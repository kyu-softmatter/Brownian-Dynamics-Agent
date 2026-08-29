"""soft-r3 완화시간 — 시드를 늘려 `τ(A=1)` vs `τ(A=0.1)` 을 판정한다.

usage:
  python scripts/soft2d_relax_seeds.py --power        # 설계 검정력만 (실행 안 함)
  python scripts/soft2d_relax_seeds.py               # 사전등록 k 로 실행 + 판정
  python scripts/soft2d_relax_seeds.py --analyze-only

## 왜 별도 스크립트인가

`soft2d_time_series.py` 를 고치면 그 파일의 `driver_hash` 가 바뀐다. 그 해시는
**이미 커밋된 런의 봉인된 `03_spec.yaml` 에 박혀 있다** — 고치면 그 런의 산출물이
가리키는 드라이버와 디스크의 파일이 어긋난다. 그래서 건드리지 않고 새로 만든다.

기하·게이트·러너는 그 스크립트에서 **import 한다** (같은 수를 쓴다는 보장).
그러므로 `provenance(driver=[이 파일, 그 파일])` 로 **둘 다** 해싱한다 —
하나만 잡으면 "code_hash + driver_hash 가 전부를 덮는다"가 거짓이 된다.

## 이 스크립트가 조심하는 것

**① 표본 수를 미리 고정한다.** 유의해질 때까지 시드를 늘리는 것은 optional stopping
   이다. `k` 는 `examples/soft-r3-relax-seeds/prediction.yaml` 에 봉인되어 있고
   이 스크립트는 그 파일에서 읽는다 — 코드에 다시 쓰지 않는다.

**② 오차는 시드 앙상블에서 온다.** `curve_fit` 공분산은 *평균 곡선 하나*의 적합
   불확실성이다. `bootstrap_relaxation_over_seeds` 로 시드를 복원추출해 다시 적합한다.

**③ "구별 안 됨"을 "같다"로 쓰지 않는다.** 판정 문구를 예측 파일이 소유한다.
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from simbot.analysis.structure import (bootstrap_relaxation_over_seeds,  # noqa: E402
                                       fit_relaxation, hex_order_series)
from simbot.estimators import seeds_for_target_sigma                    # noqa: E402
from simbot.io import RunDir, provenance, write_seal                    # noqa: E402
from simbot.policy import load_policy                                   # noqa: E402
from simbot.run import Soft2DRunConfig, run_soft2d                      # noqa: E402

#  ★ 기하·게이트·상수를 재구현하지 않는다 — 같은 수를 쓴다는 보장이 필요하다
import soft2d_time_series as TS                                         # noqa: E402

SRC = REPO / "examples" / "soft-r3-relax-seeds"
AMPLITUDES = (0.1, 1.0)          # ★ 판정 대상 두 개만. A=10 은 이미 4.7σ 로 갈렸다
DRIVERS = [Path(__file__), Path(TS.__file__)]

#  순차 설계의 단계들. **k 는 코드가 아니라 봉인된 예측 파일이 소유한다** —
#  코드에 다시 쓰면 두 곳이 갈라지고, 그러면 "미리 고정했다"를 증명할 수 없다.
STAGES = {
    1: {"pred": "prediction.yaml",
        "run_id": "2026-07-29_soft-r3-relax-seeds"},
    2: {"pred": "prediction_stage2.yaml",
        "run_id": "2026-07-29_soft-r3-relax-seeds-stage2"},
}


def stage_from_argv() -> int:
    if "--stage" in sys.argv:
        s = int(sys.argv[sys.argv.index("--stage") + 1])
        if s not in STAGES:
            raise SystemExit(f"⛔ 단계 {s} 없음. 있는 것: {sorted(STAGES)}")
        return s
    return 1


def load_prediction(stage: int) -> dict:
    return yaml.safe_load((SRC / STAGES[stage]["pred"]).read_text())


def power_table(pred: dict, stage: int) -> dict:
    """사전등록된 설계 검정력을 **다시 계산해서** 예측 파일과 일치하는지 확인한다."""
    dp = pred["design_power"]
    r3 = seeds_for_target_sigma(diff=dp["prior_diff_tau_d"],
                                se_diff=dp["prior_se_diff_tau_d"],
                                k_current=dp["prior_k"], n_sigma=3.0)
    k = int(dp["k_chosen"])
    sigma_at_k = r3["sigma_now"] * np.sqrt(k / dp["prior_k"])

    print(f"## 설계 검정력 — **실행 전에** 계산했다  (단계 {stage})\n")
    print(f"  선행 관측 (k={dp['prior_k']}):  diff = {dp['prior_diff_tau_d']:+.5f} "
          f"± {dp['prior_se_diff_tau_d']:.5f} τ_d = {r3['sigma_now']:.2f}σ")
    print(f"  3σ 에 필요한 k = {r3['k_needed_int']}  (t={r3['t_quantile']:.3f})")
    print(f"  ⇒ 사전등록 k = {k}  →  예상 {sigma_at_k:.2f}σ (차이가 참일 때)")

    #  예측 파일과 어긋나면 멈춘다 — 봉인된 설계와 코드가 갈라지면 안 된다
    if int(dp["k_for_3sigma"]) != r3["k_needed_int"]:
        raise SystemExit(
            f"⛔ 봉인된 k_for_3sigma={dp['k_for_3sigma']} 와 재계산 "
            f"{r3['k_needed_int']} 가 다르다 — 설계가 코드와 갈라졌다")
    if stage == 1:
        sigma_at_64 = r3["sigma_now"] * np.sqrt(64 / dp["prior_k"])
        print(f"  ★ k=64 로는 {sigma_at_64:.2f}σ — 3σ 에 도달하지 못한다")
        if abs(float(dp["k_at_64_sigma"]) - sigma_at_64) > 0.01:
            raise SystemExit(f"⛔ k_at_64_sigma 불일치: {dp['k_at_64_sigma']} vs "
                             f"{sigma_at_64:.2f}")
    else:
        pm = pred["stage1_postmortem"]
        print(f"\n  1단계 부검: 예상 {pm['predicted_sigma']}σ → 관측 "
              f"{pm['observed_sigma']}σ")
        print(f"    ① 차이가 {pm['cause_1_biased_diff']['change_pct']:+.1f} % 줄었다 "
              f"(저시드 편향)")
        print(f"    ② curve_fit SE 가 "
              f"{pm['cause_2_underestimated_se']['se_ratio_bootstrap_over_fit_A1']}"
              f"배 과소추정")
        if abs(float(dp["k_chosen_expected_sigma"]) - sigma_at_k) > 0.02:
            raise SystemExit(f"⛔ k_chosen_expected_sigma 불일치")
        print(f"  ★ 3단계는 없다 (no_stage_3="
              f"{pred['decision_rule'].get('no_stage_3')})")
    print("  ✅ 봉인된 설계 검정력이 재계산과 일치한다")
    return {"k": k, "sigma_expected_at_k": float(sigma_at_k),
            "k_for_3sigma": r3["k_needed_int"], "t_quantile": r3["t_quantile"],
            "stage": stage}


def _one(args) -> dict:
    cfg_dict, outdir = args
    return run_soft2d(Soft2DRunConfig(**cfg_dict), outdir=Path(outdir))


def run_batch(rd: RunDir, rows: list[dict], k: int, policy) -> dict:
    seeds = tuple(range(5, 5 + k))          # 선행 런과 같은 시작점 (5–20 을 포함)
    outroot = rd.path / "raw_early"
    outroot.mkdir(parents=True, exist_ok=True)
    jobs = []
    for r in rows:
        if r["amplitude"] not in AMPLITUDES:
            continue
        for s in seeds:
            label = f"A{r['amplitude']:g}_s{s}"
            cfg = Soft2DRunConfig(
                amplitude=r["amplitude"], n_particles=TS.N_PARTICLES,
                init=TS.INIT, box_shape=TS.BOX_SHAPE, r_min=TS.R_MIN,
                nlist_buffer=TS.NLIST_BUFFER, min_sep_init=TS.MIN_SEP_INIT,
                dt_star=r["dt_star"], equil_tau=0.0, prod_tau=TS.EARLY_TAU,
                n_frames=TS.EARLY_FRAMES, seed=s, label=label)
            jobs.append((asdict(cfg), str(outroot / label)))

    kk = policy.concurrency("default")
    print(f"\n## S5 — {len(jobs)} 런 (A {len(AMPLITUDES)}개 × 시드 {k}개, 동시 {kk})")
    t0 = time.perf_counter()
    done, failed = [], []
    with ProcessPoolExecutor(max_workers=kk) as ex:
        futs = {ex.submit(_one, j): j[0]["label"] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            label = futs[fut]
            try:
                out = fut.result()
                done.append({"label": label, "wall_s": out["wall_s"],
                             "fails": out["guards"]["failures"]})
            except Exception as e:                       # noqa: BLE001
                failed.append({"label": label, "error": repr(e)})
                print(f"  ⛔ {label}: {e!r}")
            if i % 100 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} …")
    wall = time.perf_counter() - t0
    n_fail_guard = sum(1 for d in done if d["fails"])
    print(f"  배치 wall {wall:.1f} s · 실패 {len(failed)} · 가드위반 {n_fail_guard}")
    return {"done": len(done), "failed": failed, "batch_wall_s": wall,
            "concurrency": kk, "n_jobs": len(jobs), "seeds": list(seeds),
            "guard_violations": n_fail_guard}


def analyze(rd: RunDir, pred: dict) -> dict:
    rule = pred["decision_rule"]
    root = rd.path / "raw_early"
    out: dict = {}
    for A in AMPLITUDES:
        dirs = sorted(p for p in root.glob(f"A{A:g}_s*")
                      if (p / "samples.npz").exists())
        if not dirs:
            raise SystemExit(f"⛔ A={A} 런이 없다")
        curves, t_ref = [], None
        for d in dirs:
            z = np.load(d / "samples.npz")
            cfg = json.loads((d / "manifest.json").read_text())["config"]
            stride = int(z["stride"][0])
            n = z["traj"].shape[0]
            t = np.concatenate([[0.0],
                                np.arange(1, n + 1) * stride * cfg["dt_star"]])
            frames = np.concatenate(
                [z["init_pos"][None].astype(np.float32), z["traj"]])
            s = hex_order_series(frames, Lx=float(z["box"][0]),
                                 Ly=float(z["box"][1]), t_star=t,
                                 coord_range=(3, 12))
            curves.append(s.defect_fraction)
            t_ref = t if t_ref is None else t_ref
        mat = np.array(curves)
        tail = mat[:, t_ref >= 0.5 * t_ref[-1]]
        noise = float(mat.mean(axis=0)[t_ref >= 0.5 * t_ref[-1]].std(ddof=1))
        fit = fit_relaxation(t_ref, mat.mean(axis=0), noise=noise)
        boot = bootstrap_relaxation_over_seeds(
            t_ref, mat, n_resample=int(rule["n_resample"]),
            seed=int(rule["bootstrap_seed"]), noise=noise)
        out[f"A{A:g}"] = {
            "amplitude": A, "n_seeds": len(dirs),
            "tau": fit.tau, "tau_se_fit": fit.tau_se,
            "tau_se_bootstrap": boot["tau_se_bootstrap"],
            "tau_ci95": list(boot["tau_ci95"]),
            "se_ratio_bootstrap_over_fit": boot.get(
                "se_ratio_bootstrap_over_fit"),
            "bootstrap_converged": boot["n_converged"],
            "relax_amplitude": fit.amplitude,
            "relax_converged": bool(fit.converged), "relax_note": fit.note,
            "r_squared": fit.r_squared, "defect_frame_sd": noise,
            "defect_at_t0": float(mat[:, 0].mean()),
            "defect_tail_mean": float(tail.mean()),
        }
        print(f"  A={A:<5g} k={len(dirs):<4d} τ = {fit.tau:.5f}  "
              f"SE_fit {fit.tau_se:.5f} · SE_boot {boot['tau_se_bootstrap']:.5f} "
              f"(비 {boot.get('se_ratio_bootstrap_over_fit', float('nan')):.2f})  "
              f"폭/잡음 {abs(fit.amplitude)/noise:.1f}×")

    a, b = out["A0.1"], out["A1"]
    for tag, key in (("bootstrap", "tau_se_bootstrap"), ("fit", "tau_se_fit")):
        se = float(np.hypot(a[key], b[key]))
        diff = b["tau"] - a["tau"]
        out[f"comparison_{tag}"] = {
            "diff": diff, "se_diff": se,
            "sigma": abs(diff) / se if se else float("nan")}
    return out


def main() -> int:
    policy = load_policy()
    stage = stage_from_argv()
    pred = load_prediction(stage)
    geo = TS.geometry()
    thresholds, rows = TS.gate_table(policy)
    power = power_table(pred, stage)
    if "--power" in sys.argv:
        return 0

    rd = RunDir.create(REPO / "runs", STAGES[stage]["run_id"])
    if "--analyze-only" not in sys.argv:
        shutil.copy2(SRC / STAGES[stage]["pred"], rd.file("prediction"))
        rd.write_json("prediction_json", pred)
        rd.write_json("spec", {
            "source": "scripts/soft2d_relax_seeds.py",
            "stage": stage,
            "provenance": provenance(DRIVERS),
            "parent_run": pred["parent_run"],
            "amplitudes": list(AMPLITUDES), "k_seeds": power["k"],
            "geometry": geo, "thresholds": thresholds,
            "gates": [r for r in rows if r["amplitude"] in AMPLITUDES],
            "early_tau": TS.EARLY_TAU, "early_frames": TS.EARLY_FRAMES,
            "design_power": power})
        seal = write_seal(rd)
        print(f"\n  🔒 봉인 {seal.name} — "
              f"{len(seal.read_text().splitlines())}개 문서 (실행 전)")
        batch = run_batch(rd, rows, power["k"], policy)
        rd.write_json("manifest", batch)

    print("\n## S7 — 완화시간 (부트스트랩 오차)")
    res = analyze(rd, pred)
    res["_provenance_at_analysis"] = provenance(DRIVERS)
    res["_design_power"] = power
    rd.write_json("metrics", res)

    rule = pred["decision_rule"]
    thr = float(rule["threshold_sigma"])
    cb, cf = res["comparison_bootstrap"], res["comparison_fit"]
    print(f"\n## 판정 — 사전등록 규칙 (SE 출처: {rule['se_source']}, 문턱 {thr:g}σ)\n")
    print(f"  {'SE 출처':<12} {'diff':>10} {'SE_diff':>10} {'σ':>7}")
    print("  " + "-" * 42)
    for tag, c in (("bootstrap ★", cb), ("curve_fit", cf)):
        print(f"  {tag:<12} {c['diff']:>+10.5f} {c['se_diff']:>10.5f} "
              f"{c['sigma']:>7.2f}")
    verdict = (rule["verdict_if_above"] if cb["sigma"] > thr
               else rule["verdict_if_below"])
    print(f"\n  → {verdict}")
    print(f"  (예상했던 것: {power['sigma_expected_at_k']:.2f}σ · "
          f"실제: {cb['sigma']:.2f}σ)")
    print("\n  proposed_by: agent · confirmed_by: null")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
