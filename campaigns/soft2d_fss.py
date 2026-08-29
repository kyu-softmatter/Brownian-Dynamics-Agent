"""soft-r3 유한크기 사다리 (S17) — `ψ₆ ~ N^{-p}` 의 **형태**를 검증한다.

usage:
  python scripts/soft2d_fss.py --gates
  python scripts/soft2d_fss.py
  python scripts/soft2d_fss.py --analyze-only

## 왜 새 런인가 — `r_cut` 을 고정해야 한다

`runs/2026-07-29_soft-r3-nconv` 는 두 점(`N=100, 256`)으로 `p` 를 얻었다. 형태를
검증하려면 점이 셋 이상이어야 하는데, 상자를 줄이면 **`r_cut` 도 함께 줄어든다**:

    N     L*   r_cut   βU(r_cut) at A=10
   64   8.00   3.820   0.1794
  100  10.00   4.800   0.0904
  256  16.00   7.740   0.0216
  400  20.00   9.700   0.0110

⇒ 자연 `r_cut` 으로 사다리를 만들면 `ψ₆(N)` 기울기에 **절단오차 8배 변화가 섞인다.**
그러면 측정한 `p` 가 유한크기 효과인지 절단오차 추세인지 구별할 수 없다.

**해법: `r_cut = 3.80` 을 전 `N` 에 고정한다.** 가장 작은 상자(`N=64`, `L/2 = 4.0`)가
허용하는 값이다. 절단오차는 `0.182 kT` 로 커지지만 **모든 `N` 에서 동일**하므로
거짓 기울기를 만들 수 없다. `S16` 이 절단오차 4.2배 변화에도 관측량이 `3σ` 이내로
같음을 보였으므로 이 선택은 방어된다 — 그리고 이 런의 `N=256` 점이 자연 `r_cut`
런과 일치하는지가 그것을 한 번 더 검사한다 (`8배` 여행).

## 런 길이

`prod_tau = 30 τ_d`, 분석창 `[20, 30]`. `τ_relax ≈ 0.098 τ_d` (§8.4) 이므로
`20 τ_d` 는 완화시간의 **200배**다 — 충분하다. `80 τ_d` 를 쓰면 `A=10`·`N=400` 이
예산(`600 s/런`)을 넘는다.
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

from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,   # noqa: E402
                                       LIQUID_EXPONENT_P, hex_order_series,
                                       psi6_finite_size_exponent, zahn_phase)
from simbot.analysis.trap import aggregate_seeds                    # noqa: E402
from simbot.build import square_box_for                             # noqa: E402
from simbot.io import RunDir, provenance, write_seal                # noqa: E402
from simbot.nondim import dt_max_force, dt_max_thermal              # noqa: E402
from simbot.policy import load_policy                               # noqa: E402
from simbot.run import (Soft2DRunConfig, measure_max_force_soft2d,  # noqa: E402
                        run_soft2d)
from simbot.viz import FigureSet, plot_finite_size_scaling          # noqa: E402

import soft2d_time_series as TS                                    # noqa: E402

RUN_ID = "2026-07-29_soft-r3-fss"
SRC = REPO / "examples" / "soft-r3-fss"
DRIVERS = [Path(__file__), Path(TS.__file__)]

N_LADDER = (64, 144, 256, 400)
AMPLITUDES = (0.1, 1.0, 10.0)
#  ★ **전 N 에서 초기배치가 성공하는 시드**로 선별했다 (paired 설계).
#    `min_sep = 0.8 d` 기각표집은 시드마다 실패한다 — 실측 성공률 (시드 31–90):
#      N=64 98.3 % · N=144 95.0 % · N=256 100 % · N=400 100 %
#    (60표본이라 N 의존성은 유의하지 않다. 시드 31 은 N=144 에서 실패했다.)
#    시드가 N 마다 다르면 ψ₆(N) 비교에 다른 초기배치 앙상블이 섞인다 → 짝지어야 한다.
SEEDS = (32, 33, 34, 35)
R_CUT_FIXED = 3.80                 # ★ 전 N 공통. N=64 의 L/2 = 4.0 이 허용하는 값
PROD_TAU, N_FRAMES = 30.0, 300     # stride 0.1 tau_d
WINDOW = (20.0, 30.0)
CHI2_MAX = 3.0                     # 형태 판독 문턱 (사전등록)


def gate_table(policy) -> tuple[dict, list[dict]]:
    ts = policy.timestep
    d_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    d_F = float(ts.get("max_force_displacement_sigma", 0.005))
    dt_th = dt_max_thermal(d_th, sigma=1.0, D0=1.0)
    rows = []
    for A in AMPLITUDES:
        for N in N_LADDER:
            half = square_box_for(N) / 2.0
            if R_CUT_FIXED + TS.NLIST_BUFFER > half:
                raise SystemExit(f"⛔ N={N}: r_cut+buffer > L/2 = {half}")
            probe = Soft2DRunConfig(
                amplitude=A, n_particles=N, init=TS.INIT,
                box_shape=TS.BOX_SHAPE, r_cut=R_CUT_FIXED, r_min=TS.R_MIN,
                nlist_buffer=TS.NLIST_BUFFER, min_sep_init=TS.MIN_SEP_INIT,
                seed=SEEDS[0], dt_star=1e-6)
            m = measure_max_force_soft2d(probe)
            dt_f = dt_max_force(d_F, sigma=1.0, gamma=1.0,
                                max_force=m["max_force_star"])
            act = {k: v for k, v in (("thermal_displacement", dt_th),
                                     ("force_displacement", dt_f))
                   if v is not None}
            dom = min(act, key=act.get)
            dt = act[dom]
            rows.append({"amplitude": A, "n_particles": N, "dt_star": dt,
                         "dominant_gate": dom, "steps": int(round(PROD_TAU / dt)),
                         "max_force_star": m["max_force_star"],
                         "r_cut": m["r_cut"],
                         "beta_u_at_rcut": m["beta_u_at_rcut"],
                         "Lx": m["Lx"]})
    return {"delta_thermal": d_th, "delta_force": d_F,
            "r_cut_fixed": R_CUT_FIXED}, rows


def print_gates(rows: list[dict], policy) -> None:
    lam = float(policy.get("hardware.throughput_particle_steps_per_s", 6.3e6))
    k = policy.concurrency("default")
    eff = policy.efficiency(k)
    print(f"## 게이트 — `r_cut = {R_CUT_FIXED}` 전 N 고정 · prod {PROD_TAU:g} τ_d "
          f"/ 창 {WINDOW}\n")
    hdr = (f"{'A':>6} {'N':>5} {'L*':>6} {'r_cut':>7} {'βU(rc)':>8} "
           f"{'max|F*|':>9} {'지배':<22} {'dt*':>10} {'steps':>9} {'예상':>7}")
    print(hdr); print("-" * len(hdr))
    worst, bu = 0.0, set()
    for r in rows:
        wall = r["n_particles"] * r["steps"] / (lam * eff) * 3.4
        worst = max(worst, wall)
        bu.add(round(r["beta_u_at_rcut"], 6))
        print(f"{r['amplitude']:>6g} {r['n_particles']:>5d} {r['Lx']:>6.2f} "
              f"{r['r_cut']:>7.3f} {r['beta_u_at_rcut']:>8.4f} "
              f"{r['max_force_star']:>9.2f} {r['dominant_gate']:<22} "
              f"{r['dt_star']:>10.3g} {r['steps']:>9d} {wall:>6.0f}s")
    #  ★ 절단오차가 A 마다 하나씩(N 에 무관) 이어야 한다 — 그것이 이 설계의 요점
    per_A = {A: {round(r["beta_u_at_rcut"], 6) for r in rows
                 if r["amplitude"] == A} for A in AMPLITUDES}
    for A, s in per_A.items():
        if len(s) != 1:
            raise SystemExit(f"⛔ A={A}: βU(r_cut) 이 N 마다 다르다 {s} — "
                             f"r_cut 고정이 깨졌다")
    print(f"\n  ✅ βU(r_cut) 이 각 A 에서 N 에 무관하다: "
          f"{ {A: s.pop() for A, s in per_A.items()} }")
    budget = policy.wall_budget_s
    n = len(rows) * len(SEEDS)
    print(f"  런 {n}개 (A {len(AMPLITUDES)} × N {len(N_LADDER)} × 시드 "
          f"{len(SEEDS)}) · 동시 {k} · 최장 런 {worst:.0f} s "
          f"{'≤' if worst <= budget else '>'} 예산 {budget:.0f} s")
    if worst > budget:
        raise SystemExit("⛔ 예산 초과 — 실행하지 않고 보고한다")


def _one(args) -> dict:
    cfg, outdir = args
    return run_soft2d(Soft2DRunConfig(**cfg), outdir=Path(outdir))


def run_batch(rd: RunDir, rows: list[dict], policy) -> dict:
    jobs = []
    for r in rows:
        for s in SEEDS:
            label = f"A{r['amplitude']:g}_N{r['n_particles']}_s{s}"
            cfg = Soft2DRunConfig(
                amplitude=r["amplitude"], n_particles=r["n_particles"],
                init=TS.INIT, box_shape=TS.BOX_SHAPE, r_cut=R_CUT_FIXED,
                r_min=TS.R_MIN, nlist_buffer=TS.NLIST_BUFFER,
                min_sep_init=TS.MIN_SEP_INIT, dt_star=r["dt_star"],
                equil_tau=0.0, prod_tau=PROD_TAU, n_frames=N_FRAMES,
                seed=s, label=label)
            jobs.append((asdict(cfg), str(rd.raw / label)))
    k = policy.concurrency("default")
    print(f"\n## S5 — {len(jobs)} 런 (동시 {k})")
    t0 = time.perf_counter()
    done, failed = [], []
    with ProcessPoolExecutor(max_workers=k) as ex:
        futs = {ex.submit(_one, j): j[0]["label"] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            label = futs[fut]
            try:
                out = fut.result()
                done.append({"label": label, "wall_s": out["wall_s"],
                             "fails": out["guards"]["failures"]})
                if out["guards"]["failures"]:
                    print(f"  ⚠ {label}: {out['guards']['failures']}")
            except Exception as e:                          # noqa: BLE001
                failed.append({"label": label, "error": repr(e)})
                print(f"  ⛔ {label}: {e!r}")
            if i % 12 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} …")
    wall = time.perf_counter() - t0
    print(f"  배치 wall {wall:.1f} s · 실패 {len(failed)}")
    return {"done": len(done), "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "seeds": list(SEEDS),
            "n_ladder": list(N_LADDER), "r_cut_fixed": R_CUT_FIXED}


def analyze(rd: RunDir) -> dict:
    lo, hi = WINDOW
    out: dict = {}
    for A in AMPLITUDES:
        Ns, psi, se, loc, defect, defect_se = [], [], [], [], [], []
        for N in N_LADDER:
            dirs = sorted(p for p in rd.raw.glob(f"A{A:g}_N{N}_s*")
                          if (p / "samples.npz").exists())
            if not dirs:
                raise SystemExit(f"⛔ A={A} N={N} 런이 없다")
            per = []
            for d in dirs:
                z = np.load(d / "samples.npz")
                cfg = json.loads((d / "manifest.json").read_text())["config"]
                stride = int(z["stride"][0])
                t = np.arange(1, z["traj"].shape[0] + 1) * stride * cfg["dt_star"]
                s = hex_order_series(z["traj"], Lx=float(z["box"][0]),
                                     Ly=float(z["box"][1]), t_star=t,
                                     coord_range=(3, 12))
                m = (t >= lo) & (t < hi)
                per.append({"g": float(s.psi6_global[m].mean()),
                            "l": float(s.psi6_local[m].mean()),
                            "d": float(s.defect_fraction[m].mean())})
            ag = {k_: aggregate_seeds([p[k_] for p in per])
                  for k_ in ("g", "l", "d")}
            Ns.append(N); psi.append(ag["g"].mean); se.append(ag["g"].se)
            loc.append(ag["l"].mean)
            defect.append(ag["d"].mean); defect_se.append(ag["d"].se)
        fit = psi6_finite_size_exponent(np.array(Ns, dtype=float),
                                        np.array(psi), np.array(se))
        out[f"A{A:g}"] = {
            "amplitude": A, "zahn": zahn_phase(A), "n_ladder": Ns,
            "psi6_global": psi, "psi6_global_se": se, "psi6_local": loc,
            "defect_fraction": defect, "defect_fraction_se": defect_se,
            "n_seeds": len(dirs),
            "finite_size": {
                "p": fit.p, "p_se": fit.p_se, "eta6": fit.eta6,
                "eta6_se": fit.eta6_se, "reading": fit.reading,
                "n_points": fit.n_points, "chi2_reduced": fit.chi2_reduced,
                "residuals": list(fit.residuals), "amplitude": fit.amplitude,
                "form_verdict": fit.form_verdict(CHI2_MAX),
                "form_is_testable": fit.form_is_testable},
        }
        print(f"  A={A:<5g} p = {fit.p:.3f}±{fit.p_se:.3f}  η₆ = {fit.eta6:.2f}"
              f"±{fit.eta6_se:.2f}  χ²/dof = {fit.chi2_reduced:.2f}  "
              f"{fit.reading}\n           {fit.form_verdict(CHI2_MAX)}")
    return out


def figures(rd: RunDir, metrics: dict) -> FigureSet:
    fs = FigureSet(rd.figs)
    data, local = {}, {}
    for A in AMPLITUDES:
        m = metrics[f"A{A:g}"]
        f = m["finite_size"]
        from simbot.analysis.structure import FiniteSizeExponent
        data[A] = {"N": m["n_ladder"], "psi6": m["psi6_global"],
                   "se": m["psi6_global_se"],
                   "fit": FiniteSizeExponent(
                       p=f["p"], p_se=f["p_se"], eta6=f["eta6"],
                       eta6_se=f["eta6_se"], n_points=f["n_points"],
                       reading=f["reading"], chi2_reduced=f["chi2_reduced"],
                       residuals=tuple(f["residuals"]),
                       amplitude=f["amplitude"])}
        local[A] = {"N": m["n_ladder"], "defect": m["defect_fraction"],
                    "defect_se": m["defect_fraction_se"]}
    plot_finite_size_scaling(
        fs, data, local_data=local, exponent_liquid=LIQUID_EXPONENT_P,
        exponent_hexatic=KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
        name="01_fss_ladder")
    fs.skip("voronoi", "결함의 성격은 부모 런이 이미 냈다")
    fs.skip("early_transient", "완화는 시드 1513개로 이미 닫았다 (§8.4)")
    return fs


def check(pred: dict, metrics: dict) -> list[dict]:
    """봉인 예측 대조. 측정값은 `metrics` 에서만 읽는다."""
    got: dict = {}
    for A in AMPLITUDES:
        m = metrics[f"A{A:g}"]
        f = m["finite_size"]
        got[f"chi2_reduced__A{A:g}"] = f["chi2_reduced"]
        got[f"psi6_exponent_p__A{A:g}"] = f["p"]
        #  N=256 점을 앞 런(자연 r_cut)과 대조 — 사다리에서 뽑아낸다
        i = m["n_ladder"].index(256)
        got[f"psi6_global__A{A:g}__N256"] = m["psi6_global"][i]
    f10 = metrics["A10"]["finite_size"]
    got["eta6_minus_3sigma__A10"] = f10["eta6"] - 3.0 * f10["eta6_se"]

    rows = []
    for it in pred["items"]:
        q = it["quantity"]
        v = got.get(q)
        tol = str(it["tolerance"]).strip()
        if v is None or not np.isfinite(v):
            verdict = "NOT_EVALUATED"
        elif tol.startswith(">"):
            verdict = "PASS" if v > float(tol[1:]) else "FAIL"
        elif tol.startswith("<"):
            verdict = "PASS" if v < float(tol[1:]) else "FAIL"
        else:
            verdict = ("PASS" if abs(v - float(it["value"]))
                       <= float(tol.lstrip("±")) else "FAIL")
        rows.append({"quantity": q, "predicted": it["value"],
                     "tolerance": tol, "measured": v, "verdict": verdict,
                     "discriminates": it.get("discriminates", "")})
    return rows


def residual_diagnosis(metrics: dict, A: float) -> dict:
    """`χ²` 초과가 **휘어서**인가 **오차막대가 좁아서**인가.

    ★ 휘어 있으면 잔차 부호가 단조 패턴을 만든다 (예: `−,+,+,−` 가 아니라
      `+,+,−,−`). 흩어져 있으면 부호가 자주 바뀐다.
      `χ² ∝ 1/SE²` 이므로 SE 를 2배 과소추정하면 `χ²` 가 4배가 된다 —
      4시드 SE 는 자체 불확실성이 41 % 다
      ([[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]]).
    """
    m = metrics[f"A{A:g}"]
    f = m["finite_size"]
    r = np.asarray(f["residuals"])
    se = np.asarray(m["psi6_global_se"])
    y = np.asarray(m["psi6_global"])
    z = r / (se / y)                              # 로그 잔차 / 로그 오차
    signs = np.sign(r)
    flips = int(np.sum(signs[:-1] != signs[1:]))
    #  SE 를 몇 배로 늘리면 χ²/dof = 1 이 되는가
    inflate = float(np.sqrt(f["chi2_reduced"])) if np.isfinite(
        f["chi2_reduced"]) else float("nan")
    return {"amplitude": A, "z_residuals": [float(x) for x in z],
            "sign_flips": flips, "n_points": int(r.size),
            "relative_se": [float(x) for x in se / y],
            "se_inflation_for_chi2_unity": inflate,
            "reading": ("휘어 있다 (형태 문제)" if flips <= 1
                        else "흩어져 있다 (오차막대 과소추정 쪽)")}


def main() -> int:
    policy = load_policy()
    thresholds, rows = gate_table(policy)
    print_gates(rows, policy)
    if "--gates" in sys.argv:
        return 0

    pred = yaml.safe_load((SRC / "prediction.yaml").read_text())
    rd = RunDir.create(REPO / "runs", RUN_ID)
    if "--analyze-only" not in sys.argv:
        shutil.copy2(SRC / "prediction.yaml", rd.file("prediction"))
        rd.write_json("prediction_json", pred)
        rd.write_json("spec", {"source": "scripts/soft2d_fss.py",
                               "provenance": provenance(DRIVERS),
                               "n_ladder": list(N_LADDER),
                               "amplitudes": list(AMPLITUDES),
                               "seeds": list(SEEDS),
                               "r_cut_fixed": R_CUT_FIXED,
                               "prod_tau": PROD_TAU, "n_frames": N_FRAMES,
                               "window": list(WINDOW), "chi2_max": CHI2_MAX,
                               "thresholds": thresholds, "gates": rows})
        seal = write_seal(rd)
        print(f"\n  🔒 봉인 {seal.name} — "
              f"{len(seal.read_text().splitlines())}개 문서 (실행 전)")
        rd.write_json("manifest", run_batch(rd, rows, policy))

    print("\n## S7 — 유한크기 사다리 (4점)")
    metrics = analyze(rd)
    metrics["_provenance_at_analysis"] = provenance(DRIVERS)
    fs = figures(rd, metrics)
    rd.write("figures", fs.figures_md())

    checks = check(pred, metrics)
    diag = {f"A{A:g}": residual_diagnosis(metrics, A) for A in AMPLITUDES}
    metrics["_checks"] = checks
    metrics["_residual_diagnosis"] = diag
    rd.write_json("metrics", metrics)

    print("\n## 예측 대조 (판정은 제안이다 — confirmed_by: null)")
    hdr = f"{'항목':<32} {'예측':>10} {'허용':>10} {'측정':>10} 판정"
    print(hdr); print("-" * (len(hdr) + 6))
    for c in checks:
        print(f"{c['quantity']:<32} {float(c['predicted']):>10.4g} "
              f"{c['tolerance']:>10} {c['measured']:>10.4g} {c['verdict']}")
    n_p = sum(1 for c in checks if c["verdict"] == "PASS")
    print(f"\n  PASS {n_p} · FAIL {len(checks) - n_p}")

    fails = [c for c in checks if c["verdict"] == "FAIL"]
    if fails:
        print("\n## FAIL 진단 — 휘어서인가 오차막대가 좁아서인가")
        for c in fails:
            if not c["quantity"].startswith("chi2_reduced__"):
                continue
            A = float(c["quantity"].split("__A")[1])
            d = diag[f"A{A:g}"]
            print(f"  A={A:g}: 잔차/σ = " +
                  " ".join(f"{x:+.2f}" for x in d["z_residuals"]))
            print(f"    부호 전환 {d['sign_flips']}/{d['n_points']-1} → "
                  f"**{d['reading']}**")
            print(f"    상대 SE = " +
                  " ".join(f"{x:.1%}" for x in d["relative_se"]))
            print(f"    SE 를 {d['se_inflation_for_chi2_unity']:.2f}배로 늘리면 "
                  f"χ²/dof = 1 이 된다 (4시드 SE 의 자체 불확실성은 41 %)")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
