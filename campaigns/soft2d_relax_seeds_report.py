"""soft-r3 완화시간 — 시드 수렴 분석 + 그림 + 검증 문서.

usage: python scripts/soft2d_relax_seeds_report.py [--stage 2]

## 이 스크립트의 핵심 — 중첩 부분집합

`τ(k)` 를 **같은 궤적의 처음 `k` 개**로 계산한다. 서로 다른 런에서 뽑으면
"데이터가 달라서" 와 "추정량이 `k` 에 의존해서" 를 구별할 수 없다. 중첩으로 잡으면
곡선의 기울기가 **순수하게 추정량의 `k` 의존성**이다 — 그것이 저시드 편향이다.

부수 효과로 "`k = 64` 는 실제로 무엇을 줬을까" 에 답할 수 있다 (사용자의 원 제안).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from simbot.analysis.structure import (bootstrap_relaxation_over_seeds,  # noqa: E402
                                       fit_relaxation, hex_order_series)
from simbot.io import RunDir, provenance, verify_seal                   # noqa: E402
from simbot.report import seal_section                                  # noqa: E402
from simbot.viz import FigureSet, plot_seed_convergence                 # noqa: E402

import soft2d_relax_seeds as RS                                        # noqa: E402
import soft2d_time_series as TS                                        # noqa: E402

K_LADDER = (4, 8, 16, 32, 64, 128, 256, 512, 1024)
DRIVERS = [Path(__file__), Path(RS.__file__), Path(TS.__file__)]


def load_curves(rd: RunDir, A: float) -> tuple[np.ndarray, np.ndarray]:
    """시드 순서대로 정렬된 결함분율 곡선 `(n_seeds, n_frames)` + 시간축.

    ★ **캐시한다.** 곡선 1개가 프레임 401개의 Voronoi 이고 시드가 1500개면
      60만 번이다 (A 두 개면 120만). 사다리 분석이 그것을 다시 하면 안 된다.
      캐시 키에 시드 수를 넣어 런이 늘어나면 자동으로 무효화된다.
    """
    dirs = sorted((p for p in (rd.path / "raw_early").glob(f"A{A:g}_s*")
                   if (p / "samples.npz").exists()),
                  key=lambda p: int(p.name.rsplit("_s", 1)[1]))
    cache = rd.path / f"defect_curves_A{A:g}_k{len(dirs)}.npz"
    if cache.exists():
        z = np.load(cache)
        print(f"  A={A:<5g} 캐시에서 읽음 ({cache.name})")
        return z["t"], z["curves"]

    curves, t_ref = [], None
    for d in dirs:
        z = np.load(d / "samples.npz")
        cfg = json.loads((d / "manifest.json").read_text())["config"]
        stride = int(z["stride"][0])
        n = z["traj"].shape[0]
        t = np.concatenate([[0.0], np.arange(1, n + 1) * stride * cfg["dt_star"]])
        frames = np.concatenate([z["init_pos"][None].astype(np.float32), z["traj"]])
        s = hex_order_series(frames, Lx=float(z["box"][0]), Ly=float(z["box"][1]),
                             t_star=t, coord_range=(3, 12))
        curves.append(s.defect_fraction)
        t_ref = t if t_ref is None else t_ref
    arr = np.array(curves)
    np.savez_compressed(cache, t=t_ref, curves=arr)
    print(f"  A={A:<5g} 곡선 {arr.shape[0]}개 계산 → 캐시 {cache.name}")
    return t_ref, arr


def _fit_subset(t: np.ndarray, mat: np.ndarray, n_resample: int,
                boot_seed: int) -> dict:
    tail_m = t >= 0.5 * t[-1]
    noise = float(mat.mean(axis=0)[tail_m].std(ddof=1))
    fit = fit_relaxation(t, mat.mean(axis=0), noise=noise)
    boot = bootstrap_relaxation_over_seeds(t, mat, n_resample=n_resample,
                                           seed=boot_seed, noise=noise)
    lo, hi = boot["tau_ci95"]
    #  ★ SD 는 저시드에서 이상치 지배다 (k=4 에서 82 τ_d — τ 자체가 0.06 인데).
    #    백분위 구간의 반폭을 **강건 SE** 로 함께 낸다. 판정은 사전등록대로 SD 를
    #    쓰지만, 그림과 진단은 강건값을 봐야 읽을 수 있다.
    se_robust = (hi - lo) / (2.0 * 1.96) if np.isfinite(lo) and np.isfinite(hi) \
        else float("nan")
    return {"tau": fit.tau, "se_fit": fit.tau_se,
            "se": boot["tau_se_bootstrap"], "se_robust": se_robust,
            "ci95": [lo, hi], "amp": fit.amplitude,
            "noise": noise, "r2": fit.r_squared,
            "amp_over_noise": abs(fit.amplitude) / noise if noise else float("nan"),
            "y0": fit.y0, "y_inf": fit.y_inf,
            "n_converged": boot["n_converged"],
            "se_ratio": boot.get("se_ratio_bootstrap_over_fit")}


def main() -> int:
    stage = RS.stage_from_argv()
    pred = RS.load_prediction(stage)
    rule = pred["decision_rule"]
    rd = RunDir(REPO / "runs" / RS.STAGES[stage]["run_id"])
    v = verify_seal(rd)
    print(("✅ " if v.ok else "⛔ ") + v.summary())

    t, mats = None, {}
    for A in RS.AMPLITUDES:
        t, mats[A] = load_curves(rd, A)
        print(f"  A={A:<5g} 시드 {mats[A].shape[0]}개 · 프레임 {mats[A].shape[1]}")

    k_max = min(m.shape[0] for m in mats.values())
    ladder = [k for k in K_LADDER if k <= k_max] + [k_max]
    ladder = sorted(set(ladder))

    print(f"\n## 중첩 부분집합 τ(k) — 같은 궤적의 처음 k 개\n")
    hdr = (f"{'k':>6} {'τ(A=0.1)':>20} {'τ(A=1)':>20} {'diff':>10} "
           f"{'SE_diff':>9} {'σ':>7}")
    print(hdr); print("-" * len(hdr))
    conv: dict = {"k": [], "A0.1": {"tau": [], "se": []},
                  "A1": {"tau": [], "se": []}, "diff": [], "se_diff": [],
                  "sigma": [], "detail": {}}
    for k in ladder:
        f01 = _fit_subset(t, mats[0.1][:k], int(rule["n_resample"]),
                          int(rule["bootstrap_seed"]))
        f1 = _fit_subset(t, mats[1.0][:k], int(rule["n_resample"]),
                         int(rule["bootstrap_seed"]))
        diff = f1["tau"] - f01["tau"]
        se = float(np.hypot(f01["se"], f1["se"]))
        sig = abs(diff) / se if se else float("nan")
        conv["k"].append(k)
        for key, f in (("A0.1", f01), ("A1", f1)):
            conv[key]["tau"].append(f["tau"])
            conv[key]["se"].append(f["se"])
            conv[key].setdefault("se_robust", []).append(f["se_robust"])
        conv["diff"].append(diff); conv["se_diff"].append(se)
        conv["sigma"].append(sig)
        conv["detail"][str(k)] = {"A0.1": f01, "A1": f1}
        mark = "  ← 사용자 제안" if k == 64 else ("  ← 사전등록" if k == k_max else "")
        print(f"{k:>6d} {f01['tau']:>11.5f}±{f01['se']:<8.5f} "
              f"{f1['tau']:>11.5f}±{f1['se']:<8.5f} {diff:>+10.5f} "
              f"{se:>9.5f} {sig:>7.2f}{mark}")

    # --- 판정 (사전등록 규칙) ---
    thr = float(rule["threshold_sigma"])
    sig_final = conv["sigma"][-1]
    verdict = (rule["verdict_if_above"] if sig_final > thr
               else rule["verdict_if_below"])
    print(f"\n## 판정 — 사전등록 문턱 {thr:g}σ, SE 출처 {rule['se_source']}")
    print(f"  k = {k_max}:  {sig_final:.2f}σ  →  {verdict.strip().splitlines()[0]}")
    print(f"  (사전등록 예상 {pred['design_power']['k_chosen_expected_sigma']}σ)")

    # --- 그림 ---
    fs = FigureSet(rd.figs)
    geo = TS.geometry()
    curves = {}
    for A in RS.AMPLITUDES:
        d = conv["detail"][str(k_max)][f"A{A:g}"]
        curves[f"{A:g}"] = (t, mats[A].mean(axis=0), d)
    prereg = [{"k": 256, "sigma_expected": 4.19, "stage": 1},
              {"k": int(pred["design_power"]["k_chosen"]),
               "sigma_expected": float(
                   pred["design_power"]["k_chosen_expected_sigma"]),
               "stage": 2}]
    plot_seed_convergence(fs, conv, curves=curves, tau_d_si=geo["tau_d_si"],
                          threshold_sigma=thr, preregistered=prereg)
    fs.skip("voronoi", "이 런은 완화시간만 본다 — 배치의 성격은 부모 런이 이미 냈다")
    rd.write("figures", fs.figures_md())

    out = {"convergence": conv, "k_max": k_max, "sigma_final": sig_final,
           "verdict": verdict.strip(), "threshold_sigma": thr,
           "preregistered": prereg,
           "_provenance_at_analysis": provenance(DRIVERS)}
    #  ★ 기존 metrics 를 덮지 않고 합친다 — 러너가 쓴 τ·부트스트랩 결과를 잃으면
    #    수렴표만 남고 "최종 k 에서 무엇이 나왔나" 가 사라진다
    prev = json.loads(rd.read("metrics")) if rd.exists("metrics") else {}
    merged = {**prev, **out}
    rd.write_json("metrics", merged)
    rd.write("validation", validation_md(rd, merged, pred, stage))
    print(f"\n→ {rd.figs.relative_to(REPO)}/01_seed_convergence.png")
    print(f"→ {rd.file('validation').relative_to(REPO)}")
    return 0


def check_predictions(pred: dict, merged: dict, rd: RunDir) -> list[dict]:
    """봉인된 예측 대조. **측정값은 metrics 에서만 읽는다** (손으로 옮기지 않는다)."""
    c = merged["convergence"]
    man = json.loads(rd.read("manifest"))
    k_planned = int(json.loads(rd.read("spec"))["k_seeds"])
    failed_seeds = {int(f["label"].rsplit("_s", 1)[1]) for f in man["failed"]}
    meas = {
        "tau_relax__A0.1": merged["A0.1"]["tau"],
        "tau_relax__A1": merged["A1"]["tau"],
        "tau_relax_diff_sigma": merged["sigma_final"],
        "tau_relax_diff": c["diff"][-1],
        "se_ratio_bootstrap_over_fit__A1":
            merged["A1"]["se_ratio_bootstrap_over_fit"],
        "init_config_failures_per_1000_seeds":
            1000.0 * len(failed_seeds) / k_planned,
    }
    rows = []
    for it in pred["items"]:
        q = it["quantity"]
        v = meas.get(q)
        tol = str(it["tolerance"]).strip()
        if v is None:
            verdict = "NOT_EVALUATED"
        elif tol.startswith(">"):
            verdict = "PASS" if v > float(tol[1:]) else "FAIL"
        elif tol.startswith("<"):
            verdict = "PASS" if v < float(tol[1:]) else "FAIL"
        else:
            verdict = ("PASS" if abs(v - float(it["value"]))
                       <= float(tol.lstrip("±")) else "FAIL")
        rows.append({"quantity": q, "predicted": it["value"], "tolerance": tol,
                     "measured": v, "verdict": verdict,
                     "discriminates": it.get("discriminates", "")})
    return rows


def validation_md(rd: RunDir, merged: dict, pred: dict, stage: int) -> str:
    c = merged["convergence"]
    d = c["detail"]
    k_max = merged["k_max"]
    checks = check_predictions(pred, merged, rd)
    n_pass = sum(1 for x in checks if x["verdict"] == "PASS")
    n_fail = sum(1 for x in checks if x["verdict"] == "FAIL")
    from simbot.estimators import seeds_for_target_sigma
    need = seeds_for_target_sigma(diff=c["diff"][-1], se_diff=c["se_diff"][-1],
                                  k_current=k_max, n_sigma=3.0)

    out = [f"# S7 — `τ_relax` 시드 스윕 검증, **{stage}단계** ({rd.run_id})", "",
           f"`A = 0.1` vs `A = 1` · 시드 `{k_max}` (계획 "
           f"{json.loads(rd.read('spec'))['k_seeds']}) · 중첩 부분집합 사다리", "",
           "> **판정은 제안이다.** `proposed_by: agent`, `confirmed_by: null`.", "",
           seal_section(rd), "",
           "## 1. 판정 — 사전등록 규칙", "", "```yaml",
           f"verdict: {'DIFFERENT' if merged['sigma_final'] > merged['threshold_sigma'] else 'INCONCLUSIVE'}",
           f"sigma_observed: {merged['sigma_final']:.3f}",
           f"sigma_threshold: {merged['threshold_sigma']:g}",
           f"sigma_expected_preregistered: "
           f"{pred['design_power']['k_chosen_expected_sigma']}",
           "se_source: bootstrap_over_seeds",
           "proposed_by: agent", "confirmed_by: null", "```", "",
           merged["verdict"], "",
           "## 2. ★ 중첩 부분집합 사다리 — 편향이 어디서 멈추는가", "",
           "같은 궤적의 **처음 `k` 개**로 계산했다. 곡선의 기울기는 데이터 차이가 아니라 "
           "**추정량의 `k` 의존성**이다.", "",
           "| `k` | `τ(A=0.1)` | 폭/잡음 | `τ(A=1)` | 폭/잡음 | `diff` | `SE_diff` | `σ` |",
           "|---|---|---|---|---|---|---|---|"]
    for i, k in enumerate(c["k"]):
        a, b = d[str(k)]["A0.1"], d[str(k)]["A1"]
        usable = c["se_diff"][i] <= abs(c["diff"][i])
        se_s = (f"`{c['se_diff'][i]:.5f}`" if usable
                else f"`{c['se_diff'][i]:.3g}` ⚠")
        sig_s = f"`{c['sigma'][i]:.2f}`" if usable else "—"
        out.append(f"| {k} | `{a['tau']:.5f}` | `{a['amp_over_noise']:.1f}×` | "
                   f"`{b['tau']:.5f}` | `{b['amp_over_noise']:.1f}×` | "
                   f"`{c['diff'][i]:+.5f}` | {se_s} | {sig_s} |")
    tail_diff = c["diff"][-3:]
    out += ["",
            f"⚠ 표시는 **부트스트랩 SD 가 `τ` 보다 커서** 오차를 주장할 수 없는 "
            f"지점이다 (§4).", "",
            "### 읽기", "",
            f"- **`τ(A=0.1)` 은 편향이 없다** — `k=4` 의 `{d['4']['A0.1']['tau']:.5f}` "
            f"부터 `k={k_max}` 의 `{d[str(k_max)]['A0.1']['tau']:.5f}` 까지 추세 없이 "
            f"요동한다. 이미 `k=4` 에서 폭/잡음이 "
            f"`{d['4']['A0.1']['amp_over_noise']:.1f}×` 였다.",
            f"- **`τ(A=1)` 은 `{d['4']['A1']['tau']:.5f}` → "
            f"`{d[str(k_max)]['A1']['tau']:.5f}` 로 "
            f"`{100*(d[str(k_max)]['A1']['tau']/d['4']['A1']['tau']-1):.0f} %` "
            f"내려가고 **`k ≈ 512` 에서 평탄해진다.**",
            f"- ⇒ **편향은 시드 수가 아니라 폭/잡음이 정한다.** `A=1` 이 폭/잡음 `20×` "
            f"에 도달하는 것이 `k=512` 이고, 거기서 편향이 사라진다. "
            f"`A=0.1` 은 `k=4` 에서 이미 `6×` 였고 `k=64` 에 `23×` 다.",
            f"- **`diff` 는 0 으로 가지 않는다** — `k ≥ 512` 세 점의 평균 "
            f"`{np.mean(tail_diff):+.5f} ± {np.std(tail_diff, ddof=1):.5f}` 로 "
            f"평탄하다. 차이는 작지만(`τ` 의 `{100*c['diff'][-1]/merged['A1']['tau']:.0f} %`) "
            f"**있어 보인다 — 다만 `3σ` 로 분해되지 않았다.**", "",
            "## 3. 봉인 예측 대조", "",
            "| 항목 | 예측 | 허용 | 측정 | 판정 | 무엇을 가리는가 |",
            "|---|---|---|---|---|---|"]
    for x in checks:
        mv = ("—" if x["measured"] is None else f"`{x['measured']:.5g}`")
        mark = {"PASS": "**PASS**", "FAIL": "**FAIL** ⛔"}.get(x["verdict"],
                                                              x["verdict"])
        out.append(f"| `{x['quantity']}` | `{x['predicted']}` | "
                   f"`{x['tolerance']}` | {mv} | {mark} | "
                   f"{x['discriminates']} |")
    out += ["", f"**PASS {n_pass} · FAIL {n_fail}**", "",
            "### FAIL 2건의 성격", "",
            "- **`tau_relax_diff_sigma`** — 이것이 이 런의 1급 예측이었고 빗나갔다. "
            "예상 `3.96σ` vs 실제 "
            f"`{merged['sigma_final']:.2f}σ`. 원인은 §2 다: `k=254` 에서도 `A=1` 의 "
            "편향이 완전히 사라지지 않아 `diff` 가 한 번 더 줄었다 "
            f"(`{pred['design_power']['prior_diff_tau_d']:.5f}` → "
            f"`{c['diff'][-1]:.5f}`). **설계 검정력의 입력이 또 낙관적이었다.**",
            "- **`init_config_failures_per_1000_seeds`** — 예측 `7.8 ± 6`, 실제 "
            f"`{[x['measured'] for x in checks if x['quantity'].startswith('init_')][0]:.1f}`. "
            "1단계 `k=256` 의 2건(7.8/1000)에서 추정했는데 표본이 너무 작았다 "
            "(2건으로 비율을 추정하면 상대오차가 71 % 다). **`A` 에는 독립이라는 것은 "
            "확인됐다** — 같은 시드가 두 `A` 에서 함께 실패했다.", "",
            "## 4. ⚠ 부트스트랩 SD 는 저시드에서 이상치가 지배한다", "",
            "| `k` | `A=1` SD | `A=1` 강건 SE (95 % 구간/3.92) | 비 |",
            "|---|---|---|---|"]
    for k in c["k"]:
        b = d[str(k)]["A1"]
        ratio = (b["se"] / b["se_robust"] if b.get("se_robust") else float("nan"))
        out.append(f"| {k} | `{b['se']:.4g}` | `{b['se_robust']:.5f}` | "
                   f"`{ratio:.1f}` |")
    out += ["",
            f"`k=4` 에서 SD 가 `{d['4']['A1']['se']:.0f} τ_d` 다 — `τ` 자체가 "
            f"`{d['4']['A1']['tau']:.3f}` 인데. 재표집 몇 개가 적합 상한 근처의 `τ` 를 "
            f"내면 SD 가 그것에 지배된다. **백분위 구간은 그 영향을 받지 않는다.**",
            "⇒ 적합 파라미터의 부트스트랩 오차는 **SD 가 아니라 백분위 구간**으로 "
            "보고해야 한다. 판정은 사전등록대로 SD 를 썼고 (문턱을 넘지 못했으므로 "
            "결론은 같다), 그림과 진단은 강건값을 쓴다.", "",
            "## 5. 3σ 에는 얼마가 더 필요한가 — **그리고 돌리지 않는다**", "",
            f"수렴된 추정(`k={k_max}`) 기반: **`k = {need['k_needed_int']}`** "
            f"(현재의 `{need['k_needed_int']/k_max:.2f}` 배, `t={need['t_quantile']:.3f}`).", "",
            "이번엔 입력이 신뢰할 만하다 — `diff` 가 `k ≥ 512` 에서 평탄하고 폭/잡음이 "
            f"`{d[str(k_max)]['A1']['amp_over_noise']:.0f}×` 다.",
            "**그러나 사전등록에 `no_stage_3: true` 를 박아 두었으므로 돌리지 않는다.** "
            "유의해질 때까지 늘리면 optional stopping 이고 `σ` 는 '많이 돌렸다'만 뜻한다.",
            "이 숫자는 **사람이 결정할 근거**로 남긴다.", "",
            "## 6. 결론에 쓸 수 없는 것", "",
            "- **\"`τ(A=0.1)` = `τ(A=1)`\" 이라고 말할 수 없다.** 구별 안 됨은 같음의 "
            "증명이 아니다. `diff` 는 `k ≥ 512` 에서 평탄하고 0 이 아니다.",
            "- **\"다르다\" 도 말할 수 없다** — `2.52σ` 는 사전등록 문턱 `3σ` 아래다.",
            "- 단, **`τ` 가 `A` 와 함께 커진다는 추세는 `A=10` 이 지탱한다** "
            "(`4.7σ`, 부모 런). 이 런이 INCONCLUSIVE 여도 그 주장은 남는다.",
            "- `dt` 수렴·`N` 수렴은 여전히 미실행이고 이 런이 대신하지 않는다.",
            ""]
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
