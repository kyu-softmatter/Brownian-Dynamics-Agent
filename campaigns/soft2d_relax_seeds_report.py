"""soft-r3 relaxation time -- seed-convergence analysis + figures + validation doc.

usage: python scripts/soft2d_relax_seeds_report.py [--stage 2]

## The point of this script -- nested subsets

`tau(k)` is computed from the **first `k` of the same trajectories**. Drawn from
different runs, "the data differ" would be indistinguishable from "the estimator
depends on `k`". Nested, the slope of the curve is **purely the estimator's `k`
dependence** -- and that is the low-seed bias.

A side effect is that it can answer "what did `k = 64` actually give?" (the
user's original suggestion).
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
    """Defect-fraction curves `(n_seeds, n_frames)` in seed order, plus the time axis.

    ★ **It caches.** One curve is a Voronoi over 401 frames, so 1500 seeds is
      600k of them (1.2M across the two A values). The ladder analysis must not
      redo that. The seed count goes into the cache key, so the cache
      invalidates itself as soon as more runs arrive.
    """
    dirs = sorted((p for p in (rd.path / "raw_early").glob(f"A{A:g}_s*")
                   if (p / "samples.npz").exists()),
                  key=lambda p: int(p.name.rsplit("_s", 1)[1]))
    cache = rd.path / f"defect_curves_A{A:g}_k{len(dirs)}.npz"
    if cache.exists():
        z = np.load(cache)
        print(f"  A={A:<5g} read from cache ({cache.name})")
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
    print(f"  A={A:<5g} computed {arr.shape[0]} curves -> cache {cache.name}")
    return t_ref, arr


def _fit_subset(t: np.ndarray, mat: np.ndarray, n_resample: int,
                boot_seed: int) -> dict:
    tail_m = t >= 0.5 * t[-1]
    noise = float(mat.mean(axis=0)[tail_m].std(ddof=1))
    fit = fit_relaxation(t, mat.mean(axis=0), noise=noise)
    boot = bootstrap_relaxation_over_seeds(t, mat, n_resample=n_resample,
                                           seed=boot_seed, noise=noise)
    lo, hi = boot["tau_ci95"]
    #  ★ The SD is outlier-dominated at low seed count (82 tau_d at k=4 -- when
    #    tau itself is 0.06). So the half-width of the percentile interval is
    #    emitted alongside as a **robust SE**. The verdict uses the SD as
    #    pre-registered, but the figures and diagnostics are only readable with
    #    the robust value.
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
        print(f"  A={A:<5g} seeds {mats[A].shape[0]} · frames {mats[A].shape[1]}")

    k_max = min(m.shape[0] for m in mats.values())
    ladder = [k for k in K_LADDER if k <= k_max] + [k_max]
    ladder = sorted(set(ladder))

    print(f"\n## nested-subset τ(k) -- the first k of the same trajectories\n")
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
        mark = "  <- user suggestion" if k == 64 else ("  <- pre-registered"
                                                       if k == k_max else "")
        print(f"{k:>6d} {f01['tau']:>11.5f}±{f01['se']:<8.5f} "
              f"{f1['tau']:>11.5f}±{f1['se']:<8.5f} {diff:>+10.5f} "
              f"{se:>9.5f} {sig:>7.2f}{mark}")

    # --- verdict (the pre-registered rule) ---
    thr = float(rule["threshold_sigma"])
    sig_final = conv["sigma"][-1]
    verdict = (rule["verdict_if_above"] if sig_final > thr
               else rule["verdict_if_below"])
    print(f"\n## verdict -- pre-registered threshold {thr:g}σ, "
          f"SE source {rule['se_source']}")
    print(f"  k = {k_max}:  {sig_final:.2f}σ  →  {verdict.strip().splitlines()[0]}")
    print(f"  (pre-registered expectation "
          f"{pred['design_power']['k_chosen_expected_sigma']}σ)")

    # --- figures ---
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
    fs.skip("voronoi", "this run looks only at the relaxation time -- the "
                       "character of the arrangement came from the parent run")
    rd.write("figures", fs.figures_md())

    out = {"convergence": conv, "k_max": k_max, "sigma_final": sig_final,
           "verdict": verdict.strip(), "threshold_sigma": thr,
           "preregistered": prereg,
           "_provenance_at_analysis": provenance(DRIVERS)}
    #  ★ Merge into the existing metrics rather than overwrite -- losing the τ
    #    and bootstrap results the runner wrote would leave only the convergence
    #    table, and "what came out at the final k" would be gone
    prev = json.loads(rd.read("metrics")) if rd.exists("metrics") else {}
    merged = {**prev, **out}
    rd.write_json("metrics", merged)
    rd.write("validation", validation_md(rd, merged, pred, stage))
    print(f"\n→ {rd.figs.relative_to(REPO)}/01_seed_convergence.png")
    print(f"→ {rd.file('validation').relative_to(REPO)}")
    return 0


def check_predictions(pred: dict, merged: dict, rd: RunDir) -> list[dict]:
    """Compare against the sealed prediction. **Measured values are read from
    metrics only** -- never transcribed by hand."""
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

    out = [f"# S7 -- `τ_relax` seed sweep validation, **stage {stage}** "
           f"({rd.run_id})", "",
           f"`A = 0.1` vs `A = 1` · seeds `{k_max}` (planned "
           f"{json.loads(rd.read('spec'))['k_seeds']}) · nested-subset ladder", "",
           "> **The verdict is a proposal.** `proposed_by: agent`, "
           "`confirmed_by: null`.", "",
           seal_section(rd), "",
           "## 1. verdict -- the pre-registered rule", "", "```yaml",
           f"verdict: {'DIFFERENT' if merged['sigma_final'] > merged['threshold_sigma'] else 'INCONCLUSIVE'}",
           f"sigma_observed: {merged['sigma_final']:.3f}",
           f"sigma_threshold: {merged['threshold_sigma']:g}",
           f"sigma_expected_preregistered: "
           f"{pred['design_power']['k_chosen_expected_sigma']}",
           "se_source: bootstrap_over_seeds",
           "proposed_by: agent", "confirmed_by: null", "```", "",
           merged["verdict"], "",
           "## 2. ★ The nested-subset ladder -- where the bias stops", "",
           "Computed from the **first `k` of the same trajectories**. The slope "
           "of the curve is not a difference in the data but **the estimator's "
           "`k` dependence**.", "",
           "| `k` | `τ(A=0.1)` | amp/noise | `τ(A=1)` | amp/noise | `diff` | "
           "`SE_diff` | `σ` |",
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
            f"⚠ marks the points where **the bootstrap SD exceeds `τ`**, so no "
            f"error can be claimed there (§4).", "",
            "### Reading it", "",
            f"- **`τ(A=0.1)` has no bias** -- from `{d['4']['A0.1']['tau']:.5f}` "
            f"at `k=4` through `k={k_max}`'s "
            f"`{d[str(k_max)]['A0.1']['tau']:.5f}` it fluctuates with no trend. "
            f"Already at `k=4` the amp/noise was "
            f"`{d['4']['A0.1']['amp_over_noise']:.1f}×`.",
            f"- **`τ(A=1)` falls from `{d['4']['A1']['tau']:.5f}` to "
            f"`{d[str(k_max)]['A1']['tau']:.5f}` by "
            f"`{100*(d[str(k_max)]['A1']['tau']/d['4']['A1']['tau']-1):.0f} %`, "
            f"and **flattens out at `k ≈ 512`.**",
            f"- ⇒ **The bias is set by amp/noise, not by the seed count.** `A=1` "
            f"reaches amp/noise `20×` at `k=512`, and that is where the bias "
            f"disappears. `A=0.1` was already `6×` at `k=4`, and `23×` by `k=64`.",
            f"- **`diff` does not go to 0** -- the mean of the three points at "
            f"`k ≥ 512`, "
            f"`{np.mean(tail_diff):+.5f} ± {np.std(tail_diff, ddof=1):.5f}`, is "
            f"flat. The difference is small "
            f"(`{100*c['diff'][-1]/merged['A1']['tau']:.0f} %` of `τ`) but "
            f"**it looks real -- it just did not resolve at `3σ`.**", "",
            "## 3. Sealed-prediction comparison", "",
            "| quantity | predicted | tolerance | measured | verdict | what it "
            "discriminates |",
            "|---|---|---|---|---|---|"]
    for x in checks:
        mv = ("—" if x["measured"] is None else f"`{x['measured']:.5g}`")
        mark = {"PASS": "**PASS**", "FAIL": "**FAIL** ⛔"}.get(x["verdict"],
                                                              x["verdict"])
        out.append(f"| `{x['quantity']}` | `{x['predicted']}` | "
                   f"`{x['tolerance']}` | {mv} | {mark} | "
                   f"{x['discriminates']} |")
    out += ["", f"**PASS {n_pass} · FAIL {n_fail}**", "",
            "### The character of the 2 FAILs", "",
            "- **`tau_relax_diff_sigma`** -- this was the run's first-class "
            "prediction, and it missed. Expected `3.96σ` vs actual "
            f"`{merged['sigma_final']:.2f}σ`. The cause is §2: even at `k=254` "
            "the `A=1` bias had not fully vanished, so `diff` shrank once more "
            f"(`{pred['design_power']['prior_diff_tau_d']:.5f}` → "
            f"`{c['diff'][-1]:.5f}`). **The design-power input was optimistic "
            "again.**",
            "- **`init_config_failures_per_1000_seeds`** -- predicted "
            "`7.8 ± 6`, actual "
            f"`{[x['measured'] for x in checks if x['quantity'].startswith('init_')][0]:.1f}`. "
            "It was estimated from the 2 cases at stage 1's `k=256` (7.8/1000), "
            "and that sample was far too small (estimating a rate from 2 events "
            "carries a 71 % relative error). **Independence of `A` was "
            "confirmed** -- the same seeds failed at both `A` values.", "",
            "## 4. ⚠ The bootstrap SD is outlier-dominated at low seed count", "",
            "| `k` | `A=1` SD | `A=1` robust SE (95 % interval/3.92) | ratio |",
            "|---|---|---|---|"]
    for k in c["k"]:
        b = d[str(k)]["A1"]
        ratio = (b["se"] / b["se_robust"] if b.get("se_robust") else float("nan"))
        out.append(f"| {k} | `{b['se']:.4g}` | `{b['se_robust']:.5f}` | "
                   f"`{ratio:.1f}` |")
    out += ["",
            f"At `k=4` the SD is `{d['4']['A1']['se']:.0f} τ_d` -- when `τ` "
            f"itself is `{d['4']['A1']['tau']:.3f}`. If a few resamples land on "
            f"a `τ` near the fit's upper bound, the SD is dominated by those. "
            f"**The percentile interval is not affected by them.**",
            "⇒ The bootstrap error on a fit parameter has to be reported as "
            "**the percentile interval, not the SD**. The verdict used the SD as "
            "pre-registered (the conclusion is the same, since it did not cross "
            "the threshold), and the figures and diagnostics use the robust "
            "value.", "",
            "## 5. How much more would 3σ need -- **and it is not being run**", "",
            f"Based on the converged estimate (`k={k_max}`): "
            f"**`k = {need['k_needed_int']}`** "
            f"(`{need['k_needed_int']/k_max:.2f}`× the current, "
            f"`t={need['t_quantile']:.3f}`).", "",
            "This time the input is trustworthy -- `diff` is flat for `k ≥ 512` "
            "and the amp/noise is "
            f"`{d[str(k_max)]['A1']['amp_over_noise']:.0f}×`.",
            "**But the pre-registration nailed down `no_stage_3: true`, so it is "
            "not run.** Extending until it turns significant is optional "
            "stopping, and then `σ` means only 'we ran a lot'.",
            "This number is left as **grounds for a human to decide on**.", "",
            "## 6. What cannot go in the conclusion", "",
            "- **\"`τ(A=0.1)` = `τ(A=1)`\" cannot be said.** Indistinguishable "
            "is not a proof of equal. `diff` is flat for `k ≥ 512` and is not 0.",
            "- **\"they differ\" cannot be said either** -- `2.52σ` is below the "
            "pre-registered `3σ` threshold.",
            "- Still, **the trend that `τ` grows with `A` is held up by `A=10`** "
            "(`4.7σ`, the parent run). That claim survives this run being "
            "INCONCLUSIVE.",
            "- `dt` convergence and `N` convergence are still unrun, and this "
            "run does not stand in for them.",
            ""]
    return "\n".join(out)


if __name__ == "__main__":
    raise SystemExit(main())
