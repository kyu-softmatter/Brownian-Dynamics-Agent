"""soft-r3 `A` sweep analysis -- seed-ensemble error bars plus a comparison against
Zahn's phase diagram.

usage: python scripts/soft2d_sweep_analyze.py <sweep_dir> [--json out.json]

## What this script is careful about

**(1) The seed ensemble is the honest error estimate.** Frame-to-frame spread includes
   temporal correlation, so it is not a statistical error. `aggregate_seeds` is used
   to take the SE of the seed means.

**(2) It does not report the defect fraction alone.** `A=100` from a random start
   gives `{5:2%, 6:96%, 7:2%}` -- equal numbers of 5s and 7s, the signature of
   **dislocations** -- while `A=0.1` spreads over 3-12, a **liquid**. The same
   fraction can mean entirely different physics, so the coordination distribution and
   the `5-7 symmetry` are reported alongside it.

**(3) It does not stand in for an equilibration verdict.** It compares `psi6` over the
   first and second halves and **only reports whether there is drift**. Saying
   "equilibrium was reached" would require a threshold.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.analysis.structure import (hex_order, rdf, structure_factor,
                                       zahn_phase)
from simbot.analysis.trap import aggregate_seeds


def load(d: Path) -> dict:
    man = json.loads((d / "manifest.json").read_text())
    z = np.load(d / "samples.npz")
    return {"man": man, "traj": z["traj"], "energy": z["energy"],
            "Lx": float(z["box"][0]), "Ly": float(z["box"][1])}


def per_run(r: dict) -> dict:
    """Measurements for one run. Uses the second half only, and also reports the
    difference against the first half."""
    traj, Lx, Ly = r["traj"], r["Lx"], r["Ly"]
    n = len(traj)
    first, second = traj[: n // 2], traj[n // 2:]
    h2 = hex_order(second, Lx=Lx, Ly=Ly)
    h1 = hex_order(first, Lx=Lx, Ly=Ly)
    s = structure_factor(second, Lx=Lx, Ly=Ly, n_max=18)
    g = rdf(second, Lx=Lx, Ly=Ly, bins=200)
    hist = h2.coordination_hist
    n5, n7 = hist.get(5, 0.0), hist.get(7, 0.0)
    return {
        "psi6_global": h2.psi6_global, "psi6_local": h2.psi6_local_mean,
        "defect_fraction": h2.defect_fraction,
        "sixfold": s.sixfold_modulation, "first_peak_g": g.first_peak_g,
        "first_peak_r": g.first_peak_r,
        "energy_pp": float(r["energy"][n // 2:].mean()) / r["man"]["config"]["n_particles"],
        # ★ dislocation signature: are the 5s and 7s equal in number?
        "coord_5": n5, "coord_7": n7,
        "five_seven_balance": (abs(n5 - n7) / max(n5 + n7, 1e-12)),
        "coord_spread": float(len([k for k, v in hist.items() if v > 0.005])),
        "coordination_hist": {int(k): float(v) for k, v in sorted(hist.items())},
        # drift (an equilibration diagnostic -- NOT a verdict)
        "psi6_drift": h2.psi6_global - h1.psi6_global,
        "min_separation": r["man"]["guards"]["min_separation"],
        "min_sep_over_r_min": r["man"]["guards"]["min_separation_over_r_min"],
        "force_displacement": r["man"]["guards"]["force_displacement_star"],
        "wall_s": r["man"]["wall_s"],
    }


def main() -> int:
    root = Path(sys.argv[1])
    dirs = sorted(p for p in root.glob("A*_*_s*") if (p / "samples.npz").exists())
    if not dirs:
        print(f"⛔ no runs under {root}", file=sys.stderr)
        return 2

    groups: dict[tuple[float, str], list[dict]] = {}
    for d in dirs:
        cfg = json.loads((d / "manifest.json").read_text())["config"]
        key = (float(cfg["amplitude"]), cfg["init"])
        groups.setdefault(key, []).append(per_run(load(d)))

    print(f"# soft-r3 A sweep -- {len(dirs)} runs / {len(groups)} conditions\n")
    metrics = ("psi6_global", "psi6_local", "defect_fraction", "sixfold",
               "energy_pp", "first_peak_g")
    out: dict = {}

    hdr = f"{'A':>6} {'init':<7} {'n':>2} " + " ".join(
        f"{m[:11]:>17}" for m in metrics)
    print(hdr)
    print("-" * len(hdr))
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        runs = groups[key]
        agg = {m: aggregate_seeds([r[m] for r in runs]) for m in metrics}
        out[f"A{A:g}_{init}"] = {
            "amplitude": A, "init": init, "n_seeds": len(runs),
            "zahn": zahn_phase(A),
            **{m: {"mean": agg[m].mean, "se": agg[m].se, "spread": agg[m].spread,
                   "values": agg[m].values} for m in metrics},
            "coordination_hist_seed1": runs[0]["coordination_hist"],
            "five_seven_balance": [r["five_seven_balance"] for r in runs],
            "psi6_drift": [r["psi6_drift"] for r in runs],
            "min_sep_over_r_min": min(r["min_sep_over_r_min"] for r in runs),
            "wall_s_total": sum(r["wall_s"] for r in runs),
        }
        cells = " ".join(f"{agg[m].mean:>9.4f}±{agg[m].se:<7.4f}" for m in metrics)
        print(f"{A:>6g} {init:<7} {len(runs):>2} {cells}")

    # --- comparison against Zahn's phase diagram ---
    print("\n## comparison against Zahn's phase diagram  [cited, not reproduced]")
    print(f"{'A':>6} {'Gamma':>9} {'Zahn predicts':<18} {'init':<7} "
          f"{'psi6 global':>17} {'6-fold modulation':>17} {'our reading':<12}")
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        o = out[f"A{A:g}_{init}"]
        p6, s6 = o["psi6_global"], o["sixfold"]
        # The reading is **observation-based** only; Zahn's boundaries are not used
        if p6["mean"] > 0.7 and s6["mean"] > 0.7:
            verdict = "crystal-like"
        elif s6["mean"] > 0.3:
            verdict = "6-fold residual"
        else:
            verdict = "isotropic-like"
        print(f"{A:>6g} {o['zahn']['gamma']:>9.2f} {o['zahn']['phase_zahn']:<18} "
              f"{init:<7} {p6['mean']:>9.4f}±{p6['se']:<7.4f} "
              f"{s6['mean']:>9.4f}±{s6['se']:<7.4f} {verdict:<12}")

    # --- initial-condition dependence ---
    print("\n## initial-condition dependence (hex vs random) -- against the seed "
          "ensemble")
    print(f"{'A':>6} {'psi6 hex':>17} {'psi6 random':>17} {'diff':>10} "
          f"{'combined SE':>9} {'sigma':>7}  verdict")
    for A in sorted({k[0] for k in groups}):
        h, r = out.get(f"A{A:g}_hex"), out.get(f"A{A:g}_random")
        if not (h and r):
            continue
        dh, dr = h["psi6_global"], r["psi6_global"]
        diff = dh["mean"] - dr["mean"]
        se = float(np.hypot(dh["se"] or 0.0, dr["se"] or 0.0))
        sig = abs(diff) / se if se > 0 else float("inf")
        # ★ within 3 sigma is "indistinguishable" -- NOT proof that they are equal
        mark = "indistinguishable" if sig < 3 else "❗significant difference"
        print(f"{A:>6g} {dh['mean']:>9.4f}±{dh['se']:<7.4f} "
              f"{dr['mean']:>9.4f}±{dr['se']:<7.4f} {diff:>+10.4f} {se:>9.4f} "
              f"{sig:>6.2f}σ  {mark}")

    # --- the character of the defects ---
    print("\n## the character of the defects -- the fraction alone cannot "
          "distinguish them")
    print(f"{'A':>6} {'init':<7} {'defect fraction':>17} {'5-7 imbalance':>11} "
          f"{'coord kinds':>10}  coordination distribution (seed 1)")
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        o = out[f"A{A:g}_{init}"]
        df = o["defect_fraction"]
        bal = float(np.mean(o["five_seven_balance"]))
        hist = o["coordination_hist_seed1"]
        kinds = len([k for k, v in hist.items() if v > 0.005])
        top = " ".join(f"{k}:{v:.2f}" for k, v in hist.items() if v > 0.005)
        note = "dislocations (5-7 pairs)" if bal < 0.25 and df["mean"] < 0.2 else ""
        print(f"{A:>6g} {init:<7} {df['mean']:>9.4f}±{df['se']:<7.4f} "
              f"{bal:>11.3f} {kinds:>10d}  {top}  {note}")

    # --- diagnostics (not verdicts) ---
    print("\n## diagnostics -- reporting facts, not verdicts")
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        o = out[f"A{A:g}_{init}"]
        drift = np.array(o["psi6_drift"])
        print(f"  A={A:<6g} {init:<7} psi6 drift (2nd-1st) = "
              f"{drift.mean():+.4f} ± {drift.std(ddof=1)/np.sqrt(drift.size):.4f}"
              f"   min_sep/r_min = {o['min_sep_over_r_min']:.2f}"
              f"   wall {o['wall_s_total']:.0f}s")
    print("\n  ⚠ Drift being indistinguishable from 0 does NOT mean 'equilibrium was")
    print("    reached' -- it means this run length cannot resolve that drift.")

    if "--json" in sys.argv:
        p = Path(sys.argv[sys.argv.index("--json") + 1])
        p.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=float))
        print(f"\n→ {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
