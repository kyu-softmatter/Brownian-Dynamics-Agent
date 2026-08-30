"""`trap-drag` velocity-sweep re-analysis -- **with seed-ensemble error bars**.

## Why this is being reissued

The previously published v-sweep results (CLAUDE.md: yield force
`F(v->0) = 35-45 kT/d`, "11-20 sigma from 0", shear thinning `F/gamma*v`
16.9 -> 1.87, defect count independent of v) used **a single run per velocity** and
took **block SEM** as the error. But it has since been measured that in a seed
ensemble, block SEM underestimates the spread between realisations -- at
gamma*v=24.21 by **2.35x for dU/particle and 1.37x for F_drag**
(`trap_drag_ensemble.py`). ⚠️ The factor **differs per observable**. The widely
quoted 2.35 is the dU/particle figure; the drag force treated here was 1.37 --
the two must not be interchanged. Either way, the old significances (11-20 sigma)
are overconfident.

Now that all 7 velocities x 9 seeds exist, the numbers are reissued using **the
spread between realisations**. Both errors are plotted side by side, so how the
underestimation factor varies with velocity is also visible.

⚠️ This script **reissues numbers** and does not change any physical
   interpretation. `f_drag_kT_per_d` and the rest are taken from the values already
   in each run's `metrics.json["result"]` -- recomputing them here could drift away
   from the case code.

## The stance on yield force

`F(v->0)` is an **extrapolation**, so it depends on the functional form. Rather than
picking one, all three are reported: the lowest-velocity measurement, a log-linear
extrapolation, and a Herschel-Bulkley fit. If the three disagree, that disagreement
IS the conclusion (the principle: do not assert, show what the answer depends on).

    $PY scratch/trap_drag_vsweep.py
    $PY scratch/trap_drag_vsweep.py --include-legacy   # also include old-code runs (not recommended)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import nondim as ND  # noqa: E402

OUT = ROOT / "runs" / "trap-drag-2d-hex300__ENSEMBLE"


def is_current(p: Path) -> bool:
    """Was this run produced by the current code?

    The discriminator is `step_drift_max_sigma`, present only after the L4
    measurement was wired in.
    """
    try:
        m = json.loads((p / "metrics.json").read_text())
        return "step_drift_max_sigma" in m.get("numerics", {})
    except Exception:
        return False


def collect(include_legacy: bool) -> dict:
    """velocity -> list of runs.

    Grouped by `v_star` (the reduced gamma*v), which is more reliable than the tag
    string.
    """
    by_v: dict[float, list] = {}
    skipped = []
    for p in sorted((ROOT / "runs").glob("trap-drag-2d-hex300__tr0.117647*")):
        if not (p / "metrics.json").exists() or not (p / "spec.json").exists():
            continue
        if "smoke" in p.name:
            continue
        if not include_legacy and not is_current(p):
            skipped.append(p.name)
            continue
        try:
            spec = ND.load(p / "spec.json")
            m = json.loads((p / "metrics.json").read_text())
            r = m["result"]
        except Exception as e:
            skipped.append(f"{p.name} (read failed: {e})")
            continue
        # Exclude variant runs whose warm-up/equilibration/relaxation phase settings
        # differ from the defaults (w2e3r10 and similar)
        tag = p.name.split("__")[1]
        if "-w" in tag:
            skipped.append(f"{p.name} (variant phase settings)")
            continue
        by_v.setdefault(float(spec.params["v_star"]), []).append({
            "name": p.name,
            "seed": int(spec.numerics["seed"]),
            "v_star": float(spec.params["v_star"]),
            "F": float(r["f_drag_kT_per_d"]),
            "F_blocksem": float(r["f_drag_sem"]),
            "stokes": float(r["f_stokes_bare"]),
            "n_def": float(r["n_def"]["driven"]),
            "psi6": float(r["psi6"]["driven"]),
            "rec": float(r["relax_fit_defects"]["recovered_frac"]),
            "dev": float(r["dev_max_d"]),
        })
    return by_v, skipped


def stats(xs):
    a = np.asarray([x for x in xs if np.isfinite(x)], float)
    if a.size == 0:
        return dict(mean=np.nan, std=np.nan, sem=np.nan, n=0)
    return dict(mean=a.mean(), std=a.std(ddof=1) if a.size > 1 else np.nan,
                sem=(a.std(ddof=1) / np.sqrt(a.size)) if a.size > 1 else np.nan, n=a.size)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-legacy", action="store_true")
    a = ap.parse_args()

    by_v, skipped = collect(a.include_legacy)
    if skipped:
        print(f"⚠ excluded {len(skipped)} run(s) (old code or variant settings) -- "
              f"never dropped silently:")
        for s in skipped[:10]:
            print(f"    - {s}")
        if len(skipped) > 10:
            print(f"    ... and {len(skipped)-10} more")
    if len(by_v) < 3:
        print(f"not enough velocity points ({len(by_v)})", file=sys.stderr)
        return 1

    vs = sorted(by_v)
    print("=" * 104)
    print(f"trap-drag v-sweep re-analysis -- {len(vs)} velocities, seed ensemble")
    print("=" * 104)
    print(f"{'gam.v':>8}{'seeds':>5}{'<F>':>9}{'ensSEM':>11}{'blkSEM':>10}"
          f"{'underest':>9}{'F/gam.v':>8}{'defect':>7}{'psi6':>7}{'recov%':>7}")
    rows = []
    for v in vs:
        g = by_v[v]
        F = stats([x["F"] for x in g])
        blk = float(np.mean([x["F_blocksem"] for x in g]))
        ratio = F["std"] / blk if blk else np.nan
        nd = stats([x["n_def"] for x in g])
        ps = stats([x["psi6"] for x in g])
        rc = stats([x["rec"] for x in g])
        rows.append(dict(v=v, F=F, blk=blk, ratio=ratio, nd=nd, psi6=ps, rec=rc,
                         thin=F["mean"] / v))
        print(f"{v:>8.2f}{F['n']:>5}{F['mean']:>9.1f}{F['sem']:>11.1f}{blk:>10.1f}"
              f"{ratio:>9.2f}{F['mean']/v:>8.2f}{nd['mean']:>7.2f}{ps['mean']:>7.3f}"
              f"{100*rc['mean']:>7.1f}")

    # ── block-SEM underestimation factor ───────────────────────────────
    rr = np.array([r["ratio"] for r in rows], float)
    rr = rr[np.isfinite(rr)]
    print("-" * 104)
    # ⚠️ Compare per observable, exactly: the single-v ensemble (gamma*v=24.21) gave
    #    underestimation factors of **2.35x for dU/particle and 1.37x for F_drag**.
    #    What is measured here is F_drag, so it must be compared against 1.37.
    #    Comparing it against 2.35 would be mixing two different observables.
    print(f"block-SEM underestimation factor (F_drag): median {np.median(rr):.2f}x, "
          f"range {rr.min():.2f}-{rr.max():.2f}x")
    print(f"  reference -- single-v ensemble (gamma*v=24.21): F_drag 1.37x, "
          f"dU/particle 2.35x (it differs per observable)")

    # ── yield force -- three extrapolations ────────────────────────────
    V = np.array([r["v"] for r in rows])
    FM = np.array([r["F"]["mean"] for r in rows])
    FE = np.array([r["F"]["sem"] for r in rows])
    print()
    print("yield force F(v->0) -- **an extrapolation, so it depends on the "
          "functional form. All three are reported**")
    lo = int(np.argmin(V))
    print(f"  (a) lowest-v measured   at gamma*v={V[lo]:.2f}, "
          f"F = {FM[lo]:.1f} ± {FE[lo]:.1f} kT/d"
          f"   -> {FM[lo]/FE[lo]:.1f} sigma from 0   "
          f"(no extrapolation, most defensible)")

    # (b) log-linear: F = a + b*ln(gamma*v)  -- the three lowest velocities
    k = min(3, len(V))
    idx = np.argsort(V)[:k]
    A = np.vstack([np.ones(k), np.log(V[idx])]).T
    w = 1.0 / np.where(FE[idx] > 0, FE[idx], np.nan)
    good = np.isfinite(w)
    if good.sum() >= 2:
        coef, *_ = np.linalg.lstsq(A[good] * w[good, None], FM[idx][good] * w[good], rcond=None)
        # v->0 sends ln -> -inf and diverges; the conclusion is that this form
        # cannot define a yield force at all
        print(f"  (b) log-linear fit      "
              f"F = {coef[0]:.1f} + {coef[1]:.1f}*ln(gamma*v)  "
              f"-> **diverges** as v->0 (undefined). A log form cannot give a "
              f"yield force")

    # (c) Herschel-Bulkley: F = F_y + c*(gamma*v)^n -- whole range, n scanned on a grid
    best = None
    for n_ in np.linspace(0.1, 1.5, 141):
        X = np.vstack([np.ones(len(V)), V ** n_]).T
        ww = 1.0 / np.where(FE > 0, FE, np.nan)
        m_ = np.isfinite(ww)
        if m_.sum() < 3:
            continue
        c, res, *_ = np.linalg.lstsq(X[m_] * ww[m_, None], FM[m_] * ww[m_], rcond=None)
        pred = X @ c
        chi2 = float(np.nansum(((FM - pred) * ww) ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, n_, c)
    if best:
        chi2, n_, c = best
        dof = max(1, int(np.isfinite(FE).sum()) - 3)
        # Uncertainty on F_y: the weighted least-squares covariance
        X = np.vstack([np.ones(len(V)), V ** n_]).T
        ww = 1.0 / np.where(FE > 0, FE, np.nan)
        m_ = np.isfinite(ww)
        XtX = (X[m_] * ww[m_, None]).T @ (X[m_] * ww[m_, None])
        cov = np.linalg.inv(XtX) * max(1.0, chi2 / dof)
        sF = float(np.sqrt(cov[0, 0]))
        print(f"  (c) Herschel-Bulkley    F = F_y + c*(gamma*v)^n,  "
              f"n={n_:.2f} optimal "
              f"-> F_y = {c[0]:.1f} ± {sF:.1f} kT/d   "
              f"-> {abs(c[0])/sF:.1f} sigma from 0"
              f"   (χ²/dof = {chi2/dof:.2f})")

    print()
    print("shear thinning F/gamma*v")
    print(f"  γv={V[lo]:.2f}: {FM[lo]/V[lo]:.2f}  →  γv={V[-1]:.2f}: {FM[-1]/V[-1]:.2f}"
          f"   ({(FM[lo]/V[lo])/(FM[-1]/V[-1]):.1f}x reduction)")

    nd_m = np.array([r["nd"]["mean"] for r in rows])
    nd_e = np.array([r["nd"]["sem"] for r in rows])
    print()
    print("defect count -- does it depend on v?")
    print(f"  range {nd_m.min():.2f}-{nd_m.max():.2f}, "
          f"mean ensemble SEM {np.nanmean(nd_e):.2f}")
    spread = (nd_m.max() - nd_m.min()) / np.nanmean(nd_e)
    print(f"  max-min = {spread:.1f}*SEM -> "
          + ("**v-dependent**" if spread > 3
             else "indistinguishable from v-independent"))

    _plot(rows, V, FM, FE, best)
    _save(rows, best)
    return 0


def _plot(rows, V, FM, FE, best) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]     # * labels in English (CLAUDE.md)
    matplotlib.rcParams["axes.unicode_minus"] = False

    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    blk = np.array([r["blk"] for r in rows])

    a = ax[0, 0]
    a.errorbar(V, FM, yerr=FE, fmt="o-", capsize=4, lw=1.8, color="tab:red",
               label="ensemble SEM (between seeds)")
    a.errorbar(V, FM, yerr=blk, fmt="none", capsize=8, lw=1.0, color="tab:blue",
               alpha=.8, label="block SEM (single run) -- underestimates")
    a.plot(V, V, "k--", lw=1.2, label=r"bare Stokes $\gamma v$")
    if best:
        _, n_, c = best
        vv = np.geomspace(V.min() * .5, V.max() * 1.2, 100)
        a.plot(vv, c[0] + c[1] * vv ** n_, ":", color="darkgreen", lw=1.6,
               label=f"HB fit (n={n_:.2f}, $F_y$={c[0]:.0f})")
        a.axhline(c[0], color="darkgreen", ls="-.", lw=1.0, alpha=.6)
    a.set(xscale="log", xlabel=r"$\gamma v$ [kT/d]", ylabel=r"$\langle F_x\rangle$ [kT/d]",
          title="(1) drag force vs velocity -- the two errors compared")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[0, 1]
    a.plot(V, FM / V, "o-", lw=1.8, color="tab:purple")
    a.axhline(1, color="k", ls="--", lw=1.2, label="bare Stokes (=1)")
    a.set(xscale="log", yscale="log", xlabel=r"$\gamma v$ [kT/d]",
          ylabel=r"$F/\gamma v$",
          title="(2) shear thinning -- effective drag multiplier")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[1, 0]
    nd = np.array([r["nd"]["mean"] for r in rows])
    nde = np.array([r["nd"]["sem"] for r in rows])
    a.errorbar(V, nd, yerr=nde, fmt="s-", capsize=4, lw=1.8, color="tab:orange")
    a.set(xscale="log", xlabel=r"$\gamma v$ [kT/d]",
          ylabel="defect count while driven",
          title="(3) defect count -- v dependence")
    a.grid(alpha=.3, which="both")

    a = ax[1, 1]
    ratio = np.array([r["ratio"] for r in rows])
    a.plot(V, ratio, "D-", lw=1.8, color="tab:green")
    a.axhline(1, color="k", ls="--", lw=1.2, label="agreement (=1)")
    a.axhline(1.37, color="tab:red", ls=":", lw=1.4,
              label=r"single-v measurement, $F_{drag}$ 1.37x")
    a.set(xscale="log", xlabel=r"$\gamma v$ [kT/d]",
          ylabel="std across realisations / block SEM",
          title="(4) block-SEM underestimation factor")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    fig.tight_layout()
    p = OUT / "vsweep_ensemble.png"
    fig.savefig(p, dpi=130)
    print(f"\nfigure: {p}")


def _save(rows, best) -> None:
    out = {"schema": "bdbot.trap_drag_vsweep/0.2",
           "note": "Seed-ensemble error bars. The error is the spread between "
                   "realisations (ensemble SEM), NOT block SEM.",
           "points": [{"v_star": r["v"], "n_seeds": r["F"]["n"],
                       "F_mean": r["F"]["mean"], "F_ensemble_sem": r["F"]["sem"],
                       "F_realization_std": r["F"]["std"], "F_block_sem_mean": r["blk"],
                       "block_sem_underestimate": r["ratio"],
                       "F_over_gamma_v": r["thin"],
                       "n_def_mean": r["nd"]["mean"], "n_def_sem": r["nd"]["sem"],
                       "psi6_mean": r["psi6"]["mean"],
                       "recovered_frac_mean": r["rec"]["mean"]} for r in rows]}
    if best:
        chi2, n_, c = best
        out["herschel_bulkley"] = {"n": n_, "F_yield": c[0], "c": c[1]}
    p = OUT / "vsweep_ensemble.json"
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"data: {p}")


if __name__ == "__main__":
    sys.exit(main())
