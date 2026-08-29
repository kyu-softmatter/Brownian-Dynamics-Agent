"""`trap-drag` ensemble average -- collect several runs differing only in seed and
look at the spread between realisations.

**Why**: "dU = 4.7*n_def", "defects vanish at 24.7 tau_int" and "the relaxation is a
plateau plus a sharp drop" all came from a SINGLE realisation. Dislocation creation
is stochastic, so the timing of the drop will vary from seed to seed, and if it
does, **the drop gets smeared out in the ensemble average** and a smooth curve
results. Distinguishing those two possibilities is the purpose of this script.

Checked alongside it: does a single run's block SEM agree with the **standard
deviation between realisations**? Block SEM only corrects for correlation within one
trajectory, so if it misses a 'slow degree of freedom' -- such as a different
dislocation arrangement in each realisation -- it underestimates the error. Only an
ensemble can confirm that.

    $PY scratch/trap_drag_ensemble.py [--glob 'trap-drag-2d-hex300__tr0.117647-s*']
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from bdbot import nondim as ND, stats as ST  # noqa: E402

# ★ These MUST match cases/trap_drag_2d.py exactly -- they filter the stored `phase`
#   column, and a mismatch silently selects nothing rather than raising.
PH_EQ, PH_DRAG, PH_RELAX = "equilibrium", "drag", "relaxation"


def load(d: Path):
    spec = ND.load(d / "spec.json")
    z = np.load(d / "observables.npz")
    r = json.loads((d / "metrics.json").read_text())["result"]
    ti = spec.reduced("times", "tau_int")
    t = z["_t_step"] * float(spec.numerics["dt_star"]) / ti
    return spec, z, r, t, z["phase"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="trap-drag-2d-hex300__tr0.117647-s*")
    ap.add_argument("--out", default="runs/trap-drag-2d-hex300__ENSEMBLE")
    ap.add_argument("--include-legacy", action="store_true",
                    help="also include runs with no step-resolution measurement "
                         "(pre-2026-08-05); excluded by default")
    args = ap.parse_args()

    dirs = sorted(p for p in (ROOT / "runs").glob(args.glob)
                  if (p / "observables.npz").exists() and (p / "metrics.json").exists())

    # ★ **Legacy runs with the same tag but a different hash match the same glob.**
    #   Both `tr0.117647-s1__<old hash>` and `tr0.117647-s1__<new hash>` match, so
    #   left alone this **mixes old-code and new-code runs into one average and
    #   counts the seed twice.**
    #   The discriminator is the presence of `step_drift_max_sigma` -- that is
    #   exactly "was run by the current code".
    #   ⚠️ Never drop them silently; always print what was excluded.
    def is_current(p: Path) -> bool:
        try:
            m = json.loads((p / "metrics.json").read_text())
            return "step_drift_max_sigma" in m.get("numerics", {})
        except Exception:
            return False

    if not args.include_legacy:
        legacy = [p for p in dirs if not is_current(p)]
        dirs = [p for p in dirs if is_current(p)]
        if legacy:
            print(f"⚠ excluded {len(legacy)} legacy run(s) (no step-resolution "
                  f"measurement = old code). To include them, --include-legacy:")
            for p in legacy:
                print(f"    - {p.name}")
    if len(dirs) < 2:
        print(f"not enough runs ({len(dirs)})", file=sys.stderr)
        return 1
    print(f"ensemble over {len(dirs)} runs: "
          f"{', '.join(p.name.split('__')[1] for p in dirs)}")

    rows, rel_u, rel_n, drag_f = [], [], [], []
    t_rel = None
    for d in dirs:
        spec, z, r, t, ph = load(d)
        N = spec.params["N"]
        rel = ph == PH_RELAX
        tr = t[rel] - t[rel][0]
        if t_rel is None or len(tr) < len(t_rel):
            t_rel = tr
        rel_u.append((z["u_pair"][rel] - r["U_pair_eq"]) * N)     # ΔU_total(t)
        rel_n.append(z["n_def"][rel].astype(float))
        nd = z["n_def"][rel]
        t0 = tr[np.argmax(nd == 0)] if (nd == 0).any() else np.nan
        dmask = ph == PH_DRAG
        idx = np.flatnonzero(dmask)
        ss = idx[len(idx) // 2:]                                   # steady state = the second half
        fx = -spec.params["k_star"] * z["dx_probe"][ss]
        drag_f.append(fx)
        # slope of dU vs n_def (drag + relaxation)
        m = dmask | rel
        x = z["n_def"][m].astype(float)
        y = (z["u_pair"][m] - r["U_pair_eq"]) * N
        slope = np.linalg.lstsq(np.vstack([x, np.ones_like(x)]).T, y, rcond=None)[0][0]
        rows.append(dict(seed=int(spec.numerics["seed"]), dU=r["dU_drive"],
                         dU_sem=r["dU_sem"], F=float(fx.mean()),
                         F_sem=float(ST.block_sem(fx)),
                         nd_drag=float(z["n_def"][ss].mean()),
                         nd_max=float(z["n_def"][dmask].max()),
                         t_zero=float(t0), slope=float(slope),
                         psi6=r["psi6"]["driven"]))

    n = min(len(a) for a in rel_u)
    U = np.array([a[:n] for a in rel_u])         # (M, n)
    D = np.array([a[:n] for a in rel_n])
    t_rel = t_rel[:n]
    M = len(dirs)

    print("=" * 96)
    print(f"ensemble -- {M} seeds, 2 lattice periods (3.22 d), "
          f"N={rows[0].get('N', 306)}")
    print("=" * 96)
    print(f"{'seed':>5}{'dU/part':>10}{'(runSEM)':>10}{'F_drag':>9}{'(runSEM)':>10}"
          f"{'def drag':>10}{'max':>6}{'t(def=0)':>11}{'slope':>8}{'psi6':>7}")
    for r_ in rows:
        tz = f"{r_['t_zero']:.1f}" if np.isfinite(r_["t_zero"]) else "never"
        print(f"{r_['seed']:>5}{r_['dU']:>10.4f}{r_['dU_sem']:>10.4f}{r_['F']:>9.1f}"
              f"{r_['F_sem']:>10.1f}{r_['nd_drag']:>10.2f}{r_['nd_max']:>6.0f}"
              f"{tz:>11}{r_['slope']:>8.2f}{r_['psi6']:>7.3f}")

    def summ(key):
        v = np.array([r_[key] for r_ in rows], dtype=float)
        v = v[np.isfinite(v)]
        return v.mean(), v.std(ddof=1), v.std(ddof=1) / np.sqrt(len(v)), len(v)

    print("-" * 96)
    print(f"{'':>5}{'ens mean':>14}{'std across real':>12}{'ens SEM':>12}"
          f"{'mean run SEM':>13}{'std/runSEM':>14}")
    for key, lab, semkey in (("dU", "dU/particle", "dU_sem"),
                             ("F", "F_drag", "F_sem")):
        m_, s_, e_, k_ = summ(key)
        run_sem = np.mean([r_[semkey] for r_ in rows])
        print(f"{lab:>5}{m_:>14.4f}{s_:>12.4f}{e_:>12.4f}{run_sem:>13.4f}"
              f"{s_/run_sem:>14.2f}")
    for key, lab in (("nd_drag", "defects (driven)"), ("t_zero", "t(defects=0)"),
                     ("slope", "slope")):
        m_, s_, e_, k_ = summ(key)
        print(f"{lab:>5}{m_:>14.3f}{s_:>12.3f}{e_:>12.3f}"
              + (f"      ({k_}/{M} valid)" if k_ < M else ""))

    # ── figures ───────────────────────────────────────────────────────
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, NullFormatter
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    plain = FuncFormatter(lambda v, _: f"{v:g}")

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # (1) dU(t), individual + ensemble
    a = ax[0, 0]
    for u in U:
        a.plot(t_rel, u, "-", lw=.6, color="0.75")
    mu, se = U.mean(0), U.std(0, ddof=1) / np.sqrt(M)
    a.plot(t_rel, mu, "-", lw=2.2, color="tab:red", label=f"ensemble mean (M={M})")
    a.fill_between(t_rel, mu - se, mu + se, color="tab:red", alpha=.25, label="±SEM")
    a.axhline(0, color="green", ls="--", lw=1.2, label="equilibrium")
    a.set(xlabel=r"$t/\tau_{int}$ after relaxation starts",
          ylabel=r"$\Delta U_{total}$ [kT]",
          title="(1) relaxation -- individual realisations (grey) vs ensemble mean")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # (2) log-log -- does the sharp drop survive the average?
    a = ax[0, 1]
    pos = t_rel > 0
    edges = np.geomspace(t_rel[pos][0], t_rel[-1], 26)
    ctr, mm, ss = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (t_rel >= lo) & (t_rel < hi)
        if s.sum() >= 1:
            ctr.append(np.sqrt(lo * hi)); mm.append(U[:, s].mean()); ss.append(U[:, s].mean(1).std(ddof=1) / np.sqrt(M))
    ctr, mm, ss = map(np.array, (ctr, mm, ss))
    ok = mm > 0
    a.errorbar(ctr[ok], mm[ok], yerr=ss[ok], fmt="o-", ms=4, lw=1.5, capsize=2,
               color="tab:red", label="ensemble dU")
    for u in U:
        ub = np.array([u[(t_rel >= lo) & (t_rel < hi)].mean()
                       for lo, hi in zip(edges[:-1], edges[1:])
                       if ((t_rel >= lo) & (t_rel < hi)).sum() >= 1])
        m2 = ub > 0
        a.loglog(ctr[m2], ub[m2], "-", lw=.6, color="0.8", zorder=0)
    a.set_xscale("log"); a.set_yscale("log")
    for p_, c_ in ((-0.5, "tab:orange"), (-1.0, "tab:blue")):
        a.loglog(ctr[ok], mm[ok][0] * (ctr[ok] / ctr[ok][0]) ** p_, "--", lw=1,
                 color=c_, alpha=.7, label=fr"$t^{{{p_:g}}}$")
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_formatter(plain); axis.set_minor_formatter(NullFormatter())
    a.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$\Delta U_{total}$ [kT]",
          title="(2) log-log -- does the individual drop survive averaging?")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    # (3) defects
    a = ax[1, 0]
    for dd in D:
        a.plot(t_rel, dd, "-", lw=.6, color="0.8")
    mud, sed = D.mean(0), D.std(0, ddof=1) / np.sqrt(M)
    a.plot(t_rel, mud, "-", lw=2.2, color="tab:purple",
           label=f"ensemble mean (M={M})")
    a.fill_between(t_rel, mud - sed, mud + sed, color="tab:purple", alpha=.25)
    tz = [r_["t_zero"] for r_ in rows if np.isfinite(r_["t_zero"])]
    for x in tz:
        a.axvline(x, color="tab:green", lw=.8, alpha=.6)
    a.set(xlabel=r"$t/\tau_{int}$ after relaxation starts",
          ylabel="defect count (z!=6)",
          title=f"(3) defect relaxation -- green lines = each realisation reaching 0 "
                f"({len(tz)}/{M})")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # (4) spread across realisations
    a = ax[1, 1]
    dU = np.array([r_["dU"] for r_ in rows])
    F = np.array([r_["F"] for r_ in rows])
    a.errorbar(dU, F, xerr=[r_["dU_sem"] for r_ in rows],
               yerr=[r_["F_sem"] for r_ in rows], fmt="o", ms=7, capsize=3,
               color="tab:blue", label="individual runs (bars = each run's own SEM)")
    a.errorbar([dU.mean()], [F.mean()], xerr=[dU.std(ddof=1) / np.sqrt(M)],
               yerr=[F.std(ddof=1) / np.sqrt(M)], fmt="*", ms=20, capsize=4,
               color="crimson", label="ensemble mean ±SEM")
    a.axhline(F.mean(), color="crimson", ls=":", lw=1)
    a.axvline(dU.mean(), color="crimson", ls=":", lw=1)
    a.set(xlabel=r"$\Delta U$ [kT/particle]", ylabel=r"$F_{drag}$ [kT/d]",
          title="(4) spread across realisations -- are the per-run error bars honest?")
    a.legend(fontsize=8); a.grid(alpha=.3)

    fig.suptitle(f"trap-drag-2d-hex300 ensemble ({M} seeds, 2 lattice periods "
                 f"= 3.22 d)", fontsize=12)
    fig.tight_layout()
    p = out / "ensemble.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    (out / "ensemble.json").write_text(json.dumps(
        {"n_runs": M, "runs": rows,
         "t_rel": t_rel.tolist(), "dU_mean": U.mean(0).tolist(),
         "dU_sem": (U.std(0, ddof=1) / np.sqrt(M)).tolist(),
         "ndef_mean": D.mean(0).tolist()}, indent=2, ensure_ascii=False))
    # ══ force ensemble ═════════════════════════════════════════════════
    #   Single-sample noise on F in one run is k*.l_k = 246 kT/d, larger than the
    #   signal (~98). Stacking 8 runs reduces it by sqrt(8) = 2.8x, and **folding
    #   onto the lattice period** gives 2 periods x 8 runs = 16 periods, matching
    #   the statistics of a single 17-period run.
    Fx, Fy, PHz, Feq = [], [], [], []
    for d in dirs:
        spec, z, r, t, ph = load(d)
        k = spec.params["k_star"]; v = spec.params["v_star"]; a_nn = spec.params["a_nn_star"]
        dm = ph == PH_DRAG
        Fx.append(-k * z["dx_probe"][dm]); Fy.append(-k * z["dy_probe"][dm])
        td = t[dm] - t[dm][0]
        tiv = spec.reduced("times", "tau_int")
        PHz.append(((v * td * tiv) % a_nn) / a_nn)      # phase within a lattice period
        Feq.append(-k * z["dx_probe"][ph == PH_EQ])     # ★ control: trap stopped, so
                                                        #   <F> must be 0
    nf = min(len(a) for a in Fx)
    FX = np.array([a[:nf] for a in Fx]); FY = np.array([a[:nf] for a in Fy])
    tdrag = t_rel[:0]  # placeholder
    spec0, z0, r0, t0, ph0 = load(dirs[0])
    td0 = (t0[ph0 == PH_DRAG] - t0[ph0 == PH_DRAG][0])[:nf]
    v_star = spec0.params["v_star"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    w = max(9, nf // 40); kern = np.ones(w) / w
    for f in FX:
        a.plot(td0, np.convolve(f, kern, "same"), "-", lw=.6, color="0.8")
    mu, se = FX.mean(0), FX.std(0, ddof=1) / np.sqrt(M)
    a.plot(td0, np.convolve(mu, kern, "same"), "-", lw=2.2, color="tab:red",
           label=f"ensemble $F_x$ (M={M}, {w}-point moving average)")
    a.fill_between(td0, np.convolve(mu - se, kern, "same"),
                   np.convolve(mu + se, kern, "same"), color="tab:red", alpha=.25)
    a.axhline(v_star, color="k", ls="--", lw=1.3,
              label=f"bare Stokes $\\gamma v$ = {v_star:.1f}")
    a.axhline(FX.mean(), color="darkred", ls=":", lw=1.4,
              label=f"⟨$F_x$⟩ = {FX.mean():.1f} ± {FX.mean(1).std(ddof=1)/np.sqrt(M):.1f}")
    a.set(xlabel=r"$t/\tau_{int}$ after the drag starts", ylabel=r"$F_x$ [kT/d]",
          title="(1) force on the probe -- individual (grey) vs ensemble")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[0, 1]
    allp = np.concatenate([p_[:nf] for p_ in PHz]); allf = FX.ravel()
    nb = 16; edges = np.linspace(0, 1, nb + 1); ctr = .5 * (edges[:-1] + edges[1:])
    mu2 = np.array([allf[(allp >= l) & (allp < h)].mean() for l, h in zip(edges[:-1], edges[1:])])
    se2 = np.array([allf[(allp >= l) & (allp < h)].std() /
                    max(1, np.sqrt(((allp >= l) & (allp < h)).sum()))
                    for l, h in zip(edges[:-1], edges[1:])])
    a.errorbar(np.r_[ctr, ctr + 1], np.r_[mu2, mu2], yerr=np.r_[se2, se2], fmt="o-",
               ms=5, lw=1.5, capsize=3, color="tab:red",
               label=f"ensemble ({M} runs x 2 periods = {2*M} periods)")
    a.axhline(np.nanmean(mu2), color="darkred", ls=":", lw=1.3,
              label=f"mean {np.nanmean(mu2):.0f}")
    a.axhline(v_star, color="k", ls="--", lw=1.2, label="bare Stokes")
    a.axvline(1, color="0.7", lw=.8)
    amp2 = (np.nanmax(mu2) - np.nanmin(mu2)) / 2
    a.set(xlabel="phase within the lattice period (x mod $a_{NN}$)/$a_{NN}$",
          ylabel=r"$F_x$ [kT/d]",
          title=f"(2) folded force -- modulation ±{amp2:.0f} kT/d "
                f"({100*amp2/abs(np.nanmean(mu2)):.0f}%)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[1, 0]
    muy = FY.mean(0); sey = FY.std(0, ddof=1) / np.sqrt(M)
    a.plot(td0, np.convolve(muy, kern, "same"), "-", lw=2, color="tab:blue",
           label=r"ensemble $F_y$")
    a.fill_between(td0, np.convolve(muy - sey, kern, "same"),
                   np.convolve(muy + sey, kern, "same"), color="tab:blue", alpha=.25)
    a.axhline(0, color="k", ls="--", lw=1.2)
    fy_m = FY.mean(); fy_e = FY.mean(1).std(ddof=1) / np.sqrt(M)
    a.set(xlabel=r"$t/\tau_{int}$ after the drag starts", ylabel=r"$F_y$ [kT/d]",
          title=f"(3) transverse force -- must be 0 (sanity): "
                f"{fy_m:+.1f} ± {fy_e:.1f} = {abs(fy_m)/fy_e:.1f}σ")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[1, 1]
    FEQ = np.concatenate(Feq)
    for lab, arr, c in (("equilibrium (trap stopped)", FEQ, "tab:green"),
                        ("drag", allf, "tab:red")):
        a.hist(arr, bins=70, histtype="step", lw=1.6, color=c, density=True,
               label=f"{lab}: ⟨F⟩={arr.mean():+.1f}")
    a.axvline(0, color="k", ls="--", lw=1)
    a.axvline(v_star, color="0.4", ls=":", lw=1.4,
              label=f"bare Stokes {v_star:.0f}")
    eq_e = FEQ.std() / np.sqrt(len(FEQ))
    a.set(xlabel=r"$F_x$ [kT/d]", ylabel="probability density",
          title=f"(4) control -- equilibrium <F> = {FEQ.mean():+.1f} ± {eq_e:.1f} "
                f"({abs(FEQ.mean())/eq_e:.1f} sigma, must be 0)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    fig.suptitle(f"trap-drag ensemble -- force on the probe ({M} seeds)", fontsize=12)
    fig.tight_layout()
    pf = out / "ensemble_force.png"; fig.savefig(pf, dpi=140); plt.close(fig)
    print(f"\nforce ensemble: <F_x> = {FX.mean():.2f} ± {FX.mean(1).std(ddof=1)/np.sqrt(M):.2f} kT/d"
          f"  (bare Stokes {v_star:.2f}, +{100*(FX.mean()/v_star-1):.0f}%)")
    print(f"  F_y = {fy_m:+.2f} ± {fy_e:.2f} ({abs(fy_m)/fy_e:.1f}σ)   "
          f"equilibrium control F_x = {FEQ.mean():+.2f} ± {eq_e:.2f} "
          f"({abs(FEQ.mean())/eq_e:.1f} sigma)")
    print(f"  lattice modulation ±{amp2:.1f} kT/d = "
          f"{100*amp2/abs(np.nanmean(mu2)):.0f}% of the mean")
    print(f"\n{p}\n{pf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
