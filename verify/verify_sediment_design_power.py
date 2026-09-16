"""Design power for the sedimentation campaign — from the Carnahan–Starling
profiles alone, using NO run data.

    ./bin/py verify/verify_sediment_design_power.py

★ Why this exists, and why it exists late. `prediction.yaml` revisions 1 and 2
stated their gates as fixed percentages — `l_g_fitted_tail` within 5 % of 13.375,
`Z_dilute_tail` within 5 % of 1, `D_xy` within 10 % — and never computed the
noise those percentages were supposed to sit outside of. A tolerance is only a
gate if it is wider than its own statistic's noise. Measured here:

    l_g_fitted_tail at the primary arm:  BIAS +10.6 %, 1 sigma 5.9 %
                                         against a 5 % tolerance.

So that gate's band was narrower than its own 1 sigma, and the sealed decision
rule says a failed correctness gate means "report IMPLEMENTATION FAILURE and no
physics". The campaign was about 82 % likely, per seed, to call a working
implementation broken.

⚠ **No production run had completed when this was written.** P1 was stopped at
20 % for an unrelated defect (the sealed plan required a first-half/second-half
profile comparison the output could not produce) and no run directory holds a
`metrics.json`. The numbers below are properties of the DESIGN. Revising a sealed
document on them is what the pre-registration window is for; revising it after
seeing an outcome would not be.

## How it works

The two competing hypotheses are made concrete as two profiles:

    CS(phi_nominal)  -- `cs_hydrostatic_profile(scale=1)`,      phi(0) = 0.10306
    CS(phi_eff)      -- `..., scale=(d_BH/d)**3 = 0.74073`,     phi(0) = 0.11237

Poisson counts are drawn on each, and **the case's own estimators are imported**
and run on them. A power calculation done against a second copy of the estimator
measures the copy.

Two quantities come out for every observable:

| | |
|---|---|
| **bias and sigma** | generated from the hypothesis the observable is compared against. Says whether a tolerance is a gate or a coin flip |
| **separation** | `|mean_eff − mean_nom| / sigma`. Says whether the observable can tell the two hypotheses apart at all. Under 2 sigma it cannot, whatever its tolerance |

⚠ The first version generated from `scale = 1` only and compared against
CS(phi_eff), so the primary statistic read +6.49 % and looked like a catastrophic
estimator bias. It is not — it is the MAPPING GAP, which is the thing the
campaign measures. Generating from the matching hypothesis is what makes a null a
null, and that mistake is why `expected_counts` takes `scale` explicitly.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
import cases.sediment_3d as SED  # noqa: E402

L_G, LZ, PHI = 13.375, 125.0, 0.017

#: (d_BH/d)**3 at eps_wca = 1 kT, sigma_LJ = 2**(-1/6). The mapping under test.
SCALE_EFF = 0.9048 ** 3

N_TRIAL = 800

#: (label, L_xy/d, T_obs in tau_sed, samples per tau_sed)
ARMS = (("cs_eff", 12.0, 4.0, 40), ("cs_nominal", 12.0, 4.0, 40),
        ("convergence", 12.0, 8.0, 40), ("finite_size", 24.0, 3.0, 40))


def expected_counts(lxy: float, frames: int, *, scale: float):
    """Counts per bin a CS equilibrium profile would deposit. `scale` picks WHICH
    hypothesis generates the data: 1 is "no mapping", `SCALE_EFF` is "mapping
    required"."""
    z, phi = SED.cs_hydrostatic_profile(L_G, LZ, PHI, scale=scale)
    mid = np.arange(SED.BIN_D / 2, LZ, SED.BIN_D)
    vbin = lxy * lxy * SED.BIN_D
    return mid, np.interp(mid, z, phi) * vbin / (math.pi / 6.0) * frames


def observables(rng, mid, lam, frames, lxy, n=N_TRIAL) -> dict:
    """Every candidate observable, over `n` synthetic realisations."""
    vbin = lxy * lxy * SED.BIN_D
    hi = min(6.0 * L_G, LZ - 2.0)
    keys = ("l_g_fitted_tail", "l_g_tail_se", "l_g_fitted_dense", "Z_dilute_tail",
            "phi_wall_measured", "Z_dev_wmean_vs_CS_eff", "Z_dev_wmean_vs_CS_nominal",
            "Z_dev_max_vs_CS_eff")
    out = {k: [] for k in keys}
    for _ in range(n):
        counts = rng.poisson(lam)
        phi_z = counts / frames * (math.pi / 6.0) / vbin
        z_eos = SED.eos_from_profile(phi_z, SED.BIN_D, L_G)
        ft = SED.fit_decay_length(mid, counts, SED.TAIL_LO, min(SED.TAIL_HI, LZ - 2))
        fd = SED.fit_decay_length(mid, counts, SED.FIT_LO, hi)
        if ft:
            out["l_g_fitted_tail"].append(ft["l_g"])
            out["l_g_tail_se"].append(ft["se"])
        if fd:
            out["l_g_fitted_dense"].append(fd["l_g"])
        dil = (mid >= SED.TAIL_LO) & (counts >= 50) & np.isfinite(z_eos)
        if dil.any():
            out["Z_dilute_tail"].append(float(np.nanmean(z_eos[dil])))
        out["phi_wall_measured"].append(float(phi_z[mid < 1.5].max()))
        ok = (mid >= SED.FIT_LO) & (counts >= 200) & np.isfinite(z_eos)
        if ok.any():
            w = counts[ok].astype(float)
            ce = SED.z_carnahan_starling(phi_z[ok] * SCALE_EFF)
            cn = SED.z_carnahan_starling(phi_z[ok])
            de = 100.0 * (z_eos[ok] - ce) / ce
            out["Z_dev_wmean_vs_CS_eff"].append(float((w * de).sum() / w.sum()))
            out["Z_dev_wmean_vs_CS_nominal"].append(
                float((w * 100.0 * (z_eos[ok] - cn) / cn).sum() / w.sum()))
            out["Z_dev_max_vs_CS_eff"].append(float(np.max(np.abs(de))))
    return {k: np.array(v) for k, v in out.items() if len(v)}


def apparent_l_g(lo: float, hi: float, *, scale: float) -> float:
    """The decay length a fit over `[lo, hi]` returns on the CONTINUOUS generating
    profile — i.e. with no counting noise and no count floor. Useful only as a
    contrast: it is what revision 2's P1b prediction was computed from, and the
    Monte Carlo shows the count floor moves it by 13 %."""
    z, phi = SED.cs_hydrostatic_profile(L_G, LZ, PHI, scale=scale)
    m = (z >= lo) & (z <= hi) & (phi > 1e-12)
    return float(-1.0 / np.polyfit(z[m], np.log(phi[m]), 1)[0])


#: (observable, the value it is compared against, the hypothesis that generates
#:  its null, whether it is expected to discriminate)
TARGETS = (
    ("l_g_fitted_tail", L_G, "eff", False),
    ("l_g_fitted_dense", None, "eff", True),
    ("Z_dilute_tail", 1.0, "eff", False),
    ("phi_wall_measured", None, "eff", True),
    ("Z_dev_wmean_vs_CS_eff", 0.0, "eff", True),
    ("Z_dev_max_vs_CS_eff", None, "eff", False),
)


def main() -> int:
    rng = np.random.default_rng(20260916)
    print(f"{N_TRIAL} realisations per arm per hypothesis, "
          f"scale_eff = {SCALE_EFF:.6f}")
    for label, s in (("CS(phi_nominal)", 1.0), ("CS(phi_eff)", SCALE_EFF)):
        z, phi = SED.cs_hydrostatic_profile(L_G, LZ, PHI, scale=s)
        print(f"  {label:16s} phi(0) = {phi[0]:.5f}   continuous-profile "
              f"l_g over [2,80] = {apparent_l_g(2.0, 80.0, scale=s):.3f} d")

    problems = []
    for arm, lxy, tau, spt in ARMS:
        frames = int(tau * spt)
        gen = {}
        for tag, s in (("eff", SCALE_EFF), ("nom", 1.0)):
            mid, lam = expected_counts(lxy, frames, scale=s)
            gen[tag] = observables(rng, mid, lam, frames, lxy)

        print(f"\n── {arm}: L_xy = {lxy:g} d, T_obs = {tau:g} tau_sed, "
              f"{frames} frames " + "─" * 18)
        print(f"{'observable':28s} {'from CS(eff)':>19s} {'from CS(nom)':>19s} "
              f"{'sep':>7s} {'bias':>8s} {'3 sigma':>8s}")
        for name, target, null, should_sep in TARGETS:
            if name not in gen["eff"] or name not in gen["nom"]:
                continue
            e, n_ = gen["eff"][name], gen["nom"][name]
            sd = max(e.std(), n_.std())
            sep = abs(e.mean() - n_.mean()) / sd if sd > 0 else float("inf")
            ref = e if null == "eff" else n_
            if target is None:
                bias_s, t3 = "  (n/a)", f"{3*ref.std():7.4f}"
            elif target == 0.0:
                bias_s = f"{ref.mean():+7.3f}"
                t3 = f"{3*ref.std():7.3f}"
            else:
                bias_s = f"{100*(ref.mean()/target - 1):+6.2f}%"
                t3 = f"{300*ref.std()/target:6.2f}%"
            print(f"{name:28s} {e.mean():11.4f} ±{e.std():6.4f} "
                  f"{n_.mean():11.4f} ±{n_.std():6.4f} {sep:6.1f}σ {bias_s:>8s} "
                  f"{t3:>8s}")
            if should_sep and sep < 3.0:
                problems.append(f"{arm}/{name}: only {sep:.1f} sigma of "
                                f"separation -- it cannot decide the question")
            if target not in (None, 0.0):
                b = abs(ref.mean() / target - 1) * 100
                s3 = 300 * ref.std() / target
                if b > s3:
                    problems.append(
                        f"{arm}/{name}: bias {b:.1f} % exceeds 3 sigma = "
                        f"{s3:.1f} %, so no symmetric tolerance around {target} "
                        f"is both passable and informative")

        if "l_g_tail_se" in gen["eff"]:
            se, sd = gen["eff"]["l_g_tail_se"].mean(), gen["eff"]["l_g_fitted_tail"].std()
            print(f"{'  (tail fit reports its own se)':28s} "
                  f"reported {se:.4f} d   actual sd {sd:.4f} d   "
                  f"ratio {se/sd:.3f}")
            if not 0.85 < se / sd < 1.15:
                problems.append(f"{arm}: the tail fit's reported se is "
                                f"{se/sd:.2f}x its actual spread")

    print("\n" + "=" * 78)
    for p in problems:
        print(f"  x {p}")
    print(f"\n{len(problems)} design problem(s)" if problems else
          "\nevery gate's tolerance clears its own noise, and every observable "
          "the verdict rests on separates the two hypotheses")
    return 1 if problems else 0


if __name__ == "__main__":
    raise SystemExit(main())
