"""S2 -- derive the time-resolved sweep's predictions **from the prior run's measured
values** and seal them.

usage: python scripts/soft2d_time_series_predict.py [--out examples/.../prediction.yaml]

## Why a generator script, rather than writing it by hand

A prediction's `tolerance` **must come from the theoretical statistical error**
(`tests/conftest.py` rule 1). Reading the prior run's seed SE off the screen and
copying it makes it a hand calculation at that moment, and the prediction goes stale
silently when the value changes. This script reads `runs/<prior>/metrics.json` and
computes `3*sqrt(2)*SE`:

  . the difference between two 4-seed means has
    `SE_diff = sqrt(SE1^2 + SE2^2) ~ sqrt(2)*SE`
  . 3 sigma is used as the tolerance -> `tolerance = 3*sqrt(2)*SE`

## What this script is careful about

**The prior and new runs use different seeds** (prior `1-4`, new `5-8`). With the same
seeds, HOOMD `Brownian` reproduces bit-for-bit, which would turn the late-window
comparison into an **arithmetic identity** (`conftest` rule 3: a measurement that does
not fluctuate is an identity). The seeds have to differ for this to be a real test.

**The comparison window differs per `A`.** The prior runs had different total lengths
(`A <= 1` ran `30 tau_d`, `A = 10` ran `80 tau_d`). Their second halves are
`[20,30]` and `[60,80]` respectively, so those are the windows compared -- using one
window name for both would compare different times.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.analysis.structure import zahn_phase                    # noqa: E402
from simbot.build import (box_si_for_coverage,                      # noqa: E402
                          coverage_from_sigma_over_d)
from simbot.units import scales_soft2d, water_viscosity_si          # noqa: E402

PRIOR = REPO / "runs" / "2026-07-29_soft-r3-2d-A-sweep" / "metrics.json"
OUT_DEFAULT = REPO / "examples" / "soft-r3-time-resolved" / "prediction.yaml"

SIGMA_SI = 5.0e-6
D_OVER_SIGMA = 3.0
COVERAGE_MAX = 0.10
T_SI = 298.15
N_PARTICLES = 100
AMPLITUDES = (0.1, 1.0, 10.0)
NEW_SEEDS = (5, 6, 7, 8)

#  The times the prior runs' second halves covered [tau_d] -- total length differed
#  per A
PRIOR_WINDOW = {0.1: (20.0, 30.0), 1.0: (20.0, 30.0), 10.0: (60.0, 80.0)}
#  Whole-trajectory minimum separation in the prior runs (in d) -- it gets smaller as
#  the frame count grows
PRIOR_MIN_SEP = {0.1: 0.2172, 1.0: 0.4586, 10.0: 0.7164}

N_SIGMA = 3.0


def tol_from_se(se: float, n_sigma: float = N_SIGMA) -> float:
    """An `n_sigma` tolerance on the difference between two 4-seed means."""
    return float(n_sigma * np.sqrt(2.0) * se)


def main() -> int:
    out_path = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv \
        else OUT_DEFAULT
    prior = json.loads(PRIOR.read_text())

    geom = box_si_for_coverage(n_particles=N_PARTICLES, sigma_si=SIGMA_SI,
                               coverage_max=COVERAGE_MAX,
                               d_over_sigma_round=D_OVER_SIGMA)
    scales = scales_soft2d(d_si=geom["d_si"], sigma_si=SIGMA_SI, T_si=T_SI)
    eta_si, extrapolated = water_viscosity_si(T_SI)
    assert not extrapolated, "outside the water viscosity table, downgrade provenance"

    items: list[dict] = []

    # --- P1/P2/P7: does the late window agree with the prior run (independent
    #     seeds)? -------------------
    for A in AMPLITUDES:
        key = f"A{A:g}_random_sqbox"
        p = prior[key]
        lo, hi = PRIOR_WINDOW[A]
        for quantity, field, unit in (
                ("psi6_global", "psi6_global", ""),
                ("defect_fraction", "defect_fraction", ""),
                ("energy_per_particle", "energy_pp", "kT")):
            m, se = p[field]["mean"], p[field]["se"]
            items.append({
                "quantity": f"{quantity}__A{A:g}__t{lo:g}_{hi:g}",
                "value": float(m),
                "tolerance": f"±{tol_from_se(se):.6g}",
                "basis": (f"the prior run runs/2026-07-29_soft-r3-2d-A-sweep gave "
                          f"{key} = {m:.6g} ± {se:.4g} (SE) as the 4-seed mean over "
                          f"its second half. Tolerance = 3*sqrt(2)*SE, i.e. 3 sigma "
                          f"on the difference between two 4-seed means"),
                "discriminates": ("whether sampling continuously from t=0 gives the "
                                  "same answer as 'discard equilibration, keep the "
                                  "second half'. A difference means either the "
                                  "transient contaminated the late window or the two "
                                  "analysis paths have diverged"),
                "competing_value": None,
                **({"unit": unit} if unit else {}),
            })

    # --- P3: the number of coordination kinds -- the character of the defects -------
    for A in AMPLITUDES:
        kinds = int(prior[f"A{A:g}_random_sqbox"]["coord_kinds"])
        items.append({
            "quantity": f"coord_kinds__A{A:g}",
            "value": kinds,
            "tolerance": "±0",
            "basis": (f"the prior run A{A:g}_random_sqbox's coordination histogram has "
                      f"{kinds} kinds with a fraction > 0.5 %. Being an integer "
                      f"quantity, zero error is required"),
            "discriminates": ("the **character** of the defects. Equal fractions with "
                              "different numbers of kinds is different physics -- a "
                              "liquid (6 kinds) vs dislocations only (3 kinds). "
                              "Card §8.2"),
            "competing_value": None,
        })

    # --- P4: the initial condition's min_sep=0.8 d shell fills in ------------------
    #  `random_2d_snapshot` enforces `min_sep = 0.8 d` by rejection sampling, so at
    #  t=0, g(r < 0.8 d) = 0 holds **exactly**. An equilibrium liquid does not.
    #  Time for free diffusion to fill the gap 0.8 -> 0.5 d:
    #  MSD* = 4t* => t* ~ 0.3^2/4 = 0.0225.
    fill_t = float(0.3**2 / 4.0)
    items.append({
        "quantity": "rdf_at_0.5d__A0.1__first_window",
        "value": ">0",
        "tolerance": ">0",
        "basis": (f"the initial placement enforces min_sep = 0.8 d by rejection "
                  f"sampling, so g(r<0.8d) = 0 holds exactly at t=0. This is an "
                  f"**initial-condition artefact** and free diffusion fills it: "
                  f"MSD* = 4t* => a gap of 0.3 d takes t* ~ {fill_t:.4g} tau_d. "
                  f"The first window (0-20 tau_d) is {20/fill_t:.0f}x longer"),
        "discriminates": ("whether the initial condition's excluded-volume shell "
                          "actually relaxes. If it survives, the sampling is trapped "
                          "in the initial condition"),
        "competing_value": None,
    })

    # --- P5: reference-disc overlap -- a prediction that can only be stated once
    #     sigma is attached -------------
    sigma_over_d = geom["sigma_over_d"]
    for A in AMPLITUDES:
        prior_ms = PRIOR_MIN_SEP[A]
        overlaps = prior_ms < sigma_over_d
        items.append({
            "quantity": f"min_separation_over_sigma__A{A:g}",
            "value": f"<{1.0:g}" if overlaps else f">{1.0:g}",
            "tolerance": f"<{1.0:g}" if overlaps else f">{1.0:g}",
            "basis": (f"the prior run's min separation {prior_ms:.4g} d = "
                      f"{prior_ms / sigma_over_d:.4g} σ (σ/d = {sigma_over_d:.6g}). "
                      f"the new run has 2x the frames and "
                      f"{80.0 / (PRIOR_WINDOW[A][1]):.3g}x the length, so the extreme "
                      f"value gets smaller -> the overlap verdict holds or "
                      f"strengthens"),
            "discriminates": ("whether the picture of '100 discs of 5 µm' is "
                              "physically consistent at this A. A/r^3 has no hard "
                              "core, so the model permits overlap -- if the picture "
                              "does not hold, that is the result"),
            "competing_value": None,
            "unit": "sigma",
        })

    # --- P6: the ordering of arrival times -- the transition point is slow ----------
    items.append({
        "quantity": "t90_ordering",
        "value": "t90(A=10) > t90(A=1) and t90(A=10) > 1 tau_d",
        "tolerance": "qualitative (ordering)",
        "basis": ("A=10 gives Gamma = pi^{3/2}*10 = 55.68, which is -7.0 % from Zahn's "
                  "transition point 59.88 [cited, not reproduced] -> critical slowing "
                  "is expected. A <= 1 is deep liquid, and a random initial placement "
                  "is already close to a liquid configuration, so psi6 should be "
                  "nearly flat"),
        "discriminates": "whether relaxation slows near the transition (critical slowing)",
        "competing_value": "arrival times at all three A are O(1) tau_d, indistinguishable",
    })

    doc = {
        "items": items,
        "regimes": {
            "amplitudes": list(AMPLITUDES),
            "gamma_zahn": {f"A{A:g}": zahn_phase(A)["gamma"] for A in AMPLITUDES},
            "phase_zahn": {f"A{A:g}": zahn_phase(A)["phase_zahn"]
                           for A in AMPLITUDES},
            "sigma_si": SIGMA_SI,
            "d_si": geom["d_si"],
            "L_si": geom["L_si"],
            "L_star": geom["L_star"],
            "coverage": geom["coverage"],
            "coverage_max": COVERAGE_MAX,
            "sigma_over_d": geom["sigma_over_d"],
            "tau_d_si": scales.time_si,
            "D0_si": scales.diffusivity_si,
            "eta_si": eta_si,
            "T_si": T_SI,
            "seeds": list(NEW_SEEDS),
            "note": ("the reduced physics is determined by A alone (n* = 1, no hard "
                     "core). Coverage changes only tau_d's value in seconds and the "
                     "validity of the point-particle idealisation -- its sensitivity "
                     "to psi6, g(r) and the defects is exactly 0"),
        },
        "alternatives": [
            "the seeds were made disjoint from the prior run (1-4 -> 5-8). With the "
            "same seeds, HOOMD Brownian reproduces bit-for-bit and the late-window "
            "comparison becomes an arithmetic identity.",
            "INCONCLUSIVE expected -- 'equilibrium was reached' cannot be decided by "
            "this run. In the prior run A=10's psi6 drift was -0.013 ± 0.008 "
            "(1.6 sigma), indistinguishable from 0 -- and indistinguishable means 'not "
            "visible at this length', not 'equilibrium'.",
            "INCONCLUSIVE expected -- the psi6 difference between A=0.1 and A=1 "
            "(0.049 vs 0.063) is below the finite-size floor (1/sqrt(N) = 0.1). Their "
            "orientational order can only be described as 'absent in both'.",
            "N=100 does not satisfy A=10's r_cut requirement (N >= 252, card §9). "
            "beta*U(r_cut) = 0.09 kT remains in the A=10 result -- an N convergence "
            "check is outside this run's scope.",
            "reference discs overlapping at A=0.1 is a property of the model, not a "
            "code error. Removing the overlap requires dropping coverage below 3.71 %, "
            "which leaves the reduced results unchanged and only alters tau_d.",
        ],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "# S2 PREDICTION -- write the answer down before simulating  "
        "(auto-generated: "
        "scripts/soft2d_time_series_predict.py)\n"
        "#\n"
        "# ⚠ This file is sealed before S5 runs (SEALED.sha256). Do not edit after "
        "running.\n"
        "# The tolerances were derived from the prior run's seed SE (3*sqrt(2)*SE) -- "
        "not tailored to the observed values.\n"
        "#\n"
        + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8")
    print(f"-> {out_path.relative_to(REPO)}  ({len(items)} items)")
    for it in items:
        print(f"   {it['quantity']:<46} {str(it['value']):>12}  {it['tolerance']}")
    print(f"\n  d = {geom['d_si']*1e6:.3f} µm · L = {geom['L_si']*1e6:.2f} µm · "
          f"coverage = {geom['coverage']:.4%} · τ_d = {scales.time_si:.1f} s "
          f"= {scales.time_si/60:.2f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
