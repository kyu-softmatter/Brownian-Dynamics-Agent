"""S2 -- derive the `N` convergence run's predictions **from the measured `N=100`
values** and seal them.

usage: python scripts/soft2d_nconv_predict.py

## This run's first-class question is NOT the truncation error

`N=100 -> 256` lowers `beta*U(r_cut)` from `0.0904` to `0.0216 kT` (the open item in
card §9). But **the finite-size exponent of `psi6` is worth more.**

`|⟨ψ₆⟩| ~ N^{-p}`, `η₆ = 4p`:
  . liquid (`g6` decays exponentially)  -> `p = 1/2`
  . hexatic (`g6 ~ r^{-eta6}`)          -> `p = eta6/4 <= 1/16`
                                           (KTHNY boundary `eta6 = 1/4`)
  . crystal                             -> `p -> 0`

`A=10` was read as "not a crystal, but a dislocation-rich hexatic-like state" from
`psi6 = 0.248` (card §8.3), and **the exponent tests that reading.** If `p ~ 0.5`
then that `0.248` is merely a finite-size floor; if `p` is small it is genuine
quasi-long-range orientational order.

=> This is the route to the `eta6` that Zahn reproduction condition §6-3 requires,
**without fitting `g6(r)`.**
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,  # noqa: E402
                                       LIQUID_EXPONENT_P, zahn_phase)
from simbot.build import box_si_for_coverage                        # noqa: E402
from simbot.units import scales_soft2d                              # noqa: E402

PARENT = REPO / "runs" / "2026-07-29_soft-r3-time-resolved" / "metrics.json"
OUT = REPO / "examples" / "soft-r3-nconv" / "prediction.yaml"

N_REF, N_NEW = 100, 256
AMPLITUDES = (0.1, 1.0, 10.0)
SEEDS = (21, 22, 23, 24)          # ★ kept disjoint from the parent run (5-8) and the
                                  # relaxation sweep (5-1540)
SIGMA_SI, D_OVER_SIGMA, T_SI = 5.0e-6, 3.0, 298.15

#  Gates measured at N=256 (scripts/soft2d_nconv.py recomputes and compares the same
#  values)
BETA_U_RCUT = {0.1: 0.0002, 1.0: 0.0022, 10.0: 0.0216}
N_SIGMA = 3.0


def main() -> int:
    parent = json.loads(PARENT.read_text())
    ratio = N_NEW / N_REF
    ln_ratio = float(np.log(ratio))
    items: list[dict] = []

    for A in AMPLITUDES:
        ref = parent[f"A{A:g}"]
        p6, p6se = ref["psi6_global"]["mean"], ref["psi6_global"]["se"]

        # --- (1) psi6 finite-size prediction: both extremes are written down ---
        liquid_pred = p6 * ratio ** (-LIQUID_EXPONENT_P)
        hexatic_pred = p6 * ratio ** (-KTHNY_ETA6_HEXATIC_LIQUID / 4.0)
        #  Tolerance: the parent SE scaled, then 3*sqrt(2) (assuming the new SE is
        #  comparable)
        tol = 3.0 * np.sqrt(2.0) * p6se * ratio ** (-LIQUID_EXPONENT_P)
        items.append({
            "quantity": f"psi6_global__A{A:g}__N256",
            "value": float(liquid_pred),
            "tolerance": f"±{tol:.6g}",
            "basis": (
                f"liquid hypothesis `|<psi6>| ~ N^-1/2`: {p6:.6g} x "
                f"({ratio:g})^-0.5 = "
                f"{liquid_pred:.6g}. **The competing hypothesis (hexatic, "
                f"eta6=1/4 -> p=1/16) gives "
                f"{hexatic_pred:.6g}**, and the gap between them is "
                f"{abs(hexatic_pred - liquid_pred):.6g} = "
                f"{abs(hexatic_pred - liquid_pred) / tol:.1f}x the tolerance -- so "
                f"they are distinguishable"),
            "discriminates": (
                "whether A=%g's psi6 is a finite-size floor or genuine "
                "quasi-long-range order" % A),
            "competing_value": float(hexatic_pred),
        })

        # --- (2) the finite-size exponent itself ---
        items.append({
            "quantity": f"psi6_exponent_p__A{A:g}",
            "value": LIQUID_EXPONENT_P,
            "tolerance": "±0.15",
            "basis": (
                f"`p = -dln|<psi6>|/dlnN`. A liquid gives 0.5, the KTHNY hexatic "
                f"boundary gives {KTHNY_ETA6_HEXATIC_LIQUID / 4.0:.4f}. The tolerance "
                f"±0.15 is the scale of error propagation for a two-point estimate "
                f"(parent SE {p6se:.4g}, relative {p6se / p6:.1%}, "
                f"/ln{ratio:g}={ln_ratio:.4f})"),
            "discriminates": f"the phase reading for A={A:g} (eta6 = 4p)",
            "competing_value": KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
        })

        # --- (3) local quantities must be independent of N ---
        for q, field in (("defect_fraction", "defect_fraction"),
                         ("psi6_local", "psi6_local")):
            m, se = ref[field]["mean"], ref[field]["se"]
            items.append({
                "quantity": f"{q}__A{A:g}__N256",
                "value": float(m),
                "tolerance": f"±{N_SIGMA * np.sqrt(2.0) * se:.6g}",
                "basis": (
                    f"**a local quantity, so it must be independent of `N`.** "
                    f"Observed at `N=100`: {m:.6g} ± {se:.4g}. Tolerance "
                    f"3*sqrt(2)*SE. A disagreement means `N=100` carried a "
                    f"finite-size effect"),
                "discriminates": f"whether {q} had already converged at N=100",
                "competing_value": None,
            })

    # --- (4) truncation error (card §9's original open item) ---
    for A in AMPLITUDES:
        items.append({
            "quantity": f"beta_u_at_rcut__A{A:g}__N256",
            "value": BETA_U_RCUT[A],
            "tolerance": "±0.0005",
            "basis": (
                f"`r_cut = 0.98·L/2 − buffer`, `L* = √256 = 16` → `r_cut = 7.740`. "
                f"`beta*U(r_cut) = {A:g}/7.740^3`. Deterministic, so the only error is "
                f"rounding. At `N=100` it was `r_cut = 4.800`"),
            "discriminates": "whether card §9's truncation error actually decreased",
            "competing_value": None,
        })

    # --- (5) 5-7 imbalance and coordination kinds (the character of the defects) ---
    items.append({
        "quantity": "coord_kinds_aggregate__A10__N256",
        "value": int(parent["A10"]["coord_kinds_aggregate"]),
        "tolerance": "±0",
        "basis": (
            f"{parent['A10']['coord_kinds_aggregate']} kinds in the `N=100` aggregate. "
            f"An **aggregate estimator** is used -- a per-frame threshold depends on "
            f"`N` (one particle is 1 % at `N=100` but 0.39 % at `N=256`). "
            f"Basis: findings/fraction-threshold-flips-meaning-*.md"),
        "discriminates": "whether the character of the defects is independent of N",
        "competing_value": None,
    })

    geo_ref = box_si_for_coverage(n_particles=N_REF, sigma_si=SIGMA_SI,
                                  coverage_max=0.10,
                                  d_over_sigma_round=D_OVER_SIGMA)
    geo_new = box_si_for_coverage(n_particles=N_NEW, sigma_si=SIGMA_SI,
                                  coverage_max=0.10,
                                  d_over_sigma_round=D_OVER_SIGMA)
    sc = scales_soft2d(d_si=geo_new["d_si"], sigma_si=SIGMA_SI, T_si=T_SI)

    doc = {
        "question": (
            "from N=100 -> 256: (1) does psi6 fall as 1/sqrt(N) (i.e. is it a "
            "finite-size floor)? (2) do the local quantities (defect fraction, local "
            "psi6) stay unchanged? (3) does the truncation error decrease?"),
        "parent_run": "runs/2026-07-29_soft-r3-time-resolved",
        "card_open_item": "card §9 (the N >= 252 requirement) . §10 (the g6 exponent "
                          "eta6 unimplemented)",
        "items": items,
        "regimes": {
            "n_ref": N_REF, "n_new": N_NEW, "n_ratio": ratio,
            "seeds": list(SEEDS),
            "amplitudes": list(AMPLITUDES),
            "gamma_zahn": {f"A{A:g}": zahn_phase(A)["gamma"] for A in AMPLITUDES},
            "L_star_ref": geo_ref["L_star"], "L_star_new": geo_new["L_star"],
            "coverage": geo_new["coverage"],
            "note_coverage": (
                "★ Coverage is independent of N (with n* = 1 the area per particle is "
                "exactly d^2). Only the box grows, as sqrt(N) -- L = 150 -> 240 µm"),
            "L_si_ref": geo_ref["L_si"], "L_si_new": geo_new["L_si"],
            "tau_d_si": sc.time_si,
            "exponent_liquid": LIQUID_EXPONENT_P,
            "exponent_hexatic_boundary": KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
            "eta6_hexatic_boundary": KTHNY_ETA6_HEXATIC_LIQUID,
        },
        "alternatives": [
            "★ The most important competing hypothesis: that A=10's psi6 = 0.248 is a "
            "**finite-size floor**. If p ~ 0.5 comes out, card §8.3's 'hexatic-like' "
            "reading must be weakened -- it would mean 0.248 was merely the floor at "
            "N=100.",
            "Two points cannot verify the power-law **form**. This extracts an "
            "exponent only, and psi6_finite_size_exponent's n_points carries that "
            "limitation. Discussing the form needs three or more (e.g. N=64, 144, "
            "400).",
            "A=0.1 and A=1 should give p ~ 0.5 (deep liquid). If those also come out "
            "small, the exponent estimate itself is wrong -- so **A <= 1 serves as the "
            "control**.",
            "Energy per particle cannot be compared directly across N -- "
            "power_law_table shifts U(r_cut) to 0, so changing r_cut changes the "
            "shift. It is recorded but not raised as a prediction item.",
            "The coverage control (3.71 %) is **not run.** The reduced config is "
            "bit-identical, making it an arithmetic identity -- pinned by "
            "tests/test_s7_structure.py::"
            "test_coverage_does_not_touch_the_reduced_config_at_all.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# S2 PREDICTION -- N convergence + the psi6 finite-size exponent  "
        "(auto-generated: "
        "scripts/soft2d_nconv_predict.py)\n"
        "#\n# ⚠ This file is sealed before S5 runs (SEALED.sha256). Do not edit after "
        "running.\n"
        "# The tolerances were derived by scaling the parent run's (N=100) seed SE.\n"
        "#\n"
        + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8")
    print(f"-> {OUT.relative_to(REPO)}  ({len(items)} items)\n")
    print(f"{'item':<40} {'predicted':>12} {'competing':>12} {'tolerance':>12}")
    print("-" * 80)
    for it in items:
        cv = it.get("competing_value")
        cvs = f"{cv:.5g}" if isinstance(cv, (int, float)) else "—"
        v = it["value"]
        vs = f"{v:.5g}" if isinstance(v, (int, float)) else str(v)
        print(f"{it['quantity']:<40} {vs:>12} {cvs:>12} {it['tolerance']:>12}")
    print(f"\n  N: {N_REF} → {N_NEW} · L*: {geo_ref['L_star']:.0f} → "
          f"{geo_new['L_star']:.0f} · L: {geo_ref['L_si']*1e6:.0f} → "
          f"{geo_new['L_si']*1e6:.0f} µm . coverage {geo_new['coverage']:.4%} "
          f"(unchanged)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
