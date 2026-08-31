"""S2 -- derive the finite-size ladder (S17) predictions from the two earlier runs and
seal them.

usage: python scripts/soft2d_fss_predict.py

## What this run tests is **the method**

`runs/2026-07-29_soft-r3-nconv` obtained `p` from two points (`N=100, 256`) and
claimed "`A=10` is not hexatic" at `3.5 sigma`. Two points **assume** a power law.
This run tests that assumption with four -- if the form is wrong, `eta6 = 1.46` is
unusable too.

Three things change at the same time, deliberately, so that this is an **independent
measurement**:
  . `r_cut` fixed at 3.80 (same for every N) vs the earlier runs' natural 4.80/7.74
  . `prod_tau` 30 tau_d vs 80 . window [20,30] vs [60,80]
  . seeds 32-35 vs 21-24
=> getting the same `p` means **the method is robust**; getting a different one means
  narrowing down which axis caused it.
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
                                       LIQUID_EXPONENT_P)

NCONV = REPO / "runs" / "2026-07-29_soft-r3-nconv" / "metrics.json"
OUT = REPO / "examples" / "soft-r3-fss" / "prediction.yaml"
N_LADDER = (64, 144, 256, 400)
AMPLITUDES = (0.1, 1.0, 10.0)
CHI2_MAX = 3.0


def main() -> int:
    nc = json.loads(NCONV.read_text())
    items: list[dict] = []

    # --- (1) form verification: the reason this run exists --------------------
    for A in AMPLITUDES:
        items.append({
            "quantity": f"chi2_reduced__A{A:g}",
            "value": 1.0,
            "tolerance": f"<{CHI2_MAX:g}",
            "basis": (
                f"`chi^2/dof` of a straight-line fit in log-log (dof = 4-2 = 2). If "
                f"the power law holds and the error bars are honest, the expectation "
                f"is 1. `< {CHI2_MAX:g}` is pre-registered as the pass criterion -- "
                f"above it, the data must not be summarised by a single exponent"),
            "discriminates": "whether psi6(N) really is a power law (the premise of "
                             "the two-point estimate)",
            "competing_value": None,
        })

    # --- (2) exponent: the liquid control plus a re-measurement of A=10 --------
    for A in AMPLITUDES:
        f = nc[f"A{A:g}"]["finite_size"]
        if A <= 1.0:
            items.append({
                "quantity": f"psi6_exponent_p__A{A:g}",
                "value": LIQUID_EXPONENT_P,
                "tolerance": "±0.12",
                "basis": (
                    f"deep in the liquid, so `p = 1/2`. The two-point run gave "
                    f"{f['p']:.3f} ± {f['p_se']:.3f}, and four points should shrink "
                    f"the SE. The tolerance ±0.12 is about 1.5x the two-point SE "
                    f"({f['p_se']:.3f}) -- **kept tight because this is the "
                    f"control**"),
                "discriminates": "whether the exponent estimate is healthy "
                                 "(the control)",
                "competing_value": KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
            })
        else:
            items.append({
                "quantity": f"psi6_exponent_p__A{A:g}",
                "value": float(f["p"]),
                "tolerance": f"±{3.0 * np.sqrt(2.0) * f['p_se']:.4g}",
                "basis": (
                    f"re-measures the two-point run's {f['p']:.3f} ± "
                    f"{f['p_se']:.3f}. `r_cut`, run length and seeds all differ, so "
                    f"this is an **independent measurement**. Tolerance 3*sqrt(2)*SE"),
                "discriminates": (
                    "whether the two-point estimate is robust to r_cut, run length "
                    "and seeds"),
                "competing_value": LIQUID_EXPONENT_P,
            })

    # --- (3) does A=10 still reject hexatic? ----------------------------------
    items.append({
        "quantity": "eta6_minus_3sigma__A10",
        "value": 0.5,
        "tolerance": f">{KTHNY_ETA6_HEXATIC_LIQUID:g}",
        "basis": (
            f"if `eta6 - 3*SE` exceeds the hexatic upper bound "
            f"{KTHNY_ETA6_HEXATIC_LIQUID:g}, hexatic is rejected at `3 sigma`. In the "
            f"two-point run it was "
            f"{nc['A10']['finite_size']['eta6']:.2f} − 3×"
            f"{nc['A10']['finite_size']['eta6_se']:.2f} = "
            f"{nc['A10']['finite_size']['eta6'] - 3*nc['A10']['finite_size']['eta6_se']:.2f} "
            f". The predicted 0.5 is a placeholder meaning 'comfortably exceeds'"),
        "discriminates": "whether A=10's hexatic rejection survives at four points",
        "competing_value": None,
    })

    # --- (4) is psi6 unchanged across an 8.4-fold r_cut excursion? (S16 extended) ---
    for A in AMPLITUDES:
        ref = nc[f"A{A:g}"]
        m, se = ref["psi6_global"]["mean"], ref["psi6_global"]["se"]
        items.append({
            "quantity": f"psi6_global__A{A:g}__N256",
            "value": float(m),
            "tolerance": f"±{3.0 * np.sqrt(2.0) * se:.4g}",
            "basis": (
                f"re-measures the same `N=256` at `r_cut = 3.80` (the earlier run used "
                f"7.740 -- an **8.4-fold difference in truncation error**, "
                f"beta*U: 0.0216 -> 0.182). The earlier run gave "
                f"{m:.5g} ± {se:.4g}. Agreement means S16 (that truncation error does "
                f"not bias the observables) holds across a larger excursion too. The "
                f"run length and seeds differ as well"),
            "discriminates": "whether an 8.4-fold truncation error biases psi6",
            "competing_value": None,
        })

    doc = {
        "question": (
            "whether psi6(N) really is a power law (4 points) . whether that exponent "
            "is robust to r_cut, run length and seeds . whether A=10's hexatic "
            "rejection survives"),
        "parent_runs": ["runs/2026-07-29_soft-r3-nconv",
                        "runs/2026-07-29_soft-r3-time-resolved"],
        "card_open_item": "card §8.5 S17 (power-law form verification)",
        "design": {
            "n_ladder": list(N_LADDER),
            "lever_arm_ln": float(np.log(N_LADDER[-1] / N_LADDER[0])),
            "amplitudes": list(AMPLITUDES),
            "seeds": [32, 33, 34, 35],
            "seed_screening": (
                "min_sep = 0.8 d rejection sampling fails for some seeds (measured "
                "success rates over seeds 31-90: N=64 98.3 %, N=144 95.0 %, "
                "N=256/400 100 %). Seeds that succeed at every N were selected, "
                "making this a **paired design** -- different seeds per N would mix "
                "different initial-placement ensembles into the psi6(N) comparison"),
            "r_cut_fixed": 3.80,
            "r_cut_rationale": (
                "the natural r_cut grows with N, changing the truncation error 8-fold "
                "(N=64: beta*U=0.179 -> N=400: 0.011), which mixes a truncation trend "
                "into the psi6(N) slope. Fixing it at 3.80 -- the value the smallest "
                "box (N=64, L/2=4.0) permits -- raises the truncation error to 0.182 "
                "but makes it **identical at every N**, so it cannot manufacture a "
                "false slope"),
            "prod_tau": 30.0, "window": [20.0, 30.0],
            "prod_tau_rationale": (
                "since tau_relax ~ 0.098 tau_d (§8.4), 20 tau_d is 200x the "
                "relaxation time. Using 80 tau_d would push A=10 at N=400 past the "
                "budget (600 s/run)"),
            "chi2_max": CHI2_MAX,
        },
        "items": items,
        "alternatives": [
            "★ The most informative outcome is the form BREAKING (chi^2/dof > 3). "
            "That would invalidate claiming 'A=10 is not hexatic' from a single "
            "exponent, and card §8.5 would have to be weakened again. That "
            "possibility is left open.",
            "At N=64, L* = 8, so r_cut = 3.80 is 95 % of L/2 -- almost no "
            "minimum-image margin. If that point dominates the residuals it must be "
            "revisited with three points (144, 256, 400).",
            "A=10's p may come out different from the two-point run. Since r_cut, run "
            "length and seeds all change at once, a disagreement **cannot be "
            "attributed to any one axis from this run alone** -- that would require "
            "runs reverting one axis at a time.",
            "This run does not sweep the hexatic window (A = 10.03-10.75). Verifying "
            "the method comes first; only once the form is confirmed can that "
            "window's eta6 be trusted.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# S2 PREDICTION -- finite-size ladder (S17, form verification)  "
        "(auto-generated: "
        "scripts/soft2d_fss_predict.py)\n#\n"
        "# ⚠ This file is sealed before S5 runs (SEALED.sha256). Do not edit after "
        "running.\n#\n"
        + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8")
    print(f"-> {OUT.relative_to(REPO)}  ({len(items)} items)\n")
    print(f"{'item':<34} {'predicted':>10} {'tolerance':>12} {'competing':>10}")
    print("-" * 70)
    for it in items:
        cv = it.get("competing_value")
        print(f"{it['quantity']:<34} {it['value']:>10.5g} "
              f"{it['tolerance']:>12} "
              f"{(f'{cv:.5g}' if isinstance(cv, (int, float)) else '—'):>10}")
    print(f"\n  ladder {list(N_LADDER)} . lever arm ln(400/64) = "
          f"{np.log(400/64):.3f} · dof = {len(N_LADDER)-2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
