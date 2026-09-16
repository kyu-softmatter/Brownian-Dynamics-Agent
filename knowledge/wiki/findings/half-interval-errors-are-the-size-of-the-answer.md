---
type: finding
author: agent
drafted: 2026-09-16
confirmed_by:
status: proposed
question: "When a quantity is obtained by integrating or inverting a binned
  profile, where does the error come from -- and how big is it compared with the
  threshold the answer is judged against?"
answer: "From half-interval truncation, three times in one case, and each time it
  was the same size as the whole decision band. Reading an EOS from a 0.5 d
  histogram with a left-Riemann sum was +1.65 % against a 2 % threshold; the
  trapezoid rule leaves (bin/l_g)^2/8 = 0.017 %. Sampling from the same profile
  with a cumsum CDF left the mean 0.44 % low with KS ON the 95 % critical value.
  The fix is the same in all three: integrate with the trapezoid rule, and derive
  the residual as a LAW so the test fails when the bin changes."
cites:
  - cases/sediment_3d.py
  - tests/test_sediment_eos.py
  - knowledge/wiki/systems/passive-sphere--sedimentation.md
affects_docs: [docs/05-pitfalls.md]
---

# Half-interval errors are the size of the answer

## The question

`sediment-pmma-3d` reads the equation of state out of a density profile by the
hydrostatic relation

```
Z(z) = (1 / (l_g phi(z))) * integral_z^Lz phi(z') dz'
```

and compares it against Carnahan–Starling with a **2 % band** inherited from
[[wca-reproduces-carnahan-starling]]. The profile is a histogram with
`bin = 0.5 d` and `l_g = 13.375 d`. So the integral is a sum over bins, and the
only question is where each bin's contribution is taken from.

## What was measured

**① the EOS integral — a left-Riemann sum is `+1.65 %`**

`np.cumsum(phi[::-1])[::-1]` integrates from each bin's **lower edge**, which
adds the whole bin containing `z` instead of half of it. Against the exact
finite-cell ideal gas `Z = 1 − exp(−(L_z − z)/l_g)`:

| what the integral takes | error at `bin/l_g = 0.0374` |
|---|---|
| the whole bin containing `z` (cumsum) | **+1.65 %** = `bin/(2 l_g)` |
| **half** of it, plus the whole bins above | **+0.017 %** = `(bin/l_g)² / 8` |

1.65 % against a 2 % band. It would not have crashed, would not have looked
wrong, and would have been reported as a deviation from Carnahan–Starling.

**② the initial-condition sampler — a cumsum CDF is `−0.44 %` on the mean**

The same error class, one layer up. Building the CDF of `phi(z)` with
`np.cumsum(w)` rather than `np.cumsum(0.5*(w[1:]+w[:-1])*dz)` drops the
half-interval at every step:

| CDF | KS `D` at `N = 200,000` | 95 % critical | mean |
|---|---|---|---|
| cumsum | 0.00291 | 0.00304 | **0.44 % low** |
| trapezoid | 0.00238 | 0.00304 | — |
| trapezoid, `N = 1e6` | **0.000961** (p = 0.31) | 0.00136 | — |

★ The cumsum version **passed** its KS test — at 96 % of the critical value.
Marginal is not correct, and at `N = 1e6` it would have failed. A test whose
verdict depends on how many samples you happened to draw is not a test.

**③ the decay-length fit — not half-interval, but the same shape of mistake**

`l_g` fitted over a window where `Z > 1` is **not** `l_g`: integrating the CS
hydrostatic ODE over `h ∈ [2, 80] d` gives an apparent `15.442 d` against an
imposed `13.375`, `+15.5 %`. Over the tail `[55, 100] d` it gives `13.50 d`,
`+0.9 %`. The number is window-dependent, and the window was not part of the
claim. → [[an-observable-without-its-window-is-not-an-observable]]

## Why this is a finding and not a bug report

Three independent half-interval errors in one case, each the size of the answer.
The common cause is not carelessness — each individual expression is the obvious
one. It is that **a binned profile has no single value at a point**, so every
reduction of it has to say which part of the bin it used, and `cumsum` is the
expression that silently answers "all of it".

## How to apply

- **Integrate a histogram with the trapezoid rule.** If a `cumsum` appears in an
  integral over bins, the half-interval is already wrong; the only question is
  whether it is large enough to notice.
- **Assert the residual as a LAW, not a number.** `tests/test_sediment_eos.py`
  asserts the error is `(bin/l_g)²/8` across three decades of `bin/l_g`. A test
  pinned to `0.017 %` at the production bin would pass with a first-order
  integrator at a smaller bin and say nothing.
- ⚠ **A KS test sitting at 96 % of its critical value is a failure you have not
  paid for yet.** Re-run it at 5× the samples before believing it.
- **State the window with the observable.** An `l_g` with no window is not a
  measurement of `l_g`.

## See also

[[wca-reproduces-carnahan-starling]] ·
`knowledge/wiki/systems/passive-sphere--sedimentation.md` ·
[[ks-test-needs-independent-samples]]
