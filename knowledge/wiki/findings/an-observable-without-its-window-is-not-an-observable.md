---
type: finding
author: agent
drafted: 2026-09-16
confirmed_by:
status: proposed
question: "A decay length fitted from a profile was going to be an
  `implementation_check` against the imposed value. Is it one?"
answer: "No -- not without the fit window, and the window changes the answer by
  more than the tolerance. The CS hydrostatic profile at l_g = 13.375 d fits
  15.442 d over the dense window [2, 80] (+15.5 %) and 13.50 d over the tail
  [55, 100] (+0.9 %). A 5 % check against 13.375 therefore FAILS on a correct run
  in the first window and passes in the second. The role and the window are one
  decision, not two."
cites:
  - cases/sediment_3d.py
  - tests/test_sediment_eos.py
  - campaigns/sediment_preregistration/prediction.yaml
affects_docs: [docs/02-verification.md]
---

# An observable without its window is not an observable — it is a rule 7′ error waiting

## The question

`sediment-pmma-3d` imposes `l_g/d = 13.375` as a uniform force `F_z = −kT/l_g`
and fits `l_g` back out of `ln n(h)`. The simulation imposes the number, so
recovering it looks like the definition of an `implementation_check`: a mismatch
is a bug (CLAUDE.md rule 7′).

## What was measured

It is an `implementation_check` **only in the dilute tail.** Integrating the
hydrostatic equation `d(φZ)/dz = −φ/l_g` with `Z = Z_CS` — i.e. the profile the
interacting system actually has — and fitting the result:

| fit window `h/d` | apparent `l_g/d` | vs the imposed 13.375 | a 5 % check |
|---|---|---|---|
| `[2, 80]` (dense) | **15.442** | **+15.5 %** | **FAILS** |
| `[55, 100]` (tail) | **13.50** | +0.9 % | passes |

The cause is not subtle once stated: `Z > 1` extends the column, so the profile
is **not** `exp(−z/l_g)` where the excluded volume matters. `l_g` is exactly
recoverable where `φ → 0` and nowhere else. The non-interacting capability check
reached `−0.014 %` (0.02 σ) at this `l_g` — that number is real, and it is a
measurement of the *point-particle* limit.

So the original single observable was a check that **a correct run would have
failed**, which is the failure rule 7′ exists to name: calling a physical result
a bug. It was caught by integrating the ODE before sealing, not by a run.

## What the fix looks like

One estimator, two windows, two roles:

| observable | window | role | prediction | source of the prediction |
|---|---|---|---|---|
| `l_g_fitted_tail` | `[55, 100] d` | `implementation_check` | 13.375, 5 % | the imposed force |
| `l_g_fitted_dense` | `[2, 6 l_g] d` | `hypothesis` | 15.442, 8 % | the CS hydrostatic ODE |

★ The second one is *more* informative than the first, because the dense-window
decay length and the weighted-mean `Z` deviation are two independent reductions
of one profile and must agree in sign. A single check against 13.375 would have
thrown that away and called it a bug.

## How to apply

- **Write the window next to the role, in the same line of the prediction.** If
  the role is `implementation_check`, the question to answer out loud is *over
  what range is the imposed value the correct prediction* — and if the answer is
  "not all of it", the observable splits.
- **Get the prediction for a non-trivial window by integrating the model, not by
  reasoning about it.** The `+15.5 %` here is not estimable by eye, and the wrong
  intuition ("excluded volume is a small correction at `φ = 0.017`") is exactly
  what makes it dangerous: the *mean* `φ` is 0.017, the wall value is 0.10.
- ⚠ **A tolerance cannot rescue a wrong prediction.** Widening the 5 % to 20 %
  would have made the dense-window check pass — and made it unable to detect
  anything, which is the sibling failure
  [[wide-tolerance-hides-significant-deviation]].

## See also

[[half-interval-errors-are-the-size-of-the-answer]] ·
[[wide-tolerance-hides-significant-deviation]] ·
`knowledge/wiki/systems/passive-sphere--sedimentation.md`
