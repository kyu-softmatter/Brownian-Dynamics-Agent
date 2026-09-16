---
type: finding
author: agent
drafted: 2026-09-16
confirmed_by:
status: proposed
question: "A pre-registration states its correctness gates as percentages --
  `within 5 % of 13.375`. Is a percentage a gate?"
answer: "Only if it is wider than the statistic's own noise, and that has to be
  COMPUTED. Measured on the sedimentation campaign's sealed revision 2: of three
  gates, TWO had a tolerance inside their own 1 sigma. The worst was a 5 % band
  on a statistic with a +10.4 % bias and a 6.1 % sigma -- about 82 % likely PER
  SEED to fail on a correct run, with a sealed decision rule that reads a failed
  gate as IMPLEMENTATION FAILURE. A third number, a `hypothesis` prediction, was
  wrong by 31 sigma because it had been computed for an idealisation of the
  estimator rather than for the estimator."
cites:
  - verify/verify_sediment_design_power.py
  - campaigns/sediment_preregistration/prediction.yaml
  - knowledge/wiki/findings/an-observable-without-its-window-is-not-an-observable.md
affects_docs: [docs/02-verification.md, docs/05-pitfalls.md]
---

# A tolerance narrower than its own sigma is not a gate — it is a coin flip with a verdict attached

## The question

Pre-registration makes a tolerance load-bearing: the number is fixed before the
data exists, and a sealed decision rule turns it into a verdict. The
sedimentation campaign's sealed revision 2 said

> "CORRECTNESS GATES first, and they can end it: `l_g_fitted_tail` within 5 % of
> 13.375, `Z_dilute_tail` within 5 % of 1, `D_xy` within 10 % of 1. If any fails,
> report IMPLEMENTATION FAILURE and no physics."

Three percentages, each chosen because it *sounded* generous next to the
expected physics. None of them had been compared against the noise of the
statistic it was applied to.

## What was measured

`verify/verify_sediment_design_power.py`: draw Poisson counts on the two candidate
profiles and run **the case's own estimators** on them, 800 realisations per arm.

| gate | tolerance as sealed | bias | 1 σ | tol/σ | P(fail \| correct), per seed |
|---|---|---|---|---|---|
| `l_g_fitted_tail` | 5 % | **+10.4 %** | 6.1 % | **0.8** | **~82 %** |
| `Z_dilute_tail` | 5 % | −2.0 % | 2.9 % | 1.7 | ~7 % |
| `D_xy` | 10 % | — | not calibrated | — | — |

And one `hypothesis` prediction:

| | sealed prediction | what the estimator returns | |
|---|---|---|---|
| `l_g_fitted_dense` | 15.442 d | **17.554 d** | **31 σ** |

## The three causes, which are different from each other

**① The estimator was biased, not just noisy.** `E[ln n] ≠ ln E[n]` for Poisson
counts, so a weighted least-squares fit on `ln(counts)` has its slope flattened
where counts are small — which inflates a decay length. ⚠ Raising the count floor
makes it **worse**, because it truncates the window and collapses the lever arm:
floor 20 gives 40 bins over 22.8 d and σ = 6.1 %; floor 50 gives 15 bins over
8.7 d and σ = 30 %. There was no floor that rescued it, so the observable was
demoted to a `measurement`.

**② The prediction was computed for an idealisation of the estimator.** 15.442 d
is what a fit over the *continuous* profile on the *nominal* window returns. The
estimator never sees that: its `≥ 20 count` floor truncates the window, and a
shorter window sits where the profile is flatter, which raises the apparent decay
length by 13 %. This is
[[an-observable-without-its-window-is-not-an-observable]] one level down — the
window is not what the plan says, it is what the count floor leaves.

**③ Two of the gates could never have supported the verdict anyway.** The same
calculation gives each observable's *separation* — how many σ apart the two
competing hypotheses are in it:

| observable | separation | so it is |
|---|---|---|
| `Z_dev_wmean_vs_CS_eff` | **15.6 σ** | a discriminator |
| `l_g_fitted_dense` | **15.0 σ** | a discriminator |
| `phi_wall_measured` | **4.7 σ** | a discriminator |
| `Z_dilute_tail` | 0.4 σ | a correctness gate, and nothing else |
| `l_g_fitted_tail` | 0.1 σ | a correctness gate, and nothing else |
| `Z_dev_max_vs_CS_eff` | 0.9 σ | **nothing** — 3 σ is 9.6 %, five times the threshold |

★ A near-zero separation is **correct** for a correctness gate: the dilute tail is
ideal-gas under either hypothesis, so a statistic there cannot be contaminated by
the physics question. The point is that the two roles are distinguishable by
measurement, and a pre-registration that does not compute the separation cannot
tell which of its numbers is which.

## How to apply

- **Before sealing a tolerance, generate synthetic data from the hypothesis the
  observable is compared against, and run the real estimator on it.** Not a
  reimplementation — import the case's function, or the calculation measures the
  copy.
- **Report `tol/σ`.** Under 3, the gate will fire on correct runs. That single
  ratio is what revisions 1 and 2 were missing.
- **Report the separation too**, for every observable the verdict rests on. Under
  2 σ it cannot decide, whatever its tolerance.
- ⚠ **Generate from the MATCHING hypothesis.** The first version of this
  calculation generated from CS(φ_nominal) and compared against CS(φ_eff), read
  +6.49 %, and looked like a catastrophic estimator bias. It was the *mapping
  gap* — the thing being measured. A null has to be a null.
- **Say which tolerances the calculation did NOT calibrate.** Here `D_xy` and
  `profile_halves_chi2_nu` are sanity bounds: the Monte Carlo generates a static,
  temporally uncorrelated profile, so it says nothing about a lateral MSD or about
  frame-to-frame correlation. Naming them is what stops them looking as measured
  as the other four.

## The part worth keeping

Writing the two hypotheses down as two **profiles**, in order to generate
synthetic data from each, is what turned the campaign from calibrated into
discriminating: a stationary profile *is* an equilibrium profile, so starting
seeds at each candidate and asking which one does not move answers the question
directly. That was a by-product of doing the power calculation, not the reason for
doing it.

## See also

[[an-observable-without-its-window-is-not-an-observable]] ·
[[half-interval-errors-are-the-size-of-the-answer]] ·
[[wide-tolerance-hides-significant-deviation]] ·
[[low-seed-pilots-give-optimistic-design-power]] ·
[[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]]
