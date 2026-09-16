---
type: finding
author: agent
drafted: 2026-09-16
confirmed_by:
status: proposed
question: "A Monte Carlo that drew every frame as independent said the primary
  statistic's sigma was 0.38 %. The measured integrated autocorrelation time of
  the mean particle height was 86-132 frames out of 478, i.e. N_eff = 3.6-5.5.
  Does that 10x inflation apply to the primary statistic?"
answer: "NO, and acting on it would have killed an affordable campaign. `h_mean`
  tracks the SLOW global drift, whose timescale is tau_fall = L_z l_g/D_0 =
  9.35 tau_sed; the dense-region statistics that carry the verdict decorrelate
  far faster. The measured block SEM of the primary statistic is 0.399 % against
  the Monte Carlo's 0.384 % -- a ratio of 1.04, not 10. ⚠ But the SAME campaign
  then measured a 5.4 % seed-to-seed spread on that statistic where the
  within-run block SEM was 0.4-0.85 %, so the block SEM understates the error
  too. Neither shortcut works: the error bar that matters is the one over
  independent seeds."
cites:
  - cases/sediment_3d.py
  - verify/verify_sediment_design_power.py
  - knowledge/wiki/systems/passive-sphere--sedimentation.md
affects_docs: [docs/02-verification.md]
---

# A correlation time measured on one observable is not a correlation time for another

## What happened

Three estimates of the same uncertainty, in the order they were obtained, all for
`Z_dev_wmean_vs_CS_eff` on a 4 `τ_sed` sedimentation run:

| | estimate | what it assumes |
|---|---|---|
| ① Monte Carlo on a static profile | **0.384 %** | every frame is an independent Poisson draw |
| ② `√τ_int` inflation of ①, with `τ_int` from `h_mean` | **~3.8 %** | the mean height's correlation time is the statistic's |
| ③ block SEM over 16 blocks of the run | **0.399 %** | the blocks are independent |
| ④ spread over independent seeds | **~5.4 %** (2 seeds) | nothing — the seeds really are independent |

★ ② was **wrong**, and it was nearly acted on. Applied to the campaign it turned
15.6 σ of discrimination into 1.7 σ, which would have made the experiment
unaffordable and the conclusion INCONCLUSIVE by construction.

The reason is physical and was already written in this pair's own card:
`τ_fall = L_z l_g/D_0 = 9.35 τ_sed` is the drift time across the cell, nine times
the slowest *equilibrium* mode. `⟨h⟩` is exactly the coordinate that relaxes on
`τ_fall`, so its autocorrelation measures the slowest thing in the problem. The
equation of state is read from the **dense region near the wall**, which
re-equilibrates locally and much faster.

⚠ And ③ is not safe either. The measured `τ_int` of 86–132 frames exceeds the
block length of 478/16 ≈ 30 frames, so the blocks are correlated and the block
SEM is biased low — which ④ then demonstrated directly, at about 6×.

## How to apply

- **Do not transfer a `τ_int` between observables.** Measure it on the observable
  you are putting an error bar on, or use a method that does not need one.
- **Block-average with `block_length > τ_int` or say that you have not.** Sixteen
  blocks looks like `ν = 15`; at `τ_int ≈ 3 ×` the block length it is worth
  nearer 5. Record the number of blocks *and* the correlation time beside every
  block SEM, so the reader can see which regime it is in.
- **The seed-to-seed spread is the error bar that matters**, because independent
  seeds are the only thing here that is actually independent. Report it next to
  the within-run estimate, and when they disagree by 6× report the disagreement
  rather than picking one (`A1`: *"if two pieces of evidence disagree, stop"*).
- **Ask which mode an observable couples to before trusting its noise estimate.**
  In an external field there is more than one relaxation time, and the card's
  gate table is where they are written down.
- ⚠ **A statistic that is cheap to compute from a single run is not cheap to put
  an error bar on.** This is the concrete cost of `A2`, and it is the reason the
  seed count in a campaign cannot be chosen from a static-profile calculation.

## What this cost, and what it bought

Two production runs and one pilot, and a campaign stopped three times. It bought
per-block accumulation in the case, so every profile statistic now carries an
error bar measured from its own run — and the knowledge that that error bar is a
lower bound until the seed-to-seed spread confirms it.

## See also

[[a-tolerance-narrower-than-its-own-sigma-is-not-a-gate]] ·
[[half-interval-errors-are-the-size-of-the-answer]] ·
[[ks-test-needs-independent-samples]] ·
[[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]] ·
`knowledge/wiki/systems/passive-sphere--sedimentation.md`
