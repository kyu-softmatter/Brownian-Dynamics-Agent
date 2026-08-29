# S6 · S7 — figures · validation

> **The engine does not decide. You propose; a human confirms.**

## S6 — figures

`cli.py run` produces them automatically. Your job is to **read them, and dig in
when something looks wrong.**

### The caption convention — enforced by code

`simbot.viz.FigureSet.save()` **requires** two fields:

| Field | What |
|---|---|
| `caption` | the text under the figure. **Must contain the numbers and the verdict** |
| `shows` | **"what is this figure trying to show"** — different from the caption |

If `caption` is "MSD curve", then `shows` is "whether the MSD follows a single
exponential". An empty string raises — this is **blocked at creation**, not
checked afterwards.

### A figure that was not drawn leaves a reason

`FigureSet.skip(name, reason)`. You cannot skip without one. It appears in the
final table of the figures document.

★ Distinguish the *kinds* of reason:
- **not defined** — with no pair interaction, `g(r)` carries no structural
  information
- **not implemented** — there is a pair interaction, but the structure analysis
  does not exist yet
- **environment missing** — a 3D snapshot needs `fresnel` and it is not installed

If all three render as the same "no figure", you cannot tell what to fix.

### Dual axes

Every time and length axis carries **both the dimensionless and the SI value**
(`viz.add_si_axis`). With only `t/τ_trap` you cannot answer "how many ms"; with
only `t [ms]` you cannot compare against another paper.

### Text inside a figure is English

matplotlib's default font has no Hangul glyphs, so they render as tofu (□). And
the fonts that *do* have Hangul are missing `−` (U+2212) and `ŷ` (U+0177) —
measured. **Do not fix it by switching fonts.** Write axes, legends, titles and
annotations in English from the start, and confirm zero `missing from font`
warnings. `test_s6_viz.py` checks the source.

### The correlated-sample trap — it applies to figures too

A position-distribution figure uses **independent frames only**
(`decorrelation_tau=2.0`). Using every frame inflates the sample count, which
falsely shrinks the error bar on the kurtosis and puts a number in the caption
that differs from the verified value in `metrics.json`.

The same trap was hit in a KS test:
[`ks-test-needs-independent-samples`](../../../../knowledge/wiki/findings/ks-test-needs-independent-samples.md)

### Animations

Put `kT=0` (the mode shape) and `kT>0` (how far thermal noise buries it)
**side by side**. That shows an SNR problem faster than any number. A cheap
large-`dt` animation is fine, but **it must be labelled as not being the
production measurement.**

---

## S7 — validation

### 0. The seal comes first

`validate_run` verifies the seal, and **if it is broken it does not build the
comparison table.**

```
★ seal violation — modified ['02_prediction.md']. The prediction may have been
  edited after the run, so no comparison table is produced.
```

What to do then: **revert the prediction, or start a new run.** Do not route
around it.

### 1. Three verdicts

| verdict | When |
|---|---|
| `PASS` | the deviation is inside the tolerance band **and** a verdict was possible |
| `FAIL` | the deviation is outside the band |
| **`INCONCLUSIVE`** | ① `SE > tolerance half-width` → this measurement cannot decide this band<br>② design power `< 1σ` → indistinguishable from the competing hypothesis<br>③ there is no statistical error → a number without an error bar cannot enter a conclusion |

### 2. `INCONCLUSIVE` is not a failure

**It is a fact.** "This measurement cannot decide" is a result, and it comes with
the sample multiple that would be required.

On the first end-to-end run, 2 of 9 were `INCONCLUSIVE` and **both were foreseen
in the prediction document.** Writing them up as `PASS` would be luck, not
verification.

What you judge: **was it foreseen, or is it a design mistake?**

| | Response |
|---|---|
| the prediction document says "INCONCLUSIVE expected" | fix it as a fact and move on. State that the conclusion does not depend on it |
| it was not foreseen | decide whether the tolerance was unrealistically tight or the seeds were too few |
| the conclusion depends on this item | raise the seeds (`escalate_to: 8`) or redesign into conditions where power exists |

> ★ **Do not demand a `3σ` rejection where the power cannot produce `3σ`** — that
> is an unachievable assertion. In that regime, fix `INCONCLUSIVE` as a fact.

### 3. `PASS ⚑` — what a wide tolerance hid

A `PASS` whose deviation exceeds `3σ` gets flagged by the judge.

**Two causes produce the same symptom, and you need to know which:**

| Cause | Prescription |
|---|---|
| the tolerance is far wider than the statistical precision | narrow it (or write down why it is wide) |
| **a known bias was left out of the prediction** | fix the prediction |

Real case: predicting the MSD plateau as exactly `2d` gave `3.54σ`. Since
`plateau = 2d⟨x*²⟩`, the EM bias multiplies in. Full account:
[`wide-tolerance-hides-significant-deviation`](../../../../knowledge/wiki/findings/wide-tolerance-hides-significant-deviation.md)

### 4. A `FAIL` requires a cause category

```python
validate_run(pred, meas, rundir=rd, causes={"var_x": "numerical"})
```

| Category | Meaning |
|---|---|
| `numerical` | `dt` too large, integrator, insufficient convergence |
| `modeling` | the potential, the approximation, or neglecting HI is inappropriate |
| `interpretation` | the S1 reading was wrong (arrows, dimension, units) |
| `analysis` | the measurement code is wrong ← **the most dangerous** |
| `environment` | package, parser, platform |

★ **Suspect `analysis` first.** On the first end-to-end run a `FAIL` came out
that was not physics at all — a KS test had been applied to correlated samples.
**A wrong FAIL is as bad as a wrong PASS** — it sends you chasing physics that is
not there.

⚠️ **But check the `role` before starting.** A `hypothesis`-role mismatch is not
a defect to be eliminated — it is the result, and running the elimination on it
will manufacture a spurious "cause." Only `implementation_check` mismatches are
bugs. → [docs/02](../../../../docs/02-verification.md#4--the-role-of-a-prediction-decides-what-a-mismatch-means)

Elimination order:
1. Does the self-consistency check pass (do two independent routes agree)?
2. Does the statistic fluctuate (`guards.assert_statistic_fluctuates`)?
3. Are the samples independent (were correlated frames treated as iid)?
4. Are the units and dimensions right?
5. And only then suspect the physics

### 5. The strongest check that was not in the prediction — self-consistency

An item like `plateau_over_2d_var` is **the ratio of two independent routes**
(time series vs snapshot). A departure from `1.0` is **the analysis code**, not
the physics.

On the first end-to-end run it was `0.99961` — agreement to 0.04 %. **That is the
strongest evidence available for ruling out an analysis error.** Put one such item
into every new system.

### 6. Do not delete the footnotes

`ValidationReport.notes()` keeps the footnotes on `PASS` items in the report too.
This is where remarks like "this is not an independent check" live, and dropping
them inflates the conclusion.

Example: in a pure harmonic trap `kT_conf_star` is **algebraically identical** to
`⟨x²⟩`. A stationarity figure **shows** that with a residual of `4.4e-16` — a
measurement, not a claim.

### 7. The verdict is a proposal

```yaml
verdict_overall: PASS_WITH_INCONCLUSIVE
proposed_by: agent
confirmed_by: null            # <- awaiting human confirmation
```

**No code path exists that fills `confirmed_by`** (pinned by a test). You do not
fill it either. A verdict with `confirmed_by: null` does not enter the benchmark
ledger's aggregate.

### 8. What you may and may not say

| ✅ | ❌ |
|---|---|
| "`⟨x²⟩ = 416.58 ± 1.85 nm²`, agreeing with the prediction `414.19` at `1.29σ`" | "it matches well" |
| "It looks like a PASS — on the basis of …, though … concerns me" | "it is verified" |
| "It fluctuates without a trend (see the stationarity figure)" | "it reached equilibrium" |
| "The simulation is accurate about itself" | "it agrees with real water" |

The last row matters most. On the first end-to-end run
`τ_fit/τ = 0.9998 ± 0.0008` (0.08 %), but **the correspondence error between the
model and real water is `~5 %`** (Basset, Faxén, gravity).
⇒ **You cannot claim 0.6 % precision against an experiment.**
Full account:
[`stokes-drag-corrections`](../../../../knowledge/wiki/concepts/stokes-drag-corrections.md)
