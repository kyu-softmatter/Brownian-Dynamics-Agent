# 02 · Verification — how you verify when there is no grader

Many computational fields have a **grader**. The answer is right or wrong, the
bound holds or it does not, the score reproduces. In those fields "verification"
means "ask the grader."

**There is no grader here.** No ground truth, no scalar to maximize. And the
harder problem: **the simulation always produces some number.** `g(r)` always
plots. MSD always looks like a line. Not diverging is not the same as being
right, and you cannot tell a wrong `g(r)` from a right one by eye.

So verification has to be assembled. This document is the assembly.

---

## 1 · Two axes, and the dangerous one is not the obvious one

| Axis | Asks | Where |
|---|---|---|
| **Process verification** | did the pipeline run correctly — transcription, plausibility, scales, consistency, stability, outliers | §2 (V1–V6) |
| **Result verification** | is the number that came out actually right | §3 (V7, four layers of evidence) |

They fail differently. Process fails **quietly** — it simulates the wrong system
perfectly. Result fails **plausibly**. The first is far more dangerous, because
**no amount of result verification catches it.** Four independent layers of
evidence can all agree and all be agreeing about a different system than the one
you meant.

---

## 2 · Process verification — a seven-rung ladder

Each rung catches a **different kind** of failure. They pass in order, and
verifying an upper rung while a lower one is broken is meaningless.

| | Check | Asks | Stage | On failure |
|---|---|---|---|---|
| **V1** | **Fidelity** | was the sketch / speech / text transcribed correctly | S1 → gate 1 | re-run S1 |
| **V2** | Physical plausibility | does this problem make physical sense | S1, S3 | `BLOCKED_INPUT` |
| **V3** | Scales and parameters | are the scales and parameters appropriate | S3, S4 | back to S2 / warn |
| **V4** | **Cross-stage consistency** | does anything contradict an earlier stage | pre-flight, S5 | back to S5 |
| **V5** | Numerical stability | does it survive the run | S5, S7 | repair |
| **V6** | **Outliers** | is there an outlier inside the result | S7 | investigate |
| **V7** | Physical interpretation | is the number actually right | S7, S8 | → §3 |

### V1 Fidelity — checked by back-translation

The naive check is "zero invented values." That blocks *hallucination*; it does
not check whether the input was transcribed. This failure passes it cleanly:

> Human: "500 nm **silica**" → agent: `material: polystyrene`,
> `confidence: 0.9`, `unknowns: []`
>
> **Nothing was invented. `unknowns` is empty. Every existing check passes.**

So the spec is **back-translated** into prose and put next to the original for a
human to approve. What the human reads is not YAML but **their own words handed
back** — far easier than scanning YAML for an omission, and a thing humans are
actually good at.

| Check | By |
|---|---|
| every value stated in the source appears in the spec — **zero omissions** | code |
| any spec value not in the source is `assumed: true` with a provenance | code |
| a value read off an image records **where** (`provenance: "scale bar 1 µm"`) | code |
| **does the back-translation mean the same thing as the original** | **human** |

### V4 Cross-stage consistency

Every other rung is a within-stage check. V4 is the exception, and the failure
it catches is real: **each stage passes individually and the combination is
wrong.**

| Check | Earlier | Later | Example conflict |
|---|---|---|---|
| dimension | S1 `spec.dim` | S5 box | 3D spec, 2D box |
| observation target vs run length | S1 goal | S5 `steps` | "observe gelation" but shorter than `τ_gel` |
| pre-registration vs plan | S2 | S5 | "expect crystallization" but `φ` is outside the crystalline region |
| time step | S3 `Δt/τ_D` | S5 `dt` | S5 changes `dt` and **invalidates the S3 gate** |
| box vs cutoff | S1 `N`,`φ` → `L` | S5 `r_cut` | `r_cut > L/2` — minimum image violated |
| potential vs system kind | S1 system | S5 pair | declared hard-sphere, only attraction applied |
| budget vs tier | `budget` | S5 tier | production estimate exceeds the remaining budget |

**This table getting longer is good news.** Each new row is a combination that
actually bit once, and will not bite again.

### V6 Outliers — four kinds, and aggregation destroys them

The standard diagnostics all look at **system aggregates** (mean `kT`, overall
MSD). Outliers live below that level.

| Kind | What | Method |
|---|---|---|
| **ensemble** | one seed out of five disagrees | seed-to-seed spread vs block error |
| time series | energy or pressure jumps at one frame | MAD-based robust z-score |
| particle | one particle displaces anomalously (suspect a missed neighbour-list rebuild) | tail of the displacement distribution |
| spatial | part of the box has a different density | density by sub-volume |

**Ensemble outliers matter most** — they bear directly on reproducibility, and
the other three are usually its *cause*. If one seed is anomalous and you report
the mean, **that mean is not physics; it is an average over an accident.**

This is not hypothetical here. See
[04 Cases §`trap-drag`](04-cases.md#trap-drag-2d-hex300--single-run-error-bars-were-wrong-in-both-directions):
two published conclusions were **reversed** by re-running 7 velocities × 9 seeds,
and the cause was a single run plus an underestimated error bar.

---

## 3 · Result verification — four layers of evidence

**V7.** Meaningful only for a result that cleared the six rungs above.

| Layer | Asks | Example |
|---|---|---|
| **① self-consistency** | does the simulation contradict itself | dilute tracer `D_msd` = `kT/γ` (±2%) · measured `kT` = target `kT` (±5%) |
| **② analytic limit** | go to a limit whose answer is known — do you get it | ideal gas as `φ→0` · `MSD = 2·dim·Dt` in free BD |
| **③ literature benchmark** | does a value someone else measured reproduce | Carnahan–Starling `Z(φ)` · hard-sphere `φ_freeze = 0.494` |
| **④ independent method** | does a different route give the same answer | BD vs HPMC · the same `B₂*` from different potentials |

> **Running the same code twice is one piece of evidence.** So is re-running with
> a different seed. Evidence has to differ in **kind**.

### Evidence grade — every result wears a badge

| Grade | Condition |
|---|---|
| `certified` | ③ literature + ② analytic limit + ④ independent method |
| `verified` | ① self-consistency + one of ②–④ |
| `plausible` | ① self-consistency only |
| `unverified` | nothing to compare against |

**Do not hide `unverified`.** Simulating a new system that is in no benchmark is
legitimate; calling it "verified" is not. The badge goes in the report.

### Disagreement protocol

| Situation | Action |
|---|---|
| two layers disagree | **it is a bug. Stop.** Do not proceed until it is found |
| all agree | trust it, proceed |
| all disagree | the definitions differ. Re-check **the non-dimensionalization and the observable definitions** first |

### Literature is verification infrastructure

Literature plays the grader role we do not have. So collecting papers is not
a side activity — **it is the verification infrastructure.** A paper used only
as reading context throws away half its value; extracted into a machine-readable
benchmark and run as a regression test, it becomes a grader:

```yaml
# knowledge/wiki/benchmarks/benchmarks.yaml
- id: carnahan_starling_hs_eos
  system: hard_sphere_3d
  input: {phi: 0.30}
  observable: compressibility_factor_Z
  expected: 4.577            # Z = (1+φ+φ²−φ³)/(1−φ)³
  tolerance_rel: 0.02
  cost: cheap
  evidence_layer: 3
  source: "Carnahan & Starling, J. Chem. Phys. 51, 635 (1969)"
```

⚠️ **The literature base here is narrow**: 38 of 42 distillations are the
group's own published work ([NOTICE §1](../NOTICE.md)). Any sentence of the form
*"the literature says…"* in this repository should be read against that.

---

## 4 · The role of a prediction decides what a mismatch means

This is the piece `simulation_auto` added, and it is the one that keeps
verification from destroying discovery. Every comparison carries a **role**
([`bdbot/metrics.py`](../bdbot/metrics.py) `ROLES`):

| Role | Where the prediction comes from | A mismatch means |
|---|---|---|
| `implementation_check` | derived **from the model I implemented** | **a bug** → fix it |
| `hypothesis` | an assumption the simulation does **not** impose — continuum, dilute limit, effective medium, literature | **a result** → report it |
| `measurement` | no prediction | the simulation is the answer |

Default is `measurement` (no verdict). Only `implementation_check` mismatches
are FAIL.

**Design consequence.** When specifying a case, write down separately *what the
simulation imposes* and *what the theory adds*. If the second list is empty, the
case can validate but can never discover. `abp-rod-2d-run-flip` was exactly
that: five predictions, all `implementation_check`, **zero hypotheses** — and it
produced zero findings, as designed and as measured.

Conversely `chain-bend-2d-dlvo` was built so its central prediction was a
`hypothesis` — and reporting the mismatch *was* the result.

---

## 5 · Three verdicts, not two

| Verdict | Meaning | Recorded where |
|---|---|---|
| `PASS` | inside tolerance, proceed | journal (automatic) |
| **`PASS-with-doubt`** | **inside tolerance but suspicious** | journal + verification ledger |
| `FAIL` | outside tolerance | journal + repair |

**The middle one is the point.** With only PASS/FAIL, the suspicious has nowhere
to go, so it disappears — it passed, so nobody looks again, and later when the
result is strange there is no way to retrace *what was marginal*.

Attach `PASS-with-doubt` when: the tolerance was passed **narrowly** (1.9 % on a
±2 % criterion) · the threshold's own basis is weak (convention rather than
literature) · it passed but **differed from expectation** · the check was
**skipped** for want of data or cost.

### Doubt lowers the evidence grade

Without an incentive nobody records doubt. So they are wired together:

| Condition | Grade ceiling |
|---|---|
| any `PASS-with-doubt` | **`certified` unreachable** — ceiling is `verified` |
| doubt in V1 (fidelity) or V4 (consistency) | ceiling is `plausible` — if the *system* is in doubt the numbers are meaningless |
| any unresolved `FAIL` | `unverified` |

For this to be an incentive to record doubt honestly rather than hide it, a low
grade has to be treated as **an accurate description of state, not a
punishment**. `plausible` is not a bad result; it is an honest one.

### A reason must be a cause, not a symptom

| ✗ symptom | ✓ cause |
|---|---|
| "it diverged" | "WCA's `r⁻¹³` core makes `F·dt/γ` blow up in the overdamped scheme" |
| "seed 3 is weird" | "only seed 3 kept an overlap in the initial placement, so energy spiked in the first 100 steps" |
| "`kT` is off" | "`γ` went in as SI while everything else was reduced — units mixed" |

If you cannot state the cause, **record that fact**: `cause: unknown` + the
observation + the next candidate to investigate. Three `cause: unknown` entries
on one theme is a signal that the theme is worth investigating.

---

## 6 · The failure mode this document exists to prevent

A checker that is not wired up cannot be wrong out loud. Two measured instances:

**`step_health()` never ran — in all 81 runs.**
[`bdbot/health.py`](../bdbot/health.py) documents `step_health()` as "the core of
this module." A name mismatch meant it never executed: `run.Guard` computed
`dt·|F|max` into `l4`, while the health tool looked for
`numerics["step_rms_sigma"]`. Finding nothing, it printed `"not measured"` and
**returned HEALTHY anyway**. `82/82 HEALTHY` read like coverage. Three lessons:

1. **Silence is not success.** The count of *unmeasured* runs is now printed
   separately, and that HEALTHY is explicitly labelled as "no divergence, no
   stall, no collapse" and **not** "`dt` is small enough."
2. **Measure from forces.** `dt·|F|max/γ` *is* `dt/τ_fast` with nothing to
   subtract. A position-difference measure has to subtract thermal noise
   `√(2·dim·dt)` in quadrature, and when drift is 0.5 % of thermal noise, finite
   sampling **clips the drift to zero** — measured.
3. **Keep the worst value.** Only the last guard sample used to be stored.
   Measured: peak force 1062.9 vs last 244.2 kT/σ — **4.4×**. Stability is
   decided by the worst case, not the mean.

**The pre-run gate rejected 80 of 83 specs, none of them real failures.**
`health.gate()` tested `verdict != "PASS"`, so it refused `"PASS (3 warnings)"`.
Statistics and finite-size issues are defined as ⚠ (not ❌), and `run.execute`
was reading `startswith("FAIL")` correctly. **The two disagreed and nobody knew,
because `execute` never called `gate()`** — an unwired checker can be wrong
indefinitely. The gate now blocks on three things only: hash mismatch, `FAIL`,
and an L3 integrity **error**. Warnings and thin margins are **shown but not
blocking** — and they are always printed, since a gate that passes silently is
not a gate.

⛔ **And the coverage could not be backfilled.** Both routes were measured
closed: re-running produces a *different* `run_id` because it is content-
addressed (all 79 unmeasured runs failed to reproduce), and replaying forces
from the stored trajectory is **invalid under time-dependent driving** — for
`trap-drag`, replay gave `|F|max = 17041` against the guard's true 1062.9, a
**16.03×** overestimate, because the trap anchor is pinned at `t=0`. So the
coverage denominator is *"runs executed after the measurement was wired,"* and
legacy runs are marked **not retroactively knowable**.
