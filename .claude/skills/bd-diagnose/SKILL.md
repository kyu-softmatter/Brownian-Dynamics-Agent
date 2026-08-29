---
name: bd-diagnose
description: Diagnoses a broken or suspicious Brownian Dynamics run. Use when the simulation died with NaN, when a result differs from the literature or an analytic solution by orders of magnitude, when a prediction FAILed, or when an MSD or RDF has the wrong shape. Follows the diagnostic paths in knowledge/wiki/findings in order — it gets stronger as failures accumulate. To design and run something new use bd-pipeline; to record a diagnosis as knowledge use bd-knowledge.
---

# bd-diagnose — finding the cause of a broken run

> **A wrong `FAIL` is as bad as a wrong `PASS`.** It sends you chasing physics
> that is not there. So **suspect the analysis code first.**

## 0. Five categories of cause — decide which one before anything else

| Category | Meaning | Frequency (measured, first end-to-end run) |
|---|---|---|
| **`analysis`** | the measurement or statistics code is wrong | **1/4** ← most dangerous |
| `numerical` | `dt` too large, integrator, insufficient convergence | 1/4 |
| `interpretation` | the S1 reading was wrong (dimension, units, arrows) | 0/4 |
| `modeling` | the potential or approximation is inappropriate | 0/4 |
| `environment` | package, parser, platform | 2/4 |

**Suspect physics last.** Of the four failures on the first end-to-end run,
**zero** were physics.

⚠️ **First, check the prediction's `role`.** A `hypothesis`-role mismatch is not
a defect to be diagnosed away — it is the result, and running this elimination on
it will produce a spurious "cause." Only `implementation_check` mismatches are
bugs. See [docs/02](../../../docs/02-verification.md#4--the-role-of-a-prediction-decides-what-a-mismatch-means).

---

## 1. Elimination order — follow it as written

### ① Does the statistic fluctuate?

```python
from simbot.guards import assert_statistic_fluctuates
assert_statistic_fluctuates(samples, "name")
```

**A "measurement" that does not fluctuate is an arithmetic identity.** This
happened on 2026-07-28: subtracting the mean from displacements and then
measuring the cross-correlation gives `cross/auto = −1/(n−1)`
**identically**, and the standard deviation over 200 repetitions was `6.7e-20`.
**The result was plausible enough that it nearly passed.**

### ② Do two independent routes agree? (self-consistency)

Measure the same quantity two ways and compare:

| Quantity | Route A | Route B | Must agree |
|---|---|---|---|
| `⟨x²⟩` | snapshot variance | MSD plateau / 2d | `1.0` |
| `D` | MSD slope | `kT/γ` | at short times only |

If they disagree it is **the analysis code**, not the physics.

### ③ Are the samples independent?

**A KS or χ² test on correlated samples always rejects.** And it rejects *more*
confidently as `n` grows, so the intuition "many samples, therefore trustworthy"
works **exactly backwards**.

Full account:
[`ks-test-needs-independent-samples`](../../../knowledge/wiki/findings/ks-test-needs-independent-samples.md)

The check: is the frame interval longer than the correlation time (`~2τ`)?

```python
frames_per_tau = 1.0 / (frame_interval_steps * dt_star)
step = ceil(2.0 * frames_per_tau)      # use only frames spaced this far apart
```

### ④ Units and dimensions

| Suspicion | The check |
|---|---|
| a diameter went into `γ = 6πηa` | the timescale is wrong by exactly 2× |
| the timescale was taken as `τ_D` | `roundtrip_errors` exceeds `1e-3` |
| `kT` and `ε` were both set to 1 | `T* = 1` has been pinned |
| 3D Stokes used in 2D | if intended, say so; otherwise a bug |

```bash
<PY> -c "
from simbot.spec import SystemSpec
from simbot.nondim import roundtrip_errors
print(roundtrip_errors(SystemSpec.load('<spec.yaml>')))"
```

### ⑤ Numerics — `dt` and the guards

```bash
<PY> cli.py converge <spec.yaml>
```

If the answer changes at `dt/2`, `dt` is too large. If it does not change, that
means **this error bar cannot see the difference** — it is not proof that `dt` is
fine.

Guard results: the run manifest's `batch`, and `guards` in each replica manifest.

### ⑥ And only then suspect the physics

- Is this a regime where neglecting hydrodynamic interactions breaks down?
- Is constant-`γ` Stokes insufficient? At `a ≳ 5 μm`, Basset is `4.3 %` and
  Faxén is `2.3×`
- Is the overdamped approximation valid (`τ_i/τ_process ≪ 1`)?

---

## 2. First suspicion by symptom — query the knowledge base

```bash
ls knowledge/wiki/findings/
grep -rl "<symptom keyword>" knowledge/wiki/
```

| Symptom | Look at first |
|---|---|
| NaN / blow-up | `dt` too large, initial overlap. Which frame did the finiteness guard catch? |
| per-step displacement exceeds `√3σ_step` | **impossible** — HOOMD's noise is uniform. Suspect the measurement code |
| distribution tails differ from Gaussian | **this is normal.** `findings/hoomd-brownian-scheme-and-noise.md` |
| kurtosis is exactly 3.000 | **that is the suspicious case.** It should come out as `3 − 1.2 dt*` |
| `⟨x²⟩` is high by `dt*/2` against the prediction | **normal.** Euler–Maruyama bias |
| KS p = 0.0000 | correlated samples. Go to ③ |
| the `dt` gate passed but something is wrong | the displacement gate is 1086× loose in a trap system |
| cost estimate is strange, or a TypeError | YAML 1.1 exponent notation (`5e-3` is a string) |
| KeyError trying to run a system with no card | **that is correct behaviour.** Make the card first |
| a run exists but `status` reports zero | `result.txt` is written by the case script, not the engine |

## 3. Diagnostic tools

```bash
# seal status
<PY> -c "
from simbot.io import RunDir, verify_seal
print(verify_seal(RunDir('runs_s1s8/<id>')).summary())"

# re-check the gates
<PY> -c "
from simbot.spec import SystemSpec, validate
r = validate(SystemSpec.load('runs_s1s8/<id>/03_spec.yaml'))
print(r.table()); print(r.problems)"

# look at the trajectory directly
<PY> -c "
from simbot.analysis.trap import load_run
d = load_run('runs_s1s8/<id>/raw/<label>')
print({k: (v.shape, v.dtype) for k, v in d.items()})
print('var range', d['indep_var'].min(), d['indep_var'].max())"

# regenerate only the figures, without re-running
<PY> cli.py resume runs_s1s8/<id>

# the L4 numerical-health verdict on a bdbot run
<PY> -m bdbot.cli health runs/<run_id>
```

## 4. When the diagnosis is done — record it

**A failure is a deliverable.** Delete it and half the evidence about whether
this agent works at all disappears.

```
knowledge/wiki/findings/<slug>.md              you found the cause
knowledge/wiki/findings/dead-end-<slug>.md     this route is blocked
```

A `dead-end` records the **cause** in `why_it_failed`. "It diverged" is a
symptom; **"WCA's `r⁻¹³` core made `F·dt/γ` explode under overdamped dynamics"**
is a cause. If you cannot state one, write `cause: unknown` plus the observation
and the next candidate — three of those on one theme is a signal worth
investigating.

**If "prevention: none" keeps getting written, that is the signal to build a
guard.**

## 5. What you do not do

- Do not say "fixed" without finding the cause. Fix `INCONCLUSIVE` as a fact
- A symptom disappearing is not a cause disappearing — if you cannot write down
  what was fixed and why, it is not fixed
- Do not "solve" a problem by turning off a gate. When you do turn one off,
  **the reason has to be in the card**
