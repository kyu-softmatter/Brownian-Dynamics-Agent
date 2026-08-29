# S2 — the prediction protocol

> **Write the answer down before the simulation runs.** This stage is what
> underwrites the project's scientific honesty. A prediction written afterwards is
> not a prediction; it is an explanation.

## 0. Get numbers by calling

```python
from simbot.spec import SystemSpec, derive
from simbot.estimators import harmonic_trap, euler_maruyama_trap_variance_bias
from simbot.analysis.trap import em_uniform_noise_excess_kurtosis
```

**No hand calculation.** Not even `4.14e-21 × 2` — code has no arithmetic error
and leaves a record.

Generating the prediction file *with code* is safest. The example was made that
way:
[`examples/trap-2d-5um/prediction.yaml`](../../../../examples/trap-2d-5um/prediction.yaml)

## 1. The four parts of one prediction

```yaml
- quantity: var_x_star            # must match the measurement name **exactly**
  value: 1.0025062656641603       # core function output. Do not truncate
  tolerance: ±1%                  # outside this is a FAIL
  basis: 'Euler-Maruyama stationary variance 1/(1-dt*/2). estimators.euler_maruyama_...'
  role: implementation_check       # ★ what a mismatch MEANS
  discriminates: 'is the integrator scheme EM or exact'
  competing_value: 1.0            # ★ the competing hypothesis — used for the power calculation
```

If `quantity` differs from the measurement name, `validate_run` reports **"no
corresponding measurement"**.

### ★ `role` decides what a mismatch means

| Role | Where the prediction comes from | A mismatch means |
|---|---|---|
| `implementation_check` | derived from the model **we** implemented | **a bug** → fix it |
| `hypothesis` | an assumption the simulation does **not** impose — continuum, dilute limit, effective medium, a paper | **a result** → report it |
| `measurement` | no prediction | the simulation is the answer |

**If every prediction in a case is `implementation_check`, that case can validate
but can never discover.** Say so here, at S2, rather than discovering it at S8.
One case in this repository is exactly that and is kept as the negative example.

### Do not truncate `value`

Writing `1.00251` gives a design power of `0.5628σ`; full precision
`1.0025063` gives `0.5618σ`. **Use the analytic value as it comes out.**

## 2. `tolerance` — the formats, and the trap

| Format | Meaning | Example |
|---|---|---|
| `±X%` | relative band on the predicted value | `±1.5%` |
| `±X` | absolute band | `±0.03` |
| `>X` / `p>X` / `R^2>X` | one-sided lower bound | `>0.99`, `p>0.05` |
| `<X` | one-sided upper bound | `<1e-12` |

Parser: `simbot.validate.parse_tolerance`. **It raises if it cannot read the
string** — passing silently would let that item through with no verdict at all.

### ❌ A wide tolerance so that any result PASSes

**Forbidden, and subject to review.** The judge now catches it: a `PASS` whose
deviation exceeds `3σ` is marked `PASS ⚑` and raised as a problem.

★ This actually fired. Predicting the MSD plateau as exactly `2d = 4.0` came out
`3.54σ` off *inside* a `±2%` band. The cause was **leaving a known bias out of
the prediction** — `plateau = 2d⟨x*²⟩`, so the EM bias multiplies in.
Full account:
[`wide-tolerance-hides-significant-deviation`](../../../../knowledge/wiki/findings/wide-tolerance-hides-significant-deviation.md)

⇒ **Fold known systematic bias into the prediction.** Record not the ideal value
but the value that *should come out of this scheme*.

## 3. `competing_value` — compute the power in advance

Recording a competing hypothesis lets the judge compute the **design power**:

```
power = |prediction − competing| / SE
```

Below `1σ`, that measurement cannot distinguish the two hypotheses, so it is
`INCONCLUSIVE`.

### ★ Compute this at the prediction stage

If the SE for four seeds is knowable in advance
(`estimators.samples_for_variance_precision`), then **you know beforehand which
predictions will be undecidable.** Two items on the first end-to-end run were,
and "INCONCLUSIVE expected" was written into the prediction document — so when
the result came in, it read as a foreseen limitation rather than a failure.

**Knowing in advance that something is undecidable gives two options:**
1. Proceed and fix `INCONCLUSIVE` as a fact (if the conclusion does not depend
   on it)
2. Change to conditions where power exists — verify an integrator at a
   **deliberately large `dt*`** (`3.8σ` at `dt*=2e-2`, `0.56σ` at `5e-3`)

> **You cannot verify the integrator at the production `dt*`.** Shrinking `dt*`
> shrinks the bias, and shrinks the verifiability with it.

⚠️ And a related trap from the same family: **testing an integrator assumption at
finite temperature can have no power at all.** Comparing overdamped against
`Langevin(kT=0)` across 7 frequencies bounded the difference at 0.159 %, whereas
the thermal comparison **could not have excluded a 47 % effect**. Use a `kT=0`
deterministic difference, where noise is zero and transients cancel as a common
mode.

## 4. `alternatives` — how the prediction could be wrong, and the signal

```yaml
alternatives:
  - 'If dt* is large, <x*^2> comes out high by 1/(1-dt*/2) — this is normal.'
  - 'Kurtosis is 3 - 1.2 dt*, not exactly 3 (uniform noise). Exactly 3.000 is the suspicious case.'
  - 'If tau_trap is off, that signals a wrong eta assumption or a wrong reading of a.
     <x^2> is unaffected, so the two can be separated and identified.'
```

**Writing down in advance what a wrong prediction would look like makes causal
inference at S7 far faster.**

## 5. `regimes` — dimensionless groups and distance from a boundary

```python
from simbot.nondim import groups
groups(spec)     # only what is computable comes out. Do not add what is not
```

Record **how far each dimensionless group sits from its regime boundary.**
`k*_σ = 2.4e5` is 5.4 decades from the boundary (`~1`) — you have to write down
"extremely far" to be entitled to say "a small error in this value does not
change the conclusion."

Near a boundary (`Pe = 45` with MIPS at `Pe_c ≈ 40–60`), **ask the user** — that
is a legitimate use of the question budget.

★ A convention trap here, execution-verified: the numerator of Pe/Wi is
`|E| = √(2E:E)`, **not** `|∇u|` — vorticity cannot change an isotropic
equilibrium. Using `√(E:E)` is off by `√2`, and a pure rotational flow has
`|∇u| ≠ 0` but **Pe = 0**.

## 6. Sealing

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>
```

The CLI writes `SEALED.sha256` **before running**. It is standard `sha256sum`
format, so it verifies without our code:

```bash
shasum -a 256 -c SEALED.sha256
```

**If the prediction is edited after sealing, S7 does not build the comparison
table.** Do not look for a way around this — if one exists, the seal guarantees
nothing.

Running without `--prediction` makes the CLI warn:
> ⚠️ no prediction file — proceeding without the S7 comparison. There is nothing
> to seal, so there is no device preventing post-hoc rationalization.

**That is fine while exploring. It is not fine for a run that produces a
conclusion.**

⚠️ **Know which engine you are on.** Sealing lives in `cli.py`, which has a
runner for one card only. All 8 cases run through `bdbot`, which has **no
sealing** — so for those, S2 is currently a discipline rather than a mechanism.
Say that in the conclusion instead of implying a seal existed.
→ [docs/00 §5](../../../../docs/00-merge-decisions.md#5--known-seams)

## 7. Gates

- ≥1 quantitative prediction, each with `tolerance`, `basis` and a `role`
  (`Prediction.problems()` checks the first two)
- If the drawing has a user-drawn expected curve, state the agreement or
  disagreement with it
- Is every known systematic bias reflected in the predicted value?
- Is at least one prediction a `hypothesis`? If not, say explicitly that this
  case cannot discover anything
