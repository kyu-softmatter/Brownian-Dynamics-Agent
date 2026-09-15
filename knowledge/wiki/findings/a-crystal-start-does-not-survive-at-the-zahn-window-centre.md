---
type: finding
author: agent
drafted: 2026-09-15
confirmed_by:
system: soft-repulsive-2d
dynamics: equilibrium-structure
id: a-crystal-start-does-not-survive-at-the-zahn-window-centre
kind: findings
tags: [soft-r3, 2d-melting, hysteresis, initial-condition, pre-registration, equivalence]
created: 2026-09-15
updated: 2026-09-15
confidence: high
supersedes: []
cause_class: interpretation
stage: S7
question: "Was the archived isotropic-liquid reading at Gamma_Zahn = 57.9 a supercooling upper bound, i.e. does a crystal survive there?"
answer: "No. All 5 crystal seeds melt within about 25 tau_B and reach psi6_global ~ 0.13, equivalent to a random start inside +-0.10 (TOST, 0.028). The sealed prediction was HYSTERESIS and it was wrong."
reproduced: yes
runs:
  - soft-r3-2d-A-sweep__A34.938-eq0-rc7.8-hex20x20-s20260914__5e60542041
  - soft-r3-2d-A-sweep__A34.938-eq0-rc7.8-box20x20-s20260914__fd904e5c12
  - soft-r3-2d-A-sweep__A106.224-eq0-rc7.8-hex20x20-s20260914__a688d9aa11
  - soft-r3-2d-A-sweep__A10-eq0-rc7.8-hex20x20-s20260914__1c35e8f260
cites:
  - campaigns/s30_preregistration/prediction.yaml
  - campaigns/s30_preregistration/analysis_plan.yaml
  - knowledge/wiki/findings/order-parameter-magnitude-cannot-identify-a-phase.md
  - knowledge/wiki/findings/box-shape-confounds-initial-condition-comparison.md
  - knowledge/wiki/findings/coarse-sampling-hides-the-whole-transient.md
  - knowledge/source/papers/1999-zahn-two-stage-melting-2d.md
---

# A crystal start does not survive at the centre of Zahn's window — the supercooling argument does not bite at `Γ_Z = 57.9`

> **The sealed prediction was `(n_H, n_R) = (5, 0)` — HYSTERESIS. The measurement
> is `(0, 0)`.** All five crystal seeds melted. The two arms are *equivalent*
> inside ±0.10 in `ψ₆^global` by a two-one-sided-tests result of **0.028**, which
> is a positive statement and not a failure to reject.

## What was asked

`runs_s1s8/2026-07-29_soft-r3-hexwin` measured `Γ_Zahn = 57.9` — the centre of
Zahn's hexatic window `55.87–59.88` — from **random** starts at
`N = 256/576/1024` and read isotropic liquid at every point. Correction ③ of
[[order-parameter-magnitude-cannot-identify-a-phase]] argued that reading is an
**upper bound** rather than the boundary: near a first-order transition a random
start can sit behind a nucleation barrier arbitrarily long, and a supercooled
liquid is also steady, so a random start cannot locate the boundary. It
concluded that bracketing it requires **starting from a crystal and seeing
whether it melts**, filed as S30 and never run — no run in this repository had
ever started from anything but a random configuration.

That argument is this repository's own, and it had never been measured.

## The design, and the two findings it obeys rather than rediscovers

| | |
|---|---|
| primary statistic | `(n_survive_H, n_survive_R)` — the count of seeds in **each** arm with `ψ₆^global` over the last quarter above `PSI6_CRYSTAL_FLOOR = 0.5`. A pair of integers, no variance estimate |
| arms | H = perfect 20×20 hexagonal lattice · R = RSA, **in the same commensurate box, with the same seeds** |
| controls | C1 at `Γ_Z = 176` must stay ordered · C2 at `Γ_Z = 16.6` must melt |
| window | `100 τ_B`, **sampled from `t = 0`** |
| seeds | 5 per arm |

★ **The box is the same in both arms because that is already a finding.**
[[box-shape-confounds-initial-condition-comparison]] measured an `S(k)` six-fold
modulation splitting **50×** and an energy splitting **12.8σ** when a hexagonal
and a random start were compared while the box shape changed with them — at a
coupling where both are deep liquids. So `init` and the box are two knobs here:
`cases/soft_r3_2d.py` takes `--init {rsa,hex}` and `--box {square,hex}` and
**refuses `--init hex --box square`**, while `--box hex --init rsa` is the
matched random arm and is required. `simbot.run` had reached the same conclusion
independently and carries the same refusal.

## What was measured

```
seed        H late    > floor      R late    > floor
20260914    0.1319      no         0.1652      no
20260915    0.1106      no         0.1386      no
20260916    0.1296      no         0.1303      no
20260917    0.1249      no         0.1265      no
20260918    0.1363      no         0.1341      no

(n_survive_H, n_survive_R) = (0, 0)

Delta = -0.0123   sigma = 0.0081   nu = 6.82   |t| = 1.51   p = 0.18
TOST  |Delta| + t(nu,0.95)*sigma = 0.0278  <  delta = 0.10   -> FIRES
```

| id | sealed | measured | verdict |
|---|---|---|---|
| **P1** `hypothesis` | `(5, 0)` → HYSTERESIS | `(0, 0)` | **NO HYSTERESIS DETECTED** |
| G1 gate | `ψ₆ > 0.5`, defect `< 0.10` | `0.8961`, `0.0000` | PASS |
| G2 gate | `ψ₆ < 0.30`, defect `> 0.20` | `0.0531`, `0.4637` | PASS |
| P2 `implementation_check` | defect(R) `= 0.285 ± 5 %` | `0.2899` | PASS |
| P3 `implementation_check` | energy identity `≤ 2 %` | `0.0213 %` | PASS |

⚠ **`Δ` is negative.** The crystal arm ends *slightly more disordered* than the
random arm. It is not significant (`p = 0.18`) and it is inside the equivalence
band, so the correct statement is "indistinguishable" — recorded here so that
the sign is not quietly dropped.

**Both controls passing is what makes P1 mean anything.** C1 keeps `ψ₆ = 0.896`
with **zero** Voronoi defects, so the protocol does not melt whatever it is
handed; C2 melts a crystal to `0.053` with defect fraction `0.464`, so it is
capable of melting one. Without the pair, `(0, 0)` would be consistent with a
broken protocol.

## What this says, and what it does not

**Says.** The answer at `Γ_Z = 57.87` does not depend on where you start, on
this timescale, at this resolution. So the archived isotropic-liquid reading is
**robust to the initial condition**, and correction ③'s supercooling argument —
which is sound in general — **does not bite at this coupling**. The existing
conclusion is strengthened, not overturned.

**Does not say.** Nothing about hexatic: Zahn's crystal→hexatic boundary is
defined by the *dynamic* Lindemann parameter `γ_M` crossing `0.033`, which this
repository does not implement, and `ψ₆` magnitude plus defect fraction cannot
separate hexatic from crystal. Zahn is why this coupling was chosen and nothing
more. Nothing about where the equilibrium boundary is: both arms give one-sided
bounds and neither is converged. Nothing about `η₆`: one system size.

## ★ The whole separation is in the first quarter

```
                Q1       Q2       Q3       Q4       (25 tau_B each)
H arm mean     0.259    0.138    0.153    0.127
R arm mean     0.132    0.138    0.137    0.139
difference    +0.127   ~0       +0.016   -0.012
```

`ψ₆^global` is **exactly 1.000000** for the perfect lattice
(`tests/test_lattice.py`), `0.65` at the first sample (`0.25 τ_B`), below the
crystal floor by `t ≈ 5–10 τ_B`, and merged with the random arm by `t ≈ 25 τ_B`.

⚠ **The case default `eq_frac = 0.2` would have discarded almost all of it.**
`bdbot.run` builds the equilibration phase with `collect=False`, so at the
default the first **20.00 `τ_B`** of a `100 τ_B` window carries no sample — and
that is where the melt happens. The transient would have been invisible and the
crystal arm would have looked liquid from its first recorded frame. `--eq-frac 0`
was added for this campaign.

⚠ **And the onset is still unresolved.** `1.000 → 0.65` happens inside one
sample interval. That is [[coarse-sampling-hides-the-whole-transient]]
recurring: that finding measured `stride = 0.2 τ_d` against a relaxation time of
`0.031–0.098 τ_d` and concluded the system looked equilibrated from the first
frame. Here the stride is `0.25 τ_B` and the melt onset is faster than it. **The
melting *rate* is not measured by this campaign** — only that melting completes
well inside the window.

## Prevention

- `--eq-frac 0` on `cases/soft_r3_2d.py`; `n_eq = 0` drops the phase rather than
  running a zero-step one
- `tests/test_soft_r3_init.py::test_a_run_with_no_equilibration_phase_writes_its_artefacts`
  — because the first attempt at this campaign **crashed after printing PASS**:
  `make_plots` read `res["eq_trace"]`, which does not exist when `n_eq = 0`, so
  `result.txt` was never written, the exit code was 1, and the campaign driver
  read a complete and correct run as a failure and stopped. It cost the campaign
  its first run
- `bdbot/lattice.py` + `tests/test_lattice.py` — the lattice is verified to give
  `ψ₆ = 1.000000` and zero defects before it is used as an initial condition
- to measure the melting *rate*, a run with a stride below `0.05 τ_B` over the
  first `5 τ_B` is needed. Not run

## Provenance of the claim

The prediction, the decision rule, the equivalence margin and the analysis plan
were written and **cryptographically sealed before any run existed**
(`runs/*/SEALED.sha256` over `prediction.yaml` + `analysis_plan.yaml`, verified
by `bdbot.runcard.verify_seal` at the top of the analysis and by CI thereafter),
and every run was gated at execution by `require_seal=True` and
`require_approval=True`. Three revisions of the prediction were refused by six
adversarial reviews before sealing; two of the refusals were self-contradictions
between the two documents, and one — `analysis_plan.yaml` not being valid YAML —
had stood since the first draft because nothing parsed it.

⚠ One deviation to record. The drift gate was applied **per run**, where the
sealed text says "for each arm separately (H, R, C1)" and an arm is five runs.
That is an ambiguity in the sealed document. `H1` fires under the per-run
reading. Computed both ways: `Δ = −0.0123`, TOST `0.0278` with all five seeds;
`Δ = −0.0136`, TOST `0.0302` with `H1` excluded. **The verdict is unchanged.**

## See also

[[order-parameter-magnitude-cannot-identify-a-phase]] ·
[[box-shape-confounds-initial-condition-comparison]] ·
[[coarse-sampling-hides-the-whole-transient]] ·
[[low-seed-pilots-give-optimistic-design-power]] ·
[[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]] ·
`knowledge/wiki/concepts/system-size-is-never-chosen-directly.md`
