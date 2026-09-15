---
type: concept
author: agent
drafted: 2026-09-10
confirmed_by:
cites:
  - docs/history/2026-07_bd_agent_00_decision_log.ko.md
  - bdbot/physical.py
  - verify/verify_dim_gate.py
  - verify/verify_abp_rod_dimension_impact.py
  - tests/test_dim_gate.py
  - knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md
  - knowledge/source/papers/2022-modica-porous-media-active-diffusion.md
  - knowledge/source/papers/2020-takatori-motility-induced-buckling.md
id: dimensionality-has-no-default
kind: models
tags: [dimensionality, 2d, 3d, scope, decision, provenance, D9]
created: 2026-09-10
updated: 2026-09-10
runs: []
confidence: medium
supersedes: []
closes: [D9]
---

## Summary

**There is no default dimension.** Every case declares `dim.basis` — one of
`given` · `required` · `inherited` · `sufficient` — and `given` takes precedence
over the other three.

This closes **`D9 · 차원 (2D / 3D)`**
([decision log](../../../docs/history/2026-07_bd_agent_00_decision_log.ko.md)),
which had stood at `상태: OPEN` with `기본값: 3D 기본, 2D는 스펙 플래그로 지원`
and an empty `결정일`. The history file is a verbatim provenance record and is
**not edited**; this page is where the decision now lives, the same way
[`no-hydrodynamics`](no-hydrodynamics.md) carries `D11`.

**What forced the change:** the stated default was 3D and **7 of the 8 cases are
2D**, so every one of them departs from it — and, measured 2026-09-10, **none of
the 8 records why.** A default that is contradicted by 100 % of its instances is
not a default; it is an unfilled field wearing one.

## The four bases

| `basis` | It means | Required alongside | Expires? |
|---|---|---|---|
| `given` | The input — sketch, paper, or a human asked directly — **states** it. It is an input, not a modelling choice, so it is not re-derived | `source_kind` (`artefact` \| `human`), plus `source_file` **which must exist** or `confirmed_by` — see below | no |
| `required` | The question does not exist in another dimension. 3D is a **different case**, not a better one | — | no |
| `inherited` | Matched to an existing case so the two are comparable | `compared_with` naming the case(s) | **yes** — if the comparison target changes or is dropped |
| `sufficient` | The degrees of freedom separate and the observable is one of the separated components. 3D would give the same answer; 2D is only cheaper | `checked_observables` | **yes** — when an observable is added |

**Precedence: `given` beats the rest.** If the source states it, the agent does
not re-derive it (rule 3: never invent; rule 5: transcribe before you interpret).

### `sufficient` is narrower than it looks

Two conditions, both required:

1. the Hamiltonian separates per axis — no pair interactions, no bonds, and a
   separable external field
2. the observable is one of those separated components — **not** a radial
   quantity and **not** a counting quantity

"The observable's formula does not contain `dim`" is **not enough**. `chain-bend`'s
`K′` contains no `dim`, yet in 3D there is a second transverse direction sharing
the same bonds, so the chain is a different mechanical object. Under conditions
1–2, `sufficient` holds for exactly **1 of the 8 cases**.

Observables that **do** carry the dimension: `⟨r²⟩ = dim·kT/k` · `MSD = 2·dim·D·t`
· `φ` (`Nπσ²/4L²` vs `Nπσ³/6L³`) · coordination number, loop count, percolation,
`d_f` · `D_r` and HOOMD's `rotational_diffusion` (**off by 2× in 3D** — skill
`bd-hoomd` trap 14).

### ★ `given` has to say WHICH input, and the two kinds owe different things

`basis: given` asserts that the input states the dimension. Until 2026-09-15 the
gate discharged that with the mere **presence** of a `source` field, so a file
and a half-remembered conversation satisfied it identically — the same
"presence is not validity" shape `params.blockers()` had for `derived`.

It was found by auditing which `structure` claims an artefact could settle.
Of the three `given` cases:

| case | what `source` points at | settled by |
|---|---|---|
| `network` | `observation.yaml` A4, quoted verbatim, with `confirmed_by: user` | reading the file |
| `trap-2d-5um` | the sketch's own `U(r) = ½k_t r²,  r = √(x²+y²)` | looking at the image |
| `chain-bend-2d-oscill` | **an exchange dated 2026-09-10 that appears nowhere in that case's intake** — `observation.yaml` never resolves dimensionality at all | only by asking |

The third was put to the user on 2026-09-15 and **confirmed**, so the `given`
stands. What did not stand was the check. Its own `source` prose reads *"The
sketch itself is silent"* — the record openly said no artefact stated it, and
nothing read that.

So the author now **declares** the kind and each kind owes a different thing:

| `source_kind` | means | owes | checked how |
|---|---|---|---|
| `artefact` | a file in the record states it | `source_file` | **the file must exist**, resolved against the case directory then the repo root |
| `human` | a person stated it, in conversation | `confirmed_by` | naming who, and when |

⚠ **Declared, not inferred.** The first attempt scanned the `source` string for
`sketch` / `observation.yaml` / `.jpeg` and labelled exactly the one
conversational case `artefact`, because its text contains the word *sketch*
while saying the sketch says nothing. Parse, do not grep — the same lesson
`A4`'s grep taught, one layer up.

⚠ And the existence check **cannot run on the spec path**. `run.execute()`
rebuilds a `PhysicalSystem` around the spec's `system` document, whose `path`
does not point at the case directory. There it emits a **warning** saying so
(`recorded, NOT verified`) rather than an error that would refuse every run or a
silence that would let the run path print "structure OK" over an unverifiable
claim. `bdbot.cli system check` is the full gate. Same shape as `_phi_closure`,
for the same reason.

## Anti-patterns

**★ A flat drawing is not evidence of 2D.** A hand sketch is *always* drawn on
flat paper, so its planarity carries essentially no information — and this
project has already been wrong that way. `network`'s intake leaned 2D on the
grounds *"(1) all 6 project cases are 2D (2) the drawing is planar"*, honestly
flagged the choice as convenience, and the record now reads
**`resolution: 3D (user-specified 2026-08-06, "physical"). ★ my lean (2D) was wrong`**.

**Precedent without a named comparison is inertia, not a basis.** Ground (1)
above fails for the same reason: nothing was being compared against those six
cases. The discriminator is whether `compared_with` can be filled in.
`chain-relax-2d-dlvo ← chain-bend-2d-dlvo` (its rule-8 static control) is a
legitimate `inherited`; *"we have always done 2D"* is not.

**A case name is not a decision.** Seven of the eight case labels contain `2d`
(`abp-rod-2d-run-flip`, `trap-2d-5um`, …). The one that does not is `network` —
and it is **the only case in which the dimension was ever argued**. n = 8, so
this is a correlation and not a demonstrated cause, but the mechanism is easy to
state: once the name asserts a dimension, the field looks already filled.

## Evidence

**The literature layer is thinner than "ABP is usually 2D" suggests.** Of the 7
ABP-related distillations in [`source/papers/`](../../source/papers/):

| Distillation | Dimension | How it is known |
|---|---|---|
| `2024-quah-graybox-abp-mpc-repo` | 2D | **the actual code** — `box = [width, plate_gap, 0,0,0,0]  # Lz=0 → 2D` |
| `2022-modica-porous-media-active-diffusion` | 2D | stated **with a physical reason** — silica beads sediment to the chamber floor |
| `2020-takatori-motility-induced-buckling` | 2D bidisperse | `engine:` field; "modelling bacterial monolayers as 2D fluid films" |
| `2025-quah-continuum-closures-active-control` | 2D | inherited via its `code_repo` pointer to the graybox repo |
| `2024-quah-mpc-noninteracting-abp` | **not recorded** | the "N, box size, dimension" row is blank |
| `2024-quah-nn-augmented-mpc` | **not recorded** | same row, blank |
| `2024-cheon-motility-partitioning-atps` | **not recorded** | — |

So: 3 state it directly, 1 inherits it, **3 do not record it at all.** The claim
is supported, but by four papers rather than by a field-wide convention — and
⚠️ 38 of this repository's 42 distillations are the group's own published work,
so layer ③ is narrower here than it looks.

**`2022-modica` noticed this exact contradiction and wrote it down:**
*"2D다. 우리 `D9` 기본값은 3D — 이 벤치마크를 쓰려면 2D 경로가 필요"* — the
benchmark is 2D while our own default was 3D.

**Takatori is the precedent for inferring out-of-plane physics from a 2D run.**
The paper's subject is a *three-dimensional* transition — buckling of a bacterial
monolayer out of its plane — and the simulation is 2D. Relevant wherever a 2D
choice is questioned on the grounds that the real system can leave the plane
(`chain-bend`, where 6 of 9 beads are held by DLVO bonds alone and nothing
confines their `z`).

**The choice propagates into the physics before anyone checks it.**
`abp-rod-2d-run-flip` already carries `gamma_bar_2d` — *"2D in-plane harmonic
mean `2/(1/ζ∥+1/ζ⊥)`"* — and `D_bar = kT/gamma_bar_2d`, plus hard constraint H1
*"MSAD ✓ exact (depends only on `γ_r,z`; **2D has one rotation axis**)"*. Moving
to 3D changes how the Perrin friction is averaged and takes the rotation axes
from one to three.

**`network` shows what the impact analysis should look like**, and it was done
*before* deciding: units change (2D `G*` is `[N/m]`, 3D is `[Pa]` — which decides
whether literature comparison is possible at all) and cost scales as
`N ∝ (L/d)^dim`.

## Assignment for the 8 cases

| Case | `dim` | `basis` | Held by |
|---|:---:|---|---|
| `network` | 3 | `given` | user-specified 2026-08-06, "physical" |
| `trap-2d-5um` | 2 | `given` | the sketch's `U(r)=½k_t r²`, `r=√(x²+y²)` fixes it. (`sufficient` also holds — no interactions, separable trap, per-component observables — and is recorded as corroboration) |
| `chain-bend-2d-oscill` | 2 | `given` | user-specified 2026-09-10: optical tweezers. ⚠️ residual: only 3 of ~9 beads are trapped; the rest are held by DLVO bonds alone, so out-of-plane motion of the free beads is not directly suppressed |
| `soft-r3-2d-A-sweep` | 2 | `required` | hexatic / KTHNY does not exist in 3D |
| `trap-drag-2d-hex300` | 2 | `required` | hexagonal lattice, `ψ₆`, commensurability |
| `chain-bend-2d-dlvo` | 2 | `inherited` | `compared_with: [chain-bend-2d-oscill]` — same three-point-bending geometry |
| `chain-relax-2d-dlvo` | 2 | `inherited` | `compared_with: [chain-bend-2d-dlvo]` — its rule-8 static control; a 3D run would not be a control |
| `abp-rod-2d-run-flip` | 2 | `inherited` | `compared_with: [2024-quah-graybox-abp-mpc-repo]` — the lab ABP standard, `Lz=0` in the published code. ⚠️ **partial**: the target is spheres and this is a prolate ellipsoid (aspect ratio 4), so this basis **expires** if that comparison is dropped. Decided 2026-09-10 |

All eight are assigned. The impact numbers behind `abp-rod`'s entry are in
[`verify/verify_abp_rod_dimension_impact.py`](../../../verify/verify_abp_rod_dimension_impact.py)
(4/4), which recomputes the case's own `gamma_bar_2d` from `ζ∥`/`ζ⊥` first and
refuses to report anything if that does not reproduce:

| | 2D (as built) | 3D | ratio |
|---|---:|---:|---:|
| `γ̄` [kg/s] | 7.2122e-9 | 7.5271e-9 | 1.044 |
| `D_t` [µm²/s] | 0.5743 | 0.5503 | 0.958 |
| **`D_eff`** [µm²/s] | **4.866** | **2.729** | **0.56×** |
| MSD prefactor `2·d·D_eff` | 19.46 | 16.37 | 0.84× |

**Rotation dominates, not friction.** Perrin moves `γ̄` by 4.4 %; the `(d−1)`
director factor (trap 14) halves `τ_eff` from 0.343 s to 0.261 s, and the two
together move `D_eff` by a factor 1.8. Under the *weakest* alternative
assumption (`Λ_r` unchanged) it is still 0.70×, so **`D_eff` moves by 30–44 %
whichever way the one unverified input goes** — the `τ_eff` composition is
marked `not_verified` in the case's own model_notes (rule 7). The MSD prefactor
looks mild at 0.84× only because the `2d` factor grows while `D_eff` shrinks;
the two effects are cancelling, not absent. This case measures MSD and MSAD, so
the dimension is load-bearing for its stated goal.

## Enforcement

This decision has teeth, which is the difference between it and the three
practices CLAUDE.md rule 10 records as written-and-never-enforced.

| | |
|---|---|
| the field | `structure.dim` in `system.yaml` — `value`, `basis`, `alternatives`, `what_would_change`, plus whatever the basis owes |
| the check | [`bdbot/physical.check_dim`](../../../bdbot/physical.py), called from `validate()` — **not** from a sibling tool, because `health.gate()` was reachable only from `tools/health.py` and no run ever gated itself |
| what fails | the section absent · `basis` not one of the four · a required field missing or blank · `structure.dim.value` disagreeing with the hashed `dimensions` · `given` without a `source_kind`, or with one whose obligation is unmet, or naming a `source_file` that does not exist |
| the adversarial test | [`verify/verify_dim_gate.py`](../../../verify/verify_dim_gate.py) — 28/28, breaking one field at a time, and requiring the eight **real** cases to pass. [`tests/test_dim_gate.py`](../../../tests/test_dim_gate.py) adds the `source_kind` obligations: 6 break cases, 3 mutations of the branch, all caught |
| run_id | unaffected. 278/278 archived specs keep their hash, `structure` is in `runid.DOC_KEYS`, and `dimensions` 2 → 3 still re-ids |

Both expiring bases announce themselves on every read: `inherited` prints its
`compared_with` target, `sufficient` prints how many observables it was checked
against. That warning is the only thing that makes an expiry visible, since an
expired basis looks exactly like a valid one.

## Scope / limits

**A `given` from a human overrides any of the other three, at any time.** That
is what precedence means, and `network` is the worked example: the intake leaned
2D, the user said 3D, and the user won. So `abp-rod`'s `inherited` is not a
lock — an instruction to make it 3D supersedes it. What the basis buys is not
authority over the human; it is that **the cost of the override is already
written down** (the table above, and the list below), so the instruction is
given with the consequences visible rather than discovered afterwards.

**What this page does not settle.**

- **`boundary`.** This page settles the dimension; system size is settled in
  [system-size-is-never-chosen-directly](system-size-is-never-chosen-directly.md).
  The boundary has no field and does not need one: all 8 cases are periodic,
  HOOMD's box always is, and "open" is already expressed as a large box plus
  `declare_absent("box", reason)` in the 3 chain cases. A wall would be
  *physics*, not a numerical convenience, and enters as an interaction (rule 9).
- **The other 17 `OPEN` decisions.** The log holds 38 decisions, **18 of them
  `OPEN`**, and there is no live register saying which have since been settled
  elsewhere. This page closes one and inherits that discoverability problem.

**When a basis expires**, the case is not wrong — it is *unjustified*, which is a
different state and has to be re-derived rather than assumed:

- `inherited` — the comparison target moved, or the comparison was dropped
- `sufficient` — an observable was added. Re-check conditions 1–2 against the
  **new** list, not the one in `checked_observables`

## S7 diagnostic signals

Put the dimension on the candidate-cause list when:

1. a rotational quantity (`D_r`, MSAD, `τ_r`) is off by a factor near **2** —
   `rotational_diffusion` is the director decay rate itself and is off by 2× in
   3D (trap 14)
2. a modulus cannot be compared against the literature because the units do not
   match — 2D gives `[N/m]`, 3D gives `[Pa]`
3. a friction or diffusion coefficient disagrees at the ~10 % level in a case
   using Perrin factors — the in-plane harmonic mean and the 3D average differ
4. a structural or connectivity observable (coordination, loops, percolation,
   `d_f`) is compared against a source whose dimension was never recorded

See also: [no-hydrodynamics](no-hydrodynamics.md) · skill `bd-physics` §4 ·
skill `bd-hoomd` trap 14 · [D9 in the decision log](../../../docs/history/2026-07_bd_agent_00_decision_log.ko.md)
