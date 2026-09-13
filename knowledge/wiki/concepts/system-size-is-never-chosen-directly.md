---
type: concept
author: agent
drafted: 2026-09-10
confirmed_by:
cites:
  - bdbot/physical.py
  - tests/test_size_gate.py
  - knowledge/wiki/findings/order-parameter-magnitude-cannot-identify-a-phase.md
  - campaigns/soft2d_fss.py
  - campaigns/soft2d_nconv.py
  - intake/soft-r3-2d-A-sweep/system.yaml
  - intake/trap-drag-2d-hex300/system.yaml
  - intake/trap-2d-5um/system.yaml
id: system-size-is-never-chosen-directly
kind: models
tags: [system-size, N, box, phi, finite-size, replication, provenance]
created: 2026-09-10
updated: 2026-09-10
runs: []
confidence: medium
supersedes: []
---

## Summary

**`particle.count` was carrying three different quantities under one name**, and
in none of the 8 cases was `N` chosen directly — something else was chosen and
`N` followed. Every case now declares `structure.size.role` and, more
importantly, `fixed_by`: *what* it followed from.

| role | `N` is | Raising it | 2026-09-10 |
|---|---|---|---|
| `object` | the size of the **thing** — a 25-bead chain | makes **a different object** | 3 |
| `system` | how much of an **infinite system** is simulated | is the **same material in a bigger box** | 3 |
| `replication` | a **sampling multiplier**; the physical system is smaller | only shrinks the error bar | 2 |

The three carry **opposite obligations**, which is why one field could not serve
all three. The same `converge --only N_double`:

- `replication` — a pass is reassurance ✓
- `object` — a category error. `n=50` is not the converged limit of `n=25`, it
  is a different chain
- `system` — a pass may be a **false negative**: the two sizes may simply have
  been too close to see the scaling
  ([findings/order-parameter-magnitude…](../findings/order-parameter-magnitude-cannot-identify-a-phase.md))

## `phi` links N and L but fixes neither

For a `system`, the packing fraction is one equation in two unknowns:

```
phi = N pi d^2 / (4 L^2)     (2D)
```

So a **second constraint** pins one of them and `phi` gives the other, and *that
second constraint* is what has to be recorded — writing down `phi` alone does
not say what set the size. `soft-r3` went from the sketch's `N=100` to `N=400`
because of `r_c < L/2`, **not** because `phi` moved.

★ **And the direction is not fixed.** Two of the three `system` cases pin `L`
and let `N` follow; `trap-drag` pins `N` and lets `L` follow, because a
hexagonal lattice has to close on itself — `N = n_x·n_y` is an integer pair.
Hence the `pins` field — which was called `follows` and held the *pinned* quantity, i.e. the opposite of what the name read, until adversarial review caught it.

| case | phi | second constraint | pins | follows from phi |
|---|---:|---|---|---|
| `soft-r3-2d-A-sweep` | 0.35 | `r_c < L/2` (minimum image) | L | N: 100 → **400** |
| `network` | 0.10 | measurement radius, `L/2 = 10 d` for stress propagation; cost caps it | L | N = φ(L/d)³·6/π = **1528** |
| `trap-drag-2d-hex300` | 0.35 | commensurability, `N = n_x·n_y = 17×18` | **N** | L_x, L_y |

## `replication` — the case that motivated the split

`trap-2d-5um`'s sketch has **one** particle. `N = 1000` is a device, and the
file says so in prose — *"a simulation choice … to get a statistics multiple …
the sketch has 1"* — while filing 1000 under `particle.count`, beside diameter
and density, where it reads as a property of the physical system. It is not one.

Recording it as `replication` with `n_physical × n_replicas = N` makes the claim
falsifiable and keeps three things honest: `phi` computed here would be
meaningless; `N=1000` (trap) is not comparable to `N=400` (soft-r3); and the box
is a **container for replicas**, so its checks are about leakage, not finite
size.

Measured for `trap-2d-5um`: at N = 1000 the statistical error on `<x^2>` is
**±0.0481 %** (`runs/trap-2d-5um__a5ef4f45d589`).
⚠️ Not "against goal.yaml's 5 % target" — that target is on **`f_c`**, the
`answering_quantity`, a different observable (`f_c` itself came in at +1.17 %).
`<x^2>` carries no stated *percentage* target; what the case does declare for it
is a window-length criterion, `T_obs/τ_k ≥ 100` (soft), plus a `target_bias` of
0.1 % that sets `dt` rather than the error bar. Corrected 2026-09-11.
⚠️ An earlier draft of this page quoted **±0.594 %** beside N = 1000; that is the
**N = 200 smoke run** — wrong by 12×, caught in review. Anchors on a `1 d` lattice, `L = 32 d`,
`l_k = 20.4 nm = 0.00407 d` ⟹ anchor→edge = **123 l_k**, anchor→anchor =
**246 l_k**. Overwhelmingly safe.

⚠️ **Recorded as a question, not a finding.** The check in place is `2 l_k / L`
— the fluctuation against the whole box — while the distance a replica would
actually have to cross is anchor→edge = `0.5 d`. If that reading is right the
check is 64× looser than the constraint it stands for; but the check's own basis
line says it guards the box, not the lattice cell, so the two may be guarding
different things. **Not verified.** Both pass here, by 123×.

⚠️ And `trap-2d-5um`'s `N=1000` was never inverted from a target. Its record is
`confirmed_by: run` — a round number, validated afterwards when the four
observables matched the analytic solution. `abp-rod` then inherited it
(*"1-A와 같은 방식"*). The enum names that state `confirmed_by_run` rather than
dressing it as `target_precision`.

## `object` — why it is not a finite-size effect

A chain of `n` beads has an `n → ∞` limit, so one could call finite `n` a
finite-size artefact. This repository does not, and the reason is a result:
`chain-bend-2d-dlvo` reports `K′` falling monotonically toward 0 across
`n = 5…25` as **evidence for its conclusion**. Relabelling that as an artefact
to extrapolate away would file a finding as a failure — the exact failure mode
rule 7′ exists to prevent.

`chain-bend-2d-oscill`'s `n = 25` shows the third kind of reason, which is
neither statistics nor a thermodynamic limit: **design power**. `delta_max ~ n²`,
so at `n = 11` the maximum deflection is 74 nm ≈ 3 `l_k` and the bending signal
sits inside the trap's thermal noise; `n = 25` opens the window to 429 nm.

## Enforcement

| | |
|---|---|
| field | `structure.size` — `role`, `n`, `fixed_by`, `what_would_change`, plus what the role owes |
| check | [`bdbot.physical.check_size`](../../../bdbot/physical.py), called from `validate()` — **and from `run.execute()`**, because `validate()` is reachable only from `bdbot.cli system check` / `status`, so the gate blocked no run until adversarial review said so. Same shape as its two siblings: `require_structure` defaults to False (the 254 archived specs predate the block), but a block that IS present is **always** validated |
| per role | `system` → `phi`, `pins`, `finite_size`, `finite_size_status` · `replication` → `n_physical`, `n_replicas` (product cross-checked against `particle.count`; a non-numeric operand is reported, not skipped) |
| `fixed_by` | validated **against the role** — `minimum_image` is meaningless for an `object` |
| ⚠ what the run gate does NOT check | it validates the block against the spec's **document** (`system.particle.count`) and prints a **warning** when the N the run integrates (`spec.params`) differs — not an error, because `bdbot.smoke.PROFILES` shrinks N deliberately (`trap-2d-5um N=200`, `verified: True`), so 3 of the 278 archived specs disagree by design and refusing them would refuse the repo's own smoke contract. **Two checks do not run there at all**: the phi closure and the `stated_in_source` cross-check need the Provenanced leaves and the L0 file, which the spec-side path does not carry. Both now emit a **warning** saying so — they used to be an invisible `info` and a silent `return []`, so the run path printed "structure OK -- dim/size validated" over a 350 %-packing phi. `bdbot.cli system check` is the full gate |
| tests | [`tests/test_size_gate.py`](../../../tests/test_size_gate.py) — **92** as of 2026-09-13. Unwiring `check_size` fails 44 · `check_dim` 18 · the `execute()` gate 1. That last one is deliberately a single test, and it is the non-vacuous one: the source-grep beside it survives the gate being turned off (measured), so what catches a disablement is calling `execute()` with a broken block and requiring the refusal *before* `build_fn` |
| override | a stated N may be overridden — `soft-r3` had to — but not silently: `stated_in_source` ≠ `n` makes `approved_by` required |

## Scope / limits

**The `finite_size` obligation is met in 0 of 3** — one partial, two not started.
⚠️ Separately, and often confused with it: the **phi closure** runs and agrees on
**2 of the 3** `system` cases (`soft-r3`, `trap-drag` — the latter with the
product of its rectangular box's axes). Only `network` cannot be closed, because
it is 3D and compressed, so no box length is recorded; it says so as a warning.
An earlier note here said all three could not be closed.
The gate branches on a status word (`done`/`partial`/`not_done`), not on a
substring of the prose, which failed in both directions. It warns rather than blocks —
whether an N-sweep was run is a property of the campaign, not of this file — but
it is said on every read, because an unrecorded gap and a recorded one must not
look alike.

- `soft-r3` — **partial**. The ladders exist (`campaigns/soft2d_fss.py`,
  `soft2d_nconv.py`) but **A = 100, the production point this `n = 400` was
  chosen for, has no N-sweep**, and the file's own convergence item CV3 never
  ran. This page first said `done`, reading the existence of a ladder as
  coverage of the point — corrected in review.
  ⚠️ A second correction in the same block: the `r_c` 5a→7a convergence check was
  run **at** this production point, `Γ = 29.7`, which the run measures as a
  hexagonal **crystal** (`state_predicted = 결정`, ψ₆ = 0.885) — an earlier draft
  called it "deep in the liquid" by comparing this file's `Γ ≡ U(a_mean)/kT`
  against Zahn's `Γ = π^1.5 βU(d)`. In this convention the hexatic window is
  **10.03–10.75**, so `Γ = 29.7` is 2.8× *above* it
- `network` — not done; stage 1 is incomplete (1 run)
- `trap-drag` — not done, and this is the sharp one: **81 run directories, all
  at `N = 306`**, reporting defect count, `psi_6` and `F_drag`, which are
  collective. This case has already had **two** conclusions reversed — defect
  count read as v-independent, then found to be a non-monotonic hump peaking at
  γv ≈ 24; and a 67 % recovery fraction that was one realization, whose ensemble
  mean peaks at 36 % — and both were reversed by adding **seeds, not sizes**

**Not settled here.** `boundary` has no field: all 8 cases are periodic, HOOMD's
box always is, and "open" is already expressed as a large box plus
`declare_absent("box", reason)` in the 3 chain cases. `md.external.wall` (LJ/Gaussian/Morse/…) applied to a
`hoomd.wall` geometry (`Plane`/`Sphere`/`Cylinder`) exists and would be *physics* — a substrate or a
channel — not a numerical convenience, so it enters as an interaction (rule 9)
if an experiment has one.

See also: [dimensionality-has-no-default](dimensionality-has-no-default.md) ·
[order-parameter-magnitude-cannot-identify-a-phase](../findings/order-parameter-magnitude-cannot-identify-a-phase.md)
· skill `bd-physics` §4
