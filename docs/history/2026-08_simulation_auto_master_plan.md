# Brownian Dynamics Simulation Bot — Master Plan

> Written 2026-08-03 · last updated **2026-08-04** · **version v0.5**
> Target engine: [HOOMD-blue v7.1.1](https://hoomd-blue.readthedocs.io/en/v7.1.1/)
> Execution environment: macOS (Apple Silicon, 10 cores / 16 GB), CPU only, local
> v0.1 → v0.2: multimodal intake · a literature knowledge base · an explicit non-dimensionalization layer ·
> raw-data storage and re-analysis · a human approval gate · a post-mortem learning loop
> v0.2 → v0.3: dimensions-first promoted to a hard invariant · Claude Code as the agent runtime ·
> the need for a scope extension discovered
> v0.3 → v0.4: **opening the scope + the physics module registry** (§5.6) · **Phase 0 complete** (capabilities measured,
> 9 traps) · **Phase 8-min complete** (CLAUDE.md + 2 skills) · **Phase 1-A complete** (4 observables
> matching their analytic solutions) · **the isolation checks reclassified into model/integration/geometry/statistics**
> (§6.4 — the `τ_p/dt` error corrected) · a quantitative rule that **inverts dt from the bias** ·
> **the KB starts from `record.json`** (§7.0)
> v0.4 → v0.5: **principle 8** (verification ≠ agreement) · **principle 9** (verifying an independent element alone) +
> **9.1 no existing theory in a combination** + **9.2 pre-registering the prediction** · the principles reordered 1→9 ·
> `bdbot/` 16 modules + the CLI + L3 `NondimSpec` complete · traps 13 and 14 (ABP) · KB 37 entries

---

## 0. The one-line definition

**An LLM agent that interprets a physical system from a hand-drawn sketch, a note or the literature, proposes and non-dimensionalizes parameters on the strength of accumulated knowledge, runs a Brownian dynamics simulation, and feeds even the successes and failures of that back as knowledge.**

The crux is **the closed loop**. The bot has to get smarter the more it is used.

```
       ┌─────────────────── the knowledge base (KB) ───────────────────┐
       │  literature distillations · dimensionless-number coordinates · our runs' successes and failures    │
       └──▲──────────────────────────────────────────────┬────┘
          │ (reference: the grounds for a parameter proposal)      (feedback: the post-mortem) │
          │                                              ▼
  a sketch/note/picture → reading the physical system → [human approval] → non-dimensionalization → the simulation → analysis
                     (SI, dimensional)              (based on the scale ledger)
```

**Dimensions come first.** Every system is fixed as an SI-unit physical system first, the main length,
time and energy scales present in the system are enumerated in a ledger, and the reference is chosen from
among them to non-dimensionalize.
There is no bypass route that starts from dimensionless values ([principle 3](#principle-3--dimensions-come-first-non-dimensionalization-comes-after--an-invariant)).

---

## 0.1 The current status (2026-08-04) — this section is the single source of truth for the status

| Stage | Status | The substance |
|---|---|---|
| **0** a general environment + the HOOMD capabilities measured | ✅ | `environment.yml` · [`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md) · 15 APIs working |
| **8-min** knowledge capture | ✅ | [`CLAUDE.md`](../../CLAUDE.md) 9 absolute rules (2 of them 7 — the `7`/`7'` convention) · **3** skills (`bd-hoomd` `bd-physics` `bd-intake`) |
| **1-A** `trap-2d-5um` end to end | ✅ | 4 observables matching their analytic solutions · 2 runs |
| **1-B** `soft-r3-2d-A-sweep` end to end | ✅ | 8 runs · 5 verifications + 2 convergence checks |
| **1-C** abstraction → the `bdbot/` package | ✅ | **16 modules** + the CLI (`status/intake/interactions/system/nondim/run`) |
| **the front end** the intake schema + checks | ✅ | adversarial checks 27/27 |
| **L3** `NondimSpec` (the only L2↔L4 contract) | ✅ | `specs/` 3 entries · adversarial checks 33/33 |
| **KB** `record.json` + `kb/entries/` | ✅ running | **37 entries** (8 with a run · 29 knowledge without one) |
| **L4** the numerical-health verdict | ✅ | [`bdbot/health.py`](../../bdbot/health.py) · [`tools/health.py`](../../tools/health.py) · `cli health` · adversarial 31/31 · runs 33/33 HEALTHY |
| 1-D two case scripts (`chain-bend` · `trap-drag`) | ⬅️ **next** | — |
| 1-D three cases · 2 · 8-rest · 9 | waiting/blocked | §0.2 |

**The cases, measured** (`$PY -m bdbot.cli status`)

| Case | L0 | L2 | Script | Runs |
|---|---|---|---|---|
| `trap-2d-5um` | READY | READY | ✅ | 2 |
| `soft-r3-2d-A-sweep` | READY | READY | ✅ | 8 |
| `abp-rod-2d-run-flip` | READY | READY | ✅ | 3 |
| `chain-bend-2d-oscill` | READY | READY | ❌ | 0 |
| `trap-drag-2d-hex300` | READY | READY | ❌ | 0 |

**All 7 verification scripts PASS**: `verify_{1c_equivalence, bdbot, intake_guards,
l3_spec_gaps, nondim_guards, pair_table, skill_snippets, health}`

**What does not exist yet** (in the plan but unimplemented — deliberately):
the `PhysicsModule` registry (§5.6 is design only) · the L4 execution layer · sweep orchestration ·
hooks/subagents/slash commands (Phase 8-rest) · the SQLite KB (§7.1) · raw-data Tier C

---

## 0.2 What remains — measured by "does a sketch carry over well into a simulation" (2026-08-04)

**3 of the 5 sketches** have passed the whole route. The L2 blocks of §0.2-A are **all resolved**
(`bdbot.cli status`: all 5 are L0 and L2 READY).

| Sketch | Intake | The physical system | Non-dimensionalization | Execution | Runs |
|---|---|---|---|---|---|
| `trap-2d-5um` | ✅ | ✅ | ✅ | ✅ | 2 |
| `soft-r3-2d-A-sweep` | ✅ | ✅ | ✅ | ✅ | 8 |
| `abp-rod-2d-run-flip` | ✅ | ✅ | ✅ | ✅ | 3 |
| `chain-bend-2d-oscill` | ✅ | ✅ | ⏸ | ⏸ | 0 — **no script** |
| `trap-drag-2d-hex300` | ✅ | ✅ | ⏸ | ⏸ | 0 — **no script** |

### A. ~~What is blocked right now~~ → **resolved** (2026-08-04)

The three cases' physical-system gaps (`U_ij` · the pair potential · `R`) were filled in by human confirmation
and all became L2 READY. `abp-rod` went as far as execution (3 runs).

**One judgment does remain** — `abp-rod`'s **anisotropic translational friction was concluded to be impossible in BD**
(no HI. `γ⊥/γ∥ = 1.000000` measured). The isotropic average `γ̄` is in use, and that limit
remains in the short-time MSD (§20 A · [`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md) §5).

### B. What remains in the pipeline itself

> **The direction (user instruction 2026-08-04)**: rather than rigorous simulation experiments and verification,
> concentrate first on **building the front end (L0→L3)**. That is, complete "throw one sketch at it and does the
> physical system and the non-dimensionalization come out by procedure" first, and defer formalizing the golden
> tests (Phase 2) and the raw data and sweeps.
>
> **The direction reinforced (user instruction 2026-08-04, during the L3 work)** ⭐️: reduce the weight of the
> scientific and physical verification **after** a run — it differs per system and the system may not be one that has
> been reported (that is, in the literature). After a run, look only at **numerical errors**: divergence · NaN/Inf ·
> convergence to a strange value. Instead, **extracting the parameters from the sketch, doing the
> non-dimensionalization properly and setting the simulation up well** is far more important.
>
> → So L4 is built not as a "physics verifier" but as **a numerical-health judge**. No more analytic-comparison or
> literature-comparison machinery gets built (the `role` system of `bdbot.metrics` is already sufficient).
> → The investment in rigour is concentrated on the L0→L3 side. The `validate()` layer of `bdbot.nondim` (ledger
> completeness · whether a dimensionless number really is that ratio · invertibility) is the direct product of this instruction (§6.4).

**✅ The front end (completed 2026-08-04)**

| # | What | Artifact | Result |
|---|---|---|---|
| 1 | the `Observation` schema + `intake check` | [`bdbot/intake.py`](../../bdbot/intake.py) | the schema was **derived from the actual usage frequency across 5 files**. The `resolution` key is in all 47/47 and the value is often `null` — that is the device of §8.3. All 5 have 0 errors |
| 2 | making the `PhysicalSystem` loader common | [`bdbot/physical.py`](../../bdbot/physical.py) | the `derived_from` invariant · the tier gate (§12.4) · **derived values recomputed and verified** · the cross-check that "if L0 is BLOCKED, L2 cannot exist" |
| 3 | skill `bd-intake` | [`.claude/skills/bd-intake/SKILL.md`](../../.claude/skills/bd-intake/SKILL.md) | 8 rules against inventing. Distinguishing `physical`/`choice` · self-consistency inversion · handling dimensional ambiguity |
| 4 | the `bdbot` CLI | [`bdbot/cli.py`](../../bdbot/cli.py) | `status`·`intake`·`system`·`nondim`·`run`. It exposes the verdict as exit codes 0/1/2/3 (the surface a hook will read) |
| 5 | **L3 `NondimSpec` + `specs/`** | [`bdbot/nondim.py`](../../bdbot/nondim.py) | what made ③④⑤ of §6 checkable. **The only contract** between L2 and L4. 3 cases migrated · adversarial checks 30/30 · a 16-key schema common to the three cases |

**What the tools caught at L3** (§6.4 / `scratch/verify_l3_spec_gaps.py`)
- ⭐️ **`run_id` did not cover the physical system** — the 1-B spec had no physical system, so changing `d` 5µm→0.5µm and
  `η` by 62× (a 16.1× difference in τ_B) left run_id identical. `prepare_outdir` skipped it as "an already completed run"
  and **reported the old system's result as the new system's**. → `system` was put into the hash
- **whether a dimensionless number really is that ratio could not be checked** → `Group(num,den)` now points at the ledger
  and `validate()` recomputes and compares. **1 label error detected immediately on the first real use**:
  `U(d)` had been written as "A kT" when it is actually `(A+ε_WCA)kT` (at A=100 it is 101, a 1% discrepancy)
- **the reference scale was not in the ledger** — the report writes the reference length as `d`, while the ellipsoid case's
  ledger had only `d_eq`. The check caught the symbol mismatch
- **a scale missing from the ledger is a check that does not run** — 1-B had put `dt` and `T_obs` outside the ledger so
  they did not appear in the timescale alignment table. → the 4 required **roles** are enforced

**What the tools caught at the front end** (all of it came out of running the tools on real data)
- **the verdict was wrong** — 2 cases that had already been completed came out BLOCKED, and the cause was that
  `missing_required` **mixed physical unknowns with simulation choices**. → `kind: physical|choice` added
- **the checker crashed** — changing a unit to `furlong^2` gave a pint exception instead of an error report.
  → caught by building `verify_intake_guards.py` (27 adversarial tests)
- **`run_id` reacted to a documentation edit** — merely adding `derived_from` changed the 1-A run_id.
  → hash only the physics fields with `runid.physics_only()`. Re-ran and confirmed 96 fields identical
- **`derived_from` was only in a comment** — a provenance a machine cannot read. Promoted to a field

**What remains**

| # | What | Why it is deferred | Estimate |
|---|---|---|---|
| ~~4~~ | ✅ **the L4 numerical-health judge — done** | `Guard` (immediate runtime halt) · `judge_series` (divergence, NaN, freezing, collapse) · **`step_health` (the L3 ledger feedback)** · `gate` (before execution, rejecting a touched spec). Not a physics verifier | — |
| **4b** | the `chain-bend` · `trap-drag` case scripts ⬅️ **next** | L2 is READY and there is no script. **Split as principle 9.1 prescribes** — for `chain-bend`, ① the equilibrium chain (driving OFF) → ② add the driving | half a day each |
| 5 | Phase 2 formalizing the golden tests | the verification is already running as the 4 `scratch/verify_*.py`. `pytest` formalization comes after the front end is stable | 1 day |
| 6 | Phase 8-rest, hooks + slash commands | the CLI (#4) has to exist for there to be anything to enforce | 1 day |
| 7 | Phase 3 raw data Tier B~D | Tier A (GSD positions) alone has been sufficient so far | 1.5 days |
| 8 | the rest of Phase 4 | the LLM narrative · rule promotion (insufficient samples) · the approval ledger | 0.5 days |
| 9 | Phase 5 the SQLite KB + literature distillation | literature: **0 papers**. The "Eric Furst paper" of the `chain` sketch is the first seed candidate | 2.5 days |
| 10 | Phase 6 parallel sweeps | in 1-B, 7 runs were launched by hand | 1.5 days |

### C. What remains unconfirmed (enumerated honestly)

- `soft-r3`'s **reading that `A` is dimensionless** was inferred backwards from the sketch's goal (`system.yaml.not_verified`)
- **the finite-size CV3**: only `N=400` was run. The `N=900` comparison was not done
- **the initial-condition dependence**: only one RSA random placement. No comparison against a run started from a crystal
- **0 literature comparisons**: the literature `Γ` value of the 2D `r⁻³` melting transition was not checked
- **§20 A the anisotropic translational friction**: the premise of `abp-rod`. Option A (the isotropic average) looks like it can proceed, but it is unsettled
- ⭐️ **hypothesis verification is 0 across the whole project**: counting the observable roles gives `implementation_check` 10 ·
  `measurement` 9 · **`hypothesis` 0** (34 more unmarked in old runs). Everything done so far is a proof that
  "the pipeline is right", and **physical discoveries are 0.**
  That is normal at the pipeline-building stage, and **principle 9.1 and the user instruction above ("after a run, numerical
  health only") point at the same conclusion** — do not hold an existing theory up against a combined result. But
  if a discovery is ever wanted, the "assumptions the theory adds" must not be left empty when designing a case (principle 8).
- **`abp-rod`'s 5 implementation checks are still circular**: `derivation` was filled in to get past the principle 9.1
  gate, but its content is entirely "does the value I put in come back out". The gate **exposes the circularity
  but cannot block it** — that case needs a separate design of `hypothesis` observables.
- **The `PhysicsModule` registry is design only** (§5.6). Principle 9's `standalone_check()` is
  therefore unimplemented too. With 3 cases the module boundary is not yet settled.
- **Trap 14 (the 3D `rotational_diffusion` convention) is harmless for 2D cases only**. The moment 3D active
  work starts, it bites.
- What was **deliberately not lifted** into `bdbot`: the equilibrium indicator · the observables · the verification strategy ·
  the choice of the governing timescale · the initial placement · the sampling loop. If it comes up again in a third case, it goes in then

---

## 1. The settled scope

| Item | Decision |
|---|---|
| autonomy | **Claude Code is the agent runtime** (natural language + images → a simulation), **with human approval at every stage at first** |
| the execution environment | a local MacBook, CPU only (10 physical cores, 16 GB RAM) |
| the target systems | ⭐️ **not restricted.** No particular physical system is hardcoded; they are composed and extended through a **physics module registry**. If new physics is needed, one module file is added (with no core modification) → **[§5.6](#56-physicsmodule--the-unit-of-extension-)** |
| the outputs | a GSD trajectory (raw data included) + an observables table + automatic analysis plots |
| **knowledge** | **a literature-distillation KB + a KB of our own runs' experience, used as the grounds for a parameter proposal** |
| **units** | **dimensions first (a hard invariant). Every system is fixed as an SI physical system first, then non-dimensionalized by the main length, time and energy scales. Bidirectional conversion is mandatory.** |
| **raw data** | **positions, orientations, forces and torques stored hierarchically, re-analysed on demand** |

### Non-scope (not done now)
HPC/SLURM · multi-GPU · a web UI · autonomous parameter optimization (active learning) · polymers/rigid bodies/patchy particles

---

## 2. The core design principles

### Principle 1 — the LLM does not write code. It writes **the spec and the knowledge** ⭐️
```
❌ natural language/pictures → the LLM generates a HOOMD script → exec()
✅ natural language/pictures → the LLM generates a PhysicalSystem(JSON) → deterministic non-dimensionalization → SimSpec → the builder
```
Reproducibility (`run_id = sha256(SimSpec)`), verifiability, debuggability, safety and cost all come from here.

### Principle 2 — **every number has a provenance** ⭐️ (new)
Whether a parameter came from the KB or was estimated by the LLM, **where it came from is tracked**.

```python
particle_diameter = Provenanced(
    value=1.2 * ureg.micrometer,
    source="user_sketch:2026-08-03/note1.png#annotation_3",
    confidence=1,              # confirmed by a human
)
solvent_viscosity = Provenanced(
    value=1.0e-3 * ureg.Pa * ureg.s,
    source="kb:paper/10.1103-PhysRevLett.110.238301#table1",
    confidence=2,              # extracted from a paper, unverified
)
```
- The moment the LLM distils a paper, **a hallucination risk** appears. Provenance tracking is the only line of defence.
- A lineage table of "where each value of this spec came from" is always attached to the spec report.

### Principle 3 — **dimensions come first. Non-dimensionalization comes after.** ⭐️⭐️ (an invariant)

It is this project's **hard invariant**. No bypass route is provided.

```
        every system must pass through this order — no exceptions

PhysicalSystem (SI, pint Quantity, every field dimensional)
      │
      │  ① write the scale ledger — enumerate every length, time and energy scale present in the system
      │  ② choose the reference scales — explicitly designate the main length σ*, the main time τ* and the main energy E*
      │  ③ derive the dimensionless numbers — derive every dimensionless number as "the ratio of two scales" (the interpretation follows)
      │  ④ the scale-separation check — verify that a scale you mean to neglect really is separated
      ▼
SimSpec (reduced units) + DimensionlessReport + ScaleLedger
      │
      ▼  the simulation → dimensionless results
      │  ⑤ redimensionalize() — the inverse transform
      ▼
results in physical units (µm²/s, Pa, s, ...)
```

**Three enforced conditions:**

1. **A `SimSpec` cannot exist on its own.** Every `SimSpec` has a `derived_from: PhysicalSystemRef`, and
   the `PhysicalSystem` and `ScaleLedger` that produced it are stored alongside.
   There is **no** route to writing a dimensionless spec by hand — even if you want to start from dimensionless
   values, the physical system those values imply has to be stated first (§6.6, the reverse-construction route).

2. **The reference scales are not decided automatically.** Which length and which time were taken as the reference
   is recorded in `ScaleLedger.reference` **together with the grounds for the choice**, and it becomes the subject of human check #3.

3. **Without the inverse transform it cannot go into the KB.** Dimensionless results alone cannot be compared against
   the literature or experiment. `observables.parquet` **always** stores the dimensionless value and the physical-unit value **as a pair**.

> What observing this principle buys: a dimensionless number becomes "the ratio of two physical scales" rather than
> "a number from somewhere", so the physical interpretation follows automatically, and the scale-separation check (§6.4)
> — a powerful automatic verification — becomes possible.

### Principle 4 — a human approves, and the approval history is the grounds for autonomy (new)
At first a human confirms at each of `Intake → PhysicalSystem → SimSpec`.
The approve/reject/amend history is left in an **approval ledger**.
Later, "this type was approved unamended 30 times out of 30 → an automation candidate" is judged from data.

### Principle 5 — failures fast and loud, and **recorded**
Verification before execution / monitoring during execution / a **structured post-mortem** after execution. A failure experience becomes the next verification's rule.

### Principle 6 — face the MacBook CPU reality
A HOOMD CPU run is **one run = one core**. Parallelism comes from "running several runs at once".
The default is 8 workers, with a realistic scale of N = 10³~2×10⁴ and 10⁶~10⁸ steps.

---

### Principle 7 — **verify a physical claim before stating it** ⭐️ (new, v0.4)

In this project there have been **two** cases of asserting something by reasoning and having the measurement overturn it:

| Claim | Measured |
|---|---|
| "the trap verification passes" (only k=10 checked) | **+1856%** at k=2 — a missing minimum image (§11 trap 7) |
| "bind it as a rigid body and the friction becomes anisotropic" | `γ⊥/γ∥ = 1.000000` — BD has no HI |

Both were **plausible and both were wrong.** A claim about HOOMD behaviour or a physics result is confirmed
**by execution**, not by reasoning. If it could not be confirmed, write "not confirmed"
— §10's `record.json.not_verified` field enforces that.

### Principle 8 — **"verification" does not mean agreeing with the prior hypothesis** ⭐️⭐️ (new, 2026-08-04)

The reason to compute the systems in `intake/` is that **they may differ from the prior hypothesis**.
And yet writing the verdict logic as "differs from the prediction, FAIL" **calls a discovery a failure.**

Attach a role to every comparison (`bdbot.metrics.ROLES`, `judge()`):

| Role | Where the prediction comes from | A mismatch |
|---|---|---|
| `implementation_check` | **derived analytically from the model implemented** | **a bug** → FAIL |
| `hypothesis` | an assumption the simulation does **not** impose (continuum, dilute limit, effective medium, a literature model) | **a result** — not a FAIL |
| `measurement` | none | the simulation is the answer |

**It actually bit.** `abp-rod`'s 5 predictions (`τ_eff`, `D_eff`, `D̄`, `τ`, the tumble frequency) agreed to
within 0.66%, but all five were derivable from the model I implemented (the active force + the rotational-diffusion
updater + the Poisson tumble + BD). Agreement is code verification and not a physical discovery
— it is nearly circular. **That case's hypothesis verification was 0.**

So at the case-design stage two lists are written separately:

```
the assumptions the simulation imposes   overdamped BD · isotropic translational friction · independent rotational diffusion · a Poisson tumble · dilute
the assumptions the theory adds          ← if this is empty there is no discovery
```

**The roles are fixed before looking at the results.** Change a role after seeing the results (demoting a
`hypothesis` mismatch to a `measurement`, or promoting a coincidental agreement to an `implementation_check` success)
and this principle is neutered. For the timing and configuration-level rules see **principles 9.1 and 9.2**.

This principle also changes the meaning of §12 (the verification layers) and §16 Phase 2 (the golden tests) —
**a golden test covers `implementation_check` only, and a physical discovery needs a separate design.**

### Principle 9 — **isolate independent elements and verify them one at a time** ⭐️⭐️ (new, 2026-08-04)

If a system has several mutually independent elements, first run the **minimal configuration with only each
element on**, confirm it is right on its own, and then combine.

**Why — three reasons, and the first is decisive**

1. **The isolated configuration often has an analytic solution. Combined, it does not.**
   Isolation is sometimes not a debugging convenience but **the only way to make a ground truth**.
2. **When it is wrong, the culprit is identified.** Wrong in a combination and you have to suspect N things at once,
   and you cannot even tell a standalone bug from an interaction effect.
3. **What appears only in the combination can be identified as a real interaction.**
   Without the standalone verification, `cross_check` (§5.6) has no meaning.

**The evidence — actually, in this project**

| Standalone verification | The analytic solution | Result |
|---|---|---|
| `external.harmonic_trap` (interactions OFF, N independent particles) | `⟨x²⟩=kT/k` · `τ=γ/k` · a Lorentzian | all 4 observables agree |
| `shape.rigid_rod` friction (a constant force → the terminal velocity) | `v=F/γ` | `γ⊥/γ∥ = 1.000000` — **the expectation overturned** |
| the integrator bias (a linear system) | `bias=(dt/τ)/2` | 0.07% against theory |
| `active.abp` (pair, trap and shape OFF) | `D_eff = D_t + v₀²/(d·Λ)` | **2 traps discovered** ↓ |

What came out of the `active.abp` standalone verification (`scratch/standalone_abp_diffusion.py`):

- **If the active force is 0 the rotational diffusion does not operate at all** — Λ=0 (in all 4 combinations).
- **HOOMD's `rotational_diffusion` is the director decay rate Λ itself.**
  In both 2D and 3D, `Λ/D_r = 1.00`. Not the standard theory's `(d−1)D_r` but **off by a factor of 2 in 3D.**
- So the `D_eff = D_t + v₀²/[2(d−1)D_r]` that was in §17 **was only coincidentally right in 2D**
  (a +29~31% discrepancy measured in 3D). The correct form is `D_t + v₀²/(d·Λ)` — measured to within 1.5%.

**Neither would have been caught inside a combined simulation.** Mix the active force, the interactions and the trap
and they all get buried in one MSD curve. In fact at first only `D_eff` was measured and it was "is the expression wrong",
and **only after measuring Λ and v₀ separately** did it emerge that the cause was on the rotational-diffusion side.

**The procedure**

1. When adding a module M, first run **the minimal case with only M on**
2. Compare against the analytic solution or a known limit. If the comparison fails, **split further**
   (ABP → Λ separately, v₀ separately)
3. It has to pass before going into a combination
4. If a combination fails, **first look at whether the relevant modules' standalone verifications are current** (a regression)

`standalone_check()` goes on `PhysicsModule` (§5.6) — the minimal configuration · the prediction · the tolerance.

### 9.1 **The roles differ per stage** ⭐️⭐️ — do not hold an existing theory up against a combination

The two stages of principle 9 have different epistemological status. Mix them and principle 8 gets broken.

| Stage | Role | Why | What a mismatch means |
|---|---|---|---|
| **a module alone** | `implementation_check` | in the minimal configuration **the existing theory is my model**. Free ABP's `D_eff` is derived from the equation I implemented | **a bug** |
| **the whole combination** | `measurement` or `hypothesis` | a combination usually has **no** analytic solution — that is why it is being simulated | **a result** |

**Why holding an existing theory up against a combined result is dangerous:**

- if it matches you end up saying "it is verified", when in fact you have not confirmed that the theory was not imposed
- if it does not match you end up saying "that theory does not hold in this region", and that is unfalsifiable
- **either way nothing is learned, and at worst a discovery gets erased**

To use an existing theory at the combination stage, register it as a `hypothesis` and **report a mismatch as a result**.
To use `composite` + `implementation_check`, **why that analytic expression is also derivable for the combination**
(`derivation`) must be written down — usually it is a limit (dilute, linear response, short time).

### 9.2 The prediction is fixed **before the results are seen** ⭐️

Fit the interpretation after seeing the results and there is no telling whether it was fitted or was right.

**A structural guarantee** (rather than self-declaration): the prediction function **does not take the simulation
results as an argument.** The `analytic(lg)` of `cases/*.py` is already in that form — it takes only the scale ledger,
so it has no way to depend on the results. It is dumped as `predictions.json` before execution so it can be audited.

**If the interpretation has to change after seeing the results** — that can happen. But conditions attach:

> The new interpretation has to be independently confirmed by **an observable other than the one that failed**.
> And the fact that the original prediction failed is left in the record.

**A real case (2026-08-04, the `active.abp` standalone verification)** — why this rule is needed:

| Stage | What happened |
|---|---|
| pre-registration | `D_eff = D_t + v₀²/[d(d−1)D_r]` (derived before execution) |
| result | **a +90% discrepancy in 3D — the registered prediction failed** |
| tracing the cause | `Λ` (the director decay rate) measured independently **from the autocorrelation, not the MSD** → `Λ = D_r` |
| re-prediction | `D_eff = D_t + v₀²/(d·Λ)` → −1.5~−0.6% across 4 conditions |

The grounds that it is not curve fitting: ① `Λ` came from **a different observable** ② after putting it in there are 0 free parameters
③ it matched 4 conditions simultaneously (had it been fitted, one point would be 0% and the rest worse).

**The real lesson is elsewhere.** What was wrong was not the physics expression but **an unregistered parameter mapping**.
It was registered as `D_eff = f(D_r)` while the actual structure is `D_eff = f(Λ)`, with one
**invisible assumption** slipped in: that `Λ = (d−1)D_r`. HOOMD had `Λ = D_r`.

> **When registering a prediction, register the parameter-mapping assumptions with it.**
> "The `D_r` of this expression is the same as the tool's `rotational_diffusion`" is an assumption too.

**The limits — honestly**

- Being right alone does not make the combination right. `cross_check` is needed separately.
- **Some things are inseparable in principle.** `shape` has no defined friction without `medium`.
  The "minimal configuration" has to be defined by following the dependency DAG of §5.6.
- **Some combinations become physically meaningless when separated.** Trap ① above is exactly that case —
  trying to turn off the active force and see only the rotational diffusion made HOOMD turn the rotational diffusion off.

**The roadmap impact** — two of the remaining cases put 2 modules in at once:

| Case | Modules added | Standalone-ness | Proposal |
|---|---|---|---|
| `trap-2d-5um` ✅ | external | fully standalone | — |
| `soft-r3` ✅ | pair | the pair only | — |
| `trap-drag` | external (verified) + driving | effectively driving alone | as it is |
| `chain-bend` | bonded + driving | **2 at once** | ① the equilibrium chain (driving OFF) → ② add the driving |
| `abp-rod` | shape + active | **2 at once** | ① a spherical ABP ✅ **done** → ② add the shape |

## 3. The architecture — the closed loop

```
┌── L0 INTAKE (multimodal) ────────────────────────────────────────┐
│  a sketch image / a hand note / a whiteboard photo / a paper PDF / natural language      │
│  → Claude Code Read tool → Observation (structured) ─ [human check #1] ►│
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L1 KNOWLEDGE ────────────────────────────────────────────────┐
│  KB search: similar systems · the dimensionless-coordinate neighbourhood · our past runs      │
│  proposing the missing parameters on the strength of the results (with the source and confidence stated)              │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L2 PHYSICAL SYSTEM (dimensional, SI) ──────────────────────────┐
│  σ=1.2µm, T=298K, η=1mPa·s, φ=0.6, v₀=8µm/s, ...              │
│  a Provenance attached to each field   ─── [human check #2, approval per field] ──► │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L3 NONDIMENSIONALIZATION ────────────────────────────────────┐
│  ① the scale ledger  ② choose the reference (d, τ_B, kT)  ③ derive as ratios  ④ the separation check │
│  → SimSpec + DimensionlessReport + ScaleLedger                 │
│  Pe=38.4, φ=0.60, D_r*=3.0, dt*=5e-5 ── [human check #3] ───────► │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L4 VALIDATION ── physics, numerics, cost + KB-based warnings ───────────────┐
│  "in the past, the Pe>80, dt=1e-4 combination diverged 3 times"                     │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L5 BUILD → L6 RUN (a process pool of 8) → L7 RAW DATA STORE ───────┐
│  GSD trajectory (positions/orientations/forces/torques, hierarchical frequencies) + HDF5 log   │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L8 ANALYZE (on demand) → L9 PLOT / literature comparison ─────────────────┐
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L10 POST-MORTEM ── structuring the success/failure factors ────────────────────┐
│  → KB feedback (improving the next proposal and the verification rules)  ──────► back to L1  │
└────────────────────────────────────────────────────────────────┘
```

**Important**: L3~L9 have to work on their own through the CLI / the Python API with no LLM. The LLM handles only the "interpretation and judgment" parts of L0~L2 and L10. Only that way can the pipeline be verified without an LLM.

---

## 4. The directory structure

### 4.1 Current (2026-08-04) — what actually exists, driven by the cases

```
simulation_auto/
├── CLAUDE.md                      the invariants (always loaded)
├── mater_plan.md · environment.yml
├── .claude/skills/{bd-hoomd,bd-physics}/SKILL.md
├── docs/hoomd_capabilities.md     HOOMD capability matrix, measured
├── intake/<case>/                 the sketch + observation.yaml + system.yaml
├── bdbot/                         ⭐️ Phase 1-C output — **only what appeared twice**
│   ├── units.py                   a single pint registry
│   ├── provenance.py              Provenanced (value + source + tier)
│   ├── materials.py               γ=3πηd · D_t=kT/γ · τ_B · m · τ_p
│   ├── scales.py                  ScaleLedger + thermal reference
│   ├── checks.py                  Check(model/integration/geometry/statistics) + dt = 10⁻²·γ/(local stiffness)
│   ├── report.py                  DimensionlessReport renderer
│   ├── runid.py                   content addressing + preventing a re-run
│   ├── metrics.py                 metrics.json schema
│   ├── stats.py                   block averaging · the autocorrelation correction · the unbiased autocorrelation
│   └── sim.py                     2D frame · BD integrator · GSD · seeds · minimum image
├── cases/{trap_2d_5um,soft_r3_2d}.py   the end-to-end scripts (only the case-specific physics is left)
├── tools/{postmortem,kb}.py       the post-mortem + KB queries
├── scratch/*.py                   11 verification scripts (reproducible)
└── runs/<run_id>/                 report·result·metrics·record·traj·plots
```

`bdbot` **has no CLI yet** — the cases import `bdbot` directly. The CLI is right to build together with the
Phase 8-rest hooks, which need something to enforce (`bdbot run`).

**What was deliberately not included**: the equilibrium indicator · the observables · the verification strategy · the choice of the governing timescale ·
the initial placement · the sampling loop. They differed per case (skill `bd-physics` §6.3).

### 4.2 The target (the form it will converge to after Phase 1-C)

The below is **the design intent** and not the current implementation state. Refer to it when extracting the
common denominator of 1-A/1-B in Phase 1-C, but build **only what actually appeared twice**.

```
simulation_auto/
├── mater_plan.md
├── environment.yml · pyproject.toml
│
├── bdbot/
│   ├── intake/                    # ⭐️ L0
│   │   ├── vision.py              # an image → Observation (Claude vision)
│   │   ├── pdf.py                 # a paper PDF → distillation candidates
│   │   ├── observation.py         # Observation schema
│   │   └── review.py              # the human-check UI (a terminal form)
│   │
│   ├── knowledge/                 # ⭐️ L1, L10
│   │   ├── store.py               # SQLite + FTS5
│   │   ├── schema.py              # KnowledgeEntry, Claim, Source, Provenance
│   │   ├── distill.py             # a paper → KnowledgeEntry (LLM)
│   │   ├── search.py              # keyword + dimensionless-coordinate neighbourhood search
│   │   └── feedback.py            # a run post-mortem → fed back into the KB
│   │
│   ├── physical/                  # ⭐️ L2
│   │   ├── system.py              # PhysicalSystem (pint units attached)
│   │   ├── materials.py           # a constants library: water/glycerol viscosity, silica/PS density and so on
│   │   └── provenance.py
│   │
│   ├── nondim/                    # ⭐️ L3
│   │   ├── scales.py              # the scale-selection strategy
│   │   ├── forward.py             # PhysicalSystem → SimSpec
│   │   ├── inverse.py             # dimensionless results → physical units
│   │   ├── groups.py              # Pe, φ, Re, Pé_r, λ_D/σ ... dimensionless numbers computed
│   │   └── report.py              # DimensionlessReport
│   │
│   ├── spec/                      # SimSpec (dimensionless)
│   ├── validate/                  # L4 (+ KB linked warnings)
│   ├── build/                     # L5 (HOOMD assembly)
│   ├── run/                       # L6 (the pool, the queue, guards, checkpoints)
│   ├── rawdata/                   # ⭐️ L7
│   │   ├── policy.py              # the storage tier policy (tier A/B/C)
│   │   ├── writers.py             # GSD/Burst/per-particle loggers
│   │   ├── index.py               # the frame index (fast random access)
│   │   └── loader.py              # lazy loading, memory mapping, slicing
│   │
│   ├── analyze/                   # L8 (a registry + a cache)
│   ├── plot/                      # L9 (+ a literature-comparison overlay)
│   ├── postmortem/                # ⭐️ L10
│   │   ├── taxonomy.py            # the failure classification system
│   │   ├── diagnose.py            # automatic diagnosis (based on numerical indicators)
│   │   └── narrate.py             # LLM summary → KnowledgeEntry
│   │
│   ├── approval/                  # ⭐️ the approval ledger
│   │   ├── ledger.py
│   │   └── gates.py
│   │
│   ├── distill_batch.py           # (optional) bulk distillation via the anthropic Batch API — §15.9
│   └── cli.py                     # ⭐️ the only entry point. Claude Code calls only this.
│
├── CLAUDE.md                      # ⭐️ the always-loaded invariants (§15.4)
├── .claude/                       # ⭐️ Claude Code integration surface (§15.3)
│   ├── settings.json              #   permissions + hooks
│   ├── skills/                    #   bd-physics, bd-hoomd, bd-intake, bd-distill
│   ├── agents/                    #   bd-distiller, bd-analyst, bd-reviewer
│   ├── commands/                  #   /bd-intake, /bd-spec, /bd-run, ...
│   └── hooks/                     #   guard_invariant, guard_separation, guard_cost
│
├── kb/
│   ├── knowledge.db               # SQLite (FTS5)
│   ├── papers/                    # the original PDFs
│   └── figures/                   # the extracted figures
│
├── intake/<case>/                 # the sketches/notes (thrown in by the user)
│   ├── sketch_01.jpeg
│   ├── observation.yaml           # L0 output (settled after human confirmation)
│   └── system.yaml                # L2 PhysicalSystem (settled after human confirmation)
├── specs/                         # L3 output (generated by bdbot nondim, never written by hand)
├── runs/<run_id>/                 # the raw data + the analysis results
└── tests/
```

> That **`cli.py` is the only entry point** matters. It has to work identically outside a Claude Code session
> (cron, a script, another person), and that is the defence against the single risk §15.10 pointed out.

---

## 5. The data model (5 core objects)

### 5.1 `Observation` — the L0 output
```python
class Observation(BaseModel):
    """What was read out of a sketch/note/picture. The LLM fills it and a human confirms."""
    source_files: list[str]
    raw_transcription: str              # the letters read from the image, verbatim
    system_guess: str                   # "2D active colloid monolayer"
    entities: list[Entity]              # particles, walls, fields, arrows …
    stated_quantities: list[StatedQuantity]  # "d ≈ 1 µm", "Pe 20~100"
    stated_goals: list[str]             # "MIPS phase boundary, find it"
    ambiguities: list[str]              # ⭐️ LLM states what it could not be sure of
    unread_regions: list[str]           # ⭐️ state what could not be read
```
`ambiguities` and `unread_regions` are **required fields**. A device for making it say when it does not know.

### 5.2 `KnowledgeEntry` — the KB's unit record
```python
class Source(BaseModel):
    kind: Literal["paper", "book", "our_run", "user_input", "handbook"]
    doi: str | None; title: str | None; authors: list[str]; year: int | None
    locator: str | None                 # "Fig.3b", "Table 1", "p.4 eq.(7)"
    local_path: str | None

class Claim(BaseModel):
    statement: str                      # "MIPS onset near Pe≈35 at φ=0.6 (2D ABP)"
    dimensionless_coords: dict[str, float]   # {"Pe": 35, "phi": 0.6}
    kind: Literal["parameter", "phase_boundary", "scaling", "method_note", "pitfall"]

class KnowledgeEntry(BaseModel):
    id: str
    system_tags: list[str]              # ["2D", "ABP", "WCA", "monodisperse"]
    source: Source
    claims: list[Claim]
    physical_params: dict               # the dimensional values (with the unit string included)
    dimensionless_params: dict          # Pe, φ, D_r*, ...
    numerics: dict                      # the dt, N, integrator and box that paper used
    confidence: Literal[0, 1, 2, 3]     # the table below
    extracted_by: str                   # "claude-opus-5 / distill_prompt_v3"
    verified_by: str | None; verified_at: datetime | None
    notes: str
```

**The confidence tier** ⭐️
| tier | Meaning | Usable alone in a production run |
|---|---|---|
| 0 | entered directly by a human / a property handbook | ✅ |
| 1 | extracted from the literature + human-verified | ✅ |
| 2 | extracted from the literature, unverified (the LLM alone) | ⚠️ needs approval |
| 3 | induced from our own simulation | ⚠️ needs approval |

The validator always requires human approval for "a spec composed only of tier 2 or lower values".

### 5.3 `PhysicalSystem` — the dimensional physical system
```python
class PhysicalSystem(BaseModel):
    label: str
    dimensions: Literal[2, 3]

    # the particles
    particle_diameter: Provenanced[Quantity]        # µm
    particle_density: Provenanced[Quantity] | None  # kg/m³
    n_particles: Provenanced[int] | None
    area_or_volume_fraction: Provenanced[float]

    # the medium
    temperature: Provenanced[Quantity]              # K
    solvent_viscosity: Provenanced[Quantity]        # Pa·s
    # → γ = 3πηd (Stokes) derived, or specified directly
    drag_coefficient: Provenanced[Quantity] | None

    # the interactions
    interaction: Provenanced[InteractionSpec]       # WCA/LJ/Yukawa + ε(J or kT), κ(1/m)

    # active
    self_propulsion_speed: Provenanced[Quantity] | None   # µm/s
    rotational_diffusion: Provenanced[Quantity] | None    # 1/s

    # the observation goal
    target_observables: list[str]
    target_physical_time: Provenanced[Quantity]     # s (how long to watch for)
```
Every field is `Provenanced` → the value + the source + the confidence.
Units are enforced with `pint` to block unit mistakes.

### 5.4 `SimSpec` — the dimensionless spec (straight into HOOMD)
```python
class SimSpec(BaseModel):
    schema_version: str = "0.2"
    label: str

    # ⭐️ the invariant: a dimensionless spec is always derived from a physical system
    derived_from: PhysicalSystemRef          # the hash + the storage path
    scale_ledger: ScaleLedger                # the reference scales + the full scale list
    dimensionless: dict[str, float]          # {Pe, phi, D_r_star, dt_star, Re, St, ...}

    box: BoxSpec
    types: list[ParticleType]
    modules: list[ModuleSpec]                # ⭐️ the combination of physics modules (§5.6)
    integrator: IntegratorSpec               # all in reduced units
    run: RunSpec
    raw_data: RawDataPolicy
    observables: list[ObservableSpec]
```
A `SimSpec` with no `derived_from` is **rejected by the builder**. To start from dimensionless values,
the physical system has to be stated first, through the reverse-construction route of §6.6.

### 5.5 `RunRecord` + `PostMortem`
```python
class PostMortem(BaseModel):
    run_id: str
    outcome: Literal["success", "partial", "failure"]
    failure_modes: list[FailureMode]    # the taxonomy below
    diagnostics: dict                   # the energy drift, the equilibrium verdict, the finite-size indicators …
    narrative: str                      # LLM summary
    lessons: list[Claim]                # ⭐️ KB feedback entries
    dimensionless_coords: dict          # which parameter region it was in
```

---

### 5.6 `PhysicsModule` — the unit of extension ⭐️⭐️

**Since the scope was decided not to be restricted, physical systems are not hardcoded but composed from modules.**
If new physics is needed, **one module file is added** and registered in the registry. The core is not touched.

The crux is that one module **contributes to all 7 layers by itself**. Only then does attaching new physics bring
the scale ledger, the non-dimensionalization, the separation checks and the verification along automatically (principle 3 holds under extension too).

```python
class ModuleContribution(BaseModel):
    scales:  dict[str, Quantity]          # ① the characteristic scales to add to the ledger (§6.1)
    groups:  list[DimensionlessGroup]     # ② the dimensionless numbers this module creates (§6.3)
    checks:  list[SeparationCheck]        # ③ the separation checks this module requires (§6.4)
    reduced: dict[str, float]             # ④ the dimensionless parameters

class PhysicsModule(ABC):
    kind:     ClassVar[str]               # "external.harmonic_trap"
    requires: ClassVar[set[str]] = set()  # the modules depended on (e.g. {"shape.*"})
    PhysicalParams: ClassVar[type[BaseModel]]   # the dimensional parameter schema
    ReducedParams:  ClassVar[type[BaseModel]]   # the dimensionless parameter schema

    @abstractmethod
    def contribute(self, phys, ctx: ScaleContext) -> ModuleContribution: ...  # L3
    @abstractmethod
    def build(self, sim, spec, ctx: BuildContext) -> None: ...                # L5

    def standalone_check(self) -> StandaloneCheck | None: ...  # ⭐️ the principle 9 standalone verification
    def cross_check(self, others) -> list[Issue]: return []   # L4 combination check
    def default_observables(self) -> list[str]: return []     # L8
    def periodic_safe(self) -> bool: return True              # ⭐️ trap 7 (§11)

class StandaloneCheck(BaseModel):
    """The minimal configuration with only this module on + the analytic solution + the tolerance (principle 9)."""
    minimal_spec: dict          # every other module OFF
    predictions: dict           # observable → the analytic expression
    tolerance_pct: float
    script: str | None = None   # the path of the reproduction script
```

#### The resolution order (the dependency DAG)

```
① medium.*   →  η, T, kT, ρ_f
② shape.*    →  the friction tensor γ∥, γ⊥, γ_r  →  D_t, D_r     ← the sphere/ellipsoid branch point
③ the rest   →  each contributes its scales, dimensionless numbers and checks
④ choose the reference scales (§6.2)
⑤ non-dimensionalize + the whole separation check (§6.4)
```

Why `shape` has to be resolved first: **`D_r = 3D_t/d²` holds only for a sphere.** That is not a core formula but
a contribution of the `shape.sphere` module. `shape.ellipsoid` contributes the Perrin friction factors instead.
Set up this way, a case like `abp-rod` comes in with no core modification.

> ⚠️ **The measured result (2026-08-03)**: whether spheres are bound as a rigid body (`constrain.Rigid`) or made into
> a bead chain with bonds, **the translational friction is isotropic** (`γ⊥/γ∥ = 1.000000`, on both routes). A rod's
> `γ⊥/γ∥ → 2` is **an effect of hydrodynamic interactions (HI)** and BD has no HI — it is not a limitation of HOOMD but
> **a property of the BD model itself**. Rotational friction works correctly via the `gamma_r` tensor (a ratio of 0.2000 measured).
> → the `shape.*` module declares `translational_friction: "isotropic_average" | "anisotropic"`, and
> the validator raises the warning "short-time MSD anisotropy is the goal but an isotropic approximation is in use".
> The grounds and the options: [`docs/hoomd_capabilities.md` §5.1–5.4](../../docs/hoomd_capabilities.md)

#### What a module contributes — examples

| Module | Scales added | Dimensionless numbers | Separation checks |
|---|---|---|---|
| `shape.sphere` | `d` | — | — |
| `shape.ellipsoid` | `a` (the semi-major axis) `b` (the semi-minor axis) | `p = a/b`, the aspect ratio | ⚠️ anisotropic translational diffusion **unsupported in BD** (below) |
| `pair.wca` | `r_c` | `ε*` | `r_c < L/2` |
| `pair.table` | `r_c`, the potential scale | `A*` | `dt/τ_int ≤ 1e-2` |
| `external.harmonic_trap` | `ℓ_k=√(kT/k)`, `τ_k=γ/k` | `k* = k d²/kT` | `dt/τ_k ≤ 1e-2`, `ℓ_k < L/2` |
| `bonded.bond_harmonic` | `ℓ_b=√(kT/k_b)`, `τ_b=γ/k_b` | `k_b*` | `dt/τ_b ≤ 1e-2` |
| `bonded.angle_harmonic` | — | `κ* = κ_bend/kT` (the persistence length) | `ℓ_p^chain ≤ L/4` |
| `active.abp` | `ℓ_p=v₀/D_r`, `τ_v`, `τ_r` | `Pe`, `D_r*` | `dt·D_r ≤ 1e-2`, `ℓ_p ≤ L/4` |
| `active.run_and_flip` | `τ_flip` | `Pe`, `τ_flip/τ_B` | `dt/τ_flip ≤ 1e-2` |
| `driving.oscillate` | `τ_ω = 1/ω` | `De = τ_relax·ω` (Deborah) | `dt/τ_ω ≤ 1e-2`, `T_obs ≫ τ_ω` |

> Each row means that "attach new physics and the scale ledger, the dimensionless numbers and the separation checks grow **automatically**".
> Principle 3's benefit holds under extension unchanged.

#### `ModuleSpec` — its representation inside the spec

```python
class ModuleSpec(BaseModel):
    kind: str                    # the registry key. "external.harmonic_trap"
    params: dict                 # that module's ReducedParams
    targets: str = "all"         # the target filter (a hoomd.filter expression)
```

The registry is a decorator registration in `bdbot/modules/__init__.py`:
```python
@register("external.harmonic_trap")
class HarmonicTrap(PhysicsModule): ...
```

#### The current status — measurement-based

Which modules can be built rests not on guesswork but on the survey result:
**[`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md)** — measured on hoomd 7.1.0,
15 APIs verified working, and all 5 `intake/` cases confirmed implementable.

---

## 6. The non-dimensionalization engine (L3) — in detail

It implements principle 3's 5 stages (① the scale ledger → ② choose the reference → ③ derive as ratios → ④ the separation check → ⑤ the inverse transform).

### 6.1 ① The Scale Ledger ⭐️

The first stage of non-dimensionalization is **enumerating, without omission, every characteristic scale present in the
system in SI units**. The order of computing the dimensionless numbers first and inferring the scales afterwards is forbidden.

```python
class ScaleLedger(BaseModel):
    lengths: dict[str, Quantity]     # name → the SI value
    times:   dict[str, Quantity]
    energies: dict[str, Quantity]
    reference: ReferenceScales       # the set chosen as the reference + the grounds for the choice
    separations: list[SeparationCheck]
```

**The scales that must be computed and entered in the ledger** (recorded explicitly as `None` if not applicable):

| Kind | Name | Definition | Meaning |
|---|---|---|---|
| **length** | `d` | the particle diameter | the base length |
| | `a_mean` | `ρ^(-1/dim)` = the mean interparticle spacing | the dense/dilute verdict |
| | `r_c` | the interaction cutoff | the reach of the force |
| | `λ_D` | the Debye length (a charged system) | electrostatic screening |
| | `ℓ_p` | `v₀/D_r` (active) | the orientational persistence length |
| | `L` | one side of the box | the system size |
| | `ξ` | the correlation length (post hoc, a measured value) | the finite-size verdict |
| **time** | `τ_p` | `m/γ` | momentum relaxation (inertia) |
| | `dt` | the integration time step | the numerical resolution |
| | `τ_v` | `d/v₀` (active) | the advective/ballistic time |
| | `τ_int` | `d²γ/ε` | the interaction response time |
| | `τ_r` | `1/D_r` | the orientational correlation time |
| | `τ_B` | `d²/D_t` | the diffusive (Brownian) time |
| | `T_obs` | the production physical time | the observation window |
| **energy** | `k_BT` | the thermal energy | the fluctuation scale |
| | `ε` | the interaction depth/strength | the coupling strength |
| | `f_a·d` | the self-propulsion work | the active driving |

The basic property relations (provided by `physical/materials.py`, all dimensional):
```
γ    = 3πηd                      (Stokes drag, a sphere)
D_t  = k_BT/γ                    (Stokes–Einstein)
D_r  = k_BT/(πηd³) = 3D_t/d²     (Stokes–Einstein–Debye, a sphere)
m    = ρ_p (π/6) d³              (the particle mass)
v₀   = f_a/γ                     (the self-propulsion speed)
```

### 6.2 ② Choosing the reference scales — explicit, with the grounds recorded

```python
class ReferenceScales(BaseModel):
    length: tuple[str, Quantity]     # ("d", 1.2 µm)
    time:   tuple[str, Quantity]     # ("tau_B", 3.96 s)
    energy: tuple[str, Quantity]     # ("kT", 4.11e-21 J)
    rationale: str                   # why this was chosen (the subject of a human check)
```

| Strategy | Length | Time | Energy | When |
|---|---|---|---|---|
| `thermal` **(default)** | `d` | `τ_B = d²/D_t` | `k_BT` | thermal fluctuation dominates — colloids, ABPs |
| `interaction` | `d` | `τ_int = d²γ/ε` | `ε` | the interaction dominates — strong coupling, low temperature |
| `active` | `ℓ_p` | `τ_r` | `f_a·ℓ_p` | activity dominates — a very large Pe |
| `custom` | user-specified | | | when following a literature convention |

Under the `thermal` strategy the HOOMD spec becomes `σ=1, kT=1, γ=1` → `τ_B = 1` and `D_t = 1`.

⚠️ **Change the reference and the same physical system has entirely different dimensionless numbers.** When comparing
dimensionless numbers against the literature, which reference the other paper used has to be checked, and the reference
is recorded in the `KnowledgeEntry` too.

### 6.3 ③ A dimensionless number is derived as "the ratio of two scales" ⭐️
Rather than throwing a number, **which two scales the ratio is of** is recorded with it. The physical interpretation follows.

| Dimensionless number | The scale ratio | Expression | Physical meaning |
|---|---|---|---|
| `φ` | the occupied volume / the total volume | `N v_p / V` | the packing |
| `Pe` | `τ_B / τ_v` | `v₀d/D_t = f_a d/k_BT` | advection vs diffusion |
| `D_r*` | `τ_B / τ_r` | `D_r τ_B` (spherical Stokes → 3) | rotational vs translational diffusion |
| `ℓ_p/d` | `ℓ_p / d` | `Pe / D_r*` | how many particles the persistence length is |
| `ε*` | `ε / k_BT` | | the coupling strength vs thermal fluctuation |
| `κ*` | `d / λ_D` | `κd` | screening vs the particle size |
| `dt*` | `dt / τ_B` | | the numerical resolution |
| `St` | `τ_p / τ_B` | `m/(γτ_B)` | inertia vs diffusion |
| `Re` | inertia vs viscosity (the fluid) | `ρ_f v₀ d / η` | the fluid inertia |
| `N*` | `L / d` | `L/d` | how many particles the system is |
| `T*` | `T_obs / τ_B` | | how many diffusion times the observation window is |

Each dimensionless number in `groups.py` is stored as `DimensionlessGroup(name, value, numerator_scale, denominator_scale, expression, interpretation)`, so the report automatically explains "what against what".

### 6.4 ④ The scale-separation check ⭐️⭐️ — a hard gate

The real value of non-dimensionalization is here. **Whether a scale you decided to neglect really is separated** is verified automatically.
Sort the ledger's time scales by magnitude and a violation becomes visible.

**The time-scale rule (unified into one)**: the integration step has to be **at or below 1% of the fastest physical
timescale**, and the inertial relaxation has to be **at least 100× faster than the integration step**. Every threshold
is unified at `10⁻²` so it is easy to remember.

| Check | Condition | Equivalent expression | On violation |
|---|---|---|---|
| **neglecting inertia (model)** | `τ_p / τ_dyn ≤ 10⁻²` | `τ_dyn` = the fastest timescale of interest | ❌ hard — overdamped BD invalid, Langevin needed. **Not compared against `dt`** — BD has no inertia so `dt≫τ_p` is unnecessary (measured: 0.38% accurate even at `τ_p/dt=4000`) |
| **advection resolved** | `dt / τ_v ≤ 10⁻²` | `v₀ dt ≤ 0.01 d` | ❌ hard — moving more than 1% of a diameter in one step |
| **rotation resolved** | `dt / τ_r ≤ 10⁻²` | `dt·D_r ≤ 0.01` | ❌ hard — the ABP orientational dynamics collapses |
| **interaction resolved** | `dt / τ_int ≤ 10⁻²` | `dt·ε/(d²γ) ≤ 0.01` | ❌ hard — the force integration is inaccurate |
| **low Reynolds** | `Re ≤ 10⁻³` | | ❌ hard — the fluid inertia cannot be neglected |
| **the cutoff vs the box** | `r_c < L/2` | | ❌ hard — a minimum-image violation |
| **finite size (the persistence length)** | `ℓ_p ≤ L/4` | | ⚠️ a warning — an active artefact |
| **finite size (the correlation length)** | `ξ ≤ L/4` (measured post hoc) | | ⚠️ a warning — a finite-size artefact |
| **screening** | `λ_D ≤ L/4` | | ⚠️ a warning — the screening is incomplete |
| **observation sufficient** | `T_obs ≥ 10² · max(τ_B, τ_r)` | | ⚠️ a warning — insufficient statistics |
| **insufficient margin** | within 1/5 of the limit on any hard check | | ⚠️ a warning — no room to raise the parameter |
| **dense/dilute** | report `a_mean / d` | | ℹ️ information — the regime verdict |

> The Brownian integrator has an O(δt) error (§11 trap 2), so passing the ceilings above does not guarantee **accuracy**.
> A paper-grade result needs a separate convergence confirmation with `dt` halved (§17 the convergence study).

A list of `SeparationCheck(name, ratio, threshold, verdict, message)` is stored in the `ScaleLedger`, and
if even one is ❌ **the validator refuses to run**.

### 6.4b `dt` is inverted from the bias ⭐️ (v0.4, confirmed by measurement)

"Small enough" is vague. **In a linear system the systematic bias can be computed exactly.**
The Euler–Maruyama discrete steady-state variance of a harmonic trap:

```
⟨x²⟩_discrete = (kT/k) / (1 − h/2),   h ≡ dt/τ    ⟹   the relative bias ≈ h/2
```

**Half of `dt/τ` is the systematic bias.** That HOOMD's Brownian follows this law was measured
(`scratch/dt_convergence.py`, N=2000):

| `dt/τ` | The measured bias | Theory | Difference |
|---|---|---|---|
| 0.10 | 5.262% | 5.263% | −0.002% |
| 0.05 | 2.539% | 2.564% | −0.025% |
| 0.02 | 1.079% | 1.010% | +0.069% |
| 0.01 | 0.569% | 0.503% | +0.067% |

| Target accuracy | `dt/τ_fast` | Relative cost |
|---|---|---|
| 1% | 2e-2 | 5× |
| **0.5%** | **1e-2** ← the §6.4 hard gate | 10× |
| 0.1% | 2e-3 | 50× |
| 0.025% | 5e-4 | 200× |

→ **The §6.4 hard gate of `1e-2` means "a 0.5% bias".**
A nonlinear system has no closed form so it needs a convergence study, but the O(dt) scaling holds, so
a two-point Richardson extrapolation works well. Details: skill `bd-physics` §1.2.

### 6.5 `DimensionlessReport` — the subject of human check #3

```
System: 2D active colloid  (label: abp_silica_1p2um)
Reference scales: length=d, time=τ_B, energy=kT   [strategy: thermal]
  rationale: a thermal-fluctuation-dominated region (Pe~40), and the literature convention is the τ_B reference
══════════════════════════════════════════════════════════════════════

INPUT (dimensional, SI)
  d   = 1.20 µm      [sketch:note1.png#ann3, tier 1]
  T   = 298 K        [user, tier 0]
  η   = 1.00 mPa·s   [handbook:water@25C, tier 0]
  ρ_p = 2000 kg/m³   [kb:silica, tier 0]
  φ   = 0.600        [user, tier 0]
  N   = 4000         [user, tier 0]
  v₀  = 11.6 µm/s    [kb:10.1103/PhysRevLett.110.238301#fig2, tier 2]  ⚠ unverified

DERIVED (dimensional)
  γ   = 1.131e-8 kg/s        D_t = 0.364 µm²/s       D_r = 0.758 1/s
  m   = 1.81e-15 kg          L   = 86.8 µm           ℓ_p = 15.4 µm

SCALE LEDGER
  lengths   d=1.20µm  <  a_mean=1.37µm  <  r_c=1.35µm  <  ℓ_p=15.4µm  <  L=86.8µm
  times     τ_p=1.6e-7s  <  dt=1.98e-4s  <  τ_v=0.103s  <  τ_r=1.32s
            <  τ_B=3.96s  <  T_obs=3960s
  energies  kT=4.11e-21 J  =  ε=4.11e-21 J  <  f_a·d=1.58e-19 J

DIMENSIONLESS GROUPS
  φ      = 0.600      the occupied-volume ratio
  Pe     = 38.4       τ_B/τ_v      advection vs diffusion
  D_r*   = 3.00       τ_B/τ_r      rotational vs translational  (the Stokes prediction 3.00 ✓)
  ℓ_p/d  = 12.8       Pe/D_r*      the persistence length = 12.8 particles
  ε*     = 1.00       ε/kT         WCA
  dt*    = 5.0e-5     dt/τ_B
  L/d    = 72.4       the box = 72 particles
  T*     = 1000       the observation window = 1000 τ_B
  St     = 4.0e-8     τ_p/τ_B
  Re     = 1.4e-5     the fluid inertia

SCALE SEPARATION CHECKS                        value      limit    margin
  ✓ inertia neglected  τ_p/dt        =  8.1e-4   ≤ 1e-2     12.4×
  ✓ advection resolved dt/τ_v        =  1.9e-3   ≤ 1e-2      5.2×   ← the tightest
  ✓ rotation resolved  dt/τ_r        =  1.5e-4   ≤ 1e-2     66×
  ✓ interaction resolved dt/τ_int    =  5.0e-5   ≤ 1e-2    200×
  ✓ low Reynolds       Re            =  1.4e-5   ≤ 1e-3     71×
  ✓ cutoff             r_c/L         =  0.016    <  0.5      31×
  ⚠ finite size        ℓ_p/L         =  0.178    ≤ 0.25      1.4×   ← insufficient margin
  ✓ observation enough T_obs/τ_B     =  1000     ≥ 100       10×
  ℹ packing            a_mean/d      =  1.14              → a dense system (contact-dominated)

VERDICT: PASS (2 warnings)
  ⚠ ℓ_p/L = 0.178 — the persistence length is 1/5.6 of the box. Only a 1.4× margin to the limit (1/4).
     → N raised to 16000 gives L=174µm, ℓ_p/L=0.089 (a 2.8× margin). The wall-clock time is 4×.
  ⚠ v₀ is tier 2 (an unverified literature value) — human confirmation needed.
  ℹ Pe: if it is to be raised, advection resolution (a 5.2× margin) binds first. At Pe=200, make dt* 1e-5.

RESOURCE
  production 2.0e7 steps = 1000 τ_B = 66 min (physical) ≈ 47 min (wall clock)
  raw data: A 0.5GB + B 0.2GB + D 0.05GB = 0.75 GB
```

This report is the subject of **human check #3**. Rather than showing only the dimensionless numbers, it presents the
scale ledger and the separation checks alongside so a human can judge "does this make physical sense".

### 6.6 ⑤ The inverse transform and reverse construction

**The inverse transform (results → physical units)** — always performed:
```python
D_eff = D_eff_star * (sigma**2 / tau_B)      # 1.83 → 0.666 µm²/s
P     = P_star * (kT / sigma**dim)           # dimensionless pressure → Pa
t     = t_star * tau_B                       # steps → seconds
```
`observables.parquet` always stores `y` (dimensionless) and `y_physical`+`y_unit_si` as a pair.

**Reverse construction (dimensionless → a physical system)** — the **only** route when you want to start from dimensionless values:
```python
# "Pe=40, φ=0.6 2D ABP, run it for me"  ← there is no physical system
system = PhysicalSystem.from_dimensionless(
    groups={"Pe": 40, "phi": 0.6, "D_r_star": 3.0},
    anchors={                                    # ⭐️ the anchors must be specified
        "d": 1.0 * ureg.micrometer,              #    (otherwise the physical system is undetermined)
        "T": 298 * ureg.kelvin,
        "eta": 1.0 * ureg.mPa * ureg.s,
    },
    note="literature-convention coordinates. The anchors are referenced to a typical silica colloid.",
)
```
Without the anchors the physical system is undetermined, so it is **an error**. The anchor defaults are provided by
`materials.py`, but the fact that they are "arbitrarily chosen anchors" is marked with `Provenanced.confidence = 3`
and stated in the report.
This way, even while working dimensionlessly, the scale-separation checks (§6.4) still apply.

---

## 7. The knowledge base (L1 / L10) — in detail

### 7.0 The current status — starting from `record.json` (2026-08-03) ⭐️

**The SQLite KB has not been built yet.** With 1 run and 0 papers it is premature abstraction.
Instead, **only the data format is fixed first** and things are flowed through it:

```
cases/*.py            → runs/<id>/metrics.json   (the machine-readable result)
tools/postmortem.py   → runs/<id>/record.json    (automatic diagnosis + a tier3 KB entry)
tools/kb.py           → list | query | lessons   (glob + filters)
```

This way, when Phase 5 arrives, **there is already a history to feed back.**
Once the runs pass 100 or the literature comes in, it moves to §7.1's SQLite then (with the format unchanged).

`record.json` is a reduced version of §5.2's `KnowledgeEntry` — `system_tags`,
`dimensionless`, `observables` (measured vs analytic), `outcome`, `failure_modes`,
`not_verified`, `lessons` (tier 3).

**The `not_verified` field matters.** It leaves "what was not directly confirmed" explicit
(for example `dt_convergence_direct`). A device to keep a success record from turning into overconfidence.

### 7.1 The store choice (in Phase 5)
It starts with **SQLite + FTS5**. The reasons:
- a local single machine with entries in the hundreds to thousands → a vector DB is overkill
- dimensionless-coordinate search (`Pe BETWEEN 30 AND 50 AND phi BETWEEN 0.55 AND 0.65`) is the crux, and for that a **structured SQL query** is more accurate than an embedding
- full-text search alongside, via FTS5
- once the entries pass a few thousand, add an embedding column then (the schema is opened in advance)

### 7.2 The three search modes
| Mode | Query | Use |
|---|---|---|
| **coordinate search** | a dimensionless-number neighbourhood | "prior work near Pe≈40, φ≈0.6" |
| **tag search** | matching `system_tags` | "2D ABP WCA" |
| **full-text search** | FTS5 | "MIPS onset" |

The three results are merged, sorted by confidence tier, and injected into the LLM context.

### 7.3 The paper distillation pipeline
```
PDF → (a) text/table extraction  → LLM distillation → a KnowledgeEntry draft
   → (b) figure extraction       → vision  ↗       ↓
                                     [human verification] → confidence 2 → 1
```
- The distillation prompt is **version controlled** (`distill_prompt_v3`). Change the prompt and re-distillation is possible.
- The LLM is instructed strongly to **"use null if it is not in the paper"**, and made to fill in the `locator` (the figure/table number) without fail.
- The verification UI: put the original snippet alongside each entry and ✓/✗/amend.

### 7.4 Feeding back our own runs
`PostMortem.lessons` is converted into `Claim`s and goes in as tier 3 entries.
```
"dt*=1e-4 at Pe>80 diverged (3/3 runs). Safe: dt* ≤ 2e-5."
   kind = "pitfall",  coords = {"Pe": 80, "dt_star": 1e-4}
```
→ promoted to an **automatic warning rule** at the next verification (with a promotion to a hard error proposed if the same failure occurs 3 or more times).

---

## 8. Intake (L0) — interpreting sketches, notes and pictures

### 8.1 The input forms
| Form | Handling |
|---|---|
| a hand-drawn sketch (a photo/scan) | Claude vision → `Observation` |
| hand-note text | vision (handwriting) or plain text |
| a whiteboard photo | vision, several images supported |
| a capture of a paper figure | vision + a contextual query |
| a paper PDF | text and figures separated, then distilled (7.3) |
| natural language | as it is |

The technique: an `{"type": "image", "source": {"type": "base64", ...}}` block in `client.messages.create()`.
Several images go as several image blocks in one message. A PDF goes as a `{"type": "document"}` block (a base64 PDF) or via the Files API.

### 8.2 The interpretation protocol (suppressing hallucination)
What the LLM is required to do:
1. **First copy it down verbatim** (`raw_transcription`) — transcription before interpretation
2. Then structure it (`entities`, `stated_quantities`)
3. **Fill in `ambiguities` and `unread_regions` without fail** — an empty list is suspicious
4. Never invent a value **absent** from the sketch → leave a gap as `null` and fill it in at L1 (the KB) with grounds

### 8.3 Human check #1 — approval per field
A free-form summary is not shown with "is this right?". It is presented as **a per-item form**:
```
[1] system_guess: "2D active colloidal monolayer"        [✓ approve] [✗] [amend: ___]
[2] particle diameter ≈ 1 µm    (source: the figure's top left, "d~1um")  [✓] [✗] [amend: ___]
[3] Pe range 20–100             (source: the arrow label on the figure's right)  [✓] [✗] [amend: ___]
[!] ambiguity: unclear whether the arrow is the self-propulsion direction or the shear direction    [interpretation: ___]
[!] unread: the handwriting at the bottom right is illegible                          [enter: ___]
```
The approval result is recorded in the approval ledger.

---

## 9. The raw-data storage strategy (L7) ⭐️

"Store the particle motion and the forces and take them out later" — do that as it stands and the disk blows up.
At N=10⁴ and 10⁸ steps, storing the positions and forces every step = **tens of TB**.
→ Solved with a **hierarchical storage policy**.

### 9.1 The storage tiers

| Tier | Content | Frequency | Size (referenced to N=10⁴, 10⁷ steps) | Default | **Implementation** |
|---|---|---|---|---|---|
| **A** | positions + orientations (all particles) | every `10⁴` steps | ~0.5 GB | always ON | ✅ **used in 1-A** |
| **B** | positions + orientations + **velocities/forces/torques** | every `10⁵` steps | ~0.2 GB | ON by default | only confirmed working |
| **C** | a high-frequency burst (a short stretch) | `10` steps × 1000 frames | ~2.5 GB | on request | only confirmed working |
| **D** | tracked particles (a subset, high frequency) | 100 of them, every `10` steps | ~0.05 GB | ON by default | only confirmed working |
| **L** | the global scalar log | every `10³` steps | < 10 MB | always ON | only confirmed working |

> All 5 tiers were **confirmed working** (Phase 0, `scratch/smoke.py` 15/15). 1-A had no interparticle
> interactions so only Tier A was needed. Tier B is turned on for the first case where the forces have
> meaning (1-B soft-r3).

**Tier C is the key idea**: `hoomd.write.Burst` holds the most recent N frames in memory in a sliding window and writes to disk only when `.dump()` is called. → It can store "only around the moment of an event of interest, at high resolution".
For example the moment a cluster forms, or the moment the energy changes abruptly.

**Tier D**: track a few particles at high frequency all the way through → sufficient for the MSD, the velocity autocorrelation and single-particle trajectory statistics. There is almost never a need to store every particle at high frequency.

### 9.2 The HOOMD implementation mapping
| Tier | Implementation |
|---|---|
| A | `hoomd.write.GSD(trigger=Periodic(1e4), dynamic=['property'])` |
| B | `GSD(trigger=Periodic(1e5), dynamic=['property','momentum'])` + per-particle forces added to the `logger` → the GSD `log/` namespace |
| C | `hoomd.write.Burst(max_burst_size=1000, trigger=Periodic(10))` + a custom Action calling `.dump()` when the condition is met |
| D | `GSD(filter=hoomd.filter.Tags([...]), trigger=Periodic(10))` |
| L | `hoomd.logging.Logger` + `hoomd.write.HDF5Log` |

Per-particle force logging:
```python
logger = hoomd.logging.Logger(categories=['particle'])
logger.add(lj, quantities=['forces', 'energies'])
logger.add(active, quantities=['forces', 'torques'])
gsd_b = hoomd.write.GSD(filename='traj_forces.gsd', trigger=Periodic(int(1e5)),
                        mode='xb', logger=logger, dynamic=['property','momentum'])
```

### 9.3 The policy is stated in the spec
```python
class RawDataPolicy(BaseModel):
    tier_a_every: int = 10_000
    tier_b_every: int | None = 100_000       # None turns it off
    tier_c: BurstPolicy | None = None        # conditional high frequency
    tier_d_n_tracers: int = 100
    tier_d_every: int = 10
    log_every: int = 1_000
    estimated_bytes: int                     # computed by the validator, shown to the human
```
The validator computes the expected volume and **requires human approval above 5 GB**.

### 9.4 On-demand re-analysis
- **A frame index**: `index.parquet` generated at the end of a run (frame number ↔ step ↔ file offset ↔ time)
- **A lazy loader**: `loader.frames(run_id, tier='A', steps=slice(1e6, 2e6))` → reads only the frames needed
- **An observable cache**: a computed observable is stored in `observables.parquet` and the same request is not recomputed
- **An agent tool**: `query_raw_data(run_id, what, when, who)` — for example "the force distribution of the particles inside the cluster over the 3e6~4e6 step stretch"

### 9.5 The retention policy
- a completed run: delete Tiers C/D after 90 days (keeping only A/B/L)
- a failed run: delete the trajectory after 30 days (the `PostMortem` and the logs are kept permanently)
- confirm with `bdbot gc --dry-run` and then clean up

---

## 10. The post-mortem & learning loop (L10) ⭐️

> **Implementation status (2026-08-03)**: [`tools/postmortem.py`](../../tools/postmortem.py) implements the §10.1
> taxonomy and the §10.2 automatic diagnosis and emits `record.json` (a tier 3 KB entry).
> Queries go through [`tools/kb.py`](../../tools/kb.py). **What does not exist yet**: the LLM narrative (§10.3) and
> rule promotion (§10.4). With 1 run there is no sample from which to discuss promotion.
>
> An actual output example is `runs/trap-2d-5um__70b9394e7310/record.json`.

### 10.1 The failure taxonomy
Left as free text it does not accumulate. It is mapped onto **a predefined classification**.

| Code | Meaning | The automatic detection indicator |
|---|---|---|
| `NUM_DIVERGE` | numerical divergence | NaN/Inf, PE >100× the initial |
| `NUM_DRIFT` | an energy/temperature drift | a significant PE trend after steady state |
| `EQ_INSUFFICIENT` | insufficient equilibration | the first-half/second-half block averages disagree |
| `STAT_INSUFFICIENT` | insufficient statistics | an observable's error bar > the threshold |
| `FINITE_SIZE` | a finite-size effect | the correlation length > L/4 |
| `WRONG_REGIME` | the target phenomenon did not appear | the target observable is outside the expected range |
| `RESOURCE` | time/disk exceeded | the runner's record |
| `SPEC_ERROR` | the spec itself is wrong | a human verdict |
| `SUCCESS` | success | everything above passes |

### 10.2 The automatic diagnosis (without an LLM)
Always run at the end of a run:
- **the equilibrium verdict**: divide the trajectory into 5 blocks and test the stationarity of the PE and pressure block averages
- **the energy drift**: the significance of a linear-regression slope
- **finite size**: the g(r) correlation length vs L/4
- **statistical sufficiency**: produce error bars by the block-average method and check the relative error
- **goal attainment**: whether the `target_observables` are inside the expected range

### 10.3 The LLM narrative + extracting lessons
The automatic diagnosis result + the spec + the dimensionless coordinates are given to the LLM to produce:
1. a natural-language narrative (`narrative`)
2. **reusable lessons** (`lessons: list[Claim]`) — the dimensionless coordinates must be included
3. a proposal for the next attempt (an amended-parameter proposal)

### 10.4 Rule promotion
When the same failure repeats in the same coordinate region:
```
3 repeats → propose automatically generating a validator warning rule (human-approved)
5 repeats → propose promotion to a hard error
```
The rule is generated as code in `validate/learned_rules.py`, and a human reviews and merges it.
**The LLM does not write a verification rule directly as code.** It only fills the parameters into a template.

### 10.5 Record the success factors too
Recording only failures is biased. A successful run too:
```
"φ=0.6, Pe=40, dt*=5e-5, N=4000, prod=1000τ_B → MIPS observed, statistics sufficient, 47 min"
   kind = "parameter", tier = 3
```
→ proposed as **a verified starting point** on a similar request.

---

## 11. The HOOMD API mapping (documentation v7.1.1 / **the installed build 7.1.0** — the APIs surveyed are identical)

| SimSpec | HOOMD v7.1.1 |
|---|---|
| device | `hoomd.device.CPU()` |
| seed | `hoomd.Simulation(device=dev, seed=spec.run.seed)` |
| box + init | build a `gsd.hoomd.Frame` → `sim.create_state_from_snapshot(...)` |
| pair `wca` | `md.pair.LJ(nlist=cell, default_r_cut=2**(1/6)*σ, mode='shift')` |
| pair `lj` | `md.pair.LJ(nlist=cell, default_r_cut=r_cut, mode='shift')` |
| pair `yukawa` | `md.pair.Yukawa(nlist=cell, default_r_cut=r_cut)` |
| nlist | `md.nlist.Cell(buffer=0.4)` |
| the BD integration | `md.methods.Brownian(filter=All(), kT=kT, default_gamma=γ)` |
| the integrator | `md.Integrator(dt=dt, methods=[bd], forces=[...])` |
| the ABP force | `md.force.Active(filter=All())`; `active.active_force['A'] = (f_a, 0, 0)` |
| the ABP rotational diffusion | `active.create_diffusion_updater(trigger, rotational_diffusion=D_r)` |
| the trajectory | `hoomd.write.GSD(...)` |
| a high-frequency burst | `hoomd.write.Burst(...)` + `.dump()` |
| the log | `hoomd.logging.Logger` + `hoomd.write.HDF5Log` |
| thermodynamic quantities | `md.compute.ThermodynamicQuantities(filter=All())` |
| restart | `GSD(mode='wb', truncate=True, dynamic=[...])` |
| overlap removal | a short `md.minimize.FIRE(...)` → then swap in the main integrator |

### ⚠️ The trap list that must be observed
1. **There is no dedicated WCA class** → `md.pair.LJ` + `r_cut=2^(1/6)σ` + `mode='shift'`. (`ForceShiftedLJ` is not WCA)
2. **Brownian has an O(δt) error** — the documentation states it. It needs a far smaller dt than Langevin. The default is `dt* = 1e-4`, reduced further when the forces are strong.
3. **The ABP rotation goes in the updater, not the integrator**:
   ```python
   integrator.integrate_rotational_dof = False   # ← mandatory
   sim.operations.updaters.append(
       active.create_diffusion_updater(trigger, rotational_diffusion=D_r))
   ```
   Leave it `True` and inertial rotation mixes in and it is no longer an ABP.
4. **2D**: `Lz=0`, `dimensions=2`. Initialize the orientation quaternion so that it rotates about the z axis only.
5. **BD is overdamped** — the velocity has no physical meaning. `thermalize_particle_momenta` is unnecessary. A velocity-based MSD is forbidden.
6. **`r_cut < L/2`** — otherwise a minimum-image violation.
7. ⭐️ **An external force + a periodic boundary: leave the minimum image out and it is quietly wrong.** A trap confirmed by measurement.
   In a trap towards a fixed anchor, using `d = pos - anchor` as it stands means that the moment the particle wraps across the box
   the distance jumps by L and it receives an enormous restoring force **in the wrong direction**. **It does not blow up, it is quietly wrong** —
   with a strong trap the correct value comes out and the weaker it gets the larger the error (k=10 +0.2% ✓ / k=2 **+1856%** ✗).
   ```python
   d = pos - anchors[tags]
   d -= L * np.round(d / L)      # ← this one line. It applies to every external.* module
   ```
   → `PhysicsModule.periodic_safe()` is made a required declared field, and the validator checks it.
8. **A NaN in the z component of the minimum image in 2D** — leaving the box length of a non-periodic axis as `inf` gives
   `inf * round(0/inf) = nan`. A NaN enters the force array and the runtime guard raises a false positive. Mask to the periodic axes only.
9. **`write.Burst` needs `write_at_start=True` for a new file** — without it,
   `RuntimeError: Must set write_at_start to write to a new file.`

> Traps 7~9 were discovered by measurement. The full survey and verification results: **[`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md)**

---

## 12. The verification layers (L4)

### 12.0 Structural checks (the invariants) ⭐️
- if `SimSpec.derived_from` is empty, **reject** — a dimensionless spec with no physical system cannot run (principle 3)
- if there is no `ScaleLedger`, or the reference scales' `rationale` is empty, **reject**
- confirm that every field of the `PhysicalSystem` passes the `pint` dimensional check (a unit mismatch = immediate rejection)

### 12.1 The scale-separation hard gate ⭐️
The ❌ items of §6.4 are the validator's hard errors. **All unified on the `10⁻²` criterion**:
`τ_p/τ_dyn ≤ 10⁻²` (the model) · `dt/τ_v ≤ 10⁻²` · `dt/τ_r ≤ 10⁻²` · `dt/τ_int ≤ 10⁻²` · `Re ≤ 10⁻³` · `r_c < L/2`

The error message always includes **the margin to the limit and a concrete amendment**:
```
❌ advection resolution violated: dt/τ_v = 3.2e-2 (the limit is 1e-2, exceeded 3.2×)
   cause: at Pe=200 τ_v=0.0198s while dt=6.3e-4s
   fix: lower dt* 5e-5 → 1.5e-5 (the wall-clock time increases 3.3×)
```

### 12.1b The margin warning
A warning if any hard check is within 1/5 of the limit — since raising the parameter even a little would violate it,
if a sweep is planned it is checked in advance across the whole sweep range.

### 12.2 Spec range checks
the `φ` ceiling (2D 0.9, 3D 0.74) · `100 ≤ N ≤ 10⁵` · no duplicate specification of Pe/force · the write-frequency floor · the disk ceiling

### 12.2b Physical consistency warnings
A warning if `D_r*` deviates far from the Stokes prediction (a sphere = 3) — if deliberate, a reason is required.
If `a_mean/d` is near 1 it is classified as a dense system, if ≫1 as a dilute one, and tagged for the later interpretation.

### 12.3 KB-based learned warnings ⭐️
```
⚠ a past warning: in the Pe=90, dt*=1e-4 region it diverged 3 times (run_ids: a3f9…, b71c…, e024…)
             recommended: dt* ≤ 2e-5
⚠ a literature comparison: Redner et al. 2013 used N≥10⁴ in this region (currently N=2000)
             → a finite-size effect is possible
```

### 12.4 The provenance-confidence check ⭐️
A spec composed only of tier 2 or lower values → human approval enforced.

### 12.5 The runtime guard
Check for NaN/Inf and a PE explosion every `10⁴` steps → halt + `status: "diverged"` + preserve the last snapshot.

---

## 13. Analysis & observables (L8)

### 13.1 The observable list
**structure**: `rdf` g(r) · `sq` S(q) · `psi6` the hexagonal order (2D) · `voronoi_density`
**dynamics**: `msd` · `diffusion_coefficient` · `fskt` · `vacf` (using Tier D)
**active**: `cluster_size` · `local_density_hist` (→ the MIPS bimodality) · `polar_order` · `giant_number_fluctuations`
**force statistics** ⭐️ (using Tier B/C): `force_distribution` · `stress_per_particle` · `virial_pressure_decomposition`
**thermodynamics**: `potential_energy` · `pressure` · `pressure_tensor`

### 13.1b The output formats (settled 2026-08-03)

| File | Content | Status |
|---|---|---|
| `report.txt` | the DimensionlessReport (the scale table + the separation checks + the margins) | ✅ |
| `result.txt` | + the analytic comparison after the inverse transform to physical units | ✅ |
| `metrics.json` | **the machine-readable** result — the coordinates, the observables, the checks | ✅ |
| `record.json` | the automatic diagnosis + `outcome` + `lessons` + `not_verified` | ✅ |
| `observables.npz` | the raw P(x), C(t) and PSD | ✅ |
| `observables.png` | a 4-panel plot (measured vs analytic) | ✅ |
| `traj_A.gsd` | the raw trajectory (Tier A) | ✅ |
| `observables.parquet` | the long-format table (§13.2) | ❌ unimplemented — being substituted by npz |

### 13.2 The output schema (the target)
```
run_id | observable | x | y | y_err | x_unit_reduced | y_unit_reduced
       | x_physical | y_physical | x_unit_si | y_unit_si | metadata(json)
```
⭐️ **The dimensionless value and the physical-unit value stored simultaneously** (principle 3).

### 13.3 The standard plots
1. energy/pressure vs time (**always**, the grounds for the equilibrium verdict)
2. g(r) / MSD (log-log, with a slope-1 guide) / S(q)
3. ABP: the local-density histogram, the cluster distribution, the final snapshot
4. ⭐️ **the literature-comparison plot**: overlay the KB's literature values as points at the same dimensionless coordinates

### 13.3b Checking a result is **graphs + animations by default** ⭐️⭐️ (2026-08-05)

> A result is not reported in words and tables alone. Since the usual way a human checks a result is
> **pictures and video**, producing both and attaching them to the report is the default behaviour.

**Graphs — "what is wrong" has to be visible.** Put the measurement and the prediction on the same axes and
overlay the check thresholds, the analytic solution and the literature values as lines. The model: [`scratch/viz_chain_bend.py`](../../verify/viz_chain_bend.py)
— its 6 panels correspond one each to `chain-bend`'s 5 spec errors, and in each panel the "wrong value" and the
"right value" are in the same picture (for example the ω-independent straight line the SNR check drew vs the curve that actually falls).

**Animations — how the system actually moves.** In particular, putting `kT=0` (deterministic, the mode shape) and
`kT>0` (thermal fluctuation) **side by side** shows an SNR problem faster than a number. In `chain-bend`,
two panels one above the other revealed at a glance that "the driven response is buried inside the ±ℓ_k band".
An animation may be made **cheaply with a large `dt`** (about 10% of the stability limit). But **it must be stated
that it is not the production measurement**, and it has to be confirmed that the static equilibrium quantities (ℓ_k and so on) are correct at that `dt` too.

**⚠️ Write graph labels in English, not Korean** (measured). matplotlib's default `DejaVu Sans` has no
Hangul glyphs so the labels all become `□`, and switching to a font that has Hangul loses the symbols instead —
`AppleGothic` · `Apple SD Gothic Neo` · `NanumGothic`, all three, **lack** `−` (U+2212) and
`ŷ` (U+0177). The only one with both Hangul and the symbols is `Arial Unicode MS`, so
switching fonts cannot fix it completely (`axes.unicode_minus = False` also only fixes the tick labels and
**the `−` in a string literal still breaks**). → **Do not try to match the font; write the axes, legends, titles
and annotations in English from the start.** After generating, be sure to confirm that the `missing from font`
warnings are **0** — ignore a warning and you end up presenting as a result a picture a human cannot read.
The rule is kept even when the font is not a problem — this field's physics and statistics terms (`de-correlation
time`, `shear thinning`, `yield force` and so on) are often awkward or unestablished in Korean,
so the English conveys the meaning more precisely and concisely.

ffmpeg is **not** in this environment → a GIF is written with `PillowWriter` (`FFMpegWriter.isAvailable()` = False).

---

## 14. Execution / job management (L6)

- **A process pool of 8** (10 cores − 2). A worker is a fully independent process.
- **A SQLite job queue** (`runs.db`): run_id, status, pid, progress, error
- **Checkpoints**: `restart.gsd` overwritten with truncate, `bdbot resume <run_id>`
- **The state machine**: `queued → running → {completed | failed | diverged | interrupted}`
- **Sweeps**: expand `SweepSpec.axes` → automatic deduplication by the run_id hash

---

## 15. The Claude Code integration layer ⭐️ (changed wholesale in v0.3)

**Claude Code is the agent runtime.** We do not build an agent; we build only the
**deterministic engine (`bdbot`)** that sits beneath it and the **integration surface (`.claude/`)**.

### 15.1 The division of roles — what disappears

| Plan element | v0.2 (our own implementation) | v0.3 (Claude Code) |
|---|---|---|
| the agent loop | `tool_runner()` + `@beta_tool` | **built in** |
| model calls, caching, streaming | the anthropic SDK directly | **built in** |
| reading a sketch (vision) | a base64 image block | **the Read tool — free** |
| reading a paper PDF | a document block / the Files API | **the Read tool (the `pages` argument) — free** |
| managing the conversation history | implemented by hand | **built in** |
| the 10 tools | `@beta_tool` functions | **the `bdbot` CLI (Bash)** |
| injecting domain knowledge | the system prompt + `cache_control` | **a Skill + `CLAUDE.md`** |
| the approval gates | a branch inside the tool | **Hooks + AskUserQuestion + permissions** |
| the workflow norms | prompt instructions | **Slash commands** |
| separating specialist roles | (none) | **Subagents (context isolation)** |
| polling a long-running job | implemented by hand | **Bash `run_in_background`** |

**Deleted**: most of `bdbot/agent/` · `bdbot chat` · the mandatory `anthropic` dependency
**The remaining use for anthropic**: **batch** distillation of papers only (§15.9)

### 15.2 The integration method: **CLI + Skill** (not MCP)

| Method | Advantages | Disadvantages | Verdict |
|---|---|---|---|
| **CLI + Skill** | no server needed · being file-based it can be inspected, diffed and version-controlled · permission control via a `Bash(bdbot:*)` allow-list · a human uses it identically | the schema is not exposed directly to the model (→ the Skill compensates) | ✅ **adopted** |
| an MCP server | the type schema exposed · structured returns | a server process to keep alive · overkill for a local Python package · hard to debug | thinly on top of the CLI later, if needed |
| files only + prompts | minimal | verification and gates cannot be enforced | ❌ |

**The core grounds**: the outputs of this pipeline are all files (`system.yaml` → `spec.json` → `traj.gsd` → `observables.parquet`). Claude Code's main business is file work, so
`Write(system.yaml)` → `Bash(bdbot nondim system.yaml)` → `Read(report.txt)` is the most natural.
On top of that, every intermediate output remains as a file a human can open — which fits principle 2 (provenance tracking) well.

### 15.3 The `.claude/` structure

```
.claude/
├── settings.json                    # the permission allow-list + hook registration
├── skills/
│   ├── bd-physics/SKILL.md          # the unit system, the scale ledger, the non-dimensionalization convention (§2, §6)
│   ├── bd-hoomd/SKILL.md            # HOOMD v7.1.1 API mapping + the 6 traps (§11)
│   ├── bd-intake/SKILL.md           # the sketch interpretation protocol (§8.2)
│   └── bd-distill/SKILL.md          # the paper distillation protocol (§7.3)
├── agents/
│   ├── bd-distiller.md              # 1 paper → a KnowledgeEntry
│   ├── bd-analyst.md                # a run result → a PostMortem
│   └── bd-reviewer.md               # an adversarial physical review of the spec
├── commands/
│   ├── bd-intake.md                 # /bd-intake <folder>
│   ├── bd-spec.md                   # /bd-spec <observation>
│   ├── bd-run.md · bd-sweep.md
│   ├── bd-analyze.md · bd-postmortem.md
│   └── bd-distill.md
└── hooks/
    ├── guard_invariant.py           # ⭐️ derived_from absent, block the run
    ├── guard_separation.py          # ⭐️ block a failing separation check
    ├── guard_cost.py                # require confirmation when the cost is exceeded
    └── log_approval.py              # record in the approval ledger

CLAUDE.md                            # the always-loaded invariants (kept short)
```

### 15.4 `CLAUDE.md` — the always-loaded invariants only

A Skill loads when needed, but `CLAUDE.md` is **always** in the context. It is kept short and holds
**absolute rules** only. The details are deferred to the Skills.

```markdown
# bdbot — Brownian dynamics simulation pipeline

## Absolute rules
1. Dimensions come first. Every system is fixed as an SI-unit PhysicalSystem first and
   non-dimensionalized through the scale ledger. There is no route that starts from dimensionless values.
2. SimSpec is never written by hand. Only what `bdbot nondim` generated gets run.
3. Attach a provenance (tier 0–3) to every parameter. If it is not in the KB, say it is not.
   Mark an estimate as an estimate. Never invent one.
4. HOOMD scripts are never written directly. Use only the `bdbot` CLI.
5. When reading a sketch, first transcribe it verbatim, and state without fail
   what could not be read and what is ambiguous.

## Detailed knowledge
- non-dimensionalization and the scale ledger → skill `bd-physics`
- HOOMD API·traps → skill `bd-hoomd`
- the workflow → `/bd-intake`, `/bd-spec`, `/bd-run`, `/bd-postmortem`
```

### 15.5 Hooks — **the harness** enforces it, not the prompt ⭐️⭐️

The biggest gain of this change. "Ask" for a hard invariant like principle 3 in a prompt and
the model can break it. **A PreToolUse hook is executed by the harness**, so the model cannot route around it.

| Hook | Event | Action |
|---|---|---|
| `guard_invariant` | PreToolUse `Bash(bdbot run*)` | **reject** if the spec has no `derived_from`/`scale_ledger` |
| `guard_separation` | PreToolUse `Bash(bdbot run*)` | **reject** if the §6.4 hard gate fails, and return an amendment |
| `guard_cost` | PreToolUse `Bash(bdbot run*\|sweep*)` | require user confirmation when the estimated cost is exceeded |
| `log_approval` | PostToolUse | record approvals/rejections/amendments in the ledger |
| `capture_postmortem` | Stop | notify if a completed run has no post-mortem run on it |

When a hook rejects, it **returns a concrete amendment on stderr** so Claude Code can fix it itself:
```
❌ blocked: spec has no `derived_from`.
   Create a PhysicalSystem first: bdbot init-system --template abp_2d > system.yaml
```

> The hook configuration goes in `.claude/settings.json`. That file is safest written with the `update-config` skill.

### 15.6 Subagents — the purpose is context isolation

| Subagent | Role | Why isolated |
|---|---|---|
| `bd-distiller` | 1 paper → a `KnowledgeEntry` | **the paper's full text does not contaminate the main conversation.** 5 papers can run as 5 at once |
| `bd-analyst` | a run result → a `PostMortem` | the bulk numerical output does not enter the main context |
| `bd-reviewer` | an adversarial review of the spec | **an independent perspective** — it looks only at the physical validity without knowing the context that made the spec |

`bd-reviewer` automates the "agent evaluation" of §17. Running "try to refute this spec's physical problems"
in an independent context before submitting the spec avoids the bias of reviewing your own work.

### 15.7 Slash commands — fixing the workflow

| Command | What it does |
|---|---|
| `/bd-intake <folder>` | read every sketch/note in the folder → an `Observation` YAML → human confirmation |
| `/bd-spec <observation>` | search the KB → propose a `PhysicalSystem` (with provenance attached) → human confirmation → `bdbot nondim` → the report |
| `/bd-run <spec>` | verify → estimate the cost → approve → run in the background |
| `/bd-sweep <sweep>` | expand the sweep → separation checks across the whole range → approve → enqueue |
| `/bd-analyze <run_id>` | compute the observables → plot → compare against the literature |
| `/bd-postmortem <run_id>` | the `bd-analyst` subagent → feed back into the KB |
| `/bd-distill <pdf...>` | the `bd-distiller` subagents in parallel → the verification queue |

Fix the workflow **as a command** rather than as a prompt norm and the order is observed.

### 15.8 Remapping the approval gates

| Gate | The v0.2 implementation | The v0.3 implementation |
|---|---|---|
| #1 Observation | a form inside the tool | AskUserQuestion + a YAML file review |
| #2 PhysicalSystem | a per-field form | a YAML file review (amendments confirmed by diff) + AskUserQuestion |
| #3 SimSpec/non-dimensionalization | checking the report | Read the report file + AskUserQuestion |
| #4 cost | a branch inside the tool | **the `guard_cost` hook** (harness-enforced) |
| #5 promoting a knowledge tier | a tool | AskUserQuestion |
| #6 promoting a rule | a tool | a human reviews it like a PR (a file diff) |

In the early stages, using **plan mode** so that a human approves the whole plan before execution is also valid.

### 15.9 The only place the `anthropic` SDK remains — batch distillation

Distilling 50 papers conversationally is slow and expensive in tokens. Bulk distillation goes in a separate script:
- **the Message Batches API** (`client.messages.batches.create`) — a **50% discount** against the standard price, up to 100 thousand requests
- receive Pydantic-validated results with `client.messages.parse(output_format=KnowledgeEntry)`
- the model `claude-opus-5`

It is kept as an **optional dependency** (`bdbot[distill]`) in `pyproject.toml`. Small-volume distillation is handled
in Claude Code with `/bd-distill`, so it is not normally needed.

### 15.10 The effect of this change

| Item | The change |
|---|---|
| the development schedule | Phase 7+8, 4.5 days → **1.5 days** (total 18~20 days → **15~17 days**) |
| the amount of code | ~800 lines of `bdbot/agent/` deleted, ~300 lines of `.claude/` added |
| the enforceability of the invariants | asking in a prompt → **harness-enforced** (stronger) |
| vision/PDF | implemented by hand → free |
| human intervention | a separate UI needed → the Claude Code conversation is the UI |
| the risk | automation does not work outside a Claude Code session → **all the more reason the CLI has to work on its own** |

---

## 16. The development roadmap

The **definition of done (DoD)** is stated for each Phase. It has to pass before the next.

### ✅ Phase 0 — a general environment & API demonstration **(complete 2026-08-03)**
```bash
conda env update -f environment.yml -n simulation_bot --prune
```
- [x] the environment built — extending the existing `simulation_bot` (hoomd **7.1.0** CPU, gsd 5.0.1, freud 3.5.0,
      numpy 2.5.1, + pyarrow/pydantic/pint/typer/rich added). Pinned with `environment.yml`.
- [x] **an exhaustive HOOMD capability survey** (`scratch/survey.py`) — 28 isotropic pairs, 17 anisotropic pairs,
      bonds/angles/dihedrals, frictional contact, many-body, meshes, long range, manifolds, rigid bodies, HPMC, MPCD confirmed
- [x] **15 APIs verified working** (`scratch/smoke.py`) — **15/15 PASS**.
      All 5 raw-data tiers (§9) confirmed working
- [x] **the harmonic-trap golden physics verification** (`scratch/golden_trap.py`) — `⟨x²⟩=kT/k`
      reproduced to within 0.6% over a 10× range of k, with a coefficient of variation of `⟨x²⟩·k` of **0.28%**
- [x] **3 new traps discovered** (§11 traps 7~9) — trap 7 in particular is of the quietly-wrong kind
- [x] the results documented → **[`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md)**
- **DoD met**: confirmed by measurement that all 5 `intake/` cases are implementable.
  No wrong entry in the §11 mapping table (3 traps added). The BD integrator's accuracy confirmed.

> **A side benefit**: it proved in advance that `trap-2d-5um` stands up as a golden physics test → Phase 2 shortened.

### Phase 1 — **case-driven incremental construction** ⭐️ (the method changed in v0.4)

> **Why it changed**: the previous plan was a 2-day framework-first design that built the whole `ScaleLedger`,
> 4 reference strategies, 11 dimensionless numbers and 8 separation checks **first**. But non-dimensionalization
> differs per system, and which scales are actually needed is only known by driving a case end to end. Generalize
> by guesswork and you build code nobody uses while the thing actually needed is missing.
>
> **The changed method**: drive one case to the end → find the commonality in the next case → abstract **then**.
> Each stage is half a day, and at the end of every stage **the physics is verified together with a human**.

#### The procedural rules
1. Choose one case (the easy ones first)
2. Read the sketch → a draft physical system → **confirm together**
3. Write that system's scale table by hand → **confirm together** ← the non-dimensionalization knowledge accumulates here
4. Non-dimensionalize → run → invert → compare against the golden/literature values → **review the result together**
5. Record what was learned (the scale table, the check thresholds, the traps) → on to the next case
6. **Only after two cases are finished** is the common structure extracted and abstracted

> Number 3 is the crux. "Which timescales this system has and what to take as the reference" is
> a physics judgment rather than code, and until knowledge accumulates a human has to look at it too.

#### ✅ Phase 1-A · `trap-2d-5um` end to end **(complete 2026-08-03)**
- [x] [`intake/trap-2d-5um/observation.yaml`](../../intake/trap-2d-5um/observation.yaml) — the sketch transcribed, 3 ambiguities resolved
- [x] [`intake/trap-2d-5um/system.yaml`](../../intake/trap-2d-5um/system.yaml) — the physical system (SI), a tier per field
- [x] [`cases/trap_2d_5um.py`](../../cases/trap_2d_5um.py) — the scale table → the separation checks → non-dimensionalization →
      execution → the inverse transform → the analytic comparison. **Specific to this case** (not a general framework)
- [x] the 4 separation checks classified as **model/integration/geometry/statistics**, with the margin displayed
- [x] `run_id` content addressing + preventing a re-run (§14)
- **DoD met**: a 5 µm, water, 300K YAML → results in physical units → **all 4 observables agreeing with their analytic solutions**
  (`⟨x²⟩` +0.02%, `σ` +0.01%, `τ` −0.25%, `f_c` +1.17%, `S(0)` −0.33%)
- **A discovery**: subtracting the sample mean in the displacement autocorrelation puts `τ` off by −7.75% → not subtracting gives −0.26%
  (recorded in bd-physics §5.1)
- **Output**: `runs/trap-2d-5um__70b9394e7310/` (the report, trajectory, observables and plots)

#### ✅ Phase 1-B · `soft-r3-2d-A-sweep` end to end **(complete 2026-08-04)**
- [x] `intake/soft-r3-2d-A-sweep/{observation,system}.yaml` — the density, the core, N and the anchors fixed by a human
- [x] [`cases/soft_r3_2d.py`](../../cases/soft_r3_2d.py) — `pair.Table` (r⁻³) + a WCA core
- [x] **7 runs**: an amplitude sweep of 4 (A=0.1·1·10·100) + 1 dilute-limit verification + 2 convergence checks
- [x] 5 verifications — a direct two-particle comparison **0.000%** · energy consistency **+0.00~0.67%** ·
      the hexagonal NN distance **+0.45%** · the dilute limit (including O(ρ)) RMS **2.43%** · dt halved **−0.004%**
- [x] the comparison verdict against 1-A complete → skill `bd-physics` **§6.3**
- **DoD met**: fixed in a table what appeared twice (the reference scales, the property relations, the use of `τ_p`,
  `dt=10⁻²·γ/the local stiffness`, the check classification, the report, run_id, the metrics schema) and what
  differs per case (the equilibrium indicator, the observables, the verification strategy, the governing timescale).
- **The decision**: the `ScaleLedger` abstraction gets introduced in 1-C. Except for §6.3's "do not make common".

**What was newly learned in 1-B**
- **3 traps** (skill `bd-hoomd` 10·11·12): the `pair.Table` grid is `endpoint=False`
  (otherwise the force is −1.65% off) · at `r<r_min` the force and energy are **0** (dangerous for a diverging potential) ·
  `seed` is truncated to 16 bits so **a different seed becomes the same run**
- **`A` alone is not the control parameter.** `Γ = A(d/a_mean)³` sets the structure.
  Without the density the amplitude sweep is not defined (the sketch had no density)
- **`Γ_max = N^{3/2}u_c/8`** — solving the absolute-kT cutoff criterion together with the minimum image puts
  an `A`-independent ceiling on the achievable coupling strength. The sketch's N=100 gives Γ≲1.25
- **The 2D `r⁻³` converges only in structure.** From `r_c` 5a→7a, ψ₆, NN and the coordination number are unchanged
  to within 0.15%, while the absolute `⟨U⟩/N` is +7.5% (because the tail falls off only as `1/r_c`)
- **Weak coupling is numerically more expensive.** Γ=0.03 is 12.8× Γ=30 (the WCA core sets dt)
- **Two assumptions of the post-mortem tool were specific to 1-A** — see the Phase 4 item below

#### ✅ Phase 1-C · abstraction **(complete 2026-08-04)**
Only **what actually appeared twice** in 1-A/1-B was made common → the 10 `bdbot/` modules (§4.1).
- [x] `ScaleLedger` + `thermal_reference` — only the scales actually used
- [x] `Check` — the 4-way model/integration/geometry/statistics classification + hard/soft + the margin
- [x] `dt = 10⁻²·γ/(the local stiffness)` — the trap's `γ/k` and the pair's `γ/U''(r_min)` were the same formula
- [x] `report.render` · `runid` · `metrics` · `stats` · `sim` · `materials` · `provenance`
- [x] the two cases refactored to run through the same route. Only the specific physics is left in the cases
      (`trap_2d_5um.py` 610→436 lines, `soft_r3_2d.py` 755→703 lines)
- **DoD met — confirmed by measurement**: [`scratch/verify_1c_equivalence.py`](../../verify/verify_1c_equivalence.py)
  - **all 8 `run_id`s preserved** (the spec hash byte-identical)
  - 1-A re-run: metrics **77 fields identical**, 0 violations. The 4 observables agreeing to 6 digits
  - 1-B A=100 re-run: metrics **124 fields identical**, 0 violations
  - only 3 permitted changes:
    ① the wall clock ② the block SEM changed by 1.2e-6 through a float32→float64 upcast (**the more accurate side**;
    the original trace was float32) ③ **unifying the sign convention of the energy-consistency `err_pct` — a bug fix**
    (the 1-B original had this one row as `(predicted−measured)/|measured|`, so its sign was opposite both to the
    other rows in the same file and to 1-A. The magnitude and the verdict are unchanged)
- A side verification: [`scratch/verify_bdbot.py`](../../verify/verify_bdbot.py) — whether the common modules reproduce
  the cases' measured values (5 properties · the dt convention · the hard/soft verdict · the autocorrelation correction · the trap defences)

- **`SimSpec` + `derived_from` + the Builder/Runner were not built.** What the two cases
  shared was **the convention** "a spec dict → sha256 → run_id", and `runid.py` is
  sufficient for that. The pydantic schema and the builder have not yet appeared twice even once.

#### Phase 1-D onwards · half a day per case added
The order `chain-bend` → `trap-drag` → `abp-rod`. Each case adds one module,
and each time the scale ledger grows naturally (the structure by which a §5.6 module contributes).

**Total estimate: 2~2.5 days** (previously 3.5). More importantly, **a working result comes out every half day**.

### Phase 2 — physics verification (the golden tests) (1.5 days) ⭐️
**Prove the physics is right before the LLM.**

> ⚠️ **Principle 8**: the below are all `implementation_check`s — the prediction is derived from the model
> implemented, so a mismatch is a bug. **A physical discovery does not come out here.**
> For a discovery, a case that tests an assumption the simulation does not impose has to be designed separately
> (for example `soft-r3`'s melting Γ, or whether `trap-drag`'s effective viscosity holds).
- [ ] `test_free_particle_msd`: `⟨Δr²⟩ = 2dD_t t`, recovering D to within 5%
- [ ] `test_diffusion_from_gamma`: change γ → confirm `D = kT/γ`
- [ ] `test_wca_rdf_dilute`: the dilute-limit g(r) → `exp(-βU(r))`
- [x] `test_abp_effective_diffusion`: **`D_eff = D_t + v₀²/(d·Λ)`**, `Λ` = the director decay rate.
      ~~`v₀²/(2(d−1)D_r)`~~ is only coincidentally right in 2D (a +29~31% error measured in 3D).
      HOOMD has `Λ = rotational_diffusion` (in both 2D and 3D). → `scratch/standalone_abp_diffusion.py`
- [x] `test_nondim_roundtrip`: the physical → dimensionless → physical identity — **a round-trip error ≤ 1.2e-16**
      (4 kinds: length², time, the diffusion coefficient, energy). `NondimSpec.physical(v, L=,T=,E=)` is called
      with the dimensional exponents and divided back by the reference scales to confirm the identity. A spec save→loaded
      gives the same value (`scratch/verify_nondim_guards.py` ⑥⑦). It can be used as it is when moved to pytest.
- [ ] `test_scale_invariance` ⭐️: **non-dimensionalize the same physical system with the `thermal` reference and with the
      `interaction` reference, run each → and the results inverted back to physical units agree**. The core test proving the non-dimensionalization is correct.
- [ ] `test_separation_gates`: whether each separation check catches an artificially violating case exactly (8 of them)
      — the L3 layer (ledger completeness, ratio consistency, invertibility) is complete with the 30 checks of `verify_nondim_guards.py`.
      The physical separation checks (the §6.4 gates) side is still outstanding.
- [ ] `test_reproducibility`: the same seed → an identical trajectory
- **DoD**: `pytest tests/test_golden_physics.py` all passing.

### Phase 3 — the raw-data tiers + on-demand analysis (2 days)
- [ ] implement the Tier A/B/C/D/L writers
- [ ] the frame index + the lazy loader + slicing
- [ ] the observable registry + the cache (13.1)
- [ ] the standard plot set
- [ ] CLI: `bdbot analyze <run_id> --obs msd,rdf`, `bdbot raw <run_id> --steps 1e6:2e6 --what forces`
- **DoD**: the force distribution of an arbitrary stretch can be extracted from a completed run within 30 seconds.

### 🟡 Phase 4 — automatic diagnosis + the post-mortem (**the automatic diagnosis and record.json are complete**, 0.5 days for the rest)
Complete: the failure taxonomy · the automatic diagnosis (equilibrium, drift, statistics, bias) · emitting `record.json` · KB queries
→ [`tools/postmortem.py`](../../tools/postmortem.py), [`tools/kb.py`](../../tools/kb.py)

**What was fixed in 1-B — 3 diagnostics were assumptions specific to 1-A** (all of them gave a *wrong verdict*):
1. **The equilibrium indicator was 'the displacement from the anchor'** → specific to a bound system (a trap). In a
   diffusive system the displacement grows without bound so it is always `EQ_INSUFFICIENT`. → the case declares the
   indicator with `metrics.equilibration`
2. **The drift t-test had no autocorrelation correction** → a run whose whole-range change was **−0.026%** of the mean
   was judged `NUM_DRIFT` at `t=−3.3`. `n_eff` is now obtained with the Sokal automatic window to inflate the SE, and
   **the significance and the magnitude are looked at together** (below 0.5% is not treated as a failure).
   → it was **the same kind of mistake** as 1-A's "error bars by block averaging" lesson
3. **Every separation-check failure was treated as hard** → bd-physics §4 defines statistics and finite size as
   ⚠ warnings. In 1-A everything passed so the distinction never surfaced
4. The 0.5% statistics target was hardcoded (referenced to 1-A's `⟨x²⟩`) → the case declares it with `numerics.stat_target_pct`
Remaining: the LLM narrative (§10.3) · rule promotion (§10.4, insufficient samples) · the approval ledger
- [ ] the failure taxonomy + the automatic diagnoser (10.2)
- [ ] generating the `PostMortem` (including the LLM narrative)
- [ ] the approval ledger
- **DoD**: a run with dt deliberately raised is automatically classified `NUM_DIVERGE`, and a run with a short equilibration is caught as `EQ_INSUFFICIENT`.

### Phase 5 — the knowledge base (2.5 days)
- [ ] the SQLite + FTS5 schema, the three search modes
- [ ] the paper distillation pipeline (text + figures)
- [ ] the verification UI
- [ ] the post-mortem → the KB feedback
- [ ] the KB-based verification warnings (12.3)
- [ ] the literature-comparison plot
- **DoD**: 3 ABP papers are distilled in, and a search for "Pe≈40, φ≈0.6" returns the relevant entries. Our runs' failures accumulate as tier 3 entries.

### Phase 6 — parallel sweeps (1.5 days)
- [ ] the SQLite job queue + a process pool of 8
- [ ] `SweepSpec` expansion + deduplication + the sweep plot (a phase-diagram heatmap)
- **DoD**: a 30-run Pe×φ sweep runs on 8 workers and a heatmap comes out.

### Phase 7 — the intake schema + the approval procedure (0.5 days) ⬇️ reduced
Since vision is provided free by Claude Code, what we build is **only the schema and the settling procedure**.
- [ ] the `Observation` schema + YAML serialization
- [ ] `bdbot intake init <folder>` — generate an empty `observation.yaml` template
- [ ] `bdbot intake check <folder>` — reject if `ambiguities`/`unread_regions` are unfilled
- [ ] skill `bd-intake` — the interpretation protocol (§8.2: transcription first, ambiguities stated, gaps as null)
- **DoD**: `/bd-intake intake/abp-rod-2d-run-flip` → the sketch is read and an `observation.yaml` is generated,
  the ambiguous items are explicitly enumerated, and a human confirms and amends them.

### Phase 8 — the Claude Code integration surface

Since Claude Code is the agent runtime (§15), what we build is **knowledge capture + enforcement devices**.
It splits into two lumps — the first has value right now, and the second is only meaningful once the CLI exists.

#### ✅ Phase 8-min · knowledge capture **(complete 2026-08-03)**
When a session ends, everything learned in it disappears. Preventing that was the purpose.
- [x] [`CLAUDE.md`](../../CLAUDE.md) — 6 invariants (54 lines). Number 6 is "verify a physical claim before stating it"
- [x] skill [`bd-hoomd`](../../.claude/skills/bd-hoomd/SKILL.md) (407 lines) — 9 traps (★ 5 of them of the
      quietly-wrong kind) + the hard constraint (translational anisotropy impossible, with measured grounds) + 14 verified snippets + an API reference
- [x] skill [`bd-physics`](../../.claude/skills/bd-physics/SKILL.md) (209 lines) — the 5-stage procedure ·
      the unit conventions · dimensionless number = a scale ratio · the separation checks (all 10⁻²) · **a scale table per case** (the trap verified)
- [x] `scratch/verify_skill_snippets.py` — **a regression test that extracts the documentation's code and runs it**
- **DoD met**: the trap code extracted from the skill document runs as it is with a +0.46% error and 0 NaNs.
  Syntax checks 14/14. Break the documentation and the test breaks.
- **Why it was done**: two things were wrong in this session — the missing minimum image (+1856% at k=2) and
  the misjudgment about rigid-body anisotropy. Not recorded, the same mistake gets made in the next session.

#### Phase 8-rest · the enforcement devices + the workflow (1 day) ← **after the CLI exists**
- [ ] the 4 hooks: `guard_invariant`, `guard_separation`, `guard_cost`, `log_approval`
- [ ] `.claude/settings.json` — the permission allow-list (`Bash(bdbot:*)`),
      writing to `specs/` and `runs/` forbidden (the second line of defence, appendix B.6)
- [ ] the 3 subagents: `bd-distiller`, `bd-analyst`, `bd-reviewer`
- [ ] the 7 slash commands
- [ ] skill `bd-intake`, `bd-distill`
- **DoD**: attempting `bdbot run` with a spec that has no `derived_from` gets **blocked by the hook**, which presents an amendment.
  From one sketch folder, `/bd-intake` → `/bd-spec` → `/bd-run` → `/bd-postmortem` flows all the way through.
- **Why later**: a hook is a device that intercepts `bdbot run`. With no `bdbot run` there is nothing to enforce.

### Phase 9 — (optional) extensions
automating the convergence study · automatic report generation · GPU support · active learning · **turning it into an MCP server (if needed)** ·
an anisotropic translational friction module (§20 question 10 option B)

**Total estimate: 13~15 days** (full time). **Phase 1 (driving a case) → 2 (the golden verification) is the critical path.**

### 16.5 The recommended execution order

```
✅ Phase 0       the general environment + the HOOMD capability survey + 15 APIs verified        complete
✅ Phase 8-min   CLAUDE.md + bd-hoomd/bd-physics skill              complete
✅ Phase 1-A     trap-2d-5um end to end — 4 observables matching their analytic solutions           complete
   Phase 1-B     soft-r3 end to end → compare the commonalities (judged together)            half a day   ⭐️ ← next
   Phase 1-C     abstract only what actually appeared twice + the Builder/Runner     half a day
   Phase 2       formalizing the golden tests (the trap already passes)              1 day
   Phase 1-D~    chain-bend → trap-drag → abp-rod                half a day each
   Phase 3       the raw-data Tiers B~D + on-demand analysis            1.5 days
   Phase 4       the rest of the post-mortem (the LLM narrative · rule promotion)              0.5 days  ← the automatic diagnosis is complete
   Phase 7       the intake schema                                   0.5 days
   Phase 8-rest  the hooks + the workflow + the subagents                       1 day
   Phase 5       the knowledge base                                        2.5 days
   Phase 6       parallel sweeps                                         1.5 days
```

**Phase 1 produces a working result every half day.** Since it drives a case rather than building the framework first,
changing direction midway throws little away.

---

## 17. The test strategy

| Layer | Content | Status (2026-08-03) |
|---|---|---|
| **golden physics** | proving the integrator's accuracy against an analytic solution | ✅ the trap `⟨x²⟩=kT/k` with an error <0.6% across 4 values of k, `⟨x²⟩·k` varying 0.28% · the 4 observables of 1-A |
| **verifying the documentation's code** ⭐️ | extract the skill snippets and run them | ✅ `scratch/verify_skill_snippets.py` (syntax 14/14 + the trap physics) |
| **the dt bias law** ⭐️ | confirming `bias = (dt/τ)/2` | ✅ `scratch/dt_convergence.py`, 4 points, to within 0.07% of theory |
| **automatic diagnosis** | the consistency of the equilibrium, drift, statistics and bias checks | ✅ `tools/postmortem.py` |
| **reproducibility** | a fixed seed → an identical result |  ✅ the 1-A re-run identical to the digit |
| **the non-dimensionalization round trip** | the physical → dimensionless → physical identity | partial — 1-A round-trips but `test_scale_invariance` (comparing 2 references) is unwritten |
| **the convergence study** ⭐️ | unchanged at half `dt` and double `N` | ❌ **not done** — stated in 1-A's `record.json.not_verified` |
| regression | a fixed-seed short run → a snapshot diff | ❌ |
| the schema | the SimSpec round trip, run_id stability | partial — the run_id content addressing works |
| the verification logic | a wrong spec → the exact error | partial — the separation-check gate actually blocked a smoke run |
| **distillation accuracy** ⭐️ | a 5-paper answer set | ❌ (0 papers) |
| agent evaluation | 10 requests, human-scored | ❌ |

**Fast iteration**: the case script has `--smoke` (N=200 short) · `--report` (the report only) · `--force` (re-run).
`--smoke` also has to pass the separation checks — the initial configuration actually got blocked by the statistics check.

---

## 18. Risks & mitigations

| Risk | Mitigation |
|---|---|
| the HOOMD API differing from expectation | demonstrate it all in Phase 0 before proceeding |
| **the LLM hallucinating a paper value** ⭐️ | the source (locator) mandatory · the tier system · human verification · a distillation-accuracy test |
| **misreading a sketch** ⭐️ | transcription first · `ambiguities` enforced · approval per field |
| the LLM producing a plausible but wrong physics spec | the verification layers + the golden tests + hardcoding the relations |
| dt being large and quietly inaccurate | a hard ceiling + runtime monitoring + a convergence-study utility |
| **the raw data exploding the disk** ⭐️ | tiered storage + estimating the volume in advance + an approval gate + a GC policy |
| **the KB getting contaminated with rubbish** ⭐️ | the tier system · the verification workflow · marking an entry `retracted` when refuted (not deleting it, preserving the history) |
| applying BD to a system that violates the BD assumptions | automatic `Re`/`St` checks |
| the MacBook CPU limit | cost estimation + an approval gate + a smoke mode |
| analysing without knowing the equilibration is insufficient | **always** plot the energy time series + the automatic equilibrium verdict |
| the API cost | prompt caching · short structured output · the local pipeline works without an LLM |

---

## 19. The technology stack

| Purpose | Choice | Reason |
|---|---|---|
| the simulation | HOOMD-blue 7.1.1 (`*cpu*`) | specified |
| the schema | Pydantic v2 | CLI validation + YAML serialization + (for batch distillation) directly into `messages.parse()` |
| **units** | **pint** | blocking unit mistakes in dimensional computation |
| the trajectory | GSD (+ `Burst`) | HOOMD native |
| the log | HDF5 (`hoomd.write.HDF5Log`) | large time series |
| the analysis | freud + numpy/scipy | it pairs well with HOOMD |
| tables | pandas + pyarrow (Parquet) | comparing sweeps |
| **the KB** | **SQLite + FTS5** | a vector DB is overkill at local scale, and SQL is accurate for coordinate search |
| the job queue | SQLite | no server needed |
| plots | matplotlib | minimal dependencies |
| **the agent runtime** | **Claude Code** | specified — vision, PDF, conversation and approval are all built in |
| **the integration method** | **CLI + Skill + Hooks** (not MCP) | §15.2 |
| calling the LLM directly | `anthropic` — **an optional dependency** `bdbot[distill]` | for batch distillation of papers only (§15.9) |
| the CLI | typer | a subcommand structure + `--format json` |
| tests | pytest | — |

---

## 20. Open questions

### Resolved (as of v0.4)
| # | Question | Conclusion |
|---|---|---|
| 3 | the sketch format | a photo of paper and pen, readable with the Claude Code Read tool. Low information density (2~3 parameters) |
| 4 | the default reference scales | fixed at `thermal` (σ=d, E=kT, τ=τ_B). **The unit system is not changed per system** |
| 4b | the separation-check thresholds | all unified at `10⁻²`. `dt/τ=1e-2` ⟺ a 0.5% bias (confirmed by measurement) |
| 6 | the interface | the Claude Code runtime + the CLI/Skills/Hooks (not MCP) |
| 9 | the scope | not restricted. Extended through the physics module registry (§5.6) |
| — | the medium | water (Newtonian). Viscoelasticity as a separate case (a `medium.*` module) |
| — | when the KB starts | from `record.json`; SQLite once there are 100 runs or the literature comes in (§7.0) |

### The remaining questions

**A. Anisotropic translational friction** — detailed in §10 below. A precondition of the `abp-rod` case. **A decision is needed.**

**B. When to do the convergence study** — 1-A did not re-run at half `dt` and substituted a comparison against the
bias law's prediction (stated in `record.json.not_verified`). That is justified for a linear system, but a nonlinear
system (from 1-B soft-r3 onwards) has no closed form so **a direct convergence study is needed**. Should it go into
1-B as a standard procedure, or be split off as a `tools/converge.py` utility?

**C. The observable priority** — the 5 to implement first in Phase 3.
Recommended: `msd`, `rdf`, `local_density_hist`, `cluster_size`, `force_distribution`.
(In 1-A, `P(x)`, `⟨x²⟩`, `C(t)` and `PSD` are already implemented inside the case code — candidates to be made common in 1-C)

**D. The paper sources** — feeding PDFs in directly vs fetching by arXiv/DOI. The latter needs a network-access policy.
If you designate the 5~10 papers to go in as the KB seed, the distillation prompt gets tuned to that style.
(currently 0 papers — not urgent until Phase 5)

**E. Notebook support** — if you are going to manipulate the analysis results in Jupyter, `bdbot.api` has to be
designed as first class during the Phase 1-C abstraction.

**F. The `mater_plan.md` filename** — it looks like a typo for `master_plan.md` but it is being kept as the name you
asked for. Should it be changed?

### A. Anisotropic translational friction — how should this be done? (**a decision is needed**)

The earlier judgment that binding as a rigid body produces anisotropy **was wrong on measurement** (the §5.6 note,
[`docs/hoomd_capabilities.md` §5.1](../../docs/hoomd_capabilities.md)). Both the rigid body and the bead chain give `γ⊥/γ∥ = 1.000000`.
The anisotropy is an effect of HI, and BD has no HI.

The actual impact on `abp-rod-2d-run-flip`:

| Observable | Accuracy |
|---|---|
| **MSAD** | ✓ exact (depends only on `γ_r,z`; being 2D there is only one rotation axis) |
| **the long-time MSD** (t ≫ τ_r) | ✓ exact (depends only on the isotropic average `γ̄`) |
| **the short-time MSD** (t < τ_r) | ✗ the anisotropy is lost |

Run-and-flip holds its orientation for a long time between flips, so the short-time anisotropy really can be visible.
Whether the sketch's "measure MSD, MSAD" was aiming at that signal is the criterion for the judgment.

| Option | Method | Cost | Result |
|---|---|---|---|
| **A (the recommended default)** | compute the isotropic average `γ̄` from the aspect ratio with Perrin/slender-body. Rotation exactly, via the `γ_r` tensor. The module declares the limitation | low | the MSAD and the long-time MSD exact |
| **B** | impose the anisotropy with a custom integrator (`md.half_step_hook`). The deterministic term works as a corrective force but **the noise term breaks fluctuation-dissipation** — it needs careful design | high | exact |
| **C** | introduce HI with `hoomd.mpcd` | very high | exact, effectively a different project |

**Starting with A and leaving B as a separate module** is what is recommended — being a module structure, attaching it
later does not touch the core. But if the short-time MSD anisotropy is **the core observable** of this research,
B has to go into Phase 1b, so please tell me the sketch's intent.


Looking at the 5 cases that came into `intake/`, the scope set in §1 (spherical colloids + spherical ABPs) **cannot handle a single one of them.**

| Case | The physics needed | The current scope | The HOOMD equivalent |
|---|---|---|---|
| `trap-2d-5um` | **an external harmonic trap** | ❌ | `md.force.Custom` or `md.external.field` |
| `trap-drag-2d-hex300` | a trap + **moving driving**, a hexagonal initial arrangement | ❌ | the above + moving the trap centre with `hoomd.variant` |
| `chain-bend-2d-oscill` | **bonds + a bending stiffness + time-dependent driving** | ❌ (out of scope) | `md.bond.FENE/Harmonic` + `md.angle` + `variant` |
| `abp-rod-2d-run-flip` | **an ellipsoid (non-spherical)** + **run-and-flip** (a discrete event) | ❌ | a rigid body or an anisotropic potential + **a custom updater** |
| `soft-r3-2d-A-sweep` | **a soft `r⁻³` potential**, an amplitude sweep | ❌ | `md.pair.Table` or custom |

`abp-rod` in particular affects the plan's physics formulas directly:
- **§6.1's `D_r = 3D_t/d²` holds only for a sphere.** An ellipsoid needs the Perrin friction factors, its
  translational diffusion is **anisotropic** (major/minor axis) and its rotational diffusion differs per axis too.
- **Run-and-flip is not an ABP.** An ABP is continuous rotational diffusion (`create_diffusion_updater`), whereas
  run-and-flip is **a discrete event** (a 180° reversal by a Poisson process). It cannot be done with HOOMD's built-ins and
  needs a custom updater.
- Whether the sketch's `τ_R = 0.5 s` is the rotational diffusion time or the flip interval is **ambiguous** (a typical `ambiguities` item).

**The proposal**: redefine the §1 scope as the 6 **physics modules** below and make each a plugin.
Designing `SimSpec` as a module combination covers all 5 cases and leaves extension open.

| Module | Content | Priority |
|---|---|---|
| `M1 pair` | WCA/LJ/Yukawa/**Table (an arbitrary r⁻ⁿ)** | essential |
| `M2 external` | **the harmonic trap**, a moving trap, walls | high (2 of the 5) |
| `M3 bonded` | bonds + **the bending angle** (a chain) | high (1) |
| `M4 active` | continuous ABP + **discrete run-and-tumble/flip** | high (1) |
| `M5 anisotropic` | **ellipsoids/rods** (Perrin friction, anisotropic diffusion) | medium (1 — the hardest) |
| `M6 driving` | **time-dependent driving** (`hoomd.variant`: oscillation, a ramp) | high (2) |

This **affects the schedule**: Phase 1b (the Builder) grows 1.5 days → 3~4 and M5 (non-spherical) may need its own Phase.

Three options:
- **(a)** extend the scope to cover all 5 cases — honest but +4~6 days on the schedule
- **(b)** easiest first: M1+M2+M6 first (the 2 trap cases), M3/M4/M5 later — a first result within 2 weeks
- **(c)** pick 1 representative and drive it end to end (a vertical slice), then extend — the fastest pipeline verification

**My recommendation is (c).** `trap-2d-5um` is the simplest (1 particle + a harmonic trap) and has an analytic solution, so
**it can be used directly as a golden physics test** — the position distribution of a BD particle in a trap is exactly
`P(x) ∝ exp(-kx²/2k_BT)` with `⟨x²⟩ = k_BT/k` and a relaxation time `τ = γ/k`. Verifying the whole route —
non-dimensionalization, verification, raw data, analysis, the post-mortem, the KB feedback — at minimal cost and then
attaching M3/M4/M5 is the safe order.

---

## Appendix A — the Phase 0 smoke test script

```python
# scratch/hello_bd.py — API demonstration, a minimal example
import itertools, math
import numpy as np, gsd.hoomd, hoomd

N_SIDE, PHI, KT, GAMMA, DT = 40, 0.5, 1.0, 1.0, 1e-4
N = N_SIDE ** 2
L = math.sqrt(N * math.pi / (4 * PHI))            # 2D area fraction → box length (σ=1)

a = L / N_SIDE
pos = np.array([[(i + .5) * a - L / 2, (j + .5) * a - L / 2, 0.]
                for i, j in itertools.product(range(N_SIDE), repeat=2)])

frame = gsd.hoomd.Frame()
frame.particles.N = N
frame.particles.position = pos
frame.particles.orientation = [(1, 0, 0, 0)] * N
frame.particles.typeid = [0] * N
frame.particles.types = ['A']
frame.configuration.box = [L, L, 0, 0, 0, 0]      # Lz=0 → 2D
frame.configuration.dimensions = 2

sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_snapshot(frame)

cell = hoomd.md.nlist.Cell(buffer=0.4)
lj = hoomd.md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode='shift')  # = WCA
lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)

bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=GAMMA)
integrator = hoomd.md.Integrator(dt=DT, methods=[bd], forces=[lj])
integrator.integrate_rotational_dof = False       # ABP extension: stated explicitly in preparation
sim.operations.integrator = integrator

# --- Tier A: positions/orientations ---
sim.operations.writers.append(
    hoomd.write.GSD(filename='traj_A.gsd', trigger=hoomd.trigger.Periodic(10_000),
                    mode='xb', dynamic=['property']))

# --- Tier B: per-particle forces (must be demonstrated in Phase 0) ---
plog = hoomd.logging.Logger(categories=['particle'])
plog.add(lj, quantities=['forces', 'energies'])
sim.operations.writers.append(
    hoomd.write.GSD(filename='traj_B.gsd', trigger=hoomd.trigger.Periodic(50_000),
                    mode='xb', logger=plog, dynamic=['property', 'momentum']))

# --- Tier L: global scalars ---
thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)
glog = hoomd.logging.Logger()
glog.add(thermo, quantities=['potential_energy', 'pressure'])
sim.operations.writers.append(
    hoomd.write.Table(trigger=hoomd.trigger.Periodic(10_000), logger=glog))

sim.run(100_000)
```

> The purpose of Phase 0 is to confirm **that this script runs as it is on v7.1.1**.
> If even one API signature differs, the §11 mapping table and the builder design get fixed then.
> In particular, storing per-particle forces through the `logger=` argument is the premise of Tier B, so verify it without fail.

---

## Appendix B — a sketch of the Claude Code integration outputs

### B.1 The `bdbot` CLI surface (everything Claude Code calls)

```
bdbot init-system --template <name> [--from-observation FILE]  → system.yaml skeleton
bdbot nondim system.yaml [--strategy thermal] [-o specs/x.json] → SimSpec + report
bdbot validate specs/x.json [--format json]                     → verification report (pass/fail by exit code)
bdbot estimate specs/x.json                                     → wall-time / disk estimated
bdbot run specs/x.json [--background]                           → run_id
bdbot sweep sweeps/y.yaml                                       → sweep_id (the whole range pre-checked)
bdbot status [run_id|sweep_id] [--format json]                  → progress/status
bdbot resume <run_id>                                           → restart from the checkpoint
bdbot analyze <run_id> --obs msd,rdf [--format json]            → observables.parquet
bdbot raw <run_id> --what forces --steps 3e6:4e6 [--who tracers] → Parquet extract + summary statistics
bdbot plot <run_id|sweep_id> --kind standard|phase|literature   → PNG paths
bdbot redimensionalize <run_id> --obs msd                       → values in physical units
bdbot kb search --tags 2D,ABP --coords Pe=40,phi=0.6 [--format json]
bdbot kb add entry.yaml   |   bdbot kb verify <id>   |   bdbot kb conflicts
bdbot postmortem <run_id> [--format json]                       → automatic diagnosis (without an LLM)
bdbot intake init <folder>  |  bdbot intake check <folder>
bdbot gc --dry-run
```

**The design rules**
- `--format json` on every command — so Claude Code can parse it easily
- pass/fail by exit code — used by the hooks for the verdict
- an error message **always** includes a concrete amendment (the §12.1 example)
- it works identically when a human uses it directly (no Claude Code dependency)

### B.2 The hook — enforcing a hard invariant (`.claude/hooks/guard_invariant.py`)

```python
#!/usr/bin/env python3
"""PreToolUse hook: block `bdbot run` on specs that bypass the dimensional layer.

Principle 3 (dimensions first) is enforced by the harness, not by a prompt.
stdin receives the hook input (JSON); when blocking it writes the reason + the amendment to stderr and exits 2.
"""
import json, re, sys
from pathlib import Path

payload = json.load(sys.stdin)
cmd = payload.get("tool_input", {}).get("command", "")

m = re.search(r"\bbdbot\s+run\s+(\S+)", cmd)
if not m:
    sys.exit(0)                                   # an unrelated command -- pass

spec_path = Path(m.group(1))
if not spec_path.exists():
    sys.exit(0)                                   # whether it exists is the CLI's judgment

spec = json.loads(spec_path.read_text())
missing = [k for k in ("derived_from", "scale_ledger") if not spec.get(k)]
if missing:
    print(
        f"BLOCKED: {spec_path} is missing {missing}.\n"
        "A principle 3 (dimensions first) violation -- a dimensionless spec must be derived from a PhysicalSystem.\n"
        "Fix: bdbot init-system --template <name> > system.yaml\n"
        f"      bdbot nondim system.yaml -o {spec_path}",
        file=sys.stderr,
    )
    sys.exit(2)                                   # 2 = blocked, and stderr is delivered to Claude

ledger = spec["scale_ledger"]
failed = [c for c in ledger.get("separations", [])
          if c["verdict"] == "fail"]
if failed:
    lines = "\n".join(
        f"  - {c['name']}: {c['ratio']:.2e} (limit {c['threshold']:.0e}) -- {c['message']}"
        for c in failed
    )
    print(f"BLOCKED: {len(failed)} scale-separation check failure(s)\n{lines}", file=sys.stderr)
    sys.exit(2)

sys.exit(0)
```

> The hook registration goes in `.claude/settings.json`. That file is safer written with the
> `update-config` skill than edited by hand.

### B.3 A slash command (`.claude/commands/bd-intake.md`)

```markdown
---
description: read a sketch folder, write an Observation and get human confirmation
---

The folder given as an argument: $ARGUMENTS

1. `bdbot intake init $ARGUMENTS` makes an `observation.yaml` skeleton.
2. Read every image in the folder with the Read tool.
3. skill `bd-intake` protocol, followed:
   - first transcribe **exactly as visible** into `raw_transcription`. Do not interpret.
   - then fill in `entities` and `stated_quantities`.
   - never invent a value **absent** from the sketch. Leave it `null`.
   - `ambiguities` and `unread_regions` are filled in without fail. If they are empty, look again.
4. `observation.yaml` is written.
5. `bdbot intake check $ARGUMENTS` checks the required fields.
6. Get confirmation from the user **per item** (AskUserQuestion). For ambiguous items in particular,
   present the interpretation candidates as options.
7. Update `observation.yaml` with the approved content.

The next step is `/bd-spec $ARGUMENTS`, but do not move on to it automatically now.
```

### B.4 A subagent (`.claude/agents/bd-reviewer.md`)

```markdown
---
name: bd-reviewer
description: adversarially reviews a simulation spec. Called before submitting a spec.
tools: Read, Bash, Grep
---

You are a sceptical simulation-physics reviewer. You are **not the person who made** this spec, and
you do not know their intent. You judge from the spec file and the dimensionless report alone.

Your task is not approval but **refutation**. Check the following in order:

1. Are the dimensionless numbers in a region where the target phenomenon can be observed? Compare against the literature
   (`bdbot kb search`).
2. Is there any scale-separation check with a margin under 5×? Could that
   distort the result?
3. Is the box sufficient relative to the length scales of interest (the correlation length, the persistence length)?
4. Is the observation time at least 100× the slowest timescale?
5. Does any property value at tier 2 or lower (unverified) govern the conclusion?
6. Is using Brownian dynamics for this system justified? (Re, St, the overdamped assumption)

For each item, judge it as one of **no problem / a concern / critical**, and for a concern or a critical
present a concrete amendment. If you are not sure, say "not sure".
Do not invent a problem you did not find.
```

### B.5 A skill (`.claude/skills/bd-physics/SKILL.md` — an excerpt)

```markdown
---
name: bd-physics
description: |
  Brownian dynamics simulation: the unit system, the scale ledger and the non-dimensionalization convention.
  Read it when defining a physical system, proposing parameters, interpreting a dimensionless number,
  or judging a scale separation.
---

# The dimensions-first workflow

... (the whole §6.1 scale list table of the master plan) ...
... (§6.2 reference-selection strategy table) ...
... (§6.3 dimensionless number = the ratio of two scales table) ...
... (§6.4 separation-check threshold table) ...

## Absolutely forbidden
- Do not fix the dimensionless values first and infer the physical system afterwards.
- Do not apply the spherical Stokes relation (D_r = 3D_t/d^2) to a non-spherical particle.
- Do not invent a property value that is in neither the sketch nor the literature.
```

### B.6 The permission configuration (the essentials of `.claude/settings.json`)

```json
{
  "permissions": {
    "allow": [
      "Bash(bdbot:*)",
      "Read(//Users/kyuhwan/Desktop/simulation_auto/**)",
      "Write(//Users/kyuhwan/Desktop/simulation_auto/intake/**)",
      "Write(//Users/kyuhwan/Desktop/simulation_auto/kb/**)"
    ]
  },
  "hooks": {
    "PreToolUse": [
      {
        "matcher": "Bash",
        "hooks": [
          {"type": "command", "command": ".claude/hooks/guard_invariant.py"},
          {"type": "command", "command": ".claude/hooks/guard_cost.py"}
        ]
      }
    ]
  }
}
```

Reading is allowed across the whole project, but **writing is restricted to `intake/` and `kb/`**.
`specs/` and `runs/` are writable only by the `bdbot` CLI, which structurally prevents a hand-made spec from mixing in
(the second line of defence enforcing principle 3).
