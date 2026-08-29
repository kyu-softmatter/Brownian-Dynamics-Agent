# Brownian-Dynamics Agent

An agent that turns a hand-drawn sketch into a defended Brownian-dynamics
result.

It reads a physical system out of a sketch, a note or a paper; fixes it in SI
units with a provenance on every number; derives a dimensionless specification;
runs it in HOOMD-blue; and files what it learned — **including the failures** —
back into a knowledge base the next run queries first.

The thing it is built to resist is specific. **A Brownian-dynamics simulation
always produces a number.** `g(r)` always plots, MSD always looks like a line,
and not diverging is not the same as being right. Every gate in this repository
exists because a plausible-but-wrong result got through once.

> **Public repository.** No copyrighted PDF, trajectory binary or unpublished lab
> asset is in any commit — the history starts from the merged tree, so there is
> nothing to scrub. See [NOTICE](NOTICE.md) for what was held back and where its
> content lives instead.
>
> **Merged from three predecessors** on 2026-08-28
> (`BD_agent` → `Simulation_bot` → `simulation_auto`). What came from where, and
> which seams are still open, is in
> [`docs/00-merge-decisions.md`](docs/00-merge-decisions.md).

---

## Architecture

Every stage either **reads** evidence out of the knowledge base or **writes**
evidence back into it. `R` marks a read, `W` marks a write.

```text
      INPUT     a hand sketch · a note · a paper · a text description
                (6 sketches are in intake/ — they are the actual inputs)
                                    |
                                    v
  +-------------------------------------------------------------------+
  |  S1  READ                       transcribe BEFORE you interpret   |
  |  observation / inference / assumption, split three ways.           |
  |  Anything absent from the sketch stays `null`.                     |
  |                                                                    |
  |  R  knowledge/wiki/systems/   is there a card for this system?      |
  |                               THE CARD OWNS the scales and gates    |
  |  GATE  dimension . boundary . driving fixed;                       |
  |        is the `question` falsifiable?                              |
  |        back-translate the spec to prose -- a human approves THAT,  |
  |        not the YAML, because reading your own words back is a      |
  |        thing humans are good at                                    |
  +---------------------------------+---------------------------------+
                                    v
  +-------------------------------------------------------------------+
  |  S2  PREDICT           write the answer down BEFORE running        |
  |  prediction.yaml  ->  SEALED.sha256                                |
  |                                                                    |
  |  Sealing is structural, not disciplinary: settings.json refuses to  |
  |  edit a sealed document, and a broken seal means the comparison     |
  |  table is not built at all.                                        |
  |                                                                    |
  |  R  wiki/benchmarks/   systems whose answer is already known       |
  |  GATE  >=1 quantitative prediction, each with a tolerance, a basis |
  |        and a ROLE -- and the design power to know, in advance,     |
  |        which items CANNOT be decided by this design                |
  +---------------------------------+---------------------------------+
                                    v
  +-------------------------------------------------------------------+
  |  S3  SPECIFY                    SI first. Always SI first.        |
  |  Every number carries provenance + tier + derived_from.            |
  |                                                                    |
  |  R  knowledge/    material properties, past decisions, 44 paper    |
  |                   and book distillations                          |
  |  GATE  no empty field . derived values recomputed and matched      |
  |                                                                    |
  |  BLOCKED here is a SUCCESS -- it names the ONE missing input.      |
  |  Inventing a number to reach READY is the only failure.            |
  +---------------------------------+---------------------------------+
                                    v
  +-------------------------------------------------------------------+
  |  S4  NON-DIMENSIONALIZE            *** the only contract ***       |
  |  specs/<run_id>.json  --  the ONLY thing passed downstream.        |
  |  The health layer never imports case code. It reads this.          |
  |                                                                    |
  |  GATE  round-trip < 1e-12 . ledger complete (4 required roles) .   |
  |        every dimensionless group really IS a ratio of two ledger   |
  |        entries . the inverse-transform anchor holds                |
  |                                                                    |
  |  run_id = hash(physics fields only). Editing a comment must not    |
  |  invalidate a run -- and changing d by 10x must not leave the id   |
  |  alone. Both directions have bitten. Both are now guarded.         |
  +---------------------------------+---------------------------------+
                                    |
           +------------------------+------------------------+
           v  a hard gate fails                              v  advances
  +----------------------------+            +----------------------------+
  |  REFUSED / BLOCKED         |            |  S5  RUN                   |
  |  names the missing input   |            |  @RUN.builder -> Build     |
  |  and what would supply it  |            |  L4 health: NaN, drift,    |
  |                            |            |  configurational temp,     |
  |  Warnings and thin margins |            |  step resolution           |
  |  are SHOWN BUT DO NOT      |            |                            |
  |  BLOCK -- a gate that      |            |  A case supplies ONLY      |
  |  passes silently is not a  |            |  build(spec). Loops,       |
  |  gate, and one that        |            |  guards and metrics.json   |
  |  refuses 80 of 83 specs    |            |  are common code.          |
  |  with 0 real failures is   |            |                            |
  |  worse than none. Both     |            |  Read skill bd-hoomd       |
  |  actually happened.        |            |  BEFORE writing a line.    |
  +-------------+--------------+            +-------------+--------------+
                |                                         |
                v                                         v
      back to the researcher              +----------------------------+
      -- no result is claimed             |  S6 VISUALIZE              |
                                          |  an uncaptioned figure     |
                                          |  CANNOT BE CREATED         |
                                          |  kT=0 and kT>0 side by     |
                                          |  side, because that shows  |
                                          |  an SNR problem faster     |
                                          |  than any number           |
                                          +-------------+--------------+
                                                        v
                                          +----------------------------+
                                          |  S7 VALIDATE               |
                                          |  PASS / FAIL /             |
                                          |  ** INCONCLUSIVE **        |
                                          |                            |
                                          |  The role decides what a   |
                                          |  mismatch MEANS:           |
                                          |   implementation_check     |
                                          |     -> a bug, fix it       |
                                          |   hypothesis               |
                                          |     -> A RESULT, report it |
                                          |   measurement              |
                                          |     -> the sim is the      |
                                          |        answer              |
                                          +-------------+--------------+
                                                        v
                                          +----------------------------+
                                          |  S8 CONCLUDE               |
                                          |  answer S1's question.     |
                                          |  Separate the assumptions. |
                                          |  State the confidence.     |
                                          |                            |
                                          |  W  the post-mortem, the   |
                                          |     finding, the dead-end, |
                                          |     the tooling lesson     |
                                          +-------------+--------------+
                                                        |
                                                        v
  +===================================================================+
  |  KNOWLEDGE BASE                                       knowledge/  |
  +===================================================================+
  |  wiki/systems/      (system x target dynamics) cards. The card     |
  |                     owns the non-dimensionalization and the        |
  |                     gates -- not the pipeline, not the case. 11    |
  |                                                                    |
  |  wiki/findings/     Q->A, and dead-ends. A dead end is worth as    |
  |                     much as a success, and must state a CAUSE      |
  |                     rather than a symptom. 23                      |
  |                                                                    |
  |  wiki/benchmarks/   systems with known answers, run AS REGRESSION  |
  |                     TESTS. Literature is the grader we don't have. |
  |                     Used only as reading, a paper wastes half its  |
  |                     value. 5                                       |
  |                                                                    |
  |  source/papers/     42 distillations -- equations converted into   |
  |  source/books/       our conventions, with what we could and       |
  |                      could NOT reproduce. 2 books, 56/56 claims    |
  |                      re-derived numerically                        |
  |                                                                    |
  |  entries/           126 tool-written entries. 44 of them are       |
  |                     `tooling` -- which is the honest measure of    |
  |                     how much of this is fighting instruments       |
  |                     rather than physics.                           |
  |                                                                    |
  |  runs/**/record.json   227 post-mortems. A run is not finished     |
  |                        when it exits; it is finished when its      |
  |                        post-mortem exists.                        |
  +---------------------------------+---------------------------------+
                                    |
                                    +--> R  feeds the next S1 and S3
```

### When the knowledge base changes

| Moment | Direction | What moves |
|---|---|---|
| S1 begins | **R** | is there a `wiki/systems/` card? If so it owns the scales and the gates |
| S2 sets a tolerance | **R** | `wiki/benchmarks/` — systems whose answer is known |
| S3 needs a parameter | **R** | material properties, past decisions, distillations — with tier and provenance |
| a run finishes | **W** | `record.json`, one lesson per run |
| **a `hypothesis` mismatches** | **W** | `wiki/findings/` — **this is a result, not a failure** |
| a route turns out not to work | **W** | a `dead-end` page, so the next attempt queries it first |
| a tool bites | **W** | `entries/` with `origin: tooling` |
| a paper or book is distilled | **W** | `source/` — with the numbers we could *not* reproduce stated too |
| **a judgment is made in conversation** | **W** | captured out of chat into a durable entry → [03](docs/03-knowledge-base.md#5--capturing-judgment-out-of-conversation), *the real purpose of this project* |

---

## Current status

**The engine is complete and eight cases have run; the discipline is not yet
wired to the engine.**

| | |
|---|---|
| Cases | **8**, all `READY` at L0 · L2 · L3. Six have produced runs |
| Runs | **278** specs · 261 run directories · **254** with `metrics.json` · 227 post-mortems |
| Code | `bdbot/` 21 modules (L0→L7) · `simbot/` 19 modules (S2/S6/S7/S8) · 8 case scripts · 74 verification scripts |
| Tests | **572 pass**, 2 skipped, ~11 s |
| Knowledge | 46 wiki pages · 42 paper + 2 book distillations · 126 entries |
| Agent layer | 6 skills · 9 model-tiered subagents · 4 rules |

The headline scientific result is
[`chain-bend-2d-dlvo`](docs/04-cases.md#chain-bend-2d-dlvo--the-central-result):
145 runs, 3 independent driving protocols, concluding that a colloidal chain held
together by DLVO central forces alone has **no bending stiffness and is
rheologically invisible** — attaching the chain changes the system's dissipation
by **0.996×**, i.e. not at all. That was the prediction the source paper itself
made before introducing JKR adhesion to explain its data; this project executed
the prediction side. Three candidate ways to recover stiffness without adhesion
(sliding friction, rolling resistance, DLVO tension) were then quantified and
excluded, which narrows *why* adhesion was needed.

**What is blocking progress is a seam, not missing code.** `bdbot` runs the
physics; `simbot` holds the prediction-sealing and the `INCONCLUSIVE` verdict;
**they do not call each other.** So the strongest discipline in the repository —
write the answer down first, hash it, refuse to build the comparison table if the
hash breaks — has never been applied to the results above. Those results are
defensible on their numbers and *not* defensible against the charge of having
been interpreted after the fact. Closing that is item 1 on the
[roadmap](docs/06-roadmap.md#7--what-would-most-improve-the-science-in-order).

**Where this is going.** The longer-term goal is to join this agent to
[**agentic-microscope**](https://github.com/kyu-softmatter/agentic-microscope) —
the same architecture pointed at the instrument instead of the integrator. One
decides what the system does, the other decides what the microscope can actually
record, and today they are consulted separately and can silently contradict.
Joined, a simulation would state the precision an experiment needs to see an
effect, a measurement would close the assumptions this repository currently
carries as tier-1 *choices* (`T = 300 K` foremost), and experiment would become
the fifth layer of evidence that
[02](docs/02-verification.md#3--result-verification--four-layers-of-evidence)
deliberately left unadopted.

The sharpest case for it is this repository's own central result. At the bead
diameter `chain-bend-2d-dlvo` actually used, d = 1.47 µm, its 22.3 σ separation
is a transverse displacement of **166.8 nm** without adhesion against
**9.39 nm** with JKR. Deciding *which model* means resolving the 157.5 nm
difference — a 15.7× margin against a 10 nm target precision, comfortable.
Deciding *whether the adhesive branch is separable from zero* means resolving
9.39 nm — **0.94×**, just under. The first question is settled by physics; the
second is settled by photon count and frame count. **Neither repository can
answer the second alone**, and today it is answered by consulting them separately
and trusting that the two `d` mean the same thing in the same units.

**Neither is finished, and coupling two moving targets would be a mistake** — so
it is future work, with a stated order of preconditions: sealing first, then a
shared quantity vocabulary, and never letting a simulated number set
`evidence: measured`.
→ [06 §6](docs/06-roadmap.md#6--future-work--joining-this-agent-to-the-microscope-agent)
· [the mirror of it, from the instrument's side](https://github.com/kyu-softmatter/agentic-microscope#future-work--joining-this-agent-to-the-simulation-agent)

Read [the pitfalls](docs/05-pitfalls.md) before implementing anything.

---

## What a gate looks like

Two properties matter more than the verdict. Unedited output:

```console
$ python -m bdbot.cli health --gate specs/chain-bend-2d-oscill__w85__d7c5a778ddba.json

gate — chain-bend-2d-oscill__w85__d7c5a778ddba   (L3 verdict: PASS (3 warnings))
  ⚠ soft warning [model] note: tau_p/tau_fast = 0.597 (limit 0.01)
        — does not block, but it is a statistics / finite-size limitation
  ⚠ soft warning [model] *angle force valid   min|theta-pi| = 7.26e-05
        (limit 0.00141421) — does not block, but it is a …
  ⚠ soft warning [statistics] quasi-static reached De(w_min) = 0.993
        (limit 0.1) — does not block, but it is a …
  ⚠ thin margin [model] linear elasticity    a/delta_max — only 2.1x to the limit
  ⚠ thin margin [integration] fastest mode resolved dt/tau_fast — only 1.0x to the limit
  OK — cleared to run
```

Read that as: the fastest mode is **not overdamped** (`τ_p/τ_fast = 0.597`
against a 0.01 criterion) · `angle.Harmonic`'s force is invalid here
(`min|θ−π| = 7.26e-05`, below its `1.41e-03` clamp) · the frequency sweep never
reaches the quasi-static plateau (`De(ω_min) = 0.993` against 0.1) · and two
margins are thin, one of them **exactly on** its threshold.

**It passes, and it says what is wrong anyway.** All five lines are real
problems: the fastest mode is outside the integrator's assumption, HOOMD's
bending force is invalid in this regime, and the frequency sweep does not reach
the limit the sketch asked about. None of them is a *hard* failure, so none of
them blocks — because a gate that refuses everything is worse than no gate. This
one once **refused 80 of 83 specs with zero real failures among them**, and
nobody noticed for weeks, because the runner never called the gate. An unwired
checker cannot be wrong out loud.

**A margin of 1.0× is reported as 1.0×.** `dt/τ_fast` is sitting exactly on its
threshold. Rounding that to "pass" is how a `dt` that is nominally fine produces
a result that is quietly wrong, so the ratio is printed rather than the verdict.

The refusal that matters most in this repository was of a different kind: the
`angle.Harmonic` line above was once a **hard** check that stopped the case
entirely. `md.angle.Harmonic` clamps `sin θ` at √2×10⁻³, and below that its force
becomes quadratic instead of linear — up to **96 % wrong** — while **its energy
stays 0.0000 % correct**, so no energy check can find it. All 23 angles of that
case's response profile were inside the broken region. The check was not deleted
when the case moved to a `force.Custom` implementation; it was **scoped to the
implementation** (`hard=(BENDING_IMPL == "angle_harmonic")`), so it returns the
moment anyone switches back. A constraint of an implementation should not be
recorded as a constraint of the system — and should not be thrown away either.

---

## Document map

| Document | Contents |
|---|---|
| [00 Merge decisions](docs/00-merge-decisions.md) | What came from which predecessor and why, per the three merge criteria. **The known seams**, stated rather than smoothed |
| [01 Architecture](docs/01-architecture.md) | The four layers, the eight stages with their gates, the L0–L7 module map, how to run it |
| [02 Verification](docs/02-verification.md) | ★ **How you verify when there is no grader.** The 7-rung process ladder, the 4 layers of evidence, evidence grades, `PASS-with-doubt`, and the roles that keep verification from destroying discovery |
| [03 Knowledge base](docs/03-knowledge-base.md) | The six kinds of knowledge and their contracts · provenance and tiers · content-addressed `run_id` · capturing judgment out of conversation |
| [04 Cases](docs/04-cases.md) | ★ **What has actually been run.** Eight cases, the numbers, and the several occasions on which a stated conclusion was later reversed |
| [05 Pitfalls](docs/05-pitfalls.md) | Traps that are silently wrong — in HOOMD, in the tooling, and in the reasoning |
| [06 Roadmap](docs/06-roadmap.md) | Measured status, what is blocking, what is deliberately not done, and what would most improve the science |
| [HOOMD capabilities](docs/hoomd_capabilities.md) | Measured capability matrix — what is actually buildable |
| [History](docs/history/) | The three predecessors' design documents, verbatim in Korean, as the provenance record |

**Code**

| Module | Layer | Status |
|---|---|---|
| [`bdbot/`](bdbot/) | L2 engine, L0→L7 | 21 modules. `nondim.py` is the single contract; `health.py` the numerical verdict; `run.py` the `@RUN.builder` assembly registry that all 8 cases use. Front end does not import `hoomd`, so specifying is fast. Promotion rule is only ever *"has it appeared twice?"* |
| [`simbot/`](simbot/) | L2 pipeline half, S2/S6/S7/S8 | 19 modules. Prediction sealing, `PASS/FAIL/INCONCLUSIVE` with design power, `REPORT.md` generation, figure generation that cannot emit an uncaptioned figure. ⚠️ **one runner only** — see the seam above |
| [`cases/`](cases/) | case physics | 8 scripts. Each supplies only `build(spec) -> Build`; the loops, guards and storage are common |
| [`campaigns/`](campaigns/) | sweep analyses | 18 scripts from the 3,856-run `soft-r3` campaign — finite-size scaling, hexatic window, seed sweeps |
| [`verify/`](verify/) | executable claims | 74 scripts. Every physics claim in the docs traces to one. A checker that has not been deliberately broken is not a checker |
| [`tools/`](tools/) | knowledge | `kb.py` (query/add), `postmortem.py` (a run is finished when this has run), `health.py` |
| [`.claude/`](.claude/) | L1 agent layer | 6 skills — 3 mutually-exclusive orchestrators (`bd-pipeline`, `bd-diagnose`, `bd-knowledge`) and 3 domain references the pipeline reads at a stage (`bd-intake`, `bd-physics`, `bd-hoomd`) · 9 subagents, model-tiered by whether the task needs judgment or is mechanical · 4 rules, each born from a dated accident |
| [`knowledge/`](knowledge/) | L3 | ⚠️ **two unmerged schemas** — `wiki/` Markdown and `entries/` JSON, read by different tools |
| [`intake/`](intake/), [`specs/`](specs/), [`runs/`](runs/), [`figures/`](figures/) | L4 | the inputs, the 278 contracts, the text-only run ledger, and the curated result figures |

```bash
python -m bdbot.cli status
```

---

## Why this is not a wrapper around HOOMD

**1 · Dimensions come first, and there is no path that skips them.**
Every system is fixed in SI with a provenance and a tier on each number, then
passed through a scale ledger. You cannot start from `Pe=40, φ=0.6` — the tool
will refuse, because a dimensionless number without an anchor cannot be inverted
back to an experiment. This is not pedantry: `T = 300 K` in this repository is
labelled tier 1 but is actually a *choice* inherited from a sketch with no
temperature, and because water's viscosity is 2.06 %/K sensitive, that mislabel
carries a −5 % to −14 % error into every timescale downstream. **The tier field
is the only reason that was findable at all.**

**2 · A mismatch with theory is a result, not a bug — but only if you said so
first.** Every prediction carries a role. If it was derived from the model you
implemented, a mismatch is a bug. If it comes from an assumption the simulation
does *not* impose — continuum, dilute limit, effective medium, a paper — a
mismatch is **the finding**. Without that distinction, a project that computes
systems *because they might differ from the standard picture* ends up filing its
discoveries as failures. The corollary is a design test: if the list of
assumptions the theory adds is empty, the case can validate but can never
discover. One case here is kept precisely as that negative example.

**3 · Macroscopic answers come from running the interactions, not from
inverting a formula.** If an interaction is physically needed it goes in as a
pair potential, bond or angle. And the reverse is equally forbidden: you do not
get the modulus of a particle-resolved system by substituting into a continuum
relation. `GSER  G* = K*/(6πa)` was **invalid** for a bead chain, because what a
bead feels is not a continuous medium but two neighbouring beads. A book
distillation later showed the stronger version: bulk stress presumes a volume
average over many particles *and* force-free particles, so a trapped single chain
breaks both — **our `K′` and `K″` are stiffnesses, and renaming them `G` would
make them wrong.**

---

## Sources

| Source | Status |
|---|---|
| 6 hand sketches (the actual inputs) | in [`intake/`](intake/), downscaled |
| HOOMD-blue 7.1.0, CPU, no MPI, no GPU | installed; capability matrix measured |
| 42 paper distillations | in-repo. ⚠️ 38 of 42 are the group's own published work, so *"the literature says…"* is narrower here than it looks |
| 2 book distillations (Leal 2026; Welty 5th ed.) | in-repo, 56/56 claims re-derived numerically |
| Copyrighted PDFs | **not published** — see [NOTICE](NOTICE.md) |
| Trajectories (542 MB of `.gsd`/`.npz`) | **not published** — regenerate from `spec + seed` |
| The experiment's actual temperature | **not obtained.** The one-number fix that would remove a −14 % worst case from every timescale |

**No licence is granted.** Published for reading; all rights reserved.
