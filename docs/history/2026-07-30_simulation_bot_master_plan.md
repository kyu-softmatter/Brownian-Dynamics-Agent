# Simulation Bot — Master Plan

> A Brownian Dynamics simulation agent. From one hand drawing to a verified conclusion.
>
> Status: **design settled v0.1 / the deterministic core complete / the first physics campaign run to the end** · last modified 2026-07-30
>
> **Scale** commits 28 · `simbot` 19 modules · `scripts` 18 · tests **562** ·
> systems cards 11 · findings 23 · production runs **3,856** (batch wall 85 min)

---

## The progress dashboard

**Notation convention** — `[O]` done and verified · `[~]` in progress · `[X]` not started · `[-]` out of scope for this version
> This dashboard is the single source of truth for progress. When work finishes, **update here before the code**.

### Milestones
| | Milestone | Status |
|---|---|---|
| `[O]` | **M0** the skeleton (env · directories · documents · knowledge · tests) | done |
| `[O]` | **M1+M2** the vertical slice — **S1→S8 all the way through, from a hand drawing** | **done.** With a harmonic trap instead of free diffusion. Of 9 predictions, 7 PASS / 2 INCONCLUSIVE / **0 FAIL** · 16 runs, wall **10.6 s** |
| `[O]` | **M2.5** the deterministic core — **making the run reproducible in code** | **done (2026-07-28).** `spec` `nondim` `io` `validate` `report` `policy` `session` `viz` + `cli.py`. The hand-drawing run reproduces with one line, `cli.py run` (**2.3 s**, bit-identical to the first run) |
| `[O]` | **M2.7** the L1 agent layer (`.claude/`) | **done (2026-07-28).** 3 skills + 5 reference documents + 9 subagents + `settings.json`. The structure is watched by `test_agent_layer.py` (64 tests) |
| `[~]` | **M3** extending the physics domains A→D | the trap branch + **the 2D soft-repulsive system run to the end** (§4.1). Domains B, C and D have cards only |
| `[O]` | **M4** strengthening verification | **done (2026-07-29).** the dt ladder · the seed ensemble · self-consistency · discriminating power · `converge` + **bootstrap · a finite-size scan · pre-computed design power · sequential pre-registration** |
| `[O]` | **M4.5** the first physics campaign — **soft repulsion in 2D** | **done (2026-07-29~30).** 6 runs · 3,856 production runs. The card went `draft` → 31 benchmarks. §4.1 |
| `[X]` | **M5** extending the input modalities | hand drawings done → next is comparison against experimental data (needed to close A2) |

### The result of the first full run (run `2026-07-28_trap-2d-5um_2dfb9d`)
| | |
|---|---|
| input | 1 hand drawing (2D optical tweezers, `R=5 μm`, `k=10 pN/μm`, `T=300 K`) |
| the settled answer | `⟨x²⟩ = 416.58 ± 1.85 nm²` · `τ_trap = 8.0567 ± 0.0300 ms` · `f_c = 19.7 Hz` |
| the highest precision | `τ_fit/τ = 0.9998 ± 0.0008` (**0.08 %**) · MSD `R² = 0.99998` |
| total computation | **10.6 s** (16 runs, 8 concurrent) |
| errors caught | 4 — an S0 knowledge value · the S3 truncation solver · an S5 GC detach · **an S7 correlated-sample KS false rejection** |

### M0 details
| | Item | Notes |
|---|---|---|
| `[O]` | create the `simulation_bot` conda env | python 3.12 · hoomd 7.1.0 · gsd 5.0.1 · numpy · scipy |
| `[O]` | review local CPU parallelism + measure the optimal core count | §7.3 · `wiki/findings/local-cpu-parallelism.md` |
| `[O]` | `master_plan.md` | this document |
| `[O]` | the progress dashboard | this section |
| `[O]` | `config/run_policy.yaml` — propose the optimal conditions + a human override | §7.4 |
| `[O]` | **`knowledge/` — ported from BD_agent** | §10. `source/papers` 42 papers · `wiki/{systems 11, findings 23, benchmarks 5, concepts 2}` (as of 2026-07-30) + the `wiki/CLAUDE.md` contract |
| `[O]` | the `CLAUDE.md` project conventions | centred on the first principle (propose, but ask when you do not know) |
| `[O]` | the sensitivity-analysis design (S7b) | §11 |
| `[O]` | the model-tiering design | §12 |
| `[O]` | store the first example input | `inputs/trap-2d-5um/sketch_01.jpeg` + sha256 |
| `[O]` | `environment.yml` + `wiki/techniques/env-log.md` | §7.1. The 11 latest packages compatible with HOOMD 7.1 |
| `[O]` | the `passive-sphere--harmonic-trap` card | the (system × dynamics) pair of the first example drawing |
| `[O]` | decide the scope of the CLI and session layers | **adopted in full** (session + converge + params) |
| `[O]` | git init + the initial commit | `4ac2a53`, 78 files. Independent of BD_agent. Currently 28 commits |
| `[O]` | **S1 reading + S2 prediction (the first hand drawing)** | `runs/2026-07-28_trap-2d-5um_2dfb9d/` |
| `[O]` | 4 `simbot/` modules | `units` `estimators` `forces` `guards` — all tested |
| `[O]` | **the test suite (per stage)** | **562 passing / 1 skip / 94 s.** `pytest -m "not slow"` is 14 s |
| `[O]` | measure and settle the HOOMD scheme | EM settled · the noise is uniformly distributed · the particles are independent · the configurational temperature |
| `[O]` | **the `simbot/` deterministic core** | `units` `estimators` `forces` `guards` `build` `run` `cutoff` `analysis/trap` `spec` `nondim` `io` `validate` `report` `policy` `session` `viz` — **16 modules, all tested** |
| `[O]` | **`cli.py`** | `run` `resume` `converge` `params` `calibrate`, all confirmed working |
| `[O]` | **`examples/trap-2d-5um/`** | a machine-readable `spec.yaml` + `prediction.yaml`. Reproduces the 10 derived values of the hand-written first run |
| `[O]` | **the `.claude/` agent layer** | 3 skills · 5 reference documents · 9 subagents · `settings.json`. §12.3 |
| `[O]` | **`simbot/viz.py`** | 5 figures generated automatically. **The caption and `shows` are forced at generation time**, and a skipped figure requires a reason |
| `[~]` | `simbot/analysis/` | `trap` and `structure` **done** (`structure` was opened by the campaign — RDF, `ψ₆`, Voronoi, `S(k)`, time-resolved, finite-size, bootstrap; 37 tests). Remaining: `msd` `microrheo` `active` `equilibration` — **when the system appears.** This principle was right: `structure` was made after the verification runs existed, and so the tests caught 3 bugs |

### Test status — **562 passing / 1 skip / 94 s** (`-m "not slow"` is 14 s)
| File | Target | Count |
|---|---|---|
| `test_s0_units.py` | units, constants, the scale round trip | 28 |
| `test_s2_estimators.py` | analytic identities, limits, the cost model | 32 |
| `test_s3_cutoff.py` | the `r_cut` proposal (WCA/LJ/Yukawa/Morse) | 36 |
| `test_s3_spec.py` | provenance, gates, derived-value regression, the YAML round trip | 47 |
| `test_s4_nondim.py` | the card scales, the round trip <1e-12, the dt constraints, the policy | 55 |
| `test_s5_forces.py` | `HarmonicTrap` force/energy/numerical derivative | 9 (1 skip) |
| `test_s5_guards.py` | the configurational temperature + guard firing | 15 |
| `test_s5_scheme.py` | B1 · B2 · B5 · B7 · B9 + reproducibility | 8 |
| **`test_s6_viz.py`** | **caption enforcement** · the twin axis · the skip reason · independent frames | **34** |
| **`test_s5_pair.py`** | **the pair-interaction runner** · the box shape · the Table potential | **42** |
| **`test_s7_structure.py`** | **time-resolved · the finite-size exponent · bootstrap · the phase reading** | **37** |
| `test_s7_validate.py` | tolerance parsing, the verdict, the power, the seal | 53 |
| `test_s8_io.py` | the hash, the run directory, **the seal**, provenance | 31 |
| `test_s8_report.py` | REPORT.md — does the bad news not fall out | 33 |
| `test_cli_session.py` | the session being append-only, the budget gate, `session run`, the full run | 54 |
| **`test_agent_layer.py`** | **the `.claude/` structure** — frontmatter, links, tiering, permissions | **64** |

### The pipeline stage implementation (§2)
| | Stage | Core module | Status |
|---|---|---|---|
| `[O]` | S1 Intake | — (LLM) | the protocol: `.claude/skills/bd-pipeline/references/s1_intake_drawing.md` |
| `[O]` | S2 Predict | `estimators.py` · `spec.Prediction` | the prediction YAML + the seal. A 9-item example |
| `[O]` | S3 Specify | `spec.py` `units.py` `cutoff.py` | provenance enforced · gates declared · derived values recomputed and compared |
| `[O]` | S4 Nondim | `nondim.py` `policy.py` | the scales per card · the round trip `1.6e-16` · 4 dt constraints |
| `[O]` | S5 Run | `build.py` `forces.py` `run.py` `guards.py` | concurrent batch execution + recording failed runs |
| `[O]` | S6 Visualize | `viz.py` | 5 figures + `06_figures.md`. A figure with no caption **cannot be made** |
| `[O]` | S7 Validate | `analysis/trap.py` `validate.py` | PASS/FAIL/**INCONCLUSIVE** + the design power |
| `[O]` | S8 Conclude | `report.py` `io.py` | `REPORT.md` generated automatically. **Narrating** the conclusion is the agent |

### The physics domains (§4)
| | Domain | Reference case | Regression criterion | Status |
|---|---|---|---|---|
| `[~]` | **A** free and confined BD | `examples/trap-2d-5um/` | equipartition `⟨x²⟩=kT/k` | **the confined branch only.** ★ **free diffusion `D*=1.00±0.03` is still missing** (Q9) |
| `[X]` | **B** microrheology | — | the Newtonian-fluid limit `G''=ηω` | no card |
| `[X]` | **C** active matter (ABP) | — | the ABP MSD analytic expression | the card is `draft`, the scale rules are unimplemented |
| `[X]` | **D** aggregating colloids | — | `B₂` vs `βU_min` | no card |
| `[O]` | **E** **2D soft repulsion `A/r³`** ★new | `runs/2026-07-29_soft-r3-*` (6 runs) | a perfect lattice `ψ₆=1` · the liquid exponent `p=1/2` | **the most developed domain.** 31 benchmarks · §4.1 |

> ★ **Domain E is this repository's largest body of knowledge** — 6 runs · 3,856 production runs ·
> 31 benchmarks · 6 methodology findings. The full text is in **§4.1**.

> ★ **Domain A's basic verification is empty.** The trap branch is verified down to equipartition, the relaxation
> time and the EM bias, but **free diffusion `D* = 1`** has never once been measured. It was M1's original DoD and
> was skipped because the hand drawing was a trap (§8). As a result:
> - the `brownian` route (`σ`, `τ_D`) of `CARD_SCALE_RULES` **has never been executed end to end**
> - **the `dt` displacement gate has bound 0 runs** (it is off in a trap system)
> - there is no verification target from which to build `analysis/msd.py`
>
> **Correction (2026-07-28, the 3 items above re-verified):** the second item is true **only as far as the pipeline goes**.
> On the route where `scripts/chain_bend.py:113` re-implements the same thresholds (`0.03`/`0.005`),
> the displacement gate **has already bound** — `runs/chain-bend/smoke/batch.json`:
> `binding: "force"`, `dt_force = 4.82e-6` vs `dt_diffusion = 4.5e-4`,
> measured `max|F*| = 1037.7` (the force constraint reduced `dt` 100-fold).
> ⇒ `simbot.nondim.choose_dt` still has 0 bindings, and **the gate logic has been duplicated in two places.**
> And `choose_dt`'s displacement constraint has `active=has_pair`, `has_pair = bool(spec.pair)`, while
> `SystemSpec` has no bond or angle field → **put a bond-only system on the pipeline and the gate goes quietly off.**
> The measurement above came out of exactly the system that gate is needed for.
> (The first and third items remain true on re-verification: the callers of `scales_brownian` are only `nondim.py:68` and
> `test_s0_units.py` · `simbot/analysis/` has only `trap.py`)

### Input modalities (§8 M5)
| | Modality | Status |
|---|---|---|
| `[O]` | **a photo of a hand drawing** ← the v1 goal | ✅ **run to the end.** The reading protocol + 1 real user drawing |
| `[~]` | a text description | it works if you write the spec YAML directly. There is no natural-language → spec route |
| `[X]` | an experimental screen or video | M5. **Needed to close A2 (`a` vs `R`) by measurement** |
| `[X]` | a paper PDF | M5 (`bd-lit-distill` does distil) |
| `[X]` | voice | M5 |

### Whether the additional proposals were adopted (§9) — the details are in §9
| | Proposal | Decision |
|---|---|---|
| `[O]` | 1. git init | ✅ 28 commits |
| `[O]` | 2. **sealing the prediction** | ✅ `io.py`. `shasum -c` compatible · if it breaks, no comparison table gets made |
| `[X]` | 3. a pilot run | unimplemented. Declared in the policy as `mandatory: true` but `cli.py` does not run it |
| `[O]` | 4. enforcing the unit suffix (`_si`/`_star`) | ✅ `Quantity.si` rejects strings and bools + a round-trip test |
| `[~]` | 5. a question budget | the convention is in `CLAUDE.md` and the skills. **There is no code enforcement** |
| `[O]` | 6. the **`bd-diagnose` skill** | ✅ `.claude/skills/bd-diagnose/` |
| `[~]` | 7. a parameter sweep (`sweep: [...]`) | **not in `spec.yaml` yet.** But the 6 `scripts/soft2d_*.py` actually ran `A`, `N` and seed sweeps (3,856 runs) — **the pattern is established, so there is now material to generalize** |
| `[~]` | 8. a run cache | the material exists (`spec.hash()`, `completed_stages()`). **There is no `spec_hash` lookup** |
| `[X]` | 9. an HTML report | unimplemented. Now that there are figures there is a benefit |
| `[X]` | 10. **a hand-drawing guide** | unwritten. 2 ambiguities in the first drawing + **the blank unit of `r` in `soft-r3`** (it has to be read by the Zahn convention `d = n^{-1/2}` for `Γ = π^{3/2}A` to hold — read differently and the physics changes). The grounds have grown |
| `[-]` | 11. a Langevin fallback | to be decided after review. The `overdamped` gate catches a violation and only advises |
| `[-]` | 12. an HI approximation (RPY) | out of scope for v1 |
| `[-]` | 13. direct comparison against experimental data | to be decided after review. Needed to close A2 |

---

## 0. The goal and the scope

### 0.1 The one-line definition
A Claude Code native agent that interprets material the user provides (v1: **a hand drawing**) and autonomously
carries a Brownian Dynamics simulation from **design → prediction → execution → verification → conclusion**,
accumulating the grounds for its judgments in a knowledge base along the way.

### 0.2 The settled design decisions
| Item | Decision | Notes |
|---|---|---|
| the chatbot runtime | **Claude Code native** | this conversation window is the chatbot. No API key needed |
| the physics engine | **HOOMD-Blue 7.1.0** | `md.methods.Brownian` (overdamped Langevin) |
| the execution environment | **local CPU alone** | Apple Silicon, no GPU. The N ≲ 10⁴ scale |
| the conda env | **`simulation_bot`** (new) | packages accumulated in stages, with a history in `knowledge/env_log.md` |
| the v1 input | **a photo of a hand drawing** | later extended to experimental video → paper PDFs → voice → text |
| the physics domains | all 4 (§4) | free/confined BD · microrheology · active matter · aggregating colloids |

### 0.3 Explicit non-scope (not done in v1)
- hydrodynamic interactions (HI): Oseen/RPY/Stokesian dynamics **absent**. The free-draining approximation.
  → the conditions under which this approximation breaks are stated in `knowledge/models/no_hydrodynamics.md`.
- GPU / MPI / cluster submission
- reactions (chemistry), flow-field coupling (CFD), full electromagnetic coupling

### 0.4 Design principles (running through every stage)
1. **The agent judges, and the core computes.**
   The LLM computing a number "in its head" is forbidden. Every number is obtained by calling a `simbot/` function.
   The LLM's role is to decide and record *which model to use, what values to assume, and why*.
2. **Seal the prediction first.**
   The S2 prediction document is sealed with a hash before S5 execution. Post-hoc rationalization is blocked structurally.
3. **Every number has a provenance.**
   One of `from_drawing` / `from_paper` / `from_knowledge` / `assumed` / `derived`.
4. **Treat units as a type.**
   The variable-name suffixes `_si` (physical units) vs `_star` (dimensionless) are enforced. Mixing them is caught by tests.
5. **Failures are outputs too.**
   A simulation that exploded, a wrong prediction, a misread drawing → all recorded in `knowledge/failures/`.
6. **Reproducibility is the default.**
   One run has to be fully restorable from `spec + seed + env + code hash` alone.

---

## 1. The system architecture

### 1.1 The 4 layers

```
┌──────────────────────────────────────────────────────────────┐
│ L1  the agent layer —  .claude/skills/ , CLAUDE.md            │
│     Claude Code performs it directly. Interpretation, judgment, reasoning, questions, records. │
│     It never does a deterministic computation; it calls L2.    │
└───────────────┬──────────────────────────────────────────────┘
                │ calls (python -m simbot.… / import)
┌───────────────▼──────────────────────────────────────────────┐
│ L2  the deterministic core —  simbot/                         │
│     Pure Python. No LLM. Exhaustively verified by pytest.      │
│     Unit conversion, non-dimensionalization, analytic solutions, HOOMD execution, analysis, plots, reports. │
└───────────────┬──────────────────────────────────────────────┘
                │ read/write
┌───────────────▼──────────────┐  ┌───────────────────────────┐
│ L3  the knowledge layer — knowledge/ │  │ L4  outputs — runs/<run_id>/│
│     System archetypes, parameter │  │     Every per-stage artifact, │
│     grounds, modelling grounds, failure │  │     the trajectory, figures, the report. │
│     cases, verification benchmarks. │  │     Self-contained, reproducible. │
│     ← accumulated, under version control │  │     ← a .gitignore target │
└──────────────────────────────┘  └───────────────────────────┘
```

**L3 is this project's real asset.** The L2 code can be rewritten, but
"why was η=1.2 mPa·s used for this material", "why did it blow up at dt=5e-5" disappear unless accumulated.

### 1.2 The directory structure

> Notation: `[O]` exists and is tested · `[X]` planned only · `[~]` partially implemented
> (measured 2026-07-28. This section has to match the actual tree — if it does not, it is not a design document but fiction)

```
Simulation_bot/
├── master_plan.md            [O] this document. The single source of truth for the whole design
├── CLAUDE.md                 [O] the project conventions the agent reads every session
├── environment.yml           [O] simulation_bot env reproduction
├── pyproject.toml            [O] pytest configuration (markers, paths)
├── cli.py                    [O] run · resume · converge · params · calibrate
├── README.md                 [X] usage for humans — not there yet
│
├── .claude/                  [O] L1 agent layer. §12.3–12.4
│   ├── README.md                 [O] the composition + the tiering table + the Q6 decision grounds
│   ├── skills/
│   │   ├── bd-pipeline/          [O] [main] the S1→S8 orchestrator
│   │   │   ├── SKILL.md              stages, gates, a checklist of prohibitions
│   │   │   └── references/       [O] s1_intake_drawing ★ · s2_prediction
│   │   │                             s3_s5_execute · s6_s7_validate · s8_knowledge
│   │   ├── bd-diagnose/          [O] diagnosing a blown-up run (the exclusion order)
│   │   └── bd-knowledge/         [O] knowledge/ search, add, tidy
│   ├── agents/                   [O] 9 subagents. Tiered by the model: frontmatter
│   └── settings.json             [O] interpreter allow-list + **refusing to edit a sealed document**
│
├── simbot/                   ← L2 deterministic core. 0 lines of LLM
│   ├── units.py              [O] physical constants, Scales, the per-card scale factory
│   ├── spec.py               [O] Quantity/SystemSpec/Prediction + the gate checks
│   ├── nondim.py             [O] per-card non-dimensionalization, dimensionless numbers, the 4 dt constraints
│   ├── policy.py             [O] run_policy.yaml loader + a deep merge of overrides
│   ├── estimators.py         [O] analytic solutions and scalings (the S2 prediction engine)
│   ├── cutoff.py             [O] r_cut proposal (WCA/LJ/Yukawa/Morse)
│   ├── build.py              [~] the trap snapshot only. Overlap removal is not there yet
│   ├── forces.py             [~] HarmonicTrap (md.force.Custom) only
│   ├── run.py                [O] the trap BD runner + concurrent batch execution
│   ├── guards.py             [O] NaN/displacement/configurational-temperature/fluctuation checks
│   ├── analysis/
│   │   ├── trap.py           [O] MSD fitting, seed aggregation, distribution tests
│   │   ├── msd.py            [X] the general MSD, D, the block-average error
│   │   ├── structure.py      [X] RDF, S(q), clusters, the density profile
│   │   ├── microrheo.py      [X] GSER → G'(ω), G''(ω)
│   │   ├── active.py         [X] ABP MSD crossover, MIPS verdict
│   │   └── equilibration.py  [X] the verdict on reaching equilibrium
│   ├── validate.py           [O] prediction vs measurement, PASS/FAIL/INCONCLUSIVE, the power
│   ├── report.py             [O] REPORT.md generation
│   ├── session.py            [O] the session state (turns append-only), set = an estimate only
│   ├── io.py                 [O] run directory, hashes, **the seal**
│   └── viz.py                [X] ★ absent. The S6 figures live one-off in scripts/
│
├── scripts/                  [~] one-off — code that should move into simbot
│   ├── trap_batch.py             → run.run_trap_batch absorbed it (a duplicate)
│   └── trap_analyze.py           → viz.py does not exist, so it still lives here
│
├── examples/                 [O] machine-readable reference cases
│   └── trap-2d-5um/
│       ├── spec.yaml             S3 specification (18 provenance fields)
│       └── prediction.yaml       S2 prediction, 9 items (the seal target)
│
├── knowledge/                ← L3 knowledge base. The schema is §10 (BD_agent is authoritative)
│   ├── source/papers/        [O] 42 literature distillations + an INDEX
│   └── wiki/                 [O] CLAUDE.md contract + systems 5 · findings 13
│       │                         concepts 3 · benchmarks 2 · questions 2
│       ├── systems/              ★ (system × target dynamics) cards — they own the non-dimensionalization and the gates
│       ├── findings/             Q→A + dead-end
│       ├── concepts/  techniques/  benchmarks/  questions/
│
├── inputs/                   [O] the raw material the user provides (gitignored, only .sha256 tracked)
│   └── <topic>/…
│
├── sessions/                 [O] the session state (gitignored)
│   └── <session_id>/session.yaml + spec_turnNN.yaml
│
├── config/run_policy.yaml    [O] the resource, tier and dt policy. A human's overrides: take precedence
│
├── tests/                    [O] pytest — 373 passing / 1 skip / 28 s
│   └── test_s0_units · test_s2_estimators · test_s3_cutoff · test_s3_spec
│       test_s4_nondim · test_s5_{forces,guards,scheme} · test_s7_validate
│       test_s8_{io,report} · test_cli_session
│
└── runs/                     [O] outputs (gitignored, but .md/.json/.yaml tracked)
    └── 2026-07-28_trap-2d-5um_2dfb9d/     the first hand-drawing full run (human + agent)
        2026-07-28_cli-e2e-test/           the same system reproduced with cli.py (2.3 s)
```

### 1.3 The run directory convention (self-sufficiency)

```
runs/<run_id>/
├── 00_input/            a copy of the raw material (the hand-drawing photo and so on) + sha256
├── 01_intake.md         observations/inferences/assumptions recorded separately
├── 01_intake.json
├── 02_prediction.md     ⚠ sealed. No modification after S5
├── 02_prediction.json
├── 03_spec.yaml         the full specification in physical units
├── 03_spec_rationale.md the source and the grounds for each value
├── 04_reduced.yaml      the dimensionless specification + the inverse-transform coefficients
├── 04_nondim.md         the conversion table
├── 05_run_manifest.json the code hash, the env hash, the seed, the HOOMD version, the wall time
├── traj.gsd             the trajectory
├── thermo.h5            the thermodynamic log
├── 06_figures.md        the figure list + captions
├── figs/*.png
├── 07_validation.md     the prediction vs measurement comparison table
├── metrics.json         the measured values + errors
├── 08_conclusion.md     the conclusion
└── REPORT.md            the whole summary (the final output a human reads)
```

`run_id = <ISO time>_<slug>_<the first 6 characters of the spec hash>`

---

## 2. The 8-stage pipeline — in detail

Each stage is defined as **input → processing → output → gate (the pass condition) → failure mode**.
Failing the gate means not moving to the next stage but reporting to the user or going back to the previous stage.

---

### S1. Intake — interpreting the material

**Input** the hand-drawing photo in `inputs/<topic>/` (+ the user's verbal explanation)

**Processing**
1. Read the image and list the following:
   - particles: the count (exact/approximate), size differences, kinds (distinguished by colour, hatching or labels), specially marked individuals (a probe and so on)
   - boundaries: the box border, walls (solid lines), periodic boundaries (dashed lines/arrows), slit/cylinder/sphere shapes, the dimensionality (2D/3D)
   - arrows: position, direction, length → **list the candidates for which it is: force / velocity / flow / time progression**
   - text: numbers, units, symbols (η, T, k, φ, v₀…), axis labels, captions
   - graphs: a hand-drawn expected curve (if present, grounds for comparison with the S2 prediction)
2. **The 3-way separation of observation/inference/assumption** — this is S1's core output.

   | Grade | Definition | Example |
   |---|---|---|
   | `observation` | read directly from the drawing | "about 30 particles", "solid-line walls on the left and right" |
   | `inference` | derived from the drawing + physics knowledge | "solid-line walls + dashed above and below → a slit geometry, periodic in y and z" |
   | `assumption` | absent from the drawing, filled in by me | "the medium is water, η=1.0 mPa·s, T=298 K" |

   Each entry gets a `confidence: high/medium/low` and a one-line justification.
3. **gaps** — list the information essential to the simulation that is absent, and decide how to handle it:
   `ask_user` (within the question budget) / `fill_from_knowledge` / `assume_and_flag` / `sweep` (a parameter scan)
4. Search `knowledge/systems/` → if a similar archetype exists, reuse it.

**Hand-drawing-specific rules** (`references/s1_intake_drawing.md`)
- **Do not trust the absolute sizes in a hand drawing.** What is trusted is ① the topology (what is inside/next to what)
  ② the ratios (particle:box ≈ 1:20) ③ the counts ④ the symmetry ⑤ the stated numbers and units.
- The absolute values of arrow thickness and length are meaningless. Use relative comparison only.
- For an ambiguous element, rather than interpreting arbitrarily, **state 2~3 candidates** and predict in S2 how the
  results differ per candidate.
  → The user then becomes able to tell immediately which interpretation is right.
- If the drawing has no scale information at all, leave φ (the volume fraction) as a free parameter and mark it as a sweep candidate.

**Output** `01_intake.md`, `01_intake.json`

**Gate**
- The required fields fixed: the spatial dimension `d`, the number of particle species, the boundary condition, whether there is driving/activity, and what the question is
- Is the `question` field in a falsifiable form (e.g. "how much does the diffusion slow down" ✅ / "what happens" ❌)

**Failure modes** → `knowledge/failures/intake_*.md`
- An arrow was read as a force when it was actually a velocity field
- A 2D drawing was read as a 2D simulation when it was actually a cross-section of a 3D one
- The hand-drawing particle count was mistaken for the actual N (the drawing is a sketch; N is the number statistically needed)

---

### S2. Predict — reasoning out the expected result

> **Write down the answer before the simulation.** This stage is what secures this project's scientific honesty.

**Input** `01_intake.json`, `knowledge/`

**Processing**
1. **Identify the governing physics** — which forces/timescales compete. What determines the result.
2. **Estimate the order of magnitude of the dimensionless numbers** — φ, Pe, κσ, T*=kT/ε, k σ²/kT, D_r τ_B … (`simbot.nondim`)
   Which regime each dimensionless number is in, and how far it is from the regime boundary.
3. **Quantitative prediction** — call the analytic solutions/scalings of `simbot.estimators`. For example:
   - Stokes–Einstein `D₀ = k_BT / 6πηa`
   - free diffusion `⟨Δr²⟩ = 2d D t`
   - optical tweezers `⟨Δr²⟩(t) = (2d k_BT/k)(1−e^{−t/τ_k})`, `τ_k = γ/k`, equipartition `⟨x²⟩ = k_BT/k`
   - sedimentation `ρ(z) ∝ e^{−z/ℓ_g}`, `ℓ_g = k_BT/(Δρ V g)`
   - ABP `⟨Δr²⟩ = 2dD_t t + (2v₀²/λ²)(λt − 1 + e^{−λt})`, `λ = (d−1)D_r`
     `D_eff = D_t + v₀²/[d(d−1)D_r]`
   - the long-time diffusion reduction of a dense system, the MIPS phase boundary and so on are literature correlations (with the source stated)
4. **Sealed in a falsifiable form** — each prediction must have these 4 elements:

   ```yaml
   - quantity: D_long / D_0
     value: 0.42
     tolerance: "±25%"          # outside this is a FAIL
     basis: "Batchelor + φ=0.35 semi-dilute correction, knowledge/validation/dense_diffusion.md"
     discriminates: "HI being ignored -- whether that is justified"
   ```
5. **Alternative scenarios** — write down in advance the ways the prediction could be wrong, and the signal that would then appear.
   (e.g. "if dt is too large D gets overestimated → if D changes on a half-dt re-run it is a numerical problem")

**Output** `02_prediction.md`, `02_prediction.json`
→ the sha256 of both files is recorded in `05_run_manifest.json` to **seal** them. S7 verifies this hash.

**Gate**
- ≥ 1 quantitative prediction, each with a `tolerance` and a verdict criterion
- If the drawing has an expected curve the user drew, state the agreement or disagreement with it

**Failure modes** → `knowledge/failures/prediction_*.md`
- Getting the order of magnitude of a dimensionless number wrong and misjudging the regime
- Applying a literature correlation outside its range of applicability
- Setting the tolerance so wide that any result PASSes (neutering the verification) — **forbidden, a review target**

---

### S3. Specify — making the system concrete (in physical units)

**Input** `01_intake.json`, `knowledge/parameters/`, `knowledge/systems/`

**Processing** Fill in a complete `SystemSpec` (in SI units). No empty field may remain.

| Group | Fields |
|---|---|
| geometry | `dim`, `box_si` (Lx,Ly,Lz), `boundary` (pbc/wall/slit/cylinder/sphere) |
| particles | `species[]`: `name, N, radius_si, mass_si, density_si, charge, active` |
| medium | `T_si` (K), `eta_si` (Pa·s), `rho_fluid_si`, `epsilon_r` |
| friction | `gamma_si` = 6πηa (Stokes) or an explicit value; whether a near-wall correction applies |
| interactions | `pair[]`: the potential per type pair (WCA/LJ/Yukawa/Morse/DLVO) + parameters + `r_cut` |
| external | `traps[]`, `gravity`, `shear_rate`, `field` |
| activity | `v0_si`, `D_r_si` or `tau_r_si` |
| time | `t_total_si`, `t_equil_si`, `dump_interval_si`, `thermo_interval_si` |
| numerics | `seed`, `dt_policy` |

**A provenance is mandatory on every field:**
```yaml
eta_si:
  value: 0.890e-3          # 298.15 K. 20 °C value (1.002e-3) -- do not confuse
  unit: "Pa*s"
  provenance: from_knowledge
  source: "knowledge/parameters/water_298k.md"
  note: "298.15 K pure water, IAPWS"
```

**Automatic physical-validity checks** (`simbot.spec.validate`)
- φ < 0.64 (RCP), and in 2D φ_A < 0.9
- the Reynolds number `Re = ρ v a/η ≪ 1` (the BD premise)
- the inertial timescale `τ_i = m/γ ≪ dt` (the overdamped premise) — advise Langevin on a violation
- the Debye length vs the interparticle distance, consistent
- `r_cut` < L/2 (the minimum image)
- for an active system: `v₀ τ_r` (the persistence length) vs the box size

**Output** `03_spec.yaml`, `03_spec_rationale.md`

**Gate** every field filled + a provenance present + all the validity checks passing (or an explicit exception approved)

**Failure modes** → `knowledge/failures/spec_*.md`
- Computing γ as 6πηd instead of 6πηa (confusing the radius and the diameter) — **the most common mistake**
- Using 3D Stokes friction in a 2D simulation (state it if deliberate)
- The box being smaller than the persistence length, giving finite-size artefacts

---

### S4. Nondimensionalize — non-dimensionalization

> The detailed convention is §5. Here it is only the pipeline view.

**Input** `03_spec.yaml`

**Processing**
1. Choose the 3 reference scales (default: `L*=σ`, `E*=k_BT`, `T*=τ_B=σ²/D₀`)
2. Convert every SI value → a dimensionless value, and store the inverse-transform coefficients
3. Compute the whole set of dimensionless numbers and tabulate them
4. **Choose dt** — adopt the minimum of §5.4's 4 constraints, and record the grounds
5. The round-trip conversion test: the relative error of `to_reduced → to_si` < 1e-12

**Output** `04_reduced.yaml`, `04_nondim.md` (the conversion table: quantity | SI | dimensionless | the inverse-transform coefficient)

**Gate** the round-trip error passing + all the dt constraints satisfied + the dimensionless-number table complete

**Failure modes** → `knowledge/failures/nondim_*.md`
- Confusing HOOMD's time unit (τ_LJ = σ√(m/ε)) with the Brownian time (τ_B = σ²/D₀) → **fatal**
- Omitting the dimensional factor (2d) in the definition of τ_B in 2D
- Setting both kT and ε to 1 and thereby fixing T*=1 (a bug unless deliberate)

---

### S5. Run — executing HOOMD-Blue

**Input** `04_reduced.yaml`

**Processing**
1. **The initial placement** (`simbot.build`) — lattice/random/non-overlapping insertion. If there is overlap, remove it with a soft pushoff (a gradual increase in σ).
2. **The pilot run** (active by default) — run only 0.5% of the main run's steps to confirm
   ① no guard violations ② the expected wall time ③ the validity of dt. If there is a problem, halt here.
3. **Forces and the integrator** (`simbot.forces`, `simbot.run`)
   - `md.methods.Brownian(filter, kT, default_gamma)` — D = kT/γ
   - active: `md.force.Active` + `md.update.ActiveRotationalDiffusion`
   - walls: `md.external.wall.{LJ,Morse,Yukawa,ForceShiftedLJ}` + `hoomd.wall.{Plane,Sphere,Cylinder}`
   - optical tweezers / harmonic confinement: an `md.force.Custom` subclass (not built in)
   - gravity / an electric field: `md.force.Constant` or `md.external.field.Electric`
4. **Equilibration → production**, separated. The equilibration verdict is `simbot.analysis.equilibration`.
5. **Runtime guards** (`simbot.guards`) — checked every `thermo_interval`, halting immediately on a violation and saving a diagnosis:
   - the position and force of a NaN/Inf
   - the maximum displacement per step > 0.1σ (a sign of an explosion).
     ⚠ **Do not assume the displacement distribution is Gaussian** — HOOMD's noise is uniformly distributed (`max/σ = √3`)
   - wall penetration (a particle outside the boundary)
   - **the configurational temperature** `kT_conf = ⟨|∇U|²⟩/⟨∇²U⟩` matching the input `kT`
     — measured `1.00382 ± 0.00480` (a harmonic trap). The forces are already computed so the extra cost is almost nil
   - a divergence of the pressure or the potential energy

   > ❌ **Deleted: "the kinetic-energy temperature deviates from the target kT"** — it was a guard that could not work.
   > HOOMD's `Brownian` does not integrate the velocity but **draws it from the target distribution every step.**
   > `kinetic_temperature` merely echoes the input back and cannot deviate systematically.
   > Grounds: [`findings/hoomd-brownian-scheme-and-noise.md`](../../knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md)
   >
   > **`kT` itself is still a first-class input** — `U/kT` sets the Boltzmann weight and it enters every
   > dimensionless number (`k* = kσ²/kT`, `T* = kT/ε`, `Pe = Fσ/kT`). What is unusable is
   > only reading the kinetic energy back, and the configurational temperature takes that place.
6. **Logging** — the GSD trajectory (positions, orientations, image flags), `ThermodynamicQuantities` → HDF5.
   Reproducibility: the seed, the HOOMD version, the `simbot` code hash, the env hash, the spec hash, the wall time.

**Output** `traj.gsd`, `thermo.h5`, `05_run_manifest.json`

**Gate** completion + no guard violation + equilibrium reached + the production stretch ≥ 10× the required correlation time

**Failure modes** → `knowledge/failures/run_*.md`
- In an ABP, not setting `moment_inertia` to 0 and using `Brownian`, so the rotational degree of freedom is double-integrated
- `r_cut` > L/2 violates the minimum image
- dt being large so the WCA core is penetrated (the particles overlap and then it explodes)
- Starting without removing the initial overlap → diverging in 1 step

---

### S6. Visualize — visualization

**Input** `traj.gsd`, `thermo.h5`

**Processing**
1. **The mandatory diagnostic set** (domain-independent, always generated)
   - the MSD log–log + a slope-1 reference line
   - the RDF g(r)
   - the thermo time series (PE, P, the temperature estimate)
   - the per-step displacement distribution (a post-hoc guard check)
   - a triptych of initial/middle/final snapshots
2. **Domain-specific plots** — see the §4 table
3. **Animation** — matplotlib for 2D, fresnel ray tracing for 3D. GIF/MP4.
   (ffmpeg not installed → GIF first, MP4 after installing ffmpeg if needed)
4. **Dual labelling of every axis in dimensionless values and physical units** (`t/τ_B` above, `t [ms]` below)

**Output** `figs/*.png`, `figs/*.gif`, `06_figures.md` (each figure's caption = what the figure is trying to show)

**Gate** the axis labels, units and captions all present. A figure with no caption is not accepted as an output.

---

### S7. Analyze & Validate — analysis and verification

**Input** `traj.gsd`, `02_prediction.json` (the seal hash verified), `04_reduced.yaml`

**Processing**
1. **Seal verification** — confirm that the sha256 of `02_prediction.*` matches the manifest. On a mismatch, **halt and report**.
2. **Quantitative measurement** (`simbot.analysis`) — every measured value accompanied by a statistical error
   - MSD fitting → D (a block average + a bootstrap confidence interval)
   - structure: g(r), S(q), the cluster size distribution, the density profile
   - microrheology: GSER (Mason) → G'(ω), G''(ω)
   - active: the time of the MSD crossover, D_eff, the MIPS verdict indicator
3. **The prediction vs measurement comparison table** — `PASS` / `FAIL` / `INCONCLUSIVE` per entry
   `INCONCLUSIVE` is when the statistical error is larger than the tolerance so no verdict is possible → a longer run is needed
4. **Numerical convergence checks**
   - the dt ladder: re-run at dt/2 and see whether the change in the measured value is < the statistical error (default: done for 1 key indicator)
   - finite size: re-run with N or L multiplied by √2 (on request)
   - finite time: whether the measured values from the first and second halves of the trajectory agree
5. **Comparison against known analytic solutions (sanity)** — always run
   - for dimensionless free diffusion, `D* = 1.00 ± stat`
   - for harmonic confinement, equipartition `⟨x*²⟩ = 1/k*`
   - for an ABP, superposition with the analytic MSD curve
6. **A cause hypothesis for every FAIL item** — 4 classes: `numerical` / `modeling` / `interpretation` / `analysis`

**Output** `07_validation.md`, `metrics.json`

**Gate** the sanity checks passing + a cause hypothesis present for every FAIL

---

### S8. Conclude — summarizing the conclusion + committing the knowledge

**Input** the outputs of every earlier stage

**Processing**
1. **Question → answer** — answer S1's `question` directly. One paragraph.
2. **3 lines of grounds** — which measurements support that answer.
3. **Confidence and limits** — ignoring HI, finite N/t, the dt convergence, the influence of the assumptions.
4. **Proposing the next experiment** — 1~2 minimum-cost experiments that would firm up or refute this conclusion.
5. **Updating knowledge/ (mandatory, cannot be omitted)**
   - `systems/`: a summary of this system archetype (or an update of the existing entry)
   - `parameters/`: the grounds for a newly settled parameter
   - `models/`: the grounds for a modelling decision (why this potential/approximation)
   - `failures/`: every failure that occurred (the trivial ones too)
   - `validation/`: the benchmark values confirmed this time → candidates for promotion to a pytest regression test

**Output** `08_conclusion.md`, `REPORT.md`, N knowledge entries

**Gate** at least 1 knowledge entry added or updated. `REPORT.md` self-contained, including links to the figures and numbers.

---

## 3. The data model

**Implementation complete** (2026-07-28). Dataclass-based, no pydantic. Where it diverged from the design is marked ★.

```python
# simbot/spec.py
Quantity:        value, unit, provenance, basis, confidence, ambiguity,
                 sensitivity, affects[], written_by
                 ★ Provenance was not put in a separate dataclass — with only 5 fields,
                   another wrapper layer makes the YAML two levels deeper and unwritable by hand
Species:         name, n_simulated, radius_si, density_si, n_physical, charge, active
Geometry:        dim, boundary, box_si | box_over_ref
                 ★ box_over_ref added — a trap system gives the box as a multiple of ℓ_trap
Medium:          T_si, eta_si, rho_fluid_si, species
Friction:        model, gamma_si, wall_correction, note
PairInteraction: type_a, type_b, potential, params, r_cut_si
ExternalField:   kind, params{Quantity}, implementation, note
                 ★ traps/gravity/shear/active were merged into one rather than kept as separate fields —
                   adding a field per kind of external makes the schema differ per system
Timing:          equil_in_tau, prod_in_tau, sample_interval_in_tau, target_precision
Numerics:        dt_star, seed_base, n_seeds, integrator, scheme, noise_distribution
Gate:            status(required|pass|fail|off|applicable|unknown), reason
                 ★ new — the card switches the gates on and off. off requires a reason
SystemSpec:      card, question, geometry, species[], medium, friction, pair[],
                 external[], timing, numerics, gates{}, tier, notes[]
PredictionItem:  quantity, value, tolerance, basis, discriminates, unit,
                 competing_value, note        ★ competing_value added — for the power calculation
Prediction:      items[], regimes, alternatives[]
                 ★ sealed_hash was not put here — the seal is a file hash, so
                   a document cannot contain its own hash. io.write_seal owns it

# simbot/nondim.py
Scales:          length_si, energy_si, time_si, origin     (units.py)
DtConstraint:    name, dt_si_max, active, basis, off_reason     ★ new
DtChoice:        dt_si, dt_star, dominant, constraints[], logged{}   ★ new
ReducedSpec:     card, scales, dim, n_particles, box_star, kT_star, gamma_star,
                 D_star, sigma_star, dt_star, dt_dominant, k_star,
                 equil_steps, prod_steps, sample_interval_steps, groups{}, logged{}
                 ★ sigma_star is not 1 — in the trap card it is 491.358

# simbot/validate.py
Tolerance:       kind(relative|absolute|lower_bound|upper_bound), magnitude, text
Measurement:     quantity, value, stat_err, method, n_samples, spread, unit
ValidationRow:   quantity, predicted, measured, tolerance, verdict, stat_err,
                 deviation, deviation_rel, sigma, design_power,
                 samples_needed_for_3sigma, cause_class, reason, note, flags[]
                 ★ design_power·samples_needed_for_3sigma·flags added

# simbot/io.py
RunDir:          path + RUN_LAYOUT key access
SealVerdict:     ok, changed[], missing[], unsealed[], entries{}   ★ new
RunManifest:     ★ dataclass not used; build_manifest() → dict instead.
                 The fields include an env version table and so are variable; a fixed schema gets in the way
```

**Invariants** (enforced by tests)
- a `*_si` field is always a `Quantity` (with a unit). A `*_star` field is always a pure float.
  → `Quantity.si` rejects strings and bools (`test_s3_spec.py`)
- the `SystemSpec` → YAML → `SystemSpec` round-trip error is **0** (the derived values bit-identical)
- the `SystemSpec` → `ReducedSpec` → SI round-trip error is **< 1e-12** (measured `1.6e-16`)
- every `Quantity` has a `provenance`. **It is written in the file even when it is the default** —
  omit `assumed` and "I assumed it" becomes indistinguishable from "I forgot to write it"
- a field with `provenance ∈ {inference, assumed}` cannot be written by a cheap model (§12.2)

---

## 4. Physics domain coverage

| # | Domain | Model | The HOOMD construction | Key dimensionless numbers | The verification criterion (analytic/literature) |
|---|---|---|---|---|---|
| **A** | free and confined BD | point/spherical particles, WCA excluded volume | `Brownian` + `pair.LJ` (WCA mode) + `external.wall` | φ, `kσ²/kT`, `σ/ℓ_g` | free diffusion `D*=1` · equipartition `⟨x²⟩=kT/k` · the sedimentation exponential `ℓ_g` · the slit plateau MSD |
| **B** | microrheology | a probe + a background medium/network; passive (thermal fluctuation) and active (tweezer towing) | `Brownian` + `force.Custom` (the trap) + the background pair | `Pe = Fσ/kT`, `kσ²/kT`, `ωτ_B` | the Newtonian-fluid limit: `G'=0, G''=ηω` · GSER round-trip consistency · the trap corner frequency `f_c = k/2πγ` |
| **C** | active matter | ABP: self-propulsion + rotational diffusion | `Brownian` (`moment_inertia=0`) + `force.Active` + `update.ActiveRotationalDiffusion` | `Pe = v₀/(σD_r)` or `v₀τ_r/σ`, φ | the single-particle ABP MSD analytic expression · `D_eff = D_t + v₀²/[d(d−1)D_r]` · the MIPS phase boundary (2D: Pe≳40–60, φ≳0.4) |
| **D** | aggregating colloids | Yukawa (DLVO) / Morse attraction + WCA | `pair.{Yukawa,Morse,DLVO,ExpandedLJ}` | `T*=kT/ε`, `κσ`, `βU_min`, φ | the two-body binding probability vs `βU_min` · the second virial coefficient `B₂` at low φ · the initial slope of the Smoluchowski aggregation rate · the RDF contact peak |
| **E** ★ | **2D soft repulsion** | point particles `U/kT = A/r³` (**no** hard core) | `Brownian` + `pair.Table` (a power-law table) | **`Γ = π^{3/2}A`** · `n* ≡ 1` · `βU(r_cut)` · **`η₆ = 4p`** | a perfect lattice `ψ₆ = 1` · defects `0` · the liquid exponent `p = 1/2` · the Zahn phase diagram `[source, unreproduced]` |

Each domain keeps **1 verified reference case** in `examples/`, and its numbers are fixed in `tests/` as a regression test.

### 4.1 ★ The first physics campaign — the 2D soft-repulsive system (2026-07-29~30)

Starting from one hand drawing, `soft-r3-2d-A-sweep`, **6 runs · 3,856 production runs**.
The card [`soft-repulsive-2d--equilibrium-structure`](../../knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md)
went from `draft` to this repository's largest body of knowledge, with 31 benchmarks (S1–S31).

| Run | What was asked | Runs | Result |
|---|---|---|---|
| `soft-r3-2d-A-sweep` | an `A` sweep, 4 seeds | 40 | discovered that the box shape perturbs the initial-condition comparison |
| `soft-r3-time-resolved` | **when** does the structure form | 60 | `τ_relax = 0.03–0.10 τ_d` · the physical scale (`σ=5 µm`) attached for the first time |
| `soft-r3-relax-seeds` (2 stages) | `τ(A=1)` vs `τ(A=0.1)` | 3,584 | **INCONCLUSIVE** (short of the pre-registered `3σ`) |
| `soft-r3-nconv` | `N` convergence + the `ψ₆` finite size | 12 | **`A=10` is not hexatic** (16/16 PASS) |
| `soft-r3-fss` | verifying the **form** of the power law | 48 | the form confirmed for `A=0.1` and `A=10` (`χ²/dof` `0.55` and `0.58`) |
| `soft-r3-hexwin` | the Zahn hexatic window | 72 | all 3 points inside the window are **isotropic liquid** · the truncation error biases the exponent by `2.9σ` |

**The physical conclusion (the current best values)**

| `A` | `Γ` | `η₆` | Phase | Grounds |
|---|---|---|---|---|
| 0.1 | 0.56 | `2.03 ± 0.08` | isotropic liquid | `p = 0.508 ± 0.020` (`0.4σ` from the liquid value `0.5`) |
| 1 | 5.57 | `1.86 ± 0.05` | isotropic liquid | the form verification is unresolved, SE-limited |
| 10 | 55.68 | **`2.06 ± 0.23`** | isotropic liquid | the best value at `r_cut=7.80`. Hexatic rejected at `7.9σ` |
| 10.2–10.8 | 56.8–60.1 | `1.5–2.2` | isotropic liquid | **inside the Zahn window and yet there is no hexatic** |
| 31.6 | 176.0 | `0.27 ± 0.25` | crystal | defects `0.020` |

⇒ **the liquid→crystal bracket is `A = 10.8–31.6`. No hexatic was observed.**
⚠ But **the boundary cannot be found from a random start** (supercooling). It is a value bracketed from above only.

### ★ The 6 **methodology** findings this campaign produced — they apply to other cards too

| finding | Rule |
|---|---|
| [[coarse-sampling-hides-the-whole-transient]] | `stride ≲ τ_relax/5`. Estimate the relaxation time **first** |
| [[fraction-threshold-flips-meaning-between-per-frame-and-aggregate]] | if the fraction threshold is below `1/N` there is no threshold |
| [[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]] | to erect a `3σ` from a 4-seed SE you need `t(3)=5.84` |
| [[low-seed-pilots-give-optimistic-design-power]] | compute the sample count **in advance**. The bias indicator is **width/noise ≳ 20** |
| [[order-parameter-magnitude-cannot-identify-a-phase]] | neither the magnitude nor the slope is enough on its own — look at both |
| [[provenance-must-have-one-definition-and-three-capture-points]] | capture it at three points: the seal, the trajectory and the analysis |

**These six substantively revise `master_plan`'s §11 (sensitivity) and §S2 (prediction)** — see §8.1.

---

## 5. The non-dimensionalization convention (§S4 in detail)

### 5.1 The reference scales
| Symbol | Choice | Reason |
|---|---|---|
| length `L*` | `σ` = the representative particle **diameter** | the natural scale of the interaction range and the excluded volume. (Not the radius — the most common bug is confusing them) |
| energy `E*` | `k_BT` | BD is driven by thermal fluctuation. It leaves `T*=kT/ε` as a free dimensionless number |
| time `T*` | `τ_B = σ²/D₀` | the time to diffuse its own diameter. `D₀ = k_BT/γ₀`, `γ₀ = 6πη(σ/2)` |

Under this choice the HOOMD input is fixed at **`σ*=1, kT*=1, γ*=1 ⟹ D₀*=1, τ_B*=1`**.
So **1 HOOMD time unit = 1 τ_B**. The physical time conversion is `t_si = t_star · τ_B,si`.

> ⚠ Do not confuse this with the default time unit of the HOOMD documentation, `τ_LJ = σ√(m/ε)`.
> BD being overdamped, the mass does not enter the dynamics, so `m*=1` is set and the timescale is fixed by `γ`.
> This convention broken, every time is quietly wrong. `tests/test_nondim.py` watches it.

### 5.2 The main conversions
```
σ_si   = 2 a_si                          the representative diameter
γ₀_si  = 6 π η_si a_si                   Stokes drag (a sphere in an infinite medium)
D₀_si  = k_B T_si / γ₀_si
τ_B_si = σ_si² / D₀_si
F* = F_si σ_si / k_BT_si                 force
k* = k_si σ_si² / k_BT_si                the spring constant
v* = v_si σ_si / D₀_si                   velocity (= Pe)
D_r* = D_r,si · τ_B_si                   rotational diffusion
ω* = ω_si · τ_B_si                       angular frequency
G* = G_si σ_si³ / k_BT_si                the modulus
```

### 5.3 The list of dimensionless numbers (`simbot.nondim.groups`)
| Symbol | Definition | Physical meaning | The regime boundary |
|---|---|---|---|
| φ | `N v_p / V` | the volume (area) fraction | dilute <0.05, dense >0.3, RCP 0.64 |
| `Pe_F` | `F σ / k_BT` | the driving force vs thermal motion | the transition is at ~1 |
| `Pe_a` | `v₀ / (σ D_r)` | the active persistence | MIPS ≳ 40 (2D) |
| `T*` | `k_BT / ε` | thermal vs the attraction | aggregation ≲ 0.3 |
| `κσ` | screening | the electric double-layer thickness | long-range <1, short-range >5 |
| `k*` | `k σ²/k_BT` | the confinement strength | strong confinement ≫1 |
| `σ/ℓ_g` | the gravitational Pe | sedimentation vs diffusion | the transition is at ~1 |
| `Re` | `ρ v a / η` | inertia (verifying the BD premise) | must be ≪1 |
| `τ_i/τ_B` | `m/(γ τ_B)` | the validity of overdamped | must be ≪1 |

### 5.4 The dt selection rule (`simbot.nondim.choose_dt`) — implemented

The BD update: `Δr = (F/γ) Δt + √(2 D₀ Δt) ξ`

**★ The constraints are computed in SI and divided by the card's time scale at the end.** Written directly in
dimensionless units, you would have to work out "which time is the `*` of `Δt*`" every time, and in a system where
`τ_D` and `τ_trap` differ by 240 thousand times that mistake passes quietly.

The **minimum** of the constraints is adopted, and which constraint dominated is recorded:

| # | Constraint | Expression (SI) | The default target | When it turns on |
|---|---|---|---|---|
| 1 | thermal displacement | `√(2 D₀ Δt) ≤ δ_th σ` | `δ_th = 0.03` (per component) | when there is **a partner to overlap with ∪ a partner to bond with** |
| 2 | force displacement | `max\|F\| Δt/γ ≤ δ_F σ` | `δ_F = 0.005` | as above + `max\|F\|` measured |
| 3 | **stiffness stability** | `Δt ≤ s · 2γ/λ_max` | `s = 0.2` | **when there are bonds or angles** |
| 4 | the shortest relaxation time | `Δt ≤ ζ · min(τ_trap, 1/D_r, …)` | `ζ = 0.01` | when there is confinement or activity |
| 5 | active displacement | `v₀ Δt ≤ δ_a σ` | `δ_a = 0.01` | when there is active driving |
| 6 | **the accuracy target** | `Δt ≤ 2b/(1+b) · τ_trap` | the target bias `b` | a harmonic trap + `b` stated |

`λ_max = 4 k_bond + 16 k_angle/b²` (derived by `derive()`). `1·2·5` are **accuracy**, `3` is
**stability**, and `4·6` are **observability** — the three kinds cannot be merged into one.

**If not a single constraint is active, an exception is thrown.** A run going out on the defaults with no grounds is the worst case.

**The gate expressions are owned by `nondim.dt_max_{thermal,force,active,stability}`** — they are unit-independent
functions, so `choose_dt` calls them in SI and a script running in reduced units calls the same functions with `σ=γ=D₀=1`.
The thresholds all come from `config/run_policy.yaml` §timestep. If a script rewrites the numbers, fixing the policy
does not follow through — on 2026-07-28 `scripts/chain_bend.py` actually did that.

#### ★ The displacement gate is not universal (measured 2026-07-28)

Putting the four constraints on the same axis for the first hand-drawing system:

| Constraint | The `Δt*` ceiling (in units of `τ_trap`) |
|---|---|
| thermal displacement | **`108.6`** |
| the relaxation time | **`0.01`** ← dominant |

**A ratio of 1086×.** Turn on the displacement gate alone and it blocks nothing — because the reference length is `σ`
while the distance the particle explores is `ℓ_trap = σ/491`.
The 3 measurements in
[`dt-gate-should-be-displacement-based`](../../knowledge/wiki/findings/dt-gate-should-be-displacement-based.md)
were **all pair-interaction systems**, and that is the range of applicability of that conclusion.
The full text: [`displacement-gate-is-1000x-loose-for-traps`](../../knowledge/wiki/findings/displacement-gate-is-1000x-loose-for-traps.md)

⇒ The two gates do not compete but complement. **Which one turns on is decided by the card.**

#### ★ The displacement gate is only an accuracy gate (measured 2026-07-28, a bonded system)

Applying the displacement gate as it stands to a straight chain selected `Δt* = 4.5×10⁻⁴` and the chain blew up.
The gate had been **neutered** — a perfectly straight chain is a **stationary point** with `max|F*| = 0`, so
the force gate gives no ceiling. It blows up at `kT = 0` too, so it is not a stochastic phenomenon either.

| `k_bond*` | `k_angle*` | `Δt_crit` measured (bisection) | the `2/λ_max*` lower bound | measured/bound |
|---|---|---|---|---|
| `10⁶` | `10⁴` | `1.00e-6` (N=5) · `5.87e-7` (N=9) | `4.81e-7` | `2.08` · **`1.22`** |
| `10³` | `10³` | `1.84e-4` (N=5) · `1.48e-4` (N=9) | `1.00e-4` | `1.84` · `1.48` |

In all 10 combinations **the bound is below the measurement** (ratios `1.22–2.80`) → it can be used as a gate.
`s = 0.2` is a `6–14`× margin against the worst case, `1.22`.
The full text: [`dt-gate-needs-a-stability-term-for-stiff-bonds`](../../knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md)

**`hard_floor` applies to the accuracy constraints only.** The stability ceiling of `k_bond* = 10⁶`, `9.6×10⁻⁸`, is
below the floor `10⁻⁷`, and yet the measurement showed that system stable out to `Δt* = 10⁻⁶` — the floor was the wrong one.
A stability ceiling is not negotiable (it is not a demand to lower it but a demand to lower `k_bond`),
so rather than rejecting, it is left in `logged["dt_star_below_hard_floor"]`.
⇒ Unresolved: `hard_floor` is in **the card's reference time unit**, so its meaning differs per card. It is the
"universal convention" trap already discarded for the `dt/τ_D` gate, and if it is a cost check then estimating the `Λ` budget is what should do it.

`max|F|` is obtained by computing the actual forces in the initial arrangement (**estimating is forbidden**). **"I measured
and it was 0" and "I have not measured yet" are not written as the same sentence** — the former is physics (a stationary point),
the latter a procedural violation, and `DtChoice.table()` displays them differently.
`Δt/τ_D` is **recorded only** in `logged`. It is not used as a gate.

---

## 6. The `knowledge/` accumulation schema

> ## ⚠️ This section is obsolete — §10 is authoritative
>
> §6.1~6.4 below is the schema designed on the morning of 2026-07-28. On the same day, while porting
> `BD_agent/knowledge/`, **that contract was judged better and adopted as authoritative** (§10.2).
> The actual directories are `source/papers` + `wiki/{systems,findings,concepts,techniques,
> benchmarks,questions}`, and the contract is owned by
> [`knowledge/wiki/CLAUDE.md`](../../knowledge/wiki/CLAUDE.md).
>
> **Why it is kept rather than deleted:** §10.4 cites "the case where the port refuted what had been written today",
> and if the refuted original disappears that record becomes unreadable.
> **Do not read the below as implementation guidance.**

### 6.1 The structure (obsolete — see §10)
```
knowledge/
├── INDEX.md              the whole list + a one-line summary (the agent reads it first every time)
├── systems/              system archetypes: what was simulated
├── parameters/           the source and grounds of parameter values (η, T, ε, dt, r_cut…)
├── models/               the grounds for modelling decisions (why WCA, why ignore HI)
├── failures/             failure cases: symptom → cause → prescription
├── validation/           verification benchmark values (regression candidates for pytest)
└── env_log.md            the package accumulation history
```

### 6.2 The entry format (common to every file)
```markdown
---
id: dense-diffusion-hardsphere
kind: validation            # systems | parameters | models | failures | validation
tags: [diffusion, hard-sphere, dense, phi]
created: 2026-07-28
updated: 2026-07-28
runs: [2026-07-28T14-30_free-diff_a1b2c3]   # the run that produced this knowledge
confidence: medium          # high | medium | low
supersedes: []              # the id of the earlier entry this replaced
---

## Summary
One paragraph. What this entry claims.

## Grounds
Data/literature/expressions. Numbers with their errors.

## Range of applicability / limits
When this knowledge may be trusted and when not.

## References
Literature, run links, related entries [[another-id]]
```

### 6.3 The `failures/`-specific format (the most important)
```markdown
## Symptom
What was observed. The error message, the shape of the graph, the strange number.

## The diagnostic path
What was suspected and how it was excluded. (Meet the same symptom next time and follow this order)

## The root cause
The class: numerical | modeling | interpretation | analysis | environment

## The prescription
The concrete fix. A code/configuration diff.

## Preventing recurrence
The test / guard / documentation added. If none, write "none".
```

### 6.4 The knowledge accumulation rules
- **The pipeline cannot terminate at S8 without a knowledge update.** (a gate)
- If a result contradicts an existing entry, make a new entry and link it with `supersedes`. **Do not quietly overwrite the existing entry.**
- Among the `validation/` entries, the reproducible numbers get promoted into `tests/test_knowledge_regression.py` → so that a code change automatically shows whether the knowledge breaks.
- `INDEX.md` is updated automatically when an entry is added (`simbot.io.reindex_knowledge`).

---

## 7. The environment and package accumulation

### 7.1 env: `simulation_bot` (newly created)
The existing `hoomd_slit` is not touched. It starts as a new env, adds only what is needed in stages, and
**records "why it was needed" in `knowledge/env_log.md` at every addition**.

| Stage | Packages | Purpose | When |
|---|---|---|---|
| **1. the core** | `python=3.12`, `hoomd`, `gsd`, `numpy`, `scipy` | BD execution, trajectory I/O, numerics | now |
| **2. analysis and plots** | `matplotlib`, `freud`, `pandas`, `h5py` | RDF/S(q)/clusters, plots, logs | at the start of S6/S7 |
| **3. input handling** | `pillow` | hand-drawing image metadata and preprocessing | at the start of S1 |
| **4. development** | `pytest`, `pyyaml` | tests, spec serialization | when the core is implemented |
| **5. 3D rendering** | `fresnel` | ray tracing 3D snapshots | when a 3D case appears |
| **6. video** | `ffmpeg` | MP4 animation | when GIF is not enough |

`environment.yml` is updated as the stages rise, and pinned with `conda env export --from-history`.
### 7.2 The execution convention
- Every Python execution uses the absolute path of the `simulation_bot` env's interpreter
  (`conda activate` is unreliable in a non-interactive shell)
- The interpreter path is recorded in `CLAUDE.md` so the agent can refer to it every time

---

### 7.3 Local CPU parallel computing — a measured review

> Measured 2026-07-28 · machine Apple M4 (`Mac16,12`), **4 P-cores + 6 E-cores = 10 cores**, 16 GB
> The benchmark kernel: 3D WCA + `md.methods.Brownian`, φ=0.30, `dt=1e-4`, a Cell nlist (buffer 0.3)
> The raw data: `knowledge/parameters/local_cpu_parallelism.md`

#### 7.3.1 The summary of conclusions (read this first)

| Finding | Grounds |
|---|---|
| **HOOMD 7.1 is completely single-threaded in this environment.** | measured `CPU util = 1.00x` (for all of N=500~32000) |
| **There is no CPU parallelization route at all.** HOOMD v3+ removed TBB, and the remaining route is MPI domain decomposition, but — | `hoomd.version.mpi_enabled = False` |
| **conda-forge has no MPI build for **either** osx-arm64 or linux-64.** | searching the whole build string gives 0 `mpi` matches. All are `cpu*`/`gpu*` |
| ⇒ **Strong scaling (one simulation across several cores) is impossible.** If MPI is needed, a source build is mandatory. | |
| ⇒ **Throughput parallelization (running independent runs at once) is the only route, and it actually suits this project better.** | §7.4 |

#### 7.3.2 Single-process throughput — almost independent of N

| N | TPS (steps/s) | particle-steps/s | CPU util |
|---|---|---|---|
| 500 | 13 313 | 6.66e6 | 1.00x |
| 2 000 | 3 210 | 6.42e6 | 1.00x |
| 4 000 | 1 629 | 6.52e6 | 1.00x |
| 8 000 | 785 | 6.28e6 | 1.00x |
| 32 000 | 189 | 6.05e6 | 1.00x |

**The throughput constant `Λ ≈ 6.3 × 10⁶ particle-steps/s` (1 P-core).**
Raising N by 64× degrades it by only 9 % → the Cell list works well and the working set fits in cache.
With this one constant the wall time of every run can be predicted: `wall ≈ N × steps / Λ`.

#### 7.3.3 Concurrent-execution scaling (N=4000, 12 000 steps, synchronized start)

| k | total TPS | speedup | efficiency per core | TPS per process | the gain per additional one |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 629 | 1.00x | 100 % | 1 629 | — |
| 2 | 3 090 | 1.90x | 95 % | 1 545 | +1 461 |
| 3 | 4 524 | 2.78x | 93 % | 1 508 | +1 434 |
| **4** | **6 027** | **3.70x** | **93 %** | 1 507 | +1 503 |
| 5 | 6 304 | 3.87x | 77 % | 1 261 | +277 ← **the cliff** |
| 6 | 6 804 | 4.18x | 70 % | 1 134 | +500 |
| **8** | **8 041** | **4.94x** | 62 % | 1 005 | +309 each |
| **10** | **8 906** | **5.47x** | 55 % | 891 | +432 each |
| 12 | 8 663 | 5.32x | 44 % | 722 | **−243 (a regression)** |

**How to read it**
1. **k ≤ 4 is nearly linear (93 %).** The stretch where one lands on each of the 4 P-cores.
2. **The cliff at k=5** — from the fifth it starts sharing P-cores. Given that the per-process TPS scatter is
   only `1.02x`, macOS does **not pin** processes to P/E but **time-shares and migrates** them.
   So there is no "one slow run delaying the batch" problem — a property favourable to batch scheduling.
3. **An E-core's effective performance is about 1/3 of a P-core's.**
   4 P × 1 507 = 6 027, and k=10 totals 8 906 ⇒ the 6 E contribute 2 879 ⇒ 480 TPS each ⇒ `480/1507 = 32 %`.
4. **k=12 is slower than k=10 (−2.7 %).** Oversubscription is a net loss. **The ceiling is the core count, 10.**

#### 7.3.4 The recommended concurrency

| Scenario | **k** | total TPS | vs the maximum | efficiency per core | Why adopted |
|---|:---:|---:|---:|---:|---|
| interactive (the user is working at the same time) | **4** | 6 027 | 68 % | 93 % | P-cores only. Keeps the machine responsive, the highest per-core efficiency |
| **the default (recommended)** | **8** | 8 041 | **90 %** | 62 % | gets 90 % of the maximum throughput while **leaving 2 cores spare** — the agent itself, freud/numpy analysis (Accelerate BLAS is multithreaded), plotting, the OS |
| batch (unattended) | **10** | 8 906 | 100 % | 55 % | the absolute maximum. Only when there is nothing else |
| — | 12+ | 8 663 | 97 % | 44 % | **a measured regression. Forbidden** |

**Why the default is 8 and not 10**: the extra gain from using 10 is +11 %, and the price is that the cores on which
the analysis, the plotting and the agent itself would run disappear. The pipeline does not only run simulations but
keeps the S6 and S7 analysis going alongside, so 8 has the higher effective throughput.

---

### 7.4 The optimal simulation conditions (a proposal) — a human may override

> The below are **values derived from** the §7.3 measurements. They are not settled values.
> Where a human changes them: **one place, the `overrides:` block of `config/run_policy.yaml`.**
> The agent presents the condition table at the start of every run and gets approval, but if there are `overrides` it takes those first.

#### 7.4.1 The core policy: **error bars are free**

Since HOOMD is single-threaded and there are 10 cores, **the cost of running 4 independent seeds at once is the same as 1**
(k≤4 is 93 % efficient). Therefore:

> **The default policy: every production run executes ≥ 4 seeds concurrently.**
> A single run with 1 seed is forbidden (because it becomes a result with no error bar).

This is the most important design conclusion on this hardware. Four short runs are almost always better than one long run.

#### 7.4.2 The execution tiers

`Λ = 6.3e6 particle-steps/s`, `η_k` = the efficiency per process (§7.3.3), referenced to `dt* = 5e-5 τ_B`.
`wall ≈ N × steps / (Λ · η_k)` , `t_total = steps × dt*`

| Tier | Purpose | N | steps | `t_total` | concurrent k | wall per process | total batch wall |
|---|---|---:|---:|---:|:---:|---:|---:|
| **T0** `smoke` | confirm only that the code runs | 200 | 2e3 | 0.1 τ_B | 1 | **0.06 s** | instant |
| **T1** `pilot` | verify the guards and dt + predict the wall (S5-2 **mandatory**) | the same as production | 0.5 % of production | — | 1 | **≤30 s** | ≤30 s |
| **T2** `explore` | regime exploration / a parameter sweep | 1 000 | 4e5 | 20 τ_B | **8** | 103 s | **an 8-point sweep in 1.7 min** |
| **T3** `production` | the main measurement (the default tier) | 4 000 | 1e6 | 50 τ_B | **4** | 11.4 min | **4 seeds in 11.4 min** |
| **T4** `long` | GSER · MIPS · a dense system at long times | 4 000 | 1e7 | 500 τ_B | 2–4 | 1.9 h | **user approval mandatory** |

**The tier selection rule (applied automatically by the agent)**
- A system seen for the first time → always the order `T0 → T1 → T2`. T3 after confirming the regime from the T2 result.
- The `t_total` requirement is set by the quantity being measured:
  - the diffusion coefficient `D` → `t_total ≳ 10 τ_B` (T2 is sufficient too)
  - a dense system's long-time `D_L` → `t_total ≳ 100 τ_B` (the T3 floor)
  - **GSER `G'(ω), G''(ω)` → the MSD needs 3~4 decades → `t_total ≳ 100 τ_B` (T3 is marginal, T4 recommended)**
  - **MIPS coarsening → `t_total ~ 10³ τ_B` + `N ≳ 10⁴` → T4 only**
- `wall_time_budget_s = 600` (default). If it is expected to be exceeded, **do not run and report to the user**.

#### 7.4.3 The N selection logic (a proposal)
- The statistical error is `1/√(N × the number of independent time origins)`. **In a dilute system the particles are independent of each other**,
  so raising N and lengthening the time are statistically equivalent → **and since the wall time is proportional to `N × steps`, it is indifferent.**
  → In that case **raising the number of seeds is the only real gain** (it also exposes systematic error).
- Structure and phase behaviour are different: `L > 2 × the correlation length` is a hard condition.
  Since `L = (N v_p/φ)^{1/3}`, a lower φ needs fewer N for the same L.
- An active system requires `L > 2 × v₀τ_r` (the persistence length).

#### 7.4.4 The points of human intervention (only 3)

| When | Where | What |
|---|---|---|
| a permanent change | `config/run_policy.yaml` → `overrides:` | the tier defaults, k, the wall budget. Write the reason alongside and S8 reflects it in `knowledge/` |
| for one run only | answer the condition table the agent presents at the end of S4 | "T2 instead of T3", "make it N=8000" |
| during execution | when a budget-exceeded report arrives | approve / reduce / halt |

Otherwise the agent **decides for itself by the rules above and records the grounds** (the question budget, §9-5).

---

## 8. The development roadmap

| | Milestone | Content | Definition of done (DoD) | Actual |
|---|---|---|---|---|
| `[O]` | **M0** the skeleton | env, directories, the 3 documents, the knowledge schema | the initial `pytest` passing, the documents existing | ✅ |
| `[O]` | **M1+M2** the vertical slice | S1→S8 **all the way** through | `REPORT.md` generated + verified against an analytic solution | ✅ **straight to the hand drawing.** A harmonic trap instead of free diffusion — the drawing was a trap. 7 PASS / 2 INCONCLUSIVE / 0 FAIL |
| `[O]` | **M2.5** the deterministic core | making the run reproducible **in code** | the full run with one line, `cli.py run <spec>` | ✅ **2.3 s, bit-identical to the first run.** 373 tests |
| `[O]` | **M2.7** the L1 agent layer | `.claude/skills/` + `agents/` (§12.3–12.4) | a new hand drawing → the skill delegates S1→S8 | ✅ 3 skills · 5 reference documents · 9 subagents. **The delegation has not been tested yet** (a second drawing is needed) |
| `[~]` | **M3** domain extension | A→D in sequence, each with an `examples/` reference case | the per-domain verification criterion PASSing | A (the trap) only. **Trap+WCA is on hold at the user's request** |
| `[~]` | **M4** strengthening verification | the dt ladder, finite size, bootstrap, promotion to regression | the `INCONCLUSIVE` verdict logic working | ✅ INCONCLUSIVE, the power and `converge` work. Remaining: bootstrap, finite size |
| `[X]` | **M5** input extension | experimental video → paper PDFs → voice → text | reproducing M1 with each modality | hand drawings only |

**Currently: M2.7 done → next is §8.1.**

### Why M1 went differently from this design (a record)

The DoD was "run through with free diffusion first, `D* = 1.00 ± 0.03`". In reality **the user's drawing was
optical tweezers, and free diffusion is not observed in that system** (with `τ_D/τ_trap = 2.4e5`, the
free-diffusion stretch is below `10⁻⁵ τ_trap`).

Building a free-diffusion case first and then moving to the hand drawing would have made M1 **work that gets thrown away**.
Instead it went straight to the trap card, and as a result the `(system × target dynamics)` card system was
justified by measurement — it confirmed that a universal `τ_D` convention produces `Δt = 12 τ_trap` in this system.

⇒ **The free-diffusion regression case (`D* = 1.00 ± 0.03`) does not exist yet.** §8.1-①.

---

## 8.1 What to do next — **updated 2026-07-30**

> The criterion for judging: **what is blocking something else**, and **what could be quietly wrong.**

### ① Bracketing the phase boundary by melting from the crystal (S30) — **the physics top priority** (~1 hour)

`soft-r3-hexwin` read all three points inside the Zahn window as isotropic liquid, but **the boundary cannot be
found from a random start** — because of the nucleation barrier near a first-order transition it stays
arbitrarily long supercooled. Looking like a steady state does not help the distinction either.

⇒ **Start from a hexagonal crystal and see whether it melts.** Random brackets from above and the crystal from below,
and the truth is between. A hexagonally commensurate box has to be used, and since `ψ₆` is independent of the box
shape the comparison is valid ([[box-shape-confounds-initial-condition-comparison]]).
**Without this, "there is no hexatic" cannot be claimed.**

### ② Re-measuring `A = 13.3` (S31) — the bracket's lower bound (~15 min)

The `A=13.3` crystal verdict in §8.7 was measured at `βU(r_cut) = 0.24 kT`. Since the truncation error moved the
exponent by `2.9σ` (§8.8), **whether it crystallizes may be affected too.** Re-measuring at `r_cut = 7.80`
narrows the bracket considerably from `10.8–31.6` (193 %).

### ③ Fixing the cost model — **the policy tells a falsehood** (30 min)

`hexwin` was estimated at 22 minutes vs **54 minutes measured**. The overhead factor `3.4` was obtained at `N=400`,
and the frame-extraction cost grows with `N`. `estimate_wall_time_s` is **2.5× optimistic** for
large-`N` batches — the budget gate then passes a run it should not pass.

### ④ Revising the §11 sensitivity analysis with the campaign results (1 hour)

The 6 methodology findings of §4.1 substantively revise §11. In particular:
- put **the pre-computed design power** (`seeds_for_target_sigma`) into the §11.5 verdict rules
- enforce **the `t(ν)` correction** in the tolerance derivation (currently it is the normal quantile)
- add **the sequential-design pre-registration** convention to §S2 (nail `no_stage_N` down in advance)

### ⑤ The free-diffusion regression case (Q9) — **still empty** (half a day)

The `brownian` scale route has still never been executed end to end. The soft-repulsive system used a
separate route, `Soft2DRunConfig`, so **this debt has not been paid.**

### ⑥ On hold — awaiting the user's judgment

| | Why it is on hold |
|---|---|
| the `τ_relax` stage 3 (`k = 2145`) | `no_stage_3: true` was nailed into the pre-registration. Running it requires **a human decision** |
| the `A = 100` finite-size ladder | **impossible on this machine** — a single `N=1024` run is 3.5× the budget |
| trap + WCA (M3) | on hold at the user's request (2026-07-28). The queue's `trap-drag-2d-hex300` requires it |
| the 2 remaining hand drawings in the queue | `abp-rod-2d-run-flip` · `trap-drag-2d-hex300` |

---

## 8.2 ⚠️ Obsolete — the queue reading as of 2026-07-28 (kept as a record)

> The below is the triage from when 4 hand drawings came into the queue. **`soft-r3-2d-A-sweep` was run to the end and
> `chain-bend-2d-oscill` is closed too** (before commit `fefd5c9`). Two remain.
> The reason the original is not deleted is that "what was judged to be the bottleneck, and was that right" has to
> be readable — **the judgment that the pair runner was the bottleneck was right.**

### (old) What to do next — the priorities and the grounds (as of 2026-07-28)
> The criterion for judging: **what is blocking something else**, and **what could be quietly wrong.**

### ★ The 4 hand drawings that came into the queue (2026-07-28, added by the user)

4 new hand drawings and 2 papers came into `inputs/`. **None of the four runs with the current code.**
What each requires (a reading triage — the formal S1 is done at the start):

| Drawing | System | Domain | The stated values | What is needed | Why it is blocked |
|---|---|---|---|---|---|
| **`soft-r3-2d-A-sweep`** | `U/kT = A/r³` soft repulsion, a 2D square | **A/D** structure | `A = 0.1, 1, 10, 100` · `N=100` · `Lx=Ly` | an `r⁻³` pair potential · a **sweep** · RDF · Voronoi | no pair runner · `sweep:` unsupported · no `analysis/structure.py` |
| **`trap-drag-2d-hex300`** | dragging 1 probe through a 2D hexagonal equilibrium | **B** active microrheology | `N≈300` · `R=5 μm` · `k_t=10 pN/μm` · `v=0.5 μm/s` | pair interactions · **a moving trap** `r_trap(t)=r₀+vt` · a drag measurement | the above + `HarmonicTrap`'s centre is fixed |
| **`chain-bend-2d-oscill`** | bending a chain held by optical tweezers as `y=a sin(ωt)` | **B** active microrheology | `k_t=10 pN/μm` · `R=5 μm` | **a bond potential** · an oscillating trap · GSER → `G'(ω)`, `G''(ω)` | ★ **`U_ij` is blank in the drawing** ("see the Eric Furst paper") · no `analysis/microrheo.py` |
| **`abp-rod-2d-run-flip`** | a single active ellipsoid, run-and-flip | **C** active | `τ_R=0.5 s` · `v ≤ 5 μm/s` | an anisotropic particle · the orientational degree of freedom · a 180° flip · **MSAD** | the `abp` card is `draft` · the `active_run_length` scale is `NotImplementedError` · anisotropic drag |

**The common denominator settles the answer:**
- **3 of the 4 require pair interactions** → the pair runner is the top-priority bottleneck
- **2 of the 4 require a time-dependent trap** (dragging, oscillation) → `center(t)` on `HarmonicTrap`
- The 2 papers (`PhysRevLett.94.138301`, `la7023617`) appear to have been added to fill in `chain-bend`'s
  blank `U_ij`. Installing `pypdf` is for the same reason (env-log stage 3)

⇒ **The trap+WCA the user had put on hold is now on the queue's critical path** (`trap-drag-2d-hex300`).

### ⓪ The pair-interaction runner + `soft-r3-2d-A-sweep` — **the queue's bottleneck** (1 day)

★ **Raised to top priority after reading the queue.** It is the point where 3 of the 4 are blocked, and
`soft-r3` is **the cleanest entrance** among them — pure pair interactions with no trap and no driving.

What this one thing solves at once:

| | How |
|---|---|
| the pair runner (the bottleneck common to 3 drawings) | add `RUNNERS["...equilibrium-structure"]` to `run.py` → **the dispatch design gets tested for the first time** |
| the `dt` **displacement gate** | once there are pair interactions the gate turns on → **it constrains `dt` for the first time** (half of ①'s concern is solved here) |
| the `brownian` scale route (`σ`, `τ_D`) | this system's card uses it → **executed end to end for the first time** (the other half of ①) |
| `sweep:` support (§9-7) | `A = 0.1, 1, 10, 100` is **stated** in the drawing. All four have to be run for an answer |
| `analysis/structure.py` | RDF · Voronoi · ψ₆. **A verification target exists for the first time** (`freud` is installed) |
| the `max\|F\|` measurement route | the force-displacement constraint escapes `n/a` (§5.4 "estimating is forbidden") |

**It is verifiable:** at `A=100` it should crystallize hexagonally, and at `A=0.1` it is nearly an ideal gas.
2D melting has a distillation in `knowledge/source/papers/1999-zahn-two-stage-melting-2d.md`.
With `N=100` the cost is low too.

⚠ `r⁻³` is not a built-in HOOMD potential → confirming first which of `md.pair.Mie(n=3, m=0)` or `pair.Table`
is correct comes first.

### ① The free-diffusion regression case — **a verification debt** (half a day)

> ⓪ absorbs a substantial part of this item (the `brownian` scales · the displacement gate · the runner branch).
> Even so, **the exact-solution comparison itself**, `D* = 1.00 ± 0.03`, remains separate — with pair interactions
> `D` is no longer 1 and so a pure free-diffusion case is needed separately. It gets attached as ⓪'s smoke run.

`D* = 1.00 ± 0.03`. It was M1's original DoD and it does not exist yet. Because it does not, **three things
remain unverified:**

| Unverified | Why it is dangerous |
|---|---|
| the `brownian` route of `CARD_SCALE_RULES` (`σ`, `τ_D`) | **it has never once been executed end to end.** The card system's central claim ("the scales differ per system") is only half tested |
| the `dt` **displacement gate** | it is off in a trap system → **0 runs where it actually constrained `dt`.** There are unit tests only |
| the `RUNNERS` branch | with only 1 entry, the runner dispatch design is untested |

What is needed: a free-BD runner in `run.py` (no trap force) + registering the `passive-sphere--transport`
card + `analysis/msd.py` (a **verification target** exists here for the first time).

**The analytic solution is exact (`D* = 1`) and the cost is low.** If it is wrong it shows up immediately.

### ② A hand-drawing guide — **the cheapest leverage** (30 min)

`docs/drawing_guide.md`. The first drawing produced 2 ambiguities and **both could have been eliminated by this guide:**

| Ambiguity | What the guide would require |
|---|---|
| A1 (2D vs a 3D cross-section) | "please write the axes and the dimensionality" |
| A2 (`a` vs `R` — is it the radius or the diameter) | "on a number, the unit and **what size it is**" |

A2 **could not be closed by measurement because there is no experimental `f_c`** — it was closed by parsimony
with only the refutation condition recorded. One more character on the drawing and the ambiguity would never have arisen.

It raises ③'s success rate directly.

### ③ The remaining 3 drawings in the queue — the order and the reasons

Once ⓪ is done the pair runner exists, so the following order is natural:

| Order | Drawing | Reused from ⓪ | What is newly needed |
|---|---|---|---|
| 1 | **`trap-drag-2d-hex300`** | the pair runner · RDF | a time-dependent trap centre · a drag measurement. **The trap+WCA the user had put on hold acquires a purpose here** |
| 2 | **`chain-bend-2d-oscill`** | the trap (from ⇧) | a bond potential · GSER. **The 2 papers have to be read first to fill in the `U_ij` blank** (`bd-lit-distill`) |
| 3 | **`abp-rod-2d-run-flip`** | — | an anisotropic particle · the orientational degree of freedom · MSAD. **The largest new physics** |

What each drawing tests (apart from the physics):
- **the S1 reading protocol** — do the rules induced from the first drawing hold for a new drawing
- **skill delegation** — does `.claude/` actually work (Q10; right now it is **untested**)
- **handling a pair with no card** — does `nondim` throw and make you create the card first
  (all 4 have no card → this route fires 4 times)

### ④ The pilot run — **a mismatch between the policy and the code** (1 hour)

`run_policy.yaml` declares `pilot: {mandatory: true}` and `cli.py` does not run it. **Right now the policy file
tells a falsehood.** The trap system had a wall of 2 s so there was no need, and so it was deferred, but the
moment runs become minutes long (trap+WCA, `N=8000`) it acquires meaning.

The minimum action: implement it, or set `mandatory: false` in the policy and write the reason.

### On hold — awaiting the user's judgment

| | Why it is on hold |
|---|---|
| trap + WCA (M3) | **on hold at the user's request** (2026-07-28). The point at which `kT_conf` becomes an independent check and `r_cut` acquires meaning |
| the 5 `analysis/` modules | there is neither a caller nor a verification run. ① opens `msd.py` |
| comparison against experimental data | needed to close A2 by measurement. There are no measured values |

---

## 9. Additional proposals (miscellaneous)

Items proposed at the user's request. Whether they are adopted is confirmed separately.

### Strongly recommended
1. `[O]` **git init** — ✅ done. 28 commits, since `4ac2a53`. A repo independent of BD_agent.
2. `[O]` **prediction sealing** — ✅ `simbot/io.py`.
   `SEALED.sha256` is in the standard `sha256sum` format, so **it verifies with `shasum -c` without our code** —
   because the trustworthiness of the seal must not depend on our code.
   If the seal breaks, `validate_run` **does not make the comparison table.**
   ★ Not in the design: documents made **after** the seal are reported separately as `unsealed`.
   Because "the sealed files passed" does not mean "the prediction was sealed".
3. `[X]` **the pilot run** — unimplemented. It is declared in `run_policy.yaml` as `mandatory: true`
   but `cli.py` does not run it yet. **In a trap system the wall was 2 s so there was
   no need, and so it was deferred.** It becomes necessary for the first time with trap+WCA.
4. `[O]` **enforcing the unit suffix** (`_si` / `_star`) — ✅ `Quantity.si` rejects strings
   and bools, and the round-trip test catches scale confusion (`test_s4_nondim.py`).
5. `[~]` **the question budget** — the convention is in `CLAUDE.md` and **a human keeps it.** There is no
   enforcement device in the code (with no agent layer there is not yet a point at which to enforce it either). In M2.7 the skill owns it.

### Recommended
6. `[O]` **the `bd-diagnose` skill** — ✅ `.claude/skills/bd-diagnose/SKILL.md`.
   It holds the exclusion order (statistic fluctuation → self-consistency → sample independence → units → numerics → **then physics**).
   The grounds are the measurement that 0 of the first full run's 4 problems were physics problems.
7. `[X]` **parameter sweep support** — `sweep: [...]` is unsupported in the `spec`.
   `cli.py converge` shaking dt, N and the seed is **a convergence check** and not a regime map.
8. `[~]` **a cache** — the material exists (`spec.hash()`, `RunDir.completed_stages()`,
   `cli.py resume` reusing a completed run). **There is no lookup that finds an existing run
   by `spec_hash`** — `resume` has to be handed the directory directly.
9. `[X]` **`REPORT.md` → HTML** — unimplemented. **Now that there are 5 figures there is a benefit**
   (currently the `figs/*.png` links are relative, so moving the report alone breaks the figures).
10. `[X]` **a hand-drawing guide** (`docs/drawing_guide.md`) — unwritten.
    **The first drawing's reading produced 2 ambiguities (A1 the dimensionality, A2 `a` vs `R`) and both would
    have been absent had this guide existed.** The unresolved item with the largest benefit per cost.

### To be decided after review
11. **A Langevin fallback** — BD is unsuitable for a system where `τ_i/τ_B ≪ 1` breaks (small particles, low viscosity).
    Whether to detect it automatically and switch to `md.methods.Langevin`, or only warn and halt.
12. **Introducing an HI approximation** — out of scope for v1, but it becomes necessary if dense-system results
    systematically fail to match the literature.
    The minimal introduction: the Rotne–Prager mobility matrix + Cholesky (limited to N ≲ 10³, implemented outside HOOMD).
13. **Direct comparison against experimental data** — if the user gives an experimental MSD or video, analyse it
    through the same pipeline as the simulation and compare side by side. An extension of S7.

---

## 10. Porting the knowledge layer — absorbing `BD_agent/knowledge/`

`/Users/kyuhwan/Desktop/BD_agent/knowledge/` already has a considerably mature knowledge contract.
**It is better than the schema I newly designed in §6. It is adopted as authoritative and §6 is obsoleted.**

### 10.1 What is there
```
BD_agent/knowledge/
├── raw/lab/          the original PDFs + LaTeX source tarballs (~20 papers, gitignored)
├── source/
│   ├── papers/       ★ 42 per-paper distillation .md files + INDEX.md  ← an immediate asset
│   └── lab/          unpublished lab assets
└── wiki/
    ├── CLAUDE.md     ★ the knowledge contract (the frontmatter is the machine contract, the prose is the grounds)
    ├── systems/      ★ (system × target dynamics) cards + _TEMPLATE.md + _index.md
    ├── concepts/  techniques/  benchmarks/  findings/  questions/
```
The 42 distillations are Choi, Gubbala, Arnold, Kim, Xu, Cheon, Takatori, Barakat, Quah, Modica and others —
**the user's lab's actual paper corpus.** An asset that cannot be created from scratch.

### 10.2 The design being adopted (where it beats my §6)

| Their design | Why it is better |
|---|---|
| **the 3 layers `raw / source / wiki`** | copyright (the original), the authoritative record (the distillation) and interpretation (the synthesis) are separated. My flat structure cannot do this |
| **the `(system, target dynamics)` pair card** ★ | **it refutes §5** — see 10.3 below |
| **`reproduced: yes/no/partial`** | stronger than my `verified`. "It was published in a paper ≠ it runs in our code" |
| **the `[source]` vs `[source, unreproduced]` notation** | it forces verification grounds and factual records to be distinguished in the report |
| **splitting the folders by publication status** (`papers/` vs `lab/`) | because the publication boundary is at folder granularity. Mix lab papers in with unpublished material and they drop out together at publication |
| **`precedence L0–L3`** | "a lower L wins on what is true, and a higher L wins on what to do about it" |
| closing `questions/` with a `status` rather than deleting them | the list of unresolved problems remains |
| `dead-end-<slug>.md` | a blocked path is an asset too |

### 10.3 ★ This port **amends §5** — the non-dimensionalization is not fixed by the system alone

What their `wiki/CLAUDE.md` shows by measurement:

| System × target dynamics | Reference length | Reference time | `kT` |
|---|---|---|---|
| ABP × control | **the run length `ℓ`** | **`τ_r = 1/D_r`** | **a derived quantity** |
| brush colloids × non-equilibrium contact | `σ` | `τ_D = σ²/D` | an input |
| a passive tracer × transport | `σ` | `τ_D` | an input |

> **§5 nailed `(σ, k_BT, τ_B)` down as a universal convention. That is wrong for domain C (active matter).**
> In an ABP the natural units are the run length and the rotational relaxation time, and `kT` becomes a derived quantity.
>
> **The amendment:** §5 is demoted to **the default convention of domains A, B and D.**
> The owner of the non-dimensionalization convention is the `wiki/systems/<system>--<dynamics>.md` card.
> **On meeting a pair with no card, ad hoc non-dimensionalization is forbidden** — create a `status: draft` card from `_TEMPLATE.md` first.

The gates also differ per pair (their table):

| Gate | passive spheres × equilibrium structure | ABP × a dense collective |
|---|---|---|
| the equilibration verdict | ✅ valid | ❌ **meaningless** — an active system never reaches thermal equilibrium |
| `D_msd = kT/γ` | ✅ holds | ⚠️ does not hold — `D_eff = D_t + v₀²τ_r/2` |
| the advective displacement `v₀Δt/σ` | not applicable | ✅ mandatory |

### 10.4 The case where the port **refuted what had been written today** (already amended)

`wiki/findings/dt-gate-should-be-displacement-based.md`:

| | The value I wrote today | Their measured grounds | Action |
|---|---|---|---|
| the `dt` ceiling | `hard_ceiling: 1e-4` | it **rejects 2** of 3 runs that actually ran (the preceding slit 1e-3, the Quah code 1.67e-4) | ✅ amended to `4.5e-4` (`config/run_policy.yaml`) |
| the displacement convention | `√(2d·Δt)` (the total displacement in d dimensions) | the lab practice is `√(2Δt)` (per component) | ✅ stated as `per_component` + the conversion `×√d` recorded |
| the displacement threshold | `0.02σ` | measured 0.006 / 0.018 / 0.045σ → `0.03σ` recommended | ✅ amended to `0.03σ`, the default dt kept at 0.010σ |

**The lesson:** set a gate without the literature or a measurement and it rejects configurations that actually work.
That alone is sufficient reason to do the port.

### 10.5 The measured prior additionally gained (`findings/lab-bd-conventions.md`)

| Item | The lab's measured value | The implication for our policy |
|---|---|---|
| the engine | the BD papers are **all HOOMD-blue** (Quah: 3.8.1) | we are on 7.1.0 → **an API port is needed** |
| the execution hardware | **GPU-accelerated** in many cases. Xu 2024 does `8×10⁸` steps | not reproducible on the M4 CPU (~35 hours at `N=1000`). The §7.4 budget gate must filter it |
| excluded volume | WCA is standard. Quah is stable even at **`ε/kT = 500`** | a strong WCA is not itself dangerous. The danger is **its coupling with the displacement** |
| suppressing crystallization | Takatori: a diameter ratio of **1.4**, a mole fraction of **2/3:1/3** (`φ`≤0.83) | at high `φ` a **forced-bidispersity gate** is needed. Without it a crystal gets misjudged as a "glass" |
| statistics | Xu 2023: **an average over 20 realizations** | our T3 default of 4 seeds is weaker than the lab practice → §10.6 unresolved |
| the definition of `φ` | `nπσ²` (the diameter) vs `n̄πa²/4` — **it differs per paper** | **it can be quietly 4× wrong.** `simbot` always records the diameter/radius convention when computing φ |

### 10.6 The porting plan and what is unresolved

| | Step | Notes |
|---|---|---|
| `[O]` | port the `wiki/CLAUDE.md` contract | ✅ carried over as `knowledge/wiki/CLAUDE.md` |
| `[O]` | the 42 `source/papers/` + `INDEX.md` | ✅ copied. An immediate asset |
| `[X]` | copy the `wiki/systems/` cards + `_TEMPLATE.md` + `_index.md` | the grounds for the §5 amendment |
| `[X]` | copy `wiki/findings/` and `benchmarks/` | `benchmarks.yaml` is the grounds for the pytest regression |
| `[X]` | rewrite the 3 entries made today into the new schema | `water_298k` → `wiki/concepts/`, `no_hydrodynamics` → `wiki/concepts/`, `local_cpu_parallelism` → `wiki/techniques/` |
| `[X]` | whether to copy `raw/` needs deciding | 2.3 MB. Being a gitignore target, a reference to the original location may be sufficient |
| `[X]` | **unresolved:** the default seed count, 4 vs the lab practice of 20 | at k=8 it is 8 seeds in 1.7 min (T2). How many at T3? → needs the user's judgment |
| `[X]` | **unresolved:** write the HOOMD 3.8.1 → 7.1.0 API port table | the `gamma.default` syntax and so on |

> ⚠️ The port is **an adoption, not a copy.** Their `master_plan.md` and `docs/00_decision_log.md` also have to be
> read and the decision history (D1–D16 and so on) inherited. They have not been read yet.

---

## 11. Sensitivity analysis — S7b

> Something that was in none of the pipeline's stages. **Inserted right after the S7 verification and right before the S8 conclusion.**

### 11.1 Why it is needed — it connects directly to provenance
S1 and S3 fill values absent from the drawing with `provenance: assumed`. If the conclusion depends on those assumptions,
**it is not a conclusion but a guess.** The sensitivity analysis answers this question:

> "If a value I filled in arbitrarily were wrong, would the conclusion change?"

**The automatic connection rule:** every field with `provenance: assumed` is a sensitivity-analysis candidate.
A human does not pick the candidates — the spec designates them itself.

### 11.2 ★ The sensitivity is computed **in the space of dimensionless numbers**
Measuring the sensitivity in raw SI parameters is a waste. `η, T, a` enter only through `D₀` and `τ_B`, so
shaking the three separately is **shaking the same direction three times**.

> **The rule:** the sensitivity is computed against the dimensionless-number ledger of §5 (or the relevant systems card).
> When `m` SI parameters reduce to `n` dimensionless numbers, the number of runs drops `2m → 2n`, and usually `n ≪ m`.

### 11.3 The 4 stages (cheapest first)

| Stage | Method | Cost | When |
|---|---|---|---|
| **A. regime proximity** | how far each dimensionless number is from a regime boundary. `d = \|log(X/X_c)\|` | **0 runs** | always, **at S4** |
| **B. local first order (OAT)** | one dimensionless number at a time `×2, ÷2` → the dimensionless sensitivity index `S_i = ∂lnQ/∂lnX_i` | `2n` runs (T2) | the default. Always |
| **C. second-order interactions** | only the top 2~3 with a strong `S_i`, combined on a grid | `~9` runs (T2) | when B shows a sign of nonlinearity |
| **D. global (Sobol/LHS)** | sampling the whole assumption box | `≥64` runs | only when C shows a strong interaction. User approval |

**A needs no runs and may be the most important.** If `Pe = 45` and the MIPS boundary is `Pe_c ≈ 40–60`, then
that system is by definition maximally sensitive — and you can know it before running.

### 11.4 The cost on our hardware — effectively free
On the §7.3 measurements, 1 T2 run (`N=1000`, 4e5 steps) = 103 s, with 8 concurrent possible:

| The number of dimensionless numbers `n` | OAT runs | Batches (k=8) | **total wall** |
|---|---|---|---|
| 2 | 4 | 1 | **1.7 min** |
| 4 | 8 | 1 | **1.7 min** |
| 8 | 16 | 2 | **3.4 min** |

> **Conclusion: there is no reason to omit the sensitivity analysis.** Up to 4 dimensionless numbers it finishes at
> 1/7 the cost of a single production run. **It is enabled by default.**

### 11.5 The verdict rules

| `\|S_i\|` | Interpretation | Action |
|---|---|---|
| `> 1` | the conclusion **depends strongly** on this assumption | a warning at the top of the report. The assumption has to be narrowed (a literature lookup / a question to the user / a conditional conclusion) |
| `0.2 – 1` | a moderate dependence | stated in the report |
| `< 0.2` | **irrelevant** | **report explicitly** that "this assumption does not change the conclusion" — that is a result too |

**Whether the sign and magnitude of `S_i` agree with the prediction (S2) is also checked.** A disagreement is a signal of a defect in the model understanding.

### 11.6 The outputs
- `07b_sensitivity.md` — the sensitivity table + the regime proximity + the verdict
- `figs/tornado.png` — a tornado plot (horizontal bars in descending `|S_i|`)
- The "confidence and limits" of `08_conclusion.md` **has to cite** this result (a gate)

### 11.7 Failure modes
- computing the sensitivity in SI parameters and wasting runs on duplicate directions
- a `×2, ÷2` perturbation crossing a regime boundary so `S_i` loses meaning (→ shrink the perturbation to `±20%` and retry)
- T2's statistical error being larger than `S_i` so everything is `INCONCLUSIVE` (→ raise the seeds or promote to T3)

---

## 12. Model tiering — where to spend expensive reasoning

> The principle: **extraction is cheap, interpretation is expensive.** And **computation is code, not the LLM.**

### 12.1 The allocation

| Stage / task | Model | Grounds |
|---|---|---|
| **S1 hand-drawing interpretation** — judging the geometry, boundaries and dimensionality, the meaning of arrows, generating ambiguity candidates | **Opus 5** | multimodal + physics reasoning. **Get this wrong and everything after it is wrong.** The most expensive error point |
| S1 extraction — reading text, numbers and labels, indexing EXIF, resolution and files | **Haiku 4.5** | structured extraction. No reasoning needed |
| **S2 prediction reasoning** | **Opus 5** | a scientific claim that gets sealed |
| S3 filling in the spec (a knowledge lookup + applying rules) | **Sonnet 5** | little room for judgment |
| S3 YAML serialization and provenance tidying | **Haiku 4.5** | structured |
| S4 non-dimensionalization | **code** (`simbot.nondim`) + a Sonnet 5 review | the LLM does not make the numbers |
| S5 execution | **code only** | no LLM |
| S6 figure generation / captions | code / **Sonnet 5** | |
| **S7 the verdict + the cause hypothesis** | **Opus 5** | causal inference. Misjudging a FAIL's cause is the most expensive |
| S7b sensitivity interpretation | **Sonnet 5** | the numbers are code, the interpretation only |
| **S8 the conclusion** | **Opus 5** | the final scientific claim |
| S8 drafting the knowledge entries | **Sonnet 5** | filling in a template |
| **failure diagnosis** (`bd-diagnose`) | **Opus 5** | hypothesis generation and elimination reasoning |
| literature distillation (writing a `source/` entry) | **Sonnet 5** | in bulk. But **an equation transformation gets an Opus review** |
| bulk literature scanning, bibliographic extraction, updating the `INDEX` | **Haiku 4.5** | structured and in bulk |

### 12.2 The safeguard — keeping a cheap model from making a physics judgment

> **The rule: a field whose `provenance` is `inference` or `assumed` may only be written by Opus.**
> `observation` and `derived` may be filled by a cheap model.

Why this rule is good: it sits directly on top of the existing provenance schema, and it is **mechanically checkable**.
`simbot.spec.validate` requires a `written_by` field on each field and catches violations.

### 12.3 Implementation — ✅ done (2026-07-28)
Defined by the `model:` frontmatter of `.claude/agents/*.md`, and the `bd-pipeline` skill delegates per stage.
**Whether the allocation table and the actual `model:` agree is checked by `tests/test_agent_layer.py`.**

| Subagent | model | Responsibility |
|---|---|---|
| `bd-intake-extract` | haiku | S1 extraction |
| `bd-intake-interpret` | opus | S1 interpretation |
| `bd-predict` | opus | S2 |
| `bd-spec` | sonnet | S3 |
| `bd-validate` | opus | the S7 verdict |
| `bd-conclude` | opus | S8 |
| `bd-lit-distill` | sonnet | literature distillation |
| `bd-lit-scan` | haiku | literature scanning and indexing |
| `bd-diagnose` | opus | failure diagnosis |

> The purpose is not cost saving but **the allocation of speed and quality**.
> Running the input scan on Haiku finishes listing 10 hand drawings in a few seconds, and that saving
> can be spent on the S1 interpretation and the S7 verdict.

### 12.4 The implementation result (2026-07-28) — the reference documents number 5 and not 8

**A Q6 decision.** The design had 8 per-stage reference documents. But after the deterministic core was completed,
**S3, S4 and S5 became one line, `cli.py run`**, and putting a separate document on each would create 3 thin files
saying only "call this function".

**The documents go where the content is:**

| Document | Why it is independent | Size |
|---|---|---|
| `s1_intake_drawing.md` | **the only stage that cannot be expressed in code.** The most expensive error point | the thickest |
| `s2_prediction.md` | the seal, tolerance and power discipline. Sloppy here and the verification is neutered | |
| `s3_s5_execute.md` | all three are a `cli.py` call + reading the gate → **merged, the flow becomes visible** | |
| `s6_s7_validate.md` | the figures and the verdict serve the same judgment (what is anomalous) | |
| `s8_knowledge.md` | narrating the conclusion + the knowledge contract | |

**The skill does not rewrite the physics.** It cites the `simbot` functions and the `knowledge/wiki/` cards.
The skill layer's only unique content is **the S1 hand-drawing reading protocol**.

The structure is watched by [`tests/test_agent_layer.py`](../../tests/test_agent_layer.py) (64 tests) —
frontmatter validity, link integrity, **agreement between the §12.1 allocation table and `model:`**, the cheap models'
permission boundary being stated, and `settings.json`'s refusal to edit a sealed document.

★ **`settings.json` refuses `Edit` on a sealed document:**

```json
"deny": ["Edit(./runs/**/02_prediction.md)", "Edit(./runs/**/01_intake.md)",
         "Edit(./runs/**/SEALED.sha256)", "Bash(conda activate:*)"]
```

The prediction is **generated** by `cli.py` in Python, and what is blocked is the agent afterwards
**fixing it** by editing the text. The seal verification catches it after the fact, but **it is better to make it impossible in the first place.**

The details: [`.claude/README.md`](../../.claude/README.md)

---

## 13. Open questions

### Closed

| # | Question | The answer (2026-07-28) |
|---|---|---|
| Q1 | should git be initialized? | ✅ **done.** 28 commits (2026-07-30) |
| Q2 | what is the subject of the first hand drawing? | ✅ **2D optical tweezers** (`R=5 μm`, `k=10 pN/μm`, `T=300 K`). The text free-diffusion stage was skipped (§8) |
| Q3 | the report language? | ✅ **a Korean body + English technical terms.** But **the text in matplotlib figures is in English** (the default font has no Hangul glyphs) |
| Q4 | the allowed wall-time ceiling per run? | ✅ **10 minutes.** `budget.wall_time_per_run_s: 600`. When `cli.py` expects it to be exceeded it **does not run and halts** (fixed by a test) |
| Q5 | is 2D supported as first class? | ✅ **yes.** The first system was 2D and the `Lz=0` route is verified (`plateau = 2d`, the `⟨r²⟩` ratio `3/2` confirmed to 0.3 %) |

### Newly opened

| # | Question | Status |
|---|---|---|
| Q6 | **at what granularity should the `.claude/` skills be split** | ✅ **closed — 3 skills + 5 reference documents** (§12.4) |
| Q7 | why the trap kernel's throughput is 1.2–1.4× the baseline | open → [`questions/trap-kernel-throughput-vs-wca-baseline.md`](../../knowledge/wiki/questions/trap-kernel-throughput-vs-wca-baseline.md) |
| Q8 | **should `scripts/trap_batch.py` be deleted** — `run.run_trap_batch` does the same thing | unresolved. Kept because it is the reproduction route of the first full run |
| **Q9** | **when should the free-diffusion regression case be made** | **still open** → §8.1-⑤. The soft-repulsive system used the separate `Soft2DRunConfig` route, so this debt was **not paid** |
| Q10 | does skill delegation actually work | ✅ **closed (2026-07-29).** The `soft-r3` campaign's 6 runs were completed with the `bd-pipeline` skill |
| Q11 | should `pilot: {mandatory: true}` be implemented, or the policy fixed | unresolved. **Right now the policy file tells a falsehood** |

### Newly opened on 2026-07-30

| # | Question | Status |
|---|---|---|
| **Q12** | **can a phase boundary be found from a random start** | ❌ **no — closed.** Because of supercooling it gives only an upper bound. It has to be melted from the crystal (§8.1-①) |
| **Q13** | may the truncation-error tolerance be erected on a **value** criterion | ❌ **no — closed.** `βU(r_cut)` keeping the value within `3σ` still **biases the exponent by `2.9σ`**. A separate exponent-based tolerance is needed |
| **Q14** | the overhead factor of `estimate_wall_time_s` | open → §8.1-③. At large `N` it is **2.5× optimistic** (`hexwin` 22 min expected vs 54 min measured) |
| **Q15** | can a `χ²` form test be done with 4 seeds | ❌ **no — closed.** Since `χ² ∝ 1/SE²`, the 41 % uncertainty of a 4-seed SE is amplified 4× |
| **Q16** | is `A = 13.3` a crystal under a proper truncation too | open → §8.1-② (S31) |
| **Q17** | should `master_plan` §11 (sensitivity) and §S2 (prediction) be revised with the campaign results | open → §8.1-④. The 6 findings demand a substantive revision |

---

*This document is the single source of truth for the design. On a design change, fix this document before the code.*
*An obsoleted section is not deleted but marked `⚠️ obsolete` — if the refuted original disappears, that refutation becomes unreadable.*
