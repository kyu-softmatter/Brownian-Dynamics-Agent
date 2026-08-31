# 00 · Decision Log

> **What this document is for**
> It is the one place where everything in the agent design that "has to be decided but has not been yet" is collected, and filled in one at a time.
> So that work can proceed with decisions still deferred, every item has a **recommended default** entered in advance.
> A default can be overturned at any time — when it is, change `status` to `DECIDED` and update the documents in the `impact` cell.

**Status legend**

| Status | Meaning |
|---|---|
| `OPEN` | Not yet decided by a human. **In progress with the default applied** |
| `DECIDED` | Explicitly decided by a human. The date and the reason are recorded |
| `SUPERSEDED` | An earlier decision was overturned. The new item's number is left behind |

**Update rules**
1. When a decision changes, leave **why it changed** in `grounds`. Leave only the outcome and the same argument gets repeated in 6 months.
2. Fix the documents in the `impact` cell in the same commit.
3. When the code and this document disagree, **this document is the intent** — fix the code.
   (inheriting the rule from the preceding project's `~/Research/MD_particle/brownian_slit_sim/docs/model_assumptions.md`)

---

## A. Architecture

### D1 · What form to build the agent in
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | **Hybrid + four layers** — the deterministic core `bdkit/` + a thin LLM layer `agent/` + a behaviour-rule layer `.claude/rules/` + a knowledge layer `knowledge/`<br>※ `master_plan.md` §4 counts the same structure as **"three layers"** — because `bdkit/` and `agent/` share one contract document (the root `CLAUDE.md`). **Counted by contract document it is three; counted by code boundary it is four.** Not an inconsistency but a difference in what is being counted |
| **Options** | (a) hybrid (b) Claude Code only (c) a standalone Python package (Agent SDK) |
| **Grounds** | Keeping the core independently verifiable with `pytest`, without an LLM, is the top priority. If, when the agent is wrong, you cannot tell *whether the physics is wrong or the LLM is*, debugging is impossible. The seams for moving to (c) later (`Runner`, schema-based I/O) are left in place.<br>**Why it was settled 2026-07-27:** a **behaviour-rule layer** and a **knowledge-compounding layer** were added to the initial design (core + LLM, two layers). The former leaves "why is it done this way" outside the code, and the latter makes knowledge survive the end of a session. Without those two layers the agent starts from scratch every time. |
| **Impact** | `master_plan.md` §4, `01_agent_architecture.md`, the whole repository structure |
| **Decided** | **2026-07-27** |

### D2 · Where the simulation runs
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | **Local execution + a `Runner` abstraction to keep the seam for cluster extension.** v1 implements only `LocalRunner`; `SlurmRunner` is an interface definition only |
| **Options** | (a) the local M4 only (b) local+HPC (c) HPC-centred |
| **Grounds** | User-fixed (2026-07-27): *"the simulation and the visualization will run on my local computer, but there is a possibility of extending to a cluster later."* The only compute resource at present is the M4 (10 cores, 16GB, no CUDA). The HOOMD build is also `gpu_enabled=False`, `mpi_enabled=False`. |
| **Derived constraints** | The things that are expensive to fix later are **observed from v1 onwards** — ① no absolute paths ② `simulate.py` is self-contained ③ no hardcoded device (`make_device(spec)`) ④ resumable from a checkpoint ⑤ results communicate through files only (no stdout parsing). Details in `master_plan.md` §8 |
| **Impact** | `master_plan.md` §8, `01` (S7 EXECUTE), `05` (the cost gate), `10_roadmap.md` |
| **Decided** | **2026-07-27** |

### D3 · The scope of physics v1 will cover
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **BD + spherical colloids in 3D**. Potentials: WCA / Yukawa / Morse / DLVO. **HPMC is v2** |
| **Options** | (a) BD spherical 3D (b) BD + HPMC (c) an extension of the preceding slit project |
| **Grounds** | MC's verification logic (acceptance-ratio tuning, detailed balance, the equilibration verdict) is entirely different from BD's, so in practice there would be two pipelines. Starting narrow and completing S1–S12 once takes priority. |
| **Impact** | `02`, `04`, `05`, `07`, `09`, all of them |
| **Decided** | — |

### D4 · The level of automation / where the human approval gates go
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **2 gates** — ① fixing the spec (after S2) ② immediately before production starts (during S7). Between them, rule-based and automatic |
| **Options** | (a) 2 gates (b) one at every stage (c) fully automatic |
| **Grounds** | Gate ① prevents the failure of "simulating the wrong system" and gate ② the failure of "casually starting a multi-day job". The rest of the automation has a large benefit and a small risk. |
| **Impact** | `01` (the state machine), `06` (the budget) |
| **Decided** | — |

### D5 · Support for voice (recorded) input
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **Optional in v1**. If it becomes necessary, local transcription with `faster-whisper` (no API needed, works on the M4) |
| **Grounds** | S1 stands up with text, images and PDFs alone. Voice only adds one more transcription step, so adding it later does not change the structure. |
| **Impact** | `02_system_spec.md` (S1 input channels) |
| **Decided** | — |

### D6 · Version control (git)
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | **`git init` — a private repository.** `outputs/` and `knowledge/raw/` in `.gitignore`. The commit hash is recorded in `run_state.yaml` provenance. Publication later, after tidying (`D27`) |
| **Options** | (a) no git (b) private git (c) public git from the start |
| **Grounds** | `~/Desktop/BD_agent` is not currently a git repository. Neither was the preceding project. **To claim reproducibility you have to be able to answer "exactly what code produced this result".** All the more so in an agent whose auto-repair loop changes parameters.<br>**Why it was settled 2026-07-27:** user-fixed — *"develop privately → tidy up and publish later"*. (c) was not chosen because `knowledge/source/lab/` (a senior's simulations) and intermediate research results may be unpublished. But **since eventual publication is the goal, commit hygiene is kept to public standards from the start** — cheaper than erasing history later. |
| **Impact** | `master_plan.md` §0·§4 (the publication boundary table), `01` (the provenance block), `10_roadmap.md` |
| **Decided** | **2026-07-27** |

---

## B. Physics · numerics

### D7 · The reference time unit
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **`τ_D = σ²/D₀ = σ²γ/kT`** (the Brownian time to diffuse across one diameter) |
| **Options** | (a) `τ_D` (b) the MD convention `σ√(m/ε)` |
| **Grounds** | In overdamped BD the mass effectively never appears. Using the MD convention would make a physically meaningless number the reference. |
| **Correction 2026-07-27** | The original grounds also said *"`Δt/τ_D` then serves directly as the integration stability indicator"*, and **that part was wrong.** The stability indicator is the **displacement per step**, `√(2DΔt)/σ`. Across 3 measurements `Δt/τ_D` spanned a factor of 50 while the displacement spanned 7 — because the displacement is the square root of `Δt/τ_D`. **`τ_D` is kept as the reference time unit** but taken out of the gate. Grounds: `knowledge/wiki/findings/dt-gate-should-be-displacement-based.md` |
| **Impact** | all of `03_units_nondim.md`, `05` (the dt gate) |
| **Decided** | — |

### D8 · The default excluded-volume potential
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) — **WCA.** Decided by measurement |
| **Decision** | **`WCA`** (`hoomd.md.pair.LJ` + `r_cut = 2^{1/6}σ` + `mode="shift"`). `system.interaction: wca` · tuning in the `core` block. **The `ε` default is left at 1 k_BT, but 10 is used in dense systems (`φ_eff ≳ 0.3`)** — the table below is the grounds |
| **How it was decided** | `Z(φ)` was measured the way card §8 prescribed in advance and compared against Carnahan–Starling. The distillation is [`wca-reproduces-carnahan-starling`](../../knowledge/wiki/findings/wca-reproduces-carnahan-starling.md), and the 14-run sweep is in `notes/simbot/eos_wca/`.<br>**Keeping to the order changed the conclusion** — first the estimator was checked in the dilute regime against `1+B₂ρ+B₃ρ²` (the exact value). With `B₂` alone, +0.42% (4σ) remained at `φ`=0.02, which could be read as *"the estimator is wrong"*; adding the third-order term made the residual vanish to **+0.01%**. Only **after** confirming the measurement path was right did it go to the literature (the discipline of `SD12`).<br>The `φ_eff` mapping uses `B₂` (the `B₂*` axis of card §4) — feeding a `φ` measured with `σ` straight into CS is 16% wrong at `φ`=0.45 |
| **Measured — deviation from CS** | `ε`=1: −0.32% (`φ_eff`0.11) · −0.51% (0.21) · **−1.84%** (0.32) · −5.09% (0.42) · −8.43% (0.47)<br>`ε`=10: **+0.53%** (0.13) · **+0.56%** (0.38)<br>`ε`=0.5: −0.16% (0.10) · −1.99% (0.29) · −8.76% (0.44)<br>**`ε`=10 is 3× more accurate at `φ`=0.3 and costs 1.5×** (wall clock 46.2 s vs 30.4 s, `dt` 1.25e-4 vs 1.5e-4). That `ε`=1's deviation is **negative and monotonically increasing** is the signature of a soft core — as the density rises the particles dig in further and cannot be captured by a single effective diameter |
| **What became of the earlier warning** | The 2026-07-27 hypothesis (*"WCA's danger is not the potential but its coupling with the displacement per step"*) **was right.** Choosing `dt` on a force-displacement criterion (`plan.wca_force_limited_dt`) gives 0 box escapes up to `ε`=10 · `φ`=0.45. It is the same reason the lab code's `ε/kT`=500 runs stably. `D8` turned out not to be a question of *"which potential"* but of *"was the displacement gate observed"* |
| **The limits — honestly** | **Of the three candidates only WCA was measured.** The harmonic core and Wang–Frenkel were not, so this is not *"the best fit"* but **"the first one that fits well enough"**. WCA went in as the baseline first because it is the only candidate with a standard mapping that connects to the hard-sphere literature (CS, `φ_freeze`), and coming in within 2% weakened the reason to measure more. If the 5~8% in dense systems becomes a problem, the other two get measured then |
| **A side benefit** | In the run whose `φ_eff` crossed the melting point (`ε`=10 · `φ`=0.45 → `φ_eff`=0.565) the CS deviation was −29%, and **that was the result of extrapolating CS outside the fluid.** In the same run `D_fit` dropped by two orders to 0.0073 (its neighbour: 0.3507), so **an independent observable confirmed the phase transition** (`A1`). `analyze_eos` now attaches `cs_valid: false` with a reason — reporting a deviation from an extrapolation as a mismatch misreads physics as a defect of the method |
| **Default** | ~~undecided — to be fixed in stages C/D.~~ The candidates had been: bounded harmonic (Table) / WCA / Wang–Frenkel |
| **Grounds** | The warning the preceding project left: *"WCA's r⁻¹³ core is dangerous in the overdamped case. The displacement in one step is `F·dt/γ`, so even a small overlap flings the particle out of the box"* (`src/forces.py:117`). On the other hand, reproducing hard-sphere phase behaviour means it must not be too soft. **Decided by measuring against a literature benchmark (Carnahan–Starling).** |
| **Hypothesis 2026-07-27** | The lab's public code `graybox_abp_mpc` **runs stably with a very strong WCA, `ε/kT = 500`.** Its displacement per step is 0.018σ against the preceding project's 0.045σ, 2.5× smaller.<br>→ The hypothesis that **WCA's danger is not the potential itself but its coupling with the displacement per step**. Keep the displacement below 0.02σ and even a strong WCA may be stable. If true, `D8` turns from "which potential" into "was the displacement gate observed".<br>**A reproduction experiment is needed** — a 2D sweep of `ε/kT` × displacement. Grounds: `knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md` |
| **Impact** | `simbot/spec.py` (`INTERACTIONS`, `Core`) · `simbot/build.py` (`make_core`) · `simbot/plan.py` (`wca_force_limited_dt`) · `simbot/eos.py` (new) · `simbot/analyze.py` (`analyze_eos`) · `tests/test_eos.py` 22 cases · `benchmarks.yaml`'s `carnahan_starling_hs_eos` **can be unblocked** |
| **Decided** | **2026-07-28** |

### D9 · Dimensionality (2D / 3D)
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **3D by default**, with 2D supported by a spec flag (`hoomd.Box(Lx,Ly,Lz=0)`) |
| **Grounds** | The preceding project was 2D, and interesting physics like the hexatic ψ₆ and 2D melting is in 2D. As long as the code does not hardcode the dimensionality, the cost is almost nil. |
| **Impact** | `02`, `04`, `07` |
| **Decided** | — |

### D10 · Polydispersity
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **Unsupported in v1 (monodisperse)**. But the fields are opened in the spec schema in advance |
| **Grounds** | Real colloids are always polydisperse, and **polydispersity is essential in a dense system for suppressing crystallization to see glass and gel states**. But the v1 benchmark (hard-sphere phase behaviour) is a monodisperse reference value, so it has to be matched monodisperse first. |
| **Impact** | `02`, `09` |
| **Decided** | — |

### D11 · Hydrodynamic interactions (HI)
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **Ignored (free-draining)**. The approximation is stated explicitly in the report |
| **Grounds** | HOOMD's BD/Langevin does not include HI. It can be quantitatively wrong for sedimentation, shear and dense systems, so **writing "we did not do this" in the report rather than hiding it** matters. (inheriting the preceding project's *"What this is NOT"* convention) |
| **Impact** | `03`, `12` (the report caveat section) |
| **Decided** | — |

---

## C. Software

### D12 · The configuration file format
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **YAML → a dataclass per section** (inheriting the preceding project's `src/config.py` pattern). No pydantic |
| **Grounds** | It is an already-proven pattern and familiar to the user. What matters is not the format but two habits: ① `validate() -> list[str]` to **collect everything rather than stopping at the first error** ② unknown keys flow into an `extra` bucket so that future fields do not raise. |
| **Impact** | `02`, `bdkit/config` |
| **Decided** | — |

### D13 · Parameter sweeps / provenance management
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **v1 is directory-based** (`outputs/<run_id>/`). **signac is v2** |
| **Grounds** | `signac` is exactly the right tool for "parameter space + resume + history" (the same group as HOOMD), but it has a conceptual overhead. Better to bring it in when sweeps actually start to hurt in v1. |
| **Impact** | `01` (the artifact layout), `10_roadmap.md` |
| **Decided** | — |

### D14 · The particle renderer
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | **`fresnel` + matplotlib. OVITO is not used.** The roles split — **matplotlib for 2D systems, fresnel only for 3D.** A 2D monolayer has no depth to shade, so ray tracing gives nothing, and matplotlib is better at exact radii, scalar colour encoding and vector output (PDF) |
| **How it was settled** | The default (OVITO first) was **overturned.** The criterion was *"one headless render on macOS"* and `fresnel 0.13.8` passed it with `pathtrace` (320×320, 60 spheres, embree CPU). OVITO failed at two steps: ① conda-forge's `ovito` is **the GUI app only**, so it installs `bin/ovito` and no module in `site-packages` (the absence of `py312` in the build string was the tell) ② the real Python module (OVITO's own channel, the `py312` build) **downgrades `tbb` 2023.0.0 → 2022.3.0 and unlinks and relinks `fresnel`.** `hoomd 7.1.0` is not unlinked, so it is left sitting on a changed `tbb` — the kind of thing that shakes the reproducibility baseline, so it was stopped. `--freeze-installed` does not prevent transitive dependencies |
| **The price** | OVITO's analysis modifiers (CNA, clustering) are lost. **There is no substantive loss** — `freud 3.5.0` has `Hexatic`, `Voronoi`, `StaticStructureFactorDirect` and `DiffractionPattern`, all of them |
| **Fallback** | None. If fresnel fails it drops to a matplotlib scatter plot (which is the default even now) |
| **Conditions for revisiting** | When fresnel's sphere rendering becomes insufficient for dense 3D systems. Even then OVITO goes in a **separate env** and is handed only the GSD — `hoomd_slit` is not touched |
| **Impact** | `08_visualization.md`, `master_plan.md` §5(S7.5)·§13·§14-B, `simbot/viz.py`, `simbot/environment.yaml` |
| **Decided** | **2026-07-28** |

### D15 · The trajectory storage strategy
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **Structural observables need GSD → store GSD. But minimize the `dynamic` fields** and split the scalars off into `HDF5Log` |
| **Grounds** | The preceding project got 138MB→<1MB by accumulating in memory, but that was possible because the observable was only `(x,y)`. `g(r)`, `S(k)`, `ψ₆` and Voronoi need the full coordinates per frame. Instead, the storage frequency and the fields are reduced. |
| **Impact** | `05` (the disk budget), `07` |
| **Decided** | — |

### D16 · How error bars are produced
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **block averaging** (Flyvbjerg–Petersen 1989) by default, with the autocorrelation time `τ_ac` alongside |
| **Options** | (a) block averaging (b) bootstrap (c) the standard error of an independent-seed ensemble |
| **Grounds** | A standard error that ignores the time-series correlation **underestimates the error by factors**. This is the most common error in reporting simulation results. The rule: **do not produce a number without an error bar.** |
| **Impact** | `07_observables.md`, `05` (the N_eff gate) |
| **Decided** | — |

### D17 · The conda environment
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | A new env **`bd_agent`**, Python 3.12 (the installed HOOMD 7.1.0 is a `cpu_py312` build) |
| **Grounds** | Do not contaminate the existing `hoomd_slit` (py3.12.13). Commit `environment.yml` to make it reproducible. |
| **Impact** | `README.md`, task #8 |
| **Decided** | — |

---

## D. Operations · safety

### D18 · The default budgets (the ceiling on the auto-repair loop)
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | `max_total_walltime = 6h` · `max_repair_iterations = 8` · `max_disk_gb = 20` · `max_llm_calls_per_run = 100` |
| **Grounds** | Without a ceiling on the auto-repair loop the agent quietly burns days. **The numbers may be arbitrary but they must exist.** On exhaustion, escalate to a human. |
| **Impact** | `01` (run_state.budget), `06_repair_policy.md` |
| **Decided** | — |

### D19 · The range within which the LLM may make decisions
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | The LLM **proposes, diagnoses and narrates only**. Numerical computation, the PASS/FAIL verdict and the final fixing of parameters are deterministic code |
| **Grounds** | `decision_journal.jsonl`'s `actor` field records `rule` / `llm` / `human` separately. Only then can "how many times did an LLM judgment enter this result" be counted after the fact. **Overturn this item and the whole architecture changes.** |
| **Impact** | all of `01` |
| **Decided** | — |

---

## E. The knowledge layer · the rule layer (added 2026-07-27)

### D20 · How a senior's lab simulations go into the knowledge layer
| | |
|---|---|
| **Status** | **`DECIDED`** · **active: no** (2026-07-27 — no senior's code obtained) |
| **Decision** | Add **`kind: lab`** to the source kinds. `knowledge/source/lab/<year>-<author>-<slug>.md`. The required fields include **`reproduced: yes\|no\|partial`** |
| **Activation condition** | When a senior's code is **actually obtained**. Until then `master_plan.md` §7 keeps this section shrunk to 3~4 lines — spending page space on an asset you do not have blurs what to do first. On obtaining it, add `source/lab/` to `.gitignore` pre-emptively (the `D27` boundary table) |
| **Grounds** | User-fixed (2026-07-27): *"the wiki will also include simulations written by seniors in the lab."*<br>**A paper is "a published result" and a senior's code is "a parameter set that actually ran."** Papers often do not state `dt` or the number of equilibration steps, so as a parameter prior the latter is far more valuable. In exchange, **unverified practice comes along with it.** |
| **Discipline** | **Do not cite a parameter marked `reproduced: no` as if it were literature grounds.** Until reproduced it is a factual record that "this is what was done", not grounds that "this is right". Let that distinction collapse and the wiki becomes a rumour store rather than verification infrastructure. |
| **Impact** | `knowledge/wiki/CLAUDE.md`, `master_plan.md` §7, `09_literature.md` |
| **Decided** | **2026-07-27** |

### D21 · How many behaviour rules to write to begin with
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | **4 in v1** — `axioms` · `deterministic-core` · `overdamped-stability` · `verify-against-literature`. The rest are written **after actually experiencing that failure**, citing the incident |
| **Grounds** | **Rules are not written in advance. They are written after being burned, citing the incident.** If "why does this rule exist" does not have a real incident with a date, a path and a cost attached, the rule soon becomes a ritual observed by people who do not know the reason. And then the grounds for judging whether a rule should be retired when circumstances change are gone. |
| **Exception** | `overdamped-stability` already has a real incident — `~/Research/MD_particle/brownian_slit_sim/src/forces.py:117`: WCA's `r⁻¹³` core flung a particle out of the box in the overdamped case |
| **Rule candidates** | `wall-hit-escalation` · `cycle-discipline` · `failure-is-a-finding` · `wiki-first-lookup` · `cost-gate` · `ask-the-question-first` · `error-bars-or-silence` · `compute-router` |
| **hooks** | The same principle. **After confirming a rule actually was not observed**, attach it in `warn` mode, and if that still does not work raise it to `block` (v2) |
| **Impact** | `.claude/rules/`, `.claude/CLAUDE.md`, `master_plan.md` §9 |
| **Decided** | **2026-07-27** |

### D22 · Whether to adopt a multi-perspective reconnaissance panel
| | |
|---|---|
| **Status** | **`DECIDED`** — **not adopted** in v1, held as a v3 candidate |
| **What was considered** | A reconnaissance pattern that dispatches 3–5 personas of differing perspectives in parallel per task, having them **write questions rather than answers** |
| **Grounds** | This is a device suited to **divergent reconnaissance of an unsolved problem**. When you do not know where to look, extracting questions from several angles is valuable.<br>Our v1 task — "simulate a given system correctly" — is **convergent**, so the correct path is largely fixed. The panel's benefit is unclear and it would likely only add cost and complexity. |
| **Conditions for revisiting** | When the scope widens to exploring the phase behaviour of a new system, or to tasks where which observable to look at is itself unclear |
| **Impact** | `master_plan.md` §11 (v3) |
| **Decided** | **2026-07-27** |

---

---

## F. Reflecting the user's situation (2026-07-27, second pass)

`master_plan.md` was written borrowing einstein's structure, but on confirming that **the premises that justified
that structure (an external grader, 23 problems, 380+ accumulated cycles) do not exist for BD_agent**, this is the
result of readjusting to the user's situation. All 5 below are planted in `master_plan.md`'s body as `> ❓ undecided · Dnn` markers.

### D23 · How the `agent/` layer calls the LLM
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-27) |
| **Decision** | **(a) Call the Anthropic SDK directly.** Enforce the output schema with tool-use. Only S2 ELICIT, which needs conversational round trips, is left room to be moved out to (b) later. Calls happen from the single place `agent/llm.py` (enforced by `tests/test_invariants.py`) |
| **How it was settled** | **It was overturned to (b) once and came back.** The environment had neither `anthropic` nor an API key, only the `claude` CLI, so *"(b) works today"* was the judgment — but on actually calling it: `Your organization has disabled Claude subscription access for Claude Code · Use an Anthropic API key instead`. **A key is needed either way, so (b)'s only advantage disappeared.**<br>The lesson: **do not decide from an environment survey alone. You have to call it once.** "It is installed" and "it works" were different propositions, and the cost of distinguishing them was one call. Deferring `D23` to "when S1 is actually implemented" was right.<br>Grounds: [`knowledge/wiki/findings/d23-sdk-backend.md`](../../knowledge/wiki/findings/d23-sdk-backend.md) |
| **Options** | (a) `client.messages.create()` + tool-use to enforce a JSON Schema (b) a Claude Code skill + a `claude -p` subprocess (c) the Claude Agent SDK |
| **Grounds** | Of the 6 stages an LLM touches (S1·S2·S4·S5·S9·S12), **5 are "a single call with a fixed schema"**, so (a)'s ease of pytest is a direct benefit. Only S2, which needs conversational round trips, is moved out to (b). The full comparison table is in `master_plan.md` §4.5.<br>Adopting (b) wholesale is einstein's way, but there the essence was divergent exploration in which the LLM had to grep the wiki for itself. For us the context each stage needs is deterministically fixed. |
| **Cost of overturning** | **Medium.** The prompt, the output schema and the validator are reused independently of the method, and only the call wrapper (`agent/llm.py`) is replaced. That is why it is fine to defer and start |
| **Decided** | **2026-07-27** — settled while implementing the S1 INTAKE sketch path |
| **Impact** | `master_plan.md` §4.5, `01` (the §1 layer diagram), all of `agent/` |
| **Decided** | — |

### D24 · Whether to admit experimental data as a fifth evidence layer
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **Not admitted in v1.** But `evidence_layer: 5` is reserved in the `benchmarks.yaml` schema |
| **Grounds** | The research topic includes microrheology and tracking, so measured trajectories (trackpy) may exist. This belongs to none of §6's four layers and is **an independent measured oracle**, the most valuable evidence there is in a domain with no grader.<br>**Why it is deferred anyway:** an experiment-simulation mismatch has too many candidate causes — ignoring HI (`D11`) · polydispersity (`D10`) · tracking error · the system simply being different. To become an evidence layer there first has to be **a rule for interpreting a mismatch**. Without one, nothing beyond "it does not match" is possible |
| **When it gets decided** | When a microrheology system is first handled (after M2) |
| **Impact** | `master_plan.md` §6, `07_observables.md`, `09_literature.md` |
| **Decided** | — |

### D25 · The unit of a "cycle" + the outer autonomous loop
| | |
|---|---|
| **Status** | `OPEN` (the unit) · **`DECIDED`** (the loop — not adopted in v1) |
| **Decision (the loop)** | **No outer autonomous loop is built in v1.** The pipeline is always started by a human with `bd-agent new`. There is no queue, no scheduler and no layer that "decides for itself what to do next" |
| **Default (the unit)** | **One human-started run = 1 cycle.** `run_id` ↔ one row of `cycle-log.md`, 1:1 |
| **Grounds (the loop)** | User-fixed (2026-07-27). In a domain with a weak verification oracle an autonomous loop carries a large risk of **running all night in the wrong direction** — with no grader there is also no way for it to learn it was wrong. einstein could bear an autonomous loop because it had an external grader, the arena.<br>**Derived:** step ⑤ (the human) of the escalation ladder is not a waiting state but **the end of the run (`ESCALATED`)** |
| **Why the unit still catches** | If one system is run three times changing only a parameter, is that 1 cycle or 3? §10's "first-try pass rate" changes completely on this |
| **Conditions for revisiting (the loop)** | When the 2 gates actually become tedious, and `benchmarks.yaml` becomes dense enough that the automatic verdict is trustworthy |
| **When it gets decided (the unit)** | When the first row of `cycle-log.md` is actually written (at the end of M1) |
| **Impact** | `master_plan.md` §5·§9·§10·§11·§12, `01` (the §2 entry point) |
| **Decided** | **2026-07-27** (the loop part) |

### D26 · How many papers to distil to begin with, and in what order
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | **Before M1, only 2~3 on free BD / Stokes–Einstein.** The rest as and when the relevant system comes up |
| **Grounds** | The knowledge assets currently in hand are **a bundle of research-topic paper PDFs and nothing else** (no senior's code obtained → `D20`). There are four topics — depletion gels · charged/DLVO · microrheology · dense systems/the glass transition — and these become the 4 indices of `knowledge/wiki/systems/`.<br>**Why not distil everything in advance:** filling the wiki in advance is the same failure mode as writing rules in advance (`D21`). **A distillation that is never used is never verified** — it sits in the wiki wrong and contaminates the next run.<br>Stage 1's goal is not "a sufficient wiki" but **getting `wiki-first` to hit for real, even once** |
| **When it gets decided** | At the start of M1 |
| **Impact** | `master_plan.md` §7, `09_literature.md` |
| **Decided** | — |

### D27 · The separation procedure at publication time
| | |
|---|---|
| **Status** | `OPEN` |
| **Default** | Maintain `master_plan.md` §4's **publication boundary table** throughout development → **one audit immediately before publication**: check the copyright of `source/papers/` + search the whole history for sensitive paths + re-initialize into a new repository if necessary |
| **Grounds** | `D6` chose private, but **eventual publication is the goal**. Judging what goes out at publication time means digging through history. Nailing it down now at folder granularity makes that work disappear.<br>The most dangerous path is `knowledge/source/lab/` — added to `.gitignore` pre-emptively the moment it is obtained |
| **Cost of overturning** | **Low now and high later.** That is the reason for keeping the boundary table now |
| **When it gets decided** | When publication is actually resolved on (expected after M4) |
| **Impact** | `master_plan.md` §4, `.gitignore` |
| **Decided** | — |

### D28 · Adding pipeline stages (S2.5 PREREGISTER · S7.5 EYEBALL) + making the checks cross-cutting
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | 12 → **14 stages.** ① **S2.5 PREREGISTER** — promote pre-registration from an appendage of S1 to an independent stage, written **twice: qualitative v0 (right after gate 1) → quantitative v1 (right after S4)** ② **S7.5 EYEBALL** — place the eyeball check of 3 low-resolution snapshots **before S8** ③ define the numerical-instability and physical-nonsense checks not as a separate stage but as **the exit check of every stage** (`master_plan.md` §14-E) |
| **Grounds** | While mapping the user's own description of their workflow in 9 steps (2026-07-27) onto the existing 12 stages, **three places disagreed, and the disagreeing side was right.**<br>**① Pre-registration:** the user picked out "drawing the expected result + presenting the physical grounds" as an independent step. In fact this is the **only** device by which the pipeline prevents post-hoc interpretation, so it cannot be an appendage. Writing it twice resolves a dilemma — written before non-dimensionalization the grounds are weak, written after it is no longer a pre-registration because the system has already been examined. Keeping both makes **the v0→v1 change itself data** (evidence that non-dimensionalization and the literature corrected the intuition).<br>**② The eyeball check:** the user put visualization (6) **before** analysis (7). The original plan was `S10 ANALYZE → S11 VISUALIZE`. Analysing a crystallized system as a glass / a cluster spanning the box / overlapping particles / an ongoing phase separation — all of these are failures that take **1 second by eye and are hard by number**. S11 (the final render for the report) stays as it is and one cheap check is added before it.<br>**③ Cross-cutting:** the user said "review **at each step**". The original S8 DIAGNOSE runs only once, *after* execution. Far more physical nonsense is caught before execution, and far more cheaply (`φ>0.64` · a `τ_B/τ_D` violation · `r_cut > L/2` · overlap in the initial arrangement). |
| **Cost** | 2 stages added = 9 implementation items (`14-B`). S7.5 finishes in a few seconds of rendering, so the runtime burden is negligible |
| **Impact** | `master_plan.md` §5·§14, `01` (§2 the state machine·§3 artifacts·§5 the layout) |
| **Decided** | **2026-07-27** |

### D29 · Splitting verification into two axes (process V1~V7 / results, 4 layers) + a 3-way verdict
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | ① Reorganize §6 into two axes, **process verification (§6-A, V1~V7)** and **result verification (§6-B, the four layers of evidence)**<br>② Create three layers the original plan did not have — **V1 fidelity** (back-translation), **V4 inter-stage consistency** (an S6 cross-check), **V6 Outlier** (4 kinds)<br>③ Move the verdict from 2-way (`PASS`/`FAIL`) to **3-way** (`PASS` / **`PASS-with-doubt`** / `FAIL`)<br>④ **A doubt lowers the ceiling of the evidence grade** (the §6-C table)<br>⑤ Make gate 1's identity concrete: "spec approval" → **"back-translation comparison approval"**<br>⑥ S8 diagnosis, 5 categories → **6 categories** (outlier added) |
| **Grounds** | Mapping the user's description of verification in 7 steps (2026-07-27) onto the document showed that **§6 covered result verification only.** Process verification was merely scattered across §14-E as flat exit checks, with no system to it.<br>**The two fail in different ways.** The process fails *quietly* (simulating the wrong system perfectly) and the result fails *plausibly*. **The former is not caught however well the results are verified** — even if all four layers of evidence agree, if it was a different system to begin with then it is a wrong answer on which everything agrees.<br>**V1 was an especially large hole.** The original S1's check was only "0 arbitrarily generated values", and that prevents *hallucination*, not *fidelity*. `"500 nm silica"` → `material: polystyrene, confidence: 0.9, unknowns: []` is not an invented value, so **it passes every check.** Back-translation comparison is the standard technique for catching this class, and it is something far better done than a human skimming YAML.<br>**V4:** every other check is internal to a stage. There was nowhere to catch the case where each stage passes individually but the combination is wrong (S5 changing `dt` and thereby invalidating the S3 gate, and so on).<br>**V6:** S8's 5 categories are all aggregates, so **an outlier disappears on aggregation.** In particular an ensemble outlier (one seed jumping) bears directly on reproducibility, and looking only at the mean mistakes an averaged accident for physics. |
| **Why a 3-way verdict** | With only `PASS`/`FAIL`, **anything doubtful has nowhere to go and disappears** — it passed, so nobody looks again, and when the result later seems odd there is no way to retrace "what was it that caught". Tying `PASS-with-doubt` to a grade ceiling is **to create an incentive to record it.** But it only works if a low grade is treated not as a *punishment* but as *an accurate description of the state* |
| **Cost** | 44 check items (`14-E`) — +25 from the existing 19. Most are threshold comparisons so the implementation is light, and **only V1 back-translation and V6 ensemble are substantive work** |
| **Impact** | `master_plan.md` §5·§6·§6-A·§6-B·§6-C·§9·§14, `01`(the gate 1 definition) |
| **Decided** | **2026-07-27** |

### D30 · The S1 reading model — is strong reasoning needed
| | |
|---|---|
| **Status** | **`DECIDED`** (provisional — before measurement) |
| **Decision** | **`claude-opus-5`**. Overridable with the environment variable `BD_AGENT_LLM_MODEL` |
| **Grounds** | At first `claude-sonnet-5` was set on the grounds of "a limited extraction task", but **that framing itself was wrong.** Looking at a real sketch (`tests/test_image*.jpeg`), the work needed is this.<br>① reading handwritten formulas — the subscript `k_t`, the fraction `½`, the Greek `µ`, roots<br>② **judging a drawing-label conflict** — 13 circles vs `100 particles`. Knowing which side is authoritative<br>③ recognizing that `A/r³ ≫ k_BT` is not a value but **a regime declaration**<br>④ judging that `find final configuration by minimizing U_tot` is not a description of the system but **the objective**, and that this task is energy minimization rather than BD dynamics<br>⑤ **the restraint not to fill in a blank** — the whole V1 layer exists to prevent this failure<br>All of it is strong reasoning. And it is **a single call at the very front of the pipeline**, so being wrong here contaminates everything below, and the cost is a few cents a call. Economizing here is a false economy. |
| **Unresolved** | **This is still an argument, not a measurement.** It falls short of the "do not guess, measure" principle that `D8` and `D14` established. It gets settled by attaching a per-model score on synthetic sketches with an answer key (`tests/fixtures/sketch_walls.truth.json`). After an API key is obtained. |
| **Impact** | `agent/llm.py`, `cli.py`, `master_plan.md` §4.5 |
| **Decided** | **2026-07-27** |

### D31 · The scope of the reading schema — does it hold formulas
| | |
|---|---|
| **Status** | **`DECIDED`** |
| **Decision** | **It does.** Schema v2 — `medium` · `interactions[]` · `external_potential[]` · `objective` · `relations[]` created. A potential is held **split** into `raw` (the original text) / `form` (normalized) / `params` (a value and unit per symbol). A drawing-label conflict is `count` (authoritative) + `drawn_count` (diagnostic) + `F11` (recorded only) |
| **Grounds** | v1 assumed the sketch was "a concept diagram of the system's topology", but **the actual drawing was a hand-written problem statement** — half of it is formulas. The `F8` omission check detected this design defect by itself (the `300` of `T = 300 K` had nowhere to go).<br>**`objective` was the biggest hole.** Unable to hold the objective, S5 PLAN plans time integration, an equilibration verdict and an MSD for an energy-minimization task. Everything passes and everything is meaningless — the "the process fails quietly" failure §6-A speaks of.<br>Why `raw`/`form` are split: handwritten formulas are frequently misread (`k_t`↔`k_B`, `r³`↔`r⁵`), so without the original text there is no retracing. `F12` enforces it.<br>Grounds: [`knowledge/wiki/findings/sketch-schema-v2-equations.md`](../../knowledge/wiki/findings/sketch-schema-v2-equations.md) |
| **What was kept** | `F2` stands — fill these values with `sketch:visual` and it still FAILs. **"If it is written, read it; if it is not, leave it empty."** Creating a field does not mean guessing is permitted |
| **Impact** | `bdkit/reading/sketch.py`, `fidelity.py` (F11, F12, F13 created), `backtranslate.py`, `agent/s1_intake/` |
| **Decided** | **2026-07-27** |

### D32 · Absorb `Simulation_bot` + resolve the `spec` name collision
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | ① `Simulation_bot/simbot/` → **`BD_agent/simbot/`**. One entry point, `cli.py` (`intake`, `check`, `elicit`, `compose`, `run`, `batch`)<br>② **`bdkit/spec/` → `bdkit/reading/`** — the contents are all sketch reading (`SketchReading`, `fidelity`, `gaps`, `elicit`), so the name did not match reality, and it collided with `simbot/spec.py` (the simulation configuration). There is now only one `spec` in the repository<br>③ **Unify the blank verdict** — delete `_ELICIT_TABLE` from `simbot/intake.blanks()` and delegate to `bdkit/reading/gaps.detect_gaps()`. `blanks()` only transforms the shape<br>④ **Lifted** the low-confidence confirmation check from `simbot/intake.py` up into `gaps.py` — a feature that was nearly lost in the delegation |
| **Grounds** | User-fixed (2026-07-28): *"the bot seems to work well; it would be good to integrate it into the existing system."* It is the answer to `SQ5` and the decision that opens `Simulation_bot/PLAN.md` §6's *"do not modify BD_agent"*.<br>**The two halves were already speaking the same language** — `readings/*.reading.json` is exactly the `SketchReading` v2 format, so the integration was shallow. The only real duplication was **the blank verdict**, and that is ③. |
| **How the verification-layer collision was resolved** | `PLAN.md` §2 deliberately deleted V1~V7 and the evidence grades, and `master_plan.md` §6 has them as its core. **`simbot` does not make claims, so it is not a subject of the V-ladder** — the V layers attach to *claims*, and the engine only produces numbers and pictures. At the moment verification is needed, the ladder wraps around it. Both live as they are |
| **A new invariant** | `SD2` as a test — `simbot/` does not import `agent` or `anthropic`, and apart from `intake.py` does not import `bdkit` either (the engine does not know the reading layer. That is what keeps alive the path of running from YAML with no drawing) |
| **Impact** | `simbot/` `bdkit/reading/` `cli.py` `tests/test_invariants.py` `docs/11_simbot.md`, 17 files of imports |
| **Decided** | **2026-07-28** |

### D33 · How the simulation interpreter is determined
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | **Search for it. Do not put the path in the code.** In order: `--python` → `$BD_AGENT_SIM_PYTHON` → `simbot/environment.local.yaml` (gitignored) → a search by conda env name → the current interpreter. If all of it fails, throw **with four ways to fix it** written out.<br>If necessary `cli.py` **re-executes itself** (`os.execv`) — because `--python` only sets the P5 child, while P6 and P7 run in the parent, so `gsd` and `numpy` are needed in the parent too |
| **Grounds** | User request (2026-07-28): *"can the default system be set to hoomd_slit? Separately for when it is distributed later."*<br>Writing `/opt/homebrew/.../envs/hoomd_slit/bin/python` into the code makes it **a repository that only runs on this machine**. The same logic by which `D2` prevented a hardcoded device with `make_device(spec)`.<br>**If, in a chatbot, you have to remember "which command with which python", that is a defect of the tool.** So the automatic switch was attached too |
| **The distribution boundary** | `simbot/environment.yaml` (**committed**) writes down only *what is needed* — `env_name` · `required` · `optional`. `simbot/environment.local.yaml` (**gitignored**) is *where that is on this machine*. Only the former is distributed, and the receiving side mostly just runs, via the search or a single env |
| **Bootstrap discipline** | `simbot/env.py` **uses the standard library only.** At first it imported `pyyaml` and died with `ModuleNotFoundError: yaml` — when running on an interpreter that lacks pyyaml is precisely this module's reason for existing. **A bootstrap cannot depend on the thing it is trying to find for you.** `tests/test_sim_env.py` enforces it by AST |
| **Impact** | `simbot/env.py` (new) · `simbot/environment.yaml` (new) · `cli.py` (the `env` command + the automatic switch) · `simbot/commands.py` · `.gitignore` |
| **Decided** | **2026-07-28** |

### D34 · The machine performance profile — measured per machine, not committed
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | `simbot/machine_profile.yaml` (**committed**) holds only *how it is measured* and a conservative fallback. The measurement is `machine_profile.local.json` (**gitignored**). On a new computer, run `python cli.py calibrate` once **before the first run** (quick ~2 min / `--full` ~8 min).<br>Write a **machine fingerprint** (the CPU name, the core count) into the profile alongside, and **ignore it** if it does not match |
| **Grounds** | `machine_profile.yaml` was committed holding the M4's number (10.4 M p-steps/s). Get it on another machine and it estimates cost with that number — **the same problem** as `D33`'s decision not to commit the interpreter path.<br>**Somebody else's profile is worse than no profile.** With none, it goes to a conservative default and a warning appears; with somebody else's, a plausible number comes out and nobody looks. That is why a fingerprint mismatch takes precedence over the fallback |
| **Parallel measurement (2026-07-28)** | **One run is single-threaded** — CPU utilization 0.99 cores. `mpi_enabled=False` · `gpu_enabled=False` · `device.CPU()` has no thread argument. **This is not something code can do anything about.**<br>K independent runs at once: K=2 efficiency 98% · **K=4 efficiency 81% (the knee)** · K=8 43% · K=10 38%. Maximum speedup 3.8×.<br>**The M4 has 4 performance + 6 efficiency cores**, so trusting `os.cpu_count()` (10) **overestimates by 2.6×** |
| **Derived** | The worker count is set not from the core count but from **the knee of the efficiency curve** (`recommended_workers()`, the largest `k` with efficiency ≥ 0.6). The throughput uses the measurement point at **the nearest `N`** rather than the overall mean — because the throughput depends on `N`. If it was extrapolated, that fact is left in the provenance string |
| **Impact** | `simbot/machine.py` (new) · `simbot/calibrate.py` (new) · `simbot/plan.py` (`load_profile` delegation) · `cli.py` (`calibrate`) · `.gitignore` · `tests/test_machine_profile.py` |
| **Decided** | **2026-07-28** |

---

### D35 · Verdict thresholds are born from human review — a partial answer to `D4`
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | **Early on a human reviews, and once measurements accumulate the threshold is drawn from that distribution and automated.** The tool judges **only observables with a registered threshold**, and passes the rest to a human as before. The registry is `knowledge/wiki/benchmarks/benchmarks.yaml` |
| **Options** | (a) do not judge (as until now) (b) attach a few thresholds for partial automation (c) automate the whole verification ladder |
| **Grounds** | User decision: *"given the nature of experiments there will always be error, and the threshold range will probably keep changing with the system. How about a human reviewing early on and automating once experience accumulates"*.<br>**This sets the path (b)→(c).** The reason the tool has not judged until now is that there were no thresholds (`CLAUDE.md`: *"it is a statement you can only make if you have a threshold, and that threshold is not yet in this repository"*), and inventing a threshold stamps **a plausible and wrong pass** — the most expensive failure there is in a domain with a weak verification oracle (`master_plan` §6).<br>**One key requirement:** for a review record to become a threshold it has to be left **per observable, together with the number of times it was seen**. Attaching "OK" to one run cannot make a threshold — you need a record that `D_fit` was measured three times and came out 0.9911 three times before you can write "within 1%".<br>There is already material in that direction: free BD `D_fit`, 3 runs agreeing to 6 digits (§7 #4·#16·#21) · the trap pentamer within 0.6σ (Step 1.5) · the `A`=100 ensemble scatter (§7 #31). **The first candidates for automation are these three** — attach a threshold starting from what you know the distribution of |
| **Impact** | `knowledge/wiki/benchmarks/benchmarks.yaml` (the registry) · `docs/11_simbot.md` §2's *"does not judge"* becomes **"judges only what is registered"** · connects to the gate discussion of `D4` (the level of automation) |
| ~~**What stays unresolved**~~ | ~~the format and location of the review record~~ · ~~how many entries before registering~~ → **both settled** (below) |
| **① Location — one ledger in the wiki** | `knowledge/wiki/benchmarks/reviews.jsonl` (**committed**). The reason it is not the run folder is decisive: **`runs/` is gitignored.** A review verdict is *knowledge*, so put it there and one cleanup destroys it, and it was never committed. It goes in the same folder as the threshold registry (`benchmarks.yaml`).<br>**The reason it is JSONL is concurrent appends** — `cli.py batch` runs with several workers, each a different process, and a read-modify-write YAML loses updates. Appending one line is atomic via `O_APPEND` (the same reason as `machine_measurements.local.jsonl`).<br>**Append-only** — to overturn a verdict, insert a new line and let the aggregation take the latest as authoritative. Delete and the evidence of *"how the threshold was born"* is gone |
| **② How many — not runs but 3 independent configurations** | The prohibition in `A1` (*"counting runs that differ only in seed or box size as different evidence"*) is **enforced by code**: the sha256 of the spec with `run.seed` removed is the fingerprint, and the same fingerprint counts as one. So running 5 seeds gives `n_independent`=1 and no proposal comes out (`tests/test_review.py::test_five_seed_runs_do_not_make_three_independent_evidence`).<br>3 is a convention and there is one ground for it — estimating a scatter needs at least three. The proposed tolerance is **the largest observed deviation × 1.5**: use the maximum as it stands and half of the next runs fall out (the maximum of a sample is not the maximum of the distribution).<br>**It does not register automatically.** `--propose` even produces the YAML fragment to attach to `benchmarks.yaml`, but attaching it is the human's job |
| **③ A session cannot make a threshold out of its own verdict** | This is the third problem, which surfaced during implementation. A session (the LLM in the conversation) reading the numbers and **proposing** a verdict is helpful, but counting that stamps *a pass no human ever looked at*. So `reviewer` is recorded and **the aggregation separates them** — `propose()` counts only `human` by default, and what was left with `--as-session` appears only in the preview. If a human enters a verdict on the same `(run_id, observable)`, that overwrites as authoritative |
| **The first return — a defect was caught in the first run** | The table said to compare `free_bd`'s `D_fit` against `D_theory`=1 but **with excluded volume 1 is not the expected value** (it is the non-interacting limit). Without a gate the ledger **accumulates −52% as a mismatch** at `φ`=0.3, and drawing a threshold from that material gives a 78% tolerance. It was blocked with `invalid_if="core_kind"` and nailed down with 2 tests — **the failure this layer was built to prevent came out while building this layer** |
| **Impact** | `simbot/review.py` (new) · `simbot/commands.py` (`cmd_review`) · `cli.py` (`review`) · `tests/test_review.py` 24 cases · `knowledge/wiki/benchmarks/reviews.jsonl` (new, committed) |
| **Decided** | **2026-07-28** |

### D36 · The stage is inside Claude Code — and distribution assumes the same
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | **Used inside Claude Code only.** No separate UI or headless entry point is built. Distribution to the Park lab assumes the same, so what is needed is **an installation procedure + a few example configs**, not a new interface |
| **Options** | (a) inside Claude Code only (b) support users who only use a terminal CLI |
| **Grounds** | User decision: *"I am thinking of using it only inside Claude Code"*. It extends `SD11` (the stage is Claude Code) and `SD22` (the session does the reading) to distribution, and those two already stood on this premise — **that the LLM doing the reading is inside the session** is a pillar of the design.<br>The price of not choosing (b): the receiving side needs Claude Code to start from a drawing or natural language. But **since the path of starting from a finished config is open with the CLI alone** (`run`, `batch`, `session set`), the price is not large |
| **Impact** | The distribution requirements narrow to *"example configs + installation documentation"*. Running `cli.py calibrate` (`D34`) once before the first run enters the procedure |
| **Decided** | **2026-07-28** |

### D37 · One authoritative record of progress — `master_plan` §14's checklist was folded
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | **The authoritative record of engine progress is [`docs/11_simbot.md`](11_simbot.md) §5·§7 and nowhere else.** `master_plan` §14's `14-A` (11 infrastructure items) and `14-B` (72 per-stage items) were deleted and replaced with a pointer. What was kept: `14-C` (verification infrastructure), `14-D` (knowledge and rules), `14-E` (the 7-layer verification ladder, 44 items) |
| **Options** | (a) re-map the 144 items onto the current structure (b) **fold them and make §5 authoritative** (c) leave them alone |
| **Grounds** | User decision. The reason (a) was not chosen is cost — each of the 144 rows needs **three judgments**: ① is it done (the artifact cell points at a nonexistent path like `bdkit/{spec,units,…}/`, so you have to open the code to know) ② is it still wanted (a good many are layers removed in the chatbot pivot, so an `X` is a *deletion* rather than "to do") ③ where is it now. Every row needs its grounds looked at, so 2~5 hours, and it does not parallelize.<br>**That the return is small matters more.** With progress recorded in two places they will inevitably diverge, and in fact this session too updated only §5 without touching §14. The denominator `2 / 144` **was false** — things functionally finished were sitting there as `X`.<br>**Leaving a removed layer as an `X` is especially bad** — dozens of nonexistent backlog items appear, and then nobody trusts the list at all |
| **Why the rest was kept** | `14-E` is **not yet started and exists nowhere else.** And `D35` (human review → threshold registration → automatic verdict) is the entry path to this layer, so it just came alive. `14-C` and `14-D` are only 22 rows, so they were re-marked with their actual state (individual benchmark items are better tracked by `benchmarks.yaml` together with `blocked_by`) |
| **Impact** | `master_plan.md` 1417 → 1249 lines. *"Engine progress is not here"* was nailed into the §14 header |
| **Decided** | **2026-07-28** |

### D38 · The parameter-grounds ledger — leave the value and **the reason** on the same line
| | |
|---|---|
| **Status** | **`DECIDED`** (2026-07-28) |
| **Decision** | Per run, leave **every decision in the same shape** in `runs/<id>/parameters.yaml` — `path` · `value` · `decided_by` (user/reading/rule/default) · `why` · `rule` · `source` · `alternatives` · `advice`. Reading across is `python cli.py params`. The implementation is `simbot/rationale.py` |
| **Grounds** | User decision: *"choosing the modelling and the parameters for a simulation well is what matters for verification through simulation. So that this modelling and parameter selection can be done appropriately and well, it is important to leave a record of the system and of the reason and value of the parameter modelling. What a human used to do by accumulating experience — choosing appropriate parameters — learned from a large body of material to make and carry out efficient proposals."*<br>**It was not that the grounds were absent but that they were scattered and in different shapes** — `spec.provenance` (only values coming from the reading) · `run_plan.notes` (prose, so not aggregatable) · `dt_constraint` (a single string) · `decision_log` (repository-level, so not this run's values). They have to be in the same shape to read **across several runs**, and being able to read across is what answers *"what did we choose before in a similar system"*.<br>**Separating out `decided_by` is the crux of this ledger.** A default 0.35 and a human-entered 0.35 **have the same value and are indistinguishable by eye.** Unable to distinguish them, the corpus accumulates **the false practice** that *"in this system we use 0.35"*. So `spec.specified_paths` (the paths the user actually wrote) is newly carried around.<br>`alternatives` widens `SD9`'s *"also write down the value the rule you did not adopt would have given"* to every parameter |
| **The first return** | Reading 4 runs' ledgers across immediately exposed **a parameter set nobody had ever chosen** — `init.placement` (lattice) · `init.min_separation_over_d` (0.35) · **`core.epsilon_kT` (1.0)**. The third is especially heavy: it is the value that sets how hard the core is, and in every run it was a dataclass default — the same spot as the reason `D8` (the excluded-volume potential) is open |
| **The boundary** | **The engine does not make the proposal.** Accumulating the corpus and answering queries is the engine's job, and *"this value would be good for this system"* is said by the session (the LLM in the conversation) (`SD22` · `A4`). If the engine also proposed, then when a proposal is wrong **there is no telling whether the record was wrong or the judgment was** |
| **No retroactive fill** | No ledger is made for the previous 29 runs. Without `specified_paths` there is no restoring *"what did a human write and what was a default"*, and guessing to fill it quietly contaminates the corpus. **The corpus accumulates from now on** |
| **Impact** | `simbot/rationale.py` (new) · `simbot/spec.py` (`specified_paths`) · `simbot/commands.py` (save after P3 · `cmd_params`) · `cli.py` (`params`) · `tests/test_rationale.py` (14 cases) |
| **Decided** | **2026-07-28** |

---

## Not yet itemized (to be promoted to D32+ later)

- When the LLM fills a value absent from the spec with a literature default, should human confirmation be forced beyond `assumed: true`
- The retention period / automatic cleanup policy for the artifacts of a failed run
- The format of a meta-report comparing several runs
- How BD and the pipeline get shared when HPMC is introduced
- **The wiki promotion criterion** — what should be the threshold for raising a finding to a concept (3 citations? appearing in 2 different systems?). Only that promotion is a human approval has been settled
- **How an `unverified` result is displayed in the report** — a badge? A warning box? How conspicuous should it be
- **To what level reproduction of a senior's code should be required** — the criterion for a verdict of `reproduced: partial` (1 observable matching? all of them?)
- **The set of subcommands of the `bd-agent` CLI** — what is needed beyond `new` / `resume` (`status`? `ls`? `report`?)

> **Promoted:** "comparison against experimental data" → `D24` · "the unit of a cycle" → `D25` (both 2026-07-27)
