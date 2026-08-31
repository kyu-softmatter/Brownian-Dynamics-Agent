# BD_agent — the master plan

> Reading this one document should answer three things: **(a) what is being built (b) why it is built that way
> (c) what to do next**. If it does not, that is a defect of this document.
>
> Last updated: 2026-07-27 · status: **M0 (design)** · repository: **private** (to be published after tidying → `D27`)

---

## 0. The next three things

Opening a session, this is all you need to look at. When these three are done, M0 is done.

| | To do | Why now | Resolves |
|---|---|---|---|
| **1** | `git init` (private) + `.gitignore` + `environment.yml` | the premise of provenance. Once the auto-repair loop starts changing parameters, *"exactly what code produced this result"* cannot be answered | `D6` `D17` |
| **2** | measure `machine_profile.yaml` — free BD, `N` = 1k / 5k / 10k | the §8 cost gate does not work without this file. Without the gate the agent casually starts a multi-day job | `A3` |
| **3** | the `bdkit/` skeleton + the benchmark `free_bd_stokes_einstein` passing — **0 lines of LLM** | half of M1's physics |  `A4` |

> **Number 3 is the crux.** Build the agent layer first and, when a result is wrong, there is no telling whether the
> physics is wrong or the LLM is. **Prove the core is right on its own, then put the LLM on top.** This is the real bridge M0→M1.

The detailed items and their O/X progress are in **[§14 the implementation checklist](#14-the-implementation-checklist)**. §11 is *when*, §14 is *what*.

---

## 1. One sentence

**An agent that turns a colloidal system described in natural language and pictures into a verified Brownian
Dynamics simulation, and accumulates the knowledge gained along the way so the next simulation can use it.**

Voice input is optional in v1 (`D5`), and MC/HPMC is v2 (`D3`). Putting them in a one-sentence definition blurs the v1 scope, so they are left out.

---

## 2. Why — the three problems of doing BD/MC by hand

### 2-a. Parameter selection is tacit knowledge

What to set `dt` to, where to cut `r_cut`, how long to run the equilibration — these are mostly **values learned from
a senior or imitated from a paper**. It is rare to be able to explain why that value, and there is no procedure for
confirming it is still valid when moved to a new system.

### 2-b. Verification is unsystematic

A simulation **always produces some number.** A `g(r)` always gets drawn, and an MSD always looks like a straight line.
The problem is that **you cannot tell by eye whether it is the right `g(r)`**. Not having diverged
does not mean being correct.

### 2-c. Failures are not recorded

"I tried these parameters and it did not work" mostly stays in that person's head alone. Six months later the same
person, or the next student, hits the same wall again. **Failures disappear and only successes remain in papers.**

### What this agent is trying to do

| Problem | Response |
|---|---|
| tacit knowledge | attach a **dimensionless-number basis** and a **literature source** to every parameter |
| unsystematic verification | run **literature benchmarks as regression tests**. Attach an **evidence grade** to every result |
| disappearing failures | leave a failure as a `dead-end` page, and make **the next attempt look it up first** |

---

## 3. Design principles

> ### The LLM proposes, and deterministic code decides.

| What the LLM does | What the LLM never does |
|---|---|
| natural language, pictures, voice → extract a draft spec | numerical computation |
| literature comparison, presenting a parameter prior | the verification pass/fail verdict |
| triaging the cause of a verification failure | fixing a parameter finally |
| the report narrative | generating a physical conclusion |

The 10 general principles for building an agent are in [`docs/01_agent_architecture.md`](../../docs/history/2026-07_bd_agent_01_agent_architecture.ko.md) §0.
Of those, the three that weigh most in this project:

1. **Keep state in files** (`run_state.yaml`) — conversation context disappears and cannot be inspected
2. **Journal every decision, distinguishing down to the `actor`** (`rule` / `llm` / `human`) — it has to be possible
   to count afterwards "how many times did an LLM judgment enter this result"
3. **A verdict criterion has to be a number** — "it looks plausible" is not a gate

---

## 4. Three layers, three contracts

The repository is split into three layers, and each layer gets **its own contract document**. Because the rate of
change and the verification method differ per layer — code is verified by tests, knowledge by citation, rules by an incident record.

| Layer | Path | Contract document | Content |
|---|---|---|---|
| **Code** | `bdkit/` `agent/` `tests/` | `CLAUDE.md` (root) | the deterministic core + a thin LLM layer |
| **Knowledge** | `knowledge/` | `knowledge/wiki/CLAUDE.md` | distillations of papers and seniors' code + a synthesized wiki |
| **Rules** | `.claude/rules/` | `.claude/CLAUDE.md` | behaviour rules, one topic per file |

> It is the same structure `D1` calls "four layers". Counting the deterministic core `bdkit/` and the LLM layer
> `agent/` separately gives four; bundling them as the Code layer that shares one contract document gives three.
> **This document counts by contract document — three.**

### The publication boundary — fixed from now

The repository is **developed privately and published later after tidying** (`D27`). Judging what goes out at
publication time means digging through history, so **it is nailed down now at folder granularity.**

| Path | At publication | Reason |
|---|---|---|
| `bdkit/` `agent/` `tests/` `docs/` `.claude/` | **goes out** | the tool itself is what is being published |
| `knowledge/wiki/` | **goes out** | synthesized knowledge — only citations, no originals |
| `knowledge/source/papers/` | **needs judgment** | 1:1 distillations. Check the copyright boundary before publication. **Our lab's published papers are here too** — `lab_authored: true` |
| `knowledge/source/lab/` | **does not go out** | **unpublished** lab assets only — code, notes, unpublished parameters (once assets are obtained) |
| `knowledge/raw/` | **does not go out** | the original PDFs and code |
| `outputs/` | **does not go out** | run results. The volume is also a problem |

In v1 `.gitignore` blocks only `outputs/` and `knowledge/raw/`. `source/lab/` gets added pre-emptively **the moment a
senior's code is actually obtained** — once committed it has to be erased from history.

> ❓ **Undecided · `D27`** — how is the separation at publication time to be carried out?
> **The current default:** keep to the boundary table above, then **one audit immediately before publication** —
> check the copyright of `source/papers/` + search the whole history for sensitive paths + re-initialize into a new repository if necessary
> **When it gets decided:** when publication is actually resolved on (expected after M4)
> **Cost of overturning:** **low, now.** Later it is high — that is the reason for keeping the boundary table now

> **What divides them is not the author but whether it is published.** Put our lab's papers in the same folder as
> unpublished assets and, at publication, **the distillations of papers that already have a DOI and are already
> public drop out along with them.** What needs protecting is not the author but
> *the fact that it is not yet published*. Lab authorship is marked not by the folder but by `lab_authored: true`.
> Details in [`knowledge/wiki/CLAUDE.md`](../../knowledge/wiki/CLAUDE.md).

### The invariant of the code layer

```
bdkit/     the deterministic core — 0 LLM dependencies.  grep -r "anthropic\|claude" bdkit/ must come back empty
agent/     LLM layer — three pieces: prompt + output schema + validator
```

Without this separation, when a result is wrong **there is no ever knowing whether the physics is wrong or the LLM is.**
In simulation that is fatal — because a plausible-but-wrong `g(r)` is indistinguishable by eye.

### The rule layer is born from incidents

> **Rules are not written in advance. They are written after being burned, citing the incident.**

If "why does this rule exist" does not have **a real incident with a date and a cost** attached, the rule soon becomes
a ritual observed by people who do not know the reason. And then the grounds for judging whether a rule should be
retired when circumstances change are gone.

**v1 starts with 4** (§9). One of them already has a real incident. The rest get written after being experienced.

The rule file template (no frontmatter, 25–80 lines):

```markdown
# <slug> — <a one-line summary>

<the rule body. Imperative>

**Why (the triggering incident):** <a real incident with a date, a path and a cost>
**How to apply:** <a checkable list of actions>
**Anti-patterns explicitly forbidden:** <named failure modes>

See also: [other-rule](other-rule.md)
```

---

## 4.5. How the `agent/` layer calls the LLM — a comparison table

§4 defines `agent/` only as "three pieces: prompt + output schema + validator". **The actual calling method is not yet
decided.** The table below is the basis for the decision; it proceeds on the default and gets settled in M1.

| Axis | (a) the Anthropic SDK directly | (b) Claude Code + `claude -p` | (c) the Claude Agent SDK |
|---|---|---|---|
| where the code lives | `agent/*.py` → `client.messages.create()` | `.claude/skills/*.md` + an `agent/headless.py` wrapper | an SDK session inside `agent/` |
| structured output | **enforced** as a JSON Schema by tool-use | requested in the prompt + parsed and retried by hand | enforced by a tool definition |
| `pytest` | **easy** — 1 function, inject a mock | hard — a subprocess, non-deterministic | medium |
| file exploration | none (I supply the context) | **yes** — the LLM greps the wiki for itself | yes |
| conversational | no | **very well** | yes |
| cost tracking | the API `usage` field | the `--output-format json` envelope | the SDK |
| dependency | `anthropic` | needs the Claude Code CLI installed | `claude-agent-sdk` |
| precedent | — | `jmsung/einstein` uses this way | — |

### Suitability per stage

The LLM only touches six places (`docs/01` §1).

| Stage | Character | Suits |
|---|---|---|
| S1 INTAKE | multimodal, single call, fixed schema | **(a)** |
| S2 ELICIT | a back-and-forth conversation with a human | **(b)** |
| S4 LIT-GROUND | wiki search + comparison | (b), or (a) + a search tool |
| S5 PLAN | 1 proposal, fixed schema | **(a)** |
| S9 REPAIR triage | 1 diagnosis | **(a)** |
| S12 REPORT | generating a narrative | **(a)** |

> ✅ **Settled · `D23` — (a) the SDK directly** (2026-07-27, during the S1 implementation)
> The output schema is **structurally enforced** with tool-use. Calls happen from the single place `agent/llm.py`
> (enforced by `tests/test_invariants.py`). Only S2 ELICIT, which needs conversational round trips, is left room to be moved out to (b) later.
>
> **It was overturned to (b) once and came back.** The environment had neither `anthropic` nor an API key, only the
> `claude` CLI, so *"(b) works today"* was the judgment — but on actually calling it:
> `Your organization has disabled Claude subscription access for Claude Code · Use an Anthropic API key instead`.
> **A key is needed either way, so (b)'s only advantage disappeared.**
>
> The lesson — **do not decide from an environment survey alone. You have to call it once.** "It is installed" and
> "it works" were different propositions, and the cost of distinguishing them was one call. Deferring `D23` to
> "when it is actually implemented" was right.
> How it went: [`findings/d23-sdk-backend.md`](../../knowledge/wiki/findings/d23-sdk-backend.md)

---

## 5. The pipeline — 14 stages

A state machine on top of `run_state.yaml`. It is not free conversation. Details in
[`docs/01_agent_architecture.md`](../../docs/history/2026-07_bd_agent_01_agent_architecture.ko.md) §2,
and the implementation checklist is **[§14](#14-the-implementation-checklist)**.

### Who starts it — a human

**There is no outer autonomous loop in v1.** No queue, no scheduler, and no layer that "decides for itself what to do next".

```
human:  bd-agent new "<a natural-language description>" [--image ...] [--pdf ...]
   → outputs/<run_id>/ created → proceeds from S1 → stops at a gate
   → the human approves → bd-agent resume outputs/<run_id> → continues
```

In a domain with a weak verification oracle (§6), an autonomous loop carries a large risk of **running all night in
the wrong direction**. With no grader there is also no way for it to learn it was wrong. The loop gets attached after
the gates actually become tedious (`D25`).

> ❓ **Undecided · `D25`** — what is the unit of the "cycle" that one row of `cycle-log.md` counts?
> **The current default:** **one human-started run = 1 cycle.** `run_id` ↔ one log row, 1:1
> **Why it catches:** if one system is run three times changing only a parameter, is that 1 cycle or 3?
> §10's "first-try pass rate" changes completely on this
> **When it gets decided:** when the first row is actually written (at the end of M1)

```
S1 INTAKE → S2 ELICIT →🚦gate 1→ S2.5 PREREGISTER(v0 qualitative)
   → S3 NONDIM → S4 LIT-GROUND → S2.5′ PREREGISTER(v1 quantitative)
   → S5 PLAN → S6 PREFLIGHT
   → S7 EXECUTE(smoke → pilot →🚦gate 2→ production)
   → S7.5 EYEBALL → S8 DIAGNOSE ⇄ S9 REPAIR
   → S10 ANALYZE → S11 VISUALIZE → S12 REPORT
                                        ↑ S8=REDESIGN then back to S5
```

| Stage | Owner | The core pass criterion |
|---|---|---|
| S1 INTAKE | LLM | **0 arbitrarily generated values** — if unknown, into `unknowns[]` · **classify the `(system, target dynamics)` pair** |
| S2 ELICIT | LLM+human | 🚦 **gate 1** — approval of the **back-translation comparison** (§6-A V1) |
| **S2.5 PREREGISTER** | human+LLM | nail down the expected result and its **physical grounds**. Qualitative (v0) → quantitative (v1) |
| S3 NONDIM | code | **based on the system-dynamics card** · round-trip error `<1e-12` · passing the gates the card specifies |
| S4 LIT-GROUND | LLM | judge whether the dimensionless numbers are inside or outside the literature range (a warning, not a block) |
| S5 PLAN | LLM proposes → code validates | 3 tiers defined + **the estimated time < the budget** |
| S6 PREFLIGHT | code | comparison against the real API + a 0-step dry run + **the V4 cross-check**. **Collect all errors** |
| S7 EXECUTE | Runner | 🚦 **gate 2** — immediately before production |
| **S7.5 EYEBALL** | code → human | 3 low-resolution snapshots. The things that are **hard to catch by number and take 1 second by eye** |
| S8 DIAGNOSE | code (+LLM triage) | **6 categories** (§6-B) — the sixth is outlier |
| S9 REPAIR | the rule table → LLM → human | within budget |
| S10 ANALYZE | code (freud) | **no number without an error bar** |

### ★ The non-dimensionalization and the gates differ per `(system, target dynamics)` pair

**One non-dimensionalization convention cannot be forced on every run.** Even inside the same lab the reference units split three ways.

| System · target dynamics | Reference length | Reference time | `kT` |
|---|---|---|---|
| ABP × a dense collective | **the run length `ℓ`** | **`τ_r = 1/D_r`** | **a derived quantity** (`= D_r ζ_r`) |
| brush colloids × non-equilibrium contact | `σ` | `τ_D = σ²/D` | an input |
| passive spheres × equilibrium structure | `σ` | `τ_D` | an input |

So **S1 classifies the pair, and S3 looks up that pair's card and non-dimensionalizes.**
The cards: [`knowledge/wiki/systems/`](../../knowledge/wiki/systems/_index.md) · the contract: [`knowledge/wiki/CLAUDE.md`](../../knowledge/wiki/CLAUDE.md)

**The gates are also switched on and off per pair** — this is the card's core utility.

| Gate | passive spheres × equilibrium structure | ABP × a dense collective |
|---|---|---|
| equilibration (`pymbar`) | ✅ valid | ❌ **meaningless** — an active system never reaches thermal equilibrium |
| self-consistency `D_msd = kT/γ` | ✅ holds | ⚠️ **does not hold** — `D_eff = D_t + U₀²τ_r/2` |
| advective displacement `u₀Δt/σ` | not applicable | ✅ mandatory |

> Proceed without a card and the non-dimensionalization gets done ad hoc, and that is
> **a recurrence of the "parameter selection is tacit knowledge" §2-a pointed out.**

> ⚠️ **S3's `Δt` gate was replaced on 2026-07-27.** It was originally `Δt/τ_D ≤ 1e-4`, and it was confirmed to
> **reject 2 of 3** BD simulations that actually ran and even produced papers — the preceding slit project (`1.0e-3`) and
> the lab's public code `graybox_abp_mpc` (`1.67e-4`). Measured as displacement per step, all three fall within 0.006–0.045σ.
> Grounds: [`findings/dt-gate-should-be-displacement-based.md`](../../knowledge/wiki/findings/dt-gate-should-be-displacement-based.md)
>
> **The threshold `0.03σ` is provisional, coming from 3 samples.** It has to be settled by measuring a `dt` sweep.
> `Δt/τ_D` drops out of the gate but **stays recorded in the dimensionless ledger** — it is needed for literature comparison.
| S11 VISUALIZE | code | the headless render succeeds |
| S12 REPORT | LLM narrative + code assembly | **a comparison table of the pre-registered hypothesis vs the measurement** |

### S2.5 PREREGISTER — pre-registration as an independent stage (new 2026-07-27)

It was originally an appendage of S1 (`hypothesis.yaml`). **It is promoted to an independent stage** — because
nailing down the expected result is not an ancillary task but the only device by which this pipeline prevents post-hoc interpretation.

**It is written twice.** Written once, there is a dilemma over "when to write it" — written before
non-dimensionalization the grounds are weak, and written after it is no longer a pre-registration because the system has already been examined.

| | When | Content | Where the grounds come from |
|---|---|---|---|
| **v0** | right after gate 1 | **qualitative** — which phase? What trend? What would be surprising? | intuition, experience, similar systems |
| **v1** | right after S4 | **quantitative** — numbers + tolerance ranges (the position of `g(r)`'s first peak, the MSD slope, `Z(φ)`) | dimensionless numbers + the literature |

**Both are kept.** If the prediction changed v0 → v1, **that change is itself recorded** — it is evidence that
the non-dimensionalization and the literature actually corrected my intuition, and if it did not change that is
information too. The S12 report carries **v1 vs the measurement** as a comparison table and leaves v0 in an appendix.

> A large difference from the expectation is itself **a signal to investigate**. Pre-registration's second role is as a verification oracle.

### S7.5 EYEBALL — visualization before analysis (new 2026-07-27)

The original order was `S10 ANALYZE → S11 VISUALIZE`. **It is inverted.** There are things that are caught in 1 second
by looking at the trajectory, and those are mostly **hard to catch by number.**

| 1 second by eye | By number |
|---|---|
| it crystallized and is being analysed as a glass | easy to miss looking only at `g(r)` |
| a cluster clumps across the box | you only know by running a separate finite-size check |
| the particles are overlapping (the potential is too soft) | you only know by looking at the minimum-distance histogram |
| a phase separation is in progress and equilibrium is assumed | `τ_ac` merely grows longer; no FAIL appears |

S7.5 is **a cheap and fast eyeball check** — 3 snapshots (initial, middle, final), low resolution, a few seconds to
render. The result feeds into S8 DIAGNOSE as an input. **The final render for the report is still S11**, and there
fresnel (3D) or matplotlib (2D) is used (`D14`).

---

## 6. The verification philosophy ★ — how to verify in a domain with no right answer

Many computational fields have **a grader**. An answer is right or wrong, a bound holds or does not, a score
reproduces. In such a domain "verification" is just "asking the grader".

**We have no grader.** There is no right answer, and no scalar score to maximize. As §2-b said,
a simulation **always produces some number** — not having diverged does not mean being correct.

So the verification **has to be assembled by hand.**

### Verification has two axes — process and result

The original §6 dealt only with **result verification** — "is the number that came out right". But in real work what
goes wrong more often is **the process**. If a different system was simulated to begin with, that number is not even worth verifying.

| Axis | What it asks | Where |
|---|---|---|
| **process verification** | did the pipeline run properly — transcription, validity, scale, consistency, stability, outliers | **§6-A** (V1~V6) |
| **result verification** | is the number that came out actually right | **§6-B** (V7 = the four layers of evidence) |

The two **fail in different ways.** The process fails *quietly* — it simulates the wrong system perfectly — and
the result fails *plausibly*. The former is far more dangerous. **Because it is not caught however well the results
are verified.** Even if all four layers of evidence agree, if it was a different system to begin with then it is a wrong answer on which everything agrees.

---

## 6-A. Process verification — a 7-layer ladder

Each layer catches **a different kind of failure**. They have to be passed in order, and verifying an upper layer
while a lower one is broken is meaningless.

| | Verification | What it asks | Stage | On failure |
|---|---|---|---|---|
| **V1** | **fidelity** | were the words, pictures and text transcribed properly | S1 → 🚦gate 1 | re-run S1 |
| **V2** | physical validity | does this problem make physical sense | S1 · S2.5 | `BLOCKED_INPUT` |
| **V3** | non-dimensionalization, parameters | are the scales and the parameters appropriate | S3 · S4 · S5 | back to S2 / a warning |
| **V4** | **inter-stage consistency** | is there anything conflicting with an earlier stage | S6 PREFLIGHT | back to S5 |
| **V5** | numerical stability | does it not blow up while running | S7 · S8 | S9 REPAIR |
| **V6** | **Outlier** | are there outliers inside the result | S8 | S9 / investigate |
| **V7** | physical interpretation | is the number that came out actually right | S10 · S12 | **to §6-B** |

The three in bold, **V1 · V4 · V6**, are the ones the original plan did not have. The rest bundle into layers what
was scattered across §14-E as exit checks. The implementation items are in §14-E.

> These seven layers are **a 1:1 transposition** of the 7 verification steps described on 2026-07-27. The order is unchanged too.
> The eighth thing mentioned, *"analysing and recording the reason for uncertainties and for what went wrong"*, is not a
> layer but **a discipline that attaches to every layer**, so it was split off into §6-C.

### V1 fidelity — checked by back-translation

The original S1's check was only "**0 arbitrarily generated values**". That is a device against *hallucination*, not
a device for seeing **whether it was transcribed properly**. A failure like this passes straight through.

> Human: "500 nm **silica**" → agent: `material: polystyrene`, `confidence: 0.9`, `unknowns: []`
>
> **It is not an invented value. `unknowns` is empty too. Every current check passes.**

So **back-translation** goes in.

```
the original (human) ──S1──▶ spec.yaml ──back-translate──▶ a natural-language re-description
                                              │
                        🚦gate 1: put it side by side with the original and let a human approve
```

Gate 1's identity is made concrete, from "spec approval" to **"back-translation comparison approval"**. What the human
sees has to be not YAML but **their own words come back as sentences** — far easier than skimming YAML looking for an
omission, and something a human is actually good at.

| Check | Owner |
|---|---|
| are all the values stated in the original present in the spec — **0 omissions** | code |
| if the spec has a value the original does not, `assumed: true` + a source is mandatory | code |
| for a value read from an image, **where** it was read from (`provenance: "the scale bar, 1 µm"`) | code |
| **does the back-translation mean the same as the original** | **the human** (gate 1) |

### V4 inter-stage consistency — the cross-check

V1~V3 and V5~V6 are all **within-stage** checks. Only V4 is different — it looks at whether a later stage
contradicts what an earlier one fixed. **The case where each stage individually passes and yet the combination is wrong** is real.

S6 PREFLIGHT is extended to check them all at once immediately before execution.

| Check | Earlier | Later | Example conflict |
|---|---|---|---|
| dimensionality | S1 `spec.dim` | S5 box | a 3D spec with a 2D box |
| the observation goal vs the run length | S1's goal | S5 `steps` | "observe gelation" but shorter than `τ_gel` |
| the pre-registration vs the plan | S2.5 v1 | S5 | "crystallization expected" but `φ` is outside the crystallization region |
| the time step | S3 `Δt/τ_D` | S5 `dt` | S5 changes `dt` and **invalidates the S3 gate** |
| the box vs the cutoff | S1 `N`, `φ` → `L` | S5 `r_cut` | `r_cut > L/2` — a minimum-image violation |
| the potential vs the kind of system | S1's system | S5 pair | it was said to be hard spheres but only an attraction is applied |
| the budget vs the tiers | `budget` | S5's 3 tiers | the production estimate exceeds the remaining budget |

> **This table getting longer is good news.** Every row added means it is a combination that actually bit once
> (the logic of `D21`), and that it will not bite again.

### V6 Outlier — four kinds

S8's existing diagnostics all look at **the whole system's aggregates** (the mean `kT`, the overall MSD…). An outlier
is at a different level — **it disappears on aggregation.**

| Kind | What | Method |
|---|---|---|
| **ensemble** | 1 of 5 seeds gives a different result | the seed-to-seed scatter vs the block error |
| time series | an energy or pressure jump at a particular frame | a MAD-based robust z-score |
| particle | just one particle with an anomalous displacement (a suspected nlist miss) | the tail of the displacement distribution |
| spatial | only part of the box has a different density | density by volume partition |

**The ensemble outlier matters most.** It bears directly on reproducibility, and the other three are usually its *cause*.
If one seed jumps and it is passed over on the strength of the mean alone, **that mean is not physics but an averaged accident.**

---

## 6-B. Result verification — the four layers of evidence

**V7.** From here on is the original §6. It is only meaningful for a result that has passed the six layers above.

### The four layers of evidence — they have to be **different in kind** from each other

| Layer | What it asks | Example |
|---|---|---|
| **① self-consistency** | does the simulation not contradict itself | a dilute tracer's `D_msd` = `kT/γ` (±2%) · the measured `kT` = the target `kT` (±5%) |
| **② the analytic limit** | going to a limit whose answer is known, does that answer come out | an ideal gas at `φ→0` · `MSD = 6Dt` in free BD |
| **③ the literature benchmark** | does a value others measured reproduce | Carnahan–Starling `Z(φ)` · hard-sphere `φ_freeze=0.494` |
| **④ an independent method** | does the same answer come out by a different route | BD vs HPMC · the same `B2*` with different potentials |

> **Running the same code twice is one piece of evidence.** A re-run with only the seed changed is the same.
> Evidence has to differ in **kind**.

### A fifth-layer candidate — experimental data

Since the research topic includes **microrheology/tracking** (§7), measured trajectories may be to hand. This
belongs to none of the four layers above — it is neither a literature value nor an analytic limit. **In a domain with
no grader, an independent measured oracle is the most valuable evidence there is.**

| Layer | What it asks | Example |
|---|---|---|
| **⑤ experimental comparison** | does it match a value measured in a real system | a trackpy MSD vs the simulated MSD · a measured `g(r)` |

> ❓ **Undecided · `D24`** — should ⑤ be admitted as a formal evidence layer?
> **The current default:** **not admitted in v1.** But `evidence_layer: 5` is reserved in the `benchmarks.yaml` schema
> **Why it is deferred:** an experiment-simulation comparison has too many candidate causes when it disagrees
> (ignoring HI (`D11`) · polydispersity (`D10`) · tracking error · the system simply being different). To become an
> evidence layer there first has to be **a rule for interpreting a disagreement**
> **When it gets decided:** when a microrheology system is first handled (after M2)

### The disagreement protocol

| Situation | Action |
|---|---|
| two pieces of evidence disagree | **it is a bug. Stop.** Do not proceed until it is found |
| all three agree | trust it. Proceed |
| all three disagree | the definitions differ. **Re-check the non-dimensionalization and the observable definitions first** |

### Evidence grades — attach a badge to every result

| Grade | Condition |
|---|---|
| `certified` | ③ the literature benchmark + ② the analytic limit + ④ an independent method |
| `verified` | ① self-consistency + one of ②~④ |
| `plausible` | ① self-consistency only |
| `unverified` | there is nothing to compare against |

**Do not hide an `unverified`.** Simulating a new system that is in no benchmark is legitimate, but calling it
"verified" is not. It is stated in the report as a badge.

> The table above is the **acquisition condition**, and one more **ceiling** applies on top of it — a
> `PASS-with-doubt` coming out of process verification (§6-A) lowers the ceiling of the grade. Details in §6-C. In short,
> **even with all four layers of evidence collected, if the system itself is doubtful it cannot be `certified`.**

### The literature = verification infrastructure

**The literature stands in** for the grader role we do not have. So collecting papers is not an ancillary task but
**the verification infrastructure itself**. Using a paper only as *reading context* throws away half its value —
only by **extracting it into a machine-readable benchmark table and running it as a regression test** does it become a grader.

```yaml
# knowledge/wiki/benchmarks/benchmarks.yaml
- id: carnahan_starling_hs_eos
  system: hard_sphere_3d
  input: {phi: 0.30}
  observable: compressibility_factor_Z
  expected: 4.577            # Z = (1+φ+φ²−φ³)/(1−φ)³
  tolerance_rel: 0.02
  cost: cheap                # 10 cores, tens of seconds
  evidence_layer: 3
  source: "Carnahan & Starling, J. Chem. Phys. 51, 635 (1969)"
```

`pytest tests/test_benchmarks.py` **actually runs a short simulation and compares against the literature value.**
It is the only trustworthy way of confirming that the pipeline is physically right.

**The v1 seed benchmarks** (what is affordable on 10 cores first)

| Priority | Benchmark | Cost |
|---|---|---|
| 1 | free BD `MSD = 6Dt`, `D = kT/γ` — within 1% | a few seconds |
| 2 | **Carnahan–Starling** `Z(φ)`, φ=0.1–0.4 — cheap and highly discriminating | tens of seconds |
| 3 | hard sphere `φ_freeze=0.494` / `φ_melt=0.545` (Hoover–Ree 1968) | minutes |
| 4 | the LJ equation of state — the Johnson–Zollweg–Gubbins table state points | minutes |
| 5 | Noro–Frenkel `B2*` — the analytic `B2` vs a numerical integration | seconds |
| 6 | HPMC hard-disk acceptance + EOS | v2 |

### S8's 6 diagnostic categories

| Category | Indicator | Layer |
|---|---|---|
| stability | no NaN or inf · the maximum displacement per step `< 0.1σ` · 0 box escapes | V5 |
| thermodynamics | the measured `kT` vs the target (±5%) · the pressure · the potential-energy drift | V5 |
| self-consistency | a dilute tracer's `D_msd` vs `kT/γ` (±2%) | V7 ① |
| equilibration | `pymbar.detect_equilibration` · `τ_ac` · `N_eff = N/(2τ_ac)` | V5 |
| finite size | the shift in an observable between `L` and `1.5L` is within tol | V5 |
| **Outlier** | the 4 kinds: ensemble · time series · particle · spatial (§6-A) | **V6** |

The first five look at **aggregates**, and only the sixth looks at **the distribution and the individuals**. All five
can pass and yet one seed can be off on its own, and that gets buried in the mean.

---

## 6-C. The uncertainty ledger — what passed but is doubtful

What comes out of verification is not two things but **three**.

| Verdict | Meaning | Where it is left |
|---|---|---|
| `PASS` | inside the criterion. Move on | the journal (automatic) |
| **`PASS-with-doubt`** | **it passed the criterion but is doubtful** | the journal + **the verification ledger** |
| `FAIL` | outside the criterion | the journal + S9 REPAIR |

**The middle one is the crux.** With only `PASS`/`FAIL`, anything doubtful has nowhere to go and disappears —
it passed, so nobody looks again, and when the result later seems odd there is no way to retrace *"what was it that caught"*.

Cases for attaching `PASS-with-doubt`:

- it passed the criterion **by a hair** (`1.9%` against a `±2%` criterion)
- the grounds for the criterion itself are weak (the threshold came from practice rather than the literature)
- it passed but **differs from the expectation** — a discrepancy with pre-registration v1 is a doubt even on a pass
- the check was **skipped** (because the data is absent, or for cost)

### Where it is recorded — it differs by level

```
every verdict              → decision_journal.jsonl   (automatic. rule_id + observation + threshold + verdict)
PASS-with-doubt            → docs/agent/verification-ledger.md   (a human-readable summary. Append-only)
an unresolved FAIL         → knowledge/wiki/findings/dead-end-<slug>.md   (a cause analysis)
a recurring doubt/FAIL pattern → docs/agent/wall-ledger.md   (grep it before a new attempt)
```

### The reason has to be the cause, not the symptom

§7's `dead-end` discipline is **extended to verification generally.** Write only the symptom in the verdict's reason
and it is useless next time.

| ✗ symptom | ✓ cause |
|---|---|
| "it diverged" | "WCA's `r⁻¹³` core made `F·dt/γ` explode in the overdamped case" |
| "seed 3 is anomalous" | "only seed 3 had overlap left in the initial arrangement, giving an energy spike in the first 100 steps" |
| "`kT` does not match" | "`γ` was entered in SI and the rest in reduced, so the units got mixed" |

If the cause cannot be written, **write that fact itself** — `cause: unknown` + the observed values + the next
investigation candidates. Three `cause: unknown` entries accumulating is a signal that it is a pattern worth investigating.

### A doubt cuts the evidence grade

With no incentive to record, nobody records. So **they are connected.**

| Condition | Maximum grade |
|---|---|
| there is at least one `PASS-with-doubt` | **`certified` impossible** — `verified` is the ceiling |
| a doubt in V1 (fidelity) or V4 (consistency) | **`plausible` is the ceiling** — if the system itself is doubtful the numbers are meaningless |
| there is an unresolved `FAIL` | **`unverified`** |

> For this table to be **an incentive to write a doubt down honestly rather than to hide it**, a low grade has to be
> treated not as a punishment but as **an accurate description of the state**. `plausible` is not a bad result but an honest one.
> As written in §12 — **do not hide an `unverified`.**

---

## 7. Knowledge compounding

### The pipeline

```
raw/        gitignored. PDF·original-code local cache
  ↓  distillation (human-approved)
source/     in-git. 1 .md per original. The provenance is in the frontmatter
  ├── papers/  ★ papers and arXiv        ← v1's only active seed
  └── lab/       a senior's lab simulations   ← assets not obtained. Inactive (D20)
  ↓  synthesis
wiki/
  ├── concepts/     WHAT-IS  — dimensionless numbers, phase behaviour, kinds of potential
  ├── techniques/   HOW-TO   — the equilibration verdict, error bars, initial placement, rendering
  ├── systems/      an index per system (settled as the 4 below)
  ├── benchmarks/   ★ benchmarks.yaml + a grounds page per entry
  ├── findings/     Q→A + citations.  Includes dead-end-<slug>.md
  └── questions/    what has no answer yet. Not deleted but closed with a status
```

### `systems/` — settled as 4

Matched 1:1 to the topics of the paper bundle in hand. The folders and empty index files are made in M0, and
**the content is filled in when that system is actually handled.**

| File | System | Representative verification candidate |
|---|---|---|
| `depletion-gel.md` | depletion attraction · gelation · arrested phase separation | the Lu–Zaccarelli gel boundary · the `g(r)` contact peak |
| `charged-dlvo.md` | charged colloids · Yukawa · the screening length | the Robbins–Kremer–Grest phase diagram |
| `microrheology.md` | MSD · GSER · the viscoelastic modulus | **measured tracking trajectories** (`D24`) |
| `dense-glass.md` | dense systems · the glass transition · crystallization · ψ₆ | Carnahan–Starling · `φ_freeze`=0.494 |

### Paper distillation — v1's only seed

`source/papers/` is the only knowledge asset currently in hand. The procedure is nailed down so it can be executed.

```
knowledge/raw/papers/<year>-<author>-<slug>.pdf          (gitignored)
  ↓  1 per paper. Human-approved
knowledge/source/papers/<year>-<author>-<slug>.md
     required frontmatter: doi · year · system · extracted_benchmarks[] · extracted_parameters[]
     body: only what we can use
  ↓  only when a benchmark was extracted
knowledge/wiki/benchmarks/<id>.md  +  benchmarks.yaml entry
```

**"Only what we can use" is the crux.** The output of a distillation is not a paper summary but these three:
① a reproducible number ② a parameter value and its grounds ③ the clues for carrying it over to our system. If none
of the three comes out, it is not yet time to distil that paper.

> ❓ **Undecided · `D26`** — how many papers to distil to begin with, and in what order?
> **The current default:** before M1, **only 2~3 on free BD / Stokes–Einstein.** The rest as and when the relevant
> system comes up
> **Why not fill it in advance:** filling the wiki in advance is the same failure mode as writing rules in advance (`D21`).
> **A distillation that is never used is never verified** — it sits in the wiki wrong and contaminates the next run
> **When it gets decided:** at the start of M1

### The wiki bootstrap order

In an empty wiki, `wiki-first` is meaningless. What to start with is fixed in advance.

| Step | When | What gets filled |
|---|---|---|
| 1 | M0~M1 | 2~3 paper distillations · 2~3 free-BD-related `concepts` · 1~2 `benchmarks` |
| 2 | during M1 | **only what was actually experienced**, as `findings`. Dead ends included |
| 3 | M2+ | filling in `systems/` as the systems widen |

Step 1's goal is not "a sufficient wiki" but **getting `wiki-first` to hit for real, even once**.
A rule that never once hits is soon ignored.

### A senior's lab simulations — inactive (revived when assets are obtained)

**A paper is "a published result", and a senior's code is "a parameter set that actually ran".** Papers often do not
state `dt` or the number of equilibration steps, so as a parameter prior the latter is far more valuable. In exchange,
**unverified practice comes along with it.**

**No senior's code is currently obtained.** The frontmatter schema and the discipline are settled in `D20`, so when
assets actually appear this section gets revived then. On obtaining them, add `source/lab/` to `.gitignore` immediately (§4, the publication boundary).

> Only one piece of discipline is left here: **do not cite a parameter marked `reproduced: no` as if it were literature grounds.**
> Until reproduced it is a factual record that "this is what was done", not grounds that "this is right". Let that
> distinction collapse and the wiki becomes a rumour store rather than verification infrastructure.

### Authorship symmetry

Both humans and the agent may write the wiki, but it is marked honestly in the frontmatter.

```yaml
type: concept | technique | finding | question | system | benchmark
author: agent | human | hybrid
drafted: YYYY-MM-DD
confirmed_by: human          # optional. After a human has reviewed it
cites: [paths]
```

**The `author: agent` ratio is itself a self-improvement indicator.**
**Promotion (finding → concept) is always human-approved** — it prevents the agent from inflating concepts on its own.

### What causes a wiki page to be created

| Trigger | Result |
|---|---|
| a question arose and the wiki has no answer | `questions/<date>-<slug>.md` (status: open) |
| the question was answered | `findings/<slug>.md` + the original question set to `status: answered` |
| **an approach failed** | `findings/dead-end-<slug>.md` — write down **why it did not work** |
| a benchmark was added | `benchmarks/<id>.md` + a `benchmarks.yaml` entry |
| a finding was cited 3 or more times | promoted to a concept (**human-approved**) |

### Failures are findings

Without a `dead-end` page you hit the same wall again. The minimal form:

```yaml
type: finding
subtype: dead-end
author: agent
drafted: YYYY-MM-DD
system: <in which system>
what_was_tried: <what>
why_it_failed: <why it failed — the cause, not the symptom>
evidence: <the measured values, the log path>
what_to_try_instead: <the next candidates>
```

If `why_it_failed` is "it diverged", that is a symptom and not a cause. **"WCA's `r⁻¹³` core made `F·dt/γ` explode in
the overdamped case"** is the cause.

---

## 8. Compute routing — local first, cluster-extensible

**Decision `D2` (settled 2026-07-27):** simulation and visualization run locally. The seam for cluster extension is kept.

### The current resources

| | |
|---|---|
| machine | Apple M4 · 10 CPU cores · 16GB unified memory |
| HOOMD | 7.1.0 `cpu_py312` — `gpu_enabled=False`, `mpi_enabled=False` |
| implication | **no CUDA.** N~10³–10⁴ and medium-length runs are realistic. The cost gate is mandatory |

### The workload routing matrix

| Workload | Character | v1 (local) | v2 (cluster) |
|---|---|---|---|
| smoke (N≤500, ≤10⁴ steps) | seconds | always local | local |
| pilot (N~2k, ~10⁵ steps) | minutes | local | local |
| production 3D (N~10⁴, ≥10⁷ steps) | hours to days | at the local limit | **the cluster's first priority** |
| a parameter sweep (K independent runs) | embarrassingly parallel | as many as there are cores | **the cluster wins overwhelmingly** |
| an ensemble (k independent seeds) | parallel | multiprocess | the cluster |
| the finite-size check (`L`, `1.5L`) | 2× the cost | local | local |

### The cost estimate is based on measurement

Nothing is guessed. **It is measured once on the machine and stored.**

```
bdkit/run/machine_profile.yaml    # particle-steps/sec, per potential and per N
```

S5 PLAN's cost gate reads this profile to estimate the wall-clock time and does not run if it exceeds the budget.
**Without this gate the agent casually starts a multi-day job.**

> **The measurement comes before M1** (§0's number 2). Without the profile, `A3` (the cost gate) physically cannot
> work, and then M1's S5 and S7 are meaningless wholesale.

### The Runner interface — the seam

```
Runner.estimate(run_plan) -> CostEstimate
Runner.submit(run_plan)   -> Handle
Runner.poll(handle)       -> Status
Runner.fetch(handle)      -> paths
```

`LocalRunner` (v1) / `SlurmRunner` (v2). The pipeline knows only `Runner` and does not know the backend.

### What has to be observed **now** for cluster extension

The things that get expensive to fix later. Observed from v1 onwards.

1. **No absolute paths** — every path is relative to the run directory
2. **`simulate.py` is self-contained** — it takes one argument, a config path
3. **Do not hardcode the device** — `make_device(spec)` instead of calling `hoomd.device.CPU()` directly
4. **Resumable from a checkpoint** — mandatory because of SLURM time limits
5. **Results communicate through files only** — no stdout parsing

---

## 9. The autonomy boundary

### Human approval gates — 2, permanent

| Gate | Position | The failure it blocks |
|---|---|---|
| 🚦 **gate 1** | after S2 — approval of the **back-translation comparison** | the failure of **simulating the wrong system** (§6-A V1) |
| 🚦 **gate 2** | S7, immediately before production | the failure of **casually starting a multi-day job** |

No other gates are placed. With many gates the benefit of automation disappears.

### The budget — it must exist

```yaml
max_total_walltime_s: 21600      # 6 hours
max_repair_iterations: 8
max_disk_gb: 20
max_llm_calls: 100
```

The numbers may be arbitrary. **It must not be absent.** On exhaustion, `ESCALATED`.

### The escalation ladder — the human is last

```
① wall-ledger grep  →  ② wiki lookup  →  ③ re-apply the rule table
                    →  ④ LLM triage  →  ⑤ human = the end of the run
```

**Asking the human first is an anti-pattern.** Skip the ledger and the wiki and the compounding breaks —
you end up grinding at a wall that has already been solved.

**On reaching ⑤ the run ends as `ESCALATED`.** With no outer loop (§5), a state of waiting for a human does not
exist. On ending it leaves a summary of **the symptom · everything that was tried · the next candidates**, adds one
row to the `wall-ledger`, and brings the process down. After the human judges, it re-enters with `bd-agent resume`.

### The L0 axis (axioms) — not changeable without human approval

| | |
|---|---|
| **A1** | Every numerical claim is verified by **3 pieces of evidence of different kinds**. If two disagree, that number is fake |
| **A2** | **Do not produce a number without an error bar.** Block averaging + `τ_ac` alongside |
| **A3** | The **cost gate** has to be passed before running. If the estimated time > the budget, running is forbidden |
| **A4** | `bdkit/` does not call an LLM. `grep -r "anthropic\|claude" bdkit/` must come back empty |

### The 4 v1 rules

| Rule | The incident that produced it |
|---|---|
| `axioms.md` | the L0 invariants |
| `deterministic-core.md` | enforcing A4. So as to distinguish a physics error from an LLM error |
| `overdamped-stability.md` | **there is a real incident** — the preceding project's `~/Research/MD_particle/brownian_slit_sim/src/forces.py:117`: WCA's `r⁻¹³` core flung a particle out of the box in the overdamped case |
| `verify-against-literature.md` | §6's 4-layer evidence system |

**Rule candidates** (written after actually being experienced): `wall-hit-escalation` · `cycle-discipline` ·
`failure-is-a-finding` · `wiki-first-lookup` · `cost-gate` · `ask-the-question-first` ·
`error-bars-or-silence` · `compute-router`

### Hooks — v2

They are raised in the order `rule → hook(warn) → hook(block)`. They get attached
**after confirming a rule actually was not observed** — discipline that is observed merely by being written in a
document is rare, and discipline that actually is observed usually has a deterministic enforcement device attached.

Candidates: warn if a heavy simulation command runs 3 times in a row without a wiki lookup (`PostToolUse(Bash)`) ·
block termination if a cycle leaves no log row + finding (`Stop`).

One principle to observe when attaching a gate:
**pass on an infrastructure error (fail-open).** A gate must not lock a session in.

---

## 10. Success indicators — turning "self-improvement" into a verifiable claim

Without indicators, "self-improvement" is a slogan.

### The operating policy — record everything, analyse later

Since the pipeline is started by a human (§5), cycles accumulate slowly. **The trend of the indicators over the
first dozen or so cycles is almost noise.** So the two are separated.

| | When | What |
|---|---|---|
| **recording** | from now on, every cycle | all 7 below as columns of `cycle-log.md`. The computational cost is near 0 |
| **analysis** | **after 20 cycles** | trends, correlations, the dashboard. Draw a graph before that and noise gets read as signal |

What is not recorded now cannot be filled in retroactively. Conversely, analysing now has nothing to gain.
`docs/agent/metrics.md` is created at the end of M1.

| Indicator | Definition | What it asks |
|---|---|---|
| **the first-try pass rate** | the fraction of runs that passed S6 preflight without modification | is parameter selection improving |
| **the mean repair iterations** | the number of S9 iterations to a PASS | is it learning from failures |
| **the benchmark pass rate** | the fraction of `benchmarks.yaml` passing | is the physics right |
| **the wiki reuse rate** | the fraction of new runs that cite an existing wiki page | is the knowledge actually used |
| **the `author: agent` ratio** | the share of wiki pages authored by the agent | is the agent contributing |
| **the number of human interventions** | the number of escalations outside the 2 gates | is autonomy increasing |
| **the evidence-grade distribution** | the `certified` / `verified` / `plausible` / `unverified` ratios | is the verification level rising |

> **Record the changes that move an indicator down too.** Keep only what improved and it is not an indicator but a brochure.

### The cycle log

```
docs/agent/cycle-log.md      exactly 1 row per cycle. Failures included      (created at the end of M1)
docs/agent/wall-ledger.md    append-only where it got stuck. Always grep before a new attempt   (created at the end of M1)
```

The column definitions — the 7 indicators above become the columns as they are. Nailed down now so as not to
redesign them when the first row is written.

```
| # | date | run_id | system | tier | first-try | repair | benchmark | wiki cites | agent-authored | escalation | evidence grade | note |
```

| Column | Value | Notes |
|---|---|---|
| `first-try` | `PASS` / `FAIL` | did it pass S6 preflight without modification |
| `repair` | an integer | the number of S9 iterations to a PASS |
| `benchmark` | in the form `3/4` | the number of `benchmarks.yaml` entries passing among those run in this run |
| `wiki cites` | an integer | the number of existing wiki pages this run cited. **0 means the compounding is not happening** |
| `agent-authored` | an integer | the number of pages newly written as `author: agent` in this run |
| `escalation` | an integer | human interventions **excluding** the 2 gates |
| `evidence grade` | `certified`\|`verified`\|`plausible`\|`unverified` | the final result's badge |

**No cherry-picking. No editing past rows** — a correction is made by adding a `corrigendum:` row.

---

## 11. Milestones

| | Name | Definition of done |
|---|---|---|
| **M0** | the skeleton | the checklist below. **← current** |
| **M1** | the first full run | the checklist below. One free-BD system through S1→S12, **the whole range** (14 stages) |
| **M2** | physics extension | the excluded-volume potential settled (`D8`). **Carnahan–Starling passing**. The `g(r)`, `S(k)` and `ψ₆` observables + error bars |
| **M3** | auto-repair | deliberately inject a failure → the rule table recovers. The `wall-ledger` works. The first `dead-end` page |
| **M4** | compounding | 20+ wiki pages (paper distillations + our own findings). The indicators reach 20 cycles → the first trend analysis |
| v2 | extension | HPMC · the SLURM Runner · hooks · signac |
| v3 | candidates | a multi-perspective reconnaissance panel · comparison against experimental data — see below |

**M1 matters most.** Going around once all the way, however shallowly, is better than going deep and getting cut off in the middle.

This table fixes **when**. **What** gets built is broken into items in [§14 the implementation checklist](#14-the-implementation-checklist),
and which milestone each item belongs to is in the `M` column there.

### The M0 definition-of-done checklist

Do not write all the documents before starting. `docs/02`~`10` get **made as thin stubs first and filled in while
implementing** — a compromise that keeps the whole picture visible without blocking the start of the code.

```
□ docs/02~10 as 1-page stubs each
    (a table of contents + "the 3 questions this document answers" only. The body is left empty)
□ git init (private) + .gitignore(outputs/, knowledge/raw/) + environment.yml
□ conda env `bd_agent` built → hoomd, freud, fresnel, gsd and pymbar import successfully
□ machine_profile.yaml measured (free BD, N = 1k / 5k / 10k)
□ create the repository skeleton:
    bdkit/{spec,units,plan,build,run,diagnose,repair,analyze,viz,report}/
    agent/  knowledge/{raw,source,wiki}/  tests/  outputs/
□ 3 contract documents for the 3 layers: CLAUDE.md · .claude/CLAUDE.md · knowledge/wiki/CLAUDE.md
□ .claude/rules/ 4 files (D21)
```

### The M1 definition-of-done checklist

"One free-BD system run end to end", in a checkable form. **All of it has to be true for it to be M1.**

```
□ bd-agent new "1000 spheres of 500 nm silica floating in water, no interaction" to start
□ S1~S12 all 14 stages pass with no human intervention outside the 2 gates (S2.5 and S7.5 included)
□ pytest tests/test_benchmarks.py::test_free_bd_stokes_einstein passes (within 1% error)
□ report.html — an MSD graph with error bars + an evidence-grade badge, a single self-contained file
□ decision_journal.jsonl — actor-wise (rule/llm/human) aggregation comes out
□ cycle-log.md 1 row + at least 1 finding (a dead end is allowed)
□ resuming the same run_id does not recompute the completed stages
□ grep -r "anthropic\|claude" bdkit/ comes back empty (A4)
```

The last line matters especially. **That `bdkit/` is right on its own without an LLM** is this architecture's only
means of debugging (§4), and once contaminated it is hard to undo.

### Considered but not adopted in v1 — and why

> Silently omitting what was decided against means repeating the same argument in 6 months. It is stated.

| Not adopted | Why |
|---|---|
| **a multi-perspective reconnaissance panel** (dispatching several personas in parallel to write questions from differing angles) | It is a device suited to **divergent reconnaissance of an unsolved problem**. Our v1 task — "simulate a given system correctly" — is **convergent**, so the correct path is largely fixed and the benefit is unclear. **Conditions for revisiting:** when the scope widens to an exploratory task where which observable to look at is itself unclear. → `D22` |
| **writing several rules in advance** | A rule with no incident record becomes a ritual. **v1 has 4.** → `D21` |
| **hooks** | Raised in the order `warn` → `block` after confirming a rule actually was not observed. v2 |
| **signac** | A parameter-space management tool. Brought in when sweeps actually hurt. → `D13` |
| **a Streamlit/Panel dashboard** | The outputs scatter and do not get archived. A single HTML report is sufficient |
| **an outer autonomous loop** (a layer that decides its own next task and runs all night) | In a domain with a weak verification oracle (§6) there is a large risk of **running all night in the wrong direction**. With no grader there is also no way for it to learn it was wrong. **Conditions for revisiting:** when the 2 gates actually become tedious and the benchmarks become dense enough that the automatic verdict is trustworthy. → `D25` |
| **activating the senior's-code layer** | The schema and the discipline are settled in `D20`, but **there are no assets obtained yet.** Spending page space on an asset you do not have blurs what to do first. Revived the moment they are obtained. → `D20` |
| **indicator trend analysis and a dashboard** | Cycles accumulate slowly, so the first dozen or so are noise. **Record everything from now on**, analyse after 20 cycles (§10). Cycles being slow is because there is no outer loop → `D25` |

---

## 12. Explicit non-goals — What this is NOT

> Inheriting the preceding project's *"What this is NOT"* convention. If what is not being done is not written down,
> the same argument gets repeated later.

- **Hydrodynamic interactions (HI) are not done.** HOOMD's BD/Langevin does not include HI.
  It can be quantitatively wrong for sedimentation, shear and dense systems. The approximation is **stated** in the report (`D11`).
- **It is not a tool for discovering new physics.** It is a tool for simulating existing physics
  **correctly and reproducibly**. Discovery happens on top of that, by a human.
- **It does not replace a human's physical judgment.** The 2 gates are permanent. "Fully automatic" is not a goal.
- **It does not call a system that is in no benchmark "verified".** It attaches an `unverified` badge.
- **It does not cite a senior's parameter marked `reproduced: no` as grounds.** It is a factual record, not grounds.
- **It does not treat quantum or electronic structure.** The classical colloidal scale only.
- **It does not build a Streamlit/Panel dashboard in v1.** The outputs scatter and do not get archived.
- **It does not build an outer autonomous loop in v1.** The pipeline is always started by a human (§5).
- **It does not fill the wiki in advance in v1.** A distillation that is never used is never verified (§7 `D26`).

---

## 13. Open questions

The full list and the current defaults are in **[`docs/00_decision_log.md`](../../docs/history/2026-07_bd_agent_00_decision_log.ko.md)**.
In the document body they are planted inside the relevant section as `> ❓ undecided · Dnn` markers — so that they are
visible right where a decision is needed. **Every marker has a `when it gets decided`.** An undecided item with no
date stays undecided forever.

The four heaviest right now:

| ID | Question | Status | When it gets decided |
|---|---|---|---|
| `D3` | the v1 physics scope — whether HPMC is included | `OPEN` — the default is BD only. MC's verification logic differs, so the pipeline becomes two | at the start of M2 |
| ~~`D8`~~ | ~~the default excluded-volume potential (harmonic / WCA / Wang–Frenkel)~~ | **`DECIDED`** (2026-07-28) — **WCA.** As foretold, it was settled by **measuring** against Carnahan–Starling (14 runs). `ε`=1 is within 2% for `φ_eff`≲0.32, and `ε`=10 is +0.56% at 0.38. The distillation is [`wca-reproduces-carnahan-starling`](../../knowledge/wiki/findings/wca-reproduces-carnahan-starling.md).<br>**The limit: of the three candidates only WCA was measured** — it is not *"the best fit"* but "the first one that fits well enough" | ~~M2~~ closed |
| ~~`D14`~~ | ~~the renderer, OVITO vs fresnel~~ | **`DECIDED`** (2026-07-28) — **fresnel + matplotlib.** matplotlib for 2D, fresnel only for 3D. OVITO downgrades `tbb` and endangers `hoomd` | ~~before starting M1's S11~~ closed |
| `D24` | whether to admit experimental data as a fifth evidence layer | `OPEN` — the default: not admitted in v1, the schema reserved only. Adopting it adds one **independent oracle** | at the start of a microrheology system |

`D6` (git) is **`DECIDED`** — settled as a private repository (§0's number 1). `D25` (the cycle unit), `D26` (the
distillation order) and `D27` (the publication separation procedure) are markers in §5·§7·§4 respectively.

---

## 14. The implementation checklist — **engine progress is not here**

> **Tidied 2026-07-28 (`D37`).** This section was originally a 144-item progress checklist and said
> `2 / 144`. That denominator **was false** — it was written for the structure before `D32` (absorbing the engine),
> so its artifact cells assume 10 `bdkit/{spec,units,plan,…}/` subpackages while the actual work went to
> flat `simbot/` modules, and things functionally finished were sitting there as `X`.
>
> With progress recorded in two places they **inevitably diverge.** So the authoritative record of progress was reduced to one.

### Where to look

| Progress of what | Authoritative record | Why there |
|---|---|---|
| **the engine (P1~P7 · 7 stages)** | [`docs/11_simbot.md`](docs/11_simbot.md) **§5 the stage table + §7 the progress record** | for each stage, *what it is trying to confirm, why this order, and the measured result* are together. It carries more information than a checkbox, and §7 being append-only it can be looked back over |
| **benchmarks (literature comparison)** | [`knowledge/wiki/benchmarks/benchmarks.yaml`](../../knowledge/wiki/benchmarks/benchmarks.yaml) | each entry carries its expected value, tolerance and source, and **for what cannot be done the reason is written in `blocked_by`.** A checkbox cannot hold "why did it not work" |
| **decisions** | [`docs/00_decision_log.md`](../../docs/history/2026-07_bd_agent_00_decision_log.ko.md) (`D`) · `docs/11_simbot.md` §7 (`SD`) | |
| **the verification ladder (V1~V7)** | **`14-E` below** | not yet started and existing nowhere else. `D35` (human review → threshold registration → automatic verdict) is the entry path to this layer |

### What was deleted — `14-A` infrastructure · `14-B` per-stage implementation (83 items)

**Deleted does not mean "not done".** Most of it is finished, and the way it finished differed from the table's assumptions.
To check where an item went, look at `docs/11_simbot.md` §5. Three representative mismatches:

| The table's assumption | Reality |
|---|---|
| `bdkit/state.py` updates `run_state.yaml` atomically | `simbot/state.py` does it (`SD24`, §7 #37). `bd-agent resume` is `cli.py resume` |
| `agent/s1_intake` reads natural language and makes the spec | **it was decided not to put an LLM path in the engine** (`SD22`). The session does the reading and the engine only assembles |
| 14 stages (including the 2 gates, S2.5, S8 and S9) | **reduced to 7** in the chatbot pivot (`docs/11_simbot.md` §2). The gates, the diagnosis and the auto-repair are removed layers, so they were a *deletion* rather than an `X` |

The last row is the real reason this section was deleted — **leave a removed layer as an `X` and dozens of
nonexistent backlog items appear.** To revive them, start from `D35` and `14-E`.

### `14-0` the 9-step workflow ↔ pipeline mapping is moved to §5 to be read there

The 9-step correspondence table that used to be here (text and pictures → systematization … the nonsense review)
retains value as **a record of the design intent**. But since the pipeline shrank 14 → 7 stages, the correspondence changed too:

```
1 text·pictures → systematization  →  reading (the session) + compose   SD22 · Step 6
2 expected result + physical grounds  →  (removed — S2.5 pre-registration)       11_simbot §2
3 non-dimensionalization + grounding papers  →  P2 NONDIM               SD9 (the reference differs per system)
4 building the simulation · timescales  →  P3 PLAN + P4 BUILD           SD8 · overdamped-stability
5 execution                →  P5 RUN
6 visualization            →  the figures of P7 REPORT                  the 3 snapshots were kept
7 analysis · comparison with the expectation  →  P6 ANALYZE             it does not judge (opened by D35)
8 the final physical and numerical interpretation  →  the human         11_simbot §2
9 the per-stage instability review  →  cross-cutting → 14-E
```

★ The original judgment that number 9 is not a stage but **a pass condition of every stage** remains valid.

---

### 14-C. Verification infrastructure

> **The status of individual benchmarks is not counted here** — [`benchmarks.yaml`](../../knowledge/wiki/benchmarks/benchmarks.yaml)
> carries the expected value, tolerance and source per entry, and for what cannot be done **the reason it is blocked**
> is written in `blocked_by`. A checkbox cannot hold that reason. Below counts **the infrastructure** only.

| | Item | Artifact | M |
|---|---|---|---|
| `O` | the pytest skeleton — **pure logic that runs without HOOMD** | `tests/` 15 files · 333 passing | M0 |
| `O` | the `benchmarks.yaml` schema + a loader + marker separation (`-m benchmark`) | `knowledge/wiki/benchmarks/` · `tests/test_benchmarks.py` | M1 |
| `O` | take slow benchmarks out of the default suite — **a test that does not run is not a test** | `pyproject.toml` `addopts` | M1 |
| `~` | running the registered benchmarks — free BD and the trap family run, the hard-sphere family is `blocked` | `benchmarks.yaml` 12 entries | M1~M2 |
| `X` | **estimator error regression** — measuring the estimator on synthetic data, at the same grade as a benchmark | `tests/test_simbot_estimators.py` (exists) → register it | M2 |
| `X` | the review-ledger → threshold-registration route (`D35`) | undecided — `D35`'s unresolved item | M2 |
| `–` | state-machine transition tests (mock artifacts, 1:1 with the 14-stage transition table) | — | the stages shrank to 7 so the transition table itself does not exist |
| `–` | journal completeness — an aggregation by `actor` (rule/llm/human) | — | the engine has no LLM so it is all `rule` (`SD22`) |
| `–` | LLM extraction accuracy — comparison against a golden `SystemSpec` | `tests/test_compose.py` | since the session does the reading, **assembly equivalence is tested instead of reading accuracy** (§7 #34) |
| `–` | ⑥ HPMC hard-disk acceptance + EOS | — | v2 (`D3`) |
| `X` | `reproduce <paper-slug>` — extract a spec from a distilled paper and run and compare automatically | — | v2 |

> The last row corresponds to workflow 1's *"later, verifying by simulation while checking archives or the papers of
> well-known researchers"*. **The substance of the `S0` long-term goal is this one row** — extracting a spec from a
> distillation, running the pipeline and comparing against the literature value. It is an automated form of §6's
> ③ literature-benchmark layer, and it only stands with `D35` (threshold registration) and
> `14-E` (the ladder) in front of it.

### 14-D. The knowledge and rule layers

| | Item | Artifact | M |
|---|---|---|---|
| `O` | the `knowledge/{raw,source,wiki}/` folder structure | — | M0 |
| `O` | paper distillation — **41 papers** (`D26`) | `knowledge/source/papers/` | M1 |
| `O` | 6 `systems/` cards | `knowledge/wiki/systems/` | M0 |
| `O` | the first `findings/` — 6 | `knowledge/wiki/findings/` | M1 |
| `O` | benchmark grounds pages | `benchmarks/candidates.md` · `choi2020-interfacial-rdf.md` | M1 |
| `O` | **the first `wiki-first` hit** — the excluded-volume verdict cited the `φ` ledger of the `passive-sphere` card and the CS expression (`SD25`, 2026-07-28) | `simbot/excluded.py` | M1 |
| `X` | `concepts/` — **it is empty.** Candidates: `τ_D` · overdamped · block averaging · `Z(φ)` | `knowledge/wiki/concepts/` | M1 |
| `X` | `techniques/` — it is empty | `knowledge/wiki/techniques/` | M1 |
| `X` | **the parameter-grounds ledger** — which value was chosen and why, machine-readable per run (`D38`) | `runs/<id>/parameters.yaml` | M1 |
| `X` | `cycle-log.md` + the 13-column definition (§10) | `docs/agent/` | M1 |
| `X` | `metrics.md` + the first trend analysis after reaching 20 cycles | `docs/agent/` | M4 |
| `–` | 9 `docs/02`~`10` stubs | — | they were not written and were not needed either. It runs on the three `docs/00`, `01` and `11` |

### 14-E. Cross-cutting — the 7-layer verification ladder

**The implementation of §6-A's V1~V7.** Not a stage but **a pass condition of every stage.** S8 DIAGNOSE runs only
once, *after* execution, whereas physical nonsense is **far more often caught before execution, and far more cheaply.**
**V1 · fidelity** — were the words, pictures and text transcribed properly (S1 → 🚦gate 1) ★new

| | Check | Artifact | M |
|---|---|---|---|
| `X` | **generate the back-translation** — `spec.yaml` → a natural-language re-description | `agent/s1_intake/backtranslate.py` | M1 |
| `X` | the omission check — are all the values stated in the original present in the spec | `bdkit/reading/fidelity.py` | M1 |
| `X` | force `assumed: true` + a source on a value absent from the original | `bdkit/reading/fidelity.py` | M1 |
| `X` | the gate 1 output — **the original vs the back-translation side by side** | `bdkit/cli.py` | M1 |
| `X` | `provenance` mandatory for a value read from an image (a scale bar and so on) | `bdkit/reading/` | M2 |

**V2 · physical validity** — does this problem make sense (S1 · S2.5)

| | Check | Kind | M |
|---|---|---|---|
| `X` | `φ` exceeds the geometric ceiling (RCP `0.64` / FCC `0.74`) | physics | M1 |
| `X` | `T`, `η` or the radius outside the physical range (sign, magnitude) | physics | M1 |
| `X` | outside the colloidal scale (a warning if outside 1 nm ~ 10 µm) | physics | M2 |
| `X` | the system itself is contradictory (calling it hard spheres while specifying only an attraction) | physics | M2 |
| `X` | the S2.5 expectation contradicts the spec (crystallization expected at `φ=0.05`) | physics | M2 |

**V3 · non-dimensionalization · parameters** — are the scales appropriate (S3 · S4 · S5)

| | Check | Kind | M |
|---|---|---|---|
| `X` | the unit round-trip error `≥1e-12` | numerics | M1 |
| `X` | `τ_B/τ_D ≥ 1e-3` → **the overdamped assumption breaks** | physics | M1 |
| `X` | `Δt/τ_D > 1e-4` | numerics | M1 |
| `X` | a dimensionless number outside the literature range (a warning, not a block) | physics | M1 |

**V4 · inter-stage consistency** — does it conflict with an earlier stage (S6 PREFLIGHT) ★new

| | Check (earlier → later) | Kind | M |
|---|---|---|---|
| `X` | dimensionality: S1 `spec.dim` → S5 box | consistency | M1 |
| `X` | the time step: S3 `Δt/τ_D` → S5 `dt` (**preventing gate invalidation**) | consistency | M1 |
| `X` | the box vs the cutoff: S1 `N`·`φ`→`L` → S5 `r_cut > L/2` | numerics | M1 |
| `X` | the budget vs the tiers: `budget` → the S5 production estimate | resources | M1 |
| `X` | the observation goal vs the run length: S1's goal → S5 `steps` | physics | M2 |
| `X` | the pre-registration vs the plan: S2.5 v1 → S5's parameters | physics | M2 |
| `X` | the potential vs the kind of system: S1's system → S5 pair | physics | M2 |

**V5 · numerical stability** — does it not blow up while running (S6 · S7 · S8)

| | Check | Kind | M |
|---|---|---|---|
| `X` | the particles already overlap in the initial arrangement (S6) | numerics | M1 |
| `X` | a HOOMD API signature mismatch (S6) | code | M1 |
| `X` | NaN or inf occurring | numerics | M1 |
| `X` | the maximum displacement per step `≥0.1σ` — **`F·dt/γ` exploding** (`overdamped-stability`) | numerics | M1 |
| `X` | a box escape · the measured `kT` outside ±5% of the target · a PE drift | physics | M1 |
| `X` | by eye (S7.5): crystallization · a cluster crossing the box · overlap · a phase separation in progress | physics | M1 |
| `X` | equilibration not reached — `N_eff` insufficient | statistics | M2 |

**V6 · Outlier** — what disappears on aggregation (S8) ★new

| | Check | Method | M |
|---|---|---|---|
| `X` | **ensemble** — one seed alone gives a different result | the seed-to-seed scatter vs the block error | M2 |
| `X` | time series — an energy or pressure jump at a particular frame | a MAD-based robust z-score | M2 |
| `X` | particle — just one particle with an anomalous displacement (a suspected nlist miss) | the tail of the displacement distribution | M2 |
| `X` | spatial — only part of the box has a different density | density by volume partition | M3 |

**V7 · physical interpretation** — is the number that came out right (S10 · S12) → §6-B

| | Check | Kind | M |
|---|---|---|---|
| `X` | a dilute tracer's `D_msd` vs `kT/γ` outside ±2% (evidence ①) | self-consistency | M1 |
| `X` | **a number was produced without an error bar** (an `A2` violation) | honesty | M1 |
| `X` | judge how many of the 4 evidence layers were obtained → produce the evidence grade | honesty | M1 |
| `X` | **the disagreement protocol** — halt if two pieces of evidence disagree | honesty | M2 |
| `X` | the evidence grade is `unverified` and there is no badge | honesty | M1 |

**Recording** — uncertainties and causes (§6-C) ★new

| | Item | Artifact | M |
|---|---|---|---|
| `X` | **the 3-way verdict** — `PASS` / `PASS-with-doubt` / `FAIL` | `bdkit/verify/verdict.py` | M1 |
| `X` | every verdict → the journal (`rule_id`+`observation`+`threshold`+`verdict`) | `bdkit/journal.py` | M1 |
| `X` | `PASS-with-doubt` → `verification-ledger.md` (append-only) | `docs/agent/` | M1 |
| `X` | **apply the doubt → evidence-grade ceiling** (the §6-C table) | `bdkit/report/` | M1 |
| `X` | **the cause, not the symptom** in the verdict's reason — the `cause` field mandatory | `bdkit/verify/` | M1 |
| `X` | detect 3 accumulated `cause: unknown` → promote to the `wall-ledger` | `bdkit/verify/` | M3 |
| `X` | an unresolved `FAIL` → an automatic `findings/dead-end-<slug>.md` stub | `bdkit/verify/` | M3 |

**Four shared disciplines**

1. **Do not quietly swallow a failure and move to the next stage.** A stage that passed without a check is marked in
   the report with an `UNVERIFIED` badge.
2. **Every verdict into the journal** — `rule_id` + `observation` (the measured value + the threshold) + `verdict`.
   It has to be possible to count afterwards "at which check did it catch".
3. **A verdict criterion has to be a number.** The only rows above with no threshold are V1's back-translation
   comparison and V5's by-eye check, and that is deliberate — because **things only a human can catch are real**.
   Instead, the fact that a human looked, and their judgment, are left in the journal.
4. **Pass what is doubtful, but record it.** Without `PASS-with-doubt`, anything doubtful has nowhere to go and
   disappears, and when the result later seems odd it cannot be retraced (§6-C).

---

## Appendix — references

| | |
|---|---|
| the preceding project | `~/Research/MD_particle/brownian_slit_sim` — the config pattern · the verification ladder · overdamped lore |
| the engines | HOOMD-blue 7.1.0 · freud 3.5.0 · fresnel 0.13.8 · gsd 5.0.1 · pymbar 4.2.0 |
| the agent referenced | [`jmsung/einstein`](https://github.com/jmsung/einstein) — the source of the 3-layer structure, knowledge compounding and triple verification |

### The document index

In M0, `02`~`10` are **made as stubs first** (§11). The statuses below get updated when the stubs are created.

| | Document | Status |
|---|---|---|
| `00` | [the decision log](../../docs/history/2026-07_bd_agent_00_decision_log.ko.md) | **written** |
| `01` | [the agent architecture](../../docs/history/2026-07_bd_agent_01_agent_architecture.ko.md) | **written** |
| `02` | the system spec (the `SystemSpec` schema · the S1/S2 inputs) | stub planned |
| `03` | units and non-dimensionalization (the `τ_D` reference, `UnitMap`) | stub planned |
| `04` | potentials and HOOMD script generation | stub planned |
| `05` | the run plan · the cost gate · preflight · the diagnostic indicators | stub planned |
| `06` | the repair rule table | stub planned |
| `07` | observables and error bars | stub planned |
| `08` | visualization and the report | stub planned |
| `09` | the literature and benchmarks | stub planned |
| `10` | the roadmap | stub planned |
| `11` | [the simulation engine plan](docs/11_simbot.md) | **written** — `simbot/`'s plan, progress record and decisions (`SD1`~`SD11`) |
| — | `docs/agent/metrics.md` · `cycle-log.md` · `wall-ledger.md` | created at the end of M1 |
