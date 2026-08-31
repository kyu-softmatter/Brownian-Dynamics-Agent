# 01 · Agent architecture

> This document defines **what the agent is made of, in what order it does what, and where it stores state**.
> The physics and numerics of individual stages are covered in `02`~`09`.

---

## 0. First: general principles for building an agent

When building an agent for the first time, most failures come not from the physics but from the **structure**.
The 10 below apply independently of the domain, and every design decision in this project derives from them.

| # | Principle | Its implementation in this project |
|---|---|---|
| 1 | **Keep state in files, not in the conversation.** Conversation context disappears and cannot be inspected | `run_state.yaml` |
| 2 | **Do not put an LLM inside a numerical loop.** It is non-deterministic, untestable, and quietly wrong | `bdkit/` has 0 LLM dependencies |
| 3 | **Force structured output at every boundary.** The bugs start the moment you parse free text | every LLM output passes JSON Schema validation |
| 4 | **Journal every decision, down to who made it (actor).** The only means of debugging an agent | `decision_journal.jsonl` |
| 5 | **Fix the budget first.** An infinite loop is an agent's default failure mode | `run_state.budget` (D18) |
| 6 | **Run the small tier first.** The first run must not be the expensive one | smoke → pilot → production |
| 7 | **A verdict criterion has to be a number.** "It looks plausible" is not a gate | `\|D_msd/(kT/γ) − 1\| < 0.02` |
| 8 | **Make it say when it does not know.** An LLM's biggest failure in scientific work is quietly inventing a plausible parameter | `unknowns[]` + `confidence` + `assumed: true` |
| 9 | **Make it possible to go back.** Checkpoints, idempotent stages, resume | per-stage artifacts + `resume` |
| 10 | **Put the human gates immediately before the cost jumps.** Nowhere else | 2 gates (D4) |

> **The one line that matters: the LLM proposes, and deterministic code decides.**

---

## 1. Layer structure

```
┌─────────────────────────────────────────────────────────┐
│ user        text · sketch · photo · recording · paper PDF │
└───────────────────────────┬─────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────┐
│ agent/     LLM layer                                     │
│            extraction · literature check · failure triage · report narrative │
│            → must return structured output (schema) only  │
└───────────────────────────┬─────────────────────────────┘
                            ↓  JSON Schema validation (the boundary)
┌─────────────────────────────────────────────────────────┐
│ bdkit/     the deterministic core ── LLM deps 0, verified by pytest alone │
│            computation · verification · execution · diagnosis · analysis · visualization │
└───────────────────────────┬─────────────────────────────┘
                            ↓
             HOOMD-blue 7.1.0 · freud · OVITO/fresnel
```

### `bdkit/` — the deterministic core

| Module | Responsibility | Related document |
|---|---|---|
| `spec/` | the SystemSpec schema and validation. **Collects all errors** | `02` |
| `units/` | `UnitMap` (SI ↔ reduced), the dimensionless-number ledger. **Reads the system-dynamics card to fix the reference units** | `03` |
| `plan/` | generating and validating a `RunPlan`, cost estimation | `05` |
| `build/` | generating the HOOMD script + static preflight checks | `04`, `05` |
| `run/` | the `Runner` interface · `LocalRunner` (v1) · `SlurmRunner` (v2) | D2 |
| `diagnose/` | stability, thermodynamic, self-consistency, equilibration and finite-size indicators. **Card §7 switches them on and off** | `05` |
| `repair/` | the failure → action rule table | `06` |
| `analyze/` | freud observables + block-averaged error bars | `07` |
| `viz/` | particle rendering + matplotlib plots | `08` |
| `report/` | assembling a single HTML | `08` |

**The invariant: nowhere in `bdkit/` is an LLM called.** `grep -r "anthropic\|claude" bdkit/` must come back empty.

### `agent/` — the LLM layer

Three pieces per stage: `prompt + output schema + validator`.
The stages an LLM touches are only these six: **S1, S2, S4, S5 (proposals only), S9 (when no rule hits), S12**.

**The calling method is `D23` — `DECIDED` (2026-07-27): the Anthropic SDK directly**, from the single place `agent/llm.py`.
The output schema is **structurally enforced** with tool-use (principle 3). Only S2 ELICIT, which needs conversational
round trips, is left room to be moved out to Claude Code later. The comparison table is in [`master_plan.md`](../../docs/history/2026-08_simulation_auto_master_plan.ko.md) §4.5,
and how it was settled is in [`findings/d23-sdk-backend.md`](../../knowledge/wiki/findings/d23-sdk-backend.md).

Of the three pieces, **the prompt, the output schema and the validator are reused independently of the calling method** — change
the method and only `agent/llm.py` is replaced. This seam is kept.

---

## 2. The state machine

**The entry point is a human** (`D25`). There is no outer autonomous loop, queue or scheduler in v1.

```
human:  bd-agent new "<natural-language description>" [--image ...] [--pdf ...]     → creates outputs/<run_id>/
human:  bd-agent resume outputs/<run_id> [--from S5_PLAN]         → resumes after gate approval
```

```
S1 INTAKE → S2 ELICIT →🚦gate 1→ S2.5 PREREGISTER(v0)
    → S3 NONDIM → S4 LIT-GROUND → S2.5′ PREREGISTER(v1) → S5 PLAN → S6 PREFLIGHT
    → S7 EXECUTE(smoke → pilot →🚦gate 2→ production)
    → S7.5 EYEBALL → S8 DIAGNOSE ⇄ S9 REPAIR
    → S10 ANALYZE → S11 VISUALIZE → S12 REPORT → DONE
                                                        ↑ S8=REDESIGN then back to S5
```

| State | Owner | On success → | On failure → |
|---|---|---|---|
| `S1 INTAKE` | LLM | `S2` | `BLOCKED_INPUT` (contradictory input) |
| ↳ | | **classify the `(system, target dynamics)` pair** → look up the card. On a classification failure, `unknowns[]` | |
| `S2 ELICIT` | LLM + human | `GATE_SPEC` | `S1` |
| `🚦GATE_SPEC` | human | `S2.5` | `S2` (revision requested) |
| `S2.5 PREREGISTER` | human + LLM | `S3` | `S2` (the expectation contradicts the spec) |
| `S3 NONDIM` | code | `S4` | `S2` (overdamped is unsuitable → negotiate an engine change) |
| ↳ | | **card-based** — apply `wiki/systems/<system>--<dynamics>.md` §3·§4·§7. If there is no card, create one | |
| `S4 LIT-GROUND` | LLM | `S2.5′` (warnings pass) | — (does not block) |
| `S2.5′ PREREGISTER` | human + LLM | `S5` | — (does not block. It can proceed with v0 alone) |
| `S5 PLAN` | LLM proposes → code validates | `S6` | `S5` (re-propose, deducted from the budget) |
| `S6 PREFLIGHT` | code | `S7` | `S5` |
| `S7 EXECUTE` | Runner | `S7.5` | `S7.5` (a crash is a diagnosis target too) |
| `🚦GATE_PROD` | human | `S7:production` | `S5` \| `DONE` (finish at pilot) |
| `S7.5 EYEBALL` | code → human | `S8` | `S8` (a visual anomaly is a diagnosis input too) |
| `S8 DIAGNOSE` | code + LLM triage | `S10` | `S9`(REPAIR) \| `S5`(REDESIGN) |
| `S9 REPAIR` | rule → LLM → human | `S6` (re-run after re-checking) | `ESCALATED` (budget exhausted) |
| `S10 ANALYZE` | code | `S11` | `S9` (`N_eff` insufficient → extend the run) |
| `S11 VISUALIZE` | code | `S12` | `S11` (renderer fallback) |
| `S12 REPORT` | LLM + code | `DONE` | — |

**Why S2.5 and S7.5 were added (`D28`, 2026-07-27):** pre-registration is the only device that prevents post-hoc
interpretation, so it has to be a stage rather than an appendage; and the eyeball check has to be *before* the analysis —
because failures like crystallization, clustering and overlap are real, and they take **1 second by eye and are hard by
number**. The details are in [`master_plan.md`](../../docs/history/2026-08_simulation_auto_master_plan.ko.md) §5.
`S2.5′` is not a new state but **the second pass through S2.5** (qualitative v0 → quantitative v1).

**Terminal states:** `DONE` · `ESCALATED` (needs human intervention) · `BLOCKED_INPUT` · `ABORTED` (stopped by the user)

**The idempotency rule:** each stage produces the same output given the same input artifacts (LLM stages at temperature 0 with a fixed seed).
So it can be resumed from any stage.

---

## 3. `run_state.yaml`

The single source of truth for one run. **Updated atomically after every stage transition** (a temporary file → `os.replace`).

```yaml
schema_version: 1
run_id: "2026-07-27T11-17-03Z__silica-depletion-gel"
state: S8_DIAGNOSE
tier: pilot                       # smoke | pilot | production
created_at: "2026-07-27T11:17:03Z"
updated_at: "2026-07-27T11:42:18Z"

stage_history:
  - {stage: S1_INTAKE,  entered: "...", exited: "...", result: PASS}
  - {stage: S2_ELICIT,  entered: "...", exited: "...", result: PASS}
  - {stage: GATE_SPEC,  entered: "...", exited: "...", result: APPROVED, by: human}
  - {stage: S2_5_PREREGISTER, entered: "...", exited: "...", result: PASS, version: v0}
  - {stage: S3_NONDIM,  entered: "...", exited: "...", result: PASS}

gates:
  spec:       {status: approved, at: "2026-07-27T11:25:00Z", by: human}
  production: {status: pending}

artifacts:                        # relative paths. Absent means not yet generated
  spec_draft:   artifacts/spec_draft.yaml
  spec:         artifacts/spec.yaml
  hypothesis_v0: artifacts/hypothesis-v0.yaml   # S2.5  qualitative. Basis = intuition, experience
  hypothesis_v1: artifacts/hypothesis-v1.yaml   # S2.5′ quantitative. Basis = dimensionless numbers, literature
  reduced_spec: artifacts/reduced_spec.yaml
  unit_map:     artifacts/unit_map.yaml
  dimensionless: artifacts/dimensionless.yaml
  grounding:    artifacts/grounding.md
  run_plan:     artifacts/run_plan.yaml
  cost_estimate: artifacts/cost_estimate.yaml
  preflight:    artifacts/preflight_report.md
  sim_script:   artifacts/simulate.py
  eyeball:      figures/eyeball/                # S7.5 3 snapshots (low resolution)
  diagnosis:    artifacts/diagnosis.yaml

budget:                           # D18
  max_total_walltime_s: 21600
  max_repair_iterations: 8
  max_disk_gb: 20
  max_llm_calls: 100
  spent_walltime_s: 412
  repair_iterations_used: 1
  disk_used_gb: 0.8
  llm_calls_used: 14

provenance:                       # the crux of reproducibility. See D6
  bdkit_version: "0.1.0"
  git_commit: "TBD"               # D6 undecided, so "no-vcs"
  hoomd_version: "7.1.0"
  hoomd_gpu_enabled: false
  freud_version: "TBD"
  python: "3.12.13"
  platform: "macOS-15-arm64"
  master_seed: 12345
```

---

## 4. `decision_journal.jsonl`

**The only means of debugging the agent.** One line = one decision. Append-only.

The `actor` field distinguishes `rule` (a deterministic rule) / `llm` (a model judgment) / `human` (a person).
Only with that can one count, after the fact, **"how many times did an LLM judgment enter this result"**.

```jsonl
{"ts":"...","stage":"S1_INTAKE","actor":"llm","action":"extract_field","field":"particle_radius_m","value":5.0e-7,"confidence":0.6,"provenance":"estimated from the scale bar '1 µm' in the photo","note":"confidence<0.8 → also entered in unknowns[]"}
{"ts":"...","stage":"GATE_SPEC","actor":"human","action":"approve","note":"confirmed the radius is 500nm"}
{"ts":"...","stage":"S3_NONDIM","actor":"rule","rule_id":"U03_overdamped_validity","observation":{"tau_B_over_tau_D":2.1e-7},"verdict":"PASS","threshold":1e-3}
{"ts":"...","stage":"S8_DIAGNOSE","actor":"rule","rule_id":"G02_max_displacement","observation":{"max_disp_per_step_sigma":0.34,"threshold":0.10},"verdict":"FAIL"}
{"ts":"...","stage":"S9_REPAIR","actor":"rule","rule_id":"R02_reduce_dt","action":{"param":"run_plan.dt","from":1.0e-4,"to":5.0e-5},"iteration":1,"rationale":"the rule table hit"}
{"ts":"...","stage":"S9_REPAIR","actor":"llm","action":"triage","observation":"no rule hit: the nlist rebuild frequency is anomalous","proposal":{"param":"run_plan.nlist_buffer","from":0.4,"to":0.8},"iteration":3,"rationale":"..."}
{"ts":"...","stage":"S9_REPAIR","actor":"human","action":"escalate_ack","note":"budget exhausted; decided to swap the potential for a harmonic one instead of dt"}
```

**Required fields:** `ts` `stage` `actor` `action`
**Recommended fields:** `rule_id` (actor=rule) · `observation` (the measured value + the threshold) · `rationale` · `iteration`

---

## 5. The run directory layout

```
outputs/<run_id>/
├── run_state.yaml                  ← the single source of truth
├── decision_journal.jsonl          ← append-only
├── artifacts/
│   ├── spec_draft.yaml  spec.yaml  hypothesis-v0.yaml  hypothesis-v1.yaml
│   ├── reduced_spec.yaml  unit_map.yaml  dimensionless.yaml
│   ├── grounding.md
│   ├── run_plan.yaml  cost_estimate.yaml
│   ├── preflight_report.md
│   ├── simulate.py                 ← the generated HOOMD script (a human must be able to read it and run it by hand)
│   ├── diagnosis.yaml
│   └── observables.parquet  observables_summary.csv
├── raw/
│   ├── smoke/       trajectory.gsd  log.h5
│   ├── pilot/       trajectory.gsd  log.h5  checkpoint.gsd
│   └── production/  trajectory.gsd  log.h5  checkpoint.gsd
├── figures/                        ← PNG (embedded in the report as base64)
│   └── eyeball/                    ← S7.5 3 low-resolution snapshots. For the eyeball check before analysis
└── report.html                     ← the final product, a self-contained single file
```

**`artifacts/simulate.py` must always be directly runnable by a human.** It is the escape hatch for when you do not
trust the agent, and the fastest way to check what the agent actually did.

The `run_id` convention: `<UTC ISO8601, colons→hyphens>__<slug>` — lexical order = chronological order.

---

## 6. The resume convention

```
bd-agent resume outputs/<run_id> [--from S5_PLAN]
```

1. Read `run_state.yaml` and restore `state`
2. Without `--from`, start from the last incomplete stage
3. Artifacts that already exist are not recomputed (override with `--force`)
4. `S7 EXECUTE` continues from `checkpoint.gsd` if one exists
5. The fact of the resume is itself journaled

---

## 7. Failure-handling principles

| Failure type | Handling |
|---|---|
| **a verification failure** (an expected one) | the rule table (`06`) → automatic correction within budget |
| **no rule hits** | 1 LLM triage → apply the proposal after passing it through the deterministic validator |
| **budget exhausted** | `ESCALATED`. Present to the human a summary of **the symptom · everything that was tried · the next candidates** |
| **a crash** (an unexpected one) | report the stack trace + `run_state` + the last 20 journal lines together. The state is preserved |
| **contradictory input** | `BLOCKED_INPUT`. Return pointing at the contradiction (e.g. "φ=0.7 with hard spheres — physically impossible") |

**What will never be done:** quietly swallowing a failure and moving to the next stage.
A stage that passed without verification is marked in the report with an `UNVERIFIED` badge.

---

## 8. How this architecture is verified

| Target | Method |
|---|---|
| `bdkit/` pure logic | `pytest tests/` — works without HOOMD (inheriting the 33-test convention of the preceding project) |
| physical consistency | `pytest tests/test_benchmarks.py` — run a short simulation → compare against literature values (`09`) |
| the state machine | tests that transition each stage with mock artifacts. 1:1 with the transition table (§2) |
| LLM extraction accuracy | ~10 natural-language prompts → compare against expected `SystemSpec` golden files |
| journal completeness | aggregate `actor` after a run finishes. Is the number of LLM judgments in the expected range |
