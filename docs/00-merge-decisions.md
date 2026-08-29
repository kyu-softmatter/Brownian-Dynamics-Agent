# 00 · Merge decisions — what came from where

This repository is the merge of three private predecessors, assembled
2026-08-28. They are not three alternatives: they are three **generations**, and
each generation stopped being worked on when the next began.

| Generation | Last commit | Scale at handover | What it was best at |
|---|---|---|---|
| **`BD_agent`** | 2026-07-28 | 14-stage pipeline, 3 layers, 81 py files | **Verification philosophy.** The only one that asked "how do you verify in a domain with no grader" and answered it structurally |
| **`Simulation_bot`** | 2026-07-30 | S1→S8 pipeline, 19 modules, **562 tests**, 3,856 production runs | **Agent orchestration.** The only one with a subagent roster, model tiering, prediction sealing and a PASS/FAIL/**INCONCLUSIVE** verdict |
| **`simulation_auto`** | 2026-08-17 | `bdbot` 21 modules, 8 cases, 278 specs, 254 runs, 353 KB entries | **Physics that actually ran.** The only one that took real systems to defended quantitative conclusions |

The merge followed three stated criteria, in this order:

1. **Adopt whichever was designed with the broader scope.**
2. **Carry whichever was developed more concretely.**
3. **Prefer the more recent.**

The criteria mostly do not conflict, because the generations specialized in
different directions. Where they did conflict, criterion order decided, and this
document records it.

---

## 1 · Execution engine — `bdbot/` wins outright

`simulation_auto`'s [`bdbot/`](../bdbot/) (21 modules) is both the most recent
and the most developed engine, and it is the only one that reached a **content-
addressed run** with a numerical-health verdict. All three criteria point the
same way. It is the live engine.

What that buys, concretely: an `L3` contract
([`bdbot/nondim.py`](../bdbot/nondim.py)) that is the *only* thing passed
between the physical system and the numerics, so the health layer never imports
case code; a builder registry (`@RUN.builder`) that every one of the 8 cases
uses; and a `run_id` that hashes the physics fields only, so editing a comment
does not invalidate a run.

`BD_agent/simbot/` and `BD_agent/bdkit/` were **dropped**. Everything they did
that survived had already been rewritten better twice.

## 2 · Pipeline stages — `simbot/` fills a real hole in `bdbot/`

Here the criteria conflict, and criterion 1 wins.

`bdbot` is deep on **L0→L7** (intake → non-dimensionalize → gate → run →
store). It has no notion of **predicting before running**, no verdict
vocabulary, and no report generator. `Simulation_bot`'s
[`simbot/`](../simbot/) has exactly those, and they are tested (562 tests):

| Module | What `bdbot` lacks | Why it matters |
|---|---|---|
| `estimators.py` | analytic-solution / scaling engine for S2 | you cannot pre-register a prediction you cannot compute |
| `io.py` | **prediction sealing** (`SEALED.sha256`) | blocks post-hoc rationalization *structurally*, not by discipline |
| `validate.py` | PASS / FAIL / **INCONCLUSIVE** + design power | "not measurable by this design" is a distinct verdict from "wrong" |
| `report.py` | `REPORT.md` generation | tested for whether the bad news can go missing |
| `session.py`, `policy.py` | append-only session state, resource policy | a sweep that costs 25 days needs a budget gate before it starts |
| `analysis/structure.py` | RDF, `ψ₆`, Voronoi, `S(k)`, finite-size scan, bootstrap | 37 tests, written *after* the campaign existed — which is why it caught 3 bugs |
| `viz.py` | figure generation that **refuses to emit an uncaptioned figure** | 34 tests |

**So both are carried, side by side, unmerged.** `bdbot/` is the engine;
`simbot/` is the S2/S7/S8 half. This is honest rather than tidy — see
[Known seams](#5--known-seams) below.

## 3 · Agent layer — merged, because the two sets do different jobs

`.claude/` combines all three generations:

| From | What | Why it won |
|---|---|---|
| `Simulation_bot` | skills `bd-pipeline` (+ 5 stage references), `bd-diagnose`, `bd-knowledge`; **9 subagents** with model tiering; `settings.json` that refuses to edit a sealed document | criterion 1 — the only orchestration design that exists |
| `simulation_auto` | skills `bd-hoomd` (20 measured traps), `bd-physics` (758 lines), `bd-intake` | criteria 2+3 — far more concrete domain knowledge, and newer |
| `BD_agent` | `.claude/rules/` — 4 rules | criterion 1 — the only generation with a rules layer, and its "rules are born from accidents, not written in advance" convention |

The two skill classes are **not** interchangeable, and conflating them was a
real merge defect: the orchestration skills are mutually exclusive routers, so
they must cross-reference each other; the domain skills are reference material
the pipeline *reads at a stage*. Immediately after merging, nothing pointed at
`bd-hoomd` from the pipeline — so S5 would write HOOMD code without reading the
trap list. `tests/test_agent_layer.py::test_pipeline_points_at_every_domain_skill`
now guards that.

## 4 · Knowledge base — taxonomy from one, contents from all

Criterion 1 gives the **shape**, criterion 2 gives the **contents**.

- **Shape**: the `wiki/{systems, findings, benchmarks, concepts, questions,
  techniques}` + `source/{papers, books}` taxonomy, from
  `BD_agent` → `Simulation_bot`. Six kinds of knowledge, each with its own
  schema. `Simulation_bot` had already absorbed `BD_agent`'s wiki and was a
  strict superset (11 vs 6 systems, 23 vs 7 findings), so there was nothing to
  reconcile.
- **Contents**: `Simulation_bot`'s wiki (46 pages) + 42 paper distillations,
  plus `simulation_auto`'s **126 flat `kb/entries/*.json`** and **227 run
  `record.json`** files, plus the 2 book distillations.

⚠️ **These two knowledge schemas are not unified.** `knowledge/wiki/` is
human-written Markdown with frontmatter contracts;
[`knowledge/entries/`](../knowledge/entries/) is tool-written JSON keyed by
`origin` (`intake` · `paper` · `tooling` · `method` · `handbook`). They are
queried by different tools (`bd-knowledge` skill vs `tools/kb.py`). Unifying
them is the largest outstanding piece of debt in this repository — see
[06 Roadmap](06-roadmap.md).

## 5 · Known seams

Recording these because a merge that claims to be seamless is lying.

| Seam | State | Consequence |
|---|---|---|
| **Two engines** | `bdbot/` (live, 8 cases) and `simbot/` (S2/S7/S8, tested) do not call each other | a case run through `bdbot` gets a health verdict but **no sealed prediction**; a case run through `cli.py` gets sealing but cannot use `bdbot`'s cases |
| **Two knowledge schemas** | `wiki/` Markdown vs `entries/` JSON | a lesson can be filed in either and found by only one |
| **`dt` gate logic duplicated** | `simbot.nondim.choose_dt` and `campaigns/chain_bend.py:113` reimplement the same thresholds | already caused a silent miss: `choose_dt`'s displacement gate keys off `bool(spec.pair)`, and `SystemSpec` has no bond/angle field, so **a bond-only system silently turns the gate off** — and the measurement that needed it (`max|F*| = 1037.7`, `dt` cut 100×) came from exactly such a system |
| **`result.txt` is written by the case script, not by `bdbot.run`** | `chain-bend-2d-oscill` (7 runs with metrics) and `trap-drag-2d-hex300` (80 runs with metrics) report **0 runs** in `bdbot.cli status` | the counting convention makes real runs invisible; it once caused a "clean up incomplete runs" pass to **delete 6 completed runs** |
| **Language** | README, `NOTICE`, `CLAUDE.md` and `docs/0*` are English; skills, knowledge base and `docs/history/` are Korean | a reader of the English surface cannot read the evidence layer without translation |

Fixed during the merge rather than carried: `cases/network_3d.py` existed and
had produced runs but was missing from `bdbot/cli.py`'s `CASE_SCRIPTS`, so
`bdbot.cli run network` refused with "no end-to-end script" and the status table
showed `—`. Now registered.

## 6 · What was dropped entirely

| Dropped | Why |
|---|---|
| `BD_agent/simbot/`, `BD_agent/bdkit/`, `BD_agent/agent/` | superseded twice; `agent/llm.py` called the Anthropic API directly, which the Claude-Code-native design replaced |
| `BD_agent/notes/`, `BD_agent/readings/`, `BD_agent/configs/` | scratch work from a 2-day generation, all of it re-derived later |
| `Simulation_bot/sessions/`, `BD_agent/outputs/` | machine-local state, already gitignored in their own repos |
| 3,866 per-replica `manifest.json` files | `Simulation_bot`'s own `.gitignore` documents the measurement: identical provenance fields, only `run_hash`/`seed`/`wall_s` differ, and those recompute from `spec + seed` |
| All PDFs, all `.gsd`/`.npz` | see [NOTICE](../NOTICE.md) |

The three predecessor repositories still exist locally and were not modified by
this merge. Their design documents are carried verbatim in
[`docs/history/`](history/) as the provenance record for every decision above.
