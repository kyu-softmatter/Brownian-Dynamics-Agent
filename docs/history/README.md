# Predecessor design documents — verbatim, in Korean

These are the design documents of the three private repositories this project
was merged from, carried **untranslated and otherwise unedited except for one
mechanical change: relative link paths were re-based onto this repository's
layout** (71 of them), because these files used to sit at their own repo roots
and now sit two levels down. No prose was altered. They are the provenance record
for every decision in [`00-merge-decisions.md`](../00-merge-decisions.md).

| File | Generation | Lines | What it is authoritative for |
|---|---|---:|---|
| [`2026-07-28_bd_agent_master_plan.ko.md`](2026-07-28_bd_agent_master_plan.ko.md) | `BD_agent` | 1,255 | **Verification philosophy** — the 7-rung process ladder, the 4 evidence layers, the uncertainty ledger, autonomy boundaries, the public-boundary table |
| [`2026-07_bd_agent_00_decision_log.ko.md`](2026-07_bd_agent_00_decision_log.ko.md) | `BD_agent` | 470 | numbered decisions `D1`…`D32`, each with what would reverse it and at what cost |
| [`2026-07_bd_agent_01_agent_architecture.ko.md`](2026-07_bd_agent_01_agent_architecture.ko.md) | `BD_agent` | 287 | 10 general principles for building agents; the LLM-layer comparison table |
| [`2026-07-30_simulation_bot_master_plan.ko.md`](2026-07-30_simulation_bot_master_plan.ko.md) | `Simulation_bot` | 1,632 | **The S1→S8 pipeline** stage by stage, with gates and failure modes · non-dimensionalization conventions · the `knowledge/` schema · sensitivity analysis (S7b) · model tiering · the progress dashboard |
| [`2026-08_simulation_auto_master_plan.ko.md`](2026-08_simulation_auto_master_plan.ko.md) | `simulation_auto` | 2,305 | **The L0–L10 layer design** · the non-dimensionalization engine in detail · raw-data storage strategy · the HOOMD API map · the closed learning loop |
| [`2026-08_simulation_auto_CLAUDE.ko.md`](2026-08_simulation_auto_CLAUDE.ko.md) | `simulation_auto` | 724 | **The case-by-case results log** as it stood on 2026-08-06 — the primary source behind [`04-cases.md`](../04-cases.md) |

---

**Four links still dangle**, pointing at files that were genuinely dropped in the
merge: `docs/11_simbot.md` (3×) and `other-rule.md` (a template placeholder that
never existed). They are left broken rather than redirected, because a link to a
document that no longer exists is itself part of the record.

---

## Why these are not translated

The stated language policy for this repository is English for the authored
surface — README, `NOTICE`, `CLAUDE.md`, `docs/0*`. These files are exempt on
purpose.

They are not reference documentation; they are **a corrections trail**. Their
value is concentrated in passages that record being wrong: *"correction
(2026-07-28, re-verified the three items above): the second is true only for the
pipeline…"*, *"⚠️ discarded — the queue reading as of 2026-07-28, kept as a
record"*, *"this table getting longer is good news."* Those passages are precise
about who believed what, when, and what measurement changed it. Paraphrasing
them into English would smooth exactly the seams that make them useful, and a
smoothed history is worse than no history.

So they stay as written, and the English documents cite them by section rather
than restating them.

**They may also be stale.** Where a number in these files disagrees with
[`06-roadmap.md`](../06-roadmap.md), the roadmap is the newer measurement — for
instance `Simulation_bot`'s dashboard reports 562 tests, and this merged tree
measures 572 with 6 of them repaired during the merge.

## Reading order, if you read them

1. `bd_agent_master_plan` §6 — verification with no grader. The most transferable
   part of the whole project, and independent of Brownian dynamics.
2. `simulation_bot_master_plan` §2 — the eight stages, each defined as
   *input → processing → artifact → gate → failure mode*.
3. `simulation_auto_CLAUDE` — what actually happened, case by case, including the
   several occasions on which a stated conclusion was later reversed.
