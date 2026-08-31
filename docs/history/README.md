# Predecessor design documents — translated, structure verified

These are the design documents of the three private repositories this project
was merged from. They arrived **untranslated and otherwise unedited except for one
mechanical change: relative link paths were re-based onto this repository's
layout** (71 of them), because these files used to sit at their own repo roots
and now sit two levels down. They are the provenance record
for every decision in [`00-merge-decisions.md`](../00-merge-decisions.md).

**They were translated to English on 2026-08-31** at the user's instruction, with
the correction markers preserved explicitly. What that means concretely is below.

| File | Generation | Lines | What it is authoritative for |
|---|---|---:|---|
| [`2026-07-28_bd_agent_master_plan.md`](2026-07-28_bd_agent_master_plan.md) | `BD_agent` | 1,255 | **Verification philosophy** — the 7-rung process ladder, the 4 evidence layers, the uncertainty ledger, autonomy boundaries, the public-boundary table |
| [`2026-07_bd_agent_00_decision_log.md`](2026-07_bd_agent_00_decision_log.md) | `BD_agent` | 470 | numbered decisions `D1`…`D32`, each with what would reverse it and at what cost |
| [`2026-07_bd_agent_01_agent_architecture.md`](2026-07_bd_agent_01_agent_architecture.md) | `BD_agent` | 287 | 10 general principles for building agents; the LLM-layer comparison table |
| [`2026-07-30_simulation_bot_master_plan.md`](2026-07-30_simulation_bot_master_plan.md) | `Simulation_bot` | 1,632 | **The S1→S8 pipeline** stage by stage, with gates and failure modes · non-dimensionalization conventions · the `knowledge/` schema · sensitivity analysis (S7b) · model tiering · the progress dashboard |
| [`2026-08_simulation_auto_master_plan.md`](2026-08_simulation_auto_master_plan.md) | `simulation_auto` | 2,305 | **The L0–L10 layer design** · the non-dimensionalization engine in detail · raw-data storage strategy · the HOOMD API map · the closed learning loop |
| [`2026-08_simulation_auto_CLAUDE.md`](2026-08_simulation_auto_CLAUDE.md) | `simulation_auto` | 724 | **The case-by-case results log** as it stood on 2026-08-06 — the primary source behind [`04-cases.md`](../04-cases.md) |

---

**Four links still dangle**, pointing at files that were genuinely dropped in the
merge: `docs/11_simbot.md` (3×) and `other-rule.md` (a template placeholder that
never existed). They are left broken rather than redirected, because a link to a
document that no longer exists is itself part of the record.

---

## How they were translated, and why that was the hard part

These are not reference documentation; they are **a corrections trail**. Their
value is concentrated in passages that record being wrong: *"correction
(2026-07-28, re-verified the three items above): the second is true only for the
pipeline…"*, *"⚠️ discarded — the queue reading as of 2026-07-28, kept as a
record"*, *"this table getting longer is good news."* Those passages are precise
about who believed what, when, and what measurement changed it. Paraphrase them
loosely and you smooth exactly the seams that make them useful, and a smoothed
history is worse than no history. **That was the risk, and it is what the
verification was built against.**

Every file was checked by
[`verify/verify_markdown_translation_safety.py`](../../verify/verify_markdown_translation_safety.py),
which compares OLD against NEW and reports every difference that is *not* prose:

- the sequence of heading levels
- every line of every fenced code block, byte for byte (a line that carried
  Hangul is a code line with a translated comment, so only its code prefix is
  required to survive)
- the link and image TARGETS in order; in-document anchors are retargeted with
  their headings and required to **resolve** to a heading in the same file
- each table's row count and per-row cell count
- **the multiset of numeric literals** — this is the one that matters. Every
  number in these documents is a measurement, and a translation that perturbs one
  is not a translation
- the census of correction markers — ★ ⭐️ ⚠️ ⛔ ✅ ❌ ⟹ → and the circled digits
  — plus the count of `**` emphasis delimiters, so an emphasis cannot be lost

It found **206 defects** across the 7 files (these 6 plus
[`../hoomd_capabilities.md`](../hoomd_capabilities.md)) and every one was fixed.
The pattern was consistent: arrows prosified into "from X to Y", `·`-separated
section lists rewritten as English comma lists (which turns the token `5` into
`5,`), and ASCII directory trees where translating a trailing comment moved the
code that precedes it. None of those change the meaning, which is exactly why a
census catches them and reading does not.

Two numerals are accepted by name rather than matched, because a Korean myriad
form (`10만`, `330만`) has no digit-preserving English rendering: the acceptance
is recorded in the verify command, not absorbed silently.

**They may still be stale.** Where a number in these files disagrees with
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
