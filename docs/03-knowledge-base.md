# 03 · Knowledge base

The premise: **`bdbot/` and `simbot/` can be rewritten; `knowledge/` cannot.**
*Why this material got `η = 1.2 mPa·s`* and *why `dt = 5e-5` exploded* disappear
unless they are accumulated, and re-deriving them costs the same weeks it cost
the first time.

Three things go wrong when knowledge is not accumulated deliberately:

- **Parameter choice is tacit.** What `dt`, where to cut `r_cut`, how long to
  equilibrate — mostly learned from a senior student or copied from a paper.
  It is rare to be able to say *why* that value, and there is no procedure for
  checking whether it still holds in a new system.
- **Verification is unsystematic.** The simulation always produces a number.
- **Failure is not recorded.** "I tried these parameters and it didn't work"
  usually lives in one person's head. Six months later the same person, or the
  next student, hits the same wall. **Failures vanish and only successes reach
  papers.**

| Problem | Response |
|---|---|
| tacit parameters | every parameter carries a **dimensionless-group rationale** and a **literature source** |
| unsystematic verification | literature benchmarks run **as regression tests**; every result wears an **evidence grade** |
| vanishing failure | a failure becomes a `dead-end` page, and the next attempt **queries it first** |

---

## 1 · Two schemas, not yet one

⚠️ This is the largest piece of debt in the repository. The merge carried two
knowledge stores with different shapes, different writers and different readers.

| | [`knowledge/wiki/`](../knowledge/wiki/) + [`source/`](../knowledge/source/) | [`knowledge/entries/`](../knowledge/entries/) |
|---|---|---|
| Form | Markdown + YAML frontmatter | flat JSON, one file per claim |
| Written by | a human or an agent, deliberately | `tools/kb.py add` |
| Read by | skill `bd-knowledge`, and by reading | `tools/kb.py query` |
| Size | 46 wiki pages · 42 paper + 2 book distillations | 126 entries |
| Keyed on | kind (systems / findings / …) | `origin` × `kind` |
| Strength | contracts, cross-links, human-legible rationale | cheap to append, so it actually gets appended |

A lesson filed in one is not found by a tool reading the other. Unifying them is
on the [roadmap](06-roadmap.md); the honest interim rule is **query both**.

---

## 2 · The wiki taxonomy — six kinds, because they are verified differently

Layers differ in how fast they change and how they are checked: **code is
verified by tests, knowledge by citation, rules by an accident report.**

| Kind | Count | What | Contract |
|---|---:|---|---|
| **`systems/`** | 11 | ★ **(system × target dynamics)** cards | **the card owns the non-dimensionalization and the gates.** Not the pipeline, not the case script |
| **`findings/`** | 23 | Q→A, and **dead-ends** | summary · basis · scope/limits · references |
| **`benchmarks/`** | 5 | systems whose answer is known, run as regression tests | input · observable · expected · tolerance · `evidence_layer` · source |
| **`concepts/`** | 3 | a physical concept in this project's conventions | — |
| **`techniques/`** | 2 | method notes, e.g. the environment log | — |
| **`questions/`** | 2 | open questions, with what would close them | — |
| **`source/papers/`** | 42 | 1:1 paper distillations | citation · doi · `read_depth` · `provides` · `used_by` · `lab_authored` |
| **`source/books/`** | 2 | book distillations | same, with `[BOOK]`/`[DERIV]`/`[OURS]` check labels |

**A `systems/` card owning the scales is the load-bearing decision.** If the
pipeline owned them, every new system would mean editing the pipeline; if the
case script owned them, two cases on the same physics could silently disagree.
The card is the single place a scale choice can be argued with.

A `findings/` page has five sections and the fifth is the point:

```markdown
## Summary        one paragraph — what is now known
## Basis          the measurement or derivation, with numbers
## Scope / limits where this stops being true
## References    links to runs, specs, papers
## Prevention     what changes so this does not recur
```

A **dead-end** page is a finding whose answer is "this route does not work."
Recording it is worth as much as recording a success, and it must state a
**cause, not a symptom** — see
[02 §5](02-verification.md#a-reason-must-be-a-cause-not-a-symptom).

---

## 3 · The `entries/` store — what actually gets written

126 entries, by origin. The distribution is itself a finding:

| `origin` | Count | What it means |
|---|---:|---|
| `method` | 44 | how to do something, learned by doing it |
| **`tooling`** | **44** | **a tool bit us** |
| `handbook` | 25 | a number or a formula from a reference work |
| `intake` | 10 | how to read a sketch |
| `paper` | 3 | a claim from a paper |

⚠️ **44 tooling entries against 3 paper entries is not a good sign about the
tools.** More time went into "our own machinery silently misled us" than into
"the literature said something." The tooling entries are, individually, why the
gates in [02](02-verification.md) exist; collectively they are an argument for
the two-engine seam being paid down rather than lived with.

```bash
$PY tools/kb.py list
$PY tools/kb.py query --tags dlvo --origin tooling --kind pitfall
$PY tools/kb.py lessons
$PY tools/kb.py add --origin tooling --kind pitfall \
    --source "bdbot/metrics.py#build" --claim "…"
```

Plus **227 `record.json`** post-mortems under `runs/`, written by
`tools/postmortem.py`. A run is not finished when it exits; it is finished when
its post-mortem exists.

```bash
$PY tools/postmortem.py runs/<run_id> --lesson "lesson::pitfall::coord=value"
```

---

## 4 · Provenance and tiers — a number without a source is not a number

Every parameter carries where it came from. The vocabulary:

`from_drawing` · `from_paper` · `from_knowledge` · `assumed` · `derived`

and a **tier**, which is what makes an inherited value auditable. The failure
this catches is subtle and it happened: **`T = 300 K` is recorded as tier 1
(measured) but it is a choice** — it was inherited from `trap-2d-5um`, whose
sketch contains no temperature. Water viscosity is 2.06 %/K sensitive, so that
mislabel propagates a −4 % to −14 % error into every `τ_B` downstream. The tier
field is the only reason that was findable at all.

Rules that fall out of this, and are enforced by code:

- **`BLOCKED` is not a failure.** It means the missing input has been narrowed to
  one thing. Forcing it to `READY` by inventing a number is the failure.
- **`null` for anything absent from the sketch.** Skill `bd-intake`,
  eight anti-invention rules, derived from reading five real sketches.
- **A derived value is recomputed and compared**, never trusted as written.

⚠️ **`run_id` must hash the physics and nothing else** —
[`bdbot/runid.py`](../bdbot/runid.py) `DOC_KEYS`. Both directions have bitten:

- Adding a `derived_from` field once **changed a run's id** although no physics
  changed, invalidating a completed run.
- Worse in reverse: the `soft-r3` spec contained **no physical system at all**,
  so changing `d` 5 µm→0.5 µm and `η` by 62× (a 16.1× change in τ_B) left the
  `run_id` **identical** — and the old result was reported as the new system's
  result. `verify/verify_l3_spec_gaps.py` guards both directions.

---

## 5 · Capturing judgment out of conversation

**This is the real purpose of the project**, and it is the part that has no
tooling yet.

Most of what is worth keeping is said in passing while doing something else —
"that's the wrong timescale," "we already tried that and it collapsed," "that
threshold is convention, not literature." It is not a run artefact, so
`record.json` has no slot for it, and it is gone when the session ends.

The interim mechanism is manual and it works:

```bash
$PY tools/kb.py add --origin intake|paper|tooling|method|handbook \
    --kind pitfall --source "file#anchor" --claim "…"
```

`origin: handbook` has a convention worth copying: the distillation lives in
`knowledge/source/books/` and the entry's `source` reads
`distillation#section ← [short-name] p.page`, so a claim can be walked back to
the page it came from without the book being in the repository.

**Rules are born from accidents, not written in advance.** If a rule has no
dated, costed accident attached, it becomes a ritual nobody understands, and then
there is no basis for deciding whether a changed situation retires it.
[`.claude/rules/`](../.claude/rules/) holds four, each 25–80 lines, one topic per
file.
