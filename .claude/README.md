# `.claude/` — the L1 agent layer

The core computes; this layer judges. **There is no number-producing code here.**

## Contents

```
.claude/
├── settings.json                  permissions — allow the interpreter, forbid editing a sealed document
├── rules/                         4 rules, each born from a dated accident
│   ├── axioms.md                  the four that need human approval to change
│   ├── deterministic-core.md      the core does not call an LLM
│   ├── overdamped-stability.md    dt is set by force, not by a timescale ratio
│   └── verify-against-literature.md   literature is cited, not remembered
├── skills/
│   ├── bd-pipeline/               [main] the S1→S8 orchestrator
│   │   ├── SKILL.md               stages, gates, and the checklist of prohibitions
│   │   └── references/
│   │       ├── s1_intake_drawing.md   ★ reading a hand drawing — the only content unique to the skill layer
│   │       ├── s2_prediction.md       prediction, sealing, power
│   │       ├── s3_s5_execute.md       specify, non-dimensionalize, run (mostly CLI calls)
│   │       ├── s6_s7_validate.md      figures and verdicts
│   │       └── s8_knowledge.md        conclusion and the knowledge commit
│   ├── bd-diagnose/SKILL.md       diagnosing a broken run (elimination order)
│   ├── bd-knowledge/SKILL.md      searching, adding and tidying knowledge
│   ├── bd-intake/SKILL.md         reading sketches — 8 anti-invention rules
│   ├── bd-physics/SKILL.md        scale tables, non-dimensionalization, inversion
│   └── bd-hoomd/SKILL.md          20 measured HOOMD traps + verified snippets
└── agents/                        9, tiered by the `model:` frontmatter
```

Two skill classes live here and they are **not** interchangeable.
`bd-pipeline`, `bd-diagnose` and `bd-knowledge` are mutually exclusive routers,
so they must cross-reference each other or the wrong procedure runs.
`bd-intake`, `bd-physics` and `bd-hoomd` are reference material the pipeline
**reads at a stage** — they do not self-trigger, so the pipeline has to point at
them. `tests/test_agent_layer.py` guards both properties; the second check was
empty immediately after the 2026-08-28 merge, which would have meant S5 writing
HOOMD code without reading the trap list.

## Why the references are split into five (decided 2026-07-28)

The original design had eight, one per stage. But **once the deterministic core
was finished, S3 · S4 · S5 collapsed into a single CLI invocation**, and giving
each its own document produces three thin files that say only "call this
function."

So instead, **the documents live where the content is:**

| Document | Why it stands alone |
|---|---|
| `s1_intake_drawing.md` | **The only stage that cannot be expressed as code.** The most expensive place to be wrong |
| `s2_prediction.md` | Sealing, tolerance and power discipline. Sloppiness here disables verification entirely |
| `s3_s5_execute.md` | All three are a CLI call plus reading a gate. Merged, the flow is visible |
| `s6_s7_validate.md` | Figures and verdicts serve the same judgment — *what looks wrong* |
| `s8_knowledge.md` | Writing the conclusion plus the knowledge contract |

## Model tiering

> Principle: **extraction is cheap, interpretation is expensive. And calculation
> belongs to code, not to a model.**

| Agent | model | Responsibility |
|---|---|---|
| `bd-intake-extract` | haiku | S1 text, number and file extraction |
| **`bd-intake-interpret`** | **opus** | S1 physical interpretation — wrong here means wrong everywhere |
| **`bd-predict`** | **opus** | S2, the claim that gets sealed |
| `bd-spec` | sonnet | S3 rule application |
| **`bd-validate`** | **opus** | S7 verdict and causal inference |
| **`bd-conclude`** | **opus** | S8 final claim |
| `bd-lit-distill` | sonnet | literature distillation (equation transforms go to Opus for review) |
| `bd-lit-scan` | haiku | bibliography and INDEX bulk work |
| **`bd-diagnose`** | **opus** | failure diagnosis |

**The purpose is not cost saving but allocating speed against quality.**

<a name="authority-boundary"></a>
### The authority boundary — mechanically checkable

> **Only Opus may write a field whose `provenance` is `inference` or `assumed`.**

`observation` · `derived` · `rule` · `from_knowledge` may be filled by a cheaper
model. The boundary is stated in the instructions of `bd-spec` (sonnet),
`bd-intake-extract` (haiku), `bd-lit-distill` (sonnet) and `bd-lit-scan` (haiku),
`tests/test_agent_layer.py` asserts that each of them states it, and
`simbot.spec.Quantity.problems()` checks `written_by` at runtime.

The reason this is a hard boundary and not a preference: an `assumed` value is a
physical assumption entering the system without a human ever seeing it. Every
number downstream inherits it, and nothing in the pipeline can distinguish it
later from a measured one — which is exactly the failure the `tier` field exists
to make findable (`T = 300 K` is the live example; see
[docs/03](../docs/03-knowledge-base.md#4--provenance-and-tiers--a-number-without-a-source-is-not-a-number)).

## The two decisions in `settings.json`

**① Allow the interpreter's absolute path.** `conda activate` is unreliable in a
non-interactive shell, so it is on the **deny** list — to stop it being used by
accident.

**② Deny `Edit` on a sealed document.**

```json
"deny": ["Edit(./runs/**/02_prediction.md)", "Edit(./runs/**/01_intake.md)", ...]
```

The prediction is **written** by code, and the agent is blocked from
**amending** it by text edit afterwards. This closes the easiest route to post-hoc
rationalization structurally. Seal verification (`SEALED.sha256`) would catch it
after the fact, but **not being able to do it in the first place is better.**

## What this layer does not do

- Does not produce numbers — every one is the return value of a core call
- Does not fill `confirmed_by` — humans only
- Does not read the predecessor master plans in full (`docs/history/`, 6,700
  lines) — only the linked section
- Does not restate physics — it cites `knowledge/wiki/` and the system card
