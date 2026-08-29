# 01 · Architecture

Two things are layered here, and confusing them makes the repository look
redundant when it is not.

- **Layers (L1–L4)** are *who is allowed to do what*. This is a permission
  structure: the agent judges, the deterministic core computes, the knowledge
  base remembers, the run directory records.
- **Stages (S1–S8, and inside them L0–L7)** are *what happens in what order*.
  This is a pipeline with gates.

A layer never becomes a stage and a stage never becomes a layer. S5 (run) is
executed *by* L2 code *under* L1 judgment, writing into L4, reading from L3.

---

## 1 · The four layers

```text
┌──────────────────────────────────────────────────────────────────┐
│ L1  AGENT LAYER — .claude/ , CLAUDE.md                           │
│     Claude Code itself. Reading, judging, reasoning, asking,      │
│     recording. Performs NO deterministic computation; calls L2.   │
│                                                                   │
│     6 skills · 9 subagents (model-tiered) · 4 rules               │
└──────────────────────────┬───────────────────────────────────────┘
                           │  python -m bdbot.cli … / import
┌──────────────────────────▼───────────────────────────────────────┐
│ L2  DETERMINISTIC CORE — bdbot/ (engine) , simbot/ (S2·S7·S8)     │
│     Pure Python. No LLM. Exhaustively covered by pytest.          │
│     Units, non-dimensionalization, analytic solutions, HOOMD       │
│     execution, analysis, plotting, reports.                       │
│                                                                   │
│     grep -r "anthropic\|claude" bdbot/ must come back empty.      │
└───────────┬──────────────────────────────────┬───────────────────┘
            │ read / write                     │ write
┌───────────▼──────────────────────┐ ┌─────────▼────────────────────┐
│ L3  KNOWLEDGE — knowledge/       │ │ L4  ARTIFACTS — runs/<id>/    │
│     System cards, parameter       │ │     Every stage's artifact,   │
│     provenance, modelling         │ │     trajectory, figures,      │
│     rationale, failures,          │ │     report. Self-contained    │
│     verification benchmarks.      │ │     and reproducible.         │
│     ← accumulates, version-       │ │     ← text tracked,           │
│       controlled                  │ │       binaries gitignored     │
└──────────────────────────────────┘ └──────────────────────────────┘
```

**L3 is the actual asset.** L2 can be rewritten. *Why this material got
`η = 1.2 mPa·s`* and *why `dt = 5e-5` exploded* disappear unless they are
accumulated.

Why the L1/L2 split is load-bearing and not architectural taste: without it,
when a result is wrong you cannot tell whether the physics was wrong or the
model reading it was. In this domain that is fatal, because **a
plausible-but-wrong `g(r)` is indistinguishable by eye** — the simulation always
produces *some* number, and not diverging is not the same as being right.

---

## 2 · The pipeline

```text
   INPUT   a hand sketch · a note · a paper · a text description
     │
     ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S1  READ            transcribe → observation / inference / assumption   │
│                     three-way split, plus ambiguity candidates          │
│     L0  intake/<case>/observation.yaml                                  │
│     R   knowledge/wiki/systems/  — does a card for this system exist?   │
│     GATE  dimension · boundary · driving fixed;                         │
│           is the `question` falsifiable?                                │
│           back-translate the spec to prose and have a human approve it  │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S2  PREDICT         write the answer down BEFORE running, and seal it    │
│     simbot.estimators → prediction.yaml → SEALED.sha256                 │
│     GATE  ≥1 quantitative prediction, each with `tolerance` + `basis`   │
│           + a ROLE (implementation_check / hypothesis / measurement)     │
│                                                                          │
│     Sealing is structural, not disciplinary: settings.json refuses to    │
│     edit a sealed document, and a broken seal means the comparison       │
│     table is not built at all.                                          │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S3  SPECIFY         SI physical system, every number with a provenance   │
│     L2  bdbot.physical · bdbot.provenance · simbot.spec · simbot.cutoff │
│     R   knowledge/  — parameters, material properties, past decisions   │
│     GATE  no empty field · tier + derived_from present ·                │
│           derived values recomputed and matched                         │
│                                                                          │
│     BLOCKED here is a success: it names the ONE missing input.          │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S4  NON-DIMENSIONALIZE                        ★ the L2↔L4 contract      │
│     L3  bdbot.nondim → specs/<run_id>.json     — the ONLY thing passed  │
│     GATE  round-trip error < 1e-12 · dt constraints recorded ·          │
│           ledger complete (4 required roles) ·                          │
│           every dimensionless group really is a ratio of two ledger     │
│           entries · the inverse-transform anchor holds                  │
│                                                                          │
│     Self-sufficiency is decided by `nondim show <run_id>`: if the whole │
│     report draws from the spec alone, the spec is self-sufficient.      │
│     The health layer NEVER imports case code — it reads only this.      │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S5  RUN             seed ensemble, batched                              │
│     L5  @RUN.builder("<case>") → Build      (assembly contract)         │
│     L6  bdbot.run.execute()                 (equilibrate + production)  │
│     L7  bdbot.run.execute()                 (artifact storage)          │
│     L4  bdbot.health — Guard · judge · step_health                      │
│     GATE  pre-run gate: hash mismatch · FAIL · L3 integrity ERROR       │
│           in-run guards: NaN · displacement · configurational temp      │
│                                                                          │
│     A case supplies ONLY build(spec) -> Build. The equilibration and    │
│     production loops, guard calls and metrics.json are common code.     │
│     Read skill `bd-hoomd` before writing a line of this.               │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S6  VISUALIZE       the mandatory diagnostic set                        │
│     simbot.viz → figs/ + 06_figures.md                                  │
│     GATE  every figure has a caption and a `shows` field —              │
│           an uncaptioned figure CANNOT BE CREATED (34 tests)            │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S7  VALIDATE        prediction vs measurement                           │
│     simbot.validate · simbot.analysis · bdbot.metrics                   │
│     VERDICT  PASS / FAIL / INCONCLUSIVE                                 │
│     GATE  the seal verifies, or no comparison table is produced         │
│                                                                          │
│     INCONCLUSIVE is a first-class verdict: "this design had no power    │
│     to decide this" is different from "the prediction was wrong", and   │
│     design power is computed in S2 so you know in advance which items   │
│     will land there.                                                    │
│                                                                          │
│     A `hypothesis` mismatch is NOT a failure — it is the result.        │
└──────────────────────────────┬──────────────────────────────────────────┘
                               ▼
┌─────────────────────────────────────────────────────────────────────────┐
│ S8  CONCLUDE        answer S1's question; separate the assumptions      │
│     simbot.report → REPORT.md · tools/postmortem.py → record.json       │
│     W   knowledge/  — findings, benchmarks, dead-ends, tooling lessons  │
│     GATE  the conclusion answers the S1 question directly;              │
│           confidence and limits are stated                              │
└──────────────────────────────┬──────────────────────────────────────────┘
                               │
                               ▼
                    ╔═══════════════════════════════════════════════════╗
                    ║  KNOWLEDGE BASE                       knowledge/  ║
                    ╠═══════════════════════════════════════════════════╣
                    ║  wiki/systems/      (system × target dynamics)    ║
                    ║                     cards. THE CARD OWNS the      ║
                    ║                     non-dimensionalization and    ║
                    ║                     the gates — 11 cards          ║
                    ║  wiki/findings/     Q→A, and dead-ends — 23       ║
                    ║  wiki/benchmarks/   systems with known answers,    ║
                    ║                     run as regression tests — 5   ║
                    ║  wiki/concepts/ techniques/ questions/            ║
                    ║  source/papers/     42 distillations              ║
                    ║  source/books/      2 book distillations          ║
                    ║  entries/           126 tool-written JSON          ║
                    ║  runs/**/record.json  227 post-mortems            ║
                    ╚═══════════════════════════════════════════════════╝
                               │
                               └──▶ R  feeds the next S1 and S3
```

### When the knowledge base is read and written

| Moment | Direction | What moves |
|---|---|---|
| S1 begins | **R** | is there a `wiki/systems/` card for this system? If so it owns the scales and the gates |
| S3 needs a parameter | **R** | material properties, past decisions, `source/` distillations — with tier and provenance |
| S2 sets a tolerance | **R** | `wiki/benchmarks/` — systems whose answer is known |
| a run finishes | **W** | `record.json` post-mortem, one lesson per run |
| a `hypothesis` mismatches | **W** | `wiki/findings/` — this is a result, not a failure |
| a run fails or a dead end is hit | **W** | `wiki/findings/` dead-end page, so the next attempt queries it first |
| a tool bites | **W** | `entries/` with `origin: tooling` — 44 of these exist, which is not a good sign about the tools |
| a paper or book is distilled | **W** | `source/papers/`, `source/books/` — with the numbers we could and could not reproduce |

---

## 3 · The L-layers inside S3–S7

`bdbot`'s internal numbering, which appears throughout the code and the Korean
design documents:

| | Module | Job |
|---|---|---|
| **L0** | [`bdbot/intake.py`](../bdbot/intake.py) | `Observation` schema + checks, derived from 5 real sketches |
| **L2** | [`bdbot/physical.py`](../bdbot/physical.py) | `PhysicalSystem` loader; tier, `derived_from`, recomputation of derived values |
| **L3** | [`bdbot/nondim.py`](../bdbot/nondim.py) | ⭐ `NondimSpec` — the only contract between L2 and L4 |
| **L4** | [`bdbot/health.py`](../bdbot/health.py) | numerical-health verdict: `Guard`, `judge`, `step_health` |
| **L5** | [`bdbot/run.py`](../bdbot/run.py) | `@RUN.builder` registry + `Build` dataclass |
| **L6** | [`bdbot/run.py`](../bdbot/run.py) | `execute()` — equilibration and production loops |
| **L7** | [`bdbot/run.py`](../bdbot/run.py) | `execute()` — artifact storage, `metrics.json` |

Supporting modules: `units` (one pint registry — mixing registries makes pint
refuse), `provenance`, `materials` (`γ=3πηd`, `D_t=kT/γ`, `τ_B=d²/D_t`),
`pairpot` (the numbers that set `dt`), `scales` (`ScaleLedger`), `checks`,
`traps` (one harmonic trap expresses fixed, constant-velocity and oscillatory
driving), `sim`, `lockin` (complex stiffness `K*(ω)`), `report`, `runid`,
`metrics`, `stats` (block averaging with autocorrelation correction),
`interactions`, `cli`.

**What is deliberately NOT in `bdbot/`**, because it has appeared only once or
differs per system: equilibrium criteria · observables · verification strategy ·
choice of governing timescale · initial placement · sampling loop. Those stay in
the case scripts. The promotion rule is only ever *"has it appeared twice?"*

---

## 4 · Running it

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```

The front end enforces the order — you cannot skip a stage by accident:

```bash
$PY -m bdbot.cli status                      # where am I blocked
$PY -m bdbot.cli intake init  intake/<case>  # template → fill from the image
$PY -m bdbot.cli intake check intake/<case>  # FAIL / BLOCKED / READY
$PY -m bdbot.cli system check intake/<case>  # tier · derived_from · recompute
$PY -m bdbot.cli nondim spec  <case>         # L3 → specs/<run_id>.json (no run)
$PY -m bdbot.cli nondim show  <run_id>       # reproduce the report from the spec
$PY -m bdbot.cli health --gate specs/<run_id>.json   # pre-run gate
$PY -m bdbot.cli run <case>                  # L5→L7
$PY -m bdbot.cli health runs/<run_id>        # L4 verdict
```

`bdbot.cli` does not import `hoomd` or `freud`, so the front end is fast when
you are only reading and specifying. `bdbot.run` defers its `hoomd` import into
`execute()` for the same reason.

The S1→S8 half is driven from [`cli.py`](../cli.py):

```bash
$PY cli.py run      examples/trap-2d-5um/spec.yaml   # spec → REPORT.md
$PY cli.py resume   runs_s1s8/<id>                   # pick up a dead run
$PY cli.py converge examples/trap-2d-5um/spec.yaml   # shake dt · N · init
$PY cli.py calibrate                                 # measure this machine
```

⚠️ **This half has a runner for exactly one card**
(`passive-sphere--harmonic-trap`). Any other card needs a runner written into
`simbot.run` — it will say so rather than silently running your system through
the trap runner. That narrowness is the other face of the two-engine seam in
[00 §5](00-merge-decisions.md#5--known-seams): the 8 cases that actually ran
physics all went through `bdbot`, which has no sealing.

Tests:

```bash
$PY -m pytest -q -m "not slow"
```

589 tests collected; **572 pass** in ~11 s with 2 skipped, and `-m "not slow"`
deselects 15 long ones. Measured in this merged tree on 2026-08-28 — the merge
broke 6 of them (all in `test_agent_layer.py`, which guards the `.claude/`
structure) and those were fixed rather than deleted.
