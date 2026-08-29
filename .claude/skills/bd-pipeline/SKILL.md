---
name: bd-pipeline
description: Runs the S1 to S8 pipeline from one piece of source material (a hand drawing, text, data) to a defended Brownian Dynamics conclusion. Use when the user hands over a drawing or photo and says "simulate this", "run this system", "measure the MSD", or when designing, running and verifying a BD / colloid / optical-trap / ABP system. To diagnose a run that already happened use bd-diagnose; to search or extend knowledge use bd-knowledge.
---

# bd-pipeline — S1 → S8

> **You judge; the core computes.** Do not produce a number in your head.
> Every value must be the return value of a function call.

## Read first, every time

1. [`CLAUDE.md`](../../../CLAUDE.md) — the absolute rules, the verdict
   conventions, the run environment
2. [`knowledge/wiki/systems/_index.md`](../../../knowledge/wiki/systems/_index.md)
   — is there a card for this system?
3. If there is, that card — **the card owns the non-dimensionalization and the
   gates**

Do not read: the predecessor master plans in `docs/history/` in full. Follow a
link to the section you need. The current design is
[`docs/01-architecture.md`](../../../docs/01-architecture.md).

### Domain skills — read them **before writing code**, at each stage

These three do not compete with the pipeline. The pipeline decides *when* to do
what; they decide *how* to do it without being silently wrong. They arrived with
the 2026-08-28 merge.

| Stage | Skill to read | Why |
|---|---|---|
| **S1** read | **`bd-intake`** | transcribe first · state ambiguity · absent means `null`. 8 anti-invention rules |
| **S3·S4** specify, non-dimensionalize | **`bd-physics`** | how to write a scale table · choosing the reference scale · scale-separation checks · inversion |
| **S5** run | **`bd-hoomd`** | 20 HOOMD traps. **Several produce wrong results with no error at all** (+1856 % error / `angle.Harmonic`'s force up to 96 % wrong while its energy is exact) |

⚠️ **If S5 writes even one new line of HOOMD, read `bd-hoomd` first.**
This is not advice — every trap on that list was actually walked into, and some of
them cannot be caught by an energy check.

## Interpreter

```
/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```

`conda activate` is unreliable in a non-interactive shell. **Use the absolute
path.**

---

## Stages and gates

Each stage **does not advance if its gate does not pass.** Go back, or report to
the user.

| Stage | What it does | Tool | Gate |
|---|---|---|---|
| **S1** read | source → observation / inference / assumption, plus ambiguity candidates | you (Opus) | dimension, boundary, driving fixed; is the `question` falsifiable |
| **S2** predict | write the answer down **before** running, and seal it | `simbot.estimators` | ≥1 quantitative prediction, each with `tolerance`, `basis` and a **role** |
| **S3** specify | `spec.yaml` with provenance attached | `simbot.spec` | no empty field, plausibility checks pass |
| **S4** non-dimensionalize | card scales plus `dt` | `simbot.nondim`, `bdbot.nondim` | round-trip error `< 1e-12`, `dt` constraints recorded |
| **S5** run | seed-ensemble batch | `simbot.run`, `bdbot.run` | completes, no guard violation |
| **S6** figures | the mandatory diagnostic set | `simbot.viz` | captions and dual axes all present |
| **S7** validate | verify the seal → comparison table → **propose** a verdict | `simbot.validate` | seal intact, every FAIL has a cause category |
| **S8** conclude | answer the question, commit knowledge | you (Opus) | at least one knowledge entry |

**S3–S8 are one command:**

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>
```

What you do by hand is **writing S1 and S2**, and **interpreting S7 and writing
S8**.

⚠️ **Know which engine you are on.** `cli.py` has the sealing and the
`INCONCLUSIVE` verdict but a runner for **one card only**
(`passive-sphere--harmonic-trap`). `bdbot` runs all 8 cases and has the L4 health
verdict but **no sealing**. They do not call each other. If your system is one of
the 8 cases, you are on `bdbot` and S2's seal is currently manual —
say so in the conclusion rather than implying it was sealed.
→ [docs/00 §5](../../../docs/00-merge-decisions.md#5--known-seams)

---

## Procedure

### 1. S1 — read the source → the intake document

Protocol: [`references/s1_intake_drawing.md`](references/s1_intake_drawing.md)
← **read it**

The core of it: **do not trust absolute sizes in a hand drawing.** Trust
topology, ratios, counts, symmetry and **explicitly written numbers**. For an
ambiguous element, **state 2–3 candidates** instead of picking one, and predict
how the result differs between them.

### 2. Check for a card

```bash
ls knowledge/wiki/systems/
```

- there is a card → follow its reference units and its gate table exactly
- **there is no card → improvised non-dimensionalization is forbidden.** Make a
  `status: draft` card from `_TEMPLATE.md` first and register it in
  `simbot/nondim.py::CARD_SCALE_RULES`

### 3. S2 — predict → `prediction.yaml`

Protocol: [`references/s2_prediction.md`](references/s2_prediction.md)

Get the numbers by **calling** `simbot.estimators`. Example:
[`examples/trap-2d-5um/prediction.yaml`](../../../examples/trap-2d-5um/prediction.yaml)

★ Give every prediction a **role**. If all of them are `implementation_check`,
this case can validate but cannot discover — say so up front rather than
discovering it at S8.

### 4. S3 — specify → `spec.yaml`

Protocol: [`references/s3_s5_execute.md`](references/s3_s5_execute.md)

Copying the example and editing it is fastest:
[`examples/trap-2d-5um/spec.yaml`](../../../examples/trap-2d-5um/spec.yaml)

**Every value gets a `provenance` and a `basis`.** Check:

```bash
<PY> -c "
from simbot.spec import SystemSpec, validate
r = validate(SystemSpec.load('<spec.yaml>'))
print(r.table()); print(r.problems or 'no convention violations')"
```

### 5. Check the cost → run

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>
```

The CLI **stops before running** on a budget overrun or a gate violation. Use
`--force` only when the user explicitly asks.

To shake a parameter and see only the cost, **without running**:

```bash
<PY> -m simbot.session new <spec.yaml>
<PY> -m simbot.session set numerics.dt_star=2.5e-3 species.0.n_simulated=4000
```

On the `bdbot` side, the pre-run gate is separate and **shows warnings without
blocking**:

```bash
<PY> -m bdbot.cli health --gate specs/<run_id>.json
```

### 6. Interpret S7

Protocol: [`references/s6_s7_validate.md`](references/s6_s7_validate.md)

The CLI **proposes** the comparison table and the verdict. Your job:
- was the `INCONCLUSIVE` **foreseen**, or is it a design mistake?
- the cause category of a `FAIL`: `numerical` / `modeling` / `interpretation` /
  `analysis` — and **check the role first**, because a `hypothesis` mismatch is
  a result, not a fault
- if a `PASS ⚑` appears (a significant deviation still inside tolerance), check
  **whether the prediction was sloppy**

### 7. S8 conclude → the conclusion document plus knowledge

Protocol: [`references/s8_knowledge.md`](references/s8_knowledge.md)

`REPORT.md` is generated by the CLI. **What you write is the answer to the
question, the causal hypotheses, and the next experiment.**

---

## Never do this

| ❌ | Why |
|---|---|
| calculate a number in your head | you will never know whether the physics was wrong or you were |
| fill `confirmed_by` | humans only. Otherwise a pass gets stamped that nobody looked at |
| put a number without an error bar into a conclusion | a single-seed production run is forbidden |
| say "verified" or "reached equilibrium" | those need a threshold to be sayable |
| cite a literature value from memory | cite the distillation in `knowledge/source/papers/` |
| use `reproduced: no` as a basis | that is a record of fact, not a basis. Mark it `[source, not reproduced]` |
| improvise non-dimensionalization for a card-less pair | `nondim.py` raises — do not route around it |
| call a `hypothesis` mismatch a failure | that is the result. Reporting it *is* the job |
| ask for an API key | the LLM doing the reading is already in this session |

## Question budget — three per round

**If the conversation becomes twenty questions, the tool has failed.** Ask only
when one of these three holds:

1. the value flips the conclusion (expected sensitivity `|S| > 1`)
2. it is near a regime boundary (where the interpretation forks)
3. not even the order of magnitude is known (no basis in knowledge or literature)

Otherwise fill it as `provenance: assumed`, **show the proposal table**, and
proceed. A proposal always carries three parts: **value · basis · confidence**.
And tell the user **where to change it**.

## When it fails

Switch to the `bd-diagnose` skill. And **record the failure** — delete it and half
the evidence about whether this agent works at all disappears.
