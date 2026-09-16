# 08 · Microscope link — survey of the actual gap

**Status: survey only. Nothing is wired, no code changed.** This records what
was found on 2026-09-15 by reading both sides, against the plan in
[06 §6](06-roadmap.md). Its purpose is to replace the roadmap's precondition
table — written earlier — with what the trees actually contain today.

Sources read: this repo at `microscope-link-survey` (= `main`, `a18e171`), and
the read-only mirror `~/.claude/knowledge/microscope/` (synced 2026-09-15 18:22,
`kb/` + `docs/` only — no code on that side is mirrored here).

---

## 1 · Corrections to the roadmap's precondition table

Three rows of [06 §6 "What has to be true first"](06-roadmap.md) are stale. The
direction of every correction is the same: **the gap is narrower than recorded,
but a new collision has appeared that the table does not name.**

| Roadmap row | What the tree actually shows |
|---|---|
| "the 8 `bdbot` cases do not go through" sealing | **1 of 8 now does.** `bdbot/runcard.py` has `write_seal` / `SEALED_DOCS = ("prediction.yaml", "analysis_plan.yaml")`, and `bdbot/runid.py:26` refuses to run without a seal. But `grep -l runcard cases/*.py` returns **`cases/soft_r3_2d.py` only**. 21 `SEALED.sha256` exist on disk and every one is soft-r3 or a `runs_s1s8/` run |
| "Computed values have a provenance kind of their own — they do not, on either side" | **That side already writes `evidence: computed`** — 3 uses, all in `kb/calibrations/illumination-power.yaml`. See §2: this makes the situation *worse*, not better |
| "τ_c · ℓ_c have somewhere to live over there — not yet, `kb/samples/` arrives with its Phase 4" | `kb/samples/` still does not exist (`kb/` = calibrations, decisions, expertise, literature, plans, sessions, systems). But **the schema for the slot is already written** — `docs/02-knowledge-base.md:360` specifies `characteristic_scales.time` / `.length` with `evidence`, `method`, `model`, `inputs`, `measured_by`, `review_after`. The container is missing; the shape is not |

The other four rows stand as written.

---

## 2 · The one thing the roadmap did not foresee: two name collisions

The roadmap's stated hazard is *"a simulated number must never set
`evidence: measured`."* Correct, and easy to enforce. **The live risk is not
that. It is that a simulated number can already land in two existing slots that
mean something else, and nothing would flag it.**

**Collision A — `evidence: computed`.** On the microscope side this means
*arithmetic over already-measured kb values*. `illumination-power.yaml:270`
computes an illuminated area from a measured pixel size and a measured sensor
format. That is a measurement, propagated. A τ_c from a BD run is **not** that —
it is a model output whose falsifier is this repo's gate verdict. Written into
the same token, the two become indistinguishable at the point of reading, and
the distinction the roadmap wants to preserve is lost before any gate runs.

Note also that `computed` is **off-enum**: `docs/02-knowledge-base.md:366`
sanctions `measured | assumed` only, and the comment on that line says *"the
advances verdict looks only at this"*. So three values already sit outside the
vocabulary the verdict is specified against. Whether `Verdict.advances` treats
`computed` as not-`measured` (safe) or crashes is **not determinable from the
mirror** — that side's code is not synced here. → open question Q1.

**Collision B — `method: calculation`.** The same schema line offers
`measurement | calculation | literature | expert-judgment`, and its own worked
example uses `calculation` for a Stokes–Einstein one-liner. A BD result would
also be `calculation`. So the field that ought to carry *"how much work stands
behind this number"* cannot separate a back-of-envelope from a sealed, gated,
eight-stage run. **This is the field the interface actually needs, and it is
already occupied by a coarser meaning.**

Both collisions are cheap to fix now and expensive later — they are *additions*
to a vocabulary today, and *migrations* once values exist in the field.

---

## 3 · What each side has, in the same table

| | this repo | agentic-microscope |
|---|---|---|
| provenance carrier | `bdbot/provenance.py` → `Provenanced(value, source, tier)` | per-field YAML keys: `evidence`, `method`, `model`, `inputs`, `measured_by` |
| confidence scale | `tier: 0..3` (0 given/handbook · 1 literature+verification · 2 literature unverified · 3 assumption) | **no `tier` field.** The word appears 12 times in `kb/`, every one of them prose — and `kb/literature/README.md:32` says it outright: *"This repository has two evidence tiers, `measured` and `assumed`"*. The scale **is** the `evidence` token |
| units | `pint` quantities, SI (`bdbot/units.py`) | SI in the key name (`value_s`, `value_m`, `power_at_sample_mw`) |
| serialisation | `{value, unit, source, tier}` YAML node (`provenance.load_node`) | flat YAML, unit fixed by the key suffix |
| refusal | `BLOCKED`, naming the one missing input | `BLOCKED` / `Verdict.advances = false` |
| sealing | `SEALED.sha256`, plain `sha256sum` format, both packages | not applicable |

**The roadmap says "both already carry a tier and a falsifier on every stored
value". Half of that is wrong: only this side has a tier.** That side's
equivalent is the `evidence` token plus `review_after`. So the vocabulary job is
not "agree on a common serialisation for a tier both already have" — it is
**decide whether a tier crosses at all**, and if so which side grows a field.

The unit convention is the genuinely easy half, as recorded: both are SI, and
`value_s` / `value_m` maps onto a pint quantity without loss in either
direction.

---

## 4 · Where a number would actually cross

Only one slot is specified on both sides today:

```
bdbot case → τ_c, ℓ_c  ──▶  kb/samples/<system>.md
                              characteristic_scales.time.value_s
                              characteristic_scales.length.value_m
```

Consumers of that slot, per [06 §6](06-roadmap.md): **G8** (motion-blur ceiling,
wants `D` or τ_c), **G5** (wants ℓ_c + task kind), **G11** (target error),
**G14** (trap stiffness κ).

**But this repo does not currently emit either quantity under those names.**
`grep -rn "tau_c\|ell_c" bdbot/ cases/` finds only `tau_chain_diff` in
`cases/chain_bend_dlvo_2d.py` and `tau_chain` in `cases/chain_bend_2d.py` — both
case-local variables inside a ledger, neither a named export. `bdbot/report.py`
has exactly two public functions (`render`, `kv`) and both render text for
humans. **There is no machine-readable export surface at all**, and
`bdbot/cli.py` has no subcommand that emits one.

So the wiring the roadmap calls "small once those two exist" is, concretely,
three things and only the third is small:

1. name τ_c and ℓ_c as first-class per-case outputs (they exist as ledger
   entries, under different names, per case)
2. give `bdbot` a machine-readable export (there is none today)
3. write the file on the other side

---

## 5 · Open questions — not answerable from this tree

- **Q1.** Does `Verdict.advances` reject `evidence: computed` as not-`measured`,
  or does an off-enum token fall through? That side's code is not in the mirror.
  ⚠️ This is the unwired-checker shape this project keeps finding — an off-enum
  value sitting in three files while the verdict is specified against two tokens
- **Q2.** Were the 3 `evidence: computed` uses a deliberate vocabulary extension
  or a local improvisation? If deliberate, there is a decision record on that
  side that this survey has not read
- **Q3.** Is `soft-r3` sealing meant to generalise to the other 7 cases as-is,
  or is `SEALED_DOCS = (prediction.yaml, analysis_plan.yaml)` specific to the
  goal-first pipeline ([07](07-goal-first-pipeline.md))?

---

## 6 · What this survey does not change

The roadmap's stated **order** survives it: sealing first, then vocabulary, then
wiring. Nothing found here argues for reordering. What changes is the size of
step 1 (7 cases, not 8) and the content of step 2 — which is now known to
include **two collisions and a missing tier**, not just "a common serialisation".
