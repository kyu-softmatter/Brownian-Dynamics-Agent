---
name: bd-knowledge
description: Searches, extends and tidies the knowledge/ knowledge base. Use when asked "what is the basis for this parameter", "find me this value in the literature", "record this", "make a card", or when creating a new (system x dynamics) card, finding or benchmark. To run a simulation use bd-pipeline; to diagnose a failure use bd-diagnose.
---

# bd-knowledge — searching, extending and tidying knowledge

Contract: [`knowledge/wiki/CLAUDE.md`](../../../knowledge/wiki/CLAUDE.md)
← **read this first**

> **This wiki is not a general literature store.** It does two jobs only:
> ① supply the **verification oracle** for `pytest` regressions, and ② supply a
> **parameter dictionary with provenance**. A page that serves only as reading
> context has a weak claim to be here.

## Structure

```
knowledge/
├── source/papers/   42 per-paper distillations + INDEX.md   (original PDFs are gitignored; see NOTICE.md)
├── source/books/    2 book distillations
├── entries/         126 tool-written JSON entries (tools/kb.py)
└── wiki/
    ├── systems/     ★ (system × target dynamics) cards — these OWN the non-dimensionalization and the gates
    ├── findings/    Q→A plus dead-end-<slug>.md
    ├── concepts/    WHAT-IS (dimensionless groups, phase behaviour, potentials)
    ├── techniques/  HOW-TO (equilibration criteria, error bars, environment history)
    ├── benchmarks/  the verification oracle
    └── questions/   not answered yet. Never deleted — closed via `status`
```

⚠️ **`wiki/` and `entries/` are two unmerged schemas.** `wiki/` is human-written
Markdown with frontmatter contracts; `entries/` is tool-written JSON keyed by
`origin` × `kind`. They are read by different tools, so a lesson filed in one is
invisible to a reader of the other. Until they are unified, **query both.**
→ [docs/03](../../../docs/03-knowledge-base.md#1--two-schemas-not-yet-one)

## Searching

```bash
ls knowledge/wiki/systems/                      # is there a card for this system
grep -rl "<keyword>" knowledge/wiki/
grep -rl "<keyword>" knowledge/source/papers/
sed -n '1,20p' knowledge/source/papers/INDEX.md

<PY> tools/kb.py query --tags <tag> --origin tooling --kind pitfall
<PY> tools/kb.py lessons
```

**Read the frontmatter first** — `status`, `reproduced` and `confirmed_by` decide
whether that page can be used as a basis at all.

## Citation discipline — break it and the wiki becomes a rumour store

| Marking | Meaning | When |
|---|---|---|
| `[source]` | **verified basis** | a literature benchmark, or `reproduced: yes` |
| `[source, not reproduced]` | **a record of fact** | `reproduced: no` — consult it, but never use it in a verification claim |

- **Do not cite a literature value from memory.** Cite the distillation in
  `source/papers/`. "I saw it in training data" is not a basis
- **Do not use `reproduced: no` as a basis.** That a paper printed a number does
  not mean the number works in our code
- ⚠️ 38 of the 42 distillations are the group's own published work, so the
  literature layer is narrower than it looks. Weight it accordingly
- ⚠️ Bibliographic data here is largely unchecked — one distillation carries an
  explicit `verified: false` and 41 carry no `verified` field at all. The physics
  was checked; the volume-and-page line mostly was not. Confirm before a citation
  reaches a manuscript

## Creating a new (system × target dynamics) card

**If you meet a pair with no card, improvised non-dimensionalization is
forbidden.** `simbot.nondim` raises.

```bash
cp knowledge/wiki/systems/_TEMPLATE.md \
   knowledge/wiki/systems/<system>--<dynamics>.md
```

frontmatter:

```yaml
type: system
system: passive-sphere | abp | attractive-colloid | brush-colloid | interfacial-colloid
dynamics: equilibrium-structure | transport | harmonic-trap | dense-collective | coarsening
status: draft            # draft -> usable -> validated. Promote only on measured evidence
```

**What the card owns** (and is therefore written nowhere else):

| § | Content | Why the card owns it |
|---|---|---|
| 3 | **reference units** (length, energy, time) | the same system needs different ones depending on the target dynamics |
| 4 | the dimensionless-group ledger | only the ones meaningful for this pair |
| 6 | observables | |
| 7 | **which gates apply — on and off** | apply an equilibration criterion to an active system and it can never pass |
| 8 | benchmarks | candidates for promotion into `pytest` |
| 10 | remaining blanks | a closed item is never deleted, only re-marked |

Making a card also means **registering it in code**:

```python
# simbot/nondim.py
CARD_SCALE_RULES = {
    "<system>--<dynamics>": "harmonic_trap" | "brownian" | "active_run_length",
}
```

Without the registration `scales_for` raises — **that is the intended
behaviour.** (The same shape of gap bit the merge: `cases/network_3d.py` existed
and had produced runs, but was missing from `bdbot/cli.py`'s `CASE_SCRIPTS`, so
the CLI reported "no end-to-end script" for a case that had one.)

## A new finding

Use `findings/_TEMPLATE.md`. **The diagnostic path is this document's central
value** — next time the same symptom appears, someone follows that order exactly.

**A recurring "prevention: none" is the signal to build a guard.**

## A new benchmark → promotion into pytest

```markdown
| # | Benchmark | Expected | Status |
|---|---|---|---|
| B1 | dimensionless equipartition, EM bias | `⟨x*²⟩ = 1/(1−dt*/2)` | `[O]` measured `1.01041±0.00264` (`0.1σ`) |
```

When promoting into `tests/`:

```python
@pytest.mark.benchmark
def test_B1(...):
    """Take the tolerance from the **theoretical statistical error**.
    Tailoring it to the observed value is post-hoc rationalization, not
    verification."""
```

★ **When comparing against a documented value, take the tolerance from the
significant figures printed in the document.** A single global constant lets the
loose case excuse the strict one.

★ **Reject the competing hypothesis too.** "It agrees with the prediction" is
weak on its own. But do not demand a `3σ` rejection where the design power cannot
produce `3σ` — in that regime, fix `INCONCLUSIVE` as a fact.

★ **Give the benchmark a `role`.** A benchmark whose prediction comes from the
model we implemented is an `implementation_check`, and a mismatch is a bug. One
that comes from an assumption the simulation does not impose is a `hypothesis`,
and a mismatch is **a result to report**. Filing the second kind as a failure is
how a project files its discoveries as failures.

## Authorship symmetry and promotion

```yaml
author: agent | human | hybrid
confirmed_by:              # leave empty. A human fills it after review
```

- **The `author: agent` fraction is itself a self-improvement metric**
- **Promotion from `finding` to `concept` always needs human approval.** This
  stops an agent inflating its own observations into concepts
- The quality bar is the same regardless of author: the basis traces to
  `source/` or to established literature, the citation names a concrete file or
  URL, and the claim is checkable

## When a contradiction appears

**Never silently overwrite an existing entry.** Create a new one and link it with
`supersedes: [<previous id>]`. If it only narrows the scope, that is a
cross-link rather than a `supersedes` —
[`displacement-gate-is-1000x-loose-for-traps`](../../../knowledge/wiki/findings/displacement-gate-is-1000x-loose-for-traps.md)
**narrowed the scope of**
[`dt-gate-should-be-displacement-based`](../../../knowledge/wiki/findings/dt-gate-should-be-displacement-based.md)
rather than rejecting it.

## When the environment changes

If you add a package, record **why it was needed** in
[`techniques/env-log.md`](../../../knowledge/wiki/techniques/env-log.md) and
update `environment.yml` too.
