# Brownian-Dynamics Agent — project contract

Read this every session. It is the contract, not the design document — the
design is in [`docs/`](docs/), the results are in
[`docs/04-cases.md`](docs/04-cases.md), and the accumulated judgment is in
[`knowledge/`](knowledge/).

**One line.** Read a physical system out of a hand sketch, a note or a paper;
non-dimensionalize it; run it in HOOMD-blue; and feed both the successes and the
failures back as knowledge the next run can query.

---

## Run environment

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```

conda env `simulation_bot` · HOOMD-blue **7.1.0** (CPU, **no MPI**, no GPU) ·
macOS arm64. Defined in [environment.yml](environment.yml).

`conda activate` is unreliable in a non-interactive shell. **Use the absolute
path.**

---

## The division of labour

> **You judge. The code computes.**

| The agent does | The agent never does |
|---|---|
| Read a sketch, a note, a paper into a draft spec | Arithmetic in its head |
| Choose which model applies, and record why | Decide whether a check passed |
| Triage a failed verification into causes | Fix a parameter by fiat |
| Write the narrative of a report | Originate a physical conclusion |

Every number in a report has to be the return value of a function in
[`bdbot/`](bdbot/) or [`simbot/`](simbot/). This split is not tidiness: without
it, when a result is wrong you cannot tell whether the physics was wrong or the
model reading it was. A plausible-but-wrong `g(r)` is indistinguishable by eye,
which is what makes that ambiguity fatal here.

---

## Absolute rules

**1 · Dimensions come first.** Fix every system as an SI physical system, then
pass it through a scale table to non-dimensionalize. There is no path that
starts from dimensionless values. If you must start from one, state the anchors
(particle diameter, temperature, viscosity) explicitly.

**2 · Never hand-write a dimensionless spec.** Only run what was derived from a
physical system. `specs/` and `runs/` are written **by tools only**
(`bdbot.cli nondim spec`). A hand-edited spec is caught by
`LoadedSpec.verify_hash()` — `run_id` is the hash of the content, so editing the
content breaks it.

**3 · Every parameter carries a provenance.** Sketch, literature, handbook or
estimate — say which. If you do not know, say you do not know. **Never invent a
value.** `BLOCKED` is not a failure: it means you have narrowed the missing
input to one thing. Forcing it to `READY` by making a number up is the failure.

**4 · Never write HOOMD from memory.** Several traps produce wrong results with
no error at all. Read skill **`bd-hoomd`** *before* writing the code.

**5 · When reading a sketch, transcribe before you interpret.** Write down what
is visibly there, then interpret. State what you could not read and what is
ambiguous. Anything absent from the sketch stays `null`. Skill **`bd-intake`**.

**6 · Verify a physical claim before stating it.** This project has been wrong
repeatedly: a missing minimum image gave **+1856 %**; "bind it rigid and it
becomes anisotropic" measured **1.000000**; a `pair.Table` grid convention put
the force off by **−1.65 %**; a cage-fluctuation formula was missing a
`√2·a_NN`. Claims about HOOMD behaviour and about physics formulas are settled
by execution, not by reasoning.

**7 · Isolate independent elements and verify them one at a time.** If a system
has several mutually independent elements, first run the **minimal
configuration with only one of them on**, compare against the analytic
solution, and only then combine. Often only the isolated configuration *has* an
analytic solution. If the comparison fails, **split further**.

**7′ · "Verified" does not mean "agrees with the prior hypothesis."** ⭐
The reason to compute these systems is that they may differ from the standard
picture. Call a mismatch a failure and you will be calling discoveries
failures. So every comparison carries a **role**
([`bdbot/metrics.py`](bdbot/metrics.py) `ROLES`):

| Role | Where the prediction comes from | A mismatch means |
|---|---|---|
| `implementation_check` | **Derived from the model I implemented** | **A bug** → fix it |
| `hypothesis` | An assumption the simulation does **not** impose (continuum, dilute limit, effective medium, literature) | **A result** → report it |
| `measurement` | No prediction | The simulation is the answer |

When designing a case, write down separately *what the simulation imposes* and
*what the theory adds*. If the second list is empty the case can only validate,
never discover — `abp-rod` was exactly that (5 predictions, all
`implementation_check`, zero hypotheses).

*(Two rules numbered 7 is convention, not an error — cite the isolation one as
"rule 7" and this one as "rule 7′". There are already six references to them.)*

**8 · Build the static system first; add motion afterwards.** ⭐ Before any
time-dependent term (driving, oscillation, dragging), build the static
configuration — equilibrium arrangement, bond structure, energy minimum — and
finish verifying *that*. Static problems have analytic solutions and exact
minimizations, so they can be ground truth; dynamic responses usually cannot.
Noise is zero, so systematic error is not buried. And the static structure is
the initial condition of the dynamic measurement, so inverting the order mixes
two errors inseparably. Store the static stage as its **own run** so it can be
re-run alone.

**9 · Build systems out of direct inter-particle interactions — never
substitute a macroscopic general formula.** If an interaction is physically
needed, put it in as a pair potential, bond or angle that actually acts between
particles. **The reverse is equally forbidden**: do not shortcut to a
macroscopic property (modulus, diffusivity, stress) of an already
particle-resolved system by substituting into a continuum, mean-field or
literature relation (GSER, effective-medium theory, continuum beam formulas).
Those relations have preconditions — continuous medium, dilute limit, rigid
joints — and the systems here break them often. **Get the macroscopic answer by
running the inter-particle interactions**, not by inverting a formula.

Worked example: `GSER  G* = K*/(6πa)` was **invalid** for `chain-bend`, because
what the bead feels is not a continuous medium but **two neighbouring beads**.

---

## Working practice

- **Case-driven.** Do not build the framework first. Drive one case end to end,
  and abstract only what has appeared **twice**. What is deliberately still
  duplicated: equilibrium criteria, observables, verification strategy, choice
  of the governing timescale, initial placement, sampling loop.
- **Verify together.** Scale tables and non-dimensionalization results get
  looked at with a human. Do not settle a physics judgment alone until the
  knowledge base has depth.
- **Keep what you learned even when there is no run.** Sketch readings,
  literature distillations and tooling lessons attach to no run:
  ```bash
  $PY tools/kb.py add --origin intake|paper|tooling|method|handbook \
      --kind pitfall --source "file#anchor" --claim "..."
  ```
- **Give every observable a role.** Default is `measurement` (no verdict). Only
  `implementation_check` mismatches are FAIL; `hypothesis` mismatches are
  **reported as results**.
- **Show results as graphs and animations.** ⭐ Do not report in prose and
  tables alone — the default is to produce both and attach them.
  - **Graphs**: measurement vs prediction on the same axes, with check
    thresholds, analytic solutions and literature values overlaid. What is
    *wrong* should be visible.
  - **Animations**: `kT=0` (mode shape) and `kT>0` (how far thermal noise
    buries it) **side by side** shows an SNR problem faster than any number.
    Cheap large-`dt` animations are fine but **must be labelled as not being
    the production measurement**.
  - ⚠️ **Label graphs in English, not Korean.** matplotlib's default
    `DejaVu Sans` has no Hangul, so labels render as `□`; the fonts that do
    have Hangul (`AppleGothic`, `Apple SD Gothic Neo`, `NanumGothic`) are
    missing `−` (U+2212) and `ŷ` (U+0177) — measured. Do not try to fix it by
    switching fonts. Write axes, legends, titles and annotations in English
    from the start, and confirm zero `missing from font` warnings.
- **After a run finishes, always run the post-mortem.** Do not declare success
  or failure — decide it by measurement:
  ```bash
  $PY tools/postmortem.py runs/<run_id> --lesson "lesson::pitfall::coord=value"
  ```
- **Leave verification scripts in [`verify/`](verify/).** They must be
  reproducible.
- **When you build a checker, deliberately break it and see.** "Silently
  passing" and "not checking" are different things —
  `verify_intake_guards.py` actually caught a crash bug this way.
- **Found a new trap? Add it to skill `bd-hoomd`** and leave a reproduction
  script.

---

## Where to read what

| When you need | Read |
|---|---|
| The architecture and the stage contracts | [`docs/01-architecture.md`](docs/01-architecture.md) |
| How verification is assembled | [`docs/02-verification.md`](docs/02-verification.md) |
| Reading a sketch or an image (**always first**) | skill **`bd-intake`** |
| Writing HOOMD code (**always first**) | skill **`bd-hoomd`** |
| Defining a physical system, proposing parameters, non-dimensionalizing | skill **`bd-physics`** |
| Running the S1→S8 pipeline | skill **`bd-pipeline`** |
| Diagnosing a broken run | skill **`bd-diagnose`** |
| Searching or extending the knowledge base | skill **`bd-knowledge`** |
| Using shared code | [`bdbot/__init__.py`](bdbot/__init__.py) — what exists and what deliberately **does not** |
| What is buildable at all | [`docs/hoomd_capabilities.md`](docs/hoomd_capabilities.md) |
| Defining a stress / modulus / rheology observable | [`knowledge/source/books/leal_microstructural_rheology.md`](knowledge/source/books/leal_microstructural_rheology.md) |
| A citable material property | [`knowledge/source/books/welty_transport.md`](knowledge/source/books/welty_transport.md) |
| Past runs and lessons | `$PY tools/kb.py query --tags ... --origin ... --kind ...` |
