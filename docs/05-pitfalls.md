# 05 · Pitfalls

The full, maintained list is skill
[`bd-hoomd`](../.claude/skills/bd-hoomd/SKILL.md) — 20 traps, each with a
reproduction script. **Read that skill before writing HOOMD code**, not after.
This page is the argument for why that rule exists, plus the traps that are
about *this project* rather than about HOOMD.

**★ marks a trap that is silently wrong — no error, no divergence, a plausible
number.** Those are the dangerous ones, and most of the list is starred.

---

## 1 · Why "just be careful" does not work here

Four measured examples, all of which passed some verification before being
caught.

**★ External force + periodic boundary → minimum image is mandatory.**
Using `d = pos - anchor` directly means that the instant a particle wraps, the
distance jumps by `L` and it feels an enormous force in the *wrong direction*.
**It does not blow up.** And the error depends on trap stiffness:

| `k` | without minimum image | with |
|---|---|---|
| 2 | **+1856 %** ✗ | +0.38 % ✓ |
| 5 | +344 % ✗ | +0.56 % ✓ |
| 10 | +0.16 % ✓ | −0.02 % ✓ |

Test only the stiff condition and you will believe it passed. **Verify in the
weak condition.**

**★ `pair.Table`'s grid is `endpoint=False`.** The documentation says the
implicit `r` values are `numpy.linspace(r_min, r_cut, len(U), endpoint=False)`.
Build the table with `endpoint=True` and the whole table shifts, so **the force
is quietly wrong** — and worse near the cutoff:

| separation | `endpoint=False` | `endpoint=True` |
|---|---|---|
| 0.70 | +0.000 % | −0.572 % |
| 1.50 | +0.000 % | −1.329 % |
| 2.90 | +0.000 % | **−1.646 %** |

**★★ `md.angle.Harmonic`'s force is up to 96 % wrong while its energy is exactly
right.** It converts torque to Cartesian force through `1/sin θ` and clamps
`sin θ` at √2×10⁻³, so below that the force becomes `∝ κ(θ−π)²` — quadratic, not
linear. With `t0=π` the equilibrium itself is at `sinθ=0`, so it bites hardest on
**stiff, nearly-straight chains**. Energy is **0.0000 % correct throughout**.
This one broke a whole case ([04
`chain-bend-2d-oscill`](04-cases.md#chain-bend-2d-oscill--the-hoomd-bug-and-the-way-around-it)),
and no energy-based check could have found it.

**★ `update.BoxResize` collapses bonds in a narrow well.** It affinely scales
coordinates (error 8.9e-16), which shortens *already-bonded* pairs too. Past
0.703 % linear strain per trigger the pair is pushed inside the barrier and
**collapses irreversibly** into the primary minimum (0.40 % holds, 0.80 %
collapses). Measured before running, which is why the `network` compression uses
0.4 % over 178 steps.

Others in the same family: `pair.Table` gives **zero** force and energy for
`r < r_min`, so with a diverging potential a particle that once penetrates
`r_min` simply stays overlapped, without blowing up. There is **no WCA class** —
`ForceShiftedLJ` is a different potential, and the name will mislead you.
ABP rotation must go through an updater with `integrate_rotational_dof = False`,
or inertial rotation mixes in and **it quietly stops being ABP**. `Brownian`'s
`velocity` field is not zero but **uncorrelated thermal noise**, so
velocity-based MSD or VACF is meaningless in overdamped dynamics.

---

## 2 · Project-level traps

These are not HOOMD's fault. They are ours, and they cost more.

### An unwired checker cannot be wrong out loud

`step_health()` — documented in its own module as "the core of this module" —
**never ran, in all 81 runs**, because of a name mismatch, and returned HEALTHY
anyway. Separately, the pre-run gate **falsely rejected 80 of 83 specs**, with
zero real failures among them, and nobody noticed because `execute()` never
called `gate()`. Full account in
[02 §6](02-verification.md#6--the-failure-mode-this-document-exists-to-prevent).

**Rule: `N/N HEALTHY` is not coverage.** Print the count of *unmeasured* runs
separately, and state what the verdict does and does not cover.

### A crash while writing artefacts makes a passing run look failed

The first attempt at the S30 campaign printed

```
VERDICT: ✓ PASS
```

and then died. `make_plots` read `observables.npz["eq_trace"]`, which does not
exist when a run has no equilibration phase — and a melting measurement needs
`--eq-frac 0`, because that phase is built with `collect=False` and would leave
the first 20 `τ_B` of a 100 `τ_B` window unsampled.

The physics was complete and correct. But `result.txt` is written **after**
`execute()` returns, and `RID.prepare_outdir` treats its presence as the
completion marker, so:

| | |
|---|---|
| the run | finished, `metrics.json` written, every check passed |
| `result.txt` | never written |
| exit code | 1 |
| the campaign driver | read a correct run as a failed one and stopped |
| cost | 1 of 12 runs, and the whole queue behind it |

**Rule: the completion marker must not be downstream of an optional artefact.**
Two independent things go wrong together here — `result.txt` being written by
the case script rather than the engine (already filed as a seam, below), and a
figure panel that assumes a phase exists. Either alone is survivable; together
they convert a `PASS` into a `FAIL` with no diagnostic anywhere near the cause.

⚠ And the traceback is the *lucky* version. Had `make_plots` swallowed the
`KeyError` instead, the run would have been reported as complete with a panel
silently missing — which is the same defect as the family above.

Prevention: `tests/test_soft_r3_init.py::test_a_run_with_no_equilibration_phase_writes_its_artefacts`
asserts `result.txt`, `metrics.json`, `observables.npz` and `observables.png`
all exist after an `--eq-frac 0` run, with a sibling test requiring the default
path to still produce an `eq_trace` so the branch stays a branch. Unguarding the
read fails the first and not the second (measured).

### A gate that checks presence, not validity

Rule 10's manifest refuses a number with no value, no unit or no provenance.
[`bdbot/params.py`](../bdbot/params.py)'s docstring also said, from the day it
was written, that `derived` "is checked against a recomputation … the same
invariant one layer up." **No such check existed for twelve days.**
`blockers()` verified that each `REQUIRED_KEYS` name was *present* and never
that its value was *real*, so this cleared the gate:

```python
m.add("derived", "gamma",   0.0, "1", "6 pi eta a -- in the ledger, recomputed by L3")
m.add("derived", "D_t",     0.0, "1", "kT/gamma -- in the ledger, recomputed by L3")
m.add("derived", "tau_gov", 0.0, "1", "tau_B = d^2/D_t -- in the ledger, recomputed by L3")
```

`blockers()` returned `[]`. Three physically impossible zeros — no friction, no
diffusion, no timescale — in the first manifest anyone built, written for the
run that was meant to *demonstrate* rule 10. **The provenance string was the
confession** ("in the ledger") and nothing read it.

That docstring opens by listing three practices that were written and never
enforced — `bd-intake` §2.1's empty-goal blocker, `A4`'s grep, `health.gate()`
— and calls the record *three for three*. It was four, and the fourth was in
the same file.

**Rule: a required field is a schema check, not a physics check.** If a value
follows from other values in the same document, recompute it. If it does not
follow from anything, it is a choice and must say so — `unknown(name, who)` is
what "it is in the ledger" should have used.

#### …and the fix over-corrected within the hour

The recomputation was written with `rtol = 1e-6` and with `tau_gov` compared
against `tau_B`. It immediately **refused a correct manifest** —
`trap-2d-5um`'s card, on two counts:

| declared | recomputed | verdict | what was actually wrong |
|---|---|---|---|
| `gamma = 4.0102e-8 kg/s` | `4.010243e-8` | −0.00107 %, **blocked** | the tolerance. The value is right to the five figures it was *written* to |
| `tau_gov = 4.010e-3 s` | `tau_B = 242.051 s` | −100 %, **blocked** | the check. `tau_gov` is *which timescale governs* — here `τ_k = γ/k`, and the provenance says so in words |

Both are the same mistake in opposite directions: **a comparison's tolerance
comes from how the value was written, and a field's meaning comes from what it
was for.** Rule 10 asks for "γ, `D_t`, `τ_B`, and **which timescale governs**"
— four things, and the fourth is a case-dependent choice, not a synonym for the
third. `τ_B` is recomputed; `tau_gov` is checked for being a positive time.

A gate that refuses a correct answer is worse than no gate, and this repository
had already measured that once at 80 rejections out of 83. The second version
was caught by the mutation harness in the same session that introduced it
([`verify/verify_gates_bite.py`](../verify/verify_gates_bite.py)), which is the
only reason it is a footnote rather than a fourth entry.

### A check whose success is indistinguishable from its own failure

The family below is about checks that go blind. This one is worse: it inverts.

`grep -c` exits **1** when the count is zero. For a translation sweep, zero
matches *is* the success condition — so wiring it into an `&&` chain produces a
harness that acts on failure and refuses to act on success. Measured:

| file state | the check says | `grep -c` exit | `… && git commit` |
|---|---|---|---|
| no Korean left | **passed** | 1 | **stops — the commit never runs** |
| Korean still there | **failed** | 0 | continues — **the commit runs** |

This happened on 2026-08-29. The chain printed a plausible `0`, the commit step
silently never executed, and the next command pushed an unrelated session's HEAD.
The count was correct and clearly displayed; nothing downstream could act on it.

**Rule: never put `grep -c` in an `&&` chain.** Use `[ "$(grep -c PAT f)" -eq 0 ]`,
or `! grep -q PAT f`. And **confirm a commit happened by reading `git log`** rather
than inferring it from surrounding output that still looks like success.

★ Note where the defect sat: not in the check, but in **the harness around it**.
Same for the two below — `chk_doc` computed the right value and asserted on the
wrong proposition; the refactor audit resolved the right names and was scoped so
it saw none of them. A correct measurement wired to something that cannot act on
the answer is the recurring shape, and it is invisible to the measurement itself.

### A verifier that iterates its own vocabulary instead of the artefact

★★ **`simbot.io.verify_seal` hashed zero documents on all 21 archived seal
directories.** Measured 2026-09-15 by instrumenting `sha256_file` and counting
calls. Of those 21, **12 returned `ok=True`** while printing
*"봉인 검증 통과 — 2개 문서 실행 후 미변경"* — a pass claiming two documents were
verified unchanged, having hashed none. The other 9 returned `ok=False` with the
documents listed as `unsealed`, i.e. *"never sealed"*, while their content was
byte-for-byte intact.

The cause is one line of loop:

```python
for stage in stages:              # SEALED_STAGES: 02_prediction.md, 01_intake.md, 03_spec.yaml
    p = rundir.file(stage)
    rel = _seal_relpath(p)        # where the document is NOW
    if rel not in sealed:         # the seal records where it WAS
        if p.exists(): unsealed.append(...)
        continue                  # ← content never hashed
```

The verifier iterated **the filename vocabulary this module happens to know**
rather than **the seal's own entries**. Two independent consequences:

- `bdbot.runcard` seals `prediction.yaml` + `analysis_plan.yaml`. None of the
  three `SEALED_STAGES` files exists in those runs and none of their paths is a
  seal key, so every branch `continue`d, all three problem lists stayed empty,
  and `ok = not (changed or missing or unsealed)` evaluated to **`True` on an
  empty check**.
- The 2026-08-28 merge renamed `runs/` → `runs_s1s8/`. The seals still record
  `runs/<id>/…`, so the current repo-relative path is not a key → `unsealed`.
  **"we never sealed it" and "we sealed it and then renamed the directory"
  became the same output** — the two states the entire mechanism exists to
  distinguish.

★ **And it could not see tampering.** Because the branch taken was `unsealed`
rather than `changed`, the content was never read, so appending a falsified
prediction to one of those 18 documents produced a **bit-identical**
`SealVerdict`. Demonstrated, not argued.

**Why nothing caught it for 18 days.** Every unit test builds its run in
`tmp_path` and seals it in place, which is the one configuration where the stage
path and the seal key coincide. The single test that touched the archive —
`test_verify_seal_on_real_first_run` — probed
`runs/2026-07-28_trap-2d-5um_2dfb9d`, **a path absent from every commit in the
history** (the run has been under `runs_s1s8/` since the initial commit), so it
skipped on every platform with the reason *"runs/ is gitignored — not in this
checkout"*, which is also false: `git check-ignore runs` exits 1 and 1311 files
under `runs/` are tracked. A silent skip, on a false premise, at the centre of
the pre-registration machinery.

**Two implementations, complementary blind spots, and the broken one was wired.**
`bdbot.runcard.verify_seal` was correct throughout — it iterates the seal's
entries, has a rename fallback, and separates `[warn]` from hard problems. It is
called from `run.execute`. The broken one is what `cli.py`, `simbot/report.py`
and `simbot/validate.py` use. The repository deliberately keeps two (documented:
`bdbot` cannot import `simbot`), so the fix adds a test asserting **the two
agree on every archived seal** — if they may not be unified, they must not
diverge.

**The recorded path is provenance, not a resolution strategy.** CI's bash
preferred the recorded path and fell back to the sibling. That is unsafe in a way
the rename hid: a **copied** run directory carries the original's paths, so its
seal would verify the original's documents and pass while its own went unchecked.
A seal can only ever cover documents in its own directory (`write_seal` cannot
reach outside it), so the fixed rule is: **verify the sibling; report the
recorded path not matching as drift.** Both checkers now follow it and report
the same 42 documents and the same 18 drifts.

**And the external verification command did not work.** Reports advertise
`shasum -a 256 -c SEALED.sha256` under *"verified without this code"* — the
seal's only independent check. It resolves in **no** working directory: from the
run directory the repo-relative entries become `runs/<id>/runs/<id>/…`; from the
repo root the seal file is not there. Two committed reports carry it directly
beside *"✅ 봉인 검증 통과"*. It is generated by `report.seal_check_command` now,
and two tests **execute** its output over all 21 seals — one asserting it passes,
one asserting it reports `FAILED` on an edited document, because a command that
only ever passes is indistinguishable from `true`.

**The general form:** *iterate the artefact you are checking, not the names you
expect it to have* — and make a pass state **how much it checked**. `ok` now
requires a non-empty `verified` list, and the count is in the summary. After:
21/21 verify, 42 documents hashed, 18 drifts reported, both implementations
agree 21/21, and five gate mutants are CAUGHT
(`verify/verify_gates_bite.py`, `SEAL-*`).

### Output that looks like a finding

The unwired-checker failure above has a family. Each member produces something
that *reads* as a result while carrying no information, and the tell is always
the same: **the check cannot distinguish the case it is testing from the case
where it did not run.**

| Shape | What it cannot tell apart | Found |
|---|---|---|
| a doc-scraper whose regex stops matching | `0 matches` from `passed` | 2026-08-29 — a transcription typo (`0.8580` for `0.8598`) passed 56/56 because nothing tied the computed value to the printed one |
| `chk(want=0.0, rtol=1.0)` | a 100 % tolerance from the strictest possible | the bound is `rtol × max(|want|, 1e-300)`, so `want = 0` collapses it to `1e-300` whatever `rtol` says. Passed only because the quantities were bit-exactly zero |
| a refactor audit reading AST definitions | **moved** from **deleted** | flagged 5 files, of which 4 were `def`s that had legitimately become re-export `import`s — the merge working as designed |
| an interpolation quoted without its method | `+1.03 %` from `+2.91 %` | log-linear vs linear on the same 20 K table |

**Rule: make the ambiguous call impossible rather than auditing the call sites.**
`chk` now *fails* when `want == 0` and no `atol` is given, because a relative
tolerance on zero is not a loose bound, it is an incoherent request. Doing that
surfaced a fourth call site that enumerating the known three had missed. The same
move worked twice more the same day — an accessor that names the cause instead of
raising `KeyError` found five sites, and resolving module names *at runtime*
instead of parsing the AST told moved from deleted. Three independent instances,
one lesson: **fix the contract and it finds what a sweep does not.**

### When a copied constant is legitimate

The rule against duplicated tables has one real exception, and it is worth
stating precisely, because a de-duplication pass spawns duplicates while it runs
— which happened in both directions on 2026-08-29, including a **fourth** copy of
the water-viscosity table added *inside the verifier written to fix a
table-divergence bug*, and a copy inside the test written to assert there was
only one copy.

A copy is legitimate only if all three hold:

1. it is **labelled frozen**, with the commit it snapshots;
2. it is **read by nothing but the comparison**; and
3. it is **wrong-by-construction the moment someone updates it.**

That third clause is the whole point: you cannot prove a move changed nothing by
importing the thing you moved. If a physical value legitimately changes, the
comparison *should* fail and the change should be argued, not absorbed.
→ KB `tooling__a-de-duplication-pass-spawns-duplicates-while-it`

### Search-delimited edits need an asserted end marker

A mechanical replacement that locates its **start** by matching text and its
**end** by searching for the next plausible line will, when the end search
overruns, delete everything in between and leave valid syntax behind. Measured:
one such edit replaced **41 lines instead of 2**, removing a reader function
whole. `ast.parse` caught that one because the result happened not to parse; the
dangerous version is the one that does.

**Rule: assert `count == 1` on the exact text being replaced, and assert the end
marker as well as the start.** A start-only match will happily eat the rest of
the file.

### A checker scoped by a bare `git diff` goes blind exactly when it matters

`git diff --name-only` lists **unstaged** changes only. A verifier scoped that way
audits nothing the moment its subject is staged — which is to say, at the instant
someone is ready to commit it. Measured on 2026-08-29: the refactor audit in
[`verify/verify_merge_equivalence.py`](../verify/verify_merge_equivalence.py)
covered **zero files and printed `PASS`**, and it was noticed only because a human
happened to see the empty output. That is not a mechanism.

It took three attempts to scope correctly, and each wrong answer was found by
running it rather than by reading it:

| scope | what went wrong |
|---|---|
| `git diff --name-only` | unstaged only — audited nothing once its subject was staged |
| `git diff --name-only HEAD` | **failed on a clean clone**, since a fresh checkout has nothing uncommitted |
| "whatever is uncommitted" | in a shared worktree this picked up nine `verify/*.py` files another session was mid-translation on — and because those are *scripts*, importing them **ran** them. One went looking for a `.gsd` that does not exist |

The third is the worst of the three. A check that fails on clone is merely wrong;
a check that executes code it was never pointed at has side effects on work that
isn't its business.

**Rule: scope a checker to a pinned pair of commits — never to a range, and never
to the working tree.** A range like `BASELINE..HEAD` drifts as other sessions
commit; measured here, it already annexed three files belonging to an unrelated
translation pass. Take the *file list* from the pinned pair and the *content* from
the working tree, so an uncommitted fix on top still gets checked. Empty list is a
failure in every mode. In a tree with concurrent sessions, "what changed" is not a
question about you.

### The index is shared state with no audit trail

When more than one agent works in one worktree, `.git/index` is the piece of state
they all share and none can reconstruct the history of. On 2026-08-29 25 paths
turned up staged and **no session that could be asked had run `git add`**.

The tempting evidence is `.git/index`'s mtime. It does not work: the mtime records
the most recent index *write*, and every commit rewrites the index, so by the time
anyone looks it has usually been overwritten by unrelated activity. Reading a
staging time off it produced a confident and wrong answer here — including a
provenance claim in a commit message that the mtime could not actually support.

**Rule: never infer authorship or timing from the index.** Scope checkers to
`HEAD`, commit with explicit pathspecs so a stray staged entry cannot ride along,
and if it matters who changed something, ask the other sessions rather than
measuring the filesystem.

The same day produced a second version of this from the other direction: a broken
`&&` chain (above) fell through to a `git push`, which published **another
session's HEAD commit** — its own committed work, so nothing was corrupted, but
not the pushing session's to publish. Both incidents share a cause that has
nothing to do with the files anyone edited:

> **Working-tree files are per-session. The index and the refs are not.**

Edits can be partitioned by agreement — *you take `docs/`, I take `verify/`* — and
that works. Staging, committing, and pushing cannot be partitioned that way,
because they act on state every session shares and none can attribute afterwards.
So: pass explicit pathspecs to `commit` and not just to `add`; read `git log`
before and after; and treat `push` as publishing whatever HEAD is, not whatever
you just wrote.

### `result.txt` is written by the case script, not by the engine

Consequences, all of which actually happened: `bdbot.cli status` counts the run
as zero; completed work is not skipped on re-run; and a "clean up incomplete
runs" pass **deleted 6 completed runs**. 87 real runs are currently invisible in
`status` for this reason ([04 Cases](04-cases.md)).

### Editing code during a batch

11 runs died with `NameError`. `xargs` prints `done` even for a crashed child, so
`48/48 done` looked fine while only 37 `metrics.json` existed. **Count the
artifacts, not the exit lines.**

### Adding a spec field re-ids every run

One field addition changed the `run_id` of all 137 `chain-bend-2d-dlvo` runs at
once. **Aggregate by tag, never by `run_id`.** And the reverse hazard is worse:
a spec that omits the physical system entirely will keep the same `run_id`
across a 16× change in `τ_B` — see [03
§4](03-knowledge-base.md#4--provenance-and-tiers--a-number-without-a-source-is-not-a-number).

### A content hash over a `pow` result is not portable

★ `run_id` is the sha256 of the spec's JSON, which serialises floats at full
`repr` precision. That makes the identity of a run depend on **every bit** of
every derived float — including bits that IEEE-754 does not pin down.

`sqrt` is required to be correctly rounded, so it is bit-identical everywhere.
**`pow` is not.** Measured on `soft-r3`'s `Gamma = A / a_mean³`:

```
a_mean       = sqrt(pi/(4*0.35)) = 1.4979969134027407     identical on both
a_mean**3    = 3.361497213033026        Apple libm,  osx-arm64
a_mean*a*a   = 3.3614972130330254       glibc + exact mults, linux-64
Gamma(A=100) = 29.748648790272707   vs  29.74864879027271   -- 1 ULP
run_id       = ...__A100__30caa5c9e0 vs ...__A100__079a25f073
```

So **the same physical system has two identities**, and which one you get
depends on the machine. It was invisible for six weeks because every run was
computed on one laptop; CI run 35026488825 was the first time the suite ran the
derivation on `linux-64`, and it surfaced as one failing test asserting a
literal digest. 104 of 296 specs carry a hashed `params.Gamma`, and **all 104
move** under a 1-ULP shift (measured).

**What is *not* broken:** `LoadedSpec.verify_hash()` re-hashes the **stored**
content, so it is portable and passed on `linux-64`. Rule 2's hand-edit
detector still works. Only *re-deriving* a spec from the physics is
platform-dependent — which is what "reproduce this run elsewhere" means.

**And it cannot be fixed by normalising the payload.** Rounding `Gamma` to 12
significant figures makes the derivation portable and renames the archive —
including `runs/soft-r3-2d-A-sweep__A100__30caa5c9e0`, which the **sealed**
`campaigns/s30_preregistration/prediction.yaml` cites by name and whose sha256
is locked into 12 `SEALED.sha256` files. Renaming means editing a sealed
pre-registration after the fact, which is the one thing pre-registration exists
to prevent. The mutation is in
`tests/test_soft_r3_init.py::test_the_ulp_tolerance_still_refuses_what_it_exists_to_catch`
and it is CAUGHT, so this trade-off is measured rather than asserted.

**How to live with it.** A run_id assertion must be written against the payload,
not the digest string: `test_the_archived_run_id_does_not_move` now accepts any
digest reachable by a one-ULP shift of one float leaf of the archived payload,
and refuses everything else — which still catches the two mutations that
originally survived this file (`init`, `n_x`/`n_y` written on the default path)
plus a 10th-digit change in `phi`. If the archive is ever re-identified — which
can only happen at a campaign boundary, where re-sealing is legitimate —
normalise **per field**, at a precision each field's conditioning justifies —
not "the payload at 15 digits", which is what this section said for one commit.

⚠ **That number does not generalise, measured.** 15 significant figures unifies
`params.Gamma` (0 of 104 specs disagree) and does **not** unify
`params.k_bond_star`, hashed in 186 specs. `cases/chain_bend_dlvo_2d.py`'s
`find_well` takes a central second difference with `dh = h_min*1e-4`, so the
numerator loses about eight digits to cancellation and one ULP in a single
`U_star` evaluation becomes **2.96e-9 relative** in the hashed value:
`1042362.8817700658` against `1042362.8848514813`, which agree at 9 significant
figures and disagree at 12 and at 15. Same defect class as the one this section
is about — a number measured on one field, restated about "the archive".

**The general form:** a content address is only as portable as the least
portable operation upstream of it. Hash rounded values, or accept that the
address is machine-local.

### `dt` candidate lists that omit a stiffness

With `--kt-scale` near 200 the trap becomes the *fastest* mode, and `dt` had not
been recalculated. Whenever a knob can reorder the timescales, the `dt`
derivation has to be re-run — the candidate list is not static.

### `metrics` merge semantics

`MET.build(extra=...)` merges with `m.update(extra)`, so `finalize()`'s `result`
lands at top-level `metrics["result"]`, **not** `metrics["extra"]["result"]`.
Get it wrong and the value is silently `nan`.

### Returning full arrays from `sample()`

If `sample()` hands back the whole N-particle array every sample,
`observables.npz` goes 448 KB → **148 MB (330×)**. Accumulate derived quantities
in a closure; store raw only the subset you need.

### Defining "well escape" by `U = 0`

The DLVO outer branch **asymptotes** to `U→0⁻` and never crosses it, so solving
`U=0` returns `nan`. Bond rupture is decided by the **maximum tensile force**
(`F_max = 810.4 kT/d` at h=14.6 nm), not by an energy crossing.

### Korean labels in matplotlib

The default `DejaVu Sans` has no Hangul, so labels render as `□`. Fonts that do
have Hangul (`AppleGothic`, `Apple SD Gothic Neo`, `NanumGothic`) are missing
`−` (U+2212) and `ŷ` (U+0177) — measured. **Do not fix it by switching fonts.**
Write axes, legends, titles and annotations in English from the start; confirm
zero `missing from font` warnings.

---

## 3 · Traps in reasoning, not in code

These cost the most, because no test catches them.

**Comparing integrators at finite temperature has no power.** Testing whether
`τ_p/τ_fast = 0.60` (fastest mode not overdamped, ζ=0.65) contaminated the
observable: `OverdampedViscous` vs `Langevin(kT=0)` across all 7 ω differed by at
most **0.159 %**. The thermal comparison (Brownian vs Langevin, kT=1) **could not
have excluded a 47 % effect** — because `|ŷ|/ℓ_k < 1`. **Test integrator
assumptions with a `kT=0` deterministic difference**, where noise is zero and
transients cancel as a common mode.

**A metric must be tested against both extremes before it is trusted.** The
first shape metric built for `chain-bend-2d-dlvo` had **0.1σ** of discriminating
power on configurations that were obviously different by eye.

**A metric's discriminating power depends on the protocol.** Bow separated DLVO
from JKR at **22.3σ** under a soft trap and at **1.4×** under a stiff one. The
rule: **free deformation → measure shape; imposed deformation → measure force.**
The earlier "bow is the best discriminant" claim is true *only* in the soft-trap
regime, and stating it unconditionally would have been wrong.

**Improving the apparatus can make the measurement worse.** Raising `k_t` fixes
tracking, but past ×300 the bracket in `K′ = k_t(ŷ_c/ŷ − 1)` shrinks to 0.09 —
**better tracking, worse measurement.**

**"Too expensive to run" is usually a claim about implementation, not physics.**
`md.force.Custom` calls Python every step: **26× slower**. Swapping in a ghost
particle plus `bond.Harmonic(r0=0)` — exactly `½k r²`, on the compiled path —
turned a 25-day sweep into **1.16 days**. Measure which kind of cost it is first.

**A harmonic approximation is not the prediction.** Predicting bond-stretch
variance as `kT/k_bond` disagreed with measurement at 6.67σ. The DLVO secondary
minimum is asymmetric, and Boltzmann-integrating the basin gives **4.57×** the
harmonic value — that is the real prediction, and with it the discrepancy
converged to 2.15σ as sampling grew. Same family as soft-r3's Einstein-cage
approximation failing under anharmonicity.

**A checker you have not tried to break is not a checker.** "Silently passing"
and "not checking" are different states. `verify/verify_intake_guards.py` caught
a real crash bug precisely because it was written to break things.

---

## 4 · When you find a new one

1. Write a reproduction script into [`verify/`](../verify/).
2. Add it to skill `bd-hoomd`, and decide whether it earns a ★ — *does it fail
   without an error?*
3. File a KB entry with `origin: tooling` and a **cause, not a symptom**.

There are **57** `tooling` entries of 145 (measured 2026-09-15: 57 tooling · 50
method · 25 handbook · 10 intake · 3 paper). That number is the honest measure of
how much of this work is fighting the instruments rather than the physics — and
it had drifted from 48, which is the same class of defect as the counts in the
README.
