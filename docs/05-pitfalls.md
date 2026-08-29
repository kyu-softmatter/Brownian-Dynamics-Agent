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

There are 44 `tooling` entries. That number is the honest measure of how much of
this work is fighting the instruments rather than the physics.
