---
name: bd-intake
description: |
  The interpretation protocol for reading a physical system out of a sketch, a
  handwritten note, a whiteboard photo or a paper. Transcribe first, state the
  ambiguity, leave what is absent as null. Read this when creating or amending an
  observation.yaml from an image, when extracting parameters from a sketch, or
  when deciding what is missing. These rules came out of actually reading five
  sketches — the purpose is to suppress invention.
---

# The intake protocol (L0)

> The rules in this document came from **actually reading five sketches**. Three
> of the five were blocked by a single parameter, and **pinning that down to
> exactly one thing** was what intake achieved. The checking is mechanical:
> `$PY -m bdbot.cli intake check <folder>` (`bdbot/intake.py`).

## 0. The absolute order

```
① transcribe        write down exactly what is visible. Do not interpret.
② structure         split into entities, values, goals. Only what is in the sketch.
③ state ambiguity   record everything you are not certain of, with its impact and your lean.
④ state unread      record what you could not read. Writing "none" requires a basis.
⑤ state missing     values absent from the sketch that have to be supplied. **Do not invent them.**
```

**Do not skip ①.** Interpreting before transcribing makes you read things that
are not in the sketch. This actually happened: because of a folder name, the
master plan had a case recorded as "bonds + **bending stiffness** + driving", and
transcription showed **there was no bending in the sketch at all.** Transcribe
what is on the paper, not what is in the folder name, the expectation, or the
plan.

## 1. Transcription rules

- Mark the layout: `[left figure]`, `[right text]`, `[bottom]`
- Equations as they appear. Do not "tidy" `U(r) = ½ k_t r²` into `U = 0.5*k*r**2`
- Record struck-through and overwritten text too — it may be a cancelled
  definition (one sketch had `"Δr … = optical center"` struck out with a line)
- Mark an illegible character as `[illegible]` and add it to `unread_regions`
- **Count the particles in the figure** — but whether the count becomes `N`
  **depends on the count**:

  | Particles in the figure | Status of the drawn count |
  |---|---|
  | **fewer than 10** | **that is the real `N`.** Count exactly — "8 circles", not "several circles" |
  | **10 or more** | **treat it as schematic.** `N` is proposed separately, from the observable |

  The threshold of 10 is a **rule of thumb**, not a sharp value (confirmed with
  the user 2026-08-06). The basis: the scale at which a person can draw a system
  **one particle at a time** is roughly single digits, and past that they start
  drawing "lots, like this".

  ⚠️ **Counting is still worth doing above 10** — it just does not become `N`.
  The count and the **connectivity** become the **baseline for what the
  generation protocol has to build.** A real case: counting 21 circles showed
  the topology was a **tree (zero loops)**, and since the goal was `G′(ω)` and a
  tree has no closed path to carry stress, **the goal was undefined** — which is
  what justified the "schematic" call and led to the requirement that compression
  create loops. Without counting, none of that would have been visible.

  **The count goes in `stated_quantities` (a fact); the proposed `N` goes in
  `missing_required` as `kind: choice` (my judgment)** — do not mix them.

  The criterion is not only the count: also ask **whether the target observable
  is even defined at that `N`.** If it demands a **bulk or ensemble** quantity —
  `G′(ω)`, the stress tensor, an rdf — the figure is schematic even with few
  particles; conversely if the goal is that object's own response (one chain's
  `K′`, one trap's PSD), the drawn count is the real `N` however many there are.
  **When the count and the observable disagree, the observable wins.**

  When proposing `N`, the basis comes **from the observable**, not from a feel for
  "lots": statistical error `∝1/√N` against the target precision · how many
  correlation lengths must fit in the box · the minimum size at which percolation
  or loops exist · and **cost** (a sweep in this project went from 25 days to
  1.16 days). State the upper bound too.

## 2. What counts as an `ambiguity`

Something is ambiguous if **it can be read more than one way and the choice
changes the result.** Each entry gets:

| Key | What |
|---|---|
| `id` | `A1`, `B2` … case prefix plus number |
| `issue` | what is ambiguous (one sentence) |
| `impact` | **how this choice changes the result.** Write it as a number |
| `lean` | which way you lean, plus the basis. `null` if none |
| `resolution` | **leave it `null`.** That is the human's slot |

`impact` is what makes the entry worth anything. Writing only "this is ambiguous"
gives a person nothing to decide with. A good one: *"at d=1µm, ℓ_p/d=2.5 so
activity is visible in the MSD; at d=5µm it is 0.5 and nearly invisible."*

**Split `impact` per observable.** Not "affected / unaffected" but which
observable, by how much. That is what determines **what work can proceed while
the confirmation is pending.**

> `trap-2d-5um` A1: whether `R=5µm` is a radius or a diameter was ambiguous →
> `⟨x²⟩=kT/k` is independent of `d`, so **impact zero**, and only `τ_k` changes by
> exactly 2×. That is why the golden-test design was finished before the
> confirmation arrived.

⚠️ **Do not gain confidence from having listed several bases in `lean`.** In that
A1, three bases were given; one of them ("read as the trap region it would be
30,000 kT, which is meaningless") really did exclude one option, but the
remaining two (radius vs diameter) were not separated at all. **The stronger the
exclusion argument, the more confident you become even when the remaining options
are barely distinguished.** The answer was diameter, and `lean` was wrong.

⚠️ **Be suspicious of zero ambiguities.** That is rare in a hand-drawn sketch.
The checker warns.

### 2.1 An empty `stated_goals` is also a blocker ⭐️

**Not knowing what will be measured means you cannot set `T_obs`, the sample
interval, or the success criterion.**

| Case | Goal in the sketch |
|---|---|
| `abp-rod-2d-run-flip` | "measure MSD, MSAD" ✓ |
| `trap-2d-5um` | **absent** — settled by asking the user |

Only after deciding to look at the PSD did `T_obs = 2000 τ_k` (securing
`1/T ≪ f_c`) and the sample interval `τ_k/10` (`f_Nyq ≫ f_c`) get determined. Not
knowing the goal, they would have been chosen differently, and the PSD would have
turned out unobtainable after the whole thing had run.

If `stated_goals: []`, do not push it into `choice` — **ask the user.**

## 3. Missing values — separate `physical` from `choice` ⭐️

Each `missing_required` entry gets a `kind`:

| `kind` | What | Does it block L2? |
|---|---|---|
| `physical` (default) | a property of the system. A human must supply it, or it must be found in the KB | ❌ **it blocks** |
| `choice` | a simulation choice — box size, observation window, sample count | does not block |

**Without this distinction the verdict is wrong.** When the tool was first run
over five files, `trap-2d-5um` and `soft-r3` — both already completed end to end
— came back BLOCKED. The cause was that `L` (box) and `T_obs` (observation
window) sat in the same list as genuinely unknown physics.

An unstated `kind` is treated as `physical` (conservative).

**Splitting `physical` once more reduces the confirmation cost** — does the value
enter the result, or is it used only by a check?

| Example | Where it is used | Safe at tier 3? |
|---|---|---|
| `eta` (viscosity) | `γ → τ_k, D_t` — **directly sets the result** | ❌ needs confirmation |
| `rho_p` (density) | `τ_p` — **model-validity check only** | ✅ zero effect on a BD result |

`trap-2d-5um`'s `rho_p` was arbitrarily assumed as silica (tier 3), but BD does
not use mass, so the result does not change. Writing "…used only by a check, no
effect on the result" into `note` lets a person set their confirmation priority
immediately.

### When filling with an assumption

```yaml
- symbol: eta
  kind: physical
  what: "solvent viscosity"
  assumed_value: 0.851
  assumed_unit: mPa*s
  confidence: 1          # ★ mandatory whenever assumed_value is present; the checker rejects otherwise
  note: "water@300K handbook. Medium not stated in this sketch. Confirmed as water in case 1-A -> inherited"
  resolution: null
```

**An `assumed_value` without a `confidence` is a value with no provenance**
(violating rule 3 in CLAUDE.md). The checker raises it as an error.

tier: `0` directly given or handbook · `1` literature plus verification, or **an
inherited convention that was confirmed** · `2` literature, unverified · `3`
arbitrary assumption

⚠️ **Tier 1 by inheritance is the dangerous case, and it bit this project.**
`T = 300 K` is recorded as tier 1 in every case and is in fact a *choice*
inherited from a sketch with no temperature. Water's viscosity is 2.06 %/K
sensitive, so at 298 K `η` is off by −4 % and at 293 K by −14 %, and every
timescale follows. Inheriting is legitimate; recording it as if it were measured
is not.

### If you cannot fill it, leave `null` and write `what`

```yaml
- symbol: U_ij
  what: "chain bond potential — the elasticity in G'(omega) comes from here. Blank in the sketch (C1)"
  assumed_value: null
  resolution: null
```

Without `what` you cannot tell a person **what is needed and why**. The checker
warns.

## 4. Inheriting a value from another sketch in the same notebook

The five sketches were drawn consecutively in one notebook, and **four of them did
not state the medium or the temperature.** The user confirmed "water, 300 K" for
case 1-A, so that convention was carried forward — but:

- mark it `tier: 1` (not direct confirmation but **inherited convention**)
- state in `note`: "not stated in this sketch, inheriting the value confirmed in
  case N-A"
- remember the inheritance **can be wrong.** Applying the `R=5µm` convention to
  `abp-rod` disagrees with that sketch's own `τ_R=0.5s` by 160× (so that one case
  probably has different particles)

## 5. Do not simply trust the sketch's values — check self-consistency ⭐️

When the values you read are linked by a physical relation, that **is a
cross-check.**

A real case: reading `τ_R = 0.5 s` as a rotational diffusion time and inverting
`τ_r = πηd³/kT` gives `d = 0.918 µm`. Applying another sketch's `R = 5 µm`
convention makes the Stokes prediction `80.7 s` — off by **160×**. ⇒ the
conclusion is that this case has ~1 µm particles and that `τ_R` is the rotational
diffusion time. **Two ambiguities (D1 size, D2 the meaning of `τ_R`) resolved each
other.**

Record such an inversion in `ambiguity.evidence`. It is an inference, so
`resolution` still stays `null`.

## 6. Check dimensional consistency first

`soft-r3`'s `U_ij/kT = A/r³` — whether `A` is dimensionless or carries `µm³`
split the physics by 125× (at `d=5µm`). When you see a coefficient with no unit
in the sketch:

1. List **every** dimensionally consistent reading
2. Compute the dimensionless groups under each reading
3. Choose the reading that **does not contradict the sketch's stated goal** (and
   state that this is an inference)

In `soft-r3`, if `A` carried `µm³` the entire sweep would sit below `kT`, so
**nothing** the sketch asked for — rdf, Voronoi structure — would appear. Only the
dimensionless reading is consistent with the goal.

## 7. Do not

- Record something from a folder name, a plan or an expectation as if it were read
  from the sketch
- Skip transcription and start interpreting
- Fill a material property absent from the sketch with no basis (an
  `assumed_value` without a `confidence`)
- Fill `resolution` yourself (that is the slot for the first human confirmation)
- Defer an ambiguity as "we can decide later" — writing `impact` as a number is
  what makes the priority visible
- **Omit** the `ambiguities: []` / `unread_regions: []` keys entirely (the checker
  rejects that)

## 8. Procedure

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python

$PY -m bdbot.cli intake init  intake/<case>     # generate the template
#   -> read the image and fill in steps 1 through 5
$PY -m bdbot.cli intake check intake/<case>     # schema plus readiness verdict
#   -> FAIL:    a schema error. Fix it
#   -> BLOCKED: schema intact, unresolved physical gap -> **ask the user**
#   -> READY:   the L2 system file can be written
$PY -m bdbot.cli system check  intake/<case>    # L2 checks (tier, derived_from, recomputation)
```

`BLOCKED` is not a failure. It means **you pinned down exactly what is missing**,
and that is where intake's job ends. Do not invent something to force it to
READY.

### 8.1 Right after `READY` — draw the scale table and feed it back ⭐️

**There are things intake alone cannot see.** Only the scale table reveals them.

`trap-2d-5um` looked entirely clean at intake (all three values tier 0). Only
after drawing the scale table:

```
tau_p = 3.26 µs   tau_k = 4.01 ms  ★   tau_B = 242 s
```

`τ_B` is 242 seconds while the trap catches the particle in 4 ms, so **free
diffusion is never realized.** That changed the `dt` choice, the observation
window and the reporting units (had `τ_B` been used as the reference, one step
would have been 6× the relaxation time — skill `bd-physics`).

So there is **one feedback step** between the observation and the system file:

1. Draw the scale table
2. **See what the governing timescale actually is** — it may not be `τ_B`
3. If that changes any `choice` entry in the observation (`T_obs`, the sample
   interval, `L`), go back and fix it
