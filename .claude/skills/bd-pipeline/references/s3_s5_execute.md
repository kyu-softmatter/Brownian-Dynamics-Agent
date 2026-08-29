# S3 · S4 · S5 — specify · non-dimensionalize · run

> These three stages are **almost entirely code**. Your job is to write
> `spec.yaml` and to read the gate results. `cli.py run` drives S3→S8 in one go.

## S3 — `spec.yaml`

### The fastest route: copy the example

```bash
cp examples/trap-2d-5um/spec.yaml <working dir>/spec.yaml
```

[`examples/trap-2d-5um/spec.yaml`](../../../../examples/trap-2d-5um/spec.yaml) is
a verified specification with 18 provenance fields filled in.

### Every value gets a `provenance` and a `basis`

```yaml
eta_si:
  value: 8.5566e-4
  unit: Pa*s
  provenance: assumed        # the number is accurate; "it is water" is the assumption
  basis: IAPWS table interpolation (not extrapolation). knowledge/wiki/concepts/water-298k.md
  confidence: medium
  affects: [tau_trap, D0]    # no effect on <x^2>
```

| `provenance` | When | Who may write it |
|---|---|---|
| `observation` / `from_drawing` | read directly off the source | anyone (including Haiku) |
| `derived` | computed from other fields | anyone |
| `rule` | derived from policy | anyone |
| `from_knowledge` / `from_paper` | from the wiki or a distillation | anyone |
| **`inference`** | derived from the source plus physics knowledge | **Opus only** |
| **`assumed`** | absent from the source, so supplied | **Opus only** |
| `user` | a human specified it for this run | marked automatically by code |
| `measured` | an experimental value | — |

`inference` and `assumed` **require** a `confidence`. Without one, `validate`
catches it. The authority boundary is stated in
[`.claude/README.md`](../../../README.md#authority-boundary).

⚠️ **`assumed` values inherit a false authority downstream.** `T = 300 K` is
recorded as tier 1 across every case here and is actually a *choice* inherited
from a sketch with no temperature — worth −4.4 % to −15.0 % on every timescale,
because water's viscosity is 2.06 %/K sensitive. Mark the tier honestly; the
tier field is the only thing that makes such a value findable later.

### ⚠ The YAML 1.1 trap — exponent notation

```yaml
value: 5e-3      # ❌ this is a STRING! No decimal point in the mantissa
value: 5.0e-3    # ✅ float
value: 6.3e6     # ❌ string. No sign on the exponent
value: 6.3e+6    # ✅ float
```

**This happened twice** (a throughput constant, and `session set`). Full account:
[`yaml-scientific-notation-parsed-as-string`](../../../../knowledge/wiki/findings/yaml-scientific-notation-parsed-as-string.md)

Writing via `simbot.spec.dump_yaml` is safe (a Python `float` always emits a
decimal point). **Hand-written YAML is the dangerous case.**

### Gate declarations — the card turns them on and off

```yaml
gates:
  equipartition: {status: required, reason: the first-class gate of this card}
  step_displacement_vs_sigma:
    status: off
    reason: no pair interaction, so overlap cannot occur    # ★ `off` requires a reason
```

- Only **registered** gate names (`simbot.spec.KNOWN_GATES`). A typo becomes
  **a check that never runs once** → `validate` rejects it
- Which gates turn on is decided by the card
- `required` is **a declaration, not a result**. It shows in the S3 report as
  "awaiting S7"

### Do not store derived values

Do not write derived values like `kT_si` or `tau_trap_si` into `spec.yaml` —
`simbot.spec.derive()` computes them. If they are written, `validate` recomputes
and compares to catch a mismatch, which is the only way a hand-edited derived
value gets caught.

### Check

```bash
<PY> -c "
from simbot.spec import SystemSpec, validate
r = validate(SystemSpec.load('<spec.yaml>'))
print(r.table())
print()
print('convention violations:', r.problems or 'none')
print('awaiting S7:', [c.name for c in r.deferred()])"
```

---

## S4 — non-dimensionalization

**All automatic.** Your job is to read the result.

```bash
<PY> -c "
from simbot.spec import SystemSpec
from simbot.nondim import reduce_spec, roundtrip_errors, nondim_table
sp = SystemSpec.load('<spec.yaml>')
r = reduce_spec(sp)
print(nondim_table(sp, r))
print('roundtrip error:', max(roundtrip_errors(sp, r).values()))
print('dt* =', r.dt_star, 'dominant constraint:', r.dt_dominant)"
```

On the `bdbot` side the same stage produces the L3 contract, and the spec checks
itself:

```bash
<PY> -m bdbot.cli nondim spec <case>       # -> specs/<run_id>.json
<PY> -m bdbot.cli nondim show <run_id>     # reproduce the report from the spec alone
```

★ **Self-sufficiency is decided by `nondim show`.** If the whole report draws
from the spec alone, the spec is self-sufficient — and the health layer never
imports case code, so this matters.

### Three things to check

1. **Where the scales came from** — the card must appear, e.g.
   `scales_harmonic_trap: (l_trap, kT, tau_trap)`
2. **Round-trip error `< 1e-12`** — larger is not an arithmetic slip but a
   **convention violation** (e.g. dividing by `τ_D` and inverting with `τ_trap`)
3. **The dominant constraint** — which `dt` constraint won. If it reads as an
   explicit value, a human set it

### `dt` constraints turn on differently per system

| Constraint | Turns on when |
|---|---|
| thermal displacement `√(2D₀Δt) ≤ 0.03σ` | **only with a pair interaction** |
| force displacement | pair interaction plus a **measured** `max|F|` (never estimated) |
| relaxation time `Δt ≤ 0.01 τ` | confinement or activity present |
| active displacement | active driving present |
| accuracy target | harmonic trap with a stated target bias |

★ **The displacement gate is not universal.** In a trap system the displacement
bound is **1086×** looser than the relaxation-time bound — the gate stops
nothing. Full account:
[`displacement-gate-is-1000x-loose-for-traps`](../../../../knowledge/wiki/findings/displacement-gate-is-1000x-loose-for-traps.md)

`dt/τ_D` is **recorded only**. Used as a gate it rejects runs that reached
publication.

⚠️ **And the candidate list is not static.** `choose_dt`'s displacement gate keys
off `bool(spec.pair)`, and the spec has no bond/angle field — so **a bond-only
system silently turns the gate off**, and the measurement that needed it
(`max|F*| = 1037.7`, `dt` cut 100×) came from exactly such a system. Whenever a
knob can reorder the timescales, re-derive `dt`.

### No card means an exception

```
KeyError: no scale rule registered for card 'colloid--new-thing'.
Improvised non-dimensionalization is forbidden — make a draft card from
_TEMPLATE.md first and register it in CARD_SCALE_RULES
```

**Do not route around it; make the card.** A result produced without a card
cannot be reproduced later.

---

## S5 — run

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>       # S1->S8 engine
<PY> -m bdbot.cli run <case>                                     # the 8-case engine
```

⚠️ **Read `bd-hoomd` before writing any new HOOMD.** 20 traps, several of which
produce wrong results with no error at all — including one where the force is up
to 96 % wrong while the energy stays exact, so no energy check finds it.

### Conditions under which the CLI stops **before** running

| Condition | Message |
|---|---|
| an S3 gate failed | `S3 gates not passed` |
| round-trip error `≥ 1e-12` | `S4 roundtrip gate violated — the scale convention is inconsistent` |
| seeds `< 4` | `a production run without error bars is forbidden` |
| budget overrun expected | `budget overrun expected — reporting instead of running` |

**Use `--force` only when the user explicitly asks.** If a gate was bypassed,
record that fact in the report.

On the `bdbot` side the pre-run gate blocks on three things only — hash mismatch,
a hard `FAIL`, and an L3 integrity error — and **shows warnings and thin margins
without blocking**. That is deliberate: it once rejected 80 of 83 specs with zero
real failures among them, and nobody noticed because the runner never called it.
An unwired checker cannot be wrong out loud.

### The tier ladder

For a card being run for the first time, the CLI warns:

> ⚠️ this is this card's first run. Policy does not allow skipping the ladder
> `[smoke, pilot, explore]` — this run becomes its first rung. **Do not cite the
> result as production**

### Error bars are free

At `k ≤ 4` the efficiency is 93 %, so **four seeds cost about the same as one**.
Policy enforces a minimum of four. **Four short runs beat one long run.**

⚠️ And four is a floor, not a target. In a system that generates stochastic
defects, a single run's block SEM underestimated the true spread by **1.09–2.28×**
depending on the observable and the velocity — and two published conclusions were
reversed by re-running with 9 seeds. If the system can produce discrete events,
size the ensemble to the events, not to the policy minimum.

### Failed runs are not discarded

If part of a batch dies it is recorded in the run manifest's `batch.failed` and
the CLI reports it. If it drops out silently, an error bar labelled "four seeds"
is actually three.

⚠️ Related: **editing code during a batch** killed 11 runs with `NameError`, and
`xargs` prints `done` even for a crashed child — `48/48 done` looked fine while
only 37 `metrics.json` existed. **Count the artifacts, not the exit lines.**

### Shaking a parameter — without running

```bash
<PY> -m simbot.session new <spec.yaml>
<PY> -m simbot.session set numerics.dt_star=2.5e-3 species.0.n_simulated=4000
<PY> -m simbot.session show
```

`set` **only estimates cost**. Running is `cli.py run`. A changed value is marked
`provenance: user`, and **the previous value and its original basis stay in the
basis field.**

### Convergence check

```bash
<PY> cli.py converge <spec.yaml>
```

Runs `dt/2`, `dt×2`, `N×2` and a seed shift, and judges **against the statistical
error**. Within `3σ` means "not distinguishable" — **not a proof of equality, but
a statement that this error bar cannot see the difference.**
