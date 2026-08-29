# overdamped-stability — `dt` is set by force, not by a timescale ratio

The position update in overdamped BD is `Δx = (F/γ)Δt + √(2DΔt)·ξ`. There is no
inertia, so **a particle under a large force simply travels that distance in one
step** — nothing carries it back. So `Δt` is chosen from **displacement per
step**, not as "some fraction of `τ_D`".

```
rms thermal displacement   √(2·dim·D·Δt) / ℓ  ≤ 0.03
force displacement         max|F|·Δt/γ  / ℓ   ≤ 0.03   (when excluded volume
                                                        or strong repulsion is present)
```

`ℓ` is **that system's reference length**, and it differs per system. For free BD
it is `σ`; for a harmonic trap it is `ℓ_trap = √(k_BT/k_t)`; for a repulsive pair
system it is the mean spacing `d`.

Do not use `Δt/τ_D` as the criterion: displacement goes as the **square root** of
`Δt`, so when the dimensionality changes the same `Δt/τ_D` gives a different
displacement.

**Why (the triggering incident):** this actually happened in a predecessor
project — WCA's `r⁻¹³` core **threw a particle out of the box** under overdamped
dynamics. The `Δt` chosen from thermal displacement was statistically fine, but
in one step where two particles happened to come close, `F` exploded, and with
no inertia the particle simply left. The symptom was not NaN but a **quiet box
escape**, so the log looked normal.

The second instance was a trap system (2026-07-28). `τ_D/τ_trap = 2.41e5` — five
orders of magnitude. Carrying over the free-BD `dt` rule based on `τ_D` gives
`Δt = 54 τ_trap`, which, divergence or not, **cannot see the trap at all.**
Not blowing up is not the same as being right.

**How to apply:**
- When adding a new system, **fix the reference length and time first.** The
  combination that makes `kT = D = γ = 1` in reduced units is that system's
  intrinsic length and the time to diffuse across it
- Read which constraint actually bound `dt`. That record is the only thing that
  answers "why this `dt`" after the fact
- **Also record what the rule you did not adopt would have given.** It is a cheap
  way to make a silent error visible
- When adding strong repulsion, check the minimum separation of the initial
  placement. Start overlapped and the first step throws the particle out
- ⚠️ The `dt` candidate list is not static. When a knob can reorder the
  timescales, re-derive it — with `--kt-scale` near 200 the trap becomes the
  *fastest* mode and `dt` had not been recomputed

**Anti-patterns explicitly forbidden:**
- **Measuring every system by `Δt/τ_D` alone** — change the dimensionality or the
  system kind and the same ratio gives a different displacement
- **"It didn't blow up, so it's fine"** — a `Δt` larger than the relaxation time
  cannot see the system even when it does not diverge
- **Watching only for NaN** — the symptom of this failure is a box escape, not NaN

See also: [axioms](axioms.md) ·
[docs/05-pitfalls.md](../../docs/05-pitfalls.md)
