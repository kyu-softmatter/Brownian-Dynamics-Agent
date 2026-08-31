# S1 — reading: the **delta** against the drawing

Source [`inputs/soft-r3-2d-A-sweep/sketch_01.jpeg`](../../inputs/soft-r3-2d-A-sweep/sketch_01.jpeg)
· card [`soft-repulsive-2d--equilibrium-structure`](../../knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md)

> The reading of the drawing itself is already in card §1. **This document records
> only where the user's directives override the drawing.** An overridden item is
> `provenance: user`.

## 1. Taken straight from the drawing

| item | value | provenance |
|---|---|---|
| interaction | `U/kT = A/r³` (no hard core) | from_drawing |
| dimension | 2D | from_drawing |
| `N` | `100` | from_drawing |
| boundary | periodic | from_drawing |
| initial condition | irregular scatter → `random` | inference |
| question | "final configuration?" → `rdf`, `voronoi plot`, structure analysis | from_drawing |

## 2. What the user's directives overrode ★

| # | item | drawing | **this run** | provenance |
|---|---|---|---|---|
| D1 | `A` list | `0.1, 1, 10, 100` | **`0.1, 1, 10`** (`100` excluded) | user |
| D2 | box shape | `L_x = L_y` | **fixed square** (`hex_commensurate` forbidden) | user |
| D3 | box size | absent | **`σ = 5 µm` disc coverage `< 10 %`** | user |
| D4 | observation | the final arrangement | **track `g(r)`, Voronoi and the defect structure over time** | user |

### D2 is not incidental, it is required

Card §3 · [[box-shape-confounds-initial-condition-comparison]] — changing the box
shape changes `r_cut` (derived from `= min(L)/2`) and the allowed `k` grid together,
and the `S(k)` 6-fold modulation came out **28x** apart. Fixing it square removes
that confounder. Which is why **only `init = random` is run** -- a `hex` initial
placement in a square box has the periodic boundary cut the lattice and create
artificial defects (`build_soft2d` raises).

### D3 — what the coverage does and does not set ★★

`σ` **does not enter the dynamics.** There is no hard core and `n* = 1` holds by
definition, so **the dimensionless physics is fully determined by `A` alone.** There
are exactly two things `σ` sets:

1. **The time scale** — `τ_d = d²/D₀`, `D₀ = kT/(3πησ)`. Without `σ` there is no way
   to attach seconds.
2. **The validity of the point-particle idealization** — how crowded 100 discs of
   5 µm actually are.

⇒ **Changing the coverage does not move `ψ₆`, `g(r)` or the defects by a single
digit.** What changes is the axis labels and `τ_d` in seconds. The sensitivity is
`|S| = 0`, so under CLAUDE.md §the question rule it was decided without asking.

| value | basis | confidence |
|---|---|---|
| `σ = 5.0 µm` | user directive | high |
| **`d/σ = 3`** → `d = 15.0 µm` | `φ < 10 %` requires `d/σ > 2.8025`. `3` is the cleanest value above it (`simbot.build.box_si_for_coverage`) | high |
| **`L = 150.0 µm`** | `L = √N · d`, `N = 100` | high (derived) |
| **`φ = 8.7266 %`** | `φ = (π/4)(σ/d)² = π/36` | high (derived) |
| `T = 298.15 K` · `η = 0.890 mPa·s` | `knowledge/wiki/concepts/water-298k.md`, IAPWS 298.15 K | high |
| **`τ_d = 2292.4 s = 38.21 min`** | `d²/D₀`, `D₀ = 0.09815 µm²/s` (`simbot.units.scales_soft2d`) | high (derived) |

**Where to change it** — `SIGMA_SI` and `D_OVER_SIGMA` in
`scripts/soft2d_time_series.py`.

### ⚠️ A side effect of D3 — the reference discs overlap at `A = 0.1`

Converting the measured minimum separation of the 40 prior runs
(`runs/2026-07-29_soft-r3-2d-A-sweep`) into `σ`:

| `A` | min separation | in `σ` | overlap? | to remove the overlap |
|---|---|---|---|---|
| **0.1** | `0.2172 d` | **`0.652 σ`** | **overlaps** | `φ < 3.71 %` |
| 1 | `0.4586 d` | `1.376 σ` | none | `φ < 16.5 %` |
| 10 | `0.7164 d` | `2.149 σ` | none | `φ < 40.3 %` |

`A/r³` has no hard core, so **the model permits this overlap** — that is faithful to
the drawing. But at `A = 0.1` the picture of "5 µm discs" does not hold physically:
`βU(0.652 σ) = 0.1/0.652³ = 0.36 kT`, a barrier thermal fluctuation crosses easily.
⇒ **The `A = 0.1` result may only be read as that of a point-particle
soft-repulsive system.** This is sealed into S2 in falsifiable form (P6).

### D4 — what time resolution requires

The 40 prior runs **threw the equilibration window away** and analysed only the
second half of production. That answers "what did it become" but **does not answer
"when did it become that."**

⇒ This run sets `equil_tau = 0` and **samples the whole thing from `t = 0`.**
A total of `80 τ_d` (the same total as the prior run) · `400` frames
(`stride = 0.2 τ_d`). The `60–80 τ_d` window corresponds to the prior run's "second
half of production", so it is **a direct cross-check.**

## 3. Gates — card §7 as written

| gate | on/off | in this run |
|---|---|---|
| `force_displacement` (`max\|F\|Δt ≤ 0.005`) | ✅ **dominant** | `dt` is set from a **measurement** on the initial placement, per `A` |
| `thermal_displacement` (`√(2Δt) ≤ 0.03`) | ✅ | the reference length is `d` |
| `r_cut + buffer ≤ L/2` | ✅ | `L* = 10` → `r_cut = 4.80`, buffer `0.1` |
| `min_separation > r_min` | ✅ | `r_min = 0.1` |
| `packing_fraction` | ❌ off | undefined with no hard core (the reference-disc coverage is **a record, not a gate**) |
| equilibration verdict | ⚠️ **report only** | there is no threshold, so the drift is measured and the verdict is only proposed |
| initial-condition independence | ❌ off in this run | D2 means `random` only. The prior run already showed agreement to within `0.8σ` for `A ≤ 10` |

## 4. Remaining ambiguities

- `A = 10` is `−0.3 %` from Zahn's transition point `Γ = 59.88` (`Γ = 55.68`). Being
  **on a regime boundary**, it is expected to be the slowest under time resolution --
  that is the content of P3.
- `N = 100` does not satisfy `A = 10`'s `r_cut` requirement (`N ≥ 252`, card §9).
  `βU(r_cut) = 0.09 kT`. **An `N`-convergence check is outside this run's scope** --
  it is left open.
