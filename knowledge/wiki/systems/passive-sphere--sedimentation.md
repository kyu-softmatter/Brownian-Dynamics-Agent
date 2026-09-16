---
type: system
author: agent
drafted: 2026-09-16
confirmed_by:
system: passive-sphere
dynamics: sedimentation
status: draft
cites:
  - knowledge/source/papers/2014-newman-yethiraj-sedimentation-equilibrium.md
  - verify/verify_sedimentation_wall.py
  - .claude/rules/overdamped-stability.md
---

# `passive-sphere` × `sedimentation`

> **Why this card exists:** the non-dimensionalization and the key parameters are
> **not fixed by the system alone.** The same colloidal sphere needs a different
> reference time depending on whether equilibrium structure or transport is the
> target, and the `Δt` gate, the observables and the benchmarks all follow from
> that choice. **S3 NONDIM reads this card to decide how to non-dimensionalize.**

⚠️ **This is the first card written in English.** The other ten are Korean, and
the divergence is deliberate rather than accidental: `bdbot/`, `cases/` and
`tools/` are already English, the repository is public, and the section numbering
and frontmatter keys here are identical to `_TEMPLATE.md` so nothing that parses
this file can tell the difference.

---

## 1. The system — what is in the box

| item | value |
|---|---|
| particle kind | one passive sphere species, sterically stabilised |
| dimensionality | **3D** |
| interaction | WCA (an excluded-volume stand-in for a hard sphere) |
| polydispersity | monodisperse (the source reports < 5 %) |
| **boundary** | **periodic in x, y · a repulsive wall at the bottom in z** ★ |
| solvent | implicit, free-draining (no hydrodynamic interactions) |
| external field | **uniform gravity**, `F_z = −kT/l_g`, via `md.force.Constant` |

★ The boundary is what makes this pair new. Every other card here is fully
periodic; this one has a wall, and `hoomd.wall` appeared nowhere in this
repository before 2026-09-16.

## 2. The target dynamics — what is being looked at

| | |
|---|---|
| question | does the equilibrium density profile `Φ(z)` reproduce the barometric law in the dilute limit, and does its local equation of state follow Carnahan–Starling? |
| equilibrium / non-equilibrium | **equilibrium** — sedimentation equilibrium is a true Boltzmann state in an external potential |
| is there a steady state? | yes, and it is the equilibrium state; there is no flux at steady state |

## 3. Reference units ★ — valid for this pair only

| dimension | choice | why this one |
|---|---|---|
| length | `d = 2a` (the **diameter**) | excluded volume sets both the structure and the EOS, and `Φ` is defined through `d³`. The gravitational length `l_g` is the *other* length and becomes a dimensionless group, not the unit — it varies by orders of magnitude across the source's three solvent mixtures while `d` does not |
| energy | `kT` | equilibrium in an external potential; `l_g` is itself defined by `kT` |
| time | `τ_d = d²/D₀`, `D₀ = kT/(3πηd)` | same as `passive-sphere--equilibrium-structure`, so the two cards share a clock |

So the scale rule is **`brownian`**, registered in
`simbot/nondim.py::CARD_SCALE_RULES`.

**Derived quantities and the identity check**

```
gamma  = 3 pi eta d                       [kg/s]
D_0    = kT / gamma                       [m^2/s]
tau_d  = d^2 / D_0                        [s]
l_g    = kT / ( (4/3) pi a^3 drho g )     [m]        <- note a, not d
v_sed  = kT / (gamma l_g) = D_0 / l_g     [m/s]      <- the drift speed
tau_settle = L_z / v_sed = L_z l_g / D_0  [s]        <- ★ NOT l_g^2/D_0

identity: v_sed * l_g / D_0 == 1          (exact, by construction)
identity: tau_settle / tau_d == (L_z/d)(l_g/d)
```

## 4. Dimensionless ledger ★

| symbol | definition | meaning | source range |
|---|---|---|---|
| `Φ` | `N π d³ / (6 V)` | volume fraction | 0.017 (mean) |
| `l_g*` | `l_g / d` | gravitational length in diameters | **13.375** (0.8 µm, 70:30) · 10.9 (1.0 µm, 60:40) · ≫100 (50:50) |
| `F_g*` | `1 / l_g*` | the gravity force, in `kT/d` | 0.0748 |
| `L_z*` | `L_z / d` | box height | must satisfy `exp(−L_z*/l_g*)` ≪ the smallest `Φ` resolved |
| `ε_w*`, `σ_w*` | wall depth and range, in `kT` and `d` | bottom wall | this repo: 20, 0.3 |

> **Convention traps — where this goes quietly wrong:**
>
> ① **`l_g` is defined with the radius `a`, our length unit is the diameter `d`.**
> `l_g/a` and `l_g/d` differ by exactly 2. The source itself carries both
> conventions: its text writes `Pe_g = 2a/l_g` while its Table 1 column header
> reads `Pe = a/l_g`, and the tabulated values match the second. **Never quote a
> `Pe_g` without saying which.** This card uses `l_g* = l_g/d` and nothing else.
>
> ② `Φ` uses `d³`. A `Φ` computed from `a³` is 8× wrong and still looks like a
> volume fraction.
>
> ③ `Φ(z)` is a **local** volume fraction in a slab, so it depends on the slab
> thickness. State it. The EOS comparison is only as good as that binning.

## 5. Key parameters — measured values from the literature

| parameter | value | source |
|---|---|---|
| `2a` | 0.8 µm (also 1.0 µm) | [Newman & Yethiraj 2014](../../source/papers/2014-newman-yethiraj-sedimentation-equilibrium.md) |
| `l_g` | 10.7 µm | same, Table 1 |
| `D₀` | 0.44 µm²/s | same |
| `Φ` | 0.017 | same |
| `ρ_PMMA` | 1259 kg/m³ | same, Table 1 |
| `T` | **296 ± 2 K, derived** — the paper does not state it; `l_g`, `a` and `Δρ` close on it to 1 % | same, §3 of the distillation |
| `η` | 1.22 mPa·s, derived from `D₀` by Stokes–Einstein | same |

`reproduced: no` on that source, so these are **records of fact, not evidence**
(`.claude/rules/verify-against-literature.md`). Mark them `[source, not
reproduced]` in any report until a run of ours reproduces them.

## 6. Observables

| observable | how the error is obtained | note |
|---|---|---|
| `Φ(z)` | Poisson counting per slab, `√n/n`, over blocks | the primary |
| `l_g` fitted | weighted least squares on `ln n(z)`; the **slope** is fitted, so the origin and the wall layer cannot shift it | the `implementation_check` |
| `Z(Φ) = βP/ρ` local EOS | from the profile by the hydrostatic relation `dP/dz = −ρ kT/l_g` | the `hypothesis` vs Carnahan–Starling |
| horizontal MSD | `bdbot.stats.block_sem` | `D_msd = kT/γ` holds in x, y |
| vertical MSD | — | ⚠ **saturates at ~`l_g²`**; do not fit `D` to it |

## 7. Applicable gates — which of the general ones fit this pair

| gate | applies | why |
|---|---|---|
| equilibration (`pymbar`) | ✅ | a true equilibrium state |
| `D_msd = kT/γ` | ⚠ **x, y only** | the z direction is confined by gravity, so its MSD plateaus |
| advective displacement `v Δt / d` | ✅ **required** | gravity is a constant drift, `v = D₀/l_g` |
| `r_c < L/2` minimum image | ✅ | as always |
| ★ **the wall binds `dt`** | ✅ **required** | measured: a WCA wall at `σ_w = 0.5 d`, `ε_w = 1 kT` reaches `F = 1517 kT/d` at `h = 0.4 d`, so `dt_max_force = 2e-5`. Use a **bounded** wall (`Gaussian`, `F_max = ε_w/(σ_w√e)`) and the escape mode closes by construction |
| ★ **box height** | ✅ **required** | `exp(−L_z*/l_g*)` is the weight at the periodic lid. At `L_z = 15 l_g` it is 3e-7; at `6 l_g` it is 2.5e-3 and contaminates the third decade of the profile |
| ★ **settling time** | ✅ **required** | equilibration is `L_z l_g / D₀`, **not** `l_g²/D₀`. Measured: a uniform start run for `8 l_g²/D₀` in a box of height `15 l_g` fitted `l_g = 11.5 d` against an imposed 4.0 — **+187 %**, purely from stopping before the particles had fallen |

⚠ The last three did not exist before this card. They came out of
`verify/verify_sedimentation_wall.py` failing twice, and each one is a number
that was measured rather than reasoned.

## 8. Benchmarks

| # | layer | claim | threshold |
|---|---|---|---|
| B1 | ① analytic | the dilute profile is `n(z) ∝ exp(−z/l_g)`, **exact** for point particles | the fitted `l_g` within 4σ of the imposed one; measured 13.179 ± 0.127 against 13.375 (−1.46 %, 1.54σ) with no pair potential |
| B2 | ② self-consistency | `v_sed l_g / D₀ = 1`; the fitted `l_g` is independent of `dt` | < 3 % between two `dt` |
| B3 | ③ literature | the local EOS follows Carnahan–Starling — *"excellent agreement"* per the source, **with no percentage given**, so our own threshold has to be chosen and declared | to be set in the prediction document |

★ B3's role is `hypothesis`, not `implementation_check`: the simulation imposes a
**WCA core, not a hard sphere**, and imposes no equation of state at all. A
mismatch is a result about how well a soft core stands in for a hard one at this
`Φ`, and [[wca-reproduces-carnahan-starling]] is this repository's own prior
measurement of exactly that substitution.

⚠ **What this pair cannot do.** The source's headline dynamical result,
`D = D₀(1 + K_L Φ)` with `K_L = −2.80`, is a **hydrodynamic** correction.
Free-draining BD with a scalar `γ` has no hydrodynamic interactions, so it cannot
produce `K_L` — stated here so that it is not adopted as a target and then
"discovered" to fail. That is rule 7′ applied before the run rather than after.

## See also

[[wca-reproduces-carnahan-starling]] ·
`knowledge/wiki/systems/passive-sphere--equilibrium-structure.md` ·
[overdamped-stability](../../../.claude/rules/overdamped-stability.md)
