---
type: source
kind: paper
lab_authored: false
title: Clusters in sedimentation equilibrium for an experimental hard-sphere-plus-dipolar Brownian colloidal system
authors:
  - "Newman HD"
  - "Yethiraj A"
year: 2014
journal: "arXiv preprint (no journal reference shown on the arXiv record)"
doi: 10.48550/arXiv.1412.3190
source_url: "https://arxiv.org/abs/1412.3190"
si_available: false
access: open
engine: "experiment (confocal microscopy, z-stacks)"
reproduced: no
parameters_extracted: yes
affiliation: Memorial University of Newfoundland (Yethiraj group)
ingested_at: 2026-09-16
ingested_by: agent
tags:
  - "sedimentation"
  - "gravity"
  - "hard-sphere"
  - "equation-of-state"
  - "carnahan-starling"
  - "dipolar"
  - "self-diffusion"
---
# Clusters in sedimentation equilibrium (Newman & Yethiraj 2014)

`arXiv:1412.3190v3 [cond-mat.soft]`, submitted 2014-12-10. Read from the arXiv
HTML; **the PDF was not obtained**, so there is no `raw_file` and figures were not
digitized — every number below is from the text or a table the text renders.

## Why this is in the wiki

★ **The first source in this repository for a system driven by gravity.** Every
existing case drives with a trap, a moving trap, an active force or a bond;
`grep -rl "hoomd.wall"` over the repository returned **nothing** before
2026-09-16 and `md.external.field` has no gravity class, so sedimentation was not
a buildable system here at all
([`verify/verify_sedimentation_wall.py`](../../../verify/verify_sedimentation_wall.py)
closes that).

What makes it worth building against, rather than just reading:

| evidence layer | what this paper supplies |
|---|---|
| ① analytic limit | the dilute profile is the barometric law `n(z) ∝ exp(−z/l_g)`, **exact** for point particles, and `l_g` is measured |
| ② self-consistency | `l_g`, `a` and `Δρ` are given independently and **close on each other** (§3) |
| ③ literature benchmark | *"The zero field equation of state (EOS) is hard sphere without any re-scaling of particle size"*, compared to **Carnahan–Starling** — for which this repository already holds a benchmark ([[wca-reproduces-carnahan-starling]]) |

Three layers, three different kinds, which is what `A1` asks for.

## 1 · The system

Sterically stabilised PMMA spheres in a decalin/tetrachloroethylene (TCE)
mixture, imaged by confocal microscopy in sedimentation equilibrium, with an
optional 1 MHz AC field that induces dipoles. **Our target is the zero-field
case**, where the paper states the system is hard-sphere.

| | this paper | our nearest existing case |
|---|---|---|
| dimensionality | 3D | `soft-r3-2d`, `trap-drag-2d` are 2D; `network_3d` is 3D |
| drive | **gravity + a bottom wall** | trap / moving trap / active force |
| pair interaction | hard sphere (zero field) | WCA, `A/r³`, DLVO |
| observable | density profile `Φ(z)`, EOS, long-time self-diffusion | `g(r)`, `ψ₆`, MSD, `τ_ac` |

## 2 · Parameters, as stated

| quantity | value | where |
|---|---|---|
| particle | PMMA, sterically stabilised | text |
| diameter `2a` | **0.8 µm** and 1.0 µm | text, Table 1 |
| polydispersity | *"less than 5 %"* in size | text |
| `ρ_PMMA` | **1259 kg/m³** | Table 1 |
| solvent | cis-trans-decalin + TCE, 50:50 / 60:40 / 70:30 by volume | text |
| `ρ_solvent` | 1113.8 … 1259 kg/m³ across the mixtures | Table 1 |
| `l_g` | **10.7 µm** (0.8 µm particles, 70:30) · 10.9 µm (1.0 µm, 60:40) · ≫100 µm (50:50) | Table 1 |
| `Φ` (mean) | **0.017** | text |
| `D₀` (bulk) | **0.44 µm²/s** | text |
| `D(Φ)` | `D = D₀(1 + K_L Φ)`, `K_L = −2.80 ± 0.05` | text |
| dipolar strength | `Λ = πε₀ε_f β² a³ E₀² / (2 k_B T)`, `β = (ε_p−ε_f)/(ε_p+2ε_f)` | text |
| `Λ` range | 0.0127 … 2.626 over `E` = 166.7 … 1666.7 V/mm | Table 2 |
| cell | 100 µm sample thickness; z-stacks of 77 × 77 × 108 µm at 29 fps | text |

`l_g` is defined in the text as the height at which the buoyant potential energy
equals `k_B T`:

    l_g = k_B T / ( (4/3) π a³ Δρ g )

## 3 · ★ The temperature is not stated — and the paper's own numbers pin it

**The paper gives no temperature and no solvent viscosity.** That would normally
block a dimensional intake (rule 1, rule 3). It does not have to, because the
three quantities it *does* give are not independent:

```
l_g = 10.7 µm · a = 0.4 µm · Δρ = 1259 − 1113.8 = 145.2 kg/m³
  ->  buoyant weight (4/3)πa³Δρg = 3.8173e-16 N
  ->  T = l_g × weight / k_B = 295.84 K = 22.69 °C
```

A room temperature, to about 1 %. The reverse check brackets it: 20 °C requires
`Δρ = 143.88` (`ρ_f = 1115.1`) and 25 °C requires `Δρ = 146.33` (`ρ_f = 1112.7`),
and the stated `ρ_f = 1113.8` sits between them. So `T ≈ 296 ± 2 K` is
**derived from the paper, not assumed** — `provenance: derived`, not `assumed`.

⚠ And the same calculation settles a transcription question. Written without the
`π` — `4Δρga³/3`, which is how one automated reading of the paper rendered it —
the identical inputs give `T = 94 K`, i.e. −179 °C. Not a temperature. The `π`
belongs there, and the check is what says so rather than a recollection of the
formula.

Viscosity then follows from the **measured** `D₀` by Stokes–Einstein:
`η = k_B T / (6π a D₀)` = **1.22 mPa·s** at 20 °C, 1.24 mPa·s at 25 °C. That is
the right order for a decalin/TCE mixture and slightly low, which is what a
`D₀` measured in a finite cell would give.

★ **But none of this is needed for the reduced system.** `(a, l_g, Φ)` fix the
dimensionless problem and `τ_d = a²/D₀` fixes the clock, and all four are
measured. `T` and `η` are needed only to report an SI time, so the run is
anchored in measurement even where the paper is silent.

## 4 · ⚠ Two definitions of `Pe_g`, differing by a factor 2

The text defines

> *"the ratio of gravitational potential energy to thermal energy, given by
> `Pe_g = 2a/l_g << 0.1`"*

while **Table 1's column header reads `Pe = a/l_g`**, and the tabulated values
(0.037 … 0.046) are reproduced by `a/l_g` (0.4/10.7 = 0.0374), not by `2a/l_g`
(0.0748).

Recorded rather than resolved, per rule 5. It does not touch anything downstream
here: `l_g` and `a` are both stated directly, and `Pe_g` is a label derived from
them. **Do not cite a `Pe_g` from this paper without saying which definition.**

## 5 · What can be checked against this paper

| # | claim | our role for it | how |
|---|---|---|---|
| 1 | `Φ(z) ∝ exp(−z/l_g)` in the dilute limit | `implementation_check` | exact for point particles; a fitted decay length that misses `l_g` is a bug |
| 2 | the zero-field EOS is hard-sphere, Carnahan–Starling | `hypothesis` | the simulation imposes a WCA core, **not** a hard sphere, and imposes no EOS — so a mismatch is a result |
| 3 | `K_L = −2.80 ± 0.05` in `D = D₀(1 + K_L Φ)` | `hypothesis`, and **not reachable** | `K_L` is a *hydrodynamic* correction; BD with a scalar `γ` has no hydrodynamic interactions, so the simulation cannot produce it. Stated here so it is not mistaken for a target |

★ Row 3 is the honest limit of this comparison and belongs in the prediction
document, not in a footnote: the paper's own headline dynamical result is
**outside** what free-draining Brownian dynamics can reproduce. Row 2 is the
discovery-capable comparison; row 1 is the correctness check that licenses it.

## 6 · What the paper does not give

- **temperature** — derived in §3 instead
- **solvent viscosity** — derived from `D₀` in §3 instead
- **particle count** in the observation volume
- **tabulated `Φ(z)`** — the profiles are figures only, and were **not digitized**
  (see `docs/05-pitfalls.md` on a citation to a digitizer that was never
  committed; this distillation does not repeat that)
- any percentage for the *"excellent agreement"* with Carnahan–Starling — so the
  agreement is a qualitative statement in the source, and our own comparison has
  to set its own threshold

## See also

`knowledge/wiki/findings/wca-reproduces-carnahan-starling.md` ·
[`verify/verify_sedimentation_wall.py`](../../../verify/verify_sedimentation_wall.py)
