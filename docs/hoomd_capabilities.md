# HOOMD capability matrix — the physics modules we can build

> Surveyed 2026-08-03 · **hoomd 7.1.0** (the installed build; the docs are v7.1.1 — a patch difference, and the APIs surveyed are identical)
> Environment: macOS 26.5 arm64, the `simulation_bot` conda env
> Build: `gpu=False`, `mpi=False`, `compile flags: DOUBLE[SINGLE]`
> Method: introspection of the installed build (`scratch/survey.py`) + live smoke tests (`scratch/smoke.py`)
>
> **This document is measurement, not guesswork.** The basis for the module-registry design of master plan §5.6.

---

## 0. Summary

- **Leaving the scope open was right.** HOOMD is far broader than expected:
  28 isotropic pair potentials, 17 anisotropic (non-spherical) pair potentials, bonds/angles/dihedrals, frictional contact, many-body potentials,
  meshes (membranes and vesicles), long-range electrostatics, curved-surface constraints (manifolds), rigid bodies, HPMC, and MPCD.
- **All 5 cases in `intake/` are implementable** (§3).
- **There is one hole**: an optical tweezer (a harmonic trap) is not built in → implement it directly with `md.force.Custom`.
  Confirmed working, and it passes a golden physics verification too (§4).
- **The most important constraint**: **in BD, a rod's anisotropic translational friction cannot be produced from geometry alone** — bind it
  as a rigid body or as a bead chain and either way `γ⊥/γ∥ = 1.000000`. Anisotropy is an effect of hydrodynamic interactions (HI), and BD has
  no HI. That is not a limitation of HOOMD but a property of the model itself (§5.1). Rotational friction does
  work correctly as a tensor (§5.2). Otherwise: no MPI (a single core), no GPU.

---

## 1. Live verification results (15 APIs)

| # | Item | Result | Corresponding design |
|---|---|---|---|
| 1 | BD + WCA 2D | ✓ | §11, appendix A |
| 2 | `GSD(logger=)` storing per-particle forces | ✓ a `(64,3)` array confirmed | §9.2 Tier B |
| 3 | `write.Burst` sliding window + `dump()` | ✓ (but `write_at_start=True` is required) | §9.2 Tier C |
| 4 | a `md.force.Custom` harmonic trap | ✓ | `external.harmonic_trap` |
| 5 | `variant.Ramp` / `Cycle` time-dependent driving | ✓ | `driving.*` |
| 6 | a `bond.Harmonic` + `angle.Harmonic` chain | ✓ | `bonded.*` |
| 7 | `force.Active` + `create_diffusion_updater` | ✓ | `active.abp` |
| 8 | a custom Action updater (run-and-flip) | ✓ 5 flips (3.6 expected) | `active.run_and_flip` |
| 9 | `pair.Table` for an arbitrary r⁻ⁿ | ✓ | `pair.table` |
| 10 | `pair.aniso.GayBerne` + per-axis `gamma_r` | ✓ aspect 3.0 | `shape.ellipsoid` |
| 11 | `methods.OverdampedViscous` (deterministic) | ✓ | for golden tests |
| 12 | `write.HDF5Log` global scalars | ✓ | §9.2 Tier L |
| 13 | restart GSD (`truncate=True`) | ✓ 1 frame retained + reloaded | §14 checkpoints |
| 14 | runtime monitoring from a custom Action | ✓ | §12.5 runtime guards |
| 15 | `GSD(filter=Tags)` a tracked particle subset | ✓ 8/64 | §9.2 Tier D |

**All 5 raw-data tiers (§9) were confirmed working.**

### Harmonic-trap golden physics verification (`scratch/golden_trap.py`)

The equilibrium distribution of a BD particle in a 2D harmonic trap is `⟨x²⟩ = k_BT/k`, with relaxation time `τ = γ/k`.
N=400, 340,000 steps each:

| k | τ=γ/k | dt | ⟨x²⟩ measured | kT/k predicted | error |
|---|---|---|---|---|---|
| 2.0 | 0.500 | 2.5e-4 | 0.50188 ± 0.0016 | 0.50000 | +0.38% |
| 5.0 | 0.200 | 1.0e-4 | 0.20112 ± 0.0006 | 0.20000 | +0.56% |
| 10.0 | 0.100 | 5.0e-5 | 0.09998 ± 0.0003 | 0.10000 | −0.02% |
| 20.0 | 0.050 | 2.5e-5 | 0.04993 ± 0.0001 | 0.05000 | −0.14% |

`⟨x²⟩·k = [1.0038, 1.0056, 0.9998, 0.9986]` — change k by 10× and it stays at kT (coefficient of variation **0.28%**).

→ **The BD integrator is accurate, and `trap-2d-5um` can be adopted as a golden physics test.**

---

## 2. ⚠️ Newly discovered traps (to be added to the master plan §11 trap list)

### Trap 7 — an external force + a periodic boundary: leave the minimum image out and it is quietly wrong ⭐️

When writing a trap towards a fixed anchor, using `d = pos - anchor` as it stands means that the moment the particle
wraps across the box the distance jumps by L and it receives an enormous restoring force **in the wrong direction**.

The symptom is nasty — **it does not blow up, it is quietly wrong.** With a strong trap (a low probability of reaching the boundary)
the correct value comes out, and the weaker it gets the larger the error:

| k | without the minimum image | with the minimum image |
|---|---|---|
| 2.0 | +1856% ✗ | +0.38% ✓ |
| 5.0 | +344% ✗ | +0.56% ✓ |
| 10.0 | +0.16% ✓ | −0.02% ✓ |

Had only k=10 been tested, it would have been mistaken for a pass. **It applies to every `external.*` module**, and
the validator needs a rule that "an external-force module has to declare whether it applies the minimum image".

```python
d = pos - anchors[tags]
d -= L * np.round(d / L)      # ← this one line
```

### Trap 8 — a NaN in the z component of the minimum image in 2D

In the code above, leaving the box length in z as `np.inf` gives `inf * round(0/inf) = inf*0 = nan`, so
the z force becomes NaN. Being 2D, HOOMD ignores z so the result was correct, but **a NaN entering the force array
is itself dangerous** (the §12.5 runtime guard can raise a false positive). Apply it with a mask, to the periodic dimensions only.

### Trap 9 — `write.Burst` needs `write_at_start=True` for a new file

Without it, `RuntimeError: Must set write_at_start to write to a new file.`
The consequence is that one extra initial frame goes into the file (a buffer of 10 + 1 at the start = 11).

---

## 3. Mapping the 5 `intake/` cases (all implementable)

| Case | Physics needed | HOOMD equivalent | Difficulty |
|---|---|---|---|
| `trap-2d-5um` | a harmonic trap | `md.force.Custom` (not built in) ✓verified | low |
| `trap-drag-2d-hex300` | a trap + moving driving | `force.Custom` + `variant.Ramp/Cycle` ✓verified | low |
| `chain-bend-2d-oscill` | bonds + bending + oscillation | `bond.Harmonic` + `angle.Harmonic` + `variant.Cycle` ✓verified | low |
| `soft-r3-2d-A-sweep` | a soft r⁻³ | `pair.Table` (or `pair.Mie`) ✓verified | low |
| `abp-rod-2d-run-flip` | an ellipsoid + run-and-flip | `pair.aniso.GayBerne`/`ALJ` + a custom updater ✓verified | **medium** |

`abp-rod` is the only awkward one — because of the anisotropic translational friction problem (§5.1~5.4).
But being 2D, rotation has only the one z axis, and the MSAD and the long-time MSD are reproduced exactly (§5.3).

---

## 4. The full capability list (measured)

### Isotropic pair potentials — `md.pair` (28)
`Buckingham` `DLVO` `DPD` `DPDConservative` `DPDLJ` `Ewald` `ExpandedGaussian`
`ExpandedLJ` `ExpandedMie` `ForceShiftedLJ` `Fourier` `Gaussian` **`LJ`** `LJ0804`
`LJ1208` `LJGauss` `Mie` `Moliere` `Morse` `OPP` `ReactionField` `TWF` **`Table`**
`WangFrenkel` **`Yukawa`** `ZBL` `Zetterling`

> There is still no dedicated WCA class → `LJ(r_cut=2^(1/6)σ, mode='shift')` (§11 trap 1)

### Anisotropic pairs (non-spherical particles) — `md.pair.aniso` (17)
**`GayBerne`** (ellipsoids) · **`ALJ`** (anisotropic LJ for polyhedra/ellipsoids) · `Dipole` · `YLZ` ·
`Patchy` `PatchyLJ` `PatchyGaussian` `PatchyMie` `PatchyYukawa`
`PatchyExpandedLJ` `PatchyExpandedGaussian` `PatchyExpandedMie`

### Frictional contact — `md.pair.friction` ✓ measured (2026-08-06)
`FrictionLJCoulomb` `FrictionLJCoulombNewton` `FrictionLJLinear`
(Hofmann et al. 2025, [arXiv:2507.16388](https://doi.org/10.48550/arXiv.2507.16388))

`params = {epsilon, sigma, kT, gamma_f?, kappa_f?}` · the normal is fixed to WCA ·
the tangential force = `w(r)·f(u)`, `w(r) = −dU_WCA/dr`, `u` = the relative surface velocity at the contact point ·
it produces a force **and a torque** together. Measured (`scratch/verify_pair_friction.py`, 23/23, kT=0, V=1, r*=0.95):

| Model | `f(u)` | F_tan | residual slip u/V=1e−4 |
|---|---|---|---|
| `FrictionLJLinear` | `γ_f·u` | 6.1926 | 6.19e−4 (∝u) |
| `FrictionLJCoulomb` | `κ_f` | 18.578 | **18.578 (does not shrink)** |
| `FrictionLJCoulombNewton` | `min[γ_f u, κ_f]` | 1.0 | 1e−4 (∝u) |

**Three structural properties** (all measured):
- **Exactly 0 outside the cutoff** — `w(r)=0`. The DLVO secondary minimum (`r*=1.00759`) is outside it, so `F_tan=|τ|=0.000e+00`.
- **Dissipative** — even with a tangential displacement, the force is 0 if `u=0` (the normal WCA is still alive at 6.1926).
  It contributes to `K″` only and 0 to `K′` → it vanishes in the quasi-static limit. There is no static-friction (stick) state.
- **Exempt for slip-free rolling** — at `ω=V/2R` all three are exactly 0. It cannot prevent bending (= rolling).

> ⛔ **It is not a substitute for the JKR bending stiffness.** What is needed is a tangential spring with history + a rolling resistance
> (Cundall–Strack / Dominik–Tielens), and **HOOMD does not have it** → `force.Custom(aniso=True)`.
> ★ Read traps 19 and 20 (the lever arm `R = diameter/2` · `Brownian` velocities) in skill `bd-hoomd` first.

### Bonded interactions
| Kind | Classes |
|---|---|
| `md.bond` | **`Harmonic`** **`FENEWCA`** `Table` `Tether` |
| `md.angle` | **`Harmonic`** `CosineSquared` `Table` |
| `md.dihedral` | `OPLS` `Periodic` `Table` |
| `md.improper` | `Harmonic` `Periodic` |
| `md.special_pair` | `Coulomb` `LJ` |

### External forces
| Kind | Classes | Notes |
|---|---|---|
| `md.external.field` | `Electric` `Magnetic` `Periodic` | **no harmonic trap** |
| `md.external.wall` | `LJ` `Gaussian` `Morse` `Mie` `Yukawa` `ForceShiftedLJ` | + `hoomd.wall`: `Plane` `Sphere` `Cylinder` |
| `md.force` | **`Custom`** **`Constant`** **`Active`** `ActiveOnManifold` | `Custom` is the main route for extensions |

### Integration methods — `md.methods`
**`Brownian`** `Langevin` **`OverdampedViscous`** `ConstantVolume` `ConstantPressure`
`DisplacementCapped`
`md.methods.rattle`: the methods above + **manifold constraints** (BD on a curved surface)
`md.methods.thermostats`: `Berendsen` `Bussi` `MTTK`

### Updaters and constraints
`md.update`: **`ActiveRotationalDiffusion`** `ReversePerturbationFlow` `ZeroMomentum` `Mesh*`
`md.constrain`: **`Rigid`** `Distance`
`hoomd.update`: **`BoxResize`** `CustomUpdater` `FilterUpdater` `RemoveDrift`

#### `BoxResize` measured (`scratch/verify_3d_boxresize.py`, 2026-08-06) ⭐️

```
hoomd.update.BoxResize(trigger, box: variant.box.BoxVariant, filter=All())
hoomd.variant.box.Interpolate(initial_box, final_box, variant)      # variant: 0→1 scalar
hoomd.variant.box.InverseVolumeRamp(initial_box, final_volume, t_start, t_ramp)
```

| What was checked | Result |
|---|---|
| does it affinely scale the coordinates (`r → s·r`) | **yes — maximum error 8.9e-16** (exactly as the docs claim) |
| does it stop exactly at the target size | yes (`Lx` error < 1e-6) |
| does it coexist with `Brownian` + `pair.Table` + WCA | yes. No crash, a finite PE, the cell list updated automatically |
| `r_cut < L/2` after compression | has to be checked separately (compression can break this condition — trap 6) |

★ **The affine scaling is itself the trap** — compression **also shrinks the bond length of pairs that are already bonded.**
In a narrow well (the DLVO secondary minimum), if one trigger's strain exceeds the well width the pair is
pushed inside the barrier and **collapses irreversibly into the primary minimum (contact)**. For the measured threshold → §5.5 below.

### Many-body, long-range and meshes
`md.many_body`: `Tersoff` `SquareDensity` `RevCross`
`md.long_range.pppm`: `Coulomb`
`md.mesh`: `bending.Helfrich` `bending.BendingRigidity` `bond.*` `conservation.{Area,Volume,TriangleArea}`
> Membrane and vesicle simulations are possible

### Time-dependent driving — `hoomd.variant`
`Constant` **`Ramp`** **`Cycle`** `Power`
`variant.box`: `Interpolate` `InverseVolumeRamp`

### Output — `hoomd.write`
**`GSD`** **`Burst`** **`HDF5Log`** `Table` `DCD` `CustomWriter`

### Filters — `hoomd.filter`
`All` **`Tags`** `Type` `Rigid` `Union` `Intersection` `SetDifference` `Null` `CustomFilter`

### Extensions along other axes
`hoomd.hpmc` — hard-particle Monte Carlo (`Ellipsoid` `ConvexPolyhedron` `Sphinx` … + `external.Harmonic`)
`hoomd.mpcd` — multiparticle collision dynamics (hydrodynamic coupling)
`md.alchemy` — alchemical transformations
`md.tune` / `hoomd.tune` — auto-tuning (`NeighborListBuffer` `ParticleSorter` `GridOptimizer`)

---

## 5. Constraints (must be reflected in the module design)

| Constraint | Content | Impact |
|---|---|---|
| **translational friction is a scalar — no way around it** ⭐️ | Confirmed by measurement in §5.1 below. Neither a rigid body nor a bead chain can produce anisotropy | anisotropic translational diffusion has to be **imposed explicitly**. The `shape.*` module declares the limitation |
| **rotational friction is a tensor — works correctly** ✓ | `default_gamma_r=(x,y,z)` confirmed by measurement (§5.1) | the MSAD and the rotational dynamics can be reproduced exactly |
| **no MPI** | `mpi_enabled=False` | one run = one core, settled. Parallelism only in "run several runs at once" (§2 principle 6 still holds) |
| **no GPU** | `gpu_enabled=False` | there is no CUDA build for osx-arm64. The scale ceiling stands |
| **no built-in harmonic trap** | not in `md.external.field` | `force.Custom` is mandatory + trap 7 (the minimum image) must be observed |
| **`hpmc.external.Harmonic` cannot be used** | it is HPMC (Monte Carlo) only | not compatible with the MD/BD path |

### 5.1 Anisotropic translational friction — the measured result ⭐️

**Conclusion: in BD, a rod's anisotropic translational friction does not come out of geometry alone. This is not a limitation of HOOMD but
a property of the BD model itself.**

Both routes were measured (`scratch/rigid_rod_friction.py`, `scratch/anisotropy_probe.py`).
A constant force is applied to a rod held at a fixed orientation and `γ = F/v` is measured directly from the terminal velocity (`OverdampedViscous`, no noise):

| Configuration | γ∥ | γ⊥ | γ⊥/γ∥ | Verdict |
|---|---|---|---|---|
| **rigid body** (`constrain.Rigid`, 5 spheres) | 1.00000 | 1.00000 | **1.000000** | isotropic |
| **free-draining bead rod** (5 spheres tied by bonds+angles) | 5.00000 | 5.00000 | **1.000000** | isotropic |
| theoretical expectation (slender body, with HI) | — | — | → 2 | — |

A further check: giving the rigid body's constituent particle type `gamma['A'] = 10` does not change the translational velocity at all
(v∥ = 1.000000, identical). **A rigid body's translational drag is determined entirely by the single scalar `gamma` of its central particle.**

**Why**: in the free-draining approximation each bead feels an independent Stokes drag, so the
total drag is `N·γ_bead`, independent of direction. A real rod's `γ⊥/γ∥ → 2` is due to **the hydrodynamic interaction (HI) by which the beads
shield each other's flow field**, and the shielding is stronger for axial alignment, giving `γ∥ < γ⊥`.
**In BD, which has no HI, no amount of geometric refinement produces anisotropy.**

The measured `γ = 5.000 = N·γ_bead` agrees with that explanation exactly.

### 5.2 The rotational friction tensor — works correctly ✓

Specifying `default_gamma_r=(x,y,z)` per axis works perfectly. The rotation angle measured for a constant torque:

| `gamma_r` | rotation angle for τ∥x | rotation angle for τ∥z | ratio z/x | expected |
|---|---|---|---|---|
| (1, 1, 1) | 2.00000 rad | 2.00000 rad | 1.0000 | 1.0 ✓ |
| (1, 1, 5) | 2.00000 rad | 0.40000 rad | **0.2000** | 0.2 ✓ |

→ **The MSAD and the rotational dynamics can be reproduced exactly.**

### 5.3 The actual impact on `abp-rod-2d-run-flip`

Being 2D, rotation has only the one z axis — `gamma_r` anisotropy is not needed in the first place.

| Observable | Depends on | Accuracy |
|---|---|---|
| **MSAD** | `γ_r,z` only | ✓ exact |
| **MSD (long time, t ≫ τ_r)** | the isotropic average `γ̄` only (the direction averages out) | ✓ exact |
| **MSD (short time, t < τ_r)** | `γ∥`, `γ⊥` individually | ✗ **anisotropy lost** |

Run-and-flip holds its orientation for a long time between flips, so the short-time anisotropy really can be observed.
If that signal is the goal, one of the options below is needed.

### 5.4 The options

| Option | Method | Cost | Accuracy |
|---|---|---|---|
| **A. the isotropic average γ̄ (the recommended default)** | compute `γ̄ = 3/(1/γ∥ + 2/γ⊥)` from the aspect ratio using Perrin/slender-body theory. Rotation exactly, via the `γ_r` tensor. The module declares the limitation | low | the long-time MSD and the MSAD are exact, there is no short-time anisotropy |
| **B. impose the anisotropy with a custom integrator** | the deterministic term is possible as a corrective force, but **the noise term is the problem** — HOOMD generates isotropic noise, so the fluctuation-dissipation theorem breaks. Needs the `md.half_step_hook` route | high | exact |
| **C. introduce HI via MPCD** | `hoomd.mpcd` (multiparticle collision dynamics) — HI arises naturally | very high | exact, a project along a different axis |

→ **A is left as the default, and the `shape.ellipsoid` module declares
`translational_friction: "isotropic_average"`**.
The validator can raise a warning that "the goal is short-time MSD anisotropy but an isotropic approximation is in use".
If B becomes necessary it is added as a separate module (the advantage of the module structure).

### 5.5 Box compression breaks the bonds of a narrow well — the measured threshold ⭐️

`scratch/verify_3d_boxresize.py` (2026-08-06). The `network` case has an "aggregate then squeeze"
protocol, so this was measured before running it. Units d=kT=γ=1.

The DLVO ledger (inherited from `chain-bend-2d-dlvo`): the secondary minimum `h*=0.007593`, the barrier `h*=0.000508`,
the bond length `ℓ*=1.007593`, `k_bond*=1.042e6 kT/d²`, `τ_bond*=9.594e-7 τ_B`.

**The analytic prediction**: if an affine step pushes the pair inside the barrier it cannot come back, so
```
ε_crit = (h_min* − barrier_h*) / ℓ*  =  0.703 %  (linear strain per trigger)
```

**Measured** (kT=0 deterministic, trigger interval = 21 τ_bond of relaxation):

| ε/trigger | final `h*` | Result |
|---|---|---|
| 0.40 % | 0.007591 | the secondary minimum holds ✓ |
| 0.80 % | **−0.042257** | collapse — contact/overlap |
| 2.00 % | −0.042519 | collapse |

**The measured threshold of 0.40~0.80% brackets the predicted 0.703%.** After collapse, `h*≈−0.042` (a 4.2% overlap) is
the seat where the vdW divergence and the WCA core balance, and it is **irreversible** — a reversible secondary-minimum bond
turns into a permanent contact.

**The cost that follows**: φ 0.02→0.10 is a total linear strain of 41.5%, so dividing by ε=0.4% gives
**at least 134 stages**. Making each stage quasi-static (≫τ_B) makes it ~400 τ_B, and
at `dt = 1e-2·τ_bond = 9.6e-9 τ_B` that is 4.2e10 steps — against the
**6384 steps/s** measured at N=1528 that is **1800 hours**. → Quasi-static compression is impossible at this dt.
Fast compression (with only τ_bond-scale relaxation between stages) is possible at 134×2000 = 2.7e5 steps ≈ **42 seconds**, but
the structure is then **the affinely compressed initial arrangement**, with no diffusion intervening — which is different physics.

---

## 6. How to reproduce

```bash
conda env update -f environment.yml -n simulation_bot --prune
CONDA=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
$CONDA scratch/survey.py     # the capability survey (introspection of the installed build)
$CONDA scratch/smoke.py      # 15 APIs, live
$CONDA scratch/golden_trap.py         # harmonic-trap golden physics verification (~2 min)
$CONDA scratch/rigid_rod_friction.py  # rigid-rod anisotropy (comes out isotropic)
$CONDA scratch/anisotropy_probe.py    # gamma_r anisotropy + a free-draining bead rod
```
