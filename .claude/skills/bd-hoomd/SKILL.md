---
name: bd-hoomd
description: |
  HOOMD-blue v7.1 API mapping, execution-verified code snippets, and 20 traps
  found by measurement. Read this before writing any Brownian dynamics
  simulation code — several of the traps produce silently wrong results with no
  error at all (+1856% error / for a nearly straight chain, angle.Harmonic's
  **force** is up to 96% wrong while its energy stays exact, so an energy check
  cannot find it / BoxResize compression irreversibly collapses bonds in a
  narrow well). Applies to writing HOOMD code, implementing traps, active
  particles, chains and non-spherical particles, **3D and box compression
  (gelation)**, trajectory and force storage, restarts, and runtime monitoring.
---

# HOOMD-blue v7.1 — verified usage

> Every snippet and every number in this document **was actually executed against
> the installed build and confirmed.**
> Evidence: the scripts in `verify/`, and
> [docs/hoomd_capabilities.md](../../../docs/hoomd_capabilities.md).
> Do not add a guess. A new fact goes in here only after execution confirms it.

## Environment

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```
hoomd **7.1.0** · `gpu=False` · **`mpi=False`** · macOS arm64 · double precision

**With no MPI, one run = one core.** Parallelism comes only from running several
runs at once.

---

# ⚠️ Traps — check before writing code

**★ marks the ones that are silently wrong with no error.** Those are the most
dangerous.

### ★ 1. External force plus periodic boundaries → minimum image is mandatory

In a force toward a fixed anchor, using `d = pos - anchor` directly means that
the instant a particle wraps across the box, the distance jumps by L and it
receives **an enormous force in the opposite direction. It does not blow up.**
A stiff trap gives the right answer and it gets worse as the trap weakens:

| k | without minimum image | with |
|---|---|---|
| 2 | **+1856%** ✗ | +0.38% ✓ |
| 5 | +344% ✗ | +0.56% ✓ |
| 10 | +0.16% ✓ | −0.02% ✓ |

Test only the stiff condition and you will believe it passed. **Verify in the
weak condition.**

### ★ 2. Brownian carries O(δt) error

The documentation states this. It needs a much smaller `dt` than
Langevin/ConstantVolume. With σ=kT=γ=1, start from `dt ≈ 1e-4` and shrink further
when the forces are strong (large Pe, large ε, large k).
**Passing the upper bound does not guarantee accuracy** — a publication-grade
result needs a separate `dt`-halved convergence check.

### ★ 3. ABP rotation goes through an updater, not the integrator

```python
integrator.integrate_rotational_dof = False        # <- mandatory
sim.operations.updaters.append(
    active.create_diffusion_updater(trigger, rotational_diffusion=D_r))
```
Leave it `True` and inertial rotation mixes in, so it stops being an ABP. It
quietly becomes different physics.

### ★ 4. There is no dedicated WCA class

```python
md.pair.LJ(nlist=cell, default_r_cut=2**(1/6)*sigma, mode='shift')   # = WCA
```
`ForceShiftedLJ` is **not** WCA (it is a different potential). Do not be misled
by the name.

### ★ 5. BD is overdamped — velocity has no physical meaning

`thermalize_particle_momenta()` is unnecessary. **No velocity-based MSD or VACF.**
Use position differences only.

### 6. `r_cut < L/2`

Otherwise the minimum image convention is violated. This bites often in a small
box.

### 7. NaN in the z component of a 2D minimum image

Setting a non-periodic axis's period to `np.inf` gives
`inf * round(0/inf) = inf*0 = nan`. A NaN enters the force array and the runtime
guard false-positives. **Mask the periodic axes only** (snippet below).

### 8. `write.Burst` needs `write_at_start=True` for a new file

Without it: `RuntimeError: Must set write_at_start to write to a new file.`
The result file then contains one extra initial frame (buffer 10 → file 11).

### 9. 2D means `Lz=0` **and** `dimensions=2`

Set the box z to 0 and `configuration.dimensions = 2`. Both are required.

### ★ 10. `pair.Table`'s grid is `endpoint=False`

The documentation: *"implicitly defined r values are those returned by
`numpy.linspace(r_min, r_cut, len(U), endpoint=False)`"*.
Build it with `endpoint=True` and the whole table shifts, so **the force is
silently wrong.**
Measured (`verify/verify_pair_table.py`, `U=A/r³`, r_min=0.5, r_cut=3.0, 200 bins):

| sep | `endpoint=False` | `endpoint=True` |
|---|---|---|
| 0.70 | +0.000% | −0.572% |
| 1.50 | +0.000% | −1.329% |
| 2.90 | +0.000% | **−1.646%** |

`endpoint=False` is **exact to 0.000%**. The error grows as you approach the
cutoff.

### ★ 11. `pair.Table` gives force and energy **exactly 0** for `r < r_min`

The documentation says `F=0, U=0 for r < r_min`, and **for a diverging potential
that is a trap.** The repulsion disappears, so once a particle penetrates
`r_min` it simply stays overlapped. It does not blow up.

Measured (`A/r³`, r_min=0.5): sep=0.51 → F=443.9 (normal) / sep=0.49 →
**F=0.0, U=0.0**.

→ Two layers of defence: put **the excluded-volume core in as a separate force**
to prevent `r_min` ever being reached (WCA and Table sum exactly — measured error
0.000%), and **monitor the minimum neighbour distance** during the run.

### ★ 12. `seed` is truncated to 16 bits — two different seeds can be the same run

Give `Simulation(seed=...)` a value above 65535 and it is truncated to
`seed % 65536` with a warning. The warning appears but **the run proceeds**, so
spacing the seeds of repeated runs by 65536 produces **completely identical
trajectories** — and you average them believing they are independent samples.

Measured (free BD, N=40, 2000 steps, max coordinate difference):

| seed pair | difference |
|---|---|
| `20260803` vs `10179` (= mod 65536) | **0.000e+00** |
| `20260803` vs `20260803+65536` | **0.000e+00** |
| `20260803` vs `10180` | 9.97 (properly different) |

→ Use **small consecutive integers** (1, 2, 3, …) for repeated runs.

---

### ★ 13. With `active_force = 0`, `ActiveRotationalDiffusion` **does not run**

The orientation never decorrelates — the `⟨n(0)·n(t)⟩` decay rate is
**Λ = 0.0000** (measured across all four combinations).

```python
act.active_force["A"] = (0.0, 0.0, 0.0)     # <- this turns rotational diffusion off too
```

**So an attempt to turn the active force off and look at the rotational
statistics alone is blocked.** To verify orientational diffusion by itself, leave
a small active force on and measure the director
(`verify/standalone_abp_diffusion.py`).

### ★ 14. `rotational_diffusion` **is the director decay rate Λ itself** — off by 2× in 3D

Give HOOMD `rotational_diffusion=D_r` and the measured result is:

| dim | `Λ/D_r` | `Λ/[(d−1)D_r]` (standard theory) |
|---|---|---|
| 2 | **1.00** | 1.00 ✓ |
| 3 | **1.00** | **0.50** ✗ |

So `⟨n(0)·n(t)⟩ = exp(−D_r·t)` holds in **both 2D and 3D**.
Standard spherical rotational diffusion is `exp(−(d−1)D_r^phys·t)`, so **to
reproduce a physical `D_r^phys` in 3D you must pass
`rotational_diffusion = 2·D_r^phys`.** 2D is unchanged.

The consequence — free ABP effective diffusion:
```
D_eff = D_t + v₀² / (d · Λ)          Λ = HOOMD's rotational_diffusion
```
Measured to within 1.5% (2D and 3D, two values of v₀). The commonly used
`v₀²/[2(d−1)D_r]` **coincides only in 2D, by accident**, and is +29–31% wrong in
3D.

> Every case here is 2D, so this has no effect today. **It bites the moment
> anyone does 3D active matter.**

---

### ★★ 15. `angle.Harmonic` gets the energy right and **the force wrong** on a nearly straight chain

`md.angle.Harmonic` uses `1/sin θ` to convert torque into Cartesian force, and
**clamps `sin θ` from below**. Measured and pinned:
**SMALL = 1.414217e-03 = √2×10⁻³** (s.d. 1.4e-7, across six amplitudes).

```
sin θ < SMALL  ->  the force is scaled by exactly  sinθ/SMALL
```

With `t0 = π` (a straight chain is the equilibrium), **the equilibrium itself is
the singularity at sin θ = 0**, and a stiff chain is always nearly straight. So
this trap bites **hardest on a rigid filament.** Since `sin θ ≈ |θ−π|`, the force
goes as `∝ κ(θ−π)·(θ−π)/SMALL = κ(θ−π)²/SMALL` — **quadratic, not linear.** The
chain becomes far **softer and nonlinear** than it should be, and it worsens the
straighter it gets.

Measured (κ_θ = 1.39e6 kT, n=5, smooth bending mode):

| max\|θ−π\| | energy error | force ratio (vs correct) |
|---|---|---|
| 1.757e-02 | 0.0000% | **1.000000** ✓ |
| 1.757e-03 | 0.0000% | **1.000000** ✓ |
| 8.787e-04 | 0.0000% | 0.621320 ✗ |
| 2.929e-04 | 0.0000% | 0.207107 ✗ |
| 5.858e-05 | 0.0000% | **0.041421** ✗ |

**★ The energy is exactly 0.0000% correct across the whole range.** So
**verifying by energy passes and leaves the force wrong.** This project actually
did that: `verify/verify_angle_matrix.py` measured `k = 2U/δ²` and **passed at
0.55%**, missing the force bug (energy had been chosen deliberately to dodge a
particle-ordering problem, and that choice concealed the trap).

> **Lesson: verifying a potential by energy alone is not verifying the force.**
> A new potential must have its **force compared directly against the numerical
> gradient of the energy.**

**There is no good workaround** (both excluded by measurement or derivation):
- Giving `angle.Table` the same `U` and `tau` is also wrong (force ratio 3.27 —
  the table resolution Δθ=3.1e-3 cannot represent `|θ−π|~3e-4`). Widening it
  takes the same `1/sin θ` path.
- `angle.CosineSquared`: `U = ½k(cosθ − cosθ₀)²`. At `θ₀=π`,
  `cosθ+1 ≈ (π−θ)²/2`, so `U ≈ k(π−θ)⁴/8` — **quartic.** It cannot produce a
  harmonic bending stiffness.

→ **Responses**: ① keep the system inside `min|θ−π| > 1.41e-3` (every angle! not
just the maximum — the end angles of a chain are an order of magnitude smaller
than the middle) ② implement the bending directly with `force.Custom` (exact but
26× slower — trap 16) ③ **if the regime is linear, do not use MD; solve it
analytically.**

Reproduce: `verify/verify_angle_force_small_theta.py` · the elimination path:
`verify/diagnose_chain_bend_28pct.py` (rules out dt, nonlinearity, the x degrees
of freedom, and energy in turn)

### ★ 16. `force.Custom` runs Python every step — 26× slower on a small system

`set_forces` is called every step. Measured (n=25 chain, 2D, CPU):

| Configuration | steps/s |
|---|---|
| `bond` + `angle` only (compiled) | 61,264 |
| **+ a `force.Custom` trap** | **2,339** (26× slower) |
| **+ a ghost-particle `bond.Harmonic(r0=0)` trap** | **55,551** |

`bond.Harmonic(r0=0)` is `U = ½k r²`, **functionally identical to a harmonic
trap.** Exclude the ghost particle from the integrator's `filter` and it does not
move, giving a fixed trap → entirely on the compiled path.
For time-dependent driving, move the ghost with a `CustomUpdater` (it does not
have to be every step — see trap 17).

→ **Before concluding "too expensive to run", measure whether the cost is
implementation or physics.**
Reproduce: `verify/bench_chain_bend.py`

### ★ 17. Moving the anchor every U steps makes the drive a **zero-order hold**

Updating the trap centre every `U` steps with a `CustomUpdater` turns the driving
sine into a staircase, so the fundamental shrinks by `sinc(ωΔt/2)` and lags by
`ωΔt/2` (`Δt = U·dt`).
Measured `|ŷ_c|/a = 0.98999`, phase `−0.2522 rad`; ZOH prediction `0.99040`,
`−0.2404 rad`.

→ **Do not use the nominal amplitude in a response-function estimator.** Measure
the anchor (ghost) position **as well** and use the measured phasor: the ZOH
attenuation then cancels exactly between numerator and denominator. Using the
nominal value, `K′` at De=10 came out **−6559 (wrong even in sign, 236% error)**;
with the measured phasor, 5863 (21%).
Reproduce: `verify/verify_chain_bend_gates.py --gate lockin`

### ★ 18. `update.BoxResize` **affinely scales** coordinates — bonds in a narrow well break

It scales `r → s·r` as documented (**measured error 8.9e-16**). So compression
**also shortens the bond length of pairs that are already bonded.** When the well
is narrow, as in a DLVO secondary minimum, the instant one trigger's strain
exceeds the well width the pair is pushed inside the barrier and **collapses
irreversibly into the primary minimum (contact).** It does not blow up — the bond
type just quietly changes.

```
ε_crit = (h_min − h_barrier) / ℓ          # allowed linear strain per trigger
```

Measured (DLVO: h_min*=0.007593, barrier*=0.000508, ℓ*=1.007593 → predicted
**0.703%**):

| ε per trigger | final h* | Outcome |
|---|---|---|
| 0.40% | 0.007591 | holds ✓ |
| 0.80% | **−0.042257** | collapsed (4.2% overlap) |
| 2.00% | −0.042519 | collapsed |

→ Before compressing, **compute ε_crit from the ledger and set the step count by
dividing by it.** A total strain of 41.5% (φ 0.02→0.10) split at 0.4% is
**134 steps**.
⚠️ **The step count is set by physics; the relaxation time per step is set by the
purpose** — quasi-static (≫τ_B) costs 1800 hours, relaxing only on the τ_bond
scale costs 42 seconds, and **the structures differ** (a gel that diffusion took
part in, versus an affinely compressed initial configuration). Always state which
one it is.
Also, compression can break `r_cut < L/2`, so **re-check in the final box**
(trap 6).

Reproduce: `verify/verify_3d_boxresize.py` · details:
docs/hoomd_capabilities.md

### ★ 19. `md.pair.friction`'s lever arm `R` comes from `particles.diameter/2`, not `sigma`

The three `FrictionLJ*` variants that produce tangential force and torque take a
`sigma` in `params`, but the contact-point lever arm comes from **the particle
property `diameter`.** `diameter` defaults to **1.0**, so with `sigma≠1` it is
silently inconsistent. Measured (inverting `|τ|/F_tan` for R):

| `particles.diameter` | measured `R` | |
|---|---|---|
| default (1.0), `sigma=0.8909` | **0.5000000000** | **1.1225× too large** vs σ/2 ✗ |
| set explicitly `= sigma` | 0.4454493591 | `= σ/2` exactly ✓ |

There is no error and the force magnitude looks plausible. **Only the torque and
the no-slip rolling condition are wrong.**
→ When using this force, always set `frame.particles.diameter` alongside it.

### ★★ 20. `Brownian`'s `velocity` is not 0 but **uncorrelated thermal noise**

This is a concrete consequence of trap 5 ("BD is overdamped — velocity has no
physical meaning"), but the dangerous direction is the opposite one.
`methods.Brownian` does not set the velocity to zero; it **redraws it from a
Maxwell–Boltzmann distribution every step** (measured `⟨v²⟩ = 2.771` vs
`3kT/m = 3.000`, ratio 0.924, N=80). And that value is **unrelated** to the actual
displacement rate — measured over one free-BD step:

```
stored v        = [-0.634, -0.306, -0.625]
actual Δx/dt    = [14.36, -95.99, -21.53]      ratio = [-0.044, +0.003, +0.029]
```
The overdamped displacement is `√(2Dδt) ∝ δt^½`, so `Δx/dt` diverges as `δt→0`,
while `v` is a separately drawn number.

⟹ **Combining a velocity-dependent pair force (`md.pair.friction`, the DPD
family) with `methods.Brownian` means the friction sees thermal noise, not the
actual relative motion.** The force is not zero but **garbage**, so "a result came
out" will never catch it.

Reproduce: `verify/verify_pair_friction.py` (23/23) · the three variants'
structural properties (dissipative, exempt under rolling, zero outside the
cutoff) are in
[docs/hoomd_capabilities.md](../../../docs/hoomd_capabilities.md)

### Moving to 3D (measured and confirmed) ✓

`network` is this project's first 3D case. What was confirmed:

```python
f.configuration.box = [L, L, L, 0, 0, 0]      # 3D uses Lz=L (2D uses Lz=0 -- trap 9)
f.configuration.dimensions = 3
```
Free 3D BD reproduces `⟨r²⟩ = 6·D·t` (**−1.40%, SEM 2.81%**), with each of the
three axes giving `⟨Δx²⟩=2Dt` (x/y/z = 1.023/0.938/0.997 — no axis dropped).
⚠️ Build unwrapped coordinates by adding `snapshot.particles.image × L`.
⚠️ 3D active matter hits trap 14 (`rotational_diffusion` factor of 2) — this case
has no activity.

---

# Hard constraints (no workaround — measured and confirmed)

### Translational friction is scalar only. Anisotropy cannot be produced

`Brownian(default_gamma=float)` — a **scalar** per type. Only rotation takes a
`default_gamma_r=(x,y,z)` tensor.

Both routes were measured and **both are isotropic**:

| Configuration | γ⊥/γ∥ |
|---|---|
| rigid body via `constrain.Rigid` (5 spheres) | **1.000000** |
| free-draining bead rod (bonds + angles, 5 spheres) | **1.000000** |
| theory (slender body, with HI) | → 2 |

Even giving `gamma['A']=10` to the constituent particles of a rigid body changes
the translational velocity not at all.
**A rigid body's translational drag is the one scalar gamma of its centre
particle.**

**Why**: in free draining, each bead independently feels Stokes drag, so the total
drag is `N·γ_bead` regardless of direction. A real rod's `γ⊥/γ∥→2` is a
**hydrodynamic interaction (HI)** effect, and BD has no HI.
**This is a property of the BD model itself, not a HOOMD limitation.** Do not try
to work around it with geometry.

→ Response: use the isotropic mean `γ̄` and have the module declare the
limitation. Details in `docs/hoomd_capabilities.md`.

### The rotational friction tensor works correctly ✓

| `gamma_r` | rotation angle, τ∥x | rotation angle, τ∥z | ratio |
|---|---|---|---|
| (1,1,1) | 2.00000 | 2.00000 | 1.0000 |
| (1,1,5) | 2.00000 | 0.40000 | **0.2000** |

MSAD and rotational dynamics reproduce exactly.

### There is no built-in harmonic trap

`md.external.field` has only `Electric`, `Magnetic` and `Periodic`.
`hoomd.hpmc.external.Harmonic` is Monte-Carlo only and cannot be used in MD/BD.
→ Implement it directly with `md.force.Custom` (snippet below).

---

# Verified snippets

## A 2D frame + BD + WCA (the basic form)

```python
import itertools, math
import numpy as np, gsd.hoomd, hoomd, hoomd.md as md

def frame_2d(n_side=40, phi=0.5, sigma=1.0):
    """A 2D square lattice at area fraction phi."""
    N = n_side ** 2
    L = math.sqrt(N * math.pi * sigma**2 / (4 * phi))
    a = L / n_side
    pos = np.array([[(i + .5)*a - L/2, (j + .5)*a - L/2, 0.]
                    for i, j in itertools.product(range(n_side), repeat=2)])
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = pos
    f.particles.orientation = [(1, 0, 0, 0)] * N        # needed for active / non-spherical
    f.particles.typeid = [0] * N
    f.particles.types = ['A']
    f.configuration.box = [L, L, 0, 0, 0, 0]            # Lz=0 -> 2D (trap 9)
    f.configuration.dimensions = 2
    return f, N, L

f, N, L = frame_2d()
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_snapshot(f)

cell = md.nlist.Cell(buffer=0.4)
lj = md.pair.LJ(nlist=cell, default_r_cut=2**(1/6), mode='shift')   # WCA (trap 4)
lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)

bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
integrator.integrate_rotational_dof = False
sim.operations.integrator = integrator
sim.run(100_000)
```

## A harmonic trap (with minimum image — resolves traps 1 and 7)

```python
class HarmonicTrap(md.force.Custom):
    """A harmonic trap pulling toward a per-particle anchor.

    Minimum image is applied to the periodic axes only.
    """

    def __init__(self, k, anchors, box_L, dimensions=2):
        super().__init__(aniso=False)
        self.k = float(k)
        self.anchors = np.asarray(anchors, dtype=float)
        # wrap the periodic axes only. 2D's z is left 0 to exclude it
        # (using inf gives inf*0=nan -- trap 7)
        self.period = np.array([box_L, box_L, box_L if dimensions == 3 else 0.0])

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)      # <- tag indexing is mandatory
            d = np.array(snap.particles.position, copy=True) - self.anchors[tags]
            m = self.period > 0
            d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
            arr.force[:] = -self.k * d
            arr.potential_energy[:] = 0.5 * self.k * (d ** 2).sum(axis=1)

    def displacements(self, state):
        d = np.array(state.get_snapshot().particles.position) - self.anchors
        m = self.period > 0
        d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
        return d
```

**`tags` indexing is mandatory** — HOOMD's `ParticleSorter` rearranges memory
order for cache efficiency, so the local snapshot's order is not tag order.

**Verification**: at k=2 (the most fragile), the `⟨x²⟩` error is −0.64%, zero NaNs
in the force array, and the z force is exactly 0. Full verification in
`verify/golden_trap.py` (within 5% at k=2, 5, 10 and 20; `⟨x²⟩·k` varies by
0.28%).

## A moving trap / time-dependent driving

```python
ramp = hoomd.variant.Ramp(A=0.0, B=5.0, t_start=0, t_ramp=1000)
cyc  = hoomd.variant.Cycle(A=-1.0, B=1.0, t_start=0, t_A=10, t_AB=100, t_B=10, t_BA=100)
# inside set_forces, move the anchor as center = f(ramp(timestep))
```
Confirmed: `ramp(0,500,1000) = 0.00, 2.50, 5.00` /
`cycle(0,60,160) = -1.00, 0.00, 0.20`

⚠️ If you move the anchor every `U` steps rather than every step, the drive
becomes a zero-order hold — see trap 17, and measure the anchor position rather
than using the nominal amplitude.

## ABP (continuous rotational diffusion)

```python
active = md.force.Active(filter=hoomd.filter.All())
active.active_force['A'] = (f_a, 0.0, 0.0)      # the particle's local frame
active.active_torque['A'] = (0.0, 0.0, 0.0)

integrator = md.Integrator(dt=dt, methods=[bd], forces=[lj, active])
integrator.integrate_rotational_dof = False      # <- trap 3
sim.operations.integrator = integrator
sim.operations.updaters.append(
    active.create_diffusion_updater(trigger=hoomd.trigger.Periodic(1),
                                    rotational_diffusion=D_r))
```

## run-and-flip / run-and-tumble (discrete events)

**Different from ABP.** This is a Poisson process, not continuous rotational
diffusion. There is no built-in feature, so a custom Action:

```python
class RunAndFlip(hoomd.custom.Action):
    """Flip the orientation by 180 degrees as a Poisson process."""
    def __init__(self, rate, dt, seed=7):
        self.p = rate * dt
        self.rng = np.random.default_rng(seed)
        self.n_flips = 0

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            q = np.array(snap.particles.orientation, copy=True)
            flip = self.rng.random(len(q)) < self.p
            self.n_flips += int(flip.sum())
            # 180-degree rotation about z
            q[flip] = np.column_stack([-q[flip,3], -q[flip,2], q[flip,1], q[flip,0]])
            snap.particles.orientation[:] = q

sim.operations.updaters.append(hoomd.update.CustomUpdater(
    action=RunAndFlip(rate=2.0, dt=dt), trigger=hoomd.trigger.Periodic(1)))
```
`rate·dt ≪ 1` is required for the Poisson approximation to hold
(`dt/τ_flip ≤ 1e-2`).

## A chain (bonds + bending)

```python
f.bonds.N = M - 1
f.bonds.types = ['backbone']
f.bonds.typeid = [0] * (M - 1)
f.bonds.group = [[i, i+1] for i in range(M-1)]
f.angles.N = M - 2
f.angles.types = ['bend']
f.angles.typeid = [0] * (M - 2)
f.angles.group = [[i, i+1, i+2] for i in range(M-2)]

bond = md.bond.Harmonic();  bond.params['backbone'] = dict(k=100.0, r0=1.0)
angle = md.angle.Harmonic(); angle.params['bend'] = dict(k=10.0, t0=math.pi)
```

> ⚠️ **Read trap 15 before using `t0=math.pi`.** On a nearly straight chain
> (`min|θ−π| < 1.41e-3`) **the energy is right and the force is wrong.** A stiff
> chain is always in that regime.

`md.bond.FENEWCA`, `md.angle.CosineSquared` and `md.bond.Table` also exist.

## An arbitrary potential (r⁻ⁿ and so on)

```python
r_min, r_cut, nbins = 0.5, 3.0, 200
r = np.linspace(r_min, r_cut, nbins, endpoint=False)   # * endpoint=False (trap 10)
U = A / r**3;  F = 3*A / r**4
U = U - A / r_cut**3                             # zero at the cutoff. Use r_cut, not r[-1]
tab = md.pair.Table(nlist=cell, default_r_cut=r_cut)
tab.params[('A','A')] = dict(r_min=r_min, U=U, F=F)

# * the force becomes 0 for r < r_min (trap 11) -> add a separate core to prevent reaching it
wca = md.pair.LJ(nlist=cell, default_r_cut=2**(1/6), mode='shift')
wca.params[('A','A')] = dict(epsilon=1.0, sigma=1.0)
integrator = md.Integrator(dt=dt, methods=[bd], forces=[tab, wca])   # the two forces sum
```
Watch the sign of `F = -dU/dr`. Shift `U` but do **not** shift `F`.
**Verification** (`verify/verify_pair_table.py`): the force and energy of the
combination above agree with the analytic solution to **0.000%**
(sep = 0.90 to 2.00, including the WCA region).

## Non-spherical particles

```python
gb = md.pair.aniso.GayBerne(nlist=cell, default_r_cut=4.0)
gb.params[('A','A')] = dict(epsilon=1.0, lperp=0.5, lpar=1.5)   # aspect = 3
bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                         default_gamma=1.0,               # translation is scalar only (see constraints)
                         default_gamma_r=(0.5, 0.5, 2.0)) # rotation can be per-axis
integrator.integrate_rotational_dof = True
```
`md.pair.aniso.ALJ` and `md.constrain.Rigid` also work. **But neither gives
translational anisotropy.**

## The five storage tiers

```python
# Tier A -- positions / orientations (low frequency, all particles)
hoomd.write.GSD(filename='traj_A.gsd', trigger=hoomd.trigger.Periodic(10_000),
                mode='xb', dynamic=['property'])

# Tier B -- per-particle forces / energies (medium frequency)
plog = hoomd.logging.Logger(categories=['particle'])
plog.add(lj, quantities=['forces', 'energies'])
hoomd.write.GSD(filename='traj_B.gsd', trigger=hoomd.trigger.Periodic(100_000),
                mode='xb', logger=plog, dynamic=['property', 'momentum'])
# reading it: frame.log['particles/md/pair/LJ/forces']  -> shape (N, 3)

# Tier C -- a sliding window, flushed to disk only when an event fires
burst = hoomd.write.Burst(filename='burst.gsd', trigger=hoomd.trigger.Periodic(10),
                          mode='xb', max_burst_size=1000,
                          write_at_start=True)      # <- trap 8
# ... when the condition is met:  burst.dump()

# Tier D -- a few tracer particles at high frequency
hoomd.write.GSD(filename='tracers.gsd', trigger=hoomd.trigger.Periodic(10),
                mode='xb', filter=hoomd.filter.Tags(list(range(100))),
                dynamic=['property'])

# Tier L -- global scalars
thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)
glog = hoomd.logging.Logger(categories=['scalar', 'sequence'])
glog.add(thermo, quantities=['potential_energy', 'pressure', 'kinetic_temperature'])
hoomd.write.HDF5Log(trigger=hoomd.trigger.Periodic(1000), filename='log.h5',
                    logger=glog, mode='x')
```
Call `writer.flush()` to be certain it reached the file (especially in a test).

## Restarting

```python
hoomd.write.GSD(filename='restart.gsd', trigger=hoomd.trigger.Periodic(100_000),
                mode='wb', truncate=True)          # always keep exactly one frame
# restore:
sim2.create_state_from_gsd(filename='restart.gsd')
```

## A runtime guard (NaN / energy blow-up)

```python
class Guard(hoomd.custom.Action):
    def __init__(self, thermo):
        self.thermo = thermo
    def act(self, timestep):
        pe = self.thermo.potential_energy
        if pe is not None and not math.isfinite(pe):
            raise RuntimeError(f"non-finite PE at step {timestep}")

sim.operations.writers.append(hoomd.write.CustomWriter(
    action=Guard(thermo), trigger=hoomd.trigger.Periodic(10_000)))
```

## Removing initial overlap

```python
fire = md.minimize.FIRE(dt=1e-4, force_tol=1e-2, angmom_tol=1e-2, energy_tol=1e-7,
                        methods=[...], forces=[lj])
sim.operations.integrator = fire
while not fire.converged:
    sim.run(1000)
sim.operations.integrator = integrator          # swap back to the real integrator
```

---

# API quick reference (measured signatures)

```
md.methods.Brownian(filter, kT, default_gamma=1.0, default_gamma_r=(1.0,1.0,1.0))
md.methods.Langevin(filter, kT, tally_reservoir_energy=False, default_gamma=..., default_gamma_r=...)
md.methods.OverdampedViscous(filter, default_gamma=1.0, default_gamma_r=(1.0,1.0,1.0))
md.force.Active(filter)                  # .active_force[type], .active_torque[type]
md.force.Custom(aniso=False)             # implement set_forces(self, timestep)
md.force.Constant(filter)                # .constant_force[type], .constant_torque[type]
md.update.ActiveRotationalDiffusion(trigger, active_force, rotational_diffusion)
md.pair.LJ(nlist, default_r_cut=None, default_r_on=0.0, mode='none', tail_correction=False)
md.pair.Table(nlist, default_r_cut=None)
md.nlist.Cell(buffer=...)                # Stencil and Tree also exist
md.minimize.FIRE(dt, force_tol, angmom_tol, energy_tol, ...)
hoomd.Simulation(device, seed=None)
hoomd.write.GSD(trigger, filename, filter=All(), mode='ab', truncate=False,
                dynamic=None, logger=None, precision='single')
hoomd.write.Burst(trigger, filename, ..., max_burst_size=-1, write_at_start=False,
                  clear_whole_buffer_after_dump=True)
hoomd.variant.Ramp(A, B, t_start, t_ramp)
hoomd.variant.Cycle(A, B, t_start, t_A, t_AB, t_B, t_BA)
hoomd.filter.Rigid(flags=('center',))    # 'center' | 'constituent' | 'free'
```

`OverdampedViscous` is **overdamped without thermal noise** — use it for
deterministic verification (`v=F/γ`, relaxation `τ=γ/k`). ★ And use it, or
`Langevin(kT=0)`, whenever you are testing an integrator assumption: a thermal
comparison can lack the power to exclude even a 47 % effect.

The full list (28 isotropic pair potentials, 17 anisotropic, bond/angle/dihedral,
friction, many-body, mesh, long-range, manifold, HPMC, MPCD):
[docs/hoomd_capabilities.md](../../../docs/hoomd_capabilities.md)

---

# When you find a new fact

1. Leave a reproduction script in `verify/`
2. Add it to this document (if it is a trap, decide whether it earns a ★ — does
   it fail without an error?)
3. Record the measured numbers in `docs/hoomd_capabilities.md`
4. File a KB entry with `origin: tooling` and a **cause, not a symptom**

**Do not write a guess into this document.** Everything here must have been
confirmed by execution.
