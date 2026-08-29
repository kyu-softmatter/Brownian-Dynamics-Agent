---
name: bd-physics
description: |
  The dimensions-first workflow for Brownian dynamics — unit conventions, how to
  write a scale table for each system, non-dimensionalization, scale-separation
  checks, and inversion. Read this when defining a physical system, proposing
  parameters, choosing a reference scale, interpreting a dimensionless group, or
  judging whether a dt, a box size or a run time is reasonable. Applies to any
  work that carries parameters from a sketch or the literature into a simulation.
---

# The dimensions-first workflow

> **This document grows as cases accumulate.** Each time a new system is driven
> end to end, its scale table is added to section 6 below.
> Do not build a general framework first — generalize only what has appeared
> twice.

## 0. The absolute order

```
① enumerate scales   write down every length, time and energy in the system, in SI, exhaustively
② choose references  pick which of them are the references, and record why
③ derive as ratios   make every dimensionless group "a ratio of two scales" (the interpretation follows)
④ separation check   confirm that the scale you intend to ignore really is separated
⑤ invert             carry the result back into physical units
```

**Fixing the dimensionless groups first and inferring the scales afterwards is
forbidden.** If you must start from a dimensionless value, state the **anchors**
(particle diameter, temperature, viscosity) explicitly and record that those
anchors were chosen arbitrarily.

**Step ② is a physics judgment.** Confirm it with a human until the knowledge base
has depth.

---

## 1. Unit conventions

The default strategy is `thermal` (thermal-fluctuation dominated — colloids, ABPs,
traps):

| Reference | Choice | In HOOMD |
|---|---|---|
| length | particle diameter `d` | `σ = 1` |
| energy | `k_BT` | `kT = 1` |
| time | `τ_B = d²/D_t = 3πηd³/k_BT` | `γ = 1`, `τ_B = 1`, `D_t = 1` |

Other strategies may become necessary (never used yet — add them when they are):
`interaction` (referenced on `ε`, strong coupling) · `active` (on `ℓ_p`, `τ_r`,
very high Pe) · `custom` (a literature convention).

⚠️ **A different reference gives the same system entirely different dimensionless
groups.** When comparing against the literature, check which reference they used.

Do the arithmetic **with units attached, via `pint`**. Strip units only at the
non-dimensionalization step.

### 1.1 Three decisions that must not be conflated ⭐️

| Decision | Varies per system? | What sets it |
|---|---|---|
| **the unit system** (`σ`, `kT`, `τ_B` = 1) | ❌ **fixed** | convention. Comparability with the literature |
| **`dt`** | ✅ **per system** | the **fastest** relevant timescale in that system |
| **reporting units** | ✅ by interest | the scale at which the physics reads (rescaled afterwards) |

**Do not change the unit system per system.** The same system would end up with
entirely different dimensionless groups depending on the reference, making
literature comparison impossible. Fix it at `thermal` (σ, kT, τ_B).

**`dt` cannot be unified.** Fixing the dimensionless `dt* = dt/τ_B` puts systems
up to four orders of magnitude apart. The actual values:

| System | fastest scale / `τ_B` | `dt*` needed for `dt/τ_fast=1e-2` |
|---|---|---|
| harmonic trap (`k*`=6.0e4) | `τ_k/τ_B` = 1.7e-5 | **1.7e-7** |
| ABP (Pe=40) | `τ_v/τ_B` = 2.5e-2 | 2.5e-4 |

`dt*` has to differ by 1500×. Unify on the trap value and the ABP wastes 1500×;
unify on the ABP value and the trap gets `dt/τ_k = 6`, i.e. **one step is 6× the
relaxation time** — completely broken.

**Rescale for reporting.** Plotting a trap result against `t/τ_k` makes the
physics read. That is changing the axis, not the unit system.

### 1.2 How small should `dt` be? — invert it from the bias ⭐️

"Small enough" is vague. **In a linear system the systematic bias can be computed
exactly.**

The Euler–Maruyama discrete stationary variance of a harmonic trap:
```
x_{n+1} = x_n(1 − h) + √(2D·dt)·ξ ,   h ≡ dt/τ_k
⟹  ⟨x²⟩_discrete = (kT/k) / (1 − h/2)
⟹  relative bias = 1/(1 − h/2) − 1  ≈  h/2
```

**So half of `dt/τ` *is* the systematic bias.** Measured and confirmed that HOOMD's
Brownian follows this law (`verify/dt_convergence.py`, N=2000, 1000 samples):

| `dt/τ_k` | bias, measured | bias, theory | difference |
|---|---|---|---|
| 0.10 | 5.262% | 5.263% | −0.002% |
| 0.05 | 2.539% | 2.564% | −0.025% |
| 0.02 | 1.079% | 1.010% | +0.069% |
| 0.01 | 0.569% | 0.503% | +0.067% |

(the +0.07% on the last two is at the level of 1 SEM of statistical error)

**A practical table — invert `dt` from the accuracy you want:**

| Target accuracy | `dt/τ_fast` | Relative cost |
|---|---|---|
| 5% | 1e-1 | 1× |
| 1% | 2e-2 | 5× |
| 0.5% | 1e-2 ← the hard gate | 10× |
| 0.1% | 2e-3 | 50× |
| 0.025% | 5e-4 | 200× |

Cost ∝ `1/dt`. **The hard gate of `1e-2` in section 4 means "0.5% bias."**

⚠️ This closed form holds **only for a linear system.** Nonlinear systems such as
WCA or ABP have no analytic expression, so do a **convergence study** (halve `dt`
and check the observable does not move). The O(dt) scaling does survive, so
two-point Richardson extrapolation works well.

⚠️ And `dt` is not decided once. Whenever a knob can reorder the timescales,
re-derive it — with a trap stiffness scaled by ~200 the trap becomes the *fastest*
mode, and `dt` had not been recomputed.

---

## 2. Derived material properties (for a sphere)

```
γ    = 3πηd                     Stokes drag
D_t  = k_BT / γ                 Stokes-Einstein
τ_B  = d² / D_t                 Brownian (diffusion) time
m    = ρ_p (π/6) d³             particle mass
τ_p  = m / γ                    momentum relaxation (inertia)
v₀   = f_a / γ                  self-propulsion speed
```

**True for a sphere only** — do not use these for other shapes:
```
D_r  = k_BT/(πηd³) = 3 D_t/d²   Stokes-Einstein-Debye (sphere)
```
Ellipsoids and rods need Perrin friction factors. And **translational anisotropy
is not reproducible in BD at all** (no HI — see the hard constraints in skill
`bd-hoomd`).

Density:
```
2D:  φ = N π σ² / (4 L²)        3D:  φ = N π σ³ / (6 L³)
```

### 2.1 The provenance of `η` and `T` — they must be a matched pair ⭐️

`γ`, `D_t` and `τ_B` all hang on `η` and `T`, and **there is no theory for the
viscosity of water** (the main source of liquid-viscosity knowledge is experiment
— [W] p.104). Do not try to compute it from a formula; read it from a table:

- The citable table:
  [welty_transport.md](../../../knowledge/source/books/welty_transport.md)
  (Welty 5th ed., Appendix I). It supports our `η = 0.851 mPa·s @300 K` to
  **+1.03%**.
- ⚠️ **Interpolate the table log-linearly.** Linear interpolation at 300 K is off
  by **+2.91%** (the curvature of `η(T)` is already significant over a 20 K
  interval). **The interpolation method changes the answer.**
- ⚠️⚠️ **Water's viscosity is `2.06 %/K` sensitive to temperature.** The
  `T = 300 K` used across our cases is **an inherited choice** (the originating
  sketch had no temperature) and it is recorded as tier 1. If the truth is 298 K
  that is −4.4%, at 293.15 K it is −15.0%, and **`τ_B` follows directly.** If a new
  case's sketch states a temperature, **re-read `η` — do not inherit.**
- ⚠️ **Check that `kT` and `η` use the same `T`.** If they disagree, `τ_B` is
  wrong twice over.

---

## 3. A dimensionless group is a ratio of two scales

Do not write only the number — write **which two scales it is the ratio of.**
The physical interpretation follows.

| Group | Ratio | Expression | Meaning |
|---|---|---|---|
| `φ` | occupied / total volume | — | packing |
| `Pe` | `τ_B / τ_v` | `v₀d/D_t = f_a d/k_BT` | advection vs diffusion |
| `D_r*` | `τ_B / τ_r` | `D_r τ_B` (sphere → 3) | rotational vs translational diffusion |
| `ℓ_p/d` | `ℓ_p / d` | `Pe / D_r*` | persistence length, in particle diameters |
| `ε*` | `ε / k_BT` | — | binding vs thermal fluctuation |
| `k*` | trap stiffness | `k d² / k_BT` | trap vs thermal fluctuation |
| `dt*` | `dt / τ_B` | — | numerical resolution |
| `L/d` | — | — | box, in particle diameters |
| `T*` | `T_obs / τ_B` | — | observation window |
| `St` | `τ_p / τ_B` | `m/(γτ_B)` | inertia vs diffusion |
| `Re` | — | `ρ_f v₀ d / η` | fluid inertia |

### 3.1 ⭐️ `Pe` names two different numbers — and the shear one has `|E|` on top

The `Pe = v₀d/D_t` in the table above is the **self-propulsion (advective)
Peclet**, and the `Pe` of rheology is a different number. **They share a name, so
always record which one the spec means** — storing a dimensionless group as a bare
number means nobody can catch "the name is X but the value is something else."

| Name | Definition | Where |
|---|---|---|
| advective Peclet | `v₀d/D_t` | self-propulsion, dragged driving |
| **orientational Peclet = Weissenberg** | **`|E| / D_r`** | how much the flow distorts the orientation distribution ([L] p.106) |

⚠️⚠️ **`|E| = √(2 E:E)`, not `|∇u|`** (execution-verified —
`verify/verify_book_claims.py`):

- **Vorticity merely rotates the particle uniformly and cannot change an
  isotropic equilibrium distribution.** Only the rate-of-strain component sets
  the orientation ([L] p.106). **A pure rotational flow has `|∇u| ≠ 0` but
  `Pe = 0`.**
- Using `√(E:E)` is wrong by **√2 = 1.414×**. The convention that makes
  `√(2E:E) = γ̇` in simple shear is the correct one.
- For a **near-sphere (`r−1 ≪ 1`)** the correct scale is **`|E|(r−1)/D_r`** —
  writing plain `|E|/D_r` overestimates it.

**`Wi` versus `De`** ([L] (1.1)): the same ratio `τ_relax/τ_flow`, but it is
called **`Wi` for steady flow and `De` for unsteady (oscillatory) flow.** SAOS
needs **both** — `Wi` for whether the amplitude is inside linear response, `De`
for whether the frequency is near the relaxation time ([L] p.181).

**The orientational relaxation time is `1/(6 D_r)`** — not `1/D_r` and not
`1/(2D_r)`. The loss modulus of a dilute axisymmetric suspension goes as
`De/(36+De²)`, so the pole is at `De = 6`, and it is **exactly a single Maxwell
mode with `λ = 1/(6D_r)`** (deviation <1e-12, execution-verified).
⚠️ **The 6 is the 3D `l=2` (stress) mode coefficient**, so pasting it onto a
first moment or onto a 2D angular correlation (`abp-rod`'s `⟨cos Δθ⟩`) is wrong.

→ Details plus the stress tensor, the Kramers form, and semi-dilute rods:
[leal_microstructural_rheology.md](../../../knowledge/source/books/leal_microstructural_rheology.md)

---

## 4. Scale-separation checks — all against `10⁻²`

There are **two kinds** of check. Do not conflate them.

- **Model validity** — "is BD (overdamped) the right model for this system?"
  Independent of dt.
- **Integration resolution** — "is the chosen model being integrated finely
  enough?" *The integration step must be at most 1% of the fastest physical
  timescale.*

| Kind | Check | Condition | On violation |
|---|---|---|---|
| **model** | inertia negligible | `τ_p / τ_dyn ≤ 10⁻²` | ❌ overdamped BD invalid → Langevin needed |
| **integration** | advection resolved | `dt / τ_v ≤ 10⁻²` | ❌ moves more than 1% of a diameter per step |
| **integration** | rotation resolved | `dt / τ_r ≤ 10⁻²` | ❌ orientational dynamics collapses |
| **integration** | interaction resolved | `dt / τ_int ≤ 10⁻²` | ❌ inaccurate force integration |
| **integration** | trap resolved | `dt / τ_k ≤ 10⁻²` | ❌ trap relaxation unresolved |
| **model** | low Reynolds | `Re ≤ 10⁻³` | ❌ fluid inertia cannot be neglected |
| **geometry** | cutoff | `r_c < L/2` | ❌ minimum image violated |
| **geometry** | finite size | `ℓ_p, ξ ≤ L/4` | ⚠️ artefacts |
| **statistics** | enough observation | `T_obs ≥ 10² · max(τ)` | ⚠️ insufficient statistics |
| — | thin margin | a hard check is within 1/5 of its limit | ⚠️ no room to raise a parameter |

`τ_dyn` = **the fastest physical timescale of interest** in this system (`τ_k` for
a trap, `τ_r` for an ABP, and so on).

> ⚠️ **Do not compare `τ_p` against `dt`.** BD has no inertia at all — HOOMD's
> Brownian is `dr/dt = (F_C + F_R)/γ`, and **mass does not appear in the
> equation.** There is no inertia left to resolve, so `dt ≫ τ_p` is not required.
> Demonstrated: `verify/golden_trap.py` violates that condition by five orders of
> magnitude (`τ_p/dt = 4000`, mass 1, γ=1, dt=2.5e-4) and the result was accurate
> to 0.38%.
>
> The question `τ_p` answers is **"is BD an admissible model for this system?"**
> If `τ_p ≈ τ_dyn`, the real particle moves ballistically and BD cannot produce
> that, and **no choice of dt fixes it.** Switch to Langevin.

**Always report the margin alongside.** Pass/fail alone does not tell you "if I
double Pe, which check binds first?"

> ⚠️ Passing the upper bound **does not guarantee accuracy.** Brownian carries
> O(δt) error, so a publication-grade result needs a separate `dt`-halved
> convergence check.

⚠️ **And a check that is never wired up cannot be wrong out loud.** In this
project the step-resolution check silently never ran across 81 runs because of a
name mismatch, and the pre-run gate rejected 80 of 83 specs with zero real
failures — see [docs/02](../../../docs/02-verification.md#6--the-failure-mode-this-document-exists-to-prevent).

---

## 5. Inversion

The simulation is dimensionless; comparison against literature and experiment is
in physical units. **Store both.**

```python
D_eff = D_eff_star * sigma**2 / tau_B     # → µm²/s
x2    = x2_star * sigma**2                # → µm²
t     = step * dt_star * tau_B            # → s
P     = P_star * kT / sigma**dim          # → Pa
```

---

## 5.1 Observable-extraction traps (found by measurement)

### ★ Do not subtract the sample mean from a displacement autocorrelation

The **displacement** from an anchor or a trap centre has a true mean of exactly
zero. Doing the habitual `x -= x.mean()` makes `C(t)` decay **systematically
faster** by `O(τ/T_obs)`, and the `τ` extracted from it comes out small. This is
the silently-wrong kind.

Measured (`cases/trap_2d_5um.py`, `T_obs = 120 τ_k`):

| Treatment | error in the fitted τ |
|---|---|
| sample mean subtracted | **−7.75%** ✗ |
| not subtracted | **−0.26%** ✓ |

```python
x = trace.astype(np.float64)          # <- no mean subtraction (the true mean is 0)
F = np.fft.rfft(x, n=nfft, axis=0)
ac = np.fft.irfft(F * np.conj(F), n=nfft, axis=0)[:n_t]
ac /= np.arange(n_t, 0, -1)[:, None, None]     # unbiased (divide by the overlapping-pair count)
```

For a quantity whose true mean is unknown (an absolute position, say) you must
subtract it. **A displacement is not such a quantity.**

### Extract `τ` with an exponential curve_fit, not a log-linear one

Fitting a straight line to `log C(t)` over-weights the noise in the tail.
Fit `A·exp(−t/τ)` directly over `[0, ~3τ]`.

### Error bars on `⟨x²⟩` come from block averaging

Time-series samples are correlated within `2τ`. A naive SEM underestimates the
error. Split the samples into ~20 blocks and use the standard error of the block
means.

⚠️ **And block SEM itself underestimates when the system produces discrete
stochastic events.** Measured across a velocity sweep, the block SEM was
**1.09–2.28×** too small depending on the observable and the velocity — and two
conclusions were reversed once a 9-seed ensemble replaced single runs. Sizing the
ensemble to the events, not to a policy minimum, is the fix.

### PSD normalization

Use the one-sided density so that `∫₀^∞ S(f) df = ⟨x²⟩`
(`scipy.signal.welch(..., 'density')`). The analytic solution for an OU process:
```
S(f) = 4⟨x²⟩τ / (1 + (f/f_c)²),     f_c = 1/(2πτ) = k/(2πγ)
S(0) = 4⟨x²⟩τ = 4γk_BT/k²
```

---

## 6. Per-case scale tables (verified ones only)

### 6.1 A single particle in a harmonic trap — `trap-*` ✅ verified

The simplest case, and **it has an analytic solution**, so it serves as the
golden test.

**Scales**
| Kind | Name | Definition |
|---|---|---|
| length | `d` | particle diameter (the reference) |
| | `ℓ_k = √(k_BT/k)` | thermal fluctuation width inside the trap |
| | `L` | box |
| time | `τ_p = m/γ` | inertia (to be discarded) |
| | `dt` | integration step |
| | `τ_k = γ/k` | **trap relaxation — the governing timescale of this system** |
| | `τ_B = d²/D_t` | diffusion |
| | `T_obs` | observation window |
| energy | `k_BT` | the reference |
| | `k d²` | trap stiffness (the work to pull by one diameter) |

**Dimensionless groups**: `k* = k d²/k_BT` · `ℓ_k/d = 1/√k*` · `dt/τ_k` ·
`τ_k/τ_B = 1/k*`

**Checks**: integration `dt/τ_k ≤ 1e-2` · model `τ_p/τ_k ≤ 1e-2` · geometry
`ℓ_k < L/2` (the weaker the trap, the more vulnerable to trap 1)

**Analytic solution (the golden test)**
```
per degree of freedom  ⟨x²⟩ = k_BT / k       relaxation time  τ = γ/k
2D total               ⟨r²⟩ = 2 k_BT / k     distribution     P(x) ∝ exp(-k x²/2k_BT)
```

**Measured verification** (`verify/golden_trap.py`, N=400, 340k steps):
| k | ⟨x²⟩ measured | predicted | error |
|---|---|---|---|
| 2 | 0.50188 ± 0.0016 | 0.5 | +0.38% |
| 5 | 0.20112 ± 0.0006 | 0.2 | +0.56% |
| 10 | 0.09998 ± 0.0003 | 0.1 | −0.02% |
| 20 | 0.04993 ± 0.0001 | 0.05 | −0.14% |

`⟨x²⟩·k` stays at 1.0 across a 10× range in k (coefficient of variation
**0.28%**).

**End-to-end result** (`cases/trap_2d_5um.py`, N=1000, 2000 τ_k, 1M steps, about
3 minutes). Physical-system YAML (SI) → scale table → separation checks →
non-dimensionalization → HOOMD run → inversion to physical units:

| Observable | Measured | Analytic | Error |
|---|---|---|---|
| `⟨x²⟩` | 4.14293e-4 µm² | 4.14195e-4 µm² | **+0.02%** |
| `σ = √⟨x²⟩` | 20.3542 nm | 20.3518 nm | **+0.01%** |
| `τ` (from the C(t) fit) | 4.00022 ms | 4.01024 ms | **−0.25%** |
| `f_c` (from the PSD fit) | 40.1516 Hz | 39.6871 Hz | **+1.17%** |
| `S(0)` (from the PSD fit) | 6.62187e-6 µm²/Hz | 6.64409e-6 µm²/Hz | **−0.33%** |

Statistical error ±0.091%, expected systematic bias +0.100% (`dt/τ_k = 2e-3`) —
the measured +0.02% is within 1 SEM.

Confirmed from the plots:
- `P(x)` is Gaussian **over six decades (±5σ)**
- `C(t)` is a single exponential over more than two decades. Slightly below at
  long lag — tail noise of the unbiased estimator
- The PSD rises above the Lorentzian for `f > ~8 f_c` — **aliasing from discrete
  sampling** (expected; components above Nyquist fold back in). No effect on the
  `f_c` fit

**What this case taught**
- A weak trap (small `k*`) is the most vulnerable to trap 1 (minimum image).
  **Verify in the weak condition.**
- `dt = τ_k/2000`, equilibration `20τ_k`, production `150τ_k` gives errors within
  0.6%.
- ⭐️ **Do not conflate the three decisions** (they were confused here and then
  sorted out — see section 1.1): the unit system stays **fixed** on `τ_B`, `dt` is
  **set** by `τ_k`, and the reporting is **rescaled** into units of `τ_k`.
- **Some systems make `τ_B` meaningless.** Here `τ_B = 242 s` while the trap
  catches the particle in `4 ms`, so free diffusion is never realized. Keep `τ_B`
  as the unit-system reference (fixed), but make every physical judgment in
  `τ_k`.
- Do not subtract the sample mean from a displacement autocorrelation →
  section 5.1
- ⭐️ **This case is where we found that the `τ_p` check must compare against
  `τ_dyn`, not against `dt`.** In the real case (d=5µm in water, k=10pN/µm),
  `τ_p=3.26µs` and the recommended `dt=2.0µs`, so **no dt existed at all** that
  satisfied the old criterion `τ_p/dt ≤ 1e-2`
  (`100τ_p ≤ dt ≤ 0.01τ_k` ⟹ `τ_p/τ_k ≤ 1e-4`, whereas the truth is 8.1e-4).
  Under the correct criterion `τ_p/τ_k = 8.1e-4 ≤ 1e-2` it passes with 12× margin.

**The real physical numbers** (d=5µm silica, water@300K, k=10pN/µm —
`verify/trap_scales.py`):
| Quantity | Value | Note |
|---|---|---|
| `γ` | 4.01e-8 kg/s | |
| `D_t` | 0.103 µm²/s | |
| `τ_p` | 3.26 µs | for the model check |
| **`τ_k`** | **4.01 ms** | **the governing timescale** |
| `τ_B` | 242 s | never realized, because of the trap |
| `ℓ_k` | 20.35 nm | **independent of `d`** |
| `k*` | 6.04e4 | very stiff |

`ℓ_k/d = 4.1e-3` — the particle moves only 1/250 of its own size. Normal for
optical tweezers.

### 6.2 Soft repulsive pair plus an excluded-volume core — `soft-r3-*` ✅ verified

`U/kT = A(d/r)³ + WCA(σ=d, ε=kT)`, 2D periodic, fixed φ. No analytic solution.

**Start with what was expected and wrong**: the general form `τ_int = d²γ/ε`
**cannot be used.** A soft potential has a different stiffness at every `r`, so
`τ_int` is not a constant.

**Scales**
| Kind | Name | Definition | Note |
|---|---|---|---|
| length | `d` | particle diameter (the reference) | |
| | `r_min` | **closest approach distance** ★ | dt comes from here |
| | `a_mean = ρ^(-1/2)` | mean spacing | `φ = (π/4)(d/a)²` |
| | `a_NN = √(2/√3)·a_mean` | **hexagonal nearest-neighbour distance** | = 1.07457 a_mean. Not `a_mean` |
| | `r_c` | cutoff | set it as **a multiple of a_mean** (below) |
| | `L = a_mean√N` | box | |
| time | `τ_p = m/γ` | inertia (model check) | |
| | `τ_int(r) = γ/U''(r)` | **local interaction relaxation** ★ | same structure as a trap's `τ_k=γ/k` |
| | `τ_B = d²/D_t` | diffusion = **the governing timescale** | unlike case 6.1, here it really does govern |
| energy | `k_BT` | the reference | |
| | `Γ = U(a_mean)/kT` | **coupling strength** ★★ | the real control parameter |
| | `A` | the `r⁻³` amplitude | the value written on the sketch |
| | `U(d) = (A+ε)kT` | contact binding | ⚠️ **not `A·kT`** (below) |

> ⚠️ **`A ≠ U(d)/kT`.** At `r=d` the WCA core contributes exactly `ε`
> (`4ε(1−1)+ε`), so `U(d) = (A+ε)kT`. With `A=100` that is 101 — off by 1%. The
> ledger had this entry written as "= A kT" and the L3 integrity check caught it
> (section 6.4). `Γ` is unaffected: `a_mean = 1.498d > 2^(1/6)d` is outside the
> WCA cutoff, so `Γ = A(d/a)³` is exact.

**Dimensionless groups**: `Γ = A(d/a_mean)³` · `φ` · `a_mean/d` · `L/d` ·
`r_c/a_mean` · `dt/τ_int(r_min)`

> ⭐️ **`A` alone is not the control parameter.** What sets the structure is
> `Γ = A(d/a_mean)³`. The same `A=100` gives Γ=57.9 at `a=1.2d` and Γ=1.26 at
> `a=4.3d` — crystal versus fluid. Without a density, an amplitude sweep is
> **undefined.**

#### How to set `r_min` (dt hangs on it)

Take the **smaller** of two criteria.

```
(a) pair criterion    the r where U(r) = u_max·kT.  u_max = ln(number of pair samples) ~ 12  (10^6 samples)
(b) vibration crit.   a_NN - 3σ_bond ,  σ_bond = √2·√(kT/k_cage) ,  k_cage = 3[U''(a_NN)+U'(a_NN)/a_NN]
                      * if the Lindemann ratio σ_bond/a_NN > 0.15 the cage has melted, so do not use (b)
```
Then `dt = 10⁻²·τ_int(r_min)`.

**This was wrong three times and fixed by measurement** (do not repeat them):
1. Setting `u_max = 5` put the measured minimum distance **4.3σ inside** the
   prediction. With 10⁶ pair samples the Boltzmann tail reaches
   `βU ≈ ln(10⁶) ≈ 14` → use `u_max = 12`.
2. The cage was evaluated at `a_mean`. **The hexagonal nearest-neighbour distance
   is `a_NN = 1.0746 a_mean`**, and since `U'' ∝ r⁻⁵` that is a 41% difference.
3. `σ_bond` was written as `√(kT/k_cage)`. That is the **per-component** rms
   (`u₁`); the bond-length fluctuation is `√2·u₁`. Both the Lindemann verdict and
   `r_min` come out wrong.
   (Errors 1 and 3 partially cancelled, so `r_min` happened to look close — do
   not rely on a coincidence.)

⚠️ **In weak coupling, dt is set by the WCA core, not by `r⁻³`.** Measured: Γ=0.03
was **12.8× more expensive** than Γ=30 (`dt/τ_B` 2.83e-6 vs 3.62e-5). The
physically dullest point is the numerically hardest. Strong coupling is caged, so
close approach is itself restricted.

#### The cutoff — an absolute kT criterion does not work ⭐️

In 2D the `r⁻³` tail energy falls off only as `1/r_c`:
energy outside `r_c` / nearest-neighbour energy = `2πa/(3r_c)` → still **42%** at
`r_c=5a`. Setting it by an absolute criterion (`U(r_c) ≤ 0.01 kT`) collides with
minimum image:

```
solving r_c = (A/u_c)^(1/3) d  together with  r_c < L/2 = (√N/2)a  gives
      Γ_max = N^(3/2)·u_c / 8        <- independent of A
      N=100, u_c=0.01 -> Γ_max = 1.25 (weakly correlated fluid only)
```
→ **Set the cutoff as a multiple of `a_mean` and check convergence.** Measured
(`r_c` 5a→7a, Γ=29.7):

| Quantity | Change |
|---|---|
| ψ₆ · NN distance · 6-coordinated fraction | **unchanged within 0.15%** |
| absolute `⟨U⟩/N` | **+7.5%** |

**The structure converges and the thermodynamics does not.** Forces from the far
field cancel by symmetry, but the energy accumulates. If pressure or energy is
being reported, a tail correction is needed.

#### Verification strategy when there is no analytic solution (the biggest
difference from case 6.1)

| Check | Result | What it catches |
|---|---|---|
| direct two-particle comparison (`verify/verify_pair_table.py`) | **0.000%** | the potential / table implementation |
| energy consistency `⟨U⟩/N` vs `(ρ/2)∫U g(r) 2πr dr` | **+0.00 to +0.67%** (7 runs) | consistency of the measured g(r) with the force sum |
| hexagonal NN distance vs `√(2/√3)·a_mean` | **+0.45%** | the crystal structure (a parameter-free prediction) |
| dilute-limit `g(r)` vs `e^{-βU}` | RMS **6.30%** ✗ | — |
| dilute-limit `g(r)` vs `e^{-βU}[1+ρ∫f f]` | RMS **2.43%** ✓ | the equation of state |

⚠️ **Do not compare a dilute `g(r)` directly against `e^{-βU}`.** Even at
φ=0.01, A=10 the O(ρ) cluster term contributes +5 to +14%. A repulsive system has
`f<0` so `∫f f>0` → **`g > e^{-βU}`** (the sign comes out right too).
"Is the potential implemented correctly" is answered by the two-particle
comparison; "is the statistical mechanics right" needs the first-order term
included.

#### Measured results (d=5µm, water@300K, φ=0.35, N=400, `r_c`=5a, T_obs=100 τ_B)

| A | Γ | ψ₆ | 6-coord | NN/d | σ_NN/NN | State | dt/τ_B | wall clock |
|---|---|---|---|---|---|---|---|---|
| 0.1 | 0.030 | 0.435 | 0.42 | 1.681 | 0.279 | fluid | 2.83e-6 | 138 min |
| 1 | 0.298 | 0.435 | 0.43 | 1.680 | 0.259 | fluid | 3.14e-6 | 128 min |
| 10 | 2.97 | 0.458 | 0.53 | 1.661 | 0.182 | fluid | 1.37e-5 | 44 min |
| 100 | 29.7 | **0.885** | **0.987** | 1.617 | 0.066 | **hexagonal crystal** | 3.62e-5 | 20 min |

- **The crystallization transition lies between Γ 3 and 30.** (Not compared
  against the literature at the time — the KB had zero papers then. A later
  campaign did compare it, and found that a truncation error had biased an
  exponent by 2.9σ and that a run labelled "hexatic" was in fact a crystal.)
- **Γ=0.03 and Γ=0.30 are indistinguishable** (ψ₆ is 0.4347 for both). A
  four-decade sweep in A produces only three states at φ=0.35 — the weak-coupling
  end is simply a WCA fluid.
- Einstein cage prediction `σ_bond/a_NN`=0.056 vs measured 0.066 (+18%). This is
  anharmonicity, so it is used **only as a regime indicator, not as a quantitative
  prediction** (in 2D the absolute `u_rms` diverges logarithmically). ★ The same
  trap recurs: a naive harmonic prediction for bond-stretch variance in a DLVO
  well was off by **4.57×** until the whole basin was Boltzmann-integrated.

### 6.2b Hexagonal lattice plus a moving trap — `trap-drag-2d-hex300` 🔶 L3 only

`soft-r3` (6.2) with one `trap` (6.1) laid on top. What is new is that **there are
two stiffnesses.** Script: `cases/trap_drag_2d.py`

**What this system brought for the first time**

| Item | Content |
|---|---|
| two competing stiffnesses | `τ_k = γ/k_t` = 4.01 ms vs `τ_int = γ/U''(r_min)` = 876 ms → **the trap is 218× faster** |
| `dt` decision | the **faster** of the two. `dt = 10⁻²·τ_k` = 40.1 µs. Pair resolution has 218× margin |
| advective timescale | `τ_v = d/v_x` = 10 s — a moving boundary condition. A new `dt/τ_v` check appears |
| wake healing | periodic box, so the probe **returns into its own wake**. `v·τ_int/L_x` = 3.2e-3 (313× margin) |
| box | ★ **commensurate hexagonal, therefore not square** |

**⭐️ This system's real gate is not the separation checks but the statistics.**
All seven hard checks pass, and running it as planned does not produce the value
you want:

```
Δr_ss = γv/k_t = 2.005 nm     <- signal (the probe's bare Stokes lag)
ℓ_k   = √(kT/k_t) = 20.35 nm  <- noise (thermal fluctuation inside the trap)
SNR = 0.0985                   the signal is 1/10 of the noise
```
One traverse (`T_obs = L_x/v_x` = 274 s) gives `T_obs/2τ_k` = 34,119 independent
samples, so the precision on `⟨F_drag⟩` is `1/(SNR√n)` = **5.5%** against a 2%
target. Reaching 2% takes **8 traverses.** And the lattice period recurs only
**17 times**, so if there is stick-slip modulation that binds first.
→ The design intent is to use **the strain field of all 306 lattice sites** as the
primary observable (N× the statistics).

> ⚠️ **`Δr_ss = γv/k_t` is the probe's bare Stokes lag** — the value you get with
> no lattice at all. The microrheological signal is the **excess on top of it**,
> and there is no prediction for it (`measurement`). Do not lose track of what the
> numerator is when calling something "signal to noise".

> ⭐️ **The trap overwhelms the pair interaction** (`k*` = 6.04e4 vs `Γ` = 29.7).
> The probe is effectively a **constant-velocity boundary condition** pinned by the
> trap, not a probe that "applies a force and watches the response."

> ⚠️ **Two of this case's published conclusions were later reversed**, and not by
> a code change: a 7-velocity × 9-seed re-run showed the defect count is *not*
> v-independent (it peaks non-monotonically) and the recovery fraction is not
> monotonic. The cause was a single run plus a block SEM that underestimated the
> spread by 1.09–2.28×. See [docs/04](../../../docs/04-cases.md#trap-drag-2d-hex300--single-run-error-bars-were-wrong-in-both-directions).

#### ⭐️ A commensurate hexagonal lattice — the box aspect ratio is not chosen, **commensurability sets it**

The sketch says "start from hexagonal equilibrium", so the initial lattice must be
commensurate with the periodic box. Incommensurate, and defects are injected at
the seam — and this system's observables (the lattice strain field, ψ₆) are
precisely sensitive to that.

```
L_x = n_x · a_NN              (row direction = nearest-neighbour direction)
L_y = n_y · (√3/2) · a_NN     (row spacing)
N   = n_x · n_y ,  n_y even    <- the staggered rows have a 2-row period
φ is preserved identically:  φ = πd²/(4 a_mean²),  L_x L_y = N a_mean²
```

**`L = a_mean√N` (square) is not commensurate.** φ comes out right but the lattice
does not match at the seam.

| Candidate | N | ΔN | aspect ratio | lattice periods n_x |
|---|---|---|---|---|
| **17×18** ← adopted | 306 | +2.0% | 1.091 | **17** |
| 16×18 | 288 | −4.0% | 1.026 | 16 |
| 15×20 | 300 | 0% | 0.866 | 15 |

⚠️ **Drag along the nearest-neighbour (row) direction.** Then the potential's
spatial period is `a_NN` and the number of periods per traverse is `n_x`. Drag
**perpendicular** to the rows and the period becomes 2 rows = `√3 a_NN`, so the
same box gives *fewer* periods (17 → 10).

> ⚠️ **A check that recomputes φ from the lattice and compares does not run** —
> `L_x` and `L_y` are derived from φ, so it passes identically. The comparison that
> actually means something is **against the `box_length_*` recorded by L2** (the
> same attitude as recomputing derived values). The former was written first and
> an adversarial check caught it.

### 6.2c Three-point bending chain plus oscillatory driving — `chain-bend-2d-oscill` 🔶 L3 only

`G'(ω)` and `G''(ω)` of a colloidal bead chain. The physics came from the
literature ([P1] Pantina & Furst PRL 94 138301 / [P2] Langmuir 24 1141).
Script: `cases/chain_bend_2d.py`

**★ The bond is not a pair potential** — adhesive contact (JKR) plus tangential
bending stiffness:
```
EI = κ₀ a³/3        κ_θ = EI/ℓ  (ℓ = d)          <- the k of angle.Harmonic
κ_end    = 24 EI/L³   <- the paper's definition (end-force basis). The value to compare against
κ_center = 48 EI/L³   <- what the driving trap feels
δ_max = M_c L²/(12 EI)  <- the linear-elastic amplitude limit, M<M_c
```
⚠️ **Mixing the force definitions puts you off the published value by exactly
2×.**
⚠️ The discrete-to-continuum mapping error **depends on n**: −0.35% at n=25,
−0.08% at n=51.

**⭐️ The mode that sets `dt` is not the one being observed** — a first in this
project.

| | Value | Role |
|---|---|---|
| `τ_fast = γ/λ_max` | 0.279 µs | **sets dt.** Largest eigenvalue of the stiffness matrix |
| `τ_chain = γ/κ_center` | 1.27 ms | **what we want to measure.** The G'(ω) band |
| ratio | **2.2e-4** | 4570× — the cost is fixed independently of the band of interest |

→ At the lowest ω (85 rad/s, 100 cycles) that is **2.65e9 steps**. The cause is
that `κ_θ = 1.39e6 kT`, so the chain is thermally completely rigid. Options:
(a) accept it (b) lower κ₀ (surfactant, [P2] Fig. 4) (c) measure only high ω and
replace low ω with the quasi-static limit. Cost ∝ 1/ω, so the low-ω end dominates.

⚠️ **And the governing timescale here was itself derived wrong at first.** Using
`τ_chain = γ/κ_center` rather than `τ_max = γ/λ_min` was off by **9.18×**, which
propagated into De, the equilibration length and the SNR. Setting `dt` from
`λ_max` **without also looking at `λ_min`** is the trap — both come free from the
same eigendecomposition, and it is `λ_min` that governs.

**`λ_max` is the larger of the two blocks — do not add them.** In a straight chain
the stretching (x) and bending (y) blocks **decouple** at linear order. Adding
them overestimates `λ_max` by 18% and only shrinks dt.
```
bending    4.231e-2 N/m  (= 15.86 κ_θ/ℓ², n=25)   * this one is faster
stretching 7.639e-3 N/m  (= 4 k_b)
```

**⭐️ The fastest mode is not overdamped** — `τ_p/τ_fast = 0.60`, violating the
1e-2 criterion by 60×. BD treats that mode as overdamped, so **the dynamics in
that band is wrong and no choice of dt fixes it** (section 4 — BD has no inertia
at all). It sits 4570× away from the band of interest (`τ_chain`), so it likely
does not affect `G'(ω)`.

✅ **And that was then measured**: comparing `OverdampedViscous` against
`Langevin(kT=0)` at all 7 frequencies bounded the difference at **0.159%**. BD is
admissible here. ★ The thermal comparison (Brownian vs Langevin at kT=1)
**lacked the power to exclude even a 47% effect**, because `|ŷ|/ℓ_k < 1` — which
is why an integrator assumption must be tested with a `kT=0` deterministic
difference.

The model check used for the verdict is against the **fastest scale of
interest**: `τ_p/τ_chain` = 1.3e-4 (76× margin).

**The amplitude is squeezed from both sides** (a shape no other case had):
```
ℓ_k (20.4 nm)  ≪  a (200 nm)  <  δ_max (429 nm)
   below: buried in thermal noise      above: the bond slips and a harmonic angle is invalid ([P2])
   SNR = 9.83                          linear margin 2.1x <- the tightest hard check in this system
```

**There are no periodic boundaries** → `declare_absent("box", reason)`. What plays
the role of the geometric limit is not the box but `δ_max`. A scale that does not
exist must be left empty, not invented (rule 3).

### 6.3 Comparing two cases — what really appeared twice

| | `trap-2d-5um` | `soft-r3-2d` | Common? |
|---|---|---|---|
| reference scales | `d`, `kT`, `τ_B` (thermal) | same | ✅ **appeared twice** |
| derived properties | `γ=3πηd`, `D_t=kT/γ` | same | ✅ **twice** |
| use of `τ_p` | model validity only (never compared to dt) | same | ✅ **twice** |
| `dt` decision | `10⁻²·τ_k`, `τ_k=γ/k` | `10⁻²·τ_int(r_min)`, `τ_int=γ/U''` | ✅ **twice** — both are `γ/(local stiffness)` |
| check categories | model / integration / geometry / statistics | same | ✅ **twice** |
| hard vs soft | all passed, so no distinction needed | a statistics check binds → distinction needed | ⭐️ discovered in the second case |
| role of `τ_B` | only the reference; `τ_k` governs (4ms vs 242s) | reference **and** governing | same `d`, same `τ_B`, opposite meanings |
| equilibrium indicator | anchor displacement `⟨Δr²⟩` (bound system) | `⟨U⟩/N` (diffusive system) | ❌ **differs per case** |
| verification basis | 4 analytic solutions | 5 limits, consistency and convergence checks | ❌ **differs** |
| governing length | `ℓ_k=√(kT/k)`, independent of `d` | `a_mean`, tied to `d` | ❌ differs |
| geometry check | always passes (`ℓ_k`=20nm) | **a real gate** (it sets the Γ ceiling) | ❌ differs |

**✅ What was promoted into `bdbot/`** (only what appeared twice):

| Module | Content |
|---|---|
| `units` | one pint registry (mixing registries makes pint refuse) |
| `provenance` | `Provenanced` (value + source + tier) |
| `materials` | `γ=3πηd` · `D_t=kT/γ` · `τ_B=d²/D_t` · `m` · `τ_p` |
| `scales` | `ScaleLedger` (length/time/energy + reference + basis) · `thermal_reference` |
| `checks` | `Check(kind,name,value,limit,op,note,hard)` · `verdict` · **`dt = 10⁻²·γ/(local stiffness)`** |
| `report` | the `DimensionlessReport` renderer (the case supplies the input, derived and run-plan blocks) |
| `runid` | content addressing plus re-run prevention |
| `metrics` | the `metrics.json` schema (the sole input to the post-mortem) |
| `stats` | block averaging · autocorrelation correction · unbiased autocorrelation |
| `sim` | 2D frame · BD integrator · GSD · seeds · minimum image · WCA |

⭐️ **The most important finding**: the two cases' `dt` decisions were **the same
formula** — `γ/k` for the trap, `γ/U''(r_min)` for the soft pair. Both are
`γ/(local stiffness)` with `dt = 10⁻²·τ`. Only the names differed (`τ_k`,
`τ_int`).

**What was *not* promoted** (it differed per case): the equilibrium indicator (the
case declares it) · observables · verification strategy · choice of the governing
timescale · initial placement · the sampling loop.

**Equivalence check** (`verify/verify_1c_equivalence.py`): all 8 `run_id`s
preserved, and re-running the first case gave 77 identical metric fields with
zero violations. Always run this after a refactor.

---

## 6.4 The L3 artefact `NondimSpec` — what made non-dimensionalization **checkable** ⭐️

`bdbot/nondim.py`. The **only contract** between L2 (SI) and L4 (execution),
stored as a single `specs/<run_id>.json`. Three cases had each been hand-building
their own `spec = {...}` with mutually different schemas (only three keys in
common). Four measured defects:

| Defect | Symptom | Fix |
|---|---|---|
| ⭐️ `run_id` did not cover the physical system | one spec contained **no physical system at all**, so changing `d` 5µm→0.5µm and `η` by 62× (a 16.1× change in τ_B) left the run_id identical → mistaken for a completed run of a different system, and **an old result was reported as the new system's result** | hash target = `{system(physics_only), params, numerics}` |
| the spec alone could not be inverted | no SI values for σ, τ, kT (violating section 5) | put the three anchors into `back_transform` as SI floats |
| dimensionless groups could not be checked | nobody could catch "the name is `dt/τ_int` but the value is something else" | `Group(num=(cat,symbol), den=(cat,symbol))` → **recomputed and compared** against the ledger |
| a scale missing from the ledger | one case put `dt` and `T_obs` outside the ledger, so they were invisible in the timescale ordering | enforce four required **roles** (`box`, `dt`, `observation`, `inertia`) |

**Separate the two layers** — the names are similar enough to confuse:

| Layer | What it asks | Where |
|---|---|---|
| `checks` (section 4) | is BD valid for this system, and is it integrated finely enough — **physics** | `verdict()` |
| `validate()` | was the non-dimensionalization done correctly — ledger completeness, ratio consistency, invertibility | `NondimSpec.validate()` |

**The ledger separates the symbol from the description.** The old keys were a
single blob like `"d        particle diameter (reference)"`, which a machine
cannot use. Now it is `lg.add_time("dt", dt, "integration step", role="dt")`,
accessed via `lg.get("times","dt")`. Required entries are enforced **by role, not
by symbol** — for an ellipsoid case the reference length is `d_eq`, not `d`.

A scale that does not exist is emptied **explicitly** with
`declare_absent(role, reason)` (the same attitude as rule 3).

**How self-sufficiency is decided**: `bdbot.cli nondim show <run_id>` redraws the
report from **the spec alone.** If it all draws, the spec is self-sufficient; if
not, something is missing from it. L4 uses only this path (it does not import case
scripts).

**Adversarial checks** `verify/verify_nondim_guards.py` (33/33) — deliberately
broken to test: shift a ratio by ×1.41 and by ×(1+1e-6), remove a required role,
reference a nonexistent symbol, change the physical system, hand-edit the spec,
set `dt*`=0, use a dimensionally inconsistent ratio. Floating-point error (1e-12)
must pass (to prevent false positives). The regression test is
`verify/verify_l3_spec_gaps.py`.

⚠️ **The `dimensionless` keys in `metrics.json` are symbols** (`"k*"`,
`"dt/tau_int"`). They used to be the report's display strings, and such keys
cannot be queried in what is the post-mortem's only input. Use `groups_dict()`
for humans and `metrics_dict()` for machines.

---

## 7. Starting a new case — the procedure

1. Collect the **dimensional values** from the sketch or the literature. Absent
   means `null`. Do not invent.
2. Look at the existing cases in section 6 and write down **which scales this
   system adds.**
3. Draft the scale table and **confirm it with a human.** ← do not skip
4. Choose the references and **record the basis.** The governing timescale may not
   be `τ_B` (see section 6.1).
5. Non-dimensionalize → separation checks (with margins) → run → invert.
6. Compare against an analytic solution or a published value if one exists.
   Otherwise verify against a known result in a limit (dilute, low Pe, …).
7. **Add what you learned to section 6.** Especially "the scale unique to this
   system" and "the condition that failed."

---

## 7.5 Separate the two kinds of verification ⭐️⭐️

**The reason to compute these systems is that they may differ from the standard
picture.** So writing the verdict as "differs from the prediction, therefore
FAIL" means **calling a discovery a failure.**

| Role | Where the prediction comes from | A mismatch means |
|---|---|---|
| `implementation_check` | **analytically derived from the model I implemented** | **a bug** → FAIL |
| `hypothesis` | an assumption the simulation does **not** impose (continuum, dilute limit, effective medium, a literature model) | **a result** — not a FAIL |
| `measurement` | no prediction | — |

Declare it with `bdbot.metrics.observable(..., role=...)`, and
`bdbot.metrics.judge()` judges each role differently. The default is the most
conservative, `measurement`.

### A real case — `abp-rod` was a validation case, not a discovery case

All five of its predictions — `1/τ_eff = D_r + 1/τ_tumble`,
`D_eff = D̄ + v²τ_eff/2` and so on — **followed from the model I had put in.**
Agreement to within 0.66% means the code is right and tests nothing about the
physics. It is very nearly circular.

**When designing a case, write down two separate lists:**
```
assumptions the simulation imposes   overdamped BD . isotropic translational friction .
                                     independent rotational diffusion . Poisson tumbling . dilute
assumptions the theory adds          (empty)   <- if this is empty there is no discovery
```
Making `abp-rod` a discovery case requires adding something the simulation does
not impose: finite density (interactions), or confinement (walls, channels).
**Translational anisotropy is impossible in BD** (no HI).

By contrast, among what has already run, the `hypothesis`/`measurement` items:
- `soft-r3`'s **melting transition at Γ 3–30** — there was no prediction. The
  simulation is the answer.
- `soft-r3`'s energy consistency and hexagonal NN distance — those were
  `implementation_check`.
- `trap-2d-5um`'s `⟨x²⟩=kT/k` — `implementation_check` (which is what makes it a
  good golden test).
- `chain-bend-2d-dlvo`'s `K′ = 0` — a genuine `hypothesis`, derived from the
  paper's own argument rather than from our implementation, and reporting the
  agreement **was** the result.

---

## 8. Do not

- Fix dimensionless values first and infer the physical system afterwards
- Apply the spherical Stokes relation (`D_r = 3D_t/d²`) to a non-spherical particle
- Invent a material property that is not in the sketch or the literature (if you
  do not know, say so, and mark it as an anchor)
- Say accuracy is guaranteed because the separation checks passed (convergence is
  a separate check)
- Verify only in the strong condition and call it a pass (the lesson of
  section 6.1)
- Trust a check you have not wired up and watched fire — a step-resolution check
  silently never ran across 81 runs, and a pre-run gate rejected 80 of 83 specs
  with zero real failures
