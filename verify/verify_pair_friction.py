"""Measure what the three `md.pair.friction` models actually DO in this project's
systems.

Background -- `chain-bend-2d-dlvo`'s JKR control applies the bending stiffness
directly with `force.Custom`. On the question "can friction be implemented
directly", the API survey (docs/hoomd_capabilities.md, frictional contact) already
established that HOOMD 7.1 **has** `md.pair.friction` (Hofmann et al. 2025,
arXiv:2507.16388) -- but **this project had never once run it**. So it is measured
rather than ruled out by reasoning (rule 6).

Five claims to test --
  C1  the friction weight w(r) = -dU_WCA/dr is 0 outside the WCA cutoff.
      Under this case's core convention (SIGMA_CORE_STAR = 2^(-1/6), so the cutoff
      ends exactly at r*=1), the DLVO secondary minimum at r* = 1+h_min = 1.00759
      lies **outside** the cutoff, so the friction is identically 0.
  C2  (a power control) inside the cutoff it is NOT 0. If C2 came out 0, C1 would
      prove nothing.
  C3  all three models are **dissipative** -- with a tangential displacement present
      but zero relative velocity, the force is 0. So there is no elastic storage
      (K'), and no state corresponding to static friction (stick).
  C4  they are **exempt from rolling** -- if the relative surface velocity at the
      contact point is u=0, the force is 0. The bond deformation in three-point
      bending is rolling, not sliding, so this friction cannot resist bending.
  C5  `methods.Brownian` (overdamped) does not use velocity as a physical quantity
      -> so what does this force actually see under BD?

    $PY scratch/verify_pair_friction.py
"""
from __future__ import annotations

import math
import sys

import gsd.hoomd
import hoomd
import hoomd.md as md
import hoomd.md.pair.friction as FR
import numpy as np

# ── This case's actual convention (cases/chain_bend_dlvo_2d.py:598,
#    specs/...jkr...json) ──
SIGMA_CORE_STAR = 2 ** (-1.0 / 6.0)          # so the WCA cutoff ends at r*=1
                                             # (surface contact)
R_CUT_WCA = SIGMA_CORE_STAR * 2 ** (1 / 6)   # = 1.0 (exactly)
H_MIN_STAR = 0.00759259035993831             # DLVO secondary-minimum surface gap h/d
R_WELL = 1.0 + H_MIN_STAR                    # centre-to-centre r* of the secondary
                                             # minimum
R_RAD = SIGMA_CORE_STAR / 2                  # the particle radius the friction model
                                             # uses (ASSUMED = sigma/2; tested at T0)

MODELS = {
    "LJLinear": (FR.FrictionLJLinear, dict(gamma_f=1.0)),
    "LJCoulomb": (FR.FrictionLJCoulomb, dict(kappa_f=3.0)),
    "LJCoulombNewton": (FR.FrictionLJCoulombNewton, dict(gamma_f=1.0, kappa_f=3.0)),
}

results: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((name, ok, detail))
    print(f"  {'✓' if ok else '✗'} {name}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════════════════
def make_frame(sep: float, v_rel: float, omega_z: float, *, mass=1.0, diameter=None):
    """Two spheres a centre-to-centre distance `sep` apart.

    Placed along x, with relative sliding v_rel along y, both spinning about z at
    omega_z.
    """
    I = 0.4 * mass * R_RAD ** 2                     # a homogeneous sphere
    f = gsd.hoomd.Frame()
    f.particles.N = 2
    f.particles.types = ["A"]
    f.particles.typeid = [0, 0]
    f.particles.position = [[-sep / 2, 0, 0], [sep / 2, 0, 0]]
    f.particles.orientation = [(1, 0, 0, 0)] * 2
    f.particles.velocity = [[0, -v_rel / 2, 0], [0, v_rel / 2, 0]]
    f.particles.mass = [mass, mass]
    f.particles.moment_inertia = [(I, I, I)] * 2
    # With q=(1,0,0,0), angmom = 2*q(x)(0, I*omega) = (0,0,0, 2*I*omega_z)
    f.particles.angmom = [(0, 0, 0, 2 * I * omega_z)] * 2
    if diameter is not None:
        f.particles.diameter = [diameter, diameter]
    f.configuration.box = [20, 20, 20, 0, 0, 0]
    f.configuration.dimensions = 3
    return f


def probe(model: str, sep: float, v_rel: float, omega_z: float = 0.0, *,
          kT: float = 0.0, method: str = "nve", n_steps: int = 0, diameter=None):
    """Measure force and torque. Returns (F_radial, F_tangential, |torque|) for
    particle 1."""
    cls, extra = MODELS[model]
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(make_frame(sep, v_rel, omega_z, diameter=diameter))
    cell = md.nlist.Cell(buffer=0.4)
    fr = cls(nlist=cell, default_r_cut=R_CUT_WCA)
    fr.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR, kT=kT, **extra)
    if method == "nve":
        meth = md.methods.ConstantVolume(filter=hoomd.filter.All())
    else:
        meth = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-6, methods=[meth], forces=[fr])
    integ.integrate_rotational_dof = True
    sim.operations.integrator = integ
    sim.run(n_steps if n_steps else 0)
    F = np.array(fr.forces)          # (2,3)
    T = np.array(fr.torques)         # (2,3)
    vel = np.array(sim.state.get_snapshot().particles.velocity)
    return dict(F_rad=float(F[1, 0]), F_tan=float(np.hypot(F[1, 1], F[1, 2])),
                torque=float(np.linalg.norm(T[1])), F=F, T=T, vel=vel)


# ══════════════════════════════════════════════════════════════════════════
print("=" * 84)
print("T0 . convention arithmetic -- in this case, does the DLVO secondary minimum "
      "lie outside the WCA cutoff?")
print("=" * 84)
print(f"  SIGMA_CORE_STAR = 2^(-1/6) = {SIGMA_CORE_STAR:.10f}")
print(f"  WCA cutoff r_cut = sigma*2^(1/6) = {R_CUT_WCA:.15f}")
print(f"  DLVO secondary min r* = 1 + h_min = {R_WELL:.10f}   "
      f"(h_min = {H_MIN_STAR*1470:.2f} nm)")
check("r_cut == 1.0 (ends exactly at surface contact)", abs(R_CUT_WCA - 1.0) < 1e-12,
      f"r_cut − 1 = {R_CUT_WCA - 1.0:+.3e}")
check("C1 premise: the secondary minimum is outside the cutoff", R_WELL > R_CUT_WCA,
      f"r_well/r_cut = {R_WELL / R_CUT_WCA:.6f}  "
      f"(margin {(R_WELL-R_CUT_WCA)*1470:.1f} nm)")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T1 . API -- instantiate all three + their parameter keys")
print("=" * 84)
for name, (cls, extra) in MODELS.items():
    try:
        c = hoomd.md.nlist.Cell(buffer=0.4)
        o = cls(nlist=c, default_r_cut=R_CUT_WCA)
        o.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR, kT=0.0, **extra)
        keys = sorted(o.params[("A", "A")].keys())
        check(f"{name} instantiates", True, f"params = {keys}")
    except Exception as e:                                    # noqa: BLE001
        check(f"{name} instantiates", False, f"{type(e).__name__}: {e}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T2/T3 . ★ the core test -- outside the cutoff (secondary min) vs inside "
      "(control). Sliding speed V=1.0, kT=0")
print("=" * 84)
V = 1.0
print(f"  {'model':>16} | {'r*':>8} {'F_radial':>12} {'F_tangent':>12} {'|torque|':>12}")
print("  " + "-" * 70)
for name in MODELS:
    out_well = probe(name, R_WELL, V)
    out_ctrl = probe(name, 0.95, V)
    for lbl, o in (("secondary-min", out_well), ("control", out_ctrl)):
        r = R_WELL if lbl == "secondary-min" else 0.95
        print(f"  {name if lbl=='secondary-min' else '':>16} | {r:8.5f} "
              f"{o['F_rad']:12.5g} "
              f"{o['F_tan']:12.5g} {o['torque']:12.5g}   {lbl}")
    check(f"C1 {name}: friction is 0 at the secondary minimum",
          out_well["F_tan"] == 0.0 and out_well["torque"] == 0.0,
          f"F_tan={out_well['F_tan']:.3e}, |τ|={out_well['torque']:.3e}")
    check(f"C2 {name}: friction != 0 in the control (r*=0.95) [statistical power]",
          out_ctrl["F_tan"] > 1e-12 and out_ctrl["torque"] > 1e-12,
          f"F_tan={out_ctrl['F_tan']:.5g}, |τ|={out_ctrl['torque']:.5g}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T4 . C3 dissipativeness -- tangential DISPLACEMENT present, VELOCITY zero: "
      "is there elastic storage?")
print("=" * 84)
print("  Same r*=0.95; only the sliding speed changes from V=1 to 0. If it were "
      "elastic, a force would remain.")
print(f"  {'model':>16} | {'V=1 F_tan':>12} {'V=0 F_tan':>12} {'V=0 |τ|':>12}")
print("  " + "-" * 60)
for name in MODELS:
    moving = probe(name, 0.95, 1.0)
    still = probe(name, 0.95, 0.0)
    print(f"  {name:>16} | {moving['F_tan']:12.5g} {still['F_tan']:12.5g} "
          f"{still['torque']:12.5g}")
    check(f"C3 {name}: tangential force is 0 at rest (no elastic storage)",
          still["F_tan"] == 0.0 and still["torque"] == 0.0,
          f"F_tan={still['F_tan']:.3e} -- dissipative, confirmed. The normal WCA is "
          f"still alive at {still['F_rad']:.5g}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T5a . ★ what contact radius R does the model use? -- measure it, do not "
      "assume it")
print("=" * 84)
print("  |torque| = R*F_tan, so measuring the ratio gives R. **Assuming it was "
      "sigma/2 is what made T5 fail.**")
o = probe("LJLinear", 0.95, V)
R_MEAS = o["torque"] / o["F_tan"]
print(f"  |torque|/F_tan = {R_MEAS:.10f}     (had it been sigma/2 = {R_RAD:.10f}, "
      f"that is the value it should have given)")
check("★ R is **NOT** sigma/2 -- it ignores sigma", abs(R_MEAS - R_RAD) > 1e-3,
      f"R_meas − σ/2 = {R_MEAS - R_RAD:+.4e}")

o_d = probe("LJLinear", 0.95, V, diameter=2 * R_RAD)
R_d = o_d["torque"] / o_d["F_tan"]
print(f"  setting particles.diameter = sigma gives |torque|/F_tan = {R_d:.10f}")
check("★ R = particles.diameter/2 (independent of sigma)", abs(R_d - R_RAD) < 1e-12,
      f"diameter=sigma -> R = {R_d:.10f} = sigma/2 ✓. The default diameter=1.0 is "
      f"why R came out as 0.5")
check("★★ TRAP: setting only sigma silently misaligns the lever arm",
      abs(R_MEAS / R_RAD - 1.0) > 0.1,
      f"sigma={SIGMA_CORE_STAR:.4f} but diameter defaults to 1.0 -> R is "
      f"{R_MEAS/R_RAD:.4f}x too large. No error raised")

print()
print("=" * 84)
print("T5b . C4 rolling exemption -- re-impose the no-slip condition using the "
      "MEASURED R")
print("=" * 84)
omega_roll = V / (2 * R_MEAS)
print(f"  With omega_z = V/(2R) = {omega_roll:.6f} (R={R_MEAS}, measured), the "
      f"contact point has u=0.")
print(f"  {'model':>16} | {'sliding F_tan':>14} {'rolling F_tan':>12} {'ratio':>11}")
print("  " + "-" * 60)
for name in MODELS:
    slide = probe(name, 0.95, V, 0.0)
    roll = probe(name, 0.95, V, omega_roll)
    ratio = roll["F_tan"] / slide["F_tan"] if slide["F_tan"] > 0 else float("nan")
    print(f"  {name:>16} | {slide['F_tan']:14.5g} {roll['F_tan']:12.5g} {ratio:11.3e}")
    check(f"C4 {name}: friction vanishes under no-slip rolling", ratio < 1e-9,
          f"rolling/sliding = {ratio:.3e}")

print()
print("  ⚠ sensitivity to residual slip -- Coulomb produces w(r)*kappa_f "
      "independently of u.")
print(f"  {'u/V':>8} | " + " ".join(f"{n:>16}" for n in MODELS))
print("  " + "-" * 62)
for frac in (1.0, 1e-2, 1e-4, 1e-8):
    om = (V * (1 - frac)) / (2 * R_MEAS)          # residual slip u = frac*V
    row = [probe(n, 0.95, V, om)["F_tan"] for n in MODELS]
    print(f"  {frac:8.0e} | " + " ".join(f"{x:16.6g}" for x in row))
print("  -> Linear and CoulombNewton vanish in proportion to u, but **Coulomb "
      "delivers full force even at u=1e-8**.")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T6 . C5 BD coupling -- is methods.Brownian's velocity a quantity friction "
      "may legitimately read?")
print("=" * 84)
o0 = probe("LJLinear", 0.95, 1.0, method="brownian", n_steps=0)
o1 = probe("LJLinear", 0.95, 1.0, method="brownian", n_steps=1)
o5 = probe("LJLinear", 0.95, 1.0, method="brownian", n_steps=50)
for lbl, o in (("0 steps", o0), ("1 step", o1), ("50 steps", o5)):
    print(f"  {lbl:>8}: |v| = {np.linalg.norm(o['vel'], axis=1)}   F_tan = {o['F_tan']:.5g}")
check("C5 Brownian does NOT leave velocity at 0 (friction does see something)",
      float(np.abs(o5["vel"]).max()) > 0.0,
      f"after 50 steps max|v| = {np.abs(o5['vel']).max():.4g} -- the prediction (0) "
      f"was wrong")

print()
print("  So what IS that velocity? -- compared directly against dx/dt in free BD "
      "with no forces.")
DT, MASS, KT = 1e-4, 1.0, 1.0
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=3)
f = make_frame(5.0, 0.0, 0.0)                       # far apart, so no interaction
sim.create_state_from_snapshot(f)
bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=1.0)
integ = md.Integrator(dt=DT, methods=[bd], forces=[])
integ.integrate_rotational_dof = False
sim.operations.integrator = integ
sim.run(1)
x0 = np.array(sim.state.get_snapshot().particles.position, dtype=float)
v_before = np.array(sim.state.get_snapshot().particles.velocity, dtype=float)
sim.run(1)
x1 = np.array(sim.state.get_snapshot().particles.position, dtype=float)
dx_dt = (x1 - x0) / DT
print(f"    stored v        = {v_before[0]}")
print(f"    actual dx/dt    = {dx_dt[0]}")
print(f"    ratio (v / dx/dt) = {v_before[0] / dx_dt[0]}")
check("C5' Brownian's velocity is unrelated to the actual displacement rate",
      float(np.abs(v_before[0] / dx_dt[0]).max()) < 0.5
      or float(np.abs(v_before[0] / dx_dt[0]).min()) > 2.0,
      "-> the overdamped displacement goes as sqrt(2*D*dt) ~ dt^(1/2), while v is "
      "drawn separately")

vs = []
for s in range(40):
    sm = hoomd.Simulation(device=hoomd.device.CPU(), seed=s + 1)
    sm.create_state_from_snapshot(make_frame(5.0, 0.0, 0.0))
    b = md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=1.0)
    it = md.Integrator(dt=DT, methods=[b], forces=[])
    it.integrate_rotational_dof = False
    sm.operations.integrator = it
    sm.run(2)
    vs.append(np.array(sm.state.get_snapshot().particles.velocity, dtype=float))
vs = np.concatenate(vs)
v2 = float((vs ** 2).sum(axis=1).mean())
print(f"    ⟨v²⟩ = {v2:.4f}   vs  Maxwell–Boltzmann 3kT/m = {3*KT/MASS:.4f}"
      f"   (ratio {v2/(3*KT/MASS):.3f}, N={len(vs)})")
check("C5'' velocity is redrawn from Maxwell-Boltzmann every step",
      abs(v2 / (3 * KT / MASS) - 1.0) < 0.25,
      f"<v^2>/(3kT/m) = {v2/(3*KT/MASS):.3f} -- so what the friction sees is "
      "**uncorrelated thermal noise, not real relative motion**")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
n_pass = sum(1 for _, ok, _ in results if ok)
print(f"{n_pass}/{len(results)} PASS")
for name, ok, detail in results:
    if not ok:
        print(f"  ✗ {name}   {detail}")
print("=" * 84)
sys.exit(0 if n_pass == len(results) else 1)
