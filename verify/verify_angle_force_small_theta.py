"""★★ Trap reproduction -- `md.angle.Harmonic(t0=pi)` gets the ENERGY right and the
**FORCE wrong** for a **nearly straight** chain.

How it was found: `chain-bend-2d-oscill`'s K*(omega) disagreed with the
linear-response prediction by 28% at high frequency. Narrowing it down a rung at a
time (`scratch/diagnose_chain_bend_28pct.py`) --
  . dt discretisation: 0.002% against the exact z-domain solution -> not it
  . nonlinearity, x degrees of freedom, bond stretching: exact minimization with
    scipy agrees with the linear model to 0.32% -> not it
  . energy: agrees with HOOMD to **exactly 0.000%** (in both configurations)
  . force: ✗ **completely different in the nearly straight configuration** <- cause

Why it is silently wrong -- an angle force uses `1/sin theta` to carry the torque
into coordinates. With t0=pi the equilibrium itself is the singular point
sin theta = 0, and a stiff chain (kappa_theta = 1.4e6 kT) is **always** nearly
straight. So this trap bites **hardest** precisely for a rigid filament. The energy
is computed directly from arccos and is correct, so **verifying with energy passes
and leaves the force wrong** -- which is exactly why
`scratch/verify_angle_matrix.py`, measuring via 2U/delta^2, passed at 0.55%.

Lesson: **verifying a potential by its energy alone is not verifying its force.**

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_angle_force_small_theta.py
"""
from __future__ import annotations

import math

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

K_ANGLE = 1.3912297767209476e6      # chain-bend's kappa_theta (in kT) -- tested with the real case value
N = 5                                # the smallest chain that can be checked by hand
TOL = 1e-3                           # a force ratio deviating from 1 by more than this counts
                           # as wrong


def energy_exact(r, k=K_ANGLE, t0=math.pi):
    v1 = r[:-2] - r[1:-1]
    v2 = r[2:] - r[1:-1]
    c = np.sum(v1 * v2, axis=1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
    return 0.5 * k * np.sum((np.arccos(np.clip(c, -1, 1)) - t0) ** 2)


def force_exact(r, k=K_ANGLE, t0=math.pi):
    """Analytic gradient of the energy (Richardson-extrapolated finite difference).

    Force = -grad U.
    """
    G = np.zeros_like(r)
    h0 = 1e-6
    for i in range(r.shape[0]):
        for d in range(2):
            e = []
            for h in (h0, h0 / 2):
                a = r.copy(); a[i, d] += h
                b = r.copy(); b[i, d] -= h
                e.append((energy_exact(a, k, t0) - energy_exact(b, k, t0)) / (2 * h))
            G[i, d] = (4 * e[1] - e[0]) / 3
    return -G


def hoomd_angle(r, kind="harmonic", k=K_ANGLE, t0=math.pi):
    n = r.shape[0]
    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.position = np.column_stack([r, np.zeros(n)])
    f.particles.typeid = [0] * n
    f.particles.types = ["A"]
    f.configuration.box = [400.0, 400.0, 0, 0, 0, 0]     # Lz=0 -> 2D (trap 9)
    f.configuration.dimensions = 2
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)

    if kind == "harmonic":
        ang = md.angle.Harmonic()
        ang.params["bend"] = dict(k=k, t0=t0)
    else:                                                 # tabulated angle -- the same U(theta)
                                                 # supplied directly
        nb = 1000
        th = np.linspace(0, math.pi, nb, endpoint=False)
        ang = md.angle.Table(width=nb)
        ang.params["bend"] = dict(U=0.5 * k * (th - t0) ** 2, tau=-k * (th - t0))
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0)
    integ = md.Integrator(dt=1e-14, methods=[ov], forces=[ang])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(0)
    F = np.array(ang.forces)[:, :2].copy()
    with sim.state.cpu_local_snapshot as s:
        tags = np.array(s.particles.tag, copy=True)
    out = np.zeros_like(F)
    out[tags] = F                                         # local order -> tag order
    return out, float(ang.energy)


def bent(scale, n=N, ell=1.0):
    """A smooth bending mode of amplitude `scale`. As scale -> 0, sin theta -> 0."""
    y = scale * np.sin(np.pi * np.arange(n) / (n - 1))
    return np.column_stack([np.arange(n) * ell, y])


def theta_max(r):
    v1 = r[:-2] - r[1:-1]
    v2 = r[2:] - r[1:-1]
    c = np.sum(v1 * v2, axis=1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
    return float(np.max(np.abs(np.arccos(np.clip(c, -1, 1)) - math.pi)))


def main() -> int:
    print("=" * 96)
    print("md.angle.Harmonic(t0=pi) -- energy right, force wrong for a nearly "
          "straight chain")
    print("=" * 96)
    print(f"kappa_theta = {K_ANGLE:.4e} (chain-bend's real value) . n = {N}")
    print(f"{'max|th-pi|':>11}{'sin th':>11}{'E err%':>10}{'F ratio(max)':>15}"
          f"{'mid bead |F|':>14}{'mid exact':>12}{'verdict':>7}")
    print("-" * 96)
    first_bad = None
    for scale in (3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1.5e-3, 1e-3, 7e-4, 5e-4, 3e-4, 1e-4):
        r = bent(scale)
        tm = theta_max(r)
        Fh, Eh = hoomd_angle(r)
        Fa = Ea = None
        Fa = force_exact(r)
        Ea = energy_exact(r)
        j = int(np.argmax(np.abs(Fa[:, 1])))
        ratio = Fh[j, 1] / Fa[j, 1]
        mid = N // 2
        ok = abs(ratio - 1.0) < TOL
        if not ok and first_bad is None:
            first_bad = (tm, math.sin(math.pi - tm))
        print(f"{tm:>11.3e}{math.sin(math.pi-tm):>11.3e}{100*(Eh/Ea-1):>10.4f}"
              f"{ratio:>15.6f}{abs(Fh[mid,1]):>14.4e}{abs(Fa[mid,1]):>12.4e}"
              f"{'✓' if ok else '✗':>7}")
    print("-" * 96)
    if first_bad:
        print(f"★ where the force starts to break: max|theta-pi| ~ {first_bad[0]:.2e}"
              f"  (sin theta ~ {first_bad[1]:.2e})")
        print("  -> safe rule: the force is only right if "
              "**max|theta-pi| >~ 3e-3 (sin theta >~ 3e-3) is maintained**.")
    else:
        print("the force is right throughout -- the trap does not reproduce on this "
              "install")

    # Alternative (1) angle.Table -- does supplying the same U(theta) directly take
    # the same code path?
    print("\n" + "-" * 96)
    print("alternatives considered")
    print("-" * 96)
    r = bent(3e-4)
    Fa = force_exact(r)
    j = int(np.argmax(np.abs(Fa[:, 1])))
    try:
        Ft, Et = hoomd_angle(r, kind="table")
        print(f"  angle.Table  (same U and torque, tabulated)  "
              f"F ratio = {Ft[j,1]/Fa[j,1]:.6f}"
              f"   {'✓ a viable workaround' if abs(Ft[j,1]/Fa[j,1]-1)<1e-2 else '✗ same problem'}")
    except Exception as e:
        print(f"  angle.Table  → {type(e).__name__}: {str(e)[:70]}")
    # Alternative (2) CosineSquared is **not** harmonic at t0=pi (ruled out
    # analytically)
    print("  angle.CosineSquared: U = 0.5*k*(cos th - cos th0)^2. At th0=pi, "
          "cos th + 1 ~ (pi - th)^2/2, so")
    print("    U ~ k(pi - th)^4/8 -- **quartic**. It cannot produce a harmonic "
          "bending stiffness -> not a substitute (ruled out analytically)")
    print("=" * 96)
    return 0 if first_bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
