"""Compare the bending stiffness matrix against **HOOMD's `angle.Harmonic` itself**.

★ Why this is needed -- all that had been verified so far was
`discrete matrix vs beam formula 48EI/L^3` (-0.35%), and that is **analysis against
analysis**. Whether `cases/chain_bend_2d.bending_matrix` produces the same physics as
HOOMD's `md.angle.Harmonic(k=kappa_theta, t0=pi)` had **never once been checked.**
Rules 4 and 6 point at exactly this situation (do not write it from intuition;
settle it by execution).

Setup: the static three-point bending problem. Beads 0, 12 and 24 are pinned at
y = 0, delta, 0 (excluded from integration), the free beads are relaxed, and the
bending + bond force on bead 12 is read to give k_center = -F_y/delta.
No traps are used, so no trap stiffness is mixed in -- this sees the chain's
stiffness alone.

Three things compared:
  (1) the discrete matrix (the static solution from bending_matrix)
  (2) the continuum beam 48EI/L^3
  (3) HOOMD, measured                            <- this is the new measurement

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_angle_matrix.py
"""
from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from cases.chain_bend_2d import bending_matrix, sweep_specs, trapped_indices  # noqa: E402

TOL = 1.0          # tolerance [%] -- the discrete matrix and HOOMD must be the same physics


def matrix_static(n, kth, ell, delta, fixed):
    """Static three-point bending solved with the discrete matrix.

    k_center = 2U/delta^2 (the same definition as three_point_bending).
    """
    A = bending_matrix(n, kth, ell)
    free = [i for i in range(n) if i not in fixed]
    y = np.zeros(n)
    y[fixed] = [0.0, delta, 0.0]
    y[free] = np.linalg.solve(A[np.ix_(free, free)], -A[np.ix_(free, fixed)] @ y[fixed])
    U = 0.5 * y @ A @ y
    return 2.0 * U / delta ** 2, y


def hoomd_static(n, kth, k_bond, ell, delta, fixed, dt, n_steps):
    """Solve the same problem in HOOMD.

    The pinned beads are given type 'F' and excluded from integration.
    """
    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    pos[fixed[1]][1] = delta                      # push the centre by delta
    typeid = [1 if i in fixed else 0 for i in range(n)]

    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "F"]
    f.configuration.box = [40.0 * ell, 40.0 * ell, 0, 0, 0, 0]   # Lz=0 -> 2D (trap 9)
    f.configuration.dimensions = 2
    f.bonds.N = n - 1
    f.bonds.types = ["backbone"]
    f.bonds.typeid = [0] * (n - 1)
    f.bonds.group = np.array([[i, i + 1] for i in range(n - 1)])
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=k_bond, r0=ell)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=kth, t0=math.pi)
    # Overdamped with no thermal noise -- relax to the static solution
    # (skill bd-hoomd)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Type(["A"]), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[ov], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(n_steps)

    # ★ The energy is a scalar so it has no particle-ordering problem (the
    #   per-particle force array may not be in tag order because of ParticleSorter
    #   -- skill bd-hoomd). So it is measured as k_center = 2U/delta^2.
    U = float(angle.energy) + float(bond.energy)
    snap = sim.state.get_snapshot()               # a global snapshot IS in tag order
    y = np.array(snap.particles.position)[:, 1]
    return 2.0 * U / delta ** 2, y, (float(angle.energy), float(bond.energy))


def main() -> int:
    # ★ Deliberately does not depend on a spec file -- when a hard check fails, the
    #   absence of a spec is the correct behaviour
    p = sweep_specs()[0]["params"]
    n = int(p["n_beads"])
    kth = float(p["kappa_theta_star"])
    k_bond = float(p["k_bond_star"])
    ell = float(p["L_chain_star"]) / (n - 1)
    L = float(p["L_chain_star"])
    delta = float(p["amp_star"])
    fixed = trapped_indices(n)
    EI = kth * ell
    beam = 48 * EI / L ** 3

    lam_max = float(np.linalg.eigvalsh(bending_matrix(n, kth, ell))[-1])
    dt = 0.2 / lam_max                    # 10% of the explicit-Euler stability limit (2/lambda)
    n_steps = int(round(20.0 / (0.2 * 526.4 / lam_max)))   # ≈ 20 τ_max

    print("=" * 88)
    print("bending stiffness: discrete matrix vs continuum beam vs "
          "**HOOMD angle.Harmonic, measured**")
    print("=" * 88)
    print(f"n={n}  kappa_theta*={kth:.4e}  l={ell:.3f}  L={L:.1f}  "
          f"delta={delta:.5f}  pinned={fixed}")
    print(f"lambda_max={lam_max:.4e}  dt={dt:.3e}  steps={n_steps:,}\n")

    k_mat, y_mat = matrix_static(n, kth, ell, delta, fixed)
    k_hd, y_hd, (U_ang, U_bond) = hoomd_static(n, kth, k_bond, ell, delta, fixed,
                                               dt, n_steps)

    print(f"(1) discrete matrix      k_center = {k_mat:10.2f}")
    print(f"(2) continuum beam 48EI/L^3  k_center = {beam:10.2f}   "
          f"(vs (1): {100*(k_mat/beam-1):+.2f}%)")
    print(f"(3) HOOMD, measured      k_center = {k_hd:10.2f}   (vs (1): "
          f"{100*(k_hd/k_mat-1):+.2f}%   vs (2): {100*(k_hd/beam-1):+.2f}%)")
    print(f"\nHOOMD energy breakdown: bending {U_ang:.4f} + bond {U_bond:.4e} kT"
          f"  (bond share {100*U_bond/(U_ang+U_bond):.3f}%)")
    print(f"largest difference in the displacement profile (matrix vs HOOMD) = "
          f"{100*np.max(np.abs(y_hd-y_mat))/delta:.3f}% of δ")

    err = 100 * abs(k_hd / k_mat - 1)
    ok = err < TOL
    print("-" * 88)
    if ok:
        print(f"✓ PASS -- matrix and HOOMD agree to within {err:.2f}%. "
              f"bending_matrix represents angle.Harmonic correctly")
    else:
        print(f"✗ FAIL -- matrix and HOOMD differ by {err:.2f}%. They are **not the "
              f"same physics** -- every quantity derived from the matrix "
              f"(lambda_max, tau_fast, dt, the K*(omega) prediction) is affected")
    print("=" * 88)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
