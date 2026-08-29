"""(1) Does per-axis anisotropy in gamma_r actually work?
   (2) Does a free-draining bead chain give anisotropic translational friction?
"""
import numpy as np, gsd.hoomd, hoomd, hoomd.md as md

# ── (1) check that per-axis rotational drag gamma_r works ───────────────
print("=" * 84)
print("(1) per-axis anisotropy of gamma_r -- angular velocity under constant torque "
      "(Langevin, kT=0 -> deterministic)")
print("=" * 84)

def spin(gamma_r, torque_axis, steps=2000, dt=1e-3):
    f = gsd.hoomd.Frame()
    f.particles.N = 1
    f.particles.types = ["R"]
    f.particles.typeid = [0]
    f.particles.position = [[0, 0, 0]]
    f.particles.orientation = [[1, 0, 0, 0]]
    f.particles.moment_inertia = [[1.0, 1.0, 1.0]]
    f.particles.mass = [1.0]
    f.configuration.box = [50, 50, 50, 0, 0, 0]
    f.configuration.dimensions = 3
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=3)
    sim.create_state_from_snapshot(f)

    const = md.force.Constant(filter=hoomd.filter.All())
    const.constant_force["R"] = (0, 0, 0)
    const.constant_torque["R"] = tuple(float(x) for x in torque_axis)

    # Brownian with kT=0 -> no noise, dq/dt = torque/gamma_r
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0,
                             default_gamma=1.0, default_gamma_r=tuple(gamma_r))
    integ = md.Integrator(dt=dt, methods=[bd], forces=[const],
                          integrate_rotational_dof=True)
    sim.operations.integrator = integ
    sim.run(steps)
    q = np.array(sim.state.get_snapshot().particles.orientation[0], dtype=float)
    # rotation angle = 2*arccos(w)
    return 2 * np.arccos(np.clip(abs(q[0]), -1, 1))

for gr, label in [((1.0, 1.0, 1.0), "isotropic (1,1,1)"),
                  ((1.0, 1.0, 5.0), "z axis 5x (1,1,5)")]:
    ax = spin(gr, (1.0, 0, 0))
    az = spin(gr, (0, 0, 1.0))
    print(f"  gamma_r={gr:} [{label}]")
    print(f"     torque||x -> angle {ax:.5f} rad     "
          f"torque||z -> angle {az:.5f} rad     "
          f"ratio z/x = {az/ax if ax else 0:.4f}")
print("  -> making the z component 5x should slow z-axis rotation to 1/5 "
      "(ratio ~ 0.2)")

# ── (2) free-draining bead chain: does anisotropy appear? ───────────────
print()
print("=" * 84)
print("(2) free-draining bead rod (not rigid; each bead has independent drag) -- "
      "does anisotropy appear?")
print("=" * 84)

def bead_rod_gamma(force_dir, n_beads=5, k_bond=2000.0, steps=4000, dt=1e-4):
    f = gsd.hoomd.Frame()
    f.particles.N = n_beads
    f.particles.types = ["A"]
    f.particles.typeid = [0] * n_beads
    offs = (np.arange(n_beads) - (n_beads - 1) / 2) * 1.0
    f.particles.position = [[float(o), 0.0, 0.0] for o in offs]
    f.particles.mass = [1.0] * n_beads
    f.configuration.box = [200, 200, 200, 0, 0, 0]
    f.configuration.dimensions = 3
    f.bonds.N = n_beads - 1
    f.bonds.types = ["b"]
    f.bonds.typeid = [0] * (n_beads - 1)
    f.bonds.group = [[i, i + 1] for i in range(n_beads - 1)]
    f.angles.N = n_beads - 2
    f.angles.types = ["a"]
    f.angles.typeid = [0] * (n_beads - 2)
    f.angles.group = [[i, i + 1, i + 2] for i in range(n_beads - 2)]

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=5)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic(); bond.params["b"] = dict(k=k_bond, r0=1.0)
    ang = md.angle.Harmonic(); ang.params["a"] = dict(k=k_bond, t0=np.pi)
    const = md.force.Constant(filter=hoomd.filter.All())
    # F/n per bead, so the total force on the whole rod is F
    const.constant_force["A"] = tuple(float(x) / n_beads for x in force_dir)
    const.constant_torque["A"] = (0, 0, 0)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=dt, methods=[ov],
                                              forces=[bond, ang, const])
    com0 = np.array(sim.state.get_snapshot().particles.position).mean(axis=0)
    sim.run(steps)
    com1 = np.array(sim.state.get_snapshot().particles.position).mean(axis=0)
    v = float(np.dot(com1 - com0, np.array(force_dir, float)) / (steps * dt))
    return 1.0 / v if abs(v) > 1e-14 else float("inf")   # F_total=1, so gamma=1/v

gp = bead_rod_gamma((1.0, 0, 0))
gt = bead_rod_gamma((0, 1.0, 0))
print(f"  γ∥ = {gp:.5f}   γ⊥ = {gt:.5f}   γ⊥/γ∥ = {gt/gp:.6f}")
print(f"  (5 beads x gamma_bead=1 -> total gamma = 5 expected)")
print()
print("  -> in free draining each bead feels Stokes drag independently, so")
print("    the total drag = N*gamma_bead, independent of direction. Anisotropy is an "
      "effect of hydrodynamic interaction (HI),")
print("    and in BD without HI, geometry alone does not produce anisotropy.")
