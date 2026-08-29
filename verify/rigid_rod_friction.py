"""Is the translational friction of a rigid rod actually anisotropic? -- the decisive
measurement.

The claim: "bind several spheres into a rod with md.constrain.Rigid and the friction
becomes anisotropic naturally."
The test: hold the rod's orientation fixed, apply a constant force along the axis
      (parallel) and perpendicular to it in turn, and measure the terminal velocity.
      Under OverdampedViscous (no thermal noise) v = F/gamma, so gamma = F/v comes
      out directly.

      gamma_par != gamma_perp -> anisotropic (the claim holds)
      gamma_par == gamma_perp -> isotropic (the claim is wrong, need an alternative)

Theory (prolate spheroid, Perrin): gamma_perp/gamma_par ~ 1.4-2.0 for aspect ratio
p = 2-10.
Slender-body limit: gamma_perp/gamma_par -> 2 as p -> infinity.
"""
import math
import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

N_BEADS = 5          # number of spheres making up the rod
BEAD_SEP = 1.0       # sphere spacing (body frame)
GAMMA_CENTER = 1.0   # the scalar gamma assigned to the central particle
FORCE = 1.0
L_BOX = 200.0
N_STEPS = 2000
DT = 1e-3


def build(force_dir):
    """One rigid rod. Constant force along force_dir. Rotation held fixed."""
    f = gsd.hoomd.Frame()
    f.particles.N = 1
    f.particles.types = ["R", "A"]          # R = centre (the rod), A = constituent sphere
    f.particles.typeid = [0]
    f.particles.position = [[0.0, 0.0, 0.0]]
    f.particles.orientation = [[1.0, 0.0, 0.0, 0.0]]   # rod axis = body x axis = lab x axis
    f.particles.moment_inertia = [[0.0, 100.0, 100.0]]
    f.particles.body = [0]                  # the central particle takes its own tag
    f.particles.mass = [float(N_BEADS)]
    f.configuration.box = [L_BOX, L_BOX, L_BOX, 0, 0, 0]
    f.configuration.dimensions = 3

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(f)

    # Rigid-body definition: place N_BEADS spheres along the body x axis
    offs = (np.arange(N_BEADS) - (N_BEADS - 1) / 2) * BEAD_SEP
    rigid = md.constrain.Rigid()
    rigid.body["R"] = {
        "constituent_types": ["A"] * N_BEADS,
        "positions": [[float(o), 0.0, 0.0] for o in offs],
        "orientations": [[1.0, 0.0, 0.0, 0.0]] * N_BEADS,
    }
    rigid.create_bodies(sim.state)

    const = md.force.Constant(filter=hoomd.filter.Rigid(("center",)))
    const.constant_force["R"] = tuple(float(x) * FORCE for x in force_dir)
    const.constant_torque["R"] = (0.0, 0.0, 0.0)
    const.constant_force["A"] = (0.0, 0.0, 0.0)
    const.constant_torque["A"] = (0.0, 0.0, 0.0)

    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Rigid(("center",)),
                                      default_gamma=GAMMA_CENTER)
    integ = md.Integrator(dt=DT, methods=[ov], forces=[const],
                          rigid=rigid, integrate_rotational_dof=False)
    sim.operations.integrator = integ
    return sim


def terminal_gamma(force_dir):
    sim = build(force_dir)
    p0 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
    sim.run(N_STEPS)
    p1 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
    disp = p1 - p0
    t = N_STEPS * DT
    v = float(np.dot(disp, np.array(force_dir, dtype=float)) / t)
    return FORCE / v if abs(v) > 1e-14 else float("inf"), disp, v


print("=" * 84)
print("rigid-rod translational friction anisotropy check")
print(f"rod: {N_BEADS} spheres, spacing {BEAD_SEP}, along the body x axis  |  "
      f"central gamma={GAMMA_CENTER} (scalar)")
print("=" * 84)

g_par, d_par, v_par = terminal_gamma((1, 0, 0))     # along the axis (parallel)
g_perp, d_perp, v_perp = terminal_gamma((0, 1, 0))  # perpendicular

print(f"  parallel   displacement={np.array2string(d_par, precision=5)}  "
      f"v={v_par:.6f}  gamma_par ={g_par:.6f}")
print(f"  perpendicular displacement={np.array2string(d_perp, precision=5)}  "
      f"v={v_perp:.6f}  gamma_perp ={g_perp:.6f}")
print()
ratio = g_perp / g_par
print(f"  γ⊥/γ∥ = {ratio:.6f}")
print()
if abs(ratio - 1.0) < 0.02:
    print("  ✗ isotropic. Even bound as a rigid body, the translational friction is "
          "set by the single scalar gamma of the central particle.")
    print("    -> the claim that 'binding it rigid gives anisotropy naturally' does "
          "NOT hold in HOOMD.")
else:
    print(f"  ✓ anisotropic (gamma_perp/gamma_par = {ratio:.3f}). Theory expects "
          f"1.4-2.0 (slender body -> 2).")
print("=" * 84)

# Aside: does giving the constituent particles a gamma change anything?
print()
print("secondary check: does a separate gamma on constituent type 'A' affect the "
      "translation at all?")
sim = build((1, 0, 0))
ov = sim.operations.integrator.methods[0]
ov.gamma["A"] = 10.0          # large drag on the constituent particles
p0 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
sim.run(N_STEPS)
p1 = np.array(sim.state.get_snapshot().particles.position[0], dtype=float)
v2 = float((p1 - p0)[0] / (N_STEPS * DT))
print(f"  after setting gamma['A']=10  v_par={v2:.6f}  "
      f"(previously v_par={v_par:.6f})")
print("  -> if the value is unchanged, the constituent particles' gamma plays no "
      "part in translation at all.")
