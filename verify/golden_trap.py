"""Golden check for the harmonic trap -- the version with the minimum-image
convention applied.

In an earlier attempt <x^2> blew up as k got smaller (as the trap got weaker). The
cause was the periodic boundary: when a particle wraps across the box, its distance
to the fixed anchor jumps by L, and the trap applies an enormous restoring force in
the wrong direction. Applying the minimum image to the displacement resolves it.

-> This trap applies to every external.* module (added to the trap list in
   masterplan §11).
"""
import math, tempfile
from pathlib import Path
import numpy as np, gsd.hoomd, hoomd, hoomd.md as md


def lattice_frame(N, L):
    n = int(math.ceil(math.sqrt(N)))
    a = L / n
    pos = [[(i % n + .5) * a - L / 2, (i // n + .5) * a - L / 2, 0.] for i in range(N)]
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = pos
    f.particles.typeid = [0] * N
    f.particles.types = ["A"]
    f.configuration.box = [L, L, 0, 0, 0, 0]
    f.configuration.dimensions = 2
    return f


class HarmonicTrap(md.force.Custom):
    """Harmonic trap pulling each particle to its own anchor.

    The minimum image is applied to the displacement.
    """

    def __init__(self, k, anchors, box_L):
        super().__init__(aniso=False)
        self.k = float(k)
        self.anchors = np.asarray(anchors, dtype=float)
        self.L = np.array([box_L, box_L, np.inf])      # z is infinite because this is 2D (no wrapping)

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)
            pos = np.array(snap.particles.position, copy=True)
            d = pos - self.anchors[tags]
            d -= self.L * np.round(d / self.L)          # <- minimum image
            arr.force[:] = -self.k * d
            arr.potential_energy[:] = 0.5 * self.k * (d ** 2).sum(axis=1)

    def displacements(self, state):
        snap = state.get_snapshot()
        d = np.array(snap.particles.position) - self.anchors
        d -= self.L * np.round(d / self.L)
        return d[:, :2]


kT, gamma, N, L = 1.0, 1.0, 400, 60.0
print("=" * 84)
print("Harmonic trap golden check (minimum image applied):  "
      "<x^2> = kT/k,   tau_relax = gamma/k")
print(f"kT={kT}  γ={gamma}  N={N}  L={L}   [2D]")
print("=" * 84)

rows = []
for k in (2.0, 5.0, 10.0, 20.0):
    tau = gamma / k
    dt = tau / 2000
    n_eq = int(20 * tau / dt)
    n_samp, gap = 300, int(0.5 * tau / dt)

    f = lattice_frame(N, L)
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=int(k * 17) + 1)
    sim.create_state_from_snapshot(f)
    anchors = np.array(f.particles.position)
    trap = HarmonicTrap(k, anchors, L)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=kT, default_gamma=gamma)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[trap])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    sim.run(n_eq)
    samples = []
    for _ in range(n_samp):
        sim.run(gap)
        samples.append((trap.displacements(sim.state) ** 2).mean(axis=0))
    s = np.array(samples)
    mx, my = s[:, 0].mean(), s[:, 1].mean()
    mean = 0.5 * (mx + my)
    pred = kT / k
    err = 100 * (mean - pred) / pred
    sem = s.mean(axis=1).std(ddof=1) / math.sqrt(n_samp)
    ok = abs(err) < 5
    rows.append((k, tau, dt, mx, my, mean, sem, pred, err, ok))
    print(f"  k={k:5.1f}  τ={tau:.4f}  dt={dt:.2e}  steps={n_eq + n_samp * gap:>7,}")
    print(f"     <x^2>={mx:.5f}  <y^2>={my:.5f}   mean={mean:.5f} ± {sem:.5f}")
    print(f"     predicted={pred:.5f}   error {err:+6.2f}%   "
          f"{'✓ PASS' if ok else '✗ FAIL'}")

print("=" * 84)
n_ok = sum(r[-1] for r in rows)
print(f"{n_ok}/{len(rows)} PASS")

# Scaling check: <x^2>*k must equal kT independently of k
prod = np.array([r[5] * r[0] for r in rows])
print(f"\n<x^2>*k = {np.array2string(prod, precision=4)}   "
      f"(all must equal kT={kT})")
print(f"coefficient of variation {100 * prod.std() / prod.mean():.2f}%  "
      f"-> {'✓ scaling holds' if prod.std() / prod.mean() < 0.05 else '✗ scaling broken'}")
