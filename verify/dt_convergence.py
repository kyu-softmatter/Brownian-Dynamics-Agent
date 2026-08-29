"""dt convergence study -- a quantitative answer to "how small must dt/tau be?".

The harmonic trap is a linear system, so the discrete stationary variance of
Euler-Maruyama can be obtained analytically:

    x_{n+1} = x_n (1 - h) + sqrt(2 D dt) ξ,      h ≡ dt/τ_k
    stationary:  <x^2> = (kT/k) / (1 - h/2)
    => relative bias = 1/(1 - h/2) - 1 ~ h/2

In other words **half of dt/tau IS the systematic bias**. That is what the
"dt/tau <= 1e-2" rule actually is. This checks whether HOOMD's Brownian really
follows that law.
"""
import math
import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

KT, GAMMA, K = 1.0, 1.0, 5.0
TAU_K = GAMMA / K                      # 0.2
PRED = KT / K                          # 0.2
N = 2000
N_SAMPLES = 1000


class Trap(md.force.Custom):
    def __init__(self, k, anchors, box_L, dim=2):
        super().__init__(aniso=False)
        self.k = float(k)
        self.anchors = np.asarray(anchors, float)
        self.period = np.array([box_L, box_L, box_L if dim == 3 else 0.0])

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)
            d = np.array(snap.particles.position, copy=True) - self.anchors[tags]
            m = self.period > 0
            d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
            arr.force[:] = -self.k * d
            arr.potential_energy[:] = 0.5 * self.k * (d ** 2).sum(axis=1)

    def disp(self, state):
        d = np.array(state.get_snapshot().particles.position) - self.anchors
        m = self.period > 0
        d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
        return d[:, :2]


def run(h):
    """Measure <x^2> at h = dt/tau_k."""
    dt = h * TAU_K
    L = 200.0
    n = int(math.ceil(math.sqrt(N)))
    a = L / n
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = [[(i % n + .5) * a - L/2, (i // n + .5) * a - L/2, 0.]
                            for i in range(N)]
    f.particles.typeid = [0]*N
    f.particles.types = ["A"]
    f.configuration.box = [L, L, 0, 0, 0, 0]
    f.configuration.dimensions = 2

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1234)
    sim.create_state_from_snapshot(f)
    anchors = np.array(f.particles.position)
    trap = Trap(K, anchors, L)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=GAMMA)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[trap])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    gap = max(1, int(round(2 * TAU_K / dt)))       # sample interval 2 tau_k
    sim.run(int(20 * TAU_K / dt))                  # equilibrate for 20 tau_k
    s = []
    for _ in range(N_SAMPLES):
        sim.run(gap)
        s.append((trap.disp(sim.state) ** 2).mean())
    s = np.array(s)
    return s.mean(), s.std(ddof=1) / math.sqrt(len(s)), sim.timestep


print("=" * 88)
print(f"dt convergence study -- harmonic trap,  kT={KT} gamma={GAMMA} k={K}  ->  "
      f"tau_k={TAU_K}  <x^2>_exact={PRED}")
print(f"N={N} particles, {N_SAMPLES} samples (interval 2 tau_k)")
print("=" * 88)
print(f"{'dt/tau_k':>8} {'dt':>10} {'<x^2> meas':>14} {'±SEM':>9} "
      f"{'bias meas':>10} {'bias theory':>10} {'diff':>9}  {'steps':>9}")
print("-" * 88)

rows = []
for h in (0.1, 0.05, 0.02, 0.01):
    m, sem, steps = run(h)
    bias = 100 * (m - PRED) / PRED
    bias_th = 100 * (1 / (1 - h / 2) - 1)
    rows.append((h, bias, bias_th))
    print(f"{h:>8.3f} {h*TAU_K:>10.2e} {m:>14.6f} {sem:>9.6f} "
          f"{bias:>9.3f}% {bias_th:>9.3f}% {bias-bias_th:>+8.3f}%  {steps:>9,}")

print("-" * 88)
meas = np.array([r[1] for r in rows])
theo = np.array([r[2] for r in rows])
print(f"\nmeasured/theory ratio: {np.array2string(meas/theo, precision=3)}")
ok = np.allclose(meas, theo, rtol=0.25)
print(f"{'✓ Euler-Maruyama bias law confirmed -- bias ~ (dt/tau)/2' if ok else '✗ does not match the law'}")
print()
print("Practical rule (read straight off this table):")
for h, want in [(0.1, None), (0.02, None), (0.01, None), (0.002, None), (0.0005, None)]:
    print(f"   dt/tau = {h:<7.4g} -> systematic bias ~ {100*(1/(1-h/2)-1):.3f}%")
