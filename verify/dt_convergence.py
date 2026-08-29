"""dt 수렴 연구 — "dt/τ를 얼마나 작게?"에 대한 정량적 답.

조화 트랩은 선형계라 Euler-Maruyama의 이산 정상분산을 해석적으로 구할 수 있다:

    x_{n+1} = x_n (1 - h) + sqrt(2 D dt) ξ,      h ≡ dt/τ_k
    정상상태:  <x²> = (kT/k) / (1 - h/2)
    ⟹ 상대 편향 = 1/(1 - h/2) - 1 ≈ h/2

즉 **dt/τ 의 절반이 곧 계통 편향**이다. 이게 "dt/τ ≤ 1e-2" 규칙의 정체.
HOOMD의 Brownian이 실제로 이 법칙을 따르는지 확인한다.
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
    """h = dt/tau_k 에서 <x²> 측정."""
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

    gap = max(1, int(round(2 * TAU_K / dt)))       # 표본 간격 2 τ_k
    sim.run(int(20 * TAU_K / dt))                  # 평형화 20 τ_k
    s = []
    for _ in range(N_SAMPLES):
        sim.run(gap)
        s.append((trap.disp(sim.state) ** 2).mean())
    s = np.array(s)
    return s.mean(), s.std(ddof=1) / math.sqrt(len(s)), sim.timestep


print("=" * 88)
print(f"dt 수렴 연구 — 조화 트랩,  kT={KT} γ={GAMMA} k={K}  →  τ_k={TAU_K}  <x²>_exact={PRED}")
print(f"N={N} 입자, 표본 {N_SAMPLES}개 (간격 2τ_k)")
print("=" * 88)
print(f"{'dt/τ_k':>8} {'dt':>10} {'<x²> 측정':>14} {'±SEM':>9} "
      f"{'편향 측정':>10} {'편향 이론':>10} {'차이':>9}  {'스텝':>9}")
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
print(f"\n측정/이론 비: {np.array2string(meas/theo, precision=3)}")
ok = np.allclose(meas, theo, rtol=0.25)
print(f"{'✓ Euler-Maruyama 편향 법칙 확인 — 편향 ≈ (dt/τ)/2' if ok else '✗ 법칙과 불일치'}")
print()
print("실용 규칙 (이 표에서 직접 읽음):")
for h, want in [(0.1, None), (0.02, None), (0.01, None), (0.002, None), (0.0005, None)]:
    print(f"   dt/τ = {h:<7.4g} → 계통 편향 ≈ {100*(1/(1-h/2)-1):.3f}%")
