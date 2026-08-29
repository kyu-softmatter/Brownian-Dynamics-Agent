"""조화 트랩 골든 검증 — 최소 이미지 규약 적용본.

앞선 시도에서 k가 작을수록(트랩이 약할수록) <x²>가 폭증했다. 원인은 주기경계:
입자가 박스를 넘어 wrap되면 고정 앵커까지의 거리가 L만큼 점프하고, 트랩이
거대한 잘못된 방향의 복원력을 준다. 변위에 최소 이미지를 적용하면 해소된다.

→ 이 함정은 external.* 모듈 전체에 해당한다 (마스터플랜 §11 함정 목록에 추가).
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
    """입자별 앵커로 끌어당기는 조화 트랩. 변위에 최소 이미지 적용."""

    def __init__(self, k, anchors, box_L):
        super().__init__(aniso=False)
        self.k = float(k)
        self.anchors = np.asarray(anchors, dtype=float)
        self.L = np.array([box_L, box_L, np.inf])      # z는 2D라 무한대(래핑 없음)

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)
            pos = np.array(snap.particles.position, copy=True)
            d = pos - self.anchors[tags]
            d -= self.L * np.round(d / self.L)          # ← 최소 이미지
            arr.force[:] = -self.k * d
            arr.potential_energy[:] = 0.5 * self.k * (d ** 2).sum(axis=1)

    def displacements(self, state):
        snap = state.get_snapshot()
        d = np.array(snap.particles.position) - self.anchors
        d -= self.L * np.round(d / self.L)
        return d[:, :2]


kT, gamma, N, L = 1.0, 1.0, 400, 60.0
print("=" * 84)
print("조화 트랩 골든 검증 (최소 이미지 적용):  <x²> = kT/k,   τ_relax = γ/k")
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
    print(f"     <x²>={mx:.5f}  <y²>={my:.5f}   평균={mean:.5f} ± {sem:.5f}")
    print(f"     예측={pred:.5f}   오차 {err:+6.2f}%   {'✓ PASS' if ok else '✗ FAIL'}")

print("=" * 84)
n_ok = sum(r[-1] for r in rows)
print(f"{n_ok}/{len(rows)} PASS")

# 스케일링 검증: <x²>·k 가 k에 무관하게 kT여야 함
prod = np.array([r[5] * r[0] for r in rows])
print(f"\n<x²>·k = {np.array2string(prod, precision=4)}   (전부 kT={kT}이어야 함)")
print(f"변동계수 {100 * prod.std() / prod.mean():.2f}%  "
      f"→ {'✓ 스케일링 성립' if prod.std() / prod.mean() < 0.05 else '✗ 스케일링 깨짐'}")
