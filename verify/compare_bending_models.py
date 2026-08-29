"""조화 굽힘(JKR 선형화) vs 구름 저항 — 같은 기하·같은 구동에서 동적 응답 비교.

정적으로는 이미 결론이 났다 (`scratch/verify_rolling_contact.py`, 24/24):
  · 접선 스프링만  → 굽힘강성 **정확히 0** (δ×100 에서 잔차가 1/10⁴ 로 → 반올림)
  · 구름 저항만    → κ_θ,eff = ½k_rR² 로 조화굽힘과 **1e−5 이내 일치** (방향 이완 극한)
  · 방향을 **얼리면** 6·22·66배 뻣뻣해진다 (곡률이 아니라 절대 회전을 벌함)

⟹ 남은 질문은 하나다: **이 계의 구동 주파수에서 방향 자유도가 이완할 시간이 있는가.**
   있으면 두 모델은 같은 물리이고, 없으면 구름 모델이 훨씬 뻣뻣해진다.
   교차 주파수 예측: `ω_c = 1/τ_rot`,  `τ_rot = γ_r/(k_r R²)`.

★ kT=0 결정론적으로 잰다 — 이 프로젝트의 교훈("적분기 가정을 시험할 때는 kT=0").
   잡음이 0 이라 4주기면 락인이 깨끗하고, 두 모델의 과도가 공통모드로 상쇄된다.
★ 규칙 7(격리): DLVO 표 대신 **우물 곡률과 같은 강성의 방사 본드**를 쓴다. 굽힘 모델
   두 개의 차이만 남기기 위해서다 (DLVO 의 비선형은 이 질문과 무관하다).

    $PY scratch/compare_bending_models.py            # 스윕 + 그래프
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scratch"))

from bdbot import lockin as LI                                          # noqa: E402
from rolling_contact import (k_roll_from_kappa_theta,                   # noqa: E402
                             make_rolling_force)

# ── chain-bend-2d-dlvo (n=9, ω=3000, a=632nm, --jkr) 스펙의 실제 값 ──────────
H_MIN_STAR = 0.00759259035993831
ELL = 1.0 + H_MIN_STAR
K_BOND = 1042362.8817700658          # DLVO 2차극소 곡률 [kT/d²] — 방사 강성
KAPPA_THETA = 1391229.7767209478     # [kT] — κ₀=64 mN/m → EI/ℓ
# ★ 트랩을 방사 본드만큼 뻣뻣하게 쓴다 (생산값 5217 이 아니다).
#   왜: 생산 트랩의 완화 모드가 τ = γ/k_t = 1.9e−4 로 **최고 ω 주기의 3만 배**라,
#   주기 단위로 정착시키면 과도가 안 빠져 락인이 오염된다 (첫 스윕에서 K″<0, 비가
#   0.44↔1.39 로 요동쳤다 — 물리가 아니라 미정착이었다). 여기서 묻는 것은 굽힘 모델
#   두 개의 차이뿐이므로 트랩은 빠른 경계조건이면 된다 (규칙 7 격리).
K_T = K_BOND                         # 트랩 강성 [kT/d²]
AMP = 0.05                           # 구동 진폭 [d] — kT=0 이라 선형응답이면 값 무관
R_C = 0.5                            # 접촉 반경 = d/2
GAMMA_R = 4 * R_C ** 2 / 3           # γ_r/γ_t = 8πηa³/6πηa = 4a²/3  (a=R_C, γ_t=1)
K_ROLL = k_roll_from_kappa_theta(KAPPA_THETA, R_C)
TAU_ROT = GAMMA_R / (K_ROLL * R_C ** 2)
N_BEADS = 5
UPDATE_EVERY = 1                     # 유령을 매 스텝 옮긴다 (ZOH 감쇠 최소화 — 함정 17)
OUT = ROOT / "runs" / "_bending_model_compare"


def bending_matrix(n, kappa_theta, ell):
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def make_harmonic_bending(A, n_real):
    """조화 굽힘 F_y = −A y (chain-bend-2d-dlvo 의 --jkr 구현과 동일)."""
    class Bending(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=False)
            self.A = np.ascontiguousarray(A, dtype=float)
            self.n = int(n_real)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)
                pos = np.array(snap.particles.position, copy=True)
                m = tags < self.n
                y = np.zeros(self.n)
                y[tags[m]] = pos[m, 1]
                fy = -(self.A @ y)
                arr.force[:] = 0.0
                arr.potential_energy[:] = 0.0
                arr.force[m, 1] = fy[tags[m]]
                arr.potential_energy[m] = -0.5 * y[tags[m]] * fy[tags[m]]
    return Bending()


class ClampAndDrive(hoomd.custom.Action):
    """★ 변형률 제어 — 양끝 비드의 (x,y) 를 고정하고 중앙 비드의 y 를 직접 강제한다.

    왜 트랩을 안 쓰는가: 트랩 구동은 컴플라이언스 때문에 ωτ_trap≫1 에서 추종률이
    1% 아래로 떨어지고, `K* = k_t·ŷ_c/ŷ − k_t` 가 큰 수의 비가 되어 조건수가 무너진다
    (첫 스윕 실측 — 추종 0.01~0.00001, K′ 이 음수로 뒤집혔다). 위치 강제는 추종이
    정의상 100% 라 **모든 ω 에서 조건이 같다**.
    ★ 비드는 적분기 filter 에 **남겨 둔다** — 그래야 회전 자유도가 계속 적분된다
    (G4 에서 끝 입자의 방향 이완이 강성의 일부였다). 위치만 매 스텝 덮어쓴다.
    """

    def __init__(self, clamp_tags, clamp_xy, drive_tag, amp, omega, dt):
        self.tags = np.asarray(clamp_tags, int)
        self.xy = np.asarray(clamp_xy, float)
        self.drive_tag, self.amp, self.omega, self.dt = int(drive_tag), amp, omega, dt

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            for k, tg in enumerate(self.tags):
                idx = snap.particles.rtag[tg]
                snap.particles.position[idx, 0] = self.xy[k, 0]
                snap.particles.position[idx, 1] = (
                    self.amp * math.sin(self.omega * timestep * self.dt)
                    if tg == self.drive_tag else self.xy[k, 1])


def build(model: str, omega: float, dt: float, n=N_BEADS):
    mid = n // 2
    clamped = [0, mid, n - 1]
    pos0 = np.zeros((n, 3))
    pos0[:, 0] = (np.arange(n) - (n - 1) / 2) * ELL
    quat0 = np.tile(np.array([1.0, 0, 0, 0]), (n, 1))

    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.types = ["A"]
    f.particles.typeid = [0] * n
    f.particles.position = pos0
    f.particles.orientation = [(1, 0, 0, 0)] * n
    f.particles.moment_inertia = [(1.0, 1.0, 1.0)] * n
    f.particles.diameter = [1.0] * n            # bd-hoomd 함정 19 (여기선 안 쓰지만 습관)
    f.bonds.N = n - 1
    f.bonds.types = ["radial"]
    f.bonds.typeid = [0] * (n - 1)
    f.bonds.group = [[i, i + 1] for i in range(n - 1)]
    L = 8.0 * n * ELL
    f.configuration.box = [L, L, 0, 0, 0, 0]
    f.configuration.dimensions = 2

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["radial"] = dict(k=K_BOND, r0=ELL)
    forces = [bond]
    if model == "harmonic":
        forces.append(make_harmonic_bending(bending_matrix(n, KAPPA_THETA, ELL), n))
    elif model == "rolling":
        forces.append(make_rolling_force([[i, i + 1] for i in range(n - 1)], pos0, quat0,
                                         R_C, K_ROLL, 0.0, n))
    elif model != "none":
        raise ValueError(model)

    # kT=0 결정론적 과감쇠 (OverdampedViscous — 잡음 없음)
    meth = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0,
                                        default_gamma_r=(GAMMA_R,) * 3)
    integ = md.Integrator(dt=dt, methods=[meth], forces=forces)
    integ.integrate_rotational_dof = (model == "rolling")
    sim.operations.integrator = integ
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=ClampAndDrive(clamped, pos0[clamped, :2], mid, AMP, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    return sim, mid, forces


def timescales(model: str, n=N_BEADS):
    """이 계의 시간척도. 가장 빠른 것이 dt 를, 가장 느린 것이 **정착 시간**을 정한다."""
    A = bending_matrix(n, KAPPA_THETA, ELL)
    mid = n // 2
    free = [i for i in range(n) if i not in (0, mid, n - 1)]
    lam = np.linalg.eigvalsh(A)
    lam_ff = np.linalg.eigvalsh(A[np.ix_(free, free)])       # 자유 비드의 완화
    fast = [1.0 / K_BOND, 1.0 / lam.max()]
    slow = [1.0 / max(lam_ff.min(), 1e-30)]
    if model == "rolling":
        fast.append(TAU_ROT)
        slow.append(TAU_ROT)
    return min(fast), max(slow)


def run_one(model: str, omega: float, *, n_cycles=4, samples_per_cycle=48,
            dt_div=200, settle_taus=15):
    """★ 변형률 제어로 K* 를 잰다 — 중앙 비드의 y 를 강제하고 **그 비드가 받는 힘**을 잰다.

        K* = −F̂ / ŷ      (F 는 시료 힘만: 방사 본드 + 굽힘. 용매 항력은 안 들어감)

    트랩 컴플라이언스가 없어 추종이 정의상 100% 이고 조건수가 ω 에 무관하다.
    """
    tau_fast, tau_slow = timescales(model)
    period = 2 * math.pi / omega
    dt = min(tau_fast / dt_div, period / 2000)
    settle_steps = int(max(settle_taus * tau_slow, 2 * period) / dt)
    n_meas = int(n_cycles * period / dt)
    sim, mid, forces = build(model, omega, dt)
    t0 = time.time()
    sim.run(settle_steps)
    ts, ym, fm = [], [], []
    n_chunk = max(1, n_meas // (n_cycles * samples_per_cycle))
    done = 0
    while done < n_meas:
        sim.run(min(n_chunk, n_meas - done))
        done += min(n_chunk, n_meas - done)
        p = np.array(sim.state.get_snapshot().particles.position)
        ts.append((settle_steps + done) * dt)
        ym.append(float(p[mid, 1]))
        fm.append(float(sum(np.array(f.forces)[mid, 1] for f in forces)))
    wall = time.time() - t0
    ts, ym, fm = np.array(ts), np.array(ym), np.array(fm)
    nb = min(8, max(2, len(ts) // 12))
    by = LI.lockin_blocks(ts, ym, omega, n_blocks=nb)
    bf = LI.lockin_blocks(ts, fm, omega, n_blocks=nb)
    Kb = -bf / by
    K, Ksem = LI.agg(Kb)
    yh, _ = LI.agg(by)
    return dict(K_re=K.real, K_im=K.imag, K_sem=Ksem, follow=abs(yh) / AMP,
                y_hat=abs(yh), steps=settle_steps + n_meas, dt=dt, wall=wall,
                rate=(settle_steps + n_meas) / max(wall, 1e-9))


# ══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    print("=" * 88)
    print("조화 굽힘(JKR 선형화) vs 구름 저항 — kT=0 결정론적 ω 스윕")
    print("=" * 88)
    print(f"  n = {N_BEADS},  κ_θ = {KAPPA_THETA:.6g} kT,  k_r = 2κ_θ/R² = {K_ROLL:.6g}")
    print(f"  γ_r = 4a²/3 = {GAMMA_R:.6f}   →   τ_rot = γ_r/(k_rR²) = {TAU_ROT:.4e}")
    print(f"  ★ 예측 교차 주파수  ω_c = 1/τ_rot = {1/TAU_ROT:.4e}  (무차원)")
    print(f"     생산 런의 ω* = 18453  →  ω τ_rot = {18453*TAU_ROT:.3e}  "
          f"(≪1 이면 두 모델이 같아야 한다)")
    print()

    # 정적 극한의 절대 기준값 (MD 가 이걸 재현해야 한다)
    A = bending_matrix(N_BEADS, KAPPA_THETA, ELL)
    mid = N_BEADS // 2
    fx, fr = [0, mid, N_BEADS - 1], [i for i in range(N_BEADS) if i not in (0, mid, N_BEADS - 1)]
    yfx = np.array([0.0, 1.0, 0.0])
    yv = np.zeros(N_BEADS); yv[fx] = yfx
    yv[fr] = np.linalg.solve(A[np.ix_(fr, fr)], -A[np.ix_(fr, fx)] @ yfx)
    K_STATIC = float(yv @ A @ yv)
    print(f"  ★ 정적 극한 기준값 (선형응답 정확해)  K′(ω→0) = {K_STATIC:.6g} kT/d²")
    print()

    # ω=1e4·1e5 는 뺐다 — 주기가 길어 dt(빠른 모드가 정함) 대비 스텝이 70분/점이고,
    # 준정적 극한은 이미 위의 해석적 K_STATIC 으로 정확히 갖고 있다.
    omegas = [3e5, 1e6, 3e6, 1e7, 3e7, 1e8]
    rows = []
    print(f"  {'ω*':>10} {'ωτ_rot':>10} | {'조화 K′':>13} {'구름 K′':>13} {'구름/조화':>10} "
          f"| {'조화/정적':>9} | {'steps/s':>8}")
    print("  " + "-" * 96)
    for om in omegas:
        h = run_one("harmonic", om)
        r = run_one("rolling", om)
        rows.append(dict(omega=om, harmonic=h, rolling=r, K_static=K_STATIC))
        print(f"  {om:10.3g} {om*TAU_ROT:10.3e} | {h['K_re']:13.6g} {r['K_re']:13.6g} "
              f"{r['K_re']/h['K_re']:10.5f} | {h['K_re']/K_STATIC:9.5f} "
              f"| {r['rate']:8.0f}")
        (OUT / "sweep.json").write_text(json.dumps(rows, indent=1, default=float))

    print()
    print("  ★ dt 수렴 확인 (ω*=1e7, dt 절반):")
    for m in ("harmonic", "rolling"):
        a = run_one(m, 1e7, dt_div=200)
        b = run_one(m, 1e7, dt_div=400)
        print(f"    {m:9s}  K′ {a['K_re']:.6g} → {b['K_re']:.6g}   "
              f"변화 {100*(b['K_re']/a['K_re']-1):+.3f}%")
    print()
    print(f"  → {OUT/'sweep.json'}")
