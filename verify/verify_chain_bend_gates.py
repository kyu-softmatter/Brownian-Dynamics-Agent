"""`chain-bend-2d-oscill` 실행 전 관문 2종. 둘 다 프로덕션 스윕 **전에** 통과해야 한다.

────────────────────────────────────────────────────────────────────────────
관문 A (`--gate lockin`) — 관측 코드를 해석해에 대조 (mater_plan 원칙 9 / 규칙 7)
────────────────────────────────────────────────────────────────────────────
이 케이스의 산출물 전체가 **새로 쓰는 락인 추출 코드**에 얹혀 있다. 사슬을 빼고
비드 1개 + 구동 트랩(k_t) + 정적 유령 스프링(k_s) 으로 바꾸면 답이 닫힌 형태로 나온다:

    ŷ = k_t a / (k_t + k_s + iωγ)          (비드 응답)
    K*_sample = k_s  (실수, ω 무관)         ← 추정량이 돌려줘야 하는 값

★ 이 관문이 1차 시도에서 **FAIL 했고**, 그게 이 관문을 만든 이유다. 원인은 물리가
아니라 구동이었다 — 유령을 U 스텝마다 옮기면 구동이 **영차 유지(zero-order hold)**
가 되어 기본파가 sinc(ωΔt/2) 로 줄고 위상이 ωΔt/2 만큼 늦는다 (Δt = U·dt).
실측 De=10 에서 |ŷ_c|/a = 0.98999, 위상 −0.2522 rad / ZOH 예측 0.99040, −0.2404 rad.

교훈은 **공칭 진폭 a 를 추정량에 쓰지 말라**는 것이다. 유령 위치를 같이 재서
**측정된 위상자 ŷ_c** 를 쓰면 ZOH 감쇠가 분자·분모에서 정확히 상쇠된다 — 비드는
공칭 사인이 아니라 유령이 실제로 있는 곳에 반응하기 때문이다. 아래 표는 두 추정량을
나란히 찍어 그것을 보인다 (공칭은 De 와 함께 무너지고 측정은 평평하다).

  ① 유령 위치 락인 → ZOH 예측과 일치하는가        [구동을 정량적으로 이해했는가]
  ② 비드 위치 락인 → 해석해 ŷ                    [BD + 트랩 + 락인]
  ③ K* 추정량 (측정 ŷ_c) → (k_s, 0)              [추정량 전체]
  ④ 통계 10배 → 오차가 1/√N 로 줄어드는가        [편향과 잡음의 분리]
①이 없으면 조용히 틀린다 — 유령이 안 움직여도 런은 완주하고 K* 만 엉뚱해진다.

생산 런의 dt(4.53e-10) 에서는 U=100 이어도 ωΔt = 2.19e-3 rad (De=10) 이라 ZOH 자체가
무해하다. 관문 A 는 dt 가 220배 커서 효과가 증폭된 것이고, 그래서 **추정량의 취약점을
드러내는 데 오히려 유용한 조건**이다 (함정 1 의 교훈 — 약한 조건으로 검증한다).

────────────────────────────────────────────────────────────────────────────
관문 B (`--gate inertia`) — τ_p/τ_fast = 0.60 을 논증하지 말고 측정 (규칙 6)
────────────────────────────────────────────────────────────────────────────
최속 굽힘 모드는 과감쇠가 아니다 (감쇠비 ζ = γ/2√(mλ_max) = 0.65 < 1). BD 는 그
모드를 과감쇠로 강제하므로 그 대역의 동역학이 틀린다 — 어떤 dt 로도 안 고쳐진다.
"관측 대역과 4570배 떨어져 있어 영향 없을 것"은 지금까지 **확인되지 않은 추론**이다.

같은 파라미터로 Brownian(관성 없음) vs Langevin(관성 있음) 을 돌려 측정 대역의
K*(ω) 를 비교한다. 일치하면 그 추론이 **측정된 여유**가 된다. 오염이 가장 클 곳은
τ_fast 에 가장 가까운 **최고 ω** 이므로 De = 10 과 4.7 에서 본다.

────────────────────────────────────────────────────────────────────────────
    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_chain_bend_gates.py --gate lockin
    $PY scratch/verify_chain_bend_gates.py --gate inertia --method bd       --de 10
    $PY scratch/verify_chain_bend_gates.py --gate inertia --method langevin --de 10
    $PY scratch/verify_chain_bend_gates.py --gate inertia --collect
"""
from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
import sys as _sys; _sys.path.insert(0, str(ROOT))
OUT = ROOT / "scratch" / "_gates"
UPDATE_EVERY = 100          # 구동 유령을 옮기는 주기. ZOH 를 남겨 두고 추정량으로 상쇠한다
SAMPLES_PER_CYCLE = 20      # 스펙과 같은 표본 밀도 (2000표본/100주기)
N_SIGMA = 3.0               # 해석해와 "구분 안 됨" 판정 기준
KAPPA_CENTER = 2.0 * 2415.33   # κ_center* = 2 κ_end* (스펙 원장 kappa_end_d2/kT)
TAU_P = 2.70624e-8          # 스펙 원장 tau_p/tau_B
TAU_CHAIN = 2.07011e-4      # 스펙 원장 tau_chain/tau_B


def load_specs() -> list[dict]:
    """ω 오름차순. ★ 파일명 알파벳 정렬은 ω 순서가 아니다 (w1737 < w85)."""
    specs = [json.loads(Path(p).read_text())
             for p in glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__*.json"))]
    if not specs:
        raise SystemExit("specs/chain-bend-2d-oscill__*.json 이 없습니다")
    return sorted(specs, key=lambda s: s["params"]["omega_star"])


# ════════════════════════════════════════════════════════════════════════
# 락인 — 블록별 위상자를 돌려준다 (오차막대를 블록 산포에서 만든다)
# ════════════════════════════════════════════════════════════════════════
# ★ 락인 추정량은 `bdbot/lockin.py` 로 올렸습니다 (관문 A + 생산 런에서 두 번 씀).
#   본문은 이 파일에서 검증된 것을 그대로 옮긴 것이고, 수치 동일성을 대조했습니다.
from bdbot.lockin import agg, k_star, lockin_blocks  # noqa: E402


# ════════════════════════════════════════════════════════════════════════
# 표본 수집 + 구동
# ════════════════════════════════════════════════════════════════════════
class Sampler(hoomd.custom.Action):
    """비드와 구동 유령의 y 를 함께 남긴다 — 유령까지 재야 구동을 검사·보정할 수 있다."""

    def __init__(self, bead_tag: int, ghost_tag: int, dt: float):
        self.bead_tag, self.ghost_tag, self.dt = int(bead_tag), int(ghost_tag), float(dt)
        self.t: list[float] = []
        self.y_bead: list[float] = []
        self.y_ghost: list[float] = []

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            tags = np.array(snap.particles.tag, copy=True)
            pos = np.array(snap.particles.position, copy=True)
            self.t.append(timestep * self.dt)
            self.y_bead.append(float(pos[np.flatnonzero(tags == self.bead_tag)[0], 1]))
            self.y_ghost.append(float(pos[np.flatnonzero(tags == self.ghost_tag)[0], 1]))


class MoveGhost(hoomd.custom.Action):
    """구동 유령을 y = y0 + a sin(ωt) 로 옮긴다. 컴파일 경로를 비워 두는 것이 목적."""

    def __init__(self, ghost_tag: int, y0: float, amp: float, omega: float, dt: float):
        self.ghost_tag, self.y0 = int(ghost_tag), float(y0)
        self.amp, self.omega, self.dt = float(amp), float(omega), float(dt)

    def act(self, timestep):
        y = self.y0 + self.amp * math.sin(self.omega * timestep * self.dt)
        with self._state.cpu_local_snapshot as snap:
            tags = np.array(snap.particles.tag, copy=True)
            loc = np.flatnonzero(tags == self.ghost_tag)
            if len(loc):
                snap.particles.position[loc[0], 1] = y      # 2-인덱스로 확실히 써넣는다


def attach(sim, dt, bead_tag, ghost_tag, amp, omega, sample_every):
    sampler = Sampler(bead_tag, ghost_tag, dt)
    sim.operations.writers.append(hoomd.write.CustomWriter(
        action=sampler, trigger=hoomd.trigger.Periodic(int(sample_every))))
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=MoveGhost(ghost_tag, 0.0, amp, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    return sampler


def run_window(sim, dt, bead_tag, ghost_tag, amp, omega, n_cycles, n_eq=0):
    """구동을 켜고 n_eq 스텝 지난 뒤 n_cycles 주기를 수집한다."""
    period = 2.0 * math.pi / omega
    spc = max(SAMPLES_PER_CYCLE, int(round(period / dt)))
    sample_every = max(1, spc // SAMPLES_PER_CYCLE)
    smp = attach(sim, dt, bead_tag, ghost_tag, amp, omega, sample_every)
    if n_eq:
        sim.run(int(n_eq))
        smp.t.clear(); smp.y_bead.clear(); smp.y_ghost.clear()
    sim.run(int(round(n_cycles * period / dt)))
    return (np.array(smp.t), np.array(smp.y_bead), np.array(smp.y_ghost))


# ════════════════════════════════════════════════════════════════════════
# 관문 A — 비드 1개 + 구동 트랩 + 정적 유령 스프링
# ════════════════════════════════════════════════════════════════════════
def build_single(k_t: float, k_s: float, dt: float):
    """tag 0 = 비드, 1 = 구동 유령(k_t), 2 = 정적 유령(k_s). 유령은 적분하지 않는다."""
    f = gsd.hoomd.Frame()
    f.particles.N = 3
    f.particles.position = np.zeros((3, 3))
    f.particles.typeid = [0, 1, 1]
    f.particles.types = ["A", "G"]
    f.configuration.box = [40.0, 40.0, 0, 0, 0, 0]          # Lz=0 → 2D (함정 9)
    f.configuration.dimensions = 2
    f.bonds.N = 2
    f.bonds.types = ["trap", "spring"]
    f.bonds.typeid = [0, 1]
    f.bonds.group = np.array([[0, 1], [0, 2]])
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=3)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["trap"] = dict(k=k_t, r0=0.0)               # U = ½k r² = 조화 트랩
    bond.params["spring"] = dict(k=k_s, r0=0.0)
    bd = md.methods.Brownian(filter=hoomd.filter.Type(["A"]), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[bond])
    integ.integrate_rotational_dof = False                  # BD 는 과감쇠 (함정 5)
    return sim, integ


def measure_single(k_t, k_s, amp, omega, dt, n_cycles, n_blocks=10):
    sim, integ = build_single(k_t, k_s, dt)
    sim.operations.integrator = integ
    t, yb, yg = run_window(sim, dt, 0, 1, amp, omega, n_cycles,
                           n_eq=int(round(0.1 * n_cycles * 2 * math.pi / omega / dt)))
    yb_b = lockin_blocks(t, yb, omega, n_blocks=n_blocks)
    yg_b = lockin_blocks(t, yg, omega, n_blocks=n_blocks)
    h3_b = lockin_blocks(t, yb, omega, harmonic=3, n_blocks=n_blocks)
    K_b = np.array([k_star(y, g, k_t, omega) for y, g in zip(yb_b, yg_b)])
    Kn_b = np.array([k_star(y, complex(amp, 0.0), k_t, omega) for y in yb_b])
    return dict(y=agg(yb_b), g=agg(yg_b), K=agg(K_b), Kn=agg(Kn_b),
                h3=abs(agg(h3_b)[0]), n=len(t))


def gate_lockin() -> int:
    specs = load_specs()
    k_t = float(specs[0]["params"]["k_t_star"])
    amp = float(specs[0]["params"]["amp_star"])
    k_s, gamma, n_cycles = KAPPA_CENTER, 1.0, 100
    dt = 1e-3 * gamma / (k_t + k_s)                          # BD 는 O(δt) — 넉넉히 (함정 2)

    print("=" * 104)
    print("관문 A — 락인 + K* 추정량을 해석해에 대조 (비드 1개 + 구동 트랩 + 정적 스프링)")
    print("=" * 104)
    print(f"k_t = {k_t:.2f}   k_s = κ_center* = {k_s:.2f}   a = {amp:.5f}   "
          f"γ = {gamma}   dt = {dt:.3e}   {n_cycles}주기")
    print(f"기대: K* = ({k_s:.2f}, 0)  실수, ω 무관.   판정 = 해석해와 {N_SIGMA:.0f}σ 안\n")
    print(f"{'De':>6} │{'|ŷ_c|/a':>9}{'ZOH예측':>9} │{'ŷ 오차%':>9} │"
          f"{'K′(공칭a)':>11}{'오차%':>8} │{'K′(측정ŷ_c)':>13}{'±σ':>8}{'오차%':>8}"
          f"{'K″':>9}{'±σ':>7} │{'3차/1차':>9} {'':>3}")
    print("-" * 104)

    fails, worst = [], 0.0
    for sp in specs:
        omega, de = float(sp["params"]["omega_star"]), float(sp["params"]["De"])
        m = measure_single(k_t, k_s, amp, omega, dt, n_cycles)
        (y_hat, _), (g_hat, _), (K, Ks), (Kn, _) = m["y"], m["g"], m["K"], m["Kn"]

        x = omega * dt * UPDATE_EVERY / 2.0
        zoh = math.sin(x) / x if x else 1.0
        y_ex = k_t * g_hat / complex(k_t + k_s, omega * gamma)   # 측정 구동 기준 해석해
        err_y = 100.0 * abs(y_hat - y_ex) / abs(y_ex)
        err_K = 100.0 * abs(K.real - k_s) / k_s
        err_Kn = 100.0 * abs(Kn.real - k_s) / k_s
        ok = (abs(K.real - k_s) <= N_SIGMA * Ks) and (abs(K.imag) <= N_SIGMA * Ks)
        worst = max(worst, err_K)
        if not ok:
            fails.append(round(de, 2))
        print(f"{de:>6.2f} │{abs(g_hat) / amp:>9.5f}{zoh:>9.5f} │{err_y:>9.3f} │"
              f"{Kn.real:>11.1f}{err_Kn:>8.2f} │{K.real:>13.1f}{Ks:>8.1f}{err_K:>8.2f}"
              f"{K.imag:>9.1f}{Ks:>7.1f} │{m['h3'] / abs(y_hat):>9.1e} "
              f"{'✓' if ok else '✗':>3}")

    print("-" * 104)
    print(f"측정 ŷ_c 추정량 최대 오차 {worst:.2f}%   "
          f"{'✓ 모든 ω 에서 해석해와 3σ 안' if not fails else f'✗ FAIL at De={fails}'}")

    # ④ 편향 vs 잡음 — 통계를 10배로 늘려 오차가 1/√N 로 줄면 편향이 아니다
    print("\n④ 편향 검사 — 통계 10배 (De≈1). 편향이 없으면 오차와 σ 가 함께 √10 배 줄어든다")
    sp1 = min(specs, key=lambda s: abs(s["params"]["De"] - 1.0))
    om1 = float(sp1["params"]["omega_star"])
    print(f"{'주기수':>8}{'K′':>12}{'±σ':>9}{'|K′−k_s|':>11}{'그 배수 σ':>11}{'K″':>10}")
    prev = None
    for nc, nb in ((100, 10), (1000, 10)):
        m = measure_single(k_t, k_s, amp, om1, dt, nc, n_blocks=nb)
        K, Ks = m["K"]
        d = abs(K.real - k_s)
        print(f"{nc:>8}{K.real:>12.1f}{Ks:>9.1f}{d:>11.1f}{d / Ks:>11.2f}{K.imag:>10.1f}")
        prev = (d, Ks) if prev is None else prev
    print(f"  (기대: 100→1000 주기에서 σ 가 약 √10 = 3.16 배 감소)")
    print("=" * 104)
    return 0 if not fails else 1


# ════════════════════════════════════════════════════════════════════════
# 관문 B — 사슬 전체, Brownian vs Langevin
# ════════════════════════════════════════════════════════════════════════
def build_chain(sp: dict, mass: float):
    p = sp["params"]
    n = int(p["n_beads"])
    trapped = sorted(int(t) for t in p["trapped"])
    ell = float(p["L_chain_star"]) / (n - 1)
    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:                                        # 유령을 비드 위에 겹쳐 둔다
        pos.append(list(pos[g]))
        typeid.append(1)

    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.particles.mass = [mass] * n + [1.0] * len(trapped)
    f.configuration.box = [4.0 * float(p["L_chain_star"])] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2

    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(trapped)]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(trapped)
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=5)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=float(p["k_bond_star"]), r0=ell)
    bond.params["trap"] = dict(k=float(p["k_t_star"]), r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=float(p["kappa_theta_star"]), t0=math.pi)
    mid = trapped[len(trapped) // 2]
    return sim, [bond, angle], mid, n + trapped.index(mid)


def gate_inertia_one(de_target: float, method: str, n_cycles: int, n_eq_tau: float,
                     eq_steps: int = 0) -> dict:
    specs = load_specs()
    sp = min(specs, key=lambda s: abs(s["params"]["De"] - de_target))
    p, nu = sp["params"], sp["numerics"]
    omega, amp = float(p["omega_star"]), float(p["amp_star"])
    k_t, dt, gamma = float(p["k_t_star"]), float(nu["dt_star"]), 1.0
    # 관성이 있는 쪽만 질량을 준다. kT=0 인 두 방법은 **결정론적** — 통계오차가 0 이라
    # 관성 항만 따로 떼어 볼 수 있다 (열적 비교는 검정력이 없었다).
    inertial = method in ("langevin", "lang0")
    mass = TAU_P * gamma if inertial else 0.0

    sim, forces, bead_tag, ghost_tag = build_chain(sp, mass if mass else 1.0)
    filt = hoomd.filter.Type(["A"])
    integ_m = {
        "bd": lambda: md.methods.Brownian(filter=filt, kT=1.0, default_gamma=gamma),
        "langevin": lambda: md.methods.Langevin(filter=filt, kT=1.0, default_gamma=gamma),
        # 열잡음 없는 과감쇠 — skill bd-hoomd 가 결정론적 검증용으로 지정한 적분기
        "ov": lambda: md.methods.OverdampedViscous(filter=filt, default_gamma=gamma),
        "lang0": lambda: md.methods.Langevin(filter=filt, kT=0.0, default_gamma=gamma),
    }[method]()
    integ = md.Integrator(dt=dt, methods=[integ_m], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    if method == "langevin":
        sim.state.thermalize_particle_momenta(filter=filt, kT=1.0)

    # ★ TAU_CHAIN(= γ/κ_center) 은 이 계의 최장 이완시간이 **아니다** — (A+T) 의 최저
    # 고유값에서 나오는 τ_max = γ/λ_min 이 9.18배 길다. 그것으로 평형화해야 한다.
    n_eq = int(eq_steps) if eq_steps else int(round(n_eq_tau * TAU_CHAIN / dt))
    t, yb, yg = run_window(sim, dt, bead_tag, ghost_tag, amp, omega, n_cycles, n_eq=n_eq)

    # ★ 블록은 **최소 한 주기**를 담아야 한다. 반주기 블록으로 락인하면 블록값이
    # 무의미해지고 SEM 이 물리와 무관한 수가 된다 (전체창 추정치는 영향 없음).
    nb = max(1, min(6, n_cycles))
    yb_b = lockin_blocks(t, yb, omega, n_blocks=nb)
    yg_b = lockin_blocks(t, yg, omega, n_blocks=nb)
    h3_b = lockin_blocks(t, yb, omega, harmonic=3, n_blocks=nb)
    K_b = np.array([k_star(y, g, k_t, omega, gamma, mass) for y, g in zip(yb_b, yg_b)])
    (y_hat, y_sem), (g_hat, _), (K, K_sem) = agg(yb_b), agg(yg_b), agg(K_b)

    return dict(method=method, de=float(p["De"]), omega=omega, mass=mass,
                n_eq=n_eq, n_cycles=n_cycles, n_samples=len(t),
                drive_abs=abs(g_hat) / amp, y_abs=abs(y_hat), y_sem=y_sem,
                K_re=K.real, K_im=K.imag, K_sem=K_sem,
                h3_rel=abs(agg(h3_b)[0]) / abs(y_hat))


def gate_inertia_collect() -> int:
    rows = [json.loads(f.read_text()) for f in sorted(OUT.glob("inertia_*.json"))]
    if not rows:
        raise SystemExit(f"{OUT}/inertia_*.json 이 없습니다 — 먼저 각 구성을 돌리세요")

    print("=" * 96)
    print("관문 B — τ_p/τ_fast = 0.60 의 오염을 측정 (Brownian vs Langevin, 같은 파라미터)")
    print("=" * 96)
    print(f"κ_center* = {KAPPA_CENTER:.1f}   ζ_fast = 0.65 (< 1 → 최속 모드는 실제로 저감쇠)")
    print("Langevin 은 관성이 있고 BD 는 없다. 측정 대역에서 K* 가 같으면 BD 로 재도 된다.\n")
    print(f"{'De':>6} {'방법':<10}{'|ŷ_c|/a':>10}{'K′':>11}{'±σ':>9}{'K″':>11}{'±σ':>9}"
          f"{'3차/1차':>10}{'표본':>7}")
    print("-" * 96)
    for r in sorted(rows, key=lambda x: (-x["de"], x["method"])):
        print(f"{r['de']:>6.2f} {r['method']:<10}{r['drive_abs']:>10.6f}"
              f"{r['K_re']:>11.1f}{r['K_sem']:>9.1f}{r['K_im']:>11.1f}{r['K_sem']:>9.1f}"
              f"{r['h3_rel']:>10.1e}{r['n_samples']:>7}")
    print("-" * 96)

    verdict = 0
    for de in sorted({r["de"] for r in rows}, reverse=True):
        pair = {r["method"]: r for r in rows if r["de"] == de}
        if len(pair) < 2:
            print(f"De={de:.2f}: 짝이 없어 비교 못 함 ({list(pair)})")
            verdict = 1
            continue
        b, l = pair["bd"], pair["langevin"]
        for lbl, kb, kl in (("K′", b["K_re"], l["K_re"]), ("K″", b["K_im"], l["K_im"])):
            d, sig = abs(kb - kl), math.hypot(b["K_sem"], l["K_sem"])
            nsig = d / sig if sig > 0 else float("inf")
            rel = 100.0 * d / max(abs(kb), abs(kl), 1e-30)
            ok = nsig <= N_SIGMA
            verdict |= 0 if ok else 1
            print(f"De={de:>5.2f} {lbl}: |BD−Langevin| = {d:9.2f} = {nsig:5.2f}σ "
                  f"({rel:6.2f}%)  {'✓ 구분 안 됨' if ok else '✗ 유의한 차이'}")
    print("-" * 96)
    print("✓ PASS — 관성이 측정 대역의 K* 를 바꾸지 않는다. BD 로 스윕해도 된다."
          if verdict == 0 else
          "✗ FAIL — 관성이 측정 대역에 들어온다. BD 로는 이 대역을 못 잰다.")
    print("=" * 96)
    return verdict


def gate_det_collect() -> int:
    """결정론적(kT=0) 비교 — 통계오차가 0 이므로 관성 항만 정확히 떼어 본다.

    부수 산출: 잡음 없는 K*(ω) 곡선. 저주파 극한이 κ_center 로 가는지가
    implementation_check 이고, 곡선의 **모양**이 hypothesis (단일 Maxwell vs 모드 스펙트럼).
    """
    rows = [json.loads(f.read_text()) for f in sorted(OUT.glob("det_*.json"))]
    if not rows:
        raise SystemExit(f"{OUT}/det_*.json 이 없습니다")

    print("=" * 100)
    print("관문 B′ — 결정론적 비교 (OverdampedViscous vs Langevin kT=0). 통계오차 0")
    print("=" * 100)
    print(f"κ_center* = {KAPPA_CENTER:.1f}   m* = {TAU_P:.3e}   ζ_fast = 0.65 (최속 모드는 저감쇠)")
    print("열적 비교는 |ŷ|/ℓ_k < 1 때문에 검정력이 없었다. kT=0 이면 그 문제가 사라진다.\n")
    print(f"{'De':>7} │{'K′(과감쇠)':>13}{'K′(관성)':>13}{'차이%':>8} │"
          f"{'K″(과감쇠)':>13}{'K″(관성)':>13}{'차이%':>8} │{'K′/κ_c':>8}")
    print("-" * 100)

    worst, pairs = 0.0, 0
    for de in sorted({r["de"] for r in rows}):
        pr = {r["method"]: r for r in rows if r["de"] == de}
        if len(pr) < 2:
            print(f"{de:>7.3f} │ 짝 없음 ({list(pr)})")
            continue
        o, l = pr["ov"], pr["lang0"]
        pairs += 1
        d_re = 100.0 * abs(o["K_re"] - l["K_re"]) / max(abs(o["K_re"]), 1e-30)
        d_im = 100.0 * abs(o["K_im"] - l["K_im"]) / max(abs(o["K_im"]), 1e-30)
        worst = max(worst, d_re, d_im)
        print(f"{de:>7.3f} │{o['K_re']:>13.1f}{l['K_re']:>13.1f}{d_re:>8.3f} │"
              f"{o['K_im']:>13.1f}{l['K_im']:>13.1f}{d_im:>8.3f} │"
              f"{o['K_re'] / KAPPA_CENTER:>8.3f}")
    print("-" * 100)
    ok = worst < 1.0
    print(f"짝 {pairs}개 · 관성 유무의 최대 차이 = {worst:.3f}%")
    print("✓ PASS — 관성 항이 측정 대역의 K*(ω) 를 1% 미만으로 바꾼다. BD 로 스윕해도 된다."
          if ok else
          f"✗ FAIL — 관성이 최대 {worst:.1f}% 바꾼다. BD 로는 이 대역을 못 잰다.")
    print("  (이 검정은 관성 항 자체만 본다. 열적 링잉 × 비선형 결합은 덮지 못한다 —")
    print("   max|θ| = 2.8e-3 rad 라 작을 것으로 보이지만 별도 확인이 필요하다.)")
    print("=" * 100)
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", choices=["lockin", "inertia", "det"], required=True)
    ap.add_argument("--method", choices=["bd", "langevin", "ov", "lang0"])
    ap.add_argument("--de", type=float, default=10.0)
    ap.add_argument("--cycles", type=int, default=60)
    ap.add_argument("--eq-tau", type=float, default=5.0)
    ap.add_argument("--eq-steps", type=int, default=0,
                    help="평형화 스텝을 절대값으로 지정 (τ_max 기준으로 줄 때)")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()

    if a.gate == "lockin":
        return gate_lockin()
    if a.collect:
        return gate_det_collect() if a.gate == "det" else gate_inertia_collect()
    if not a.method:
        raise SystemExit("--method bd|langevin|ov|lang0 이 필요합니다")
    OUT.mkdir(parents=True, exist_ok=True)
    res = gate_inertia_one(a.de, a.method, a.cycles, a.eq_tau, a.eq_steps)
    pre = ("deq" if a.eq_steps else "det") if a.gate == "det" else "inertia"
    (OUT / f"{pre}_de{res['de']:.3f}_{a.method}.json").write_text(json.dumps(res, indent=2))
    print(f"[de{res['de']:.2f} {a.method}] K* = ({res['K_re']:.1f}, {res['K_im']:.1f}) "
          f"± {res['K_sem']:.1f}   유령 {res['drive_abs']:.6f}   표본 {res['n_samples']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
