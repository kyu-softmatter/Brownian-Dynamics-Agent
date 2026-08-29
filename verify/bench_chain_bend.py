"""`chain-bend-2d-oscill` 의 실제 스텝률을 잰다 — L4 비용 판단의 근거.

L3 스펙은 최저 ω 에서 2.65e9 스텝을 요구한다. "비싸다"는 말만으로는 결정할 수 없어서
초당 스텝을 실측한다. 세 구성을 비교한다:

  A  bond + angle 만                          — 컴파일된 힘만. 하한(가능한 최고 속도)
  B  + md.force.Custom 트랩 (매 스텝 파이썬)   — 현재 설계 (trap-2d/trap-drag 방식)
  C  + 유령입자에 bond.Harmonic(r0=0) 트랩     — 컴파일 경로. 구동 앵커만 updater 로 이동

C의 착안: bond.Harmonic(r0=0) 은 U = ½k r² 로 조화 트랩과 정확히 같다. 유령입자를
적분기 filter 에서 빼면 움직이지 않으므로 고정 트랩이 된다. 구동 트랩은 유령을
CustomUpdater 로 옮기는데, ω dt = 2.4e-7 이라 100스텝마다 옮겨도 위상 오차가
2.4e-5 주기다 — 매 스텝 파이썬을 호출할 이유가 없다.

파라미터는 specs/chain-bend-2d-oscill__w85__*.json (최저 ω, 가장 비싼 점) 에서 읽는다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/bench_chain_bend.py
"""
from __future__ import annotations

import glob
import json
import math
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH_STEPS = 200_000          # 실측 구간. 워밍업 후
WARMUP = 5_000
UPDATE_EVERY = 100             # C 구성에서 구동 앵커를 옮기는 주기


def load_spec() -> dict:
    """최저 ω 스펙 = 가장 비싼 점."""
    cands = sorted(glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__w85__*.json")))
    if not cands:
        raise SystemExit("specs/chain-bend-2d-oscill__w85__*.json 이 없습니다")
    return json.loads(Path(cands[0]).read_text())


def build_frame(n: int, L_chain: float, *, ghosts: list[int] | None = None):
    """직선 사슬 (x축). ghosts 가 주어지면 그 비드마다 유령입자를 하나 더 붙인다."""
    ell = L_chain / (n - 1)
    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    types = ["A"]
    typeid = [0] * n
    n_real = n
    if ghosts:
        types.append("G")
        for g in ghosts:
            pos.append(list(pos[g]))        # 유령은 초기 위치에 겹쳐 둔다 (Δr=0)
            typeid.append(1)

    box_L = 4.0 * L_chain                   # 사슬보다 훨씬 크게 → 래핑이 물리에 안 닿는다
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = types
    f.configuration.box = [box_L, box_L, 0, 0, 0, 0]     # Lz=0 → 2D (함정 9)
    f.configuration.dimensions = 2

    f.bonds.N = n - 1
    f.bonds.types = ["backbone"]
    f.bonds.typeid = [0] * (n - 1)
    bond_group = [[i, i + 1] for i in range(n - 1)]
    if ghosts:
        f.bonds.types = ["backbone", "trap"]
        for j, g in enumerate(ghosts):
            bond_group.append([g, n_real + j])
        f.bonds.N = len(bond_group)
        f.bonds.typeid = [0] * (n - 1) + [1] * len(ghosts)
    f.bonds.group = np.array(bond_group)

    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])
    return f, n_real, box_L, ell


class CustomTrap(md.force.Custom):
    """구성 B — 매 스텝 파이썬. trap-2d-5um / trap-drag 와 같은 방식."""

    def __init__(self, k, trapped, anchors, amp, omega, dt, drive_row):
        super().__init__(aniso=False)
        self.k = float(k)
        self.trapped = np.asarray(trapped, dtype=int)
        self.anchors = np.asarray(anchors, dtype=float)
        self.amp, self.omega, self.dt = float(amp), float(omega), float(dt)
        self.drive_row = int(drive_row)

    def set_forces(self, timestep):
        anc = self.anchors.copy()
        anc[self.drive_row, 1] += self.amp * math.sin(self.omega * timestep * self.dt)
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)
            pos = np.array(snap.particles.position, copy=True)
            arr.force[:] = 0.0
            arr.potential_energy[:] = 0.0
            for row, tg in enumerate(self.trapped):     # tag 인덱싱 필수
                loc = np.flatnonzero(tags == tg)
                d = pos[loc] - anc[row]
                arr.force[loc] = -self.k * d
                arr.potential_energy[loc] = 0.5 * self.k * (d ** 2).sum(axis=1)


class MoveGhost(hoomd.custom.Action):
    """구성 C — 구동 유령입자만 옮긴다. UPDATE_EVERY 스텝마다."""

    def __init__(self, ghost_tag, y0, amp, omega, dt):
        self.ghost_tag = int(ghost_tag)
        self.y0, self.amp = float(y0), float(amp)
        self.omega, self.dt = float(omega), float(dt)

    def act(self, timestep):
        y = self.y0 + self.amp * math.sin(self.omega * timestep * self.dt)
        with self._state.cpu_local_snapshot as snap:
            tags = np.array(snap.particles.tag, copy=True)
            loc = np.flatnonzero(tags == self.ghost_tag)
            if len(loc):
                snap.particles.position[loc[0]][1] = y


def make_sim(variant: str, p: dict, nu: dict):
    n = int(p["n_beads"])
    trapped = [int(t) for t in p["trapped"]]
    k_t, k_b = float(p["k_t_star"]), float(p["k_bond_star"])
    kappa = float(p["kappa_theta_star"])
    amp, omega = float(p["amp_star"]), float(p["omega_star"])
    dt = float(nu["dt_star"])
    ghosts = trapped if variant == "C" else None

    f, n_real, box_L, ell = build_frame(n, float(p["L_chain_star"]), ghosts=ghosts)
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)

    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=k_b, r0=ell)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=kappa, t0=math.pi)
    forces = [bond, angle]

    if variant == "C":
        # bond.Harmonic(r0=0) = ½k r²  → 조화 트랩과 동일. 컴파일 경로.
        bond.params["trap"] = dict(k=k_t, r0=0.0)
        integrated = hoomd.filter.Type(["A"])        # 유령("G")은 적분 안 함 → 고정
    else:
        integrated = hoomd.filter.All()
        if variant == "B":
            anchors = np.array([[(t - (n - 1) / 2) * ell, 0.0, 0.0] for t in trapped])
            drive_row = trapped.index(sorted(trapped)[len(trapped) // 2])
            forces.append(CustomTrap(k_t, trapped, anchors, amp, omega, dt, drive_row))

    bd = md.methods.Brownian(filter=integrated, kT=1.0, default_gamma=1.0)
    integrator = md.Integrator(dt=dt, methods=[bd], forces=forces)
    integrator.integrate_rotational_dof = False      # BD 는 과감쇠 (함정 5)
    sim.operations.integrator = integrator

    if variant == "C":
        mid = sorted(trapped)[len(trapped) // 2]
        ghost_tag = n_real + trapped.index(mid)
        sim.operations.updaters.append(hoomd.update.CustomUpdater(
            action=MoveGhost(ghost_tag, 0.0, amp, omega, dt),
            trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    return sim


def main() -> int:
    spec = load_spec()
    p, nu = spec["params"], spec["numerics"]
    total = int(nu["n_prod"]) + int(nu["n_eq"])
    print("=" * 78)
    print("chain-bend-2d-oscill — 스텝률 실측 (최저 ω, 가장 비싼 점)")
    print("=" * 78)
    print(f"n_beads={p['n_beads']}  dt*={nu['dt_star']:.3e}  "
          f"요구 스텝(이 ω) = {total:,}")
    print(f"실측 구간 {BENCH_STEPS:,} 스텝 (워밍업 {WARMUP:,})\n")
    print(f"{'구성':<44}{'steps/s':>12}{'이 ω 소요':>14}")
    print("-" * 78)

    labels = {
        "A": "A  bond + angle 만 (컴파일 힘만)",
        "B": "B  + force.Custom 트랩 (매 스텝 파이썬)",
        "C": f"C  + 유령 bond 트랩 (updater {UPDATE_EVERY}스텝)",
    }
    rates = {}
    for v in ("A", "B", "C"):
        sim = make_sim(v, p, nu)
        sim.run(WARMUP)
        t0 = time.perf_counter()
        sim.run(BENCH_STEPS)
        el = time.perf_counter() - t0
        rate = BENCH_STEPS / el
        rates[v] = rate
        days = total / rate / 86400
        span = f"{days:.1f} 일" if days >= 1 else f"{total / rate / 3600:.1f} 시간"
        print(f"{labels[v]:<44}{rate:>12,.0f}{span:>14}")

    print("-" * 78)
    print(f"파이썬 트랩의 대가  B/A = {rates['B'] / rates['A']:.3f}  "
          f"(느려짐 {rates['A'] / rates['B']:.1f}×)")
    print(f"유령 트랩의 회수    C/B = {rates['C'] / rates['B']:.2f}× 빠름")

    # 스윕 전체 비용. ★ 파일명 알파벳 정렬은 ω 순서가 아니다 (w1737 < w85) — 값으로 정렬한다
    allspecs = [json.loads(Path(s).read_text())
                for s in glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__*.json"))]
    steps = sorted(s["numerics"]["n_prod"] + s["numerics"]["n_eq"] for s in allspecs)
    tot, worst = sum(steps), steps[-1]      # 최저 ω = 최다 스텝 = 병렬 시 벽시계
    print(f"\n스윕 {len(allspecs)}점 합계 = {tot:,} 스텝 "
          f"(최다 = {worst:,})")
    for v in ("B", "C"):
        ser, par = tot / rates[v] / 86400, worst / rates[v] / 3600
        print(f"  구성 {v}: 직렬 {ser:6.1f} 일   7런 병렬 → 벽시계 {par:6.1f} 시간")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
