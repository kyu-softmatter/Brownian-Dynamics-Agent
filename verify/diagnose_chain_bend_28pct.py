"""chain-bend 의 **28% 미해명 불일치**를 격리로 좁힌다 (규칙 7).

배경: 소각 선형응답 (iωγI + A + T)ŷ = k_t a e_mid 로 계산한 K*(ω) 가 평형화를 충분히 준
HOOMD 결정론 실측과 고주파에서 최대 28% 어긋난다. 이미 확인된 것 —
  · 굽힘 행렬 ↔ HOOMD angle.Harmonic : 0.55% 일치 (정적, 트랩 없음, 강체 고정)
  · λ_max : 스펙과 일치
  · 유령 트랩 : 관문 A 에서 7개 ω 전부 해석해와 3σ 안
  · 과도 : 20 τ_max 로 수렴 (σ/K = 0.02%)
남은 후보를 사다리로 하나씩 떼어 본다.

  ① dt (명시적 오일러 이산화)  — **시뮬레이션 없이** z-영역 정확해로 판정
     명시적 오일러 y_{n+1} = (I − (A+T)dt/γ) y_n + (k_t dt/γ) u_n 의 정상응답은
     (e^{iωdt} I − M) ŷ = B â 로 **닫힌 형태**다. 연속 해와 비교하면 dt 기여가 정확히 나온다.
  ② 정적 + 트랩              — 트랩까지 켠 정적 강성을 HOOMD 로 직접 (kT=0, 큰 dt)
     ①에서 검증된 것은 '강체 고정, 트랩 없음' 이었다. 트랩이 들어오면 달라지는가?
  ③ 동적 모드 형태            — 25개 비드 전부의 복소 위상자 ŷ_i 를 재서 모델과 비교
     불일치가 **공간적으로 어디서** 생기는지 본다 (구동 비드만? 끝? 전체 스케일?)

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/diagnose_chain_bend_28pct.py --stage 12      # ①② (빠름)
    $PY scratch/diagnose_chain_bend_28pct.py --stage 3 --de 91.8
"""
from __future__ import annotations

import argparse
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
GATES = ROOT / "scratch" / "_gates"
sys.path.insert(0, str(ROOT))
from cases.chain_bend_2d import bending_matrix, sweep_specs, trapped_indices   # noqa: E402


def setup():
    sp = sweep_specs()                      # 스펙 파일 없이도 동작 (하드 검사 실패가 정상 상태)
    p = sp[0]["params"]
    n = int(p["n_beads"])
    kth = float(p["kappa_theta_star"]); kt = float(p["k_t_star"])
    kb = float(p["k_bond_star"])
    ell = float(p["L_chain_star"]) / (n - 1)
    amp = float(p["amp_star"]); dt = float(sp[0]["numerics"]["dt_star"])
    idx = trapped_indices(n); mid = idx[len(idx) // 2]
    A = bending_matrix(n, kth, ell)
    T = np.zeros((n, n))
    for e in idx:
        T[e, e] += kt
    return dict(specs=sp, n=n, kth=kth, kt=kt, kb=kb, ell=ell, amp=amp, dt=dt,
                idx=idx, mid=mid, A=A, T=T, AT=A + T)


def K_continuous(S, omega, gamma=1.0):
    n, mid = S["n"], S["mid"]
    y = np.linalg.solve(1j * omega * gamma * np.eye(n) + S["AT"], S["kt"] * np.eye(n)[mid])
    return S["kt"] / y[mid] - S["kt"] - 1j * omega * gamma, y


def K_euler(S, omega, dt, gamma=1.0):
    """명시적 오일러의 **정확한** 정상응답. y_n = Im[ŷ e^{iωn dt}] 를 넣으면
    (e^{iωdt} I − M) ŷ = B â,  M = I − (A+T)dt/γ,  B = (k_t dt/γ) e_mid."""
    n, mid = S["n"], S["mid"]
    M = np.eye(n) - S["AT"] * dt / gamma
    B = (S["kt"] * dt / gamma) * np.eye(n)[mid]
    y = np.linalg.solve(math.e ** (1j * omega * dt) * np.eye(n) - M, B)
    return S["kt"] / y[mid] - S["kt"] - 1j * omega * gamma, y


def stage1(S):
    print("=" * 92)
    print("① dt (명시적 오일러) — 시뮬레이션 없이 z-영역 정확해로 판정")
    print("=" * 92)
    print(f"dt = {S['dt']:.4e}   dt·λ_max/γ = {S['dt']*np.linalg.eigvalsh(S['AT'])[-1]:.4f}")
    print(f"{'De':>8}{'K′ 연속':>12}{'K′ 오일러':>12}{'차이%':>8}"
          f"{'K″ 연속':>12}{'K″ 오일러':>12}{'차이%':>8}")
    print("-" * 92)
    worst = 0.0
    for sp in S["specs"]:
        om = float(sp["params"]["omega_star"]); de = float(sp["params"]["De"])
        Kc, _ = K_continuous(S, om)
        Ke, _ = K_euler(S, om, S["dt"])
        dre = 100 * (Ke.real / Kc.real - 1); dim = 100 * (Ke.imag / Kc.imag - 1)
        worst = max(worst, abs(dre), abs(dim))
        print(f"{de:>8.3f}{Kc.real:>12.1f}{Ke.real:>12.1f}{dre:>8.3f}"
              f"{Kc.imag:>12.1f}{Ke.imag:>12.1f}{dim:>8.3f}")
    print("-" * 92)
    print(f"최대 차이 {worst:.3f}%  →  "
          f"{'dt 는 원인이 아니다 (28% 를 설명 못 함)' if worst < 1 else 'dt 가 기여한다'}")
    print("=" * 92)
    return worst


def stage2(S):
    """트랩까지 켠 정적 강성을 HOOMD 로. kT=0 · 큰 dt · 유령을 δ 만큼 옮겨 고정."""
    n, kt, ell, mid = S["n"], S["kt"], S["ell"], S["mid"]
    delta = S["amp"]
    lam = np.linalg.eigvalsh(S["AT"])
    dt = 0.2 / lam[-1]
    n_steps = int(round(25.0 / (lam[0] * dt)))          # 25 τ_max

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in S["idx"]:
        pos.append(list(pos[g])); typeid.append(1)
    pos[n + S["idx"].index(mid)][1] = delta              # 구동 유령만 δ 로 옮겨 고정

    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.configuration.box = [40.0 * n] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(S["idx"])]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(S["idx"])
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=2)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=S["kb"], r0=ell)
    bond.params["trap"] = dict(k=kt, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=S["kth"], t0=math.pi)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Type(["A"]), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[ov], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(n_steps)

    y = np.array(sim.state.get_snapshot().particles.position)[:n, 1]
    K_hd = kt * (delta / y[mid] - 1.0)                   # 정적: K = k_t(y_c/y_mid − 1)
    y_mod = np.linalg.solve(S["AT"], kt * delta * np.eye(n)[mid])
    K_mod = kt * (delta / y_mod[mid] - 1.0)

    print("=" * 92)
    print("② 정적 + 트랩 — 트랩을 켠 정적 강성 (①의 검증은 '강체 고정·트랩 없음' 이었다)")
    print("=" * 92)
    print(f"dt={dt:.3e}  스텝={n_steps:,}  δ={delta:.5f}")
    print(f"모델  K_static = {K_mod:10.2f}")
    print(f"HOOMD K_static = {K_hd:10.2f}   ({100*(K_hd/K_mod-1):+.2f}%)")
    print(f"y_mid: 모델 {y_mod[mid]:.6f} · HOOMD {y[mid]:.6f}  "
          f"({100*(y[mid]/y_mod[mid]-1):+.3f}%)")
    print(f"변위 프로파일 최대 차이 = {100*np.max(np.abs(y-y_mod))/delta:.3f}% of δ")
    ok = abs(K_hd / K_mod - 1) < 0.02
    print("-" * 92)
    print(f"→ {'정적은 트랩까지 켜도 일치한다. 불일치는 **동역학**이다' if ok else '정적부터 어긋난다 — 트랩/경계 쪽 문제'}")
    print("=" * 92)
    return K_hd, K_mod, y, y_mod


def stage3(S, de_target, n_cycles=3, eq_tau=12.0):
    """동적 모드 형태 — 25개 비드 전부의 복소 위상자를 재서 모델과 비교."""
    sp = min(S["specs"], key=lambda s: abs(s["params"]["De"] - de_target))
    p = sp["params"]
    om = float(p["omega_star"]); n, mid, ell = S["n"], S["mid"], S["ell"]
    dt = S["dt"]; amp = S["amp"]; kt = S["kt"]
    lam = np.linalg.eigvalsh(S["AT"])

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in S["idx"]:
        pos.append(list(pos[g])); typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.configuration.box = [40.0 * n] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(S["idx"])]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(S["idx"])
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=4)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=S["kb"], r0=ell)
    bond.params["trap"] = dict(k=kt, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=S["kth"], t0=math.pi)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Type(["A"]), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[ov], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    ghost_mid = n + S["idx"].index(mid)
    UPD = 100

    class Move(hoomd.custom.Action):
        def act(self, ts):
            y = amp * math.sin(om * ts * dt)
            with self._state.cpu_local_snapshot as s:
                tg = np.array(s.particles.tag, copy=True)
                loc = np.flatnonzero(tg == ghost_mid)
                if len(loc):
                    s.particles.position[loc[0], 1] = y

    class All(hoomd.custom.Action):
        def __init__(self): self.t, self.Y, self.g = [], [], []
        def act(self, ts):
            with self._state.cpu_local_snapshot as s:
                tg = np.array(s.particles.tag, copy=True)
                q = np.array(s.particles.position, copy=True)
                o = np.argsort(tg)
                self.t.append(ts * dt); self.Y.append(q[o][:n, 1].copy())
                self.g.append(q[o][ghost_mid, 1])

    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=Move(), trigger=hoomd.trigger.Periodic(UPD)))
    period = 2 * math.pi / om
    spc = max(20, int(round(period / dt)))
    smp = All()
    sim.operations.writers.append(hoomd.write.CustomWriter(
        action=smp, trigger=hoomd.trigger.Periodic(max(1, spc // 20))))
    n_eq = int(round(eq_tau / (lam[0] * dt)))
    sim.run(n_eq)
    smp.t.clear(); smp.Y.clear(); smp.g.clear()
    sim.run(int(round(n_cycles * period / dt)))

    t = np.array(smp.t); Y = np.array(smp.Y); g = np.array(smp.g)
    def lk(s):
        return complex(2 * np.mean(s * np.sin(om * t)), 2 * np.mean(s * np.cos(om * t)))
    y_hd = np.array([lk(Y[:, i]) for i in range(n)])
    g_hat = lk(g)
    y_mod = np.linalg.solve(1j * om * np.eye(n) + S["AT"], kt * g_hat * np.eye(n)[mid])

    print("=" * 92)
    print(f"③ 동적 모드 형태 — De={float(p['De']):.2f}  (25비드 전부의 복소 위상자)")
    print("=" * 92)
    print(f"평형화 {n_eq:,} 스텝 = {eq_tau:g} τ_max · 표본 {len(t)}")
    print(f"{'비드':>5}{'|ŷ| 모델':>12}{'|ŷ| HOOMD':>12}{'비':>8}"
          f"{'위상 모델':>11}{'위상 HOOMD':>12}{'차[rad]':>10}")
    print("-" * 92)
    for i in range(n):
        if i in S["idx"] or i % 4 == 0:
            mark = " ←트랩" if i in S["idx"] else ""
            print(f"{i:>5}{abs(y_mod[i]):>12.6f}{abs(y_hd[i]):>12.6f}"
                  f"{abs(y_hd[i])/abs(y_mod[i]):>8.3f}"
                  f"{np.angle(y_mod[i]):>11.4f}{np.angle(y_hd[i]):>12.4f}"
                  f"{np.angle(y_hd[i])-np.angle(y_mod[i]):>10.4f}{mark}")
    print("-" * 92)
    rat = np.abs(y_hd) / np.abs(y_mod)
    print(f"진폭비 범위 {rat.min():.3f} ~ {rat.max():.3f}   "
          f"{'★ 전 비드가 같은 비율 → 전체 스케일 문제' if rat.max()/rat.min() < 1.05 else '★ 비드마다 다름 → 모드 형태가 다르다'}")
    Kh = kt * g_hat / y_hd[mid] - kt - 1j * om
    Km = kt * g_hat / y_mod[mid] - kt - 1j * om
    print(f"K* 모델 ({Km.real:.0f}, {Km.imag:.0f})  ·  HOOMD ({Kh.real:.0f}, {Kh.imag:.0f})"
          f"  → K′ {100*(Kh.real/Km.real-1):+.1f}%  K″ {100*(Kh.imag/Km.imag-1):+.1f}%")
    print("=" * 92)
    np.savez(GATES / f"modeshape_de{float(p['De']):.1f}.npz",
             y_hd=y_hd, y_mod=y_mod, g_hat=g_hat, omega=om, de=float(p["De"]))
    return y_hd, y_mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="12")
    ap.add_argument("--de", type=float, default=91.8)
    ap.add_argument("--eq-tau", type=float, default=12.0)
    a = ap.parse_args()
    S = setup()
    if "1" in a.stage:
        stage1(S)
    if "2" in a.stage:
        stage2(S)
    if "3" in a.stage:
        stage3(S, a.de, eq_tau=a.eq_tau)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
