"""굽힘 강성 행렬을 **HOOMD 의 `angle.Harmonic` 자체와** 대조한다.

★ 왜 필요한가 — 지금까지 검증된 것은 `이산 행렬 vs 빔 공식 48EI/L³` (−0.35%) 뿐이고,
그건 **해석 대 해석**이다. `cases/chain_bend_2d.bending_matrix` 가 HOOMD 의
`md.angle.Harmonic(k=κ_θ, t0=π)` 와 같은 물리를 내는지는 **한 번도 확인되지 않았다.**
규칙 4·6 이 정확히 이 상황을 가리킨다 (감으로 쓰지 않는다 · 실행으로 확인한다).

구성: 3점 굽힘의 정적 문제. 비드 0·12·24 를 y = 0·δ·0 에 고정(적분에서 제외)하고
자유 비드를 완화시킨 뒤, 비드 12 에 걸린 굽힘+본드 힘을 읽어 k_center = −F_y/δ 를 얻는다.
트랩을 쓰지 않으므로 트랩 강성이 섞이지 않는다 — 순수하게 사슬의 강성만 본다.

대조 대상 3종:
  ① 이산 행렬 (bending_matrix 로 푼 정적해)
  ② 연속체 빔 48EI/L³
  ③ HOOMD 실측                                  ← 이것이 새로 재는 것

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_angle_matrix.py
"""
from __future__ import annotations

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
sys.path.insert(0, str(ROOT))
from cases.chain_bend_2d import bending_matrix, sweep_specs, trapped_indices  # noqa: E402

TOL = 1.0          # 허용 오차 [%] — 이산 행렬과 HOOMD 가 같은 물리여야 한다


def matrix_static(n, kth, ell, delta, fixed):
    """이산 행렬로 푼 정적 3점 굽힘. k_center = 2U/δ² (three_point_bending 과 같은 정의)."""
    A = bending_matrix(n, kth, ell)
    free = [i for i in range(n) if i not in fixed]
    y = np.zeros(n)
    y[fixed] = [0.0, delta, 0.0]
    y[free] = np.linalg.solve(A[np.ix_(free, free)], -A[np.ix_(free, fixed)] @ y[fixed])
    U = 0.5 * y @ A @ y
    return 2.0 * U / delta ** 2, y


def hoomd_static(n, kth, k_bond, ell, delta, fixed, dt, n_steps):
    """HOOMD 로 같은 문제를 푼다. 고정 비드는 타입 'F' 로 두고 적분에서 뺀다."""
    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    pos[fixed[1]][1] = delta                      # 중앙을 δ 만큼 밀어 둔다
    typeid = [1 if i in fixed else 0 for i in range(n)]

    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "F"]
    f.configuration.box = [40.0 * ell, 40.0 * ell, 0, 0, 0, 0]   # Lz=0 → 2D (함정 9)
    f.configuration.dimensions = 2
    f.bonds.N = n - 1
    f.bonds.types = ["backbone"]
    f.bonds.typeid = [0] * (n - 1)
    f.bonds.group = np.array([[i, i + 1] for i in range(n - 1)])
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=k_bond, r0=ell)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=kth, t0=math.pi)
    # 열잡음 없는 과감쇠 — 정적 해로 완화시킨다 (skill bd-hoomd)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Type(["A"]), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[ov], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(n_steps)

    # ★ 에너지는 스칼라라 입자 순서 문제가 없다 (per-particle force 배열은 ParticleSorter
    #   때문에 tag 순서가 아닐 수 있다 — skill bd-hoomd). k_center = 2U/δ² 로 잰다.
    U = float(angle.energy) + float(bond.energy)
    snap = sim.state.get_snapshot()               # 전역 스냅샷은 tag 순서다
    y = np.array(snap.particles.position)[:, 1]
    return 2.0 * U / delta ** 2, y, (float(angle.energy), float(bond.energy))


def main() -> int:
    # ★ 스펙 파일에 의존하지 않는다 — 하드 검사 실패 시 스펙이 없는 것이 정상이다
    p = sweep_specs()[0]["params"]
    n = int(p["n_beads"])
    kth = float(p["kappa_theta_star"])
    k_bond = float(p["k_bond_star"])
    ell = float(p["L_chain_star"]) / (n - 1)
    L = float(p["L_chain_star"])
    delta = float(p["amp_star"])
    fixed = trapped_indices(n)
    EI = kth * ell
    beam = 48 * EI / L ** 3

    lam_max = float(np.linalg.eigvalsh(bending_matrix(n, kth, ell))[-1])
    dt = 0.2 / lam_max                    # 명시적 오일러 안정 한계(2/λ)의 10%
    n_steps = int(round(20.0 / (0.2 * 526.4 / lam_max)))   # ≈ 20 τ_max

    print("=" * 88)
    print("굽힘 강성: 이산 행렬 vs 연속체 빔 vs **HOOMD angle.Harmonic 실측**")
    print("=" * 88)
    print(f"n={n}  κ_θ*={kth:.4e}  ℓ={ell:.3f}  L={L:.1f}  δ={delta:.5f}  고정={fixed}")
    print(f"λ_max={lam_max:.4e}  dt={dt:.3e}  스텝={n_steps:,}\n")

    k_mat, y_mat = matrix_static(n, kth, ell, delta, fixed)
    k_hd, y_hd, (U_ang, U_bond) = hoomd_static(n, kth, k_bond, ell, delta, fixed,
                                               dt, n_steps)

    print(f"① 이산 행렬          k_center = {k_mat:10.2f}")
    print(f"② 연속체 빔 48EI/L³  k_center = {beam:10.2f}   (① 대비 {100*(k_mat/beam-1):+.2f}%)")
    print(f"③ HOOMD 실측         k_center = {k_hd:10.2f}   (① 대비 "
          f"{100*(k_hd/k_mat-1):+.2f}%   ② 대비 {100*(k_hd/beam-1):+.2f}%)")
    print(f"\nHOOMD 에너지 분해: 굽힘 {U_ang:.4f} + 본드 {U_bond:.4e} kT"
          f"  (본드 비중 {100*U_bond/(U_ang+U_bond):.3f}%)")
    print(f"변위 프로파일 최대 차이 (행렬 vs HOOMD) = "
          f"{100*np.max(np.abs(y_hd-y_mat))/delta:.3f}% of δ")

    err = 100 * abs(k_hd / k_mat - 1)
    ok = err < TOL
    print("-" * 88)
    if ok:
        print(f"✓ PASS — 행렬과 HOOMD 가 {err:.2f}% 안에서 같다. bending_matrix 는 "
              f"angle.Harmonic 을 옳게 표현한다")
    else:
        print(f"✗ FAIL — 행렬과 HOOMD 가 {err:.2f}% 다르다. 둘은 **같은 물리가 아니다** — "
              f"행렬로 유도한 λ_max·τ_fast·dt·K*(ω) 예측이 전부 영향을 받는다")
    print("=" * 88)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
