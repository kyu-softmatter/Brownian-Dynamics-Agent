"""★★ 함정 재현 — `md.angle.Harmonic(t0=π)` 은 **거의 곧은** 사슬에서 에너지는 맞고 **힘이 틀린다.**

발견 경로: `chain-bend-2d-oscill` 의 K*(ω) 가 선형응답 예측과 고주파에서 28% 어긋났다.
사다리로 좁혔더니 (`scratch/diagnose_chain_bend_28pct.py`) —
  · dt 이산화 : z-영역 정확해로 0.002% → 무관
  · 비선형·x자유도·본드신축 : scipy 로 정확 최소화 → 선형 모델과 0.32% 일치 (무관)
  · 에너지 : HOOMD 와 **정확히 0.000% 일치** (두 배치에서)
  · 힘 : ✗ **거의 곧은 배치에서 완전히 다르다** ← 원인

왜 조용히 틀리는가 — 각도 힘은 토크를 좌표로 옮길 때 `1/sin θ` 를 쓴다. t0=π 이면
평형 자체가 sin θ = 0 인 특이점이고, 뻣뻣한 사슬(κ_θ = 1.4e6 kT)은 **항상** 거의 곧다.
즉 이 함정은 강체 필라멘트에서 **가장 심하게** 작동한다. 에너지는 arccos 로 바로 계산되어
맞으므로, **에너지로 검증하면 통과하고 힘은 틀린 채로 남는다** —
`scratch/verify_angle_matrix.py` 가 2U/δ² 로 재서 0.55% 로 통과한 이유가 정확히 이것이다.

교훈: **퍼텐셜을 에너지로만 검증하면 힘을 검증한 것이 아니다.**

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_angle_force_small_theta.py
"""
from __future__ import annotations

import math

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

K_ANGLE = 1.3912297767209476e6      # chain-bend 의 κ_θ (kT 단위) — 실제 케이스 값으로 시험
N = 5                                # 손으로 확인 가능한 최소 사슬
TOL = 1e-3                           # 힘 비가 1 에서 이만큼 벗어나면 틀린 것으로 본다


def energy_exact(r, k=K_ANGLE, t0=math.pi):
    v1 = r[:-2] - r[1:-1]
    v2 = r[2:] - r[1:-1]
    c = np.sum(v1 * v2, axis=1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
    return 0.5 * k * np.sum((np.arccos(np.clip(c, -1, 1)) - t0) ** 2)


def force_exact(r, k=K_ANGLE, t0=math.pi):
    """에너지의 해석적 기울기 (리처드슨 외삽 유한차분). 힘 = −∇U."""
    G = np.zeros_like(r)
    h0 = 1e-6
    for i in range(r.shape[0]):
        for d in range(2):
            e = []
            for h in (h0, h0 / 2):
                a = r.copy(); a[i, d] += h
                b = r.copy(); b[i, d] -= h
                e.append((energy_exact(a, k, t0) - energy_exact(b, k, t0)) / (2 * h))
            G[i, d] = (4 * e[1] - e[0]) / 3
    return -G


def hoomd_angle(r, kind="harmonic", k=K_ANGLE, t0=math.pi):
    n = r.shape[0]
    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.position = np.column_stack([r, np.zeros(n)])
    f.particles.typeid = [0] * n
    f.particles.types = ["A"]
    f.configuration.box = [400.0, 400.0, 0, 0, 0, 0]     # Lz=0 → 2D (함정 9)
    f.configuration.dimensions = 2
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)

    if kind == "harmonic":
        ang = md.angle.Harmonic()
        ang.params["bend"] = dict(k=k, t0=t0)
    else:                                                 # 표 각도 — 같은 U(θ) 를 직접 준다
        nb = 1000
        th = np.linspace(0, math.pi, nb, endpoint=False)
        ang = md.angle.Table(width=nb)
        ang.params["bend"] = dict(U=0.5 * k * (th - t0) ** 2, tau=-k * (th - t0))
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0)
    integ = md.Integrator(dt=1e-14, methods=[ov], forces=[ang])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(0)
    F = np.array(ang.forces)[:, :2].copy()
    with sim.state.cpu_local_snapshot as s:
        tags = np.array(s.particles.tag, copy=True)
    out = np.zeros_like(F)
    out[tags] = F                                         # 로컬순서 → tag 순서
    return out, float(ang.energy)


def bent(scale, n=N, ell=1.0):
    """진폭 scale 의 매끄러운 굽힘 모드. scale → 0 이면 sin θ → 0."""
    y = scale * np.sin(np.pi * np.arange(n) / (n - 1))
    return np.column_stack([np.arange(n) * ell, y])


def theta_max(r):
    v1 = r[:-2] - r[1:-1]
    v2 = r[2:] - r[1:-1]
    c = np.sum(v1 * v2, axis=1) / (np.linalg.norm(v1, axis=1) * np.linalg.norm(v2, axis=1))
    return float(np.max(np.abs(np.arccos(np.clip(c, -1, 1)) - math.pi)))


def main() -> int:
    print("=" * 96)
    print("md.angle.Harmonic(t0=π) — 거의 곧은 사슬에서 에너지는 맞고 힘이 틀린다")
    print("=" * 96)
    print(f"κ_θ = {K_ANGLE:.4e} (chain-bend 실제 값) · n = {N}")
    print(f"{'max|θ−π|':>11}{'sin θ':>11}{'E 오차%':>10}{'힘비(최대성분)':>15}"
          f"{'중앙비드 |F|':>14}{'중앙 정답':>12}{'판정':>7}")
    print("-" * 96)
    first_bad = None
    for scale in (3e-1, 1e-1, 3e-2, 1e-2, 3e-3, 1.5e-3, 1e-3, 7e-4, 5e-4, 3e-4, 1e-4):
        r = bent(scale)
        tm = theta_max(r)
        Fh, Eh = hoomd_angle(r)
        Fa = Ea = None
        Fa = force_exact(r)
        Ea = energy_exact(r)
        j = int(np.argmax(np.abs(Fa[:, 1])))
        ratio = Fh[j, 1] / Fa[j, 1]
        mid = N // 2
        ok = abs(ratio - 1.0) < TOL
        if not ok and first_bad is None:
            first_bad = (tm, math.sin(math.pi - tm))
        print(f"{tm:>11.3e}{math.sin(math.pi-tm):>11.3e}{100*(Eh/Ea-1):>10.4f}"
              f"{ratio:>15.6f}{abs(Fh[mid,1]):>14.4e}{abs(Fa[mid,1]):>12.4e}"
              f"{'✓' if ok else '✗':>7}")
    print("-" * 96)
    if first_bad:
        print(f"★ 힘이 깨지기 시작하는 지점: max|θ−π| ≈ {first_bad[0]:.2e}  (sin θ ≈ {first_bad[1]:.2e})")
        print("  → 안전 규칙: **max|θ−π| ≳ 3e-3 (sin θ ≳ 3e-3) 를 유지**해야 힘이 맞는다.")
    else:
        print("힘이 전 구간에서 맞다 — 이 설치본에서는 함정이 재현되지 않는다")

    # 대안 ① angle.Table — 같은 U(θ) 를 직접 줘도 같은 경로를 타는가?
    print("\n" + "-" * 96)
    print("대안 검토")
    print("-" * 96)
    r = bent(3e-4)
    Fa = force_exact(r)
    j = int(np.argmax(np.abs(Fa[:, 1])))
    try:
        Ft, Et = hoomd_angle(r, kind="table")
        print(f"  angle.Table  (같은 U·τ 를 표로)  힘비 = {Ft[j,1]/Fa[j,1]:.6f}"
              f"   {'✓ 우회 가능' if abs(Ft[j,1]/Fa[j,1]-1)<1e-2 else '✗ 같은 문제'}")
    except Exception as e:
        print(f"  angle.Table  → {type(e).__name__}: {str(e)[:70]}")
    # 대안 ② CosineSquared 는 t0=π 에서 조화가 **아니다** (해석적으로 배제)
    print("  angle.CosineSquared: U = ½k(cosθ − cosθ₀)². θ₀=π 에서 cosθ+1 ≈ (π−θ)²/2 이므로")
    print("    U ≈ k(π−θ)⁴/8 — **4차**다. 조화 굽힘 강성을 못 낸다 → 대체 불가 (해석적으로 배제)")
    print("=" * 96)
    return 0 if first_bad else 1


if __name__ == "__main__":
    raise SystemExit(main())
