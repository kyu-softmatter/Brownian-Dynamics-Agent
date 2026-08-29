"""`scratch/rolling_contact.py` 의 구름 저항 + 접선 스프링 모델을 정적으로 검증한다.

원칙 8 — 정적인 계를 먼저 세우고 움직임은 그 다음에 얹는다. 여기서 확정할 것:
  G0  ★ 해석적 힘·토크 == 에너지의 수치 기울기 (bd-hoomd 함정 15: 에너지만 맞고 힘이 틀리는 사례)
  G1  기준(곧은) 상태에서 U=0, F=0, τ=0
  G2  쌍 전체 강체회전에서 U=0  ← 골든 테스트
  G3  ★ 접선 스프링만이면 사슬 굽힘강성 = 0  (예측 P3, 실행 전 고정)
  G4  ★ 구름 저항만이면 κ_θ,eff = ½k_rR²  (예측 P4) — 조화굽힘 행렬과 직접 대조
  G5  방향이 얼면 곡률이 아니라 절대 회전을 벌한다 (예측 P5)
  G6  HOOMD `force.Custom(aniso=True)` 가 넘파이 참조 모델과 같은 힘·토크를 낸다

    $PY scratch/verify_rolling_contact.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from rolling_contact import (RollingContact, k_roll_from_kappa_theta,  # noqa: E402
                             kappa_theta_eff, q_from_axis_angle, q_rotate)

RNG = np.random.default_rng(20260806)
results: list[tuple[str, bool, str]] = []


def check(name, ok, detail=""):
    results.append((name, ok, detail))
    print(f"  {'✓' if ok else '✗'} {name}" + (f"   {detail}" if detail else ""))


# ── 계 설정 (chain-bend-2d-dlvo 의 실제 값) ────────────────────────────────
H_MIN_STAR = 0.00759259035993831
ELL = 1.0 + H_MIN_STAR
R_C = 0.5                                  # 접촉 반경(입자 반경, d* 단위)
KAPPA_THETA = 1391229.7767209478           # specs/...jkr...json  params.kappa_theta_star
K_ROLL = k_roll_from_kappa_theta(KAPPA_THETA, R_C)


def straight_chain(n, ell=ELL):
    pos = np.zeros((n, 3))
    pos[:, 0] = np.arange(n) * ell
    quat = np.tile(np.array([1.0, 0, 0, 0]), (n, 1))
    return pos, quat


def bonds_of(n):
    return [[i, i + 1] for i in range(n - 1)]


def zrot(theta):
    """2D — z축 회전 쿼터니언 배열."""
    return q_from_axis_angle([0, 0, 1.0], np.asarray(theta, float))


# ══════════════════════════════════════════════════════════════════════════
print("=" * 84)
print("G0 · ★ 해석적 힘·토크 vs 에너지의 수치 기울기 (함정 15 규율)")
print("=" * 84)
n = 4
pos0, quat0 = straight_chain(n)
mdl = RollingContact(bonds_of(n), pos0, quat0, R_C, k_roll=1.7, k_slide=2.3)

pos = pos0 + RNG.normal(0, 0.02, (n, 3))
pos[:, 2] = 0.0
quat = zrot(RNG.normal(0, 0.05, n))

F_an, T_an = mdl.force_torque(pos, quat)
h = 1e-7
F_num = np.zeros_like(F_an)
for i in range(n):
    for c in range(3):
        p = pos.copy(); p[i, c] += h; Up = mdl.energy(p, quat)
        p = pos.copy(); p[i, c] -= h; Um = mdl.energy(p, quat)
        F_num[i, c] = -(Up - Um) / (2 * h)
T_num = np.zeros_like(T_an)
for i in range(n):
    for c, ax in enumerate(np.eye(3)):
        q = quat.copy(); q[i] = _qm = None  # placeholder replaced below
        def rot_by(eps):
            qq = quat.copy()
            dq = q_from_axis_angle(ax, [eps])[0]
            w1, v1 = dq[0], dq[1:]
            w2, v2 = qq[i, 0], qq[i, 1:]
            qq[i] = np.concatenate([[w1 * w2 - v1 @ v2], w1 * v2 + w2 * v1 + np.cross(v1, v2)])
            return mdl.energy(pos, qq)
        T_num[i, c] = -(rot_by(h) - rot_by(-h)) / (2 * h)

fe = np.abs(F_an - F_num).max() / max(np.abs(F_num).max(), 1e-30)
te = np.abs(T_an - T_num).max() / max(np.abs(T_num).max(), 1e-30)
print(f"  최대 |F| = {np.abs(F_num).max():.6g}   최대 |τ| = {np.abs(T_num).max():.6g}")
check("힘이 −∂U/∂r 과 일치", fe < 1e-6, f"상대 최대차 {fe:.3e}")
check("토크가 −∂U/∂θ 와 일치", te < 1e-6, f"상대 최대차 {te:.3e}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("G1/G2 · 기준 상태와 강체회전 — 골든 테스트")
print("=" * 84)
pos0, quat0 = straight_chain(5)
m = RollingContact(bonds_of(5), pos0, quat0, R_C, K_ROLL, K_ROLL)
U0 = m.energy(pos0, quat0)
F0, T0 = m.force_torque(pos0, quat0)
check("G1 기준 상태 U=0, F=0, τ=0", U0 == 0.0 and np.abs(F0).max() == 0.0
      and np.abs(T0).max() == 0.0,
      f"U={U0:.3e}, max|F|={np.abs(F0).max():.3e}, max|τ|={np.abs(T0).max():.3e}")

for ang in (0.1, 0.7, 2.0):
    c, s = np.cos(ang), np.sin(ang)
    Rz = np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    pr = pos0 @ Rz.T
    qr = np.array([[np.cos(ang / 2), 0, 0, np.sin(ang / 2)]] * 5)
    Ur = m.energy(pr, qr)
    # 허용오차는 **에너지 스케일 k_rR² 대비 기계정밀도**로 잡는다. U 는 O(k_rR²) 항들의
    # 상쇄로 만들어지므로 절대 잔차가 k_rR²·1e-16 규모인 것이 정상이다.
    scale = K_ROLL * R_C ** 2
    check(f"G2 강체회전 {ang} rad 에서 U=0", abs(Ur) / scale < 1e-13,
          f"U/(k_rR²) = {Ur/scale:+.2e}   (U = {Ur:.3e}, 스케일 {scale:.4g})")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("G3/G4 · ★ 3점 굽힘 강성 — 접선 스프링만 vs 구름 저항만")
print("=" * 84)
print("  양끝 y=0, 중앙 y=δ 고정. 나머지 y 와 **모든 방향 θ** 를 이완시켜 U 를 최소화하고")
print("  k_center = 2U/δ² 를 잰다. 조화굽힘(bending_matrix) 예측과 대조한다.")


def bending_matrix(n, kappa_theta, ell):
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def harmonic_3point(n, kappa_theta, ell, delta=1.0):
    A = bending_matrix(n, kappa_theta, ell)
    c = n // 2
    fixed, free = [0, c, n - 1], [i for i in range(n) if i not in (0, c, n - 1)]
    yf = np.array([0.0, delta, 0.0])
    y = np.zeros(n)
    y[fixed] = yf
    y[free] = np.linalg.solve(A[np.ix_(free, free)], -A[np.ix_(free, fixed)] @ yf)
    return 2 * (0.5 * y @ A @ y) / delta ** 2


def rolling_3point(n, k_roll, k_slide, delta, ell=ELL, relax_orientation=True):
    pos0, quat0 = straight_chain(n, ell)
    mdl = RollingContact(bonds_of(n), pos0, quat0, R_C, k_roll, k_slide)
    c = n // 2
    free_y = [i for i in range(n) if i not in (0, c, n - 1)]

    def unpack(x):
        y = np.zeros(n)
        y[c] = delta
        y[free_y] = x[:len(free_y)]
        th = x[len(free_y):] if relax_orientation else np.zeros(n)
        p = pos0.copy()
        p[:, 1] = y
        return p, zrot(th)

    # ★ 해석적 기울기를 준다 — 수치 기울기로는 U 가 O(k_rR²) 상쇄라 최소화기가 못 내려간다
    #   (첫 시도에서 접선-only 의 U 가 **음수**로 나왔다 = 미수렴의 증거).
    #   dU/dy_i = −F_i,y,  dU/dθ_i = −τ_i,z
    def fun_jac(x):
        p, q = unpack(x)
        F, T = mdl.force_torque(p, q)
        g = -F[free_y, 1]
        if relax_orientation:
            g = np.concatenate([g, -T[:, 2]])
        return mdl.energy(p, q), g

    n_var = len(free_y) + (n if relax_orientation else 0)
    res = minimize(fun_jac, np.zeros(n_var), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=100000, maxfun=100000, ftol=0.0, gtol=0.0))
    return 2 * res.fun / delta ** 2, res


DELTA = 1e-4 * ELL
print()
print(f"  δ = {DELTA:.3e} d,  R = {R_C},  κ_θ = {KAPPA_THETA:.6g} kT  →  "
      f"k_r = 2κ_θ/R² = {K_ROLL:.6g}")
print(f"  {'n':>4} {'접선만 k':>14} {'구름만 k':>14} {'조화굽힘 k':>14} {'구름/조화':>11}")
print("  " + "-" * 62)
for nb in (3, 5, 9, 15):
    k_slide_only, _ = rolling_3point(nb, 0.0, K_ROLL, DELTA)
    k_roll_only, _ = rolling_3point(nb, K_ROLL, 0.0, DELTA)
    k_harm = harmonic_3point(nb, kappa_theta_eff(K_ROLL, R_C), ELL, DELTA)
    print(f"  {nb:4d} {k_slide_only:14.6g} {k_roll_only:14.6g} {k_harm:14.6g} "
          f"{k_roll_only/k_harm:11.6f}")
    # ★ 정규화는 k_harm 이 아니라 **접촉 강성 스케일 k·R²** 로 한다 — k_harm 은 n 에 따라
    #   1/L³ 로 줄어드는데 반올림 바닥은 안 줄어서, 비를 쓰면 n 이 커질수록 거짓 실패한다.
    scale = K_ROLL * R_C ** 2
    check(f"G3 n={nb}: 접선 스프링만 → 굽힘강성 ≈ 0",
          abs(k_slide_only) / scale < 1e-5,
          f"k_slide/(k_sR²) = {k_slide_only/scale:+.2e}  vs  "
          f"k_roll/(k_rR²) = {k_roll_only/scale:.4e}")
    check(f"G4 n={nb}: 구름만 == 조화굽힘(κ_θ=½k_rR²)",
          abs(k_roll_only / k_harm - 1.0) < 2e-3,
          f"비 = {k_roll_only/k_harm:.8f}")

print()
print("  ★ 판별 — 접선-only 의 잔차가 진짜 강성인가 반올림인가:")
print("    진짜 2차 강성이면 k = 2U/δ² 가 δ 에 **무관**하고, 상쇄 반올림이면 **k ∝ 1/δ²** 다.")
print(f"    {'δ/d':>10} {'접선만 k':>14} {'구름만 k':>14}")
print("    " + "-" * 40)
ref_s = ref_r = None
for mult in (1, 10, 100):
    ks, _ = rolling_3point(9, 0.0, K_ROLL, DELTA * mult)
    kr, _ = rolling_3point(9, K_ROLL, 0.0, DELTA * mult)
    print(f"    {DELTA*mult:10.2e} {ks:14.6g} {kr:14.6g}")
    if ref_s is None:
        ref_s, ref_r = abs(ks), kr
    else:
        check(f"G3′ δ×{mult}: 접선 잔차가 1/δ² 로 줄어듦 (= 반올림)",
              abs(ks) < ref_s / (mult ** 2) * 5,
              f"|k_slide| {ref_s:.3g} → {abs(ks):.3g} (예상 ~{ref_s/mult**2:.3g})")
        check(f"G4′ δ×{mult}: 구름 강성은 δ 에 무관 (= 진짜 2차 강성)",
              abs(kr / ref_r - 1.0) < 1e-4, f"k_roll 비 = {kr/ref_r:.8f}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("G5 · P5 방향이 얼면 다른 물리가 된다 (곡률 → 절대 회전)")
print("=" * 84)
print(f"  {'n':>4} {'θ 이완':>14} {'θ 고정':>14} {'비':>10}")
print("  " + "-" * 46)
for nb in (5, 9, 15):
    k_free, _ = rolling_3point(nb, K_ROLL, 0.0, DELTA, relax_orientation=True)
    k_frozen, _ = rolling_3point(nb, K_ROLL, 0.0, DELTA, relax_orientation=False)
    print(f"  {nb:4d} {k_free:14.6g} {k_frozen:14.6g} {k_frozen/k_free:10.4f}")
    check(f"G5 n={nb}: θ 고정이 더 뻣뻣", k_frozen > k_free * 1.5,
          f"{k_frozen/k_free:.3f}배 — 준정적↔고주파에서 갈린다")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("G6 · HOOMD force.Custom(aniso=True) ↔ 넘파이 참조 모델")
print("=" * 84)
import gsd.hoomd                                                       # noqa: E402
import hoomd                                                           # noqa: E402
import hoomd.md as md                                                  # noqa: E402
from rolling_contact import make_rolling_force                         # noqa: E402

nb = 6
pos0, quat0 = straight_chain(nb)
pos = pos0 + RNG.normal(0, 0.01, (nb, 3)); pos[:, 2] = 0.0
quat = zrot(RNG.normal(0, 0.03, nb))

f = gsd.hoomd.Frame()
f.particles.N = nb
f.particles.types = ["A"]
f.particles.typeid = [0] * nb
f.particles.position = pos
f.particles.orientation = quat
f.particles.moment_inertia = [(1.0, 1.0, 1.0)] * nb
f.configuration.box = [60, 60, 0, 0, 0, 0]
f.configuration.dimensions = 2
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_snapshot(f)
rf = make_rolling_force(bonds_of(nb), pos0, quat0, R_C, 1.7, 2.3, nb)
bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0,
                         default_gamma_r=(1.0, 1.0, 1.0))
integ = md.Integrator(dt=1e-9, methods=[bd], forces=[rf])
integ.integrate_rotational_dof = True
sim.operations.integrator = integ
sim.run(0)

mdl = RollingContact(bonds_of(nb), pos0, quat0, R_C, 1.7, 2.3)
F_ref, T_ref = mdl.force_torque(pos, quat)
F_h, T_h = np.array(rf.forces), np.array(rf.torques)
ef = np.abs(F_h - F_ref).max() / max(np.abs(F_ref).max(), 1e-30)
et = np.abs(T_h - T_ref).max() / max(np.abs(T_ref).max(), 1e-30)
check("HOOMD 힘 == 참조 모델", ef < 1e-12, f"상대 최대차 {ef:.3e}")
check("HOOMD 토크 == 참조 모델", et < 1e-12, f"상대 최대차 {et:.3e}")
check("HOOMD 총 에너지 == 참조 모델",
      abs(rf.energy - mdl.energy(pos, quat)) < 1e-9 * max(1.0, abs(rf.energy)),
      f"{rf.energy:.10g} vs {mdl.energy(pos, quat):.10g}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
n_pass = sum(1 for _, ok, _ in results if ok)
print(f"{n_pass}/{len(results)} PASS")
for name, ok, detail in results:
    if not ok:
        print(f"  ✗ {name}   {detail}")
print("=" * 84)
sys.exit(0 if n_pass == len(results) else 1)
