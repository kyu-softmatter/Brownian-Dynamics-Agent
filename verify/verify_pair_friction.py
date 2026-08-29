"""`md.pair.friction` 3종이 이 프로젝트의 계에서 실제로 무엇을 하는지 실측한다.

배경 — `chain-bend-2d-dlvo` 의 JKR 대조군은 굽힘강성을 `force.Custom` 으로 직접 건다.
"마찰을 직접 구현할 수 없나"에 대해 HOOMD 7.1 에 `md.pair.friction` (Hofmann et al. 2025,
arXiv:2507.16388) 이 **있다**는 것까지는 API 서베이(docs/hoomd_capabilities.md §마찰 접촉)로
알고 있었지만 **이 프로젝트에서 한 번도 실행한 적이 없다**. 추론으로 배제하지 않고 잰다
(규칙 6).

검사할 주장 5개 —
  C1  마찰 가중 w(r) = −dU_WCA/dr 는 WCA 컷오프 밖에서 0 이다.
      이 케이스의 코어 규약(SIGMA_CORE_STAR = 2^(−1/6) → 컷오프가 정확히 r*=1)에서
      DLVO 2차극소 r* = 1+h_min = 1.00759 는 **컷오프 밖**이므로 마찰이 항등적으로 0.
  C2  (검정력 대조군) 컷오프 **안**에서는 0 이 아니다. C2 가 0 이면 C1 은 아무것도 증명 못 한다.
  C3  세 모델 전부 **산일성**이다 — 접선 변위가 있어도 상대속도가 0 이면 힘이 0.
      즉 탄성 저장(K′)이 없다. 정지마찰(stick)에 해당하는 상태가 없다.
  C4  **구름(rolling)에는 면제**된다 — 접촉점 상대 표면속도 u=0 이면 힘이 0.
      3점 굽힘의 결합 변형은 미끄러짐이 아니라 구름이므로, 이 마찰로는 굽힘을 막지 못한다.
  C5  `methods.Brownian`(과감쇠)은 속도를 물리량으로 쓰지 않는다 → BD 에서 이 힘이
      무엇을 보는가.

    $PY scratch/verify_pair_friction.py
"""
from __future__ import annotations

import math
import sys

import gsd.hoomd
import hoomd
import hoomd.md as md
import hoomd.md.pair.friction as FR
import numpy as np

# ── 이 케이스의 실제 규약 (cases/chain_bend_dlvo_2d.py:598, specs/...jkr...json) ──
SIGMA_CORE_STAR = 2 ** (-1.0 / 6.0)          # WCA 컷오프가 r*=1(표면 접촉)에서 끝나도록
R_CUT_WCA = SIGMA_CORE_STAR * 2 ** (1 / 6)   # = 1.0 (정확히)
H_MIN_STAR = 0.00759259035993831             # DLVO 2차극소 표면간극 h/d
R_WELL = 1.0 + H_MIN_STAR                    # 2차극소의 중심간 거리 r*
R_RAD = SIGMA_CORE_STAR / 2                  # 마찰 모델이 쓰는 입자 반경 (=σ/2 로 가정, T0에서 검사)

MODELS = {
    "LJLinear": (FR.FrictionLJLinear, dict(gamma_f=1.0)),
    "LJCoulomb": (FR.FrictionLJCoulomb, dict(kappa_f=3.0)),
    "LJCoulombNewton": (FR.FrictionLJCoulombNewton, dict(gamma_f=1.0, kappa_f=3.0)),
}

results: list[tuple[str, bool, str]] = []


def check(name: str, ok: bool, detail: str = "") -> None:
    results.append((name, ok, detail))
    print(f"  {'✓' if ok else '✗'} {name}" + (f"   {detail}" if detail else ""))


# ══════════════════════════════════════════════════════════════════════════
def make_frame(sep: float, v_rel: float, omega_z: float, *, mass=1.0, diameter=None):
    """중심간 거리 sep 인 두 구. x축에 놓고 y로 상대 미끄러짐 v_rel, 둘 다 ω_z 로 회전."""
    I = 0.4 * mass * R_RAD ** 2                     # 균질 구
    f = gsd.hoomd.Frame()
    f.particles.N = 2
    f.particles.types = ["A"]
    f.particles.typeid = [0, 0]
    f.particles.position = [[-sep / 2, 0, 0], [sep / 2, 0, 0]]
    f.particles.orientation = [(1, 0, 0, 0)] * 2
    f.particles.velocity = [[0, -v_rel / 2, 0], [0, v_rel / 2, 0]]
    f.particles.mass = [mass, mass]
    f.particles.moment_inertia = [(I, I, I)] * 2
    # q=(1,0,0,0) 이면 angmom = 2·q⊗(0, I·ω) = (0,0,0, 2·I·ω_z)
    f.particles.angmom = [(0, 0, 0, 2 * I * omega_z)] * 2
    if diameter is not None:
        f.particles.diameter = [diameter, diameter]
    f.configuration.box = [20, 20, 20, 0, 0, 0]
    f.configuration.dimensions = 3
    return f


def probe(model: str, sep: float, v_rel: float, omega_z: float = 0.0, *,
          kT: float = 0.0, method: str = "nve", n_steps: int = 0, diameter=None):
    """힘·토크를 잰다. 반환: (F_radial, F_tangential, |torque|) — 입자 1 기준."""
    cls, extra = MODELS[model]
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(make_frame(sep, v_rel, omega_z, diameter=diameter))
    cell = md.nlist.Cell(buffer=0.4)
    fr = cls(nlist=cell, default_r_cut=R_CUT_WCA)
    fr.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR, kT=kT, **extra)
    if method == "nve":
        meth = md.methods.ConstantVolume(filter=hoomd.filter.All())
    else:
        meth = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-6, methods=[meth], forces=[fr])
    integ.integrate_rotational_dof = True
    sim.operations.integrator = integ
    sim.run(n_steps if n_steps else 0)
    F = np.array(fr.forces)          # (2,3)
    T = np.array(fr.torques)         # (2,3)
    vel = np.array(sim.state.get_snapshot().particles.velocity)
    return dict(F_rad=float(F[1, 0]), F_tan=float(np.hypot(F[1, 1], F[1, 2])),
                torque=float(np.linalg.norm(T[1])), F=F, T=T, vel=vel)


# ══════════════════════════════════════════════════════════════════════════
print("=" * 84)
print("T0 · 규약 산술 — 이 케이스에서 DLVO 2차극소가 WCA 컷오프 밖인가")
print("=" * 84)
print(f"  SIGMA_CORE_STAR = 2^(-1/6) = {SIGMA_CORE_STAR:.10f}")
print(f"  WCA 컷오프 r_cut = σ·2^(1/6) = {R_CUT_WCA:.15f}")
print(f"  DLVO 2차극소   r* = 1 + h_min = {R_WELL:.10f}   (h_min = {H_MIN_STAR*1470:.2f} nm)")
check("r_cut == 1.0 (표면 접촉에서 정확히 끝남)", abs(R_CUT_WCA - 1.0) < 1e-12,
      f"r_cut − 1 = {R_CUT_WCA - 1.0:+.3e}")
check("C1 전제: 2차극소가 컷오프 밖", R_WELL > R_CUT_WCA,
      f"r_well/r_cut = {R_WELL / R_CUT_WCA:.6f}  (여유 {(R_WELL-R_CUT_WCA)*1470:.1f} nm)")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T1 · API — 3종 인스턴스화 + 파라미터 키")
print("=" * 84)
for name, (cls, extra) in MODELS.items():
    try:
        c = hoomd.md.nlist.Cell(buffer=0.4)
        o = cls(nlist=c, default_r_cut=R_CUT_WCA)
        o.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR, kT=0.0, **extra)
        keys = sorted(o.params[("A", "A")].keys())
        check(f"{name} 인스턴스화", True, f"params = {keys}")
    except Exception as e:                                    # noqa: BLE001
        check(f"{name} 인스턴스화", False, f"{type(e).__name__}: {e}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T2/T3 · ★ 핵심 — 컷오프 밖(2차극소) vs 안(대조군). 미끄러짐 속도 V=1.0, kT=0")
print("=" * 84)
V = 1.0
print(f"  {'model':>16} | {'r*':>8} {'F_radial':>12} {'F_tangent':>12} {'|torque|':>12}")
print("  " + "-" * 70)
for name in MODELS:
    out_well = probe(name, R_WELL, V)
    out_ctrl = probe(name, 0.95, V)
    for lbl, o in (("2차극소", out_well), ("대조군", out_ctrl)):
        r = R_WELL if lbl == "2차극소" else 0.95
        print(f"  {name if lbl=='2차극소' else '':>16} | {r:8.5f} {o['F_rad']:12.5g} "
              f"{o['F_tan']:12.5g} {o['torque']:12.5g}   {lbl}")
    check(f"C1 {name}: 2차극소에서 마찰 0",
          out_well["F_tan"] == 0.0 and out_well["torque"] == 0.0,
          f"F_tan={out_well['F_tan']:.3e}, |τ|={out_well['torque']:.3e}")
    check(f"C2 {name}: 대조군(r*=0.95)에서 마찰 ≠ 0 [검정력]",
          out_ctrl["F_tan"] > 1e-12 and out_ctrl["torque"] > 1e-12,
          f"F_tan={out_ctrl['F_tan']:.5g}, |τ|={out_ctrl['torque']:.5g}")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T4 · C3 산일성 — 접선 '변위'는 있고 '속도'가 0 이면? (탄성 저장이 있는가)")
print("=" * 84)
print("  같은 r*=0.95, 미끄러짐 속도만 V=1 → 0 으로 바꾼다. 탄성이면 힘이 남아야 한다.")
print(f"  {'model':>16} | {'V=1 F_tan':>12} {'V=0 F_tan':>12} {'V=0 |τ|':>12}")
print("  " + "-" * 60)
for name in MODELS:
    moving = probe(name, 0.95, 1.0)
    still = probe(name, 0.95, 0.0)
    print(f"  {name:>16} | {moving['F_tan']:12.5g} {still['F_tan']:12.5g} "
          f"{still['torque']:12.5g}")
    check(f"C3 {name}: 정지 상태에서 접선력 0 (탄성 저장 없음)",
          still["F_tan"] == 0.0 and still["torque"] == 0.0,
          f"F_tan={still['F_tan']:.3e} — 산일성 확인. 법선 WCA 는 {still['F_rad']:.5g} 로 살아있다")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T5a · ★ 모델이 쓰는 접촉 반경 R 은 무엇인가 — 가정하지 말고 잰다")
print("=" * 84)
print("  |τ| = R·F_tan 이므로 비를 재면 R 이 나온다. σ/2 일 것이라고 **가정했다가 T5 가 실패했다**.")
o = probe("LJLinear", 0.95, V)
R_MEAS = o["torque"] / o["F_tan"]
print(f"  |τ|/F_tan = {R_MEAS:.10f}     (σ/2 = {R_RAD:.10f} 이었다면 그 값이 나왔어야)")
check("★ R 은 σ/2 가 **아니다** — sigma 를 무시한다", abs(R_MEAS - R_RAD) > 1e-3,
      f"R_meas − σ/2 = {R_MEAS - R_RAD:+.4e}")

o_d = probe("LJLinear", 0.95, V, diameter=2 * R_RAD)
R_d = o_d["torque"] / o_d["F_tan"]
print(f"  particles.diameter = σ 로 주면 |τ|/F_tan = {R_d:.10f}")
check("★ R = particles.diameter/2 (sigma 와 독립)", abs(R_d - R_RAD) < 1e-12,
      f"diameter=σ → R = {R_d:.10f} = σ/2 ✓. 기본 diameter=1.0 이라 R 이 0.5 로 나왔던 것")
check("★★ 함정: sigma 만 설정하면 지렛대가 조용히 어긋난다",
      abs(R_MEAS / R_RAD - 1.0) > 0.1,
      f"σ={SIGMA_CORE_STAR:.4f} 인데 diameter 기본값 1.0 → R 이 {R_MEAS/R_RAD:.4f} 배 과대. "
      "에러 없음")

print()
print("=" * 84)
print("T5b · C4 구름 면제 — 잰 R 로 무슬립 조건을 다시 건다")
print("=" * 84)
omega_roll = V / (2 * R_MEAS)
print(f"  ω_z = V/(2R) = {omega_roll:.6f} (R={R_MEAS} 측정값) 이면 접촉점 u=0.")
print(f"  {'model':>16} | {'미끄러짐 F_tan':>14} {'구름 F_tan':>12} {'감소비':>11}")
print("  " + "-" * 60)
for name in MODELS:
    slide = probe(name, 0.95, V, 0.0)
    roll = probe(name, 0.95, V, omega_roll)
    ratio = roll["F_tan"] / slide["F_tan"] if slide["F_tan"] > 0 else float("nan")
    print(f"  {name:>16} | {slide['F_tan']:14.5g} {roll['F_tan']:12.5g} {ratio:11.3e}")
    check(f"C4 {name}: 무슬립 구름에서 마찰이 사라짐", ratio < 1e-9,
          f"구름/미끄러짐 = {ratio:.3e}")

print()
print("  ⚠ 잔류 슬립에 대한 민감도 — Coulomb 은 u 에 무관하게 w(r)·κ_f 를 낸다.")
print(f"  {'u/V':>8} | " + " ".join(f"{n:>16}" for n in MODELS))
print("  " + "-" * 62)
for frac in (1.0, 1e-2, 1e-4, 1e-8):
    om = (V * (1 - frac)) / (2 * R_MEAS)          # 잔류 슬립 u = frac·V
    row = [probe(n, 0.95, V, om)["F_tan"] for n in MODELS]
    print(f"  {frac:8.0e} | " + " ".join(f"{x:16.6g}" for x in row))
print("  → Linear/CoulombNewton 은 u 에 비례해 사라지지만 **Coulomb 은 u=1e−8 에서도 만힘**이다.")

# ══════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("T6 · C5 BD 결합 — methods.Brownian 의 velocity 는 마찰이 봐도 되는 양인가")
print("=" * 84)
o0 = probe("LJLinear", 0.95, 1.0, method="brownian", n_steps=0)
o1 = probe("LJLinear", 0.95, 1.0, method="brownian", n_steps=1)
o5 = probe("LJLinear", 0.95, 1.0, method="brownian", n_steps=50)
for lbl, o in (("0 스텝", o0), ("1 스텝", o1), ("50 스텝", o5)):
    print(f"  {lbl:>8}: |v| = {np.linalg.norm(o['vel'], axis=1)}   F_tan = {o['F_tan']:.5g}")
check("C5 Brownian 이 velocity 를 0 으로 두지 않는다 (마찰이 무언가를 본다)",
      float(np.abs(o5["vel"]).max()) > 0.0,
      f"50 스텝 후 max|v| = {np.abs(o5['vel']).max():.4g} — 예측(0)이 틀렸다")

print()
print("  그럼 그 velocity 는 무엇인가 — 힘 없는 자유 BD 로 Δx/dt 와 직접 대조한다.")
DT, MASS, KT = 1e-4, 1.0, 1.0
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=3)
f = make_frame(5.0, 0.0, 0.0)                       # 멀리 떨어뜨려 상호작용 없음
sim.create_state_from_snapshot(f)
bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=1.0)
integ = md.Integrator(dt=DT, methods=[bd], forces=[])
integ.integrate_rotational_dof = False
sim.operations.integrator = integ
sim.run(1)
x0 = np.array(sim.state.get_snapshot().particles.position, dtype=float)
v_before = np.array(sim.state.get_snapshot().particles.velocity, dtype=float)
sim.run(1)
x1 = np.array(sim.state.get_snapshot().particles.position, dtype=float)
dx_dt = (x1 - x0) / DT
print(f"    저장된 v      = {v_before[0]}")
print(f"    실제 Δx/dt    = {dx_dt[0]}")
print(f"    비 (v / Δx/dt) = {v_before[0] / dx_dt[0]}")
check("C5′ Brownian 의 velocity 는 실제 변위율과 무관하다",
      float(np.abs(v_before[0] / dx_dt[0]).max()) < 0.5
      or float(np.abs(v_before[0] / dx_dt[0]).min()) > 2.0,
      "→ 과감쇠 변위는 √(2Dδt) ∝ δt^½ 인데 v 는 별도로 뽑힌 값이다")

vs = []
for s in range(40):
    sm = hoomd.Simulation(device=hoomd.device.CPU(), seed=s + 1)
    sm.create_state_from_snapshot(make_frame(5.0, 0.0, 0.0))
    b = md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=1.0)
    it = md.Integrator(dt=DT, methods=[b], forces=[])
    it.integrate_rotational_dof = False
    sm.operations.integrator = it
    sm.run(2)
    vs.append(np.array(sm.state.get_snapshot().particles.velocity, dtype=float))
vs = np.concatenate(vs)
v2 = float((vs ** 2).sum(axis=1).mean())
print(f"    ⟨v²⟩ = {v2:.4f}   vs  Maxwell–Boltzmann 3kT/m = {3*KT/MASS:.4f}"
      f"   (비 {v2/(3*KT/MASS):.3f}, N={len(vs)})")
check("C5″ velocity 는 매 스텝 Maxwell–Boltzmann 에서 새로 뽑힌 값",
      abs(v2 / (3 * KT / MASS) - 1.0) < 0.25,
      f"⟨v²⟩/(3kT/m) = {v2/(3*KT/MASS):.3f} — 즉 마찰이 보는 것은 "
      "**실제 상대운동이 아니라 무상관 열잡음**이다")

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
