"""콜로이드 접착 접촉의 두 채널(법선·접선)을 HOOMD 파라미터로 환산.

목적: 나중에 이 계를 바로 세울 수 있도록 SI ↔ 축약 단위 표를 낸다.
출처는 전부 knowledge/source/papers/2005-pantina-furst-bending-coefficient.md
(계수는 poppler 렌더링으로 확정한 `κ₀ = 24 EI/a³`).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from simbot.units import (kT_si, stokes_drag_si, stokes_einstein_D_si,
                          water_viscosity_si)

# ── 문헌값 ────────────────────────────────────────────────────────────────
T = 298.15
A = 1.47e-6 / 2            # 입자 반지름 [m]
SIG = 2 * A                # 결합 길이 = 지름 = 축약 길이 단위
E, NU = 3.1e9, 0.4         # PMMA
A_C = 40e-9                # JKR 접촉반경 (저염 근사)
KAPPA0 = {"10 mM": 64e-3, "250 mM": 0.21, "500 mM": 0.64}   # N/m
M_C = 35e-12 * 1e-6        # 35 pN*um -> N*m
SLIP = 32e-9               # 재배열 슬립 길이 [m]
C_BEAM = 24.0              # κ₀ = 24 EI/a³  (C6 확정)

eta, _ = water_viscosity_si(T)
kT = kT_si(T)
gam = stokes_drag_si(eta, A)
D0 = stokes_einstein_D_si(T, gam)
tau_D = SIG ** 2 / D0
F_UNIT = kT / SIG                    # 힘 단위
K_UNIT = kT / SIG ** 2               # 스프링 상수 단위
E_UNIT = kT                          # 에너지·토크 단위

print("=" * 78)
print("축약 단위 (σ = 2a, kT, τ_D = σ²/D₀)")
print("=" * 78)
print(f"  σ = 2a      = {SIG*1e6:.4f} μm       a = {A*1e6:.4f} μm")
print(f"  kT          = {kT:.4e} J")
print(f"  γ           = {gam:.4e} kg/s        D₀ = {D0*1e12:.4f} μm²/s")
print(f"  τ_D         = {tau_D*1e3:.4f} ms")
print(f"  kT/σ        = {F_UNIT*1e12:.4f} pN            ← 힘 단위")
print(f"  kT/σ²       = {K_UNIT*1e6:.6f} pN/μm         ← 스프링 단위")

print()
print("=" * 78)
print("채널 1 — 법선 (JKR 접착 접촉).  HOOMD: md.bond.Harmonic")
print("=" * 78)
E_star = E / (2 * (1 - NU ** 2))          # 동일 구 2개
k_n = 2 * E_star * A_C                    # Hertz/JKR 접촉 법선 강성 dP/dδ ≈ 2E*a_c
print(f"  E* = E/2(1-ν²)   = {E_star/1e9:.4f} GPa")
print(f"  k_n = 2 E* a_c   = {k_n:.4f} N/m = {k_n/1e-6:.4g} pN/μm")
print(f"  k_n*             = {k_n/K_UNIT:.4e}   ← 축약")
print()
print("  ★ 이 강성을 그대로 넣을 수 없다. 안정성 게이트가")
dt_phys = 0.2 * 2.0 / (4 * k_n / K_UNIT)
print(f"     Δt* ≤ 0.2·2/(4 k_n*) = {dt_phys:.3e} τ_D 를 요구한다 → 완전히 비현실적.")
print("     findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md")
print("  ⇒ 법선은 '굽힘을 오염시키지 않을 만큼만' 단단하게: k_bond* = C · κ(N)*")

print()
print("=" * 78)
print("채널 2 — 접선 (접촉면 탄성).  HOOMD: md.angle.Harmonic(t0=π)")
print("=" * 78)
print(f"{'[MgCl₂]':>9} {'κ₀ [N/m]':>10} {'EI [N m²]':>12} {'k_θ [J]':>12} "
      f"{'k_θ* [kT/rad²]':>15} {'k_θ*/k_n*':>11}")
KTH = {}
for name, k0 in KAPPA0.items():
    EI = k0 * A ** 3 / C_BEAM
    kth = EI / SIG
    KTH[name] = kth
    print(f"{name:>9} {k0:>10.4g} {EI:>12.4e} {kth:>12.4e} "
          f"{kth/E_UNIT:>15.4e} {(kth/E_UNIT)/(k_n/K_UNIT):>11.2e}")
print()
print("  ⇒ 접선은 법선보다 10⁵~10⁶ 배 부드럽다. **법선을 10⁴배 물러뜨려도 굽힘은 안전하다.**")
print("     이것이 k_bond* = C·κ(N)* (C ~ 100) 가 정당한 이유다.")

print()
print("=" * 78)
print("항복 — M_c.  현재 모델에 대응물 없음 (각 스프링은 항복하지 않는다)")
print("=" * 78)
print(f"  M_c        = {M_C:.4e} N·m = {M_C/E_UNIT:.4g} kT")
print(f"  슬립 길이   = {SLIP*1e9:.0f} nm = {SLIP/SIG:.4f} σ   (≈ a_c = {A_C/SIG:.4f} σ)")
print()
print(f"{'[MgCl₂]':>9} {'항복 굽힘각 ε_c = M_c/k_θ':>26} {'[deg]':>8}")
for name, kth in KTH.items():
    eps = M_C / kth
    print(f"{name:>9} {eps:>26.4f} rad {np.degrees(eps):>8.3f}")
print()
print("  ⇒ 결합은 겨우 몇 도 굽으면 항복한다. 선형 G'(ω) 에는 무관하지만")
print("     항복·재배열을 보려면 angle.Table 로 ε > ε_c 에서 토크를 꺾어야 한다.")

print()
print("=" * 78)
print("바로 쓸 수 있는 런 설정 (10 mM, k_θ 를 문헌값으로 고정)")
print("=" * 78)
kth_star = KTH["10 mM"] / E_UNIT
print(f"  k_theta_star = {kth_star:.5g}      (= κ₀ a³ /(24 σ kT), 10 mM)")
print(f"{'N':>4} {'L*':>4} {'κ_sim* ':>11} {'κ_sim [pN/μm]':>14} {'rms y_c [σ]':>12} "
      f"{'k_bond*(C=100)':>15} {'Δt* (안정)':>12}")
for n in (7, 11, 15, 21, 31, 41):
    L = n - 1.0
    kap = 48 * kth_star / L ** 3
    kb = 100 * kap
    dt = 0.2 * 2.0 / (4 * kb + 16 * kth_star)
    print(f"{n:>4} {L:>4.0f} {kap:>11.4g} {kap*K_UNIT/1e-6:>14.4g} "
          f"{1/np.sqrt(kap):>12.5f} {kb:>15.4g} {dt:>12.3e}")
print()
print("  주의: k_θ* 를 문헌값(10 mM)으로 고정하면 짧은 사슬의 rms y_c 가 매우 작다.")
print("        지수 검증만 목적이면 k_θ* 를 낮춰 접근 가능한 영역에서 재는 것이 낫다")
print("        (지수 -3 은 k_θ 에 무관하다 — C1/C2 에서 그렇게 했다).")
