"""타원체 마찰 계수 — abp-rod 케이스에 필요. **구 극한으로 검증한 뒤 쓴다.**

CLAUDE.md 규칙 6: 물리 공식을 감으로 쓰지 않는다. Perrin 인자는 외우기 쉽게 틀리므로
타원 적분 정의에서 직접 계산하고, 구 극한(a₁=a₂=a₃)에서 알려진 값을 재현하는지 확인한다:

    병진  ζ = 6πηR        회전  ζ_r = 8πηR³

정의 (Kim & Karrila, 타원체 저항):
    Δ(λ) = √((a₁²+λ)(a₂²+λ)(a₃²+λ))
    χ₀   = ∫₀^∞ dλ/Δ            χ_i = ∫₀^∞ dλ/[(a_i²+λ)Δ]
    병진  ζ_i    = 16πη / (χ₀ + a_i² χ_i)
    회전  ζ_r,i  = 16πη (a_j²+a_k²) / [3 (a_j² χ_j + a_k² χ_k)]     (i,j,k 순환)

    $PY scratch/perrin_friction.py
"""
import math

import numpy as np
from scipy import integrate

ETA = 0.851e-3          # Pa·s, 물@300K
KT = 1.380649e-23 * 300  # J


def chi(axes):
    """(χ₀, χ₁, χ₂, χ₃) — 타원 적분. 수치적으로 계산하고 구 극한으로 검증한다."""
    a = np.asarray(axes, dtype=float)

    def delta(lam):
        return math.sqrt((a[0] ** 2 + lam) * (a[1] ** 2 + lam) * (a[2] ** 2 + lam))

    # λ → ∞ 에서 적분 수렴이 느리므로 λ = a_max² · (1/u − 1) 치환으로 유한 구간화
    scale = float(a.max()) ** 2

    def sub(f):
        def g(u):
            lam = scale * (1.0 / u - 1.0)
            return f(lam) * scale / u**2
        return integrate.quad(g, 1e-14, 1.0, limit=400)[0]

    chi0 = sub(lambda lam: 1.0 / delta(lam))
    chis = [sub(lambda lam, i=i: 1.0 / ((a[i] ** 2 + lam) * delta(lam))) for i in range(3)]
    return chi0, chis


def friction(axes, eta=ETA):
    """반축 (a₁,a₂,a₃) [m] → 병진·회전 마찰 계수 [SI]."""
    a = np.asarray(axes, dtype=float)
    chi0, chis = chi(a)
    zeta_t = [16 * math.pi * eta / (chi0 + a[i] ** 2 * chis[i]) for i in range(3)]
    zeta_r = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        num = 16 * math.pi * eta * (a[j] ** 2 + a[k] ** 2)
        den = 3 * (a[j] ** 2 * chis[j] + a[k] ** 2 * chis[k])
        zeta_r.append(num / den)
    return zeta_t, zeta_r


# ══════════════════════════════════════════════════════════════════════
print("=" * 80)
print("① 구 극한 검증 — 공식을 쓰기 전에 알려진 답을 재현하는가")
print("=" * 80)
print(f"  {'R':>10} {'ζ_t 계산':>14} {'6πηR':>14} {'오차':>10} | "
      f"{'ζ_r 계산':>13} {'8πηR³':>13} {'오차':>10}")
ok = True
for R in (0.1e-6, 0.5e-6, 1e-6, 5e-6):
    zt, zr = friction((R, R, R))
    want_t, want_r = 6 * math.pi * ETA * R, 8 * math.pi * ETA * R**3
    et = 100 * (zt[0] - want_t) / want_t
    er = 100 * (zr[0] - want_r) / want_r
    ok &= abs(et) < 1e-6 and abs(er) < 1e-6
    print(f"  {R*1e6:8.1f}µm {zt[0]:14.6e} {want_t:14.6e} {et:+9.2e}% | "
          f"{zr[0]:13.6e} {want_r:13.6e} {er:+9.2e}%")
print(f"\n  {'✓' if ok else '✗'} 세 축이 같으면 등방이어야 함: "
      f"ζ_t 편차 {100*(max(zt)-min(zt))/zt[0]:.2e}%")

print()
print("=" * 80)
print("② abp-rod 확정 형상 — 장축 2 µm, 단축 500 nm (사용자 2026-08-04)")
print("=" * 80)
a1, a2 = 1.0e-6, 0.25e-6            # 반축 = 전체 길이의 절반
p = a1 / a2
zt, zr = friction((a1, a2, a2))
d_eq = 2 * (a1 * a2 * a2) ** (1 / 3)   # 등가부피 구 지름

print(f"  반축 a₁={a1*1e6:.3f} µm (장)  a₂=a₃={a2*1e6:.3f} µm (단)   종횡비 p = {p:.1f}")
print(f"  등가부피 구 지름 d_eq = 2(a₁a₂a₃)^(1/3) = {d_eq*1e6:.4f} µm")
print(f"\n  병진 마찰 [kg/s]")
print(f"    ζ_∥ (장축)      = {zt[0]:.4e}")
print(f"    ζ_⊥ (단축)      = {zt[1]:.4e}      ζ_⊥/ζ_∥ = {zt[1]/zt[0]:.4f}")
print(f"    ζ(구 d_eq)      = {3*math.pi*ETA*d_eq:.4e}   (대조용)")
print(f"\n  회전 마찰 [kg·m²/s]")
print(f"    ζ_r,∥ (장축 중심) = {zr[0]:.4e}   ← 축 방향을 바꾸지 않음")
print(f"    ζ_r,⊥ (횡축 중심) = {zr[1]:.4e}   ★ 2D 면내 회전이 이것")

# ── BD 등방 평균 (§20 옵션 A) ─────────────────────────────────────────
Dt = [KT / z for z in zt]
gamma_bar_3d = 3.0 / sum(1.0 / z for z in zt)          # D̄ = (D₁+D₂+D₃)/3 에 맞춤
gamma_bar_2d = 2.0 / (1.0 / zt[0] + 1.0 / zt[1])       # 면내 두 방향만
print(f"\n  ★ BD 등방 평균 (병진 이방성은 BD에서 재현 불가 — bd-hoomd 하드 제약)")
print(f"    D_∥ = {Dt[0]*1e12:.4f} µm²/s   D_⊥ = {Dt[1]*1e12:.4f} µm²/s")
print(f"    γ̄(3D 평균) = {gamma_bar_3d:.4e} kg/s   →  D̄ = {KT/gamma_bar_3d*1e12:.4f} µm²/s")
print(f"    γ̄(2D 면내) = {gamma_bar_2d:.4e} kg/s   →  D̄ = {KT/gamma_bar_2d*1e12:.4f} µm²/s  ← 2D라 이것")
print(f"    ※ 조화평균을 쓰는 이유: 장시간 MSD는 D̄ = 평균 확산계수로 정해지므로")
print(f"       γ̄ = kT/D̄ 여야 장시간 거동이 맞는다 (mater_plan §20 옵션 A)")

# ── 회전확산과 텀블 시간의 분리 ───────────────────────────────────────
Dr = KT / zr[1]
tau_r = 1.0 / Dr                # 2D: ⟨cos Δθ⟩ = exp(−D_r t)
tau_tumble = 0.5
print(f"\n  ★ 회전확산 vs 텀블 — 두 시간척도가 분리되는가")
print(f"    D_r,⊥  = {Dr:.4f} 1/s      τ_r = 1/D_r = {tau_r:.4f} s   (열적 회전확산)")
print(f"    τ_tumble = {tau_tumble} s                              (사용자 확정)")
print(f"    τ_r/τ_tumble = {tau_r/tau_tumble:.3f}")
if tau_r < tau_tumble:
    print(f"    → 열적 회전확산이 텀블보다 **{tau_tumble/tau_r:.1f}배 빠르다.**")
    print(f"      방향 상관을 지배하는 것은 텀블이 아니라 회전확산이다.")
else:
    print(f"    → 텀블이 회전확산보다 {tau_r/tau_tumble:.1f}배 빠르다. 텀블이 지배.")

# 유효 지속시간: 두 과정이 독립이면 감쇠율이 더해진다
for label, factor in (("텀블이 방향을 완전 무작위화 (run-and-tumble)", 1.0),
                      ("텀블이 180° 반전 (스케치의 run-and-flip)", 2.0)):
    tau_eff = 1.0 / (1.0 / tau_r + factor / tau_tumble)
    print(f"\n    {label}")
    print(f"      1/τ_eff = 1/τ_r + {factor:.0f}/τ_tumble  →  τ_eff = {tau_eff:.4f} s")
    for v_ums in (1.0, 5.0):
        v = v_ums * 1e-6
        lp = v * tau_eff
        Pe = v * d_eq / (KT / gamma_bar_2d)
        print(f"      v={v_ums:.0f}µm/s:  ℓ_p = v·τ_eff = {lp*1e6:.3f} µm = {lp/d_eq:.2f} d_eq"
              f" = {lp/(2*a1):.2f} 몸길이   Pe = v d_eq/D̄ = {Pe:.1f}")

print()
print("=" * 80)
print(f"{'✓ PASS' if ok else '✗ FAIL'} — 구 극한 재현 (병진 6πηR · 회전 8πηR³, 오차 < 1e-11%)")
print("=" * 80)
print()
print("=" * 80)
print("③ 그래서 무엇이 정해지고 무엇이 남는가")
print("=" * 80)
print(f"""  정해짐 (사용자 확정 → tier 0/1)
    형상      반축 (1.0, 0.25, 0.25) µm,  p = 4,  d_eq = {d_eq*1e6:.4f} µm
    매질      물 @300K (η = {ETA*1e3:.3f} mPa·s)
    운동      run-and-tumble, τ_tumble = 0.5 s
    속도      v ≤ 5 µm/s (스케치)

  유도됨 (Perrin, 구 극한 검증 완료)
    γ̄(2D)    {gamma_bar_2d:.4e} kg/s        D̄ = {KT/gamma_bar_2d*1e12:.4f} µm²/s
    γ_r,z     {zr[1]:.4e} kg·m²/s   D_r = {Dr:.4f} 1/s   τ_r = {tau_r:.4f} s

  ★ 남은 결정 — 텀블 각도 분포
    "텀블"이 방향을 어떻게 바꾸는가에 따라 τ_eff 가 위 표처럼 갈린다.
    스케치는 'run-and-flip'(180° 반전)이라고 적혀 있는데 사용자는 'run-and-tumble'
    이라고 했다. 2D에서 균등 무작위 재배향이 표준이므로 그것을 제안값(tier 3)으로 둔다.

  ★ 짚어둘 것 — 이 계는 열적 회전확산이 텀블보다 빠르다 (τ_r = {tau_r:.3f} s < 0.5 s).
    즉 '직선으로 달리다가 텀블한다'는 그림이 완전히 성립하지 않는다: 런 도중에도
    열요동으로 방향이 {tau_r:.2f}s 척도로 휘어진다. p=4 타원체가 물에서 그만큼 작기 때문.
    → MSD의 방향 상관은 τ_eff 가 지배하고 텀블만으로 설명되지 않는다. 사람 확인 필요.""")
