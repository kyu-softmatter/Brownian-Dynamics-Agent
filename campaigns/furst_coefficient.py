"""Pantina & Furst (PRL 94, 138301) 의 굽힘 강성 계수 확정.

poppler 로 페이지를 렌더링해 **이미지로 직접 읽은** 값만 쓴다 (텍스트 층은 수학 글리프를 잃는다).

읽어낸 것 (전부 렌더링 이미지 확인):
  Eq. (1)   y(x) = -(F_bend/EI) ( L x^2/4  -  |x^3|/6 )
  kappa_0   = 3 pi a_c^4 E / (4 a^3)
  EI        = pi E a_c^4 / 4                     -> kappa_0 = 3 EI / a^3
  kappa     = F_bend / delta,  delta = 중앙 입자의 양 끝 대비 횡변위
  scaling   kappa(s) = kappa_0 (a/s)^(2+d_b),  d_b = 1 (직선 사슬)
  Fig 3(a)  x축 s/a, y축 kappa [pN/um], 파선 = 트랩 상한 40 pN/um
  Fig 3(a) inset  kappa(s)/kappa_0 vs s/a, 기울기 -3, 최좌측점 ~ 1.0e-3
  Fig 3(b)  kappa_0 (N/m) vs [MgCl2]: ~0.21 @250 mM, ~0.64 @500 mM
  본문      kappa_0 = 64 +- 0.5 mN/m @10 mM (실험), JKR 예측 80 mN/m
  본문      2a = 1.47 +- 0.1 um,  E = 3100 MPa, nu = 0.4, a_c ~ 40 nm
  Fig 1     11-입자 응집체, x 범위 -8..8 um
  Fig 2(a)  11-입자 @150 mM, F_bend ~ 20 pN 까지, delta ~ 1.2 um 에서 최대
"""
import numpy as np

A = 1.47e-6 / 2          # 입자 반지름 [m]
SIG = 2 * A              # 결합 길이 = 지름
E, NU = 3.1e9, 0.4
pN_um = 1e-6             # 1 pN/um = 1e-6 N/m

print("=" * 78)
print("1) Eq. (1) 의 정체 — 계수 4 와 6 이 무엇을 뜻하는가")
print("=" * 78)
print("  Eq.(1):    y(x) = -(F/EI)( L x^2/4 - |x|^3/6 )")
print("  캔틸레버:   y(x) = -(F/EI)( l x^2/2 - |x|^3/6 ),  길이 l, 자유단 하중 F")
print("  -> Eq.(1) 은 **l = L/2 인 캔틸레버**다. 3점 굽힘의 절반 (중앙은 기울기 0).")
print()
print("  곡률에서 모멘트: M(x) = EI y'' = -F(L/2 - x)")
print("     지점 x=L/2 에서 M=0  (단순지지 ✓),  중앙 x=0 에서 M = -F L/2")
print("  교과서 단순지지 중앙하중 P:  M_max = P L/4")
print("  ==> F_bend = P/2  —  F_bend 는 **지점 하나의 반력**이다")
print("      ('measured by the displacement of the end particles' 와 일치)")
print()
d_over = 1/16 - 1/48
print(f"  delta = |y(L/2)| = (F/EI) L^3 (1/16 - 1/48) = (F/EI) L^3 * {d_over:.10f}")
print(f"                    = F L^3 / {1/d_over:.0f} EI")
print(f"  kappa_PF = F_bend/delta = {1/d_over:.0f} EI / L^3")
print(f"  kappa_sim = P/delta     = {2/d_over:.0f} EI / L^3   (우리 정의: 총하중/처짐)")
print(f"  ==> kappa_PF = kappa_sim / 2.   우리 BD 는 {2/d_over:.0f} 을 0.9 % 내로 확인했다 (C2)")

print()
print("=" * 78)
print("2) Fig 3 의 s/a 는 반지름 기준 전체 윤곽길이인가 — 세 가지 독립 검사")
print("=" * 78)
print(f"  (i)  Fig 1: 11-입자, x 범위 -8..8 um.  s=(N-1)*2a = {10*SIG*1e6:.2f} um")
print(f"       -> s/a = {10*SIG/A:.0f}.  a=반지름·s=전체윤곽 이면 s/a = 2(N-1) ✓")
print()
KAP0 = {"10mM": 64e-3, "150mM": 0.15, "250mM": 0.21, "500mM": 0.64}   # Fig3(b)+본문
print("  (ii) Fig 3(a) 250 mM 최좌측 사각형: 그림에서 (s/a ~ 21, kappa ~ 24 pN/um)")
for name, expr, f in [("kappa_0 (a/s)^3", "(a/s)^3", lambda x: 1/x**3),
                      ("kappa_0 (2a/s)^3", "(2a/s)^3", lambda x: (2/x)**3)]:
    v = KAP0["250mM"] * f(21.0) / pN_um
    print(f"       {name:<18} -> {v:>9.2f} pN/um   {'✓ 일치' if 15 < v < 35 else '✗ 불일치'}")
print()
print("  (iii) Fig 2(a) 11-입자 @150 mM: F~20 pN, delta~1.2 um -> kappa ~ 16.7 pN/um")
for name, f in [("kappa_0 (a/s)^3", lambda x: 1/x**3),
                ("kappa_0 (2a/s)^3", lambda x: (2/x)**3)]:
    v = KAP0["150mM"] * f(20.0) / pN_um
    print(f"       {name:<18} -> {v:>9.2f} pN/um   {'✓ 일치' if 10 < v < 30 else '✗ 불일치'}")
print()
print("  (iv) inset 최좌측점 kappa/kappa_0 = 1.0e-3 -> (a/s)^3 = 1e-3 -> s/a = 10 ✓")
print("       (2a/s)^3 = 1e-3 이면 s/a = 20 인데, 그 점은 x축 '10' 눈금 바로 옆이다 ✗")
print()
print("  ==> 확정:  kappa_PF(s) = kappa_0 (a/s)^3,  a = 반지름,  s = (N-1)*2a")

print()
print("=" * 78)
print("3) 계수 모순 — 논문의 두 진술이 정확히 2^3 배 어긋난다")
print("=" * 78)
c_beam = 1 / d_over            # 24
print(f"  (A) Eq.(1) + Fig 3 규약에서:   kappa_0 = kappa_PF(s=a) = {c_beam:.0f} EI / a^3")
print(f"  (B) 논문이 적은 값:            kappa_0 = 3 pi a_c^4 E/(4a^3) = 3 EI / a^3")
print(f"  비 = {c_beam/3:.0f} = 2^3")
print()
print("  기원: 캔틸레버 관계 kappa = 3EI/l^3 을 쓰면서 l(=L/2, 반스팬) 을 s(전체 윤곽) 로")
print("        동일시하면 정확히 (L/2 -> L) 즉 2^3 = 8 배가 빠진다.")
print(f"        3 EI/(L/2)^3 = {3*8} EI/L^3 = {c_beam:.0f} EI/L^3  <- (A) 와 일치")

print()
print("=" * 78)
print("4) 그래서 우리 시뮬레이션의 k_theta 는 얼마인가 (논문 실측 kappa_0 기준)")
print("=" * 78)
print("  경로 A (권장) — Fig 3 의 **적합된** kappa_0 을 쓴다. 자료에 직접 붙어 있다.")
for name in ("10mM", "250mM", "500mM"):
    k0 = KAP0[name]
    EI = k0 * A ** 3 / c_beam            # kappa_0 = 24 EI/a^3
    kth = EI / SIG                       # EI = k_theta * b,  b = sigma
    print(f"    {name:>6}: kappa_0={k0*1e3:>6.1f} mN/m -> EI={EI:.4e} N m^2 -> "
          f"k_theta={kth:.4e} J")
print()
print("  경로 B (기각) — JKR a_c 로 EI 를 만들고 논문의 3EI/a^3 을 쓰는 길.")
A_C = 40e-9
EI_jkr = np.pi * E * A_C ** 4 / 4
print(f"    a_c=40 nm -> EI={EI_jkr:.4e} N m^2 -> 논문식 kappa_0=3EI/a^3="
      f"{3*EI_jkr/A**3*1e3:.1f} mN/m,  올바른 식 24EI/a^3={c_beam*EI_jkr/A**3*1e3:.1f} mN/m")
print("    본문은 JKR 예측 80 mN/m 이라고 적었다 -> 그 값은 3EI/a^3 규약의 값이다.")

print()
print("=" * 78)
print("5) 이전에 보고한 '11.77 배' 의 분해")
print("=" * 78)
EI_fit = KAP0["10mM"] * A ** 3 / c_beam
print(f"  이전 비교: 48 EI(a_c=40nm)/L^3   vs   kappa_0(a/s)^3, kappa_0=64 mN/m")
f1 = 2.0
f2 = EI_jkr / EI_fit
print(f"    인자 1 = 2.00      (48 = P/delta  vs  24 = F_bend/delta — 반력 규약)")
print(f"    인자 2 = {f2:.3f}     (EI from a_c=40nm  vs  EI from 적합 kappa_0 with 24)")
print(f"    곱     = {f1*f2:.2f}   <- 이전에 측정한 11.77 과 일치")

print()
print("=" * 78)
print("6) 최종 SI 표 — 확정된 계수로")
print("=" * 78)
print(f"{'N':>4} {'s [um]':>8} {'s/a':>6} {'kappa_PF [pN/um]':>18} {'kappa_sim [pN/um]':>19} "
      f"{'40 pN/um 상한?':>15}")
for n in (7, 9, 11, 15, 21, 31, 41):
    s = (n - 1) * SIG
    sa = s / A
    kpf = KAP0["10mM"] / sa ** 3
    print(f"{n:>4} {s*1e6:>8.2f} {sa:>6.0f} {kpf/pN_um:>18.3f} {2*kpf/pN_um:>19.3f} "
          f"{'측정불가' if kpf/pN_um > 40 else 'OK':>15}")
