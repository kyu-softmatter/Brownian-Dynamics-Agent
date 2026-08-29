import math
import numpy as np

kB = 1.380649e-23
T = 300.0
kT = kB * T
h_planck = 6.62607015e-34

# ── ① 독립 계산: Tabor-Winterton/단순화 Lifshitz 근사 (Israelachvili 교재 표준식) ──
# A ≈ (3/4)kT[(e1-e3)/(e1+e3)]^2  +  3h*nu_e/(16*sqrt(2)) * (n1^2-n3^2)^2/(n1^2+n3^2)^1.5
eps_PMMA = 3.0      # 핸드북: PMMA 저주파(정전) 유전상수, 통상 2.6~3.5 범위
eps_water = 78.5    # 핸드북: 물@298K (system.yaml 값과 동일)
n_PMMA = 1.49        # 핸드북: PMMA 굴절률 (589nm)
n_water = 1.33       # 핸드북: 물 굴절률
nu_e = 3.0e15        # Israelachvili 표준값: 유기물·물의 주요 UV 흡수진동수 어림 (~3x10^15 Hz)

def lifshitz_A(eps1, eps3, n1, n3, nu_e=3.0e15):
    static = 0.75 * kT * ((eps1 - eps3) / (eps1 + eps3)) ** 2
    disp = (3 * h_planck * nu_e / (16 * math.sqrt(2))) * \
           (n1**2 - n3**2)**2 / (n1**2 + n3**2)**1.5
    return static, disp, static + disp

static, disp, A_calc = lifshitz_A(eps_PMMA, eps_water, n_PMMA, n_water)
print("="*72)
print("① 독립 계산 (Tabor-Winterton 근사, 핸드북 유전율·굴절률)")
print("="*72)
print(f"  정전 항 (ν=0)   = {static:.4e} J = {static/kT:.3f} kT")
print(f"  분산 항 (UV)    = {disp:.4e} J = {disp/kT:.3f} kT")
print(f"  A_calc (합)     = {A_calc:.4e} J")
print(f"  웹검색 값 A_web = 1.05e-20 J")
print(f"  차이            = {100*(A_calc-1.05e-20)/1.05e-20:+.1f}%")

# 민감도: eps_PMMA, n_PMMA 문헌 범위
print("\n  민감도 (PMMA 유전율·굴절률 문헌 범위):")
for eps1, n1, tag in [(2.6, 1.483, "하한"), (3.0, 1.49, "중심"), (3.6, 1.50, "상한")]:
    s, d, a = lifshitz_A(eps1, eps_water, n1, n_water)
    print(f"    eps={eps1}, n={n1} ({tag}): A = {a:.3e} J")

# ── ② 폴리스티렌 계열 문헌값과 비교 (같은 자릿수인지, 이미 확인된 것 재정리) ──
print()
print("="*72)
print("② 유사 폴리머 콜로이드(폴리스티렌-물) 문헌값 (웹검색 교차확인, 2026-08-05)")
print("="*72)
lit_PS = [1.2e-20, 5.0e-21, 7.0e-21, 1.5e-21*10]  # 앞서 검색된 값들 정리(마지막은 표기 오류 방지용 예시 제외)
print("  PS-물-PS: 5e-21 ~ 1.5e-20 J (여러 DLVO 논문) — PMMA 도 유사 유전특성이라 같은 자릿수 기대")

# ── ③ 민감도: A 를 바꾸면 우물이 어떻게 변하는가 (10mM 고정) ──
print()
print("="*72)
print("③ A_H 민감도 — 우물 깊이·위치가 얼마나 흔들리는가 (I=10mM MgCl2 고정)")
print("="*72)
E_CHARGE = 1.602176634e-19
NA = 6.02214076e23
EPS0 = 8.8541878128e-12
d = 1.47e-6
a = d/2
psi0 = 40e-3
eps_r = 78.5
c = 10e-3
I_SI = 0.5*(c*4+2*c*1)*NA
kappa = math.sqrt(2*I_SI*E_CHARGE**2/(EPS0*eps_r*kB*T))

def well(A_H):
    def U(h):
        return 2*math.pi*EPS0*eps_r*a*psi0**2*math.log(1+math.exp(-kappa*h)) - A_H*a/(12*h)
    hs = np.geomspace(1e-11, 1e-7, 20000)
    Us = np.array([U(h) for h in hs])
    ibar = np.argmax(Us[:len(hs)//2])
    tail = Us[ibar:]
    imin = np.argmin(tail)
    return hs[ibar+imin], tail[imin]/kT, Us[ibar]/kT

for A_H, tag in [(0.5e-20,"PS 하한"), (1.05e-20,"채택값(웹)"), (0.94e-20,"①계산값"),
                 (1.5e-20,"PS 상한"), (2.0e-20,"보수적 상한")]:
    h_min, U_min, barrier = well(A_H)
    print(f"  A={A_H:.2e} J ({tag:10s}): 2차극소 {U_min:+7.2f} kT @ h={h_min*1e9:6.2f}nm"
          f"   장벽 {barrier:6.1f} kT")
