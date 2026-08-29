import math

# ── 상수 ──
kB = 1.380649e-23
T = 300.0
kT = kB * T
e = 1.602176634e-19
NA = 6.02214076e23
eps0 = 8.8541878128e-12
eps_r = 78.5          # 물 (298K 근처, 핸드북)
eta = 0.851e-3

# ── 논문값 ──
d = 1.47e-6           # PRL p.1
a = d/2
psi0 = 40e-3          # PRL p.1, "measured surface potential of 40±3 mV"
c_MgCl2 = 250e-3      # Langmuir 2008 "bare" 고정 조건, mol/L

# ── 이온강도 (MgCl2 → Mg2+ + 2Cl-) ──
c_Mg = c_MgCl2
c_Cl = 2*c_MgCl2
I = 0.5*(c_Mg*(2**2) + c_Cl*(1**2))   # mol/L
I_SI = I*1000*NA   # 1/m^3 (mol/L -> mol/m^3 -> 1/m^3)

kappa = math.sqrt(2*I_SI*e**2/(eps0*eps_r*kB*T))
lambda_D = 1/kappa
print(f"이온강도 I = {I:.3f} mol/L")
print(f"디바이 길이 λ_D = {lambda_D*1e9:.4f} nm   (κ = {kappa:.3e} 1/m)")
print(f"κa = {kappa*a:.1f}   (≫1 → 평판 근사/약한 곡률 근사 타당)")

# ── Hamaker 상수 (웹 검색, 출처 불확실 — tier 3 추정) ──
A_H = 1.05e-20   # J

def U_edl(h):
    """Hogg-Healy-Fuerstenau, 등전위·약한 중첩 근사, 동일 구 반지름 a."""
    return 2*math.pi*eps0*eps_r*a*psi0**2*math.log(1+math.exp(-kappa*h))

def U_vdw(h):
    """비지연 vdW, 구-구 근접 근사 (Derjaguin)."""
    return -A_H*a/(12*h)

import numpy as np
hs_nm = np.geomspace(0.05, 50, 400)   # 표면간 거리 [nm]
hs = hs_nm*1e-9
U = np.array([U_edl(h)+U_vdw(h) for h in hs])
U_kT = U/kT

imin = np.argmin(U_kT)
print(f"\n최소점: h = {hs_nm[imin]:.3f} nm, U_min = {U_kT[imin]:.2f} kT")
print(f"\n{'h[nm]':>8} {'U_edl/kT':>10} {'U_vdw/kT':>10} {'U_tot/kT':>10}")
for h_nm in [0.1, 0.2, 0.3, 0.35, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    h = h_nm*1e-9
    ue, uv = U_edl(h)/kT, U_vdw(h)/kT
    print(f"{h_nm:8.2f} {ue:10.3f} {uv:10.3f} {ue+uv:10.3f}")

# 장벽이 있는지 확인
barrier = U_kT.max()
ibar = np.argmax(U_kT)
print(f"\n장벽(있다면) = {barrier:.3f} kT @ h={hs_nm[ibar]:.3f} nm")

print("\n" + "="*70)
print("비교: 10mM MgCl2 (PRL 스윕 조건)")
print("="*70)
c_MgCl2_2 = 10e-3
I2 = 0.5*(c_MgCl2_2*4 + 2*c_MgCl2_2*1)
I2_SI = I2*1000*NA
kappa2 = math.sqrt(2*I2_SI*e**2/(eps0*eps_r*kB*T))
lambda_D2 = 1/kappa2
print(f"λ_D(10mM) = {lambda_D2*1e9:.3f} nm   κa = {kappa2*a:.1f}")

def U_edl2(h):
    return 2*math.pi*eps0*eps_r*a*psi0**2*math.log(1+math.exp(-kappa2*h))

Us = []
for h_nm in np.geomspace(0.1, 50, 300):
    h = h_nm*1e-9
    Us.append((h_nm, (U_edl2(h)+U_vdw(h))/kT))
Us = np.array(Us)
imax = np.argmax(Us[:,1])
# 장벽 이후 최소(2차극소) 찾기
tail = Us[imax:]
imin_after = np.argmin(tail[:,1])
print(f"장벽 = {Us[imax,1]:.2f} kT @ h={Us[imax,0]:.2f} nm")
print(f"2차극소 = {tail[imin_after,1]:.2f} kT @ h={tail[imin_after,0]:.2f} nm  (장벽 이후)")
for h_nm in [0.3, 0.5, 1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50]:
    h = h_nm*1e-9
    print(f"  h={h_nm:6.2f}nm  U={( U_edl2(h)+U_vdw(h))/kT:9.3f} kT")
