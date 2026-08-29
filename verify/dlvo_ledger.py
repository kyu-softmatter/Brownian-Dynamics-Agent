import math
import numpy as np

kB = 1.380649e-23
T = 300.0
kT = kB*T
e = 1.602176634e-19
NA = 6.02214076e23
eps0 = 8.8541878128e-12
eps_r = 78.5
eta = 0.851e-3
rho_p = 1180.0

d = 1.47e-6
a = d/2
psi0 = 40e-3
c_MgCl2 = 10e-3
A_H = 1.05e-20

c_Mg, c_Cl = c_MgCl2, 2*c_MgCl2
I = 0.5*(c_Mg*4 + c_Cl*1)
I_SI = I*1000*NA
kappa = math.sqrt(2*I_SI*e**2/(eps0*eps_r*kB*T))
lam_D = 1/kappa

def U_edl(h): return 2*math.pi*eps0*eps_r*a*psi0**2*math.log(1+math.exp(-kappa*h))
def U_vdw(h): return -A_H*a/(12*h)
def U_tot(h): return U_edl(h)+U_vdw(h)

hs = np.geomspace(0.05e-9, 100e-9, 20000)
Us = np.array([U_tot(h) for h in hs])
ibar = np.argmax(Us[:int(len(hs)*0.3)])   # 장벽은 앞쪽에 있음
barrier_h, barrier_U = hs[ibar], Us[ibar]
tail = Us[ibar:]
imin = np.argmin(tail)
h_min = hs[ibar+imin]
U_min = tail[imin]

# 곡률 (중심차분)
dh = h_min*1e-4
k_bond = (U_tot(h_min+dh) - 2*U_tot(h_min) + U_tot(h_min-dh))/dh**2   # d2U/dh2 [J/m^2]

print(f"λ_D = {lam_D*1e9:.4f} nm")
print(f"장벽 = {barrier_U/kT:.2f} kT @ h={barrier_h*1e9:.3f} nm")
print(f"2차극소 = {U_min/kT:.3f} kT @ h={h_min*1e9:.3f} nm")
print(f"결합 강성(곡률) k_bond = {k_bond:.4e} N/m = {k_bond*d**2/kT:.4e} kT/d^2")
print(f"결합 열요동 폭 σ_bond = sqrt(kT/k_bond) = {math.sqrt(kT/k_bond)*1e9:.3f} nm"
      f"  (우물 깊이 {abs(U_min)/kT:.2f} kT 대비 σ/h_min = {math.sqrt(kT/k_bond)/h_min:.3f})")

# ── 벌크 스케일 ──
gamma = 3*math.pi*eta*d
D_t = kT/gamma
tau_B = d**2/D_t
m = rho_p*(math.pi/6)*d**3
tau_p = m/gamma
tau_bond = gamma/k_bond
print(f"\nγ = {gamma:.4e} kg/s")
print(f"τ_B = {tau_B:.4f} s")
print(f"τ_p = {tau_p:.4e} s  (St = τ_p/τ_B = {tau_p/tau_B:.3e})")
print(f"τ_bond = γ/k_bond = {tau_bond*1e3:.4f} ms  ← 이 계의 '최속 모드' 후보 (사슬 신축)")

# 사슬 실효 길이(결합 길이 ℓ) = d + h_min (표면간극 포함)
ell = d + h_min
print(f"\n결합길이 ℓ = d + h_min = {ell*1e9:.1f} nm  (d={d*1e9:.0f}nm 의 {ell/d:.4f}배)")
for n in (5, 9, 15, 25):
    L = (n-1)*ell
    print(f"  n={n:3d}  L_chain={L*1e6:.3f} um")
