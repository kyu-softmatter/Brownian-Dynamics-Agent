import math

# ── constants ──
kB = 1.380649e-23
T = 300.0
kT = kB * T
e = 1.602176634e-19
NA = 6.02214076e23
eps0 = 8.8541878128e-12
eps_r = 78.5          # water (near 298K, handbook)
eta = 0.851e-3

# ── values from the papers ──
d = 1.47e-6           # PRL p.1
a = d/2
psi0 = 40e-3          # PRL p.1, "measured surface potential of 40±3 mV"
c_MgCl2 = 250e-3      # Langmuir 2008 "bare" fixed condition, mol/L

# ── ionic strength (MgCl2 -> Mg2+ + 2Cl-) ──
c_Mg = c_MgCl2
c_Cl = 2*c_MgCl2
I = 0.5*(c_Mg*(2**2) + c_Cl*(1**2))   # mol/L
I_SI = I*1000*NA   # 1/m^3 (mol/L -> mol/m^3 -> 1/m^3)

kappa = math.sqrt(2*I_SI*e**2/(eps0*eps_r*kB*T))
lambda_D = 1/kappa
print(f"ionic strength I = {I:.3f} mol/L")
print(f"Debye length lambda_D = {lambda_D*1e9:.4f} nm   (kappa = {kappa:.3e} 1/m)")
print(f"kappa*a = {kappa*a:.1f}   (>>1 -> the flat-plate / weak-curvature "
      f"approximation is valid)")

# ── Hamaker constant (from a web search, provenance uncertain -- tier 3 estimate) ──
A_H = 1.05e-20   # J

def U_edl(h):
    """Hogg-Healy-Fuerstenau, constant-potential weak-overlap approximation, equal
    spheres of radius a."""
    return 2*math.pi*eps0*eps_r*a*psi0**2*math.log(1+math.exp(-kappa*h))

def U_vdw(h):
    """Non-retarded vdW, sphere-sphere proximity approximation (Derjaguin)."""
    return -A_H*a/(12*h)

import numpy as np
hs_nm = np.geomspace(0.05, 50, 400)   # surface-to-surface separation [nm]
hs = hs_nm*1e-9
U = np.array([U_edl(h)+U_vdw(h) for h in hs])
U_kT = U/kT

imin = np.argmin(U_kT)
print(f"\nminimum: h = {hs_nm[imin]:.3f} nm, U_min = {U_kT[imin]:.2f} kT")
print(f"\n{'h[nm]':>8} {'U_edl/kT':>10} {'U_vdw/kT':>10} {'U_tot/kT':>10}")
for h_nm in [0.1, 0.2, 0.3, 0.35, 0.5, 0.7, 1.0, 1.5, 2.0, 3.0, 5.0, 10.0]:
    h = h_nm*1e-9
    ue, uv = U_edl(h)/kT, U_vdw(h)/kT
    print(f"{h_nm:8.2f} {ue:10.3f} {uv:10.3f} {ue+uv:10.3f}")

# check whether a barrier exists
barrier = U_kT.max()
ibar = np.argmax(U_kT)
print(f"\nbarrier (if any) = {barrier:.3f} kT @ h={hs_nm[ibar]:.3f} nm")

print("\n" + "="*70)
print("comparison: 10mM MgCl2 (the PRL sweep condition)")
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
# find the minimum after the barrier (the secondary minimum)
tail = Us[imax:]
imin_after = np.argmin(tail[:,1])
print(f"barrier = {Us[imax,1]:.2f} kT @ h={Us[imax,0]:.2f} nm")
print(f"secondary minimum = {tail[imin_after,1]:.2f} kT @ "
      f"h={tail[imin_after,0]:.2f} nm  (after the barrier)")
for h_nm in [0.3, 0.5, 1, 2, 3, 5, 8, 10, 15, 20, 30, 40, 50]:
    h = h_nm*1e-9
    print(f"  h={h_nm:6.2f}nm  U={( U_edl2(h)+U_vdw(h))/kT:9.3f} kT")
