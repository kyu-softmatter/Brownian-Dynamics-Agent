"""Convert the two channels of a colloidal adhesive contact (normal, tangential) into
HOOMD parameters.

Purpose: produce an SI <-> reduced unit table so the system can be built directly
later. Every source is
knowledge/source/papers/2005-pantina-furst-bending-coefficient.md
(the coefficient `kappa_0 = 24 EI/a^3` was settled by rendering with poppler).
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from simbot.units import (kT_si, stokes_drag_si, stokes_einstein_D_si,
                          water_viscosity_si)

# ── literature values ─────────────────────────────────────────────────────
T = 298.15
A = 1.47e-6 / 2            # particle radius [m]
SIG = 2 * A                # bond length = diameter = the reduced length unit
E, NU = 3.1e9, 0.4         # PMMA
A_C = 40e-9                # JKR contact radius (low-salt approximation)
KAPPA0 = {"10 mM": 64e-3, "250 mM": 0.21, "500 mM": 0.64}   # N/m
M_C = 35e-12 * 1e-6        # 35 pN*um -> N*m
SLIP = 32e-9               # rearrangement slip length [m]
C_BEAM = 24.0              # kappa_0 = 24 EI/a^3  (settled in C6)

eta, _ = water_viscosity_si(T)
kT = kT_si(T)
gam = stokes_drag_si(eta, A)
D0 = stokes_einstein_D_si(T, gam)
tau_D = SIG ** 2 / D0
F_UNIT = kT / SIG                    # force unit
K_UNIT = kT / SIG ** 2               # spring-constant unit
E_UNIT = kT                          # energy and torque unit

print("=" * 78)
print("reduced units (sigma = 2a, kT, tau_D = sigma^2/D_0)")
print("=" * 78)
print(f"  σ = 2a      = {SIG*1e6:.4f} μm       a = {A*1e6:.4f} μm")
print(f"  kT          = {kT:.4e} J")
print(f"  γ           = {gam:.4e} kg/s        D₀ = {D0*1e12:.4f} μm²/s")
print(f"  τ_D         = {tau_D*1e3:.4f} ms")
print(f"  kT/sigma    = {F_UNIT*1e12:.4f} pN            <- force unit")
print(f"  kT/sigma^2  = {K_UNIT*1e6:.6f} pN/um         <- spring unit")

print()
print("=" * 78)
print("channel 1 -- normal (JKR adhesive contact).  HOOMD: md.bond.Harmonic")
print("=" * 78)
E_star = E / (2 * (1 - NU ** 2))          # two identical spheres
k_n = 2 * E_star * A_C                    # Hertz/JKR normal contact stiffness
                                          # dP/ddelta ~ 2 E* a_c
print(f"  E* = E/2(1-ν²)   = {E_star/1e9:.4f} GPa")
print(f"  k_n = 2 E* a_c   = {k_n:.4f} N/m = {k_n/1e-6:.4g} pN/μm")
print(f"  k_n*             = {k_n/K_UNIT:.4e}   <- reduced")
print()
print("  ★ This stiffness cannot be used as-is. The stability gate demands")
dt_phys = 0.2 * 2.0 / (4 * k_n / K_UNIT)
print(f"     dt* <= 0.2*2/(4 k_n*) = {dt_phys:.3e} tau_D -- completely impractical.")
print("     findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md")
print("  => normal channel only as stiff as needed to not contaminate bending: k_bond* = C * kappa(N)*")

print()
print("=" * 78)
print("channel 2 -- tangential (contact-surface elasticity).  "
      "HOOMD: md.angle.Harmonic(t0=pi)")
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
print("  => tangential is 1e5-1e6x softer than normal. **Softening normal by 1e4 leaves bending safe.**")
print("     That is why k_bond* = C*kappa(N)* (C ~ 100) is justified.")

print()
print("=" * 78)
print("yield -- M_c.  No counterpart in the current model (an angle spring does not "
      "yield)")
print("=" * 78)
print(f"  M_c        = {M_C:.4e} N·m = {M_C/E_UNIT:.4g} kT")
print(f"  slip length = {SLIP*1e9:.0f} nm = {SLIP/SIG:.4f} sigma   "
      f"(~ a_c = {A_C/SIG:.4f} sigma)")
print()
print(f"{'[MgCl2]':>9} {'yield bend angle eps_c = M_c/k_theta':>26} {'[deg]':>8}")
for name, kth in KTH.items():
    eps = M_C / kth
    print(f"{name:>9} {eps:>26.4f} rad {np.degrees(eps):>8.3f}")
print()
print("  => a bond yields after only a few degrees of bending. Irrelevant to linear G'(omega), but")
print("     seeing yield/rearrangement needs angle.Table to break the torque above eps > eps_c.")

print()
print("=" * 78)
print("ready-to-use run settings (10 mM, k_theta fixed at the literature value)")
print("=" * 78)
kth_star = KTH["10 mM"] / E_UNIT
print(f"  k_theta_star = {kth_star:.5g}      (= κ₀ a³ /(24 σ kT), 10 mM)")
print(f"{'N':>4} {'L*':>4} {'κ_sim* ':>11} {'κ_sim [pN/μm]':>14} {'rms y_c [σ]':>12} "
      f"{'k_bond*(C=100)':>15} {'dt* (stable)':>12}")
for n in (7, 11, 15, 21, 31, 41):
    L = n - 1.0
    kap = 48 * kth_star / L ** 3
    kb = 100 * kap
    dt = 0.2 * 2.0 / (4 * kb + 16 * kth_star)
    print(f"{n:>4} {L:>4.0f} {kap:>11.4g} {kap*K_UNIT/1e-6:>14.4g} "
          f"{1/np.sqrt(kap):>12.5f} {kb:>15.4g} {dt:>12.3e}")
print()
print("  Caution: fixing k_theta* at the literature value (10 mM) makes rms y_c tiny")
print("        for short chains. To verify the exponent alone, lower k_theta* and")
print("        measure in an accessible range (-3 is k_theta-independent; C1/C2 did so).")
