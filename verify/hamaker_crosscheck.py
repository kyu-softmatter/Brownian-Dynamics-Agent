import math
import numpy as np

kB = 1.380649e-23
T = 300.0
kT = kB * T
h_planck = 6.62607015e-34

# ── (1) independent calculation: Tabor-Winterton / simplified Lifshitz
#        approximation (the standard form in Israelachvili's textbook) ──
# A ≈ (3/4)kT[(e1-e3)/(e1+e3)]^2  +  3h*nu_e/(16*sqrt(2)) * (n1^2-n3^2)^2/(n1^2+n3^2)^1.5
eps_PMMA = 3.0      # handbook: PMMA low-frequency (static) permittivity, usually 2.6-3.5
eps_water = 78.5    # handbook: water@298K (same value as in system.yaml)
n_PMMA = 1.49        # handbook: PMMA refractive index (589nm)
n_water = 1.33       # handbook: refractive index of water
nu_e = 3.0e15        # Israelachvili's standard value: rough main UV absorption frequency for
        # organics and water (~3e15 Hz)

def lifshitz_A(eps1, eps3, n1, n3, nu_e=3.0e15):
    static = 0.75 * kT * ((eps1 - eps3) / (eps1 + eps3)) ** 2
    disp = (3 * h_planck * nu_e / (16 * math.sqrt(2))) * \
           (n1**2 - n3**2)**2 / (n1**2 + n3**2)**1.5
    return static, disp, static + disp

static, disp, A_calc = lifshitz_A(eps_PMMA, eps_water, n_PMMA, n_water)
print("="*72)
print("(1) independent calculation (Tabor-Winterton approximation, handbook "
      "permittivities and refractive indices)")
print("="*72)
print(f"  static term (nu=0) = {static:.4e} J = {static/kT:.3f} kT")
print(f"  dispersion term (UV) = {disp:.4e} J = {disp/kT:.3f} kT")
print(f"  A_calc (sum)    = {A_calc:.4e} J")
print(f"  web-search value A_web = 1.05e-20 J")
print(f"  difference      = {100*(A_calc-1.05e-20)/1.05e-20:+.1f}%")

# Sensitivity: the literature ranges of eps_PMMA and n_PMMA
print("\n  sensitivity (literature range of PMMA permittivity and refractive index):")
for eps1, n1, tag in [(2.6, 1.483, "lower"), (3.0, 1.49, "central"),
                      (3.6, 1.50, "upper")]:
    s, d, a = lifshitz_A(eps1, eps_water, n1, n_water)
    print(f"    eps={eps1}, n={n1} ({tag}): A = {a:.3e} J")

# ── (2) compare with polystyrene-family literature values (is it the same order
#        of magnitude? -- restating what was already confirmed) ──
print()
print("="*72)
print("(2) literature values for a similar polymer colloid (polystyrene-water) "
      "(web cross-check, 2026-08-05)")
print("="*72)
lit_PS = [1.2e-20, 5.0e-21, 7.0e-21, 1.5e-21*10]  # the values found earlier (the last one is excluded, kept only to avoid a
  # transcription slip)
print("  PS-water-PS: 5e-21 to 1.5e-20 J (several DLVO papers) -- PMMA has similar "
      "dielectric properties, so the same order of magnitude is expected")

# ── (3) sensitivity: how does the well change with A (10mM fixed) ──
print()
print("="*72)
print("(3) A_H sensitivity -- how much the well depth and position move "
      "(I=10mM MgCl2 fixed)")
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

for A_H, tag in [(0.5e-20,"PS lower"), (1.05e-20,"adopted (web)"),
                 (0.94e-20,"(1) computed"),
                 (1.5e-20,"PS upper"), (2.0e-20,"conservative upper")]:
    h_min, U_min, barrier = well(A_H)
    print(f"  A={A_H:.2e} J ({tag:10s}): secondary min {U_min:+7.2f} kT @ "
          f"h={h_min*1e9:6.2f}nm"
          f"   barrier {barrier:6.1f} kT")
