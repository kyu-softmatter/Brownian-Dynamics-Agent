"""Settle the bending-stiffness coefficient of Pantina & Furst (PRL 94, 138301).

Uses only values **read directly off rendered images** produced with poppler (the text
layer loses mathematical glyphs).

What was read (all confirmed against the rendered images):
  Eq. (1)   y(x) = -(F_bend/EI) ( L x^2/4  -  |x^3|/6 )
  kappa_0   = 3 pi a_c^4 E / (4 a^3)
  EI        = pi E a_c^4 / 4                     -> kappa_0 = 3 EI / a^3
  kappa     = F_bend / delta,  delta = transverse displacement of the centre particle
                              relative to the two ends
  scaling   kappa(s) = kappa_0 (a/s)^(2+d_b),  d_b = 1 (a straight chain)
  Fig 3(a)  x axis s/a, y axis kappa [pN/um], dashed line = trap limit 40 pN/um
  Fig 3(a) inset  kappa(s)/kappa_0 vs s/a, slope -3, leftmost point ~ 1.0e-3
  Fig 3(b)  kappa_0 (N/m) vs [MgCl2]: ~0.21 @250 mM, ~0.64 @500 mM
  body      kappa_0 = 64 +- 0.5 mN/m @10 mM (experiment), JKR predicts 80 mN/m
  body      2a = 1.47 +- 0.1 um,  E = 3100 MPa, nu = 0.4, a_c ~ 40 nm
  Fig 1     11-particle aggregate, x range -8..8 um
  Fig 2(a)  11 particles @150 mM, F_bend up to ~ 20 pN, peaking at delta ~ 1.2 um
"""
import numpy as np

A = 1.47e-6 / 2          # particle radius [m]
SIG = 2 * A              # bond length = diameter
E, NU = 3.1e9, 0.4
pN_um = 1e-6             # 1 pN/um = 1e-6 N/m

print("=" * 78)
print("1) what Eq. (1) actually is -- what the coefficients 4 and 6 mean")
print("=" * 78)
print("  Eq.(1):    y(x) = -(F/EI)( L x^2/4 - |x|^3/6 )")
print("  cantilever: y(x) = -(F/EI)( l x^2/2 - |x|^3/6 ),  length l, end load F")
print("  -> Eq.(1) is **a cantilever with l = L/2**: half of three-point bending "
      "(zero slope at the centre).")
print()
print("  moment from curvature: M(x) = EI y'' = -F(L/2 - x)")
print("     M=0 at the support x=L/2 (simply supported ✓), M = -F L/2 at the centre")
print("  textbook simply-supported centre load P:  M_max = P L/4")
print("  ==> F_bend = P/2  --  F_bend is **the reaction at a single support**")
print("      (consistent with 'measured by the displacement of the end particles')")
print()
d_over = 1/16 - 1/48
print(f"  delta = |y(L/2)| = (F/EI) L^3 (1/16 - 1/48) = (F/EI) L^3 * {d_over:.10f}")
print(f"                    = F L^3 / {1/d_over:.0f} EI")
print(f"  kappa_PF = F_bend/delta = {1/d_over:.0f} EI / L^3")
print(f"  kappa_sim = P/delta     = {2/d_over:.0f} EI / L^3   "
      f"(our definition: total load / deflection)")
print(f"  ==> kappa_PF = kappa_sim / 2.   our BD confirmed {2/d_over:.0f} to within "
      f"0.9 % (C2)")

print()
print("=" * 78)
print("2) is Fig 3's s/a the full contour length over the radius? -- three "
      "independent checks")
print("=" * 78)
print(f"  (i)  Fig 1: 11 particles, x range -8..8 um.  "
      f"s=(N-1)*2a = {10*SIG*1e6:.2f} um")
print(f"       -> s/a = {10*SIG/A:.0f}.  with a=radius and s=full contour, "
      f"s/a = 2(N-1) ✓")
print()
KAP0 = {"10mM": 64e-3, "150mM": 0.15, "250mM": 0.21, "500mM": 0.64}   # Fig3(b)+body
print("  (ii) Fig 3(a), the leftmost 250 mM square: from the figure "
      "(s/a ~ 21, kappa ~ 24 pN/um)")
for name, expr, f in [("kappa_0 (a/s)^3", "(a/s)^3", lambda x: 1/x**3),
                      ("kappa_0 (2a/s)^3", "(2a/s)^3", lambda x: (2/x)**3)]:
    v = KAP0["250mM"] * f(21.0) / pN_um
    print(f"       {name:<18} -> {v:>9.2f} pN/um   "
          f"{'✓ agrees' if 15 < v < 35 else '✗ disagrees'}")
print()
print("  (iii) Fig 2(a), 11 particles @150 mM: F~20 pN, delta~1.2 um "
      "-> kappa ~ 16.7 pN/um")
for name, f in [("kappa_0 (a/s)^3", lambda x: 1/x**3),
                ("kappa_0 (2a/s)^3", lambda x: (2/x)**3)]:
    v = KAP0["150mM"] * f(20.0) / pN_um
    print(f"       {name:<18} -> {v:>9.2f} pN/um   "
          f"{'✓ agrees' if 10 < v < 30 else '✗ disagrees'}")
print()
print("  (iv) the inset's leftmost point kappa/kappa_0 = 1.0e-3 -> (a/s)^3 = 1e-3 "
      "-> s/a = 10 ✓")
print("       (2a/s)^3 = 1e-3 would give s/a = 20, but that point sits right next to "
      "the '10' tick ✗")
print()
print("  ==> settled:  kappa_PF(s) = kappa_0 (a/s)^3,  a = radius,  s = (N-1)*2a")

print()
print("=" * 78)
print("3) the coefficient contradiction -- two statements in the paper differ by "
      "exactly 2^3")
print("=" * 78)
c_beam = 1 / d_over            # 24
print(f"  (A) from the Eq.(1) + Fig 3 convention:  "
      f"kappa_0 = kappa_PF(s=a) = {c_beam:.0f} EI / a^3")
print(f"  (B) the value the paper writes:  "
      f"kappa_0 = 3 pi a_c^4 E/(4a^3) = 3 EI / a^3")
print(f"  ratio = {c_beam/3:.0f} = 2^3")
print()
print("  origin: using the cantilever relation kappa = 3EI/l^3 while identifying")
print("        l (= L/2, half-span) with s (full contour) loses exactly 2^3 = 8.")
print(f"        3 EI/(L/2)^3 = {3*8} EI/L^3 = {c_beam:.0f} EI/L^3  <- matches (A)")

print()
print("=" * 78)
print("4) so what is k_theta for our simulation (against the paper's measured "
      "kappa_0)")
print("=" * 78)
print("  route A (recommended) -- use Fig 3's **fitted** kappa_0, attached "
      "directly to the data.")
for name in ("10mM", "250mM", "500mM"):
    k0 = KAP0[name]
    EI = k0 * A ** 3 / c_beam            # kappa_0 = 24 EI/a^3
    kth = EI / SIG                       # EI = k_theta * b,  b = sigma
    print(f"    {name:>6}: kappa_0={k0*1e3:>6.1f} mN/m -> EI={EI:.4e} N m^2 -> "
          f"k_theta={kth:.4e} J")
print()
print("  route B (rejected) -- build EI from the JKR a_c and use the paper's 3EI/a^3.")
A_C = 40e-9
EI_jkr = np.pi * E * A_C ** 4 / 4
print(f"    a_c=40 nm -> EI={EI_jkr:.4e} N m^2 -> paper's form kappa_0=3EI/a^3="
      f"{3*EI_jkr/A**3*1e3:.1f} mN/m,  correct form "
      f"24EI/a^3={c_beam*EI_jkr/A**3*1e3:.1f} mN/m")
print("    the body text says JKR predicts 80 mN/m -> that value is in the 3EI/a^3 "
      "convention.")

print()
print("=" * 78)
print("5) decomposition of the previously reported 'factor 11.77'")
print("=" * 78)
EI_fit = KAP0["10mM"] * A ** 3 / c_beam
print(f"  the earlier comparison: 48 EI(a_c=40nm)/L^3   vs   kappa_0(a/s)^3, "
      f"kappa_0=64 mN/m")
f1 = 2.0
f2 = EI_jkr / EI_fit
print(f"    factor 1 = 2.00    (48 = P/delta  vs  24 = F_bend/delta -- the reaction "
      f"convention)")
print(f"    factor 2 = {f2:.3f}   (EI from a_c=40nm  vs  EI from the fitted kappa_0 "
      f"with 24)")
print(f"    product  = {f1*f2:.2f}   <- matches the 11.77 measured earlier")

print()
print("=" * 78)
print("6) final SI table -- with the settled coefficient")
print("=" * 78)
print(f"{'N':>4} {'s [um]':>8} {'s/a':>6} {'kappa_PF [pN/um]':>18} {'kappa_sim [pN/um]':>19} "
      f"{'above 40 pN/um?':>15}")
for n in (7, 9, 11, 15, 21, 31, 41):
    s = (n - 1) * SIG
    sa = s / A
    kpf = KAP0["10mM"] / sa ** 3
    print(f"{n:>4} {s*1e6:>8.2f} {sa:>6.0f} {kpf/pN_um:>18.3f} {2*kpf/pN_um:>19.3f} "
          f"{'unmeasurable' if kpf/pN_um > 40 else 'OK':>15}")
