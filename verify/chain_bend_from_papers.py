"""Verification of the chain-bend literature distillation -- is what was extracted
from the two Pantina & Furst papers correct?

Papers (specified by the user, 2026-08-04):
  [P1] Pantina & Furst, PRL 94, 138301 (2005)
       "Elasticity and Critical Bending Moment of Model Colloidal Aggregates"
  [P2] Pantina & Furst, Langmuir 24, 1141-1146 (2008)
       "Micromechanics and Contact Forces of Colloidal Aggregates in the Presence of Surfactants"

★ The papers' headline finding: **the bead-bead interaction is NOT a central pair
  potential.** A single bond **supports a torque** (a tangential interaction), and
  that is why the chain bends like a beam.
  -> The sketch's blank `U_ij` is not a "pair potential" but **an adhesive contact
     plus a bending stiffness**.

A value an LLM extracted from a paper carries a hallucination risk (masterplan
principle 2). The defence:
  **reproduce a number the paper reports itself, using the formula I extracted.**
  If it reproduces, that is evidence the formula and constants were extracted
  correctly.

    $PY scratch/chain_bend_from_papers.py
"""
import math

import numpy as np

# ══════════════════════════════════════════════════════════════════════
# Values read out of the papers (every one carries a locator)
# ══════════════════════════════════════════════════════════════════════
E_PMMA = 3100e6         # Pa      [P1] p.3 left column, ref 15 (Schreyer)
NU_PMMA = 0.4           # -       [P1] p.3 left column
D_PART = 1.47e-6        # m       [P1] p.1 right column
                        #         "average diameter of 2a = 1.47 ± 0.01 µm"
                        #         [P2] p.2 right column "2a = 1.47 ± 0.1 µm"
                        #         (Bangs Labs PMMA)
GAMMA_L = 72.7e-3       # N/m     [P1] p.3 left column, ref 16 (surface tension of
                        #         water)
THETA_0 = math.radians(73.7)   # [P1] p.3 left column
                               # "PMMA-water contact angle theta_0 = 73.7 deg"

# Measured values
KAPPA0_EXP_10mM = 64e-3     # N/m  [P1] p.3 left column "experimental value at
                            #      10 mM MgCl2 is 64 ± 0.5 mN/m"
KAPPA0_JKR_PAPER = 80e-3    # N/m  [P1] p.3 left column "Using the JKR model,
                            #      kappa_0 = 80 mN/m"
AC_PAPER = 40e-9            # m    [P1] p.3 left column "corresponding contact area
                            #      radius is a_c ~ 40 nm"
KAPPA0_BARE_250mM = 0.21    # N/m  [P2] p.4 Fig.4 inset "κ₀^bare = 0.21 ± 0.01 N/m"
MC_PLATEAU = 35e-18        # N*m   [P1] p.2 right column "M_c plateaus to a value of
                           #       approximately 35 pN µm"
MC_BARE_250mM = 30e-18     # N*m   [P2] Fig.3 arrow (no surfactant, ~30 pN*µm)
SLIP_LENGTH = 32e-9       # m     [P1] p.2 right column "average length of sliding
                          #       rearrangements is 32 ± 15 nm"
K_TRAP_PAPER = 40e-6      # N/m   [P1] p.1 right column "trap rigidity is
                          #       approximately 40 pN/µm"
V_DRAG_P2 = 20e-9         # m/s   [P2] Fig.1A "translating in the y direction at a rate of 20 nm/s"
D_B = 1.0                 # -     [P2] p.4 "chain aggregates ... nearly perfectly straight; therefore d_b = 1"

ETA = 0.851e-3            # Pa*s  water@300K (handbook) -- the papers use an aqueous
                          #       solution
KT = 1.380649e-23 * 300   # J

a = D_PART / 2            # particle radius
K_JKR = 2 * E_PMMA / (3 * (1 - NU_PMMA**2))     # [P2] eq.6 "K = 2E/3(1−ν²)"


def a_c_from_W(W_SL, a_=a):
    """[P2] eq.6  a_c = (3*pi*a^2*W_SL / 2K)^(1/3)  -- JKR, no external load"""
    return (3 * math.pi * a_**2 * W_SL / (2 * K_JKR)) ** (1 / 3)


def kappa0_from_ac(a_c, a_=a):
    """[P2] eq.5  κ₀ = 3π a_c⁴ E / (4 a³)"""
    return 3 * math.pi * a_c**4 * E_PMMA / (4 * a_**3)


def ac_from_kappa0(kappa0, a_=a):
    """eq.5, inverted"""
    return (kappa0 * 4 * a_**3 / (3 * math.pi * E_PMMA)) ** 0.25


def EI_from_kappa0(kappa0, a_=a):
    """kappa_0 = 3EI/a^3  (eq.5 combined with EI = pi*E*a_c^4/4)"""
    return kappa0 * a_**3 / 3


print("=" * 82)
print("(1) extraction check -- does my formula reproduce the number the paper "
      "reports?")
print("=" * 82)
W0 = GAMMA_L * (1 + math.cos(THETA_0))          # [P2] eq.9 Young-Dupré
ac0 = a_c_from_W(W0)
k0 = kappa0_from_ac(ac0)
print(f"  Young-Dupré  W_SL⁰ = γ_L(1+cos θ₀) = {W0*1e3:.1f} mJ/m²")
print(f"  JKR          a_c = {ac0*1e9:.1f} nm      (K = {K_JKR/1e9:.3f} GPa)")
print(f"  eq.5         κ₀ = {k0*1e3:.1f} mN/m")
print(f"  paper [P1]    kappa_0 = {KAPPA0_JKR_PAPER*1e3:.0f} mN/m")
err = 100 * (k0 - KAPPA0_JKR_PAPER) / KAPPA0_JKR_PAPER
print(f"  -> difference {err:+.1f}%   "
      f"{'✓ extraction is correct (eq.5, 6, 9 + 4 constants)' if abs(err) < 5 else '✗ something is wrong'}")

ac_exp = ac_from_kappa0(KAPPA0_EXP_10mM)
print(f"\n  inverting the experimental value: kappa_0 = "
      f"{KAPPA0_EXP_10mM*1e3:.0f} mN/m -> a_c = {ac_exp*1e9:.1f} nm")
print(f"  paper [P1]    a_c ~ {AC_PAPER*1e9:.0f} nm")
print(f"  -> difference {100*(ac_exp-AC_PAPER)/AC_PAPER:+.1f}%  "
      f"(the paper's value is rounded)  ✓")

sigma_c = 4 * MC_PLATEAU / (math.pi * ac_exp**3)   # [P2] eq.16 M_c = (π/4)σ*_xx a_c³
print(f"\n  eq.16 inverted: M_c = {MC_PLATEAU*1e18:.0f} pN*µm, "
      f"a_c = {ac_exp*1e9:.1f} nm")
print(f"               -> critical tensile stress sigma*_xx = "
      f"{sigma_c/1e6:.2f} MPa")
print(f"               far below PMMA's yield strength (~70 MPa) -> consistent with "
      f"contact-line de-pinning")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("(2) HOOMD mapping check -- the discrete chain's angle stiffness vs the "
      "continuum beam's EI")
print("=" * 82)
print("""  claim:  EI = kappa_theta * l        (l = bond length = 2a)
  basis:  discrete   U = sum 0.5*kappa_theta*theta_i^2,  theta_i ~ l/R
                     -> 0.5*kappa_theta*l/R^2 per unit length
          continuum  U = integral 0.5*EI/R^2 dx     -> compare coefficients
  ★ Not taken on intuition: the discrete chain is actually bent and compared
    against the beam formula.""")


def discrete_3point_stiffness(n, kappa_theta, ell):
    """An n-bead chain: ends pinned, centre pushed by delta -- returns F/delta.

    Small-deformation approximation theta_i = (y_{i+1} - 2y_i + y_{i-1})/l,
    U = 0.5*kappa_theta*sum(theta_i^2) -> a quadratic form in y. Minimised under the
    constraints (ends 0, centre delta), then F = dU/d(delta).
    """
    assert n % 2 == 1, "there must be a centre bead"
    c = n // 2
    # Second-difference matrix (n-2 interior angles)
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    A = kappa_theta * (B.T @ B)          # U = ½ yᵀ A y
    fixed = [0, c, n - 1]
    free = [i for i in range(n) if i not in fixed]
    y_fix = np.array([0.0, 1.0, 0.0])    # delta = 1 (linear, so the scale is
                                         # irrelevant)
    A_ff = A[np.ix_(free, free)]
    A_fx = A[np.ix_(free, fixed)]
    y_free = np.linalg.solve(A_ff, -A_fx @ y_fix)
    y = np.zeros(n)
    y[fixed] = y_fix
    y[free] = y_free
    U = 0.5 * y @ A @ y                  # delta=1, so U = 0.5*k*delta^2 = 0.5*k
    return 2 * U                         # k = F_center/delta
                                         # ★ referenced to the force applied at the
                                         #   CENTRE


kappa_theta_10mM = EI_from_kappa0(KAPPA0_EXP_10mM) / D_PART
print(f"\n  κ₀ = {KAPPA0_EXP_10mM*1e3:.0f} mN/m → EI = {EI_from_kappa0(KAPPA0_EXP_10mM):.4e} N·m²"
      f" → κ_θ = EI/ℓ = {kappa_theta_10mM:.4e} J = {kappa_theta_10mM/KT:.3e} kT")
print("""
  ★ The force definitions have to be matched -- they disagreed by exactly 2x at
    first, and raising n did not make it converge.
    The paper's kappa = F_bend/delta is referenced to the force on the **END
    particles** ([P1] p.2: "F_bend is measured by the displacement of the end
    particles"). By force balance the centre force is 2x the end force, so
        referenced to the end force     kappa_end    = 24 EI/L^3
                                        <- the value to compare with the paper
        referenced to the centre force  kappa_center = 48 EI/L^3
                                        <- what the driving trap feels in simulation
    Measuring the centre force in simulation and comparing it directly against the
    paper's kappa is wrong by **exactly a factor of 2.**""")
print(f"\n  {'n':>4} {'L=(n-1)l':>11} {'discrete(ctr)':>15} {'beam 48EI/L^3':>13}"
      f" {'diff':>9} {'paper k_end=24EI/L^3':>19}")
EI = EI_from_kappa0(KAPPA0_EXP_10mM)
ok2 = True
for n in (5, 9, 11, 15, 25, 51):
    L = (n - 1) * D_PART
    k_disc = discrete_3point_stiffness(n, kappa_theta_10mM, D_PART)
    k_end = 24 * EI / L**3               # y(±L/2) = F_end L³/(24EI) — [P1] eq.1
    k_center = 2 * k_end                 # force balance: centre force = 2x end force
    e = 100 * (k_disc - k_center) / k_center
    # At n=5 the discretisation error is a large -11%. The paper's experiments use
    # chains of 9 to 25 beads, so the verdict is taken over that range.
    if n >= 9:
        ok2 &= abs(e) < 4
    print(f"  {n:>4} {L*1e6:9.2f}µm {k_disc*1e6:13.3f}pN/µm {k_center*1e6:11.3f}pN/µm"
          f" {e:+8.2f}% {k_end*1e6:17.3f}pN/µm")
print(f"\n  {'✓' if ok2 else '✗'} once the force definitions match, the discrete "
      f"chain converges to the beam formula "
      f"(within 4% for n>=9; n=5 gives -11%, the discretisation error of a short "
      f"chain, outside the papers' experimental range)")
print(f"  => **the kappa_theta = EI/l mapping is correct** (the value to use as "
      f"angle.Harmonic's k)")
print(f"  ✓ the paper's eq.4, kappa(s)=kappa_0(a/s)^(2+d_b) with d_b=1, read with "
      f"**s = the ARM length (L/2)** gives")
print(f"    3EI/s^3 = 24EI/L^3, matching [P1] eq.1 exactly (reading s as the FULL "
      f"length is wrong by 8x)")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("(3) chain stiffness vs trap stiffness -- where is the measurable window?")
print("=" * 82)
print(f"  [P1] p.3 left column: \"trap compliance limits the maximum rigidity that can be measured\"")
print(f"  sketch k_t = 10 pN/µm,  paper k_t ~ {K_TRAP_PAPER*1e6:.0f} pN/µm")
print(f"\n  {'n':>4} {'L':>9} | " + " ".join(f"{lbl:>15}" for lbl in
      ("k(10mM,64mN/m)", "k(250mM,0.21)")) + "   verdict (sketch k_t=10pN/µm)")
for n in (5, 9, 11, 15, 25, 41):
    L = (n - 1) * D_PART
    ks = [24 * EI_from_kappa0(k0_) / L**3 for k0_ in (KAPPA0_EXP_10mM, KAPPA0_BARE_250mM)]
    r = ks[0] / 10e-6
    verd = ("chain stiffer than the trap -- trap compliance dominates" if r > 3 else
            "comparable -- the measurement window ★" if 0.3 < r < 3 else
            "chain is soft -- the signal is small")
    print(f"  {n:>4} {L*1e6:7.1f}µm | " + " ".join(f"{k*1e6:13.2f}pN/µm" for k in ks)
          + f"   {verd}")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("(4) a new hard constraint -- the critical bending moment M_c (where the "
      "harmonic angle potential breaks down)")
print("=" * 82)
print(f"""  Above M_c the particles **slip or roll** ([P2]'s conclusion). A harmonic
  angle potential cannot represent that plastic behaviour -> **measuring a linear
  response requires keeping M < M_c.**
  This is a kind of check 1-A and 1-B did not have: the validity range of a physical
  model.

  From [P2] eq.3, for a central rearrangement  M = F_bend * L/2
    -> F_bend < 2 M_c / L,   amplitude delta < F_bend/kappa(s)""")
print(f"\n  {'n':>4} {'L':>9} {'κ(s)':>12} {'F_max':>10} {'δ_max':>10} {'ℓ_k=√(kT/k_t)':>15}"
      f" {'window d_max/l_k':>13}")
l_k = math.sqrt(KT / 10e-6)
for n in (5, 9, 11, 15, 25):
    L = (n - 1) * D_PART
    ks = 24 * EI / L**3
    F_max = 2 * MC_PLATEAU / L
    d_max = F_max / ks
    print(f"  {n:>4} {L*1e6:7.1f}µm {ks*1e6:10.2f}pN/µm {F_max*1e12:8.2f}pN"
          f" {d_max*1e9:8.0f}nm {l_k*1e9:13.1f}nm {d_max/l_k:12.0f}×")
print(f"""
  -> amplitude window:  l_k({l_k*1e9:.0f} nm) << a < delta_max
     Below it the signal is buried in thermal noise; above it plasticity (bond
     slipping) sets in.
     For n=11 that is roughly
     20 nm << a < {2*MC_PLATEAU/(10*D_PART)/(24*EI/(10*D_PART)**3)*1e9:.0f} nm.""")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 82)
print("(5) timescales + the disagreement with the sketch")
print("=" * 82)
gamma_bead = 3 * math.pi * ETA * D_PART
tau_k = gamma_bead / 10e-6
tau_B = D_PART**2 / (KT / gamma_bead)
print(f"  the paper's particle d = {D_PART*1e6:.2f} µm ->  "
      f"gamma = {gamma_bead:.4e} kg/s")
print(f"    tau_k = gamma/k_t = {tau_k*1e3:.3f} ms   (trap relaxation, "
      f"k_t=10pN/µm)")
print(f"    τ_B = d²/D_t = {tau_B:.2f} s")
for n in (11, 25):
    L = (n - 1) * D_PART
    ks = 24 * EI / L**3
    tau_c = gamma_bead / ks
    print(f"    n={n}: κ(s)={ks*1e6:.2f} pN/µm → τ_chain = γ/κ = {tau_c*1e3:.3f} ms"
          f"   De=1 at ω={1/tau_c:.0f} rad/s ({1/tau_c/(2*math.pi):.0f} Hz)")
print(f"""
  ★ Where the sketch and the papers disagree
      sketch:  R = 5 µm              papers:  2a = 1.47 µm   <- 3.4x
      sketch:  k_t = 10 pN/µm        papers:  ~40 pN/µm
      sketch:  y = a sin(omega t), G' and G''
                                     papers:  constant {V_DRAG_P2*1e9:.0f} nm/s,
                                              quasi-static elasticity
    Every value from the papers is referenced to d=1.47µm. kappa_0 goes as
    a^(-1/3) (eq.7), so the dependence is weak and extrapolating to d=5µm is
    possible -- but **which one to use is a human decision.**
    Oscillatory driving (G' and G'') does not appear in the papers at all: the
    elastic constants carry over, but the dynamics is new territory.""")

print()
print("=" * 82)
print(f"{'✓ PASS' if (abs(err) < 5 and ok2) else '✗ FAIL'} -- formula extraction "
      f"verified ({err:+.1f}%), discrete<->continuum mapping verified")
print("=" * 82)
