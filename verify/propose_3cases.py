"""Derive the proposed values for the three blocked cases -- every number destined
for a system.yaml computed in one place.

User instruction (2026-08-04): "fill it in with proposals and show me the report
first."
-> Values coming from a sketch, a paper or a user confirmation are marked tier
0/1/2; values I chose myself are marked **tier 3**.

    $PY scratch/propose_3cases.py
"""
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from scipy import integrate  # noqa: E402

ETA = 0.851e-3
KT = 1.380649e-23 * 300
PI = math.pi
HEX_NN = math.sqrt(2 / math.sqrt(3))


def hr(t):
    print(f"\n{'=' * 84}\n{t}\n{'=' * 84}")


# ══════════════════════════════════════════════════════════════════════
hr("A.  abp-rod-2d-run-tumble")
# ── Ellipsoid friction (same calculation as verify/perrin_friction.py, whose
#    sphere limit is already verified) ──
def chi(axes):
    a = np.asarray(axes, float)
    dl = lambda l: math.sqrt((a[0]**2 + l) * (a[1]**2 + l) * (a[2]**2 + l))
    sc = float(a.max())**2
    sub = lambda f: integrate.quad(lambda u: f(sc * (1 / u - 1)) * sc / u**2,
                                   1e-14, 1.0, limit=400)[0]
    return sub(lambda l: 1 / dl(l)), [sub(lambda l, i=i: 1 / ((a[i]**2 + l) * dl(l)))
                                      for i in range(3)]


a1, a2 = 1.0e-6, 0.25e-6                       # user-confirmed: major axis 2µm, minor 500nm
chi0, chis = chi((a1, a2, a2))
zt = [16 * PI * ETA / (chi0 + x**2 * c) for x, c in zip((a1, a2, a2), chis)]
zr_perp = (16 * PI * ETA * (a1**2 + a2**2)) / (3 * (a1**2 * chis[0] + a2**2 * chis[1]))
d_eq = 2 * (a1 * a2 * a2) ** (1 / 3)
gamma_bar = 2.0 / (1 / zt[0] + 1 / zt[1])      # 2D in-plane harmonic mean
D_bar = KT / gamma_bar
D_r = KT / zr_perp
tau_r = 1 / D_r
TAU_TUMBLE = 0.5                                # user-confirmed
V = 5e-6                                        # the sketch's upper bound
tau_eff = 1 / (1 / tau_r + 1 / TAU_TUMBLE)      # uniform random reorientation
                                                # (proposed, tier 3)
l_p = V * tau_eff
tau_B_abp = d_eq**2 / D_bar
tau_v = d_eq / V
RHO_P = 1050.0                                  # proposed, tier 3 (polystyrene-like)
m_ell = RHO_P * (4 * PI / 3) * a1 * a2 * a2
tau_p_abp = m_ell / gamma_bar

print(f"  fixed     semi-axes (1.0, 0.25, 0.25) µm . p=4 . water@300K . "
      f"tau_tumble=0.5s . v<=5µm/s")
print(f"  derived   d_eq = {d_eq*1e6:.4f} µm   gammabar(2D) = {gamma_bar:.4e} kg/s"
      f"   D̄ = {D_bar*1e12:.4f} µm²/s")
print(f"         γ_r,z = {zr_perp:.4e} kg·m²/s   D_r = {D_r:.4f} 1/s   τ_r = {tau_r:.4f} s")
print(f"         τ_B(d_eq) = {tau_B_abp:.4f} s   τ_v = d_eq/v = {tau_v:.4f} s"
      f"   τ_p = {tau_p_abp*1e6:.4f} µs")
print(f"  proposed  tumble = uniform random reorientation (2D)  -> "
      f"tau_eff = {tau_eff:.4f} s")
print(f"            l_p = v*tau_eff = {l_p*1e6:.4f} µm = {l_p/d_eq:.3f} d_eq = "
      f"{l_p/(2*a1):.3f} body lengths")
print(f"         Pe = v d_eq/D̄ = {V*d_eq/D_bar:.3f}    D_r* = D_r τ_B = {D_r*tau_B_abp:.3f}")
N_ABP, L_ABP_D = 1000, 32.0                     # proposed, tier 3
L_abp = L_ABP_D * d_eq
dt_abp = 1e-2 * min(tau_v, tau_r, TAU_TUMBLE)
print(f"            N = {N_ABP} (independent ensemble, same approach as 1-A) . "
      f"L = {L_ABP_D:.0f} d_eq"
      f" = {L_abp*1e6:.2f} µm")
print(f"  dt        fastest relevant scale = tau_v = {tau_v:.4f} s  ->  "
      f"dt = {dt_abp*1e3:.4f} ms"
      f" = {dt_abp/tau_B_abp:.3e} τ_B")
T_ABP = 200 * max(tau_r, TAU_TUMBLE)
print(f"  T_obs  {T_ABP:.0f} s = {T_ABP/max(tau_r,TAU_TUMBLE):.0f}×max(τ_r,τ_tumble)"
      f"  →  {T_ABP/dt_abp:.2e} steps")
print(f"  checks    l_p/(L/4) = {l_p/(L_abp/4):.4f}  .  "
      f"tau_p/tau_v = {tau_p_abp/tau_v:.3e}"
      f"  ·  dt·D_r = {dt_abp*D_r:.3e}  ·  dt/τ_tumble = {dt_abp/TAU_TUMBLE:.3e}")

# ══════════════════════════════════════════════════════════════════════
hr("B.  trap-drag-2d-hex300")
d_td = 5.0e-6
gam_td = 3 * PI * ETA * d_td
Dt_td = KT / gam_td
tauB_td = d_td**2 / Dt_td
k_t = 1e-5                                      # 10 pN/µm (from the sketch)
tau_k_td = gam_td / k_t
l_k = math.sqrt(KT / k_t)
v_x = 0.5e-6
A_TD, PHI_TD, N_TD = 100.0, 0.35, 300           # A and phi proposed, tier 3 (the
                                                # combination 1-B confirmed gives a
                                                # hexagonal crystal)
a_mean_td = d_td * math.sqrt(PI / (4 * PHI_TD))
a_nn_td = HEX_NN * a_mean_td
L_td = a_mean_td * math.sqrt(N_TD)
r_c_td = 5 * a_mean_td
Gamma_td = A_TD / (a_mean_td / d_td) ** 3
# Pair local stiffness (same definition as 1-B): U''(r) = 12A kT d^3/r^5 + WCA
r_min_td = a_nn_td / d_td - 3 * math.sqrt(2 / (3 * (12 * A_TD / (a_nn_td/d_td)**5
                                                   - 3 * A_TD / (a_nn_td/d_td)**5)))
Upp_pair = 12 * A_TD / r_min_td**5              # [kT/d²]
tau_int_td = tauB_td / Upp_pair
print(f"  fixed     d=5µm (convention) . k_t=10pN/µm . v_x=0.5µm/s . N~300 "
      f"(all from the sketch)")
print(f"            pair = A/r^3 + WCA (user-confirmed 2026-08-04)")
print(f"  proposed  A = {A_TD:.0f} . phi = {PHI_TD} -> Gamma = {Gamma_td:.2f}"
      f"   (1-B confirmed this Gamma gives a hexagonal crystal, psi6=0.885)")
print(f"         a_mean = {a_mean_td*1e6:.3f} µm   a_NN = {a_nn_td*1e6:.3f} µm"
      f"   L = {L_td*1e6:.2f} µm = {L_td/d_td:.1f} d")
print(f"            r_c = 5 a_mean = {r_c_td*1e6:.2f} µm    minimum-image margin "
      f"{(L_td/2)/r_c_td:.2f}×")
print(f"  derived   gamma = {gam_td:.4e} kg/s  D_t = {Dt_td*1e12:.4f} µm^2/s  "
      f"tau_B = {tauB_td:.2f} s")
print(f"         τ_k = γ/k_t = {tau_k_td*1e3:.3f} ms   ℓ_k = {l_k*1e9:.2f} nm")
print(f"            tau_int(pair, r_min={r_min_td:.3f}d) = {tau_int_td*1e3:.3f} ms")
print(f"         τ_v = d/v_x = {d_td/v_x:.1f} s   Δr_ss = γv/k = {gam_td*v_x/k_t*1e9:.3f} nm")
dt_td = 1e-2 * min(tau_k_td, tau_int_td)
binder = "the trap, tau_k" if tau_k_td < tau_int_td else "the pair, tau_int"
print(f"  dt        ★ two stiffnesses compete -- the faster is {binder}")
print(f"         dt = 10⁻²×{binder.split()[1]} = {dt_td*1e6:.2f} µs = {dt_td/tauB_td:.3e} τ_B")
T_TD = L_td / v_x
print(f"  T_obs     box traverse = L/v = {T_TD:.0f} s  ->  {T_TD/dt_td:.2e} steps")
print(f"  checks    SNR = dr_ss/l_k = {gam_td*v_x/k_t/l_k:.4f}  <- the signal is 1/10 "
      f"of the noise (the standing warning)")

# ══════════════════════════════════════════════════════════════════════
hr("C.  chain-bend-2d-oscill")
d_cb = 1.47e-6                                  # papers [P1][P2] -- disagrees with
                                                # the sketch's 5µm
gam_cb = 3 * PI * ETA * d_cb
Dt_cb = KT / gam_cb
tauB_cb = d_cb**2 / Dt_cb
tau_k_cb = gam_cb / k_t
KAPPA0 = 64e-3                                  # [P1] 10 mM MgCl₂, tier 2
a_cb = d_cb / 2
EI = KAPPA0 * a_cb**3 / 3
kappa_theta = EI / d_cb
# ★ Taking n=11 gave delta_max = 74 nm, which nearly closed the amplitude window
#   (just above l_k=20nm).
#   delta_max = M_c L^2/(12 EI) ~ n^2, so lengthening the chain opens the window.
#   [P1] Fig.4 uses a 25-particle chain -> n=25 is proposed.
N_BEADS = 25                                    # [P1] Fig.4 (25-particle aggregate)
L_cb = (N_BEADS - 1) * d_cb
k_end = 24 * EI / L_cb**3
k_center = 2 * k_end
MC = 35e-18                                     # [P1] tier 2
# Bond (radial) stiffness: keep the thermal extension below 0.1% of the diameter
K_BOND = KT / (1e-3 * d_cb) ** 2
print(f"  papers    d = {d_cb*1e6:.2f} µm . E=3100MPa . nu=0.4 . "
      f"kappa_0={KAPPA0*1e3:.0f} mN/m"
      f" · M_c={MC*1e18:.0f} pN·µm  (tier 2)")
print(f"  derived   EI = {EI:.4e} N.m^2   kappa_theta = EI/l = {kappa_theta:.4e} J"
      f" = {kappa_theta/KT:.3e} kT")
print(f"         n={N_BEADS} → L = {L_cb*1e6:.2f} µm  κ_end = {k_end*1e6:.2f} pN/µm"
      f"  κ_center = {k_center*1e6:.2f} pN/µm")
print(f"         γ = {gam_cb:.4e} kg/s  τ_B = {tauB_cb:.3f} s  τ_k = {tau_k_cb*1e3:.3f} ms")
print(f"  proposed  k_bond = {K_BOND*1e6:.0f} pN/µm (thermal extension <= 0.1% d) "
      f". tier 3")

# ★ dt comes from the largest eigenvalue of the stiffness matrix (the fastest mode
#   sets dt)
n = N_BEADS
B = np.zeros((n - 2, n))
for i in range(n - 2):
    B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
B /= d_cb
A_bend = kappa_theta * (B.T @ B)                 # transverse bending
G = np.zeros((n - 1, n))
for i in range(n - 1):
    G[i, i], G[i, i + 1] = -1.0, 1.0
G /= 1.0
A_bond = K_BOND * (G.T @ G)                      # longitudinal stretching
lam_bend = float(np.linalg.eigvalsh(A_bend).max())
lam_bond = float(np.linalg.eigvalsh(A_bond).max())
lam = max(lam_bend, lam_bond)
tau_fast = gam_cb / lam
dt_cb = 1e-2 * tau_fast
print(f"\n  ★ dt is taken from the **largest eigenvalue of the stiffness matrix** "
      f"(the fastest mode sets it)")
print(f"    lambda_max(bending)    = {lam_bend:.4e} N/m  -> "
      f"tau = {gam_cb/lam_bend*1e6:.4f} µs")
print(f"    lambda_max(stretching) = {lam_bond:.4e} N/m  -> "
      f"tau = {gam_cb/lam_bond*1e6:.4f} µs")
print(f"    governing: {'bending' if lam_bend > lam_bond else 'stretching'}  ->  "
      f"tau_fast = {tau_fast*1e6:.4f} µs")
print(f"    dt = 10⁻² τ_fast = {dt_cb*1e9:.2f} ns = {dt_cb/tauB_cb:.3e} τ_B")
print(f"    note: the collective bending mode (kappa_center) is "
      f"{gam_cb/k_center*1e6:.1f} µs -- "
      f"{(gam_cb/k_center)/tau_fast:.0f}x slower than the fastest mode")

# Amplitude window
F_max = 2 * MC / L_cb
d_max = F_max / k_end
print(f"\n  amplitude window  l_k = {l_k*1e9:.1f} nm << a < "
      f"delta_max = {d_max*1e9:.0f} nm"
      f"   (M<M_c: F<{F_max*1e12:.2f} pN)")
print(f"  {'n':>4} {'L':>9} {'kap_end':>11} {'del_max':>9}   amplitude window")
for nb in (11, 15, 25, 41):
    Lb = (nb - 1) * d_cb
    ke = 24 * EI / Lb**3
    dm = MC * Lb**2 / (12 * EI)
    print(f"  {nb:>4} {Lb*1e6:7.1f}µm {ke*1e6:9.2f}pN/µm {dm*1e9:7.0f}nm   "
          f"{'closed (just above l_k)' if dm < 4*l_k else f'20nm << a < {dm*1e9:.0f}nm'}")
A_AMP = 200e-9
print(f"  proposed   a = {A_AMP*1e9:.0f} nm  -> a/l_k = {A_AMP/l_k:.1f} (SNR)"
      f"  a/delta_max = {A_AMP/d_max:.3f} "
      f"(linear margin {d_max/A_AMP:.1f}x) . tier 3")

# Frequency window + cost
print(f"\n  frequency  which tau is the tau in De = omega*tau? -- the trap (tau_k) "
      f"vs the chain's collective mode")
tau_chain = gam_cb / k_center
for lbl, tau_ in (("tau_k(trap)", tau_k_cb), ("tau_chain(collective)", tau_chain)):
    print(f"    {lbl:<16} = {tau_*1e3:8.4f} ms  →  De=1 at ω = {1/tau_:9.1f} rad/s"
          f" ({1/tau_/(2*PI):8.1f} Hz)")
OMEGA_LO, OMEGA_HI = 0.1 / tau_k_cb, 10 / tau_chain
print(f"  proposed   omega sweep {OMEGA_LO:.0f} to {OMEGA_HI:.0f} rad/s "
      f"(covers De 0.1-10) . tier 3")
n_cyc = 100
T_cb = n_cyc * 2 * PI / OMEGA_LO
print(f"\n  ★ cost  {n_cyc} cycles at the lowest omega = {T_cb:.2f} s  ->  "
      f"{T_cb/dt_cb:.2e} steps")
print(f"    1-B's measured throughput (N=400, 79 neighbours) was ~6600 steps/s. "
      f"Here N={N_BEADS} is")
print(f"    far smaller so it will be faster, but the step COUNT itself is "
      f"{T_cb/dt_cb:.1e}.")
print(f"    -> **using the papers' stiffness as-is makes direct BD expensive.** "
      f"Because kappa_theta = {kappa_theta/KT:.1e} kT makes")
print(f"      the chain thermally rigid, so the fastest mode -- irrelevant to the "
      f"measurement -- sets dt.")
print(f"    options: (a) accept it  (b) lower kappa_0 (surfactant condition, "
      f"[P2] Fig.4) and sweep")
print(f"             (c) measure only high omega and replace low omega with the "
      f"quasi-static limit")

hr("summary -- the tiers to record in system.yaml")
print("""  tier 0  written on the sketch      k_t . v_x . N~300 . tau_R=0.5s . v<=5µm/s
                                     . A={0.1,1,10,100}
  tier 1  user-confirmed / inherited ellipsoid shape . medium (water 300K)
                                     . d=5µm convention . pair=A/r^3+WCA
  tier 2  extracted from papers      d=1.47µm . E . nu . kappa_0=64mN/m
          (unverified)               . M_c=35pN.µm
  tier 3  ★ my proposals            tumble angle distribution . N_abp=1000 . L
                                     . A=100, phi=0.35 (trap-drag)
                                     . n_beads=25 . k_bond . amplitude 200nm
                                     . omega sweep range . rho_p""")
