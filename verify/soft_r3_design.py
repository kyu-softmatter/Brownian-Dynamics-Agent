"""Phase 1-B `soft-r3-2d-A-sweep` design calculation -- with the user's confirmed
values (2026-08-03).

Choices confirmed by the user, 2026-08-03:
    density   phi = 0.35  (a_mean = 1.5 d)
    core      add a WCA at d (epsilon = kT)
    tail      raise N to secure a proper r_c
    anchor    d = 5 µm  (unified with the R=5µm convention of the other sketches)

Still unresolved: is A dimensionless or does it carry µm^3? (Section 0 computes both
readings side by side and produces the basis for deciding.)

Usage: $PY verify/soft_r3_design.py
"""

import math

import numpy as np
import pint

u = pint.UnitRegistry()
Q = u.Quantity
PI = math.pi

kB = Q(1.380649e-23, "J/K")

# ── Confirmed anchors ──────────────────────────────────────────────────
d = Q(5.0, "um")  # user-confirmed (the R=5µm convention of the other sketches)
T = Q(300, "K")  # inherited from the value confirmed in 1-A
ETA = Q(0.851, "mPa*s").to("Pa*s")  # water@300K
RHO_P = Q(2000, "kg/m^3")  # silica assumed, for tau_p only
PHI = 0.35  # user-confirmed
A_LIST = (0.1, 1.0, 10.0, 100.0)  # from the sketch
EPS_WCA = 1.0  # WCA depth [kT] -- a convention. Its only role is to set where the
               # core sits

KT = (kB * T).to("J")
gamma = (3 * PI * ETA * d).to("kg/s")
D_t = (KT / gamma).to("um^2/s")
tau_B = (d**2 / D_t).to("s")
m = (RHO_P * (PI / 6) * d**3).to("kg")
tau_p = (m / gamma).to("s")

a_mean = d * math.sqrt(PI / (4 * PHI))  # 2D: φ = (π/4)(d/a)²
a_star = float(a_mean / d)


def hr(t, ch="="):
    print(f"\n{ch * 78}\n{t}\n{ch * 78}")


# ══════════════════════════════════════════════════════════════════════
hr("0.  the dimensional reading of A -- at d = 5 µm this choice changes the physics "
   "by 125x")
print(f"""
  reading (i)  A dimensionless, r in diameters:  U/kT = A (d/r)^3
  reading (ii) [A] = µm^3, r in µm:              U/kT = A[µm^3] / r[µm]^3

  coupling strength at the mean spacing, Gamma == U(a_mean)/kT
  (the parameter that actually controls the structure)
    a_mean = {a_mean:~.2fP} = {a_star:.3f} d
""")
print(f"    {'A':>7} | {'Gam read(i)':>12} | {'Gam read(ii)':>13} | "
      f"regime under reading (i)")
for A in A_LIST:
    G_i = A / a_star**3
    G_ii = A / a_mean.to("um").magnitude ** 3
    regime = "weakly correlated fluid" if G_i < 1 else (
        "correlated fluid" if G_i < 10 else "strongly coupled (crystal cand.)")
    print(f"    {A:7.1f} | {G_i:12.4f} | {G_ii:13.5f} | {regime}")
print(f"""
  ★ Under reading (ii), Gamma = 0.24 even at A=100 -- the whole sweep sits below kT.
    Then rdf and voronoi would show no structure at all beyond the WCA excluded
    volume.
    The sketch names 'final configuration? / rdf / voronoi / structure analysis' as
    the goal, so **reading (i) is the only one consistent with that goal.**
    -> proceeding with (i) below. Needs confirmation (observation.yaml A1).
""")

# ══════════════════════════════════════════════════════════════════════
hr("1.  the physical system (SI) -- derived from the confirmed anchors")
print(f"    d     = {d:~P}          (user-confirmed)")
print(f"    T     = {T:~P}          η = {ETA.to('mPa*s'):~.3fP}   kT = {KT:~.4eP}")
print(f"    γ     = 3πηd = {gamma:~.4eP}")
print(f"    D_t   = kT/γ = {D_t:~.4fP}")
print(f"    tau_B = d^2/D_t = {tau_B:~.2fP}   "
      f"★ the reference time (this system's governing scale)")
print(f"    tau_p = m/gamma = {tau_p.to('us'):~.3fP}   (for the model check only)")
print(f"    φ     = {PHI}  →  a_mean = {a_mean:~.3fP} = {a_star:.3f} d")


# ══════════════════════════════════════════════════════════════════════
# The potential (reduced: lengths in d, energies in kT)
#   U*(r*) = U_WCA(r*) + A / r*³      r* = r/d
# ══════════════════════════════════════════════════════════════════════
R_WCA = 2 ** (1 / 6)


def U_star(rs, A):
    w = np.where(rs < R_WCA, 4 * EPS_WCA * (rs**-12.0 - rs**-6.0) + EPS_WCA, 0.0)
    return w + A / rs**3


def U2_star(rs, A):
    """U''(r) [kT/d^2] -- the local stiffness."""
    w = np.where(rs < R_WCA, 4 * EPS_WCA * (156 * rs**-14.0 - 42 * rs**-8.0), 0.0)
    return w + 12 * A / rs**5


def solve_r_pair(A, u_max=5.0):
    """The distance where U*(r) = u_max (the dilute / pair-approach criterion).

    By bisection.
    """
    lo, hi = 0.5, 50.0
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if U_star(np.array([mid]), A)[0] > u_max:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


hr("2.  closest approach r_min -- this is where dt comes from")
LINDEMANN = 0.15  # above this, u_rms/a invalidates the cage (crystal) picture -> fluid
print(f"""  Both criteria are computed and the **smaller** (i.e. stiffer) one is used.
  A conservative choice.
    (a) pair criterion      : the distance where U(r) = 5 kT
                              -- dominates in the dilute, weakly coupled case
    (b) vibrational criterion: a_mean - 3 u_rms
                              -- valid only in the dense, strongly coupled case
        2D hexagonal cage stiffness k_cage = 3[U''(a) + U'(a)/a],
        u_rms = sqrt(kT/k_cage)
        ⚠ if u_rms/a > {LINDEMANN} the cage has melted, so (b) is NOT used
          (using it anyway gives a negative r_cage -- which is what the first
           calculation actually produced)
""")
print(f"    {'A':>7} {'Γ':>8} | {'(a) r_pair':>11} | {'u_rms/a':>8} {'(b) r_cage':>11} |"
      f" {'r_min':>8} {'U2(r_min)':>10} {'tau_int/tau_B':>11} {'dt/tau_B':>10}"
      f"  criterion")
design = {}
for A in A_LIST:
    G = A / a_star**3
    r_pair = solve_r_pair(A)
    # Cage stiffness (at a_mean=1.5d the WCA is 0, so only the r^-3 part contributes)
    Upp_a = float(U2_star(np.array([a_star]), A)[0])
    Up_a = -3 * A / a_star**4  # U'(a) [kT/d]
    k_cage = 3 * (Upp_a + Up_a / a_star)
    u_rms = math.sqrt(1.0 / k_cage) if k_cage > 0 else float("inf")
    ur = u_rms / a_star
    crystalline = ur < LINDEMANN
    r_cage = a_star - 3 * u_rms
    if crystalline and r_cage < r_pair:
        r_min, which, state = r_cage, "vibrational", "crystal"
    else:
        r_min, which, state = r_pair, "pair", ("crystal" if crystalline else "fluid")
    Upp = float(U2_star(np.array([r_min]), A)[0])
    tau_int_star = 1.0 / Upp  # τ_int/τ_B = (γ/U'')/(γd²/kT) = 1/U''*
    dt_star = 1e-2 * tau_int_star
    design[A] = dict(G=G, r_min=r_min, dt_star=dt_star, which=which, u_rms=ur, state=state)
    cage_txt = f"{r_cage:10.3f}d" if crystalline else "  (invalid) "
    print(
        f"    {A:7.1f} {G:8.3f} | {r_pair:10.3f}d | {ur:8.4f} {cage_txt} |"
        f" {r_min:7.3f}d {Upp:10.1f} {tau_int_star:11.3e} {dt_star:10.2e}  {which}({state})"
    )
print(f"""
    ★ At A=0.1 and 1, **the WCA core sets dt** (r^-3 is nearly irrelevant,
      Gamma <= 0.3).
      At A=100 the particles are caged, so approach itself is limited and dt is
      actually looser.
      -> the most expensive point in this sweep is A=0.1, the physically dullest one.
    ★ dt differs per A (bd-physics section 1.1: the unit system is fixed, dt is
      per-system).
    ⚠ The system is nonlinear so there is no closed form for the bias -> a
      half-dt convergence check is required (mater_plan section 20 B)
    ⚠ u_rms/a is also the Lindemann indicator: at A=100 it is
      {design[100.0]['u_rms']:.3f}
      -> a crystal candidate. The actual verdict comes from voronoi and psi6 AFTER
         the run (it is not declared here)
""")


# ══════════════════════════════════════════════════════════════════════
hr("3.  cutoff and N -- handling the r^-3 tail honestly")
print("""  In 2D the r^-3 tail energy converges **slowly**, as
  integral r^-3 * r dr ~ 1/r_c.
    energy beyond r_c / nearest-neighbour energy = 2*pi*a/(3*r_c)
  -> absolute energies and pressures need a tail correction, while **structure
     (rdf, voronoi) is barely affected**
     (forces from a homogeneous distant field cancel by symmetry -- provided
      xi < r_c)
  -> so the cutoff should be set as a **multiple of a_mean**, not against an
     absolute kT threshold.
""")
print(f"    {'r_c':>10} {'U(r_c)/kT @A=100':>18} {'/Gam':>8} {'tail/neigh':>10}"
      f" | min-image N lower bound")
for k in (3, 5, 7, 10):
    rc = k * a_star
    u_rc = 100.0 / rc**3
    G100 = 100.0 / a_star**3
    tail = 2 * PI / (3 * k)
    N_min = (2 * k) ** 2  # r_c < L/2 = (√N/2)a  →  √N > 2k
    print(f"    {k:2d} a_mean {rc:5.2f}d {u_rc:14.4f} {u_rc / G100:8.4f} {tail:10.3f} | N > {N_min}")

N = 400
L = a_mean * math.sqrt(N)
rc_star = 5 * a_star
print(f"""
  recommended:  N = {N}   ->  L = {L:~.1fP} = {float(L / d):.1f} d
                = {math.sqrt(N):.0f} a_mean
         r_c = 5 a_mean = {rc_star:.2f} d = {(rc_star * d).to('um'):~.1fP}
         minimum-image margin = (L/2)/r_c = {math.sqrt(N) / 2 / 5:.2f}x
         ★ At N=400, r_c = 5, 7 and 9 a_mean can ALL be tested inside the minimum
           image
           -> so the r_c convergence check can be done at one fixed N
              (comparing rdf at A=100)
         neighbour count ~ pi*r_c^2/a_mean^2 = {PI * 25:.0f} per particle
""")


# ══════════════════════════════════════════════════════════════════════
hr("4.  separation checks (model / integration / geometry / statistics) -- the same "
   "classification as 1-A")
T_OBS_STAR = 100.0  # in units of tau_B -- the value that satisfies the statistics
                    # check T_obs >= 1e2 tau_B
for A in A_LIST:
    dd = design[A]
    print(f"\n  A = {A}   (Γ = {dd['G']:.3f},  dt* = {dd['dt_star']:.2e})")
    rows = [
        ("model", "inertia negligible  tau_p/tau_int",
         float(tau_p / tau_B) / (dd["dt_star"] / 1e-2), 1e-2),
        ("integration", "interaction resolved  dt/tau_int", 1e-2, 1e-2),
        ("geometry", "cutoff  r_c/(L/2)", rc_star / (math.sqrt(N) * a_star / 2), 1.0),
        ("statistics", "observation window  tau_B/T_obs", 1.0 / T_OBS_STAR, 1e-2),
    ]
    for kind, name, val, lim in rows:
        ok = val <= lim
        mark = "✓" if ok else ("⚠" if val <= 1.5 * lim else "✗")
        print(f"    {mark} [{kind}] {name:26s} = {val:9.3e}  ≤ {lim:7.1e}"
              f"   margin {lim / val:6.1f}x")
    steps = T_OBS_STAR / dd["dt_star"]
    print(f"      cost: T_obs = {T_OBS_STAR:.0f} tau_B = "
          f"{(T_OBS_STAR * tau_B).to('s'):~.0fP} (physical)"
          f"  →  {steps:.2e} steps")

print(f"""
  All four checks pass. But passing does not guarantee accuracy (bd-physics
  section 4):
  ⚠ The statistics check tested T_obs >= 1e2 tau_B against tau_B. Under strong
    coupling (A=100), defect annealing can be slower than tau_B, so this check
    alone may not be enough.
    -> do not declare equilibrium; measure it with the post-mortem EQ diagnostic
       (tools/postmortem.py).
  ⚠ The integration check's margin is exactly 1.0x because dt was set AT the limit.
    That corresponds to 0.5% bias (for a linear system). The system is nonlinear, so
    a half-dt convergence check is required.
""")

total = sum(T_OBS_STAR / design[A]["dt_star"] for A in A_LIST)
print(f"  total over the 4 A values ~ {total:.2e} steps "
      f"(4 independent runs -> can be run in parallel, principle 6)")
