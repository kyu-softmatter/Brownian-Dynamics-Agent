"""The remaining 4 intake/ cases -- scale tables and a draft non-dimensionalization
(2026-08-03).

Follows CLAUDE.md rule 1 (dimensions first) and rule 3 (every value carries a
provenance).
  READ   = a value written on the sketch (tier 0)
  ANCHOR = a value absent from the sketch that had to be filled in (tier stated --
           mandatory, because the non-dimensionalization is otherwise undetermined)
  OPEN   = a value that cannot be filled in (not invented; for a human to confirm)

Usage: $PY scratch/intake_scales.py
"""

import math

import pint

u = pint.UnitRegistry()
Q = u.Quantity
PI = math.pi

kB = Q(1.380649e-23, "J/K")

# ══════════════════════════════════════════════════════════════════════════
# Shared anchors -- none of the 4 sketches states a medium or a temperature.
#   In 1-A (trap-2d-5um) the user confirmed "water, 300 K", so that convention is
#   carried over. They are consecutive sketches in one notebook, so reading them as
#   the same experimental system is reasonable -- but the sketches themselves give
#   no basis for it, so this is marked tier 1 (inherited convention).
# ══════════════════════════════════════════════════════════════════════════
T = Q(300, "K")  # ANCHOR tier1 -- inherited from the value confirmed in 1-A
ETA = Q(0.851, "mPa*s").to("Pa*s")  # ANCHOR tier0 (handbook water@300K) / the medium
                                    # assumption itself is tier1
RHO_P = Q(2000, "kg/m^3")  # ANCHOR tier3 -- silica assumed. Used only for tau_p
                           # (the model check)
KT = (kB * T).to("J")


def stokes(d):
    """Derive sphere-based material properties (bd-physics §2). d = diameter."""
    d = d.to("m")
    gamma = (3 * PI * ETA * d).to("kg/s")
    D_t = (KT / gamma).to("um^2/s")
    tau_B = (d**2 / D_t).to("s")
    m = (RHO_P * (PI / 6) * d**3).to("kg")
    tau_p = (m / gamma).to("s")
    # Rotation (holds for a sphere only -- do not use for non-spherical particles)
    D_r = (KT / (PI * ETA * d**3)).to("1/s")
    tau_r = (1 / D_r).to("s")
    return dict(d=d, gamma=gamma, D_t=D_t, tau_B=tau_B, m=m, tau_p=tau_p, D_r=D_r, tau_r=tau_r)


def hr(title, ch="="):
    print(f"\n{ch * 78}\n{title}\n{ch * 78}")


def chk(name, value, limit, kind="≤", note=""):
    """One separation-check line. Reports the margin alongside (bd-physics §4)."""
    ok = value <= limit if kind == "≤" else value < limit
    margin = limit / value if value > 0 else float("inf")
    mark = "✓" if ok else "✗"
    print(f"    {mark} {name:26s} = {value:9.3e}  {kind} {limit:7.1e}   margin {margin:8.1f}x  {note}")
    return ok


print(f"shared anchors:  T = {T:~P}   eta = {ETA.to('mPa*s'):~.3fP}(water@300K)   "
      f"kT = {KT:~.4eP}")
print(f"           rho_p = {RHO_P:~P} (silica assumed, for the tau_p check only)")
print("NOTE all 4 sketches omit the medium and temperature -> the 'water, 300 K' "
      "convention confirmed in 1-A is inherited")


# ══════════════════════════════════════════════════════════════════════════
# CASE 1-B ─ soft-r3-2d-A-sweep
#   READ : U_ij/kT = A/r_ij³ ,  A = 0.1, 1, 10, 100 ,  N = 100 ,  Lx = Ly, 2D
#   OPEN : the density (phi, or L, or a_mean) -- without it, A alone does not
#          determine the physics
#   OPEN : whether an excluded-volume core is present (with r^-3 alone, particles
#          interpenetrate for A <~ 1)
#   ANCHOR: d -- needed to pin down the dimensionless reading of A
# ══════════════════════════════════════════════════════════════════════════
hr("CASE 1-B  soft-r3-2d-A-sweep   <- the mater_plan §16 Phase 1-B target")

d_soft = Q(1.0, "um")  # ANCHOR tier3 (reasoning below)
s = stokes(d_soft)
print(f"""
READ    U_ij/kT = A/r_ij³   A ∈ {{0.1, 1, 10, 100}}   N = 100   Lx = Ly (2D)
        goal: final configuration / rdf / voronoi / structure analysis
ANCHOR  d = {d_soft:~P}  (tier 3)
        ★ Is A dimensionless? For U/kT = A/r^3 to hold dimensionally, either
            (i)  read r in units of the diameter, making A dimensionless
                 ->  U/kT = A (d/r)^3
            (ii) read r in µm, making [A] = µm^3
                 ->  U/kT = A[µm^3]/r[µm]^3
          With the d = 1 µm anchor, (i) and (ii) are **numerically identical**.
          -> choosing the anchor at 1 µm dissolves the ambiguity (which is why
             1 µm was chosen)
OPEN    the density (phi / L / a_mean) is not stated -- A alone does not determine
        the structure (demonstrated below)
OPEN    whether an excluded-volume core (WCA or similar) is present is not stated
""")

print(f"DERIVED (d = {d_soft:~P}, water@300K)")
print(f"    γ   = {s['gamma']:~.4eP}      D_t = {s['D_t']:~.4fP}")
print(f"    τ_B = d²/D_t = {s['tau_B']:~.4fP}   ★ the governing timescale of this system (contrast 1-A, where it was tau_k)")
print(f"    τ_p = m/γ    = {s['tau_p'].to('us'):~.3fP}  (for the model check only)")

# ── A length scale unique to this system: the distance at which U(r) = c*kT ──
print(f"\nnew length scale -- r_u(A,c): the distance where U(r)/kT = c, "
      f"= A^(1/3) c^(-1/3) d")
print(f"    {'A':>6} | {'r(U=kT)':>9} {'r(U=5kT)':>9} {'r(U=0.01kT)':>12}   reading")
for A in (0.1, 1.0, 10.0, 100.0):
    r1 = A ** (1 / 3)
    r5 = (A / 5) ** (1 / 3)
    rc = (A / 0.01) ** (1 / 3)
    tag = "interpenetrates without a core" if r5 < 1.0 else "contact blocked by r^-3 alone"
    print(f"    {A:6.1f} | {r1:8.3f}d {r5:8.3f}d {rc:11.2f}d   {tag}")
print("    -> for A <~ 5, r^-3 repulsion alone cannot prevent overlap "
      "(r(U=5kT) < d).")
print("      Whether an excluded-volume core exists must be confirmed -- if not, it "
      "has to be stated as a 'point-particle model'.")

# ── The real control parameter: the coupling strength at the mean spacing ──
print(f"\ndimensionless groups -- the structure is set by A combined with the "
      f"density, not by A alone")
print(f"    Γ ≡ U(a_mean)/kT = A (d/a_mean)³      a_mean = ρ^(-1/2) = L/√N  (2D)")
print(f"    φ = (π/4)(d/a_mean)²                  L = a_mean √N")
print(f"\n    {'a/d':>5} {'φ':>7} | " + " ".join(f"{'Γ(A=' + str(A) + ')':>12}" for A in (0.1, 1, 10, 100)))
for a_over_d in (1.2, 1.5, 2.0, 3.0, 4.3):
    phi = (PI / 4) / a_over_d**2
    row = " ".join(f"{A / a_over_d**3:12.3f}" for A in (0.1, 1, 10, 100))
    print(f"    {a_over_d:5.1f} {phi:7.3f} | {row}")
print("    -> at the same A, Gamma varies by more than 10x with the density. "
      "Without the density the sweep is not defined.")

# ── Cutoff vs minimum image: the limit of N=100 (trap 7 territory) ─────────
print(f"\ngeometry check -- r^-3 has a long tail. Can r_c < L/2 hold while still "
      f"respecting a truncation error u_c?")
print(f"    r_c = (A/u_c)^(1/3) d   ,   L/2 = (√N/2) a_mean")
print(f"    upper bound on the coupling that satisfies both at once:  "
      f"Gamma_max = N^(3/2) u_c / 8   (independent of A!)")
N_soft = 100
PHI_MAX = 0.9069  # 2D hexagonal close packing -- nothing can be denser
A_MIN = math.sqrt(PI / (4 * PHI_MAX))  # lower bound on a_mean/d = 0.934
for u_c in (0.01, 0.1):
    G_max = N_soft**1.5 * u_c / 8
    print(f"\n    u_c = {u_c:<5} (potential remaining at truncation)   N = {N_soft} ->  "
          f"Gamma_max(cutoff) = {G_max:.3f}")
    for A in (0.1, 1.0, 10.0, 100.0):
        rc = (A / u_c) ** (1 / 3)  # in units of d
        a_need = max(2 * rc / math.sqrt(N_soft), A_MIN)  # cutoff constraint OR
                                                        # close-packing constraint
        phi_max = (PI / 4) / a_need**2
        G_at = A / a_need**3
        bind = "cutoff" if 2 * rc / math.sqrt(N_soft) > A_MIN else "close packing"
        print(
            f"      A={A:6.1f}  r_c = {rc:6.2f}d  → a_mean ≥ {a_need:5.2f}d "
            f"(φ ≤ {phi_max:.3f})  Γ ≤ {G_at:7.3f}   <- {bind} binds"
        )
    for G_t in (10, 60):
        N_need = (8 * G_t / u_c) ** (2 / 3)
        print(f"      reaching Gamma = {G_t:3d} requires N >= {N_need:8.0f}")
print("\n    ★ N = 100 (the sketch) can only reach Gamma <~ 1. To see 'the final "
      "structure at strong coupling'")
print("      you must raise N (hundreds to thousands), accept a larger truncation "
      "error, or use r_c = L/2 plus a tail correction.")
print("      -> trap 7 (minimum-image violation, the historical +1856% error) is "
      "exactly this point.")

# ── dt: worked back from the curvature (= the local stiffness) ────────────
print(f"\nintegration resolution -- tau_int is defined from the local curvature "
      f"(the same structure as a trap's tau_k = gamma/k)")
print(f"    U'' (r) = 12 A kT d³/r⁵   →   τ_int(r) = γ/U'' = τ_B (r/d)⁵ / (12A)")
print(f"    dt <= 1e-2 tau_int(r_min) ,  r_min = the closest approach distance")
print(f"\n    {'A':>6} | {'no core: r_min(U=5kT)':>22} {'dt/tau_B':>10} {'dt [ms]':>9} |"
      f" {'core d: r_min=d':>16} {'dt/tau_B':>10} {'dt [ms]':>9}")
tau_B_ms = s["tau_B"].to("ms").magnitude
for A in (0.1, 1.0, 10.0, 100.0):
    rmin_nc = (A / 5) ** (1 / 3)
    dt_nc = 1e-2 * rmin_nc**5 / (12 * A)
    dt_c = 1e-2 * 1.0 / (12 * A)
    print(
        f"    {A:6.1f} | {rmin_nc:20.3f}d {dt_nc:10.2e} {dt_nc * tau_B_ms:9.4f} |"
        f" {1.0:15.1f}d {dt_c:10.2e} {dt_c * tau_B_ms:9.4f}"
    )
print("    ★ Without a core, **a SMALL A is numerically harder** (deeper "
      "interpenetration -> the stiffness spikes).")
print("      With a core it is the opposite: large A is harder. Which branch applies "
      "depends on a physics choice.")
print(f"    for reference: one thermal step sqrt(2 D_t dt) <= 0.01 r  ->  "
      f"dt/tau_B <= 5e-5 (r/d)^2")

# ── Cost ────────────────────────────────────────────────────────────────
print(f"\ncost (assuming T_obs = 100 tau_B = "
      f"{100 * s['tau_B'].to('s').magnitude:.0f} s)")
for A, dt_star in ((1.0, 1e-2 * 1.0 / 12), (100.0, 1e-2 / 1200)):
    steps = 100 / dt_star
    print(f"    A = {A:5.1f} (core d)  dt* = {dt_star:.2e}  ->  {steps:.2e} steps")


# ══════════════════════════════════════════════════════════════════════════
# CASE ─ trap-drag-2d-hex300
#   READ : N ~ 300, R = 5 µm, k_t = 10 pN/µm, v_x = 0.5 µm/s, hexagonal initial
#          placement
#          U_trap = ½ k_t (Δr)², Δr = |r_trap − r_i|, r_trap(t) = r_trap(0) + v t
#   OPEN : the pair potential is not stated -- there is no force on the sketch that
#          could produce a 'hexagonal equilibrium'
#   OPEN : lattice spacing / box / observation goal
# ══════════════════════════════════════════════════════════════════════════
hr("CASE  trap-drag-2d-hex300")

d_dr = Q(5.0, "um")  # READ 'R = 5µm' + the 'diameter' convention confirmed in 1-A
k_t = Q(10, "pN/um").to("N/m")
v_x = Q(0.5, "um/s").to("m/s")
t = stokes(d_dr)
tau_k = (t["gamma"] / k_t).to("s")
l_k = ((KT / k_t) ** 0.5).to("nm")
dr_ss = (t["gamma"] * v_x / k_t).to("nm")  # steady-state lag under constant-speed
                                           # dragging
snr = float((dr_ss / l_k).to("dimensionless"))
Pe_drag = float((v_x * d_dr / t["D_t"]).to("dimensionless"))
tau_v = (d_dr / v_x).to("s")

print(f"""
READ    N ~ 300   R = 5 µm   k_t = 10 pN/µm   v_x = 0.5 µm/s   hexagonal initial
        placement
        U_trap = 0.5 k_t (dr)^2 ,  r_trap(t) = r_trap(0) + v_x t
        (only ONE trap moves)
ANCHOR  d = 5 µm -- 'R=5µm' read as a diameter (the convention the user confirmed
        in 1-A, tier 1)
OPEN    ★ There is no pair potential. The sketch has no repulsion that would hold a
        hexagonal lattice together.
        (Is it the same A/r^3 system as soft-r3? Or is the lattice pinned by an
         array of traps?)
OPEN    lattice spacing a / box L / what is to be measured (drag force? lattice
        deformation? defect creation?)
""")

print("DERIVED / scale ledger (smallest first)")
print(f"    τ_p = {t['tau_p'].to('us'):~.3fP}   (model check)")
print(f"    dt  ≤ 10⁻² τ_k = {(0.01 * tau_k).to('us'):~.2fP}")
print(f"    τ_k = γ/k_t = {tau_k.to('ms'):~.3fP}   ★ the fastest physical time -> it sets dt")
print(f"    τ_B = {t['tau_B']:~.1fP}   (never realised, because the trap holds it -- same as 1-A)")
print(f"    τ_v = d/v_x = {tau_v:~.1fP}   ★ the dragging time -> it sets T_obs")
print(f"    lengths:  l_k = sqrt(kT/k) = {l_k:~.2fP}  <  d = {d_dr:~P}")
print(f"           Δr_ss = γv/k_t = {dr_ss:~.3fP}   <- the deterministic lag the dragging produces")

print("\ndimensionless groups")
print(f"    k*  = k d²/kT      = {float((k_t * d_dr**2 / KT).to('')):9.3e}   trap vs thermal fluctuation (very stiff)")
print(f"    Pe  = v d/D_t      = {Pe_drag:9.3f}     advection vs diffusion")
print(f"    τ_k/τ_v            = {float((tau_k / tau_v).to('')):9.3e}   trap relaxation vs dragging (extreme separation)")
print(f"    ★ SNR = Δr_ss/ℓ_k  = {snr:9.4f}     drag signal vs in-trap thermal fluctuation")
print(f"      = γv/√(k kT)  →  SNR ∝ v/√k")

v_snr1 = (((k_t * KT) ** 0.5) / t["gamma"]).to("um/s")
print(f"\n    ★★ Under these conditions the deterministic lag (={dr_ss:~.2fP}) is "
      f"1/10 of the thermal fluctuation (={l_k:~.2fP}).")
print(f"       the speed at which SNR = 1:  v = sqrt(k kT)/gamma = {v_snr1:~.2fP}   "
      f"(10x the sketch value)")
n_ind = (1 / (0.01 * snr)) ** 2
T_need = (2 * tau_k * n_ind).to("s")
print(f"       measuring dr_ss to 1% accuracy needs ~{n_ind:.2e} independent samples")
print(f"         -> T_obs ~ 2 tau_k x samples = {T_need:~.0fP} = "
      f"{T_need.to('hour'):~.2fP} (physical time)")
print(f"         -> with dt = 1e-2 tau_k that is "
      f"{float((T_need / (0.01 * tau_k)).to('')):.2e} steps")
print("       (for the single dragged particle. Using the response of all 300 "
      "lattice particles reduces it)")

print("\ncost -- the time to cross the box once (per assumed lattice spacing a)")
for a_over_d in (1.5, 2.0, 3.0):
    a = a_over_d * d_dr
    L = (a * math.sqrt(300 * math.sqrt(3) / 2)).to("um")  # hexagonal: area per
                                                          # particle (sqrt(3)/2)a^2
    t_cross = (L / v_x).to("s")
    steps = float((t_cross / (0.01 * tau_k)).to(""))
    print(
        f"    a = {a_over_d:.1f}d = {a:~.1fP}  →  L = {L:~.0fP}   crossing {t_cross:~.0fP}"
        f"   =  {steps:.2e} steps   (ℓ_k/L = {float((l_k / L).to('')):.1e})"
    )


# ══════════════════════════════════════════════════════════════════════════
# CASE ─ chain-bend-2d-oscill
#   READ : a bead chain, trapped by optical tweezers, k_t = 10 pN/µm, R = 5 µm,
#          y = a sin(omega t) (oscillation in y), goal G' and G'' vs omega,
#          "see Eric Furst's papers"
#   OPEN : ★ U_ij (the chain bonding potential) is left **blank** on the sketch
#   OPEN : a, the omega range, the number of beads, which beads are trapped
# ══════════════════════════════════════════════════════════════════════════
hr("CASE  chain-bend-2d-oscill")

c = stokes(Q(5.0, "um"))
tau_k_c = (c["gamma"] / k_t).to("s")
l_k_c = ((KT / k_t) ** 0.5).to("nm")
omega_1 = (1 / tau_k_c).to("1/s")

print(f"""
READ    a bead chain (~6 circles + "..."), some marked with x = "trapped by optical
        tweezers"
        U = ½ k_t (Δr)²   k_t = 10 pN/µm   R = 5 µm
        y = a sin(omega t)  "oscillation in y-dir"   goal: G' & G'' = ?
                                                     (a graph vs omega)
        "particles are connected by U_ij = ~~~"   <- **left BLANK**
        "Eric Furst, see the papers"              <- a literature pointer
                                                     (candidate KB seed)
OPEN    ★ U_ij is not stated -> the chain's elastic and relaxation spectrum cannot
        be computed (the core physics of this case)
OPEN    amplitude a, frequency range omega, bead count N, which beads are trapped
        (the x marks appear to be on beads 1, 4 and the ends)
OPEN    the bending stiffness (the folder is named 'chain-bend') is absent from the
        sketch -- no angle potential is stated
""")

print("the part that CAN be computed (trap + drive only)")
print(f"    τ_k = γ/k_t = {tau_k_c.to('ms'):~.3fP}     ℓ_k = {l_k_c:~.2fP}")
print(f"    new timescale tau_omega = 1/omega  ->  dimensionless De = omega tau_k "
      f"(Deborah)")
f_1 = (omega_1 / (2 * PI)).to("Hz").magnitude
print(f"    ★ the frequency at which De = 1:  omega = 1/tau_k = "
      f"{omega_1.magnitude:.1f} rad/s  (f = {f_1:.1f} Hz)")
print(f"      the G'/G'' crossover appears near here -> the omega sweep should cover "
      f"{omega_1.magnitude / 100:.1f} to "
      f"{omega_1.magnitude * 10:.0f} rad/s or so")
print(f"\n    the amplitude window (two dimensionless groups squeeze it from both "
      f"sides)")
print(f"      a/l_k >> 1  (to get the signal above the thermal fluctuation)   "
      f"l_k = {l_k_c:~.2fP}")
print(f"      a/d   << 1  (to stay in linear response)                        "
      f"d   = 5 µm")
for a_nm in (20, 100, 200, 500):
    a_q = Q(a_nm, "nm")
    print(
        f"      a = {a_q:~4.0fP}  →  a/ℓ_k = {float((a_q / l_k_c).to('')):6.2f} (SNR)"
        f"   a/d = {float((a_q / Q(5.0, 'um')).to('')):7.4f}"
    )
print("      -> a ~ 100-500 nm is the window satisfying both at once")

print("\n★ A point to flag in the physical reading (rule 6: verify before stating "
      "-- this is a property of the model)")
print("    In BD with a Newtonian solvent (water), **the solvent has G' = 0 and "
      "G'' = eta*omega** (no elasticity).")
print("    So the G' and G'' obtained here are NOT the 'medium' but **the elasticity "
      "of the chain (its bonds)**.")
print("    If the intent was to measure the medium's viscoelasticity, BD cannot do "
      "it and a medium.* module would be required")
print("    (the same item as the 'non-Newtonian medium' follow-up request in 1-A's "
      "observation.yaml).")

print("\ncost (dt = 1e-2 tau_k, observing 100 cycles at the lowest omega)")
for de in (0.1, 1.0):
    omega = de / tau_k_c
    period = (2 * PI / omega).to("s")
    T_obs = 100 * period
    print(
        f"    De = {de:4.1f}  ω = {omega.to('1/s'):~.1fP}  period {period.to('ms'):~.1fP}"
        f"  T_obs = {T_obs:~.2fP}  →  {float((T_obs / (0.01 * tau_k_c)).to('')):.2e} steps"
    )


# ══════════════════════════════════════════════════════════════════════════
# CASE ─ abp-rod-2d-run-flip
#   READ : ellipsoid "active" particle, run-and-flip motion
#          τ_R : rotation time 0.5 s ,  v : speed  v ≤ 5 µm/s
#          measure MSD, MSAD          the drawing marks dimensions "2R" and "R"
#                                     (no numbers)
#   OPEN : the value of R, the aspect ratio, N, the box, and (below) the flip rate
# ══════════════════════════════════════════════════════════════════════════
hr("CASE  abp-rod-2d-run-flip")

tau_R = Q(0.5, "s")
v_a = Q(5.0, "um/s")

print(f"""
READ    "ellipsoid, 'active' particle    run-and-flip motion"
        τ_R : rotation time = 0.5 s        v : speed  v ≤ 5 µm/s
        "measure MSD. MSAD."
        drawing: an ellipse + a dimension line "2R" (reading somewhat uncertain)
                 + an arrow below with "R"
OPEN    ★ R has no number (the other 3 sketches all write 'R = 5 µm')
OPEN    aspect ratio p = a/b, particle count N, box L
""")

# ── Reading tau_R as a rotational-diffusion time works the size back out
#    (sphere approximation) ──
d_impl = ((KT * tau_R / (PI * ETA)) ** (1 / 3)).to("um")
print("★ Reading tau_R = 0.5 s as a 'rotational diffusion time' works the size back "
      "out (sphere approximation)")
print(f"    τ_r = 1/D_r = πηd³/kT   →   d = (kT τ_R/πη)^(1/3) = {d_impl:~.3fP}")
print("    -> in water@300K, tau_R=0.5s matches a sphere of diameter ~0.92 µm "
      "exactly.")
print("      (an ellipsoid has larger rotational friction, so its equivalent "
      "diameter is slightly smaller)")
print("    => applying the 'R = 5 µm' convention here makes tau_R disagree with the "
      "sketch value (table below)")

print(f"\ndimensionless groups per size candidate (v = {v_a:~P}, "
      f"l_p = v*tau_R = {(v_a * tau_R).to('um'):~.2fP}, which is fixed)")
print(f"    {'d':>8} {'D_t':>10} {'τ_B':>9} {'τ_v=d/v':>9} {'Stokes τ_r':>11} "
      f"{'Pe=vd/D_t':>10} {'ℓ_p/d':>7} {'D_r*':>6}")
for d_um in (0.5, d_impl.magnitude, 1.0, 5.0, 10.0):
    a = stokes(Q(d_um, "um"))
    Pe = float((v_a * a["d"] / a["D_t"]).to(""))
    l_p = (v_a * tau_R).to("um")
    print(
        f"    {d_um:6.2f}µm {a['D_t'].magnitude:8.4f}µm²/s {a['tau_B'].to('s').magnitude:8.2f}s"
        f" {(a['d'] / v_a).to('s').magnitude:8.3f}s {a['tau_r'].to('s').magnitude:10.3f}s"
        f" {Pe:10.2f} {float((l_p / a['d']).to('')):7.2f} {float((a['D_r'] * a['tau_B']).to('')):6.2f}"
    )
print("    -> taking d = 5 µm gives l_p/d = 0.5: it reorients before travelling even "
      "half a body length (activity is barely visible)")
print("      at d ~ 1 µm, l_p/d ~ 2.5 and Pe ~ 9 -- a clear active regime appears in "
      "the MSD")
print("      ★ the reading self-consistent with the sketch's tau_R is d ~ 1 µm")

print("\n★ Mapping observables to parameters -- one more parameter is missing here")
print("    An ellipsoid's 'body orientation' n and its 'propulsion direction' +-n "
      "are distinct.")
print("      . MSAD  <- rotational diffusion of the body (tau_R)   ... on the sketch")
print("      . MSD   <- propulsion speed v + **the flip rate tau_flip**"
      "   ... tau_flip is absent")
print("    Either reading leaves one parameter missing:")
print("      (a) tau_R = the rotational diffusion time -> tau_flip is not stated "
      "[consistent with the size inversion -> leaning this way]")
print("      (b) tau_R = the flip interval -> rotational diffusion is not stated, "
      "which makes the MSAD trivial")
print("    velocity correlation for run-and-flip: <v(0).v(t)> = "
      "v^2 exp(-2t/tau_flip)  (Poisson +- reversal)")
print("      -> the effective persistence time is tau_flip/2, and that sets the MSD "
      "crossover")
print("    NOTE this narrows the ambiguity of mater_plan §20 question 10 to this "
      "form (resolving it needs a human)")

print("\nhard constraint (measured, bd-hoomd): BD has no HI, so **the translational "
      "friction is isotropic**")
print("    -> MSAD ✓ / long-time MSD ✓ / **anisotropy of the short-time MSD ✗**")
print("    the sketch says only 'MSD, MSAD' and does not call for anisotropy")
print("    -> can proceed with mater_plan §20 option A (isotropic mean friction + an "
      "exact gamma_r). Needs confirmation")

hr("summary -- no case has a 'complete' non-dimensionalization. Where each is stuck",
   "-")
print("""
  soft-r3        no density (phi/L) + unknown whether a core exists -> the A sweep
                 is not defined.
                 On top of that, N=100 can only reach Gamma <~ 1 (minimum-image
                 constraint)
  trap-drag      no pair potential (no force that would build the hexagonal lattice)
                 + no lattice spacing.
                 At v=0.5µm/s the signal is 1/10 of the thermal fluctuation (a
                 measurement-design problem)
  chain-oscill   U_ij is blank + no a, no omega -> the chain relaxation spectrum
                 cannot be computed
  abp-rod        no value for R + one of tau_flip / rotational diffusion is missing
""")
