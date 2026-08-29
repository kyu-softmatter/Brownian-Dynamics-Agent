"""`chain-bend-2d-oscill` — L3 non-dimensionalization + L4 builder.

Three-point optical-tweezer bending: a colloidal bead chain is oscillated along y to
measure G'(ω) and G''(ω).
The end beads sit in fixed traps (force sensors); the centre bead's trap is driven as
`y = a sin(ωt)`.

The physics comes from two papers the user specified (the sketch's `U_ij` was a squiggle
= blank):
  [P1] Pantina & Furst, PRL 94, 138301 (2005)
  [P2] Pantina & Furst, Langmuir 24, 1141-1146 (2008)
★ The bond is **not a pair potential** — it is an adhesive contact (JKR) plus a
tangential bending stiffness κ_θ = EI/ℓ.
The extraction is verified in `scratch/chain_bend_from_papers.py` (reproduces the
papers' κ₀ to +1.6%; the discrete<->continuum mapping is
−0.35% at n=25 — the −0.08% that system.yaml cites is the n=51 row of that
convergence table).

**What is fundamentally different from the preceding three cases**:

  ① What sets `dt` is **not what is being observed.** The fastest bending mode, coming
     from the largest eigenvalue of the stiffness matrix (τ_fast = 0.279 µs), sets dt —
     but what we want is the collective bending at τ_chain (1.27 ms), a factor 4570
     away. In 1-A, 1-B and trap-drag the scale that set dt was the scale of interest,
     or close to it.
  ② ★ The fastest mode is not overdamped (τ_p/τ_fast = 0.60, ζ=0.65). BD treats that
     mode as overdamped, so the dynamics in that band is wrong and **no choice of dt
     fixes it**.
     ✔ **Resolved** — comparing `OverdampedViscous` against `Langevin(kT=0)` at all
     seven ω, the difference in K*(ω) is **at most 0.159%**. No effect on the observed
     band, measured
     (`scratch/verify_chain_bend_gates.py --gate det --collect`).
  ③ There are no periodic boundaries -> the `box` role is stated explicitly via
     `declare_absent`.
  ④ The amplitude is squeezed from **three sides**: ℓ_k << a < δ_max **and**
     min|θ−π| > SMALL.
     The upper limit is M_c (beyond it the bond slips), the lower one is SNR, and ⑤ is
     the new constraint.

⚠️ **⑤ This case cannot be run with `angle.Harmonic`** (2026-08-05).
   HOOMD clamps `sin θ` at `ANGLE_SIN_SMALL` (=√2×10⁻³), so for a nearly straight chain
   **only the force** is shrunk (the energy stays exact -> an energy check cannot catch
   it). In this system **all 23 angles** of the response profile are in that regime
   (min|θ−π| = 7.3e-5 = SMALL/19).
   ★ **This does NOT block the case.** The Korean original said `--spec` writes nothing
   and `--run` is refused; that was true only while the bending force was
   `angle.Harmonic`. `BENDING_IMPL` is now `custom_linear` (a `force.Custom` that
   computes F = −A y directly, see `make_bending_force`), the hard check is conditional
   on the implementation, and measured on 2026-08-29 this label carries 15 specs with 8
   of their run directories present. Corrected rather than translated forward.
   The L4 construction itself (ghost traps, drive, lock-in) is verified by two gates and
   **the builder is complete**.
   Details and the way around it: `assert_angle_force_valid()`, skill bd-hoomd trap 15.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/chain_bend_2d.py --report              # ω = ω_min (the most expensive point)
    $PY cases/chain_bend_2d.py --omega 7853 --report
    $PY cases/chain_bend_2d.py --sweep --spec        # the whole sweep as specs
    $PY cases/chain_bend_2d.py --run                 # L4 run (see ⑤ — allowed under custom_linear)
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bdbot import checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import lockin as LI, run as RUN, sim as SIM  # noqa: E402
from bdbot import nondim as ND, scales as SC  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from bdbot.units import Q  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
N_CYCLES = 100.0            # observation window = how many cycles at the lowest ω
N_SWEEP = 7                 # number of ω sweep points (log-spaced)
AC_GATE = 1e-1              # ★proposed: may the contact be treated as a point hinge — upper bound on a_c/d
THETA_GATE = 1e-1           # ★proposed: small-angle linearization limit of the harmonic angle [rad]

# ★★ A hard HOOMD constant (established by measurement,
#    `scratch/verify_angle_force_small_theta.py`).
# md.angle.Harmonic uses 1/sin θ when moving torque into coordinates, and **clamps**
# sin θ at this value.
# For sin θ < SMALL the force is shrunk by sinθ/SMALL -> force ∝ κ(θ−π)²/SMALL, i.e.
# **quadratic** rather than linear, so
# the chain comes out far softer and more nonlinear than it really is. **The energy is
# exact** (0.000%) — which is why
# verifying with energy passes and leaves the force wrong. With t0=π the equilibrium
# itself has sin θ=0, so
# a stiff chain is always in this regime. Measured SMALL = 1.414217e-03
# (standard deviation 1.4e-7) = √2×10⁻³.
ANGLE_SIN_SMALL = 1.414214e-3


# ════════════════════════════════════════════════════════════════════════
# ① The physical system (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    P = load_node
    con = raw["interactions"][0]
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "n": int(raw["geometry"]["n_beads"]),
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "E": P(raw["particle"]["youngs_modulus"]),
        "kappa_0": P(con["kappa_0"]),
        "a_c": P(con["contact_radius"]),
        "k_bond": P(con["bond_stiffness"]),
        "M_c": P(con["critical_moment"]),
        "k_t": P(raw["external"]["stiffness"]),
        "amp": P(raw["external"]["amplitude"]),
        "omega_range": [float(x) for x in raw["external"]["omega_range"]["value"]],
        "n_trapped": int(raw["external"]["n_trapped"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② Numerics specific to this system — the stiffness matrix of the discrete chain
#    ★ Not written from intuition. Same construction as
#      scratch/chain_bend_from_papers.py, which reproduces the papers' values; here the
#      **eigenvalues** are used as well (there, only the stiffness was).
# ════════════════════════════════════════════════════════════════════════
def bending_matrix(n: int, kappa_theta: float, ell: float) -> np.ndarray:
    """Quadratic form A (along y) of U = ½ κ_θ Σ θ_i², θ_i = (y_{i+1} − 2y_i + y_{i−1})/ℓ.

    Small-deformation approximation. `A = κ_θ Bᵀ B`, where B is the (n−2)×n second-
    difference matrix.
    """
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def bond_matrix(n: int, k_bond: float) -> np.ndarray:
    """Quadratic form of bond stretching (along x) — the path-graph Laplacian times k_b.

    ★ For a straight chain, stretching (x) and bending (y) **decouple at linear order**.
      So the largest eigenvalue is taken separately from the two blocks and the larger
      one is used. Adding them into a single matrix mixes decoupled degrees of freedom,
      overestimates λ_max, and makes dt smaller than it needs to be.
    """
    G = np.zeros((n, n))
    for i in range(n - 1):
        G[i, i] += 1.0
        G[i + 1, i + 1] += 1.0
        G[i, i + 1] -= 1.0
        G[i + 1, i] -= 1.0
    return k_bond * G


def trapped_indices(n: int) -> list[int]:
    """[P2] Fig.1A — the two ends (fixed, force sensors) plus the centre (driven): 3."""
    return [0, n // 2, n - 1]


def three_point_bending(n: int, kappa_theta: float, ell: float, delta: float):
    """Ends fixed, centre pushed by delta: returns (centre stiffness F/δ, max bond
    angle |θ|).

    Same construction as
    `scratch/chain_bend_from_papers.discrete_3point_stiffness`
    (already compared against the beam formula in that convergence table —
    **n=25 gives −0.35%**, n=51 gives −0.08%),
    but here the angle is returned as well — the angle is needed to see whether the
    small-angle linearization of the harmonic angle potential holds.
    """
    A = bending_matrix(n, kappa_theta, ell)
    fixed = [0, n // 2, n - 1]
    free = [i for i in range(n) if i not in fixed]
    y_fix = np.array([0.0, delta, 0.0])
    y = np.zeros(n)
    y[fixed] = y_fix
    y[free] = np.linalg.solve(A[np.ix_(free, free)], -A[np.ix_(free, fixed)] @ y_fix)
    U = 0.5 * y @ A @ y
    k_center = 2 * U / delta**2                      # U = ½ k δ²
    theta = np.abs(np.diff(y, n=2)) / ell            # |θ_i|
    return k_center, float(theta.max())


def trapped_stiffness_matrix(n: int, kappa_theta: float, ell: float, k_t: float):
    """Stiffness matrix of bending plus the three traps, and the driven bead index.

    ★ The ends are **traps, not rigid clamps** (finite stiffness k_t). The difference is
      more than 30%, so a value derived with rigid boundary conditions must not be
      called "what the trap feels".
    """
    A = bending_matrix(n, kappa_theta, ell)
    idx = trapped_indices(n)
    for i in idx:
        A[i, i] += k_t
    return A, idx[1]


def driven_static_stiffness(n: int, kappa_theta: float, ell: float, k_t: float) -> float:
    """The static stiffness K(ω->0) the driving trap **actually** feels.

    Solve (A_bend + T) y = k_t y_c e_mid and take K = k_t(y_c/y_mid − 1). Because the
    end traps have finite stiffness the whole chain shifts, so this is **smaller** than
    the 48EI/L³ of the rigid-clamp assumption.
    Measured comparison: at the lowest ω (De_true≈1) HOOMD gives 0.95x this value
    (`scratch/verify_chain_bend_gates.py --gate det`).
    """
    A, mid = trapped_stiffness_matrix(n, kappa_theta, ell, k_t)
    e = np.zeros(n)
    e[mid] = k_t
    return float(k_t * (1.0 / np.linalg.solve(A, e)[mid] - 1.0))


def driven_response(n: int, kappa_theta: float, ell: float, k_t: float,
                    gamma: float, omega: float, amp: float) -> float:
    """Response amplitude |y_hat(ω)| of the driven bead — small-angle linear response
    (iωγI + A + T) y_hat = k_t a e_mid.

    ★ Used to measure SNR. The a/ℓ_k the spec used to check is the **drive amplitude**,
      whereas the observable is the **response** — the chain and the drag resist, so
      |y_hat| << a, and at high frequency it drops below ℓ_k.

    ✔ **Cause identified (2026-08-05).** The old 28% disagreement with the HOOMD
      measurement was not this prediction being wrong but **HOOMD being wrong** —
      `md.angle.Harmonic` clamps sin θ at `ANGLE_SIN_SMALL` and so shrinks the force for
      a nearly straight chain. All 23 angles of this system are in that regime.
      This model agrees with an exact nonlinear minimization (scipy, all 2n coordinates)
      to **0.32%**, so it is the correct one.
      Trail: `scratch/diagnose_chain_bend_28pct.py`,
      `scratch/verify_angle_force_small_theta.py`
    """
    A, mid = trapped_stiffness_matrix(n, kappa_theta, ell, k_t)
    e = np.zeros(n, dtype=complex)
    e[mid] = k_t * amp
    M = 1j * omega * gamma * np.eye(n) + A
    y = np.linalg.solve(M, e)
    return float(abs(y[mid])), np.abs(y)


def response_angles(prof: np.ndarray, ell: float) -> tuple[float, float]:
    """(max, min) of the bond-angle deviation |θ−π| over the response amplitude profile.

    ★★ Used for the `ANGLE_SIN_SMALL` verdict. For the force to be right **every** angle
    must be above it — looking only at the maximum is not enough. The angles at the
    chain ends are more than an order of magnitude smaller than at the centre.
    """
    th = np.abs(np.diff(prof, n=2)) / ell
    return float(th.max()), float(th.min())


# ════════════════════════════════════════════════════════════════════════
# ③ The scale ledger (bd-physics §0 ①②)
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, omega: float, *, dt_scale=1.0, n_cycles=N_CYCLES) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B = b["kT"], b["gamma"], b["tau_B"]
    k_t = sys_["k_t"].value.to("N/m")
    amp = sys_["amp"].value.to("m")
    n = sys_["n"]
    a_rad = d / 2                                    # particle radius (the papers' a)

    # ── Contact -> bending stiffness ([P2] eq.5 + the discrete mapping) ──────
    EI = (sys_["kappa_0"].value.to("N/m") * a_rad**3 / 3).to("N*m^2")
    kappa_theta = (EI / d).to("J")                   # κ_θ = EI/ℓ,  ℓ = d
    L_chain = (n - 1) * d                            # contour length
    # ★ Mixing the two force definitions disagrees with the papers by exactly a factor
    #   of 2 (chain_bend_from_papers §②).
    kappa_end = (24 * EI / L_chain**3).to("N/m")     # the papers' definition (end-force based)
    # ★ 48EI/L³ is the three-point bending value for **rigidly clamped ends**. This
    #   system's ends are traps (k_t), so it is not what the driving trap feels — that
    #   is kappa_drive below (a 32% difference).
    kappa_center = 2 * kappa_end                     # beam formula 48EI/L^3 (rigid-clamp assumption)
    M_c = sys_["M_c"].value.to("N*m")
    delta_max = (M_c * L_chain**2 / (12 * EI)).to("m")    # linear-elastic limit amplitude, M<M_c

    # ── Timescales ────────────────────────────────────────────────────────
    tau_k = C.relaxation_time(gamma, k_t)                       # trap
    tau_chain = C.relaxation_time(gamma, kappa_center)          # gamma/kappa_center (beam-formula based)
    k_b = sys_["k_bond"].value.to("N/m")

    # ★ Fastest mode — largest eigenvalue of the stiffness matrix. Stretching (x) and
    #   bending (y) decouple, so they are taken separately.
    idx = trapped_indices(n)
    kth = float(kappa_theta.magnitude)
    ell = float(d.magnitude)
    Ay = bending_matrix(n, kth, ell)
    Ax = bond_matrix(n, float(k_b.magnitude))
    for i in idx:                                    # a trap adds k_t on the diagonal
        Ay[i, i] += float(k_t.magnitude)
        Ax[i, i] += float(k_t.magnitude)
    ev_bend = np.linalg.eigvalsh(Ay)
    lam_bend = float(ev_bend[-1])
    lam_bond = float(np.linalg.eigvalsh(Ax)[-1])
    lam_max = max(lam_bend, lam_bond)
    tau_fast = C.relaxation_time(gamma, Q(lam_max, "N/m"))
    tau_bond = C.relaxation_time(gamma, Q(lam_bond, "N/m"))

    # ★★ Smallest eigenvalue -> the **longest** relaxation time. This used not to be
    #    computed, and τ_chain = γ/κ_center was used as the governing scale — but τ_max
    #    is 9.18x longer. As a result (a) De was underestimated by 9.18x so the ω sweep
    #    never entered the quasi-static region at all, and (b) equilibration
    #    (20 τ_chain = 2.2 τ_max) was insufficient so K* was wrong by up to 21%.
    #    Confirmed by measurement:
    #    `scratch/verify_chain_bend_gates.py --gate det --eq-steps` (raising
    #    equilibration to 10 τ_max shrinks the block scatter 1000x and moves K' by 21%).
    lam_min = float(ev_bend[0])
    tau_max = C.relaxation_time(gamma, Q(lam_min, "N/m"))
    # Does the low-frequency end of the **whole** sweep reach quasi-static? Judged on
    # the lower bound of the range, not this spec's ω (looking at De per point lets the
    # lowest point pass and hides the range problem).
    de_lo = float((Q(min(sys_["omega_range"]), "1/s") * tau_max).to("dimensionless").magnitude)
    kappa_drive = driven_static_stiffness(n, kth, ell, float(k_t.magnitude))
    y_resp, y_prof = driven_response(n, kth, ell, float(k_t.magnitude),
                                     float(gamma.to("kg/s").magnitude), omega,
                                     float(amp.magnitude))
    th_hi, th_lo = response_angles(y_prof, ell)      # ★★ for the HOOMD angle-force validity verdict

    dt = dt_scale * C.dt_from_gate(tau_fast)         # the fastest mode sets dt
    tau_w = Q(1.0 / omega, "s")                      # drive (De = tau_chain/tau_w)
    tau_period = Q(2 * math.pi / omega, "s")
    T_obs = (n_cycles * tau_period).to("s")

    l_k = (kT / k_t) ** 0.5
    k_center_disc, theta_max = three_point_bending(n, kth, ell, float(amp.magnitude))

    lg = SC.ScaleLedger()
    lg.add_length("l_k", l_k.to("m"), "sqrt(kT/k_t) trap fluctuation (noise floor)", star=True)
    lg.add_length("y_resp", Q(y_resp, "m"),
                  "★ response amplitude |y_hat(omega)| of the driven bead — the "
                  "observable. This, not a, is the numerator of SNR"
                  " (linear-response estimate; the old 28% gap against HOOMD at high "
                  "frequency is explained — see driven_response)", star=True)
    lg.add_length("a_c", sys_["a_c"].value.to("m"), "JKR contact radius (point-hinge assumption)")
    lg.add_length("a", amp, "drive amplitude", star=True)
    lg.add_length("delta_max", delta_max, "M_c linear-elastic limit amplitude", star=True)
    lg.add_length("d", d, "bead diameter = bond length l (reference)")
    lg.add_length("L_chain", L_chain.to("m"), "chain contour length (n-1)d")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_fast", tau_fast, "gamma/lambda_max fastest bending mode — sets dt", star=True)
    lg.add_time("tau_bond", tau_bond, "gamma/lambda_max(stretch) bond stretching")
    lg.add_time("tau_w", tau_w, f"1/omega drive (omega = {omega:.0f} rad/s)")
    lg.add_time("tau_k", tau_k, "gamma/k_t trap")
    lg.add_time("tau_chain", tau_chain, "gamma/kappa_center (beam-formula based — NOT the governing scale)")
    lg.add_time("tau_max", tau_max,
                "★★ gamma/lambda_min longest relaxation time — **the governing scale**. "
                "De and equilibration must use this",
                star=True)
    lg.add_time("tau_period", tau_period, "2*pi/omega drive period")
    lg.add_time("tau_B", tau_B, "d^2/D_t diffusion (reference)")
    lg.add_time("T_obs", T_obs, f"observation window ({n_cycles:g} cycles)", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("k_t_d2", (k_t * d**2).to("J"), "k_t*d^2 trap stiffness")
    lg.add_energy("kappa_end_d2", (kappa_end * d**2).to("J"), "kappa_end d^2 chain stiffness (the papers' definition)")
    lg.add_energy("kappa_drive_d2", Q(kappa_drive, "N/m").to("N/m") * d**2,
                  "★ kappa_drive d^2 — the static stiffness the driving trap "
                  "**actually** feels "
                  "(the end traps have finite stiffness; 0.68x the rigid-clamp "
                  "48EI/L^3)", star=True)
    lg.add_energy("k_b_d2", (k_b * d**2).to("J"), "k_b d^2 bond stretching stiffness")
    lg.add_energy("M_c", M_c, "critical bending moment (moment has energy dimensions)")
    lg.add_energy("kappa_theta", kappa_theta, "bond-angle stiffness kappa_theta = EI/l", star=True)
    # ★ There are no periodic boundaries — a single chain, with traps fixing its
    #   position. Left empty rather than invented.
    lg.declare_absent(
        "box",
        "No periodic boundaries (geometry.periodic=false). A single chain whose position "
        "is fixed by traps, so "
        "there is no box in the physics to act as the denominator of a minimum-image or "
        "finite-size check. L4 only has to make the HOOMD box "
        "larger than the chain's maximum extent plus the WCA cutoff, and that value does "
        "not change the physics. "
        "What plays the role of the geometric limit here is not the box but "
        "delta_max (M<M_c).")

    lg.derived = dict(gamma=gamma, D_t=b["D_t"], m=b["m"], kT=kT, d=d, tau_B=tau_B,
                      EI=EI, kappa_theta=kappa_theta, kappa_end=kappa_end,
                      kappa_center=kappa_center, delta_max=delta_max, M_c=M_c,
                      L_chain=L_chain.to("m"), l_k=l_k.to("m"), k_t=k_t, k_b=k_b,
                      tau_k=tau_k, tau_chain=tau_chain, tau_fast=tau_fast,
                      tau_max=tau_max, lam_min=lam_min, kappa_drive=kappa_drive,
                      y_resp=y_resp, de_lo=de_lo, th_hi=th_hi, th_lo=th_lo,
                      tau_bond=tau_bond, tau_w=tau_w, tau_period=tau_period,
                      dt=dt, T_obs=T_obs, omega=omega, n=n, amp=amp,
                      lam_bend=lam_bend, lam_bond=lam_bond, lam_max=lam_max,
                      k_center_disc=k_center_disc, theta_max=theta_max,
                      trapped=trapped_indices(n))
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ This system's timescales span four decades, and the "
        "scale that sets dt (tau_fast, the fastest bending mode) is **not what is being "
        "observed** — what we want is tau_chain. "
        "So the cost is fixed independently of the observed band.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ④ Dimensionless groups + ⑤ separation checks (bd-physics §3, §4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg):
    n = sys_["n"]
    D = lg.derived
    r = lg.ratio
    n = sys_["n"]

    groups = [
        ND.Group("kappa_theta/kT", r("energies", "kappa_theta", "kT"),
                 ("energies", "kappa_theta"), ("energies", "kT"), "EI/(ℓ kT)",
                 "★ the chain is thermally completely stiff (fluctuations cannot bend it)"),
        ND.Group("M_c/kT", r("energies", "M_c", "kT"),
                 ("energies", "M_c"), ("energies", "kT"), "",
                 "the moment at which the bond slips (validity range of the harmonic angle)"),
        ND.Group("k*", r("energies", "k_t_d2", "kT"),
                 ("energies", "k_t_d2"), ("energies", "kT"), "k_t d²/kT",
                 "trap vs thermal fluctuation"),
        ND.Group("kappa_end/k_t", r("energies", "kappa_end_d2", "k_t_d2"),
                 ("energies", "kappa_end_d2"), ("energies", "k_t_d2"), "",
                 "★ chain vs trap stiffness — they must be comparable to fall in the measurement window"),
        ND.Group("k_b/k_t", r("energies", "k_b_d2", "k_t_d2"),
                 ("energies", "k_b_d2"), ("energies", "k_t_d2"), "",
                 "bond stretching is effectively inextensible (192x stiffer than the trap)"),
        ND.Group("a/l_k", r("lengths", "a", "l_k"), ("lengths", "a"), ("lengths", "l_k"),
                 "", "★ SNR — amplitude vs in-trap thermal fluctuation"),
        ND.Group("a/delta_max", r("lengths", "a", "delta_max"),
                 ("lengths", "a"), ("lengths", "delta_max"), "",
                 "★ linear-elastic margin (M<M_c). Above 1 the harmonic angle is invalid"),
        ND.Group("a/d", r("lengths", "a", "d"), ("lengths", "a"), ("lengths", "d"), "",
                 "amplitude vs bead"),
        ND.Group("a_c/d", r("lengths", "a_c", "d"), ("lengths", "a_c"), ("lengths", "d"), "",
                 "is the contact a point hinge — the premise of the JKR derivation"),
        ND.Group("L_chain/d", r("lengths", "L_chain", "d"),
                 ("lengths", "L_chain"), ("lengths", "d"), "n−1", "chain length in beads"),
        ND.Group("n_beads", float(n), None, None, "", "bead count (input, ★proposed)"),
        ND.Group("De", r("times", "tau_max", "tau_w"),
                 ("times", "tau_max"), ("times", "tau_w"), "ω τ_max",
                 "★★ Deborah — **based on the longest relaxation time**. Quasi-static limit is De << 1"),
        ND.Group("De_chain_old", r("times", "tau_chain", "tau_w"),
                 ("times", "tau_chain"), ("times", "tau_w"), "ω τ_chain",
                 "★ the old definition (based on gamma/kappa_center). Underestimates De by 9.18x — kept only for the record"),
        ND.Group("tau_max/tau_chain", r("times", "tau_max", "tau_chain"),
                 ("times", "tau_max"), ("times", "tau_chain"), "κ_center/λ_min",
                 "★★ the governing scale is this much longer than the beam formula — the factor the old spec missed"),
        ND.Group("kappa_drive/k_t", r("energies", "kappa_drive_d2", "k_t_d2"),
                 ("energies", "kappa_drive_d2"), ("energies", "k_t_d2"), "",
                 "★ chain vs trap — the value that accounts for the trap boundary (this, not kappa_end/k_t)"),
        ND.Group("y_resp/l_k", r("lengths", "y_resp", "l_k"),
                 ("lengths", "y_resp"), ("lengths", "l_k"), "",
                 "★★ the real SNR — is the response larger than the thermal fluctuation? a/l_k overestimates this by up to 60x"),
        ND.Group("De_trap", r("times", "tau_k", "tau_w"),
                 ("times", "tau_k"), ("times", "tau_w"), "ω τ_k", "trap reference"),
        ND.Group("tau_fast/tau_chain", r("times", "tau_fast", "tau_chain"),
                 ("times", "tau_fast"), ("times", "tau_chain"), "",
                 "★ scale-separation width — the mode that sets dt is this much faster than the mode of interest"),
        ND.Group("dt/tau_fast", r("times", "dt", "tau_fast"),
                 ("times", "dt"), ("times", "tau_fast"), "", "integration resolution (governing scale)"),
        ND.Group("n_cycles", r("times", "T_obs", "tau_period"),
                 ("times", "T_obs"), ("times", "tau_period"), "", "cycles observed"),
        ND.Group("St", r("times", "tau_p", "tau_B"),
                 ("times", "tau_p"), ("times", "tau_B"), "tau_p/tau_B", "inertia vs diffusion"),
    ]
    checks = [
        C.Check("model", "inertia negligible   tau_p/tau_max", r("times", "tau_p", "tau_max"),
                C.GATE, "<=",
                "tau_dyn = the governing scale of the observed band = tau_max (longest "
                "relaxation). That is the band in which G'(omega) is measured"),
        C.Check("model", "note: tau_p/tau_fast", r("times", "tau_p", "tau_fast"), C.GATE, "<=",
                "★ The fastest bending mode is not overdamped "
                "(zeta = gamma/2*sqrt(m*lambda_max) = 0.65 < 1). BD treats that "
                "mode as overdamped, so the dynamics in that band is wrong and **no "
                "choice of dt fixes it**. ✔ **Measured** — with the same parameters, "
                "OverdampedViscous vs "
                "Langevin(kT=0) compared at all seven omega gives a K*(omega) "
                "difference of at most **0.159%** "
                "(De_old=10). No effect on the observed band, confirmed "
                "(`scratch/verify_chain_bend_gates.py --gate det --collect`). "
                "thermal ringing combined with nonlinear coupling is NOT covered by "
                "this test", hard=False),
        C.Check("model", "linear elasticity    a/delta_max", r("lengths", "a", "delta_max"), 1.0, "<=",
                f"★ M < M_c. Beyond it the bond slips or rolls ([P2]'s conclusion) and "
                f"the harmonic angle "
                f"potential becomes invalid. delta_max = M_c L^2/(12EI) = "
                f"{D['delta_max'].to('nm'):~.0fP}"),
        C.Check("model", "small-angle linear   max|theta| [rad]", D["theta_max"], THETA_GATE, "<=",
                f"The harmonic angle is a small-angle approximation. This is the maximum "
                f"bond angle from solving the discrete three-point bending directly at "
                f"a={D['amp'].to('nm'):~.0fP}. "
                f"The upper bound {THETA_GATE:g} rad is ★proposed", hard=False),
        C.Check("model", "point contact        a_c/d", r("lengths", "a_c", "d"), AC_GATE, "<=",
                f"kappa_0 = 3*pi*a_c^4*E/(4a^3) presupposes a_c << a. The upper bound "
                f"{AC_GATE:g} is ★proposed",
                hard=False),
        C.Check("model", "*angle force valid   min|theta-pi|", D["th_lo"], ANGLE_SIN_SMALL, ">=",
                f"★★ **A hard HOOMD constraint.** md.angle.Harmonic clamps sin theta at "
                f"{ANGLE_SIN_SMALL:.3e}, and below that the force is shrunk by "
                f"sin(theta)/SMALL "
                f"(force ∝ kappa*(theta−pi)^2 — quadratic, not linear). "
                f"**The energy is exact** (0.000%) "
                f"-> verifying with energy passes while only the force is wrong. "
                f"**Every** angle of the response profile "
                f"must lie above it. Measured and reproducible: "
                f"scratch/verify_angle_force_small_theta.py. "
                f"In this system max|theta−pi|={D['th_hi']:.2e} but the minimum is "
                f"{D['th_lo']:.2e}, so "
                f"all {n-2} angles are in the broken regime. "
                + ("-> cannot be run with angle.Harmonic (hard)"
                   if BENDING_IMPL == "angle_harmonic" else
                   f"-> ★ bending is currently implemented as `{BENDING_IMPL}` "
                   f"(a force.Custom computing F = −A y "
                   "directly), so **this constraint does not apply.** The check is not "
                   "deleted but kept for reference — reverting the implementation to "
                   "angle.Harmonic makes it hard again immediately (BENDING_IMPL)."),
                # ★ This is a constraint of **the implementation**, not of the system.
                #   So it is not deleted but made conditional on the implementation. The
                #   validity of the custom force is judged by 'small-angle linear' above.
                hard=(BENDING_IMPL == "angle_harmonic")),
        C.Check("integration", "fastest mode resolved dt/tau_fast", r("times", "dt", "tau_fast"),
                C.GATE, "<=",
                f"largest eigenvalue of the stiffness matrix "
                f"lambda_max = {D['lam_max']:.4e} N/m (bending "
                f"{D['lam_bend']:.3e} vs stretching {D['lam_bond']:.3e} — the larger). "
                "Miss this and it blows up"),
        C.Check("integration", "stretch resolved     dt/tau_bond", r("times", "dt", "tau_bond"),
                C.GATE, "<=", "bond stretching mode (slower than bending, so there is plenty of margin)"),
        C.Check("integration", "drive resolved       dt/tau_w", r("times", "dt", "tau_w"), C.GATE, "<=",
                f"drive omega = {D['omega']:.0f} rad/s resolved"),
        C.Check("statistics", "SNR   |y_hat(w)|/l_k", r("lengths", "y_resp", "l_k"), 3.0, ">=",
                "★★ The **response** amplitude must exceed the thermal fluctuation for "
                "phase extraction to work. The old check put "
                "the drive amplitude a in the numerator, returned an omega-independent "
                "9.83 and passed — but the real SNR "
                "falls with omega and is below 1 at high frequency (measured 0.165 at "
                "De_old=10). "
                "It was the same hole as trap-drag's 'the hard checks pass but the "
                "statistics do not'. "
                "⚠ |y_hat| is a linear-response estimate; the old 28% gap against HOOMD "
                "at high frequency is now explained (see driven_response)",
                hard=False),
        C.Check("statistics", "drive amplitude      a/l_k", r("lengths", "a", "l_k"), 3.0, ">=",
                "is the drive larger than the thermal fluctuation — necessary but "
                "**not sufficient** "
                "(the |y_hat|/l_k check above is the real verdict)", hard=False),
        C.Check("statistics", "quasi-static reached De(w_min)", D["de_lo"], 0.1, "<=",
                "★★ De = omega*tau_max. The sweep must include the quasi-static limit "
                "(De << 1) for the elastic plateau of K' "
                "to be visible. Under the old definition (tau_chain) it looked like "
                "0.1-10 was covered, but in reality "
                "De ≈ 1-92, so it **never enters** the plateau region — the saturation "
                "curve the sketch asked for cannot be produced by this sweep. The omega "
                "range is system.yaml tier 3 "
                "(user-approved), so it is not changed here but exposed by the check",
                hard=False),
        C.Check("statistics", "cycles observed      T_obs/(2pi/w)", r("times", "T_obs", "tau_period"),
                N_CYCLES, ">=", "cycles used for the phase average", hard=False),
    ]
    return groups, checks


def report_blocks(sys_, lg, n_eq, n_prod):
    D = lg.derived
    inp = [R.kv(k, f"{sys_[k].value:~.4gP}", sys_[k].tier, sys_[k].source[:44])
           for k in ("d", "T", "eta", "kappa_0", "M_c", "k_t", "amp", "k_bond")]
    inp += [R.kv("n", f"{sys_['n']}", 3, "★proposed [P1] Fig.4 — at n=11 the amplitude window closes"),
            R.kv("omega", f"{D['omega']:.0f} rad/s", 3, "★proposed — one point of the De 0.1-10 sweep")]
    der = [
        f"  EI = {D['EI']:~.4eP}   κ_θ = EI/ℓ = {D['kappa_theta']:~.4eP}"
        f" = {lg.ratio('energies', 'kappa_theta', 'kT'):.3e} kT",
        f"  κ_end = {D['kappa_end'].to('pN/um'):~.3fP} (papers' definition, end force)"
        f"   κ_center = {D['kappa_center'].to('pN/um'):~.3fP} (beam 48EI/L^3, rigid-clamp assumption)",
        f"  ★ κ_drive = {Q(D['kappa_drive'], 'N/m').to('pN/um'):~.3fP}"
        f" = κ_center × {D['kappa_drive']/float(D['kappa_center'].magnitude):.3f}"
        f"  <- what the driving trap **actually** feels (the ends are traps, so finite"
        f" stiffness)",
        f"  ★★ λ_min = {D['lam_min']:.4e} N/m → τ_max = {D['tau_max'].to_compact():~.4gP}"
        f" = {float(D['tau_max']/D['tau_chain']):.2f} × τ_chain  <- **the governing scale**",
        f"  ★ Mixing the force definitions disagrees with the papers by exactly 2x"
        f" (confirmed in scratch/chain_bend_from_papers.py §②)",
        f"  discrete 3-point bending, solved directly ="
        f" {Q(D['k_center_disc'], 'N/m').to('pN/um'):~.3fP}"
        f"  vs  beam 48EI/L^3 = {D['kappa_center'].to('pN/um'):~.3fP}"
        f"  ({100*(D['k_center_disc']/float(D['kappa_center'].magnitude)-1):+.2f}%)",
        f"  lambda_max = {D['lam_max']:.4e} N/m  (bending {D['lam_bend']:.3e}, "
        f"stretching {D['lam_bond']:.3e} — the two blocks decouple for a straight "
        f"chain)",
        f"  δ_max = M_c L²/(12EI) = {D['delta_max'].to('nm'):~.0fP}"
        f"   ℓ_k = {D['l_k'].to('nm'):~.2fP}   -> amplitude window l_k << a < delta_max",
        f"  trapped beads = {D['trapped']}  ([P2] Fig.1A: force sensors at both ends"
        f" + driven centre)",
    ]
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}"
        f"  = {lg.ratio('times', 'dt', 'tau_B'):.3e} τ_B   (set by the fastest bending mode)",
        f"  ω       = {D['omega']:.0f} rad/s"
        f"  = {D['omega']/(2*math.pi):.1f} Hz   →  De = ω τ_max ="
        f" {lg.ratio('times', 'tau_max', 'tau_w'):.3f}"
        f"   (old definition omega*tau_chain ="
        f" {lg.ratio('times', 'tau_chain', 'tau_w'):.3f})",
        f"  SNR     = |ŷ|/ℓ_k = {lg.ratio('lengths', 'y_resp', 'l_k'):.3f}"
        f"   (drive a/l_k = {lg.ratio('lengths', 'a', 'l_k'):.2f} — this is that with"
        f" the numerator replaced by the response)",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}  ({N_CYCLES:g} cycles)",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}   × n={sys_['n']}",
        f"  ⚠ Cost: prod {n_prod:,} steps. kappa_theta = "
        f"{lg.ratio('energies', 'kappa_theta', 'kT'):.2e} kT makes the chain stiff, so"
        f" the fastest mode — **unrelated to the measurement** — sets dt.",
        f"    Cost ∝ 1/omega, so the low-omega end dominates. Options — (a) accept it"
        f" (b) lower kappa_0 (surfactant, [P2] Fig.4)",
        f"    (c) measure only high omega and replace low omega with the quasi-static"
        f" limit.",
    ]
    return inp, der, plan


def sweep_specs(dt_scale: float = 1.0, cycles: float = N_CYCLES, samples: int = 2000):
    """Build the spec documents for the 7 sweep points **in memory** and return them
    (writes no file).

    ★ If a hard check fails, `--spec` writes nothing to `specs/` (correct behaviour,
    rule 2).
    But the verification, diagnosis and visualization scripts need the parameters — if
    they depended on the spec file, a failed check would also mean "no figure either".
    Hence this one access point.
    """
    import argparse as _ap
    sys_ = load_system(ROOT / "intake" / "chain-bend-2d-oscill" / "system.yaml")
    args = _ap.Namespace(dt_scale=dt_scale, cycles=cycles, samples=samples)
    lo, hi = sys_["omega_range"]
    out = []
    for om in np.geomspace(lo, hi, N_SWEEP):
        _, spec, _, _, _, _ = build_spec(sys_, float(om), args)
        out.append({"params": spec.params, "numerics": spec.numerics})
    return sorted(out, key=lambda s: s["params"]["omega_star"])


# ════════════════════════════════════════════════════════════════════════
# main — through L3
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, omega, args):
    lg = build_ledger(sys_, omega, dt_scale=args.dt_scale, n_cycles=args.cycles)
    D = lg.derived
    dt = lg.get("times", "dt")
    # Equilibration: relax the chain **before** driving. ★★ The scale is **tau_max**
    # (the longest relaxation).
    # It used to be 20 tau_chain, which is only 2.2 tau_max, leaving an 11% residual
    # transient and making K* wrong by up to 21% (measured: raising equilibration to
    # 10 tau_max shrank the block scatter 1000x).
    n_eq = int(round(20 * float((D["tau_max"] / dt).to(""))))
    n_prod = int(round(float((D["T_obs"] / dt).to(""))))
    sample_every = max(1, n_prod // args.samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(sys_, lg)      # ★ built once — the spec and the report
    tag = f"w{omega:.0f}"                          #    must see the same object or they diverge
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"
    if args.cycles != N_CYCLES:
        tag += f"-nc{args.cycles:g}"

    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_beads": sys_["n"], "trapped": D["trapped"],
                "L_chain_star": lg.ratio("lengths", "L_chain", "d"),
                # All reduced stiffnesses are in kT/d^2 (the same number as
                # ledger energy / kT)
                "kappa_theta_star": lg.ratio("energies", "kappa_theta", "kT"),
                "k_t_star": lg.ratio("energies", "k_t_d2", "kT"),
                "k_bond_star": lg.ratio("energies", "k_b_d2", "kT"),
                "amp_star": lg.ratio("lengths", "a", "d"),
                "omega_star": float((Q(omega, "1/s") * D["tau_B"]).to("dimensionless").magnitude),
                "De": lg.ratio("times", "tau_max", "tau_w"),
                "De_chain_old": lg.ratio("times", "tau_chain", "tau_w"),
                "kappa_drive_star": lg.ratio("energies", "kappa_drive_d2", "kT"),
                "snr_response": lg.ratio("lengths", "y_resp", "l_k"),
                "M_c_star": lg.ratio("energies", "M_c", "kT"),
                "n_trapped": sys_["n_trapped"],
                "bending_impl": BENDING_IMPL},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "dt_over_tau_fast": args.dt_scale * C.GATE,
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": 20260804},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_eq, n_prod


def emit(sys_, omega, args) -> int:
    lg, spec, groups, checks, n_eq, n_prod = build_spec(sys_, omega, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_eq, n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  ω={omega:.0f} rad/s   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)

    if spec.errors:
        print(f"\n❌ {len(spec.errors)} L3 integrity error(s) — the "
              f"non-dimensionalization does not hold.")
        return 1
    if verdict == "FAIL":
        print("\n❌ A hard separation check failed — no spec is written.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 spec: {p.relative_to(ROOT)}")
        return 0

    # ── L4 — run by **re-reading** the on-disk spec (that is when the hash check
    #    fires) ──
    outdir = ROOT / "runs" / run_id
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.txt").write_text(report)
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="the L3 report only")
    ap.add_argument("--spec", action="store_true", help="the L3 spec -> specs/<run_id>.json")
    ap.add_argument("--omega", type=float, default=None,
                    help="drive angular frequency [rad/s]. Default is the lowest sweep "
                         "value (the most expensive point)")
    ap.add_argument("--sweep", action="store_true", help="the whole omega sweep")
    ap.add_argument("--cycles", type=float, default=N_CYCLES, help="cycles observed")
    ap.add_argument("--force", action="store_true", help="re-run a completed run")
    ap.add_argument("--dt-scale", type=float, default=1.0, help="dt multiplier (for a convergence check)")
    ap.add_argument("--samples", type=int, default=2000, help="total number of samples, NOT samples per cycle")
    ap.add_argument("--run", action="store_true",
                    help="run L4. ★ Refused only if BENDING_IMPL is angle_harmonic "
                         "(see the module docstring ⑤)")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/chain-bend-2d-oscill/system.yaml")
    lo, hi = sys_["omega_range"]

    if args.sweep:
        omegas = list(np.geomspace(lo, hi, N_SWEEP))
    else:
        omegas = [args.omega if args.omega is not None else lo]

    if not (args.report or args.spec or args.run):
        print("choose what to do -- `--report`, `--spec` or `--run`")
        return 3

    rc = 0
    for i, om in enumerate(omegas):
        if i:
            print()
        rc |= emit(sys_, float(om), args)
    return rc




# ════════════════════════════════════════════════════════════════════════
# L4 — build the system from the spec alone (bdbot.run runs it, bdbot.health judges)
#
# The construction is exactly what the two gates in
# `scratch/verify_chain_bend_gates.py` verified:
#   · trap = **a ghost particle + a harmonic bond** (r0=0). No custom force needed.
#   · drive = an updater moving the centre ghost as `y = a sin(omega*t)`
#   · estimation = `bdbot.lockin` (gate A compares it against the analytic solution
#     to within 3 sigma)
#   · no pair force — at kappa_theta=1.4e6 kT the chain cannot bend far enough to
#     self-contact
#     (the amplitude-window upper bound delta_max guarantees it). The minimum
#     non-bonded distance is monitored instead.
#
# ★★ The estimator does **not use the nominal amplitude** — the ghost y is measured
#    alongside and the measured phasor is used
#    (bdbot/lockin.py docstring: using the nominal value gets even the sign wrong at
#    De=92).
# ════════════════════════════════════════════════════════════════════════
BENDING_IMPL = "custom_linear"   # "custom_linear" | "angle_harmonic"
# ★ Which bending implementation is in use. The sin(theta) clamp constraint of
#   `angle.Harmonic` is a constraint of **that implementation**,
#   not of the system, so the hard check is made conditional on the implementation
#   (narrowed in scope rather than deleted).
UPDATE_EVERY = 100      # ghost move period. The ZOH is left in and cancelled by the measured phasor


def make_bending_force(A, n_beads):
    """Linearized bending implemented directly as `md.force.Custom`.

    **Because `angle.Harmonic` cannot be used here.**

    ★★ Why not the built-in — `md.angle.Harmonic` clamps `sin θ` at SMALL=1.414e-3,
       and below that the force is shrunk by `sinθ/SMALL` (force ∝ κ(θ−π)² —
       quadratic, not linear). **The energy is exact to 0.000%**, so an energy check
       cannot catch it.
       In this system **all 23** angles are in the broken regime
       (max|θ−π| = 9.4e-4 < 1.414e-3).
       Raising the amplitude cannot escape it either — the angle profile spans a factor
       20 and δ_max blocks the way.
       `angle.Table` has the same problem; `CosineSquared` is quartic at θ₀=π and is
       ruled out.
       Reproduction: `scratch/verify_angle_force_small_theta.py`

    What is implemented is **exactly the linearized form the L3 ledger uses**:
        U = ½ κ_θ Σ_i θ_i²,  θ_i = (y_{i+1} − 2y_i + y_{i−1})/ℓ  →  F_y = −A y,  F_x = 0
    `A = κ_θ BᵀB` IS `bending_matrix()`, so **model and implementation agree exactly**
    (λ_max, τ_fast and driven_static_stiffness all come out of this same A).
    B annihilates constants and linear terms, so it is **invariant under translation
    and (linearized) rotation**.

    ⚠️ At large deformation this differs from a true angle potential. Here
       max|θ| ~ 1e-3, deep inside the linear regime, and L3's 'small-angle linear'
       check gates exactly that.
    """
    import hoomd.md as md

    class Bending(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=False)
            self.A = np.ascontiguousarray(A, dtype=float)
            self.n = int(n_beads)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)   # ★ tag indexing is mandatory
                pos = np.array(snap.particles.position, copy=True)
                m = tags < self.n                                # exclude ghosts
                y = np.zeros(self.n)
                y[tags[m]] = pos[m, 1]
                fy = -(self.A @ y)
                arr.force[:] = 0.0
                arr.potential_energy[:] = 0.0
                arr.force[m, 1] = fy[tags[m]]
                # Distributed per particle so the total is 1/2 y^T A y (the split
                # itself is a convention)
                arr.potential_energy[m] = -0.5 * y[tags[m]] * fy[tags[m]]

    return Bending()


def _move_ghost_action(ghost_tag, amp, omega, dt):
    import hoomd

    class MoveGhost(hoomd.custom.Action):
        def act(self, timestep):
            y = amp * math.sin(omega * timestep * dt)
            with self._state.cpu_local_snapshot as snap:
                tags = np.array(snap.particles.tag, copy=True)
                loc = np.flatnonzero(tags == ghost_tag)
                if len(loc):
                    snap.particles.position[loc[0], 1] = y
    return MoveGhost()


def _chain_frame(n, trapped, ell, L_chain):
    import gsd.hoomd

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:                                  # ghosts are placed on top of their beads
        pos.append(list(pos[g]))
        typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.particles.mass = [1.0] * len(pos)                # BD does not use mass (trap 5)
    f.configuration.box = [4.0 * L_chain] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(trapped)]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(trapped)
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])
    return f


def assert_angle_force_valid(params: dict) -> tuple[float, float]:
    """★★ Check that `angle.Harmonic`'s force is in its valid regime; refuse if not.

    Makes **the same judgment as the L3 hard check, a second time in the builder**
    (defence in depth). Why:
    the L3 check stops a spec from being written and `RUN.execute` refuses a FAILed
    spec, but the paths that call the builder **directly** (tests, manual calls, old
    specs) are stopped by nothing. And
    this failure is not a crash but **silently wrong physics** — the chain comes out
    softer than it really is.
    """
    n = int(params["n_beads"])
    ell = float(params["L_chain_star"]) / (n - 1)
    _, prof = driven_response(n, float(params["kappa_theta_star"]), ell,
                              float(params["k_t_star"]), 1.0,
                              float(params["omega_star"]), float(params["amp_star"]))
    th_hi, th_lo = response_angles(prof, ell)
    if th_lo < ANGLE_SIN_SMALL:
        raise ValueError(
            f"angle.Harmonic's force is outside its valid regime — refusing to run.\n"
            f"  min|θ−π| = {th_lo:.3e} < SMALL = {ANGLE_SIN_SMALL:.3e}  "
            f"(max|θ−π| = {th_hi:.3e})\n"
            f"  HOOMD clamps sin(theta) at SMALL, so below it the force is shrunk by\n"
            f"  sin(theta)/SMALL (force ∝ kappa*(theta−pi)^2 — quadratic, not linear).\n"
            f"  **The energy stays exact** ->\n"
            f"  verifying with energy passes while only the force is wrong. The chain\n"
            f"  comes out softer than it really is.\n"
            f"  Reproduction: scratch/verify_angle_force_small_theta.py "
            f"(skill bd-hoomd trap 15)\n"
            f"  Ways out: (1) implement bending directly with force.Custom (exact but\n"
            f"      26x slower)\n"
            f"      (2) lower kappa_0 to soften the chain ([P2], surfactant) -> theta\n"
            f"      grows\n"
            f"      (3) this regime is linear, so solve it analytically without MD\n"
            f"      (0.32% against the exact minimization)")
    return th_hi, th_lo


@RUN.builder("chain-bend-2d-oscill")
def build(spec, outdir=None) -> RUN.Build:
    # ★ Defence in depth — refuse if the implementation the spec declares differs from
    #   the one this builder actually uses. Running an `angle.Harmonic` spec with the
    #   custom force is silently different physics.
    _impl = spec.params.get("bending_impl", "angle_harmonic")
    if _impl != BENDING_IMPL:
        raise ValueError(f"the spec declares bending implementation '{_impl}' but the "
                         f"builder is "
                         f"'{BENDING_IMPL}'. Regenerate the spec.")
    if _impl == "angle_harmonic":
        assert_angle_force_valid(spec.params)
    # ★ `assert_angle_force_valid` is called above **only when the implementation is
    #   angle_harmonic**. The unconditional call that used to be here was removed —
    #   the custom force has no sin(theta) clamp, so applying that verdict as-is would
    #   reject a valid configuration.
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_beads"])
    trapped = sorted(int(t) for t in P["trapped"])
    L_chain = float(P["L_chain_star"])
    ell = L_chain / (n - 1)
    k_t, k_b, kth = float(P["k_t_star"]), float(P["k_bond_star"]), float(P["kappa_theta_star"])
    amp, omega = float(P["amp_star"]), float(P["omega_star"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    mid = trapped[len(trapped) // 2]
    ghost_mid = n + trapped.index(mid)

    sim = SIM.make_sim(_chain_frame(n, trapped, ell, L_chain), seed=seed)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=k_b, r0=ell)
    bond.params["trap"] = dict(k=k_t, r0=0.0)          # trap = a harmonic bond to the ghost
    # ★ angle.Harmonic gets the force wrong in this regime (see make_bending_force)
    angle = make_bending_force(bending_matrix(n, kth, ell), n)
    # ★ Ghosts are not integrated — the updater writes their position directly. Only
    #   the beads get BD.
    bd = md.methods.Brownian(filter=hoomd.filter.Type(["A"]), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[bond, angle])
    integ.integrate_rotational_dof = False             # trap 3
    sim.operations.integrator = integ
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=_move_ghost_action(ghost_mid, amp, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, int(Nm["n_prod"]) // 200))

    def pe_per_particle():
        return float(np.array(bond.energies).sum() + np.array(angle.energies).sum()) / n

    def sample(timestep, phase):
        snap = sim.state.get_snapshot()
        p = np.array(snap.particles.position, dtype=float)
        # ★ The ghost y is measured **alongside** — without it the ZOH passes straight
        #   through as a K* error
        nb = p[:n, :2]
        sep = np.linalg.norm(nb[:, None, :] - nb[None, :, :], axis=-1)
        np.fill_diagonal(sep, np.inf)
        for i in range(n - 1):                          # exclude bonded neighbours (this watches for overlap)
            sep[i, i + 1] = sep[i + 1, i] = np.inf
        return {"t": timestep * dt, "y_bead": float(p[mid, 1]),
                "y_ghost": float(p[ghost_mid, 1]),
                "y_end0": float(p[trapped[0], 1]), "y_end1": float(p[trapped[-1], 1]),
                "min_sep_nonbonded": float(sep.min()),
                "max_theta": float(np.abs(np.diff(p[:n, 1], n=2)).max() / ell)}

    def finalize(cols):
        t, yb, yg = cols["t"], cols["y_bead"], cols["y_ghost"]
        blocks_y = LI.lockin_blocks(t, yb, omega, n_blocks=10)
        blocks_g = LI.lockin_blocks(t, yg, omega, n_blocks=10)
        Kb = np.array([LI.k_star(a, b, k_t, omega) for a, b in zip(blocks_y, blocks_g)])
        K, Ksem = LI.agg(Kb)
        yh, _ = LI.agg(blocks_y)
        gh, _ = LI.agg(blocks_g)
        # Also computed with the nominal amplitude, to **leave the difference on
        # record** (after-the-fact evidence of the ZOH trap)
        K_nom, _ = LI.agg(np.array([LI.k_star(a, complex(amp), k_t, omega) for a in blocks_y]))
        # Static-limit prediction — **derived from the model I implemented**, hence
        # implementation_check
        K_static = driven_static_stiffness(n, kth, ell, k_t)
        de = float(P["De"])
        obs = [
            MET.observable("K_prime", float(K.real), K_static if de < 1.5 else None,
                           "kT/d^2", role="implementation_check" if de < 1.5 else "measurement",
                           tol_pct=15.0,
                           note=f"storage stiffness K'(omega). De={de:.3f}. "
                                f"Static-limit prediction "
                                f"{K_static:.4g} (driven_static_stiffness — includes "
                                f"the finite trap stiffness). "
                                f"Compared only for De<1.5"),
            MET.observable("K_doubleprime", float(K.imag), None, "kT/d^2",
                           role="measurement", note="loss stiffness K''(omega)"),
            MET.observable("K_sem", Ksem, None, "kT/d^2", role="measurement",
                           note="scatter over 10 blocks (the larger of real and "
                                "imaginary)"),
            MET.observable("y_response", float(abs(yh)), None, "d", role="measurement",
                           note=f"response amplitude |y_hat|. Drive "
                                f"|y_hat_c|={abs(gh):.5g} (nominal {amp:g})"),
        ]
        return {"observables": obs,
                "extra": {"K_prime": float(K.real), "K_doubleprime": float(K.imag),
                          "K_sem": Ksem, "K_static_pred": K_static,
                          "K_prime_nominal_amp": float(K_nom.real),
                          "nominal_vs_measured_pct":
                              100 * (K_nom.real - K.real) / abs(K.real) if K.real else None,
                          "y_resp_abs": float(abs(yh)), "drive_abs": float(abs(gh)),
                          "drive_over_nominal": float(abs(gh) / amp),
                          "De": de, "omega_star": omega,
                          "min_sep_nonbonded": float(cols["min_sep_nonbonded"].min()),
                          "max_theta": float(cols["max_theta"].max())}}

    every = max(1, int(Nm["sample_every"]))
    return RUN.Build(
        sim=sim, forces=[bond, angle], n_particles=n,
        sample=sample, pe_per_particle=pe_per_particle, sample_every=every,
        phases=[RUN.Phase("warm-up", int(Nm["n_eq"]), collect=False,
                          note="drive ON, samples discarded (20 tau_max — 2.2 tau_max "
                               "left an 11% transient)"),
                RUN.Phase("production", int(Nm["n_prod"]), every,
                          note=f"lock-in collection, omega*={omega:.4g}, "
                               f"De={P['De']:.3f}")],
        tags=["2D", "chain", "bending", "angle_harmonic", "oscillatory_drive",
              "microrheology", "lockin", "newtonian"],
        physical={"n_beads": n, "De": float(P["De"]), "omega_star": omega,
                  "kappa_theta_star": kth, "k_t_star": k_t, "amp_star": amp,
                  **{k: v for k, v in spec.raw["back_transform"].items()
                     if isinstance(v, float)}},
        finalize=finalize)


if __name__ == "__main__":
    sys.exit(main())
