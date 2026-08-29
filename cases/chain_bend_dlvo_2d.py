"""`chain-bend-2d-dlvo` -- the alternative-hypothesis branch of
chain-bend-2d-oscill.

Same three-point-bending, oscillatory geometry (fixed traps at both ends plus an
oscillating central trap), but the bead bond is made from **the secondary minimum
of a DLVO central pair potential** rather than [P1][P2]'s JKR adhesive contact plus
bending stiffness. **There is no explicit angular (bending) potential -- that
absence is the hypothesis.**

[P1] Pantina & Furst, PRL 94, 138301 (2005), p.2 left column predicts it
explicitly: such assumptions have been made based on DLVO interactions between
particles, which are centro-symmetric; if particles did undergo free rotations, we
would expect the aggregates to respond to the bending moment by forming a
trianglelike structure. The paper found that prediction contradicted by experiment
and introduced JKR adhesive contact. This case runs the prediction side (rules 6
and 7' -- a mismatch is a discovery, not a failure).

** The structural fact, derived before running
   (`intake/chain-bend-2d-dlvo/observation.yaml` G1):
in a straight chain, if a pure central-force bond sits at its **natural length
(U'=0)**, the bond energy varies as O(y^4) in the transverse displacement y -- the
coefficient of the O(y^2) term is U'(ell), which is zero where the force is zero.
So **the linear bending stiffness is exactly 0.** There is nothing in this system
corresponding to chain-bend-2d-oscill's bending_matrix(). Under three-point
bending, the "trianglelike" local kinking [P1] described is more likely than smooth
beam curvature -- and that is confirmed at L4 rather than assumed.

Physical parameter provenance: [P1] gives d=1.47um and psi0=40mV; MgCl2 10mM
(ionic strength); Hamaker A=1.05e-20 J (from a web search, tier 3 -- the primary
source was not confirmed, and that uncertainty is stated).
DLVO curve: HHF (constant-potential, weak-overlap) EDL + Derjaguin non-retarded
vdW. Computed in `verify/dlvo_ledger.py` (barrier 382kT @ 0.49nm, secondary minimum
-11.7kT @ 11.2nm).

Potential implementation: **a global pair potential, not a bond list**
(`md.pair.Table`, applied to every neighbour pair -- the same pattern as
`soft-r3-2d-A-sweep`). A "bond" arises because adjacent beads sit in each other's
secondary minimum; it is not a predeclared bond topology. That IS the hypothesis
that a chain can be made from a central pair potential alone. The cutoff
(h <~ 60nm) is far shorter than the non-adjacent bead spacing (~3um), so in the
normal state only nearest neighbours interact -- but if the chain folds, distant
beads interact naturally too, and that is also real physics, so it is not
prevented.
"""
from __future__ import annotations

import argparse
import json
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
N_CYCLES = 50.0
NA = 6.02214076e23
E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
R_WCA = 2 ** (1 / 6)          # WCA cutoff (in sigma=d) -- the same convention as bdbot/pairpot.py
CUTOFF_H_STAR = 0.06          # table cutoff -- surface gap h/d. 8x outside the
                              # secondary minimum (~0.0076)


# ════════════════════════════════════════════════════════════════════════
# 1. the physical system (SI)
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
        "psi0": P(raw["particle"]["surface_potential"]),
        "n_list": [int(x) for x in raw["particle"]["count"]["value"]],
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "eps_r": P(raw["medium"]["relative_permittivity"]),
        "ionic_strength": P(raw["medium"]["ionic_strength"]),
        "A_H": P(con["hamaker_constant"]),
        # * the control -- used only with `--jkr`. interactions[1] is OFF by default.
        "kappa_theta": P(raw["interactions"][1]["angle_stiffness"]),
        "k_t": P(raw["external"]["stiffness"]),
        "amp_range": [float(x) for x in raw["external"]["amplitude"]["value"]],
        "omega_range": [float(x) for x in raw["external"]["omega_range"]["value"]],
        "n_trapped": int(raw["external"]["n_trapped"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# 2. DLVO(h) -- determined by three reduced variables in h*=h/d
#    (kappa_star, edl_amp, vdw_amp).
#    * Not written from memory -- this is the reduced form of exactly the expression
#    (HHF constant-potential weak-overlap EDL + Derjaguin non-retarded vdW) that was
#    first computed and verified in SI in verify/dlvo_ledger.py.
# ════════════════════════════════════════════════════════════════════════
def dlvo_reduced_params(sys_: dict) -> dict:
    d = sys_["d"].value.to("m").magnitude
    a = d / 2
    T = sys_["T"].value.to("K").magnitude
    kT = 1.380649e-23 * T
    eps_r = sys_["eps_r"].value.to("dimensionless").magnitude
    psi0 = sys_["psi0"].value.to("V").magnitude
    A_H = sys_["A_H"].value.to("J").magnitude
    c_salt = sys_["ionic_strength"].value.to("mol/m^3").magnitude   # MgCl2 molarity
    I_SI = 0.5 * (c_salt * (2 ** 2) + 2 * c_salt * (1 ** 2)) * NA   # 1/m^3 (ionic strength, MgCl2 -> Mg2+ + 2Cl-)
    kappa = math.sqrt(2 * I_SI * E_CHARGE ** 2 / (EPS0 * eps_r * kT))
    a_star = 0.5
    return {
        "a_star": a_star,
        "kappa_star": kappa * d,
        "edl_amp": 2 * math.pi * EPS0 * eps_r * a * psi0 ** 2 / kT,   # the coefficient in front of U_edl*
        "vdw_amp": A_H / (12 * kT),                                    # U_vdw* = -vdw_amp*a_star/h*
        "kT": kT, "d": d,
    }


def U_star(h_star, p: dict):
    """U(h)/kT for h*=h/d>0. EDL (HHF constant-potential, weak-overlap) plus vdW
    (Derjaguin, non-retarded).
    """
    h_star = np.asarray(h_star, dtype=float)
    u_edl = p["edl_amp"] * np.log1p(np.exp(-p["kappa_star"] * h_star))
    u_vdw = -p["vdw_amp"] * p["a_star"] / h_star
    return u_edl + u_vdw


def F_h_star(h_star, p: dict):
    """-dU*/dh* -- the force in the direction of increasing h* (positive = repulsive).
    The analytic derivative.
    """
    h_star = np.asarray(h_star, dtype=float)
    k = p["kappa_star"]
    dU_edl = p["edl_amp"] * (-k) * np.exp(-k * h_star) / (1 + np.exp(-k * h_star))
    dU_vdw = p["vdw_amp"] * p["a_star"] / h_star ** 2
    return -(dU_edl + dU_vdw)


def find_well(p: dict) -> dict:
    """The barrier and secondary-minimum positions and depths, plus the well curvature
    (= the bond radial stiffness, in kT/d^2), via bisection and central differences.
    """
    hs = np.geomspace(1e-4, CUTOFF_H_STAR * 3, 4000)
    Us = U_star(hs, p)
    ibar = int(np.argmax(Us[: int(len(hs) * 0.5)]))
    barrier_h, barrier_U = float(hs[ibar]), float(Us[ibar])
    tail = Us[ibar:]
    imin = int(np.argmin(tail))
    h_min, U_min = float(hs[ibar + imin]), float(tail[imin])
    dh = h_min * 1e-4
    k_bond_star = (U_star(h_min + dh, p) - 2 * U_star(h_min, p) + U_star(h_min - dh, p)) / dh ** 2
    return {"barrier_h": barrier_h, "barrier_U": float(barrier_U), "h_min": h_min,
            "U_min": float(U_min), "k_bond_star": float(k_bond_star)}


def trapped_indices(n: int) -> list[int]:
    return [0, n // 2, n - 1]


def _bow(y: np.ndarray) -> float:
    """Bow -- the maximum deviation from the chord joining the two ends [d]. Pure
    bending, with rigid translation and tilt removed.

    * Why this is needed: this system's traps have finite stiffness (k_t), so the
      whole chain gets pushed. Measuring raw y lets that rigid-body motion mask the
      bending signal. Measured (n=9, omega=3000, a=632nm):
          full y range   DLVO 0.303 d  vs JKR 0.090 d  ->  3.4x
          bow            DLVO 0.1175 d vs JKR 0.0060 d ->  **19.6x**
      Same data, and the discriminating power differs by 5.8x.

    WARNING: this holds **only under a soft trap (free deformation).** Under
      position-forced driving, the bow ratio collapses from 15.4x to 1.4x -- when the
      deformation is imposed, shape carries no information. Free deformation ->
      measure shape; imposed deformation -> measure force.
    """
    n = len(y)
    base = np.linspace(0.0, 1.0, n) * (y[-1] - y[0]) + y[0]
    return float(np.abs(y - base).max())


# ════════════════════════════════════════════════════════════════════════
# 3. the scale ledger
# ════════════════════════════════════════════════════════════════════════
def bending_matrix(n: int, kappa_theta: float, ell: float) -> np.ndarray:
    """The quadratic form A = kappa_theta * B^T B of
    U = (1/2)*kappa_theta*sum(theta_i^2), theta_i=(y_{i+1}-2y_i+y_{i-1})/ell.

    * **The same expression** as chain-bend-2d-oscill, where the paper's kappa_0 was
      reproduced and the discrete-to-continuum mapping verified. Used here only for
      the control (`--jkr`).
    """
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def build_ledger(sys_, n: int, omega: float, amp_nm: float, *, dt_scale=1.0,
                  n_cycles=N_CYCLES, jkr: bool = False,
                  kt_scale: float = 1.0) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B = b["kT"], b["gamma"], b["tau_B"]
    # * kt_scale -- deliberately stiffen the trap to push the system toward a
    #   **position-control condition.**
    #   Why it is needed: with the default k_t (sketch tier 0, 10 pN/um) the driven
    #   bead does not follow the trap position at all (measured tracking |y_hat|/a --
    #   JKR 3.5%, DLVO 29%). The chain and the viscosity are stronger than the trap,
    #   so the trap simply stretches and the bead is not pulled along.
    #   Measured design table (n=9, omega=3000, JKR):
    #       scale       1   10   30  100  300  1000
    #       tracking 0.035 0.33 0.59 0.82 0.93 0.98
    #       dt cost   1.00 1.00 1.00 1.01 1.02 1.08 (JKR) / 1.0 1.0 - 2.0 6.0 20 (DLVO)
    #   * Use 100x: tracking 0.82 at almost no cost. From 300x the cancellation in the
    #     K' estimator becomes dangerous -- K' = k_t(y_hat_c/y_hat - 1), and
    #     (y_hat_c/y_hat - 1) shrinks to 0.09, amplifying the relative noise (0.027 at
    #     1000x). **Better tracking makes the measurement worse.**
    #   WARNING on physical honesty: 100x = 1000 pN/um, which is 25x stiffer than a
    #     real optical trap ([P1] ~40 pN/um). This is **a numerical experiment in the
    #     position-control limit**, not the system in the sketch.
    #     (Same context as trap-drag recording that a stiff trap = a position-control
    #     condition.)
    k_t = sys_["k_t"].value.to("N/m") * float(kt_scale)
    amp = Q(amp_nm, "nm").to("m")

    p = dlvo_reduced_params(sys_)
    w = find_well(p)
    ell = d * (1 + w["h_min"])                       # natural bond length (centre to centre, U'=0)
    L_chain = (n - 1) * ell
    k_bond = Q(w["k_bond_star"], "dimensionless") * kT / d ** 2   # radial (stretch) stiffness
    sigma_bond = (kT / k_bond) ** 0.5

    tau_bond = C.relaxation_time(gamma, k_bond)      # DLVO bond stretch
    tau_k = C.relaxation_time(gamma, k_t)

    # * The control (`--jkr`): with the bending term on, the largest eigenvalue of the
    #   stiffness matrix may set dt. That is what happened in chain-bend-2d-oscill
    #   (tau_fast=0.28us dominated) -- so here too dt is taken from **whichever is
    #   faster.** Otherwise it diverges silently.
    kappa_theta = None
    tau_bend_fast = None
    if jkr:
        kappa_theta = sys_["kappa_theta"].value.to("J")
        kth_star = float((kappa_theta / kT).to("dimensionless").magnitude)
        A_bend = bending_matrix(n, kth_star, float((ell / d).to("dimensionless").magnitude))
        for i in trapped_indices(n):
            A_bend[i, i] += float((k_t * d ** 2 / kT).to("dimensionless").magnitude)
        lam_max_star = float(np.linalg.eigvalsh(A_bend)[-1])          # in kT/d^2
        tau_bend_fast = C.relaxation_time(gamma, Q(lam_max_star, "dimensionless") * kT / d ** 2)
    # ** dt is set by **the fastest mode** -- all three candidates go in.
    #   This was a hole: tau_k (the trap) used to be omitted and only tau_bond and
    #   tau_bend were considered. Stiffening the trap with --kt-scale makes k_t exceed
    #   k_bond (around 200x, where k_bond*=1.04e6), at which point the trap becomes the
    #   fastest mode while dt stays put -- and dt becomes **silently insufficient.**
    #   (The DLVO branch was especially exposed, since with no bending matrix tau_bend
    #   is None.)
    cands = [tau_bond, tau_k] + ([tau_bend_fast] if tau_bend_fast is not None else [])
    tau_fast = min(cands, key=lambda q: float(q.to("s").magnitude))
    dt = dt_scale * C.dt_from_gate(tau_fast)

    # ** The transverse (bending) linear stiffness is structurally 0 (natural-length
    #    equilibrium plus central forces) -- no invented substitute scale is used.
    #    Instead the whole-chain shape relaxation is roughly estimated as a "contour
    #    length diffusion time" (a *proposed upper bound, not an exact Rouse
    #    spectrum).
    tau_chain_diffusion = (L_chain ** 2 / b["D_t"]).to("s")

    tau_w = Q(1.0 / omega, "s")
    tau_period = Q(2 * math.pi / omega, "s")
    T_obs = (n_cycles * tau_period).to("s")

    lg = SC.ScaleLedger()
    lg.add_length("sigma_bond", sigma_bond.to("m"), "bond radial thermal width sqrt(kT/k_bond)", star=True)
    lg.add_length("h_min", Q(w["h_min"], "dimensionless") * d, "secondary-minimum position (surface gap)", star=True)
    lg.add_length("a", amp, "drive amplitude", star=True)
    lg.add_length("d", d, "bead diameter")
    lg.add_length("ell", ell.to("m"), "natural bond length (centre to centre, d+h_min)")
    lg.add_length("L_chain", L_chain.to("m"), "chain contour length (n-1)*ell")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_bond", tau_bond, "* gamma/k_bond bond stretch -- the only linear stiffness mode. It sets dt",
                star=True)
    lg.add_time("tau_k", tau_k, "gamma/k_t trap")
    lg.add_time("tau_w", tau_w, f"1/omega drive (omega={omega:.0f} rad/s)")
    lg.add_time("tau_period", tau_period, "2*pi/omega drive period")
    lg.add_time("tau_chain_diff", tau_chain_diffusion,
                "*proposed: L_chain^2/D_t -- a rough upper bound on chain shape relaxation "
                "(with zero bending stiffness, deriving an exact Rouse spectrum needs "
                "separate verification; this is an estimate)", star=True)
    lg.add_time("tau_B", tau_B, "d^2/D_t diffusion (reference)")
    lg.add_time("T_obs", T_obs, f"observation window ({n_cycles:g} cycles)", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("k_t_d2", (k_t * d ** 2).to("J"), "k_t*d^2 trap stiffness")
    lg.add_energy("k_bond_d2", (k_bond * d ** 2).to("J"), "k_bond*d^2 bond radial stiffness", star=True)
    lg.add_energy("well_depth", Q(-w["U_min"], "dimensionless") * kT,
                  "|secondary-minimum depth| -- the thermal scale at which bonding can be reversible", star=True)
    lg.add_time("tau_fast", tau_fast, "the fastest mode, which sets dt", role="", star=True)
    lg.declare_absent(
        "box",
        "no periodic boundaries (one chain, with traps fixing its position). Same reason "
        "as chain-bend-2d-oscill.")
    if jkr:
        lg.add_energy("kappa_theta", kappa_theta,
                      "* the control: JKR tangential bending stiffness kappa_theta = EI/ell ([P1][P2])", star=True)
        lg.add_time("tau_bend_fast", tau_bend_fast,
                    "the fastest mode from the largest eigenvalue of the bending stiffness matrix -- it can set dt", star=True)
    else:
        lg.declare_absent(
            "bending_stiffness",
            "** structurally absent (no invented substitute value goes here) -- for a "
            "straight chain with pure central forces and bonds at their natural length "
            "(U'=0), the linear (O(y^2)) restoring force against transverse "
            "displacement is exactly 0, because the O(y^2) coefficient is U'(ell) and "
            "that is zero at the natural length. There is nothing in this system "
            "corresponding to chain-bend-2d-oscill's bending_matrix/lambda_min -- that "
            "absence IS the hypothesis (G1), confirmed from the shape of the L4 "
            "trajectory (smooth curvature vs local kinking). "
            "* Turning on `--jkr` makes it the control.")
    lg.derived = dict(gamma=gamma, D_t=b["D_t"], m=b["m"], kT=kT, d=d, tau_B=tau_B,
                      ell=ell.to("m"), L_chain=L_chain.to("m"), k_bond=k_bond,
                      sigma_bond=sigma_bond.to("m"), tau_bond=tau_bond, tau_k=tau_k,
                      tau_w=tau_w, tau_period=tau_period, tau_chain_diff=tau_chain_diffusion,
                      dt=dt, T_obs=T_obs, omega=omega, n=n, amp=amp,
                      jkr=jkr, kappa_theta=kappa_theta, tau_fast=tau_fast,
                      tau_bend_fast=tau_bend_fast,
                      k_bond_star=w["k_bond_star"], h_min_star=w["h_min"],
                      U_min_star=w["U_min"], barrier_star=w["barrier_U"],
                      reduced=p, trapped=trapped_indices(n))
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " * Unlike chain-bend-2d-oscill, this system has no "
        "bending linear stiffness at all, structurally -- so the relation between the "
        "scale that sets dt (tau_bond, bond stretch) and what is being measured (chain "
        "shape relaxation, tau_chain_diff) is not fixed in advance. It is confirmed by "
        "running L4.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 4. dimensionless groups + checks
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg, n):
    D = lg.derived
    r = lg.ratio

    groups = [
        ND.Group("well_depth/kT", r("energies", "well_depth", "kT"),
                 ("energies", "well_depth"), ("energies", "kT"), "",
                 "* bond depth -- on the scale of kT it is reversible (easily broken)"),
        ND.Group("k_bond_star", r("energies", "k_bond_d2", "kT"),
                 ("energies", "k_bond_d2"), ("energies", "kT"), "k_bond d²/kT",
                 "bond radial (stretch) stiffness -- NOT a bending stiffness"),
        ND.Group("k_t/k_bond", r("energies", "k_t_d2", "k_bond_d2"),
                 ("energies", "k_t_d2"), ("energies", "k_bond_d2"), "",
                 "trap vs bond stiffness"),
        ND.Group("sigma_bond/h_min", r("lengths", "sigma_bond", "h_min"),
                 ("lengths", "sigma_bond"), ("lengths", "h_min"), "",
                 "bond thermal width / well position -- close to 1 means the bond is unstable"),
        ND.Group("a/sigma_bond", r("lengths", "a", "sigma_bond"),
                 ("lengths", "a"), ("lengths", "sigma_bond"), "",
                 "* drive amplitude vs bond thermal fluctuation -- a proxy for SNR"),
        ND.Group("a/L_chain", r("lengths", "a", "L_chain"),
                 ("lengths", "a"), ("lengths", "L_chain"), "", "amplitude vs chain length"),
        ND.Group("n_beads", float(n), None, None, "", "bead count (input, *suggested sweep)"),
        ND.Group("De_bond", r("times", "tau_bond", "tau_w"),
                 ("times", "tau_bond"), ("times", "tau_w"), "ω τ_bond",
                 "Deborah number referenced on bond stretch (expected very small -- tau_bond is extremely fast)"),
        ND.Group("De_chain_diff", r("times", "tau_chain_diff", "tau_w"),
                 ("times", "tau_chain_diff"), ("times", "tau_w"), "ω τ_chain_diff",
                 "*proposed: Deborah referenced on chain shape relaxation (an estimate)"),
        ND.Group("De_trap", r("times", "tau_k", "tau_w"),
                 ("times", "tau_k"), ("times", "tau_w"), "ω τ_k", "trap reference"),
        ND.Group("tau_bond/tau_chain_diff", r("times", "tau_bond", "tau_chain_diff"),
                 ("times", "tau_bond"), ("times", "tau_chain_diff"), "",
                 "* the scale-separation span (the mode that sets dt vs the shape-relaxation estimate)"),
        ND.Group("dt/tau_bond", r("times", "dt", "tau_bond"),
                 ("times", "dt"), ("times", "tau_bond"), "", "integration resolution"),
        ND.Group("n_cycles", r("times", "T_obs", "tau_period"),
                 ("times", "T_obs"), ("times", "tau_period"), "", "cycles observed"),
        ND.Group("St", r("times", "tau_p", "tau_B"),
                 ("times", "tau_p"), ("times", "tau_B"), "tau_p/tau_B", "inertia vs diffusion"),
    ]
    checks = [
        C.Check("model", "note: tau_p/tau_bond", r("times", "tau_p", "tau_bond"),
                C.GATE, "<=",
                "* the bond well is deep and narrow, so tau_bond approaches tau_p "
                "(exceeding the inertia-negligible criterion by ~2.8x, independent of "
                "n, omega and amplitude -- a property of the bond physics itself). "
                "chain-bend-2d-oscill had the same class of violation (zeta=0.65, far "
                "worse than here) and it was verified harmless (0.16%) by comparing "
                "OverdampedViscous against Langevin(kT=0). **That comparison has not "
                "been run for this system yet.** It is treated as soft while unverified, "
                "and will be verified if needed after the structural smoke test", hard=False),
        C.Check("model", "bond stability     sigma_bond/h_min", r("lengths", "sigma_bond", "h_min"),
                0.5, "<=",
                "** if the bond thermal width exceeds half the well position, the bond "
                "is unstable enough to keep breaking thermally -- that may itself be the "
                "result, but it is flagged in advance",
                hard=False),
        C.Check("integration", "bond stretch resolved dt/tau_bond", r("times", "dt", "tau_bond"),
                C.GATE, "<=", "the DLVO bond stretch mode. Miss it and it diverges"),
        C.Check("integration", "fastest mode resolved dt/tau_fast", r("times", "dt", "tau_fast"),
                C.GATE, "<=",
                "* in the control (--jkr) the largest eigenvalue of the bending "
                "stiffness matrix can be faster -- dt is taken from that (which is what "
                "happened in chain-bend-2d-oscill)"),
        C.Check("integration", "drive resolved       dt/tau_w", r("times", "dt", "tau_w"), C.GATE, "<=",
                f"resolves the drive omega = {lg.derived['omega']:.0f} rad/s"),
        C.Check("statistics", "SNR (bond)         a/sigma_bond", r("lengths", "a", "sigma_bond"), 3.0, ">=",
                "the drive amplitude must exceed the bond thermal fluctuation for the signal to clear the noise", hard=False),
        C.Check("statistics", "cycles observed      T_obs/(2pi/w)", r("times", "T_obs", "tau_period"),
                N_CYCLES, ">=", "cycles used for the phase average", hard=False),
    ]
    return groups, checks


def report_blocks(sys_, lg, n_eq, n_prod, n):
    D = lg.derived
    inp = [R.kv("d", f"{sys_['d'].value:~.4gP}", sys_["d"].tier, sys_["d"].source[:44]),
           R.kv("psi0", f"{sys_['psi0'].value:~.4gP}", sys_["psi0"].tier, sys_["psi0"].source[:44]),
           R.kv("I(MgCl2)", f"{sys_['ionic_strength'].value:~.4gP}",
                sys_["ionic_strength"].tier, sys_["ionic_strength"].source[:44]),
           R.kv("A_H", f"{sys_['A_H'].value:~.4gP}", sys_["A_H"].tier, sys_["A_H"].source[:44]),
           R.kv("n", f"{n}", 3, "*suggested sweep"),
           R.kv("omega", f"{D['omega']:.0f} rad/s", 3, "*suggested sweep"),
           R.kv("amp", f"{D['amp'].to('nm'):~.1fP}", 3, "*suggested sweep")]
    der = [
        f"  lambda_D (radius-independent) -- see verify/dlvo_ledger.py, outside build_ledger",
        f"  bond: barrier {D['barrier_star']:.2f} kT   secondary minimum {D['U_min_star']:.3f} kT"
        f" @ h={D['h_min_star']*float(D['d'].to('nm').magnitude):.2f} nm",
        f"  k_bond = {D['k_bond'].to('pN/um'):~.4fP} = {D['k_bond_star']:.4e} kT/d²"
        f"   σ_bond = {D['sigma_bond'].to('nm'):~.3fP}",
        f"  ell (natural length) = {D['ell'].to('nm'):~.2fP}   L_chain = {D['L_chain'].to('um'):~.3fP}",
        f"  ** linear bending stiffness = 0 (structural, declare_absent) -- there is no "
        f"value in this system corresponding to chain-bend-2d-oscill's lambda_min",
        f"  tau_bond = {D['tau_bond'].to_compact():~.4gP}  (the only mode that sets dt)",
        f"  tau_chain_diff (estimate) = {D['tau_chain_diff'].to_compact():~.4gP}"
        f" = {float(D['tau_chain_diff']/D['tau_bond']):.3e} × τ_bond",
    ]
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}  = {lg.ratio('times','dt','tau_B'):.3e} τ_B",
        f"  ω       = {D['omega']:.0f} rad/s  →  De_bond = {lg.ratio('times','tau_bond','tau_w'):.3e}"
        f"   De_chain_diff (estimate) = {lg.ratio('times','tau_chain_diff','tau_w'):.3f}",
        f"  SNR     = a/σ_bond = {lg.ratio('lengths','a','sigma_bond'):.3f}",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}  ({N_CYCLES:g} cycles)",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}   × n={n}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# main — L3
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, n, omega, amp_nm, args):
    lg = build_ledger(sys_, n, omega, amp_nm, dt_scale=args.dt_scale, n_cycles=args.cycles,
                      jkr=args.jkr, kt_scale=args.kt_scale)
    D = lg.derived
    dt = lg.get("times", "dt")
    # ** Equilibrating for 5x tau_chain_diff (the estimate) makes n_eq explode into
    #   tens of billions of steps (tau_chain_diff/tau_bond ~ 1e7-1e8). With no bending
    #   stiffness the very meaning of "fully equilibrated" is unclear, and the goal
    #   right now is **structural confirmation** (a smoke test), not statistical
    #   convergence -- so it is set cheaply as a multiple of the local (bond)
    #   relaxation. To measure whole-chain shape relaxation, raise --eq-scale and
    #   verify separately (user-confirmed: the smoke test comes first).
    n_eq = int(round(args.eq_scale * float((D["tau_fast"] / dt).to(""))))
    n_prod = int(round(float((D["T_obs"] / dt).to(""))))
    sample_every = max(1, n_prod // args.samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(sys_, lg, n)
    tag = f"n{n}-w{omega:.0f}-a{amp_nm:.0f}"
    if args.jkr:
        tag += "-jkr"
    if args.kt_scale != 1.0:
        tag += f"-kt{args.kt_scale:g}"
    if args.drive_mode != "trap":
        tag += f"-{args.drive_mode}"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"

    p = D["reduced"]
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_beads": n, "trapped": D["trapped"],
                "kappa_star": p["kappa_star"], "edl_amp": p["edl_amp"], "vdw_amp": p["vdw_amp"],
                "a_star": p["a_star"], "cutoff_h_star": CUTOFF_H_STAR,
                "k_t_star": lg.ratio("energies", "k_t_d2", "kT"),
                "amp_star": lg.ratio("lengths", "a", "d"),
                "omega_star": float((Q(omega, "1/s") * D["tau_B"]).to("dimensionless").magnitude),
                "De_bond": lg.ratio("times", "tau_bond", "tau_w"),
                "well_depth_star": D["U_min_star"], "h_min_star": D["h_min_star"],
                "k_bond_star": D["k_bond_star"],
                # * the control switch -- it has to be in the run_id hash for the two
                #   branches to be distinct
                "jkr": bool(args.jkr),
                "drive_mode": args.drive_mode,
                "kappa_theta_star": (lg.ratio("energies", "kappa_theta", "kT")
                                     if args.jkr else 0.0),
                "n_trapped": sys_["n_trapped"]},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": args.seed},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_eq, n_prod


def emit(sys_, n, omega, amp_nm, args) -> int:
    lg, spec, groups, checks, n_eq, n_prod = build_spec(sys_, n, omega, amp_nm, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_eq, n_prod, n)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  n={n} ω={omega:.0f} a={amp_nm:.0f}nm"
              f"   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)

    if spec.errors:
        print(f"\nx {len(spec.errors)} L3 integrity error(s).")
        return 1
    if verdict == "FAIL":
        print("\nx a hard separation check failed -- not writing the spec.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 spec: {p.relative_to(ROOT)}")
        return 0

    outdir = ROOT / "runs" / run_id
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.txt").write_text(report)
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    verdict_txt = RUN.render_verdict(v)
    print(verdict_txt)

    # ** `result.txt` must be written -- it is the project's **completion marker**
    #   (the other three cases all write it). `bdbot.run.execute` does not write it;
    #   it is **the case script's responsibility.** Omit it and three things break
    #   silently:
    #     1. `bdbot.cli status` counts the run as zero (cli.py looks at result.txt)
    #     2. `runid.prepare_outdir` cannot recognize completion and keeps re-running
    #        the same run
    #     3. an "incomplete cleanup" pass **deletes a completed run** -- 6 were
    #        actually destroyed that way
    #   (Added late; the 137 earlier runs were backfilled with
    #   verify/backfill_result_txt.py.)
    if v["status"] != "skipped":
        obs_lines = []
        try:
            mj = json.loads((outdir / "metrics.json").read_text())
            for o in mj.get("observables", []):
                m = o.get("measured")
                p_ = o.get("predicted")
                tail = f"   (predicted {p_:.6g})" if isinstance(p_, (int, float)) else ""
                obs_lines.append(f"  {o['name']:<22} {m:.6g}{tail}" if m is not None
                                 else f"  {o['name']:<22} —")
        except Exception as e:                      # leave the marker even if the result cannot be read
            obs_lines.append(f"  (could not read metrics.json: {e})")
        result = "\n".join(["=" * 84, f"RESULT -- {run_id}", "=" * 84,
                            *obs_lines, "=" * 84, verdict_txt])
        (outdir / "result.txt").write_text(report + "\n" + result)
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n", type=int, default=None, help="bead count (default: the whole list in system.yaml)")
    ap.add_argument("--omega", type=float, default=None, help="rad/s (default: the lowest value in the range)")
    ap.add_argument("--amp", type=float, default=None, help="nm (default: the median of the range)")
    ap.add_argument("--cycles", type=float, default=N_CYCLES)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dt-scale", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eq-scale", type=float, default=200.0,
                    help="equilibration = this value x tau_fast/dt (default 200 -- local relaxation "
                         "only. * It does NOT cover whole-chain shape relaxation "
                         "(tau_chain_diff); this is for the structural smoke test")
    ap.add_argument("--drive-mode", choices=("trap", "position"), default="trap",
                    help="trap = move the trap centre (experimental, with compliance) / "
                         "position = force the driven bead's y directly (strain control, the "
                         "rheological standard). position yields K_transfer, **a different "
                         "quantity from trap's K'** -- do not compare them directly")
    ap.add_argument("--kt-scale", type=float, default=1.0,
                    help="trap stiffness multiplier. * At the default 1 the driven bead follows the "
                         "trap only 3.5-29%%. At 100 it is 82%% (at almost no cost). The design "
                         "table is in the build_ledger docstring")
    ap.add_argument("--jkr", action="store_true",
                    help="* the control: add JKR tangential bending stiffness (kappa_theta=EI/ell) on "
                         "top of DLVO. Geometry, traps, DLVO and the seed stay identical and only "
                         "the bending term is switched on, for a direct comparison")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/chain-bend-2d-dlvo/system.yaml")
    ns = [args.n] if args.n is not None else sys_["n_list"]
    omega = args.omega if args.omega is not None else sys_["omega_range"][0]
    amp = args.amp if args.amp is not None else math.sqrt(
        sys_["amp_range"][0] * sys_["amp_range"][1])

    if not (args.report or args.spec or args.run):
        print("choose what to do -- `--report`, `--spec` or `--run`")
        return 3

    rc = 0
    for i, n in enumerate(ns):
        if i:
            print()
        rc |= emit(sys_, n, omega, amp, args)
    return rc



# ════════════════════════════════════════════════════════════════════════
# L4 -- the HOOMD builder
#
# ** The WCA core must NOT use the particle diameter (sigma=d) -- WCA is repulsive
#    across all of r < 2^(1/6)*sigma, so with sigma=d the repulsion reaches
#    r=1.122d (= h=180nm) and completely tramples our secondary minimum
#    (h=11.16nm << 180nm). Instead sigma_c = d*2^(-1/6), so the WCA cutoff ends
#    **exactly at r=d** -- then WCA repels only for r<d (particle overlap) and is
#    exactly 0 for r>=d (surface gap h>=0, where DLVO lives). No overlap and no
#    double counting.
# ════════════════════════════════════════════════════════════════════════
SIGMA_CORE_STAR = 2 ** (-1.0 / 6.0)     # so the WCA cutoff ends at r*=1 (=d, surface contact)


def build_table_arrays(P: dict, r_min_star: float, r_cut_star: float, nbins: int = 8000):
    """The (U, F) arrays for `md.pair.Table`. r* = 1+h*, h*>=0.

    endpoint=False (bd-hoomd trap 10).
    """
    r = np.linspace(r_min_star, r_cut_star, nbins, endpoint=False)
    h = r - 1.0
    h = np.maximum(h, 1e-6)             # r_min_star is >= 1+eps, so this is effectively moot
    U = U_star(h, P)
    F = F_h_star(h, P)                  # increasing r* = increasing h* (same direction), so use it as-is
    U_cut = float(U_star(max(r_cut_star - 1.0, 1e-6), P))
    U = U - U_cut                       # zero at the cutoff (the bd-hoomd snippet convention)
    return r, U, F


def make_frame(n: int, ell_star: float, trapped: list[int], box_star: float,
               skip_trap: int | None = None):
    """skip_trap: do not attach a trap bond to this bead (the driven bead in position
    mode).
    """
    import gsd.hoomd
    pos = [[(i - (n - 1) / 2) * ell_star, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:
        pos.append(list(pos[g]))
        typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.particles.mass = [1.0] * len(pos)
    f.configuration.box = [box_star] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    # Trap bonds only. No backbone bond (that is the hypothesis). * In position
    # driving the driven bead is not integrated, so a trap on it is physically
    # meaningless and the stretched bond contaminates pe_per_particle.
    grp = [[g, n + j] for j, g in enumerate(trapped) if g != skip_trap]
    f.bonds.N = len(grp)
    f.bonds.types = ["trap"]
    f.bonds.typeid = [0] * len(grp)
    f.bonds.group = np.array(grp)
    return f


def make_bending_force(A, n_beads):
    """Implement the linearized bending U = (1/2)*kappa_theta*sum(theta_i^2) directly
    with `md.force.Custom` (F_y = -A y).

    ** Why not `md.angle.Harmonic`: it clamps sin theta at SMALL=1.414e-3, which
       shrinks **only the force** on a nearly straight chain (the energy is 0.000%
       accurate, so an energy check does not catch it). This was established by
       measurement in chain-bend-2d-oscill (bd-hoomd trap 15), and this system is in
       the same regime. The implementation verified there is carried over unchanged --
       the model (bending_matrix) and the implementation match exactly.
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
                tags = np.array(snap.particles.tag, copy=True)   # * tag indexing is mandatory
                pos = np.array(snap.particles.position, copy=True)
                m = tags < self.n                                # exclude the ghosts
                y = np.zeros(self.n)
                y[tags[m]] = pos[m, 1]
                fy = -(self.A @ y)
                arr.force[:] = 0.0
                arr.potential_energy[:] = 0.0
                arr.force[m, 1] = fy[tags[m]]
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


UPDATE_EVERY = 50


@RUN.builder("chain-bend-2d-dlvo")
def build(spec, outdir=None) -> RUN.Build:
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_beads"])
    trapped = sorted(int(t) for t in P["trapped"])
    k_t = float(P["k_t_star"])
    amp, omega = float(P["amp_star"]), float(P["omega_star"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    mid = trapped[len(trapped) // 2]
    ghost_mid = n + trapped.index(mid)

    reduced = {"kappa_star": P["kappa_star"], "edl_amp": P["edl_amp"],
              "vdw_amp": P["vdw_amp"], "a_star": P["a_star"]}
    h_min_star = float(P["h_min_star"])
    ell_star = 1.0 + h_min_star                      # natural bond length (centre to centre, r*=1+h_min)
    r_cut_star = 1.0 + float(P["cutoff_h_star"])
    r_min_star = 1.0 + 1e-6
    box_star = 4.0 * (n - 1) * ell_star

    _pos_drive = str(P.get("drive_mode", "trap")) == "position"
    sim = SIM.make_sim(make_frame(n, ell_star, trapped, box_star,
                                  skip_trap=mid if _pos_drive else None), seed=seed)

    cell = md.nlist.Cell(buffer=0.2)
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_star, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_star, U=U_arr, F=F_arr)
    tab.params[("A", "G")] = dict(r_min=r_min_star, U=U_arr * 0, F=F_arr * 0)
    tab.params[("G", "G")] = dict(r_min=r_min_star, U=U_arr * 0, F=F_arr * 0)
    tab.r_cut[("A", "G")] = tab.r_cut[("G", "G")] = r_min_star     # effectively off

    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
    wca.params[("A", "G")] = wca.params[("G", "G")] = dict(epsilon=0.0, sigma=SIGMA_CORE_STAR)

    bond = md.bond.Harmonic()
    bond.params["trap"] = dict(k=k_t, r0=0.0)          # the trap = a harmonic bond to a ghost. No backbone

    # * The control -- JKR bending stiffness. `angle.Harmonic` is **not** used: on a
    #   nearly straight chain the sin-theta clamp makes the force silently wrong
    #   (bd-hoomd trap 15, established by measurement in chain-bend-2d-oscill). The
    #   **same** linearized bending is applied directly via force.Custom.
    forces = [tab, wca, bond]
    bend = None
    if bool(P.get("jkr", False)):
        kth_star = float(P["kappa_theta_star"])
        bend = make_bending_force(bending_matrix(n, kth_star, ell_star), n)
        forces.append(bend)

    # ** Two driving modes (rheologically different quantities)
    #   "trap"     -- move the trap centre as y=a*sin(omega*t) and let the bead be
    #                pulled along. Closer to the experiment (optical tweezers), but
    #                **trap compliance** means the bead does not follow the commanded
    #                position (measured tracking: JKR 3.5% / DLVO 29%).
    #                -> K* = k_t(y_hat_c/y_hat - 1)  [the stiffness the driving trap feels]
    #   "position" -- force the driven bead's y **directly** (strain control). There is
    #                no compliance at all, so the deformation is imposed exactly, and it
    #                matches the standard rheological protocol (impose a strain, measure
    #                a stress). [P1][P2] also, experimentally, move the centre and use
    #                **the force on the end beads** as the sensor.
    #                -> K*_transfer = k_t_sensor*y_hat_end / y_hat_mid  [transfer stiffness]
    #   WARNING: the two K* are **not the same quantity** -- the first is a
    #      drive-point stiffness, the second a transfer function. Do not compare their
    #      values directly; read each as a DLVO vs JKR contrast.
    pos_drive = str(P.get("drive_mode", "trap")) == "position"
    if pos_drive:
        # exclude the driven bead from integration -> no Brownian motion, and the
        # updater writes its position directly
        bd_filter = hoomd.filter.SetDifference(hoomd.filter.Type(["A"]),
                                               hoomd.filter.Tags([mid]))
    else:
        bd_filter = hoomd.filter.Type(["A"])
    bd = md.methods.Brownian(filter=bd_filter, kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    # position mode moves **the bead itself**; trap mode moves **the ghost**
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=_move_ghost_action(mid if pos_drive else ghost_mid, amp, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, int(Nm["n_prod"]) // 200))

    def pe_per_particle():
        return float(sum(np.array(f.energies).sum() for f in forces)) / n

    def sample(timestep, phase):
        snap = sim.state.get_snapshot()
        p = np.array(snap.particles.position, dtype=float)
        nb = p[:n, :2]
        sep = np.linalg.norm(nb[:, None, :] - nb[None, :, :], axis=-1)
        np.fill_diagonal(sep, np.inf)
        nn_sep = np.array([sep[i, i + 1] for i in range(n - 1)])   # NN spacing (bond monitoring)
        return {"t": timestep * dt, "y_bead": float(p[mid, 1]),
                "y_ghost": float(p[ghost_mid, 1]),
                "y_end0": float(p[trapped[0], 1]), "y_end1": float(p[trapped[-1], 1]),
                "min_sep_all": float(sep.min()), "nn_sep_max": float(nn_sep.max()),
                "nn_sep_min": float(nn_sep.min()),
                # ** bow -- the maximum deviation from the chord joining the two ends.
                #   **Pure bending deformation, with rigid translation and rotation
                #   removed.** Why it is this case's central discriminant: measured over
                #   the full y range, DLVO and JKR differ by 3.4x; measured as bow, by
                #   19.6x -- because most of what a JKR chain does is translate as a
                #   whole (the traps have finite stiffness), not bend.
                #   shape_localization failed to distinguish the two systems at all
                #   (1.481 vs 1.470), being dominated by thermal noise. Bow replaces it.
                "bow": float(_bow(p[:n, 1])),
                "shape_y": p[:n, 1].copy()}       # * the full y profile -- for judging smooth curvature vs local kinking

    def finalize(cols):
        t, yb, yg = cols["t"], cols["y_bead"], cols["y_ghost"]
        blocks_y = LI.lockin_blocks(t, yb, omega, n_blocks=min(10, max(2, len(t) // 20)))
        blocks_g = LI.lockin_blocks(t, yg, omega, n_blocks=min(10, max(2, len(t) // 20)))
        Kb = np.array([LI.k_star(a, b, k_t, omega) for a, b in zip(blocks_y, blocks_g)])
        K, Ksem = LI.agg(Kb)
        yh, _ = LI.agg(blocks_y)
        gh, _ = LI.agg(blocks_g)

        # ** Under position driving, emit **the transfer stiffness** separately.
        #   The deformation is imposed exactly (zero compliance), so what corresponds to
        #   a stress is the force on the end sensor beads, F_end = k_t*y_end
        #   ([P1][P2]'s experimental protocol).
        #       K_transfer = k_t·⟨ŷ_end⟩ / ŷ_mid
        #   WARNING: **a different quantity** from trap mode's K' (a drive-point
        #      stiffness). Do not compare their values directly.
        pos_mode = str(P.get("drive_mode", "trap")) == "position"
        K_tr = None
        if pos_mode:
            ye = 0.5 * (np.asarray(cols["y_end0"], float) + np.asarray(cols["y_end1"], float))
            be = LI.lockin_blocks(t, ye - ye.mean(), omega,
                                  n_blocks=min(10, max(2, len(t) // 20)))
            eh, esem = LI.agg(be)
            mh, _ = LI.agg(blocks_y)          # the driven bead = the imposed deformation (measured phasor)
            K_tr = (k_t * eh / mh) if abs(mh) > 0 else complex("nan")

        # * The prediction differs per branch -- and it is fixed before the result is seen.
        if bool(P.get("jkr", False)):
            # the control: with the bending term on, the static stiffness the driving
            # trap feels is the prediction.
            #   (A_bend + T) y = k_t y_c e_mid  →  K = k_t(y_c/y_mid − 1)
            #   the same expression as chain-bend-2d-oscill.driven_static_stiffness.
            A_ = bending_matrix(n, float(P["kappa_theta_star"]), ell_star)
            for i in trapped:
                A_[i, i] += k_t
            e_ = np.zeros(n)
            e_[mid] = k_t
            K_pred = float(k_t * (1.0 / np.linalg.solve(A_, e_)[mid] - 1.0))
            K_derivation = ("with the bending term present, the static limit of linear "
                            "response gives the prediction: solve "
                            "(A_bend+T)y = k_t*y_c*e_mid for K = k_t(y_c/y_mid - 1). "
                            "It is not 48EI/L^3, because the ends are finite-stiffness "
                            "traps rather than rigid clamps. Valid only at low frequency "
                            "(quasi-static) -- at large De, viscosity mixes in.")
            K_note = (f"* the control (JKR bending ON). Static-limit prediction {K_pred:.4g} "
                      f"kT/d^2. Compared against the DLVO-only branch (prediction 0) at "
                      f"identical geometry and identical seed")
        else:
            K_pred = 0.0
            K_derivation = ("for a central force U(r) with the bond at its natural length "
                            "ell (U'(ell)=0), the bond energy against transverse "
                            "displacements y_i, y_{i+1} is "
                            "U(ell)+U'(ell)(y_i-y_{i+1})^2/(2*ell)+O(y^4) -- and since "
                            "U'(ell)=0 the O(y^2) term vanishes. It is a local symmetry "
                            "argument independent of the bead count and the bond "
                            "topology, so it carries over unchanged to the full "
                            "combination (the chain). Limit conditions: small y, and the "
                            "traps neither stretching nor compressing the chain.")
            K_note = ("* prediction 0 -- G1 (a straight chain with pure central forces at "
                      "natural-length equilibrium has exactly zero linear bending "
                      "stiffness). A significantly non-zero value means tension or "
                      "finite-deformation effects have mixed in")

        shapes = np.stack(cols["shape_y"])          # (n_samples, n)
        # ** The G1 check: smooth curvature spreads the profile's second derivative
        #   gently, whereas "trianglelike" buckling should make the angle (the second
        #   derivative) spike locally at one particular bond.
        d2 = np.diff(shapes, n=2, axis=1)            # (n_samples, n-2)
        d2_mean = np.abs(d2).mean(axis=0)
        kurt_like = float(d2_mean.max() / (d2_mean.mean() + 1e-30))   # near 1 = smooth, large = locally concentrated

        # ** bow -- this case's central discriminant (see the _bow docstring).
        #   Two versions are emitted separately:
        #     bow_rms    the time average -- **thermal fluctuation AND drive** mixed
        #     bow_drive  only the drive-frequency component, extracted by lock-in --
        #                **bending due to the drive alone**
        #   The latter is the real response. Reading only the former makes a thermally
        #   floppy chain (DLVO) look like it "bends a lot", which is not an elastic
        #   response to the drive.
        bow = np.asarray(cols["bow"], dtype=float)
        bow_rms = float(np.sqrt((bow ** 2).mean()))
        bb = LI.lockin_blocks(t, bow - bow.mean(), omega,
                              n_blocks=min(10, max(2, len(t) // 20)))
        bow_hat, bow_sem = LI.agg(bb)
        bow_drive = float(abs(bow_hat))

        # ** Under position driving, K' (the driving-trap stiffness) is **undefined** --
        #   there is no trap. Computing it anyway gives y_hat_c=0 (the ghost does not
        #   move) and therefore K' = k_t(0/y_hat - 1) = -k_t, **a meaningless constant**
        #   (measured -5217.11 = -k_t*, with K_sem=1.2e-12, so even the spread is 0).
        #   The number looks plausible enough to be misread as a real measurement, so it
        #   is not emitted at all. The observable in position mode is K_transfer_*.
        obs = ([] if pos_mode else [
            MET.observable("K_prime", float(K.real), K_pred, "kT/d^2",
                           role="implementation_check", sigma=Ksem, tol_sigma=3.0,
                           derivation=K_derivation,
                           note=K_note),
            MET.observable("K_doubleprime", float(K.imag), None, "kT/d^2", role="measurement"),
            MET.observable("K_sem", Ksem, None, "kT/d^2", role="measurement"),
        ]) + [
            MET.observable("y_response", float(abs(yh)), None, "d", role="measurement"),
            *([MET.observable("K_transfer_prime", float(K_tr.real), None, "kT/d^2",
                              role="measurement",
                              note="* the real part of the transfer stiffness under position (strain) "
                                   "control = k_t*y_hat_end/y_hat_mid. The deformation is imposed "
                                   "exactly, so there is no trap compliance (the same "
                                   "configuration as [P1][P2]'s experimental protocol: move the "
                                   "centre, measure the force at the end sensors). WARNING: a "
                                   "different quantity from trap mode's K'"),
               MET.observable("K_transfer_dprime", float(K_tr.imag), None, "kT/d^2",
                              role="measurement", note="imaginary part of the transfer stiffness")]
              if K_tr is not None else []),
            MET.observable("bow_rms", bow_rms, None, "d", role="measurement",
                           note="** bow RMS -- the maximum deviation from the chord joining the "
                                "two ends. Pure bending deformation with rigid translation "
                                "removed (the traps have finite stiffness so the chain gets "
                                "pushed as a whole, and measuring raw y lets that mask the "
                                "bending signal). Measured discriminating power: 3.4x over the "
                                "full y range vs 19.6x for bow (DLVO vs JKR). "
                                "* But this value mixes thermal fluctuation with the drive -- "
                                "for the drive response alone use bow_drive below"),
            MET.observable("bow_drive", bow_drive, None, "d", role="measurement",
                           sigma=bow_sem,
                           note="* bow with only the drive-frequency component extracted by "
                                "lock-in = **the bending response due to the drive.** Reading "
                                "bow_rms alone makes a thermally floppy chain look like it "
                                "'bends a lot', which is not an elastic response"),
            MET.observable("shape_localization", kurt_like, None, "dimensionless",
                           role="measurement",
                           note="WARNING: **measured to have no discriminating power** (JKR 1.481 "
                                "vs DLVO 1.470 -- obviously different by eye, identical by this "
                                "metric), because thermal noise dominates the second derivative "
                                "in both. Kept but must not be used for a verdict -- bow_drive "
                                "replaces it. Defined as max|theta''|/mean|theta''|"),
        ]
        return {"observables": obs,
                "extra": {**({} if pos_mode else
                              {"K_prime": float(K.real), "K_doubleprime": float(K.imag),
                               "K_sem": Ksem}),
                          "y_resp_abs": float(abs(yh)),
                          "drive_abs": float(abs(gh)), "omega_star": omega,
                          **({"K_transfer_prime": float(K_tr.real),
                              "K_transfer_dprime": float(K_tr.imag)} if K_tr is not None else {}),
                          "drive_mode": str(P.get("drive_mode", "trap")),
                          "bow_rms": bow_rms, "bow_drive": bow_drive,
                          "bow_drive_sem": float(bow_sem),
                          "bow_max": float(bow.max()),
                          "shape_localization": kurt_like,
                          "shape_profile_mean": shapes.mean(axis=0).tolist(),
                          "min_sep_all": float(cols["min_sep_all"].min()),
                          "nn_sep_max": float(cols["nn_sep_max"].max()),
                          "nn_sep_min": float(cols["nn_sep_min"].min())}}

    every = max(1, int(Nm["sample_every"]))
    return RUN.Build(
        sim=sim, forces=forces, n_particles=n,
        sample=sample, pe_per_particle=pe_per_particle, sample_every=every,
        phases=[RUN.Phase("warm-up", int(Nm["n_eq"]), collect=False,
                          note="drive ON, local (bond) relaxation only -- not whole-chain shape relaxation (smoke)"),
                RUN.Phase("production", int(Nm["n_prod"]), every,
                          note=f"ω*={omega:.4g}")],
        tags=["2D", "chain", "dlvo", "pair_table", "oscillatory_drive", "no_bending",
              "hypothesis_test", "structural"],
        physical={"n_beads": n, "omega_star": omega, "amp_star": amp, "k_t_star": k_t,
                  "well_depth_star": float(P["well_depth_star"]),
                  "k_bond_star": float(P["k_bond_star"])},
        finalize=finalize)


if __name__ == "__main__":
    sys.exit(main())
