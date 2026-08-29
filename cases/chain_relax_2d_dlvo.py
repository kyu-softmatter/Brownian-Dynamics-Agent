"""`chain-relax-2d-dlvo` -- the undriven (free-relaxation) partner of
chain-bend-2d-dlvo.

Same physics (PMMA beads, a DLVO secondary-minimum central-force bond, no explicit
bending and no friction), but with **the trap and the oscillation removed
entirely** -- the simplest undriven form of "a chain with attraction only and no
inter-particle friction" (CLAUDE.md rule 8: build the static system first, add
motion afterwards). chain-bend-2d-dlvo skipped that static stage and went straight
to the oscillation experiment; this case fills it in.

Two experiments (--init):
  straight  thermalize a straight chain and measure the **radial (bond-stretch)**
            thermal fluctuation. An equipartition golden test of the same kind as
            the trap's <x^2>=kT/k -- it confirms that the DLVO table potential was
            ported into this case correctly.
            WARNING: the prediction is NOT kT*/k_bond*. The harmonic approximation
            underestimated the true value by 4.6x, because the well is asymmetric
            and softer on the outside. See bond_variance_boltzmann below.

  kink      give exactly one central bond a precise turn angle dPhi (every other
            bond perfectly straight, every bond starting exactly at its natural
            length -- pure bending perturbed with no stretch signal), release it,
            and see whether the bow recovers elastically or disperses diffusively.
            G1 (the transverse linear bending stiffness is exactly zero) had only
            been derived algebraically and never confirmed directly without driving.
            There is no prior quantitative prediction, so this is a measurement.

The DLVO expressions and the table-potential implementation are **not redefined**;
they are imported from chain_bend_dlvo_2d as-is (the second case to do so,
following the precedent `network` already set). Whether to promote them into
`bdbot/` is decided when a third case appears -- the "has it appeared twice"
principle in CLAUDE.md.
    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/chain_relax_2d_dlvo.py --init straight --report
    $PY cases/chain_relax_2d_dlvo.py --init kink --kink-angle 0.3 --report
    $PY cases/chain_relax_2d_dlvo.py --init straight --smoke --run     # a quick sanity check
    $PY cases/chain_relax_2d_dlvo.py --init kink --run
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bdbot import checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, scales as SC, sim as SIM, stats as ST  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from bdbot.units import Q  # noqa: E402

# * The DLVO expressions are used **as-is** from chain-bend-2d-dlvo, where they were
#   verified in SI -- they are not written twice (the same precedent as `network`).
from chain_bend_dlvo_2d import (  # noqa: E402
    CUTOFF_H_STAR, SIGMA_CORE_STAR, build_table_arrays, dlvo_reduced_params, find_well,
    F_h_star, U_star,
)

ROOT = Path(__file__).resolve().parent.parent


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
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# 2. geometry -- shape descriptors invariant under rotation and translation
#    ** chain-bend-2d-dlvo's `_bow(y)` used lab-frame y directly, because the trap
#    fixed the orientation. This case has no trap, so the whole chain rotates and
#    translates freely -- the descriptor MUST be rotated into the chain's own body
#    frame (the axis joining the two ends). Otherwise "bow" is really just measuring
#    thermal rotational drift.
# ════════════════════════════════════════════════════════════════════════
def bond_vectors(pos: np.ndarray) -> np.ndarray:
    return pos[1:] - pos[:-1]


def bend_angles(pos: np.ndarray) -> np.ndarray:
    """The local turn angle dtheta_i of the bond direction at each interior bead,
    (n-2,). Zero means that point is locally straight.

    Invariant under rotation and translation -- it looks only at the **difference**
    between consecutive bond directions. Equivalent at small angles to the discrete
    curvature theta_i=(y_{i+1}-2y_i+y_{i-1})/ell used by `bending_matrix()` in
    chain_bend_dlvo_2d.py (both are first-order discrete curvature) -- redefined via
    turn angles here because with no trap the y-based definition cannot be used.
    """
    bv = bond_vectors(pos)
    ang = np.arctan2(bv[:, 1], bv[:, 0])
    dtheta = np.diff(ang)
    return (dtheta + np.pi) % (2 * np.pi) - np.pi


def bow_metrics(pos: np.ndarray) -> tuple[float, float]:
    """Align the axis joining the two ends with x', then measure the deviation from
    that axis (max, rms) [d].
    """
    d = pos[-1] - pos[0]
    L_ee = float(np.hypot(d[0], d[1]))
    if L_ee < 1e-9:
        return 0.0, 0.0
    u = d / L_ee
    rel = pos - pos[0]
    yp = -rel[:, 0] * u[1] + rel[:, 1] * u[0]        # transverse component in the body frame.
                                                 # yp[0]=yp[-1]=0 by construction
    return float(np.abs(yp).max()), float(np.sqrt(np.mean(yp ** 2)))


def min_nnn_gap_star(pos: np.ndarray) -> float:
    """The minimum surface gap h*=r*-1 over pairs with |i-j|>=2.

    Detects premature non-adjacent bonding (the "trianglelike" structure).
    """
    n = len(pos)
    if n < 4:
        return float("inf")
    diff = pos[:, None, :] - pos[None, :, :]
    dist = np.hypot(diff[..., 0], diff[..., 1])
    idx = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    mask = idx >= 2
    return float(dist[mask].min()) - 1.0


def kink_positions(n: int, ell_star: float, kink_angle: float) -> np.ndarray:
    """Give exactly one central bond the turn angle `kink_angle`; every other bond is
    perfectly straight.

    Every bond starts at exactly `ell_star` -- so no radial (stretch) signal is mixed
    in and only the bending degrees of freedom are perturbed. `kink_angle=0` gives a
    straight chain (= --init straight).
    """
    mid = n // 2
    pos = np.zeros((n, 2))
    th_l, th_r = -kink_angle / 2.0, kink_angle / 2.0
    for k in range(1, n - mid):
        pos[mid + k] = pos[mid + k - 1] + ell_star * np.array([math.cos(th_r), math.sin(th_r)])
    for k in range(1, mid + 1):
        pos[mid - k] = pos[mid - k + 1] - ell_star * np.array([math.cos(th_l), math.sin(th_l)])
    pos -= pos.mean(axis=0)
    return pos


def bond_variance_boltzmann(p: dict, w: dict, cutoff_h_star: float, nbins: int = 400_000):
    """The **true** (anharmonic) thermal-equilibrium variance of the bond stretch
    (h-h_min), integrated numerically inside the secondary-minimum well
    (h>barrier_h) -- a basin-restricted Boltzmann average at kT*=1.

    ** The harmonic approximation `1/k_bond_star` sees only the local curvature at
    the bottom of the well. This well is **asymmetric** -- steep on the inside
    (barrier side) and much softer on the outside (h->inf, U->0-) -- so the harmonic
    approximation underestimates the true value. Measured (smoke run, n=9): 2.9x the
    harmonic prediction, and 4.6x for the full-well numerical integral. The same
    class of trap as soft-r3's "the Einstein cage approximation breaks under
    anharmonicity" (bd-physics 6.2). Here it is **not demoted to a regime
    indicator** -- the integral gives an exact prediction and it is used as a real
    golden test, because this system does not require the harmonic approximation
    (free relaxation, no trap; soft-r3 could not avoid it because its dt design
    needed it).

    Why the integration range is not opened past the barrier: the barrier is 416 kT,
    so it can **never** be crossed on the simulation timescale. Integrating over the
    full range (including the primary minimum at h->0) makes the distribution
    collapse entirely to h->0 because of the vdW divergence (measured). What the
    simulation actually samples is **inside the secondary-minimum well only**, so the
    integral must be restricted to that basin.
    """
    h = np.linspace(w["barrier_h"], cutoff_h_star, nbins)
    U = U_star(h, p)
    wt = np.exp(-(U - U.min()))
    Z = np.trapezoid(wt, h)
    mean_h = float(np.trapezoid(h * wt, h) / Z)
    var_h = float(np.trapezoid((h - mean_h) ** 2 * wt, h) / Z)
    return var_h, mean_h


# ════════════════════════════════════════════════════════════════════════
# 3. the scale ledger
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, n: int, kink_angle: float, *, dt_scale=1.0,
                 T_obs_tau: float, eq_scale: float = 200.0) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B, tau_p = b["kT"], b["gamma"], b["tau_B"], b["tau_p"]

    p = dlvo_reduced_params(sys_)
    w = find_well(p)
    ell = d * (1 + w["h_min"])
    k_bond = Q(w["k_bond_star"], "dimensionless") * kT / d ** 2
    sigma_bond = (kT / k_bond) ** 0.5
    L_chain = (n - 1) * ell

    # ** The bond stretch is this system's ONLY stiffness mode (the bending stiffness
    #   is structurally 0 and there is no trap) -- unlike chain-bend-2d-dlvo, there is
    #   no need to pick the fastest among several candidates.
    tau_bond = C.relaxation_time(gamma, k_bond)
    dt = dt_scale * C.dt_from_gate(tau_bond)
    T_obs = Q(T_obs_tau, "dimensionless") * tau_bond

    # Safety of the kink initial condition -- every bond starts at exactly ell so the
    # "stretch" is 0, but the kink itself can narrow the next-nearest-neighbour
    # (two-bonds-away) gap. Checked directly from the geometry.
    ell_star = float((ell / d).to("dimensionless").magnitude)
    nnn_gap0 = min_nnn_gap_star(kink_positions(n, ell_star, kink_angle)) if kink_angle else float("inf")

    # How far past h_min a bond has to stretch to reach the maximum tensile force
    # (F_max) -- an informational anchor. (The kink was built to produce no stretch
    # signal so it is not actually used here, but it stays in the ledger as a
    # reference scale for dt and safety discussions.)
    hs = np.linspace(w["h_min"], CUTOFF_H_STAR, 20_000)
    Fs = F_h_star(hs, p)
    F_max = float(-Fs.min())

    # ** The golden test's real predicted value -- not the harmonic approximation
    #   (1/k_bond_star) but the anharmonic Boltzmann integral over the whole well
    #   (basin-restricted; see bond_variance_boltzmann above).
    dl_var_boltz, dl_mean_boltz_h = bond_variance_boltzmann(p, w, CUTOFF_H_STAR)

    lg = SC.ScaleLedger()
    lg.add_length("sigma_bond", sigma_bond.to("m"), "bond-stretch thermal width sqrt(kT/k_bond)", star=True)
    lg.add_length("h_min", Q(w["h_min"], "dimensionless") * d, "secondary-minimum position (surface gap)")
    lg.add_length("d", d, "bead diameter")
    lg.add_length("ell", ell.to("m"), "natural bond length (centre to centre, d+h_min)")
    lg.add_length("L_chain", L_chain.to("m"), "chain contour length (n-1)*ell")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_bond", tau_bond,
               "** gamma/k_bond bond stretch -- this system's **only** stiffness mode. It sets dt",
               star=True)
    lg.add_time("tau_B", tau_B, "d^2/D_t diffusion (reference)")
    lg.add_time("T_obs", T_obs, "observation window (multiples of tau_bond -- only the "
               "local mode is needed; whole-chain shape relaxation "
               "tau_chain_diff is not)", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("k_bond_d2", (k_bond * d ** 2).to("J"), "k_bond*d^2 bond-stretch stiffness", star=True)
    lg.add_energy("well_depth", Q(-w["U_min"], "dimensionless") * kT, "|secondary-minimum depth|")
    lg.declare_absent(
        "box",
        "no periodic boundaries (one chain, and with no trap it relaxes and rotates "
        "freely). HOOMD's frame needs a formal box, so it is set large -- 4x the "
        "chain's own maximum extent, the same margin as chain-bend-2d-dlvo -- purely "
        "to avoid interacting with its own periodic image. It is not a physically "
        "meaningful confinement scale.")
    lg.declare_absent(
        "bending_stiffness",
        "** structurally absent (G1, derived in chain-bend-2d-dlvo). Confirming that "
        "structural fact by execution, in the minimal undriven configuration, is the "
        "whole point of this case. No invented substitute scale is put here.")
    lg.derived = dict(gamma=gamma, kT=kT, d=d, tau_B=tau_B, ell=ell.to("m"),
                      L_chain=L_chain.to("m"), ell_star=ell_star, k_bond=k_bond,
                      k_bond_star=w["k_bond_star"], sigma_bond=sigma_bond.to("m"),
                      tau_bond=tau_bond, dt=dt, T_obs=T_obs, reduced=p,
                      h_min_star=w["h_min"], well_star=w["U_min"], barrier_star=w["barrier_U"],
                      F_max_star=F_max, nnn_gap0_star=nnn_gap0, kink_angle=kink_angle, n=n,
                      dl_var_boltz=dl_var_boltz, dl_mean_boltz_h=dl_mean_boltz_h)
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " * This system has neither a bending stiffness nor a "
        "trap -- the only stiffness mode is the bond stretch (tau_bond), which makes "
        "the dt and equilibration verdicts simpler than in chain-bend-2d-dlvo.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 4. dimensionless groups + separation checks
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(lg, n, init, kink_angle):
    D = lg.derived
    r = lg.ratio

    groups = [
        ND.Group("k_bond_star", D["k_bond_star"], ("energies", "k_bond_d2"),
                 ("energies", "kT"), "k_bond d^2/kT", "bond-stretch stiffness -- the only harmonic mode"),
        ND.Group("well_depth/kT", r("energies", "well_depth", "kT"), ("energies", "well_depth"),
                 ("energies", "kT"), "", "bond depth -- reversible if it is on the scale of kT"),
        ND.Group("sigma_bond/d", r("lengths", "sigma_bond", "d"), ("lengths", "sigma_bond"),
                 ("lengths", "d"), "", "bond-stretch thermal fluctuation (the predicted golden test)"),
        ND.Group("n_beads", float(n), None, None, "", "chain length (input)"),
        ND.Group("kink_angle_rad", float(kink_angle), None, None, "",
                 "initial kink angle (input; 0 for straight)"),
        ND.Group("dt/tau_bond", r("times", "dt", "tau_bond"), ("times", "dt"),
                 ("times", "tau_bond"), "", "integration resolution -- the only stiffness mode"),
        ND.Group("T_obs/tau_bond", r("times", "T_obs", "tau_bond"), ("times", "T_obs"),
                 ("times", "tau_bond"), "", "observation window (for local-mode statistics)"),
        ND.Group("St", r("times", "tau_p", "tau_bond"), ("times", "tau_p"),
                 ("times", "tau_bond"), "tau_p/tau_bond", "inertia vs bond stretch"),
    ]
    checks = [
        C.Check("model", "note: tau_p/tau_bond", r("times", "tau_p", "tau_bond"), C.GATE, "<=",
              "the same unverified state as chain-bend-2d-dlvo (same particles, same "
              "bond) -- not re-verified here. When that case verifies it against "
              "OverdampedViscous, this case is verified with it", hard=False),
        C.Check("integration", "bond stretch resolved dt/tau_bond", r("times", "dt", "tau_bond"), C.GATE, "<=",
              "this system's only stiffness mode. Miss it and it diverges"),
        C.Check("statistics", "observation sufficient T_obs/tau_bond", r("times", "T_obs", "tau_bond"),
              1000.0, ">=", "the minimum multiple needed for local (bond) "
              "equipartition statistics -- whole-chain shape relaxation "
              "(tau_chain_diff) is not required", hard=False),
    ]
    if kink_angle:
        checks.append(C.Check(
            "geometry", "kink NNN gap (initial) vs cutoff", D["nnn_gap0_star"], CUTOFF_H_STAR, ">=",
            f"if the two-bonds-away bead is already inside the DLVO cutoff at the "
            f"moment of release (premature non-adjacent bonding), the design intent of "
            f"'perturbing pure bending only' is broken", hard=False))
    return groups, checks


def report_blocks(sys_, lg, n, init, kink_angle, n_steps):
    D = lg.derived
    inp = [R.kv("d", f"{sys_['d'].value:~.4gP}", sys_["d"].tier, sys_["d"].source[:44]),
           R.kv("psi0", f"{sys_['psi0'].value:~.4gP}", sys_["psi0"].tier, sys_["psi0"].source[:44]),
           R.kv("I(MgCl2)", f"{sys_['ionic_strength'].value:~.4gP}",
                sys_["ionic_strength"].tier, sys_["ionic_strength"].source[:44]),
           R.kv("n", f"{n}", 3, "matched to the value inherited from chain-bend-2d-dlvo"),
           R.kv("init", init, 3, "straight = thermalize a straight chain / kink = release and relax"),
           R.kv("kink_angle", f"{kink_angle:.3f} rad", 3, "*proposed (observation.yaml R1)")]
    der = [
        f"  bond: secondary minimum {D['well_star']:.3f} kT @ h_min*={D['h_min_star']:.5f}"
        f"   barrier {D['barrier_star']:.1f} kT",
        f"  k_bond = {D['k_bond'].to('pN/um'):~.4fP} = {D['k_bond_star']:.4e} kT/d²"
        f"   sigma_bond (harmonic) = {D['sigma_bond'].to('nm'):~.4fP}",
        f"  * bond-stretch golden-test prediction (anharmonic Boltzmann integral) = {D['dl_var_boltz']:.4e} d^2"
        f"   ({D['dl_var_boltz']*D['k_bond_star']:.2f}x the harmonic approximation -- because "
        f"the well is asymmetric and softer on the outside; the same class of "
        f"anharmonic correction as soft-r3's Einstein cage)",
        f"  ell (natural length) = {D['ell'].to('nm'):~.2fP}   L_chain = {D['L_chain'].to('um'):~.3fP}",
        f"  F_max (the well's maximum tensile force) = {D['F_max_star']:.1f} kT/d -- for "
        f"reference, since the kink produces no stretch signal (it starts exactly at ell)",
        f"  ** linear bending stiffness = 0 (structural, declare_absent) -- confirming this directly is what this run is for",
    ]
    if kink_angle:
        der.append(f"  kink NNN gap (initial, h*) = {D['nnn_gap0_star']:+.4f}"
                   f"   (cutoff {CUTOFF_H_STAR:.2f} -- larger means no premature non-adjacent bonding)")
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}  = {lg.ratio('times','dt','tau_B'):.3e} τ_B",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}  = {lg.ratio('times','T_obs','tau_bond'):.3g} τ_bond",
        f"  steps   = {n_steps:,}   × n={n}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# main — L3
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, n, init, kink_angle, args):
    T_obs_tau = args.tobs if args.tobs is not None else (3.0e3 if init == "kink" else 2.0e4)
    lg = build_ledger(sys_, n, kink_angle, dt_scale=args.dt_scale, T_obs_tau=T_obs_tau,
                      eq_scale=args.eq_scale)
    D = lg.derived
    dt = lg.get("times", "dt")

    if init == "kink":
        n_eq = 0
        n_release = int(round(float((D["T_obs"] / dt).to(""))))
        sample_every = max(1, n_release // args.samples)
        n_release = (n_release // sample_every) * sample_every
        n_prod = n_release
    else:
        n_eq = int(round(args.eq_scale * float((D["tau_bond"] / dt).to(""))))
        n_prod = int(round(float((D["T_obs"] / dt).to(""))))
        sample_every = max(1, n_prod // args.samples)
        n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(lg, n, init, kink_angle)
    tag = f"n{n}-{init}"
    if init == "kink":
        tag += f"-a{kink_angle:.3f}"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"

    p = D["reduced"]
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_beads": n, "init": init, "kink_angle": kink_angle,
                "kappa_star": p["kappa_star"], "edl_amp": p["edl_amp"], "vdw_amp": p["vdw_amp"],
                "a_star": p["a_star"], "cutoff_h_star": CUTOFF_H_STAR,
                "h_min_star": D["h_min_star"], "k_bond_star": D["k_bond_star"],
                "well_star": D["well_star"],
                # * The true (anharmonic) golden-test prediction, integrated once at
                #   L3 -- L4 reads this number as-is (the spec is the only contract; it
                #   is never re-derived).
                "dl_var_boltz": D["dl_var_boltz"]},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": args.seed},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_eq, n_prod


def emit(sys_, n, init, kink_angle, args) -> int:
    lg, spec, groups, checks, n_eq, n_prod = build_spec(sys_, n, init, kink_angle, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n, init, kink_angle, n_eq + n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  n={n} init={init}"
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

    # ** result.txt -- the completion marker. Without it, `status` counts the run as
    #   zero and an "incomplete cleanup" pass deletes a completed run. This is the
    #   case script's responsibility, not the engine's.
    if v["status"] != "skipped":
        obs_lines = []
        try:
            mj = json.loads((outdir / "metrics.json").read_text())
            for o in mj.get("observables", []):
                m = o.get("measured")
                p_ = o.get("predicted")
                tail = f"   (predicted {p_:.6g})" if isinstance(p_, (int, float)) else ""
                obs_lines.append(f"  {o['name']:<28} {m:.6g}{tail}" if m is not None
                                 else f"  {o['name']:<28} —")
        except Exception as e:
            obs_lines.append(f"  (could not read metrics.json: {e})")
        result = "\n".join(["=" * 84, f"RESULT -- {run_id}", "=" * 84,
                            *obs_lines, "=" * 84, verdict_txt])
        (outdir / "result.txt").write_text(report + "\n" + result)
        make_plots(sys_, lg, n, init, kink_angle, outdir)
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


# ════════════════════════════════════════════════════════════════════════
# L4 -- the HOOMD builder. With no trap and no driving this is far simpler than
# chain-bend-2d-dlvo: no ghost particle, no CustomUpdater, no force.Custom. Just N
# real beads + Table (DLVO) + WCA (core) + Brownian over all particles.
# ════════════════════════════════════════════════════════════════════════
@RUN.builder("chain-relax-2d-dlvo")
def build(spec, outdir=None) -> RUN.Build:
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_beads"])
    init, kink_angle = str(P["init"]), float(P["kink_angle"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])

    reduced = {"kappa_star": P["kappa_star"], "edl_amp": P["edl_amp"],
              "vdw_amp": P["vdw_amp"], "a_star": P["a_star"]}
    h_min_star = float(P["h_min_star"])
    ell_star = 1.0 + h_min_star
    r_cut_star = 1.0 + float(P["cutoff_h_star"])
    r_min_star = 1.0 + 1e-6
    box_star = 4.0 * max(1, n - 1) * ell_star           # a generous box -- purely to avoid the self periodic image

    pos0 = kink_positions(n, ell_star, kink_angle)
    sim = SIM.make_sim(SIM.frame_2d(pos0, box_star), seed=seed)

    cell = md.nlist.Cell(buffer=0.2)
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_star, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_star, U=U_arr, F=F_arr)

    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)

    integ, bd = SIM.attach_brownian(sim, dt, [tab, wca])           # no friction term -- plain BD
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, n_prod // 200))
    L = box_star

    def unwrapped_xy():
        # * `sim.state.get_snapshot()` already gathers in tag order (unlike
        #   force.Custom's cpu_local_snapshot) -- no tag re-indexing is needed (the
        #   same pattern as xy() in soft_r3_2d.py). Tag indexing is for local
        #   snapshots only (see the bd-hoomd trap).
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, dtype=float)[:, :2]
        img = np.array(snap.particles.image, dtype=float)[:, :2]
        return pos + img * L

    def pe_pp():
        return float(np.array(tab.energies).sum() + np.array(wca.energies).sum()) / n

    def sample(timestep, phase):
        xy = unwrapped_xy()
        bow_max, bow_rms = bow_metrics(xy)
        dtheta = bend_angles(xy)
        bl = np.hypot(*bond_vectors(xy).T)
        return {"bow_max": bow_max, "bow_rms": bow_rms,
                "dtheta_var": float(np.mean(dtheta ** 2)) if len(dtheta) else 0.0,
                "dl_var": float(np.mean((bl - ell_star) ** 2)),
                "min_sep": float(bl.min()) if len(bl) else np.nan,
                "nnn_gap": min_nnn_gap_star(xy),
                "dtheta": dtheta, "bond_len": bl}

    def finalize(cols):
        n_s = len(cols["pe"])
        dtheta_all = np.concatenate(cols["dtheta"]) if n_s else np.zeros(0)
        dl_all = np.concatenate([bl - ell_star for bl in cols["bond_len"]]) if n_s else np.zeros(0)
        dtheta_var = float(np.mean(dtheta_all ** 2)) if dtheta_all.size else float("nan")
        dl_var = float(np.mean(dl_all ** 2)) if dl_all.size else float("nan")
        dl_var_sem = ST.block_sem(np.array([float(np.mean((bl - ell_star) ** 2))
                                            for bl in cols["bond_len"]])) if n_s else float("nan")
        k_bond_star = float(P["k_bond_star"])
        dl_var_pred = float(P["dl_var_boltz"])         # * the anharmonic Boltzmann integral -- see below

        obs = []
        # -- implementation_check -- the radial (bond-stretch) equipartition golden test
        # ** The prediction is **not** the harmonic approximation (kT*/k_bond*). A
        #   smoke run measured the harmonic approximation underestimating the true
        #   value by 4.6x (the well is asymmetric and softer on the outside -- the same
        #   trap as soft-r3's Einstein cage approximation). Instead it uses the value
        #   `bond_variance_boltzmann()` in build_ledger() obtains by numerically
        #   integrating over the whole well (the basin inside the barrier) -- still
        #   **a prediction that follows directly from the same model this case builds
        #   (DLVO+WCA)**, so implementation_check is the right role: it simply does not
        #   assume harmonicity, and is the exact result for this potential. Its job is
        #   to confirm that the DLVO table potential was ported into this case
        #   correctly (chain-bend-2d-dlvo mixed in a trap and driving, so it never
        #   measured this degree of freedom in isolation).
        obs.append(MET.observable(
            "bond-stretch equipartition <dl*^2>", dl_var, dl_var_pred, "d²",
            "boltzmann_integral", role="implementation_check", scope="module", tol_pct=8.0,
            sigma=dl_var_sem if dl_var_sem > 0 else None,
            note=f"prediction = the anharmonic Boltzmann integral over the basin (exact). "
                 f"For reference, the harmonic approximation "
                 f"kT*/k_bond*={1.0 / k_bond_star:.3e} is only "
                 f"1/{dl_var_pred * k_bond_star:.2f} of this value -- the harmonic "
                 f"approximation must not be used as the prediction (regime indicator only)"))
        # -- measurement -- the bending (angular) degrees of freedom. G1 says only
        #    that the linear stiffness is 0; it does not predict the magnitude (higher
        #    order terms dominate). No prior prediction, so measurement.
        obs.append(MET.observable(
            "bending-angle thermal fluctuation <dtheta^2>", dtheta_var, None, "rad²", "none",
            role="measurement",
            note="a direct consequence of G1 (linear stiffness = 0) -- the fluctuation width "
                 "itself has no prediction (second-order and higher terms plus excluded "
                 "volume dominate). Assuming equipartition via kappa_theta,eff=1/<dtheta^2> gives "
                 f"{(1.0 / dtheta_var if dtheta_var > 0 else float('nan')):.3g} kT/rad² — "
                 "a reference number only -- it does not mean the mode is harmonic"))

        min_sep = float(np.min(cols["min_sep"]))
        nnn_min = float(np.min(cols["nnn_gap"]))
        post_checks = [
            C.Check("geometry", "table lower-bound margin r_table_min/min_sep",
                  (1.0 + 1e-6) / min_sep, 1.0, "<=",
                  f"pair.Table trap 11: the force is 0 for r<r_min. Measured minimum bond length {min_sep:.4f}"),
            C.Check("geometry", "min NNN gap (whole run) vs cutoff", nnn_min, CUTOFF_H_STAR, ">=",
                  "did a two-bonds-away bead ever come inside the DLVO cutoff during "
                  "the run (if so, a candidate for 'trianglelike' local folding -- an "
                  "observation, not a failure)",
                  hard=False),
        ]

        bow_final = float(np.mean(cols["bow_rms"][-max(1, n_s // 10):])) if n_s else float("nan")
        bow_final_sem = ST.block_sem(cols["bow_rms"][-max(1, n_s // 10):]) if n_s >= 8 else float("nan")
        post_dicts = [{**c.as_dict("post_run"), "note": c.note} for c in post_checks]
        return {"observables": obs,
                "extra": {"dtheta_var": dtheta_var, "dl_var": dl_var, "dl_var_sem": dl_var_sem,
                          "dl_var_pred": dl_var_pred, "k_bond_star": k_bond_star,
                          "ell_star": ell_star, "min_sep": min_sep, "nnn_min": nnn_min,
                          "bow_initial": float(cols["bow_rms"][0]) if n_s else float("nan"),
                          "bow_final_rms": bow_final, "bow_final_sem": bow_final_sem,
                          "post_checks": post_dicts,
                          "post_checks_ok": all(c.ok for c in post_checks)},
                "arrays": {**{k: cols[k] for k in
                              ("bow_max", "bow_rms", "min_sep", "nnn_gap", "pe")},
                          "dl_flat": dl_all, "dtheta_flat": dtheta_all}}

    phases = None
    if init == "kink":
        phases = [RUN.Phase("release", n_prod, sample_every=sample_every, collect=True,
                            expect_steady=False,
                            note="recorded from the moment of release -- the bow changing IS the physics "
                                 "(it is the observable), and running a drift check here would "
                                 "call a discovery a warning (rule 7')")]

    return RUN.Build(
        sim=sim, forces=[tab, wca], n_particles=n,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq if init != "kink" else 0, n_prod=n_prod,
        sample_every=sample_every, phases=phases or [],
        tags=["2D", "dlvo_secondary_minimum", "WCA_core", "no_bending", "no_friction",
             "no_drive", "chain", init],
        physical={"n_beads": n, "init": init, "kink_angle": kink_angle},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# visualization -- results are shown as graphs (CLAUDE.md)
# ════════════════════════════════════════════════════════════════════════
def make_plots(sys_, lg, n, init, kink_angle, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]     # * labels in English (CLAUDE.md)
    matplotlib.rcParams["axes.unicode_minus"] = False

    res = np.load(outdir / "observables.npz")
    m = json.loads((outdir / "metrics.json").read_text())
    D = lg.derived
    dt_over_tau_bond = lg.ratio("times", "dt", "tau_bond")     # * this system's natural time axis
    t = np.arange(len(res["bow_rms"])) * (m["numerics"].get("sample_every", 1)) * dt_over_tau_bond

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    ax[0, 0].plot(t, res["bow_rms"], "-", lw=1.0, label="bow rms (body frame)")
    ax[0, 0].plot(t, res["bow_max"], "-", lw=0.6, alpha=0.6, label="bow max")
    ax[0, 0].set(xlabel=r"$t / \tau_{bond}$", ylabel="bow [d]",
                title=f"Shape relaxation — n={n}, init={init}"
                      + (f", kink={kink_angle:.2f} rad" if kink_angle else ""))
    ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

    ax[0, 1].plot(t, res["pe"], "-", lw=0.8)
    ax[0, 1].set(xlabel=r"$t / \tau_{bond}$", ylabel=r"$\langle U \rangle / N$ [kT]",
                title="Potential energy per particle")
    ax[0, 1].grid(alpha=.3)

    extra = m.get("result", {})
    dtheta_var = extra.get("dtheta_var", float("nan"))
    dl_var = extra.get("dl_var", float("nan"))
    dl_var_pred = extra.get("dl_var_pred", float("nan"))
    k_bond_star = extra.get("k_bond_star", 1.0)
    if "dl_flat" in res:
        ax[1, 0].hist(res["dl_flat"], bins=60, density=True, alpha=.45,
                      color="tab:blue", label="measured (all bonds x samples)")
    harm_sigma = math.sqrt(1.0 / k_bond_star)
    boltz_sigma = math.sqrt(dl_var_pred) if dl_var_pred > 0 else harm_sigma
    xs2 = np.linspace(-5 * boltz_sigma, 5 * boltz_sigma, 400)
    ax[1, 0].plot(xs2, np.exp(-xs2 ** 2 / (2 * boltz_sigma ** 2)) / math.sqrt(2 * math.pi) / boltz_sigma,
                 "k-", lw=1.4, label=f"predicted (basin Boltzmann integral)  σ={boltz_sigma:.4f} d")
    ax[1, 0].plot(xs2, np.exp(-xs2 ** 2 / (2 * harm_sigma ** 2)) / math.sqrt(2 * math.pi) / harm_sigma,
                 "k:", lw=1.0, alpha=.6, label=f"harmonic approx (NOT the prediction)  σ={harm_sigma:.4f} d")
    ax[1, 0].axvline(0, color="gray", lw=.5)
    ax[1, 0].set(xlabel=r"$\delta\ell^* = \ell^* - \ell^*_{nat}$  [d]", ylabel="density",
                title=f"Bond-length golden test — measured Var={dl_var:.3e} d²"
                      f" vs predicted {dl_var_pred:.3e} d²")
    ax[1, 0].legend(fontsize=7); ax[1, 0].grid(alpha=.3)

    ax[1, 1].plot(t, res["min_sep"], "-", lw=.8, label="min bond length")
    ax[1, 1].axhline(1.0 + D["h_min_star"], color="tab:green", ls=":", label="natural length")
    ax[1, 1].plot(t, res["nnn_gap"] + 1.0, "-", lw=.6, alpha=.7, label="min NNN separation")
    ax[1, 1].axhline(1.0 + float(CUTOFF_H_STAR), color="tab:red", ls=":", label="DLVO cutoff")
    ax[1, 1].set(xlabel=r"$t / \tau_{bond}$", ylabel="r* [d]",
                title="Bond safety / non-adjacent folding")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3)

    # * The bending (angular) degrees of freedom -- exactly the observable G1 says it
    # has no prediction for. With no prior prediction, no reference line is drawn here
    # (it is a measurement). Instead, for reference, it is shown on the same angular
    # scale so it can be compared visually against the bond stretch's golden-test
    # sigma (harmonic): if the bending fluctuation width is overwhelmingly larger than
    # the stretch fluctuation width (a direct comparison is impossible -- radians vs a
    # dimensionless length), that is itself a qualitative signal that there is no
    # linear stiffness -- having no lower bound means it is not held by any scale the
    # system defines.
    if "dtheta_flat" in res and res["dtheta_flat"].size:
        ax[0, 2].hist(res["dtheta_flat"], bins=60, density=True, color="tab:orange", alpha=.6)
    ax[0, 2].axvline(0, color="gray", lw=.5)
    ax[0, 2].set(xlabel=r"$d\theta_i$ [rad] (local bend, all internal beads x samples)",
                ylabel="density",
                title=f"Bend-angle fluctuation — measured $\\langle d\\theta^2\\rangle$="
                      f"{dtheta_var:.3e} rad$^2$  (G1: no prior prediction)")
    ax[0, 2].grid(alpha=.3)

    ax[1, 2].axis("off")
    ax[1, 2].text(0.02, 0.95,
                  "G1 check (this run)\n"
                  "─────────────────\n"
                  f"bond stretch  Var={dl_var:.3e} d²\n"
                  f"  predicted   ={dl_var_pred:.3e} d²  (basin Boltzmann)\n"
                  f"  err         ={100*(dl_var/dl_var_pred-1):+.1f}%"
                  f"  (implementation_check)\n\n"
                  f"bend angle    Var={dtheta_var:.3e} rad²\n"
                  f"  no prior prediction (measurement)\n"
                  f"  naive equipartition kappa_eff="
                  f"{(1.0/dtheta_var if dtheta_var > 0 else float('nan')):.3g} kT/rad²\n"
                  f"  (reference number only — not a claim of harmonicity)\n\n"
                  + (f"bow: initial={extra.get('bow_initial', float('nan')):.4f} d, "
                     f"final(rms, last decile)={extra.get('bow_final_rms', float('nan')):.4f}"
                     f" +/- {extra.get('bow_final_sem', float('nan')):.4f} d\n"
                     if init == "kink" else ""),
                  transform=ax[1, 2].transAxes, fontsize=9, va="top", family="monospace")

    fig.suptitle(f"{sys_['label']}  n={n}  init={init}", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", choices=("straight", "kink"), required=True)
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n", type=int, default=None, help="number of beads (default: system.yaml, n=9)")
    ap.add_argument("--kink-angle", type=float, default=0.30, help="rad (--init kink only)")
    ap.add_argument("--tobs", type=float, default=None,
                    help="observation window (multiples of tau_bond). Default: straight=2e4, kink=3e3")
    ap.add_argument("--cycles", type=float, default=None)   # unused -- kept for interface symmetry
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dt-scale", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eq-scale", type=float, default=200.0,
                    help="(--init straight) equilibration = this value x tau_bond/dt")
    ap.add_argument("--smoke", action="store_true", help="for a quick sanity check -- steps heavily reduced")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/chain-relax-2d-dlvo/system.yaml")
    n = args.n if args.n is not None else sys_["n_list"][0]
    kink_angle = args.kink_angle if args.init == "kink" else 0.0

    if args.smoke:
        args.tobs = min(args.tobs or 50.0, 50.0)
        args.samples = min(args.samples, 100)
        args.eq_scale = min(args.eq_scale, 20.0)

    if not (args.report or args.spec or args.run):
        print("choose what to do -- `--report`, `--spec` or `--run`")
        return 3

    return emit(sys_, n, args.init, kink_angle, args)


if __name__ == "__main__":
    sys.exit(main())
