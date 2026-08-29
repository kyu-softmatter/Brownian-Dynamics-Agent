"""`network` -- a 3D colloidal network. **Stage 1: compression gelation**
(CLAUDE.md rule 8).

The sketch (intake/network/) says only "colloidal network / identical bead /
DLVO, & JKR / stress propagation / x(t) = A sin(omega*t)" and **contains no numbers
at all.** Every material property is inherited from `chain-bend-2d-dlvo` (see the
tiers in system.yaml).

What this file does -- **only up to building the network**:
    1. scatter into a 3D box (phi0=0.02) without overlap
    2. aggregate via the DLVO secondary minimum (pair potential only -- no topology
       is declared)
    3. compress the box isotropically over 178 stages to reach phi=0.10
    4. relax again after compression and measure the **structure**:
       z, independent loops, free-end fraction, percolation, d_f, g(r), and the
       fraction of collapsed bonds

The driving (x(t)=A sin(omega*t)) and G'(omega) are **stage 2** and are not here.
Why the split:
  * rule 8: build the static system first, add motion afterwards
  * HOOMD's bond/angle are **static topology** and cannot be declared during
    aggregation. The contact list has to be extracted and the topology frozen after
    gelation before stage 2 exists
  * JKR bending is stage 2 for the same reason
    (`interactions[1].enabled_stage: 2`)

** The central constraint on compression, measured before running
   (verify/verify_3d_boxresize.py):
    `update.BoxResize` **affinely scales** the coordinates (error 8.9e-16), so an
    already-bonded pair's bond length shrinks with it. Once the linear strain per
    trigger exceeds (h_min - h_barrier)/ell = 0.703%, the pair is pushed inside the
    barrier and **collapses irreversibly into the primary minimum (contact).**
    Measured: 0.40% holds (h*=0.007591) / 0.80% collapses (h*=-0.042257).
    -> Use 0.4% per stage. `Interpolate` is linear in L, so the worst case is the
      last stage, which makes the stage count **178**, not 134
      (system.yaml geometry.n_stages).

Usage:
    $PY cases/network_3d.py --report                 # the L3 report only
    $PY cases/network_3d.py --spec                   # save the spec (does not run)
    $PY cases/network_3d.py --n 512 --stage-tau 1e-3 # a pilot run
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
from bdbot import nondim as ND, run as RUN, scales as SC, sim as SIM        # noqa: E402
from bdbot.provenance import load_node                                      # noqa: E402
from bdbot.units import Q                                                   # noqa: E402

# The DLVO expressions are used **as-is** from chain-bend-2d-dlvo, where they were
# verified in SI -- they are not written twice
from chain_bend_dlvo_2d import (                                            # noqa: E402
    SIGMA_CORE_STAR, build_table_arrays, dlvo_reduced_params, find_well, U_star,
)

ROOT = Path(__file__).resolve().parent.parent
CUTOFF_H_STAR = 0.06          # pair-table cutoff (surface gap h/d). 8x outside the
                              # secondary minimum (0.0076)
EPS_PER_STAGE = 0.004         # allowed linear strain per compression trigger
                              # (the measured safe line)
R_WCA = 2 ** (1 / 6)


# ════════════════════════════════════════════════════════════════════════
# 1. the physical system (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    P = load_node
    dlvo, jkr = raw["interactions"][0], raw["interactions"][1]
    g = raw["geometry"]
    return {
        "label": raw["label"], "dim": int(raw["dimensions"]),
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "psi0": P(raw["particle"]["surface_potential"]),
        "n_list": [int(x) for x in raw["particle"]["count"]["value"]],
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "eps_r": P(raw["medium"]["relative_permittivity"]),
        "ionic_strength": P(raw["medium"]["ionic_strength"]),
        "A_H": P(dlvo["hamaker_constant"]),
        "kappa_theta": P(jkr["angle_stiffness"]),
        "phi0": float(g["volume_fraction_initial"]["value"]),
        "phi1": float(g["volume_fraction_final"]["value"]),
        "n_stages": int(g["n_stages"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


def bond_edge_h_star(p: dict, w: dict) -> float:
    """The bonding boundary -- the h* on the outer flank where U returns to **half**
    the well depth.

    An arbitrary cutoff would make z, the loop count and the free-end fraction wobble
    with that choice. Half the well depth is a value that comes out of the ledger, so
    it is reproducible across cases.
    """
    h = np.geomspace(w["h_min"], CUTOFF_H_STAR, 200_000)
    U = U_star(h, p)
    return float(h[int(np.argmin(np.abs(U - w["U_min"] / 2.0)))])


def box_star(n: int, phi: float) -> float:
    """The cube edge in units of d, from phi. N*(pi/6)*d^3 / L^3 = phi."""
    return (n * math.pi / (6.0 * phi)) ** (1.0 / 3.0)


# ════════════════════════════════════════════════════════════════════════
# 2. the scale ledger
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, n: int, *, dt_scale=1.0, stage_tau=1e-3,
                 agg_tau=0.2, post_tau=0.2, init="scatter", max_coord=4,
                 n_seeds=1, loop_bias=1) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B, tau_p = b["kT"], b["gamma"], b["tau_B"], b["tau_p"]

    p = dlvo_reduced_params(sys_)
    w = find_well(p)
    h_edge = bond_edge_h_star(p, w)
    ell = d * (1.0 + w["h_min"])
    k_bond = Q(w["k_bond_star"], "dimensionless") * kT / d ** 2
    sigma_bond = (kT / k_bond) ** 0.5
    lam_D = (Q(1.0, "dimensionless") / (p["kappa_star"] / d)).to("nm")

    tau_bond = C.relaxation_time(gamma, k_bond)
    # * The bond escape time -- the well is 11.7 kT, so bonding is **reversible.**
    #   The reference time for quasi-staticity of the compression is this, not tau_B
    #   (Kramers estimate: prefactor tau_bond, exponent |U_min|/kT).
    tau_esc = tau_bond * math.exp(abs(w["U_min"]))

    # dt is set by the fastest mode. In stage 1 the only candidate is bond stretch
    # (bending and the trap are stage 2).
    # WARNING: the seat a collapsed pair settles into (h*<0, the WCA-vdW balance) is
    #    stiffer -- so if collapse occurs, dt has to be revisited, which is why
    #    crushed_bond_fraction is emitted as an observable.
    tau_fast = tau_bond
    dt = dt_scale * C.dt_from_gate(tau_fast)

    # * `init="sprout"` -- build the network directly at the target phi, so there is
    #   neither aggregation nor compression. L0 = L1 and the stage count is 0. (A user
    #   suggestion: if compression is hard, sprout outward in random directions from a
    #   seed particle. Compression was not in fact hard, but aggregation barely
    #   progressed at phi0=0.02 -- measured -- so this route is better.)
    sprout = (init == "sprout")
    L1 = box_star(n, sys_["phi1"])
    L0 = L1 if sprout else box_star(n, sys_["phi0"])
    n_stages = 0 if sprout else sys_["n_stages"]
    r_cut = 1.0 + CUTOFF_H_STAR
    # physical time per compression stage (stage_tau is in units of tau_B) -> the
    # observation window
    T_stage = Q(stage_tau, "dimensionless") * tau_B
    T_compress = n_stages * T_stage
    T_obs = (Q(0.0 if sprout else agg_tau, "dimensionless") * tau_B + T_compress
             + Q(post_tau, "dimensionless") * tau_B)

    lg = SC.ScaleLedger(ref=SC.thermal_reference(
        d, kT, tau_B,
        rationale="the d, kT, tau_B thermal references. In HOOMD this becomes "
                             "sigma=kT=gamma=1, and compression changes only the box."))
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("ell", ell, "DLVO bond length (centre to centre, natural)", star=True)
    lg.add_length("h_min", d * w["h_min"], "secondary-minimum surface gap")
    lg.add_length("h_edge", d * h_edge, "bonding boundary (U=well/2)")
    lg.add_length("sigma_bond", sigma_bond.to("nm"), "bond thermal fluctuation width")
    lg.add_length("lambda_D", lam_D, "Debye length")
    lg.add_length("r_cut", d * r_cut, "pair-table cutoff")
    lg.add_length("L0", d * L0, f"initial box (phi={sys_['phi0']:.3f})")
    # * role="box" takes the box **after** compression -- the worst case for the
    #   geometry check (r_cut <= L/2)
    lg.add_length("L", d * L1, f"box edge (after compression, phi={sys_['phi1']:.3f})", role="box")

    lg.add_time("tau_B", tau_B, "Brownian time (reference)")
    lg.add_time("tau_p", tau_p, "momentum relaxation", role="inertia")
    lg.add_time("tau_bond", tau_bond, "* DLVO bond stretch -- the fastest mode", star=True)
    lg.add_time("tau_esc", tau_esc, "* bond escape (structural relaxation of a reversible gel)", star=True)
    lg.add_time("T_stage", T_stage, "physical time of one compression stage")
    lg.add_time("T_compress", T_compress if n_stages else dt,
                f"the whole compression ({n_stages} stages)" if n_stages else "no compression (sprout)")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("T_obs", T_obs, "observation window (aggregation + compression + post-relaxation)", role="observation")

    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("well_depth", Q(-w["U_min"], "dimensionless") * kT, "secondary-minimum depth (sign inverted)")
    lg.add_energy("barrier", Q(w["barrier_U"], "dimensionless") * kT, "DLVO barrier")
    lg.add_energy("k_bond_d2", k_bond * d ** 2, "bond stretch stiffness x d^2")

    total_strain = 1.0 - L1 / L0
    # with no compression the affine strain is 0, so there is no risk of breaking a bond at all
    eps_max = 0.0 if n_stages == 0 else (L0 / L1 - 1.0) / n_stages
    eps_crit = (w["h_min"] - w["barrier_h"]) / (1.0 + w["h_min"])

    lg.derived.update(
        reduced=p, well=w, h_edge=h_edge, ell=ell, k_bond=k_bond, tau_fast=tau_fast,
        tau_esc=tau_esc, tau_B=tau_B, T_obs=T_obs, T_stage=T_stage,
        L0_star=L0, L1_star=L1, r_cut_star=r_cut, r_bond_star=1.0 + h_edge,
        total_strain=total_strain, eps_max=eps_max, eps_crit=eps_crit, init=init,
        sprout=sprout, max_coord=max_coord, n_seeds=n_seeds, loop_bias=loop_bias,
        n_stages=n_stages, stage_tau=stage_tau,
        agg_tau=0.0 if sprout else agg_tau, post_tau=post_tau,
        U_min_star=w["U_min"], barrier_star=w["barrier_U"], h_min_star=w["h_min"],
        k_bond_star=w["k_bond_star"])
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 3. dimensionless groups + checks
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg, n):
    D = lg.derived
    r = lg.ratio
    groups = [
        ND.Group("phi_final", sys_["phi1"], None, None,
                 "N(pi/6)d^3/L^3", "volume fraction -- a sparse gel"),
        ND.Group("well_depth/kT", r("energies", "well_depth", "kT"),
                 ("energies", "well_depth"), ("energies", "kT"), "",
                 "* 11.7 -> bonding is reversible (tau_esc = 0.115 tau_B)"),
        ND.Group("barrier/kT", r("energies", "barrier", "kT"),
                 ("energies", "barrier"), ("energies", "kT"), "",
                 "the transition to the primary minimum is effectively impossible"),
        ND.Group("k_bond_star", r("energies", "k_bond_d2", "kT"),
                 ("energies", "k_bond_d2"), ("energies", "kT"), "",
                 "bond stretch stiffness -- it sets dt"),
        ND.Group("tau_bond/tau_B", r("times", "tau_bond", "tau_B"),
                 ("times", "tau_bond"), ("times", "tau_B"), "",
                 "* six decades of scale separation -- the source of the dt cost"),
        ND.Group("tau_esc/tau_B", r("times", "tau_esc", "tau_B"),
                 ("times", "tau_esc"), ("times", "tau_B"), "",
                 "bond rearrangement / diffusion -- the criterion for quasi-static compression"),
        ND.Group("T_stage/tau_esc", r("times", "T_stage", "tau_esc"),
                 ("times", "T_stage"), ("times", "tau_esc"), "",
                 "* >=1 is quasi-static compression, <<1 is affine compression (different physics)"),
        ND.Group("ell/d", r("lengths", "ell", "d"),
                 ("lengths", "ell"), ("lengths", "d"), "", "bond length / diameter"),
        ND.Group("L/d", r("lengths", "L", "d"),
                 ("lengths", "L"), ("lengths", "d"), "", "box after compression (effective propagation radius L/2)"),
        ND.Group("eps_max", D["eps_max"], None, None,
                 "(L0/L1-1)/n_stages", "* maximum linear strain per compression stage"),
    ]
    checks = [
        C.Check("integration", "dt / tau_fast", lg.ratio("times", "dt", "tau_bond"), C.GATE,
                note="Brownian is O(dt) -- the fastest mode is bond stretch"),
        # * Treating this as soft follows the established convention from
        #   `chain-bend-2d-dlvo` -- same particles, same DLVO bond, so this ratio is
        #   **literally the same value (0.0282)** and that case already documented the
        #   same judgment. The check is kept and reported rather than removed.
        C.Check("model", "note: tau_p/tau_bond", lg.ratio("times", "tau_p", "tau_bond"),
                C.GATE, "<=",
                "* the bond well is deep and narrow, so tau_bond approaches tau_p "
                "(exceeding the inertia-negligible criterion by 2.82x -- a property of "
                "**the bond physics itself**, independent of N, phi and compression). "
                "But zeta = (1/2)(tau_p/tau_bond)^(-1/2) = 2.98, so the overdamped "
                "condition (zeta>1) is satisfied with margin. In chain-bend-2d-oscill a "
                "far worse zeta=0.65 (tau_p/tau_fast=0.60) was compared "
                "OverdampedViscous vs Langevin(kT=0) and measured to affect the "
                "observables by at most 0.159% -- this case is 21x more favourable. "
                "WARNING: that comparison was not actually run for this system "
                "(recorded in not_verified)",
                hard=False),
        C.Check("model", "bond stability  sigma_bond/h_min",
                lg.ratio("lengths", "sigma_bond", "h_min"), 0.5, "<=",
                "* if the bond thermal fluctuation width exceeds half the well "
                "position, bonds keep breaking thermally -- that may itself be the "
                "result, but it is flagged in advance (the same check as chain-bend)",
                hard=False),
        C.Check("geometry", "r_cut / (L/2)",
                lg.ratio("lengths", "r_cut", "L") * 2.0, 1.0,
                note="* judged against the box **after** compression (bd-hoomd trap 6)"),
        C.Check("geometry", "lambda_D / (L/4)",
                lg.ratio("lengths", "lambda_D", "L") * 4.0, 1.0,
                note="does the double layer fit in the box"),
        C.Check("geometry", "ell / (L/2)", lg.ratio("lengths", "ell", "L") * 2.0, 1.0,
                note="the bond length is inside minimum image"),
        # ** this case's central check -- it forces the measured threshold into the design
        C.Check("integration", "eps_max / eps_crit", D["eps_max"] / D["eps_crit"], 1.0,
                note="* the condition under which compression does not break a bond. "
                     "Measured: 0.40% holds / 0.80% collapses"
                     + ("  (sprout -- no compression, so affine strain is 0)" if D["sprout"] else "")),
        C.Check("statistics", "N", float(n), 512.0, op=">=", hard=False,
                note="structural statistics -- 512 for the pilot, 1528 for production"),
        C.Check("statistics", "T_stage / tau_bond",
                lg.ratio("times", "T_stage", "tau_bond") if not D["sprout"] else 1e9,
                10.0, op=">=", hard=False,
                note="the time a bond has to relax per stage. Below 1, the collapse verdict is meaningless"
                     + ("  (sprout -- not applicable)" if D["sprout"] else "")),
        C.Check("finite-size", "L / ell", lg.ratio("lengths", "L", "ell"), 8.0, op=">=",
                hard=False, note="how many bond lengths the box spans -- does the mesh fit"),
    ]
    return groups, checks


def report_blocks(sys_, lg, n, n_agg, n_stage, n_post):
    D = lg.derived
    inp = [R.kv("d", f"{sys_['d'].value:~.4gP}", sys_["d"].tier, sys_["d"].source[:44]),
           R.kv("T", f"{sys_['T'].value:~.4gP}", sys_["T"].tier, "water"),
           R.kv("eta", f"{sys_['eta'].value:~.4gP}", sys_["eta"].tier, "water@300K"),
           R.kv("A_H", f"{sys_['A_H'].value:~.4gP}", sys_["A_H"].tier, sys_["A_H"].source[:44]),
           R.kv("I", f"{sys_['ionic_strength'].value:~.4gP}",
                sys_["ionic_strength"].tier, "MgCl2"),
           R.kv("N", str(n), 3, "*proposed -- the 21 in the figure are schematic (A6)"),
           R.kv("phi", f"{sys_['phi0']:.3f} → {sys_['phi1']:.3f}", 1, "compression gelation")]
    der = [f"  bond: barrier {D['barrier_star']:.2f} kT   secondary minimum {D['U_min_star']:.3f} kT"
           f" @ h*={D['h_min_star']:.6f}",
           f"  bonding boundary r*<{D['r_bond_star']:.6f} (U=well/2, h*={D['h_edge']:.6f})",
           f"  k_bond* {D['k_bond_star']:.4g} kT/d²  → τ_bond {D['tau_fast']:~.4gP}",
           f"  τ_esc {D['tau_esc']:~.4gP} = {float((D['tau_esc']/D['tau_B']).to('')):.4f} τ_B"
           f"  <- the structural relaxation time of a reversible gel",
           f"  box {D['L0_star']:.3f} d -> {D['L1_star']:.3f} d"
           f"   (total linear strain {D['total_strain']*100:.2f}%)",
           (f"  topology generator **sprout** -- no compression (affine strain 0, zero risk of "
            f"bond collapse). The {D['eps_crit']*100:.4f}% threshold does not apply"
            if D["eps_max"] == 0 else
            f"  compression in {D['n_stages']} stages . max strain per stage {D['eps_max']*100:.4f}%"
            f"  (threshold {D['eps_crit']*100:.4f}%, margin {D['eps_crit']/D['eps_max']:.2f}x)")]
    if D["sprout"]:
        plan = [f"  1. sprouting (not MD)   seeds -> random directions, placed exactly at ell, "
                f"max_coord={D['max_coord']}, seeds={D['n_seeds']}, loop_bias={D['loop_bias']}",
                f"  2. relax   {n_post:>12,} steps  ({D['post_tau']:.4g} tau_B, phi={sys_['phi1']:.3f})",
                f"  total      {n_post:>12,} steps",
                "  * no compression -> affine strain 0, zero risk of bond collapse",
                "  WARNING: this topology is **imposed**, not produced by DLVO dynamics -- "
                "it has to be compared against the compression route before topology "
                "dependence is visible (rule 7')"]
    else:
        plan = [f"  1. aggregate {n_agg:>12,} steps  ({D['agg_tau']:.4g} tau_B, phi={sys_['phi0']:.3f} fixed)",
                f"  2. compress  {n_stage*D['n_stages']:>12,} steps  "
                f"({D['n_stages']} × {n_stage:,} = {D['stage_tau']*D['n_stages']:.4g} τ_B)",
                f"  3. relax     {n_post:>12,} steps  ({D['post_tau']:.4g} tau_B, phi={sys_['phi1']:.3f})",
                f"  total        {n_agg+n_stage*D['n_stages']+n_post:>12,} steps",
                f"  ★ T_stage/τ_esc = {D['stage_tau']/float((D['tau_esc']/D['tau_B']).to('')):.4g}"
                f"  ({'quasi-static' if D['stage_tau'] >= float((D['tau_esc']/D['tau_B']).to('')) else 'affine compression -- diffusion does not participate'})"]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# L3 -- the spec
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, n, args):
    lg = build_ledger(sys_, n, dt_scale=args.dt_scale, stage_tau=args.stage_tau,
                      agg_tau=args.agg_tau, post_tau=args.post_tau,
                      init=args.init, max_coord=args.max_coord,
                      n_seeds=args.n_seeds, loop_bias=args.loop_bias)
    D = lg.derived
    dt = lg.get("times", "dt")
    per = lambda tau: max(1, int(round(float((Q(tau, "dimensionless") * D["tau_B"] / dt).to("")))))
    n_agg = 0 if D["sprout"] else per(args.agg_tau)
    n_stage, n_post = per(args.stage_tau), per(args.post_tau)
    n_prod = n_agg + n_stage * D["n_stages"] + n_post
    sample_every = max(1, n_prod // args.samples)

    groups, checks = analyze_scales(sys_, lg, n)
    # * the tag must separate different physics -- init changes the topology
    #   generator, so it has to be included
    if D["sprout"]:
        tag = (f"N{n}-sprout-mc{args.max_coord}-sd{args.n_seeds}"
               f"-lb{args.loop_bias}-po{args.post_tau:g}")
    else:
        tag = f"N{n}-st{args.stage_tau:g}-ag{args.agg_tau:g}"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"
    p = D["reduced"]
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_particles": n, "phi0": sys_["phi0"], "phi1": sys_["phi1"],
                "L0_star": D["L0_star"], "L1_star": D["L1_star"],
                "n_stages": D["n_stages"], "eps_max": D["eps_max"], "eps_crit": D["eps_crit"],
                "kappa_star": p["kappa_star"], "edl_amp": p["edl_amp"],
                "vdw_amp": p["vdw_amp"], "a_star": p["a_star"],
                "cutoff_h_star": CUTOFF_H_STAR,
                "h_min_star": D["h_min_star"], "h_edge_star": D["h_edge"],
                "r_bond_star": D["r_bond_star"],
                "well_depth_star": D["U_min_star"], "barrier_star": D["barrier_star"],
                "k_bond_star": D["k_bond_star"],
                "tau_esc_star": float((D["tau_esc"] / D["tau_B"]).to("")),
                "stage_tau": args.stage_tau, "agg_tau": D["agg_tau"],
                "post_tau": args.post_tau,
                # * the topology generator -- it fixes the physics, so it must be in the hash
                "init": D["init"], "max_coord": int(args.max_coord),
                "n_seeds": int(args.n_seeds), "loop_bias": int(args.loop_bias),
                "stage": 1, "jkr": False, "drive": "none"},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "n_eq": 0, "n_prod": n_prod, "n_agg": n_agg,
                  "n_stage": n_stage, "n_post": n_post,
                  "n_samples": args.samples, "sample_every": sample_every,
                  "seed": args.seed},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_agg, n_stage, n_post


def emit(sys_, n, args) -> int:
    lg, spec, groups, checks, n_agg, n_stage, n_post = build_spec(sys_, n, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n, n_agg, n_stage, n_post)
    report, verdict = R.render(
        title=f"DimensionlessReport -- {sys_['label']} stage 1 (gelation)  N={n}"
              f"  T_stage={args.stage_tau:g}τ_B   run_id={run_id}",
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

    # `result.txt` is the project's completion marker, and writing it is the case
    # script's responsibility -- without it `status` counts the run as zero and a
    # cleanup pass deletes a completed run
    if v["status"] != "skipped":
        lines = []
        try:
            mj = json.loads((outdir / "metrics.json").read_text())
            for o in mj.get("observables", []):
                m, pr = o.get("measured"), o.get("predicted")
                tail = f"   (predicted {pr:.6g})" if isinstance(pr, (int, float)) else ""
                lines.append(f"  {o['name']:<26} {m:.6g}{tail}" if m is not None
                             else f"  {o['name']:<26} —")
        except Exception as e:                        # noqa: BLE001
            lines.append(f"  (could not read metrics.json: {e})")
        (outdir / "result.txt").write_text(
            report + "\n" + "\n".join(["=" * 84, f"RESULT -- {run_id}", "=" * 84,
                                       *lines, "=" * 84, verdict_txt]))
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


# ════════════════════════════════════════════════════════════════════════
# L4 -- structural analysis (pure functions, testable independently of a run)
# ════════════════════════════════════════════════════════════════════════
def contacts(pos: np.ndarray, L: float, r_bond: float):
    """The list of bonded pairs under PBC minimum image, plus each pair's surface gap.

    O(N^2) -- adequate for N <~ 2000.
    """
    n = len(pos)
    d = pos[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d ** 2).sum(-1))
    iu = np.triu_indices(n, 1)
    rr = r[iu]
    m = rr < r_bond
    return np.column_stack([iu[0][m], iu[1][m]]), rr[m], rr


def topology(n: int, pairs: np.ndarray):
    """Coordination number, connected components, independent loops (Betti-1), free ends."""
    deg = np.zeros(n, dtype=int)
    adj = [[] for _ in range(n)]
    for i, j in pairs:
        deg[i] += 1
        deg[j] += 1
        adj[i].append(j)
        adj[j].append(i)
    seen = np.zeros(n, dtype=bool)
    comps, sizes = 0, []
    for s in range(n):
        if seen[s]:
            continue
        comps += 1
        cnt, stack = 0, [s]
        while stack:
            u = stack.pop()
            if seen[u]:
                continue
            seen[u] = True
            cnt += 1
            stack.extend(v for v in adj[u] if not seen[v])
        sizes.append(cnt)
    e = len(pairs)
    return dict(z=float(2 * e / n), n_links=e, n_components=comps,
                loops=int(e - n + comps), dangling=float((deg == 1).sum() / n),
                isolated=float((deg == 0).sum() / n),
                largest_cluster=float(max(sizes) / n) if sizes else 0.0)


def percolates(pos: np.ndarray, L: float, pairs: np.ndarray) -> float:
    """The PBC spanning verdict -- walk the bonds looking for a loop whose accumulated
    image offset is non-zero.

    BFS while summing minimum-image displacements; if an already-visited node is
    reached with a **different accumulated displacement**, that component wraps the
    box and joins onto itself = percolation.
    Returns: the number of percolating axes / 3.
    """
    n = len(pos)
    adj = [[] for _ in range(n)]
    for i, j in pairs:
        dv = pos[j] - pos[i]
        dv -= L * np.round(dv / L)
        adj[i].append((j, dv))
        adj[j].append((i, -dv))
    axes = np.zeros(3, dtype=bool)
    off = np.full((n, 3), np.nan)
    for s in range(n):
        if not np.isnan(off[s, 0]):
            continue
        off[s] = 0.0
        stack = [s]
        while stack:
            u = stack.pop()
            for v, dv in adj[u]:
                cand = off[u] + dv
                if np.isnan(off[v, 0]):
                    off[v] = cand
                    stack.append(v)
                else:
                    w = np.abs(cand - off[v])
                    axes |= w > 0.5 * L          # the wrapped axis
    return float(axes.sum() / 3.0)


def fractal_dimension(pos: np.ndarray, L: float, rng) -> float:
    """Mass-radius scaling d_f: N(r) proportional to r^d_f.

    Averaged over several centres, r in [1.5d, L/4].
    """
    rs = np.geomspace(1.5, max(L / 4.0, 2.0), 12)
    cent = pos[rng.choice(len(pos), size=min(40, len(pos)), replace=False)]
    d = cent[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d ** 2).sum(-1))
    cnt = np.array([(r < rr).sum(axis=1).mean() for rr in rs])
    ok = cnt > 1
    if ok.sum() < 4:
        return float("nan")
    return float(np.polyfit(np.log(rs[ok]), np.log(cnt[ok]), 1)[0])


def rdf(r_all: np.ndarray, n: int, L: float, r_max: float, nbins: int = 200):
    """g(r), from the full array of pair distances. r_max <= L/2."""
    h, edges = np.histogram(r_all, bins=nbins, range=(0.8, r_max))
    mid = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * math.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho = n / L ** 3
    return mid, h / (0.5 * n * rho * shell)


# ════════════════════════════════════════════════════════════════════════
# L4 -- the HOOMD builder
# ════════════════════════════════════════════════════════════════════════
def sprout_network(n: int, L: float, ell: float, rng, *, max_coord: int = 4,
                   n_seeds: int = 1, loop_bias: int = 1,
                   min_center: float | None = None,
                   r_bond: float = 1.0180308, tries: int = 3000):
    """* Build the network directly by **sprouting** in random directions from seed
    particles (a user suggestion).

    It removes the two weaknesses of the compression route:
      * aggregation barely progresses at phi0=0.02 (measured: <r_nn> 1.7519 ->
        1.7534, i.e. it **increased**)
      * affine compression risks breaking bonds (the eps_crit threshold) -- here
        there is no compression at all
    Bonds are placed **exactly at the natural length ell**, so the initial stress is
    zero, and it is built directly at the target phi.

    WARNING: **state the epistemic difference** (rule 7'). The compression route's
    topology is produced by DLVO dynamics; this route's topology is **imposed by me.**
    So a structure obtained this way is not a `hypothesis` but **an input** -- to see
    how much stage 2's G'(omega) depends on the topology, **both** generators have to
    be run and compared. That comparison is the real value of this method.

    * On loops: branching from a single seed makes the topology a **tree**, which is
    exactly the failure mode of the sketch figure (zero loops -> G'(omega->0)->0,
    observation A6). So a new particle is placed **as long as it does not overlap** --
    when branches approach each other a bond forms automatically and **a loop
    closes.** The loop count is not forced; it is measured as a result.

    `max_coord`: a particle that reaches this coordination is no longer chosen as a
                 parent (Eden-growth style). Small gives a sparse network with long
                 branches; large gives a dense one.
    `min_center`: the minimum allowed centre-to-centre distance. **The default is ell
        (the bond length).**
        WARNING: with 1.0 (contact) the measured min_sep dropped to 1.00018 -- that is
        **inside** the DLVO barrier (h*=0.000508), so such a pair collapses
        irreversibly into the primary minimum (contact) on the first relaxation step.
        That would be committing **exactly the same accident** that eps_crit prevents
        on the compression route, but through the initial placement instead. With ell,
        every pair sits at or beyond the well minimum, the initial stress is zero, and
        no pair is past the barrier.
    """
    if min_center is None:
        min_center = ell
    pos = np.zeros((n, 3))
    coord = np.zeros(n, dtype=int)
    # * Use several seeds -- with one, growth is a single blob and **no load path
    #   wrapping the periodic boundary appears** (measured: at N=512 all 14 loops had
    #   winding 0 and percolation was 0.00). Scattering the seeds lets branches from
    #   different blobs meet across the boundary and form wrapping loops.
    k = 0
    while k < n_seeds and k < n:
        c = rng.uniform(-L / 2, L / 2, 3) if k else np.zeros(3)
        if k:
            dd = pos[:k] - c
            dd -= L * np.round(dd / L)
            if np.sqrt((dd ** 2).sum(-1)).min() < min_center:
                continue
        pos[k] = c
        k += 1
    stall = 0
    while k < n and stall < tries:
        cand_parents = np.flatnonzero(coord[:k] < max_coord)
        if len(cand_parents) == 0:               # all saturated -- start a new seed
            c = rng.uniform(-L / 2, L / 2, 3)
            dd = pos[:k] - c
            dd -= L * np.round(dd / L)
            if np.sqrt((dd ** 2).sum(-1)).min() < min_center:
                stall += 1
                continue
            pos[k] = c
            k += 1
            stall = 0
            continue
        par = int(rng.choice(cand_parents))
        # * Try loop_bias candidate directions and pick the position that **forms the
        #   most bonds.**
        #   Why this is needed (measured): the window in which a bond forms is
        #   ell=1.00759 to r_bond=1.01803 -- only **0.0104 d wide** -- so a single
        #   random direction forms just the one bond to its parent and has essentially
        #   zero chance of meeting another branch. The topology then approaches a tree
        #   (14 loops out of 512) and **no wrapping load path appears** (percolation
        #   0.00). MD relaxation cannot fix it either: the well is short and does not
        #   pull anything beyond 0.018 d (measured: <r_nn> went 1.0072 -> 1.0069 during
        #   relaxation, with the topology essentially unchanged).
        #   loop_bias=1 reproduces the old behaviour exactly.
        best = None
        for _ in range(max(1, loop_bias)):
            v = rng.normal(size=3)
            v /= np.linalg.norm(v)
            c = pos[par] + ell * v
            c -= L * np.round(c / L)             # PBC
            dd = pos[:k] - c
            dd -= L * np.round(dd / L)
            r = np.sqrt((dd ** 2).sum(-1))
            if r.min() < min_center:             # overlap -> discard this candidate
                continue
            nb = int((r < r_bond).sum())
            if best is None or nb > best[0]:
                best = (nb, c, r)
                if nb >= max_coord:              # no room left to improve
                    break
        if best is None:
            stall += 1
            continue
        _, c, r = best
        pos[k] = c
        # give coordination to every neighbour that actually bonds at this position
        # -> loops get counted
        coord[k] = int((r < r_bond).sum())
        coord[:k][r < r_bond] += 1
        k += 1
        stall = 0
    if k < n:
        raise RuntimeError(f"sprouting failed: only {k}/{n} placed (phi too high, or "
                           f"max_coord={max_coord} too small)")
    return pos


def scatter_no_overlap(n: int, L: float, rng, min_sep: float = 1.25, tries: int = 400):
    """Scatter with no overlap and no immediate bonding (min_sep > r_bond is required
    for the initial z=0).
    """
    pos = np.empty((n, 3))
    k = 0
    for _ in range(tries * n):
        if k == n:
            break
        c = rng.uniform(-L / 2, L / 2, 3)
        if k:
            d = pos[:k] - c
            d -= L * np.round(d / L)
            if np.sqrt((d ** 2).sum(-1)).min() < min_sep:
                continue
        pos[k] = c
        k += 1
    if k < n:
        raise RuntimeError(f"scattering failed: {k}/{n} (phi too high, or min_sep too large)")
    return pos


@RUN.builder("network")
def build(spec, outdir=None) -> RUN.Build:
    import gsd.hoomd
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_particles"])
    L0, L1 = float(P["L0_star"]), float(P["L1_star"])
    n_stages = int(P["n_stages"])
    r_bond = float(P["r_bond_star"])
    r_cut = 1.0 + float(P["cutoff_h_star"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_agg, n_stage, n_post = int(Nm["n_agg"]), int(Nm["n_stage"]), int(Nm["n_post"])
    rng = np.random.default_rng(seed)

    init = str(P.get("init", "scatter"))
    if init == "sprout":
        # * built directly at the target phi -> no aggregation and no compression
        #   (L0 == L1 is enforced by the spec)
        pos0 = sprout_network(n, L1, 1.0 + float(P["h_min_star"]), rng,
                              max_coord=int(P.get("max_coord", 4)),
                              n_seeds=int(P.get("n_seeds", 1)),
                              loop_bias=int(P.get("loop_bias", 1)),
                              r_bond=r_bond)
    else:
        pos0 = scatter_no_overlap(n, L0, rng, min_sep=max(1.25, r_bond * 1.1))

    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.position = pos0
    f.particles.typeid = [0] * n
    f.particles.types = ["A"]
    f.particles.mass = [1.0] * n
    f.configuration.box = [L0, L0, L0, 0, 0, 0]      # * 3D: Lz=L (2D uses 0 -- trap 9)
    f.configuration.dimensions = 3
    sim = SIM.make_sim(f, seed=seed)

    cell = md.nlist.Cell(buffer=0.2)
    reduced = {k: P[k] for k in ("kappa_star", "edl_amp", "vdw_amp", "a_star")}
    r_min = 1.0 + 1e-4
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min, r_cut)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut)
    tab.params[("A", "A")] = dict(r_min=r_min, U=U_arr, F=F_arr)
    # * sigma_c = d*2^(-1/6) so the WCA cutoff ends exactly at r=d (the same
    #   convention as chain-bend)
    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * R_WCA, mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
    forces = [tab, wca]

    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    # ** Compression -- runs only after the aggregation phase ends. `Interpolate` is
    #   linear in L, so the relative strain grows toward the end -> the stage count is
    #   178, not 134 (see system.yaml).
    #   With `init="sprout"` it is already at the target phi, so there is no
    #   compression (n_stages=0).
    if n_stages > 0:
        box_var = hoomd.variant.box.Interpolate(
            initial_box=hoomd.Box(Lx=L0, Ly=L0, Lz=L0),
            final_box=hoomd.Box(Lx=L1, Ly=L1, Lz=L1),
            variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=n_agg,
                                       t_ramp=n_stage * n_stages))
        sim.operations.updaters.append(hoomd.update.BoxResize(
            trigger=hoomd.trigger.Periodic(period=n_stage, phase=n_agg), box=box_var))
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, (n_agg + n_stage * n_stages + n_post) // 200))

    # ** `<U>/N` **cannot be used** as the equilibrium/FROZEN indicator here.
    #   Early in aggregation no pair is inside the cutoff (1.06 d), so PE is **exactly
    #   0**, and health reads that as "there is thermal noise and it does not change"
    #   -> FROZEN. It is a false positive, and an expensive one: `run.execute` calls
    #   `finalize()` only when `verdict == OK`, so **the observables vanish entirely**
    #   (measured in a smoke run -- 0 observables, result None).
    #   Instead it follows the `abp-rod` precedent: use a quantity that **can never be
    #   exactly constant** while the system is alive -- the mean nearest-neighbour
    #   distance <r_nn> [d] over 64 reference particles.
    #     * aggregation -> decreases (bonds forming)   * compression -> decreases
    #     * post-relaxation -> steady state
    #     * if motion stops it pins exactly -> FROZEN fires **correctly**
    #   WARNING: the health report labels this value "<U>/N [kT]" -- the unit is
    #      actually d. PE/N itself is kept separately as a sample (`pe_energy`).
    n_ref = min(64, n)
    ref = np.arange(n_ref)          # get_snapshot() is in tag order (unlike local_snapshot)

    def r_nn_mean():
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, copy=True)
        L = float(snap.configuration.box[0])
        dd = pos[ref][:, None, :] - pos[None, :, :]
        dd -= L * np.round(dd / L)
        r = np.sqrt((dd ** 2).sum(-1))
        r[ref, ref] = np.inf                       # exclude self
        return float(r.min(axis=1).mean())

    def pe_energy():
        return float(sum(np.array(fo.energies).sum() for fo in forces)) / n

    # * Accumulate derived quantities in the closure and keep only the raw arrays that
    #   are needed -- returning the full N-array every sample bloats observables.npz by
    #   hundreds of times (448KB -> 148MB in chain-bend)
    series: dict = {k: [] for k in
                    ("t", "L", "phi", "z", "loops", "dangling", "isolated",
                     "largest_cluster", "n_components", "percolation", "d_f",
                     "crushed", "min_sep", "h_mean", "pe_energy", "r_nn")}
    last: dict = {}

    def sample(timestep, phase):
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, copy=True)
        L = float(snap.configuration.box[0])
        pairs, h_pairs, r_all = contacts(pos, L, r_bond)
        top = topology(n, pairs)
        phi = n * math.pi / 6.0 / L ** 3
        crushed = float((h_pairs < 1.0).sum() / max(len(h_pairs), 1))   # r<d ⟺ h<0
        perc = percolates(pos, L, pairs) if len(pairs) else 0.0
        df = fractal_dimension(pos, L, rng) if top["largest_cluster"] > 0.1 else float("nan")
        row = dict(t=timestep * dt, L=L, phi=phi, percolation=perc, d_f=df,
                   crushed=crushed, min_sep=float(r_all.min()),
                   h_mean=float(h_pairs.mean() - 1.0) if len(h_pairs) else float("nan"),
                   pe_energy=pe_energy(), r_nn=r_nn_mean(),
                   **{k: top[k] for k in ("z", "loops", "dangling", "isolated",
                                          "largest_cluster", "n_components")})
        for k, v in row.items():
            series[k].append(float(v))
        last.clear()
        last.update(row, pos=pos, pairs=pairs, L=L)
        return {k: row[k] for k in ("z", "loops", "phi", "percolation", "crushed",
                                    "min_sep", "pe_energy", "r_nn")}

    def finalize(_res):
        pos, L = last["pos"], last["L"]
        _, _, r_all = contacts(pos, L, r_bond)
        mid, g = rdf(r_all, n, L, min(L / 2.0, 4.0))
        # * The first peak is measured **separately and finely.** The g(r) above has
        #   0.016 d bins and cannot resolve the secondary minimum (h_min=0.0076 d),
        #   which would make the comparison against the predicted 1.00759 meaningless.
        fine_r, fine_g = rdf(r_all, n, L, 1.10, nbins=400)     # bin width 5e-4 d
        peak = float(fine_r[int(np.argmax(fine_g))]) if np.any(fine_g > 0) else float("nan")
        obs = [
            MET.observable(
                "phi_final", last["phi"], predicted=float(P["phi1"]),
                role="implementation_check", tol_pct=1.0,
                note="did the compression reach the target volume fraction exactly",
                derivation="phi = N(pi/6)d^3/L^3 is a **geometric identity**, so it holds for the "
                           "combination too -- this is not a physical prediction but a wiring "
                           "check on whether BoxResize reached the target box L1 (measured: "
                           "Lx error <1e-6, verify/verify_3d_boxresize.py). A mismatch is a "
                           "compression wiring bug."),
            MET.observable(
                "crushed_bond_fraction", last["crushed"], predicted=0.0,
                sigma=1.0 / max(len(last["pairs"]), 1),
                role="implementation_check",
                note="* the fraction of bonds compression crushed into the primary minimum. Must be 0",
                derivation="the design satisfies eps_max=0.399% < eps_crit=0.703%, so no affine "
                           "compression step can push a bonded pair inside the barrier. That "
                           "threshold was **measured on a single pair** (0.40% holds / 0.80% "
                           "collapses, h*=-0.042). Why it extends to the combination: affine "
                           "scaling applies the same relative strain independently to every "
                           "pair, and the collapse verdict depends only on that pair's h. "
                           "WARNING: a percolated network can stiffen and impede relaxation, "
                           "so a non-zero value reads as **a breakdown of that assumption, "
                           "not a bug** -- reduce eps."),
            MET.observable(
                "rdf_first_peak", peak, predicted=float(1.0 + P["h_min_star"]),
                role="implementation_check", tol_pct=3.0,
                note="the first g(r) peak must be the DLVO secondary-minimum bond distance",
                derivation="the pair table was built directly from U_star, so the two-particle "
                           "free-energy minimum is h_min by definition. In the **dilute limit** "
                           "(phi=0.1 with coordination z <~ 3, so the many-body compression on "
                           "any one bond is small) the mode of g(r) lands there. A shift of up "
                           "to 3% from many-body effects is allowed -- hence tol_pct=3."),
            MET.observable("coordination_number", last["z"], role="measurement",
                           note="coordination number z -- the schematic tree in the figure had 1.905"),
            MET.observable("independent_loops", last["loops"], role="measurement",
                           note="independent loop count (Betti-1). The figure had 0 (a tree)"),
            MET.observable("dangling_fraction", last["dangling"], role="measurement",
                           note="free-end fraction. The figure had 6/21 = 0.286"),
            MET.observable("largest_cluster_fraction", last["largest_cluster"],
                           role="measurement"),
            MET.observable("percolation", last["percolation"], role="measurement",
                           note="number of percolating axes / 3"),
            MET.observable("fractal_dimension", last["d_f"], role="measurement",
                           note="mass-radius scaling. DLCA is ~1.8; a dense mesh tends to 3"),
            MET.observable("min_separation", last["min_sep"], role="measurement",
                           note="minimum centre-to-centre distance [d]. Below 1 means overlap"),
            MET.observable("pe_per_particle", last["pe_energy"], role="measurement",
                           note="<U>/N [kT]. For cross-checking consistency against z*well_depth/2"),
            MET.observable("r_nn_mean", last["r_nn"], role="measurement",
                           note="mean nearest-neighbour distance over the 64 reference particles [d] -- the FROZEN indicator"),
        ]
        arrays = {f"series_{k}": np.asarray(v) for k, v in series.items()}
        arrays["rdf_r"], arrays["rdf_g"] = mid, g
        arrays["rdf_fine_r"], arrays["rdf_fine_g"] = fine_r, fine_g
        arrays["final_positions"] = pos
        arrays["final_pairs"] = (last["pairs"] if len(last["pairs"])
                                else np.zeros((0, 2), dtype=int))
        # * `run.execute` wraps `finalize`'s `extra` as `extra={"result": extra}` and
        #   puts it into `metrics["result"]` **as-is**. So wrapping it in
        #   `{"result": ...}` again here would nest it one level deeper as
        #   `metrics["result"]["result"]` -- confirmed by measurement. Keep it **flat.**
        return {"observables": obs, "arrays": arrays,
                "extra": {
                    "n_particles": n, "L_final": L, "r_bond_star": r_bond,
                    "n_links": int(len(last["pairs"])),
                    "eps_max": float(P["eps_max"]), "eps_crit": float(P["eps_crit"]),
                    "frozen_indicator": "r_nn_mean [d] -- <U>/N is exactly 0 early in "
                                        "aggregation, which false-positives FROZEN, and then "
                                        "finalize is skipped and the observables vanish "
                                        "entirely (the abp-rod precedent)",
                }}

    if n_stages > 0:
        phases = [RUN.Phase("aggregation", n_agg, expect_steady=False, note="phi0 fixed, clusters forming"),
                  RUN.Phase("compression", n_stage * n_stages, expect_steady=False,
                            note=f"isotropic compression in {n_stages} stages"),
                  RUN.Phase("post-relaxation", n_post, expect_steady=True, note="the structure settles at phi1")]
    else:
        phases = [RUN.Phase("relaxation", n_post, expect_steady=True,
                            note="relax the sprouted network at phi1 (no compression)")]
    return RUN.Build(
        sim=sim, forces=forces, n_particles=n, sample=sample,
        pe_per_particle=r_nn_mean, n_eq=0,       # * <r_nn> -- see the comment above
        n_prod=n_agg + n_stage * n_stages + n_post,
        sample_every=int(Nm["sample_every"]), phases=phases,
        gsd_path=(Path(outdir) / "traj_A.gsd") if outdir else None,
        tags=["3d", "dlvo", "gelation", "stage1",
              "compression" if n_stages > 0 else "sprout"],
        physical={"phi0": float(P["phi0"]), "phi1": float(P["phi1"]),
                  "n_stages": n_stages, "r_bond_star": r_bond},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="report only")
    ap.add_argument("--spec", action="store_true", help="save the spec, do not run")
    ap.add_argument("--n", type=int, default=None, help="number of particles (default: the first value in system.yaml)")
    ap.add_argument("--stage-tau", type=float, default=1e-3,
                    help="physical time of one compression stage [tau_B]. Close to tau_esc (0.115) is quasi-static")
    ap.add_argument("--agg-tau", type=float, default=0.2, help="aggregation time before compression [tau_B]")
    ap.add_argument("--post-tau", type=float, default=0.2, help="relaxation time after compression [tau_B]")
    ap.add_argument("--init", choices=("scatter", "sprout"), default="scatter",
                    help="scatter = scatter, aggregate, then compress (DLVO dynamics creates the "
                         "topology) . sprout = sprout from seeds directly at the target phi "
                         "(the topology is imposed)")
    ap.add_argument("--max-coord", type=int, default=4,
                    help="the maximum coordination at which sprout still picks a particle as a parent (smaller is sparser)")
    ap.add_argument("--n-seeds", type=int, default=8,
                    help="number of sprout seeds. Scattering them lets growth fronts meet across the "
                         "boundary and form wrapping loops (= percolation). Measured optimum 8 (N=512)")
    ap.add_argument("--loop-bias", type=int, default=100,
                    help="number of candidate directions tried; the position forming the most bonds "
                         "is chosen. 1 = no bias (close to a tree). Measured at 100: loops "
                         "15 -> 190, z 2.05 -> 2.74")
    ap.add_argument("--dt-scale", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1, help="* must be <65536 (bd-hoomd trap 12)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/network/system.yaml")
    n = args.n or sys_["n_list"][0]
    return emit(sys_, n, args)


if __name__ == "__main__":
    raise SystemExit(main())
