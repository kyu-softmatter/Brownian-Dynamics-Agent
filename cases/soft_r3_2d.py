"""soft-r3-2d-A-sweep -- the second end-to-end case.

Physical system (SI) -> scale table -> separation checks ->
non-dimensionalization -> run -> inversion -> verification.

**The shared parts were promoted into `bdbot/`.** What remains here is unique to
this system: the r^-3 table potential plus a WCA core, inverting dt from the
closest approach distance, RSA initial placement, the structural observables
(g(r), psi_6, Voronoi), three verification routes, and the plots.

What differs from the trap case: **there is no analytic solution.** Three checks
stand in for it:
  1. dilute limit    g(r) -> exp(-beta*U)   is the potential implemented right (--dilute)
  2. energy consistency  <U>/N  vs  (rho/2)*integral U g(r) 2*pi*r dr
  3. minimum neighbour distance monitoring -- the defence against pair.Table
     trap 11 (force is 0 for r<r_min)

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/soft_r3_2d.py --A 100            # one A (an independent run, so parallelizable)
    $PY cases/soft_r3_2d.py --A 1 --smoke      # short
    $PY cases/soft_r3_2d.py --dilute           # the dilute-limit verification run
    $PY cases/soft_r3_2d.py --A 100 --report   # report only
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
from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC  # noqa: E402
from bdbot import sim as SIM, stats as ST  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
# * The potential numerics (U, U'', r_min) were promoted into bdbot when
#   `trap-drag` started using the same ones. Redefining them here would let the two
#   cases' dt diverge (see bdbot/pairpot.py).
from bdbot.pairpot import HEX_NN, R_WCA, U2_star, U_star, approach_distance  # noqa: E402,F401

ROOT = Path(__file__).resolve().parent.parent
R_TABLE_MIN = 0.5          # the pair.Table lower bound (same as trap-drag). The basis
                           # of the trap-11 defence


# ════════════════════════════════════════════════════════════════════════
# 1. the physical system (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    P = load_node
    r3 = raw["interactions"][0]
    wca = raw["interactions"][1]
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "phi": float(raw["geometry"]["area_fraction"]["value"]),
        "A_list": list(r3["amplitude_A"]["value"]),
        "r_c": P(r3["cutoff"]),
        "wca_eps_kT": float(wca["epsilon"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# 2. the scale table (bd-physics section 0, step 1)
#    Against the trap case: the lengths go from 3 to 5, and tau_int is a function of
#    r rather than a constant.
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, A, N, phi, r_c_star, dt_scale=1.0, T_obs_tau=None) -> SC.ScaleLedger:
    """The ledger. * `dt` and `T_obs` go in too -- so they are visible in the
    timescale ordering.

    `dt` comes from tau_int(r_min) and `T_obs` is a multiple of tau_B. Both derive
    from values in this ledger, so keeping them outside halves the ledger's purpose
    of making a separation violation visible.
    """
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, D_t, tau_B, m = b["kT"], b["gamma"], b["D_t"], b["tau_B"], b["m"]

    a_star = math.sqrt(math.pi / (4 * phi))             # a_mean/d
    L_star = a_star * math.sqrt(N)
    r_min_star, crit, u_rms_rel, state = approach_distance(A, a_star, sys_["wca_eps_kT"])
    # tau_int = gamma/U''(r_min) -- the same structure as the trap's tau_k = gamma/k
    # (bdbot.checks.relaxation_time). The stiffness is a dimensionless value in units
    # of kT/d^2, so it is converted through tau_B here.
    tau_int = (tau_B / float(U2_star(r_min_star, A, sys_["wca_eps_kT"]))).to("s")

    dt = dt_scale * 1e-2 * tau_int                       # the default sits exactly on the hard gate
    T_obs = (T_obs_tau if T_obs_tau is not None
             else float(sys_["numerics"]["production_tau_B"])) * tau_B

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("r_min", r_min_star * d, "closest approach distance", star=True)
    lg.add_length("a_mean", a_star * d, "mean spacing")
    lg.add_length("r_c", r_c_star * d, "cutoff")
    lg.add_length("L", L_star * d, "box", role="box")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_int", tau_int, "gamma/U''(r_min) interaction", star=True)
    lg.add_time("tau_B", tau_B, "d^2/D_t diffusion (reference)")
    lg.add_time("T_obs", T_obs, "observation window", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("U_a", (float(U_star(a_star, A, sys_["wca_eps_kT"])) * kT).to("J"),
                  "U(a_mean) mean-spacing coupling = Gamma*kT")
    # * U(d) = A + eps_WCA -- **not A*kT.** At r=d the WCA core contributes exactly
    #   epsilon (4*eps*(1-1)+eps). The ledger used to label this "= A kT" and it was
    #   off by 1% at A=100. The L3 integrity check (groups.A != U_d/kT = 101) caught
    #   the mislabel.
    lg.add_energy("U_d", (float(U_star(1.0, A, sys_["wca_eps_kT"])) * kT).to("J"),
                  "U(d) contact coupling = (A+eps_WCA)*kT")
    lg.derived = dict(gamma=gamma, D_t=D_t, m=m, kT=kT, d=d, tau_B=tau_B,
                      a_star=a_star, L_star=L_star, r_min_star=r_min_star,
                      crit=crit, u_rms_rel=u_rms_rel, state=state, tau_int=tau_int,
                      dt=dt, T_obs=T_obs)
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " Unlike the trap case, here tau_B really is the "
        "governing timescale, and dt is set from tau_int(r_min).")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 3. dimensionless groups + 4. separation checks (bd-physics sections 3, 4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg, A, phi, r_c_star):
    """The dimensionless groups and the separation checks.

    dt and T_obs are read from the ledger (no longer passed as arguments).
    """
    D = lg.derived
    f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)
    a_star, L_star = D["a_star"], D["L_star"]
    Gamma = float(U_star(a_star, A, sys_["wca_eps_kT"]))
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    tau_p = lg.get("times", "tau_p")

    # * Gamma is the real control parameter, not A alone (bd-physics 6.2).
    #   It is compared against the ledger's U_a/kT, so if the A(d/a)^3 computation is
    #   wrong, validate() catches it.
    groups = [
        ND.Group("Gamma", Gamma, ("energies", "U_a"), ("energies", "kT"),
                 "U(a_mean)/kT", "coupling vs thermal fluctuation *"),
        # A is the amplitude of the r^-3 term (an input) and is **not** a ratio of two
    # ledger entries -- U(d)/kT = A+eps.
        ND.Group("A", A, None, None, "", "the r^-3 amplitude (input, from the sketch)"),
        ND.Group("U(d)/kT", float(U_star(1.0, A, sys_["wca_eps_kT"])),
                 ("energies", "U_d"), ("energies", "kT"), "(A+ε_WCA)", "contact coupling"),
        ND.Group("phi", phi, None, None, "", "packing"),
        ND.Group("a_mean/d", a_star, ("lengths", "a_mean"), ("lengths", "d"), "", "mean spacing"),
        ND.Group("L/d", L_star, ("lengths", "L"), ("lengths", "d"), "", "box size"),
        ND.Group("r_c/d", r_c_star, ("lengths", "r_c"), ("lengths", "d"), "", "cutoff"),
        ND.Group("r_c/a_mean", r_c_star / a_star, ("lengths", "r_c"), ("lengths", "a_mean"),
                 "", "cutoff (in neighbour shells)"),
        ND.Group("dt/tau_int", f(dt / D["tau_int"]), ("times", "dt"), ("times", "tau_int"),
                 "", "integration resolution"),
        ND.Group("T_obs/tau_B", f(T_obs / D["tau_B"]), ("times", "T_obs"), ("times", "tau_B"),
                 "", "observation window"),
        ND.Group("St", f(tau_p / D["tau_B"]), ("times", "tau_p"), ("times", "tau_B"),
                 "tau_p/tau_B", "inertia vs diffusion"),
    ]
    checks = [
        C.Check("model", "inertia negligible   tau_p/tau_int",
              f(tau_p / D["tau_int"]), C.GATE, "<=",
              "is overdamped BD valid. Independent of dt (bd-physics section 4)"),
        C.Check("integration", "interaction resolved dt/tau_int", f(dt / D["tau_int"]), C.GATE, "<=",
              f"tau_int = gamma/U''(r_min={D['r_min_star']:.3f}d), by the {D['crit']} criterion. "
              f"Bias for a linear system ~ {C.bias_from_dt(dt, D['tau_int']):.3f}% -- this is "
              f"nonlinear, so a convergence check is separate"),
        C.Check("geometry", "cutoff             r_c/(L/2)", r_c_star / (L_star / 2), 1.0, "<=",
              "minimum image (bd-hoomd trap 6). A past violation gave +1856%"),
        C.Check("geometry", "core margin        r_table_min/r_min",
              R_TABLE_MIN / D["r_min_star"], 1.0, "<=",
              "pair.Table trap 11: the force is 0 for r<r_min. Is the approach distance above the table's lower bound"),
        C.Check("statistics", "observation window T_obs/tau_B", f(T_obs / D["tau_B"]), 100.0, ">=",
              "statistics for structural relaxation, referenced on tau_B -- defect annealing at strong coupling can be slower",
              hard=False),
    ]
    return groups, checks, Gamma


def report_blocks(sys_, lg, A, phi, N, n_eq, n_prod):
    """The case-specific blocks of the report (the shared frame is bdbot.report.render)."""
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    inp = [R.kv(key, f"{sys_[key].value:~.4gP}", sys_[key].tier, sys_[key].source[:46])
           for key in ("d", "T", "eta", "rho_p")]
    inp += [
        R.kv("A", f"{A}", 0, "sketch 'A = 0.1, 1, 10, 100' (read as dimensionless -- needs confirmation)"),
        R.kv("phi", f"{phi}", 1, "confirmed by the user (not stated in the sketch)"),
        R.kv("N", f"{N}", 1, "sketch says 100 -> raised for minimum image"),
    ]
    der = [f"  {k_:<8} = {D[k_].to_compact():~.4gP}" for k_ in ("gamma", "D_t", "m")]
    der += [
        f"  state estimate = {D['state']}  (Lindemann sigma_bond/a_NN = {D['u_rms_rel']:.4f}, limit 0.15)",
        f"  a_NN (hexagonal prediction) = {HEX_NN*D['a_star']:.4f} d   <- NOT a_mean",
    ]
    plan = [
        f"  dt      = {dt.to_compact():~.4gP}  = {float((dt/D['tau_B']).to('')):.3e} τ_B",
        f"  T_obs   = {T_obs.to_compact():~.4gP}  = {float((T_obs/D['tau_B']).to('')):.1f} τ_B",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}   × N={N}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# 6. execution (dimensionless units)
# ════════════════════════════════════════════════════════════════════════
def rsa_positions(N, L, min_sep, rng):
    """Random sequential adsorption -- starting from a lattice biases the structural result."""
    pos = np.empty((N, 2))
    n = 0
    tries = 0
    while n < N:
        tries += 1
        if tries > 400 * N:
            raise RuntimeError(f"RSA failed: {n}/{N} (min_sep={min_sep} is too large)")
        p = rng.uniform(-L / 2, L / 2, 2)
        if n:
            dr = pos[:n] - p
            dr -= L * np.round(dr / L)
            if (dr**2).sum(axis=1).min() < min_sep**2:
                continue
        pos[n] = p
        n += 1
    return pos


RDF_BINS = 300             # freud g(r) bin count -- an output setting, not physics, so
                           # it does not go into the spec


@RUN.builder("soft-r3-2d-A-sweep")
def build(spec, outdir=None) -> RUN.Build:
    """Spec -> HOOMD system. r^-3 soft repulsion plus a WCA core, RSA initial placement.

    * The three verification routes (energy consistency, hexagonal NN distance, core
      margin) are computed by `finalize` here. This case has no analytic solution, so
      they are a different layer from `checks` (the hard/soft separation checks) --
      the post-hoc guards originally called `post_checks` have no slot yet in
      `bdbot.run`'s checks schema (it gets promoted when a second case needs it), so
      they go into `extra.post_checks` via `Check.as_dict()`.
    """
    import freud
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    N, A = int(P["N"]), float(P["A"])
    phi, r_c_star, eps = float(P["phi"]), float(P["r_c_star"]), float(P["wca_eps"])
    a_star = math.sqrt(math.pi / (4 * phi))
    L_star = a_star * math.sqrt(N)
    r_min_star = approach_distance(A, a_star, eps)[0]
    dt_star, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])

    np_seed, _ = SIM.resolve_seed(seed)
    rng = np.random.default_rng(np_seed)
    pos = rsa_positions(N, L_star, 1.0, rng)          # * starting from a lattice biases the structure
    sim = SIM.make_sim(SIM.frame_2d(pos, L_star), seed=seed)

    cell = md.nlist.Cell(buffer=0.4)
    # the r^-3 tail -- pair.Table. * endpoint=False (trap 10), shifted at the cutoff
    nbins = 1000
    rr = np.linspace(R_TABLE_MIN, r_c_star, nbins, endpoint=False)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_c_star)
    tab.params[("A", "A")] = dict(r_min=R_TABLE_MIN, U=A / rr**3 - A / r_c_star**3,
                                  F=3 * A / rr**4)
    # the excluded-volume core -- a separate force. The first line of defence against
    # trap 11 (force 0 for r<r_min)
    wca = SIM.wca(cell, epsilon=eps, sigma=1.0)                        # trap 4
    SIM.attach_brownian(sim, dt_star, [tab, wca])
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, max(1, n_prod // 200))

    box = freud.box.Box.square(L_star)
    rdf = freud.density.RDF(bins=RDF_BINS, r_max=min(r_c_star, L_star / 2 - 1e-6), r_min=0.3)
    hexatic = freud.order.Hexatic(k=6, weighted=True)
    voro = freud.locality.Voronoi()
    # * Keeps the same memory footprint as the original implementation -- the
    #   coordination-number histogram accumulates immediately in this closure
    #   (`sample()` could return 13 bins each time, but normalizing once as the
    #   original did is symmetric with finalize, so it stays that way).
    coord_hist = np.zeros(13)

    def xy():
        return np.array(sim.state.get_snapshot().particles.position, dtype=float)

    def pe_pp():
        return float(np.array(tab.energies).sum() + np.array(wca.energies).sum()) / N

    def sample(timestep, phase):
        p = xy()
        rdf.compute((box, p), reset=False)
        vn = voro.compute((box, p)).nlist
        counts = np.asarray(vn.neighbor_counts)
        coord_hist[:] += np.bincount(np.clip(counts, 0, 12), minlength=13)
        dists = np.asarray(vn.distances)
        psi6 = float(np.abs(hexatic.compute((box, p), neighbors=vn).particle_order).mean())
        return {"psi6": psi6, "min_sep": float(dists.min()),
                "bond_mean": float(dists.mean()), "bond_std": float(dists.std())}

    def finalize(cols):
        psi6 = float(cols["psi6"].mean())
        psi6_sem = ST.block_sem(cols["psi6"])
        min_sep = float(cols["min_sep"].min())
        bond_mean = float(cols["bond_mean"].mean())
        bond_std = float(cols["bond_std"].mean())
        pe_mean = float(cols["pe"].mean())
        pe_sem = ST.block_sem(cols["pe"])
        rdf_r, rdf_g = np.array(rdf.bin_centers), np.array(rdf.rdf)
        pe_rdf = energy_from_rdf(rdf_r, rdf_g, A, eps, phi, r_c_star)

        # -- verification -- no analytic solution, so via limits and consistency
        obs = []
        # The energy identity is an **exact expression** that holds whatever g(r) is
        # (no approximation) -- so it is an implementation_check for a crystal and for
        # a fluid alike. A mismatch means a g(r) binning bug or a potential-definition
        # bug. That is why the derivation carries over even at scope=composite
        # (rule 7).
        obs.append(MET.observable(
            "energy consistency <U>/N", pe_mean, pe_rdf, "kT", "consistency",
            role="implementation_check", scope="composite", tol_pct=2.0,
            note="consistency: the HOOMD force sum vs (rho/2)*integral U g(r) 2*pi*r dr",
            derivation="<U>/N = (rho/2)*integral U(r) g(r) 2*pi*r dr is an identity that holds "
                       "for any g(r) (it is not a mean-field or dilute approximation) -- it "
                       "carries over unchanged to a crystal or a fluid"))
        # Assuming a perfect hexagonal crystal, by contrast, is an assumption the
        # simulation does not impose -- whether the r^-3 potential actually condenses
        # into that structure is this case's question, so it is a hypothesis.
        a_nn_pred = HEX_NN * a_star
        if psi6 > 0.6:      # a prediction that only means anything for a crystal
            obs.append(MET.observable(
                "hexagonal NN distance", bond_mean, a_nn_pred, "d", "lattice",
                role="hypothesis", tol_pct=2.0,
                note="lattice: a_NN = sqrt(2/sqrt(3))*a_mean (perfect hexagonal). Applied only when psi_6>0.6"))

        # -- post-hoc guards -- did the run stay inside its design assumptions
        post_checks = [
            C.Check("geometry", "table lower-bound margin r_tab_min/min_sep", R_TABLE_MIN / min_sep, 1.0, "<=",
                  f"pair.Table trap 11: the force is 0 for r<{R_TABLE_MIN}d. Measured minimum {min_sep:.3f}d"),
            C.Check("integration", "design r_min respected r_min/min_sep",
                  r_min_star / min_sep, 1.15, "<=",
                  f"dt was set from the local stiffness at r_min={r_min_star:.3f}d. "
                  f"Measured minimum {min_sep:.3f}d -- if it is well inside, revisit dt"),
        ]

        _, _, u_rms_rel, state = approach_distance(A, a_star, eps)
        final_xy = xy()[:, :2]
        post_dicts = [{**c.as_dict("post_run"), "note": c.note} for c in post_checks]
        return {"observables": obs,
                "extra": {"psi6": psi6, "psi6_sem": psi6_sem,
                          "nn_distance_d": bond_mean, "nn_std_rel": bond_std / bond_mean,
                          "min_sep_d": min_sep, "Gamma": float(U_star(a_star, A, eps)),
                          "coord_hist": list(map(float, coord_hist / coord_hist.sum())),
                          "state_predicted": state, "u_rms_rel_einstein": float(u_rms_rel),
                          "pe_mean": pe_mean, "pe_sem": pe_sem, "pe_rdf": pe_rdf,
                          "post_checks": post_dicts,
                          "post_checks_ok": all(c.ok for c in post_checks)},
                "arrays": {"rdf_r": rdf_r, "rdf_g": rdf_g, "final_xy": final_xy,
                          "coord_hist": coord_hist / coord_hist.sum()}}

    return RUN.Build(
        sim=sim, forces=[tab, wca], n_particles=N,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq, n_prod=n_prod, sample_every=sample_every,
        tags=["2D", "soft_repulsion", "r^-3", "WCA_core", "newtonian",
             "pair_interaction", "structure"],
        physical={"N": N, "A": A, "phi": phi, "r_c_star": r_c_star, "L_star": L_star},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# 7. verification -- no analytic solution, so via limits and consistency
# ════════════════════════════════════════════════════════════════════════
def energy_from_rdf(r, g, A, eps, phi, r_c_star):
    """<U>/N = (rho/2) * integral U(r) g(r) 2*pi*r dr   (2D).

    rho = 4*phi/(pi*d^2) -> dimensionless rho* = 4*phi/pi
    """
    rho_star = 4 * phi / math.pi
    m = r < r_c_star
    Ur = U_star(r[m], A, eps) - A / r_c_star**3          # the actual shifted potential
    integ = np.trapezoid(Ur * g[m] * 2 * math.pi * r[m], r[m])
    return 0.5 * rho_star * integ


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", type=float, default=None, help="the r^-3 amplitude (dimensionless)")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="write the L3 spec to specs/<run_id>.json and exit (does not run)")
    ap.add_argument("--dilute", action="store_true", help="the dilute-limit verification run")
    ap.add_argument("--tobs", type=float, default=None, help="override T_obs/tau_B")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--N", type=int, default=None, help="override the particle count")
    ap.add_argument("--samples", type=int, default=400, help="number of samples in the production phase")
    ap.add_argument("--dt-scale", type=float, default=1.0,
                    help="dt multiplier, for the convergence check -- 0.5 halves it")
    ap.add_argument("--rc-shells", type=float, default=5.0,
                    help="r_c = (this value) * a_mean, for the cutoff-convergence check")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/soft-r3-2d-A-sweep/system.yaml")
    num = sys_["numerics"]

    if args.dilute:
        # The dilute-limit check: g(r) -> exp(-beta*U). Per-bin statistics are what
        # matter, so N and T_obs are set large. The first attempt (N=200, 60 tau_B)
        # had fewer than 6 expected pairs for r<1.5d, which made the comparison
        # impossible.
        A, phi, N = 10.0, 0.01, 800
        r_c_star = 8.0
        T_obs_tau = args.tobs or 200.0
        tag = "dilute"
    else:
        if args.A is None:
            ap.error("--A or --dilute is required")
        A, phi, N = args.A, sys_["phi"], sys_["N"]
        a_star = math.sqrt(math.pi / (4 * phi))
        r_c_star = args.rc_shells * a_star
        T_obs_tau = args.tobs or float(num["production_tau_B"])
        tag = f"A{A:g}"
    if args.N:
        N = args.N
    if args.smoke:
        N, T_obs_tau = min(N, 144), min(T_obs_tau, 3.0)
        tag += "-smoke"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"
    if args.rc_shells != 5.0:
        tag += f"-rc{args.rc_shells:g}"

    lg = build_ledger(sys_, A, N, phi, r_c_star, args.dt_scale, T_obs_tau)
    D = lg.derived
    tau_B, tau_int = D["tau_B"], D["tau_int"]
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")

    n_eq = int(round(float((0.2 * T_obs / dt).to(""))))   # 20% of the observation window to equilibration
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, n_prod // args.samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, Gamma = analyze_scales(sys_, lg, A, phi, r_c_star)
    # * `system` goes into the spec. The old spec had no physical system, so changing
    #   d, eta or rho_p left the run_id identical (even at a 16x difference in tau_B)
    #   and it was mistaken for a completed run of a different system
    #   (`verify/verify_l3_spec_gaps.py`, defect 1).
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"A": A, "phi": phi, "N": N, "r_c_star": r_c_star,
                "wca_eps": sys_["wca_eps_kT"], "Gamma": Gamma},
        numerics={"dt_star": float((dt / tau_B).to("")),
                  "dt_over_tau_int": args.dt_scale * 1e-2,
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": 20260803},
        tag=tag, nhex=10)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, A, phi, N, n_eq, n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  A={A}   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)
    if spec.errors:
        print(f"\nx {len(spec.errors)} L3 integrity error(s) -- the non-dimensionalization does not hold.")
        return 1
    if verdict == "FAIL":
        print("\nx a hard separation check failed -- not running.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 spec: {p.relative_to(ROOT)}")
        return 0

    # -- L4 -- read the spec back off disk and run it (bdbot.run; the hash check fires there)
    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    if v["status"] == "skipped":
        return 0
    if v["status"] != RUN.OK:
        return 1
    # * Write report.txt AFTER execute() -- writing it first would let prepare_outdir delete it.
    (outdir / "report.txt").write_text(report)

    # -- verification -- print what finalize() already computed into metrics.json as a table
    m = json.loads((outdir / "metrics.json").read_text())
    obs_out = m["observables"]
    res_extra = m.get("result", {})
    post_checks = res_extra.get("post_checks", [])
    psi6, psi6_sem = res_extra["psi6"], res_extra["psi6_sem"]
    bond_mean = res_extra["nn_distance_d"]
    bond_std = bond_mean * res_extra["nn_std_rel"]
    min_sep = res_extra["min_sep_d"]
    pe_mean, pe_sem = res_extra["pe_mean"], res_extra["pe_sem"]
    coord_hist = res_extra["coord_hist"]

    lines = ["", "=" * R.W, f"RESULT -- {sys_['label']} A={A} (Gamma={Gamma:.4f})", "=" * R.W]
    lines.append(f"{'verification (measured vs predicted)':<38}{'measured':>14}{'predicted':>14}{'diff':>10}   verdict")
    for o in obs_out:
        ok = o["err_pct"] is not None and abs(o["err_pct"]) < (o["tol_pct"] or 2.0)
        lines.append(f"{o['name']:<26}{o['measured']:>14.6g}{o['predicted']:>14.6g}"
                     f"{o['err_pct']:>+9.2f}%   {'✓' if ok else '✗'}   [{o['unit']}]")
        lines.append(f"    {o['note']}")
    lines.append("")
    lines.append(f"{'post-hoc guard':<38}{'value':>14}{'limit':>14}{'margin':>10}   verdict")
    for c in post_checks:
        lines.append(f"{c['name']:<26}{c['value']:>14.4g}{c['limit']:>14.4g}"
                     f"{c['margin']:>9.2f}×   {'✓' if c['ok'] else '✗'}")
        lines.append(f"    {c['note']}")
    all_ok = (all(abs(o["err_pct"]) < (o["tol_pct"] or 2.0) for o in obs_out
                 if o["err_pct"] is not None)
             and res_extra.get("post_checks_ok", True))

    lines += ["", "OBSERVABLES (structure)",
              f"  ⟨U⟩/N        = {pe_mean:.5f} ± {pe_sem:.5f} kT",
              f"  psi_6 (Voronoi-weighted) = {psi6:.4f} +/- {psi6_sem:.4f}",
              f"  NN distance   = {bond_mean:.4f} d   (s.d. {bond_std:.4f} d"
              f" = {100*bond_std/bond_mean:.2f}%)",
              f"  min neighbour = {min_sep:.4f} d   (minimum over all samples)",
              f"  Voronoi coordination distribution: " +
              "  ".join(f"{i}:{coord_hist[i]:.3f}" for i in range(4, 10)
                        if coord_hist[i] > 0.005),
              "",
              f"  Einstein cage prediction sigma_bond/a_NN = {D['u_rms_rel']:.4f}"
              f"  vs measured {bond_std/bond_mean:.4f}"
              f"  ({100*(bond_std/bond_mean - D['u_rms_rel'])/D['u_rms_rel']:+.1f}%)",
              "    note: the absolute u_rms of a 2D crystal diverges logarithmically",
              "      (Mermin-Wagner). What is finite is the relative (NN) fluctuation. The",
              "      Einstein approximation assumes harmonicity and uncorrelated neighbour",
              "      displacements, so the measurement can exceed it by the",
              "      anharmonicity (r^-3 is softer on the outside). Used as a regime",
              "      indicator, not a quantitative prediction.",
              ]

    # inversion (physical units)
    d_um = float(D["d"].to("um").magnitude)
    lines += ["",
              "INVERSION (physical units)",
              f"  a_mean = {D['a_star']*d_um:.3f} µm      L = {D['L_star']*d_um:.1f} µm",
              f"  NN distance = {bond_mean*d_um:.3f} µm      tau_B = {float(tau_B.to('s').magnitude):.2f} s",
              f"  dt = {float(dt.to('ms').magnitude):.4f} ms   "
              f"T_obs = {float(T_obs.to('s').magnitude):.0f} s = {T_obs_tau:.0f} τ_B",
              "=" * R.W,
              f"VERDICT: {'✓ PASS' if all_ok else '✗ FAIL'}",
              "=" * R.W]
    result = "\n".join(lines)
    print(result)
    (outdir / "result.txt").write_text(report + "\n" + result)

    make_plots(sys_, lg, A, phi, r_c_star, Gamma, outdir, args.dilute)
    print("\n".join(RID.list_artifacts(outdir, ROOT)))
    return 0 if all_ok else 1


def make_plots(sys_, lg, A, phi, r_c_star, Gamma, outdir, dilute):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    res = np.load(outdir / "observables.npz")
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    eps = sys_["wca_eps_kT"]

    # 1. g(r)  (+ compare against exp(-beta*U) in the dilute-limit run)
    r, g = res["rdf_r"], res["rdf_g"]
    ax[0, 0].plot(r, g, "-", lw=1.4, label="measured g(r)")
    if dilute:
        Ur = U_star(r, A, eps) - A / r_c_star**3
        ax[0, 0].plot(r, np.exp(-Ur), "--", lw=1.8, label=r"dilute limit $e^{-\beta U}$")
    ax[0, 0].axhline(1, color="k", lw=.5, alpha=.5)
    ax[0, 0].axvline(lg.derived["a_star"], ls=":", c="gray", label=r"$a_{mean}$")
    ax[0, 0].set(xlabel="r / d", ylabel="g(r)",
                 title=f"1. radial distribution function  (A={A:g}, Gamma={Gamma:.3f})",
                 xlim=(0.5, min(r.max(), 6 * lg.derived["a_star"])))
    ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

    # 2. final configuration + Voronoi coordination
    xy = res["final_xy"]
    L = lg.derived["L_star"]
    ax[0, 1].plot(xy[:, 0], xy[:, 1], "o", ms=max(1.5, 260 / math.sqrt(len(xy)) / 4))
    ax[0, 1].set(xlim=(-L / 2, L / 2), ylim=(-L / 2, L / 2), aspect="equal",
                 xlabel="x / d", ylabel="y / d", title="2. final configuration")
    ax[0, 1].grid(alpha=.2)

    # 3. equilibration + the energy time series
    eq = res["eq_trace"]
    n_eq_pts = len(eq)
    ax[1, 0].plot(np.arange(n_eq_pts), eq[:, 1], "-", label="equilibration")
    ax[1, 0].plot(np.linspace(n_eq_pts, n_eq_pts + 20, len(res["pe"])), res["pe"],
                  "-", lw=.8, alpha=.8, label="production")
    ax[1, 0].set(xlabel="segment (20 equilibration + production)", ylabel="⟨U⟩/N [kT]",
                 title="3. potential energy -- equilibrated?")
    ax[1, 0].legend(); ax[1, 0].grid(alpha=.3)

    # 4. psi_6 time series + coordination distribution
    ax[1, 1].plot(res["psi6"], "-", lw=.9)
    ax[1, 1].set(xlabel="sample", ylabel=r"$|\psi_6|$", ylim=(0, 1),
                 title=f"4. hexagonal order  <psi_6>={res['psi6'].mean():.3f}")
    ax[1, 1].grid(alpha=.3)
    a2 = ax[1, 1].twinx()
    ch = res["coord_hist"]
    a2.bar(np.arange(13), ch, alpha=.25, color="tab:orange", width=.6)
    a2.set_ylabel("Voronoi coordination distribution", color="tab:orange")
    a2.set_xlim(-0.5, 12.5)

    for a in ax.ravel():
        a.title.set_fontsize(10)
    fig.suptitle(f"{sys_['label']}  A={A:g}  φ={phi}  N={len(xy)}", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
