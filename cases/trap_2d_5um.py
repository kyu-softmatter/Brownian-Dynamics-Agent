"""trap-2d-5um -- the first end-to-end case.

Physical system (SI) -> scale table -> separation checks -> non-dimensionalization
-> run -> inversion -> comparison against the analytic solution.

**The shared parts were promoted into `bdbot/`.** What remains here is unique to
this system: the harmonic trap force, the four anchor-displacement observables, the
analytic solution, the equilibrium indicator (anchor displacement), and the plots.
What to promote and what to leave was decided by the comparison table in skill
`bd-physics` section 6.3.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/trap_2d_5um.py              # full run (~3 min)
    $PY cases/trap_2d_5um.py --smoke      # short (~20 s)
    $PY cases/trap_2d_5um.py --report     # report only, does not run
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

from bdbot import Q, checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC  # noqa: E402
from bdbot import sim as SIM, stats as ST, traps as TR  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


# ════════════════════════════════════════════════════════════════════════
# 1. the physical system (SI) -- read from YAML with units attached
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": load_node(raw["particle"]["diameter"]),
        "rho_p": load_node(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": load_node(raw["medium"]["temperature"]),
        "eta": load_node(raw["medium"]["viscosity"]),
        "k_t": load_node(raw["external"]["stiffness"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# 2. the scale table (bd-physics section 0, step 1)
#    * case-specific: the l_k and tau_k the trap creates. The governing timescale is
#      tau_k, not tau_B.
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_: dict, L_box, num: dict) -> SC.ScaleLedger:
    """The ledger. * `dt` and `T_obs` are set here too -- both derive from tau_k, so
    they are part of the ledger.

    main() used to compute dt and T_obs separately and keep them outside the ledger.
    Then neither appeared in the timescale ordering table, which halves the ledger's
    purpose of making a separation violation visible.
    (bdbot.scales.MANDATORY_ROLES).
    """
    d = sys_["d"].value.to("m")
    k = sys_["k_t"].value.to("N/m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma = b["kT"], b["gamma"]
    tau_k = C.relaxation_time(gamma, k)
    dt = C.dt_from_bias(tau_k, num["target_bias"])       # inverted from the bias (bd-physics 1.2)
    T_obs = num["production_tau"] * tau_k

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("l_k", ((kT / k) ** 0.5).to("m"), "sqrt(kT/k) trap fluctuation width")
    lg.add_length("L", L_box.to("m"), "box", role="box")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_k", tau_k, "gamma/k trap relaxation", star=True)
    lg.add_time("tau_B", b["tau_B"], "d^2/D_t diffusion (reference)")
    lg.add_time("T_obs", T_obs, "observation window", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("k_d2", (k * d**2).to("J"), "k*d^2 trap stiffness")
    lg.derived = {"gamma": gamma, "D_t": b["D_t"], "m": b["m"], "kT": kT, "k": k, "d": d,
                  "tau_k": tau_k, "dt": dt, "T_obs": T_obs}
    lg.ref = SC.thermal_reference(
        d, kT, b["tau_B"],
        SC.THERMAL_RATIONALE + " The governing timescale of this system is tau_k -- "
        "tau_B=242s is never realized because the trap catches the particle in 4ms. "
        "Results are rescaled into units of tau_k for reporting.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 3. dimensionless groups + 4. separation checks (bd-physics sections 3, 4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg):
    """The dimensionless groups and the separation checks.

    dt and T_obs are read from the ledger (no longer passed as arguments).
    """
    d, kT, k = lg.derived["d"], lg.derived["kT"], lg.derived["k"]
    tau_k, tau_p = lg.get("times", "tau_k"), lg.get("times", "tau_p")
    tau_B, dt, T_obs = lg.get("times", "tau_B"), lg.get("times", "dt"), lg.get("times", "T_obs")
    l_k, L = lg.get("lengths", "l_k"), lg.get("lengths", "L")

    f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)
    k_star = f(k * d**2 / kT)
    # Where num/den are attached, `NondimSpec.validate()` recomputes the ratio from
    # the ledger and compares. k* is not a ratio of two ledger entries but an energy
    # ratio (k*d^2/kT), so it is expressed via ledger symbols.
    groups = [
        ND.Group("k*", k_star, ("energies", "k_d2"), ("energies", "kT"),
                 "k d^2/kT", "trap vs thermal fluctuation"),
        ND.Group("l_k/d", f(l_k / d), ("lengths", "l_k"), ("lengths", "d"),
                 "1/sqrt(k*)", "fluctuation width vs particle"),
        ND.Group("tau_k/tau_B", f(tau_k / tau_B), ("times", "tau_k"), ("times", "tau_B"),
                 "1/k*", "trap relaxation vs diffusion"),
        ND.Group("dt/tau_k", f(dt / tau_k), ("times", "dt"), ("times", "tau_k"),
                 "", "integration resolution"),
        ND.Group("T_obs/tau_k", f(T_obs / tau_k), ("times", "T_obs"), ("times", "tau_k"),
                 "", "observation window"),
    ]
    checks = [
        C.Check("model", "inertia negligible  tau_p/tau_k", f(tau_p / tau_k), C.GATE, "<=",
                "is overdamped BD valid. Independent of dt (bd-physics section 4)"),
        C.Check("integration", "trap resolved      dt/tau_k", f(dt / tau_k), C.GATE, "<=",
                f"bias ~ (dt/tau)/2 = {C.bias_from_dt(dt, tau_k):.3f}% (bd-physics 1.2)"),
        C.Check("geometry", "fluctuation vs box 2l_k/L", f(2 * l_k / L), 0.5, "<=",
                "minimum-image safety margin (bd-hoomd traps 1, 6)"),
        C.Check("statistics", "observation window T_obs/tau_k", f(T_obs / tau_k), 100.0, ">=",
                "steady-state statistical sufficiency", hard=False),
    ]
    return groups, checks


def report_blocks(sys_, lg, n_steps):
    """The case-specific blocks of the report (the shared frame is bdbot.report.render)."""
    tau_k, dt = lg.get("times", "tau_k"), lg.get("times", "dt")
    T_obs = lg.get("times", "T_obs")
    inp = [R.kv(key, f"{sys_[key].value:~.4gP}", sys_[key].tier, sys_[key].source[:44], val_w=20)
           for key in ("d", "T", "eta", "k_t", "rho_p")]
    inp.append(f"  N      = {sys_['N']}")
    der = [f"  {k_:<6} = {lg.derived[k_].to_compact():~.4gP}"
           for k_ in ("gamma", "D_t", "m", "kT")]
    plan = [
        f"  dt      = {dt.to_compact():~.4gP}   (= {float((dt / tau_k).to('')):.2e} τ_k)",
        f"  T_obs   = {T_obs.to_compact():~.4gP}",
        f"  steps   = {n_steps:,}   × N={sys_['N']}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# 5. the analytic solution (the golden-test basis) -- * case-specific. The
#    soft-repulsion case has none.
# ════════════════════════════════════════════════════════════════════════
def analytic(lg):
    kT, k, gamma = lg.derived["kT"], lg.derived["k"], lg.derived["gamma"]
    tau_k = C.relaxation_time(gamma, k)
    x2 = (kT / k).to("um^2")
    return {
        "x2": x2,                                      # <x^2> per degree of freedom
        "sigma": (x2 ** 0.5).to("nm"),                 # the width of P(x)
        "tau_k": tau_k,
        "f_c": (1 / (2 * math.pi * tau_k)).to("Hz"),   # PSD corner
        "S0": (4 * x2 * tau_k).to("um^2/Hz"),          # the PSD f->0 limit
    }


# ════════════════════════════════════════════════════════════════════════
# 6. L4 -- build the system from the spec alone (bdbot.run runs it and judges)
#    * case-specific: a fixed harmonic trap plus 4 observables from the anchor
#      displacement
# ════════════════════════════════════════════════════════════════════════
def analytic_star(k_star: float) -> dict:
    """The closed form in dimensionless units (sigma=d, E=kT, tau=tau_B). The same
    expressions as `analytic(lg)` -- only k* is needed.

    If the two places (`analytic(lg)` for the plots and `finalize` for the observable
    comparison) each built their own values they could diverge, so they are collected
    here.
    """
    return {"x2_star": 1.0 / k_star, "tau_k_star": 1.0 / k_star,
            "fc_star": k_star / (2 * math.pi), "S0_star": 4.0 / k_star ** 2}


@RUN.builder("trap-2d-5um")
def build(spec, outdir=None) -> RUN.Build:
    """Spec -> HOOMD system. The trap is `bdbot.traps.make_trap` (no velocity and no
    drive = fixed).

    * This case used to carry its own `HarmonicTrap` -- that class is exactly what
      `traps.py` records as "promoted after three appearances". The minimum-image and
      tag-indexing conventions are identical (the same `period` computation as
      `minimum_image` in bdbot/sim.py), so swapping it in here does not change the
      result -- verified by re-running and comparing observables to 15 decimal
      places.
    """
    P, Nm = spec.params, spec.numerics
    N, L_star, k_star = int(P["N"]), float(P["L_star"]), float(P["k_star"])
    dt_star = float(Nm["dt_star"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])
    seed = int(Nm["seed"])

    n_side = int(math.ceil(math.sqrt(N)))
    a = L_star / n_side
    pos0 = np.array([[(i % n_side + .5) * a - L_star / 2,
                      (i // n_side + .5) * a - L_star / 2, 0.0] for i in range(N)])
    sim = SIM.make_sim(SIM.frame_2d(pos0, L_star), seed=seed)

    trap = TR.make_trap(k_star, pos0, L_star, dt_star=dt_star)
    SIM.attach_brownian(sim, dt_star, [trap])
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, max(1, n_prod // 200))

    def pe_pp():
        return float(np.array(trap.energies).sum()) / N

    # -- sample accumulators -- keeps the same memory footprint as the original
    #    implementation.
    #   If `sample()` hands the full N-particle array back in `cols` every sample,
    #   RUN.execute stores all of it in observables.npz (stacked as (n_samp,N,2)).
    #   At 1000 particles x 2e4 samples that is 148 MB -- the original was 448 KB
    #   (derived quantities only). So P(x) and <x^2> accumulate immediately in this
    #   closure, and only the subset needed for C(t) and the PSD (n_trace of them) is
    #   filled into a pre-allocated array -- neither passes through `cols`, so
    #   neither is stored.
    n_samp = n_prod // sample_every if sample_every else 0
    n_trace = min(250, N)
    trace = np.empty((max(n_samp, 1), n_trace, 2), dtype=np.float32)
    sum_x2 = np.zeros(2)
    per_sample_x2 = np.empty(max(n_samp, 1))         # for the block SEM (<d^2> over all N)
    hist_edges = np.linspace(-6 / math.sqrt(k_star), 6 / math.sqrt(k_star), 121)
    hist = np.zeros(len(hist_edges) - 1)
    i_sample = [0]

    def sample(timestep, phase):
        dxy = trap.displacement(sim.state, timestep)[:, :2]
        i = i_sample[0]
        if i < n_samp:
            trace[i] = dxy[:n_trace]
        sum_x2[:] += (dxy ** 2).mean(axis=0)
        if i < n_samp:
            per_sample_x2[i] = float((dxy ** 2).mean())
        hist[:] += np.histogram(dxy.ravel(), bins=hist_edges)[0]
        i_sample[0] = i + 1
        return {}

    def finalize(cols):
        from scipy import optimize, signal

        tr = trace[:i_sample[0]]                    # for C(t) and the PSD (the original used this subset too)
        dt_sample_star = dt_star * sample_every

        x2 = float(sum_x2.sum() / (2 * i_sample[0]))  # <x^2> per degree of freedom, over all N
        x2_sem = ST.block_sem(per_sample_x2[:i_sample[0]], 20)

        # P(x) -- the histogram over all N and all samples (accumulated)
        centers = 0.5 * (hist_edges[1:] + hist_edges[:-1])
        px = hist / (hist.sum() * (hist_edges[1] - hist_edges[0]))

        # C(t) -- the sample mean is NOT subtracted (bd-physics section 5.1)
        ac = ST.autocorr_unbiased(tr)
        t = np.arange(len(ac)) * dt_sample_star
        tau_guess = 1.0 / k_star
        fit_n = min(max(10, int(3 * tau_guess / dt_sample_star)), len(ac))
        popt, _ = optimize.curve_fit(lambda tt, A, tau: A * np.exp(-tt / tau),
                                     t[:fit_n], ac[:fit_n],
                                     p0=[ac[0], tau_guess], maxfev=20000)
        tau_fit = float(abs(popt[1]))

        # PSD (one-sided density; ∫S df = variance)
        x = tr[:, :, 0].astype(np.float64).T
        fs = 1.0 / dt_sample_star
        nper = min(len(tr), 4096)
        f, S = signal.welch(x, fs=fs, nperseg=nper, axis=-1, detrend="constant")
        S = S.mean(axis=0)
        psd_f, psd = f[1:], S[1:]

        def lorentz(ff, S0, fc):
            return S0 / (1 + (ff / fc) ** 2)

        try:
            popt, _ = optimize.curve_fit(lorentz, psd_f, psd,
                                         p0=[psd[0], 1.0 / (2 * math.pi / k_star)],
                                         maxfev=20000)
            psd_S0, psd_fc = float(popt[0]), float(abs(popt[1]))
        except Exception as e:
            psd_S0, psd_fc = float("nan"), float("nan")
            print(f"    (PSD fit failed: {e})")

        # -- comparison against the analytic solution
        #    (role: implementation_check, scope: module -- bd-physics section 7.5)
        #   A single trap plus BD is the closed form of exactly the model the
        #   simulation solves. There is no combination assumption.
        ana = analytic_star(k_star)
        # (display name, measured*, predicted*, L, T, unit). The display names and
        # units are kept exactly as the original implementation had them -- the
        # post-mortem and the report are read by a human through these strings.
        rows = [
            ("⟨x²⟩", x2, ana["x2_star"], 2, 0, "µm²"),
            ("σ = √⟨x²⟩", math.sqrt(x2), math.sqrt(ana["x2_star"]), 1, 0, "nm"),
            ("tau (from the C(t) fit)", tau_fit, ana["tau_k_star"], 0, 1, "ms"),
            ("f_c (from the PSD fit)", psd_fc, ana["fc_star"], 0, -1, "Hz"),
            ("S(0) (from the PSD fit)", psd_S0, ana["S0_star"], 2, 1, "µm²/Hz"),
        ]
        obs = []
        for name, meas_star, pred_star, L, T, unit in rows:
            meas = float(spec.physical(meas_star, L=L, T=T).to(unit).magnitude)
            pred = float(spec.physical(pred_star, L=L, T=T).to(unit).magnitude)
            obs.append(MET.observable(
                name, meas, pred, unit, "analytic", role="implementation_check",
                scope="module", tol_pct=5.0,
                note="a single harmonic trap plus BD -- the closed form of an OU process (no combination assumption)"))

        bias_pct = 50.0 * dt_star * k_star     # bias ≈ (dt/τ_k)/2, τ_k* = 1/k*
        return {"observables": obs,
                "extra": {"x2_star": x2, "x2_sem_pct": 100 * x2_sem / x2 if x2 else None,
                          "tau_fit_star": tau_fit, "psd_fc_star": psd_fc,
                          "psd_S0_star": psd_S0, "bias_predicted_pct": bias_pct,
                          "stat_target_pct": 0.5},
                "arrays": {"px_centers": centers, "px": px, "ac_t": t, "ac": ac,
                          "psd_f": psd_f, "psd": psd}}

    return RUN.Build(
        sim=sim, forces=[trap], n_particles=N,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq, n_prod=n_prod, sample_every=sample_every,
        tags=["2D", "harmonic_trap", "newtonian", "single_particle", "no_pair_interaction"],
        physical={"N": N, "k_star": k_star, "L_star": L_star},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="write the L3 spec to specs/<run_id>.json and exit (does not run)")
    ap.add_argument("--force", action="store_true", help="re-run even if a result with the same run_id exists")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/trap-2d-5um/system.yaml")
    num = sys_["numerics"]
    if args.smoke:
        sys_["N"], num = 200, dict(num, production_tau=120, equilibration_tau=10)

    # box: lattice spacing = 1d, overwhelmingly larger than the fluctuation width (~0.004d)
    n_side = int(math.ceil(math.sqrt(sys_["N"])))
    L_star = float(n_side)
    lg = build_ledger(sys_, Q(L_star, "dimensionless") * sys_["d"].value, num)

    tau_k, dt = lg.get("times", "tau_k"), lg.get("times", "dt")
    tau_B, T_obs = lg.ref["time"][1], lg.get("times", "T_obs")
    n_eq = int(round(float((num["equilibration_tau"] * tau_k / dt).to(""))))
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, int(round(float((tau_k / num["samples_per_tau"] / dt).to("")))))
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(sys_, lg)
    # * The L3 artefact. The `run_id` hash covers only {system, params, numerics},
    #   with physics_only applied (if editing a comment or a source invalidates a run,
    #   content addressing is useless).
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"N": sys_["N"], "L_star": L_star,
                "k_star": float((lg.derived["k"] * lg.derived["d"] ** 2
                                 / lg.derived["kT"]).to(""))},
        numerics={"dt_star": float((dt / tau_B).to("")),
                  "dt_over_tau_k": float((dt / tau_k).to("")),
                  "n_eq": n_eq, "n_prod": n_prod, "sample_every": sample_every,
                  "seed": 20260803},
        nhex=12)
    run_id = spec.run_id()

    # L3 integrity (ledger completeness, is each group really that ratio) -- a
    # different layer from the physics checks.
    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}   run_id={run_id}",
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

    # -- L4 -- read the spec back off disk and run it (bdbot.run; the hash check
    #    fires at that point)
    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    if v["status"] == "skipped":
        return 0
    if v["status"] != RUN.OK:
        return 1
    # * Write report.txt AFTER execute() -- execute() clears outdir when it starts,
    #   so writing it first would delete it (prepare_outdir removes everything except
    #   result.txt and record.json).
    (outdir / "report.txt").write_text(report)

    # -- the analytic comparison -- print what finalize() already computed into
    #    metrics.json as a table (this is the case's central deliverable, so it also
    #    goes to the console and result.txt)
    m = json.loads((outdir / "metrics.json").read_text())
    obs_out = m["observables"]
    res_extra = m.get("result", {})
    lines = ["", "=" * R.W, "RESULT -- inverted to physical units, compared against the analytic solution", "=" * R.W,
             f"{'observable':<26}{'measured':>18}{'analytic':>18}{'error':>10}   verdict"]
    all_ok = True
    for o in obs_out:
        ok = o["err_pct"] is not None and abs(o["err_pct"]) < (o["tol_pct"] or 5.0)
        all_ok &= ok
        lines.append(f"{o['name']:<18}{o['measured']:>18.6g}{o['predicted']:>18.6g}"
                     f"{o['err_pct']:>+9.2f}%   {'✓' if ok else '✗'}")
    lines += ["",
              f"  (expected systematic bias on <x^2> = (dt/tau_k)/2 = "
              f"+{res_extra.get('bias_predicted_pct', float('nan')):.3f}%  — bd-physics §1.2)",
              f"  (statistical error on <x^2> +/-{res_extra.get('x2_sem_pct', float('nan')):.3f}%)",
              "=" * R.W,
              f"VERDICT: {'PASS -- all 4 observables match the analytic solution' if all_ok else 'FAIL'}",
              "=" * R.W]
    result = "\n".join(lines)
    print(result)
    (outdir / "result.txt").write_text(report + "\n" + result)

    ana = analytic(lg)
    make_plots(ana, lg, outdir)
    print("\n".join(RID.list_artifacts(outdir, ROOT)))
    return 0 if all_ok else 1


def make_plots(ana, lg, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

    obs = np.load(outdir / "observables.npz")
    d, tau_B = lg.derived["d"], lg.get("times", "tau_B")
    tau_k = ana["tau_k"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    xc = obs["px_centers"] * float(d.to("nm").magnitude)
    sig = float(ana["sigma"].to("nm").magnitude)
    ax[0, 0].plot(xc, obs["px"] / float(d.to("nm").magnitude), "o", ms=3, label="measured")
    xx = np.linspace(xc.min(), xc.max(), 400)
    ax[0, 0].plot(xx, np.exp(-xx**2 / (2 * sig**2)) / (sig * math.sqrt(2 * math.pi)),
                  "-", lw=2, label=f"Gaussian σ={sig:.2f} nm")
    ax[0, 0].set(xlabel="x [nm]", ylabel="P(x) [1/nm]", title="1. position distribution")
    ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

    ax[0, 1].semilogy(xc, obs["px"] / float(d.to("nm").magnitude), "o", ms=3)
    ax[0, 1].semilogy(xx, np.exp(-xx**2 / (2 * sig**2)) / (sig * math.sqrt(2 * math.pi)),
                      "-", lw=2)
    ax[0, 1].set(xlabel="x [nm]", ylabel="P(x)",
                 title="2. position distribution (log) -- Gaussian into the tails?")
    ax[0, 1].grid(alpha=.3)

    t_ms = obs["ac_t"] * float(tau_B.to("ms").magnitude)
    c_um2 = obs["ac"] * float((d**2).to("um^2").magnitude)
    n = min(len(t_ms), int(6 * float(tau_k.to("ms").magnitude) / max(t_ms[1], 1e-12)))
    ax[1, 0].semilogy(t_ms[1:n], c_um2[1:n], "o", ms=3, label="measured")
    tt = np.linspace(0, t_ms[n - 1], 300)
    ax[1, 0].semilogy(tt, float(ana["x2"].to("um^2").magnitude)
                      * np.exp(-tt / float(tau_k.to("ms").magnitude)),
                      "-", lw=2, label=f"exp(−t/τ), τ={float(tau_k.to('ms').magnitude):.3f} ms")
    ax[1, 0].set(xlabel="t [ms]", ylabel="⟨x(0)x(t)⟩ [µm²]", title="3. position autocorrelation")
    ax[1, 0].legend(); ax[1, 0].grid(alpha=.3)

    f_hz = obs["psd_f"] / float(tau_B.to("s").magnitude)
    s_phys = obs["psd"] * float((d**2 * tau_B).to("um^2*s").magnitude)
    ax[1, 1].loglog(f_hz, s_phys, "-", lw=1, alpha=.7, label="measured")
    S0, fc = float(ana["S0"].to("um^2/Hz").magnitude), float(ana["f_c"].to("Hz").magnitude)
    ax[1, 1].loglog(f_hz, S0 / (1 + (f_hz / fc) ** 2), "-", lw=2,
                    label=f"Lorentzian, $f_c$={fc:.1f} Hz")
    ax[1, 1].axvline(fc, ls="--", c="k", alpha=.5)
    ax[1, 1].set(xlabel="f [Hz]", ylabel="S(f) [µm²/Hz]", title="4. power spectrum")
    ax[1, 1].legend(); ax[1, 1].grid(alpha=.3, which="both")

    for a in ax.ravel():
        a.title.set_fontsize(10)
    fig.suptitle("trap-2d-5um -- measured vs analytic (physical units)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
