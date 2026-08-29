"""abp-rod-2d-run-flip -- a run-and-tumble ellipsoid.

Physical system (SI) -> scale table -> separation checks ->
non-dimensionalization -> run -> inversion -> comparison against the analytic
solution.

The shared parts are in `bdbot/`. What remains here is unique to this system:
ellipsoidal friction (Perrin), a custom run-and-tumble updater, the MSD and MSAD
observables, and two analytic predictions:

    orientational correlation   <cos dTheta(t)> = exp(-t/tau_eff),
                                1/tau_eff = D_r + 1/tau_tumble
    effective diffusion         D_eff = D_bar + v^2*tau_eff/2   (2D)

Both follow from assuming that tumbling and thermal rotational diffusion are
**independent processes** -- which is what this run verifies.

WARNING: **this case is kept as a negative example of case design.** All five of
its predictions are `implementation_check` -- they follow from the model that was
implemented -- so it validates the code and discovers nothing. See
docs/02-verification.md section 4 and bd-physics section 7.5.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/abp_rod_2d.py              # full (~140k steps x N=1000)
    $PY cases/abp_rod_2d.py --smoke      # short
    $PY cases/abp_rod_2d.py --report     # report only
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bdbot import Q, checks as C, metrics as MET, physical as PH, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC  # noqa: E402
from bdbot import sim as SIM  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
CASE = ROOT / "intake/abp-rod-2d-run-flip"
f_ = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)

# An orientation jump larger than this between samples is counted as a tumble
# (thermal rotational diffusion alone cannot jump that far within one sample
# interval -- Dr*sample_dt << 1). Computed once at save time and put into the npz,
# so make_plots and a later notebook see the same verdict (the estimate is not
# scattered across several places).
TUMBLE_JUMP_RAD = 1.0


def detect_tumbles(theta, threshold=TUMBLE_JUMP_RAD):
    """theta: (n_t, n_particles). Returns (n_t-1, n_particles) bool -- the per-interval
    tumble estimate.
    """
    dth = np.abs((np.diff(theta, axis=0) + math.pi) % (2 * math.pi) - math.pi)
    return dth > threshold


# ════════════════════════════════════════════════════════════════════════
# 1. the physical system -- read and checked by bdbot.physical (tier, derived_from, shape)
# ════════════════════════════════════════════════════════════════════════
def node(s, *path):
    cur = s.raw
    for k in path:
        cur = cur[k]
    u = cur.get("unit") or "dimensionless"
    return Q(cur["value"], u)


# ════════════════════════════════════════════════════════════════════════
# 2. the scale table -- * case-specific: being an ellipsoid gives three lengths
#    (major, minor, equivalent), and tumbling enters the times
# ════════════════════════════════════════════════════════════════════════
def build_ledger(s, dt=None, T_obs=None):
    """The ledger. * `dt` and `T_obs` go in too (bdbot.scales.MANDATORY_ROLES).

    This system takes dt **directly** from `numerics.dt` -- unlike the trap case
    (inverted from the bias) and the soft-repulsion case (inverted from tau_int),
    here it is a human-chosen value. It has to be in the ledger for its position in
    the timescale ordering to be visible.
    """
    d = node(s, "particle", "diameter").to("m")                 # equivalent-volume sphere diameter (reference)
    semi = node(s, "particle", "semi_axes").to("m")
    kT = node(s, "derived_scales", "kT").to("J") if "kT" in s.raw["derived_scales"] \
        else (Q(1.380649e-23, "J/K") * node(s, "medium", "temperature")).to("J")
    gbar = node(s, "friction", "gamma_bar_2d").to("kg/s")
    gr = node(s, "friction", "gamma_rot_z").to("kg*m^2/s")
    v = node(s, "active", "speed").to("m/s")
    tau_tumble = node(s, "active", "tumble_interval").to("s")
    L = node(s, "geometry", "box_length").to("m")

    D_bar = (kT / gbar).to("m^2/s")
    tau_B = (d**2 / D_bar).to("s")
    D_r = (kT / gr).to("1/s")
    tau_r = (1 / D_r).to("s")
    # * the independence assumption: the decay rates add
    tau_eff = (1 / (D_r + 1 / tau_tumble)).to("s")
    l_p = (v * tau_eff).to("m")
    tau_v = (d / v).to("s")
    rho = node(s, "particle", "density").to("kg/m^3")
    m = (rho * (4 * math.pi / 3) * semi[0] * semi[1] * semi[2]).to("kg")
    tau_p = (m / gbar).to("s")

    lg = SC.ScaleLedger()
    lg.add_length("2a_minor", (2 * semi[1]).to("m"), "minor axis length")
    lg.add_length("d_eq", d, "equivalent-volume sphere diameter (reference)")
    lg.add_length("l_p", l_p, "v*tau_eff persistence length", star=True)
    lg.add_length("2a_major", (2 * semi[0]).to("m"), "major axis length")
    lg.add_length("L", L, "box", role="box")
    lg.add_time("tau_p", tau_p, "m/gamma_bar inertia", role="inertia")
    if dt is not None:
        lg.add_time("dt", dt, "integration step (numerics.dt, human-chosen)", role="dt")
    lg.add_time("tau_v", tau_v, "d_eq/v advection", star=True)
    lg.add_time("tau_eff", tau_eff, "orientational correlation", star=True)
    lg.add_time("tau_tumble", tau_tumble, "tumble interval")
    lg.add_time("tau_r", tau_r, "1/D_r thermal rotational diffusion")
    lg.add_time("tau_B", tau_B, "d_eq^2/D_bar diffusion (reference)")
    if T_obs is not None:
        lg.add_time("T_obs", T_obs, "observation window", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("fa_d", (gbar * v * d).to("J"), "f_a*d_eq self-propulsion work")
    lg.derived = dict(d=d, kT=kT, gbar=gbar, gr=gr, v=v, D_bar=D_bar, tau_B=tau_B,
                      D_r=D_r, tau_r=tau_r, tau_eff=tau_eff, tau_tumble=tau_tumble,
                      tau_v=tau_v, tau_p=tau_p, l_p=l_p, L=L, semi=semi)
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " The reference length is the equivalent-volume sphere "
        "diameter d_eq -- an ellipsoid has separate major and minor axes but the "
        "non-dimensionalization reference must be one length. Translational friction "
        "is the isotropic mean gamma_bar, forced by the BD constraint (the harmonic "
        "mean, so that the long-time MSD is set by D_bar).",
        length_symbol="d_eq")                     # * must match the ledger symbol
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 3. dimensionless groups + 4. separation checks
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(s, lg, N):
    """The dimensionless groups and the separation checks. dt and T_obs come from the
    ledger.
    """
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    Pe = f_(D["v"] * D["d"] / D["D_bar"])
    # Pe **equals** tau_B/tau_v (v*d/D_bar = (d^2/D_bar)/(d/v)) -- cross-checked
    # through the ledger.
    groups = [
        ND.Group("Pe", Pe, ("times", "tau_B"), ("times", "tau_v"),
                 "v d_eq/D̄ = τ_B/τ_v", "advection vs diffusion"),
        ND.Group("D_r*", f_(D["D_r"] * D["tau_B"]), ("times", "tau_B"), ("times", "tau_r"),
                 "D_r*tau_B", "rotational vs translational"),
        ND.Group("l_p/d_eq", f_(D["l_p"] / D["d"]), ("lengths", "l_p"), ("lengths", "d_eq"),
                 "", "persistence length"),
        ND.Group("p", f_(D["semi"][0] / D["semi"][1]), ("lengths", "2a_major"),
                 ("lengths", "2a_minor"), "a_maj/a_min", "aspect ratio"),
        ND.Group("tau_tumble/tau_r", f_(D["tau_tumble"] / D["tau_r"]),
                 ("times", "tau_tumble"), ("times", "tau_r"), "", "* tumbling vs rotational diffusion"),
        ND.Group("L/d_eq", f_(D["L"] / D["d"]), ("lengths", "L"), ("lengths", "d_eq"),
                 "", "box"),
        ND.Group("T_obs/tau_eff", f_(T_obs / D["tau_eff"]), ("times", "T_obs"),
                 ("times", "tau_eff"), "", "observation window"),
        ND.Group("St", f_(D["tau_p"] / D["tau_B"]), ("times", "tau_p"), ("times", "tau_B"),
                 "tau_p/tau_B", "inertia vs diffusion"),
    ]
    ck = [
        C.Check("model", "inertia negligible   tau_p/tau_v", f_(D["tau_p"] / D["tau_v"]), C.GATE, "<=",
                "tau_dyn = the fastest scale of interest = tau_v (advection). BD validity, independent of dt"),
        C.Check("integration", "advection resolved   dt/tau_v", f_(dt / D["tau_v"]), C.GATE, "<=",
                f"moves {100*f_(dt/D['tau_v']):.1f}% of d_eq per step"),
        C.Check("integration", "rotation resolved   dt*D_r", f_(dt * D["D_r"]), C.GATE, "<=",
                "resolves thermal rotational diffusion"),
        C.Check("integration", "tumbling resolved   dt/tau_tumble", f_(dt / D["tau_tumble"]), C.GATE, "<=",
                "the Poisson approximation (p = dt/tau_tumble << 1)"),
        C.Check("geometry", "finite size        l_p/(L/4)", f_(D["l_p"] / (D["L"] / 4)), 1.0, "<=",
                "the persistence length within 1/4 of the box (active artefacts)", hard=False),
        C.Check("statistics", "observation window T_obs/tau_eff", f_(T_obs / D["tau_eff"]), 100.0, ">=",
                f"referenced on the orientational correlation time. {N} independent particles multiply the statistics", hard=False),
    ]
    return groups, ck


def report_blocks(s, lg, n_eq, n_prod, N, sample_every):
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    inp = [
        R.kv("semi_axes", "(1.0,0.25,0.25) µm", 1, "confirmed by the user", val_w=20),
        R.kv("medium", "water @300K", 1, "user-confirmed (medium) + inherited from the trap case (temperature)", val_w=20),
        R.kv("v", f"{D['v'].to('um/s'):~.3gP}", 0, "sketch 'v <= 5 µm/s', the upper bound", val_w=20),
        R.kv("tau_tumble", f"{D['tau_tumble']:~.3gP}", 1, "confirmed by the user", val_w=20),
        R.kv("N", str(N), 3, "*proposed: an independent ensemble", val_w=20),
        R.kv("tumble", "uniform random (2D)", 3, "*proposed (the sketch says 'flip' = 180 deg)", val_w=20),
    ]
    der = [
        f"  gamma_bar(2D) = {D['gbar']:~.4eP}   D_bar = {D['D_bar'].to('um^2/s'):~.4fP}  (Perrin, verified in the sphere limit)",
        f"  γ_r,z   = {D['gr']:~.4eP}   D_r = {D['D_r']:~.4fP}   τ_r = {D['tau_r']:~.4fP}",
        f"  * the real friction is anisotropic (zeta_perp/zeta_par=1.287) but BD is isotropic only -- the short-time MSD anisotropy is lost",
        f"  tau_eff = 1/(D_r + 1/tau_tumble) = {D['tau_eff']:~.4fP}   <- assumes independent processes",
    ]
    plan = [
        f"  dt      = {dt.to('ms'):~.4gP}  = {f_(dt/D['tau_B']):.3e} τ_B",
        f"  T_obs   = {T_obs:~.4gP} = {f_(T_obs/D['tau_eff']):.0f} τ_eff",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}  ({n_prod//sample_every:,} samples)  x N={N}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# 5. the analytic solution -- * case-specific. 2D run-and-tumble plus thermal
#    rotational diffusion
# ════════════════════════════════════════════════════════════════════════
def analytic(lg):
    D = lg.derived
    tau_eff, v, D_bar = D["tau_eff"], D["v"], D["D_bar"]
    return {
        "tau_eff": tau_eff,                              # the <cos dTheta> decay time
        "D_eff": (D_bar + v**2 * tau_eff / 2).to("um^2/s"),   # 2D long-time
        "D_bar": D_bar.to("um^2/s"),                     # short-time
        "ratio": f_((D_bar + v**2 * tau_eff / 2) / D_bar),
    }


def msd_analytic(t, D_bar, v, tau):
    """2D: MSD = 4D̄t + 2v²τ²(t/τ − 1 + e^{−t/τ})"""
    x = t / tau
    return 4 * D_bar * t + 2 * v**2 * tau**2 * (x - 1 + np.exp(-x))


def analytic_star(Pe, Dr_star, tau_tumble_star) -> dict:
    """The same expressions as `analytic(lg)`, in dimensionless units
    (sigma=d_eq, tau=tau_B). D_bar*=1 is true by definition -- D_bar = kT/gamma_bar
    IS the reference diffusion (tau_B := d^2/D_bar), so it is an identity rather than
    a separate derivation.
    """
    tau_eff_star = 1.0 / (Dr_star + 1.0 / tau_tumble_star)
    D_bar_star = 1.0
    D_eff_star = D_bar_star + Pe ** 2 * tau_eff_star / 2
    return {"tau_eff_star": tau_eff_star, "D_bar_star": D_bar_star,
            "D_eff_star": D_eff_star, "ratio": D_eff_star / D_bar_star}


# ════════════════════════════════════════════════════════════════════════
# 6. L4 -- build the system from the spec alone (bdbot.run runs it and judges)
#    * case-specific: the active force plus a thermal rotational-diffusion updater
#      plus a run-and-tumble updater
# ════════════════════════════════════════════════════════════════════════
@RUN.builder("abp-rod-2d-run-flip")
def build(spec, outdir=None) -> RUN.Build:
    """Spec -> HOOMD system.

    * This system has neither a pair interaction nor a trap, so there is no potential
      energy for `pe_per_particle` to point at (the active force is
      non-conservative). Instead it uses `<cos dTheta>` (the orientational
      correlation between samples), which the original implementation already used as
      its equilibrium indicator -- if rotational diffusion dies it pins at 1.0
      (exactly constant -> FROZEN), and when healthy every call sees a different noise
      realization so it can never be exactly constant.
    """
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    N, Pe = int(P["N"]), float(P["Pe"])
    Dr_star = float(P["Dr_star"])
    tau_tumble_star, L_star = float(P["tau_tumble_star"]), float(P["L_star"])
    dt_star, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])

    rng = np.random.default_rng(seed)
    n_side = int(math.ceil(math.sqrt(N)))
    a = L_star / n_side
    pos = np.array([[(i % n_side + .5) * a - L_star / 2,
                     (i // n_side + .5) * a - L_star / 2, 0.0] for i in range(N)])
    fr = SIM.frame_2d(pos, L_star, orientation=True)
    th0 = rng.uniform(0, 2 * math.pi, N)
    fr.particles.orientation = np.c_[np.cos(th0 / 2), np.zeros(N), np.zeros(N),
                                     np.sin(th0 / 2)]
    sim = SIM.make_sim(fr, seed=seed)

    active = md.force.Active(filter=hoomd.filter.All())
    active.active_force["A"] = (Pe, 0.0, 0.0)      # dimensionless: f_a* = Pe (particle local frame)
    active.active_torque["A"] = (0.0, 0.0, 0.0)
    SIM.attach_brownian(sim, dt_star, [active])   # gamma=1, kT=1 (the thermal convention)

    # thermal rotational diffusion -- via an updater, not the integrator (bd-hoomd trap 3)
    sim.operations.updaters.append(
        active.create_diffusion_updater(trigger=hoomd.trigger.Periodic(1),
                                        rotational_diffusion=Dr_star))

    # * run-and-tumble: **uniformly random reorientation** as a Poisson process.
    #   Different from bd-hoomd's run-and-flip snippet -- that is a 180-degree flip,
    #   this is fully random.
    class RunAndTumble(hoomd.custom.Action):
        def __init__(self, rate, dt, seed=7):
            self.p = rate * dt                      # tumble probability per step (must be << 1)
            self.rng = np.random.default_rng(seed)
            self.n_tumbles = 0
            self.n_steps = 0

        def act(self, timestep):
            self.n_steps += 1
            with self._state.cpu_local_snapshot as snap:
                q = np.array(snap.particles.orientation, copy=True)
                hit = self.rng.random(len(q)) < self.p
                k = int(hit.sum())
                if k:
                    self.n_tumbles += k
                    th = self.rng.uniform(0, 2 * math.pi, k)
                    q[hit] = np.c_[np.cos(th / 2), np.zeros(k), np.zeros(k), np.sin(th / 2)]
                    snap.particles.orientation[:] = q

    tumble = RunAndTumble(1.0 / tau_tumble_star, dt_star, seed=seed + 1)
    sim.operations.updaters.append(
        hoomd.update.CustomUpdater(action=tumble, trigger=hoomd.trigger.Periodic(1)))
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, max(1, n_prod // 200))

    # -- sample accumulators -- keeps the same memory footprint as the original.
    #   If sample() hands the full N-particle positions and orientations back in cols
    #   every sample, they stack as (n_samp,N,...) into observables.npz. The original
    #   also computed the MSD and the correlation functions over all N, but **stored
    #   only 8 representative particles** (the trajectory is a first-class
    #   deliverable). To keep that asymmetry, the full arrays live in this closure and
    #   only what goes to the npz is sliced in finalize.
    n_samp = n_prod // sample_every if sample_every else 0
    xy = np.empty((max(n_samp, 1), N, 2), dtype=np.float64)     # unwrapped positions (for the MSD)
    th_arr = np.empty((max(n_samp, 1), N), dtype=np.float64)    # orientation angles
    i_sample = [0]
    prev_theta = [None]

    def pe_pp():
        """<cos dTheta> against the previous call -- this system's only signal that must
        always change.
        """
        sn = sim.state.get_snapshot()
        q = np.array(sn.particles.orientation)
        theta = 2 * np.arctan2(q[:, 3], q[:, 0])
        if prev_theta[0] is None:
            prev_theta[0] = theta
            return 1.0
        val = float(np.cos(theta - prev_theta[0]).mean())
        prev_theta[0] = theta
        return val

    def sample(timestep, phase):
        sn = sim.state.get_snapshot()
        p = np.array(sn.particles.position)[:, :2]
        img = np.array(sn.particles.image)[:, :2]
        i = i_sample[0]
        if i < n_samp:
            xy[i] = p + img * L_star                        # unwrap -- remove the periodic wrapping
            q = np.array(sn.particles.orientation)
            th_arr[i] = 2 * np.arctan2(q[:, 3], q[:, 0])     # rotation angle about z
        i_sample[0] = i + 1
        return {}

    def finalize(cols):
        from scipy import optimize

        n = i_sample[0]
        XY, TH = xy[:n], th_arr[:n]
        dt_sample = dt_star * sample_every
        v_star = Pe

        n_t = len(XY)
        lags = np.unique(np.round(np.logspace(0, math.log10(n_t // 2), 40)).astype(int))
        msd, cth, msad = [], [], []
        for L in lags:
            d = XY[L:] - XY[:-L]
            msd.append(float((d ** 2).sum(axis=2).mean()))
            dth = TH[L:] - TH[:-L]
            cth.append(float(np.cos(dth).mean()))
            # MSAD: tumbling produces jumps up to pi, which makes unwrapping
            #       ambiguous -> fold into [-pi,pi] and label it a "folded MSAD".
            #       The orientational correlation C(t) is the unambiguous observable.
            msad.append(float((((dth + math.pi) % (2 * math.pi) - math.pi) ** 2).mean()))
        t = lags * dt_sample
        msd, cth, msad = np.array(msd), np.array(cth), np.array(msad)

        # * displacement correlation C(t) = <dr(t0).dr(t0+t)> / <dr^2>
        #   BD is overdamped, so **velocity has no physical meaning** (bd-hoomd
        #   trap 5) -> no velocity autocorrelation. Instead, correlate displacements
        #   over a fixed interval delta -- the legitimate discrete-time counterpart,
        #   which for run-and-tumble must decay as exp(-t/tau_eff) (one more
        #   independent check).
        dr = XY[1:] - XY[:-1]                       # (n_t-1, N, 2) consecutive displacements
        dlags = np.unique(np.round(np.logspace(
            0, math.log10(max(2, len(dr) // 3)), 28)).astype(int))
        norm = float((dr * dr).sum(axis=2).mean())
        cdr = np.array([float((dr[L:] * dr[:-L]).sum(axis=2).mean()) / norm
                        for L in dlags])
        disp_corr_t = dlags * dt_sample

        # tau_eff: fit <cos dTheta> to an exponential (the sample mean is NOT
        # subtracted -- the true mean is not 0)
        m = cth > 0.05
        if m.sum() >= 4:
            popt, _ = optimize.curve_fit(lambda tt, A, tau: A * np.exp(-tt / tau),
                                         t[m], cth[m],
                                         p0=[1.0, t[m][len(t[m]) // 2]], maxfev=20000)
            tau_eff_fit, cos_A = float(abs(popt[1])), float(popt[0])
        else:
            tau_eff_fit = cos_A = float("nan")

        # * The MSD is fitted with **the full analytic expression, two free
        #   parameters** (D_bar and tau). Fitting the short-time and long-time regimes
        #   separately put D_bar +37% off, from active contamination.
        def msd_model(tt, D_bar, tau):
            x = tt / tau
            return 4 * D_bar * tt + 2 * v_star ** 2 * tau ** 2 * (x - 1 + np.exp(-x))

        tau0 = tau_eff_fit if np.isfinite(tau_eff_fit) else t[len(t) // 3]
        try:
            popt, _ = optimize.curve_fit(msd_model, t, msd, p0=[1.0, tau0],
                                         sigma=msd, absolute_sigma=False, maxfev=40000)
            D_bar_fit, tau_msd_fit = float(popt[0]), float(abs(popt[1]))
            msd_resid_pct = float(100 * np.sqrt(np.mean(
                ((msd - msd_model(t, *popt)) / msd) ** 2)))
        except Exception as e:
            D_bar_fit = tau_msd_fit = msd_resid_pct = float("nan")
            print(f"    (MSD fit failed: {e})")
        tail = t > 5 * tau0
        D_eff_fit = (float(np.polyfit(t[tail], msd[tail], 1)[0] / 4)
                    if tail.sum() >= 3 else float("nan"))
        head = t < 0.1 * tau0
        D_short_naive = (float(np.polyfit(t[head], msd[head], 1)[0] / 4)
                        if head.sum() >= 3 else float("nan"))
        if head.sum():
            th_ = float(t[head][-1]); xh = th_ / tau0
            head_contam_pct = float(
                100 * 2 * v_star ** 2 * tau0 ** 2 * (xh - 1 + math.exp(-xh))
                / (4 * D_bar_fit * th_))
        else:
            head_contam_pct = float("nan")

        rate_measured = tumble.n_tumbles / (tumble.n_steps * dt_star * N)
        ana = analytic_star(Pe, Dr_star, tau_tumble_star)

        # * Assign the roles honestly (bdbot.metrics.ROLES). All five follow
        #   analytically from the model that was implemented (active force +
        #   rotational-diffusion updater + Poisson tumbling + BD) -- so they are
        #   implementation_check. Agreement = code verification. This case has zero
        #   hypothesis checks, and that is the point of keeping it.
        ROLE, TOL = "implementation_check", 5.0
        rows = [
            ("tau_eff (from the <cos dTheta> fit)", tau_eff_fit, ana["tau_eff_star"], 0, 1, "s", ROLE, TOL,
             "1/tau_eff = D_r + 1/tau_tumble -- it follows because the two processes were **implemented** as independent"),
            ("D_eff (MSD long-time)", D_eff_fit, ana["D_eff_star"], 2, -1, "um^2/s", ROLE, TOL,
             "D_eff = D_bar + v^2*tau_eff/2 -- a consequence of the same model"),
            ("D_bar (full MSD fit)", D_bar_fit, ana["D_bar_star"], 2, -1, "um^2/s", ROLE, TOL,
             "kT/gamma_bar -- does the value put into the integrator come back out (the most circular of the five)"),
            ("tau (full MSD fit)", tau_msd_fit, ana["tau_eff_star"], 0, 1, "s", ROLE, TOL,
             "do the <cos dTheta> and MSD routes give the same tau -- internal consistency"),
            ("tumble rate", rate_measured, 1.0 / tau_tumble_star, 0, -1, "1/s", ROLE, 3.0,
             "does the updater produce the rate it was given"),
        ]
        obs = []
        for name, meas_star, pred_star, L, T, u, role, tol, note in rows:
            meas = float(spec.physical(meas_star, L=L, T=T).to(u).magnitude)
            pred = float(spec.physical(pred_star, L=L, T=T).to(u).magnitude)
            # scope=composite + implementation_check requires a derivation (rule 7).
            # The note IS the derivation, and the fact that all of them amount to
            # "does what I put in come back out" is exactly the circularity -- zero
            # hypothesis checks.
            obs.append(MET.observable(name, meas, pred, u, "analytic_from_model",
                                      role=role, tol_pct=tol, note=note,
                                      scope="composite", derivation=note))

        n_tr = min(8, N)
        traj_theta = TH[:, :n_tr]
        tumble_mask = detect_tumbles(traj_theta)
        cos_series = np.cos(TH[1:] - TH[:-1]).mean(axis=1)
        return {"observables": obs,
                "extra": {"tumble_rate_star": rate_measured, "n_tumbles": tumble.n_tumbles,
                          "D_eff_over_Dbar": D_eff_fit / D_bar_fit,
                          "ratio_model": ana["ratio"], "msd_resid_pct": msd_resid_pct,
                          "D_short_naive_star": D_short_naive,
                          "head_contam_pct": head_contam_pct,
                          "n_particles_traced": n_tr, "n_particles_total": N,
                          "sample_dt_star": dt_sample,
                          "tumble_jump_threshold_rad": TUMBLE_JUMP_RAD},
                "arrays": {"t": t, "msd": msd, "cos_theta": cth, "msad_folded": msad,
                          "cos_theta_series": cos_series,
                          "disp_corr_t": disp_corr_t, "disp_corr": cdr,
                          "traj_xy": XY[:, :n_tr, :], "traj_theta": traj_theta,
                          "tumble_mask": tumble_mask,
                          "final_xy": XY[-1], "final_theta": TH[-1]}}

    return RUN.Build(
        sim=sim, forces=[active], n_particles=N,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq, n_prod=n_prod, sample_every=sample_every,
        tags=["2D", "active", "run_and_tumble", "ellipsoid", "newtonian",
             "no_pair_interaction", "MSD"],
        physical={"N": N, "Pe": Pe, "Dr_star": Dr_star, "tau_tumble_star": tau_tumble_star},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="write the L3 spec to specs/<run_id>.json and exit")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--N", type=int, default=None)
    ap.add_argument("--tobs", type=float, default=None, help="T_obs/τ_eff")
    args = ap.parse_args()

    s = PH.load(CASE)
    if s.errors:
        print(PH.render_check(s))
        return 1
    lg0 = build_ledger(s)                      # tau_eff is needed before dt and T_obs can be set
    D = lg0.derived
    N = args.N or int(node(s, "particle", "count").magnitude)
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    if args.tobs:
        T_obs = args.tobs * D["tau_eff"]
    if args.smoke:
        N, T_obs = min(N, 200), 30 * D["tau_eff"]
    lg = build_ledger(s, dt=dt, T_obs=T_obs)   # the ledger is now complete (every required role filled)
    D = lg.derived

    n_prod = int(round(f_(T_obs / dt)))
    n_eq = max(1, int(round(f_(5 * D["tau_eff"] / dt))))     # equilibrate for 5 orientational correlation times
    sample_every = max(1, n_prod // 3000)
    n_prod = (n_prod // sample_every) * sample_every

    groups, ck = analyze_scales(s, lg, N)
    # * L4 (build()) reads only the spec -- Dr_star, tau_tumble_star and L_star have
    #   to go in here too, or the system cannot be built without re-reading the case
    #   YAML (the L2<->L4 contract is the spec and nothing else).
    spec = ND.NondimSpec(
        case=s.label, system=s.raw, reference=lg.ref, ledger=lg, groups=groups, checks=ck,
        params={"N": N, "Pe": groups[0].value,
                "Dr_star": f_(D["D_r"] * D["tau_B"]),
                "tau_tumble_star": f_(D["tau_tumble"] / D["tau_B"]),
                "L_star": f_(D["L"] / D["d"])},
        numerics={"dt_star": f_(dt / D["tau_B"]), "dt_over_tau_B": f_(dt / D["tau_B"]),
                  "n_eq": n_eq, "n_prod": n_prod, "sample_every": sample_every,
                  "seed": 20260804},
        tag="smoke" if args.smoke else None, nhex=12)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(s, lg, n_eq, n_prod, N, sample_every)
    report, verdict = R.render(
        title=f"DimensionlessReport — {s.label}   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=ck,
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

    # -- the analytic comparison -- print what finalize() already put into metrics.json
    m = json.loads((outdir / "metrics.json").read_text())
    obs_out = m["observables"]
    res_extra = m.get("result", {})
    ana = analytic(lg)

    verdict_o, bad_impl, dev_hypo, meas = MET.judge(obs_out)
    lines = ["", "=" * R.W, "RESULT -- inverted to physical units and compared", "=" * R.W,
             f"{'observable':<32}{'measured':>15}{'predicted':>15}{'error':>9}  role        verdict"]
    for o in obs_out:
        mark = "✓" if o not in bad_impl and o not in dev_hypo else "✗"
        r_short = {"implementation_check": "impl-check", "hypothesis": "hypothesis",
                   "measurement": "measurement"}[o["role"]]
        lines.append(f"{o['name']:<24}{o['measured']:>15.6g}{o['predicted']:>15.6g}"
                     f"{o['err_pct']:>+8.2f}%  {r_short:<10}  {mark}")
    d_um, tau_B_s = f_(D["d"] / Q(1, "um")), f_(D["tau_B"] / Q(1, "s"))
    d_short = res_extra["D_short_naive_star"]
    lines += ["",
              f"  D_eff/D_bar  measured {res_extra['D_eff_over_Dbar']:.3f}"
              f"  vs model {res_extra['ratio_model']:.3f}   (how much activity multiplied the diffusion)",
              f"  full-MSD fit residual RMS {res_extra['msd_resid_pct']:.2f}%",
              (f"  ! naive short-time slope D_bar={d_short*d_um**2/tau_B_s:.4f}"
               f" µm^2/s -- active contamination in that window is {res_extra['head_contam_pct']:.0f}%. Do not use it"
               if math.isfinite(d_short)
               else "  (the sample interval exceeds 0.1 tau_eff, so there is no naive short-time window at all)"),
              "",
              "* what this run confirmed and what it did not",
              "  confirmed    : the implementation matches the intended model (active force,",
              "                 rotational diffusion, tumbling, BD integration). All five",
              "                 predictions derive from that model, so agreement is **code",
              "                 verification.**",
              "  NOT confirmed: the actual physics of this system. Every prediction rests on",
              "                 the simulation's own assumptions, so nothing that could have",
              "                 come out differently was tested -- zero hypothesis checks.",
              "  A hypothesis needs an assumption the simulation does not impose. For example:"
              "    - raise the density to add interactions -> does the dilute D_eff formula break"
              "    - add confinement (walls, a channel) -> does accumulation or upstream swimming appear"
              "    - translational anisotropy -> **impossible in BD** (no HI). Needs another method",
              "=" * R.W,
              f"VERDICT: {verdict_o}",
              "=" * R.W]
    all_ok = not bad_impl
    result = "\n".join(lines)
    print(result)
    (outdir / "result.txt").write_text(report + "\n" + result)

    make_plots(lg, ana, outdir)
    print("\n".join(RID.list_artifacts(outdir, ROOT)))
    return 0 if all_ok else 1


def make_plots(lg, ana, outdir):
    """1. single-particle trajectories  2. MSD  3. two correlation functions
    4. folded MSAD
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

    D = lg.derived
    tau_B = f_(D["tau_B"] / Q(1, "s"))
    d_um = f_(D["d"] / Q(1, "um"))
    te = f_(ana["tau_eff"] / Q(1, "s"))
    Db = f_(ana["D_bar"] / Q(1, "um^2/s"))
    v_um = f_(D["v"] / Q(1, "um/s"))
    z = np.load(outdir / "observables.npz")
    obs = z
    res_extra = json.loads((outdir / "metrics.json").read_text())["result"]
    msd_resid_pct = res_extra["msd_resid_pct"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9.5))

    # -- 1. single-particle trajectories
    tr = z["traj_xy"] * d_um            # (n_t, n_tr, 2) [µm], unwrapped
    tumble_mask = z["tumble_mask"]      # the tumble verdict fixed at save time (detect_tumbles) -- never re-estimated
    dts = float(res_extra["sample_dt_star"]) * tau_B
    cmap = plt.cm.turbo(np.linspace(0.08, 0.92, tr.shape[1]))
    for k in range(tr.shape[1]):
        xy = tr[:, k, :] - tr[0, k, :]
        ax[0, 0].plot(xy[:, 0], xy[:, 1], "-", lw=0.9, color=cmap[k], alpha=.85)
        jump = np.where(tumble_mask[:, k])[0]
        ax[0, 0].plot(xy[jump + 1, 0], xy[jump + 1, 1], ".", ms=3.5,
                      color=cmap[k], alpha=.7)
    ax[0, 0].plot(0, 0, "k+", ms=9, mew=1.6)
    ax[0, 0].set(xlabel="x [µm]", ylabel="y [µm]", aspect="equal",
                 title=f"1. {tr.shape[1]} single-particle trajectories (dots = estimated tumbles)")
    ax[0, 0].grid(alpha=.25)
    ax[0, 0].text(.02, .98, f"ℓ_p = {f_(D['l_p']/Q(1,'um')):.2f} µm\n"
                            f"sample interval {dts*1e3:.0f} ms",
                  transform=ax[0, 0].transAxes, va="top", fontsize=8,
                  bbox=dict(fc="w", alpha=.7, lw=0))

    # ── ② MSD ────────────────────────────────────────────────────────
    t_s = obs["t"] * tau_B
    ax[0, 1].loglog(t_s, obs["msd"] * d_um**2, "o", ms=4.5, label="measured")
    tt = np.logspace(math.log10(t_s[0]), math.log10(t_s[-1]), 300)
    ax[0, 1].loglog(tt, msd_analytic(tt, Db, v_um, te), "-", lw=2,
                    label="analytic (2D run-and-tumble)")
    ax[0, 1].loglog(tt, 4 * Db * tt, "--", lw=1.2, c="gray",
                    label=r"short-time $4\bar{D}t$")
    ax[0, 1].loglog(tt, 4 * f_(ana["D_eff"] / Q(1, "um^2/s")) * tt, ":", lw=1.4,
                    c="tab:red", label=r"long-time $4D_{eff}t$")
    ax[0, 1].axvline(te, ls=":", c="k", alpha=.6)
    ax[0, 1].text(te, ax[0, 1].get_ylim()[0] * 2, r" $\tau_{eff}$", fontsize=9)
    ax[0, 1].set(xlabel="t [s]", ylabel=r"MSD [µm²]",
                 title=f"2. MSD -- residual RMS {msd_resid_pct:.2f}%")
    ax[0, 1].legend(fontsize=8); ax[0, 1].grid(alpha=.3, which="both")

    # -- 3. two correlation functions
    ax[1, 0].semilogy(t_s, np.clip(obs["cos_theta"], 1e-3, None), "o", ms=4.5,
                      label=r"orientation $\langle\cos\Delta\theta\rangle$")
    dct = z["disp_corr_t"] * tau_B
    ax[1, 0].semilogy(dct, np.clip(z["disp_corr"], 1e-3, None), "s", ms=4.5,
                      mfc="none", label=r"displacement $\langle\Delta r\cdot\Delta r\rangle$ (normalized)")
    ax[1, 0].semilogy(tt, np.exp(-tt / te), "-", lw=2,
                      label=f"exp(−t/τ_eff), τ={te:.3f} s")
    ax[1, 0].set(xlabel="t [s]", ylabel="correlation (normalized)", ylim=(1e-3, 2),
                 xlim=(0, min(6 * te, t_s[-1])),
                 title="3. correlation functions -- both routes must give the same tau_eff")
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.3)
    ax[1, 0].text(.98, .95, "* BD is overdamped, so velocity has no meaning (trap 5)\n"
                            "   -> a **displacement** correlation, not a velocity autocorrelation",
                  transform=ax[1, 0].transAxes, ha="right", va="top", fontsize=7.5,
                  bbox=dict(fc="w", alpha=.75, lw=0))

    # -- 4. folded MSAD
    ax[1, 1].semilogx(t_s, obs["msad_folded"], "o", ms=4.5, label="measured (folded)")
    ax[1, 1].axhline(math.pi**2 / 3, ls="--", c="gray",
                     label=r"uniform-distribution limit $\pi^2/3$")
    ax[1, 1].axvline(te, ls=":", c="k", alpha=.6)
    ax[1, 1].set(xlabel="t [s]", ylabel=r"folded MSAD [rad$^2$]",
                 title="4. folded MSAD -- tumble jumps make unwrapping ambiguous")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3, which="both")

    for a in ax.ravel():
        a.title.set_fontsize(10)
    fig.suptitle(f"abp-rod-2d -- run-and-tumble ellipsoid (2µm x 500nm, water 300K, "
                 f"Pe={f_(D['v']*D['d']/D['D_bar']):.2f})", fontsize=12.5)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
