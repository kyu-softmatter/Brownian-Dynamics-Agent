"""The S2 prediction engine — analytic solutions and scaling. 0 lines of LLM.

Every function here produces its answer **before the simulation runs.**
S7 compares these values against the measurement.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

from .units import K_B, kT_si, stokes_drag_si, stokes_einstein_D_si


# =============================================================================
# harmonic trap (optical tweezers) — domains A/B
# =============================================================================
@dataclass
class HarmonicTrapPrediction:
    dim: int
    # inputs (SI)
    T_si: float
    eta_si: float
    radius_si: float
    k_si: float
    # derived (SI)
    kT_si: float
    gamma_si: float
    D0_si: float
    sigma_si: float
    tau_trap_si: float
    tau_D_si: float
    l_trap_si: float
    corner_freq_si: float
    var_per_component_si: float
    msd_plateau_si: float
    rms_radial_si: float
    # dimensionless
    k_star_sigma: float          # k sigma^2 / kT  — trap stiffness (sigma-based)
    l_trap_over_sigma: float
    tau_sep: float               # tau_D / tau_trap = k_star_sigma
    notes: list[str] = field(default_factory=list)

    def msd_si(self, t_si: float) -> float:
        """<|r(t)-r(0)|^2> = 2 d (kT/k) (1 - exp(-t/tau_trap))."""
        return self.msd_plateau_si * (1.0 - math.exp(-t_si / self.tau_trap_si))


def harmonic_trap(
    *, T_si: float, eta_si: float, radius_si: float, k_si: float, dim: int = 2
) -> HarmonicTrapPrediction:
    """Overdamped BD prediction for one spherical particle held in a harmonic trap.

    U(r) = 1/2 k r^2,  r^2 = sum over `dim` components

    Analytic solutions (exact, no approximation):
      per-component equipartition  <x^2> = kT/k              (dim-independent)
      radial square                <r^2> = dim * kT/k
      MSD                          <|dr(t)|^2> = 2 dim (kT/k)(1 - e^{-t/tau})
      relaxation time              tau_trap = gamma/k
      corner frequency             f_c = k/(2 pi gamma)
      confinement length           l_trap = sqrt(kT/k)
    """
    kT = kT_si(T_si)
    gamma = stokes_drag_si(eta_si, radius_si)
    D0 = stokes_einstein_D_si(T_si, gamma)
    sigma = 2.0 * radius_si

    tau_trap = gamma / k_si
    tau_D = sigma**2 / D0
    l_trap = math.sqrt(kT / k_si)
    var_1c = kT / k_si                       # per-component <x^2>
    msd_plateau = 2.0 * dim * var_1c         # MSD(t->inf) = 2 <r^2>
    rms_radial = math.sqrt(dim * var_1c)

    k_star_sigma = k_si * sigma**2 / kT

    notes: list[str] = []
    if k_star_sigma > 1e3:
        notes.append(
            f"k*_sigma = {k_star_sigma:.3g} >> 1 — a very strong trap. "
            f"The particle moves only {l_trap/sigma:.2e} of its own diameter. "
            f"Excluded volume is irrelevant."
        )
    if tau_D / tau_trap > 10:
        notes.append(
            f"tau_D/tau_trap = {tau_D/tau_trap:.3g}. Take tau_D as the reference "
            f"time and dt exceeds the relaxation time — tau_trap must be used."
        )
    return HarmonicTrapPrediction(
        dim=dim, T_si=T_si, eta_si=eta_si, radius_si=radius_si, k_si=k_si,
        kT_si=kT, gamma_si=gamma, D0_si=D0, sigma_si=sigma,
        tau_trap_si=tau_trap, tau_D_si=tau_D, l_trap_si=l_trap,
        corner_freq_si=k_si / (2.0 * math.pi * gamma),
        var_per_component_si=var_1c, msd_plateau_si=msd_plateau,
        rms_radial_si=rms_radial,
        k_star_sigma=k_star_sigma, l_trap_over_sigma=l_trap / sigma,
        tau_sep=tau_D / tau_trap, notes=notes,
    )


# =============================================================================
# systematic error of the numerical scheme — predictable, therefore predicted
# =============================================================================
# --- Euler-Maruyama variance bias: the single source of truth is `bdbot.dt` ---
#  ★ Merged 2026-08-29. The same physics existed in two places:
#      simbot.estimators.euler_maruyama_trap_variance_bias   1/(1-dt*/2)-1   exact
#      bdbot.checks.bias_from_dt                             dt*/2           1st order
#    At dt* = 1e-2 that is 0.5025 % vs 0.5000 % -- **they agree** to the precision
#    each one states. So this was never the kind of duplicate that surfaces as a
#    wrong answer; it is the kind that drifts apart silently the day one side gets
#    "improved". The equations therefore live once, in `bdbot/dt.py`, under two
#    explicitly different names (`em_variance_bias` vs `em_variance_bias_linearized`).
#    ⚠ `bdbot.checks` **keeps the first-order form on purpose.**
#      `cases/trap_2d_5um.py:77` picks its `dt` from it and `run_id` is the hash of
#      the spec content, so switching to the exact form would break the hashes of
#      runs that already exist. Who ran which:
#        exact       -> the S2 trap prediction documents and `campaigns/trap_batch.py`
#                       (via `choose_dt(..., target_em_bias=...)` below)
#        linearized  -> `cases/trap_2d_5um.py`, and `bias_from_dt` reporting in
#                       `trap_drag_2d`, `soft_r3_2d`
#      The gap is exactly `b` (measured); `bdbot.dt.em_bias_form_gap` computes it.
#
#  This value is a **known error**, not a defect: in S7 the measured <x^2> is
#  supposed to come out this much high, and its absence signals something else wrong.
#
#  Scheme:  x_{n+1} = x_n (1 - dt*) + sqrt(2 dt*) xi     (units of tau_trap, D* = 1)
#  Steady:  Var* = 2 dt* / (2 dt* - dt*^2) = 1/(1 - dt*/2)
#  Exact:   Var* = 1
#  =>       relative bias = (dt*/2)/(1 - dt*/2)  ~ dt*/2
from bdbot.dt import (em_variance_bias as euler_maruyama_trap_variance_bias,  # noqa: E402
                      dt_star_for_em_bias as dt_star_for_trap_bias)


# =============================================================================
# validity of the overdamped and continuum approximations (the S3 gate)
# =============================================================================
@dataclass
class ValidityChecks:
    tau_inertial_si: float
    inertial_ratio: float        # tau_inertial / tau_process — << 1 for overdamped
    reynolds: float              # << 1 for Stokes
    passed: bool
    failures: list[str] = field(default_factory=list)


def overdamped_validity(
    *, gamma_si: float, mass_si: float, tau_process_si: float,
    velocity_scale_si: float, radius_si: float, eta_si: float, rho_fluid_si: float,
    inertial_tol: float = 1e-2, reynolds_tol: float = 1e-2,
) -> ValidityChecks:
    """Check the premises of BD (overdamped + Stokes)."""
    tau_i = mass_si / gamma_si
    ratio = tau_i / tau_process_si
    Re = rho_fluid_si * velocity_scale_si * radius_si / eta_si
    fails: list[str] = []
    if ratio > inertial_tol:
        fails.append(f"inertial timescale ratio {ratio:.3g} > {inertial_tol:g} — "
                     f"Langevin needs review")
    if Re > reynolds_tol:
        fails.append(f"Reynolds {Re:.3g} > {reynolds_tol:g} — Stokes drag in doubt")
    return ValidityChecks(tau_i, ratio, Re, not fails, fails)


# =============================================================================
# statistical precision — how many samples are needed
# =============================================================================
def samples_for_variance_precision(rel_err: float) -> float:
    """Independent samples needed to hold a variance estimate's relative standard
    error at or below rel_err.

    Variance estimator of a Gaussian sample: SE(s^2)/s^2 = sqrt(2/(n-1)) ≈ sqrt(2/n)
    """
    return 2.0 / rel_err**2


def seeds_for_target_sigma(*, diff: float, se_diff: float, k_current: int,
                           n_sigma: float = 3.0, t_correction: bool = True) -> dict:
    """Seeds needed to resolve an observed difference at `n_sigma`. **Computed
    before running.**

    ★ Skip this calculation and you end up demanding an unachievable test
      (CLAUDE.md §the 4 test rules ④: "do not demand a 3σ rejection where the
      design power cannot produce 3σ"). And raising the seed count until it turns
      significant is **optional stopping**, which ruins the p value — this is for
      **fixing the target `k` in advance**.

    Assumes `SE ∝ 1/√k`, because that is how the standard error of a seed mean
    shrinks.

    Args:
        t_correction: with `True`, convert `n_sigma` into a `t(ν)` quantile.
            When `k` is small the `SE` itself swings a lot, so the normal quantile
            underestimates
            (basis: findings/tolerance-from-a-4-seed-se-is-not-a-3-sigma-test.md).
    """
    if se_diff <= 0.0:
        raise ValueError(f"se_diff = {se_diff} — must be positive")
    if k_current < 2:
        raise ValueError(f"k_current = {k_current} — at least 2")
    from scipy.stats import t as _t

    sigma_now = abs(diff) / se_diff

    def needed(target_quantile: float) -> float:
        #  SE ∝ 1/√k  ⇒  k_needed = k_current · (target/current σ)²
        return k_current * (target_quantile / sigma_now) ** 2

    #  the normal baseline
    k_normal = needed(n_sigma)
    out = {"sigma_now": sigma_now, "k_current": k_current,
           "n_sigma": n_sigma, "k_needed_normal": k_normal}

    if t_correction:
        #  the t quantile depends on ν = k−1 and ν depends on k → iterate to
        #  self-consistency
        two_sided_p = 2.0 * (1.0 - _norm_cdf(n_sigma))
        k = max(k_normal, 3.0)
        for _ in range(60):
            q = float(_t.ppf(1.0 - two_sided_p / 2.0, df=max(k - 1.0, 1.0)))
            k_new = needed(q)
            if abs(k_new - k) < 1e-9:
                k = k_new
                break
            k = k_new
        out.update({"t_quantile": q, "k_needed": k,
                    "two_sided_p": two_sided_p})
    else:
        out["k_needed"] = k_normal
    out["k_needed_int"] = int(math.ceil(out["k_needed"]))
    out["se_diff_target"] = abs(diff) / out.get("t_quantile", n_sigma)
    return out


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def trap_run_length(
    *, n_particles: int, rel_err_target: float, tau_trap_si: float,
    decorrelation_in_tau: float = 2.0,
) -> dict[str, float]:
    """Run length needed with N non-interacting particles in the same trap.

    The particles do not interact, so **one snapshot is N independent samples**
    (the Barakat 2022 approach:
    knowledge/source/papers/2022-barakat-enhanced-dispersion-harmonic-traps.md).
    Independence along time is secured roughly every 2 tau_trap.
    """
    n_needed = samples_for_variance_precision(rel_err_target)
    n_timepoints = max(1.0, n_needed / n_particles)
    t_total_tau = n_timepoints * decorrelation_in_tau
    return {
        "independent_samples_needed": n_needed,
        "independent_timepoints_needed": n_timepoints,
        "t_total_in_tau_trap": t_total_tau,
        "t_total_si": t_total_tau * tau_trap_si,
    }


# =============================================================================
# cost estimate — measured constants from
# knowledge/wiki/findings/local-cpu-parallelism.md
# =============================================================================
THROUGHPUT_PARTICLE_STEPS_PER_S = 6.3e6
EFFICIENCY_BY_K = {1: 1.0, 2: 0.948, 3: 0.926, 4: 0.925, 5: 0.774,
                   6: 0.696, 8: 0.617, 10: 0.547, 12: 0.443}


def estimate_wall_time_s(n_particles: int, n_steps: int, concurrency: int = 1) -> float:
    """Wall time of a single process [s]. Per-process time when k run concurrently.

    ⚠ This model was measured at N >= 500. At small N the fixed per-step overhead
      dominates, so it underestimates. Measure small N with a pilot run.
    """
    eta = EFFICIENCY_BY_K.get(concurrency)
    if eta is None:
        ks = sorted(EFFICIENCY_BY_K)
        lo = max([k for k in ks if k <= concurrency], default=ks[0])
        eta = EFFICIENCY_BY_K[lo]
    return n_particles * n_steps / (THROUGHPUT_PARTICLE_STEPS_PER_S * eta)
