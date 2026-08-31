"""S7 — analysis of the harmonic-trap system. It measures only; it does not judge.

The verdict is **proposed** by `simbot.validate`; a human confirms it
(CLAUDE.md §verdicts).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit
from scipy.stats import kstest, norm


# =============================================================================
def msd_model(t: np.ndarray, plateau: float, tau: float) -> np.ndarray:
    """Harmonic-trap MSD solution: 2d(kT/k)(1 - e^{-t/tau}). plateau = 2d reduced."""
    return plateau * (1.0 - np.exp(-t / tau))


@dataclass
class MSDFit:
    plateau: float
    plateau_se: float
    tau: float
    tau_se: float
    r_squared: float
    n_points: int


def fit_msd(lags_tau: np.ndarray, msd: np.ndarray, dim: int) -> MSDFit:
    """Fit an MSD curve to `plateau (1 - exp(-t/tau))`. Reduced units."""
    t = np.asarray(lags_tau, dtype=np.float64)
    y = np.asarray(msd, dtype=np.float64)
    p0 = (2.0 * dim, 1.0)
    popt, pcov = curve_fit(msd_model, t, y, p0=p0, maxfev=20000)
    resid = y - msd_model(t, *popt)
    ss_res = float(np.sum(resid**2))
    ss_tot = float(np.sum((y - y.mean()) ** 2))
    se = np.sqrt(np.diag(pcov))
    return MSDFit(float(popt[0]), float(se[0]), float(popt[1]), float(se[1]),
                  1.0 - ss_res / ss_tot, t.size)


# =============================================================================
def em_uniform_noise_excess_kurtosis(dt_star: float) -> float:
    """★ The **stationary excess kurtosis** of HOOMD `Brownian` + a harmonic trap
    -- it is known analytically.

    The scheme is AR(1):  `x_{n+1} = a x_n + sqrt(2 dt) U`,  `a = 1 - dt*`.
    For iid noise of excess kurtosis `g_u`, the stationary excess kurtosis is
        `g_x = g_u (1-a^2)/(1+a^2)  ~  g_u * dt*`      (dt* << 1)
    HOOMD's noise is **uniform**, so `g_u = -1.2` (findings §2).
        `g_x ~ -1.2 dt*`
    ⇒ at `dt* = 5e-3` the kurtosis is `3 - 0.006 = 2.994`.

    **Getting exactly 3.000 would be the suspicious outcome.** The position
    distribution approaches Gaussian by the CLT over the relaxation time, but a
    residual non-Gaussianity of order `dt*` remains.
    """
    a = 1.0 - dt_star
    return -1.2 * (1.0 - a * a) / (1.0 + a * a)


@dataclass
class DistributionCheck:
    ks_stat: float
    ks_p: float
    kurtosis: float
    kurtosis_predicted: float
    kurtosis_se: float
    n_independent_frames: int
    n_effective_samples: int


def check_position_distribution(traj: np.ndarray, *, dt_star: float,
                                frame_interval_steps: int,
                                decorrelation_tau: float = 2.0
                                ) -> DistributionCheck:
    """Is the **shape** of the equilibrium position Gaussian? The width is P1's
    job (`<x^2>`), separately.

    ⚠ **Only independent frames are used.** If the frame interval is shorter than
      `2 tau_trap` the samples are correlated, KS's null distribution is wrong,
      and it **falsely rejects**. That actually happened on 2026-07-28: the last
      20 frames (a 1.6 tau_trap window) were treated as 50,000 independent
      samples and gave p = 0.0000. Ignore the correlation and KS always rejects.

    ⚠ The sample is **normalized by the measured standard deviation.** Catch the
      width difference too (the 0.25 % EM bias) and it stops being a shape test
      and becomes a width test — which is P1's job.
    """
    frames_per_tau = 1.0 / (frame_interval_steps * dt_star)
    step = max(1, int(math.ceil(decorrelation_tau * frames_per_tau)))
    frames = traj[::step]                                     # independent frames only
    x = frames.reshape(-1).astype(np.float64)
    x = (x - x.mean()) / x.std(ddof=1)                        # look at the shape only

    ks = kstest(x, norm(loc=0.0, scale=1.0).cdf)
    kurt = float(np.mean(x**4))
    return DistributionCheck(
        ks_stat=float(ks.statistic), ks_p=float(ks.pvalue), kurtosis=kurt,
        kurtosis_predicted=3.0 + em_uniform_noise_excess_kurtosis(dt_star),
        kurtosis_se=float(np.sqrt(24.0 / x.size)),
        n_independent_frames=frames.shape[0], n_effective_samples=x.size,
    )


# =============================================================================
@dataclass
class SeedAggregate:
    """Seed-ensemble aggregate. **The seed-to-seed spread is the honest error.**"""

    mean: float
    se: float                 # standard error of the seed mean
    spread: float             # seed-to-seed standard deviation
    n_seeds: int
    values: list[float]


def aggregate_seeds(values: list[float]) -> SeedAggregate:
    a = np.asarray(values, dtype=np.float64)
    if a.size < 2:
        return SeedAggregate(float(a.mean()), float("nan"), float("nan"), a.size,
                             a.tolist())
    return SeedAggregate(float(a.mean()), float(a.std(ddof=1) / np.sqrt(a.size)),
                         float(a.std(ddof=1)), a.size, a.tolist())


def load_run(outdir: Path) -> dict:
    d = np.load(Path(outdir) / "samples.npz")
    return {k: d[k] for k in d.files}
