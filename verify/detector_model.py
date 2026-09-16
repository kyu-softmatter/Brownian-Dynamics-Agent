"""The camera layer r2 left in `gaps[1]` -- exposure integration and localisation noise.

A BD trajectory has no detector. A camera has two, and both act on `var(x)`,
which is what the equipartition stiffness is read from:

  * **motion blur** -- a frame is the position *averaged over the exposure*, not
    sampled at an instant. That is a low-pass filter, so it DEFLATES var(x).
  * **localisation noise** -- the centroid of a finite-photon image lands
    `epsilon` away from the true centre. Independent of the motion, so it
    INFLATES var(x) by exactly `epsilon**2`.

This is post-processing over sampled positions, not a change of physics: nothing
here touches the equation of motion. It lives in `verify/` and not in `bdbot/`
because it has appeared once (CLAUDE.md: promote at the second appearance).

Everything is in units of the trap relaxation time `tau` and the thermal
amplitude `sigma = sqrt(kT/k)`, so it applies to any harmonic trap.

    $PY verify/verify_detector_model.py     # the check against known answers
"""
from __future__ import annotations

import math

import numpy as np

__all__ = ["ou_exact", "blur_variance_factor", "blur_deflation_exact",
           "blur_deflation_free_particle", "apply_detector", "lorentzian_fc"]


# ---------------------------------------------------------------- generator
def ou_exact(n: int, dt_over_tau: float, sigma: float, rng, n_paths: int = 1):
    """`n` samples of an Ornstein-Uhlenbeck process, EXACT at any step size.

        x_{k+1} = x_k e^{-a} + sigma sqrt(1 - e^{-2a}) xi,   a = dt/tau

    Exact means the sampled series has the right variance and the right
    autocorrelation for *any* `dt_over_tau` -- there is no Euler-Maruyama bias to
    confuse with the detector bias being measured. That is the whole reason the
    check below uses this instead of a BD run: it isolates one element (rule 7).
    """
    a = float(dt_over_tau)
    rho = math.exp(-a)
    s = sigma * math.sqrt(1.0 - rho * rho)
    xi = rng.standard_normal((n, n_paths))
    x = np.empty((n, n_paths))
    x[0] = sigma * rng.standard_normal(n_paths)          # start in equilibrium
    for k in range(1, n):
        x[k] = rho * x[k - 1] + s * xi[k]
    return x


# ---------------------------------------------------------------- analytics
def blur_variance_factor(u: float) -> float:
    """`var(boxcar average over u=t_exp/tau) / var(x)`, exact for an OU process.

        2 (u - 1 + e^{-u}) / u^2

    -> 1 - u/3 as u -> 0, so the deflation is `D t_exp / 3` and NOT
    `2 D t_exp / 3`. See `blur_deflation_free_particle`.
    """
    u = float(u)
    if u == 0.0:
        return 1.0
    return 2.0 * (u - 1.0 + math.exp(-u)) / (u * u)


def blur_deflation_exact(u: float) -> float:
    """Fractional deflation of var(x) by exposure, exact. `1 - factor`."""
    return 1.0 - blur_variance_factor(u)


def blur_deflation_free_particle(u: float) -> float:
    """`2 D t_exp / 3` expressed as a fraction of `var(x) = kT/k`, i.e. `2u/3`.

    This is the FREE-PARTICLE MSD blur coefficient (Savin & Doyle; the term that
    appears as `-2 D t_exp / 3` in a measured MSD). It is the value AM's r1
    `camera_model` assumption quotes against a trapped bead's variance. Carried
    here so the two can be printed side by side rather than argued about -- for a
    trapped particle the exact coefficient is `u/3`, so this is 2x.
    """
    return 2.0 * float(u) / 3.0


# ----------------------------------------------------------------- detector
def apply_detector(x, dt_over_tau: float, t_exp_over_tau: float,
                   frame_over_tau: float, epsilon_over_sigma: float, rng):
    """Turn a densely sampled trajectory into camera frames.

    `x` is `(n_steps, n_paths)` sampled every `dt_over_tau`. One frame every
    `frame_over_tau`; within it the first `t_exp_over_tau` is averaged (the
    exposure), then `epsilon` of Gaussian noise is added per frame per path.

    Returns `(frames, n_in_exposure)`. Raises if the exposure is not resolved --
    a boxcar over one sample is not an exposure, it is the instant it started,
    and that failure is silent unless it refuses.
    """
    n_exp = int(round(t_exp_over_tau / dt_over_tau))
    n_frame = int(round(frame_over_tau / dt_over_tau))
    if n_exp < 2:
        raise ValueError(
            f"exposure spans {n_exp} sample(s) at dt/tau={dt_over_tau:g}: the "
            f"boxcar would be a point sample. Need dt <= t_exp/2 "
            f"(t_exp/tau={t_exp_over_tau:g}).")
    if n_frame < n_exp:
        raise ValueError(f"frame period ({n_frame} steps) is shorter than the "
                         f"exposure ({n_exp} steps).")
    n_full = (len(x) // n_frame) * n_frame
    blk = x[:n_full].reshape(-1, n_frame, x.shape[1])
    frames = blk[:, :n_exp, :].mean(axis=1)
    if epsilon_over_sigma:
        frames = frames + epsilon_over_sigma * rng.standard_normal(frames.shape)
    return frames, n_exp


# ------------------------------------------------------------------ readout
def lorentzian_fc(x, fs, nperseg: int | None = None):
    """`f_c` from a one-sided Welch PSD, fitted with `S0/(1+(f/f_c)^2)`.

    Deliberately the SAME estimator and the same fit form as
    `cases/trap_2d_5um.py` `finalize()`, so a difference between a detectored and
    an undetectored series is the detector and not a second estimator.
    """
    from scipy import optimize, signal

    xx = np.asarray(x, dtype=np.float64).T                # (n_paths, n_samples)
    nper = min(xx.shape[-1], nperseg or 4096)
    f, S = signal.welch(xx, fs=fs, nperseg=nper, axis=-1, detrend="constant")
    S = S.mean(axis=0) if S.ndim > 1 else S
    f, S = f[1:], S[1:]
    popt, _ = optimize.curve_fit(lambda ff, S0, fc: S0 / (1 + (ff / fc) ** 2),
                                 f, S, p0=[S[0], fs / 20.0], maxfev=40000)
    return float(abs(popt[1])), float(popt[0]), f, S
