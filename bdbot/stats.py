"""Time-series statistics -- tools that all follow from one fact: **samples are
correlated.**

The same mistake was made twice in this project:
  case 1: the `<x^2>` error bar used a naive SEM -> underestimated. Fixed with
          block averaging.
  case 2: the drift t-test used a naive SE -> a run whose total change was
          -0.026% came out at t=-3.3 (a false positive). Fixed by obtaining
          `n_eff` from a Sokal automatic window and inflating the SE.

WARNING: block averaging is itself not the end of it. Measured across a velocity
sweep, the block SEM underestimated the true realization-to-realization spread by
**1.09-2.28x** depending on the observable and the velocity -- and two published
conclusions were reversed once a 9-seed ensemble replaced single runs. In a system
that produces discrete stochastic events, size the ensemble to the events.

Collected here so that `tools/postmortem.py` and the case scripts use the same
implementation.
"""
from __future__ import annotations

import math

import numpy as np


def block_sem(x, n_blocks: int = 20) -> float:
    """Standard error of the block means. The honest error bar for a correlated
    time series.
    """
    x = np.asarray(x, dtype=float)
    if len(x) < n_blocks * 2:
        n_blocks = max(2, len(x) // 2)
    blocks = np.array_split(x, n_blocks)
    means = np.array([b.mean() for b in blocks])
    return float(means.std(ddof=1) / math.sqrt(len(means)))


def tau_int(y, c: float = 5.0) -> float:
    """Integrated autocorrelation time (in samples, Sokal automatic window).

    `n_eff = n/(2*tau+1)`.
    """
    y = np.asarray(y, dtype=float)
    n = len(y)
    if n < 8 or y.std() == 0:
        return 0.0
    y = y - y.mean()
    nfft = 1 << (2 * n - 1).bit_length()
    F = np.fft.rfft(y, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    ac /= ac[0]
    tau, k = 0.0, 1
    while k < n:
        tau += ac[k]
        if k >= c * (2 * tau + 1):
            break
        k += 1
    return max(0.0, float(tau))


def n_eff(y, c: float = 5.0) -> float:
    return len(y) / (2 * tau_int(y, c) + 1)


def autocorr_unbiased(trace):
    """FFT autocorrelation, normalized by the overlapping-pair count.

    `trace` (n_t, ...) -> (n_t,)

    * **The sample mean is NOT subtracted.** The displacement from an anchor or a
      trap centre has a true mean of exactly zero, and subtracting the sample mean
      makes `C(t)` decay systematically faster by `O(tau/T_obs)` (measured: the
      fitted tau came out -7.75% versus -0.26% without subtracting). See skill
      `bd-physics` section 5.1.
      For a quantity whose true mean is unknown you must subtract it -- a
      displacement is not such a quantity.
    """
    x = np.asarray(trace, dtype=np.float64)
    n_t = x.shape[0]
    nfft = 1 << (2 * n_t - 1).bit_length()
    F = np.fft.rfft(x, n=nfft, axis=0)
    ac = np.fft.irfft(F * np.conj(F), n=nfft, axis=0)[:n_t]
    shape = (n_t,) + (1,) * (x.ndim - 1)
    ac /= np.arange(n_t, 0, -1).reshape(shape)
    return ac.reshape(n_t, -1).mean(axis=1)


def stationarity(series, steps=None) -> dict:
    """First-half/second-half z plus a linear-trend t, autocorrelation-corrected,
    reporting the drift **magnitude** alongside.

    Judging on significance alone false-positives a long, well-converged run --
    the magnitude (`drift_span_rel_pct`) has to be read with it.
    """
    y = np.asarray(series, dtype=float)
    n = len(y)
    steps = np.arange(n, dtype=float) if steps is None else np.asarray(steps, dtype=float)
    tau = tau_int(y)
    infl = math.sqrt(2 * tau + 1)

    half = n // 2
    a, b = y[:half], y[half:]
    pooled = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) * infl
    z = float((b.mean() - a.mean()) / pooled) if pooled > 0 else 0.0

    x = (steps - steps.mean()) / (steps.std() or 1.0)
    slope, icept = np.polyfit(x, y, 1)
    resid = y - (slope * x + icept)
    se = resid.std(ddof=2) / math.sqrt(n) * infl
    span = float(slope * (x.max() - x.min()))
    mean = float(y.mean())
    return {"equilibrium_z": z,
            "trend_t": float(slope / se) if se > 0 else 0.0,
            "first_half": float(a.mean()), "second_half": float(b.mean()),
            "tau_int_samples": tau, "n_eff": float(n / (2 * tau + 1)),
            "drift_span": span,
            "drift_span_rel_pct": 100 * span / abs(mean) if mean else None}


__all__ = ["block_sem", "tau_int", "n_eff", "autocorr_unbiased", "stationarity"]
