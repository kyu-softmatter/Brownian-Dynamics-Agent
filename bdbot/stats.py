"""시계열 통계 — **표본은 상관되어 있다**는 사실 하나에서 나온 도구들.

이 프로젝트에서 같은 실수를 두 번 했습니다:
  1-A: `⟨x²⟩` 오차막대를 naive SEM으로 → 과소평가. 블록 평균으로 고침.
  1-B: 드리프트 t 검정을 naive SE로 → 전 구간 변화 −0.026%인 런이 t=−3.3 (오탐).
       Sokal 자동창으로 `n_eff`를 구해 SE를 팽창시켜 고침.

`tools/postmortem.py`와 케이스 스크립트가 같은 구현을 쓰도록 여기로 모았습니다.
"""
from __future__ import annotations

import math

import numpy as np


def block_sem(x, n_blocks: int = 20) -> float:
    """블록 평균의 표준오차. 상관된 시계열의 정직한 오차막대."""
    x = np.asarray(x, dtype=float)
    if len(x) < n_blocks * 2:
        n_blocks = max(2, len(x) // 2)
    blocks = np.array_split(x, n_blocks)
    means = np.array([b.mean() for b in blocks])
    return float(means.std(ddof=1) / math.sqrt(len(means)))


def tau_int(y, c: float = 5.0) -> float:
    """적분 자기상관 시간 (표본 단위, Sokal 자동창). `n_eff = n/(2τ+1)`."""
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
    """FFT 자기상관, 겹치는 쌍 개수로 정규화. `trace` (n_t, ...) → (n_t,)

    ★ **표본평균을 빼지 않습니다.** 앵커/트랩 중심으로부터의 변위는 참 평균이 정확히 0이고,
      표본평균을 빼면 `C(t)`가 `O(τ/T_obs)`만큼 체계적으로 빨리 감쇠합니다
      (1-A 실측: 피팅된 τ가 −7.75% vs 안 빼면 −0.26%). skill `bd-physics` §5.1.
      참 평균을 모르는 양이라면 빼야 합니다 — 변위는 아닙니다.
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
    """전반/후반 z + 선형 추세 t. 자기상관 보정 포함, 드리프트 **크기**도 함께 보고.

    유의성만으로 판정하면 잘 수렴한 긴 런이 오탐됩니다 — 크기(`drift_span_rel_pct`)를
    함께 봐야 합니다.
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
