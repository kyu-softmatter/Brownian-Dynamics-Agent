"""S7 — 조화 트랩 계의 분석. 측정만 하고 판정하지 않는다.

판정은 `simbot.validate` 가 **제안**하고, 확정은 사람이 한다 (CLAUDE.md §판정).
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
    """조화 트랩 MSD 해석해:  2d(kT/k)(1 - e^{-t/tau}).  축약 단위에서 plateau = 2d."""
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
    """MSD 곡선을 `plateau (1 - exp(-t/tau))` 로 피팅. 축약 단위."""
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
    """★ HOOMD `Brownian` + 조화 트랩의 **정상분포 초과첨도** — 해석적으로 알려져 있다.

    스킴은 AR(1) 이다:  `x_{n+1} = a x_n + sqrt(2 dt) U`,  `a = 1 - dt*`.
    iid 노이즈의 초과첨도 `g_u` 에 대해 정상분포의 초과첨도는
        `g_x = g_u (1-a^2)/(1+a^2)  ~  g_u * dt*`      (dt* << 1)
    HOOMD 노이즈는 **균일분포**이므로 `g_u = -1.2` (findings §2).
        `g_x ~ -1.2 dt*`
    ⇒ `dt* = 5e-3` 에서 첨도 = `3 - 0.006 = 2.994`.

    **정확히 3.000 이 나오면 오히려 이상하다.** 위치 분포는 완화시간에 걸쳐
    CLT 로 Gaussian 에 접근하지만 `dt*` 차수의 잔여 비가우시안성이 남는다.
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
    """평형 위치의 **형태**가 Gaussian 인가. 폭은 P1(`<x^2>`)이 따로 본다.

    ⚠ **독립 프레임만 쓴다.** 프레임 간격이 `2 tau_trap` 보다 짧으면 표본이 상관되어
      KS 의 귀무분포가 틀리고 **거짓 기각**이 난다. 2026-07-28 에 실제로 그랬다:
      마지막 20 프레임(= 1.6 tau_trap 구간)을 50,000 개 독립표본으로 취급해
      p = 0.0000 을 얻었다. 상관을 무시하면 KS 는 항상 기각한다.

    ⚠ 표본을 **측정된 표준편차로 규격화**한다. 폭 차이(EM 편향 0.25 %)까지 KS 로
      잡으면 형태 검정이 아니라 폭 검정이 된다 — 그건 P1 의 일이다.
    """
    frames_per_tau = 1.0 / (frame_interval_steps * dt_star)
    step = max(1, int(math.ceil(decorrelation_tau * frames_per_tau)))
    frames = traj[::step]                                     # 독립 프레임만
    x = frames.reshape(-1).astype(np.float64)
    x = (x - x.mean()) / x.std(ddof=1)                        # 형태만 본다

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
    """시드 앙상블 집계. **시드간 산포가 오차의 정직한 추정치다.**"""

    mean: float
    se: float                 # 시드 평균의 표준오차
    spread: float             # 시드간 표준편차
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
