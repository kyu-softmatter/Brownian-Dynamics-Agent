"""런타임 가드 — 시뮬레이션이 조용히 틀리는 것을 잡는다.

핵심 원칙: **가드는 계통적으로 이탈할 수 있는 양을 봐야 한다.**
HOOMD `Brownian`의 운동에너지 온도는 매 스텝 목표 분포에서 뽑히므로 계통 이탈이
불가능하다 → 가드로 쓸 수 없다. 배위 온도가 그 자리를 대신한다.
근거: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

import numpy as np


# =============================================================================
# 배위 온도 — BD 의 진짜 온도계
# =============================================================================
def configurational_temperature(forces: np.ndarray, laplacian_U_total: float) -> float:
    """kT_conf = <|grad U|^2> / <laplacian U>   (Rugh / Butler 배위 온도)

    위치와 힘만 쓴다 — 속도를 쓰지 않으므로 BD 에서 유효하다.

    Args:
        forces: (N, 3) 또는 (N, d) 힘 배열. `F = -grad U` 이므로 `|F|^2 = |grad U|^2`
        laplacian_U_total: 입자 전체에 대한 `sum_i laplacian_i U` 의 앙상블 평균.
            포텐셜마다 해석적으로 다르므로 호출자가 공급한다.
            조화 트랩 `U = 1/2 k r^2` (활성축 d개): 입자당 `d*k`

    Returns:
        kT_conf. 입력 kT 와 일치해야 한다.

    ⚠ **순수 조화 트랩에서는 `<x^2> = kT/k` 검사와 대수적으로 동일하다** — 새 정보가 아니다.
      독립적인 검사가 되는 것은 쌍 상호작용이 있을 때다 (grad U 가 이웃에서도 오므로).
    """
    if laplacian_U_total <= 0:
        raise ValueError(f"laplacian_U_total must be > 0, got {laplacian_U_total}")
    return float(np.mean(np.sum(np.asarray(forces, dtype=np.float64) ** 2, axis=1))
                 / laplacian_U_total)


# =============================================================================
# 스텝당 변위 — 폭발 감지
# =============================================================================
@dataclass
class DisplacementReport:
    max_over_sigma: float
    rms_over_sigma: float
    max_over_rms: float
    n_exceeding: int
    passed: bool
    note: str = ""


def check_step_displacements(
    dr: np.ndarray, sigma: float, max_frac: float = 0.10
) -> DisplacementReport:
    """스텝당 변위가 sigma 의 max_frac 을 넘는 입자가 있는가.

    ⚠ **변위 분포를 Gaussian 으로 가정하지 말 것.** HOOMD `Brownian` 의 노이즈는
      균일분포이므로 성분별 `max/sigma_step = sqrt(3) = 1.732` 가 구조적 상한이다
      (Gaussian 이면 상한이 없다). 따라서 "n-sigma 이상이면 이상" 류의 판정은
      이 엔진에서 물리적 의미가 다르다 — **절대 변위로 판정한다.**
      근거: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
    """
    dr = np.asarray(dr, dtype=np.float64)
    mag = np.linalg.norm(dr, axis=1)
    rms = float(np.sqrt(np.mean(mag**2)))
    mx = float(mag.max()) if mag.size else 0.0
    n_bad = int(np.count_nonzero(mag > max_frac * sigma))
    note = ""
    if n_bad:
        note = (f"{n_bad}개 입자가 스텝당 {max_frac:.3g} sigma 초과 "
                f"(최대 {mx/sigma:.4g} sigma) — dt 과대 또는 초기 겹침 의심")
    return DisplacementReport(
        max_over_sigma=mx / sigma, rms_over_sigma=rms / sigma,
        max_over_rms=(mx / rms if rms > 0 else float("nan")),
        n_exceeding=n_bad, passed=(n_bad == 0), note=note,
    )


# =============================================================================
# 유한성 · 경계
# =============================================================================
def check_finite(**arrays: np.ndarray) -> tuple[bool, list[str]]:
    """NaN/Inf 검사. 위반한 배열 이름과 개수를 돌려준다."""
    fails = []
    for name, a in arrays.items():
        a = np.asarray(a)
        n_bad = int(np.count_nonzero(~np.isfinite(a)))
        if n_bad:
            fails.append(f"{name}: {n_bad}개 비유한값")
    return (not fails), fails


def check_inside_box(positions: np.ndarray, box_lengths, dims: int = 3,
                    tol: float = 1e-9) -> tuple[bool, int]:
    """입자가 박스 안에 있는가 (벽 관통 감지). PBC 계에서는 래핑되므로 항상 통과한다."""
    pos = np.asarray(positions, dtype=np.float64)[:, :dims]
    half = np.asarray(box_lengths, dtype=np.float64)[:dims] / 2.0
    outside = np.any(np.abs(pos) > half + tol, axis=1)
    n = int(np.count_nonzero(outside))
    return (n == 0), n


# =============================================================================
# 결합 계 전용 — 결합 길이가 살아 있는가
# =============================================================================
def check_bond_lengths(positions: np.ndarray, bonds, target: float,
                       tol: float = 0.05, dims: int = 3
                       ) -> tuple[bool, dict]:
    """결합 길이가 목표값 근처인가.

    **`check_finite` 만으로는 결합 계의 폭발을 못 잡는다.** `constrain.Distance` 를
    `Brownian` 과 쓰면 결합 길이가 `1.0 -> 5.8e7` 이 되는데 전부 유한하므로
    `check_finite` 는 통과한다. 폭발한 사슬에서 `kappa` 를 재면 그럴듯한 `s^-3` 이
    나올 수 있다.
    근거: findings/dead-end-distance-constraint-with-brownian.md

    Args:
        bonds: (n_bonds, 2) 입자 인덱스
        target: 목표 결합 길이 (축약 단위)
        tol: 허용 상대 오차. 굽힘 측정에서는 신축이 굽힘을 오염시키지 않아야 하므로
            느슨하게 잡지 말 것.

    Returns:
        (ok, 상세). 상세에는 max/mean 상대 편차와 최악 결합 인덱스가 들어간다.
    """
    pos = np.asarray(positions, dtype=np.float64)[:, :dims]
    b = np.asarray(bonds, dtype=int)
    d = np.linalg.norm(pos[b[:, 1]] - pos[b[:, 0]], axis=1)
    rel = np.abs(d - target) / target
    imax = int(np.argmax(rel))
    info = {
        "max_rel_dev": float(rel[imax]),
        "mean_rel_dev": float(rel.mean()),
        "worst_bond": imax,
        "worst_length": float(d[imax]),
        "target": float(target),
        "tol": float(tol),
        "n_violating": int(np.count_nonzero(rel > tol)),
    }
    return bool(rel.max() <= tol), info


# =============================================================================
# 통계량이 요동하는지 — "항등식을 측정으로 착각하는" 실패 방지
# =============================================================================
def assert_statistic_fluctuates(samples, name: str = "statistic",
                                min_rel_std: float = 1e-12) -> None:
    """독립 표본에서 얻은 통계량이 실제로 요동하는지 확인한다.

    요동하지 않는 "측정값"은 측정이 아니라 **산술 항등식**이다.
    2026-07-28 에 실제로 겪었다: 변위에서 평균을 뺀 뒤 교차상관을 재면
    `cross/auto = -1/(n-1)` 이 항등적으로 나오고, 200회 반복의 표준편차가 `6.7e-20`
    이었다. 결과가 그럴듯해서 통과할 뻔했다.
    근거: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md §3

    Raises:
        AssertionError: 상대 표준편차가 min_rel_std 미만이면
    """
    s = np.asarray(samples, dtype=np.float64)
    if s.size < 2:
        raise ValueError("need >= 2 samples")
    scale = max(abs(float(np.mean(s))), 1e-300)
    rel = float(np.std(s, ddof=1)) / scale
    if rel < min_rel_std:
        raise AssertionError(
            f"{name}: 상대 표준편차 {rel:.3e} < {min_rel_std:.3e} — "
            f"통계량이 요동하지 않는다. 측정이 아니라 산술 항등식일 가능성이 크다."
        )
