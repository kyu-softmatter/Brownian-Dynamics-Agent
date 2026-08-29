"""소프트 반발 페어 퍼텐셜의 수치 — `U`, `U''`, 최근접 접근거리 (skill `bd-physics` §6.2).

**두 번 나와서 올렸습니다** (CLAUDE.md 추상화 규칙). 1-B `soft-r3-2d-A-sweep` 가 유일한
사용자였는데 `trap-drag-2d-hex300` 이 같은 페어 퍼텐셜(`A(d/r)³ + WCA`)에 같은 밀도로
같은 육방 격자를 만들면서 `r_min → τ_int → dt` 경로를 그대로 씁니다.

여기 있는 함수 3개는 **`dt` 를 정하는 물리**입니다. 값이 조용히 틀리면 적분이 조용히
틀립니다. 1-B에서 세 번 고쳤고 (`approach_distance` 도크스트링 참조), 그 교정이 한 곳에만
있어야 두 케이스가 갈라지지 않습니다.

⚠️ **`bdbot.interactions` 와 다른 층입니다.** 그쪽은 "무엇을 고를까"의 카탈로그(텍스트
   메타데이터)이고, 여기는 고른 것의 수치입니다.
"""
from __future__ import annotations

import math

import numpy as np

R_WCA = 2 ** (1 / 6)                 # WCA 컷오프 (σ 단위) — LJ 최소점
# 2D 육방격자: 입자당 면적 = (√3/2) a_NN²  ⟹  a_mean ≡ ρ^(-1/2) = √(√3/2)·a_NN
# ★ 최근접이웃 거리는 a_mean 이 아니다. 1-B에서 a_mean으로 잡았다가 스모크 측정에서 교정.
HEX_NN = math.sqrt(2 / math.sqrt(3))          # a_NN / a_mean = 1.07457


def U_star(rs, A, eps=1.0):
    """무차원 퍼텐셜 U/kT, 길이는 d 단위. WCA 코어 + A/r³ (컷오프 시프트는 별도).

    ⚠️ `U_star(1.0, A)` 는 `A` 가 **아니라** `A + eps` 입니다 — r=d 에서 WCA 코어가
       정확히 ε 을 냅니다 (4ε(1−1)+ε). A=100 이면 101 로 1% 어긋납니다.
    """
    rs = np.asarray(rs, dtype=float)
    w = np.where(rs < R_WCA, 4 * eps * (rs**-12.0 - rs**-6.0) + eps, 0.0)
    return w + A / rs**3


def U2_star(rs, A, eps=1.0):
    """U''(r) [kT/d²] — 국소 강성. τ_int = γ/U'' 의 분모."""
    rs = np.asarray(rs, dtype=float)
    w = np.where(rs < R_WCA, 4 * eps * (156 * rs**-14.0 - 42 * rs**-8.0), 0.0)
    return w + 12 * A / rs**5


def approach_distance(A, a_star, eps=1.0, u_max=12.0, lindemann=0.15):
    """최근접 접근거리 r_min* [d]. 두 기준 중 작은 쪽. dt는 여기서 나온다.

    (a) 쌍 기준   U(r) = u_max kT.
        u_max=12 인 이유: 결합 표본이 ~10⁶개(400 표본 × 400 입자 × 6 이웃)라
        Boltzmann 꼬리의 극단값이 βU ≈ ln(10⁶) ≈ 14 근처까지 간다.
        u_max=5로 잡았더니 스모크 측정 최소거리가 그보다 훨씬 안쪽이었다.
    (b) 진동 기준 a_NN − 3 σ_bond   (케이지가 살아 있을 때만 유효)
        ★ 이웃거리는 a_mean 이 아니라 a_NN = 1.07457 a_mean (육방).
        ★ 결합길이 요동 σ_bond = √2 · u₁ (u₁ = 성분당 rms). √2를 빠뜨리기 쉽다.
        판정에 쓰는 Lindemann 지표는 σ_bond/a_NN (2D에서 유한한 것은 상대 요동).

    반환: (r_min*, 어느 기준, Lindemann σ_bond/a_NN, 상태 추정)
    """
    lo, hi = 0.4, 60.0
    for _ in range(200):                       # 이분법
        mid = 0.5 * (lo + hi)
        if float(U_star(mid, A, eps)) > u_max:
            lo = mid
        else:
            hi = mid
    r_pair = 0.5 * (lo + hi)

    a_nn = HEX_NN * a_star
    k1 = 3 * (float(U2_star(a_nn, A, eps)) + (-3 * A / a_nn**4) / a_nn)   # 성분당 케이지 강성
    if k1 <= 0:
        return r_pair, "쌍", float("nan"), "유체"
    sigma_bond = math.sqrt(2.0 / k1)
    lind = sigma_bond / a_nn
    r_cage = a_nn - 3 * sigma_bond
    if lind < lindemann and r_cage < r_pair:
        return r_cage, "진동", lind, "결정"
    return r_pair, "쌍", lind, ("결정" if lind < lindemann else "유체")


def a_mean_star(phi: float) -> float:
    """평균 간격 a_mean/d = √(π/4φ)  (2D, φ = Nπd²/4L²)."""
    return math.sqrt(math.pi / (4 * phi))


__all__ = ["R_WCA", "HEX_NN", "U_star", "U2_star", "approach_distance", "a_mean_star"]
