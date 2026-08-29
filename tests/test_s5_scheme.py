"""S5 — HOOMD 적분기·RNG 벤치마크 회귀. [slow]

`knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md` 의 측정을 고정한다.
HOOMD 버전을 올리거나 forces.py 를 고쳤을 때 조용히 깨지는 것을 잡는다.

허용오차는 **이론 통계오차의 배수**로 쓴다 (관측값에 맞추지 않는다).
seed 는 고정한다 — HOOMD `Brownian` 은 같은 seed 에서 비트 단위 재현된다.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from simbot.estimators import euler_maruyama_trap_variance_bias
from simbot.guards import assert_statistic_fluctuates

from .conftest import TRAP_DIM, TRAP_N, se_of_mean, sigma_away

pytestmark = pytest.mark.slow

N_SIGMA = 4.0        # 4σ — 고정 seed 이므로 flaky 하지 않고, 이론 오차에서 나온 값


def _equilibrium_variance(sim, positions_of, dt_star, n_blocks=25):
    """평형 후 성분별 <x^2> 를 블록별로 수집한다."""
    sim.run(int(20 / dt_star))                       # 평형화 20 tau_trap
    stride = max(1, int(2.0 / dt_star))              # 2 tau_trap 마다 → 시간 독립
    blocks = []
    for _ in range(n_blocks):
        sim.run(stride)
        blocks.append(float(np.mean(positions_of(sim) ** 2)))
    return blocks


# =============================================================================
# B1 — 적분 스킴은 Euler-Maruyama 다
# =============================================================================
@pytest.mark.benchmark
@pytest.mark.parametrize("dt_star,n_blocks,require_rejection", [
    # dt* 가 클수록 편향(dt*/2)이 커져 판별력이 생긴다.
    # require_rejection 은 관측이 아니라 **설계 검정력**에서 결정된다:
    #   기대 판별력 = 편향 / SE.  3σ 를 못 만드는 설계에서 3σ 를 요구하면
    #   달성 불가능한 assert 가 된다 (2026-07-28 에 실제로 그렇게 썼다).
    (2.0e-2, 25, True),    # 편향 1.01 % → 기대 ~3.8σ.  기각 요구
    (1.0e-2, 25, False),   # 편향 0.50 % → 기대 ~1.8σ.  INCONCLUSIVE 를 사실로 고정
])
def test_B1_brownian_is_euler_maruyama(trap_sim_factory, positions_of, dt_star,
                                      n_blocks, require_rejection):
    """★ <x*^2> = 1/(1-dt*/2). 판별력이 있는 곳에서만 경쟁 가설을 기각한다.

    ★ **적분기 검증은 일부러 큰 dt* 에서 해야 한다.** dt* 를 줄이면 편향이 줄지만
      검증 가능성도 함께 줄어든다 — 프로덕션 dt*(5e-3) 에서는 적분기를 검증할 수 없다.
      근거: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md §1
    """
    blocks = []
    for seed in (1, 2, 3, 4):
        sim, _ = trap_sim_factory(dt_star, seed)
        blocks += _equilibrium_variance(sim, positions_of, dt_star, n_blocks=n_blocks)

    assert_statistic_fluctuates(blocks, "block <x^2>")

    measured = float(np.mean(blocks))
    se = se_of_mean(blocks)
    bias = euler_maruyama_trap_variance_bias(dt_star)
    em_pred = 1.0 + bias

    dev_em = sigma_away(measured, em_pred, se)
    dev_exact = sigma_away(measured, 1.0, se)
    power = bias / se                      # 설계 검정력: EM 과 exact 를 가를 수 있는 σ

    assert dev_em < N_SIGMA, (
        f"dt*={dt_star:g}: 측정 {measured:.5f}±{se:.5f} 가 EM 예측 {em_pred:.5f} 에서 "
        f"{dev_em:.1f}σ 벗어났다")

    if require_rejection:
        assert power > 3.0, (
            f"dt*={dt_star:g}: 설계 검정력이 {power:.1f}σ 뿐이다 — n_blocks 를 늘려야 "
            f"이 케이스가 의미를 갖는다 (필요 배수 {(3.0/power)**2:.1f}x)")
        assert dev_exact > 3.0, (
            f"dt*={dt_star:g}: exact 스킴 가설(1.0)을 {dev_exact:.1f}σ 로만 기각한다")
    else:
        # 판별 불가 구간임을 **사실로 고정한다.** 여기서 3σ 기각이 가능해졌다면
        # 통계가 개선된 것이므로 require_rejection=True 로 승격해야 한다.
        assert power < 3.0, (
            f"dt*={dt_star:g}: 검정력이 {power:.1f}σ 로 올라갔다 — "
            f"require_rejection=True 로 승격하라")


@pytest.mark.benchmark
def test_B5_em_bias_halves_when_dt_halves(trap_sim_factory, positions_of):
    """★ 1차 스킴 확인: dt* 를 절반으로 하면 편향도 절반 (2차면 1/4)."""
    bias = {}
    for dt_star in (2.0e-2, 1.0e-2):
        blocks = []
        for seed in (11, 12, 13, 14):
            sim, _ = trap_sim_factory(dt_star, seed)
            blocks += _equilibrium_variance(sim, positions_of, dt_star)
        bias[dt_star] = (float(np.mean(blocks)) - 1.0, se_of_mean(blocks))

    b_hi, se_hi = bias[2.0e-2]
    b_lo, se_lo = bias[1.0e-2]
    ratio = b_hi / b_lo
    # 비의 오차 전파
    se_ratio = abs(ratio) * math.hypot(se_hi / b_hi, se_lo / b_lo)

    assert sigma_away(ratio, 2.0, se_ratio) < N_SIGMA, (
        f"편향 비 {ratio:.2f}±{se_ratio:.2f} 가 1차 예측 2.0 에서 벗어났다")
    assert abs(ratio - 4.0) > 2 * se_ratio, (
        f"2차 스킴 가설(비=4)을 기각하지 못한다 (비 {ratio:.2f}±{se_ratio:.2f})")


# =============================================================================
# B7 — N개 비상호작용 입자는 독립표본이다
# =============================================================================
@pytest.mark.benchmark
def test_B7_particles_are_independent_samples(trap_sim_factory, positions_of):
    """★ Var_t(mean_i d_i) * N / <Var_i(d_i)> = 1 이어야 한다.

    "1% 정밀도 1.3초" 라는 비용 주장 전체가 이 등식에 서 있다.
    깨지면 오차막대가 실제보다 작아져 **거짓 정밀도를 보고**하게 된다.

    ⚠ 이 검사에서 변위의 평균을 빼면 안 된다 — sum(d)=0 이 항등적으로 성립해
      교차상관이 -1/(n-1) 로 고정되고 측정이 무의미해진다 (findings §3).
    """
    dt_star, M = 5e-3, 2000
    sim, _ = trap_sim_factory(dt_star, seed=101, with_trap=False)  # 순수 노이즈
    sim.run(10)

    means, variances = [], []
    prev = positions_of(sim)
    for _ in range(M):
        sim.run(1)
        cur = positions_of(sim)
        d = (cur - prev)[:, 0]                 # 평균을 빼지 않는다
        means.append(float(d.mean()))
        variances.append(float(d.var(ddof=1)))
        prev = cur

    assert_statistic_fluctuates(means, "mean_i d_i")
    assert_statistic_fluctuates(variances, "Var_i d_i")

    ratio = float(np.var(means, ddof=1)) * TRAP_N / float(np.mean(variances))
    se_ratio = ratio * math.sqrt(2.0 / (M - 1))

    assert sigma_away(ratio, 1.0, se_ratio) < N_SIGMA, (
        f"독립성 비 {ratio:.4f}±{se_ratio:.4f} != 1 → 실효표본 {TRAP_N/ratio:.0f}/{TRAP_N}")


# =============================================================================
# B9 — 노이즈는 균일분포다 (Gaussian 이 아니다)
# =============================================================================
@pytest.mark.benchmark
def test_B9_noise_is_uniform_not_gaussian(trap_sim_factory, positions_of):
    """★ 첨도 1.80 (균일) vs 3.00 (Gaussian), max/sigma = sqrt(3).

    2차 모멘트는 영향 없지만 **꼬리가 잘려 있어 장벽 넘기 속도가 틀린다.**
    이 사실이 잊히면 탈출률 계산에서 조용히 틀린다.
    """
    dt_star, n = 1e-3, 20000
    sim, _ = trap_sim_factory(dt_star, seed=99, n=n, with_trap=False)
    q0 = positions_of(sim)[:, 0]
    sim.run(1)
    step = positions_of(sim)[:, 0] - q0

    sd = float(step.std())
    # 분산은 요동-소산으로 맞아야 한다: sqrt(2 D* dt*), D*=1
    assert sd == pytest.approx(math.sqrt(2 * dt_star), rel=5e-3)

    kurt = float(np.mean((step / sd) ** 4))
    assert kurt == pytest.approx(1.80, abs=0.05), f"첨도 {kurt:.4f} — 균일분포는 1.80"
    assert abs(kurt - 3.0) > 1.0, "Gaussian 가설(3.00)이 기각되어야 한다"

    max_over_sd = float(np.abs(step).max() / sd)
    assert max_over_sd == pytest.approx(math.sqrt(3), rel=0.02), (
        f"max/sigma = {max_over_sd:.4f} — 균일분포의 구조적 상한은 sqrt(3)=1.7321")


# =============================================================================
# 재현성 — 통계 테스트의 전제
# =============================================================================
def test_same_seed_reproduces_bitwise(trap_sim_factory, positions_of):
    """★ 고정 seed 통계 테스트가 정당한 이유. 깨지면 위 테스트들이 flaky 해진다."""
    outs = []
    for _ in range(2):
        sim, _ = trap_sim_factory(5e-3, seed=42, n=200)
        sim.run(500)
        outs.append(positions_of(sim))
    assert np.array_equal(outs[0], outs[1]), "같은 seed 인데 결과가 다르다"


def test_different_seeds_diverge(trap_sim_factory, positions_of):
    """seed 가 실제로 효과가 있는가 — 없으면 '독립 시드 4개'가 거짓이 된다."""
    outs = []
    for seed in (42, 43):
        sim, _ = trap_sim_factory(5e-3, seed=seed, n=200)
        sim.run(500)
        outs.append(positions_of(sim))
    assert not np.allclose(outs[0], outs[1])


# =============================================================================
# 물리 정합성 — 차원 스케일링 (B2/B3)
# =============================================================================
@pytest.mark.benchmark
def test_B2_radial_variance_scales_with_dimension(trap_sim_factory, positions_of):
    """<r^2>(3D)/<r^2>(2D) = 1.5 — 01_intake.md 모호성 A1 의 판별자."""
    dt_star = 1.0e-2
    res = {}
    for dim in (2, 3):
        vals = []
        for seed in (21, 22):
            sim, _ = trap_sim_factory(dt_star, seed, dim=dim)
            sim.run(int(20 / dt_star))
            stride = int(2.0 / dt_star)
            for _ in range(25):
                sim.run(stride)
                p = positions_of(sim, dim=dim)
                vals.append(float(np.mean(np.sum(p**2, axis=1))))   # <r^2>
        res[dim] = (float(np.mean(vals)), se_of_mean(vals))

    (r2_2d, se2), (r2_3d, se3) = res[2], res[3]
    ratio = r2_3d / r2_2d
    se_ratio = ratio * math.hypot(se3 / r2_3d, se2 / r2_2d)
    assert sigma_away(ratio, 1.5, se_ratio) < N_SIGMA, (
        f"<r^2> 비 {ratio:.4f}±{se_ratio:.4f} != 1.5")
