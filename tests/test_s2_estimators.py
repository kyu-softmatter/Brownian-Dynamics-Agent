"""S2 — 예측 엔진. 해석해가 정말 해석해인지 검사한다.

여기 있는 값들이 `runs/*/02_prediction.md` 에 봉인되므로, **예측 엔진이 틀리면
검증 전체가 무의미해진다.** 그래서 항등식과 극한으로 이중 확인한다.
"""
from __future__ import annotations

import math

import pytest

from simbot.estimators import (
    EFFICIENCY_BY_K, THROUGHPUT_PARTICLE_STEPS_PER_S, dt_star_for_trap_bias,
    estimate_wall_time_s, euler_maruyama_trap_variance_bias, harmonic_trap,
    overdamped_validity, samples_for_variance_precision, trap_run_length,
)
from simbot.units import kT_si, stokes_drag_si, stokes_einstein_D_si, water_viscosity_si

# 첫 손그림의 계 — 이 값들이 02_prediction.md 에 봉인되어 있다
SKETCH = dict(T_si=300.0, eta_si=8.5566e-4, radius_si=5e-6, k_si=1e-5)


# =============================================================================
# 조화 트랩 해석해 — 항등식으로 검사
# =============================================================================
@pytest.mark.parametrize("dim", [2, 3])
def test_equipartition_per_component_is_kT_over_k_regardless_of_dim(dim):
    """★ <x_i^2> = kT/k 는 차원에 무관하다. 그래서 차원을 판별하지 못한다."""
    p = harmonic_trap(dim=dim, **SKETCH)
    assert p.var_per_component_si == pytest.approx(p.kT_si / p.k_si, rel=1e-14)


@pytest.mark.parametrize("dim", [2, 3])
def test_radial_variance_scales_with_dim(dim):
    """<r^2> = d * kT/k. 이것이 차원 판별자다 (모호성 A1)."""
    p = harmonic_trap(dim=dim, **SKETCH)
    assert p.rms_radial_si**2 == pytest.approx(dim * p.kT_si / p.k_si, rel=1e-14)


def test_dimension_discriminator_ratio_is_exactly_three_halves():
    """★ 01_intake.md 모호성 A1 의 판별자: <r^2>(3D)/<r^2>(2D) = 3/2 정확히."""
    p2 = harmonic_trap(dim=2, **SKETCH)
    p3 = harmonic_trap(dim=3, **SKETCH)
    assert (p3.rms_radial_si**2) / (p2.rms_radial_si**2) == pytest.approx(1.5, rel=1e-14)
    assert p3.msd_plateau_si / p2.msd_plateau_si == pytest.approx(1.5, rel=1e-14)


@pytest.mark.parametrize("dim", [2, 3])
def test_msd_plateau_is_twice_radial_variance(dim):
    """MSD(t->inf) = <|r(inf)-r(0)|^2> = 2<r^2>  (두 위치가 독립이므로)."""
    p = harmonic_trap(dim=dim, **SKETCH)
    assert p.msd_plateau_si == pytest.approx(2 * dim * p.var_per_component_si, rel=1e-14)


def test_msd_reaches_plateau_at_long_times():
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.msd_si(50 * p.tau_trap_si) == pytest.approx(p.msd_plateau_si, rel=1e-9)


def test_msd_short_time_limit_recovers_free_diffusion():
    """★ 비자명한 극한 검사: t << tau_trap 에서 MSD -> 2 d D0 t.

    (2d kT/k)(1 - e^{-t/tau}) ~ (2d kT/k)(t k/gamma) = 2 d (kT/gamma) t = 2 d D0 t
    트랩 해석해와 Stokes-Einstein 이 서로 정합함을 확인한다.
    """
    p = harmonic_trap(dim=2, **SKETCH)
    t = 1e-4 * p.tau_trap_si
    free = 2 * p.dim * p.D0_si * t
    assert p.msd_si(t) == pytest.approx(free, rel=1e-3)


def test_relaxation_time_and_corner_frequency_are_consistent():
    """f_c = 1/(2 pi tau_trap). 광집게 실험이 f_c 를 재므로 중요하다."""
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.tau_trap_si == pytest.approx(p.gamma_si / p.k_si, rel=1e-14)
    assert p.corner_freq_si == pytest.approx(1 / (2 * math.pi * p.tau_trap_si), rel=1e-14)


def test_confinement_length_is_sqrt_kT_over_k():
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.l_trap_si == pytest.approx(math.sqrt(p.kT_si / p.k_si), rel=1e-14)
    # <x^2> = l_trap^2 — 같은 것을 두 이름으로 부르고 있는지 확인
    assert p.l_trap_si**2 == pytest.approx(p.var_per_component_si, rel=1e-14)


def test_time_scale_separation_equals_k_star_sigma():
    """tau_D/tau_trap = k sigma^2/kT. 카드 선택의 근거가 되는 항등식."""
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.tau_sep == pytest.approx(p.k_star_sigma, rel=1e-12)
    assert p.l_trap_over_sigma == pytest.approx(1 / math.sqrt(p.k_star_sigma), rel=1e-12)


def test_estimator_warns_when_tau_D_would_be_wrong_time_unit():
    """강한 트랩에서 카드 경고가 실제로 발동하는가."""
    p = harmonic_trap(dim=2, **SKETCH)
    joined = " ".join(p.notes)
    assert "tau_trap" in joined and "tau_D" in joined
    assert "배제부피" in joined


# =============================================================================
# B6 — SI 참조값. 02_prediction.md 에 봉인된 숫자를 고정한다
# =============================================================================
@pytest.mark.benchmark
def test_B6_sketch_reference_values():
    """조화 트랩 카드 B6.  a=5um, k=10 pN/um, T=300K, 물.

    이 숫자가 바뀌면 runs/2026-07-28_trap-2d-5um_2dfb9d/02_prediction.md 의
    봉인된 예측이 무효가 된다.
    """
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.kT_si * 1e21 == pytest.approx(4.14195, rel=1e-5)     # pN*nm
    assert p.gamma_si == pytest.approx(8.0644e-8, rel=1e-4)       # kg/s
    assert p.D0_si * 1e12 == pytest.approx(0.05136, rel=1e-3)     # um^2/s
    assert p.tau_trap_si * 1e3 == pytest.approx(8.064, rel=1e-3)  # ms
    assert p.tau_D_si == pytest.approx(1947.0, rel=1e-3)          # s
    assert p.l_trap_si * 1e9 == pytest.approx(20.35, rel=1e-3)    # nm
    assert p.corner_freq_si == pytest.approx(19.735, rel=1e-3)    # Hz
    assert p.k_star_sigma == pytest.approx(2.4143e5, rel=1e-3)
    assert p.var_per_component_si * 1e18 == pytest.approx(414.19, rel=1e-4)  # nm^2
    assert p.msd_plateau_si * 1e18 == pytest.approx(1656.78, rel=1e-4)       # nm^2


def test_sketch_values_reproduce_from_water_table():
    """SKETCH 의 eta 가 units 모듈의 물 표에서 온 값과 일치하는가.

    두 곳에 하드코딩된 상수가 어긋나는 것을 막는다.
    """
    eta, extrap = water_viscosity_si(300.0)
    assert extrap is False
    assert eta == pytest.approx(SKETCH["eta_si"], rel=1e-4)


# =============================================================================
# Euler-Maruyama 계통편향
# =============================================================================
@pytest.mark.parametrize("dt_star,expected_pct", [
    (2.0e-2, 1.0101), (1.0e-2, 0.5025), (5.0e-3, 0.2506),
    (2.5e-3, 0.1252), (1.0e-3, 0.0500),
])
def test_em_bias_known_values(dt_star, expected_pct):
    """Var* = 1/(1-dt*/2). 이 값들이 dt 래더 검증의 기준선이다."""
    assert euler_maruyama_trap_variance_bias(dt_star) * 100 == pytest.approx(
        expected_pct, rel=2e-3)


def test_em_bias_is_first_order_in_dt():
    """★ 1차 스킴임을 고정: dt 를 절반으로 하면 편향도 (거의) 절반."""
    b1 = euler_maruyama_trap_variance_bias(1e-2)
    b2 = euler_maruyama_trap_variance_bias(5e-3)
    assert b1 / b2 == pytest.approx(2.0, rel=5e-3)
    # 2차 스킴이라면 비가 4가 되어야 한다 — 그 가설을 기각
    assert abs(b1 / b2 - 4.0) > 1.0


def test_em_bias_leading_order_matches_dt_over_two():
    for dt in (1e-3, 1e-4):
        assert euler_maruyama_trap_variance_bias(dt) == pytest.approx(dt / 2, rel=1e-2)


def test_em_bias_rejects_unstable_dt():
    """dt* >= 2 는 결정론 항 (1-dt*) 이 |.|>1 이 되어 발산한다."""
    for bad in (0.0, -1e-3, 2.0, 3.0):
        with pytest.raises(ValueError):
            euler_maruyama_trap_variance_bias(bad)


def test_dt_star_for_bias_is_inverse_of_bias():
    for target in (1e-4, 1e-3, 2.5e-3, 1e-2):
        dt = dt_star_for_trap_bias(target)
        assert euler_maruyama_trap_variance_bias(dt) == pytest.approx(target, rel=1e-12)


# =============================================================================
# 전제 검사 (과감쇠 · Stokes)
# =============================================================================
def test_sketch_system_passes_overdamped_and_stokes():
    p = harmonic_trap(dim=2, **SKETCH)
    V = 4 / 3 * math.pi * SKETCH["radius_si"] ** 3
    c = overdamped_validity(
        gamma_si=p.gamma_si, mass_si=1050 * V, tau_process_si=p.tau_trap_si,
        velocity_scale_si=p.l_trap_si / p.tau_trap_si,
        radius_si=SKETCH["radius_si"], eta_si=SKETCH["eta_si"], rho_fluid_si=996.5)
    assert c.passed, c.failures
    assert c.inertial_ratio == pytest.approx(8.454e-4, rel=1e-2)
    assert c.reynolds == pytest.approx(1.470e-5, rel=1e-2)


def test_overdamped_check_fails_for_heavy_particle():
    """가드가 실제로 발동하는가 — 통과만 테스트하면 가드가 죽어도 모른다."""
    c = overdamped_validity(
        gamma_si=1e-8, mass_si=1.0, tau_process_si=1e-3,
        velocity_scale_si=1.0, radius_si=1e-6, eta_si=1e-3, rho_fluid_si=1000.0)
    assert not c.passed
    assert any("관성" in f for f in c.failures)
    assert any("Reynolds" in f for f in c.failures)


# =============================================================================
# 통계 정밀도 · 비용 모델
# =============================================================================
def test_samples_for_precision_follows_inverse_square():
    """SE(s^2)/s^2 = sqrt(2/n)  =>  n = 2/rel^2"""
    assert samples_for_variance_precision(0.01) == pytest.approx(20000)
    assert samples_for_variance_precision(0.005) == pytest.approx(80000)
    # 정밀도를 2배로 올리면 표본은 4배
    assert (samples_for_variance_precision(0.005)
            / samples_for_variance_precision(0.01)) == pytest.approx(4.0)


def test_trap_run_length_matches_prediction_document():
    """02_prediction.md §6 의 비용표를 고정한다."""
    p = harmonic_trap(dim=2, **SKETCH)
    r = trap_run_length(n_particles=1000, rel_err_target=0.01,
                        tau_trap_si=p.tau_trap_si)
    assert r["independent_samples_needed"] == pytest.approx(20000)
    assert r["t_total_in_tau_trap"] == pytest.approx(40.0)
    steps = r["t_total_in_tau_trap"] / 5e-3
    assert steps == pytest.approx(8000)
    assert estimate_wall_time_s(1000, int(steps), 1) == pytest.approx(1.27, rel=0.05)


def test_wall_time_scales_linearly_in_work():
    assert estimate_wall_time_s(1000, 1000) == pytest.approx(
        estimate_wall_time_s(500, 2000), rel=1e-12)
    assert estimate_wall_time_s(1000, 2000) == pytest.approx(
        2 * estimate_wall_time_s(1000, 1000), rel=1e-12)


def test_wall_time_uses_measured_efficiency_and_penalizes_concurrency():
    """실측 효율표를 쓰는가. k=8 은 k=1 보다 프로세스당 느려야 한다."""
    w1 = estimate_wall_time_s(1000, 10000, 1)
    w8 = estimate_wall_time_s(1000, 10000, 8)
    assert w8 > w1
    assert w8 / w1 == pytest.approx(1 / EFFICIENCY_BY_K[8], rel=1e-12)


def test_wall_time_falls_back_for_unmeasured_concurrency():
    """표에 없는 k 는 아래쪽 실측값으로 보수적으로 근사한다."""
    assert estimate_wall_time_s(100, 100, 7) == pytest.approx(
        estimate_wall_time_s(100, 100, 6), rel=1e-12)


def test_throughput_constant_matches_measurement():
    """knowledge/wiki/findings/local-cpu-parallelism.md 의 실측 상수."""
    assert THROUGHPUT_PARTICLE_STEPS_PER_S == pytest.approx(6.3e6, rel=1e-9)
    assert EFFICIENCY_BY_K[4] == pytest.approx(0.925, abs=1e-3)
    assert EFFICIENCY_BY_K[12] < EFFICIENCY_BY_K[10], "k=12 는 회귀 — 측정 사실"
