"""S0 — 단위·상수·무차원화 척도. 모든 단계의 기반.

이 파일이 잡으려는 버그: **조용히 틀리는 종류.**
반지름/직경 혼동, kT 와 eta 혼동, 왕복 변환 손실 — 발산하지 않으므로
테스트가 없으면 영영 모른다.
"""
from __future__ import annotations

import math

import pytest

from simbot.units import (
    K_B, Scales, kT_si, scales_brownian, scales_harmonic_trap,
    stokes_drag_si, stokes_einstein_D_si, water_density_si, water_viscosity_si,
)


# --- 상수 -------------------------------------------------------------------
def test_boltzmann_is_si_2019_exact_value():
    assert K_B == 1.380649e-23


@pytest.mark.parametrize("T_si,expected_pN_nm", [
    (293.15, 4.047373),   # 20 C  — 정정 2026-07-28: 4.047872 로 잘못 적었던 값
    (298.15, 4.116405),   # 25 C
    (300.00, 4.141947),   # 300 K — 첫 손그림
])
def test_kT_in_pN_nm(T_si, expected_pN_nm):
    """1 pN*nm = 1e-21 J. 광집게 문헌이 이 단위를 쓴다."""
    assert kT_si(T_si) * 1e21 == pytest.approx(expected_pN_nm, rel=1e-6)


# --- 물 물성 ----------------------------------------------------------------
def test_water_viscosity_at_table_points_is_exact():
    """표에 있는 점은 보간하지 않는다."""
    eta, extrap = water_viscosity_si(298.15)
    assert eta == 0.8900e-3
    assert extrap is False


def test_water_viscosity_20C_and_25C_differ_by_11_percent():
    """★ 가장 흔한 혼동: 1.002 mPa*s 는 20 C 값이고 25 C 는 0.890 mPa*s.

    이 차이(11 %)가 D0 와 모든 시간척도에 그대로 전파된다.
    """
    eta20, _ = water_viscosity_si(293.15)
    eta25, _ = water_viscosity_si(298.15)
    assert eta20 == pytest.approx(1.0016e-3)
    assert eta25 == pytest.approx(0.8900e-3)
    assert (eta20 - eta25) / eta25 == pytest.approx(0.1254, rel=0.01)


def test_water_viscosity_interpolates_monotonically():
    ts = [293.15, 295.0, 298.15, 300.0, 303.15, 306.0, 308.15]
    etas = [water_viscosity_si(t)[0] for t in ts]
    assert all(a > b for a, b in zip(etas, etas[1:])), "점도는 온도에 단조감소"


def test_water_viscosity_at_300K_matches_recorded_value():
    """첫 손그림의 T = 300 K. 예측 문서에 봉인된 값과 일치해야 한다."""
    eta, extrap = water_viscosity_si(300.0)
    assert extrap is False
    assert eta == pytest.approx(8.5566e-4, rel=1e-4)


def test_water_viscosity_flags_extrapolation():
    """보간 범위(293-308 K) 밖은 외삽임을 알려야 한다 — provenance 를 낮추기 위해."""
    _, extrap_lo = water_viscosity_si(280.0)
    _, extrap_hi = water_viscosity_si(330.0)
    _, extrap_in = water_viscosity_si(300.0)
    assert extrap_lo is True and extrap_hi is True and extrap_in is False


def test_water_density_is_near_1000():
    rho, _ = water_density_si(298.15)
    assert 990 < rho < 1000


# --- Stokes 항력: 반지름 vs 직경 -------------------------------------------
def test_stokes_drag_uses_radius_not_diameter():
    """★ `master_plan` §S3 이 "가장 흔한 실수"로 지목한 버그를 고정한다.

    gamma = 6*pi*eta*a  (a = 반지름).  직경을 넣으면 gamma 가 2배가 되고
    D0 가 절반이 되어 **모든 시간척도가 2배 틀린다.** 발산하지 않는다.
    """
    eta, a = 1e-3, 5e-6
    gamma = stokes_drag_si(eta, a)
    assert gamma == pytest.approx(6 * math.pi * eta * a)
    # 직경을 넣으면 정확히 2배 — 이 관계를 테스트로 못박아 둔다
    assert stokes_drag_si(eta, 2 * a) == pytest.approx(2 * gamma)


def test_stokes_einstein_reference_case_1um_sphere_in_water_25C():
    """문헌 관례값: 물 25 C 에서 지름 1 um 구는 D ~ 0.49 um^2/s.

    knowledge/wiki/concepts/water-298k.md 의 참조 케이스.
    """
    T, a = 298.15, 0.5e-6
    eta, _ = water_viscosity_si(T)
    gamma = stokes_drag_si(eta, a)
    D0 = stokes_einstein_D_si(T, gamma)
    assert D0 * 1e12 == pytest.approx(0.4909, rel=2e-3)   # um^2/s
    tau_B = (2 * a) ** 2 / D0
    assert tau_B == pytest.approx(2.037, rel=2e-3)        # s


# --- Scales 왕복 ------------------------------------------------------------
ALL_KINDS = ["length", "energy", "time", "force", "stiffness", "velocity",
             "diffusivity", "rate", "modulus_3d", "area", "volume"]


@pytest.mark.parametrize("kind", ALL_KINDS)
def test_scales_roundtrip_is_lossless(kind):
    """SI -> star -> SI 왕복 오차가 부동소수점 한계 이내여야 한다.

    `master_plan` §S4 게이트: 상대오차 < 1e-12
    """
    s = Scales(length_si=2.0352e-8, energy_si=4.1419e-21, time_si=8.0644e-3)
    value_si = 3.14159e-7
    back = s.to_si(s.to_star(value_si, kind), kind)
    assert back == pytest.approx(value_si, rel=1e-14)


def test_scales_rejects_unknown_kind():
    s = Scales(length_si=1.0, energy_si=1.0, time_si=1.0)
    with pytest.raises(KeyError, match="unknown scale kind"):
        s.to_star(1.0, "wibble")


def test_derived_scales_are_dimensionally_consistent():
    s = Scales(length_si=3.0, energy_si=7.0, time_si=11.0)
    assert s.force_si == pytest.approx(7.0 / 3.0)
    assert s.stiffness_si == pytest.approx(7.0 / 9.0)
    assert s.velocity_si == pytest.approx(3.0 / 11.0)
    assert s.diffusivity_si == pytest.approx(9.0 / 11.0)
    assert s.rate_si == pytest.approx(1.0 / 11.0)
    assert s.modulus_3d_si == pytest.approx(7.0 / 27.0)


# --- 카드별 척도의 정의 불변식 ---------------------------------------------
def test_harmonic_trap_scales_normalize_D_and_k_to_exactly_one():
    """★ 조화 트랩 카드 §3 의 핵심 주장.

    (l_trap, kT, tau_trap) 을 고르면 무차원 운동방정식이 파라미터 없이
    `dr*/dt* = -r* + sqrt(2) xi` 로 정규화된다  =>  D* = 1, k* = 1, <x*^2> = 1.
    """
    T, eta, a, k = 300.0, 8.5566e-4, 5e-6, 1e-5
    gamma = stokes_drag_si(eta, a)
    s = scales_harmonic_trap(k_si=k, T_si=T, gamma_si=gamma)

    D0 = stokes_einstein_D_si(T, gamma)
    assert s.to_star(D0, "diffusivity") == pytest.approx(1.0, rel=1e-14)
    assert s.to_star(k, "stiffness") == pytest.approx(1.0, rel=1e-14)
    # 등분배 <x^2> = kT/k 도 무차원으로 정확히 1
    assert s.to_star(kT_si(T) / k, "area") == pytest.approx(1.0, rel=1e-14)


def test_brownian_scales_normalize_D_and_tauD_to_exactly_one():
    """수동 구형 × 수송 카드: (sigma, kT, tau_D) 를 고르면 D* = 1, tau_D* = 1."""
    T, eta, a = 298.15, 0.8900e-3, 0.5e-6
    gamma = stokes_drag_si(eta, a)
    sigma = 2 * a
    s = scales_brownian(sigma_si=sigma, T_si=T, gamma_si=gamma)

    D0 = stokes_einstein_D_si(T, gamma)
    assert s.to_star(D0, "diffusivity") == pytest.approx(1.0, rel=1e-14)
    assert s.to_star(sigma**2 / D0, "time") == pytest.approx(1.0, rel=1e-14)


def test_two_cards_give_time_scales_separated_by_k_star_sigma():
    """★ 카드 체계가 필요한 이유를 수치로 고정한다.

    tau_D / tau_trap = k sigma^2 / kT = k*_sigma
    첫 손그림에서 이 값이 2.41e5 다  =>  tau_D 기준 dt 는 완화시간의 12배가 된다.
    """
    T, eta, a, k = 300.0, 8.5566e-4, 5e-6, 1e-5
    gamma = stokes_drag_si(eta, a)
    sigma = 2 * a
    s_bd = scales_brownian(sigma_si=sigma, T_si=T, gamma_si=gamma)
    s_tr = scales_harmonic_trap(k_si=k, T_si=T, gamma_si=gamma)

    ratio = s_bd.time_si / s_tr.time_si
    k_star_sigma = k * sigma**2 / kT_si(T)
    assert ratio == pytest.approx(k_star_sigma, rel=1e-12)
    assert ratio == pytest.approx(2.414e5, rel=1e-3)

    # tau_D 기준 dt* = 5e-5 를 쓰면 완화시간의 몇 배가 되는가
    dt_in_tau_trap = 5e-5 * ratio
    assert dt_in_tau_trap == pytest.approx(12.07, rel=1e-2)
    assert dt_in_tau_trap > 1.0, "보편 규약이 이 계에서 완화시간을 넘는다는 사실을 고정"
