"""S3 — 절단 반경 제안. HOOMD 불필요 (순수 수치), 전부 빠름.

핵심 검증: **`r_cut` 이 관습이 아니라 명시된 오차 허용치에서 나오는가**,
그리고 **희박·장거리에서는 늘리고 조밀에서는 예산으로 제한하는가.**
"""
from __future__ import annotations

import math

import pytest

from simbot.cutoff import (
    CONVENTIONAL_LJ_RCUT, DEFAULT_NEIGHBOR_BUDGET, LAMBDA_BASELINE_PHI,
    LAMBDA_BASELINE_RCUT, TOLERANCE_PRESETS, WCA_RCUT_OVER_SIGMA_LJ,
    barker_henderson_diameter, cutoff_from_tolerance, lj_cutoff, morse_cutoff,
    neighbor_list_buffer, neighbors_per_particle, pair_cost_vs_lambda_baseline,
    propose_cutoff, wca_cutoff, yukawa_cutoff,
)


# =============================================================================
# WCA — 정의상 고정
# =============================================================================
def test_wca_cutoff_is_the_lj_minimum_exactly():
    p = wca_cutoff()
    assert p.r_cut_over_sigma == pytest.approx(2 ** (1 / 6), rel=1e-15)
    assert p.r_cut_over_sigma == pytest.approx(1.1224620483, rel=1e-9)
    assert p.exact is True
    assert p.criterion == "potential_definition"
    # 최솟점이므로 힘과 포텐셜이 정확히 0 — 절단 인공물이 없다
    assert p.beta_U_at_cut == 0.0 and p.beta_F_sigma_at_cut == 0.0


def test_wca_force_vanishes_at_its_cutoff():
    """★ r_cut = 2^{1/6} sigma 에서 LJ 힘이 정확히 0 임을 직접 확인."""
    r = WCA_RCUT_OVER_SIGMA_LJ
    f = 24.0 * (2.0 * r ** -13 - r ** -7)          # -dU/dr, eps=sigma=1
    assert f == pytest.approx(0.0, abs=1e-14)


def test_wca_hoomd_recipe_records_the_shift_requirement():
    """랩 코드가 빠뜨린 mode="shift" 를 레시피가 명시하는가."""
    p = wca_cutoff()
    assert 'mode="shift"' in p.hoomd_note
    assert "2**(1/6)" in p.hoomd_note or "2^{1/6}" in p.rationale


def test_propose_cutoff_short_circuits_for_exact_potentials():
    """WCA 는 phi 나 예산과 무관하게 값이 바뀌지 않는다."""
    for phi, n in [(0.0, 1), (0.4, 10000)]:
        p = propose_cutoff(wca_cutoff(), phi=phi, n_particles=n)
        assert p.r_cut_over_sigma == pytest.approx(WCA_RCUT_OVER_SIGMA_LJ, rel=1e-15)
        assert p.criterion == "potential_definition"


# =============================================================================
# Barker-Henderson — 지식 베이스와의 독립 대조
# =============================================================================
@pytest.mark.benchmark
def test_barker_henderson_reproduces_knowledge_base_value():
    """★ knowledge/wiki/findings/wca-reproduces-carnahan-starling.md 는
    `beta*eps = 1` 에서 `d_eff = 1.017 sigma` 라고 적고 있다.

    이 테스트가 그것을 독립적으로 재현한다 (수치 적분).
    """
    d = barker_henderson_diameter(beta_epsilon=1.0)
    assert d == pytest.approx(1.017, rel=3e-3), f"d_BH = {d:.5f}, 지식값 1.017"


def test_barker_henderson_increases_with_stiffness_toward_r_min():
    """eps 가 커지면 유효 경구가 딱딱해져 r_min 에 접근한다 (넘지는 않는다)."""
    ds = [barker_henderson_diameter(be) for be in (0.5, 1.0, 2.0, 10.0, 500.0)]
    assert all(a < b for a, b in zip(ds, ds[1:])), "beta*eps 에 단조증가"
    assert ds[-1] < WCA_RCUT_OVER_SIGMA_LJ, "r_min 을 넘을 수 없다"
    assert ds[-1] == pytest.approx(WCA_RCUT_OVER_SIGMA_LJ, rel=1e-2), "강한 극한에서 접근"


def test_barker_henderson_rejects_nonpositive_epsilon():
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError):
            barker_henderson_diameter(bad)


# =============================================================================
# 허용치 기반 해 — 프리셋과 지배 제약
# =============================================================================
def test_tolerance_presets_are_ordered_from_loose_to_strict():
    order = ["convention", "structure", "thermodynamics", "precision"]
    us = [TOLERANCE_PRESETS[k][0] for k in order]
    assert all(a > b for a, b in zip(us, us[1:])), "허용치가 단조 감소해야 한다"


def test_convention_preset_reproduces_the_literature_2p5_sigma():
    """★ `convention` 프리셋이 문헌 관습값 2.5 sigma 를 재현하는가.

    재현하지 못하면 프리셋의 허용치가 관습과 무관한 값이라는 뜻이다.
    """
    p = lj_cutoff(1.0, preset="convention")
    assert p.r_cut_over_sigma == pytest.approx(CONVENTIONAL_LJ_RCUT, rel=0.05)


def test_stricter_preset_gives_larger_cutoff():
    rs = [lj_cutoff(1.0, preset=k).r_cut_over_sigma
          for k in ("convention", "structure", "thermodynamics", "precision")]
    assert all(a < b for a, b in zip(rs, rs[1:]))


def test_lj_cutoff_reports_residual_at_the_convention():
    """관습값을 쓰면 무엇이 남는지 근거에 적어야 한다 — 판단의 재료."""
    p = lj_cutoff(1.0, preset="thermodynamics")
    assert "2.5 sigma" in p.rationale
    assert p.r_cut_over_sigma == pytest.approx(3.984, rel=1e-3)
    assert p.beta_U_at_cut == pytest.approx(1e-3, rel=1e-6)


def test_cutoff_records_which_constraint_dominated():
    p = lj_cutoff(1.0)
    assert p.criterion in ("potential_tolerance", "force_tolerance")
    assert "지배" in p.rationale


def test_force_tolerance_can_dominate_when_made_strict():
    """★ 두 제약 중 힘이 지배하는 경우도 실제로 발생하는가.

    한쪽 분기만 테스트하면 다른 분기가 죽어도 모른다.
    """
    p = lj_cutoff(1.0, beta_u_tol=1.0, beta_f_sigma_tol=1e-6)
    assert p.criterion == "force_tolerance"
    assert p.beta_F_sigma_at_cut == pytest.approx(1e-6, rel=1e-5)


def test_unknown_preset_raises():
    with pytest.raises(KeyError, match="unknown preset"):
        lj_cutoff(1.0, preset="vibes")


def test_long_range_that_never_meets_tolerance_is_flagged():
    """r_max 까지 가도 허용치에 못 들면 경고해야 한다 (조용히 잘라내면 안 된다)."""
    p = cutoff_from_tolerance(
        potential="coulomb-like", beta_u=lambda r: 100.0 / r,
        beta_f_sigma=lambda r: 100.0 / r ** 2, r_start=1.0, r_max=10.0)
    assert p.criterion == "r_max_reached"
    assert p.warnings and "장거리" in p.warnings[0]


# =============================================================================
# Yukawa — 스크리닝이 사거리를 정한다
# =============================================================================
def test_yukawa_cutoff_grows_as_screening_weakens():
    """★ 계면 콜로이드 카드의 함정: kappa*sigma -> 0 이면 사거리가 폭발한다."""
    rs = [yukawa_cutoff(10.0, ks).r_cut_over_sigma for ks in (10, 5, 2, 1, 0.3, 0.1)]
    assert all(a < b for a, b in zip(rs, rs[1:])), "kappa 감소 → r_cut 증가"
    assert rs[-1] > 50, f"kappa*sigma=0.1 에서 r_cut = {rs[-1]:.1f} sigma"


def test_yukawa_cutoff_scales_like_inverse_kappa():
    """r_cut ~ ln(beta eps / tol)/kappa 이므로 kappa 를 절반으로 하면 대략 2배."""
    r1 = yukawa_cutoff(10.0, 1.0).r_cut_over_sigma
    r2 = yukawa_cutoff(10.0, 0.5).r_cut_over_sigma
    assert 1.7 < r2 / r1 < 2.4


def test_minimum_image_violation_is_detected():
    """★ 계면 콜로이드 실사례 재현: 박스가 상호작용 사거리를 담지 못한다."""
    p = yukawa_cutoff(beta_epsilon=1e5, kappa_sigma=0.1)
    ok, msg = p.box_ok(130.0)
    assert not ok
    assert "minimum image 위반" in msg
    assert f"{p.min_box_over_sigma:.4g}" in msg
    # 충분히 큰 박스면 통과
    assert p.box_ok(p.min_box_over_sigma)[0] is True


# =============================================================================
# Morse
# =============================================================================
def test_morse_cutoff_grows_as_range_alpha_decreases():
    rs = [morse_cutoff(5.0, a).r_cut_over_sigma for a in (20, 10, 5, 3)]
    assert all(a < b for a, b in zip(rs, rs[1:]))


def test_morse_cutoff_exceeds_the_well_position():
    p = morse_cutoff(beta_D0=5.0, alpha_sigma=10.0, r0_over_sigma=1.0)
    assert p.r_cut_over_sigma > 1.0


# =============================================================================
# 이웃 수 — 비용의 전부
# =============================================================================
@pytest.mark.parametrize("dim,coef", [(3, 8.0), (2, 4.0)])
def test_neighbor_count_formula(dim, coef):
    phi, r = 0.3, 2.5
    expected = coef * phi * r ** dim
    assert neighbors_per_particle(r, phi, dim) == pytest.approx(expected, rel=1e-14)


def test_neighbor_count_is_zero_for_single_particle_limit():
    assert neighbors_per_particle(10.0, 0.0, 3) == 0.0


def test_conventional_lj_neighbor_count_is_in_the_expected_range():
    """phi=0.3, r_cut=2.5 sigma 는 3D LJ 액체의 전형값 ~40개여야 한다."""
    nb = neighbors_per_particle(CONVENTIONAL_LJ_RCUT, 0.3, 3)
    assert 30 < nb < 50, f"{nb:.1f}"


def test_lambda_baseline_cost_is_unity_at_measurement_conditions():
    """★ 처리량 상수가 측정된 조건에서 비용 배수가 1 이어야 한다."""
    c = pair_cost_vs_lambda_baseline(LAMBDA_BASELINE_RCUT, LAMBDA_BASELINE_PHI, 3)
    assert c == pytest.approx(1.0, rel=1e-12)


def test_lambda_baseline_was_measured_in_a_cheap_regime():
    """★ Lambda 는 WCA(이웃 3.4개)에서 측정됐다 — 관습 LJ 는 그보다 10배 비싸다.

    이 사실을 잊으면 wall time 예측이 낙관적으로 틀린다.
    """
    nb_base = neighbors_per_particle(LAMBDA_BASELINE_RCUT, LAMBDA_BASELINE_PHI, 3)
    assert nb_base == pytest.approx(3.4, rel=0.05)
    cost_conv = pair_cost_vs_lambda_baseline(CONVENTIONAL_LJ_RCUT, 0.3, 3)
    assert cost_conv > 9.0, f"관습 LJ 비용 {cost_conv:.1f}x"


# =============================================================================
# ★ 최종 제안 — 사용자 규칙 (관습 우선, 희박·장거리는 늘림)
# =============================================================================
def test_single_particle_makes_cutoff_irrelevant():
    p = propose_cutoff(lj_cutoff(1.0), phi=0.0, n_particles=1, dim=2)
    assert p.criterion == "single_particle_free"
    assert any("1개" in w for w in p.warnings)


def test_dilute_system_gets_the_full_tolerance_cutoff():
    """★ 사용자 규칙: 희박하면 관습보다 길어도 된다 — 이웃이 없으니 공짜."""
    tol = lj_cutoff(1.0, preset="thermodynamics")
    p = propose_cutoff(tol, phi=0.01, n_particles=100, dim=3)
    assert p.r_cut_over_sigma == pytest.approx(tol.r_cut_over_sigma, rel=1e-12)
    assert p.criterion == "dilute_or_longrange_affordable"
    assert neighbors_per_particle(p.r_cut_over_sigma, 0.01, 3) < 10


def test_dense_system_is_limited_by_neighbor_budget():
    """조밀계에서는 예산이 걸리고, **남은 오차를 한계로 기록**해야 한다."""
    tol = lj_cutoff(1.0, preset="thermodynamics")
    p = propose_cutoff(tol, phi=0.40, n_particles=4000, dim=3)
    assert p.criterion == "neighbor_budget_limited"
    assert p.r_cut_over_sigma < tol.r_cut_over_sigma
    nb = neighbors_per_particle(p.r_cut_over_sigma, 0.40, 3)
    assert nb == pytest.approx(DEFAULT_NEIGHBOR_BUDGET, rel=1e-9)
    assert any("알려진 한계" in w for w in p.warnings)
    assert math.isnan(p.beta_U_at_cut), "줄였으면 원래 허용치를 주장할 수 없다"


def test_budget_limited_never_goes_below_convention():
    """예산이 아무리 빡빡해도 관습값 아래로는 내리지 않는다."""
    tol = lj_cutoff(1.0, preset="precision")
    p = propose_cutoff(tol, phi=0.6, n_particles=10000, dim=3, neighbor_budget=1.0)
    assert p.r_cut_over_sigma >= CONVENTIONAL_LJ_RCUT


def test_tolerance_below_convention_is_taken_as_free_win():
    """강한 스크리닝이면 허용치 해가 관습보다 작다 — 그대로 쓴다."""
    tol = yukawa_cutoff(beta_epsilon=10.0, kappa_sigma=5.0)
    assert tol.r_cut_over_sigma < CONVENTIONAL_LJ_RCUT
    p = propose_cutoff(tol, phi=0.3, n_particles=2000, dim=3)
    assert p.criterion == "tolerance_below_convention"
    assert p.r_cut_over_sigma == pytest.approx(tol.r_cut_over_sigma, rel=1e-12)


def test_long_range_dilute_is_affordable_but_may_break_the_box():
    """★ 희박 + 장거리: 이웃 예산은 통과하지만 minimum image 가 걸릴 수 있다."""
    tol = yukawa_cutoff(beta_epsilon=10.0, kappa_sigma=0.3)
    p = propose_cutoff(tol, phi=0.001, n_particles=2000, dim=3, box_over_sigma=40.0)
    assert p.criterion == "dilute_or_longrange_affordable"
    assert any("minimum image" in w for w in p.warnings)


def test_proposal_always_reports_cost_and_neighbor_count():
    """근거에 비용과 이웃 수가 항상 들어가야 한다 — 예산 판단의 재료."""
    p = propose_cutoff(lj_cutoff(1.0), phi=0.2, n_particles=1000, dim=3)
    assert "이웃 수" in p.rationale
    assert "비용" in p.rationale


# =============================================================================
# 이웃 리스트 버퍼
# =============================================================================
def test_neighbor_buffer_grows_with_step_displacement():
    b1 = neighbor_list_buffer(1.1225, 0.005)
    b2 = neighbor_list_buffer(1.1225, 0.03)
    assert b2 > b1


def test_neighbor_buffer_is_clamped_to_practical_range():
    tiny = neighbor_list_buffer(2.5, 1e-6)
    huge = neighbor_list_buffer(2.5, 10.0)
    assert tiny == pytest.approx(0.1 * 2.5)
    assert huge == pytest.approx(0.5 * 2.5)


def test_neighbor_buffer_rejects_nonpositive_displacement():
    with pytest.raises(ValueError):
        neighbor_list_buffer(2.5, 0.0)
