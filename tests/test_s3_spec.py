"""S3 — 시스템 명세 데이터 모델 + 타당성 검사.

핵심 회귀: **손으로 쓴 첫 런의 파생값을 코드가 재현하는가.**
`runs/2026-07-28_trap-2d-5um_2dfb9d/03_spec.yaml` 의 `derived:` 블록은 사람이
`simbot` 출력을 옮겨 적은 것이다. 그 값이 정답지이고, 여기서 어긋나면 둘 중
하나가 틀렸다는 뜻이다 — 어느 쪽이든 알아야 한다.
"""
from __future__ import annotations

import math

import pytest

from simbot import spec as S
from simbot.spec import (Friction, Gate, Geometry, Medium, Numerics, PredictionItem,
                         Prediction, Q, Quantity, Species, SystemSpec, Timing,
                         ExternalField, PairInteraction)

EXAMPLE_SPEC = S.Path("examples/trap-2d-5um/spec.yaml")

# 손으로 쓴 03_spec.yaml §derived — **문서에 적힌 문자열 그대로** 옮긴다.
# 문자열로 두는 이유: 허용오차를 값마다 그 값의 유효숫자에서 뽑기 위해서다.
# `6.817e-6` 은 4자리이므로 반올림 반폭이 7.3e-5 이고, `4.14195e-21` 은 6자리이므로
# 7.2e-7 이다. 전역 상수 하나로 잡으면 느슨한 쪽이 엄격한 값을 봐주게 된다.
HAND_DERIVED: dict[str, str] = {
    "kT_si": "4.14195e-21",
    "D0_si": "5.1361e-14",
    "gamma_si": "8.0644e-8",
    "sigma_si": "1.0e-5",
    "tau_trap_si": "8.0644e-3",
    "tau_D_si": "1.9470e3",
    "l_trap_si": "2.0352e-8",
    "corner_freq_si": "19.735",
    "tau_inertial_si": "6.817e-6",
    "k_star_sigma": "2.4143e5",
}


def rounding_halfwidth(written: str) -> float:
    """문서에 `written` 으로 적힌 값의 반올림 반폭 (상대값).

    유효숫자 `n` 자리로 적힌 값 `v` 는 마지막 자리에서 ±0.5 만큼 잘렸을 수 있다.
    이것이 문서와 재계산을 비교할 때 **이론적으로 정당한 유일한 허용오차**다.
    """
    mant, _, exp = written.lower().partition("e")
    digits = len(mant.replace("-", "").replace(".", "").lstrip("0")) or 1
    v = abs(float(written))
    # 마지막 유효숫자의 크기 = v 의 10^(floor(log10 v) - (n-1))
    last_place = 10.0 ** (math.floor(math.log10(v)) - (digits - 1))
    return 0.5 * last_place / v


@pytest.fixture
def trap_spec() -> SystemSpec:
    if not EXAMPLE_SPEC.exists():
        pytest.skip("examples/trap-2d-5um/spec.yaml 없음")
    return SystemSpec.load(EXAMPLE_SPEC)


# =============================================================================
# Quantity — provenance 규약
# =============================================================================
def test_quantity_rejects_unknown_provenance():
    with pytest.raises(ValueError, match="provenance"):
        Quantity(value=1.0, provenance="vibes")


def test_quantity_rejects_unknown_confidence():
    with pytest.raises(ValueError, match="confidence"):
        Quantity(value=1.0, basis="x", confidence="pretty sure")


def test_quantity_without_basis_is_a_problem():
    assert any("basis" in p for p in Quantity(value=1.0).problems())


def test_assumed_without_confidence_is_a_problem():
    q = Quantity(value=1.0, provenance="assumed", basis="물이라고 가정")
    assert any("confidence" in p for p in q.problems())


def test_derived_does_not_require_confidence():
    q = Quantity(value=1.0, provenance="derived", basis="sigma = 2a")
    assert q.problems() == []


def test_cheap_model_cannot_write_inference_field():
    """master_plan §12.2 — inference/assumed 는 Opus 만."""
    q = Quantity(value=2, provenance="inference", basis="z 축이 없다",
                 confidence="medium", written_by="haiku")
    assert any("§12.2" in p for p in q.problems())


def test_cheap_model_may_write_observation_field():
    q = Quantity(value=1, provenance="observation", basis="원이 1개",
                 written_by="haiku")
    assert q.problems() == []


def test_si_property_refuses_non_numeric():
    with pytest.raises(TypeError):
        Q("periodic", provenance="rule", basis="x").si


def test_si_property_refuses_bool():
    """`True` 는 파이썬에서 1 이다 — 단위 산술에 흘러들면 조용히 틀린다."""
    with pytest.raises(TypeError):
        Q(True, provenance="rule", basis="x").si


# =============================================================================
# 게이트
# =============================================================================
def test_gate_rejects_unknown_status():
    with pytest.raises(ValueError):
        Gate(status="probably")


def test_unknown_gate_name_is_a_problem(trap_spec):
    """오타 난 게이트 이름은 '한 번도 실행되지 않는 검사'다."""
    trap_spec.gates["equiparition"] = Gate("required")     # 오타
    assert any("등록되지 않은 게이트" in p for p in S.validate(trap_spec).problems)


def test_gate_off_without_reason_is_a_problem(trap_spec):
    trap_spec.gates["equilibration_detection"] = Gate("off", "")
    assert any("이유가 없다" in p for p in S.validate(trap_spec).problems)


def test_required_gate_is_deferred_not_passed(trap_spec):
    """★ 선언은 결과가 아니다. `required` 를 pass 로 찍으면 아무도 보지 않은 합격."""
    rep = S.validate(trap_spec)
    names = {c.name: c.status for c in rep.checks}
    assert names["equipartition"] == "declared"
    assert names["em_bias_reproduced"] == "declared"
    assert "equipartition" in [c.name for c in rep.deferred()]


def test_computed_result_fills_declared_gate(trap_spec):
    """카드가 선언하고 S3 이 계산할 수 있는 게이트는 실제 판정이 들어간다."""
    rep = S.validate(trap_spec)
    c = next(c for c in rep.checks if c.name == "overdamped")
    assert c.status == "pass" and c.value is not None


def test_no_duplicate_check_rows(trap_spec):
    rep = S.validate(trap_spec)
    names = [c.name for c in rep.checks]
    assert len(names) == len(set(names))


# =============================================================================
# 파생값 — 첫 런 회귀
# =============================================================================
@pytest.mark.benchmark
def test_derive_reproduces_hand_written_first_run(trap_spec):
    """손문서의 파생값 10개를 전부 재현한다.

    허용오차는 **각 값이 문서에 적힌 유효숫자**에서 나온다 (conftest 규칙 1 —
    관측값에 맞춰 재단하지 않는다). 재계산이 반올림 반폭 안에 들어오면
    "문서와 코드가 같은 계산을 했다"가 증명된다.
    """
    d = S.derive(trap_spec)
    for key, written in HAND_DERIVED.items():
        hand = float(written)
        tol = rounding_halfwidth(written)
        rel = abs(d[key] - hand) / abs(hand)
        assert rel <= tol, (f"{key}: 재계산 {d[key]:.6g} vs 문서 {written} "
                            f"(상대차 {rel:.2e} > 반올림 반폭 {tol:.2e})")


def test_var_per_component_is_kT_over_k(trap_spec):
    """등분배 `⟨x²⟩ = kT/k` — 이 계의 1급 게이트가 보는 양."""
    d = S.derive(trap_spec)
    assert d["var_per_component_si"] == pytest.approx(d["kT_si"] / d["k_si"], rel=1e-15)
    assert d["msd_plateau_si"] == pytest.approx(
        2 * trap_spec.geometry.d * d["var_per_component_si"], rel=1e-15)


def test_rounding_halfwidth_matches_significant_figures():
    """허용오차 계산기 자체를 검증한다 — 틀리면 위 회귀가 무의미해진다."""
    assert rounding_halfwidth("6.817e-6") == pytest.approx(0.5e-3 / 6.817, rel=1e-9)
    assert rounding_halfwidth("4.14195e-21") == pytest.approx(0.5e-5 / 4.14195, rel=1e-9)
    assert rounding_halfwidth("19.735") == pytest.approx(0.5e-3 / 19.735, rel=1e-9)
    assert rounding_halfwidth("1.0e-5") == pytest.approx(0.5e-1 / 1.0, rel=1e-9)


def test_derive_gamma_uses_radius_not_diameter(trap_spec):
    """가장 흔한 실수: 6πηa 에 직경을 넣기. 넣으면 모든 시간척도가 2배 틀린다."""
    d = S.derive(trap_spec)
    a = trap_spec.primary.radius_si.si
    expected = 6.0 * math.pi * trap_spec.medium.eta_si.si * a
    assert d["gamma_si"] == pytest.approx(expected, rel=1e-12)
    assert d["gamma_si"] != pytest.approx(2 * expected, rel=1e-3)


def test_sigma_is_twice_the_radius(trap_spec):
    assert trap_spec.primary.sigma_si == pytest.approx(
        2 * trap_spec.primary.radius_si.si, rel=1e-15)


def test_reference_length_is_l_trap_for_trap_card(trap_spec):
    """★ 트랩 계의 기준 길이는 σ 가 아니라 ℓ_trap 이다 (카드 §3)."""
    d = S.derive(trap_spec)
    assert S.reference_length_si(trap_spec, d) == pytest.approx(d["l_trap_si"])
    assert d["l_trap_si"] / d["sigma_si"] < 1e-2      # 자기 지름의 0.2 %


def test_time_scale_separation_is_large(trap_spec):
    """이 계가 τ_D 규약을 기각하는 이유 — 실측 2.41e5 배."""
    d = S.derive(trap_spec)
    assert d["tau_D_si"] / d["tau_trap_si"] == pytest.approx(2.4143e5, rel=1e-3)


def test_derive_omits_trap_fields_when_no_trap(trap_spec):
    trap_spec.external = []
    d = S.derive(trap_spec)
    assert "tau_trap_si" not in d and "tau_D_si" in d


def test_derive_omits_inertia_without_density(trap_spec):
    trap_spec.primary.density_si = None
    d = S.derive(trap_spec)
    assert "tau_inertial_si" not in d
    assert next(c for c in S.validate(trap_spec).checks
                if c.name == "overdamped").status == "na"


# =============================================================================
# 파생값 대조 — 손으로 고친 숫자를 잡는다
# =============================================================================
HAND_DERIVED_F = {k: float(v) for k, v in HAND_DERIVED.items()}


def test_stored_derived_mismatch_is_caught(trap_spec):
    """★ 2026-07-28 의 kT 4번째 자리 오류 부류를 잡는 검사."""
    bad = dict(HAND_DERIVED_F, kT_si=4.1420e-21 * 1.01)      # 1 % 틀린 값
    rep = S.validate(trap_spec, stored_derived=bad)
    assert not rep.ok
    assert any("kT_si" in p for p in rep.problems)


def test_stored_derived_match_passes(trap_spec):
    rep = S.validate(trap_spec, stored_derived=HAND_DERIVED_F)
    assert rep.ok, rep.problems
    c = next(c for c in rep.checks if c.name == "derived_consistency")
    assert c.status == "pass"
    # 실제로 비교가 일어났는지 — 0건 비교로 통과하면 공허한 테스트다
    assert f"{len(HAND_DERIVED_F)}개 대조" in c.detail


def test_zero_comparisons_is_not_reported_as_pass(trap_spec):
    """★ 대조 0건을 pass 로 보고하면 '검사했다'가 거짓이 된다."""
    rep = S.validate(trap_spec, stored_derived={"unknown_key": 1.0,
                                                "note": "무한매질 Stokes"})
    c = next(c for c in rep.checks if c.name == "derived_consistency")
    assert c.status == "na"
    assert c.value == 0.0 and "대조 불가" in c.detail


def test_stored_derived_ignores_non_numeric(trap_spec):
    rep = S.validate(trap_spec, stored_derived={"note": "무한매질 Stokes"})
    assert rep.ok, rep.problems


# =============================================================================
# 타당성 검사
# =============================================================================
def test_example_spec_has_no_problems(trap_spec):
    rep = S.validate(trap_spec)
    assert rep.problems == []
    assert rep.failed() == []


def test_overdamped_fails_for_heavy_particle(trap_spec):
    trap_spec.primary.density_si = Q(1e9, "kg/m^3", "assumed", "비현실적으로 무거움",
                                     confidence="low")
    c = next(c for c in S.validate(trap_spec).checks if c.name == "overdamped")
    assert c.status == "fail" and "Langevin" in c.detail


def test_box_too_small_fails(trap_spec):
    trap_spec.geometry.box_over_ref = Q(3.0, "l_trap", "rule", "일부러 작게")
    c = next(c for c in S.validate(trap_spec).checks
             if c.name == "box_much_larger_than_l_trap")
    assert c.status == "fail"


def test_packing_fraction_is_off_without_pair_interactions(trap_spec):
    """★ 같은 트랩 안의 비상호작용 복제에 φ 를 걸면 φ=4741 로 통과 불가가 된다 —
    존재하지 않는 문제다. 값은 보여주되 판정하지 않는다."""
    c = next(c for c in S.validate(trap_spec).checks if c.name == "packing_fraction")
    assert c.status == "off"
    assert c.value is not None and c.value > 1.0        # 값 자체는 보고된다


def test_r_cut_gate_catches_minimum_image_violation(trap_spec):
    d = S.derive(trap_spec)
    box_min = 200.0 * d["l_trap_si"]
    trap_spec.pair = [PairInteraction(
        "probe", "probe", "wca",
        r_cut_si=Q(box_min, "m", "rule", "L/2 를 넘김 — 최소이미지 위반"))]
    trap_spec.gates.pop("r_cut_le_half_box")
    c = next(c for c in S.validate(trap_spec).checks if c.name == "r_cut_le_half_box")
    assert c.status == "fail"


def test_pair_without_r_cut_is_a_problem(trap_spec):
    trap_spec.pair = [PairInteraction("probe", "probe", "wca")]
    assert any("r_cut" in p for p in S.validate(trap_spec).problems)


def test_missing_box_is_a_problem(trap_spec):
    trap_spec.geometry.box_over_ref = None
    trap_spec.geometry.box_si = None
    assert any("박스" in p for p in S.validate(trap_spec).problems)


def test_explicit_box_si_is_used(trap_spec):
    trap_spec.geometry.box_over_ref = None
    trap_spec.geometry.box_si = Q([4.07e-6, 4.07e-6, 0.0], "m", "rule", "명시 박스")
    assert trap_spec.box_lengths_si()[0] == pytest.approx(4.07e-6)


# =============================================================================
# YAML 왕복 — 조용한 단위 손실을 잡는다
# =============================================================================
def test_yaml_roundtrip_is_exact(trap_spec):
    back = SystemSpec.from_yaml(trap_spec.to_yaml())
    assert S.to_dict(back) == S.to_dict(trap_spec)


def test_yaml_roundtrip_preserves_derived_bit_for_bit(trap_spec):
    """왕복 후 파생값이 **정확히** 같아야 한다 (master_plan §S4 게이트: < 1e-12)."""
    back = SystemSpec.from_yaml(trap_spec.to_yaml())
    d1, d2 = S.derive(trap_spec), S.derive(back)
    assert set(d1) == set(d2)
    for k in d1:
        assert abs(d1[k] - d2[k]) <= 1e-12 * abs(d1[k]), k


def test_hash_is_stable_across_roundtrip(trap_spec):
    assert SystemSpec.from_yaml(trap_spec.to_yaml()).hash() == trap_spec.hash()


def test_bonds_and_angles_survive_yaml_roundtrip(trap_spec):
    """결합·각의 `k_si` 는 차원이 다르다 (N/m vs J/rad²) — 왕복에서 섞이면 λ_max 가 틀린다."""
    from simbot.spec import AngleInteraction, BondInteraction
    trap_spec.bonds = [BondInteraction(params={
        "k_si": Q(1.0e-3, "N/m", "rule", "결합 스프링"),
        "r0_si": Q(1.0e-5, "m", "rule", "결합 길이 = σ")})]
    trap_spec.angles = [AngleInteraction(params={
        "k_si": Q(4.14e-17, "J/rad^2", "rule", "각 스프링"),
        "t0_rad": Q(math.pi, "rad", "rule", "직선 사슬")})]
    back = SystemSpec.from_yaml(trap_spec.to_yaml())
    assert S.to_dict(back) == S.to_dict(trap_spec)
    assert back.bond_stiffness_si() == trap_spec.bond_stiffness_si()
    assert back.angle_stiffness_si() == trap_spec.angle_stiffness_si()
    assert S.derive(back)["lambda_max_si"] == S.derive(trap_spec)["lambda_max_si"]


def test_empty_bonds_do_not_change_the_spec_hash(trap_spec):
    """★ 필드를 추가해도 결합이 없는 계의 해시는 그대로여야 한다 —
    아니면 기존 런의 `run_id`·캐시 키가 전부 갈린다."""
    text = trap_spec.to_yaml()
    assert "bonds:" not in text and "angles:" not in text


def test_hash_changes_with_physics(trap_spec):
    before = trap_spec.hash()
    trap_spec.medium.T_si = Q(310.0, "K", "from_drawing", "다른 온도")
    assert trap_spec.hash() != before


def test_yaml_always_writes_provenance(trap_spec):
    """`assumed` 가 기본값이라고 생략하면 '가정'과 '적기를 잊음'이 구별되지 않는다."""
    text = trap_spec.to_yaml()
    assert text.count("provenance:") == len(S._iter_quantities(trap_spec))
    assert "provenance: assumed" in text


def test_bare_value_loads_as_assumed_with_no_basis():
    """provenance 없이 맨 값을 쓰면 규약 위반으로 흘러가야 한다 (조용히 통과 금지)."""
    q = S.from_dict(Quantity, 300.0)
    assert q.provenance == "assumed"
    assert any("basis" in p for p in q.problems())


# =============================================================================
# 예측 (S2)
# =============================================================================
def test_empty_prediction_is_a_problem():
    assert any("0개" in p for p in Prediction(items=[]).problems())


def test_prediction_item_requires_tolerance_and_basis():
    it = PredictionItem(quantity="D*", value=1.0, tolerance="", basis="")
    probs = it.problems()
    assert any("tolerance" in p for p in probs) and any("basis" in p for p in probs)


def test_valid_prediction_has_no_problems():
    p = Prediction(items=[PredictionItem(
        quantity="var_x_star", value=1.00251, tolerance="±1%",
        basis="EM 편향 1/(1-dt*/2)", discriminates="적분기 스킴",
        competing_value=1.0)])
    assert p.problems() == []
