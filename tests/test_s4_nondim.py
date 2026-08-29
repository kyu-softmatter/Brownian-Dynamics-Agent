"""S4 — 무차원화 · 무차원수 · `dt` 선택.

두 가지를 지킨다:

1. **왕복 오차 < 1e-12** (master_plan §S4 게이트). 무차원화는 나눗셈 한 번이므로
   큰 오차는 계산 실수가 아니라 **규약 위반**의 신호다 (τ_D 로 나누고 τ_trap 으로
   되돌리는 식).

2. **`dt` 게이트가 계마다 다르게 켜진다.** 트랩 계에서 변위 게이트만 믿으면
   `Δt*` 상한이 완화시간 한계의 1000배가 된다 — 발산하지 않으므로 더 위험하다.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from simbot import nondim as N
from simbot.policy import Policy, deep_merge, load_policy
from simbot.spec import Q, SystemSpec, derive

EXAMPLE_SPEC = Path("examples/trap-2d-5um/spec.yaml")


@pytest.fixture
def trap_spec() -> SystemSpec:
    if not EXAMPLE_SPEC.exists():
        pytest.skip("examples/trap-2d-5um/spec.yaml 없음")
    return SystemSpec.load(EXAMPLE_SPEC)


@pytest.fixture
def policy() -> Policy:
    return load_policy()


# =============================================================================
# 카드가 척도를 소유한다
# =============================================================================
def test_trap_card_uses_l_trap_and_tau_trap(trap_spec):
    d = derive(trap_spec)
    sc = N.scales_for(trap_spec, d)
    assert sc.length_si == pytest.approx(d["l_trap_si"], rel=1e-15)
    assert sc.time_si == pytest.approx(d["tau_trap_si"], rel=1e-15)
    assert sc.energy_si == pytest.approx(d["kT_si"], rel=1e-15)


def test_trap_card_does_not_use_tau_D(trap_spec):
    """★ τ_D 를 골랐다면 시간 척도가 24만 배 어긋난다."""
    d = derive(trap_spec)
    sc = N.scales_for(trap_spec, d)
    assert sc.time_si / d["tau_D_si"] < 1e-4


def test_transport_card_uses_sigma_and_tau_D(trap_spec):
    """같은 계라도 목적동역학이 수송이면 척도가 (σ, τ_D) 다."""
    trap_spec.card = "passive-sphere--transport"
    d = derive(trap_spec)
    sc = N.scales_for(trap_spec, d)
    assert sc.length_si == pytest.approx(d["sigma_si"], rel=1e-15)
    assert sc.time_si == pytest.approx(d["tau_D_si"], rel=1e-15)


def test_unregistered_card_refuses_improvised_nondim(trap_spec):
    """★ 카드가 없는 쌍에 즉흥 무차원화 금지 (CLAUDE.md)."""
    trap_spec.card = "colloid--something-new"
    with pytest.raises(KeyError, match="즉흥 무차원화는 금지"):
        N.scales_for(trap_spec)


def test_error_message_points_at_the_template(trap_spec):
    trap_spec.card = "unknown--pair"
    with pytest.raises(KeyError, match="_TEMPLATE"):
        N.scales_for(trap_spec)


def test_abp_card_is_declared_unimplemented(trap_spec):
    """구현 안 된 것을 조용히 다른 척도로 대체하면 안 된다."""
    trap_spec.card = "abp--dense-collective"
    with pytest.raises(NotImplementedError):
        N.scales_for(trap_spec)


# =============================================================================
# 카드 단위의 정규화 — 정확히 1 이어야 한다
# =============================================================================
def test_trap_reduced_units_are_exactly_unity(trap_spec):
    """(ℓ_trap, kT, τ_trap) 하에서 `k* = D* = γ* = kT* = 1` 이 **정의상** 성립한다.

    1 에서 벗어나면 척도 정의가 틀렸다 — 부동소수점 오차(≲1e-15) 외의 편차는 버그다.
    """
    r = N.reduce_spec(trap_spec)
    for name, value in (("k_star", r.k_star), ("D_star", r.D_star),
                        ("gamma_star", r.gamma_star), ("kT_star", r.kT_star)):
        assert abs(value - 1.0) < 1e-14, f"{name} = {value!r}"


def test_equation_of_motion_has_no_free_parameters(trap_spec):
    """카드 §3 의 주장: `dr*/dt* = -r* + √2 ξ`. 파라미터가 없다."""
    r = N.reduce_spec(trap_spec)
    assert r.k_star == pytest.approx(1.0, abs=1e-14)
    assert r.D_star == pytest.approx(1.0, abs=1e-14)


def test_sigma_star_is_huge_in_trap_units(trap_spec):
    """★ 입자 지름이 기준길이의 491배다 — 입자는 자기 크기의 0.2 %만 움직인다."""
    r = N.reduce_spec(trap_spec)
    assert r.sigma_star == pytest.approx(491.358, rel=1e-4)
    assert 1.0 / r.sigma_star < 0.01


def test_step_counts_come_from_timing(trap_spec):
    r = N.reduce_spec(trap_spec)
    assert r.equil_steps == 2000          # 10 τ / 5e-3
    assert r.prod_steps == 8000           # 40 τ / 5e-3
    assert r.sample_interval_steps == 400  # 2 τ / 5e-3


# =============================================================================
# 왕복 — master_plan §S4 게이트
# =============================================================================
def test_roundtrip_error_below_gate(trap_spec):
    errs = N.roundtrip_errors(trap_spec)
    assert errs, "왕복 대상이 0개면 검사한 것이 없다"
    for name, e in errs.items():
        assert e < 1e-12, f"{name}: 상대오차 {e:.2e}"


def test_roundtrip_covers_every_scale_kind(trap_spec):
    errs = N.roundtrip_errors(trap_spec)
    for key in ("kT", "D0", "sigma", "gamma", "k", "box_x"):
        assert key in errs, key


def test_roundtrip_catches_wrong_time_scale(trap_spec):
    """규약을 어기면 왕복이 깨지는지 — 게이트가 실제로 무언가를 잡는지 확인."""
    d = derive(trap_spec)
    r = N.reduce_spec(trap_spec)
    from simbot.units import Scales
    broken = Scales(length_si=r.scales.length_si, energy_si=r.scales.energy_si,
                    time_si=d["tau_D_si"], origin="일부러 틀린 τ_D")
    object.__setattr__(r, "scales", broken)
    errs = N.roundtrip_errors(trap_spec, r)
    assert max(errs.values()) > 1e-3, "틀린 척도인데 왕복이 통과했다 — 게이트가 무력하다"


# =============================================================================
# dt 선택 — 계마다 켜지는 제약이 다르다
# =============================================================================
def test_relaxation_constraint_dominates_in_trap(trap_spec, policy):
    ch = N.choose_dt(trap_spec, policy=policy)
    assert ch.dominant == "relaxation_time"
    assert ch.dt_star == pytest.approx(0.01, rel=1e-9)


def test_em_bias_target_can_take_over(trap_spec, policy):
    """정확도 목표를 명시하면 그것이 지배한다 — 첫 런의 dt* = 5e-3 을 재현."""
    ch = N.choose_dt(trap_spec, policy=policy, target_em_bias=0.0025)
    assert ch.dominant == "em_bias_target"
    assert ch.dt_star == pytest.approx(5.0e-3, rel=5e-3)


@pytest.mark.benchmark
def test_displacement_gate_is_1000x_too_loose_for_a_trap(trap_spec, policy):
    """★ 변위 게이트만 믿으면 트랩 계에서 `Δt*` 상한이 완화시간 한계의 1000배다.

    `dt-gate-should-be-displacement-based` 는 **쌍 상호작용 계**에서 도출된 결론이고,
    트랩 계에서는 구속 완화시간이 진짜 제약이다. 두 게이트는 경쟁하지 않고
    보완한다 — 그것이 (계 × 목적동역학) 카드가 게이트 표를 소유하는 이유다.
    """
    ch = N.choose_dt(trap_spec, policy=policy)
    by = {c.name: c for c in ch.constraints}
    thermal = by["thermal_displacement"]
    relax = by["relaxation_time"]
    assert not thermal.active, "쌍 상호작용이 없는데 변위 게이트가 켜져 있다"
    ratio = thermal.dt_si_max / relax.dt_si_max
    assert ratio > 1e3, f"변위/완화 상한 비 = {ratio:.3g} (1000배 초과를 기대)"


def test_displacement_gate_turns_on_with_pair_interactions(trap_spec, policy):
    from simbot.spec import PairInteraction
    trap_spec.pair = [PairInteraction("probe", "probe", "wca",
                                      r_cut_si=Q(1.1e-5, "m", "rule", "WCA 2^(1/6)σ"))]
    ch = N.choose_dt(trap_spec, policy=policy)
    by = {c.name: c for c in ch.constraints}
    assert by["thermal_displacement"].active


def test_force_constraint_is_na_without_measured_force(trap_spec, policy):
    """★ `max|F|` 는 실제 힘 계산으로 얻어야 한다 — 추정 금지 (§5.4)."""
    from simbot.spec import PairInteraction
    trap_spec.pair = [PairInteraction("probe", "probe", "wca",
                                      r_cut_si=Q(1.1e-5, "m", "rule", "x"))]
    ch = N.choose_dt(trap_spec, policy=policy)
    force = next(c for c in ch.constraints if c.name == "force_displacement")
    assert force.dt_si_max is None
    assert "추정 금지" in force.basis


def test_force_constraint_applies_when_force_is_measured(trap_spec, policy):
    """힘 제약이 완화시간 제약을 이기는 문턱을 해석적으로 계산해서 넘긴다.

    `0.005 σ γ / F < ζ τ_trap`  ⟺  `F > 0.005 σ γ / (ζ τ_trap)`.
    관측값에 맞춰 고른 숫자가 아니라 제약식에서 나온 문턱이다.
    """
    from simbot.spec import PairInteraction
    trap_spec.pair = [PairInteraction("probe", "probe", "wca",
                                      r_cut_si=Q(1.1e-5, "m", "rule", "x"))]
    d = derive(trap_spec)
    ts = policy.timestep
    threshold = (ts["max_force_displacement_sigma"] * d["sigma_si"] * d["gamma_si"]
                 / (ts["relaxation_safety_factor"] * d["tau_trap_si"]))
    ch = N.choose_dt(trap_spec, policy=policy, max_force_si=2.0 * threshold)
    assert ch.dominant == "force_displacement"
    # 문턱 아래면 완화시간이 다시 이긴다 — 문턱이 실재함을 보인다
    ch2 = N.choose_dt(trap_spec, policy=policy, max_force_si=0.5 * threshold)
    assert ch2.dominant == "relaxation_time"


def test_dt_over_tau_D_is_logged_but_not_a_gate(trap_spec, policy):
    """`dt/τ_D` 는 기록만 한다 — 게이트로 쓰면 논문까지 나온 런을 기각한다."""
    ch = N.choose_dt(trap_spec, policy=policy)
    assert "dt_over_tau_D" in ch.logged
    assert not any(c.name == "dt_over_tau_D" for c in ch.constraints)


def test_universal_tau_D_convention_would_be_catastrophic(trap_spec):
    """★ 첫 런의 실측 재현: τ_D 기준 `dt* = 5e-5` 는 `12 τ_trap` 이 된다."""
    d = derive(trap_spec)
    dt_si = 5.0e-5 * d["tau_D_si"]
    assert dt_si / d["tau_trap_si"] == pytest.approx(12.07, rel=1e-2)


def test_choose_dt_refuses_when_no_constraint_applies(trap_spec, policy):
    """제약이 하나도 없으면 dt 를 고를 근거가 없다 — 조용히 기본값을 쓰지 않는다."""
    trap_spec.external = []                       # 트랩 제거 → 완화시간 없음
    trap_spec.card = "passive-sphere--transport"  # 쌍 상호작용도 없음
    with pytest.raises(ValueError, match="활성 제약이 하나도 없다"):
        N.choose_dt(trap_spec, policy=policy)


def test_dt_table_marks_the_dominant_constraint(trap_spec, policy):
    ch = N.choose_dt(trap_spec, policy=policy, target_em_bias=0.0025)
    table = ch.table()
    assert "**←지배**" in table
    assert table.count("←지배") == 1


# =============================================================================
# 결합 계 — 변위 게이트가 꺼지면 안 되고, 변위만으로는 부족하다
# =============================================================================
#  ★ 이 절이 없어서 두 가지가 조용히 통과했다:
#    ① `active=bool(spec.pair)` 라서 결합만 있는 계에서 변위 게이트가 꺼졌다
#    ② 곧은 사슬은 max|F*| = 0 인 정류점이라 힘 게이트도 무력하다 (kT=0 에서도 터진다)
#  근거: knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
def _chain(trap_spec, *, k_bond_star: float | None = 1.0e6,
           k_angle_star: float | None = 1.0e4, keep_trap: bool = False):
    """축약 강성 목표를 SI 로 되돌려 결합·각을 붙인다.

    `k_bond* = k_bond σ²/kT` · `k_angle* = k_angle/kT` (각은 J/rad² 이라 길이를 안 쓴다).
    결합 길이는 `σ` 로 둔다 — finding 의 실측 표가 `b* = 1` 에서 얻어졌다.
    """
    from simbot.spec import AngleInteraction, BondInteraction
    d = derive(trap_spec)
    kT, sigma = d["kT_si"], d["sigma_si"]
    if not keep_trap:
        trap_spec.external = []                       # 트랩 제거 → 완화시간 제약 없음
        trap_spec.card = "passive-sphere--transport"   # (σ, τ_D) 척도
    if k_bond_star is not None:
        trap_spec.bonds = [BondInteraction(params={
            "k_si": Q(k_bond_star * kT / sigma**2, "N/m", "rule", f"k* = {k_bond_star:g}"),
            "r0_si": Q(sigma, "m", "rule", "결합 길이 = σ")})]
    if k_angle_star is not None:
        trap_spec.angles = [AngleInteraction(params={
            "k_si": Q(k_angle_star * kT, "J/rad^2", "rule", f"k* = {k_angle_star:g}"),
            "t0_rad": Q(3.141592653589793, "rad", "rule", "직선 사슬")})]
    return trap_spec


def test_lambda_max_combines_bond_and_angle_stiffness(trap_spec):
    """`λ_max ≈ 4k_bond + 16k_angle/b²`. 각의 J/rad² 을 b² 로 나누지 않으면 차원이 틀린다."""
    spec = _chain(trap_spec, k_bond_star=1.0e6, k_angle_star=1.0e4)
    d = derive(spec)
    expect = 4.0 * d["k_bond_si"] + 16.0 * d["k_angle_si"] / d["bond_length_si"] ** 2
    assert d["lambda_max_si"] == pytest.approx(expect, rel=1e-15)
    # 축약으로 되돌리면 finding 의 λ_max* = 4k_b* + 16k_a* 와 같아야 한다
    lam_star = d["lambda_max_si"] * d["sigma_si"] ** 2 / d["kT_si"]
    assert lam_star == pytest.approx(4.0e6 + 16.0e4, rel=1e-12)


def test_displacement_gates_turn_on_for_a_bond_only_system(trap_spec, policy):
    """★ 쌍 상호작용이 없어도 결합 상대가 있으면 변위 게이트는 켜져야 한다."""
    spec = _chain(trap_spec)
    assert not spec.pair and spec.has_neighbor_interaction
    ch = N.choose_dt(spec, policy=policy)
    by = {c.name: c for c in ch.constraints}
    assert by["thermal_displacement"].active, "결합 계인데 확산 변위 게이트가 꺼졌다"
    assert by["force_displacement"].active


def test_stability_gate_dominates_for_stiff_bonds(trap_spec, policy):
    """지배 제약과 그 값이 **식에서** 나온다: `Δt* = s·2/(4k_b* + 16k_a*)`."""
    spec = _chain(trap_spec, k_bond_star=1.0e6, k_angle_star=1.0e4)
    ch = N.choose_dt(spec, policy=policy)
    assert ch.dominant == "stiff_stability"
    s = policy.timestep["stability_safety_factor"]
    assert ch.dt_star == pytest.approx(s * 2.0 / (4.0e6 + 16.0e4), rel=1e-12)
    # 선형 안정 한계는 dt/τ_stiff = 2 다. 지배하면 정확히 2s 여야 한다.
    assert ch.logged["dt_over_tau_stiff"] == pytest.approx(2.0 * s, rel=1e-12)


def test_force_gate_is_powerless_for_a_straight_chain(trap_spec, policy):
    """★ 실측 힘이 정확히 0 인 정류점 — 힘 게이트가 `None` 이어도 dt 는 정해져야 한다.

    이전 코드는 이 계에서 활성 제약을 하나도 못 찾아 `dt` 선택을 거부했고,
    그래서 스크립트가 게이트를 따로 재구현했다.
    """
    spec = _chain(trap_spec)
    ch = N.choose_dt(spec, policy=policy, max_force_si=0.0)
    by = {c.name: c for c in ch.constraints}
    assert by["force_displacement"].dt_si_max is None, "힘 0 인데 상한이 생겼다"
    assert by["force_displacement"].active                     # 무력하지만 켜져 있어야 한다
    assert ch.dominant == "stiff_stability"
    # "재봤더니 0" 과 "아직 안 재봤다" 가 표에서 구별되어야 한다
    assert "무력" in by["force_displacement"].basis
    assert "추정 금지" in N.choose_dt(spec, policy=policy).constraints[1].basis


def test_stiff_bonds_are_not_hidden_by_the_trap_relaxation_gate(trap_spec, policy):
    """★ 트랩 + 강성 결합 (큐의 `trap-drag-2d-hex300` 모양).

    문턱을 제약식에서 계산해서 양쪽으로 넘는다 — `s·2γ/λ = ζ·τ_trap` 인 지점.
    관측값에 맞춘 숫자가 아니다.
    """
    d0 = derive(trap_spec)
    ts = policy.timestep
    lam_crit = (ts["stability_safety_factor"] * 2.0 * d0["gamma_si"]
                / (ts["relaxation_safety_factor"] * d0["tau_trap_si"]))
    kT, sigma = d0["kT_si"], d0["sigma_si"]

    stiff = _chain(trap_spec, k_bond_star=2.0 * lam_crit / 4.0 * sigma**2 / kT,
                   k_angle_star=None, keep_trap=True)
    ch = N.choose_dt(stiff, policy=policy)
    assert ch.dominant == "stiff_stability", "결합 강성이 트랩 완화시간에 가려졌다"
    by = {c.name: c for c in ch.constraints}
    assert by["relaxation_time"].dt_si_max / ch.dt_si == pytest.approx(2.0, rel=1e-9)

    soft = _chain(trap_spec, k_bond_star=0.5 * lam_crit / 4.0 * sigma**2 / kT,
                  k_angle_star=None, keep_trap=True)
    assert N.choose_dt(soft, policy=policy).dominant == "relaxation_time"


def test_stability_gate_is_exempt_from_the_accuracy_floor(trap_spec, policy):
    """★ `k_bond* = 1e6` 의 안정성 상한 `9.6e-8` 은 `hard_floor = 1e-7` 아래다.

    floor 로 기각하면 finding 이 실제로 완주시킨 런을 막는다. 안정성 상한은 협상 대상이
    아니므로 통과시키되 **floor 아래라는 사실을 기록**한다 (비용 지레는 `k_bond`).
    정확도 제약이 floor 를 깨는 경우는 여전히 기각해야 한다.
    """
    spec = _chain(trap_spec, k_bond_star=1.0e6, k_angle_star=1.0e4)
    ch = N.choose_dt(spec, policy=policy)
    assert ch.dt_star < policy.timestep["hard_floor"]
    assert ch.logged["dt_star_below_hard_floor"] == pytest.approx(ch.dt_star, rel=1e-15)

    # 정확도 쪽으로 floor 를 깨면 기각된다 — 면제가 안정성에만 걸려 있음을 보인다
    d = derive(spec)
    huge = (policy.timestep["max_force_displacement_sigma"] * d["sigma_si"]
            * d["gamma_si"] / (1e-3 * policy.timestep["hard_floor"] * d["tau_D_si"]))
    with pytest.raises(ValueError, match="hard_floor"):
        N.choose_dt(spec, policy=policy, max_force_si=huge)


#  finding 의 이분법 실측표 (노이즈 0 · eps=1e-3 · 4000 step · 결합 1±5 % 를 "안정"으로)
STABILITY_MEASURED = [
    (1.0e6, 1.0e4, 5, 1.00e-6), (1.0e6, 1.0e4, 9, 5.87e-7),
    (1.0e5, 1.0e4, 5, 1.00e-5), (1.0e5, 1.0e4, 9, 5.87e-6),
    (1.0e4, 1.0e4, 5, 1.84e-5), (1.0e4, 1.0e4, 9, 1.48e-5),
    (1.0e4, 1.0e3, 5, 1.00e-4), (1.0e4, 1.0e3, 9, 5.87e-5),
    (1.0e3, 1.0e3, 5, 1.84e-4), (1.0e3, 1.0e3, 9, 1.48e-4),
]


@pytest.mark.benchmark
@pytest.mark.parametrize("k_bond_star,k_angle_star,n,dt_crit", STABILITY_MEASURED,
                         ids=[f"kb{a:g}_ka{b:g}_N{c}" for a, b, c, _ in STABILITY_MEASURED])
def test_stability_gate_stays_below_measured_critical_dt(
        trap_spec, policy, k_bond_star, k_angle_star, n, dt_crit):
    """★ 게이트가 실측 임계 `Δt` 아래에 있고, 근거 없이 보수적이지도 않다.

    안전계수 `0.2` 의 유일한 정당화가 이 표다 (`reproduced: yes`). 여유가 1 미만이면
    터지는 런을 통과시키고, 15 배를 넘으면 스텝 수를 근거 없이 늘린다.
    finding 이 주장한 범위는 `6–14` 배이고 `(1e5,1e4)`·`(1e4,1e3)` 행이 상한 `14.0` 이다.
    `n` 은 게이트에 들어가지 않는다 — λ_max 근사가 사슬 길이에 무관함을 이 표가 확인한다.
    """
    spec = _chain(trap_spec, k_bond_star=k_bond_star, k_angle_star=k_angle_star)
    ch = N.choose_dt(spec, policy=policy)
    assert ch.dominant == "stiff_stability"
    margin = dt_crit / ch.dt_star
    assert margin > 1.0, f"게이트가 실측 임계값 위다 (여유 {margin:.2f}배)"
    assert 6.0 <= margin <= 15.0, f"여유 {margin:.2f}배 — finding 의 6–14배 범위를 벗어났다"


# =============================================================================
# 무차원수
# =============================================================================
def test_groups_reproduce_first_run_values(trap_spec):
    g = N.groups(trap_spec)
    assert g["k_star_sigma"] == pytest.approx(2.4143e5, rel=1e-3)
    assert g["l_trap_over_sigma"] == pytest.approx(2.0352e-3, rel=1e-3)
    assert g["tau_D_over_tau_trap"] == pytest.approx(2.4143e5, rel=1e-3)
    assert g["reynolds"] == pytest.approx(1.470e-5, rel=1e-2)
    assert g["tau_inertial_over_ref_time"] == pytest.approx(8.454e-4, rel=1e-3)


def test_phi_absent_without_pair_interactions(trap_spec):
    """★ 계산 불가능한 것을 0 으로 채우면 '이 값은 0 이다'로 읽힌다."""
    assert "phi" not in N.groups(trap_spec)


def test_k_star_sigma_and_k_star_are_different_numbers(trap_spec):
    """`kσ²/kT` 와 카드 단위 `k*` 를 혼동하면 5 decade 틀린다."""
    g = N.groups(trap_spec)
    assert g["k_star"] == pytest.approx(1.0, abs=1e-14)
    assert g["k_star_sigma"] > 1e5


# =============================================================================
# 정책
# =============================================================================
def test_deep_merge_preserves_siblings():
    base = {"tiers": {"production": {"N": 4000, "steps": 1_000_000}}}
    out = deep_merge(base, {"tiers": {"production": {"N": 8000}}})
    assert out["tiers"]["production"] == {"N": 8000, "steps": 1_000_000}


def test_deep_merge_replaces_lists_wholesale():
    """사다리를 **줄이려는** override 가 늘어나면 안 된다."""
    base = {"budget": {"required_tier_ladder": ["smoke", "pilot", "explore"]}}
    out = deep_merge(base, {"budget": {"required_tier_ladder": ["pilot"]}})
    assert out["budget"]["required_tier_ladder"] == ["pilot"]


def test_deep_merge_does_not_mutate_base():
    base = {"a": {"b": 1}}
    deep_merge(base, {"a": {"b": 2}})
    assert base == {"a": {"b": 1}}


def test_policy_loads_measured_constants(policy):
    assert policy.get("hardware.throughput_particle_steps_per_s") == pytest.approx(6.3e6)
    assert policy.seeds_default == 4


def test_policy_rejects_number_parsed_as_string(tmp_path):
    """★ YAML 1.1: `6.3e6` 은 문자열이다. 실제로 처리량 상수가 그랬다 (2026-07-28)."""
    p = tmp_path / "policy.yaml"
    p.write_text("hardware:\n  throughput_particle_steps_per_s: 6.3e6\n",
                 encoding="utf-8")
    with pytest.raises(ValueError, match="문자열로 파싱"):
        load_policy(p)


def test_real_policy_has_no_string_numbers(policy):
    """정책 파일 전체를 감시한다 — 새 값이 들어올 때 이 테스트가 잡는다."""
    from simbot.policy import _find_numeric_strings
    assert _find_numeric_strings(policy.raw) == []


def test_concurrency_respects_hard_max(policy):
    """k=12 에서 총처리량 회귀가 측정됐다 — 넘으면 안 된다."""
    assert policy.concurrency("batch") <= policy.get("concurrency.hard_max")


def test_concurrency_clamps_absurd_override(tmp_path):
    p = tmp_path / "policy.yaml"
    p.write_text("concurrency:\n  default: 8\n  hard_max: 10\n"
                 "overrides:\n  concurrency:\n    default: 64\n", encoding="utf-8")
    assert load_policy(p).concurrency() == 10


def test_efficiency_table_is_a_step_function(policy):
    assert policy.efficiency(4) == pytest.approx(0.925)
    assert policy.efficiency(7) == pytest.approx(policy.efficiency(6))
    assert policy.efficiency(1) == pytest.approx(1.0)


def test_unknown_tier_raises_with_options(policy):
    with pytest.raises(KeyError, match="있는 것"):
        policy.tier("turbo")


def test_overridden_paths_are_reported(tmp_path):
    """`params` 명령의 ⚠(아무도 안 고른 기본값) 표시가 이 목록에 의존한다."""
    p = tmp_path / "policy.yaml"
    p.write_text("seeds:\n  default: 4\n"
                 "overrides:\n  seeds:\n    default: 8\n  _why: 'INCONCLUSIVE 보강'\n",
                 encoding="utf-8")
    pol = load_policy(p)
    assert pol.seeds_default == 8
    assert pol.overridden_paths == ["seeds.default"]      # `_why` 는 제외


def test_real_policy_has_no_overrides_yet(policy):
    """현재 정책에 사람 override 가 없음을 고정 — 생기면 이 테스트가 알려준다."""
    assert policy.overridden_paths == []


def test_cost_constants_agree_between_code_and_policy(policy):
    """★ 처리량 상수가 `estimators.py` 와 `run_policy.yaml` 두 곳에 있다.

    어긋나면 어느 경로로 추정했는지에 따라 비용이 달라지고, 예산 게이트가
    조용히 다른 문턱을 쓰게 된다.
    """
    from simbot.estimators import EFFICIENCY_BY_K, THROUGHPUT_PARTICLE_STEPS_PER_S
    assert THROUGHPUT_PARTICLE_STEPS_PER_S == pytest.approx(
        policy.get("hardware.throughput_particle_steps_per_s"), rel=1e-12)
    yaml_eff = {int(k): float(v)
                for k, v in policy.get("hardware.efficiency_by_k").items()}
    assert EFFICIENCY_BY_K == pytest.approx(yaml_eff)
