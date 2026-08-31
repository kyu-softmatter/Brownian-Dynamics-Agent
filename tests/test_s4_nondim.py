"""S4 — non-dimensionalization · dimensionless groups · choosing `dt`.

Two things are enforced:

1. **round-trip error < 1e-12** (the master_plan §S4 gate). Non-dimensionalization
   is one division, so a large error is not an arithmetic slip but a signal of a
   **broken convention** (dividing by τ_D and converting back with τ_trap, say).

2. **The `dt` gates turn on differently per system.** Trust only the displacement
   gate in a trap system and the `Δt*` cap comes out 1000x the relaxation-time
   limit — which is more dangerous because it does not diverge.
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
        pytest.skip("examples/trap-2d-5um/spec.yaml is absent")
    return SystemSpec.load(EXAMPLE_SPEC)


@pytest.fixture
def policy() -> Policy:
    return load_policy()


# =============================================================================
# the card owns the scales
# =============================================================================
def test_trap_card_uses_l_trap_and_tau_trap(trap_spec):
    d = derive(trap_spec)
    sc = N.scales_for(trap_spec, d)
    assert sc.length_si == pytest.approx(d["l_trap_si"], rel=1e-15)
    assert sc.time_si == pytest.approx(d["tau_trap_si"], rel=1e-15)
    assert sc.energy_si == pytest.approx(d["kT_si"], rel=1e-15)


def test_trap_card_does_not_use_tau_D(trap_spec):
    """★ Had τ_D been chosen, the time scale would be off by 240,000x."""
    d = derive(trap_spec)
    sc = N.scales_for(trap_spec, d)
    assert sc.time_si / d["tau_D_si"] < 1e-4


def test_transport_card_uses_sigma_and_tau_D(trap_spec):
    """For the same system, if the target dynamics is transport the scales are
    (σ, τ_D)."""
    trap_spec.card = "passive-sphere--transport"
    d = derive(trap_spec)
    sc = N.scales_for(trap_spec, d)
    assert sc.length_si == pytest.approx(d["sigma_si"], rel=1e-15)
    assert sc.time_si == pytest.approx(d["tau_D_si"], rel=1e-15)


def test_unregistered_card_refuses_improvised_nondim(trap_spec):
    """★ No improvised non-dimensionalization for an unregistered pair
    (CLAUDE.md)."""
    trap_spec.card = "colloid--something-new"
    with pytest.raises(KeyError, match="Improvised"):
        N.scales_for(trap_spec)


def test_error_message_points_at_the_template(trap_spec):
    trap_spec.card = "unknown--pair"
    with pytest.raises(KeyError, match="_TEMPLATE"):
        N.scales_for(trap_spec)


def test_abp_card_is_declared_unimplemented(trap_spec):
    """Something unimplemented must not be quietly substituted with other scales."""
    trap_spec.card = "abp--dense-collective"
    with pytest.raises(NotImplementedError):
        N.scales_for(trap_spec)


# =============================================================================
# normalization in card units — must be exactly 1
# =============================================================================
def test_trap_reduced_units_are_exactly_unity(trap_spec):
    """Under (ℓ_trap, kT, τ_trap), `k* = D* = γ* = kT* = 1` holds **by definition.**

    A departure from 1 means the scale definition is wrong -- anything beyond
    floating-point error (≲1e-15) is a bug.
    """
    r = N.reduce_spec(trap_spec)
    for name, value in (("k_star", r.k_star), ("D_star", r.D_star),
                        ("gamma_star", r.gamma_star), ("kT_star", r.kT_star)):
        assert abs(value - 1.0) < 1e-14, f"{name} = {value!r}"


def test_equation_of_motion_has_no_free_parameters(trap_spec):
    """Card §3's claim: `dr*/dt* = -r* + √2 ξ`. No parameters left."""
    r = N.reduce_spec(trap_spec)
    assert r.k_star == pytest.approx(1.0, abs=1e-14)
    assert r.D_star == pytest.approx(1.0, abs=1e-14)


def test_sigma_star_is_huge_in_trap_units(trap_spec):
    """★ The particle diameter is 491x the reference length — the particle moves
    only 0.2 % of its own size."""
    r = N.reduce_spec(trap_spec)
    assert r.sigma_star == pytest.approx(491.358, rel=1e-4)
    assert 1.0 / r.sigma_star < 0.01


def test_step_counts_come_from_timing(trap_spec):
    r = N.reduce_spec(trap_spec)
    assert r.equil_steps == 2000          # 10 τ / 5e-3
    assert r.prod_steps == 8000           # 40 τ / 5e-3
    assert r.sample_interval_steps == 400  # 2 τ / 5e-3


# =============================================================================
# the round-trip — the master_plan §S4 gate
# =============================================================================
def test_roundtrip_error_below_gate(trap_spec):
    errs = N.roundtrip_errors(trap_spec)
    assert errs, "0 round-trip targets means nothing was checked"
    for name, e in errs.items():
        assert e < 1e-12, f"{name}: relative error {e:.2e}"


def test_roundtrip_covers_every_scale_kind(trap_spec):
    errs = N.roundtrip_errors(trap_spec)
    for key in ("kT", "D0", "sigma", "gamma", "k", "box_x"):
        assert key in errs, key


def test_roundtrip_catches_wrong_time_scale(trap_spec):
    """Does breaking the convention break the round-trip -- i.e. does the gate
    actually catch anything."""
    d = derive(trap_spec)
    r = N.reduce_spec(trap_spec)
    from simbot.units import Scales
    broken = Scales(length_si=r.scales.length_si, energy_si=r.scales.energy_si,
                    time_si=d["tau_D_si"], origin="deliberately wrong τ_D")
    object.__setattr__(r, "scales", broken)
    errs = N.roundtrip_errors(trap_spec, r)
    assert max(errs.values()) > 1e-3, \
        "the scales are wrong and the round-trip passed — the gate is toothless"


# =============================================================================
# choosing dt — different constraints turn on per system
# =============================================================================
def test_relaxation_constraint_dominates_in_trap(trap_spec, policy):
    ch = N.choose_dt(trap_spec, policy=policy)
    assert ch.dominant == "relaxation_time"
    assert ch.dt_star == pytest.approx(0.01, rel=1e-9)


def test_em_bias_target_can_take_over(trap_spec, policy):
    """State an accuracy target and it dominates — reproduces the first run's
    dt* = 5e-3."""
    ch = N.choose_dt(trap_spec, policy=policy, target_em_bias=0.0025)
    assert ch.dominant == "em_bias_target"
    assert ch.dt_star == pytest.approx(5.0e-3, rel=5e-3)


@pytest.mark.benchmark
def test_displacement_gate_is_1000x_too_loose_for_a_trap(trap_spec, policy):
    """★ Trust only the displacement gate and a trap system's `Δt*` cap is 1000x
    the relaxation-time limit.

    `dt-gate-should-be-displacement-based` is a conclusion drawn from **pair-
    interacting systems**; in a trap system the confinement relaxation time is the
    real constraint. The two gates do not compete, they complement -- which is why
    the (system × target dynamics) card owns the gate table.
    """
    ch = N.choose_dt(trap_spec, policy=policy)
    by = {c.name: c for c in ch.constraints}
    thermal = by["thermal_displacement"]
    relax = by["relaxation_time"]
    assert not thermal.active, \
        "no pair interaction and yet the displacement gate is on"
    ratio = thermal.dt_si_max / relax.dt_si_max
    assert ratio > 1e3, \
        f"displacement/relaxation cap ratio = {ratio:.3g} (expected over 1000x)"


def test_displacement_gate_turns_on_with_pair_interactions(trap_spec, policy):
    from simbot.spec import PairInteraction
    trap_spec.pair = [PairInteraction("probe", "probe", "wca",
                                      r_cut_si=Q(1.1e-5, "m", "rule", "WCA 2^(1/6)σ"))]
    ch = N.choose_dt(trap_spec, policy=policy)
    by = {c.name: c for c in ch.constraints}
    assert by["thermal_displacement"].active


def test_force_constraint_is_na_without_measured_force(trap_spec, policy):
    """★ `max|F|` has to come from an actual force computation — estimating is
    forbidden (§5.4)."""
    from simbot.spec import PairInteraction
    trap_spec.pair = [PairInteraction("probe", "probe", "wca",
                                      r_cut_si=Q(1.1e-5, "m", "rule", "x"))]
    ch = N.choose_dt(trap_spec, policy=policy)
    force = next(c for c in ch.constraints if c.name == "force_displacement")
    assert force.dt_si_max is None
    assert "estimating is forbidden" in force.basis


def test_force_constraint_applies_when_force_is_measured(trap_spec, policy):
    """Crosses the threshold where the force constraint beats the relaxation-time
    one, computed analytically.

    `0.005 σ γ / F < ζ τ_trap`  ⟺  `F > 0.005 σ γ / (ζ τ_trap)`.
    A threshold from the constraint equation, not a number picked to fit the
    observation.
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
    # below the threshold relaxation wins again — showing the threshold is real
    ch2 = N.choose_dt(trap_spec, policy=policy, max_force_si=0.5 * threshold)
    assert ch2.dominant == "relaxation_time"


def test_dt_over_tau_D_is_logged_but_not_a_gate(trap_spec, policy):
    """`dt/τ_D` is only recorded — used as a gate it rejects runs that made it into
    papers."""
    ch = N.choose_dt(trap_spec, policy=policy)
    assert "dt_over_tau_D" in ch.logged
    assert not any(c.name == "dt_over_tau_D" for c in ch.constraints)


def test_universal_tau_D_convention_would_be_catastrophic(trap_spec):
    """★ Reproduces the first run's measurement: `dt* = 5e-5` on the τ_D basis
    comes out as `12 τ_trap`."""
    d = derive(trap_spec)
    dt_si = 5.0e-5 * d["tau_D_si"]
    assert dt_si / d["tau_trap_si"] == pytest.approx(12.07, rel=1e-2)


def test_choose_dt_refuses_when_no_constraint_applies(trap_spec, policy):
    """With no constraint at all there is no basis for choosing dt — it does not
    quietly fall back to a default."""
    trap_spec.external = []                       # remove the trap → no relaxation
    trap_spec.card = "passive-sphere--transport"  # no pair interaction either
    with pytest.raises(ValueError, match="not a single active constraint"):
        N.choose_dt(trap_spec, policy=policy)


def test_dt_table_marks_the_dominant_constraint(trap_spec, policy):
    ch = N.choose_dt(trap_spec, policy=policy, target_em_bias=0.0025)
    table = ch.table()
    assert "**←dominant**" in table
    assert table.count("←dominant") == 1


# =============================================================================
# bonded systems — the displacement gate must not turn off, and displacement alone
# is not enough
# =============================================================================
#  ★ Without this section two things passed quietly:
#    ① `active=bool(spec.pair)` turned the displacement gate off in a bond-only
#       system
#    ② a straight chain is a stationary point with max|F*| = 0, so the force gate
#       is toothless too (it blows up even at kT=0)
#  Basis: knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
def _chain(trap_spec, *, k_bond_star: float | None = 1.0e6,
           k_angle_star: float | None = 1.0e4, keep_trap: bool = False):
    """Convert reduced stiffness targets back to SI and attach bonds and angles.

    `k_bond* = k_bond σ²/kT` · `k_angle* = k_angle/kT` (an angle is J/rad², so it
    uses no length). The bond length is set to `σ` -- the finding's measured table
    was obtained at `b* = 1`.
    """
    from simbot.spec import AngleInteraction, BondInteraction
    d = derive(trap_spec)
    kT, sigma = d["kT_si"], d["sigma_si"]
    if not keep_trap:
        trap_spec.external = []                       # remove the trap → no
                                                      # relaxation constraint
        trap_spec.card = "passive-sphere--transport"   # (σ, τ_D) scales
    if k_bond_star is not None:
        trap_spec.bonds = [BondInteraction(params={
            "k_si": Q(k_bond_star * kT / sigma**2, "N/m", "rule", f"k* = {k_bond_star:g}"),
            "r0_si": Q(sigma, "m", "rule", "bond length = σ")})]
    if k_angle_star is not None:
        trap_spec.angles = [AngleInteraction(params={
            "k_si": Q(k_angle_star * kT, "J/rad^2", "rule", f"k* = {k_angle_star:g}"),
            "t0_rad": Q(3.141592653589793, "rad", "rule", "a straight chain")})]
    return trap_spec


def test_lambda_max_combines_bond_and_angle_stiffness(trap_spec):
    """`λ_max ≈ 4k_bond + 16k_angle/b²`. Fail to divide the angle's J/rad² by b² and
    the dimensions are wrong."""
    spec = _chain(trap_spec, k_bond_star=1.0e6, k_angle_star=1.0e4)
    d = derive(spec)
    expect = 4.0 * d["k_bond_si"] + 16.0 * d["k_angle_si"] / d["bond_length_si"] ** 2
    assert d["lambda_max_si"] == pytest.approx(expect, rel=1e-15)
    # converted back to reduced it must equal the finding's λ_max* = 4k_b* + 16k_a*
    lam_star = d["lambda_max_si"] * d["sigma_si"] ** 2 / d["kT_si"]
    assert lam_star == pytest.approx(4.0e6 + 16.0e4, rel=1e-12)


def test_displacement_gates_turn_on_for_a_bond_only_system(trap_spec, policy):
    """★ With no pair interaction but a bonded partner, the displacement gate must
    still be on."""
    spec = _chain(trap_spec)
    assert not spec.pair and spec.has_neighbor_interaction
    ch = N.choose_dt(spec, policy=policy)
    by = {c.name: c for c in ch.constraints}
    assert by["thermal_displacement"].active, \
        "a bonded system with the diffusive displacement gate off"
    assert by["force_displacement"].active


def test_stability_gate_dominates_for_stiff_bonds(trap_spec, policy):
    """The dominant constraint and its value come **from the equation**:
    `Δt* = s·2/(4k_b* + 16k_a*)`."""
    spec = _chain(trap_spec, k_bond_star=1.0e6, k_angle_star=1.0e4)
    ch = N.choose_dt(spec, policy=policy)
    assert ch.dominant == "stiff_stability"
    s = policy.timestep["stability_safety_factor"]
    assert ch.dt_star == pytest.approx(s * 2.0 / (4.0e6 + 16.0e4), rel=1e-12)
    # the linear stability limit is dt/τ_stiff = 2. When it dominates, exactly 2s.
    assert ch.logged["dt_over_tau_stiff"] == pytest.approx(2.0 * s, rel=1e-12)


def test_force_gate_is_powerless_for_a_straight_chain(trap_spec, policy):
    """★ A stationary point where the measured force is exactly 0 — dt must still be
    determined even with the force gate at `None`.

    The earlier code found no active constraint at all for this system and refused
    to choose `dt`, which is why a script re-implemented the gates separately.
    """
    spec = _chain(trap_spec)
    ch = N.choose_dt(spec, policy=policy, max_force_si=0.0)
    by = {c.name: c for c in ch.constraints}
    assert by["force_displacement"].dt_si_max is None, \
        "the force is 0 and yet a cap appeared"
    assert by["force_displacement"].active            # toothless but must stay on
    assert ch.dominant == "stiff_stability"
    # "measured it and got 0" and "have not measured it" must be distinguishable
    assert "toothless" in by["force_displacement"].basis
    assert "estimating is forbidden" in \
        N.choose_dt(spec, policy=policy).constraints[1].basis


def test_stiff_bonds_are_not_hidden_by_the_trap_relaxation_gate(trap_spec, policy):
    """★ Trap + stiff bonds (the shape of the queue's `trap-drag-2d-hex300`).

    The threshold is computed from the constraint equation and crossed from both
    sides -- the point where `s·2γ/λ = ζ·τ_trap`. Not a number fitted to the
    observation.
    """
    d0 = derive(trap_spec)
    ts = policy.timestep
    lam_crit = (ts["stability_safety_factor"] * 2.0 * d0["gamma_si"]
                / (ts["relaxation_safety_factor"] * d0["tau_trap_si"]))
    kT, sigma = d0["kT_si"], d0["sigma_si"]

    stiff = _chain(trap_spec, k_bond_star=2.0 * lam_crit / 4.0 * sigma**2 / kT,
                   k_angle_star=None, keep_trap=True)
    ch = N.choose_dt(stiff, policy=policy)
    assert ch.dominant == "stiff_stability", \
        "the bond stiffness was hidden by the trap relaxation time"
    by = {c.name: c for c in ch.constraints}
    assert by["relaxation_time"].dt_si_max / ch.dt_si == pytest.approx(2.0, rel=1e-9)

    soft = _chain(trap_spec, k_bond_star=0.5 * lam_crit / 4.0 * sigma**2 / kT,
                  k_angle_star=None, keep_trap=True)
    assert N.choose_dt(soft, policy=policy).dominant == "relaxation_time"


def test_stability_gate_is_exempt_from_the_accuracy_floor(trap_spec, policy):
    """★ At `k_bond* = 1e6` the stability cap `9.6e-8` is below `hard_floor = 1e-7`.

    Rejecting on the floor blocks a run the finding actually completed. A stability
    cap is not negotiable, so it passes but **the fact that it is below the floor is
    recorded** (the cost lever is `k_bond`). An accuracy constraint breaking the
    floor still has to be rejected.
    """
    spec = _chain(trap_spec, k_bond_star=1.0e6, k_angle_star=1.0e4)
    ch = N.choose_dt(spec, policy=policy)
    assert ch.dt_star < policy.timestep["hard_floor"]
    assert ch.logged["dt_star_below_hard_floor"] == pytest.approx(ch.dt_star, rel=1e-15)

    # breaking the floor on the accuracy side is rejected — showing the exemption
    # applies to stability only
    d = derive(spec)
    huge = (policy.timestep["max_force_displacement_sigma"] * d["sigma_si"]
            * d["gamma_si"] / (1e-3 * policy.timestep["hard_floor"] * d["tau_D_si"]))
    with pytest.raises(ValueError, match="hard_floor"):
        N.choose_dt(spec, policy=policy, max_force_si=huge)


#  The finding's bisection measurement table (noise 0 · eps=1e-3 · 4000 steps ·
#  bond length within 1±5 % counted as "stable")
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
    """★ The gate sits below the measured critical `Δt`, and is not conservative
    without a basis.

    This table is the only justification for the safety factor `0.2`
    (`reproduced: yes`). A margin under 1 would pass a run that blows up; over 15x
    would inflate the step count with no basis. The range the finding claims is
    `6–14`x, with the `(1e5,1e4)` and `(1e4,1e3)` rows at the `14.0` upper end.
    `n` does not enter the gate -- this table confirms the λ_max approximation is
    independent of chain length.
    """
    spec = _chain(trap_spec, k_bond_star=k_bond_star, k_angle_star=k_angle_star)
    ch = N.choose_dt(spec, policy=policy)
    assert ch.dominant == "stiff_stability"
    margin = dt_crit / ch.dt_star
    assert margin > 1.0, f"the gate is above the measured critical value " \
                         f"(margin {margin:.2f}x)"
    assert 6.0 <= margin <= 15.0, \
        f"margin {margin:.2f}x — outside the finding's 6–14x range"


# =============================================================================
# dimensionless groups
# =============================================================================
def test_groups_reproduce_first_run_values(trap_spec):
    g = N.groups(trap_spec)
    assert g["k_star_sigma"] == pytest.approx(2.4143e5, rel=1e-3)
    assert g["l_trap_over_sigma"] == pytest.approx(2.0352e-3, rel=1e-3)
    assert g["tau_D_over_tau_trap"] == pytest.approx(2.4143e5, rel=1e-3)
    assert g["reynolds"] == pytest.approx(1.470e-5, rel=1e-2)
    assert g["tau_inertial_over_ref_time"] == pytest.approx(8.454e-4, rel=1e-3)


def test_phi_absent_without_pair_interactions(trap_spec):
    """★ Fill something uncomputable with 0 and it reads as 'this value is 0'."""
    assert "phi" not in N.groups(trap_spec)


def test_k_star_sigma_and_k_star_are_different_numbers(trap_spec):
    """Confuse `kσ²/kT` with the card-unit `k*` and it is off by 5 decades."""
    g = N.groups(trap_spec)
    assert g["k_star"] == pytest.approx(1.0, abs=1e-14)
    assert g["k_star_sigma"] > 1e5


# =============================================================================
# policy
# =============================================================================
def test_deep_merge_preserves_siblings():
    base = {"tiers": {"production": {"N": 4000, "steps": 1_000_000}}}
    out = deep_merge(base, {"tiers": {"production": {"N": 8000}}})
    assert out["tiers"]["production"] == {"N": 8000, "steps": 1_000_000}


def test_deep_merge_replaces_lists_wholesale():
    """An override meant to **shorten** the ladder must not lengthen it."""
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
    """★ YAML 1.1: `6.3e6` is a string. The throughput constant actually was one
    (2026-07-28)."""
    p = tmp_path / "policy.yaml"
    p.write_text("hardware:\n  throughput_particle_steps_per_s: 6.3e6\n",
                 encoding="utf-8")
    with pytest.raises(ValueError, match="parsed as a string"):
        load_policy(p)


def test_real_policy_has_no_string_numbers(policy):
    """Watches the whole policy file — this test catches a new value as it lands."""
    from simbot.policy import _find_numeric_strings
    assert _find_numeric_strings(policy.raw) == []


def test_concurrency_respects_hard_max(policy):
    """A total-throughput regression was measured at k=12 — it must not be exceeded."""
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
    with pytest.raises(KeyError, match="available"):
        policy.tier("turbo")


def test_overridden_paths_are_reported(tmp_path):
    """The `params` command's ⚠ (a default nobody picked) depends on this list."""
    p = tmp_path / "policy.yaml"
    p.write_text("seeds:\n  default: 4\n"
                 "overrides:\n  seeds:\n    default: 8\n"
                 "  _why: 'reinforcing an INCONCLUSIVE'\n",
                 encoding="utf-8")
    pol = load_policy(p)
    assert pol.seeds_default == 8
    assert pol.overridden_paths == ["seeds.default"]      # `_why` excluded


def test_real_policy_has_no_overrides_yet(policy):
    """Pins that the current policy has no human override — this test tells you when
    one appears."""
    assert policy.overridden_paths == []


def test_cost_constants_agree_between_code_and_policy(policy):
    """★ The throughput constant lives in two places: `estimators.py` and
    `run_policy.yaml`.

    Let them drift and the cost depends on which path did the estimating, and the
    budget gate quietly uses a different threshold.
    """
    from simbot.estimators import EFFICIENCY_BY_K, THROUGHPUT_PARTICLE_STEPS_PER_S
    assert THROUGHPUT_PARTICLE_STEPS_PER_S == pytest.approx(
        policy.get("hardware.throughput_particle_steps_per_s"), rel=1e-12)
    yaml_eff = {int(k): float(v)
                for k, v in policy.get("hardware.efficiency_by_k").items()}
    assert EFFICIENCY_BY_K == pytest.approx(yaml_eff)
