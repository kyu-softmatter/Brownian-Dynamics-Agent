"""S3 — proposing a cutoff radius. No HOOMD needed (pure numerics), all fast.

The central check: **does `r_cut` come from a stated error tolerance rather than
from convention**, and **is it extended for dilute or long-range systems and
limited by the budget for dense ones.**
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
# WCA — fixed by definition
# =============================================================================
def test_wca_cutoff_is_the_lj_minimum_exactly():
    p = wca_cutoff()
    assert p.r_cut_over_sigma == pytest.approx(2 ** (1 / 6), rel=1e-15)
    assert p.r_cut_over_sigma == pytest.approx(1.1224620483, rel=1e-9)
    assert p.exact is True
    assert p.criterion == "potential_definition"
    # at the minimum the force and potential are exactly 0 — no truncation artefact
    assert p.beta_U_at_cut == 0.0 and p.beta_F_sigma_at_cut == 0.0


def test_wca_force_vanishes_at_its_cutoff():
    """★ Confirm directly that the LJ force is exactly 0 at r_cut = 2^{1/6} sigma."""
    r = WCA_RCUT_OVER_SIGMA_LJ
    f = 24.0 * (2.0 * r ** -13 - r ** -7)          # -dU/dr, eps=sigma=1
    assert f == pytest.approx(0.0, abs=1e-14)


def test_wca_hoomd_recipe_records_the_shift_requirement():
    """Does the recipe state the mode="shift" the lab code left out?"""
    p = wca_cutoff()
    assert 'mode="shift"' in p.hoomd_note
    assert "2**(1/6)" in p.hoomd_note or "2^{1/6}" in p.rationale


def test_propose_cutoff_short_circuits_for_exact_potentials():
    """WCA's value does not change with phi or with the budget."""
    for phi, n in [(0.0, 1), (0.4, 10000)]:
        p = propose_cutoff(wca_cutoff(), phi=phi, n_particles=n)
        assert p.r_cut_over_sigma == pytest.approx(WCA_RCUT_OVER_SIGMA_LJ, rel=1e-15)
        assert p.criterion == "potential_definition"


# =============================================================================
# Barker-Henderson — an independent comparison against the knowledge base
# =============================================================================
@pytest.mark.benchmark
def test_barker_henderson_reproduces_knowledge_base_value():
    """★ knowledge/wiki/findings/wca-reproduces-carnahan-starling.md records
    `d_eff = 1.017 sigma` at `beta*eps = 1`.

    This test reproduces that independently (numerical integration).
    """
    d = barker_henderson_diameter(beta_epsilon=1.0)
    assert d == pytest.approx(1.017, rel=3e-3), \
        f"d_BH = {d:.5f}, knowledge value 1.017"


def test_barker_henderson_increases_with_stiffness_toward_r_min():
    """As eps grows the effective hard sphere stiffens and approaches r_min (without
    exceeding it)."""
    ds = [barker_henderson_diameter(be) for be in (0.5, 1.0, 2.0, 10.0, 500.0)]
    assert all(a < b for a, b in zip(ds, ds[1:])), \
        "monotonically increasing in beta*eps"
    assert ds[-1] < WCA_RCUT_OVER_SIGMA_LJ, "cannot exceed r_min"
    assert ds[-1] == pytest.approx(WCA_RCUT_OVER_SIGMA_LJ, rel=1e-2), \
        "approaches it in the strong limit"


def test_barker_henderson_rejects_nonpositive_epsilon():
    for bad in (0.0, -1.0):
        with pytest.raises(ValueError):
            barker_henderson_diameter(bad)


# =============================================================================
# The tolerance-based solution — presets and the dominant constraint
# =============================================================================
def test_tolerance_presets_are_ordered_from_loose_to_strict():
    order = ["convention", "structure", "thermodynamics", "precision"]
    us = [TOLERANCE_PRESETS[k][0] for k in order]
    assert all(a > b for a, b in zip(us, us[1:])), \
        "the tolerances have to decrease monotonically"


def test_convention_preset_reproduces_the_literature_2p5_sigma():
    """★ Does the `convention` preset reproduce the literature's 2.5 sigma?

    Failing to would mean the preset's tolerance is a value unrelated to the
    convention.
    """
    p = lj_cutoff(1.0, preset="convention")
    assert p.r_cut_over_sigma == pytest.approx(CONVENTIONAL_LJ_RCUT, rel=0.05)


def test_stricter_preset_gives_larger_cutoff():
    rs = [lj_cutoff(1.0, preset=k).r_cut_over_sigma
          for k in ("convention", "structure", "thermodynamics", "precision")]
    assert all(a < b for a, b in zip(rs, rs[1:]))


def test_lj_cutoff_reports_residual_at_the_convention():
    """The basis has to record what is left at the conventional value — the
    material for a judgment."""
    p = lj_cutoff(1.0, preset="thermodynamics")
    assert "2.5 sigma" in p.rationale
    assert p.r_cut_over_sigma == pytest.approx(3.984, rel=1e-3)
    assert p.beta_U_at_cut == pytest.approx(1e-3, rel=1e-6)


def test_cutoff_records_which_constraint_dominated():
    p = lj_cutoff(1.0)
    assert p.criterion in ("potential_tolerance", "force_tolerance")
    assert "dominates" in p.rationale


def test_force_tolerance_can_dominate_when_made_strict():
    """★ Does the case where the force constraint dominates actually occur?

    Test only one branch and you will not notice when the other one dies.
    """
    p = lj_cutoff(1.0, beta_u_tol=1.0, beta_f_sigma_tol=1e-6)
    assert p.criterion == "force_tolerance"
    assert p.beta_F_sigma_at_cut == pytest.approx(1e-6, rel=1e-5)


def test_unknown_preset_raises():
    with pytest.raises(KeyError, match="unknown preset"):
        lj_cutoff(1.0, preset="vibes")


def test_long_range_that_never_meets_tolerance_is_flagged():
    """Still outside tolerance at r_max has to warn (never truncate quietly)."""
    p = cutoff_from_tolerance(
        potential="coulomb-like", beta_u=lambda r: 100.0 / r,
        beta_f_sigma=lambda r: 100.0 / r ** 2, r_start=1.0, r_max=10.0)
    assert p.criterion == "r_max_reached"
    assert p.warnings and "long-range" in p.warnings[0]


# =============================================================================
# Yukawa — the screening sets the range
# =============================================================================
def test_yukawa_cutoff_grows_as_screening_weakens():
    """★ The interfacial-colloid card's trap: as kappa*sigma -> 0 the range
    explodes."""
    rs = [yukawa_cutoff(10.0, ks).r_cut_over_sigma for ks in (10, 5, 2, 1, 0.3, 0.1)]
    assert all(a < b for a, b in zip(rs, rs[1:])), \
        "kappa down → r_cut up"
    assert rs[-1] > 50, f"at kappa*sigma=0.1, r_cut = {rs[-1]:.1f} sigma"


def test_yukawa_cutoff_scales_like_inverse_kappa():
    """r_cut ~ ln(beta eps / tol)/kappa, so halving kappa roughly doubles it."""
    r1 = yukawa_cutoff(10.0, 1.0).r_cut_over_sigma
    r2 = yukawa_cutoff(10.0, 0.5).r_cut_over_sigma
    assert 1.7 < r2 / r1 < 2.4


def test_minimum_image_violation_is_detected():
    """★ Reproduces the real interfacial-colloid case: the box cannot contain the
    interaction range."""
    p = yukawa_cutoff(beta_epsilon=1e5, kappa_sigma=0.1)
    ok, msg = p.box_ok(130.0)
    assert not ok
    assert "minimum image violation" in msg
    assert f"{p.min_box_over_sigma:.4g}" in msg
    # a large enough box passes
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
# Neighbour count — the entire cost
# =============================================================================
@pytest.mark.parametrize("dim,coef", [(3, 8.0), (2, 4.0)])
def test_neighbor_count_formula(dim, coef):
    phi, r = 0.3, 2.5
    expected = coef * phi * r ** dim
    assert neighbors_per_particle(r, phi, dim) == pytest.approx(expected, rel=1e-14)


def test_neighbor_count_is_zero_for_single_particle_limit():
    assert neighbors_per_particle(10.0, 0.0, 3) == 0.0


def test_conventional_lj_neighbor_count_is_in_the_expected_range():
    """phi=0.3, r_cut=2.5 sigma should be ~40, the typical 3D LJ-liquid value."""
    nb = neighbors_per_particle(CONVENTIONAL_LJ_RCUT, 0.3, 3)
    assert 30 < nb < 50, f"{nb:.1f}"


def test_lambda_baseline_cost_is_unity_at_measurement_conditions():
    """★ At the condition the throughput constant was measured under, the cost
    multiple has to be 1."""
    c = pair_cost_vs_lambda_baseline(LAMBDA_BASELINE_RCUT, LAMBDA_BASELINE_PHI, 3)
    assert c == pytest.approx(1.0, rel=1e-12)


def test_lambda_baseline_was_measured_in_a_cheap_regime():
    """★ Lambda was measured on WCA (3.4 neighbours) — conventional LJ is 10x more
    expensive than that.

    Forget this and the wall-time estimate is optimistically wrong.
    """
    nb_base = neighbors_per_particle(LAMBDA_BASELINE_RCUT, LAMBDA_BASELINE_PHI, 3)
    assert nb_base == pytest.approx(3.4, rel=0.05)
    cost_conv = pair_cost_vs_lambda_baseline(CONVENTIONAL_LJ_RCUT, 0.3, 3)
    assert cost_conv > 9.0, f"conventional LJ cost {cost_conv:.1f}x"


# =============================================================================
# ★ The final proposal — the user's rule (convention first, extend for dilute and
# long-range)
# =============================================================================
def test_single_particle_makes_cutoff_irrelevant():
    p = propose_cutoff(lj_cutoff(1.0), phi=0.0, n_particles=1, dim=2)
    assert p.criterion == "single_particle_free"
    assert any("1 particle" in w for w in p.warnings)


def test_dilute_system_gets_the_full_tolerance_cutoff():
    """★ The user's rule: when dilute it may go longer than the convention — with
    no neighbours it is free."""
    tol = lj_cutoff(1.0, preset="thermodynamics")
    p = propose_cutoff(tol, phi=0.01, n_particles=100, dim=3)
    assert p.r_cut_over_sigma == pytest.approx(tol.r_cut_over_sigma, rel=1e-12)
    assert p.criterion == "dilute_or_longrange_affordable"
    assert neighbors_per_particle(p.r_cut_over_sigma, 0.01, 3) < 10


def test_dense_system_is_limited_by_neighbor_budget():
    """In a dense system the budget binds, and **the residual error has to be
    recorded as a limit.**"""
    tol = lj_cutoff(1.0, preset="thermodynamics")
    p = propose_cutoff(tol, phi=0.40, n_particles=4000, dim=3)
    assert p.criterion == "neighbor_budget_limited"
    assert p.r_cut_over_sigma < tol.r_cut_over_sigma
    nb = neighbors_per_particle(p.r_cut_over_sigma, 0.40, 3)
    assert nb == pytest.approx(DEFAULT_NEIGHBOR_BUDGET, rel=1e-9)
    assert any("known limit" in w for w in p.warnings)
    assert math.isnan(p.beta_U_at_cut), \
        "having reduced it, the original tolerance cannot be claimed"


def test_budget_limited_never_goes_below_convention():
    """However tight the budget, it never drops below the conventional value."""
    tol = lj_cutoff(1.0, preset="precision")
    p = propose_cutoff(tol, phi=0.6, n_particles=10000, dim=3, neighbor_budget=1.0)
    assert p.r_cut_over_sigma >= CONVENTIONAL_LJ_RCUT


def test_tolerance_below_convention_is_taken_as_free_win():
    """With strong screening the tolerance solution is below the convention — use
    it as is."""
    tol = yukawa_cutoff(beta_epsilon=10.0, kappa_sigma=5.0)
    assert tol.r_cut_over_sigma < CONVENTIONAL_LJ_RCUT
    p = propose_cutoff(tol, phi=0.3, n_particles=2000, dim=3)
    assert p.criterion == "tolerance_below_convention"
    assert p.r_cut_over_sigma == pytest.approx(tol.r_cut_over_sigma, rel=1e-12)


def test_long_range_dilute_is_affordable_but_may_break_the_box():
    """★ Dilute + long-range: the neighbour budget passes but minimum image may
    still bind."""
    tol = yukawa_cutoff(beta_epsilon=10.0, kappa_sigma=0.3)
    p = propose_cutoff(tol, phi=0.001, n_particles=2000, dim=3, box_over_sigma=40.0)
    assert p.criterion == "dilute_or_longrange_affordable"
    assert any("minimum image" in w for w in p.warnings)


def test_proposal_always_reports_cost_and_neighbor_count():
    """The basis must always carry the cost and the neighbour count — the material
    for a budget judgment."""
    p = propose_cutoff(lj_cutoff(1.0), phi=0.2, n_particles=1000, dim=3)
    assert "neighbours" in p.rationale
    assert "cost" in p.rationale


# =============================================================================
# The neighbour-list buffer
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
