"""S3 — the system-specification data model + validity checks.

The central regression: **does the code reproduce the hand-written derived values
of the first run.** The `derived:` block in
`runs/2026-07-28_trap-2d-5um_2dfb9d/03_spec.yaml` was transcribed by a person from
`simbot`'s output. Those values are the answer key, and a mismatch here means one
of the two is wrong — either way that has to be known.
"""
from __future__ import annotations

import math

import pytest

from simbot import spec as S
from simbot.spec import (Friction, Gate, Geometry, Medium, Numerics, PredictionItem,
                         Prediction, Q, Quantity, Species, SystemSpec, Timing,
                         ExternalField, PairInteraction)

EXAMPLE_SPEC = S.Path("examples/trap-2d-5um/spec.yaml")

# The hand-written 03_spec.yaml §derived — copied **exactly as the strings appear
# in the document**. Why keep them as strings: so the tolerance for each value can
# be taken from that value's own significant figures. `6.817e-6` has 4 figures, so
# its rounding half-width is 7.3e-5; `4.14195e-21` has 6, so 7.2e-7. One global
# constant would let the loose end excuse the strict value.
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
    """The rounding half-width (relative) of a value written as `written`.

    A value `v` written to `n` significant figures may have been truncated by ±0.5
    in its last digit. That is **the only theoretically justified tolerance** when
    comparing a document against a recomputation.
    """
    mant, _, exp = written.lower().partition("e")
    digits = len(mant.replace("-", "").replace(".", "").lstrip("0")) or 1
    v = abs(float(written))
    # size of the last significant figure = v's 10^(floor(log10 v) - (n-1))
    last_place = 10.0 ** (math.floor(math.log10(v)) - (digits - 1))
    return 0.5 * last_place / v


@pytest.fixture
def trap_spec() -> SystemSpec:
    if not EXAMPLE_SPEC.exists():
        pytest.skip("examples/trap-2d-5um/spec.yaml is absent")
    return SystemSpec.load(EXAMPLE_SPEC)


# =============================================================================
# Quantity — the provenance convention
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
    q = Quantity(value=1.0, provenance="assumed", basis="assumed to be water")
    assert any("confidence" in p for p in q.problems())


def test_derived_does_not_require_confidence():
    q = Quantity(value=1.0, provenance="derived", basis="sigma = 2a")
    assert q.problems() == []


def test_cheap_model_cannot_write_inference_field():
    """master_plan §12.2 — inference/assumed is Opus-only."""
    q = Quantity(value=2, provenance="inference", basis="there is no z axis",
                 confidence="medium", written_by="haiku")
    assert any("§12.2" in p for p in q.problems())


def test_cheap_model_may_write_observation_field():
    q = Quantity(value=1, provenance="observation", basis="one circle",
                 written_by="haiku")
    assert q.problems() == []


def test_si_property_refuses_non_numeric():
    with pytest.raises(TypeError):
        Q("periodic", provenance="rule", basis="x").si


def test_si_property_refuses_bool():
    """`True` is 1 in Python — let it into unit arithmetic and it is quietly wrong."""
    with pytest.raises(TypeError):
        Q(True, provenance="rule", basis="x").si


# =============================================================================
# gates
# =============================================================================
def test_gate_rejects_unknown_status():
    with pytest.raises(ValueError):
        Gate(status="probably")


def test_unknown_gate_name_is_a_problem(trap_spec):
    """A mistyped gate name is 'a check that never runs'."""
    trap_spec.gates["equiparition"] = Gate("required")     # a typo
    assert any("unregistered gate name" in p
               for p in S.validate(trap_spec).problems)


def test_gate_off_without_reason_is_a_problem(trap_spec):
    trap_spec.gates["equilibration_detection"] = Gate("off", "")
    assert any("with no reason" in p for p in S.validate(trap_spec).problems)


def test_required_gate_is_deferred_not_passed(trap_spec):
    """★ A declaration is not a result. Stamping `required` as pass is a pass
    nobody looked at."""
    rep = S.validate(trap_spec)
    names = {c.name: c.status for c in rep.checks}
    assert names["equipartition"] == "declared"
    assert names["em_bias_reproduced"] == "declared"
    assert "equipartition" in [c.name for c in rep.deferred()]


def test_computed_result_fills_declared_gate(trap_spec):
    """A gate the card declares and S3 can compute gets a real verdict."""
    rep = S.validate(trap_spec)
    c = next(c for c in rep.checks if c.name == "overdamped")
    assert c.status == "pass" and c.value is not None


def test_no_duplicate_check_rows(trap_spec):
    rep = S.validate(trap_spec)
    names = [c.name for c in rep.checks]
    assert len(names) == len(set(names))


# =============================================================================
# derived values — the first-run regression
# =============================================================================
@pytest.mark.benchmark
def test_derive_reproduces_hand_written_first_run(trap_spec):
    """Reproduces all 10 derived values from the hand-written document.

    The tolerance comes from **each value's own significant figures as written**
    (conftest rule 1 — never cut to fit the observed value). If the recomputation
    lands inside the rounding half-width, "the document and the code did the same
    calculation" is proved.
    """
    d = S.derive(trap_spec)
    for key, written in HAND_DERIVED.items():
        hand = float(written)
        tol = rounding_halfwidth(written)
        rel = abs(d[key] - hand) / abs(hand)
        assert rel <= tol, (f"{key}: recomputed {d[key]:.6g} vs document {written} "
                            f"(relative diff {rel:.2e} > rounding half-width "
                            f"{tol:.2e})")


def test_var_per_component_is_kT_over_k(trap_spec):
    """Equipartition `⟨x²⟩ = kT/k` — the quantity this system's first-class gate
    looks at."""
    d = S.derive(trap_spec)
    assert d["var_per_component_si"] == pytest.approx(d["kT_si"] / d["k_si"], rel=1e-15)
    assert d["msd_plateau_si"] == pytest.approx(
        2 * trap_spec.geometry.d * d["var_per_component_si"], rel=1e-15)


def test_rounding_halfwidth_matches_significant_figures():
    """Verifies the tolerance calculator itself — wrong, and the regression above
    means nothing."""
    assert rounding_halfwidth("6.817e-6") == pytest.approx(0.5e-3 / 6.817, rel=1e-9)
    assert rounding_halfwidth("4.14195e-21") == pytest.approx(0.5e-5 / 4.14195, rel=1e-9)
    assert rounding_halfwidth("19.735") == pytest.approx(0.5e-3 / 19.735, rel=1e-9)
    assert rounding_halfwidth("1.0e-5") == pytest.approx(0.5e-1 / 1.0, rel=1e-9)


def test_derive_gamma_uses_radius_not_diameter(trap_spec):
    """The most common mistake: passing the diameter to 6πηa. Do it and every
    timescale is wrong by 2x."""
    d = S.derive(trap_spec)
    a = trap_spec.primary.radius_si.si
    expected = 6.0 * math.pi * trap_spec.medium.eta_si.si * a
    assert d["gamma_si"] == pytest.approx(expected, rel=1e-12)
    assert d["gamma_si"] != pytest.approx(2 * expected, rel=1e-3)


def test_sigma_is_twice_the_radius(trap_spec):
    assert trap_spec.primary.sigma_si == pytest.approx(
        2 * trap_spec.primary.radius_si.si, rel=1e-15)


def test_reference_length_is_l_trap_for_trap_card(trap_spec):
    """★ A trap system's reference length is ℓ_trap, not σ (card §3)."""
    d = S.derive(trap_spec)
    assert S.reference_length_si(trap_spec, d) == pytest.approx(d["l_trap_si"])
    assert d["l_trap_si"] / d["sigma_si"] < 1e-2      # 0.2 % of its own diameter


def test_time_scale_separation_is_large(trap_spec):
    """Why this system rejects the τ_D convention — a measured factor of 2.41e5."""
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
# comparing derived values — catching a hand-edited number
# =============================================================================
HAND_DERIVED_F = {k: float(v) for k, v in HAND_DERIVED.items()}


def test_stored_derived_mismatch_is_caught(trap_spec):
    """★ The check that catches the 2026-07-28 class of error: kT wrong in the 4th
    digit."""
    bad = dict(HAND_DERIVED_F, kT_si=4.1420e-21 * 1.01)      # a value 1 % off
    rep = S.validate(trap_spec, stored_derived=bad)
    assert not rep.ok
    assert any("kT_si" in p for p in rep.problems)


def test_stored_derived_match_passes(trap_spec):
    rep = S.validate(trap_spec, stored_derived=HAND_DERIVED_F)
    assert rep.ok, rep.problems
    c = next(c for c in rep.checks if c.name == "derived_consistency")
    assert c.status == "pass"
    # did a comparison actually happen — passing with 0 comparisons is an empty test
    assert f"{len(HAND_DERIVED_F)} compared" in c.detail


def test_zero_comparisons_is_not_reported_as_pass(trap_spec):
    """★ Reporting 0 comparisons as pass makes 'it was checked' a false statement."""
    rep = S.validate(trap_spec, stored_derived={"unknown_key": 1.0,
                                                "note": "unbounded-medium Stokes"})
    c = next(c for c in rep.checks if c.name == "derived_consistency")
    assert c.status == "na"
    assert c.value == 0.0 and "uncomparable" in c.detail


def test_stored_derived_ignores_non_numeric(trap_spec):
    rep = S.validate(trap_spec, stored_derived={"note": "unbounded-medium Stokes"})
    assert rep.ok, rep.problems


# =============================================================================
# validity checks
# =============================================================================
def test_example_spec_has_no_problems(trap_spec):
    rep = S.validate(trap_spec)
    assert rep.problems == []
    assert rep.failed() == []


def test_overdamped_fails_for_heavy_particle(trap_spec):
    trap_spec.primary.density_si = Q(1e9, "kg/m^3", "assumed", "unrealistically heavy",
                                     confidence="low")
    c = next(c for c in S.validate(trap_spec).checks if c.name == "overdamped")
    assert c.status == "fail" and "Langevin" in c.detail


def test_box_too_small_fails(trap_spec):
    trap_spec.geometry.box_over_ref = Q(3.0, "l_trap", "rule", "deliberately small")
    c = next(c for c in S.validate(trap_spec).checks
             if c.name == "box_much_larger_than_l_trap")
    assert c.status == "fail"


def test_packing_fraction_is_off_without_pair_interactions(trap_spec):
    """★ Gating on φ for non-interacting replicas in the same trap gives φ=4741 and
    a can-never-pass verdict — a problem that does not exist. Show the value, do
    not judge it."""
    c = next(c for c in S.validate(trap_spec).checks if c.name == "packing_fraction")
    assert c.status == "off"
    assert c.value is not None and c.value > 1.0        # the value is still reported


def test_r_cut_gate_catches_minimum_image_violation(trap_spec):
    d = S.derive(trap_spec)
    box_min = 200.0 * d["l_trap_si"]
    trap_spec.pair = [PairInteraction(
        "probe", "probe", "wca",
        r_cut_si=Q(box_min, "m", "rule", "exceeds L/2 — minimum-image violation"))]
    trap_spec.gates.pop("r_cut_le_half_box")
    c = next(c for c in S.validate(trap_spec).checks if c.name == "r_cut_le_half_box")
    assert c.status == "fail"


def test_pair_without_r_cut_is_a_problem(trap_spec):
    trap_spec.pair = [PairInteraction("probe", "probe", "wca")]
    assert any("r_cut" in p for p in S.validate(trap_spec).problems)


def test_missing_box_is_a_problem(trap_spec):
    trap_spec.geometry.box_over_ref = None
    trap_spec.geometry.box_si = None
    assert any("box" in p for p in S.validate(trap_spec).problems)


def test_explicit_box_si_is_used(trap_spec):
    trap_spec.geometry.box_over_ref = None
    trap_spec.geometry.box_si = Q([4.07e-6, 4.07e-6, 0.0], "m", "rule",
                                  "an explicit box")
    assert trap_spec.box_lengths_si()[0] == pytest.approx(4.07e-6)


# =============================================================================
# the YAML round-trip — catching a silent unit loss
# =============================================================================
def test_yaml_roundtrip_is_exact(trap_spec):
    back = SystemSpec.from_yaml(trap_spec.to_yaml())
    assert S.to_dict(back) == S.to_dict(trap_spec)


def test_yaml_roundtrip_preserves_derived_bit_for_bit(trap_spec):
    """After a round-trip the derived values must be **exactly** the same
    (master_plan §S4 gate: < 1e-12)."""
    back = SystemSpec.from_yaml(trap_spec.to_yaml())
    d1, d2 = S.derive(trap_spec), S.derive(back)
    assert set(d1) == set(d2)
    for k in d1:
        assert abs(d1[k] - d2[k]) <= 1e-12 * abs(d1[k]), k


def test_hash_is_stable_across_roundtrip(trap_spec):
    assert SystemSpec.from_yaml(trap_spec.to_yaml()).hash() == trap_spec.hash()


def test_bonds_and_angles_survive_yaml_roundtrip(trap_spec):
    """A bond's and an angle's `k_si` have different dimensions (N/m vs J/rad²) —
    mixed in the round-trip, λ_max comes out wrong."""
    from simbot.spec import AngleInteraction, BondInteraction
    trap_spec.bonds = [BondInteraction(params={
        "k_si": Q(1.0e-3, "N/m", "rule", "bond spring"),
        "r0_si": Q(1.0e-5, "m", "rule", "bond length = σ")})]
    trap_spec.angles = [AngleInteraction(params={
        "k_si": Q(4.14e-17, "J/rad^2", "rule", "angle spring"),
        "t0_rad": Q(math.pi, "rad", "rule", "a straight chain")})]
    back = SystemSpec.from_yaml(trap_spec.to_yaml())
    assert S.to_dict(back) == S.to_dict(trap_spec)
    assert back.bond_stiffness_si() == trap_spec.bond_stiffness_si()
    assert back.angle_stiffness_si() == trap_spec.angle_stiffness_si()
    assert S.derive(back)["lambda_max_si"] == S.derive(trap_spec)["lambda_max_si"]


def test_empty_bonds_do_not_change_the_spec_hash(trap_spec):
    """★ Adding a field must leave the hash of a bond-free system unchanged —
    otherwise every existing run's `run_id` and cache key shifts."""
    text = trap_spec.to_yaml()
    assert "bonds:" not in text and "angles:" not in text


def test_hash_changes_with_physics(trap_spec):
    before = trap_spec.hash()
    trap_spec.medium.T_si = Q(310.0, "K", "from_drawing", "a different temperature")
    assert trap_spec.hash() != before


def test_yaml_always_writes_provenance(trap_spec):
    """Omit `assumed` because it is the default and 'assumed' becomes
    indistinguishable from 'forgot to write it'."""
    text = trap_spec.to_yaml()
    assert text.count("provenance:") == len(S._iter_quantities(trap_spec))
    assert "provenance: assumed" in text


def test_bare_value_loads_as_assumed_with_no_basis():
    """A bare value with no provenance must flow through as a convention violation
    (never a silent pass)."""
    q = S.from_dict(Quantity, 300.0)
    assert q.provenance == "assumed"
    assert any("basis" in p for p in q.problems())


# =============================================================================
# predictions (S2)
# =============================================================================
def test_empty_prediction_is_a_problem():
    assert any("0 quantitative predictions" in p
               for p in Prediction(items=[]).problems())


def test_prediction_item_requires_tolerance_and_basis():
    it = PredictionItem(quantity="D*", value=1.0, tolerance="", basis="")
    probs = it.problems()
    assert any("tolerance" in p for p in probs) and any("basis" in p for p in probs)


def test_valid_prediction_has_no_problems():
    p = Prediction(items=[PredictionItem(
        quantity="var_x_star", value=1.00251, tolerance="±1%",
        basis="EM bias 1/(1-dt*/2)", discriminates="integrator scheme",
        competing_value=1.0)])
    assert p.problems() == []
