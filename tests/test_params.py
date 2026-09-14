"""Rule 10 — the parameter manifest.

*"Every time you formulate a problem, lay out all the numbers you are using ...
before you start the run, and ask for permission."*

Structure: does it fire · is it wired · does turning it off change the verdict ·
is its scope narrow on purpose.
"""
from __future__ import annotations

from pathlib import Path

import pytest

from bdbot import params as PRM
from bdbot import runcard as RC

ROOT = Path(__file__).resolve().parent.parent


def water_probe(approved=None) -> PRM.Manifest:
    """The manifest the `bead-water-3d` run should have presented first.

    Every value is the one actually used, and every provenance is real.
    """
    m = PRM.Manifest("bead-water-3d", approved_by=approved)
    m.add("geometry", "d", 1.0, "um", "choice: probe size for this validation", 3)
    m.add("geometry", "radius", 0.5, "um",
          "d/2 -- GSER's `a` is the RADIUS (distillation trap 1)", 3)
    m.add("geometry", "L", 0.0, "1",
          "no box: force-free system has no interactions, so no PBC", 0,
          note="also removes the unwrap trap sim.unwrap guards")
    m.add("geometry", "dim", 3, "1",
          "GSER presumes a 3D spherical probe; 2D is refused", 0)
    m.add("geometry", "N", 120, "1", "choice: independent probes for statistics", 3)
    m.add("energy", "T", 300.0, "K", "tier 0, inherited from trap-2d-5um's sketch", 1)
    m.add("energy", "kT", 1.0, "1", "reduced convention: kT = 1", 0)
    m.add("medium", "eta", 0.851e-3, "Pa*s",
          "Welty anchor -- the value the 8 archived cases ran", 0,
          note="IAPWS gives 0.85566; the 0.548 % gap is deliberately unreconciled")
    m.not_applicable("interaction",
                     "free probe: no trap, no pair potential, no bonds. refbd "
                     "raises if any force is passed")
    m.add("numerics", "dt", 1.5e-4, "tau_B",
          "displacement gate: sqrt(2*dim*D*dt)/sigma <= 0.03 at dim=3", 0)
    m.add("numerics", "n_eq", 0, "1", "no equilibration: r(0)=0 is the exact "
          "initial condition for free diffusion", 0)
    m.add("numerics", "n_prod", 60_000, "1", "choice: gives T_obs = 9 tau_B", 3)
    m.add("numerics", "sample_every", 1, "1",
          "every step -- sample_every=10 costs a full decade of GSER band", 0)
    m.add("numerics", "seed", 20260902, "1", "choice", 3)
    m.add("derived", "gamma", 1.0, "1", "reduced: gamma = 1", 0)
    m.add("derived", "D_t", 1.0, "1", "reduced: D = kT/gamma = 1", 0)
    m.add("derived", "tau_gov", 1.0, "tau_B",
          "tau_B governs: no trap and no bond, so tau_B is the only timescale", 0)
    m.add("derived", "tau_B", 1.936, "s", "d^2/D_t at d=1um, water 300K", 0)
    return m


# ── 1. does it fire ────────────────────────────────────────────────────────

def test_a_complete_approved_manifest_is_ready():
    m = water_probe(approved="Saksham 2026-09-02")
    assert m.blockers() == [], m.blockers()
    assert m.ready is True


def test_unapproved_blocks_and_no_code_sets_it():
    m = water_probe()
    assert m.approved_by is None
    assert m.ready is False
    assert any("approved_by is unset" in b and "rule 10" in b.lower()
               for b in m.blockers()), m.blockers()


@pytest.mark.parametrize("category", PRM.CATEGORIES)
def test_every_category_must_be_present_or_declared_na(category):
    m = water_probe(approved="x")
    m.params.pop(category, None)
    m.absent.pop(category, None)
    b = m.blockers()
    assert any(category in x for x in b), (category, b)


@pytest.mark.parametrize("cat,key", [(c, k) for c, ks in PRM.REQUIRED_KEYS.items()
                                     for k in ks])
def test_every_required_key_is_required(cat, key):
    m = water_probe(approved="x")
    if cat in m.absent:
        pytest.skip(f"{cat} is N/A in this manifest")
    m.params[cat].pop(key, None)
    b = m.blockers()
    # Either blocker is correct: removing the only required key of a category
    # empties it, and "empty and not declared N/A" is the stronger message.
    assert any(key in x for x in b) or any(f"{cat!r} is empty" in x for x in b), \
        (cat, key, b)


# ── 2. a number cannot sneak in without value / unit / provenance ─────────

def test_a_bare_number_is_rejected():
    with pytest.raises(ValueError, match="no unit"):
        PRM.Param("d", 1.0, "", "somewhere")
    with pytest.raises(ValueError, match="no provenance"):
        PRM.Param("d", 1.0, "um", "")
    with pytest.raises(ValueError, match="no value"):
        PRM.Param("d", None, "um", "sketch")
    with pytest.raises(ValueError, match="tier"):
        PRM.Param("d", 1.0, "um", "sketch", tier=7)


def test_dimensionless_must_say_so_explicitly():
    """`unit='1'` is how a dimensionless quantity declares itself. Blank must mean
    forgotten, or the two states are indistinguishable."""
    p = PRM.Param("k_star", 6.04e4, "1", "k d^2/kT", 0)
    assert p.unit == "1"
    assert "6.04e+04" in p.line() or "60400" in p.line()


# ── 3. "if applicable" means stated, not omitted ──────────────────────────

def test_not_applicable_requires_a_reason():
    m = PRM.Manifest("x")
    with pytest.raises(ValueError, match="requires a reason"):
        m.not_applicable("interaction", "")
    m.not_applicable("interaction", "free probe, no forces at all")
    assert "interaction" in m.absent


def test_a_category_cannot_be_both_na_and_populated():
    m = PRM.Manifest("x")
    m.not_applicable("interaction", "free probe")
    with pytest.raises(ValueError, match="not applicable"):
        m.add("interaction", "k", 1.0, "1", "somewhere")

    m2 = PRM.Manifest("y")
    m2.add("interaction", "k", 1.0, "1", "sketch")
    with pytest.raises(ValueError, match="already carries"):
        m2.not_applicable("interaction", "changed my mind")


def test_an_unknown_number_blocks_and_names_its_supplier():
    """Rule 3: BLOCKED is correct, inventing the value is the failure."""
    m = water_probe(approved="x")
    m.unknown("lambda_p", "the experimentalist -- which Boger fluid?")
    b = m.blockers()
    assert any("lambda_p is UNKNOWN" in x and "experimentalist" in x for x in b), b
    assert not m.ready
    with pytest.raises(ValueError, match="who or what"):
        m.unknown("eta", "")


# ── 4. is it wired ────────────────────────────────────────────────────────

def test_a_run_card_without_a_manifest_is_not_sealable():
    card = RC.RunCard(case="c", run_id="r", question="Does X hold?",
                      cost_line="200 s", cannot_decide=("nothing",),
                      predictions=(("f_c", 1.0, 5.0, "implementation_check"),),
                      analysis=(("P1", "f", "d", "implementation_check"),))
    assert card.manifest is None
    assert any("rule 10" in b for b in card.blockers()), card.blockers()
    assert "NOT SEALABLE" in card.render()


def test_a_run_card_with_an_unapproved_manifest_is_not_sealable():
    card = RC.RunCard(case="bead-water-3d", run_id="r", question="Does X hold?",
                      cost_line="200 s", cannot_decide=("nothing",),
                      predictions=(("f_c", 1.0, 5.0, "implementation_check"),),
                      analysis=(("P1", "f", "d", "implementation_check"),),
                      manifest=water_probe())
    assert any("approved_by is unset" in b for b in card.blockers())


def test_a_complete_card_shows_the_numbers_before_the_design():
    """★ Rule 10's ordering: the numbers are what the reader is approving, so they
    come before the design and the cost."""
    card = RC.RunCard(case="bead-water-3d", run_id="r", question="Does X hold?",
                      cost_line="200 s vs 600 s", cannot_decide=("nothing",),
                      predictions=(("eta_star", 0.1061, 1.0, "implementation_check"),),
                      analysis=(("P1", "GSER", "is eta* recovered",
                                 "implementation_check"),),
                      manifest=water_probe(approved="Saksham 2026-09-02"))
    assert card.blockers() == [], card.blockers()
    text = card.render()
    assert "PARAMETER MANIFEST" in text
    assert text.index("PARAMETER MANIFEST") < text.index("COST")
    for name in ("radius", "kT", "dt", "n_prod", "seed", "tau_gov"):
        assert name in text, name
    assert "NOT APPLICABLE" in text          # the interaction category


def test_execute_checks_the_manifest_from_inside_itself():
    """Same lesson as the seal: a gate on an optional path is not a gate."""
    src = (ROOT / "bdbot" / "run.py").read_text()
    body = src.split("def execute(")[1].split("\ndef ")[0]
    assert "params.json" in body
    assert "require_approval" in body
    assert "diff_against_spec" in body


# ── 5. ★ approving one set of numbers and running another ─────────────────

def test_a_manifest_that_disagrees_with_the_spec_is_caught():
    """The failure that makes approval meaningful. Approving `dt = 1.5e-4` and
    then running `dt = 1.5e-3` must not pass silently."""
    m = water_probe(approved="Saksham")
    assert PRM.diff_against_spec(m, {"dt": 1.5e-4, "n_prod": 60_000}) == []

    drift = PRM.diff_against_spec(m, {"dt": 1.5e-3, "n_prod": 60_000})
    assert drift and "dt" in drift[0]
    assert "not the numbers about to run" in drift[0]

    drift2 = PRM.diff_against_spec(m, {"dt": 1.5e-4, "n_prod": 600_000})
    assert drift2 and "n_prod" in drift2[0]


def test_star_suffixed_spec_keys_are_matched():
    """Specs write `dt_star`; the manifest writes `dt`. They must still compare."""
    m = water_probe(approved="x")
    assert PRM.diff_against_spec(m, {"dt_star": 1.5e-4}) == []
    assert PRM.diff_against_spec(m, {"dt_star": 3.0e-4})


# ── 6. round trip and tier accounting ─────────────────────────────────────

def test_round_trip_through_json(tmp_path):
    m = water_probe(approved="Saksham 2026-09-02")
    p = m.write(tmp_path)
    assert p.name == "params.json"
    back = PRM.load(tmp_path)
    assert back.approved_by == "Saksham 2026-09-02"
    assert back.ready
    assert back.absent["interaction"].startswith("free probe")
    assert back.get("geometry", "radius").value == 0.5
    assert back.tier_counts() == m.tier_counts()


def test_tier_three_values_are_counted_and_flagged():
    """Tier 3 is an arbitrary assumption. The count is the honest measure of how
    much of a run is choice rather than measurement."""
    m = water_probe(approved="x")
    tc = m.tier_counts()
    assert tc.get(3, 0) >= 4          # d, radius, N, n_prod, seed are choices
    assert "arbitrary" in m.render()


def test_rule_10_is_in_the_project_contract():
    """A rule that lives only in a module docstring is the failure mode this rule
    was written to avoid. It has to be in CLAUDE.md."""
    md = (ROOT / "CLAUDE.md").read_text()
    assert "10 ·" in md
    assert "params.py" in md
    for word in ("radius", "box size", "trap stiffness", "provenance"):
        assert word in md, word


# ── 6. presence is not validity — the `derived` recomputation ──────────────
#
# ★ `bdbot/params.py`'s docstring said, from the day it was written, that
#   `derived` "is checked against a recomputation, because `physical.verify()`
#   already catches a hand-edited derived value and this is the same invariant
#   one layer up." No such check existed for twelve days. `blockers()` checked
#   that each `REQUIRED_KEYS` name was PRESENT and never that its value was
#   REAL, so this cleared the gate:
#
#       m.add("derived", "gamma",   0.0, "1", "6 pi eta a -- in the ledger")
#       m.add("derived", "D_t",     0.0, "1", "kT/gamma -- in the ledger")
#       m.add("derived", "tau_gov", 0.0, "1", "tau_B -- in the ledger")
#
#   Three physically impossible zeros -- no friction, no diffusion, no timescale
#   -- in the first manifest anyone built, written for the run that was meant to
#   demonstrate rule 10. The docstring's own list of practices-that-were-never-
#   enforced ("three for three") had a fourth entry in the same file.

#: The truth, at d = 5 um, T = 300 K, eta = 0.851 mPa*s -- `soft-r3`'s numbers.
#: Cross-checked against the closed forms by hand, independently of
#: `materials.sphere_bulk`, and they agree to 15 significant figures:
#:   gamma = 3 pi eta d  = 3*pi*0.851e-3*5e-6      = 4.0102430223073706e-08
#:   D_t   = kT/gamma    = 1.380649e-23*300/gamma  = 1.0328418943590235e-13
#:   tau_B = d^2/D_t                               = 242.05059977272592
#: The last-digit differences are float association order.
TRUE_GAMMA = 4.01024302230737e-08       # kg/s
TRUE_D_T = 1.0328418943590236e-13       # m^2/s
TRUE_TAU_B = 242.05059977272577         # s


def si_probe(approved="a human", **derived) -> PRM.Manifest:
    """A manifest whose inputs are SI, so `derived` is checkable by recomputation.
    `soft-r3`'s numbers: d = 5 um, T = 300 K, water at 0.851 mPa*s."""
    m = PRM.Manifest("si-probe", approved_by=approved)
    m.add("geometry", "d", 5.0, "um", "sketch", 0)
    m.add("geometry", "L", 29.96, "um", "a_mean*sqrt(N)", 1)
    m.add("geometry", "dim", 2, "1", "2D", 0)
    m.add("geometry", "N", 400, "1", "the spec", 1)
    m.add("energy", "T", 300.0, "K", "sketch", 0)
    m.add("energy", "kT", 1.0, "1", "reduced convention", 0)
    m.add("medium", "eta", 0.851, "mPa*s", "handbook, water at 300 K", 0)
    m.not_applicable("interaction", "checked elsewhere in this probe")
    m.add("numerics", "dt", 4.915e-05, "1", "dt/tau_B", 0)
    m.add("numerics", "n_prod", 2034400, "1", "T_obs/dt", 1)
    m.add("numerics", "seed", 20260914, "1", "fixed", 0)
    m.add("derived", "gamma", derived.get("gamma", TRUE_GAMMA),
          derived.get("gamma_u", "kg/s"), "6 pi eta a", 1)
    m.add("derived", "D_t", derived.get("D_t", TRUE_D_T),
          derived.get("D_t_u", "m^2/s"), "kT/gamma", 1)
    m.add("derived", "tau_gov", derived.get("tau_gov", TRUE_TAU_B),
          derived.get("tau_gov_u", "s"),
          "tau_B governs here: no trap and no bond", 1)
    m.add("derived", "tau_B", derived.get("tau_B", TRUE_TAU_B),
          derived.get("tau_B_u", "s"), "d^2/D_t", 1)
    return m


def test_the_si_probe_is_correct_to_begin_with():
    """Guard on the guard: every break below starts from this manifest, so a
    'caught' result proves nothing unless this one is clean."""
    assert si_probe().blockers() == []


def test_the_measured_defect_three_placeholder_zeros():
    """★ The exact manifest that cleared the gate on 2026-09-14."""
    m = si_probe(gamma=0.0, gamma_u="1", D_t=0.0, D_t_u="1",
                 tau_gov=0.0, tau_gov_u="1", tau_B=0.0, tau_B_u="1")
    b = m.blockers()
    for name in ("gamma", "D_t", "tau_gov", "tau_B"):
        assert any(f"derived.{name}" in x for x in b), (name, b)
    #  gamma / D_t / tau_B are refuted by the reduced-convention rule (they are
    #  1 by construction there); tau_gov by positivity, because a governing
    #  timescale is a case-dependent choice and has no single right value.
    assert any("by construction" in x for x in b), b
    assert any("cannot be zero or negative" in x for x in b), b


@pytest.mark.parametrize("secs", [0.0, -1.0])
def test_a_governing_timescale_must_be_positive(secs):
    m = si_probe(tau_gov=secs, tau_gov_u="s")
    b = [x for x in m.blockers() if "tau_gov" in x]
    assert b and "cannot be zero or negative" in b[0], b


def test_tau_gov_is_not_required_to_equal_tau_B():
    """★★ The fix that this check needed within the hour of being written.

    `tau_gov` is "which timescale GOVERNS" (CLAUDE.md rule 10), and it is a
    case-dependent choice. `trap-2d-5um`'s real card declares
    `tau_gov = 4.010e-3 s` with the provenance "tau_k = gamma/k governs, NOT
    tau_B = 242 s" -- right by five orders of magnitude. The first version of
    `check_derived` compared `tau_gov` against tau_B and refused that manifest
    at -100 %. A gate that refuses a correct answer is worse than no gate; this
    repository's own pre-run gate once rejected 80 of 83 specs with zero real
    failures among them.
    """
    m = si_probe(tau_gov=4.010e-3, tau_gov_u="s")      # a trap, not diffusion
    assert [x for x in m.blockers() if "tau_gov" in x] == [], m.blockers()
    #  ...and tau_B, which IS derivable, is still checked in the same manifest
    m2 = si_probe(tau_gov=4.010e-3, tau_gov_u="s", tau_B=TRUE_TAU_B * 1.05)
    assert [x for x in m2.blockers() if "derived.tau_B" in x], m2.blockers()


def test_the_real_trap_card_manifest_is_accepted():
    """The regression, on the actual numbers rather than a mutation of mine.
    `trap-2d-5um`'s manifest is written to 4-5 significant figures throughout,
    which is why DERIVED_RTOL is 1e-3 and not 1e-6."""
    m = PRM.Manifest("trap-2d-5um", approved_by="a human")
    m.add("geometry", "d", 5.0, "um", "sketch", 0)
    m.add("geometry", "L", 160.0, "um", "n_side * d", 0)
    m.add("geometry", "dim", 2, "1", "sketch: a 2D trap", 0)
    m.add("geometry", "N", 1000, "1", "choice", 3)
    m.add("energy", "T", 300.0, "K", "sketch", 0)
    m.add("energy", "kT", 4.142e-21, "J", "k_B T at 300 K", 0)
    m.add("medium", "eta", 0.851e-3, "Pa*s", "Welty, water at 300 K", 0)
    m.add("interaction", "k_t", 10.0, "pN/um", "sketch", 0)
    m.add("numerics", "dt", 8.02e-6, "s", "0.1 % EM bias", 0)
    m.add("numerics", "n_prod", 1_000_000, "1", "T_obs / dt", 0)
    m.add("numerics", "seed", 20260803, "1", "choice", 3)
    m.add("derived", "gamma", 4.0102e-8, "kg/s", "3 pi eta d", 0)
    m.add("derived", "D_t", 0.1033, "um**2/s", "kT/gamma", 0)
    m.add("derived", "tau_gov", 4.010e-3, "s",
          "tau_k = gamma/k governs, NOT tau_B = 242 s", 0)
    assert m.blockers() == [], m.blockers()


@pytest.mark.parametrize("key,unit,truth", [
    ("gamma", "kg/s", TRUE_GAMMA),
    ("D_t", "m^2/s", TRUE_D_T),
    ("tau_B", "s", TRUE_TAU_B),
])
def test_a_typo_in_a_recomputable_derived_value_is_caught(key, unit, truth):
    """1 % -- small enough to miss by eye, ten times the tolerance. The tight
    comparison lives one layer down in `physical.verify()`, where both sides
    are full precision; here the tolerance is set by how the value was written
    (four to five significant figures across the real manifests)."""
    m = si_probe(**{key: truth * 1.01, f"{key}_u": unit})
    b = [x for x in m.blockers() if f"derived.{key}" in x]
    assert b and "+1" in b[0], b


def test_the_recomputation_is_not_vacuous():
    """...and it accepts the truth. Without this the test above would pass for a
    check that rejects every value."""
    assert [x for x in si_probe().blockers() if "derived." in x] == []


def test_a_derived_value_in_the_wrong_dimension_is_refused_not_converted():
    m = si_probe(gamma=TRUE_GAMMA, gamma_u="m")
    b = [x for x in m.blockers() if "derived.gamma" in x]
    assert b and "not convertible" in b[0], b


@pytest.mark.parametrize("missing", ["d", "T", "eta"])
def test_a_missing_input_is_reported_rather_than_skipped(missing):
    """An unrunnable check that returns [] is indistinguishable from a check
    that ran and passed -- this repository's most-repeated defect
    (docs/05 section 2). It must say it could not run."""
    cat = {"d": "geometry", "T": "energy", "eta": "medium"}[missing]
    m = si_probe()
    del m.params[cat][missing]
    b = [x for x in m.blockers() if "cannot be recomputed" in x]
    assert b and f"{cat}.{missing}" in b[0], b


def test_the_reduced_convention_is_checked_not_exempted():
    """`gamma = 1` with a reduced unit is correct BY CONSTRUCTION, and 0 is not.
    `bead-water-3d`'s real manifest is of this kind, so exempting reduced units
    would have left the original defect open on the very manifest that shipped."""
    assert water_probe("a human").blockers() == []
    m = water_probe("a human")
    m.params["derived"]["gamma"] = PRM.Param("gamma", 0.0, "1", "reduced", 0)
    b = m.blockers()
    assert b and "by construction" in b[0], b


def test_declaring_derived_not_applicable_skips_the_check_explicitly():
    """The escape hatch is `not_applicable`, which raises on an empty reason --
    the same shape as `scales.declare_absent`. Silence is not an escape hatch."""
    m = si_probe()
    m.params.pop("derived")
    m.not_applicable("derived", "a probe with no medium: nothing to derive")
    assert m.blockers() == []
    with pytest.raises(ValueError):
        PRM.Manifest("x").not_applicable("derived", "   ")


def test_the_docstring_no_longer_claims_an_unimplemented_check():
    """The paragraph promising the recomputation stood for twelve days with
    nothing behind it. If the check is ever removed, the docstring must not
    quietly go back to claiming it."""
    src = Path(PRM.__file__).read_text()
    assert "Four for four" in src, \
        "the correction record was removed from bdbot/params.py"
    assert "def check_derived" in src, \
        "check_derived is gone but the docstring still promises a recomputation"
