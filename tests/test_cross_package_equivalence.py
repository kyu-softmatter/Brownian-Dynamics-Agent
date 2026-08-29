"""`bdbot` and `simbot` must not drift apart where they overlap.

The 2026-08-29 de-duplication pass merged what could be merged and left three
things deliberately unmerged, because merging them would change numbers that
existing runs are hashed against. This file is what keeps both halves honest:

  MERGED, so assert identity        -- the same function object reached from both
                                       packages. If someone re-forks it, this fails
  UNMERGED, so assert the gap       -- pin the measured difference. If someone
                                       "fixes" it, this fails and they have to read
                                       why it was left alone

The second kind is the unusual one and the reason this file exists. A duplicate
that agrees today is not safe; it is one well-meaning edit away from disagreeing,
and in this domain a 0.5 % disagreement in `eta` is invisible in the output.

Named by: `bdbot/materials.py`, `bdbot/runid.py`, `simbot/io.py`.
"""
from __future__ import annotations

import json

import pytest

import bdbot.dt as BDT
import bdbot.health as BH
import simbot.nondim as SN
import simbot.units as SU
from bdbot import constants as BC
from bdbot import materials as BM
from bdbot.runid import spec_hash
from bdbot.units import Q
from simbot.estimators import (dt_star_for_trap_bias,
                               euler_maruyama_trap_variance_bias)
from simbot.io import sha256_payload

T_GRID = [293.15, 295.0, 298.15, 300.0, 303.15, 306.0, 308.15]


# ═══════════════════════════════════════════════════════════════════════════
# MERGED — assert identity, not equality
# ═══════════════════════════════════════════════════════════════════════════
def test_boltzmann_constant_is_one_object():
    assert SU.K_B == BC.K_B
    assert BM.kB.to("J/K").magnitude == BC.K_B


@pytest.mark.parametrize("name", ["dt_max_thermal", "dt_max_force",
                                  "dt_max_active", "dt_max_stability"])
def test_dt_gates_are_the_same_function(name):
    """Not "give the same answer" -- literally the same object.

    Before the merge these existed twice with **different criteria** (displacement
    in `simbot`, timescale ratio in `bdbot.checks`), which
    `.claude/rules/overdamped-stability.md` forbids.
    """
    assert getattr(SN, name) is getattr(BDT, name)


@pytest.mark.parametrize("name", ["configurational_temperature", "check_finite",
                                  "check_inside_box", "check_bond_lengths",
                                  "check_step_displacements",
                                  "assert_statistic_fluctuates",
                                  "DisplacementReport"])
def test_guard_primitives_are_the_same_object(name):
    import simbot.guards as SG
    assert getattr(SG, name) is getattr(BH, name)


def test_em_bias_is_the_same_function():
    assert euler_maruyama_trap_variance_bias is BDT.em_variance_bias
    assert dt_star_for_trap_bias is BDT.dt_star_for_em_bias


@pytest.mark.parametrize("T", T_GRID)
def test_water_viscosity_agrees_across_packages(T):
    """One table. `bdbot` returns pint, `simbot` returns a float; same number."""
    si, extrapolated = SU.water_viscosity_si(T)
    assert not extrapolated
    assert BM.water_viscosity(Q(T, "K")).to("Pa*s").magnitude == si


@pytest.mark.parametrize("T", T_GRID)
def test_water_density_agrees_across_packages(T):
    si, _ = SU.water_density_si(T)
    assert BM.water_density(Q(T, "K")).to("kg/m^3").magnitude == si


@pytest.mark.parametrize("T,d_um", [(293.15, 1.0), (298.15, 0.5), (300.0, 5.0),
                                    (303.15, 2.0)])
def test_drag_and_diffusion_agree_across_conventions(T, d_um):
    """⚠ The factor-2 trap, pinned.

    `bdbot.materials.sphere_drag` takes a **diameter** (`3*pi*eta*d`),
    `simbot.units.stokes_drag_si` takes a **radius** (`6*pi*eta*a`). Same relation,
    and the formulas are deliberately NOT routed through a shared kernel -- pint
    round-tripping on a quantity that feeds `run_id` hashes and a 1e-12 round-trip
    check is a last-bit risk for no gain. This test is what makes that safe.
    """
    eta_si, _ = SU.water_viscosity_si(T)
    d_si = d_um * 1e-6

    gamma_bd = BM.sphere_drag(Q(eta_si, "Pa*s"), Q(d_si, "m")).to("kg/s").magnitude
    gamma_sb = SU.stokes_drag_si(eta_si, d_si / 2.0)
    assert gamma_bd == pytest.approx(gamma_sb, rel=1e-15)

    kT_bd = BM.thermal_energy(Q(T, "K")).to("J").magnitude
    assert kT_bd == pytest.approx(SU.kT_si(T), rel=1e-15)

    D_bd = BM.diffusion(Q(kT_bd, "J"), Q(gamma_bd, "kg/s")).to("m^2/s").magnitude
    assert D_bd == pytest.approx(SU.stokes_einstein_D_si(T, gamma_sb), rel=1e-15)


def test_radius_diameter_confusion_would_be_caught():
    """Deliberately break it and see: feeding a diameter where a radius belongs
    must show up as a factor of exactly 2, not as a rounding difference."""
    eta_si, _ = SU.water_viscosity_si(298.15)
    d_si = 1e-6
    correct = SU.stokes_drag_si(eta_si, d_si / 2.0)
    wrong = SU.stokes_drag_si(eta_si, d_si)
    assert wrong / correct == pytest.approx(2.0, rel=1e-15)


def test_water_viscosity_dict_is_a_view_over_the_table():
    """`bdbot.materials.WATER_VISCOSITY` must not be an independent third copy."""
    assert BM.WATER_VISCOSITY[298].to("Pa*s").magnitude == BC.WATER_ETA_SI[298.15]
    assert BM.WATER_VISCOSITY[293].to("Pa*s").magnitude == BC.WATER_ETA_SI[293.15]
    # 300 K is the separately-sourced anchor, NOT the table value -- see below.
    assert BM.WATER_VISCOSITY[300].to("Pa*s").magnitude == BC.WATER_ETA_SOURCED_SI[300.0]


# ═══════════════════════════════════════════════════════════════════════════
# UNMERGED — pin the measured gap so nobody "fixes" it silently
# ═══════════════════════════════════════════════════════════════════════════
def test_water_viscosity_provenance_gap_at_300K_is_pinned():
    """Two anchors, both in use, 0.545 % apart. **Do not reconcile this in code.**

    The 8 `bdbot` cases ran 0.851e-3 (written into their spec files and therefore
    into their `run_id`s); the `simbot` S2 documents sealed the interpolated
    8.5566e-4. Changing either breaks a content hash or a `SEALED.sha256`.
    `.claude/rules/verify-against-literature.md`: do not average a disagreement.
    """
    gap = BC.water_viscosity_provenance_gap(300.0)
    assert gap is not None
    assert gap["eta_table_si"] == pytest.approx(8.5566e-4, rel=1e-4)
    assert gap["eta_sourced_si"] == 0.851e-3
    assert gap["rel_gap"] == pytest.approx(0.00548, abs=1e-5)


def test_the_sealed_300K_value_is_unchanged():
    """The S2 prediction documents contain this number. If it moves, a seal breaks."""
    eta, extrapolated = SU.water_viscosity_si(300.0)
    assert not extrapolated
    assert eta == pytest.approx(8.5566e-4, rel=1e-4)


def test_temperature_uncertainty_dwarfs_the_anchor_disagreement():
    """The proportion that keeps the 0.545 % from being over-analysed.

    `T = 300 K` is tier 1 across the cases but was never measured -- it was
    inherited from a sketch that stated no temperature. One kelvin is worth more
    than four times the gap being argued about.
    """
    per_K = BC.water_viscosity_sensitivity_per_K(300.0)
    gap = BC.water_viscosity_provenance_gap(300.0)["rel_gap"]
    assert per_K == pytest.approx(0.02202, abs=1e-5)
    assert per_K / gap > 4.0


def test_the_temperature_convention_is_298_15_not_298():
    """⚠ `25 C = 298.15 K`. Pinning the 0.15 K, because it cost two sessions a round.

    Both this file's author and a peer session published a viscosity sensitivity
    computed at `298.00 K` in the same exchange in which one of them diagnosed that
    exact slip in the other's work. It is worth 0.3 percentage points -- small
    enough to survive review, large enough to matter next to a 0.545 % argument.

    At 298.15 K the IAPWS table has a **direct row**, so the interpolation axis
    vanishes; at 298.00 K it does not, and the answer becomes interpolation-
    dependent. That is the practical reason the convention matters.
    """
    #  ⚠ `water_eta_row` rather than `WATER_ETA_SI[298.15]`: a missing row must
    #    **report** ("no direct row at 298.0 K; nearest 298.15") instead of raising
    #    a bare `KeyError` that says nothing about why. A peer session measured that
    #    this exact pattern made an adversarial fork die rather than report.
    #  ⚠ And no Kelvin literals -- `BC.celsius()` is the only place the 273.15 lives.
    assert BC.T_25C == 298.15 and BC.T_20C == 293.15
    eta_25 = BC.water_eta_row(BC.T_25C)            # raises with a diagnosis, not KeyError
    eta_20 = BC.water_eta_row(BC.T_20C)

    at_25C, extrap = BC.water_viscosity_si(BC.T_25C)
    at_298_00, _ = BC.water_viscosity_si(BC.KELVIN_0C + 24.85)   # the wrong convention
    assert not extrap
    assert at_25C == eta_25                        # the row, exactly
    assert at_298_00 != at_25C                     # 298.00 interpolates

    base = BC.water_viscosity_si(300.0)[0]
    assert (at_25C - base) / base == pytest.approx(0.0401, abs=5e-5)
    assert (at_298_00 - base) / base == pytest.approx(0.0440, abs=5e-5)

    #  and what 0.851e-3 is worth against the truth at each convention
    anchor = BC.WATER_ETA_SOURCED_SI[300.0]
    assert (anchor - at_25C) / at_25C == pytest.approx(-0.0438, abs=5e-5)
    assert (anchor - eta_20) / eta_20 == pytest.approx(-0.1504, abs=5e-5)


def test_a_coarse_table_should_not_be_interpolated_when_a_fine_one_exists():
    """The other half of the same lesson, as a measurement.

    Log-interpolating the Welty 20 K table to 25 C gives 893.157 uPa*s against the
    IAPWS **row** value of 890.000 -- 0.35 %. Linear gives 906.738, i.e. 1.9 %. Both
    are avoidable: `WATER_ETA_SI` has the row. The four labelling choices together
    span -4.38 % to -6.41 %, which is why this module refuses to quote a
    sensitivity without naming all four.
    """
    #  ⚠ The rows come from `bdbot.constants`, NOT copied into this test body. An
    #    earlier version of this function copied the three Welty rows inline -- i.e.
    #    a de-duplication test that was itself a duplicate. A peer session hit the
    #    same thing from the other side, adding a fourth copy of the water table to
    #    a verifier while fixing a bug caused by copies of that table diverging.
    welty = BC.WATER_ETA_HANDBOOK_20K_ROWS_SI
    welty_log = BC.interp_table_log(welty, BC.T_25C)[0]
    welty_lin = BC.interp_table(welty, BC.T_25C)[0]
    row = BC.water_eta_row(BC.T_25C)

    assert welty_log == pytest.approx(893.157e-6, rel=1e-5)
    assert welty_lin == pytest.approx(906.738e-6, rel=1e-5)
    assert abs(welty_log - row) / row == pytest.approx(0.00355, abs=5e-5)
    assert abs(welty_lin - row) / row == pytest.approx(0.01881, abs=5e-5)
    #  the spread across all four labelling choices, at ONE temperature
    anchor = BC.WATER_ETA_SOURCED_SI[300.0]
    spreads = [(anchor - e) / e for e in (row, welty_log, welty_lin)]
    assert min(spreads) == pytest.approx(-0.0615, abs=5e-4)
    assert max(spreads) == pytest.approx(-0.0438, abs=5e-4)


def test_sensitivity_is_a_segment_value_not_a_table_constant():
    """⚠ Pins which segment each figure comes from.

    Three sessions computed this quantity on 2026-08-29 and got three answers,
    entirely from leaving the table and the baseline implicit. `2.062 %/K` is a
    value of **both** tables (the 303-308 K segment here, the 293-313 K segment of
    the Welty 20 K table), which is exactly why the confusion was hard to see.
    """
    per_segment = {T: BC.water_viscosity_sensitivity_per_K(T)
                   for T in (295.0, 300.0, 305.0)}
    assert per_segment[295.0] == pytest.approx(0.02363, abs=1e-5)
    assert per_segment[300.0] == pytest.approx(0.02202, abs=1e-5)
    assert per_segment[305.0] == pytest.approx(0.02062, abs=1e-5)
    # it varies by >14 % across the table -- so it is not a constant to quote bare
    assert max(per_segment.values()) / min(per_segment.values()) > 1.14


def test_the_handbook_anchor_is_internally_consistent():
    """The typo that a 56/56-passing verifier could not see.

    `welty_transport.md` printed `0.8580 mPa*s` next to its own `+1.03 %`, and
    +1.03 % of 0.851 is 0.8598. The check asserted only that the *computed* value
    was within 1.5 % of 0.851 -- true at 1.03 % and at 0.82 % alike -- and never
    compared it to the *printed* one. Fixed upstream in `20da7a9`; pinned here so
    the internal consistency is what gets checked, not the tolerance.
    """
    log_linear = BC.WATER_ETA_HANDBOOK_20K_SI["log_linear"]
    linear = BC.WATER_ETA_HANDBOOK_20K_SI["linear"]
    anchor = BC.WATER_ETA_SOURCED_SI[300.0]
    assert (log_linear - anchor) / anchor == pytest.approx(0.0103, abs=5e-5)
    assert (linear - anchor) / anchor == pytest.approx(0.0291, abs=5e-5)
    #  and the value 0.8580 would NOT have been consistent with +1.03 %
    assert (8.580e-4 - anchor) / anchor == pytest.approx(0.0082, abs=5e-5)


@pytest.mark.parametrize("payload,expect_same", [
    ({"d": 5e-6, "eta": 8.51e-4, "T": 300}, True),
    ({"N": 512, "k_star": 1.0}, True),
    ({"d": 5e-6, "src": "물@300K 핸드북"}, False),
    ({"점도": 8.51e-4}, False),
])
def test_the_two_spec_hashes_agree_on_ascii_and_differ_otherwise(payload, expect_same):
    """⚠ `bdbot.runid.spec_hash` and `simbot.io.sha256_payload` are NOT the same.

    The single difference is `ensure_ascii`. They agree on every ASCII payload,
    which is why the divergence never surfaced -- and this repository is full of
    Korean strings. Deliberately unmerged: 263 run directories and 279 spec files
    are named by the `bdbot` one.

    What holds the ASCII agreement together is `physics_only`, which strips
    `source`, `note`, `description` and the rest of `DOC_KEYS` -- exactly the
    prose-carrying fields -- before `spec_hash` ever sees a spec.
    """
    a = spec_hash(payload, 12)
    b = sha256_payload(payload)[:12]
    assert (a == b) is expect_same
    if not expect_same:
        # and the cause is exactly ensure_ascii, not something else
        assert a == __import__("hashlib").sha256(
            json.dumps(payload, sort_keys=True, default=str).encode()
        ).hexdigest()[:12]


def test_em_bias_linearized_and_exact_relationship_is_pinned():
    """`bdbot.checks.dt_from_bias` keeps the first-order inverse on purpose.

    `cases/trap_2d_5um.py` picks its `dt` from it and `run_id` hashes the spec, so
    the exact form would rename existing runs.

    The gap is **exactly the target bias `b`** -- measured here, not asserted. A
    first pass at this docstring said `b/(1+b)`, which is wrong by `b^2` and would
    never have shown up in any output. That is the whole argument for pinning an
    "obviously equivalent" pair with a test instead of a comment.
    """
    for b in (1e-4, 1e-3, 1e-2, 2.5e-3):
        exact = BDT.dt_star_for_em_bias(b)
        lin = BDT.dt_star_for_em_bias_linearized(b)
        assert lin == pytest.approx(2.0 * b, rel=1e-15)
        assert BDT.em_bias_form_gap(b) == pytest.approx(abs(lin - exact) / exact, rel=1e-12)
        assert BDT.em_bias_form_gap(b) == pytest.approx(b, rel=1e-12)
        assert BDT.em_bias_form_gap(b) != pytest.approx(b / (1.0 + b), rel=1e-9)
    # the exact form round-trips; the linearized one does not, and by how much
    assert BDT.em_variance_bias(BDT.dt_star_for_em_bias(1e-3)) == pytest.approx(1e-3, rel=1e-12)
    assert BDT.em_variance_bias(1e-2) == pytest.approx(0.0050251, rel=1e-5)
    assert BDT.em_variance_bias_linearized(1e-2) == 0.005


def test_pint_wrappers_reproduce_the_old_inline_formulas():
    """`bdbot.checks`' dt helpers now call `bdbot.dt`. Bit-identical, by execution.

    `bias_from_dt` went from `50.0*r` to `100.0*(r/2)`, and `dt_from_bias` from
    `2*bias*tau` to `(2.0*bias)*tau`. Both are the same real number, but "the same
    real number" is not the same as "the same float" -- so measure it.
    """
    import bdbot.checks as BCK

    tau = Q(1.234567e-3, "s")
    for bias in (1e-4, 1e-3, 5e-3, 1e-2):
        assert BCK.dt_from_bias(tau, bias).to("s").magnitude == (2 * bias * tau).to("s").magnitude
    for dt_over_tau in (1e-4, 1e-3, 1e-2, 3.7e-3):
        dt = dt_over_tau * tau
        #  ⚠ Compare against the OLD formula on the **same** ratio, not against a
        #    literal. `(dt/tau)` through pint is not bit-exactly `dt_over_tau`, and
        #    a first pass at this test compared to `50.0*dt_over_tau` and "failed"
        #    on a pint artefact rather than on a real change.
        r = float((dt / tau).to("dimensionless").magnitude)
        assert BCK.bias_from_dt(dt, tau) == 50.0 * r
    assert BCK.dt_from_gate(tau).to("s").magnitude == (BCK.GATE * tau).to("s").magnitude
    assert BCK.GATE is BDT.GATE


@pytest.mark.parametrize("r", [1e-8, 1e-4, 1e-3, 1e-2, 3.7e-3, 0.1, 1 / 3, 0.999])
def test_the_rewritten_bias_expression_is_bit_identical(r):
    """`50.0*r` became `100.0*(r/2.0)`. Same real number -- but that is not the same
    claim as the same float, so it is measured over a probe set (200k random draws
    found zero disagreements) rather than argued from binary exactness."""
    assert 100.0 * (r / 2.0) == 50.0 * r
