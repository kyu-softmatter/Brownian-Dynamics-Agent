"""S7 — proposing a verdict.

## The central test in this file

`test_reproduces_first_run_verdicts` — it **reproduces the first completed run's 9
verdicts exactly** (7 PASS · 2 INCONCLUSIVE · 0 FAIL). If the verdict a human
reached by hand and the verdict the code reaches diverge, one of the two is wrong,
and either way that has to be known.

The input numbers come from `runs/2026-07-28_trap-2d-5um_2dfb9d/metrics.json`.
`runs/` is gitignored, so rather than depending on that file **the values are pinned
here** -- a regression test has to run from any checkout.
"""
from __future__ import annotations

import math

import pytest

from simbot import validate as V
from simbot.io import RunDir, write_seal
from simbot.spec import Prediction, PredictionItem
from simbot.validate import FAIL, INCONCLUSIVE, PASS, Measurement, compare

# =============================================================================
# parsing a tolerance
# =============================================================================
@pytest.mark.parametrize("text,kind,mag", [
    ("±1.5%", "relative", 0.015),
    ("± 2 %", "relative", 0.02),
    ("+/-0.03", "absolute", 0.03),
    ("±1e-3", "absolute", 1e-3),
    (">0.99", "lower_bound", 0.99),
    ("p>0.05", "lower_bound", 0.05),
    ("R^2>0.99", "lower_bound", 0.99),
    ("<1e-12", "upper_bound", 1e-12),
    (">=3", "lower_bound", 3.0),
])
def test_parse_tolerance_forms(text, kind, mag):
    t = V.parse_tolerance(text)
    assert t.kind == kind
    assert t.magnitude == pytest.approx(mag)


def test_unparseable_tolerance_raises():
    """★ Skip an unreadable tolerance quietly and that item passes with no verdict."""
    with pytest.raises(ValueError, match="cannot be read"):
        V.parse_tolerance("roughly right is fine")


def test_zero_tolerance_raises():
    """A half-width of 0 makes every result FAIL — that is not verification."""
    with pytest.raises(ValueError, match="0 or less"):
        V.parse_tolerance("±0")


def test_relative_half_width_scales_with_prediction():
    t = V.parse_tolerance("±1.5%")
    assert t.half_width(414.19) == pytest.approx(414.19 * 0.015)
    assert t.half_width(None) is None


# =============================================================================
# measurements — a number with no error cannot go in a conclusion
# =============================================================================
def test_measurement_without_stat_err_is_a_problem():
    assert any("no statistical error" in p
               for p in Measurement("D", 1.0).problems())


def test_zero_stat_err_flags_identity_risk():
    """A 'measurement' that does not fluctuate is an arithmetic identity
    (an actual case, 2026-07-28)."""
    assert any("identity" in p for p in Measurement("x", 1.0, stat_err=0.0).problems())


def test_nonfinite_measurement_is_a_problem():
    assert any("non-finite" in p
               for p in Measurement("x", math.nan, stat_err=0.1).problems())


def test_verdict_is_inconclusive_without_error_bar():
    item = PredictionItem("D", 1.0, "±3%", "Stokes-Einstein")
    row = compare(item, Measurement("D", 1.0))
    assert row.verdict == INCONCLUSIVE and "unbarred number" in row.reason


# =============================================================================
# verdict logic
# =============================================================================
def test_pass_when_inside_band_and_decidable():
    item = PredictionItem("var_x", 414.19, "±1.5%", "kT/k")
    row = compare(item, Measurement("var_x", 416.58, stat_err=1.85))
    assert row.verdict == PASS
    assert row.deviation_rel == pytest.approx(0.00577, rel=1e-2)


def test_fail_when_outside_band():
    item = PredictionItem("var_x", 414.19, "±1.5%", "kT/k")
    row = compare(item, Measurement("var_x", 460.0, stat_err=1.85))
    assert row.verdict == FAIL


def test_inconclusive_when_se_exceeds_tolerance():
    """★ If SE exceeds the band, this measurement cannot decide this tolerance.

    Writing that as PASS is luck, not verification -- it landed inside not because
    the band is wide but because **the measurement cannot resolve that band.**
    """
    item = PredictionItem("var_x_star", 1.00251, "±0.1%", "EM bias")
    row = compare(item, Measurement("var_x_star", 1.00577, stat_err=0.00446))
    assert row.verdict == INCONCLUSIVE
    assert "cannot decide this band" in row.reason
    assert "sample multiple needed" in row.reason


def test_required_sample_multiple_is_the_variance_ratio():
    item = PredictionItem("x", 1.0, "±0.1%", "b")
    row = compare(item, Measurement("x", 1.0, stat_err=0.01))
    # SE/half = 0.01/0.001 = 10 → 100x
    assert "100" in row.reason


def test_inside_band_but_low_power_is_inconclusive():
    """★ CLAUDE.md's 4 statistics rules — do not demand 3σ where the power cannot
    produce 3σ.

    Landing inside the band still discriminates nothing if it cannot be told apart
    from the competing hypothesis. `INCONCLUSIVE` is pinned as a fact.
    """
    # The predicted value uses the full-precision 1/(1−dt*/2) — a truncated value
    # gives a power of 0.5628 and disagrees with the document's 0.5618
    item = PredictionItem("var_x_star", 1.0025063, "±1%", "EM bias",
                          competing_value=1.0)      # the exact scheme
    row = compare(item, Measurement("var_x_star", 1.0057652, stat_err=0.0044614))
    assert row.verdict == INCONCLUSIVE
    assert row.design_power == pytest.approx(0.5618, rel=1e-3)
    assert "a foreseen limit, not a failure" in row.reason


def test_high_power_gives_pass():
    """Raise dt* and discriminating power appears — the first run's `dt*=2e-2`."""
    item = PredictionItem("var_x_star", 1.01010, "±1%", "EM bias",
                          competing_value=1.0)
    row = compare(item, Measurement("var_x_star", 1.01041, stat_err=0.00264))
    assert row.verdict == PASS
    assert row.design_power > 3.0


def test_samples_needed_for_3sigma_is_reported():
    item = PredictionItem("x", 1.00251, "±1%", "b", competing_value=1.0)
    row = compare(item, Measurement("x", 1.0, stat_err=0.00446))
    assert row.samples_needed_for_3sigma == pytest.approx((3 / 0.5617) ** 2, rel=1e-2)


# =============================================================================
# one-sided bounds
# =============================================================================
def test_lower_bound_pass():
    item = PredictionItem("msd_r_squared", ">0.99", ">0.99", "single-exponential form")
    row = compare(item, Measurement("msd_r_squared", 0.999977, stat_err=1e-5))
    assert row.verdict == PASS


def test_lower_bound_fail():
    item = PredictionItem("ks_p", "p>0.05", "p>0.05", "Gaussian")
    row = compare(item, Measurement("ks_p", 0.0000, stat_err=0.01))
    assert row.verdict == FAIL


def test_bound_too_close_to_call_is_inconclusive():
    item = PredictionItem("ks_p", "p>0.05", "p>0.05", "Gaussian")
    row = compare(item, Measurement("ks_p", 0.055, stat_err=0.02))
    assert row.verdict == INCONCLUSIVE
    assert "indistinguishable from the bound" in row.reason


def test_upper_bound_pass():
    item = PredictionItem("roundtrip_err", "<1e-12", "<1e-12", "the S4 gate")
    row = compare(item, Measurement("roundtrip_err", 1.6e-16, stat_err=1e-18))
    assert row.verdict == PASS


def test_non_numeric_prediction_with_two_sided_band_is_inconclusive():
    item = PredictionItem("shape", "single exponential", "±5%", "shape")
    assert compare(item, Measurement("shape", 1.0, stat_err=0.01)).verdict == INCONCLUSIVE


# =============================================================================
# ★ reproducing the first completed run's 9 verdicts
# =============================================================================
#  Values from runs/2026-07-28_trap-2d-5um_2dfb9d/{02_prediction.md, metrics.json}.
#  runs/ is gitignored, so they are pinned here.
FIRST_RUN = [
    # (prediction item, measurement, the verdict recorded in the document)
    (PredictionItem("var_x_2d_nm2", 414.19, "±1.5%", "⟨x²⟩ = kT/k (exact)"),
     Measurement("var_x_2d_nm2", 416.5826, stat_err=1.8479, unit="nm^2"), PASS),
    (PredictionItem("var_x_3d_nm2", 414.19, "±1.5%", "⟨x²⟩ = kT/k, dim-independent"),
     Measurement("var_x_3d_nm2", 415.3363, stat_err=1.3685, unit="nm^2"), PASS),
    (PredictionItem("var_x_star", 1.0025063, "±1%", "EM bias 1/(1−dt*/2)"),
     Measurement("var_x_star", 1.0057652, stat_err=0.0044614), PASS),
    (PredictionItem("var_r_2d_nm2", 828.39, "±1.5%", "⟨r²⟩ = d·kT/k"),
     Measurement("var_r_2d_nm2", 833.1652, stat_err=3.6958, unit="nm^2"), PASS),
    (PredictionItem("msd_plateau_2d_star", 4.0, "±2%", "plateau = 2d"),
     Measurement("msd_plateau_2d_star", 4.0214904, stat_err=0.0060703), PASS),
    (PredictionItem("tau_trap_ms", 8.0644, "±5%", "τ = γ/k"),
     Measurement("tau_trap_ms", 8.0566504, stat_err=0.0299748, unit="ms"), PASS),
    (PredictionItem("msd_r_squared", ">0.99", ">0.99", "single-exponential form"),
     Measurement("msd_r_squared", 0.9999765, stat_err=1.0e-5), PASS),
    (PredictionItem("ks_p", "p>0.05", "p>0.05", "Gaussian position distribution"),
     Measurement("ks_p", 0.29, stat_err=0.05), PASS),
    # P8 — insufficient power. Document: 0.56σ
    (PredictionItem("em_bias_reproduced", 1.0025063, "±0.1%",
                    "EM bias reproduced",
                    competing_value=1.0, discriminates="integrator scheme"),
     Measurement("em_bias_reproduced", 1.0057652, stat_err=0.0044614), INCONCLUSIVE),
    # P9 — the dt ladder. Document: 0.24σ
    (PredictionItem("em_bias_halved", 1.0012516, "±0.1%",
                    "dt halved → bias halved",
                    competing_value=1.0, discriminates="dt-linearity of the bias"),
     Measurement("em_bias_halved", 1.0016875, stat_err=0.0051896), INCONCLUSIVE),
]


@pytest.mark.benchmark
@pytest.mark.parametrize("item,meas,expected", FIRST_RUN,
                         ids=[i.quantity for i, _, _ in FIRST_RUN])
def test_reproduces_first_run_verdicts(item, meas, expected):
    """★ The code reproduces the verdict a human reached by hand."""
    row = compare(item, meas)
    assert row.verdict == expected, row.reason


@pytest.mark.benchmark
def test_first_run_verdict_counts():
    """The document's tally: 8 PASS · 2 INCONCLUSIVE · 0 FAIL (P1 and P3 counted
    per item)."""
    verdicts = [compare(i, m).verdict for i, m, _ in FIRST_RUN]
    assert verdicts.count(FAIL) == 0
    assert verdicts.count(INCONCLUSIVE) == 2
    assert verdicts.count(PASS) == len(FIRST_RUN) - 2


@pytest.mark.benchmark
def test_p8_design_power_matches_recorded_value():
    """Reproduces the recorded power `0.56σ` and `29x for 3σ`."""
    item, meas, _ = FIRST_RUN[8]
    row = compare(item, meas)
    assert row.design_power == pytest.approx(0.5618, rel=2e-3)
    assert row.samples_needed_for_3sigma == pytest.approx(28.5, rel=0.05)


# =============================================================================
# sealing — a broken seal produces no comparison table
# =============================================================================
@pytest.fixture
def sealed_rundir(tmp_path):
    rd = RunDir.create(tmp_path, "r1")
    rd.write("intake", "# S1\n")
    rd.write("prediction", "# S2\nvar_x = 414.19 nm^2 ±1.5%\n")
    rd.write("spec", "card: passive-sphere--harmonic-trap\n")
    write_seal(rd)
    return rd


def _prediction():
    return Prediction(items=[PredictionItem("var_x", 414.19, "±1.5%", "kT/k")])


def _measurements():
    return {"var_x": Measurement("var_x", 416.58, stat_err=1.85)}


def test_validate_run_passes_with_intact_seal(sealed_rundir):
    rep = V.validate_run(_prediction(), _measurements(), rundir=sealed_rundir)
    assert rep.seal.ok
    assert rep.verdict_overall == PASS
    assert rep.problems == []


def test_broken_seal_produces_no_comparison_table(sealed_rundir):
    """★ If the prediction may have been edited after seeing the result, the
    comparison table is not a verification."""
    sealed_rundir.write("prediction",
                        "# S2\nvar_x = 416 nm^2 ±0.5%  (edited after the result)\n")
    rep = V.validate_run(_prediction(), _measurements(), rundir=sealed_rundir)
    assert rep.rows == []
    assert rep.verdict_overall == "SEAL_BROKEN"
    assert any("seal violation" in p for p in rep.problems)


def test_missing_measurement_is_reported(sealed_rundir):
    """A sealed prediction cannot be quietly dropped."""
    rep = V.validate_run(_prediction(), {}, rundir=sealed_rundir)
    assert any("no measurement matches the prediction" in p for p in rep.problems)


def test_unsealed_extra_measurement_is_reported(sealed_rundir):
    meas = {**_measurements(), "kurtosis": Measurement("kurtosis", 2.9968,
                                                       stat_err=0.012)}
    rep = V.validate_run(_prediction(), meas, rundir=sealed_rundir)
    assert any("unsealed measurements" in p for p in rep.problems)


def test_empty_prediction_is_reported(sealed_rundir):
    rep = V.validate_run(Prediction(items=[]), {}, rundir=sealed_rundir)
    assert any("0개" in p for p in rep.problems)


# =============================================================================
# the verdict block — code never fills in confirmed_by
# =============================================================================
def test_yaml_block_never_confirms():
    """★ If this test does not pass, a pass stamp gets applied that no human saw."""
    rep = V.ValidationReport(rows=[compare(*FIRST_RUN[0][:2])])
    block = rep.yaml_block()
    assert "confirmed_by: null" in block
    assert "proposed_by: agent" in block
    assert "confirmed_by: agent" not in block


def test_no_public_api_sets_confirmed_by():
    """No code path that fills in `confirmed_by` may exist."""
    src = (V.__file__)
    text = open(src, encoding="utf-8").read()
    for bad in ("confirmed_by =", 'confirmed_by": "', "confirmed_by: agent"):
        assert bad not in text, bad


def test_overall_verdict_prefers_fail_over_inconclusive():
    rows = [compare(*FIRST_RUN[8][:2]),                       # INCONCLUSIVE
            compare(PredictionItem("x", 1.0, "±1%", "b"),
                    Measurement("x", 2.0, stat_err=0.001))]   # FAIL
    assert V.ValidationReport(rows=rows).verdict_overall == FAIL


def test_overall_verdict_marks_inconclusive():
    rep = V.ValidationReport(rows=[compare(*FIRST_RUN[0][:2]),
                                   compare(*FIRST_RUN[8][:2])])
    assert rep.verdict_overall == "PASS_WITH_INCONCLUSIVE"
    assert rep.count(PASS) == 1 and rep.count(INCONCLUSIVE) == 1


def test_fail_without_cause_class_is_a_problem(sealed_rundir):
    """A FAIL must carry one of the 4 cause classes (S7 gate)."""
    pred = Prediction(items=[PredictionItem("var_x", 414.19, "±1.5%", "kT/k")])
    meas = {"var_x": Measurement("var_x", 600.0, stat_err=1.85)}
    rep = V.validate_run(pred, meas, rundir=sealed_rundir)
    assert rep.verdict_overall == FAIL
    assert any("no cause class" in p for p in rep.problems)


def test_fail_with_cause_class_has_no_problem(sealed_rundir):
    pred = Prediction(items=[PredictionItem("var_x", 414.19, "±1.5%", "kT/k")])
    meas = {"var_x": Measurement("var_x", 600.0, stat_err=1.85)}
    rep = V.validate_run(pred, meas, rundir=sealed_rundir,
                         causes={"var_x": "numerical"})
    assert rep.problems == []


def test_table_and_reasons_render(sealed_rundir):
    rep = V.validate_run(_prediction(), _measurements(), rundir=sealed_rundir)
    assert "| quantity | predicted (sealed) |" in rep.table()
    assert "PASS" in rep.table()
    assert rep.reasons() == "_none — every item PASS._"


def test_reasons_explain_inconclusive():
    rep = V.ValidationReport(rows=[compare(*FIRST_RUN[8][:2])])
    assert "INCONCLUSIVE" in rep.reasons()
    assert "a foreseen limit" in rep.reasons()
