"""S8 — generating `REPORT.md`.

The report is an **audit record**, not a summary. So what this file checks is not
"does it come out pretty" but **does any bad news go missing**:

- when the seal is broken, is the comparison table withheld and a warning put on top
- do the `INCONCLUSIVE` items and their reasons survive
- is `git_dirty` reported
- is an uncaptioned figure flagged
- is `confirmed_by: null` present
"""
from __future__ import annotations

from pathlib import Path

import pytest

from simbot import report as R
from simbot.io import RunDir, build_manifest, write_seal
from simbot.spec import Prediction, PredictionItem, SystemSpec
from simbot.validate import Measurement, validate_run

EXAMPLE_SPEC = Path("examples/trap-2d-5um/spec.yaml")


@pytest.fixture
def spec() -> SystemSpec:
    if not EXAMPLE_SPEC.exists():
        pytest.skip("examples/trap-2d-5um/spec.yaml is absent")
    return SystemSpec.load(EXAMPLE_SPEC)


@pytest.fixture
def rundir(tmp_path, spec) -> RunDir:
    rd = RunDir.create(tmp_path, "2026-07-28_trap-2d-5um_c59e93")
    rd.write("spec", spec.to_yaml())
    rd.write("intake", "# S1 INTAKE\n")
    rd.write("prediction", "# S2 PREDICTION\n")
    rd.write("conclusion", "# S8 CONCLUSION\n")
    write_seal(rd)
    return rd


def _prediction() -> Prediction:
    return Prediction(items=[
        PredictionItem("var_x_2d_nm2", 414.19, "±1.5%", "⟨x²⟩ = kT/k (exact)"),
        PredictionItem("em_bias_reproduced", 1.0025063, "±0.1%",
                       "EM bias reproduced",
                       competing_value=1.0, discriminates="integrator scheme"),
    ])


def _measurements() -> dict[str, Measurement]:
    return {
        "var_x_2d_nm2": Measurement("var_x_2d_nm2", 416.5826, stat_err=1.8479,
                                    unit="nm^2", n_samples=4),
        "em_bias_reproduced": Measurement("em_bias_reproduced", 1.0057652,
                                          stat_err=0.0044614, n_samples=4),
    }


@pytest.fixture
def full_inputs(rundir, spec) -> R.ReportInputs:
    rep = validate_run(_prediction(), _measurements(), rundir=rundir)
    man = build_manifest(run_id=rundir.run_id, spec_hash=spec.hash(),
                         seed=[1, 2, 3, 4], rundir=rundir)
    return R.ReportInputs(spec=spec, validation=rep, manifest=man,
                          wall_s=10.6, n_runs=16)


# =============================================================================
# whole-report render
# =============================================================================
def test_report_renders_and_writes(rundir, full_inputs):
    p = R.write_report(rundir, full_inputs)
    assert p.exists()
    text = p.read_text(encoding="utf-8")
    assert text.startswith("# REPORT — `2026-07-28_trap-2d-5um_c59e93`")


def test_report_never_confirms(rundir, full_inputs):
    """★ A human confirmation stamp must never appear in the report."""
    text = R.render(rundir, full_inputs)
    assert "confirmed_by: null" in text
    assert "confirmed_by: agent" not in text


def test_report_states_verdict_is_a_proposal(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "The verdict is a **proposal**" in text


def test_report_is_self_contained_on_all_stages(rundir, full_inputs):
    """§1-§6 must all be present for a person to judge from this file alone."""
    text = R.render(rundir, full_inputs)
    for head in ("## 1. Verdict summary",
                 "## 2. System specification and gates (S3)",
                 "## 3. Non-dimensionalization (S4)", "## 4. Figures (S6)",
                 "## 5. Reproducibility", "## 6. Agent-written documents"):
        assert head in text, head


def test_report_survives_missing_pieces(rundir):
    """Empty inputs still render, and what is absent comes out as 'none'."""
    text = R.render(rundir, R.ReportInputs())
    assert "no verdict" in text
    assert "_no S7 verdict._" in text
    assert "_no S3 check results._" in text


# =============================================================================
# sealing — a broken seal withholds the comparison table
# =============================================================================
def test_intact_seal_shows_external_verification_command(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "shasum -a 256 -c SEALED.sha256" in text


def test_broken_seal_suppresses_the_comparison_table(rundir, full_inputs):
    """★ A comparison table under a broken seal looks like a verification and is
    not one."""
    rundir.write("prediction", "# S2 (edited after seeing the result)\n")
    text = R.render(rundir, full_inputs)
    assert "⛔ seal violation" in text
    assert "The comparison table was **not generated**" in text
    assert "| quantity | predicted (sealed) |" not in text


def test_broken_seal_warning_appears_before_results(rundir, full_inputs):
    rundir.write("prediction", "# edited\n")
    text = R.render(rundir, full_inputs)
    assert text.index("⛔ seal violation") < text.index("## 1. Verdict summary")


# =============================================================================
# INCONCLUSIVE does not disappear
# =============================================================================
def test_inconclusive_and_its_reason_are_reported(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "INCONCLUSIVE" in text
    assert "design power" in text
    assert "28.5x the samples needed" in text or "28.5" in text


def test_inconclusive_is_framed_as_a_fact_not_a_failure(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "not a failure" in text


def test_counts_are_reported(rundir, full_inputs):
    assert "2 items: 1 PASS · 1 INCONCLUSIVE · 0 FAIL." in R.render(rundir,
                                                                   full_inputs)


def test_validation_problems_are_surfaced(rundir, spec):
    """A measurement with no error bar must surface in the report."""
    meas = {"var_x_2d_nm2": Measurement("var_x_2d_nm2", 416.6),   # no stat_err
            "em_bias_reproduced": Measurement("em_bias_reproduced", 1.0057, stat_err=0.004)}
    rep = validate_run(_prediction(), meas, rundir=rundir)
    text = R.render(rundir, R.ReportInputs(spec=spec, validation=rep))
    assert "problems with the verification procedure" in text
    assert "no statistical error" in text


# =============================================================================
# reproducibility
# =============================================================================
def test_dirty_git_is_reported_as_a_reproducibility_limit(rundir, spec):
    man = {"run_id": "r", "git_rev": "abc1234", "git_dirty": True, "env": {}}
    text = R.reproducibility_section(man)
    assert "uncommitted changes present" in text
    assert "does not reproduce from `git_rev` alone" in text


def test_clean_git_is_reported_as_clean():
    text = R.reproducibility_section({"git_rev": "abc", "git_dirty": False, "env": {}})
    assert "clean" in text


def test_unknown_git_state_is_not_claimed_clean():
    """Writing 'undecidable' as clean falsely claims reproducibility."""
    text = R.reproducibility_section({"git_rev": "?", "git_dirty": None, "env": {}})
    assert "undecidable" in text
    assert "clean" not in text


def test_missing_manifest_refuses_to_claim_reproducibility():
    assert "cannot be claimed" in R.reproducibility_section(None)


def test_manifest_records_versions(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "hoomd" in text and "7.1.0" in text


# =============================================================================
# gates
# =============================================================================
def test_deferred_gates_are_listed(rundir, full_inputs):
    """A person has to know how many gates S3 could not decide."""
    text = R.render(rundir, full_inputs)
    assert "gates S7 has to decide" in text
    assert "`equipartition`" in text


def test_off_gates_show_their_reason(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "no pair interaction, so overlap cannot occur" in text


def test_spec_problems_are_surfaced(rundir, spec):
    from simbot.spec import Q
    spec.medium.rho_fluid_si = Q(996.5, "kg/m^3", "assumed",
                                 "has a basis but no confidence")
    text = R.gates_section(R.validate_spec(spec))
    assert "convention violation" in text and "confidence" in text


# =============================================================================
# the non-dimensionalization section
# =============================================================================
def test_nondim_section_reports_roundtrip_gate(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "round-trip error" in text
    assert "✅ pass" in text


def test_nondim_section_separates_logged_from_gated(rundir, full_inputs):
    """The report must show that `dt/τ_D` is a record, not a gate."""
    text = R.render(rundir, full_inputs)
    assert "for the record (not a gate" in text
    assert "dt_over_tau_D" in text


def test_nondim_section_names_the_scale_origin(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "scales_harmonic_trap: (l_trap, kT, tau_trap)" in text


def test_nondim_section_empty_without_spec(rundir):
    assert "no S4 non-dimensionalization results" in R.nondim_section(None, None)


def test_unregistered_card_does_not_crash_the_report(rundir, spec):
    """The report must still come out for a card with no scale rule (only that
    section is empty)."""
    spec.card = "brand--new-pair"
    text = R.render(rundir, R.ReportInputs(spec=spec))
    assert "REPORT" in text
    assert "_no S4 non-dimensionalization results._" in text


# =============================================================================
# figures — an uncaptioned figure is not an artefact
# =============================================================================
def test_figure_without_caption_is_flagged(rundir, full_inputs):
    (rundir.figs / "01_msd.png").write_bytes(b"\x89PNG\r\n")
    text = R.render(rundir, full_inputs)
    assert "no caption" in text
    assert "§S6 gate" in text


def test_figure_with_caption_is_embedded(rundir, full_inputs):
    (rundir.figs / "01_msd.png").write_bytes(b"\x89PNG\r\n")
    full_inputs.figures = {"01_msd.png": "MSD and the solution 2d(1−e^{−t/τ})"}
    text = R.render(rundir, full_inputs)
    assert "![MSD and the solution" in text
    assert "no caption" not in text


def test_no_figures_says_so(rundir, full_inputs):
    assert "_no figures._" in R.render(rundir, full_inputs)


# =============================================================================
# cost
# =============================================================================
def test_cost_section_reports_per_run_average():
    assert "mean per run `0.66 s`" in R.cost_section(10.6, 16)


def test_missing_cost_says_so():
    assert "no compute-cost record" in R.cost_section(None, None)


# =============================================================================
# quoting the agent's documents — it does not write them instead
# =============================================================================
def test_missing_agent_document_is_named_not_invented(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "`07_validation.md` is missing" in text
    assert "has to be written by the agent" in text


def test_existing_agent_documents_are_linked(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "[`08_conclusion.md`](08_conclusion.md)" in text


def test_report_states_the_division_of_labour(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "does not write them instead" in text
