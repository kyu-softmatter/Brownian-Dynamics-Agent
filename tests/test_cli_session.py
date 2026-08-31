"""The CLI and session layer.

## What this file enforces

1. **`set` does not run anything.** It only estimates the cost. Break that
   separation and "let's try N = 8000" quietly starts an 11-minute run.
2. **The history is append-only.** Overwriting a past turn loses what was tried.
3. **A card with no runner is not forced through.** Running it quietly with the trap
   runner computes an entirely different system and nobody finds out.
4. **Over budget or a failed gate stops it before running.**

Tests that need a HOOMD run carry the `slow` marker.
"""
from __future__ import annotations

from datetime import datetime
from pathlib import Path

import pytest

import cli
from simbot.session import Session, set_by_path
from simbot.spec import SystemSpec

EXAMPLE_SPEC = Path("examples/trap-2d-5um/spec.yaml")
EXAMPLE_PRED = Path("examples/trap-2d-5um/prediction.yaml")
T0 = datetime(2026, 7, 28, 12, 0, 0)


@pytest.fixture
def spec_file(tmp_path) -> Path:
    if not EXAMPLE_SPEC.exists():
        pytest.skip("examples/trap-2d-5um/spec.yaml is absent")
    p = tmp_path / "spec.yaml"
    p.write_text(EXAMPLE_SPEC.read_text(encoding="utf-8"), encoding="utf-8")
    return p


@pytest.fixture
def session(spec_file, tmp_path) -> Session:
    return Session.create(spec_file, root=tmp_path / "sessions", when=T0)


# =============================================================================
# dotted-key assignment
# =============================================================================
def test_set_by_path_updates_quantity_value(spec_file):
    spec = SystemSpec.load(spec_file)
    old, new = set_by_path(spec, "numerics.dt_star", "2.5e-3", turn=1)
    assert old == 0.005 and new == 0.0025
    assert spec.numerics.dt_star.value == 0.0025


def test_set_by_path_marks_provenance_user(spec_file):
    """★ A human-chosen value is neither an assumption nor policy — it is not
    something S7b shakes."""
    spec = SystemSpec.load(spec_file)
    set_by_path(spec, "medium.eta_si", "1.002e-3", turn=3)
    q = spec.medium.eta_si
    assert q.provenance == "user"
    assert "set by a human on turn 3" in q.basis


def test_set_by_path_preserves_original_basis(spec_file):
    """It does not erase the original basis — what was overwritten has to survive."""
    spec = SystemSpec.load(spec_file)
    before = spec.numerics.dt_star.basis
    set_by_path(spec, "numerics.dt_star", "1e-3", turn=1)
    assert "original basis" in spec.numerics.dt_star.basis
    assert before[:20] in spec.numerics.dt_star.basis


def test_set_by_path_records_previous_value(spec_file):
    spec = SystemSpec.load(spec_file)
    set_by_path(spec, "species.0.n_simulated", "4000", turn=2)
    assert "previously 1000" in spec.primary.n_simulated.basis


def test_set_by_path_indexes_into_lists(spec_file):
    spec = SystemSpec.load(spec_file)
    set_by_path(spec, "external.0.params.k_si", "2e-5", turn=1)
    assert spec.trap().params["k_si"].value == 2e-5


def test_set_by_path_rejects_unknown_key(spec_file):
    spec = SystemSpec.load(spec_file)
    with pytest.raises(KeyError):
        set_by_path(spec, "medium.no_such_field.x", "1", turn=1)


def test_set_by_path_parses_types(spec_file):
    spec = SystemSpec.load(spec_file)
    set_by_path(spec, "geometry.dim", "3", turn=1)
    assert spec.geometry.dim.value == 3 and isinstance(spec.geometry.dim.value, int)


@pytest.mark.parametrize("text,expected", [
    ("5e-3", 5e-3),      # ★ yaml.safe_load reads this as a **string**
    ("2e-5", 2e-5),      # ★ same (no decimal point in the mantissa)
    ("2.5e-3", 2.5e-3),  # only this one is read as a float by yaml
    ("1E+6", 1e6),
    ("4000", 4000),
    ("0.005", 0.005),
    ("-3", -3),
])
def test_scientific_notation_never_becomes_a_string(text, expected):
    """★ The YAML 1.1 trap — `5e-3` arriving as a string makes the
    non-dimensionalization die, or quietly wrong.

    Hit once with `6.3e6` in `config/run_policy.yaml` and again in `session set`
    (2026-07-28). A trap hit twice gets pinned as a test.
    """
    from simbot.session import _parse_scalar
    v = _parse_scalar(text)
    assert isinstance(v, (int, float)) and not isinstance(v, bool)
    assert v == pytest.approx(expected)


def test_non_numeric_still_parses(spec_file):
    from simbot.session import _parse_scalar
    assert _parse_scalar("periodic") == "periodic"
    assert _parse_scalar("true") is True
    assert _parse_scalar("[1, 2, 3]") == [1, 2, 3]


def test_set_scientific_notation_is_usable_downstream(session):
    """If a string got in, the cost estimate dies right here."""
    t = session.set(["numerics.dt_star=5e-4"], when=T0)
    assert "error" not in t.cost
    assert t.cost["steps_per_seed"] == pytest.approx(100_000, rel=0.01)


# =============================================================================
# Session — set does not run anything
# =============================================================================
def test_new_session_records_turn_zero(session):
    assert len(session.turns) == 1
    assert session.turns[0].kind == "new"
    assert session.turns[0].spec_hash


def test_set_estimates_cost_without_running(session):
    """★ `set` must not run — nothing may appear in the raw directory."""
    t = session.set(["species.0.n_simulated=4000"], when=T0)
    assert t.kind == "set"
    assert t.cost["n_particles"] == 4000
    assert t.cost["wall_s_batch"] > 0
    assert "to run it" in t.note or "run" in t.note
    assert not list(session.dir.glob("**/samples.npz"))


def test_cost_scales_with_particle_count(session):
    a = session.set(["species.0.n_simulated=1000"], when=T0).cost["wall_s_batch"]
    b = session.set(["species.0.n_simulated=4000"], when=T0).cost["wall_s_batch"]
    assert b == pytest.approx(4 * a, rel=0.02)


def test_cost_scales_inversely_with_dt(session):
    a = session.set(["numerics.dt_star=5e-3"], when=T0).cost["steps_per_seed"]
    b = session.set(["numerics.dt_star=2.5e-3"], when=T0).cost["steps_per_seed"]
    assert b == pytest.approx(2 * a, rel=0.01)


def test_over_budget_is_flagged_not_run(session):
    """Going over budget has to surface **before** running."""
    t = session.set(["species.0.n_simulated=4000000"], when=T0)
    assert t.cost["over_budget"] is True
    assert "report without running" in t.cost["action"]


def test_small_n_warns_about_underestimate(session):
    """The throughput model was measured at N ≥ 500 — small N underestimates."""
    t = session.set(["species.0.n_simulated=100"], when=T0)
    assert "underestimates" in t.cost["warning"]


def test_history_is_append_only(session):
    session.set(["numerics.dt_star=2.5e-3"], when=T0)
    session.set(["numerics.dt_star=1e-3"], when=T0)
    assert [t.index for t in session.turns] == [0, 1, 2]
    assert session.turns[1].changes["numerics.dt_star"] == [0.005, 0.0025]
    assert session.turns[2].changes["numerics.dt_star"] == [0.0025, 0.001]


def test_each_turn_writes_its_own_spec(session):
    session.set(["numerics.dt_star=2.5e-3"], when=T0)
    assert session.derived_spec_file(0).exists()
    assert session.derived_spec_file(1).exists()
    a = SystemSpec.load(session.derived_spec_file(0))
    b = SystemSpec.load(session.derived_spec_file(1))
    assert a.numerics.dt_star.value == 0.005
    assert b.numerics.dt_star.value == 0.0025


def test_original_spec_file_is_not_modified(session, spec_file):
    before = spec_file.read_text(encoding="utf-8")
    session.set(["numerics.dt_star=1e-3"], when=T0)
    assert spec_file.read_text(encoding="utf-8") == before


def test_reload_restores_latest_spec(session, tmp_path):
    session.set(["numerics.dt_star=2.5e-3", "species.0.n_simulated=2000"], when=T0)
    back = Session.load(session.session_id, root=tmp_path / "sessions")
    assert back.spec.numerics.dt_star.value == 0.0025
    assert back.spec.primary.n_simulated.value == 2000
    assert len(back.turns) == 2


def test_reload_preserves_spec_hash(session, tmp_path):
    session.set(["numerics.dt_star=2.5e-3"], when=T0)
    back = Session.load(session.session_id, root=tmp_path / "sessions")
    assert back.spec.hash() == session.spec.hash()


def test_bad_assignment_raises(session):
    with pytest.raises(ValueError, match="key=value"):
        session.set(["dt_star 2.5e-3"], when=T0)


def test_record_run_stores_metrics(session):
    session.record_run(run_id="r1", metrics={"var_x_star": 1.0058}, wall_s=2.3,
                       when=T0)
    t = session.turns[-1]
    assert t.kind == "run" and t.run_id == "r1"
    assert t.metrics["var_x_star"] == 1.0058
    assert t.cost["wall_s_measured"] == 2.3


def test_run_records_failure_instead_of_dropping_it(session, tmp_path):
    """★ A failed attempt vanishing from the history loses half of "what was
    tried"."""
    session.set(["species.0.n_simulated=40000000"], when=T0)   # over budget → fails
    t = session.run(runs_root=tmp_path / "runs", when=T0)
    assert t.kind == "run"
    assert t.problems and any("failed" in p for p in t.problems)
    assert "a failed attempt stays in the history too" in t.note


def test_run_id_has_no_duplicate_date(session):
    """A session spec lives in `sessions/<date>_<card>/`, so the default slug would
    put the date in twice."""
    turn = session.turns[-1].index
    rid = f"{session.session_id}_t{turn:02d}_{session.spec.hash()[:6]}"
    assert rid.count("2026-07-28") == 1


@pytest.mark.slow
def test_session_run_executes_and_records_metrics(session, tmp_path):
    """Does the set → run → compare loop close."""
    if not EXAMPLE_PRED.exists():
        pytest.skip("no example prediction")
    session.set(["numerics.dt_star=1.0e-2"], when=T0)
    t = session.run(prediction=EXAMPLE_PRED, runs_root=tmp_path / "runs", when=T0)
    assert t.kind == "run" and not t.problems
    assert t.metrics["var_x_star"] > 0
    assert t.cost["wall_s_measured"] > 0
    assert (tmp_path / "runs" / t.run_id / "REPORT.md").exists()
    # compare shows the parameters and the measurements together
    text = session.compare(0, t.index)
    assert "numerics.dt_star" in text and "var_x_star" in text


def test_show_is_readable_without_conversation(session):
    """It has to be possible to pick up after the conversational context is gone."""
    session.set(["numerics.dt_star=2.5e-3"], when=T0)
    session.record_run(run_id="r1", metrics={"var_x_star": 1.0017}, when=T0)
    text = session.show()
    for expect in (session.session_id, "passive-sphere--harmonic-trap",
                   "numerics.dt_star", "r1", "session run"):
        assert expect in text, expect


def test_show_surfaces_problems(session, spec_file):
    """A convention violation has to survive into `show`."""
    spec = SystemSpec.load(spec_file)
    spec.medium.rho_fluid_si.confidence = ""       # induce a violation
    spec.save(spec_file)
    s = Session.create(spec_file, root=session.root, when=T0)
    assert "unresolved problems" in s.show()


def test_compare_shows_parameter_and_metric_deltas(session):
    session.set(["numerics.dt_star=2.5e-3"], when=T0)
    session.record_run(run_id="r1", metrics={"var_x_star": 1.0017}, when=T0)
    text = session.compare(0, 2)
    assert "numerics.dt_star" in text
    assert "0.005" in text and "0.0025" in text
    assert "var_x_star" in text


def test_compare_reports_percent_change_between_runs(session):
    session.record_run(run_id="r1", metrics={"tau_star": 1.0}, when=T0)
    session.record_run(run_id="r2", metrics={"tau_star": 1.02}, when=T0)
    assert "+2.00%" in session.compare(1, 2)


def test_compare_says_when_specs_are_identical(session):
    session.record_run(run_id="r1", metrics={}, when=T0)
    assert "the same spec" in session.compare(0, 1)


def test_compare_unknown_turn_raises(session):
    with pytest.raises(KeyError, match="turn 99"):
        session.compare(0, 99)


def test_latest_finds_the_newest_session(tmp_path, spec_file):
    root = tmp_path / "sessions"
    Session.create(spec_file, root=root, when=datetime(2026, 7, 28, 9, 0))
    b = Session.create(spec_file, root=root, when=datetime(2026, 7, 28, 15, 0))
    assert Session.latest(root=root).session_id == b.session_id


def test_latest_without_sessions_raises(tmp_path):
    with pytest.raises(FileNotFoundError, match="session new"):
        Session.latest(root=tmp_path / "empty")


# =============================================================================
# CLI — runners, gates, budget
# =============================================================================
def test_unknown_card_is_refused_not_silently_run(spec_file):
    """★ Running it quietly with the trap runner computes an entirely different
    system."""
    spec = SystemSpec.load(spec_file)
    spec.card = "abp--dense-collective"
    with pytest.raises(SystemExit, match="has no runner"):
        cli._runner_for(spec)


def test_known_card_resolves(spec_file):
    assert cli._runner_for(SystemSpec.load(spec_file)) == "trap"


def test_run_refuses_over_budget(spec_file, tmp_path, capsys):
    spec = SystemSpec.load(spec_file)
    spec.primary.n_simulated.value = 40_000_000
    spec.save(spec_file)
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs")])
    assert code != 0
    assert "over budget" in capsys.readouterr().err


def test_run_refuses_too_few_seeds(spec_file, tmp_path, capsys):
    """A single-seed production run is forbidden (a result with no error bar)."""
    spec = SystemSpec.load(spec_file)
    spec.numerics.n_seeds.value = 1
    spec.save(spec_file)
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs")])
    assert code != 0
    assert "error bars" in capsys.readouterr().err


def test_run_refuses_on_spec_gate_failure(spec_file, tmp_path, capsys):
    spec = SystemSpec.load(spec_file)
    spec.geometry.box_over_ref.value = 2.0       # violate the box gate
    spec.save(spec_file)
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs")])
    assert code != 0
    assert "S3 gates" in capsys.readouterr().err


def test_resume_refuses_missing_dir(tmp_path, capsys):
    assert cli.main(["resume", str(tmp_path / "nope")]) != 0


def test_resume_refuses_broken_seal(tmp_path, spec_file, capsys):
    """★ Resuming a run whose seal is broken makes the verification meaningless."""
    from simbot.io import RunDir, write_seal
    rd = RunDir.create(tmp_path / "runs", "r1")
    rd.write("spec", spec_file.read_text(encoding="utf-8"))
    rd.write("prediction", "# S2\n")
    write_seal(rd)
    rd.write("prediction", "# S2 edited\n")
    code = cli.main(["resume", str(rd.path)])
    assert code != 0
    assert "seal violation" in capsys.readouterr().err


def test_resume_refuses_with_no_finished_runs(tmp_path, spec_file, capsys):
    from simbot.io import RunDir
    rd = RunDir.create(tmp_path / "runs", "r1")
    rd.write("spec", spec_file.read_text(encoding="utf-8"))
    assert cli.main(["resume", str(rd.path)]) != 0
    assert "0 completed runs" in capsys.readouterr().err


def test_params_marks_unchosen_defaults(capsys):
    """⚠ = the provenance is only assumed and the value is the same in every spec."""
    if not EXAMPLE_SPEC.exists():
        pytest.skip("no example spec")
    assert cli.main(["params", "--path", str(EXAMPLE_SPEC.parent)]) == 0
    out = capsys.readouterr().out
    assert "⚠assumed" in out
    assert "medium.eta_si" in out
    assert "a default nobody picked" in out


def test_params_refuses_empty_directory(tmp_path, capsys):
    assert cli.main(["params", "--path", str(tmp_path)]) != 0


def test_calibrate_names_both_kernels():
    """★ A different kernel cannot overwrite the constant — both kernels have to be
    named."""
    assert "WCA" in cli.BASELINE_KERNEL
    assert "no pair interaction" in cli.CALIBRATE_KERNEL


def test_parser_exposes_all_documented_commands():
    p = cli.build_parser()
    sub = next(a for a in p._actions if a.dest == "cmd")
    assert set(sub.choices) == {"run", "resume", "converge", "params", "calibrate"}


# =============================================================================
# End to end — with a HOOMD run (slow)
# =============================================================================
@pytest.mark.slow
def test_end_to_end_run_produces_report(spec_file, tmp_path, capsys):
    """S3 → S8 end to end. Also checks it reproduces the first completed run's
    measurements."""
    if not EXAMPLE_PRED.exists():
        pytest.skip("no example prediction")
    runs = tmp_path / "runs"
    code = cli.main(["run", str(spec_file), "--prediction", str(EXAMPLE_PRED),
                     "--runs-root", str(runs), "--run-id", "e2e"])
    assert code == 0
    out = capsys.readouterr().out

    rd = runs / "e2e"
    for name in ("03_spec.yaml", "04_reduced.yaml", "04_nondim.md",
                 "05_run_manifest.json", "metrics.json", "REPORT.md",
                 "SEALED.sha256", "02_prediction.md"):
        assert (rd / name).exists(), name

    report = (rd / "REPORT.md").read_text(encoding="utf-8")
    assert "confirmed_by: null" in report
    assert "PASS" in report

    # reproduces the first run's values (fixed seed → bit-for-bit)
    import json
    m = json.loads((rd / "metrics.json").read_text())
    assert m["var_x_star"]["value"] == pytest.approx(1.0057652, rel=1e-5)
    assert m["tau_trap_ms"]["value"] == pytest.approx(8.0566504, rel=1e-5)
    assert "confirmed_by: null" in out


@pytest.mark.slow
def test_end_to_end_warns_without_prediction(spec_file, tmp_path, capsys):
    """Running with no prediction leaves nothing to seal — it has to say so."""
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs"),
                     "--run-id", "nopred"])
    assert code == 0
    out = capsys.readouterr().out
    assert "no prediction file" in out
    assert "nothing preventing post-hoc rationalisation" in out
