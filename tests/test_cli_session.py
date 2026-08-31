"""CLI · 세션 층.

## 이 파일이 지키는 것

1. **`set` 은 실행하지 않는다.** 비용 추정만 한다. 이 분리가 깨지면
   "N 을 8000 으로 해보자"가 11분짜리 런을 조용히 시작한다.
2. **이력은 append-only.** 과거 턴을 덮어쓰면 무엇을 시도했는지가 사라진다.
3. **러너 없는 카드를 억지로 돌리지 않는다.** 조용히 트랩 러너로 돌리면
   전혀 다른 계를 계산하고 그 사실을 아무도 모른다.
4. **예산·게이트를 넘기면 실행 전에 멈춘다.**

HOOMD 실행이 필요한 테스트는 `slow` 마커를 단다.
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
        pytest.skip("examples/trap-2d-5um/spec.yaml 없음")
    p = tmp_path / "spec.yaml"
    p.write_text(EXAMPLE_SPEC.read_text(encoding="utf-8"), encoding="utf-8")
    return p


@pytest.fixture
def session(spec_file, tmp_path) -> Session:
    return Session.create(spec_file, root=tmp_path / "sessions", when=T0)


# =============================================================================
# dotted-key 설정
# =============================================================================
def test_set_by_path_updates_quantity_value(spec_file):
    spec = SystemSpec.load(spec_file)
    old, new = set_by_path(spec, "numerics.dt_star", "2.5e-3", turn=1)
    assert old == 0.005 and new == 0.0025
    assert spec.numerics.dt_star.value == 0.0025


def test_set_by_path_marks_provenance_user(spec_file):
    """★ 사람이 정한 값은 가정도 정책도 아니다 — S7b 가 흔들 대상이 아니다."""
    spec = SystemSpec.load(spec_file)
    set_by_path(spec, "medium.eta_si", "1.002e-3", turn=3)
    q = spec.medium.eta_si
    assert q.provenance == "user"
    assert "턴 3 에서 사람이 지정" in q.basis


def test_set_by_path_preserves_original_basis(spec_file):
    """원래 근거를 지우지 않는다 — 무엇을 덮어썼는지 남아야 한다."""
    spec = SystemSpec.load(spec_file)
    before = spec.numerics.dt_star.basis
    set_by_path(spec, "numerics.dt_star", "1e-3", turn=1)
    assert "원래 근거" in spec.numerics.dt_star.basis
    assert before[:20] in spec.numerics.dt_star.basis


def test_set_by_path_records_previous_value(spec_file):
    spec = SystemSpec.load(spec_file)
    set_by_path(spec, "species.0.n_simulated", "4000", turn=2)
    assert "이전 1000" in spec.primary.n_simulated.basis


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
    ("5e-3", 5e-3),      # ★ yaml.safe_load 는 이것을 **문자열**로 읽는다
    ("2e-5", 2e-5),      # ★ 같음 (만티사에 소수점이 없다)
    ("2.5e-3", 2.5e-3),  # 이것만 yaml 이 float 로 읽는다
    ("1E+6", 1e6),
    ("4000", 4000),
    ("0.005", 0.005),
    ("-3", -3),
])
def test_scientific_notation_never_becomes_a_string(text, expected):
    """★ YAML 1.1 함정 — `5e-3` 이 문자열로 들어가면 무차원화가 죽거나 조용히 틀린다.

    `config/run_policy.yaml` 의 `6.3e6` 에서 한 번, `session set` 에서 또 한 번
    겪었다 (2026-07-28). 두 번 겪은 함정은 테스트로 고정한다.
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
    """문자열이 들어갔다면 여기서 비용 추정이 죽는다."""
    t = session.set(["numerics.dt_star=5e-4"], when=T0)
    assert "error" not in t.cost
    assert t.cost["steps_per_seed"] == pytest.approx(100_000, rel=0.01)


# =============================================================================
# 세션 — set 은 실행하지 않는다
# =============================================================================
def test_new_session_records_turn_zero(session):
    assert len(session.turns) == 1
    assert session.turns[0].kind == "new"
    assert session.turns[0].spec_hash


def test_set_estimates_cost_without_running(session):
    """★ `set` 이 실행하면 안 된다 — raw 디렉터리에 아무것도 생기지 않아야 한다."""
    t = session.set(["species.0.n_simulated=4000"], when=T0)
    assert t.kind == "set"
    assert t.cost["n_particles"] == 4000
    assert t.cost["wall_s_batch"] > 0
    assert "실행하려면" in t.note or "실행" in t.note
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
    """예산 초과는 **실행 전에** 드러나야 한다."""
    t = session.set(["species.0.n_simulated=4000000"], when=T0)
    assert t.cost["over_budget"] is True
    assert "실행하지 않고" in t.cost["action"]


def test_small_n_warns_about_underestimate(session):
    """처리량 모델은 N ≥ 500 에서 실측됐다 — 작은 N 은 과소추정이다."""
    t = session.set(["species.0.n_simulated=100"], when=T0)
    assert "과소추정" in t.cost["warning"]


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
    with pytest.raises(ValueError, match="키=값"):
        session.set(["dt_star 2.5e-3"], when=T0)


def test_record_run_stores_metrics(session):
    session.record_run(run_id="r1", metrics={"var_x_star": 1.0058}, wall_s=2.3,
                       when=T0)
    t = session.turns[-1]
    assert t.kind == "run" and t.run_id == "r1"
    assert t.metrics["var_x_star"] == 1.0058
    assert t.cost["wall_s_measured"] == 2.3


def test_run_records_failure_instead_of_dropping_it(session, tmp_path):
    """★ 실패한 시도가 이력에서 사라지면 '무엇을 시도했는가'의 절반이 없어진다."""
    session.set(["species.0.n_simulated=40000000"], when=T0)   # 예산 초과 → 실패
    t = session.run(runs_root=tmp_path / "runs", when=T0)
    assert t.kind == "run"
    assert t.problems and any("실패" in p for p in t.problems)
    assert "실패한 시도도 이력에 남는다" in t.note


def test_run_id_has_no_duplicate_date(session):
    """세션 spec 은 `sessions/<날짜>_<카드>/` 에 있어 기본 슬러그가 날짜를 두 번 넣는다."""
    turn = session.turns[-1].index
    rid = f"{session.session_id}_t{turn:02d}_{session.spec.hash()[:6]}"
    assert rid.count("2026-07-28") == 1


@pytest.mark.slow
def test_session_run_executes_and_records_metrics(session, tmp_path):
    """set → run → compare 루프가 닫히는가."""
    if not EXAMPLE_PRED.exists():
        pytest.skip("예시 예측 없음")
    session.set(["numerics.dt_star=1.0e-2"], when=T0)
    t = session.run(prediction=EXAMPLE_PRED, runs_root=tmp_path / "runs", when=T0)
    assert t.kind == "run" and not t.problems
    assert t.metrics["var_x_star"] > 0
    assert t.cost["wall_s_measured"] > 0
    assert (tmp_path / "runs" / t.run_id / "REPORT.md").exists()
    # compare 가 파라미터와 측정값을 함께 보여준다
    text = session.compare(0, t.index)
    assert "numerics.dt_star" in text and "var_x_star" in text


def test_show_is_readable_without_conversation(session):
    """대화 컨텍스트가 날아가도 이어받을 수 있어야 한다."""
    session.set(["numerics.dt_star=2.5e-3"], when=T0)
    session.record_run(run_id="r1", metrics={"var_x_star": 1.0017}, when=T0)
    text = session.show()
    for expect in (session.session_id, "passive-sphere--harmonic-trap",
                   "numerics.dt_star", "r1", "session run"):
        assert expect in text, expect


def test_show_surfaces_problems(session, spec_file):
    """규약 위반이 있으면 show 에 남아야 한다."""
    spec = SystemSpec.load(spec_file)
    spec.medium.rho_fluid_si.confidence = ""       # 위반 유발
    spec.save(spec_file)
    s = Session.create(spec_file, root=session.root, when=T0)
    assert "미해결 문제" in s.show()


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
    assert "같은 spec" in session.compare(0, 1)


def test_compare_unknown_turn_raises(session):
    with pytest.raises(KeyError, match="턴 99"):
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
# CLI — 러너·게이트·예산
# =============================================================================
def test_unknown_card_is_refused_not_silently_run(spec_file):
    """★ 조용히 트랩 러너로 돌리면 전혀 다른 계를 계산한다."""
    spec = SystemSpec.load(spec_file)
    spec.card = "abp--dense-collective"
    with pytest.raises(SystemExit, match="러너가 없다"):
        cli._runner_for(spec)


def test_known_card_resolves(spec_file):
    assert cli._runner_for(SystemSpec.load(spec_file)) == "trap"


def test_run_refuses_over_budget(spec_file, tmp_path, capsys):
    spec = SystemSpec.load(spec_file)
    spec.primary.n_simulated.value = 40_000_000
    spec.save(spec_file)
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs")])
    assert code != 0
    assert "예산 초과" in capsys.readouterr().err


def test_run_refuses_too_few_seeds(spec_file, tmp_path, capsys):
    """시드 1개짜리 프로덕션 런은 금지 (오차 막대 없는 결과)."""
    spec = SystemSpec.load(spec_file)
    spec.numerics.n_seeds.value = 1
    spec.save(spec_file)
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs")])
    assert code != 0
    assert "오차 막대" in capsys.readouterr().err


def test_run_refuses_on_spec_gate_failure(spec_file, tmp_path, capsys):
    spec = SystemSpec.load(spec_file)
    spec.geometry.box_over_ref.value = 2.0       # 박스 게이트 위반
    spec.save(spec_file)
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs")])
    assert code != 0
    assert "S3 게이트" in capsys.readouterr().err


def test_resume_refuses_missing_dir(tmp_path, capsys):
    assert cli.main(["resume", str(tmp_path / "nope")]) != 0


def test_resume_refuses_broken_seal(tmp_path, spec_file, capsys):
    """★ 봉인이 깨진 런을 이어받으면 검증이 무의미해진다."""
    from simbot.io import RunDir, write_seal
    rd = RunDir.create(tmp_path / "runs", "r1")
    rd.write("spec", spec_file.read_text(encoding="utf-8"))
    rd.write("prediction", "# S2\n")
    write_seal(rd)
    rd.write("prediction", "# S2 고침\n")
    code = cli.main(["resume", str(rd.path)])
    assert code != 0
    assert "seal violation" in capsys.readouterr().err


def test_resume_refuses_with_no_finished_runs(tmp_path, spec_file, capsys):
    from simbot.io import RunDir
    rd = RunDir.create(tmp_path / "runs", "r1")
    rd.write("spec", spec_file.read_text(encoding="utf-8"))
    assert cli.main(["resume", str(rd.path)]) != 0
    assert "완주한 런이 0개" in capsys.readouterr().err


def test_params_marks_unchosen_defaults(capsys):
    """⚠ = provenance 가 assumed 뿐이고 모든 spec 에서 같은 값."""
    if not EXAMPLE_SPEC.exists():
        pytest.skip("예시 spec 없음")
    assert cli.main(["params", "--path", str(EXAMPLE_SPEC.parent)]) == 0
    out = capsys.readouterr().out
    assert "⚠assumed" in out
    assert "medium.eta_si" in out
    assert "아무도 고르지 않은 기본값" in out


def test_params_refuses_empty_directory(tmp_path, capsys):
    assert cli.main(["params", "--path", str(tmp_path)]) != 0


def test_calibrate_names_both_kernels():
    """★ 커널이 다르면 상수를 덮어쓸 수 없다 — 두 커널이 다 적혀야 한다."""
    assert "WCA" in cli.BASELINE_KERNEL
    assert "쌍 상호작용 없음" in cli.CALIBRATE_KERNEL


def test_parser_exposes_all_documented_commands():
    p = cli.build_parser()
    sub = next(a for a in p._actions if a.dest == "cmd")
    assert set(sub.choices) == {"run", "resume", "converge", "params", "calibrate"}


# =============================================================================
# 전체 관통 — HOOMD 실행 (slow)
# =============================================================================
@pytest.mark.slow
def test_end_to_end_run_produces_report(spec_file, tmp_path, capsys):
    """S3 → S8 관통. 첫 완주의 측정값을 재현하는지도 함께 본다."""
    if not EXAMPLE_PRED.exists():
        pytest.skip("예시 예측 없음")
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

    # 첫 완주 값 재현 (seed 고정 → 비트 단위 재현)
    import json
    m = json.loads((rd / "metrics.json").read_text())
    assert m["var_x_star"]["value"] == pytest.approx(1.0057652, rel=1e-5)
    assert m["tau_trap_ms"]["value"] == pytest.approx(8.0566504, rel=1e-5)
    assert "confirmed_by: null" in out


@pytest.mark.slow
def test_end_to_end_warns_without_prediction(spec_file, tmp_path, capsys):
    """예측 없이 돌리면 봉인할 것이 없다 — 그 사실을 말해야 한다."""
    code = cli.main(["run", str(spec_file), "--runs-root", str(tmp_path / "runs"),
                     "--run-id", "nopred"])
    assert code == 0
    out = capsys.readouterr().out
    assert "예측 파일이 없다" in out
    assert "사후합리화를 막을 장치가 없다" in out
