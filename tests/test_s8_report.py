"""S8 — `REPORT.md` 생성.

리포트는 요약이 아니라 **감사 기록**이다. 그래서 이 파일이 확인하는 것은
"예쁘게 나오는가"가 아니라 **나쁜 소식이 빠지지 않는가**다:

- 봉인이 깨졌으면 대조표를 싣지 않고 맨 위에 경고가 나오는가
- `INCONCLUSIVE` 와 그 이유가 남는가
- `git_dirty` 가 보고되는가
- 캡션 없는 그림이 표시되는가
- `confirmed_by: null` 이 있는가
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
        pytest.skip("examples/trap-2d-5um/spec.yaml 없음")
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
        PredictionItem("var_x_2d_nm2", 414.19, "±1.5%", "⟨x²⟩ = kT/k (정확해)"),
        PredictionItem("em_bias_reproduced", 1.0025063, "±0.1%", "EM 편향 재현",
                       competing_value=1.0, discriminates="적분기 스킴"),
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
# 전체 렌더
# =============================================================================
def test_report_renders_and_writes(rundir, full_inputs):
    p = R.write_report(rundir, full_inputs)
    assert p.exists()
    text = p.read_text(encoding="utf-8")
    assert text.startswith("# REPORT — `2026-07-28_trap-2d-5um_c59e93`")


def test_report_never_confirms(rundir, full_inputs):
    """★ 리포트에 사람 확정 도장이 찍히면 안 된다."""
    text = R.render(rundir, full_inputs)
    assert "confirmed_by: null" in text
    assert "confirmed_by: agent" not in text


def test_report_states_verdict_is_a_proposal(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "판정은 **제안**이다" in text


def test_report_is_self_contained_on_all_stages(rundir, full_inputs):
    """§1~§6 이 전부 있어야 사람이 이 파일 하나로 판단할 수 있다."""
    text = R.render(rundir, full_inputs)
    for head in ("## 1. 판정 요약", "## 2. 시스템 명세와 게이트 (S3)",
                 "## 3. 무차원화 (S4)", "## 4. 그림 (S6)",
                 "## 5. 재현 가능성", "## 6. 에이전트가 쓴 문서"):
        assert head in text, head


def test_report_survives_missing_pieces(rundir):
    """입력이 비어도 렌더는 되고, 없는 것은 '없음'으로 나온다."""
    text = R.render(rundir, R.ReportInputs())
    assert "판정 없음" in text
    assert "_S7 판정 없음._" in text
    assert "_S3 검사 결과 없음._" in text


# =============================================================================
# 봉인 — 깨졌으면 대조표를 싣지 않는다
# =============================================================================
def test_intact_seal_shows_external_verification_command(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "shasum -a 256 -c SEALED.sha256" in text


def test_broken_seal_suppresses_the_comparison_table(rundir, full_inputs):
    """★ 봉인이 깨진 상태의 대조표는 검증처럼 보이지만 검증이 아니다."""
    rundir.write("prediction", "# S2 (결과 보고 나서 고침)\n")
    text = R.render(rundir, full_inputs)
    assert "⛔ seal violation" in text
    assert "대조표는 **생성하지 않았다**" in text
    assert "| 양 | 예측 (봉인) |" not in text


def test_broken_seal_warning_appears_before_results(rundir, full_inputs):
    rundir.write("prediction", "# 고침\n")
    text = R.render(rundir, full_inputs)
    assert text.index("⛔ seal violation") < text.index("## 1. 판정 요약")


# =============================================================================
# INCONCLUSIVE 가 사라지지 않는다
# =============================================================================
def test_inconclusive_and_its_reason_are_reported(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "INCONCLUSIVE" in text
    assert "설계 검정력" in text
    assert "표본 28.5배 필요" in text or "28.5" in text


def test_inconclusive_is_framed_as_a_fact_not_a_failure(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "실패가 아니다" in text


def test_counts_are_reported(rundir, full_inputs):
    assert "2개 중 1 PASS · 1 INCONCLUSIVE · 0 FAIL." in R.render(rundir, full_inputs)


def test_validation_problems_are_surfaced(rundir, spec):
    """오차 막대 없는 측정이 리포트에서 드러나야 한다."""
    meas = {"var_x_2d_nm2": Measurement("var_x_2d_nm2", 416.6),   # stat_err 없음
            "em_bias_reproduced": Measurement("em_bias_reproduced", 1.0057, stat_err=0.004)}
    rep = validate_run(_prediction(), meas, rundir=rundir)
    text = R.render(rundir, R.ReportInputs(spec=spec, validation=rep))
    assert "검증 절차의 문제" in text
    assert "통계오차가 없다" in text


# =============================================================================
# 재현 가능성
# =============================================================================
def test_dirty_git_is_reported_as_a_reproducibility_limit(rundir, spec):
    man = {"run_id": "r", "git_rev": "abc1234", "git_dirty": True, "env": {}}
    text = R.reproducibility_section(man)
    assert "커밋되지 않은 변경 있음" in text
    assert "재현되지 않는다" in text


def test_clean_git_is_reported_as_clean():
    text = R.reproducibility_section({"git_rev": "abc", "git_dirty": False, "env": {}})
    assert "clean" in text


def test_unknown_git_state_is_not_claimed_clean():
    """판정 불가를 clean 으로 적으면 재현 가능성을 거짓 주장한다."""
    text = R.reproducibility_section({"git_rev": "?", "git_dirty": None, "env": {}})
    assert "판정 불가" in text
    assert "clean" not in text


def test_missing_manifest_refuses_to_claim_reproducibility():
    assert "주장할 수 없다" in R.reproducibility_section(None)


def test_manifest_records_versions(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "hoomd" in text and "7.1.0" in text


# =============================================================================
# 게이트
# =============================================================================
def test_deferred_gates_are_listed(rundir, full_inputs):
    """S3 에서 판정할 수 없는 게이트가 몇 개인지 사람이 알아야 한다."""
    text = R.render(rundir, full_inputs)
    assert "S7 이 판정해야 하는 게이트" in text
    assert "`equipartition`" in text


def test_off_gates_show_their_reason(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "쌍 상호작용이 없어 겹침이 발생할 수 없다" in text


def test_spec_problems_are_surfaced(rundir, spec):
    from simbot.spec import Q
    spec.medium.rho_fluid_si = Q(996.5, "kg/m^3", "assumed", "근거는 있으나 신뢰도 없음")
    text = R.gates_section(R.validate_spec(spec))
    assert "규약 위반" in text and "confidence" in text


# =============================================================================
# 무차원화 절
# =============================================================================
def test_nondim_section_reports_roundtrip_gate(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "왕복 오차" in text
    assert "✅ 통과" in text


def test_nondim_section_separates_logged_from_gated(rundir, full_inputs):
    """`dt/τ_D` 는 기록이지 게이트가 아니라는 것이 리포트에 드러나야 한다."""
    text = R.render(rundir, full_inputs)
    assert "기록용 (게이트 아님" in text
    assert "dt_over_tau_D" in text


def test_nondim_section_names_the_scale_origin(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "scales_harmonic_trap: (l_trap, kT, tau_trap)" in text


def test_nondim_section_empty_without_spec(rundir):
    assert "없음" in R.nondim_section(None, None)


def test_unregistered_card_does_not_crash_the_report(rundir, spec):
    """척도 규칙이 없는 카드여도 리포트는 나와야 한다 (그 절만 빈다)."""
    spec.card = "brand--new-pair"
    text = R.render(rundir, R.ReportInputs(spec=spec))
    assert "REPORT" in text
    assert "_S4 무차원화 결과 없음._" in text


# =============================================================================
# 그림 — 캡션 없는 그림은 산출물이 아니다
# =============================================================================
def test_figure_without_caption_is_flagged(rundir, full_inputs):
    (rundir.figs / "01_msd.png").write_bytes(b"\x89PNG\r\n")
    text = R.render(rundir, full_inputs)
    assert "캡션 없음" in text
    assert "§S6 게이트" in text


def test_figure_with_caption_is_embedded(rundir, full_inputs):
    (rundir.figs / "01_msd.png").write_bytes(b"\x89PNG\r\n")
    full_inputs.figures = {"01_msd.png": "MSD 와 해석해 2d(1−e^{−t/τ})"}
    text = R.render(rundir, full_inputs)
    assert "![MSD 와 해석해" in text
    assert "캡션 없음" not in text


def test_no_figures_says_so(rundir, full_inputs):
    assert "_그림 없음._" in R.render(rundir, full_inputs)


# =============================================================================
# 비용
# =============================================================================
def test_cost_section_reports_per_run_average():
    assert "런당 평균 `0.66 s`" in R.cost_section(10.6, 16)


def test_missing_cost_says_so():
    assert "기록 없음" in R.cost_section(None, None)


# =============================================================================
# 에이전트 문서 인용 — 대신 쓰지 않는다
# =============================================================================
def test_missing_agent_document_is_named_not_invented(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "`07_validation.md` 가 없다" in text
    assert "에이전트가 써야 한다" in text


def test_existing_agent_documents_are_linked(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "[`08_conclusion.md`](08_conclusion.md)" in text


def test_report_states_the_division_of_labour(rundir, full_inputs):
    text = R.render(rundir, full_inputs)
    assert "대신 쓰지 않는다" in text
