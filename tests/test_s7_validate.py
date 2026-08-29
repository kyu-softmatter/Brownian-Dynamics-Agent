"""S7 — 판정 제안.

## 이 파일의 중심 테스트

`test_reproduces_first_run_verdicts` — 첫 완주의 **9개 판정을 그대로 재현한다**
(7 PASS · 2 INCONCLUSIVE · 0 FAIL). 사람이 손으로 내린 판정과 코드의 판정이
갈리면 둘 중 하나가 틀렸고, 어느 쪽이든 알아야 한다.

입력 숫자는 `runs/2026-07-28_trap-2d-5um_2dfb9d/metrics.json` 에서 온다.
`runs/` 는 gitignore 대상이므로 그 파일에 의존하지 않고 **값을 여기 박아둔다** —
회귀 테스트는 체크아웃 어디서나 돌아야 한다.
"""
from __future__ import annotations

import math

import pytest

from simbot import validate as V
from simbot.io import RunDir, write_seal
from simbot.spec import Prediction, PredictionItem
from simbot.validate import FAIL, INCONCLUSIVE, PASS, Measurement, compare

# =============================================================================
# tolerance 파싱
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
    """★ 읽지 못한 tolerance 를 조용히 넘기면 그 항목은 판정 없이 통과한다."""
    with pytest.raises(ValueError, match="읽을 수 없다"):
        V.parse_tolerance("대충 맞으면 됨")


def test_zero_tolerance_raises():
    """반폭 0 은 어떤 결과든 FAIL 로 만든다 — 검증이 아니다."""
    with pytest.raises(ValueError, match="0 이하"):
        V.parse_tolerance("±0")


def test_relative_half_width_scales_with_prediction():
    t = V.parse_tolerance("±1.5%")
    assert t.half_width(414.19) == pytest.approx(414.19 * 0.015)
    assert t.half_width(None) is None


# =============================================================================
# 측정값 — 오차 없는 수는 결론에 못 쓴다
# =============================================================================
def test_measurement_without_stat_err_is_a_problem():
    assert any("통계오차가 없다" in p
               for p in Measurement("D", 1.0).problems())


def test_zero_stat_err_flags_identity_risk():
    """요동하지 않는 '측정값'은 산술 항등식이다 (2026-07-28 실제 사례)."""
    assert any("항등식" in p for p in Measurement("x", 1.0, stat_err=0.0).problems())


def test_nonfinite_measurement_is_a_problem():
    assert any("비유한값" in p
               for p in Measurement("x", math.nan, stat_err=0.1).problems())


def test_verdict_is_inconclusive_without_error_bar():
    item = PredictionItem("D", 1.0, "±3%", "Stokes-Einstein")
    row = compare(item, Measurement("D", 1.0))
    assert row.verdict == INCONCLUSIVE and "오차 막대" in row.reason


# =============================================================================
# 판정 논리
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
    """★ SE 가 대역보다 크면 이 측정은 이 tolerance 를 판정할 수 없다.

    이걸 PASS 로 쓰면 검증이 아니라 요행이다 — 대역이 넓어서가 아니라
    **측정이 그 대역을 해상하지 못해서** 우연히 안에 들어온 것이다.
    """
    item = PredictionItem("var_x_star", 1.00251, "±0.1%", "EM 편향")
    row = compare(item, Measurement("var_x_star", 1.00577, stat_err=0.00446))
    assert row.verdict == INCONCLUSIVE
    assert "판정할 수 없다" in row.reason
    assert "필요 표본 배수" in row.reason


def test_required_sample_multiple_is_the_variance_ratio():
    item = PredictionItem("x", 1.0, "±0.1%", "b")
    row = compare(item, Measurement("x", 1.0, stat_err=0.01))
    # SE/half = 0.01/0.001 = 10 → 100배
    assert "100" in row.reason


def test_inside_band_but_low_power_is_inconclusive():
    """★ CLAUDE.md 통계 4규칙 — 검정력이 3σ 를 못 만드는 곳에서 3σ 를 요구하지 않는다.

    대역 안에 들어왔더라도 경쟁 가설과 구별되지 않으면 그 항목은 아무것도
    판별하지 못했다. `INCONCLUSIVE` 를 사실로 고정한다.
    """
    # 예측값은 해석식 1/(1−dt*/2) 의 전정밀도 값을 쓴다 — 절단값을 쓰면
    # 검정력이 0.5628 로 나오고 문서의 0.5618 과 어긋난다
    item = PredictionItem("var_x_star", 1.0025063, "±1%", "EM 편향",
                          competing_value=1.0)      # exact 스킴
    row = compare(item, Measurement("var_x_star", 1.0057652, stat_err=0.0044614))
    assert row.verdict == INCONCLUSIVE
    assert row.design_power == pytest.approx(0.5618, rel=1e-3)
    assert "예견된 한계이지 실패가 아니다" in row.reason


def test_high_power_gives_pass():
    """dt* 를 키우면 판별력이 생긴다 — 첫 런의 `dt*=2e-2` 상황."""
    item = PredictionItem("var_x_star", 1.01010, "±1%", "EM 편향",
                          competing_value=1.0)
    row = compare(item, Measurement("var_x_star", 1.01041, stat_err=0.00264))
    assert row.verdict == PASS
    assert row.design_power > 3.0


def test_samples_needed_for_3sigma_is_reported():
    item = PredictionItem("x", 1.00251, "±1%", "b", competing_value=1.0)
    row = compare(item, Measurement("x", 1.0, stat_err=0.00446))
    assert row.samples_needed_for_3sigma == pytest.approx((3 / 0.5617) ** 2, rel=1e-2)


# =============================================================================
# 단측 경계
# =============================================================================
def test_lower_bound_pass():
    item = PredictionItem("msd_r_squared", ">0.99", ">0.99", "단일지수 형태")
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
    assert "구별되지 않는다" in row.reason


def test_upper_bound_pass():
    item = PredictionItem("roundtrip_err", "<1e-12", "<1e-12", "S4 게이트")
    row = compare(item, Measurement("roundtrip_err", 1.6e-16, stat_err=1e-18))
    assert row.verdict == PASS


def test_non_numeric_prediction_with_two_sided_band_is_inconclusive():
    item = PredictionItem("shape", "single exponential", "±5%", "형태")
    assert compare(item, Measurement("shape", 1.0, stat_err=0.01)).verdict == INCONCLUSIVE


# =============================================================================
# ★ 첫 완주 9개 판정 재현
# =============================================================================
#  runs/2026-07-28_trap-2d-5um_2dfb9d/{02_prediction.md, metrics.json} 의 값.
#  runs/ 는 gitignore 대상이라 여기 박아둔다.
FIRST_RUN = [
    # (예측 항목, 측정, 문서에 기록된 verdict)
    (PredictionItem("var_x_2d_nm2", 414.19, "±1.5%", "⟨x²⟩ = kT/k (정확해)"),
     Measurement("var_x_2d_nm2", 416.5826, stat_err=1.8479, unit="nm^2"), PASS),
    (PredictionItem("var_x_3d_nm2", 414.19, "±1.5%", "⟨x²⟩ = kT/k, dim 무관"),
     Measurement("var_x_3d_nm2", 415.3363, stat_err=1.3685, unit="nm^2"), PASS),
    (PredictionItem("var_x_star", 1.0025063, "±1%", "EM 편향 1/(1−dt*/2)"),
     Measurement("var_x_star", 1.0057652, stat_err=0.0044614), PASS),
    (PredictionItem("var_r_2d_nm2", 828.39, "±1.5%", "⟨r²⟩ = d·kT/k"),
     Measurement("var_r_2d_nm2", 833.1652, stat_err=3.6958, unit="nm^2"), PASS),
    (PredictionItem("msd_plateau_2d_star", 4.0, "±2%", "plateau = 2d"),
     Measurement("msd_plateau_2d_star", 4.0214904, stat_err=0.0060703), PASS),
    (PredictionItem("tau_trap_ms", 8.0644, "±5%", "τ = γ/k"),
     Measurement("tau_trap_ms", 8.0566504, stat_err=0.0299748, unit="ms"), PASS),
    (PredictionItem("msd_r_squared", ">0.99", ">0.99", "단일지수 형태"),
     Measurement("msd_r_squared", 0.9999765, stat_err=1.0e-5), PASS),
    (PredictionItem("ks_p", "p>0.05", "p>0.05", "위치분포 Gaussian"),
     Measurement("ks_p", 0.29, stat_err=0.05), PASS),
    # P8 — 판별력 부족. 문서: 0.56σ
    (PredictionItem("em_bias_reproduced", 1.0025063, "±0.1%", "EM 편향 재현",
                    competing_value=1.0, discriminates="적분기 스킴"),
     Measurement("em_bias_reproduced", 1.0057652, stat_err=0.0044614), INCONCLUSIVE),
    # P9 — dt 래더. 문서: 0.24σ
    (PredictionItem("em_bias_halved", 1.0012516, "±0.1%", "dt 절반 → 편향 절반",
                    competing_value=1.0, discriminates="편향의 dt 선형성"),
     Measurement("em_bias_halved", 1.0016875, stat_err=0.0051896), INCONCLUSIVE),
]


@pytest.mark.benchmark
@pytest.mark.parametrize("item,meas,expected", FIRST_RUN,
                         ids=[i.quantity for i, _, _ in FIRST_RUN])
def test_reproduces_first_run_verdicts(item, meas, expected):
    """★ 사람이 손으로 내린 판정을 코드가 재현한다."""
    row = compare(item, meas)
    assert row.verdict == expected, row.reason


@pytest.mark.benchmark
def test_first_run_verdict_counts():
    """문서의 집계: 8 PASS · 2 INCONCLUSIVE · 0 FAIL (P1·P3 을 항목별로 분리한 셈)."""
    verdicts = [compare(i, m).verdict for i, m, _ in FIRST_RUN]
    assert verdicts.count(FAIL) == 0
    assert verdicts.count(INCONCLUSIVE) == 2
    assert verdicts.count(PASS) == len(FIRST_RUN) - 2


@pytest.mark.benchmark
def test_p8_design_power_matches_recorded_value():
    """문서에 기록된 판별력 `0.56σ` 와 `3σ 에 29배` 를 재현한다."""
    item, meas, _ = FIRST_RUN[8]
    row = compare(item, meas)
    assert row.design_power == pytest.approx(0.5618, rel=2e-3)
    assert row.samples_needed_for_3sigma == pytest.approx(28.5, rel=0.05)


# =============================================================================
# 봉인 — 깨지면 대조표를 만들지 않는다
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
    """★ 예측이 결과를 보고 수정됐을 수 있으면 대조표는 검증이 아니다."""
    sealed_rundir.write("prediction", "# S2\nvar_x = 416 nm^2 ±0.5%  (결과 보고 고침)\n")
    rep = V.validate_run(_prediction(), _measurements(), rundir=sealed_rundir)
    assert rep.rows == []
    assert rep.verdict_overall == "SEAL_BROKEN"
    assert any("봉인 위반" in p for p in rep.problems)


def test_missing_measurement_is_reported(sealed_rundir):
    """봉인된 예측을 조용히 빼놓을 수 없다."""
    rep = V.validate_run(_prediction(), {}, rundir=sealed_rundir)
    assert any("대응하는 측정이 없다" in p for p in rep.problems)


def test_unsealed_extra_measurement_is_reported(sealed_rundir):
    meas = {**_measurements(), "kurtosis": Measurement("kurtosis", 2.9968,
                                                       stat_err=0.012)}
    rep = V.validate_run(_prediction(), meas, rundir=sealed_rundir)
    assert any("봉인되지 않은 측정" in p for p in rep.problems)


def test_empty_prediction_is_reported(sealed_rundir):
    rep = V.validate_run(Prediction(items=[]), {}, rundir=sealed_rundir)
    assert any("0개" in p for p in rep.problems)


# =============================================================================
# 판정 블록 — confirmed_by 는 코드가 채우지 않는다
# =============================================================================
def test_yaml_block_never_confirms():
    """★ 이 테스트가 통과하지 않으면 사람이 보지 않은 합격 도장이 찍힌다."""
    rep = V.ValidationReport(rows=[compare(*FIRST_RUN[0][:2])])
    block = rep.yaml_block()
    assert "confirmed_by: null" in block
    assert "proposed_by: agent" in block
    assert "confirmed_by: agent" not in block


def test_no_public_api_sets_confirmed_by():
    """`confirmed_by` 를 채우는 코드 경로가 존재하지 않아야 한다."""
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
    """FAIL 에는 원인 4분류 중 하나가 붙어야 한다 (S7 게이트)."""
    pred = Prediction(items=[PredictionItem("var_x", 414.19, "±1.5%", "kT/k")])
    meas = {"var_x": Measurement("var_x", 600.0, stat_err=1.85)}
    rep = V.validate_run(pred, meas, rundir=sealed_rundir)
    assert rep.verdict_overall == FAIL
    assert any("원인 분류가 없다" in p for p in rep.problems)


def test_fail_with_cause_class_has_no_problem(sealed_rundir):
    pred = Prediction(items=[PredictionItem("var_x", 414.19, "±1.5%", "kT/k")])
    meas = {"var_x": Measurement("var_x", 600.0, stat_err=1.85)}
    rep = V.validate_run(pred, meas, rundir=sealed_rundir,
                         causes={"var_x": "numerical"})
    assert rep.problems == []


def test_table_and_reasons_render(sealed_rundir):
    rep = V.validate_run(_prediction(), _measurements(), rundir=sealed_rundir)
    assert "| 양 | 예측 (봉인) |" in rep.table()
    assert "PASS" in rep.table()
    assert rep.reasons() == "_없음 — 전 항목 PASS._"


def test_reasons_explain_inconclusive():
    rep = V.ValidationReport(rows=[compare(*FIRST_RUN[8][:2])])
    assert "INCONCLUSIVE" in rep.reasons()
    assert "예견된 한계" in rep.reasons()
