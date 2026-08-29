"""S7 — 봉인된 예측 vs 측정. **판정을 제안하고, 확정하지 않는다.**

## 이 모듈이 절대 하지 않는 것

`confirmed_by` 를 채우지 않는다. 그 필드는 사람만 쓴다 (CLAUDE.md §판정).
`confirmed_by: null` 인 판정은 벤치마크 원장 집계에 들어가지 않는다 — 그것이 세어지면
사람이 한 번도 보지 않은 합격 도장이 찍힌다.

## `INCONCLUSIVE` 는 실패가 아니다

첫 완주에서 예측 9개 중 2개가 `INCONCLUSIVE` 였고, **둘 다 예측 문서를 쓸 때 이미
판별력 부족이 예견되어 있었다.** 통계오차가 tolerance 보다 크면 그 측정은 그
tolerance 를 판정할 수 없다 — 그걸 `PASS` 라고 쓰면 검증이 아니라 요행이다.

두 가지를 따로 계산한다:

| | 질문 | 실패 시 |
|---|---|---|
| **판정 가능성** | `SE` 가 tolerance 대역보다 작은가 | `INCONCLUSIVE` |
| **설계 검정력** | 경쟁 가설과의 차이가 `SE` 의 몇 배인가 | `INCONCLUSIVE` + 필요 표본 배수 보고 |

★ 검정력이 `3σ` 를 못 만드는 곳에서 `3σ` 기각을 요구하지 않는다 — 달성 불가능한
assert 가 된다. 그 구간은 `INCONCLUSIVE` 를 **사실로 고정**한다 (CLAUDE.md 통계 4규칙).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from .io import RunDir, SealVerdict, verify_seal
from .spec import Prediction, PredictionItem

PASS, FAIL, INCONCLUSIVE = "PASS", "FAIL", "INCONCLUSIVE"

# FAIL 원인 4분류 (master_plan §S7-6). FAIL 에는 반드시 하나가 붙어야 한다.
CAUSE_CLASSES: frozenset[str] = frozenset({
    "numerical", "modeling", "interpretation", "analysis", "environment"})


# =============================================================================
# tolerance 파싱
# =============================================================================
@dataclass(frozen=True)
class Tolerance:
    """`"±1.5%"` 같은 문자열의 기계 표현.

    kind:
      relative     — 예측값의 ±X %
      absolute     — ±X (측정 단위)
      lower_bound  — 측정값이 이 값보다 커야 한다 (`R² > 0.99`, `p > 0.05`)
      upper_bound  — 측정값이 이 값보다 작아야 한다
    """

    kind: str
    magnitude: float
    text: str

    def half_width(self, predicted: float | None) -> float | None:
        """판정 대역의 반폭. 단측 경계에는 대역이 없으므로 `None`."""
        if self.kind == "absolute":
            return self.magnitude
        if self.kind == "relative":
            if predicted is None:
                return None
            return abs(predicted) * self.magnitude
        return None


_TOL_RE = re.compile(r"""
    ^\s*
    (?:(?P<name>[A-Za-z_][\w^]*)\s*)?           # p, R2, R^2 …
    (?P<op>±|\+/-|>=|<=|>|<)\s*
    (?P<num>[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)
    \s*(?P<pct>%)?\s*$
""", re.VERBOSE)


def parse_tolerance(text: str) -> Tolerance:
    """tolerance 문자열 → `Tolerance`. 못 읽으면 **예외** (조용히 넘기지 않는다).

    tolerance 를 읽지 못하고 넘어가면 그 항목은 판정 없이 통과한다 — 검증을
    무력화하는 가장 쉬운 방법이다 (master_plan §S2 실패모드).
    """
    m = _TOL_RE.match(str(text))
    if not m:
        raise ValueError(
            f"tolerance {text!r} 를 읽을 수 없다. 허용 형식: "
            f"'±1.5%', '±0.03', '>0.99', 'p>0.05', '<1e-3'")
    num = float(m.group("num"))
    op = m.group("op")
    if op in ("±", "+/-"):
        if num <= 0:
            raise ValueError(f"tolerance {text!r}: 반폭이 0 이하다 — "
                             f"어떤 결과든 FAIL 이 된다")
        return Tolerance("relative" if m.group("pct") else "absolute",
                         num / 100.0 if m.group("pct") else num, str(text))
    kind = "lower_bound" if op in (">", ">=") else "upper_bound"
    return Tolerance(kind, num / 100.0 if m.group("pct") else num, str(text))


# =============================================================================
# 측정값
# =============================================================================
@dataclass
class Measurement:
    """측정 1건. **`stat_err` 없는 값은 결론에 쓸 수 없다** (CLAUDE.md).

    `stat_err=None` 은 "오차를 모른다"이고, 그 상태에서는 판정이 나가지 않는다.
    시드 1개짜리 프로덕션 런이 여기서 막힌다.
    """

    quantity: str
    value: float
    stat_err: float | None = None
    method: str = ""
    n_samples: int | None = None
    spread: float | None = None          # 시드간 산포 — 계통오차 지표
    unit: str = ""

    def problems(self) -> list[str]:
        out = []
        if self.stat_err is None:
            out.append(f"{self.quantity}: 통계오차가 없다 — 시드 1개 런이거나 "
                       f"오차를 계산하지 않았다. 결론에 쓸 수 없다")
        elif self.stat_err <= 0:
            out.append(f"{self.quantity}: stat_err = {self.stat_err} — 요동하지 않는 "
                       f"'측정값'은 산술 항등식일 가능성이 크다 "
                       f"(guards.assert_statistic_fluctuates 참조)")
        if not math.isfinite(self.value):
            out.append(f"{self.quantity}: 측정값이 비유한값이다")
        return out


# =============================================================================
# 판정 1행
# =============================================================================
@dataclass
class ValidationRow:
    quantity: str
    predicted: float | str
    measured: float
    tolerance: str
    verdict: str
    stat_err: float | None = None
    deviation: float | None = None        # 절대 편차
    deviation_rel: float | None = None    # 상대 편차
    sigma: float | None = None            # 편차 / SE
    design_power: float | None = None     # |예측 − 경쟁| / SE
    samples_needed_for_3sigma: float | None = None
    cause_class: str = ""
    reason: str = ""
    note: str = ""
    flags: list[str] = field(default_factory=list)

    def problems(self) -> list[str]:
        out = []
        if self.verdict == FAIL and self.cause_class not in CAUSE_CLASSES:
            out.append(f"{self.quantity}: FAIL 인데 원인 분류가 없다 — "
                       f"{sorted(CAUSE_CLASSES)} 중 하나가 필요하다 (S7 게이트)")
        if "significant_deviation_within_tolerance" in self.flags:
            out.append(
                f"{self.quantity}: **PASS 이지만 편차가 {self.sigma:.2f}σ** 다 — "
                f"tolerance 대역 안이지만 통계적으로 유의한 어긋남이다. "
                f"tolerance 가 너무 넓어 검증이 무력화됐는지, 아니면 예측에 넣지 않은 "
                f"알려진 편향이 있는지 확인할 것 (master_plan §S2 실패모드)")
        return out


def compare(item: PredictionItem, meas: Measurement, *,
            power_threshold: float = 1.0,
            significance_sigma: float = 3.0) -> ValidationRow:
    """예측 1개 vs 측정 1개 → 판정 **제안** 1행.

    Args:
        power_threshold: 설계 검정력이 이 값 미만이면 `INCONCLUSIVE`.
            기본 `1.0` — 경쟁 가설과의 차이가 `1σ` 도 안 되면 그 측정은 두 가설을
            구별하지 못한다. `3.0` 으로 올리면 더 보수적이다.
        significance_sigma: PASS 인데도 편차가 이만큼 유의하면 표시한다.
            ★ 넓은 tolerance 는 어떤 결과든 PASS 로 만든다 (§S2 실패모드). 대역 안에
            들어왔다는 것과 예측이 맞았다는 것은 다른 말이다 — 3σ 어긋난 PASS 는
            "예측에 넣지 않은 편향이 있다"는 신호다.
    """
    tol = parse_tolerance(item.tolerance)
    se = meas.stat_err
    row = ValidationRow(quantity=item.quantity, predicted=item.value,
                        measured=meas.value, tolerance=item.tolerance,
                        verdict=INCONCLUSIVE, stat_err=se, note=item.note)

    # --- 오차를 모르면 판정하지 않는다 ---
    if se is None:
        row.reason = "통계오차 없음 — 오차 막대 없는 수를 결론에 쓰지 않는다"
        return row

    # --- 단측 경계 (R² > 0.99, p > 0.05) ---
    if tol.kind in ("lower_bound", "upper_bound"):
        bound = tol.magnitude
        row.deviation = meas.value - bound
        row.sigma = abs(row.deviation) / se if se > 0 else None
        ok = meas.value > bound if tol.kind == "lower_bound" else meas.value < bound
        if se > 0 and abs(meas.value - bound) < se:
            row.verdict = INCONCLUSIVE
            row.reason = (f"경계 {bound:g} 에서 {abs(row.deviation):.3g} 떨어져 있고 "
                          f"SE = {se:.3g} — 경계와 구별되지 않는다")
        else:
            row.verdict = PASS if ok else FAIL
            row.reason = (f"측정 {meas.value:.6g} "
                          f"{'>' if tol.kind == 'lower_bound' else '<'} {bound:g}"
                          + ("" if ok else " 위반"))
        return row

    # --- 양측 대역 ---
    if not isinstance(item.value, (int, float)):
        row.reason = f"예측값이 수치가 아니다 ({item.value!r}) — 양측 대역을 쓸 수 없다"
        return row

    pred = float(item.value)
    half = tol.half_width(pred)
    row.deviation = meas.value - pred
    row.deviation_rel = row.deviation / pred if pred else None
    row.sigma = abs(row.deviation) / se if se > 0 else None

    # ① 판정 가능성 — SE 가 대역보다 크면 이 측정은 이 tolerance 를 판정할 수 없다
    if half is not None and se > half:
        row.verdict = INCONCLUSIVE
        row.reason = (f"SE = {se:.4g} > tolerance 반폭 {half:.4g} — "
                      f"이 측정은 이 대역을 판정할 수 없다. "
                      f"필요 표본 배수 ≈ {(se / half) ** 2:.3g}×")
    else:
        inside = half is not None and abs(row.deviation) <= half
        row.verdict = PASS if inside else FAIL
        row.reason = (f"편차 {row.deviation:+.4g}"
                      + (f" ({row.deviation_rel:+.3%})" if row.deviation_rel else "")
                      + f", 대역 ±{half:.4g}")

    # ② 설계 검정력 — 경쟁 가설과 구별되는가
    if item.competing_value is not None and se > 0:
        gap = abs(pred - float(item.competing_value))
        row.design_power = gap / se
        row.samples_needed_for_3sigma = (3.0 / row.design_power) ** 2 \
            if row.design_power > 0 else float("inf")
        if row.design_power < power_threshold:
            row.verdict = INCONCLUSIVE
            row.reason = (
                f"설계 검정력 {row.design_power:.2f}σ < {power_threshold:g}σ — "
                f"경쟁 가설({item.competing_value:g})과 구별되지 않는다. "
                f"3σ 판별에 표본 {row.samples_needed_for_3sigma:.3g}배 필요. "
                f"★ 예견된 한계이지 실패가 아니다")

    # ③ 넓은 tolerance 로 가려진 유의한 어긋남
    if (row.verdict == PASS and row.sigma is not None
            and row.sigma > significance_sigma):
        row.flags.append("significant_deviation_within_tolerance")
    return row


# =============================================================================
# 리포트
# =============================================================================
@dataclass
class ValidationReport:
    rows: list[ValidationRow]
    seal: SealVerdict | None = None
    problems: list[str] = field(default_factory=list)
    sanity: list[ValidationRow] = field(default_factory=list)

    # --- 집계 ---
    def count(self, verdict: str) -> int:
        return sum(1 for r in self.all_rows() if r.verdict == verdict)

    def all_rows(self) -> list[ValidationRow]:
        return [*self.rows, *self.sanity]

    @property
    def verdict_overall(self) -> str:
        if self.seal is not None and not self.seal.ok:
            return "SEAL_BROKEN"
        if self.count(FAIL):
            return "FAIL"
        if self.count(INCONCLUSIVE):
            return "PASS_WITH_INCONCLUSIVE"
        return PASS

    def yaml_block(self) -> str:
        """판정 블록. ★ `confirmed_by` 는 **항상** `null` 이다 — 코드가 채우지 않는다."""
        return (
            "```yaml\n"
            f"verdict_overall: {self.verdict_overall}\n"
            "proposed_by: agent\n"
            "confirmed_by: null            # ← 사람 확정 대기\n"
            f"pass: {self.count(PASS)}\n"
            f"inconclusive: {self.count(INCONCLUSIVE)}\n"
            f"fail: {self.count(FAIL)}\n"
            "```")

    def table(self) -> str:
        rows = ["| 양 | 예측 (봉인) | 측정 | tolerance | 편차 | 검정력 | verdict |",
                "|---|---|---|---|---|---|---|"]
        mark = {PASS: "**PASS**", FAIL: "❌ **FAIL**",
                INCONCLUSIVE: "⚠ **INCONCLUSIVE**"}
        for r in self.all_rows():
            suffix = " ⚑" if r.flags else ""
            pred = f"`{r.predicted:.6g}`" if isinstance(r.predicted, (int, float)) \
                else f"`{r.predicted}`"
            meas = f"`{r.measured:.6g}`"
            if r.stat_err is not None:
                meas += f" ± `{r.stat_err:.3g}`"
            dev = ""
            if r.deviation_rel is not None:
                dev = f"{r.deviation_rel:+.2%}"
            elif r.deviation is not None:
                dev = f"{r.deviation:+.3g}"
            if r.sigma is not None:
                dev += f" ({r.sigma:.2f}σ)"
            power = f"`{r.design_power:.2f}σ`" if r.design_power is not None else "—"
            rows.append(f"| `{r.quantity}` | {pred} | {meas} | `{r.tolerance}` "
                        f"| {dev} | {power} | {mark.get(r.verdict, r.verdict)}{suffix} |")
        if any(r.flags for r in self.all_rows()):
            rows.append("")
            rows.append("⚑ = PASS 이지만 편차가 통계적으로 유의하다 — 아래 참조.")
        return "\n".join(rows)

    def reasons(self) -> str:
        """`INCONCLUSIVE`·`FAIL` 의 이유를 풀어 쓴다. PASS 는 생략."""
        out = []
        for r in self.all_rows():
            if r.verdict == PASS:
                continue
            head = f"- **`{r.quantity}`** — {r.verdict}"
            if r.cause_class:
                head += f" (`{r.cause_class}`)"
            out.append(f"{head}: {r.reason}")
        return "\n".join(out) if out else "_없음 — 전 항목 PASS._"

    def notes(self) -> str:
        """예측 문서에 달린 부기. **PASS 라도 남긴다** — 여기 "이건 독립 검사가
        아니다" 같은 한계가 들어 있고, 그게 사라지면 결론이 과대해진다."""
        out = [f"- **`{r.quantity}`**: {r.note}" for r in self.all_rows() if r.note]
        return "\n".join(out) if out else ""


def validate_run(prediction: Prediction, measurements: dict[str, Measurement], *,
                 rundir: RunDir | None = None, power_threshold: float = 1.0,
                 causes: dict[str, str] | None = None) -> ValidationReport:
    """봉인 검증 → 예측 대조 → 판정 제안.

    ⚠ **봉인이 깨졌으면 판정을 만들지 않는다.** 예측이 결과를 보고 나서 수정됐을
      가능성이 있는 상태에서 대조표를 그리면, 그 표는 검증처럼 보이지만 검증이 아니다.

    Args:
        measurements: `quantity` → 측정값. 예측에 대응하는 측정이 없으면 문제로 보고.
        causes: FAIL 항목의 원인 분류 (에이전트가 공급). 없으면 게이트가 잡는다.
    """
    problems: list[str] = list(prediction.problems())
    seal = verify_seal(rundir) if rundir is not None else None

    if seal is not None and not seal.ok:
        return ValidationReport(
            rows=[], seal=seal,
            problems=[f"★ 봉인 위반 — {seal.summary()}. "
                      f"예측이 실행 후 수정됐을 수 있으므로 **대조표를 만들지 않는다** "
                      f"(master_plan §S7-1)."] + problems)

    rows: list[ValidationRow] = []
    for item in prediction.items:
        meas = measurements.get(item.quantity)
        if meas is None:
            problems.append(f"예측 {item.quantity!r} 에 대응하는 측정이 없다 — "
                            f"봉인된 예측을 조용히 빼놓을 수 없다")
            continue
        problems += meas.problems()
        row = compare(item, meas, power_threshold=power_threshold)
        if causes and item.quantity in causes:
            row.cause_class = causes[item.quantity]
        rows.append(row)
        problems += row.problems()

    extra = sorted(set(measurements) - {i.quantity for i in prediction.items})
    if extra:
        problems.append(f"봉인되지 않은 측정 {extra} — 사후에 추가된 항목은 "
                        f"'예측 대조'가 아니라 '관찰'로 따로 보고할 것")

    return ValidationReport(rows=rows, seal=seal, problems=problems)
