"""S8 — `REPORT.md` 생성. LLM 0줄.

**리포트는 요약이 아니라 감사 기록이다.** 사람이 이 파일 하나만 읽고
"이 결론을 믿어도 되는가"에 답할 수 있어야 한다. 그래서 좋은 소식만 담지 않는다:

- 봉인 상태 (깨졌으면 **맨 위에**)
- `INCONCLUSIVE` 와 그 이유
- 재현 가능성 (`git_dirty` 포함)
- 아직 판정되지 않은 게이트
- `confirmed_by: null`

`simbot` 이 만들 수 없는 것 — 질문에 대한 답, 원인 가설, 다음 실험 제안 — 은
에이전트가 쓴 `08_conclusion.md` 를 **인용**한다. 지어내지 않는다.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .io import RUN_LAYOUT, RunDir, verify_seal
from .nondim import ReducedSpec, nondim_table, reduce_spec, roundtrip_errors
from .spec import SpecReport, SystemSpec, validate as validate_spec
from .validate import ValidationReport

ROUNDTRIP_GATE = 1e-12


@dataclass
class ReportInputs:
    """리포트가 필요로 하는 것. 없는 것은 `None` 으로 두면 그 절이 '없음'으로 나온다."""

    spec: SystemSpec | None = None
    spec_report: SpecReport | None = None
    reduced: ReducedSpec | None = None
    validation: ValidationReport | None = None
    manifest: dict | None = None
    figures: dict[str, str] | None = None      # 파일명 → 캡션
    wall_s: float | None = None
    n_runs: int | None = None


def _fmt(x: float | None, spec: str = ".6g") -> str:
    return "—" if x is None else f"`{x:{spec}}`"


# =============================================================================
# 절 단위 렌더러 — 각각 독립적으로 테스트된다
# =============================================================================
def seal_section(rundir: RunDir) -> str:
    v = verify_seal(rundir)
    if v.ok:
        return (f"**봉인 검증** ✅ {v.summary()}\n\n"
                f"검증 명령 (이 코드 없이도 확인된다):\n\n"
                f"```bash\nshasum -a 256 -c {RUN_LAYOUT['seal']}\n```")
    return ("> ## ⛔ seal violation\n>\n"
            f"> {v.summary()}\n>\n"
            "> **예측이 실행 후 수정됐을 수 있다.** 아래 대조표는 검증으로 읽으면 안 된다.\n"
            "> 근거: master_plan §S7-1.")


def reproducibility_section(manifest: dict | None) -> str:
    if not manifest:
        return "_manifest 없음 — 재현 가능성을 주장할 수 없다._"
    rows = ["| 항목 | 값 |", "|---|---|"]
    for key, label in (("run_id", "run_id"), ("spec_hash", "spec 해시"),
                       ("code_hash", "코드 해시"), ("git_rev", "git 커밋"),
                       ("env_hash", "env 해시"), ("seed", "시드")):
        if key in manifest:
            rows.append(f"| {label} | `{manifest[key]}` |")
    env = manifest.get("env", {})
    for k in ("hoomd", "numpy", "scipy", "python"):
        if k in env:
            rows.append(f"| {k} | `{env[k]}` |")

    dirty = manifest.get("git_dirty")
    if dirty is True:
        rows.append("| **git 상태** | ⚠️ **커밋되지 않은 변경 있음** — "
                    "`git_rev` 만으로 이 런이 재현되지 않는다 |")
    elif dirty is False:
        rows.append("| git 상태 | ✅ clean — `git_rev` 로 코드가 특정된다 |")
    else:
        rows.append("| git 상태 | ❔ 판정 불가 |")
    return "\n".join(rows)


def gates_section(spec_report: SpecReport | None) -> str:
    if spec_report is None:
        return "_S3 검사 결과 없음._"
    out = [spec_report.table()]
    deferred = spec_report.deferred()
    if deferred:
        out.append("")
        out.append(f"⏳ **S7 이 판정해야 하는 게이트 {len(deferred)}개** — "
                   f"S3 에서는 계산할 수 없는 양이다: "
                   + ", ".join(f"`{c.name}`" for c in deferred))
    if spec_report.problems:
        out.append("")
        out.append("### ⚠️ 규약 위반")
        out += [f"- {p}" for p in spec_report.problems]
    return "\n".join(out)


def nondim_section(spec: SystemSpec | None, reduced: ReducedSpec | None) -> str:
    if spec is None or reduced is None:
        return "_S4 무차원화 결과 없음._"
    errs = roundtrip_errors(spec, reduced)
    worst = max(errs.values()) if errs else 0.0
    gate = "✅ 통과" if worst < ROUNDTRIP_GATE else "❌ **위반**"
    out = [f"기준 척도: **{reduced.scales.origin}**", "",
           nondim_table(spec, reduced), "",
           f"**왕복 오차** 최대 `{worst:.2e}` (게이트 `< {ROUNDTRIP_GATE:g}`) — {gate}", "",
           f"`Δt*` = `{reduced.dt_star:.6g}` · 지배 제약 **{reduced.dt_dominant}**"]
    if reduced.logged:
        out.append("")
        out.append("기록용 (게이트 아님 — 다른 논문과 비교할 때만 쓴다):")
        out += [f"- `{k}` = `{v:.4g}`" for k, v in reduced.logged.items()]
    if reduced.groups:
        out.append("")
        out.append("| 무차원수 | 값 |")
        out.append("|---|---|")
        out += [f"| `{k}` | `{v:.6g}` |" for k, v in reduced.groups.items()]
    return "\n".join(out)


def validation_section(validation: ValidationReport | None) -> str:
    if validation is None:
        return "_S7 판정 없음._"
    out = [validation.table(), ""]
    n_pass = validation.count("PASS")
    n_inc = validation.count("INCONCLUSIVE")
    n_fail = validation.count("FAIL")
    out.append(f"**{len(validation.all_rows())}개 중 {n_pass} PASS · "
               f"{n_inc} INCONCLUSIVE · {n_fail} FAIL.**")
    if n_inc or n_fail:
        out += ["", "### PASS 아닌 항목", validation.reasons()]
    if n_inc:
        out += ["", "> `INCONCLUSIVE` 는 실패가 아니다 — 통계오차가 tolerance 보다 "
                    "커서 **판정이 불가능**하다는 사실이다. 필요한 표본 배수가 위에 있다."]
    notes = validation.notes()
    if notes:
        out += ["", "### 부기 (PASS 항목 포함)", notes,
                "", "> 예측 문서에 적힌 한계다. PASS 라고 사라지지 않는다 — "
                    "여기 '이건 독립 검사가 아니다' 류가 들어 있고, 그게 빠지면 "
                    "결론이 과대해진다."]
    if validation.problems:
        out += ["", "### ⚠️ 검증 절차의 문제", *[f"- {p}" for p in validation.problems]]
    return "\n".join(out)


def figures_section(rundir: RunDir, figures: dict[str, str] | None) -> str:
    """캡션 없는 그림은 산출물로 인정하지 않는다 (master_plan §S6 게이트)."""
    files = sorted(p.name for p in rundir.figs.glob("*.png")) if rundir.figs.exists() \
        else []
    if not files:
        return "_그림 없음._"
    caps = figures or {}
    out = []
    for f in files:
        cap = caps.get(f)
        if cap:
            # alt 텍스트는 한 줄로 자른다 — 여러 줄 캡션을 `![...]` 안에 넣으면
            # 마크다운 이미지 문법이 깨진다. 전체 캡션은 그림 아래로.
            alt = " ".join(cap.split())
            if len(alt) > 90:
                alt = alt[:89].rstrip() + "…"
            out.append(f"### {f}\n\n![{alt}](figs/{f})\n\n{cap}")
        else:
            out.append(f"### {f}\n\n![{f}](figs/{f})\n\n"
                       f"⚠️ **캡션 없음** — 무엇을 보이려는 그림인지 적어야 한다 "
                       f"(§S6 게이트).")
    return "\n\n".join(out)


def cost_section(wall_s: float | None, n_runs: int | None) -> str:
    if wall_s is None:
        return "_계산 비용 기록 없음._"
    parts = [f"**총 계산시간 `{wall_s:.1f} s`**"]
    if n_runs:
        parts.append(f"{n_runs}런")
        parts.append(f"런당 평균 `{wall_s / n_runs:.2f} s`")
    return " · ".join(parts)


def _excerpt(rundir: RunDir, stage: str, title: str) -> str:
    """에이전트가 쓴 문서를 **인용**한다. 없으면 없다고 적는다 — 지어내지 않는다."""
    if not rundir.exists(stage):
        return f"_`{RUN_LAYOUT[stage]}` 가 없다 — {title}은 에이전트가 써야 한다._"
    return f"→ [`{RUN_LAYOUT[stage]}`]({RUN_LAYOUT[stage]})"


# =============================================================================
# 조립
# =============================================================================
def render(rundir: RunDir, inputs: ReportInputs) -> str:
    spec = inputs.spec
    sr = inputs.spec_report
    if sr is None and spec is not None:
        sr = validate_spec(spec)
    reduced = inputs.reduced
    if reduced is None and spec is not None:
        try:
            reduced = reduce_spec(spec)
        except Exception:                      # 척도 미등록 카드 등 — 절만 비운다
            reduced = None

    v = inputs.validation
    headline = v.verdict_overall if v is not None else "판정 없음"

    parts = [
        f"# REPORT — `{rundir.run_id}`",
        "",
        f"**판정 (제안)** `{headline}` · **확정 대기** (`confirmed_by: null`)",
    ]
    if spec is not None:
        parts += ["", f"**질문** {spec.question.strip()}",
                  f"**카드** [`{spec.card}`]"
                  f"(../../knowledge/wiki/systems/{spec.card}.md)"]
    parts += ["", cost_section(inputs.wall_s, inputs.n_runs), "",
              seal_section(rundir), "", "---", ""]

    # 봉인이 깨졌으면 대조표를 리포트에 싣지 않는다
    seal_ok = verify_seal(rundir).ok

    parts += ["## 1. 판정 요약", ""]
    if not seal_ok:
        parts += ["대조표는 **생성하지 않았다** — 봉인이 깨졌기 때문이다.", ""]
    else:
        parts += [validation_section(v), ""]
    if v is not None:
        parts += [v.yaml_block(), "",
                  "> 판정은 **제안**이다. 사람이 `confirmed_by` 를 채우기 전까지 "
                  "벤치마크 원장 집계에 들어가지 않는다 (CLAUDE.md §판정).", ""]

    parts += ["---", "", "## 2. 시스템 명세와 게이트 (S3)", "",
              gates_section(sr), "",
              "---", "", "## 3. 무차원화 (S4)", "",
              nondim_section(spec, reduced), "",
              "---", "", "## 4. 그림 (S6)", "",
              figures_section(rundir, inputs.figures), "",
              "---", "", "## 5. 재현 가능성", "",
              reproducibility_section(inputs.manifest), "",
              "---", "", "## 6. 에이전트가 쓴 문서", "",
              f"- 판독 (S1) {_excerpt(rundir, 'intake', '판독')}",
              f"- 예측 (S2, 봉인) {_excerpt(rundir, 'prediction', '예측')}",
              f"- 검증 서술 (S7) {_excerpt(rundir, 'validation', '검증 서술')}",
              f"- 결론 (S8) {_excerpt(rundir, 'conclusion', '결론')}",
              "",
              "> `simbot` 은 수치를 만들고, 질문에 대한 답·원인 가설·다음 실험 제안은 "
              "에이전트가 쓴다. 이 리포트는 그 문서를 **인용**하며 대신 쓰지 않는다.",
              ""]
    return "\n".join(parts)


def write_report(rundir: RunDir, inputs: ReportInputs) -> Path:
    return rundir.write("report", render(rundir, inputs))
