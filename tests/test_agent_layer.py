"""L1 에이전트 층 (`.claude/`) 의 구조 검사.

## 왜 마크다운을 테스트하는가

이 층은 코드가 아니라 문서다. 그런데 **frontmatter 는 기계가 읽는 계약**이고,
링크가 깨지면 에이전트가 프로토콜을 못 찾고, 티어링 경계가 빠지면 값싼 모델이
물리 판단을 하게 된다. 전부 조용히 실패하는 종류다.

## 검사하는 것

1. frontmatter 유효 (`name`·`description`·`model`), 이름과 파일명 일치
2. 상대 링크가 실제로 존재
3. **모델 티어링** — §12.2 배분표와 일치, 저가 모델에 권한 경계 명시
4. **스킬이 물리를 다시 적지 않는다** — 코어 함수를 인용하는가
5. `settings.json` 이 봉인 문서 편집을 거부하는가
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
import yaml

CLAUDE = Path(".claude")
pytestmark = pytest.mark.skipif(not CLAUDE.exists(), reason=".claude/ 없음")

SKILLS = sorted(CLAUDE.glob("skills/*/SKILL.md"))
AGENTS = sorted(CLAUDE.glob("agents/*.md"))

# 2026-08-28 병합으로 스킬이 두 부류가 됐다. 섞어서 한 규칙으로 재면 안 된다 —
# 오케스트레이션 스킬끼리는 **서로 배타적**이라 잘못 트리거되면 엉뚱한 절차가
# 돌아가므로 상호참조가 필수다. 도메인 스킬은 파이프라인이 단계마다 *불러 읽는*
# 참조물이라 배타적이지 않고, 대신 파이프라인이 그것들을 가리켜야 한다.
ORCH_SKILLS = {"bd-pipeline", "bd-diagnose", "bd-knowledge"}
DOMAIN_SKILLS = {"bd-intake", "bd-physics", "bd-hoomd"}

# master_plan §12.1 배분표. 여기서 벗어나면 설계와 코드가 갈린 것이다.
EXPECTED_MODELS = {
    "bd-intake-extract": "haiku",
    "bd-intake-interpret": "opus",
    "bd-predict": "opus",
    "bd-spec": "sonnet",
    "bd-validate": "opus",
    "bd-conclude": "opus",
    "bd-lit-distill": "sonnet",
    "bd-lit-scan": "haiku",
    "bd-diagnose": "opus",
}


def frontmatter(path: Path) -> dict:
    m = re.match(r"^---\n(.*?)\n---\n", path.read_text(encoding="utf-8"), re.S)
    assert m, f"{path}: frontmatter 가 없다 — 기계가 읽는 계약이다"
    return yaml.safe_load(m.group(1))


# =============================================================================
# 존재
# =============================================================================
def test_all_six_skills_exist():
    assert {p.parent.name for p in SKILLS} == ORCH_SKILLS | DOMAIN_SKILLS


def test_all_nine_agents_exist():
    assert {p.stem for p in AGENTS} == set(EXPECTED_MODELS)


def test_pipeline_has_all_reference_docs():
    refs = {p.name for p in (CLAUDE / "skills/bd-pipeline/references").glob("*.md")}
    assert refs == {"s1_intake_drawing.md", "s2_prediction.md", "s3_s5_execute.md",
                    "s6_s7_validate.md", "s8_knowledge.md"}


# =============================================================================
# frontmatter
# =============================================================================
@pytest.mark.parametrize("path", SKILLS, ids=lambda p: p.parent.name)
def test_skill_frontmatter(path):
    fm = frontmatter(path)
    assert fm["name"] == path.parent.name, "name 이 디렉터리명과 달라 호출되지 않는다"
    assert len(fm["description"]) > 40, "description 이 짧으면 트리거되지 않는다"


@pytest.mark.parametrize("path", AGENTS, ids=lambda p: p.stem)
def test_agent_frontmatter(path):
    fm = frontmatter(path)
    assert fm["name"] == path.stem
    assert fm["model"] in ("opus", "sonnet", "haiku", "inherit")
    assert len(fm["description"]) > 40


@pytest.mark.parametrize("path", AGENTS, ids=lambda p: p.stem)
def test_agent_model_matches_the_tiering_table(path):
    """★ master_plan §12.1 배분표와 일치하는가.

    어긋나면 설계 문서와 실제 위임이 갈린 것이다 — 어느 쪽이 맞는지 결정해야 한다.
    """
    fm = frontmatter(path)
    assert fm["model"] == EXPECTED_MODELS[path.stem]


def test_orch_skill_descriptions_say_when_not_to_use():
    """배타적인 스킬끼리는 서로를 가리켜야 잘못된 스킬이 트리거되지 않는다."""
    for path in SKILLS:
        if path.parent.name not in ORCH_SKILLS:
            continue
        desc = frontmatter(path)["description"]
        others = ORCH_SKILLS - {path.parent.name}
        assert any(o in desc for o in others), f"{path.parent.name}: 다른 스킬 언급 없음"


def test_pipeline_points_at_every_domain_skill():
    """도메인 스킬은 스스로 트리거되지 않는다 — 파이프라인이 불러야 읽힌다.

    병합 직후 이것이 실제로 비어 있었다: S5 가 HOOMD 코드를 쓰는데 함정 20개를
    담은 bd-hoomd 를 아무도 가리키지 않았다. 이 검사가 그 공백을 지킨다.
    """
    body = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    missing = sorted(s for s in DOMAIN_SKILLS if s not in body)
    assert not missing, f"bd-pipeline 이 가리키지 않는 도메인 스킬: {missing}"


# =============================================================================
# 링크
# =============================================================================
def _links(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    return [m.group(1).split("#")[0]
            for m in re.finditer(r"\[[^\]]+\]\((\.\.?/[^)#]*)\)", text)]


@pytest.mark.parametrize("path", sorted(CLAUDE.rglob("*.md")),
                         ids=lambda p: str(p.relative_to(CLAUDE)))
def test_relative_links_resolve(path):
    """깨진 링크는 에이전트가 프로토콜을 못 찾는다는 뜻이다."""
    for link in _links(path):
        assert (path.parent / link).resolve().exists(), f"{path}: {link}"


def test_pipeline_skill_links_every_reference():
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    for ref in (CLAUDE / "skills/bd-pipeline/references").glob("*.md"):
        assert ref.name in text, f"SKILL.md 가 {ref.name} 를 인용하지 않는다"


# =============================================================================
# 티어링 안전장치 — 저가 모델이 물리 판단을 하지 못하게
# =============================================================================
@pytest.mark.parametrize(
    "path", [p for p in AGENTS if EXPECTED_MODELS[p.stem] != "opus"],
    ids=lambda p: p.stem)
def test_cheap_agents_state_their_authority_boundary(path):
    """★ `provenance` 가 `inference`/`assumed` 인 필드는 Opus 만 (§12.2).

    저가 모델의 지시문에 이 경계가 없으면 그 모델이 물리 가정을 확정하게 된다.
    """
    text = path.read_text(encoding="utf-8")
    assert "§12.2" in text, "권한 경계 근거(§12.2) 인용 없음"
    assert "inference" in text and "assumed" in text
    assert "Opus" in text, "누구에게 넘겨야 하는지 명시 없음"


def test_extract_agent_cannot_write_files():
    """추출 에이전트에 Write 를 주면 spec 을 직접 고칠 수 있다."""
    tools = frontmatter(CLAUDE / "agents/bd-intake-extract.md")["tools"]
    assert "Write" not in tools


# =============================================================================
# 스킬이 물리를 다시 적지 않는다
# =============================================================================
def test_skills_delegate_to_the_core_rather_than_restating_physics():
    """스킬은 `simbot` 함수와 `cli.py` 를 인용한다 — 숫자를 만들지 않는다."""
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "cli.py run" in text
    assert "simbot" in text
    assert "숫자를 머리로" in text            # 금지 규약이 명시돼야 한다


def test_pipeline_skill_forbids_confirming_verdicts():
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "confirmed_by" in text
    assert "사람만" in text or "사람이 확정" in text


def test_pipeline_skill_names_the_interpreter_absolutely():
    """`conda activate` 는 non-interactive shell 에서 불안정하다."""
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python" in text


def test_pipeline_skill_states_the_question_budget():
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "질문 예산" in text and "3개" in text


def test_s1_reference_is_the_substantive_one():
    """S1 판독은 코드로 표현할 수 없는 유일한 단계 — 문서가 가장 두꺼워야 한다."""
    refs = {p.name: len(p.read_text(encoding="utf-8"))
            for p in (CLAUDE / "skills/bd-pipeline/references").glob("*.md")}
    assert refs["s1_intake_drawing.md"] > 3000


def test_s1_reference_warns_against_trusting_absolute_size():
    text = (CLAUDE / "skills/bd-pipeline/references/s1_intake_drawing.md").read_text(
        encoding="utf-8")
    assert "절대 크기를 신뢰하지 않는다" in text
    assert "후보" in text                     # 모호성 후보 나열 규약


def test_diagnose_skill_suspects_analysis_first():
    """★ 첫 완주 4건 중 물리 문제는 0건이었다."""
    text = (CLAUDE / "skills/bd-diagnose/SKILL.md").read_text(encoding="utf-8")
    assert "analysis" in text
    assert "먼저 의심" in text
    assert text.index("analysis") < text.index("modeling")


def test_knowledge_skill_states_the_citation_discipline():
    text = (CLAUDE / "skills/bd-knowledge/SKILL.md").read_text(encoding="utf-8")
    assert "미재현" in text
    assert "기억으로 인용" in text


# =============================================================================
# settings.json
# =============================================================================
def test_settings_is_valid_json():
    json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))


def test_settings_allows_the_project_interpreter():
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    assert any("simulation_bot/bin/python" in a for a in s["permissions"]["allow"])


def test_settings_denies_conda_activate():
    """CLAUDE.md: non-interactive shell 에서 불안정하다."""
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    assert any("conda activate" in d for d in s["permissions"]["deny"])


def test_settings_denies_editing_sealed_documents():
    """★ 사후합리화의 가장 쉬운 경로를 구조적으로 닫는다.

    봉인 검증이 사후에 잡지만, 애초에 못 하게 하는 것이 낫다.
    """
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    deny = " ".join(s["permissions"]["deny"])
    for sealed in ("02_prediction.md", "01_intake.md", "SEALED.sha256"):
        assert sealed in deny, f"{sealed} 편집이 허용돼 있다"


def test_settings_does_not_deny_reading_inputs():
    """S1 이 손그림을 읽어야 한다 — 여기를 막으면 파이프라인이 죽는다."""
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    deny = " ".join(s["permissions"]["deny"])
    assert "Read(./inputs" not in deny
    assert "Read(./knowledge" not in deny


# =============================================================================
# 문서화된 결정
# =============================================================================
def test_readme_records_the_granularity_decision():
    """참조문서를 8개→5개로 줄인 이유가 기록돼 있어야 한다 (master_plan Q6)."""
    text = (CLAUDE / "README.md").read_text(encoding="utf-8")
    assert "Q6" in text
    assert "5개" in text


def test_readme_lists_the_tiering_table():
    text = (CLAUDE / "README.md").read_text(encoding="utf-8")
    for name in EXPECTED_MODELS:
        assert name in text, name
