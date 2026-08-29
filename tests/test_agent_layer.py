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
    """Only Opus may write a field whose `provenance` is `inference`/`assumed`.

    Without this boundary in a cheap model's instructions, that model ends up
    settling a physical assumption -- and every number downstream inherits it
    with nothing able to tell it later from a measured one.

    The citation used to be the predecessor master plan's section 12.2. That
    document now lives in docs/history/, so the canonical statement moved to
    .claude/README.md#authority-boundary (2026-08-28 merge). The assertion is
    kept, not dropped: an unsourced boundary becomes a rule nobody can argue
    with when circumstances change.
    """
    text = path.read_text(encoding="utf-8")
    assert "authority-boundary" in text, (
        "no citation of the authority boundary's basis "
        "(.claude/README.md#authority-boundary)")
    assert "inference" in text and "assumed" in text
    assert "Opus" in text, "does not say who to hand it to"


def test_extract_agent_cannot_write_files():
    """추출 에이전트에 Write 를 주면 spec 을 직접 고칠 수 있다."""
    tools = frontmatter(CLAUDE / "agents/bd-intake-extract.md")["tools"]
    assert "Write" not in tools


# =============================================================================
# 스킬이 물리를 다시 적지 않는다
# =============================================================================
def test_skills_delegate_to_the_core_rather_than_restating_physics():
    """The skill cites core functions and the CLI -- it does not produce numbers.

    The last assertion is the load-bearing one: without the prohibition stated
    in the skill, an agent that can do arithmetic will do arithmetic, and then
    a wrong result cannot be attributed to the physics or to the model.
    """
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "cli.py run" in text
    assert "simbot" in text
    assert "in your head" in text, "the prohibition on mental arithmetic is missing"


def test_pipeline_skill_forbids_confirming_verdicts():
    """Only a human fills `confirmed_by`, and the skill has to say so.

    Otherwise a pass gets stamped that nobody looked at, and the whole verdict
    chain becomes decoration.
    """
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "confirmed_by" in text
    assert "humans only" in text, "does not say that only a human may confirm"


def test_pipeline_skill_names_the_interpreter_absolutely():
    """`conda activate` is unreliable in a non-interactive shell."""
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python" in text


def test_pipeline_skill_states_the_question_budget():
    """If the conversation becomes twenty questions, the tool has failed."""
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    assert "Question budget" in text and "three per round" in text


def test_s1_reference_is_the_substantive_one():
    """S1 판독은 코드로 표현할 수 없는 유일한 단계 — 문서가 가장 두꺼워야 한다."""
    refs = {p.name: len(p.read_text(encoding="utf-8"))
            for p in (CLAUDE / "skills/bd-pipeline/references").glob("*.md")}
    assert refs["s1_intake_drawing.md"] > 3000


def test_s1_reference_warns_against_trusting_absolute_size():
    """The two rules that stop a sketch reading from inventing physics.

    People draw the box small and the particles large, so the drawn ratio
    overestimates phi almost every time; and an ambiguity resolved by picking one
    candidate silently is an assumption nobody can audit later.
    """
    text = (CLAUDE / "skills/bd-pipeline/references/s1_intake_drawing.md").read_text(
        encoding="utf-8")
    assert "Do not trust absolute sizes" in text
    assert "candidates" in text            # the ambiguity-enumeration convention


def test_diagnose_skill_suspects_analysis_first():
    """Of the 4 failures on the first end-to-end run, zero were physics.

    The ordering assertion is the real content: if `modeling` ever comes before
    `analysis` in this document, the elimination order has been inverted and the
    skill will send its reader chasing physics that is not there.
    """
    text = (CLAUDE / "skills/bd-diagnose/SKILL.md").read_text(encoding="utf-8")
    assert "analysis" in text
    assert "suspect the analysis code first" in text.lower()
    assert text.index("analysis") < text.index("modeling")


def test_knowledge_skill_states_the_citation_discipline():
    """Both halves of the discipline, or the wiki becomes a rumour store."""
    text = (CLAUDE / "skills/bd-knowledge/SKILL.md").read_text(encoding="utf-8")
    assert "not reproduced" in text, "the [source, not reproduced] marking is missing"
    assert "from memory" in text, "the ban on citing from memory is missing"


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
    """The reason the references went from 8 to 5 must stay recorded.

    A design that collapsed three stages into one CLI call is the kind of thing
    that looks like an oversight later. The count and the reason both have to be
    present, or someone re-splits them.
    """
    text = (CLAUDE / "README.md").read_text(encoding="utf-8")
    assert "split into five" in text, "the granularity decision is not recorded"
    assert "eight" in text, "does not say what it was reduced from"


def test_readme_lists_the_tiering_table():
    text = (CLAUDE / "README.md").read_text(encoding="utf-8")
    for name in EXPECTED_MODELS:
        assert name in text, name
