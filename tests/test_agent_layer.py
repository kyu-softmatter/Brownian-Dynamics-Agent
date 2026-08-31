"""Structural checks on the L1 agent layer (`.claude/`).

## Why markdown gets tested

This layer is documentation, not code. But **the frontmatter is a contract a machine
reads**, a broken link means the agent cannot find the protocol, and a missing
tiering boundary means a cheap model makes a physics judgment. All of them fail
silently.

## What is checked

1. valid frontmatter (`name`, `description`, `model`), with the name matching the
   filename
2. relative links actually resolve
3. **model tiering** — matches the §12.2 allocation table, with permission
   boundaries stated for the cheap models
4. **a skill does not re-write the physics** — does it cite the core function
5. does `settings.json` refuse edits to a sealed document
"""
from __future__ import annotations

import json
import re
from pathlib import Path

import pytest
import yaml

CLAUDE = Path(".claude")
pytestmark = pytest.mark.skipif(not CLAUDE.exists(), reason=".claude/ is absent")

SKILLS = sorted(CLAUDE.glob("skills/*/SKILL.md"))
AGENTS = sorted(CLAUDE.glob("agents/*.md"))

# The 2026-08-28 merge split the skills into two kinds, and they must not be measured
# by one rule. The orchestration skills are **mutually exclusive** -- trigger the
# wrong one and an entirely different procedure runs -- so cross-references are
# mandatory. The domain skills are references the pipeline *pulls in* at each stage,
# so they are not exclusive; instead the pipeline has to point at them.
ORCH_SKILLS = {"bd-pipeline", "bd-diagnose", "bd-knowledge"}
DOMAIN_SKILLS = {"bd-intake", "bd-physics", "bd-hoomd"}

# The master_plan §12.1 allocation table. A departure means the design and the code
# have drifted apart.
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
    assert m, f"{path}: no frontmatter — it is a contract a machine reads"
    return yaml.safe_load(m.group(1))


# =============================================================================
# Existence
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
    assert fm["name"] == path.parent.name, \
        "the name differs from the directory name, so it never gets invoked"
    assert len(fm["description"]) > 40, \
        "too short a description does not trigger"


@pytest.mark.parametrize("path", AGENTS, ids=lambda p: p.stem)
def test_agent_frontmatter(path):
    fm = frontmatter(path)
    assert fm["name"] == path.stem
    assert fm["model"] in ("opus", "sonnet", "haiku", "inherit")
    assert len(fm["description"]) > 40


@pytest.mark.parametrize("path", AGENTS, ids=lambda p: p.stem)
def test_agent_model_matches_the_tiering_table(path):
    """★ Does it match the master_plan §12.1 allocation table.

    A mismatch means the design document and the actual delegation have drifted --
    and which one is right has to be decided.
    """
    fm = frontmatter(path)
    assert fm["model"] == EXPECTED_MODELS[path.stem]


def test_orch_skill_descriptions_say_when_not_to_use():
    """Mutually exclusive skills have to point at each other, or the wrong one gets
    triggered."""
    for path in SKILLS:
        if path.parent.name not in ORCH_SKILLS:
            continue
        desc = frontmatter(path)["description"]
        others = ORCH_SKILLS - {path.parent.name}
        assert any(o in desc for o in others), \
            f"{path.parent.name}: no mention of the other skills"


def test_pipeline_points_at_every_domain_skill():
    """A domain skill does not trigger itself — it only gets read when the pipeline
    pulls it in.

    Right after the merge this was actually empty: S5 writes HOOMD code and nothing
    pointed at bd-hoomd, which holds 20 traps. This check guards that gap.
    """
    body = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    missing = sorted(s for s in DOMAIN_SKILLS if s not in body)
    assert not missing, f"domain skills bd-pipeline does not point at: {missing}"


# =============================================================================
# Links
# =============================================================================
def _links(path: Path) -> list[str]:
    text = path.read_text(encoding="utf-8")
    return [m.group(1).split("#")[0]
            for m in re.finditer(r"\[[^\]]+\]\((\.\.?/[^)#]*)\)", text)]


@pytest.mark.parametrize("path", sorted(CLAUDE.rglob("*.md")),
                         ids=lambda p: str(p.relative_to(CLAUDE)))
def test_relative_links_resolve(path):
    """A broken link means the agent cannot find the protocol."""
    for link in _links(path):
        assert (path.parent / link).resolve().exists(), f"{path}: {link}"


def test_pipeline_skill_links_every_reference():
    text = (CLAUDE / "skills/bd-pipeline/SKILL.md").read_text(encoding="utf-8")
    for ref in (CLAUDE / "skills/bd-pipeline/references").glob("*.md"):
        assert ref.name in text, f"SKILL.md does not cite {ref.name}"


# =============================================================================
# The tiering safeguard — keeping a cheap model out of physics judgments
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
    """Give the extraction agent Write and it can edit the spec directly."""
    tools = frontmatter(CLAUDE / "agents/bd-intake-extract.md")["tools"]
    assert "Write" not in tools


# =============================================================================
# A skill does not re-write the physics
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
    """The S1 reading is the one stage that cannot be expressed as code — its
    documentation has to be the thickest."""
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
    """CLAUDE.md: it is unreliable in a non-interactive shell."""
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    assert any("conda activate" in d for d in s["permissions"]["deny"])


def test_settings_denies_editing_sealed_documents():
    """★ Structurally closes the easiest route to post-hoc rationalisation.

    Seal verification catches it after the fact, but preventing it in the first place
    is better.
    """
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    deny = " ".join(s["permissions"]["deny"])
    for sealed in ("02_prediction.md", "01_intake.md", "SEALED.sha256"):
        assert sealed in deny, f"editing {sealed} is permitted"


def test_settings_does_not_deny_reading_inputs():
    """S1 has to read the hand sketch — block this and the pipeline dies."""
    s = json.loads((CLAUDE / "settings.json").read_text(encoding="utf-8"))
    deny = " ".join(s["permissions"]["deny"])
    assert "Read(./inputs" not in deny
    assert "Read(./knowledge" not in deny


# =============================================================================
# Documented decisions
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
