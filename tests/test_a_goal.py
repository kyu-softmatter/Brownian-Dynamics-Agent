"""Block A — `bdbot/goal.py`, the goal artefact.

Named `test_a_goal.py` rather than `test_s1_*`: block A runs *before* S1, because
the goal is fixed before the sketch is read. The filename is the pipeline order.

Structure mirrors `tests/test_s1_intake_goal.py`: does it fire, is it wired,
does turning it off change the verdict, and is its scope narrow on purpose.
"""
from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from bdbot import goal as G
from bdbot import intake as I

ROOT = Path(__file__).resolve().parent.parent

VALID = {
    "schema": G.SCHEMA,
    "case": "synthetic",
    "asked_by": "experimentalist",
    "question": "Does the bead's PSD corner frequency recover the trap stiffness?",
    "answering_quantity": {"symbol": "f_c", "unit": "Hz", "why": "f_c = k/(2 pi gamma)"},
    "decisive_precision": {
        "value": 5.0, "basis": "the 5 % tolerance this case already carries",
        "what_would_change_my_mind": "f_c off by more than 5 %, or a non-Gaussian histogram"},
    "physics_that_matters": [{"symbol": "k_t", "why": "sets f_c directly"}],
    "analysis_implied": ["PSD with a Lorentzian fit"],
    "provenance": {"source": "conversation 2026-09-02", "confirmed_by": None},
}


@pytest.fixture
def write_goal(tmp_path):
    def go(raw, name="goal.yaml"):
        if raw is None:                       # write nothing -> ABSENT
            return G.load(tmp_path)
        (tmp_path / name).write_text(yaml.safe_dump(raw, allow_unicode=True))
        return G.load(tmp_path)
    go.dir = tmp_path
    return go


# -- 1. the happy path, and the DRAFT/CONFIRMED distinction ------------------

def test_a_valid_goal_has_no_errors_and_is_known(write_goal):
    g = write_goal(VALID)
    assert g.errors == [], [str(i) for i in g.errors]
    assert g.known is True
    assert g.confirmed is False               # confirmed_by is null
    assert G.status(write_goal.dir)[0] == "DRAFT"
    assert G.blocks(write_goal.dir) is False, "a DRAFT must not block"


def test_confirmed_by_flips_draft_to_confirmed(write_goal):
    raw = copy.deepcopy(VALID)
    raw["provenance"]["confirmed_by"] = "Kyu Hwan 2026-09-02"
    g = write_goal(raw)
    assert g.confirmed is True
    assert G.status(write_goal.dir)[0] == "CONFIRMED"


def test_a_missing_file_is_absent_and_blocks(write_goal):
    g = write_goal(None)
    assert g.exists is False
    assert g.known is False
    verdict, reason = G.status(write_goal.dir)
    assert verdict == "ABSENT"
    assert "goal init" in g.errors[0].msg     # tells you how to fix it
    assert G.blocks(write_goal.dir) is True


# -- 2. every required field actually required -------------------------------

@pytest.mark.parametrize("field", G.REQUIRED_TOP)
def test_each_required_field_is_required(write_goal, field):
    raw = copy.deepcopy(VALID)
    raw.pop(field)
    g = write_goal(raw)
    assert any(i.where == field and i.level == "error" for i in g.issues), \
        f"dropping {field!r} produced no error: {[str(i) for i in g.issues]}"


@pytest.mark.parametrize("question, ok", [
    ("Does the bead's PSD recover the trap stiffness?", True),
    ("", False),
    ("   ", False),
    ("MSD", False),                            # a label, not a question
    ("measure MSD, MSAD", True),               # the shortest real one, 17 chars
])
def test_the_question_must_be_a_question(write_goal, question, ok):
    raw = copy.deepcopy(VALID)
    raw["question"] = question
    g = write_goal(raw)
    has_err = any(i.where == "question" and i.level == "error" for i in g.issues)
    assert has_err is (not ok), [str(i) for i in g.issues]


def test_the_falsifiability_gate_is_enforced(write_goal):
    """★ `what_would_change_my_mind` is the gate S1 has always stated in prose.
    This is the first time it is checked by code anywhere in the repository.
    """
    raw = copy.deepcopy(VALID)
    raw["decisive_precision"]["what_would_change_my_mind"] = ""
    g = write_goal(raw)
    assert any("what_would_change_my_mind" in i.where and i.level == "error"
               for i in g.issues), [str(i) for i in g.issues]
    assert G.status(write_goal.dir)[0] == "FAIL"


def test_a_null_precision_is_allowed_but_needs_a_basis(write_goal):
    """Often the decisive precision is unknown until the design exists. That is
    permitted -- silently omitting the reason is not."""
    raw = copy.deepcopy(VALID)
    raw["decisive_precision"]["value"] = None
    raw["decisive_precision"]["basis"] = "not known until the trap stiffness is fixed"
    assert write_goal(raw).errors == []

    raw["decisive_precision"]["basis"] = ""
    g = write_goal(raw)
    assert any(i.where == "decisive_precision" and i.level == "error"
               for i in g.issues), [str(i) for i in g.issues]


@pytest.mark.parametrize("key", ["physics_that_matters", "analysis_implied"])
def test_the_two_lists_may_not_be_empty(write_goal, key):
    """Empty defeats the point: `physics_that_matters` is why stating the goal
    first helps at all, and `analysis_implied` is how the goal reaches the end."""
    raw = copy.deepcopy(VALID)
    raw[key] = []
    g = write_goal(raw)
    assert any(i.where == key and i.level == "error" for i in g.issues)


def test_physics_items_need_a_reason_not_just_a_symbol(write_goal):
    raw = copy.deepcopy(VALID)
    raw["physics_that_matters"] = [{"symbol": "k_t", "why": ""}]
    g = write_goal(raw)
    assert any(i.level == "error" and "why" in i.msg for i in g.issues)


# -- 3. ★ break it deliberately -----------------------------------------------

def test_disabling_the_gate_restores_a_clean_goal(write_goal, monkeypatch):
    """CLAUDE.md: *"When you build a checker, deliberately break it and see."*
    Forcing `validate` to return nothing must flip FAIL to DRAFT -- otherwise
    something other than this validator is producing the verdict.
    """
    raw = copy.deepcopy(VALID)
    raw["question"] = ""                       # a real, detectable error
    assert G.status(write_goal.dir) or True    # (dir has no file yet)
    g = write_goal(raw)
    assert g.errors, "precondition: the validator is catching the empty question"
    assert G.status(write_goal.dir)[0] == "FAIL"

    monkeypatch.setattr(G, "validate", lambda g: [])
    assert G.status(write_goal.dir)[0] == "DRAFT", (
        "disabling validate() did not change the verdict -- status() is not "
        "actually driven by the validator")


def test_status_is_the_single_source_of_the_verdict(write_goal):
    """`blocks()` must be derivable from `status()` and never drift from it."""
    for raw, expect in [(VALID, False),
                        ({**VALID, "question": ""}, True),
                        (None, True)]:
        write_goal(raw) if raw is not None else None
        if raw is None:
            import shutil
            shutil.rmtree(write_goal.dir)
            write_goal.dir.mkdir()
        assert G.blocks(write_goal.dir) is expect
        assert (G.status(write_goal.dir)[0] in ("ABSENT", "FAIL")) is expect


# -- 4. narrow on purpose: the sketch may legitimately say nothing -----------

def test_an_empty_stated_goals_no_longer_blocks_l0():
    """★ The correction. `stated_goals: []` records that the SKETCH is silent,
    which rule 5 requires. Blocking on it punished honest transcription; the
    goal now lives in goal.yaml instead.

    Regression lock against re-adding the block to `intake`.
    """
    obs = I.load(ROOT / "intake/trap-2d-5um")
    assert obs.open_goal is True, "this sketch does state no goal"
    assert obs.ready_for_system is True, (
        "a silent sketch must not block L0 -- see bdbot/intake.py "
        "GOAL_ABSENT_NOTE for why this was reverted")
    assert "stated_goals" not in obs.blockers


def test_intake_still_reports_the_silent_sketch():
    """Not blocking is not the same as not mentioning. The fact is still surfaced,
    with a pointer to where the goal belongs."""
    obs = I.load(ROOT / "intake/trap-2d-5um")
    text = I.render_check(obs)
    assert "STATES NO GOAL" in text
    assert "goal.yaml" in text
    assert "VERDICT: READY" in text


# -- 5. the real tree --------------------------------------------------------

def test_trap_2d_5um_has_a_draft_goal():
    """Backfilled 2026-09-02 by transcription from bd-intake section 2.1.
    It is DRAFT until a human sets `confirmed_by`, and DRAFT does not block.
    """
    d = ROOT / "intake/trap-2d-5um"
    g = G.load(d)
    assert g.exists and g.errors == [], [str(i) for i in g.errors]
    assert g.answering_symbol == "f_c"
    assert G.status(d)[0] == "DRAFT"
    assert G.blocks(d) is False


def test_every_case_folder_is_accounted_for():
    """Each of the 8 cases is ABSENT, DRAFT, FAIL or CONFIRMED -- never something
    else. Guards against `status()` growing a fourth state nobody handles."""
    cases = sorted(p.parent for p in ROOT.glob("intake/*/observation.yaml"))
    assert len(cases) == 8
    for d in cases:
        v, reason = G.status(d, I.load(d))
        assert v in ("ABSENT", "DRAFT", "FAIL", "CONFIRMED"), (d.name, v)
        assert reason, f"{d.name}: a verdict with no reason"


def test_the_template_it_writes_passes_its_own_required_field_check(tmp_path):
    """A template that fails its own schema teaches the wrong shape. It must have
    every required KEY present -- while still failing on empty CONTENT, because an
    untouched template is not a goal.
    """
    made, _ = G.init_template(tmp_path, case="synthetic")
    assert made
    g = G.load(tmp_path)
    missing_key = [i for i in g.errors if i.msg == "required field missing"]
    assert missing_key == [], [str(i) for i in missing_key]
    assert g.errors, "an untouched template must not pass as a real goal"
