"""S1 / L0 — what `stated_goals` does and does not block.

**This file records a reversal, so read the history before changing it.**

  2026-09-02 (morning)  `stated_goals: []` was made a BLOCKER here. The rule had
                        been written in `bd-intake` SKILL.md §2.1 and in
                        `intake.py`'s own TEMPLATE, and enforced in neither, while
                        2 of 8 real intakes carried an empty list and had run
                        anyway (4 and 81 run directories).

  2026-09-02 (later)    ★ Reverted. `stated_goals` records what the SKETCH says,
                        and rule 5 requires that anything absent from the sketch
                        stay null — so an empty list on a silent sketch is a
                        *correct transcription*, and blocking on it punished
                        rule 5. `trap-2d-5um` is exactly that case: `bd-intake`
                        §2.1 records its goal as "absent — settled by asking the
                        user", so the goal existed with nowhere to live.
                        The fix is block A, `bdbot/goal.py` + `tests/test_a_goal.py`.

So what is asserted here now is the **narrow** claim: a silent sketch is reported
and does not block, physical gaps still do, and the two verdicts stay separate.
The blocking half lives in `test_a_goal.py`.

Still the first pytest coverage `bdbot/intake.py` has had — before this the only
checks on the L0 validator were in `verify/verify_intake_guards.py`, which CI does
not run.
"""
from __future__ import annotations

import copy
from pathlib import Path

import pytest
import yaml

from bdbot import goal as G
from bdbot import intake as I
from bdbot import physical as P

ROOT = Path(__file__).resolve().parent.parent
DONOR = ROOT / "intake/soft-r3-2d-A-sweep"

# The two cases whose sketch states no goal. Not blocked for it — see the header.
SILENT_SKETCH_CASES = ("trap-2d-5um", "trap-drag-2d-hex300")


@pytest.fixture
def donor_obs() -> dict:
    return yaml.safe_load((DONOR / "observation.yaml").read_text())


@pytest.fixture
def write_obs(tmp_path):
    """Synthetic input on purpose. `verify_intake_guards.py` records why: an
    adversarial test pinned to a real case silently disables itself the moment
    that case is resolved, and `abp-rod` actually did that.
    """
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")

    def go(raw: dict) -> I.Observation:
        (tmp_path / "observation.yaml").write_text(
            yaml.safe_dump(raw, allow_unicode=True))
        return I.load(tmp_path)

    go.dir = tmp_path
    return go


# -- 1. open_goal detects a silent sketch, and that is all it does -----------

@pytest.mark.parametrize("goals, silent", [
    (["measure the PSD and fit f_c"], False),
    (["a", "b"], False),
    ([], True),
    (None, True),
])
def test_open_goal_detects_a_silent_sketch(donor_obs, write_obs, goals, silent):
    raw = copy.deepcopy(donor_obs)
    raw["stated_goals"] = goals
    obs = write_obs(raw)
    assert obs.open_goal is silent
    assert not obs.errors, [str(i) for i in obs.errors]


def test_a_silent_sketch_does_not_block_l0(donor_obs, write_obs):
    """★ The reversal, as a regression lock. Do not re-add this block here."""
    raw = copy.deepcopy(donor_obs)
    raw["stated_goals"] = []
    obs = write_obs(raw)

    assert obs.open_goal is True
    assert obs.ready_for_system is True
    assert obs.blockers == []
    assert "VERDICT: READY" in I.render_check(obs)


def test_a_silent_sketch_is_still_reported(donor_obs, write_obs):
    """Not blocking is not the same as not mentioning. The fact stays visible and
    points at where the goal belongs — otherwise reverting the block would have
    thrown away the finding along with it."""
    raw = copy.deepcopy(donor_obs)
    raw["stated_goals"] = []
    text = I.render_check(write_obs(raw))
    assert "STATES NO GOAL" in text
    assert "goal.yaml" in text
    assert "rule 5" in text


def test_a_missing_key_is_still_a_schema_error(donor_obs, write_obs):
    """`stated_goals` is in REQUIRED_TOP. Absent-key and present-but-empty are
    different: the first is dishonest bookkeeping, the second can be the truth."""
    raw = copy.deepcopy(donor_obs)
    raw.pop("stated_goals")
    obs = write_obs(raw)
    assert any(i.where == "stated_goals" and i.level == "error" for i in obs.issues)


# -- 2. blockers stays the single source ------------------------------------

def test_blockers_is_the_single_source_of_the_blocking_set(donor_obs, write_obs):
    """Regression lock on the bug this work started from: `ready_for_system` had
    zero call sites while four consumers each re-derived the blocking set from
    `open_missing` alone, so a blocker could be added to one and not the other.
    """
    for goals in ([], None, ["real goal"]):
        raw = copy.deepcopy(donor_obs)
        raw["stated_goals"] = goals
        obs = write_obs(raw)
        assert obs.ready_for_system == (not obs.blockers)


def test_a_physical_gap_does_block(donor_obs, write_obs):
    """The thing L0 *is* allowed to block on, so the tests above cannot pass by
    the blocking mechanism being broken outright."""
    raw = copy.deepcopy(donor_obs)
    raw["missing_required"].append({
        "symbol": "made_up_param", "kind": "physical",
        "what": "synthetic unresolved physical gap",
        "assumed_value": None, "resolution": None})
    obs = write_obs(raw)
    assert obs.ready_for_system is False
    assert obs.blockers == ["made_up_param"]


# -- 3. the two verdicts are separate, and L2 survives a goal-only problem ---

def test_l0_and_goal_verdicts_are_independent(donor_obs, tmp_path):
    """A silent sketch with no goal.yaml: L0 READY, goal ABSENT. `bdbot.cli status`
    shows both columns; neither answers for the other."""
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    raw = copy.deepcopy(donor_obs)
    raw["stated_goals"] = []
    (tmp_path / "observation.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))

    obs = I.load(tmp_path)
    assert obs.ready_for_system is True
    assert G.status(tmp_path, obs)[0] == "ABSENT"
    assert G.blocks(tmp_path, obs) is True


def test_a_goal_problem_does_not_invalidate_a_settled_l2(donor_obs, tmp_path):
    """`physical.validate()` errors when L0 has unresolved *physical* gaps, on the
    grounds that *"a value may have been invented somewhere (rule 3)."* A goal
    problem implies no invented value, so it must not cascade — 85 existing runs
    across the two silent-sketch cases depend on their L2 staying valid.
    """
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    raw = copy.deepcopy(donor_obs)
    raw["stated_goals"] = []
    (tmp_path / "observation.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))
    (tmp_path / "system.yaml").write_text((DONOR / "system.yaml").read_text())

    s = P.load(tmp_path)
    assert not s.errors, [str(i) for i in s.errors]


def test_a_physical_gap_still_does_invalidate_l2(donor_obs, tmp_path):
    """The counterpart. Also re-asserts the invariant whose own assertion in
    `verify_intake_guards.py` was filtering for a Korean string that
    `physical.py`'s English message never contained — it reported "not caught"
    while the checker was catching it (pre-existing at 4b4503a, fixed 2026-09-02).
    """
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    raw = copy.deepcopy(donor_obs)
    raw["missing_required"].append({
        "symbol": "made_up_param", "kind": "physical",
        "what": "synthetic unresolved physical gap",
        "assumed_value": None, "resolution": None})
    (tmp_path / "observation.yaml").write_text(yaml.safe_dump(raw, allow_unicode=True))
    (tmp_path / "system.yaml").write_text((DONOR / "system.yaml").read_text())

    s = P.load(tmp_path)
    assert [i for i in s.errors
            if i.where == "derived_from" and "made_up_param" in i.msg], \
        [str(i) for i in s.errors]


# -- 4. the real tree -------------------------------------------------------

@pytest.mark.parametrize("case", SILENT_SKETCH_CASES)
def test_the_two_silent_sketches_are_ready_at_l0(case):
    obs = I.load(ROOT / "intake" / case)
    assert obs.open_goal is True, "this sketch does state no goal"
    assert obs.ready_for_system is True, "a silent sketch must not block L0"
    assert not obs.errors, [str(i) for i in obs.errors]


def test_all_eight_cases_are_ready_at_l0():
    """L0 was READY for all 8 at HEAD and must stay that way — this whole change
    should have moved nothing at L0."""
    cases = sorted(p.parent for p in ROOT.glob("intake/*/observation.yaml"))
    assert len(cases) == 8
    for d in cases:
        obs = I.load(d)
        assert obs.ready_for_system, f"{d.name}: {obs.blockers}"
