"""The `kind` vocabulary must keep reading the frozen archive.

`bdbot.checks` classifies a check as hard or soft from its `kind` string. That
vocabulary was translated Korean -> English in the code, but `specs/` and `runs/`
are frozen, content-addressed archives written by the older code, so almost every
archived check still carries the Korean spelling.

Measured on 2026-08-29, before `canon_kind` existed:

    2,350 archived checks, of which 2,333 carry a Korean `kind`
      738 archived SOFT checks classified as HARD
      390 of those were failing, so re-reading those specs flipped the verdict
          from "PASS (N warnings)" to FAIL and `health.gate()` blocked the run

That is the same false-rejection bug `bdbot.health.gate`'s own docstring records
(80 of 83 specs falsely rejected), reintroduced by a vocabulary change rather
than by a predicate change.

`run_id` was never affected, because `kind` is not part of the hash payload. That
is exactly what made it easy to miss: **"not hash-covered" and "not load-bearing"
are different properties**, and only the first was checked when the vocabulary
moved.

These tests are the wiring. An unwired fix is the failure this repository keeps
finding, so they run in the normal suite and do not depend on any particular run
or spec existing.
"""
from __future__ import annotations

import glob
import json

from bdbot import checks as C

# The five kinds the older code wrote into the archive, and what each means now.
LEGACY_TO_CURRENT = {
    "통계": "statistics",
    "적분": "integration",
    "모델": "model",
    "기하": "geometry",
    "유한크기": "finite-size",
}


def test_canon_kind_maps_every_legacy_spelling():
    for legacy, current in LEGACY_TO_CURRENT.items():
        assert C.canon_kind(legacy) == current, legacy


def test_canon_kind_passes_current_spellings_through_unchanged():
    for current in LEGACY_TO_CURRENT.values():
        assert C.canon_kind(current) == current, current


def test_canon_kind_passes_an_unknown_kind_through_rather_than_dropping_it():
    # A kind nobody has seen must survive, not become None or "".
    assert C.canon_kind("something-new") == "something-new"


def test_legacy_soft_kinds_are_still_soft():
    """The regression itself: a Korean soft kind must not read as hard."""
    for legacy in ("통계", "유한크기"):
        assert C.soft(legacy) is True, (
            f"{legacy!r} classified as HARD -- re-reading an archived spec would "
            f"flip its verdict to FAIL and gate() would block the run")


def test_legacy_hard_kinds_are_still_hard():
    """The opposite direction has to hold too, or the fix would pass everything."""
    for legacy in ("적분", "모델", "기하"):
        assert C.soft(legacy) is False, legacy


def test_current_vocabulary_classification_is_unchanged():
    assert C.soft("statistics") is True
    assert C.soft("finite-size") is True
    assert C.soft("integration") is False
    assert C.soft("model") is False
    assert C.soft("geometry") is False


def test_every_kind_in_the_archive_is_recognised():
    """No archived `kind` may fall through `canon_kind` unrecognised.

    This is the check that would catch a *sixth* legacy spelling nobody mapped.
    It is skipped rather than failed when specs/ is absent, so a fresh clone
    without the archive does not see a spurious failure.
    """
    known = set(LEGACY_TO_CURRENT) | set(LEGACY_TO_CURRENT.values())
    seen: set[str] = set()
    for path in glob.glob("specs/*.json"):
        try:
            raw = json.loads(open(path, encoding="utf-8").read())
        except Exception:
            continue
        for ch in raw.get("checks", []):
            k = ch.get("kind")
            if isinstance(k, str):
                seen.add(k)
    if not seen:
        import pytest
        pytest.skip("no specs/ archive present")
    unknown = seen - known
    assert not unknown, (
        f"archived kind(s) not in the vocabulary: {sorted(unknown)} -- add them to "
        f"bdbot.checks.LEGACY_KINDS or they will be misclassified")


def test_no_archived_soft_check_reads_as_hard():
    """The end-to-end assertion, measured against the real archive.

    Before the fix this counted 738. It must be 0.
    """
    misclassified = 0
    for path in glob.glob("specs/*.json"):
        try:
            raw = json.loads(open(path, encoding="utf-8").read())
        except Exception:
            continue
        for ch in raw.get("checks", []):
            k = ch.get("kind")
            if C.canon_kind(k) in C.SOFT_KINDS and not C.soft(k):
                misclassified += 1
    assert misclassified == 0, (
        f"{misclassified} archived soft checks classify as hard")
