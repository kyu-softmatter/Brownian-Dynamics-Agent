"""Every count claimed in the documentation is re-measured, not maintained.

★ Why this file exists. The counts in `README.md` and `docs/` drifted three times
and were corrected by hand three times -- and the correcting commits themselves
introduced drift: `8bb504f` re-measured "every count in the README", changed 11 of
them, turned 42 paper distillations into 43 by counting a generated index, and
left `56/56` two lines below the line it edited. `aeb4a4b` fixed the tooling-entry
count 48 -> 57 and the next commit in the same session made it 58.

`CLAUDE.md` rule 10 records this repository's own tally: practices without gates
are three for three against. A count corrected by hand is a practice. This is the
gate.

The registry lives in `verify/verify_counts.py` so it is runnable on its own and
carries a `--fix`; this file is what makes CI fail.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "verify_counts", ROOT / "verify" / "verify_counts.py")
VC = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(VC)


@pytest.fixture(scope="module")
def measured() -> dict:
    return {name: fn() for name, fn in VC.COUNTERS.items()}


def test_every_claim_pattern_matches_exactly_once():
    """★ Guard on the gate, and the first thing to check.

    A regex that matches nothing reports no drift, which is indistinguishable
    from a claim that holds. Every pattern must occur exactly once in its file --
    a pattern that has rotted is an error, never a pass.
    """
    import re
    bad = []
    for claim in VC.CLAIMS:
        text = (ROOT / claim.file).read_text()
        n = len(re.findall(claim.pattern, text))
        if n != 1:
            bad.append(f"{claim.file}: matched {n}x — {claim.pattern}")
    assert not bad, "\n".join(bad)


def test_the_registry_covers_more_than_a_token_number_of_claims():
    """A registry can be emptied and this file would still pass. Assert the
    scale: the documentation carries ~95 checkable numbers across 8 files."""
    assert len(VC.CLAIMS) >= 50, len(VC.CLAIMS)
    files = {c.file for c in VC.CLAIMS}
    assert len(files) >= 7, sorted(files)
    assert "README.md" in files and "NOTICE.md" in files


@pytest.mark.parametrize("claim_index", range(len(VC.CLAIMS)),
                         ids=[f"{c.file}:{'+'.join(c.counters)}" for c in VC.CLAIMS])
def test_the_documented_count_is_what_is_measured(claim_index, measured):
    """One test per claim, so a failure names the file and the quantity.

    `$PY verify/verify_counts.py --fix` rewrites the numbers; do not hand-edit
    them, because that is the practice this gate replaces.
    """
    import re
    claim = VC.CLAIMS[claim_index]
    text = (ROOT / claim.file).read_text()
    hits = list(re.finditer(claim.pattern, text))
    assert len(hits) == 1, f"pattern matched {len(hits)}x"
    got = [int(g) for g in hits[0].groups()]

    if claim.kind == "sum":
        want = measured[claim.counters[0]]
        assert sum(got) == want, (
            f"{claim.file}: {' + '.join(map(str, got))} = {sum(got)}, "
            f"but {claim.counters[0]} measures {want}")
        return

    for value, name in zip(got, claim.counters):
        assert value == measured[name], (
            f"{claim.file}: {name} is claimed as {value}, measured {measured[name]}"
            f" — run `python verify/verify_counts.py --fix`")


def test_the_merge_record_is_deliberately_not_a_live_claim():
    """`docs/00-merge-decisions.md` records the three predecessor repositories as
    they stood at the 2026-08-28 merge -- "`bdbot` 21 modules, 8 cases, 278
    specs, 254 runs". Those are historical snapshots; re-measuring them would
    destroy the record. This asserts the intent, so nobody "fixes" them and
    nobody adds them to the registry by accident.
    """
    assert not any(c.file == "docs/00-merge-decisions.md" for c in VC.CLAIMS)
    text = (ROOT / "docs/00-merge-decisions.md").read_text()
    assert "21 modules" in text, "the merge-era snapshot was overwritten"
    assert "562 tests" in text, "the merge-era snapshot was overwritten"
