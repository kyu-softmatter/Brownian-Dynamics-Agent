"""Every repo-relative link in the documentation resolves.

★ Why this file exists. 39 links were broken across 17 files (measured
2026-09-15), and two of them were claims rather than navigation: a benchmark page
and its system card both cited `docs/tools/digitize_fig4b.py` as the script that
produced the digitized `g(r)` tolerances, and said the reading was reproducible
with it. That file, and the whole `docs/tools/` directory, has never existed in
this repository — so a verification tolerance rested on an artefact nobody can
open. Seven more still pointed at `inputs/`, renamed to `intake/` at the
2026-08-28 merge, one of them inside live report-generating code.

The registry of categories lives in `verify/verify_links.py` so it runs on its
own; this file is what makes CI fail.
"""
from __future__ import annotations

import importlib.util
import pathlib

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "verify_links", ROOT / "verify" / "verify_links.py")
VL = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(VL)


@pytest.fixture(scope="module")
def scanned():
    return VL.scan()


def test_the_scan_actually_sees_the_documentation(scanned):
    """Guard on the gate. A scan that finds nothing reports nothing broken, which
    is indistinguishable from a repository whose links all resolve."""
    tally, _ = scanned
    total = sum(tally.values())
    assert total >= VL.MIN_LINKS, f"only {total} links found"
    assert tally["ok"] > 0.8 * total, tally


def test_no_documentation_link_is_broken(scanned):
    tally, broken = scanned
    assert not broken, "\n".join(f"{f} -> {t}" for f, t, _ in broken)


def test_the_allowed_categories_are_each_small_and_accounted_for(scanned):
    """★ Three categories are allowed to be unresolvable, and each could be used
    to hide a real break. Their sizes are asserted so growth is visible.

    `unpublished` — per-run figures and trajectories, which `NOTICE.md` says are
    not published. `frozen` — the containing file is under seal, so the rename
    cannot be applied to it without breaking the seal; that is sealing working.
    `historical` — `docs/history/` holds copies of the predecessor repositories'
    documents, whose references point into repositories that no longer exist.
    """
    tally, _ = scanned
    assert tally["unpublished"] <= 40, tally["unpublished"]
    assert tally["frozen"] <= 5, tally["frozen"]
    assert tally["historical"] <= 8, tally["historical"]


def test_the_frozen_category_is_derived_from_the_seals_not_from_a_list():
    """A hand-written exemption list is what drifted everywhere else in this
    repository. `frozen` must come from the seal files themselves."""
    sealed = VL.sealed_documents()
    assert len(sealed) >= 60, len(sealed)          # 42 documents + 21 seals
    assert all("/" in s for s in sealed)
    assert any(s.endswith("SEALED.sha256") for s in sealed)
    assert any(s.endswith("01_intake.md") for s in sealed)


def test_docs_tools_is_still_absent_so_the_correction_stays_true():
    """Two pages now state that `docs/tools/digitize_fig4b.py` is not in the
    repository. If someone commits it, those corrections become false and should
    be reverted — so this fails rather than letting them rot in the other
    direction."""
    assert not (ROOT / "docs" / "tools").exists(), (
        "docs/tools/ exists now — the 'not reproducible' notes on "
        "knowledge/wiki/benchmarks/choi2020-interfacial-rdf.md and its system "
        "card should be revisited")


def test_each_allowed_link_really_belongs_to_the_category_it_claims():
    """★★ The category caps above are not enough, and a mutation proved it: a
    change that classified EVERY unresolvable link as `unpublished` kept all
    three caps satisfied and the suite passed. A category is only an exemption if
    membership in it is checkable, so each one is re-derived here independently.
    """
    import subprocess

    _, _, rows = VL.scan(with_rows=True)
    assert rows, "nothing unresolvable at all — the scan is not working"
    sealed = VL.sealed_documents()

    for kind, rel_file, target, rel_target in rows:
        if kind == "unpublished":
            rc = subprocess.run(["git", "check-ignore", "-q", rel_target],
                                cwd=ROOT).returncode
            assert rc == 0, (
                f"{rel_file} -> {target} is counted as deliberately unpublished, "
                f"but `git check-ignore` does not match {rel_target}")
        elif kind == "frozen":
            assert rel_file in sealed, (
                f"{rel_file} is counted as frozen by a seal but no seal covers it")
        elif kind == "historical":
            assert rel_file.startswith(VL.HISTORICAL_PREFIXES), rel_file


def test_every_unpublished_target_is_a_figure_or_a_trajectory():
    """Narrower still: the `.gitignore` rules that may excuse a link are the
    per-run figure and trajectory rules, not any rule at all. `knowledge/raw/`
    being ignored must not become a way to link at nothing."""
    _, _, rows = VL.scan(with_rows=True)
    allowed = (".png", ".gsd", ".npz", ".h5", ".parquet", ".pdf")
    bad = [(f, t) for kind, f, t, _ in rows
           if kind == "unpublished" and not t.endswith(allowed)]
    assert not bad, bad
