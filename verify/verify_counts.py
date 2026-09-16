"""Every count claimed in the documentation, re-measured.

    $PY verify/verify_counts.py          # check; non-zero exit on any drift
    $PY verify/verify_counts.py --fix    # rewrite the claimed numbers in place
    $PY verify/verify_counts.py -v       # show every claim, not only the drifted

**Why this file exists.** The counts in `README.md` and `docs/` have drifted
three times and been corrected by hand three times:

  · 2026-09-15 `8bb504f` re-measured "every count in the README" and changed 11
    of them -- and introduced one (42 -> 43 paper distillations, by counting the
    generated `INDEX.md`) while leaving `56/56` untouched two lines below the one
    it edited
  · 2026-09-15 `aeb4a4b` corrected the tooling-entry count 48 -> 57, and the
    next commit in the same session made it 58 by adding an entry
  · the architecture diagrams still carried the 2026-08-28 merge-era numbers, so
    `README.md` stated the entry count as both 145 and 126

This repository's own record on practices-without-gates is documented as three
for three against (`CLAUDE.md` rule 10). A count corrected by hand is a
practice. So the counts are now derived, and a drift is a test failure.

**What makes the number mean anything.** A regex that matches nothing reports no
drift, which is indistinguishable from a claim that holds -- `docs/05-pitfalls.md`
section 2, *"a check whose success is indistinguishable from its own failure."*
So every claim asserts its pattern occurs **exactly once** in its file before any
comparison, and a pattern that does not is an **ERROR**, never a pass. That is
the same guard `verify/verify_gates_bite.py` puts on its anchors.

**What is deliberately NOT in here.** `docs/00-merge-decisions.md` records the
three predecessor repositories as they stood at the 2026-08-28 merge -- "`bdbot`
21 modules, 8 cases, 278 specs, 254 runs", "`Simulation_bot` 19 modules, 562
tests". Those are historical snapshots, not live claims, and updating them would
destroy the record. Anything absent from `CLAIMS` is by construction not a live
claim; if you add a count to the docs, add it here.
"""
from __future__ import annotations

import argparse
import collections
import json
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
PY = sys.executable


# ── the counters: one rule each, applied uniformly ─────────────────────────
#
# ⚠ Two rules that are easy to get wrong, and did get wrong:
#   * modules are counted RECURSIVELY. `simbot/` has an `analysis/` subpackage,
#     so `simbot/*.py` is 16 while the documented (and correct) number is 19.
#   * the artefact counts are of what a FRESH CLONE has, i.e. git-tracked. Two
#     `runs/` directories here hold only `observables.npz`/`traj_A.gsd` with no
#     `metrics.json` and no `record.json`; they are not in the ledger and not in
#     the repository, so the run-directory count is 273, not the 275 on this disk.

def _files(pattern: str) -> int:
    return len([p for p in ROOT.glob(pattern) if p.is_file()])


def _modules(package: str) -> int:
    return len([p for p in (ROOT / package).rglob("*.py")
                if "__pycache__" not in p.parts])


def _tracked(pattern: str) -> list[str]:
    out = subprocess.run(["git", "ls-files", pattern], cwd=ROOT,
                         capture_output=True, text=True, check=True).stdout
    return [line for line in out.split("\n") if line.strip()]


def _tracked_run_dirs() -> int:
    return len({p.split("/")[1] for p in _tracked("runs/*/*") if p.count("/") >= 2})


def _papers() -> int:
    """`INDEX.md` is generated -- its own header says so, and counts 40 of the 42
    files beside it. `bd-knowledge/SKILL.md` writes the basis down as
    "42 per-paper distillations + INDEX.md", i.e. the index sits alongside the
    count, not inside it. Counting it is what turned 42 into 43."""
    return _files("knowledge/source/papers/*.md") - 1


def _entries_by_origin() -> collections.Counter:
    c = collections.Counter()
    for p in (ROOT / "knowledge/entries").glob("*.json"):
        c[json.loads(p.read_text()).get("origin")] += 1
    return c


def _wiki_pages() -> int:
    return len([p for p in (ROOT / "knowledge/wiki").rglob("*.md")])


def _lab_authored() -> int:
    return sum(1 for p in (ROOT / "knowledge/source/papers").glob("*.md")
               if p.name != "INDEX.md" and "lab_authored: true" in p.read_text())


def _papers_field(field: str, *, value: str | None = None, absent: bool = False) -> int:
    """Count distillations by a frontmatter field. `absent=True` counts the ones
    that do not carry the field at all, which `NOTICE.md` treats as operationally
    the same as `verified: false`."""
    n = 0
    for p in (ROOT / "knowledge/source/papers").glob("*.md"):
        if p.name == "INDEX.md":
            continue
        text = p.read_text()
        head = text.split("---")[1] if text.startswith("---") else ""
        m = re.search(rf"^{field}:[ \t]*(\S.*)$", head, re.M)
        got = m.group(1).strip() if m else None
        if got in ("null", '""', "''", "~"):
            got = None
        if absent:
            n += got is None
        elif value is None:
            n += got is not None
        else:
            n += got == value
    return n


def _papers_with_doi() -> int:
    n = 0
    for p in (ROOT / "knowledge/source/papers").glob("*.md"):
        if p.name == "INDEX.md":
            continue
        text = p.read_text()
        head = text.split("---")[1] if text.startswith("---") else ""
        m = re.search(r"^doi:[ \t]*(\S.*)$", head, re.M)
        if m and m.group(1).strip() not in ("null", '""', "''", "~"):
            n += 1
    return n


def _tests_collected() -> int:
    """`pytest --collect-only -q`, 1.7 s. The suite's own `N passed, M skipped`
    cannot be derived without running it -- skips are decided at runtime -- so the
    claims are checked as a SUM against this instead. That catches the drift that
    actually happens (the suite grew) without a second full run."""
    out = subprocess.run([PY, "-m", "pytest", "--collect-only", "-q"], cwd=ROOT,
                         capture_output=True, text=True).stdout
    m = re.search(r"(\d+) tests? collected", out)
    if not m:
        raise RuntimeError(f"could not read the collected count from:\n{out[-800:]}")
    return int(m.group(1))


def _book_claims() -> int:
    """`verify/verify_book_claims.py`, 0.4 s. The claim was written as 56/56 and
    the script has printed 64/64 since 2026-08-29 (`4e52fa0`); the +8 is 2 added
    `[OURS]` and 6 added `[DOC]` checks, the latter added precisely because a
    transcription typo said "passed 56/56"."""
    out = subprocess.run([PY, "verify/verify_book_claims.py"], cwd=ROOT,
                         capture_output=True, text=True).stdout
    hits = re.findall(r"(\d+)/(\d+) PASS", out)
    if not hits:
        raise RuntimeError(f"no 'N/N PASS' line in:\n{out[-800:]}")
    total, of = hits[-1]
    if total != of:
        raise RuntimeError(f"verify_book_claims.py is failing: {total}/{of}")
    return int(of)


def _sealed_documents() -> int:
    return sum(len([q for q in p.read_text().splitlines() if q.strip()])
               for p in ROOT.glob("runs*/*/SEALED.sha256"))


def _index_stated() -> int:
    """What the generated `knowledge/source/papers/INDEX.md` says about itself.
    The README cites the gap between this and the file count as a finding of the
    librarian agent, so both halves of that sentence have to stay true."""
    m = re.search(r"항목 수 \| \*\*(\d+)\*\*",
                  (ROOT / "knowledge/source/papers/INDEX.md").read_text())
    if not m:
        raise RuntimeError("INDEX.md no longer states its own item count")
    return int(m.group(1))


COUNTERS: dict[str, callable] = {
    "specs":            lambda: _files("specs/*.json"),
    "run_dirs":         _tracked_run_dirs,
    "metrics":          lambda: len(_tracked("runs/*/metrics.json")),
    "postmortems":      lambda: len(_tracked("runs/*/record.json")),
    "verify_scripts":   lambda: _files("verify/*.py"),
    "bdbot_modules":    lambda: _modules("bdbot"),
    "simbot_modules":   lambda: _modules("simbot"),
    "case_scripts":     lambda: _files("cases/*.py"),
    "campaign_scripts": lambda: _files("campaigns/*.py"),
    "tests_collected":  _tests_collected,
    "wiki_pages":       _wiki_pages,
    "findings":         lambda: _files("knowledge/wiki/findings/*.md"),
    "concepts":         lambda: _files("knowledge/wiki/concepts/*.md"),
    "system_cards":     lambda: _files("knowledge/wiki/systems/*.md"),
    "benchmarks":       lambda: _files("knowledge/wiki/benchmarks/*.md"),
    "papers":           _papers,
    "books":            lambda: _files("knowledge/source/books/*.md"),
    "distillations":    lambda: _papers() + _files("knowledge/source/books/*.md"),
    "papers_with_doi":  _papers_with_doi,
    "papers_no_doi":    lambda: _papers() - _papers_with_doi(),
    "lab_authored":     _lab_authored,
    "papers_raw_file":  lambda: _papers_field("raw_file"),
    "papers_unverified": lambda: _papers_field("verified", value="false"),
    "papers_no_verified": lambda: _papers_field("verified", absent=True),
    "entries":          lambda: _files("knowledge/entries/*.json"),
    "entries_tooling":  lambda: _entries_by_origin()["tooling"],
    "entries_method":   lambda: _entries_by_origin()["method"],
    "entries_handbook": lambda: _entries_by_origin()["handbook"],
    "entries_intake":   lambda: _entries_by_origin()["intake"],
    "entries_paper":    lambda: _entries_by_origin()["paper"],
    "book_claims":      _book_claims,
    "sealed_documents": _sealed_documents,
    "seal_dirs":        lambda: len(list(ROOT.glob("runs*/*/SEALED.sha256"))),
    "index_stated":     _index_stated,
    "oscill_specs":     lambda: _files("specs/chain-bend-2d-oscill__*.json"),
    "oscill_runs":      lambda: len({p.split("/")[1] for p in _tracked("runs/*/*")
                                     if p.count("/") >= 2
                                     and p.split("/")[1].startswith("chain-bend-2d-oscill")}),
}


# ── the claims ─────────────────────────────────────────────────────────────
#
# `kind`:
#   "equal" -> group i must equal COUNTERS[counters[i]]
#   "sum"   -> the groups must SUM to COUNTERS[counters[0]]
#
# Patterns are anchored on surrounding prose rather than on line numbers,
# because line numbers move and a stale line number silently matches nothing.

class Claim:
    def __init__(self, file, pattern, counters, kind="equal", why=""):
        self.file, self.pattern, self.kind, self.why = file, pattern, kind, why
        self.counters = (counters,) if isinstance(counters, str) else tuple(counters)


C = Claim
CLAIMS: list[Claim] = [
    # ── README.md: the badge paragraph and the headline table ──────────────
    C("README.md", r"\*\*(\d+) passed, (\d+) skipped\*\* on HOOMD-blue",
      "tests_collected", kind="sum"),
    C("README.md", r"plus \*\*(\d+) sealed", "sealed_documents"),
    C("README.md", r"\*\*(\d+) specifications · (\d+) run directories\*\* \| (\d+) carrying",
      ("specs", "run_dirs", "metrics")),
    C("README.md", r"\| \*\*(\d+) post-mortems\*\* \|", "postmortems"),
    C("README.md", r"\| \*\*(\d+) verification scripts\*\* \|", "verify_scripts"),
    C("README.md", r"\| \*\*(\d+) tests\*\* \| (\d+) skipped", "tests_collected", kind="sum"),
    C("README.md",
      r"(\d+) wiki pages — (\d+) system cards · (\d+) findings · (\d+) concepts — plus "
      r"(\d+) paper and (\d+) book distillations and (\d+) tool-written entries",
      ("wiki_pages", "system_cards", "findings", "concepts", "papers", "books", "entries")),

    # ── README.md: the ASCII stage and knowledge diagrams ─────────────────
    C("README.md", r"past decisions, (\d+) paper", "distillations"),
    C("README.md", r"rather than a symptom\. (\d+) ", "findings"),
    C("README.md", r"value\. (\d+) ", "benchmarks"),
    C("README.md", r"source/papers/     (\d+) distillations", "papers"),
    C("README.md", r"(\d+) books, (\d+)/(\d+) claims", ("books", "book_claims", "book_claims")),
    C("README.md", r"entries/           (\d+) tool-written entries\. (\d+) of them are",
      ("entries", "entries_tooling")),
    C("README.md", r"record\.json   (\d+) post-mortems", "postmortems"),

    # ── README.md: prose and the closing tables ───────────────────────────
    C("README.md", r"There are (\d+) worked cases in \[`cases/`\]\(cases/\) and (\d+) modules in",
      ("case_scripts", "bdbot_modules")),
    C("README.md",
      r"\| Runs \| \*\*(\d+)\*\* specs · (\d+) run directories · \*\*(\d+)\*\* with "
      r"`metrics\.json` · (\d+) post-mortems \|",
      ("specs", "run_dirs", "metrics", "postmortems")),
    C("README.md",
      r"`bdbot/` (\d+) modules \(L0→L7\) · `simbot/` (\d+) modules \(S2/S6/S7/S8\) · "
      r"(\d+) case scripts · (\d+) verification scripts",
      ("bdbot_modules", "simbot_modules", "case_scripts", "verify_scripts")),
    C("README.md", r"\| Tests \| \*\*(\d+) pass\*\*, (\d+) skipped", "tests_collected", kind="sum"),
    C("README.md",
      r"\| Knowledge \| (\d+) wiki pages · (\d+) paper \+ (\d+) book distillations · (\d+) entries \|",
      ("wiki_pages", "papers", "books", "entries")),
    C("README.md", r"states (\d+) entries where (\d+) files", ("index_stated", "papers"),
      why="both halves of the librarian-agent finding have to stay true"),
    C("README.md", r"L2 engine, L0→L7 \| (\d+) modules", "bdbot_modules"),
    C("README.md", r"S2/S6/S7/S8 \| (\d+) modules", "simbot_modules"),
    C("README.md", r"case physics \| (\d+) scripts", "case_scripts"),
    C("README.md", r"sweep analyses \| (\d+) scripts", "campaign_scripts"),
    C("README.md", r"executable claims \| (\d+) scripts", "verify_scripts"),
    C("README.md", r"the (\d+) contracts", "specs"),
    C("README.md", r"\| (\d+) paper distillations \| in-repo", "papers"),
    C("README.md", r"⚠️ \*\*(\d+)\*\* of them are the group's own", "lab_authored"),
    C("README.md", r"38 `true`, (\d+) `false`, (\d+) with no field, totalling the (\d+)",
      ("papers_unverified", "papers_no_verified", "papers")),
    C("README.md", r"tabulates separately as (\d+) / (\d+)",
      ("papers_raw_file", "index_stated"),
      why="INDEX.md is generated and its own count is stale on purpose -- the "
          "README cites that gap as a finding, so both numbers must stay true"),
    C("README.md", r"\| (\d+) book distillations \(Leal 2026; Welty 5th ed\.\) \| in-repo, "
                   r"(\d+)/(\d+) claims re-derived",
      ("books", "book_claims", "book_claims")),

    # ── NOTICE.md ─────────────────────────────────────────────────────────
    C("NOTICE.md", r"(\d+) distillations in \[", "papers"),
    C("NOTICE.md", r"(\d+) of the (\d+)\ncarry a DOI", ("papers_with_doi", "papers")),
    C("NOTICE.md", r"(\d+)\ncarry no DOI field", "papers_no_doi"),
    C("NOTICE.md", r"(\d+) claims taken out\nof the two books", "book_claims"),
    C("NOTICE.md", r"\((\d+)/(\d+) pass\)", ("book_claims", "book_claims")),
    C("NOTICE.md", r"\*\*(\d+)\*\* of the (\d+) are marked `lab_authored: true`",
      ("lab_authored", "papers")),
    C("NOTICE.md", r"largely unchecked\.\*\* (\d+)\ncarries an explicit `verified: false`",
      "papers_unverified"),
    C("NOTICE.md", r"the other (\d+) carry no `verified` field", "papers_no_verified"),

    # ── docs/01-architecture.md ───────────────────────────────────────────
    C("docs/01-architecture.md", r"Q→A, and dead-ends — (\d+)", "findings"),
    C("docs/01-architecture.md", r"entries/           (\d+) tool-written JSON", "entries"),
    C("docs/01-architecture.md", r"record\.json  (\d+) post-mortems", "postmortems"),

    # ── docs/03-knowledge-base.md ─────────────────────────────────────────
    C("docs/03-knowledge-base.md",
      r"\| Size \| (\d+) wiki pages · (\d+) paper \+ (\d+) book distillations \| (\d+) entries \|",
      ("wiki_pages", "papers", "books", "entries")),
    C("docs/03-knowledge-base.md", r"\| \*\*`findings/`\*\* \| (\d+) \|", "findings"),
    C("docs/03-knowledge-base.md", r"(\d+) entries, by origin", "entries"),
    C("docs/03-knowledge-base.md", r"\| \*\*`tooling`\*\* \| \*\*(\d+)\*\* \|", "entries_tooling"),
    C("docs/03-knowledge-base.md",
      r"\*\*(\d+) tooling entries against (\d+) paper entries",
      ("entries_tooling", "entries_paper")),
    C("docs/03-knowledge-base.md", r"Plus \*\*(\d+) `record\.json`\*\*", "postmortems"),

    # ── docs/06-roadmap.md ────────────────────────────────────────────────
    C("docs/06-roadmap.md", r"✅ (\d+) modules, L0→L7", "bdbot_modules"),
    C("docs/06-roadmap.md",
      r"🔶 (\d+) wiki pages · (\d+) distillations · (\d+) entries · (\d+) post-mortems",
      ("wiki_pages", "distillations", "entries", "postmortems")),
    C("docs/06-roadmap.md",
      r"✅ (\d+) specs · (\d+) run directories · (\d+) with `metrics\.json`",
      ("specs", "run_dirs", "metrics")),

    # ── docs/04-cases.md ──────────────────────────────────────────────────
    C("docs/04-cases.md", r"(\d+) KB `handbook` entries, (\d+)/(\d+) checks",
      ("entries_handbook", "book_claims", "book_claims")),
    C("docs/04-cases.md", r"\*\*L3-only\*\* — (\d+) specs, (\d+) run",
      ("oscill_specs", "oscill_runs")),

    # ── docs/05-pitfalls.md ───────────────────────────────────────────────
    C("docs/05-pitfalls.md", r"There are \*\*(\d+)\*\* `tooling` entries of (\d+)",
      ("entries_tooling", "entries")),
    C("docs/05-pitfalls.md",
      r"(\d+) tooling · (\d+)\nmethod · (\d+) handbook · (\d+) intake · (\d+) paper",
      ("entries_tooling", "entries_method", "entries_handbook",
       "entries_intake", "entries_paper")),
]


# ── the check ──────────────────────────────────────────────────────────────
def run(fix: bool = False, verbose: bool = False) -> int:
    measured: dict[str, int] = {}
    for name, fn in COUNTERS.items():
        measured[name] = fn()

    rows, errors, drift = [], [], []
    edits: dict[str, list[tuple[int, int, str]]] = collections.defaultdict(list)

    for claim in CLAIMS:
        path = ROOT / claim.file
        text = path.read_text()
        hits = list(re.finditer(claim.pattern, text))

        if len(hits) != 1:
            errors.append(f"{claim.file}: pattern matched {len(hits)} times, "
                          f"expected exactly 1 — {claim.pattern!r}")
            continue

        m = hits[0]
        got = [int(g) for g in m.groups()]

        if claim.kind == "sum":
            want_total = measured[claim.counters[0]]
            ok = sum(got) == want_total
            label = f"{' + '.join(map(str, got))} = {sum(got)}"
            rows.append((claim.file, claim.counters[0], label, want_total, ok))
            if not ok:
                drift.append(f"{claim.file}: {' + '.join(map(str, got))} "
                             f"= {sum(got)}, measured {want_total} "
                             f"({claim.counters[0]})")
                if fix:
                    # only the first group is adjustable; the second is the
                    # runtime skip count, which this file cannot derive
                    errors.append(f"{claim.file}: --fix cannot repair a `sum` "
                                  f"claim — run the suite and set both numbers")
            continue

        assert len(got) == len(claim.counters), (claim.file, claim.pattern)
        for i, (g, name) in enumerate(zip(got, claim.counters)):
            want = measured[name]
            ok = g == want
            rows.append((claim.file, name, str(g), want, ok))
            if not ok:
                drift.append(f"{claim.file}: {name} claimed {g}, measured {want}")
                if fix:
                    s, e = m.span(i + 1)
                    edits[claim.file].append((s, e, str(want)))

    if fix and edits:
        for f, spans in edits.items():
            path = ROOT / f
            text = path.read_text()
            for s, e, new in sorted(spans, reverse=True):
                text = text[:s] + new + text[e:]
            path.write_text(text)
            print(f"  rewrote {len(spans)} number(s) in {f}")

    width = max(len(r[1]) for r in rows) if rows else 10
    if verbose:
        print(f"{'file':<28} {'counter':<{width}} {'claimed':>10} {'measured':>9}")
        print("-" * (52 + width))
    for f, name, got, want, ok in rows:
        if verbose or not ok:
            mark = " " if ok else "  <-- DRIFT"
            print(f"{f:<28} {name:<{width}} {got:>10} {want:>9}{mark}")

    print()
    print(f"{len(CLAIMS)} claims, {len(rows)} numbers checked, "
          f"{len(drift)} drifted, {len(errors)} unmatched")
    for e in errors:
        print(f"  ERROR  {e}")
    if fix and drift and not errors:
        print("  re-run without --fix to confirm")
    return 1 if (drift or errors) else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--fix", action="store_true",
                    help="rewrite the claimed numbers to the measured values")
    ap.add_argument("-v", "--verbose", action="store_true",
                    help="print every claim, not only the drifted ones")
    a = ap.parse_args()
    return run(fix=a.fix, verbose=a.verbose)


if __name__ == "__main__":
    raise SystemExit(main())
