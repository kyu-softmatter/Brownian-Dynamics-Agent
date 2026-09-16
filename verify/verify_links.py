"""Every repo-relative link in a tracked Markdown file, resolved.

    $PY verify/verify_links.py          # non-zero exit if anything is broken
    $PY verify/verify_links.py -v       # print the allowed categories too

**Why this file exists.** 39 links were broken across 17 files, measured
2026-09-15, and two of them were load-bearing rather than cosmetic:

  · `knowledge/wiki/benchmarks/choi2020-interfacial-rdf.md` and its system card
    both cited `docs/tools/digitize_fig4b.py` as the script that produced the
    digitized `g(r)` tolerances and said the reading was *"재현 가능"*
    (reproducible with it). That file — and the whole `docs/tools/` directory —
    has never been in this repository. A verification tolerance rested on an
    artefact that does not exist.
  · seven links still pointed at `inputs/`, which was renamed to `intake/` at the
    2026-08-28 merge. One of them was live report-generating code, so every
    report it wrote carried the dead path forward.

A link is the cheapest possible claim to check and the easiest to let rot, which
is exactly the shape this repository keeps recording.

**Four categories, and only one of them is a failure.** Lumping them would make
the gate either useless or unpassable:

| category | meaning | gates? |
|---|---|---|
| `ok` | the target exists | — |
| `unpublished` | the target matches a `.gitignore` rule | no — `NOTICE.md` says trajectories and per-run figures are not published |
| `frozen` | the **containing file is under seal** | no — editing it would break the seal, and `.claude/settings.json` denies the edit |
| `historical` | the containing file is under `docs/history/` | no — copies of the predecessor repositories' own documents, whose references point into repositories that no longer exist |
| `broken` | anything else | **yes** |

`frozen` is derived from the seals themselves, not from a list. Two archived
`01_intake.md` documents link to `inputs/…` and are sealed, so the rename cannot
be applied to them; that is a consequence of sealing working, not a defect.

**The guard that makes the number mean anything.** A scan that finds no links
reports no breakage. The link count is asserted to be in the hundreds, so an
empty or mis-globbed scan fails instead of passing.
"""
from __future__ import annotations

import argparse
import collections
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
LINK = re.compile(r"\[[^\]]*\]\(([^)\s]+)\)")
EXTERNAL = ("http://", "https://", "mailto:", "#")
HISTORICAL_PREFIXES = ("docs/history/",)

#: A scan that matches nothing cannot fail. 502 repo-relative links were found on
#: 2026-09-15; the floor is well below that so ordinary editing does not trip it,
#: and well above zero so a broken glob does.
MIN_LINKS = 400


def _tracked_markdown() -> list[pathlib.Path]:
    out = subprocess.run(["git", "ls-files", "-z", "*.md"], cwd=ROOT,
                         capture_output=True, check=True).stdout
    return [ROOT / f.decode("utf-8", "surrogateescape")
            for f in out.split(b"\0") if f]


def sealed_documents() -> set[str]:
    """Repo-relative paths of every document covered by a tracked seal, plus the
    seals. Derived from the seal files, because a hand-written list is what
    drifts."""
    out: set[str] = set()
    for seal in sorted(ROOT.glob("runs*/*/SEALED.sha256")):
        out.add(seal.relative_to(ROOT).as_posix())
        for line in seal.read_text().splitlines():
            if not line.strip():
                continue
            _, _, rel = line.partition("  ")
            sibling = seal.parent / pathlib.PurePath(rel.strip()).name
            out.add(sibling.relative_to(ROOT).as_posix())
    return out


def _ignored(paths: list[pathlib.Path]) -> set[str]:
    """One `git check-ignore` call for the whole batch; per-link calls made this
    the slowest part of the check by an order of magnitude."""
    if not paths:
        return set()
    payload = "\0".join(p.as_posix() for p in paths)
    r = subprocess.run(["git", "check-ignore", "--stdin", "-z"], cwd=ROOT,
                       input=payload, capture_output=True, text=True)
    return {line for line in r.stdout.split("\0") if line}


def scan(with_rows: bool = False):
    """`(tally, broken)`, or `(tally, broken, classified)` with `with_rows`.

    `classified` is `[(category, file, target, resolved_relpath)]` for every link
    that does not resolve. It exists because a size cap on the allowed categories
    is not enough: a mutation that classified EVERY unresolvable link as
    `unpublished` kept all three caps satisfied and passed. The classification
    itself has to be checkable, so `tests/test_links.py` re-derives it.
    """
    return _scan(with_rows)


def _scan(with_rows: bool = False):
    sealed = sealed_documents()
    rows: list[tuple[str, str, pathlib.Path]] = []
    for path in _tracked_markdown():
        rel_file = path.relative_to(ROOT).as_posix()
        for m in LINK.finditer(path.read_text(errors="replace")):
            target = m.group(1)
            if target.startswith(EXTERNAL):
                continue
            bare = target.split("#")[0]
            if not bare:
                continue
            rows.append((rel_file, target, (path.parent / bare)))

    missing = [p for _, _, p in rows if not p.exists()]
    ignored = _ignored(missing)

    tally: collections.Counter = collections.Counter()
    broken: list[tuple[str, str, str]] = []
    classified: list[tuple[str, str, str, str]] = []
    for rel_file, target, resolved in rows:
        if resolved.exists():
            tally["ok"] += 1
            continue
        try:
            rel_target = resolved.resolve().relative_to(ROOT).as_posix()
        except ValueError:
            rel_target = resolved.as_posix()
        if rel_target in ignored or resolved.as_posix() in ignored:
            kind = "unpublished"
        elif rel_file in sealed:
            kind = "frozen"
        elif rel_file.startswith(HISTORICAL_PREFIXES):
            kind = "historical"
        else:
            kind = "broken"
            broken.append((rel_file, target, rel_target))
        tally[kind] += 1
        classified.append((kind, rel_file, target, rel_target))
    if with_rows:
        return tally, broken, classified
    return tally, broken


def run(verbose: bool = False) -> int:
    tally, broken = scan()
    total = sum(tally.values())

    if total < MIN_LINKS:
        print(f"ERROR: only {total} repo-relative links found (floor {MIN_LINKS}) "
              f"— the scan is not seeing the documentation")
        return 1

    for rel_file, target, resolved in broken:
        print(f"BROKEN  {rel_file}\n          -> {target}")
    if verbose:
        for k in ("unpublished", "frozen", "historical"):
            print(f"  allowed: {tally[k]:>3}  {k}")

    print(f"\n{total} repo-relative links · {tally['ok']} resolve · "
          f"{tally['unpublished']} unpublished · {tally['frozen']} frozen by a seal · "
          f"{tally['historical']} historical · {tally['broken']} broken")
    return 1 if tally["broken"] else 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("-v", "--verbose", action="store_true")
    return run(verbose=ap.parse_args().verbose)


if __name__ == "__main__":
    raise SystemExit(main())
