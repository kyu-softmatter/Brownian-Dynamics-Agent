#!/usr/bin/env python
"""Compare two markdown files and report every difference that is NOT prose.

The analogue of `verify_translation_safety.py` (token stream, for `.py`) and
`verify_yaml_translation_safety.py` (structure, for `.yaml`), for the third and
largest surface: `docs/` and `knowledge/` are ~10,000 lines of markdown.

What a translation is allowed to change: **the words**. What it must not change:

  headings      the sequence of ATX heading levels (#, ##, ...) -- not their text
  code fences   the fence infostring and every line of every fenced block, byte
                for byte. A fenced line that carried Hangul in OLD is a code line
                with a translated comment, so only the **prefix up to the first
                Hangul character** is required to survive -- `x = 1  # KO` may
                become `x = 1  # a comment`, but not `x = 2  # a comment`.
  links         every link/image TARGET, in order (the label is prose, the URL is not)
  tables        for each table, the row count and the per-row cell count
  numbers       the multiset of numeric literals -- this is the one that matters.
                CLAUDE.md's whole premise is that a number in a document is the
                return value of a function; a translation that perturbs one is not
                a translation. Percent signs, exponents and separators are kept
                attached so that 1e-2, 1.65% and 1,255 stay distinguishable.
  markers       the multiset of correction/emphasis markers the user asked to be
                preserved explicitly -- ★ ⭐ ⚠️ ⛔ ✅ ❌ ⟹ → and the circled
                digits, plus the count of **bold** runs.

**The policy this checker assumes: preserve every non-alphabetic glyph exactly.**
Calibrating it against the four markdown files this repository already translated
(`.claude/skills/bd-{hoomd,physics,diagnose}/SKILL.md`,
`examples/soft-r3-time-resolved/08_conclusion.md`) produced 77 findings, and the
ones that were not real content changes were all one thing: those translations
ASCII-ified typography -- `# <- ` for `# ← `, `Stokes-Einstein` for
`Stokes–Einstein`, `->` for `→`. Weakening the comparison to tolerate that would
also hide a real edit inside a code fence, so the checker stays strict and the
translation keeps the glyph. That is also what the user asked for on
`docs/history/`: the correction markers are preserved explicitly.

Exit 0 when the only differences are prose; 1 otherwise. Every finding prints the
path and the first differing item, never a bare "FAIL".

    $PY verify/verify_markdown_translation_safety.py OLD NEW
    $PY verify/verify_markdown_translation_safety.py --selftest
"""
import re
import sys
from collections import Counter
from pathlib import Path

KO = re.compile(r"[ᄀ-ᇿ㄰-㆏가-힯]")
FENCE = re.compile(r"^(\s*)(```+|~~~+)(.*)$")
HEADING = re.compile(r"^(#{1,6})\s")
# a link or image target: [label](target) / ![label](target). Nested parens in the
# target are not markdown-legal, so a non-greedy stop at ) is correct here.
LINK = re.compile(r"!?\[[^\]]*\]\(([^)]*)\)")
# numeric literals, sign and exponent and trailing % kept attached
NUM = re.compile(r"(?<![\w.])[+\-−]?\d[\d,]*(?:\.\d+)?(?:[eE][+\-]?\d+)?%?")
MARKERS = "★⭐⚠⛔✅❌⟹→" + "".join(
    chr(c) for c in range(0x2460, 0x2474))
BOLD = re.compile(r"\*\*[^*\n]+\*\*")


def parse(text: str) -> dict:
    """Everything about a markdown file that a translation must not change."""
    headings, fences, tables, fence_ko = [], [], [], 0
    in_fence, fence_close, cur_fence = False, "", []
    cur_table = []
    for line in text.split("\n"):
        m = FENCE.match(line)
        if in_fence:
            if m and m.group(2).startswith(fence_close):
                fences.append(cur_fence)
                in_fence, cur_fence = False, []
            else:
                # a fenced line that carried Hangul is a code line whose comment
                # (or string) is prose. Keep the code that precedes the Hangul as
                # a required prefix; the rest is free to be translated.
                m2 = KO.search(line)
                if m2:
                    fence_ko += 1
                    cur_fence.append(("KO-PREFIX", line[:m2.start()]))
                else:
                    cur_fence.append(line)
            continue
        if m:
            in_fence, fence_close = True, m.group(2)
            cur_fence = [("FENCE-OPEN", m.group(1), m.group(2), m.group(3).strip())]
            continue
        h = HEADING.match(line)
        if h:
            headings.append(len(h.group(1)))
        # a table row: starts and ends with | after stripping
        s = line.strip()
        if s.startswith("|") and s.endswith("|") and len(s) > 1:
            cur_table.append(s.count("|"))
        elif cur_table:
            tables.append(cur_table)
            cur_table = []
    if in_fence:                      # unterminated fence: keep what we have
        fences.append(cur_fence)
    if cur_table:
        tables.append(cur_table)

    # links/numbers/markers are read from the WHOLE file, fences included:
    # a URL or a number inside a code block is exactly as load-bearing.
    return {
        "headings": headings,
        "fences": fences,
        "tables": tables,
        "links": LINK.findall(text),
        "numbers": Counter(NUM.findall(text)),
        "markers": Counter(ch for ch in text if ch in MARKERS),
        "bold": len(BOLD.findall(text)),
        "_fence_ko": fence_ko,
    }


def _seq_diff(name: str, a: list, b: list) -> list:
    if a == b:
        return []
    out = [f"{name}: {len(a)} -> {len(b)}" if len(a) != len(b) else f"{name}: same length, differing content"]
    for i, (x, y) in enumerate(zip(a, b)):
        if x != y:
            out.append(f"    first difference at index {i}:")
            out.append(f"      OLD {x!r}")
            out.append(f"      NEW {y!r}")
            break
    return out


def _fence_diff(name: str, a: list, b: list) -> list:
    """Like _seq_diff, but a ("KO-PREFIX", prefix) entry in OLD only requires that
    the NEW line at that position still starts with `prefix`."""
    if len(a) != len(b):
        return [f"{name}: {len(a)} -> {len(b)} line(s)"]
    for i, (x, y) in enumerate(zip(a, b)):
        if isinstance(x, tuple) and x and x[0] == "KO-PREFIX":
            if isinstance(y, str) and y.startswith(x[1]):
                continue
            return [f"{name}: translated code comment at line {i} lost its code",
                    f"      OLD must start with {x[1]!r}",
                    f"      NEW {y!r}"]
        if x != y:
            return [f"{name}: line {i} differs",
                    f"      OLD {x!r}",
                    f"      NEW {y!r}"]
    return []


def compare(old: str, new: str) -> list:
    A, B = parse(old), parse(new)
    p = []
    p += _seq_diff("heading levels", A["headings"], B["headings"])
    if len(A["fences"]) != len(B["fences"]):
        p.append(f"fenced code blocks: {len(A['fences'])} -> {len(B['fences'])}")
    else:
        for i, (fa, fb) in enumerate(zip(A["fences"], B["fences"])):
            p += _fence_diff(f"code fence #{i}", fa, fb)
    p += _seq_diff("table shapes", A["tables"], B["tables"])
    p += _seq_diff("link targets", A["links"], B["links"])
    for label, key in (("numbers", "numbers"), ("markers", "markers")):
        if A[key] != B[key]:
            lost = A[key] - B[key]
            gained = B[key] - A[key]
            p.append(f"{label} changed: lost {dict(lost)} gained {dict(gained)}")
    if A["bold"] != B["bold"]:
        p.append(f"bold runs: {A['bold']} -> {B['bold']}")
    return p


SELFTESTS = [
    # (label, old, new, expect_problem)
    ("identical", "# A\n\ntext 1.5\n", "# A\n\ntext 1.5\n", False),
    ("prose only", "# A\n\n한국어 1.5\n", "# A\n\nEnglish 1.5\n", False),
    ("heading level changed", "# A\n\n## B\n", "# A\n\n### B\n", True),
    ("heading dropped", "# A\n\n## B\n", "# A\n", True),
    ("number changed", "value 1.65%\n", "value 1.66%\n", True),
    ("number dropped", "a 42 b 7\n", "a 42 b\n", True),
    ("exponent changed", "dt = 1e-2\n", "dt = 1e-3\n", True),
    ("thousands sep kept", "1,255 lines\n", "1,255 rows\n", False),
    ("link target changed", "[x](a.md)\n", "[y](b.md)\n", True),
    ("link label only", "[한국어](a.md)\n", "[English](a.md)\n", False),
    ("code line changed", "```py\nx = 1\n```\n", "```py\nx = 2\n```\n", True),
    ("code comment translated", "```py\nx = 1  # 주석\n```\n",
     "```py\nx = 1  # a comment\n```\n", False),
    ("code changed under a translated comment", "```py\nx = 1  # 주석\n```\n",
     "```py\nx = 2  # a comment\n```\n", True),
    ("string translated in a fence", "```py\nprint(\"한국어\")\n```\n",
     "```py\nprint(\"English\")\n```\n", False),
    ("call renamed under a translated string", "```py\nprint(\"한국어\")\n```\n",
     "```py\nputs(\"English\")\n```\n", True),
    ("fence infostring changed", "```py\nx = 1\n```\n", "```sh\nx = 1\n```\n", True),
    ("fence dropped", "```py\nx = 1\n```\n", "x = 1\n", True),
    ("table column dropped", "| a | b |\n|---|---|\n| 1 | 2 |\n",
     "| a |\n|---|\n| 1 |\n", True),
    ("table cell prose", "| 한국어 | b |\n|---|---|\n", "| English | b |\n|---|---|\n", False),
    ("marker dropped", "★ important\n", "important\n", True),
    ("warning marker dropped", "⚠ careful\n", "careful\n", True),
    ("bold count changed", "**a** and **b**\n", "**a** and b\n", True),
]


def selftest() -> int:
    bad = 0
    for label, old, new, expect in SELFTESTS:
        got = compare(old, new)
        ok = bool(got) == expect
        print(f"  {'PASS' if ok else 'FAIL'}  {label:28} "
              f"expected={'problem' if expect else 'clean'} got={len(got)} finding(s)")
        if not ok:
            bad += 1
            for line in got:
                print(f"          {line}")
    print(f"selftest: {len(SELFTESTS) - bad}/{len(SELFTESTS)} passed")
    return 1 if bad else 0


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    if len(sys.argv) != 3:
        print(__doc__)
        return 2
    old_p, new_p = Path(sys.argv[1]), Path(sys.argv[2])
    old, new = old_p.read_text(encoding="utf-8"), new_p.read_text(encoding="utf-8")
    problems = compare(old, new)
    ko_left = sum(1 for ln in new.split("\n") if KO.search(ln))
    fence_ko = parse(old)["_fence_ko"]
    if problems:
        print(f"DIFF {new_p}: {len(problems)} structural difference(s)")
        for line in problems:
            print("  " + line)
        return 1
    print(f"OK   {new_p}: structure identical, only prose differs "
          f"(fenced lines with Hangul in OLD: {fence_ko}), Hangul lines left {ko_left}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
