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
  links         every link/image TARGET, in order (the label is prose, the URL is
                not). One exception: an IN-DOCUMENT anchor (`](#...)`) is derived
                from the heading text, which IS prose, so translating a heading
                must retarget it. Those are not compared as strings; instead each
                is required to RESOLVE to a heading in the same file -- a check
                the string comparison could not make. An anchor already broken in
                OLD is reported as such rather than counted against NEW.
  tables        for each table, the row count and the per-row cell count
  numbers       the multiset of numeric literals -- this is the one that matters.
                CLAUDE.md's whole premise is that a number in a document is the
                return value of a function; a translation that perturbs one is not
                a translation. Percent signs, exponents and separators are kept
                attached so that 1e-2, 1.65% and 1,255 stay distinguishable.
                One exception, and it is measured rather than assumed: a numeral
                IMMEDIATELY ADJACENT to a Hangul character on either side is part
                of a Korean word, not a measurement -- the 2 of `2차극소` is the
                "second" of "secondary minimum", the 3 of `3종` is "three", the 1
                of `1개` is "one", the 1 of `게이트1` is the gate's number. Those become English words, so they are censused
                separately and EXEMPTED in both directions: such a numeral may
                vanish (`2차극소` -> "secondary minimum") or persist (`28종` ->
                "28 isotropic"), and neither is a finding. What is still a
                finding: a numeral in OLD that was not Hangul-attached and is now
                missing, or a numeral in NEW that OLD's exempt pool cannot
                account for. Calibrating on docs/hoomd_capabilities.md, 9 of the
                11 number findings were exactly this and 2 were real bugs I had
                introduced.
                One case the census cannot settle by itself: Korean writes
                magnitudes as `330만` (3.3 million) and `2,085만` (20.85 million),
                where the VALUE survives translation but the digits cannot. There
                is no digit-preserving English form, so such a rescale has to be
                accepted BY NAME with `--allow-gained 3.3,20.85`. It is not
                absorbed automatically: an automatic budget would silently swallow
                whichever unexplained numeral sorted first, which on
                docs/history/2026-07-30_simulation_bot_master_plan.md was two
                real defects rather than the rescale. The count of myriad forms
                available in OLD is reported so the claim can be checked.
  markers       the multiset of correction/emphasis markers the user asked to be
                preserved explicitly -- ★ ⭐ ⚠️ ⛔ ✅ ❌ ⟹ → and the circled
                digits, plus the number of `**` emphasis delimiters. The
                delimiters are counted rather than the runs because a bold span
                may cross a line break, and a translation re-wraps lines -- on
                docs/history/2026-08_simulation_auto_CLAUDE.md that alone
                moved a run count from 356 to 362 with no emphasis added or lost.

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
# Korean magnitude words. `330만` is 3.3 million and `2,085만` is 20.85 million:
# the DIGITS cannot survive translation, only the value. Each occurrence buys the
# translation one otherwise-unexplained numeral (see _number_census).
KO_SCALE = "만억조천백십"
FENCE = re.compile(r"^(\s*)(```+|~~~+)(.*)$")
HEADING = re.compile(r"^(#{1,6})\s")
# a link or image target: [label](target) / ![label](target). Nested parens in the
# target are not markdown-legal, so a non-greedy stop at ) is correct here.
LINK = re.compile(r"!?\[[^\]]*\]\(([^)]*)\)")
HEAD_TEXT = re.compile(r"^#{1,6}\s+(.*?)\s*$")
# GitHub's heading -> anchor rule, near enough: lowercase, drop everything that is
# not a word character / space / hyphen, then spaces -> hyphens.
_SLUG_DROP = re.compile(r"[^\w\s-]", re.UNICODE)


def slug(heading_text: str) -> str:
    t = _SLUG_DROP.sub("", heading_text.strip().lower())
    return re.sub(r"\s+", "-", t)
# numeric literals, sign and exponent and trailing % kept attached
NUM = re.compile(r"(?<![\w.])[+\-−]?\d[\d,]*(?:\.\d+)?(?:[eE][+\-]?\d+)?%?")
# the same, without the lookbehind, so that a numeral GLUED to a Hangul word
# (`게이트1`) is findable at all -- NUM cannot see it, because Hangul is \w
NUM_ANY = re.compile(r"[+\-−]?\d[\d,]*(?:\.\d+)?(?:[eE][+\-]?\d+)?%?")
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
    heading_slugs = {slug(m.group(1)) for line in text.split("\n")
                     if (m := HEAD_TEXT.match(line))}
    all_links = LINK.findall(text)
    return {
        "heading_slugs": heading_slugs,
        "anchors": [t for t in all_links if t.startswith("#")],
        "headings": headings,
        "fences": fences,
        "tables": tables,
        "links": [t for t in all_links if not t.startswith("#")],
        **_number_census(text),
        "markers": Counter(ch for ch in text if ch in MARKERS),
        "bold": text.count("**"),
        "_fence_ko": fence_ko,
    }


def _number_census(text: str) -> dict:
    """Split the numerals into the ones a translation must preserve and the ones
    it is free to turn into words.

    A numeral adjacent to Hangul on EITHER side belongs to a Korean word, so it
    is exempt. A numeral in the middle of a longer token (the `0` of `7.1.0`) is
    neither -- it is not a separate number and is dropped from both.
    """
    census, exempt, scaled = Counter(), Counter(), 0
    for m in NUM_ANY.finditer(text):
        before = text[m.start() - 1] if m.start() else ""
        after = text[m.end():m.end() + 1]
        if (before and KO.match(before)) or KO.match(text, m.end()):
            exempt[m.group(0)] += 1
            if after in KO_SCALE:
                scaled += 1          # a myriad form: the value survives, the digits do not
        elif NUM.match(text, m.start()):
            census[m.group(0)] += 1
        # else: a digit run inside a longer token -- not a number of its own
    return {"numbers": census, "_ko_numerals": exempt, "_ko_scaled": scaled}


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


def compare(old: str, new: str, allow_gained=(), allow_lost=()) -> list:
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
    if len(A["anchors"]) != len(B["anchors"]):
        p.append(f"in-document anchors: {len(A['anchors'])} -> {len(B['anchors'])}")
    bad_old = [a for a in A["anchors"] if a[1:] not in A["heading_slugs"]]
    bad_new = [a for a in B["anchors"] if a[1:] not in B["heading_slugs"]]
    # only anchors that NEW broke count; ones already broken in OLD are noted
    if len(bad_new) > len(bad_old):
        p.append(f"in-document anchors that resolve to no heading: "
                 f"{bad_new} (already broken in OLD: {bad_old})")
    # numbers: OLD's Hangul-attached numerals are exempt in BOTH directions
    new_all = B["numbers"] + B["_ko_numerals"]
    lost = A["numbers"] - new_all
    gained = new_all - A["numbers"] - A["_ko_numerals"]
    for tok in allow_gained:
        if gained.get(tok):
            gained[tok] -= 1
            if not gained[tok]:
                del gained[tok]
    for tok in allow_lost:
        if lost.get(tok):
            lost[tok] -= 1
            if not lost[tok]:
                del lost[tok]
    if lost or gained:
        p.append(f"numbers changed: lost {dict(lost)} gained {dict(gained)}")
    if A["markers"] != B["markers"]:
        p.append(f"markers changed: lost {dict(A['markers'] - B['markers'])} "
                 f"gained {dict(B['markers'] - A['markers'])}")
    if A["bold"] != B["bold"]:
        p.append(f"** emphasis delimiters: {A['bold']} -> {B['bold']}")
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
    ("korean ordinal becomes a word", "2차극소 유지\n", "the secondary minimum holds\n", False),
    # both numerals here are Hangul-attached, so 21 persisting and 3 becoming
    # "three" are both allowed -- this is the exemption working, not a miss
    ("korean counters, one kept one worded", "원 21개 · 3종 전부\n",
     "21 circles, all three\n", False),
    ("latin-suffixed numeral still checked", "10mM salt\n", "11mM salt\n", True),
    ("numeral, space, hangul still checked", "N=1528 에서\n", "at N=1527\n", True),
    ("korean-attached numeral may persist", "등방 페어 28종\n", "28 isotropic pairs\n", False),
    ("korean-attached numeral may vanish", "등방 페어 28종\n", "isotropic pairs\n", False),
    ("but an unaccounted new numeral is a finding", "등방 페어 28종\n",
     "28 isotropic pairs over 3 runs\n", True),
    ("and it is only exempt up to its count", "2차극소 하나\n", "2 and 2 minima\n", True),
    ("hangul-PREFIXED numeral may appear", "게이트1 통과\n", "gate 1 passed\n", False),
    ("mid-token digit is not a number", "hoomd 7.1.0\n", "hoomd 7.1.0 build\n", False),
    ("a myriad rescale is NOT absorbed automatically", "330만 스텝\n",
     "3.3 million steps\n", True),
    ("but a version bump is", "hoomd 7.1.0\n", "hoomd 7.2.0\n", True),
    ("link target changed", "[x](a.md)\n", "[y](b.md)\n", True),
    ("link label only", "[한국어](a.md)\n", "[English](a.md)\n", False),
    ("anchor retargeted with its heading", "## 한국어 제목\n\n[x](#한국어-제목)\n",
     "## English Heading\n\n[x](#english-heading)\n", False),
    ("anchor left pointing at the old slug", "## 한국어 제목\n\n[x](#한국어-제목)\n",
     "## English Heading\n\n[x](#한국어-제목)\n", True),
    ("an anchor already broken in OLD is not charged to NEW",
     "## 한국어 제목\n\n[x](#nope)\n", "## English Heading\n\n[x](#nope)\n", False),
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
    ("bold across a line break is still counted", "**한국어\n둘째줄**\n",
     "**English on one line**\n", False),
    ("but losing that bold is a finding", "**한국어\n둘째줄**\n",
     "English on one line\n", True),
]


ALLOW_TESTS = [
    ("--allow-gained accepts by name", "330만 스텝\n", "3.3 million steps\n",
     ("3.3",), (), False),
    ("--allow-gained accepts only the named token", "330만 스텝\n",
     "3.3 million steps over 12 runs\n", ("3.3",), (), True),
    ("--allow-lost accepts by name", "42 runs\n", "many runs\n", (), ("42",), False),
]


def selftest() -> int:
    bad = 0
    for label, old, new, ag, al, expect in ALLOW_TESTS:
        got = compare(old, new, ag, al)
        ok = bool(got) == expect
        print(f"  {'PASS' if ok else 'FAIL'}  {label:28} "
              f"expected={'problem' if expect else 'clean'} got={len(got)} finding(s)")
        if not ok:
            bad += 1
            for line in got:
                print(f"          {line}")
    for label, old, new, expect in SELFTESTS:
        got = compare(old, new)
        ok = bool(got) == expect
        print(f"  {'PASS' if ok else 'FAIL'}  {label:28} "
              f"expected={'problem' if expect else 'clean'} got={len(got)} finding(s)")
        if not ok:
            bad += 1
            for line in got:
                print(f"          {line}")
    total = len(SELFTESTS) + len(ALLOW_TESTS)
    print(f"selftest: {total - bad}/{total} passed")
    return 1 if bad else 0


def main() -> int:
    if "--selftest" in sys.argv:
        return selftest()
    def _csv(flag):
        for i, a in enumerate(sys.argv):
            if a == flag and i + 1 < len(sys.argv):
                return tuple(t for t in sys.argv[i + 1].split(",") if t)
        return ()
    allow_gained, allow_lost = _csv("--allow-gained"), _csv("--allow-lost")
    skip, positional = set(), []
    for i, a in enumerate(sys.argv[1:], 1):
        if a in ("--allow-gained", "--allow-lost"):
            skip.add(i + 1)
            continue
        if i in skip or a.startswith("--"):
            continue
        positional.append(a)
    if len(positional) != 2:
        print(__doc__)
        return 2
    old_p, new_p = Path(positional[0]), Path(positional[1])
    old, new = old_p.read_text(encoding="utf-8"), new_p.read_text(encoding="utf-8")
    problems = compare(old, new, allow_gained, allow_lost)
    ko_left = sum(1 for ln in new.split("\n") if KO.search(ln))
    A = parse(old)
    fence_ko = A["_fence_ko"]
    ko_num = sum(A["_ko_numerals"].values())
    if problems:
        print(f"DIFF {new_p}: {len(problems)} structural difference(s)")
        for line in problems:
            print("  " + line)
        return 1
    print(f"OK   {new_p}: structure identical, only prose differs "
          f"(fenced lines with Hangul in OLD: {fence_ko}, "
          f"Hangul-attached numerals in OLD: {ko_num}, of them myriad forms: "
          f"{A['_ko_scaled']}), Hangul lines left {ko_left}"
          + (f", numerals accepted by name: gained={list(allow_gained)} "
             f"lost={list(allow_lost)}" if (allow_gained or allow_lost) else ""))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
