#!/usr/bin/env python
"""Prove that translating a Python file changed ONLY comments and string contents.

The Korean -> English migration edits hundreds of files mechanically. `ast.parse`
is not enough: it accepts a file where a number, a variable name or an f-string
field moved. This does the stronger check.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python

    $PY verify/verify_translation_safety.py --selftest          # prove it has power
    $PY verify/verify_translation_safety.py before.py after.py  # compare two files
    $PY verify/verify_translation_safety.py --apply t.py subs   # translate + check

**How it works.** Both files are tokenised and every token is compared, with the
*contents* of string tokens blanked to a single opaque marker. If the resulting
token streams are identical, no executable code was altered -- only comments and
string contents, which is exactly what a translation is allowed to touch.

Runs of *adjacent* string tokens collapse to one marker, because Python
concatenates adjacent literals implicitly, so re-wrapping one long string as two
shorter ones is the same expression -- and re-wrapping is unavoidable when an
English line grows past the column limit. Strings separated by an operator (two
arguments, say) keep that operator between them and are still compared
positionally, so this does not weaken the check.

**Two rules for writing the replacement text**, both learned by having the checker
reject the same mistake repeatedly across ~45 files:

  1. **Never remove a `;` statement separator.** Splitting
     `ax.set_xlabel(...); ax.set_ylabel(...)` onto two lines is a code change.
     Rejected four times before it became a habit; keep the original line shape.
  2. **Never introduce a line break between two adjacent `{...}` fields** in an
     f-string. Adjacent string literals concatenate implicitly, so re-wrapping is
     normally fine -- but splitting *between two replacement fields* inserts an
     f-string boundary where there was none. Keep the original split points.

Also: reordering two f-string fields is rejected, and correctly so -- it silently
swaps which value prints where.

**Why the self-test matters.** The first version of this checker "passed" its own
adversarial test, and the reason was that the `sed` used to inject the bug never
matched -- so PASS meant "nothing was tested". That is the same defect family this
project keeps finding in its own checkers:

  - `chk_doc` computed the right viscosity and asserted nothing about the value
    printed in the document it was checking (verify_book_claims.py, 2026-08-29)
  - `chk(want=0.0, rtol=1.0)` reads as a 100% tolerance and evaluates as
    "exactly zero", because the tolerance is rtol*max(|want|, 1e-300)
  - an audit scoped by bare `git diff` audits nothing once another session stages
    a file, and still reports PASS (verify_merge_equivalence.py)
  - a syntax tally that counted failures by substring-matching its own message,
    which would silently read 0 if the message were translated and the filter
    were not (verify_skill_snippets.py)
  - and the one that had already gone wrong in the tree: verify_bdbot.py asserted
    `v == '<the Korean verdict string>'` against a verdict that bdbot/checks.py
    had rendered in
    English since c0074a2, so it FAILed for a day unnoticed -- it is not in the
    pytest suite, so the suite stayed green. Fixed by asserting the proposition
    (no hard failure, exactly one soft failure) instead of the wording.

In every case the check could not distinguish "I looked and it was fine" from "I
never looked". So `--selftest` asserts each mutation *exists* before asking the
checker to catch it, and fails loudly if a mutation could not be constructed.
"""
from __future__ import annotations

import argparse
import io
import re
import sys
import token
import tokenize
from pathlib import Path

# ⛔ The one deliberately Korean line in this file, and it cannot be otherwise:
#    this IS the character class that detects Korean. Same load-bearing-Korean
#    rule as bdbot/interactions.py -- do not "finish the translation" here.
KO = re.compile(r"[ㄱ-ㆎ가-힣]")

# Structural tokens that carry no meaning for this comparison. COMMENT is skipped
# because comments are exactly what we permit to change.
_SKIP = {token.COMMENT, token.NL, token.NEWLINE, token.INDENT, token.DEDENT,
         token.ENCODING, token.ENDMARKER}

_STRINGY = {token.STRING}
for _name in ("FSTRING_START", "FSTRING_MIDDLE", "FSTRING_END"):
    if hasattr(token, _name):
        _STRINGY.add(getattr(token, _name))


def code_tokens(src: str) -> list[tuple[str, str | None]]:
    """The executable token stream, with string contents blanked."""
    out: list[tuple[str, str | None]] = []
    for tok in tokenize.generate_tokens(io.StringIO(src).readline):
        if tok.type in _SKIP:
            continue
        if tok.type in _STRINGY:
            if not (out and out[-1] == ("STR", None)):
                out.append(("STR", None))
        else:
            out.append((token.tok_name[tok.type], tok.string))
    return out


def compare(before: str, after: str) -> list[str]:
    """Return a list of problems; empty means only comments/strings changed."""
    ta, tb = code_tokens(before), code_tokens(after)
    if ta == tb:
        return []
    n = min(len(ta), len(tb))
    i = next((k for k in range(n) if ta[k] != tb[k]), n)
    return [f"code tokens differ at index {i} ({len(ta)} vs {len(tb)} tokens)\n"
            f"    before: {ta[max(0, i - 3):i + 3]}\n"
            f"    after : {tb[max(0, i - 3):i + 3]}"]


def check_files(a: Path, b: Path) -> int:
    before, after = a.read_text(encoding="utf-8"), b.read_text(encoding="utf-8")
    problems = compare(before, after)
    left = sum(1 for line in after.split("\n") if KO.search(line))
    if problems:
        print(f"FAIL {b}")
        for p in problems:
            print("  " + p)
        return 1
    print(f"OK   {b}: {len(code_tokens(after))} code tokens identical, "
          f"Hangul lines left {left}")
    return 0


# ════════════════════════════════════════════════════════════════════════
# Self-test -- deliberately break it and confirm it notices
# ════════════════════════════════════════════════════════════════════════
_SAMPLE = '''\
"""Docstring with a comment-like sentence."""
import math

SCALE = 3.25                      # a numeric literal
NAME = "unchanged"


def area(radius, height=1.0):
    """Compute something."""
    base = math.pi * radius ** 2   # a comment
    total = base * height
    print(f"radius {radius:.2f} -> area {total:.4f}")
    return total
'''


def _mutations() -> dict[str, tuple[str, bool]]:
    """name -> (mutated source, must_be_rejected)."""
    s = _SAMPLE
    out: dict[str, tuple[str, bool]] = {}

    m = re.search(r"^SCALE = (3\.25)", s, re.M)
    out["number"] = (s[:m.start(1)] + "9.75" + s[m.end(1):], True)

    m = re.search(r"^(\s+)(base)( = math)", s, re.M)
    out["varname"] = (s[:m.start(2)] + "basis" + s[m.end(2):], True)

    m = re.search(r"base \* height", s)
    out["operator"] = (s[:m.start()] + "base / height" + s[m.end():], True)

    lines = s.split("\n")
    i = next(k for k, l in enumerate(lines) if l.strip().startswith("total = "))
    out["dropped_line"] = ("\n".join(lines[:i] + lines[i + 1:]), True)

    m = re.search(r'\{total:\.4f\}', s)
    out["fstring_field"] = (s[:m.start()] + "{base:.4f}" + s[m.end():], True)

    m = re.search(r'"unchanged"', s)
    out["extra_argument"] = (s[:m.start()] + '"unchanged", "extra"' + s[m.end():], True)

    # Legitimate translation edits -- these MUST be accepted.
    out["comment_translated"] = (
        s.replace("# a numeric literal", "# EIN Zahlenliteral")
         .replace("# a comment", "# EIN Kommentar"), False)
    out["string_rewrapped"] = (
        s.replace('"""Compute something."""',
                  '"""Compute something, at greater length so it must wrap."""'), False)
    m = re.search(r'f"radius \{radius:\.2f\} -> area \{total:\.4f\}"', s)
    out["string_split_in_two"] = (
        s[:m.start()] + 'f"radius {radius:.2f} -> "\n              '
        'f"area {total:.4f}"' + s[m.end():], False)
    return out


def selftest() -> int:
    fails: list[str] = []
    muts = _mutations()
    print(f"self-test: {len(muts)} cases "
          f"({sum(1 for _, r in muts.values() if r)} must be rejected, "
          f"{sum(1 for _, r in muts.values() if not r)} must be accepted)")
    for name, (mutated, must_reject) in sorted(muts.items()):
        # ★ Assert the mutation actually changed the file BEFORE asking the
        #   checker about it. The first version of this self-test was worthless
        #   because a mutation silently failed to apply and the checker's "OK"
        #   was read as a pass.
        if mutated == _SAMPLE:
            fails.append(f"{name}: mutation did not change the source -- "
                         f"this case tested NOTHING")
            print(f"  ✗ {name:22s} mutation was a no-op")
            continue
        rejected = bool(compare(_SAMPLE, mutated))
        ok = rejected == must_reject
        if not ok:
            fails.append(f"{name}: expected "
                         f"{'reject' if must_reject else 'accept'}, "
                         f"got {'reject' if rejected else 'accept'}")
        verb = "rejected" if rejected else "accepted"
        print(f"  {'✓' if ok else '✗'} {name:22s} {verb}"
              f"{'' if ok else '   <-- WRONG'}")
    print()
    if fails:
        print(f"✗ FAIL — {len(fails)}/{len(muts)} self-test cases wrong")
        for f in fails:
            print("    " + f)
        return 1
    print(f"✓ PASS — {len(muts)}/{len(muts)}. The checker rejects code changes, "
          f"accepts comment and string changes.")
    return 0


# ════════════════════════════════════════════════════════════════════════
# Applier -- exactly-once substitutions, refusing to write unless safe
# ════════════════════════════════════════════════════════════════════════
def parse_subs(text: str) -> tuple[list[dict], list[str]]:
    """Parse the marker format: --OLD / <text> / --NEW / <text>, then --KEEP."""
    pairs: list[dict] = []
    keep: list[str] = []
    cur: dict | None = None
    mode: str | None = None
    for line in text.split("\n"):
        if line == "--OLD":
            if cur:
                pairs.append(cur)
            cur, mode = {"old": [], "new": []}, "old"
        elif line == "--NEW" and cur is not None:
            mode = "new"
        elif line == "--KEEP":
            if cur:
                pairs.append(cur)
                cur = None
            mode = "keep"
        elif mode == "keep":
            if line.strip():
                keep.append(line.strip())
        elif cur is not None and mode:
            cur[mode].append(line)
    if cur:
        pairs.append(cur)
    for p in pairs:
        for k in ("old", "new"):
            # A heredoc leaves a trailing empty line before the next marker;
            # without dropping it the LAST pair silently gains a newline, which
            # inserts a stray blank line into the target (this happened once).
            if p[k] and p[k][-1] == "":
                p[k] = p[k][:-1]
            p[k] = "\n".join(p[k])
    return pairs, keep


def apply_subs(target: Path, subsfile: Path) -> int:
    original = target.read_text(encoding="utf-8")
    pairs, keep = parse_subs(subsfile.read_text(encoding="utf-8"))
    src, errs = original, []
    for i, p in enumerate(pairs, 1):
        count = src.count(p["old"])
        if count != 1:
            errs.append(f"substitution #{i}: OLD occurs {count} times "
                        f"(needs exactly 1)\n    {p['old'][:110]!r}")
            continue
        src = src.replace(p["old"], p["new"], 1)
    if errs:
        print(f"ABORT {target} — {len(errs)} problem(s), nothing written:")
        for e in errs:
            print("  " + e)
        return 1

    stragglers = [ln for ln in src.split("\n")
                  if KO.search(ln) and not any(k in ln for k in keep)]
    if stragglers:
        print(f"ABORT {target} — {len(stragglers)} Hangul line(s) left, "
              f"nothing written:")
        for ln in stragglers[:12]:
            print(f"    {ln.strip()[:110]!r}")
        return 1

    problems = compare(original, src)
    if problems:
        print(f"ABORT {target} — code would change, nothing written:")
        for p in problems:
            print("  " + p)
        return 1

    target.write_text(src, encoding="utf-8")
    print(f"WROTE {target}: {len(pairs)} substitutions, code tokens identical"
          + (f", {len(keep)} deliberate Hangul kept" if keep else ""))
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--selftest", action="store_true",
                    help="break the checker on purpose and confirm it notices")
    ap.add_argument("--apply", nargs=2, metavar=("TARGET", "SUBS"),
                    help="apply a substitution file, refusing to write if unsafe")
    ap.add_argument("files", nargs="*", metavar="BEFORE AFTER",
                    help="two paths to compare")
    args = ap.parse_args()

    if args.selftest:
        return selftest()
    if args.apply:
        return apply_subs(Path(args.apply[0]), Path(args.apply[1]))
    if len(args.files) == 2:
        return check_files(Path(args.files[0]), Path(args.files[1]))
    ap.print_help()
    return 2


if __name__ == "__main__":
    sys.exit(main())
