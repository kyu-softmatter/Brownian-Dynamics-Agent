"""Which of bdbot's error paths has anything ever been seen to fire?

The general form that came out of the bridge collaboration (2026-09-16), applied
to this repository: **enumerate every error-producing site and list the ones with
neither a negative test nor a live firing.** Applied to the bridge it returned
four rules -- R0, R1, R9, R10 -- all of which passed once fixtures existed, so
none had been dead; nobody had ever watched them run. R1's error path in
particular had never executed, because every fixture tripped a provenance rule
and none tripped the schema.

CLAUDE.md already carries the practice ("when you build a checker, deliberately
break it and see") and three instances of what happens when nobody does. This is
the procedure form of it, which is the part that does not need a reader to
notice: discovery needs a reader, coverage does not.

    $PY -m pytest tests/ -q -p covplug     # writes /tmp/bdcov.json (see verify/_covplug.py)
    $PY verify/verify_error_paths_fire.py
"""
from __future__ import annotations
import ast, json, pathlib, sys
from collections import Counter

ROOT = pathlib.Path(__file__).resolve().parent.parent
COV = pathlib.Path("/tmp/bdcov.json")

# Gates whose failure is load-bearing: a run is refused, a seal is broken, a
# number is rejected. An unfired error path here is the unwired-checker shape.
LOAD_BEARING = {"params", "runcard", "runid", "physical", "metrics", "checks",
                "health", "intake", "goal", "nondim", "scales", "design"}


def error_lines(path: pathlib.Path) -> dict[int, str]:
    """{lineno: short description} for every raise / Issue('error') / .err() site."""
    src = path.read_text()
    tree = ast.parse(src)
    lines = src.splitlines()
    out = {}
    for n in ast.walk(tree):
        hit = isinstance(n, ast.Raise) or (
            isinstance(n, ast.Call) and getattr(n.func, "attr", None) in ("err",)) or (
            isinstance(n, ast.Call) and getattr(n.func, "id", None) == "Issue"
            and n.args and isinstance(n.args[0], ast.Constant) and n.args[0].value == "error")
        if not hit:
            continue
        ln = n.lineno
        out[ln] = lines[ln - 1].strip()[:84]
    return out


def main():
    if not COV.exists():
        print(f"no {COV}. Run the suite under the plugin first:\n"
              f"  $PY -m pytest tests/ -q -p covplug")
        return 2
    cov = {pathlib.Path(k).resolve(): set(v) for k, v in json.loads(COV.read_text()).items()}

    total = unfired = 0
    rows = []
    for f in sorted((ROOT / "bdbot").glob("*.py")):
        el = error_lines(f)
        if not el:
            continue
        seen = cov.get(f.resolve(), set())
        for ln, txt in sorted(el.items()):
            total += 1
            if ln not in seen:
                unfired += 1
                rows.append((f.stem, ln, txt))

    print("=" * 96)
    print("bdbot error paths never executed by the test suite")
    print("=" * 96)
    print(f"  {total} error sites, {unfired} never fired ({100*unfired/total:.0f} %)")

    by_mod = Counter(r[0] for r in rows)
    print(f"\n  {'module':<16} {'unfired':>8}   load-bearing?")
    for mod, k in by_mod.most_common():
        print(f"  {mod:<16} {k:8d}   {'YES' if mod in LOAD_BEARING else ''}")

    print(f"\n--- the load-bearing ones, which are the unwired-checker candidates ---")
    n = 0
    for mod, ln, txt in rows:
        if mod in LOAD_BEARING:
            n += 1
            print(f"  bdbot/{mod}.py:{ln}  {txt}")
    print(f"\n  {n} unfired error paths in gate modules.")
    print("  An unfired path is not necessarily dead -- the bridge's four all passed once")
    print("  fixtures existed. It means nobody has ever seen it run, which is the state in")
    print("  which `chain`, `A4`'s grep and `health.gate()` all sat while looking fine.")
    (ROOT / "verify/_out/error_paths.json").write_text(json.dumps(
        {"total": total, "unfired": unfired,
         "rows": [{"module": m, "line": l, "src": t} for m, l, t in rows]}, indent=1))
    print(f"\nwrote verify/_out/error_paths.json")
    return 0


if __name__ == "__main__":
    sys.exit(main())
