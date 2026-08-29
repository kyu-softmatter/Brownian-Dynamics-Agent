"""The knowledge base -- reads both branches together.

    runs/*/record.json          run experience (written by tools/postmortem.py)
    knowledge/entries/*.json    * knowledge with no run -- sketch readings,
                                literature distillations, tooling lessons

A SQLite + FTS5 store is the plan, but it is overkill below 100 runs.
**Fix the data format first** and move to a DB when it hurts.

* Why `entries/` is needed: only lessons that came out of runs were being stored.
  But **knowledge produced by the front end (sketch -> physical system) is reusable
  too, and forgetting it is a loss.** For example: "inheriting the R=5um convention
  into another sketch disagreed with abp-rod's tau_R by 160x"; "a paper extraction
  is verified by reproducing the paper's own numbers"; "patching YAML with a regex
  silently edits the wrong node."
  Those have no run, so there was no slot for them in record.json.

WARNING: this is one of **two unmerged knowledge schemas.** `knowledge/wiki/` is
human-written Markdown read by the `bd-knowledge` skill; this store is
tool-written JSON read by this file. A lesson filed in one is invisible to a
reader of the other -- query both. See docs/03-knowledge-base.md.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY tools/kb.py list
    $PY tools/kb.py query --tags 2D,harmonic_trap
    $PY tools/kb.py query --kind intake
    $PY tools/kb.py lessons
    $PY tools/kb.py add --claim "..." --kind pitfall --origin intake \
                        --source "intake/abp-rod/observation.yaml#D1" --tags abp,convention
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


# 2026-08-29 (merge fix): this pointed at `kb/entries`, which the merge renamed
# to `knowledge/entries`. The result was not an error -- `kb.py list` simply
# reported "run-less knowledge 0" for 126 existing entries. A silently empty
# read is the same failure mode as an unwired checker.
ENTRY_DIR = ROOT / "knowledge" / "entries"
SCHEMA_ENTRY = "bdbot.kb_entry/0.1"
# Source kinds for non-run knowledge
ORIGINS = ("intake", "paper", "tooling", "method", "handbook", "user_input")


def load_all():
    """Return run records and run-less entries merged into the same shape."""
    out = []
    for f in sorted((ROOT / "runs").glob("*/record.json")):
        try:
            r = json.loads(f.read_text())
            r.setdefault("origin", "our_run")
            out.append(r)
        except Exception as e:
            print(f"  (skipped {f}: {e})", file=sys.stderr)
    for f in sorted(ENTRY_DIR.glob("*.json")):
        try:
            out.append(json.loads(f.read_text()))
        except Exception as e:
            print(f"  (skipped {f}: {e})", file=sys.stderr)
    return out


def add_entry(claim, kind, origin, source, tags, coords, not_verified) -> Path:
    """Write one run-less knowledge entry as a single file.

    Uses the same field names as record.json.
    """
    if origin not in ORIGINS:
        raise SystemExit(f"origin must be one of {ORIGINS} (got: {origin})")
    ENTRY_DIR.mkdir(parents=True, exist_ok=True)
    slug = "".join(c if c.isalnum() or c in "-_" else "-" for c in claim.lower())[:48]
    slug = "-".join(x for x in slug.split("-") if x) or "entry"
    path = ENTRY_DIR / f"{origin}__{slug}.json"
    n = 2
    while path.exists():
        path = ENTRY_DIR / f"{origin}__{slug}-{n}.json"
        n += 1
    rec = {
        "schema": SCHEMA_ENTRY,
        "run_id": None,                 # there is no run -- that is why this entry exists
        "origin": origin,
        "case": None,
        "outcome": None,
        "system_tags": tags,
        "dimensionless": coords,
        "source": source,
        "lessons": [{"claim": claim, "kind": kind, "coords": coords}],
        "not_verified": not_verified,
    }
    path.write_text(json.dumps(rec, indent=2, ensure_ascii=False))
    return path


def match_coord(rec, expr):
    """'key=lo:hi' or 'key=value' (10% tolerance). The key matches partially."""
    key, rng = expr.split("=", 1)
    key = key.strip()
    hits = [(k, v) for k, v in rec.get("dimensionless", {}).items() if key in k]
    if not hits:
        return False
    _, val = hits[0]
    if ":" in rng:
        lo, hi = (float(x) for x in rng.split(":"))
        return lo <= val <= hi
    t = float(rng)
    return abs(val - t) <= 0.1 * abs(t)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    q = sub.add_parser("query")
    q.add_argument("--tags", default="")
    q.add_argument("--coord", action="append", default=[])
    q.add_argument("--outcome", default="")
    q.add_argument("--kind", default="", help="lesson kind (pitfall/method_note/...)")
    q.add_argument("--origin", default="", help=f"source kind {ORIGINS + ('our_run',)}")
    sub.add_parser("lessons").add_argument("--origin", default="")
    a_ = sub.add_parser("add", help="add a run-less knowledge entry (sketch reading, literature, tooling lesson)")
    a_.add_argument("--claim", required=True)
    a_.add_argument("--kind", default="method_note",
                    help="parameter | phase_boundary | scaling | method_note | pitfall")
    a_.add_argument("--origin", required=True, help=f"{ORIGINS}")
    a_.add_argument("--source", default="", help="where it came from (file#anchor, or a paper locator)")
    a_.add_argument("--tags", default="")
    a_.add_argument("--coord", action="append", default=[],
                    help="dimensionless coordinates 'name=value' (for search)")
    a_.add_argument("--not-verified", action="append", default=[])
    a = ap.parse_args()

    if a.cmd == "add":
        coords = {}
        for c in a.coord:
            k, v = c.split("=", 1)
            coords[k.strip()] = float(v)
        path = add_entry(a.claim, a.kind, a.origin, a.source,
                         [t.strip() for t in a.tags.split(",") if t.strip()],
                         coords, a.not_verified)
        print(f"→ {path.relative_to(ROOT)}")
        return 0

    recs = load_all()
    if not recs:
        print("no records. Use tools/postmortem.py (for a run) or kb.py add (for run-less knowledge).")
        return 1

    def rid(r):
        return r.get("run_id") or f"[{r.get('origin', '?')}] {r.get('source') or '—'}"

    if a.cmd == "list":
        print(f"{'source':<44}{'origin':<11}{'outcome':<10}{'tags'}")
        print("-" * 104)
        for r in sorted(recs, key=lambda x: (x.get("origin") != "our_run", rid(x))):
            print(f"{rid(r)[:43]:<44}{r.get('origin', '?'):<11}"
                  f"{str(r.get('outcome') or '—'):<10}"
                  f"{','.join((r.get('system_tags') or [])[:4])}")
        n_run = sum(1 for r in recs if r.get("origin") == "our_run")
        print(f"\n{len(recs)} total  (runs {n_run} . run-less knowledge {len(recs) - n_run})")
        return 0

    if a.cmd == "lessons":
        n = 0
        for r in recs:
            if a.origin and r.get("origin") != a.origin:
                continue
            for l_ in r.get("lessons", []):
                n += 1
                c = f"   coords={l_['coords']}" if l_.get("coords") else ""
                tier = l_.get("tier", r.get("tier", 3))
                print(f"[tier{tier} {l_['kind']}] {l_['claim']}{c}")
                print(f"      ← {rid(r)}"
                      + (f" ({r['outcome']})" if r.get("outcome") else
                         f" · origin={r.get('origin')}"))
        print(f"\n{n} total")
        return 0

    sel = recs
    if a.tags:
        want = {t.strip() for t in a.tags.split(",")}
        sel = [r for r in sel if want <= set(r.get("system_tags") or [])]
    for c in a.coord:
        sel = [r for r in sel if match_coord(r, c)]
    if a.outcome:
        sel = [r for r in sel if r.get("outcome") == a.outcome]
    if a.origin:
        sel = [r for r in sel if r.get("origin") == a.origin]
    if a.kind:
        sel = [r for r in sel
               if any(l_.get("kind") == a.kind for l_ in r.get("lessons", []))]

    for r in sel:
        print("=" * 78)
        print(f"{rid(r)}   [{r.get('outcome') or r.get('origin')}]")
        if r.get("system_tags"):
            print(f"  tags: {', '.join(r['system_tags'])}")
        if r.get("dimensionless"):
            print("  dimensionless groups:")
            for k, v in r["dimensionless"].items():
                print(f"    {k:<36} {v:.4g}")
        if r.get("observables"):
            print("  observables (measured vs predicted):")
            for o in r["observables"]:
                pv = o.get("predicted")
                if pv is None:
                    print(f"    {o['name']:<16} {o['measured']:>13.6g} {o['unit']:<8} (no prediction)")
                else:
                    print(f"    {o['name']:<16} {o['measured']:>13.6g} vs {pv:>13.6g} "
                          f"{o['unit']:<8} {o['err_pct']:+7.2f}%")
        if r.get("lessons"):
            print("  lessons:")
            for l_ in r["lessons"]:
                print(f"    [{l_['kind']}] {l_['claim']}")
        if r.get("not_verified"):
            print("  not verified:")
            for nv in r["not_verified"]:
                print(f"    · {nv}")
    print(f"\n{len(sel)}/{len(recs)} matched")
    return 0


if __name__ == "__main__":
    sys.exit(main())
