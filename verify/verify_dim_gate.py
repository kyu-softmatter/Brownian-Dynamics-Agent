"""Adversarial test of the `structure.dim` gate — does it bite, and does it let
a correct answer through?

`knowledge/wiki/concepts/dimensionality-has-no-default.md` closes `D9`. This
script is what stops that decision from being the fourth entry in rule 10's list
of practices that were written and never enforced. Two failure modes have to be
excluded, and this repository has hit both:

  · **an unwired checker cannot be wrong out loud** — the step-resolution check
    silently never ran across 81 runs
  · **a gate that refuses everything is worse than none** — the pre-run gate once
    rejected 80 of 83 specs with zero real failures among them

Section ② feeds it broken input one field at a time. Section ③ requires the eight
**real** cases to pass — they are the strongest possible "does it let a correct
answer through", because they were not written to satisfy the checker. Section ①
is the regression that matters most to the archive.

    $PY verify/verify_dim_gate.py
"""
import copy
import hashlib
import json
import pathlib
import sys
import tempfile

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import physical as P  # noqa: E402
from bdbot import runid as RID  # noqa: E402

DONOR = ROOT / "intake/soft-r3-2d-A-sweep"

results = []


def report(ok, label, detail=""):
    results.append(bool(ok))
    print(f"  {'OK ' if ok else 'XX '} {label:<50}{detail}")


# ══════════════════════════════════════════════════════════════════════
print("=" * 78)
print("(1) run_id stability -- the archive must not be renamed")
print("=" * 78)
OLD = frozenset(k for k in RID.DOC_KEYS if k != "structure")


def _strip(node, keys):
    if isinstance(node, dict):
        return {k: _strip(v, keys) for k, v in node.items() if k not in keys}
    if isinstance(node, list):
        return [_strip(v, keys) for v in node]
    return node


def _h(d):
    return hashlib.sha256(
        json.dumps(d, sort_keys=True, default=str).encode()).hexdigest()[:12]


specs = sorted((ROOT / "specs").glob("*.json"))
moved = [f.name for f in specs
         if _h(_strip(json.load(open(f)), OLD)) != _h(_strip(json.load(open(f)),
                                                             RID.DOC_KEYS))]
report(specs and not moved, f"{len(specs)} archived specs keep their hash",
       f"{len(moved)} moved" if moved else "0 moved")

#  ...and the exclusion is not vacuous: a doc that DOES carry the key must change
probe = {"dimensions": 2, "structure": {"dim": {"value": 2}}}
report(_h(_strip(probe, OLD)) != _h(_strip(probe, RID.DOC_KEYS)),
       "  -> and the exclusion is real, not vacuous",
       "a structure key does change the stripped form")

#  ...while the physics value stays hashed: 2 -> 3 must still re-id
report(_h(_strip({"dimensions": 2}, RID.DOC_KEYS))
       != _h(_strip({"dimensions": 3}, RID.DOC_KEYS)),
       "  -> and `dimensions` 2 -> 3 still re-ids the run", "both directions held")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("(2) break it one field at a time -- every rule must fire")
print("=" * 78)
base_sys = yaml.safe_load((DONOR / "system.yaml").read_text())
GOOD = copy.deepcopy(base_sys["structure"])          # the donor's real, passing block

tmp = pathlib.Path(tempfile.mkdtemp())
(tmp / "sketch_01.jpeg").write_bytes(b"")
(tmp / "observation.yaml").write_text((DONOR / "observation.yaml").read_text())


def check(mutate, label, where_contains, want_error=True):
    d = copy.deepcopy(base_sys)
    d["structure"] = copy.deepcopy(GOOD)
    mutate(d)
    (tmp / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    hits = [i for i in P.load(tmp).errors if where_contains in i.where]
    report(bool(hits) == want_error, label, (hits[0].msg[:40] if hits else "NOT CAUGHT"))


check(lambda d: d.pop("structure"), "no `structure:` section at all", "structure.dim")
check(lambda d: d["structure"].pop("dim"), "`structure.dim` missing", "structure.dim")
check(lambda d: d["structure"]["dim"].pop("basis"), "`basis` missing", "structure.dim.basis")
check(lambda d: d["structure"]["dim"].update(basis="because"),
      "`basis` not one of the four", "structure.dim.basis")
check(lambda d: d["structure"]["dim"].pop("what_would_change"),
      "* `what_would_change` missing", "structure.dim.what_would_change")
check(lambda d: d["structure"]["dim"].update(what_would_change="   "),
      "* `what_would_change` present but blank", "structure.dim.what_would_change")
check(lambda d: d["structure"]["dim"].update(alternatives=[]),
      "`alternatives` empty list", "structure.dim.alternatives")
check(lambda d: d["structure"]["dim"].update(value=3),
      "* dim says 3 while `dimensions:` says 2", "structure.dim.value")

print()
print("  -- the per-basis evidence each owes --")
check(lambda d: d["structure"]["dim"].update(basis="given"),
      "basis `given` without `source`", "structure.dim.source")
check(lambda d: d["structure"]["dim"].update(basis="inherited"),
      "basis `inherited` without `compared_with`", "structure.dim.compared_with")
check(lambda d: d["structure"]["dim"].update(basis="sufficient"),
      "basis `sufficient` without `checked_observables`",
      "structure.dim.checked_observables")
check(lambda d: d["structure"]["dim"].update(basis="required"),
      "basis `required` owes nothing extra -> no error",
      "structure.dim", want_error=False)

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("(3) the eight real cases pass -- they were not written for the checker")
print("=" * 78)
cases = sorted(p.parent for p in ROOT.glob("intake/*/system.yaml"))
report(len(cases) == 8, f"{len(cases)} cases on disk", "expected 8")
dist, expiring = {}, []
for c in cases:
    s = P.load(c)
    se = [i for i in s.errors if i.where.startswith("structure")]
    b = ((s.raw.get("structure") or {}).get("dim") or {}).get("basis")
    dist[b] = dist.get(b, 0) + 1
    if any("EXPIRES" in i.msg for i in s.issues):
        expiring.append(c.name)
    report(not se, f"{c.name} -> {s.dim}D  basis={b}",
           f"{len(se)} error(s)" if se else "clean")

print()
report(sum(dist.values()) == 8 and None not in dist,
       "every case declares a basis", " . ".join(f"{k}={v}" for k, v in sorted(dist.items())))
report(len(expiring) == 3,
       "the expiring bases say so, every read",
       f"{len(expiring)}: {', '.join(expiring)}"[:56])

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("(4) the gate is WIRED -- `validate` calls it, not a sibling tool")
print("=" * 78)
d = copy.deepcopy(base_sys)
d.pop("structure")
(tmp / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
report(any(i.where.startswith("structure") for i in P.load(tmp).errors),
       "reached through physical.load -> validate", "not only via check_dim()")

print()
print("=" * 78)
n_ok = sum(results)
print(f"{n_ok}/{len(results)} passed")
print("=" * 78)
sys.exit(0 if n_ok == len(results) else 1)
