"""Mutation test at the **gate** level — for each enforcement point, remove it
and count how many tests notice.

    $PY verify/verify_gates_bite.py            # full suite per mutant (~20 min)
    $PY verify/verify_gates_bite.py --fast     # declared subset per mutant
    $PY verify/verify_gates_bite.py --only L2-dim-wired

**What this measures, precisely.** The column is *how many tests fail when the
gate is deleted* — i.e. **suite coverage of the gate**. It is NOT whether the
gate is correct. A gate can be perfectly written and covered by nothing, and
that is the failure this repository keeps hitting:

  · `health.gate()` is reachable only from `tools/health.py`; `run.execute` never
    called it, so **no run has ever gated itself** — and because it was unwired,
    the fact that it tested `verdict != "PASS"` (80 false rejections out of 83
    specs) went unnoticed for a month. An unwired checker cannot be wrong out loud
  · `A4` was a `grep` that matched the prose explaining the rule (7 false hits)
  · `bd-intake` section 2.1's empty-goal blocker was written twice, enforced zero
    times, and walked past by 2 of 8 cases which then produced 85 runs

So `expect="uncaught"` rows are not bugs in this script. They are the measured,
named gaps, and they are the reason the table exists.

**The guard that makes the number mean anything.** A `sed` that matched nothing
produces "0 tests failed", which is indistinguishable from a gate the suite
ignores — docs/05 section 2, *"a check whose success is indistinguishable from
its own failure."* So every mutation asserts its anchor occurs **exactly once**
before it is applied, and asserts the file actually changed after. A mutation
that cannot be applied is an ERROR, never a 0.

Mutates the working tree in place and restores with `git checkout --`. Refuses
to start if any target file is already dirty.
"""
from __future__ import annotations

import argparse
import json
import pathlib
import re
import subprocess
import sys
import time

ROOT = pathlib.Path(__file__).resolve().parent.parent

# ── the enforcement points ─────────────────────────────────────────────────
# `kind`:   replace | inject
# `expect`: caught   -> at least one test must fail. A survivor is a finding.
#           uncaught -> documented gap. Becoming caught is also reported, because
#                       a stale table is the thing this file exists to prevent.
MUTANTS = [
    dict(id="L2-dim-wired", file="bdbot/physical.py", kind="replace",
         anchor="    out += check_dim(s)", into="    pass  # MUTANT",
         gate="physical.validate calls check_dim", expect="caught",
         fast=["tests/test_dim_gate.py"]),
    dict(id="L2-size-wired", file="bdbot/physical.py", kind="replace",
         anchor="    out += check_size(s)", into="    pass  # MUTANT",
         gate="physical.validate calls check_size", expect="caught",
         fast=["tests/test_size_gate.py"]),
    dict(id="RUN-structure-validate", file="bdbot/run.py", kind="replace",
         anchor='        _errs = [i for i in _iss if i.level == "error"]',
         into="        _errs = []  # MUTANT",
         gate="execute() refuses a structure block that does not validate",
         expect="caught", fast=["tests/test_size_gate.py", "tests/test_dim_gate.py"]),
    dict(id="RUN-structure-required", file="bdbot/run.py", kind="replace",
         anchor="    elif require_structure:", into="    elif False:  # MUTANT",
         gate="require_structure=True refuses a spec with no block",
         expect="caught", fast=["tests/test_size_gate.py"]),
    dict(id="RUN-spec-hash", file="bdbot/run.py", kind="replace",
         anchor="    if not ok_hash:", into="    if False:  # MUTANT",
         gate="execute() refuses a hand-edited spec (rule 2)",
         expect="caught", fast=["tests/test_run_gates.py"]),
    dict(id="RUN-verdict-fail", file="bdbot/run.py", kind="replace",
         anchor='    if spec.verdict.startswith("FAIL"):',
         into="    if False:  # MUTANT",
         gate="execute() refuses a spec whose L3 verdict is FAIL",
         expect="caught", fast=["tests/test_run_gates.py"]),
    dict(id="RUN-approval-missing", file="bdbot/run.py", kind="replace",
         anchor="    elif require_approval:", into="    elif False:  # MUTANT",
         gate="rule 10: require_approval=True refuses a run with no params.json",
         expect="caught", fast=["tests/test_run_gates.py"]),
    dict(id="RUN-approval-unsigned", file="bdbot/run.py", kind="replace",
         anchor="        if require_approval and not man.approved_by:",
         into="        if False:  # MUTANT",
         gate="rule 10: an unapproved params.json does not run",
         expect="caught", fast=["tests/test_run_gates.py"]),
    dict(id="RUN-params-drift", file="bdbot/run.py", kind="replace",
         anchor="        if drift:", into="        if False:  # MUTANT",
         gate="approving one set of numbers and running another",
         expect="caught", fast=["tests/test_run_gates.py"]),
    dict(id="SEAL-broken", file="bdbot/runcard.py", kind="replace",
         anchor="    if not ok:", into="    if False:  # MUTANT",
         gate="a broken seal is a hard stop", expect="caught",
         fast=["tests/test_run_gates.py"]),
    dict(id="SEAL-missing", file="bdbot/runcard.py", kind="replace",
         anchor="        if require:", into="        if False:  # MUTANT",
         gate="require_seal=True refuses an unsealed run", expect="caught",
         fast=["tests/test_run_gates.py"]),
    dict(id="A4-core-imports-llm", file="bdbot/checks.py", kind="inject",
         anchor="import anthropic  # MUTANT",
         gate="A4: the deterministic core carries no LLM dependency",
         expect="caught", fast=["tests/test_invariants.py"]),
    dict(id="RULE10-derived-recomputed", file="bdbot/params.py", kind="replace",
         anchor="            out += self.check_derived()",
         into="            pass  # MUTANT",
         gate="rule 10: a `derived` value must follow from d, T and eta",
         expect="caught",
         fast=["tests/test_params.py"]),
    dict(id="PRESERVE-gate-inputs", file="bdbot/runid.py", kind="replace",
         anchor="""    "record.json",             # a lesson must outlive the run artefacts
    "params.json",             # rule 10: the approved numbers
    SEAL_FILE,                 # the pre-registration""",
         into="    \"record.json\",  # MUTANT",
         gate="prepare_outdir keeps the human-authored gate inputs",
         expect="caught",
         fast=["tests/test_run_gates.py"]),
    dict(id="HEALTH-gate-wired", file="bdbot/health.py", kind="replace",
         anchor="    problems = []\n    ok, want = spec.verify_hash()",
         into="    return []  # MUTANT\n    problems = []\n    ok, want = spec.verify_hash()",
         gate="health.gate() blocks a bad spec",
         expect="uncaught",
         why="documented in docs/05 section 2 and in health.gate's own docstring: "
             "reachable only from tools/health.py, so run.execute never calls it",
         fast=["tests/test_blocks_cdfgk.py"]),
]


def sh(*args, **kw):
    return subprocess.run(args, cwd=ROOT, capture_output=True, text=True, **kw)


def snapshot(paths) -> dict:
    """Read every target file once, up front.

    ★ Restoring with `git checkout --` was wrong twice over: it required a
      clean tree, so the harness could not verify a fix before that fix was
      committed — exactly when you most want to run it — and it would have
      silently reverted uncommitted work to HEAD. Snapshotting the bytes has
      neither problem and drops the git dependency.
    """
    return {f: (ROOT / f).read_text() for f in sorted(set(paths))}


def apply_mutant(m) -> None:
    """Apply, or raise. Never silently no-op -- see the module docstring."""
    p = ROOT / m["file"]
    before = p.read_text()
    if m["kind"] == "inject":
        after = m["anchor"] + "\n" + before
    else:
        n = before.count(m["anchor"])
        if n != 1:
            raise RuntimeError(
                f"{m['id']}: anchor occurs {n}x in {m['file']}, need exactly 1. "
                f"A mutation that cannot be placed must not be reported as 0.")
        after = before.replace(m["anchor"], m["into"])
    if after == before:
        raise RuntimeError(f"{m['id']}: file unchanged after mutation")
    p.write_text(after)


def restore(snap: dict) -> None:
    for f, text in snap.items():
        p = ROOT / f
        if p.read_text() != text:
            p.write_text(text)


#: ⚠ Parse pytest's FINAL SUMMARY LINE only. Summing every `N failed` / `N
#: error` across the whole capture double-counts: the A4 mutant printed "1
#: error in 1.41s" plus an `ERROR tests/...` line and was reported as **2**.
#: A paper table with an inflated count is the defect this file is about.
SUMMARY = re.compile(r"(\d+) (failed|error|errors)\b")


def run_tests(targets) -> tuple[int, str]:
    cmd = [sys.executable, "-m", "pytest", "-q", "-p", "no:cacheprovider"]
    cmd += list(targets) if targets else []
    r = sh(*cmd)
    lines = [ln for ln in (r.stdout or "").splitlines() if ln.strip()]
    last = lines[-1] if lines else ""
    n = sum(int(m.group(1)) for m in SUMMARY.finditer(last))
    return n, last.strip()


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--fast", action="store_true",
                    help="run each mutant's declared test subset, not the suite")
    ap.add_argument("--only", action="append", default=[])
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    sel = [m for m in MUTANTS if not a.only or m["id"] in a.only]
    snap = snapshot([m["file"] for m in sel])

    scope = "declared subset" if a.fast else "full suite"
    print(f"{len(sel)} gates - scope: {scope} - {sys.executable}")
    print()

    rows, bad = [], []
    for i, m in enumerate(sel, 1):
        t0 = time.time()
        print(f"[{i}/{len(sel)}] {m['id']:24s} ", end="", flush=True)
        try:
            apply_mutant(m)
            n, last = run_tests(m["fast"] if a.fast else None)
            status = "ERROR" if False else None
        except RuntimeError as exc:
            restore(snap)
            print(f"CANNOT APPLY - {exc}")
            rows.append(dict(id=m["id"], gate=m["gate"], caught=None,
                             expect=m["expect"], note=str(exc)))
            bad.append(m["id"])
            continue
        finally:
            restore(snap)

        dt = time.time() - t0
        caught = n > 0
        ok = caught == (m["expect"] == "caught")
        mark = "ok " if ok else "!! "
        print(f"{n:4d} tests fail  ({dt:.0f}s)  {mark}{last}")
        if not ok:
            bad.append(m["id"])
        rows.append(dict(id=m["id"], gate=m["gate"], failures=n, caught=caught,
                         expect=m["expect"], agrees=ok, why=m.get("why"),
                         seconds=round(dt, 1)))

    print()
    print("| gate | removing it fails | expected | |")
    print("|---|---:|---|---|")
    for r in rows:
        n = "n/a" if r.get("failures") is None else f"**{r['failures']}**"
        exp = r["expect"]
        mk = "✓" if r.get("agrees") else "✗"
        print(f"| {r['gate']} | {n} | {exp} | {mk} |")

    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(
            dict(scope=scope, rows=rows), indent=1, ensure_ascii=False))
        print(f"\nwrote {a.json}")

    print()
    if bad:
        print(f"DISAGREES WITH THE RECORD: {bad}")
        print("A survivor of an `expect=caught` gate is an unenforced gate.")
        print("An `expect=uncaught` gate that is now caught means this table is "
              "stale -- update it, that is good news.")
        return 1
    print(f"all {len(rows)} gates agree with the record")
    return 0


if __name__ == "__main__":
    sys.exit(main())
