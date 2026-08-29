"""Phase 1-C equivalence regression test -- did the refactor change any result?

Phase 1-C's DoD: **1-A and 1-B run through the same code path, and the results are
identical to before.**
`verify/ref_1b/*.metrics.json` are snapshots of the results from **before** the
refactor (Phase 1-B).

⛔ **THE BASELINE IS MISSING, SO THIS TEST CANNOT RUN** (measured 2026-08-29).
   `verify/ref_1b/` is empty and `git log --all -- '*/ref_1b/*'` returns nothing:
   the snapshots were never committed, so no clone has ever been able to run this.
   It used to fail with a bare FileNotFoundError; it now says so explicitly.

   ⚠️ **Do NOT regenerate the snapshots from the current code.** A snapshot taken
   now would make this a comparison of current-against-current, which passes
   unconditionally and guards nothing -- converting a visibly broken check into a
   silently passing one, which is strictly worse. The baseline is only meaningful
   if it predates the 1-C refactor. Either recover those two files from outside the
   repository, or retire this script; do not manufacture agreement.

Three things are checked:
  (1) run_id  -- the spec hash. If equal, the ledger -> checks -> dt path is
      byte-identical.
      ⚠️ 1-A's run_id **was deliberately changed once, in earlier work** (see
         EXPECTED_RUN_ID below). Cause: the spec hashed the whole of system.yaml, so
         merely adding a `derived_from` field invalidated it. Fixed to hash only the
         physics fields (bdbot.runid.physics_only). All 96 physics values were
         confirmed unchanged.
  (2) physics values -- do the observables, dimensionless groups and check values
      match the snapshot?
  (3) schema -- state which fields 1-C added (deletions are not permitted)

    $PY verify/verify_1c_equivalence.py
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
REF = ROOT / "verify/ref_1b"

TARGETS = [
    ("1-A trap-2d-5um", "cases/trap_2d_5um.py", [], "trap-2d-5um.metrics.json"),
    ("1-B soft-r3 A=100", "cases/soft_r3_2d.py", ["--A", "100"], "soft-r3-A100.metrics.json"),
]

# Cases where a run_id differing from the snapshot is **correct** -- a reason has to
# be written down before it is allowed to pass.
EXPECTED_RUN_ID = {
    "trap-2d-5um__70b9394e7310": (
        "trap-2d-5um__a5ef4f45d589",
        "(1) earlier work #2 changed the hash basis from 'the whole YAML' to "
        "'physics fields only' "
        "(-> a49f2508556b). Invalidating a run merely by adding a `derived_from` "
        "field is something "
        "which defeats the purpose of content addressing. "
        "(2) the L3 work unified the spec schema across the three cases "
        "(`bdbot.nondim.NondimSpec`, -> d724e8d507cc). The hash payload went from a "
        "flat dict to {system, params, numerics} -- **the same physics in a "
        "different arrangement**. There was no defect in 1-A itself; this is a side "
        "effect of unifying the schema. "
        "(3) the L5 migration (2026-08-05) added `numerics.seed` to the spec "
        "(-> a5ef4f45d589) -- `RUN.builder`'s `build(spec)` does not re-read the case "
        "YAML, so seed=20260803, which previously lived only inside main(), had to go "
        "into the spec for build() to see it. The re-run was confirmed to match the "
        "old run to 15 decimal places on all 5 observables."),
    "soft-r3-2d-A-sweep__A100__27f70deab9": (
        "soft-r3-2d-A-sweep__A100__30caa5c9e0",
        "⭐️ **A DEFECT FIX.** The 1-B spec had **no physical system in it** -- "
        "changing d from 5µm to 0.5µm, eta by 62x, or rho_p left the run_id "
        "unchanged (even though tau_B differed by 16.1x). It would be mistaken for a "
        "completed run of a different system, and the old results reported as the new "
        "system's. Fixed by putting `system` (physics_only) into the hash "
        "(reproduction: `verify/verify_l3_spec_gaps.py`). "
        "The 7 existing runs are preserved as legacy (user decision 2026-08-04) -- "
        "their physics and results are valid; only their directory names are not "
        "reproducible under the new convention. "
        "The L5 migration (2026-08-05) did not change any params/numerics key, so "
        "this run_id is unchanged -- `30caa5c9e0` was re-run through `RUN.builder` "
        "and the observables matched to 15 decimal places (energy consistency "
        "105.510722899358, hexagonal NN distance 1.6170154574513436)."),
}

# ★ The L5 migration (2026-08-05, bdbot.run's @RUN.builder contract) changed the
#   shape of metrics.json. chain-bend and trap-drag were already in that shape, and
#   this script had never compared those two, so it went unnoticed -- it first
#   surfaced when trap-2d-5um and soft-r3 were migrated. The entries below are either
#   **a move rather than a loss of information** (RENAMED, where the value is still
#   compared) or a deliberate omission from this schema (checks[] accepts only the
#   design and post_run phases -- post_checks stays in result; prediction_source was
#   superseded by the richer role/scope/derivation/note).
#   The physics values themselves were separately confirmed to match the old run to
#   15 decimal places (see the conversation record) -- so "no violation" here does
#   NOT mean the schema is the same. It means **every deleted field is accounted
#   for**.
ALLOWED_DELETION_PREFIXES = (
    "physical.",             # now holds only the star parameters (the trap-drag
                              # convention). The SI values can be reconstructed from
                              # the spec at any time via spec.physical().
    "checks[5]", "checks[6]",  # post_run checks -> moved to result.post_checks[]
    "numerics.bias_predicted_pct", "numerics.dt_over_tau_B",
    "numerics.primary_sem_pct", "numerics.stat_target_pct",
    # soft-r3: observables with no prediction (psi6, U_per_N) now live only in
    # result.psi6 / result.pe_mean -- observables[] keeps just the two
    # 'measured vs predicted' entries (no information is lost).
    "observables[2]", "observables[3]",
)
ALLOWED_DELETION_SUFFIXES = (
    ".prediction_source",    # -> role/scope/derivation/note (a richer structure)
)

# Fields the refactor is allowed to change
TOLERATED = {
    "run_id": "judged separately in (1) via EXPECTED_RUN_ID (not double-judged here)",
    "wall_seconds": "wall clock -- differs on every run",
    "steps_per_second": "derived from the wall clock",
    "numerics.x2_sem_pct": "block SEM upcast to float64 (1.2e-6, the more accurate "
                           "of the two)",
    "numerics.primary_sem_pct": "same reason",
    "observables[0].err_pct": (
        "sign-convention unification -- **a bug fix**. In the 1-B original this one "
        "row alone used (predicted - measured)/|measured|, giving it the opposite "
        "sign from the other rows in the same file ((measured - reference)/reference) "
        "and from 1-A. The denominator also changed from |measured| to |predicted|, "
        "so the magnitude shifts by the relative difference between the two "
        "(8.7e-5). The |error| < 2% verdict is unchanged."),
    "observables[0].prediction_source": (
        "the 1-B original packed the whole explanatory sentence into this one field. "
        "The L5 schema separates a short tag (source) from the full sentence (note), "
        "role and derivation -- this run was produced by code from **before** that "
        "separation existed, so its source stayed 'none' (the default), while the "
        "current code passes source='consistency' (the full sentence is present "
        "verbatim in observables[0].note -- no information lost, the value did "
        "migrate)."),
    "observables[1].prediction_source": "same reason -- the current code passes "
                                        "source='lattice'.",
}
# Schema fields **deliberately added** in 1-C and the work before it. Deletions are
# still not permitted.
EXPECTED_NEW_FIELDS = {
    "schema": "the metrics schema version (bdbot.metrics/0.2)",
    "checks[].hard": "the hard/soft distinction -- 1-A passed everything, so the "
                     "distinction never showed",
    "checks[].phase": "design / post_run",
    "equilibration.source": "the case declares its equilibration indicator (1-A's "
                            "anchor displacement cannot be used for a diffusive "
                            "system)",
    "equilibration.series_key": "same reason",
    "equilibration.label": "same reason",
    "numerics.dt_over_tau_B": "dt against the reference time is recorded too",
    "numerics.primary_sem_pct": "the primary indicator differs per case "
                                "(x2_sem_pct is 1-A only)",
    "numerics.stat_target_pct": "the case declares its statistical target",
    "checks[].op": ("the comparison direction (<=/>=). Needed when re-reading an L3 "
                    "spec -- without it a '>=' check (the observation window, for "
                    "instance) is restored as '<=' and the verdict inverts "
                    "(bdbot.nondim.load)."),
}

# Entries where only the sign convention changed: if the magnitude moves by more than
# the relative difference between the two values, treat it as a real regression
SIGN_ONLY = {"observables[0].err_pct": 1e-3}
RTOL = 1e-9

# ★ Fields whose **name only** changed (old key -> new key). The value is still
#   compared -- skipping the check because the name changed would be "silently
#   passing".
#   The L3 work changed the keys of `metrics.dimensionless` from report display
#   strings to symbols. metrics.json is postmortem's only input, and a key like
#   `'k*     = k d²/kT   트랩 vs 열요동'` could not be queried.
#
# ⛔ **THE KEYS ON THE LEFT ARE DELIBERATELY KOREAN AND MUST STAY THAT WAY.**
#   They are not prose: they are the literal keys inside archived metrics.json
#   files. Measured 2026-08-29: 9 archived runs carry 81 such Korean keys under
#   `dimensionless`, and runs/ is a frozen, non-reproducible archive.
#   Translating a key here would stop it matching, so the field would be reported as
#   DELETED -- a false violation. Same load-bearing-Korean rule as
#   bdbot/interactions.py and bdbot/runid.py:84.
RENAMED = {
    "dimensionless.k*     = k d²/kT   트랩 vs 열요동": "dimensionless.k*",
    "dimensionless.l_k/d  = 1/√k*     요동폭 vs 입자": "dimensionless.l_k/d",
    "dimensionless.tau_k/tau_B = 1/k*": "dimensionless.tau_k/tau_B",
    "dimensionless.dt/tau_k           적분 해상": "dimensionless.dt/tau_k",
    "dimensionless.T_obs/tau_k        관측창": "dimensionless.T_obs/tau_k",
    "dimensionless.Gamma  = U(a_mean)/kT   결합 vs 열요동 ★": "dimensionless.Gamma",
    "dimensionless.A      = U(d)/kT        접촉 결합": "dimensionless.A",
    "dimensionless.phi                     밀집도": "dimensionless.phi",
    "dimensionless.a_mean/d                평균간격": "dimensionless.a_mean/d",
    "dimensionless.L/d                     박스 크기": "dimensionless.L/d",
    "dimensionless.r_c/d                   컷오프": "dimensionless.r_c/d",
    "dimensionless.r_c/a_mean              컷오프(이웃 껍질 수)": "dimensionless.r_c/a_mean",
    "dimensionless.dt/tau_int              적분 해상": "dimensionless.dt/tau_int",
    "dimensionless.T_obs/tau_B             관측창": "dimensionless.T_obs/tau_B",
    "dimensionless.St     = tau_p/tau_B    관성 vs 확산": "dimensionless.St",
    "numerics.x2_sem_pct": "result.x2_sem_pct",
    # L5 migration (2026-08-05) -- per-case extra information moved from the
    # top-level `structure` to the `result` that RUN.execute() fills. The value is
    # still compared as-is.
    "structure.psi6": "result.psi6",
    "structure.psi6_sem": "result.psi6_sem",
    "structure.nn_distance_d": "result.nn_distance_d",
    "structure.nn_std_rel": "result.nn_std_rel",
    "structure.min_sep_d": "result.min_sep_d",
    "structure.Gamma": "result.Gamma",
    "structure.coord_hist": "result.coord_hist",
    "structure.state_predicted": "result.state_predicted",
    "structure.u_rms_rel_einstein": "result.u_rms_rel_einstein",
}


def flatten(d, pre=""):
    out = {}
    for k, v in d.items():
        p = f"{pre}{k}"
        if isinstance(v, dict):
            out.update(flatten(v, p + "."))
        elif isinstance(v, list) and v and isinstance(v[0], dict):
            for i, x in enumerate(v):
                out.update(flatten(x, f"{p}[{i}]."))
        else:
            out[p] = v
    return out


def run_id_of(script, extra) -> str | None:
    out = subprocess.run([PY, str(ROOT / script)] + extra + ["--report"],
                         capture_output=True, text=True, cwd=ROOT).stdout
    m = re.search(r"run_id=(\S+)", out)
    return m.group(1) if m else None


def compare(label, ref_path: Path, new_path: Path) -> tuple[bool, list[str]]:
    ref, new = json.loads(ref_path.read_text()), json.loads(new_path.read_text())
    fr, fn = flatten(ref), flatten(new)
    lines, bad = [], []
    same = 0
    renamed = 0
    for k in sorted(fr):
        key = k
        if key not in fn and key in RENAMED and RENAMED[key] in fn:
            key = RENAMED[k]          # name-only change -- the value is compared below
            renamed += 1
        if key not in fn:
            if k.startswith(ALLOWED_DELETION_PREFIXES) or k.endswith(ALLOWED_DELETION_SUFFIXES):
                lines.append(f"    ~ deleted (permitted): {k}")
            else:
                bad.append(f"    ✗ deleted field: {k}")
            continue
        a, b = fr[k], fn[key]
        if isinstance(a, (int, float)) and not isinstance(a, bool) \
                and isinstance(b, (int, float)) and not isinstance(b, bool):
            if a == b or abs(a - b) <= RTOL * max(abs(a), abs(b), 1e-300):
                same += 1
            elif k in TOLERATED:
                tol = SIGN_ONLY.get(k)
                if tol is not None and abs(abs(a) - abs(b)) > tol * max(abs(a), abs(b)):
                    bad.append(f"    ✗ {k}: not just a sign change -- the magnitude "
                               f"differs too "
                               f"({a:.12g} → {b:.12g})")
                else:
                    lines.append(f"    ~ {k}: {a:.12g} → {b:.12g}   ({TOLERATED[k]})")
            else:
                rel = 100 * (b - a) / abs(a) if a else float("nan")
                bad.append(f"    ✗ {k}: {a!r} → {b!r}  ({rel:+.4g}%)")
        elif a == b:
            same += 1
        elif k in TOLERATED:
            lines.append(f"    ~ {k}  ({TOLERATED[k]})")
        else:
            bad.append(f"    ✗ {k}:\n        before: {str(a)[:80]}\n"
                       f"        after:  {str(b)[:80]}")
    added = sorted(set(fn) - set(fr) - set(RENAMED.values()))
    print(f"  identical {same} . renamed only {renamed} . permitted changes "
          f"{len(lines)} . added {len(added)} . violations {len(bad)}")
    for ln in lines:
        print(ln)
    for ln in bad:
        print(ln)
    if added:
        groups = sorted({re.sub(r"\[\d+\]", "[]", a) for a in added})
        print(f"    + added fields ({len(added)}): " + ", ".join(groups))
    return not bad, bad


def main() -> int:
    print("=" * 84)
    print("Phase 1-C equivalence check -- pre-refactor (1-B snapshot) vs current code")
    print("=" * 84)
    ok_all = True
    for label, script, extra, ref_name in TARGETS:
        print(f"\n■ {label}   ({script} {' '.join(extra)})")
        rid = run_id_of(script, extra)
        ref_path = REF / ref_name
        if not ref_path.exists():
            # ⛔ Say what is wrong. This used to raise a bare FileNotFoundError, which
            #    reads like a bug in the script rather than an absent baseline.
            print(f"  ⛔ BASELINE MISSING: {ref_path.relative_to(ROOT)}")
            print(f"     This snapshot was never committed, so this regression test "
                  f"has no baseline")
            print(f"     and cannot run -- for anyone, on any clone.")
            print(f"     ⚠️ Do NOT regenerate it from the current code. A snapshot "
                  f"taken now compares")
            print(f"        current against current, passes unconditionally, and "
                  f"guards nothing -- turning")
            print(f"        a visibly broken check into a silently passing one. "
                  f"Recover the pre-1-C file,")
            print(f"        or retire this script. Do not manufacture agreement.")
            ok_all = False
            continue
        ref = json.loads(ref_path.read_text())
        exp = EXPECTED_RUN_ID.get(ref["run_id"])
        if rid == ref["run_id"]:
            id_ok, note = True, "✓ identical"
        elif exp and rid == exp[0]:
            id_ok, note = True, "✓ an intended change"
        else:
            id_ok, note = False, "✗ differs (no reason recorded)"
        ok_all &= id_ok
        print(f"  ① run_id  {rid}")
        print(f"     snapshot {ref['run_id']}   {note}")
        if exp and rid == exp[0]:
            print(f"     reason: {exp[1]}")
        new_path = ROOT / "runs" / (rid or "") / "metrics.json"
        if not new_path.exists():
            print(f"  (2) no re-run result -- run it first:  "
                  f"$PY {script} {' '.join(extra)} --force")
            ok_all = False
            continue
        print("  (2)(3) metrics comparison")
        ok, _ = compare(label, ref_path, new_path)
        ok_all &= ok
    print()
    print("=" * 84)
    print("✓ PASS -- the refactor did not change any result" if ok_all
          else "✗ FAIL -- check the violations above")
    print("=" * 84)
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
