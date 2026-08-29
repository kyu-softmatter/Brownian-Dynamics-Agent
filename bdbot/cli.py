"""The `bdbot` CLI -- the single entry point to the front end (L0 -> L3).

It has to behave identically **outside** a Claude Code session: cron, a script, or
another person must get the same result from the same command, and these commands
are also what a hook would intercept.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python

    $PY -m bdbot.cli status                     pipeline progress for every case
    $PY -m bdbot.cli intake init  <folder>      an observation.yaml template
    $PY -m bdbot.cli intake check <folder>      L0 schema plus readiness verdict
    $PY -m bdbot.cli system check <folder>      L2 schema, tiers, recomputed derived values
    $PY -m bdbot.cli nondim report <case>       the L3 report (does not run anything)
    $PY -m bdbot.cli nondim spec   <case>       the L3 spec -> specs/<run_id>.json
    $PY -m bdbot.cli nondim show   <run_id>     redraw the report from the spec alone, and verify the hash
    $PY -m bdbot.cli nondim list                list specs/
    $PY -m bdbot.cli run <case> [-- ...]        run (delegated to the case script)

Exit codes: 0 ok . 1 FAIL (schema error) . 2 BLOCKED (unresolved gap) . 3 usage error
-> a script or a hook can read the verdict as a code.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from . import intake as _intake
from . import interactions as _inter
from . import physical as _physical

ROOT = Path(__file__).resolve().parent.parent
INTAKE_DIR = ROOT / "intake"

EXIT_OK, EXIT_FAIL, EXIT_BLOCKED, EXIT_USAGE = 0, 1, 2, 3

# case -> script. `L3_ONLY` means it goes as far as non-dimensionalization
# (--report/--spec) but has no L4 (execution) -- `status` shows "L3" instead of "O"
# to distinguish it from an end-to-end script. The point of that table is to show
# precisely **what is missing**, not what is present.
CASE_SCRIPTS = {
    "trap-2d-5um": "cases/trap_2d_5um.py",
    "soft-r3-2d-A-sweep": "cases/soft_r3_2d.py",
    "abp-rod-2d-run-flip": "cases/abp_rod_2d.py",
    "trap-drag-2d-hex300": "cases/trap_drag_2d.py",
    "chain-bend-2d-oscill": "cases/chain_bend_2d.py",
    "chain-bend-2d-dlvo": "cases/chain_bend_dlvo_2d.py",
    "chain-relax-2d-dlvo": "cases/chain_relax_2d_dlvo.py",
    # 2026-08-28 (merge): `cases/network_3d.py` existed and had produced runs, but
    # was missing from this registry alone, so `bdbot.cli run network` refused with
    # "no end-to-end script". The `-` in the status table's script column did not
    # mean the script was absent either.
    "network": "cases/network_3d.py",
}
L3_ONLY = frozenset({"trap-drag-2d-hex300", "chain-bend-2d-oscill"})


def _resolve(folder: str) -> Path:
    p = Path(folder)
    if not p.exists():
        alt = INTAKE_DIR / folder
        if alt.exists():
            return alt
        cand = [d.name for d in sorted(INTAKE_DIR.iterdir()) if d.is_dir()]
        print(f"path not found: {folder}\ncases: {', '.join(cand)}", file=sys.stderr)
        raise SystemExit(EXIT_USAGE)
    return p


def _cases() -> list[Path]:
    if not INTAKE_DIR.exists():
        return []
    return [d for d in sorted(INTAKE_DIR.iterdir())
            if d.is_dir() and (d / "observation.yaml").exists()]


# ══════════════════════════════════════════════════════════════════════
def cmd_intake_init(args) -> int:
    p = Path(args.folder) if "/" in args.folder else INTAKE_DIR / args.folder
    ok, msg = _intake.init_template(p, force=args.force)
    print(msg)
    if ok:
        print("\nnext: read the image and start by filling in the transcription "
              "(skill `bd-intake` step 1).")
        print(f"      then check:  $PY -m bdbot.cli intake check {p}")
    return EXIT_OK if ok else EXIT_USAGE


def cmd_intake_check(args) -> int:
    obs = _intake.load(_resolve(args.folder))
    print(_intake.render_check(obs))
    if obs.errors:
        return EXIT_FAIL
    unspec = [m for m in obs.open_missing
              if _inter.looks_like_interaction(str(m.get("symbol", "")))]
    if unspec:
        syms = ", ".join(str(m.get("symbol")) for m in unspec)
        print(f"\nthe interaction is unspecified ({syms}).")
        print(f"  -> to see the standard candidates:  "
              f"$PY -m bdbot.cli intake suggest {args.folder}")
    return EXIT_BLOCKED if obs.open_missing else EXIT_OK


def cmd_intake_suggest(args) -> int:
    """When the sketch has no interaction, recommend and **ask.**

    The decision is the human's.
    """
    obs = _intake.load(_resolve(args.folder))
    if obs.errors:
        print(_intake.render_check(obs))
        return EXIT_FAIL
    print(_inter.render_suggestion(obs))
    return EXIT_BLOCKED if obs.open_missing else EXIT_OK


def cmd_interactions_list(args) -> int:
    print(_inter.render_catalog())
    return EXIT_OK


def cmd_system_check(args) -> int:
    s = _physical.load(_resolve(args.folder))
    print(_physical.render_check(s))
    return EXIT_FAIL if s.errors else EXIT_OK


def cmd_status(args) -> int:
    """L0 -> L2 -> L3/run progress in one table. The screen for seeing what is
    blocked where.

    WARNING: the run count only counts runs that have a `result.txt`, and that file
    is written by the *case script*, not by `bdbot.run`. A case that never added
    that line reports 0 runs while its metrics.json files sit on disk.
    """
    rows = []
    for d in _cases():
        obs = _intake.load(d)
        l0 = "FAIL" if obs.errors else ("BLOCKED" if obs.open_missing else "READY")
        has_sys = (d / "system.yaml").exists()
        if has_sys:
            s = _physical.load(d)
            l2 = "FAIL" if s.errors else "READY"
        else:
            l2 = "—"
        script = CASE_SCRIPTS.get(d.name)
        mark = ("L3" if d.name in L3_ONLY else "O") if script else "—"
        runs = sorted(ROOT.glob(f"runs/{d.name}__*"))
        runs = [r for r in runs if (r / "result.txt").exists()]
        # L3 artefacts -- written by `nondim spec`. An earlier stage than a run, so
        # counted separately.
        specs = sorted((ROOT / "specs").glob(f"{d.name}__*.json")) \
            if (ROOT / "specs").exists() else []
        blockers = ", ".join(m.get("symbol", "?") for m in obs.open_missing) or "—"
        rows.append((d.name, l0, l2, mark, len(specs), len(runs), blockers))

    w = max([len(r[0]) for r in rows] + [8])
    print("=" * (w + 62))
    print("bdbot status -- from sketch to run")
    print("=" * (w + 62))
    print(f"{'case':<{w}}  {'L0':<8}{'L2':<8}{'script':<9}{'specs':>6}{'runs':>6}   blocking gaps")
    print("-" * (w + 62))
    for name, l0, l2, sc, ns, nr, blk in rows:
        print(f"{name:<{w}}  {l0:<8}{l2:<8}{sc:<9}{ns:>5}{nr:>5}   {blk[:30]}")
    print("-" * (w + 62))
    n_ready = sum(1 for r in rows if r[1] == "READY")
    n_spec = sum(1 for r in rows if r[4] > 0)
    n_run = sum(1 for r in rows if r[5] > 0)
    print(f"{len(rows)} cases . L0 READY {n_ready} . cases with an L3 spec {n_spec} . "
          f"cases with runs {n_run}")
    print("script: O = end-to-end (through L4) . L3 = non-dimensionalization only (--report/--spec)")
    if any(r[1] == "FAIL" for r in rows):
        print("\nFAIL cases: run `intake check <case>` to see the schema errors.")
    blocked = [r[0] for r in rows if r[1] == "BLOCKED"]
    if blocked:
        print(f"\nBLOCKED: {len(blocked)} -- stopped rather than inventing a value the "
              f"sketch does not contain.")
        print("  a human must supply it, or it must be found in the KB. "
              "`intake check <case>` has the list.")
    print("=" * (w + 58))
    return EXIT_OK


def _dispatch(case: str, extra: list[str]) -> int:
    script = CASE_SCRIPTS.get(case)
    if not script:
        print(f"'{case}' has no end-to-end script yet.\n"
              f"available: {', '.join(CASE_SCRIPTS)}\n"
              f"(for the rest, the physical system is not settled -- see `bdbot.cli status`)",
              file=sys.stderr)
        return EXIT_USAGE
    cmd = [sys.executable, str(ROOT / script), *extra]
    print(f"→ {' '.join(cmd[1:])}\n", flush=True)
    return subprocess.call(cmd, cwd=ROOT)


def cmd_nondim_report(args) -> int:
    """The L3 report. **Delegated to the case script** -- the scale ledger differs
    per system.

    The shared parts (ledger structure, check categories, the renderer) live in
    `bdbot.scales`/`checks`/`report`; which scales go in is decided by the case
    (skill `bd-physics` section 6.3).
    """
    return _dispatch(args.case, ["--report", *args.extra])


def cmd_nondim_spec(args) -> int:
    """Write the L3 spec to `specs/<run_id>.json` (does not run) -- the object of
    the third human check.
    """
    return _dispatch(args.case, ["--spec", *args.extra])


def cmd_nondim_show(args) -> int:
    """Redraw the report from **the stored spec alone** -- with no case script.

    This is where a human checks the path L4 will use. If everything shows here the
    spec is self-sufficient; if not, something is missing from it. A hand-edited
    spec is caught by the hash verification.
    """
    from . import nondim as _nd
    p = Path(args.spec)
    if not p.exists():
        alt = ROOT / "specs" / args.spec
        alt = alt if alt.exists() else ROOT / "specs" / f"{args.spec}.json"
        if not alt.exists():
            have = sorted(x.name for x in (ROOT / "specs").glob("*.json")) \
                if (ROOT / "specs").exists() else []
            print(f"spec not found: {args.spec}\n"
                  + ("available:\n  " + "\n  ".join(have) if have else
                     "specs/ is empty -- create one first with `nondim spec <case>`."),
                  file=sys.stderr)
            return EXIT_USAGE
        p = alt
    try:
        spec = _nd.load(p)
    except Exception as e:
        print(f"cannot read the spec: {e}", file=sys.stderr)
        return EXIT_FAIL
    print(spec.render())
    ok, _ = spec.verify_hash()
    n_err = sum(1 for i in spec.raw.get("l3_issues", []) if i["level"] == "error")
    if not ok:
        print("\nthis spec was hand-edited -- the run_id does not match its contents (rule 2).")
        return EXIT_FAIL
    return EXIT_FAIL if (n_err or spec.verdict.startswith("FAIL")) else EXIT_OK


def cmd_nondim_list(args) -> int:
    from . import nondim as _nd
    d = ROOT / "specs"
    files = sorted(d.glob("*.json")) if d.exists() else []
    if not files:
        print("specs/ is empty -- create one with `nondim spec <case>`.")
        return EXIT_OK
    print(f"{'run_id':<44}{'verdict':<22}hash")
    print("-" * 78)
    for f in files:
        try:
            s = _nd.load(f)
            ok, _ = s.verify_hash()
            print(f"{s.run_id:<44}{s.verdict:<22}{'ok' if ok else 'HAND-EDITED'}")
        except Exception as e:
            print(f"{f.name:<44}{'read failed':<22}{e}")
    return EXIT_OK


def cmd_health(args) -> int:
    """L4 -- the numerical-health verdict (`bdbot.health`).

    Not a physics verifier. Divergence, NaN/Inf, stalling, collapse, and the
    comparison against the L3 ledger -- nothing else.
    """
    import subprocess
    cmd = [sys.executable, str(ROOT / "tools" / "health.py")]
    if args.all:
        cmd.append("--all")
    elif args.gate:
        cmd += ["--gate", args.target]
    elif args.target:
        cmd.append(args.target)
    else:
        cmd.append("--all")
    return subprocess.call(cmd)


def cmd_run(args) -> int:
    return _dispatch(args.case, list(args.extra))


# ══════════════════════════════════════════════════════════════════════
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="bdbot", description="the Brownian dynamics pipeline -- front-end entry point",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="exit codes: 0 ok . 1 FAIL . 2 BLOCKED . 3 usage error")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("status", help="progress for every case").set_defaults(fn=cmd_status)

    p_in = sub.add_parser("intake", help="L0 intake")
    s_in = p_in.add_subparsers(dest="sub", required=True)
    q = s_in.add_parser("init", help="create an observation.yaml template")
    q.add_argument("folder")
    q.add_argument("--force", action="store_true")
    q.set_defaults(fn=cmd_intake_init)
    q = s_in.add_parser("check", help="schema plus readiness verdict")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_intake_check)
    q = s_in.add_parser("suggest", help="recommendations and questions for an unspecified interaction")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_intake_suggest)
    q = s_in.add_parser("list", help="alias for status")
    q.set_defaults(fn=cmd_status)

    q = sub.add_parser("interactions", help="the colloidal interaction catalogue")
    q.set_defaults(fn=cmd_interactions_list)

    p_sy = sub.add_parser("system", help="the L2 physical system")
    s_sy = p_sy.add_subparsers(dest="sub", required=True)
    q = s_sy.add_parser("check", help="schema, tiers, recomputed derived values")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_system_check)

    p_nd = sub.add_parser("nondim", help="L3 non-dimensionalization")
    s_nd = p_nd.add_subparsers(dest="sub", required=True)
    q = s_nd.add_parser("report", help="the report (does not run anything)")
    q.add_argument("case")
    q.add_argument("extra", nargs=argparse.REMAINDER)
    q.set_defaults(fn=cmd_nondim_report)
    q = s_nd.add_parser("spec", help="write the L3 spec to specs/<run_id>.json")
    q.add_argument("case")
    q.add_argument("extra", nargs=argparse.REMAINDER)
    q.set_defaults(fn=cmd_nondim_spec)
    q = s_nd.add_parser("show", help="redraw the report from the stored spec alone (no case code)")
    q.add_argument("spec", help="a filename under specs/, or a run_id")
    q.set_defaults(fn=cmd_nondim_show)
    q = s_nd.add_parser("list", help="list specs/ and verify hashes")
    q.set_defaults(fn=cmd_nondim_list)

    q = sub.add_parser("health", help="the L4 numerical-health verdict (not physics verification)")
    q.add_argument("target", nargs="?", help="runs/<run_id>, or a spec path for --gate")
    q.add_argument("--all", action="store_true", help="every completed run")
    q.add_argument("--gate", action="store_true", help="the pre-run gate only (a spec path)")
    q.set_defaults(fn=cmd_health)

    q = sub.add_parser("run", help="run a case (delegated to its script)")
    q.add_argument("case")
    q.add_argument("extra", nargs=argparse.REMAINDER)
    q.set_defaults(fn=cmd_run)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
