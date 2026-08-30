"""Re-run the `trap-drag` ensemble **with the current code**, including the
step-resolution measurement.

## Why re-run at all

The 79 legacy runs do not carry `step_drift_max_sigma` (the input to the L4->L3
feedback), and it **cannot be measured retroactively** -- `run_id` is a content
address, so re-running produces a new run under a different id, and the GSD-replay
workaround is invalid under time-dependent driving (a moving trap), overestimating
by 16x (`verify/probe_gsd_replay.py`). So a fresh ensemble is the only way to fill
the coverage.

## Design (reconstructed from the existing 63 runs)

    traverse = 0.117647           one box traverse (T_obs multiplier). Same as the
                                  existing ensemble
    v        = 0.05 0.1 0.5 1.5 4 12 32   µm/s  (the 7 points in L2
                                  `external.drag_velocity`)
    seed     = default (20260804) + 1..8   -> 9 runs per velocity, 63 in total

The existing 63 runs totalled **72.5 hours** of wall-clock. Hence:

  . **longest first** (v=0.05 is 1.9h per run, v=0.5 is 0.3h). Running the short ones
    first leaves only long ones at the end, so cores idle and the makespan grows
    (LPT scheduling).
  . Concurrency defaults to 6. This machine has 4 performance + 6 efficiency cores,
    so going above 6 makes runs steal performance cores from each other. HOOMD is
    single-threaded per run (no MPI).
  . **Restartable**: completed runs are skipped. Interrupting it is safe -- it
    resumes. (See `measured()` below for what "completed" actually means here; it is
    NOT result.txt.)

    $PY scratch/run_trap_drag_ensemble.py --jobs 6
    $PY verify/run_trap_drag_ensemble.py --dry-run     # plan only
    $PY verify/run_trap_drag_ensemble.py --wait-for chain_bend_2d
                                                      # wait for that to finish first
"""
from __future__ import annotations

import argparse
import json
import re
import subprocess
import sys
import time
from concurrent.futures import FIRST_COMPLETED, ThreadPoolExecutor, wait
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
CASE = ROOT / "cases" / "trap_drag_2d.py"
LOG = ROOT / "verify" / "trap_drag_ensemble_run.log"

TRAVERSE = 0.117647
VELOCITIES = [0.05, 0.1, 0.5, 1.5, 4.0, 12.0, 32.0]
SEEDS = [None, 1, 2, 3, 4, 5, 6, 7, 8]        # None = the case's default seed
                                              # (20260804)

# Wall-clock per run [h] measured from the existing ensemble -- used only to order
# longest-first, so accuracy does not matter
COST_H = {0.05: 1.86, 0.1: 1.33, 1.5: 1.17, 4.0: 1.15, 12.0: 1.14, 32.0: 1.13, 0.5: 0.29}


def jobs_plan() -> list[dict]:
    out = []
    for v in sorted(VELOCITIES, key=lambda x: -COST_H.get(x, 1.0)):    # longest first
        for s in SEEDS:
            args = ["--traverse", repr(TRAVERSE), "--v", f"{v:g}"]
            if s is not None:
                args += ["--seed", str(s)]
            out.append({"v": v, "seed": s, "args": args, "cost_h": COST_H.get(v, 1.0)})
    return out


def spec_of(job: dict) -> tuple[str, Path] | tuple[None, None]:
    """Pin down run_id up front via `--spec` (without running).

    Used to decide whether a run is already complete.
    """
    p = subprocess.run([PY, str(CASE), *job["args"], "--spec"],
                       capture_output=True, text=True, cwd=ROOT)
    m = re.search(r"(specs/\S+\.json)", p.stdout)
    if not m:
        return None, None
    sp = ROOT / m.group(1)
    return sp.stem, sp


def already_done(run_id: str) -> bool:
    """Completion test. ⚠️ **Do NOT test this with `result.txt`.**

    `trap-drag` goes through the `RUN.execute` path, which **does not write**
    `result.txt` (it leaves only l4.json, metrics.json, observables.npz, spec.json
    and traj_A.gsd). Testing on result.txt broke two things:
      (1) this function was always False -> it reported "0 runs already complete",
          and because the **run_id matched a legacy run**,
          `runid.prepare_outdir` kept only `record.json` and deleted and overwrote
          the rest (33 legacy runs lost -- 2026-08-06).
      (2) restarts stopped being idempotent, so finished runs were re-run from
          scratch.
    The real completion signal for this case is `metrics.json`, and since the point
    of this ensemble is whether the **measurement** happened too, both are checked.
    """
    return (ROOT / "runs" / run_id / "metrics.json").exists() and measured(run_id)


def measured(run_id: str) -> bool:
    f = ROOT / "runs" / run_id / "metrics.json"
    if not f.exists():
        return False
    try:
        return "step_drift_max_sigma" in json.loads(f.read_text()).get("numerics", {})
    except Exception:
        return False


def log(msg: str) -> None:
    line = f"[{time.strftime('%H:%M:%S')}] {msg}"
    print(line, flush=True)
    with LOG.open("a") as fh:
        fh.write(line + "\n")


def run_one(job: dict) -> dict:
    t0 = time.time()
    label = f"v={job['v']:g}" + (f" s={job['seed']}" if job["seed"] is not None
                                 else " s=default")
    p = subprocess.run([PY, str(CASE), *job["args"]],
                       capture_output=True, text=True, cwd=ROOT)
    dt = time.time() - t0
    rid = job.get("run_id") or "?"
    ok = p.returncode == 0 and measured(rid)
    tail = "\n".join(p.stdout.strip().splitlines()[-4:]) if p.returncode else ""
    if p.returncode != 0:
        tail = "\n".join((p.stdout + p.stderr).strip().splitlines()[-6:])
    log(f"{'✓' if ok else '✗'} {label:<16} {dt/3600:5.2f}h  {rid}"
        + ("" if ok else f"\n    rc={p.returncode} measured={measured(rid)}\n"
                         f"    {tail}"))
    return {**job, "ok": ok, "wall_h": dt / 3600, "rc": p.returncode}


def wait_for(pattern: str, poll: int = 60) -> None:
    """Wait until processes matching the pattern are gone.

    Avoids fighting another session's work for cores.

    ★ **Must not count itself.** This script's own command line contains
      `--wait-for <pattern>`, so `pgrep -f <pattern>` matches its own process ->
      it would wait forever. (The same trap was hit once with a monitor in this
      session.) Its own PID is excluded.
    """
    import os
    me = {os.getpid(), os.getppid()}
    first = True
    while True:
        r = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True)
        pids = {int(x) for x in r.stdout.split() if x.strip().isdigit()} - me
        if not pids:
            log(f"'{pattern}' confirmed finished -- starting the ensemble")
            return
        if first:
            log(f"'{pattern}': {len(pids)} still running -- waiting for them "
                f"(every {poll}s) "
                f"pids={sorted(pids)[:8]}")
            first = False
        time.sleep(poll)


def procs_matching(pattern: str) -> set[int]:
    """PIDs of **other** processes matching the pattern.

    Excludes this process, whose own command line contains the pattern.
    """
    import os
    r = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True)
    return ({int(x) for x in r.stdout.split() if x.strip().isdigit()}
            - {os.getpid(), os.getppid()})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6, help="concurrency (4 performance + 6 efficiency cores)")
    ap.add_argument("--jobs-max", type=int, default=None,
                    help="raise concurrency to this once --ramp-when-clear is clear")
    ap.add_argument("--ramp-when-clear", default=None,
                    help="raise concurrency to --jobs-max once processes matching "
                         "this pattern are gone (starts immediately rather than "
                         "waiting -- shares cores first, then takes all of them)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--wait-for", default=None, help="wait for processes matching this pattern to finish first")
    a = ap.parse_args()

    plan = jobs_plan()
    log("=" * 76)
    log(f"trap-drag ensemble -- {len(plan)} runs (traverse={TRAVERSE}, "
        f"{len(VELOCITIES)} velocities x {len(SEEDS)} seeds)")
    log(f"expected total {sum(j['cost_h'] for j in plan):.1f}h based on the existing "
        f"ensemble . concurrency {a.jobs} -> makespan roughly "
        f"{sum(j['cost_h'] for j in plan)/a.jobs:.1f}h")
    log("=" * 76)

    # Pin down run_id first (generating a spec is cheap, it is not a run) -> so
    # already-finished runs can be skipped
    todo, skip = [], []
    for j in plan:
        rid, sp = spec_of(j)
        if rid is None:
            log(f"✗ spec generation failed -- v={j['v']} seed={j['seed']} (skipped)")
            continue
        j["run_id"] = rid
        (skip if already_done(rid) else todo).append(j)
    log(f"{len(todo)} run(s) to do . {len(skip)} already complete")
    for j in skip:
        log(f"  . skipped {j['run_id']}  "
            f"(measured={'Y' if measured(j['run_id']) else 'N'})")

    if a.dry_run:
        log("--dry-run -- stopping here")
        for j in todo:
            log(f"  -> {j['run_id']}  expected {j['cost_h']:.2f}h")
        return 0
    if not todo:
        log("nothing to do")
        return 0

    if a.wait_for:
        wait_for(a.wait_for)

    t0 = time.time()
    done, failed = [], []
    # ★ Scheduling is done by hand so concurrency can be **raised while running**.
    #   ThreadPoolExecutor cannot change max_workers after construction. This starts
    #   now, sharing cores with another session's work (chain-bend), and raises
    #   concurrency for the remaining queue once that finishes -- restarting later to
    #   raise it would **throw away runs in progress** (up to 1.9h each).
    limit = a.jobs
    ramped = a.ramp_when_clear is None or a.jobs_max is None
    queue, running = list(todo), {}
    pool = ThreadPoolExecutor(max_workers=max(a.jobs, a.jobs_max or a.jobs))
    n_done = 0
    while queue or running:
        if not ramped and not procs_matching(a.ramp_when_clear):
            limit = a.jobs_max
            ramped = True
            log(f"'{a.ramp_when_clear}' confirmed finished -> raising concurrency "
                f"{a.jobs} -> {limit}")
        while queue and len(running) < limit:
            j = queue.pop(0)
            running[pool.submit(run_one, j)] = j
        if not running:
            continue
        fin, _ = wait(list(running), return_when=FIRST_COMPLETED, timeout=30)
        for f in fin:
            r = f.result()
            del running[f]
            n_done += 1
            (done if r["ok"] else failed).append(r)
            el = (time.time() - t0) / 3600
            log(f"    progress {n_done}/{len(todo)} . elapsed {el:.2f}h . "
                f"ok {len(done)} failed {len(failed)} . concurrency {limit} . "
                f"queued {len(queue)}")
    pool.shutdown()

    log("=" * 76)
    log(f"done -- ok {len(done)} . failed {len(failed)} . "
        f"total {(time.time()-t0)/3600:.2f}h")
    for r in failed:
        log(f"  ✗ v={r['v']:g} seed={r['seed']} rc={r['rc']}")
    log("next: $PY -m bdbot.cli health --all      (check step coverage)")
    log("      $PY verify/trap_drag_ensemble.py   (re-analyse the ensemble)")
    log("=" * 76)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
