"""`trap-drag` 앙상블을 **현재 코드로** 새로 돌린다 — 스텝 해상 측정을 포함해서.

## 왜 새로 돌리나

레거시 런 79개는 `step_drift_max_sigma`(L4→L3 되먹임의 입력)를 갖고 있지 않고,
**소급 측정이 불가능**합니다 — `run_id` 가 콘텐츠 주소라 재실행하면 다른 id 의 새 런이
생기고, GSD 재생 우회는 시간 의존 구동(이동 트랩)에서 16배 과대로 무효입니다
(`scratch/probe_gsd_replay.py`). 그래서 커버리지를 채우는 유일한 길은 새 앙상블입니다.

## 설계 (기존 63런에서 복원)

    traverse = 0.117647           박스 횡단 1회 (T_obs 배율). 기존 앙상블과 동일
    v        = 0.05 0.1 0.5 1.5 4 12 32   µm/s  (L2 `external.drag_velocity` 7점)
    seed     = 기본(20260804) + 1…8       → 속도당 9런, 전체 63런

기존 63런의 벽시계 합계는 **72.5시간**이었습니다. 그래서:

  · **긴 것부터** 넣습니다 (v=0.05 가 런당 1.9h, v=0.5 는 0.3h). 짧은 것을 먼저 돌리면
    마지막에 긴 것만 남아 코어가 놀고 makespan 이 늘어납니다 (LPT 스케줄링).
  · 병렬도는 기본 6. 이 기계는 성능코어 4 + 효율코어 6 이라 6을 넘기면 성능코어를
    서로 빼앗습니다. HOOMD 는 런당 단일 스레드입니다 (MPI 없음).
  · **재시작 가능**: 완료된 런(`result.txt` 존재)은 건너뜁니다. 중단해도 이어서 돌립니다.

    $PY scratch/run_trap_drag_ensemble.py --jobs 6
    $PY scratch/run_trap_drag_ensemble.py --dry-run     # 계획만
    $PY scratch/run_trap_drag_ensemble.py --wait-for chain_bend_2d   # 끝나길 기다린 뒤 시작
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
SEEDS = [None, 1, 2, 3, 4, 5, 6, 7, 8]        # None = 케이스 기본 시드(20260804)

# 기존 앙상블에서 측정한 런당 벽시계 [h] — 긴 것부터 넣기 위한 순서용 (정확도 불필요)
COST_H = {0.05: 1.86, 0.1: 1.33, 1.5: 1.17, 4.0: 1.15, 12.0: 1.14, 32.0: 1.13, 0.5: 0.29}


def jobs_plan() -> list[dict]:
    out = []
    for v in sorted(VELOCITIES, key=lambda x: -COST_H.get(x, 1.0)):    # 긴 것부터
        for s in SEEDS:
            args = ["--traverse", repr(TRAVERSE), "--v", f"{v:g}"]
            if s is not None:
                args += ["--seed", str(s)]
            out.append({"v": v, "seed": s, "args": args, "cost_h": COST_H.get(v, 1.0)})
    return out


def spec_of(job: dict) -> tuple[str, Path] | tuple[None, None]:
    """`--spec` 로 run_id 를 먼저 확정한다 (실행하지 않음). 완료 여부 판정에 씀."""
    p = subprocess.run([PY, str(CASE), *job["args"], "--spec"],
                       capture_output=True, text=True, cwd=ROOT)
    m = re.search(r"(specs/\S+\.json)", p.stdout)
    if not m:
        return None, None
    sp = ROOT / m.group(1)
    return sp.stem, sp


def already_done(run_id: str) -> bool:
    """완료 판정. ⚠️ **`result.txt` 로 판정하면 안 됩니다.**

    `trap-drag` 는 `RUN.execute` 경로라 `result.txt` 를 **쓰지 않습니다**
    (l4.json·metrics.json·observables.npz·spec.json·traj_A.gsd 만 남깁니다).
    result.txt 로 판정했다가 두 가지가 깨졌습니다:
      ① 이 함수가 항상 False → "이미 완료 0런" 으로 보고했고, 실제로는 **레거시 런과
         run_id 가 같아서** `runid.prepare_outdir` 가 `record.json` 만 남기고 나머지를
         지우고 덮어썼습니다 (레거시 33런 손실 — 2026-08-06).
      ② 재시작이 멱등하지 않아 끝난 런을 처음부터 다시 돌립니다.
    이 케이스의 진짜 완료 신호는 `metrics.json` 이고, **측정까지** 됐는지가 이 앙상블의
    목적이므로 둘을 같이 봅니다.
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
    label = f"v={job['v']:g}" + (f" s={job['seed']}" if job["seed"] is not None else " s=기본")
    p = subprocess.run([PY, str(CASE), *job["args"]],
                       capture_output=True, text=True, cwd=ROOT)
    dt = time.time() - t0
    rid = job.get("run_id") or "?"
    ok = p.returncode == 0 and measured(rid)
    tail = "\n".join(p.stdout.strip().splitlines()[-4:]) if p.returncode else ""
    if p.returncode != 0:
        tail = "\n".join((p.stdout + p.stderr).strip().splitlines()[-6:])
    log(f"{'✓' if ok else '✗'} {label:<16} {dt/3600:5.2f}h  {rid}"
        + ("" if ok else f"\n    rc={p.returncode} 측정={measured(rid)}\n    {tail}"))
    return {**job, "ok": ok, "wall_h": dt / 3600, "rc": p.returncode}


def wait_for(pattern: str, poll: int = 60) -> None:
    """해당 패턴의 프로세스가 사라질 때까지 기다린다 (다른 세션 작업과 코어 다툼 방지).

    ★ **자기 자신을 세면 안 됩니다.** 이 스크립트의 명령줄에 `--wait-for <pattern>` 이
      들어 있어서 `pgrep -f <pattern>` 이 자기 프로세스를 잡습니다 → 영원히 기다립니다.
      (같은 함정을 이 세션에서 monitor 로 한 번 겪었습니다.) 자기 PID 를 뺍니다.
    """
    import os
    me = {os.getpid(), os.getppid()}
    first = True
    while True:
        r = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True)
        pids = {int(x) for x in r.stdout.split() if x.strip().isdigit()} - me
        if not pids:
            log(f"'{pattern}' 종료 확인 — 앙상블 시작")
            return
        if first:
            log(f"'{pattern}' {len(pids)}개 실행 중 — 끝나길 기다립니다 ({poll}s 간격) "
                f"pids={sorted(pids)[:8]}")
            first = False
        time.sleep(poll)


def procs_matching(pattern: str) -> set[int]:
    """패턴에 맞는 **남의** 프로세스 PID. 자기 자신은 뺀다 (명령줄에 패턴이 들어 있다)."""
    import os
    r = subprocess.run(["pgrep", "-f", pattern], capture_output=True, text=True)
    return ({int(x) for x in r.stdout.split() if x.strip().isdigit()}
            - {os.getpid(), os.getppid()})


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--jobs", type=int, default=6, help="동시 실행 수 (성능코어 4 + 효율코어 6)")
    ap.add_argument("--jobs-max", type=int, default=None,
                    help="--ramp-when-clear 가 비면 이 값까지 올린다")
    ap.add_argument("--ramp-when-clear", default=None,
                    help="이 패턴의 프로세스가 사라지면 병렬도를 --jobs-max 로 올린다 "
                         "(기다리지 않고 바로 시작한다 — 코어를 나눠 쓰다가 이후 전부 쓴다)")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--wait-for", default=None, help="이 패턴의 프로세스가 끝나길 먼저 기다림")
    a = ap.parse_args()

    plan = jobs_plan()
    log("=" * 76)
    log(f"trap-drag 앙상블 — {len(plan)}런 (traverse={TRAVERSE}, 속도 {len(VELOCITIES)} × "
        f"시드 {len(SEEDS)})")
    log(f"기존 앙상블 기준 예상 합계 {sum(j['cost_h'] for j in plan):.1f}h · "
        f"병렬 {a.jobs} → makespan 대략 {sum(j['cost_h'] for j in plan)/a.jobs:.1f}h")
    log("=" * 76)

    # run_id 를 먼저 확정 (스펙 생성은 실행이 아니라 싸다) → 이미 끝난 것 건너뛰기
    todo, skip = [], []
    for j in plan:
        rid, sp = spec_of(j)
        if rid is None:
            log(f"✗ 스펙 생성 실패 — v={j['v']} seed={j['seed']} (건너뜀)")
            continue
        j["run_id"] = rid
        (skip if already_done(rid) else todo).append(j)
    log(f"실행 대상 {len(todo)}런 · 이미 완료 {len(skip)}런")
    for j in skip:
        log(f"  · 건너뜀 {j['run_id']}  (측정={'O' if measured(j['run_id']) else 'X'})")

    if a.dry_run:
        log("--dry-run — 여기서 멈춥니다")
        for j in todo:
            log(f"  → {j['run_id']}  예상 {j['cost_h']:.2f}h")
        return 0
    if not todo:
        log("할 일이 없습니다")
        return 0

    if a.wait_for:
        wait_for(a.wait_for)

    t0 = time.time()
    done, failed = [], []
    # ★ 병렬도를 **실행 중에 올릴 수 있게** 직접 스케줄합니다. ThreadPoolExecutor 는
    #   max_workers 를 나중에 못 바꿉니다. 다른 세션 작업(chain-bend)과 코어를 나눠 쓰며
    #   지금 시작하고, 그쪽이 끝나면 남은 큐를 전부 쓰는 병렬도로 올립니다 —
    #   나중에 재시작해서 올리면 **진행 중인 런이 버려집니다**(런당 최대 1.9h).
    limit = a.jobs
    ramped = a.ramp_when_clear is None or a.jobs_max is None
    queue, running = list(todo), {}
    pool = ThreadPoolExecutor(max_workers=max(a.jobs, a.jobs_max or a.jobs))
    n_done = 0
    while queue or running:
        if not ramped and not procs_matching(a.ramp_when_clear):
            limit = a.jobs_max
            ramped = True
            log(f"'{a.ramp_when_clear}' 종료 확인 → 병렬도 {a.jobs} → {limit} 로 올립니다")
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
            log(f"    진행 {n_done}/{len(todo)} · 경과 {el:.2f}h · 성공 {len(done)} "
                f"실패 {len(failed)} · 병렬 {limit} · 대기 {len(queue)}")
    pool.shutdown()

    log("=" * 76)
    log(f"완료 — 성공 {len(done)} · 실패 {len(failed)} · 총 {(time.time()-t0)/3600:.2f}h")
    for r in failed:
        log(f"  ✗ v={r['v']:g} seed={r['seed']} rc={r['rc']}")
    log("다음: $PY -m bdbot.cli health --all      (스텝 커버리지 확인)")
    log("      $PY scratch/trap_drag_ensemble.py  (앙상블 재분석)")
    log("=" * 76)
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
