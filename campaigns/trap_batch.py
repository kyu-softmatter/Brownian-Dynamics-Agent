"""첫 손그림 배치 실행 — 모호성 A1(차원) 판별 + dt 래더를 한 번에.

    d ∈ {2, 3}  ×  dt* ∈ {5e-3, 2.5e-3}  ×  seed ∈ {1,2,3,4}  =  16 런

동시 실행 k=8 (config/run_policy.yaml 기본값). HOOMD 는 단일스레드이므로
독립 런 동시 실행이 유일한 병렬화 경로다.

usage:
    python scripts/trap_batch.py <run_dir>
    python scripts/trap_batch.py <run_dir> --worker <dim> <dt_star> <seed> <out>
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from dataclasses import asdict
from pathlib import Path

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

PY = sys.executable
DIMS = (2, 3)
DTS = (5.0e-3, 2.5e-3)
SEEDS = (1, 2, 3, 4)
CONCURRENCY = 8


def worker(dim: int, dt_star: float, seed: int, out: Path) -> None:
    from simbot.run import TrapRunConfig, run_trap
    cfg = TrapRunConfig(dim=dim, dt_star=dt_star, seed=seed,
                        label=f"d{dim}_dt{dt_star:g}_s{seed}")
    r = run_trap(cfg, outdir=out)
    print(json.dumps({
        "label": cfg.label, "dim": dim, "dt_star": dt_star, "seed": seed,
        "var_c": r.var_per_component_star, "var_c_se": r.var_per_component_se,
        "var_r": r.var_radial_star, "var_r_se": r.var_radial_se,
        "kT_conf": r.kT_conf_star, "kT_conf_se": r.kT_conf_se,
        "wall_s": r.wall_s, "guards_ok": r.guards["finite"],
        "max_disp": r.guards["max_step_displacement_l_trap"],
        "n_snap": r.n_independent_snapshots,
        "outdir": str(out),
    }))


def main() -> None:
    run_dir = Path(sys.argv[1])
    if "--worker" in sys.argv:
        i = sys.argv.index("--worker")
        worker(int(sys.argv[i + 1]), float(sys.argv[i + 2]),
               int(sys.argv[i + 3]), Path(sys.argv[i + 4]))
        return

    jobs = [(d, dt, s) for d in DIMS for dt in DTS for s in SEEDS]
    print(f"# {len(jobs)} 런, 동시 {CONCURRENCY}", flush=True)
    t0 = time.perf_counter()
    results, pending = [], list(jobs)
    running: list[tuple[subprocess.Popen, tuple]] = []

    while pending or running:
        while pending and len(running) < CONCURRENCY:
            d, dt, s = pending.pop(0)
            out = run_dir / "raw" / f"d{d}_dt{dt:g}_s{s}"
            p = subprocess.Popen(
                [PY, __file__, str(run_dir), "--worker", str(d), str(dt), str(s),
                 str(out)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
            running.append((p, (d, dt, s)))
        done = [(p, j) for p, j in running if p.poll() is not None]
        for p, j in done:
            running.remove((p, j))
            so, se = p.communicate()
            line = next((l for l in so.splitlines() if l.strip().startswith("{")), None)
            if line is None:
                print(f"  !! {j} 실패\n{se[-1500:]}", flush=True)
                continue
            r = json.loads(line)
            results.append(r)
            print(f"  [{len(results):2d}/{len(jobs)}] {r['label']:<18} "
                  f"<x*^2>={r['var_c']:.5f}±{r['var_c_se']:.5f}  "
                  f"kT_conf={r['kT_conf']:.5f}  {r['wall_s']:.1f}s", flush=True)
        if running:
            time.sleep(0.2)

    wall = time.perf_counter() - t0
    run_dir.mkdir(parents=True, exist_ok=True)
    (run_dir / "05_batch_results.json").write_text(
        json.dumps({"jobs": results, "batch_wall_s": wall,
                    "concurrency": CONCURRENCY}, indent=2))
    print(f"\n# 배치 완료: {wall:.1f}s (동시 {CONCURRENCY})", flush=True)


if __name__ == "__main__":
    main()
