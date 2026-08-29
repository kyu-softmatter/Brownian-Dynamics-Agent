"""GSD 재생으로 max|F| 를 측정할 수 있는가 — **가드 진값이 있는 런에서 교정한다.**

레거시 런 79개의 스텝 커버리지를 채우려면 원래 벽시계 80시간을 다시 써야 합니다.
그런데 궤적(GSD)이 저장돼 있으니, 프레임마다 힘만 다시 계산하면 `dt·|F|max` 를
**재시뮬레이션 없이** 얻을 수 있습니다 (프레임 200~5667개 × 힘 1회).

⚠️ 다만 GSD 프레임은 가드 표본보다 드뭅니다 — 최댓값을 **과소평가**할 수 있습니다.
   그래서 먼저 **둘 다 있는 런**에서 재생값 vs 가드 진값을 대조합니다.
   교정 없이 쓰면 "측정했다"면서 틀린 값을 넣게 됩니다 (규칙 6).

    $PY scratch/probe_gsd_replay.py runs/<run_id>
"""
from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from bdbot import nondim as ND, run as RUN  # noqa: E402


def _import_case_builders() -> None:
    """`@RUN.builder` 데코레이터가 등록되게 케이스 모듈을 임포트한다."""
    import importlib
    for mod in ("trap_drag_2d", "chain_bend_2d", "trap_2d_5um", "abp_rod_2d", "soft_r3_2d"):
        try:
            importlib.import_module(mod)
        except Exception as e:
            print(f"  (케이스 {mod} 임포트 실패: {type(e).__name__}: {e})")


def replay(run_dir: Path, max_frames: int | None = None, verbose: bool = True) -> dict:
    import gsd.hoomd

    spec = ND.load(run_dir / "spec.json")
    build_fn = RUN.get_builder(spec.case)
    dt_star = float(spec.numerics["dt_star"])

    gsd_files = sorted(run_dir.glob("*.gsd"))
    if not gsd_files:
        raise FileNotFoundError(f"{run_dir} 에 GSD 가 없습니다")

    with tempfile.TemporaryDirectory() as td:
        b = build_fn(spec, Path(td))
        sim = b.sim
        n_build = int(sim.state.N_particles)

        with gsd.hoomd.open(str(gsd_files[0]), "r") as traj:
            n_frames = len(traj)
            idx = range(n_frames) if max_frames is None else \
                np.unique(np.linspace(0, n_frames - 1, max_frames).astype(int))
            f_max = 0.0
            n_used = 0
            per_frame = []
            for i in idx:
                fr = traj[int(i)]
                pos = np.asarray(fr.particles.position, dtype=float)
                if pos.shape[0] != n_build:
                    raise ValueError(
                        f"입자 수 불일치: GSD {pos.shape[0]} vs build {n_build} "
                        f"— 이 런은 재생할 수 없습니다")
                snap = sim.state.get_snapshot()
                snap.particles.position[:] = pos
                if getattr(fr.particles, "orientation", None) is not None \
                        and snap.particles.orientation is not None:
                    o = np.asarray(fr.particles.orientation, dtype=float)
                    if o.shape == snap.particles.orientation.shape:
                        snap.particles.orientation[:] = o
                sim.state.set_snapshot(snap)
                sim.run(0)                        # 힘 재계산 (상태는 전진하지 않음)
                fm = 0.0
                for f in b.forces:
                    arr = np.asarray(f.forces, dtype=float)
                    if arr.size:
                        fm = max(fm, float(np.abs(arr).max()))
                per_frame.append(fm)
                f_max = max(f_max, fm)
                n_used += 1

    return {"case": spec.case, "run_id": spec.run_id, "dt_star": dt_star,
            "n_frames_total": n_frames, "n_frames_used": n_used,
            "f_max_gsd": f_max, "drift_gsd": dt_star * f_max,
            "per_frame": np.asarray(per_frame)}


def main() -> int:
    if len(sys.argv) < 2:
        print(__doc__)
        return 3
    run_dir = Path(sys.argv[1])
    if not run_dir.is_absolute():
        run_dir = ROOT / run_dir
    _import_case_builders()

    print("=" * 78)
    print(f"GSD 재생 교정 — {run_dir.name}")
    print("=" * 78)

    r = replay(run_dir)
    m = json.loads((run_dir / "metrics.json").read_text())
    num = m.get("numerics", {})
    truth_drift = num.get("step_drift_max_sigma")
    truth_f = num.get("f_max_kT_per_sigma")

    print(f"  케이스        {r['case']}")
    print(f"  dt*           {r['dt_star']:.6e}")
    print(f"  GSD 프레임    {r['n_frames_used']} / {r['n_frames_total']}")
    print(f"  재생 |F|max   {r['f_max_gsd']:.6g} kT/σ")
    print(f"  재생 drift    {r['drift_gsd']:.6e}")
    if truth_drift is None:
        print("  가드 진값     없음 (이 런은 교정에 못 씀 — 진값이 있는 런으로 하세요)")
        return 0
    print(f"  가드 |F|max   {truth_f:.6g} kT/σ   ← 진값")
    print(f"  가드 drift    {truth_drift:.6e}   ← 진값")
    ratio = r["drift_gsd"] / truth_drift if truth_drift else float("nan")
    print()
    print(f"  재생/가드 = {ratio:.4f}×   "
          + ("(재생이 과소평가)" if ratio < 0.95 else
             "(일치)" if ratio <= 1.05 else "(재생이 더 큼 — 가드가 놓친 순간을 잡음)"))
    pf = r["per_frame"]
    print(f"  프레임별 |F| 분포: 중앙 {np.median(pf):.4g} · 90% {np.percentile(pf, 90):.4g} "
          f"· 최대 {pf.max():.4g}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
