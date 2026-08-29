"""Can max|F| be measured by replaying a GSD? -- **calibrated on a run that has the
guard's ground truth.**

Filling in the step coverage of the 79 legacy runs would mean spending the original
80 hours of wall clock again. But the trajectories (GSD) are stored, so recomputing
only the forces frame by frame gives `dt*|F|max` **without re-simulating**
(200-5667 frames x one force evaluation each).

⚠️ The catch is that GSD frames are sparser than the guard's samples, so the maximum
   can be **underestimated**. So the replayed value is first compared against the
   guard's ground truth on a run that has **both**.
   Using it without that calibration would mean entering a wrong number while
   claiming to have measured it (rule 6).

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
    """Import the case modules so their `@RUN.builder` decorators register."""
    import importlib
    for mod in ("trap_drag_2d", "chain_bend_2d", "trap_2d_5um", "abp_rod_2d", "soft_r3_2d"):
        try:
            importlib.import_module(mod)
        except Exception as e:
            print(f"  (case {mod} failed to import: {type(e).__name__}: {e})")


def replay(run_dir: Path, max_frames: int | None = None, verbose: bool = True) -> dict:
    import gsd.hoomd

    spec = ND.load(run_dir / "spec.json")
    build_fn = RUN.get_builder(spec.case)
    dt_star = float(spec.numerics["dt_star"])

    gsd_files = sorted(run_dir.glob("*.gsd"))
    if not gsd_files:
        raise FileNotFoundError(f"no GSD in {run_dir}")

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
                        f"particle-count mismatch: GSD {pos.shape[0]} vs "
                        f"build {n_build} -- this run cannot be replayed")
                snap = sim.state.get_snapshot()
                snap.particles.position[:] = pos
                if getattr(fr.particles, "orientation", None) is not None \
                        and snap.particles.orientation is not None:
                    o = np.asarray(fr.particles.orientation, dtype=float)
                    if o.shape == snap.particles.orientation.shape:
                        snap.particles.orientation[:] = o
                sim.state.set_snapshot(snap)
                sim.run(0)                        # recompute forces (the state does not advance)
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
    print(f"GSD replay calibration -- {run_dir.name}")
    print("=" * 78)

    r = replay(run_dir)
    m = json.loads((run_dir / "metrics.json").read_text())
    num = m.get("numerics", {})
    truth_drift = num.get("step_drift_max_sigma")
    truth_f = num.get("f_max_kT_per_sigma")

    print(f"  case          {r['case']}")
    print(f"  dt*           {r['dt_star']:.6e}")
    print(f"  GSD frames    {r['n_frames_used']} / {r['n_frames_total']}")
    print(f"  replay |F|max {r['f_max_gsd']:.6g} kT/sigma")
    print(f"  replay drift  {r['drift_gsd']:.6e}")
    if truth_drift is None:
        print("  guard truth   none (this run cannot calibrate -- use one that has it)")
        return 0
    print(f"  guard |F|max  {truth_f:.6g} kT/sigma   <- ground truth")
    print(f"  guard drift   {truth_drift:.6e}   <- ground truth")
    ratio = r["drift_gsd"] / truth_drift if truth_drift else float("nan")
    print()
    print(f"  replay/guard = {ratio:.4f}x   "
          + ("(replay underestimates)" if ratio < 0.95 else
             "(agrees)" if ratio <= 1.05 else
             "(replay is larger -- it caught an instant the guard missed)"))
    pf = r["per_frame"]
    print(f"  per-frame |F| distribution: median {np.median(pf):.4g}, "
          f"90th {np.percentile(pf, 90):.4g} "
          f", max {pf.max():.4g}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
