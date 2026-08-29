"""Backfill `bow` into completed chain-bend-2d-dlvo runs, computed after the fact
from the GSD trajectory.

★ Why post-hoc computation is valid here: bow is a geometric quantity fixed by
  positions alone, so it is reproduced exactly from the trajectory.
  (It would have been impossible for a force-based quantity -- see the invalid
   GSD-replay case in CLAUDE.md: under time-dependent driving, recomputing forces
   pins the trap anchor at t=0 and overestimates by 16x. Bow has no such trap.)

⚠️ Limitation: GSD is sparsely sampled (~200 frames per run vs 2000 original
  samples). So the value filled in here goes in under a **different name**,
  `bow_rms_gsd` -- mixing it with a new run's `bow_rms` (all samples) would mean
  comparing different resolutions under one name. The lock-in (bow_drive) is not
  backfilled, because too few frames make it unreliable.

    $PY scratch/backfill_bow.py [--dry]
"""
import json
import sys
from pathlib import Path

import gsd.hoomd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))
from chain_bend_dlvo_2d import _bow  # noqa: E402

dry = "--dry" in sys.argv
n_done = n_skip = 0
for d in sorted((ROOT / "runs").glob("chain-bend-2d-dlvo__*")):
    mp, gp, sp = d / "metrics.json", d / "traj_A.gsd", d / "spec.json"
    if not (mp.exists() and gp.exists() and sp.exists()):
        n_skip += 1
        continue
    met = json.loads(mp.read_text())
    names = {o["name"] for o in met["observables"]}
    if "bow_rms_gsd" in names:
        n_skip += 1
        continue
    n = int(json.loads(sp.read_text())["params"]["n_beads"])
    traj = gsd.hoomd.open(str(gp), "r")
    bows = np.array([_bow(traj[i].particles.position[:n, 1]) for i in range(len(traj))])
    rms = float(np.sqrt((bows ** 2).mean()))
    met["observables"].append({
        "name": "bow_rms_gsd", "measured": rms, "predicted": None, "unit": "d",
        "err_pct": None, "err_sigma": None, "sigma": None,
        "prediction_source": "none", "role": "measurement", "scope": "composite",
        "derivation": "", "tol_pct": None, "tol_sigma": None,
        "note": f"bow RMS -- computed after the fact from {len(traj)} GSD frames "
                f"(scratch/backfill_bow.py). Sparser than a new run's bow_rms "
                f"(all samples)"})
    if not dry:
        mp.write_text(json.dumps(met, indent=2, ensure_ascii=False))
    n_done += 1
    print(f"  {d.name[:62]:<62} bow_rms_gsd={rms:.5f} d  ({len(traj)}f)")
print(f"\nfilled {n_done} · skipped {n_skip}"
      f"{' (dry-run, nothing written)' if dry else ''}")
