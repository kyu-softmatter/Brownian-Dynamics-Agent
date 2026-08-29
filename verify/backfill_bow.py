"""완료된 chain-bend-2d-dlvo 런에 `bow`(굽음) 를 GSD 궤적에서 사후 계산해 채운다.

★ 왜 사후계산이 가능한가: 굽음은 위치만으로 정해지는 기하량이라 궤적에서 정확히 재현된다.
  (힘 기반 량이었다면 불가능했다 — CLAUDE.md 의 GSD 재생 무효 사례 참조: 시간 의존 구동에서
   힘을 재계산하면 트랩 앵커가 t=0 에 고정돼 16배 과대평가된다. 굽음은 그 함정이 없다.)

⚠️ 한계: GSD 는 표본이 성기다(런당 ~200 프레임 vs 원본 표본 2000개). 그래서 여기서 채운
  값은 `bow_rms_gsd` 로 **이름을 다르게** 넣는다 — 신규 런의 `bow_rms`(전 표본)와 섞으면
  분해능이 다른 값을 같은 이름으로 비교하게 된다. 락인(bow_drive)은 프레임이 적어 신뢰도가
  낮으므로 채우지 않는다.

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
        "note": f"굽음 RMS — GSD 궤적 {len(traj)}프레임에서 사후계산 "
                f"(scratch/backfill_bow.py). 신규 런의 bow_rms(전 표본)보다 성기다"})
    if not dry:
        mp.write_text(json.dumps(met, indent=2, ensure_ascii=False))
    n_done += 1
    print(f"  {d.name[:62]:<62} bow_rms_gsd={rms:.5f} d  ({len(traj)}f)")
print(f"\n채움 {n_done}건 · 건너뜀 {n_skip}건{' (dry-run, 안 씀)' if dry else ''}")
