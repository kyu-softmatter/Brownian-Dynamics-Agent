"""`chain-bend-2d-dlvo` 의 완료 런에 `result.txt`(완료 마커)를 소급해서 채운다.

★ 왜 필요했나: `result.txt` 는 이 프로젝트의 **완료 마커**인데 `bdbot.run.execute` 가
  아니라 **케이스 스크립트가** 쓴다 (soft_r3_2d·abp_rod_2d·trap_2d_5um 셋 다 그렇다).
  `chain_bend_dlvo_2d.py` 를 새로 쓰면서 이걸 빠뜨려서 셋이 조용히 깨졌다:
    ① `bdbot.cli status` 가 137런을 **0개**로 셌다 (cli.py:133)
    ② `runid.prepare_outdir` 가 완료를 못 알아봐 같은 런을 계속 재실행했다
    ③ "미완료 정리" 를 result.txt 기준으로 돌렸다가 **완료 런 6개를 지웠다**
  케이스 쪽은 고쳤고(2026-08-06), 그 전에 끝난 런은 이 스크립트로 채운다.

    $PY scratch/backfill_result_txt.py [--dry]
"""
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
dry = "--dry" in sys.argv

done = skip = 0
for d in sorted((ROOT / "runs").glob("chain-bend-2d-dlvo__*")):
    mp, rp = d / "metrics.json", d / "result.txt"
    if not mp.exists() or rp.exists():
        skip += 1
        continue
    mj = json.loads(mp.read_text())
    lines = []
    for o in mj.get("observables", []):
        m, p_ = o.get("measured"), o.get("predicted")
        tail = f"   (예측 {p_:.6g})" if isinstance(p_, (int, float)) else ""
        lines.append(f"  {o['name']:<22} {m:.6g}{tail}" if m is not None
                     else f"  {o['name']:<22} —")
    l4 = d / "l4.json"
    verdict = ""
    if l4.exists():
        try:
            verdict = f"L4: {json.loads(l4.read_text()).get('verdict', '?')}"
        except Exception:
            pass
    body = "\n".join(["=" * 84, f"결과 — {d.name}", "=" * 84, *lines,
                      "=" * 84, verdict,
                      "(★ result.txt 는 scratch/backfill_result_txt.py 로 소급 생성. "
                      "원본 리포트는 이 런에 없다 — 케이스가 당시 안 썼다)"])
    report = (d / "report.txt").read_text() if (d / "report.txt").exists() else ""
    if not dry:
        rp.write_text((report + "\n" if report else "") + body)
    done += 1

print(f"채움 {done}건 · 건너뜀 {skip}건{' (dry-run)' if dry else ''}")
