"""Retroactively write `result.txt` (the completion marker) into completed
`chain-bend-2d-dlvo` runs.

★ Why this was needed: `result.txt` is this project's **completion marker**, but it
  is written by **the case script**, not by `bdbot.run.execute` (all three of
  soft_r3_2d, abp_rod_2d and trap_2d_5um do it themselves).
  Writing `chain_bend_dlvo_2d.py` omitted it, and three things broke silently:
    (1) `bdbot.cli status` counted 137 runs as **0** (cli.py:133)
    (2) `runid.prepare_outdir` could not recognise completion and kept re-running
        the same run
    (3) a "clean up the incomplete runs" pass keyed on result.txt **deleted 6
        completed runs**
  The case side is fixed (2026-08-06); runs that finished before that are filled in
  by this script.

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
        tail = f"   (prediction {p_:.6g})" if isinstance(p_, (int, float)) else ""
        lines.append(f"  {o['name']:<22} {m:.6g}{tail}" if m is not None
                     else f"  {o['name']:<22} —")
    l4 = d / "l4.json"
    verdict = ""
    if l4.exists():
        try:
            verdict = f"L4: {json.loads(l4.read_text()).get('verdict', '?')}"
        except Exception:
            pass
    body = "\n".join(["=" * 84, f"result — {d.name}", "=" * 84, *lines,
                      "=" * 84, verdict,
                      "(★ this result.txt was generated retroactively by "
                      "scratch/backfill_result_txt.py. The original report does not "
                      "exist for this run -- the case did not write one at the time)"])
    report = (d / "report.txt").read_text() if (d / "report.txt").exists() else ""
    if not dry:
        rp.write_text((report + "\n" if report else "") + body)
    done += 1

print(f"filled {done} · skipped {skip} runs{' (dry-run)' if dry else ''}")
