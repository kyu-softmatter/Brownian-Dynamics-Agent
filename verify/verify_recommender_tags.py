#!/usr/bin/env python
"""Snapshot bdbot.interactions' recommender output for every intake case.

The recommender keys on **both English and Korean** keywords, because
`intake/*/observation.yaml` was written in Korean. Removing the Korean keywords is
only safe if the tags and the top recommendation are unchanged on every real case --
so the check is a before/after comparison of exactly that, not an inspection of the
keyword list.

    $PY recbase.py            # print + write baseline.json
    $PY recbase.py --check    # compare against baseline.json
"""
import glob
import json
import sys
from pathlib import Path

import yaml

HERE = str(Path(__file__).with_name("recommender_tags_baseline.json"))
sys.path.insert(0, ".")
from bdbot.interactions import infer_tags, recommend    # noqa: E402


def snapshot() -> dict:
    out = {}
    for p in sorted(glob.glob("intake/*/observation.yaml")):
        raw = yaml.safe_load(open(p, encoding="utf-8"))
        scored, tags = recommend(raw, top=3)
        out[p.split("/")[1]] = {
            "tags": sorted(infer_tags(raw)),
            "top": [it.key for it, _, _ in scored],
            "scores": [round(s, 2) for _, s, _ in scored],
        }
    return out


def main() -> int:
    cur = snapshot()
    if "--check" in sys.argv:
        base = json.load(open(HERE))
        bad = [k for k in sorted(set(base) | set(cur)) if base.get(k) != cur.get(k)]
        for k in sorted(cur):
            mark = "✗" if k in bad else "✓"
            print(f"  {mark} {k:24} tags={cur[k]['tags']}")
            if k in bad:
                print(f"      was {base.get(k)}")
                print(f"      now {cur[k]}")
        print(f"{len(cur)} cases, {len(bad)} changed")
        return 1 if bad else 0
    json.dump(cur, open(HERE, "w"), indent=1, ensure_ascii=False)
    for k in sorted(cur):
        print(f"  {k:24} tags={cur[k]['tags']}")
        print(f"  {'':24} top={cur[k]['top']} scores={cur[k]['scores']}")
    print(f"baseline written: {len(cur)} cases")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
