"""지식베이스 — 두 갈래를 함께 읽는다.

    runs/*/record.json    런 경험 (tools/postmortem.py 가 만든다)
    kb/entries/*.json     ★ 런이 없는 지식 — 스케치 해석·문헌 증류·도구 교훈

마스터플랜 §7은 SQLite+FTS5를 계획하지만, 런 100개 미만에서는 과잉이다.
**데이터 형식만 먼저 고정**하고 아파지면 DB로 옮긴다 (§7.0).

★ `kb/entries/` 가 왜 필요한가 (2026-08-04 추가):
  런에서 나온 교훈만 저장하고 있었다. 그런데 **앞단(스케치 → 물리계)에서 나오는 지식도
  재사용 가능하고 잊으면 손해다.** 예: "R=5µm 관례를 다른 스케치에 승계했더니 abp-rod에서
  τ_R과 160배 어긋났다", "논문 추출은 논문 자체 수치를 재현해 검증한다",
  "정규식으로 YAML을 패치하면 조용히 엉뚱한 노드를 고친다".
  이런 것들은 런이 없어서 record.json 에 넣을 자리가 없었다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY tools/kb.py list
    $PY tools/kb.py query --tags 2D,harmonic_trap
    $PY tools/kb.py query --kind intake
    $PY tools/kb.py lessons
    $PY tools/kb.py add --claim "..." --kind pitfall --origin intake \
                        --source "intake/abp-rod/observation.yaml#D1" --tags abp,convention
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent


ENTRY_DIR = ROOT / "kb" / "entries"
SCHEMA_ENTRY = "bdbot.kb_entry/0.1"
# 런이 아닌 지식의 출처 종류 (마스터플랜 §5.2 Source.kind 와 정렬)
ORIGINS = ("intake", "paper", "tooling", "method", "handbook", "user_input")


def load_all():
    """런 기록 + 런 없는 엔트리를 같은 모양으로 합쳐 돌려준다."""
    out = []
    for f in sorted((ROOT / "runs").glob("*/record.json")):
        try:
            r = json.loads(f.read_text())
            r.setdefault("origin", "our_run")
            out.append(r)
        except Exception as e:
            print(f"  (건너뜀 {f}: {e})", file=sys.stderr)
    for f in sorted(ENTRY_DIR.glob("*.json")):
        try:
            out.append(json.loads(f.read_text()))
        except Exception as e:
            print(f"  (건너뜀 {f}: {e})", file=sys.stderr)
    return out


def add_entry(claim, kind, origin, source, tags, coords, not_verified) -> Path:
    """런 없는 지식 엔트리를 파일 하나로 남긴다. record.json 과 같은 필드 이름을 쓴다."""
    if origin not in ORIGINS:
        raise SystemExit(f"origin 은 {ORIGINS} 중 하나여야 합니다 (받은 값: {origin})")
    ENTRY_DIR.mkdir(parents=True, exist_ok=True)
    slug = "".join(c if c.isalnum() or c in "-_" else "-" for c in claim.lower())[:48]
    slug = "-".join(x for x in slug.split("-") if x) or "entry"
    path = ENTRY_DIR / f"{origin}__{slug}.json"
    n = 2
    while path.exists():
        path = ENTRY_DIR / f"{origin}__{slug}-{n}.json"
        n += 1
    rec = {
        "schema": SCHEMA_ENTRY,
        "run_id": None,                 # 런이 없다 — 그게 이 엔트리의 존재 이유
        "origin": origin,
        "case": None,
        "outcome": None,
        "system_tags": tags,
        "dimensionless": coords,
        "source": source,
        "lessons": [{"claim": claim, "kind": kind, "coords": coords}],
        "not_verified": not_verified,
    }
    path.write_text(json.dumps(rec, indent=2, ensure_ascii=False))
    return path


def match_coord(rec, expr):
    """'키=lo:hi' 또는 '키=값' (10% 허용). 키는 부분 일치."""
    key, rng = expr.split("=", 1)
    key = key.strip()
    hits = [(k, v) for k, v in rec.get("dimensionless", {}).items() if key in k]
    if not hits:
        return False
    _, val = hits[0]
    if ":" in rng:
        lo, hi = (float(x) for x in rng.split(":"))
        return lo <= val <= hi
    t = float(rng)
    return abs(val - t) <= 0.1 * abs(t)


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    q = sub.add_parser("query")
    q.add_argument("--tags", default="")
    q.add_argument("--coord", action="append", default=[])
    q.add_argument("--outcome", default="")
    q.add_argument("--kind", default="", help="교훈 종류 (pitfall/method_note/...)")
    q.add_argument("--origin", default="", help=f"출처 종류 {ORIGINS + ('our_run',)}")
    sub.add_parser("lessons").add_argument("--origin", default="")
    a_ = sub.add_parser("add", help="런 없는 지식 엔트리 추가 (스케치 해석·문헌·도구 교훈)")
    a_.add_argument("--claim", required=True)
    a_.add_argument("--kind", default="method_note",
                    help="parameter | phase_boundary | scaling | method_note | pitfall")
    a_.add_argument("--origin", required=True, help=f"{ORIGINS}")
    a_.add_argument("--source", default="", help="어디서 나왔는가 (파일#앵커 · 논문 locator)")
    a_.add_argument("--tags", default="")
    a_.add_argument("--coord", action="append", default=[],
                    help="무차원 좌표 '이름=값' (검색용)")
    a_.add_argument("--not-verified", action="append", default=[])
    a = ap.parse_args()

    if a.cmd == "add":
        coords = {}
        for c in a.coord:
            k, v = c.split("=", 1)
            coords[k.strip()] = float(v)
        path = add_entry(a.claim, a.kind, a.origin, a.source,
                         [t.strip() for t in a.tags.split(",") if t.strip()],
                         coords, a.not_verified)
        print(f"→ {path.relative_to(ROOT)}")
        return 0

    recs = load_all()
    if not recs:
        print("기록이 없습니다. tools/postmortem.py (런) 또는 kb.py add (런 없는 지식) 를 쓰세요.")
        return 1

    def rid(r):
        return r.get("run_id") or f"[{r.get('origin', '?')}] {r.get('source') or '—'}"

    if a.cmd == "list":
        print(f"{'출처':<44}{'origin':<11}{'outcome':<10}{'tags'}")
        print("-" * 104)
        for r in sorted(recs, key=lambda x: (x.get("origin") != "our_run", rid(x))):
            print(f"{rid(r)[:43]:<44}{r.get('origin', '?'):<11}"
                  f"{str(r.get('outcome') or '—'):<10}"
                  f"{','.join((r.get('system_tags') or [])[:4])}")
        n_run = sum(1 for r in recs if r.get("origin") == "our_run")
        print(f"\n총 {len(recs)}건  (런 {n_run} · 런 없는 지식 {len(recs) - n_run})")
        return 0

    if a.cmd == "lessons":
        n = 0
        for r in recs:
            if a.origin and r.get("origin") != a.origin:
                continue
            for l_ in r.get("lessons", []):
                n += 1
                c = f"   coords={l_['coords']}" if l_.get("coords") else ""
                tier = l_.get("tier", r.get("tier", 3))
                print(f"[tier{tier} {l_['kind']}] {l_['claim']}{c}")
                print(f"      ← {rid(r)}"
                      + (f" ({r['outcome']})" if r.get("outcome") else
                         f" · origin={r.get('origin')}"))
        print(f"\n총 {n}건")
        return 0

    sel = recs
    if a.tags:
        want = {t.strip() for t in a.tags.split(",")}
        sel = [r for r in sel if want <= set(r.get("system_tags") or [])]
    for c in a.coord:
        sel = [r for r in sel if match_coord(r, c)]
    if a.outcome:
        sel = [r for r in sel if r.get("outcome") == a.outcome]
    if a.origin:
        sel = [r for r in sel if r.get("origin") == a.origin]
    if a.kind:
        sel = [r for r in sel
               if any(l_.get("kind") == a.kind for l_ in r.get("lessons", []))]

    for r in sel:
        print("=" * 78)
        print(f"{rid(r)}   [{r.get('outcome') or r.get('origin')}]")
        if r.get("system_tags"):
            print(f"  tags: {', '.join(r['system_tags'])}")
        if r.get("dimensionless"):
            print("  무차원수:")
            for k, v in r["dimensionless"].items():
                print(f"    {k:<36} {v:.4g}")
        if r.get("observables"):
            print("  관측량 (측정 vs 예측):")
            for o in r["observables"]:
                pv = o.get("predicted")
                if pv is None:
                    print(f"    {o['name']:<16} {o['measured']:>13.6g} {o['unit']:<8} (예측 없음)")
                else:
                    print(f"    {o['name']:<16} {o['measured']:>13.6g} vs {pv:>13.6g} "
                          f"{o['unit']:<8} {o['err_pct']:+7.2f}%")
        if r.get("lessons"):
            print("  교훈:")
            for l_ in r["lessons"]:
                print(f"    [{l_['kind']}] {l_['claim']}")
        if r.get("not_verified"):
            print("  미검증:")
            for nv in r["not_verified"]:
                print(f"    · {nv}")
    print(f"\n{len(sel)}/{len(recs)} 건 일치")
    return 0


if __name__ == "__main__":
    sys.exit(main())
