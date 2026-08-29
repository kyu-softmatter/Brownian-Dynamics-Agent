"""L5 이관이 `trap-drag` 의 관측량 정의를 바꿨는가 — 옛 런 vs 새 런 대조.

## 왜 이 대조가 가능한가

v-스윕 재실행(2026-08-06)에서 대부분의 레거시 런은 **run_id 가 같아 덮였습니다**.
그런데 `v=0.5` 그룹만은 살아남았습니다 — L2 `system.yaml` 의 `external.drag_velocity` 가
단일값 `0.5` → 7점 리스트로 바뀌면서 **스펙 해시가 달라졌기** 때문입니다.
그래서 이 그룹만 옛/새 세대가 **같은 시드로 공존**하고, 유일한 대조 창입니다.

## 무엇을 시험하는가

`params`·`numerics` 가 0 차이면 **같은 물리 · 같은 시드 · 같은 dt** 입니다.
그러면 궤적이 같아야 하고, 관측량도 같아야 합니다. 다르면 **관측량 정의가 바뀐 것**입니다.
(정의가 같다면 v-스윕 정정이 코드 변경 탓이 아니라는 근거가 됩니다.)

⚠️ 이 대조는 **v=0.5 에만** 가능합니다. 다른 속도의 레거시는 덮여 사라졌습니다.

    $PY scratch/verify_trap_drag_generations.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
RTOL = 1e-12
PASS, FAIL = [], []


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASS if ok else FAIL).append(name)
    print(f"  {'✓' if ok else '✗'} {name}" + (f"\n        {detail}" if detail else ""))


def is_new(m: dict) -> bool:
    return "step_drift_max_sigma" in m.get("numerics", {})


def flat(d, pre=""):
    out = {}
    items = d.items() if isinstance(d, dict) else enumerate(d)
    for k, v in items:
        p = f"{pre}{k}"
        if isinstance(v, (dict, list)) and not isinstance(v, str):
            out.update(flat(v, p + "."))
        else:
            out[p] = v
    return out


def numeric_leaves(d, pre=""):
    out = {}
    for k, v in d.items():
        p = f"{pre}{k}"
        if isinstance(v, dict):
            out.update(numeric_leaves(v, p + "."))
        elif isinstance(v, (int, float)) and not isinstance(v, bool):
            out[p] = float(v)
    return out


def main() -> int:
    print("=" * 80)
    print("trap-drag 세대 대조 — L5 이관이 관측량 정의를 바꿨는가 (v=0.5 그룹)")
    print("=" * 80)

    pairs: dict[int, dict] = {}
    for p in sorted((ROOT / "runs").glob("trap-drag-2d-hex300__tr0.117647-s*")):
        if not (p / "metrics.json").exists() or not (p / "spec.json").exists():
            continue
        m = json.loads((p / "metrics.json").read_text())
        spec = json.loads((p / "spec.json").read_text())
        pairs.setdefault(int(spec["numerics"]["seed"]), {})[
            "new" if is_new(m) else "old"] = (p, m, spec)

    both = {s: g for s, g in pairs.items() if len(g) == 2}
    check("옛/새 세대가 짝을 이루는 시드가 있다", len(both) >= 2,
          f"짝 {len(both)}개 (시드 {sorted(both)})")
    if not both:
        print("\n대조 불가 — 레거시가 남아 있지 않습니다.")
        print("=" * 80)
        return 1

    # ── ① 해시 대상 스펙 차이가 '문서'뿐인가 ──────────────────────────
    print("\n[①] 스펙 차이 — params·numerics 가 같아야 '같은 물리'다")
    s0 = both[sorted(both)[0]]
    for sect, must_match in (("params", True), ("numerics", True), ("system", False)):
        fo, fn = flat(s0["old"][2].get(sect, {})), flat(s0["new"][2].get(sect, {}))
        diff = [k for k in set(fo) | set(fn) if fo.get(k, "∅") != fn.get(k, "∅")]
        if must_match:
            check(f"{sect} 차이 0건", not diff, f"차이: {sorted(diff)[:6]}")
        else:
            keys = sorted({k.split(".")[0] + "." + k.split(".")[1]
                           for k in diff if "." in k})
            check(f"{sect} 차이는 drag_velocity 선언뿐",
                  all("drag_velocity" in k for k in diff) if diff else True,
                  f"{len(diff)}건, 최상위: {keys[:4]}")

    # ── ② 관측량 스칼라 ───────────────────────────────────────────────
    print("\n[②] `metrics.result` 스칼라 — 정의가 같으면 비트 동일해야 한다")
    worst = {}
    for s in sorted(both):
        lo = numeric_leaves(both[s]["old"][1].get("result", {}))
        ln = numeric_leaves(both[s]["new"][1].get("result", {}))
        common = set(lo) & set(ln)
        for k in common:
            a, b = lo[k], ln[k]
            if not (math.isfinite(a) and math.isfinite(b)):
                continue
            rel = abs(a - b) / max(abs(a), abs(b), 1e-30)
            if k not in worst or rel > worst[k][0]:
                worst[k] = (rel, a, b, s)
        if s == sorted(both)[0]:
            check("옛/새 result 필드 집합이 같다", set(lo) == set(ln),
                  f"공통 {len(common)} · 옛만 {sorted(set(lo)-set(ln))[:4]} · "
                  f"새만 {sorted(set(ln)-set(lo))[:4]}")
    bad = {k: v for k, v in worst.items() if v[0] > RTOL}
    check(f"스칼라 {len(worst)}개가 전 시드에서 비트 동일 (<{RTOL:g})", not bad,
          "불일치: " + ", ".join(f"{k}={v[0]:.2e}" for k, v in list(bad.items())[:5])
          if bad else f"{len(worst)}개 × 시드 {len(both)}개 전부 일치")

    # ── ③ 스텝 수 ─────────────────────────────────────────────────────
    print("\n[③] 스텝 수 — 관측량이 같으면 적분도 같아야 한다")
    steps_ok, speed = True, []
    for s in sorted(both):
        o, n = both[s]["old"][1], both[s]["new"][1]
        so, sn = o["numerics"].get("steps_done"), n["numerics"].get("steps_done")
        steps_ok &= (so == sn)
        if o.get("wall_seconds") and n.get("wall_seconds"):
            speed.append(o["wall_seconds"] / n["wall_seconds"])
    check("스텝 수가 전 시드에서 동일", steps_ok,
          f"steps_done = {both[sorted(both)[0]]['old'][1]['numerics'].get('steps_done'):,}")
    if speed:
        print(f"        (벽시계 배속 {min(speed):.2f}~{max(speed):.2f}× — "
              f"관측량은 같고 속도만 달라졌다)")

    # ── ④ observables.npz 배열 ────────────────────────────────────────
    print("\n[④] `observables.npz` 배열 — 스칼라만 같고 배열이 다를 수도 있다")
    arr_bad, n_arr, n_seed = [], 0, 0
    for s in sorted(both):
        po, pn = both[s]["old"][0], both[s]["new"][0]
        if not ((po / "observables.npz").exists() and (pn / "observables.npz").exists()):
            continue
        n_seed += 1
        with np.load(po / "observables.npz") as zo, np.load(pn / "observables.npz") as zn:
            if set(zo.files) != set(zn.files):
                arr_bad.append(f"s{s}: 키 불일치")
                continue
            for k in sorted(set(zo.files)):
                a, b = zo[k], zn[k]
                if a.shape != b.shape:
                    arr_bad.append(f"s{s}.{k}: shape {a.shape}≠{b.shape}")
                    continue
                if not np.issubdtype(a.dtype, np.number):
                    continue
                n_arr += 1
                d = np.abs(a.astype(float) - b.astype(float))
                sc = float(np.nanmax(np.abs(a.astype(float)))) or 1.0
                if float(np.nanmax(d) if d.size else 0.0) / sc > RTOL:
                    arr_bad.append(f"s{s}.{k}: 상대차 {np.nanmax(d)/sc:.2e}")
    check(f"배열 {n_arr}개 (시드 {n_seed}개분) 전부 일치", not arr_bad,
          "; ".join(arr_bad[:5]) if arr_bad else f"불일치 0")

    print("\n" + "=" * 80)
    print(f"통과 {len(PASS)} · 실패 {len(FAIL)}")
    for f in FAIL:
        print(f"   ✗ {f}")
    if not FAIL:
        print("✓ PASS — L5 이관은 관측량 정의를 바꾸지 않았다.")
        print("  → v-스윕 정정(결함 v-의존성·회복률 비단조)은 코드 변경이 아니라")
        print("     단일 런 + 과소평가된 오차 때문이다.")
        print("  ⚠️ 이 대조는 v=0.5 에만 가능하다 — 다른 속도의 레거시는 덮여 사라졌다.")
    print("=" * 80)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
