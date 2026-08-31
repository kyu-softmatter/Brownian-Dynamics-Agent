"""Did the L5 migration change `trap-drag`'s observable definitions? -- old runs
compared against new.

## Why this comparison is possible at all

In the v-sweep re-run (2026-08-06) most legacy runs were **overwritten, because their
run_id was unchanged**. Only the `v=0.5` group survived -- because L2
`system.yaml`'s `external.drag_velocity` changed from the single value `0.5` to a
7-point list, which **changed the spec hash**.
So only in that group do the old and new generations **coexist under the same seed**,
and it is the only window for this comparison.

## What is being tested

If `params` and `numerics` differ in 0 places, then it is **the same physics, the same
seed and the same dt**. The trajectory must then be the same, and so must the
observables. If they differ, **an observable definition changed**.
(If the definitions agree, that is grounds for saying the v-sweep correction was not
caused by a code change.)

⚠️ This comparison is possible **only for v=0.5**. The legacy runs at other velocities
were overwritten and are gone.

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
    print("trap-drag generation comparison -- did the L5 migration change any "
          "observable definition? (the v=0.5 group)")
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
    check("there are seeds where the old and new generations pair up", len(both) >= 2,
          f"{len(both)} pairs (seeds {sorted(both)})")
    if not both:
        print("\ncomparison impossible -- no legacy runs remain.")
        print("=" * 80)
        return 1

    # ── (1) is the hashed spec difference confined to 'documentation'? ────────
    print("\n[1] spec differences -- params and numerics must agree for it to be "
          "'the same physics'")
    s0 = both[sorted(both)[0]]
    for sect, must_match in (("params", True), ("numerics", True), ("system", False)):
        fo, fn = flat(s0["old"][2].get(sect, {})), flat(s0["new"][2].get(sect, {}))
        diff = [k for k in set(fo) | set(fn) if fo.get(k, "∅") != fn.get(k, "∅")]
        if must_match:
            check(f"{sect}: 0 differences", not diff,
                  f"differences: {sorted(diff)[:6]}")
        else:
            keys = sorted({k.split(".")[0] + "." + k.split(".")[1]
                           for k in diff if "." in k})
            check(f"{sect}: the only difference is the drag_velocity declaration",
                  all("drag_velocity" in k for k in diff) if diff else True,
                  f"{len(diff)} differences, top-level: {keys[:4]}")

    # ── (2) observable scalars ──────────────────────────────────────
    print("\n[2] `metrics.result` scalars -- if the definitions agree these must be "
          "bit-identical")
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
            check("the old and new result field sets are the same", set(lo) == set(ln),
                  f"common {len(common)} . old only {sorted(set(lo)-set(ln))[:4]} . "
                  f"new only {sorted(set(ln)-set(lo))[:4]}")
    bad = {k: v for k, v in worst.items() if v[0] > RTOL}
    check(f"{len(worst)} scalars bit-identical across every seed (<{RTOL:g})", not bad,
          "mismatched: " + ", ".join(f"{k}={v[0]:.2e}"
                                     for k, v in list(bad.items())[:5])
          if bad else f"all {len(worst)} x {len(both)} seeds agree")

    # ── (3) step counts ─────────────────────────────────────────────
    print("\n[3] step counts -- if the observables agree, the integration must too")
    steps_ok, speed = True, []
    for s in sorted(both):
        o, n = both[s]["old"][1], both[s]["new"][1]
        so, sn = o["numerics"].get("steps_done"), n["numerics"].get("steps_done")
        steps_ok &= (so == sn)
        if o.get("wall_seconds") and n.get("wall_seconds"):
            speed.append(o["wall_seconds"] / n["wall_seconds"])
    check("step counts identical across every seed", steps_ok,
          f"steps_done = {both[sorted(both)[0]]['old'][1]['numerics'].get('steps_done'):,}")
    if speed:
        print(f"        (wall-clock speedup {min(speed):.2f}-{max(speed):.2f}x -- "
              f"the observables are unchanged, only the speed differs)")

    # ── (4) observables.npz arrays ──────────────────────────────────
    print("\n[4] `observables.npz` arrays -- the scalars could agree while the arrays "
          "differ")
    arr_bad, n_arr, n_seed = [], 0, 0
    for s in sorted(both):
        po, pn = both[s]["old"][0], both[s]["new"][0]
        if not ((po / "observables.npz").exists() and (pn / "observables.npz").exists()):
            continue
        n_seed += 1
        with np.load(po / "observables.npz") as zo, np.load(pn / "observables.npz") as zn:
            if set(zo.files) != set(zn.files):
                arr_bad.append(f"s{s}: key mismatch")
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
                    arr_bad.append(f"s{s}.{k}: rel diff {np.nanmax(d)/sc:.2e}")
    check(f"all {n_arr} arrays (over {n_seed} seeds) agree", not arr_bad,
          "; ".join(arr_bad[:5]) if arr_bad else f"0 mismatches")

    print("\n" + "=" * 80)
    print(f"passed {len(PASS)} . failed {len(FAIL)}")
    for f in FAIL:
        print(f"   ✗ {f}")
    if not FAIL:
        print("✓ PASS -- the L5 migration did not change any observable definition.")
        print("  -> the v-sweep correction came from single runs plus underestimated")
        print("     errors, not from a code change.")
        print("  ⚠️ Possible only for v=0.5; legacy runs at other v were overwritten.")
    print("=" * 80)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
