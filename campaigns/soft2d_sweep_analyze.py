"""soft-r3 `A` 스윕 분석 — 시드 앙상블 오차막대 + Zahn 상도 대조.

usage: python scripts/soft2d_sweep_analyze.py <sweep_dir> [--json out.json]

## 이 스크립트가 조심하는 것

**① 시드 앙상블이 오차의 정직한 추정치다.** 프레임간 산포는 시간 상관을 포함하므로
   통계오차가 아니다. `aggregate_seeds` 로 시드 평균의 SE 를 쓴다.

**② 결함 분율만 보고하지 않는다.** `A=100` 무작위 시작은 `{5:2%, 6:96%, 7:2%}` 로
   5·7 이 동수인 **전위** 서명이고, `A=0.1` 은 3~12 로 퍼진 **액체**다. 같은 분율
   숫자가 전혀 다른 물리를 뜻한다 → 배위수 분포와 `5-7 대칭성`을 함께 낸다.

**③ 평형화 판정을 대신하지 않는다.** 전반/후반 절반의 `ψ₆` 를 비교해 **표류가 있는지
   보고만** 한다. "평형에 도달했다"는 임계값이 있어야 할 수 있는 말이다.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.analysis.structure import (hex_order, rdf, structure_factor,
                                       zahn_phase)
from simbot.analysis.trap import aggregate_seeds


def load(d: Path) -> dict:
    man = json.loads((d / "manifest.json").read_text())
    z = np.load(d / "samples.npz")
    return {"man": man, "traj": z["traj"], "energy": z["energy"],
            "Lx": float(z["box"][0]), "Ly": float(z["box"][1])}


def per_run(r: dict) -> dict:
    """런 1개의 측정값. 후반 절반만 쓰고, 전반과의 차이도 함께 낸다."""
    traj, Lx, Ly = r["traj"], r["Lx"], r["Ly"]
    n = len(traj)
    first, second = traj[: n // 2], traj[n // 2:]
    h2 = hex_order(second, Lx=Lx, Ly=Ly)
    h1 = hex_order(first, Lx=Lx, Ly=Ly)
    s = structure_factor(second, Lx=Lx, Ly=Ly, n_max=18)
    g = rdf(second, Lx=Lx, Ly=Ly, bins=200)
    hist = h2.coordination_hist
    n5, n7 = hist.get(5, 0.0), hist.get(7, 0.0)
    return {
        "psi6_global": h2.psi6_global, "psi6_local": h2.psi6_local_mean,
        "defect_fraction": h2.defect_fraction,
        "sixfold": s.sixfold_modulation, "first_peak_g": g.first_peak_g,
        "first_peak_r": g.first_peak_r,
        "energy_pp": float(r["energy"][n // 2:].mean()) / r["man"]["config"]["n_particles"],
        # ★ 전위 서명: 5 와 7 이 동수인가
        "coord_5": n5, "coord_7": n7,
        "five_seven_balance": (abs(n5 - n7) / max(n5 + n7, 1e-12)),
        "coord_spread": float(len([k for k, v in hist.items() if v > 0.005])),
        "coordination_hist": {int(k): float(v) for k, v in sorted(hist.items())},
        # 표류 (평형화 진단 — 판정 아님)
        "psi6_drift": h2.psi6_global - h1.psi6_global,
        "min_separation": r["man"]["guards"]["min_separation"],
        "min_sep_over_r_min": r["man"]["guards"]["min_separation_over_r_min"],
        "force_displacement": r["man"]["guards"]["force_displacement_star"],
        "wall_s": r["man"]["wall_s"],
    }


def main() -> int:
    root = Path(sys.argv[1])
    dirs = sorted(p for p in root.glob("A*_*_s*") if (p / "samples.npz").exists())
    if not dirs:
        print(f"⛔ {root} 에 런이 없다", file=sys.stderr)
        return 2

    groups: dict[tuple[float, str], list[dict]] = {}
    for d in dirs:
        cfg = json.loads((d / "manifest.json").read_text())["config"]
        key = (float(cfg["amplitude"]), cfg["init"])
        groups.setdefault(key, []).append(per_run(load(d)))

    print(f"# soft-r3 A 스윕 — {len(dirs)} 런 / {len(groups)} 조건\n")
    metrics = ("psi6_global", "psi6_local", "defect_fraction", "sixfold",
               "energy_pp", "first_peak_g")
    out: dict = {}

    hdr = f"{'A':>6} {'init':<7} {'n':>2} " + " ".join(
        f"{m[:11]:>17}" for m in metrics)
    print(hdr)
    print("-" * len(hdr))
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        runs = groups[key]
        agg = {m: aggregate_seeds([r[m] for r in runs]) for m in metrics}
        out[f"A{A:g}_{init}"] = {
            "amplitude": A, "init": init, "n_seeds": len(runs),
            "zahn": zahn_phase(A),
            **{m: {"mean": agg[m].mean, "se": agg[m].se, "spread": agg[m].spread,
                   "values": agg[m].values} for m in metrics},
            "coordination_hist_seed1": runs[0]["coordination_hist"],
            "five_seven_balance": [r["five_seven_balance"] for r in runs],
            "psi6_drift": [r["psi6_drift"] for r in runs],
            "min_sep_over_r_min": min(r["min_sep_over_r_min"] for r in runs),
            "wall_s_total": sum(r["wall_s"] for r in runs),
        }
        cells = " ".join(f"{agg[m].mean:>9.4f}±{agg[m].se:<7.4f}" for m in metrics)
        print(f"{A:>6g} {init:<7} {len(runs):>2} {cells}")

    # --- Zahn 상도 대조 ---
    print("\n## Zahn 상도 대조  [출처, 미재현]")
    print(f"{'A':>6} {'Γ':>9} {'Zahn 예측':<18} {'init':<7} "
          f"{'ψ₆ 전역':>17} {'6겹변조':>17} {'우리 판독':<12}")
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        o = out[f"A{A:g}_{init}"]
        p6, s6 = o["psi6_global"], o["sixfold"]
        # 판독은 **관측 기반**으로만. Zahn 경계를 쓰지 않는다
        if p6["mean"] > 0.7 and s6["mean"] > 0.7:
            verdict = "결정-유사"
        elif s6["mean"] > 0.3:
            verdict = "6겹 잔류"
        else:
            verdict = "등방-유사"
        print(f"{A:>6g} {o['zahn']['gamma']:>9.2f} {o['zahn']['phase_zahn']:<18} "
              f"{init:<7} {p6['mean']:>9.4f}±{p6['se']:<7.4f} "
              f"{s6['mean']:>9.4f}±{s6['se']:<7.4f} {verdict:<12}")

    # --- 초기조건 의존성 ---
    print("\n## 초기조건 의존성 (hex vs random) — 시드 앙상블 대비")
    print(f"{'A':>6} {'ψ₆ hex':>17} {'ψ₆ random':>17} {'차이':>10} {'결합SE':>9} {'σ':>7}  판정")
    for A in sorted({k[0] for k in groups}):
        h, r = out.get(f"A{A:g}_hex"), out.get(f"A{A:g}_random")
        if not (h and r):
            continue
        dh, dr = h["psi6_global"], r["psi6_global"]
        diff = dh["mean"] - dr["mean"]
        se = float(np.hypot(dh["se"] or 0.0, dr["se"] or 0.0))
        sig = abs(diff) / se if se > 0 else float("inf")
        # ★ 3σ 이내면 "구별 안 됨" — 같다는 증명이 아니다
        mark = "구별 안 됨" if sig < 3 else "❗유의한 차이"
        print(f"{A:>6g} {dh['mean']:>9.4f}±{dh['se']:<7.4f} "
              f"{dr['mean']:>9.4f}±{dr['se']:<7.4f} {diff:>+10.4f} {se:>9.4f} "
              f"{sig:>6.2f}σ  {mark}")

    # --- 결함의 성격 ---
    print("\n## 결함의 성격 — 분율만으로는 구별되지 않는다")
    print(f"{'A':>6} {'init':<7} {'결함 분율':>17} {'5-7 불균형':>11} "
          f"{'배위수 종류':>10}  배위수 분포 (시드 1)")
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        o = out[f"A{A:g}_{init}"]
        df = o["defect_fraction"]
        bal = float(np.mean(o["five_seven_balance"]))
        hist = o["coordination_hist_seed1"]
        kinds = len([k for k, v in hist.items() if v > 0.005])
        top = " ".join(f"{k}:{v:.2f}" for k, v in hist.items() if v > 0.005)
        note = "전위(5-7 쌍)" if bal < 0.25 and df["mean"] < 0.2 else ""
        print(f"{A:>6g} {init:<7} {df['mean']:>9.4f}±{df['se']:<7.4f} "
              f"{bal:>11.3f} {kinds:>10d}  {top}  {note}")

    # --- 진단 (판정 아님) ---
    print("\n## 진단 — 판정이 아니라 사실 보고")
    for key in sorted(groups, key=lambda k: (k[0], k[1])):
        A, init = key
        o = out[f"A{A:g}_{init}"]
        drift = np.array(o["psi6_drift"])
        print(f"  A={A:<6g} {init:<7} ψ₆ 표류(후−전) = "
              f"{drift.mean():+.4f} ± {drift.std(ddof=1)/np.sqrt(drift.size):.4f}"
              f"   최소분리/r_min = {o['min_sep_over_r_min']:.2f}"
              f"   wall {o['wall_s_total']:.0f}s")
    print("\n  ⚠ 표류가 0 과 구별되지 않는다는 것이 '평형에 도달했다'는 뜻은 아니다 —")
    print("    이 런 길이로는 그 표류를 볼 수 없다는 뜻이다.")

    if "--json" in sys.argv:
        p = Path(sys.argv[sys.argv.index("--json") + 1])
        p.write_text(json.dumps(out, indent=2, ensure_ascii=False, default=float))
        print(f"\n→ {p}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
