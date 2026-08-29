"""`chain-bend` ω 스윕 → K'(ω)·K''(ω) 그림.

De≈1 에서 `driven_static_stiffness` 와 대조하는 것이 이 케이스의 `implementation_check`
입니다 — 그 예측은 **내가 구현한 모델(강성 행렬 A + 유한강성 트랩)에서 유도**된 것이므로
불일치는 발견이 아니라 버그입니다 (CLAUDE.md 규칙 7').

    $PY scratch/chain_bend_sweep_plot.py
"""
from __future__ import annotations
import json, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path[:0] = [str(ROOT), str(ROOT / "cases")]
from bdbot import nondim as ND  # noqa: E402

def main() -> int:
    import matplotlib; matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, NullFormatter
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    pl = FuncFormatter(lambda x, _: f"{x:g}")

    rows = []
    for d in sorted((ROOT / "runs").glob("chain-bend-2d-oscill__*")):
        if not (d / "metrics.json").exists():
            continue
        m = json.loads((d / "metrics.json").read_text())
        r = m.get("result")
        if not r:
            continue
        rows.append(r)
    if not rows:
        print("완료된 런이 없습니다", file=sys.stderr); return 1
    rows.sort(key=lambda r: r["De"])
    De = np.array([r["De"] for r in rows])
    Kp = np.array([r["K_prime"] for r in rows])
    Kpp = np.array([r["K_doubleprime"] for r in rows])
    Ke = np.array([r["K_sem"] for r in rows])
    Ks = rows[0]["K_static_pred"]

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
    a = ax[0]
    a.errorbar(De, Kp, yerr=Ke, fmt="o-", ms=7, capsize=3, color="tab:blue", label="K′ (저장)")
    a.errorbar(De, Kpp, yerr=Ke, fmt="s--", ms=6, capsize=3, color="tab:red", label="K″ (손실)")
    a.axhline(Ks, color="k", ls=":", lw=1.6, label=f"정적 예측 {Ks:.0f}")
    a.axvline(1, color="0.6", lw=1, ls=":")
    a.set_xscale("log"); a.set_yscale("log")
    for x in (a.xaxis, a.yaxis): x.set_major_formatter(pl); x.set_minor_formatter(NullFormatter())
    a.set(xlabel="De = ω τ_max", ylabel="K [kT/d²]", title="① 복소 강성")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[1]
    lo = De < 1.5
    if lo.any():
        rel = 100 * (Kp[lo] - Ks) / Ks
        a.errorbar(De[lo], rel, yerr=100 * Ke[lo] / Ks, fmt="o-", ms=8, capsize=4, color="tab:green")
    a.axhline(0, color="k", ls="--", lw=1.4)
    a.axhspan(-15, 15, color="green", alpha=.10, label="허용 ±15%")
    a.set_xscale("log"); a.xaxis.set_major_formatter(pl); a.xaxis.set_minor_formatter(NullFormatter())
    a.set(xlabel="De", ylabel="(K′ − 예측)/예측  [%]",
          title="② ★ implementation_check — De<1.5 정적 극한")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[2]
    nom = np.array([r.get("nominal_vs_measured_pct") or np.nan for r in rows])
    a.plot(De, np.abs(nom), "o-", ms=7, color="tab:purple")
    a.set_xscale("log"); a.set_yscale("log")
    for x in (a.xaxis, a.yaxis): x.set_major_formatter(pl); x.set_minor_formatter(NullFormatter())
    a.set(xlabel="De", ylabel="|공칭 − 측정| / |측정|  [%]",
          title="③ 공칭 진폭을 썼다면 생겼을 오차 (ZOH)")
    a.grid(alpha=.3, which="both")

    fig.suptitle(f"chain-bend-2d-oscill — K*(ω) · {len(rows)}점 (10주기)", fontsize=12)
    fig.tight_layout()
    p = ROOT / "runs/chain-bend-2d-oscill__SWEEP"; p.mkdir(parents=True, exist_ok=True)
    fig.savefig(p / "kstar.png", dpi=145); plt.close(fig)
    print(f"{'De':>9}{'K′':>12}{'±':>9}{'K″':>12}{'예측대비%':>11}")
    for r_, k1, k2, e in zip(rows, Kp, Kpp, Ke):
        print(f"{r_['De']:>9.3f}{k1:>12.1f}{e:>9.1f}{k2:>12.1f}"
              f"{100*(k1-Ks)/Ks:>11.1f}")
    print(f"\n{p/'kstar.png'}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
