"""`chain-bend` omega sweep -> K'(omega) and K''(omega) figure.

Comparing against `driven_static_stiffness` at De ~ 1 is this case's
`implementation_check` -- that prediction is **derived from the model I implemented**
(the stiffness matrix A plus finite-stiffness traps), so a mismatch is a bug, not a
discovery (CLAUDE.md rule 7').

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
        print("no completed runs", file=sys.stderr); return 1
    rows.sort(key=lambda r: r["De"])
    De = np.array([r["De"] for r in rows])
    Kp = np.array([r["K_prime"] for r in rows])
    Kpp = np.array([r["K_doubleprime"] for r in rows])
    Ke = np.array([r["K_sem"] for r in rows])
    Ks = rows[0]["K_static_pred"]

    fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.8))
    a = ax[0]
    a.errorbar(De, Kp, yerr=Ke, fmt="o-", ms=7, capsize=3, color="tab:blue", label="K' (storage)")
    a.errorbar(De, Kpp, yerr=Ke, fmt="s--", ms=6, capsize=3, color="tab:red", label="K'' (loss)")
    a.axhline(Ks, color="k", ls=":", lw=1.6, label=f"static prediction {Ks:.0f}")
    a.axvline(1, color="0.6", lw=1, ls=":")
    a.set_xscale("log"); a.set_yscale("log")
    for x in (a.xaxis, a.yaxis): x.set_major_formatter(pl); x.set_minor_formatter(NullFormatter())
    a.set(xlabel="De = omega tau_max", ylabel="K [kT/d^2]",
          title="(1) complex stiffness")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[1]
    lo = De < 1.5
    if lo.any():
        rel = 100 * (Kp[lo] - Ks) / Ks
        a.errorbar(De[lo], rel, yerr=100 * Ke[lo] / Ks, fmt="o-", ms=8, capsize=4, color="tab:green")
    a.axhline(0, color="k", ls="--", lw=1.4)
    a.axhspan(-15, 15, color="green", alpha=.10, label="tolerance ±15%")
    a.set_xscale("log"); a.xaxis.set_major_formatter(pl); a.xaxis.set_minor_formatter(NullFormatter())
    a.set(xlabel="De", ylabel="(K' - predicted)/predicted  [%]",
          title="(2) ★ implementation_check -- De<1.5 static limit")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[2]
    nom = np.array([r.get("nominal_vs_measured_pct") or np.nan for r in rows])
    a.plot(De, np.abs(nom), "o-", ms=7, color="tab:purple")
    a.set_xscale("log"); a.set_yscale("log")
    for x in (a.xaxis, a.yaxis): x.set_major_formatter(pl); x.set_minor_formatter(NullFormatter())
    a.set(xlabel="De", ylabel="|nominal - measured| / |measured|  [%]",
          title="(3) the error that using the nominal amplitude would have caused (ZOH)")
    a.grid(alpha=.3, which="both")

    fig.suptitle(f"chain-bend-2d-oscill -- K*(omega), {len(rows)} points (10 cycles)",
                 fontsize=12)
    fig.tight_layout()
    p = ROOT / "runs/chain-bend-2d-oscill__SWEEP"; p.mkdir(parents=True, exist_ok=True)
    fig.savefig(p / "kstar.png", dpi=145); plt.close(fig)
    print(f"{'De':>9}{'K1':>12}{'±':>9}{'K2':>12}{'vs pred %':>11}")
    for r_, k1, k2, e in zip(rows, Kp, Kpp, Ke):
        print(f"{r_['De']:>9.3f}{k1:>12.1f}{e:>9.1f}{k2:>12.1f}"
              f"{100*(k1-Ks)/Ks:>11.1f}")
    print(f"\n{p/'kstar.png'}")
    return 0

if __name__ == "__main__":
    sys.exit(main())
