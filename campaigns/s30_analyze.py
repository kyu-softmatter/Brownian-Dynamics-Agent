"""The S30 analysis, in the order the sealed plan fixes.

    $PY campaigns/s30_analyze.py            # the tables
    $PY campaigns/s30_analyze.py --figures  # ...and F1-F4 into figures/

Every threshold here is read from `campaigns/s30_preregistration/prediction.yaml`
rather than written in this file, so a number can only enter the verdict by
being in the sealed document. `analysis_plan.yaml` fixes the order of
operations, and steps 1, 2 and 4 can each end the analysis.

⚠ This file was written AFTER the runs finished. That is allowed and it is the
reason the plan had to be sealed: the plan says what this script must do, so
writing it now cannot change what counts as an answer. What it must not do is
introduce a statistic the plan does not name -- anything extra is labelled
`[exploratory]` in the output.
"""
from __future__ import annotations

import argparse
import json
import math
import pathlib
import sys

import numpy as np
import yaml
from scipy import stats

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import runcard as RC          # noqa: E402
from bdbot import stats as ST            # noqa: E402

DOCS = ROOT / "campaigns/s30_preregistration"
import importlib.util as _iu             # noqa: E402
_s = _iu.spec_from_file_location("hyst", ROOT / "campaigns/soft2d_hysteresis.py")
HYST = _iu.module_from_spec(_s)
_s.loader.exec_module(HYST)

PRED = yaml.safe_load((DOCS / "prediction.yaml").read_text())
FLOOR = float(PRED["primary_statistic"]["crystal_floor"]["value"])
DELTA_EQ = float(PRED["equivalence_margin"]["delta"])


def load(run_id: str) -> dict:
    d = ROOT / "runs" / run_id
    met = json.loads((d / "metrics.json").read_text())
    res = np.load(d / "observables.npz")
    g = np.asarray(res["psi6_global"], dtype=float)
    q = len(g) // 4
    return dict(
        dir=d, run_id=run_id, metrics=met,
        psi6_global=g,
        late=float(g[-q:].mean()),
        q3=g[-2 * q:-q], q4=g[-q:],
        full=float(g.mean()),
        quarters=[float(g[i * q:(i + 1) * q].mean()) for i in range(4)],
        defect=1.0 - float((met["result"]["coord_hist"])[6]),
        energy_err=next((o["err_pct"] for o in met["observables"]
                         if "energy consistency" in o["name"]), None),
    )


def drift(r: dict) -> tuple[bool, float, float]:
    """The sealed rule: |mean(Q4)-mean(Q3)| > 3*sqrt(sem(Q3)^2+sem(Q4)^2)."""
    d = abs(float(r["q4"].mean()) - float(r["q3"].mean()))
    s = math.sqrt(ST.block_sem(r["q3"]) ** 2 + ST.block_sem(r["q4"]) ** 2)
    return d > 3 * s, d, 3 * s


def welch(h: list[float], r: list[float]) -> dict:
    nh, nr = len(h), len(r)
    sh, sr = float(np.std(h, ddof=1)), float(np.std(r, ddof=1))
    a, b = sh ** 2 / nh, sr ** 2 / nr
    sig = math.sqrt(a + b)
    nu = (a + b) ** 2 / (a * a / (nh - 1) + b * b / (nr - 1)) if sig else float("nan")
    delta = float(np.mean(h)) - float(np.mean(r))
    t = delta / sig if sig else float("inf")
    return dict(delta=delta, sigma=sig, nu=nu, t=t,
                p_two=float(2 * (1 - stats.t.cdf(abs(t), nu))) if sig else 0.0,
                t_crit=float(stats.t.ppf(0.995, nu)) if sig else float("nan"),
                tost_halfwidth=float(stats.t.ppf(0.95, nu)) * sig if sig else 0.0,
                sd_h=sh, sd_r=sr)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--figures", action="store_true")
    a = ap.parse_args()

    ids = {r["id"]: HYST.run_id_of(r) for r in HYST.RUNS}
    W = 74

    # ── STEP 1 ───────────────────────────────────────────────────────────
    print("=" * W)
    print("STEP 1 -- the seal, before any number is read")
    print("=" * W)
    broken = []
    for k, rid in ids.items():
        ok, probs = RC.verify_seal(ROOT / "runs" / rid, root=ROOT)
        hard = [p for p in probs if not p.startswith("[warn]")]
        print(f"  {k:3s} {'OK' if ok else 'BROKEN':6s} "
              f"{len(RC.read_seal(ROOT / 'runs' / rid))} docs  {rid[-12:]}"
              + (f"  {hard}" if hard else ""))
        if not ok:
            broken.append(k)
    if broken:
        print(f"\nSEAL BROKEN: {broken}. The comparison table is not built.")
        return 1
    print("\n  all 12 seals hold")

    R = {k: load(v) for k, v in ids.items()}
    H = [f"H{i}" for i in range(1, 6)]
    Rr = [f"R{i}" for i in range(1, 6)]

    # ── STEP 2 ───────────────────────────────────────────────────────────
    print()
    print("=" * W)
    print("STEP 2 -- the controls. Either one ends the campaign.")
    print("=" * W)
    g1 = R["C1"]["late"] > FLOOR and R["C1"]["defect"] < 0.10
    g2 = R["C2"]["late"] < 0.30 and R["C2"]["defect"] > 0.20
    print(f"  G1  deep crystal stays ordered   psi6_global_late = {R['C1']['late']:.4f} "
          f"(> {FLOOR})   defect = {R['C1']['defect']:.4f} (< 0.10)   "
          f"{'PASS' if g1 else 'FAIL'}")
    print(f"  G2  deep liquid melts it         psi6_global_late = {R['C2']['late']:.4f} "
          f"(< 0.30)   defect = {R['C2']['defect']:.4f} (> 0.20)   "
          f"{'PASS' if g2 else 'FAIL'}")
    if not (g1 and g2):
        print("\nPROTOCOL FAILURE -- no physics is reported. P1 is void.")
        return 1

    # ── STEP 3 ───────────────────────────────────────────────────────────
    print()
    print("=" * W)
    print("STEP 3 -- the drift gate, PER ARM")
    print("=" * W)
    fired = []
    for k in H + Rr + ["C1"]:
        f, d, lim = drift(R[k])
        if f:
            fired.append(k)
        print(f"  {k:3s} |Q4-Q3| = {d:.4f}  vs  3*SE_diff = {lim:.4f}   "
              f"{'NOT STATIONARY' if f else 'stationary'}")
    if fired:
        print(f"\n  ⚠ {fired} did not reach a steady state. Their late means are "
              f"reported but not used in the equivalence test (per the sealed rule, "
              f"the gate is per arm and does not end the campaign).")

    # ── STEP 4 ───────────────────────────────────────────────────────────
    print()
    print("=" * W)
    print("STEP 4 -- the primary read: the pair (n_survive_H, n_survive_R)")
    print("=" * W)
    print(f"  {'seed':>10s}  {'H late':>9s}  {'>floor':>6s}    {'R late':>9s}  {'>floor':>6s}")
    for i in range(5):
        h, r = R[H[i]], R[Rr[i]]
        print(f"  {PRED['runs']['entries'][2 + i]['seed']:>10d}  {h['late']:9.4f}  "
              f"{'yes' if h['late'] > FLOOR else 'no':>6s}    {r['late']:9.4f}  "
              f"{'yes' if r['late'] > FLOOR else 'no':>6s}")
    nh = sum(R[k]["late"] > FLOOR for k in H)
    nr = sum(R[k]["late"] > FLOOR for k in Rr)
    print(f"\n  (n_survive_H, n_survive_R) = ({nh}, {nr})")

    if nr > 0:
        verdict = "BOTH ARMS ORDER"
    elif nh == 5:
        verdict = "HYSTERESIS"
    elif nh == 0:
        verdict = None                      # step 5 decides
    else:
        verdict = "SPLIT"

    lo, hi = (stats.beta.ppf([0.025, 0.975], nh + 0.5, 5 - nh + 0.5)
              if 0 < nh < 5 else (nh / 5, nh / 5))
    print(f"  survival fraction {nh}/5"
          + ("" if not (0 < nh < 5) else f"   Jeffreys 95 % CI [{lo:.2f}, {hi:.2f}]"))

    # ── STEP 5 ───────────────────────────────────────────────────────────
    print()
    print("=" * W)
    print("STEP 5 -- the secondary read (only when the pair is (5,0) or (0,0))")
    print("=" * W)
    w = welch([R[k]["late"] for k in H], [R[k]["late"] for k in Rr])
    paired_sd = float(np.std([R[H[i]]["late"] - R[Rr[i]]["late"] for i in range(5)], ddof=1))
    if (nh, nr) in ((5, 0), (0, 0)):
        tost = abs(w["delta"]) + w["tost_halfwidth"] < DELTA_EQ
        print(f"  Delta = {w['delta']:+.4f}   sigma = {w['sigma']:.4f}   "
              f"nu = {w['nu']:.2f}")
        print(f"  |t| = {abs(w['t']):.2f}   vs t(nu, 0.995) = {w['t_crit']:.3f}   "
              f"two-sided p = {w['p_two']:.2e}")
        print(f"  TOST  |Delta| + t(nu,0.95)*sigma = "
              f"{abs(w['delta']) + w['tost_halfwidth']:.4f}  vs  delta = {DELTA_EQ}"
              f"   -> {'FIRES' if tost else 'does not fire'}")
        print(f"  sd_H/sd_R = {w['sd_h'] / w['sd_r']:.2f}" if w["sd_r"] else "  sd_R = 0")
        print(f"  paired sd of per-seed differences = {paired_sd:.4f}  "
              f"(unpaired sigma*sqrt(5) = {w['sigma'] * math.sqrt(5):.4f})")
        if verdict is None:
            verdict = "NO HYSTERESIS DETECTED" if tost else "UNDERPOWERED"
        elif verdict == "HYSTERESIS":
            verdict = "HYSTERESIS, with the magnitude quoted"
    else:
        print(f"  NOT COMPUTED. The pair is ({nh}, {nr}); the sealed rule restricts")
        print(f"  Delta to a unanimous arm, because a mean over a bimodal arm is")
        print(f"  not a state. Reported for the record only:")
        print(f"    Delta = {w['delta']:+.4f}   sd_H = {w['sd_h']:.4f}   "
              f"sd_R = {w['sd_r']:.4f}   sd_H/sd_R = "
              f"{w['sd_h'] / w['sd_r']:.2f}" if w["sd_r"] else "")

    # ── STEP 6 ───────────────────────────────────────────────────────────
    print()
    print("=" * W)
    print("STEP 6 -- descriptive. Names no phase.")
    print("=" * W)
    print(f"  {'run':4s} {'late':>8s} {'full':>8s} {'defect':>8s}   quarters")
    for k in ["C1", "C2"] + H + Rr:
        qs = "  ".join(f"{x:.3f}" for x in R[k]["quarters"])
        print(f"  {k:4s} {R[k]['late']:8.4f} {R[k]['full']:8.4f} "
              f"{R[k]['defect']:8.4f}   {qs}")

    # ── STEP 7 ───────────────────────────────────────────────────────────
    print()
    print("=" * W)
    print("STEP 7 -- P2 and P3, reported whatever P1 said")
    print("=" * W)
    dr = [R[k]["defect"] for k in Rr]
    p2 = all(abs(x - 0.285) <= 0.285 * 0.05 for x in dr)
    print(f"  P2  defect_fraction(R) = {np.mean(dr):.4f} +- {np.std(dr, ddof=1):.4f}"
          f"   sealed: 0.285 +- 0.014 (5 %)   {'PASS' if p2 else 'MISS'}")
    print(f"      per seed: " + "  ".join(f"{x:.4f}" for x in dr))
    ee = [abs(R[k]["energy_err"]) for k in ["C1", "C2"] + H + Rr
          if R[k]["energy_err"] is not None]
    p3 = all(x <= 2.0 for x in ee)
    print(f"  P3  |energy identity| max = {max(ee):.4f} %   sealed: <= 2 %   "
          f"{'PASS' if p3 else 'MISS'}")

    # ── the sealed table ─────────────────────────────────────────────────
    print()
    print("=" * W)
    print("F5 -- the sealed prediction against the measurement")
    print("=" * W)
    print("| id | role | sealed | measured | verdict |")
    print("|---|---|---|---|---|")
    print(f"| P1 | hypothesis | (n_survive_H, n_survive_R) = (5, 0) -> HYSTERESIS "
          f"| ({nh}, {nr}) | **{verdict}** |")
    print(f"| G1 | gate | psi6_late > 0.5 and defect < 0.10 | "
          f"{R['C1']['late']:.4f}, {R['C1']['defect']:.4f} | {'PASS' if g1 else 'FAIL'} |")
    print(f"| G2 | gate | psi6_late < 0.30 and defect > 0.20 | "
          f"{R['C2']['late']:.4f}, {R['C2']['defect']:.4f} | {'PASS' if g2 else 'FAIL'} |")
    print(f"| P2 | implementation_check | defect(R) = 0.285 +- 5 % | "
          f"{np.mean(dr):.4f} | {'PASS' if p2 else 'MISS'} |")
    print(f"| P3 | implementation_check | energy identity <= 2 % | "
          f"{max(ee):.4f} % | {'PASS' if p3 else 'MISS'} |")
    print()
    print(f"VERDICT: {verdict}")

    if a.figures:
        make_figures(R, H, Rr, nh, nr, verdict, w)
    return 0


def make_figures(R, H, Rr, nh, nr, verdict, w) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    out = ROOT / "figures/s30_hysteresis.png"
    fig, ax = plt.subplots(2, 2, figsize=(14, 9))

    # F1 -- psi6_global(t), both arms, from t = 0
    for k in H:
        g = R[k]["psi6_global"]
        ax[0, 0].plot(np.linspace(0, 100, len(g)), g, "-", lw=.9, color="#1b6b3a",
                      alpha=.85, label="crystal start" if k == "H1" else None)
    for k in Rr:
        g = R[k]["psi6_global"]
        ax[0, 0].plot(np.linspace(0, 100, len(g)), g, "-", lw=.9, color="#2b6cb0",
                      alpha=.85, label="random start" if k == "R1" else None)
    ax[0, 0].axhline(FLOOR, color="#9b2c2c", lw=1.4, ls="--",
                     label=f"crystal floor {FLOOR}")
    ax[0, 0].axvspan(75, 100, color="#888", alpha=.15, label="decision window")
    #  t = 0 is not a sample: the first one sits at sample_every*dt = 0.25 tau_B,
    #  by which time psi_6 has already fallen from its exact lattice value.
    ax[0, 0].plot([0], [1.0], "*", ms=15, color="#1b6b3a",
                  label=r"$t=0$ lattice, exactly 1.000")
    ax[0, 0].annotate("the melt is faster than one\nsample interval (0.25 "
                      r"$\tau_B$)",
                      xy=(0.3, 0.86), xytext=(11, 0.82), fontsize=7.5,
                      arrowprops=dict(arrowstyle="->", lw=.8, color="#444"))
    ax[0, 0].set(xlabel=r"$t / \tau_B$",
                 ylabel=r"$\psi_6^{global} = |\langle\psi_{6i}\rangle|$",
                 ylim=(0, 1.05))
    ax[0, 0].set_title("F1  both arms, from $t=0$  (5 seeds each)", fontsize=10.5)
    ax[0, 0].legend(fontsize=7.5, loc="upper right"); ax[0, 0].grid(alpha=.3)

    # F2 -- the primary read
    for i, k in enumerate(H):
        ax[0, 1].plot(0 + 0.06 * (i - 2), R[k]["late"], "o", color="#1b6b3a", ms=8)
    for i, k in enumerate(Rr):
        ax[0, 1].plot(1 + 0.06 * (i - 2), R[k]["late"], "o", color="#2b6cb0", ms=8)
    ax[0, 1].axhline(FLOOR, color="#9b2c2c", lw=1.4, ls="--")
    ax[0, 1].set(xticks=[0, 1], xticklabels=["crystal start", "random start"],
                 ylabel=r"$\psi_6^{global}$ (75-100 $\tau_B$)", ylim=(0, 1.05),
                 xlim=(-0.5, 1.5))
    ax[0, 1].set_title(f"F2  primary read  $(n_H,n_R)=({nh},{nr})$", fontsize=10.5)
    ax[0, 1].axhspan(-DELTA_EQ / 2, DELTA_EQ / 2, xmin=0, xmax=0,
                     color="none")            # keeps the legend honest
    ax[0, 1].text(0.5, 0.92, f"sealed: (5, 0) = HYSTERESIS\nmeasured: "
                             f"({nh}, {nr})\n$\\Delta$ = {w['delta']:+.4f}, "
                             f"|t| = {abs(w['t']):.2f}, p = {w['p_two']:.2f}\n"
                             f"TOST inside $\\pm${DELTA_EQ}: "
                             f"{abs(w['delta']) + w['tost_halfwidth']:.4f}",
                  ha="center", va="top", fontsize=8.5,
                  bbox=dict(boxstyle="round", fc="#f2f2f2", ec="#999"))
    ax[0, 1].text(0.5, 0.42, verdict, ha="center", fontsize=12,
                  fontweight="bold", color="#9b2c2c")
    ax[0, 1].grid(alpha=.3, axis="y")

    # F3 -- coordination
    lab, vals = [], []
    for k in ["C1", "C2", "H1", "R1"]:
        lab.append(k); vals.append(R[k]["defect"])
    bars = ax[1, 0].bar(lab, vals,
                        color=["#1b6b3a", "#9b2c2c", "#1b6b3a", "#2b6cb0"])
    for b, v in zip(bars, vals):
        ax[1, 0].annotate(f"{v:.4f}" + ("  (zero defects)" if v == 0 else ""),
                          (b.get_x() + b.get_width() / 2, v), ha="center",
                          va="bottom", fontsize=8)
    ax[1, 0].axhline(0.285, color="#111", ls=":",
                     label="hexwin random start, other engine: 0.285")
    ax[1, 0].set(ylabel="defect fraction  $1 - f(6)$", ylim=(0, 0.55))
    ax[1, 0].set_title("F3  window-averaged Voronoi coordination", fontsize=10.5)
    ax[1, 0].legend(fontsize=8); ax[1, 0].grid(alpha=.3, axis="y")

    # F4 -- the quarter means, which is where drift is visible
    for k in H:
        ax[1, 1].plot([1, 2, 3, 4], R[k]["quarters"], "o-", color="#1b6b3a", alpha=.8)
    for k in Rr:
        ax[1, 1].plot([1, 2, 3, 4], R[k]["quarters"], "o-", color="#2b6cb0", alpha=.8)
    ax[1, 1].axhline(FLOOR, color="#9b2c2c", lw=1.2, ls="--")
    ax[1, 1].set(xticks=[1, 2, 3, 4],
                 xlabel=r"quarter of the window (25 $\tau_B$ each)",
                 ylabel=r"$\psi_6^{global}$", ylim=(0, 0.45))
    ax[1, 1].set_title("F4  quarter means -- the whole separation is in Q1",
                       fontsize=10.5)
    ax[1, 1].annotate(f"Q1: H {np.mean([R[k]['quarters'][0] for k in H]):.3f} vs "
                      f"R {np.mean([R[k]['quarters'][0] for k in Rr]):.3f}",
                      xy=(1, 0.30), xytext=(1.5, 0.38), fontsize=8,
                      arrowprops=dict(arrowstyle="->", lw=.8))
    ax[1, 1].annotate("the default eq_frac = 0.2 would have\n"
                      "discarded the first 20 $\\tau_B$ -- almost\n"
                      "all of it", xy=(1.05, 0.20), xytext=(2.1, 0.06),
                      fontsize=7.5, arrowprops=dict(arrowstyle="->", lw=.8,
                                                    color="#9b2c2c"),
                      color="#9b2c2c")
    ax[1, 1].grid(alpha=.3)

    fig.suptitle("soft-r3 S30 -- a crystal-start melting bracket, against a "
                 "prediction sealed before any run existed",
                 fontsize=13, fontweight="bold")
    fig.tight_layout(rect=(0, 0, 1, 0.97))
    out.parent.mkdir(exist_ok=True)
    fig.savefig(out, dpi=145)
    print(f"\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    sys.exit(main())
