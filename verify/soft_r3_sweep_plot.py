"""Phase 1-B sweep summary -- answers the sketch's question
("final configuration? rdf / voronoi / structure").

    $PY scratch/soft_r3_sweep_plot.py
    → runs/soft-r3-2d-A-sweep__SUMMARY/{summary.png, summary.md}
"""
import glob
import json
import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"   # keep mathtext on a font
                                                         # carrying minus and micro

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "runs" / "soft-r3-2d-A-sweep__SUMMARY"
OUT.mkdir(parents=True, exist_ok=True)
HEX_NN = math.sqrt(2 / math.sqrt(3))

runs = []
for tag in ("A0.1", "A1", "A10", "A100"):
    d = glob.glob(str(ROOT / f"runs/soft-r3-2d-A-sweep__{tag}__*"))
    assert len(d) == 1, (tag, d)
    m = json.load(open(d[0] + "/metrics.json"))
    z = np.load(d[0] + "/observables.npz")
    runs.append(dict(tag=tag, m=m, z=z,
                     A=m["dimensionless"]["A      = U(d)/kT        접촉 결합"],
                     G=m["structure"]["Gamma"]))

a_star = math.sqrt(math.pi / (4 * runs[0]["m"]["physical"]["phi"]))
colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(runs)))

fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(3, 4, height_ratios=[1.15, 1.0, 0.95], hspace=0.55, wspace=0.32)

# ── (1) g(r), all of it ───────────────────────────────────────────────
ax = fig.add_subplot(gs[0, :2])
for r_, c in zip(runs, colors):
    ax.plot(r_["z"]["rdf_r"], r_["z"]["rdf_g"], lw=1.5, color=c,
            label=f"A={r_['A']:g}  Γ={r_['G']:.2f}")
ax.axvline(HEX_NN * a_star, ls=":", c="k", alpha=.6)
ax.text(HEX_NN * a_star, ax.get_ylim()[1] * .93, " $a_{NN}$(hex)", fontsize=8)
ax.axhline(1, c="k", lw=.5, alpha=.4)
ax.set(xlabel="r / d", ylabel="g(r)", xlim=(0.8, 6.5),
       title="(1) radial distribution function -- raising the coupling Gamma makes "
             "the crystal peaks appear")
ax.legend(fontsize=8); ax.grid(alpha=.3)

# ── (2) psi6 and the 6-fold fraction vs Gamma ─────────────────────────
ax = fig.add_subplot(gs[0, 2:])
G = [r_["G"] for r_ in runs]
p6 = [r_["m"]["structure"]["psi6"] for r_ in runs]
p6e = [r_["m"]["structure"]["psi6_sem"] for r_ in runs]
c6 = [r_["m"]["structure"]["coord_hist"][6] for r_ in runs]
ax.errorbar(G, p6, yerr=p6e, marker="o", lw=1.6, capsize=3, label=r"$|\psi_6|$")
ax.plot(G, c6, marker="s", lw=1.6, ls="--", label="Voronoi 6-fold fraction")
ax.set(xscale="log", xlabel=r"$\Gamma = U(a_{mean})/k_BT$", ylabel="order parameter",
       ylim=(0, 1.05),
       title="(2) hexagonal order -- the transition lies between Gamma 3 and 30")
ax.axvspan(3, 30, alpha=.12, color="tab:red")
ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
for g_, y_ in zip(G, p6):
    ax.annotate(f"{y_:.3f}", (g_, y_), textcoords="offset points", xytext=(4, 6), fontsize=7)

# ── (3) four final configurations ─────────────────────────────────────
for i, (r_, c) in enumerate(zip(runs, colors)):
    ax = fig.add_subplot(gs[1, i])
    xy = r_["z"]["final_xy"]
    L = r_["m"]["dimensionless"]["L/d                     박스 크기"]
    ax.plot(xy[:, 0], xy[:, 1], "o", ms=3.0, color=c, mec="none")
    ax.set(xlim=(-L / 2, L / 2), ylim=(-L / 2, L / 2), aspect="equal",
           title=f"A={r_['A']:g}  Γ={r_['G']:.2f}", xticks=[], yticks=[])
    ax.set_xlabel(f"$\\psi_6$={r_['m']['structure']['psi6']:.3f}", fontsize=9)
fig.text(0.5, 0.655, "(3) final configurations (phi=0.35, N=400, T_obs=100 tau_B)",
         ha="center", fontsize=10.5)

# ── (4) Voronoi coordination distribution ─────────────────────────────
ax = fig.add_subplot(gs[2, :2])
w = 0.2
for j, (r_, c) in enumerate(zip(runs, colors)):
    ch = np.array(r_["m"]["structure"]["coord_hist"])
    ks = np.arange(3, 11)
    ax.bar(ks + (j - 1.5) * w, ch[3:11], width=w, color=c, label=f"Γ={r_['G']:.2f}")
ax.set(xlabel="Voronoi coordination", ylabel="fraction",
       title="(4) coordination distribution -- the crystal concentrates at 6")
ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")

# ── (5) convergence check ─────────────────────────────────────────────
ax = fig.add_subplot(gs[2, 2:])
base = json.load(open(glob.glob(str(ROOT / "runs/*A100__*/metrics.json"))[0]))
labels, dt_ch, rc_ch = [], [], []
for key, get in (("⟨U⟩/N", lambda mm: [x for x in mm["observables"]
                                       if x["name"].startswith("에너지")][0]["measured"]),
                 (r"$\psi_6$", lambda mm: mm["structure"]["psi6"]),
                 ("NN distance", lambda mm: mm["structure"]["nn_distance_d"]),
                 ("6-fold", lambda mm: mm["structure"]["coord_hist"][6])):
    labels.append(key)
    for tag, acc in (("A100-dt0.5", dt_ch), ("A100-rc7", rc_ch)):
        c = json.load(open(glob.glob(str(ROOT / f"runs/*{tag}*/metrics.json"))[0]))
        acc.append(100 * (get(c) - get(base)) / abs(get(base)))
x = np.arange(len(labels))
ax.bar(x - 0.19, dt_ch, 0.36, label="CV1  dt halved")
ax.bar(x + 0.19, rc_ch, 0.36, label=r"CV2  $r_c$ 5a→7a")
ax.axhline(0, c="k", lw=.6)
ax.set(xticks=x, xticklabels=labels, ylabel="change [%]",
       title="(5) convergence -- structure invariant, only the absolute energy "
             "depends on $r_c$")
ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")
for xi, v in zip(x - 0.19, dt_ch):
    ax.annotate(f"{v:+.3f}", (xi, v), textcoords="offset points", xytext=(0, 3 if v >= 0 else -11),
                ha="center", fontsize=7)
for xi, v in zip(x + 0.19, rc_ch):
    ax.annotate(f"{v:+.2f}", (xi, v), textcoords="offset points", xytext=(0, 3 if v >= 0 else -11),
                ha="center", fontsize=7)

fig.suptitle("Phase 1-B  soft-r3-2d-A-sweep  --  $U/k_BT = A(d/r)^3$ + WCA core, "
             "d=5 um silica / water 300 K, phi=0.35", fontsize=12.5)
fig.savefig(OUT / "summary.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# ── summary table ─────────────────────────────────────────────────────
lines = ["# Phase 1-B `soft-r3-2d-A-sweep` -- sweep results", "",
         "Physical system: `U/kT = A(d/r)^3` + WCA(sigma=d, eps=kT), 2D periodic, "
         "d=5 µm, water@300 K,",
         "φ=0.35 (a_mean=1.498 d), N=400, r_c=5 a_mean, T_obs=100 τ_B (τ_B=242.05 s)", "",
         "| A | Gamma=U(a_mean)/kT | psi6 | 6-fold | NN/d | sigma_NN/NN | state | "
         "dt/tau_B | wall-clock |",
         "|---|---|---|---|---|---|---|---|---|"]
for r_ in runs:
    st = r_["m"]["structure"]
    lines.append(f"| {r_['A']:g} | {r_['G']:.4f} | {st['psi6']:.4f} | "
                 f"{st['coord_hist'][6]:.3f} | {st['nn_distance_d']:.3f} | "
                 f"{st['nn_std_rel']:.4f} | "
                 f"{'hexagonal crystal' if st['psi6'] > 0.6 else 'fluid'} | "
                 f"{r_['m']['numerics']['dt_over_tau_B']:.2e} | "
                 f"{r_['m']['wall_seconds']/60:.0f} min |")
lines += ["", "## The answer to the sketch's 'final configuration?'", "",
          "- **Gamma <~ 3 is fluid, Gamma ~ 30 is a hexagonal crystal.** The "
          "transition lies between Gamma 3 and 30 (narrowing it further needs a "
          "denser A sweep).",
          "- **A=0.1 and A=1 are indistinguishable** (psi6 = 0.4347 for both). At "
          "phi=0.35 a four-decade A sweep",
          "  produces only 3 states -- the weak-coupling end is a WCA fluid where "
          "r^-3 is irrelevant.",
          "- The crystal's nearest-neighbour distance agrees with the hexagonal "
          "prediction `a_NN = sqrt(2/sqrt(3))*a_mean` to **+0.45%**.", "",
          "## Verification", "",
          "| check | result |", "|---|---|",
          "| direct two-particle comparison (`verify/verify_pair_table.py`) | 0.000% |",
          "| energy consistency `<U>/N` vs `(rho/2) int U g(r) 2 pi r dr` | "
          "+0.00 to +0.67% (7 runs) |",
          "| hexagonal NN distance vs a parameter-free prediction | +0.45% |",
          "| dilute limit `g(r)` vs `e^{-beta U}[1+rho int f f]` | "
          "RMS 2.43% (6.30% using 0th order alone) |",
          "| CV1 dt halved | `<U>/N` -0.004%, psi6 +0.13% |",
          "| CV2 `r_c` 5a->7a | structure invariant within 0.15%, absolute `<U>/N` "
          "**+7.5%** |", "",
          "## What was NOT checked", "",
          "- the dimensionless reading of `A` was inferred backwards from the "
          "sketch's goal (`system.yaml.not_verified`)",
          "- finite size: only N=400 was run. No N=900 comparison (CV3 incomplete)",
          "- no literature comparison (0 papers in the KB). The literature Gamma for "
          "the 2D r^-3 melting transition is unconfirmed",
          "- initial-condition dependence: only one RSA random placement. No "
          "comparison against a run started from a crystal",
          ]
(OUT / "summary.md").write_text("\n".join(lines))
print(f"→ {OUT.relative_to(ROOT)}/summary.png")
print(f"→ {OUT.relative_to(ROOT)}/summary.md")
