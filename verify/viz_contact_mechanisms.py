"""Four-way comparison figure of contact mechanisms.

Which of them produces a bending stiffness and which cannot.

Panel (a) md.pair.friction is identically 0 at the DLVO secondary minimum, and even
          in contact it is exempt from rolling
      (b) ★ the discriminator -- amplitude dependence of the stiffness: DLVO tension
          gives K proportional to delta^2, real bending is delta-independent
      (c) the exact equivalence between rolling resistance and a harmonic bend
          (kappa_theta,eff = 0.5*k_r*R^2), and the orientation-freezing effect
      (d) a kT=0 omega sweep -- the frequency at which the two models diverge

⚠ All labels are in English (CLAUDE.md -- a Hangul matplotlib font has no
  minus sign or y-hat glyph).

    $PY scratch/viz_contact_mechanisms.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                        # noqa: E402
import numpy as np                                                     # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scratch"))
sys.path.insert(0, str(ROOT / "cases"))

from dlvo_tension_stiffening import (ELL, K_BOND, KAPPA_THETA,          # noqa: E402
                                     k_tension_pred, static_3point)
from rolling_contact import RollingContact, k_roll_from_kappa_theta     # noqa: E402
from scipy.optimize import minimize                                     # noqa: E402

R_C = 0.5
K_ROLL = k_roll_from_kappa_theta(KAPPA_THETA, R_C)
H_MIN_STAR = 0.00759259035993831
OUT = ROOT / "runs" / "_bending_model_compare"

C_BEND, C_ROLL, C_DLVO, C_REF = "#1f77b4", "#d62728", "#2ca02c", "#888888"
fig, axes = plt.subplots(2, 2, figsize=(13.5, 10.6))
fig.suptitle("What can and cannot produce chain bending stiffness from contact mechanics",
             fontsize=13, fontweight="bold")

# ══ (a) pair.friction ═════════════════════════════════════════════════════
ax = axes[0, 0]
SIG = 2 ** (-1 / 6)
r = np.linspace(0.86, 1.06, 600)
w = np.where(r < 1.0, 24 * 1.0 / r * (2 * (SIG / r) ** 12 - (SIG / r) ** 6), 0.0)
ax.plot(r, np.maximum(w, 0), color=C_DLVO, lw=2,
        label=r"friction weight $w(r)=-dU_{WCA}/dr$")
ax.set_ylim(-2.5, 72)
ax.axvline(1.0, color="k", ls="--", lw=1.2)
ax.text(0.996, 46, "WCA cutoff\n$r^*=1$", ha="right", fontsize=8.5)
ax.axvline(1 + H_MIN_STAR, color=C_ROLL, lw=2.2)
ax.text(1.011, 38, "DLVO secondary\nminimum\n$r^*=1.0076$\n"
        r"$\Rightarrow w \equiv 0$", color=C_ROLL, fontsize=8.5, va="top")
ax.set_xlabel(r"center-to-center distance $r^*$  [$d$]")
ax.set_ylabel(r"friction weight  [$k_BT/d$]")
ax.set_title("(a) md.pair.friction is inert at the DLVO bond", fontsize=10.5)
ax.legend(fontsize=8.5, loc="upper right")
axi = ax.inset_axes([0.13, 0.30, 0.36, 0.34])
meas = {"LJLinear": (6.1926, 0.0), "LJCoulomb": (18.578, 0.0), "CoulombNewton": (1.0, 0.0)}
xs = np.arange(3)
axi.bar(xs - 0.19, [v[0] for v in meas.values()], 0.36, color="#777", label="sliding")
axi.bar(xs + 0.19, [v[1] for v in meas.values()], 0.36, color=C_ROLL, label="rolling")
axi.set_xticks(xs); axi.set_xticklabels(["Lin", "Coul", "CN"], fontsize=7)
axi.set_ylabel(r"$F_{tan}$", fontsize=7.5); axi.tick_params(labelsize=7)
axi.set_title("in contact ($r^*$=0.95)", fontsize=7.5)
axi.legend(fontsize=6.5, loc="upper left")

# ══ (b) ★ amplitude discriminator ═════════════════════════════════════════
ax = axes[0, 1]
N = 9
deltas = np.geomspace(1e-3, 1.0, 11)
k_none = np.array([static_3point(N, d, bend="none") for d in deltas])
k_harm = np.array([static_3point(N, d, bend="harmonic") for d in deltas])
k_roll = np.array([static_3point(N, d, bend="rolling") for d in deltas])
ax.loglog(deltas, k_none, "o-", color=C_DLVO, lw=2, ms=5,
          label="DLVO only (central force)")
ax.loglog(deltas, k_tension_pred(N, deltas), "--", color="k", lw=1.3,
          label=r"prediction  $4k_{bond}\delta^2/((n{-}1)^3\ell^2)$")
ax.loglog(deltas, k_harm, "s-", color=C_BEND, lw=2, ms=5, label="JKR harmonic bending")
ax.loglog(deltas, k_roll, "^--", color=C_ROLL, lw=1.6, ms=5, label="rolling resistance")
ax.annotate(r"slope $=1.9999$" "\n" r"$\Rightarrow$ linear stiffness is 0",
            xy=(0.03, k_tension_pred(N, 0.03)), xytext=(0.0016, 3e2),
            arrowprops=dict(arrowstyle="->", color=C_DLVO), color=C_DLVO, fontsize=9)
ax.annotate("flat: true 2nd-order potential", xy=(0.02, k_harm[3]),
            xytext=(0.0016, 2.6e4), color=C_BEND, fontsize=9,
            arrowprops=dict(arrowstyle="->", color=C_BEND))
ax.set_xlabel(r"imposed mid-bead displacement  $\delta$  [$d$]")
ax.set_ylabel(r"3-point stiffness  $k=2U/\delta^2$  [$k_BT/d^2$]")
ax.set_title(r"(b) The discriminator is $\delta$-scaling, not magnitude  [n=9]", fontsize=10.5)
ax.legend(fontsize=8, loc="lower right")
ax.grid(alpha=0.25, which="both")

# ══ (c) rolling ↔ harmonic equivalence ════════════════════════════════════
ax = axes[1, 0]


def three_point_rolling(n, delta, relax):
    pos0 = np.zeros((n, 3)); pos0[:, 0] = np.arange(n) * ELL
    quat0 = np.tile(np.array([1.0, 0, 0, 0]), (n, 1))
    rc = RollingContact([[i, i + 1] for i in range(n - 1)], pos0, quat0, R_C, K_ROLL, 0.0)
    c = n // 2
    fy = [i for i in range(n) if i not in (0, c, n - 1)]

    def unpack(v):
        y = np.zeros(n); y[c] = delta; y[fy] = v[:len(fy)]
        th = v[len(fy):] if relax else np.zeros(n)
        p = pos0.copy(); p[:, 1] = y
        cq, sq = np.cos(th / 2), np.sin(th / 2)
        return p, np.stack([cq, np.zeros(n), np.zeros(n), sq], axis=1)

    def fj(v):
        p, q = unpack(v)
        F, T = rc.force_torque(p, q)
        g = -F[fy, 1]
        if relax:
            g = np.concatenate([g, -T[:, 2]])
        return rc.energy(p, q), g

    nv = len(fy) + (n if relax else 0)
    res = minimize(fj, np.zeros(nv), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=100000, maxfun=100000, ftol=0.0, gtol=0.0))
    return 2 * res.fun / delta ** 2


def harm_3p(n, kt):
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ELL
    A = kt * (B.T @ B)
    c = n // 2
    fx, fr = [0, c, n - 1], [i for i in range(n) if i not in (0, c, n - 1)]
    yf = np.array([0.0, 1.0, 0.0]); y = np.zeros(n); y[fx] = yf
    y[fr] = np.linalg.solve(A[np.ix_(fr, fr)], -A[np.ix_(fr, fx)] @ yf)
    return float(y @ A @ y)


ns = [5, 7, 9, 11, 15]
D = 1e-4
kr_free = np.array([three_point_rolling(n, D, True) for n in ns])
kr_froz = np.array([three_point_rolling(n, D, False) for n in ns])
kh = np.array([harm_3p(n, KAPPA_THETA) for n in ns])
ax.semilogy(ns, kh, "s-", color=C_BEND, lw=2, ms=7, label="harmonic bending  $\\kappa_\\theta$")
ax.semilogy(ns, kr_free, "o--", color=C_ROLL, lw=1.8, ms=8, mfc="none",
            label=r"rolling, orientations relaxed  ($\kappa_{\theta,eff}=\frac{1}{2}k_rR^2$)")
ax.semilogy(ns, kr_froz, "^:", color="#ff7f0e", lw=1.8, ms=7,
            label="rolling, orientations frozen")
for n, a, b in zip(ns, kr_free, kr_froz):
    ax.annotate(f"{b/a:.0f}x", (n, b), textcoords="offset points", xytext=(4, 4),
                fontsize=8, color="#ff7f0e")
ax.set_xlabel("chain length  $n$  [beads]")
ax.set_ylabel(r"3-point stiffness  [$k_BT/d^2$]")
ax.set_title("(c) Rolling resistance == harmonic bending (if orientations relax)",
             fontsize=10.5)
ax.legend(fontsize=8, loc="upper right")
ax.grid(alpha=0.25, which="both")
ax.text(0.03, 0.06, f"max deviation (relaxed vs harmonic): "
        f"{100*np.abs(kr_free/kh-1).max():.4f}%", transform=ax.transAxes,
        fontsize=8.5, color=C_ROLL,
        bbox=dict(fc="white", ec=C_ROLL, alpha=0.85, boxstyle="round,pad=0.3"))

# ══ (d) ω sweep ═══════════════════════════════════════════════════════════
ax = axes[1, 1]
TAU_ROT = (4 * R_C ** 2 / 3) / (K_ROLL * R_C ** 2)
sw = OUT / "sweep.json"
if sw.exists():
    rows = json.loads(sw.read_text())
    rows = [r for r in rows if "K_static" in r]
else:
    rows = []
if rows:
    om = np.array([r["omega"] for r in rows])
    Kh = np.array([r["harmonic"]["K_re"] for r in rows])
    Kr = np.array([r["rolling"]["K_re"] for r in rows])
    Ks = rows[0]["K_static"]
    ax.axhline(Ks, color=C_REF, ls="--", lw=1.3,
               label=f"static limit (exact) = {Ks:.4g}")
    ax.loglog(om, Kh, "s-", color=C_BEND, lw=2, ms=6, label="harmonic bending")
    ax.loglog(om, Kr, "o--", color=C_ROLL, lw=1.8, ms=7, mfc="none",
              label="rolling resistance")
    ax.axvline(1 / TAU_ROT, color="#ff7f0e", lw=1.5, ls=":")
    ax.annotate(r"$\omega_c=1/\tau_{rot}$" "\n" f"= {1/TAU_ROT:.2e}",
                xy=(1 / TAU_ROT, 0.30), xycoords=("data", "axes fraction"),
                xytext=(6, 0), textcoords="offset points",
                color="#ff7f0e", fontsize=8.5, va="center")
    ax.axvline(18453, color="k", lw=1.5, ls="-.")
    ax.text(18453 * 1.15, Ks * 1.6, "production\n$\\omega^*$ = 18453", fontsize=8.5)
    # The ratio of the two models is only meaningful where bending dominates --
    # above that it is radial-bond / drag saturation
    bend_dom = om[Kh / Ks <= 2.0]
    if len(bend_dom):
        ax.axvspan(om.min() * 0.5, bend_dom.max() * 1.7, color="#2ca02c", alpha=0.07)
        ax.text(om.min() * 0.7, Kh.max() * 0.55, "bending-dominated\n($K'/K'_{static}\\leq2$)",
                fontsize=8, color="#2ca02c")
    ax.set_xlabel(r"drive frequency  $\omega^*$  [$1/\tau$]")
    ax.set_ylabel(r"$K'$  [$k_BT/d^2$]")
    ax.legend(fontsize=8, loc="lower right")
    ax.grid(alpha=0.25, which="both")
    axr = ax.inset_axes([0.13, 0.63, 0.40, 0.32])
    axr.semilogx(om, Kr / Kh, "o-", color=C_ROLL, ms=4, lw=1.5)
    axr.axhline(1.0, color="k", lw=0.9, ls="--")
    axr.axvline(18453, color="k", lw=1.0, ls="-.")
    axr.set_ylabel("rolling/harmonic", fontsize=7); axr.tick_params(labelsize=6.5)
    axr.set_xlim(ax.get_xlim()); axr.grid(alpha=0.25)
else:
    ax.text(0.5, 0.5, "sweep.json not ready", ha="center", va="center",
            transform=ax.transAxes, fontsize=11, color="#999")
ax.set_title("(d) kT=0 deterministic sweep: where the two models diverge", fontsize=10.5)

fig.tight_layout(rect=(0, 0, 1, 0.965))
fig.subplots_adjust(hspace=0.30, wspace=0.24)
OUT.mkdir(parents=True, exist_ok=True)
p = OUT / "contact_mechanisms.png"
fig.savefig(p, dpi=155)
print(f"→ {p}")
