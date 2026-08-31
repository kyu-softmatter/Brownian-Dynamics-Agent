"""Chain-bending sweep results: figures plus SI conversion.

All figure text is English (CLAUDE.md).
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from simbot.units import (kT_si, stokes_drag_si, water_viscosity_si)

ROOT = Path("runs/chain-bend")
K_THETA = 1.0e4


def load(tag):
    p = ROOT / tag / "batch.json"
    if not p.exists():
        return []
    return [j for j in json.loads(p.read_text())["jobs"] if not j.get("skipped")]


st = load("explore-static")
th = load("explore-thermal")

S = {j["config"]["n_particles"]: j for j in st}
T = defaultdict(list)
for j in th:
    T[j["config"]["n_particles"]].append(j)

Ns = sorted(S)
Ls = np.array([n - 1.0 for n in Ns], float)
k_st = np.array([S[n]["kappa_star"] for n in Ns])
k_pr = np.array([S[n]["kappa_pred_star"] for n in Ns])
tau = np.array([S[n]["tau_relax_fit_star"] for n in Ns])
taub = 1.0 / k_pr

Nt = sorted(T)
Lt = np.array([n - 1.0 for n in Nt], float)
k_th = np.array([np.mean([x["kappa_star"] for x in T[n]]) for n in Nt])
e_th = np.array([np.std([x["kappa_star"] for x in T[n]], ddof=1) / np.sqrt(len(T[n]))
                 for n in Nt])
k_th_pr = np.array([T[n][0]["kappa_pred_star"] for n in Nt])

fig, ax = plt.subplots(1, 3, figsize=(15.5, 4.6))

# --- (a) kappa vs L -----------------------------------------------------
a = ax[0]
Lf = np.logspace(np.log10(Ls.min() * 0.85), np.log10(Ls.max() * 1.18), 60)
a.plot(Lf, 48 * K_THETA / Lf ** 3, "k-", lw=1.4,
       label=r"beam theory  $48\,k_\theta/L^3$")
a.plot(Ls, k_st, "o", ms=8, mfc="none", mew=1.8, color="#1f77b4",
       label="static (deterministic, $k_BT=0$)")
a.errorbar(Lt, k_th, yerr=e_th, fmt="s", ms=6, color="#d62728", capsize=3,
           label="thermal (4 seeds, $1/\\langle y_c^2\\rangle$)")
a.set_xscale("log"); a.set_yscale("log")
a.set_xlabel(r"span  $L^* = (N-1)\,\sigma$")
a.set_ylabel(r"point stiffness  $\kappa^*\;[k_BT/\sigma^2]$")
a.set_title("(a) chain bending stiffness vs length")
a.legend(fontsize=8.5, frameon=False)
a.grid(alpha=.25, which="both")

p3, a3 = np.polyfit(np.log(Ls[3:]), np.log(k_st[3:]), 1)
a.text(.04, .06, f"static fit ($N\\geq{Ns[3]}$):  $p = {p3:.3f}$",
       transform=a.transAxes, fontsize=9,
       bbox=dict(fc="white", ec="0.7", alpha=.9))

# --- (b) ratio ----------------------------------------------------------
b = ax[1]
b.axhline(1.0, color="k", lw=1.2)
b.axhspan(0.99, 1.01, color="k", alpha=.07, label=r"$\pm 1\,\%$")
b.plot(Ls, k_st / k_pr, "o-", ms=7, mfc="none", mew=1.8, color="#1f77b4",
       label="static")
b.errorbar(Lt, k_th / k_th_pr, yerr=e_th / k_th_pr, fmt="s-", ms=6,
           color="#d62728", capsize=3, label="thermal")
b.set_xscale("log")
b.set_xlabel(r"span  $L^*$")
b.set_ylabel(r"$\kappa_{\rm meas}/\kappa_{\rm pred}$")
b.set_title("(b) ratio to beam theory")
b.legend(fontsize=9, frameon=False, loc="lower right")
b.grid(alpha=.25, which="both")
b.set_ylim(0.84, 1.06)

# --- (c) relaxation time ------------------------------------------------
c = ax[2]
r = tau / taub
c.plot(Ls, r, "o", ms=8, mfc="none", mew=1.8, color="#2ca02c",
       label=r"measured  $\tau_{\rm relax}/\tau_\kappa$")
sl = np.polyfit(Ls, r, 1)[0]
c.plot(Lf, sl * Lf, "k--", lw=1.2, label=f"linear in $L$, slope {sl:.3f}")
# theory: tau_1 = gamma L^4 / (pi^4 EI),  EI = k_theta b,  tau_kappa = L^3/(48 k_theta)
c.plot(Lf, (48 / np.pi ** 4) * Lf, "-", color="#ff7f0e", lw=1.4,
       label=r"$48/\pi^4 \cdot L$  (fundamental mode)")
c.set_xlabel(r"span  $L^*$")
c.set_ylabel(r"$\tau_{\rm relax}\,/\,(\gamma/\kappa)$")
c.set_title(r"(c) relaxation time:  $\tau \propto L^4$, not $L^3$")
c.legend(fontsize=8.5, frameon=False)
c.grid(alpha=.25)

fig.tight_layout()
outd = ROOT / "figs"
outd.mkdir(parents=True, exist_ok=True)
fig.savefig(outd / "chain_bend_length_sweep.png", dpi=155)
print(f"-> {outd/'chain_bend_length_sweep.png'}")

# =========================================================================
print()
print("=" * 74)
print("measurement summary (reduced units)")
print("=" * 74)
print(f"{'N':>4} {'L*':>5} {'static kappa*':>14} {'ratio':>8} | "
      f"{'thermal kappa*':>18} {'ratio':>16}")
for i, n in enumerate(Ns):
    if n in T:
        j = Nt.index(n)
        t = f"{k_th[j]:>11.5g}+-{e_th[j]:<6.3g} {k_th[j]/k_th_pr[j]:>8.4f}+-{e_th[j]/k_th_pr[j]:<6.4f}"
    else:
        t = " " * 18 + "        —"
    print(f"{n:>4} {Ls[i]:>5.0f} {k_st[i]:>14.5g} {k_st[i]/k_pr[i]:>8.4f} | {t}")

print()
print(f"static exponent (N>={Ns[3]}):  p = {p3:.4f}   (textbook -3)")
print(f"static coefficient:       exp(a)/(48 k_theta) = {np.exp(a3)/(48*K_THETA):.4f}")
print(f"tau_relax/tau_kappa slope: {sl:.4f}   "
      f"(fundamental mode 48/pi^4 = {48/np.pi**4:.4f})")

# =========================================================================
print()
print("=" * 74)
print("SI conversion -- what the papers' values give directly  "
      "[cited, not reproduced]")
print("=" * 74)
T_K = 298.15
eta, _ = water_viscosity_si(T_K)
A = 1.47e-6 / 2
SIG = 2 * A
kT = kT_si(T_K)
gam = stokes_drag_si(eta, A)
E_PMMA, A_C = 3.1e9, 40e-9              # Pa, m  [cited, not reproduced] PRL
EI = np.pi * E_PMMA * A_C ** 4 / 4.0     # the papers' circular-contact EI
KAPPA0 = 64e-3                           # N/m @10 mM  [cited, not reproduced]
print(f"  EI = pi E a_c^4/4 = {EI:.4e} N*m^2     (E={E_PMMA/1e9:g} GPa, a_c={A_C*1e9:g} nm)")
print(f"  k_theta = EI/sigma = {EI/SIG:.4e} J = {EI/SIG/kT:.4e} kT/rad^2")
print()
print(f"{'N':>4} {'L[um]':>8} {'48EI/L^3 [pN/um]':>18} "
      f"{'kappa_0(a/s)^3 [pN/um]':>24} {'ratio':>8}")
for n in (7, 11, 15, 21, 31, 41):
    L = (n - 1) * SIG
    beam = 48 * EI / L ** 3
    paper = KAPPA0 / (2.0 * (n - 1)) ** 3          # the papers' form (a = sigma/2)
    print(f"{n:>4} {L*1e6:>8.2f} {beam*1e6:>18.3f} {paper*1e6:>24.3f} {beam/paper:>8.2f}")
print()
print("  ★ The two differ by a constant factor. The papers' kappa_0 coefficient")
print("    cannot be pinned down from the PDF (two-column typesetting corrupts")
print("    coefficients, env-log stage 3). Exponent -3 agrees; coefficient open.")
