"""Two panels: what the exposure formula gets wrong, and what f_c's error really is.

    $PY verify/plot_detector_and_fc_bias.py
"""
from __future__ import annotations
import json, math, sys, warnings
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "verify"))
from detector_model import (apply_detector, blur_deflation_exact,          # noqa: E402
                            blur_deflation_free_particle, ou_exact)

warnings.filterwarnings("error", message=".*missing from font.*")

U_THREAD = 0.027881          # t_exp/tau at t_exp = 0.45 ms, tau = 16.14 ms
res = json.loads((ROOT / "verify/_out/fc_estimator_bias.json").read_text())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.5, 4.8))

# ---- panel A: blur coefficient -------------------------------------------
uu = np.logspace(-2.2, 0.7, 200)
ax1.plot(uu, 100 * np.array([blur_deflation_exact(u) for u in uu]), "-", lw=2,
         color="#1f77b4", label=r"exact OU:  $1-2(u-1+e^{-u})/u^2$")
ax1.plot(uu, 100 * np.array([blur_deflation_free_particle(u) for u in uu]), "--", lw=2,
         color="#d62728", label=r"free-particle MSD term:  $2u/3$  (AM r1)")

rng = np.random.default_rng(4711)
mu, mm, me = [], [], []
for u in (0.027881, 0.1, 0.3, 1.0, 3.0):
    n_exp, nb = 32, 6
    dt = u / n_exp
    n_frame = max(n_exp, int(round(5.0 / dt)))
    n = n_frame * 400
    r = []
    for _ in range(nb):
        x = ou_exact(n, dt, 1.0, rng, 16)
        blurred, _ = apply_detector(x, dt, u, n_frame * dt, 0.0, rng)
        inst = x[:(len(x)//n_frame)*n_frame].reshape(-1, n_frame, x.shape[1])[:, 0, :]
        r.append(float((blurred**2).mean() / (inst**2).mean()))
    mu.append(u); mm.append(100*(1-np.mean(r)))
    me.append(100*np.std(r, ddof=1)/math.sqrt(nb))
ax1.errorbar(mu, mm, yerr=me, fmt="o", ms=7, color="k", zorder=5, capsize=3,
             label="measured (exact-OU + boxcar)")
ax1.axvline(U_THREAD, color="grey", ls=":", lw=1.2)
ax1.annotate(f"this thread\n$t_{{exp}}$=0.45 ms, $\\tau$=16.14 ms\n"
             f"exact 0.92 %  vs  AM 1.86 %",
             xy=(U_THREAD, 1.0), xytext=(0.012, 12), fontsize=8.5,
             arrowprops=dict(arrowstyle="->", lw=1, color="grey"))
ax1.set_xscale("log"); ax1.set_yscale("log")
ax1.set_xlabel(r"$u = t_{exp}/\tau_k$")
ax1.set_ylabel(r"deflation of $\mathrm{var}(x)$ by exposure  [%]")
ax1.set_title("A · The blur coefficient for a TRAPPED bead is $u/3$, not $2u/3$")
ax1.legend(fontsize=8.5, loc="upper left"); ax1.grid(alpha=.3, which="both")

# ---- panel B: f_c estimator bias -----------------------------------------
spt = sorted(int(k) for k in res["vs_sampling"])
b = [res["vs_sampling"][str(k)][0] for k in spt]
e = [res["vs_sampling"][str(k)][1] for k in spt]
ax2.errorbar(spt, b, yerr=e, fmt="o-", ms=7, lw=1.8, color="#1f77b4", capsize=3,
             label="estimator bias on EXACT OU (answer known)")
ax2.axhline(res["run_err_pct"], color="#d62728", ls="--", lw=2,
            label=f"archived run reported  {res['run_err_pct']:+.2f} %")
ax2.axhline(0, color="k", lw=.8)
ax2.axhspan(-5, 5, color="green", alpha=.06)
ax2.text(2.2, 4.3, "run's 5 % tolerance band", fontsize=8, color="green")
ax2.axvline(10, color="grey", ls=":", lw=1.2)
ax2.annotate("BD convention\nand the run", xy=(10, b[spt.index(10)]),
             xytext=(11, 5.6), fontsize=8.5,
             arrowprops=dict(arrowstyle="->", lw=1, color="grey"))
ax2.set_xscale("log")
ax2.set_xticks(spt); ax2.set_xticklabels([str(s) for s in spt])
ax2.minorticks_off()
ax2.set_xlabel(r"samples per $\tau_k$")
ax2.set_ylabel(r"bias on fitted $f_c$  [%]")
ax2.set_title(r"B · The run's $+1.17\,\%$ on $f_c$ is the estimator, not the physics")
ax2.legend(fontsize=8.5, loc="upper center"); ax2.grid(alpha=.3, which="both")
ax2.set_ylim(-1, 9)

fig.suptitle("Detector layer and Lorentzian-fit bias — synthetic OU, answers known in "
             "closed form   (verify/verify_*.py)", fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.95))
out = ROOT / "verify/_out/detector_and_fc_bias.png"
fig.savefig(out, dpi=150)
print(f"wrote {out.relative_to(ROOT)}")
