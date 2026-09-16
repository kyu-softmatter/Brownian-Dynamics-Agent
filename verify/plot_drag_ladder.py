"""Two panels: the per-rung error the plan assumes, and what it does to h0."""
from __future__ import annotations
import json, sys, warnings
from pathlib import Path
import numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
warnings.filterwarnings("error", message=".*missing from font.*")

ROOT = Path(__file__).resolve().parent.parent
R = json.loads((ROOT / "verify/_out/drag_ladder.json").read_text())
S = json.loads((ROOT / "verify/_out/drag_ladder_fit_sweep.json").read_text())

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12.6, 4.9))

# ---- A: per-rung sigma on gamma -------------------------------------------
hs = sorted(float(k) for k in R["rungs"])
sg = [R["rungs"][str(h)]["sigma_gamma_pct"] for h in hs]
sa = [R["rungs"][str(h)]["sigma_alpha_pct"] for h in hs]
ss = [R["rungs"][str(h)]["sigma_slope_pct"] for h in hs]
bg = [R["rungs"][str(h)]["bias_gamma_pct"] for h in hs]
ax1.plot(hs, sg, "o-", lw=2, ms=8, color="#1f77b4", label=r"measured  $\sigma_\gamma$ (camera on)")
ax1.plot(hs, sa, "s--", lw=1.4, ms=6, color="#7f7f7f",
         label=r"from $\alpha=k_BT/\mathrm{var}(x)$ alone")
ax1.plot(hs, ss, "^--", lw=1.4, ms=6, color="#2ca02c", label=r"from the drag slope alone")
ax1.axhline(3.0, color="#d62728", ls="--", lw=2,
            label="assumed in the plan:  3 %")
ax1.plot(hs, bg, "v:", lw=1.6, ms=7, color="#9467bd",
         label=r"camera BIAS on $\gamma$ (not scatter)")
ax1.axhline(0, color="k", lw=.7)
ax1.set_xlabel(r"rung height  $h$  [$\mu$m]")
ax1.set_ylabel(r"per-rung error on $\gamma$  [%]")
ax1.set_title(r"A · The per-rung error is 7–9.5 %, not 3 % — and it is $\alpha$, not the slope")
ax1.legend(fontsize=8, loc="center right"); ax1.grid(alpha=.3)
ax1.set_ylim(-9, 13)

# ---- B: sigma_h0 against the two stated thresholds -------------------------
h0s = sorted(float(k) for k in S)
sd = [S[str(k)]["sd"] for k in h0s]
sde = [S[str(k)]["sd_err"] for k in h0s]
x = np.arange(len(h0s))
ax2.bar(x, sd, yerr=sde, width=.5, color="#1f77b4", capsize=5,
        label=r"measured  $\sigma_{h_0}$  (n = 192 ladders each)")
ax2.axhline(0.195, color="#d62728", ls="--", lw=2)
ax2.text(len(h0s) - .45, 0.205, "claimed  ±0.195 µm", color="#d62728", fontsize=9, ha="right")
ax2.axhline(0.40, color="#ff7f0e", ls="-.", lw=2)
ax2.text(len(h0s) - .45, 0.415,
         "r1's own falsifier:  \"scatter exceeds ~0.4 µm → the ladder\ndoes not produce a trapping height\"",
         color="#ff7f0e", fontsize=8, ha="right", va="bottom")
ax2.set_xticks(x); ax2.set_xticklabels([f"{v:g}" for v in h0s])
ax2.set_xlabel(r"injected true offset  $h_0$  [$\mu$m]")
ax2.set_ylabel(r"$\sigma_{h_0}$  recovered by the six-rung fit  [$\mu$m]")
ax2.set_title(r"B · The six-rung fit returns $h_0$ to 0.43 µm — 2.2× the target")
ax2.legend(fontsize=8.5, loc="upper left"); ax2.grid(alpha=.3, axis="y")
ax2.set_ylim(0, 0.62)

fig.suptitle("Drag-slope ladder — per-rung precision measured, propagated into $h_0$   "
             "(exact OU + camera layer; the wall is synthetic and not simulated)", fontsize=10.5)
fig.tight_layout(rect=(0, 0, 1, 0.94))
out = ROOT / "verify/_out/drag_ladder.png"
fig.savefig(out, dpi=150)
print(f"wrote {out.relative_to(ROOT)}")
