"""Single-particle trajectories, densely sampled -- to see the run-and-tumble
structure by eye.

The production run (cases/abp_rod_2d.py) samples every 71 ms while tau_eff = 343 ms,
so a run segment is captured by only 5 samples. Here the sample interval is set to dt
so the straight runs and the abrupt tumbles are both drawn. This is **for looking**,
not for statistics (the observables come from the production run).

    $PY scratch/abp_traj_zoom.py
"""
import math, sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "cases"))
from bdbot import Q, physical as PH
import abp_rod_2d as CB

s = PH.load(ROOT / "intake/abp-rod-2d-run-flip")
lg = CB.build_ledger(s); D = lg.derived
f_ = CB.f_
Pe = f_(D["v"] * D["d"] / D["D_bar"])
dt_star = f_(CB.node(s, "numerics", "dt").to("s") / D["tau_B"])
tau_B = f_(D["tau_B"] / Q(1, "s")); d_um = f_(D["d"] / Q(1, "um"))
te = f_(D["tau_eff"] / Q(1, "s"))

N, n_prod = 6, 8000          # sample interval = dt (every step) -> 8000 samples = 12.7 s
res = CB.simulate(N, Pe, f_(D["D_r"] * D["tau_B"]), f_(D["tau_tumble"] / D["tau_B"]),
                  dt_star, n_eq=500, n_prod=n_prod, sample_every=1,
                  L_star=f_(D["L"] / D["d"]), seed=20260804, gsd_path=None,
                  progress=False)
print(f"  {res['xy'].shape[0]} samples, interval {dt_star*tau_B*1e3:.2f} ms "
      f"(tau_eff = {te*1e3:.0f} ms), {res['n_tumbles']} tumbles")

import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
fig, ax = plt.subplots(1, 2, figsize=(13, 5.6))
xy = res["xy"] * d_um; th = res["theta"]
t = np.arange(len(xy)) * dt_star * tau_B

cmap = plt.cm.turbo(np.linspace(.08, .92, N))
for k in range(N):
    p = xy[:, k, :] - xy[0, k, :]
    ax[0].plot(p[:, 0], p[:, 1], "-", lw=1.1, color=cmap[k], alpha=.9)
    d_th = np.abs((np.diff(th[:, k]) + math.pi) % (2 * math.pi) - math.pi)
    j = np.where(d_th > 1.5)[0]
    ax[0].plot(p[j + 1, 0], p[j + 1, 1], "o", ms=5, mfc="none", mec=cmap[k], mew=1.4)
ax[0].plot(0, 0, "k+", ms=10, mew=1.8)
ax[0].set(xlabel="x [µm]", ylabel="y [µm]", aspect="equal",
          title=f"{N} single-particle trajectories -- sample interval "
                f"{dt_star*tau_B*1e3:.1f} ms (circles = tumbles)")
ax[0].grid(alpha=.25)

# Zoom on one particle, together with its orientation angle
k = 0
p = xy[:, k, :] - xy[0, k, :]
sc = ax[1].scatter(p[:, 0], p[:, 1], c=t, s=5, cmap="viridis")
d_th = np.abs((np.diff(th[:, k]) + math.pi) % (2 * math.pi) - math.pi)
j = np.where(d_th > 1.5)[0]
ax[1].plot(p[j + 1, 0], p[j + 1, 1], "o", ms=9, mfc="none", mec="crimson", mew=1.8,
           label=f"{len(j)} tumbles")
ax[1].set(xlabel="x [µm]", ylabel="y [µm]", aspect="equal",
          title=f"one particle, zoomed -- straight runs + abrupt reorientation")
ax[1].legend(fontsize=9); ax[1].grid(alpha=.25)
plt.colorbar(sc, ax=ax[1], label="t [s]")
fig.suptitle(f"abp-rod-2d trajectories (2µm x 500nm ellipsoid, Pe={Pe:.2f}, "
             f"ℓ_p={f_(D['l_p']/Q(1,'um')):.2f} µm, τ_eff={te*1e3:.0f} ms)", fontsize=12)
fig.tight_layout()
out = ROOT / "runs/abp-rod-2d-run-flip__74fb2d81066a/trajectory_zoom.png"
fig.savefig(out, dpi=140); print(f"→ {out.relative_to(ROOT)}")
