"""단일 입자 궤적을 촘촘히 — run-and-tumble 구조를 눈으로 확인.

본 런(cases/abp_rod_2d.py)의 표본 간격은 71 ms 인데 τ_eff = 343 ms 라 런 구간이
5 표본밖에 안 잡힌다. 여기서는 표본 간격을 dt 로 두어 직선 런 + 급격한 텀블을 그린다.
통계용이 아니라 **보기용**이다 (관측량은 본 런에서 나온다).

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

N, n_prod = 6, 8000          # 표본 간격 = dt (매 스텝) → 8000 표본 = 12.7 s
res = CB.simulate(N, Pe, f_(D["D_r"] * D["tau_B"]), f_(D["tau_tumble"] / D["tau_B"]),
                  dt_star, n_eq=500, n_prod=n_prod, sample_every=1,
                  L_star=f_(D["L"] / D["d"]), seed=20260804, gsd_path=None,
                  progress=False)
print(f"  표본 {res['xy'].shape[0]}개 · 간격 {dt_star*tau_B*1e3:.2f} ms "
      f"(τ_eff = {te*1e3:.0f} ms) · 텀블 {res['n_tumbles']}회")

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
          title=f"단일 입자 궤적 {N}개 — 표본 간격 {dt_star*tau_B*1e3:.1f} ms "
                f"(원 = 텀블)")
ax[0].grid(alpha=.25)

# 한 입자를 확대: 방향각과 함께
k = 0
p = xy[:, k, :] - xy[0, k, :]
sc = ax[1].scatter(p[:, 0], p[:, 1], c=t, s=5, cmap="viridis")
d_th = np.abs((np.diff(th[:, k]) + math.pi) % (2 * math.pi) - math.pi)
j = np.where(d_th > 1.5)[0]
ax[1].plot(p[j + 1, 0], p[j + 1, 1], "o", ms=9, mfc="none", mec="crimson", mew=1.8,
           label=f"텀블 {len(j)}회")
ax[1].set(xlabel="x [µm]", ylabel="y [µm]", aspect="equal",
          title=f"입자 1개 확대 — 직선 런 + 급격한 방향 전환")
ax[1].legend(fontsize=9); ax[1].grid(alpha=.25)
plt.colorbar(sc, ax=ax[1], label="t [s]")
fig.suptitle(f"abp-rod-2d 궤적 (2µm × 500nm 타원체, Pe={Pe:.2f}, "
             f"ℓ_p={f_(D['l_p']/Q(1,'um')):.2f} µm, τ_eff={te*1e3:.0f} ms)", fontsize=12)
fig.tight_layout()
out = ROOT / "runs/abp-rod-2d-run-flip__74fb2d81066a/trajectory_zoom.png"
fig.savefig(out, dpi=140); print(f"→ {out.relative_to(ROOT)}")
