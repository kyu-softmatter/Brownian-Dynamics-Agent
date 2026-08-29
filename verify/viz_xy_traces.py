import json
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import gsd.hoomd
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

n = 9
D = "runs/chain-bend-2d-dlvo__n9-w3000-a1470-kt100__d040a3311aab"
J = "runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr-kt100__759c7ac313b7"
TRAPPED = [0, 4, 8]      # 0·8 = 고정(힘센서), 4 = 구동

def load(d):
    t = gsd.hoomd.open(f"{d}/traj_A.gsd", "r")
    p = np.array([t[i].particles.position[:n, :2] for i in range(len(t))])
    s = json.load(open(f"{d}/spec.json"))
    nprod, dt = s["numerics"]["n_prod"], s["numerics"]["dt_star"]
    om = s["params"]["omega_star"]
    # 프레임 → 구동 주기 단위 시간
    tt = np.linspace(0, nprod * dt, len(p)) * om / (2 * np.pi)
    return p, tt

pD, tD = load(D)
pJ, tJ = load(J)

fig, ax = plt.subplots(2, 2, figsize=(15, 9), sharex="col")
cmap = plt.cm.viridis(np.linspace(0, 0.92, n))

for col, (p, tt, name) in enumerate([(pD, tD, "DLVO-only"), (pJ, tJ, "DLVO + JKR 굽힘")]):
    # ── y(t) ──
    a = ax[0, col]
    for i in range(n):
        lw, ls = (2.0, "-") if i in TRAPPED else (1.0, "-")
        a.plot(tt, p[:, i, 1], color=cmap[i], lw=lw, ls=ls,
               label=f"{i}{' (구동)' if i == 4 else ' (센서)' if i in TRAPPED else ''}")
    a.axhline(0, color="gray", lw=0.5)
    a.set_ylabel("y [d]")
    a.set_title(f"{name} — y 방향 (구동 방향)", fontsize=11)
    a.legend(fontsize=7, ncol=3, loc="upper right")
    a.grid(alpha=0.25)

    # ── x(t), 초기위치 대비 ──
    a = ax[1, col]
    dx = p[:, :, 0] - p[0, :, 0]
    for i in range(n):
        lw = 2.0 if i in TRAPPED else 1.0
        a.plot(tt, dx[:, i], color=cmap[i], lw=lw)
    a.axhline(0, color="gray", lw=0.5)
    a.set_xlabel("시간 [구동 주기]"); a.set_ylabel("Δx [d]  (초기위치 대비)")
    a.set_title(f"{name} — x 방향 (구동과 수직)\n"
                f"전범위 {np.ptp(dx):.4f} d", fontsize=11)
    a.grid(alpha=0.25)

# 같은 행은 같은 축으로 (직접 비교 가능하게)
for r in range(2):
    lo = min(ax[r, 0].get_ylim()[0], ax[r, 1].get_ylim()[0])
    hi = max(ax[r, 0].get_ylim()[1], ax[r, 1].get_ylim()[1])
    ax[r, 0].set_ylim(lo, hi); ax[r, 1].set_ylim(lo, hi)

fig.suptitle("입자별 x·y 시간 궤적 — n=9, ω=3000 rad/s, a=1470nm(=1d), k_t×100, seed=1\n"
             "굵은 선 = 트랩된 비드 (0·8 고정 센서, 4 구동) · 색 = 비드 번호(왼→오)",
             fontsize=12)
fig.tight_layout()
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/xy_traces.png"
fig.savefig(out, dpi=140); print("saved", out)

# 수치 요약
for p, name in [(pD, "DLVO"), (pJ, "JKR")]:
    dx = p[:, :, 0] - p[0, :, 0]
    print(f"{name}: y 진폭(비드4) {np.ptp(p[:,4,1])/2:.4f} d | "
          f"y 진폭(비드2) {np.ptp(p[:,2,1])/2:.4f} d | "
          f"Δx 전범위 {np.ptp(dx):.5f} d | Δx RMS {dx.std():.5f} d")
