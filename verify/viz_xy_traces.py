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
TRAPPED = [0, 4, 8]      # 0 and 8 = fixed (force sensors), 4 = driven

def load(d):
    t = gsd.hoomd.open(f"{d}/traj_A.gsd", "r")
    p = np.array([t[i].particles.position[:n, :2] for i in range(len(t))])
    s = json.load(open(f"{d}/spec.json"))
    nprod, dt = s["numerics"]["n_prod"], s["numerics"]["dt_star"]
    om = s["params"]["omega_star"]
    # frame -> time in units of the drive period
    tt = np.linspace(0, nprod * dt, len(p)) * om / (2 * np.pi)
    return p, tt

pD, tD = load(D)
pJ, tJ = load(J)

fig, ax = plt.subplots(2, 2, figsize=(15, 9), sharex="col")
cmap = plt.cm.viridis(np.linspace(0, 0.92, n))

for col, (p, tt, name) in enumerate([(pD, tD, "DLVO-only"),
                                     (pJ, tJ, "DLVO + JKR bending")]):
    # ── y(t) ──
    a = ax[0, col]
    for i in range(n):
        lw, ls = (2.0, "-") if i in TRAPPED else (1.0, "-")
        a.plot(tt, p[:, i, 1], color=cmap[i], lw=lw, ls=ls,
               label=f"{i}{' (driven)' if i == 4 else ' (sensor)' if i in TRAPPED else ''}")
    a.axhline(0, color="gray", lw=0.5)
    a.set_ylabel("y [d]")
    a.set_title(f"{name} -- y direction (the drive direction)", fontsize=11)
    a.legend(fontsize=7, ncol=3, loc="upper right")
    a.grid(alpha=0.25)

    # ── x(t), relative to the initial position ──
    a = ax[1, col]
    dx = p[:, :, 0] - p[0, :, 0]
    for i in range(n):
        lw = 2.0 if i in TRAPPED else 1.0
        a.plot(tt, dx[:, i], color=cmap[i], lw=lw)
    a.axhline(0, color="gray", lw=0.5)
    a.set_xlabel("time [drive periods]"); a.set_ylabel("dx [d]  (from initial position)")
    a.set_title(f"{name} -- x direction (perpendicular to the drive)\n"
                f"full range {np.ptp(dx):.4f} d", fontsize=11)
    a.grid(alpha=0.25)

# Same row shares an axis, so the two are directly comparable
for r in range(2):
    lo = min(ax[r, 0].get_ylim()[0], ax[r, 1].get_ylim()[0])
    hi = max(ax[r, 0].get_ylim()[1], ax[r, 1].get_ylim()[1])
    ax[r, 0].set_ylim(lo, hi); ax[r, 1].set_ylim(lo, hi)

fig.suptitle("per-particle x and y time traces -- n=9, omega=3000 rad/s, "
             "a=1470nm(=1d), k_t x100, seed=1\n"
             "thick lines = trapped beads (0 and 8 fixed sensors, 4 driven), "
             "colour = bead index (left to right)",
             fontsize=12)
fig.tight_layout()
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/xy_traces.png"
fig.savefig(out, dpi=140); print("saved", out)

# Numerical summary
for p, name in [(pD, "DLVO"), (pJ, "JKR")]:
    dx = p[:, :, 0] - p[0, :, 0]
    print(f"{name}: y amplitude (bead 4) {np.ptp(p[:,4,1])/2:.4f} d | "
          f"y amplitude (bead 2) {np.ptp(p[:,2,1])/2:.4f} d | "
          f"dx full range {np.ptp(dx):.5f} d | dx RMS {dx.std():.5f} d")
