import json, glob
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np
import gsd.hoomd

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

def collect(pat):
    Kp, Kpp, y, pred = [], [], [], None
    for d in sorted(glob.glob(pat)):
        try: met = json.load(open(f"{d}/metrics.json"))
        except FileNotFoundError: continue
        o = {x["name"]: x for x in met["observables"]}
        Kp.append(o["K_prime"]["measured"]); Kpp.append(o["K_doubleprime"]["measured"])
        y.append(o["y_response"]["measured"]); pred = o["K_prime"]["predicted"]
    return np.array(Kp), np.array(Kpp), np.array(y), pred

# DLVO-only: only the cycles=20 runs (n_prod >= 600000)
Kp_d, Kpp_d, y_d = [], [], []
for d in sorted(glob.glob("runs/chain-bend-2d-dlvo__n9-w3000-a632__*")):
    spec = json.load(open(f"{d}/spec.json"))
    if spec["numerics"]["n_prod"] < 600000: continue
    met = json.load(open(f"{d}/metrics.json"))
    o = {x["name"]: x for x in met["observables"]}
    Kp_d.append(o["K_prime"]["measured"]); Kpp_d.append(o["K_doubleprime"]["measured"])
    y_d.append(o["y_response"]["measured"])
Kp_d, Kpp_d, y_d = np.array(Kp_d), np.array(Kpp_d), np.array(y_d)
Kp_j, Kpp_j, y_j, pred_j = collect("runs/chain-bend-2d-dlvo__n9-w3000-a632-jkr__*")

def ms(a): return a.mean(), a.std(ddof=1)/np.sqrt(len(a))

fig = plt.figure(figsize=(15, 5.2))
gs = fig.add_gridspec(1, 3, width_ratios=[1, 1, 1.25])

# ── (1) K' comparison (log) ──
ax = fig.add_subplot(gs[0])
labels = ["DLVO-only", "DLVO+JKR bending"]
mK = [ms(Kp_d)[0], ms(Kp_j)[0]]; eK = [ms(Kp_d)[1], ms(Kp_j)[1]]
mL = [ms(Kpp_d)[0], ms(Kpp_j)[0]]; eL = [ms(Kpp_d)[1], ms(Kpp_j)[1]]
x = np.arange(2); w = 0.34
ax.bar(x-w/2, mK, w, yerr=eK, capsize=5, color="#1f77b4", label="K' (storage)")
ax.bar(x+w/2, mL, w, yerr=eL, capsize=5, color="#ff7f0e", alpha=0.8, label="K'' (loss)")
ax.axhline(pred_j, color="green", ls="--", lw=1.6,
           label=f"JKR static-limit prediction {pred_j:.0f}")
ax.set_yscale("log"); ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("K [kT/d^2]  (log)")
ax.set_title(f"K' differs by a factor {mK[1]/mK[0]:.0f}\n"
             f"DLVO {mK[0]:.0f}±{eK[0]:.0f}  vs  JKR {mK[1]:.0f}±{eK[1]:.0f}", fontsize=10)
for xi, v in zip(x-w/2, mK): ax.text(xi, v*1.35, f"{v:.0f}", ha="center", fontsize=8)
ax.legend(fontsize=8, loc="upper left"); ax.grid(alpha=0.3, axis="y")

# ── (2) response amplitude (against the linear-response prediction) ──
ax = fig.add_subplot(gs[1])
my = [ms(y_d)[0], ms(y_j)[0]]; ey = [ms(y_d)[1], ms(y_j)[1]]
ax.bar(x, my, 0.5, yerr=ey, capsize=5, color=["#1f77b4", "#d62728"])
ax.axhline(0.01506, color="green", ls="--", lw=1.6, label="JKR linear-response prediction 0.01506 d")
ax.axhline(0.4300, color="gray", ls=":", lw=1.4, label="drive amplitude 0.43 d")
ax.set_xticks(x); ax.set_xticklabels(labels, fontsize=9)
ax.set_ylabel("|ŷ| [d]")
ax.set_title(f"response amplitude -- JKR is {100*(my[1]-0.01506)/0.01506:+.1f}% "
             f"against prediction\n"
             "the stiffer it is, the less it moves however hard the trap pulls",
             fontsize=10)
ax.legend(fontsize=8); ax.grid(alpha=0.3, axis="y")

# ── (3) instantaneous chain shapes (NOT time-averaged) ──
ax = fig.add_subplot(gs[2])
n = 9
for pat, col, lab in [("runs/chain-bend-2d-dlvo__n9-w3000-a632__afe1dd36499b", "#1f77b4", "DLVO-only"),
                      (sorted(glob.glob("runs/chain-bend-2d-dlvo__n9-w3000-a632-jkr__*"))[0], "#d62728", "DLVO+JKR")]:
    traj = gsd.hoomd.open(f"{pat}/traj_A.gsd", "r")
    idx = np.linspace(len(traj)//2, len(traj)-1, 6, dtype=int)
    for k, i in enumerate(idx):
        p = traj[int(i)].particles.position[:n]
        ax.plot(p[:,0], p[:,1], "o-", color=col, alpha=0.45, ms=4, lw=1.2,
                label=lab if k == 0 else None)
ax.axhline(0, color="gray", lw=0.5)
ax.set_xlabel("x [d]"); ax.set_ylabel("y [d]")
ax.set_title("six instantaneous chain shapes each (late production)\n"
             "JKR (red) is nearly straight -- too stiff to bend\n"
             "DLVO (blue) is heavily disordered by thermal noise", fontsize=10)
ax.legend(fontsize=8)

fig.suptitle("JKR control -- same geometry, same seeds, only the bending term "
             "toggled ON/OFF (n=9, omega=3000 rad/s, a=632nm, 6 seeds)",
             fontsize=11)
fig.tight_layout()
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/jkr_compare.png"
fig.savefig(out, dpi=140); print("saved", out)
