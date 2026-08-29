import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.animation as animation
import gsd.hoomd, glob
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

n = 9
D = "runs/chain-bend-2d-dlvo__n9-w3000-a1470__f4bef323fba8"          # DLVO, seed=1
J = sorted(glob.glob("runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr__*"))[0]

def load(d):
    t = gsd.hoomd.open(f"{d}/traj_A.gsd", "r")
    return np.array([t[i].particles.position[:n, :2] for i in range(len(t))])

pd_, pj = load(D), load(J)
NF = min(len(pd_), len(pj)); pd_, pj = pd_[:NF], pj[:NF]
x = pd_[0][:, 0]

def bow(p):
    y = p[:, :, 1]
    base = np.linspace(0, 1, n)[None, :] * (y[:, -1:] - y[:, :1]) + y[:, :1]
    return np.abs(y - base).max(axis=1)
bd_, bj = bow(pd_), bow(pj)

ys = (min(pd_[:,:,1].min(), pj[:,:,1].min())-0.05, max(pd_[:,:,1].max(), pj[:,:,1].max())+0.05)
yz = (pj[:,:,1].min()-0.01, pj[:,:,1].max()+0.01)

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.8))
cfg = [(axes[0], pd_, "#1f77b4", ys,
        f"DLVO-only (no bending)\nbow RMS = {np.sqrt((bd_**2).mean()):.4f} d"),
       (axes[1], pj, "#d62728", ys,
        f"DLVO + JKR bending  ★same y axis\n"
        f"bow RMS = {np.sqrt((bj**2).mean()):.4f} d"),
       (axes[2], pj, "#d62728", yz,
        f"JKR zoom ({(ys[1]-ys[0])/(yz[1]-yz[0]):.0f}x)\n"
        f"bends with a smooth beam-like curvature")]
arts = []
for ax, dat, col, yl, title in cfg:
    ln, = ax.plot([], [], "o-", color=col, ms=8, lw=2)
    ax.set_xlim(x.min()-0.4, x.max()+0.4); ax.set_ylim(*yl)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("x [d]"); ax.set_ylabel("y [d]"); ax.set_title(title, fontsize=10)
    arts.append((ln, dat))
txt = axes[0].text(0.02, 0.93, "", transform=axes[0].transAxes, fontsize=9)

fig.suptitle("a = 1470 nm = 1d (maximum amplitude), omega = 3000 rad/s, n = 9 "
             "-- only the bending term toggled ON/OFF", fontsize=11)
fig.tight_layout()

def init():
    for ln, _ in arts: ln.set_data([], [])
    return [a[0] for a in arts] + [txt]

def update(i):
    for ln, dat in arts: ln.set_data(dat[i][:,0], dat[i][:,1])
    txt.set_text(f"frame {i}/{NF}  (20 drive periods)")
    return [a[0] for a in arts] + [txt]

ani = animation.FuncAnimation(fig, update, frames=range(NF), init_func=init, blit=True, interval=50)
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/a1470_jkr_vs_dlvo.gif"
ani.save(out, writer=animation.PillowWriter(fps=20))
print("saved", out)
print(f"bow RMS: DLVO {np.sqrt((bd_**2).mean()):.5f} d   "
      f"JKR {np.sqrt((bj**2).mean()):.5f} d"
      f"  -> {np.sqrt((bd_**2).mean())/np.sqrt((bj**2).mean()):.1f}x")
print(f"full y range: DLVO {np.ptp(pd_[:,:,1]):.4f}   JKR {np.ptp(pj[:,:,1]):.4f}"
      f"  -> {np.ptp(pd_[:,:,1])/np.ptp(pj[:,:,1]):.1f}x")
