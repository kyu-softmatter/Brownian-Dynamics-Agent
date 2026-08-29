import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import matplotlib.animation as animation
import gsd.hoomd
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

n = 9
D = "runs/chain-bend-2d-dlvo__n9-w3000-a632__afe1dd36499b"
J = "runs/chain-bend-2d-dlvo__n9-w3000-a632-jkr__0247bac45e98"

def load(d):
    t = gsd.hoomd.open(f"{d}/traj_A.gsd", "r")
    return np.array([t[i].particles.position[:n, :2] for i in range(len(t))])

pd_, pj = load(D), load(J)
NF = min(len(pd_), len(pj))
pd_, pj = pd_[:NF], pj[:NF]
x = pd_[0][:, 0]

ylim_shared = (min(pd_[:,:,1].min(), pj[:,:,1].min()) - 0.02,
               max(pd_[:,:,1].max(), pj[:,:,1].max()) + 0.02)
ylim_zoom = (pj[:,:,1].min() - 0.005, pj[:,:,1].max() + 0.005)

fig, axes = plt.subplots(1, 3, figsize=(15.5, 4.6))
cfg = [
    (axes[0], pd_, "#1f77b4", ylim_shared,
     "DLVO-only (굽힘항 없음)\n열요동으로 사슬이 크게 흐트러진다"),
    (axes[1], pj, "#d62728", ylim_shared,
     "DLVO + JKR 굽힘  ★같은 y축\n거의 직선 — 트랩이 끌어도 안 휜다"),
    (axes[2], pj, "#d62728", ylim_zoom,
     f"JKR 확대 (y축 {(ylim_shared[1]-ylim_shared[0])/(ylim_zoom[1]-ylim_zoom[0]):.0f}배 확대)\n"
     "움직이긴 한다 — 매끈한 빔형 곡률로"),
]
arts = []
for ax, dat, col, yl, title in cfg:
    ln, = ax.plot([], [], "o-", color=col, ms=8, lw=2)
    ax.set_xlim(x.min()-0.4, x.max()+0.4); ax.set_ylim(*yl)
    ax.axhline(0, color="gray", lw=0.5)
    ax.set_xlabel("x [d]"); ax.set_ylabel("y [d]")
    ax.set_title(title, fontsize=10)
    arts.append((ln, dat))
txt = axes[0].text(0.02, 0.93, "", transform=axes[0].transAxes, fontsize=9)

fig.suptitle("같은 기하·같은 시드·같은 구동 (n=9, ω=3000 rad/s, a=632nm) — 굽힘항만 ON/OFF",
             fontsize=11)
fig.tight_layout()

def init():
    for ln, _ in arts: ln.set_data([], [])
    return [a[0] for a in arts] + [txt]

def update(i):
    for ln, dat in arts:
        ln.set_data(dat[i][:, 0], dat[i][:, 1])
    txt.set_text(f"frame {i}/{NF}  (20 구동주기)")
    return [a[0] for a in arts] + [txt]

ani = animation.FuncAnimation(fig, update, frames=range(NF), init_func=init,
                              blit=True, interval=50)
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/jkr_vs_dlvo_motion.gif"
ani.save(out, writer=animation.PillowWriter(fps=20))
print("saved", out)
print(f"y범위: DLVO {pd_[:,:,1].ptp():.4f} d   JKR {pj[:,:,1].ptp():.4f} d  → {pd_[:,:,1].ptp()/pj[:,:,1].ptp():.1f}배")
