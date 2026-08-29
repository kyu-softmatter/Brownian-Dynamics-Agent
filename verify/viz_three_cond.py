import pickle
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

S="/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad"
R = pickle.load(open(f"{S}/three_cond.pkl","rb"))
labs = [o["lab"] for o in R]
x = np.arange(3); w = 0.35

fig, ax = plt.subplots(2, 2, figsize=(14, 9.5))

# ── (1) tracking ratio ──
a = ax[0,0]
yd = [100*o["yD"][0] for o in R]; yj = [100*o["yJ"][0] for o in R]
a.bar(x-w/2, yd, w, color="#1f77b4", label="DLVO-only")
a.bar(x+w/2, yj, w, color="#d62728", label="DLVO+JKR")
for xi,v in zip(x-w/2,yd): a.text(xi,v+2,f"{v:.0f}%",ha="center",fontsize=9)
for xi,v in zip(x+w/2,yj): a.text(xi,v+2,f"{v:.0f}%",ha="center",fontsize=9)
a.axhline(100, color="gray", ls=":", lw=1.2)
a.set_xticks(x); a.set_xticklabels(labs, fontsize=9); a.set_ylim(0,118)
a.set_ylabel("tracking ratio |y_hat|/a  [%]")
a.set_title("(1) tracking ratio by drive protocol\nthe default trap follows only "
            "3.5-27% of the command", fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── (2) how many sigma DLVO sits from 0 (a protocol-independent measure) ──
a = ax[0,1]
z = [o["z"] for o in R]
cols = ["#ff7f0e" if v>2 else "#2ca02c" for v in z]
a.bar(x, z, 0.5, color=cols)
for xi,v in zip(x,z): a.text(xi, v+0.08, f"{v:.2f}σ", ha="center", fontsize=10)
a.axhline(2, color="crimson", ls="--", lw=1.6, label="2 sigma -- below this it is indistinguishable from 0")
a.set_xticks(x); a.set_xticklabels(labs, fontsize=9); a.set_ylim(0, 3.0)
a.set_ylabel("distance of the DLVO stiffness from 0 [sigma]")
a.set_title("★ (2) the better the tracking, the closer DLVO converges to 0\n"
            "2.28 -> 0.84 -> 0.66 sigma  (tracking 27% -> 100% -> 100%)\n"
            "the default trap's 'finite K'' was an artefact of failing to track",
            fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── (3) JKR/DLVO ratio (protocol-independent, log) ──
a = ax[1,0]
ratio = [abs(o["J"][0]/o["D"][0]) for o in R]
a.bar(x, ratio, 0.5, color="#6a51a3")
for xi,v,o in zip(x,ratio,R):
    a.text(xi, v*1.4, f"{v:,.0f}x", ha="center", fontsize=10)
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
a.set_ylim(100, 3e5)
a.set_ylabel("|JKR| / |DLVO|  (log)")
a.set_title("(3) JKR dominates under all three protocols\n"
            "400x -> 12,000x -> 2,500x\n"
            "⚠ the denominator is close to 0 so the ratio itself is unstable -- "
            "read the order of magnitude only", fontsize=10.5)
a.grid(alpha=0.3, axis="y", which="both")

# ── (4) bow -- its discriminating power inverts with the protocol ──
a = ax[1,1]
bd = [o["bD"][0] for o in R]; bj = [o["bJ"][0] for o in R]
a.bar(x-w/2, bd, w, color="#1f77b4", label="DLVO-only")
a.bar(x+w/2, bj, w, color="#d62728", label="DLVO+JKR")
for xi,(d_,j_) in zip(x, zip(bd,bj)):
    a.text(xi, max(d_,j_)*1.35, f"{d_/j_:.1f}x", ha="center", fontsize=10,
           color="darkgreen")
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
a.set_ylabel("bow RMS [d]  (log)"); a.set_ylim(5e-3, 3)
a.set_title("(4) the discriminating power of bow (shape) depends on the protocol\n"
            "15.4x -> 1.4x -> 1.9x\n"
            "★ once the deformation is imposed, shape carries no information -- "
            "look at the force", fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y", which="both")

fig.suptitle("all three drive protocols -- n=9, omega=3000 rad/s, a=1470nm(=1d), "
             "6 seeds\n"
             "the two trap protocols give K' (stiffness at the drive point), "
             "position gives K_transfer (transfer stiffness) -- "
             "different quantities, so never compare absolute values across "
             "conditions; only the DLVO/JKR contrast within each condition is valid",
             fontsize=11.5)
fig.tight_layout()
out=f"{S}/three_conditions.png"
fig.savefig(out, dpi=140); print("saved", out)
