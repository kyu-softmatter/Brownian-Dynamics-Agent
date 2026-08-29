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
D = pickle.load(open(f"{S}/w10_rows.pkl","rb"))
nn  = np.array([r[0] for r in D["n"]]);  Kn = np.array([r[1] for r in D["n"]]);  Kne = np.array([r[2] for r in D["n"]])
Bn  = np.array([r[3] for r in D["n"]]);  Bne = np.array([r[4] for r in D["n"]])
aa  = np.array([r[0] for r in D["amp"]])/1470.0
Ka  = np.array([r[1] for r in D["amp"]]); Kae = np.array([r[2] for r in D["amp"]])
Ba  = np.array([r[3] for r in D["amp"]]); Bae = np.array([r[4] for r in D["amp"]])

JKR_K, JKR_B = 43064.0, 0.006387     # the omega=3000 control (reference line)

fig, ax = plt.subplots(2, 2, figsize=(13.5, 9))

# ① K' vs n
a0 = ax[0,0]
a0.errorbar(nn, Kn, yerr=Kne, fmt="o-", capsize=4, color="#1f77b4", lw=2, ms=7)
a0.axhline(0, color="gray", lw=1.3, label="G1 prediction: K'=0")
a0.set_xlabel("n (bead count)"); a0.set_ylabel("K' [kT/d^2]")
a0.set_title("(1) the longer the chain, the closer K' converges to 0\n"
             "76.9 -> 12.7 (n=5->25), a 6x reduction -- longer chains are closer to "
             "a free hinge", fontsize=10)
for x,y,e in zip(nn,Kn,Kne): a0.annotate(f"{y:.0f}±{e:.0f}",(x,y),xytext=(6,7),
                                          textcoords="offset points",fontsize=8)
a0.legend(fontsize=9); a0.grid(alpha=0.3)

# (2) K' vs amplitude
a1 = ax[0,1]
a1.errorbar(aa, Ka, yerr=Kae, fmt="s-", capsize=4, color="#2ca02c", lw=2, ms=7)
a1.axhline(0, color="gray", lw=1.3, label="G1 prediction: K'=0")
a1.set_xlabel("a/d (drive amplitude/diameter)"); a1.set_ylabel("K' [kT/d^2]")
a1.set_title("(2) the larger the amplitude, the closer K' converges to 0 (and the "
             "error shrinks with it)\n"
             "81±31 -> 9.2±2.6 -- the more precise it gets, the closer to 0",
             fontsize=10)
for x,y,e in zip(aa,Ka,Kae): a1.annotate(f"{y:.0f}±{e:.0f}",(x,y),xytext=(6,7),
                                          textcoords="offset points",fontsize=8)
a1.legend(fontsize=9); a1.grid(alpha=0.3)

# (3) bow vs n / amplitude  (JKR reference line)
a2 = ax[1,0]
a2.errorbar(nn, Bn, yerr=Bne, fmt="o-", capsize=4, color="#1f77b4", lw=2, ms=7, label="DLVO-only (ω=10)")
a2.axhline(JKR_B, color="#d62728", ls="--", lw=1.8, label=f"JKR control {JKR_B:.4f} d (omega=3000)")
a2.set_xlabel("n (bead count)"); a2.set_ylabel("bow RMS [d]")
a2.set_yscale("log")
a2.set_title(f"(3) bow -- the DLVO chain bends {Bn.mean()/JKR_B:.0f}x more than JKR\n"
             "longer chains bend more (0.536->0.728). Direct evidence that there is "
             "no bending stiffness", fontsize=10)
a2.legend(fontsize=8); a2.grid(alpha=0.3)

# (4) summary -- order-of-magnitude comparison of K'
a3 = ax[1,1]
cats = ["DLVO\nn=25\n(longest chain)", "DLVO\na=d\n(max amplitude)",
        "DLVO\nn=9,a=632", "JKR bending ON\n(omega=3000)"]
vals = [Kn[-1], Ka[-1], Kn[1], JKR_K]
errs = [Kne[-1], Kae[-1], Kne[1], 2516.0]
cols = ["#1f77b4","#2ca02c","#7f7f7f","#d62728"]
a3.bar(range(4), vals, yerr=errs, capsize=5, color=cols)
a3.set_yscale("log"); a3.set_xticks(range(4)); a3.set_xticklabels(cats, fontsize=8)
a3.set_ylabel("K' [kT/d^2]  (log)")
a3.axhline(1.391e6, color="k", ls=":", lw=1.5,
           label="JKR bending stiffness kappa_theta* = 1.4e6")
a3.set_title(f"(4) order-of-magnitude comparison -- the cleanest DLVO condition "
             f"(n=25) is 1/{JKR_K/Kn[-1]:.0f} of JKR\n"
             "under no condition does it come near beam elasticity", fontsize=10)
for i,(v,e) in enumerate(zip(vals,errs)): a3.text(i, v*1.5, f"{v:.0f}", ha="center", fontsize=8)
a3.legend(fontsize=8); a3.grid(alpha=0.3, axis="y")

fig.suptitle("chain-bend-2d-dlvo -- omega=10 rad/s (quasi-static), full sweep, "
             "6-seed ensemble, 48 runs total",
             fontsize=12)
fig.tight_layout()
out=f"{S}/w10_results.png"; fig.savefig(out, dpi=140); print("saved", out)
