import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

S = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad"
R = np.load(f"{S}/corr_res.npy", allow_pickle=True).item()

fig = plt.figure(figsize=(15, 8.6))
gs = fig.add_gridspec(2, 3, height_ratios=[1, 0.95], hspace=0.42, wspace=0.32)

# ── top row: three displacement correlation matrices ──
mats = [("dlvo_full", "DLVO-only  (n=9)",
         "no bending -> each bead moves on its own"),
        ("jkr_full",  "DLVO+JKR  (n=9)",
         "bending present -> the chain moves as a whole"),
        ("n25_full",  "DLVO-only  (n=25)",
         "still uncorrelated even for a long chain")]
for k, (key, title, sub) in enumerate(mats):
    ax = fig.add_subplot(gs[0, k])
    C = R[key]["C"]
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    off = C[np.triu_indices_from(C, k=1)]
    ax.set_title(f"{title}\n{sub}\noff-diagonal mean = {off.mean():+.3f}", fontsize=9.5)
    ax.set_xlabel("bead j"); ax.set_ylabel("bead i")
    fig.colorbar(im, ax=ax, fraction=0.046)

# ── bottom left: tangent correlation ──
ax = fig.add_subplot(gs[1, :2])
for key, col, lab in [("dlvo_full", "#1f77b4", "DLVO-only n=9"),
                      ("jkr_full", "#d62728", "DLVO+JKR n=9"),
                      ("n25_full", "#2ca02c", "DLVO-only n=25")]:
    c = R[key]["c"]
    ax.plot(np.arange(len(c)), c, "o-", color=col, label=lab, ms=5)
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_xlabel("bond separation m"); ax.set_ylabel("<t_i . t_{i+m}>")
ax.set_ylim(0.6, 1.03)
ax.set_title("tangent-tangent correlation -- weakly discriminating in this system\n"
             "both chains are nearly straight along x (displacement << bond length), "
             "so cos(theta) ~ 1.\n"
             "n=25 drops to 0.80 at m=1, which shows the absence of bending, but is "
             "flat beyond that -- because of the trap constraint",
             fontsize=9.5)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# ── bottom right: summary bars ──
ax = fig.add_subplot(gs[1, 2])
labs = ["DLVO\nn=9", "JKR\nn=9", "DLVO\nn=25"]
vals = [R[k]["C"][np.triu_indices_from(R[k]["C"], k=1)].mean()
        for k in ("dlvo_full", "jkr_full", "n25_full")]
cols = ["#1f77b4", "#d62728", "#2ca02c"]
ax.bar(range(3), vals, color=cols)
ax.axhline(0, color="k", lw=1)
ax.set_xticks(range(3)); ax.set_xticklabels(labs, fontsize=9)
ax.set_ylabel("off-diagonal mean of the displacement correlation")
ax.set_ylim(-0.15, 1.0)
for i, v in enumerate(vals):
    ax.text(i, v + (0.04 if v > 0 else -0.08), f"{v:+.3f}", ha="center", fontsize=9)
ax.set_title("★ the real discriminator\none bending term takes it from -0.02 to +0.91",
             fontsize=9.5)
ax.grid(alpha=0.3, axis="y")

fig.suptitle("particle position correlations -- displacement correlation matrix "
             "<dy_i dy_j> and tangent correlation <t_i . t_j> "
             "(omega=3000, a=1d / for n=25, omega=10, a=632nm)", fontsize=11.5)
out = f"{S}/correlations.png"
fig.savefig(out, dpi=140, bbox_inches="tight"); print("saved", out)
