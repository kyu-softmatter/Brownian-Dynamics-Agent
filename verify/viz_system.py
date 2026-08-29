import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

OM = 18453.1                       # omega*gamma (reduced, gamma*=1) = drag dissipation of one bead
D = dict(tot=1.838e4, chain=-56.96, store=641.2, Kp_lock=-9.203)
J = dict(tot=7.559e4, chain=5.727e4, store=1.276e5, Kp_lock=1.128e5)

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.2))

# ── (1) dissipation breakdown (stacked bars) ──
a = ax[0]
labs = ["bead alone\n(no chain)", "DLVO chain", "JKR chain"]
drag = [OM, OM, OM]
chain = [0.0, D["chain"], J["chain"]]
x = np.arange(3)
a.bar(x, drag, 0.55, color="#9ecae1",
      label="the driven bead's own solvent drag (omega*gamma)")
a.bar(x, chain, 0.55, bottom=drag, color="#d62728",
      label="dissipation carried by the chain")
for i, (dr, ch) in enumerate(zip(drag, chain)):
    tot = dr + ch
    a.text(i, tot * 1.04, f"{tot:.3g}\n({tot/OM:.2f}×)", ha="center", fontsize=9)
a.set_ylabel("K'' [kT/d^2]  (dissipation per cycle is proportional to this)")
a.set_ylim(0, 9.5e4)
a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
a.set_title("★ (1) dissipation breakdown for the whole system\n"
            "attaching the DLVO chain leaves the dissipation **the same as the bead "
            "alone** (0.996x)\n"
            "JKR is 4.10x -- the chain carries 75.8% of the dissipation", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── (2) cross-check via two independent routes ──
a = ax[1]
w = 0.35; x = np.arange(2)
tot = [D["tot"], J["tot"]]
recon = [D["chain"] + OM, J["chain"] + OM]
a.bar(x - w/2, tot, w, color="#2ca02c",
      label="(1) energy integral  cyclic-int F.dy/(pi|y_hat|^2)")
a.bar(x + w/2, recon, w, color="#ff7f0e",
      label="(2) lock-in K''_chain + omega*gamma")
for xi, v in zip(x - w/2, tot): a.text(xi, v*1.03, f"{v:.4g}", ha="center", fontsize=8)
for xi, v in zip(x + w/2, recon): a.text(xi, v*1.03, f"{v:.4g}", ha="center", fontsize=8)
for i, (t_, r_) in enumerate(zip(tot, recon)):
    a.text(i, max(t_, r_)*1.14, f"{100*(t_-r_)/t_:+.2f}%", ha="center", fontsize=10,
           color="darkgreen", weight="bold")
a.set_xticks(x); a.set_xticklabels(["DLVO chain", "JKR chain"], fontsize=9)
a.set_ylabel("K″ [kT/d²]"); a.set_ylim(0, 9.5e4)
a.set_title("(2) the two routes agree with no free parameters\n"
            "whole-system energy integral vs drive-point lock-in + drag\n"
            "DLVO −0.11% · JKR −0.17%", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── (3) storage: whole system vs drive point ──
a = ax[2]
x = np.arange(2)
sys_ = [D["store"], J["store"]]
loc = [abs(D["Kp_lock"]), abs(J["Kp_lock"])]
a.bar(x - w/2, sys_, w, color="#6a51a3",
      label="whole-system storage K'_sys (the 2*omega component of U)")
a.bar(x + w/2, loc, w, color="#1f77b4", label="drive-point |K'| (lock-in)")
for xi, v in zip(x - w/2, sys_): a.text(xi, v*1.5, f"{v:.4g}", ha="center", fontsize=8)
for xi, v in zip(x + w/2, loc): a.text(xi, v*1.5, f"{v:.3g}", ha="center", fontsize=8)
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(["DLVO chain", "JKR chain"], fontsize=9)
a.set_ylabel("K' [kT/d^2]  (log)"); a.set_ylim(1, 1e6)
a.set_title("(3) storage: whole system != drive point\n"
            "JKR 1.13x (bending stores across the whole chain)\n"
            "★ DLVO stores 641 even though its drive-point K' ~ 0 -- bond stretching",
            fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y", which="both")

fig.suptitle("system-level rheology -- the energy budget of the whole system "
             "(no per-bead fitting)\n"
             "n=9, omega=3000 rad/s, a=1d, k_t x100, 6 seeds", fontsize=12)
fig.tight_layout()
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/system_rheology.png"
fig.savefig(out, dpi=140); print("saved", out)
