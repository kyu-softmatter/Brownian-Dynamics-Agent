import pickle
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name
        break
plt.rcParams["axes.unicode_minus"] = False

SCRATCH = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad"
omega_rows = pickle.load(open(f"{SCRATCH}/omega_sweep_rows.pkl", "rb"))
amp_rows = pickle.load(open(f"{SCRATCH}/amp_sweep_rows.pkl", "rb"))

fig, axes = plt.subplots(1, 2, figsize=(13.5, 5.5))

# ── (1) omega sweep ──
ax = axes[0]
w = np.array([r[0] for r in omega_rows], dtype=float)
Kp = np.array([r[1] for r in omega_rows]); Kp_sem = np.array([r[2] for r in omega_rows])
Kpp = np.array([r[3] for r in omega_rows]); Kpp_sem = np.array([r[4] for r in omega_rows])
n_seeds = [r[5] for r in omega_rows]
ax.errorbar(w, Kp, yerr=Kp_sem, fmt="o-", capsize=4, color="#1f77b4", label="K'(ω)")
ax.errorbar(w, Kpp, yerr=Kpp_sem, fmt="s--", capsize=4, color="#ff7f0e", alpha=0.7, label="K''(ω)")
ax.axhline(0, color="gray", lw=1.2, label="G1 prediction: K'=0")
ax.set_xscale("log")
ax.set_xlabel("ω [rad/s]"); ax.set_ylabel("K [kT/d²]")
ax.set_title("omega sweep (n=9, a=632nm)\nclosest to 0 and significant (4-5 sigma) at "
             "low frequency (10-100)\n"
             "-> no elastic plateau in the quasi-static limit (consistent with G1); "
             "increasingly viscoelastic at higher frequency")
for wi, ni, kpi in zip(w, n_seeds, Kp):
    ax.annotate(f"N={ni}", (wi, kpi), textcoords="offset points", xytext=(6, 8), fontsize=7,
                color="gray")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# ── (2) amplitude sweep ──
ax = axes[1]
a = np.array([r[0] for r in amp_rows], dtype=float) / 1470.0    # a/d
Kp2 = np.array([r[1] for r in amp_rows]); Kp2_sem = np.array([r[2] for r in amp_rows])
n2 = [r[3] for r in amp_rows]
ax.errorbar(a, Kp2, yerr=Kp2_sem, fmt="o-", capsize=4, color="#2ca02c")
ax.axhline(0, color="gray", lw=1.2, label="G1 prediction: K'=0")
ax.set_xlabel("a/d (drive amplitude / particle diameter)"); ax.set_ylabel("K' [kT/d^2]")
ax.set_title("amplitude sweep (n=9, omega=3000 rad/s)\nas the amplitude grows K' "
             "shrinks toward 0 and the SEM shrinks too\n"
             "(small amplitude = noise dominated, large amplitude = precisely close "
             "to 0) -- consistent with G1")
for ai, ni, kpi, semi in zip(a, n2, Kp2, Kp2_sem):
    ax.annotate(f"N={ni}", (ai, kpi), textcoords="offset points", xytext=(6, 8), fontsize=7,
                color="gray")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.suptitle("chain-bend-2d-dlvo -- K'(omega) and K'(a) for a chain bonded by DLVO "
             "alone (n=9 fixed)", fontsize=11)
fig.tight_layout()
out = f"{SCRATCH}/dlvo_final_sweeps.png"
fig.savefig(out, dpi=140)
print("saved", out)
