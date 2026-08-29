import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

S = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad"
z = np.load(f"{S}/mode_profile.npz")
s, A, Ae, Ap = z["s_idx"], z["A_m"], z["A_e"], z["Apred"]
phi, phip = np.degrees(z["phi_m"]), np.degrees(z["phipred"])
mid = int(z["mid"])
TRAP = [0, 4, 8]

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.0))

# ── (1) amplitude profile ──
a = ax[0]
a.plot(s, Ap, "-", color="gray", lw=6, alpha=0.4, label="linear-response prediction")
a.errorbar(s, A, yerr=Ae, fmt="o", color="#d62728", ms=8, capsize=4,
           label="measured (6 seeds)")
for i in TRAP:
    a.plot(s[i], A[i], "s", ms=14, mfc="none", mec="k", mew=1.6,
           label="trapped bead" if i == 0 else None)
a.set_xlabel("s = i - mid  (distance from the driven bead, in bonds)")
a.set_ylabel("A [d]")
a.set_title(f"(1) amplitude profile -- **+0.01 to 0.04%** against prediction\n"
            f"driven {A[mid]:.3f} d -> ends {A[0]:.3f} d (9.3x attenuation)",
            fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# ── (2) phase profile ──
a = ax[1]
a.plot(s, phip, "-", color="gray", lw=6, alpha=0.4,
       label="linear-response prediction")
a.errorbar(s, phi, yerr=np.degrees(z["phi_e"]), fmt="o", color="#1f77b4", ms=8,
           capsize=4, label="measured (6 seeds)")
for i in TRAP:
    a.plot(s[i], phi[i], "s", ms=14, mfc="none", mec="k", mew=1.6)
a.axhline(0, color="gray", lw=0.6)
a.set_xlabel("s = i - mid"); a.set_ylabel("phi - phi_driven  [deg]")
a.set_title("(2) phase lag -- within **0.06 deg** of prediction\n"
            "the further from the centre, the further behind (viscous propagation). "
            "The ends jump because of the traps", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# ── (3) shape of the attenuation (log) ──
a = ax[2]
a.semilogy(np.abs(s), A, "o", color="#d62728", ms=8, label="measured")
a.semilogy(np.abs(s), Ap, "-", color="gray", lw=5, alpha=0.4,
           label="linear-response prediction")
# Exponential fit over the free range (|s|=0..3) only -- |s|=4 is trap-dominated
# and is looked at separately
m = np.abs(s) <= 3
k = np.polyfit(np.abs(s)[m], np.log(A[m]), 1)[0]
xs = np.linspace(0, 4, 50)
a.semilogy(xs, A[mid] * np.exp(k * xs), "--", color="#2ca02c", lw=1.6,
           label=f"exp fit |s|<=3: decay length {-1/k:.2f} bonds")
a.plot(4, A[0], "s", ms=14, mfc="none", mec="k", mew=1.6, label="trapped ends")
a.set_xlabel("|s|  (distance from the centre)"); a.set_ylabel("A [d]  (log)")
a.set_title("(3) the attenuation is not a single exponential\n"
            "gentle for |s|<=3, then a sharp drop at |s|=4 (the traps) -- the "
            "boundary condition dominates", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, which="both")

fig.suptitle("per-bead A*sin(omega*t+phi) fit for the JKR chain -- n=9, "
             "omega=3000 rad/s, a=1d, k_t x100, 6 seeds\n"
             "thick grey line = analytic linear response "
             "(i*omega*gamma*I + A_bend + T)y_hat = k_t*a*e_mid  -- no free parameters",
             fontsize=11)
fig.tight_layout()
out = f"{S}/mode_profile.png"
fig.savefig(out, dpi=140); print("saved", out)
print(f"exponential decay length (|s|<=3): {-1/k:.3f} bonds = {-1/k*1.0076:.3f} d")
