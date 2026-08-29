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

# ── ① omega 스윕 ──
ax = axes[0]
w = np.array([r[0] for r in omega_rows], dtype=float)
Kp = np.array([r[1] for r in omega_rows]); Kp_sem = np.array([r[2] for r in omega_rows])
Kpp = np.array([r[3] for r in omega_rows]); Kpp_sem = np.array([r[4] for r in omega_rows])
n_seeds = [r[5] for r in omega_rows]
ax.errorbar(w, Kp, yerr=Kp_sem, fmt="o-", capsize=4, color="#1f77b4", label="K'(ω)")
ax.errorbar(w, Kpp, yerr=Kpp_sem, fmt="s--", capsize=4, color="#ff7f0e", alpha=0.7, label="K''(ω)")
ax.axhline(0, color="gray", lw=1.2, label="G1 예측: K'=0")
ax.set_xscale("log")
ax.set_xlabel("ω [rad/s]"); ax.set_ylabel("K [kT/d²]")
ax.set_title("ω 스윕 (n=9, a=632nm)\n저주파(10~100)에서 0에 가장 가깝고 유의미(4~5σ)\n"
             "→ 준정적 극한에서 탄성 고원 없음(G1과 일치), 고주파일수록 점탄성적 증가")
for wi, ni, kpi in zip(w, n_seeds, Kp):
    ax.annotate(f"N={ni}", (wi, kpi), textcoords="offset points", xytext=(6, 8), fontsize=7,
                color="gray")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# ── ② amplitude 스윕 ──
ax = axes[1]
a = np.array([r[0] for r in amp_rows], dtype=float) / 1470.0    # a/d
Kp2 = np.array([r[1] for r in amp_rows]); Kp2_sem = np.array([r[2] for r in amp_rows])
n2 = [r[3] for r in amp_rows]
ax.errorbar(a, Kp2, yerr=Kp2_sem, fmt="o-", capsize=4, color="#2ca02c")
ax.axhline(0, color="gray", lw=1.2, label="G1 예측: K'=0")
ax.set_xlabel("a/d (구동진폭 / 입자지름)"); ax.set_ylabel("K' [kT/d²]")
ax.set_title("진폭 스윕 (n=9, ω=3000 rad/s)\n진폭이 커질수록 K'는 0쪽으로 줄고 SEM도 줄어듦\n"
             "(작은 진폭=노이즈 지배, 큰 진폭=정밀하게 0에 근접) — G1과 일치")
for ai, ni, kpi, semi in zip(a, n2, Kp2, Kp2_sem):
    ax.annotate(f"N={ni}", (ai, kpi), textcoords="offset points", xytext=(6, 8), fontsize=7,
                color="gray")
ax.legend(fontsize=8); ax.grid(alpha=0.3)

fig.suptitle("chain-bend-2d-dlvo — DLVO만으로 연결한 사슬의 K'(ω), K'(a) (n=9 고정)", fontsize=11)
fig.tight_layout()
out = f"{SCRATCH}/dlvo_final_sweeps.png"
fig.savefig(out, dpi=140)
print("saved", out)
