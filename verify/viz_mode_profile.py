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

# ── ① 진폭 프로파일 ──
a = ax[0]
a.plot(s, Ap, "-", color="gray", lw=6, alpha=0.4, label="선형응답 예측")
a.errorbar(s, A, yerr=Ae, fmt="o", color="#d62728", ms=8, capsize=4, label="측정 (시드 6)")
for i in TRAP:
    a.plot(s[i], A[i], "s", ms=14, mfc="none", mec="k", mew=1.6,
           label="트랩된 비드" if i == 0 else None)
a.set_xlabel("s = i − mid  (구동 비드로부터의 거리 [결합])")
a.set_ylabel("A [d]")
a.set_title(f"① 진폭 프로파일 — 예측과 **+0.01~0.04%**\n"
            f"구동 {A[mid]:.3f} d → 양끝 {A[0]:.3f} d (9.3배 감쇠)", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# ── ② 위상 프로파일 ──
a = ax[1]
a.plot(s, phip, "-", color="gray", lw=6, alpha=0.4, label="선형응답 예측")
a.errorbar(s, phi, yerr=np.degrees(z["phi_e"]), fmt="o", color="#1f77b4", ms=8,
           capsize=4, label="측정 (시드 6)")
for i in TRAP:
    a.plot(s[i], phi[i], "s", ms=14, mfc="none", mec="k", mew=1.6)
a.axhline(0, color="gray", lw=0.6)
a.set_xlabel("s = i − mid"); a.set_ylabel("φ − φ_구동  [°]")
a.set_title("② 위상 지연 — 예측과 **0.06° 이내**\n"
            "중심에서 멀수록 뒤처진다 (점성 전파). 양끝은 트랩 때문에 급증", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3)

# ── ③ 감쇠 형태 (로그) ──
a = ax[2]
a.semilogy(np.abs(s), A, "o", color="#d62728", ms=8, label="측정")
a.semilogy(np.abs(s), Ap, "-", color="gray", lw=5, alpha=0.4, label="선형응답 예측")
# 자유 구간(|s|=0..3)만 지수 적합 — |s|=4 는 트랩이 지배해 따로 본다
m = np.abs(s) <= 3
k = np.polyfit(np.abs(s)[m], np.log(A[m]), 1)[0]
xs = np.linspace(0, 4, 50)
a.semilogy(xs, A[mid] * np.exp(k * xs), "--", color="#2ca02c", lw=1.6,
           label=f"지수적합 |s|≤3: 감쇠길이 {-1/k:.2f} 결합")
a.plot(4, A[0], "s", ms=14, mfc="none", mec="k", mew=1.6, label="트랩된 양끝")
a.set_xlabel("|s|  (중심으로부터 거리)"); a.set_ylabel("A [d]  (로그)")
a.set_title("③ 감쇠는 단일 지수가 아니다\n"
            "|s|≤3 은 완만, |s|=4(트랩)에서 급락 — 경계조건이 지배", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, which="both")

fig.suptitle("JKR 사슬의 비드별 A·sin(ωt+φ) 피팅 — n=9, ω=3000 rad/s, a=1d, k_t×100, 시드 6개\n"
             "회색 굵은 선 = 해석적 선형응답 (iωγI + A_bend + T)ŷ = k_t·a·e_mid  — 자유 파라미터 없음",
             fontsize=11)
fig.tight_layout()
out = f"{S}/mode_profile.png"
fig.savefig(out, dpi=140); print("saved", out)
print(f"지수 감쇠길이(|s|<=3): {-1/k:.3f} 결합 = {-1/k*1.0076:.3f} d")
