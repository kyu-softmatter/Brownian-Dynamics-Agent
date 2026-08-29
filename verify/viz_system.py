import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm
import numpy as np

for f in fm.fontManager.ttflist:
    if "Arial Unicode" in f.name:
        plt.rcParams["font.family"] = f.name; break
plt.rcParams["axes.unicode_minus"] = False

OM = 18453.1                       # ωγ (무차원, γ*=1) = 비드 하나의 항력 소산
D = dict(tot=1.838e4, chain=-56.96, store=641.2, Kp_lock=-9.203)
J = dict(tot=7.559e4, chain=5.727e4, store=1.276e5, Kp_lock=1.128e5)

fig, ax = plt.subplots(1, 3, figsize=(15.5, 5.2))

# ── ① 소산 분해 (누적 막대) ──
a = ax[0]
labs = ["비드 혼자\n(사슬 없음)", "DLVO 사슬", "JKR 사슬"]
drag = [OM, OM, OM]
chain = [0.0, D["chain"], J["chain"]]
x = np.arange(3)
a.bar(x, drag, 0.55, color="#9ecae1", label="구동 비드 자신의 용매 항력 (ωγ)")
a.bar(x, chain, 0.55, bottom=drag, color="#d62728", label="사슬이 담당하는 소산")
for i, (dr, ch) in enumerate(zip(drag, chain)):
    tot = dr + ch
    a.text(i, tot * 1.04, f"{tot:.3g}\n({tot/OM:.2f}×)", ha="center", fontsize=9)
a.set_ylabel("K″ [kT/d²]  (한 주기 소산 ∝ 이 값)")
a.set_ylim(0, 9.5e4)
a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
a.set_title("★ ① 계 전체의 소산 분해\n"
            "DLVO 사슬을 붙여도 소산이 **비드 혼자일 때와 같다** (0.996배)\n"
            "JKR 은 4.10배 — 사슬이 소산의 75.8% 담당", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── ② 두 경로 교차검증 ──
a = ax[1]
w = 0.35; x = np.arange(2)
tot = [D["tot"], J["tot"]]
recon = [D["chain"] + OM, J["chain"] + OM]
a.bar(x - w/2, tot, w, color="#2ca02c", label="① 에너지 적분  ∮F·dy/(π|ŷ|²)")
a.bar(x + w/2, recon, w, color="#ff7f0e", label="② 락인 K″_chain + ωγ")
for xi, v in zip(x - w/2, tot): a.text(xi, v*1.03, f"{v:.4g}", ha="center", fontsize=8)
for xi, v in zip(x + w/2, recon): a.text(xi, v*1.03, f"{v:.4g}", ha="center", fontsize=8)
for i, (t_, r_) in enumerate(zip(tot, recon)):
    a.text(i, max(t_, r_)*1.14, f"{100*(t_-r_)/t_:+.2f}%", ha="center", fontsize=10,
           color="darkgreen", weight="bold")
a.set_xticks(x); a.set_xticklabels(["DLVO 사슬", "JKR 사슬"], fontsize=9)
a.set_ylabel("K″ [kT/d²]"); a.set_ylim(0, 9.5e4)
a.set_title("② 두 경로가 자유 파라미터 없이 일치\n"
            "계 전체 에너지 적분 ↔ 구동점 락인 + 항력\n"
            "DLVO −0.11% · JKR −0.17%", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── ③ 저장: 계 전체 vs 구동점 ──
a = ax[2]
x = np.arange(2)
sys_ = [D["store"], J["store"]]
loc = [abs(D["Kp_lock"]), abs(J["Kp_lock"])]
a.bar(x - w/2, sys_, w, color="#6a51a3", label="계 전체 저장 K′_sys (U 의 2ω)")
a.bar(x + w/2, loc, w, color="#1f77b4", label="구동점 |K′| (락인)")
for xi, v in zip(x - w/2, sys_): a.text(xi, v*1.5, f"{v:.4g}", ha="center", fontsize=8)
for xi, v in zip(x + w/2, loc): a.text(xi, v*1.5, f"{v:.3g}", ha="center", fontsize=8)
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(["DLVO 사슬", "JKR 사슬"], fontsize=9)
a.set_ylabel("K′ [kT/d²]  (로그)"); a.set_ylim(1, 1e6)
a.set_title("③ 저장은 계 전체 ≠ 구동점\n"
            "JKR 1.13배 (굽힘이 사슬 전체에 저장)\n"
            "★ DLVO 는 구동점 K′≈0 인데도 641 저장 — 결합 신축", fontsize=10)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y", which="both")

fig.suptitle("시스템 수준 유변학 — 계 전체의 에너지 수지 (개별 비드 피팅 없음)\n"
             "n=9, ω=3000 rad/s, a=1d, k_t×100, 시드 6개", fontsize=12)
fig.tight_layout()
out = "/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/system_rheology.png"
fig.savefig(out, dpi=140); print("saved", out)
