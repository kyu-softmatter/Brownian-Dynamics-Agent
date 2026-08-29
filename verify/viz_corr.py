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

# ── 상단: 변위 상관행렬 3개 ──
mats = [("dlvo_full", "DLVO-only  (n=9)", "굽힘 없음 → 비드끼리 제각각"),
        ("jkr_full",  "DLVO+JKR  (n=9)",  "굽힘 있음 → 사슬이 통째로 움직임"),
        ("n25_full",  "DLVO-only  (n=25)", "긴 사슬에서도 상관 없음")]
for k, (key, title, sub) in enumerate(mats):
    ax = fig.add_subplot(gs[0, k])
    C = R[key]["C"]
    im = ax.imshow(C, cmap="RdBu_r", vmin=-1, vmax=1)
    off = C[np.triu_indices_from(C, k=1)]
    ax.set_title(f"{title}\n{sub}\n비대각 평균 = {off.mean():+.3f}", fontsize=9.5)
    ax.set_xlabel("비드 j"); ax.set_ylabel("비드 i")
    fig.colorbar(im, ax=ax, fraction=0.046)

# ── 하단 좌: 접선 상관 ──
ax = fig.add_subplot(gs[1, :2])
for key, col, lab in [("dlvo_full", "#1f77b4", "DLVO-only n=9"),
                      ("jkr_full", "#d62728", "DLVO+JKR n=9"),
                      ("n25_full", "#2ca02c", "DLVO-only n=25")]:
    c = R[key]["c"]
    ax.plot(np.arange(len(c)), c, "o-", color=col, label=lab, ms=5)
ax.axhline(1.0, color="gray", lw=0.8, ls=":")
ax.set_xlabel("결합 간격 m"); ax.set_ylabel("<t_i . t_{i+m}>")
ax.set_ylim(0.6, 1.03)
ax.set_title("접선-접선 상관 — 이 계에서는 판별력이 약하다\n"
             "두 사슬 다 x축을 따라 거의 곧아서(변위 ≪ 결합길이) cos θ ≈ 1 이 되어버린다.\n"
             "n=25 는 m=1 에서 0.80 으로 떨어져 굽힘 없음이 보이지만, 그 뒤 평평 — 트랩 구속 때문",
             fontsize=9.5)
ax.legend(fontsize=8); ax.grid(alpha=0.3)

# ── 하단 우: 요약 막대 ──
ax = fig.add_subplot(gs[1, 2])
labs = ["DLVO\nn=9", "JKR\nn=9", "DLVO\nn=25"]
vals = [R[k]["C"][np.triu_indices_from(R[k]["C"], k=1)].mean()
        for k in ("dlvo_full", "jkr_full", "n25_full")]
cols = ["#1f77b4", "#d62728", "#2ca02c"]
ax.bar(range(3), vals, color=cols)
ax.axhline(0, color="k", lw=1)
ax.set_xticks(range(3)); ax.set_xticklabels(labs, fontsize=9)
ax.set_ylabel("변위 상관 비대각 평균")
ax.set_ylim(-0.15, 1.0)
for i, v in enumerate(vals):
    ax.text(i, v + (0.04 if v > 0 else -0.08), f"{v:+.3f}", ha="center", fontsize=9)
ax.set_title("★ 진짜 판별량\n굽힘항 하나로 −0.02 → +0.91", fontsize=9.5)
ax.grid(alpha=0.3, axis="y")

fig.suptitle("입자 위치 상관 — 변위 상관행렬 <dy_i dy_j> 와 접선 상관 <t_i . t_j> "
             "(ω=3000·a=1d / n=25는 ω=10·a=632nm)", fontsize=11.5)
out = f"{S}/correlations.png"
fig.savefig(out, dpi=140, bbox_inches="tight"); print("saved", out)
