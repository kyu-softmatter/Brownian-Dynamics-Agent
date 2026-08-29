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

# ── ① 추종률 ──
a = ax[0,0]
yd = [100*o["yD"][0] for o in R]; yj = [100*o["yJ"][0] for o in R]
a.bar(x-w/2, yd, w, color="#1f77b4", label="DLVO-only")
a.bar(x+w/2, yj, w, color="#d62728", label="DLVO+JKR")
for xi,v in zip(x-w/2,yd): a.text(xi,v+2,f"{v:.0f}%",ha="center",fontsize=9)
for xi,v in zip(x+w/2,yj): a.text(xi,v+2,f"{v:.0f}%",ha="center",fontsize=9)
a.axhline(100, color="gray", ls=":", lw=1.2)
a.set_xticks(x); a.set_xticklabels(labs, fontsize=9); a.set_ylim(0,118)
a.set_ylabel("추종률 |ŷ|/a  [%]")
a.set_title("① 구동 프로토콜별 추종률\n기본 트랩은 명령의 3.5~27%만 따라간다", fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── ② DLVO 가 0 에서 몇 σ (프로토콜 무관 지표) ──
a = ax[0,1]
z = [o["z"] for o in R]
cols = ["#ff7f0e" if v>2 else "#2ca02c" for v in z]
a.bar(x, z, 0.5, color=cols)
for xi,v in zip(x,z): a.text(xi, v+0.08, f"{v:.2f}σ", ha="center", fontsize=10)
a.axhline(2, color="crimson", ls="--", lw=1.6, label="2σ — 이 아래면 0 과 구별 안 됨")
a.set_xticks(x); a.set_xticklabels(labs, fontsize=9); a.set_ylim(0, 3.0)
a.set_ylabel("DLVO 강성이 0 에서 떨어진 정도 [σ]")
a.set_title("★ ② 추종이 좋아질수록 DLVO 는 0 에 수렴\n"
            "2.28σ → 0.84σ → 0.66σ  (추종 27% → 100% → 100%)\n"
            "기본 트랩의 '유한한 K′' 는 추종 실패가 만든 인공물이었다", fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y")

# ── ③ JKR/DLVO 비 (프로토콜 무관, 로그) ──
a = ax[1,0]
ratio = [abs(o["J"][0]/o["D"][0]) for o in R]
a.bar(x, ratio, 0.5, color="#6a51a3")
for xi,v,o in zip(x,ratio,R):
    a.text(xi, v*1.4, f"{v:,.0f}배", ha="center", fontsize=10)
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
a.set_ylim(100, 3e5)
a.set_ylabel("|JKR| / |DLVO|  (로그)")
a.set_title("③ 세 프로토콜 모두 JKR 이 압도\n"
            "400배 → 12,000배 → 2,500배\n"
            "⚠ 분모가 0 에 가까워 비 자체는 불안정 — 자릿수만 볼 것", fontsize=10.5)
a.grid(alpha=0.3, axis="y", which="both")

# ── ④ 굽음 — 판별력이 프로토콜에 따라 뒤집힌다 ──
a = ax[1,1]
bd = [o["bD"][0] for o in R]; bj = [o["bJ"][0] for o in R]
a.bar(x-w/2, bd, w, color="#1f77b4", label="DLVO-only")
a.bar(x+w/2, bj, w, color="#d62728", label="DLVO+JKR")
for xi,(d_,j_) in zip(x, zip(bd,bj)):
    a.text(xi, max(d_,j_)*1.35, f"{d_/j_:.1f}배", ha="center", fontsize=10, color="darkgreen")
a.set_yscale("log"); a.set_xticks(x); a.set_xticklabels(labs, fontsize=9)
a.set_ylabel("굽음 RMS [d]  (로그)"); a.set_ylim(5e-3, 3)
a.set_title("④ 굽음(모양)의 판별력은 프로토콜에 좌우된다\n"
            "15.4배 → 1.4배 → 1.9배\n"
            "★ 변형이 강제되면 모양은 정보가 아니다 — 힘을 봐야 한다", fontsize=10.5)
a.legend(fontsize=8); a.grid(alpha=0.3, axis="y", which="both")

fig.suptitle("3 구동 프로토콜 종합 — n=9, ω=3000 rad/s, a=1470nm(=1d), 시드 6개\n"
             "trap 2종은 K′(구동점 강성), position 은 K_transfer(전달강성) — "
             "다른 양이므로 조건 간 절대값 비교 금지, 각 조건 안의 DLVO/JKR 대비만 유효",
             fontsize=11.5)
fig.tight_layout()
out=f"{S}/three_conditions.png"
fig.savefig(out, dpi=140); print("saved", out)
