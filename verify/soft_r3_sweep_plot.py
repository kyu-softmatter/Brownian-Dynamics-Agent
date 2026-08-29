"""Phase 1-B 스윕 요약 — 스케치의 질문("final configuration? rdf / voronoi / structure")에 답한다.

    $PY scratch/soft_r3_sweep_plot.py
    → runs/soft-r3-2d-A-sweep__SUMMARY/{summary.png, summary.md}
"""
import glob
import json
import math
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False
matplotlib.rcParams["mathtext.fontset"] = "dejavusans"   # 한글 폰트에 −·µ 글리프가 없음

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "runs" / "soft-r3-2d-A-sweep__SUMMARY"
OUT.mkdir(parents=True, exist_ok=True)
HEX_NN = math.sqrt(2 / math.sqrt(3))

runs = []
for tag in ("A0.1", "A1", "A10", "A100"):
    d = glob.glob(str(ROOT / f"runs/soft-r3-2d-A-sweep__{tag}__*"))
    assert len(d) == 1, (tag, d)
    m = json.load(open(d[0] + "/metrics.json"))
    z = np.load(d[0] + "/observables.npz")
    runs.append(dict(tag=tag, m=m, z=z,
                     A=m["dimensionless"]["A      = U(d)/kT        접촉 결합"],
                     G=m["structure"]["Gamma"]))

a_star = math.sqrt(math.pi / (4 * runs[0]["m"]["physical"]["phi"]))
colors = plt.cm.viridis(np.linspace(0.05, 0.9, len(runs)))

fig = plt.figure(figsize=(14, 9))
gs = fig.add_gridspec(3, 4, height_ratios=[1.15, 1.0, 0.95], hspace=0.55, wspace=0.32)

# ── ① g(r) 전체 ────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[0, :2])
for r_, c in zip(runs, colors):
    ax.plot(r_["z"]["rdf_r"], r_["z"]["rdf_g"], lw=1.5, color=c,
            label=f"A={r_['A']:g}  Γ={r_['G']:.2f}")
ax.axvline(HEX_NN * a_star, ls=":", c="k", alpha=.6)
ax.text(HEX_NN * a_star, ax.get_ylim()[1] * .93, " $a_{NN}$(육방)", fontsize=8)
ax.axhline(1, c="k", lw=.5, alpha=.4)
ax.set(xlabel="r / d", ylabel="g(r)", xlim=(0.8, 6.5),
       title="① 동경분포함수 — 결합세기 Γ 를 올리면 결정 피크가 선다")
ax.legend(fontsize=8); ax.grid(alpha=.3)

# ── ② ψ₆ 와 6배위 비율 vs Γ ────────────────────────────────────────────
ax = fig.add_subplot(gs[0, 2:])
G = [r_["G"] for r_ in runs]
p6 = [r_["m"]["structure"]["psi6"] for r_ in runs]
p6e = [r_["m"]["structure"]["psi6_sem"] for r_ in runs]
c6 = [r_["m"]["structure"]["coord_hist"][6] for r_ in runs]
ax.errorbar(G, p6, yerr=p6e, marker="o", lw=1.6, capsize=3, label=r"$|\psi_6|$")
ax.plot(G, c6, marker="s", lw=1.6, ls="--", label="Voronoi 6배위 비율")
ax.set(xscale="log", xlabel=r"$\Gamma = U(a_{mean})/k_BT$", ylabel="질서 지표",
       ylim=(0, 1.05), title="② 육방 질서 — 전이는 Γ 3~30 사이")
ax.axvspan(3, 30, alpha=.12, color="tab:red")
ax.legend(fontsize=8); ax.grid(alpha=.3, which="both")
for g_, y_ in zip(G, p6):
    ax.annotate(f"{y_:.3f}", (g_, y_), textcoords="offset points", xytext=(4, 6), fontsize=7)

# ── ③ 최종 배치 4장 ────────────────────────────────────────────────────
for i, (r_, c) in enumerate(zip(runs, colors)):
    ax = fig.add_subplot(gs[1, i])
    xy = r_["z"]["final_xy"]
    L = r_["m"]["dimensionless"]["L/d                     박스 크기"]
    ax.plot(xy[:, 0], xy[:, 1], "o", ms=3.0, color=c, mec="none")
    ax.set(xlim=(-L / 2, L / 2), ylim=(-L / 2, L / 2), aspect="equal",
           title=f"A={r_['A']:g}  Γ={r_['G']:.2f}", xticks=[], yticks=[])
    ax.set_xlabel(f"$\\psi_6$={r_['m']['structure']['psi6']:.3f}", fontsize=9)
fig.text(0.5, 0.655, "③ 최종 배치 (phi=0.35, N=400, T_obs=100 tau_B)",
         ha="center", fontsize=10.5)

# ── ④ Voronoi 배위수 분포 ──────────────────────────────────────────────
ax = fig.add_subplot(gs[2, :2])
w = 0.2
for j, (r_, c) in enumerate(zip(runs, colors)):
    ch = np.array(r_["m"]["structure"]["coord_hist"])
    ks = np.arange(3, 11)
    ax.bar(ks + (j - 1.5) * w, ch[3:11], width=w, color=c, label=f"Γ={r_['G']:.2f}")
ax.set(xlabel="Voronoi 배위수", ylabel="분율", title="④ 배위수 분포 — 결정은 6에 몰린다")
ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")

# ── ⑤ 수렴 확인 ────────────────────────────────────────────────────────
ax = fig.add_subplot(gs[2, 2:])
base = json.load(open(glob.glob(str(ROOT / "runs/*A100__*/metrics.json"))[0]))
labels, dt_ch, rc_ch = [], [], []
for key, get in (("⟨U⟩/N", lambda mm: [x for x in mm["observables"]
                                       if x["name"].startswith("에너지")][0]["measured"]),
                 (r"$\psi_6$", lambda mm: mm["structure"]["psi6"]),
                 ("NN 거리", lambda mm: mm["structure"]["nn_distance_d"]),
                 ("6배위", lambda mm: mm["structure"]["coord_hist"][6])):
    labels.append(key)
    for tag, acc in (("A100-dt0.5", dt_ch), ("A100-rc7", rc_ch)):
        c = json.load(open(glob.glob(str(ROOT / f"runs/*{tag}*/metrics.json"))[0]))
        acc.append(100 * (get(c) - get(base)) / abs(get(base)))
x = np.arange(len(labels))
ax.bar(x - 0.19, dt_ch, 0.36, label="CV1  dt 절반")
ax.bar(x + 0.19, rc_ch, 0.36, label=r"CV2  $r_c$ 5a→7a")
ax.axhline(0, c="k", lw=.6)
ax.set(xticks=x, xticklabels=labels, ylabel="변화 [%]",
       title="⑤ 수렴 — 구조는 불변, 절대 에너지만 $r_c$ 의존")
ax.legend(fontsize=8); ax.grid(alpha=.3, axis="y")
for xi, v in zip(x - 0.19, dt_ch):
    ax.annotate(f"{v:+.3f}", (xi, v), textcoords="offset points", xytext=(0, 3 if v >= 0 else -11),
                ha="center", fontsize=7)
for xi, v in zip(x + 0.19, rc_ch):
    ax.annotate(f"{v:+.2f}", (xi, v), textcoords="offset points", xytext=(0, 3 if v >= 0 else -11),
                ha="center", fontsize=7)

fig.suptitle("Phase 1-B  soft-r3-2d-A-sweep  —  $U/k_BT = A(d/r)^3$ + WCA 코어, "
             "d=5 um 실리카 / 물 300 K, phi=0.35", fontsize=12.5)
fig.savefig(OUT / "summary.png", dpi=140, bbox_inches="tight")
plt.close(fig)

# ── 요약 표 ────────────────────────────────────────────────────────────
lines = ["# Phase 1-B `soft-r3-2d-A-sweep` — 스윕 결과", "",
         "물리계: `U/kT = A(d/r)³` + WCA(σ=d, ε=kT), 2D 주기, d=5 µm, 물@300 K,",
         "φ=0.35 (a_mean=1.498 d), N=400, r_c=5 a_mean, T_obs=100 τ_B (τ_B=242.05 s)", "",
         "| A | Γ=U(a_mean)/kT | ψ₆ | 6배위 | NN/d | σ_NN/NN | 상태 | dt/τ_B | 벽시계 |",
         "|---|---|---|---|---|---|---|---|---|"]
for r_ in runs:
    st = r_["m"]["structure"]
    lines.append(f"| {r_['A']:g} | {r_['G']:.4f} | {st['psi6']:.4f} | "
                 f"{st['coord_hist'][6]:.3f} | {st['nn_distance_d']:.3f} | "
                 f"{st['nn_std_rel']:.4f} | {'육방 결정' if st['psi6'] > 0.6 else '유체'} | "
                 f"{r_['m']['numerics']['dt_over_tau_B']:.2e} | "
                 f"{r_['m']['wall_seconds']/60:.0f}분 |")
lines += ["", "## 답 — 스케치의 'final configuration?'", "",
          "- **Γ ≲ 3 유체, Γ ≈ 30 육방 결정.** 전이는 Γ 3~30 사이 (더 좁히려면 A 스윕을 촘촘히).",
          "- **A=0.1 과 A=1 은 구별되지 않는다** (ψ₆ 둘 다 0.4347). φ=0.35에서 4자릿수 A 스윕이",
          "  만드는 상태는 3개뿐 — 약결합 쪽은 WCA 유체이고 r⁻³는 무의미하다.",
          "- 결정의 최근접거리는 육방 예측 `a_NN = √(2/√3)·a_mean` 과 **+0.45%** 일치.", "",
          "## 검증", "",
          "| 검증 | 결과 |", "|---|---|",
          "| 2입자 직접 대조 (`scratch/verify_pair_table.py`) | 0.000% |",
          "| 에너지 일관성 `⟨U⟩/N` vs `(ρ/2)∫U g(r) 2πr dr` | +0.00 ~ +0.67% (7런) |",
          "| 육방 NN 거리 vs 파라미터 없는 예측 | +0.45% |",
          "| 희박극한 `g(r)` vs `e^{-βU}[1+ρ∫f f]` | RMS 2.43% (0차만 쓰면 6.30%) |",
          "| CV1 dt 절반 | `⟨U⟩/N` −0.004%, ψ₆ +0.13% |",
          "| CV2 `r_c` 5a→7a | 구조 0.15% 이내 불변, 절대 `⟨U⟩/N` **+7.5%** |", "",
          "## 확인하지 않은 것", "",
          "- `A`의 무차원 해석은 스케치 목표에서 역추론한 것 (`system.yaml.not_verified`)",
          "- 유한크기: N=400 하나만 돌렸다. N=900 대조 미실시 (CV3 미완)",
          "- 문헌 대조 없음 (KB 문헌 0편). 2D r⁻³ 녹음 전이의 문헌 Γ 값 미확인",
          "- 초기조건 의존성: RSA 랜덤 배치 하나만. 결정에서 시작한 런과 대조 안 함",
          ]
(OUT / "summary.md").write_text("\n".join(lines))
print(f"→ {OUT.relative_to(ROOT)}/summary.png")
print(f"→ {OUT.relative_to(ROOT)}/summary.md")
