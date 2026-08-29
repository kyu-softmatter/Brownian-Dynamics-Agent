"""`chain-bend-2d-dlvo` 첫 런 결과 확인 — 그래프 + 애니메이션.

**결과를 말로만 보고하지 않는다** (작업 관행). G1 가설(직선사슬 + 순수 중심력 DLVO +
자연장 평형 ⇒ 선형 굽힘강성이 정확히 0)을 이 런이 지지하는지 눈으로 확인한다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/viz_chain_bend_dlvo.py [run_id]
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "scratch" / "_viz"
OUT.mkdir(parents=True, exist_ok=True)

for _f in ("Arial Unicode MS", "Apple SD Gothic Neo", "AppleGothic", "NanumGothic"):
    if _f in {f.name for f in matplotlib.font_manager.fontManager.ttflist}:
        plt.rcParams["font.family"] = _f
        break
plt.rcParams["axes.unicode_minus"] = False

C_MEAS, C_ZERO, C_GOOD, C_BAD = "#1f77b4", "#d62728", "#2ca02c", "#ff7f0e"

RUN_ID = sys.argv[1] if len(sys.argv) > 1 else "chain-bend-2d-dlvo__n9-w1000-a100__cbe816dbec24"
RUN = ROOT / "runs" / RUN_ID


def load():
    m = json.loads((RUN / "metrics.json").read_text())
    npz = np.load(RUN / "observables.npz", allow_pickle=True)
    return m, npz


def make_figure(m, npz):
    n = int(m["physical"]["n_beads"])
    kp = next(o for o in m["observables"] if o["name"] == "K_prime")
    shape_mean = np.array(m["result"]["shape_profile_mean"])
    t, yb, yg = npz["t"], npz["y_bead"], npz["y_ghost"]
    shape_y = npz["shape_y"]

    fig, ax = plt.subplots(1, 3, figsize=(15, 4.4))
    fig.suptitle(f"chain-bend-2d-dlvo — 첫 런 (n={n}, DLVO 중심력만, 굽힘항 없음)  "
                 f"run_id={RUN_ID}", fontsize=11, y=1.03)

    # ① 평균 형태 프로파일 — 매끈한 곡률(포물선형) vs 국소 꺾임(삼각형형)
    a = ax[0]
    idx = np.arange(n)
    a.plot(idx, shape_mean * 1e3, "o-", color=C_MEAS, label="측정 <y_i> (시간평균)")
    # 매끈한 빔형 곡률이라면 이차함수로 잘 맞아야 한다 — 참고선으로 겹친다
    coef = np.polyfit(idx, shape_mean, 2)
    a.plot(idx, np.polyval(coef, idx) * 1e3, "--", color=C_GOOD, alpha=0.7,
           label="2차 다항 최소적합 (매끈하면 이 근처)")
    a.set_xlabel("비드 인덱스 i")
    a.set_ylabel("<y_i>  [d, ×1e-3]")
    a.set_title(f"형태 프로파일 (shape_localization={m['result']['shape_localization']:.3f}"
                f"  ·  1=매끈, ≫1=국소꺾임)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # ② 구동 vs 응답 시계열 (마지막 몇 주기)
    a = ax[1]
    T_period = 2 * np.pi / m["physical"]["omega_star"] / (t[1] - t[0]) * (t[1] - t[0])
    mask = t > t.max() - 6 * (2 * np.pi / (m["physical"]["omega_star"] / m["numerics"].get("dt_star", 1)))
    # 안전하게: 그냥 마지막 400 샘플만
    sl = slice(max(0, len(t) - 400), len(t))
    a.plot(t[sl] * 1e3, yg[sl], color="#888", lw=1, label="구동 (유령, y_ghost)")
    a.plot(t[sl] * 1e3, yb[sl], color=C_MEAS, lw=1.3, label="중앙 비드 응답 (y_bead)")
    a.set_xlabel("t  [ms]")
    a.set_ylabel("y  [d]")
    a.set_title("구동 vs 중앙 비드 응답 (마지막 400 샘플)")
    a.legend(fontsize=8)
    a.grid(alpha=0.3)

    # ③ G1 검사: K' 측정 vs 예측 0 (σ 기준)
    a = ax[2]
    meas, sig = kp["measured"], kp["sigma"]
    z = kp["err_sigma"]
    a.errorbar([0], [meas], yerr=[sig], fmt="o", color=C_MEAS, capsize=6, ms=9,
               label=f"측정 K′ = {meas:.1f} ± {sig:.1f} kT/d²")
    a.axhline(0, color=C_ZERO, ls="--", label="예측 (G1: 선형 굽힘강성 = 0)")
    ok = abs(z) < (kp.get("tol_sigma") or 3.0)
    a.set_title(f"G1 판정: {'✓ 일치' if ok else '✗ 불일치'} ({z:+.2f}σ, 기준 ±3σ)",
               color=C_GOOD if ok else C_BAD)
    a.set_xlim(-1, 1)
    a.set_xticks([])
    a.set_ylabel("K′  [kT/d²]")
    a.legend(fontsize=8, loc="upper right")
    a.grid(alpha=0.3)

    fig.tight_layout()
    p = OUT / "chain_bend_dlvo_results.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


def make_animation(npz):
    """사슬 형태의 시간 변화 — 삼각형형 좌굴이 실제로 보이는지 눈으로.

    ★ 실제 생산 데이터에서 뽑은 것 (2000 샘플 중 처음 300개, ~3주기) — 별도 저비용
    런이 아니라 이 런의 진짜 궤적이다.
    """
    shape_y = npz["shape_y"][:300]
    t = npz["t"][:300]
    n = shape_y.shape[1]
    idx = np.arange(n)

    fig, ax = plt.subplots(figsize=(5, 4))
    ymax = np.abs(shape_y).max() * 1.15 + 1e-6
    line, = ax.plot([], [], "o-", color=C_MEAS, lw=2, ms=6)
    ax.set_xlim(-0.5, n - 0.5)
    ax.set_ylim(-ymax, ymax)
    ax.set_xlabel("비드 인덱스 i")
    ax.set_ylabel("y  [d]")
    title = ax.set_title("")
    ax.grid(alpha=0.3)

    def update(fr):
        line.set_data(idx, shape_y[fr])
        title.set_text(f"chain-bend-2d-dlvo  t={t[fr]*1e3:.3f} ms  (생산 궤적)")
        return line, title

    anim = FuncAnimation(fig, update, frames=len(t), interval=60, blit=False)
    p = OUT / "chain_bend_dlvo_motion.gif"
    anim.save(p, writer=PillowWriter(fps=16))
    plt.close(fig)
    return p


def main():
    m, npz = load()
    png = make_figure(m, npz)
    gif = make_animation(npz)
    print(f"그래프: {png.relative_to(ROOT)}")
    print(f"애니메이션: {gif.relative_to(ROOT)}  (이 런의 실제 생산 궤적, 처음 ~3주기)")
    missing = {w for w in ("missing from font",) if False}
    return 0


if __name__ == "__main__":
    sys.exit(main())
