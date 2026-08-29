"""`trap-drag` 궤적 → 애니메이션 GIF. 결함을 색으로 구분해 전위 생성·소멸을 보게.

에너지 그래프(평탄역 + 급락)와 결함 곡선이 **같은 모양**이라는 것이 이 케이스의 결론인데,
숫자로만 보면 "왜 그런지"가 안 보입니다. 탐침이 격자를 밀고 지나가며 5–7 쌍을 만들고,
멈춘 뒤 그것들이 서로 만나 사라지는 과정을 눈으로 보는 것이 목적입니다.

색: 회색 = 6배위(정상) · 주황 = 5배위 · 보라 = 7배위 · 빨강 = 그 외
    파랑 테두리 = 끌리는 탐침 · × = 트랩 중심

    $PY scratch/trap_drag_movie.py runs/<run_id> [--stride 2] [--fps 20]
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from bdbot import nondim as ND  # noqa: E402

COLOR = {5: "tab:orange", 6: "0.72", 7: "tab:purple"}
PH = ("예열", "평형", "끌기", "이완")


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("rundir")
    ap.add_argument("--stride", type=int, default=2, help="프레임 건너뛰기")
    ap.add_argument("--fps", type=int, default=18)
    ap.add_argument("--out", default="trajectory.gif")
    # ★ 빠른 v 는 끌기가 전체의 4% 뿐이라(v=32: 5667프레임 중 209개) 전체를 균등
    #   추출하면 정작 볼 구간이 스쳐 지나갑니다. 프레임 범위를 잘라 씁니다.
    ap.add_argument("--frames", default=None, help="시작:끝 (프레임 인덱스)")
    args = ap.parse_args()

    import freud
    import gsd.hoomd
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.animation import FuncAnimation, PillowWriter
    from matplotlib.collections import EllipseCollection
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    d = Path(args.rundir)
    spec = ND.load(d / "spec.json")
    P, Nm = spec.params, spec.numerics
    n_x, n_y, N = int(P["n_x"]), int(P["n_y"]), int(P["N"])
    Lx, Ly, a_nn = P["Lx_star"], P["Ly_star"], P["a_nn_star"]
    v_star, dt = P["v_star"], float(Nm["dt_star"])
    n_w, n_e, n_d = int(Nm["n_warm"]), int(Nm["n_equil"]), int(Nm["n_prod"])
    tau_int = spec.reduced("times", "tau_int")

    import trap_drag_2d as TD
    pos0 = TD.hex_lattice(n_x, n_y, a_nn)
    probe = (n_y // 2) * n_x + n_x // 2
    T1, T2 = (n_w + n_e) * dt, (n_w + n_e + n_d) * dt

    def phase_of(step):
        return PH[0] if step < n_w else PH[1] if step < n_w + n_e else \
            PH[2] if step < n_w + n_e + n_d else PH[3]

    traj = gsd.hoomd.open(str(d / "traj_A.gsd"), "r")
    lo, hi = 0, len(traj)
    if args.frames:
        a_, b_ = args.frames.split(":"); lo, hi = int(a_ or 0), int(b_ or len(traj))
    frames = list(range(lo, min(hi, len(traj)), max(1, args.stride)))
    box = freud.box.Box(Lx=Lx, Ly=Ly, is2D=True)
    voro = freud.locality.Voronoi()

    # 에너지 시계열 (아래 띠)
    z = np.load(d / "observables.npz")
    res = json.loads((d / "metrics.json").read_text())["result"]
    t_u = z["_t_step"] * dt / tau_int
    u = z["u_pair"]
    u_eq = res["U_pair_eq"]

    fig = plt.figure(figsize=(8.4, 7.4), dpi=85)
    gs = fig.add_gridspec(2, 1, height_ratios=[3.1, 1.0], hspace=0.28)
    ax, axu = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])

    ax.set(xlim=(-Lx / 2, Lx / 2), ylim=(-Ly / 2, Ly / 2), aspect="equal",
           xlabel="x / d", ylabel="y / d")
    ec = EllipseCollection(np.ones(N), np.ones(N), np.zeros(N), units="xy",
                           offsets=np.zeros((N, 2)), offset_transform=ax.transData,
                           facecolors=["0.72"] * N, edgecolors="0.35", linewidths=.4)
    ax.add_collection(ec)
    ring = ax.scatter([], [], s=140, facecolors="none", edgecolors="tab:blue", lw=2.2,
                      zorder=5)
    cross, = ax.plot([], [], "x", color="k", ms=11, mew=2.5, zorder=6)
    ttl = ax.set_title("")

    axu.plot(t_u, u, "-", lw=.7, color="0.6")
    axu.axhline(u_eq, color="green", ls="--", lw=1.1, label="평형")
    mark, = axu.plot([], [], "o", color="crimson", ms=7, zorder=5)
    vline = axu.axvline(t_u[0], color="crimson", lw=1)
    axu.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$\langle U\rangle/N$ [kT]",
            xlim=(t_u.min(), t_u.max()))
    axu.legend(fontsize=8, loc="upper left"); axu.grid(alpha=.3)

    from matplotlib.lines import Line2D
    ax.legend(handles=[Line2D([], [], marker="o", ls="", mfc=COLOR[6], mec="0.35",
                              label="6배위 (정상)"),
                       Line2D([], [], marker="o", ls="", mfc=COLOR[5], mec="0.35",
                              label="5배위"),
                       Line2D([], [], marker="o", ls="", mfc=COLOR[7], mec="0.35",
                              label="7배위"),
                       Line2D([], [], marker="o", ls="", mfc="none", mec="tab:blue",
                              mew=2, label="탐침"),
                       Line2D([], [], marker="x", ls="", color="k", label="트랩 중심")],
              fontsize=8, loc="upper right", ncol=2, framealpha=.92)

    def draw(fi):
        fr = traj[fi]
        p = np.array(fr.particles.position, dtype=float)
        step = fr.configuration.step
        vn = voro.compute((box, p)).nlist
        zc = np.asarray(vn.neighbor_counts)
        ec.set_offsets(p[:, :2])
        ec.set_facecolors([COLOR.get(int(c), "tab:red") for c in zc])
        ring.set_offsets(p[probe:probe + 1, :2])
        # 트랩 중심 (조각별 선형) — 박스 밖으로 나가면 되감아 그린다
        cx = pos0[probe, 0] + v_star * (min(max(step * dt, T1), T2) - T1)
        cx = (cx + Lx / 2) % Lx - Lx / 2
        cross.set_data([cx], [pos0[probe, 1]])
        tt = step * dt / tau_int
        nd = int((zc != 6).sum())
        ttl.set_text(f"[{phase_of(step)}]   t = {tt:6.1f} $\\tau_{{int}}$   "
                     f"결함 {nd:2d}개   (5배위 {int((zc==5).sum())} · "
                     f"7배위 {int((zc==7).sum())})")
        j = int(np.searchsorted(t_u, tt))
        if 0 < j < len(t_u):
            mark.set_data([t_u[j]], [u[j]])
        vline.set_xdata([tt, tt])
        return ec, ring, cross, ttl, mark, vline

    fig.suptitle(f"{spec.label}   Γ={P['Gamma']:.1f}  φ={P['phi']}  "
                 f"{n_x}×{n_y}={N}   끌기 {n_d*dt*v_star/a_nn:.1f}격자", fontsize=11)
    anim = FuncAnimation(fig, draw, frames=frames, blit=False)
    out = d / args.out
    anim.save(str(out), writer=PillowWriter(fps=args.fps))
    plt.close(fig)
    print(f"{out}  ({out.stat().st_size/1e6:.1f} MB, {len(frames)}프레임)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
