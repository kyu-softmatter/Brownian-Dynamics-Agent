"""애니메이션 — chain-relax-2d-dlvo 킹크 방출을 kT=0(결정론) vs kT=1(실제 열적) 나란히.

CLAUDE.md 작업 관행: "kT=0(모드 형태)과 kT>0(열요동에 묻히는 정도)을 나란히 두면 SNR
문제가 숫자보다 빨리 보인다." G1(선형 굽힘강성=0)이 맞다면 kT=0 패널은 거의 안 움직여야
한다 — 복원할 힘이 없기 때문이다(정확히 자연장에서 시작해 신장 신호도 없다).

★ 이건 생산 측정이 아니다 — 시각화 전용. kT=1 궤적은 기존 프로덕션 런의 GSD를 그대로
재사용하고(재실행 안 함), kT=0 만 새로 짧게 돌린다(같은 dt·같은 초기조건·같은 힘).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/anim_chain_relax_kink.py
"""
import json
import math
import sys
from pathlib import Path

import gsd.hoomd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from chain_bend_dlvo_2d import SIGMA_CORE_STAR, build_table_arrays  # noqa: E402
from chain_relax_2d_dlvo import kink_positions  # noqa: E402
from bdbot import sim as SIM  # noqa: E402

KINK_RUN = ROOT / "runs/chain-relax-2d-dlvo__n9-kink-a0.300__1e7b282680f2"


def rotate_to_body_frame(pos):
    """양끝을 잇는 축을 x'축으로 — chain_relax_2d_dlvo.bow_metrics 와 같은 변환."""
    d = pos[-1] - pos[0]
    L = float(np.hypot(d[0], d[1]))
    if L < 1e-9:
        return pos - pos.mean(axis=0)
    u = d / L
    rel = pos - pos[0]
    xp = rel[:, 0] * u[0] + rel[:, 1] * u[1]
    yp = -rel[:, 0] * u[1] + rel[:, 1] * u[0]
    xp -= xp.mean()
    return np.stack([xp, yp], axis=1)


def load_kT1_frames():
    spec = json.loads((KINK_RUN / "spec.json").read_text())
    with gsd.hoomd.open(str(KINK_RUN / "traj_A.gsd"), "r") as f:
        L = float(f[0].configuration.box[0])
        frames = []
        for fr in f:
            pos = np.array(fr.particles.position, dtype=float)[:, :2]
            img = np.array(fr.particles.image, dtype=float)[:, :2]
            frames.append(pos + img * L)
    return np.array(frames), spec


def run_kT0(spec, n_frames_target):
    import hoomd
    import hoomd.md as md

    P, Nm = spec["params"], spec["numerics"]
    n, kink_angle = int(P["n_beads"]), float(P["kink_angle"])
    h_min_star, dt = float(P["h_min_star"]), float(Nm["dt_star"])
    ell_star = 1.0 + h_min_star
    r_cut_star = 1.0 + float(P["cutoff_h_star"])
    r_min_star = 1.0 + 1e-6
    box_star = 4.0 * (n - 1) * ell_star
    reduced = {k: P[k] for k in ("kappa_star", "edl_amp", "vdw_amp", "a_star")}

    pos0 = kink_positions(n, ell_star, kink_angle)
    sim = SIM.make_sim(SIM.frame_2d(pos0, box_star), seed=int(Nm["seed"]))

    cell = md.nlist.Cell(buffer=0.2)
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_star, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_star, U=U_arr, F=F_arr)
    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)

    method = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[method], forces=[tab, wca])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    n_prod = int(Nm["n_prod"])
    period = max(1, n_prod // n_frames_target)
    frames = [np.array(sim.state.get_snapshot().particles.position, dtype=float)[:, :2].copy()]
    for _ in range(n_prod // period):
        sim.run(period)
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, dtype=float)[:, :2]
        img = np.array(snap.particles.image, dtype=float)[:, :2]
        frames.append(pos + img * box_star)
    return np.array(frames)


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.animation as anim
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    frames_kT1, spec = load_kT1_frames()
    print(f"kT=1 (실제 런) 프레임 {len(frames_kT1)}개 로드")
    frames_kT0 = run_kT0(spec, n_frames_target=len(frames_kT1) - 1)
    print(f"kT=0 (결정론, 신규 짧은 런) 프레임 {len(frames_kT0)}개 생성")

    n_show = min(len(frames_kT0), len(frames_kT1))
    body_kT0 = np.array([rotate_to_body_frame(p) for p in frames_kT0[:n_show]])
    body_kT1 = np.array([rotate_to_body_frame(p) for p in frames_kT1[:n_show]])

    period = int(spec["numerics"]["n_prod"]) // (n_show - 1)
    # dt/tau_bond = C.GATE = 0.01 정확히(설계 상수, build_ledger 의 dt_from_gate) —
    # tau_B 를 거칠 필요 없이 프레임 간 스텝수만으로 τ_bond 단위 시간이 나온다.
    t_per_frame_tau_bond = period * 0.01

    xmax = float(np.max(np.abs(np.concatenate([body_kT0[:, :, 0], body_kT1[:, :, 0]])))) * 1.1
    ymax = max(0.35, float(np.max(np.abs(np.concatenate([body_kT0[:, :, 1], body_kT1[:, :, 1]])))) * 1.3)

    fig, (axL, axR) = plt.subplots(1, 2, figsize=(13, 3.4))
    for ax, title in ((axL, "kT = 0 (deterministic — no thermal noise)"),
                      (axR, "kT = 1 (real production run)")):
        ax.set_xlim(-xmax, xmax); ax.set_ylim(-ymax, ymax)
        ax.set_aspect("equal"); ax.grid(alpha=.3)
        ax.set_xlabel("body-frame x [d]"); ax.set_ylabel("body-frame y [d]")
        ax.set_title(title, fontsize=10)
        ax.axhline(0, color="gray", lw=.5)

    lineL, = axL.plot([], [], "o-", color="tab:blue", ms=6, lw=1.5)
    lineR, = axR.plot([], [], "o-", color="tab:red", ms=6, lw=1.5)
    txtL = axL.text(0.02, 0.95, "", transform=axL.transAxes, va="top", fontsize=9)
    txtR = axR.text(0.02, 0.95, "", transform=axR.transAxes, va="top", fontsize=9)
    fig.suptitle("chain-relax-2d-dlvo — kink release (n=9, angle=0.30 rad)\n"
                "G1 check: no linear bending stiffness -> kT=0 panel should barely move",
                fontsize=11, y=0.98)
    fig.subplots_adjust(top=0.72, bottom=0.15, left=0.06, right=0.98, wspace=0.18)

    def update(i):
        lineL.set_data(body_kT0[i, :, 0], body_kT0[i, :, 1])
        lineR.set_data(body_kT1[i, :, 0], body_kT1[i, :, 1])
        t = i * t_per_frame_tau_bond
        txtL.set_text(f"t = {t:,.0f} $\\tau_{{bond}}$\nbow_rms = {np.sqrt(np.mean(body_kT0[i,:,1]**2)):.4f} d")
        txtR.set_text(f"t = {t:,.0f} $\\tau_{{bond}}$\nbow_rms = {np.sqrt(np.mean(body_kT1[i,:,1]**2)):.4f} d")
        return lineL, lineR, txtL, txtR

    ani = anim.FuncAnimation(fig, update, frames=n_show, interval=80, blit=False)
    out = ROOT / "runs" / KINK_RUN.name / "kink_release_kT0_vs_kT1.gif"
    ani.save(out, writer="pillow", fps=12)
    plt.close(fig)
    print("저장:", out)


if __name__ == "__main__":
    main()
