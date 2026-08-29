"""Animation -- chain-relax-2d-dlvo kink release, kT=0 (deterministic) side by side
with kT=1 (the real thermal case).

CLAUDE.md working practice: "putting kT=0 (the mode shape) next to kT>0 (how far
thermal noise buries it) shows an SNR problem faster than any number." If G1 (linear
bending stiffness = 0) is right, the kT=0 panel should barely move -- there is no
restoring force (it starts at exactly the natural length, so there is no stretching
signal either).

★ This is NOT a production measurement -- visualization only. The kT=1 trajectory
reuses the GSD of an existing production run as-is (nothing is re-run); only kT=0 is
run fresh and short (same dt, same initial condition, same forces).

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
    """Take the axis joining the two ends as x'.

    The same transform as chain_relax_2d_dlvo.bow_metrics.
    """
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
    print(f"kT=1 (the real run): loaded {len(frames_kT1)} frames")
    frames_kT0 = run_kT0(spec, n_frames_target=len(frames_kT1) - 1)
    print(f"kT=0 (deterministic, fresh short run): generated "
          f"{len(frames_kT0)} frames")

    n_show = min(len(frames_kT0), len(frames_kT1))
    body_kT0 = np.array([rotate_to_body_frame(p) for p in frames_kT0[:n_show]])
    body_kT1 = np.array([rotate_to_body_frame(p) for p in frames_kT1[:n_show]])

    period = int(spec["numerics"]["n_prod"]) // (n_show - 1)
    # dt/tau_bond = C.GATE = 0.01 exactly (a design constant, dt_from_gate in
    # build_ledger) -- so time in units of tau_bond follows from the step count
    # between frames alone, with no need to go via tau_B.
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
    print("saved:", out)


if __name__ == "__main__":
    main()
