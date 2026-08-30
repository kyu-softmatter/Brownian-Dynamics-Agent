"""Exploration -- drive the centre bead of the chain-relax-2d-dlvo chain directly as
y(t) = A sin(omega*t).

Takes chain-relax-2d-dlvo (no friction, DLVO attraction only, no explicit bending,
**and both ends free with no traps either**), adds a single local drive, and asks
"how far does that drive propagate along the chain?".
The forcing is the same as chain-bend-2d-dlvo's `--drive-mode position` (the bead
position is written directly rather than moving a trap, reusing
`_move_ghost_action`), but **the ends are NOT held by traps** -- so no trap
compliance absorbs the displacement. By G1 (central forces + natural bond length =>
transverse bond energy is O(y^4)), at small amplitude the neighbouring bonds stretch
only as O(y^2) -- which predicts that **at small amplitude the drive should barely
propagate at all** (a prior prediction, checkable with `--amp`).

★★ Amplitude safety limit -- with no trap compliance, 100% of it goes into bond
   stretching:
    induced stretch ~ amp^2 / (2*ell)   (geometry, O(y^2))
    it must stay inside the F_max margin (h_min -> h_min+3.46nm) or the bond fails
    amp=50nm -> stretch 0.84nm (safe) /
    amp=632nm (the chain-bend-2d-dlvo convention) -> 134.8nm (fails)
-> the default amplitude was lowered to **50nm**. A larger amplitude is available via
--amp, but it accepts the possibility of bond failure -- which is itself a result.

★ This is a visualization/exploration script -- it does NOT pass through the L3 spec
or the health gates (it is not a production run).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/drive_chain_relax_center.py --omega 3000 --amp 50 --cycles 20
"""
import argparse
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from chain_bend_dlvo_2d import (  # noqa: E402
    SIGMA_CORE_STAR, build_table_arrays, _move_ghost_action, find_well, dlvo_reduced_params,
    F_h_star,
)
from chain_relax_2d_dlvo import kink_positions, bend_angles, bow_metrics  # noqa: E402
from bdbot import materials as M, sim as SIM  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from bdbot.units import Q  # noqa: E402
import yaml  # noqa: E402


def load_physics():
    raw = yaml.safe_load((ROOT / "intake/chain-bend-2d-dlvo/system.yaml").read_text())
    P = load_node
    sys_ = {"d": P(raw["particle"]["diameter"]), "T": P(raw["medium"]["temperature"]),
           "eta": P(raw["medium"]["viscosity"]), "rho_p": P(raw["particle"]["density"]),
           "eps_r": P(raw["medium"]["relative_permittivity"]), "psi0": P(raw["particle"]["surface_potential"]),
           "ionic_strength": P(raw["medium"]["ionic_strength"]),
           "A_H": P(raw["interactions"][0]["hamaker_constant"])}
    return sys_


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=9)
    ap.add_argument("--omega", type=float, default=3000.0, help="rad/s")
    ap.add_argument("--amp", type=float, default=50.0, help="nm")
    ap.add_argument("--cycles", type=float, default=20.0)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eq-tau-bond", type=float, default=200.0)
    ap.add_argument("--kT", type=float, default=1.0, help="0.0 = deterministic (OverdampedViscous), "
                    "1.0 = really thermal (Brownian) -- run twice to compare kT=0 "
                    "against kT=1 side by side")
    ap.add_argument("--vdw-amp", type=float, default=None,
                    help="rescale the vdW amplitude (reduced, = A_H/(12kT)). Default "
                         "is the literature value (0.2113). At 0.2000, "
                         "U(h*=0.1)=-1.00kT (barrier and secondary minimum confirmed "
                         "to survive: barrier=427.6kT, well=-10.98kT). "
                         "edl_amp and kappa_star (ionic strength) are left alone")
    ap.add_argument("--r-cut", type=float, default=None,
                    help="rescale the table cutoff r*=r/d (default 1.06). Use 2.0 to "
                         "reach out to h=1d")
    ap.add_argument("--r-min", type=float, default=None,
                    help="table lower bound r* (default 1.000001; trap 11: F=0 for "
                         "r<r_min). ★ Required when used with a barrier-free vdw_amp "
                         "(e.g. 2.0) -- the force there diverges as h->0, so r_min "
                         "must sit slightly outside physical contact (1.0) for dt to "
                         "stay finite (e.g. 1.01)")
    args = ap.parse_args()

    sys_ = load_physics()
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    gamma, tau_B = b["gamma"], b["tau_B"]
    p = dlvo_reduced_params(sys_)
    if args.vdw_amp is not None:
        p = dict(p, vdw_amp=args.vdw_amp)
    w = find_well(p)
    print(f"DLVO: vdw_amp={p['vdw_amp']:.4f}  barrier_h*={w['barrier_h']:.5f} "
         f"({w['barrier_U']:.1f}kT)  h_min*={w['h_min']:.5f}  well={w['U_min']:.2f}kT  "
         f"k_bond*={w['k_bond_star']:.4e}")

    no_barrier_mode = args.r_min is not None
    if no_barrier_mode:
        # ★★ With no barrier (and hence no secondary minimum), find_well()'s
        #   k_bond_star is meaningless -- dt is instead back-computed from the maximum
        #   force at the table's lower bound (r_min), with a 1% safety margin, in the
        #   same spirit as bd-hoomd Guard's step_disp_max convention. r_min sits
        #   slightly outside physical contact (r*=1) so the true h->0 divergence is
        #   never met; that region is covered instead by the "dead zone" of the table
        #   plus WCA (the table gives F=0 for r<r_min, and WCA only repels for r<1).
        r_min_star = args.r_min
        F_at_rmin = abs(float(F_h_star(r_min_star - 1.0, p)))
        dt = 0.01 / F_at_rmin
        ell_star = r_min_star
        print(f"barrier-free mode -- at r_min*={r_min_star:.4f}, "
              f"|F|={F_at_rmin:.2f} kT/d, "
              f"dt(1%% margin)={dt:.4e} tau_B")
    else:
        h_min_star = w["h_min"]
        ell_star = 1.0 + h_min_star
        k_bond = Q(w["k_bond_star"], "dimensionless") * b["kT"] / d ** 2
        tau_bond = (gamma / k_bond).to("s")
        dt_star = 1e-2                                  # dt/tau_bond, same as this
                                                        # system's gate
        dt = dt_star * float((tau_bond / tau_B).to("dimensionless").magnitude)  # dt*=dt/tau_B

    omega_star = float((Q(args.omega, "1/s") * tau_B).to("dimensionless").magnitude)
    amp_star = float((Q(args.amp, "nm") / d).to("dimensionless").magnitude)

    tau_period_star = 2 * math.pi / omega_star     # the period, in tau_B
    period_steps = int(round(tau_period_star / dt))
    n_cycles = args.cycles
    n_prod = int(round(n_cycles * tau_period_star / dt))
    frames_per_cycle = 24
    n_frames = int(round(n_cycles * frames_per_cycle))
    capture_every = max(1, n_prod // n_frames)
    # Keep the ZOH ratio (bd-hoomd trap 17) UPDATE_EVERY/period_steps below 1% -- dt
    # varies enormously between systems (barrier-free mode is 100x larger), so a fixed
    # 50 flattens the sine into a staircase when period_steps is small (measured: at
    # period=332, a value of 50 updates only 6-7 times per cycle).
    UPDATE_EVERY = max(1, min(50, period_steps // 100))

    n, mid = args.n, args.n // 2
    print(f"n={n} mid={mid}  omega={args.omega:.0f} rad/s (omega*={omega_star:.4e})  "
         f"amp={args.amp:.0f} nm (amp*={amp_star:.4e})")
    print(f"dt*={dt:.4e}  period={period_steps:,} steps  n_prod={n_prod:,} steps "
         f"({n_cycles:g} cycles)  capture every {capture_every} steps -> {n_prod//capture_every} frames")

    import hoomd
    import hoomd.md as md

    r_cut_star = args.r_cut if args.r_cut is not None else 1.0 + 0.06
    r_min_used = args.r_min if args.r_min is not None else 1.0 + 1e-6
    box_star = max(4.0 * (n - 1) * ell_star, 4.0 * r_cut_star)
    pos0 = kink_positions(n, ell_star, 0.0)             # straight, starting at the
                                                        # natural length (or r_min)
    sim = SIM.make_sim(SIM.frame_2d(pos0, box_star), seed=args.seed)

    cell = md.nlist.Cell(buffer=0.2)
    reduced = {k: p[k] for k in ("kappa_star", "edl_amp", "vdw_amp", "a_star")}
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_used, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_used, U=U_arr, F=F_arr)
    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)

    # ── (1) short equilibration -- all particles free, no drive yet. Meaningless at
    #        kT=0, so skipped there ──
    n_eq = 0
    if args.kT > 0:
        bd_all = md.methods.Brownian(filter=hoomd.filter.All(), kT=args.kT, default_gamma=1.0)
        integ = md.Integrator(dt=dt, methods=[bd_all], forces=[tab, wca])
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        if no_barrier_mode:
            # There is no tau_bond here (no secondary minimum, hence no k_bond_star)
            # -- the number of steps needed to settle inside the "dead zone" at the
            # r_min boundary (table F=0, WCA F=0) is just taken as a multiple of dt
            # (fixed, ★an exploratory guess -- NOT a rigorous timescale).
            # ★ Keep it SHORT -- the point is only to damp the local jitter of the
            #   initial placement (exactly on the r_min boundary).
            #   Run it long (say 2e4) and a chain with no bending stiffness has
            #   already rearranged into a random coil in the meantime (measured: bead
            #   positions spread over ±400nm -- far wider than an individual bond's
            #   0.01d dead zone, due to accumulated multi-joint rotation), which
            #   breaks the drive experiment's premise of "starting from straight".
            n_eq = int(args.eq_tau_bond * 2)
        else:
            n_eq = int(round(args.eq_tau_bond
                            * (tau_bond / (dt * tau_B)).to("dimensionless").magnitude))
        sim.run(n_eq)
    print(f"equilibration of {n_eq:,} steps done (to settle ⟨U⟩/N, no drive)"
          + (" -- skipped, kT=0" if args.kT == 0 else ""))

    # ── (2) remove the centre bead from integration and write its position directly ──
    bd_filter = hoomd.filter.SetDifference(hoomd.filter.All(), hoomd.filter.Tags([mid]))
    if args.kT > 0:
        bd = md.methods.Brownian(filter=bd_filter, kT=args.kT, default_gamma=1.0)
    else:
        bd = md.methods.OverdampedViscous(filter=bd_filter, default_gamma=1.0)
    integ2 = md.Integrator(dt=dt, methods=[bd], forces=[tab, wca])
    integ2.integrate_rotational_dof = False
    sim.operations.integrator = integ2
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=_move_ghost_action(mid, amp_star, omega_star, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))

    # ── (3) drive phase -- capture frames ─────────────────────────────
    frames = [np.array(sim.state.get_snapshot().particles.position, dtype=float)[:, :2].copy()]
    done = 0
    while done < n_prod:
        chunk = min(capture_every, n_prod - done)
        sim.run(chunk)
        done += chunk
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, dtype=float)[:, :2]
        img = np.array(snap.particles.image, dtype=float)[:, :2]
        frames.append(pos + img * box_star)
    frames = np.array(frames)
    print(f"drive done -- captured {len(frames)} frames")

    t_star = np.arange(len(frames)) * capture_every * dt          # in tau_B
    t_over_period = t_star / tau_period_star
    y_all = frames[:, :, 1] - frames[:1, :, 1].mean()              # rough centring,
                                                                   # informational only

    out_dir = ROOT / "runs" / "_scratch_drive_chain_relax"
    out_dir.mkdir(exist_ok=True)
    tag = f"n{n}-w{args.omega:g}-a{args.amp:g}-kT{args.kT:g}"
    if no_barrier_mode:
        tag += f"-vdw{p['vdw_amp']:g}-rc{r_cut_star:g}-rm{r_min_used:g}"
    np.savez_compressed(out_dir / f"{tag}.npz", frames=frames, t_over_period=t_over_period,
                       omega=args.omega, amp_nm=args.amp, mid=mid, ell_star=ell_star,
                       box_star=box_star)
    make_plots(frames, t_over_period, mid, n, args, out_dir, tag)


def make_plots(frames, t_over_period, mid, n, args, out_dir, tag):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.animation as anim
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    d_nm = 1470.0
    y_nm = frames[:, :, 1] * d_nm

    # (1) y_i(t) -- per bead, fading in colour from the centre outwards
    fig1, ax1 = plt.subplots(figsize=(10, 5))
    cmap = plt.cm.viridis
    for i in range(n):
        dist = abs(i - mid)
        ax1.plot(t_over_period, y_nm[:, i], lw=1.3 if dist == 0 else 0.9,
                 color=cmap(1 - dist / max(1, mid)), alpha=1.0 if dist == 0 else 0.85,
                 label=f"bead {i}" + ("  (driven)" if i == mid else f"  (|Δi|={dist})"))
    ax1.set(xlabel="t / period", ylabel="y [nm]",
           title=f"chain-relax-2d-dlvo, driven center bead — n={n}, "
                 f"ω={args.omega:g} rad/s, A={args.amp:g} nm, kT={args.kT:g}")
    ax1.legend(fontsize=7, ncol=2)
    ax1.grid(alpha=.3)
    fig1.tight_layout()
    fig1.savefig(out_dir / f"{tag}_y_of_t.png", dpi=140)
    plt.close(fig1)

    # (2) amplitude propagation -- each bead's y amplitude (steady-state RMS*sqrt(2))
    #     against |i-mid|
    steady = t_over_period > t_over_period.max() * 0.4
    amp_meas = y_nm[steady].std(axis=0) * math.sqrt(2)
    fig2, ax2 = plt.subplots(figsize=(6, 4.5))
    dist = np.abs(np.arange(n) - mid)
    ax2.semilogy(dist, amp_meas, "o-")
    ax2.axhline(args.amp, color="gray", ls=":", lw=1, label=f"driven amplitude ({args.amp:g} nm)")
    ax2.set(xlabel="|bead index - driven bead|", ylabel="oscillation amplitude [nm] (log)",
           title="Amplitude vs distance from drive\n(no bending stiffness -> expect fast decay)")
    ax2.legend(fontsize=8); ax2.grid(alpha=.3, which="both")
    fig2.tight_layout()
    fig2.savefig(out_dir / f"{tag}_propagation.png", dpi=140)
    plt.close(fig2)
    print("propagated amplitude [nm] vs |di|:",
          dict(zip(dist.tolist(), np.round(amp_meas, 3).tolist())))

    # (3) animation -- in the body frame (origin at the driven bead rather than the
    #     ends), with labels
    xmax = float(np.abs(frames[:, :, 0] - frames[:, mid:mid+1, 0]).max()) * 1.15
    ymax = max(args.amp / d_nm * 1.4, float(np.abs(y_nm).max()) / d_nm * 1.15)
    fig3, ax3 = plt.subplots(figsize=(7, 3.6))
    ax3.set_xlim(-xmax, xmax); ax3.set_ylim(-ymax, ymax)
    ax3.set_aspect("equal"); ax3.grid(alpha=.3)
    ax3.set_xlabel("x - x_driven [d]"); ax3.set_ylabel("y [d]")
    ax3.set_title(f"n={n}, ω={args.omega:g} rad/s, A={args.amp:g} nm, kT={args.kT:g}")
    line, = ax3.plot([], [], "o-", color="tab:purple", ms=6, lw=1.5)
    driven_pt, = ax3.plot([], [], "o", color="tab:red", ms=9)
    txt = ax3.text(0.02, 0.92, "", transform=ax3.transAxes, fontsize=9)

    def update(i):
        xr = frames[i, :, 0] - frames[i, mid, 0]
        yr = frames[i, :, 1]
        line.set_data(xr, yr)
        driven_pt.set_data([xr[mid]], [yr[mid]])
        txt.set_text(f"t = {t_over_period[i]:.2f} periods")
        return line, driven_pt, txt

    ani = anim.FuncAnimation(fig3, update, frames=len(frames), interval=60, blit=False)
    ani.save(out_dir / f"{tag}_anim.gif", writer="pillow", fps=15)
    plt.close(fig3)
    print("saved:", out_dir / f"{tag}_y_of_t.png", "/",
          out_dir / f"{tag}_propagation.png",
         "/", out_dir / f"{tag}_anim.gif")


if __name__ == "__main__":
    main()
