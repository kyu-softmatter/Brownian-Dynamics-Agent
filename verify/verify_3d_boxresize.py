#!/usr/bin/env python
"""Capability survey before running the `network` case -- 3D BD + `update.BoxResize`
(CLAUDE.md rules 4 and 7).

Why this is needed: all six existing cases in this project are **2D**, and box
compression (`hoomd.update.BoxResize`) had **never been used here** (no entry in
skill `bd-hoomd` or docs/hoomd_capabilities.md -> intake/network/observation.yaml
N4 and N6).

Five checks --
  (1) 3D free diffusion alone:  <r^2> = 6*D*t   (isolation check. N4)
  (2) does BoxResize really **affinely scale** particle coordinates?
      (measuring the documented claim)
  (3) do BoxResize + Brownian + pair.Table + WCA run together?
      (cell list, finiteness)
  (4) ★ **the threshold at which compression breaks a DLVO bond** -- this IS the
      design number for the gelation protocol.
      prediction: linear strain per trigger
      eps_crit = (h_min* - barrier_h*)/l*
      (if an affine step pushes a bond inside the barrier it falls to the primary
       minimum, irreversibly)
  (5) cost: steps/s for the gelation configuration and the wall-clock that follows
      (N3 -- measure it before deciding)

Units: d=1, kT=1, gamma=1  =>  D_t=1, tau_B = d^2/D_t = 1.
The reduced DLVO parameters are taken directly from
`cases/chain_bend_dlvo_2d.py` (so the same physics is not written twice -- those
expressions were verified in SI in that case).

Run:  $PY verify/verify_3d_boxresize.py [--quick]
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

import gsd.hoomd                                        # noqa: E402
import hoomd                                            # noqa: E402
import hoomd.md as md                                   # noqa: E402

from chain_bend_dlvo_2d import (                         # noqa: E402
    SIGMA_CORE_STAR, CUTOFF_H_STAR,
    build_table_arrays, dlvo_reduced_params, find_well, load_system,
)

R_WCA = 2 ** (1 / 6)
PASS, FAIL = [], []


def check(ok: bool, label: str, detail: str = "") -> None:
    (PASS if ok else FAIL).append(label)
    print(f"  {'✓' if ok else '✗'} {label}" + (f"   {detail}" if detail else ""))


def cpu(seed: int) -> hoomd.Simulation:
    # ★ seed < 65536 (bd-hoomd trap 12 -- it is truncated to 16 bits)
    return hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)


def frame_3d(pos, L, types=("A",), typeid=None):
    f = gsd.hoomd.Frame()
    pos = np.asarray(pos, dtype=float)
    f.particles.N = len(pos)
    f.particles.position = pos
    f.particles.typeid = [0] * len(pos) if typeid is None else list(typeid)
    f.particles.types = list(types)
    f.particles.mass = [1.0] * len(pos)
    f.configuration.box = [L, L, L, 0, 0, 0]      # ★ 3D: Lz=L (2D uses Lz=0 -- trap 9)
    f.configuration.dimensions = 3
    return f


def unwrapped(sim, L=None):
    s = sim.state.get_snapshot()
    pos = np.array(s.particles.position, copy=True)
    img = np.array(s.particles.image, copy=True)
    box = s.configuration.box
    Ls = np.array([box[0], box[1], box[2]], dtype=float) if L is None else np.full(3, L)
    return pos + img * Ls


# ═══════════════════════════════════════════════════════════════════════
# (1) 3D free diffusion -- <r^2> = 6 D t   (analytic. isolation check)
# ═══════════════════════════════════════════════════════════════════════
def check_free_diffusion_3d(n_part=800, dt=1e-4, n_steps=20_000):
    print("\n(1) 3D free diffusion (no interactions) -- <r^2> = 6*D*t")
    L = 40.0
    rng = np.random.default_rng(3)
    pos = rng.uniform(-L / 2, L / 2, size=(n_part, 3))
    sim = cpu(11)
    sim.create_state_from_snapshot(frame_3d(pos, L))
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    r0 = unwrapped(sim)
    sim.run(n_steps)
    r1 = unwrapped(sim)
    t = dt * n_steps
    msd = float(((r1 - r0) ** 2).sum(axis=1).mean())
    pred = 6.0 * 1.0 * t                         # D = kT/γ = 1
    rel = msd / pred - 1.0
    sem = float(((r1 - r0) ** 2).sum(axis=1).std(ddof=1) / math.sqrt(n_part)) / pred
    check(abs(rel) < 4 * sem + 0.02,
          "3D free diffusion agrees with 6Dt",
          f"measured {msd:.5f} / predicted {pred:.5f} = {1+rel:.4f}  "
          f"({rel*100:+.2f}%, SEM {sem*100:.2f}%)")

    # Per-component isotropy -- catches the mistake of dropping one axis when
    # moving to 3D
    per_axis = ((r1 - r0) ** 2).mean(axis=0)
    iso = per_axis / (2.0 * t)
    check(np.all(np.abs(iso - 1) < 0.12),
          "each of the three axes gives <dx^2>=2Dt (no axis dropped)",
          f"x/y/z = {iso[0]:.3f} / {iso[1]:.3f} / {iso[2]:.3f}")
    return msd, pred


# ═══════════════════════════════════════════════════════════════════════
# (2) does BoxResize affinely scale the coordinates?
# ═══════════════════════════════════════════════════════════════════════
def check_boxresize_affine():
    print("\n(2) BoxResize -- does it affinely scale particle coordinates? "
          "(measuring the documented claim)")
    L0, s = 10.0, 0.5
    pos = np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, -2.5], [0.0, 0.0, 0.0]])
    sim = cpu(12)
    sim.create_state_from_snapshot(frame_3d(pos, L0))
    # The updater with no integrator -- isolates BoxResize as the only possible
    # cause of a coordinate change
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-12, methods=[bd], forces=[])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    box0 = hoomd.Box(Lx=L0, Ly=L0, Lz=L0)
    box1 = hoomd.Box(Lx=L0 * s, Ly=L0 * s, Lz=L0 * s)
    var = hoomd.variant.box.Interpolate(
        initial_box=box0, final_box=box1,
        variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=0, t_ramp=10))
    sim.operations.updaters.append(
        hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(1), box=var))
    sim.run(11)

    got = np.array(sim.state.get_snapshot().particles.position, copy=True)
    box = sim.state.get_snapshot().configuration.box
    check(abs(box[0] - L0 * s) < 1e-9, "the box shrank to the target size",
          f"Lx {box[0]:.6f} (target {L0*s:.6f})")
    err = np.abs(got - pos * s).max()
    check(err < 1e-9, "★ coordinates are scaled exactly affinely (r -> s*r)",
          f"max error {err:.3e} -- meaning **bond lengths shrink along with it**")
    return err


# ═══════════════════════════════════════════════════════════════════════
# DLVO table + WCA core (same convention as chain-bend-2d-dlvo)
# ═══════════════════════════════════════════════════════════════════════
def dlvo_forces(nlist, P, extra_types=()):
    r_min = 1.0 + 1e-4
    r_cut = 1.0 + CUTOFF_H_STAR
    r, U, F = build_table_arrays(P, r_min, r_cut)
    tab = md.pair.Table(nlist=nlist, default_r_cut=r_cut)
    wca = md.pair.LJ(nlist=nlist, default_r_cut=SIGMA_CORE_STAR * R_WCA, mode="shift")
    for a in ("A",) + tuple(extra_types):
        for b in ("A",) + tuple(extra_types):
            tab.params[(a, b)] = dict(r_min=r_min, U=U, F=F)
            wca.params[(a, b)] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
    return [tab, wca], r_cut


# ═══════════════════════════════════════════════════════════════════════
# (3) do BoxResize + Brownian + the pair forces run together?
# ═══════════════════════════════════════════════════════════════════════
def check_boxresize_with_pair(P):
    print("\n(3) BoxResize + Brownian + pair.Table + WCA (cell list, finiteness)")
    n_side, L0 = 6, 18.0
    a = L0 / n_side
    pos = np.array([[(i + .5) * a - L0 / 2, (j + .5) * a - L0 / 2, (k + .5) * a - L0 / 2]
                    for i in range(n_side) for j in range(n_side) for k in range(n_side)])
    sim = cpu(13)
    sim.create_state_from_snapshot(frame_3d(pos, L0))
    nl = md.nlist.Cell(buffer=0.2)
    forces, r_cut = dlvo_forces(nl, P)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-7, methods=[bd], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    L1 = 12.0
    var = hoomd.variant.box.Interpolate(
        initial_box=hoomd.Box(Lx=L0, Ly=L0, Lz=L0),
        final_box=hoomd.Box(Lx=L1, Ly=L1, Lz=L1),
        variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=0, t_ramp=20_000))
    sim.operations.updaters.append(
        hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(100), box=var))
    try:
        sim.run(20_100)
        ran = True
    except Exception as e:                                     # noqa: BLE001
        ran = False
        print(f"      crash: {type(e).__name__}: {e}")
    check(ran, "no crash during compression")
    if not ran:
        return None
    pe = thermo.potential_energy
    box = sim.state.get_snapshot().configuration.box
    check(pe is not None and math.isfinite(pe), "the potential energy is finite",
          f"PE = {pe:.4f} kT (N={len(pos)})")
    check(abs(box[0] - L1) < 1e-6, "compression stops exactly at the target",
          f"Lx {box[0]:.6f} (target {L1})")
    check(r_cut < box[0] / 2, "r_cut < L/2 still holds after compression (trap 6)",
          f"r_cut {r_cut:.4f} < L/2 {box[0]/2:.4f}")
    return pe


# ═══════════════════════════════════════════════════════════════════════
# (4) ★ the threshold at which compression breaks a DLVO bond
# ═══════════════════════════════════════════════════════════════════════
def check_crush_threshold(P, W, quick=False):
    print("\n(4) ★ the threshold at which compression breaks a DLVO "
          "secondary-minimum bond")
    h_min, h_bar = W["h_min"], W["barrier_h"]
    ell = 1.0 + h_min
    eps_crit = (h_min - h_bar) / ell
    tau_bond = 1.0 / W["k_bond_star"]
    print(f"   ledger: h_min*={h_min:.6f}  barrier*={h_bar:.6f}  l*={ell:.6f}")
    print(f"         k_bond*={W['k_bond_star']:.4g} kT/d²   τ_bond*={tau_bond:.4g} τ_B")
    print(f"   predicted threshold eps_crit = (h_min*-barrier*)/l* = "
          f"{eps_crit:.6f}  ({eps_crit*100:.3f}% per trigger)")

    dt = 1e-8                       # dt/τ_bond ≈ 0.0104
    T = 2000                        # trigger interval -> relaxation time
                                    # T*dt = 21 tau_bond
    L0 = 6.0
    s_tot = 0.95                    # 5% total linear compression
    eps_list = [0.002, 0.004, 0.006, 0.008, 0.012, 0.020]
    if quick:
        eps_list = [0.004, 0.008, 0.020]

    print(f"   dt={dt:g} (dt/tau_bond={dt/tau_bond:.4f}) . trigger interval {T} steps"
          f" (= {T*dt/tau_bond:.0f} tau_bond of relaxation)")
    print("   eps/trigger   final h*       verdict")
    rows = []
    for eps in eps_list:
        t_ramp = max(int(round((1.0 - s_tot) * T / eps)), T)
        sim = cpu(14)
        sim.create_state_from_snapshot(
            frame_3d([[-ell / 2, 0, 0], [ell / 2, 0, 0]], L0))
        nl = md.nlist.Cell(buffer=0.2)
        forces, _ = dlvo_forces(nl, P)
        # kT=0, deterministic -- so noise cannot blur the threshold
        # (bd-physics: test an integrator assumption at kT=0)
        bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
        integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        var = hoomd.variant.box.Interpolate(
            initial_box=hoomd.Box(Lx=L0, Ly=L0, Lz=L0),
            final_box=hoomd.Box(Lx=L0 * s_tot, Ly=L0 * s_tot, Lz=L0 * s_tot),
            variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=0, t_ramp=t_ramp))
        sim.operations.updaters.append(
            hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(T), box=var))
        sim.run(t_ramp + T)
        p = np.array(sim.state.get_snapshot().particles.position, copy=True)
        sep = float(np.linalg.norm(p[1] - p[0]))
        h_end = sep - 1.0
        survived = h_end > h_bar * 3          # well outside the barrier = still in
                                              # the secondary minimum
        rows.append((eps, h_end, survived))
        print(f"   {eps*100:6.2f}%   {h_end:12.6f}   "
              f"{'secondary min held' if survived else '★ collapsed (primary min/contact)'}")

    ok = [e for e, _, s in rows if s]
    bad = [e for e, _, s in rows if not s]
    check(bool(ok) and bool(bad),
          "the threshold lies inside the swept range (both outcomes observed)",
          f"held <={max(ok)*100:.2f}% . collapsed >={min(bad)*100:.2f}%"
          if ok and bad else "")
    if ok and bad:
        lo, hi = max(ok), min(bad)
        inside = lo <= eps_crit <= hi or abs(eps_crit - lo) / eps_crit < 0.5
        check(inside,
              "★ the measured threshold is consistent with the analytic eps_crit",
              f"measured {lo*100:.2f}-{hi*100:.2f}% . predicted {eps_crit*100:.3f}%")
    return eps_crit, rows


# ═══════════════════════════════════════════════════════════════════════
# (5) cost -- steps/s for the gelation configuration -> wall-clock
# ═══════════════════════════════════════════════════════════════════════
def check_cost(P, W, n_part=1528, phi0=0.02, phi1=0.10, bench_steps=3000):
    print("\n(5) cost -- steps/s for the gelation configuration and the wall-clock "
          "(N3: measure before deciding)")
    tau_bond = 1.0 / W["k_bond_star"]
    dt = 1e-2 * tau_bond                          # design convention
                                                  # dt/tau_fast = 1e-2
    out = {}
    for tag, phi in (("initial (dilute)", phi0), ("final (compressed)", phi1)):
        L = (n_part * math.pi / (6.0 * phi)) ** (1 / 3)
        n_side = int(math.ceil(n_part ** (1 / 3)))
        a = L / n_side
        pos = np.array([[(i + .5) * a - L / 2, (j + .5) * a - L / 2, (k + .5) * a - L / 2]
                        for i in range(n_side) for j in range(n_side)
                        for k in range(n_side)])[:n_part]
        sim = cpu(15)
        sim.create_state_from_snapshot(frame_3d(pos, L))
        nl = md.nlist.Cell(buffer=0.2)
        forces, _ = dlvo_forces(nl, P)
        bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
        integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        sim.run(200)                                          # warm-up (cell list,
                                                              # autotuner)
        t0 = time.perf_counter()
        sim.run(bench_steps)
        rate = bench_steps / (time.perf_counter() - t0)
        out[tag] = (L, rate)
        print(f"   {tag:12s} φ={phi:.3f}  L/d={L:6.2f}  {rate:9.0f} steps/s")

    rate_min = min(r for _, r in out.values())
    print(f"\n   dt = 1e-2·τ_bond = {dt:.4g} τ_B   (τ_bond*={tau_bond:.4g})")
    for tau_target in (1.0, 5.0, 10.0):
        steps = tau_target / dt
        hours = steps / rate_min / 3600
        print(f"   gelation {tau_target:4.1f} tau_B  ->  {steps:.3g} steps  ->  "
              f"{hours:8.1f} hours/seed (at the lowest steps/s)")
    check(True, "cost measured (decide from the table above)",
          f"lowest {rate_min:.0f} steps/s")
    return out, dt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="reduce the number of sweep points in (4)")
    args = ap.parse_args()

    print("=" * 78)
    print("network capability survey -- 3D BD + BoxResize   "
          "(rule 4: never from intuition)")
    print("=" * 78)
    print(f"hoomd {hoomd.version.version} . CPU . units d=kT=gamma=1 (tau_B=1)")

    sys_ = load_system(ROOT / "intake/chain-bend-2d-dlvo/system.yaml")
    P = dlvo_reduced_params(sys_)
    W = find_well(P)
    print(f"DLVO ledger inherited: barrier {W['barrier_U']:.2f} kT @ "
          f"h*={W['barrier_h']:.5f} . "
          f"secondary min {W['U_min']:.3f} kT @ h*={W['h_min']:.5f}")

    check_free_diffusion_3d()
    check_boxresize_affine()
    check_boxresize_with_pair(P)
    check_crush_threshold(P, W, quick=args.quick)
    check_cost(P, W)

    print("\n" + "=" * 78)
    print(f"{'✓ PASS' if not FAIL else '✗ FAIL'} -- "
          f"{len(PASS)}/{len(PASS)+len(FAIL)} OK")
    for f in FAIL:
        print(f"   failed: {f}")
    print("=" * 78)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
