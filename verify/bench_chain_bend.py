"""Measure the actual step rate of `chain-bend-2d-oscill` -- the basis for judging
the L4 cost.

The L3 spec demands 2.65e9 steps at the lowest omega. "It is expensive" is not
something you can decide on, so the steps per second are measured. Three
configurations are compared:

  A  bond + angle only                        -- compiled forces only. The lower
                                                 bound (the best possible speed)
  B  + an md.force.Custom trap (Python every step)
                                              -- the current design (the
                                                 trap-2d / trap-drag approach)
  C  + a bond.Harmonic(r0=0) trap to a ghost particle
                                              -- the compiled path. Only the driven
                                                 anchor is moved, by an updater

The idea behind C: bond.Harmonic(r0=0) is U = 0.5*k*r^2, exactly a harmonic trap.
Leave the ghost particle out of the integrator's filter and it does not move, giving
a fixed trap. The driven trap moves its ghost with a CustomUpdater -- and since
omega*dt = 2.4e-7, moving it every 100 steps costs a phase error of only 2.4e-5 of a
cycle. There is no reason to call into Python every step.

The parameters are read from specs/chain-bend-2d-oscill__w85__*.json (the lowest
omega, the most expensive point).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/bench_chain_bend.py
"""
from __future__ import annotations

import glob
import json
import math
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
BENCH_STEPS = 200_000          # the measured interval, after warm-up
WARMUP = 5_000
UPDATE_EVERY = 100             # how often configuration C moves the driven anchor


def load_spec() -> dict:
    """The lowest-omega spec = the most expensive point."""
    cands = sorted(glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__w85__*.json")))
    if not cands:
        raise SystemExit("specs/chain-bend-2d-oscill__w85__*.json not found")
    return json.loads(Path(cands[0]).read_text())


def build_frame(n: int, L_chain: float, *, ghosts: list[int] | None = None):
    """A straight chain along x.

    If `ghosts` is given, one extra ghost particle is attached per listed bead.
    """
    ell = L_chain / (n - 1)
    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    types = ["A"]
    typeid = [0] * n
    n_real = n
    if ghosts:
        types.append("G")
        for g in ghosts:
            pos.append(list(pos[g]))        # the ghost is placed on top of the initial position (dr=0)
            typeid.append(1)

    box_L = 4.0 * L_chain                   # much larger than the chain -> wrapping never reaches the physics
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = types
    f.configuration.box = [box_L, box_L, 0, 0, 0, 0]     # Lz=0 -> 2D (trap 9)
    f.configuration.dimensions = 2

    f.bonds.N = n - 1
    f.bonds.types = ["backbone"]
    f.bonds.typeid = [0] * (n - 1)
    bond_group = [[i, i + 1] for i in range(n - 1)]
    if ghosts:
        f.bonds.types = ["backbone", "trap"]
        for j, g in enumerate(ghosts):
            bond_group.append([g, n_real + j])
        f.bonds.N = len(bond_group)
        f.bonds.typeid = [0] * (n - 1) + [1] * len(ghosts)
    f.bonds.group = np.array(bond_group)

    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])
    return f, n_real, box_L, ell


class CustomTrap(md.force.Custom):
    """Configuration B -- Python every step. The same approach as
    trap-2d-5um / trap-drag."""

    def __init__(self, k, trapped, anchors, amp, omega, dt, drive_row):
        super().__init__(aniso=False)
        self.k = float(k)
        self.trapped = np.asarray(trapped, dtype=int)
        self.anchors = np.asarray(anchors, dtype=float)
        self.amp, self.omega, self.dt = float(amp), float(omega), float(dt)
        self.drive_row = int(drive_row)

    def set_forces(self, timestep):
        anc = self.anchors.copy()
        anc[self.drive_row, 1] += self.amp * math.sin(self.omega * timestep * self.dt)
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)
            pos = np.array(snap.particles.position, copy=True)
            arr.force[:] = 0.0
            arr.potential_energy[:] = 0.0
            for row, tg in enumerate(self.trapped):     # tag indexing is mandatory
                loc = np.flatnonzero(tags == tg)
                d = pos[loc] - anc[row]
                arr.force[loc] = -self.k * d
                arr.potential_energy[loc] = 0.5 * self.k * (d ** 2).sum(axis=1)


class MoveGhost(hoomd.custom.Action):
    """Configuration C -- move only the driven ghost, every UPDATE_EVERY steps."""

    def __init__(self, ghost_tag, y0, amp, omega, dt):
        self.ghost_tag = int(ghost_tag)
        self.y0, self.amp = float(y0), float(amp)
        self.omega, self.dt = float(omega), float(dt)

    def act(self, timestep):
        y = self.y0 + self.amp * math.sin(self.omega * timestep * self.dt)
        with self._state.cpu_local_snapshot as snap:
            tags = np.array(snap.particles.tag, copy=True)
            loc = np.flatnonzero(tags == self.ghost_tag)
            if len(loc):
                snap.particles.position[loc[0]][1] = y


def make_sim(variant: str, p: dict, nu: dict):
    n = int(p["n_beads"])
    trapped = [int(t) for t in p["trapped"]]
    k_t, k_b = float(p["k_t_star"]), float(p["k_bond_star"])
    kappa = float(p["kappa_theta_star"])
    amp, omega = float(p["amp_star"]), float(p["omega_star"])
    dt = float(nu["dt_star"])
    ghosts = trapped if variant == "C" else None

    f, n_real, box_L, ell = build_frame(n, float(p["L_chain_star"]), ghosts=ghosts)
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)

    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=k_b, r0=ell)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=kappa, t0=math.pi)
    forces = [bond, angle]

    if variant == "C":
        # bond.Harmonic(r0=0) = 0.5*k*r^2 -> identical to a harmonic trap, on the
        # compiled path.
        bond.params["trap"] = dict(k=k_t, r0=0.0)
        integrated = hoomd.filter.Type(["A"])        # the ghosts ("G") are not integrated -> fixed
    else:
        integrated = hoomd.filter.All()
        if variant == "B":
            anchors = np.array([[(t - (n - 1) / 2) * ell, 0.0, 0.0] for t in trapped])
            drive_row = trapped.index(sorted(trapped)[len(trapped) // 2])
            forces.append(CustomTrap(k_t, trapped, anchors, amp, omega, dt, drive_row))

    bd = md.methods.Brownian(filter=integrated, kT=1.0, default_gamma=1.0)
    integrator = md.Integrator(dt=dt, methods=[bd], forces=forces)
    integrator.integrate_rotational_dof = False      # BD is overdamped (trap 5)
    sim.operations.integrator = integrator

    if variant == "C":
        mid = sorted(trapped)[len(trapped) // 2]
        ghost_tag = n_real + trapped.index(mid)
        sim.operations.updaters.append(hoomd.update.CustomUpdater(
            action=MoveGhost(ghost_tag, 0.0, amp, omega, dt),
            trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    return sim


def main() -> int:
    spec = load_spec()
    p, nu = spec["params"], spec["numerics"]
    total = int(nu["n_prod"]) + int(nu["n_eq"])
    print("=" * 78)
    print("chain-bend-2d-oscill -- measured step rate (lowest omega, the most "
          "expensive point)")
    print("=" * 78)
    print(f"n_beads={p['n_beads']}  dt*={nu['dt_star']:.3e}  "
          f"steps required at this omega = {total:,}")
    print(f"measured over {BENCH_STEPS:,} steps (warm-up {WARMUP:,})\n")
    print(f"{'configuration':<44}{'steps/s':>12}{'time at this w':>14}")
    print("-" * 78)

    labels = {
        "A": "A  bond + angle only (compiled forces only)",
        "B": "B  + force.Custom trap (Python every step)",
        "C": f"C  + ghost bond trap (updater every {UPDATE_EVERY} steps)",
    }
    rates = {}
    for v in ("A", "B", "C"):
        sim = make_sim(v, p, nu)
        sim.run(WARMUP)
        t0 = time.perf_counter()
        sim.run(BENCH_STEPS)
        el = time.perf_counter() - t0
        rate = BENCH_STEPS / el
        rates[v] = rate
        days = total / rate / 86400
        span = f"{days:.1f} d" if days >= 1 else f"{total / rate / 3600:.1f} h"
        print(f"{labels[v]:<44}{rate:>12,.0f}{span:>14}")

    print("-" * 78)
    print(f"cost of the Python trap  B/A = {rates['B'] / rates['A']:.3f}  "
          f"({rates['A'] / rates['B']:.1f}x slower)")
    print(f"recovered by the ghost trap  C/B = {rates['C'] / rates['B']:.2f}x faster")

    # Cost of the whole sweep. ★ Alphabetical filename order is NOT omega order
    # (w1737 sorts before w85) -- sort by the value.
    allspecs = [json.loads(Path(s).read_text())
                for s in glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__*.json"))]
    steps = sorted(s["numerics"]["n_prod"] + s["numerics"]["n_eq"] for s in allspecs)
    tot, worst = sum(steps), steps[-1]      # lowest omega = most steps = the wall clock when run in parallel
    print(f"\nsweep of {len(allspecs)} points totals {tot:,} steps "
          f"(largest = {worst:,})")
    for v in ("B", "C"):
        ser, par = tot / rates[v] / 86400, worst / rates[v] / 3600
        print(f"  config {v}: serial {ser:6.1f} d   7 runs in parallel -> "
              f"wall clock {par:6.1f} h")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
