"""Harmonic bending (linearized JKR) vs rolling resistance -- dynamic response
compared at identical geometry and identical driving.

Statically this is already settled (`verify/verify_rolling_contact.py`, 24/24):
  . tangential spring only -> bending stiffness **exactly 0** (at delta x100 the
    residual falls by 1/10^4, i.e. it was rounding)
  . rolling resistance only -> kappa_theta,eff = 0.5*k_r*R^2, agreeing with harmonic
    bending to **within 1e-5** (in the orientation-relaxed limit)
  . **freeze** the orientations and it stiffens by 6x, 22x, 66x (it penalises
    absolute rotation, not curvature)

=> One question remains: **at this system's drive frequency, do the orientational
   degrees of freedom have time to relax?**
   If they do the two models are the same physics; if not, the rolling model is far
   stiffer.
   Predicted crossover: `omega_c = 1/tau_rot`,  `tau_rot = gamma_r/(k_r R^2)`.

★ Measured at kT=0, deterministically -- this project's lesson ("test an integrator
   assumption at kT=0"). With zero noise, 4 cycles give a clean lock-in and the two
   models' transients cancel as common mode.
★ Rule 7 (isolation): instead of the DLVO table, a **radial bond whose stiffness
   equals the well curvature** is used, so only the difference between the two
   bending models remains (DLVO's nonlinearity is irrelevant to this question).

    $PY scratch/compare_bending_models.py            # sweep + figures
"""
from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "verify"))

from bdbot import lockin as LI                                          # noqa: E402
from rolling_contact import (k_roll_from_kappa_theta,                   # noqa: E402
                             make_rolling_force)

# ── Actual values from the chain-bend-2d-dlvo spec
#    (n=9, omega=3000, a=632nm, --jkr) ──────────
H_MIN_STAR = 0.00759259035993831
ELL = 1.0 + H_MIN_STAR
K_BOND = 1042362.8817700658          # DLVO secondary-minimum curvature [kT/d^2]
                                     # -- the radial stiffness
KAPPA_THETA = 1391229.7767209478     # [kT] — κ₀=64 mN/m → EI/ℓ
# ★ The trap is made as stiff as the radial bond (NOT the production value of 5217).
#   Why: the production trap's relaxation mode is tau = gamma/k_t = 1.9e-4, which is
#   **30,000x the period of the highest omega**, so settling in units of periods never
#   drains the transient and the lock-in is contaminated (the first sweep gave K''<0
#   and a ratio oscillating 0.44<->1.39 -- not physics, just unsettled). The only
#   question here is the difference between two bending models, so the trap need only
#   be a fast boundary condition (rule 7 isolation).
K_T = K_BOND                         # trap stiffness [kT/d^2]
AMP = 0.05                           # drive amplitude [d] -- at kT=0, irrelevant if
                                     # the response is linear
R_C = 0.5                            # contact radius = d/2
GAMMA_R = 4 * R_C ** 2 / 3           # γ_r/γ_t = 8πηa³/6πηa = 4a²/3  (a=R_C, γ_t=1)
K_ROLL = k_roll_from_kappa_theta(KAPPA_THETA, R_C)
TAU_ROT = GAMMA_R / (K_ROLL * R_C ** 2)
N_BEADS = 5
UPDATE_EVERY = 1                     # move the ghost every step (minimises ZOH
                                     # attenuation -- trap 17)
OUT = ROOT / "runs" / "_bending_model_compare"


def bending_matrix(n, kappa_theta, ell):
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def make_harmonic_bending(A, n_real):
    """Harmonic bending F_y = -A y (identical to chain-bend-2d-dlvo's --jkr)."""
    class Bending(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=False)
            self.A = np.ascontiguousarray(A, dtype=float)
            self.n = int(n_real)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)
                pos = np.array(snap.particles.position, copy=True)
                m = tags < self.n
                y = np.zeros(self.n)
                y[tags[m]] = pos[m, 1]
                fy = -(self.A @ y)
                arr.force[:] = 0.0
                arr.potential_energy[:] = 0.0
                arr.force[m, 1] = fy[tags[m]]
                arr.potential_energy[m] = -0.5 * y[tags[m]] * fy[tags[m]]
    return Bending()


class ClampAndDrive(hoomd.custom.Action):
    """★ Strain control -- pin the end beads' (x,y) and force the centre bead's y
    directly.

    Why not a trap: trap driving suffers compliance, so above omega*tau_trap >> 1 the
    tracking ratio falls below 1% and `K* = k_t*y_c_hat/y_hat - k_t` becomes a ratio
    of large numbers with a collapsing condition number (measured in the first sweep:
    tracking 0.01 down to 0.00001, with K' flipping negative). Position forcing has
    100% tracking by definition, so **the conditioning is the same at every omega**.
    ★ The bead is **left in** the integrator's filter, so its rotational degree of
    freedom keeps being integrated (in G4 the end particles' orientational relaxation
    was part of the stiffness). Only the position is overwritten each step.
    """

    def __init__(self, clamp_tags, clamp_xy, drive_tag, amp, omega, dt):
        self.tags = np.asarray(clamp_tags, int)
        self.xy = np.asarray(clamp_xy, float)
        self.drive_tag, self.amp, self.omega, self.dt = int(drive_tag), amp, omega, dt

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            for k, tg in enumerate(self.tags):
                idx = snap.particles.rtag[tg]
                snap.particles.position[idx, 0] = self.xy[k, 0]
                snap.particles.position[idx, 1] = (
                    self.amp * math.sin(self.omega * timestep * self.dt)
                    if tg == self.drive_tag else self.xy[k, 1])


def build(model: str, omega: float, dt: float, n=N_BEADS):
    mid = n // 2
    clamped = [0, mid, n - 1]
    pos0 = np.zeros((n, 3))
    pos0[:, 0] = (np.arange(n) - (n - 1) / 2) * ELL
    quat0 = np.tile(np.array([1.0, 0, 0, 0]), (n, 1))

    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.types = ["A"]
    f.particles.typeid = [0] * n
    f.particles.position = pos0
    f.particles.orientation = [(1, 0, 0, 0)] * n
    f.particles.moment_inertia = [(1.0, 1.0, 1.0)] * n
    f.particles.diameter = [1.0] * n            # bd-hoomd trap 19 (unused here, but
                                                # kept as habit)
    f.bonds.N = n - 1
    f.bonds.types = ["radial"]
    f.bonds.typeid = [0] * (n - 1)
    f.bonds.group = [[i, i + 1] for i in range(n - 1)]
    L = 8.0 * n * ELL
    f.configuration.box = [L, L, 0, 0, 0, 0]
    f.configuration.dimensions = 2

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["radial"] = dict(k=K_BOND, r0=ELL)
    forces = [bond]
    if model == "harmonic":
        forces.append(make_harmonic_bending(bending_matrix(n, KAPPA_THETA, ELL), n))
    elif model == "rolling":
        forces.append(make_rolling_force([[i, i + 1] for i in range(n - 1)], pos0, quat0,
                                         R_C, K_ROLL, 0.0, n))
    elif model != "none":
        raise ValueError(model)

    # kT=0 deterministic overdamped (OverdampedViscous -- no noise)
    meth = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0,
                                        default_gamma_r=(GAMMA_R,) * 3)
    integ = md.Integrator(dt=dt, methods=[meth], forces=forces)
    integ.integrate_rotational_dof = (model == "rolling")
    sim.operations.integrator = integ
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=ClampAndDrive(clamped, pos0[clamped, :2], mid, AMP, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    return sim, mid, forces


def timescales(model: str, n=N_BEADS):
    """This system's timescales.

    The fastest sets dt; the slowest sets the **settling time**.
    """
    A = bending_matrix(n, KAPPA_THETA, ELL)
    mid = n // 2
    free = [i for i in range(n) if i not in (0, mid, n - 1)]
    lam = np.linalg.eigvalsh(A)
    lam_ff = np.linalg.eigvalsh(A[np.ix_(free, free)])       # relaxation of the free
                                                             # beads
    fast = [1.0 / K_BOND, 1.0 / lam.max()]
    slow = [1.0 / max(lam_ff.min(), 1e-30)]
    if model == "rolling":
        fast.append(TAU_ROT)
        slow.append(TAU_ROT)
    return min(fast), max(slow)


def run_one(model: str, omega: float, *, n_cycles=4, samples_per_cycle=48,
            dt_div=200, settle_taus=15):
    """★ Measure K* under strain control -- force the centre bead's y and measure
    **the force on that bead**.

        K* = -F_hat / y_hat    (F is the sample force only: radial bond + bending.
                                Solvent drag is excluded)

    With no trap compliance, tracking is 100% by definition and the condition number
    is independent of omega.
    """
    tau_fast, tau_slow = timescales(model)
    period = 2 * math.pi / omega
    dt = min(tau_fast / dt_div, period / 2000)
    settle_steps = int(max(settle_taus * tau_slow, 2 * period) / dt)
    n_meas = int(n_cycles * period / dt)
    sim, mid, forces = build(model, omega, dt)
    t0 = time.time()
    sim.run(settle_steps)
    ts, ym, fm = [], [], []
    n_chunk = max(1, n_meas // (n_cycles * samples_per_cycle))
    done = 0
    while done < n_meas:
        sim.run(min(n_chunk, n_meas - done))
        done += min(n_chunk, n_meas - done)
        p = np.array(sim.state.get_snapshot().particles.position)
        ts.append((settle_steps + done) * dt)
        ym.append(float(p[mid, 1]))
        fm.append(float(sum(np.array(f.forces)[mid, 1] for f in forces)))
    wall = time.time() - t0
    ts, ym, fm = np.array(ts), np.array(ym), np.array(fm)
    nb = min(8, max(2, len(ts) // 12))
    by = LI.lockin_blocks(ts, ym, omega, n_blocks=nb)
    bf = LI.lockin_blocks(ts, fm, omega, n_blocks=nb)
    Kb = -bf / by
    K, Ksem = LI.agg(Kb)
    yh, _ = LI.agg(by)
    return dict(K_re=K.real, K_im=K.imag, K_sem=Ksem, follow=abs(yh) / AMP,
                y_hat=abs(yh), steps=settle_steps + n_meas, dt=dt, wall=wall,
                rate=(settle_steps + n_meas) / max(wall, 1e-9))


# ══════════════════════════════════════════════════════════════════════════
if __name__ == "__main__":
    OUT.mkdir(parents=True, exist_ok=True)
    print("=" * 88)
    print("harmonic bending (linearized JKR) vs rolling resistance -- "
          "kT=0 deterministic omega sweep")
    print("=" * 88)
    print(f"  n = {N_BEADS},  κ_θ = {KAPPA_THETA:.6g} kT,  k_r = 2κ_θ/R² = {K_ROLL:.6g}")
    print(f"  γ_r = 4a²/3 = {GAMMA_R:.6f}   →   τ_rot = γ_r/(k_rR²) = {TAU_ROT:.4e}")
    print(f"  ★ predicted crossover  omega_c = 1/tau_rot = {1/TAU_ROT:.4e}  (reduced)")
    print(f"     the production run's omega* = 18453  ->  "
          f"omega*tau_rot = {18453*TAU_ROT:.3e}  "
          f"(if <<1 the two models must agree)")
    print()

    # Absolute reference for the static limit (MD has to reproduce this)
    A = bending_matrix(N_BEADS, KAPPA_THETA, ELL)
    mid = N_BEADS // 2
    fx, fr = [0, mid, N_BEADS - 1], [i for i in range(N_BEADS) if i not in (0, mid, N_BEADS - 1)]
    yfx = np.array([0.0, 1.0, 0.0])
    yv = np.zeros(N_BEADS); yv[fx] = yfx
    yv[fr] = np.linalg.solve(A[np.ix_(fr, fr)], -A[np.ix_(fr, fx)] @ yfx)
    K_STATIC = float(yv @ A @ yv)
    print(f"  ★ static-limit reference (exact linear response)  "
          f"K'(omega->0) = {K_STATIC:.6g} kT/d^2")
    print()

    # omega=1e4 and 1e5 were dropped -- their periods are long enough that, against a
    # dt set by the fast mode, each point costs 70 minutes, and the quasi-static limit
    # is already held exactly by the analytic K_STATIC above.
    omegas = [3e5, 1e6, 3e6, 1e7, 3e7, 1e8]
    rows = []
    print(f"  {'omega*':>10} {'om.tau_rot':>10} | {'harmonic K1':>13} "
          f"{'rolling K1':>13} {'roll/harm':>10} "
          f"| {'harm/static':>9} | {'steps/s':>8}")
    print("  " + "-" * 96)
    for om in omegas:
        h = run_one("harmonic", om)
        r = run_one("rolling", om)
        rows.append(dict(omega=om, harmonic=h, rolling=r, K_static=K_STATIC))
        print(f"  {om:10.3g} {om*TAU_ROT:10.3e} | {h['K_re']:13.6g} {r['K_re']:13.6g} "
              f"{r['K_re']/h['K_re']:10.5f} | {h['K_re']/K_STATIC:9.5f} "
              f"| {r['rate']:8.0f}")
        (OUT / "sweep.json").write_text(json.dumps(rows, indent=1, default=float))

    print()
    print("  ★ dt convergence check (omega*=1e7, dt halved):")
    for m in ("harmonic", "rolling"):
        a = run_one(m, 1e7, dt_div=200)
        b = run_one(m, 1e7, dt_div=400)
        print(f"    {m:9s}  K′ {a['K_re']:.6g} → {b['K_re']:.6g}   "
              f"change {100*(b['K_re']/a['K_re']-1):+.3f}%")
    print()
    print(f"  → {OUT/'sweep.json'}")
