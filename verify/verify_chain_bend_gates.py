"""Two pre-run gates for `chain-bend-2d-oscill`. Both must pass **before** the
production sweep.

────────────────────────────────────────────────────────────────────────────
Gate A (`--gate lockin`) -- compare the measurement code against an analytic
solution (mater_plan principle 9 / rule 7)
────────────────────────────────────────────────────────────────────────────
Everything this case produces rests on **newly written lock-in extraction code**.
Remove the chain and replace it with one bead + a driving trap (k_t) + a static
ghost spring (k_s), and the answer comes out in closed form:

    y_hat = k_t a / (k_t + k_s + i*omega*gamma)   (the bead response)
    K*_sample = k_s  (real, omega-independent)    <- what the estimator must return

★ This gate **FAILED on the first attempt**, and that is why it exists. The cause was
not the physics but the drive -- moving the ghost only every U steps makes the drive
a **zero-order hold**, which attenuates the fundamental by sinc(omega*dt_h/2) and
retards the phase by omega*dt_h/2 (dt_h = U*dt).
Measured at De=10: |y_hat_c|/a = 0.98999, phase -0.2522 rad, against the ZOH
prediction of 0.99040 and -0.2404 rad.

The lesson is **do not use the nominal amplitude a in the estimator**. Measure the
ghost's position too and use the **measured phasor y_hat_c**, and the ZOH attenuation
cancels exactly between numerator and denominator -- because the bead responds to
where the ghost actually is, not to the nominal sine. The table below prints both
estimators side by side to show it (the nominal one collapses with De; the measured
one stays flat).

  (1) lock-in on the ghost position -> does it match the ZOH prediction?
      [is the drive quantitatively understood]
  (2) lock-in on the bead position -> the analytic y_hat
      [BD + trap + lock-in]
  (3) the K* estimator (using the measured y_hat_c) -> (k_s, 0)
      [the estimator as a whole]
  (4) 10x the statistics -> does the error fall as 1/sqrt(N)?
      [separating bias from noise]
Without (1) it fails silently -- if the ghost never moves, the run still completes
and only K* comes out wrong.

At the production run's dt (4.53e-10), even U=100 gives omega*dt_h = 2.19e-3 rad at
De=10, so the ZOH is harmless there. Gate A's dt is 220x larger, which amplifies the
effect -- making it **a more useful condition for exposing the estimator's weak
point** (the lesson of trap 1: verify under the weak condition).

────────────────────────────────────────────────────────────────────────────
Gate B (`--gate inertia`) -- measure tau_p/tau_fast = 0.60 rather than argue about
it (rule 6)
────────────────────────────────────────────────────────────────────────────
The fastest bending mode is not overdamped (damping ratio
zeta = gamma/2*sqrt(m*lambda_max) = 0.65 < 1). BD forces that mode to be overdamped,
so the dynamics in that band is wrong -- and no choice of dt fixes it.
"It sits 4570x away from the observed band so it will not matter" was, until now,
**an unverified inference**.

Run Brownian (no inertia) against Langevin (with inertia) at identical parameters and
compare K*(omega) in the measured band. If they agree, that inference becomes **a
measured margin**. Contamination is largest where omega is closest to tau_fast, i.e.
at the **highest omega**, so De = 10 and 4.7 are the ones examined.

────────────────────────────────────────────────────────────────────────────
    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_chain_bend_gates.py --gate lockin
    $PY scratch/verify_chain_bend_gates.py --gate inertia --method bd       --de 10
    $PY scratch/verify_chain_bend_gates.py --gate inertia --method langevin --de 10
    $PY scratch/verify_chain_bend_gates.py --gate inertia --collect
"""
from __future__ import annotations

import argparse
import glob
import json
import math
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
import sys as _sys; _sys.path.insert(0, str(ROOT))
OUT = ROOT / "verify" / "_gates"
UPDATE_EVERY = 100          # how often the driving ghost is moved. The ZOH is left in
                            # and cancelled by the estimator
SAMPLES_PER_CYCLE = 20      # the same sample density as the spec
                            # (2000 samples / 100 cycles)
N_SIGMA = 3.0               # the threshold for "indistinguishable from the analytic
                            # solution"
KAPPA_CENTER = 2.0 * 2415.33   # kappa_center* = 2 kappa_end* (the spec ledger's
                               # kappa_end_d2/kT)
TAU_P = 2.70624e-8          # the spec ledger's tau_p/tau_B
TAU_CHAIN = 2.07011e-4      # the spec ledger's tau_chain/tau_B


def load_specs() -> list[dict]:
    """Ascending in omega. ★ Alphabetical filename order is NOT omega order
    (w1737 sorts before w85)."""
    specs = [json.loads(Path(p).read_text())
             for p in glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__*.json"))]
    if not specs:
        raise SystemExit("no specs/chain-bend-2d-oscill__*.json found")
    return sorted(specs, key=lambda s: s["params"]["omega_star"])


# ════════════════════════════════════════════════════════════════════════
# Lock-in -- returns a phasor per block (the error bars come from the block scatter)
# ════════════════════════════════════════════════════════════════════════
# ★ The lock-in estimator was promoted to `bdbot/lockin.py` (used twice: gate A and
#   the production run). The body was moved across verbatim from what this file
#   verified, and numerical identity was checked.
from bdbot.lockin import agg, k_star, lockin_blocks  # noqa: E402


# ════════════════════════════════════════════════════════════════════════
# Sample collection + driving
# ════════════════════════════════════════════════════════════════════════
class Sampler(hoomd.custom.Action):
    """Record the y of both the bead and the driving ghost.

    The ghost has to be measured too, or the drive cannot be checked or corrected.
    """

    def __init__(self, bead_tag: int, ghost_tag: int, dt: float):
        self.bead_tag, self.ghost_tag, self.dt = int(bead_tag), int(ghost_tag), float(dt)
        self.t: list[float] = []
        self.y_bead: list[float] = []
        self.y_ghost: list[float] = []

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            tags = np.array(snap.particles.tag, copy=True)
            pos = np.array(snap.particles.position, copy=True)
            self.t.append(timestep * self.dt)
            self.y_bead.append(float(pos[np.flatnonzero(tags == self.bead_tag)[0], 1]))
            self.y_ghost.append(float(pos[np.flatnonzero(tags == self.ghost_tag)[0], 1]))


class MoveGhost(hoomd.custom.Action):
    """Move the driving ghost as y = y0 + a sin(omega*t).

    The point is to leave the compiled path free of Python.
    """

    def __init__(self, ghost_tag: int, y0: float, amp: float, omega: float, dt: float):
        self.ghost_tag, self.y0 = int(ghost_tag), float(y0)
        self.amp, self.omega, self.dt = float(amp), float(omega), float(dt)

    def act(self, timestep):
        y = self.y0 + self.amp * math.sin(self.omega * timestep * self.dt)
        with self._state.cpu_local_snapshot as snap:
            tags = np.array(snap.particles.tag, copy=True)
            loc = np.flatnonzero(tags == self.ghost_tag)
            if len(loc):
                snap.particles.position[loc[0], 1] = y   # write via a 2-index so it
                                                         # definitely lands


def attach(sim, dt, bead_tag, ghost_tag, amp, omega, sample_every):
    sampler = Sampler(bead_tag, ghost_tag, dt)
    sim.operations.writers.append(hoomd.write.CustomWriter(
        action=sampler, trigger=hoomd.trigger.Periodic(int(sample_every))))
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=MoveGhost(ghost_tag, 0.0, amp, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    return sampler


def run_window(sim, dt, bead_tag, ghost_tag, amp, omega, n_cycles, n_eq=0):
    """Turn the drive on, let n_eq steps pass, then collect n_cycles cycles."""
    period = 2.0 * math.pi / omega
    spc = max(SAMPLES_PER_CYCLE, int(round(period / dt)))
    sample_every = max(1, spc // SAMPLES_PER_CYCLE)
    smp = attach(sim, dt, bead_tag, ghost_tag, amp, omega, sample_every)
    if n_eq:
        sim.run(int(n_eq))
        smp.t.clear(); smp.y_bead.clear(); smp.y_ghost.clear()
    sim.run(int(round(n_cycles * period / dt)))
    return (np.array(smp.t), np.array(smp.y_bead), np.array(smp.y_ghost))


# ════════════════════════════════════════════════════════════════════════
# Gate A -- one bead + a driving trap + a static ghost spring
# ════════════════════════════════════════════════════════════════════════
def build_single(k_t: float, k_s: float, dt: float):
    """tag 0 = the bead, 1 = the driving ghost (k_t), 2 = the static ghost (k_s).

    Ghosts are not integrated.
    """
    f = gsd.hoomd.Frame()
    f.particles.N = 3
    f.particles.position = np.zeros((3, 3))
    f.particles.typeid = [0, 1, 1]
    f.particles.types = ["A", "G"]
    f.configuration.box = [40.0, 40.0, 0, 0, 0, 0]          # Lz=0 -> 2D (trap 9)
    f.configuration.dimensions = 2
    f.bonds.N = 2
    f.bonds.types = ["trap", "spring"]
    f.bonds.typeid = [0, 1]
    f.bonds.group = np.array([[0, 1], [0, 2]])
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=3)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["trap"] = dict(k=k_t, r0=0.0)     # U = 0.5*k*r^2 = a harmonic trap
    bond.params["spring"] = dict(k=k_s, r0=0.0)
    bd = md.methods.Brownian(filter=hoomd.filter.Type(["A"]), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[bond])
    integ.integrate_rotational_dof = False       # BD is overdamped (trap 5)
    return sim, integ


def measure_single(k_t, k_s, amp, omega, dt, n_cycles, n_blocks=10):
    sim, integ = build_single(k_t, k_s, dt)
    sim.operations.integrator = integ
    t, yb, yg = run_window(sim, dt, 0, 1, amp, omega, n_cycles,
                           n_eq=int(round(0.1 * n_cycles * 2 * math.pi / omega / dt)))
    yb_b = lockin_blocks(t, yb, omega, n_blocks=n_blocks)
    yg_b = lockin_blocks(t, yg, omega, n_blocks=n_blocks)
    h3_b = lockin_blocks(t, yb, omega, harmonic=3, n_blocks=n_blocks)
    K_b = np.array([k_star(y, g, k_t, omega) for y, g in zip(yb_b, yg_b)])
    Kn_b = np.array([k_star(y, complex(amp, 0.0), k_t, omega) for y in yb_b])
    return dict(y=agg(yb_b), g=agg(yg_b), K=agg(K_b), Kn=agg(Kn_b),
                h3=abs(agg(h3_b)[0]), n=len(t))


def gate_lockin() -> int:
    specs = load_specs()
    k_t = float(specs[0]["params"]["k_t_star"])
    amp = float(specs[0]["params"]["amp_star"])
    k_s, gamma, n_cycles = KAPPA_CENTER, 1.0, 100
    dt = 1e-3 * gamma / (k_t + k_s)              # BD is O(dt) -- keep it generous
                                                 # (trap 2)

    print("=" * 104)
    print("Gate A -- lock-in + the K* estimator against the analytic solution "
          "(one bead + driving trap + static spring)")
    print("=" * 104)
    print(f"k_t = {k_t:.2f}   k_s = κ_center* = {k_s:.2f}   a = {amp:.5f}   "
          f"gamma = {gamma}   dt = {dt:.3e}   {n_cycles} cycles")
    print(f"expected: K* = ({k_s:.2f}, 0), real and omega-independent.   "
          f"verdict = within {N_SIGMA:.0f} sigma of the analytic solution\n")
    print(f"{'De':>6} |{'|y_c|/a':>9}{'ZOH pred':>9} |{'y err%':>9} |"
          f"{'K1(nom a)':>11}{'err%':>8} |{'K1(meas y_c)':>13}{'±sig':>8}{'err%':>8}"
          f"{'K2':>9}{'±sig':>7} |{'3rd/1st':>9} {'':>3}")
    print("-" * 104)

    fails, worst = [], 0.0
    for sp in specs:
        omega, de = float(sp["params"]["omega_star"]), float(sp["params"]["De"])
        m = measure_single(k_t, k_s, amp, omega, dt, n_cycles)
        (y_hat, _), (g_hat, _), (K, Ks), (Kn, _) = m["y"], m["g"], m["K"], m["Kn"]

        x = omega * dt * UPDATE_EVERY / 2.0
        zoh = math.sin(x) / x if x else 1.0
        y_ex = k_t * g_hat / complex(k_t + k_s, omega * gamma)   # the analytic
                                                                # solution referenced
                                                                # to the MEASURED drive
        err_y = 100.0 * abs(y_hat - y_ex) / abs(y_ex)
        err_K = 100.0 * abs(K.real - k_s) / k_s
        err_Kn = 100.0 * abs(Kn.real - k_s) / k_s
        ok = (abs(K.real - k_s) <= N_SIGMA * Ks) and (abs(K.imag) <= N_SIGMA * Ks)
        worst = max(worst, err_K)
        if not ok:
            fails.append(round(de, 2))
        print(f"{de:>6.2f} │{abs(g_hat) / amp:>9.5f}{zoh:>9.5f} │{err_y:>9.3f} │"
              f"{Kn.real:>11.1f}{err_Kn:>8.2f} │{K.real:>13.1f}{Ks:>8.1f}{err_K:>8.2f}"
              f"{K.imag:>9.1f}{Ks:>7.1f} │{m['h3'] / abs(y_hat):>9.1e} "
              f"{'✓' if ok else '✗':>3}")

    print("-" * 104)
    print(f"worst error of the measured-y_c estimator: {worst:.2f}%   "
          f"{'✓ within 3 sigma of the analytic solution at every omega' if not fails else f'✗ FAIL at De={fails}'}")

    # (4) bias vs noise -- raise the statistics 10x; if the error falls as 1/sqrt(N)
    #     it is not a bias
    print("\n(4) bias check -- 10x the statistics (De ~ 1). With no bias, both the "
          "error and sigma fall by sqrt(10)")
    sp1 = min(specs, key=lambda s: abs(s["params"]["De"] - 1.0))
    om1 = float(sp1["params"]["omega_star"])
    print(f"{'cycles':>8}{'K1':>12}{'±sig':>9}{'|K1-k_s|':>11}{'in sig':>11}{'K2':>10}")
    prev = None
    for nc, nb in ((100, 10), (1000, 10)):
        m = measure_single(k_t, k_s, amp, om1, dt, nc, n_blocks=nb)
        K, Ks = m["K"]
        d = abs(K.real - k_s)
        print(f"{nc:>8}{K.real:>12.1f}{Ks:>9.1f}{d:>11.1f}{d / Ks:>11.2f}{K.imag:>10.1f}")
        prev = (d, Ks) if prev is None else prev
    print(f"  (expected: from 100 to 1000 cycles, sigma falls by about "
          f"sqrt(10) = 3.16)")
    print("=" * 104)
    return 0 if not fails else 1


# ════════════════════════════════════════════════════════════════════════
# Gate B -- the full chain, Brownian vs Langevin
# ════════════════════════════════════════════════════════════════════════
def build_chain(sp: dict, mass: float):
    p = sp["params"]
    n = int(p["n_beads"])
    trapped = sorted(int(t) for t in p["trapped"])
    ell = float(p["L_chain_star"]) / (n - 1)
    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:                            # place each ghost on top of its bead
        pos.append(list(pos[g]))
        typeid.append(1)

    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.particles.mass = [mass] * n + [1.0] * len(trapped)
    f.configuration.box = [4.0 * float(p["L_chain_star"])] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2

    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(trapped)]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(trapped)
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=5)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=float(p["k_bond_star"]), r0=ell)
    bond.params["trap"] = dict(k=float(p["k_t_star"]), r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=float(p["kappa_theta_star"]), t0=math.pi)
    mid = trapped[len(trapped) // 2]
    return sim, [bond, angle], mid, n + trapped.index(mid)


def gate_inertia_one(de_target: float, method: str, n_cycles: int, n_eq_tau: float,
                     eq_steps: int = 0) -> dict:
    specs = load_specs()
    sp = min(specs, key=lambda s: abs(s["params"]["De"] - de_target))
    p, nu = sp["params"], sp["numerics"]
    omega, amp = float(p["omega_star"]), float(p["amp_star"])
    k_t, dt, gamma = float(p["k_t_star"]), float(nu["dt_star"]), 1.0
    # Only the inertial side gets a mass. With kT=0 both methods are
    # **deterministic** -- the statistical error is 0, so the inertial term can be
    # isolated on its own (the thermal comparison had no statistical power).
    inertial = method in ("langevin", "lang0")
    mass = TAU_P * gamma if inertial else 0.0

    sim, forces, bead_tag, ghost_tag = build_chain(sp, mass if mass else 1.0)
    filt = hoomd.filter.Type(["A"])
    integ_m = {
        "bd": lambda: md.methods.Brownian(filter=filt, kT=1.0, default_gamma=gamma),
        "langevin": lambda: md.methods.Langevin(filter=filt, kT=1.0, default_gamma=gamma),
        # Overdamped with no thermal noise -- the integrator skill bd-hoomd
        # designates for deterministic verification
        "ov": lambda: md.methods.OverdampedViscous(filter=filt, default_gamma=gamma),
        "lang0": lambda: md.methods.Langevin(filter=filt, kT=0.0, default_gamma=gamma),
    }[method]()
    integ = md.Integrator(dt=dt, methods=[integ_m], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    if method == "langevin":
        sim.state.thermalize_particle_momenta(filter=filt, kT=1.0)

    # ★ TAU_CHAIN (= gamma/kappa_center) is **NOT** this system's longest relaxation
    # time -- tau_max = gamma/lambda_min, from the smallest eigenvalue of (A+T), is
    # 9.18x longer. Equilibration must use that one.
    n_eq = int(eq_steps) if eq_steps else int(round(n_eq_tau * TAU_CHAIN / dt))
    t, yb, yg = run_window(sim, dt, bead_tag, ghost_tag, amp, omega, n_cycles, n_eq=n_eq)

    # ★ A block must contain **at least one full cycle**. Locking in on half-cycle
    # blocks makes the per-block values meaningless and turns the SEM into a number
    # unrelated to the physics (the whole-window estimate is unaffected).
    nb = max(1, min(6, n_cycles))
    yb_b = lockin_blocks(t, yb, omega, n_blocks=nb)
    yg_b = lockin_blocks(t, yg, omega, n_blocks=nb)
    h3_b = lockin_blocks(t, yb, omega, harmonic=3, n_blocks=nb)
    K_b = np.array([k_star(y, g, k_t, omega, gamma, mass) for y, g in zip(yb_b, yg_b)])
    (y_hat, y_sem), (g_hat, _), (K, K_sem) = agg(yb_b), agg(yg_b), agg(K_b)

    return dict(method=method, de=float(p["De"]), omega=omega, mass=mass,
                n_eq=n_eq, n_cycles=n_cycles, n_samples=len(t),
                drive_abs=abs(g_hat) / amp, y_abs=abs(y_hat), y_sem=y_sem,
                K_re=K.real, K_im=K.imag, K_sem=K_sem,
                h3_rel=abs(agg(h3_b)[0]) / abs(y_hat))


def gate_inertia_collect() -> int:
    rows = [json.loads(f.read_text()) for f in sorted(OUT.glob("inertia_*.json"))]
    if not rows:
        raise SystemExit(f"no {OUT}/inertia_*.json found -- run each configuration "
                         f"first")

    print("=" * 96)
    print("Gate B -- measure the contamination from tau_p/tau_fast = 0.60 "
          "(Brownian vs Langevin, identical parameters)")
    print("=" * 96)
    print(f"kappa_center* = {KAPPA_CENTER:.1f}   zeta_fast = 0.65 "
          f"(< 1 -> the fastest mode really is underdamped)")
    print("Langevin has inertia and BD does not. If K* agrees in the measured band, "
          "BD is a valid instrument for it.\n")
    print(f"{'De':>6} {'method':<10}{'|y_c|/a':>10}{'K1':>11}{'±sig':>9}{'K2':>11}{'±sig':>9}"
          f"{'3rd/1st':>10}{'samp':>7}")
    print("-" * 96)
    for r in sorted(rows, key=lambda x: (-x["de"], x["method"])):
        print(f"{r['de']:>6.2f} {r['method']:<10}{r['drive_abs']:>10.6f}"
              f"{r['K_re']:>11.1f}{r['K_sem']:>9.1f}{r['K_im']:>11.1f}{r['K_sem']:>9.1f}"
              f"{r['h3_rel']:>10.1e}{r['n_samples']:>7}")
    print("-" * 96)

    verdict = 0
    for de in sorted({r["de"] for r in rows}, reverse=True):
        pair = {r["method"]: r for r in rows if r["de"] == de}
        if len(pair) < 2:
            print(f"De={de:.2f}: unpaired, cannot compare ({list(pair)})")
            verdict = 1
            continue
        b, l = pair["bd"], pair["langevin"]
        for lbl, kb, kl in (("K′", b["K_re"], l["K_re"]), ("K″", b["K_im"], l["K_im"])):
            d, sig = abs(kb - kl), math.hypot(b["K_sem"], l["K_sem"])
            nsig = d / sig if sig > 0 else float("inf")
            rel = 100.0 * d / max(abs(kb), abs(kl), 1e-30)
            ok = nsig <= N_SIGMA
            verdict |= 0 if ok else 1
            print(f"De={de:>5.2f} {lbl}: |BD−Langevin| = {d:9.2f} = {nsig:5.2f}σ "
                  f"({rel:6.2f}%)  "
                  f"{'✓ indistinguishable' if ok else '✗ a significant difference'}")
    print("-" * 96)
    print("✓ PASS -- inertia does not change K* in the measured band. The sweep may "
          "use BD."
          if verdict == 0 else
          "✗ FAIL -- inertia reaches into the measured band. BD cannot measure it.")
    print("=" * 96)
    return verdict


def gate_det_collect() -> int:
    """Deterministic (kT=0) comparison -- the statistical error is 0, so the inertial
    term is isolated exactly.

    Side product: a noise-free K*(omega) curve. Whether the low-frequency limit goes
    to kappa_center is the implementation_check; the **shape** of the curve is a
    hypothesis (a single Maxwell mode vs a mode spectrum).
    """
    rows = [json.loads(f.read_text()) for f in sorted(OUT.glob("det_*.json"))]
    if not rows:
        raise SystemExit(f"no {OUT}/det_*.json found")

    print("=" * 100)
    print("Gate B' -- deterministic comparison (OverdampedViscous vs Langevin kT=0). "
          "Statistical error 0")
    print("=" * 100)
    print(f"kappa_center* = {KAPPA_CENTER:.1f}   m* = {TAU_P:.3e}   "
          f"zeta_fast = 0.65 (the fastest mode is underdamped)")
    print("The thermal comparison had no statistical power because |y_hat|/l_k < 1. "
          "With kT=0 that problem disappears.\n")
    print(f"{'De':>7} |{'K1(overdamp)':>13}{'K1(inertia)':>13}{'diff%':>8} |"
          f"{'K2(overdamp)':>13}{'K2(inertia)':>13}{'diff%':>8} |{'K1/kap_c':>8}")
    print("-" * 100)

    worst, pairs = 0.0, 0
    for de in sorted({r["de"] for r in rows}):
        pr = {r["method"]: r for r in rows if r["de"] == de}
        if len(pr) < 2:
            print(f"{de:>7.3f} | unpaired ({list(pr)})")
            continue
        o, l = pr["ov"], pr["lang0"]
        pairs += 1
        d_re = 100.0 * abs(o["K_re"] - l["K_re"]) / max(abs(o["K_re"]), 1e-30)
        d_im = 100.0 * abs(o["K_im"] - l["K_im"]) / max(abs(o["K_im"]), 1e-30)
        worst = max(worst, d_re, d_im)
        print(f"{de:>7.3f} │{o['K_re']:>13.1f}{l['K_re']:>13.1f}{d_re:>8.3f} │"
              f"{o['K_im']:>13.1f}{l['K_im']:>13.1f}{d_im:>8.3f} │"
              f"{o['K_re'] / KAPPA_CENTER:>8.3f}")
    print("-" * 100)
    ok = worst < 1.0
    print(f"{pairs} pairs . largest difference with/without inertia = {worst:.3f}%")
    print("✓ PASS -- the inertial term changes K*(omega) in the measured band by "
          "less than 1%. The sweep may use BD."
          if ok else
          f"✗ FAIL -- inertia changes it by up to {worst:.1f}%. BD cannot measure "
          f"this band.")
    print("  (this test looks at the inertial term alone. It does NOT cover "
          "thermal ringing combined with")
    print("   nonlinear coupling -- max|theta| = 2.8e-3 rad suggests that is small, "
          "but needs a separate check.)")
    print("=" * 100)
    return 0 if ok else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--gate", choices=["lockin", "inertia", "det"], required=True)
    ap.add_argument("--method", choices=["bd", "langevin", "ov", "lang0"])
    ap.add_argument("--de", type=float, default=10.0)
    ap.add_argument("--cycles", type=int, default=60)
    ap.add_argument("--eq-tau", type=float, default=5.0)
    ap.add_argument("--eq-steps", type=int, default=0,
                    help="give the equilibration steps as an absolute number "
                         "(when specifying them in units of tau_max)")
    ap.add_argument("--collect", action="store_true")
    a = ap.parse_args()

    if a.gate == "lockin":
        return gate_lockin()
    if a.collect:
        return gate_det_collect() if a.gate == "det" else gate_inertia_collect()
    if not a.method:
        raise SystemExit("--method bd|langevin|ov|lang0 is required")
    OUT.mkdir(parents=True, exist_ok=True)
    res = gate_inertia_one(a.de, a.method, a.cycles, a.eq_tau, a.eq_steps)
    pre = ("deq" if a.eq_steps else "det") if a.gate == "det" else "inertia"
    (OUT / f"{pre}_de{res['de']:.3f}_{a.method}.json").write_text(json.dumps(res, indent=2))
    print(f"[de{res['de']:.2f} {a.method}] K* = ({res['K_re']:.1f}, {res['K_im']:.1f}) "
          f"± {res['K_sem']:.1f}   ghost {res['drive_abs']:.6f}   "
          f"samples {res['n_samples']}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
