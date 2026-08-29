"""Narrow chain-bend's **unexplained 28% discrepancy** by isolation (rule 7).

Background: K*(omega) from the small-angle linear response
(i*omega*gamma*I + A + T) y_hat = k_t a e_mid disagrees by up to 28% at high
frequency with a well-equilibrated deterministic HOOMD measurement. Already
established --
  . bending matrix vs HOOMD angle.Harmonic: agree to 0.55% (static, no traps,
    rigidly clamped)
  . lambda_max: matches the spec
  . ghost traps: within 3 sigma of the analytic solution at all 7 omega, in gate A
  . transient: converged by 20 tau_max (sigma/K = 0.02%)
The remaining candidates are peeled off one rung at a time.

  (1) dt (explicit-Euler discretisation) -- settled **without simulating**, using
      the exact z-domain solution. The steady response of explicit Euler
      y_{n+1} = (I - (A+T)dt/gamma) y_n + (k_t dt/gamma) u_n is the **closed form**
      (e^{i*omega*dt} I - M) y_hat = B a_hat. Comparing with the continuous solution
      gives the dt contribution exactly.
  (2) static + traps -- the static stiffness with the traps on, measured directly in
      HOOMD (kT=0, large dt). What (1) verified was 'rigidly clamped, no traps'.
      Does it change once the traps are in?
  (3) dynamic mode shape -- measure the complex phasor y_hat_i of all 25 beads and
      compare against the model, to see **where spatially** the disagreement arises
      (only the driven bead? the ends? an overall scale?)

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/diagnose_chain_bend_28pct.py --stage 12      # (1) and (2) (fast)
    $PY scratch/diagnose_chain_bend_28pct.py --stage 3 --de 91.8
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import sys
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
GATES = ROOT / "verify" / "_gates"
sys.path.insert(0, str(ROOT))
from cases.chain_bend_2d import bending_matrix, sweep_specs, trapped_indices   # noqa: E402


def setup():
    sp = sweep_specs()                      # works without a spec file (a failed hard check is the normal
                      # state here)
    p = sp[0]["params"]
    n = int(p["n_beads"])
    kth = float(p["kappa_theta_star"]); kt = float(p["k_t_star"])
    kb = float(p["k_bond_star"])
    ell = float(p["L_chain_star"]) / (n - 1)
    amp = float(p["amp_star"]); dt = float(sp[0]["numerics"]["dt_star"])
    idx = trapped_indices(n); mid = idx[len(idx) // 2]
    A = bending_matrix(n, kth, ell)
    T = np.zeros((n, n))
    for e in idx:
        T[e, e] += kt
    return dict(specs=sp, n=n, kth=kth, kt=kt, kb=kb, ell=ell, amp=amp, dt=dt,
                idx=idx, mid=mid, A=A, T=T, AT=A + T)


def K_continuous(S, omega, gamma=1.0):
    n, mid = S["n"], S["mid"]
    y = np.linalg.solve(1j * omega * gamma * np.eye(n) + S["AT"], S["kt"] * np.eye(n)[mid])
    return S["kt"] / y[mid] - S["kt"] - 1j * omega * gamma, y


def K_euler(S, omega, dt, gamma=1.0):
    """The **exact** steady response of explicit Euler.

    Substituting y_n = Im[y_hat e^{i*omega*n*dt}] gives
    (e^{iωdt} I − M) ŷ = B â,  M = I − (A+T)dt/γ,  B = (k_t dt/γ) e_mid."""
    n, mid = S["n"], S["mid"]
    M = np.eye(n) - S["AT"] * dt / gamma
    B = (S["kt"] * dt / gamma) * np.eye(n)[mid]
    y = np.linalg.solve(math.e ** (1j * omega * dt) * np.eye(n) - M, B)
    return S["kt"] / y[mid] - S["kt"] - 1j * omega * gamma, y


def stage1(S):
    print("=" * 92)
    print("(1) dt (explicit Euler) -- settled by the exact z-domain solution, "
          "without simulating")
    print("=" * 92)
    print(f"dt = {S['dt']:.4e}   dt·λ_max/γ = {S['dt']*np.linalg.eigvalsh(S['AT'])[-1]:.4f}")
    print(f"{'De':>8}{'K1 cont':>12}{'K1 Euler':>12}{'diff%':>8}"
          f"{'K2 cont':>12}{'K2 Euler':>12}{'diff%':>8}")
    print("-" * 92)
    worst = 0.0
    for sp in S["specs"]:
        om = float(sp["params"]["omega_star"]); de = float(sp["params"]["De"])
        Kc, _ = K_continuous(S, om)
        Ke, _ = K_euler(S, om, S["dt"])
        dre = 100 * (Ke.real / Kc.real - 1); dim = 100 * (Ke.imag / Kc.imag - 1)
        worst = max(worst, abs(dre), abs(dim))
        print(f"{de:>8.3f}{Kc.real:>12.1f}{Ke.real:>12.1f}{dre:>8.3f}"
              f"{Kc.imag:>12.1f}{Ke.imag:>12.1f}{dim:>8.3f}")
    print("-" * 92)
    print(f"largest difference {worst:.3f}%  ->  "
          f"{'dt is NOT the cause (it cannot explain 28%)' if worst < 1 else 'dt does contribute'}")
    print("=" * 92)
    return worst


def stage2(S):
    """The static stiffness with the traps on, in HOOMD.

    kT=0, large dt, and the ghost displaced by delta and held there.
    """
    n, kt, ell, mid = S["n"], S["kt"], S["ell"], S["mid"]
    delta = S["amp"]
    lam = np.linalg.eigvalsh(S["AT"])
    dt = 0.2 / lam[-1]
    n_steps = int(round(25.0 / (lam[0] * dt)))          # 25 τ_max

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in S["idx"]:
        pos.append(list(pos[g])); typeid.append(1)
    pos[n + S["idx"].index(mid)][1] = delta              # only the driven ghost is displaced by delta and held

    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.configuration.box = [40.0 * n] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(S["idx"])]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(S["idx"])
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=2)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=S["kb"], r0=ell)
    bond.params["trap"] = dict(k=kt, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=S["kth"], t0=math.pi)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Type(["A"]), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[ov], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(n_steps)

    y = np.array(sim.state.get_snapshot().particles.position)[:n, 1]
    K_hd = kt * (delta / y[mid] - 1.0)                   # static: K = k_t(y_c/y_mid - 1)
    y_mod = np.linalg.solve(S["AT"], kt * delta * np.eye(n)[mid])
    K_mod = kt * (delta / y_mod[mid] - 1.0)

    print("=" * 92)
    print("(2) static + traps -- the static stiffness with traps on ((1) verified "
          "'rigidly clamped, no traps')")
    print("=" * 92)
    print(f"dt={dt:.3e}  steps={n_steps:,}  delta={delta:.5f}")
    print(f"model  K_static = {K_mod:10.2f}")
    print(f"HOOMD K_static = {K_hd:10.2f}   ({100*(K_hd/K_mod-1):+.2f}%)")
    print(f"y_mid: model {y_mod[mid]:.6f}, HOOMD {y[mid]:.6f}  "
          f"({100*(y[mid]/y_mod[mid]-1):+.3f}%)")
    print(f"largest difference in the displacement profile = "
          f"{100*np.max(np.abs(y-y_mod))/delta:.3f}% of delta")
    ok = abs(K_hd / K_mod - 1) < 0.02
    print("-" * 92)
    print(f"-> {'the static case agrees even with the traps on, so the disagreement is **dynamical**' if ok else 'it already disagrees statically -- a trap/boundary problem'}")
    print("=" * 92)
    return K_hd, K_mod, y, y_mod


def stage3(S, de_target, n_cycles=3, eq_tau=12.0):
    """Dynamic mode shape -- measure the complex phasor of all 25 beads and compare
    against the model."""
    sp = min(S["specs"], key=lambda s: abs(s["params"]["De"] - de_target))
    p = sp["params"]
    om = float(p["omega_star"]); n, mid, ell = S["n"], S["mid"], S["ell"]
    dt = S["dt"]; amp = S["amp"]; kt = S["kt"]
    lam = np.linalg.eigvalsh(S["AT"])

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in S["idx"]:
        pos.append(list(pos[g])); typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.configuration.box = [40.0 * n] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(S["idx"])]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(S["idx"])
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=4)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=S["kb"], r0=ell)
    bond.params["trap"] = dict(k=kt, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=S["kth"], t0=math.pi)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.Type(["A"]), default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[ov], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    ghost_mid = n + S["idx"].index(mid)
    UPD = 100

    class Move(hoomd.custom.Action):
        def act(self, ts):
            y = amp * math.sin(om * ts * dt)
            with self._state.cpu_local_snapshot as s:
                tg = np.array(s.particles.tag, copy=True)
                loc = np.flatnonzero(tg == ghost_mid)
                if len(loc):
                    s.particles.position[loc[0], 1] = y

    class All(hoomd.custom.Action):
        def __init__(self): self.t, self.Y, self.g = [], [], []
        def act(self, ts):
            with self._state.cpu_local_snapshot as s:
                tg = np.array(s.particles.tag, copy=True)
                q = np.array(s.particles.position, copy=True)
                o = np.argsort(tg)
                self.t.append(ts * dt); self.Y.append(q[o][:n, 1].copy())
                self.g.append(q[o][ghost_mid, 1])

    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=Move(), trigger=hoomd.trigger.Periodic(UPD)))
    period = 2 * math.pi / om
    spc = max(20, int(round(period / dt)))
    smp = All()
    sim.operations.writers.append(hoomd.write.CustomWriter(
        action=smp, trigger=hoomd.trigger.Periodic(max(1, spc // 20))))
    n_eq = int(round(eq_tau / (lam[0] * dt)))
    sim.run(n_eq)
    smp.t.clear(); smp.Y.clear(); smp.g.clear()
    sim.run(int(round(n_cycles * period / dt)))

    t = np.array(smp.t); Y = np.array(smp.Y); g = np.array(smp.g)
    def lk(s):
        return complex(2 * np.mean(s * np.sin(om * t)), 2 * np.mean(s * np.cos(om * t)))
    y_hd = np.array([lk(Y[:, i]) for i in range(n)])
    g_hat = lk(g)
    y_mod = np.linalg.solve(1j * om * np.eye(n) + S["AT"], kt * g_hat * np.eye(n)[mid])

    print("=" * 92)
    print(f"(3) dynamic mode shape -- De={float(p['De']):.2f}  "
          f"(complex phasors of all 25 beads)")
    print("=" * 92)
    print(f"equilibration {n_eq:,} steps = {eq_tau:g} tau_max, {len(t)} samples")
    print(f"{'bead':>5}{'|y| model':>12}{'|y| HOOMD':>12}{'ratio':>8}"
          f"{'ph model':>11}{'ph HOOMD':>12}{'diff[rad]':>10}")
    print("-" * 92)
    for i in range(n):
        if i in S["idx"] or i % 4 == 0:
            mark = " <-trap" if i in S["idx"] else ""
            print(f"{i:>5}{abs(y_mod[i]):>12.6f}{abs(y_hd[i]):>12.6f}"
                  f"{abs(y_hd[i])/abs(y_mod[i]):>8.3f}"
                  f"{np.angle(y_mod[i]):>11.4f}{np.angle(y_hd[i]):>12.4f}"
                  f"{np.angle(y_hd[i])-np.angle(y_mod[i]):>10.4f}{mark}")
    print("-" * 92)
    rat = np.abs(y_hd) / np.abs(y_mod)
    print(f"amplitude-ratio range {rat.min():.3f} to {rat.max():.3f}   "
          f"{'★ every bead has the same ratio -> an overall scale problem' if rat.max()/rat.min() < 1.05 else '★ it varies per bead -> the mode SHAPE differs'}")
    Kh = kt * g_hat / y_hd[mid] - kt - 1j * om
    Km = kt * g_hat / y_mod[mid] - kt - 1j * om
    print(f"K* model ({Km.real:.0f}, {Km.imag:.0f})  .  "
          f"HOOMD ({Kh.real:.0f}, {Kh.imag:.0f})"
          f"  → K′ {100*(Kh.real/Km.real-1):+.1f}%  K″ {100*(Kh.imag/Km.imag-1):+.1f}%")
    print("=" * 92)
    np.savez(GATES / f"modeshape_de{float(p['De']):.1f}.npz",
             y_hd=y_hd, y_mod=y_mod, g_hat=g_hat, omega=om, de=float(p["De"]))
    return y_hd, y_mod


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="12")
    ap.add_argument("--de", type=float, default=91.8)
    ap.add_argument("--eq-tau", type=float, default=12.0)
    a = ap.parse_args()
    S = setup()
    if "1" in a.stage:
        stage1(S)
    if "2" in a.stage:
        stage2(S)
    if "3" in a.stage:
        stage3(S, a.de, eq_tau=a.eq_tau)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
