"""Golden test for the moving trap -- is `bdbot.traps.make_trap(velocity=...)` right?

CLAUDE.md rule 7: **verify independent elements one at a time.**
`trap-drag-2d-hex300` combines two elements -- a moving trap and a soft hexagonal
lattice -- and the combination has no analytic solution.
**The minimal configuration with only the moving trap on** does.

An overdamped particle dragged at constant velocity:
    γẋ = −k(x − vt) + ξ ,   u ≡ x − vt
    ⟹ γu̇ = −k u − γv + ξ
    => <u> = -gamma*v/k    (it lags. The NEGATIVE sign is the point -- the combined
                            code got this wrong)
       Var(u) = kT/k       ★ **dragging does not change the variance** (drift moves
                            the mean only)
       Var(y) = kT/k

In reduced units (gamma=kT=1): <dx> = -v*/k*, Var = 1/k*.

Two things are targeted:
  (1) **sign and magnitude** -- the combined code wrote `F_drag = +k<dx>` and produced
      -493% against bare Stokes. The sign convention is pinned down here.
  (2) **minimum image (trap 1)** -- the trap centre is `anchor0 + v*t`, so it crosses
      the box **over and over**. With a fixed trap only the particle wrapped; here
      **the anchor recedes without bound.**
      Wrong wrapping breaks silently for a weak trap -- which is why weak k is used.

    $PY scratch/verify_moving_trap.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import sim as SIM, traps as TR  # noqa: E402

L = 40.0
N = 1000                 # independent particles -- N times the statistics
                         # (no interactions)
SEED = 20260804


def run(k_star, v_star, n_tau=200.0, n_samples=200, dt_over_tau=1e-2, seed=SEED):
    """Drag N independent particles, each with its own moving trap. (a dict of
    measurements)

    ★ `seed` MUST differ between combinations. With gamma=kT=1 and dt = h*tau_k, the
      dynamics of u in units of tau_k is **identical regardless of k**, so reusing one
      seed makes all four rows the same trajectory -- not four independent checks but
      one check printed four times. (That actually happened, and it was noticed
      because all four rows agreed to the decimal at +0.15%.)
      ⚠️ Trap 12: HOOMD seeds are truncated to 16 bits, so small consecutive integers
      are used.
    """
    tau_k = 1.0 / k_star                      # γ/k, γ=1
    dt = dt_over_tau * tau_k
    n_steps = int(round(n_tau * tau_k / dt))
    every = max(1, n_steps // n_samples)

    rng = np.random.default_rng(seed)
    pos = rng.uniform(-L / 2, L / 2, (N, 2))
    sim = SIM.make_sim(SIM.frame_2d(pos, L), seed=seed)

    vel = np.zeros((N, 3))
    vel[:, 0] = v_star
    trap = TR.make_trap(k_star, pos, L, dt_star=dt, velocity=vel)
    SIM.attach_brownian(sim, dt, [trap])     # ★ no pair force -- the minimal
                                             # configuration, trap only

    sim.run(int(20 * tau_k / dt))             # equilibrate for 20 tau_k
    dx, dy = [], []
    for _ in range(n_samples):
        sim.run(every)
        d = trap.displacement(sim.state, sim.timestep)
        dx.append(d[:, 0].copy())
        dy.append(d[:, 1].copy())
    dx, dy = np.array(dx), np.array(dy)

    n_ind = N * n_tau / 2.0                   # correlation time 2 tau_k -> number of
                                              # independent samples
    return {
        "mean_x": float(dx.mean()), "sem_x": float(dx.std() / math.sqrt(n_ind)),
        "var_x": float(dx.var()), "var_y": float(dy.var()),
        "traverse": v_star * n_tau * tau_k / L,      # how many times the trap crossed
                                                     # the box
    }


print("=" * 92)
print("moving-trap golden test -- <dx> = -v/k, Var = 1/k "
      "(dragging leaves the variance unchanged)")
print("=" * 92)
print(f"  N={N} independent particles . box L={L:g} . no pair force . "
      f"dt = 1e-2 tau_k")
print("  ★ The verdict is **how many sigma against its own statistical error (SEM)**, "
      "not a fixed % --")
print("    in a system with SNR=0.0985, 'within 2%' is a demand the statistics never "
      "permitted.")
print(f"\n{'k*':>7}{'v*':>7}{'SNR':>6}{'cross':>6} | "
      f"{'<dx> meas':>13}{'pred -v/k':>12}{'error':>8}{'sigma':>7} | "
      f"{'Var_x.k':>9}{'Var_y.k':>9}{'theory':>8}")
rows, ok = [], True
# ★ The weaker the trap, the more vulnerable to a minimum-image defect
#   (trap 1: +1856% at k=2).
#   SNR = v/sqrt(k) is kept near 1 so the statistics work out. Only the last row uses
#   the real case values.
VAR_TH = 1.0 / (1.0 - 1e-2 / 2)          # discrete-OU stationary variance bias
                                         # (bd-physics §1.2)
for i, (k_star, v_star) in enumerate(((2.0, 1.5), (5.0, 2.5), (10.0, 3.5),
                                      (60358.0, 24.205))):
    r = run(k_star, v_star, seed=101 + i)          # ★ a different seed per row
                                                   # (see the docstring above)
    pred = -v_star / k_star
    err = 100 * (r["mean_x"] - pred) / abs(pred)
    nsig = abs(r["mean_x"] - pred) / r["sem_x"]
    vx, vy = r["var_x"] * k_star, r["var_y"] * k_star
    # relative statistical error on a variance ~ sqrt(2/n_indep)
    v_tol = 4 * math.sqrt(2 / (N * 200.0 / 2))
    good = nsig < 3.0 and abs(vx - VAR_TH) < v_tol and abs(vy - VAR_TH) < v_tol
    ok &= good
    rows.append((k_star, err, nsig, vx, vy, good))
    print(f"{k_star:>7g}{v_star:>7g}{v_star/math.sqrt(k_star):>6.2f}"
          f"{r['traverse']:>6.1f} | {r['mean_x']:>13.6f}{pred:>12.6f}{err:>+7.2f}%"
          f"{nsig:>6.1f}σ | {vx:>9.4f}{vy:>9.4f}{VAR_TH:>8.4f}   {'✓' if good else '✗'}")

# ── dt convergence -- is the residual systematic error O(dt) discretisation, or a
#    bug? ────────────
#   Above, k=2, 5 and 10 all showed a **constant** systematic error of +0.15%. The
#   size is small but the constancy is what is worrying -- halve dt, and if the error
#   halves too it is discretisation; if it does not move, it is a convention bug.
#   Settled by running rather than guessing (CLAUDE.md rule 6).
#   ⭐️ **The variance** is where the dt check belongs -- NOT the mean. The mean of a
#      discrete OU process is `<u> = -v*gamma/k`, **exact regardless of dt** (derived
#      above), and the bias enters the variance alone as `1/(1-h/2)`.
#      The variance is also sharper: its relative statistical error is
#      sqrt(2/n) ~ 0.45%, against 0.28% of signal for the mean.
print("\n  dt convergence (k*=5, v*=2.5) -- does the variance bias follow the "
      "theoretical 1/(1-h/2)?")
print(f"    {'dt/tau_k':>8}{'Var_x.k':>10}{'theory':>9}{'diff':>9}{'stat err':>9}")
conv_ok = True
for j, hh in enumerate((4e-2, 2e-2, 1e-2, 5e-3)):
    r = run(5.0, 2.5, dt_over_tau=hh, seed=211 + j)
    th = 1.0 / (1.0 - hh / 2)
    got = r["var_x"] * 5
    tol = 3 * math.sqrt(2 / (N * 200.0 / 2))
    good = abs(got - th) < tol
    conv_ok &= good
    print(f"    {hh:>8.0e}{got:>10.4f}{th:>9.4f}{got-th:>+9.4f}{tol:>9.4f}   "
          f"{'✓' if good else '✗'}")
ok &= conv_ok
print(f"    -> {'✓ the variance bias matches theory -- dragging leaves the variance at kT/k' if conv_ok else '✗ it disagrees with theory'}")

print(f"""
  Criteria: <dx> within 3 sigma of prediction, and the Var bias within statistical
  error of the theoretical 1/(1-h/2).
  ⚠️ The first version used a fixed criterion like 'error < 2%' AND **the same seed on
     every row**, and all four rows came out at +0.15%, agreeing to the decimal -- not
     four independent checks but one check printed four times. Halving dt did not move
     it either, so it was nearly misdiagnosed as a bug.
  ★ The last row uses the real case values (k*=6.04e4, v*=24.2, SNR=0.0985) -- the
    signal is 1/10 of the noise, so it is only barely visible after averaging
    N={N} particles. The statistics problem L3 predicted shows up here too.""")
print("=" * 92)
print(f"{'✓ PASS' if ok else '✗ FAIL'} -- moving trap "
      f"{sum(1 for r in rows if r[4])}/{len(rows)}")
print("=" * 92)
sys.exit(0 if ok else 1)
