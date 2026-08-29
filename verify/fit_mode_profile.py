"""Per-bead A*sin(omega*t+phi) fit for the JKR chain -> amplitude and phase profile
versus distance from the centre.

Each bead's y(t) is lock-in detected at the drive frequency to get the complex phasor
y_hat_i; then A_i=|y_hat_i| and phi_i=arg(y_hat_i) are tabulated against the distance
s=|i-mid| from the centre (the driven bead).

★ **There is an analytic prediction** -- linear response:
  (i*omega*gamma*I + A_bend + T) y_hat = k_t*a*e_mid
  In this system the bending matrix A_bend and the traps T are fully known, so y_hat
  can be solved exactly. Comparing with the measurement verifies the implementation
  and the physics **bead by bead** (so far only the driven bead had been checked --
  |y_hat_mid| was +1.3% against prediction).

  Phases are referenced to the driven bead (phi_mid == 0). Absolute phase is
  contaminated by the ZOH and by sampling delay, but **the difference cancels**
  (same logic as bd-hoomd trap 17).

    $PY scratch/fit_mode_profile.py
"""
import glob
import json
import sys
from pathlib import Path

import gsd.hoomd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))
from bdbot import lockin as LI            # noqa: E402
from chain_bend_dlvo_2d import bending_matrix, trapped_indices  # noqa: E402

PAT = "runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr-kt100__*"


def phasors_from_npz(d):
    """Per-bead phasor from every sample in observables.npz.

    Denser sampling than GSD.
    """
    z = np.load(Path(d) / "observables.npz", allow_pickle=True)
    keys = set(z.files)
    if "shape_y" not in keys:
        return None, None
    t = np.asarray(z["t"], dtype=float)
    ys = np.asarray(z["shape_y"], dtype=float)      # (T, n)
    return t, ys


def phasors_from_gsd(d, n):
    tr = gsd.hoomd.open(str(Path(d) / "traj_A.gsd"), "r")
    ys = np.array([tr[i].particles.position[:n, 1] for i in range(len(tr))])
    s = json.loads((Path(d) / "spec.json").read_text())
    t = np.linspace(0, s["numerics"]["n_prod"] * s["numerics"]["dt_star"], len(ys))
    return t, ys


rows = []
spec0 = None
for d in sorted(glob.glob(PAT)):
    if not (Path(d) / "metrics.json").exists():
        continue
    s = json.loads((Path(d) / "spec.json").read_text())
    spec0 = spec0 or s
    n = int(s["params"]["n_beads"])
    om = float(s["params"]["omega_star"])
    t, ys = phasors_from_npz(d)
    src = "npz"
    if t is None:
        t, ys = phasors_from_gsd(d, n)
        src = "gsd"
    ph = []
    for i in range(n):
        blocks = LI.lockin_blocks(t, ys[:, i] - ys[:, i].mean(), om,
                                  n_blocks=min(10, max(2, len(t) // 20)))
        h, _ = LI.agg(blocks)
        ph.append(h)
    rows.append(np.array(ph))
    print(f"  {Path(d).name[-12:]}  seed={s['numerics']['seed']}  sample source={src} ({len(t)})")

Z = np.array(rows)                      # (n_seeds, n)
n = Z.shape[1]
mid = trapped_indices(n)[1]
# Phasors relative to the driven bead (cancels the ZOH and sampling delay in the
# absolute phase)
Zrel = Z / Z[:, [mid]]

A = np.abs(Z)
A_m, A_e = A.mean(0), A.std(0, ddof=1) / np.sqrt(len(A))
phi = np.angle(Zrel)
phi_m = np.angle(Zrel.mean(0))
phi_e = np.abs(Zrel).std(0, ddof=1) / np.sqrt(len(Zrel)) / np.abs(Zrel.mean(0)).clip(1e-30)

# ── analytic linear-response prediction ──
P = spec0["params"]
kth, k_t = float(P["kappa_theta_star"]), float(P["k_t_star"])
ell = 1.0 + float(P["h_min_star"])
amp, om = float(P["amp_star"]), float(P["omega_star"])
Ab = bending_matrix(n, kth, ell)
for i in trapped_indices(n):
    Ab[i, i] += k_t
e = np.zeros(n, dtype=complex); e[mid] = k_t * amp
ypred = np.linalg.solve(1j * om * np.eye(n) + Ab, e)
Apred = np.abs(ypred)
phipred = np.angle(ypred / ypred[mid])

s_idx = np.arange(n) - mid
print()
print("=" * 88)
print(f"JKR per-bead A*sin(omega*t+phi) fit  ({len(Z)} seeds, n={n}, "
      f"driven = bead {mid})")
print("=" * 88)
print(f"{'bead':>4} {'s=i-mid':>8} {'A meas [d]':>18} {'A pred':>10} {'diff%':>8} "
      f"{'phi meas [deg]':>13} {'phi pred':>9} {'diff':>7}")
for i in range(n):
    dA = 100 * (A_m[i] - Apred[i]) / Apred[i]
    dphi = np.degrees(np.angle(np.exp(1j * (phi_m[i] - phipred[i]))))
    print(f"{i:>4} {s_idx[i]:>8} {A_m[i]:>11.5f}±{A_e[i]:<6.4f} {Apred[i]:>10.5f} "
          f"{dA:>+8.2f} {np.degrees(phi_m[i]):>13.2f} {np.degrees(phipred[i]):>9.2f} {dphi:>+7.2f}")

np.savez("/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/"
         "7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/mode_profile.npz",
         A_m=A_m, A_e=A_e, phi_m=phi_m, phi_e=phi_e,
         Apred=Apred, phipred=phipred, s_idx=s_idx, mid=mid, n=n)
print("\nsaved: mode_profile.npz")
