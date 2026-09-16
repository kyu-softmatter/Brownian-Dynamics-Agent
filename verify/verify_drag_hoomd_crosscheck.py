"""HOOMD cross-check: is the exact-OU reduction the right model of the BD system?

`verify_drag_ladder.py` measures the per-rung precision of the drag-slope route on
exact Ornstein-Uhlenbeck trajectories. That is a *statistics* argument and it
assumes the BD system reduces to OU. This run checks that assumption against the
actual integrator, which is a different kind of evidence (A1): the OU result is
the analytic layer, this is self-consistency against the thing being modelled.

Approved parameter manifest: verify/_out/params_trap_drag_ladder.txt
No camera here -- the detector is a post-processing layer and is checked
separately in verify_detector_model.py. What is checked is BD == OU.

    $PY verify/verify_drag_hoomd_crosscheck.py
"""
from __future__ import annotations
import json, math, sys, time
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from bdbot import sim as SIM, traps as TR          # noqa: E402

K_STAR   = 21221.4                  # k_t d^2 / kT, from the manifest
DT_STAR  = 9.42443e-08              # dt/tau_k = 0.002
N_EQ     = 10_000
N_PROD   = 161_436                  # T_obs/tau_k = 322.9
SAMPLE   = 50                       # 10 samples per tau_k
L_STAR   = 32.0
N        = 1000
SEED     = 11                       # bd-hoomd trap 12: small integer
V_FRACS  = (0.25, 0.50, 0.75, 1.00)
X_EQ_MAX = 0.075                    # 0.15a/d


def main():
    n_side = int(math.ceil(math.sqrt(N)))
    a = L_STAR / n_side
    pos0 = np.array([[(i % n_side + .5) * a - L_STAR / 2,
                      (i // n_side + .5) * a - L_STAR / 2, 0.0] for i in range(N)])
    # four velocity groups, 250 beads each -- one run, four velocities
    grp = np.arange(N) % len(V_FRACS)
    v_of = np.array([f * X_EQ_MAX * K_STAR for f in V_FRACS])   # v* = x_eq* k*
    vel = np.zeros((N, 3)); vel[:, 0] = v_of[grp]

    sim = SIM.make_sim(SIM.frame_2d(pos0, L_STAR), seed=SEED)
    trap = TR.make_trap(K_STAR, pos0, L_STAR, dt_star=DT_STAR, velocity=vel)
    SIM.attach_brownian(sim, DT_STAR, [trap])

    t0 = time.time()
    sim.run(N_EQ)
    n_s = N_PROD // SAMPLE
    acc_x = np.zeros(N); acc_x2 = np.zeros(N)
    for i in range(n_s):
        sim.run(SAMPLE)
        dx = trap.displacement(sim.state, sim.timestep)[:, 0]
        acc_x += dx; acc_x2 += dx * dx
    wall = time.time() - t0

    mean = acc_x / n_s
    var = acc_x2 / n_s - mean ** 2
    if not np.isfinite(mean).all():
        raise RuntimeError("non-finite displacement -- run is invalid")

    print("=" * 84)
    print("HOOMD cross-check of the OU reduction -- trap + constant-velocity drive")
    print("=" * 84)
    print(f"  k* = {K_STAR}  dt/tau_k = {DT_STAR*K_STAR:.4f}  T_obs/tau_k = "
          f"{N_PROD*DT_STAR*K_STAR:.1f}  N = {N}  wall = {wall/60:.1f} min")
    print(f"\n  {'v*':>10} {'<x>/d measured':>22} {'analytic -v*/k*':>16} {'err':>9}"
          f"  {'var(x)/d^2':>13} {'analytic 1/k*':>14} {'err':>9}")
    out = {"k_star": K_STAR, "dt_over_tau_k": DT_STAR * K_STAR,
           "T_obs_over_tau_k": N_PROD * DT_STAR * K_STAR, "wall_s": wall, "groups": []}
    for g, v in enumerate(v_of):
        m = mean[grp == g]; vv = var[grp == g]
        pred_m = -v / K_STAR; pred_v = 1.0 / K_STAR
        em = 100 * (m.mean() / pred_m - 1); ev = 100 * (vv.mean() / pred_v - 1)
        sem_m = m.std(ddof=1) / math.sqrt(len(m)); sem_v = vv.std(ddof=1) / math.sqrt(len(vv))
        out["groups"].append(dict(v_star=float(v), n=int(len(m)),
                                  mean=float(m.mean()), mean_sem=float(sem_m),
                                  mean_pred=float(pred_m), mean_err_pct=float(em),
                                  var=float(vv.mean()), var_sem=float(sem_v),
                                  var_pred=float(pred_v), var_err_pct=float(ev)))
        print(f"  {v:10.1f} {m.mean():+14.6e}+/-{sem_m:.0e} {pred_m:+16.6e} {em:+8.3f} %"
              f"  {vv.mean():13.6e} {pred_v:14.6e} {ev:+8.3f} %")

    # the slope, which is what the ladder route actually uses
    A = np.vstack([v_of, np.ones_like(v_of)]).T
    y = np.array([o["mean"] for o in out["groups"]])
    slope = np.linalg.lstsq(A, y, rcond=None)[0][0]
    out["slope_measured"] = float(slope); out["slope_pred"] = float(-1.0 / K_STAR)
    out["slope_err_pct"] = float(100 * (slope / (-1.0 / K_STAR) - 1))
    print(f"\n  slope d<x>/dv  = {slope:+.6e}   analytic -1/k* = {-1/K_STAR:+.6e}"
          f"   err {out['slope_err_pct']:+.3f} %")
    p = ROOT / "verify/_out/drag_hoomd_crosscheck.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
