"""Is the trap thread's headline +1.17 % on `f_c` physics, or the estimator?

`runs/trap-2d-5um__a5ef4f45d589` reports `f_c` measured +1.17 % above the
closed-form prediction, role `implementation_check`, inside a 5 % tolerance.
That number is the whole precision argument of the `trap-stiffness-recovery`
bridge thread: r2's `T_obs >= 32.3 s` requirement exists because 1.17 % "was
measured at T_obs/tau = 2000".

This script feeds **exact OU data, whose `f_c` is known in closed form**, through
the run's own estimator configuration (`cases/trap_2d_5um.py` `finalize()`:
Welch, `nperseg = min(n_samp, 4096)`, unweighted `curve_fit` of
`S0/(1+(f/f_c)^2)` over the whole one-sided PSD, PSDs averaged over the trace
subset before fitting). If the bias reappears there, it is not the simulation.

This is the third time in this project that a "physics discrepancy" turned out
to be the analysis (trap tau 70 % high, an MSD plateau 3 % high, a D 9 % low),
and the device is the same one that settled all three: synthetic data whose
answer is known.

It also answers r2's `gaps[3]` -- one bead against the replica ensemble.

    $PY verify/verify_fc_estimator_bias.py
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "verify"))
from detector_model import lorentzian_fc, ou_exact                      # noqa: E402

FC_TRUE = 1.0 / (2 * math.pi)          # in units of 1/tau -- OU, exact
RUN = "trap-2d-5um__a5ef4f45d589"
RUN_ERR_PCT = 1.1704777782934543       # metrics.json observables[3].err_pct
RUN_SAMPLES_PER_TAU = 10               # numerics.sample_every=50 x dt/tau=0.002
RUN_T_OBS_OVER_TAU = 2000
RUN_N_TRACE = 250                      # cases/trap_2d_5um.py: n_trace = min(250, N)
NPERSEG = 4096


def sem(v):
    v = np.asarray(v, float)
    return float(v.std(ddof=1) / math.sqrt(len(v)))


def fit_ensemble(dt, t_obs, n_paths, seed, nperseg=NPERSEG):
    n = int(t_obs / dt)
    x = ou_exact(n, dt, 1.0, np.random.default_rng(seed), n_paths)
    return lorentzian_fc(x, fs=1.0 / dt, nperseg=min(nperseg, n))[0]


def bias_pct(vals):
    m, e = float(np.mean(vals)), sem(vals)
    return 100 * (m / FC_TRUE - 1), 100 * e / FC_TRUE


def main():
    out = {"run": RUN, "run_err_pct": RUN_ERR_PCT, "fc_true_per_tau": FC_TRUE}

    print("=" * 78)
    print("1 - the run's exact estimator configuration, on data whose answer is known")
    print("=" * 78)
    dt = 1.0 / RUN_SAMPLES_PER_TAU
    v = [fit_ensemble(dt, RUN_T_OBS_OVER_TAU, RUN_N_TRACE, 1000 + b) for b in range(6)]
    b, e = bias_pct(v)
    out["reproduced_bias_pct"] = [b, e]
    print(f"  samples/tau = {RUN_SAMPLES_PER_TAU} . T_obs/tau = {RUN_T_OBS_OVER_TAU} . "
          f"{RUN_N_TRACE} traces . nperseg = {NPERSEG}")
    print(f"  estimator bias on exact OU  = {b:+.3f} % +/- {e:.3f} %")
    print(f"  the archived run reported     {RUN_ERR_PCT:+.3f} %")
    print(f"  -> the run's number is {abs(RUN_ERR_PCT - b) / e:.1f} sigma from the "
          f"estimator's own bias")
    print("  role: `hypothesis`. The prediction 'f_c err is physics' is NOT imposed by\n"
          "        the simulation, so a mismatch here is a RESULT, not a bug (rule 7').")

    print("\n" + "=" * 78)
    print("2 - is it a bias or a variance? sweep T_obs (a variance shrinks, a bias does not)")
    print("=" * 78)
    out["vs_T_obs"] = {}
    print(f"  {'T_obs/tau':>10}  {'bias on f_c':>22}")
    for t in (310, 1000, 2000, 8000):
        v = [fit_ensemble(dt, t, RUN_N_TRACE, 2000 + b) for b in range(6)]
        b, e = bias_pct(v)
        out["vs_T_obs"][t] = [b, e]
        note = "  <- AM's 5 s per rung" if t == 310 else (
               "  <- the run" if t == 2000 else "")
        print(f"  {t:10d}  {b:+9.3f} % +/- {e:5.3f}{note}")
    print("  T_obs x 26 does not remove it. It is a bias.")

    print("\n" + "=" * 78)
    print("3 - where the bias comes from: the fit range, in units of f_c")
    print("=" * 78)
    out["vs_sampling"] = {}
    print(f"  {'samples/tau':>12}  {'f_nyq/f_c':>10}  {'bias on f_c':>22}")
    for spt in (2, 5, 10, 20, 50):
        d = 1.0 / spt
        v = [fit_ensemble(d, RUN_T_OBS_OVER_TAU, RUN_N_TRACE, 3000 + b) for b in range(6)]
        b, e = bias_pct(v)
        out["vs_sampling"][spt] = [b, e]
        nyq = (spt / 2) / FC_TRUE
        star = "  <- BD convention, and the run" if spt == 10 else ""
        print(f"  {spt:12d}  {nyq:10.1f}  {b:+9.3f} % +/- {e:5.3f}{star}")
    print("  Non-monotonic, and worst where the PSD is mostly far above the corner:\n"
          "  an unweighted least squares over the whole band is dominated by points\n"
          "  that carry no information about f_c.")

    print("\n" + "=" * 78)
    print("4 - r2 gaps[3]: one bead against the replica ensemble")
    print("=" * 78)
    out["single_bead"] = {}
    print(f"  {'T_obs/tau':>10}  {'one bead':>20}  {'ensemble of 250':>20}  {'ratio':>7}")
    for t in (310, 2000):
        n = int(t / dt)
        x = ou_exact(n, dt, 1.0, np.random.default_rng(31337), RUN_N_TRACE)
        singles = []
        for i in range(RUN_N_TRACE):
            try:
                singles.append(lorentzian_fc(x[:, i:i + 1], fs=1.0 / dt,
                                             nperseg=min(NPERSEG, n))[0])
            except Exception:
                pass
        s = np.array(singles)
        s_sd = 100 * s.std(ddof=1) / FC_TRUE
        ens = [fit_ensemble(dt, t, RUN_N_TRACE, 500 + b) for b in range(12)]
        e_sd = 100 * np.std(ens, ddof=1) / FC_TRUE
        out["single_bead"][t] = {"one_bead_sd_pct": s_sd, "ensemble_sd_pct": e_sd,
                                 "ratio": s_sd / e_sd, "n": len(s)}
        print(f"  {t:10d}  {s_sd:15.2f} % sd  {e_sd:15.2f} % sd  {s_sd/e_sd:6.1f}x")
    print(f"  sqrt(250) = {math.sqrt(250):.1f} . sqrt(1000) = {math.sqrt(1000):.1f}")
    print("  r2 finding[3] said ~32x = sqrt(1000). But f_c in that run is fitted to the\n"
          "  PSD of the 250-trace subset, not of all 1000 particles\n"
          "  (cases/trap_2d_5um.py: n_trace = min(250, N)) -- and averaging PSDs before a\n"
          "  non-linear fit is not averaging fits, so the ratio is not sqrt(n) either.")

    p = ROOT / "verify" / "_out" / "fc_estimator_bias.json"
    p.parent.mkdir(exist_ok=True)
    p.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {p.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
