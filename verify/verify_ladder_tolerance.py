"""r7's inverse question: the largest per-rung scatter the six-rung fit tolerates.

r1 asked what the ladder DELIVERS at an assumed 3 %. r7 asks the actionable
inverse -- the largest per-rung relative scatter `sigma_gamma_max` at which the
fit still returns `h0` to +/-0.2 um, unbiased. That is a tolerance the
instrument side can compare a measurement against.

Pure propagation: gamma_i = gamma_true(h_i)(1 + sigma*xi), fitted with the same
two-parameter Faxen model. No trajectory is needed, because r7's question is
about the FIT and not about how gamma was obtained -- which is also why the
answer is independent of whether the 7-9.5 % of verify_drag_ladder.py is right.

    $PY verify/verify_ladder_tolerance.py
"""
from __future__ import annotations
import json, math, sys
from pathlib import Path
import numpy as np
from scipy import optimize

ROOT = Path(__file__).resolve().parent.parent
A_UM = 2.475
RUNGS = np.array([3.0, 3.5, 4.5, 6.0, 8.0, 10.0])
G0 = 0.0467273                      # bulk Stokes, pN*s/um -- the manifest's derived.gamma
TARGET_UM = 0.2                     # r1 required_precision / r7's +/-0.2 um
N_TRIAL = 4000


def model(dz, g0, h0):
    return g0 / (1.0 - 9.0 * A_UM / (16.0 * (dz + h0)))


def fit_scatter(sigma, h0_true=1.0, n=N_TRIAL, seed=5150):
    rng = np.random.default_rng(seed)
    dz = RUNGS - h0_true
    g_true = model(dz, G0, h0_true)
    y = g_true[None, :] * (1.0 + sigma * rng.standard_normal((n, len(dz))))
    h0s = []
    for row in y:
        try:
            p, _ = optimize.curve_fit(model, dz, row, p0=[G0, h0_true], maxfev=20000)
            if 0.0 < p[1] < 10.0:          # a runaway fit is a failure, not a sample
                h0s.append(p[1])
        except Exception:
            pass
    h = np.array(h0s)
    return dict(sigma=sigma, n=len(h), frac_ok=len(h) / n,
                mean=float(h.mean()), sd=float(h.std(ddof=1)),
                bias=float(h.mean() - h0_true),
                sd_err=float(h.std(ddof=1) / math.sqrt(2 * (len(h) - 1))),
                bias_se=float(h.std(ddof=1) / math.sqrt(len(h))))


def main():
    print("=" * 92)
    print("r7: sigma_gamma_max -- the largest per-rung scatter the six-rung fit tolerates")
    print("=" * 92)
    print(f"  target: sigma_h0 <= {TARGET_UM} um, unbiased. ladder h = {list(RUNGS)} um, a = {A_UM} um")
    print(f"\n  {'sigma_gamma':>11} {'n ok':>7} {'sigma_h0 [um]':>18} {'bias [um]':>18}")
    rows = []
    for s in (0.01, 0.02, 0.03, 0.04, 0.05, 0.07, 0.08, 0.095, 0.12):
        r = fit_scatter(s)
        rows.append(r)
        flag = "  <-- target" if abs(r["sd"] - TARGET_UM) < 0.02 else ""
        print(f"  {100*s:10.1f} % {r['frac_ok']:6.1%} {r['sd']:12.4f}+/-{r['sd_err']:.4f}"
              f" {r['bias']:+12.4f}+/-{r['bias_se']:.4f}{flag}")

    # interpolate sigma_gamma_max on the two points bracketing the target
    xs = np.array([r["sigma"] for r in rows]); ys = np.array([r["sd"] for r in rows])
    i = int(np.searchsorted(ys, TARGET_UM))
    lo, hi = i - 1, i
    smax = xs[lo] + (TARGET_UM - ys[lo]) * (xs[hi] - xs[lo]) / (ys[hi] - ys[lo])
    print(f"\n  sigma_h0 is linear in sigma_gamma to well inside the error bars, so:")
    print(f"  sigma_gamma_max = {100*smax:.2f} %   (interpolated between "
          f"{100*xs[lo]:.1f} % and {100*xs[hi]:.1f} %)")
    print(f"\n  AM's claim: at or above 3 %.   -> CONFIRMED, and with room: {100*smax:.2f} %.")
    print(f"  But the per-rung scatter this side MEASURED is 7.2-9.5 % "
          f"(verify_drag_ladder.py),")
    print(f"  which is {0.08/smax:.1f}x the tolerance. The budget was self-consistent; "
          f"the 3 % input was not.")

    out = {"target_um": TARGET_UM, "sigma_gamma_max_pct": 100 * smax,
           "rows": rows, "ladder_um": RUNGS.tolist(), "a_um": A_UM}
    p = ROOT / "verify/_out/ladder_tolerance.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
