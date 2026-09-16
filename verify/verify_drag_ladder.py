"""The drag-slope ladder: what IS the per-rung precision the plan's primary route rests on?

AM's D-5: the run does not get `gamma` from `f_c`. It gets it from
`gamma = alpha_equipartition x slope`, where `slope = d<x_eq>/dv` from the drag
segments and `alpha = kT/var(x)` from the held segment. The ~3 % per-rung figure
the whole error budget stands on is a toy-model assumption (r1's 400-realisation
numpy estimate at 3 % Gaussian per point) and **neither side has measured it**.

This measures it, and then propagates it into the six-rung Faxen fit for `h0`.

Exact-OU again, for the reason r4 and r5 used it: the sampler is unbiased at any
step size, so the detector bias and the finite-window statistics are not mixed
with Euler-Maruyama error. The companion HOOMD run is the independent check that
the OU reduction is the right model of the BD system, not a substitute for it.

⚠ The wall is NOT simulated and cannot be: `k*` carries no drag, so all six rungs
are the same dimensionless run and any wall factor reported would be the one put
in (BD r4/r5, AM D-6, independently). `gamma(h)` here is a SYNTHETIC ladder built
from AM's declared Faxen factor. What is measured is the FIT PROCEDURE.

    $PY verify/verify_drag_ladder.py
"""
from __future__ import annotations

import json, math, sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT)); sys.path.insert(0, str(ROOT / "verify"))
from detector_model import apply_detector, blur_variance_factor, ou_exact   # noqa: E402

from bdbot import Q                                                         # noqa: E402
from bdbot import constants as C                                            # noqa: E402

PI = math.pi
# ---- SI, from the approved manifest (verify/_out/params_trap_drag_ladder.txt) ----
d = Q(4.95, "um"); a = d / 2
T = Q(293.15, "K"); eta = Q(1.0016e-3, "Pa*s"); k_t = Q(3.5054, "pN/um")
kT = (Q(C.K_B, "J/K") * T).to("pN*um")
G0 = (6 * PI * eta * a).to("pN*s/um")           # BULK -- the manifest's derived.gamma
L_K = float(((kT / k_t) ** 0.5).to("nm").m)     # 33.98 nm, h-independent
TAU_K0 = float((G0 / k_t).to("s").m)
A_UM = float(a.to("um").m)

T_EXP_S = 0.45e-3
FS = 520.0
EPS_NM = 10.0
T_SEG_S = 5.0                                   # AM sequence 9b/9c: 5 s per segment
X_EQ_FRAC = 0.15                                # v chosen so x_eq = 0.15a at each rung
V_FRACS = (0.25, 0.50, 0.75, 1.00)
RUNGS = (3.0, 3.5, 4.5, 6.0, 8.0, 10.0)
H0_TRUE = 0.0                                   # the ladder is written in absolute h here;
                                                # the fit's unknown offset is injected below
N_BEADS = 64
N_LADDERS = 48                                  # independent six-rung ladders -> sigma_h0
N_EXP = 12                                      # sub-steps per exposure


def faxen(h_um: float) -> float:
    """gamma(h)/gamma_0, parallel, first order. AM's declared form. NOT derived here."""
    return 1.0 / (1.0 - 9.0 * A_UM / (16.0 * h_um))


def one_rung(h_um, rng, n_beads=N_BEADS, camera=True):
    """Return (slope, alpha, gamma) per bead at one rung. Everything in um / s / pN."""
    fax = faxen(h_um)
    tau_k = TAU_K0 * fax                                   # s
    sig_um = L_K * 1e-3                                    # thermal amplitude, um
    u = T_EXP_S / tau_k
    frame_tau = (1.0 / FS) / tau_k
    dt = u / N_EXP                                         # in units of tau_k
    n_frame = max(N_EXP, int(round(frame_tau / dt)))
    frame_eff = n_frame * dt
    n_frames = int(round((T_SEG_S / tau_k) / frame_eff))
    n = n_frames * n_frame
    eps = (EPS_NM * 1e-3) / sig_um if camera else 0.0

    def segment():
        x = ou_exact(n, dt, 1.0, rng, n_beads)             # units of sigma
        if camera:
            fr, _ = apply_detector(x, dt, u, frame_eff, eps, rng)
        else:
            fr, _ = apply_detector(x, dt, 2 * dt, frame_eff, 0.0, rng)
        del x
        return fr * sig_um                                 # -> um

    # --- drag segments: <x> at four velocities -------------------------------
    gamma_h = float(G0.m) * fax                            # pN*s/um
    x_eq_max = X_EQ_FRAC * A_UM                            # um
    v_max = x_eq_max * float(k_t.m) / gamma_h              # um/s
    vs = np.array([f * v_max for f in V_FRACS])
    xbar = np.empty((len(vs), n_beads))
    for j, v in enumerate(vs):
        xbar[j] = segment().mean(axis=0) + gamma_h * v / float(k_t.m)
    # least squares slope of <x> vs v, per bead
    vc = vs - vs.mean()
    slope = (vc[:, None] * (xbar - xbar.mean(axis=0))).sum(axis=0) / (vc ** 2).sum()

    # --- held segment: alpha = kT/var(x) -------------------------------------
    held = segment()
    alpha = float(kT.m) / held.var(axis=0, ddof=1)         # pN/um

    return slope, alpha, alpha * slope, gamma_h, u


def main():
    rng = np.random.default_rng(20260916)
    out = {"rungs": {}, "inputs": {"l_k_nm": L_K, "tau_k0_ms": TAU_K0 * 1e3,
                                   "gamma_0_pN_s_um": float(G0.m), "k_t_pN_um": float(k_t.m)}}

    print("=" * 96)
    print("1 - per-rung precision of the drag-slope route   (camera ON: exposure + epsilon)")
    print("=" * 96)
    print(f"  {'h [um]':>7} {'gamma(h)':>9} {'u=t_e/tau':>10} {'sigma_slope':>12}"
          f" {'sigma_alpha':>12} {'sigma_gamma':>12} {'bias_gamma':>11}")
    per_rung = {}
    for h in RUNGS:
        sl, al, gm, g_true, u = one_rung(h, rng)
        s_sl = 100 * sl.std(ddof=1) / sl.mean()
        s_al = 100 * al.std(ddof=1) / al.mean()
        s_gm = 100 * gm.std(ddof=1) / gm.mean()
        b_gm = 100 * (gm.mean() / g_true - 1)
        per_rung[h] = dict(sigma_slope_pct=s_sl, sigma_alpha_pct=s_al,
                           sigma_gamma_pct=s_gm, bias_gamma_pct=b_gm,
                           gamma_true=g_true, u=u)
        print(f"  {h:7.1f} {g_true:9.5f} {u:10.5f} {s_sl:11.2f} % {s_al:11.2f} %"
              f" {s_gm:11.2f} % {b_gm:+10.2f} %")
    out["rungs"] = {str(k): v for k, v in per_rung.items()}
    print(f"\n  AM's assumption, on which the whole error budget rests: 3 % per rung.")

    # --- camera off, to separate statistics from the detector -----------------
    print("\n" + "=" * 96)
    print("2 - the same, camera OFF -- which part of the error is the detector?")
    print("=" * 96)
    rng2 = np.random.default_rng(7788)
    print(f"  {'h [um]':>7} {'sigma_gamma':>12} {'bias_gamma':>11}")
    for h in (3.0, 8.0, 10.0):
        sl, al, gm, g_true, _ = one_rung(h, rng2, camera=False)
        out["rungs"][str(h)]["sigma_gamma_pct_nocam"] = 100 * gm.std(ddof=1) / gm.mean()
        out["rungs"][str(h)]["bias_gamma_pct_nocam"] = 100 * (gm.mean() / g_true - 1)
        print(f"  {h:7.1f} {out['rungs'][str(h)]['sigma_gamma_pct_nocam']:11.2f} %"
              f" {out['rungs'][str(h)]['bias_gamma_pct_nocam']:+10.2f} %")

    p = ROOT / "verify/_out/drag_ladder.json"
    p.write_text(json.dumps(out, indent=1))
    print(f"\nwrote {p.relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    sys.exit(main())


# =============================================================================
# 3 - the six-rung Faxen fit. What the per-rung error above does to h0.
# =============================================================================
def ladder_fit(h0_true_um=1.0, n_ladders=N_LADDERS, camera=True, seed=4242):
    """Recover (gamma_0, h0) from six rungs, with the MEASURED per-rung error.

    The piezo reads `dz`; the absolute height is `h = dz + h0` with `h0` unknown.
    So the six commanded `dz` are `RUNGS - h0_true`, and the fit gets `dz` and the
    measured `gamma`, never `h`.
    """
    from scipy import optimize
    rng = np.random.default_rng(seed)
    dz = np.array(RUNGS) - h0_true_um
    g0_true = float(G0.m)

    # per-bead gamma at each rung, drawn afresh for each ladder
    samples = {}
    for h in RUNGS:
        _, _, gm, g_true, _ = one_rung(h, rng, n_beads=n_ladders, camera=camera)
        samples[h] = gm

    def model(dz_, g0_, h0_):
        return g0_ / (1.0 - 9.0 * A_UM / (16.0 * (dz_ + h0_)))

    rec_g0, rec_h0 = [], []
    for i in range(n_ladders):
        y = np.array([samples[h][i] for h in RUNGS])
        try:
            popt, _ = optimize.curve_fit(model, dz, y, p0=[g0_true, h0_true_um],
                                         maxfev=40000)
        except Exception:
            continue
        rec_g0.append(popt[0]); rec_h0.append(popt[1])
    return np.array(rec_g0), np.array(rec_h0), g0_true, h0_true_um


def part3():
    print("\n" + "=" * 96)
    print("3 - the six-rung fit for h0, with the per-rung error MEASURED above")
    print("=" * 96)
    res = {}
    for camera, lab in ((True, "camera ON "), (False, "camera OFF")):
        g0, h0, g0_true, h0_true = ladder_fit(camera=camera)
        n = len(h0)
        se = lambda v: v.std(ddof=1) / math.sqrt(len(v))
        res[lab.strip()] = dict(
            n=n, h0_mean=float(h0.mean()), h0_sd=float(h0.std(ddof=1)),
            h0_bias=float(h0.mean() - h0_true),
            g0_mean=float(g0.mean()), g0_bias_pct=float(100 * (g0.mean() / g0_true - 1)),
            g0_sd_pct=float(100 * g0.std(ddof=1) / g0_true))
        print(f"  {lab}  n={n:3d}   h0 = {h0.mean():+.4f} +/- {se(h0):.4f} um "
              f"(true {h0_true:.2f})   bias {h0.mean()-h0_true:+.4f} um   "
              f"sd {h0.std(ddof=1):.4f} um")
        print(f"              gamma_0 bias {res[lab.strip()]['g0_bias_pct']:+7.2f} %   "
              f"sd {res[lab.strip()]['g0_sd_pct']:6.2f} %")
    print(f"\n  AM's claim (r1, 400-realisation numpy at 3 % per rung): h0 to +/- 0.195 um, "
          f"gamma_bulk to 2.1 %, UNBIASED.")
    p = ROOT / "verify/_out/drag_ladder_fit.json"
    p.write_text(json.dumps(res, indent=1))
    print(f"\nwrote {p.relative_to(ROOT)}")


if __name__ == "__main__":
    part3()
