"""Check the camera layer against answers that are known in closed form.

Three questions, in the order rule 7 asks them -- one element at a time:

  1. does `ou_exact` reproduce var and C(t)?            (the generator)
  2. does exposure deflate var(x) by the exact OU        (blur alone, no noise)
     factor, and is AM's `2 D t_exp / 3` the right
     coefficient for a TRAPPED bead?
  3. what do blur and localisation noise each do to      (the readout)
     `f_c` from a Lorentzian fit?

Synthetic OU, not a BD run, on purpose: the OU sampler is exact at any step
size, so there is no Euler-Maruyama bias to confuse with the detector bias being
measured. What this script does NOT do is verify the physics of the trap -- that
is `runs/trap-2d-5um__a5ef4f45d589`. It asks only whether the detector layer
does what it claims.

Numbers are the trap-stiffness-recovery thread's:
tau = 16.14 ms, sigma = 33.98 nm, t_exp = 0.45 ms, f_s = 520 Hz, epsilon = 10 nm.

    $PY verify/verify_detector_model.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "verify"))

from detector_model import (apply_detector, blur_deflation_exact,           # noqa: E402
                            blur_deflation_free_particle, blur_variance_factor,
                            lorentzian_fc, ou_exact)

# ---- the thread's physical numbers (r1 ask_simulation / r2 ask_experiment) ----
TAU_MS = 16.14          # gamma_corr/k_t, from r2 parsed_back.derived_here
SIGMA_NM = 33.98        # l_k = sqrt(kT/k_t)
T_EXP_MS = 0.45         # r1 instrument_envelope, the run's chosen exposure
FS_HZ = 520.0           # r1 instrument_envelope, hard
EPS_NM = 10.0           # r1 instrument_envelope, assumed pending P8

N_PATHS = 16
N_BATCH = 6             # independent batches -> the error bar (A2)
T_OBS_OVER_TAU = 400


def sem(v):
    v = np.asarray(v, dtype=float)
    return float(v.std(ddof=1) / math.sqrt(len(v)))


def banner(s):
    print(f"\n{'=' * 78}\n{s}\n{'=' * 78}")


# =============================================================================
# 1 - the generator itself
# =============================================================================
def check_generator():
    banner("1 - ou_exact: does the generator have the right var and C(t)?")
    rng = np.random.default_rng(20260915)
    dt_over_tau = 0.05
    n = int(T_OBS_OVER_TAU / dt_over_tau)
    v, tfit = [], []
    for _ in range(N_BATCH):
        x = ou_exact(n, dt_over_tau, 1.0, rng, N_PATHS)
        v.append(float((x ** 2).mean()))
        # C(t) by FFT, unbiased, then a 1/e crossing read off the curve
        xc = x - 0.0
        m = 1 << (2 * len(xc) - 1).bit_length()
        F = np.fft.rfft(xc, n=m, axis=0)
        ac = np.fft.irfft(F * np.conj(F), n=m, axis=0)[:len(xc)].mean(axis=1)
        ac /= np.arange(len(xc), 0, -1)
        ac /= ac[0]
        k = int(np.argmax(ac < math.exp(-1.0)))
        # linear interpolation on the crossing
        frac = (ac[k - 1] - math.exp(-1.0)) / (ac[k - 1] - ac[k])
        tfit.append((k - 1 + frac) * dt_over_tau)
    print(f"  var / sigma^2      = {np.mean(v):.5f} +/- {sem(v):.5f}   (exact 1)")
    print(f"  tau_1/e / tau      = {np.mean(tfit):.5f} +/- {sem(tfit):.5f}   (exact 1)")
    ok = abs(np.mean(v) - 1) < 4 * sem(v) + 1e-3 and abs(np.mean(tfit) - 1) < 0.02
    print(f"  -> {'PASS' if ok else 'FAIL'}")
    return ok


# =============================================================================
# 2 - blur alone. The coefficient question.
# =============================================================================
def check_blur():
    banner("2 - exposure blur: measured deflation of var(x) vs two formulas")
    print(f"  {'u=t_exp/tau':>12}  {'measured':>20}  {'exact 1-2(u-1+e^-u)/u^2':>24}"
          f"  {'AM 2u/3':>10}")
    rng = np.random.default_rng(4711)
    rows, ok = [], True
    for u in (0.027881, 0.1, 0.3, 1.0, 3.0):
        # dt fine enough that the boxcar is a boxcar; frame >> tau so frames are
        # independent and the paired ratio is not correlated across frames.
        n_exp = 32
        dt = u / n_exp
        n_frame = max(n_exp, int(round(5.0 / dt)))          # frame = 5 tau
        n = n_frame * 400
        r = []
        for _ in range(N_BATCH):
            x = ou_exact(n, dt, 1.0, rng, N_PATHS)
            blurred, _ = apply_detector(x, dt, u, n_frame * dt, 0.0, rng)
            inst = x[: (len(x) // n_frame) * n_frame].reshape(-1, n_frame, x.shape[1])[:, 0, :]
            r.append(float((blurred ** 2).mean() / (inst ** 2).mean()))
        meas = 1.0 - np.mean(r)
        err = sem(r)
        ex = blur_deflation_exact(u)
        am = blur_deflation_free_particle(u)
        n_sig_ex = abs(meas - ex) / err if err else float("inf")
        n_sig_am = abs(meas - am) / err if err else float("inf")
        rows.append((u, meas, err, ex, am, n_sig_ex, n_sig_am))
        print(f"  {u:12.6f}  {meas*100:11.4f} +/- {err*100:5.4f} %"
              f"  {ex*100:21.4f} %  {am*100:8.4f} %"
              f"     [exact {n_sig_ex:5.1f} sigma, AM {n_sig_am:6.1f} sigma]")
        if n_sig_ex > 4:
            ok = False
    print(f"\n  role: exact formula = implementation_check (it is what the code "
          f"claims to do)\n        AM's 2u/3         = hypothesis (their assumption, "
          f"not imposed here)")
    print(f"  -> {'PASS' if ok else 'FAIL'} on the implementation_check")
    return ok, rows


# =============================================================================
# 3 - what the detector does to f_c
# =============================================================================
def check_fc():
    banner("3 - f_c from a Lorentzian fit: blur down, localisation noise up")
    tau_s = TAU_MS * 1e-3
    fc_true = 1.0 / (2 * math.pi * tau_s)
    u = (T_EXP_MS * 1e-3) / tau_s
    eps = EPS_NM / SIGMA_NM
    print(f"  tau = {TAU_MS} ms -> f_c(true) = {fc_true:.4f} Hz")
    print(f"  u = t_exp/tau = {u:.6f} . epsilon/sigma = {eps:.5f} . f_s = {FS_HZ} Hz")

    n_exp = 24
    dt = u / n_exp                                   # in units of tau
    frame_tau = (1.0 / FS_HZ) / tau_s
    n_frame = int(round(frame_tau / dt))
    frame_eff = n_frame * dt
    print(f"  sampling: dt/tau = {dt:.6g} ({n_exp} steps per exposure), "
          f"frame = {n_frame} steps -> f_s(effective) = "
          f"{1.0/(frame_eff*tau_s):.2f} Hz")

    n = n_frame * 20000
    rng = np.random.default_rng(90210)
    out = {k: [] for k in ("ideal", "blur", "noise", "both")}
    for _ in range(N_BATCH):
        x = ou_exact(n, dt, 1.0, rng, N_PATHS)
        fs_frames = 1.0 / (frame_eff * tau_s)
        ideal, _ = apply_detector(x, dt, 2 * dt, frame_eff, 0.0, rng)   # ~instant
        blur, _ = apply_detector(x, dt, u, frame_eff, 0.0, rng)
        noise, _ = apply_detector(x, dt, 2 * dt, frame_eff, eps, rng)
        both, _ = apply_detector(x, dt, u, frame_eff, eps, rng)
        for k, series in (("ideal", ideal), ("blur", blur),
                          ("noise", noise), ("both", both)):
            out[k].append(lorentzian_fc(series, fs_frames)[0])

    print(f"\n  {'series':>22}  {'f_c [Hz]':>18}  {'vs true':>10}")
    res = {}
    for k, label in (("ideal", "no detector"), ("blur", "exposure only"),
                     ("noise", "epsilon only"), ("both", "exposure + epsilon")):
        m, e = float(np.mean(out[k])), sem(out[k])
        res[k] = (m, e)
        print(f"  {label:>22}  {m:10.4f} +/- {e:5.4f}  {100*(m/fc_true-1):+9.2f} %")
    print(f"\n  variance check (what equipartition alpha is read from):")
    return res, fc_true


def main():
    ok1 = check_generator()
    ok2, _ = check_blur()
    check_fc()
    banner("verdict")
    print(f"  generator   {'PASS' if ok1 else 'FAIL'}")
    print(f"  blur layer  {'PASS' if ok2 else 'FAIL'}")
    return 0 if (ok1 and ok2) else 1


if __name__ == "__main__":
    sys.exit(main())
