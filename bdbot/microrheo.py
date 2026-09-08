"""GSER — probe MSD to complex modulus. Mason's algebraic form. LLM 0 lines.

Distilled source: [`knowledge/source/papers/1995-mason-weitz-gser-microrheology.md`
](../knowledge/source/papers/1995-mason-weitz-gser-microrheology.md). Every
threshold, trap and verification value below is taken from that file, not from
memory.

## Why this module is in `bdbot/` and not where the distillation says

The distillation's header declares `used_by: [simbot/analysis/microrheo.py, …]`
and the file **did not exist** — measured 2026-09-02, `simbot/analysis/` held only
`structure.py` and `trap.py`. It cannot exist there either: `simbot` imports
`bdbot` in four modules (`units`, `guards`, `nondim`, `estimators`), so a `bdbot`
case calling an estimator in `simbot/analysis/` would be a cycle. All 8 real cases
run through `bdbot`. So the estimator lives here, beside `lockin.py` and
`stats.py`, which are the other two analysis modules.

## The formula

    |G*(w)|  =  kT / ( pi a <dr^2(1/w)> Gamma(1 + alpha(w)) )
    G'(w)    =  |G*(w)| cos(pi alpha / 2)
    G''(w)   =  |G*(w)| sin(pi alpha / 2)
    alpha(w) =  d ln <dr^2(t)> / d ln t   evaluated at t = 1/w

`a` is the probe **RADIUS** and `<dr^2>` is the **3-component** MSD. Both are the
two traps below.

## The three traps, from the distillation, each guarded here

  1. **`a` is the radius.** Passing the diameter makes every modulus exactly 2x
     wrong -- a pure factor, no shape change, so in a viscoelastic fluid it looks
     like a real modulus of the wrong size. `radius=` is a required keyword and
     `assert_radius_not_diameter()` exists to make the check explicit in a case.
  2. **The MSD must be the 3-component sum.** A 2-component measurement needs
     `<dr^2>_3D = (3/2) <dr^2>_2D`. The distillation calls this *"the most common
     error in microrheology."*
  3. ★ **GSER's derivation presumes a 3D spherical probe, so it may not be applied
     to a 2D system at all** -- not even after the 3/2 conversion. 7 of the 8 cases
     in this repository are 2D, so this refusal is the common path, and it is a
     `raise`, not a warning.

Two more, guarded the same way: `Gamma(1+alpha)` omitted costs ~10 % near
`alpha = 0.5`, and `alpha` outside `[0, 1]` is physically invalid at that
frequency and is masked out rather than reported.

## Verification values this module is tested against

Free probe in a Newtonian solvent, reduced units (`sigma = kT = gamma = 1`, so
`a = 1/2`, `D* = 1`, `eta* = 1/(3 pi)`):

    G'*(w*)  = 0                    exactly
    G''*(w*) = w* / (3 pi)          = 0.106103 w*

`tests/test_microrheo.py` recovers both from an analytic MSD before any simulation
exists. Per [`deterministic-core`](../.claude/rules/deterministic-core.md): an
estimator is tested on input whose answer is already known, because otherwise a
disagreement has two candidate causes instead of one.
"""
from __future__ import annotations

import math

import numpy as np

SCHEMA = "bdbot.microrheo/0.1"

#: GSER's derivation presumes a 3D spherical probe (distillation, trap 2).
GSER_PROBE_DIM = 3

#: Reduced-unit Newtonian slope, `eta* = 1/(3 pi)`. The regression target.
ETA_STAR_NEWTONIAN = 1.0 / (3.0 * math.pi)      # 0.10610329539459689

#: Decades trimmed from each end of the MSD lag range. The distillation:
#: "the trust band is only the interval with one decade cut from each end."
TRIM_DECADES = 1.0

#: Half-width of the local log-log window used for `alpha`. 2 -> 5 points.
ALPHA_HALF_WINDOW = 2

#: `alpha` outside this is superdiffusive or negative-slope -> masked out.
ALPHA_RANGE = (0.0, 1.0)

#: Boundary tolerance on `ALPHA_RANGE`. **Not cosmetic.** A purely viscous fluid
#: has `alpha = 1` exactly, i.e. it sits ON the upper limit, and `np.polyfit`
#: returns `1.0000000000000002` for about half the points of an exact power law.
#: With a strict `<=` the Newtonian reference case masks out 116 of 240 points and
#: therefore rejects itself -- measured 2026-09-02, caught by
#: `tests/test_microrheo.py::test_newtonian_reduced_recovers_eta_star`.
#: Same reason `checks.CMP_RTOL` exists: a threshold a real case sits exactly on
#: needs a tolerance, or the reference answer fails the check built to confirm it.
ALPHA_TOL = 1e-9

#: ★ How to handle `alpha > 1`. **This choice is not cosmetic and it has a
#: measured cost.** The distillation says discard those frequencies as
#: unphysical, which is right for a viscoelastic medium where the true alpha is
#: genuinely below 1. It is WRONG in the Newtonian limit, where the true alpha is
#: exactly 1 and therefore sits ON the boundary: noise scatters symmetrically
#: about it and the mask keeps only the low side.
#:
#: Measured 2026-09-02 on a purely viscous fluid, 100 probes x 60,000 steps:
#:     in the trust band, before masking : 5991 points, mean alpha 0.99626
#:                                          2137 above 1, 3854 below
#:     after masking alpha > 1           : 3854 points, mean alpha 0.99140
#:                                          0 above 1  <- one-sided
#: The surviving alpha is biased low by 0.0086, and via
#: `G'/G'' = cot(pi alpha/2) ~= (pi/2)(1 - alpha)` that predicts a SPURIOUS
#: `G'/G'' = 0.01351` where the true value is 0. Measured: 0.01400. So the mask
#: manufactured a 1.4 % storage modulus in a fluid that has none, and it did NOT
#: shrink with run length (0.0152 / 0.0093 / 0.0103 at 30k / 60k / 120k steps) --
#: it is a bias, not noise. It was caught by a sealed prediction of
#: `G'/G'' <= 0.01`, which failed at 0.01308.
#:
#: `clip` sets alpha>1 to 1, so those frequencies contribute `cos(pi/2) = 0`,
#: which is the correct physical answer, and the estimate stays unbiased.
ALPHA_POLICY_NOTE = "mask: discard alpha>1 (viscoelastic) | clip: alpha->1 (near-Newtonian)"


def msd_2d_to_3d(msd_2d):
    """`<dr^2>_3D = (3/2) <dr^2>_2D`. The conversion, exposed so that a case has to
    name it rather than fold a 1.5 into an expression.

    ⚠️ This makes the *number* right and does not make GSER *applicable* -- see
    trap 3 in the module docstring. `msd_to_gstar` still refuses a 2D system.
    """
    return 1.5 * np.asarray(msd_2d, dtype=float)


def assert_radius_not_diameter(radius, sigma, rtol=1e-9):
    """Raise if `radius` looks like it is actually the diameter.

    Trap 1 costs a factor of exactly 2 in every modulus with no change of shape,
    so it cannot be spotted in a plot. A case calls this once with its own
    `sigma` and the mistake becomes impossible rather than unlikely.
    """
    radius, sigma = float(radius), float(sigma)
    if radius <= 0:
        raise ValueError(f"radius must be > 0, got {radius}")
    if math.isclose(radius, sigma, rel_tol=rtol):
        raise ValueError(
            f"radius={radius} equals sigma={sigma}: that is the DIAMETER. "
            f"GSER's `a` is the radius; passing the diameter makes every modulus "
            f"exactly 2x wrong with no change of shape (distillation, trap 1). "
            f"Pass sigma/2 = {sigma / 2}.")
    return float(radius)


def alpha_loglog(t, msd, half_window: int = ALPHA_HALF_WINDOW):
    """Local logarithmic slope `d ln msd / d ln t` by a local **quadratic** fit.

    Not a finite difference: the distillation requires a local polynomial fit in
    log space because `alpha` is differentiated from noisy data, and a two-point
    slope on log-spaced lags is dominated by the spacing.

    Returns an array the same length as `t`; the outermost `half_window` points
    are fitted from a one-sided window and are trimmed away later anyway.
    """
    t = np.asarray(t, dtype=float)
    msd = np.asarray(msd, dtype=float)
    if t.ndim != 1 or t.shape != msd.shape:
        raise ValueError(f"t and msd must be 1-D and the same length, "
                         f"got {t.shape} and {msd.shape}")
    n = t.size
    if n < 2 * half_window + 3:
        raise ValueError(f"need at least {2 * half_window + 3} MSD points for a "
                         f"local quadratic fit, got {n}")
    if np.any(t <= 0):
        raise ValueError("t must be > 0 (log-log fit)")
    if np.any(msd <= 0):
        raise ValueError("msd must be > 0 (log-log fit); a zero or negative lag "
                         "value means the MSD estimator, not this function, is wrong")
    if np.any(np.diff(t) <= 0):
        raise ValueError("t must be strictly increasing")

    x, y = np.log(t), np.log(msd)
    out = np.empty(n)
    for i in range(n):
        lo, hi = max(0, i - half_window), min(n, i + half_window + 1)
        if hi - lo < 3:                       # cannot fit a quadratic
            lo, hi = max(0, min(lo, n - 3)), min(n, max(hi, lo + 3))
        dx = x[lo:hi] - x[i]
        # y = c0 + c1 dx + c2 dx^2  ->  d y / d x at dx=0 is c1
        c2, c1, c0 = np.polyfit(dx, y[lo:hi], 2)
        out[i] = c1
    return out


def _trim_mask(t, trim_decades: float):
    """Keep only lags at least `trim_decades` inside each end (distillation)."""
    t = np.asarray(t, dtype=float)
    lo = t[0] * 10.0 ** trim_decades
    hi = t[-1] / 10.0 ** trim_decades
    return (t >= lo) & (t <= hi)


def msd_to_gstar(t, msd, *, radius, kT=1.0, msd_dim=GSER_PROBE_DIM,
                 system_dim=GSER_PROBE_DIM, allow_2d=False,
                 trim_decades=TRIM_DECADES, half_window=ALPHA_HALF_WINDOW,
                 apply_gamma=True, alpha_policy="mask"):
    """GSER: probe MSD -> `G'(w)`, `G''(w)`. Mason's algebraic form.

    Parameters
    ----------
    t, msd
        Lag times and the **3-component** MSD, same length, both > 0, `t`
        strictly increasing.
    radius
        The probe **radius** `a`. Keyword-only on purpose (trap 1).
    kT
        Thermal energy in the same units as `radius`^3 * modulus.
    msd_dim
        How many components the supplied `msd` sums. If 2, it is converted with
        `msd_2d_to_3d` and that is recorded in the result.
    system_dim
        The dimensionality of the *simulation*. 2 raises unless `allow_2d`
        (trap 3): GSER's derivation presumes a 3D spherical probe.
    apply_gamma
        Include `Gamma(1 + alpha)`. Only set False to reproduce the ~10 % error
        the distillation attributes to omitting it -- `tests/` does exactly that.

    Returns
    -------
    dict with `omega`, `g_prime`, `g_double_prime`, `g_abs`, `alpha`, `valid`,
    and a `notes` list recording every conversion and refusal that applied.
    """
    notes = []
    radius = float(radius)
    if radius <= 0:
        raise ValueError(f"radius must be > 0, got {radius}")

    if int(system_dim) != GSER_PROBE_DIM and not allow_2d:
        raise ValueError(
            f"system_dim={system_dim}: GSER's derivation presumes a 3D spherical "
            f"probe, so it may not be applied to a {system_dim}D system -- not "
            f"even after the 3/2 MSD conversion (distillation, trap 2/3). "
            f"7 of the 8 cases in this repository are 2D, so this is the common "
            f"path, not an edge case. Pass allow_2d=True only with a written "
            f"justification in the run's goal or spec.")
    if int(system_dim) != GSER_PROBE_DIM:
        notes.append(f"system_dim={system_dim} accepted via allow_2d=True -- "
                     f"GSER is being applied outside its derivation")

    msd = np.asarray(msd, dtype=float)
    if int(msd_dim) == 2:
        msd = msd_2d_to_3d(msd)
        notes.append("msd converted 2-component -> 3-component with a factor 3/2")
    elif int(msd_dim) != 3:
        raise ValueError(f"msd_dim must be 2 or 3, got {msd_dim}")

    t = np.asarray(t, dtype=float)
    alpha = alpha_loglog(t, msd, half_window=half_window)

    gam = np.array([math.gamma(1.0 + a) if apply_gamma else 1.0 for a in alpha])
    if not apply_gamma:
        notes.append("Gamma(1+alpha) NOT applied -- expect ~10 % error near "
                     "alpha=0.5 (distillation, trap 3)")

    with np.errstate(divide="ignore", invalid="ignore"):
        g_abs = kT / (math.pi * radius * msd * gam)
        g_p = g_abs * np.cos(math.pi * alpha / 2.0)
        g_pp = g_abs * np.sin(math.pi * alpha / 2.0)

    omega = 1.0 / t

    lo, hi = ALPHA_RANGE
    lo_t, hi_t = lo - ALPHA_TOL, hi + ALPHA_TOL      # see ALPHA_TOL: alpha=1 is a
    n_over = int((alpha > hi_t).sum())               # real case sitting on the limit

    if alpha_policy == "clip":
        # ★ See ALPHA_POLICY_NOTE. Clip alpha>1 to 1 instead of discarding it:
        #   cos(pi/2)=0 so those frequencies contribute G'=0, which IS the correct
        #   physical answer, and the estimate stays unbiased.
        alpha = np.minimum(alpha, hi)
        gam = np.array([math.gamma(1.0 + a) if apply_gamma else 1.0 for a in alpha])
        with np.errstate(divide="ignore", invalid="ignore"):
            g_abs = kT / (math.pi * radius * msd * gam)
            g_p = g_abs * np.cos(math.pi * alpha / 2.0)
            g_pp = g_abs * np.sin(math.pi * alpha / 2.0)
        notes.append(f"alpha_policy='clip': {n_over} point(s) with alpha>1 clipped "
                     f"to 1 rather than discarded (avoids the one-sided bias)")
    elif alpha_policy != "mask":
        raise ValueError(f"alpha_policy must be 'mask' or 'clip', got {alpha_policy!r}")

    keep_alpha = ((alpha >= lo_t) & (alpha <= hi_t) if alpha_policy == "mask"
                  else (alpha >= lo_t))
    valid = _trim_mask(t, trim_decades) & keep_alpha & np.isfinite(g_abs)
    n_trim = int((~_trim_mask(t, trim_decades)).sum())
    n_alpha = int((~keep_alpha).sum())
    notes.append(f"{n_trim} point(s) outside the {trim_decades:g}-decade trust band")
    if n_alpha:
        notes.append(f"{n_alpha} point(s) masked for alpha outside [{lo}, {hi}] "
                     f"(superdiffusive or negative slope -> unphysical there)")

    # frequency runs opposite to lag time; report ascending in omega
    order = np.argsort(omega)
    return {
        "schema": SCHEMA,
        "omega": omega[order],
        "g_prime": g_p[order],
        "g_double_prime": g_pp[order],
        "g_abs": g_abs[order],
        "alpha": alpha[order],
        "valid": valid[order],
        "radius": radius,
        "kT": float(kT),
        "n_valid": int(valid.sum()),
        "band_decades": (math.log10(t[-1] / t[0]) - 2.0 * trim_decades
                         if t[0] > 0 else float("nan")),
        "notes": notes,
    }


def newtonian_msd(t, D, dim=GSER_PROBE_DIM):
    """`<dr^2> = 2 dim D t`. The analytic input the estimator is tested on."""
    return 2.0 * int(dim) * float(D) * np.asarray(t, dtype=float)


def newtonian_target(omega, eta):
    """What GSER must return for a Newtonian solvent: `(G', G'') = (0, eta w)`."""
    omega = np.asarray(omega, dtype=float)
    return np.zeros_like(omega), float(eta) * omega


def powerlaw_msd(t, amp, alpha):
    """`<dr^2> = amp t^alpha`. A single-exponent input, so `alpha` is exact and
    `G'/G'' = cot(pi alpha / 2)` is a closed-form check independent of the
    Newtonian one."""
    return float(amp) * np.asarray(t, dtype=float) ** float(alpha)


def powerlaw_target(t, amp, alpha, *, radius, kT=1.0):
    """Closed form for `powerlaw_msd`, so a test does not re-derive the estimator."""
    t = np.asarray(t, dtype=float)
    g_abs = kT / (math.pi * float(radius) * powerlaw_msd(t, amp, alpha)
                  * math.gamma(1.0 + alpha))
    return (g_abs * math.cos(math.pi * alpha / 2.0),
            g_abs * math.sin(math.pi * alpha / 2.0))


__all__ = ["SCHEMA", "GSER_PROBE_DIM", "ETA_STAR_NEWTONIAN", "TRIM_DECADES",
           "ALPHA_HALF_WINDOW", "ALPHA_RANGE", "ALPHA_TOL", "ALPHA_POLICY_NOTE",
           "msd_2d_to_3d",
           "assert_radius_not_diameter", "alpha_loglog", "msd_to_gstar",
           "newtonian_msd", "newtonian_target", "powerlaw_msd", "powerlaw_target"]
