"""GSER estimator, tested on analytic input whose answer is already known.

No HOOMD, no run, no trajectory. Per
[`deterministic-core`](../.claude/rules/deterministic-core.md): *"Test a new
analysis estimator against synthetic data whose answer is known. That is not
physics verification — it asks 'does the code answer correctly on an input whose
answer it should know'."* If the estimator cannot recover `0.106103` from a
perfect MSD, no simulation will tell us anything.

Two independent analytic inputs, so a single algebra error cannot pass both:

  Newtonian   <dr^2> = 2 dim D t   ->  G' = 0,  G'' = eta w
              reduced:                 G''* = w*/(3 pi) = 0.106103 w*
  Power law   <dr^2> = A t^alpha   ->  G'/G'' = cot(pi alpha / 2)
              at alpha = 1/2 that is exactly 1, i.e. G' = G''

Then the three traps the distillation names, each asserted to be **real** — the
wrong call must actually produce the wrong number, or the guard is decoration.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from bdbot import microrheo as MR


def log_lags(t_lo=1e-3, t_hi=1e3, n=240):
    """Log-spaced lags spanning 6 decades, so 4 survive the 1-decade trim."""
    return np.logspace(math.log10(t_lo), math.log10(t_hi), n)


# ── 1. the Newtonian limit, in reduced units ────────────────────────────────

def test_newtonian_reduced_recovers_eta_star():
    """★ The headline regression value from the distillation:
    `G''* = w*/(3 pi) = 0.106103 w*` and `G'* = 0`, for a free probe in a
    Newtonian solvent with `sigma = kT = gamma = 1` so `a = 1/2`, `D* = 1`.
    """
    t = log_lags()
    msd = MR.newtonian_msd(t, D=1.0, dim=3)          # reduced: D* = 1
    r = MR.msd_to_gstar(t, msd, radius=0.5, kT=1.0)  # a = sigma/2

    v = r["valid"]
    assert r["n_valid"] > 100, r["notes"]

    # alpha must be 1 to machine-ish precision on an exact power law
    assert np.allclose(r["alpha"][v], 1.0, atol=1e-9)

    # G'' / omega must equal eta* = 1/(3 pi)
    slope = r["g_double_prime"][v] / r["omega"][v]
    assert np.allclose(slope, MR.ETA_STAR_NEWTONIAN, rtol=1e-9), \
        f"recovered eta* = {slope[0]:.9f}, expected {MR.ETA_STAR_NEWTONIAN:.9f}"
    assert abs(slope[0] - 0.106103) < 1e-6      # the digits written in the doc

    # G' must be zero, not merely small: cos(pi/2) at alpha=1
    assert np.allclose(r["g_prime"][v], 0.0, atol=1e-12)


def test_newtonian_dimensional_recovers_eta():
    """Same check with real units, so the reduced-unit convention is not itself
    carrying the result. Water at 300 K, 1 um probe."""
    eta = 0.851e-3          # Pa*s, the Welty anchor the 8 cases ran
    kT = 1.380649e-23 * 300.0
    d = 1.0e-6
    gamma = 3.0 * math.pi * eta * d          # materials.sphere_drag: 3 pi eta d
    D = kT / gamma

    t = log_lags(1e-4, 1e2, 240)
    msd = MR.newtonian_msd(t, D=D, dim=3)
    r = MR.msd_to_gstar(t, msd, radius=d / 2, kT=kT)

    v = r["valid"]
    got = r["g_double_prime"][v] / r["omega"][v]
    assert np.allclose(got, eta, rtol=1e-9), f"recovered eta={got[0]:.6e}, want {eta:.6e}"
    _, want_pp = MR.newtonian_target(r["omega"][v], eta)
    assert np.allclose(r["g_double_prime"][v], want_pp, rtol=1e-9)


# ── 2. a second, independent analytic input ─────────────────────────────────

@pytest.mark.parametrize("alpha", [0.25, 0.5, 0.75])
def test_powerlaw_ratio_is_cot_half_pi_alpha(alpha):
    """`G'/G'' = cot(pi alpha / 2)` for a single-exponent MSD. Independent of the
    Newtonian check: it tests the phase split, which `alpha = 1` cannot (there
    `G'` is identically zero)."""
    t = log_lags()
    msd = MR.powerlaw_msd(t, amp=1.0, alpha=alpha)
    r = MR.msd_to_gstar(t, msd, radius=0.5, kT=1.0)
    v = r["valid"]

    assert np.allclose(r["alpha"][v], alpha, atol=1e-9)
    ratio = r["g_prime"][v] / r["g_double_prime"][v]
    assert np.allclose(ratio, 1.0 / math.tan(math.pi * alpha / 2.0), rtol=1e-9)

    want_p, want_pp = MR.powerlaw_target(t, 1.0, alpha, radius=0.5, kT=1.0)
    assert np.allclose(r["g_prime"][v], want_p[np.argsort(1.0 / t)][v], rtol=1e-9)


def test_alpha_half_gives_equal_moduli():
    """The cleanest single number in the power-law family: at `alpha = 1/2`,
    `cot(pi/4) = 1`, so `G' = G''` exactly."""
    t = log_lags()
    r = MR.msd_to_gstar(t, MR.powerlaw_msd(t, 1.0, 0.5), radius=0.5)
    v = r["valid"]
    assert np.allclose(r["g_prime"][v], r["g_double_prime"][v], rtol=1e-9)


# ── 3. ★ the three traps, each asserted to be REAL ──────────────────────────

def test_diameter_instead_of_radius_is_exactly_a_factor_two():
    """Trap 1. The wrong call must actually produce a wrong number, or the guard
    is decoration. And it must be a pure factor with no change of shape — which
    is precisely why it cannot be seen in a plot."""
    t = log_lags()
    msd = MR.newtonian_msd(t, D=1.0, dim=3)
    good = MR.msd_to_gstar(t, msd, radius=0.5)      # radius
    bad = MR.msd_to_gstar(t, msd, radius=1.0)       # diameter -- the mistake
    v = good["valid"]
    ratio = good["g_double_prime"][v] / bad["g_double_prime"][v]
    assert np.allclose(ratio, 2.0, rtol=1e-12), "trap 1 is not a factor of 2 here"
    # shape identical -> invisible in a log-log plot
    assert np.allclose(good["alpha"][v], bad["alpha"][v], atol=1e-12)


def test_assert_radius_not_diameter_catches_the_mistake():
    with pytest.raises(ValueError, match="DIAMETER"):
        MR.assert_radius_not_diameter(radius=1.0, sigma=1.0)
    assert MR.assert_radius_not_diameter(radius=0.5, sigma=1.0) == 0.5
    with pytest.raises(ValueError):
        MR.assert_radius_not_diameter(radius=0.0, sigma=1.0)


def test_2d_msd_conversion_is_exactly_three_halves():
    """Trap 2. A 2-component MSD converted with 3/2 must reproduce the 3-component
    answer, and the estimator must record that it did the conversion."""
    t = log_lags()
    msd3 = MR.newtonian_msd(t, D=1.0, dim=3)
    msd2 = MR.newtonian_msd(t, D=1.0, dim=2)
    assert np.allclose(MR.msd_2d_to_3d(msd2), msd3, rtol=1e-12)

    a = MR.msd_to_gstar(t, msd3, radius=0.5)
    b = MR.msd_to_gstar(t, msd2, radius=0.5, msd_dim=2)
    v = a["valid"]
    assert np.allclose(a["g_double_prime"][v], b["g_double_prime"][v], rtol=1e-12)
    assert any("3/2" in n for n in b["notes"]), b["notes"]


def test_a_2d_system_is_refused_outright():
    """★ Trap 3, and the one that matters most here: GSER's derivation presumes a
    3D spherical probe, so a 2D *system* is refused even though the MSD could be
    converted. 7 of 8 cases in this repo are 2D."""
    t = log_lags()
    msd = MR.newtonian_msd(t, D=1.0, dim=3)
    with pytest.raises(ValueError, match="3D spherical"):
        MR.msd_to_gstar(t, msd, radius=0.5, system_dim=2)

    r = MR.msd_to_gstar(t, msd, radius=0.5, system_dim=2, allow_2d=True)
    assert any("outside its derivation" in n for n in r["notes"]), r["notes"]


def test_omitting_gamma_costs_about_ten_percent_near_alpha_half():
    """Trap 4 from the distillation, quantified rather than repeated. It says
    ~10 % near alpha=0.5; `Gamma(1.5) = 0.8862`, so the error is 1/0.8862 - 1."""
    t = log_lags()
    msd = MR.powerlaw_msd(t, 1.0, 0.5)
    with_g = MR.msd_to_gstar(t, msd, radius=0.5, apply_gamma=True)
    without = MR.msd_to_gstar(t, msd, radius=0.5, apply_gamma=False)
    v = with_g["valid"]
    ratio = without["g_abs"][v] / with_g["g_abs"][v]
    expected = math.gamma(1.5)                       # 0.886227
    assert np.allclose(ratio, expected, rtol=1e-12)
    err_pct = abs(1.0 / expected - 1.0) * 100
    assert 12.0 < err_pct < 13.0, f"{err_pct:.2f} %"   # 12.84 %, the doc says ~10


# ── 4. the trust band and the alpha mask ────────────────────────────────────

def test_one_decade_is_trimmed_from_each_end():
    t = log_lags(1e-3, 1e3, 240)                     # 6 decades
    r = MR.msd_to_gstar(t, MR.newtonian_msd(t, 1.0, 3), radius=0.5)
    assert math.isclose(r["band_decades"], 4.0, abs_tol=1e-9)
    kept = 1.0 / r["omega"][r["valid"]]
    assert kept.min() >= 1e-3 * 10 - 1e-9
    assert kept.max() <= 1e3 / 10 + 1e-9
    assert any("trust band" in n for n in r["notes"])


def test_superdiffusive_points_are_masked_not_reported():
    """`alpha > 1` is unphysical for GSER at that frequency. The distillation says
    discard those frequencies; they must not silently become moduli."""
    t = log_lags()
    msd = MR.powerlaw_msd(t, 1.0, 1.6)               # ballistic-ish, alpha > 1
    r = MR.msd_to_gstar(t, msd, radius=0.5)
    assert r["n_valid"] == 0, "a superdiffusive MSD must yield no valid moduli"
    assert any("alpha outside" in n for n in r["notes"]), r["notes"]


# ── 5. ★ break the estimator and confirm the tests notice ───────────────────

def test_the_newtonian_check_is_actually_sensitive(monkeypatch):
    """CLAUDE.md: *"When you build a checker, deliberately break it and see."*
    Perturb the constant in the formula by 1 % and the headline assertion must
    fail — otherwise it was passing on tolerance, not on correctness."""
    t = log_lags()
    msd = MR.newtonian_msd(t, D=1.0, dim=3)
    r = MR.msd_to_gstar(t, msd, radius=0.5 * 1.01)   # 1 % wrong radius
    v = r["valid"]
    slope = r["g_double_prime"][v] / r["omega"][v]
    assert not np.allclose(slope, MR.ETA_STAR_NEWTONIAN, rtol=1e-9), \
        "a 1 % error in the radius did not move the recovered eta* -- the " \
        "Newtonian assertion is not testing what it claims"
    assert np.allclose(slope, MR.ETA_STAR_NEWTONIAN / 1.01, rtol=1e-9)


# ── 6. input hygiene ────────────────────────────────────────────────────────

@pytest.mark.parametrize("bad", ["nonmono", "zero_t", "neg_msd", "short", "mismatch"])
def test_bad_input_raises_rather_than_returning_nonsense(bad):
    t = log_lags(1e-3, 1e3, 40)
    msd = MR.newtonian_msd(t, 1.0, 3)
    if bad == "nonmono":
        t = t[::-1]
    elif bad == "zero_t":
        t = t.copy(); t[0] = 0.0
    elif bad == "neg_msd":
        msd = msd.copy(); msd[5] = -1.0
    elif bad == "short":
        t, msd = t[:4], msd[:4]
    elif bad == "mismatch":
        msd = msd[:-1]
    with pytest.raises(ValueError):
        MR.msd_to_gstar(t, msd, radius=0.5)


def test_radius_is_keyword_only():
    """Trap 1 again, at the API level: `a` cannot be passed positionally, so it
    cannot be confused with another length by argument order."""
    t = log_lags(1e-3, 1e3, 40)
    with pytest.raises(TypeError):
        MR.msd_to_gstar(t, MR.newtonian_msd(t, 1.0, 3), 0.5)     # positional
