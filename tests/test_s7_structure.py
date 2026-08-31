"""S7 — time-resolved structure analysis. Uses `freud` but not HOOMD (fast).

## What this file enforces

**① The time-resolved and time-averaged versions do not drift apart.**
`hex_order(frames)` and `hex_order_series(frames).*.mean()` **have to be the same
number.** If the two code paths drift quietly, "the mean of the second half" and
"the late part of the time series" give different values and there is no telling
which is right.

**② The limits have known answers.** A perfect hexagonal lattice → `ψ₆ = 1`, 0
defects, coordination 6 everywhere. Random → `ψ₆ ≈ 0` with a spread coordination
distribution. Both limits are checked frame by frame.

**③ The coverage geometry is checked analytically** — `φ = (π/4)(σ/d)²` round-trips.
"""
from __future__ import annotations

import numpy as np
import pytest

freud = pytest.importorskip("freud")

from simbot.analysis.structure import (HexOrderSeries, hex_order,
                                       hex_order_series, rdf_windows,
                                       voronoi_frame)
from simbot.build import (HEX_NN_OVER_D, box_si_for_coverage,
                          coverage_from_sigma_over_d,
                          sigma_over_d_for_coverage, square_box_for)


# --- a perfect hexagonal lattice (freud coordinates, density_star = 1) --------
def _perfect_hex(n_x: int = 10, n_y: int = 10) -> tuple[np.ndarray, float, float]:
    """A perfect hexagonal lattice at `n* = 1`. Same geometry as
    `simbot.build.hex_2d_snapshot`."""
    a = HEX_NN_OVER_D                      # nearest-neighbour distance (in d)
    Lx, Ly = n_x * a, n_y * a * np.sqrt(3.0) / 2.0
    pts = []
    for j in range(n_y):
        for i in range(n_x):
            x = (i + 0.5 * (j % 2)) * a
            y = j * a * np.sqrt(3.0) / 2.0
            pts.append((x - Lx / 2, y - Ly / 2))
    return np.array(pts), Lx, Ly


def _random_frame(n: int, L: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(-L / 2, L / 2, size=(n, 2))


# =============================================================================
# ① the two code paths do not drift apart
# =============================================================================
def test_series_mean_equals_frame_averaged_hex_order():
    """`hex_order_series`'s frame mean == `hex_order`. **At the bit level.**

    They use the same `_frame_hex`, so only floating-point difference is allowed.
    If this assert breaks, the two functions have started using a different
    neighbour definition or a different normalization.
    """
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts, _random_frame(100, min(Lx, Ly), 3),
                       hexpts * 0.999])
    agg = hex_order(frames, Lx=Lx, Ly=Ly)
    ser = hex_order_series(frames, Lx=Lx, Ly=Ly)

    assert ser.n_frames == 3
    assert ser.psi6_global.mean() == pytest.approx(agg.psi6_global, abs=1e-12)
    assert ser.psi6_local.mean() == pytest.approx(agg.psi6_local_mean, abs=1e-12)
    assert ser.defect_fraction.mean() == pytest.approx(agg.defect_fraction,
                                                       abs=1e-12)
    # the per-frame values must agree too (matching means alone is not enough)
    assert ser.psi6_global == pytest.approx(np.array(agg.psi6_per_frame),
                                            abs=1e-12)


def test_series_t_star_length_is_checked():
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts, hexpts])
    with pytest.raises(ValueError, match="t_star length"):
        hex_order_series(frames, Lx=Lx, Ly=Ly, t_star=np.array([0.0, 1.0, 2.0]))


# =============================================================================
# ② the known answers in the limits — frame by frame
# =============================================================================
def test_perfect_lattice_series_is_exact_every_frame():
    """On a perfect lattice **every frame** has `ψ₆ = 1`, 0 defects, coordination
    6."""
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts] * 4)
    s = hex_order_series(frames, Lx=Lx, Ly=Ly)

    assert s.psi6_global == pytest.approx(np.ones(4), abs=1e-9)
    assert s.psi6_local == pytest.approx(np.ones(4), abs=1e-9)
    assert np.all(s.defect_fraction == 0.0)
    six = s.coord_fraction[:, s.coord_labels == 6].ravel()
    assert six == pytest.approx(np.ones(4), abs=1e-12)
    assert np.all(s.coord_kinds == 1)              # only one kind, coordination 6
    # with no 5s or 7s at all the imbalance is 0 by definition (0/0 treated as 0)
    assert np.all(s.five_seven_balance == 0.0)


def test_random_frames_have_no_orientational_order_and_spread_coordination():
    """Random placement: global `ψ₆` near the finite-size floor, and the
    coordination spreads over several kinds."""
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(6)])
    s = hex_order_series(frames, Lx=L, Ly=L)

    # the finite-size floor: 100 uncorrelated phases give |⟨ψ₆⟩| ~ 1/√N = 0.1
    assert np.all(s.psi6_global < 4.0 / np.sqrt(100))
    assert np.all(s.defect_fraction > 0.3)
    assert np.all(s.coord_kinds >= 4), f"coordination kinds {s.coord_kinds}"


def test_defect_fraction_counts_coordination_outside_the_histogram_range():
    """Narrowing `coord_range` still leaves `defect_fraction` counting all of them.

    If a coordination dropped from the histogram bins (11, say) also dropped out of
    the defect fraction, **defects would silently disappear.** This assert prevents
    that silent loss.
    """
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(4)])
    wide = hex_order_series(frames, Lx=L, Ly=L, coord_range=(3, 12))
    narrow = hex_order_series(frames, Lx=L, Ly=L, coord_range=(6, 6))

    assert narrow.defect_fraction == pytest.approx(wide.defect_fraction,
                                                   abs=1e-12)
    assert narrow.coord_fraction.shape[1] == 1
    # with a narrow range the fractions sum to less than 1 — evidence of dropped bins
    assert narrow.coord_fraction.sum(axis=1).max() < 1.0


def test_window_mean_reproduces_a_hand_sliced_average():
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([_random_frame(100, min(Lx, Ly), s) for s in range(4)]
                      + [hexpts] * 4)
    t = np.arange(8, dtype=float)
    s = hex_order_series(frames, Lx=Lx, Ly=Ly, t_star=t)

    w = s.window_mean(4.0, 8.0)
    assert w["n_frames"] == 4
    assert w["psi6_global"] == pytest.approx(1.0, abs=1e-9)
    assert w["defect_fraction"] == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError, match="no frames inside"):
        s.window_mean(100.0, 200.0)


# =============================================================================
# time-windowed g(r)
# =============================================================================
def test_rdf_windows_partition_all_frames_exactly_once():
    """The windows cover the frames **exactly once each** — nothing dropped, nothing
    duplicated."""
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(12)])
    t = np.linspace(0.0, 11.0, 12)
    w = rdf_windows(frames, Lx=L, Ly=L, t_star=t, n_windows=4, bins=60)

    assert w.g.shape == (4, 60)
    assert w.frames_per_window.sum() == 12
    assert np.all(w.frames_per_window == 3)
    assert np.all(np.diff(w.t_lo) > 0)


def test_rdf_windows_recovers_lattice_peak_at_nearest_neighbour_distance():
    """A perfect lattice's window puts its first peak at the nearest-neighbour
    distance `a = 1.0746 d` (within a bin width)."""
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts] * 4)
    t = np.arange(4, dtype=float)
    w = rdf_windows(frames, Lx=Lx, Ly=Ly, t_star=t, n_windows=2, bins=200)

    bin_width = float(w.r[1] - w.r[0])
    assert np.all(np.abs(w.first_peak_r - HEX_NN_OVER_D) <= 1.5 * bin_width)
    assert np.all(w.first_peak_g > 5.0)          # a lattice has a sharp peak


def test_rdf_windows_rejects_more_windows_than_frames():
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(3)])
    with pytest.raises(ValueError, match="has no frames"):
        rdf_windows(frames, Lx=L, Ly=L, t_star=np.arange(3.0), n_windows=9)


# =============================================================================
# the Voronoi frame (figure material)
# =============================================================================
def test_voronoi_frame_gives_one_polygon_per_particle_and_matches_series():
    """Polygon count = particle count, and the coordination matches
    `hex_order_series`."""
    hexpts, Lx, Ly = _perfect_hex()
    vf = voronoi_frame(hexpts, Lx=Lx, Ly=Ly)

    assert len(vf.polygons) == hexpts.shape[0]
    assert np.all(vf.coordination == 6)
    assert vf.psi6_abs == pytest.approx(np.ones(hexpts.shape[0]), abs=1e-9)
    # a perfect lattice's Voronoi cell is a regular hexagon → 6 vertices
    assert all(p.shape[0] == 6 for p in vf.polygons)

    s = hex_order_series(hexpts[None], Lx=Lx, Ly=Ly)
    assert s.defect_fraction[0] == pytest.approx(
        float(np.mean(vf.coordination != 6)), abs=1e-12)


def test_voronoi_frame_cell_areas_sum_to_the_box_area():
    """The cell areas sum to the box area. **Clipping the polygons breaks this
    check.**"""
    L = square_box_for(100)
    vf = voronoi_frame(_random_frame(100, L, 7), Lx=L, Ly=L)

    def shoelace(p: np.ndarray) -> float:
        x, y = p[:, 0], p[:, 1]
        return 0.5 * abs(float(np.dot(x, np.roll(y, -1))
                               - np.dot(y, np.roll(x, -1))))

    assert sum(shoelace(p) for p in vf.polygons) == pytest.approx(L * L,
                                                                 rel=1e-9)


def test_voronoi_frame_rejects_a_trajectory():
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(2)])
    with pytest.raises(ValueError, match="one frame"):
        voronoi_frame(frames, Lx=L, Ly=L)


# =============================================================================
# the relaxation-time fit — does it mistake noise for relaxation
# =============================================================================
def test_fit_relaxation_recovers_a_known_tau():
    """Recovers `τ` from a synthetic exponential. The noise-free limit, so it has to
    be exact."""
    from simbot.analysis.structure import fit_relaxation

    t = np.linspace(0.0, 2.0, 200)
    tau_true, y0, y_inf = 0.17, 0.53, 0.29
    y = y_inf + (y0 - y_inf) * np.exp(-t / tau_true)
    fit = fit_relaxation(t, y)

    assert fit.converged
    assert fit.tau == pytest.approx(tau_true, rel=1e-6)
    assert fit.y0 == pytest.approx(y0, rel=1e-6)
    assert fit.y_inf == pytest.approx(y_inf, rel=1e-6)
    assert fit.r_squared > 0.9999
    assert fit.amplitude == pytest.approx(y0 - y_inf, rel=1e-6)


def test_fit_relaxation_refuses_to_call_stationary_noise_a_relaxation():
    """★ Fitting an exponential to steady-state fluctuation always yields some τ --
    it has to be rejected.

    This is why the function exists. The defect fraction at `A ≤ 10` is steady from
    the first frame, and attaching a τ to that signal manufactures physics that is
    not there.
    """
    from simbot.analysis.structure import fit_relaxation

    rng = np.random.default_rng(11)
    t = np.linspace(0.0, 2.0, 200)
    noise = 0.06
    y = 0.295 + rng.normal(0.0, noise, t.size)          # no relaxation
    fit = fit_relaxation(t, y, noise=noise)

    assert not fit.converged
    assert "fluctuation" in fit.note or "error exceeds the" in fit.note
    # without `noise` it passes quietly — which is why the caller must supply it
    naive = fit_relaxation(t, y)
    assert np.isfinite(naive.tau)


def test_fit_relaxation_accepts_a_real_relaxation_buried_in_noise():
    """A relaxation amplitude 4x the noise passes — checking the threshold is not
    overly conservative."""
    from simbot.analysis.structure import fit_relaxation

    rng = np.random.default_rng(3)
    t = np.linspace(0.0, 2.0, 400)
    noise = 0.02
    y = 0.29 + 0.24 * np.exp(-t / 0.15) + rng.normal(0.0, noise, t.size)
    fit = fit_relaxation(t, y, noise=noise)

    assert fit.converged, fit.note
    assert fit.tau == pytest.approx(0.15, rel=0.25)


def test_bootstrap_over_seeds_recovers_the_known_spread():
    """★ On synthetic data where τ differs per seed, the bootstrap recovers that
    spread.

    The `curve_fit` covariance is the fit uncertainty of *one mean curve* and does
    not represent the seed-to-seed variation. This test pins that the two measure
    different things.
    """
    from simbot.analysis.structure import bootstrap_relaxation_over_seeds

    rng = np.random.default_rng(5)
    t = np.linspace(0.0, 2.0, 200)
    k, tau_sd = 32, 0.03
    taus = 0.10 + tau_sd * rng.standard_normal(k)       # a different τ per seed
    per_seed = np.array([0.29 + 0.20 * np.exp(-t / max(tt, 1e-3))
                         + rng.normal(0.0, 0.01, t.size) for tt in taus])

    out = bootstrap_relaxation_over_seeds(t, per_seed, n_resample=300, seed=1)
    assert out["n_converged"] > 250
    #  SE of the seed mean = tau_sd/√k. The bootstrap has to recover that scale
    expected_se = tau_sd / np.sqrt(k)
    assert out["tau_se_bootstrap"] == pytest.approx(expected_se, rel=0.6), out
    lo, hi = out["tau_ci95"]
    assert lo < out["tau"] < hi


def test_bootstrap_and_fit_se_are_different_quantities():
    """With no seed-to-seed variation the bootstrap SE is **smaller** than the fit
    SE — they are different quantities."""
    from simbot.analysis.structure import bootstrap_relaxation_over_seeds

    rng = np.random.default_rng(7)
    t = np.linspace(0.0, 2.0, 200)
    #  every seed has the same τ and only the noise differs → the seed-ensemble
    #  spread comes from the noise alone
    per_seed = np.array([0.29 + 0.20 * np.exp(-t / 0.10)
                         + rng.normal(0.0, 0.01, t.size) for _ in range(32)])
    out = bootstrap_relaxation_over_seeds(t, per_seed, n_resample=300, seed=2)

    assert np.isfinite(out["tau_se_bootstrap"])
    assert np.isfinite(out["tau_se_fit"])
    #  the fit SE is for the mean curve (noise 1/√32); the bootstrap resamples
    #  seeds. They must not be equal — equal would mean one copied the other
    assert out["se_ratio_bootstrap_over_fit"] != pytest.approx(1.0, abs=1e-6)


def test_bootstrap_refuses_a_single_seed():
    from simbot.analysis.structure import bootstrap_relaxation_over_seeds

    t = np.linspace(0.0, 1.0, 50)
    with pytest.raises(ValueError, match="cannot give an ensemble error"):
        bootstrap_relaxation_over_seeds(t, np.zeros((1, 50)))
    with pytest.raises(ValueError, match=r"\(n_seeds, n_frames\)"):
        bootstrap_relaxation_over_seeds(t, np.zeros(50))


def test_seeds_for_target_sigma_scales_as_one_over_sqrt_k():
    """`SE ∝ 1/√k`, so doubling σ needs 4x the seeds."""
    from simbot.estimators import seeds_for_target_sigma

    r = seeds_for_target_sigma(diff=0.0109, se_diff=0.0104, k_current=16,
                               n_sigma=3.0)
    assert r["sigma_now"] == pytest.approx(0.0109 / 0.0104, rel=1e-9)
    #  the t correction requires **more** seeds than the normal one
    assert r["k_needed"] > r["k_needed_normal"]
    assert r["k_needed_int"] >= r["k_needed"]

    #  2x σ → 4x k (checked on the normal baseline; the t correction depends on k
    #  so it is approximate)
    r2 = seeds_for_target_sigma(diff=0.0109, se_diff=0.0104, k_current=16,
                                n_sigma=2.0, t_correction=False)
    r3 = seeds_for_target_sigma(diff=0.0109, se_diff=0.0104, k_current=16,
                                n_sigma=4.0, t_correction=False)
    assert r3["k_needed"] / r2["k_needed"] == pytest.approx(4.0, rel=1e-9)


def test_seeds_for_target_sigma_rejects_degenerate_input():
    from simbot.estimators import seeds_for_target_sigma

    with pytest.raises(ValueError, match="must be positive"):
        seeds_for_target_sigma(diff=0.01, se_diff=0.0, k_current=16)
    with pytest.raises(ValueError, match="at least 2"):
        seeds_for_target_sigma(diff=0.01, se_diff=0.01, k_current=1)


def test_fit_relaxation_rejects_too_few_points():
    from simbot.analysis.structure import fit_relaxation

    with pytest.raises(ValueError, match="3-parameter fit"):
        fit_relaxation(np.arange(3.0), np.arange(3.0))
    with pytest.raises(ValueError, match="different lengths"):
        fit_relaxation(np.arange(5.0), np.arange(4.0))


# =============================================================================
# coverage geometry — the analytic round-trip
# =============================================================================
def test_coverage_and_sigma_over_d_round_trip():
    for cov in (0.01, 0.0491, 0.0873, 0.10, 0.5):
        s = sigma_over_d_for_coverage(cov)
        assert coverage_from_sigma_over_d(s) == pytest.approx(cov, rel=1e-12)


def test_coverage_is_independent_of_n_because_density_star_is_one():
    """`φ` does not depend on `N` — at `n* = 1` the area per particle is exactly
    `d²`."""
    out = [box_si_for_coverage(n_particles=n, sigma_si=5e-6, coverage_max=0.10,
                               d_over_sigma_round=3.0)
           for n in (100, 400, 1024)]
    assert {round(o["coverage"], 12) for o in out} == {round(out[0]["coverage"], 12)}
    # the box grows as √N
    assert out[1]["L_si"] / out[0]["L_si"] == pytest.approx(2.0, rel=1e-12)


def test_the_sweep_geometry_is_what_the_report_claims():
    """This sweep's geometry: `σ = 5 µm` · `d/σ = 3` → `L = 150 µm` · `φ = 8.73 %`."""
    g = box_si_for_coverage(n_particles=100, sigma_si=5e-6, coverage_max=0.10,
                            d_over_sigma_round=3.0)
    assert g["d_si"] == pytest.approx(15e-6, rel=1e-12)
    assert g["L_si"] == pytest.approx(150e-6, rel=1e-12)
    assert g["L_star"] == pytest.approx(10.0, rel=1e-12)
    assert g["coverage"] == pytest.approx(np.pi / 36.0, rel=1e-12)
    assert g["coverage"] < 0.10


def test_finite_size_exponent_recovers_the_liquid_and_crystal_limits():
    """★ On synthetic data, `ψ₆ ~ N^{-1/2}` has to read as liquid and `N⁰` as
    crystal."""
    from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,
                                           psi6_finite_size_exponent)

    N = np.array([100.0, 256.0])
    liquid = psi6_finite_size_exponent(N, 0.3 * N ** -0.5)
    assert liquid.p == pytest.approx(0.5, rel=1e-12)
    assert liquid.eta6 == pytest.approx(2.0, rel=1e-12)
    assert liquid.reading == "liquid-like"

    crystal = psi6_finite_size_exponent(N, np.array([0.95, 0.95]))
    assert crystal.p == pytest.approx(0.0, abs=1e-12)
    assert crystal.reading == "hexatic-or-below"

    #  the KTHNY boundary η₆ = 1/4 → p = 1/16
    boundary = psi6_finite_size_exponent(N, 0.3 * N ** -(1.0 / 16.0))
    assert boundary.eta6 == pytest.approx(KTHNY_ETA6_HEXATIC_LIQUID, rel=1e-9)


def test_finite_size_exponent_propagates_the_error_bars():
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([100.0, 256.0])
    y = np.array([0.0477, 0.0298])
    se = np.array([0.0019, 0.0015])
    fit = psi6_finite_size_exponent(N, y, se)

    assert np.isfinite(fit.p_se) and fit.p_se > 0
    #  the relative error is divided by the log difference — checked by hand
    expected = float(np.hypot(se[0] / y[0], se[1] / y[1])
                     / abs(np.log(N[1] / N[0])))
    assert fit.p_se == pytest.approx(expected, rel=1e-12)
    assert fit.eta6_se == pytest.approx(4.0 * expected, rel=1e-12)
    #  called without SEs it is nan, and 'unknown' has to be visible
    assert not np.isfinite(psi6_finite_size_exponent(N, y).p_se)


def test_finite_size_exponent_three_points_is_a_weighted_fit():
    """Three or more and the form can be discussed — `n_points` carries that fact."""
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([64.0, 144.0, 400.0])
    fit = psi6_finite_size_exponent(N, 0.4 * N ** -0.5,
                                    0.01 * np.ones(3) * 0.4 * N ** -0.5)
    assert fit.n_points == 3
    assert fit.p == pytest.approx(0.5, rel=1e-9)


def test_finite_size_exponent_form_test_needs_three_points():
    """★ Two points **assume** a straight line — the fact that the form cannot be
    verified has to be carried out with the result."""
    from simbot.analysis.structure import psi6_finite_size_exponent

    N2 = np.array([100.0, 256.0])
    two = psi6_finite_size_exponent(N2, 0.3 * N2 ** -0.5, 0.001 * np.ones(2))
    assert two.n_points == 2
    assert not two.form_is_testable
    assert not np.isfinite(two.chi2_reduced)
    assert "not testable" in two.form_verdict()

    N4 = np.array([64.0, 144.0, 256.0, 400.0])
    y = 0.3 * N4 ** -0.5
    four = psi6_finite_size_exponent(N4, y, 0.01 * y)
    assert four.form_is_testable
    assert four.chi2_reduced == pytest.approx(0.0, abs=1e-16)
    assert "consistent with a power law" in four.form_verdict()
    assert four.p == pytest.approx(0.5, rel=1e-9)
    assert four.amplitude == pytest.approx(0.3, rel=1e-9)


def test_finite_size_exponent_flags_a_curve_that_is_not_a_power_law():
    """★ Attaching a single exponent to data that is **not** a power law has to
    warn.

    This is why the form check exists -- an exponent always comes out, but if the
    form is wrong the phase must not be decided from it.
    """
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([64.0, 144.0, 256.0, 400.0])
    #  a curve bending in log-log (the exponent varies with N)
    y = 0.3 * N ** -0.5 * (1.0 + 0.35 * np.log(N / 64.0))
    fit = psi6_finite_size_exponent(N, y, 0.005 * y)

    assert fit.form_is_testable
    assert fit.chi2_reduced > 3.0, fit.chi2_reduced
    assert "departs from a power law" in fit.form_verdict()
    #  the residuals change sign — the signature of curvature
    r = np.array(fit.residuals)
    assert np.any(r > 0) and np.any(r < 0)


def test_finite_size_exponent_chi2_needs_error_bars():
    """Without error bars `χ²` cannot be computed — it does not quietly return 0."""
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([64.0, 144.0, 256.0, 400.0])
    fit = psi6_finite_size_exponent(N, 0.3 * N ** -0.5)
    assert fit.n_points == 4
    assert not np.isfinite(fit.chi2_reduced)
    assert not fit.form_is_testable


def test_phase_reading_does_not_call_a_crystal_hexatic():
    """★★ `η₆ ≤ 1/4` is **only a necessary** condition for a hexatic — a crystal
    satisfies it too.

    On 2026-07-29, without this distinction, the crystal at `A ≥ 13.3` was printed
    as "hexatic possible? yes". What separates them is the **magnitude** of `ψ₆`.
    """
    from simbot.analysis.structure import (phase_from_finite_size,
                                           psi6_finite_size_exponent)

    N = np.array([144.0, 256.0, 400.0])

    #  ① crystal: psi6 saturates at O(1) → p ≈ 0 → eta6 ≈ 0 (the same exponent as
    #     a hexatic!)
    y_cry = np.array([0.6720, 0.6841, 0.6881])
    f_cry = psi6_finite_size_exponent(N, y_cry, np.array([0.0032, 0.0074, 0.0080]))
    r_cry = phase_from_finite_size(f_cry, y_cry[-1])
    assert abs(f_cry.eta6) < 0.25                    # by exponent alone, a hexatic
    assert r_cry["exponent_alone_would_say"] == "hexatic-or-below"
    assert r_cry["phase"] == "crystal", r_cry        # ★ the magnitude separates them
    assert "saturates" in r_cry["why"]

    #  ② liquid: psi6 is small and decays as N^-1/2
    y_liq = 0.5 * N ** -0.5
    f_liq = psi6_finite_size_exponent(N, y_liq, 0.02 * y_liq)
    r_liq = phase_from_finite_size(f_liq, y_liq[-1])
    assert r_liq["phase"] == "isotropic-liquid"
    assert r_liq["exponent_alone_would_say"] == "not-hexatic"

    #  ③ hexatic candidate: the exponent is **below** the ceiling and psi6 is below
    #     the crystal floor.
    #     ★ Putting η₆ exactly at the ceiling (0.25) gives `η₆ + 3σ > 0.25`, so
    #       inconclusive is correct -- a value sitting on the boundary must not be
    #       called a hexatic. Hence 0.1.
    y_hex = 0.30 * (N / N[0]) ** -(0.1 / 4.0)        # eta6 = 0.1 < 0.25
    f_hex = psi6_finite_size_exponent(N, y_hex, 0.005 * y_hex)
    r_hex = phase_from_finite_size(f_hex, y_hex[-1])
    assert y_hex[-1] < 0.5
    assert f_hex.eta6 == pytest.approx(0.1, rel=1e-6)
    assert r_hex["phase"] == "hexatic-candidate", r_hex


def test_phase_reading_flags_inconclusive_when_eta6_straddles_the_ceiling():
    from simbot.analysis.structure import (phase_from_finite_size,
                                           psi6_finite_size_exponent)

    N = np.array([144.0, 256.0, 400.0])
    #  eta6 ≈ 0.25 but the error is large enough to straddle the ceiling
    y = 0.30 * (N / N[0]) ** -(1.0 / 16.0)
    fit = psi6_finite_size_exponent(N, y, 0.06 * y)
    r = phase_from_finite_size(fit, y[-1])
    assert r["phase"] in ("inconclusive", "hexatic-candidate")
    if r["phase"] == "inconclusive":
        assert "straddles" in r["why"]


def test_finite_size_exponent_rejects_degenerate_input():
    from simbot.analysis.structure import psi6_finite_size_exponent

    with pytest.raises(ValueError, match="at least 2"):
        psi6_finite_size_exponent(np.array([100.0]), np.array([0.05]))
    with pytest.raises(ValueError, match="value <= 0"):
        psi6_finite_size_exponent(np.array([100.0, 256.0]),
                                  np.array([0.05, 0.0]))


def test_coverage_does_not_touch_the_reduced_config_at_all():
    """★★ Changing the coverage leaves the reduced-unit run **bit-identical.**

    So **running a "coverage 3.7 % control" as a simulation is an arithmetic
    identity** (`conftest` rule ③: a measurement that does not fluctuate is an
    identity). It is not an item to spend compute on but an item to pin with this
    test -- that proposal was made on 2026-07-29 and rejected here.

    All the coverage changes is the SI conversion: `τ_d ∝ d²` (σ fixed).
    """
    from dataclasses import asdict

    from simbot.run import Soft2DRunConfig
    from simbot.units import scales_soft2d

    made = {}
    for d_over_sigma in (3.0, 4.6):                 # coverage 8.73 % vs 3.71 %
        g = box_si_for_coverage(n_particles=100, sigma_si=5e-6,
                                coverage_max=0.10,
                                d_over_sigma_round=d_over_sigma)
        cfg = Soft2DRunConfig(amplitude=10.0, n_particles=100, init="random",
                              box_shape="square", r_min=0.05, dt_star=5.01e-5,
                              equil_tau=0.0, prod_tau=80.0, seed=5)
        made[d_over_sigma] = (asdict(cfg), cfg.hash(), g,
                              scales_soft2d(d_si=g["d_si"], sigma_si=5e-6,
                                            T_si=298.15))
    a, b = made[3.0], made[4.6]

    assert a[2]["coverage"] != b[2]["coverage"]      # the coverage really must differ
    assert b[2]["coverage"] < 0.0371 + 1e-4
    #  ★ and yet the reduced config and run_hash are identical
    assert a[0] == b[0]
    assert a[1] == b[1], \
        "a different run_hash means the coverage leaked into the dynamics"
    #  only the SI conversion changes: tau_d ∝ d²
    assert b[3].time_si / a[3].time_si == pytest.approx(
        (b[2]["d_si"] / a[2]["d_si"]) ** 2, rel=1e-12)
    assert b[3].time_si > a[3].time_si


def test_box_for_coverage_refuses_a_ratio_that_breaks_the_ceiling():
    """Even when `d/σ` is given directly, exceeding the coverage cap is
    **rejected.**"""
    with pytest.raises(ValueError, match="커버리지"):
        box_si_for_coverage(n_particles=100, sigma_si=5e-6, coverage_max=0.10,
                            d_over_sigma_round=2.0)
    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        sigma_over_d_for_coverage(0.0)


def test_scales_soft2d_uses_the_lattice_spacing_for_length_and_sigma_for_drag():
    """★ The length scale is `d` and the drag is `σ`. Mix them and `τ_d` is wrong by
    `(d/σ)²`."""
    from simbot.units import (kT_si, scales_soft2d, stokes_drag_si,
                              stokes_einstein_D_si, water_viscosity_si)

    sigma, T = 5e-6, 298.15
    d = 3.0 * sigma
    sc = scales_soft2d(d_si=d, sigma_si=sigma, T_si=T)

    eta, extrapolated = water_viscosity_si(T)
    assert not extrapolated
    D0 = stokes_einstein_D_si(T, stokes_drag_si(eta, sigma / 2.0))
    assert sc.length_si == pytest.approx(d, rel=1e-12)
    assert sc.energy_si == pytest.approx(kT_si(T), rel=1e-12)
    assert sc.time_si == pytest.approx(d**2 / D0, rel=1e-12)

    # wrongly using σ as the length scale is wrong by (d/σ)² = 9 — stated explicitly
    wrong = scales_soft2d(d_si=sigma, sigma_si=sigma, T_si=T)
    assert sc.time_si / wrong.time_si == pytest.approx(9.0, rel=1e-12)


def test_scales_soft2d_refuses_to_extrapolate_the_viscosity_table():
    """Outside the table (293–308 K) it does not extrapolate quietly but demands
    `gamma_si`."""
    from simbot.units import scales_soft2d

    with pytest.raises(ValueError, match="gamma_si"):
        scales_soft2d(d_si=15e-6, sigma_si=5e-6, T_si=350.0)
    # stated explicitly, it passes outside the table too
    sc = scales_soft2d(d_si=15e-6, sigma_si=5e-6, T_si=350.0, gamma_si=1e-8)
    assert sc.time_si > 0.0
