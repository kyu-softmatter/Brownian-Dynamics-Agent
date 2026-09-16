"""The sedimentation EOS estimator, against answers it must already know.

★ Why this file exists. The whole measurement of `sediment-pmma-3d` is `Z(phi)`
read out of the density profile. If that read is wrong, every number downstream
is wrong and **looks like physics** -- a deviation from Carnahan-Starling is
exactly what the run is trying to measure, so a bug in the estimator is
indistinguishable from a result.

`.claude/rules/deterministic-core.md`: *"Test a new analysis estimator against
synthetic data whose answer is known."*

★★ **It earned its keep immediately.** The first version of `eos_from_profile`
used a plain `cumsum`, which integrates from each bin's LOWER EDGE rather than
its midpoint and so overestimates the column above `z` by half a bin. Measured:

    l_g = 13.375, bin = 0.50  ->  Z high by +1.65 %      <- the production setting
    l_g =  4.000, bin = 1.00  ->  Z high by +12.79 %

**+1.65 % against a pre-registered 2 % band.** The run would have reported a
binning artefact as a deviation from Carnahan-Starling, and the number would have
looked entirely reasonable. Taking half of the bin containing `z` plus the whole
bins above leaves -0.22 %, and that residual turns out not to be discretisation
at all -- it is the finite cell, which has an exact answer of its own.

⚠ And two of this file's own first assertions were wrong while the code was
right: the Carnahan-Starling values at phi = 0.3 and 0.4 were written from memory
(4.04956 / 7.66049 against the correct 3.973761 / 6.925926 -- exactly what
`CLAUDE.md` forbids), and the range asserted for
`closest_sampled_separation` forgot that shifted WCA has `U = eps` at
`r = sigma_LJ`, so `U = 8 kT` at `eps = 1` is necessarily INSIDE sigma_LJ. Both
are now checked against closed forms instead of recollection.
"""
from __future__ import annotations

import importlib.util
import pathlib

import numpy as np
import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "sediment_3d", ROOT / "cases" / "sediment_3d.py")
SED = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(SED)


def _grid(lz, bin_d):
    """Bin midpoints, and `lz` SNAPPED to the grid.

    ⚠ Two of this file's rows failed until `lz` was snapped: with
    `lz = 12 * 13.375 = 160.5` and `bin = 1.0` the grid stops at 160.0 while the
    exact formula below uses 160.5, so the test compared the estimator against a
    cell half a bin taller than the one it was given. The residual jumped from
    `(bin/l_g)^2/8` to 50x that, and it was the TEST that was wrong.
    """
    nb = int(round(lz / bin_d))
    return (np.arange(nb) + 0.5) * bin_d, nb * bin_d


def _ideal_Z(mid, lz, l_g):
    """The EXACT answer for an ideal gas in a cell of finite height.

    `Z(z) = (1/(l_g n(z))) * integral_z^Lz n dz'` with `n = n0 exp(-z/l_g)` gives
    `1 - exp(-(Lz - z)/l_g)` — not 1. The difference is the mass that an infinite
    column would have above the cell, and it is not an approximation: the
    simulated cell has a wall at the top, so there genuinely is nothing above it.
    """
    return 1.0 - np.exp(-(lz - mid) / l_g)


@pytest.mark.parametrize("l_g", [4.0, 13.375, 40.0])
@pytest.mark.parametrize("bin_d", [0.25, 0.5, 1.0])
def test_an_ideal_gas_profile_gives_the_exact_finite_cell_Z(l_g, bin_d):
    """★ The load-bearing test, over three decades of `bin/l_g`. Independent of
    `l_g`, of the amplitude and of the bin width, so a factor that depends on any
    of the three is a bug this catches."""
    mid, lz = _grid(12.0 * l_g, bin_d)
    Z = SED.eos_from_profile(0.1 * np.exp(-mid / l_g), bin_d, l_g)
    want = _ideal_Z(mid, lz, l_g)
    err = float(np.max(np.abs(Z - want)))
    #  ★ The residual is not asserted as a flat tolerance -- it obeys a law.
    #    Measured over three decades of bin/l_g (0.0063 to 0.5): the maximum
    #    deviation is (bin/l_g)^2 / 8, to within a few percent. Asserting the LAW
    #    catches a change in the estimator's order of accuracy, which a flat
    #    tolerance loose enough to pass everywhere would not.
    q = bin_d / l_g
    assert err == pytest.approx(q * q / 8.0, rel=0.10), (l_g, bin_d, err, q * q / 8)


def test_deep_in_the_cell_Z_is_one():
    """The statement people actually quote. Holds only where the top wall is far,
    which is why the fit and EOS windows have an upper bound."""
    l_g, bin_d = 13.375, 0.5
    mid, lz = _grid(125.0, bin_d)
    Z = SED.eos_from_profile(0.1 * np.exp(-mid / l_g), bin_d, l_g)
    deep = (lz - mid) > 6.0 * l_g
    assert deep.any()
    assert np.allclose(Z[deep], 1.0, atol=3e-3), Z[deep].min()


def test_the_half_bin_correction_is_what_makes_that_true():
    """★ Guard on the guard. A plain `cumsum` -- the first implementation -- is
    high by `bin/(2 l_g)`. Asserted so that reverting the correction fails here
    with the number that explains why, rather than surfacing as a physics result.
    """
    l_g, bin_d = 13.375, 0.5
    mid, lz = _grid(125.0, bin_d)
    phi = 0.1 * np.exp(-mid / l_g)
    naive = np.cumsum(phi[::-1])[::-1] * bin_d / (l_g * phi)
    good = SED.eos_from_profile(phi, bin_d, l_g)
    deep = (lz - mid) > 6.0 * l_g
    assert naive[deep].mean() - 1.0 == pytest.approx(bin_d / (2 * l_g), rel=0.15)
    assert abs(good[deep].mean() - 1.0) < 0.2 * (bin_d / (2 * l_g))


def test_the_amplitude_cancels():
    l_g, bin_d = 13.375, 0.5
    mid, lz = _grid(125.0, bin_d)
    a = SED.eos_from_profile(0.01 * np.exp(-mid / l_g), bin_d, l_g)
    b = SED.eos_from_profile(0.30 * np.exp(-mid / l_g), bin_d, l_g)
    assert np.allclose(a, b, atol=1e-12)


@pytest.mark.parametrize("factor", [0.5, 2.0])
def test_a_wrong_l_g_shows_up_as_a_proportional_offset(factor):
    """If the estimator were insensitive to `l_g`, the tests above would pass for
    the wrong reason."""
    l_g, bin_d = 13.375, 0.5
    mid, lz = _grid(125.0, bin_d)
    Z = SED.eos_from_profile(0.1 * np.exp(-mid / l_g), bin_d, l_g * factor)
    deep = (lz - mid) > 6.0 * l_g
    assert np.allclose(Z[deep], 1.0 / factor, atol=5e-3), (factor, Z[deep].mean())


def test_a_flat_profile_is_not_an_ideal_gas_and_the_estimator_says_so():
    """A uniform profile has no hydrostatic balance, so `Z` must NOT come out 1 --
    it rises linearly with the column above, and with the half-bin correction the
    estimator gets that EXACTLY. Pins sensitivity to the shape, not only to the
    normalisation."""
    l_g, bin_d = 13.375, 0.5
    mid, lz = _grid(60.0, bin_d)
    Z = SED.eos_from_profile(np.full(len(mid), 0.05), bin_d, l_g)
    assert np.allclose(Z, (lz - mid) / l_g, rtol=1e-12)
    assert Z[0] > 4.0, Z[0]


def test_carnahan_starling_against_the_virial_expansion_not_from_memory():
    """`Z_CS = 1 + 4 phi + 10 phi^2 + 18 phi^3 + ...` — the hard-sphere virial
    coefficients in `phi` units. Checking the polynomial against its own
    expansion is independent; checking it against remembered values is not, and
    the first version of this test got two of those values wrong."""
    phi = np.array([1e-6, 1e-5, 1e-4])
    z = SED.z_carnahan_starling(phi)
    assert np.allclose((z - 1.0) / phi, 4.0, rtol=2e-3), (z - 1.0) / phi
    mid = np.array([0.01, 0.02])
    zm = SED.z_carnahan_starling(mid)
    assert np.allclose(zm, 1 + 4 * mid + 10 * mid ** 2, rtol=1e-3), zm
    assert SED.z_carnahan_starling(np.array([0.0]))[0] == pytest.approx(1.0)
    dense = SED.z_carnahan_starling(np.linspace(0.0, 0.45, 40))
    assert np.all(np.diff(dense) > 0)


def test_the_barker_henderson_mapping_is_the_discriminator_it_is_claimed_to_be():
    """The pre-registration rests on CS(phi_eff) and CS(phi_nominal) being
    separable. If they were not, the run could not answer its own question, so
    the separation is asserted rather than assumed."""
    from simbot.cutoff import barker_henderson_diameter
    d_bh = barker_henderson_diameter(1.0, sigma_lj=SED.SIGMA_LJ)
    assert d_bh == pytest.approx(0.9048, abs=5e-4), d_bh
    phi_wall = 0.15889
    gap = (SED.z_carnahan_starling(np.array([phi_wall * d_bh ** 3]))[0]
           / SED.z_carnahan_starling(np.array([phi_wall]))[0] - 1.0)
    assert gap == pytest.approx(-0.1705, abs=2e-3), gap
    assert abs(gap) > 8 * 0.02, "the 2 % band must be small against the gap"


def test_the_dt_gate_is_read_from_the_force_the_run_can_sample():
    """`closest_sampled_separation` is what set `dt`; a wrong root would move the
    gate silently. Solved by bisection, so it is checked against the potential —
    and against the shifted-WCA identity `U(sigma_LJ) = eps`, which is what makes
    the root's position obvious."""
    for eps in (1.0, 10.0):
        assert SED.u_wca(SED.SIGMA_LJ, eps) == pytest.approx(eps, rel=1e-12)
        assert SED.u_wca(SED.R_CUT, eps) == pytest.approx(0.0, abs=1e-12)
        r = SED.closest_sampled_separation(eps, u_kT=8.0)
        assert SED.u_wca(r, eps) == pytest.approx(8.0, rel=1e-6), (eps, r)
        assert 0.5 < r < SED.R_CUT
        #  U = 8 kT sits inside sigma_LJ exactly when eps < 8
        assert (r < SED.SIGMA_LJ) == (eps < 8.0), (eps, r)


# ── the Carnahan-Starling hydrostatic profile, and sampling from it ────────
#
# ★ Why this section exists. The primary arm's INITIAL CONDITION is this profile,
#   and that is only defensible if the profile is what it claims to be. Measured
#   on an ideal-gas start instead: 8 tau_sed moved phi(0) just 77 % of the way
#   from 0.1589 to 0.1031 and the drift check still read 6.3 sigma, so ~25
#   tau_sed per seed would be needed -- 4.7 h for the primary arm alone.

def test_the_cs_profile_satisfies_its_own_hydrostatic_equation():
    """★ The load-bearing test: check the ODE, not the integrator. Substituting
    the returned profile back into `d(phi Z)/dz = -phi/l_g` is independent of how
    it was obtained."""
    l_g, lz, phi_mean = 13.375, 125.0, 0.017
    z, phi = SED.cs_hydrostatic_profile(l_g, lz, phi_mean)
    Z = SED.z_carnahan_starling(phi)
    lhs = np.gradient(phi * Z, z)
    rhs = -phi / l_g
    m = phi > 1e-5                     # where the gradient is resolved
    assert np.allclose(lhs[m], rhs[m], rtol=2e-3), (
        float(np.max(np.abs(lhs[m] / rhs[m] - 1))))


def test_the_cs_profile_has_the_mean_volume_fraction_it_was_asked_for():
    for phi_mean in (0.005, 0.017, 0.05):
        z, phi = SED.cs_hydrostatic_profile(13.375, 125.0, phi_mean)
        assert float(np.trapezoid(phi, z)) / 125.0 == pytest.approx(phi_mean, rel=1e-6)


def test_the_cs_profile_is_less_peaked_than_the_ideal_gas_and_that_is_the_point():
    """Z > 1 extends the column, so phi(0) must come out BELOW the ideal-gas
    value. The direction is asserted because it is the prediction P5 records, and
    a sign error here would invert the whole comparison."""
    l_g, lz, phi_mean = 13.375, 125.0, 0.017
    z, phi = SED.cs_hydrostatic_profile(l_g, lz, phi_mean)
    ideal0 = phi_mean * (lz / l_g) / (1.0 - np.exp(-lz / l_g))
    assert ideal0 == pytest.approx(0.15889, abs=1e-4)
    assert phi[0] == pytest.approx(0.10306, abs=5e-4), phi[0]
    assert phi[0] < ideal0
    #  and it must OVERSHOOT the ideal gas high up: the mass has to go somewhere
    hi = z > 3.0 * l_g
    ideal_hi = ideal0 * np.exp(-z[hi] / l_g)
    assert np.all(phi[hi] > ideal_hi)


def test_the_two_fit_windows_measure_different_things():
    """★ This is the mistake the case made first. Over the DENSE window the
    CS-equilibrium profile's apparent decay length is 15.44 d, so an
    `implementation_check` against 13.375 there calls a correct result a bug. Over
    the TAIL window it is 13.50 d, +0.9 %, which a 5 % check survives."""
    l_g, lz = 13.375, 125.0
    z, phi = SED.cs_hydrostatic_profile(l_g, lz, 0.017)

    def apparent(lo, hi):
        m = (z >= lo) & (z <= hi) & (phi > 1e-12)
        return -1.0 / np.polyfit(z[m], np.log(phi[m]), 1)[0]

    dense = apparent(SED.FIT_LO, 6.0 * l_g)
    tail = apparent(SED.TAIL_LO, SED.TAIL_HI)
    assert dense == pytest.approx(15.44, abs=0.05), dense
    assert tail == pytest.approx(13.50, abs=0.05), tail
    assert abs(tail / l_g - 1) < 0.05, "the tail window must survive a 5 % check"
    assert abs(dense / l_g - 1) > 0.10, (
        "the dense window must NOT -- if it did, the split would be pointless")


def test_the_height_sampler_reproduces_the_profile_it_was_built_from():
    """Kolmogorov-Smirnov against the profile's own trapezoid CDF, which is
    normalisation-free.

    ⚠ A `cumsum` CDF -- the first implementation -- left the sampled mean 0.44 %
      low and KS sitting ON the 95 % critical value (D = 0.00291 against 0.00304
      at N = 200,000). Marginal, not correct, and the same half-interval error
      class as the half-bin in `eos_from_profile`.
    """
    from scipy import stats
    l_g, lz, phi_mean, margin = 13.375, 125.0, 0.017, 0.5
    z, phi = SED.cs_hydrostatic_profile(l_g, lz, phi_mean)
    keep = (z >= margin) & (z <= lz - margin)
    zk, w = z[keep], phi[keep]
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(zk))])
    cdf /= cdf[-1]
    n = 200_000
    h = SED.cs_height_sampler(l_g, lz, phi_mean, margin=margin)(
        np.random.default_rng(0), n)
    ks = stats.kstest(h, lambda x: np.interp(x, zk, cdf))
    assert ks.statistic < 1.36 / np.sqrt(n), (ks.statistic, 1.36 / np.sqrt(n))
    assert h.min() >= margin and h.max() <= lz - margin


# ── the stationarity split ─────────────────────────────────────────────────

def test_the_half_window_split_is_on_the_timestep_not_the_frame_count():
    """★ `bdbot.run.execute` does NOT sample every `sample_every` steps -- it caps
    the interval. Measured 2026-09-16: `sample_every = 29815` and the runner
    sampled every 10,000, so the case's `n_prod // sample_every` predicted 12
    frames where 36 arrived, and the "halves" came out 6 and 30.

    The statistic still returned an ordinary-looking 1.6356. Nothing failed, and
    `profile_halves_chi2_nu` is the whole basis of the campaign's verdict -- so
    this is the class where a check's success is indistinguishable from its
    failure. The split is now on the timestep and the balance is ASSERTED.
    """
    src = (ROOT / "cases" / "sediment_3d.py").read_text()
    #  the arithmetic that caused it must not come back
    assert "n_prod_frames" not in src, (
        "a predicted frame count is back in the case; the runner's cadence is not "
        "this case's to predict")
    assert 'split_step' in src and 'timestep <= acc["split_step"]' in src
    assert 'int(Nm["n_eq"]) + int(Nm["n_prod"]) // 2' in src


def test_the_balance_guard_raises_on_a_lopsided_split():
    """The guard is inside the run, so it is exercised here against the numbers
    that actually occurred rather than against its own source."""
    def guard(f1r, f2r, split_step=178890):
        tot = f1r + f2r
        if tot >= 8 and abs(f1r - f2r) > max(2, 0.15 * tot):
            raise RuntimeError(f"{f1r}/{f2r}")
        return True

    #  the measured failure, and the production-length version of it
    with pytest.raises(RuntimeError):
        guard(6, 30)
    with pytest.raises(RuntimeError):
        guard(80, 397)
    #  the measured fix, and an exactly even split
    assert guard(17, 19) and guard(18, 18) and guard(238, 239)
    #  a short run below the guard's floor is not rejected -- it has no halves
    #  worth testing and the statistic reports NaN from the 5-bin minimum instead
    assert guard(2, 5)


# ── the block accumulation, which is where every error bar now comes from ──

def test_block_reduction_holds_the_window_fixed_across_blocks():
    """★ A block average is ONE estimator applied to subsets. Letting each block
    re-derive its own window from its own counts makes it sixteen estimators.

    Measured: with the `>= 200` count floor applied per block, at 1/16 of the
    frames the EOS window came out empty and `Z_dev_wmean_block_sem` returned NaN
    on a short run. The masks now come from the pooled accumulation.
    """
    src = (ROOT / "cases" / "sediment_3d.py").read_text()
    body = src[src.index("def reduce_profile"):src.index("block_stats = [reduce_profile")]
    #  the pooled masks are used, and no per-block count floor is re-derived
    assert "okb = ok & np.isfinite(ze)" in body
    assert "dil = dilute & np.isfinite(ze)" in body
    assert "counts >= 200" not in body, "a per-block count floor is back"
    assert "counts >= 50" not in body, "a per-block dilute floor is back"


def test_a_single_empty_bin_does_not_void_a_whole_block():
    """⚠ One zero-count bin inside the pooled window makes `eos_from_profile` NaN
    there, and a plain `.sum()` propagates it to the entire weighted mean.
    Measured: one seed's block 0 came back NaN while its neighbours were fine, so
    the SEM was NaN from a single bin. Dropping the bin is the same estimator on
    the bins that have data; voiding the block is not."""
    src = (ROOT / "cases" / "sediment_3d.py").read_text()
    body = src[src.index("def reduce_profile"):src.index("block_stats = [reduce_profile")]
    assert "fe, fn = np.isfinite(de), np.isfinite(dn)" in body
    assert "eos_bins_dropped" in body, "the drop count must be recorded, not silent"

    #  and the arithmetic itself, on a window with one hole
    w = np.array([100.0, 200.0, 50.0, 300.0])
    dev = np.array([1.0, 2.0, np.nan, 4.0])
    fe = np.isfinite(dev)
    got = float((w[fe] * dev[fe]).sum() / w[fe].sum())
    assert got == pytest.approx((100 * 1 + 200 * 2 + 300 * 4) / 600.0)
    #  the naive version is the bug
    assert not np.isfinite(float((w * dev).sum() / w.sum()))


def test_the_block_sem_is_NaN_rather_than_zero_when_it_cannot_be_formed():
    """★ A zero error bar is worse than no error bar: it passes every sigma-based
    check. With fewer than 3 finite blocks the SEM must be NaN, and the primary
    statistic's check must then be INCONCLUSIVE rather than passed."""
    src = (ROOT / "cases" / "sediment_3d.py").read_text()
    assert 'if len(v) < 3:\n                return float("nan"), len(v)' in src
    #  and the observable refuses to invent a sigma
    assert "sigma=(sem_primary if np.isfinite(sem_primary) else None)" in src
    assert "tol_sigma=(3.0 if np.isfinite(sem_primary) else None)" in src


def test_sixteen_blocks_and_the_reason_is_recorded():
    assert SED.N_BLOCKS == 16
    doc = (ROOT / "cases" / "sediment_3d.py").read_text()
    i = doc.index("N_BLOCKS = 16")
    preamble = doc[max(0, i - 900):i]
    #  the number of blocks must carry its own justification, because the blocks
    #  are correlated at the measured tau_int and that has to be visible
    assert "nu = 15" in preamble and "correlated" in preamble
