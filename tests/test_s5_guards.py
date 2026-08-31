"""S5 — runtime guards. The pure functions are fast; only the configurational
temperature (B8) uses HOOMD.

The principle for testing a guard: **test only the passing case and you will not
notice when the guard dies.** The cases where it must fire are always tested too.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from simbot.guards import (
    assert_statistic_fluctuates, check_finite, check_inside_box,
    check_step_displacements, configurational_temperature,
)

from .conftest import TRAP_DIM, se_of_mean, sigma_away


# =============================================================================
# configurational temperature — pure functions
# =============================================================================
def test_configurational_temperature_recovers_kT_analytically():
    """In a harmonic trap, kT_conf = k^2<r^2>/(d k) = k<r^2>/d.

    With d=2, k=1, <r^2> = 2 (i.e. <x^2>=1), kT_conf = 1.
    """
    rng = np.random.default_rng(0)
    n, d, k, kT = 200000, 2, 1.0, 1.0
    # Boltzmann distribution: per-component SD sqrt(kT/k)
    pos = rng.normal(0.0, math.sqrt(kT / k), size=(n, d))
    forces = -k * pos
    kT_conf = configurational_temperature(forces, laplacian_U_total=d * k)
    assert kT_conf == pytest.approx(kT, rel=0.01)


@pytest.mark.parametrize("kT", [0.5, 1.0, 2.5])
def test_configurational_temperature_tracks_input_kT(kT):
    rng = np.random.default_rng(1)
    n, d, k = 200000, 3, 2.0
    pos = rng.normal(0.0, math.sqrt(kT / k), size=(n, d))
    kT_conf = configurational_temperature(-k * pos, laplacian_U_total=d * k)
    assert kT_conf == pytest.approx(kT, rel=0.01)


def test_configurational_temperature_rejects_bad_laplacian():
    with pytest.raises(ValueError):
        configurational_temperature(np.ones((3, 2)), laplacian_U_total=0.0)
    with pytest.raises(ValueError):
        configurational_temperature(np.ones((3, 2)), laplacian_U_total=-1.0)


# =============================================================================
# displacement per step
# =============================================================================
def test_displacement_guard_passes_for_small_steps():
    dr = np.full((100, 3), 0.001)
    r = check_step_displacements(dr, sigma=1.0, max_frac=0.10)
    assert r.passed and r.n_exceeding == 0 and r.note == ""


def test_displacement_guard_fires_on_explosion():
    """★ Does the guard actually fire?"""
    dr = np.zeros((100, 3))
    dr[7] = [0.5, 0.0, 0.0]              # 0.5 sigma — a sign of blow-up
    r = check_step_displacements(dr, sigma=1.0, max_frac=0.10)
    assert not r.passed
    assert r.n_exceeding == 1
    assert "dt too large" in r.note        # message text moved with the impl to bdbot.health
    assert r.max_over_sigma == pytest.approx(0.5)


def test_displacement_guard_reports_max_over_rms():
    """With uniform noise, max/rms must be visible as being near sqrt(3)."""
    rng = np.random.default_rng(3)
    dr = rng.uniform(-1, 1, size=(50000, 1)) * 0.01
    dr = np.hstack([dr, np.zeros((len(dr), 2))])
    r = check_step_displacements(dr, sigma=1.0, max_frac=1.0)
    assert r.passed
    assert r.max_over_rms == pytest.approx(math.sqrt(3), rel=0.02)


# =============================================================================
# finiteness · boundaries
# =============================================================================
def test_check_finite_catches_nan_and_inf():
    ok, fails = check_finite(pos=np.ones((3, 3)), force=np.ones((3, 3)))
    assert ok and not fails

    bad = np.ones((3, 3))
    bad[1, 1] = np.nan
    worse = np.ones((3, 3))
    worse[0, 0] = np.inf
    ok, fails = check_finite(pos=bad, force=worse)
    assert not ok
    assert any("pos" in f for f in fails) and any("force" in f for f in fails)


def test_check_inside_box():
    ok, n = check_inside_box(np.array([[1.0, 1.0, 0.0]]), [10, 10, 10], dims=2)
    assert ok and n == 0
    ok, n = check_inside_box(np.array([[6.0, 0.0, 0.0], [0.0, -7.0, 0.0]]),
                             [10, 10, 10], dims=2)
    assert not ok and n == 2


# =============================================================================
# the device against "mistaking an identity for a measurement"
# =============================================================================
def test_fluctuation_check_passes_for_real_statistic():
    rng = np.random.default_rng(5)
    assert_statistic_fluctuates(rng.normal(1.0, 0.1, 500), "real")


def test_fluctuation_check_catches_arithmetic_identity():
    """★ Pin a failure actually experienced on 2026-07-28 as a test.

    Subtract the mean from the displacements and then measure the cross
    correlation, and cross/auto = -1/(n-1) comes out identically -- the SD over
    200 repeats was 6.7e-20. The result was plausible enough that it nearly passed.
    """
    n = 1000
    identity_values = np.full(300, -1.0 / (n - 1))     # a non-fluctuating "measurement"
    with pytest.raises(AssertionError, match="arithmetic identity"):
        assert_statistic_fluctuates(identity_values, "cross/auto")


def test_fluctuation_check_reproduces_the_original_bug(hoomd_mod):
    """Confirm directly that subtracting the mean really does give -1/(n-1)
    (no HOOMD needed, pure arithmetic)."""
    rng = np.random.default_rng(7)
    n = 1000
    seen = []
    for _ in range(50):
        d = rng.normal(0.0, 1.0, n)
        d = d - d.mean()                              # ← this one line destroys it
        cross = (d.sum() ** 2 - np.sum(d**2)) / (n * (n - 1))
        seen.append(cross / np.mean(d**2))
    seen = np.array(seen)
    np.testing.assert_allclose(seen, -1.0 / (n - 1), rtol=1e-10)
    assert seen.std() < 1e-15, "it is an identity, so there must be no fluctuation"


# =============================================================================
# B8 — does the configurational temperature recover kT in a HOOMD run?  [slow]
# =============================================================================
@pytest.mark.slow
@pytest.mark.benchmark
def test_B8_configurational_temperature_recovers_kT_in_simulation(
        trap_sim_factory, positions_of):
    """★ BD's real thermometer. The kinetic temperature is drawn fresh every step,
    so it cannot be used.

    Harmonic trap: kT_conf = <|grad U|^2>/<lap U> = k^2<r^2>/(d k)
    """
    dt_star, k = 5e-3, 1.0
    sim, _ = trap_sim_factory(dt_star, seed=303)
    sim.run(int(20 / dt_star))

    stride = int(2.0 / dt_star)
    blocks = []
    for _ in range(40):
        sim.run(stride)
        p = positions_of(sim)
        blocks.append(configurational_temperature(k * p, laplacian_U_total=TRAP_DIM * k))

    assert_statistic_fluctuates(blocks, "kT_conf blocks")
    kT_conf = float(np.mean(blocks))
    se = se_of_mean(blocks)

    assert sigma_away(kT_conf, 1.0, se) < 4.0, (
        f"kT_conf = {kT_conf:.5f}±{se:.5f} departed from the input kT*=1.0")


@pytest.mark.slow
def test_kinetic_temperature_cannot_deviate_systematically(trap_sim_factory, hoomd_mod):
    """★ Pin as a test that this cannot be used as a guard.

    HOOMD `Brownian` does not integrate the velocity -- it **draws** it from the
    target distribution every step. So kinetic_temperature stays near kT even
    when dt* is raised 10x and the integration is a mess -- it cannot detect a
    systematic departure.
    """
    hoomd = hoomd_mod
    temps = {}
    for dt_star in (1e-3, 1e-1):          # 100x apart; the latter integrates badly
        sim, _ = trap_sim_factory(dt_star, seed=404)
        tq = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
        sim.operations.computes.append(tq)
        sim.run(200)
        temps[dt_star] = tq.kinetic_temperature

    # raise dt* 100x and the kinetic temperature stays near kT=1 → no diagnostic power
    for dt_star, T in temps.items():
        assert T == pytest.approx(1.0, rel=0.1), (
            f"at dt*={dt_star:g} the kinetic temperature is {T:.4f}")
