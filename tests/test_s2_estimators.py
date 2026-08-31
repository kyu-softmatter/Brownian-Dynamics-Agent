"""S2 — the prediction engine. Checks that the analytic solutions really are one.

The values here get sealed into `runs/*/02_prediction.md`, so **if the prediction
engine is wrong the whole verification is meaningless.** Hence the double check,
by identity and by limit.
"""
from __future__ import annotations

import math

import pytest

from simbot.estimators import (
    EFFICIENCY_BY_K, THROUGHPUT_PARTICLE_STEPS_PER_S, dt_star_for_trap_bias,
    estimate_wall_time_s, euler_maruyama_trap_variance_bias, harmonic_trap,
    overdamped_validity, samples_for_variance_precision, trap_run_length,
)
from simbot.units import kT_si, stokes_drag_si, stokes_einstein_D_si, water_viscosity_si

# the system of the first hand sketch — these values are sealed in 02_prediction.md
SKETCH = dict(T_si=300.0, eta_si=8.5566e-4, radius_si=5e-6, k_si=1e-5)


# =============================================================================
# harmonic-trap analytic solutions — checked by identity
# =============================================================================
@pytest.mark.parametrize("dim", [2, 3])
def test_equipartition_per_component_is_kT_over_k_regardless_of_dim(dim):
    """★ <x_i^2> = kT/k is dimension-independent. So it cannot discriminate dim."""
    p = harmonic_trap(dim=dim, **SKETCH)
    assert p.var_per_component_si == pytest.approx(p.kT_si / p.k_si, rel=1e-14)


@pytest.mark.parametrize("dim", [2, 3])
def test_radial_variance_scales_with_dim(dim):
    """<r^2> = d * kT/k. This is the dimension discriminator (ambiguity A1)."""
    p = harmonic_trap(dim=dim, **SKETCH)
    assert p.rms_radial_si**2 == pytest.approx(dim * p.kT_si / p.k_si, rel=1e-14)


def test_dimension_discriminator_ratio_is_exactly_three_halves():
    """★ The discriminator for 01_intake.md ambiguity A1:
    <r^2>(3D)/<r^2>(2D) = 3/2 exactly."""
    p2 = harmonic_trap(dim=2, **SKETCH)
    p3 = harmonic_trap(dim=3, **SKETCH)
    assert (p3.rms_radial_si**2) / (p2.rms_radial_si**2) == pytest.approx(1.5, rel=1e-14)
    assert p3.msd_plateau_si / p2.msd_plateau_si == pytest.approx(1.5, rel=1e-14)


@pytest.mark.parametrize("dim", [2, 3])
def test_msd_plateau_is_twice_radial_variance(dim):
    """MSD(t->inf) = <|r(inf)-r(0)|^2> = 2<r^2>  (the two positions are
    independent)."""
    p = harmonic_trap(dim=dim, **SKETCH)
    assert p.msd_plateau_si == pytest.approx(2 * dim * p.var_per_component_si, rel=1e-14)


def test_msd_reaches_plateau_at_long_times():
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.msd_si(50 * p.tau_trap_si) == pytest.approx(p.msd_plateau_si, rel=1e-9)


def test_msd_short_time_limit_recovers_free_diffusion():
    """★ A non-trivial limit check: for t << tau_trap, MSD -> 2 d D0 t.

    (2d kT/k)(1 - e^{-t/tau}) ~ (2d kT/k)(t k/gamma) = 2 d (kT/gamma) t = 2 d D0 t
    Confirms that the trap solution and Stokes-Einstein are mutually consistent.
    """
    p = harmonic_trap(dim=2, **SKETCH)
    t = 1e-4 * p.tau_trap_si
    free = 2 * p.dim * p.D0_si * t
    assert p.msd_si(t) == pytest.approx(free, rel=1e-3)


def test_relaxation_time_and_corner_frequency_are_consistent():
    """f_c = 1/(2 pi tau_trap). Matters because tweezer experiments measure f_c."""
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.tau_trap_si == pytest.approx(p.gamma_si / p.k_si, rel=1e-14)
    assert p.corner_freq_si == pytest.approx(1 / (2 * math.pi * p.tau_trap_si), rel=1e-14)


def test_confinement_length_is_sqrt_kT_over_k():
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.l_trap_si == pytest.approx(math.sqrt(p.kT_si / p.k_si), rel=1e-14)
    # <x^2> = l_trap^2 — check the two names really are the same thing
    assert p.l_trap_si**2 == pytest.approx(p.var_per_component_si, rel=1e-14)


def test_time_scale_separation_equals_k_star_sigma():
    """tau_D/tau_trap = k sigma^2/kT. The identity the card choice rests on."""
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.tau_sep == pytest.approx(p.k_star_sigma, rel=1e-12)
    assert p.l_trap_over_sigma == pytest.approx(1 / math.sqrt(p.k_star_sigma), rel=1e-12)


def test_estimator_warns_when_tau_D_would_be_wrong_time_unit():
    """Does the card's warning actually fire for a strong trap?"""
    p = harmonic_trap(dim=2, **SKETCH)
    joined = " ".join(p.notes)
    assert "tau_trap" in joined and "tau_D" in joined
    #  ★ Filter on the stable token, not the sentence. This read the Korean for
    #    "excluded volume" until 2026-08-30 and would have gone quiet the moment
    #    simbot/estimators.py was translated -- the same failure family as
    #    verify_intake_guards.py and verify_health.py. Caught here only because
    #    this file IS in the pytest suite.
    assert "excluded volume" in joined.lower()


# =============================================================================
# B6 — SI reference values. Pins the numbers sealed into 02_prediction.md
# =============================================================================
@pytest.mark.benchmark
def test_B6_sketch_reference_values():
    """Harmonic-trap card B6.  a=5um, k=10 pN/um, T=300K, water.

    Change these numbers and the sealed prediction in
    runs/2026-07-28_trap-2d-5um_2dfb9d/02_prediction.md becomes void.
    """
    p = harmonic_trap(dim=2, **SKETCH)
    assert p.kT_si * 1e21 == pytest.approx(4.14195, rel=1e-5)     # pN*nm
    assert p.gamma_si == pytest.approx(8.0644e-8, rel=1e-4)       # kg/s
    assert p.D0_si * 1e12 == pytest.approx(0.05136, rel=1e-3)     # um^2/s
    assert p.tau_trap_si * 1e3 == pytest.approx(8.064, rel=1e-3)  # ms
    assert p.tau_D_si == pytest.approx(1947.0, rel=1e-3)          # s
    assert p.l_trap_si * 1e9 == pytest.approx(20.35, rel=1e-3)    # nm
    assert p.corner_freq_si == pytest.approx(19.735, rel=1e-3)    # Hz
    assert p.k_star_sigma == pytest.approx(2.4143e5, rel=1e-3)
    assert p.var_per_component_si * 1e18 == pytest.approx(414.19, rel=1e-4)  # nm^2
    assert p.msd_plateau_si * 1e18 == pytest.approx(1656.78, rel=1e-4)       # nm^2


def test_sketch_values_reproduce_from_water_table():
    """Does SKETCH's eta match the value from the units module's water table?

    Stops two hard-coded copies of the same constant drifting apart.
    """
    eta, extrap = water_viscosity_si(300.0)
    assert extrap is False
    assert eta == pytest.approx(SKETCH["eta_si"], rel=1e-4)


# =============================================================================
# Euler-Maruyama systematic bias
# =============================================================================
@pytest.mark.parametrize("dt_star,expected_pct", [
    (2.0e-2, 1.0101), (1.0e-2, 0.5025), (5.0e-3, 0.2506),
    (2.5e-3, 0.1252), (1.0e-3, 0.0500),
])
def test_em_bias_known_values(dt_star, expected_pct):
    """Var* = 1/(1-dt*/2). These values are the baseline of the dt-ladder check."""
    assert euler_maruyama_trap_variance_bias(dt_star) * 100 == pytest.approx(
        expected_pct, rel=2e-3)


def test_em_bias_is_first_order_in_dt():
    """★ Pin that it is first order: halve dt and the bias (nearly) halves."""
    b1 = euler_maruyama_trap_variance_bias(1e-2)
    b2 = euler_maruyama_trap_variance_bias(5e-3)
    assert b1 / b2 == pytest.approx(2.0, rel=5e-3)
    # a second-order scheme would give a ratio of 4 — reject that hypothesis
    assert abs(b1 / b2 - 4.0) > 1.0


def test_em_bias_leading_order_matches_dt_over_two():
    for dt in (1e-3, 1e-4):
        assert euler_maruyama_trap_variance_bias(dt) == pytest.approx(dt / 2, rel=1e-2)


def test_em_bias_rejects_unstable_dt():
    """At dt* >= 2 the deterministic term (1-dt*) has |.|>1 and it diverges."""
    for bad in (0.0, -1e-3, 2.0, 3.0):
        with pytest.raises(ValueError):
            euler_maruyama_trap_variance_bias(bad)


def test_dt_star_for_bias_is_inverse_of_bias():
    for target in (1e-4, 1e-3, 2.5e-3, 1e-2):
        dt = dt_star_for_trap_bias(target)
        assert euler_maruyama_trap_variance_bias(dt) == pytest.approx(target, rel=1e-12)


# =============================================================================
# premise checks (overdamped · Stokes)
# =============================================================================
def test_sketch_system_passes_overdamped_and_stokes():
    p = harmonic_trap(dim=2, **SKETCH)
    V = 4 / 3 * math.pi * SKETCH["radius_si"] ** 3
    c = overdamped_validity(
        gamma_si=p.gamma_si, mass_si=1050 * V, tau_process_si=p.tau_trap_si,
        velocity_scale_si=p.l_trap_si / p.tau_trap_si,
        radius_si=SKETCH["radius_si"], eta_si=SKETCH["eta_si"], rho_fluid_si=996.5)
    assert c.passed, c.failures
    assert c.inertial_ratio == pytest.approx(8.454e-4, rel=1e-2)
    assert c.reynolds == pytest.approx(1.470e-5, rel=1e-2)


def test_overdamped_check_fails_for_heavy_particle():
    """Does the guard actually fire? Test only the passing case and you will not
    notice when it dies."""
    c = overdamped_validity(
        gamma_si=1e-8, mass_si=1.0, tau_process_si=1e-3,
        velocity_scale_si=1.0, radius_si=1e-6, eta_si=1e-3, rho_fluid_si=1000.0)
    assert not c.passed
    assert any("inertial" in f for f in c.failures)
    assert any("Reynolds" in f for f in c.failures)


# =============================================================================
# statistical precision · cost model
# =============================================================================
def test_samples_for_precision_follows_inverse_square():
    """SE(s^2)/s^2 = sqrt(2/n)  =>  n = 2/rel^2"""
    assert samples_for_variance_precision(0.01) == pytest.approx(20000)
    assert samples_for_variance_precision(0.005) == pytest.approx(80000)
    # double the precision and the sample count goes up 4x
    assert (samples_for_variance_precision(0.005)
            / samples_for_variance_precision(0.01)) == pytest.approx(4.0)


def test_trap_run_length_matches_prediction_document():
    """Pins the cost table in 02_prediction.md §6."""
    p = harmonic_trap(dim=2, **SKETCH)
    r = trap_run_length(n_particles=1000, rel_err_target=0.01,
                        tau_trap_si=p.tau_trap_si)
    assert r["independent_samples_needed"] == pytest.approx(20000)
    assert r["t_total_in_tau_trap"] == pytest.approx(40.0)
    steps = r["t_total_in_tau_trap"] / 5e-3
    assert steps == pytest.approx(8000)
    assert estimate_wall_time_s(1000, int(steps), 1) == pytest.approx(1.27, rel=0.05)


def test_wall_time_scales_linearly_in_work():
    assert estimate_wall_time_s(1000, 1000) == pytest.approx(
        estimate_wall_time_s(500, 2000), rel=1e-12)
    assert estimate_wall_time_s(1000, 2000) == pytest.approx(
        2 * estimate_wall_time_s(1000, 1000), rel=1e-12)


def test_wall_time_uses_measured_efficiency_and_penalizes_concurrency():
    """Does it use the measured efficiency table? k=8 must be slower per process
    than k=1."""
    w1 = estimate_wall_time_s(1000, 10000, 1)
    w8 = estimate_wall_time_s(1000, 10000, 8)
    assert w8 > w1
    assert w8 / w1 == pytest.approx(1 / EFFICIENCY_BY_K[8], rel=1e-12)


def test_wall_time_falls_back_for_unmeasured_concurrency():
    """A k not in the table is approximated conservatively by the value below it."""
    assert estimate_wall_time_s(100, 100, 7) == pytest.approx(
        estimate_wall_time_s(100, 100, 6), rel=1e-12)


def test_throughput_constant_matches_measurement():
    """The measured constants from knowledge/wiki/findings/local-cpu-parallelism.md."""
    assert THROUGHPUT_PARTICLE_STEPS_PER_S == pytest.approx(6.3e6, rel=1e-9)
    assert EFFICIENCY_BY_K[4] == pytest.approx(0.925, abs=1e-3)
    assert EFFICIENCY_BY_K[12] < EFFICIENCY_BY_K[10], \
        "k=12 is a regression — a measured fact"
