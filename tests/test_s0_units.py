"""S0 — units, constants and non-dimensionalization scales. The base of every stage.

The bug class this file is out to catch: **the silently wrong kind.**
Radius/diameter confusion, kT confused with eta, round-trip conversion loss --
none of them diverge, so without a test you never find out.
"""
from __future__ import annotations

import math

import pytest

from simbot.units import (
    K_B, Scales, kT_si, scales_brownian, scales_harmonic_trap,
    stokes_drag_si, stokes_einstein_D_si, water_density_si, water_viscosity_si,
)


# --- constants --------------------------------------------------------------
def test_boltzmann_is_si_2019_exact_value():
    assert K_B == 1.380649e-23


@pytest.mark.parametrize("T_si,expected_pN_nm", [
    (293.15, 4.047373),   # 20 C  — corrected 2026-07-28: had wrongly read 4.047872
    (298.15, 4.116405),   # 25 C
    (300.00, 4.141947),   # 300 K — the first hand sketch
])
def test_kT_in_pN_nm(T_si, expected_pN_nm):
    """1 pN*nm = 1e-21 J. The optical-tweezer literature uses this unit."""
    assert kT_si(T_si) * 1e21 == pytest.approx(expected_pN_nm, rel=1e-6)


# --- water properties -------------------------------------------------------
def test_water_viscosity_at_table_points_is_exact():
    """A point that is in the table is not interpolated."""
    eta, extrap = water_viscosity_si(298.15)
    assert eta == 0.8900e-3
    assert extrap is False


def test_water_viscosity_20C_and_25C_differ_by_11_percent():
    """★ The most common confusion: 1.002 mPa*s is the 20 C value; 25 C is
    0.890 mPa*s.

    That difference (11 %) propagates straight into D0 and every timescale.
    """
    eta20, _ = water_viscosity_si(293.15)
    eta25, _ = water_viscosity_si(298.15)
    assert eta20 == pytest.approx(1.0016e-3)
    assert eta25 == pytest.approx(0.8900e-3)
    assert (eta20 - eta25) / eta25 == pytest.approx(0.1254, rel=0.01)


def test_water_viscosity_interpolates_monotonically():
    ts = [293.15, 295.0, 298.15, 300.0, 303.15, 306.0, 308.15]
    etas = [water_viscosity_si(t)[0] for t in ts]
    assert all(a > b for a, b in zip(etas, etas[1:])), \
        "viscosity decreases monotonically with temperature"


def test_water_viscosity_at_300K_matches_recorded_value():
    """T = 300 K, from the first hand sketch. Must match the sealed prediction."""
    eta, extrap = water_viscosity_si(300.0)
    assert extrap is False
    assert eta == pytest.approx(8.5566e-4, rel=1e-4)


def test_water_viscosity_flags_extrapolation():
    """Outside the interpolation range (293-308 K) it must say so -- in order to
    lower the provenance."""
    _, extrap_lo = water_viscosity_si(280.0)
    _, extrap_hi = water_viscosity_si(330.0)
    _, extrap_in = water_viscosity_si(300.0)
    assert extrap_lo is True and extrap_hi is True and extrap_in is False


def test_water_density_is_near_1000():
    rho, _ = water_density_si(298.15)
    assert 990 < rho < 1000


# --- Stokes drag: radius vs diameter ----------------------------------------
def test_stokes_drag_uses_radius_not_diameter():
    """★ Pin the bug `master_plan` §S3 names as "the most common mistake".

    gamma = 6*pi*eta*a  (a = the radius).  Pass the diameter and gamma doubles,
    D0 halves, and **every timescale is wrong by 2x.** It does not diverge.
    """
    eta, a = 1e-3, 5e-6
    gamma = stokes_drag_si(eta, a)
    assert gamma == pytest.approx(6 * math.pi * eta * a)
    # the diameter gives exactly 2x — nail that relation down as a test
    assert stokes_drag_si(eta, 2 * a) == pytest.approx(2 * gamma)


def test_stokes_einstein_reference_case_1um_sphere_in_water_25C():
    """The literature's conventional value: a 1 um-diameter sphere in water at
    25 C has D ~ 0.49 um^2/s.

    The reference case in knowledge/wiki/concepts/water-298k.md.
    """
    T, a = 298.15, 0.5e-6
    eta, _ = water_viscosity_si(T)
    gamma = stokes_drag_si(eta, a)
    D0 = stokes_einstein_D_si(T, gamma)
    assert D0 * 1e12 == pytest.approx(0.4909, rel=2e-3)   # um^2/s
    tau_B = (2 * a) ** 2 / D0
    assert tau_B == pytest.approx(2.037, rel=2e-3)        # s


# --- Scales round-trip ------------------------------------------------------
ALL_KINDS = ["length", "energy", "time", "force", "stiffness", "velocity",
             "diffusivity", "rate", "modulus_3d", "area", "volume"]


@pytest.mark.parametrize("kind", ALL_KINDS)
def test_scales_roundtrip_is_lossless(kind):
    """The SI -> star -> SI round-trip error must be within the floating-point
    limit.

    `master_plan` §S4 gate: relative error < 1e-12
    """
    s = Scales(length_si=2.0352e-8, energy_si=4.1419e-21, time_si=8.0644e-3)
    value_si = 3.14159e-7
    back = s.to_si(s.to_star(value_si, kind), kind)
    assert back == pytest.approx(value_si, rel=1e-14)


def test_scales_rejects_unknown_kind():
    s = Scales(length_si=1.0, energy_si=1.0, time_si=1.0)
    with pytest.raises(KeyError, match="unknown scale kind"):
        s.to_star(1.0, "wibble")


def test_derived_scales_are_dimensionally_consistent():
    s = Scales(length_si=3.0, energy_si=7.0, time_si=11.0)
    assert s.force_si == pytest.approx(7.0 / 3.0)
    assert s.stiffness_si == pytest.approx(7.0 / 9.0)
    assert s.velocity_si == pytest.approx(3.0 / 11.0)
    assert s.diffusivity_si == pytest.approx(9.0 / 11.0)
    assert s.rate_si == pytest.approx(1.0 / 11.0)
    assert s.modulus_3d_si == pytest.approx(7.0 / 27.0)


# --- definitional invariants of each card's scales --------------------------
def test_harmonic_trap_scales_normalize_D_and_k_to_exactly_one():
    """★ The central claim of harmonic-trap card §3.

    Choose (l_trap, kT, tau_trap) and the dimensionless equation of motion
    normalizes, with no parameters left, to
    `dr*/dt* = -r* + sqrt(2) xi`  =>  D* = 1, k* = 1, <x*^2> = 1.
    """
    T, eta, a, k = 300.0, 8.5566e-4, 5e-6, 1e-5
    gamma = stokes_drag_si(eta, a)
    s = scales_harmonic_trap(k_si=k, T_si=T, gamma_si=gamma)

    D0 = stokes_einstein_D_si(T, gamma)
    assert s.to_star(D0, "diffusivity") == pytest.approx(1.0, rel=1e-14)
    assert s.to_star(k, "stiffness") == pytest.approx(1.0, rel=1e-14)
    # equipartition <x^2> = kT/k is also exactly 1 in reduced units
    assert s.to_star(kT_si(T) / k, "area") == pytest.approx(1.0, rel=1e-14)


def test_brownian_scales_normalize_D_and_tauD_to_exactly_one():
    """Passive sphere x transport card: choose (sigma, kT, tau_D) and D* = 1,
    tau_D* = 1."""
    T, eta, a = 298.15, 0.8900e-3, 0.5e-6
    gamma = stokes_drag_si(eta, a)
    sigma = 2 * a
    s = scales_brownian(sigma_si=sigma, T_si=T, gamma_si=gamma)

    D0 = stokes_einstein_D_si(T, gamma)
    assert s.to_star(D0, "diffusivity") == pytest.approx(1.0, rel=1e-14)
    assert s.to_star(sigma**2 / D0, "time") == pytest.approx(1.0, rel=1e-14)


def test_two_cards_give_time_scales_separated_by_k_star_sigma():
    """★ Pin, numerically, why the card system is needed at all.

    tau_D / tau_trap = k sigma^2 / kT = k*_sigma
    In the first hand sketch this is 2.41e5  =>  a dt chosen from tau_D comes out
    12x the relaxation time.
    """
    T, eta, a, k = 300.0, 8.5566e-4, 5e-6, 1e-5
    gamma = stokes_drag_si(eta, a)
    sigma = 2 * a
    s_bd = scales_brownian(sigma_si=sigma, T_si=T, gamma_si=gamma)
    s_tr = scales_harmonic_trap(k_si=k, T_si=T, gamma_si=gamma)

    ratio = s_bd.time_si / s_tr.time_si
    k_star_sigma = k * sigma**2 / kT_si(T)
    assert ratio == pytest.approx(k_star_sigma, rel=1e-12)
    assert ratio == pytest.approx(2.414e5, rel=1e-3)

    # if dt* = 5e-5 is taken from tau_D, how many relaxation times is that
    dt_in_tau_trap = 5e-5 * ratio
    assert dt_in_tau_trap == pytest.approx(12.07, rel=1e-2)
    assert dt_in_tau_trap > 1.0, \
        "pin the fact that the universal convention exceeds this system's " \
        "relaxation time"
