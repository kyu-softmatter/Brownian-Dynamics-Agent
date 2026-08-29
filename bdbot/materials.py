"""Derived material properties -- all dimensional (skill `bd-physics` section 2).

Only what the first two cases computed identically: gamma, D_t, tau_B, m, tau_p.

WARNING: do not mix what holds only for a sphere with what is shape-independent.
   Anything named `sphere_*` is sphere-only. Ellipsoids and rods need Perrin
   factors, and **translational anisotropy is not reproducible in BD at all**
   (see the hard constraints in skill `bd-hoomd`).
"""
from __future__ import annotations

import math

from .units import Q, kB


def thermal_energy(T):
    return (kB * T).to("J")


def sphere_drag(eta, d):
    """Stokes drag coefficient gamma = 3*pi*eta*d (sphere)."""
    return (3 * math.pi * eta * d).to("kg/s")


def diffusion(kT, gamma):
    """Stokes-Einstein D_t = kT/gamma. Shape-independent (only gamma has to be right)."""
    return (kT / gamma).to("m^2/s")


def brownian_time(d, D_t):
    """Diffusion time tau_B = d^2/D_t. The value used as the reference time."""
    return (d**2 / D_t).to("s")


def sphere_mass(rho_p, d):
    return (rho_p * (math.pi / 6) * d**3).to("kg")


def momentum_time(m, gamma):
    """Momentum relaxation tau_p = m/gamma.

    * Do not compare this against `dt`. BD has no inertia at all (skill
      `bd-physics` section 4). The question `tau_p` answers is "is BD an
      admissible model for this system", and the thing to compare against is that
      system's **fastest timescale of interest**, `tau_dyn`.
    """
    return (m / gamma).to("s")


def sphere_rotational_diffusion(kT, eta, d):
    """Stokes-Einstein-Debye D_r = kT/(pi*eta*d^3) = 3D_t/d^2.  Sphere only."""
    return (kT / (math.pi * eta * d**3)).to("1/s")


def sphere_bulk(d, T, eta, rho_p=None) -> dict:
    """The basic properties of a sphere in a Newtonian fluid, in one call.

    Both of the first two cases used exactly this bundle.
    """
    d = d.to("m")
    kT = thermal_energy(T)
    gamma = sphere_drag(eta.to("Pa*s"), d)
    D_t = diffusion(kT, gamma)
    out = {"kT": kT, "gamma": gamma, "D_t": D_t, "tau_B": brownian_time(d, D_t), "d": d}
    if rho_p is not None:
        m = sphere_mass(rho_p.to("kg/m^3"), d)
        out["m"] = m
        out["tau_p"] = momentum_time(m, gamma)
    return out


# -- material constants (handbook, tier 0) ------------------------------------
WATER_VISCOSITY = {300: Q(0.851, "mPa*s"), 298: Q(0.890, "mPa*s"), 293: Q(1.002, "mPa*s")}
DENSITY = {"silica": Q(2000, "kg/m^3"), "polystyrene": Q(1050, "kg/m^3")}

__all__ = ["thermal_energy", "sphere_drag", "diffusion", "brownian_time", "sphere_mass",
           "momentum_time", "sphere_rotational_diffusion", "sphere_bulk",
           "WATER_VISCOSITY", "DENSITY"]
