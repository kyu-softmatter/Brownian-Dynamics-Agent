"""Derived material properties -- all dimensional (skill `bd-physics` section 2).

Only what the first two cases computed identically: gamma, D_t, tau_B, m, tau_p.

WARNING: do not mix what holds only for a sphere with what is shape-independent.
   Anything named `sphere_*` is sphere-only. Ellipsoids and rods need Perrin
   factors, and **translational anisotropy is not reproducible in BD at all**
   (see the hard constraints in skill `bd-hoomd`).

This is the **pint-carrying** face of the medium properties. The float face is
`simbot.units` (`*_si`), and the numbers both of them read come from
[`bdbot.constants`](constants.py) -- one table, not two. Before that merge the two
copies had already drifted 0.545 % at 300 K.

**The formulas are deliberately NOT routed through a shared float kernel.**
`gamma = 3*pi*eta*d` here and `gamma = 6*pi*eta*a` in `simbot.units` are the same
relation in two conventions, and collapsing them would mean stripping and
re-attaching pint units on a quantity that feeds `run_id` hashes and a 1e-12
round-trip check -- a last-bit change for no gain. What kills the drift instead is
`tests/test_cross_package_equivalence.py`, which asserts the two agree exactly.
⚠ `sphere_drag` takes a **diameter**, `simbot.units.stokes_drag_si` takes a
  **radius**. That is the factor-2 trap the equivalence test exists to pin down.
"""
from __future__ import annotations

import math

from . import constants as _const
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


# -- medium properties: a pint view over the shared table ---------------------
def water_viscosity(T, *, strict: bool = True):
    """Water viscosity at `T` [pint]. Interpolated from the IAPWS table in
    `bdbot.constants`, which `simbot.units` reads too.

    `strict=True` refuses to extrapolate outside 293-308 K rather than returning a
    number whose provenance has quietly dropped from `derived` to `assumed`
    (CLAUDE.md rule 3). Pass `strict=False` to accept it deliberately.
    """
    eta, extrapolated = _const.water_viscosity_si(T.to("K").magnitude)
    if extrapolated and strict:
        raise ValueError(
            f"T = {T} is outside the water table {_const.WATER_TABLE_RANGE_K} K, so "
            f"this is an extrapolation. Pass strict=False to accept it, and drop the "
            f"provenance tier when you do.")
    return Q(eta, "Pa*s")


def water_density(T, *, strict: bool = True):
    """Water density at `T` [pint]. Same table, same extrapolation contract."""
    rho, extrapolated = _const.water_density_si(T.to("K").magnitude)
    if extrapolated and strict:
        raise ValueError(f"T = {T} is outside {_const.WATER_TABLE_RANGE_K} K "
                         f"(density table); pass strict=False to accept it.")
    return Q(rho, "kg/m^3")


# The three values the first cases were written against, kept as a dict because
# they are cited that way in the case scripts' provenance notes.
#  ★ `300 K -> 0.851 mPa*s` is NOT the IAPWS table value (that interpolates to
#    0.85566). It is the separately-sourced anchor in
#    `knowledge/source/books/welty_transport.md` section 1.1, and it is what the 8
#    `bdbot` cases actually ran -- it is written into their spec files and therefore
#    into their `run_id` hashes. `simbot`'s sealed S2 documents used the table
#    value. **Both are in use; the 0.545 % gap is reported, not resolved**
#    (`bdbot.constants.water_viscosity_provenance_gap`).
WATER_VISCOSITY = {
    300: Q(_const.WATER_ETA_SOURCED_SI[300.00], "Pa*s"),     # welty section 1.1
    298: Q(_const.WATER_ETA_SI[298.15], "Pa*s"),             # IAPWS 25 C
    293: Q(_const.WATER_ETA_SI[293.15], "Pa*s"),             # IAPWS 20 C
}
DENSITY = {"silica": Q(2000, "kg/m^3"), "polystyrene": Q(1050, "kg/m^3")}

__all__ = ["thermal_energy", "sphere_drag", "diffusion", "brownian_time", "sphere_mass",
           "momentum_time", "sphere_rotational_diffusion", "sphere_bulk",
           "water_viscosity", "water_density", "WATER_VISCOSITY", "DENSITY"]
