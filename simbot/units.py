"""Physical constants and the SI ↔ dimensionless conversion.

Convention (CLAUDE.md):
  - `*_si`   : an SI quantity (a float; the unit is pinned by the name and docs)
  - `*_star` : a dimensionless value (a bare float)
  Arithmetic mixing the two suffixes is a bug. `tests/test_units.py` watches for it.

The reference scales for non-dimensionalization are owned by the
**(system × target dynamics) card**. So this module does not impose a universal
convention -- it takes `Scales` as an argument.
  → knowledge/wiki/systems/_index.md
"""
from __future__ import annotations

from dataclasses import dataclass
import math

# --- constants and water properties: the single source of truth is `bdbot.constants` ---
#  ★ Merged 2026-08-29. This table also existed in `bdbot.materials.WATER_VISCOSITY`,
#    and the two copies **had already diverged** -- 0.851 vs 0.85566 mPa*s at 300 K,
#    a 0.545 % gap. 0.545 % in `eta` is 0.545 % in `gamma`, in `D`, and in **every
#    timescale derived from them.**
#    Values, interpolation and the extrapolation flag are unchanged and must stay
#    bit-identical: the sealed S2 documents contain 8.5566e-4. The only thing that
#    changed is that the definition now exists once.
#  ⚠ Importing `bdbot` is not expensive here: `bdbot/__init__` is lazy (PEP 562), so
#    `bdbot.constants` pulls in neither pint nor numpy (0.01 s, measured).
#    **The dependency runs one way only, `simbot -> bdbot`.** Reversed, `bdbot.cli`
#    would start pulling in matplotlib.
from bdbot.constants import (K_B, WATER_ETA_SI as _WATER_ETA_TABLE_SI,
                             WATER_RHO_SI as _WATER_RHO_TABLE_SI,
                             interp_table as _interp,
                             water_density_si, water_viscosity_si)


# --- basic relations --------------------------------------------------------
def kT_si(T_si: float) -> float:
    """Thermal energy [J]."""
    return K_B * T_si


def stokes_drag_si(eta_si: float, radius_si: float) -> float:
    """Stokes drag coefficient of a sphere, gamma = 6*pi*eta*a  [kg/s].

    ⚠ The argument is the **radius**. Pass the diameter and gamma doubles, D
       halves, and every timescale is wrong by 2x. See
       knowledge/wiki/concepts/water-298k.md.
    """
    return 6.0 * math.pi * eta_si * radius_si


def stokes_einstein_D_si(T_si: float, gamma_si: float) -> float:
    """Translational diffusion coefficient D = kT/gamma  [m^2/s]."""
    return kT_si(T_si) / gamma_si


# --- reference scales for non-dimensionalization ----------------------------
@dataclass(frozen=True)
class Scales:
    """The 3 non-dimensionalization references. Which values to pick is decided by
    the (system × target dynamics) card.

    length_si : the reference length [m]
    energy_si : the reference energy [J]
    time_si   : the reference time [s]
    origin    : where this choice came from (a card path, say)
    """

    length_si: float
    energy_si: float
    time_si: float
    origin: str = ""

    # derived scales
    @property
    def force_si(self) -> float:
        return self.energy_si / self.length_si

    @property
    def stiffness_si(self) -> float:
        """Spring-constant scale [N/m]."""
        return self.energy_si / self.length_si**2

    @property
    def velocity_si(self) -> float:
        return self.length_si / self.time_si

    @property
    def diffusivity_si(self) -> float:
        return self.length_si**2 / self.time_si

    @property
    def rate_si(self) -> float:
        """Angular frequency, rotational diffusion and the like [1/s]."""
        return 1.0 / self.time_si

    @property
    def modulus_3d_si(self) -> float:
        """Modulus scale [Pa] = energy/length^3."""
        return self.energy_si / self.length_si**3

    # --- conversion ---
    def to_star(self, value_si: float, kind: str) -> float:
        return value_si / self._scale_for(kind)

    def to_si(self, value_star: float, kind: str) -> float:
        return value_star * self._scale_for(kind)

    def _scale_for(self, kind: str) -> float:
        table = {
            "length": self.length_si,
            "energy": self.energy_si,
            "time": self.time_si,
            "force": self.force_si,
            "stiffness": self.stiffness_si,
            "velocity": self.velocity_si,
            "diffusivity": self.diffusivity_si,
            "rate": self.rate_si,
            "modulus_3d": self.modulus_3d_si,
            "area": self.length_si**2,
            "volume": self.length_si**3,
        }
        if kind not in table:
            raise KeyError(f"unknown scale kind {kind!r}; known: {sorted(table)}")
        return table[kind]


# --- per-card reference-scale factories -------------------------------------
def scales_brownian(sigma_si: float, T_si: float, gamma_si: float) -> Scales:
    """Passive sphere × transport:  (sigma, kT, tau_D = sigma^2/D0).

    Card: knowledge/wiki/systems/passive-sphere--*.md
    """
    D0 = stokes_einstein_D_si(T_si, gamma_si)
    return Scales(
        length_si=sigma_si,
        energy_si=kT_si(T_si),
        time_si=sigma_si**2 / D0,
        origin="scales_brownian: (sigma, kT, tau_D)",
    )


def scales_soft2d(d_si: float, sigma_si: float, T_si: float,
                  gamma_si: float | None = None) -> Scales:
    """2D soft-repulsive × equilibrium structure:
    (`d = n^{-1/2}`, `kT`, `τ_d = d²/D₀`).

    Card: knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md §3

    ★ **The only card whose length scale differs from the particle size.** The
      length unit is the lattice spacing `d = n^{-1/2}`, while the drag is set by
      the particle diameter `σ`: `γ = 6πη(σ/2) = 3πησ`. Treat the two as equal and
      `τ_d` is wrong by `(d/σ)²` -- **9x** at this sweep's `d/σ = 3`.

    ⚠ `σ` **does not enter the dynamics** (there is no hard core). It is used only
      for the timescale and for judging physical validity →
      `build.coverage_from_sigma_over_d`.
    """
    if gamma_si is None:
        eta_si, extrapolated = water_viscosity_si(T_si)
        if extrapolated:
            raise ValueError(
                f"T = {T_si} K is outside the water-viscosity table (293–308 K), "
                f"so it would be extrapolated. State gamma_si explicitly -- using "
                f"an extrapolated value quietly turns an 'assumed' provenance into "
                f"one pretending to be 'derived'")
        gamma_si = stokes_drag_si(eta_si, sigma_si / 2.0)
    D0 = stokes_einstein_D_si(T_si, gamma_si)
    return Scales(
        length_si=d_si,
        energy_si=kT_si(T_si),
        time_si=d_si**2 / D0,
        origin=f"scales_soft2d: (d={d_si:.4g} m, kT, tau_d=d^2/D0), "
               f"gamma from sigma={sigma_si:.4g} m",
    )


def scales_harmonic_trap(k_si: float, T_si: float, gamma_si: float) -> Scales:
    """Harmonic trap:  (l_trap = sqrt(kT/k), kT, tau_trap = gamma/k).

    Under this choice the dimensionless equation of motion normalizes to
    `dr*/dt* = -r* + sqrt(2) xi`, giving `D* = 1` and `k* = 1`.

    ★ Why tau_D must not be imposed: tau_trap/tau_D = kT/(k sigma^2) =
      1/k*_sigma. In a strong trap (k*_sigma >> 1) a dt taken from tau_D is larger
      than the relaxation time.
      Card: knowledge/wiki/systems/passive-sphere--harmonic-trap.md
    """
    return Scales(
        length_si=math.sqrt(kT_si(T_si) / k_si),
        energy_si=kT_si(T_si),
        time_si=gamma_si / k_si,
        origin="scales_harmonic_trap: (l_trap, kT, tau_trap)",
    )
