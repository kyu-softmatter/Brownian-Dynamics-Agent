"""External potentials — the ones HOOMD's built-ins cannot express.

Dimensionless convention: every parameter is taken in the reduced units (`*`) of the
relevant (system × target dynamics) card. The harmonic-trap card is based on
`l_trap`, so `k_star = 1` and `center_star` is in units of `l_trap`.
"""
from __future__ import annotations

import numpy as np
import hoomd


class HarmonicTrap(hoomd.md.force.Custom):
    """An isotropic harmonic trap:  U = 1/2 k |r - r0|^2,  F = -k (r - r0)

    HOOMD 7.1 has no built-in harmonic confinement potential (`md.external.field`
    offers only Periodic, Electric and Magnetic). A ghost particle plus
    `md.bond.Harmonic` would also work, but Custom is more direct and easier to
    inspect (user decision, 2026-07-28).

    Args:
        k_star:  the spring constant (reduced units). 1.0 under the trap card's
            convention
        center_star: the trap centre (in reduced length units)
        active_axes: the axes the trap acts on. (True, True, False) for a 2D system

    ⚠ **The box has to be much larger than the confinement length.** Wrapped
      coordinates are used here as they are -- a particle near the boundary gets
      pulled in the wrong direction. The harmonic-trap card's §7 gate
      "box >> l_trap" guarantees this (satisfied automatically since
      l_trap/sigma ~ 1e-3).
    """

    def __init__(self, k_star: float, center_star=(0.0, 0.0, 0.0),
                 active_axes=(True, True, True)):
        super().__init__()
        self.k_star = float(k_star)
        self._center = np.asarray(center_star, dtype=np.float64)
        self._mask = np.asarray(active_axes, dtype=bool)
        if self._center.shape != (3,) or self._mask.shape != (3,):
            raise ValueError("center_star and active_axes must have length 3")

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            # wrapped coordinates. See the warning above.
            dr = np.array(snap.particles.position, dtype=np.float64) - self._center
            dr[:, ~self._mask] = 0.0
            with self.cpu_local_force_arrays as arrays:
                arrays.force[:] = -self.k_star * dr
                # U = 1/2 k |dr|^2  (per particle)
                arrays.potential_energy[:] = 0.5 * self.k_star * np.sum(dr * dr, axis=1)
                # virial: W_ab = F_a * dr_b  (an external potential's contribution)
                # order convention: xx, xy, xz, yy, yz, zz
                #
                # ⚠ **UNVERIFIED.** `Force.virials` is not exposed on this path, so
                #   the test skips (tests/test_s5_forces.py). And the convention for
                #   an external field's virial is ambiguous -- a trap is not part of
                #   the system's momentum flux. Pressure is not a meaningful
                #   observable for a single confined particle, so this does not
                #   affect the pipeline, but **logging a pressure requires verifying
                #   this part first.**
                #   → registered as a knowledge/wiki/questions/ entry
                f = -self.k_star * dr
                arrays.virial[:, 0] = f[:, 0] * dr[:, 0]
                arrays.virial[:, 1] = f[:, 0] * dr[:, 1]
                arrays.virial[:, 2] = f[:, 0] * dr[:, 2]
                arrays.virial[:, 3] = f[:, 1] * dr[:, 1]
                arrays.virial[:, 4] = f[:, 1] * dr[:, 2]
                arrays.virial[:, 5] = f[:, 2] * dr[:, 2]


# =============================================================================
# Power-law repulsive pair potential — implemented with `md.pair.Table`
# =============================================================================
#  ★ Why Table
#  `U = A/r^n` is not among HOOMD's built-ins. The alternatives checked:
#    · `md.pair.Mie` — `U ∝ ε[(σ/r)^n − (σ/r)^m]`. Pure repulsion needs `m → 0`,
#      but the coefficient `(n/(n−m))(n/m)^{m/(n−m)}` diverges at `m=0`. Unusable
#    · `md.pair.Yukawa` — `exp(−κr)/r`, so only `r^{-1}`
#    · `md.pair.Table` — takes an arbitrary `U(r)` and `F(r)` as a table. **This is
#      the right one**
#  Table interpolates linearly, so the sample count sets the accuracy → a test
#  comparing against the analytic expression is mandatory (`tests/test_s5_pair.py`).
def power_law_table(hoomd, *, amplitude: float, exponent: float = 3.0,
                    r_cut: float, r_min: float = 0.2, n_points: int = 4096,
                    shift: bool = True, nlist=None, buffer: float = 0.1):
    """Builds the `U(r) = A/r^p` repulsion with `md.pair.Table` (reduced units).

    Args:
        amplitude: `A`. In reduced units `A = βU(r=1)`.
        exponent: `p`. Dipolar repulsion of paramagnetic colloids is `p = 3`.
        r_cut: the cutoff distance. **Must be `≤ L/2`** (minimum image). The caller
            guarantees it.
        r_min: where the table starts. **Come closer than this and HOOMD leaves the
            table** -- watch it with `guards.check_min_separation`.
        shift: subtract a constant so that `U(r_cut) = 0`. Removes the energy
            discontinuity.
            ⚠ **The force discontinuity remains** (subtracting a constant does not
            change the derivative). The size of the truncation error is computed by
            `pair_truncation_error()`.
        buffer: the cell-list buffer (**an absolute distance**, in `d`).
            ★ **HOOMD requires `r_cut + buffer ≤ L/2`** -- not `r_cut` alone.
            Exceed it and it is refused at runtime:
            `nlist: Simulation box is too small, the neighbor list is searching
            beyond the minimum image` (measured 2026-07-28).
            So the buffer must not be set as a fraction of `r_cut` -- that explodes
            at large `r_cut`.

    Returns:
        `(pair, info)` — `info` carries `u_shift` and the truncation error.
        Why `u_shift` is recorded: recovering the true total potential energy needs
        `U_true = U_hoomd + n_pairs_within_rcut × u_shift`.
    """
    if r_cut <= r_min:
        raise ValueError(f"r_cut({r_cut}) <= r_min({r_min})")
    nlist = nlist or hoomd.md.nlist.Cell(buffer=buffer)
    r = np.linspace(r_min, r_cut, n_points)
    u = amplitude / r**exponent
    f = exponent * amplitude / r**(exponent + 1)      # F = -dU/dr > 0 (repulsive)
    u_shift = float(amplitude / r_cut**exponent)
    if shift:
        u = u - u_shift

    pair = hoomd.md.pair.Table(nlist=nlist, default_r_cut=r_cut)
    pair.params[("A", "A")] = dict(r_min=r_min, U=u, F=f)
    info = {
        "potential": f"{amplitude:g}/r^{exponent:g}",
        "amplitude": amplitude, "exponent": exponent,
        "r_cut": r_cut, "r_min": r_min, "n_points": n_points,
        "shifted": shift, "u_shift": u_shift, "nlist_buffer": buffer,
        "r_cut_plus_buffer": r_cut + buffer,
        **pair_truncation_error(amplitude=amplitude, exponent=exponent, r_cut=r_cut),
    }
    return pair, info


def pair_truncation_error(*, amplitude: float, exponent: float, r_cut: float,
                          r_ref: float = 1.0) -> dict:
    """The error the truncation leaves. **Returns the absolute and relative values
    together.**

    ★ `simbot.cutoff`'s tolerances are on a `kT` basis (`βU ≤ 0.02`). That is right
      when `kT` is the governing scale. But in a crystal at `Γ = π^{3/2}A = 557`
      **`A·kT` is the governing scale**, and the same truncation is large on the
      `kT` basis and small on the `A` basis. Both are reported so a person chooses
      which basis to judge on.
    """
    u_cut = amplitude / r_cut**exponent
    f_cut = exponent * amplitude / r_cut**(exponent + 1)
    u_ref = amplitude / r_ref**exponent
    f_ref = exponent * amplitude / r_ref**(exponent + 1)
    return {
        "beta_u_at_rcut": u_cut,                 # kT basis (compare cutoff.py presets)
        "beta_f_at_rcut": f_cut,
        "u_rel_to_nearest": u_cut / u_ref,       # on the interaction scale
        "f_rel_to_nearest": f_cut / f_ref,
    }
