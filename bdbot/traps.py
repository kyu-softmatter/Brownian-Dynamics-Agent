"""A harmonic trap (`md.force.Custom`) -- fixed, constant-velocity and
oscillatory driving in one class.

**HOOMD has no harmonic trap** (`md.external.field` has only Electric, Magnetic
and Periodic, and `hpmc.external.Harmonic` is Monte-Carlo only). This generalizes
the verified snippet in skill `bd-hoomd`.

**Promoted after appearing three times**: `trap-2d-5um` (fixed),
`trap-drag-2d-hex300` (constant velocity) and `chain-bend-2d-oscill`
(oscillatory). All three differ only in `anchor(t)`.

    anchors(t) = anchor0 + velocity·t + drive(t)

Three of these are **silently-wrong traps**, so they are collected in one place:
  trap 1  minimum image -- without it, a weak trap is off by +1856% (with no error)
  trap 7  setting a non-periodic axis's period to inf gives
          `inf*round(0/inf) = nan`
  tags    `ParticleSorter` rearranges memory order -- a local snapshot is not in
          tag order
"""
from __future__ import annotations

import numpy as np


def _md():
    import hoomd.md as md
    return md


def make_trap(k, anchors, box, *, dt_star: float = 0.0, dims: int = 2,
              velocity=None, drive=None):
    """One harmonic trap. `k` is a scalar or a per-particle array (0 for
    untrapped particles).

    Arguments
      k         scalar or (N,) -- dimensionless stiffness [kT/d^2]
      anchors   (N,3) or (N,2) -- the anchors at t=0 (dimensionless length)
      box       scalar L or (Lx, Ly) -- the periodic lengths. 0/None makes that
                axis non-periodic
      dt_star   integration step (dimensionless). Needed only when the anchor
                moves -- t = timestep*dt_star
      velocity  (N,3) constant-velocity motion. Used by `trap-drag`
      drive     callable(t) -> (N,3) offset. Used by `chain-bend`'s a*sin(omega*t)

    * This function builds the class only when hoomd is imported -- so a front end
      like `bdbot.cli` does not pull in the heavy dependency (the convention in
      bdbot/__init__).

    WARNING: if the anchor is moved every U steps rather than every step, the
      drive becomes a zero-order hold: the fundamental shrinks by
      `sinc(omega*dt/2)` and lags by `omega*dt/2`. Measure the anchor position and
      use the measured phasor -- using the nominal amplitude put one K' off by
      236%, sign included. See `bdbot.lockin`.
    """
    md = _md()

    class HarmonicTrap(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=False)
            a = np.asarray(anchors, dtype=float)
            if a.shape[1] == 2:
                a = np.c_[a, np.zeros(len(a))]
            self.anchor0 = a
            n = len(a)
            self.k = (np.full(n, float(k)) if np.isscalar(k)
                      else np.asarray(k, dtype=float).reshape(n))
            b = (box, box) if np.isscalar(box) else tuple(box)
            # * trap 7: a non-periodic axis is left at 0 so the mask excludes it
            #   (never inf)
            self.period = np.array([float(b[0] or 0.0), float(b[1] or 0.0), 0.0])
            self.dt_star = float(dt_star)
            self.velocity = None if velocity is None else np.asarray(velocity, dtype=float)
            self.drive = drive
            self.dims = dims

        def centers(self, t: float) -> np.ndarray:
            c = self.anchor0
            if self.velocity is not None:
                c = c + self.velocity * t
            if self.drive is not None:
                c = c + np.asarray(self.drive(t), dtype=float)
            return c

        def _delta(self, pos, tags, t):
            d = pos - self.centers(t)[tags]
            m = self.period > 0                       # wrap periodic axes only (traps 1, 7)
            d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
            return d

        def set_forces(self, timestep):
            t = timestep * self.dt_star
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)     # * tag indexing is mandatory
                pos = np.array(snap.particles.position, copy=True)
                d = self._delta(pos, tags, t)
                kk = self.k[tags][:, None]
                arr.force[:] = -kk * d
                arr.potential_energy[:] = 0.5 * self.k[tags] * (d ** 2).sum(axis=1)

        def displacement(self, state, timestep) -> np.ndarray:
            """Displacement from the trap centre (N,3).

            The raw material for extracting drag force and phase.
            """
            snap = state.get_snapshot()
            pos = np.array(snap.particles.position, dtype=float)
            tags = np.arange(len(pos))       # get_snapshot is already in tag order
            return self._delta(pos, tags, timestep * self.dt_star)

    return HarmonicTrap()


__all__ = ["make_trap"]
