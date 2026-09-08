"""The HOOMD BD execution skeleton -- only the parts both of the first two cases
used identically.

Common: building a 2D frame, Simulation, Brownian + Integrator, the GSD Tier A
writer, seed handling.
**What differs per case**: the forces, the initial placement, the sampling loop ->
those stay in the case.

Every convention here was verified by measurement in skill `bd-hoomd`:
  trap 3   `integrate_rotational_dof = False`
  trap 5   BD is overdamped -- `thermalize_particle_momenta()` is unnecessary and
           velocity has no meaning
  trap 9   2D needs `Lz=0` **and** `dimensions=2`
  trap 12  the seed is truncated to 16 bits
"""
from __future__ import annotations

import numpy as np

HOOMD_SEED_MAX = 65535


def resolve_seed(seed: int) -> tuple[int, int]:
    """(numpy seed, HOOMD seed). HOOMD truncates to 16 bits, so truncate up front.

    * trap 12: passing `seed=20260803` makes HOOMD truncate it to `10179` with a
      warning. Two seeds differing by 65536 produce **the same trajectory.** numpy
      (used for the initial placement) takes the full seed, so the two are passed
      separately.
    """
    return int(seed), int(seed) & HOOMD_SEED_MAX


def frame_2d(positions, L, types=("A",), typeid=None, orientation=False,
             bonds=None, angles=None):
    """A 2D periodic frame. `positions` is (N,2) or (N,3).

    `L` is a scalar (square) or `(Lx, Ly)` -- * this became necessary when
    `trap-drag`'s commensurate hexagonal lattice required a rectangular box
    (commensurability sets the aspect ratio, it is not chosen).

    `bonds`/`angles` are (M,2)/(M,3) index arrays -- used by the chain cases.
    """
    import gsd.hoomd

    p = np.asarray(positions, dtype=float)
    if p.shape[1] == 2:
        p = np.c_[p, np.zeros(len(p))]
    n = len(p)
    Lx, Ly = (float(L), float(L)) if np.isscalar(L) else (float(L[0]), float(L[1]))
    fr = gsd.hoomd.Frame()
    fr.particles.N = n
    fr.particles.position = p
    fr.particles.typeid = [0] * n if typeid is None else list(typeid)
    fr.particles.types = list(types)
    if orientation:
        fr.particles.orientation = [(1.0, 0.0, 0.0, 0.0)] * n
    if bonds is not None:
        b = np.asarray(bonds, dtype=int)
        fr.bonds.N = len(b)
        fr.bonds.types = ["backbone"]
        fr.bonds.typeid = [0] * len(b)
        fr.bonds.group = b
    if angles is not None:
        a = np.asarray(angles, dtype=int)
        fr.angles.N = len(a)
        fr.angles.types = ["bend"]
        fr.angles.typeid = [0] * len(a)
        fr.angles.group = a
    fr.configuration.box = [Lx, Ly, 0, 0, 0, 0]   # Lz=0 -> 2D (trap 9)
    fr.configuration.dimensions = 2
    return fr


def make_sim(frame, seed: int, notice_level: int = 0):
    import hoomd

    _, hseed = resolve_seed(seed)
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=notice_level), seed=hseed)
    sim.create_state_from_snapshot(frame)
    return sim


def frame_3d(positions, L, types=("A",), typeid=None, bonds=None, angles=None):
    """A 3D periodic frame. `positions` is (N,3), `L` a scalar or `(Lx,Ly,Lz)`.

    * Promoted 2026-09-02 on its second appearance: `cases/network_3d.py:754`
      inlines this, and `bead-water-3d` needs it. Promotion rule unchanged --
      *"has it appeared twice?"*

    ⚠ The 2D/3D difference is one number and it is trap 9 in skill `bd-hoomd`:
      2D sets `box=[Lx,Ly,0,...]` **and** `dimensions=2`; 3D sets `Lz=L` and
      `dimensions=3`. Setting `Lz=0` with `dimensions=3` gives a zero-volume box.
    """
    import gsd.hoomd

    p = np.asarray(positions, dtype=float)
    if p.ndim != 2 or p.shape[1] != 3:
        raise ValueError(f"3D positions must be (N,3), got {p.shape}")
    Ls = ([float(L)] * 3 if np.isscalar(L) else [float(x) for x in L])
    if len(Ls) != 3 or any(x <= 0 for x in Ls):
        raise ValueError(f"L must give three positive lengths, got {L}")

    f = gsd.hoomd.Frame()
    f.particles.N = int(p.shape[0])
    f.particles.position = p
    f.particles.types = list(types)
    f.particles.typeid = (np.zeros(p.shape[0], dtype=int) if typeid is None
                          else np.asarray(typeid, dtype=int))
    f.particles.mass = np.ones(p.shape[0])
    f.configuration.box = [Ls[0], Ls[1], Ls[2], 0, 0, 0]
    f.configuration.dimensions = 3
    if bonds is not None:
        b = np.asarray(bonds, dtype=int)
        f.bonds.N, f.bonds.types = int(b.shape[0]), ["A-A"]
        f.bonds.typeid = np.zeros(b.shape[0], dtype=int)
        f.bonds.group = b
    if angles is not None:
        a = np.asarray(angles, dtype=int)
        f.angles.N, f.angles.types = int(a.shape[0]), ["A-A-A"]
        f.angles.typeid = np.zeros(a.shape[0], dtype=int)
        f.angles.group = a
    return f


def unwrap(positions, images, L, dims: int = 3):
    """`r_unwrapped = r_wrapped + image * L`. **Required for any MSD under PBC.**

    ⚠ This is the mirror image of the worst error in this project's history. A
      *missing* minimum image gave **+1856 %** on a pair quantity; applying the
      minimum image (or reading raw wrapped coordinates) to an MSD is the same
      class of mistake in the other direction -- the MSD silently saturates at
      `~L^2` instead of growing, which looks exactly like a plateau, i.e. exactly
      like a real physical result.

      A free probe run to 1500 tau_B wraps the box many times, so this is not an
      edge case for `bead-water-3d`; it is the common path.
    """
    p = np.asarray(positions, dtype=float)
    im = np.asarray(images, dtype=float)
    if p.shape != im.shape:
        raise ValueError(f"positions {p.shape} and images {im.shape} must match")
    Ls = np.asarray([float(L)] * 3 if np.isscalar(L) else [float(x) for x in L],
                    dtype=float)
    return (p + im * Ls)[..., :int(dims)]


def attach_brownian(sim, dt_star: float, forces, kT: float = 1.0, gamma: float = 1.0):
    """Attach the dimensionless BD integrator and return (integrator, method).

    `kT=1, gamma=1` follows from the thermal convention
    (sigma=d, E=kT, tau=tau_B -> gamma=1, D_t=1).
    """
    import hoomd
    import hoomd.md as md

    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=kT, default_gamma=gamma)
    integ = md.Integrator(dt=dt_star, methods=[bd], forces=list(forces))
    integ.integrate_rotational_dof = False          # trap 3
    sim.operations.integrator = integ
    return integ, bd


def add_trajectory_writer(sim, path, period: int):
    """The Tier A trajectory. With `path=None`, nothing is attached."""
    import hoomd

    if path is None:
        return None
    p = max(1, int(period))
    wr = hoomd.write.GSD(filename=str(path), trigger=hoomd.trigger.Periodic(p),
                         mode="xb", dynamic=["property"])
    sim.operations.writers.append(wr)
    return wr


def flush_writers(sim) -> None:
    for w in sim.operations.writers:
        if hasattr(w, "flush"):
            w.flush()


def wca(nlist, epsilon: float = 1.0, sigma: float = 1.0, types=("A", "A")):
    """WCA = LJ with cutoff `2^(1/6)*sigma` and shift (trap 4: there is no
    dedicated WCA class).
    """
    import hoomd.md as md

    lj = md.pair.LJ(nlist=nlist, default_r_cut=2 ** (1 / 6) * sigma, mode="shift")
    lj.params[types] = dict(epsilon=epsilon, sigma=sigma)
    return lj


# ── minimum image: ONE definition ────────────────────────────────────────────
#  ★ Merged 2026-08-29. This existed in **three** places: here, inlined in
#    `bdbot/traps.py::_delta`, and inlined in
#    `bdbot/health.py::measure_step_displacement`. Three copies of the single
#    correction whose absence measured **+1856 %** is not a formatting issue.
#    ⚠ And the copies were not identical: `traps.py` hardcoded `period[2] = 0`, so
#      a **3D** trap would not have wrapped z while these two do. No 3D trap
#      exists yet (`cases/network_3d.py` does not use `traps`), so nothing
#      produced a wrong number -- but the divergence was already there.
#      `period_array(L, dims=3)` now gives one answer for all three callers.
def period_array(L, dims: int = 2):
    """The periodic-length vector `(3,)`. **`0` means "do not wrap this axis".**

    `L` is a scalar or `(Lx, Ly)` (a rectangular box, as required by a
    commensurate hexagonal lattice). In 2D, z gets `0`.
    * Setting a non-periodic axis's period to `inf` gives
      `inf*round(0/inf) = nan` (trap 7). `0` plus the `period > 0` mask is the
      convention that avoids it -- do not "simplify" this to `inf`.
    """
    Lx, Ly = (L, L) if np.isscalar(L) else (L[0], L[1])
    Lz = (float(L[2]) if (not np.isscalar(L) and len(L) > 2)
          else (float(Lx) if dims == 3 else 0.0))
    return np.array([float(Lx), float(Ly), Lz])


def wrap_minimum_image(delta, period):
    """Minimum image against a **precomputed** `period` vector, in place.

    Separate from `minimum_image` because `bdbot.traps` calls this every step and
    should not rebuild the period array each time.
    """
    m = period > 0
    delta[:, m] -= period[m] * np.round(delta[:, m] / period[m])
    return delta


def minimum_image(delta, L, dims: int = 2):
    """Apply minimum image to the periodic axes only (traps 1 and 7).
    `delta` is (N,3); in 2D, z is left alone.
    """
    return wrap_minimum_image(np.asarray(delta, dtype=float).copy(),
                              period_array(L, dims))


def progress(i, total, t_elapsed, extra: str = "") -> str:
    pct = 100 * i / total if total else 0.0
    return f"    {i:>6}/{total}  ({pct:4.0f}%)  {t_elapsed:6.1f}s   {extra}"


__all__ = ["resolve_seed", "frame_2d", "frame_3d", "unwrap", "make_sim", "attach_brownian",
           "add_trajectory_writer",
           "flush_writers", "wca", "period_array", "wrap_minimum_image", "minimum_image",
           "progress", "HOOMD_SEED_MAX"]
