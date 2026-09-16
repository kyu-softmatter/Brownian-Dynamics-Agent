"""Initial placement with no overlap — for a slab, i.e. periodic in x and y and
walled in z.

**Promoted on the third occurrence.** `cases/soft_r3_2d.rsa_positions` (2D,
fully periodic) and `cases/network_3d.scatter_no_overlap` (3D, fully periodic)
are the first two, and `passive-sphere--sedimentation` is the third — but it
cannot use either, because a sedimentation cell is **periodic in x and y and
bounded by walls in z**, so the minimum image must be applied per axis and the
height is not uniform.

⚠ **The two predecessors are deliberately NOT switched over.** They consume the
RNG in their own order, and `run_id` is the content hash of the spec while the
trajectory is reproduced from `spec + seed` — so changing the stream would make
every archived run of those two cases unreproducible under its own name. Same
reasoning as `bdbot/runid.py`'s note on `spec_hash`: two implementations kept
apart on purpose, with the reason written down.

## Why this module exists at all rather than uniform random placement

`.claude/rules/overdamped-stability.md`: *"When adding strong repulsion, check the
minimum separation of the initial placement. Start overlapped and the first step
throws the particle out."* Measured here on 2026-09-16, placing 584 WCA particles
at uniform random `x, y` with exponential `z`:

    RuntimeError: Particle with unique tag 11 is no longer in the simulation box.
    Cartesian coordinates: x: 39.3173  y: 7.96552  z: -44.5194
    Local box lo: (-6, -6, -62.5)  hi: (6, 6, 62.5)

A particle 39 diameters outside a box 12 wide, on the **first few thousand
steps**, because two particles were placed on top of each other and WCA's `r⁻¹³`
core did the rest. No NaN; the run simply stopped.
"""
from __future__ import annotations

import numpy as np

#: Per-particle attempt budget before giving up. RSA in 3D jams near a volume
#: fraction of 0.38, so a failure here means the requested density is not
#: reachable by sequential addition, not that the budget was mean.
TRIES_PER_PARTICLE = 400


def exponential_height(l_g: float, lz: float, *, margin: float):
    """A height sampler for sedimentation equilibrium: `p(z) ∝ exp(−z/l_g)`.

    Starting from the analytic profile is not a shortcut around equilibration, it
    is what makes the run affordable: from a uniform start the system has to
    *fall*, which takes `L_z · l_g / D_0` — 1672 `τ_d` for this cell — while the
    slowest mode of the equilibrium state is `l_g²/D_0 = 179 τ_d`, nine times
    shorter. `verify/verify_sedimentation_wall.py` measures both, and measures
    that the uniform start does converge to the same profile.

    Args:
        l_g: gravitational length, in the same units as `lz`.
        lz: cell height. Samples are redrawn until they land inside the margins.
        margin: keep this far from both walls (one particle radius, usually).
    """
    if not (l_g > 0 and lz > 2 * margin):
        raise ValueError(f"l_g={l_g} and lz={lz} with margin={margin} leave no room")

    def sample(rng, n):
        out = np.empty(n)
        k = 0
        while k < n:
            z = margin + rng.exponential(l_g, size=max(n - k, 16))
            z = z[z < lz - margin]
            take = min(len(z), n - k)
            out[k:k + take] = z[:take]
            k += take
        return out

    return sample


def rsa_slab(n: int, *, lxy: float, lz: float, min_sep: float, rng,
             height_sampler=None, margin: float = 0.0,
             tries_per_particle: int = TRIES_PER_PARTICLE) -> np.ndarray:
    """Random sequential addition in a slab. Returns `(n, 3)` positions.

    The box is assumed **centred on the origin**, as HOOMD's is: `x, y ∈
    [−lxy/2, lxy/2]` and `z ∈ [−lz/2 + margin, lz/2 − margin]`.

    ★ The minimum image is applied **to x and y only**. Applying it to z as well
    is the error this module exists to make impossible: it would let a particle
    at the bottom of the cell be judged "too close" to one at the top, so RSA
    would refuse placements that are physically fine, and — worse — two particles
    genuinely adjacent across the (non-existent) z boundary would be accepted.

    Args:
        min_sep: reject a candidate closer than this to any placed particle.
            For a WCA core of diameter `d`, `min_sep = d` puts the closest pair at
            the cutoff where `F = 24 kT/d`; below about `0.9 d` the force passes
            `139 kT/d` and starts to set `dt` (`bdbot.dt.dt_max_force`).
        height_sampler: `f(rng, k) -> k heights`, in box coordinates measured
            from the **bottom** of the cell. `None` means uniform.
    """
    if n < 0:
        raise ValueError(f"n={n}")
    if min_sep < 0:
        raise ValueError(f"min_sep={min_sep}")
    pos = np.empty((n, 3))
    lv = np.array([float(lxy), float(lxy)])
    lo, hi = -lz / 2 + margin, lz / 2 - margin
    if hi <= lo:
        raise ValueError(f"margin={margin} leaves no room in lz={lz}")

    #  ★★ The heights are drawn ONCE, up front, and a rejection retries only
    #     `x, y` at that FIXED height.
    #
    #     ⚠ The first version redrew `z` on every rejection, and that biases the
    #     profile: rejection is most likely exactly where the density is highest,
    #     so the low-`z` end gets thinned out. Measured on 2000 particles in a
    #     24 x 24 x 125 slab with `min_sep = 1`, against an imposed
    #     `l_g = 13.375 d`:
    #
    #         fit window h in [2.0, 80.2]   l_g = 16.43 d   +22.9 %
    #         fit window h in [10.0, 80.2]  l_g = 16.65 d   +24.5 %
    #         fit window h in [20.0, 80.2]  l_g = 15.79 d   +18.0 %
    #
    #     Not a wall-layer artefact -- the dilute window is off by 18 % too. A
    #     sedimentation run started from that profile has to drift back, and
    #     drift is the slowest thing in the problem (`L_z l_g/D_0`). Fixing the
    #     height removes the bias at its source instead of paying for it in run
    #     length.
    heights = None
    if height_sampler is not None:
        heights = np.asarray(height_sampler(rng, n), dtype=float) - lz / 2
        if np.any((heights < lo - 1e-12) | (heights > hi + 1e-12)):
            raise ValueError("height_sampler returned a height outside the margins")

    k, attempts, redraws = 0, 0, 0
    budget = tries_per_particle * max(n, 1)
    while k < n:
        attempts += 1
        if attempts > budget:
            raise RuntimeError(
                f"RSA failed after {attempts} attempts: placed {k} of {n} in a "
                f"{lxy}x{lxy}x{lz} slab at min_sep={min_sep} ({redraws} height "
                f"redraws). RSA jams near a volume fraction of 0.38 in 3D, so "
                f"either the density is not reachable sequentially or min_sep is "
                f"too large.")
        xy = rng.uniform(-lv / 2, lv / 2, 2)
        z = rng.uniform(lo, hi) if heights is None else float(heights[k])
        cand = np.array([xy[0], xy[1], z])
        if k:
            dr = pos[:k] - cand
            dr[:, :2] -= lv * np.round(dr[:, :2] / lv)      # x, y only  ★
            if (dr ** 2).sum(axis=1).min() < min_sep ** 2:
                #  a height that cannot be placed laterally after many tries is
                #  redrawn rather than looped on forever -- but it is COUNTED, so
                #  the bias that the first version hid is visible if it returns.
                if heights is not None and attempts % (8 * tries_per_particle) == 0:
                    heights[k] = float(height_sampler(rng, 1)[0]) - lz / 2
                    redraws += 1
                continue
        pos[k] = cand
        k += 1
    return pos


def min_separation_slab(pos: np.ndarray, *, lxy: float) -> float:
    """The closest pair distance under the slab's boundaries. Use it to assert on
    a placement rather than trusting that RSA did what it was asked."""
    if len(pos) < 2:
        return float("inf")
    lv = np.array([float(lxy), float(lxy)])
    best = np.inf
    for i in range(len(pos) - 1):
        dr = pos[i + 1:] - pos[i]
        dr[:, :2] -= lv * np.round(dr[:, :2] / lv)
        best = min(best, float(np.sqrt((dr ** 2).sum(axis=1)).min()))
    return best


__all__ = ["exponential_height", "rsa_slab", "min_separation_slab",
           "TRIES_PER_PARTICLE"]
