"""Commensurate lattices for a periodic box.

A lattice that does not close on itself across the periodic boundary carries a
seam, and the seam is a line of defects that melts on its own. So the box aspect
ratio is **not a free choice** -- commensurability sets it, and `N` follows from
the tiling rather than being picked.

Third occurrence, hence the promotion (CLAUDE.md: abstract what has appeared
twice). `cases/trap_drag_2d.py` needed it because its observables are the
lattice deformation field and `psi_6`, which are sensitive to exactly that
defect; `cases/soft_r3_2d.py` needs it for the crystal-start arm of the melting
bracket; `simbot/build.py` has its own copy under a different parametrisation
(`(n_x, n_y)` there == `(n_x, 2*n_y)` here -- the same lattice, which is why
`n_y` must be even here).
"""
from __future__ import annotations

import math

import numpy as np

ROW = math.sqrt(3) / 2          #: row spacing / a_NN for a hexagonal lattice


def hex_lattice(n_x: int, n_y: int, a_nn: float) -> np.ndarray:
    """Commensurate hexagonal lattice (`n_x*n_y` particles), centred on the origin.

    Rows run along x; odd rows are staggered by `a_NN/2`. The box is
    `[-L_x/2, L_x/2) x [-L_y/2, L_y/2)` with `L_x = n_x a_NN` and
    `L_y = n_y (sqrt3/2) a_NN`. `n_y` must be even for the stagger to join across
    the periodic boundary.
    """
    Lx, Ly = box_for(n_x, n_y, a_nn)
    j, i = np.divmod(np.arange(n_x * n_y), n_x)
    x = (i + 0.5 * (j % 2)) * a_nn - Lx / 2
    y = j * (ROW * a_nn) - Ly / 2
    return np.c_[x, y]


def box_for(n_x: int, n_y: int, a_nn: float) -> tuple[float, float]:
    """`(L_x, L_y)` of the box that `hex_lattice(n_x, n_y, a_nn)` tiles exactly."""
    return n_x * a_nn, n_y * ROW * a_nn


def check(n_x: int, n_y: int, n: int) -> None:
    """Raise unless `(n_x, n_y)` tiles a periodic box with exactly `n` particles.

    Two separate refusals, because they fail differently: a wrong product means
    the run integrates a different number of particles than the spec declares,
    while an odd `n_y` still gives the right `N` and quietly leaves a seam.
    """
    if n_x * n_y != n:
        raise ValueError(f"commensurate hexagon broken: n_x*n_y = {n_x}*{n_y} = "
                         f"{n_x * n_y} != N = {n}")
    if n_y % 2:
        raise ValueError(f"n_y = {n_y} is odd -- the staggered rows have a period of "
                         f"2 rows, so they do not join across the periodic boundary")


def commensurate(n: int) -> tuple[int, int]:
    """The `(n_x, n_y)` tiling of exactly `n` particles whose box is closest to
    square. Raises if `n` admits none.

    `N` is **not** approximated. `simbot.build.hex_tiling_for(400)` returns
    `(18, 11)` -> 396 particles while the caller still believes 400; the 1 %
    silent drop is the failure this refuses instead.

    The box aspect is `L_x/L_y = 2 n_x / (sqrt3 n_y)`, so `n_x = n_y` gives
    `2/sqrt3 = 1.1547` -- the closest to square a square particle count allows.
    """
    best = None
    for n_y in range(2, n + 1, 2):
        if n % n_y:
            continue
        n_x = n // n_y
        aspect = 2 * n_x / (math.sqrt(3) * n_y)
        key = abs(math.log(aspect))                 # symmetric in aspect vs 1/aspect
        if best is None or key < best[0]:
            best = (key, n_x, n_y)
    if best is None:
        raise ValueError(
            f"N = {n} admits no hexagonal tiling with an even row count. "
            f"N must have an even divisor; pick a nearby N that does.")
    return best[1], best[2]


__all__ = ["ROW", "hex_lattice", "box_for", "check", "commensurate"]
