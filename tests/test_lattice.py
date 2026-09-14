"""`bdbot.lattice` — the commensurate hexagonal tiling, and the two refusals.

A lattice that does not close on itself across the periodic boundary carries a
seam, and the seam is a line of defects. For `soft-r3`'s melting bracket that
matters more than usual: the question being asked is *does the crystal melt*,
and a seam melts on its own, so a silent commensurability error would answer it
in the wrong direction.

The module was promoted out of `cases/trap_drag_2d.py` on its third occurrence.
`test_reproduces_trap_drags_own_choice` is the one that ties the promotion to
the case it came from.
"""
from __future__ import annotations

import collections
import math

import numpy as np
import pytest

from bdbot import lattice as LAT
from bdbot.pairpot import HEX_NN, a_mean_star

PHI = 0.35
A_NN = HEX_NN * a_mean_star(PHI)


# ── the tiling is exact ────────────────────────────────────────────────────

@pytest.mark.parametrize("n", [144, 256, 306, 400, 576, 700, 1024])
def test_the_box_is_filled_at_the_declared_packing(n):
    """`L_x*L_y = N a_mean^2` identically, so phi is conserved by construction.
    If this drifts, the run is at a different density than the spec says."""
    n_x, n_y = LAT.commensurate(n)
    LAT.check(n_x, n_y, n)
    lx, ly = LAT.box_for(n_x, n_y, A_NN)
    assert math.isclose(math.pi * n / (4 * lx * ly), PHI, rel_tol=1e-12)


@pytest.mark.parametrize("n", [144, 400, 576])
def test_the_lattice_is_a_defect_free_crystal(n):
    """psi_6 = 1 and every particle six-fold, measured with the same freud calls
    the case uses. This is the initial condition of the crystal arm: if it is not
    perfect, "the crystal melted" is not a measurement of anything."""
    freud = pytest.importorskip("freud")
    n_x, n_y = LAT.commensurate(n)
    lx, ly = LAT.box_for(n_x, n_y, A_NN)
    pos = LAT.hex_lattice(n_x, n_y, A_NN)
    p3 = np.c_[pos, np.zeros(len(pos))]
    box = freud.box.Box(Lx=lx, Ly=ly, is2D=True)

    psi6 = abs(freud.order.Hexatic(k=6, weighted=True)
               .compute((box, p3)).particle_order.mean())
    assert psi6 == pytest.approx(1.0, abs=1e-9), psi6

    nlist = freud.locality.Voronoi().compute((box, p3)).nlist
    counts = collections.Counter(np.asarray(nlist.neighbor_counts).tolist())
    assert set(counts) == {6}, dict(counts)
    assert np.asarray(nlist.distances).mean() == pytest.approx(A_NN, rel=1e-6)


@pytest.mark.parametrize("n", [144, 400, 576])
def test_no_pair_overlaps_across_the_periodic_boundary(n):
    """The nearest-neighbour distance must survive the minimum image. Starting
    overlapped throws a particle out of the box on the first step under
    overdamped dynamics (.claude/rules/overdamped-stability.md), and the symptom
    is a quiet box escape, not NaN."""
    freud = pytest.importorskip("freud")
    n_x, n_y = LAT.commensurate(n)
    lx, ly = LAT.box_for(n_x, n_y, A_NN)
    p3 = np.c_[LAT.hex_lattice(n_x, n_y, A_NN), np.zeros(n)]
    box = freud.box.Box(Lx=lx, Ly=ly, is2D=True)
    nlist = freud.locality.Voronoi().compute((box, p3)).nlist
    dmin = box.compute_distances(p3[nlist.query_point_indices],
                                 p3[nlist.point_indices]).min()
    assert dmin == pytest.approx(A_NN, rel=1e-5)
    assert dmin > 1.0, "closer than one diameter: the WCA core would explode"


def test_positions_lie_inside_the_box():
    n_x, n_y = LAT.commensurate(400)
    lx, ly = LAT.box_for(n_x, n_y, A_NN)
    pos = LAT.hex_lattice(n_x, n_y, A_NN)
    assert pos[:, 0].min() >= -lx / 2 - 1e-12 and pos[:, 0].max() < lx / 2
    assert pos[:, 1].min() >= -ly / 2 - 1e-12 and pos[:, 1].max() < ly / 2


# ── the refusals fire ──────────────────────────────────────────────────────

def test_a_wrong_product_is_refused():
    with pytest.raises(ValueError, match="commensurate hexagon broken"):
        LAT.check(20, 20, 399)


def test_an_odd_row_count_is_refused():
    """It still gives the right N, which is why it needs its own refusal: the
    stagger has a period of two rows, so an odd count does not join."""
    with pytest.raises(ValueError, match="odd"):
        LAT.check(20, 21, 420)


def test_a_valid_pair_is_not_refused():
    """Guard on the guard — a gate that refuses everything is worse than none."""
    LAT.check(20, 20, 400)


def test_an_untileable_n_is_refused_rather_than_approximated():
    """`simbot.build.hex_tiling_for(400)` returns (18, 11) -> 396 particles while
    the caller still believes 400. A 1 % silent drop in N is the failure this
    refuses instead."""
    with pytest.raises(ValueError, match="no hexagonal tiling"):
        LAT.commensurate(7)             # prime: no even divisor


@pytest.mark.parametrize("n", [144, 256, 306, 400, 576, 700, 1024])
def test_commensurate_never_approximates_n(n):
    n_x, n_y = LAT.commensurate(n)
    assert n_x * n_y == n
    assert len(LAT.hex_lattice(n_x, n_y, A_NN)) == n


# ── it is the same lattice the case it came from derived by hand ───────────

def test_reproduces_trap_drags_own_choice():
    """`trap-drag-2d-hex300` pins N = n_x*n_y = 17x18 = 306, derived by hand in
    that case under three constraints. `commensurate` reaches the same pair from
    the aspect ratio alone — an independent check on both."""
    assert LAT.commensurate(306) == (17, 18)


def test_the_square_count_aspect_is_two_over_root_three():
    """n_x = n_y gives L_x/L_y = 2/sqrt3 = 1.1547, which is why a commensurate
    box is never square and why the minimum image must be read off the SHORT
    side. Reading it off the area-equivalent square is 7.5 % optimistic."""
    n_x, n_y = LAT.commensurate(400)
    assert n_x == n_y
    lx, ly = LAT.box_for(n_x, n_y, A_NN)
    assert lx / ly == pytest.approx(2 / math.sqrt(3), rel=1e-12)
    #  ...and that factor is exactly HEX_NN
    l_square = a_mean_star(PHI) * math.sqrt(400)
    assert l_square / ly == pytest.approx(HEX_NN, rel=1e-12)
