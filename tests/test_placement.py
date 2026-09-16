"""`bdbot.placement` — overlap-free placement in a slab.

★ Why this module exists at all, and therefore what these tests have to pin:
`.claude/rules/overdamped-stability.md` says *"Start overlapped and the first step
throws the particle out."* Measured 2026-09-16 with 584 WCA particles placed at
uniform random x, y:

    RuntimeError: Particle with unique tag 11 is no longer in the simulation box.
    Cartesian coordinates: x: 39.3173  y: 7.96552  z: -44.5194
    Local box lo: (-6, -6, -62.5)  hi: (6, 6, 62.5)

39 diameters outside a box 12 wide, in the first few thousand steps, no NaN.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from bdbot import placement as PL

LXY, LZ, LG = 12.0, 125.0, 13.375


def test_the_requested_minimum_separation_is_achieved():
    rng = np.random.default_rng(0)
    pos = PL.rsa_slab(200, lxy=LXY, lz=LZ, min_sep=1.0, rng=rng, margin=0.5)
    assert len(pos) == 200
    assert PL.min_separation_slab(pos, lxy=LXY) >= 1.0


@pytest.mark.parametrize("min_sep", [0.5, 1.0, 1.5])
def test_the_separation_holds_at_several_radii(min_sep):
    rng = np.random.default_rng(1)
    pos = PL.rsa_slab(150, lxy=LXY, lz=LZ, min_sep=min_sep, rng=rng, margin=0.5)
    assert PL.min_separation_slab(pos, lxy=LXY) >= min_sep


def test_the_minimum_image_is_applied_to_x_and_y():
    """★ The load-bearing asymmetry, shown on two particles rather than by
    statistics. Two particles on opposite LATERAL faces are neighbours; two at
    opposite ENDS of the cell are not, because there is a wall between them.

    ⚠ The first version of this test looked for a spurious close pair among 400
    random placements and failed: at that N the tail of the cell is too sparse
    for a lateral coincidence to show up. A deterministic pair says it exactly —
    124.2 d apart in the slab, 0.8 d apart if z is wrapped.
    """
    #  bottom and top of the cell, laterally aligned
    pair = np.array([[0.0, 0.0, -LZ / 2 + 0.4],
                     [0.0, 0.0, LZ / 2 - 0.4]])
    assert PL.min_separation_slab(pair, lxy=LXY) == pytest.approx(LZ - 0.8)

    lv = np.array([LXY, LXY, LZ])
    dr = pair[1] - pair[0]
    dr -= lv * np.round(dr / lv)
    wrapped = float(np.sqrt((dr ** 2).sum()))
    assert wrapped == pytest.approx(0.8), wrapped
    assert wrapped < 1.0 < PL.min_separation_slab(pair, lxy=LXY)

    #  and laterally: opposite faces ARE within the minimum image
    lateral = np.array([[-LXY / 2 + 0.1, 0.0, 0.0],
                        [LXY / 2 - 0.1, 0.0, 0.0]])
    assert PL.min_separation_slab(lateral, lxy=LXY) == pytest.approx(0.2)

    #  RSA therefore refuses that lateral pair and allows the vertical one
    rng = np.random.default_rng(2)
    pos = PL.rsa_slab(400, lxy=LXY, lz=LZ, min_sep=1.0, rng=rng, margin=0.5)
    assert PL.min_separation_slab(pos, lxy=LXY) >= 1.0


def test_the_margin_keeps_particles_off_the_walls():
    rng = np.random.default_rng(3)
    margin = 0.5
    pos = PL.rsa_slab(300, lxy=LXY, lz=LZ, min_sep=1.0, rng=rng, margin=margin)
    h = pos[:, 2] + LZ / 2
    assert h.min() >= margin - 1e-12, h.min()
    assert h.max() <= LZ - margin + 1e-12, h.max()


def test_the_exponential_sampler_reproduces_its_decay_length():
    """The sampler is the initial condition of a sedimentation run, so a wrong
    decay length here starts the run in the wrong state and the drift back is
    the slowest thing in the problem (L_z l_g/D_0)."""
    rng = np.random.default_rng(4)
    s = PL.exponential_height(LG, LZ, margin=0.5)
    z = s(rng, 200_000)
    assert z.min() >= 0.5 and z.max() <= LZ - 0.5
    # mean of a truncated exponential, solved rather than approximated
    lo, hi = 0.5, LZ - 0.5
    a, b = (lo - 0.5) / LG, (hi - 0.5) / LG
    want = 0.5 + LG * (1 - (1 + b) * math.exp(-b)) / (1 - math.exp(-b))
    assert z.mean() == pytest.approx(want, rel=0.02), (z.mean(), want)


def test_rsa_does_not_bias_the_height_distribution():
    """★★ The strongest statement available, and it was earned the hard way.

    The first version redrew `z` on every rejection. Rejection is most likely
    exactly where the density is highest, so the low-`z` end was thinned out and
    the fitted decay length came out **+18 % to +27 %** too long — at every fit
    window, including a dilute one, so it was not a wall artefact. A run started
    from that profile has to drift back, and drift is the slowest thing in the
    problem.

    With `z` drawn once and only `x, y` retried, the placed set is the drawn set
    **exactly** — which is a stronger claim than any distributional test, so it is
    the one asserted.
    """
    s = PL.exponential_height(LG, LZ, margin=0.5)
    wanted = s(np.random.default_rng(11), 2000)
    pos = PL.rsa_slab(2000, lxy=24.0, lz=LZ, min_sep=1.0,
                      rng=np.random.default_rng(11), margin=0.5, height_sampler=s)
    got = pos[:, 2] + LZ / 2
    assert np.allclose(np.sort(wanted), np.sort(got)), (
        "RSA changed the set of heights, so it is biasing the profile again")
    assert PL.min_separation_slab(pos, lxy=24.0) >= 1.0


def test_the_log_histogram_fit_is_the_noisy_part_not_the_placement():
    """⚠ Recorded because it nearly caused a wrong conclusion. Fitting `l_g` from
    `ln n(h)` on a SINGLE placement is biased by the bin count, not by the
    placement: the same fit on the sampler's own output gives the same answer.

    Measured: +1.93 % at 2000 particles (11 usable bins) and +9.17 % at 584
    (5 bins). The production run averages ~160,000 particle-samples and the
    capability check reached 0.02 sigma, so this is a small-sample property of
    the estimator and belongs in the analysis plan's choice of window.
    """
    s = PL.exponential_height(LG, LZ, margin=0.5)
    wanted = s(np.random.default_rng(11), 2000)
    pos = PL.rsa_slab(2000, lxy=24.0, lz=LZ, min_sep=1.0,
                      rng=np.random.default_rng(11), margin=0.5, height_sampler=s)

    def fit(h):
        cnt, e = np.histogram(h, bins=np.linspace(2.0, 5 * LG, 20))
        mid = 0.5 * (e[1:] + e[:-1])
        ok = cnt >= 30
        return -1.0 / np.polyfit(mid[ok], np.log(cnt[ok]), 1)[0], int(ok.sum())

    a, na = fit(wanted)
    b, nb = fit(pos[:, 2] + LZ / 2)
    assert (a, na) == pytest.approx((b, nb)), (a, b)
    assert a == pytest.approx(LG, rel=0.10), (a, na)


def test_it_refuses_rather_than_returning_fewer_than_asked():
    """RSA jams near a volume fraction of 0.38 in 3D. A silent short return would
    change the density the run thinks it has."""
    rng = np.random.default_rng(6)
    with pytest.raises(RuntimeError, match="RSA failed"):
        PL.rsa_slab(500, lxy=4.0, lz=4.0, min_sep=1.5, rng=rng, margin=0.5)


def test_a_margin_larger_than_the_cell_is_refused():
    rng = np.random.default_rng(7)
    with pytest.raises(ValueError):
        PL.rsa_slab(10, lxy=LXY, lz=2.0, min_sep=1.0, rng=rng, margin=1.5)
    with pytest.raises(ValueError):
        PL.exponential_height(LG, 1.0, margin=1.0)


def test_placement_is_reproducible_from_the_seed():
    """`run_id` names a spec and the trajectory is reproduced from `spec + seed`,
    so the placement has to be a function of the seed alone."""
    a = PL.rsa_slab(100, lxy=LXY, lz=LZ, min_sep=1.0,
                    rng=np.random.default_rng(42), margin=0.5)
    b = PL.rsa_slab(100, lxy=LXY, lz=LZ, min_sep=1.0,
                    rng=np.random.default_rng(42), margin=0.5)
    assert np.array_equal(a, b)
    c = PL.rsa_slab(100, lxy=LXY, lz=LZ, min_sep=1.0,
                    rng=np.random.default_rng(43), margin=0.5)
    assert not np.array_equal(a, c)


def test_zero_particles_is_not_an_error():
    pos = PL.rsa_slab(0, lxy=LXY, lz=LZ, min_sep=1.0,
                      rng=np.random.default_rng(8), margin=0.5)
    assert pos.shape == (0, 3)
    assert PL.min_separation_slab(pos, lxy=LXY) == float("inf")
