"""The numerics of a soft repulsive pair potential -- `U`, `U''`, and the closest
approach distance (skill `bd-physics` section 6.2).

**Promoted after appearing twice** (the abstraction rule in CLAUDE.md).
`soft-r3-2d-A-sweep` was the only user until `trap-drag-2d-hex300` built the same
hexagonal lattice at the same density with the same pair potential
(`A(d/r)^3 + WCA`) and reused the `r_min -> tau_int -> dt` path verbatim.

The three functions here are **the physics that sets `dt`.** If a value is
silently wrong, the integration is silently wrong. It was fixed three times (see
the `approach_distance` docstring), and that correction has to live in exactly one
place or the two cases diverge.

WARNING: **this is a different layer from `bdbot.interactions`.** That one is the
   catalogue of "what should I choose" (text metadata); this one is the numerics
   of what was chosen.
"""
from __future__ import annotations

import math

import numpy as np

R_WCA = 2 ** (1 / 6)                 # WCA cutoff (in sigma) -- the LJ minimum
# 2D hexagonal lattice: area per particle = (sqrt(3)/2) a_NN^2
#   =>  a_mean == rho^(-1/2) = sqrt(sqrt(3)/2) * a_NN
# * The nearest-neighbour distance is NOT a_mean. It was first taken as a_mean and
#   corrected by a smoke measurement -- and since U'' goes as r^-5, that is a 41%
#   difference.
HEX_NN = math.sqrt(2 / math.sqrt(3))          # a_NN / a_mean = 1.07457


def U_star(rs, A, eps=1.0):
    """Dimensionless potential U/kT, lengths in units of d. A WCA core plus A/r^3
    (the cutoff shift is applied separately).

    WARNING: `U_star(1.0, A)` is **not** `A` but `A + eps` -- at r=d the WCA core
       contributes exactly epsilon (4*eps*(1-1)+eps). With A=100 that is 101, off
       by 1%. The ledger had this written as "= A kT" and the L3 integrity check
       caught it.
    """
    rs = np.asarray(rs, dtype=float)
    w = np.where(rs < R_WCA, 4 * eps * (rs**-12.0 - rs**-6.0) + eps, 0.0)
    return w + A / rs**3


def U2_star(rs, A, eps=1.0):
    """U''(r) [kT/d^2] -- the local stiffness. The denominator of tau_int = gamma/U''."""
    rs = np.asarray(rs, dtype=float)
    w = np.where(rs < R_WCA, 4 * eps * (156 * rs**-14.0 - 42 * rs**-8.0), 0.0)
    return w + 12 * A / rs**5


def approach_distance(A, a_star, eps=1.0, u_max=12.0, lindemann=0.15):
    """Closest approach distance r_min* [d]. The smaller of two criteria. dt comes
    from here.

    (a) pair criterion       the r where U(r) = u_max*kT.
        Why u_max=12: there are ~10^6 pair samples (400 samples x 400 particles x
        6 neighbours), so the extreme of the Boltzmann tail reaches
        beta*U ~ ln(10^6) ~ 14. Setting u_max=5 put the smoke-measured minimum
        distance far inside the prediction.
    (b) vibration criterion  a_NN - 3*sigma_bond   (valid only while the cage lives)
        * the neighbour distance is a_NN = 1.07457*a_mean (hexagonal), not a_mean.
        * the bond-length fluctuation is sigma_bond = sqrt(2)*u1 (u1 = per-component
          rms). The sqrt(2) is easy to drop.
        The Lindemann indicator used for the verdict is sigma_bond/a_NN (in 2D it is
        the *relative* fluctuation that is finite).

    Errors (a) and (b) partially cancelled, so r_min happened to look close -- do
    not rely on that coincidence.

    Returns: (r_min*, which criterion, Lindemann sigma_bond/a_NN, state estimate)
    """
    lo, hi = 0.4, 60.0
    for _ in range(200):                       # bisection
        mid = 0.5 * (lo + hi)
        if float(U_star(mid, A, eps)) > u_max:
            lo = mid
        else:
            hi = mid
    r_pair = 0.5 * (lo + hi)

    a_nn = HEX_NN * a_star
    k1 = 3 * (float(U2_star(a_nn, A, eps)) + (-3 * A / a_nn**4) / a_nn)   # per-component cage stiffness
    if k1 <= 0:
        return r_pair, "pair", float("nan"), "fluid"
    sigma_bond = math.sqrt(2.0 / k1)
    lind = sigma_bond / a_nn
    r_cage = a_nn - 3 * sigma_bond
    if lind < lindemann and r_cage < r_pair:
        return r_cage, "vibration", lind, "crystal"
    return r_pair, "pair", lind, ("crystal" if lind < lindemann else "fluid")


def a_mean_star(phi: float) -> float:
    """Mean spacing a_mean/d = sqrt(pi/(4*phi))  (2D, phi = N*pi*d^2/(4*L^2))."""
    return math.sqrt(math.pi / (4 * phi))


__all__ = ["R_WCA", "HEX_NN", "U_star", "U2_star", "approach_distance", "a_mean_star"]
