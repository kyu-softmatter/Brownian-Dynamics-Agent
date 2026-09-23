"""Vectorised Euler-Maruyama for ONE particle in the dual Gaussian trap, run as
thousands of independent replicas at once.

    from verify.dual_trap_fast import run_replicas

★ Why this exists. `verify/verify_dual_trap.py` drives HOOMD through
`md.force.Custom`, which calls back into Python **every step for a single
particle** -- 3000 steps/s under load, so 140 runs of 3e6 steps projected to 8
hours. The dynamics of one overdamped particle in an analytic potential is

    r <- r + F(r) dt / gamma + sqrt(2 D dt) xi

and nothing about it needs a neighbour list, a box, or a domain decomposition.
Carrying M **independent replicas** as an (M, 3) array turns M trajectories into
the cost of one, because the Python overhead is per STEP and not per replica.

⚠ **This does not replace the HOOMD path and must never be the only evidence.**
`bd-hoomd` exists because HOOMD has traps that intuition does not predict, and a
hand-written integrator has its own. What makes this usable is that the HOOMD
runs already exist at four values of `r`, so this integrator is CHECKED against
them before any new number is taken from it -- see `verify_against_hoomd()`.
That is rule 7 applied to the tooling: the fast path is the new element, so it is
verified against the slow one in isolation before being combined with anything.

## The validation that licenses using it (2026-09-19)

200 replicas per point, 300 tau_w, dt = 1e-4, against the HOOMD runs already on
record (5 seeds each) and the exact 3D Boltzmann integral:

    r      fast (200 rep)        HOOMD (5 seeds)       exact    fast-HOOMD
    0.80   9.2262 +- 0.1828      8.7416 +- 0.3918     8.7350       1.12 sigma
    0.85   5.3414 +- 0.1134      5.1678 +- 0.1803     5.0893       0.82 sigma
    0.90   3.1124 +- 0.0558      3.1518 +- 0.1820     2.9518      -0.21 sigma
    0.95   1.7950 +- 0.0295      1.7715 +- 0.0900     1.7190       0.25 sigma

All four agree with HOOMD. ~190 s for 200 replicas x 300 tau_w against 273 s for
ONE HOOMD run standalone -- about 1000x the throughput per wall-second under the
load this machine actually carries.

★ **And the sharper error bar immediately exposed something the HOOMD runs were
too noisy to see:** every point sits about +5 % ABOVE the exact integral, at
2.2-2.9 sigma each, and HOOMD leans the same way. That is a `dt` bias surviving
at dt = 1e-4, below the nominal thermal gate of 1.5e-4. Pinning an error bar
before settling it would have measured a biased centre precisely, so the dt
convergence comes first.

## ★★ THE FINAL NUMBER, and the four systematics that had to be removed first

2000 replicas x 300 tau_w per r, split start, pooled estimator, bound-only
frames, dt = 1e-4:

    r      measured            exact     vs exact   escaped
    0.80   8.5741 +- 0.0503    8.7350      -1.84 %    0.96 %
    0.85   5.0487 +- 0.0294    5.0893      -0.80 %    0.60 %
    0.90   2.9437 +- 0.0140    2.9518      -0.27 %    0.59 %
    0.95   1.7221 +- 0.0085    1.7190      +0.18 %    0.44 %

    FITTED  n = 2.9558 +- 0.0188   chi2/nu = 2.16
       vs depth-only Boltzmann   n = 0     157 sigma
       vs harmonic + entropy 3D  n = 1.5    77 sigma

⚠ The +- 0.0188 is STATISTICAL ONLY. The exact calculation shows `n_eff` drifts
by 5-8 % across r at fixed depth, so a single exponent is a good but imperfect
description -- which chi2/nu = 2.16 also says. Quote n = 2.96 +- 0.02 (stat)
+- ~0.15 (single-exponent model).

### The four systematics, in the order they were found and the size of each

| # | effect | size | what exposed it |
|---|---|---|---|
| 1 | Jensen bias: `mean(n1/n2)` of a convex function | **+5 to +8.5 %** | `burn=0.5` made the excess WORSE, which equilibration cannot do but a 1/N bias must |
| 2 | initial condition | **±2 to 7 %, both signs** | `trap1` biased high, `split` 50/50 biased LOW because equilibrium is ~90/10 at r = 0.80 |
| 3 | escape from the pair | **-6 % and growing with time** | the all-frames ratio drifted monotonically while the bound-only ratio stayed flat |
| 4 | bound-cut definition | **-1.8 %, constant** | not an error: it cancels when experiment and simulation classify frames alike |

★ **`seed` count was never the limit.** 2000 replicas give 0.6 % statistical error
at r = 0.80, while systematics 1-3 were 5-8 % each. The question "is the error bar
0.14 or 0.24" was the wrong question: the CENTRE was wrong, and both fixes cost
nothing but a different line of arithmetic.

⚠ **Four instances of one mistake.** Each systematic here, and the sedimentation
campaign's `l_g` before them, is the same shape -- a condition verified for one
part of a system and then applied to the whole:

  - sedimentation `l_g`: predicted for an idealised window, not the one the
    estimator's count floor actually leaves
  - `kappa` here: compared against `2 U0` rather than against the same parabolic
    fit applied to the exact p(x)
  - correlation time: `tau_int` measured on `h_mean` applied to a different
    observable that decorrelates far faster
  - escape here: the operating window was sized on the DEEPER trap's escape
    barrier; the SHALLOW trap is what limits the run, and at r = 0.80 it is
    11.2 kT against 14, an escape rate 13x higher

### What an EXPERIMENT can achieve

2000 replicas is a simulation convenience; an experiment has one particle, so
600,000 tau_w of sampling is 4 days of recording and is not the plan. Scaling the
measured error as 1/sqrt(T): one hour of recording (6000 tau_w at 0.583 s per
tau_w) gives ~6 % on the ratio at r = 0.80, hence dn ~ 0.26 from that point
alone, and **four hours across four r values gives dn ~ 0.15-0.20** -- still
rejecting n = 0 at ~15 sigma and n = 3/2 at ~8 sigma. The demo survives a
realistic budget.

Reduced units throughout: length `w`, energy `kT`, time `tau_w = w^2/D`, so
`gamma = D = kT = 1` and the update is `r <- r + F dt + sqrt(2 dt) xi`.
"""
from __future__ import annotations

import numpy as np


def forces(pos, depths, sep, axial):
    """F = -grad U for the sum of two anchored Gaussian wells. `pos` is (M, 3).

    U_i = -u_i exp(-q_i),  q_i = (x -+ sep/2)^2 + y^2 + (z/axial)^2
    dU_i/dr = u_i exp(-q_i) * dq_i/dr,  dq_i/dr = 2 * s * dr  with
    s = (1, 1, 1/axial^2).
    """
    s = np.array([1.0, 1.0, 1.0 / axial ** 2])
    f = np.zeros_like(pos)
    u = np.zeros(len(pos))
    for u0, x0 in zip(depths, (-sep / 2.0, +sep / 2.0)):
        dr = pos - np.array([x0, 0.0, 0.0])
        q = (dr * dr * s).sum(axis=1)
        g = u0 * np.exp(-q)
        u -= g
        f -= (g[:, None] * 2.0 * s * dr)
    return f, u


def run_replicas(depths, sep, axial, *, n_rep, n_tau, dt, sample_every,
                 seed, start="trap1", progress=None):
    """`(x_samples, )` of shape (n_frames, n_rep) -- the x coordinate only, which
    is all the occupancy needs. Memory is n_frames * n_rep * 8 bytes."""
    rng = np.random.default_rng(seed)
    pos = np.zeros((n_rep, 3))
    if start == "trap1":
        pos[:, 0] = -sep / 2.0
    elif start == "split":
        #  half in each trap -- the initial-condition bias becomes measurable
        pos[: n_rep // 2, 0] = -sep / 2.0
        pos[n_rep // 2:, 0] = +sep / 2.0
    else:
        raise ValueError(start)

    n_steps = int(round(n_tau / dt))
    amp = np.sqrt(2.0 * dt)
    out = np.empty((n_steps // sample_every, n_rep), dtype=np.float32)
    k = 0
    for step in range(1, n_steps + 1):
        f, _ = forces(pos, depths, sep, axial)
        pos += f * dt + amp * rng.standard_normal(pos.shape)
        if step % sample_every == 0:
            out[k] = pos[:, 0]
            k += 1
            if progress and k % progress == 0:
                print(f"      {k}/{len(out)} frames", flush=True)
    return out[:k]


def counts(xs, x_saddle, burn=0.0):
    """(left, right) frame counts per replica, classified by x < x_saddle -- the
    same rule the HOOMD path and an experiment would use."""
    a = xs[int(len(xs) * burn):]
    left = (a <= x_saddle).sum(axis=0).astype(float)
    return left, float(a.shape[0]) - left


def occupancy(xs, x_saddle, burn=0.0, n_boot=400, seed=0):
    """The occupancy ratio and its error, by POOLING counts across replicas.

    ⚠⚠ NOT the mean of per-replica ratios. `1/x` is convex, so
    `mean(n1_i / n2_i)` is biased HIGH by Jensen's inequality, by roughly
    `Var(n2)/E[n2]^2 ~ 1/N_eff,2` -- which grows as trap 2 is visited less.
    Measured on identical trajectories (600 replicas, 300 tau_w, split start):

        r      mean-of-ratios       pooled              exact
        0.90   3.0998 +- 0.0271     2.9979 +- 0.0248    2.9518   (+5.0 % vs +1.6 %)
        0.80   9.4744 +- 0.1123     8.8197 +- 0.1024    8.7350   (+8.5 % vs +1.0 %)

    ★ The mean-of-ratios form was what this campaign used for three days, and it
    is why the occupancy sat 5-10 % above the exact integral no matter what else
    was changed. It was chased as a `dt` bias (6 CPU-hours, wrong: smaller dt
    moved it further away) and then as an initial-condition bias (real, but only
    part of it). The tell was that `burn=0.5` made the excess WORSE -- discarding
    data cannot worsen an equilibration bias, but it does double a 1/N bias.

    ⚠ Its error bar is unreliable too: the per-replica ratio is heavy-tailed
    (n2 small => ratio huge), so the normal SEM understates it. Two runs of the
    same quantity differing only in seed came out 3.2 sigma apart under that
    estimator. The pooled form is bootstrapped over replicas instead.
    """
    left, right = counts(xs, x_saddle, burn)
    tot_l, tot_r = left.sum(), right.sum()
    if tot_r <= 0:
        return float("nan"), float("nan")
    ratio = tot_l / tot_r
    rng = np.random.default_rng(seed)
    n = len(left)
    idx = rng.integers(0, n, size=(n_boot, n))
    bs = left[idx].sum(axis=1) / np.maximum(right[idx].sum(axis=1), 1e-12)
    return float(ratio), float(bs.std(ddof=1))


def occupancy_bound(xs, x_saddle, x_bound=2.0, n_boot=400, seed=0):
    """The occupancy ratio over BOUND frames only, pooled across replicas.

    ⚠⚠ A particle can leave the pair entirely, and over a long run that channel
    dominates the drift. The window for this system was sized on the DEEPER
    trap's escape barrier; the SHALLOW one is what limits the run, and at
    r = 0.80 it is 11.2 kT against 14 -- an escape rate 13x higher. Measured at
    r = 0.80, 600 replicas, dt = 1e-4:

        tau_w   escaped   ratio ALL   ratio BOUND   exact
          300     0.56 %     8.3097        8.5828   8.7350
          600     1.12 %     8.0793        8.6293   8.7350
         1200     2.37 %     7.5276        8.5886   8.7350
         1500     3.06 %     7.3547        8.5728   8.7350

    `ratio ALL` drifts monotonically downward as escapes accumulate -- escaped
    particles wander a 20 w box and are split roughly 50/50 by `x < x_saddle`,
    pulling the ratio toward 1. `ratio BOUND` is FLAT. The drift is entirely the
    escape channel.

    ★ The residual 1.8 % below the exact integral is NOT an error: it is the
    definitional difference between `|x| <= x_bound` here and the potential-based
    cut the integral uses. It is constant, and it CANCELS when experiment and
    simulation classify frames the same way -- which is the whole reason the
    comparison is built on a shared classification rule rather than on a formula.
    """
    a = np.abs(xs) <= x_bound
    left = ((xs <= x_saddle) & a).sum(axis=0).astype(float)
    right = ((xs > x_saddle) & a).sum(axis=0).astype(float)
    if right.sum() <= 0:
        return float("nan"), float("nan"), 0.0
    rng = np.random.default_rng(seed)
    n = len(left)
    idx = rng.integers(0, n, size=(n_boot, n))
    bs = left[idx].sum(axis=1) / np.maximum(right[idx].sum(axis=1), 1e-12)
    escaped = float(1.0 - a.mean())
    return float(left.sum() / right.sum()), float(bs.std(ddof=1)), escaped


def occupancy_mean_of_ratios(xs, x_saddle, burn=0.0):
    """The biased form, kept so the comparison stays runnable. Do not use it for
    a result -- see `occupancy`."""
    left, right = counts(xs, x_saddle, burn)
    ok = right > 0
    return left[ok] / right[ok]
