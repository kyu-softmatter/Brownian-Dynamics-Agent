"""S5 — HOOMD integrator and RNG benchmark regressions. [slow]

Pins the measurements in
`knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md`. Catches a silent break
when the HOOMD version is bumped or forces.py is edited.

The tolerances are **multiples of the theoretical statistical error** (never fitted
to the observation). The seed is fixed -- HOOMD `Brownian` reproduces bit-for-bit on
the same seed.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from simbot.estimators import euler_maruyama_trap_variance_bias
from simbot.guards import assert_statistic_fluctuates

from .conftest import TRAP_DIM, TRAP_N, se_of_mean, sigma_away

pytestmark = pytest.mark.slow

N_SIGMA = 4.0        # 4σ — not flaky because the seed is fixed, and the value comes
                     # from the theoretical error


def _equilibrium_variance(sim, positions_of, dt_star, n_blocks=25):
    """Collects the per-component <x^2> in blocks after equilibration."""
    sim.run(int(20 / dt_star))                       # equilibrate 20 tau_trap
    stride = max(1, int(2.0 / dt_star))              # every 2 tau_trap → independent
    blocks = []
    for _ in range(n_blocks):
        sim.run(stride)
        blocks.append(float(np.mean(positions_of(sim) ** 2)))
    return blocks


# =============================================================================
# B1 — the integration scheme is Euler-Maruyama
# =============================================================================
@pytest.mark.benchmark
@pytest.mark.parametrize("dt_star,n_blocks,require_rejection", [
    # A larger dt* grows the bias (dt*/2), which is what gives discriminating power.
    # require_rejection is decided from the **design power**, not from observation:
    #   expected power = bias / SE.  Demanding 3σ from a design that cannot produce
    #   3σ is an unachievable assert (it was actually written that way 2026-07-28).
    (2.0e-2, 25, True),    # bias 1.01 % → expect ~3.8σ.  rejection required
    (1.0e-2, 25, False),   # bias 0.50 % → expect ~1.8σ.  pins INCONCLUSIVE as a fact
])
def test_B1_brownian_is_euler_maruyama(trap_sim_factory, positions_of, dt_star,
                                      n_blocks, require_rejection):
    """★ <x*^2> = 1/(1-dt*/2). The competing hypothesis is rejected only where there
    is power to do so.

    ★ **Verifying the integrator has to be done at a deliberately large dt*.**
      Reducing dt* reduces the bias and the verifiability along with it -- at the
      production dt* (5e-3) the integrator cannot be verified at all.
      Basis: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md §1
    """
    blocks = []
    for seed in (1, 2, 3, 4):
        sim, _ = trap_sim_factory(dt_star, seed)
        blocks += _equilibrium_variance(sim, positions_of, dt_star, n_blocks=n_blocks)

    assert_statistic_fluctuates(blocks, "block <x^2>")

    measured = float(np.mean(blocks))
    se = se_of_mean(blocks)
    bias = euler_maruyama_trap_variance_bias(dt_star)
    em_pred = 1.0 + bias

    dev_em = sigma_away(measured, em_pred, se)
    dev_exact = sigma_away(measured, 1.0, se)
    power = bias / se                      # design power: the σ separating EM from
                                           # exact

    assert dev_em < N_SIGMA, (
        f"dt*={dt_star:g}: measured {measured:.5f}±{se:.5f} departs from the EM "
        f"prediction {em_pred:.5f} by {dev_em:.1f}σ")

    if require_rejection:
        assert power > 3.0, (
            f"dt*={dt_star:g}: the design power is only {power:.1f}σ — n_blocks has "
            f"to grow for this case to mean anything "
            f"(needed multiple {(3.0/power)**2:.1f}x)")
        assert dev_exact > 3.0, (
            f"dt*={dt_star:g}: the exact-scheme hypothesis (1.0) is rejected at only "
            f"{dev_exact:.1f}σ")
    else:
        # **Pins as a fact** that this is the undecidable regime. If a 3σ rejection
        # has become possible here, the statistics improved and this case should be
        # promoted to require_rejection=True.
        assert power < 3.0, (
            f"dt*={dt_star:g}: the power rose to {power:.1f}σ — "
            f"promote it to require_rejection=True")


@pytest.mark.benchmark
def test_B5_em_bias_halves_when_dt_halves(trap_sim_factory, positions_of):
    """★ Confirms first order: halving dt* halves the bias (second order would give
    1/4)."""
    bias = {}
    for dt_star in (2.0e-2, 1.0e-2):
        blocks = []
        for seed in (11, 12, 13, 14):
            sim, _ = trap_sim_factory(dt_star, seed)
            blocks += _equilibrium_variance(sim, positions_of, dt_star)
        bias[dt_star] = (float(np.mean(blocks)) - 1.0, se_of_mean(blocks))

    b_hi, se_hi = bias[2.0e-2]
    b_lo, se_lo = bias[1.0e-2]
    ratio = b_hi / b_lo
    # error propagation on the ratio
    se_ratio = abs(ratio) * math.hypot(se_hi / b_hi, se_lo / b_lo)

    assert sigma_away(ratio, 2.0, se_ratio) < N_SIGMA, (
        f"the bias ratio {ratio:.2f}±{se_ratio:.2f} departs from the first-order "
        f"prediction 2.0")
    assert abs(ratio - 4.0) > 2 * se_ratio, (
        f"the second-order hypothesis (ratio=4) is not rejected "
        f"(ratio {ratio:.2f}±{se_ratio:.2f})")


# =============================================================================
# B7 — N non-interacting particles are independent samples
# =============================================================================
@pytest.mark.benchmark
def test_B7_particles_are_independent_samples(trap_sim_factory, positions_of):
    """★ Var_t(mean_i d_i) * N / <Var_i(d_i)> has to be 1.

    The entire cost claim of "1 % precision in 1.3 seconds" rests on this identity.
    Break it and the error bars come out smaller than they are, so **a false
    precision gets reported.**

    ⚠ The mean must not be subtracted from the displacements in this check --
      sum(d)=0 then holds identically, the cross correlation is fixed at -1/(n-1),
      and the measurement is meaningless (findings §3).
    """
    dt_star, M = 5e-3, 2000
    sim, _ = trap_sim_factory(dt_star, seed=101, with_trap=False)  # pure noise
    sim.run(10)

    means, variances = [], []
    prev = positions_of(sim)
    for _ in range(M):
        sim.run(1)
        cur = positions_of(sim)
        d = (cur - prev)[:, 0]                 # the mean is NOT subtracted
        means.append(float(d.mean()))
        variances.append(float(d.var(ddof=1)))
        prev = cur

    assert_statistic_fluctuates(means, "mean_i d_i")
    assert_statistic_fluctuates(variances, "Var_i d_i")

    ratio = float(np.var(means, ddof=1)) * TRAP_N / float(np.mean(variances))
    se_ratio = ratio * math.sqrt(2.0 / (M - 1))

    assert sigma_away(ratio, 1.0, se_ratio) < N_SIGMA, (
        f"independence ratio {ratio:.4f}±{se_ratio:.4f} != 1 → effective sample "
        f"{TRAP_N/ratio:.0f}/{TRAP_N}")


# =============================================================================
# B9 — the noise is uniform, not Gaussian
# =============================================================================
@pytest.mark.benchmark
def test_B9_noise_is_uniform_not_gaussian(trap_sim_factory, positions_of):
    """★ Kurtosis 1.80 (uniform) vs 3.00 (Gaussian), max/sigma = sqrt(3).

    The second moment is unaffected, but **the tails are cut off, so barrier-crossing
    rates come out wrong.** Forget this and an escape-rate calculation is quietly
    wrong.
    """
    dt_star, n = 1e-3, 20000
    sim, _ = trap_sim_factory(dt_star, seed=99, n=n, with_trap=False)
    q0 = positions_of(sim)[:, 0]
    sim.run(1)
    step = positions_of(sim)[:, 0] - q0

    sd = float(step.std())
    # the variance has to match fluctuation-dissipation: sqrt(2 D* dt*), D*=1
    assert sd == pytest.approx(math.sqrt(2 * dt_star), rel=5e-3)

    kurt = float(np.mean((step / sd) ** 4))
    assert kurt == pytest.approx(1.80, abs=0.05), \
        f"kurtosis {kurt:.4f} — a uniform distribution gives 1.80"
    assert abs(kurt - 3.0) > 1.0, "the Gaussian hypothesis (3.00) has to be rejected"

    max_over_sd = float(np.abs(step).max() / sd)
    assert max_over_sd == pytest.approx(math.sqrt(3), rel=0.02), (
        f"max/sigma = {max_over_sd:.4f} — a uniform distribution's structural upper "
        f"bound is sqrt(3)=1.7321")


# =============================================================================
# Reproducibility — the premise of the statistical tests
# =============================================================================
def test_same_seed_reproduces_bitwise(trap_sim_factory, positions_of):
    """★ Why fixed-seed statistical tests are legitimate. Break this and every test
    above becomes flaky."""
    outs = []
    for _ in range(2):
        sim, _ = trap_sim_factory(5e-3, seed=42, n=200)
        sim.run(500)
        outs.append(positions_of(sim))
    assert np.array_equal(outs[0], outs[1]), \
        "the same seed gave different results"


def test_different_seeds_diverge(trap_sim_factory, positions_of):
    """Does the seed actually do anything — if not, "4 independent seeds" is false."""
    outs = []
    for seed in (42, 43):
        sim, _ = trap_sim_factory(5e-3, seed=seed, n=200)
        sim.run(500)
        outs.append(positions_of(sim))
    assert not np.allclose(outs[0], outs[1])


# =============================================================================
# Physical consistency — dimensional scaling (B2/B3)
# =============================================================================
@pytest.mark.benchmark
def test_B2_radial_variance_scales_with_dimension(trap_sim_factory, positions_of):
    """<r^2>(3D)/<r^2>(2D) = 1.5 — the discriminator for 01_intake.md ambiguity A1."""
    dt_star = 1.0e-2
    res = {}
    for dim in (2, 3):
        vals = []
        for seed in (21, 22):
            sim, _ = trap_sim_factory(dt_star, seed, dim=dim)
            sim.run(int(20 / dt_star))
            stride = int(2.0 / dt_star)
            for _ in range(25):
                sim.run(stride)
                p = positions_of(sim, dim=dim)
                vals.append(float(np.mean(np.sum(p**2, axis=1))))   # <r^2>
        res[dim] = (float(np.mean(vals)), se_of_mean(vals))

    (r2_2d, se2), (r2_3d, se3) = res[2], res[3]
    ratio = r2_3d / r2_2d
    se_ratio = ratio * math.hypot(se3 / r2_3d, se2 / r2_2d)
    assert sigma_away(ratio, 1.5, se_ratio) < N_SIGMA, (
        f"<r^2> ratio {ratio:.4f}±{se_ratio:.4f} != 1.5")
