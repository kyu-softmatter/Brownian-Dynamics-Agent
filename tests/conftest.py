"""Shared fixtures · markers · statistics helpers.

# Test-file naming convention — by pipeline stage
#   test_s0_units.py        units, constants, scales (the base of every stage)
#   test_s2_estimators.py   the S2 prediction engine (analytic solutions)
#   test_s4_nondim.py       S4 non-dimensionalization round-trip (after nondim.py)
#   test_s5_forces.py       S5 force correctness
#   test_s5_scheme.py       S5 integrator and RNG benchmarks  [slow]
#   test_s7_analysis.py     S7 analysis                (after analysis/ exists)
#   test_knowledge.py       knowledge/wiki/benchmarks regression

## Rules for writing a statistical test (read them and follow them)

1. **Take the tolerance from the theoretical statistical error.** Cutting the
   tolerance to fit the observed value is not verification but post-hoc
   rationalisation (master_plan §S2 failure modes). Write it in the form
   `assert abs(measured - predicted) < n_sigma * SE`.

2. **Fix the seed.** HOOMD `Brownian` reproduces bit-for-bit on the same seed
   (confirmed 2026-07-28: max absolute difference `0.0`). So a fixed seed plus a
   theoretical tolerance is what makes a test "reproducible and honest" at once.

3. **Check that the statistic fluctuates.** A measurement that does not fluctuate
   is an arithmetic identity. Use
   `simbot.guards.assert_statistic_fluctuates`.

4. **Reject the competing hypothesis too.** "It agrees with the prediction" alone
   is weak. Also assert how many σ the alternative (an exact integration scheme,
   say) is rejected by.
"""
from __future__ import annotations

import numpy as np
import pytest

# Harmonic-trap card reduced units:  l_trap = 1, kT = 1, tau_trap = 1
#   =>  k* = 1, D* = 1, gamma* = 1
TRAP_BOX = 200.0        # >> l_trap = 1.  trap card §7 "box >> l_trap" gate
TRAP_N = 1000
TRAP_DIM = 2


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: needs a HOOMD run (seconds)")
    config.addinivalue_line("markers",
                            "benchmark: a knowledge/wiki/benchmarks regression item")


@pytest.fixture(scope="session")
def hoomd_mod():
    """Import HOOMD once per session (the import is expensive)."""
    return pytest.importorskip("hoomd")


@pytest.fixture
def trap_sim_factory(hoomd_mod):
    """Harmonic-trap simulation factory. Reduced units, 2D, N non-interacting."""
    hoomd = hoomd_mod
    from simbot.forces import HarmonicTrap

    def make(dt_star: float, seed: int, n: int = TRAP_N, with_trap: bool = True,
             dim: int = TRAP_DIM):
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)
        snap = hoomd.Snapshot()
        # Lz = 0 => 2D.  configuration.dimensions has no setter (HOOMD 7)
        lz = 0.0 if dim == 2 else TRAP_BOX
        snap.configuration.box = [TRAP_BOX, TRAP_BOX, lz, 0, 0, 0]
        snap.particles.N = n
        snap.particles.types = ["A"]
        snap.particles.position[:] = 0.0
        snap.particles.typeid[:] = 0
        snap.particles.mass[:] = 1.0
        snap.particles.moment_inertia[:] = 0.0   # rotational DOF off
        sim.create_state_from_snapshot(snap)

        trap = None
        if with_trap:
            axes = (True, True, dim == 3)
            trap = HarmonicTrap(k_star=1.0, active_axes=axes)
        bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                                      default_gamma=1.0)
        sim.operations.integrator = hoomd.md.Integrator(
            dt=dt_star, methods=[bd], forces=([trap] if trap else []))
        return sim, trap

    return make


@pytest.fixture
def positions_of():
    """Pull the active-dimension position array out of a simulation."""
    def get(sim, dim: int = TRAP_DIM):
        s = sim.state.get_snapshot()
        return np.array(s.particles.position[:, :dim], dtype=np.float64)
    return get


# --- statistics helpers ----------------------------------------------------
def se_of_mean(samples) -> float:
    """Standard error of the block means."""
    s = np.asarray(samples, dtype=np.float64)
    return float(np.std(s, ddof=1) / np.sqrt(s.size))


def sigma_away(measured: float, predicted: float, se: float) -> float:
    """How many σ away from the prediction."""
    if se <= 0:
        raise ValueError("se must be > 0 — a statistic that does not fluctuate "
                         "is an identity")
    return abs(measured - predicted) / se
