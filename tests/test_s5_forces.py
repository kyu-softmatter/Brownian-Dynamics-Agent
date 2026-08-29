"""S5 — HarmonicTrap 포스의 정확성. HOOMD 필요하지만 적분은 하지 않는다 (빠름).

힘과 에너지가 서로 정합하는지, 활성축 마스킹이 작동하는지 확인한다.
포스가 틀리면 그 뒤 모든 것이 틀리므로 여기가 첫 관문이다.
"""
from __future__ import annotations

import numpy as np
import pytest

from simbot.forces import HarmonicTrap

K = 1.0


@pytest.fixture
def trap_state(hoomd_mod):
    """알려진 위치에 입자를 놓고 힘을 계산한다 (적분 0스텝)."""
    hoomd = hoomd_mod

    # ⚠ sim 참조를 살려둬야 한다. 버리면 GC 가 detach 시키고
    #   `trap.forces` 접근이 DataAccessError 로 죽는다 (2026-07-28에 겪음).
    keepalive = []

    def build(positions, active_axes=(True, True, False), center=(0.0, 0.0, 0.0),
              k_star=K, box=100.0, two_d=True, want_virial=False):
        pos = np.asarray(positions, dtype=np.float64)
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
        snap = hoomd.Snapshot()
        snap.configuration.box = [box, box, 0.0 if two_d else box, 0, 0, 0]
        snap.particles.N = len(pos)
        snap.particles.types = ["A"]
        snap.particles.position[:] = pos
        snap.particles.mass[:] = 1.0
        snap.particles.moment_inertia[:] = 0.0
        sim.create_state_from_snapshot(snap)
        trap = HarmonicTrap(k_star=k_star, center_star=center, active_axes=active_axes)
        bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                                      default_gamma=1.0)
        sim.operations.integrator = hoomd.md.Integrator(dt=1e-6, methods=[bd],
                                                       forces=[trap])
        if want_virial:
            # virial 은 필요할 때만 계산된다 — 압력을 요구해서 켠다
            sim.operations.computes.append(
                hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All()))
        sim.run(0)          # 힘 계산만
        keepalive.append(sim)
        return sim, trap, pos

    yield build
    keepalive.clear()


def test_force_is_minus_k_times_displacement(trap_state):
    pos = [[1.0, 0.0, 0.0], [0.0, -2.0, 0.0], [3.0, 4.0, 0.0], [0.0, 0.0, 0.0]]
    sim, trap, p = trap_state(pos)
    f = np.array(trap.forces)
    expected = -K * p
    expected[:, 2] = 0.0            # z 축 비활성
    np.testing.assert_allclose(f, expected, rtol=1e-13, atol=1e-15)


def test_potential_energy_is_half_k_r_squared(trap_state):
    pos = [[3.0, 4.0, 0.0], [1.0, 0.0, 0.0], [0.0, 0.0, 0.0]]
    sim, trap, p = trap_state(pos)
    e = np.array(trap.energies)
    r2 = p[:, 0] ** 2 + p[:, 1] ** 2
    np.testing.assert_allclose(e, 0.5 * K * r2, rtol=1e-13, atol=1e-15)
    # (3,4) 는 r=5 이므로 U = 12.5
    assert e[0] == pytest.approx(12.5, rel=1e-13)


def test_force_and_energy_are_mutually_consistent(trap_state):
    """★ F = -dU/dr 을 수치미분으로 확인. 두 코드경로가 같은 포텐셜을 구현하는가."""
    r0, h = 2.3, 1e-6
    _, t_a, _ = trap_state([[r0 - h, 0.0, 0.0]])
    _, t_b, _ = trap_state([[r0 + h, 0.0, 0.0]])
    _, t_m, _ = trap_state([[r0, 0.0, 0.0]])
    dUdx = (np.array(t_b.energies)[0] - np.array(t_a.energies)[0]) / (2 * h)
    fx = np.array(t_m.forces)[0, 0]
    assert fx == pytest.approx(-dUdx, rel=1e-6)


def test_inactive_axis_receives_no_force(trap_state):
    """2D 계에서 z 방향 힘이 정확히 0이어야 한다 (아니면 입자가 평면을 벗어난다)."""
    sim, trap, _ = trap_state([[1.0, 1.0, 0.0]], active_axes=(True, True, False))
    assert np.array(trap.forces)[0, 2] == 0.0


def test_all_three_axes_active_in_3d(trap_state):
    sim, trap, p = trap_state([[1.0, 2.0, 3.0]], active_axes=(True, True, True),
                              two_d=False)
    np.testing.assert_allclose(np.array(trap.forces)[0], -K * p[0], rtol=1e-13)
    assert np.array(trap.energies)[0] == pytest.approx(0.5 * K * 14.0, rel=1e-13)


def test_offset_center_shifts_the_minimum(trap_state):
    """트랩 중심이 (2,0,0) 이면 그 점에서 힘이 0."""
    sim, trap, _ = trap_state([[2.0, 0.0, 0.0], [4.0, 0.0, 0.0]],
                              center=(2.0, 0.0, 0.0))
    f = np.array(trap.forces)
    assert f[0, 0] == pytest.approx(0.0, abs=1e-15)
    assert f[1, 0] == pytest.approx(-2.0 * K, rel=1e-13)


def test_stiffness_scales_force_linearly(trap_state):
    _, t1, _ = trap_state([[1.5, 0.0, 0.0]], k_star=1.0)
    _, t3, _ = trap_state([[1.5, 0.0, 0.0]], k_star=3.0)
    assert np.array(t3.forces)[0, 0] == pytest.approx(
        3 * np.array(t1.forces)[0, 0], rel=1e-13)


def test_virial_trace_relates_to_force_dot_displacement(trap_state):
    """비리얼 대각합 = F . dr = -k r^2 = -2U.

    virial 은 압력이 요구될 때만 계산된다 → want_virial=True 로 켠다.
    """
    sim, trap, p = trap_state([[3.0, 4.0, 0.0]], want_virial=True)
    v = np.asarray(trap.virials)
    if v.ndim != 2 or v.shape[0] != 1:
        pytest.skip(f"이 빌드에서 virial 배열을 얻지 못했다 (shape={v.shape})")
    trace = v[0, 0] + v[0, 3] + v[0, 5]          # xx + yy + zz
    assert trace == pytest.approx(-K * 25.0, rel=1e-13)
    assert trace == pytest.approx(-2 * np.asarray(trap.energies)[0], rel=1e-13)


def test_rejects_malformed_arguments():
    with pytest.raises(ValueError):
        HarmonicTrap(k_star=1.0, center_star=(0.0, 0.0))
    with pytest.raises(ValueError):
        HarmonicTrap(k_star=1.0, active_axes=(True, True))
