"""S5 — 런타임 가드. 순수 함수는 빠르고, 배위 온도(B8)만 HOOMD 를 쓴다.

가드 테스트의 원칙: **통과만 테스트하면 가드가 죽어도 모른다.**
발동해야 하는 경우도 반드시 함께 테스트한다.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from simbot.guards import (
    assert_statistic_fluctuates, check_finite, check_inside_box,
    check_step_displacements, configurational_temperature,
)

from .conftest import TRAP_DIM, se_of_mean, sigma_away


# =============================================================================
# 배위 온도 — 순수 함수
# =============================================================================
def test_configurational_temperature_recovers_kT_analytically():
    """조화 트랩에서 kT_conf = k^2<r^2>/(d k) = k<r^2>/d.

    d=2, k=1, <r^2> = 2 (즉 <x^2>=1) 이면 kT_conf = 1.
    """
    rng = np.random.default_rng(0)
    n, d, k, kT = 200000, 2, 1.0, 1.0
    # Boltzmann 분포: 성분별 표준편차 sqrt(kT/k)
    pos = rng.normal(0.0, math.sqrt(kT / k), size=(n, d))
    forces = -k * pos
    kT_conf = configurational_temperature(forces, laplacian_U_total=d * k)
    assert kT_conf == pytest.approx(kT, rel=0.01)


@pytest.mark.parametrize("kT", [0.5, 1.0, 2.5])
def test_configurational_temperature_tracks_input_kT(kT):
    rng = np.random.default_rng(1)
    n, d, k = 200000, 3, 2.0
    pos = rng.normal(0.0, math.sqrt(kT / k), size=(n, d))
    kT_conf = configurational_temperature(-k * pos, laplacian_U_total=d * k)
    assert kT_conf == pytest.approx(kT, rel=0.01)


def test_configurational_temperature_rejects_bad_laplacian():
    with pytest.raises(ValueError):
        configurational_temperature(np.ones((3, 2)), laplacian_U_total=0.0)
    with pytest.raises(ValueError):
        configurational_temperature(np.ones((3, 2)), laplacian_U_total=-1.0)


# =============================================================================
# 스텝당 변위
# =============================================================================
def test_displacement_guard_passes_for_small_steps():
    dr = np.full((100, 3), 0.001)
    r = check_step_displacements(dr, sigma=1.0, max_frac=0.10)
    assert r.passed and r.n_exceeding == 0 and r.note == ""


def test_displacement_guard_fires_on_explosion():
    """★ 가드가 실제로 발동하는가."""
    dr = np.zeros((100, 3))
    dr[7] = [0.5, 0.0, 0.0]              # 0.5 sigma — 폭발 징후
    r = check_step_displacements(dr, sigma=1.0, max_frac=0.10)
    assert not r.passed
    assert r.n_exceeding == 1
    assert "dt 과대" in r.note
    assert r.max_over_sigma == pytest.approx(0.5)


def test_displacement_guard_reports_max_over_rms():
    """균일분포 노이즈에서 max/rms 가 sqrt(3) 근처인지 볼 수 있어야 한다."""
    rng = np.random.default_rng(3)
    dr = rng.uniform(-1, 1, size=(50000, 1)) * 0.01
    dr = np.hstack([dr, np.zeros((len(dr), 2))])
    r = check_step_displacements(dr, sigma=1.0, max_frac=1.0)
    assert r.passed
    assert r.max_over_rms == pytest.approx(math.sqrt(3), rel=0.02)


# =============================================================================
# 유한성 · 경계
# =============================================================================
def test_check_finite_catches_nan_and_inf():
    ok, fails = check_finite(pos=np.ones((3, 3)), force=np.ones((3, 3)))
    assert ok and not fails

    bad = np.ones((3, 3))
    bad[1, 1] = np.nan
    worse = np.ones((3, 3))
    worse[0, 0] = np.inf
    ok, fails = check_finite(pos=bad, force=worse)
    assert not ok
    assert any("pos" in f for f in fails) and any("force" in f for f in fails)


def test_check_inside_box():
    ok, n = check_inside_box(np.array([[1.0, 1.0, 0.0]]), [10, 10, 10], dims=2)
    assert ok and n == 0
    ok, n = check_inside_box(np.array([[6.0, 0.0, 0.0], [0.0, -7.0, 0.0]]),
                             [10, 10, 10], dims=2)
    assert not ok and n == 2


# =============================================================================
# "항등식을 측정으로 착각하는" 실패 방지 장치
# =============================================================================
def test_fluctuation_check_passes_for_real_statistic():
    rng = np.random.default_rng(5)
    assert_statistic_fluctuates(rng.normal(1.0, 0.1, 500), "real")


def test_fluctuation_check_catches_arithmetic_identity():
    """★ 2026-07-28 에 실제로 겪은 실패를 테스트로 고정한다.

    변위에서 평균을 뺀 뒤 교차상관을 재면 cross/auto = -1/(n-1) 이 항등적으로
    나오고, 200회 반복의 표준편차가 6.7e-20 이었다. 결과가 그럴듯해서 통과할 뻔했다.
    """
    n = 1000
    identity_values = np.full(300, -1.0 / (n - 1))     # 요동하지 않는 "측정값"
    with pytest.raises(AssertionError, match="산술 항등식"):
        assert_statistic_fluctuates(identity_values, "cross/auto")


def test_fluctuation_check_reproduces_the_original_bug(hoomd_mod):
    """평균을 빼면 정말로 -1/(n-1) 이 나오는지 직접 확인한다 (HOOMD 불필요, 순수 산술)."""
    rng = np.random.default_rng(7)
    n = 1000
    seen = []
    for _ in range(50):
        d = rng.normal(0.0, 1.0, n)
        d = d - d.mean()                              # ← 이 한 줄이 측정을 파괴한다
        cross = (d.sum() ** 2 - np.sum(d**2)) / (n * (n - 1))
        seen.append(cross / np.mean(d**2))
    seen = np.array(seen)
    np.testing.assert_allclose(seen, -1.0 / (n - 1), rtol=1e-10)
    assert seen.std() < 1e-15, "항등식이므로 요동이 없어야 한다"


# =============================================================================
# B8 — 배위 온도가 HOOMD 시뮬레이션에서 kT 를 복원하는가  [slow]
# =============================================================================
@pytest.mark.slow
@pytest.mark.benchmark
def test_B8_configurational_temperature_recovers_kT_in_simulation(
        trap_sim_factory, positions_of):
    """★ BD 의 진짜 온도계. 운동에너지 온도는 매 스텝 뽑히므로 쓸 수 없다.

    조화 트랩: kT_conf = <|grad U|^2>/<lap U> = k^2<r^2>/(d k)
    """
    dt_star, k = 5e-3, 1.0
    sim, _ = trap_sim_factory(dt_star, seed=303)
    sim.run(int(20 / dt_star))

    stride = int(2.0 / dt_star)
    blocks = []
    for _ in range(40):
        sim.run(stride)
        p = positions_of(sim)
        blocks.append(configurational_temperature(k * p, laplacian_U_total=TRAP_DIM * k))

    assert_statistic_fluctuates(blocks, "kT_conf blocks")
    kT_conf = float(np.mean(blocks))
    se = se_of_mean(blocks)

    assert sigma_away(kT_conf, 1.0, se) < 4.0, (
        f"kT_conf = {kT_conf:.5f}±{se:.5f} 가 입력 kT*=1.0 에서 벗어났다")


@pytest.mark.slow
def test_kinetic_temperature_cannot_deviate_systematically(trap_sim_factory, hoomd_mod):
    """★ 가드로 쓸 수 없다는 것을 테스트로 고정한다.

    HOOMD `Brownian` 은 속도를 적분하지 않고 매 스텝 목표 분포에서 **뽑는다.**
    따라서 kinetic_temperature 는 dt* 를 10배 키워 적분이 엉망이 되어도
    여전히 kT 근처에 머문다 — 계통 이탈을 감지할 수 없다.
    """
    hoomd = hoomd_mod
    temps = {}
    for dt_star in (1e-3, 1e-1):          # 100배 차이. 후자는 적분 정확도가 나쁘다
        sim, _ = trap_sim_factory(dt_star, seed=404)
        tq = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
        sim.operations.computes.append(tq)
        sim.run(200)
        temps[dt_star] = tq.kinetic_temperature

    # dt* 를 100배 키워도 운동온도는 kT=1 근처에 머문다 → 진단 능력 없음
    for dt_star, T in temps.items():
        assert T == pytest.approx(1.0, rel=0.1), (
            f"dt*={dt_star:g} 에서 운동온도 {T:.4f}")
