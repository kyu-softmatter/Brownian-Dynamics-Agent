"""공용 fixture · 마커 · 통계 도우미.

# 테스트 파일 명명 규약 — 파이프라인 단계별
#   test_s0_units.py        단위·상수·척도 (모든 단계의 기반)
#   test_s2_estimators.py   S2 예측 엔진 (해석해)
#   test_s4_nondim.py       S4 무차원화 왕복        (nondim.py 구현 후)
#   test_s5_forces.py       S5 포스 정확성
#   test_s5_scheme.py       S5 적분기·RNG 벤치마크  [slow]
#   test_s7_analysis.py     S7 분석                (analysis/ 구현 후)
#   test_knowledge.py       knowledge/wiki/benchmarks 회귀

## 통계 테스트 작성 규칙 (읽고 따를 것)

1. **허용오차는 이론 통계오차에서 뽑는다.** 관측값을 보고 오차를 재단하면
   그것은 검증이 아니라 사후합리화다 (master_plan §S2 실패모드).
   `assert abs(measured - predicted) < n_sigma * SE` 형태로 쓴다.

2. **seed 를 고정한다.** HOOMD `Brownian` 은 같은 seed 에서 비트 단위로 재현된다
   (2026-07-28 확인: 최대 절대차 `0.0`). 따라서 고정 seed + 이론 허용오차 조합이
   "재현 가능하면서 정직한" 테스트를 만든다.

3. **통계량이 요동하는지 확인한다.** 요동하지 않는 측정값은 산술 항등식이다.
   `simbot.guards.assert_statistic_fluctuates` 를 쓴다.

4. **경쟁 가설을 함께 기각한다.** "예측과 맞다"만으로는 약하다.
   대안 가설(예: exact 적분 스킴)이 몇 σ 로 기각되는지도 assert 한다.
"""
from __future__ import annotations

import numpy as np
import pytest

# 조화 트랩 카드 축약 단위:  l_trap = 1, kT = 1, tau_trap = 1  =>  k* = 1, D* = 1, gamma* = 1
TRAP_BOX = 200.0        # >> l_trap = 1.  트랩 카드 §7 "박스 >> l_trap" 게이트
TRAP_N = 1000
TRAP_DIM = 2


def pytest_configure(config):
    config.addinivalue_line("markers", "slow: HOOMD 실행이 필요한 테스트 (초 단위)")
    config.addinivalue_line("markers", "benchmark: knowledge/wiki/benchmarks 회귀 항목")


@pytest.fixture(scope="session")
def hoomd_mod():
    """HOOMD 를 세션당 한 번만 import (import 비용이 크다)."""
    return pytest.importorskip("hoomd")


@pytest.fixture
def trap_sim_factory(hoomd_mod):
    """조화 트랩 시뮬레이션 팩토리. 축약 단위, 2D, 비상호작용 N개 입자."""
    hoomd = hoomd_mod
    from simbot.forces import HarmonicTrap

    def make(dt_star: float, seed: int, n: int = TRAP_N, with_trap: bool = True,
             dim: int = TRAP_DIM):
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)
        snap = hoomd.Snapshot()
        # Lz = 0 => 2D.  configuration.dimensions 에는 setter 가 없다 (HOOMD 7)
        lz = 0.0 if dim == 2 else TRAP_BOX
        snap.configuration.box = [TRAP_BOX, TRAP_BOX, lz, 0, 0, 0]
        snap.particles.N = n
        snap.particles.types = ["A"]
        snap.particles.position[:] = 0.0
        snap.particles.typeid[:] = 0
        snap.particles.mass[:] = 1.0
        snap.particles.moment_inertia[:] = 0.0   # 회전 자유도 끔
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
    """시뮬레이션에서 활성 차원의 위치 배열을 꺼낸다."""
    def get(sim, dim: int = TRAP_DIM):
        s = sim.state.get_snapshot()
        return np.array(s.particles.position[:, :dim], dtype=np.float64)
    return get


# --- 통계 도우미 -----------------------------------------------------------
def se_of_mean(samples) -> float:
    """블록 평균들의 표준오차."""
    s = np.asarray(samples, dtype=np.float64)
    return float(np.std(s, ddof=1) / np.sqrt(s.size))


def sigma_away(measured: float, predicted: float, se: float) -> float:
    """예측에서 몇 σ 떨어졌는가."""
    if se <= 0:
        raise ValueError("se must be > 0 — 요동하지 않는 통계량은 항등식이다")
    return abs(measured - predicted) / se
