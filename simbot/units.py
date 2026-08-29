"""물리 상수와 SI ↔ 무차원 변환.

규약 (CLAUDE.md):
  - `*_si`   : SI 단위 물리량 (float, 단위는 이름/문서로 고정)
  - `*_star` : 무차원 값 (순수 float)
  두 접미사를 섞은 산술은 버그다. `tests/test_units.py`가 감시한다.

무차원화의 기준 척도는 **(계 × 목적동역학) 카드**가 소유한다.
따라서 이 모듈은 보편 규약을 강요하지 않고, `Scales`를 인자로 받는다.
  → knowledge/wiki/systems/_index.md
"""
from __future__ import annotations

from dataclasses import dataclass
import math

# --- 상수 (SI 2019 정의값) ---------------------------------------------------
K_B = 1.380649e-23  # J/K, 정확값

# --- 물 물성 (IAPWS). 출처: knowledge/wiki/concepts/water-298k.md ------------
_WATER_ETA_TABLE_SI = {  # T[K] -> eta[Pa*s]
    293.15: 1.0016e-3,
    298.15: 0.8900e-3,
    303.15: 0.7972e-3,
    308.15: 0.7191e-3,
}
_WATER_RHO_TABLE_SI = {293.15: 998.21, 298.15: 997.05, 303.15: 995.65, 308.15: 994.03}


def _interp(table: dict[float, float], T_si: float) -> tuple[float, bool]:
    """표에서 선형 보간. (값, 외삽했는가)를 돌려준다."""
    ts = sorted(table)
    if T_si in table:
        return table[T_si], False
    if T_si < ts[0] or T_si > ts[-1]:
        # 외삽 — 호출자가 provenance를 낮춰야 한다
        lo, hi = (ts[0], ts[1]) if T_si < ts[0] else (ts[-2], ts[-1])
        f = (T_si - lo) / (hi - lo)
        return table[lo] + f * (table[hi] - table[lo]), True
    for lo, hi in zip(ts, ts[1:]):
        if lo <= T_si <= hi:
            f = (T_si - lo) / (hi - lo)
            return table[lo] + f * (table[hi] - table[lo]), False
    raise AssertionError("unreachable")


def water_viscosity_si(T_si: float) -> tuple[float, bool]:
    """물의 점도 [Pa*s]. (값, 외삽여부). 293-308 K 밖은 외삽이므로 신뢰도를 낮출 것."""
    return _interp(_WATER_ETA_TABLE_SI, T_si)


def water_density_si(T_si: float) -> tuple[float, bool]:
    """물의 밀도 [kg/m^3]. (값, 외삽여부)."""
    return _interp(_WATER_RHO_TABLE_SI, T_si)


# --- 기본 관계식 ------------------------------------------------------------
def kT_si(T_si: float) -> float:
    """열에너지 [J]."""
    return K_B * T_si


def stokes_drag_si(eta_si: float, radius_si: float) -> float:
    """구의 Stokes 항력계수 gamma = 6*pi*eta*a  [kg/s].

    ⚠ 인자는 **반지름**이다. 직경을 넣으면 gamma가 2배, D가 절반이 되어
       모든 시간척도가 2배 틀린다. knowledge/wiki/concepts/water-298k.md 참조.
    """
    return 6.0 * math.pi * eta_si * radius_si


def stokes_einstein_D_si(T_si: float, gamma_si: float) -> float:
    """병진 확산계수 D = kT/gamma  [m^2/s]."""
    return kT_si(T_si) / gamma_si


# --- 무차원화 기준 척도 -----------------------------------------------------
@dataclass(frozen=True)
class Scales:
    """무차원화 기준 3개. 어느 값을 고를지는 (계 × 목적동역학) 카드가 정한다.

    length_si : 기준 길이 [m]
    energy_si : 기준 에너지 [J]
    time_si   : 기준 시간 [s]
    origin    : 이 선택의 출처 (카드 경로 등)
    """

    length_si: float
    energy_si: float
    time_si: float
    origin: str = ""

    # 파생 척도
    @property
    def force_si(self) -> float:
        return self.energy_si / self.length_si

    @property
    def stiffness_si(self) -> float:
        """스프링 상수 척도 [N/m]."""
        return self.energy_si / self.length_si**2

    @property
    def velocity_si(self) -> float:
        return self.length_si / self.time_si

    @property
    def diffusivity_si(self) -> float:
        return self.length_si**2 / self.time_si

    @property
    def rate_si(self) -> float:
        """각주파수·회전확산 등 [1/s]."""
        return 1.0 / self.time_si

    @property
    def modulus_3d_si(self) -> float:
        """탄성률 척도 [Pa] = energy/length^3."""
        return self.energy_si / self.length_si**3

    # --- 변환 ---
    def to_star(self, value_si: float, kind: str) -> float:
        return value_si / self._scale_for(kind)

    def to_si(self, value_star: float, kind: str) -> float:
        return value_star * self._scale_for(kind)

    def _scale_for(self, kind: str) -> float:
        table = {
            "length": self.length_si,
            "energy": self.energy_si,
            "time": self.time_si,
            "force": self.force_si,
            "stiffness": self.stiffness_si,
            "velocity": self.velocity_si,
            "diffusivity": self.diffusivity_si,
            "rate": self.rate_si,
            "modulus_3d": self.modulus_3d_si,
            "area": self.length_si**2,
            "volume": self.length_si**3,
        }
        if kind not in table:
            raise KeyError(f"unknown scale kind {kind!r}; known: {sorted(table)}")
        return table[kind]


# --- 카드별 기준 척도 팩토리 -----------------------------------------------
def scales_brownian(sigma_si: float, T_si: float, gamma_si: float) -> Scales:
    """수동 구형 × 수송:  (sigma, kT, tau_D = sigma^2/D0).

    카드: knowledge/wiki/systems/passive-sphere--*.md
    """
    D0 = stokes_einstein_D_si(T_si, gamma_si)
    return Scales(
        length_si=sigma_si,
        energy_si=kT_si(T_si),
        time_si=sigma_si**2 / D0,
        origin="scales_brownian: (sigma, kT, tau_D)",
    )


def scales_soft2d(d_si: float, sigma_si: float, T_si: float,
                  gamma_si: float | None = None) -> Scales:
    """2D 소프트 반발계 × 평형 구조:  (`d = n^{-1/2}`, `kT`, `τ_d = d²/D₀`).

    카드: knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md §3

    ★ **길이 척도와 입자 크기가 다른 유일한 카드다.** 길이 단위는 격자 간격
      `d = n^{-1/2}` 이고, 항력은 입자 직경 `σ` 가 정한다:
      `γ = 6πη(σ/2) = 3πησ`. 둘을 같다고 두면 `τ_d` 가 `(d/σ)²` 배 틀린다 —
      이 스윕의 `d/σ = 3` 에서 **9배**다.

    ⚠ `σ` 는 **동역학에 들어가지 않는다** (경질 코어가 없다). 시간 척도와
      물리적 타당성 판단에만 쓰인다 → `build.coverage_from_sigma_over_d`.
    """
    if gamma_si is None:
        eta_si, extrapolated = water_viscosity_si(T_si)
        if extrapolated:
            raise ValueError(
                f"T = {T_si} K 는 물 점도 표(293–308 K) 밖이라 외삽이 된다. "
                f"gamma_si 를 명시하라 — 외삽값을 조용히 쓰면 provenance 가 "
                f"'derived' 인 척하는 'assumed' 가 된다")
        gamma_si = stokes_drag_si(eta_si, sigma_si / 2.0)
    D0 = stokes_einstein_D_si(T_si, gamma_si)
    return Scales(
        length_si=d_si,
        energy_si=kT_si(T_si),
        time_si=d_si**2 / D0,
        origin=f"scales_soft2d: (d={d_si:.4g} m, kT, tau_d=d^2/D0), "
               f"gamma from sigma={sigma_si:.4g} m",
    )


def scales_harmonic_trap(k_si: float, T_si: float, gamma_si: float) -> Scales:
    """조화 트랩:  (l_trap = sqrt(kT/k), kT, tau_trap = gamma/k).

    이 선택 하에 무차원 운동방정식이 `dr*/dt* = -r* + sqrt(2) xi` 로 정규화되고
    `D* = 1`, `k* = 1` 이 된다.

    ★ tau_D 를 강요하면 안 되는 이유: tau_trap/tau_D = kT/(k sigma^2) = 1/k*_sigma.
      강한 트랩(k*_sigma >> 1)에서 tau_D 기준 dt 는 완화시간보다 커진다.
      카드: knowledge/wiki/systems/passive-sphere--harmonic-trap.md
    """
    return Scales(
        length_si=math.sqrt(kT_si(T_si) / k_si),
        energy_si=kT_si(T_si),
        time_si=gamma_si / k_si,
        origin="scales_harmonic_trap: (l_trap, kT, tau_trap)",
    )
