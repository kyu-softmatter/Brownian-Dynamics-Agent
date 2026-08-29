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

# --- constants and water properties: the single source of truth is `bdbot.constants` ---
#  ★ Merged 2026-08-29. This table also existed in `bdbot.materials.WATER_VISCOSITY`,
#    and the two copies **had already diverged** -- 0.851 vs 0.85566 mPa*s at 300 K,
#    a 0.545 % gap. 0.545 % in `eta` is 0.545 % in `gamma`, in `D`, and in **every
#    timescale derived from them.**
#    Values, interpolation and the extrapolation flag are unchanged and must stay
#    bit-identical: the sealed S2 documents contain 8.5566e-4. The only thing that
#    changed is that the definition now exists once.
#  ⚠ Importing `bdbot` is not expensive here: `bdbot/__init__` is lazy (PEP 562), so
#    `bdbot.constants` pulls in neither pint nor numpy (0.01 s, measured).
#    **The dependency runs one way only, `simbot -> bdbot`.** Reversed, `bdbot.cli`
#    would start pulling in matplotlib.
from bdbot.constants import (K_B, WATER_ETA_SI as _WATER_ETA_TABLE_SI,
                             WATER_RHO_SI as _WATER_RHO_TABLE_SI,
                             interp_table as _interp,
                             water_density_si, water_viscosity_si)


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
