"""물성 유도식 — 전부 차원 있는 값 (skill `bd-physics` §2).

1-A·1-B가 똑같이 계산한 것만 넣었습니다: γ, D_t, τ_B, m, τ_p.

⚠️ 구에만 성립하는 것과 형상 무관한 것을 섞지 마세요. `sphere_*` 이름이 붙은 것은
   구 전용입니다. 타원체·막대는 Perrin 인자가 필요하고, **BD에서는 병진 이방성이
   애초에 재현되지 않습니다** (skill `bd-hoomd` 하드 제약).
"""
from __future__ import annotations

import math

from .units import Q, kB


def thermal_energy(T):
    return (kB * T).to("J")


def sphere_drag(eta, d):
    """Stokes 항력 계수 γ = 3πηd (구)."""
    return (3 * math.pi * eta * d).to("kg/s")


def diffusion(kT, gamma):
    """Stokes–Einstein D_t = kT/γ. 형상 무관 (γ만 맞으면 됨)."""
    return (kT / gamma).to("m^2/s")


def brownian_time(d, D_t):
    """확산 시간 τ_B = d²/D_t. 기준 시간으로 쓰는 값."""
    return (d**2 / D_t).to("s")


def sphere_mass(rho_p, d):
    return (rho_p * (math.pi / 6) * d**3).to("kg")


def momentum_time(m, gamma):
    """관성 이완 τ_p = m/γ.

    ★ `dt`와 비교하지 마세요. BD는 관성이 아예 없습니다 (skill `bd-physics` §4).
      `τ_p`가 답하는 질문은 "BD를 써도 되는 계인가"이고, 비교 대상은 그 계의
      **관심 최속 시간척도** `τ_dyn` 입니다.
    """
    return (m / gamma).to("s")


def sphere_rotational_diffusion(kT, eta, d):
    """Stokes–Einstein–Debye D_r = kT/(πηd³) = 3D_t/d².  ⚠️ 구에만 성립."""
    return (kT / (math.pi * eta * d**3)).to("1/s")


def sphere_bulk(d, T, eta, rho_p=None) -> dict:
    """구 + 뉴턴 유체의 기본 물성을 한 번에. 두 케이스가 똑같이 이 묶음을 씁니다."""
    d = d.to("m")
    kT = thermal_energy(T)
    gamma = sphere_drag(eta.to("Pa*s"), d)
    D_t = diffusion(kT, gamma)
    out = {"kT": kT, "gamma": gamma, "D_t": D_t, "tau_B": brownian_time(d, D_t), "d": d}
    if rho_p is not None:
        m = sphere_mass(rho_p.to("kg/m^3"), d)
        out["m"] = m
        out["tau_p"] = momentum_time(m, gamma)
    return out


# ── 물성 상수 (핸드북, tier 0) ─────────────────────────────────────────
WATER_VISCOSITY = {300: Q(0.851, "mPa*s"), 298: Q(0.890, "mPa*s"), 293: Q(1.002, "mPa*s")}
DENSITY = {"silica": Q(2000, "kg/m^3"), "polystyrene": Q(1050, "kg/m^3")}

__all__ = ["thermal_energy", "sphere_drag", "diffusion", "brownian_time", "sphere_mass",
           "momentum_time", "sphere_rotational_diffusion", "sphere_bulk",
           "WATER_VISCOSITY", "DENSITY"]
