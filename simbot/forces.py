"""외장 포텐셜 — HOOMD 내장으로 표현할 수 없는 것들.

무차원 규약: 모든 파라미터는 해당 (계 × 목적동역학) 카드의 축약 단위(`*`)로 받는다.
조화 트랩 카드에서는 `l_trap` 기준이므로 `k_star = 1`, `center_star`는 `l_trap` 단위다.
"""
from __future__ import annotations

import numpy as np
import hoomd


class HarmonicTrap(hoomd.md.force.Custom):
    """등방 조화 트랩:  U = 1/2 k |r - r0|^2,  F = -k (r - r0)

    HOOMD 7.1에는 조화 구속 내장 포텐셜이 없다 (`md.external.field`는 Periodic·
    Electric·Magnetic 뿐). ghost 입자 + `md.bond.Harmonic` 으로도 가능하지만
    Custom 이 더 직접적이고 검사하기 쉽다 (사용자 결정 2026-07-28).

    Args:
        k_star:  스프링 상수 (축약 단위). 트랩 카드 규약에서는 1.0
        center_star: 트랩 중심 (축약 길이 단위)
        active_axes: 트랩이 작용하는 축. 2D 계는 (True, True, False)

    ⚠ **박스가 구속 길이보다 훨씬 커야 한다.** 여기서는 래핑된 좌표를 그대로 쓴다 —
      입자가 경계에 가면 트랩이 잘못된 방향으로 당긴다. 조화 트랩 카드 §7의
      "박스 >> l_trap" 게이트가 이것을 보장한다 (l_trap/sigma ~ 1e-3 이라 자동 만족).
    """

    def __init__(self, k_star: float, center_star=(0.0, 0.0, 0.0),
                 active_axes=(True, True, True)):
        super().__init__()
        self.k_star = float(k_star)
        self._center = np.asarray(center_star, dtype=np.float64)
        self._mask = np.asarray(active_axes, dtype=bool)
        if self._center.shape != (3,) or self._mask.shape != (3,):
            raise ValueError("center_star and active_axes must have length 3")

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            # 래핑된 좌표. 위 경고 참조.
            dr = np.array(snap.particles.position, dtype=np.float64) - self._center
            dr[:, ~self._mask] = 0.0
            with self.cpu_local_force_arrays as arrays:
                arrays.force[:] = -self.k_star * dr
                # U = 1/2 k |dr|^2  (입자당)
                arrays.potential_energy[:] = 0.5 * self.k_star * np.sum(dr * dr, axis=1)
                # 비리얼: W_ab = F_a * dr_b  (외장 포텐셜의 기여)
                # 순서 규약: xx, xy, xz, yy, yz, zz
                #
                # ⚠ **미검증.** `Force.virials` 가 이 경로에서 노출되지 않아
                #   테스트가 skip 된다 (tests/test_s5_forces.py). 게다가 외장 장의
                #   비리얼은 규약이 모호하다 — 트랩은 계의 운동량 유속의 일부가 아니다.
                #   단일 구속 입자에서는 압력이 의미 있는 관측량이 아니므로 파이프라인에
                #   영향은 없지만, **압력을 로깅하려면 먼저 이 부분을 검증해야 한다.**
                #   → knowledge/wiki/questions/ 항목으로 등재
                f = -self.k_star * dr
                arrays.virial[:, 0] = f[:, 0] * dr[:, 0]
                arrays.virial[:, 1] = f[:, 0] * dr[:, 1]
                arrays.virial[:, 2] = f[:, 0] * dr[:, 2]
                arrays.virial[:, 3] = f[:, 1] * dr[:, 1]
                arrays.virial[:, 4] = f[:, 1] * dr[:, 2]
                arrays.virial[:, 5] = f[:, 2] * dr[:, 2]


# =============================================================================
# 거듭제곱 반발 쌍 포텐셜 — `md.pair.Table` 로 구현한다
# =============================================================================
#  ★ 왜 Table 인가
#  `U = A/r^n` 은 HOOMD 내장에 없다. 확인한 대안들:
#    · `md.pair.Mie` — `U ∝ ε[(σ/r)^n − (σ/r)^m]`. 순수 반발은 `m → 0` 인데
#      계수 `(n/(n−m))(n/m)^{m/(n−m)}` 가 `m=0` 에서 발산한다. 쓸 수 없다
#    · `md.pair.Yukawa` — `exp(−κr)/r` 로 `r^{-1}` 뿐
#    · `md.pair.Table` — 임의 `U(r)`·`F(r)` 를 표로 준다. **이것이 맞다**
#  Table 은 선형 보간이므로 표본 수가 정확도를 정한다 → 해석식과 대조하는
#  테스트가 필수다 (`tests/test_s5_pair.py`).
def power_law_table(hoomd, *, amplitude: float, exponent: float = 3.0,
                    r_cut: float, r_min: float = 0.2, n_points: int = 4096,
                    shift: bool = True, nlist=None, buffer: float = 0.1):
    """`U(r) = A/r^p` 반발을 `md.pair.Table` 로 만든다 (축약 단위).

    Args:
        amplitude: `A`. 축약 단위에서 `A = βU(r=1)` 이다.
        exponent: `p`. 상자성 콜로이드의 쌍극자 반발은 `p = 3`.
        r_cut: 절단 거리. **`≤ L/2` 여야 한다** (최소이미지). 호출자가 보장한다.
        r_min: 표의 시작. **이보다 가까워지면 HOOMD 가 표를 벗어난다** —
            `guards.check_min_separation` 으로 감시할 것.
        shift: `U(r_cut) = 0` 이 되도록 상수를 뺀다. 에너지 불연속을 없앤다.
            ⚠ **힘의 불연속은 남는다** (상수를 빼도 미분은 안 바뀐다).
            절단 오차의 크기는 `pair_truncation_error()` 로 계산한다.
        buffer: Cell 리스트 버퍼 (**절대 거리**, `d` 단위).
            ★ **HOOMD 는 `r_cut + buffer ≤ L/2` 를 요구한다** — `r_cut` 만이 아니다.
            넘으면 런타임에 거부한다:
            `nlist: Simulation box is too small, the neighbor list is searching
            beyond the minimum image` (2026-07-28 실측).
            그래서 버퍼를 `r_cut` 의 비율로 잡으면 안 된다 — 큰 `r_cut` 에서 폭발한다.

    Returns:
        `(pair, info)` — `info` 에 `u_shift` 와 절단 오차가 담긴다.
        `u_shift` 를 기록하는 이유: 총 포텐셜 에너지를 참값으로 되돌리려면
        `U_true = U_hoomd + n_pairs_within_rcut × u_shift` 가 필요하다.
    """
    if r_cut <= r_min:
        raise ValueError(f"r_cut({r_cut}) <= r_min({r_min})")
    nlist = nlist or hoomd.md.nlist.Cell(buffer=buffer)
    r = np.linspace(r_min, r_cut, n_points)
    u = amplitude / r**exponent
    f = exponent * amplitude / r**(exponent + 1)      # F = -dU/dr > 0 (반발)
    u_shift = float(amplitude / r_cut**exponent)
    if shift:
        u = u - u_shift

    pair = hoomd.md.pair.Table(nlist=nlist, default_r_cut=r_cut)
    pair.params[("A", "A")] = dict(r_min=r_min, U=u, F=f)
    info = {
        "potential": f"{amplitude:g}/r^{exponent:g}",
        "amplitude": amplitude, "exponent": exponent,
        "r_cut": r_cut, "r_min": r_min, "n_points": n_points,
        "shifted": shift, "u_shift": u_shift, "nlist_buffer": buffer,
        "r_cut_plus_buffer": r_cut + buffer,
        **pair_truncation_error(amplitude=amplitude, exponent=exponent, r_cut=r_cut),
    }
    return pair, info


def pair_truncation_error(*, amplitude: float, exponent: float, r_cut: float,
                          r_ref: float = 1.0) -> dict:
    """절단이 남기는 오차. **절대값과 상대값을 함께 준다.**

    ★ `simbot.cutoff` 의 허용치는 `kT` 기준(`βU ≤ 0.02`)이다. 그것은 `kT` 가
      지배 척도일 때 옳다. 그런데 `Γ = π^{3/2}A = 557` 인 결정에서는 **`A·kT` 가
      지배 척도**이고, 같은 절단이 `kT` 기준으로는 크고 `A` 기준으로는 작다.
      둘을 함께 보고해서 어느 기준으로 판단할지 사람이 고르게 한다.
    """
    u_cut = amplitude / r_cut**exponent
    f_cut = exponent * amplitude / r_cut**(exponent + 1)
    u_ref = amplitude / r_ref**exponent
    f_ref = exponent * amplitude / r_ref**(exponent + 1)
    return {
        "beta_u_at_rcut": u_cut,                 # kT 기준 (cutoff.py 프리셋과 비교)
        "beta_f_at_rcut": f_cut,
        "u_rel_to_nearest": u_cut / u_ref,       # 상호작용 척도 기준
        "f_rel_to_nearest": f_cut / f_ref,
    }
