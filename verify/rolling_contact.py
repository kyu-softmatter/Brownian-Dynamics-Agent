"""구름 저항 + 접선 스프링 접촉 모델 — `md.force.Custom(aniso=True)` 구현 + 넘파이 참조 모델.

`md.pair.friction`(산일성·구름 면제, `scratch/verify_pair_friction.py`)으로는 굽힘강성을
만들 수 없다는 것이 실측됐다. JKR 접착접촉이 실제로 하는 일은 **접촉점이 표면 위를
구르는 것에 저항하는 것**이므로, 그걸 직접 구현한다 (Dominik–Tielens 1995 / Wada 2007 계열).

━━ 운동학 ━━
결합 (i,j) 마다 **물체고정 접촉점 마커** `m_i`(입자 i 중심 → 접촉점, 단위벡터)를 결합
생성 시점에 저장한다. 현재 랩프레임 값은 `a_i = Rot(q_i)·m_i`.
`d = r_j − r_i`, `ℓ = |d|`, `n̂ = d/ℓ`. 기준(곧은) 상태에서 `a_i = n̂`, `a_j = −n̂`.

  구름   ξ_i = R(a_i − n̂),  ξ_j = R(a_j + n̂)      각 입자의 마커가 현재 중심선에서 벗어난 양
         U_r = ½k_r(|ξ_i|²+|ξ_j|²) = k_r R² (2 − (a_i − a_j)·n̂)      ★ 닫힌 형태
  접선   b = a_j − a_i,  ξ_s = P b·R  (P = I − n̂n̂ᵀ)
         U_s = ½k_s R² (|b|² − (b·n̂)²)

━━ 해석적 힘·토크 (에너지에서 직접 미분) ━━
  ∂n̂/∂d = P/ℓ,   τ = −a × ∂U/∂a
  구름:  w = −k_rR²(a_i−a_j);  F_j = −Pw/ℓ = −F_i
         τ_i = +k_rR²(a_i × n̂),   τ_j = −k_rR²(a_j × n̂)
  접선:  F_j = +k_sR²(b·n̂)(Pb)/ℓ = −F_i
         τ_i = +k_sR²(a_i × Pb),  τ_j = −k_sR²(a_j × Pb)

━━ ★ 종이 위 예측 (실행 전 고정 — 원칙 9.2) ━━
  P1  기준 상태에서 U=0, F=0, τ=0.
  P2  쌍 전체를 강체회전하면 U=0 (마커와 n̂ 이 함께 돈다).
  P3  **접선 스프링만으로는 사슬 굽힘강성이 0** 이다. 방향 θ 를 이완시키면
      U_s ∝ Σ(θ_i+θ_j+2φ)² 를 정확히 0 으로 만들 수 있다 (미지수가 구속보다 많다).
      ⟹ `pair.friction` 이 구름에 면제되는 것과 **같은 구조**다.
  P4  구름 저항만이면 굽힘강성이 생긴다. 방향 이완 후
      U_r = ¼ k_r R² Σ θ_bend²  ⟹  **κ_θ,eff = ½ k_r R²**
      (내부 입자가 양쪽 결합에 공유돼 좌절되는 것이 기원. 끝 입자는 완전히 이완된다)
  P5  방향이 **얼지면** 다른 물리가 된다: U_r = k_rR² Σ φ_k² — 곡률이 아니라 결합의
      절대 회전을 벌한다(강체 고정). ⟹ 준정적에서만 조화굽힘과 같고 고주파에서 갈린다.
"""
from __future__ import annotations

import numpy as np

# ──────────────────────────────────────────────────────────────────────────
# 넘파이 참조 모델 (HOOMD 없이 — 해석해·정적 최소화용)
# ──────────────────────────────────────────────────────────────────────────


def q_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """쿼터니언 q=(w,x,y,z) 로 벡터 v 회전. (..,4)·(..,3) → (..,3)"""
    w, u = q[..., :1], q[..., 1:]
    return v + 2 * np.cross(u, np.cross(u, v) + w * v)


def q_from_axis_angle(axis: np.ndarray, ang) -> np.ndarray:
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    ang = np.atleast_1d(np.asarray(ang, float))
    h = ang / 2
    return np.concatenate([np.cos(h)[:, None], np.sin(h)[:, None] * axis[None, :]], axis=1)


class RollingContact:
    """결합 리스트 기반 구름 저항 + 접선 스프링. 마커는 생성 시 기준 배치에서 고정."""

    def __init__(self, bonds, pos0, quat0, R, k_roll, k_slide):
        self.bonds = np.asarray(bonds, int)                 # (M,2)
        self.R = float(R)
        self.k_r = float(k_roll)
        self.k_s = float(k_slide)
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        d = np.asarray(pos0, float)[j] - np.asarray(pos0, float)[i]
        n = d / np.linalg.norm(d, axis=1, keepdims=True)
        qi, qj = np.asarray(quat0, float)[i], np.asarray(quat0, float)[j]
        # 마커를 **물체 프레임**으로 저장 → q0 의 역회전
        self.m_i = q_rotate(qi * np.array([1, -1, -1, -1.0]), +n)
        self.m_j = q_rotate(qj * np.array([1, -1, -1, -1.0]), -n)

    def _kin(self, pos, quat):
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        pos, quat = np.asarray(pos, float), np.asarray(quat, float)
        d = pos[j] - pos[i]
        ell = np.linalg.norm(d, axis=1, keepdims=True)
        n = d / ell
        a_i = q_rotate(quat[i], self.m_i)
        a_j = q_rotate(quat[j], self.m_j)
        return i, j, n, ell, a_i, a_j

    def energy(self, pos, quat) -> float:
        _, _, n, _, a_i, a_j = self._kin(pos, quat)
        U = self.k_r * self.R ** 2 * (2.0 - ((a_i - a_j) * n).sum(1)).sum()
        b = a_j - a_i
        bn = (b * n).sum(1)
        U += 0.5 * self.k_s * self.R ** 2 * ((b * b).sum(1) - bn ** 2).sum()
        return float(U)

    def force_torque(self, pos, quat):
        """해석적 힘·토크. 반환 (N,3)·(N,3)."""
        i, j, n, ell, a_i, a_j = self._kin(pos, quat)
        N = len(pos)
        F = np.zeros((N, 3))
        T = np.zeros((N, 3))
        R2 = self.R ** 2

        def perp(v):                                    # P v = v − (v·n̂)n̂
            return v - ((v * n).sum(1, keepdims=True)) * n

        # ── 구름 ──
        w = -self.k_r * R2 * (a_i - a_j)
        Fj = -perp(w) / ell
        np.add.at(F, j, Fj)
        np.add.at(F, i, -Fj)
        np.add.at(T, i, self.k_r * R2 * np.cross(a_i, n))
        np.add.at(T, j, -self.k_r * R2 * np.cross(a_j, n))

        # ── 접선 ──
        b = a_j - a_i
        bn = (b * n).sum(1, keepdims=True)
        Pb = perp(b)
        Fj = self.k_s * R2 * bn * Pb / ell
        np.add.at(F, j, Fj)
        np.add.at(F, i, -Fj)
        np.add.at(T, i, self.k_s * R2 * np.cross(a_i, Pb))
        np.add.at(T, j, -self.k_s * R2 * np.cross(a_j, Pb))
        return F, T


def kappa_theta_eff(k_roll: float, R: float) -> float:
    """P4 — 방향 이완 극한에서 구름 저항 ↔ 조화 굽힘 강성의 정확한 대응."""
    return 0.5 * k_roll * R ** 2


def k_roll_from_kappa_theta(kappa_theta: float, R: float) -> float:
    return 2.0 * kappa_theta / R ** 2


# ──────────────────────────────────────────────────────────────────────────
# HOOMD 힘 (aniso=True — 토크를 쓴다)
# ──────────────────────────────────────────────────────────────────────────
def make_rolling_force(bonds, pos0, quat0, R, k_roll, k_slide, n_real):
    """`md.force.Custom(aniso=True)`. 유령입자(tag >= n_real)는 건드리지 않는다."""
    import hoomd.md as md

    model = RollingContact(bonds, pos0, quat0, R, k_roll, k_slide)

    class RollingCustom(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=True)          # ★ 토크를 내려면 필수
            self.model = model
            self.n = int(n_real)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)   # ★ tag 인덱싱 (정렬 때문)
                pos = np.array(snap.particles.position, copy=True)
                quat = np.array(snap.particles.orientation, copy=True)
                m = tags < self.n
                P = np.zeros((self.n, 3))
                Q = np.tile(np.array([1.0, 0, 0, 0]), (self.n, 1))
                P[tags[m]] = pos[m]
                Q[tags[m]] = quat[m]
                F, T = self.model.force_torque(P, Q)
                U = self.model.energy(P, Q)
                arr.force[:] = 0.0
                arr.torque[:] = 0.0
                arr.potential_energy[:] = 0.0
                arr.force[m] = F[tags[m]]
                arr.torque[m] = T[tags[m]]
                arr.potential_energy[m] = U / self.n     # 총합만 의미 있음

    return RollingCustom()
