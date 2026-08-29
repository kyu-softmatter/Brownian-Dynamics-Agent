"""Rolling resistance + tangential spring contact model -- an
`md.force.Custom(aniso=True)` implementation plus a numpy reference model.

It was measured that `md.pair.friction` cannot produce a bending stiffness (it is
dissipative and exempt from rolling -- `scratch/verify_pair_friction.py`). What a JKR
adhesive contact actually does is **resist the contact point rolling over the
surface**, so that is implemented directly (in the Dominik-Tielens 1995 / Wada 2007
family).

-- Kinematics --
For each bond (i,j) a **body-fixed contact-point marker** `m_i` (unit vector from the
centre of particle i to the contact point) is stored at the moment the bond is
created. Its current lab-frame value is `a_i = Rot(q_i).m_i`.
`d = r_j - r_i`, `l = |d|`, `n = d/l`. In the reference (straight) state
`a_i = n` and `a_j = -n`.

  rolling    xi_i = R(a_i - n),  xi_j = R(a_j + n)
             how far each particle's marker has moved off the current centre line
             U_r = 0.5*k_r(|xi_i|^2+|xi_j|^2) = k_r R^2 (2 - (a_i - a_j).n)
             ★ closed form
  tangential b = a_j - a_i,  xi_s = P b * R  (P = I - n n^T)
         U_s = ½k_s R² (|b|² − (b·n̂)²)

-- Analytic force and torque (differentiated directly from the energy) --
  ∂n̂/∂d = P/ℓ,   τ = −a × ∂U/∂a
  rolling:    w = -k_r R^2 (a_i - a_j);  F_j = -P w/l = -F_i
         τ_i = +k_rR²(a_i × n̂),   τ_j = −k_rR²(a_j × n̂)
  tangential: F_j = +k_s R^2 (b.n)(P b)/l = -F_i
         τ_i = +k_sR²(a_i × Pb),  τ_j = −k_sR²(a_j × Pb)

-- ★ Predictions written down on paper (sealed before running -- principle 9.2) --
  P1  In the reference state U=0, F=0, torque=0.
  P2  Rigidly rotating the whole pair gives U=0 (the markers and n turn together).
  P3  **The tangential spring alone gives ZERO chain bending stiffness.** Relaxing
      the orientations theta lets U_s ~ sum(theta_i + theta_j + 2*phi)^2 be driven
      to exactly 0 (there are more unknowns than constraints).
      => the **same structure** as pair.friction being exempt from rolling.
  P4  Rolling resistance alone does give a bending stiffness. After relaxing the
      orientations
      U_r = ¼ k_r R² Σ θ_bend²  ⟹  **κ_θ,eff = ½ k_r R²**
      (the origin is that an interior particle is shared by two bonds and so is
       frustrated; an end particle relaxes completely)
  P5  If the orientations are **frozen** it becomes different physics:
      U_r = k_r R^2 sum(phi_k^2) -- it penalises the absolute rotation of each bond
      rather than the curvature (a rigid clamp). => it matches a harmonic bend only
      quasi-statically, and diverges from it at high frequency.
"""
from __future__ import annotations

import numpy as np

# ──────────────────────────────────────────────────────────────────────────
# numpy reference model (no HOOMD -- for the analytic solution and static
# minimization)
# ──────────────────────────────────────────────────────────────────────────


def q_rotate(q: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Rotate vector v by quaternion q=(w,x,y,z). (..,4) and (..,3) -> (..,3)"""
    w, u = q[..., :1], q[..., 1:]
    return v + 2 * np.cross(u, np.cross(u, v) + w * v)


def q_from_axis_angle(axis: np.ndarray, ang) -> np.ndarray:
    axis = np.asarray(axis, float)
    axis = axis / np.linalg.norm(axis)
    ang = np.atleast_1d(np.asarray(ang, float))
    h = ang / 2
    return np.concatenate([np.cos(h)[:, None], np.sin(h)[:, None] * axis[None, :]], axis=1)


class RollingContact:
    """Bond-list based rolling resistance + tangential spring.

    Markers are fixed from the reference configuration at construction time.
    """

    def __init__(self, bonds, pos0, quat0, R, k_roll, k_slide):
        self.bonds = np.asarray(bonds, int)                 # (M,2)
        self.R = float(R)
        self.k_r = float(k_roll)
        self.k_s = float(k_slide)
        i, j = self.bonds[:, 0], self.bonds[:, 1]
        d = np.asarray(pos0, float)[j] - np.asarray(pos0, float)[i]
        n = d / np.linalg.norm(d, axis=1, keepdims=True)
        qi, qj = np.asarray(quat0, float)[i], np.asarray(quat0, float)[j]
        # Store the markers in the **body frame** -> inverse rotation by q0
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
        """Analytic force and torque. Returns (N,3) and (N,3)."""
        i, j, n, ell, a_i, a_j = self._kin(pos, quat)
        N = len(pos)
        F = np.zeros((N, 3))
        T = np.zeros((N, 3))
        R2 = self.R ** 2

        def perp(v):                                    # P v = v − (v·n̂)n̂
            return v - ((v * n).sum(1, keepdims=True)) * n

        # ── rolling ──
        w = -self.k_r * R2 * (a_i - a_j)
        Fj = -perp(w) / ell
        np.add.at(F, j, Fj)
        np.add.at(F, i, -Fj)
        np.add.at(T, i, self.k_r * R2 * np.cross(a_i, n))
        np.add.at(T, j, -self.k_r * R2 * np.cross(a_j, n))

        # ── tangential ──
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
    """P4 -- the exact correspondence between rolling resistance and a harmonic
    bending stiffness, in the orientation-relaxed limit."""
    return 0.5 * k_roll * R ** 2


def k_roll_from_kappa_theta(kappa_theta: float, R: float) -> float:
    return 2.0 * kappa_theta / R ** 2


# ──────────────────────────────────────────────────────────────────────────
# HOOMD force (aniso=True -- it uses torques)
# ──────────────────────────────────────────────────────────────────────────
def make_rolling_force(bonds, pos0, quat0, R, k_roll, k_slide, n_real):
    """`md.force.Custom(aniso=True)`. Ghost particles (tag >= n_real) are untouched."""
    import hoomd.md as md

    model = RollingContact(bonds, pos0, quat0, R, k_roll, k_slide)

    class RollingCustom(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=True)          # ★ required in order to emit torques
            self.model = model
            self.n = int(n_real)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)   # ★ tag indexing (because of reordering)
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
                arr.potential_energy[m] = U / self.n     # only the total is meaningful

    return RollingCustom()
