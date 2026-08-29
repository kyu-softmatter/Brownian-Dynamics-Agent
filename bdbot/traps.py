"""조화 트랩 (`md.force.Custom`) — 고정 · 등속 이동 · 진동 구동을 한 클래스로.

**HOOMD에 조화 트랩이 없습니다** (`md.external.field` 는 Electric/Magnetic/Periodic 뿐이고
`hpmc.external.Harmonic` 은 몬테카를로 전용). skill `bd-hoomd` 의 검증된 스니펫을
일반화한 것입니다.

**세 번 나와서 올렸습니다**: 1-A `trap-2d-5um`(고정) · `trap-drag-2d-hex300`(등속 이동) ·
`chain-bend-2d-oscill`(진동). 세 경우가 전부 `앵커(t)` 만 다릅니다.

    anchors(t) = anchor0 + velocity·t + drive(t)

⚠️ 1-A(`cases/trap_2d_5um.py`)는 아직 자기 안에 같은 클래스를 갖고 있습니다. 옮기려면
   재실행 + `verify_1c_equivalence` 가 필요해서 이번 작업에서 건드리지 않았습니다.

세 가지가 **조용히 틀리는 함정**이라 여기 한 곳에 모읍니다:
  함정 1  최소 이미지 — 없으면 약한 트랩에서 +1856% 어긋난다 (에러가 나지 않는다)
  함정 7  비주기 축의 주기를 inf 로 두면 `inf*round(0/inf) = nan`
  tag     `ParticleSorter` 가 메모리 순서를 재배열한다 — 로컬 스냅샷은 tag 순이 아니다
"""
from __future__ import annotations

import numpy as np


def _md():
    import hoomd.md as md
    return md


def make_trap(k, anchors, box, *, dt_star: float = 0.0, dims: int = 2,
              velocity=None, drive=None):
    """조화 트랩 하나. `k` 는 스칼라이거나 입자별 배열(트랩 없는 입자는 0).

    인자
      k         스칼라 또는 (N,) — 무차원 강성 [kT/d²]
      anchors   (N,3) 또는 (N,2) — t=0 의 앵커 (무차원 길이)
      box       스칼라 L 또는 (Lx, Ly) — 주기 길이. 0/None 이면 그 축은 비주기
      dt_star   적분 스텝 (무차원). 앵커가 움직일 때만 필요 — t = timestep·dt_star
      velocity  (N,3) 등속 이동. `trap-drag` 가 쓴다
      drive     callable(t) -> (N,3) 오프셋. `chain-bend` 의 a·sin(ωt) 가 쓴다

    ★ 이 함수는 hoomd 를 임포트할 때만 클래스를 만듭니다 — `bdbot.cli` 같은 앞단이
      무거운 의존성을 끌어오지 않게 하려는 것입니다 (bdbot/__init__ 의 규약).
    """
    md = _md()

    class HarmonicTrap(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=False)
            a = np.asarray(anchors, dtype=float)
            if a.shape[1] == 2:
                a = np.c_[a, np.zeros(len(a))]
            self.anchor0 = a
            n = len(a)
            self.k = (np.full(n, float(k)) if np.isscalar(k)
                      else np.asarray(k, dtype=float).reshape(n))
            b = (box, box) if np.isscalar(box) else tuple(box)
            # ★ 함정 7: 비주기 축은 0 으로 둬서 마스크에서 빠지게 한다 (inf 금지)
            self.period = np.array([float(b[0] or 0.0), float(b[1] or 0.0), 0.0])
            self.dt_star = float(dt_star)
            self.velocity = None if velocity is None else np.asarray(velocity, dtype=float)
            self.drive = drive
            self.dims = dims

        def centers(self, t: float) -> np.ndarray:
            c = self.anchor0
            if self.velocity is not None:
                c = c + self.velocity * t
            if self.drive is not None:
                c = c + np.asarray(self.drive(t), dtype=float)
            return c

        def _delta(self, pos, tags, t):
            d = pos - self.centers(t)[tags]
            m = self.period > 0                       # 주기 축만 래핑 (함정 1·7)
            d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
            return d

        def set_forces(self, timestep):
            t = timestep * self.dt_star
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)     # ★ tag 인덱싱 필수
                pos = np.array(snap.particles.position, copy=True)
                d = self._delta(pos, tags, t)
                kk = self.k[tags][:, None]
                arr.force[:] = -kk * d
                arr.potential_energy[:] = 0.5 * self.k[tags] * (d ** 2).sum(axis=1)

        def displacement(self, state, timestep) -> np.ndarray:
            """트랩 중심으로부터의 변위 (N,3). 끌림힘·위상 추출의 원자료."""
            snap = state.get_snapshot()
            pos = np.array(snap.particles.position, dtype=float)
            tags = np.arange(len(pos))       # get_snapshot 은 tag 순으로 정렬되어 있다
            return self._delta(pos, tags, timestep * self.dt_star)

    return HarmonicTrap()


__all__ = ["make_trap"]
