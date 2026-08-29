"""HOOMD BD 실행 골격 — 두 케이스가 똑같이 쓴 부분만.

공통: 2D 프레임 만들기 · Simulation · Brownian + Integrator · GSD Tier A 라이터 · 시드 처리
**케이스마다 다른 것**: 힘(forces), 초기 배치, 표본 수집 루프 → 케이스에 남깁니다.

여기 있는 규약은 전부 skill `bd-hoomd`에서 실측 검증된 것입니다:
  함정 3  `integrate_rotational_dof = False`
  함정 5  BD는 과감쇠 — `thermalize_particle_momenta()` 불필요, 속도에 의미 없음
  함정 9  2D는 `Lz=0` + `dimensions=2` 둘 다
  함정 12 시드가 16비트로 잘림
"""
from __future__ import annotations

import numpy as np

HOOMD_SEED_MAX = 65535


def resolve_seed(seed: int) -> tuple[int, int]:
    """(numpy 시드, HOOMD 시드). HOOMD는 16비트로 잘라 쓰므로 미리 잘라 넘깁니다.

    ★ 함정 12: `seed=20260803` 을 주면 HOOMD가 경고와 함께 `10179`로 자릅니다.
      65536만큼 다른 두 시드는 **같은 궤적**이 됩니다. numpy(초기배치)는 전체 시드를
      쓰므로 둘을 분리해 넘깁니다.
    """
    return int(seed), int(seed) & HOOMD_SEED_MAX


def frame_2d(positions, L, types=("A",), typeid=None, orientation=False,
             bonds=None, angles=None):
    """2D 주기 프레임. `positions`는 (N,2) 또는 (N,3).

    `L` 은 스칼라(정사각) 또는 `(Lx, Ly)` 입니다 — ★ `trap-drag` 의 정합 육방 격자가
    직사각 박스를 요구하면서 필요해졌습니다 (정합이 종횡비를 정합니다).

    `bonds`/`angles` 는 (M,2)/(M,3) 인덱스 배열 — 사슬 케이스가 씁니다.
    """
    import gsd.hoomd

    p = np.asarray(positions, dtype=float)
    if p.shape[1] == 2:
        p = np.c_[p, np.zeros(len(p))]
    n = len(p)
    Lx, Ly = (float(L), float(L)) if np.isscalar(L) else (float(L[0]), float(L[1]))
    fr = gsd.hoomd.Frame()
    fr.particles.N = n
    fr.particles.position = p
    fr.particles.typeid = [0] * n if typeid is None else list(typeid)
    fr.particles.types = list(types)
    if orientation:
        fr.particles.orientation = [(1.0, 0.0, 0.0, 0.0)] * n
    if bonds is not None:
        b = np.asarray(bonds, dtype=int)
        fr.bonds.N = len(b)
        fr.bonds.types = ["backbone"]
        fr.bonds.typeid = [0] * len(b)
        fr.bonds.group = b
    if angles is not None:
        a = np.asarray(angles, dtype=int)
        fr.angles.N = len(a)
        fr.angles.types = ["bend"]
        fr.angles.typeid = [0] * len(a)
        fr.angles.group = a
    fr.configuration.box = [Lx, Ly, 0, 0, 0, 0]   # Lz=0 → 2D (함정 9)
    fr.configuration.dimensions = 2
    return fr


def make_sim(frame, seed: int, notice_level: int = 0):
    import hoomd

    _, hseed = resolve_seed(seed)
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=notice_level), seed=hseed)
    sim.create_state_from_snapshot(frame)
    return sim


def attach_brownian(sim, dt_star: float, forces, kT: float = 1.0, gamma: float = 1.0):
    """무차원 BD 적분기를 붙이고 (integrator, method)를 반환.

    `kT=1, gamma=1` 은 thermal 규약의 결과입니다 (σ=d, E=kT, τ=τ_B → γ=1, D_t=1).
    """
    import hoomd
    import hoomd.md as md

    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=kT, default_gamma=gamma)
    integ = md.Integrator(dt=dt_star, methods=[bd], forces=list(forces))
    integ.integrate_rotational_dof = False          # 함정 3
    sim.operations.integrator = integ
    return integ, bd


def add_trajectory_writer(sim, path, period: int):
    """Tier A 궤적 (마스터플랜 §9). `path=None` 이면 아무것도 안 붙입니다."""
    import hoomd

    if path is None:
        return None
    p = max(1, int(period))
    wr = hoomd.write.GSD(filename=str(path), trigger=hoomd.trigger.Periodic(p),
                         mode="xb", dynamic=["property"])
    sim.operations.writers.append(wr)
    return wr


def flush_writers(sim) -> None:
    for w in sim.operations.writers:
        if hasattr(w, "flush"):
            w.flush()


def wca(nlist, epsilon: float = 1.0, sigma: float = 1.0, types=("A", "A")):
    """WCA = 컷오프 `2^(1/6)σ` + shift 인 LJ (함정 4: WCA 전용 클래스는 없음)."""
    import hoomd.md as md

    lj = md.pair.LJ(nlist=nlist, default_r_cut=2 ** (1 / 6) * sigma, mode="shift")
    lj.params[types] = dict(epsilon=epsilon, sigma=sigma)
    return lj


def minimum_image(delta, L, dims: int = 2):
    """주기 축만 최소 이미지 적용 (함정 1·7). `delta` (N,3), 2D면 z는 건드리지 않음.

    `L` 은 스칼라 또는 `(Lx, Ly)` (직사각 박스 — 정합 육방 격자가 요구).
    ★ 비주기 축의 주기를 `inf`로 두면 `inf*round(0/inf) = nan` 이 됩니다 (함정 7).
    """
    d = np.asarray(delta, dtype=float).copy()
    Lx, Ly = (L, L) if np.isscalar(L) else (L[0], L[1])
    period = np.array([float(Lx), float(Ly), float(Lx) if dims == 3 else 0.0])
    m = period > 0
    d[:, m] -= period[m] * np.round(d[:, m] / period[m])
    return d


def progress(i, total, t_elapsed, extra: str = "") -> str:
    pct = 100 * i / total if total else 0.0
    return f"    {i:>6}/{total}  ({pct:4.0f}%)  {t_elapsed:6.1f}s   {extra}"


__all__ = ["resolve_seed", "frame_2d", "make_sim", "attach_brownian", "add_trajectory_writer",
           "flush_writers", "wca", "minimum_image", "progress", "HOOMD_SEED_MAX"]
