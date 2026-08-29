"""초기 배치 생성.

조화 트랩 계는 특수하다 — 쌍 상호작용이 없으므로 **겹침이라는 개념이 없다.**
전 입자를 트랩 중심에 놓고 시작해도 되고, 그게 가장 단순하다.
평형화가 그 초기조건을 지운다 (`t_equil >= 10 tau_trap`).
"""
from __future__ import annotations

import numpy as np


def trap_snapshot(hoomd, *, n: int, dim: int, box_over_l_trap: float,
                  start: str = "center", seed: int = 0):
    """비상호작용 입자 `n` 개를 조화 트랩 계에 배치한 Snapshot.

    Args:
        dim: 2 또는 3.  2D 는 `Lz = 0` 으로 지정한다 (HOOMD 7 은
            `configuration.dimensions` 에 setter 가 없다).
        start: "center" — 전부 원점. "equilibrium" — Boltzmann 분포에서 표집.
            축약 단위에서 평형 표준편차는 1 (`<x*^2> = 1`).

    ⚠ 박스는 `l_trap` 의 `box_over_l_trap` 배다. 트랩 카드 §7 게이트가
      `L >> l_trap` 을 요구한다 — `HarmonicTrap` 이 래핑된 좌표를 쓰기 때문이다.
    """
    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}")
    if box_over_l_trap < 20:
        raise ValueError(
            f"box_over_l_trap = {box_over_l_trap} 은 너무 작다. 트랩 카드 §7 게이트: "
            f"박스가 l_trap 보다 훨씬 커야 한다 (래핑된 좌표를 쓰므로)")

    L = float(box_over_l_trap)
    snap = hoomd.Snapshot()
    snap.configuration.box = [L, L, (0.0 if dim == 2 else L), 0, 0, 0]
    snap.particles.N = int(n)
    snap.particles.types = ["A"]

    pos = np.zeros((n, 3), dtype=np.float64)
    if start == "equilibrium":
        rng = np.random.default_rng(seed)
        pos[:, :dim] = rng.normal(0.0, 1.0, size=(n, dim))   # <x*^2> = 1
    elif start != "center":
        raise ValueError(f"unknown start {start!r}")

    snap.particles.position[:] = pos
    snap.particles.typeid[:] = 0
    snap.particles.mass[:] = 1.0
    snap.particles.moment_inertia[:] = 0.0        # 회전 자유도 끔 (findings §5)
    return snap


# =============================================================================
# 콜로이드 사슬 — 굽힘 강성 계
# =============================================================================
def chain_snapshot(hoomd, *, n: int, dim: int = 2, bond_length: float = 1.0,
                   box_over_length: float = 4.0, end_type: bool = True,
                   center_type: bool = True):
    """직선 사슬 `n` 개. 결합 `n-1` 개 + 각 `n-2` 개.

    축약 단위: 길이 = sigma = 2a = 결합 길이 = 1.

    Args:
        end_type: True 면 양 끝 입자를 타입 `"E"` 로 둔다. 끝을 고정하거나
            트랩을 걸 때 `hoomd.filter.Type` 으로 잡기 위한 것 —
            태그 기반 필터보다 스냅샷을 봐서 검사하기 쉽다.
        center_type: True 면 중앙 입자(홀수 `n` 일 때만)를 타입 `"C"` 로 둔다.
            3점 굽힘에서 하중을 거는 자리다. `md.force.Constant` 는 타입별로
            힘을 주므로 중앙을 타입으로 분리해야 한다.

    ⚠ **거리 구속(`constrain.Distance`)을 쓰지 말 것.** `Brownian` 과 호환되지 않고
      조용히 발산한다. 근거: findings/dead-end-distance-constraint-with-brownian.md
      결합은 `md.bond.Harmonic` 으로 하고 `k_bond* >> kappa(N)*` 를 만족시킨다.

    ⚠ 사슬은 x 축에 놓이고 중심이 원점이다. 박스는 사슬 길이의 `box_over_length` 배 —
      주기 이미지와 상호작용하지 않아야 한다 (쌍포텐셜이 없으므로 결합·각만 신경쓰면 되고,
      끝-끝이 박스를 넘어 최소이미지로 잘못 이어지는 것만 막으면 된다).
    """
    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}")
    if n < 3:
        raise ValueError(f"n = {n} — 각(angle)이 성립하려면 최소 3개 필요하다")

    import numpy as np

    contour = (n - 1) * float(bond_length)
    L = max(contour * float(box_over_length), 4.0 * bond_length)

    use_center = bool(center_type) and n % 2 == 1
    types = ["A"]
    if end_type:
        types.append("E")
    if use_center:
        types.append("C")

    snap = hoomd.Snapshot()
    snap.configuration.box = [L, L, (0.0 if dim == 2 else L), 0, 0, 0]
    snap.particles.N = int(n)
    snap.particles.types = types

    pos = np.zeros((n, 3), dtype=np.float64)
    pos[:, 0] = (np.arange(n) - (n - 1) / 2.0) * float(bond_length)
    snap.particles.position[:] = pos

    tid = np.zeros(n, dtype=int)
    if end_type:
        tid[0] = tid[-1] = types.index("E")
    if use_center:
        tid[(n - 1) // 2] = types.index("C")
    snap.particles.typeid[:] = tid
    snap.particles.mass[:] = 1.0
    snap.particles.moment_inertia[:] = 0.0

    snap.bonds.N = n - 1
    snap.bonds.types = ["b"]
    snap.bonds.group[:] = [[i, i + 1] for i in range(n - 1)]
    snap.bonds.typeid[:] = 0

    snap.angles.N = n - 2
    snap.angles.types = ["a"]
    snap.angles.group[:] = [[i, i + 1, i + 2] for i in range(n - 2)]
    snap.angles.typeid[:] = 0

    return snap


# =============================================================================
# 2D 쌍 상호작용 계 — 길이 단위는 `d = n^{-1/2}` (Zahn 규약)
# =============================================================================
#  ★ 왜 `d = n^{-1/2}` 인가 (최근접거리가 아니라)
#  Zahn 1999 의 `Γ = (μ₀/4π)M²π^{3/2}/(kT d³)` 에서 `d³ = n^{-3/2}` 이므로
#  `d ≡ n^{-1/2}` 다. 육방격자의 최근접거리는 `a = (2/√3)^{1/2} d = 1.0746 d` 로
#  7 % 다르다. 이 7 % 를 혼동하면 `Γ` 가 23 % 틀린다 (`Γ ∝ d^{-3}`).
#  근거: knowledge/source/papers/1999-zahn-two-stage-melting-2d.md §2
HEX_NN_OVER_D = (2.0 / np.sqrt(3.0)) ** 0.5        # 1.074570...


def square_box_for(n_particles: int, density_star: float = 1.0) -> float:
    """`n* = N/L²` 를 만족하는 정사각 상자 변 길이 (단위 `d`).

    `density_star = 1` 이면 정의상 `L = √N` 이다 — 길이 단위가 `n^{-1/2}` 이므로.
    """
    return float(np.sqrt(n_particles / density_star))


# --- 기준 원판 커버리지 ↔ 상자 크기 ------------------------------------------
#  ★ 이 계에는 경질 코어가 없다 (`U = A/r^p` 뿐). 따라서 `σ` 는 **동역학에 들어가지
#    않는다** — 카드 §"σ 가 없다". 그런데도 `σ` 가 필요한 이유는 둘이다:
#      ① 시간 척도.  `τ_d = d²/D₀`,  `D₀ = kT/(3πησ)` — `σ` 없이는 초를 못 붙인다.
#      ② 물리적 타당성.  "5 µm 원판 100개" 의 면적 점유율이 얼마인지가
#         점입자 이상화가 정당한가를 결정한다.
#  ⇒ `σ/d` 만이 자유도이고, **무차원 물리는 `A` 하나로 완전히 결정된다**
#    (`n* = 1` 이 정의상 성립하므로). 커버리지를 바꾸면 축의 라벨만 바뀐다.
def coverage_from_sigma_over_d(sigma_over_d: float) -> float:
    """기준 원판(직경 `σ`)의 면적 점유율.  `φ = (π/4)(σ/d)²`.

    `n* = 1` (길이 단위가 `d = n^{-1/2}`) 이므로 입자당 면적이 정확히 `d²` 다.
    따라서 `φ = (πσ²/4)/d²` 이고 `N` 과 `L` 에 **따로 의존하지 않는다.**
    """
    return float(np.pi / 4.0 * sigma_over_d**2)


def sigma_over_d_for_coverage(coverage: float) -> float:
    """주어진 커버리지를 만드는 `σ/d`. `coverage_from_sigma_over_d` 의 역함수."""
    if not 0.0 < coverage <= 1.0:
        raise ValueError(f"커버리지 {coverage} 는 (0, 1] 이어야 한다")
    return float(np.sqrt(4.0 * coverage / np.pi))


def box_si_for_coverage(*, n_particles: int, sigma_si: float, coverage_max: float,
                        d_over_sigma_round: float | None = None) -> dict:
    """"`σ` 원판 `N` 개를 커버리지 `≤ coverage_max` 로 담는 정사각 상자".

    Args:
        d_over_sigma_round: 지정하면 `d/σ` 를 **이 값으로 고정**하고 그 결과
            커버리지를 돌려준다 (깔끔한 정수비를 쓰고 싶을 때). 커버리지가
            `coverage_max` 를 넘으면 예외.

    Returns: `d_si` · `L_si` · `coverage` · `sigma_over_d` · `L_star`(= √N).
    """
    sigma_over_d_max = sigma_over_d_for_coverage(coverage_max)
    if d_over_sigma_round is None:
        d_over_sigma = 1.0 / sigma_over_d_max
    else:
        d_over_sigma = float(d_over_sigma_round)
    sigma_over_d = 1.0 / d_over_sigma
    cov = coverage_from_sigma_over_d(sigma_over_d)
    if cov > coverage_max:
        raise ValueError(
            f"d/σ = {d_over_sigma:g} 는 커버리지 {cov:.4%} > 상한 {coverage_max:.4%} "
            f"를 준다. d/σ ≥ {1 / sigma_over_d_max:.4f} 이어야 한다")
    d_si = d_over_sigma * sigma_si
    L_star = square_box_for(n_particles, 1.0)
    return {"d_si": d_si, "L_si": L_star * d_si, "coverage": cov,
            "sigma_over_d": sigma_over_d, "d_over_sigma": d_over_sigma,
            "L_star": L_star, "coverage_max": coverage_max,
            "sigma_over_d_max": sigma_over_d_max}


def random_2d_snapshot(hoomd, *, n: int, box: float | None = None,
                       box_x: float | None = None, box_y: float | None = None,
                       min_sep: float = 0.5, seed: int = 0,
                       max_tries: int = 200):
    """2D 상자에 **무작위 비중첩** 배치. 정사각(`box`) 또는 직사각(`box_x`,`box_y`).

    ★ 직사각을 받는 이유: 초기조건(hex vs random)을 비교할 때 **상자 모양이 같아야**
      한다. 다르면 `r_cut`(= min(L)/2 에서 파생)과 허용 k-벡터가 함께 바뀌어
      비교가 교란된다 (2026-07-29 실측: S(k) 6겹 변조가 50배 갈렸다).

    `min_sep` 보다 가까운 쌍이 없도록 기각표집한다. 겹침을 남기면 `A/r^p` 가
    발산해서 1 스텝에 폭발한다 (`master_plan` §S5 실패모드).

    실패하면 **예외를 던진다** — 겹친 배치를 조용히 돌려주면 나중에 원인을 못 찾는다.
    """
    if box is not None:
        box_x = box_y = box
    if box_x is None or box_y is None:
        raise ValueError("box 또는 (box_x, box_y) 를 줘야 한다")
    L = np.array([box_x, box_y])
    rng = np.random.default_rng(seed)
    pos = np.empty((n, 2))
    for i in range(n):
        for _ in range(max_tries):
            p = rng.uniform(-L / 2, L / 2)
            if i == 0:
                break
            d = pos[:i] - p
            d -= L * np.round(d / L)                     # 최소이미지 (성분별)
            if np.min(np.linalg.norm(d, axis=1)) >= min_sep:
                break
        else:
            raise RuntimeError(
                f"입자 {i}/{n} 를 min_sep={min_sep} 로 넣지 못했다 "
                f"(box={box_x:.3f}x{box_y:.3f}, 시도 {max_tries}회). 밀도가 너무 "
                f"높거나 min_sep 이 과하다 — soft pushoff 가 필요하다")
        pos[i] = p
    return _make_2d_snapshot(hoomd, pos, box_x, box_y)


def hex_2d_snapshot(hoomd, *, n_x: int, n_y: int, density_star: float = 1.0):
    """2D **육방 격자**. `N = 2 n_x n_y` (단위격자에 2개).

    상자는 격자와 정합하게 잡는다 (`Lx = n_x a`, `Ly = n_y a√3`) — 정합하지
    않으면 주기경계가 격자를 끊어 인공 결함이 생긴다.

    Returns:
        `(snapshot, info)` — `info` 에 `a`(최근접거리), 상자, 종횡비.
    """
    a = HEX_NN_OVER_D / np.sqrt(density_star)
    Lx, Ly = n_x * a, n_y * a * np.sqrt(3.0)
    basis = np.array([[0.0, 0.0], [0.5 * a, 0.5 * a * np.sqrt(3.0)]])
    cells = np.array([[i * a, j * a * np.sqrt(3.0)]
                      for i in range(n_x) for j in range(n_y)])
    pos = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 2)
    pos[:, 0] = (pos[:, 0] + Lx / 2) % Lx - Lx / 2       # 상자 중심으로
    pos[:, 1] = (pos[:, 1] + Ly / 2) % Ly - Ly / 2
    info = {"a_nn": float(a), "Lx": float(Lx), "Ly": float(Ly),
            "n_particles": int(pos.shape[0]),
            "aspect": float(Lx / Ly),
            "density_star": float(pos.shape[0] / (Lx * Ly))}
    return _make_2d_snapshot(hoomd, pos, Lx, Ly), info


def hex_tiling_for(n_target: int) -> tuple[int, int]:
    """`N = 2 n_x n_y ≈ n_target` 이면서 상자가 가장 정사각에 가까운 `(n_x, n_y)`.

    종횡비는 `Lx/Ly = n_x/(n_y√3)` 이므로 `n_x ≈ √3 n_y` 가 정사각이다.
    """
    best, best_cost = None, float("inf")
    for n_y in range(1, n_target):
        n_x = max(1, round(n_target / (2 * n_y)))
        n = 2 * n_x * n_y
        aspect = n_x / (n_y * np.sqrt(3.0))
        cost = abs(n - n_target) / n_target + abs(np.log(aspect))
        if cost < best_cost:
            best, best_cost = (n_x, n_y), cost
    return best


def _make_2d_snapshot(hoomd, pos_2d: np.ndarray, Lx: float, Ly: float):
    snap = hoomd.Snapshot()
    snap.configuration.box = [Lx, Ly, 0, 0, 0, 0]        # Lz = 0 → 2D
    snap.particles.N = int(pos_2d.shape[0])
    snap.particles.types = ["A"]
    snap.particles.position[:, :2] = pos_2d
    snap.particles.position[:, 2] = 0.0
    snap.particles.typeid[:] = 0
    snap.particles.mass[:] = 1.0
    snap.particles.moment_inertia[:] = 0.0
    return snap
