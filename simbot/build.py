"""Generating the initial placement.

The harmonic-trap system is special -- with no pair interaction **there is no such
thing as an overlap.** Every particle may start at the trap centre, and that is the
simplest choice. Equilibration erases that initial condition
(`t_equil >= 10 tau_trap`).
"""
from __future__ import annotations

import numpy as np


def trap_snapshot(hoomd, *, n: int, dim: int, box_over_l_trap: float,
                  start: str = "center", seed: int = 0):
    """A Snapshot placing `n` non-interacting particles in a harmonic-trap system.

    Args:
        dim: 2 or 3.  2D is specified as `Lz = 0` (HOOMD 7 has no setter on
            `configuration.dimensions`).
        start: "center" — all at the origin. "equilibrium" — sampled from the
            Boltzmann distribution. In reduced units the equilibrium standard
            deviation is 1 (`<x*^2> = 1`).

    ⚠ The box is `box_over_l_trap` times `l_trap`. The trap card's §7 gate requires
      `L >> l_trap`, because `HarmonicTrap` uses wrapped coordinates.
    """
    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}")
    if box_over_l_trap < 20:
        raise ValueError(
            f"box_over_l_trap = {box_over_l_trap} is too small. Trap card §7 gate: "
            f"the box has to be much larger than l_trap (wrapped coordinates)")

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
    snap.particles.moment_inertia[:] = 0.0        # rotational DOF off (findings §5)
    return snap


# =============================================================================
# Colloidal chain — the bending-stiffness system
# =============================================================================
def chain_snapshot(hoomd, *, n: int, dim: int = 2, bond_length: float = 1.0,
                   box_over_length: float = 4.0, end_type: bool = True,
                   center_type: bool = True):
    """A straight chain of `n` particles. `n-1` bonds + `n-2` angles.

    Reduced units: length = sigma = 2a = bond length = 1.

    Args:
        end_type: with True, the two end particles get type `"E"`. This is so that
            `hoomd.filter.Type` can catch them when pinning the ends or attaching a
            trap -- easier to inspect from the snapshot than a tag-based filter.
        center_type: with True, the centre particle (only for odd `n`) gets type
            `"C"`. That is where the load is applied in three-point bending.
            `md.force.Constant` applies force per type, so the centre has to be
            separated into its own type.

    ⚠ **Do not use a distance constraint (`constrain.Distance`).** It is
      incompatible with `Brownian` and diverges quietly. Basis:
      findings/dead-end-distance-constraint-with-brownian.md
      Use `md.bond.Harmonic` for the bonds and satisfy `k_bond* >> kappa(N)*`.

    ⚠ The chain lies on the x axis, centred on the origin. The box is
      `box_over_length` times the chain length -- it must not interact with its
      periodic image (there is no pair potential, so only bonds and angles matter,
      and all that has to be prevented is an end-to-end pair being wrongly joined
      across the box by minimum image).
    """
    if dim not in (2, 3):
        raise ValueError(f"dim must be 2 or 3, got {dim}")
    if n < 3:
        raise ValueError(f"n = {n} — an angle needs at least 3 particles")

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
# 2D pair-interacting systems — the length unit is `d = n^{-1/2}` (Zahn convention)
# =============================================================================
#  ★ Why `d = n^{-1/2}` and not the nearest-neighbour distance
#  In Zahn 1999's `Γ = (μ₀/4π)M²π^{3/2}/(kT d³)`, `d³ = n^{-3/2}`, so
#  `d ≡ n^{-1/2}`. A hexagonal lattice's nearest-neighbour distance is
#  `a = (2/√3)^{1/2} d = 1.0746 d`, 7 % different. Confuse that 7 % and `Γ` is
#  wrong by 23 % (`Γ ∝ d^{-3}`).
#  Basis: knowledge/source/papers/1999-zahn-two-stage-melting-2d.md §2
HEX_NN_OVER_D = (2.0 / np.sqrt(3.0)) ** 0.5        # 1.074570...


def square_box_for(n_particles: int, density_star: float = 1.0) -> float:
    """The square box edge (in `d`) satisfying `n* = N/L²`.

    At `density_star = 1`, `L = √N` by definition -- because the length unit is
    `n^{-1/2}`.
    """
    return float(np.sqrt(n_particles / density_star))


# --- reference-disc coverage ↔ box size --------------------------------------
#  ★ This system has no hard core (only `U = A/r^p`). So `σ` **does not enter the
#    dynamics** -- card §"there is no σ". There are nonetheless two reasons `σ` is
#    needed:
#      ① the time scale.  `τ_d = d²/D₀`,  `D₀ = kT/(3πησ)` — without `σ` there is
#         no way to attach seconds.
#      ② physical validity.  What area fraction "100 discs of 5 µm" occupy decides
#         whether the point-particle idealization is justified.
#  ⇒ `σ/d` is the only degree of freedom, and **the dimensionless physics is fully
#    determined by `A` alone** (because `n* = 1` holds by definition). Changing the
#    coverage changes only the axis labels.
def coverage_from_sigma_over_d(sigma_over_d: float) -> float:
    """Area fraction of the reference disc (diameter `σ`).  `φ = (π/4)(σ/d)²`.

    Since `n* = 1` (the length unit is `d = n^{-1/2}`), the area per particle is
    exactly `d²`. So `φ = (πσ²/4)/d²` and it **does not depend separately on `N`
    or `L`.**
    """
    return float(np.pi / 4.0 * sigma_over_d**2)


def sigma_over_d_for_coverage(coverage: float) -> float:
    """The `σ/d` that produces a given coverage. Inverse of
    `coverage_from_sigma_over_d`."""
    if not 0.0 < coverage <= 1.0:
        raise ValueError(f"coverage {coverage} has to be in (0, 1]")
    return float(np.sqrt(4.0 * coverage / np.pi))


def box_si_for_coverage(*, n_particles: int, sigma_si: float, coverage_max: float,
                        d_over_sigma_round: float | None = None) -> dict:
    """"A square box holding `N` discs of diameter `σ` at coverage
    `≤ coverage_max`".

    Args:
        d_over_sigma_round: given, `d/σ` is **fixed to this value** and the
            resulting coverage is returned (for when a clean integer ratio is
            wanted). Raises if the coverage exceeds `coverage_max`.

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
            f"d/σ = {d_over_sigma:g} gives coverage {cov:.4%} > the cap "
            f"{coverage_max:.4%}. d/σ ≥ {1 / sigma_over_d_max:.4f} is required")
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
    """A **random non-overlapping** placement in a 2D box. Square (`box`) or
    rectangular (`box_x`, `box_y`).

    ★ Why it accepts a rectangle: when comparing initial conditions (hex vs random)
      **the box shape has to be the same.** Different shapes change `r_cut` (derived
      from min(L)/2) and the allowed k-vectors together, which contaminates the
      comparison (measured 2026-07-29: the S(k) 6-fold modulation differed by 50x).

    Rejection sampling ensures no pair is closer than `min_sep`. Leave an overlap
    and `A/r^p` diverges and blows up in one step (`master_plan` §S5 failure modes).

    On failure it **raises** -- returning an overlapped placement quietly makes the
    cause unfindable later.
    """
    if box is not None:
        box_x = box_y = box
    if box_x is None or box_y is None:
        raise ValueError("either box or (box_x, box_y) has to be given")
    L = np.array([box_x, box_y])
    rng = np.random.default_rng(seed)
    pos = np.empty((n, 2))
    for i in range(n):
        for _ in range(max_tries):
            p = rng.uniform(-L / 2, L / 2)
            if i == 0:
                break
            d = pos[:i] - p
            d -= L * np.round(d / L)                     # minimum image, per component
            if np.min(np.linalg.norm(d, axis=1)) >= min_sep:
                break
        else:
            raise RuntimeError(
                f"could not place particle {i}/{n} at min_sep={min_sep} "
                f"(box={box_x:.3f}x{box_y:.3f}, {max_tries} tries). Either the "
                f"density is too high or min_sep is excessive — a soft pushoff is "
                f"needed")
        pos[i] = p
    return _make_2d_snapshot(hoomd, pos, box_x, box_y)


def hex_2d_snapshot(hoomd, *, n_x: int, n_y: int, density_star: float = 1.0):
    """A 2D **hexagonal lattice**. `N = 2 n_x n_y` (2 per unit cell).

    The box is made commensurate with the lattice (`Lx = n_x a`, `Ly = n_y a√3`) --
    otherwise the periodic boundary cuts the lattice and creates artificial defects.

    Returns:
        `(snapshot, info)` — `info` carries `a` (the nearest-neighbour distance),
        the box and the aspect ratio.
    """
    a = HEX_NN_OVER_D / np.sqrt(density_star)
    Lx, Ly = n_x * a, n_y * a * np.sqrt(3.0)
    basis = np.array([[0.0, 0.0], [0.5 * a, 0.5 * a * np.sqrt(3.0)]])
    cells = np.array([[i * a, j * a * np.sqrt(3.0)]
                      for i in range(n_x) for j in range(n_y)])
    pos = (cells[:, None, :] + basis[None, :, :]).reshape(-1, 2)
    pos[:, 0] = (pos[:, 0] + Lx / 2) % Lx - Lx / 2       # centre on the box
    pos[:, 1] = (pos[:, 1] + Ly / 2) % Ly - Ly / 2
    info = {"a_nn": float(a), "Lx": float(Lx), "Ly": float(Ly),
            "n_particles": int(pos.shape[0]),
            "aspect": float(Lx / Ly),
            "density_star": float(pos.shape[0] / (Lx * Ly))}
    return _make_2d_snapshot(hoomd, pos, Lx, Ly), info


def hex_tiling_for(n_target: int) -> tuple[int, int]:
    """The `(n_x, n_y)` with `N = 2 n_x n_y ≈ n_target` whose box is closest to
    square.

    The aspect ratio is `Lx/Ly = n_x/(n_y√3)`, so `n_x ≈ √3 n_y` gives a square.
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
