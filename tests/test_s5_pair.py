"""S5 — the power-law pair potential (`A/r^p`) + 2D structure analysis.

## Where the tolerances come from

Most checks in this file **have an analytic answer** -- the limit is floating-point
and interpolation error, not statistical error. So the tolerances split into `1e-12`
(analytic identities) and values derived from the interpolation grid. None are cut
to fit the observation (conftest rule 1).

## The two configurations with known answers

| configuration | `ψ₆` | defects | `S(k)` 6-fold modulation |
|---|---|---|---|
| perfect hexagonal lattice | **1** | **0** | **1** |
| random | ~0 | many | ~0 |

These two pin the analyser's upper and lower bounds.
"""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from simbot.analysis.structure import (ZAHN_GAMMA_ISO, ZAHN_GAMMA_MELT,
                                       ZAHN_GAMMA_OVER_A, amplitude_for_gamma,
                                       hex_order, min_separation, rdf,
                                       structure_factor, zahn_phase)
from simbot.build import (HEX_NN_OVER_D, hex_2d_snapshot, hex_tiling_for,
                          random_2d_snapshot, square_box_for)
from simbot.forces import pair_truncation_error


# =============================================================================
# The length convention — confuse 7 % and Γ is wrong by 23 %
# =============================================================================
def test_hex_nn_distance_over_d():
    """`a = (2/√3)^{1/2} d`. `d = n^{-1/2}`, not the nearest-neighbour distance."""
    assert HEX_NN_OVER_D == pytest.approx((2 / math.sqrt(3)) ** 0.5, rel=1e-15)
    assert HEX_NN_OVER_D == pytest.approx(1.074570, abs=1e-6)


def test_gamma_is_23_percent_wrong_if_lengths_are_confused():
    """★ This check pins the price of confusing the convention.

    `Γ ∝ d^{-3}`, so wrongly using the nearest-neighbour distance for `d` makes `Γ`
    wrong by `1.0746³ = 1.24`. The phase boundaries are `55.87 → 59.88` apart (7 %),
    so **the phase is misread.**
    """
    wrong = HEX_NN_OVER_D ** 3
    assert wrong == pytest.approx(1.2408, rel=1e-3)
    # convert A=10 wrongly and a liquid looks like a hexatic
    G_right = ZAHN_GAMMA_OVER_A * 10.0
    G_wrong = G_right * wrong
    assert G_right < ZAHN_GAMMA_ISO < G_wrong


def test_square_box_gives_unit_density():
    """The length unit is `n^{-1/2}`, so `n* = 1` holds by definition."""
    for N in (100, 300, 780):
        L = square_box_for(N)
        assert N / L**2 == pytest.approx(1.0, rel=1e-14)


def test_hex_lattice_has_unit_density(hoomd_mod):
    """The hexagonal builder has to give `n* = 1` too — a convention check."""
    for target in (100, 300):
        nx, ny = hex_tiling_for(target)
        _, info = hex_2d_snapshot(hoomd_mod, n_x=nx, n_y=ny)
        assert info["density_star"] == pytest.approx(1.0, rel=1e-12)
        assert info["a_nn"] == pytest.approx(HEX_NN_OVER_D, rel=1e-12)


def test_hex_tiling_is_near_square(hoomd_mod):
    nx, ny = hex_tiling_for(300)
    _, info = hex_2d_snapshot(hoomd_mod, n_x=nx, n_y=ny)
    assert 0.7 < info["aspect"] < 1.4


# =============================================================================
# Conversion to Zahn's phase diagram
# =============================================================================
def test_zahn_gamma_over_a_is_pi_three_halves():
    assert ZAHN_GAMMA_OVER_A == pytest.approx(math.pi ** 1.5, rel=1e-15)
    assert ZAHN_GAMMA_OVER_A == pytest.approx(5.5683, abs=1e-4)


@pytest.mark.parametrize("A,phase", [
    (0.1, "isotropic-liquid"),
    (1.0, "isotropic-liquid"),
    (10.0, "isotropic-liquid"),     # ★ Γ=55.68 < 55.87 — not a crystal
    (100.0, "crystal"),
])
def test_zahn_phase_classification(A, phase):
    assert zahn_phase(A)["phase_zahn"] == phase


def test_A10_sits_on_the_transition():
    """★ `A=10` sits within 0.3 % of the hexatic→liquid transition (distillation
    §5)."""
    z = zahn_phase(10.0)
    assert abs(z["distance_to_iso"]) < 0.005


def test_crystal_needs_A_above_ten_point_seven_five():
    """The minimum `A` needed for a hexagonal crystal. The `trap-drag` card depends
    on this value."""
    A_min = amplitude_for_gamma(ZAHN_GAMMA_MELT)
    assert A_min == pytest.approx(10.75, abs=0.01)
    assert zahn_phase(A_min * 1.001)["phase_zahn"] == "crystal"
    assert zahn_phase(A_min * 0.999)["phase_zahn"] != "crystal"


def test_zahn_values_are_marked_unreproduced():
    """A literature value with `reproduced: no` is not used as a basis
    (knowledge/wiki/CLAUDE.md)."""
    assert "not reproduced" in zahn_phase(100.0)["citation"]


# =============================================================================
# Truncation error — the absolute and the relative basis
# =============================================================================
def test_truncation_error_reports_both_scales():
    """★ The `kT` basis and the interaction (`A`) basis differ greatly — both are
    needed."""
    e = pair_truncation_error(amplitude=100.0, exponent=3.0, r_cut=4.8)
    assert e["beta_u_at_rcut"] == pytest.approx(100 / 4.8**3, rel=1e-12)
    assert e["u_rel_to_nearest"] == pytest.approx(1 / 4.8**3, rel=1e-12)
    # large on the kT basis (0.9) and small on the interaction basis (0.009)
    assert e["beta_u_at_rcut"] > 0.5
    assert e["u_rel_to_nearest"] < 0.01


def test_truncation_relative_error_is_amplitude_independent():
    """The relative error is independent of `A` — which makes it the basis that is
    comparable across an `A` sweep."""
    a = pair_truncation_error(amplitude=0.1, exponent=3.0, r_cut=4.8)
    b = pair_truncation_error(amplitude=100.0, exponent=3.0, r_cut=4.8)
    assert a["u_rel_to_nearest"] == pytest.approx(b["u_rel_to_nearest"], rel=1e-12)
    assert a["beta_u_at_rcut"] != pytest.approx(b["beta_u_at_rcut"], rel=1e-3)


def test_force_truncation_is_tighter_than_energy():
    """`F ∝ r^{-(p+1)}`, so the force dies faster — checked on the relative value."""
    e = pair_truncation_error(amplitude=100.0, exponent=3.0, r_cut=4.8)
    assert e["f_rel_to_nearest"] < e["u_rel_to_nearest"]


# =============================================================================
# The perfect hexagonal lattice — the analyser's answer key
# =============================================================================
@pytest.fixture
def perfect_hex(hoomd_mod):
    nx, ny = hex_tiling_for(300)
    snap, info = hex_2d_snapshot(hoomd_mod, n_x=nx, n_y=ny)
    pos = np.array(snap.particles.position[:, :2], dtype=np.float64)
    return pos, info


@pytest.mark.benchmark
def test_perfect_hex_has_psi6_exactly_one(perfect_hex):
    """★ The answer is 1. Only floating-point error is allowed."""
    pos, info = perfect_hex
    h = hex_order(pos, Lx=info["Lx"], Ly=info["Ly"])
    assert h.psi6_global == pytest.approx(1.0, abs=1e-9)
    assert h.psi6_local_mean == pytest.approx(1.0, abs=1e-9)


@pytest.mark.benchmark
def test_perfect_hex_has_no_defects(perfect_hex):
    pos, info = perfect_hex
    h = hex_order(pos, Lx=info["Lx"], Ly=info["Ly"])
    assert h.defect_fraction == 0.0
    assert h.coordination_hist == {6: 1.0}


@pytest.mark.benchmark
def test_perfect_hex_sixfold_modulation_is_one(perfect_hex):
    """★ `S(k)` 6-fold modulation = 1 (Bragg points). This is what separates a
    crystal from a hexatic."""
    pos, info = perfect_hex
    s = structure_factor(pos, Lx=info["Lx"], Ly=info["Ly"], n_max=20)
    assert s.sixfold_modulation == pytest.approx(1.0, abs=1e-6)


def test_perfect_hex_first_bragg_peak(perfect_hex):
    """The first Bragg peak: `|G| = 4π/(a√3)`. The tolerance comes from the k-grid
    spacing."""
    pos, info = perfect_hex
    s = structure_factor(pos, Lx=info["Lx"], Ly=info["Ly"], n_max=20, bins=120)
    expect = 4 * math.pi / (info["a_nn"] * math.sqrt(3))
    dk = s.k_radial[1] - s.k_radial[0]
    assert abs(s.first_peak_k - expect) <= dk


def test_perfect_hex_rdf_first_peak_at_nn(perfect_hex):
    pos, info = perfect_hex
    r = rdf(pos, Lx=info["Lx"], Ly=info["Ly"], bins=400)
    dr = r.r[1] - r.r[0]
    assert abs(r.first_peak_r - info["a_nn"]) <= 2 * dr


def test_perfect_hex_min_separation_is_the_nn_distance(perfect_hex):
    pos, info = perfect_hex
    assert min_separation(pos, Lx=info["Lx"], Ly=info["Ly"]) == pytest.approx(
        info["a_nn"], rel=1e-10)


# =============================================================================
# Random placement — the analyser's lower bound
# =============================================================================
@pytest.fixture
def random_config(hoomd_mod):
    L = square_box_for(300)
    snap = random_2d_snapshot(hoomd_mod, n=300, box=L, min_sep=0.4, seed=3,
                             max_tries=20000)
    return np.array(snap.particles.position[:, :2], dtype=np.float64), L


def test_random_has_low_psi6(random_config):
    pos, L = random_config
    h = hex_order(pos, Lx=L, Ly=L)
    assert h.psi6_global < 0.15
    assert h.defect_fraction > 0.4


def test_random_is_nearly_isotropic(random_config):
    """The 6-fold modulation has to be small — not exactly 0, since the sample is
    finite."""
    pos, L = random_config
    s = structure_factor(pos, Lx=L, Ly=L, n_max=20)
    assert s.sixfold_modulation < 0.25


def test_random_respects_min_sep(random_config):
    pos, L = random_config
    assert min_separation(pos, Lx=L, Ly=L) >= 0.4 - 1e-12


def test_random_builder_refuses_impossible_density(hoomd_mod):
    """It does not quietly return an overlapped placement — that causes a blow-up."""
    with pytest.raises(RuntimeError, match="could not place"):
        random_2d_snapshot(hoomd_mod, n=300, box=square_box_for(300),
                           min_sep=1.5, seed=1, max_tries=300)


# =============================================================================
# Local vs global ψ₆ — distinguishing a polycrystal
# =============================================================================
def test_local_psi6_is_at_least_global(perfect_hex, random_config):
    """`|⟨ψ₆ᵢ⟩| ≤ ⟨|ψ₆ᵢ|⟩` holds **identically**, by the triangle inequality.

    Breaking it means the computation is wrong.
    """
    for pos, box in ((perfect_hex[0], (perfect_hex[1]["Lx"], perfect_hex[1]["Ly"])),
                     (random_config[0], (random_config[1], random_config[1]))):
        h = hex_order(pos, Lx=box[0], Ly=box[1])
        assert h.psi6_global <= h.psi6_local_mean + 1e-12


def test_rotating_a_periodic_crystal_in_a_fixed_box_creates_real_defects():
    """★ A crystal **cannot be rotated** inside a periodic box — commensurability
    with the boundary breaks.

    A 30° rotation was applied at first to confirm that "local `ψ₆` is
    rotation-invariant", and it fell from `1.0` to `0.807`. That is not a bug in
    `hex_order` -- **the rotation creates real defects at the boundary** (the
    lattice period no longer matches the box period).

    Why this is pinned as a test: when handling polycrystals or crystal-plane
    orientation (oscillation card §A9) it is easy to make the mistake of "just
    rotate the lattice inside the box". **The box has to be rebuilt to match the
    lattice.**
    """
    import hoomd
    snap, info = hex_2d_snapshot(hoomd, n_x=6, n_y=3)
    pos = np.array(snap.particles.position[:, :2], dtype=np.float64)
    clean = hex_order(pos, Lx=info["Lx"], Ly=info["Ly"])
    assert clean.psi6_local_mean == pytest.approx(1.0, abs=1e-9)
    assert clean.defect_fraction == 0.0

    th = np.deg2rad(30.0)
    R = np.array([[np.cos(th), -np.sin(th)], [np.sin(th), np.cos(th)]])
    rot = pos @ R.T
    spoiled = hex_order(rot, Lx=info["Lx"], Ly=info["Ly"])
    assert spoiled.psi6_local_mean < 0.9, \
        "the rotation created no defects — that is not what was expected"
    assert spoiled.defect_fraction > 0.0


# =============================================================================
# The neighbour definition changes the answer
# =============================================================================
def test_neighbor_mode_matters_for_disordered_configs(random_config):
    """★ Voronoi and nearest6 give different answers — which one was used has to be
    recorded."""
    pos, L = random_config
    v = hex_order(pos, Lx=L, Ly=L, neighbor_mode="voronoi")
    n6 = hex_order(pos, Lx=L, Ly=L, neighbor_mode="nearest6")
    assert v.neighbor_mode == "voronoi" and n6.neighbor_mode == "nearest6"
    assert v.psi6_local_mean != pytest.approx(n6.psi6_local_mean, rel=1e-3)


def test_unknown_neighbor_mode_raises(random_config):
    pos, L = random_config
    with pytest.raises(ValueError, match="voronoi"):
        hex_order(pos, Lx=L, Ly=L, neighbor_mode="guess")


# =============================================================================
# Multiple frames
# =============================================================================
def test_multi_frame_averaging(perfect_hex):
    pos, info = perfect_hex
    frames = np.stack([pos, pos, pos])
    h = hex_order(frames, Lx=info["Lx"], Ly=info["Ly"])
    assert h.n_frames == 3
    assert h.psi6_global == pytest.approx(1.0, abs=1e-9)
    assert math.isnan(h.psi6_spread) or h.psi6_spread < 1e-12


def test_single_frame_accepts_2d_array(perfect_hex):
    pos, info = perfect_hex
    h = hex_order(pos, Lx=info["Lx"], Ly=info["Ly"])
    assert h.n_frames == 1


# =============================================================================
# Box shape — a confounder in initial-condition comparisons
# =============================================================================
#  Basis: knowledge/wiki/findings/box-shape-confounds-initial-condition-comparison.md
def test_box_shape_is_decoupled_from_initial_condition(hoomd_mod):
    """★ A `random` initial placement has to be able to use a hex-commensurate box.

    Tied together, a hex vs random comparison also changes the box aspect ratio and
    `r_cut`.
    """
    from simbot.run import Soft2DRunConfig, build_soft2d
    _, _, hex_info = build_soft2d(hoomd_mod, Soft2DRunConfig(
        amplitude=1.0, init="hex", n_particles=100))
    _, _, rnd_info = build_soft2d(hoomd_mod, Soft2DRunConfig(
        amplitude=1.0, init="random", box_shape="hex_commensurate",
        n_particles=100))
    assert rnd_info["Lx"] == pytest.approx(hex_info["Lx"], rel=1e-12)
    assert rnd_info["Ly"] == pytest.approx(hex_info["Ly"], rel=1e-12)
    assert rnd_info["r_cut"] == pytest.approx(hex_info["r_cut"], rel=1e-12)


def test_auto_box_shape_differs_by_init_and_that_is_the_trap(hoomd_mod):
    """Pins the fact that `auto` branches on the initial condition — it has to be
    stated explicitly when comparing."""
    from simbot.run import Soft2DRunConfig, build_soft2d
    _, _, h = build_soft2d(hoomd_mod, Soft2DRunConfig(init="hex", n_particles=100))
    _, _, r = build_soft2d(hoomd_mod, Soft2DRunConfig(init="random", n_particles=100))
    assert h["box_shape"] == "hex_commensurate" and r["box_shape"] == "square"
    assert h["aspect"] != pytest.approx(r["aspect"], rel=1e-3)
    assert h["r_cut"] != pytest.approx(r["r_cut"], rel=1e-3)


def test_hex_init_in_square_box_is_refused(hoomd_mod):
    """A hexagonal placement is only perfect in a commensurate box — a square one
    cuts the lattice."""
    from simbot.run import Soft2DRunConfig, build_soft2d
    with pytest.raises(ValueError, match="only perfect in a hex-commensurate box"):
        build_soft2d(hoomd_mod, Soft2DRunConfig(init="hex", box_shape="square"))


def test_unknown_box_shape_is_refused(hoomd_mod):
    from simbot.run import Soft2DRunConfig, build_soft2d
    with pytest.raises(ValueError, match="box_shape"):
        build_soft2d(hoomd_mod, Soft2DRunConfig(box_shape="hexagonalish"))


def test_box_shape_is_recorded_in_manifest_info(hoomd_mod):
    """Checking later whether 'these two runs are comparable' requires the record."""
    from simbot.run import Soft2DRunConfig, build_soft2d
    _, _, info = build_soft2d(hoomd_mod, Soft2DRunConfig(init="random"))
    assert info["box_shape"] in ("square", "hex_commensurate")


def test_random_builder_accepts_rectangular_box(hoomd_mod):
    """The minimum separation has to hold in a rectangular box too (per-component
    minimum image)."""
    snap = random_2d_snapshot(hoomd_mod, n=100, box_x=10.7457, box_y=9.3060,
                              min_sep=0.5, seed=1, max_tries=20000)
    pos = np.array(snap.particles.position[:, :2], dtype=np.float64)
    assert min_separation(pos, Lx=10.7457, Ly=9.3060) >= 0.5 - 1e-12


def test_random_builder_needs_a_box(hoomd_mod):
    with pytest.raises(ValueError, match="box"):
        random_2d_snapshot(hoomd_mod, n=10)


# =============================================================================
# Reproducibility — does the same config give the same trajectory
# =============================================================================
#  ★ This is the **only direct evidence** for the reproducibility claim.
#    `code_hash`, `git_rev` and `env_hash` record "what made it" and prove nothing
#    about "is it the same when re-run".
#    CLAUDE.md states that HOOMD `Brownian` reproduces bit-for-bit on the same seed,
#    and whether that claim actually holds can only be known by running it.
@pytest.mark.slow
def test_same_config_reproduces_the_trajectory_bit_for_bit(tmp_path):
    """The same `Soft2DRunConfig` twice → **byte-identical** trajectories.

    ★ It has to be independent of the chunk boundaries: HOOMD's counter-based RNG
      keys on `(timestep, tag, seed)`, so however many calls `run()` is split into,
      the same step gets the same random numbers. That is why moving the sampling
      window by changing `equil_tau` leaves the coordinates at the same absolute
      time identical.
    """
    from simbot.run import Soft2DRunConfig, run_soft2d

    cfg = Soft2DRunConfig(amplitude=1.0, n_particles=64, init="random",
                          box_shape="square", dt_star=4.5e-4, equil_tau=0.0,
                          prod_tau=0.2, n_frames=8, seed=17, label="repro")
    a = run_soft2d(cfg, outdir=tmp_path / "a")
    b = run_soft2d(cfg, outdir=tmp_path / "b")

    za = np.load(tmp_path / "a" / "samples.npz")
    zb = np.load(tmp_path / "b" / "samples.npz")
    for key in ("traj", "energy", "max_force", "init_pos", "box"):
        assert za[key].tobytes() == zb[key].tobytes(), f"{key} does not reproduce"

    # the derived scalars have to match bitwise too (down to accumulation order)
    assert a["energy_per_particle"] == b["energy_per_particle"]
    assert a["guards"]["min_separation"] == b["guards"]["min_separation"]
    # run_hash is a function of the config alone — the same config must give the same
    assert a["manifest"]["run_hash"] == b["manifest"]["run_hash"] == cfg.hash()


@pytest.mark.slow
def test_run_hash_changes_when_any_physical_knob_changes(tmp_path):
    """Does `run_hash` actually distinguish the parameters — if not, different runs
    look identical."""
    from simbot.run import Soft2DRunConfig

    base = Soft2DRunConfig(amplitude=1.0, n_particles=64, seed=17)
    for field, value in (("amplitude", 2.0), ("seed", 18), ("dt_star", 1e-4),
                         ("n_particles", 65), ("init", "hex"),
                         ("box_shape", "hex_commensurate"), ("r_min", 0.15),
                         ("prod_tau", 99.0), ("n_frames", 7),
                         ("min_sep_init", 0.75), ("nlist_buffer", 0.2),
                         ("exponent", 4.0), ("density_star", 0.9)):
        #  ★ Passing the same value as the default makes the check quietly toothless
        #    (this happened: `r_min=0.2` is the default, and it was misdiagnosed as
        #    "changed it and the hash is the same")
        assert getattr(base, field) != value, (
            f"{field}={value!r} equals the default — this row checks nothing")
        other = replace(base, **{field: value})
        assert other.hash() != base.hash(), f"{field} is not in run_hash"


def test_manifest_records_the_libraries_that_determine_the_numbers():
    """★ Recording `hoomd` alone is not enough — `freud` computes the Voronoi
    tessellation and `ψ₆`.

    A version bump can give a different measurement from the same trajectory. With
    no record, that cannot be known after the fact.
    """
    from simbot.io import ENV_PACKAGES, env_versions

    env = env_versions()
    for name in ("hoomd", "freud", "numpy", "scipy"):
        assert name in ENV_PACKAGES, f"{name} is not in ENV_PACKAGES"
        assert env[name] not in ("absent", "unknown"), \
            f"{name}'s version cannot be read"
