"""S5 — 거듭제곱 쌍 포텐셜 (`A/r^p`) + 2D 구조 분석.

## 허용오차의 출처

이 파일의 검사는 대부분 **해석적 정답이 있다** — 통계오차가 아니라 부동소수점·보간
오차가 한계다. 그래서 허용치를 `1e-12`(해석 항등식)와 보간 격자에서 유도한 값으로
나눠 쓴다. 관측값에 맞춰 재단하지 않는다 (conftest 규칙 1).

## 정답이 있는 두 배치

| 배치 | `ψ₆` | 결함 | `S(k)` 6겹 변조 |
|---|---|---|---|
| 완벽 육방격자 | **1** | **0** | **1** |
| 무작위 | ~0 | 많음 | ~0 |

이 둘이 분석기의 상·하한을 고정한다.
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
# 길이 규약 — 7 % 를 혼동하면 Γ 가 23 % 틀린다
# =============================================================================
def test_hex_nn_distance_over_d():
    """`a = (2/√3)^{1/2} d`. `d = n^{-1/2}` 이지 최근접거리가 아니다."""
    assert HEX_NN_OVER_D == pytest.approx((2 / math.sqrt(3)) ** 0.5, rel=1e-15)
    assert HEX_NN_OVER_D == pytest.approx(1.074570, abs=1e-6)


def test_gamma_is_23_percent_wrong_if_lengths_are_confused():
    """★ 이 검사가 규약 혼동의 대가를 고정한다.

    `Γ ∝ d^{-3}` 이므로 `d` 를 최근접거리로 잘못 쓰면 `Γ` 가 `1.0746³ = 1.24` 배
    틀린다. 상경계 간격이 `55.87 → 59.88` (7 %) 이므로 **상을 오판한다.**
    """
    wrong = HEX_NN_OVER_D ** 3
    assert wrong == pytest.approx(1.2408, rel=1e-3)
    # A=10 을 잘못 환산하면 액체가 hexatic 으로 보인다
    G_right = ZAHN_GAMMA_OVER_A * 10.0
    G_wrong = G_right * wrong
    assert G_right < ZAHN_GAMMA_ISO < G_wrong


def test_square_box_gives_unit_density():
    """길이 단위가 `n^{-1/2}` 이므로 `n* = 1` 이 정의상 성립한다."""
    for N in (100, 300, 780):
        L = square_box_for(N)
        assert N / L**2 == pytest.approx(1.0, rel=1e-14)


def test_hex_lattice_has_unit_density(hoomd_mod):
    """육방 빌더도 `n* = 1` 을 내야 한다 — 규약 검산."""
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
# Zahn 상도 환산
# =============================================================================
def test_zahn_gamma_over_a_is_pi_three_halves():
    assert ZAHN_GAMMA_OVER_A == pytest.approx(math.pi ** 1.5, rel=1e-15)
    assert ZAHN_GAMMA_OVER_A == pytest.approx(5.5683, abs=1e-4)


@pytest.mark.parametrize("A,phase", [
    (0.1, "isotropic-liquid"),
    (1.0, "isotropic-liquid"),
    (10.0, "isotropic-liquid"),     # ★ Γ=55.68 < 55.87 — 결정이 아니다
    (100.0, "crystal"),
])
def test_zahn_phase_classification(A, phase):
    assert zahn_phase(A)["phase_zahn"] == phase


def test_A10_sits_on_the_transition():
    """★ `A=10` 이 hexatic→액체 전이점에서 0.3 % 안에 있다 (증류 §5)."""
    z = zahn_phase(10.0)
    assert abs(z["distance_to_iso"]) < 0.005


def test_crystal_needs_A_above_ten_point_seven_five():
    """육방 결정에 필요한 최소 `A`. `trap-drag` 카드가 이 값에 의존한다."""
    A_min = amplitude_for_gamma(ZAHN_GAMMA_MELT)
    assert A_min == pytest.approx(10.75, abs=0.01)
    assert zahn_phase(A_min * 1.001)["phase_zahn"] == "crystal"
    assert zahn_phase(A_min * 0.999)["phase_zahn"] != "crystal"


def test_zahn_values_are_marked_unreproduced():
    """`reproduced: no` 인 문헌값을 근거처럼 쓰지 않는다 (knowledge/wiki/CLAUDE.md)."""
    assert "미재현" in zahn_phase(100.0)["citation"]


# =============================================================================
# 절단 오차 — 절대 기준과 상대 기준
# =============================================================================
def test_truncation_error_reports_both_scales():
    """★ `kT` 기준과 상호작용(`A`) 기준이 크게 다르다 — 둘 다 필요하다."""
    e = pair_truncation_error(amplitude=100.0, exponent=3.0, r_cut=4.8)
    assert e["beta_u_at_rcut"] == pytest.approx(100 / 4.8**3, rel=1e-12)
    assert e["u_rel_to_nearest"] == pytest.approx(1 / 4.8**3, rel=1e-12)
    # kT 기준으로는 크고(0.9) 상호작용 기준으로는 작다(0.009)
    assert e["beta_u_at_rcut"] > 0.5
    assert e["u_rel_to_nearest"] < 0.01


def test_truncation_relative_error_is_amplitude_independent():
    """상대 오차는 `A` 에 무관하다 — 그래서 `A` 스윕에서 비교 가능한 기준이다."""
    a = pair_truncation_error(amplitude=0.1, exponent=3.0, r_cut=4.8)
    b = pair_truncation_error(amplitude=100.0, exponent=3.0, r_cut=4.8)
    assert a["u_rel_to_nearest"] == pytest.approx(b["u_rel_to_nearest"], rel=1e-12)
    assert a["beta_u_at_rcut"] != pytest.approx(b["beta_u_at_rcut"], rel=1e-3)


def test_force_truncation_is_tighter_than_energy():
    """`F ∝ r^{-(p+1)}` 이므로 힘이 더 빨리 죽는다 — 상대값으로 확인."""
    e = pair_truncation_error(amplitude=100.0, exponent=3.0, r_cut=4.8)
    assert e["f_rel_to_nearest"] < e["u_rel_to_nearest"]


# =============================================================================
# 완벽 육방격자 — 분석기의 정답지
# =============================================================================
@pytest.fixture
def perfect_hex(hoomd_mod):
    nx, ny = hex_tiling_for(300)
    snap, info = hex_2d_snapshot(hoomd_mod, n_x=nx, n_y=ny)
    pos = np.array(snap.particles.position[:, :2], dtype=np.float64)
    return pos, info


@pytest.mark.benchmark
def test_perfect_hex_has_psi6_exactly_one(perfect_hex):
    """★ 정답이 1 이다. 부동소수점 오차만 허용."""
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
    """★ `S(k)` 6겹 변조 = 1 (Bragg 점). 이것이 결정/hexatic 을 가른다."""
    pos, info = perfect_hex
    s = structure_factor(pos, Lx=info["Lx"], Ly=info["Ly"], n_max=20)
    assert s.sixfold_modulation == pytest.approx(1.0, abs=1e-6)


def test_perfect_hex_first_bragg_peak(perfect_hex):
    """첫 Bragg: `|G| = 4π/(a√3)`. 허용치는 k 격자 간격에서 나온다."""
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
# 무작위 배치 — 분석기의 하한
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
    """6겹 변조가 작아야 한다 — 유한 표본이라 정확히 0 은 아니다."""
    pos, L = random_config
    s = structure_factor(pos, Lx=L, Ly=L, n_max=20)
    assert s.sixfold_modulation < 0.25


def test_random_respects_min_sep(random_config):
    pos, L = random_config
    assert min_separation(pos, Lx=L, Ly=L) >= 0.4 - 1e-12


def test_random_builder_refuses_impossible_density(hoomd_mod):
    """겹친 배치를 조용히 돌려주지 않는다 — 폭발의 원인이 된다."""
    with pytest.raises(RuntimeError, match="넣지 못했다"):
        random_2d_snapshot(hoomd_mod, n=300, box=square_box_for(300),
                           min_sep=1.5, seed=1, max_tries=300)


# =============================================================================
# 국소 vs 전역 ψ₆ — 다결정을 구별한다
# =============================================================================
def test_local_psi6_is_at_least_global(perfect_hex, random_config):
    """`|⟨ψ₆ᵢ⟩| ≤ ⟨|ψ₆ᵢ|⟩` 는 삼각부등식이므로 **항등적으로** 성립한다.

    깨지면 계산이 틀렸다는 뜻이다.
    """
    for pos, box in ((perfect_hex[0], (perfect_hex[1]["Lx"], perfect_hex[1]["Ly"])),
                     (random_config[0], (random_config[1], random_config[1]))):
        h = hex_order(pos, Lx=box[0], Ly=box[1])
        assert h.psi6_global <= h.psi6_local_mean + 1e-12


def test_rotating_a_periodic_crystal_in_a_fixed_box_creates_real_defects():
    """★ 주기 상자 안에서 결정을 **회전시킬 수 없다** — 경계와 정합성이 깨진다.

    처음에 "국소 `ψ₆` 는 회전 불변" 을 확인하려고 30° 회전을 걸었더니 `1.0 → 0.807`
    로 떨어졌다. `hex_order` 의 버그가 아니라 **회전이 경계에서 진짜 결함을 만든** 것이다
    (격자 주기가 상자 주기와 더 이상 맞지 않는다).

    이것을 테스트로 고정하는 이유: 다결정·결정면 방향(진동 카드 §A9)을 다룰 때
    "상자 안에서 격자를 돌리면 된다"는 실수를 하기 쉽다. **상자를 격자에 맞춰
    다시 만들어야 한다.**
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
    assert spoiled.psi6_local_mean < 0.9, "회전이 결함을 만들지 않았다 — 예상과 다르다"
    assert spoiled.defect_fraction > 0.0


# =============================================================================
# 이웃 정의가 답을 바꾼다
# =============================================================================
def test_neighbor_mode_matters_for_disordered_configs(random_config):
    """★ Voronoi 와 nearest6 은 다른 답을 준다 — 어느 것을 썼는지 기록해야 한다."""
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
# 다중 프레임
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
# 상자 모양 — 초기조건 비교의 교란 요인
# =============================================================================
#  근거: knowledge/wiki/findings/box-shape-confounds-initial-condition-comparison.md
def test_box_shape_is_decoupled_from_initial_condition(hoomd_mod):
    """★ `random` 초기배치가 육방정합 상자를 쓸 수 있어야 한다.

    묶여 있으면 hex vs random 비교가 상자 종횡비와 `r_cut` 까지 함께 바꾼다.
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
    """`auto` 가 초기조건에 따라 갈린다는 사실을 고정한다 — 비교할 때 명시해야 한다."""
    from simbot.run import Soft2DRunConfig, build_soft2d
    _, _, h = build_soft2d(hoomd_mod, Soft2DRunConfig(init="hex", n_particles=100))
    _, _, r = build_soft2d(hoomd_mod, Soft2DRunConfig(init="random", n_particles=100))
    assert h["box_shape"] == "hex_commensurate" and r["box_shape"] == "square"
    assert h["aspect"] != pytest.approx(r["aspect"], rel=1e-3)
    assert h["r_cut"] != pytest.approx(r["r_cut"], rel=1e-3)


def test_hex_init_in_square_box_is_refused(hoomd_mod):
    """육방 배치는 정합 상자에서만 완벽하다 — 정사각은 격자를 끊는다."""
    from simbot.run import Soft2DRunConfig, build_soft2d
    with pytest.raises(ValueError, match="육방정합 상자에서만"):
        build_soft2d(hoomd_mod, Soft2DRunConfig(init="hex", box_shape="square"))


def test_unknown_box_shape_is_refused(hoomd_mod):
    from simbot.run import Soft2DRunConfig, build_soft2d
    with pytest.raises(ValueError, match="box_shape"):
        build_soft2d(hoomd_mod, Soft2DRunConfig(box_shape="hexagonalish"))


def test_box_shape_is_recorded_in_manifest_info(hoomd_mod):
    """나중에 '이 두 런을 비교할 수 있는가' 를 검사하려면 기록이 있어야 한다."""
    from simbot.run import Soft2DRunConfig, build_soft2d
    _, _, info = build_soft2d(hoomd_mod, Soft2DRunConfig(init="random"))
    assert info["box_shape"] in ("square", "hex_commensurate")


def test_random_builder_accepts_rectangular_box(hoomd_mod):
    """직사각 상자에서도 최소분리를 지켜야 한다 (성분별 최소이미지)."""
    snap = random_2d_snapshot(hoomd_mod, n=100, box_x=10.7457, box_y=9.3060,
                              min_sep=0.5, seed=1, max_tries=20000)
    pos = np.array(snap.particles.position[:, :2], dtype=np.float64)
    assert min_separation(pos, Lx=10.7457, Ly=9.3060) >= 0.5 - 1e-12


def test_random_builder_needs_a_box(hoomd_mod):
    with pytest.raises(ValueError, match="box"):
        random_2d_snapshot(hoomd_mod, n=10)


# =============================================================================
# 재현성 — 같은 config 가 같은 궤적을 주는가
# =============================================================================
#  ★ 이것이 재현성 주장의 **유일한 직접 증거**다. `code_hash`·`git_rev`·`env_hash` 는
#    "무엇이 만들었는지" 를 기록할 뿐 "다시 돌리면 같은가" 를 증명하지 않는다.
#    CLAUDE.md 는 HOOMD `Brownian` 이 같은 seed 에서 비트 재현된다고 적고 있는데,
#    그 주장이 실제로 성립하는지는 실행해 봐야 안다.
@pytest.mark.slow
def test_same_config_reproduces_the_trajectory_bit_for_bit(tmp_path):
    """같은 `Soft2DRunConfig` 두 번 → 궤적이 **바이트 단위로** 같다.

    ★ 청크 경계에 무관해야 한다: HOOMD 의 counter-based RNG 는 `(timestep, tag, seed)`
      로 키를 만들므로 `run()` 을 몇 번에 나눠 부르든 같은 스텝에서 같은 난수가 나온다.
      그래서 `equil_tau` 를 바꿔 표집 구간을 옮겨도 같은 절대 시각의 좌표는 같다.
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
        assert za[key].tobytes() == zb[key].tobytes(), f"{key} 가 재현되지 않는다"

    # 파생 스칼라도 비트 단위로 같아야 한다 (부동소수점 누적 순서까지)
    assert a["energy_per_particle"] == b["energy_per_particle"]
    assert a["guards"]["min_separation"] == b["guards"]["min_separation"]
    # run_hash 는 config 만의 함수다 — 같은 config 면 같아야 한다
    assert a["manifest"]["run_hash"] == b["manifest"]["run_hash"] == cfg.hash()


@pytest.mark.slow
def test_run_hash_changes_when_any_physical_knob_changes(tmp_path):
    """`run_hash` 가 파라미터를 실제로 구별하는가 — 안 하면 다른 런이 같아 보인다."""
    from simbot.run import Soft2DRunConfig

    base = Soft2DRunConfig(amplitude=1.0, n_particles=64, seed=17)
    for field, value in (("amplitude", 2.0), ("seed", 18), ("dt_star", 1e-4),
                         ("n_particles", 65), ("init", "hex"),
                         ("box_shape", "hex_commensurate"), ("r_min", 0.15),
                         ("prod_tau", 99.0), ("n_frames", 7),
                         ("min_sep_init", 0.75), ("nlist_buffer", 0.2),
                         ("exponent", 4.0), ("density_star", 0.9)):
        #  ★ 기본값과 같은 값을 넣으면 검사가 조용히 무력해진다 (실제로 겪었다:
        #    `r_min=0.2` 가 기본값이라 "바꿨는데 해시가 같다" 로 오진했다)
        assert getattr(base, field) != value, (
            f"{field}={value!r} 가 기본값과 같다 — 이 행은 아무것도 검사하지 않는다")
        other = replace(base, **{field: value})
        assert other.hash() != base.hash(), f"{field} 가 run_hash 에 안 들어간다"


def test_manifest_records_the_libraries_that_determine_the_numbers():
    """★ `hoomd` 만 기록하면 부족하다 — `freud` 가 Voronoi·`ψ₆` 를 계산한다.

    버전이 올라가면 같은 궤적에서 다른 측정값이 나올 수 있다. 기록이 없으면
    그 사실을 사후에 알 수 없다.
    """
    from simbot.io import ENV_PACKAGES, env_versions

    env = env_versions()
    for name in ("hoomd", "freud", "numpy", "scipy"):
        assert name in ENV_PACKAGES, f"{name} 이 ENV_PACKAGES 에 없다"
        assert env[name] not in ("absent", "unknown"), f"{name} 버전을 못 읽는다"
