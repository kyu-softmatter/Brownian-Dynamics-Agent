"""S7 — 시간분해 구조 분석. `freud` 는 쓰지만 HOOMD 는 쓰지 않는다 (빠름).

## 이 파일이 지키는 것

**① 시간분해 판본과 시간평균 판본이 갈라지지 않는다.** `hex_order(frames)` 와
`hex_order_series(frames).*.mean()` 은 **같은 수여야 한다.** 두 코드 경로가
조용히 갈라지면 "후반 절반 평균" 과 "시계열 후반부" 가 다른 값을 내고, 어느 쪽이
맞는지 알 수 없게 된다.

**② 극한에서 정답이 있다.** 완벽 육방격자 → `ψ₆ = 1`, 결함 0, 배위수 전부 6.
무작위 → `ψ₆ ≈ 0`, 배위수 분포가 퍼진다. 두 극한을 프레임 단위로 검사한다.

**③ 커버리지 기하는 해석적으로 검산된다** — `φ = (π/4)(σ/d)²` 는 왕복 가능하다.
"""
from __future__ import annotations

import numpy as np
import pytest

freud = pytest.importorskip("freud")

from simbot.analysis.structure import (HexOrderSeries, hex_order,
                                       hex_order_series, rdf_windows,
                                       voronoi_frame)
from simbot.build import (HEX_NN_OVER_D, box_si_for_coverage,
                          coverage_from_sigma_over_d,
                          sigma_over_d_for_coverage, square_box_for)


# --- 완벽 육방격자 (freud 좌표계, density_star = 1) ---------------------------
def _perfect_hex(n_x: int = 10, n_y: int = 10) -> tuple[np.ndarray, float, float]:
    """`n* = 1` 인 완벽 육방격자. `simbot.build.hex_2d_snapshot` 과 같은 기하."""
    a = HEX_NN_OVER_D                      # 최근접거리 (단위 d)
    Lx, Ly = n_x * a, n_y * a * np.sqrt(3.0) / 2.0
    pts = []
    for j in range(n_y):
        for i in range(n_x):
            x = (i + 0.5 * (j % 2)) * a
            y = j * a * np.sqrt(3.0) / 2.0
            pts.append((x - Lx / 2, y - Ly / 2))
    return np.array(pts), Lx, Ly


def _random_frame(n: int, L: float, seed: int) -> np.ndarray:
    rng = np.random.default_rng(seed)
    return rng.uniform(-L / 2, L / 2, size=(n, 2))


# =============================================================================
# ① 두 코드 경로가 갈라지지 않는다
# =============================================================================
def test_series_mean_equals_frame_averaged_hex_order():
    """`hex_order_series` 의 프레임 평균 == `hex_order`. **비트 수준으로.**

    같은 `_frame_hex` 를 쓰므로 부동소수점 차이만 허용한다. 이 assert 가 깨지면
    두 함수가 다른 이웃 정의나 다른 정규화를 쓰기 시작한 것이다.
    """
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts, _random_frame(100, min(Lx, Ly), 3),
                       hexpts * 0.999])
    agg = hex_order(frames, Lx=Lx, Ly=Ly)
    ser = hex_order_series(frames, Lx=Lx, Ly=Ly)

    assert ser.n_frames == 3
    assert ser.psi6_global.mean() == pytest.approx(agg.psi6_global, abs=1e-12)
    assert ser.psi6_local.mean() == pytest.approx(agg.psi6_local_mean, abs=1e-12)
    assert ser.defect_fraction.mean() == pytest.approx(agg.defect_fraction,
                                                       abs=1e-12)
    # 프레임별 값도 일치해야 한다 (평균만 맞는 것으로는 부족하다)
    assert ser.psi6_global == pytest.approx(np.array(agg.psi6_per_frame),
                                            abs=1e-12)


def test_series_t_star_length_is_checked():
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts, hexpts])
    with pytest.raises(ValueError, match="t_star 길이"):
        hex_order_series(frames, Lx=Lx, Ly=Ly, t_star=np.array([0.0, 1.0, 2.0]))


# =============================================================================
# ② 극한의 정답 — 프레임 단위로
# =============================================================================
def test_perfect_lattice_series_is_exact_every_frame():
    """완벽 격자에서는 **모든 프레임이** `ψ₆ = 1`, 결함 0, 배위수 6 이다."""
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts] * 4)
    s = hex_order_series(frames, Lx=Lx, Ly=Ly)

    assert s.psi6_global == pytest.approx(np.ones(4), abs=1e-9)
    assert s.psi6_local == pytest.approx(np.ones(4), abs=1e-9)
    assert np.all(s.defect_fraction == 0.0)
    six = s.coord_fraction[:, s.coord_labels == 6].ravel()
    assert six == pytest.approx(np.ones(4), abs=1e-12)
    assert np.all(s.coord_kinds == 1)              # 배위수 6 한 종류뿐
    # 5·7 이 아예 없으면 불균형은 정의상 0 (0/0 을 0 으로 처리한다)
    assert np.all(s.five_seven_balance == 0.0)


def test_random_frames_have_no_orientational_order_and_spread_coordination():
    """무작위 배치: `ψ₆` 전역이 유한크기 바닥 근처, 배위수는 여러 종류로 퍼진다."""
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(6)])
    s = hex_order_series(frames, Lx=L, Ly=L)

    # 유한크기 바닥: 무상관 위상 100개면 |⟨ψ₆⟩| ~ 1/√N = 0.1 규모다
    assert np.all(s.psi6_global < 4.0 / np.sqrt(100))
    assert np.all(s.defect_fraction > 0.3)
    assert np.all(s.coord_kinds >= 4), f"배위수 종류 {s.coord_kinds}"


def test_defect_fraction_counts_coordination_outside_the_histogram_range():
    """`coord_range` 를 좁혀도 `defect_fraction` 은 전량으로 센다.

    히스토그램 칸에서 빠진 배위수(예: 11) 가 결함 분율에서도 빠지면
    **결함이 조용히 사라진다.** 그 조용한 손실을 막는 assert 다.
    """
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(4)])
    wide = hex_order_series(frames, Lx=L, Ly=L, coord_range=(3, 12))
    narrow = hex_order_series(frames, Lx=L, Ly=L, coord_range=(6, 6))

    assert narrow.defect_fraction == pytest.approx(wide.defect_fraction,
                                                   abs=1e-12)
    assert narrow.coord_fraction.shape[1] == 1
    # 좁은 범위에서는 분율 합이 1 보다 작다 — 버려진 칸이 있다는 증거
    assert narrow.coord_fraction.sum(axis=1).max() < 1.0


def test_window_mean_reproduces_a_hand_sliced_average():
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([_random_frame(100, min(Lx, Ly), s) for s in range(4)]
                      + [hexpts] * 4)
    t = np.arange(8, dtype=float)
    s = hex_order_series(frames, Lx=Lx, Ly=Ly, t_star=t)

    w = s.window_mean(4.0, 8.0)
    assert w["n_frames"] == 4
    assert w["psi6_global"] == pytest.approx(1.0, abs=1e-9)
    assert w["defect_fraction"] == pytest.approx(0.0, abs=1e-12)
    with pytest.raises(ValueError, match="프레임이 없다"):
        s.window_mean(100.0, 200.0)


# =============================================================================
# 시간창 g(r)
# =============================================================================
def test_rdf_windows_partition_all_frames_exactly_once():
    """창들이 프레임을 **정확히 한 번씩** 덮는다 — 빠짐도 중복도 없다."""
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(12)])
    t = np.linspace(0.0, 11.0, 12)
    w = rdf_windows(frames, Lx=L, Ly=L, t_star=t, n_windows=4, bins=60)

    assert w.g.shape == (4, 60)
    assert w.frames_per_window.sum() == 12
    assert np.all(w.frames_per_window == 3)
    assert np.all(np.diff(w.t_lo) > 0)


def test_rdf_windows_recovers_lattice_peak_at_nearest_neighbour_distance():
    """완벽 격자 창의 첫 봉이 최근접거리 `a = 1.0746 d` 에 온다 (빈 폭 내)."""
    hexpts, Lx, Ly = _perfect_hex()
    frames = np.stack([hexpts] * 4)
    t = np.arange(4, dtype=float)
    w = rdf_windows(frames, Lx=Lx, Ly=Ly, t_star=t, n_windows=2, bins=200)

    bin_width = float(w.r[1] - w.r[0])
    assert np.all(np.abs(w.first_peak_r - HEX_NN_OVER_D) <= 1.5 * bin_width)
    assert np.all(w.first_peak_g > 5.0)          # 격자는 봉이 날카롭다


def test_rdf_windows_rejects_more_windows_than_frames():
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(3)])
    with pytest.raises(ValueError, match="프레임이 없다"):
        rdf_windows(frames, Lx=L, Ly=L, t_star=np.arange(3.0), n_windows=9)


# =============================================================================
# Voronoi 프레임 (그림 재료)
# =============================================================================
def test_voronoi_frame_gives_one_polygon_per_particle_and_matches_series():
    """다각형 개수 = 입자 수, 배위수는 `hex_order_series` 와 일치한다."""
    hexpts, Lx, Ly = _perfect_hex()
    vf = voronoi_frame(hexpts, Lx=Lx, Ly=Ly)

    assert len(vf.polygons) == hexpts.shape[0]
    assert np.all(vf.coordination == 6)
    assert vf.psi6_abs == pytest.approx(np.ones(hexpts.shape[0]), abs=1e-9)
    # 완벽 격자의 Voronoi 셀은 정육각형 → 꼭짓점 6개
    assert all(p.shape[0] == 6 for p in vf.polygons)

    s = hex_order_series(hexpts[None], Lx=Lx, Ly=Ly)
    assert s.defect_fraction[0] == pytest.approx(
        float(np.mean(vf.coordination != 6)), abs=1e-12)


def test_voronoi_frame_cell_areas_sum_to_the_box_area():
    """셀 면적의 합 = 상자 면적. **다각형을 잘라내면 이 검사가 깨진다.**"""
    L = square_box_for(100)
    vf = voronoi_frame(_random_frame(100, L, 7), Lx=L, Ly=L)

    def shoelace(p: np.ndarray) -> float:
        x, y = p[:, 0], p[:, 1]
        return 0.5 * abs(float(np.dot(x, np.roll(y, -1))
                               - np.dot(y, np.roll(x, -1))))

    assert sum(shoelace(p) for p in vf.polygons) == pytest.approx(L * L,
                                                                 rel=1e-9)


def test_voronoi_frame_rejects_a_trajectory():
    L = square_box_for(100)
    frames = np.stack([_random_frame(100, L, s) for s in range(2)])
    with pytest.raises(ValueError, match="프레임 1개"):
        voronoi_frame(frames, Lx=L, Ly=L)


# =============================================================================
# 완화시간 적합 — 잡음을 완화로 착각하지 않는가
# =============================================================================
def test_fit_relaxation_recovers_a_known_tau():
    """합성 지수에서 `τ` 를 회수한다. 잡음 없는 극한이므로 정확해야 한다."""
    from simbot.analysis.structure import fit_relaxation

    t = np.linspace(0.0, 2.0, 200)
    tau_true, y0, y_inf = 0.17, 0.53, 0.29
    y = y_inf + (y0 - y_inf) * np.exp(-t / tau_true)
    fit = fit_relaxation(t, y)

    assert fit.converged
    assert fit.tau == pytest.approx(tau_true, rel=1e-6)
    assert fit.y0 == pytest.approx(y0, rel=1e-6)
    assert fit.y_inf == pytest.approx(y_inf, rel=1e-6)
    assert fit.r_squared > 0.9999
    assert fit.amplitude == pytest.approx(y0 - y_inf, rel=1e-6)


def test_fit_relaxation_refuses_to_call_stationary_noise_a_relaxation():
    """★ 정상상태 요동에 지수를 맞추면 항상 어떤 τ 가 나온다 — 거부해야 한다.

    이것이 이 함수의 존재 이유다. `A ≤ 10` 의 결함 분율은 첫 프레임부터 정상이고,
    그 신호에 τ 를 붙여 보고하면 없는 물리를 만들어낸다.
    """
    from simbot.analysis.structure import fit_relaxation

    rng = np.random.default_rng(11)
    t = np.linspace(0.0, 2.0, 200)
    noise = 0.06
    y = 0.295 + rng.normal(0.0, noise, t.size)          # 완화 없음
    fit = fit_relaxation(t, y, noise=noise)

    assert not fit.converged
    assert "요동" in fit.note or "오차가 값보다" in fit.note
    # noise 를 주지 않으면 조용히 통과한다 — 그래서 호출자가 반드시 줘야 한다
    naive = fit_relaxation(t, y)
    assert np.isfinite(naive.tau)


def test_fit_relaxation_accepts_a_real_relaxation_buried_in_noise():
    """완화 폭이 잡음의 4배면 통과한다 — 문턱이 과하게 보수적이지 않은지 확인."""
    from simbot.analysis.structure import fit_relaxation

    rng = np.random.default_rng(3)
    t = np.linspace(0.0, 2.0, 400)
    noise = 0.02
    y = 0.29 + 0.24 * np.exp(-t / 0.15) + rng.normal(0.0, noise, t.size)
    fit = fit_relaxation(t, y, noise=noise)

    assert fit.converged, fit.note
    assert fit.tau == pytest.approx(0.15, rel=0.25)


def test_bootstrap_over_seeds_recovers_the_known_spread():
    """★ 시드마다 τ 가 다른 합성 데이터에서 부트스트랩이 그 산포를 회수한다.

    `curve_fit` 공분산은 *평균 곡선 하나*의 적합 불확실성이므로 시드 간 변동을
    대표하지 않는다. 둘이 다른 것을 재고 있음을 이 테스트가 고정한다.
    """
    from simbot.analysis.structure import bootstrap_relaxation_over_seeds

    rng = np.random.default_rng(5)
    t = np.linspace(0.0, 2.0, 200)
    k, tau_sd = 32, 0.03
    taus = 0.10 + tau_sd * rng.standard_normal(k)       # 시드마다 다른 τ
    per_seed = np.array([0.29 + 0.20 * np.exp(-t / max(tt, 1e-3))
                         + rng.normal(0.0, 0.01, t.size) for tt in taus])

    out = bootstrap_relaxation_over_seeds(t, per_seed, n_resample=300, seed=1)
    assert out["n_converged"] > 250
    #  시드 평균의 SE = tau_sd/√k. 부트스트랩이 그 규모를 회수해야 한다
    expected_se = tau_sd / np.sqrt(k)
    assert out["tau_se_bootstrap"] == pytest.approx(expected_se, rel=0.6), out
    lo, hi = out["tau_ci95"]
    assert lo < out["tau"] < hi


def test_bootstrap_and_fit_se_are_different_quantities():
    """시드 간 변동이 없으면 부트스트랩 SE 가 적합 SE 보다 **작다** — 다른 양이다."""
    from simbot.analysis.structure import bootstrap_relaxation_over_seeds

    rng = np.random.default_rng(7)
    t = np.linspace(0.0, 2.0, 200)
    #  전 시드가 같은 τ, 잡음만 다르다 → 시드 앙상블 산포는 잡음에서만 온다
    per_seed = np.array([0.29 + 0.20 * np.exp(-t / 0.10)
                         + rng.normal(0.0, 0.01, t.size) for _ in range(32)])
    out = bootstrap_relaxation_over_seeds(t, per_seed, n_resample=300, seed=2)

    assert np.isfinite(out["tau_se_bootstrap"])
    assert np.isfinite(out["tau_se_fit"])
    #  적합 SE 는 평균 곡선(잡음 1/√32)에 대한 것이고 부트스트랩은 시드 재표집이다.
    #  같지 않아야 한다 — 같으면 한쪽이 다른 쪽을 그대로 베낀 것이다
    assert out["se_ratio_bootstrap_over_fit"] != pytest.approx(1.0, abs=1e-6)


def test_bootstrap_refuses_a_single_seed():
    from simbot.analysis.structure import bootstrap_relaxation_over_seeds

    t = np.linspace(0.0, 1.0, 50)
    with pytest.raises(ValueError, match="앙상블 오차를 낼 수 없다"):
        bootstrap_relaxation_over_seeds(t, np.zeros((1, 50)))
    with pytest.raises(ValueError, match=r"\(n_seeds, n_frames\)"):
        bootstrap_relaxation_over_seeds(t, np.zeros(50))


def test_seeds_for_target_sigma_scales_as_one_over_sqrt_k():
    """`SE ∝ 1/√k` 이므로 σ 를 2배로 올리려면 시드가 4배 필요하다."""
    from simbot.estimators import seeds_for_target_sigma

    r = seeds_for_target_sigma(diff=0.0109, se_diff=0.0104, k_current=16,
                               n_sigma=3.0)
    assert r["sigma_now"] == pytest.approx(0.0109 / 0.0104, rel=1e-9)
    #  t 보정은 정규보다 **더 많은** 시드를 요구한다
    assert r["k_needed"] > r["k_needed_normal"]
    assert r["k_needed_int"] >= r["k_needed"]

    #  2배 σ → 4배 k (정규 기준으로 확인, t 보정은 k 의존성이 있어 근사)
    r2 = seeds_for_target_sigma(diff=0.0109, se_diff=0.0104, k_current=16,
                                n_sigma=2.0, t_correction=False)
    r3 = seeds_for_target_sigma(diff=0.0109, se_diff=0.0104, k_current=16,
                                n_sigma=4.0, t_correction=False)
    assert r3["k_needed"] / r2["k_needed"] == pytest.approx(4.0, rel=1e-9)


def test_seeds_for_target_sigma_rejects_degenerate_input():
    from simbot.estimators import seeds_for_target_sigma

    with pytest.raises(ValueError, match="must be positive"):
        seeds_for_target_sigma(diff=0.01, se_diff=0.0, k_current=16)
    with pytest.raises(ValueError, match="at least 2"):
        seeds_for_target_sigma(diff=0.01, se_diff=0.01, k_current=1)


def test_fit_relaxation_rejects_too_few_points():
    from simbot.analysis.structure import fit_relaxation

    with pytest.raises(ValueError, match="3파라미터"):
        fit_relaxation(np.arange(3.0), np.arange(3.0))
    with pytest.raises(ValueError, match="길이가 다르다"):
        fit_relaxation(np.arange(5.0), np.arange(4.0))


# =============================================================================
# 커버리지 기하 — 해석적 왕복
# =============================================================================
def test_coverage_and_sigma_over_d_round_trip():
    for cov in (0.01, 0.0491, 0.0873, 0.10, 0.5):
        s = sigma_over_d_for_coverage(cov)
        assert coverage_from_sigma_over_d(s) == pytest.approx(cov, rel=1e-12)


def test_coverage_is_independent_of_n_because_density_star_is_one():
    """`φ` 는 `N` 에 의존하지 않는다 — `n* = 1` 이라 입자당 면적이 정확히 `d²` 다."""
    out = [box_si_for_coverage(n_particles=n, sigma_si=5e-6, coverage_max=0.10,
                               d_over_sigma_round=3.0)
           for n in (100, 400, 1024)]
    assert {round(o["coverage"], 12) for o in out} == {round(out[0]["coverage"], 12)}
    # 상자는 √N 으로 커진다
    assert out[1]["L_si"] / out[0]["L_si"] == pytest.approx(2.0, rel=1e-12)


def test_the_sweep_geometry_is_what_the_report_claims():
    """이 스윕의 기하: `σ = 5 µm` · `d/σ = 3` → `L = 150 µm` · `φ = 8.73 %`."""
    g = box_si_for_coverage(n_particles=100, sigma_si=5e-6, coverage_max=0.10,
                            d_over_sigma_round=3.0)
    assert g["d_si"] == pytest.approx(15e-6, rel=1e-12)
    assert g["L_si"] == pytest.approx(150e-6, rel=1e-12)
    assert g["L_star"] == pytest.approx(10.0, rel=1e-12)
    assert g["coverage"] == pytest.approx(np.pi / 36.0, rel=1e-12)
    assert g["coverage"] < 0.10


def test_finite_size_exponent_recovers_the_liquid_and_crystal_limits():
    """★ 합성 데이터: `ψ₆ ~ N^{-1/2}` 는 액체, `N⁰` 는 결정으로 읽어야 한다."""
    from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,
                                           psi6_finite_size_exponent)

    N = np.array([100.0, 256.0])
    liquid = psi6_finite_size_exponent(N, 0.3 * N ** -0.5)
    assert liquid.p == pytest.approx(0.5, rel=1e-12)
    assert liquid.eta6 == pytest.approx(2.0, rel=1e-12)
    assert liquid.reading == "liquid-like"

    crystal = psi6_finite_size_exponent(N, np.array([0.95, 0.95]))
    assert crystal.p == pytest.approx(0.0, abs=1e-12)
    assert crystal.reading == "hexatic-or-below"

    #  KTHNY 경계 η₆ = 1/4 → p = 1/16
    boundary = psi6_finite_size_exponent(N, 0.3 * N ** -(1.0 / 16.0))
    assert boundary.eta6 == pytest.approx(KTHNY_ETA6_HEXATIC_LIQUID, rel=1e-9)


def test_finite_size_exponent_propagates_the_error_bars():
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([100.0, 256.0])
    y = np.array([0.0477, 0.0298])
    se = np.array([0.0019, 0.0015])
    fit = psi6_finite_size_exponent(N, y, se)

    assert np.isfinite(fit.p_se) and fit.p_se > 0
    #  상대오차가 로그차로 나눠진다 — 손계산과 대조
    expected = float(np.hypot(se[0] / y[0], se[1] / y[1])
                     / abs(np.log(N[1] / N[0])))
    assert fit.p_se == pytest.approx(expected, rel=1e-12)
    assert fit.eta6_se == pytest.approx(4.0 * expected, rel=1e-12)
    #  SE 없이 부르면 nan 이고 '모른다'가 드러나야 한다
    assert not np.isfinite(psi6_finite_size_exponent(N, y).p_se)


def test_finite_size_exponent_three_points_is_a_weighted_fit():
    """셋 이상이면 형태를 논할 수 있다 — `n_points` 가 그 사실을 나른다."""
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([64.0, 144.0, 400.0])
    fit = psi6_finite_size_exponent(N, 0.4 * N ** -0.5,
                                    0.01 * np.ones(3) * 0.4 * N ** -0.5)
    assert fit.n_points == 3
    assert fit.p == pytest.approx(0.5, rel=1e-9)


def test_finite_size_exponent_form_test_needs_three_points():
    """★ 두 점은 직선을 **가정**한다 — 형태를 검증할 수 없다는 사실이 실려 나가야 한다."""
    from simbot.analysis.structure import psi6_finite_size_exponent

    N2 = np.array([100.0, 256.0])
    two = psi6_finite_size_exponent(N2, 0.3 * N2 ** -0.5, 0.001 * np.ones(2))
    assert two.n_points == 2
    assert not two.form_is_testable
    assert not np.isfinite(two.chi2_reduced)
    assert "검증 불가" in two.form_verdict()

    N4 = np.array([64.0, 144.0, 256.0, 400.0])
    y = 0.3 * N4 ** -0.5
    four = psi6_finite_size_exponent(N4, y, 0.01 * y)
    assert four.form_is_testable
    assert four.chi2_reduced == pytest.approx(0.0, abs=1e-16)
    assert "모순 없음" in four.form_verdict()
    assert four.p == pytest.approx(0.5, rel=1e-9)
    assert four.amplitude == pytest.approx(0.3, rel=1e-9)


def test_finite_size_exponent_flags_a_curve_that_is_not_a_power_law():
    """★ 멱함수가 **아닌** 데이터에 단일 지수를 붙이면 경고해야 한다.

    이것이 형태 검증의 존재 이유다 — 지수는 항상 나오지만 형태가 틀렸으면
    그 지수로 상을 판정하면 안 된다.
    """
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([64.0, 144.0, 256.0, 400.0])
    #  로그-로그에서 휘는 곡선 (지수가 N 과 함께 변한다)
    y = 0.3 * N ** -0.5 * (1.0 + 0.35 * np.log(N / 64.0))
    fit = psi6_finite_size_exponent(N, y, 0.005 * y)

    assert fit.form_is_testable
    assert fit.chi2_reduced > 3.0, fit.chi2_reduced
    assert "벗어난다" in fit.form_verdict()
    #  잔차가 부호를 바꾼다 — 휘어 있다는 서명
    r = np.array(fit.residuals)
    assert np.any(r > 0) and np.any(r < 0)


def test_finite_size_exponent_chi2_needs_error_bars():
    """오차막대 없이는 `χ²` 를 계산할 수 없다 — 조용히 0 을 내지 않는다."""
    from simbot.analysis.structure import psi6_finite_size_exponent

    N = np.array([64.0, 144.0, 256.0, 400.0])
    fit = psi6_finite_size_exponent(N, 0.3 * N ** -0.5)
    assert fit.n_points == 4
    assert not np.isfinite(fit.chi2_reduced)
    assert not fit.form_is_testable


def test_phase_reading_does_not_call_a_crystal_hexatic():
    """★★ `η₆ ≤ 1/4` 는 hexatic 의 **필요조건일 뿐**이다 — 결정도 만족한다.

    2026-07-29 에 이 구별 없이 `A ≥ 13.3` 의 결정을 "hexatic 가능? 예" 로 출력했다.
    갈리는 것은 `ψ₆` 의 **크기**다.
    """
    from simbot.analysis.structure import (phase_from_finite_size,
                                           psi6_finite_size_exponent)

    N = np.array([144.0, 256.0, 400.0])

    #  ① 결정: psi6 가 O(1) 로 포화 → p ≈ 0 → eta6 ≈ 0 (hexatic 과 같은 지수!)
    y_cry = np.array([0.6720, 0.6841, 0.6881])
    f_cry = psi6_finite_size_exponent(N, y_cry, np.array([0.0032, 0.0074, 0.0080]))
    r_cry = phase_from_finite_size(f_cry, y_cry[-1])
    assert abs(f_cry.eta6) < 0.25                    # 지수만 보면 hexatic 후보다
    assert r_cry["exponent_alone_would_say"] == "hexatic-or-below"
    assert r_cry["phase"] == "crystal", r_cry        # ★ 크기가 구별한다
    assert "포화" in r_cry["why"]

    #  ② 액체: psi6 작고 N^-1/2 로 줄어든다
    y_liq = 0.5 * N ** -0.5
    f_liq = psi6_finite_size_exponent(N, y_liq, 0.02 * y_liq)
    r_liq = phase_from_finite_size(f_liq, y_liq[-1])
    assert r_liq["phase"] == "isotropic-liquid"
    assert r_liq["exponent_alone_would_say"] == "not-hexatic"

    #  ③ hexatic 후보: 지수가 상한 **아래**이고 psi6 는 결정 바닥 아래.
    #     ★ η₆ 를 상한(0.25)에 정확히 놓으면 `η₆ + 3σ > 0.25` 라 inconclusive 가
    #       맞다 — 경계에 앉은 값을 hexatic 이라 부르면 안 된다. 그래서 0.1 을 쓴다.
    y_hex = 0.30 * (N / N[0]) ** -(0.1 / 4.0)        # eta6 = 0.1 < 0.25
    f_hex = psi6_finite_size_exponent(N, y_hex, 0.005 * y_hex)
    r_hex = phase_from_finite_size(f_hex, y_hex[-1])
    assert y_hex[-1] < 0.5
    assert f_hex.eta6 == pytest.approx(0.1, rel=1e-6)
    assert r_hex["phase"] == "hexatic-candidate", r_hex


def test_phase_reading_flags_inconclusive_when_eta6_straddles_the_ceiling():
    from simbot.analysis.structure import (phase_from_finite_size,
                                           psi6_finite_size_exponent)

    N = np.array([144.0, 256.0, 400.0])
    #  eta6 ≈ 0.25 이지만 오차가 커서 상한을 걸친다
    y = 0.30 * (N / N[0]) ** -(1.0 / 16.0)
    fit = psi6_finite_size_exponent(N, y, 0.06 * y)
    r = phase_from_finite_size(fit, y[-1])
    assert r["phase"] in ("inconclusive", "hexatic-candidate")
    if r["phase"] == "inconclusive":
        assert "걸친다" in r["why"]


def test_finite_size_exponent_rejects_degenerate_input():
    from simbot.analysis.structure import psi6_finite_size_exponent

    with pytest.raises(ValueError, match="최소 2개"):
        psi6_finite_size_exponent(np.array([100.0]), np.array([0.05]))
    with pytest.raises(ValueError, match="0 이하"):
        psi6_finite_size_exponent(np.array([100.0, 256.0]),
                                  np.array([0.05, 0.0]))


def test_coverage_does_not_touch_the_reduced_config_at_all():
    """★★ 커버리지를 바꿔도 축약 단위 런은 **비트 단위로 같다.**

    그래서 "커버리지 3.7 % 대조군" 을 **시뮬레이션으로 돌리는 것은 산술 항등식**이다
    (`conftest` 규칙 ③: 요동하지 않는 측정값은 항등식). 컴퓨트를 쓸 항목이 아니라
    이 테스트로 고정할 항목이다 — 2026-07-29 에 그 제안을 했다가 여기서 기각했다.

    커버리지가 바꾸는 것은 SI 환산뿐이다: `τ_d ∝ d²` (σ 고정).
    """
    from dataclasses import asdict

    from simbot.run import Soft2DRunConfig
    from simbot.units import scales_soft2d

    made = {}
    for d_over_sigma in (3.0, 4.6):                 # 커버리지 8.73 % vs 3.71 %
        g = box_si_for_coverage(n_particles=100, sigma_si=5e-6,
                                coverage_max=0.10,
                                d_over_sigma_round=d_over_sigma)
        cfg = Soft2DRunConfig(amplitude=10.0, n_particles=100, init="random",
                              box_shape="square", r_min=0.05, dt_star=5.01e-5,
                              equil_tau=0.0, prod_tau=80.0, seed=5)
        made[d_over_sigma] = (asdict(cfg), cfg.hash(), g,
                              scales_soft2d(d_si=g["d_si"], sigma_si=5e-6,
                                            T_si=298.15))
    a, b = made[3.0], made[4.6]

    assert a[2]["coverage"] != b[2]["coverage"]      # 커버리지는 실제로 달라야 한다
    assert b[2]["coverage"] < 0.0371 + 1e-4
    #  ★ 그런데 축약 config 와 run_hash 는 동일하다
    assert a[0] == b[0]
    assert a[1] == b[1], "run_hash 가 다르면 커버리지가 동역학에 새어 들어간 것이다"
    #  SI 환산만 바뀐다: tau_d ∝ d²
    assert b[3].time_si / a[3].time_si == pytest.approx(
        (b[2]["d_si"] / a[2]["d_si"]) ** 2, rel=1e-12)
    assert b[3].time_si > a[3].time_si


def test_box_for_coverage_refuses_a_ratio_that_breaks_the_ceiling():
    """`d/σ` 를 직접 주더라도 커버리지 상한을 넘으면 **거부한다.**"""
    with pytest.raises(ValueError, match="커버리지"):
        box_si_for_coverage(n_particles=100, sigma_si=5e-6, coverage_max=0.10,
                            d_over_sigma_round=2.0)
    with pytest.raises(ValueError, match=r"\(0, 1\]"):
        sigma_over_d_for_coverage(0.0)


def test_scales_soft2d_uses_the_lattice_spacing_for_length_and_sigma_for_drag():
    """★ 길이 척도는 `d`, 항력은 `σ` 다. 섞으면 `τ_d` 가 `(d/σ)²` 배 틀린다."""
    from simbot.units import (kT_si, scales_soft2d, stokes_drag_si,
                              stokes_einstein_D_si, water_viscosity_si)

    sigma, T = 5e-6, 298.15
    d = 3.0 * sigma
    sc = scales_soft2d(d_si=d, sigma_si=sigma, T_si=T)

    eta, extrapolated = water_viscosity_si(T)
    assert not extrapolated
    D0 = stokes_einstein_D_si(T, stokes_drag_si(eta, sigma / 2.0))
    assert sc.length_si == pytest.approx(d, rel=1e-12)
    assert sc.energy_si == pytest.approx(kT_si(T), rel=1e-12)
    assert sc.time_si == pytest.approx(d**2 / D0, rel=1e-12)

    # σ 를 길이 척도로 잘못 쓰면 (d/σ)² = 9 배 틀린다 — 그 크기를 명시한다
    wrong = scales_soft2d(d_si=sigma, sigma_si=sigma, T_si=T)
    assert sc.time_si / wrong.time_si == pytest.approx(9.0, rel=1e-12)


def test_scales_soft2d_refuses_to_extrapolate_the_viscosity_table():
    """표(293–308 K) 밖에서는 조용히 외삽하지 않고 `gamma_si` 를 요구한다."""
    from simbot.units import scales_soft2d

    with pytest.raises(ValueError, match="gamma_si"):
        scales_soft2d(d_si=15e-6, sigma_si=5e-6, T_si=350.0)
    # 명시하면 표 밖에서도 통과한다
    sc = scales_soft2d(d_si=15e-6, sigma_si=5e-6, T_si=350.0, gamma_si=1e-8)
    assert sc.time_si > 0.0
