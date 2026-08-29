"""S6 — 그림.

## 이 파일이 지키는 것

그림의 **모양**은 테스트하지 않는다 (픽셀 비교는 깨지기 쉽고 무엇이 틀렸는지도
알려주지 않는다). 대신 §S6 게이트를 지킨다:

1. **캡션 없는 그림은 만들 수 없다** — 사후 검사가 아니라 생성 시점 차단
2. **안 그린 그림에는 이유가 남는다** — 조용한 누락 금지
3. **모든 시간·길이 축에 이중 표기** (무차원 + SI)
4. 그림 안 텍스트는 **영문** — matplotlib 기본 폰트에 한글 글리프가 없다
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from simbot import viz


# =============================================================================
# 합성 데이터 — HOOMD 없이 도는 테스트용
# =============================================================================
@pytest.fixture
def fake_run():
    """조화 트랩 정상 상태를 흉내낸 궤적. AR(1) 로 만든다 (실제 스킴과 같은 구조)."""
    rng = np.random.default_rng(7)
    n_frames, n_part, dim, dt = 60, 200, 2, 5e-3
    stride = 16
    a = (1.0 - dt) ** stride
    traj = np.empty((n_frames, n_part, dim), dtype=np.float32)
    x = rng.normal(size=(n_part, dim))
    for f in range(n_frames):
        x = a * x + np.sqrt(1 - a**2) * rng.normal(size=(n_part, dim))
        traj[f] = x
    lags = np.unique(np.geomspace(1, n_frames - 1, 20).astype(int)) * stride
    lags_tau = lags * dt
    msd = 2 * dim * (1 - np.exp(-lags_tau))
    var = 1.0 + 0.01 * rng.normal(size=20)
    return {"traj": traj, "lags_steps": lags, "msd": msd,
            "indep_var": var, "indep_kT": var.copy()}   # 대수적으로 동일


@pytest.fixture
def runs(fake_run):
    return {f"s{i}": fake_run for i in range(1, 5)}


TRAP_SI = dict(tau_trap_si=8.0644e-3, l_trap_si=2.0352e-8)


# =============================================================================
# 캡션 강제 — 이 절이 이 파일의 존재 이유다
# =============================================================================
def test_save_refuses_empty_caption(tmp_path):
    """★ §S6 게이트: 캡션 없는 그림은 산출물로 인정하지 않는다."""
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    with pytest.raises(ValueError, match="캡션이 없다"):
        fs.save(fig, "x", caption="   ", shows="무언가")
    plt.close(fig)


def test_save_refuses_empty_shows(tmp_path):
    """'무엇을 보이려는 그림인가'에 답하지 못하면 진단에 쓸 수 없다."""
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    with pytest.raises(ValueError, match="shows"):
        fs.save(fig, "x", caption="설명", shows="")
    plt.close(fig)


def test_save_writes_png_and_records(tmp_path):
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    rec = fs.save(fig, "01_test", caption="캡션", shows="무엇을")
    assert rec.path.exists() and rec.name == "01_test.png"
    assert fs.captions == {"01_test.png": "캡션"}


def test_skip_refuses_empty_reason(tmp_path):
    """★ 빠진 그림이 '해당 없음'인지 '잊음'인지 구별되어야 한다."""
    fs = viz.FigureSet(tmp_path)
    with pytest.raises(ValueError, match="건너뛴 이유가 없다"):
        fs.skip("07_rdf.png", "")


def test_skip_records_reason(tmp_path):
    fs = viz.FigureSet(tmp_path)
    fs.skip("07_rdf.png", "쌍 상호작용 없음")
    assert fs.skipped == {"07_rdf.png": "쌍 상호작용 없음"}


# =============================================================================
# 06_figures.md
# =============================================================================
def test_figures_md_lists_captions_and_skips(tmp_path):
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    fs.save(fig, "01_a", caption="첫 그림 캡션", shows="첫 그림이 보이려는 것")
    fs.skip("02_b.png", "데이터가 없다")
    md = fs.figures_md()
    assert "# S6 FIGURES" in md
    assert "첫 그림 캡션" in md
    assert "첫 그림이 보이려는 것" in md
    assert "그리지 않은 그림과 이유" in md
    assert "데이터가 없다" in md
    assert "![첫 그림 캡션](figs/01_a.png)" in md


def test_figures_md_without_skips_omits_that_section(tmp_path):
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    fs.save(fig, "01_a", caption="c", shows="s")
    assert "그리지 않은 그림" not in fs.figures_md()


def test_figures_md_states_the_font_constraint(tmp_path):
    fs = viz.FigureSet(tmp_path)
    fs.skip("x.png", "이유")
    assert "한글 글리프" in fs.figures_md()


def test_multiline_caption_does_not_break_image_syntax(tmp_path):
    """★ 여러 줄 캡션을 `![...]` 안에 넣으면 마크다운 이미지가 깨진다."""
    import re
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    fs.save(fig, "01_a", caption="첫 줄.\n\n둘째 단락.", shows="s")
    md = fs.figures_md()
    imgs = re.findall(r"!\[([^\]]*)\]\(([^)]+)\)", md)
    assert len(imgs) == 1
    assert "\n" not in imgs[0][0]
    assert "둘째 단락." in md          # 전체 캡션은 본문으로 나간다


def test_alt_text_is_truncated_not_dropped():
    rec = viz.FigureRecord(name="x.png", caption="가" * 200, shows="s",
                           path=Path("x.png"))
    alt = viz.alt_text(rec)
    assert len(alt) <= 90 and alt.endswith("…")


def test_alt_text_keeps_short_caption_intact():
    rec = viz.FigureRecord(name="x.png", caption="짧은 캡션", shows="s",
                           path=Path("x.png"))
    assert viz.alt_text(rec) == "짧은 캡션"


# =============================================================================
# 이중 축
# =============================================================================
def test_add_si_axis_maps_reduced_to_si():
    """무차원 10 τ_trap → 80.644 ms. (보조축 한계는 draw 후에 확정된다)"""
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_xlim(0, 10)
    sec = viz.add_si_axis(ax, 8.0644e-3, "ms", axis="x", si_multiplier=1e3)
    fig.canvas.draw()
    assert sec.get_xlim()[1] == pytest.approx(80.644, rel=1e-6)
    plt.close(fig)


def test_add_si_axis_supports_y():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    ax.set_ylim(0, 4)
    sec = viz.add_si_axis(ax, 100.0, "nm²", axis="y")
    fig.canvas.draw()
    assert sec.get_ylim()[1] == pytest.approx(400.0, rel=1e-6)
    plt.close(fig)


def test_si_axis_label_carries_the_unit():
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    sec = viz.add_si_axis(ax, 1.0, "ms", axis="x", label="$t$  [ms]")
    assert "ms" in sec.get_xlabel()
    plt.close(fig)


# =============================================================================
# 개별 그림 — 캡션에 수치가 들어가는가
# =============================================================================
def test_msd_figure_reports_plateau_in_caption(tmp_path, runs):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_msd(fs, runs, dim=2, dt_star=5e-3, **TRAP_SI)
    assert "2d = 4" in rec.caption or "`2d = 4`" in rec.caption
    assert "자유확산" in rec.shows          # 단시간 극한 검증이 목적에 있어야 한다


def test_distribution_figure_warns_against_exactly_three(tmp_path, fake_run):
    """★ 첨도가 정확히 3.000 이면 오히려 이상하다 — 캡션이 이걸 말해야 한다."""
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_position_distribution(fs, fake_run["traj"], dt_star=5e-3, dim=2,
                                         l_trap_si=TRAP_SI["l_trap_si"],
                                         frame_interval_steps=16)
    assert "3.000" in rec.caption and "이상하다" in rec.caption
    assert "균일분포" in rec.caption


def test_distribution_figure_normalises_so_it_tests_shape_not_width(tmp_path, fake_run):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_position_distribution(fs, fake_run["traj"], dt_star=5e-3, dim=2,
                                         l_trap_si=TRAP_SI["l_trap_si"],
                                         frame_interval_steps=16)
    assert "형태" in rec.caption and "폭" in rec.caption


def test_distribution_uses_independent_frames_only(tmp_path, fake_run):
    """★ 상관 프레임을 전부 쓰면 표본 수가 부풀어 오차막대가 거짓으로 작아진다.

    같은 함정을 KS 검정에서 겪었다 ([[ks-test-needs-independent-samples]]).
    `frame_interval_steps=16`, `dt*=5e-3` → 프레임 간격 `0.08 τ` →
    `2 τ` 독립성에는 25 프레임마다 하나. 60 프레임 중 3개만 남는다.
    """
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_position_distribution(fs, fake_run["traj"], dt_star=5e-3, dim=2,
                                         l_trap_si=TRAP_SI["l_trap_si"],
                                         frame_interval_steps=16)
    n_total = fake_run["traj"].shape[0]
    n_used = int(__import__("re").search(r"독립 프레임 (\d+)개", rec.caption).group(1))
    assert n_used < n_total, "전체 프레임을 썼다 — 상관 표본이다"
    assert "거짓으로 작아진다" in rec.caption


def test_distribution_error_bar_shrinks_with_more_independent_frames(tmp_path,
                                                                    fake_run):
    """오차막대가 실제로 표본 수를 반영하는지 — 상수를 박아두지 않았는지 확인."""
    import re
    fs = viz.FigureSet(tmp_path)
    a = viz.plot_position_distribution(
        fs, fake_run["traj"], dt_star=5e-3, dim=2, name="a",
        l_trap_si=TRAP_SI["l_trap_si"], frame_interval_steps=16)
    b = viz.plot_position_distribution(
        fs, fake_run["traj"], dt_star=5e-3, dim=2, name="b",
        l_trap_si=TRAP_SI["l_trap_si"], frame_interval_steps=16,
        decorrelation_tau=0.1)                      # 더 많은 프레임을 독립으로 취급
    se = lambda r: float(re.search(r"± ([\d.]+)`", r.caption).group(1))
    assert se(b) < se(a)


def test_stationarity_figure_demonstrates_the_identity(tmp_path, runs):
    """★ kT_conf 가 ⟨x²⟩ 와 대수적으로 같다는 것을 **보여준다** (주장하지 않는다)."""
    import re
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_stationarity(fs, runs, dt_star=5e-3, sample_interval_steps=400,
                                tau_trap_si=TRAP_SI["tau_trap_si"])
    assert "대수적으로 동일" in rec.caption
    assert "독립 검사가 되려면 쌍 상호작용이 필요" in rec.caption
    # 잔차 최대값이 **숫자로** 캡션에 들어가야 한다 (주장이 아니라 측정값)
    assert re.search(r"최대 `\d\.\de[+-]\d+`", rec.caption), rec.caption


def test_stationarity_residual_is_machine_precision_when_identical(tmp_path, runs):
    """합성 데이터에서 두 통계량을 같게 만들었으므로 잔차가 0 이어야 한다.

    실제 런에서는 4.4e-16 (부동소수점 수준) 이 나온다 — 어느 쪽이든
    "독립 검사가 아니다"라는 결론이 그림에서 읽힌다.
    """
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_stationarity(fs, runs, dt_star=5e-3, sample_interval_steps=400,
                                tau_trap_si=TRAP_SI["tau_trap_si"])
    assert "0.0e+00" in rec.caption or "e-1" in rec.caption


def test_stationarity_shows_mentions_independence_check(tmp_path, runs):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_stationarity(fs, runs, dt_star=5e-3, sample_interval_steps=400,
                                tau_trap_si=TRAP_SI["tau_trap_si"])
    assert "독립이 아니다" in rec.shows


def test_displacement_figure_disclaims_per_step(tmp_path, fake_run):
    """★ 프레임 변위는 스텝당 변위가 아니다 — 균일 노이즈 √3 상한이 적용되지 않는다."""
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_frame_displacements(fs, fake_run["traj"], dt_star=5e-3,
                                       frame_interval_steps=16, sigma_star=491.358)
    assert "스텝당 변위가 아니다" in rec.caption
    assert "√3" in rec.caption


def test_equipartition_figure_plots_the_competing_hypothesis(tmp_path):
    """경쟁 가설(exact)을 함께 그려야 판별 가능성이 보인다."""
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_equipartition_vs_dt(fs, {5e-3: (1.00577, 0.00446),
                                            2.5e-3: (1.00169, 0.00519)})
    assert "검정력" in rec.caption
    assert "exact" in rec.shows or "구별" in rec.shows


def test_snapshots_title_carries_both_units(tmp_path, fake_run):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_snapshots(fs, fake_run["traj"], dim=2, dt_star=5e-3,
                             frame_interval_steps=16,
                             tau_trap_si=TRAP_SI["tau_trap_si"],
                             l_trap_si=TRAP_SI["l_trap_si"])
    assert "nm" in rec.caption                      # ℓ_trap 물리값
    assert "독립 복제" in rec.caption               # 서스펜션이 아니라는 경고


# =============================================================================
# 오케스트레이터
# =============================================================================
def test_trap_diagnostics_produces_the_required_set(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, **TRAP_SI)
    names = {r.name for r in fs.records}
    for expect in ("01_msd.png", "03_distribution.png", "04_stationarity.png",
                   "05_displacements.png", "06_snapshots.png"):
        assert expect in names, expect
    assert all(p.exists() for p in tmp_path.glob("*.png"))


def test_every_figure_has_nonempty_caption_and_shows(tmp_path, runs):
    """게이트 전수 검사 — 하나라도 비면 §S6 위반."""
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, **TRAP_SI)
    for r in fs.records:
        assert r.caption.strip() and r.shows.strip(), r.name


def test_single_dt_skips_bias_figure_with_reason(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, by_dt={5e-3: (1.0, 0.004)},
        **TRAP_SI)
    assert "02_equipartition_dt.png" in fs.skipped
    assert "converge" in fs.skipped["02_equipartition_dt.png"]


def test_two_dt_draws_bias_figure(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358,
        by_dt={5e-3: (1.00577, 0.00446), 2.5e-3: (1.00169, 0.00519)}, **TRAP_SI)
    assert "02_equipartition_dt.png" in {r.name for r in fs.records}


def test_no_pair_skips_rdf_because_it_is_undefined(tmp_path, runs):
    """★ 쌍 상호작용이 없으면 g(r) 이 **정의되지 않는다** — 미구현과 다른 이유다."""
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, has_pair=False, **TRAP_SI)
    assert "정의되지 않는다" in fs.skipped["07_rdf.png"]


def test_with_pair_skips_rdf_because_it_is_unimplemented(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, has_pair=True, **TRAP_SI)
    assert "아직 없다" in fs.skipped["07_rdf.png"]


def test_3d_skips_snapshots_with_reason(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=3, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, **TRAP_SI)
    assert "fresnel" in fs.skipped["06_snapshots.png"]


def test_matplotlib_uses_headless_backend():
    """Agg 가 아니면 CI·헤드리스에서 창을 띄우려다 죽는다."""
    import matplotlib
    assert matplotlib.get_backend().lower() == "agg"


def test_figure_text_is_ascii_only(tmp_path, runs):
    """★ 그림 안 텍스트가 한글이면 두부(□)로 렌더된다 — 소스에서 막는다."""
    src = Path(viz.__file__).read_text(encoding="utf-8")
    # set_xlabel/set_ylabel/set_title/label= 안의 문자열에 한글이 있으면 안 된다
    import re
    bad = []
    for m in re.finditer(r"(set_xlabel|set_ylabel|set_title|suptitle)\(([^\n]*)", src):
        if re.search(r"[가-힣]", m.group(2)):
            bad.append(m.group(0)[:70])
    assert not bad, f"그림 안 텍스트에 한글: {bad}"
