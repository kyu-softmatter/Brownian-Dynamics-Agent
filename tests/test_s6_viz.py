"""S6 — figures.

## What this file enforces

It does not test what a figure **looks like** (a pixel comparison is brittle and
tells you nothing about what went wrong). It enforces the §S6 gates instead:

1. **A figure without a caption cannot be created** — blocked at creation time, not
   checked after the fact
2. **A figure not drawn leaves a reason** — no silent omissions
3. **Every time and length axis is dual-labelled** (dimensionless + SI)
4. Text inside a figure is **English** — matplotlib's default font has no Hangul
   glyphs
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest

from simbot import viz


# =============================================================================
# Synthetic data — so the tests run without HOOMD
# =============================================================================
@pytest.fixture
def fake_run():
    """A trajectory mimicking a harmonic trap's steady state. Built as AR(1), the
    same structure as the real scheme."""
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
            "indep_var": var, "indep_kT": var.copy()}   # algebraically identical


@pytest.fixture
def runs(fake_run):
    return {f"s{i}": fake_run for i in range(1, 5)}


TRAP_SI = dict(tau_trap_si=8.0644e-3, l_trap_si=2.0352e-8)


# =============================================================================
# Captions are compulsory — this section is why the file exists
# =============================================================================
def test_save_refuses_empty_caption(tmp_path):
    """★ The §S6 gate: an uncaptioned figure is not accepted as an artefact."""
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    with pytest.raises(ValueError, match="no caption"):
        fs.save(fig, "x", caption="   ", shows="something")
    plt.close(fig)


def test_save_refuses_empty_shows(tmp_path):
    """A figure that cannot answer 'what is this meant to show' is useless for
    diagnosis."""
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    with pytest.raises(ValueError, match="shows"):
        fs.save(fig, "x", caption="a description", shows="")
    plt.close(fig)


def test_save_writes_png_and_records(tmp_path):
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    rec = fs.save(fig, "01_test", caption="a caption", shows="what for")
    assert rec.path.exists() and rec.name == "01_test.png"
    assert fs.captions == {"01_test.png": "a caption"}


def test_skip_refuses_empty_reason(tmp_path):
    """★ A missing figure has to be distinguishable between 'not applicable' and
    'forgotten'."""
    fs = viz.FigureSet(tmp_path)
    with pytest.raises(ValueError, match="no reason for skipping"):
        fs.skip("07_rdf.png", "")


def test_skip_records_reason(tmp_path):
    fs = viz.FigureSet(tmp_path)
    fs.skip("07_rdf.png", "no pair interaction")
    assert fs.skipped == {"07_rdf.png": "no pair interaction"}


# =============================================================================
# 06_figures.md
# =============================================================================
def test_figures_md_lists_captions_and_skips(tmp_path):
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    fs.save(fig, "01_a", caption="the first figure's caption",
            shows="what the first figure shows")
    fs.skip("02_b.png", "there is no data")
    md = fs.figures_md()
    assert "# S6 FIGURES" in md
    assert "the first figure's caption" in md
    assert "what the first figure shows" in md
    assert "Figures not drawn, and why" in md
    assert "there is no data" in md
    assert "![the first figure's caption](figs/01_a.png)" in md


def test_figures_md_without_skips_omits_that_section(tmp_path):
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    fs.save(fig, "01_a", caption="c", shows="s")
    assert "Figures not drawn" not in fs.figures_md()


def test_figures_md_states_the_font_constraint(tmp_path):
    fs = viz.FigureSet(tmp_path)
    fs.skip("x.png", "a reason")
    assert "Hangul glyphs" in fs.figures_md()


def test_multiline_caption_does_not_break_image_syntax(tmp_path):
    """★ A multi-line caption inside `![...]` breaks the markdown image."""
    import re
    import matplotlib.pyplot as plt
    fs = viz.FigureSet(tmp_path)
    fig, _ = plt.subplots()
    fs.save(fig, "01_a", caption="First line.\n\nSecond paragraph.", shows="s")
    md = fs.figures_md()
    imgs = re.findall(r"!\[([^\]]*)\]\(([^)]+)\)", md)
    assert len(imgs) == 1
    assert "\n" not in imgs[0][0]
    assert "Second paragraph." in md   # the full caption goes out as body text


def test_alt_text_is_truncated_not_dropped():
    rec = viz.FigureRecord(name="x.png", caption="a" * 200, shows="s",
                           path=Path("x.png"))
    alt = viz.alt_text(rec)
    assert len(alt) <= 90 and alt.endswith("…")


def test_alt_text_keeps_short_caption_intact():
    rec = viz.FigureRecord(name="x.png", caption="a short caption", shows="s",
                           path=Path("x.png"))
    assert viz.alt_text(rec) == "a short caption"


# =============================================================================
# Dual axes
# =============================================================================
def test_add_si_axis_maps_reduced_to_si():
    """Dimensionless 10 τ_trap → 80.644 ms. (The secondary axis limits are only
    fixed after a draw.)"""
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
# Individual figures — does the caption carry the numbers
# =============================================================================
def test_msd_figure_reports_plateau_in_caption(tmp_path, runs):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_msd(fs, runs, dim=2, dt_star=5e-3, **TRAP_SI)
    assert "2d = 4" in rec.caption or "`2d = 4`" in rec.caption
    assert "free diffusion" in rec.shows   # verifying the short-time limit has to be
                                           # among the purposes


def test_distribution_figure_warns_against_exactly_three(tmp_path, fake_run):
    """★ A kurtosis of exactly 3.000 would be suspicious — the caption has to say
    so."""
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_position_distribution(fs, fake_run["traj"], dt_star=5e-3, dim=2,
                                         l_trap_si=TRAP_SI["l_trap_si"],
                                         frame_interval_steps=16)
    assert "3.000" in rec.caption and "suspicious" in rec.caption
    assert "uniform" in rec.caption


def test_distribution_figure_normalises_so_it_tests_shape_not_width(tmp_path, fake_run):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_position_distribution(fs, fake_run["traj"], dt_star=5e-3, dim=2,
                                         l_trap_si=TRAP_SI["l_trap_si"],
                                         frame_interval_steps=16)
    assert "shape" in rec.caption and "width" in rec.caption


def test_distribution_uses_independent_frames_only(tmp_path, fake_run):
    """★ Using every correlated frame inflates the sample count and makes the error
    bar falsely small.

    The same trap was hit with the KS test ([[ks-test-needs-independent-samples]]).
    `frame_interval_steps=16`, `dt*=5e-3` → a frame interval of `0.08 τ` → one frame
    in 25 for `2 τ` independence. Only 3 of 60 frames survive.
    """
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_position_distribution(fs, fake_run["traj"], dt_star=5e-3, dim=2,
                                         l_trap_si=TRAP_SI["l_trap_si"],
                                         frame_interval_steps=16)
    n_total = fake_run["traj"].shape[0]
    n_used = int(__import__("re").search(r"\*\*(\d+) independent frames",
                                         rec.caption).group(1))
    assert n_used < n_total, "it used every frame — those are correlated samples"
    assert "falsely small" in rec.caption


def test_distribution_error_bar_shrinks_with_more_independent_frames(tmp_path,
                                                                    fake_run):
    """Does the error bar actually reflect the sample count — checking no constant was
    nailed in."""
    import re
    fs = viz.FigureSet(tmp_path)
    a = viz.plot_position_distribution(
        fs, fake_run["traj"], dt_star=5e-3, dim=2, name="a",
        l_trap_si=TRAP_SI["l_trap_si"], frame_interval_steps=16)
    b = viz.plot_position_distribution(
        fs, fake_run["traj"], dt_star=5e-3, dim=2, name="b",
        l_trap_si=TRAP_SI["l_trap_si"], frame_interval_steps=16,
        decorrelation_tau=0.1)                      # treat more frames as independent
    se = lambda r: float(re.search(r"± ([\d.]+)`", r.caption).group(1))
    assert se(b) < se(a)


def test_stationarity_figure_demonstrates_the_identity(tmp_path, runs):
    """★ It **shows** that kT_conf is algebraically identical to ⟨x²⟩; it does not
    assert it."""
    import re
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_stationarity(fs, runs, dt_star=5e-3, sample_interval_steps=400,
                                tau_trap_si=TRAP_SI["tau_trap_si"])
    assert "algebraically identical" in rec.caption
    assert "requires a pair interaction" in rec.caption
    # the maximum residual has to appear in the caption **as a number** (a
    # measurement, not an assertion)
    assert re.search(r"at most `\d\.\de[+-]\d+`", rec.caption), rec.caption


def test_stationarity_residual_is_machine_precision_when_identical(tmp_path, runs):
    """The two statistics were made equal in the synthetic data, so the residual has
    to be 0.

    A real run gives 4.4e-16 (floating-point level) -- either way the conclusion
    "this is not an independent check" is readable from the figure.
    """
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_stationarity(fs, runs, dt_star=5e-3, sample_interval_steps=400,
                                tau_trap_si=TRAP_SI["tau_trap_si"])
    assert "0.0e+00" in rec.caption or "e-1" in rec.caption


def test_stationarity_shows_mentions_independence_check(tmp_path, runs):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_stationarity(fs, runs, dt_star=5e-3, sample_interval_steps=400,
                                tau_trap_si=TRAP_SI["tau_trap_si"])
    assert "not** independent" in rec.shows


def test_displacement_figure_disclaims_per_step(tmp_path, fake_run):
    """★ A frame displacement is not a per-step displacement — the uniform noise's √3
    bound does not apply."""
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_frame_displacements(fs, fake_run["traj"], dt_star=5e-3,
                                       frame_interval_steps=16, sigma_star=491.358)
    assert "not the per-step displacement" in rec.caption
    assert "√3" in rec.caption


def test_equipartition_figure_plots_the_competing_hypothesis(tmp_path):
    """Drawing the competing hypothesis (exact) is what makes the decidability
    visible."""
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_equipartition_vs_dt(fs, {5e-3: (1.00577, 0.00446),
                                            2.5e-3: (1.00169, 0.00519)})
    assert "power" in rec.caption
    assert "exact" in rec.shows or "distinguishable" in rec.shows


def test_snapshots_title_carries_both_units(tmp_path, fake_run):
    fs = viz.FigureSet(tmp_path)
    rec = viz.plot_snapshots(fs, fake_run["traj"], dim=2, dt_star=5e-3,
                             frame_interval_steps=16,
                             tau_trap_si=TRAP_SI["tau_trap_si"],
                             l_trap_si=TRAP_SI["l_trap_si"])
    assert "nm" in rec.caption                      # ℓ_trap's physical value
    assert "independent replicas" in rec.caption    # the not-a-suspension warning


# =============================================================================
# Orchestrator
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
    """An exhaustive gate check — one empty entry is a §S6 violation."""
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
    """★ With no pair interaction g(r) is **not defined** — a different reason from
    'not implemented'."""
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, has_pair=False, **TRAP_SI)
    assert "is not defined" in fs.skipped["07_rdf.png"]


def test_with_pair_skips_rdf_because_it_is_unimplemented(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=2, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, has_pair=True, **TRAP_SI)
    assert "no RDF calculator" in fs.skipped["07_rdf.png"]


def test_3d_skips_snapshots_with_reason(tmp_path, runs):
    fs = viz.trap_diagnostics(
        tmp_path, runs, dim=3, dt_star=5e-3, sample_interval_steps=400,
        frame_interval_steps=16, sigma_star=491.358, **TRAP_SI)
    assert "fresnel" in fs.skipped["06_snapshots.png"]


def test_matplotlib_uses_headless_backend():
    """Anything but Agg tries to open a window and dies in CI or headless."""
    import matplotlib
    assert matplotlib.get_backend().lower() == "agg"


def test_figure_text_is_ascii_only(tmp_path, runs):
    """★ Hangul in a figure's text renders as tofu (□) — blocked at the source."""
    src = Path(viz.__file__).read_text(encoding="utf-8")
    # no string inside set_xlabel/set_ylabel/set_title/label= may contain Hangul
    import re
    bad = []
    for m in re.finditer(r"(set_xlabel|set_ylabel|set_title|suptitle)\(([^\n]*)", src):
        if re.search(r"[가-힣]", m.group(2)):
            bad.append(m.group(0)[:70])
    assert not bad, f"Hangul in figure text: {bad}"
