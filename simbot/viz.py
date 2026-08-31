"""S6 — figures. 0 lines of LLM.

## The two rules of this module

**① A figure without a caption cannot be created.** `save()` takes the caption as a
required argument and rejects an empty string. The master_plan §S6 gate says "an
uncaptioned figure is not accepted as an artefact", but checking after the fact is
already too late -- **it is blocked at creation time.**

**② A figure not drawn leaves a reason.** `FigureSet.skip()`. If there is no way to
tell whether the RDF was skipped "because there is no pair interaction" or "because
someone forgot", a missing figure becomes a silent omission. The same discipline as
requiring a reason to turn a gate off (`simbot.spec.Gate`).

## Dual axis labelling

Every time and length axis carries **the dimensionless value and the physical unit
together** (§S6-4). With only `t/τ_trap` there is no answering "so how many ms is
that", and with only `t [ms]` it cannot be compared against another paper.

## Text is English

matplotlib's default font has no Hangul glyphs (CLAUDE.md). Both the captions
(markdown) and **the text inside a figure** are English; in-figure text must never
be Korean, which `tests/test_s6_viz.py` checks by scanning the rendered artists.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")                      # headless. The import order matters
import matplotlib.pyplot as plt            # noqa: E402

from .analysis.trap import em_uniform_noise_excess_kurtosis, msd_model
from .estimators import euler_maruyama_trap_variance_bias

PLOT_STYLE = {"font.size": 9, "figure.dpi": 130, "axes.grid": True,
              "grid.alpha": 0.25, "legend.framealpha": 0.85,
              "savefig.bbox": "tight"}


# =============================================================================
# Captions are compulsory
# =============================================================================
@dataclass
class FigureRecord:
    name: str
    caption: str
    shows: str                 # "what is this figure meant to show" — not the caption
    path: Path


def alt_text(rec: FigureRecord, limit: int = 90) -> str:
    """Alt text for the markdown image — **cut to one line.**

    Putting the caption straight into `![...]` breaks the image syntax for a
    multi-line caption. The full caption goes below the figure as body text.
    """
    line = " ".join(rec.caption.split())
    return line if len(line) <= limit else line[: limit - 1].rstrip() + "…"


@dataclass
class FigureSet:
    """A set of figures. It owns the captions and the skip reasons together."""

    outdir: Path
    records: list[FigureRecord] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.outdir = Path(self.outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)

    def save(self, fig, name: str, caption: str, shows: str) -> FigureRecord:
        """Saves a figure. **Refuses when `caption` or `shows` is empty.**

        Args:
            caption: the text under the figure. It has to include the numbers and
                the verdict.
            shows: **what this figure is meant to show** (the §S6 gate). If the
                caption is "the MSD curve", shows is "does the MSD follow a single
                exponential".
        """
        if not caption.strip():
            raise ValueError(f"{name}: no caption — an uncaptioned figure is not "
                             f"accepted as an artefact (master_plan §S6 gate)")
        if not shows.strip():
            raise ValueError(f"{name}: no `shows` — a figure that cannot answer "
                             f"'what is this meant to show' is useless for diagnosis")
        stem = name if name.endswith(".png") else f"{name}.png"
        path = self.outdir / stem
        fig.savefig(path)
        plt.close(fig)
        rec = FigureRecord(name=stem, caption=caption.strip(), shows=shows.strip(),
                           path=path)
        self.records.append(rec)
        return rec

    def skip(self, name: str, reason: str) -> None:
        """A figure not drawn, and why. **Skipping without a reason is refused.**"""
        if not reason.strip():
            raise ValueError(f"{name}: no reason for skipping — a missing figure "
                             f"would be indistinguishable between 'not applicable' "
                             f"and 'forgotten'")
        self.skipped[name] = reason.strip()

    @property
    def captions(self) -> dict[str, str]:
        """The shape passed straight to `report.ReportInputs.figures`."""
        return {r.name: r.caption for r in self.records}

    def figures_md(self) -> str:
        """`06_figures.md` — the figure list + captions + what was skipped."""
        out = ["# S6 FIGURES", "",
               f"{len(self.records)} figures · {len(self.skipped)} skipped", "",
               "> Text inside a figure is English — matplotlib's default font has no "
               "Hangul glyphs.",
               "> Every time and length axis is dual-labelled with the dimensionless "
               "value and the physical unit (§S6-4).", ""]
        for r in self.records:
            out += [f"## {r.name}", "",
                    f"**What this figure is meant to show** {r.shows}", "",
                    # keep the alt text short — a multi-line caption inside `![...]`
                    # breaks the markdown image syntax
                    f"![{alt_text(r)}](figs/{r.name})", "",
                    r.caption, ""]
        if self.skipped:
            out += ["---", "", "## Figures not drawn, and why", "",
                    "| figure | reason |", "|---|---|"]
            out += [f"| `{n}` | {why} |" for n, why in sorted(self.skipped.items())]
            out += ["", "> No figure was skipped without a reason — "
                        "`FigureSet.skip()` requires one.", ""]
        return "\n".join(out)


# =============================================================================
# Dual axes — dimensionless above, SI below
# =============================================================================
def add_si_axis(ax, scale_si: float, unit: str, *, axis: str = "x",
                label: str = "", si_multiplier: float = 1.0):
    """Attaches an SI secondary axis to a dimensionless one.

    Args:
        scale_si: the SI value corresponding to dimensionless 1
            (`τ_trap = 8.064e-3 s`, say)
        unit: the display-unit string (`"ms"`, say)
        si_multiplier: SI base unit → display unit multiplier (`1e3` for s→ms)
    """
    f = scale_si * si_multiplier
    fwd, inv = (lambda v: np.asarray(v) * f), (lambda v: np.asarray(v) / f)
    if axis == "x":
        sec = ax.secondary_xaxis("top", functions=(fwd, inv))
        sec.set_xlabel(label or f"[{unit}]")
    else:
        sec = ax.secondary_yaxis("right", functions=(fwd, inv))
        sec.set_ylabel(label or f"[{unit}]")
    return sec


# =============================================================================
# The mandatory diagnostic set — harmonic trap
# =============================================================================
def plot_msd(fs: FigureSet, runs: dict, *, dim: int, dt_star: float,
             tau_trap_si: float, l_trap_si: float, name: str = "01_msd") -> FigureRecord:
    """MSD log–log + the analytic solution + a free-diffusion reference line.

    The free-diffusion reference (`2dD₀t`) is what makes **whether the short-time
    limit is right** visible -- it lets the eye separate the case where only the
    plateau matches and the short time is wrong.
    """
    fig, ax = plt.subplots(figsize=(5.4, 4.0))
    curves = []
    for i, (label, r) in enumerate(sorted(runs.items())):
        t = np.asarray(r["lags_steps"], dtype=float) * dt_star
        m = np.asarray(r["msd"], dtype=float)
        curves.append((t, m))
        ax.loglog(t, m, "o", ms=2.6, alpha=0.5, color="C0",
                  label=f"simulation ({len(runs)} seeds)" if i == 0 else None)
    t_all = curves[0][0]
    tt = np.geomspace(t_all[0], t_all[-1], 200)
    ax.loglog(tt, msd_model(tt, 2 * dim, 1.0), "k-", lw=1.7,
              label=r"analytic  $2d(1-e^{-t/\tau})$")
    ax.loglog(tt, 2 * dim * tt, "r--", lw=1.0, alpha=0.75,
              label=r"free diffusion  $2dD_0t$")
    ax.axhline(2 * dim, color="gray", ls=":", lw=1)
    ax.set_xlabel(r"lag  $t/\tau_{\rm trap}$")
    ax.set_ylabel(r"MSD  $\langle\Delta r^{*2}\rangle$   [$\ell_{\rm trap}^2$]")
    ax.set_title(f"{dim}D harmonic trap MSD   (plateau $=2d={2*dim}$)")
    ax.legend(fontsize=7, loc="lower right")
    add_si_axis(ax, tau_trap_si, "ms", axis="x", label="lag  $t$  [ms]",
                si_multiplier=1e3)
    add_si_axis(ax, l_trap_si**2 * 1e18, "nm²", axis="y",
                label=r"MSD  [nm$^2$]")
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"Harmonic-trap MSD ({dim}D, {len(runs)} seeds). It follows the "
                 f"free-diffusion slope 1 at short times and saturates at "
                 f"`2d = {2*dim}` for `t ≳ τ_trap`. The top and right axes are in "
                 f"physical units (converted by `τ_trap` and `ℓ_trap²`)."),
        shows=("does the MSD follow the single exponential `2d(1−e^{−t/τ})`, and "
               "does the short-time limit agree with free diffusion — this rules out "
               "the case where only the plateau matches"))


def plot_equipartition_vs_dt(fs: FigureSet, by_dt: dict[float, tuple[float, float]],
                             *, name: str = "02_equipartition_dt") -> FigureRecord:
    """`⟨x*²⟩` vs `dt*` + the Euler–Maruyama bias curve + an exact-scheme reference.

    ★ **The competing hypothesis (exact) is drawn too.** With only the EM curve the
      figure merely "looks right", and the fact that the gap between the two curves
      is smaller than the statistical error (= undecidable) stays invisible.
    """
    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    dts = np.geomspace(min(by_dt) / 3, max(by_dt) * 3, 200)
    ax.plot(dts, [1 + euler_maruyama_trap_variance_bias(x) for x in dts], "k-",
            lw=1.6, label=r"Euler-Maruyama  $1/(1-\Delta t^*/2)$")
    ax.axhline(1.0, color="r", ls="--", lw=1.2, label="exact scheme (no bias)")
    xs = sorted(by_dt)
    ax.errorbar(xs, [by_dt[x][0] for x in xs], yerr=[by_dt[x][1] for x in xs],
                fmt="o", color="C0", capsize=3, ms=5, label="measured (seed ensemble)")
    ax.set_xscale("log")
    ax.set_xlabel(r"$\Delta t^* = \Delta t/\tau_{\rm trap}$")
    ax.set_ylabel(r"$\langle x^{*2}\rangle$")
    ax.set_title("Equipartition and integrator bias")
    ax.legend(fontsize=7)
    fig.tight_layout()
    worst = max(by_dt, key=lambda d: abs(1 + euler_maruyama_trap_variance_bias(d) - 1))
    gap = euler_maruyama_trap_variance_bias(worst)
    se = by_dt[worst][1]
    return fs.save(
        fig, name,
        caption=(f"Equipartition and the integrator bias. Where the measured points "
                 f"sit between the EM curve (black) and exact (red) is what "
                 f"discriminates the scheme. At `dt* = {worst:g}` the gap between "
                 f"the two hypotheses is `{gap:.3%}` and the statistical error is "
                 f"`±{se:.3%}` — a power of `{gap/se:.2f}σ`."),
        shows=("does the measured `⟨x*²⟩` reproduce the Euler–Maruyama bias, and is "
               "it **distinguishable** from the exact scheme — the gap between the "
               "curves against the error-bar size *is* the power"))


def plot_position_distribution(fs: FigureSet, traj: np.ndarray, *, dt_star: float,
                               dim: int, l_trap_si: float,
                               frame_interval_steps: int,
                               decorrelation_tau: float = 2.0,
                               name: str = "03_distribution") -> FigureRecord:
    """The equilibrium position distribution vs Gaussian (log y).

    ★ Drawn on log y -- a departure in the tails is invisible on a linear axis. The
      tails are the point of interest because HOOMD's noise is uniform
      (findings §2).

    ⚠ **Only independent frames are used**, the same discipline as
      `analysis.trap.check_position_distribution`. Using every correlated frame
      inflates the sample count, makes the kurtosis error bar smaller than it is,
      and bakes a number into the caption that differs from the verified value in
      `metrics.json`. The same trap was hit with the KS test on 2026-07-28
      ([[ks-test-needs-independent-samples]]).
    """
    tr = np.asarray(traj, dtype=np.float64)
    frames_per_tau = 1.0 / (frame_interval_steps * dt_star)
    step = max(1, int(np.ceil(decorrelation_tau * frames_per_tau)))
    indep = tr[::step]
    x = indep.reshape(-1)
    x = (x - x.mean()) / x.std(ddof=1)
    kurt = float(np.mean(x**4))
    kurt_se = float(np.sqrt(24.0 / x.size))
    pred = 3.0 + em_uniform_noise_excess_kurtosis(dt_star)
    sigma_away = abs(kurt - pred) / kurt_se

    fig, ax = plt.subplots(figsize=(5.2, 3.7))
    ax.hist(x, bins=140, density=True, alpha=0.6, color="C0",
            label=f"simulation ({indep.shape[0]} independent frames)")
    xs = np.linspace(-5, 5, 400)
    ax.plot(xs, np.exp(-xs**2 / 2) / np.sqrt(2 * np.pi), "k-", lw=1.6,
            label=r"$\mathcal{N}(0,1)$  (Boltzmann)")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    ax.set_xlabel(r"$x/\sigma_{\rm measured}$   (normalised)")
    ax.set_ylabel("PDF")
    ax.set_title(f"{dim}D component positions   "
                 f"kurtosis {kurt:.4f}$\\pm${kurt_se:.4f} vs predicted {pred:.4f}")
    ax.legend(fontsize=7)
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"The equilibrium position distribution (normalized by the measured "
                 f"standard deviation — this looks at the **shape, not the width**; "
                 f"the width is checked separately by `⟨x²⟩`). Kurtosis "
                 f"`{kurt:.4f} ± {kurt_se:.4f}` against the prediction "
                 f"`3 − 1.2 dt* = {pred:.4f}` → `{sigma_away:.2f}σ`. "
                 f"★ **Exactly 3.000 would be the suspicious result** — HOOMD's "
                 f"noise is uniform, so a residual non-Gaussianity of order `dt*` "
                 f"remains. "
                 f"⚠ Only frames separated by at least "
                 f"`{decorrelation_tau:g} τ_trap` were used — **{indep.shape[0]} "
                 f"independent frames**; correlated frames inflate the sample count "
                 f"and make the error bar falsely small."),
        shows=("does the uniform noise approach a Gaussian by the CLT, and does it "
               "**reproduce the known residual deviation** — a stronger check than "
               "blindly confirming Gaussianity"))


def plot_stationarity(fs: FigureSet, runs: dict, *, dt_star: float,
                      sample_interval_steps: int, tau_trap_si: float,
                      name: str = "04_stationarity") -> FigureRecord:
    """`⟨x*²⟩` and `kT_conf` time series — confirming stationarity (equilibrium).

    This is where a thermo time series belongs. HOOMD `Brownian`'s kinetic
    temperature cannot depart systematically and is therefore useless
    (findings §1), so the **configurational temperature** and the equipartition
    statistic are what get plotted.
    """
    fig, axes = plt.subplots(2, 1, figsize=(5.8, 4.6), sharex=True,
                             gridspec_kw={"height_ratios": [3, 1]})
    resid_max = 0.0
    for i, (label, r) in enumerate(sorted(runs.items())):
        v = np.asarray(r["indep_var"], dtype=float)
        k = np.asarray(r["indep_kT"], dtype=float)
        t = np.arange(v.size) * sample_interval_steps * dt_star
        axes[0].plot(t, v, "-o", ms=2.4, lw=0.8, alpha=0.75, color=f"C{i}",
                     label=f"seed {i+1}" if i < 4 else None)
        # ★ the two statistics are overlaid — if they are algebraically identical
        #   they have to lie exactly on top of each other
        axes[0].plot(t, k, "x", ms=4, mew=0.9, alpha=0.9, color="k",
                     label=r"$kT_{\rm conf}/kT$ (all seeds)" if i == 0 else None)
        resid_max = max(resid_max, float(np.abs(k / v - 1.0).max()))
        axes[1].plot(t, k / v - 1.0, "-o", ms=2.4, lw=0.8, alpha=0.75, color=f"C{i}")

    bias = euler_maruyama_trap_variance_bias(dt_star)
    axes[0].axhline(1 + bias, color="k", ls="-", lw=1.4,
                    label=rf"EM prediction $={1+bias:.5f}$")
    axes[0].axhline(1.0, color="r", ls="--", lw=1.0, label="exact scheme $=1$")
    axes[0].set_ylabel(r"$\langle x^{*2}\rangle$  and  $kT_{\rm conf}/kT$")
    axes[0].set_title("Stationarity — and the algebraic identity of the two probes")
    axes[0].legend(fontsize=6.5, ncol=2, loc="upper left")

    axes[1].axhline(0.0, color="k", lw=1.0)
    axes[1].set_ylabel(r"$\dfrac{kT_{\rm conf}}{\langle x^{*2}\rangle}-1$", fontsize=8)
    axes[1].set_xlabel(r"$t/\tau_{\rm trap}$")
    span = max(resid_max, 1e-16) * 3
    axes[1].set_ylim(-span, span)
    axes[1].text(0.98, 0.85, f"max $|{{\\rm resid}}|$ = {resid_max:.1e}",
                 transform=axes[1].transAxes, ha="right", va="top", fontsize=7)
    add_si_axis(axes[0], tau_trap_si, "ms", axis="x", label="$t$  [ms]",
                si_multiplier=1e3)
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"**Top** `⟨x*²⟩` per independent sample time (colour) and the "
                 f"configurational temperature `kT_conf` (black ×). Fluctuating with "
                 f"no trend means equilibrium was reached. The EM bias "
                 f"`{bias:.3%}` is far below the per-sample fluctuation (of order "
                 f"`±5 %`) and is not visible here -- which is why a seed-ensemble "
                 f"mean is needed.\n\n"
                 f"**Bottom** the relative residual between the two statistics, at "
                 f"most `{resid_max:.1e}` — floating-point level. ⚠ **In a pure "
                 f"harmonic trap `kT_conf` is algebraically identical to the "
                 f"`⟨x²⟩` check** (foretold by trap card §6). This figure does not "
                 f"assert that, it **shows** it — becoming an independent check "
                 f"requires a pair interaction."),
        shows=("① is the sample stationary — a trend means insufficient "
               "equilibration ② does the statistic actually fluctuate (if not it is "
               "an arithmetic identity) ③ is `kT_conf` independent of `⟨x*²⟩` — a "
               "zero residual means it is **not** independent"))


def plot_frame_displacements(fs: FigureSet, traj: np.ndarray, *, dt_star: float,
                             frame_interval_steps: int, sigma_star: float,
                             name: str = "05_displacements") -> FigureRecord:
    """Frame-to-frame displacement distribution — an after-the-fact guard check.

    ⚠ **This is not the per-step displacement.** The trajectory is stored every
      `frame_interval_steps` steps, so the uniform noise's structural bound
      `max/σ = √3` does not apply here (a sum of many steps approaches a Gaussian by
      the CLT). The per-step check is done by the runtime guard
      (`simbot.guards.check_step_displacements`).
    """
    tr = np.asarray(traj, dtype=np.float64)
    d = np.linalg.norm(tr[1:] - tr[:-1], axis=2).reshape(-1)
    n_steps = frame_interval_steps
    expected_rms = np.sqrt(2 * tr.shape[2] * n_steps * dt_star)   # free-diffusion cap

    fig, ax = plt.subplots(figsize=(5.0, 3.6))
    ax.hist(d, bins=120, density=True, alpha=0.65, color="C2",
            label=f"measured ({n_steps} steps/frame)")
    ax.axvline(float(d.max()), color="r", ls="--", lw=1.2,
               label=f"max = {d.max():.3f}")
    ax.axvline(expected_rms, color="k", ls=":", lw=1.4,
               label=rf"free-diffusion rms $\sqrt{{2d\,n\Delta t^*}}$ = {expected_rms:.3f}")
    ax.set_xlabel(r"frame displacement  $|\Delta r|/\ell_{\rm trap}$")
    ax.set_ylabel("PDF")
    ax.set_title(f"Frame displacements   (max = {d.max()/sigma_star:.2e} $\\sigma$)")
    ax.legend(fontsize=7)
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"Frame-to-frame displacement distribution (`{n_steps}`-step "
                 f"interval). Maximum `{d.max():.3f} ℓ_trap` = "
                 f"`{d.max()/sigma_star:.2e} σ` — a factor "
                 f"{0.1/(d.max()/sigma_star):.0f} or more below the `0.1 σ` overlap "
                 f"threshold. ⚠ **This is not the per-step displacement.** The "
                 f"uniform noise's `max/σ = √3` bound does not apply here (it is a "
                 f"sum of many steps)."),
        shows=("does the displacement show any sign of blow-up, and is it inside the "
               "free-diffusion bound — confinement means it has to be smaller than "
               "free diffusion"))


def plot_snapshots(fs: FigureSet, traj: np.ndarray, *, dim: int, l_trap_si: float,
                   dt_star: float, frame_interval_steps: int, tau_trap_si: float,
                   name: str = "06_snapshots") -> FigureRecord:
    """A triptych of the initial, middle and final snapshots + `ℓ_trap` circles.

    2D only. 3D needs `fresnel`, which is not installed yet.
    """
    tr = np.asarray(traj, dtype=np.float64)
    idx = [0, tr.shape[0] // 2, tr.shape[0] - 1]
    fig, axes = plt.subplots(1, 3, figsize=(9.6, 3.4), sharex=True, sharey=True)
    lim = float(np.abs(tr[:, :, :2]).max()) * 1.1
    for ax, i in zip(axes, idx):
        p = tr[i]
        ax.scatter(p[:, 0], p[:, 1], s=3, alpha=0.45, color="C0", linewidths=0)
        for rad, ls in ((1.0, "-"), (2.0, "--"), (3.0, ":")):
            ax.add_patch(plt.Circle((0, 0), rad, fill=False, color="k", lw=1.0,
                                    ls=ls, alpha=0.7))
        t_star = i * frame_interval_steps * dt_star
        t_ms = t_star * tau_trap_si * 1e3
        # the time is dual-labelled too — a snapshot cannot carry a secondary axis,
        # so it goes in the title
        ax.set_title(f"$t = {t_star:.2f}\\,\\tau_{{\\rm trap}}$"
                     f"  $= {t_ms:.2f}$ ms")
        ax.set_xlabel(r"$x/\ell_{\rm trap}$")
        ax.set_aspect("equal")
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
    axes[0].set_ylabel(r"$y/\ell_{\rm trap}$")
    fig.suptitle(r"Snapshots with $1,2,3\,\ell_{\rm trap}$ circles   "
                 f"($\\ell_{{\\rm trap}} = {l_trap_si*1e9:.2f}$ nm)", fontsize=10)
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"Initial, middle and final snapshots. The circles are "
                 f"`1, 2, 3 ℓ_trap` (`ℓ_trap = {l_trap_si*1e9:.2f} nm`). If the "
                 f"three distributions are indistinguishable by eye, equilibrium was "
                 f"reached. The {tr.shape[1]} particles are **independent replicas** "
                 f"in the same trap, not a suspension."),
        shows=("are the particles distributed isotropically about the trap centre, "
               "and are the three times indistinguishable (= steady state)"))


# =============================================================================
# The mandatory diagnostic set — 2D soft-repulsive (time-resolved)
# =============================================================================
#  Card: knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md
#  ★ The time axis is `τ_d = d²/D₀` and the secondary axis is in **minutes**
#    (`τ_d ≈ 38 min`). In ms there would be 6 spare digits and it would be
#    unreadable.

#  Coordination colours — 5 (red) / 6 (grey) / 7 (blue) is the standard scheme for
#  dislocation pairs.
#  ★ The point is putting 5 and 7 in **contrasting** colours: a dislocation shows 5
#    and 7 in equal numbers, so whether the two colours appear in equal counts
#    becomes a signature the eye can judge.
COORD_COLORS = {3: "#4a148c", 4: "#00838f", 5: "#d32f2f", 6: "#cfcfcf",
                7: "#1565c0", 8: "#f57c00", 9: "#6d4c41", 10: "#2e7d32"}
COORD_OTHER = "#000000"


def _coord_color(z: int) -> str:
    return COORD_COLORS.get(int(z), COORD_OTHER)


def plot_structure_timeseries(fs: FigureSet, series_by_A: dict, *,
                              tau_d_si: float, energy_by_A: dict | None = None,
                              name: str = "01_structure_timeseries"
                              ) -> FigureRecord:
    """The **time trajectories** of `ψ₆`, the defect fraction, the coordination
    populations and the 5-7 imbalance.

    Args:
        series_by_A: `{A: [HexOrderSeries, ...per seed]}`. The seed ensemble is the
            error bar -- **frame-to-frame scatter carries time correlation and is
            therefore not a statistical error.**
        energy_by_A: `{A: (t_star, [energy_pp array, ...per seed])}` (optional)

    ★ What this figure answers is not "what did it become" but **"when did it become
      that"**.
    """
    As = sorted(series_by_A)
    n_panel = 4 if energy_by_A is None else 5
    fig, axes = plt.subplots(n_panel, 1, figsize=(7.6, 2.05 * n_panel),
                             sharex=True)
    cmap = {A: f"C{i}" for i, A in enumerate(As)}

    def band(ax, t, mat, color, label):
        """Seed mean ± SE. `mat` is `(n_seeds, n_frames)`."""
        m = mat.mean(axis=0)
        ax.plot(t, m, color=color, lw=1.3, label=label)
        if mat.shape[0] > 1:
            se = mat.std(axis=0, ddof=1) / np.sqrt(mat.shape[0])
            ax.fill_between(t, m - se, m + se, color=color, alpha=0.22,
                            linewidth=0)
        return m

    finals = {}
    for A in As:
        ser = series_by_A[A]
        t = ser[0].t_star
        c = cmap[A]
        band(axes[0], t, np.array([s.psi6_global for s in ser]), c, f"A = {A:g}")
        band(axes[1], t, np.array([s.defect_fraction for s in ser]), c, f"A = {A:g}")
        band(axes[2], t, np.array([s.coord_kinds.astype(float) for s in ser]), c,
             f"A = {A:g}")
        band(axes[3], t, np.array([s.five_seven_balance for s in ser]), c,
             f"A = {A:g}")
        finals[A] = float(np.mean([s.psi6_global[-1] for s in ser]))
    if energy_by_A is not None:
        for A in As:
            t, mats = energy_by_A[A]
            band(axes[4], t, np.array(mats), cmap[A], f"A = {A:g}")
        axes[4].set_ylabel(r"$U/N$  [$k_BT$]")
        axes[4].set_yscale("log")

    axes[0].set_ylabel(r"$|\langle\psi_6\rangle|$")
    axes[1].set_ylabel("defect fraction")
    axes[2].set_ylabel("coord. kinds\n(frac > 0.5%)")
    axes[3].set_ylabel(r"$|n_5-n_7|/(n_5+n_7)$")
    axes[-1].set_xlabel(r"$t/\tau_d$")
    axes[0].legend(ncol=len(As), fontsize=8, loc="upper left")
    add_si_axis(axes[0], tau_d_si, "min", axis="x", label="t [min]",
                si_multiplier=1.0 / 60.0)
    for ax in axes:
        ax.set_xlim(float(series_by_A[As[0]][0].t_star[0]),
                    float(series_by_A[As[0]][0].t_star[-1]))
    fig.suptitle(f"Time-resolved structure   "
                 f"($\\tau_d = {tau_d_si/60:.1f}$ min)", fontsize=10)
    fig.tight_layout()

    fin = " · ".join(f"`A={A:g}` → `{finals[A]:.3f}`" for A in As)
    return fs.save(
        fig, name,
        caption=(f"From the top: global `ψ₆` · Voronoi defect fraction · number of "
                 f"coordination kinds · 5-7 imbalance"
                 + (" · energy per particle" if energy_by_A is not None else "")
                 + f". The bands are the **seed-ensemble SE** (not the "
                 f"frame-to-frame scatter — that carries time correlation). `ψ₆` at "
                 f"the last frame: {fin}. The secondary axis is real time, with "
                 f"`τ_d = {tau_d_si/60:.1f} min`. "
                 f"**A 5-7 imbalance falling to 0 means the defects have organized "
                 f"into dislocation pairs**, and the number of coordination kinds "
                 f"dropping 6→3 means the liquid-like defects are gone."),
        shows=("when the structure gets made — not the final value but the time it is "
               "reached and the direction of drift. A curve that has not flattened is "
               "evidence that equilibration cannot be claimed at this run length"))


def plot_rdf_evolution(fs: FigureSet, windows_by_A: dict, *, tau_d_si: float,
                       sigma_over_d: float | None = None,
                       name: str = "02_rdf_evolution") -> FigureRecord:
    """`g(r)` per time window — one panel per `A`, one curve per window.

    Args:
        sigma_over_d: given, a vertical line is drawn at `r = σ`. That line makes
            **whether the reference discs overlap** visible to the eye (there is no
            hard core, so they can).
    """
    As = sorted(windows_by_A)
    fig, axes = plt.subplots(1, len(As), figsize=(3.5 * len(As), 3.3),
                             sharey=True)
    axes = np.atleast_1d(axes)
    overlap_note = ""
    for ax, A in zip(axes, As):
        w = windows_by_A[A]
        n_w = w.g.shape[0]
        for j in range(n_w):
            frac = j / max(n_w - 1, 1)
            ax.plot(w.r, w.g[j], lw=1.2, color=plt.cm.viridis(0.12 + 0.78 * frac),
                    label=f"{w.t_lo[j]:.0f}–{w.t_hi[j]:.0f}")
        ax.axhline(1.0, color="k", lw=0.7, ls=":", alpha=0.6)
        if sigma_over_d is not None:
            ax.axvline(sigma_over_d, color="#d32f2f", lw=1.1, ls="--", alpha=0.9)
        ax.set_title(f"$A = {A:g}$")
        ax.set_xlabel("$r/d$")
        ax.set_xlim(0.0, min(4.0, float(w.r[-1])))
        ax.legend(fontsize=7, title=r"$t/\tau_d$", title_fontsize=7)
    axes[0].set_ylabel("$g(r)$")
    if sigma_over_d is not None:
        overlap_note = (f" The red dashed line is `r = σ = {sigma_over_d:.3f} d`, "
                        f"the distance at which the reference discs touch. "
                        f"**`g(r) > 0` inside it means the 5 µm discs overlap** "
                        f"(`A/r³` has no hard core, so the model allows it).")
    fig.suptitle("RDF vs time window", fontsize=10)
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"`g(r)` per time window. The darker the window, the later the time "
                 f"(`τ_d = {tau_d_si/60:.1f} min`). Windows lying on top of each "
                 f"other mean the structure stopped changing over that stretch; "
                 f"windows still separating mean it is still in progress."
                 + overlap_note),
        shows="when `g(r)`'s first peak grows and when it stops (and whether it does)")


def plot_early_transient(fs: FigureSet, series_by_A: dict, *, tau_d_si: float,
                         fits: dict | None = None,
                         name: str = "06_early_transient") -> FigureRecord:
    """The early transient — **on a logarithmic time axis.** A stretch that coarse
    sampling cannot see.

    Args:
        fits: `{A: RelaxationFit}` (optional). Only the `converged` ones get a curve
            overlaid; a rejected one has its **rejection reason** written into the
            legend -- a relaxation that is not there is never drawn.

    ★ `t = 0` cannot go on a log axis, so **the initial placement before the first
      frame is marked with a star in the left margin.** Drop the initial placement
      and "relaxed from what" disappears.
    """
    As = sorted(series_by_A)
    fig, axes = plt.subplots(2, 1, figsize=(7.4, 5.4), sharex=True)
    for i, A in enumerate(As):
        sers = series_by_A[A]
        t = sers[0].t_star
        c = f"C{i}"
        for ax, attr, lbl in ((axes[0], "defect_fraction", "defect fraction"),
                              (axes[1], "psi6_global", r"$|\langle\psi_6\rangle|$")):
            mat = np.array([getattr(s, attr) for s in sers])
            m = mat.mean(axis=0)
            pos = t > 0
            ax.plot(t[pos], m[pos], color=c, lw=1.0, alpha=0.55)
            if mat.shape[0] > 1:
                se = mat.std(axis=0, ddof=1) / np.sqrt(mat.shape[0])
                ax.fill_between(t[pos], (m - se)[pos], (m + se)[pos], color=c,
                                alpha=0.18, linewidth=0)
            #  t=0 (the initial placement) is off the log axis — a star at the far left
            if (~pos).any():
                ax.plot([t[pos][0] * 0.55], [m[~pos][0]], marker="*",
                        markersize=11, color=c, linestyle="")
            ax.set_ylabel(lbl)
        f = (fits or {}).get(A)
        label = f"A = {A:g}"
        if f is not None and f.converged:
            tt = np.logspace(np.log10(max(t[t > 0][0], 1e-4)), np.log10(t[-1]), 200)
            axes[0].plot(tt, f.y_inf + (f.y0 - f.y_inf) * np.exp(-tt / f.tau),
                         color=c, lw=2.0, ls="--")
            label += (f"  ($\\tau$ = {f.tau:.3g} $\\tau_d$"
                      f" = {f.tau*tau_d_si:.0f} s)")
        elif f is not None:
            label += "  (no relaxation resolved)"
        axes[0].plot([], [], color=c, lw=1.6, label=label)

    axes[0].set_xscale("log")
    axes[0].legend(fontsize=8, loc="best")
    axes[1].set_xlabel(r"$t/\tau_d$   (log)")
    add_si_axis(axes[0], tau_d_si, "s", axis="x", label="t [s]")
    fig.suptitle("Early transient (log time).  "
                 r"$\star$ = initial configuration ($t=0$)", fontsize=10)
    fig.tight_layout()

    lines = []
    for A in As:
        f = (fits or {}).get(A)
        if f is None:
            continue
        if f.converged:
            lines.append(f"`A={A:g}`: `τ = {f.tau:.3g} ± {f.tau_se:.3g} τ_d` "
                         f"(`{f.tau*tau_d_si:.0f} s`), amplitude "
                         f"`{f.amplitude:+.3f}`")
        else:
            lines.append(f"`A={A:g}`: **no relaxation detected** — {f.note}")
    return fs.save(
        fig, name,
        caption=("The early transient on a **logarithmic time axis**. The stars are "
                 "the initial placement (`t = 0`, the placement whose "
                 "`min_sep = 0.8 d` was enforced by rejection sampling), and **the "
                 "stars of the three `A` coincide exactly** -- the initial placement "
                 "does not depend on `A` (same seed → same placement), so the three "
                 "curves **diverge from the same point.** The curves are "
                 "seed-ensemble means ± SE. The dashed lines are single-exponential "
                 "fits, and **nothing is drawn when the relaxation amplitude fails "
                 "to exceed twice the noise** -- because fitting an exponential to "
                 "steady-state fluctuation always yields some `τ`. "
                 + (" · ".join(lines) if lines else "")),
        shows=("does a relaxation from the initial placement to the steady state "
               "actually exist, and if so how long is it — under coarse sampling "
               "(0.2 τ_d) it already looks finished by the first frame"))


def plot_seed_convergence(fs: FigureSet, conv: dict, *, curves: dict,
                          tau_d_si: float, threshold_sigma: float = 3.0,
                          preregistered: list | None = None,
                          name: str = "01_seed_convergence") -> FigureRecord:
    """A triptych of `τ(k)` convergence, `σ(k)`, and the relaxation curve.

    Args:
        conv: `{"k": [...], "A0.1": {"tau": [...], "se": [...]}, "A1": {...},
                "diff": [...], "se_diff": [...], "sigma": [...]}`
            — these have to be **nested subsets of the same data**. Then the slope of
            the curve represents only **the estimator's `k` dependence**, not a
            difference in the data.
        curves: `{A: (t_star, mean_defect, fit)}` the relaxation curve and fit at the
            final `k`.
        preregistered: `[{"k":…, "sigma_expected":…, "stage":…}]` the pre-registered
            expectations.

    ★ The point of this figure is **the mismatch between expected and actual.**
      Plotting the pre-registered expectation alongside is what makes "the low-seed
      pilot was optimistic" visible.
    """
    k = np.asarray(conv["k"], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.7))

    # --- ① tau vs k ---------------------------------------------------------
    #  ★ The error bars are the **robust interval** (bootstrap percentile). The SD is
    #    outlier-dominated at low seed count (82 τ_d at k=4 -- τ itself is 0.06),
    #    which destroys the axis and makes the figure unreadable.
    #    Points where the SD exceeds τ are drawn as **hollow markers** to mark "no
    #    error can be claimed".
    ax = axes[0]
    tau_all = []
    for i, A in enumerate(("A0.1", "A1")):
        tau = np.asarray(conv[A]["tau"], dtype=float)
        se = np.asarray(conv[A]["se"], dtype=float)
        se_r = np.asarray(conv[A].get("se_robust", se), dtype=float)
        usable = se <= np.abs(tau)          # an SD larger than τ is meaningless
        c = f"C{i}"
        ax.errorbar(k[usable], tau[usable], yerr=se_r[usable], marker="o",
                    ms=4, lw=1.3, color=c, capsize=2, label=f"$A = {A[1:]}$")
        if (~usable).any():
            ax.errorbar(k[~usable], tau[~usable], yerr=se_r[~usable],
                        marker="o", ms=5, lw=1.0, color=c, capsize=2,
                        markerfacecolor="white", ls="--", alpha=0.8)
        ax.axhline(tau[-1], color=c, lw=0.8, ls=":", alpha=0.6)
        tau_all += list(tau)
    lo, hi = min(tau_all), max(tau_all)
    pad = 0.25 * (hi - lo)
    ax.set_ylim(lo - pad, hi + pad)         # ★ the axis is set by τ, not the error bars
    ax.set_xscale("log")
    ax.set_xlabel("seeds $k$")
    ax.set_ylabel(r"$\tau_{\rm relax}\ [\tau_d]$")
    ax.set_title(r"① $\tau(A{=}1)$ drifts until $k \approx 512$", fontsize=9)
    ax.plot([], [], marker="o", ls="--", color="0.4", markerfacecolor="white",
            ms=5, label="bootstrap SD unusable")
    ax.legend(fontsize=7, loc="upper right")
    add_si_axis(ax, tau_d_si, "s", axis="y", label=r"$\tau$ [s]")

    # --- ② sigma vs k -------------------------------------------------------
    ax = axes[1]
    sig = np.asarray(conv["sigma"])
    ax.plot(k, sig, marker="o", ms=3.5, lw=1.4, color="C3", label="observed")
    ax.axhline(threshold_sigma, color="k", lw=1.1, ls="--",
               label=f"threshold {threshold_sigma:g}$\\sigma$")
    for j, p in enumerate(preregistered or []):
        ax.plot([p["k"]], [p["sigma_expected"]], marker="*", ms=15,
                color="C2", linestyle="", zorder=4,
                label="pre-registered expectation" if j == 0 else None)
        #  an arrow joins actual to expected — the mismatch is the point of the figure
        near = int(np.argmin(np.abs(k - p["k"])))
        ax.annotate("", xy=(k[near], sig[near]),
                    xytext=(p["k"], p["sigma_expected"]),
                    arrowprops=dict(arrowstyle="->", color="C2", lw=1.0,
                                    alpha=0.65, shrinkA=6, shrinkB=4))
        ax.annotate(f"stage {p['stage']}: {p['sigma_expected']:.2f}$\\sigma$",
                    xy=(p["k"], p["sigma_expected"]), fontsize=7,
                    xytext=(-4, 7), textcoords="offset points", ha="right",
                    color="C2")
    ax.set_xscale("log")
    ax.set_xlabel("seeds $k$")
    ax.set_ylabel(r"$|\Delta\tau| / SE_{\Delta}$   [$\sigma$]")
    ax.set_title("② expected vs achieved power", fontsize=9)
    ax.set_ylim(0, max(4.6, float(np.nanmax(sig)) * 1.15))
    ax.legend(fontsize=7, loc="upper left")

    # --- ③ the relaxation curve ---------------------------------------------
    ax = axes[2]
    for i, (A, (t, mean, fit)) in enumerate(sorted(curves.items())):
        c = f"C{i}"
        pos = np.asarray(t) > 0
        ax.plot(np.asarray(t)[pos], np.asarray(mean)[pos], color=c, lw=1.0,
                alpha=0.5)
        ax.plot([np.asarray(t)[pos][0] * 0.55], [np.asarray(mean)[~pos][0]],
                marker="*", ms=11, color=c, linestyle="")
        tt = np.logspace(np.log10(np.asarray(t)[pos][0]),
                         np.log10(np.asarray(t)[-1]), 200)
        ax.plot(tt, fit["y_inf"] + (fit["y0"] - fit["y_inf"])
                * np.exp(-tt / fit["tau"]), color=c, lw=2.0, ls="--",
                label=f"$A = {A}$  ($\\tau$ = {fit['tau']:.4f} $\\tau_d$)")
    ax.set_xscale("log")
    ax.set_xlabel(r"$t/\tau_d$   (log)")
    ax.set_ylabel("defect fraction")
    ax.set_title(r"③ relaxation at final $k$  ($\star$ = $t{=}0$)", fontsize=9)
    ax.legend(fontsize=7)
    fig.tight_layout()

    k_final = int(k[-1])
    return fs.save(
        fig, name,
        caption=(f"**① The `τ` estimate moves with the seed count.** These curves "
                 f"come from **nested subsets of the same trajectories**, so the "
                 f"slope is not a difference in the data but **the estimator's `k` "
                 f"dependence** -- i.e. the low-seed bias. The error bars are the "
                 f"seed bootstrap (not the `curve_fit` covariance). "
                 f"**② The mismatch between the pre-registered expectation (stars) "
                 f"and the actual (lines).** "
                 f"**③ The defect relaxation curve and single-exponential fit at the "
                 f"final `k = {k_final}`**; the stars are the initial placement "
                 f"(`t=0`) and the two `A` diverge from the same point. "
                 f"`τ_d = {tau_d_si/60:.1f} min`."),
        shows=("why a low-seed pilot gives an optimistic power — a biased difference "
               "and an underestimated error multiply together. And whether τ "
               "separates in the final sample"))


def plot_finite_size_scaling(fs: FigureSet, data: dict, *, local_data: dict,
                             exponent_liquid: float = 0.5,
                             exponent_hexatic: float = 0.0625,
                             name: str = "01_finite_size_scaling") -> FigureRecord:
    """`ψ₆` vs `N` (log-log) + the `N`-independence of a local quantity, as a pair.

    Args:
        data: `{A: {"N": [...], "psi6": [...], "se": [...], "fit": FiniteSizeExponent}}`
        local_data: `{A: {"N": [...], "defect": [...], "defect_se": [...]}}`

    ★ The point of this figure is **the slope**, not the absolute value -- in a
      finite system even a disordered one produces a `ψ₆` of order `1/√N`. The liquid
      slope (`−1/2`) and the KTHNY hexatic-boundary slope (`−1/16`) are laid down as
      guide lines so the eye can compare.
    """
    As = sorted(data)
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 3.9))

    ax = axes[0]
    for i, A in enumerate(As):
        d = data[A]
        N = np.asarray(d["N"], dtype=float)
        y = np.asarray(d["psi6"], dtype=float)
        se = np.asarray(d["se"], dtype=float)
        c = f"C{i}"
        ax.errorbar(N, y, yerr=se, marker="o", ms=5, lw=0, color=c, capsize=3,
                    elinewidth=1.2)
        f = d.get("fit")
        lab = f"$A = {A:g}$"
        if f is not None:
            lab += (f"   $p = {f.p:.2f}$" +
                    (f" $\\pm$ {f.p_se:.2f}" if np.isfinite(f.p_se) else "") +
                    f",  $\\eta_6 = {f.eta6:.2f}$")
            #  the fitted slope
            nn = np.array([N.min() * 0.9, N.max() * 1.1])
            ax.plot(nn, y[0] * (nn / N[0]) ** (-f.p), color=c, lw=1.5, alpha=0.9)
        ax.plot([], [], marker="o", ls="-", color=c, label=lab)
        #  the two hypothesis guide lines — anchored at the first point
        for p_ref, ls, alpha in ((exponent_liquid, "--", 0.55),
                                 (exponent_hexatic, ":", 0.55)):
            nn = np.array([N.min(), N.max()])
            ax.plot(nn, y[0] * (nn / N[0]) ** (-p_ref), color=c, lw=1.0,
                    ls=ls, alpha=alpha)
    ax.plot([], [], color="0.35", ls="--", lw=1.0,
            label=f"liquid  $p = {exponent_liquid:g}$")
    ax.plot([], [], color="0.35", ls=":", lw=1.0,
            label=f"hexatic bdry  $p = {exponent_hexatic:g}$")
    ax.set_xscale("log"); ax.set_yscale("log")
    ax.set_xlabel("$N$")
    ax.set_ylabel(r"$|\langle\psi_6\rangle|$")
    ax.set_title(r"① slope $p$ discriminates the phase, not $|\psi_6|$",
                 fontsize=9)
    ax.legend(fontsize=7, loc="lower left")

    ax = axes[1]
    for i, A in enumerate(sorted(local_data)):
        d = local_data[A]
        N = np.asarray(d["N"], dtype=float)
        ax.errorbar(N, d["defect"], yerr=d["defect_se"], marker="s", ms=5,
                    lw=1.2, color=f"C{i}", capsize=3, label=f"$A = {A:g}$")
    ax.set_xscale("log")
    ax.set_xlabel("$N$")
    ax.set_ylabel("defect fraction")
    ax.set_title("② local quantity should be flat in $N$", fontsize=9)
    ax.legend(fontsize=8)
    fig.tight_layout()

    lines = []
    for A in As:
        f = data[A].get("fit")
        if f is not None:
            lines.append(f"`A={A:g}`: `p = {f.p:.3f}`"
                         + (f" ± `{f.p_se:.3f}`" if np.isfinite(f.p_se) else "")
                         + f" → `η₆ = {f.eta6:.2f}`, reading **{f.reading}**")
    return fs.save(
        fig, name,
        caption=("**① `ψ₆`'s `N` dependence.** The dashed line is the liquid slope "
                 f"(`p = {exponent_liquid:g}`, `|⟨ψ₆⟩| ~ N^-1/2`) and the dotted "
                 f"line the KTHNY hexatic-liquid boundary "
                 f"(`p = {exponent_hexatic:g}`, `η₆ = 1/4`), each anchored at that "
                 f"`A`'s first point. **The slope discriminates the phase, not the "
                 f"absolute value** -- in a finite system even a disordered one "
                 f"produces a `ψ₆` of order `1/√N`. " + " · ".join(lines) +
                 ". **② The defect fraction is a local quantity and therefore has to "
                 "be flat in `N`** -- a slope would mean there was a finite-size "
                 "effect at `N=100`. "
                 "⚠ Two points cannot verify the power-law **form**; only the "
                 "exponent was extracted."),
        shows=("is `ψ₆`'s absolute value the finite-size floor or genuine "
               "quasi-long-range orientational order — and had the local quantity "
               "already converged at N=100"))


def plot_voronoi_timelapse(fs: FigureSet, frames: list, *, t_star: list,
                           tau_d_si: float, L_si: float, amplitude: float,
                           mean_defect_fraction: float | None = None,
                           defect_frame_sd: float | None = None,
                           name: str = "03_voronoi_timelapse") -> FigureRecord:
    """A time sequence of Voronoi tilings — the cells are coloured **by
    coordination**.

    The unimplemented item from card §10 ("a `voronoi plot` was requested for the
    figures but is not in `viz.py`").

    ★ Coordination 5 (red) and 7 (blue) appearing **in pairs** is a dislocation;
      5, 7, 4 and 8 scattered and mixed is a liquid. The defect **fraction** alone
      loses this distinction (card §8.2).
    """
    if len(frames) != len(t_star):
        raise ValueError(f"{len(frames)} frames vs {len(t_star)} times — "
                         f"the counts have to match")
    n = len(frames)
    fig, axes = plt.subplots(1, n, figsize=(2.75 * n, 3.15))
    axes = np.atleast_1d(axes)
    seen: set[int] = set()
    for ax, vf, t in zip(axes, frames, t_star):
        Lx, Ly = vf.Lx, vf.Ly
        for poly, z in zip(vf.polygons, vf.coordination):
            ax.add_patch(plt.Polygon(poly, closed=True, facecolor=_coord_color(z),
                                     edgecolor="white", lw=0.45, alpha=0.95))
            seen.add(int(z))
        ax.scatter(vf.positions[:, 0], vf.positions[:, 1], s=2.2, color="k",
                   zorder=3, linewidths=0)
        n_def = int(np.sum(vf.coordination != 6))
        ax.set_title(f"$t = {t:.1f}\\,\\tau_d = {t*tau_d_si/60:.0f}$ min\n"
                     f"defects {n_def}/{vf.coordination.size}", fontsize=8)
        ax.set_xlim(-Lx / 2, Lx / 2)
        ax.set_ylim(-Ly / 2, Ly / 2)
        ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([])
        ax.grid(False)
    handles = [plt.Line2D([], [], marker="s", ls="", markersize=7,
                          markerfacecolor=_coord_color(z), markeredgecolor="0.4",
                          label=f"z = {z}") for z in sorted(seen)]
    fig.legend(handles=handles, loc="lower center", ncol=min(len(handles), 8),
               fontsize=8, frameon=False, bbox_to_anchor=(0.5, -0.04))
    ref = ""
    if mean_defect_fraction is not None:
        n_mean = mean_defect_fraction * frames[0].coordination.size
        ref = f"   ensemble mean {n_mean:.0f}"
        if defect_frame_sd is not None:
            ref += f" $\\pm$ {defect_frame_sd*frames[0].coordination.size:.0f} (per frame)"
    fig.suptitle(f"Voronoi tessellation, $A = {amplitude:g}$   "
                 f"(box {L_si*1e6:.0f} $\\mu$m){ref}", fontsize=10)
    fig.tight_layout()

    #  ★ A single frame fluctuates — a difference between panels must not be read as
    #    relaxation.
    caveat = ""
    if mean_defect_fraction is not None and defect_frame_sd is not None:
        n = frames[0].coordination.size
        caveat = (f" ⚠ **Do not read a difference in defect count between panels as "
                  f"a change in time** — in this system the per-frame defect count "
                  f"fluctuates in steady state as "
                  f"`{mean_defect_fraction*n:.0f} ± {defect_frame_sd*n:.0f}` (the "
                  f"frame standard deviation). Most of the panel-to-panel difference "
                  f"is that fluctuation. A change in time has to be judged from the "
                  f"ensemble mean in the time-series figure.")
    return fs.save(
        fig, name,
        caption=(f"A time sequence of Voronoi tilings at `A = {amplitude:g}` "
                 f"(**one seed**). The colour is the coordination `z` — **red "
                 f"`z=5`, blue `z=7`, grey `z=6`**. The box is "
                 f"`{L_si*1e6:.0f} µm` square with periodic boundaries (cells "
                 f"extending outside the box are cells open toward a periodic "
                 f"neighbour and were not clipped -- clipping makes the areas "
                 f"wrong). Red and blue **adjacent, in pairs** is a dislocation; "
                 f"several colours scattered about is a liquid." + caveat),
        shows=("the **character** of the defects — the fraction alone does not "
               "distinguish a dislocation from a liquid. This figure shows the "
               "character; the basis for a change in time is the time-series figure"))


# =============================================================================
# Orchestrator
# =============================================================================
def trap_diagnostics(outdir: Path, runs: dict, *, dim: int, dt_star: float,
                     sample_interval_steps: int, frame_interval_steps: int,
                     sigma_star: float, tau_trap_si: float, l_trap_si: float,
                     by_dt: dict[float, tuple[float, float]] | None = None,
                     has_pair: bool = False) -> FigureSet:
    """The mandatory diagnostic set for the harmonic-trap system.

    Args:
        runs: `label → the contents of samples.npz`. The whole seed ensemble.
        by_dt: `dt* → (⟨x*²⟩ mean, SE)`. Only when several `dt*` were run. Without
            it, that figure is **skipped with a reason.**
        has_pair: whether there is a pair interaction. Without one the RDF is skipped
            with a reason.
    """
    plt.rcParams.update(PLOT_STYLE)
    fs = FigureSet(outdir)
    first = runs[sorted(runs)[0]]

    plot_msd(fs, runs, dim=dim, dt_star=dt_star, tau_trap_si=tau_trap_si,
             l_trap_si=l_trap_si)

    if by_dt and len(by_dt) >= 2:
        plot_equipartition_vs_dt(fs, by_dt)
    else:
        fs.skip("02_equipartition_dt.png",
                f"only {len(by_dt or {})} `dt*` values — the bias's `dt*` dependence "
                f"cannot be plotted. Running the dt ladder with `cli.py converge` "
                f"produces it")

    plot_position_distribution(fs, first["traj"], dt_star=dt_star, dim=dim,
                               l_trap_si=l_trap_si,
                               frame_interval_steps=frame_interval_steps)
    plot_stationarity(fs, runs, dt_star=dt_star,
                      sample_interval_steps=sample_interval_steps,
                      tau_trap_si=tau_trap_si)
    plot_frame_displacements(fs, first["traj"], dt_star=dt_star,
                             frame_interval_steps=frame_interval_steps,
                             sigma_star=sigma_star)

    if dim == 2:
        plot_snapshots(fs, first["traj"], dim=dim, l_trap_si=l_trap_si,
                       dt_star=dt_star, frame_interval_steps=frame_interval_steps,
                       tau_trap_si=tau_trap_si)
    else:
        fs.skip("06_snapshots.png",
                "a 3D snapshot needs `fresnel` ray tracing, which is not installed "
                "yet (techniques/env-log.md). A 2D projection would mislead")

    if has_pair:
        fs.skip("07_rdf.png",
                "there is a pair interaction but no RDF calculator "
                "(`analysis/structure.py`) yet — it will be built together with the "
                "trap+WCA system")
    else:
        fs.skip("07_rdf.png",
                "★ with no pair interaction `g(r)` is not defined. The particles are "
                "independent replicas in the same trap, so the relative-distance "
                "distribution is merely the difference of two Gaussians and carries "
                "no structural information")
    return fs
