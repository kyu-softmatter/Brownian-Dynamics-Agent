"""S6 — 그림. LLM 0줄.

## 이 모듈의 두 규칙

**① 캡션 없는 그림은 만들 수 없다.** `save()` 가 캡션을 필수 인자로 받고 빈 문자열을
거부한다. master_plan §S6 게이트는 "캡션 없는 그림은 산출물로 인정하지 않는다"인데,
사후에 검사하면 이미 그림이 생긴 뒤다 — **생성 시점에 막는다.**

**② 안 그린 그림에는 이유가 남는다.** `FigureSet.skip()`. RDF 를 그리지 않은 것이
"쌍 상호작용이 없어서"인지 "잊어서"인지 구별되지 않으면, 빠진 그림은 조용한 누락이 된다.
게이트를 끌 때 이유를 요구하는 것(`simbot.spec.Gate`)과 같은 규율이다.

## 이중 축 표기

모든 시간·길이 축에 **무차원 값과 물리 단위를 함께** 적는다 (§S6-4).
`t/τ_trap` 만 적으면 "그래서 몇 ms 인가"에 답할 수 없고, `t [ms]` 만 적으면
다른 논문과 비교할 수 없다.

## 글자는 영문

matplotlib 기본 폰트에 한글 글리프가 없다 (CLAUDE.md). 캡션(마크다운)은 한국어,
**그림 안 텍스트는 영문**이다.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import numpy as np

import matplotlib
matplotlib.use("Agg")                      # 헤드리스. import 순서가 중요하다
import matplotlib.pyplot as plt            # noqa: E402

from .analysis.trap import em_uniform_noise_excess_kurtosis, msd_model
from .estimators import euler_maruyama_trap_variance_bias

PLOT_STYLE = {"font.size": 9, "figure.dpi": 130, "axes.grid": True,
              "grid.alpha": 0.25, "legend.framealpha": 0.85,
              "savefig.bbox": "tight"}


# =============================================================================
# 캡션 강제
# =============================================================================
@dataclass
class FigureRecord:
    name: str
    caption: str
    shows: str                 # "무엇을 보이려는 그림인가" — 캡션과 다르다
    path: Path


def alt_text(rec: FigureRecord, limit: int = 90) -> str:
    """마크다운 이미지의 alt 텍스트 — **한 줄로 자른다.**

    캡션을 그대로 `![...]` 안에 넣으면 여러 줄 캡션에서 이미지 문법이 깨진다.
    전체 캡션은 그림 아래 본문으로 나간다.
    """
    line = " ".join(rec.caption.split())
    return line if len(line) <= limit else line[: limit - 1].rstrip() + "…"


@dataclass
class FigureSet:
    """그림 모음. 캡션과 '건너뛴 이유'를 함께 소유한다."""

    outdir: Path
    records: list[FigureRecord] = field(default_factory=list)
    skipped: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        self.outdir = Path(self.outdir)
        self.outdir.mkdir(parents=True, exist_ok=True)

    def save(self, fig, name: str, caption: str, shows: str) -> FigureRecord:
        """그림을 저장한다. **캡션과 `shows` 가 비어 있으면 거부한다.**

        Args:
            caption: 그림 아래에 붙는 설명. 수치와 판정을 포함해야 한다.
            shows: 이 그림이 **무엇을 보이려는 것인가** (§S6 게이트).
                캡션이 "MSD 곡선"이면 shows 는 "MSD 가 단일지수를 따르는가"다.
        """
        if not caption.strip():
            raise ValueError(f"{name}: 캡션이 없다 — 캡션 없는 그림은 산출물로 "
                             f"인정하지 않는다 (master_plan §S6 게이트)")
        if not shows.strip():
            raise ValueError(f"{name}: `shows` 가 없다 — '무엇을 보이려는 그림인가'에 "
                             f"답하지 못하는 그림은 진단에 쓸 수 없다")
        stem = name if name.endswith(".png") else f"{name}.png"
        path = self.outdir / stem
        fig.savefig(path)
        plt.close(fig)
        rec = FigureRecord(name=stem, caption=caption.strip(), shows=shows.strip(),
                           path=path)
        self.records.append(rec)
        return rec

    def skip(self, name: str, reason: str) -> None:
        """그리지 않은 그림과 그 이유. **이유 없이 건너뛸 수 없다.**"""
        if not reason.strip():
            raise ValueError(f"{name}: 건너뛴 이유가 없다 — 빠진 그림이 "
                             f"'해당 없음'인지 '잊음'인지 구별되지 않는다")
        self.skipped[name] = reason.strip()

    @property
    def captions(self) -> dict[str, str]:
        """`report.ReportInputs.figures` 에 그대로 넘기는 형태."""
        return {r.name: r.caption for r in self.records}

    def figures_md(self) -> str:
        """`06_figures.md` — 그림 목록 + 캡션 + 건너뛴 것."""
        out = ["# S6 FIGURES", "",
               f"그림 {len(self.records)}장 · 건너뜀 {len(self.skipped)}건", "",
               "> 그림 안 텍스트는 영문이다 — matplotlib 기본 폰트에 한글 글리프가 없다.",
               "> 모든 시간·길이 축은 무차원 값과 물리 단위를 이중 표기한다 (§S6-4).", ""]
        for r in self.records:
            out += [f"## {r.name}", "",
                    f"**무엇을 보이려는 그림인가** {r.shows}", "",
                    # alt 텍스트는 짧게 — 여러 줄 캡션을 `![...]` 안에 넣으면
                    # 마크다운 이미지 문법이 깨진다
                    f"![{alt_text(r)}](figs/{r.name})", "",
                    r.caption, ""]
        if self.skipped:
            out += ["---", "", "## 그리지 않은 그림과 이유", "",
                    "| 그림 | 이유 |", "|---|---|"]
            out += [f"| `{n}` | {why} |" for n, why in sorted(self.skipped.items())]
            out += ["", "> 이유 없이 건너뛴 그림은 없다 — `FigureSet.skip()` 이 "
                        "이유를 필수로 요구한다.", ""]
        return "\n".join(out)


# =============================================================================
# 이중 축 — 무차원 위, SI 아래
# =============================================================================
def add_si_axis(ax, scale_si: float, unit: str, *, axis: str = "x",
                label: str = "", si_multiplier: float = 1.0):
    """무차원 축에 SI 보조축을 붙인다.

    Args:
        scale_si: 무차원 1 에 해당하는 SI 값 (예: `τ_trap = 8.064e-3 s`)
        unit: 표시 단위 문자열 (예: `"ms"`)
        si_multiplier: SI 기본단위 → 표시단위 배수 (초→ms 면 `1e3`)
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
# 필수 진단 세트 — 조화 트랩
# =============================================================================
def plot_msd(fs: FigureSet, runs: dict, *, dim: int, dt_star: float,
             tau_trap_si: float, l_trap_si: float, name: str = "01_msd") -> FigureRecord:
    """MSD log–log + 해석해 + 자유확산 참조선.

    자유확산 참조선(`2dD₀t`)이 있어야 **단시간 극한이 맞는지**가 보인다 — plateau 만
    맞고 단시간이 틀린 경우를 눈으로 구별할 수 있다.
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
        caption=(f"조화 트랩 MSD ({dim}D, 시드 {len(runs)}개). 단시간에 자유확산 "
                 f"기울기 1 을 따르고 `t ≳ τ_trap` 에서 `2d = {2*dim}` 로 포화한다. "
                 f"위/오른쪽 축은 물리 단위 (`τ_trap`·`ℓ_trap²` 환산)."),
        shows=("MSD 가 단일지수 `2d(1−e^{−t/τ})` 를 따르는가, 그리고 단시간 극한이 "
               "자유확산과 일치하는가 — plateau 만 맞고 단시간이 틀린 경우를 배제한다"))


def plot_equipartition_vs_dt(fs: FigureSet, by_dt: dict[float, tuple[float, float]],
                             *, name: str = "02_equipartition_dt") -> FigureRecord:
    """`⟨x*²⟩` vs `dt*` + Euler–Maruyama 편향 곡선 + exact 스킴 참조선.

    ★ **경쟁 가설(exact)을 함께 그린다.** EM 곡선만 그리면 "맞는 것처럼 보이는" 그림이
      되고, 두 곡선의 간격이 통계오차보다 작다는 사실(= 판별 불가)이 보이지 않는다.
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
        caption=(f"등분배와 적분기 편향. 측정점이 EM 곡선(검정)과 exact(빨강) 사이 어디에 "
                 f"있는지가 스킴을 판별한다. `dt* = {worst:g}` 에서 두 가설의 간격은 "
                 f"`{gap:.3%}` 이고 통계오차는 `±{se:.3%}` — 검정력 "
                 f"`{gap/se:.2f}σ`."),
        shows=("측정된 `⟨x*²⟩` 가 Euler–Maruyama 편향을 재현하는가, 그리고 exact 스킴과 "
               "**구별되는가** — 두 곡선의 간격 대비 오차막대 크기가 곧 검정력이다"))


def plot_position_distribution(fs: FigureSet, traj: np.ndarray, *, dt_star: float,
                               dim: int, l_trap_si: float,
                               frame_interval_steps: int,
                               decorrelation_tau: float = 2.0,
                               name: str = "03_distribution") -> FigureRecord:
    """평형 위치 분포 vs Gaussian (log y).

    ★ log y 로 그린다 — 꼬리에서 벗어나는 것이 선형 축에서는 안 보인다.
      HOOMD 노이즈가 균일분포라서 꼬리가 관심사다 (findings §2).

    ⚠ **독립 프레임만 쓴다.** `analysis.trap.check_position_distribution` 과 같은
      규율이다. 상관된 프레임을 전부 쓰면 표본 수가 부풀고, 첨도의 오차막대가
      실제보다 작게 나오며, `metrics.json` 의 검증된 값과 다른 숫자가 캡션에 박힌다.
      2026-07-28 에 KS 검정에서 같은 함정을 겪었다
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
        caption=(f"평형 위치 분포 (측정 표준편차로 규격화 — **폭이 아니라 형태**를 본다. "
                 f"폭은 `⟨x²⟩` 가 따로 검사한다). 첨도 `{kurt:.4f} ± {kurt_se:.4f}`, "
                 f"예측 `3 − 1.2 dt* = {pred:.4f}` → `{sigma_away:.2f}σ`. "
                 f"★ **정확히 3.000 이면 오히려 이상하다** — HOOMD 노이즈가 균일분포이므로 "
                 f"`dt*` 차수의 잔여 비가우시안성이 남는다. "
                 f"⚠ `{decorrelation_tau:g} τ_trap` 이상 떨어진 **독립 프레임 "
                 f"{indep.shape[0]}개**만 썼다 — 상관 프레임을 쓰면 표본 수가 부풀어 "
                 f"오차막대가 거짓으로 작아진다."),
        shows=("균일 노이즈가 CLT 로 Gaussian 에 접근하는가, 그리고 **알려진 잔여 편차를 "
               "재현하는가** — 무작정 Gaussian 을 확인하는 것보다 강한 검증이다"))


def plot_stationarity(fs: FigureSet, runs: dict, *, dt_star: float,
                      sample_interval_steps: int, tau_trap_si: float,
                      name: str = "04_stationarity") -> FigureRecord:
    """`⟨x*²⟩`·`kT_conf` 시계열 — 정상성(평형 도달) 확인.

    thermo 시계열의 자리다. HOOMD `Brownian` 의 운동에너지 온도는 계통 이탈이 불가능해
    쓸 수 없으므로 (findings §1), **배위 온도**와 등분배 통계량을 본다.
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
        # ★ 두 통계량을 겹쳐 그린다 — 대수적으로 같다면 완전히 겹쳐야 한다
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
        caption=(f"**위** 독립 표본 시점별 `⟨x*²⟩`(색)과 배위 온도 `kT_conf`(검정 ×). "
                 f"추세 없이 요동하면 평형에 도달했다는 뜻이다. EM 편향 `{bias:.3%}` 는 "
                 f"시점별 요동(`±5 %` 규모)보다 훨씬 작아 여기서는 보이지 않는다 — "
                 f"그래서 시드 앙상블 평균이 필요하다.\n\n"
                 f"**아래** 두 통계량의 상대 잔차. 최대 `{resid_max:.1e}` — 부동소수점 "
                 f"수준이다. ⚠ **순수 조화 트랩에서 `kT_conf` 는 `⟨x²⟩` 검사와 "
                 f"대수적으로 동일하다** (트랩 카드 §6 이 예고). 이 그림은 그것을 "
                 f"주장하지 않고 **보여준다** — 독립 검사가 되려면 쌍 상호작용이 필요하다."),
        shows=("① 표본이 정상(stationary)인가 — 추세가 있으면 평형화가 부족하다 "
               "② 통계량이 실제로 요동하는가 (요동하지 않으면 산술 항등식이다) "
               "③ `kT_conf` 가 `⟨x*²⟩` 와 독립인가 — 잔차가 0 이면 **독립이 아니다**"))


def plot_frame_displacements(fs: FigureSet, traj: np.ndarray, *, dt_star: float,
                             frame_interval_steps: int, sigma_star: float,
                             name: str = "05_displacements") -> FigureRecord:
    """프레임 간 변위 분포 — 가드 사후 확인.

    ⚠ **이것은 스텝당 변위가 아니다.** 궤적은 `frame_interval_steps` 스텝마다
      저장되므로, 균일 노이즈의 구조적 상한 `max/σ = √3` 은 여기 적용되지 않는다
      (여러 스텝의 합은 CLT 로 Gaussian 에 가까워진다). 스텝당 검사는 런타임
      가드(`simbot.guards.check_step_displacements`)가 한다.
    """
    tr = np.asarray(traj, dtype=np.float64)
    d = np.linalg.norm(tr[1:] - tr[:-1], axis=2).reshape(-1)
    n_steps = frame_interval_steps
    expected_rms = np.sqrt(2 * tr.shape[2] * n_steps * dt_star)   # 자유확산 상한

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
        caption=(f"프레임 간 변위 분포 (`{n_steps}` 스텝 간격). 최대값 "
                 f"`{d.max():.3f} ℓ_trap` = `{d.max()/sigma_star:.2e} σ` — 겹침 문턱 "
                 f"`0.1 σ` 에서 {0.1/(d.max()/sigma_star):.0f} 배 이하다. "
                 f"⚠ **스텝당 변위가 아니다.** 균일 노이즈의 `max/σ = √3` 상한은 "
                 f"여기 적용되지 않는다 (여러 스텝의 합)."),
        shows=("변위가 폭발 징후를 보이는가, 그리고 자유확산 상한 안에 있는가 — "
               "구속 때문에 자유확산보다 작아야 한다"))


def plot_snapshots(fs: FigureSet, traj: np.ndarray, *, dim: int, l_trap_si: float,
                   dt_star: float, frame_interval_steps: int, tau_trap_si: float,
                   name: str = "06_snapshots") -> FigureRecord:
    """초기·중간·최종 스냅샷 3연 + `ℓ_trap` 원.

    2D 만 그린다. 3D 는 `fresnel` 이 필요하고 아직 설치하지 않았다.
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
        # 시간도 이중 표기 — 스냅샷에는 보조축을 붙일 수 없으므로 제목에 넣는다
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
        caption=(f"초기·중간·최종 스냅샷. 원은 `1, 2, 3 ℓ_trap` "
                 f"(`ℓ_trap = {l_trap_si*1e9:.2f} nm`). 세 시점의 분포가 눈으로 "
                 f"구별되지 않으면 평형에 도달했다는 뜻이다. 입자 "
                 f"{tr.shape[1]}개는 같은 트랩 안의 **독립 복제**이지 서스펜션이 아니다."),
        shows=("입자들이 트랩 중심 주위에 등방으로 분포하는가, 그리고 세 시점이 "
               "구별되지 않는가 (= 정상 상태)"))


# =============================================================================
# 필수 진단 세트 — 2D 소프트 반발계 (시간분해)
# =============================================================================
#  카드: knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md
#  ★ 시간 축은 `τ_d = d²/D₀` 이고 보조축은 **분** 이다 (`τ_d ≈ 38 분` 규모).
#    ms 로 적으면 자리수가 6개 남아 읽을 수 없다.

#  배위수 색 — 5(빨강)/6(회색)/7(파랑) 이 전위 쌍의 표준 배색이다.
#  ★ 5 와 7 을 **대비되는 색**으로 두는 것이 요점이다: 전위는 5-7 이 동수로
#    나타나므로, 두 색이 같은 개수로 보이는지가 눈으로 판정할 수 있는 서명이 된다.
COORD_COLORS = {3: "#4a148c", 4: "#00838f", 5: "#d32f2f", 6: "#cfcfcf",
                7: "#1565c0", 8: "#f57c00", 9: "#6d4c41", 10: "#2e7d32"}
COORD_OTHER = "#000000"


def _coord_color(z: int) -> str:
    return COORD_COLORS.get(int(z), COORD_OTHER)


def plot_structure_timeseries(fs: FigureSet, series_by_A: dict, *,
                              tau_d_si: float, energy_by_A: dict | None = None,
                              name: str = "01_structure_timeseries"
                              ) -> FigureRecord:
    """`ψ₆` · 결함 분율 · 배위수 개체수 · 5-7 불균형의 **시간 궤적**.

    Args:
        series_by_A: `{A: [HexOrderSeries, ...시드별]}`. 시드 앙상블이 오차막대다 —
            **프레임 간 산포는 시간 상관을 포함하므로 통계오차가 아니다.**
        energy_by_A: `{A: (t_star, [energy_pp 배열, ...시드별])}` (선택)

    ★ 이 그림이 답하는 것은 "무엇이 되었는가" 가 아니라 **"언제 되었는가"** 다.
    """
    As = sorted(series_by_A)
    n_panel = 4 if energy_by_A is None else 5
    fig, axes = plt.subplots(n_panel, 1, figsize=(7.6, 2.05 * n_panel),
                             sharex=True)
    cmap = {A: f"C{i}" for i, A in enumerate(As)}

    def band(ax, t, mat, color, label):
        """시드 평균 ± SE. `mat` 은 `(n_seeds, n_frames)`."""
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
        caption=(f"위에서: 전역 `ψ₆` · Voronoi 결함 분율 · 배위수 종류 수 · "
                 f"5-7 불균형"
                 + (" · 입자당 에너지" if energy_by_A is not None else "")
                 + f". 띠는 **시드 앙상블 SE** (프레임 간 산포가 아니다 — 그것은 "
                 f"시간 상관을 포함한다). 마지막 프레임의 `ψ₆`: {fin}. "
                 f"보조축은 실제 시간이며 `τ_d = {tau_d_si/60:.1f} 분` 이다. "
                 f"**5-7 불균형이 0 으로 내려가면 결함이 전위 쌍으로 정리된 것이고**, "
                 f"배위수 종류 수가 6→3 으로 줄면 액체형 결함이 사라진 것이다."),
        shows=("구조가 언제 만들어지는가 — 최종값이 아니라 도달 시각과 표류 방향. "
               "평탄해지지 않은 곡선은 이 런 길이로 평형화를 주장할 수 없다는 증거다"))


def plot_rdf_evolution(fs: FigureSet, windows_by_A: dict, *, tau_d_si: float,
                       sigma_over_d: float | None = None,
                       name: str = "02_rdf_evolution") -> FigureRecord:
    """시간창별 `g(r)` — `A` 마다 패널 1개, 창마다 곡선 1개.

    Args:
        sigma_over_d: 주면 `r = σ` 에 수직선을 긋는다. **기준 원판이 겹치는지**
            눈으로 보이게 하는 선이다 (경질 코어가 아니므로 겹칠 수 있다).
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
        overlap_note = (f" 붉은 파선은 `r = σ = {sigma_over_d:.3f} d` — 기준 원판이 "
                        f"접하는 거리다. **이 안쪽에 `g(r) > 0` 이면 5 µm 원판이 "
                        f"겹친다** (`A/r³` 에는 경질 코어가 없으므로 모델은 허용한다).")
    fig.suptitle("RDF vs time window", fontsize=10)
    fig.tight_layout()
    return fs.save(
        fig, name,
        caption=(f"시간창별 `g(r)`. 창이 진해질수록 늦은 시각이다 "
                 f"(`τ_d = {tau_d_si/60:.1f} 분`). 창끼리 겹치면 그 구간에서 구조가 "
                 f"더 변하지 않는 것이고, 계속 갈라지면 아직 진행 중이다."
                 + overlap_note),
        shows="g(r) 의 첫 봉이 언제 자라고 언제 멈추는가 (그리고 멈추는가)")


def plot_early_transient(fs: FigureSet, series_by_A: dict, *, tau_d_si: float,
                         fits: dict | None = None,
                         name: str = "06_early_transient") -> FigureRecord:
    """초기 과도구간 — **로그 시간축.** 굵은 표집으로는 보이지 않는 구간이다.

    Args:
        fits: `{A: RelaxationFit}` (선택). `converged` 인 것만 곡선을 겹쳐 그리고,
            거부된 것은 범례에 **거부 사유**를 적는다 — 없는 완화를 그리지 않는다.

    ★ `t = 0` 은 로그축에 올릴 수 없으므로 **첫 프레임 앞의 초기배치는 왼쪽 여백에
      별표로** 찍는다. 초기배치를 빼면 "무엇으로부터 완화했는가" 가 사라진다.
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
            #  t=0 (초기배치) 은 로그축 밖 — 왼쪽 끝에 별표로
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
                         f"(`{f.tau*tau_d_si:.0f} s`), 폭 `{f.amplitude:+.3f}`")
        else:
            lines.append(f"`A={A:g}`: **완화 미검출** — {f.note}")
    return fs.save(
        fig, name,
        caption=("초기 과도구간을 **로그 시간축**으로. 별표는 초기배치(`t = 0`, "
                 "기각표집으로 `min_sep = 0.8 d` 를 강제한 배치) 이며 **세 `A` 의 "
                 "별표는 정확히 겹친다** — 초기배치가 `A` 에 의존하지 않으므로 "
                 "(같은 시드 → 같은 배치) 세 곡선은 **같은 점에서 갈라진다.** "
                 "곡선은 "
                 "시드 앙상블 평균 ± SE 다. 파선은 단일지수 적합이며 **완화 폭이 "
                 "잡음의 2배를 넘지 못하면 그리지 않는다** — 정상상태 요동에 지수를 "
                 "맞추면 항상 어떤 `τ` 가 나오기 때문이다. "
                 + (" · ".join(lines) if lines else "")),
        shows=("초기배치에서 정상상태로 가는 완화가 실제로 존재하는가, 그리고 "
               "존재하면 그 시간이 얼마인가 — 굵은 표집(0.2 τ_d)으로는 첫 프레임에 "
               "이미 끝나 보인다"))


def plot_seed_convergence(fs: FigureSet, conv: dict, *, curves: dict,
                          tau_d_si: float, threshold_sigma: float = 3.0,
                          preregistered: list | None = None,
                          name: str = "01_seed_convergence") -> FigureRecord:
    """`τ(k)` 수렴 · `σ(k)` · 완화 곡선 3연.

    Args:
        conv: `{"k": [...], "A0.1": {"tau": [...], "se": [...]}, "A1": {...},
                "diff": [...], "se_diff": [...], "sigma": [...]}`
            — **같은 데이터의 중첩 부분집합**이어야 한다. 그러면 곡선의 기울기가
            데이터 차이가 아니라 **추정량의 `k` 의존성**만 나타낸다.
        curves: `{A: (t_star, mean_defect, fit)}` 최종 `k` 의 완화 곡선과 적합.
        preregistered: `[{"k":…, "sigma_expected":…, "stage":…}]` 사전등록 예상치.

    ★ 이 그림의 요점은 **예상과 실제의 어긋남**이다. 사전등록 예상치를 함께 찍어야
      "저시드 예비런이 낙관적이었다"가 눈에 보인다.
    """
    k = np.asarray(conv["k"], dtype=float)
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 3.7))

    # --- ① tau vs k ---------------------------------------------------------
    #  ★ 오차막대는 **강건 구간**(부트스트랩 백분위)이다. SD 는 저시드에서 이상치가
    #    지배해서 (k=4 에서 82 τ_d — τ 자체는 0.06) 축을 파괴하고 그림을 못 읽게 한다.
    #    SD 가 τ 보다 큰 지점은 **속 빈 표식**으로 "오차를 주장할 수 없다"를 표시한다.
    ax = axes[0]
    tau_all = []
    for i, A in enumerate(("A0.1", "A1")):
        tau = np.asarray(conv[A]["tau"], dtype=float)
        se = np.asarray(conv[A]["se"], dtype=float)
        se_r = np.asarray(conv[A].get("se_robust", se), dtype=float)
        usable = se <= np.abs(tau)          # SD 가 τ 보다 크면 무의미하다
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
    ax.set_ylim(lo - pad, hi + pad)         # ★ τ 값이 정하는 축 — 오차막대가 아니다
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
        #  실제와 예상을 화살표로 잇는다 — 어긋남이 이 그림의 요점이다
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

    # --- ③ 완화 곡선 --------------------------------------------------------
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
        caption=(f"**① `τ` 추정치가 시드 수와 함께 움직인다.** 이 곡선들은 "
                 f"**같은 궤적의 중첩 부분집합**에서 나왔으므로 기울기는 데이터 차이가 "
                 f"아니라 **추정량의 `k` 의존성**이다 — 즉 저시드 편향이다. 오차막대는 "
                 f"시드 부트스트랩 (`curve_fit` 공분산이 아니다). "
                 f"**② 사전등록 예상(별)과 실제(선)의 어긋남.** "
                 f"**③ 최종 `k = {k_final}` 의 결함 완화 곡선과 단일지수 적합**; "
                 f"별표는 초기배치(`t=0`)이며 두 `A` 가 같은 점에서 갈라진다. "
                 f"`τ_d = {tau_d_si/60:.1f} 분`."),
        shows=("저시드 예비런이 왜 낙관적인 검정력을 주는가 — 편향된 차이와 "
               "과소추정된 오차가 곱해진다. 그리고 최종 표본에서 τ 가 갈리는가"))


def plot_finite_size_scaling(fs: FigureSet, data: dict, *, local_data: dict,
                             exponent_liquid: float = 0.5,
                             exponent_hexatic: float = 0.0625,
                             name: str = "01_finite_size_scaling") -> FigureRecord:
    """`ψ₆` vs `N` (log-log) + 국소량의 `N` 무관성 2연.

    Args:
        data: `{A: {"N": [...], "psi6": [...], "se": [...], "fit": FiniteSizeExponent}}`
        local_data: `{A: {"N": [...], "defect": [...], "defect_se": [...]}}`

    ★ 이 그림의 요점은 **기울기**다. 절대값이 아니다 — 유한계에서는 무질서한 계도
      `1/√N` 규모의 `ψ₆` 를 낸다. 액체 기울기(`−1/2`)와 KTHNY hexatic 경계
      기울기(`−1/16`)를 안내선으로 깔아 눈으로 비교되게 한다.
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
            #  적합 기울기
            nn = np.array([N.min() * 0.9, N.max() * 1.1])
            ax.plot(nn, y[0] * (nn / N[0]) ** (-f.p), color=c, lw=1.5, alpha=0.9)
        ax.plot([], [], marker="o", ls="-", color=c, label=lab)
        #  두 가설 안내선 — 첫 점에 고정
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
                         + f" → `η₆ = {f.eta6:.2f}`, 판독 **{f.reading}**")
    return fs.save(
        fig, name,
        caption=("**① `ψ₆` 의 `N` 의존성.** 파선은 액체 기울기 "
                 f"(`p = {exponent_liquid:g}`, `|⟨ψ₆⟩| ~ N^-1/2`), 점선은 KTHNY "
                 f"hexatic-액체 경계 (`p = {exponent_hexatic:g}`, `η₆ = 1/4`) 이며 "
                 f"각 `A` 의 첫 점에 고정했다. **절대값이 아니라 기울기가 상을 "
                 f"판별한다** — 유한계에서는 무질서한 계도 `1/√N` 규모의 `ψ₆` 를 "
                 f"낸다. " + " · ".join(lines) +
                 ". **② 결함 분율은 국소량이므로 `N` 에 평탄해야 한다** — "
                 "기울어지면 `N=100` 에 유한크기 효과가 있었다는 뜻이다. "
                 "⚠ 두 점으로는 멱함수 **형태**를 검증할 수 없다; 지수만 뽑았다."),
        shows=("ψ₆ 의 절대값이 유한크기 바닥인가 진짜 준장거리 배향 질서인가 — "
               "그리고 국소량이 N=100 에서 이미 수렴했는가"))


def plot_voronoi_timelapse(fs: FigureSet, frames: list, *, t_star: list,
                           tau_d_si: float, L_si: float, amplitude: float,
                           mean_defect_fraction: float | None = None,
                           defect_frame_sd: float | None = None,
                           name: str = "03_voronoi_timelapse") -> FigureRecord:
    """Voronoi 타일링 시간열 — 셀을 **배위수로** 칠한다.

    카드 §10 의 미구현 항목("`voronoi plot` 이 그림에 요청됐으나 `viz.py` 에 없다").

    ★ 배위수 5(빨강)/7(파랑) 이 **쌍으로** 보이면 전위이고, 5·7·4·8 이 뒤섞여
      흩어져 있으면 액체다. 결함 **분율** 하나로는 이 구별이 사라진다 (카드 §8.2).
    """
    if len(frames) != len(t_star):
        raise ValueError(f"프레임 {len(frames)}개 vs 시각 {len(t_star)}개 — "
                         f"개수가 같아야 한다")
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

    #  ★ 단일 프레임은 요동한다 — 패널 사이의 차이를 완화로 읽으면 안 된다.
    caveat = ""
    if mean_defect_fraction is not None and defect_frame_sd is not None:
        n = frames[0].coordination.size
        caveat = (f" ⚠ **패널의 결함 수 차이를 시간 변화로 읽지 말 것** — 이 계의 "
                  f"프레임별 결함 수는 정상상태에서 "
                  f"`{mean_defect_fraction*n:.0f} ± {defect_frame_sd*n:.0f}` 로 "
                  f"요동한다 (프레임 표준편차). 패널 간 차이는 대부분 그 요동이다. "
                  f"시간 변화는 시계열 그림에서 앙상블 평균으로 판단해야 한다.")
    return fs.save(
        fig, name,
        caption=(f"`A = {amplitude:g}` 의 Voronoi 타일링 시간열 (**시드 1개**). 색은 "
                 f"배위수 `z` — **빨강 `z=5`, 파랑 `z=7`, 회색 `z=6`**. 상자는 "
                 f"`{L_si*1e6:.0f} µm` 정사각, 주기경계다 (셀이 상자 밖으로 나가는 "
                 f"것은 주기 이웃 쪽으로 열린 셀이며 잘라내지 않았다 — 자르면 면적이 "
                 f"틀린다). 빨강·파랑이 **붙어서 쌍으로** 보이면 전위이고, 여러 색이 "
                 f"흩어져 있으면 액체다." + caveat),
        shows=("결함의 **성격** — 분율 숫자로는 전위와 액체가 구별되지 않는다. "
               "이 그림은 성격을 보이는 것이고, 시간 변화의 근거는 시계열 그림이다"))


# =============================================================================
# 오케스트레이터
# =============================================================================
def trap_diagnostics(outdir: Path, runs: dict, *, dim: int, dt_star: float,
                     sample_interval_steps: int, frame_interval_steps: int,
                     sigma_star: float, tau_trap_si: float, l_trap_si: float,
                     by_dt: dict[float, tuple[float, float]] | None = None,
                     has_pair: bool = False) -> FigureSet:
    """조화 트랩 계의 필수 진단 세트.

    Args:
        runs: `label → samples.npz 내용`. 시드 앙상블 전체.
        by_dt: `dt* → (⟨x*²⟩ 평균, SE)`. 여러 `dt*` 를 돌렸을 때만. 없으면
            그 그림을 **이유와 함께 건너뛴다.**
        has_pair: 쌍 상호작용 유무. 없으면 RDF 를 이유와 함께 건너뛴다.
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
                f"`dt*` 가 {len(by_dt or {})}개뿐이다 — 편향의 `dt*` 의존성을 그릴 수 "
                f"없다. `cli.py converge` 로 dt 래더를 돌리면 생긴다")

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
                "3D 스냅샷은 `fresnel` 레이트레이싱이 필요하고 아직 설치하지 않았다 "
                "(techniques/env-log.md). 2D 투영은 오해를 만든다")

    if has_pair:
        fs.skip("07_rdf.png",
                "쌍 상호작용이 있으나 RDF 계산기(`analysis/structure.py`)가 아직 없다 — "
                "트랩+WCA 계에서 함께 만든다")
    else:
        fs.skip("07_rdf.png",
                "★ 쌍 상호작용이 없어 `g(r)` 이 정의되지 않는다. 입자들은 같은 트랩 안의 "
                "독립 복제이므로 상대 거리 분포는 두 Gaussian 의 차이일 뿐이고 "
                "구조 정보가 없다")
    return fs
