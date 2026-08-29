"""S7 — 2D 구조 분석. 측정만 하고 판정하지 않는다.

`freud` 3.5 를 쓴다. RDF · 육방 order parameter `ψ₆` · Voronoi 결함 · `S(q)`.

## 이 모듈이 조심하는 것

**① 길이 단위는 `d = n^{-1/2}` 다** (최근접거리 `a = 1.0746 d` 가 아니다).
`Γ ∝ d^{-3}` 이므로 7 % 혼동이 `Γ` 를 23 % 틀리게 만든다.

**② `ψ₆` 는 이웃 정의에 의존한다.** Voronoi 이웃과 "가장 가까운 6개"는 다른 답을
준다 — 결함 근처에서 특히. **어느 것을 썼는지 함께 보고한다.**

**③ 프레임을 독립으로 취급하지 않는다.** 구조 완화는 느리다. 시드 앙상블과
프레임 간격을 호출자가 정하고, 이 모듈은 `n_frames` 를 그대로 보고한다.
근거: [[ks-test-needs-independent-samples]]
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _box(Lx: float, Ly: float):
    import freud
    return freud.box.Box(Lx=Lx, Ly=Ly, is2D=True)


def _xyz(pos_2d: np.ndarray) -> np.ndarray:
    """freud 는 2D 에서도 (N, 3) 을 요구한다 (z = 0)."""
    p = np.asarray(pos_2d, dtype=np.float64)
    if p.shape[1] == 3:
        return p
    out = np.zeros((p.shape[0], 3))
    out[:, :2] = p
    return out


# =============================================================================
# RDF
# =============================================================================
@dataclass
class RDFResult:
    r: np.ndarray
    g: np.ndarray
    first_peak_r: float
    first_peak_g: float
    n_frames: int

    def coordination(self, r_max: float, density_star: float = 1.0) -> float:
        """`r < r_max` 안의 평균 이웃 수 `n(r) = 2πn ∫ g(r) r dr` (2D)."""
        m = self.r <= r_max
        return float(2 * np.pi * density_star
                     * np.trapezoid(self.g[m] * self.r[m], self.r[m]))


def rdf(frames: np.ndarray, *, Lx: float, Ly: float, r_max: float | None = None,
        bins: int = 200) -> RDFResult:
    """여러 프레임 평균 `g(r)`. `frames` 는 `(n_frames, N, 2)`."""
    import freud
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    r_max = r_max if r_max is not None else min(Lx, Ly) / 2 * 0.99
    box = _box(Lx, Ly)
    calc = freud.density.RDF(bins=bins, r_max=r_max)
    for f in fr:
        calc.compute(system=(box, _xyz(f)), reset=False)
    g, r = np.asarray(calc.rdf), np.asarray(calc.bin_centers)
    # 첫 봉우리 — r > 0.5 d 에서 찾는다 (그 안쪽은 표본이 거의 없다)
    m = r > 0.5
    i = int(np.argmax(g[m]))
    return RDFResult(r=r, g=g, first_peak_r=float(r[m][i]),
                     first_peak_g=float(g[m][i]), n_frames=int(fr.shape[0]))


# =============================================================================
# 육방 order parameter + Voronoi 결함
# =============================================================================
@dataclass
class HexOrderResult:
    psi6_global: float              # |<psi6>| — 전역 배향 질서
    psi6_local_mean: float          # <|psi6_i|> — 국소 질서 (전역보다 항상 크거나 같다)
    psi6_per_frame: list[float]
    defect_fraction: float          # Voronoi 배위수 != 6 인 입자 분율
    coordination_hist: dict[int, float]
    neighbor_mode: str
    n_frames: int
    n_particles: int

    @property
    def psi6_spread(self) -> float:
        a = np.asarray(self.psi6_per_frame)
        return float(a.std(ddof=1)) if a.size > 1 else float("nan")


def _frame_hex(f: np.ndarray, box, hexatic, voro, neighbor_mode: str
               ) -> tuple[np.ndarray, np.ndarray]:
    """프레임 1개의 `(ψ₆ᵢ, 배위수ᵢ)`. `hex_order` 와 `hex_order_series` 가 공유한다.

    ★ 공유하는 이유: 시간분해 판본이 따로 계산하면 두 함수가 조용히 갈라진다.
      "후반 절반 평균" 과 "프레임별 시계열의 후반 평균" 은 **같은 수여야 한다.**
    """
    pts = _xyz(f)
    system = (box, pts)
    if neighbor_mode == "voronoi":
        voro.compute(system)
        hexatic.compute(system, neighbors=voro.nlist)
    elif neighbor_mode == "nearest6":
        hexatic.compute(system, neighbors={"num_neighbors": 6})
        voro.compute(system)
    else:
        raise ValueError(f"neighbor_mode {neighbor_mode!r} 는 "
                         f"'voronoi' 또는 'nearest6' 여야 한다")
    counts = np.bincount(np.asarray(voro.nlist.query_point_indices),
                         minlength=pts.shape[0])
    return np.asarray(hexatic.particle_order), counts


def hex_order(frames: np.ndarray, *, Lx: float, Ly: float,
              neighbor_mode: str = "voronoi") -> HexOrderResult:
    """`ψ₆` 와 Voronoi 결함 분율.

    Args:
        neighbor_mode: `"voronoi"` (권장) 또는 `"nearest6"`.
            ★ 둘은 다른 답을 준다 — 결함 근처에서 Voronoi 는 5·7 이웃을 그대로
            세지만 `nearest6` 은 강제로 6개를 만든다. **결함을 보려면 Voronoi 다.**

    `psi6_global = |⟨ψ₆ᵢ⟩|` 와 `psi6_local_mean = ⟨|ψ₆ᵢ|⟩` 를 **둘 다** 준다:
      · 국소만 크면 결정립이 있으나 방향이 갈린다 (다결정)
      · 둘 다 크면 단결정
      · 국소만 크고 전역이 0 이면 hexatic 후보
    """
    import freud
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    box = _box(Lx, Ly)
    hexatic = freud.order.Hexatic(k=6)
    voro = freud.locality.Voronoi()

    per_frame, local_means, coord_counts = [], [], []
    for f in fr:
        psi, counts = _frame_hex(f, box, hexatic, voro, neighbor_mode)
        per_frame.append(float(np.abs(psi.mean())))
        local_means.append(float(np.abs(psi).mean()))
        coord_counts.append(counts)

    all_counts = np.concatenate(coord_counts)
    uniq, cnt = np.unique(all_counts, return_counts=True)
    hist = {int(u): float(c / all_counts.size) for u, c in zip(uniq, cnt)}
    return HexOrderResult(
        psi6_global=float(np.mean(per_frame)),
        psi6_local_mean=float(np.mean(local_means)),
        psi6_per_frame=[float(x) for x in per_frame],
        defect_fraction=float(np.mean(all_counts != 6)),
        coordination_hist=hist, neighbor_mode=neighbor_mode,
        n_frames=int(fr.shape[0]), n_particles=int(fr.shape[1]))


# =============================================================================
# 시간분해 — 구조가 **언제** 만들어지는가
# =============================================================================
#  ★ 왜 시간 평균과 따로 두는가
#  후반 절반 평균은 "무엇이 되었는가" 에 답하고 **"언제 되었는가" 에는 답하지 않는다.**
#  같은 최종 `ψ₆` 가 (ⅰ) 1 τ_d 에 도달해 머문 것인지 (ⅱ) 아직 오르는 중인 것인지
#  구별되지 않는다 — 그 구별이 평형화 진단의 전부다.
#
#  ⚠ **프레임은 독립 표본이 아니다** (모듈 §③). 시계열의 프레임 간 산포를
#    통계오차로 쓰면 시간 상관 때문에 과소평가된다. 오차막대는 시드 앙상블에서 온다.
@dataclass
class HexOrderSeries:
    """프레임별 `ψ₆`·결함·배위수. 시간 축은 호출자가 준 `t_star`."""

    t_star: np.ndarray                  # (n_frames,)
    psi6_global: np.ndarray             # (n_frames,) |⟨ψ₆ᵢ⟩| — 프레임 내 평균의 크기
    psi6_local: np.ndarray              # (n_frames,) ⟨|ψ₆ᵢ|⟩
    defect_fraction: np.ndarray         # (n_frames,) 배위수 ≠ 6 분율
    coord_labels: np.ndarray            # (n_coord,) 배위수 값 (예: [3,4,…,10])
    coord_fraction: np.ndarray          # (n_frames, n_coord) 프레임별 분율
    neighbor_mode: str
    n_particles: int

    @property
    def n_frames(self) -> int:
        return int(self.t_star.size)

    def _col(self, z: int) -> np.ndarray:
        j = np.nonzero(self.coord_labels == z)[0]
        return self.coord_fraction[:, j[0]] if j.size else np.zeros(self.n_frames)

    @property
    def five_seven_balance(self) -> np.ndarray:
        """`|n₅ − n₇|/(n₅ + n₇)` 프레임별. **0 이면 전위(5-7 쌍) 서명이다.**

        결함 분율만으로는 전위와 액체가 구별되지 않는다 (카드 §8.2) — 이 양이
        그 구별을 시간축 위에서 유지한다.
        """
        n5, n7 = self._col(5), self._col(7)
        return np.abs(n5 - n7) / np.maximum(n5 + n7, 1e-12)

    @property
    def coord_kinds(self) -> np.ndarray:
        """프레임별 '분율 > 0.5 %' 인 배위수 종류 수. 액체는 6종, 결정은 3종."""
        return (self.coord_fraction > 0.005).sum(axis=1)

    def window_mean(self, t_lo: float, t_hi: float) -> dict:
        """`t_lo ≤ t < t_hi` 구간 평균. **후반 절반 평균과 대조하는 데 쓴다.**"""
        m = (self.t_star >= t_lo) & (self.t_star < t_hi)
        if not m.any():
            raise ValueError(f"[{t_lo}, {t_hi}) 안에 프레임이 없다")
        return {"n_frames": int(m.sum()),
                "psi6_global": float(self.psi6_global[m].mean()),
                "psi6_local": float(self.psi6_local[m].mean()),
                "defect_fraction": float(self.defect_fraction[m].mean())}


def hex_order_series(frames: np.ndarray, *, Lx: float, Ly: float,
                     t_star: np.ndarray | None = None,
                     neighbor_mode: str = "voronoi",
                     coord_range: tuple[int, int] = (3, 10)) -> HexOrderSeries:
    """`hex_order` 의 **시간분해** 판본 — 프레임을 평균하지 않고 그대로 돌려준다.

    `hex_order(frames).psi6_global` 은 이 함수의 `psi6_global.mean()` 과
    **정확히 같다** (같은 `_frame_hex` 를 쓴다). 테스트가 그것을 감시한다.

    Args:
        coord_range: 배위수 히스토그램 범위 (양끝 포함). 밖의 값은 마지막 칸에
            누적하지 않고 **버린다** — 대신 `defect_fraction` 은 전량으로 센다.
    """
    import freud
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    n_frames, N = fr.shape[0], fr.shape[1]
    t = (np.arange(n_frames, dtype=np.float64) if t_star is None
         else np.asarray(t_star, dtype=np.float64))
    if t.size != n_frames:
        raise ValueError(f"t_star 길이 {t.size} != 프레임 수 {n_frames}")

    box = _box(Lx, Ly)
    hexatic = freud.order.Hexatic(k=6)
    voro = freud.locality.Voronoi()
    labels = np.arange(coord_range[0], coord_range[1] + 1)

    pg = np.empty(n_frames)
    pl = np.empty(n_frames)
    df = np.empty(n_frames)
    cf = np.zeros((n_frames, labels.size))
    for k, f in enumerate(fr):
        psi, counts = _frame_hex(f, box, hexatic, voro, neighbor_mode)
        pg[k] = abs(psi.mean())
        pl[k] = np.abs(psi).mean()
        df[k] = float(np.mean(counts != 6))
        for j, z in enumerate(labels):
            cf[k, j] = float(np.mean(counts == z))
    return HexOrderSeries(t_star=t, psi6_global=pg, psi6_local=pl,
                          defect_fraction=df, coord_labels=labels,
                          coord_fraction=cf, neighbor_mode=neighbor_mode,
                          n_particles=int(N))


@dataclass
class RelaxationFit:
    """`y(t) = y_inf + (y_0 − y_inf) exp(−t/τ)` 적합. **판정하지 않는다.**"""

    y0: float
    y_inf: float
    tau: float
    tau_se: float
    r_squared: float
    n_points: int
    converged: bool
    note: str = ""

    @property
    def amplitude(self) -> float:
        """완화 폭 `y_0 − y_inf`. 이것이 잡음보다 작으면 `τ` 는 의미가 없다."""
        return self.y0 - self.y_inf


def fit_relaxation(t_star: np.ndarray, y: np.ndarray, *,
                   noise: float | None = None) -> RelaxationFit:
    """단일지수 완화 시간.

    ★ **`τ` 를 보고하기 전에 `amplitude` 를 잡음과 비교해야 한다.** 정상상태를
      요동하는 신호에 지수를 맞추면 항상 어떤 `τ` 가 나온다 — 그것은 완화가 아니라
      잡음이다. `noise` 를 주면 `converged` 에 그 판정을 반영한다.

    Args:
        noise: 신호의 프레임 간 표준편차. `|amplitude| < 2·noise` 면
            `converged=False` 로 두고 이유를 `note` 에 남긴다.
    """
    from scipy.optimize import curve_fit

    t = np.asarray(t_star, dtype=np.float64)
    yy = np.asarray(y, dtype=np.float64)
    if t.size != yy.size:
        raise ValueError(f"t({t.size}) 와 y({yy.size}) 길이가 다르다")
    if t.size < 4:
        raise ValueError(f"점 {t.size}개로는 3파라미터 적합을 할 수 없다")

    def model(tt, y0, y_inf, tau):
        return y_inf + (y0 - y_inf) * np.exp(-tt / np.maximum(tau, 1e-12))

    span = float(t[-1] - t[0]) or 1.0
    p0 = (float(yy[0]), float(yy[-1]), span / 5.0)
    try:
        popt, pcov = curve_fit(model, t, yy, p0=p0, maxfev=40000,
                               bounds=([-np.inf, -np.inf, 1e-9],
                                       [np.inf, np.inf, 100.0 * span]))
        ok = True
    except Exception as e:                                   # noqa: BLE001
        return RelaxationFit(y0=float(yy[0]), y_inf=float(yy[-1]),
                             tau=float("nan"), tau_se=float("nan"),
                             r_squared=float("nan"), n_points=int(t.size),
                             converged=False, note=f"적합 실패: {e!r}")
    resid = yy - model(t, *popt)
    ss_tot = float(((yy - yy.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    tau_se = float(np.sqrt(pcov[2, 2])) if np.all(np.isfinite(pcov)) else float("nan")

    note = ""
    amp = float(popt[0] - popt[1])
    if noise is not None and abs(amp) < 2.0 * noise:
        ok = False
        note = (f"완화 폭 |{amp:.4g}| < 2×잡음({noise:.4g}) — 이것은 완화가 아니라 "
                f"정상상태 요동이다. τ 를 보고하지 말 것")
    if np.isfinite(tau_se) and tau_se > abs(popt[2]):
        ok = False
        note = note or (f"τ = {popt[2]:.4g} ± {tau_se:.4g} — 오차가 값보다 크다")
    return RelaxationFit(y0=float(popt[0]), y_inf=float(popt[1]),
                         tau=float(popt[2]), tau_se=tau_se, r_squared=r2,
                         n_points=int(t.size), converged=ok, note=note)


def bootstrap_relaxation_over_seeds(
        t_star: np.ndarray, per_seed: np.ndarray, *, n_resample: int = 400,
        seed: int = 0, noise: float | None = None) -> dict:
    """`τ` 의 **시드 앙상블** 오차. `curve_fit` 공분산과 다른 것을 잰다.

    ★ `fit_relaxation` 이 주는 `tau_se` 는 *시드 평균 곡선 하나*에 대한 적합
      불확실성이다. 이 저장소의 규약은 **"시드 앙상블이 오차의 정직한 추정치"**
      (`scripts/soft2d_sweep_analyze.py` §①) 이므로, 시드를 재표집해서 다시 적합한
      분포의 산포가 옳은 오차다.

      둘이 크게 다르면 **적합 불확실성이 시드 간 변동을 대표하지 못한다**는 뜻이고,
      그때 `curve_fit` SE 로 세운 σ 거리는 믿을 수 없다.

    Args:
        per_seed: `(n_seeds, n_frames)` — 시드별 곡선. 평균하지 말고 그대로 준다.
        n_resample: 복원추출 횟수.

    Returns: `tau` (전체 평균 곡선 적합) · `tau_se_bootstrap` · `tau_ci95` ·
        `n_converged` · `tau_se_fit` (비교용).
    """
    y = np.asarray(per_seed, dtype=np.float64)
    if y.ndim != 2:
        raise ValueError(f"per_seed 는 (n_seeds, n_frames) 여야 한다 — {y.shape}")
    k = y.shape[0]
    if k < 2:
        raise ValueError(f"시드 {k}개로는 앙상블 오차를 낼 수 없다")

    full = fit_relaxation(t_star, y.mean(axis=0), noise=noise)
    rng = np.random.default_rng(seed)
    taus: list[float] = []
    for _ in range(int(n_resample)):
        idx = rng.integers(0, k, size=k)               # 복원추출
        try:
            f = fit_relaxation(t_star, y[idx].mean(axis=0))
        except Exception:                              # noqa: BLE001
            continue
        if np.isfinite(f.tau):
            taus.append(f.tau)
    arr = np.asarray(taus)
    if arr.size < 10:
        return {"tau": full.tau, "tau_se_bootstrap": float("nan"),
                "tau_ci95": (float("nan"), float("nan")),
                "n_converged": int(arr.size), "tau_se_fit": full.tau_se,
                "note": f"부트스트랩 적합이 {arr.size}회만 수렴 — 오차를 주장할 수 없다"}
    return {"tau": full.tau,
            "tau_se_bootstrap": float(arr.std(ddof=1)),
            "tau_ci95": (float(np.percentile(arr, 2.5)),
                         float(np.percentile(arr, 97.5))),
            "n_converged": int(arr.size), "tau_se_fit": full.tau_se,
            "se_ratio_bootstrap_over_fit": (float(arr.std(ddof=1) / full.tau_se)
                                            if full.tau_se else float("nan")),
            "note": ""}


@dataclass
class RDFWindows:
    """시간창별 `g(r)`. **창 안에서는 평균하고 창끼리는 섞지 않는다.**"""

    r: np.ndarray                       # (bins,)
    g: np.ndarray                       # (n_windows, bins)
    t_lo: np.ndarray                    # (n_windows,)
    t_hi: np.ndarray
    frames_per_window: np.ndarray
    first_peak_r: np.ndarray            # (n_windows,)
    first_peak_g: np.ndarray


def rdf_windows(frames: np.ndarray, *, Lx: float, Ly: float,
                t_star: np.ndarray, n_windows: int = 4,
                r_max: float | None = None, bins: int = 200) -> RDFWindows:
    """궤적을 `n_windows` 개 시간창으로 나눠 각각의 `g(r)`.

    ★ 창 1개짜리 `g(r)` 은 프레임이 적어 시끄럽다 — 창 수를 늘리면 시간해상도가
      오르고 신호대잡음이 내린다. **그 교환을 호출자가 정한다.**
    """
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    t = np.asarray(t_star, dtype=np.float64)
    if t.size != fr.shape[0]:
        raise ValueError(f"t_star 길이 {t.size} != 프레임 수 {fr.shape[0]}")
    if n_windows < 1:
        raise ValueError("n_windows ≥ 1")

    edges = np.linspace(t[0], t[-1] + (t[-1] - t[0]) * 1e-9, n_windows + 1)
    r_ref, gs, los, his, cnts, pr, pg = None, [], [], [], [], [], []
    for w in range(n_windows):
        m = (t >= edges[w]) & (t < edges[w + 1])
        if not m.any():
            raise ValueError(
                f"창 {w} [{edges[w]:.4g}, {edges[w + 1]:.4g}) 에 프레임이 없다 — "
                f"n_windows({n_windows}) 가 프레임 수({t.size}) 에 비해 크다")
        res = rdf(fr[m], Lx=Lx, Ly=Ly, r_max=r_max, bins=bins)
        r_ref = res.r if r_ref is None else r_ref
        gs.append(res.g)
        los.append(edges[w]); his.append(edges[w + 1])
        cnts.append(int(m.sum()))
        pr.append(res.first_peak_r); pg.append(res.first_peak_g)
    return RDFWindows(r=r_ref, g=np.array(gs), t_lo=np.array(los),
                      t_hi=np.array(his), frames_per_window=np.array(cnts),
                      first_peak_r=np.array(pr), first_peak_g=np.array(pg))


@dataclass
class VoronoiFrame:
    """프레임 1개의 Voronoi 타일링 — **그림용**. 측정은 `hex_order_series` 가 한다."""

    positions: np.ndarray               # (N, 2)
    polygons: list                      # 길이 N, 각 (n_v, 2)
    coordination: np.ndarray            # (N,) 정수
    psi6_abs: np.ndarray                # (N,) |ψ₆ᵢ| — 국소 육방성
    Lx: float
    Ly: float


def voronoi_frame(frame: np.ndarray, *, Lx: float, Ly: float,
                  neighbor_mode: str = "voronoi") -> VoronoiFrame:
    """Voronoi 다각형 + 배위수 + 국소 `ψ₆`. 카드 §10 의 `voronoi plot` 재료.

    ⚠ 다각형 꼭짓점은 **상자 밖으로 나갈 수 있다** (주기 이웃 쪽으로 열린 셀).
      그림에서 잘라내는 것은 축 한계가 할 일이고, 여기서 자르면 셀 면적이 틀린다.
    """
    import freud
    f = np.asarray(frame, dtype=np.float64)
    if f.ndim != 2:
        raise ValueError(f"프레임 1개 (N, 2) 를 줘야 한다 — 받은 shape {f.shape}")
    box = _box(Lx, Ly)
    hexatic = freud.order.Hexatic(k=6)
    voro = freud.locality.Voronoi()
    psi, counts = _frame_hex(f, box, hexatic, voro, neighbor_mode)
    polys = [np.asarray(p)[:, :2] for p in voro.polytopes]
    return VoronoiFrame(positions=f[:, :2], polygons=polys,
                        coordination=counts.astype(int),
                        psi6_abs=np.abs(psi), Lx=Lx, Ly=Ly)


# =============================================================================
# S(q)
# =============================================================================
@dataclass
class StructureFactorResult:
    """2D `S(k)` — **벡터 맵과 방향 평균을 함께** 준다.

    ★ 방향 평균만 보면 육방 결정과 hexatic 을 구별할 수 없다. 둘 다 같은 `|k|` 에
      첫 봉우리가 선다. 갈리는 것은 **각도 의존성**이다:
        결정   → 6겹 **점** (Bragg)
        hexatic → 6겹으로 변조된 **고리**
        액체   → 등방 고리
      그래서 `S_map` 과 `sixfold_modulation` 이 필요하다.
    """

    kx: np.ndarray                  # (n_kx,) 격자와 정합한 k 성분
    ky: np.ndarray
    S_map: np.ndarray               # (n_kx, n_ky) 벡터 S(k)
    k_radial: np.ndarray            # 방향 평균용 |k| 빈 중심
    S_radial: np.ndarray
    first_peak_k: float
    first_peak_S: float
    sixfold_modulation: float       # 첫 봉 고리에서의 6겹 세기 (0 = 등방)
    k_min: float                    # 상자가 허용하는 최소 |k| = 2π/max(L)
    n_frames: int


def structure_factor(frames: np.ndarray, *, Lx: float, Ly: float,
                     n_max: int = 24, bins: int = 120,
                     ring_width: float = 0.15) -> StructureFactorResult:
    """2D 정적 구조인자를 **직접** 계산한다.

    `S(k) = |Σ_j exp(−i k·r_j)|² / N`, `k = 2π(m/Lx, n/Ly)` — 상자와 정합한
    k-벡터만 쓰므로 **표집 오차가 없다** (근사 없는 이산 합).

    ⚠ **`freud` 를 쓰지 않는다.** `freud.diffraction.StaticStructureFactorDirect`
      는 `ValueError: 2D boxes are not currently supported` 를 던진다 (freud 3.5.0,
      2026-07-28 실측). 2D 가 이 프로젝트의 1급 지원 대상이라 직접 구현했다.

    Args:
        n_max: k 격자 반경 (정수 지수). `|k|_max ≈ 2π n_max/min(L)`
        ring_width: 6겹 변조를 재는 고리의 상대 두께 (첫 봉 `|k|` 기준)
    """
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    n_frames, N = fr.shape[0], fr.shape[1]

    m = np.arange(-n_max, n_max + 1)
    kx = 2 * np.pi * m / Lx
    ky = 2 * np.pi * m / Ly
    KX, KY = np.meshgrid(kx, ky, indexing="ij")

    S_map = np.zeros(KX.shape)
    for f in fr:
        # exp(-i k.r) 를 (n_kx, n_ky, N) 없이 계산 — 성분별로 분해
        px, py = f[:, 0], f[:, 1]
        ex = np.exp(-1j * np.outer(kx, px))          # (n_kx, N)
        ey = np.exp(-1j * np.outer(ky, py))          # (n_ky, N)
        rho = ex @ ey.T if False else np.einsum("an,bn->ab", ex, ey)
        S_map += np.abs(rho) ** 2 / N
    S_map /= n_frames

    K = np.hypot(KX, KY)
    S_map[K == 0.0] = 0.0                            # k=0 은 N (자명) — 제외

    # 방향 평균
    k_min = 2 * np.pi / max(Lx, Ly)
    k_hi = float(K.max())
    edges = np.linspace(0.0, k_hi, bins + 1)
    idx = np.digitize(K.ravel(), edges) - 1
    Sr = np.zeros(bins)
    cnt = np.zeros(bins)
    flat = S_map.ravel()
    ok = (idx >= 0) & (idx < bins) & (K.ravel() > 0)
    np.add.at(Sr, idx[ok], flat[ok])
    np.add.at(cnt, idx[ok], 1.0)
    with np.errstate(invalid="ignore", divide="ignore"):
        Sr = np.where(cnt > 0, Sr / np.maximum(cnt, 1), 0.0)
    kr = 0.5 * (edges[:-1] + edges[1:])

    valid = (kr >= k_min) & (cnt > 0)
    i = int(np.argmax(Sr[valid]))
    k1 = float(kr[valid][i])

    # 첫 봉 고리에서 6겹 변조: S(θ) 의 cos(6θ) 성분 / 평균
    ring = np.abs(K - k1) <= ring_width * k1
    if ring.sum() >= 12:
        th = np.arctan2(KY[ring], KX[ring])
        s = S_map[ring]
        c6 = float(np.abs(np.mean(s * np.exp(6j * th))) / max(np.mean(s), 1e-30))
    else:
        c6 = float("nan")

    return StructureFactorResult(
        kx=kx, ky=ky, S_map=S_map, k_radial=kr, S_radial=Sr,
        first_peak_k=k1, first_peak_S=float(Sr[valid][i]),
        sixfold_modulation=c6, k_min=float(k_min), n_frames=n_frames)


# =============================================================================
# 최소 분리 — 가드
# =============================================================================
def min_separation(frames: np.ndarray, *, Lx: float, Ly: float) -> float:
    """전 프레임에서 가장 가까운 쌍의 거리.

    `power_law_table` 의 `r_min` 보다 작아지면 **표를 벗어난다** — HOOMD 가
    조용히 외삽하거나 죽는다. 런타임 가드가 이 값을 봐야 한다.
    """
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    L = np.array([Lx, Ly])
    worst = np.inf
    for f in fr:
        d = f[:, None, :2] - f[None, :, :2]
        d -= L * np.round(d / L)
        dist = np.linalg.norm(d, axis=2)
        np.fill_diagonal(dist, np.inf)
        worst = min(worst, float(dist.min()))
    return worst


# =============================================================================
# 유한크기 스케일링 — `ψ₆` 의 `N` 의존성이 상을 판별한다
# =============================================================================
#  ★★ `|⟨ψ₆⟩|` 의 절대값은 상을 말해주지 않는다. 유한계에서는 무질서한 계도
#    `1/√N` 규모의 값을 낸다 (`N = 100` 이면 `0.1`). **갈리는 것은 `N` 의존성이다.**
#
#    `|⟨ψ₆⟩|² = N^{-2} Σ_ij ⟨ψ₆ᵢ ψ₆ⱼ*⟩` 이므로
#      · 액체 (`g₆` 지수 감쇠, 상관길이 `ξ`):  Σ ~ N ξ²  ⇒  `|⟨ψ₆⟩| ~ N^{-1/2}`
#      · hexatic (`g₆ ~ r^{-η₆}`):            Σ ~ N L^{2-η₆}  ⇒  `|⟨ψ₆⟩| ~ N^{-η₆/4}`
#      · 결정 (`g₆ → const`):                 `|⟨ψ₆⟩| ~ N⁰`
#
#    ⇒ 지수 `p ≡ -d ln|⟨ψ₆⟩| / d ln N` 에서 **`η₆ = 4p`** 다.
#      KTHNY 의 hexatic-액체 경계는 `η₆ = 1/4` → `p = 1/16 = 0.0625`.
#      `p ≈ 0.5` 는 액체, `p ≲ 0.0625` 는 hexatic 이하다.
#
#  근거: knowledge/source/papers/1999-zahn-two-stage-melting-2d.md §6-3 이 요구하는
#    `η₆` 를 **두 계 크기만으로** 얻는 경로다 (`g₆(r)` 직접 적합 없이).
KTHNY_ETA6_HEXATIC_LIQUID = 0.25        # hexatic → 등방 액체 경계
LIQUID_EXPONENT_P = 0.5                 # `|⟨ψ₆⟩| ~ N^{-1/2}`


@dataclass
class FiniteSizeExponent:
    """`|⟨ψ₆⟩| ~ N^{-p}` 의 두 점(또는 그 이상) 추정. **판정하지 않는다.**"""

    p: float
    p_se: float
    eta6: float                          # = 4p
    eta6_se: float
    n_points: int
    reading: str                         # 'liquid-like' | 'hexatic-or-below' | 'between'
    #  ★ 형태 검증 — 점 3개 이상에서만 의미가 있다. 두 점은 직선을 **가정**한다.
    chi2_reduced: float = float("nan")   # χ²/dof (로그-로그 직선 적합)
    residuals: tuple = ()                # ln y − 적합 (로그 잔차)
    amplitude: float = float("nan")       # 적합의 절편 exp(b₀)

    @property
    def is_liquid_like(self) -> bool:
        """`p` 가 `0.5` 와 `1σ` 안에서 일치하는가."""
        return abs(self.p - LIQUID_EXPONENT_P) <= self.p_se

    @property
    def form_is_testable(self) -> bool:
        """멱함수 **형태**를 논할 수 있는가 — 점 3개 이상 + `χ²` 계산됨."""
        return self.n_points >= 3 and np.isfinite(self.chi2_reduced)

    def form_verdict(self, chi2_max: float = 3.0) -> str:
        """`χ²/dof` 로 본 형태 판독. **판정이 아니라 보고다.**"""
        if not self.form_is_testable:
            return f"검증 불가 (점 {self.n_points}개) — 멱함수를 가정했다"
        if self.chi2_reduced <= chi2_max:
            return (f"멱함수와 모순 없음 (χ²/dof = {self.chi2_reduced:.2f} "
                    f"≤ {chi2_max:g})")
        return (f"⚠ 멱함수에서 벗어난다 (χ²/dof = {self.chi2_reduced:.2f} "
                f"> {chi2_max:g}) — 단일 지수로 요약하면 안 된다")


#  ★★ `η₆ ≤ 1/4` 는 hexatic 의 **필요조건이지 충분조건이 아니다.**
#  결정도 `η₆ ≈ 0 ≤ 1/4` 를 만족한다 (`ψ₆ → const` 이므로 `p = 0`).
#  갈리는 것은 **`ψ₆` 의 크기**다:
#      hexatic:  `ψ₆ → 0` 을 **느리게** (`N^{-η₆/4}`, 유한계에서는 작은 값)
#      결정:     `ψ₆ → O(1)` 상수
#  2026-07-29 에 이 구별 없이 판정해서 `A ≥ 13.3` 의 **결정을 hexatic 이라고 불렀다.**
PSI6_CRYSTAL_FLOOR = 0.5     # 이 위면 결정 후보 — 유한계 무질서 바닥(1/√N)의 5배 이상


def phase_from_finite_size(fit: FiniteSizeExponent, psi6_at_largest_N: float,
                           *, crystal_floor: float = PSI6_CRYSTAL_FLOOR) -> dict:
    """`(지수, ψ₆ 크기)` 두 축으로 상을 읽는다. **지수만으로는 부족하다.**

    Args:
        psi6_at_largest_N: 사다리에서 가장 큰 `N` 의 `|⟨ψ₆⟩|`.
        crystal_floor: 이 값을 넘으면 결정 후보.

    Returns: `phase` · `why` · `exponent_alone_would_say`
    """
    eta_hi = fit.eta6 + 3.0 * fit.eta6_se if np.isfinite(fit.eta6_se) else fit.eta6
    eta_lo = fit.eta6 - 3.0 * fit.eta6_se if np.isfinite(fit.eta6_se) else fit.eta6
    exponent_says = ("hexatic-or-below" if eta_lo <= KTHNY_ETA6_HEXATIC_LIQUID
                     else "not-hexatic")

    if psi6_at_largest_N >= crystal_floor:
        phase = "crystal"
        why = (f"`ψ₆({fit.n_points} 점 중 최대 N) = {psi6_at_largest_N:.3f}` "
               f"≥ {crystal_floor:g} — `ψ₆` 가 `O(1)` 로 포화한다. "
               f"`η₆ ≈ 0` 은 결정의 서명이기도 하다")
    elif eta_lo > KTHNY_ETA6_HEXATIC_LIQUID:
        phase = "isotropic-liquid"
        why = (f"`η₆ − 3σ = {eta_lo:.2f} > {KTHNY_ETA6_HEXATIC_LIQUID:g}` — "
               f"hexatic 상한을 넘는다. 그리고 `ψ₆ = {psi6_at_largest_N:.3f}` 가 작다")
    elif eta_hi <= KTHNY_ETA6_HEXATIC_LIQUID:
        phase = "hexatic-candidate"
        why = (f"`η₆ + 3σ = {eta_hi:.2f} ≤ {KTHNY_ETA6_HEXATIC_LIQUID:g}` 이고 "
               f"`ψ₆ = {psi6_at_largest_N:.3f}` 가 결정 바닥 {crystal_floor:g} "
               f"아래다 — 준장거리 질서 후보")
    else:
        phase = "inconclusive"
        why = (f"`η₆ = {fit.eta6:.2f} ± {fit.eta6_se:.2f}` 가 상한 "
               f"{KTHNY_ETA6_HEXATIC_LIQUID:g} 을 걸친다")
    return {"phase": phase, "why": why,
            "exponent_alone_would_say": exponent_says,
            "psi6_at_largest_N": float(psi6_at_largest_N),
            "crystal_floor": crystal_floor,
            "form_ok": bool(fit.form_is_testable
                            and fit.chi2_reduced <= 3.0)}


def psi6_finite_size_exponent(n_particles: np.ndarray, psi6: np.ndarray,
                              psi6_se: np.ndarray | None = None
                              ) -> FiniteSizeExponent:
    """`|⟨ψ₆⟩| ~ N^{-p}` 의 지수. 두 점이면 해석해, 셋 이상이면 가중 최소제곱.

    Args:
        psi6_se: 각 점의 SE. 주면 `p` 의 오차를 전파한다 (없으면 `nan`).

    ⚠ **두 점으로는 멱함수 형태를 검증할 수 없다** — 가정하고 지수만 뽑는다.
      그 사실이 `n_points` 에 드러난다. 셋 이상이어야 형태를 논할 수 있다.
    """
    N = np.asarray(n_particles, dtype=np.float64)
    y = np.asarray(psi6, dtype=np.float64)
    if N.size != y.size or N.size < 2:
        raise ValueError(f"점이 {N.size}개 — 최소 2개, 길이가 같아야 한다")
    if np.any(y <= 0):
        raise ValueError("psi6 에 0 이하가 있다 — 로그를 취할 수 없다")

    lnN, lny = np.log(N), np.log(y)
    chi2_red, resid, amp = float("nan"), (), float("nan")
    if N.size == 2:
        p = -(lny[1] - lny[0]) / (lnN[1] - lnN[0])
        if psi6_se is not None:
            se = np.asarray(psi6_se, dtype=np.float64)
            #  d(ln y) = dy/y 이므로 상대오차가 전파된다
            rel = se / y
            p_se = float(np.hypot(rel[0], rel[1]) / abs(lnN[1] - lnN[0]))
        else:
            p_se = float("nan")
        amp = float(np.exp(lny[0] + p * lnN[0]))
    else:
        #  로그 공간의 오차는 상대오차다: σ_lny = σ_y / y
        sig = (np.ones_like(y) if psi6_se is None
               else np.maximum(np.asarray(psi6_se, dtype=np.float64) / y, 1e-12))
        w = 1.0 / sig ** 2
        A = np.vstack([np.ones_like(lnN), lnN]).T
        cov = np.linalg.inv(A.T @ np.diag(w) @ A)
        beta = cov @ (A.T @ (w * lny))
        p = -float(beta[1])
        p_se = float(np.sqrt(cov[1, 1]))
        amp = float(np.exp(beta[0]))
        #  ★ 형태 검증: 잔차가 오차막대와 정합하는가 (dof = n − 2)
        r = lny - A @ beta
        dof = int(N.size - 2)
        if dof > 0 and psi6_se is not None:
            chi2_red = float(np.sum((r / sig) ** 2) / dof)
        resid = tuple(float(x) for x in r)

    p = float(p)
    #  ★ `p_se` 가 `nan` 이면 (SE 를 안 준 경우) `max(nan, x)` 가 `nan` 이 되어
    #    모든 비교가 False 로 떨어진다 — 판독이 조용히 'between' 이 된다.
    #    SE 를 모를 때는 고정 허용폭을 쓴다.
    tol = 0.05 if not np.isfinite(p_se) else max(p_se, 0.05)
    reading = ("liquid-like" if abs(p - LIQUID_EXPONENT_P) <= tol
               else "hexatic-or-below"
               if p <= KTHNY_ETA6_HEXATIC_LIQUID / 4.0 + max(tol, 0.02)
               else "between")
    return FiniteSizeExponent(p=p, p_se=p_se, eta6=4.0 * p,
                              eta6_se=4.0 * p_se, n_points=int(N.size),
                              reading=reading, chi2_reduced=chi2_red,
                              residuals=resid, amplitude=amp)


# =============================================================================
# Zahn 상도 환산 — 문헌 대조의 유일한 창
# =============================================================================
#  Γ = π^{3/2} A   (A = βU(d), d = n^{-1/2})
#  근거: knowledge/source/papers/1999-zahn-two-stage-melting-2d.md §2
ZAHN_GAMMA_OVER_A = float(np.pi ** 1.5)          # 5.568328...
ZAHN_GAMMA_MELT = 59.88                          # 결정 → hexatic
ZAHN_GAMMA_ISO = 55.87                           # hexatic → 등방액체


def zahn_phase(amplitude: float) -> dict:
    """`A` → Zahn 상도 위치. **`reproduced: no` 인 문헌값이다** (`[출처, 미재현]`)."""
    G = ZAHN_GAMMA_OVER_A * amplitude
    if G > ZAHN_GAMMA_MELT:
        phase = "crystal"
    elif G > ZAHN_GAMMA_ISO:
        phase = "hexatic"
    else:
        phase = "isotropic-liquid"
    return {
        "amplitude": amplitude, "gamma": G, "phase_zahn": phase,
        "distance_to_melt": (G - ZAHN_GAMMA_MELT) / ZAHN_GAMMA_MELT,
        "distance_to_iso": (G - ZAHN_GAMMA_ISO) / ZAHN_GAMMA_ISO,
        "citation": "[출처, 미재현] Zahn 1999 PRL 82, 2721",
    }


def amplitude_for_gamma(gamma: float) -> float:
    """역변환 — 원하는 `Γ` 를 주는 `A`."""
    return gamma / ZAHN_GAMMA_OVER_A
