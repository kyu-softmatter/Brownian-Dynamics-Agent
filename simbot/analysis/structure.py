"""S7 — 2D structure analysis. It measures only; it does not judge.

Uses `freud` 3.5. RDF · hexagonal order parameter `ψ₆` · Voronoi defects · `S(q)`.

## What this module is careful about

**① The length unit is `d = n^{-1/2}`** (not the nearest-neighbour distance
`a = 1.0746 d`). `Γ ∝ d^{-3}`, so a 7 % confusion makes `Γ` wrong by 23 %.

**② `ψ₆` depends on the neighbour definition.** Voronoi neighbours and "the nearest
6" give different answers -- especially near a defect. **Which one was used is
reported alongside.**

**③ Frames are not treated as independent.** Structural relaxation is slow. The
caller fixes the seed ensemble and the frame spacing, and this module reports
`n_frames` as it is.
Basis: [[ks-test-needs-independent-samples]]
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np


def _box(Lx: float, Ly: float):
    import freud
    return freud.box.Box(Lx=Lx, Ly=Ly, is2D=True)


def _xyz(pos_2d: np.ndarray) -> np.ndarray:
    """freud requires (N, 3) even in 2D (z = 0)."""
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
        """Mean neighbours inside `r < r_max`: `n(r) = 2πn ∫ g(r) r dr` (2D)."""
        m = self.r <= r_max
        return float(2 * np.pi * density_star
                     * np.trapezoid(self.g[m] * self.r[m], self.r[m]))


def rdf(frames: np.ndarray, *, Lx: float, Ly: float, r_max: float | None = None,
        bins: int = 200) -> RDFResult:
    """`g(r)` averaged over several frames. `frames` is `(n_frames, N, 2)`."""
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
    # the first peak — searched for at r > 0.5 d (inside that there is almost no
    # sample)
    m = r > 0.5
    i = int(np.argmax(g[m]))
    return RDFResult(r=r, g=g, first_peak_r=float(r[m][i]),
                     first_peak_g=float(g[m][i]), n_frames=int(fr.shape[0]))


# =============================================================================
# hexagonal order parameter + Voronoi defects
# =============================================================================
@dataclass
class HexOrderResult:
    psi6_global: float              # |<psi6>| — global orientational order
    psi6_local_mean: float          # <|psi6_i|> — local order (always >= global)
    psi6_per_frame: list[float]
    defect_fraction: float          # fraction with Voronoi coordination != 6
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
    """`(ψ₆ᵢ, coordinationᵢ)` for one frame. Shared by `hex_order` and
    `hex_order_series`.

    ★ Why it is shared: if the time-resolved version computed it separately the two
      functions would drift apart quietly. "The mean of the second half" and "the
      late mean of the per-frame series" **have to be the same number.**
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
        raise ValueError(f"neighbor_mode {neighbor_mode!r} must be "
                         f"'voronoi' or 'nearest6'")
    counts = np.bincount(np.asarray(voro.nlist.query_point_indices),
                         minlength=pts.shape[0])
    return np.asarray(hexatic.particle_order), counts


def hex_order(frames: np.ndarray, *, Lx: float, Ly: float,
              neighbor_mode: str = "voronoi") -> HexOrderResult:
    """`ψ₆` and the Voronoi defect fraction.

    Args:
        neighbor_mode: `"voronoi"` (recommended) or `"nearest6"`.
            ★ The two give different answers -- near a defect Voronoi counts the 5
            and 7 neighbours as they are, while `nearest6` forces 6 of them.
            **To see defects, use Voronoi.**

    **Both** `psi6_global = |⟨ψ₆ᵢ⟩|` and `psi6_local_mean = ⟨|ψ₆ᵢ|⟩` are returned:
      · local large only → there are grains but their orientations differ
        (polycrystal)
      · both large → a single crystal
      · local large with global 0 → a hexatic candidate
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
# Time-resolved — **when** the structure gets made
# =============================================================================
#  ★ Why this is kept separate from the time average
#  A second-half mean answers "what did it become" and **does not answer "when did
#  it become that."** The same final `ψ₆` is indistinguishable between (ⅰ) reached
#  at 1 τ_d and stayed there and (ⅱ) still climbing -- and that distinction is the
#  whole of an equilibration diagnosis.
#
#  ⚠ **Frames are not independent samples** (module §③). Using the frame-to-frame
#    scatter of a time series as the statistical error underestimates it, because of
#    the time correlation. Error bars come from the seed ensemble.
@dataclass
class HexOrderSeries:
    """Per-frame `ψ₆`, defects and coordination. The time axis is the caller's
    `t_star`."""

    t_star: np.ndarray                  # (n_frames,)
    psi6_global: np.ndarray             # (n_frames,) |⟨ψ₆ᵢ⟩| — magnitude of the
                                        # within-frame mean
    psi6_local: np.ndarray              # (n_frames,) ⟨|ψ₆ᵢ|⟩
    defect_fraction: np.ndarray         # (n_frames,) fraction with coordination ≠ 6
    coord_labels: np.ndarray            # (n_coord,) coordination values ([3,4,…,10])
    coord_fraction: np.ndarray          # (n_frames, n_coord) per-frame fractions
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
        """`|n₅ − n₇|/(n₅ + n₇)` per frame. **0 is the signature of dislocations
        (5-7 pairs).**

        The defect fraction alone does not distinguish a dislocation from a liquid
        (card §8.2) -- this quantity keeps that distinction along the time axis.
        """
        n5, n7 = self._col(5), self._col(7)
        return np.abs(n5 - n7) / np.maximum(n5 + n7, 1e-12)

    @property
    def coord_kinds(self) -> np.ndarray:
        """Per frame, how many coordination kinds have fraction > 0.5 %. A liquid
        has 6 kinds, a crystal 3."""
        return (self.coord_fraction > 0.005).sum(axis=1)

    def window_mean(self, t_lo: float, t_hi: float) -> dict:
        """Mean over `t_lo ≤ t < t_hi`. **Used to compare against a second-half
        mean.**"""
        m = (self.t_star >= t_lo) & (self.t_star < t_hi)
        if not m.any():
            raise ValueError(f"no frames inside [{t_lo}, {t_hi})")
        return {"n_frames": int(m.sum()),
                "psi6_global": float(self.psi6_global[m].mean()),
                "psi6_local": float(self.psi6_local[m].mean()),
                "defect_fraction": float(self.defect_fraction[m].mean())}


def hex_order_series(frames: np.ndarray, *, Lx: float, Ly: float,
                     t_star: np.ndarray | None = None,
                     neighbor_mode: str = "voronoi",
                     coord_range: tuple[int, int] = (3, 10)) -> HexOrderSeries:
    """The **time-resolved** version of `hex_order` — it returns the frames as they
    are instead of averaging them.

    `hex_order(frames).psi6_global` is **exactly equal** to this function's
    `psi6_global.mean()` (they use the same `_frame_hex`). A test watches that.

    Args:
        coord_range: the coordination histogram range (both ends inclusive). Values
            outside are **discarded** rather than accumulated into the last bin --
            `defect_fraction`, by contrast, counts all of them.
    """
    import freud
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    n_frames, N = fr.shape[0], fr.shape[1]
    t = (np.arange(n_frames, dtype=np.float64) if t_star is None
         else np.asarray(t_star, dtype=np.float64))
    if t.size != n_frames:
        raise ValueError(f"t_star length {t.size} != frame count {n_frames}")

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
    """A `y(t) = y_inf + (y_0 − y_inf) exp(−t/τ)` fit. **It does not judge.**"""

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
        """The relaxation amplitude `y_0 − y_inf`. Smaller than the noise and `τ`
        is meaningless."""
        return self.y0 - self.y_inf


def fit_relaxation(t_star: np.ndarray, y: np.ndarray, *,
                   noise: float | None = None) -> RelaxationFit:
    """Single-exponential relaxation time.

    ★ **Before reporting `τ`, the `amplitude` has to be compared against the
      noise.** Fitting an exponential to a signal fluctuating about a steady state
      always yields some `τ` -- and that is noise, not relaxation. Given `noise`,
      that judgment is reflected in `converged`.

    Args:
        noise: the signal's frame-to-frame standard deviation. If
            `|amplitude| < 2·noise`, `converged=False` and the reason is left in
            `note`.
    """
    from scipy.optimize import curve_fit

    t = np.asarray(t_star, dtype=np.float64)
    yy = np.asarray(y, dtype=np.float64)
    if t.size != yy.size:
        raise ValueError(f"t({t.size}) and y({yy.size}) have different lengths")
    if t.size < 4:
        raise ValueError(f"{t.size} points cannot support a 3-parameter fit")

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
                             converged=False, note=f"fit failed: {e!r}")
    resid = yy - model(t, *popt)
    ss_tot = float(((yy - yy.mean()) ** 2).sum())
    r2 = 1.0 - float((resid ** 2).sum()) / ss_tot if ss_tot > 0 else float("nan")
    tau_se = float(np.sqrt(pcov[2, 2])) if np.all(np.isfinite(pcov)) else float("nan")

    note = ""
    amp = float(popt[0] - popt[1])
    if noise is not None and abs(amp) < 2.0 * noise:
        ok = False
        note = (f"relaxation amplitude |{amp:.4g}| < 2×noise({noise:.4g}) — this is "
                f"steady-state fluctuation, not relaxation. Do not report τ")
    if np.isfinite(tau_se) and tau_se > abs(popt[2]):
        ok = False
        note = note or (f"τ = {popt[2]:.4g} ± {tau_se:.4g} — the error exceeds the "
                        f"value")
    return RelaxationFit(y0=float(popt[0]), y_inf=float(popt[1]),
                         tau=float(popt[2]), tau_se=tau_se, r_squared=r2,
                         n_points=int(t.size), converged=ok, note=note)


def bootstrap_relaxation_over_seeds(
        t_star: np.ndarray, per_seed: np.ndarray, *, n_resample: int = 400,
        seed: int = 0, noise: float | None = None) -> dict:
    """The **seed-ensemble** error on `τ`. It measures something different from the
    `curve_fit` covariance.

    ★ The `tau_se` from `fit_relaxation` is the fit uncertainty of *one seed-mean
      curve*. This repository's convention is that **"the seed ensemble is the
      honest error estimate"** (`scripts/soft2d_sweep_analyze.py` §①), so the
      correct error is the spread of the distribution obtained by resampling the
      seeds and refitting.

      A large difference between the two means **the fit uncertainty does not
      represent the seed-to-seed variation**, and then any σ distance built on the
      `curve_fit` SE cannot be trusted.

    Args:
        per_seed: `(n_seeds, n_frames)` — the per-seed curves. Pass them as they
            are, do not average.
        n_resample: number of resamples with replacement.

    Returns: `tau` (fit of the overall mean curve) · `tau_se_bootstrap` ·
        `tau_ci95` · `n_converged` · `tau_se_fit` (for comparison).
    """
    y = np.asarray(per_seed, dtype=np.float64)
    if y.ndim != 2:
        raise ValueError(f"per_seed must be (n_seeds, n_frames) — got {y.shape}")
    k = y.shape[0]
    if k < 2:
        raise ValueError(f"{k} seeds cannot give an ensemble error")

    full = fit_relaxation(t_star, y.mean(axis=0), noise=noise)
    rng = np.random.default_rng(seed)
    taus: list[float] = []
    for _ in range(int(n_resample)):
        idx = rng.integers(0, k, size=k)               # resample with replacement
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
                "note": f"the bootstrap fit converged only {arr.size} times — no "
                        f"error can be claimed"}
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
    """`g(r)` per time window. **Averaged within a window, never mixed across
    windows.**"""

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
    """Split the trajectory into `n_windows` time windows and take `g(r)` in each.

    ★ A single window's `g(r)` is noisy because it has few frames -- more windows
      raises the time resolution and lowers the signal-to-noise. **The caller
      decides that trade.**
    """
    fr = np.asarray(frames, dtype=np.float64)
    if fr.ndim == 2:
        fr = fr[None]
    t = np.asarray(t_star, dtype=np.float64)
    if t.size != fr.shape[0]:
        raise ValueError(f"t_star length {t.size} != frame count {fr.shape[0]}")
    if n_windows < 1:
        raise ValueError("n_windows ≥ 1")

    edges = np.linspace(t[0], t[-1] + (t[-1] - t[0]) * 1e-9, n_windows + 1)
    r_ref, gs, los, his, cnts, pr, pg = None, [], [], [], [], [], []
    for w in range(n_windows):
        m = (t >= edges[w]) & (t < edges[w + 1])
        if not m.any():
            raise ValueError(
                f"window {w} [{edges[w]:.4g}, {edges[w + 1]:.4g}) has no frames — "
                f"n_windows({n_windows}) is large relative to the frame count "
                f"({t.size})")
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
    """Voronoi tiling of one frame — **for figures**. The measuring is done by
    `hex_order_series`."""

    positions: np.ndarray               # (N, 2)
    polygons: list                      # length N, each (n_v, 2)
    coordination: np.ndarray            # (N,) integers
    psi6_abs: np.ndarray                # (N,) |ψ₆ᵢ| — local hexagonality
    Lx: float
    Ly: float


def voronoi_frame(frame: np.ndarray, *, Lx: float, Ly: float,
                  neighbor_mode: str = "voronoi") -> VoronoiFrame:
    """Voronoi polygons + coordination + local `ψ₆`. The material for card §10's
    `voronoi plot`.

    ⚠ Polygon vertices **can fall outside the box** (a cell open toward a periodic
      neighbour). Clipping for the figure is the axis limits' job; clipping here
      would make the cell areas wrong.
    """
    import freud
    f = np.asarray(frame, dtype=np.float64)
    if f.ndim != 2:
        raise ValueError(f"one frame (N, 2) is required — got shape {f.shape}")
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
    """2D `S(k)` — returns **the vector map and the angular average together.**

    ★ The angular average alone cannot distinguish a hexagonal crystal from a
      hexatic. Both put their first peak at the same `|k|`. What separates them is
      **the angular dependence**:
        crystal → 6-fold **points** (Bragg)
        hexatic → a **ring** modulated 6-fold
        liquid  → an isotropic ring
      Which is why `S_map` and `sixfold_modulation` are needed.
    """

    kx: np.ndarray                  # (n_kx,) k components commensurate with the box
    ky: np.ndarray
    S_map: np.ndarray               # (n_kx, n_ky) the vector S(k)
    k_radial: np.ndarray            # |k| bin centres for the angular average
    S_radial: np.ndarray
    first_peak_k: float
    first_peak_S: float
    sixfold_modulation: float       # 6-fold strength on the first-peak ring
                                    # (0 = isotropic)
    k_min: float                    # the smallest |k| the box allows = 2π/max(L)
    n_frames: int


def structure_factor(frames: np.ndarray, *, Lx: float, Ly: float,
                     n_max: int = 24, bins: int = 120,
                     ring_width: float = 0.15) -> StructureFactorResult:
    """Computes the 2D static structure factor **directly**.

    `S(k) = |Σ_j exp(−i k·r_j)|² / N`, `k = 2π(m/Lx, n/Ly)` — only k-vectors
    commensurate with the box are used, so **there is no sampling error** (an exact
    discrete sum, no approximation).

    ⚠ **`freud` is not used.** `freud.diffraction.StaticStructureFactorDirect`
      raises `ValueError: 2D boxes are not currently supported` (freud 3.5.0,
      measured 2026-07-28). 2D is a first-class target for this project, so this is
      implemented directly.

    Args:
        n_max: the k-grid radius (integer index). `|k|_max ≈ 2π n_max/min(L)`
        ring_width: relative thickness of the ring the 6-fold modulation is measured
            on (relative to the first peak's `|k|`)
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
        # compute exp(-i k.r) without an (n_kx, n_ky, N) array — factored per axis
        px, py = f[:, 0], f[:, 1]
        ex = np.exp(-1j * np.outer(kx, px))          # (n_kx, N)
        ey = np.exp(-1j * np.outer(ky, py))          # (n_ky, N)
        rho = ex @ ey.T if False else np.einsum("an,bn->ab", ex, ey)
        S_map += np.abs(rho) ** 2 / N
    S_map /= n_frames

    K = np.hypot(KX, KY)
    S_map[K == 0.0] = 0.0                            # k=0 is N (trivial) — excluded

    # angular average
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

    # 6-fold modulation on the first-peak ring: the cos(6θ) component of S(θ) / mean
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
# Minimum separation — a guard
# =============================================================================
def min_separation(frames: np.ndarray, *, Lx: float, Ly: float) -> float:
    """The closest pair distance over all frames.

    Fall below `power_law_table`'s `r_min` and it **leaves the table** -- HOOMD
    either extrapolates quietly or dies. A runtime guard has to watch this value.
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
# Finite-size scaling — `ψ₆`'s `N` dependence discriminates the phase
# =============================================================================
#  ★★ The absolute value of `|⟨ψ₆⟩|` says nothing about the phase. In a finite
#    system even a disordered one produces a value of order `1/√N` (`0.1` at
#    `N = 100`). **What separates them is the `N` dependence.**
#
#    Since `|⟨ψ₆⟩|² = N^{-2} Σ_ij ⟨ψ₆ᵢ ψ₆ⱼ*⟩`:
#      · liquid (`g₆` decays exponentially, correlation length `ξ`):
#        Σ ~ N ξ²  ⇒  `|⟨ψ₆⟩| ~ N^{-1/2}`
#      · hexatic (`g₆ ~ r^{-η₆}`):  Σ ~ N L^{2-η₆}  ⇒  `|⟨ψ₆⟩| ~ N^{-η₆/4}`
#      · crystal (`g₆ → const`):    `|⟨ψ₆⟩| ~ N⁰`
#
#    ⇒ from the exponent `p ≡ -d ln|⟨ψ₆⟩| / d ln N`, **`η₆ = 4p`**.
#      KTHNY's hexatic-liquid boundary is `η₆ = 1/4` → `p = 1/16 = 0.0625`.
#      `p ≈ 0.5` is a liquid; `p ≲ 0.0625` is hexatic or below.
#
#  Basis: this is the path to the `η₆` that
#    knowledge/source/papers/1999-zahn-two-stage-melting-2d.md §6-3 requires, using
#    **two system sizes alone** (no direct `g₆(r)` fit).
KTHNY_ETA6_HEXATIC_LIQUID = 0.25        # the hexatic → isotropic liquid boundary
LIQUID_EXPONENT_P = 0.5                 # `|⟨ψ₆⟩| ~ N^{-1/2}`


@dataclass
class FiniteSizeExponent:
    """A two-point (or more) estimate of `|⟨ψ₆⟩| ~ N^{-p}`. **It does not judge.**"""

    p: float
    p_se: float
    eta6: float                          # = 4p
    eta6_se: float
    n_points: int
    reading: str                         # 'liquid-like' | 'hexatic-or-below' | 'between'
    #  ★ Form verification — only meaningful with 3 or more points. Two points
    #    **assume** a straight line.
    chi2_reduced: float = float("nan")   # χ²/dof (log-log straight-line fit)
    residuals: tuple = ()                # ln y − fit (log residuals)
    amplitude: float = float("nan")       # the fit's intercept exp(b₀)

    @property
    def is_liquid_like(self) -> bool:
        """Does `p` agree with `0.5` to within `1σ`?"""
        return abs(self.p - LIQUID_EXPONENT_P) <= self.p_se

    @property
    def form_is_testable(self) -> bool:
        """Can the power-law **form** be discussed — 3+ points and `χ²` computed."""
        return self.n_points >= 3 and np.isfinite(self.chi2_reduced)

    def form_verdict(self, chi2_max: float = 3.0) -> str:
        """Reading the form via `χ²/dof`. **A report, not a verdict.**"""
        if not self.form_is_testable:
            return f"not testable ({self.n_points} points) — power law assumed"
        if self.chi2_reduced <= chi2_max:
            return (f"consistent with a power law (χ²/dof = "
                    f"{self.chi2_reduced:.2f} ≤ {chi2_max:g})")
        return (f"⚠ departs from a power law (χ²/dof = {self.chi2_reduced:.2f} "
                f"> {chi2_max:g}) — must not be summarized by a single exponent")


#  ★★ `η₆ ≤ 1/4` is a **necessary but not sufficient** condition for a hexatic.
#  A crystal also satisfies `η₆ ≈ 0 ≤ 1/4` (`ψ₆ → const`, so `p = 0`).
#  What separates them is **the magnitude of `ψ₆`**:
#      hexatic:  `ψ₆ → 0` **slowly** (`N^{-η₆/4}`, a small value in a finite system)
#      crystal:  `ψ₆ → O(1)`, constant
#  On 2026-07-29 a verdict was issued without this distinction and **called the
#  crystal at `A ≥ 13.3` a hexatic.**
PSI6_CRYSTAL_FLOOR = 0.5     # above this is a crystal candidate — at least 5x the
                             # finite-size disorder floor (1/√N)


def phase_from_finite_size(fit: FiniteSizeExponent, psi6_at_largest_N: float,
                           *, crystal_floor: float = PSI6_CRYSTAL_FLOOR) -> dict:
    """Reads the phase on two axes, `(exponent, ψ₆ magnitude)`. **The exponent
    alone is not enough.**

    Args:
        psi6_at_largest_N: `|⟨ψ₆⟩|` at the largest `N` on the ladder.
        crystal_floor: above this value it is a crystal candidate.

    Returns: `phase` · `why` · `exponent_alone_would_say`
    """
    eta_hi = fit.eta6 + 3.0 * fit.eta6_se if np.isfinite(fit.eta6_se) else fit.eta6
    eta_lo = fit.eta6 - 3.0 * fit.eta6_se if np.isfinite(fit.eta6_se) else fit.eta6
    exponent_says = ("hexatic-or-below" if eta_lo <= KTHNY_ETA6_HEXATIC_LIQUID
                     else "not-hexatic")

    if psi6_at_largest_N >= crystal_floor:
        phase = "crystal"
        why = (f"`ψ₆(largest N of {fit.n_points} points) = "
               f"{psi6_at_largest_N:.3f}` ≥ {crystal_floor:g} — `ψ₆` saturates at "
               f"`O(1)`. `η₆ ≈ 0` is also the signature of a crystal")
    elif eta_lo > KTHNY_ETA6_HEXATIC_LIQUID:
        phase = "isotropic-liquid"
        why = (f"`η₆ − 3σ = {eta_lo:.2f} > {KTHNY_ETA6_HEXATIC_LIQUID:g}` — "
               f"exceeds the hexatic upper bound. And "
               f"`ψ₆ = {psi6_at_largest_N:.3f}` is small")
    elif eta_hi <= KTHNY_ETA6_HEXATIC_LIQUID:
        phase = "hexatic-candidate"
        why = (f"`η₆ + 3σ = {eta_hi:.2f} ≤ {KTHNY_ETA6_HEXATIC_LIQUID:g}` and "
               f"`ψ₆ = {psi6_at_largest_N:.3f}` is below the crystal floor "
               f"{crystal_floor:g} — a quasi-long-range-order candidate")
    else:
        phase = "inconclusive"
        why = (f"`η₆ = {fit.eta6:.2f} ± {fit.eta6_se:.2f}` straddles the upper "
               f"bound {KTHNY_ETA6_HEXATIC_LIQUID:g}")
    return {"phase": phase, "why": why,
            "exponent_alone_would_say": exponent_says,
            "psi6_at_largest_N": float(psi6_at_largest_N),
            "crystal_floor": crystal_floor,
            "form_ok": bool(fit.form_is_testable
                            and fit.chi2_reduced <= 3.0)}


def psi6_finite_size_exponent(n_particles: np.ndarray, psi6: np.ndarray,
                              psi6_se: np.ndarray | None = None
                              ) -> FiniteSizeExponent:
    """The exponent of `|⟨ψ₆⟩| ~ N^{-p}`. Analytic for two points, weighted least
    squares for three or more.

    Args:
        psi6_se: the SE of each point. Given them, the error on `p` is propagated
            (`nan` otherwise).

    ⚠ **Two points cannot verify the power-law form** -- it is assumed and only the
      exponent is extracted. That fact shows up in `n_points`. Three or more are
      needed before the form can be discussed.
    """
    N = np.asarray(n_particles, dtype=np.float64)
    y = np.asarray(psi6, dtype=np.float64)
    if N.size != y.size or N.size < 2:
        raise ValueError(f"{N.size} points — at least 2, and the lengths must match")
    if np.any(y <= 0):
        raise ValueError("psi6 contains a value <= 0 — the log cannot be taken")

    lnN, lny = np.log(N), np.log(y)
    chi2_red, resid, amp = float("nan"), (), float("nan")
    if N.size == 2:
        p = -(lny[1] - lny[0]) / (lnN[1] - lnN[0])
        if psi6_se is not None:
            se = np.asarray(psi6_se, dtype=np.float64)
            #  d(ln y) = dy/y, so the relative error propagates
            rel = se / y
            p_se = float(np.hypot(rel[0], rel[1]) / abs(lnN[1] - lnN[0]))
        else:
            p_se = float("nan")
        amp = float(np.exp(lny[0] + p * lnN[0]))
    else:
        #  the error in log space is the relative error: σ_lny = σ_y / y
        sig = (np.ones_like(y) if psi6_se is None
               else np.maximum(np.asarray(psi6_se, dtype=np.float64) / y, 1e-12))
        w = 1.0 / sig ** 2
        A = np.vstack([np.ones_like(lnN), lnN]).T
        cov = np.linalg.inv(A.T @ np.diag(w) @ A)
        beta = cov @ (A.T @ (w * lny))
        p = -float(beta[1])
        p_se = float(np.sqrt(cov[1, 1]))
        amp = float(np.exp(beta[0]))
        #  ★ Form verification: are the residuals consistent with the error bars
        #    (dof = n − 2)
        r = lny - A @ beta
        dof = int(N.size - 2)
        if dof > 0 and psi6_se is not None:
            chi2_red = float(np.sum((r / sig) ** 2) / dof)
        resid = tuple(float(x) for x in r)

    p = float(p)
    #  ★ When `p_se` is `nan` (no SE was given), `max(nan, x)` is `nan` and every
    #    comparison falls to False -- the reading quietly becomes 'between'.
    #    With an unknown SE, a fixed tolerance is used instead.
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
# Conversion to Zahn's phase diagram — the only window onto the literature
# =============================================================================
#  Γ = π^{3/2} A   (A = βU(d), d = n^{-1/2})
#  Basis: knowledge/source/papers/1999-zahn-two-stage-melting-2d.md §2
ZAHN_GAMMA_OVER_A = float(np.pi ** 1.5)          # 5.568328...
ZAHN_GAMMA_MELT = 59.88                          # crystal → hexatic
ZAHN_GAMMA_ISO = 55.87                           # hexatic → isotropic liquid


def zahn_phase(amplitude: float) -> dict:
    """`A` → a position on Zahn's phase diagram. **These are literature values with
    `reproduced: no`** (`[source, not reproduced]`)."""
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
        "citation": "[source, not reproduced] Zahn 1999 PRL 82, 2721",
    }


def amplitude_for_gamma(gamma: float) -> float:
    """The inverse — the `A` that gives a desired `Γ`."""
    return gamma / ZAHN_GAMMA_OVER_A
