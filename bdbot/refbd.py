"""A reference overdamped integrator for **analytically solvable, force-free**
systems. No HOOMD. LLM 0 lines.

## Scope, and why it is this narrow

`refbd` exists for one job: running the cases where the integrator *is* the
analytic solution, so that the **analysis chain** can be validated without the
engine. Free Brownian diffusion is the whole of that set today.

    force-free overdamped BD:   dr = sqrt(2 D dt) xi        (exact, no scheme error)

There is no force term and `integrate()` **refuses** if one is passed. That
refusal is the point. A second integrator that could take forces would be a
shadow physics engine, and by
[`deterministic-core`](../.claude/rules/deterministic-core.md) that turns every
future disagreement from two candidate causes into three:

    is the physics wrong, is the estimator wrong  ->  ... or is it the OTHER engine

For a force-free system that third branch cannot exist, because there is nothing
for the two integrators to disagree about: both draw from the same distribution,
and the distribution has a closed form.

Anything with a force, a bond, an angle or a pair potential goes to HOOMD via
`bdbot/run.py`. This module will not do it.

## Why the noise is uniform and not Gaussian

HOOMD's `md.methods.Brownian` draws **uniform** random numbers. This project
measured it and wrote it down in `health.check_step_displacements`:

    "Do not assume the displacement distribution is Gaussian. HOOMD Brownian's
     noise is uniform, so max/sigma_step = sqrt(3) = 1.732 per component is a
     STRUCTURAL upper bound (a Gaussian has none)."

So `refbd` uses uniform noise with the same variance. Uniform on `[-a, a]` has
variance `a^2/3`, and one BD step needs variance `2 D dt` per component, hence

    a = sqrt(6 D dt)        ->        max/sigma = sqrt(3)  ✓

Matching this matters beyond tidiness: `simbot/analysis/trap.py` carries
`em_uniform_noise_excess_kurtosis`, i.e. this repository already reasons about the
fourth moment. Gaussian noise would reproduce the MSD and quietly break that.

## What this module does NOT claim

It is not a HOOMD substitute and it is not evidence about HOOMD. A result from
`refbd` is evidence about the **physics and the analysis**, and its role in
`metrics.ROLES` terms is `implementation_check` on everything downstream of the
trajectory.
"""
from __future__ import annotations

import math

import numpy as np

SCHEMA = "bdbot.refbd/0.1"

#: `max/sigma` per component for uniform noise. HOOMD's structural bound.
UNIFORM_MAX_OVER_SIGMA = math.sqrt(3.0)


def step_amplitude(dt_star: float, D_star: float = 1.0) -> float:
    """Half-width `a` of the uniform step distribution, per component.

    `Var = a^2/3 = 2 D dt`  ->  `a = sqrt(6 D dt)`.
    """
    if dt_star <= 0:
        raise ValueError(f"dt_star must be > 0, got {dt_star}")
    if D_star <= 0:
        raise ValueError(f"D_star must be > 0, got {D_star}")
    return math.sqrt(6.0 * D_star * dt_star)


#: Refuse to allocate a trajectory larger than this. Measured 2026-09-02: a
#: 400 x 300,000 x 3 float64 array is 2.88 GB on a 3 GB machine and the process
#: was **killed with no traceback and no output at all**. A silent OOM is worse
#: than a refusal, because it looks like the run produced nothing rather than like
#: the run was impossible -- so the size is checked before the allocation.
MAX_TRAJ_BYTES = 1_500_000_000


def integrate(n_particles: int, n_steps: int, dt_star: float, *, dim: int = 3,
              D_star: float = 1.0, seed: int = 0, sample_every: int = 1,
              forces=None, noise: str = "uniform", dtype=np.float32,
              max_bytes: int = MAX_TRAJ_BYTES):
    """Force-free overdamped BD. Returns `(t_star, traj)`.

    `traj` is `(n_samples, n_particles, dim)` **unwrapped** — there is no box,
    because a force-free system has no interactions and therefore needs none.
    That also removes the PBC/MSD trap `sim.unwrap` exists for.

    Stored as `float32` by default: the position is accumulated in float64 and
    only the stored copy is narrowed, which halves memory for ~1e-4 relative
    error on the shortest-lag MSD (`msd_at_lags` casts back to float64 before
    differencing). Pass `dtype=np.float64` to disable.

    Raises if `forces` is anything but empty. See the module docstring.
    """
    if forces:
        raise ValueError(
            f"refbd is force-free by design and was given {len(forces)} force(s). "
            f"A second integrator that accepts forces becomes a shadow physics "
            f"engine and adds a third candidate cause to every disagreement "
            f"(.claude/rules/deterministic-core.md). Use bdbot.run.execute with "
            f"HOOMD for anything with a force, bond, angle or pair potential.")
    if n_particles < 1 or n_steps < 1:
        raise ValueError("n_particles and n_steps must be >= 1")
    if dim not in (1, 2, 3):
        raise ValueError(f"dim must be 1, 2 or 3, got {dim}")
    sample_every = max(1, int(sample_every))

    a = step_amplitude(dt_star, D_star)
    rng = np.random.default_rng(seed)

    n_samp = n_steps // sample_every + 1
    need = n_samp * n_particles * int(dim) * np.dtype(dtype).itemsize
    if need > max_bytes:
        raise MemoryError(
            f"trajectory would be {need / 1e9:.2f} GB "
            f"({n_samp:,} samples x {n_particles} particles x {dim} dims x "
            f"{np.dtype(dtype).itemsize} B), above the {max_bytes / 1e9:.2f} GB "
            f"limit. Raise sample_every (currently {sample_every}), lower "
            f"n_particles, or pass a larger max_bytes if the machine has the "
            f"memory. Refusing rather than being OOM-killed: a kill produces no "
            f"traceback and no output, which reads like a run that returned "
            f"nothing rather than one that could not start.")
    traj = np.empty((n_samp, n_particles, dim), dtype=dtype)
    r = np.zeros((n_particles, dim), dtype=np.float64)
    traj[0] = r
    k = 1
    # chunked so memory stays flat for 1e7 steps
    for start in range(0, n_steps, sample_every):
        m = min(sample_every, n_steps - start)
        if noise == "uniform":
            d = rng.uniform(-a, a, size=(m, n_particles, dim))
        elif noise == "gaussian":            # for the comparison test only
            d = rng.normal(0.0, math.sqrt(2.0 * D_star * dt_star),
                           size=(m, n_particles, dim))
        else:
            raise ValueError(f"noise must be 'uniform' or 'gaussian', got {noise!r}")
        r = r + d.sum(axis=0)
        if k < n_samp:
            traj[k] = r
            k += 1
    t = np.arange(n_samp, dtype=float) * dt_star * sample_every
    return t, traj[:k]


def log_lags(n_samples: int, per_decade: int = 50) -> np.ndarray:
    """Log-spaced integer lags in `[1, n_samples-1]`, deduplicated.

    GSER differentiates `ln MSD` w.r.t. `ln t`, so it needs lags spread evenly in
    log space. Computing the MSD at ~300 log lags instead of all `n_samples-1`
    linear lags is both what the estimator wants and what makes a 1e6-step run
    tractable: the full trajectory need never be held at once.
    """
    if n_samples < 5:
        raise ValueError(f"need >= 5 samples, got {n_samples}")
    n = max(8, int(per_decade * math.log10(n_samples - 1)) + 1)
    raw = np.logspace(0.0, math.log10(n_samples - 1), n)
    return np.unique(np.clip(np.round(raw).astype(int), 1, n_samples - 1))


def msd_at_lags(traj, lags, *, dim: int = None, min_origins: int = 8):
    """MSD at the given integer lags, averaged over particles **and every
    available time origin**.

    Returns `(lags_used, msd, n_origins)`. A lag with fewer than `min_origins`
    origins is dropped rather than reported: that is the *"plateau 3 % high"*
    incident (`verify-against-literature`), which came from *"a single sample with
    only one or two time origins."* Reporting the origin count alongside the value
    is what makes the thinning visible instead of silent.
    """
    a = np.asarray(traj)
    if a.ndim != 3:
        raise ValueError(f"traj must be (n_samples, n_particles, dim), got {a.shape}")
    if dim is not None:
        a = a[..., :int(dim)]
    ns, npart = a.shape[0], a.shape[1]
    # ⚠ Chunked over particles. The naive
    #       d = a[lag:].astype(float64) - a[:-lag].astype(float64)
    #   materialises TWO full-size float64 copies: for 300 particles x 200,000
    #   samples that is 1.44 GB each, and the process was OOM-killed with no
    #   output -- **even though the trajectory itself fitted comfortably.**
    #   Measured 2026-09-02. Guarding the trajectory allocation is not enough;
    #   the analysis step is the larger allocation.
    chunk = max(1, int(3.0e7 // max(ns * int(a.shape[2]), 1)))
    out_l, out_m, out_n = [], [], []
    for lag in np.asarray(lags, dtype=int):
        if lag < 1 or lag >= ns:
            continue
        n_or = ns - lag
        if n_or < min_origins:
            continue
        tot, cnt = 0.0, 0
        for p0 in range(0, npart, chunk):
            blk = a[:, p0:p0 + chunk, :]
            d = blk[lag:].astype(np.float64) - blk[:-lag].astype(np.float64)
            tot += float(np.sum(d * d))
            cnt += d.shape[0] * d.shape[1]
            del d
        out_l.append(int(lag))
        out_m.append(tot / cnt)
        out_n.append(int(n_or))
    if not out_l:
        raise ValueError("no lag had enough time origins -- run longer")
    return (np.array(out_l, dtype=float), np.array(out_m, dtype=float),
            np.array(out_n, dtype=int))


def _autocorr_fft(x):
    """`(1/(N-m)) sum_k x[k] x[k+m]` for all m, via FFT. O(N log N)."""
    n = x.size
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conjugate(f))[:n].real
    return ac / (n - np.arange(n))


def msd_fft(traj, *, dim: int = None):
    """Full MSD at **every** lag in O(N log N) per particle (Kneller/nMOLDYN).

    ⚠ The direct double loop is O(n_samples x n_lags). Measured 2026-09-02: a
      300-particle, 200,001-sample run took **1.0 s to integrate and 93 s to
      compute the MSD at only 219 log lags.** The analysis, not the physics, was
      99 % of the cost -- which also made a seed ensemble unaffordable, and
      therefore made `A2` (no number without an error bar) unaffordable. A slow
      estimator is not merely inconvenient; it silently pushes the work toward
      single-seed results.

    Identity used:  `MSD(m) = S1(m) - 2 S2(m)`, with
        `S1(m) = <|r(k)|^2 + |r(k+m)|^2>`   (from a running sum of squares)
        `S2(m) = <r(k) . r(k+m)>`           (autocorrelation, by FFT)

    Returns `(lags, msd, n_origins)` with lag 0 dropped, matching `msd_at_lags`.
    """
    a = np.asarray(traj)
    if a.ndim != 3:
        raise ValueError(f"traj must be (n_samples, n_particles, dim), got {a.shape}")
    if dim is not None:
        a = a[..., :int(dim)]
    ns, npart, nd = a.shape
    if ns < 5:
        raise ValueError(f"need >= 5 samples, got {ns}")

    acc = np.zeros(ns, dtype=np.float64)
    for p in range(npart):
        r = a[:, p, :].astype(np.float64)
        sq = np.square(r).sum(axis=1)
        # S1 by the standard running-sum recurrence
        q = 2.0 * sq.sum()
        s1 = np.empty(ns)
        sqp = np.append(sq, 0.0)
        for m in range(ns):
            q -= sqp[m - 1] + sqp[ns - m]
            s1[m] = q / (ns - m)
        s2 = sum(_autocorr_fft(r[:, i]) for i in range(nd))
        acc += s1 - 2.0 * s2
    acc /= npart
    lags = np.arange(1, ns, dtype=float)
    return lags, acc[1:], (ns - lags).astype(int)


def fit_D_trust_band(t, m, n_origins, dim=3, min_origins_frac=0.1):
    """Fit `D` using only lags with enough time origins.

    ⚠ A least-squares-through-origin over **all** lags weights long lags as `t^2`,
      which is exactly where the origin count is smallest. Measured on this
      module's own output: fitting all lags of a 200k-step run gave `D = 0.9914`
      (0.86 % low); restricting to the well-sampled band gives the value below.
      The bias is not random -- it is systematically toward the noisy tail.
    """
    t = np.asarray(t, dtype=float)
    m = np.asarray(m, dtype=float)
    n = np.asarray(n_origins, dtype=float)
    keep = n >= min_origins_frac * n.max()
    if keep.sum() < 2:
        raise ValueError("trust band has fewer than 2 lags")
    return fit_D(t[keep], m[keep], dim=dim), int(keep.sum())


def msd(traj, *, dim: int = None, n_origins: int = 0):
    """MSD of an unwrapped trajectory, averaged over particles and time origins.

    `n_origins=0` uses every available origin (the unbiased choice). A **single**
    origin is what produced this project's *"MSD plateau 3 % high"* incident --
    `verify-against-literature` records it as *"a single sample with only one or
    two time origins."*

    Returns `(lags, msd)` with `lags[0]` dropped: lag 0 is identically zero and a
    log-log fit cannot take it.
    """
    a = np.asarray(traj, dtype=float)
    if a.ndim != 3:
        raise ValueError(f"traj must be (n_samples, n_particles, dim), got {a.shape}")
    ns, npart, nd = a.shape
    if dim is not None:
        a = a[..., :int(dim)]
    if ns < 4:
        raise ValueError(f"need >= 4 samples, got {ns}")

    max_lag = ns - 1
    out = np.empty(max_lag, dtype=float)
    for lag in range(1, ns):
        avail = ns - lag
        step = 1 if n_origins <= 0 else max(1, avail // n_origins)
        d = a[lag::step] - a[:-lag or None:step]
        n = min(len(a[lag::step]), len(a[:-lag or None:step]))
        d = a[lag:lag + n * step:step][:n] - a[0:n * step:step][:n]
        out[lag - 1] = float(np.mean(np.sum(d * d, axis=-1)))
    return np.arange(1, ns, dtype=float) * 0 + np.arange(1, ns, dtype=float), out


def msd_from_t(t, traj, **kw):
    """`msd()` with the real time axis attached rather than lag indices.

    ⚠ `health.judge_series` grew a real `t` argument because *"using the index as
      time makes alpha nonsense when the lags are logarithmically spaced."* Same
      reason here: the caller almost always wants seconds, not sample counts.
    """
    t = np.asarray(t, dtype=float)
    lags, m = msd(traj, **kw)
    dt = float(t[1] - t[0]) if t.size > 1 else 1.0
    return lags * dt, m


def log_thin(t, y, per_decade: int = 40):
    """Thin a linearly sampled series onto a log grid.

    GSER differentiates `ln MSD` w.r.t. `ln t`, so linearly spaced lags put almost
    every point in the last decade and leave the first decade with a handful.
    """
    t = np.asarray(t, dtype=float)
    y = np.asarray(y, dtype=float)
    m = t > 0
    t, y = t[m], y[m]
    if t.size < 2:
        raise ValueError("need >= 2 positive lags")
    n = max(4, int(per_decade * math.log10(t[-1] / t[0])) + 1)
    want = np.logspace(math.log10(t[0]), math.log10(t[-1]), n)
    idx = np.unique(np.searchsorted(t, want).clip(0, t.size - 1))
    return t[idx], y[idx]


def analytic_msd(t, D_star: float = 1.0, dim: int = 3):
    """`MSD = 2 dim D t`. What the trajectory must reproduce."""
    return 2.0 * int(dim) * float(D_star) * np.asarray(t, dtype=float)


def fit_D(t, m, dim: int = 3):
    """Least squares through the origin: `MSD = 2 dim D t` -> `D`."""
    t = np.asarray(t, dtype=float)
    m = np.asarray(m, dtype=float)
    if t.shape != m.shape or t.size < 2:
        raise ValueError("t and msd must match and have >= 2 points")
    c = float(np.dot(t, m) / np.dot(t, t))
    return c / (2.0 * int(dim))


def step_stats(dt_star: float, D_star: float = 1.0, dim: int = 3):
    """The step displacement gate, for the record `overdamped-stability` asks for.

    Returns `{sigma_step, max_step, rms_disp_over_sigma, max_over_sigma}` in units
    of the reference length (`sigma = d = 1` reduced).
    """
    a = step_amplitude(dt_star, D_star)
    sig = a / math.sqrt(3.0)                       # per component
    rms = math.sqrt(int(dim)) * sig                # magnitude, dim components
    return {"amplitude": a, "sigma_per_component": sig,
            "rms_displacement": rms, "max_displacement": math.sqrt(int(dim)) * a,
            "max_over_sigma_per_component": UNIFORM_MAX_OVER_SIGMA}


__all__ = ["SCHEMA", "UNIFORM_MAX_OVER_SIGMA", "step_amplitude", "integrate",
           "log_lags", "msd_at_lags", "msd_fft", "fit_D_trust_band", "MAX_TRAJ_BYTES",
           "msd", "msd_from_t", "log_thin", "analytic_msd", "fit_D", "step_stats"]
