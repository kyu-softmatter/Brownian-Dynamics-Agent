"""L4 -- the numerical-health judge.

**Not a physics verifier.** By instruction:

    de-emphasize post-run scientific and physical verification -- it differs per
    system, and the system may not be the one that was reported. After a run, look
    at **numerical errors only**: divergence, NaN/Inf, convergence to a strange
    value.

Rule 7 reaches the same conclusion -- do not apply a standard theory to a
combined result. So there is no analytic solution and no literature value here.
**Only the numerical properties of the time series.**

Three parts:

  1. `Guard`          in-run monitoring. On NaN/Inf/runaway, **abort immediately.**
  2. `judge(...)`     post-run time-series verdict. Divergence, stalling, collapse.
  3. `step_health()`  * **feedback into L3.** L4 measures the `dt/tau_fast` that L3
                      predicted. A discrepancy means the scale ledger has **a
                      missing timescale.**

Part 3 is the core of this module. In the dimensionless convention
(sigma=kT=gamma=1), the deterministic displacement of one step is

    drift_per_step / σ  =  F* dt* / (γ* σ)  =  dt / τ_fast

so **the measured step displacement IS `dt/tau_fast`.** Comparing it against L3's
prediction (computed from the timescales in the ledger) checks the ledger's
completeness after the fact.

WARNING: part 3 **never ran across all 81 runs** because of a name mismatch --
`run.Guard` computed `dt*|F|max` into `l4` while the health tool looked for
`numerics["step_rms_sigma"]`, found nothing, printed "not measured" and returned
HEALTHY anyway. `82/82 HEALTHY` read like coverage. Silence is not success: the
count of unmeasured runs is printed separately now, and a HEALTHY verdict means
"no divergence, no stall, no collapse" and explicitly **not** "dt is small
enough".
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# Failure classes (the numerical ones only)
NUMERIC_MODES = ("NUM_NONFINITE", "NUM_DIVERGE", "NUM_FROZEN", "NUM_COLLAPSE",
                 "NUM_STEP_TOO_COARSE", "LEDGER_INCOMPLETE")

STEP_HARD = 1e-2        # hard limit on dt/tau_fast (same value as bd-physics section 4)
LEDGER_TOL = 3.0        # above this measured/predicted ratio, suspect the ledger


# ════════════════════════════════════════════════════════════════════════
# 0. guard primitives -- **the single source of truth**
# ════════════════════════════════════════════════════════════════════════
#  ★ Merged from `simbot/guards.py` 2026-08-29. Both packages had runtime guards
#    and neither was a superset: `bdbot` had the minimum image in its displacement
#    measure and the force-vs-thermal split, `simbot` had the configurational
#    thermometer, the bond-length check and the does-it-fluctuate assertion.
#    Union, hosted here because this is the L4 layer the engine calls; the
#    `simbot.guards` names still resolve (re-export).
#
#  The organising principle, which is `simbot/guards.py`'s and worth keeping
#  verbatim: **a guard has to watch a quantity that CAN drift systematically.**
#  HOOMD `Brownian`'s kinetic temperature is redrawn from the target distribution
#  every step, so systematic drift in it is impossible -- it cannot be a guard.
#  The configurational temperature takes that seat.
#  Basis: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
def configurational_temperature(forces, laplacian_U_total: float) -> float:
    """`kT_conf = <|grad U|^2> / <laplacian U>` (Rugh / Butler configurational temp).

    Uses positions and forces only -- no velocities -- which is what makes it valid
    in BD.

    Arguments
      forces              (N,3) or (N,d). `F = -grad U`, so `|F|^2 = |grad U|^2`
      laplacian_U_total   ensemble mean of `sum_i laplacian_i U` over all
                          particles. Analytic per potential, so the caller
                          supplies it. Harmonic trap `U = k r^2/2` on `d` active
                          axes: `d*k` per particle

    Returns `kT_conf`, which must match the input `kT`.

    ⚠ **In a pure harmonic trap this is algebraically identical to the
      `<x^2> = kT/k` check** -- it is not new information there. It becomes an
      independent check once pair interactions exist, because then `grad U` also
      picks up the neighbours.
    """
    if laplacian_U_total <= 0:
        raise ValueError(f"laplacian_U_total must be > 0, got {laplacian_U_total}")
    return float(np.mean(np.sum(np.asarray(forces, dtype=np.float64) ** 2, axis=1))
                 / laplacian_U_total)


@dataclass
class DisplacementReport:
    max_over_sigma: float
    rms_over_sigma: float
    max_over_rms: float
    n_exceeding: int
    passed: bool
    note: str = ""


def check_step_displacements(dr, sigma: float, max_frac: float = 0.10) -> DisplacementReport:
    """Does any particle move more than `max_frac` of `sigma` in one step?

    ⚠ **Do not assume the displacement distribution is Gaussian.** HOOMD
      `Brownian`'s noise is uniform, so `max/sigma_step = sqrt(3) = 1.732` per
      component is a *structural* upper bound (a Gaussian has none). An "anomalous
      beyond n sigma" rule therefore means something different in this engine --
      **judge on absolute displacement.**
      Basis: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
    """
    dr = np.asarray(dr, dtype=np.float64)
    mag = np.linalg.norm(dr, axis=1)
    rms = float(np.sqrt(np.mean(mag**2)))
    mx = float(mag.max()) if mag.size else 0.0
    n_bad = int(np.count_nonzero(mag > max_frac * sigma))
    note = ""
    if n_bad:
        note = (f"{n_bad} particle(s) moved more than {max_frac:.3g} sigma in one step "
                f"(max {mx/sigma:.4g} sigma) -- suspect dt too large, or initial overlap")
    return DisplacementReport(
        max_over_sigma=mx / sigma, rms_over_sigma=rms / sigma,
        max_over_rms=(mx / rms if rms > 0 else float("nan")),
        n_exceeding=n_bad, passed=(n_bad == 0), note=note,
    )


def check_finite(**arrays) -> tuple[bool, list[str]]:
    """NaN/Inf check. Returns the offending array names and counts.

    Both `Guard` and `run.StepGuard` route their non-finite branch through this, so
    the definition of "non-finite" exists once.
    """
    fails = []
    for name, a in arrays.items():
        a = np.asarray(a)
        n_bad = int(np.count_nonzero(~np.isfinite(a)))
        if n_bad:
            fails.append(f"{name}: {n_bad} non-finite value(s)")
    return (not fails), fails


def check_inside_box(positions, box_lengths, dims: int = 3,
                     tol: float = 1e-9) -> tuple[bool, int]:
    """Are the particles inside the box (wall leakage)? Always passes under PBC,
    because the coordinates are wrapped."""
    pos = np.asarray(positions, dtype=np.float64)[:, :dims]
    half = np.asarray(box_lengths, dtype=np.float64)[:dims] / 2.0
    outside = np.any(np.abs(pos) > half + tol, axis=1)
    n = int(np.count_nonzero(outside))
    return (n == 0), n


def check_bond_lengths(positions, bonds, target: float, tol: float = 0.05,
                       dims: int = 3) -> tuple[bool, dict]:
    """Are the bond lengths still near their target?

    **`check_finite` alone cannot catch a bonded system blowing up.** Pairing
    `constrain.Distance` with `Brownian` takes the bond length from `1.0` to
    `5.8e7` -- all finite, so `check_finite` passes. And measuring `kappa` on the
    exploded chain can still return a plausible-looking `s^-3`.
    Basis: findings/dead-end-distance-constraint-with-brownian.md

    Arguments
      bonds    (n_bonds, 2) particle indices
      target   target bond length (dimensionless)
      tol      allowed relative error. Do not set this loose in a bending
               measurement -- stretch must not contaminate bend.

    Returns `(ok, detail)`; the detail carries max/mean relative deviation and the
    worst bond's index.
    """
    pos = np.asarray(positions, dtype=np.float64)[:, :dims]
    b = np.asarray(bonds, dtype=int)
    d = np.linalg.norm(pos[b[:, 1]] - pos[b[:, 0]], axis=1)
    rel = np.abs(d - target) / target
    imax = int(np.argmax(rel))
    info = {
        "max_rel_dev": float(rel[imax]),
        "mean_rel_dev": float(rel.mean()),
        "worst_bond": imax,
        "worst_length": float(d[imax]),
        "target": float(target),
        "tol": float(tol),
        "n_violating": int(np.count_nonzero(rel > tol)),
    }
    return bool(rel.max() <= tol), info


def assert_statistic_fluctuates(samples, name: str = "statistic",
                                min_rel_std: float = 1e-12) -> None:
    """Check that a statistic taken over independent samples actually fluctuates.

    A "measurement" that does not fluctuate is not a measurement -- it is an
    **arithmetic identity.** This happened for real on 2026-07-28: subtracting the
    mean from the displacements and then measuring the cross-correlation gives
    `cross/auto = -1/(n-1)` identically, and the standard deviation over 200
    repetitions was `6.7e-20`. The result looked plausible and nearly passed.
    Basis: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md section 3

    Raises `AssertionError` when the relative standard deviation is below
    `min_rel_std`.
    """
    s = np.asarray(samples, dtype=np.float64)
    if s.size < 2:
        raise ValueError("need >= 2 samples")
    scale = max(abs(float(np.mean(s))), 1e-300)
    rel = float(np.std(s, ddof=1)) / scale
    if rel < min_rel_std:
        raise AssertionError(
            f"{name}: relative std {rel:.3e} < {min_rel_std:.3e} -- the statistic "
            f"does not fluctuate. It is probably an arithmetic identity, not a "
            f"measurement.")


# ════════════════════════════════════════════════════════════════════════
# 1. in-run monitoring
# ════════════════════════════════════════════════════════════════════════
class Guard:
    """A runtime monitor, wrapped in a `hoomd.custom.Action`.

    Does not import HOOMD -- kept as pure functions so it stays testable; hoomd is
    pulled in only inside `as_action()`.

        g = Guard(box_L=32.0)
        sim.operations.writers.append(g.as_action(period=10_000, state=sim.state,
                                                  thermo=thermo))
    """

    def __init__(self, box_L: float, pe_blowup: float = 1e3):
        self.box_L = float(box_L)
        self.pe_blowup = float(pe_blowup)
        self.pe0: float | None = None
        self.n_checks = 0
        self.history: list[tuple[int, float]] = []

    def check(self, timestep: int, positions: np.ndarray, pe) -> None:
        """RuntimeError on violation. Never passes silently."""
        self.n_checks += 1
        ok, fails = check_finite(position=positions)      # one definition, section 0
        if not ok:
            n = int((~np.isfinite(np.asarray(positions))).sum())
            raise RuntimeError(f"[NUM_NONFINITE] step {timestep}: {n} non-finite position(s)")
        far = np.abs(positions).max()
        if far > 50 * self.box_L:
            raise RuntimeError(
                f"[NUM_DIVERGE] step {timestep}: |r|max={far:.3g} > 50·L={50*self.box_L:.3g}")
        if pe is not None:
            pe = float(pe)
            if not math.isfinite(pe):
                raise RuntimeError(f"[NUM_NONFINITE] step {timestep}: PE={pe}")
            self.history.append((int(timestep), pe))
            if self.pe0 is None and pe != 0:
                self.pe0 = abs(pe)
            elif self.pe0 and abs(pe) > self.pe_blowup * self.pe0:
                raise RuntimeError(
                    f"[NUM_DIVERGE] step {timestep}: PE {pe:.3g} is "
                    f"{abs(pe)/self.pe0:.0f}x the initial {self.pe0:.3g} "
                    f"(limit {self.pe_blowup:.0f}x)")

    def as_action(self, period: int, state, thermo=None):
        import hoomd

        outer = self

        class _A(hoomd.custom.Action):
            def act(self, timestep):
                snap = state.get_snapshot()
                pe = None if thermo is None else thermo.potential_energy
                outer.check(timestep, np.asarray(snap.particles.position), pe)

        return hoomd.write.CustomWriter(action=_A(), trigger=hoomd.trigger.Periodic(period))


# ════════════════════════════════════════════════════════════════════════
# 2. post-run time-series verdict
# ════════════════════════════════════════════════════════════════════════
@dataclass
class Finding:
    ok: bool
    mode: str | None
    name: str
    detail: str


@dataclass
class HealthReport:
    findings: list = field(default_factory=list)
    failure_modes: list = field(default_factory=list)
    measured: dict = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        return "HEALTHY" if not self.failure_modes else "UNHEALTHY"

    def add(self, ok, mode, name, detail):
        self.findings.append(Finding(ok, None if ok else mode, name, detail))
        if not ok and mode and mode not in self.failure_modes:
            self.failure_modes.append(mode)

    def render(self) -> str:
        L = ["=" * 78,
             "L4 numerical health  (not physics verification -- divergence, NaN, "
             "strange convergence only)", "=" * 78]
        for f in self.findings:
            L.append(f"  {'✓' if f.ok else '✗'} {f.name:<26} {f.detail}")
        L += ["", f"VERDICT: {self.verdict}"
              + (f"   modes={self.failure_modes}" if self.failure_modes else ""), "=" * 78]
        return "\n".join(L)


def judge_series(name: str, y, rep: HealthReport, *, positive: bool = False,
                 cumulative: bool = False, t=None) -> None:
    """The numerical health of one time series. **It does not judge physical
    correctness.**

    `cumulative=True` marks **a quantity for which growing is normal** (MSD, MSAD
    and other accumulated quantities).
    * Without that distinction it false-positives. This actually bit: the MSD of
      two `abp-rod` runs was judged NUM_DIVERGE for being "1010x larger in the
      second half" -- which was simply diffusion. Adversarial testing had used only
      synthetic stationary series and missed it; **the real data caught it.**

    For a cumulative quantity the check becomes a **super-ballistic** one instead:
    the log-log slope alpha (y ~ t^alpha). No overdamped or inertial dynamics can
    exceed alpha <= 2 (ballistic), so alpha > 2.5 is numerical divergence
    regardless of the physical assumptions.

    * Always pass `t` (the real time axis). Using the index as time makes
      **alpha nonsense when the lags are logarithmically spaced.** The real data
      caught this too -- `abp-rod`'s MSD has log-spaced lags, so index-based alpha
      exceeded 2.5 and was misread as divergence.
      With `t=None`, uniform spacing is assumed.
    """
    y = np.asarray(y, dtype=float)
    if y.size < 8:
        rep.add(True, None, f"{name}", f"{y.size} sample(s) -- verdict skipped")
        return

    if not np.all(np.isfinite(y)):
        rep.add(False, "NUM_NONFINITE", f"{name} finiteness",
                f"{int((~np.isfinite(y)).sum())}/{y.size} non-finite")
        return
    rep.add(True, None, f"{name} finiteness", f"all {y.size} finite")

    m = float(np.mean(y))
    scale = abs(m) if m else float(np.max(np.abs(y)) or 1.0)

    # stalling -- the integrator never ran, or it froze completely
    rel_std = float(np.std(y)) / scale
    if rel_std < 1e-12:
        rep.add(False, "NUM_FROZEN", f"{name} variation", f"relative s.d. {rel_std:.1e} -- constant")
    else:
        rep.add(True, None, f"{name} variation", f"relative s.d. {rel_std:.3g}")

    q = max(2, y.size // 4)
    if cumulative:
        # cumulative -- growing is normal. Check the ballistic bound alpha <= 2 instead.
        tt = (np.arange(1, y.size + 1, dtype=float) if t is None
              else np.asarray(t, dtype=float))
        axis = "index (uniform spacing assumed)" if t is None else "t"
        m_ = (np.abs(y) > 0) & (tt > 0) & np.isfinite(tt)
        if m_.sum() >= 8:
            alpha = float(np.polyfit(np.log(tt[m_]), np.log(np.abs(y[m_])), 1)[0])
            rep.add(alpha <= 2.5, "NUM_DIVERGE", f"{name} growth exponent",
                    f"alpha = {alpha:.3f}  (y ~ t^alpha, axis={axis}; ballistic bound 2, limit 2.5)")
        else:
            rep.add(True, None, f"{name} growth exponent", "too few positive samples -- skipped")
    else:
        # stationary quantity -- divergence if the last quarter explodes vs the first
        a, b = np.abs(y[:q]).mean(), np.abs(y[-q:]).mean()
        growth = b / a if a else float("inf")
        rep.add(growth <= 1e3, "NUM_DIVERGE", f"{name} blow-up",
                f"last/first = {growth:.4g}x" + ("" if growth <= 1e3 else " (limit 1e3)"))

    # collapse -- a quantity that must be positive converging to 0 (not applicable
    # to a cumulative quantity, which starts at 0)
    if positive and not cumulative:
        tail = float(np.abs(y[-q:]).mean())
        if tail < 1e-12 * scale:
            rep.add(False, "NUM_COLLAPSE", f"{name} collapse", f"last-quarter mean {tail:.1e} ~ 0")
        else:
            rep.add(True, None, f"{name} collapse", f"last-quarter mean {tail:.4g}")


# ════════════════════════════════════════════════════════════════════════
# 3. * feedback into L3 -- measure the predicted dt/tau_fast
# ════════════════════════════════════════════════════════════════════════
def step_health(step_disp_rms: float | None, dt_star: float, dim: int,
                predicted_dt_over_tau: float | None, rep: HealthReport,
                *, drift_direct: float | None = None) -> None:
    """**Measure** `dt/tau_fast` from the one-step displacement and compare against
    the L3 prediction.

    There are **two** measurement routes. Use (a) when it exists -- it is strictly
    better than (b).

    (a) **force-based** (`drift_direct`, measured at runtime by `bdbot.run.Guard`):

        drift = dt* · |F*|max / γ*  =  dt/τ_fast        (γ*=σ=1)

      No thermal noise is mixed in, so **there is nothing to subtract.** And it is
      the **worst value over the whole run**, which is what a stability verdict
      wants -- measured, the peak force was 1062.9 against 244.2 kT/sigma for the
      last sample, a factor of 4.4.

    (b) **position difference** (`step_disp_rms`, from two snapshots):

        Δr = (F*/γ*)·dt*  +  √(2 D* dt*)·ξ ,    D* = 1
        drift = √(max(0, ⟨Δr²⟩ − 2·dim·dt*))

      The thermal component has to be subtracted in quadrature, and when
      drift << thermal noise that is **a difference of two similar numbers** and
      cannot be trusted (which is why it falls through to "thermally dominated --
      comparison meaningless" below). Measured: when the drift is 0.5% of the
      thermal noise, finite sampling **clips the drift to zero.**
      This is the backward-compatible route for older runs that did not use
      `run.execute`.

    `predicted_dt_over_tau` is what L3 **computed from the timescales in the
    ledger.** A measurement much larger than the prediction means **there is a
    faster timescale that is not in the ledger** -- the scale table is incomplete,
    and this is the only signal L4 can hand back to the front end (L3).
    """
    if drift_direct is not None:
        drift = float(drift_direct)
        rep.measured["dt_over_tau_fast_measured"] = drift
        rep.measured["step_method"] = "force"
        detail = (f"dt/tau_fast = {drift:.3e}  (limit {STEP_HARD:.0e}, "
                  f"force-based dt*|F|max -- worst value over the whole run)")
    else:
        thermal2 = 2.0 * dim * dt_star
        meas2 = float(step_disp_rms) ** 2
        drift = math.sqrt(max(0.0, meas2 - thermal2))
        rep.measured["step_rms_sigma"] = float(step_disp_rms)
        rep.measured["thermal_rms_sigma"] = math.sqrt(thermal2)
        rep.measured["dt_over_tau_fast_measured"] = drift
        rep.measured["step_method"] = "position"
        detail = (f"dt/tau_fast = {drift:.3e}  (limit {STEP_HARD:.0e}, "
                  f"thermal part {math.sqrt(thermal2):.3e} subtracted)")

    rep.add(drift <= STEP_HARD, "NUM_STEP_TOO_COARSE", "step displacement (measured)", detail)

    if predicted_dt_over_tau is None:
        rep.add(True, None, "L3 ledger comparison", "no L3 prediction -- comparison skipped")
        return
    rep.measured["dt_over_tau_fast_predicted"] = float(predicted_dt_over_tau)
    if drift <= 0 or predicted_dt_over_tau <= 0:
        rep.add(True, None, "L3 ledger comparison",
                "drift is 0 -- thermally dominated, comparison meaningless")
        return
    ratio = drift / predicted_dt_over_tau
    rep.measured["ledger_ratio"] = ratio
    rep.add(ratio <= LEDGER_TOL, "LEDGER_INCOMPLETE", "L3 ledger completeness",
            f"measured/predicted = {ratio:.2f}x "
            + (f"-- suspect a faster timescale missing from the ledger "
               f"(limit {LEDGER_TOL:.0f}x)"
               if ratio > LEDGER_TOL else "-- the ledger holds the fastest scale"))

    # -- the other direction: an over-conservative design is a **cost** (not a
    #    failure -- information).
    # The ledger check only looks at ratio >> 1 (a missing scale). ratio << 1 is
    # healthy, but it means "dt was set smaller than necessary" = that much extra
    # wall clock. Cost has been a recurring problem in this project (a sweep went
    # from 25 days to 1.16 days).
    # WARNING: do not state this as a conclusion. (1) the guard only samples every
    #    GUARD_EVERY steps, so it can miss the maximum between samples, and (2) L3's
    #    r_min is the **design worst-case approach distance**, which this run may
    #    never have reached. So it is reported only as "worth checking".
    if ratio < 1.0 / LEDGER_TOL and ratio > 0:
        rep.measured["dt_headroom"] = 1.0 / ratio
        rep.add(True, None, "dt headroom (cost)",
                f"measured is {ratio:.3f}x the prediction -- judged by the largest "
                f"force actually encountered there is room to raise dt by about "
                f"{1.0/ratio:.0f}x (wall clock / {1.0/ratio:.0f}). "
                f"Do not state this as a conclusion: the maximum between guard "
                f"samples may have been missed, or this run may never have reached "
                f"L3's design worst-case approach distance -- verify with a "
                f"convergence check")


def measure_step_displacement(positions_t0, positions_t1, L: float, dim: int) -> float:
    """rms of the one-step displacement from two consecutive snapshots (in sigma).

    Minimum image applied (traps 1 and 7) via `bdbot.sim.minimum_image` -- one
    definition, shared with `bdbot.traps` (merged 2026-08-29).
    """
    from .sim import minimum_image           # deferred: sim is numpy-only but this
                                            # keeps health importable on its own
    d = minimum_image(np.asarray(positions_t1, float) - np.asarray(positions_t0, float),
                      L, dims=dim)
    return float(np.sqrt((d[:, :dim] ** 2).sum(axis=1).mean()))


# ════════════════════════════════════════════════════════════════════════
# extracting the L3 prediction from a spec
# ════════════════════════════════════════════════════════════════════════
def predicted_dt_over_tau(spec) -> float | None:
    """Pull L3's predicted `dt/tau_fast` out of a `LoadedSpec`'s
    integration-resolution checks.

    Check names differ per case, so take the **largest value** among
    `kind == 'integration'` (the fastest timescale produces the largest ratio).
    """
    vals = [c.value for c in getattr(spec, "checks", [])
            if getattr(c, "kind", "") == "integration" and isinstance(c.value, (int, float))]
    return max(vals) if vals else None


def gate(spec) -> list[str]:
    """The **pre-run** gate. Returns the reasons to **block** (empty list = pass).

    L4 reads only the spec — it never imports case code (the L2<->L4 contract).

    WARNING: **soft warnings do not block.** This used to test `verdict != "PASS"`,
       which rejected `"PASS (3 warnings)"` — measured over 83 specs, **80 were
       false rejections** and **zero** of them were real hard failures. Skill
       bd-physics section 4 defines statistics and finite-size issues as warnings,
       not failures. `run.execute` was reading `startswith("FAIL")` correctly, and
       the reason the two disagreed unnoticed is that **`execute` never called
       `gate()`** — an unwired checker cannot be wrong out loud.
       Only three things block: (1) hash mismatch, (2) a hard check failing
       (FAIL), (3) an L3 integrity **error**.
    """
    problems = []
    ok, want = spec.verify_hash()
    if not ok:
        problems.append(f"run_id mismatch: stored {spec.run_id} vs computed {want} "
                        f"— was the spec hand-edited? (rule 2)")
    if spec.verdict.startswith("FAIL"):
        bad = [c.name for c in spec.checks if getattr(c, "hard", True) and not _ok(c)]
        problems.append(f"L3 verdict={spec.verdict}"
                        + (f" — hard checks failed: {bad}" if bad else ""))
    errs = [i for i in spec.raw.get("l3_issues", []) if i.get("level") == "error"]
    if errs:
        problems.append(f"L3 integrity errors ({len(errs)}): {errs[:2]}")
    return problems


def gate_notes(spec) -> list[str]:
    """Things that do **not** block the gate but a human must see.

    A gate that passes silently is not a gate.
    """
    notes = []
    soft = [c for c in spec.checks if not getattr(c, "hard", True) and not _ok(c)]
    for c in soft:
        notes.append(f"soft warning [{c.kind}] {c.name.strip()} = {c.value:.3g} "
                     f"(limit {c.limit:g}) — does not block, but it is a "
                     f"statistics / finite-size limitation")
    tight = [c for c in spec.checks
             if _ok(c) and getattr(c, "hard", True) and _margin(c) < 5.0]
    for c in tight:
        notes.append(f"thin margin [{c.kind}] {c.name.strip()} — only "
                     f"{_margin(c):.1f}x to the limit")
    warn = [i for i in spec.raw.get("l3_issues", []) if i.get("level") != "error"]
    for i in warn:
        notes.append(f"L3 {i.get('level')} [{i.get('where')}] {i.get('msg')}")
    return notes


def _margin(c) -> float:
    if getattr(c, "op", "<=") == "<=":
        return c.limit / c.value if c.value else float("inf")
    return c.value / c.limit if c.limit else float("inf")


def _ok(c) -> bool:
    return c.value <= c.limit if getattr(c, "op", "<=") == "<=" else c.value >= c.limit


__all__ = ["Guard", "HealthReport", "judge_series", "step_health",
           "measure_step_displacement", "predicted_dt_over_tau", "gate", "gate_notes",
           "NUMERIC_MODES", "STEP_HARD", "LEDGER_TOL",
           # section 0 -- shared with simbot.guards
           "configurational_temperature", "DisplacementReport",
           "check_step_displacements", "check_finite", "check_inside_box",
           "check_bond_lengths", "assert_statistic_fluctuates"]
