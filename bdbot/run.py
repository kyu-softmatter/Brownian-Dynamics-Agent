r"""L4 -- the layer that runs from the spec alone, plus the **numerical-health
verdict.**

* **L4 is not a physics verifier.** After a run it looks only at **divergence,
  NaN/Inf and strange convergence.** Comparison against analytic solutions and the
  literature is already handled by the `role` system in `bdbot.metrics`, and the
  rigour is invested in the front end (L0 -> L3).

## What lives here and what stays in the case

Here (L4, system-independent):
  * `StepGuard` runtime monitoring -- non-finite PE/force, PE blow-up, an
                over-large single-step displacement (`Guard` is an alias;
                `bdbot.health.Guard` is a *different* class -- see there)
  * `judge`     post-run verdict -- nan / diverged / frozen / drifting / ok
  * `execute`   run-directory lifecycle, equilibration + production loops, sample
                collection, emitting metrics

The case (its own physics):
  * `build(spec) -> Build`  builds the HOOMD system from the spec's dimensionless
                            parameters and supplies the sampling function

**The data contract is the spec, and nothing else** -- `build` takes only a
`LoadedSpec` and never re-reads the case YAML. The **code** that builds the physics
is necessarily per-system, so it is attached through a registry
(`@builder("case-name")`). What `nondim show` guarantees self-sufficiency for is
the **data** side, and that is not broken here.

## The verdict criteria -- why these four

| Status | What was observed | Why it is a numerical problem |
|---|---|---|
| `nan` | PE or a force is non-finite | the integration broke |
| `diverged` | \|PE\|/N blew up, or a single-step displacement > `step_disp_max` | dt is too large |
| `frozen` | the variance of PE/N over the production phase is 0 | there is thermal noise and nothing moves -- the force went to zero (pair.Table trap 11 has exactly this shape) |
| `drifting` | first-half/second-half mean difference > 3 * block SEM | insufficient equilibration (a warning, not a failure) |

`frozen` came out of a real trap in this project: `pair.Table` gives force and
energy **exactly 0** for `r < r_min`, so particles quietly stop while overlapped.
It does not blow up, so a NaN check does not catch it, and PE does not explode
either. **Not changing** is the only signal.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np

from . import metrics as MET
from . import runid as RID
from . import stats as ST
from .health import check_finite

SCHEMA = "bdbot.run/0.1"

OK, DIVERGED, NAN, FROZEN = "ok", "diverged", "nan", "frozen"
DRIFT_SIGMA = 3.0            # above this many block SEMs of half-to-half difference, under-equilibrated
GUARD_EVERY = 10_000         # every 10^4 steps


# ══════════════════════════════════════════════════════════════════════
# what the case returns
# ══════════════════════════════════════════════════════════════════════
@dataclass
class Phase:
    """One phase of a run. Introduced for **systems with a driving protocol.**

    `trap-drag`'s question was "equilibrium energy vs energy while dragging vs
    relaxation after stopping", which needs three phases of different character
    inside one trajectory. Two stages (eq -> prod) cannot express that.

    * `expect_steady=False` matters. In a relaxation phase **the energy changing IS
      the physics**, so running the drift check there produces a false warning.
      The verdict is applied only to phases where a steady state is expected.
    """

    name: str
    n_steps: int
    sample_every: int = 0            # 0 means use Build.sample_every
    collect: bool = True             # False discards the samples (initial relaxation)
    expect_steady: bool = True       # False skips the drift check
    note: str = ""


@dataclass
class Build:
    """What the case builds from the spec. L4 takes only this and runs."""

    sim: Any                                   # hoomd.Simulation
    forces: list                               # for diagnostics (the force arrays are read)
    n_particles: int
    sample: Callable[..., dict]                # (timestep, phase) -> {name: value}
    pe_per_particle: Callable[[], float]       # <U>/N [kT] -- the indicator for the equilibrium and frozen verdicts
    n_eq: int = 0
    n_prod: int = 0
    sample_every: int = 1
    phases: list = field(default_factory=list)   # if empty, built from n_eq/n_prod
    gsd_path: Any = None
    tags: list = field(default_factory=list)   # the system_tags in metrics
    physical: dict = field(default_factory=dict)
    finalize: Callable[[dict], dict] | None = None    # samples -> observables and extras

    def plan(self) -> list:
        """The phases to run. A case that supplies none gets the original two:
        equilibration -> production.
        """
        if self.phases:
            return list(self.phases)
        return [Phase("equilibration", self.n_eq, GUARD_EVERY, collect=False),
                Phase("production", self.n_prod, self.sample_every)]


# ══════════════════════════════════════════════════════════════════════
# the registry -- a case registers its own builder
# ══════════════════════════════════════════════════════════════════════
BUILDERS: dict = {}


def builder(case: str):
    def deco(fn):
        BUILDERS[case] = fn
        return fn
    return deco


def get_builder(case: str):
    if case not in BUILDERS:
        raise KeyError(f"no L4 builder is registered for '{case}'. "
                       f"registered: {', '.join(sorted(BUILDERS)) or '(none)'}")
    return BUILDERS[case]


# ══════════════════════════════════════════════════════════════════════
# the runtime guard
# ══════════════════════════════════════════════════════════════════════
class Diverged(RuntimeError):
    """Raised when a guard trips. `status` distinguishes the kind."""

    def __init__(self, status: str, step: int, msg: str):
        super().__init__(f"[{status}] step {step}: {msg}")
        self.status, self.step, self.msg = status, step, msg


@dataclass
class StepGuard:
    """Every `10^4` steps, check for non-finite values, blow-up and over-large
    displacement.

    `step_disp_max` bounds **the force-driven displacement of one step**,
    `dt*|F|/gamma` (dimensionless gamma=1). BD is O(dt), so being pushed more than
    10% of a diameter in one step already makes the integration untrustworthy.
    WARNING: the thermal displacement `sqrt(2 dt)` is deliberately excluded -- that
       is physics, and what is being judged is the **force term.**

    ★ Renamed from `Guard` on 2026-08-29 because there were **two classes called
      `Guard`** in this package and they watch different things:

        bdbot.health.Guard   positions and PE -- box escape, relative PE blow-up.
                             Takes `box_L`; attachable as a hoomd action
        bdbot.run.StepGuard  PE/N and the force arrays -- absolute PE bound,
                             non-finite forces, force-driven step displacement.
                             Takes `dt_star`; driven by `execute()`

      Same-name-different-job is exactly the confusion that already cost this
      project once: `run.Guard` wrote `dt*|F|max` into `l4` while the health tool
      looked for `numerics["step_rms_sigma"]`, found nothing, and reported
      `82/82 HEALTHY` (see `bdbot/health.py`'s module docstring).
      `Guard = StepGuard` is kept below as an alias so nothing breaks.
    """

    dt_star: float
    n_particles: int
    pe_abs_max: float = 1e8          # |PE|/N [kT] -- above this it is a blow-up, not physics
    step_disp_max: float = 0.1       # [σ]
    trips: list = field(default_factory=list)
    # * The **worst value over the whole run.** Keeping only the last sample misses
    #   the worst moment in between -- measured, the peak force was 1062.9 against
    #   244.2 kT/sigma for the last sample, a factor of 4.4. A `dt/tau_fast` verdict
    #   must be made on the worst value, not the mean (it is a stability criterion).
    f_max_seen: float = 0.0
    step_disp_seen: float = 0.0

    def check(self, step: int, pe_per_n: float, forces) -> None:
        # ★ non-finite is defined once, in `bdbot.health.check_finite` (section 0).
        #   The messages stay distinct because the two guards report different
        #   quantities -- only the predicate is shared.
        if not check_finite(pe_per_n=[pe_per_n])[0]:
            raise Diverged(NAN, step, f"PE/N is non-finite ({pe_per_n!r})")
        if abs(pe_per_n) > self.pe_abs_max:
            raise Diverged(DIVERGED, step,
                           f"PE/N = {pe_per_n:.4e} kT exceeds the bound {self.pe_abs_max:.0e}")
        fmax = 0.0
        for f in forces:
            arr = np.asarray(f.forces, dtype=float)
            if arr.size == 0:
                continue
            if not check_finite(force=arr)[0]:
                raise Diverged(NAN, step, f"non-finite value in {type(f).__name__}'s force array")
            fmax = max(fmax, float(np.abs(arr).max()))
        disp = self.dt_star * fmax
        self.f_max_seen = max(self.f_max_seen, fmax)
        self.step_disp_seen = max(self.step_disp_seen, disp)
        if disp > self.step_disp_max:
            raise Diverged(DIVERGED, step,
                           f"single-step force displacement dt*|F|max = {disp:.4g} sigma "
                           f"exceeds the bound {self.step_disp_max} sigma "
                           f"(|F|max = {fmax:.4g} kT/sigma). Reduce dt")
        self.trips.append({"step": step, "pe_per_n": pe_per_n, "f_max": fmax,
                           "step_disp": disp})


# ══════════════════════════════════════════════════════════════════════
# post-run verdict -- numerical only
# ══════════════════════════════════════════════════════════════════════
def judge(pe_series, *, status: str = OK, expect_steady: bool = True,
          label: str = "production") -> dict:
    """Judge the numerical health of one `<U>/N` series.

    `frozen` is the key one -- in a system with thermal noise, PE/N being
    **exactly constant** means the force went to zero, and that does not blow up so
    NaN and blow-up checks do not catch it.

    * With `expect_steady=False` the drift check is skipped -- issuing a drift
      warning for a phase where **the energy changing IS the physics**, such as a
      relaxation phase, means calling a discovery a warning (rule 7').
    """
    s = np.asarray(pe_series, dtype=float)
    out = {"status": status, "n_samples": int(s.size), "warnings": []}
    if status != OK:
        return out
    if s.size == 0:
        out["status"] = NAN
        out["warnings"].append(f"{label} has 0 samples")
        return out
    if not np.all(np.isfinite(s)):
        out["status"] = NAN
        return out

    out["pe_mean"] = float(s.mean())
    out["pe_std"] = float(s.std())
    if s.size >= 2 and out["pe_std"] == 0.0:
        out["status"] = FROZEN
        out["warnings"].append(
            f"[{label}] <U>/N is exactly constant ({s[0]:.6g}) across {s.size} samples "
            "-- there is thermal noise and it does not change. The force may have "
            "gone to zero (pair.Table trap 11: F=0 for r<r_min)")
        return out

    # equilibration verdict -- half-to-half mean difference vs the block SEM.
    # A warning, not a failure.
    if expect_steady and s.size >= 8:
        h = s.size // 2
        d = float(s[h:].mean() - s[:h].mean())
        sem = float(ST.block_sem(s))
        out["drift"] = d
        out["drift_sem"] = sem
        out["drift_sigma"] = abs(d) / sem if sem > 0 else float("inf")
        if sem > 0 and abs(d) > DRIFT_SIGMA * sem:
            out["warnings"].append(
                f"[{label}] <U>/N is drifting: second half - first half = {d:+.4g} kT "
                f"= {abs(d)/sem:.1f} * block SEM (limit {DRIFT_SIGMA:g}). "
                f"Equilibration is insufficient")
    return out


# ══════════════════════════════════════════════════════════════════════
# execution -- one spec -> one run directory
# ══════════════════════════════════════════════════════════════════════
def execute(spec, build_fn, outdir, *, force: bool = False, progress: bool = True,
            guard_every: int = GUARD_EVERY, extra_metrics=None) -> dict:
    """Run one spec and leave a `metrics.json`. Returns the verdict dict.

    `spec` is a `bdbot.nondim.LoadedSpec` -- that is, a `specs/<run_id>.json` read
    back.

    WARNING: this does not write `result.txt`. The case script does. A case that
    omits it has its runs counted as zero by `bdbot.cli status`.
    """
    from . import sim as SIM       # pull hoomd in only here

    outdir = Path(outdir)
    go, msg = RID.prepare_outdir(outdir, force)
    if not go:
        print(msg)
        return {"status": "skipped", "run_id": spec.run_id}

    ok_hash, want = spec.verify_hash()
    if not ok_hash:
        raise ValueError(f"this spec was hand-edited -- the run_id does not match its "
                         f"contents (expected {want}). Rule 2.")
    if spec.verdict.startswith("FAIL"):
        raise ValueError(f"a spec whose L3 verdict is FAIL is not run: {spec.verdict}")

    b = build_fn(spec, outdir)          # * pass outdir, so the GSD path never enters the spec
    dt_star = float(spec.numerics["dt_star"])
    guard = StepGuard(dt_star=dt_star, n_particles=b.n_particles)
    t0 = time.time()
    status, trip = OK, None

    def _loop(n_steps, label, collect, every):
        """The loop shared by equilibration and production. Stops every `every`
        steps to run the guard and take a sample.

        * The guard period and the sampling period are made **the same** (the
          smaller of the two). Run separately, a divergence occurring between
          samples never reaches the samples.
        """
        nonlocal status, trip
        done = 0
        every = max(1, min(every, n_steps)) if n_steps else 1
        next_print = 0
        while done < n_steps:
            chunk = min(every, n_steps - done)
            b.sim.run(chunk)
            done += chunk
            pe = b.pe_per_particle()
            try:
                guard.check(b.sim.timestep, pe, b.forces)
            except Diverged as e:
                status, trip = e.status, {"step": e.step, "msg": e.msg, "phase": label}
                print(f"\n  x runtime guard: {e}", flush=True)
                return False, done
            collect(b.sim.timestep, pe)
            if progress and done >= next_print:
                print(SIM.progress(done, n_steps, time.time() - t0,
                                   f"{label}  ⟨U⟩/N={pe:9.4f}"), flush=True)
                next_print += max(1, n_steps // 10)
        return True, done

    # -- the phases, in order -------------------------------------------------
    plan = b.plan()
    samples: list = []
    pe_series: list = []          # for the verdict -- the last phase expecting a steady state
    phase_of: list = []
    eq_trace: list = []
    per_phase: dict = {}
    n_done = 0
    alive = True
    for ph in plan:
        if not alive or ph.n_steps <= 0:
            continue
        pe_ph: list = []
        print(f"\n  > {ph.name}  {ph.n_steps:,} steps"
              + (f"   {ph.note}" if ph.note else ""), flush=True)

        def collect(ts, pe, _ph=ph, _pe=pe_ph):
            _pe.append(pe)
            if _ph.collect:
                s = b.sample(ts, _ph.name)
                s["_t_step"] = ts
                samples.append(s)
                phase_of.append(_ph.name)
                pe_series.append(pe)
            else:
                eq_trace.append((ts, pe))

        every = ph.sample_every or b.sample_every
        alive, done = _loop(ph.n_steps, ph.name, collect,
                            min(every, guard_every) if ph.collect else guard_every)
        n_done += done
        per_phase[ph.name] = judge(pe_ph, status=OK if alive else status,
                                   expect_steady=ph.expect_steady, label=ph.name)
        per_phase[ph.name]["steps"] = done
        per_phase[ph.name]["expect_steady"] = ph.expect_steady

    SIM.flush_writers(b.sim)
    wall = time.time() - t0

    # the overall verdict = the worst of the phase verdicts. Drift was already
    # checked per phase, so it is skipped here.
    verdict = judge(pe_series, status=status, expect_steady=False)
    for name, pv in per_phase.items():
        if pv["status"] != OK and verdict["status"] == OK:
            verdict["status"] = pv["status"]
        verdict["warnings"].extend(pv.get("warnings", []))
    verdict["phases"] = per_phase
    verdict.update(run_id=spec.run_id, wall_seconds=wall,
                   steps_done=n_done, steps_planned=sum(p.n_steps for p in plan),
                   trip=trip, guard_samples=len(guard.trips))
    if guard.trips:
        last = guard.trips[-1]
        verdict["f_max_last"] = last["f_max"]
        verdict["step_disp_last"] = last["step_disp"]
        verdict["f_max_seen"] = guard.f_max_seen
        verdict["step_disp_seen"] = guard.step_disp_seen

    # -- artefacts ------------------------------------------------------------
    cols = {}
    if samples:
        for k in samples[0]:
            try:
                cols[k] = np.asarray([s[k] for s in samples])
            except Exception:
                pass
    cols["pe"] = np.asarray(pe_series, dtype=float)
    cols["phase"] = np.asarray(phase_of)            # so it can be sliced per phase
    if eq_trace:
        cols["eq_trace"] = np.asarray(eq_trace, dtype=float)

    obs, extra = [], {}
    if b.finalize is not None and verdict["status"] == OK:
        res = b.finalize(cols) or {}
        obs, extra = res.get("observables", []), res.get("extra", {})
        # So a case can send derived arrays (g(r) and the like) to the npz. Putting
        # them in JSON bloats metrics.json with thousands of numbers and makes it
        # hard for the post-mortem to read.
        cols.update({k: np.asarray(v) for k, v in res.get("arrays", {}).items()})

    np.savez_compressed(outdir / "observables.npz",
                        **{k: v for k, v in cols.items() if isinstance(v, np.ndarray)})

    m = MET.build(
        run_id=spec.run_id, case=spec.case,
        system_tags=list(b.tags),
        reference_scales={k: spec.raw["reference"][k]["symbol"]
                          for k in ("length", "energy", "time")},
        physical=dict(b.physical),
        dimensionless={g.name: g.value for g in spec.groups},
        checks=[(c, "design") for c in spec.checks],
        observables=obs,
        equilibration=MET.equilibration_series("pe", "⟨U⟩/N [kT]"),
        numerics={**{k: v for k, v in spec.numerics.items()
                     if isinstance(v, (int, float))},
                  "steps_done": n_done,
                  # * The **input** to the L4 -> L3 feedback. `health.step_health()`
                  #   reads this key and compares it against L3's predicted
                  #   dt/tau_fast. The guard used to compute this value and put it
                  #   only inside `l4`, under a different name from the key health
                  #   looked for (`step_rms_sigma`), so **all 81 runs read "not
                  #   measured"** -- the module's core check never ran once.
                  **({"step_drift_max_sigma": guard.step_disp_seen,
                      "f_max_kT_per_sigma": guard.f_max_seen} if guard.trips else {})},
        wall_seconds=wall,
        steps_per_second=(n_done / wall) if wall > 0 else None,
        extra={"l4": verdict, **({"result": extra} if extra else {})})
    MET.write(outdir, m)
    (outdir / "spec.json").write_text(json.dumps(spec.raw, indent=2, ensure_ascii=False))
    (outdir / "l4.json").write_text(json.dumps(verdict, indent=2, ensure_ascii=False,
                                               default=str))
    return verdict


def render_verdict(v: dict) -> str:
    """One human-readable paragraph. Shared by `cli run` and the case scripts."""
    mark = {OK: "✓", DIVERGED: "✗", NAN: "✗", FROZEN: "✗"}.get(v["status"], "?")
    L = ["", "=" * 84,
         f"L4 numerical health: {mark} {v['status'].upper()}   run_id={v.get('run_id','?')}",
         "=" * 84,
         f"  steps    {v.get('steps_done', 0):,} / {v.get('steps_planned', 0):,}"
         f"   wall {v.get('wall_seconds', 0):.1f}s",
         f"  guard    passed {v.get('guard_samples', 0)}x"
         + (f"   max force {v['f_max_last']:.4g} kT/sigma"
            f"   step displacement {v['step_disp_last']:.3e} sigma" if "f_max_last" in v else "")]
    if v.get("trip"):
        L.append(f"  x ABORTED [{v['trip']['phase']}] step {v['trip']['step']}: "
                 f"{v['trip']['msg']}")
    if v.get("phases"):
        L.append(f"  {'phase':<14}{'steps':>12}{'<U>/N [kT]':>16}{'s.d.':>11}"
                 f"{'drift':>8}   steady expected")
        for name, p in v["phases"].items():
            drift = f"{p['drift_sigma']:.1f}σ" if "drift_sigma" in p else "—"
            L.append(f"  {name:<12}{p.get('steps', 0):>12,}"
                     f"{p.get('pe_mean', float('nan')):>16.6g}"
                     f"{p.get('pe_std', float('nan')):>11.4g}{drift:>8}"
                     f"   {'yes' if p.get('expect_steady') else 'no (changing IS the physics)'}")
    elif "pe_mean" in v:
        L.append(f"  ⟨U⟩/N   {v['pe_mean']:.6g} ± {v.get('pe_std', 0):.4g} kT"
                 + (f"   drift {v['drift_sigma']:.1f} sigma" if "drift_sigma" in v else ""))
    for w in v.get("warnings", []):
        L.append(f"  ⚠ {w}")
    L.append("  note: L4 looks at numerics only -- physics comparison is the role "
             "system in metrics.observables")
    L.append("=" * 84)
    return "\n".join(L)


#  Backward-compatible alias: `run.Guard` was the name for 8 cases' worth of
#  history. Prefer `StepGuard` in new code.
Guard = StepGuard

__all__ = ["SCHEMA", "OK", "DIVERGED", "NAN", "FROZEN", "Build", "StepGuard", "Guard", "Diverged",
           "builder", "get_builder", "BUILDERS", "judge", "execute", "render_verdict"]
