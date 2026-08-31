#!/usr/bin/env python
"""Simulation Bot CLI — S3 SPEC → S8 REPORT. 0 lines of LLM.

**One command surface (merged 2026-08-29).** This module holds the S5-S8 half;
`bdbot.cli` holds L0-L4 and now exposes these same commands under `pipeline`, so
there is one entry point rather than two engines with separate front ends:

    python -m bdbot.cli pipeline run <spec.yaml>     the merged entry point
    python cli.py run <spec.yaml>                    the root shim, unchanged
    python -m simbot.cli run <spec.yaml>             this module directly

All three dispatch to the same `main()`. `bdbot.cli` imports this **lazily**, inside
the command handler, so its front end does not pull in matplotlib (`simbot.viz`) —
that property is load-bearing and `tests/test_bdbot_lazy_api.py` guards it.

    python cli.py run <spec.yaml>            one spec → REPORT.md
    python cli.py resume runs/<id>           pick up a dead run (no recomputation
                                             of completed stages)
    python cli.py converge <spec.yaml>       is the answer the same when dt, N and
                                             the initial condition are shaken
    python cli.py params [--path runs]       several runs' parameters side by side
    python cli.py calibrate                  measure this machine's throughput

## What this CLI does not do

- **It does not confirm a verdict.** Only a human writes `confirmed_by`.
- **It does not run over budget.** When the estimate exceeds it, it reports and
  stops before running.
- **It does not force an unknown card through.** With no runner, it says so.

Exactly one card currently has a runner: `passive-sphere--harmonic-trap`. Another
card needs a runner built in `simbot.run` -- it is never quietly run with the trap
runner.
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from datetime import date
from pathlib import Path

import numpy as np

#  ★ Moved from the repo root to `simbot/cli.py` on 2026-08-29, merging the two
#    CLIs (docs/00-merge-decisions.md section 5, "two engines"). The repo root is
#    now two levels up. `cli.py` at the root is a shim that re-exports this module,
#    so `python cli.py run ...` and `import cli` both still work.
REPO = Path(__file__).resolve().parent.parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from simbot import io
from simbot.analysis.trap import aggregate_seeds, fit_msd, load_run
from simbot.estimators import estimate_wall_time_s
from simbot.nondim import reduce_spec
from simbot.policy import load_policy
from simbot.report import ReportInputs, write_report
from simbot.run import TrapRunConfig, run_trap_batch
from simbot.spec import SystemSpec, derive, load_prediction, validate as validate_spec
from simbot.validate import Measurement, validate_run
from simbot.viz import trap_diagnostics

RUNNERS = {"passive-sphere--harmonic-trap": "trap"}


# =============================================================================
# Shared
# =============================================================================
def _fail(msg: str, code: int = 2) -> int:
    print(f"\n⛔ {msg}", file=sys.stderr)
    return code


def _runner_for(spec: SystemSpec) -> str:
    kind = RUNNERS.get(spec.card)
    if kind is None:
        raise SystemExit(
            f"⛔ card {spec.card!r} has no runner.\n"
            f"   cards with a runner: {sorted(RUNNERS)}\n"
            f"   running another card with the trap runner quietly computes the\n"
            f"   wrong system — build a runner in `simbot/run.py` and register it\n"
            f"   in RUNNERS.")
    return kind


def _trap_configs(spec: SystemSpec, reduced, seeds: list[int]) -> list[TrapRunConfig]:
    t = spec.timing
    return [TrapRunConfig(
        dim=reduced.dim, n_particles=reduced.n_particles, dt_star=reduced.dt_star,
        equil_tau=(t.equil_in_tau.si if t.equil_in_tau else 10.0),
        prod_tau=(t.prod_in_tau.si if t.prod_in_tau else 40.0),
        sample_interval_tau=(t.sample_interval_in_tau.si
                             if t.sample_interval_in_tau else 2.0),
        box_over_l_trap=reduced.box_star[0], k_star=reduced.k_star or 1.0,
        seed=s, label=f"d{reduced.dim}_dt{reduced.dt_star:g}_s{s}")
        for s in seeds]


def _print_cost(reduced, n_seeds: int, pol) -> tuple[float, int, bool]:
    k = min(n_seeds, pol.concurrency())
    steps = reduced.equil_steps + reduced.prod_steps
    wall = estimate_wall_time_s(reduced.n_particles, steps, k)
    over = wall > pol.wall_budget_s
    print(f"  N={reduced.n_particles} · steps/seed={steps} · {n_seeds} seeds "
          f"· concurrency {k} (efficiency {pol.efficiency(k):.3f})")
    print(f"  estimated wall ≈ {wall:.2f} s  (budget {pol.wall_budget_s:g} s)"
          + ("  ⚠️ **over**" if over else ""))
    return wall, k, over


# =============================================================================
# Measurement — batch results → Measurement
# =============================================================================
def measure_trap_batch(batch: dict, reduced, spec: SystemSpec) -> dict:
    """Extracts the measurements from a batch result. **Every value carries a
    seed-ensemble error.**"""
    jobs = batch["jobs"]
    if not jobs:
        raise ValueError("0 successful runs — there is nothing to measure")

    d = derive(spec)
    dim = reduced.dim
    var_c = aggregate_seeds([j["var_c"] for j in jobs])
    kT_conf = aggregate_seeds([j["kT_conf"] for j in jobs])

    plateaus, taus, r2s = [], [], []
    for j in jobs:
        r = load_run(Path(j["outdir"]))
        lags_tau = r["lags_steps"] * reduced.dt_star
        f = fit_msd(lags_tau, r["msd"], dim)
        plateaus.append(f.plateau)
        taus.append(f.tau)
        r2s.append(f.r_squared)
    pl, tau = aggregate_seeds(plateaus), aggregate_seeds(taus)

    # reduced → SI
    var_si = d["var_per_component_si"]
    nm2 = 1e18
    meas = {
        "var_x_star": Measurement("var_x_star", var_c.mean, var_c.se,
                                  method="per-component variance of independent "
                                         "snapshots, seed ensemble",
                                  n_samples=var_c.n_seeds, spread=var_c.spread),
        "var_x_nm2": Measurement("var_x_nm2", var_c.mean * var_si * nm2,
                                 var_c.se * var_si * nm2, unit="nm^2",
                                 method="⟨x*²⟩ × kT/k", n_samples=var_c.n_seeds),
        "var_r_nm2": Measurement("var_r_nm2", dim * var_c.mean * var_si * nm2,
                                 dim * var_c.se * var_si * nm2, unit="nm^2",
                                 method="d × ⟨x²⟩ (isotropic)",
                                 n_samples=var_c.n_seeds),
        "msd_plateau_star": Measurement("msd_plateau_star", pl.mean, pl.se,
                                        method="plateau of the MSD fit",
                                        n_samples=pl.n_seeds, spread=pl.spread),
        "tau_star": Measurement("tau_star", tau.mean, tau.se,
                                method="relaxation time of the MSD fit",
                                n_samples=tau.n_seeds,
                                spread=tau.spread),
        "tau_trap_ms": Measurement("tau_trap_ms",
                                   tau.mean * d["tau_trap_si"] * 1e3,
                                   tau.se * d["tau_trap_si"] * 1e3, unit="ms",
                                   method="τ* × τ_trap", n_samples=tau.n_seeds),
        "msd_r_squared": Measurement("msd_r_squared", float(np.mean(r2s)),
                                     float(np.std(r2s, ddof=1) / np.sqrt(len(r2s)))
                                     if len(r2s) > 1 else None,
                                     method="R² of the single-exponential fit",
                                     n_samples=len(r2s)),
        "kT_conf_star": Measurement("kT_conf_star", kT_conf.mean, kT_conf.se,
                                    method="configurational temperature "
                                           "⟨|∇U|²⟩/⟨∇²U⟩",
                                    n_samples=kT_conf.n_seeds),
    }
    # self-consistency: plateau = 2d⟨x*²⟩ — the ratio of two independent paths. A
    # strong check that was not in the prediction
    meas["plateau_over_2d_var"] = Measurement(
        "plateau_over_2d_var", pl.mean / (2 * dim * var_c.mean),
        pl.se / (2 * dim * var_c.mean),
        method="MSD (time series) / variance (snapshot) — independent-path comparison",
        n_samples=pl.n_seeds)
    return meas


# =============================================================================
# run
# =============================================================================
def cmd_run(args) -> int:
    spec_path = Path(args.spec)
    spec = SystemSpec.load(spec_path)
    _runner_for(spec)
    pol = load_policy(args.policy)

    print(f"# S3 SPEC — {spec_path}")
    print(f"  card {spec.card}")
    rep = validate_spec(spec)
    for c in rep.failed():
        print(f"  ❌ {c.name}: {c.detail}")
    for p in rep.problems:
        print(f"  ⚠️  {p}")
    if not rep.ok and not args.force:
        return _fail("the S3 gates did not pass. Fix it, or `--force` "
                     "(not recommended).")

    print("\n# S4 NONDIM")
    reduced = reduce_spec(spec, policy=pol)
    print(f"  scales {reduced.scales.origin}")
    print(f"  dt* = {reduced.dt_star:g}  (dominant constraint: "
          f"{reduced.dt_dominant})")
    errs = __import__("simbot.nondim", fromlist=["roundtrip_errors"]) \
        .roundtrip_errors(spec, reduced)
    worst = max(errs.values())
    print(f"  max round-trip error {worst:.2e}  (gate < 1e-12)"
          + ("" if worst < 1e-12 else "  ❌"))
    if worst >= 1e-12:
        return _fail("S4 round-trip gate violated — a scale convention was broken.")

    n_seeds = int(spec.numerics.n_seeds.si) if spec.numerics.n_seeds \
        else pol.seeds_default
    if n_seeds < int(pol.get("seeds.minimum", 4)) and not args.force:
        return _fail(f"{n_seeds} seeds < the minimum "
                     f"{pol.get('seeds.minimum')}. A production run without error "
                     f"bars is forbidden (CLAUDE.md).")

    print("\n# COST")
    wall_est, k, over = _print_cost(reduced, n_seeds, pol)
    if over and not args.force:
        return _fail(f"the estimate is over budget — **report without running** "
                     f"(run_policy §5 on_exceed). `--force` overrides.")

    # --- the run directory ---
    run_id = args.run_id or io.new_run_id(spec_path.parent.name or spec.card,
                                          spec.hash(), date.today())
    rd = io.RunDir.create(args.runs_root, run_id)
    rd.write("spec", spec.to_yaml())
    rd.write("reduced", _reduced_yaml(reduced))
    rd.write("nondim", __import__("simbot.nondim", fromlist=["nondim_table"])
             .nondim_table(spec, reduced))
    print(f"\n# run directory {rd.path}")

    # --- sealing the prediction ---
    prediction = None
    if args.prediction:
        src = Path(args.prediction)
        rd.write("prediction", src.read_text(encoding="utf-8"))
        prediction = load_prediction(src)
        seal = io.write_seal(rd)
        print(f"  sealed {seal.name} — {len(prediction.items)} predictions")
    else:
        print("  ⚠️  no prediction file (`--prediction`) — proceeding with no S7 "
              "comparison. With nothing to seal there is nothing preventing "
              "post-hoc rationalisation.")

    # --- the tier ladder ---
    ladder_note = _tier_ladder_check(spec, Path(args.runs_root), pol, args.force)
    if ladder_note:
        print(f"  {ladder_note}")

    # --- S5 execution ---
    print(f"\n# S5 RUN — {n_seeds} seeds, concurrency {k}")
    cfgs = _trap_configs(spec, reduced, list(range(spec.numerics.seed_base,
                                                   spec.numerics.seed_base + n_seeds)))
    t0 = time.perf_counter()

    def progress(rec, i, n):
        print(f"  [{i:2d}/{n}] {rec['label']:<20} "
              f"<x*²>={rec['var_c']:.5f}±{rec['var_c_se']:.5f} "
              f"kT_conf={rec['kT_conf']:.5f}  {rec['wall_s']:.2f}s", flush=True)

    batch = run_trap_batch(cfgs, rd.raw, concurrency=k, on_done=progress)
    wall = time.perf_counter() - t0
    rd.write_json("manifest", io.build_manifest(
        run_id=run_id, spec_hash=spec.hash(),
        seed=[c.seed for c in cfgs], rundir=rd,
        extra={"batch": {kk: vv for kk, vv in batch.items() if kk != "jobs"},
               "wall_s_measured": round(wall, 3),
               "wall_s_estimated": round(wall_est, 3)}))
    print(f"  measured wall {wall:.2f} s (estimated {wall_est:.2f} s, "
          f"ratio {wall / wall_est:.2f})")
    if batch["failed"]:
        print(f"  ❌ {len(batch['failed'])} runs failed — the seed count dropped. "
              f"See 05_run_manifest.json")
        for f in batch["failed"]:
            print(f"     {f['label']}: {f['stderr'].strip().splitlines()[-1:]}")
    for j in batch["jobs"]:
        if not j["guards"]["finite"]:
            print(f"  ❌ guard violation {j['label']}: {j['guards']['failures']}")

    # --- S6 figures ---
    print("\n# S6 VISUALIZE")
    figset = _make_figures(rd, batch, reduced, spec)
    rd.write("figures", figset.figures_md())
    print(f"  {len(figset.records)} figures → {rd.figs}")
    for r in figset.records:
        print(f"    {r.name:<26} {r.shows[:60]}")
    for n, why in sorted(figset.skipped.items()):
        print(f"    (skipped) {n:<17} {why[:60]}")

    # --- S7 measurement + verdict ---
    print("\n# S7 VALIDATE")
    meas = measure_trap_batch(batch, reduced, spec)
    rd.write_json("metrics", {k2: vars(v) for k2, v in meas.items()})
    for name, m in meas.items():
        err = f" ± {m.stat_err:.4g}" if m.stat_err is not None else " (no error)"
        print(f"  {name:<22} {m.value:>14.6g}{err}  {m.unit}")

    vrep = None
    if prediction is not None:
        vrep = validate_run(prediction, meas, rundir=rd,
                            power_threshold=args.power_threshold)
        print()
        print(vrep.table())
        print()
        print(vrep.reasons())
        for p in vrep.problems:
            print(f"  ⚠️  {p}")

    # --- S8 report ---
    write_report(rd, ReportInputs(spec=spec, spec_report=rep, reduced=reduced,
                                  validation=vrep, manifest=rd.read_json("manifest"),
                                  figures=figset.captions,
                                  wall_s=wall, n_runs=len(batch["jobs"])))
    print(f"\n# S8 REPORT → {rd.file('report')}")
    if vrep is not None:
        print(f"  verdict (proposed) {vrep.verdict_overall} · confirmed_by: null")
    print("  ★ Confirming the verdict is a human's job. The conclusion narrative "
          "(08_conclusion.md) is written by the agent.")
    return 0


def _make_figures(rd: io.RunDir, batch: dict, reduced, spec: SystemSpec):
    """Batch result → the S6 figure set. Branches per runner (only trap exists)."""
    d = derive(spec)
    runs = {j["label"]: load_run(Path(j["outdir"])) for j in batch["jobs"]}
    first = runs[sorted(runs)[0]]
    frame_interval = int(first["lags_steps"][0])      # msd_stride

    # several dt* in the same batch makes the bias figure drawable
    by_dt: dict[float, tuple[float, float]] = {}
    groups: dict[float, list[float]] = {}
    for j in batch["jobs"]:
        groups.setdefault(float(j["config"]["dt_star"]), []).append(j["var_c"])
    for dt, vals in groups.items():
        agg = aggregate_seeds(vals)
        by_dt[dt] = (agg.mean, agg.se if agg.se == agg.se else 0.0)

    return trap_diagnostics(
        rd.figs, runs, dim=reduced.dim, dt_star=reduced.dt_star,
        sample_interval_steps=reduced.sample_interval_steps,
        frame_interval_steps=frame_interval, sigma_star=reduced.sigma_star,
        tau_trap_si=d["tau_trap_si"], l_trap_si=d["l_trap_si"],
        by_dt=by_dt, has_pair=bool(spec.pair))


def _reduced_yaml(reduced) -> str:
    import yaml
    d = {k: v for k, v in vars(reduced).items() if k != "scales"}
    d["scales"] = {"length_si": reduced.scales.length_si,
                   "energy_si": reduced.scales.energy_si,
                   "time_si": reduced.scales.time_si,
                   "origin": reduced.scales.origin}
    d["inverse"] = reduced.inverse
    return yaml.dump(d, sort_keys=False, allow_unicode=True, default_flow_style=False)


def _tier_ladder_check(spec: SystemSpec, runs_root: Path, pol, force: bool) -> str:
    """`smoke → pilot → explore` cannot be skipped on a system seen for the first
    time.

    "seen for the first time" is judged by **whether a prior run of the same card
    exists.**
    """
    ladder = pol.tier_ladder()
    if not runs_root.exists():
        prior = []
    else:
        prior = [p for p in runs_root.glob("*/03_spec.yaml")
                 if f"card: {spec.card}" in p.read_text(encoding="utf-8")[:400]]
    if prior:
        return f"tier ladder: {len(prior)} prior runs of this card → ladder ok"
    return (f"⚠️  this is the card's first run. Policy does not allow skipping the "
            f"ladder {ladder} — this run becomes its first rung. Do not cite the "
            f"result as production"
            + ("" if not force else " (--force)"))


# =============================================================================
# resume
# =============================================================================
def cmd_resume(args) -> int:
    rd = io.RunDir(Path(args.rundir))
    if not rd.path.exists():
        return _fail(f"{rd.path} does not exist")
    done = rd.completed_stages()
    print(f"# resume {rd.run_id}")
    print(f"  {len(done)} completed stages: {done}")

    if not rd.exists("spec"):
        return _fail("no 03_spec.yaml — there is no basis to resume from")
    spec = SystemSpec.load(rd.file("spec"))

    if rd.exists("seal"):
        v = io.verify_seal(rd)
        print(f"  seal: {'✅' if v.ok else '❌'} {v.summary()}")
        if not v.ok:
            return _fail("seal violation — resuming makes the verification "
                         "meaningless. Revert the prediction or start a new run.")

    raw_dirs = sorted(p for p in rd.raw.glob("*") if (p / "samples.npz").exists())
    print(f"  {len(raw_dirs)} completed runs: {[p.name for p in raw_dirs]}")
    if not raw_dirs:
        return _fail("0 completed runs — start fresh with `cli.py run`")

    pol = load_policy(args.policy)
    reduced = reduce_spec(spec, policy=pol)
    batch = {"jobs": [_job_from_dir(p) for p in raw_dirs], "failed": [],
             "n_requested": len(raw_dirs), "batch_wall_s": 0.0, "concurrency": 0}

    print("\n# S6 (figures regenerated — the trajectories exist, so the run is not "
          "repeated)")
    figset = _make_figures(rd, batch, reduced, spec)
    rd.write("figures", figset.figures_md())
    print(f"  {len(figset.records)} figures · {len(figset.skipped)} skipped")

    print("\n# S7 (only what needs recomputing)")
    meas = measure_trap_batch(batch, reduced, spec)
    rd.write_json("metrics", {k: vars(v) for k, v in meas.items()})
    for name, m in meas.items():
        err = f" ± {m.stat_err:.4g}" if m.stat_err is not None else ""
        print(f"  {name:<22} {m.value:>14.6g}{err}")

    vrep = None
    if args.prediction or rd.exists("prediction_json"):
        src = Path(args.prediction) if args.prediction else rd.file("prediction_json")
        vrep = validate_run(load_prediction(src), meas, rundir=rd)
        print()
        print(vrep.table())

    man = rd.read_json("manifest") if rd.exists("manifest") else None
    write_report(rd, ReportInputs(spec=spec, reduced=reduced, validation=vrep,
                                  manifest=man, figures=figset.captions,
                                  n_runs=len(raw_dirs)))
    print(f"\n# REPORT → {rd.file('report')}")
    return 0


def _job_from_dir(p: Path) -> dict:
    man = json.loads((p / "manifest.json").read_text())
    d = load_run(p)
    var = d["indep_var"]
    n = var.size
    return {"label": p.name, "outdir": str(p),
            "var_c": float(var.mean()),
            "var_c_se": float(var.std(ddof=1) / np.sqrt(n)) if n > 1 else float("nan"),
            "kT_conf": float(d["indep_kT"].mean()),
            "kT_conf_se": 0.0, "wall_s": man["manifest"].get("wall_s", 0.0),
            "guards": man.get("guards", {"finite": True, "failures": []}),
            "n_snap": n, "config": man["config"]}


# =============================================================================
# converge
# =============================================================================
def cmd_converge(args) -> int:
    """Is the answer the same when `dt`, `N` and the initial condition are shaken.

    Convergence is judged **relative to the statistical error.** An absolute
    threshold unfairly rejects a precise run and unfairly passes a coarse one.
    """
    spec = SystemSpec.load(args.spec)
    _runner_for(spec)
    pol = load_policy(args.policy)
    base = reduce_spec(spec, policy=pol)
    n_seeds = int(spec.numerics.n_seeds.si) if spec.numerics.n_seeds \
        else pol.seeds_default

    variants = {
        "base": dict(),
        "dt_half": dict(dt_star=base.dt_star / 2),
        "dt_double": dict(dt_star=base.dt_star * 2),
        "N_double": dict(n_particles=base.n_particles * 2),
        "seed_shift": dict(seed_offset=100),
    }
    if args.only:
        variants = {k: v for k, v in variants.items() if k in set(args.only) | {"base"}}

    print(f"# converge — {len(variants)} variants × {n_seeds} seeds")
    k = min(n_seeds, pol.concurrency())
    total = sum(estimate_wall_time_s(
        v.get("n_particles", base.n_particles),
        int(round((base.equil_steps + base.prod_steps)
                  * base.dt_star / v.get("dt_star", base.dt_star))), k)
        for v in variants.values())
    print(f"  estimated total wall ≈ {total:.1f} s (budget "
          f"{pol.get('budget.wall_time_per_batch_s', 1800):g} s)")
    if total > pol.get("budget.wall_time_per_batch_s", 1800) and not args.force:
        return _fail("the batch estimate is over budget — not running.")

    outroot = Path(args.out or (REPO / "runs" / f"converge_{spec.hash()[:6]}"))
    results = {}
    for name, ov in variants.items():
        dt = ov.get("dt_star", base.dt_star)
        npart = ov.get("n_particles", base.n_particles)
        off = ov.get("seed_offset", 0)
        cfgs = [TrapRunConfig(
            dim=base.dim, n_particles=npart, dt_star=dt,
            equil_tau=base.equil_steps * base.dt_star,
            prod_tau=base.prod_steps * base.dt_star,
            sample_interval_tau=base.sample_interval_steps * base.dt_star,
            box_over_l_trap=base.box_star[0], k_star=base.k_star or 1.0,
            seed=spec.numerics.seed_base + off + i,
            label=f"{name}_s{spec.numerics.seed_base + off + i}")
            for i in range(n_seeds)]
        print(f"\n  [{name}] dt*={dt:g} N={npart} seed+{off}")
        batch = run_trap_batch(cfgs, outroot / name, concurrency=k)
        if not batch["jobs"]:
            print("    ❌ all failed")
            continue
        agg = aggregate_seeds([j["var_c"] for j in batch["jobs"]])
        taus = []
        for j in batch["jobs"]:
            r = load_run(Path(j["outdir"]))
            taus.append(fit_msd(r["lags_steps"] * dt, r["msd"], base.dim).tau)
        tagg = aggregate_seeds(taus)
        results[name] = {"var_x_star": agg, "tau_star": tagg, "dt_star": dt,
                         "n_particles": npart}
        print(f"    <x*²> = {agg.mean:.6f} ± {agg.se:.6f} · "
              f"τ* = {tagg.mean:.6f} ± {tagg.se:.6f}")

    if "base" not in results:
        return _fail("the baseline run failed — there is nothing to compare against")

    print("\n# convergence verdict (vs the baseline, in units of statistical error)")
    print(f"{'variant':<12} {'quantity':<12} {'baseline':>11} {'variant':>11} "
          f"{'diff':>10} {'combSE':>10} {'σ':>7}  verdict")
    b = results["base"]
    verdicts = []
    for name, r in results.items():
        if name == "base":
            continue
        for q in ("var_x_star", "tau_star"):
            x0, x1 = b[q], r[q]
            diff = x1.mean - x0.mean
            se = float(np.hypot(x0.se, x1.se))
            sig = abs(diff) / se if se > 0 else float("inf")
            # ★ Within 3σ is "indistinguishable" — an observation that they do not
            #   differ, not a proof that they are equal
            ok = sig < 3.0
            verdicts.append(ok)
            print(f"{name:<12} {q:<12} {x0.mean:>11.6f} {x1.mean:>11.6f} "
                  f"{diff:>+10.6f} {se:>10.6f} {sig:>6.2f}σ  "
                  f"{'indistinguishable' if ok else '❌ significant difference'}")

    outroot.mkdir(parents=True, exist_ok=True)
    (outroot / "converge.json").write_text(json.dumps(
        {n: {"var_x_star": vars(r["var_x_star"]), "tau_star": vars(r["tau_star"]),
             "dt_star": r["dt_star"], "n_particles": r["n_particles"]}
         for n, r in results.items()}, indent=2, default=float), encoding="utf-8")
    print(f"\n  → {outroot / 'converge.json'}")
    print("\n★ 'indistinguishable' is **not a proof of equality** — it means this "
          "statistical error cannot see the difference.\n  A real bias (the EM bias, "
          "say) has to be measured separately where there is power to do so.")
    return 0 if all(verdicts) else 1


# =============================================================================
# params
# =============================================================================
def cmd_params(args) -> int:
    """Several runs' parameters **side by side**. A default nobody picked gets ⚠."""
    root = Path(args.path)
    specs: dict[str, SystemSpec] = {}
    for p in sorted(root.glob("*/03_spec.yaml")):
        try:
            specs[p.parent.name] = SystemSpec.load(p)
        except Exception as e:
            print(f"  (skipped {p.parent.name}: {e})", file=sys.stderr)
    for p in sorted(root.glob("*.yaml")):
        try:
            specs[p.stem] = SystemSpec.load(p)
        except Exception:
            pass
    if not specs:
        return _fail(f"no spec found under {root} "
                     f"(looking for: */03_spec.yaml or *.yaml)")

    from simbot.spec import _iter_quantities
    rows: dict[str, dict[str, str]] = {}
    prov: dict[str, set[str]] = {}
    for name, sp in specs.items():
        for path, q in _iter_quantities(sp):
            rows.setdefault(path, {})[name] = _short_value(q.value)
            prov.setdefault(path, set()).add(q.provenance)

    names = list(specs)
    print(f"# params — {len(names)} specs × {len(rows)} fields\n")
    head = f"| {'field':<34} | " \
           + " | ".join(f"{n[:22]:<22}" for n in names) + " | prov |"
    print(head)
    print("|" + "-" * 36 + "|" + "|".join(["-" * 24] * len(names)) + "|------|")
    for path in sorted(rows):
        vals = [rows[path].get(n, "—") for n in names]
        differs = len(set(vals)) > 1
        # ⚠ = a value nobody picked (assumed, and identical across every spec)
        kinds = prov[path]
        mark = "⚠" if (kinds <= {"assumed"} and not differs) else ""
        star = "**" if differs else ""
        print(f"| {star}{path:<32}{star} | "
              + " | ".join(f"{v:<22}" for v in vals) + f" | {mark}{'/'.join(sorted(kinds))} |")
    print("\n⚠ = the provenance is only `assumed` and the value is the same in every "
          "spec — **a default nobody picked.**")
    print("**bold** = a field whose value differs between specs.")
    return 0


def _short_value(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    s = str(v)
    return s if len(s) <= 22 else s[:19] + "…"


# =============================================================================
# calibrate
# =============================================================================
# The kernel the policy constant `Λ` was measured on. A value measured on a different
# kernel is not directly comparable.
#   Basis: knowledge/wiki/findings/local-cpu-parallelism.md §2 (the measurement header)
BASELINE_KERNEL = "3D WCA + Brownian, phi=0.30, Cell nlist(buffer 0.3), dt=1e-4"
CALIBRATE_KERNEL = "2D harmonic trap + Brownian, no pair interaction, no nlist"


def cmd_calibrate(args) -> int:
    """Measures this machine's throughput.

    ★ **A different kernel does not overwrite the constant.** The policy constant
      `Λ = 6.3e6` was measured on 3D WCA + Cell nlist, and the only kernel this
      command can run is a harmonic trap with **no** pair interaction (`simbot.run`
      has no other runner). A kernel with the neighbour-list work removed being
      faster is expected, and updating the global constant with that value makes
      **the cost estimate for WCA systems optimistically wrong.**
      So it is reported as a per-kernel multiple, and no global-constant update is
      proposed.
    """
    pol = load_policy(args.policy)
    claimed = pol.get("hardware.throughput_particle_steps_per_s")
    print(f"# calibrate")
    print(f"  policy constant Λ = {claimed:.3g} particle·steps/s")
    print(f"    measured on: {BASELINE_KERNEL}")
    print(f"  measured here on: {CALIBRATE_KERNEL}")
    print(f"  ⚠ the throughput model was measured at N ≥ 500. At small N the "
          f"overhead dominates.\n")

    cases = [(int(n), int(s)) for n, s in (args.cases or [(1000, 20000),
                                                          (4000, 20000),
                                                          (8000, 10000)])]
    print(f"{'N':>7} {'steps':>8} {'wall [s]':>10} {'part·steps/s':>14} "
          f"{'vs Λ':>8}")
    measured = []
    for n, steps in cases:
        cfg = TrapRunConfig(n_particles=n, dt_star=5e-3,
                            equil_tau=0.0, prod_tau=steps * 5e-3,
                            sample_interval_tau=2.0, label=f"cal_N{n}")
        from simbot.run import run_trap
        r = run_trap(cfg)
        thr = n * steps / r.wall_s
        measured.append(thr)
        print(f"{n:>7} {steps:>8} {r.wall_s:>10.3f} {thr:>14.4g} "
              f"{thr / claimed:>7.2f}×")

    med = float(np.median(measured))
    ratio = med / claimed
    print(f"\n  median {med:.4g} particle·steps/s · vs Λ {ratio:.3f}×")
    print(f"\n  ⇒ **trap-kernel multiple = {ratio:.2f}** (no pair interaction / the "
          f"WCA baseline)")
    print(f"     This is the value to use for trap-system cost estimates. The global "
          f"Λ must not be")
    print(f"     replaced with it — WCA estimates would be {ratio:.2f}x "
          f"optimistically wrong.")
    print(f"\n  Note: the first run's 03_spec.yaml set `pair_cost_multiplier: 1.0`")
    print(f"        (no pair interaction, so the baseline as is). This measurement "
          f"shows it")
    print(f"        is not that but {ratio:.2f} — the trap system is faster than the "
          f"baseline.")

    if abs(ratio - 1) > 0.2:
        print(f"\n  ⚠️  {abs(ratio - 1):.0%} away from the baseline. Telling whether "
              f"that is **a kernel difference or\n"
              f"     a change in the machine's state** requires measuring on the "
              f"baseline's own kernel\n"
              f"     (3D WCA φ=0.30), and that runner does not exist yet "
              f"(`simbot/run.py` has only\n"
              f"     the trap runner).\n"
              f"     ⇒ What can be concluded now: **estimating trap-system cost with "
              f"Λ is {ratio:.2f}x conservative.**")
    else:
        print("\n  ✅ Within 20 % of the baseline despite the different kernel — Λ "
              "can be used as is.")
    print("\n  ★ The code never edits the policy file automatically — a measured "
          "constant is confirmed and committed by a human.")
    return 0


# =============================================================================
# Entry point
# =============================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="cli.py", description=__doc__,   # noqa: E501
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--policy", default=None, help="path to run_policy.yaml")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("run", help="one spec → REPORT.md")
    q.add_argument("spec")
    q.add_argument("--prediction", default=None,
                   help="the S2 prediction YAML (what gets sealed)")
    q.add_argument("--runs-root", default=str(REPO / "runs"))
    q.add_argument("--run-id", default=None)
    q.add_argument("--power-threshold", type=float, default=1.0)
    q.add_argument("--force", action="store_true",
                   help="ignore the gates and the budget")
    q.set_defaults(func=cmd_run)

    q = sub.add_parser("resume", help="pick up a dead run")
    q.add_argument("rundir")
    q.add_argument("--prediction", default=None)
    q.set_defaults(func=cmd_resume)

    q = sub.add_parser("converge",
                       help="is the answer the same when dt, N and the initial "
                            "condition are shaken")
    q.add_argument("spec")
    q.add_argument("--only", nargs="+", default=None,
                   choices=["dt_half", "dt_double", "N_double", "seed_shift"])
    q.add_argument("--out", default=None)
    q.add_argument("--force", action="store_true")
    q.set_defaults(func=cmd_converge)

    q = sub.add_parser("params", help="several runs' parameters side by side")
    q.add_argument("--path", default=str(REPO / "runs"))
    q.set_defaults(func=cmd_params)

    q = sub.add_parser("calibrate", help="measure this machine's throughput")
    q.add_argument("--cases", nargs="+", type=int, default=None,
                   metavar="N STEPS", help="N steps pairs, concatenated")
    q.set_defaults(func=cmd_calibrate)
    return p


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    if a.cmd == "calibrate" and a.cases:
        it = iter(a.cases)
        a.cases = list(zip(it, it))
    return a.func(a)


if __name__ == "__main__":
    raise SystemExit(main())
