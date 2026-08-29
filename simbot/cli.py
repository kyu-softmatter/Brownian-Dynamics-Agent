#!/usr/bin/env python
"""Simulation Bot CLI — S3 SPEC → S8 REPORT. LLM 0줄.

**One command surface (merged 2026-08-29).** This module holds the S5-S8 half;
`bdbot.cli` holds L0-L4 and now exposes these same commands under `pipeline`, so
there is one entry point rather than two engines with separate front ends:

    python -m bdbot.cli pipeline run <spec.yaml>     the merged entry point
    python cli.py run <spec.yaml>                    the root shim, unchanged
    python -m simbot.cli run <spec.yaml>             this module directly

All three dispatch to the same `main()`. `bdbot.cli` imports this **lazily**, inside
the command handler, so its front end does not pull in matplotlib (`simbot.viz`) —
that property is load-bearing and `tests/test_bdbot_lazy_api.py` guards it.

    python cli.py run <spec.yaml>            spec 하나 → REPORT.md
    python cli.py resume runs/<id>           죽은 런 이어받기 (완료 단계 재계산 안 함)
    python cli.py converge <spec.yaml>       dt·N·초기조건을 흔들어도 답이 같은가
    python cli.py params [--path runs]       여러 런의 파라미터를 가로로
    python cli.py calibrate                  이 기계의 처리량 실측

## 이 CLI 가 하지 않는 것

- **판정을 확정하지 않는다.** `confirmed_by` 는 사람만 쓴다.
- **예산을 넘겨 실행하지 않는다.** 초과 예상 시 실행 전에 보고하고 멈춘다.
- **모르는 카드를 억지로 돌리지 않는다.** 러너가 없으면 그렇게 말한다.

현재 러너가 있는 카드는 `passive-sphere--harmonic-trap` 하나다. 다른 카드는
`simbot.run` 에 러너를 만들어야 한다 — 조용히 트랩 러너로 돌리지 않는다.
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
# 공통
# =============================================================================
def _fail(msg: str, code: int = 2) -> int:
    print(f"\n⛔ {msg}", file=sys.stderr)
    return code


def _runner_for(spec: SystemSpec) -> str:
    kind = RUNNERS.get(spec.card)
    if kind is None:
        raise SystemExit(
            f"⛔ 카드 {spec.card!r} 에는 러너가 없다.\n"
            f"   러너가 있는 카드: {sorted(RUNNERS)}\n"
            f"   다른 카드를 트랩 러너로 돌리면 조용히 틀린 계를 계산한다 —\n"
            f"   `simbot/run.py` 에 러너를 만들고 RUNNERS 에 등록할 것.")
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
    print(f"  N={reduced.n_particles} · steps/seed={steps} · 시드 {n_seeds} "
          f"· 동시 {k} (효율 {pol.efficiency(k):.3f})")
    print(f"  추정 wall ≈ {wall:.2f} s  (예산 {pol.wall_budget_s:g} s)"
          + ("  ⚠️ **초과**" if over else ""))
    return wall, k, over


# =============================================================================
# 측정 — 배치 결과 → Measurement
# =============================================================================
def measure_trap_batch(batch: dict, reduced, spec: SystemSpec) -> dict:
    """배치 결과에서 측정값을 뽑는다. **모든 값에 시드 앙상블 오차가 붙는다.**"""
    jobs = batch["jobs"]
    if not jobs:
        raise ValueError("성공한 런이 0개다 — 측정할 것이 없다")

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

    # 무차원 → SI
    var_si = d["var_per_component_si"]
    nm2 = 1e18
    meas = {
        "var_x_star": Measurement("var_x_star", var_c.mean, var_c.se,
                                  method="독립 스냅샷 성분분산, 시드 앙상블",
                                  n_samples=var_c.n_seeds, spread=var_c.spread),
        "var_x_nm2": Measurement("var_x_nm2", var_c.mean * var_si * nm2,
                                 var_c.se * var_si * nm2, unit="nm^2",
                                 method="⟨x*²⟩ × kT/k", n_samples=var_c.n_seeds),
        "var_r_nm2": Measurement("var_r_nm2", dim * var_c.mean * var_si * nm2,
                                 dim * var_c.se * var_si * nm2, unit="nm^2",
                                 method="d × ⟨x²⟩ (등방)", n_samples=var_c.n_seeds),
        "msd_plateau_star": Measurement("msd_plateau_star", pl.mean, pl.se,
                                        method="MSD 피팅 plateau",
                                        n_samples=pl.n_seeds, spread=pl.spread),
        "tau_star": Measurement("tau_star", tau.mean, tau.se,
                                method="MSD 피팅 완화시간", n_samples=tau.n_seeds,
                                spread=tau.spread),
        "tau_trap_ms": Measurement("tau_trap_ms",
                                   tau.mean * d["tau_trap_si"] * 1e3,
                                   tau.se * d["tau_trap_si"] * 1e3, unit="ms",
                                   method="τ* × τ_trap", n_samples=tau.n_seeds),
        "msd_r_squared": Measurement("msd_r_squared", float(np.mean(r2s)),
                                     float(np.std(r2s, ddof=1) / np.sqrt(len(r2s)))
                                     if len(r2s) > 1 else None,
                                     method="단일지수 피팅 R²", n_samples=len(r2s)),
        "kT_conf_star": Measurement("kT_conf_star", kT_conf.mean, kT_conf.se,
                                    method="배위 온도 ⟨|∇U|²⟩/⟨∇²U⟩",
                                    n_samples=kT_conf.n_seeds),
    }
    # 자기일관성: plateau = 2d⟨x*²⟩ — 독립 경로 둘의 비. 예측에 없던 강한 확인
    meas["plateau_over_2d_var"] = Measurement(
        "plateau_over_2d_var", pl.mean / (2 * dim * var_c.mean),
        pl.se / (2 * dim * var_c.mean),
        method="MSD(시계열) / 분산(스냅샷) — 독립 경로 대조", n_samples=pl.n_seeds)
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
    print(f"  카드 {spec.card}")
    rep = validate_spec(spec)
    for c in rep.failed():
        print(f"  ❌ {c.name}: {c.detail}")
    for p in rep.problems:
        print(f"  ⚠️  {p}")
    if not rep.ok and not args.force:
        return _fail("S3 게이트 미통과. 고치거나 `--force` (권장하지 않음).")

    print("\n# S4 NONDIM")
    reduced = reduce_spec(spec, policy=pol)
    print(f"  척도 {reduced.scales.origin}")
    print(f"  dt* = {reduced.dt_star:g}  (지배 제약: {reduced.dt_dominant})")
    errs = __import__("simbot.nondim", fromlist=["roundtrip_errors"]) \
        .roundtrip_errors(spec, reduced)
    worst = max(errs.values())
    print(f"  왕복오차 최대 {worst:.2e}  (게이트 < 1e-12)"
          + ("" if worst < 1e-12 else "  ❌"))
    if worst >= 1e-12:
        return _fail("S4 왕복 게이트 위반 — 척도 규약이 어긋났다.")

    n_seeds = int(spec.numerics.n_seeds.si) if spec.numerics.n_seeds \
        else pol.seeds_default
    if n_seeds < int(pol.get("seeds.minimum", 4)) and not args.force:
        return _fail(f"시드 {n_seeds}개 < 최소 {pol.get('seeds.minimum')}개. "
                     f"오차 막대 없는 프로덕션 런은 금지다 (CLAUDE.md).")

    print("\n# 비용")
    wall_est, k, over = _print_cost(reduced, n_seeds, pol)
    if over and not args.force:
        return _fail(f"예산 초과 예상 — **실행하지 않고 보고한다** "
                     f"(run_policy §5 on_exceed). `--force` 로 강행 가능.")

    # --- run 디렉터리 ---
    run_id = args.run_id or io.new_run_id(spec_path.parent.name or spec.card,
                                          spec.hash(), date.today())
    rd = io.RunDir.create(args.runs_root, run_id)
    rd.write("spec", spec.to_yaml())
    rd.write("reduced", _reduced_yaml(reduced))
    rd.write("nondim", __import__("simbot.nondim", fromlist=["nondim_table"])
             .nondim_table(spec, reduced))
    print(f"\n# run 디렉터리 {rd.path}")

    # --- 예측 봉인 ---
    prediction = None
    if args.prediction:
        src = Path(args.prediction)
        rd.write("prediction", src.read_text(encoding="utf-8"))
        prediction = load_prediction(src)
        seal = io.write_seal(rd)
        print(f"  봉인 {seal.name} — 예측 {len(prediction.items)}개")
    else:
        print("  ⚠️  예측 파일이 없다 (`--prediction`) — S7 대조 없이 진행한다. "
              "봉인할 것이 없으므로 사후합리화를 막을 장치가 없다.")

    # --- 티어 사다리 ---
    ladder_note = _tier_ladder_check(spec, Path(args.runs_root), pol, args.force)
    if ladder_note:
        print(f"  {ladder_note}")

    # --- S5 실행 ---
    print(f"\n# S5 RUN — 시드 {n_seeds}개, 동시 {k}")
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
    print(f"  실측 wall {wall:.2f} s (추정 {wall_est:.2f} s, "
          f"비 {wall / wall_est:.2f})")
    if batch["failed"]:
        print(f"  ❌ 실패한 런 {len(batch['failed'])}개 — 시드 수가 줄었다. "
              f"05_run_manifest.json 참조")
        for f in batch["failed"]:
            print(f"     {f['label']}: {f['stderr'].strip().splitlines()[-1:]}")
    for j in batch["jobs"]:
        if not j["guards"]["finite"]:
            print(f"  ❌ 가드 위반 {j['label']}: {j['guards']['failures']}")

    # --- S6 그림 ---
    print("\n# S6 VISUALIZE")
    figset = _make_figures(rd, batch, reduced, spec)
    rd.write("figures", figset.figures_md())
    print(f"  그림 {len(figset.records)}장 → {rd.figs}")
    for r in figset.records:
        print(f"    {r.name:<26} {r.shows[:60]}")
    for n, why in sorted(figset.skipped.items()):
        print(f"    (건너뜀) {n:<17} {why[:60]}")

    # --- S7 측정 + 판정 ---
    print("\n# S7 VALIDATE")
    meas = measure_trap_batch(batch, reduced, spec)
    rd.write_json("metrics", {k2: vars(v) for k2, v in meas.items()})
    for name, m in meas.items():
        err = f" ± {m.stat_err:.4g}" if m.stat_err is not None else " (오차 없음)"
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

    # --- S8 리포트 ---
    write_report(rd, ReportInputs(spec=spec, spec_report=rep, reduced=reduced,
                                  validation=vrep, manifest=rd.read_json("manifest"),
                                  figures=figset.captions,
                                  wall_s=wall, n_runs=len(batch["jobs"])))
    print(f"\n# S8 REPORT → {rd.file('report')}")
    if vrep is not None:
        print(f"  판정 (제안) {vrep.verdict_overall} · confirmed_by: null")
    print("  ★ 판정 확정은 사람이 한다. 결론 서술(08_conclusion.md)은 에이전트가 쓴다.")
    return 0


def _make_figures(rd: io.RunDir, batch: dict, reduced, spec: SystemSpec):
    """배치 결과 → S6 그림 세트. 러너별로 분기한다 (트랩만 존재)."""
    d = derive(spec)
    runs = {j["label"]: load_run(Path(j["outdir"])) for j in batch["jobs"]}
    first = runs[sorted(runs)[0]]
    frame_interval = int(first["lags_steps"][0])      # msd_stride

    # 같은 배치 안에 여러 dt* 가 있으면 편향 그림을 그릴 수 있다
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
    """처음 보는 계에서 `smoke → pilot → explore` 를 건너뛸 수 없다.

    '처음 보는가'는 **같은 카드의 이전 런이 있는가**로 판단한다.
    """
    ladder = pol.tier_ladder()
    if not runs_root.exists():
        prior = []
    else:
        prior = [p for p in runs_root.glob("*/03_spec.yaml")
                 if f"card: {spec.card}" in p.read_text(encoding="utf-8")[:400]]
    if prior:
        return f"티어 사다리: 같은 카드의 이전 런 {len(prior)}개 있음 → 사다리 통과"
    return (f"⚠️  이 카드의 첫 런이다. 정책상 사다리 {ladder} 를 건너뛸 수 없다 — "
            f"이 런이 사다리의 첫 칸이 된다. 결과를 프로덕션으로 인용하지 말 것"
            + ("" if not force else " (--force)"))


# =============================================================================
# resume
# =============================================================================
def cmd_resume(args) -> int:
    rd = io.RunDir(Path(args.rundir))
    if not rd.path.exists():
        return _fail(f"{rd.path} 가 없다")
    done = rd.completed_stages()
    print(f"# resume {rd.run_id}")
    print(f"  완료된 단계 {len(done)}개: {done}")

    if not rd.exists("spec"):
        return _fail("03_spec.yaml 이 없다 — 이어받을 근거가 없다")
    spec = SystemSpec.load(rd.file("spec"))

    if rd.exists("seal"):
        v = io.verify_seal(rd)
        print(f"  봉인: {'✅' if v.ok else '❌'} {v.summary()}")
        if not v.ok:
            return _fail("봉인 위반 — 이어받으면 검증이 무의미해진다. "
                         "예측을 되돌리거나 새 run 을 시작할 것.")

    raw_dirs = sorted(p for p in rd.raw.glob("*") if (p / "samples.npz").exists())
    print(f"  완주한 런 {len(raw_dirs)}개: {[p.name for p in raw_dirs]}")
    if not raw_dirs:
        return _fail("완주한 런이 0개 — `cli.py run` 으로 새로 시작할 것")

    pol = load_policy(args.policy)
    reduced = reduce_spec(spec, policy=pol)
    batch = {"jobs": [_job_from_dir(p) for p in raw_dirs], "failed": [],
             "n_requested": len(raw_dirs), "batch_wall_s": 0.0, "concurrency": 0}

    print("\n# S6 (그림 재생성 — 궤적이 있으므로 런은 다시 돌리지 않는다)")
    figset = _make_figures(rd, batch, reduced, spec)
    rd.write("figures", figset.figures_md())
    print(f"  그림 {len(figset.records)}장 · 건너뜀 {len(figset.skipped)}건")

    print("\n# S7 (재계산 대상만)")
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
    """`dt`·`N`·초기조건을 흔들어도 답이 같은가.

    수렴 판정은 **통계오차 대비**로 한다. 절대 문턱을 쓰면 정밀한 런이 부당하게
    기각되고 거친 런이 부당하게 통과한다.
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

    print(f"# converge — {len(variants)} 변형 × 시드 {n_seeds}")
    k = min(n_seeds, pol.concurrency())
    total = sum(estimate_wall_time_s(
        v.get("n_particles", base.n_particles),
        int(round((base.equil_steps + base.prod_steps)
                  * base.dt_star / v.get("dt_star", base.dt_star))), k)
        for v in variants.values())
    print(f"  추정 총 wall ≈ {total:.1f} s (예산 "
          f"{pol.get('budget.wall_time_per_batch_s', 1800):g} s)")
    if total > pol.get("budget.wall_time_per_batch_s", 1800) and not args.force:
        return _fail("배치 예산 초과 예상 — 실행하지 않는다.")

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
            print("    ❌ 전부 실패")
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
        return _fail("기준 런이 실패했다 — 비교 대상이 없다")

    print("\n# 수렴 판정 (기준 대비, 통계오차 단위)")
    print(f"{'변형':<12} {'양':<12} {'기준':>11} {'변형':>11} {'차이':>10} "
          f"{'결합SE':>10} {'σ':>7}  판정")
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
            # ★ 3σ 이내면 "구별되지 않는다" — 같다는 증명이 아니라 다르지 않다는 관찰
            ok = sig < 3.0
            verdicts.append(ok)
            print(f"{name:<12} {q:<12} {x0.mean:>11.6f} {x1.mean:>11.6f} "
                  f"{diff:>+10.6f} {se:>10.6f} {sig:>6.2f}σ  "
                  f"{'구별 안 됨' if ok else '❌ 유의한 차이'}")

    outroot.mkdir(parents=True, exist_ok=True)
    (outroot / "converge.json").write_text(json.dumps(
        {n: {"var_x_star": vars(r["var_x_star"]), "tau_star": vars(r["tau_star"]),
             "dt_star": r["dt_star"], "n_particles": r["n_particles"]}
         for n, r in results.items()}, indent=2, default=float), encoding="utf-8")
    print(f"\n  → {outroot / 'converge.json'}")
    print("\n★ '구별 안 됨'은 **같다는 증명이 아니다** — 이 통계오차로는 차이를 "
          "볼 수 없다는 뜻이다.\n  진짜 편향(EM 편향 등)은 검정력이 있는 곳에서 따로 재야 한다.")
    return 0 if all(verdicts) else 1


# =============================================================================
# params
# =============================================================================
def cmd_params(args) -> int:
    """여러 런의 파라미터를 **가로로**. 아무도 안 고른 기본값에 ⚠ 표시."""
    root = Path(args.path)
    specs: dict[str, SystemSpec] = {}
    for p in sorted(root.glob("*/03_spec.yaml")):
        try:
            specs[p.parent.name] = SystemSpec.load(p)
        except Exception as e:
            print(f"  (건너뜀 {p.parent.name}: {e})", file=sys.stderr)
    for p in sorted(root.glob("*.yaml")):
        try:
            specs[p.stem] = SystemSpec.load(p)
        except Exception:
            pass
    if not specs:
        return _fail(f"{root} 에서 spec 을 찾지 못했다 "
                     f"(찾는 것: */03_spec.yaml 또는 *.yaml)")

    from simbot.spec import _iter_quantities
    rows: dict[str, dict[str, str]] = {}
    prov: dict[str, set[str]] = {}
    for name, sp in specs.items():
        for path, q in _iter_quantities(sp):
            rows.setdefault(path, {})[name] = _short_value(q.value)
            prov.setdefault(path, set()).add(q.provenance)

    names = list(specs)
    print(f"# params — {len(names)} spec × {len(rows)} 필드\n")
    head = f"| {'필드':<34} | " + " | ".join(f"{n[:22]:<22}" for n in names) + " | prov |"
    print(head)
    print("|" + "-" * 36 + "|" + "|".join(["-" * 24] * len(names)) + "|------|")
    for path in sorted(rows):
        vals = [rows[path].get(n, "—") for n in names]
        differs = len(set(vals)) > 1
        # ⚠ = 아무도 고르지 않은 값 (assumed 이고 전 spec 에서 동일)
        kinds = prov[path]
        mark = "⚠" if (kinds <= {"assumed"} and not differs) else ""
        star = "**" if differs else ""
        print(f"| {star}{path:<32}{star} | "
              + " | ".join(f"{v:<22}" for v in vals) + f" | {mark}{'/'.join(sorted(kinds))} |")
    print("\n⚠ = provenance 가 `assumed` 뿐이고 모든 spec 에서 같은 값 — "
          "**아무도 고르지 않은 기본값**이다.")
    print("**굵게** = spec 간 값이 다른 필드.")
    return 0


def _short_value(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    s = str(v)
    return s if len(s) <= 22 else s[:19] + "…"


# =============================================================================
# calibrate
# =============================================================================
# 정책 상수 `Λ` 가 측정된 커널. 다른 커널로 잰 값과 직접 비교할 수 없다.
#   근거: knowledge/wiki/findings/local-cpu-parallelism.md §2 (측정 헤더)
BASELINE_KERNEL = "3D WCA + Brownian, phi=0.30, Cell nlist(buffer 0.3), dt=1e-4"
CALIBRATE_KERNEL = "2D harmonic trap + Brownian, 쌍 상호작용 없음, nlist 없음"


def cmd_calibrate(args) -> int:
    """이 기계의 처리량 실측.

    ★ **커널이 다르면 상수를 덮어쓰지 않는다.** 정책 상수 `Λ = 6.3e6` 은
      3D WCA + Cell nlist 로 측정됐고, 이 명령이 돌릴 수 있는 유일한 커널은
      쌍 상호작용이 **없는** 조화 트랩이다 (`simbot.run` 에 다른 러너가 없다).
      이웃 목록 작업이 빠진 커널이 더 빠른 것은 당연하고, 그 값으로 전역 상수를
      갱신하면 **WCA 계의 비용 추정이 낙관적으로 틀린다.**
      따라서 커널별 배수로 보고하고, 전역 상수 갱신은 제안하지 않는다.
    """
    pol = load_policy(args.policy)
    claimed = pol.get("hardware.throughput_particle_steps_per_s")
    print(f"# calibrate")
    print(f"  정책 상수 Λ = {claimed:.3g} 입자·스텝/s")
    print(f"    측정 커널: {BASELINE_KERNEL}")
    print(f"  이번 측정 커널: {CALIBRATE_KERNEL}")
    print(f"  ⚠ 처리량 모델은 N ≥ 500 에서 실측됐다. 작은 N 은 오버헤드가 지배한다.\n")

    cases = [(int(n), int(s)) for n, s in (args.cases or [(1000, 20000),
                                                          (4000, 20000),
                                                          (8000, 10000)])]
    print(f"{'N':>7} {'steps':>8} {'wall [s]':>10} {'입자·스텝/s':>14} {'Λ 대비':>8}")
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
    print(f"\n  중앙값 {med:.4g} 입자·스텝/s · Λ 대비 {ratio:.3f}×")
    print(f"\n  ⇒ **트랩 커널 배수 = {ratio:.2f}** (쌍 상호작용 없음 / WCA 기준선)")
    print(f"     트랩 계 비용 추정에 쓸 값이다. 전역 Λ 를 이 값으로 바꾸면 안 된다 —")
    print(f"     WCA 계 추정이 {ratio:.2f}배 낙관적으로 틀린다.")
    print(f"\n  참고: 첫 런의 03_spec.yaml 은 `pair_cost_multiplier: 1.0` 으로 두었다")
    print(f"        (쌍 상호작용이 없으니 기준선 그대로). 이 측정은 그것이 아니라")
    print(f"        {ratio:.2f} 임을 보인다 — 트랩 계는 기준선보다 빠르다.")

    if abs(ratio - 1) > 0.2:
        print(f"\n  ⚠️  기준선과 {abs(ratio - 1):.0%} 어긋난다. 이것이 **커널 차이인지 "
              f"기계 상태 변화인지**\n"
              f"     구별하려면 기준선과 같은 커널(3D WCA φ=0.30)로 재야 하고,\n"
              f"     그 러너는 아직 없다 (`simbot/run.py` 에 트랩 러너만 있다).\n"
              f"     ⇒ 지금 결론할 수 있는 것: **트랩 계 비용은 Λ 로 추정하면 "
              f"{ratio:.2f}배 보수적이다.**")
    else:
        print("\n  ✅ 커널이 달라도 기준선과 20 % 이내 — Λ 를 그대로 써도 된다.")
    print("\n  ★ 코드가 정책 파일을 자동으로 고치지 않는다 — 측정 상수는 사람이 "
          "확인하고 커밋한다.")
    return 0


# =============================================================================
# 진입점
# =============================================================================
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(prog="cli.py", description=__doc__,   # noqa: E501
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--policy", default=None, help="run_policy.yaml 경로")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("run", help="spec 하나 → REPORT.md")
    q.add_argument("spec")
    q.add_argument("--prediction", default=None, help="S2 예측 YAML (봉인 대상)")
    q.add_argument("--runs-root", default=str(REPO / "runs"))
    q.add_argument("--run-id", default=None)
    q.add_argument("--power-threshold", type=float, default=1.0)
    q.add_argument("--force", action="store_true", help="게이트·예산을 무시한다")
    q.set_defaults(func=cmd_run)

    q = sub.add_parser("resume", help="죽은 런 이어받기")
    q.add_argument("rundir")
    q.add_argument("--prediction", default=None)
    q.set_defaults(func=cmd_resume)

    q = sub.add_parser("converge", help="dt·N·초기조건을 흔들어도 답이 같은가")
    q.add_argument("spec")
    q.add_argument("--only", nargs="+", default=None,
                   choices=["dt_half", "dt_double", "N_double", "seed_shift"])
    q.add_argument("--out", default=None)
    q.add_argument("--force", action="store_true")
    q.set_defaults(func=cmd_converge)

    q = sub.add_parser("params", help="여러 런의 파라미터를 가로로")
    q.add_argument("--path", default=str(REPO / "runs"))
    q.set_defaults(func=cmd_params)

    q = sub.add_parser("calibrate", help="이 기계의 처리량 실측")
    q.add_argument("--cases", nargs="+", type=int, default=None,
                   metavar="N STEPS", help="N steps 쌍을 이어서")
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
