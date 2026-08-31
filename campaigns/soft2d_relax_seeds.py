"""soft-r3 relaxation time -- add seeds to decide `tau(A=1)` vs `tau(A=0.1)`.

usage:
  python scripts/soft2d_relax_seeds.py --power        # design power only (does not run)
  python scripts/soft2d_relax_seeds.py               # run at the pre-registered k, then decide
  python scripts/soft2d_relax_seeds.py --analyze-only

## Why a separate script

Editing `soft2d_time_series.py` changes that file's `driver_hash`, and that hash is
**baked into the sealed `03_spec.yaml` of already-committed runs** -- editing it would
make those runs' outputs point at a driver that no longer matches the file on disk. So
it is left alone and this is written fresh.

The geometry, gates and runner are **imported** from that script (a guarantee that the
same numbers are used). Consequently `provenance(driver=[this file, that file])` hashes
**both** -- capturing only one would make "code_hash + driver_hash covers everything"
false.

## What this script is careful about

**(1) The sample size is fixed in advance.** Adding seeds until the result becomes
   significant is optional stopping. `k` is sealed in
   `examples/soft-r3-relax-seeds/prediction.yaml` and this script READS it from there
   -- it is not restated in code.

**(2) The error comes from the seed ensemble.** The `curve_fit` covariance is the fit
   uncertainty of *one mean curve*. `bootstrap_relaxation_over_seeds` resamples the
   seeds with replacement and refits.

**(3) "Indistinguishable" is never written as "the same".** The prediction file owns
   the wording of the verdict.
"""
from __future__ import annotations

import json
import shutil
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from dataclasses import asdict
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from simbot.analysis.structure import (bootstrap_relaxation_over_seeds,  # noqa: E402
                                       fit_relaxation, hex_order_series)
from simbot.estimators import seeds_for_target_sigma                    # noqa: E402
from simbot.io import RunDir, provenance, write_seal                    # noqa: E402
from simbot.policy import load_policy                                   # noqa: E402
from simbot.run import Soft2DRunConfig, run_soft2d                      # noqa: E402

#  ★ Geometry, gates and constants are NOT reimplemented -- a guarantee that the same
#    numbers are used is required
import soft2d_time_series as TS                                         # noqa: E402

SRC = REPO / "examples" / "soft-r3-relax-seeds"
AMPLITUDES = (0.1, 1.0)          # ★ only the two under decision. A=10 already
                                 # separated at 4.7 sigma
DRIVERS = [Path(__file__), Path(TS.__file__)]

#  The stages of the sequential design. **k is owned by the sealed prediction file,
#  not by the code** -- restating it in code lets the two diverge, and then "it was
#  fixed in advance" becomes unprovable.
STAGES = {
    1: {"pred": "prediction.yaml",
        "run_id": "2026-07-29_soft-r3-relax-seeds"},
    2: {"pred": "prediction_stage2.yaml",
        "run_id": "2026-07-29_soft-r3-relax-seeds-stage2"},
}


def stage_from_argv() -> int:
    if "--stage" in sys.argv:
        s = int(sys.argv[sys.argv.index("--stage") + 1])
        if s not in STAGES:
            raise SystemExit(f"⛔ no stage {s}. Available: {sorted(STAGES)}")
        return s
    return 1


def load_prediction(stage: int) -> dict:
    return yaml.safe_load((SRC / STAGES[stage]["pred"]).read_text())


def power_table(pred: dict, stage: int) -> dict:
    """**Recompute** the pre-registered design power and confirm it matches the
    prediction file."""
    dp = pred["design_power"]
    r3 = seeds_for_target_sigma(diff=dp["prior_diff_tau_d"],
                                se_diff=dp["prior_se_diff_tau_d"],
                                k_current=dp["prior_k"], n_sigma=3.0)
    k = int(dp["k_chosen"])
    sigma_at_k = r3["sigma_now"] * np.sqrt(k / dp["prior_k"])

    print(f"## design power -- computed **before running**  (stage {stage})\n")
    print(f"  prior observation (k={dp['prior_k']}):  "
          f"diff = {dp['prior_diff_tau_d']:+.5f} "
          f"± {dp['prior_se_diff_tau_d']:.5f} τ_d = {r3['sigma_now']:.2f}σ")
    print(f"  k needed for 3 sigma = {r3['k_needed_int']}  "
          f"(t={r3['t_quantile']:.3f})")
    print(f"  => pre-registered k = {k}  ->  expected {sigma_at_k:.2f} sigma "
          f"(if the difference is real)")

    #  Stop on disagreement with the prediction file -- the sealed design and the code
    #  must not diverge
    if int(dp["k_for_3sigma"]) != r3["k_needed_int"]:
        raise SystemExit(
            f"⛔ sealed k_for_3sigma={dp['k_for_3sigma']} differs from the "
            f"recomputed {r3['k_needed_int']} -- the design has diverged from the code")
    if stage == 1:
        sigma_at_64 = r3["sigma_now"] * np.sqrt(64 / dp["prior_k"])
        print(f"  ★ k=64 gives only {sigma_at_64:.2f} sigma -- it does not reach 3")
        if abs(float(dp["k_at_64_sigma"]) - sigma_at_64) > 0.01:
            raise SystemExit(f"⛔ k_at_64_sigma mismatch: {dp['k_at_64_sigma']} vs "
                             f"{sigma_at_64:.2f}")
    else:
        pm = pred["stage1_postmortem"]
        print(f"\n  stage-1 post-mortem: expected {pm['predicted_sigma']} sigma -> "
              f"observed "
              f"{pm['observed_sigma']}σ")
        print(f"    (1) the difference shrank by "
              f"{pm['cause_1_biased_diff']['change_pct']:+.1f} % "
              f"(low-seed bias)")
        print(f"    (2) the curve_fit SE was underestimated by "
              f"{pm['cause_2_underestimated_se']['se_ratio_bootstrap_over_fit_A1']}"
              f"x")
        if abs(float(dp["k_chosen_expected_sigma"]) - sigma_at_k) > 0.02:
            raise SystemExit(f"⛔ k_chosen_expected_sigma mismatch")
        print(f"  ★ there is no stage 3 (no_stage_3="
              f"{pred['decision_rule'].get('no_stage_3')})")
    print("  ✅ the sealed design power matches the recomputation")
    return {"k": k, "sigma_expected_at_k": float(sigma_at_k),
            "k_for_3sigma": r3["k_needed_int"], "t_quantile": r3["t_quantile"],
            "stage": stage}


def _one(args) -> dict:
    cfg_dict, outdir = args
    return run_soft2d(Soft2DRunConfig(**cfg_dict), outdir=Path(outdir))


def run_batch(rd: RunDir, rows: list[dict], k: int, policy) -> dict:
    seeds = tuple(range(5, 5 + k))          # same starting point as the prior run
                                            # (includes 5-20)
    outroot = rd.path / "raw_early"
    outroot.mkdir(parents=True, exist_ok=True)
    jobs = []
    for r in rows:
        if r["amplitude"] not in AMPLITUDES:
            continue
        for s in seeds:
            label = f"A{r['amplitude']:g}_s{s}"
            cfg = Soft2DRunConfig(
                amplitude=r["amplitude"], n_particles=TS.N_PARTICLES,
                init=TS.INIT, box_shape=TS.BOX_SHAPE, r_min=TS.R_MIN,
                nlist_buffer=TS.NLIST_BUFFER, min_sep_init=TS.MIN_SEP_INIT,
                dt_star=r["dt_star"], equil_tau=0.0, prod_tau=TS.EARLY_TAU,
                n_frames=TS.EARLY_FRAMES, seed=s, label=label)
            jobs.append((asdict(cfg), str(outroot / label)))

    kk = policy.concurrency("default")
    print(f"\n## S5 -- {len(jobs)} runs "
          f"({len(AMPLITUDES)} amplitudes x {k} seeds, concurrency {kk})")
    t0 = time.perf_counter()
    done, failed = [], []
    with ProcessPoolExecutor(max_workers=kk) as ex:
        futs = {ex.submit(_one, j): j[0]["label"] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            label = futs[fut]
            try:
                out = fut.result()
                done.append({"label": label, "wall_s": out["wall_s"],
                             "fails": out["guards"]["failures"]})
            except Exception as e:                       # noqa: BLE001
                failed.append({"label": label, "error": repr(e)})
                print(f"  ⛔ {label}: {e!r}")
            if i % 100 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} …")
    wall = time.perf_counter() - t0
    n_fail_guard = sum(1 for d in done if d["fails"])
    print(f"  batch wall {wall:.1f} s . failed {len(failed)} . "
          f"guard violations {n_fail_guard}")
    return {"done": len(done), "failed": failed, "batch_wall_s": wall,
            "concurrency": kk, "n_jobs": len(jobs), "seeds": list(seeds),
            "guard_violations": n_fail_guard}


def analyze(rd: RunDir, pred: dict) -> dict:
    rule = pred["decision_rule"]
    root = rd.path / "raw_early"
    out: dict = {}
    for A in AMPLITUDES:
        dirs = sorted(p for p in root.glob(f"A{A:g}_s*")
                      if (p / "samples.npz").exists())
        if not dirs:
            raise SystemExit(f"⛔ no runs for A={A}")
        curves, t_ref = [], None
        for d in dirs:
            z = np.load(d / "samples.npz")
            cfg = json.loads((d / "manifest.json").read_text())["config"]
            stride = int(z["stride"][0])
            n = z["traj"].shape[0]
            t = np.concatenate([[0.0],
                                np.arange(1, n + 1) * stride * cfg["dt_star"]])
            frames = np.concatenate(
                [z["init_pos"][None].astype(np.float32), z["traj"]])
            s = hex_order_series(frames, Lx=float(z["box"][0]),
                                 Ly=float(z["box"][1]), t_star=t,
                                 coord_range=(3, 12))
            curves.append(s.defect_fraction)
            t_ref = t if t_ref is None else t_ref
        mat = np.array(curves)
        tail = mat[:, t_ref >= 0.5 * t_ref[-1]]
        noise = float(mat.mean(axis=0)[t_ref >= 0.5 * t_ref[-1]].std(ddof=1))
        fit = fit_relaxation(t_ref, mat.mean(axis=0), noise=noise)
        boot = bootstrap_relaxation_over_seeds(
            t_ref, mat, n_resample=int(rule["n_resample"]),
            seed=int(rule["bootstrap_seed"]), noise=noise)
        out[f"A{A:g}"] = {
            "amplitude": A, "n_seeds": len(dirs),
            "tau": fit.tau, "tau_se_fit": fit.tau_se,
            "tau_se_bootstrap": boot["tau_se_bootstrap"],
            "tau_ci95": list(boot["tau_ci95"]),
            "se_ratio_bootstrap_over_fit": boot.get(
                "se_ratio_bootstrap_over_fit"),
            "bootstrap_converged": boot["n_converged"],
            "relax_amplitude": fit.amplitude,
            "relax_converged": bool(fit.converged), "relax_note": fit.note,
            "r_squared": fit.r_squared, "defect_frame_sd": noise,
            "defect_at_t0": float(mat[:, 0].mean()),
            "defect_tail_mean": float(tail.mean()),
        }
        print(f"  A={A:<5g} k={len(dirs):<4d} τ = {fit.tau:.5f}  "
              f"SE_fit {fit.tau_se:.5f} · SE_boot {boot['tau_se_bootstrap']:.5f} "
              f"(ratio {boot.get('se_ratio_bootstrap_over_fit', float('nan')):.2f})  "
              f"amplitude/noise {abs(fit.amplitude)/noise:.1f}x")

    a, b = out["A0.1"], out["A1"]
    for tag, key in (("bootstrap", "tau_se_bootstrap"), ("fit", "tau_se_fit")):
        se = float(np.hypot(a[key], b[key]))
        diff = b["tau"] - a["tau"]
        out[f"comparison_{tag}"] = {
            "diff": diff, "se_diff": se,
            "sigma": abs(diff) / se if se else float("nan")}
    return out


def main() -> int:
    policy = load_policy()
    stage = stage_from_argv()
    pred = load_prediction(stage)
    geo = TS.geometry()
    thresholds, rows = TS.gate_table(policy)
    power = power_table(pred, stage)
    if "--power" in sys.argv:
        return 0

    rd = RunDir.create(REPO / "runs", STAGES[stage]["run_id"])
    if "--analyze-only" not in sys.argv:
        shutil.copy2(SRC / STAGES[stage]["pred"], rd.file("prediction"))
        rd.write_json("prediction_json", pred)
        rd.write_json("spec", {
            "source": "scripts/soft2d_relax_seeds.py",
            "stage": stage,
            "provenance": provenance(DRIVERS),
            "parent_run": pred["parent_run"],
            "amplitudes": list(AMPLITUDES), "k_seeds": power["k"],
            "geometry": geo, "thresholds": thresholds,
            "gates": [r for r in rows if r["amplitude"] in AMPLITUDES],
            "early_tau": TS.EARLY_TAU, "early_frames": TS.EARLY_FRAMES,
            "design_power": power})
        seal = write_seal(rd)
        print(f"\n  🔒 sealed {seal.name} -- "
              f"{len(seal.read_text().splitlines())} documents (before running)")
        batch = run_batch(rd, rows, power["k"], policy)
        rd.write_json("manifest", batch)

    print("\n## S7 -- relaxation time (bootstrap errors)")
    res = analyze(rd, pred)
    res["_provenance_at_analysis"] = provenance(DRIVERS)
    res["_design_power"] = power
    rd.write_json("metrics", res)

    rule = pred["decision_rule"]
    thr = float(rule["threshold_sigma"])
    cb, cf = res["comparison_bootstrap"], res["comparison_fit"]
    print(f"\n## verdict -- the pre-registered rule "
          f"(SE source: {rule['se_source']}, threshold {thr:g} sigma)\n")
    print(f"  {'SE source':<12} {'diff':>10} {'SE_diff':>10} {'sigma':>7}")
    print("  " + "-" * 42)
    for tag, c in (("bootstrap ★", cb), ("curve_fit", cf)):
        print(f"  {tag:<12} {c['diff']:>+10.5f} {c['se_diff']:>10.5f} "
              f"{c['sigma']:>7.2f}")
    verdict = (rule["verdict_if_above"] if cb["sigma"] > thr
               else rule["verdict_if_below"])
    print(f"\n  → {verdict}")
    print(f"  (expected: {power['sigma_expected_at_k']:.2f} sigma . "
          f"actual: {cb['sigma']:.2f} sigma)")
    print("\n  proposed_by: agent · confirmed_by: null")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
