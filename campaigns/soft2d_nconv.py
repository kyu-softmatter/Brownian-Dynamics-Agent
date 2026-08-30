"""soft-r3 `N` convergence -- identify the phase from `N=256` plus the finite-size
exponent of `psi6`.

usage:
  python scripts/soft2d_nconv.py --gates        # gates and cost only (does not run)
  python scripts/soft2d_nconv.py                # all of S5 through S7
  python scripts/soft2d_nconv.py --analyze-only

## Two things closed at once

**(1) Card §9** -- `N=100` did not satisfy `A=10`'s `r_cut` requirement (`N >= 252`).
   At `N=256`, `beta*U(r_cut)` drops from `0.0904` to `0.0216 kT`.

**(2) Card §10** -- the `g6(r)` exponent `eta6` was unimplemented. It is obtained as
   `eta6 = 4p` from the `p` in `|<psi6>| ~ N^{-p}`. **Without fitting `g6(r)`** --
   from two system sizes alone.

## What this script is careful about

**(1) It does not re-run the parent run.** The `N=100` values are read from the
   committed `runs/2026-07-29_soft-r3-time-resolved/metrics.json`.

**(2) It keeps the seeds disjoint** (`21-24`). The parent run used `5-8` and the
   relaxation sweep used `5-1540`.

**(3) It counts coordination kinds with an aggregate estimator.** A per-frame threshold
   depends on `N` -- one particle is `1 %` at `N=100` but `0.39 %` at `N=256`.
   Basis: `findings/fraction-threshold-flips-meaning-between-per-frame-and-aggregate.md`
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

from simbot.analysis.structure import (hex_order_series,             # noqa: E402
                                       psi6_finite_size_exponent,
                                       structure_factor, zahn_phase)
from simbot.analysis.trap import aggregate_seeds                     # noqa: E402
from simbot.build import box_si_for_coverage                         # noqa: E402
from simbot.io import RunDir, provenance, write_seal                 # noqa: E402
from simbot.nondim import dt_max_force, dt_max_thermal               # noqa: E402
from simbot.policy import load_policy                                # noqa: E402
from simbot.run import (Soft2DRunConfig, measure_max_force_soft2d,   # noqa: E402
                        run_soft2d)
from simbot.units import scales_soft2d                               # noqa: E402
from simbot.viz import FigureSet, plot_finite_size_scaling           # noqa: E402

import soft2d_time_series as TS                                     # noqa: E402

RUN_ID = "2026-07-29_soft-r3-nconv"
SRC = REPO / "examples" / "soft-r3-nconv"
PARENT = REPO / "runs" / "2026-07-29_soft-r3-time-resolved" / "metrics.json"
DRIVERS = [Path(__file__), Path(TS.__file__)]

N_NEW, N_REF = 256, 100
AMPLITUDES = (0.1, 1.0, 10.0)
SEEDS = (21, 22, 23, 24)
TOTAL_TAU, N_FRAMES = TS.TOTAL_TAU, TS.N_FRAMES
WINDOW = {0.1: (20.0, 30.0), 1.0: (20.0, 30.0), 10.0: (60.0, 80.0)}


def geometry(n: int) -> dict:
    g = box_si_for_coverage(n_particles=n, sigma_si=TS.SIGMA_SI,
                            coverage_max=TS.COVERAGE_MAX,
                            d_over_sigma_round=TS.D_OVER_SIGMA)
    sc = scales_soft2d(d_si=g["d_si"], sigma_si=TS.SIGMA_SI, T_si=TS.T_SI)
    return {**g, "tau_d_si": sc.time_si, "n_particles": n}


def gate_table(policy) -> tuple[dict, list[dict]]:
    ts = policy.timestep
    d_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    d_F = float(ts.get("max_force_displacement_sigma", 0.005))
    dt_th = dt_max_thermal(d_th, sigma=1.0, D0=1.0)
    rows = []
    for A in AMPLITUDES:
        probe = Soft2DRunConfig(amplitude=A, n_particles=N_NEW, init=TS.INIT,
                                box_shape=TS.BOX_SHAPE, r_min=TS.R_MIN,
                                nlist_buffer=TS.NLIST_BUFFER,
                                min_sep_init=TS.MIN_SEP_INIT, seed=SEEDS[0],
                                dt_star=1e-6)
        m = measure_max_force_soft2d(probe)
        dt_f = dt_max_force(d_F, sigma=1.0, gamma=1.0,
                            max_force=m["max_force_star"])
        cands = {"thermal_displacement": dt_th, "force_displacement": dt_f}
        active = {k: v for k, v in cands.items() if v is not None}
        dom = min(active, key=active.get)
        dt = active[dom]
        rows.append({"amplitude": A, "dt_star": dt, "dominant_gate": dom,
                     "dt_max_thermal": dt_th, "dt_max_force": dt_f,
                     "steps": int(round(TOTAL_TAU / dt)),
                     "max_force_star": m["max_force_star"],
                     "r_cut": m["r_cut"], "beta_u_at_rcut": m["beta_u_at_rcut"],
                     "u_rel_to_nearest": m["u_rel_to_nearest"],
                     "Lx": m["Lx"], "Ly": m["Ly"]})
    return {"delta_thermal": d_th, "delta_force": d_F}, rows


def print_gates(geo: dict, rows: list[dict], policy, pred: dict) -> None:
    lam = float(policy.get("hardware.throughput_particle_steps_per_s", 6.3e6))
    k = policy.concurrency("default")
    eff = policy.efficiency(k)
    ref = geometry(N_REF)
    print(f"## geometry  N: {N_REF} -> {N_NEW}\n"
          f"  L*: {ref['L_star']:.0f} → {geo['L_star']:.0f}  ·  "
          f"L: {ref['L_si']*1e6:.0f} → {geo['L_si']*1e6:.0f} µm  ·  "
          f"coverage {geo['coverage']:.4%} (**independent of N -- n* = 1**)\n"
          f"  tau_d = {geo['tau_d_si']:.1f} s = {geo['tau_d_si']/60:.2f} min "
          f"(independent of N: d is set by sigma alone)\n")
    hdr = (f"{'A':>6} {'Γ':>7} {'r_cut':>7} {'βU(rc)':>8} {'N=100 βU':>9} "
           f"{'max|F*|':>9} {'governing':<22} {'dt*':>10} {'steps':>10} "
           f"{'est. wall':>9}")
    print(hdr); print("-" * len(hdr))
    sealed = {it["quantity"]: it["value"] for it in pred["items"]}
    total = 0.0
    for r in rows:
        A = r["amplitude"]
        wall = N_NEW * r["steps"] / (lam * eff) * 3.4      # frame-overhead factor
        total = max(total, wall)
        key = f"beta_u_at_rcut__A{A:g}__N256"
        ref_bu = {0.1: 0.0009, 1.0: 0.0090, 10.0: 0.0904}[A]
        print(f"{A:>6g} {zahn_phase(A)['gamma']:>7.2f} {r['r_cut']:>7.3f} "
              f"{r['beta_u_at_rcut']:>8.4f} {ref_bu:>9.4f} "
              f"{r['max_force_star']:>9.2f} {r['dominant_gate']:<22} "
              f"{r['dt_star']:>10.3g} {r['steps']:>10d} {wall:>8.0f}s")
        #  ★ stop if it disagrees with the sealed prediction
        if key in sealed and abs(r["beta_u_at_rcut"] - sealed[key]) > 5e-4:
            raise SystemExit(f"⛔ {key}: sealed {sealed[key]} vs measured "
                             f"{r['beta_u_at_rcut']:.4f} -- the design has diverged")
    print(f"\n  ✅ beta*U(r_cut) agrees with the sealed prediction")
    budget = policy.wall_budget_s
    print(f"  {len(rows)*len(SEEDS)} runs . concurrency {k} "
          f"(efficiency {eff:.2f}) . "
          f"longest run estimated {total:.0f} s "
          f"{'<=' if total <= budget else '>'} budget {budget:.0f} s")
    if total > budget:
        raise SystemExit("⛔ over budget -- reporting without running")


def _one(args) -> dict:
    cfg_dict, outdir = args
    return run_soft2d(Soft2DRunConfig(**cfg_dict), outdir=Path(outdir))


def run_batch(rd: RunDir, rows: list[dict], policy) -> dict:
    jobs = []
    for r in rows:
        for s in SEEDS:
            label = f"A{r['amplitude']:g}_N{N_NEW}_s{s}"
            cfg = Soft2DRunConfig(
                amplitude=r["amplitude"], n_particles=N_NEW, init=TS.INIT,
                box_shape=TS.BOX_SHAPE, r_min=TS.R_MIN,
                nlist_buffer=TS.NLIST_BUFFER, min_sep_init=TS.MIN_SEP_INIT,
                dt_star=r["dt_star"], equil_tau=0.0, prod_tau=TOTAL_TAU,
                n_frames=N_FRAMES, seed=s, label=label)
            jobs.append((asdict(cfg), str(rd.raw / label)))
    k = policy.concurrency("default")
    print(f"\n## S5 -- {len(jobs)} runs (concurrency {k})")
    t0 = time.perf_counter()
    done, failed = [], []
    with ProcessPoolExecutor(max_workers=k) as ex:
        futs = {ex.submit(_one, j): j[0]["label"] for j in jobs}
        for fut in as_completed(futs):
            label = futs[fut]
            try:
                out = fut.result()
                fails = out["guards"]["failures"]
                done.append({"label": label, "wall_s": out["wall_s"],
                             "fails": fails})
                print(f"  {'✅' if not fails else '⚠'} {label:<20} "
                      f"{out['wall_s']:>7.1f} s  min_sep "
                      f"{out['guards']['min_separation']:.4f} d"
                      + (f"  ⚠ {fails}" if fails else ""))
            except Exception as e:                          # noqa: BLE001
                failed.append({"label": label, "error": repr(e)})
                print(f"  ⛔ {label}: {e!r}")
    wall = time.perf_counter() - t0
    print(f"\n  batch wall {wall:.1f} s . failed {len(failed)}")
    return {"done": done, "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "seeds": list(SEEDS),
            "n_particles": N_NEW}


def analyze(rd: RunDir) -> dict:
    parent = json.loads(PARENT.read_text())
    out: dict = {}
    for A in AMPLITUDES:
        dirs = sorted(p for p in rd.raw.glob(f"A{A:g}_N{N_NEW}_s*")
                      if (p / "samples.npz").exists())
        if not dirs:
            raise SystemExit(f"⛔ no runs for A={A}")
        lo, hi = WINDOW[A]
        per_seed, agg_frac, n_tot, min_seps, late_frames = [], None, 0, [], []
        for d in dirs:
            z = np.load(d / "samples.npz")
            man = json.loads((d / "manifest.json").read_text())
            cfg = man["config"]
            stride = int(z["stride"][0])
            t = np.arange(1, z["traj"].shape[0] + 1) * stride * cfg["dt_star"]
            s = hex_order_series(z["traj"], Lx=float(z["box"][0]),
                                 Ly=float(z["box"][1]), t_star=t,
                                 coord_range=(3, 12))
            m = (t >= lo) & (t < hi)
            per_seed.append({"psi6_global": float(s.psi6_global[m].mean()),
                             "psi6_local": float(s.psi6_local[m].mean()),
                             "defect_fraction": float(s.defect_fraction[m].mean())})
            #  aggregate histogram (a per-frame threshold is N-dependent, so it is
            #  not used)
            f = s.coord_fraction[m].sum(axis=0)
            agg_frac = f if agg_frac is None else agg_frac + f
            n_tot += int(m.sum())
            min_seps.append(man["guards"]["min_separation"])
            late_frames.append(z["traj"][m])
        agg_frac = agg_frac / max(n_tot, 1)
        labels = np.arange(3, 13)
        aggd = {m_: aggregate_seeds([r[m_] for r in per_seed])
                for m_ in ("psi6_global", "psi6_local", "defect_fraction")}
        sk = structure_factor(np.concatenate(late_frames, axis=0),
                              Lx=float(z["box"][0]), Ly=float(z["box"][1]),
                              n_max=18)
        ref = parent[f"A{A:g}"]
        fit = psi6_finite_size_exponent(
            np.array([N_REF, N_NEW], dtype=float),
            np.array([ref["psi6_global"]["mean"], aggd["psi6_global"].mean]),
            np.array([ref["psi6_global"]["se"], aggd["psi6_global"].se]))
        out[f"A{A:g}"] = {
            "amplitude": A, "n_particles": N_NEW, "n_seeds": len(dirs),
            "window": [lo, hi], "zahn": zahn_phase(A),
            **{m_: {"mean": a.mean, "se": a.se, "values": a.values}
               for m_, a in aggd.items()},
            "coord_kinds_aggregate": int((agg_frac > 0.005).sum()),
            "coordination_hist_aggregate": {int(z_): float(v) for z_, v
                                            in zip(labels, agg_frac) if v > 0},
            "five_seven_balance": float(
                abs(agg_frac[labels == 5][0] - agg_frac[labels == 7][0])
                / max(agg_frac[labels == 5][0] + agg_frac[labels == 7][0], 1e-12)),
            "sixfold_modulation": sk.sixfold_modulation,
            "min_separation_d": min(min_seps),
            "reference_N100": {k_: ref[k_] for k_ in
                               ("psi6_global", "psi6_local", "defect_fraction")},
            "finite_size": {"p": fit.p, "p_se": fit.p_se, "eta6": fit.eta6,
                            "eta6_se": fit.eta6_se, "reading": fit.reading,
                            "n_points": fit.n_points},
        }
        print(f"  A={A:<5g} ψ₆ {aggd['psi6_global'].mean:.4f}±"
              f"{aggd['psi6_global'].se:.4f} (N=100: "
              f"{ref['psi6_global']['mean']:.4f}) → p = {fit.p:.3f}±{fit.p_se:.3f}"
              f"  η₆ = {fit.eta6:.2f}  **{fit.reading}**")
    return out


def figures(rd: RunDir, metrics: dict) -> FigureSet:
    fs = FigureSet(rd.figs)
    from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,
                                           LIQUID_EXPONENT_P,
                                           psi6_finite_size_exponent)
    data, local = {}, {}
    for A in AMPLITUDES:
        m = metrics[f"A{A:g}"]
        ref = m["reference_N100"]
        N = [N_REF, N_NEW]
        y = [ref["psi6_global"]["mean"], m["psi6_global"]["mean"]]
        se = [ref["psi6_global"]["se"], m["psi6_global"]["se"]]
        data[A] = {"N": N, "psi6": y, "se": se,
                   "fit": psi6_finite_size_exponent(np.array(N, dtype=float),
                                                    np.array(y), np.array(se))}
        local[A] = {"N": N,
                    "defect": [ref["defect_fraction"]["mean"],
                               m["defect_fraction"]["mean"]],
                    "defect_se": [ref["defect_fraction"]["se"],
                                  m["defect_fraction"]["se"]]}
    plot_finite_size_scaling(fs, data, local_data=local,
                             exponent_liquid=LIQUID_EXPONENT_P,
                             exponent_hexatic=KTHNY_ETA6_HEXATIC_LIQUID / 4.0)
    fs.skip("voronoi", "the parent run already characterised the defects -- this run "
                       "looks only at N dependence")
    fs.skip("early_transient", "the relaxation sweep already closed the transient with "
                               "1513 seeds")
    return fs


def check(pred: dict, metrics: dict, rows: list[dict]) -> list[dict]:
    gates = {r["amplitude"]: r for r in rows}
    got: dict = {}
    for A in AMPLITUDES:
        m = metrics[f"A{A:g}"]
        got[f"psi6_global__A{A:g}__N256"] = m["psi6_global"]["mean"]
        got[f"psi6_exponent_p__A{A:g}"] = m["finite_size"]["p"]
        got[f"defect_fraction__A{A:g}__N256"] = m["defect_fraction"]["mean"]
        got[f"psi6_local__A{A:g}__N256"] = m["psi6_local"]["mean"]
        got[f"beta_u_at_rcut__A{A:g}__N256"] = gates[A]["beta_u_at_rcut"]
    got["coord_kinds_aggregate__A10__N256"] = metrics["A10"]["coord_kinds_aggregate"]

    rows_out = []
    for it in pred["items"]:
        q = it["quantity"]
        v = got.get(q)
        tol = str(it["tolerance"]).strip()
        if v is None:
            verdict = "NOT_EVALUATED"
        elif tol.startswith(">"):
            verdict = "PASS" if v > float(tol[1:]) else "FAIL"
        elif tol.startswith("<"):
            verdict = "PASS" if v < float(tol[1:]) else "FAIL"
        else:
            verdict = ("PASS" if abs(v - float(it["value"]))
                       <= float(tol.lstrip("±")) else "FAIL")
        cv = it.get("competing_value")
        note = ""
        if isinstance(cv, (int, float)) and v is not None:
            d_pred = abs(v - float(it["value"]))
            d_comp = abs(v - float(cv))
            note = (f"competing hypothesis `{cv:.5g}` at distance `{d_comp:.4g}` "
                    f"vs "
                    f"`{d_pred:.4g}` → "
                    + ("**favours the prediction**" if d_pred < d_comp
                       else "**favours the competitor**"))
        rows_out.append({"quantity": q, "predicted": it["value"],
                         "tolerance": tol, "measured": v, "verdict": verdict,
                         "competing_value": cv, "note": note,
                         "discriminates": it.get("discriminates", "")})
    return rows_out


def main() -> int:
    policy = load_policy()
    pred = yaml.safe_load((SRC / "prediction.yaml").read_text())
    geo = geometry(N_NEW)
    thresholds, rows = gate_table(policy)
    print_gates(geo, rows, policy, pred)
    if "--gates" in sys.argv:
        return 0

    rd = RunDir.create(REPO / "runs", RUN_ID)
    if "--analyze-only" not in sys.argv:
        shutil.copy2(SRC / "prediction.yaml", rd.file("prediction"))
        rd.write_json("prediction_json", pred)
        rd.write_json("spec", {"source": "scripts/soft2d_nconv.py",
                               "provenance": provenance(DRIVERS),
                               "parent_run": pred["parent_run"],
                               "n_particles": N_NEW, "n_ref": N_REF,
                               "geometry": geo, "thresholds": thresholds,
                               "gates": rows, "seeds": list(SEEDS),
                               "total_tau": TOTAL_TAU, "n_frames": N_FRAMES,
                               "window": {str(k): v for k, v in WINDOW.items()}})
        seal = write_seal(rd)
        print(f"\n  🔒 sealed {seal.name} -- "
              f"{len(seal.read_text().splitlines())} documents (before running)")
        rd.write_json("manifest", run_batch(rd, rows, policy))

    print("\n## S7 -- finite-size analysis")
    metrics = analyze(rd)
    metrics["_provenance_at_analysis"] = provenance(DRIVERS)
    fs = figures(rd, metrics)
    rd.write("figures", fs.figures_md())
    checks = check(pred, metrics, rows)
    metrics["_checks"] = checks
    rd.write_json("metrics", metrics)

    print("\n## prediction comparison (the verdict is a proposal -- "
          "confirmed_by: null)")
    hdr = f"{'item':<38} {'predicted':>10} {'competing':>10} {'measured':>10} verdict"
    print(hdr); print("-" * (len(hdr) + 6))
    for c in checks:
        cv = c["competing_value"]
        cvs = f"{cv:.5g}" if isinstance(cv, (int, float)) else "—"
        print(f"{c['quantity']:<38} {float(c['predicted']):>10.5g} {cvs:>10} "
              f"{c['measured']:>10.5g} {c['verdict']}")
    n_p = sum(1 for c in checks if c["verdict"] == "PASS")
    print(f"\n  PASS {n_p} · FAIL {len(checks) - n_p}")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
