"""soft-r3 `A` sweep -- **time-resolved**. When do `g(r)`, the Voronoi structure
and the defect structure get made?

usage:
  python scripts/soft2d_time_series.py --gates        # gates and cost only (no run)
  python scripts/soft2d_time_series.py                # all of S5→S7
  python scripts/soft2d_time_series.py --analyze-only # re-analyse an existing raw/

## What differs from the prior run

`runs/2026-07-29_soft-r3-2d-A-sweep` **threw the equilibration window away** and
analysed only the second half of production. That answers "what did it become"
and **does not answer "when did it become that."** This run sets
`equil_tau = 0` and samples the whole thing from `t = 0`.

## What this script is careful about

**① It does not re-write the gate thresholds.** It calls
   `simbot.nondim.dt_max_*` in reduced units (`d = 1`, `D₀ = 1`, `γ = 1`), and
   reads the thresholds from `config/run_policy.yaml` -- nailed into the script
   they would not follow a policy change (precedent: 2026-07-28).

**② It actually computes `max|F|`.** Estimating it is forbidden (master_plan
   §5.4). It is measured on the random initial placement, and the maximum during
   production is recorded too -- because on a symmetric placement `max|F| = 0`,
   which makes the gate toothless (card §7).

**③ It seals **before** running.** A seal made after the run guarantees nothing.

**④ The seeds were split off from the prior run** (`1–4` → `5–8`). HOOMD
   `Brownian` reproduces bit-for-bit on the same seed, so reusing the seeds would
   turn the late-window comparison into an **arithmetic identity**.
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

from simbot.analysis.structure import (fit_relaxation, hex_order_series,  # noqa: E402
                                       rdf_windows, structure_factor,
                                       voronoi_frame, zahn_phase)
from simbot.analysis.trap import aggregate_seeds                       # noqa: E402
from simbot.build import box_si_for_coverage                           # noqa: E402
from simbot.io import RunDir, provenance, write_seal                   # noqa: E402
from simbot.nondim import dt_max_force, dt_max_thermal                 # noqa: E402
from simbot.policy import load_policy                                  # noqa: E402
from simbot.run import Soft2DRunConfig, measure_max_force_soft2d, run_soft2d  # noqa: E402
from simbot.units import scales_soft2d, water_viscosity_si             # noqa: E402
from simbot.viz import (FigureSet, plot_early_transient,               # noqa: E402
                        plot_rdf_evolution, plot_structure_timeseries,
                        plot_voronoi_timelapse)

# --- user directives (S1 delta D1–D4) ---------------------------------------
RUN_ID = "2026-07-29_soft-r3-time-resolved"
SRC = REPO / "examples" / "soft-r3-time-resolved"
AMPLITUDES = (0.1, 1.0, 10.0)          # D1 — A=100 excluded
BOX_SHAPE = "square"                   # D2 — L_x = L_y fixed
INIT = "random"                        # D2 — in a square box hex is an artefact
SIGMA_SI = 5.0e-6                      # D3
D_OVER_SIGMA = 3.0                     # D3 — coverage 8.73 % < 10 %
COVERAGE_MAX = 0.10
T_SI = 298.15
N_PARTICLES = 100
SEEDS = (5, 6, 7, 8)                   # ★ split off from the prior run (1–4)
TOTAL_TAU = 80.0                       # D4 — sample all of it (prior run's total)
N_FRAMES = 400                         # stride = 0.2 tau_d
R_MIN = 0.05                           # Table lower bound
MIN_SEP_INIT = 0.8
NLIST_BUFFER = 0.1
N_RDF_WINDOWS = 4
VORONOI_TIMES = (0.0, 1.0, 10.0, 40.0, 80.0)   # [tau_d] timelapse times (approx)

#  ★ A pass dedicated to the early transient (added 2026-07-29)
#  The main pass's stride is 80/400 = 0.2 τ_d. But the time for the initial
#  placement's excluded-volume shell (`min_sep = 0.8 d`) to fill in is, by free
#  diffusion, `t* ≈ 0.3²/4 = 0.023 τ_d` -- **9x faster than the first frame.**
#  That is, the main pass swallows the whole transient inside its first frame and
#  makes it look like "steady state from the start".
#  ⇒ Run a short, dense pass separately on the same seeds. It costs 1/40 of the
#  main pass.
EARLY_TAU = 2.0
EARLY_FRAMES = 400                             # stride = 0.005 tau_d
#  ★ 16 seeds. 12 runs took 3.9 s, so 48 runs is ~16 s -- **the regime where
#    error bars are free** (CLAUDE.md §resource policy). The A dependence of the
#    relaxation time is this run's only new claim, so that is where the seeds go.
#    The main pass cannot make the same choice: A=10 costs 65 s/run.
EARLY_SEEDS = tuple(range(5, 21))


# =============================================================================
# S4 — dt gates (reduced units: d = 1, D0 = 1, gamma = 1, kT = 1)
# =============================================================================
def gate_table(policy) -> tuple[dict, list[dict]]:
    """`dt*` and the dominant gate per `A`. **`max|F|` is actually computed.**"""
    ts = policy.timestep
    delta_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    delta_F = float(ts.get("max_force_displacement_sigma", 0.005))
    hard_floor = float(ts.get("hard_floor", 1e-7))

    dt_th = dt_max_thermal(delta_th, sigma=1.0, D0=1.0)      # reduced: d=1, D0=1
    rows = []
    for A in AMPLITUDES:
        probe = Soft2DRunConfig(
            amplitude=A, n_particles=N_PARTICLES, init=INIT, box_shape=BOX_SHAPE,
            r_min=R_MIN, nlist_buffer=NLIST_BUFFER, min_sep_init=MIN_SEP_INIT,
            seed=SEEDS[0], dt_star=1e-6)
        m = measure_max_force_soft2d(probe)
        dt_F = dt_max_force(delta_F, sigma=1.0, gamma=1.0,
                            max_force=m["max_force_star"])
        cands = {"thermal_displacement": dt_th, "force_displacement": dt_F}
        active = {k: v for k, v in cands.items() if v is not None}
        dominant = min(active, key=active.get)
        dt = active[dominant]
        if dt < hard_floor:
            raise SystemExit(f"A={A}: Δt* = {dt:.3g} < hard_floor {hard_floor:g}")
        steps = int(round(TOTAL_TAU / dt))
        rows.append({
            "amplitude": A, "gamma_zahn": zahn_phase(A)["gamma"],
            "phase_zahn": zahn_phase(A)["phase_zahn"],
            "max_force_star": m["max_force_star"],
            "dt_max_thermal": dt_th, "dt_max_force": dt_F,
            "dominant_gate": dominant, "dt_star": dt, "steps": steps,
            "r_cut": m["r_cut"], "beta_u_at_rcut": m["beta_u_at_rcut"],
            "u_rel_to_nearest": m["u_rel_to_nearest"],
            "Lx": m["Lx"], "Ly": m["Ly"],
        })
    return {"delta_thermal": delta_th, "delta_force": delta_F,
            "hard_floor": hard_floor}, rows


def geometry() -> dict:
    g = box_si_for_coverage(n_particles=N_PARTICLES, sigma_si=SIGMA_SI,
                            coverage_max=COVERAGE_MAX,
                            d_over_sigma_round=D_OVER_SIGMA)
    sc = scales_soft2d(d_si=g["d_si"], sigma_si=SIGMA_SI, T_si=T_SI)
    eta_si, extrapolated = water_viscosity_si(T_SI)
    return {**g, "tau_d_si": sc.time_si, "D0_si": sc.diffusivity_si,
            "eta_si": eta_si, "eta_extrapolated": extrapolated, "T_si": T_SI,
            "sigma_si": SIGMA_SI, "kT_si": sc.energy_si}


def print_gates(geo: dict, thresholds: dict, rows: list[dict], policy) -> float:
    lam = float(policy.get("hardware.throughput_particle_steps_per_s", 6.3e6))
    k = policy.concurrency("default")
    eff = policy.efficiency(k)
    print(f"## geometry (user directive D3)\n"
          f"  σ = {geo['sigma_si']*1e6:.2f} µm · d/σ = {geo['d_over_sigma']:.3f} → "
          f"d = {geo['d_si']*1e6:.3f} µm · L = {geo['L_si']*1e6:.2f} µm "
          f"(L* = {geo['L_star']:.1f})\n"
          f"  **coverage = {geo['coverage']:.4%}** (cap {COVERAGE_MAX:.1%}) · "
          f"σ/d = {geo['sigma_over_d']:.6g}\n"
          f"  η = {geo['eta_si']*1e3:.4f} mPa·s · D₀ = {geo['D0_si']*1e12:.5f} µm²/s · "
          f"**τ_d = {geo['tau_d_si']:.1f} s = {geo['tau_d_si']/60:.2f} min**\n"
          f"  → {TOTAL_TAU:g} τ_d = {TOTAL_TAU*geo['tau_d_si']/3600:.1f} hours "
          f"of real experiment time\n")
    print(f"## dt gates (displacement criterion, reduced units)  thresholds: "
          f"thermal {thresholds['delta_thermal']} d · "
          f"force {thresholds['delta_force']} d\n")
    hdr = (f"{'A':>6} {'Γ':>7} {'Zahn':<16} {'max|F*|':>9} {'dt_th':>10} "
           f"{'dt_F':>10} {'dominant':<22} {'dt*':>10} {'steps':>10} {'βU(rcut)':>9}")
    print(hdr); print("-" * len(hdr))
    total = 0.0
    for r in rows:
        wall = N_PARTICLES * r["steps"] / (lam * eff)
        total += wall
        print(f"{r['amplitude']:>6g} {r['gamma_zahn']:>7.2f} {r['phase_zahn']:<16} "
              f"{r['max_force_star']:>9.3f} {r['dt_max_thermal']:>10.3g} "
              f"{r['dt_max_force']:>10.3g} {r['dominant_gate']:<22} "
              f"{r['dt_star']:>10.3g} {r['steps']:>10d} {r['beta_u_at_rcut']:>9.4f}")
    n_runs = len(rows) * len(SEEDS)
    per_seed = total * len(SEEDS)
    print(f"\n## cost  {n_runs} runs = {len(rows)} A × {len(SEEDS)} seeds · "
          f"concurrency {k} (efficiency {eff:.2f})")
    print(f"  Λ = {lam:.3g} particle·steps/s → integration estimate "
          f"{per_seed:.0f} s (= {per_seed/60:.1f} min, before concurrency)")
    print(f"  ⚠ the estimate is integration only. Snapshot-extraction overhead "
          f"for {N_FRAMES} frames was large in the prior run "
          f"(A=10: integration 20 s vs measured 67 s)")
    budget = policy.wall_budget_s
    worst = max(N_PARTICLES * r["steps"] / lam for r in rows) * 3.4   # overhead factor
    if worst > budget:
        print(f"\n  ⛔ longest run estimate {worst:.0f} s > budget {budget:.0f} s/run "
              f"-- reporting without running (CLAUDE.md §resource policy)")
    else:
        print(f"  ✅ longest run estimate {worst:.0f} s ≤ budget {budget:.0f} s/run")
    return total


# =============================================================================
# S5 — run
# =============================================================================
def _one(args) -> dict:
    cfg_dict, outdir = args
    cfg = Soft2DRunConfig(**cfg_dict)
    return run_soft2d(cfg, outdir=Path(outdir),
                      extra_manifest={"label": cfg.label})


def run_batch(rd: RunDir, rows: list[dict], policy, *, early: bool = False) -> dict:
    prod_tau = EARLY_TAU if early else TOTAL_TAU
    n_frames = EARLY_FRAMES if early else N_FRAMES
    outroot = (rd.path / "raw_early") if early else rd.raw
    outroot.mkdir(parents=True, exist_ok=True)

    jobs = []
    for r in rows:
        for s in (EARLY_SEEDS if early else SEEDS):
            label = f"A{r['amplitude']:g}_s{s}"
            cfg = Soft2DRunConfig(
                amplitude=r["amplitude"], n_particles=N_PARTICLES, init=INIT,
                box_shape=BOX_SHAPE, r_min=R_MIN, nlist_buffer=NLIST_BUFFER,
                min_sep_init=MIN_SEP_INIT, dt_star=r["dt_star"],
                equil_tau=0.0,                       # ★ D4 — sample from t=0
                prod_tau=prod_tau, n_frames=n_frames, seed=s, label=label)
            jobs.append((asdict(cfg), str(outroot / label)))

    k = policy.concurrency("default")
    tag = "early-transient" if early else "main"
    print(f"\n## S5 — starting the {tag} pass, {len(jobs)} runs "
          f"(concurrency {k}, prod {prod_tau:g} τ_d / {n_frames} frames → "
          f"stride {prod_tau/n_frames:.4g} τ_d)")
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
                             "min_sep": out["guards"]["min_separation"],
                             "fails": fails})
                mark = "✅" if not fails else "⚠"
                print(f"  {mark} {label:<12} {out['wall_s']:>7.1f} s  "
                      f"min_sep {out['guards']['min_separation']:.4f} d"
                      + (f"  ⚠ {fails}" if fails else ""))
            except Exception as e:                   # noqa: BLE001
                failed.append({"label": label, "error": repr(e)})
                print(f"  ⛔ {label:<12} {e!r}")
    wall = time.perf_counter() - t0
    print(f"\n  batch wall {wall:.1f} s")
    return {"done": done, "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "early": early,
            "prod_tau": prod_tau, "n_frames": n_frames}


# =============================================================================
# S7b — the early transient: is there a relaxation **at all**, and how long
# =============================================================================
def analyze_early(rd: RunDir, geo: dict) -> dict:
    """The short dense pass. The initial placement (`t = 0`) is prepended to
    the front of the time series."""
    root = rd.path / "raw_early"
    dirs = sorted(p for p in root.glob("A*_s*") if (p / "samples.npz").exists())
    if not dirs:
        return {}

    series_by_A: dict[float, list] = {}
    fits: dict[float, object] = {}
    summary: dict[str, dict] = {}
    for d in dirs:
        r = load_run(d)
        A = float(r["cfg"]["amplitude"])
        #  ★ Prepend the initial placement as the t=0 frame -- the runner only
        #    starts sampling after the first stride
        frames = np.concatenate([r["init_pos"][None].astype(np.float32), r["traj"]])
        t = np.concatenate([[0.0], r["t_star"]])
        ser = hex_order_series(frames, Lx=r["Lx"], Ly=r["Ly"], t_star=t,
                               coord_range=(3, 12))
        series_by_A.setdefault(A, []).append(ser)

    for A, sers in sorted(series_by_A.items()):
        t = sers[0].t_star
        mat = np.array([s.defect_fraction for s in sers])
        mean = mat.mean(axis=0)
        #  noise = frame-to-frame SD over the second half (of the seed-mean curve)
        tail = mean[t >= 0.5 * t[-1]]
        noise = float(tail.std(ddof=1)) if tail.size > 2 else None
        fit = fit_relaxation(t, mean, noise=noise)
        fits[A] = fit
        psi_mat = np.array([s.psi6_global for s in sers])
        summary[f"A{A:g}"] = {
            "amplitude": A,
            "defect_at_t0": float(mean[0]),
            "defect_tail_mean": float(tail.mean()),
            "defect_frame_sd": noise,
            "psi6_at_t0": float(psi_mat[:, 0].mean()),
            "psi6_tail_mean": float(psi_mat[:, t >= 0.5 * t[-1]].mean()),
            "tau_relax_tau_d": fit.tau, "tau_relax_se": fit.tau_se,
            "tau_relax_s": fit.tau * geo["tau_d_si"],
            "relax_amplitude": fit.amplitude,
            "relax_r_squared": fit.r_squared,
            "relax_converged": bool(fit.converged),
            "relax_note": fit.note,
            "stride_tau_d": float(t[2] - t[1]),
            "n_frames": int(t.size),
        }
        print(f"  A={A:<5g} defects t=0 {mean[0]:.4f} → late {tail.mean():.4f} "
              f"(frame SD {noise:.4f}) · τ = {fit.tau:.4g} τ_d "
              f"({'accepted' if fit.converged else 'rejected'}) "
              f"{'· ' + fit.note if fit.note else ''}")
    return {"series_by_A": series_by_A, "fits": fits, "summary": summary}


# =============================================================================
# S7 — time-resolved analysis
# =============================================================================
def load_run(d: Path) -> dict:
    man = json.loads((d / "manifest.json").read_text())
    z = np.load(d / "samples.npz")
    cfg = man["config"]
    stride = int(z["stride"][0])
    n = z["traj"].shape[0]
    #  frame k is (k+1)*stride steps in -- the runner calls run(stride) first,
    #  then samples
    t_star = (np.arange(1, n + 1) * stride * cfg["dt_star"])
    return {"man": man, "cfg": cfg, "traj": z["traj"], "energy": z["energy"],
            "max_force": z["max_force"], "init_pos": z["init_pos"],
            "Lx": float(z["box"][0]), "Ly": float(z["box"][1]),
            "t_star": t_star, "stride": stride}


def analyze(rd: RunDir, geo: dict) -> dict:
    dirs = sorted(p for p in rd.raw.glob("A*_s*") if (p / "samples.npz").exists())
    if not dirs:
        raise SystemExit(f"⛔ no runs in {rd.raw}")

    by_A: dict[float, list[dict]] = {}
    for d in dirs:
        r = load_run(d)
        A = float(r["cfg"]["amplitude"])
        ser = hex_order_series(r["traj"], Lx=r["Lx"], Ly=r["Ly"],
                               t_star=r["t_star"], coord_range=(3, 12))
        by_A.setdefault(A, []).append({"run": r, "series": ser, "dir": d})
        print(f"  analysing {d.name:<12} frames {ser.n_frames} · "
              f"ψ₆(last) {ser.psi6_global[-1]:.4f}")

    #  --- late window: the same times as the prior run's second half ---
    #      the prior total length differed per A (A ≤ 1 ran 30 τ_d, A = 10 ran 80)
    prior_window = {0.1: (20.0, 30.0), 1.0: (20.0, 30.0), 10.0: (60.0, 80.0)}
    metrics: dict = {}
    windows_by_A, series_by_A, energy_by_A = {}, {}, {}

    for A, runs in sorted(by_A.items()):
        sers = [x["series"] for x in runs]
        series_by_A[A] = sers
        t = sers[0].t_star
        lo, hi = prior_window[A]

        wins = [s.window_mean(lo, hi) for s in sers]
        agg = {m: aggregate_seeds([w[m] for w in wins])
               for m in ("psi6_global", "psi6_local", "defect_fraction")}
        #  energy: per-seed mean over the same window
        e_pp = []
        for x in runs:
            r = x["run"]
            m = (r["t_star"] >= lo) & (r["t_star"] < hi)
            e_pp.append(float(r["energy"][m].mean()) / N_PARTICLES)
        agg["energy_pp"] = aggregate_seeds(e_pp)
        energy_by_A[A] = (t, [x["run"]["energy"] / N_PARTICLES for x in runs])

        #  coordination kinds and 5-7 imbalance in the final window (seed mean)
        kinds = [float(s.coord_kinds[(s.t_star >= lo) & (s.t_star < hi)].mean())
                 for s in sers]
        #  ★ Coordination-kind count the **same way** as the prior run: count the
        #    kinds whose fraction is > 0.5 % in a histogram aggregated over all
        #    frames and seeds. Counting per frame gives a different number -- at
        #    `N = 100` a single particle is already 1 %, so a per-frame threshold
        #    degenerates into "passes if it exists at all". Both estimators are
        #    reported.
        agg_frac = np.zeros(sers[0].coord_labels.size)
        n_tot = 0
        for s in sers:
            m_w = (s.t_star >= lo) & (s.t_star < hi)
            agg_frac += s.coord_fraction[m_w].sum(axis=0)
            n_tot += int(m_w.sum())
        agg_frac /= max(n_tot, 1)
        kinds_aggregate = int((agg_frac > 0.005).sum())
        bal = [float(s.five_seven_balance[(s.t_star >= lo) & (s.t_star < hi)].mean())
               for s in sers]

        #  t90 — first time the seed-mean curve reaches 90 % of the late-window mean
        psi_mat = np.array([s.psi6_global for s in sers])
        psi_mean = psi_mat.mean(axis=0)
        target = 0.9 * agg["psi6_global"].mean
        reached = np.nonzero(psi_mean >= target)[0]
        t90 = float(t[reached[0]]) if reached.size else float("nan")
        #  ★ A liquid is already above the target at t=0 (a random placement IS a
        #    liquid placement) → t90 ≈ the first frame. That fact is a result

        #  time-windowed g(r) — seeds pooled per window (frames = 4 seeds × 100)
        all_traj = np.concatenate([x["run"]["traj"] for x in runs], axis=0)
        all_t = np.concatenate([x["run"]["t_star"] for x in runs])
        order = np.argsort(all_t, kind="stable")
        windows_by_A[A] = rdf_windows(all_traj[order], Lx=runs[0]["run"]["Lx"],
                                      Ly=runs[0]["run"]["Ly"],
                                      t_star=all_t[order],
                                      n_windows=N_RDF_WINDOWS, bins=200)

        #  minimum separation (whole trajectory) · converted to σ
        min_sep = min(x["run"]["man"]["guards"]["min_separation"] for x in runs)
        #  S(k) over the late window — same box shape, so A values are comparable
        late = np.concatenate([x["run"]["traj"][(x["run"]["t_star"] >= lo)
                                               & (x["run"]["t_star"] < hi)]
                               for x in runs], axis=0)
        sk = structure_factor(late, Lx=runs[0]["run"]["Lx"],
                              Ly=runs[0]["run"]["Ly"], n_max=18)

        metrics[f"A{A:g}"] = {
            "amplitude": A, "n_seeds": len(runs),
            "zahn": zahn_phase(A),
            "window": [lo, hi],
            "dt_star": runs[0]["run"]["cfg"]["dt_star"],
            **{m: {"mean": a.mean, "se": a.se, "values": a.values}
               for m, a in agg.items()},
            "coord_kinds": {"mean": float(np.mean(kinds)), "values": kinds},
            "coord_kinds_aggregate": kinds_aggregate,
            "coordination_hist_aggregate": {
                int(z): float(v) for z, v in
                zip(sers[0].coord_labels, agg_frac) if v > 0.0},
            "five_seven_balance": {"mean": float(np.mean(bal)), "values": bal},
            "t90_tau_d": t90,
            "t90_target": target,
            "psi6_at_first_frame": float(psi_mean[0]),
            "min_separation_d": min_sep,
            "min_separation_over_sigma": min_sep / geo["sigma_over_d"],
            "sixfold_modulation": sk.sixfold_modulation,
            "rdf_first_peak_r": windows_by_A[A].first_peak_r.tolist(),
            "rdf_first_peak_g": windows_by_A[A].first_peak_g.tolist(),
            "rdf_window_t": [windows_by_A[A].t_lo.tolist(),
                             windows_by_A[A].t_hi.tolist()],
            "wall_s_total": sum(x["run"]["man"]["manifest"]["wall_s"]
                                for x in runs),
        }

        #  --- P4: does the initial-condition shell (min_sep = 0.8 d) fill in? ---
        w = windows_by_A[A]
        i_half = int(np.argmin(np.abs(w.r - 0.5)))
        metrics[f"A{A:g}"]["g_at_0.5d_by_window"] = w.g[:, i_half].tolist()
        i_08 = int(np.argmin(np.abs(w.r - 0.75)))
        metrics[f"A{A:g}"]["g_at_0.75d_by_window"] = w.g[:, i_08].tolist()

    return {"metrics": metrics, "series_by_A": series_by_A,
            "windows_by_A": windows_by_A, "energy_by_A": energy_by_A,
            "by_A": by_A}


# =============================================================================
# S6 — figures
# =============================================================================
def figures(rd: RunDir, res: dict, geo: dict, early: dict | None = None
            ) -> FigureSet:
    fs = FigureSet(rd.figs)
    plot_structure_timeseries(fs, res["series_by_A"], tau_d_si=geo["tau_d_si"],
                              energy_by_A=res["energy_by_A"])
    plot_rdf_evolution(fs, res["windows_by_A"], tau_d_si=geo["tau_d_si"],
                       sigma_over_d=geo["sigma_over_d"])

    for i, (A, runs) in enumerate(sorted(res["by_A"].items())):
        r = runs[0]["run"]                       # one seed's timelapse, not the batch
        t = r["t_star"]
        idx = [int(np.argmin(np.abs(t - tv))) for tv in VORONOI_TIMES if tv > 0]
        frames = [voronoi_frame(r["init_pos"], Lx=r["Lx"], Ly=r["Ly"])]
        times = [0.0]
        for j in idx:
            frames.append(voronoi_frame(r["traj"][j], Lx=r["Lx"], Ly=r["Ly"]))
            times.append(float(t[j]))
        #  ★ Put the per-frame fluctuation width in the caption -- so a
        #    panel-to-panel difference is not read as relaxation
        sd = float(np.mean([s.defect_fraction.std(ddof=1)
                            for s in res["series_by_A"][A]]))
        plot_voronoi_timelapse(
            fs, frames, t_star=times, tau_d_si=geo["tau_d_si"],
            L_si=geo["L_si"], amplitude=A,
            mean_defect_fraction=res["metrics"][f"A{A:g}"]["defect_fraction"]["mean"],
            defect_frame_sd=sd, name=f"0{3+i}_voronoi_A{A:g}")

    if early:
        plot_early_transient(fs, early["series_by_A"], tau_d_si=geo["tau_d_si"],
                             fits=early["fits"])
    else:
        fs.skip("06_early_transient",
                "no early-transient pass (raw_early/) — run it with `--early`")

    fs.skip("snapshots_3d", "a 2D system — 3D ray tracing (fresnel) is not installed")
    fs.skip("msd", "this card's target dynamics is the equilibrium structure. No "
                   "transport quantity was asked (§2's observable list has no MSD)")
    return fs


# =============================================================================
# S7 — seal verification + prediction comparison
# =============================================================================
def check_predictions(pred: dict, metrics: dict, geo: dict) -> list[dict]:
    rows = []
    for item in pred["items"]:
        q = item["quantity"]
        measured, note = None, ""
        if "__" in q and q.split("__")[0] in ("psi6_global", "defect_fraction",
                                              "energy_per_particle"):
            base, akey = q.split("__")[0], q.split("__")[1]
            field = {"psi6_global": "psi6_global",
                     "defect_fraction": "defect_fraction",
                     "energy_per_particle": "energy_pp"}[base]
            m = metrics[akey]
            measured = m[field]["mean"]
            #  ★ The tolerance was built from **the prior run's SE alone** (the
            #    new SE was unknowable at sealing time). Re-measured against the
            #    real SE_diff, how many σ it comes to is recorded alongside --
            #    "it exceeded the tolerance" and "it differs significantly" are
            #    two different statements.
            se_new = m[field]["se"]
            se_prior = float(str(item["tolerance"]).lstrip("±")) / (3.0 * 2 ** 0.5)
            se_diff = float(np.hypot(se_new, se_prior))
            sig = abs(measured - float(item["value"])) / se_diff if se_diff else None
            note = (f"SE_new = {se_new:.4g} ({m['n_seeds']} seeds) · "
                    f"SE_prior = {se_prior:.4g} · **actually {sig:.2f}σ**"
                    if sig is not None else f"SE = {se_new:.4g}")
        elif q.startswith("coord_kinds__"):
            akey = q.split("__")[1]
            measured = metrics[akey]["coord_kinds"]["mean"]
            note = f"per seed {metrics[akey]['coord_kinds']['values']}"
        elif q.startswith("min_separation_over_sigma__"):
            akey = q.split("__")[1]
            measured = metrics[akey]["min_separation_over_sigma"]
            note = f"= {metrics[akey]['min_separation_d']:.4g} d"
        elif q == "rdf_at_0.5d__A0.1__first_window":
            measured = metrics["A0.1"]["g_at_0.5d_by_window"][0]
            note = "first time window g(0.5 d)"
        elif q == "t90_ordering":
            t = {k: metrics[k]["t90_tau_d"] for k in metrics}
            measured = t
            note = " · ".join(f"{k}: {v:.3g} τ_d" for k, v in t.items())

        rows.append({"quantity": q, "predicted": item["value"],
                     "tolerance": item["tolerance"], "measured": measured,
                     "note": note, "basis": item["basis"],
                     "verdict": verdict_for(item, measured)})
    return rows


def verdict_for(item: dict, measured) -> str:
    """PASS/FAIL/INCONCLUSIVE. **It only proposes a verdict** — a human confirms."""
    tol, pred = str(item["tolerance"]).strip(), item["value"]
    if measured is None:
        return "NOT_EVALUATED"
    if tol.startswith("±"):
        try:
            width = float(tol[1:])
        except ValueError:
            return "NOT_EVALUATED"
        return "PASS" if abs(float(measured) - float(pred)) <= width else "FAIL"
    if tol.startswith(">"):
        return "PASS" if float(measured) > float(tol[1:]) else "FAIL"
    if tol.startswith("<"):
        return "PASS" if float(measured) < float(tol[1:]) else "FAIL"
    if isinstance(measured, dict):                     # the t90 ordering prediction
        keys = sorted(measured, key=lambda k: float(k[1:]))
        vals = [measured[k] for k in keys]
        return "PASS" if vals[-1] >= max(vals[:-1]) else "FAIL"
    return "NOT_EVALUATED"


def main() -> int:
    policy = load_policy()
    geo = geometry()
    thresholds, rows = gate_table(policy)
    print_gates(geo, thresholds, rows, policy)
    if "--gates" in sys.argv:
        return 0

    rd = RunDir.create(REPO / "runs", RUN_ID)
    analyze_only = "--analyze-only" in sys.argv
    #  --early: with the main pass already there, add only the early-transient pass
    early_only = "--early" in sys.argv

    if early_only:
        batch = run_batch(rd, rows, policy, early=True)
        #  the measured wall goes in the batch summary (`05_run_manifest.json`) --
        #  so the report need not re-read markdown to compare cost against the
        #  main pass
        if rd.exists("manifest"):
            man = json.loads(rd.file("manifest").read_text())
            man["early_batch"] = batch
            rd.write_json("manifest", man)
        #  07b is a markdown document — `write_json` would put JSON inside the `.md`
        rd.write("sensitivity", "\n".join([
            "# S7b — the early-transient pass (sampling-interval sensitivity)", "",
            f"The main pass's stride (`{TOTAL_TAU/N_FRAMES:.3g} τ_d`) is coarser "
            f"than the relaxation time, so it swallowed the transient inside the "
            f"first frame. This is the record of running a short, dense pass "
            f"separately on the same seeds.", "",
            "| item | value |", "|---|---|",
            f"| runs | {batch['n_jobs']} (`A` {len(rows)} × "
            f"{len(EARLY_SEEDS)} seeds) |",
            f"| `prod_tau` | `{batch['prod_tau']:g} τ_d` (main pass "
            f"`{TOTAL_TAU:g}`) |",
            f"| frames | `{batch['n_frames']}` → stride "
            f"`{batch['prod_tau']/batch['n_frames']:.4g} τ_d` |",
            f"| concurrency | {batch['concurrency']} |",
            f"| batch wall | `{batch['batch_wall_s']:.1f} s` |",
            f"| failed | {len(batch['failed'])} |", "",
            "> **Why the seed count went to 16** — this pass is 4 seconds for 12 "
            "runs. The `A` dependence of the relaxation time is this run's only "
            "new claim, so that is where the seeds go. The main pass cannot make "
            "the same choice, because `A=10` costs `65 s/run` "
            "(CLAUDE.md §\"error bars are free\").", "",
            "Result table and relaxation fit: "
            "[`07_validation.md`](07_validation.md) §4b · figure "
            "[`figs/06_early_transient.png`](figs/06_early_transient.png)", "",
            "Basis: [[coarse-sampling-hides-the-whole-transient]]", "",
        ]))

    if not analyze_only and not early_only:
        # --- seal: **before** the run ---
        shutil.copy2(SRC / "01_intake.md", rd.file("intake"))
        pred = yaml.safe_load((SRC / "prediction.yaml").read_text())
        rd.write("prediction", (SRC / "prediction.yaml").read_text())
        rd.write_json("prediction_json", pred)
        #  ★ Seal the driver **by hash, not by name**. With only the name, a later
        #    change to the script leaves the spec unchanged, so "what made this
        #    run" becomes unanswerable (an actual hole in this run, 2026-07-29).
        rd.write_json("spec", {"source": "scripts/soft2d_time_series.py",
                               "provenance": provenance(__file__),
                               "geometry": geo, "thresholds": thresholds,
                               "gates": rows, "seeds": list(SEEDS),
                               "total_tau": TOTAL_TAU, "n_frames": N_FRAMES,
                               "early_tau": EARLY_TAU,
                               "early_frames": EARLY_FRAMES,
                               "early_seeds": list(EARLY_SEEDS),
                               "box_shape": BOX_SHAPE, "init": INIT})
        seal = write_seal(rd)
        print(f"\n  🔒 sealed {seal.name} — "
              f"{len(seal.read_text().splitlines())} documents")
        batch = run_batch(rd, rows, policy)
        rd.write_json("manifest", batch)

    pred = yaml.safe_load(rd.file("prediction").read_text())
    print("\n## S7 — time-resolved analysis")
    res = analyze(rd, geo)

    early = None
    if (rd.path / "raw_early").exists():
        print("\n## S7b — early transient (the dense pass)")
        early = analyze_early(rd, geo) or None

    fs = figures(rd, res, geo, early)
    rd.write("figures", fs.figures_md())
    metrics = dict(res["metrics"])
    if early:
        metrics["_early_transient"] = early["summary"]
    #  ★ Provenance at analysis time. The trajectory manifest's hash is from
    #    *trajectory-generation* time, and the analysis can be re-run separately
    #    with `--analyze-only` — without this block there is no way to tell which
    #    analysis code and which freud produced metrics.json.
    metrics["_provenance_at_analysis"] = provenance(__file__)
    rd.write_json("metrics", metrics)

    checks = check_predictions(pred, res["metrics"], geo)
    rd.write_json("prediction_json", {"items": pred["items"], "checks": checks})

    print("\n## prediction comparison (the verdict is a proposal — "
          "confirmed_by: null)")
    hdr = f"{'quantity':<44} {'predicted':>14} {'measured':>14} {'tol':>12} verdict"
    print(hdr); print("-" * (len(hdr) + 8))
    for c in checks:
        mv = c["measured"]
        ms = (f"{mv:.5g}" if isinstance(mv, (int, float))
              else ("dict" if isinstance(mv, dict) else "—"))
        pv = c["predicted"]
        ps = f"{pv:.5g}" if isinstance(pv, (int, float)) else str(pv)[:14]
        print(f"{c['quantity']:<44} {ps:>14} {ms:>14} "
              f"{str(c['tolerance']):>12} {c['verdict']}")

    n_pass = sum(1 for c in checks if c["verdict"] == "PASS")
    n_fail = sum(1 for c in checks if c["verdict"] == "FAIL")
    print(f"\n  PASS {n_pass} · FAIL {n_fail} · "
          f"not evaluated {len(checks) - n_pass - n_fail}")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
