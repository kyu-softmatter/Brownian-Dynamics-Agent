"""soft-r3 finite-size ladder (S17) -- verify the **functional form** of
`psi6 ~ N^{-p}`.

usage:
  python scripts/soft2d_fss.py --gates
  python scripts/soft2d_fss.py
  python scripts/soft2d_fss.py --analyze-only

## Why a new run -- `r_cut` has to be held fixed

`runs/2026-07-29_soft-r3-nconv` obtained `p` from two points (`N=100, 256`).
Verifying the form needs three or more, but shrinking the box **shrinks `r_cut`
along with it**:

    N     L*   r_cut   βU(r_cut) at A=10
   64   8.00   3.820   0.1794
  100  10.00   4.800   0.0904
  256  16.00   7.740   0.0216
  400  20.00   9.700   0.0110

=> building the ladder with the natural `r_cut` mixes an **8-fold change in truncation
error** into the `psi6(N)` slope, and then the measured `p` cannot be told apart from
a truncation-error trend.

**Solution: fix `r_cut = 3.80` across every `N`.** That is the value the smallest box
(`N=64`, `L/2 = 4.0`) permits. The truncation error grows to `0.182 kT`, but it is
**identical at every `N`**, so it cannot manufacture a false slope. `S16` showed the
observables agree within `3 sigma` across a 4.2-fold change in truncation error, which
defends this choice -- and whether this run's `N=256` point matches the natural-`r_cut`
run checks it once more (an `8-fold` excursion).

## Run length

`prod_tau = 30 tau_d`, analysis window `[20, 30]`. Since
`tau_relax ~ 0.098 tau_d` (§8.4), `20 tau_d` is **200x** the relaxation time --
sufficient. Using `80 tau_d` would push `A=10` at `N=400` past the budget
(`600 s/run`).
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

from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,   # noqa: E402
                                       LIQUID_EXPONENT_P, hex_order_series,
                                       psi6_finite_size_exponent, zahn_phase)
from simbot.analysis.trap import aggregate_seeds                    # noqa: E402
from simbot.build import square_box_for                             # noqa: E402
from simbot.io import RunDir, provenance, write_seal                # noqa: E402
from simbot.nondim import dt_max_force, dt_max_thermal              # noqa: E402
from simbot.policy import load_policy                               # noqa: E402
from simbot.run import (Soft2DRunConfig, measure_max_force_soft2d,  # noqa: E402
                        run_soft2d)
from simbot.viz import FigureSet, plot_finite_size_scaling          # noqa: E402

import soft2d_time_series as TS                                    # noqa: E402

RUN_ID = "2026-07-29_soft-r3-fss"
SRC = REPO / "examples" / "soft-r3-fss"
DRIVERS = [Path(__file__), Path(TS.__file__)]

N_LADDER = (64, 144, 256, 400)
AMPLITUDES = (0.1, 1.0, 10.0)
#  ★ Seeds were selected so that **initial placement succeeds at every N** (a paired
#    design).
#    `min_sep = 0.8 d` rejection sampling fails for some seeds -- measured success
#    rate (seeds 31-90):
#      N=64 98.3 % · N=144 95.0 % · N=256 100 % · N=400 100 %
#    (With 60 samples the N dependence is not significant. Seed 31 failed at N=144.)
#    Different seeds per N would mix different initial-placement ensembles into the
#    psi6(N) comparison -> they must be paired.
SEEDS = (32, 33, 34, 35)
R_CUT_FIXED = 3.80                 # ★ common to every N. The value N=64's L/2 = 4.0
                                   # permits
PROD_TAU, N_FRAMES = 30.0, 300     # stride 0.1 tau_d
WINDOW = (20.0, 30.0)
CHI2_MAX = 3.0                     # threshold for reading the form (pre-registered)


def gate_table(policy) -> tuple[dict, list[dict]]:
    ts = policy.timestep
    d_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    d_F = float(ts.get("max_force_displacement_sigma", 0.005))
    dt_th = dt_max_thermal(d_th, sigma=1.0, D0=1.0)
    rows = []
    for A in AMPLITUDES:
        for N in N_LADDER:
            half = square_box_for(N) / 2.0
            if R_CUT_FIXED + TS.NLIST_BUFFER > half:
                raise SystemExit(f"⛔ N={N}: r_cut+buffer > L/2 = {half}")
            probe = Soft2DRunConfig(
                amplitude=A, n_particles=N, init=TS.INIT,
                box_shape=TS.BOX_SHAPE, r_cut=R_CUT_FIXED, r_min=TS.R_MIN,
                nlist_buffer=TS.NLIST_BUFFER, min_sep_init=TS.MIN_SEP_INIT,
                seed=SEEDS[0], dt_star=1e-6)
            m = measure_max_force_soft2d(probe)
            dt_f = dt_max_force(d_F, sigma=1.0, gamma=1.0,
                                max_force=m["max_force_star"])
            act = {k: v for k, v in (("thermal_displacement", dt_th),
                                     ("force_displacement", dt_f))
                   if v is not None}
            dom = min(act, key=act.get)
            dt = act[dom]
            rows.append({"amplitude": A, "n_particles": N, "dt_star": dt,
                         "dominant_gate": dom, "steps": int(round(PROD_TAU / dt)),
                         "max_force_star": m["max_force_star"],
                         "r_cut": m["r_cut"],
                         "beta_u_at_rcut": m["beta_u_at_rcut"],
                         "Lx": m["Lx"]})
    return {"delta_thermal": d_th, "delta_force": d_F,
            "r_cut_fixed": R_CUT_FIXED}, rows


def print_gates(rows: list[dict], policy) -> None:
    lam = float(policy.get("hardware.throughput_particle_steps_per_s", 6.3e6))
    k = policy.concurrency("default")
    eff = policy.efficiency(k)
    print(f"## gates -- `r_cut = {R_CUT_FIXED}` fixed across every N . "
          f"prod {PROD_TAU:g} tau_d / window {WINDOW}\n")
    hdr = (f"{'A':>6} {'N':>5} {'L*':>6} {'r_cut':>7} {'βU(rc)':>8} "
           f"{'max|F*|':>9} {'governing':<22} {'dt*':>10} {'steps':>9} {'est':>7}")
    print(hdr); print("-" * len(hdr))
    worst, bu = 0.0, set()
    for r in rows:
        wall = r["n_particles"] * r["steps"] / (lam * eff) * 3.4
        worst = max(worst, wall)
        bu.add(round(r["beta_u_at_rcut"], 6))
        print(f"{r['amplitude']:>6g} {r['n_particles']:>5d} {r['Lx']:>6.2f} "
              f"{r['r_cut']:>7.3f} {r['beta_u_at_rcut']:>8.4f} "
              f"{r['max_force_star']:>9.2f} {r['dominant_gate']:<22} "
              f"{r['dt_star']:>10.3g} {r['steps']:>9d} {wall:>6.0f}s")
    #  ★ There must be one truncation error per A, independent of N -- that is the
    #    whole point of this design
    per_A = {A: {round(r["beta_u_at_rcut"], 6) for r in rows
                 if r["amplitude"] == A} for A in AMPLITUDES}
    for A, s in per_A.items():
        if len(s) != 1:
            raise SystemExit(f"⛔ A={A}: beta*U(r_cut) differs across N {s} -- "
                             f"the fixed r_cut has broken")
    print(f"\n  ✅ beta*U(r_cut) is independent of N at each A: "
          f"{ {A: s.pop() for A, s in per_A.items()} }")
    budget = policy.wall_budget_s
    n = len(rows) * len(SEEDS)
    print(f"  {n} runs ({len(AMPLITUDES)} amplitudes x {len(N_LADDER)} N x "
          f"{len(SEEDS)} seeds) . concurrency {k} . longest run {worst:.0f} s "
          f"{'<=' if worst <= budget else '>'} budget {budget:.0f} s")
    if worst > budget:
        raise SystemExit("⛔ over budget -- reporting without running")


def _one(args) -> dict:
    cfg, outdir = args
    return run_soft2d(Soft2DRunConfig(**cfg), outdir=Path(outdir))


def run_batch(rd: RunDir, rows: list[dict], policy) -> dict:
    jobs = []
    for r in rows:
        for s in SEEDS:
            label = f"A{r['amplitude']:g}_N{r['n_particles']}_s{s}"
            cfg = Soft2DRunConfig(
                amplitude=r["amplitude"], n_particles=r["n_particles"],
                init=TS.INIT, box_shape=TS.BOX_SHAPE, r_cut=R_CUT_FIXED,
                r_min=TS.R_MIN, nlist_buffer=TS.NLIST_BUFFER,
                min_sep_init=TS.MIN_SEP_INIT, dt_star=r["dt_star"],
                equil_tau=0.0, prod_tau=PROD_TAU, n_frames=N_FRAMES,
                seed=s, label=label)
            jobs.append((asdict(cfg), str(rd.raw / label)))
    k = policy.concurrency("default")
    print(f"\n## S5 -- {len(jobs)} runs (concurrency {k})")
    t0 = time.perf_counter()
    done, failed = [], []
    with ProcessPoolExecutor(max_workers=k) as ex:
        futs = {ex.submit(_one, j): j[0]["label"] for j in jobs}
        for i, fut in enumerate(as_completed(futs), 1):
            label = futs[fut]
            try:
                out = fut.result()
                done.append({"label": label, "wall_s": out["wall_s"],
                             "fails": out["guards"]["failures"]})
                if out["guards"]["failures"]:
                    print(f"  ⚠ {label}: {out['guards']['failures']}")
            except Exception as e:                          # noqa: BLE001
                failed.append({"label": label, "error": repr(e)})
                print(f"  ⛔ {label}: {e!r}")
            if i % 12 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} …")
    wall = time.perf_counter() - t0
    print(f"  batch wall {wall:.1f} s . failed {len(failed)}")
    return {"done": len(done), "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "seeds": list(SEEDS),
            "n_ladder": list(N_LADDER), "r_cut_fixed": R_CUT_FIXED}


def analyze(rd: RunDir) -> dict:
    lo, hi = WINDOW
    out: dict = {}
    for A in AMPLITUDES:
        Ns, psi, se, loc, defect, defect_se = [], [], [], [], [], []
        for N in N_LADDER:
            dirs = sorted(p for p in rd.raw.glob(f"A{A:g}_N{N}_s*")
                          if (p / "samples.npz").exists())
            if not dirs:
                raise SystemExit(f"⛔ no runs for A={A} N={N}")
            per = []
            for d in dirs:
                z = np.load(d / "samples.npz")
                cfg = json.loads((d / "manifest.json").read_text())["config"]
                stride = int(z["stride"][0])
                t = np.arange(1, z["traj"].shape[0] + 1) * stride * cfg["dt_star"]
                s = hex_order_series(z["traj"], Lx=float(z["box"][0]),
                                     Ly=float(z["box"][1]), t_star=t,
                                     coord_range=(3, 12))
                m = (t >= lo) & (t < hi)
                per.append({"g": float(s.psi6_global[m].mean()),
                            "l": float(s.psi6_local[m].mean()),
                            "d": float(s.defect_fraction[m].mean())})
            ag = {k_: aggregate_seeds([p[k_] for p in per])
                  for k_ in ("g", "l", "d")}
            Ns.append(N); psi.append(ag["g"].mean); se.append(ag["g"].se)
            loc.append(ag["l"].mean)
            defect.append(ag["d"].mean); defect_se.append(ag["d"].se)
        fit = psi6_finite_size_exponent(np.array(Ns, dtype=float),
                                        np.array(psi), np.array(se))
        out[f"A{A:g}"] = {
            "amplitude": A, "zahn": zahn_phase(A), "n_ladder": Ns,
            "psi6_global": psi, "psi6_global_se": se, "psi6_local": loc,
            "defect_fraction": defect, "defect_fraction_se": defect_se,
            "n_seeds": len(dirs),
            "finite_size": {
                "p": fit.p, "p_se": fit.p_se, "eta6": fit.eta6,
                "eta6_se": fit.eta6_se, "reading": fit.reading,
                "n_points": fit.n_points, "chi2_reduced": fit.chi2_reduced,
                "residuals": list(fit.residuals), "amplitude": fit.amplitude,
                "form_verdict": fit.form_verdict(CHI2_MAX),
                "form_is_testable": fit.form_is_testable},
        }
        print(f"  A={A:<5g} p = {fit.p:.3f}±{fit.p_se:.3f}  η₆ = {fit.eta6:.2f}"
              f"±{fit.eta6_se:.2f}  χ²/dof = {fit.chi2_reduced:.2f}  "
              f"{fit.reading}\n           {fit.form_verdict(CHI2_MAX)}")
    return out


def figures(rd: RunDir, metrics: dict) -> FigureSet:
    fs = FigureSet(rd.figs)
    data, local = {}, {}
    for A in AMPLITUDES:
        m = metrics[f"A{A:g}"]
        f = m["finite_size"]
        from simbot.analysis.structure import FiniteSizeExponent
        data[A] = {"N": m["n_ladder"], "psi6": m["psi6_global"],
                   "se": m["psi6_global_se"],
                   "fit": FiniteSizeExponent(
                       p=f["p"], p_se=f["p_se"], eta6=f["eta6"],
                       eta6_se=f["eta6_se"], n_points=f["n_points"],
                       reading=f["reading"], chi2_reduced=f["chi2_reduced"],
                       residuals=tuple(f["residuals"]),
                       amplitude=f["amplitude"])}
        local[A] = {"N": m["n_ladder"], "defect": m["defect_fraction"],
                    "defect_se": m["defect_fraction_se"]}
    plot_finite_size_scaling(
        fs, data, local_data=local, exponent_liquid=LIQUID_EXPONENT_P,
        exponent_hexatic=KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
        name="01_fss_ladder")
    fs.skip("voronoi", "the parent run already characterised the defects")
    fs.skip("early_transient", "relaxation was already closed with 1513 seeds (§8.4)")
    return fs


def check(pred: dict, metrics: dict) -> list[dict]:
    """Compare against the sealed prediction. Measurements are read from `metrics`
    only."""
    got: dict = {}
    for A in AMPLITUDES:
        m = metrics[f"A{A:g}"]
        f = m["finite_size"]
        got[f"chi2_reduced__A{A:g}"] = f["chi2_reduced"]
        got[f"psi6_exponent_p__A{A:g}"] = f["p"]
        #  Compare the N=256 point against the earlier run (natural r_cut) -- pulled
        #  out of the ladder
        i = m["n_ladder"].index(256)
        got[f"psi6_global__A{A:g}__N256"] = m["psi6_global"][i]
    f10 = metrics["A10"]["finite_size"]
    got["eta6_minus_3sigma__A10"] = f10["eta6"] - 3.0 * f10["eta6_se"]

    rows = []
    for it in pred["items"]:
        q = it["quantity"]
        v = got.get(q)
        tol = str(it["tolerance"]).strip()
        if v is None or not np.isfinite(v):
            verdict = "NOT_EVALUATED"
        elif tol.startswith(">"):
            verdict = "PASS" if v > float(tol[1:]) else "FAIL"
        elif tol.startswith("<"):
            verdict = "PASS" if v < float(tol[1:]) else "FAIL"
        else:
            verdict = ("PASS" if abs(v - float(it["value"]))
                       <= float(tol.lstrip("±")) else "FAIL")
        rows.append({"quantity": q, "predicted": it["value"],
                     "tolerance": tol, "measured": v, "verdict": verdict,
                     "discriminates": it.get("discriminates", "")})
    return rows


def residual_diagnosis(metrics: dict, A: float) -> dict:
    """Is an excess `chi^2` due to **curvature** or to **error bars being too narrow**?

    ★ Curvature makes the residual signs form a monotone pattern (e.g. `+,+,-,-`
      rather than `-,+,+,-`). Scatter makes the sign flip often.
      Since `chi^2 ~ 1/SE^2`, underestimating SE by 2x inflates `chi^2` 4x -- and a
      4-seed SE has 41 % uncertainty of its own.
      ([[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]]).
    """
    m = metrics[f"A{A:g}"]
    f = m["finite_size"]
    r = np.asarray(f["residuals"])
    se = np.asarray(m["psi6_global_se"])
    y = np.asarray(m["psi6_global"])
    z = r / (se / y)                              # log residual / log error
    signs = np.sign(r)
    flips = int(np.sum(signs[:-1] != signs[1:]))
    #  By what factor must SE grow to bring chi^2/dof to 1?
    inflate = float(np.sqrt(f["chi2_reduced"])) if np.isfinite(
        f["chi2_reduced"]) else float("nan")
    return {"amplitude": A, "z_residuals": [float(x) for x in z],
            "sign_flips": flips, "n_points": int(r.size),
            "relative_se": [float(x) for x in se / y],
            "se_inflation_for_chi2_unity": inflate,
            "reading": ("curved (a problem with the form)" if flips <= 1
                        else "scattered (points to underestimated error bars)")}


def main() -> int:
    policy = load_policy()
    thresholds, rows = gate_table(policy)
    print_gates(rows, policy)
    if "--gates" in sys.argv:
        return 0

    pred = yaml.safe_load((SRC / "prediction.yaml").read_text())
    rd = RunDir.create(REPO / "runs", RUN_ID)
    if "--analyze-only" not in sys.argv:
        shutil.copy2(SRC / "prediction.yaml", rd.file("prediction"))
        rd.write_json("prediction_json", pred)
        rd.write_json("spec", {"source": "scripts/soft2d_fss.py",
                               "provenance": provenance(DRIVERS),
                               "n_ladder": list(N_LADDER),
                               "amplitudes": list(AMPLITUDES),
                               "seeds": list(SEEDS),
                               "r_cut_fixed": R_CUT_FIXED,
                               "prod_tau": PROD_TAU, "n_frames": N_FRAMES,
                               "window": list(WINDOW), "chi2_max": CHI2_MAX,
                               "thresholds": thresholds, "gates": rows})
        seal = write_seal(rd)
        print(f"\n  🔒 sealed {seal.name} -- "
              f"{len(seal.read_text().splitlines())} documents (before running)")
        rd.write_json("manifest", run_batch(rd, rows, policy))

    print("\n## S7 -- finite-size ladder (4 points)")
    metrics = analyze(rd)
    metrics["_provenance_at_analysis"] = provenance(DRIVERS)
    fs = figures(rd, metrics)
    rd.write("figures", fs.figures_md())

    checks = check(pred, metrics)
    diag = {f"A{A:g}": residual_diagnosis(metrics, A) for A in AMPLITUDES}
    metrics["_checks"] = checks
    metrics["_residual_diagnosis"] = diag
    rd.write_json("metrics", metrics)

    print("\n## prediction comparison (the verdict is a proposal -- "
          "confirmed_by: null)")
    hdr = f"{'item':<32} {'predicted':>10} {'tolerance':>10} {'measured':>10} verdict"
    print(hdr); print("-" * (len(hdr) + 6))
    for c in checks:
        print(f"{c['quantity']:<32} {float(c['predicted']):>10.4g} "
              f"{c['tolerance']:>10} {c['measured']:>10.4g} {c['verdict']}")
    n_p = sum(1 for c in checks if c["verdict"] == "PASS")
    print(f"\n  PASS {n_p} · FAIL {len(checks) - n_p}")

    fails = [c for c in checks if c["verdict"] == "FAIL"]
    if fails:
        print("\n## FAIL diagnosis -- curvature, or error bars too narrow?")
        for c in fails:
            if not c["quantity"].startswith("chi2_reduced__"):
                continue
            A = float(c["quantity"].split("__A")[1])
            d = diag[f"A{A:g}"]
            print(f"  A={A:g}: residual/sigma = " +
                  " ".join(f"{x:+.2f}" for x in d["z_residuals"]))
            print(f"    sign flips {d['sign_flips']}/{d['n_points']-1} -> "
                  f"**{d['reading']}**")
            print(f"    relative SE = " +
                  " ".join(f"{x:.1%}" for x in d["relative_se"]))
            print(f"    inflating SE by {d['se_inflation_for_chi2_unity']:.2f}x "
                  f"brings chi^2/dof to 1 (a 4-seed SE has 41 % uncertainty itself)")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
