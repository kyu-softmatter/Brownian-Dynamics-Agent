"""soft-r3 `eta6(A)` scan -- bracket the liquid->crystal transition **in our own
system**.

usage:
  python scripts/soft2d_ascan.py --gates
  python scripts/soft2d_ascan.py
  python scripts/soft2d_ascan.py --analyze-only

## ★ Why Zahn's narrow window (`A = 10.03-10.75`) was not swept

`runs/2026-07-29_soft-r3-fss` measured `eta6` at three `A`: `A=0.1` -> `2.03`,
`A=1` -> `1.86`, `A=10` -> `1.27`. `A=10` is **`-0.3 %`** from Zahn's hexatic
boundary (`Gamma=55.87`, `A=10.03`), yet its `eta6` is **5.1x** the upper bound
`0.25`. Just below the boundary it ought to be near `0.25` -- that is a tension. And
`A=100` is already a crystal (card §8.1).
=> **the transition lies between `A = 10` and `100`, and is not a narrow window.**

## Design -- sweep `A` logarithmically

`A = 10, 13.3, 17.8, 23.7, 31.6` (quarter-decade). At each `A` the exponent comes from
an `N` ladder of `64, 144, 256, 400`. `r_cut` fixed at `3.80`, seeds `32-35` -- the
S17 convention. `prod_tau = 10 tau_d`, window `[6, 10]` (60x
`tau_relax ~ 0.098 tau_d`).
`A = 10` overlaps S17 (`30 tau_d`), which **tests the effect of the shortened run
length.**

## ★★ The scan's upper limit is set by truncation error, not by budget

`beta*U(r_cut) = A/r_cut^3`, so it is **proportional** to `A`. With `r_cut` fixed at
`3.80`, `A=100` gives `1.82 kT` -- cutting the potential where it is still of thermal
scale. `r_cut` cannot be raised because the smallest box (`N=64`, `L/2 = 4.0`) bounds
it.
=> only up to `A <= 31.6` (`beta*U <= 0.58`). Beyond that **a larger box is
required.**

## ★★★ `eta6 <= 1/4` must NOT be read as hexatic

**A crystal also satisfies `eta6 ~ 0`** (`psi6 -> const`, hence `p = 0`).
The first version lacked this distinction and **reported the crystals at `A >= 13.3`
as "possibly hexatic".**
What separates them is the **magnitude** of `psi6` -- `O(1)` for a crystal, small and
slowly decaying for a hexatic.
=> `analysis.structure.phase_from_finite_size` reads both axes.

## Results summary

    A     Gamma  psi6(400)      eta6      chi2/dof  phase
   10    55.68     0.144   +1.40±0.22     0.04   isotropic liquid
   13.3  74.06     0.688   -0.19±0.04    46.9    crystal
   17.8  99.12     0.778   -0.22±0.01   106      crystal
   23.7 131.97     0.833   -0.41±0.01    24.8    crystal
   31.6 175.96     0.872   -0.25±0.01   322      crystal

Liquid->crystal lies within `A = 10-13.3` (33 % width). Zahn's crystal boundary
`Gamma = 59.88` (`A = 10.75`) is **inside that bracket**. The hexatic window (7 %
width) hides inside the bracket and is therefore **not resolved by this scan.**

The huge `chi^2` and negative `p` at `A >= 13.3` are because a power law does not hold
in a crystal (`psi6` saturates at `1`). `N=64` has `L = 8 d`, too small for a crystal
to form properly, so it comes out unusually low (`0.394` against `N=144`'s `0.672`) --
and that is what drives `p` negative. **In the crystal region this exponent is not a
physical quantity.**
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

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))
sys.path.insert(0, str(REPO / "scripts"))

from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,   # noqa: E402
                                       LIQUID_EXPONENT_P, hex_order_series,
                                       phase_from_finite_size,
                                       psi6_finite_size_exponent, zahn_phase)
from simbot.analysis.trap import aggregate_seeds                    # noqa: E402
from simbot.io import RunDir, provenance, write_seal                # noqa: E402
from simbot.nondim import dt_max_force, dt_max_thermal              # noqa: E402
from simbot.policy import load_policy                               # noqa: E402
from simbot.run import (Soft2DRunConfig, measure_max_force_soft2d,  # noqa: E402
                        run_soft2d)

import soft2d_fss as FSS                                            # noqa: E402
import soft2d_time_series as TS                                     # noqa: E402

RUN_ID = "2026-07-29_soft-r3-ascan"
SRC = REPO / "examples" / "soft-r3-ascan"
DRIVERS = [Path(__file__), Path(FSS.__file__), Path(TS.__file__)]

#  ★★ The scan's upper limit is set by **truncation error**, not by budget.
#  `beta*U(r_cut) = A/r_cut^3`, so it is **proportional** to `A`. With `r_cut` fixed
#  at 3.80:
#      A=10 → 0.18 kT · A=31.6 → 0.58 · A=56.2 → 1.02 · A=100 → 1.82
#  `A=100` cuts the potential where it is still `1.8 kT` -- the structure comes out
#  wrong. Raising `r_cut` requires a larger smallest box (`r_cut+buffer <= L/2`),
#  which means discarding the bottom of the ladder and losing the `N` lever arm.
#  => **sweep only up to `A <= 31.6`** (5 quarter-decade points). Above that a larger
#    box is needed.
#    (the budget gate also rejected A=100 at N=400, 665 s > 600 s -- the two
#     constraints point the same way)
AMPLITUDES = (10.0, 13.3, 17.8, 23.7, 31.6)       # quarter-decade, 10 -> 10^1.5
N_LADDER = FSS.N_LADDER                            # 64·144·256·400
SEEDS = FSS.SEEDS                                  # 32-35 (selected to succeed at
                                                   # every N)
R_CUT_FIXED = FSS.R_CUT_FIXED                      # 3.80
PROD_TAU, N_FRAMES = 10.0, 200                     # stride 0.05 tau_d
WINDOW = (6.0, 10.0)
CHI2_MAX = FSS.CHI2_MAX


def gate_table(policy):
    ts = policy.timestep
    dt_th = dt_max_thermal(float(ts["max_thermal_displacement_sigma"]),
                           sigma=1.0, D0=1.0)
    d_F = float(ts["max_force_displacement_sigma"])
    rows = []
    for A in AMPLITUDES:
        for N in N_LADDER:
            probe = Soft2DRunConfig(
                amplitude=A, n_particles=N, init=TS.INIT,
                box_shape=TS.BOX_SHAPE, r_cut=R_CUT_FIXED, r_min=TS.R_MIN,
                nlist_buffer=TS.NLIST_BUFFER, min_sep_init=TS.MIN_SEP_INIT,
                seed=SEEDS[0], dt_star=1e-6)
            m = measure_max_force_soft2d(probe)
            dt_f = dt_max_force(d_F, sigma=1.0, gamma=1.0,
                                max_force=m["max_force_star"])
            act = {"thermal_displacement": dt_th}
            if dt_f is not None:
                act["force_displacement"] = dt_f
            dom = min(act, key=act.get)
            dt = act[dom]
            rows.append({"amplitude": A, "n_particles": N, "dt_star": dt,
                         "dominant_gate": dom,
                         "steps": int(round(PROD_TAU / dt)),
                         "max_force_star": m["max_force_star"],
                         "r_cut": m["r_cut"],
                         "beta_u_at_rcut": m["beta_u_at_rcut"]})
    return rows


def print_gates(rows, policy) -> float:
    lam = float(policy.get("hardware.throughput_particle_steps_per_s", 6.3e6))
    k = policy.concurrency("default")
    eff = policy.efficiency(k)
    print(f"## gates -- r_cut fixed at {R_CUT_FIXED} . prod {PROD_TAU:g} tau_d / "
          f"window {WINDOW}\n")
    print(f"{'A':>7} {'Γ':>8} {'βU(rc)':>8} {'N=400 dt*':>10} "
          f"{'steps':>9} {'N=400 est':>11}")
    print("-" * 60)
    worst, total = 0.0, 0.0
    for A in AMPLITUDES:
        rs = [r for r in rows if r["amplitude"] == A]
        r4 = [r for r in rs if r["n_particles"] == 400][0]
        for r in rs:
            w = r["n_particles"] * r["steps"] / (lam * eff) * 3.4
            worst = max(worst, w)
            total += w * len(SEEDS)
        w4 = 400 * r4["steps"] / (lam * eff) * 3.4
        print(f"{A:>7g} {zahn_phase(A)['gamma']:>8.2f} "
              f"{r4['beta_u_at_rcut']:>8.4f} {r4['dt_star']:>10.3g} "
              f"{r4['steps']:>9d} {w4:>10.0f}s")
    budget = policy.wall_budget_s
    n = len(rows) * len(SEEDS)
    print(f"\n  {n} runs (A {len(AMPLITUDES)} x N {len(N_LADDER)} x seeds "
          f"{len(SEEDS)})")
    print(f"  longest run {worst:.0f} s {'<=' if worst <= budget else '>'} budget "
          f"{budget:.0f} s . total work {total/60:.0f} min -> at concurrency {k}, "
          f"about {total/k/60:.0f} min")
    if worst > budget:
        raise SystemExit("⛔ over budget -- reporting without running")
    return total


def _one(args):
    cfg, outdir = args
    return run_soft2d(Soft2DRunConfig(**cfg), outdir=Path(outdir))


def run_batch(rd: RunDir, rows, policy) -> dict:
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
            if i % 20 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} … ({time.perf_counter()-t0:.0f} s)")
    wall = time.perf_counter() - t0
    print(f"  batch wall {wall:.1f} s . failed {len(failed)}")
    return {"done": len(done), "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "seeds": list(SEEDS),
            "amplitudes": list(AMPLITUDES), "n_ladder": list(N_LADDER)}


def analyze(rd: RunDir) -> dict:
    lo, hi = WINDOW
    out: dict = {}
    for A in AMPLITUDES:
        Ns, psi, se, defect, defect_se, loc = [], [], [], [], [], []
        for N in N_LADDER:
            dirs = sorted(p for p in rd.raw.glob(f"A{A:g}_N{N}_s*")
                          if (p / "samples.npz").exists())
            if not dirs:
                raise SystemExit(f"⛔ missing A={A} N={N}")
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
            "finite_size": {
                "p": fit.p, "p_se": fit.p_se, "eta6": fit.eta6,
                "eta6_se": fit.eta6_se, "reading": fit.reading,
                "chi2_reduced": fit.chi2_reduced,
                "residuals": list(fit.residuals),
                "amplitude_fit": fit.amplitude,
                "form_verdict": fit.form_verdict(CHI2_MAX)},
        }
        #  ★ The verdict is never from the exponent alone -- a crystal also satisfies
        #    eta6 ~ 0. The magnitude of psi6 must be read alongside (on 2026-07-29,
        #    lacking this distinction, crystals were reported as "possibly hexatic").
        ph = phase_from_finite_size(fit, psi[-1])
        out[f"A{A:g}"]["phase"] = ph
        print(f"  A={A:<6g} Γ={zahn_phase(A)['gamma']:>7.2f}  "
              f"ψ₆(400)={psi[-1]:.4f}  p={fit.p:+.3f}±{fit.p_se:.3f}  "
              f"η₆={fit.eta6:+.3f}±{fit.eta6_se:.3f}  "
              f"χ²/dof={fit.chi2_reduced:>7.2f}  **{ph['phase']}**")
        if ph["exponent_alone_would_say"] == "hexatic-or-below" \
                and ph["phase"] == "crystal":
            print(f"           ⚠ by exponent alone this reads as 'at or below "
                  f"hexatic' -- the psi6 magnitude settled it as crystal")
    return out


def transition_bracket(metrics: dict) -> dict:
    """Find the `A` **interval** across which the phase changes.

    ★ An earlier version looked for where `eta6` crosses `0.25` and called that the
      "hexatic boundary". **That was wrong** -- a crystal also has `eta6 ~ 0`, so that
      crossing was the liquid->crystal boundary.
      The phase is read by `phase_from_finite_size` on two axes: the exponent and the
      `psi6` magnitude.
    """
    As = [metrics[f"A{A:g}"]["amplitude"] for A in AMPLITUDES]
    ph = [metrics[f"A{A:g}"]["phase"]["phase"] for A in AMPLITUDES]
    eta = [metrics[f"A{A:g}"]["finite_size"]["eta6"] for A in AMPLITUDES]
    psi = [metrics[f"A{A:g}"]["psi6_global"][-1] for A in AMPLITUDES]
    out = {"amplitudes": As, "phases": ph, "eta6": eta,
           "psi6_at_largest_N": psi, "ceiling": KTHNY_ETA6_HEXATIC_LIQUID}
    changes = [(As[i], As[i + 1], ph[i], ph[i + 1])
               for i in range(len(As) - 1) if ph[i] != ph[i + 1]]
    out["transitions"] = [{"A_lo": a, "A_hi": b, "from": f, "to": t,
                           "width_pct": 100.0 * (b - a) / a}
                          for a, b, f, t in changes]
    if not changes:
        out["note"] = f"the entire scan is {ph[0]} -- the transition lies outside it"
    else:
        out["note"] = " · ".join(
            f"{f} -> {t} between A = {a:g} and {b:g} ({100*(b-a)/a:.0f} % width)"
            for a, b, f, t in changes)
    #  was a hexatic observed?
    out["hexatic_observed"] = any(x == "hexatic-candidate" for x in ph)
    return out


def main() -> int:
    policy = load_policy()
    rows = gate_table(policy)
    print_gates(rows, policy)
    if "--gates" in sys.argv:
        return 0

    rd = RunDir.create(REPO / "runs", RUN_ID)
    if "--analyze-only" not in sys.argv:
        if (SRC / "prediction.yaml").exists():
            shutil.copy2(SRC / "prediction.yaml", rd.file("prediction"))
        rd.write_json("spec", {"source": "scripts/soft2d_ascan.py",
                               "provenance": provenance(DRIVERS),
                               "parent_run": "runs/2026-07-29_soft-r3-fss",
                               "amplitudes": list(AMPLITUDES),
                               "n_ladder": list(N_LADDER),
                               "seeds": list(SEEDS),
                               "r_cut_fixed": R_CUT_FIXED,
                               "prod_tau": PROD_TAU, "window": list(WINDOW),
                               "n_frames": N_FRAMES, "gates": rows,
                               "eta6_ceiling": KTHNY_ETA6_HEXATIC_LIQUID})
        seal = write_seal(rd)
        print(f"\n  🔒 sealed {seal.name} -- "
              f"{len(seal.read_text().splitlines())} documents (before running)")
        rd.write_json("manifest", run_batch(rd, rows, policy))

    print(f"\n## S7 -- eta6(A) scan (ladder {list(N_LADDER)}, window {WINDOW})")
    metrics = analyze(rd)
    cr = transition_bracket(metrics)
    metrics["_transition_bracket"] = cr
    metrics["_provenance_at_analysis"] = provenance(DRIVERS)
    rd.write_json("metrics", metrics)
    print(f"\n## phase-transition bracket")
    print(f"  {cr['note']}")
    for t in cr["transitions"]:
        g_lo = zahn_phase(t["A_lo"])["gamma"]
        g_hi = zahn_phase(t["A_hi"])["gamma"]
        print(f"  → A = {t['A_lo']:g}–{t['A_hi']:g}  (Γ = {g_lo:.1f}–{g_hi:.1f})")
        print(f"    Zahn's crystal boundary Gamma = 59.88 (A = 10.75) -> "
              f"{'**inside** the bracket' if g_lo <= 59.88 <= g_hi else 'outside it'}")
    print(f"  hexatic observed? "
          f"{'yes' if cr['hexatic_observed'] else '**no**'} -- "
          f"Zahn's window A = 10.03-10.75 (7 % width) hides inside the bracket above")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
