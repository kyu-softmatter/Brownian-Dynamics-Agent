"""soft-r3 hexatic 창 (S25) + 고-`A` 절단 해방 (S26) — **큰 사다리**로 둘을 함께.

usage:
  python scripts/soft2d_hexwin.py --gates
  python scripts/soft2d_hexwin.py
  python scripts/soft2d_hexwin.py --analyze-only

## 왜 한 런인가 — 큰 상자가 두 미결을 동시에 푼다

앞 스캔(`soft-r3-ascan`)의 사다리는 `N = 64–400` 이었고 `r_cut` 이 최소 상자에
묶여 `3.80` 이었다. 그래서 둘이 막혔다:

- **S25 (hexatic 창)** — Zahn 창 `A = 10.03–10.75` (폭 `7 %`) 이 우리 브래킷
  `A = 10–13.3` (폭 `33 %`) 안에 숨었다. `A` 를 촘촘히 하려면 `ψ₆` 오차가 작아야 하고,
  그러려면 `N` 이 커야 한다.
- **S26 (고 `A`)** — `βU(r_cut) = A/r_cut³` 이 `A` 에 비례해서 `A=100` 에서
  `1.82 kT` 였다. 퍼텐셜을 아직 열에너지 규모인 곳에서 자른다.

**사다리를 `N = 256·576·1024` 로 올리면 `r_cut = 7.80` 이 허용된다** (`N=256` 의
`L/2 = 8.0`). 그러면 절단오차가 **8.7배** 줄어든다:

    A       r_cut=3.80    r_cut=7.80
   10.4        0.190          0.022
   31.6        0.576          0.067
  100          1.822          0.211

⇒ `A = 10.0–10.8` 을 `0.2` 간격으로 훑고(S25), `A = 31.6` 을 제대로 된 절단으로
  다시 본다(S26 부분).

## 왜 `A = 100` 은 사다리로 못 보는가 (S26 이 부분만 닫히는 이유)

`dt ∝ 1/A` 이므로 `steps ∝ A` 다. `A=100`·`N=1024` 는 런 하나가 `~2100 s` 로
예산(`600 s`)의 3.5배다. **이 기계에서 `A=100` 의 유한크기 사다리는 불가능하다** —
`A = 31.6` 까지가 한계다. 그 사실을 결과로 기록한다.

## 런 길이 — `10 τ_d` 를 쓰는 근거

`soft-r3-ascan` 의 S23 이 `10 τ_d` 로 `A=10` 의 `p` 를 `30 τ_d` 판본과 `0.5σ` 이내로
재현했다. 이 런의 `A` 는 모두 `10` 근처이므로 그 검증이 직접 적용된다.
**다만 가정하지 않고 창 내부 표류를 함께 보고한다** (`[6,8]` vs `[8,10]`).

## 시드

`41–44` — `N = 256·576·1024` 전부에서 초기배치 성공 (실측 20/20).
큰 상자가 오히려 쉽다 (`min_sep = 0.8 d` 기각표집).
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
                                       hex_order_series,
                                       phase_from_finite_size,
                                       psi6_finite_size_exponent, zahn_phase)
from simbot.analysis.trap import aggregate_seeds                    # noqa: E402
from simbot.io import RunDir, provenance, write_seal                # noqa: E402
from simbot.nondim import dt_max_force, dt_max_thermal              # noqa: E402
from simbot.policy import load_policy                               # noqa: E402
from simbot.run import (Soft2DRunConfig, measure_max_force_soft2d,  # noqa: E402
                        run_soft2d)

import soft2d_time_series as TS                                     # noqa: E402

RUN_ID = "2026-07-29_soft-r3-hexwin"
SRC = REPO / "examples" / "soft-r3-hexwin"
DRIVERS = [Path(__file__), Path(TS.__file__)]

#  S25: Zahn 창 10.03–10.75 를 0.2 간격으로 · S26: 31.6 을 제대로 된 절단으로
AMPLITUDES = (10.0, 10.2, 10.4, 10.6, 10.8, 31.6)
N_LADDER = (256, 576, 1024)
SEEDS = (41, 42, 43, 44)
R_CUT_FIXED = 7.80                 # N=256 의 L/2 = 8.0 이 허용하는 값
PROD_TAU, N_FRAMES = 10.0, 200
WINDOW = (6.0, 10.0)
DRIFT_SPLIT = 8.0                  # 표류 진단: [6,8] vs [8,10]
CHI2_MAX = 3.0


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
            rows.append({"amplitude": A, "n_particles": N,
                         "dt_star": act[dom], "dominant_gate": dom,
                         "steps": int(round(PROD_TAU / act[dom])),
                         "max_force_star": m["max_force_star"],
                         "r_cut": m["r_cut"],
                         "beta_u_at_rcut": m["beta_u_at_rcut"]})
    return rows


def print_gates(rows, policy) -> None:
    lam = float(policy.get("hardware.throughput_particle_steps_per_s", 6.3e6))
    k = policy.concurrency("default")
    eff = policy.efficiency(k)
    print(f"## 게이트 — r_cut = {R_CUT_FIXED} 고정 · 사다리 {list(N_LADDER)} · "
          f"prod {PROD_TAU:g} τ_d / 창 {WINDOW}\n")
    print(f"{'A':>7} {'Γ':>8} {'βU(rc)':>8} {'rc=3.8 이면':>11} "
          f"{'dt*':>10} {'N=1024 예상':>12}")
    print("-" * 62)
    worst, total = 0.0, 0.0
    for A in AMPLITUDES:
        rs = [r for r in rows if r["amplitude"] == A]
        big = max(rs, key=lambda r: r["n_particles"])
        for r in rs:
            w = r["n_particles"] * r["steps"] / (lam * eff) * 3.4
            worst = max(worst, w)
            total += w * len(SEEDS)
        wb = big["n_particles"] * big["steps"] / (lam * eff) * 3.4
        print(f"{A:>7g} {zahn_phase(A)['gamma']:>8.2f} "
              f"{big['beta_u_at_rcut']:>8.4f} {A/3.8**3:>11.4f} "
              f"{big['dt_star']:>10.3g} {wb:>11.0f}s")
    budget = policy.wall_budget_s
    print(f"\n  런 {len(rows)*len(SEEDS)}개 (A {len(AMPLITUDES)} × N "
          f"{len(N_LADDER)} × 시드 {len(SEEDS)})")
    print(f"  최장 런 {worst:.0f} s {'≤' if worst <= budget else '>'} 예산 "
          f"{budget:.0f} s · 총 작업 {total/60:.0f} 분 → 동시 {k} 에서 "
          f"약 {total/k/60:.0f} 분")
    if worst > budget:
        raise SystemExit("⛔ 예산 초과 — 실행하지 않고 보고한다")


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
    print(f"\n## S5 — {len(jobs)} 런 (동시 {k})")
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
                print(f"  {i}/{len(jobs)} … ({time.perf_counter()-t0:.0f} s)")
    wall = time.perf_counter() - t0
    print(f"  배치 wall {wall:.1f} s · 실패 {len(failed)}")
    return {"done": len(done), "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "seeds": list(SEEDS),
            "amplitudes": list(AMPLITUDES), "n_ladder": list(N_LADDER),
            "r_cut_fixed": R_CUT_FIXED}


def analyze(rd: RunDir) -> dict:
    lo, hi = WINDOW
    out: dict = {}
    for A in AMPLITUDES:
        Ns, psi, se, loc, defect, defect_se, drift = [], [], [], [], [], [], []
        for N in N_LADDER:
            dirs = sorted(p for p in rd.raw.glob(f"A{A:g}_N{N}_s*")
                          if (p / "samples.npz").exists())
            if not dirs:
                raise SystemExit(f"⛔ A={A} N={N} 없다")
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
                m1 = (t >= lo) & (t < DRIFT_SPLIT)
                m2 = (t >= DRIFT_SPLIT) & (t < hi)
                per.append({"g": float(s.psi6_global[m].mean()),
                            "l": float(s.psi6_local[m].mean()),
                            "d": float(s.defect_fraction[m].mean()),
                            "drift": float(s.psi6_global[m2].mean()
                                           - s.psi6_global[m1].mean())})
            ag = {k_: aggregate_seeds([p[k_] for p in per])
                  for k_ in ("g", "l", "d", "drift")}
            Ns.append(N); psi.append(ag["g"].mean); se.append(ag["g"].se)
            loc.append(ag["l"].mean)
            defect.append(ag["d"].mean); defect_se.append(ag["d"].se)
            drift.append({"mean": ag["drift"].mean, "se": ag["drift"].se})
        fit = psi6_finite_size_exponent(np.array(Ns, dtype=float),
                                       np.array(psi), np.array(se))
        ph = phase_from_finite_size(fit, psi[-1])
        out[f"A{A:g}"] = {
            "amplitude": A, "zahn": zahn_phase(A), "n_ladder": Ns,
            "psi6_global": psi, "psi6_global_se": se, "psi6_local": loc,
            "defect_fraction": defect, "defect_fraction_se": defect_se,
            "psi6_drift": drift, "phase": ph,
            "finite_size": {
                "p": fit.p, "p_se": fit.p_se, "eta6": fit.eta6,
                "eta6_se": fit.eta6_se, "reading": fit.reading,
                "chi2_reduced": fit.chi2_reduced,
                "residuals": list(fit.residuals),
                "form_verdict": fit.form_verdict(CHI2_MAX)},
        }
        dmax = max(abs(d["mean"]) for d in drift)
        dsig = max(abs(d["mean"]) / d["se"] if d["se"] else 0.0 for d in drift)
        print(f"  A={A:<6g} Γ={zahn_phase(A)['gamma']:>7.2f}  "
              f"ψ₆(1024)={psi[-1]:.4f}±{se[-1]:.4f}  "
              f"η₆={fit.eta6:+.3f}±{fit.eta6_se:.3f}  "
              f"χ²/dof={fit.chi2_reduced:>7.2f}  **{ph['phase']}**")
        print(f"           표류 |ψ₆([8,10])−ψ₆([6,8])| ≤ {dmax:.4f} "
              f"({dsig:.1f}σ)  결함(1024)={defect[-1]:.4f}")
    return out


def summarise(metrics: dict) -> dict:
    As = [metrics[f"A{A:g}"]["amplitude"] for A in AMPLITUDES]
    ph = [metrics[f"A{A:g}"]["phase"]["phase"] for A in AMPLITUDES]
    eta = [metrics[f"A{A:g}"]["finite_size"]["eta6"] for A in AMPLITUDES]
    eta_se = [metrics[f"A{A:g}"]["finite_size"]["eta6_se"] for A in AMPLITUDES]
    win = [(a, p) for a, p in zip(As, ph) if 10.03 <= a <= 10.75]
    out = {"amplitudes": As, "phases": ph, "eta6": eta, "eta6_se": eta_se,
           "zahn_window": [10.03, 10.75],
           "phases_inside_zahn_window": win,
           "hexatic_observed": any(p == "hexatic-candidate" for p in ph)}
    changes = [(As[i], As[i + 1], ph[i], ph[i + 1])
               for i in range(len(As) - 1) if ph[i] != ph[i + 1]]
    out["transitions"] = [{"A_lo": a, "A_hi": b, "from": f, "to": t,
                           "width_pct": 100.0 * (b - a) / a}
                          for a, b, f, t in changes]
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
        rd.write_json("spec", {
            "source": "scripts/soft2d_hexwin.py",
            "provenance": provenance(DRIVERS),
            "parent_runs": ["runs/2026-07-29_soft-r3-ascan",
                            "runs/2026-07-29_soft-r3-fss"],
            "closes": ["S25 (hexatic 창)", "S26 부분 (A=31.6 절단 해방)"],
            "amplitudes": list(AMPLITUDES), "n_ladder": list(N_LADDER),
            "seeds": list(SEEDS), "r_cut_fixed": R_CUT_FIXED,
            "prod_tau": PROD_TAU, "window": list(WINDOW),
            "drift_split": DRIFT_SPLIT, "n_frames": N_FRAMES, "gates": rows,
            "eta6_ceiling": KTHNY_ETA6_HEXATIC_LIQUID,
            "a100_infeasible": (
                "A=100·N=1024 는 dt ∝ 1/A 때문에 런 하나가 ~2100 s = 예산의 3.5배다. "
                "이 기계에서 A=100 의 유한크기 사다리는 불가능하다 — A=31.6 이 한계"),
        })
        seal = write_seal(rd)
        print(f"\n  🔒 봉인 {seal.name} — "
              f"{len(seal.read_text().splitlines())}개 문서 (실행 전)")
        rd.write_json("manifest", run_batch(rd, rows, policy))

    print(f"\n## S7 — hexatic 창 (사다리 {list(N_LADDER)}, r_cut {R_CUT_FIXED})")
    metrics = analyze(rd)
    s = summarise(metrics)
    metrics["_summary"] = s
    metrics["_provenance_at_analysis"] = provenance(DRIVERS)
    rd.write_json("metrics", metrics)

    print(f"\n## Zahn 창 {s['zahn_window']} 안의 상")
    for a, p in s["phases_inside_zahn_window"]:
        print(f"  A = {a:g} → {p}")
    print(f"\n  hexatic 관측? {'**예**' if s['hexatic_observed'] else '**아니다**'}")
    for t in s["transitions"]:
        print(f"  전이: {t['from']} → {t['to']} 가 A = {t['A_lo']:g}–"
              f"{t['A_hi']:g} 사이 ({t['width_pct']:.0f} % 폭)")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
