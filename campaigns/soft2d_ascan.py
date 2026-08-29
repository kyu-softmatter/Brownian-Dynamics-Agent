"""soft-r3 `η₆(A)` 스캔 — 액체→결정 전이를 **우리 계에서** 브래킷한다.

usage:
  python scripts/soft2d_ascan.py --gates
  python scripts/soft2d_ascan.py
  python scripts/soft2d_ascan.py --analyze-only

## ★ 왜 Zahn 의 좁은 창(`A = 10.03–10.75`)을 훑지 않았는가

`runs/2026-07-29_soft-r3-fss` 가 `η₆` 를 세 `A` 에서 쟀다: `A=0.1` → `2.03`,
`A=1` → `1.86`, `A=10` → `1.27`. `A=10` 은 Zahn hexatic 경계(`Γ=55.87`, `A=10.03`)의
**`−0.3 %`** 인데 `η₆` 가 상한 `0.25` 의 **5.1배**다. 경계 바로 아래라면 `0.25`
근처여야 한다 — 긴장이다. 그리고 `A=100` 은 이미 결정이다 (카드 §8.1).
⇒ **전이는 `A = 10–100` 사이이고 좁은 창이 아니다.**

## 설계 — `A` 를 로그로 훑는다

`A = 10, 13.3, 17.8, 23.7, 31.6` (사분-decade). 각 `A` 마다 `N` 사다리
`64·144·256·400` 으로 지수를 얻는다. `r_cut = 3.80` 고정 · 시드 `32–35` — S17 규약.
`prod_tau = 10 τ_d`, 창 `[6, 10]` (`τ_relax ≈ 0.098 τ_d` 의 60배).
`A = 10` 이 S17(`30 τ_d`)과 겹쳐 **런 길이 단축의 영향을 검사한다.**

## ★★ 스캔 상한은 예산이 아니라 절단오차가 정한다

`βU(r_cut) = A/r_cut³` 이므로 `A` 에 **비례**한다. `r_cut = 3.80` 고정에서
`A=100` 은 `1.82 kT` — 퍼텐셜을 아직 열에너지 규모인 곳에서 자른다.
`r_cut` 은 최소 상자(`N=64`, `L/2 = 4.0`)가 상한이라 키울 수 없다.
⇒ `A ≤ 31.6` (`βU ≤ 0.58`) 까지만. 그 위는 **더 큰 상자가 필요하다.**

## ★★★ `η₆ ≤ 1/4` 를 hexatic 이라고 읽으면 안 된다

**결정도 `η₆ ≈ 0` 을 만족한다** (`ψ₆ → const` 이므로 `p = 0`).
첫 판본이 이 구별 없이 `A ≥ 13.3` 의 **결정을 "hexatic 가능" 으로 출력했다.**
갈리는 것은 `ψ₆` 의 **크기**다 — 결정은 `O(1)`, hexatic 은 작고 느리게 감소.
⇒ `analysis.structure.phase_from_finite_size` 가 두 축으로 읽는다.

## 결과 요약

    A       Γ    ψ₆(400)        η₆      χ²/dof   상
   10    55.68     0.144   +1.40±0.22     0.04   등방 액체
   13.3  74.06     0.688   −0.19±0.04    46.9    결정
   17.8  99.12     0.778   −0.22±0.01   106      결정
   23.7 131.97     0.833   −0.41±0.01    24.8    결정
   31.6 175.96     0.872   −0.25±0.01   322      결정

액체→결정이 `A = 10–13.3` (33 % 폭) 안에 있다. Zahn 결정 경계 `Γ = 59.88`
(`A = 10.75`) 가 **그 브래킷 안**이다. hexatic 창(7 % 폭)은 브래킷 안에 숨어
있어 **이 스캔으로는 분해되지 않는다.**

`A ≥ 13.3` 의 `χ²` 가 거대하고 `p` 가 음수인 것은 결정에서 멱함수가 성립하지
않기 때문이다 (`ψ₆` 가 `1` 에 포화). `N=64` 는 `L = 8 d` 로 결정이 제대로
형성되지 않아 유독 낮다 (`0.394` vs `N=144` 의 `0.672`) — 그것이 `p` 를 음수로
만든다. **결정 영역에서 이 지수는 물리량이 아니다.**
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

#  ★★ 스캔 상한은 **절단오차**가 정한다, 예산이 아니다.
#  `βU(r_cut) = A/r_cut³` 이므로 `A` 에 **비례**한다. `r_cut = 3.80` 고정에서:
#      A=10 → 0.18 kT · A=31.6 → 0.58 · A=56.2 → 1.02 · A=100 → 1.82
#  `A=100` 은 퍼텐셜을 아직 `1.8 kT` 인 곳에서 자른다 — 구조가 틀린다.
#  `r_cut` 을 키우려면 최소 상자를 키워야 하는데 (`r_cut+buffer ≤ L/2`), 그러면
#  사다리의 아래쪽 점을 버려야 하고 `N` 레버암이 사라진다.
#  ⇒ **`A ≤ 31.6` 까지만 훑는다** (사분-decade 5점). 그 위는 더 큰 상자가 필요하다.
#    (예산 게이트도 A=100·N=400 을 665 s > 600 s 로 기각했다 — 두 제약이 같은 방향)
AMPLITUDES = (10.0, 13.3, 17.8, 23.7, 31.6)       # 사분-decade, 10 → 10^1.5
N_LADDER = FSS.N_LADDER                            # 64·144·256·400
SEEDS = FSS.SEEDS                                  # 32–35 (전 N 성공 선별)
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
    print(f"## 게이트 — r_cut = {R_CUT_FIXED} 고정 · prod {PROD_TAU:g} τ_d / "
          f"창 {WINDOW}\n")
    print(f"{'A':>7} {'Γ':>8} {'βU(rc)':>8} {'N=400 dt*':>10} "
          f"{'steps':>9} {'N=400 예상':>11}")
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
    print(f"\n  런 {n}개 (A {len(AMPLITUDES)} × N {len(N_LADDER)} × 시드 "
          f"{len(SEEDS)})")
    print(f"  최장 런 {worst:.0f} s {'≤' if worst <= budget else '>'} 예산 "
          f"{budget:.0f} s · 총 작업 {total/60:.0f} 분 → 동시 {k} 에서 "
          f"약 {total/k/60:.0f} 분")
    if worst > budget:
        raise SystemExit("⛔ 예산 초과 — 실행하지 않고 보고한다")
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
            if i % 20 == 0 or i == len(jobs):
                print(f"  {i}/{len(jobs)} … ({time.perf_counter()-t0:.0f} s)")
    wall = time.perf_counter() - t0
    print(f"  배치 wall {wall:.1f} s · 실패 {len(failed)}")
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
        #  ★ 지수만으로 판정하지 않는다 — 결정도 η₆ ≈ 0 을 만족한다.
        #    ψ₆ 의 크기를 함께 봐야 한다 (2026-07-29 에 이 구별 없이 결정을
        #    "hexatic 가능" 으로 출력했다).
        ph = phase_from_finite_size(fit, psi[-1])
        out[f"A{A:g}"]["phase"] = ph
        print(f"  A={A:<6g} Γ={zahn_phase(A)['gamma']:>7.2f}  "
              f"ψ₆(400)={psi[-1]:.4f}  p={fit.p:+.3f}±{fit.p_se:.3f}  "
              f"η₆={fit.eta6:+.3f}±{fit.eta6_se:.3f}  "
              f"χ²/dof={fit.chi2_reduced:>7.2f}  **{ph['phase']}**")
        if ph["exponent_alone_would_say"] == "hexatic-or-below" \
                and ph["phase"] == "crystal":
            print(f"           ⚠ 지수만 보면 'hexatic 이하' 다 — ψ₆ 크기가 "
                  f"결정으로 갈랐다")
    return out


def transition_bracket(metrics: dict) -> dict:
    """상이 바뀌는 `A` **구간**을 찾는다.

    ★ 이전 판본은 `η₆` 가 `0.25` 를 지나는 곳을 찾아 "hexatic 경계" 라고 불렀다.
      **틀렸다** — 결정도 `η₆ ≈ 0` 이라 그 교차점은 액체→결정 경계였다.
      상은 `phase_from_finite_size` 가 (지수, `ψ₆` 크기) 두 축으로 읽는다.
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
        out["note"] = f"스캔 전체가 {ph[0]} — 전이가 스캔 밖이다"
    else:
        out["note"] = " · ".join(
            f"{f} → {t} 가 A = {a:g}–{b:g} 사이 ({100*(b-a)/a:.0f} % 폭)"
            for a, b, f, t in changes)
    #  hexatic 이 관측됐는가
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
        print(f"\n  🔒 봉인 {seal.name} — "
              f"{len(seal.read_text().splitlines())}개 문서 (실행 전)")
        rd.write_json("manifest", run_batch(rd, rows, policy))

    print(f"\n## S7 — η₆(A) 스캔 (사다리 {list(N_LADDER)}, 창 {WINDOW})")
    metrics = analyze(rd)
    cr = transition_bracket(metrics)
    metrics["_transition_bracket"] = cr
    metrics["_provenance_at_analysis"] = provenance(DRIVERS)
    rd.write_json("metrics", metrics)
    print(f"\n## 상 전이 브래킷")
    print(f"  {cr['note']}")
    for t in cr["transitions"]:
        g_lo = zahn_phase(t["A_lo"])["gamma"]
        g_hi = zahn_phase(t["A_hi"])["gamma"]
        print(f"  → A = {t['A_lo']:g}–{t['A_hi']:g}  (Γ = {g_lo:.1f}–{g_hi:.1f})")
        print(f"    Zahn 결정 경계 Γ = 59.88 (A = 10.75) → "
              f"{'브래킷 **안**' if g_lo <= 59.88 <= g_hi else '브래킷 밖'}")
    print(f"  hexatic 관측? {'예' if cr['hexatic_observed'] else '**아니다**'} — "
          f"Zahn 창 A = 10.03–10.75 (7 % 폭) 은 위 브래킷 안에 숨어 있다")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
