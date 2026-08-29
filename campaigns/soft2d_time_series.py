"""soft-r3 `A` 스윕 — **시간분해**. `g(r)`·Voronoi·결함구조가 언제 만들어지는가.

usage:
  python scripts/soft2d_time_series.py --gates        # 게이트·비용만 (실행 안 함)
  python scripts/soft2d_time_series.py                # S5→S7 전부
  python scripts/soft2d_time_series.py --analyze-only # 이미 있는 raw/ 로 재분석

## 선행 런과 무엇이 다른가

`runs/2026-07-29_soft-r3-2d-A-sweep` 은 **평형화 구간을 버리고** 프로덕션 후반
절반만 분석했다. 그것은 "무엇이 되었는가" 에 답하고 **"언제 되었는가" 에는 답하지
않는다.** 이 런은 `equil_tau = 0` 으로 `t = 0` 부터 전 구간을 표집한다.

## 이 스크립트가 조심하는 것

**① 게이트 임계값을 다시 쓰지 않는다.** `simbot.nondim.dt_max_*` 를 축약 단위
   (`d = 1`, `D₀ = 1`, `γ = 1`) 로 호출한다. 문턱은 `config/run_policy.yaml` 에서
   읽는다 — 스크립트에 박아 두면 정책이 바뀌어도 따라오지 않는다 (2026-07-28 전례).

**② `max|F|` 를 실제로 계산한다.** 추정 금지 (master_plan §5.4). 무작위 초기배치에서
   재고, 프로덕션 중 최댓값도 함께 기록한다 — 대칭 배치에서는 `max|F| = 0` 이 되어
   게이트가 무력해지기 때문이다 (카드 §7).

**③ 봉인을 실행 **전에** 한다.** 실행 후 봉인은 아무것도 보증하지 않는다.

**④ 시드를 선행 런과 갈랐다** (`1–4` → `5–8`). HOOMD `Brownian` 은 같은 시드에서
   비트 단위로 재현되므로, 같은 시드를 쓰면 후반 창 비교가 **산술 항등식**이 된다.
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

# --- 사용자 지시 (S1 델타 D1–D4) --------------------------------------------
RUN_ID = "2026-07-29_soft-r3-time-resolved"
SRC = REPO / "examples" / "soft-r3-time-resolved"
AMPLITUDES = (0.1, 1.0, 10.0)          # D1 — A=100 제외
BOX_SHAPE = "square"                   # D2 — L_x = L_y 고정
INIT = "random"                        # D2 — 정사각 상자에서 hex 는 인공 결함
SIGMA_SI = 5.0e-6                      # D3
D_OVER_SIGMA = 3.0                     # D3 — coverage 8.73 % < 10 %
COVERAGE_MAX = 0.10
T_SI = 298.15
N_PARTICLES = 100
SEEDS = (5, 6, 7, 8)                   # ★ 선행 런(1–4)과 갈랐다
TOTAL_TAU = 80.0                       # D4 — 전 구간 표집 (선행 런과 같은 총량)
N_FRAMES = 400                         # stride = 0.2 tau_d
R_MIN = 0.05                           # Table 하한
MIN_SEP_INIT = 0.8
NLIST_BUFFER = 0.1
N_RDF_WINDOWS = 4
VORONOI_TIMES = (0.0, 1.0, 10.0, 40.0, 80.0)   # [tau_d] 타임랩스 시점 (근사)

#  ★ 초기 과도구간 전용 패스 (2026-07-29 추가)
#  본 패스의 stride = 80/400 = 0.2 τ_d 다. 그런데 초기배치의 배제부피 껍질
#  (`min_sep = 0.8 d`) 이 메워지는 시간은 자유확산으로 `t* ≈ 0.3²/4 = 0.023 τ_d` —
#  **첫 프레임보다 9배 빠르다.** 즉 본 패스는 과도구간을 전부 첫 프레임 안에
#  삼켜 버리고 "처음부터 정상상태" 로 보이게 만든다.
#  ⇒ 같은 시드로 짧고 촘촘한 패스를 따로 돌린다. 비용은 본 패스의 1/40 이다.
EARLY_TAU = 2.0
EARLY_FRAMES = 400                             # stride = 0.005 tau_d
#  ★ 시드 16개. 12런 3.9 s 였으므로 48런도 ~16 s 다 — **오차막대가 공짜인 구간**
#    (CLAUDE.md §자원 정책). 완화시간의 A 의존성이 이 런의 유일한 새 주장이므로
#    거기에 시드를 쓴다. 본 패스는 A=10 이 65 s/런이라 같은 선택을 할 수 없다.
EARLY_SEEDS = tuple(range(5, 21))


# =============================================================================
# S4 — dt 게이트 (축약 단위: d = 1, D0 = 1, gamma = 1, kT = 1)
# =============================================================================
def gate_table(policy) -> tuple[dict, list[dict]]:
    """`A` 마다 `dt*` 와 지배 게이트. **`max|F|` 를 실제로 계산한다.**"""
    ts = policy.timestep
    delta_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    delta_F = float(ts.get("max_force_displacement_sigma", 0.005))
    hard_floor = float(ts.get("hard_floor", 1e-7))

    dt_th = dt_max_thermal(delta_th, sigma=1.0, D0=1.0)      # 축약: d=1, D0=1
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
    print(f"## 기하 (사용자 지시 D3)\n"
          f"  σ = {geo['sigma_si']*1e6:.2f} µm · d/σ = {geo['d_over_sigma']:.3f} → "
          f"d = {geo['d_si']*1e6:.3f} µm · L = {geo['L_si']*1e6:.2f} µm "
          f"(L* = {geo['L_star']:.1f})\n"
          f"  **커버리지 = {geo['coverage']:.4%}** (상한 {COVERAGE_MAX:.1%}) · "
          f"σ/d = {geo['sigma_over_d']:.6g}\n"
          f"  η = {geo['eta_si']*1e3:.4f} mPa·s · D₀ = {geo['D0_si']*1e12:.5f} µm²/s · "
          f"**τ_d = {geo['tau_d_si']:.1f} s = {geo['tau_d_si']/60:.2f} 분**\n"
          f"  → {TOTAL_TAU:g} τ_d = {TOTAL_TAU*geo['tau_d_si']/3600:.1f} 시간의 "
          f"실제 실험 시간\n")
    print(f"## dt 게이트 (변위 기준, 축약 단위)  문턱: 열 {thresholds['delta_thermal']} d · "
          f"힘 {thresholds['delta_force']} d\n")
    hdr = (f"{'A':>6} {'Γ':>7} {'Zahn':<16} {'max|F*|':>9} {'dt_th':>10} "
           f"{'dt_F':>10} {'지배':<22} {'dt*':>10} {'steps':>10} {'βU(rcut)':>9}")
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
    print(f"\n## 비용  런 {n_runs}개 = A {len(rows)}개 × 시드 {len(SEEDS)}개 · "
          f"동시 {k} (효율 {eff:.2f})")
    print(f"  Λ = {lam:.3g} 입자·스텝/s → 적분 예상 {per_seed:.0f} s "
          f"(= {per_seed/60:.1f} 분, 동시실행 감안 전)")
    print(f"  ⚠ 예상은 적분만이다. 프레임 {N_FRAMES}개의 스냅샷 추출 오버헤드가 "
          f"선행 런에서 컸다 (A=10: 적분 20 s vs 실측 67 s)")
    budget = policy.wall_budget_s
    worst = max(N_PARTICLES * r["steps"] / lam for r in rows) * 3.4   # 오버헤드 계수
    if worst > budget:
        print(f"\n  ⛔ 최장 런 예상 {worst:.0f} s > 예산 {budget:.0f} s/런 — "
              f"실행하지 않고 보고한다 (CLAUDE.md §자원 정책)")
    else:
        print(f"  ✅ 최장 런 예상 {worst:.0f} s ≤ 예산 {budget:.0f} s/런")
    return total


# =============================================================================
# S5 — 실행
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
                equil_tau=0.0,                       # ★ D4 — t=0 부터 표집
                prod_tau=prod_tau, n_frames=n_frames, seed=s, label=label)
            jobs.append((asdict(cfg), str(outroot / label)))

    k = policy.concurrency("default")
    tag = "초기 과도구간" if early else "본"
    print(f"\n## S5 — {tag} 패스 {len(jobs)} 런 시작 (동시 {k}, "
          f"prod {prod_tau:g} τ_d / {n_frames} 프레임 → "
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
    print(f"\n  배치 wall {wall:.1f} s")
    return {"done": done, "failed": failed, "batch_wall_s": wall,
            "concurrency": k, "n_jobs": len(jobs), "early": early,
            "prod_tau": prod_tau, "n_frames": n_frames}


# =============================================================================
# S7b — 초기 과도구간: 완화가 **있는가**, 있으면 얼마인가
# =============================================================================
def analyze_early(rd: RunDir, geo: dict) -> dict:
    """짧고 촘촘한 패스. 초기배치(`t = 0`)를 시계열 앞에 붙여서 본다."""
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
        #  ★ 초기배치를 t=0 프레임으로 앞에 붙인다 — 러너는 첫 stride 뒤부터 표집한다
        frames = np.concatenate([r["init_pos"][None].astype(np.float32), r["traj"]])
        t = np.concatenate([[0.0], r["t_star"]])
        ser = hex_order_series(frames, Lx=r["Lx"], Ly=r["Ly"], t_star=t,
                               coord_range=(3, 12))
        series_by_A.setdefault(A, []).append(ser)

    for A, sers in sorted(series_by_A.items()):
        t = sers[0].t_star
        mat = np.array([s.defect_fraction for s in sers])
        mean = mat.mean(axis=0)
        #  잡음 = 후반 절반의 프레임 간 표준편차 (시드 평균 곡선에서)
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
        print(f"  A={A:<5g} 결함 t=0 {mean[0]:.4f} → 후반 {tail.mean():.4f} "
              f"(프레임 SD {noise:.4f}) · τ = {fit.tau:.4g} τ_d "
              f"({'채택' if fit.converged else '거부'}) "
              f"{'· ' + fit.note if fit.note else ''}")
    return {"series_by_A": series_by_A, "fits": fits, "summary": summary}


# =============================================================================
# S7 — 시간분해 분석
# =============================================================================
def load_run(d: Path) -> dict:
    man = json.loads((d / "manifest.json").read_text())
    z = np.load(d / "samples.npz")
    cfg = man["config"]
    stride = int(z["stride"][0])
    n = z["traj"].shape[0]
    #  프레임 k 는 (k+1)*stride 스텝 뒤다 (러너가 먼저 run(stride) 한 뒤 표집한다)
    t_star = (np.arange(1, n + 1) * stride * cfg["dt_star"])
    return {"man": man, "cfg": cfg, "traj": z["traj"], "energy": z["energy"],
            "max_force": z["max_force"], "init_pos": z["init_pos"],
            "Lx": float(z["box"][0]), "Ly": float(z["box"][1]),
            "t_star": t_star, "stride": stride}


def analyze(rd: RunDir, geo: dict) -> dict:
    dirs = sorted(p for p in rd.raw.glob("A*_s*") if (p / "samples.npz").exists())
    if not dirs:
        raise SystemExit(f"⛔ {rd.raw} 에 런이 없다")

    by_A: dict[float, list[dict]] = {}
    for d in dirs:
        r = load_run(d)
        A = float(r["cfg"]["amplitude"])
        ser = hex_order_series(r["traj"], Lx=r["Lx"], Ly=r["Ly"],
                               t_star=r["t_star"], coord_range=(3, 12))
        by_A.setdefault(A, []).append({"run": r, "series": ser, "dir": d})
        print(f"  분석 {d.name:<12} 프레임 {ser.n_frames} · "
              f"ψ₆(마지막) {ser.psi6_global[-1]:.4f}")

    #  --- 후반 창: 선행 런의 후반 절반과 같은 시각으로 잡는다 ---
    #      선행 총 길이가 A 마다 달랐다 (A ≤ 1 은 30 τ_d, A = 10 은 80 τ_d)
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
        #  에너지: 같은 창에서 시드별 평균
        e_pp = []
        for x in runs:
            r = x["run"]
            m = (r["t_star"] >= lo) & (r["t_star"] < hi)
            e_pp.append(float(r["energy"][m].mean()) / N_PARTICLES)
        agg["energy_pp"] = aggregate_seeds(e_pp)
        energy_by_A[A] = (t, [x["run"]["energy"] / N_PARTICLES for x in runs])

        #  최종 창의 배위수 종류·5-7 불균형 (시드 평균)
        kinds = [float(s.coord_kinds[(s.t_star >= lo) & (s.t_star < hi)].mean())
                 for s in sers]
        #  ★ 선행 런과 **같은 방식**의 배위수 종류 수: 프레임·시드를 전부 집계한
        #    히스토그램에서 분율 > 0.5 % 인 종류를 센다. 프레임별로 세면 다른 수가
        #    나온다 — `N = 100` 에서 입자 1개가 이미 1 % 이므로 프레임 문턱이
        #    "존재하기만 하면 통과" 로 무력해진다. 두 추정량을 함께 보고한다.
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

        #  t90 — 후반 창 평균의 90 % 에 처음 도달하는 시각 (시드 평균 곡선에서)
        psi_mat = np.array([s.psi6_global for s in sers])
        psi_mean = psi_mat.mean(axis=0)
        target = 0.9 * agg["psi6_global"].mean
        reached = np.nonzero(psi_mean >= target)[0]
        t90 = float(t[reached[0]]) if reached.size else float("nan")
        #  ★ 액체는 t=0 부터 이미 목표 위에 있다 (무작위 배치가 이미 액체 배치다)
        #    → t90 ≈ 첫 프레임. 그 사실 자체가 결과다

        #  시간창 g(r) — 시드를 모아서 창마다 (프레임 수 = 4 시드 × 100 프레임)
        all_traj = np.concatenate([x["run"]["traj"] for x in runs], axis=0)
        all_t = np.concatenate([x["run"]["t_star"] for x in runs])
        order = np.argsort(all_t, kind="stable")
        windows_by_A[A] = rdf_windows(all_traj[order], Lx=runs[0]["run"]["Lx"],
                                      Ly=runs[0]["run"]["Ly"],
                                      t_star=all_t[order],
                                      n_windows=N_RDF_WINDOWS, bins=200)

        #  최소분리 (전 궤적) · σ 환산
        min_sep = min(x["run"]["man"]["guards"]["min_separation"] for x in runs)
        #  S(k) 후반 창 — 상자 모양이 같으므로 A 끼리 비교 가능
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

        #  --- P4: 초기조건 껍질(min_sep = 0.8 d)이 채워지는가 ---
        w = windows_by_A[A]
        i_half = int(np.argmin(np.abs(w.r - 0.5)))
        metrics[f"A{A:g}"]["g_at_0.5d_by_window"] = w.g[:, i_half].tolist()
        i_08 = int(np.argmin(np.abs(w.r - 0.75)))
        metrics[f"A{A:g}"]["g_at_0.75d_by_window"] = w.g[:, i_08].tolist()

    return {"metrics": metrics, "series_by_A": series_by_A,
            "windows_by_A": windows_by_A, "energy_by_A": energy_by_A,
            "by_A": by_A}


# =============================================================================
# S6 — 그림
# =============================================================================
def figures(rd: RunDir, res: dict, geo: dict, early: dict | None = None
            ) -> FigureSet:
    fs = FigureSet(rd.figs)
    plot_structure_timeseries(fs, res["series_by_A"], tau_d_si=geo["tau_d_si"],
                              energy_by_A=res["energy_by_A"])
    plot_rdf_evolution(fs, res["windows_by_A"], tau_d_si=geo["tau_d_si"],
                       sigma_over_d=geo["sigma_over_d"])

    for i, (A, runs) in enumerate(sorted(res["by_A"].items())):
        r = runs[0]["run"]                       # 시드 1개의 타임랩스 (배치가 아니다)
        t = r["t_star"]
        idx = [int(np.argmin(np.abs(t - tv))) for tv in VORONOI_TIMES if tv > 0]
        frames = [voronoi_frame(r["init_pos"], Lx=r["Lx"], Ly=r["Ly"])]
        times = [0.0]
        for j in idx:
            frames.append(voronoi_frame(r["traj"][j], Lx=r["Lx"], Ly=r["Ly"]))
            times.append(float(t[j]))
        #  ★ 프레임별 요동 폭을 캡션에 넣는다 — 패널 간 차이를 완화로 읽지 않게
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
                "초기 과도구간 패스(raw_early/)가 없다 — `--early` 로 돌릴 것")

    fs.skip("snapshots_3d", "2D 계다 — 3D 레이트레이싱(fresnel)은 설치하지 않았다")
    fs.skip("msd", "이 카드의 목적동역학은 평형 구조다. 수송량은 묻지 않았다 "
                   "(카드 §2 관측량 목록에 MSD 가 없다)")
    return fs


# =============================================================================
# S7 — 봉인 검증 + 예측 대조
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
            #  ★ 허용오차는 **선행 런의 SE 만으로** 세워졌다 (봉인 시점에 새 SE 를
            #    알 수 없었으므로). 실제 SE_diff 로 다시 재면 몇 σ 인지 함께 적는다 —
            #    "허용오차를 넘었다" 와 "유의하게 다르다" 는 다른 말이다.
            se_new = m[field]["se"]
            se_prior = float(str(item["tolerance"]).lstrip("±")) / (3.0 * 2 ** 0.5)
            se_diff = float(np.hypot(se_new, se_prior))
            sig = abs(measured - float(item["value"])) / se_diff if se_diff else None
            note = (f"SE_new = {se_new:.4g} (시드 {m['n_seeds']}개) · "
                    f"SE_prior = {se_prior:.4g} · **실제 {sig:.2f}σ**"
                    if sig is not None else f"SE = {se_new:.4g}")
        elif q.startswith("coord_kinds__"):
            akey = q.split("__")[1]
            measured = metrics[akey]["coord_kinds"]["mean"]
            note = f"시드별 {metrics[akey]['coord_kinds']['values']}"
        elif q.startswith("min_separation_over_sigma__"):
            akey = q.split("__")[1]
            measured = metrics[akey]["min_separation_over_sigma"]
            note = f"= {metrics[akey]['min_separation_d']:.4g} d"
        elif q == "rdf_at_0.5d__A0.1__first_window":
            measured = metrics["A0.1"]["g_at_0.5d_by_window"][0]
            note = "첫 시간창 g(0.5 d)"
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
    """PASS/FAIL/INCONCLUSIVE. **판정을 제안만 한다** — 확정은 사람이 한다."""
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
    if isinstance(measured, dict):                     # t90 순서 예측
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
    #  --early: 본 패스가 이미 있는 상태에서 초기 과도구간 패스만 추가한다
    early_only = "--early" in sys.argv

    if early_only:
        batch = run_batch(rd, rows, policy, early=True)
        #  실측 wall 은 배치 요약(`05_run_manifest.json`)에 둔다 — 리포트가 본 패스와
        #  비용을 대조할 때 마크다운을 되읽지 않게
        if rd.exists("manifest"):
            man = json.loads(rd.file("manifest").read_text())
            man["early_batch"] = batch
            rd.write_json("manifest", man)
        #  07b 는 마크다운 문서다 — `write_json` 을 쓰면 `.md` 안에 JSON 이 들어간다
        rd.write("sensitivity", "\n".join([
            "# S7b — 초기 과도구간 패스 (표집 간격 민감도)", "",
            f"본 패스의 stride(`{TOTAL_TAU/N_FRAMES:.3g} τ_d`) 가 완화시간보다 굵어서 "
            f"과도구간을 첫 프레임 안에 삼켰다. 같은 시드로 짧고 촘촘한 패스를 "
            f"따로 돌린 기록이다.", "",
            "| 항목 | 값 |", "|---|---|",
            f"| 런 수 | {batch['n_jobs']} (`A` {len(rows)}개 × 시드 "
            f"{len(EARLY_SEEDS)}개) |",
            f"| `prod_tau` | `{batch['prod_tau']:g} τ_d` (본 패스 "
            f"`{TOTAL_TAU:g}`) |",
            f"| 프레임 | `{batch['n_frames']}` → stride "
            f"`{batch['prod_tau']/batch['n_frames']:.4g} τ_d` |",
            f"| 동시 실행 | {batch['concurrency']} |",
            f"| 배치 wall | `{batch['batch_wall_s']:.1f} s` |",
            f"| 실패 | {len(batch['failed'])} |", "",
            "> **시드를 16개로 올린 이유** — 이 패스는 12런에 4초다. 완화시간의 `A` "
            "의존성이 이 런의 유일한 새 주장이므로 거기에 시드를 쓴다. 본 패스는 "
            "`A=10` 이 `65 s/런` 이라 같은 선택을 할 수 없다 "
            "(CLAUDE.md §\"오차 막대는 공짜다\").", "",
            "결과 표와 완화시간 적합: [`07_validation.md`](07_validation.md) §4b · "
            "그림 [`figs/06_early_transient.png`](figs/06_early_transient.png)", "",
            "근거: [[coarse-sampling-hides-the-whole-transient]]", "",
        ]))

    if not analyze_only and not early_only:
        # --- 봉인: 실행 **전에** ---
        shutil.copy2(SRC / "01_intake.md", rd.file("intake"))
        pred = yaml.safe_load((SRC / "prediction.yaml").read_text())
        rd.write("prediction", (SRC / "prediction.yaml").read_text())
        rd.write_json("prediction_json", pred)
        #  ★ 드라이버를 **이름이 아니라 해시로** 봉인한다. 이름만 적으면 나중에
        #    스크립트가 바뀌어도 spec 이 그대로여서 "무엇이 이 런을 만들었나" 에
        #    답할 수 없다 (2026-07-29 이 런에서 실제로 구멍이었다).
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
        print(f"\n  🔒 봉인 {seal.name} — {len(seal.read_text().splitlines())}개 문서")
        batch = run_batch(rd, rows, policy)
        rd.write_json("manifest", batch)

    pred = yaml.safe_load(rd.file("prediction").read_text())
    print("\n## S7 — 시간분해 분석")
    res = analyze(rd, geo)

    early = None
    if (rd.path / "raw_early").exists():
        print("\n## S7b — 초기 과도구간 (촘촘한 패스)")
        early = analyze_early(rd, geo) or None

    fs = figures(rd, res, geo, early)
    rd.write("figures", fs.figures_md())
    metrics = dict(res["metrics"])
    if early:
        metrics["_early_transient"] = early["summary"]
    #  ★ 분석 시점의 provenance. 궤적 manifest 의 해시는 *궤적 생성 시점*이고
    #    분석은 `--analyze-only` 로 따로 돌 수 있다 — 이 블록이 없으면
    #    metrics.json 이 어느 분석 코드·어느 freud 에서 나왔는지 알 수 없다.
    metrics["_provenance_at_analysis"] = provenance(__file__)
    rd.write_json("metrics", metrics)

    checks = check_predictions(pred, res["metrics"], geo)
    rd.write_json("prediction_json", {"items": pred["items"], "checks": checks})

    print("\n## 예측 대조 (판정은 제안이다 — confirmed_by: null)")
    hdr = f"{'항목':<44} {'예측':>14} {'측정':>14} {'허용':>12} 판정"
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
          f"미평가 {len(checks) - n_pass - n_fail}")
    print(f"\n→ {rd.path.relative_to(REPO)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
