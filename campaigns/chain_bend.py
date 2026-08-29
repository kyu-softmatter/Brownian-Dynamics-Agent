"""콜로이드 사슬 굽힘 강성 — 사슬 길이 N 스윕.

목표: `kappa(N) ~ L^-3` 의 **지수와 계수**를 확인한다.
Pantina & Furst (PRL 94, 138301) 의 `kappa = kappa_0 (a/s)^(2+d_b)`, `d_b=1` 의 BD 대응물.
**지수 -3 은 규약 무관, 계수는 규약 의존** — 논문 PDF 의 계수는 2단 조판 추출에서
깨져 신뢰할 수 없으므로(env-log 3단계) 교과서 보 이론을 기준으로 검증한다.

축약 단위: 길이 sigma = 2a = 결합 길이 = 1, 에너지 kT = 1, 시간 tau_D = sigma^2/D0
(gamma* = 1, D0* = 1).

## 예측 (단순지지 보 + 이산 사슬 매핑)

    이산 굽힘 에너지  U = 1/2 k_th sum (dtheta_i)^2
    연속 대응          EI = k_th * b,   b = 결합 길이 = 1
    중앙 점하중, 양단 단순지지:  delta = F L^3 / (48 EI)

    ->  kappa* = 48 k_theta* / (N-1)^3

양 끝 입자는 위치가 고정되지만 **각 구속은 없다**(각 삼중항의 꼭짓점은 1..N-2 뿐)
-> 단순지지(pinned)이고 clamped 가 아니다. 그래서 48 이다.

## 두 프로토콜을 함께 돌린다

- `mode="static"`  : **결정론.** kT=0, 중앙 입자에 알려진 힘 F* 를 걸고 이완 -> kappa = F/delta.
  노이즈가 없으므로 매핑(계수 48)을 날카롭게 판정한다. 통계 불필요.
- `mode="thermal"` : 요동. 중앙 입자 처짐의 등분배 -> kappa* = 1/<y_c*^2>.
  오차막대가 나오고, 나중에 G'(omega) 로 확장되는 경로다.

둘이 일치하면 매핑이 맞다. 어긋나면 어느 쪽이 틀렸는지 분리된다.

⚠ 거리 구속을 쓰지 않는다 — `Brownian` 과 호환 안 되고 조용히 발산한다.
  findings/dead-end-distance-constraint-with-brownian.md
"""
from __future__ import annotations

import json
import platform
import sys
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from simbot.build import chain_snapshot
from simbot.guards import (assert_statistic_fluctuates, check_bond_lengths,
                           check_finite)
from simbot.io import code_hash as _code_hash, git_dirty as _git_dirty, git_rev as _git_rev
from simbot.nondim import dt_max_force, dt_max_stability, dt_max_thermal
from simbot.policy import load_policy

LAMBDA = 6.3e6            # 처리량 상수 [입자*스텝/s], config/run_policy.yaml
WALL_BUDGET_S = 600.0     # CLAUDE.md: 10분/런
MAX_FRAMES = 8000         # get_snapshot 호출 비용 상한

# ★ dt 게이트의 식은 `simbot.nondim`, 문턱은 `config/run_policy.yaml` 이 소유한다.
#   이 스크립트가 `0.03`·`0.005`·`0.2` 를 다시 쓰면 정책을 고쳐도 따라오지 않는다.
#   축약 단위가 sigma = kT = gamma = D0 = 1 이므로 게이트 함수에 1.0 을 넘기면
#   그대로 축약 상한이 나온다 (게이트 식은 단위 무관).
_TS = load_policy().timestep


# =============================================================================
@dataclass
class ChainConfig:
    n_particles: int = 11             # 홀수만 (중앙 입자가 정의되어야 한다)
    dim: int = 2
    mode: str = "thermal"             # "static" | "thermal"
    k_theta_star: float = 1.0e4       # 각 스프링 [kT/rad^2]
    k_bond_ratio: float = 100.0       # k_bond* = ratio * kappa(N)*.  findings 참조:
                                      # 상수로 두면 스텝 수가 L^3 로 늘어난다
    delta_over_span: float = 0.01     # static: 목표 처짐 delta/L (선형 유지)
    dt_star: float = 0.0              # 0 이면 측정된 힘에서 게이트로 정한다
    relax_tau_bend: float = 30.0      # static: 이완 길이 [tau_bend]
    equil_tau_bend: float = 30.0      # thermal: 평형화 길이
    prod_tau_bend: float = 400.0      # thermal: 프로덕션 길이
    frames_per_tau_bend: float = 5.0  # thermal: tau_bend 당 표집 프레임
    seed: int = 1
    label: str = ""

    def __post_init__(self):
        if self.n_particles % 2 == 0:
            raise ValueError(f"n={self.n_particles} — 중앙 입자가 있어야 하므로 홀수만")
        if self.n_particles < 5:
            raise ValueError(f"n={self.n_particles} — 내부 자유도가 최소 3개 필요")
        if self.mode not in ("static", "thermal"):
            raise ValueError(f"mode={self.mode!r}")

    @property
    def span(self) -> float:
        """지점 간 거리 L* = (N-1) * 결합길이."""
        return float(self.n_particles - 1)

    @property
    def kappa_pred_star(self) -> float:
        return 48.0 * self.k_theta_star / self.span ** 3

    @property
    def tau_bend_star(self) -> float:
        """굽힘 완화시간 추정 gamma*/kappa*  (gamma* = 1)."""
        return 1.0 / self.kappa_pred_star

    @property
    def k_bond_star(self) -> float:
        """결합 스프링을 kappa(N) 에 비례시킨다 -> dt 가 L^3 로 커져 비용이 N 에 거의 무관."""
        return self.k_bond_ratio * self.kappa_pred_star

    @property
    def lambda_max_star(self) -> float:
        """강성행렬 최대고유값 추정. 1D 스프링 사슬 4k + 굽힘(4계 차분) 16k."""
        return 4.0 * self.k_bond_star + 16.0 * self.k_theta_star

    @property
    def dt_stability_star(self) -> float:
        """명시적 오일러 안정 한계 `2/lambda_max` 에 정책의 안전계수.

        실측 임계값은 이 하한의 1.22-2.80 배 -> 안전계수 0.2 에서 6-14 배 여유.
        근거: findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
        ★ 변위 게이트만으로는 못 잡는다 — 곧은 사슬은 max|F*| = 0 이라 힘 게이트가 무력하다.
        """
        return dt_max_stability(_TS["stability_safety_factor"], 1.0,
                                self.lambda_max_star)

    @property
    def load_star(self) -> float:
        """static 모드에서 걸 힘. delta = F/kappa 이므로 F = kappa * delta."""
        return self.kappa_pred_star * self.delta_over_span * self.span


# =============================================================================
def _seed_dt(cfg: ChainConfig) -> float:
    """dt 초기 추정. **안정성 항이 거의 항상 binding 이다.**

    직선 사슬에서 힘을 재면 정확히 0 이라 힘 게이트가 무력화된다
    (2026-07-28 실측: dt 가 4.5e-4 로 잡혀 kT=0 에서도 사슬이 터졌다).
    근거: findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
    """
    f = [np.sqrt(cfg.k_bond_star), np.sqrt(cfg.k_theta_star), 1.0]
    if cfg.mode == "static":
        f.append(cfg.load_star)
    proxy_force = 3.0 * max(f)          # 실측 전 임시 프록시. 예열 후 실측값으로 교체된다
    return min(_dt_gates(proxy_force, cfg).values())


def _dt_gates(max_f_star: float, cfg: ChainConfig) -> dict[str, float]:
    """게이트별 dt 상한 (축약 단위). 무력한 게이트는 **키가 빠진다.**"""
    gates = {
        "diffusion": dt_max_thermal(_TS["max_thermal_displacement_sigma"], 1.0, 1.0),
        "force": dt_max_force(_TS["max_force_displacement_sigma"], 1.0, 1.0, max_f_star),
        "stability": cfg.dt_stability_star,
    }
    return {k: v for k, v in gates.items() if v is not None}


def _dt_from_force(max_f_star: float, cfg: ChainConfig) -> tuple[float, dict]:
    """dt 게이트 셋 — 변위 2개(정확도) + 안정성 1개.

    변위 게이트는 여전히 필수지만 **불충분**하다. `dt/tau_D` 고정 게이트는 쓰지 않는다.
    식·문턱은 `simbot.nondim` + `run_policy.yaml` 이 소유한다 (복제 금지).
    """
    gates = _dt_gates(max_f_star, cfg)
    binding = min(gates, key=gates.get)
    return gates[binding], {
        "dt_diffusion_gate": gates["diffusion"],
        # 곧은 사슬에서는 max|F*| = 0 이라 힘 게이트가 **존재하지 않는다** (None).
        # 1e300 같은 수로 뭉개면 "게이트가 있었지만 느슨했다"로 잘못 읽힌다.
        "dt_force_gate": gates.get("force"),
        "dt_stability_gate": gates["stability"],
        "lambda_max_star": cfg.lambda_max_star,
        "k_bond_star": cfg.k_bond_star,
        "max_force_star_measured": float(max_f_star),
        "binding": binding}


def _integrated_act(x: np.ndarray, c_window: float = 6.0) -> float:
    """적분 자기상관시간 [프레임]. Sokal 자동창. 굽힘 모드는 느리므로 독립 가정 금지."""
    x = np.asarray(x, dtype=np.float64)
    x = x - x.mean()
    n = x.size
    f = np.fft.rfft(x, 2 * n)
    ac = np.fft.irfft(f * np.conj(f))[:n].real
    if ac[0] <= 0:
        return 1.0
    ac /= ac[0]
    tau = 0.5
    for k in range(1, n):
        tau += ac[k]
        if k >= c_window * tau:
            break
    return max(float(tau), 0.5)


def _build(hoomd, cfg: ChainConfig):
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=cfg.seed)
    snap = chain_snapshot(hoomd, n=cfg.n_particles, dim=cfg.dim)
    sim.create_state_from_snapshot(snap)

    bond = hoomd.md.bond.Harmonic()
    bond.params["b"] = dict(k=cfg.k_bond_star, r0=1.0)
    angle = hoomd.md.angle.Harmonic()
    angle.params["a"] = dict(k=cfg.k_theta_star, t0=np.pi)
    forces = [bond, angle]

    if cfg.mode == "static":
        load = hoomd.md.force.Constant(filter=hoomd.filter.Type(["C"]))
        load.constant_force["C"] = (0.0, cfg.load_star, 0.0)
        load.constant_torque["C"] = (0.0, 0.0, 0.0)
        forces.append(load)

    # 양 끝(타입 "E")은 적분하지 않는다 -> 완전 고정 지점.
    # 트랩으로 잡으면 트랩 컴플라이언스가 kappa 에 섞인다 (논문의 40 pN/um 천장 문제).
    kT = 0.0 if cfg.mode == "static" else 1.0
    bd = hoomd.md.methods.Brownian(filter=hoomd.filter.Type(["A", "C"]),
                                   kT=kT, default_gamma=1.0)
    sim.operations.integrator = hoomd.md.Integrator(
        dt=_seed_dt(cfg), methods=[bd], forces=forces)
    return sim, np.array(snap.bonds.group)


def _max_force(sim) -> float:
    return max(float(np.abs(np.array(f.forces)).max()) if np.size(f.forces) else 0.0
               for f in sim.operations.integrator.forces)


def run_chain(cfg: ChainConfig) -> dict:
    import hoomd

    t0 = time.perf_counter()
    sim, bonds = _build(hoomd, cfg)
    ic = (cfg.n_particles - 1) // 2
    tau_b = cfg.tau_bend_star

    # --- 1) 예열 후 실측 힘으로 dt 확정 ---
    dt_seed = sim.operations.integrator.dt
    warm = max(50, int(round(0.5 * tau_b / dt_seed)))
    sim.run(min(warm, 20000))
    max_f = _max_force(sim)
    # 여러 시점의 최댓값 (한 순간만 보면 과소평가한다)
    for _ in range(20):
        sim.run(max(1, warm // 20))
        max_f = max(max_f, _max_force(sim))
    dt, dt_info = _dt_from_force(max_f, cfg)
    dt_info["dt_seed"] = dt_seed
    if cfg.dt_star > 0:
        dt = cfg.dt_star
        dt_info["overridden_to"] = dt
    sim.operations.integrator.dt = dt

    def positions() -> np.ndarray:
        return np.array(sim.state.get_snapshot().particles.position,
                        dtype=np.float64)

    guard_fail: list[str] = []
    bond_info: dict = {}
    max_disp = 0.0
    out: dict = {}

    # =====================================================================
    if cfg.mode == "static":
        # **적응형 이완.** 고정 길이(30 tau_bend)로는 안 된다 — 실측 결과 참 완화시간은
        # 1/kappa_point 가 아니라 그보다 N 배쯤 길다 (사슬 전체가 끌려야 하므로).
        # 2026-07-28: 고정 길이로 돌려 N=41 에서 kappa 가 22 % 높게 나왔다.
        chunk = max(1, int(round(0.5 * tau_b / dt)))
        max_steps = int(0.8 * WALL_BUDGET_S * LAMBDA / cfg.n_particles)
        tol, need = 1e-5, 3
        curve, tsteps = [], []
        done, prev, stable = 0, 0.0, 0
        while done < max_steps and len(curve) < 4000:
            sim.run(chunk)
            done += chunk
            p = positions()
            y = float(p[ic, 1])
            curve.append(y)
            tsteps.append(done)
            ok, fails = check_finite(position=p)
            if not ok:
                guard_fail += fails
                break
            if y > 0 and abs(y - prev) / y < tol:
                stable += 1
                if stable >= need:
                    break
            else:
                stable = 0
            prev = y

        p = positions()
        ok_b, bond_info = check_bond_lengths(p, bonds, target=1.0, tol=0.05,
                                            dims=cfg.dim)
        if not ok_b:
            guard_fail.append(f"결합 길이 위반 max_rel_dev={bond_info['max_rel_dev']:.3e}")

        delta = float(p[ic, 1])
        c = np.asarray(curve, dtype=np.float64)
        # 마지막 구간의 잔여 변화율
        drift = float(abs(c[-1] - c[-2]) / max(abs(c[-1]), 1e-300)) if c.size > 1 else np.nan
        # Aitken(등간격 3점) 외삽 -> delta_inf. 이완이 덜 끝났어도 고원값을 추정한다.
        delta_inf, tau_fit = float("nan"), float("nan")
        if c.size >= 9:
            i1, i2, i3 = c.size - 1 - 2 * (c.size // 8), c.size - 1 - (c.size // 8), c.size - 1
            d1, d2, d3 = c[i1], c[i2], c[i3]
            den = d1 + d3 - 2 * d2
            if abs(den) > 1e-300:
                delta_inf = float((d1 * d3 - d2 * d2) / den)
            # delta(t) = d_inf (1 - A exp(-t/tau)) 를 로그선형 적합
            if np.isfinite(delta_inf) and delta_inf > c[-1]:
                res = delta_inf - c
                m = res > 0
                if m.sum() >= 4:
                    sl = np.polyfit(np.asarray(tsteps)[m] * dt, np.log(res[m]), 1)[0]
                    if sl < 0:
                        tau_fit = float(-1.0 / sl)

        kappa = cfg.load_star / delta if delta != 0 else float("nan")
        kappa_inf = (cfg.load_star / delta_inf
                     if np.isfinite(delta_inf) and delta_inf != 0 else float("nan"))
        out = {
            "delta_star": delta, "delta_over_span": delta / cfg.span,
            "load_star": cfg.load_star,
            "kappa_star": kappa, "kappa_star_se": 0.0,   # 결정론 — 통계오차 없음
            "delta_inf_star": delta_inf,
            "kappa_inf_star": kappa_inf,
            "ratio_inf_over_pred": kappa_inf / cfg.kappa_pred_star,
            "tau_relax_fit_star": tau_fit,
            "tau_relax_over_tau_bend": tau_fit / tau_b if np.isfinite(tau_fit) else float("nan"),
            "relax_curve": curve, "relax_tsteps": tsteps,
            "tail_drift_rel": drift,
            "converged": bool(stable >= need),
            "steps": {"relax": done, "chunk": chunk, "n_chunks": len(curve),
                      "max_steps": max_steps},
        }

    # =====================================================================
    else:
        equil_steps = max(200, int(round(cfg.equil_tau_bend * tau_b / dt)))
        prod_steps = max(1000, int(round(cfg.prod_tau_bend * tau_b / dt)))
        stride = max(1, int(round(tau_b / (cfg.frames_per_tau_bend * dt))))
        n_frames = min(MAX_FRAMES, prod_steps // stride)
        prod_steps = n_frames * stride

        est = cfg.n_particles * (equil_steps + prod_steps) / LAMBDA
        if est > WALL_BUDGET_S:
            return {"config": asdict(cfg), "skipped": True, "est_wall_s": est,
                    "reason": f"추정 wall {est:.0f}s > 예산 {WALL_BUDGET_S:.0f}s",
                    "dt_star": dt, "dt_info": dt_info,
                    "steps": {"equil": equil_steps, "prod": prod_steps}}

        sim.run(equil_steps)
        yc = np.empty(n_frames, dtype=np.float64)
        prev = positions()
        nf = 0
        for f in range(n_frames):
            sim.run(stride)
            p = positions()
            yc[f] = p[ic, 1]
            nf = f + 1
            ok, fails = check_finite(position=p)
            if not ok:
                guard_fail += [f"frame {f}: {x}" for x in fails]
                break
            max_disp = max(max_disp, float(np.abs(p - prev).max()))
            prev = p
            if f % max(1, n_frames // 20) == 0:
                ok_b, bond_info = check_bond_lengths(p, bonds, target=1.0,
                                                    tol=0.05, dims=cfg.dim)
                if not ok_b:
                    guard_fail.append(
                        f"frame {f}: 결합 길이 위반 "
                        f"max_rel_dev={bond_info['max_rel_dev']:.3e}")
                    break
        yc = yc[:nf]

        var = float(np.mean(yc ** 2))
        act = _integrated_act(yc)
        n_eff = max(2.0, yc.size / (2.0 * act))
        var_rel_se = np.sqrt(2.0 / n_eff)
        kappa = 1.0 / var if var > 0 else float("nan")

        fluct_ok, fluct_msg = True, ""
        try:
            blocks = [float(np.mean(b ** 2))
                      for b in np.array_split(yc, 20) if b.size]
            assert_statistic_fluctuates(blocks, name="<y_c^2> 블록평균")
        except (AssertionError, ValueError) as e:
            fluct_ok, fluct_msg = False, str(e)

        out = {
            "kappa_star": kappa, "kappa_star_se": kappa * var_rel_se,
            "var_yc_star": var, "rms_yc_star": float(np.sqrt(var)),
            "rms_over_span": float(np.sqrt(var)) / cfg.span,
            "act_frames": act, "n_frames": int(yc.size), "n_eff": n_eff,
            "statistic_fluctuates": fluct_ok,
            "statistic_fluctuates_msg": fluct_msg,
            "steps": {"equil": equil_steps, "prod": prod_steps, "stride": stride},
        }

    wall = time.perf_counter() - t0
    return {
        "config": asdict(cfg), "skipped": False,
        "span_L_star": cfg.span,
        "kappa_pred_star": cfg.kappa_pred_star,
        "ratio_meas_over_pred": out["kappa_star"] / cfg.kappa_pred_star,
        "dt_star": dt, "dt_info": dt_info,
        "wall_s": round(wall, 3),
        "guards": {"finite": not guard_fail, "failures": guard_fail,
                   "max_step_displacement_sigma": max_disp,
                   "bond_lengths": bond_info},
        "manifest": {
            "code_hash": _code_hash(), "git_rev": _git_rev(),
            "git_dirty": _git_dirty(),
            "hoomd_version": __import__("hoomd").version.version,
            "python": platform.python_version(), "platform": platform.platform(),
        },
        **out,
    }


# =============================================================================
def batch(configs: list[ChainConfig], concurrency: int = 8) -> dict:
    """독립 런 동시 실행. HOOMD 는 단일스레드이므로 이것이 유일한 병렬화다."""
    import subprocess

    t0 = time.perf_counter()
    pending = list(enumerate(configs))
    running: list = []
    done: list[dict] = []
    failed: list[dict] = []
    here = Path(__file__).resolve()

    while pending or running:
        while pending and len(running) < concurrency:
            idx, cfg = pending.pop(0)
            label = cfg.label or f"run{idx:03d}"
            p = subprocess.Popen(
                [sys.executable, str(here), "--worker", json.dumps(asdict(cfg))],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                cwd=str(here.parent.parent))
            running.append((p, label, cfg))
        fin = [(p, l, c) for p, l, c in running if p.poll() is not None]
        for p, label, cfg in fin:
            running.remove((p, label, cfg))
            so, se = p.communicate()
            line = next((x for x in so.splitlines() if x.strip().startswith("{")), None)
            if line is None:
                failed.append({"label": label, "config": asdict(cfg),
                               "stderr": se[-1500:]})
                print(f"  x {label}: {se.strip().splitlines()[-1] if se.strip() else 'no output'}",
                      flush=True)
                continue
            rec = json.loads(line)
            rec["label"] = label
            done.append(rec)
            if rec.get("skipped"):
                print(f"  - {label}: SKIP {rec['reason']}", flush=True)
            else:
                g = "" if rec["guards"]["finite"] and not rec["guards"]["failures"] \
                    else "  GUARD-FAIL"
                print(f"  o {label}: N={rec['config']['n_particles']:>3} "
                      f"kappa*={rec['kappa_star']:.5g}+-{rec['kappa_star_se']:.3g} "
                      f"pred={rec['kappa_pred_star']:.5g} "
                      f"ratio={rec['ratio_meas_over_pred']:.4f} "
                      f"dt={rec['dt_star']:.2e} ({rec['wall_s']:.1f}s){g}", flush=True)
        if running:
            time.sleep(0.15)

    return {"jobs": done, "failed": failed, "n_requested": len(configs),
            "batch_wall_s": round(time.perf_counter() - t0, 3),
            "concurrency": concurrency}


def main(argv: list[str]) -> int:
    if argv and argv[0] == "--worker":
        print(json.dumps(run_chain(ChainConfig(**json.loads(argv[1])))))
        return 0

    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--mode", default="static", choices=["static", "thermal", "both"])
    ap.add_argument("--tier", default="smoke", choices=["smoke", "pilot", "explore"])
    ap.add_argument("--seeds", type=int, default=1)
    ap.add_argument("--k-theta", type=float, default=1.0e4)
    ap.add_argument("--k-bond-ratio", type=float, default=100.0)
    ap.add_argument("--concurrency", type=int, default=8)
    ap.add_argument("--out", default="runs/chain-bend")
    a = ap.parse_args(argv)

    tiers = {"smoke": [5, 9], "pilot": [5, 7, 9, 11, 15],
             "explore": [5, 7, 9, 11, 15, 21, 31, 41]}
    ns = tiers[a.tier]
    modes = ["static", "thermal"] if a.mode == "both" else [a.mode]

    cfgs = []
    for m in modes:
        seeds = 1 if m == "static" else a.seeds     # 결정론은 시드가 무의미
        for n in ns:
            for s in range(1, seeds + 1):
                cfgs.append(ChainConfig(
                    n_particles=n, mode=m, seed=1000 + s,
                    label=f"{m[:4]}_N{n:03d}_s{s}",
                    k_theta_star=a.k_theta, k_bond_ratio=a.k_bond_ratio))

    print(f"tier={a.tier} mode={a.mode} N={ns} jobs={len(cfgs)} "
          f"k_theta*={a.k_theta:g} k_bond_ratio={a.k_bond_ratio:g}")
    res = batch(cfgs, concurrency=a.concurrency)
    out = Path(a.out) / f"{a.tier}-{a.mode}"
    out.mkdir(parents=True, exist_ok=True)
    (out / "batch.json").write_text(json.dumps(res, indent=2, ensure_ascii=False))
    print(f"\nbatch wall={res['batch_wall_s']}s ok={len(res['jobs'])} "
          f"failed={len(res['failed'])}\n-> {out/'batch.json'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
