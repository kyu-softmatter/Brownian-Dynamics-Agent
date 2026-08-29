"""S5 — HOOMD BD 러너. 조화 트랩 계.

산출물:
  - `samples.npz`   : 평형 위치 표본 + MSD 용 시계열 (float32)
  - `manifest.json` : 재현에 필요한 전부 (spec/예측 해시, seed, 버전, 가드 결과)

설계 원칙: **러너는 판정하지 않는다.** 수치를 내놓고 가드 위반만 보고한다.
판정은 S7 이 하고, 확정은 사람이 한다 (CLAUDE.md §판정).
"""
from __future__ import annotations

import hashlib
import json
import platform
import time
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np

from .build import trap_snapshot
from .estimators import euler_maruyama_trap_variance_bias
from .forces import HarmonicTrap
from .guards import check_finite, configurational_temperature
from .io import provenance as _provenance


@dataclass
class TrapRunConfig:
    """조화 트랩 런 1회의 완전 명세. 축약 단위 (l_trap, kT, tau_trap)."""

    dim: int = 2
    n_particles: int = 1000
    dt_star: float = 5.0e-3
    equil_tau: float = 10.0            # 평형화 길이 [tau_trap]
    prod_tau: float = 40.0             # 프로덕션 길이 [tau_trap]
    sample_interval_tau: float = 2.0   # 독립 표본 간격 [tau_trap]
    msd_frames: int = 500              # MSD 용 시계열 프레임 수
    box_over_l_trap: float = 200.0
    k_star: float = 1.0
    seed: int = 1
    label: str = ""

    @property
    def equil_steps(self) -> int:
        return int(round(self.equil_tau / self.dt_star))

    @property
    def prod_steps(self) -> int:
        return int(round(self.prod_tau / self.dt_star))

    @property
    def sample_interval_steps(self) -> int:
        return max(1, int(round(self.sample_interval_tau / self.dt_star)))

    @property
    def msd_stride(self) -> int:
        return max(1, self.prod_steps // self.msd_frames)

    def hash(self) -> str:
        payload = json.dumps(asdict(self), sort_keys=True).encode()
        return hashlib.sha256(payload).hexdigest()[:12]


@dataclass
class TrapRunResult:
    config: dict
    # 평형 표본 (독립 시점들)
    n_independent_snapshots: int
    var_per_component_star: float          # <x*^2>
    var_per_component_se: float
    var_radial_star: float                 # <r*^2>
    var_radial_se: float
    kT_conf_star: float
    kT_conf_se: float
    # MSD
    msd_lags_tau: list[float]
    msd_star: list[float]
    # 실행 정보
    wall_s: float
    guards: dict
    manifest: dict


def run_trap(cfg: TrapRunConfig, *, outdir: Path | None = None,
             extra_manifest: dict | None = None) -> TrapRunResult:
    """조화 트랩 BD 런 1회."""
    import hoomd

    t0 = time.perf_counter()
    dev = hoomd.device.CPU()
    sim = hoomd.Simulation(device=dev, seed=cfg.seed)
    sim.create_state_from_snapshot(
        trap_snapshot(hoomd, n=cfg.n_particles, dim=cfg.dim,
                      box_over_l_trap=cfg.box_over_l_trap))

    trap = HarmonicTrap(k_star=cfg.k_star,
                        active_axes=(True, True, cfg.dim == 3))
    bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                                  default_gamma=1.0)
    sim.operations.integrator = hoomd.md.Integrator(dt=cfg.dt_star, methods=[bd],
                                                   forces=[trap])

    def pos() -> np.ndarray:
        s = sim.state.get_snapshot()
        return np.array(s.particles.position[:, :cfg.dim], dtype=np.float64)

    # --- 평형화 ---
    sim.run(cfg.equil_steps)

    # --- 프로덕션: MSD 시계열 + 독립 표본 동시 수집 ---
    stride = cfg.msd_stride
    n_frames = cfg.prod_steps // stride
    traj = np.empty((n_frames, cfg.n_particles, cfg.dim), dtype=np.float32)

    indep_var, indep_kT = [], []
    indep_every = max(1, cfg.sample_interval_steps // stride)
    guard_fail: list[str] = []
    max_disp = 0.0
    prev = pos()

    for f in range(n_frames):
        sim.run(stride)
        p = pos()
        traj[f] = p.astype(np.float32)

        ok, fails = check_finite(position=p)
        if not ok:
            guard_fail += [f"frame {f}: {x}" for x in fails]
            break
        max_disp = max(max_disp, float(np.abs(p - prev).max()))
        prev = p

        if f % indep_every == 0:
            indep_var.append(float(np.mean(p**2)))          # 성분별 <x*^2>
            indep_kT.append(configurational_temperature(
                cfg.k_star * p, laplacian_U_total=cfg.dim * cfg.k_star))

    wall = time.perf_counter() - t0

    # --- 집계 ---
    def mean_se(a: list[float]) -> tuple[float, float]:
        arr = np.asarray(a, dtype=np.float64)
        if arr.size < 2:
            return float(arr.mean()), float("nan")   # 표본 1개면 오차를 주장할 수 없다
        return float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(arr.size))

    var_c, var_c_se = mean_se(indep_var)
    kTc, kTc_se = mean_se(indep_kT)
    # <r*^2> = dim * <x*^2>  (등방)
    var_r, var_r_se = cfg.dim * var_c, cfg.dim * var_c_se

    # --- MSD (다중 시간원점) ---
    max_lag = min(n_frames - 1, int(round(10.0 / (stride * cfg.dt_star))))
    lags = np.unique(np.geomspace(1, max(max_lag, 2), 40).astype(int))
    msd = []
    for lag in lags:
        d = traj[lag:] - traj[:-lag]
        msd.append(float(np.mean(np.sum(d.astype(np.float64) ** 2, axis=2))))

    manifest = {
        "run_hash": cfg.hash(),
        **_provenance(),                    # code_hash·git·env_hash·env — 한 곳에서
        #  ★ `hoomd_version`·`python` 은 `env` 안에도 있다. 남겨 두는 이유는
        #    **이미 디스크에 있는 manifest 들이 이 키를 쓰기 때문**이다 —
        #    리더(`report`·분석 스크립트)가 옛 런과 새 런을 함께 읽어야 한다.
        "hoomd_version": __import__("hoomd").version.version,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "seed": cfg.seed,
        "wall_s": round(wall, 3),
        "em_bias_expected": euler_maruyama_trap_variance_bias(cfg.dt_star),
        **(extra_manifest or {}),
    }
    guards = {
        "finite": not guard_fail,
        "failures": guard_fail,
        "max_step_displacement_l_trap": max_disp,
        "n_independent_snapshots": len(indep_var),
    }

    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(
            outdir / "samples.npz", traj=traj,
            lags_steps=lags * stride, msd=np.asarray(msd),
            indep_var=np.asarray(indep_var), indep_kT=np.asarray(indep_kT))
        (outdir / "manifest.json").write_text(
            json.dumps({"config": asdict(cfg), "manifest": manifest,
                        "guards": guards}, indent=2, ensure_ascii=False))

    return TrapRunResult(
        config=asdict(cfg),
        n_independent_snapshots=len(indep_var),
        var_per_component_star=var_c, var_per_component_se=var_c_se,
        var_radial_star=var_r, var_radial_se=var_r_se,
        kT_conf_star=kTc, kT_conf_se=kTc_se,
        msd_lags_tau=(lags * stride * cfg.dt_star).tolist(),
        msd_star=msd, wall_s=wall, guards=guards, manifest=manifest,
    )


# =============================================================================
# 배치 — 독립 런 동시 실행
# =============================================================================
#  HOOMD 는 완전 단일스레드다 (v3+ 에서 TBB 제거, MPI 빌드 없음). 따라서 **독립 런
#  동시 실행이 유일한 병렬화 경로**다. 프로세스로 띄우는 이유도 그것이다.
#  근거: knowledge/wiki/findings/local-cpu-parallelism.md
def run_trap_batch(configs: list[TrapRunConfig], outroot: Path,
                   concurrency: int = 8, *, on_done=None) -> dict:
    """`configs` 를 동시 `concurrency` 개씩 실행. 실패한 런은 버리지 않고 기록한다.

    조용히 빠진 런이 있으면 "시드 4개"라고 적힌 오차막대가 실제로는 3개짜리가 된다.
    """
    import subprocess
    import sys

    outroot = Path(outroot)
    outroot.mkdir(parents=True, exist_ok=True)
    t0 = time.perf_counter()
    pending = list(enumerate(configs))
    running: list = []
    done: list[dict] = []
    failed: list[dict] = []

    while pending or running:
        while pending and len(running) < concurrency:
            idx, cfg = pending.pop(0)
            label = cfg.label or f"run{idx:03d}"
            out = outroot / label
            p = subprocess.Popen(
                [sys.executable, "-m", "simbot.run", "--worker",
                 json.dumps(asdict(cfg)), str(out)],
                stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                cwd=str(Path(__file__).parent.parent))
            running.append((p, label, cfg))
        finished = [(p, l, c) for p, l, c in running if p.poll() is not None]
        for p, label, cfg in finished:
            running.remove((p, label, cfg))
            so, se = p.communicate()
            line = next((x for x in so.splitlines() if x.strip().startswith("{")), None)
            if line is None:
                failed.append({"label": label, "config": asdict(cfg),
                               "stderr": se[-2000:]})
                continue
            rec = json.loads(line)
            rec["label"] = label
            done.append(rec)
            if on_done:
                on_done(rec, len(done), len(configs))
        if running:
            time.sleep(0.15)

    return {"jobs": done, "failed": failed, "n_requested": len(configs),
            "batch_wall_s": round(time.perf_counter() - t0, 3),
            "concurrency": concurrency}


def _worker_main(argv: list[str]) -> int:
    cfg = TrapRunConfig(**json.loads(argv[0]))
    r = run_trap(cfg, outdir=Path(argv[1]))
    print(json.dumps({
        "config": asdict(cfg),
        "var_c": r.var_per_component_star, "var_c_se": r.var_per_component_se,
        "var_r": r.var_radial_star, "var_r_se": r.var_radial_se,
        "kT_conf": r.kT_conf_star, "kT_conf_se": r.kT_conf_se,
        "wall_s": r.wall_s, "guards": r.guards,
        "n_snap": r.n_independent_snapshots, "outdir": argv[1],
    }))
    return 0


if __name__ == "__main__":
    import sys
    if len(sys.argv) > 1 and sys.argv[1] == "--worker":
        raise SystemExit(_worker_main(sys.argv[2:]))
    raise SystemExit("usage: python -m simbot.run --worker <config-json> <outdir>")


# =============================================================================
# 2D 소프트 반발계 러너 — `A/r^p`
# =============================================================================
#  카드: knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md
#  축약 단위: 길이 `d = n^{-1/2}` · 에너지 `kT` · 시간 `tau_d = d^2/D0`
#  ⇒ kT* = 1, gamma* = 1, D* = 1, n* = 1, tau_d* = 1
@dataclass
class Soft2DRunConfig:
    """`U = A/r^p` 2D 계 런 1회의 완전 명세 (축약 단위)."""

    amplitude: float = 100.0            # A = beta*U(r=d)
    exponent: float = 3.0
    n_particles: int = 100
    density_star: float = 1.0           # 정의상 1 (길이 단위가 n^{-1/2})
    r_cut: float = 0.0                  # 0 이면 L/2 * 0.99 로 자동
    r_min: float = 0.2                  # Table 하한
    nlist_buffer: float = 0.1           # Cell 버퍼 (절대 거리). r_cut+buffer <= L/2
    init: str = "random"                # random | hex
    box_shape: str = "auto"             # auto | square | hex_commensurate
                                        # ★ auto 는 init 에 따라 갈린다 → 초기조건 비교를
                                        #   교란한다. 비교할 때는 명시할 것 (아래 참조)
    dt_star: float = 2.0e-5
    equil_tau: float = 10.0             # [tau_d]
    prod_tau: float = 20.0
    n_frames: int = 200
    min_sep_init: float = 0.5
    max_tries_init: int = 20000         # 기각표집 시도 한도. n*=1 에서 min_sep=0.8 은
                                        # 200 회로 실패한다 (실측) — 마지막 몇 개가 어렵다
    seed: int = 1
    label: str = ""

    @property
    def equil_steps(self) -> int:
        return int(round(self.equil_tau / self.dt_star))

    @property
    def prod_steps(self) -> int:
        return int(round(self.prod_tau / self.dt_star))

    @property
    def stride(self) -> int:
        return max(1, self.prod_steps // self.n_frames)

    def hash(self) -> str:
        return hashlib.sha256(
            json.dumps(asdict(self), sort_keys=True).encode()).hexdigest()[:12]


def build_soft2d(hoomd, cfg: Soft2DRunConfig):
    """스냅샷 + 쌍 포텐셜 + 상자 정보. 러너와 `max|F|` 측정이 공유한다."""
    from .build import (hex_2d_snapshot, hex_tiling_for, random_2d_snapshot,
                        square_box_for)
    from .forces import power_law_table

    # ★★ 상자 모양을 초기조건과 **분리한다** (2026-07-29).
    #  처음에는 `hex` → 육방정합 상자(종횡비 1.1547), `random` → 정사각(1.0) 으로
    #  묶어 두었다. 그러면 hex vs random 비교가 초기조건 외에
    #    ① 상자 종횡비  ② r_cut (= min(L)/2 에서 파생되므로 4.46 vs 4.80)
    #  까지 함께 바꾼다. 실측 결과 `A=0.1` 의 S(k) 6겹 변조가 hex 상자에서
    #  0.215–0.461, 정사각에서 0.004–0.008 로 **50배** 갈렸다 — 육방정합 상자에서는
    #  Bragg k-벡터가 허용 k-격자에 정확히 놓이기 때문이다 (측정량이 상자 모양에 의존).
    #  ⇒ 초기조건을 비교하려면 `box_shape` 를 **같게 고정**해야 한다.
    shape = cfg.box_shape
    if shape == "auto":
        shape = "hex_commensurate" if cfg.init == "hex" else "square"

    if shape == "hex_commensurate":
        nx, ny = hex_tiling_for(cfg.n_particles)
        _, geom = hex_2d_snapshot(hoomd, n_x=nx, n_y=ny,
                                  density_star=cfg.density_star)
        Lx, Ly = geom["Lx"], geom["Ly"]
    elif shape == "square":
        Lx = Ly = square_box_for(cfg.n_particles, cfg.density_star)
        geom = {"Lx": Lx, "Ly": Ly, "a_nn": None, "aspect": 1.0,
                "n_particles": cfg.n_particles, "density_star": cfg.density_star}
    else:
        raise ValueError(f"box_shape {shape!r} 는 'square' 또는 "
                         f"'hex_commensurate' 여야 한다")

    if cfg.init == "hex":
        if shape != "hex_commensurate":
            raise ValueError(
                "육방 초기배치는 육방정합 상자에서만 완벽하다. "
                f"box_shape={shape!r} 로는 주기경계가 격자를 끊어 인공 결함이 생긴다 "
                "(test_s5_pair.py::test_rotating_a_periodic_crystal... 참조)")
        nx, ny = hex_tiling_for(cfg.n_particles)
        snap, geom = hex_2d_snapshot(hoomd, n_x=nx, n_y=ny,
                                     density_star=cfg.density_star)
    elif cfg.init == "random":
        snap = random_2d_snapshot(hoomd, n=cfg.n_particles, box_x=Lx, box_y=Ly,
                                  min_sep=cfg.min_sep_init, seed=cfg.seed,
                                  max_tries=cfg.max_tries_init)
    else:
        raise ValueError(f"init {cfg.init!r} 는 'random' 또는 'hex' 여야 한다")
    info = {**geom, "box_shape": shape}

    # ★ HOOMD 는 `r_cut + buffer <= L/2` 를 요구한다 (r_cut 만이 아니다).
    #   버퍼 자리를 남기고 r_cut 을 정한다.
    half = min(Lx, Ly) / 2
    r_cut = cfg.r_cut or (0.98 * half - cfg.nlist_buffer)
    if r_cut + cfg.nlist_buffer > half:
        raise ValueError(
            f"r_cut({r_cut:.4f}) + buffer({cfg.nlist_buffer:.4f}) > L/2({half:.4f}) "
            f"— 최소이미지 위반. N 을 늘리거나 r_cut/buffer 를 줄일 것")
    pair, pinfo = power_law_table(hoomd, amplitude=cfg.amplitude,
                                 exponent=cfg.exponent, r_cut=r_cut,
                                 r_min=cfg.r_min, buffer=cfg.nlist_buffer)
    return snap, pair, {**info, **pinfo}


def measure_max_force_soft2d(cfg: Soft2DRunConfig) -> dict:
    """초기배치에서 **실제 힘을 계산한다.** 추정 금지 (master_plan §5.4).

    힘 변위 제약 `max|F|Δt ≤ δ_F` 를 걸려면 이 값이 필요하고, 이 계에서는
    **힘 제약이 열 변위 제약을 이긴다** (`A=100` 에서 최근접 힘이 `~225 kT/d`).
    """
    import hoomd

    snap, pair, info = build_soft2d(hoomd, cfg)
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=cfg.seed)
    sim.create_state_from_snapshot(snap)
    bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                                  default_gamma=1.0)
    sim.operations.integrator = hoomd.md.Integrator(dt=cfg.dt_star, methods=[bd],
                                                   forces=[pair])
    sim.run(0)                                   # 힘을 계산시킨다
    f = np.asarray(pair.forces)
    mag = np.linalg.norm(f[:, :2], axis=1)
    return {"max_force_star": float(mag.max()), "mean_force_star": float(mag.mean()),
            "pair_energy_star": float(pair.energy), **info}


def run_soft2d(cfg: Soft2DRunConfig, *, outdir: Path | None = None,
               extra_manifest: dict | None = None) -> dict:
    """`A/r^p` 2D 계 BD 런 1회. **판정하지 않는다** — 궤적과 가드만 돌려준다."""
    import hoomd

    t0 = time.perf_counter()
    snap, pair, info = build_soft2d(hoomd, cfg)
    Lx, Ly = info["Lx"], info["Ly"]

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=cfg.seed)
    sim.create_state_from_snapshot(snap)
    bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                                  default_gamma=1.0)
    sim.operations.integrator = hoomd.md.Integrator(dt=cfg.dt_star, methods=[bd],
                                                   forces=[pair])
    sim.run(0)
    f0 = np.linalg.norm(np.asarray(pair.forces)[:, :2], axis=1)
    max_force_initial = float(f0.max())

    def pos() -> np.ndarray:
        s = sim.state.get_snapshot()
        return np.array(s.particles.position[:, :2], dtype=np.float64)

    init_pos = pos()
    sim.run(cfg.equil_steps)

    stride = cfg.stride
    n_frames = cfg.prod_steps // stride
    traj = np.empty((n_frames, cfg.n_particles, 2), dtype=np.float32)
    energy = np.empty(n_frames)
    max_force = np.empty(n_frames)
    guard_fail: list[str] = []

    for k in range(n_frames):
        sim.run(stride)
        p = pos()
        traj[k] = p.astype(np.float32)
        energy[k] = float(pair.energy)
        max_force[k] = float(np.linalg.norm(
            np.asarray(pair.forces)[:, :2], axis=1).max())
        ok, fails = check_finite(position=p)
        if not ok:
            guard_fail += [f"frame {k}: {x}" for x in fails]
            break
    wall = time.perf_counter() - t0

    # --- 가드: Table 하한을 벗어났는가 ---
    from .analysis.structure import min_separation
    sep = min_separation(traj[:k + 1], Lx=Lx, Ly=Ly)
    if sep < cfg.r_min:
        guard_fail.append(
            f"최소분리 {sep:.4f} < r_min {cfg.r_min} — Table 표를 벗어났다. "
            f"U(r) 가 조용히 틀렸을 수 있다")

    guards = {
        "finite": not any("비유한값" in g for g in guard_fail),
        "failures": guard_fail,
        "min_separation": sep,
        "min_separation_over_r_min": sep / cfg.r_min,
        "max_force_initial": max_force_initial,
        "max_force_production": float(max_force[:k + 1].max()),
        "force_displacement_star": max_force_initial * cfg.dt_star,
        "thermal_displacement_star": float(np.sqrt(2 * cfg.dt_star)),
        "n_frames": int(k + 1),
    }
    #  ★ `hoomd` 만 기록하면 부족하다 — 이 런의 가드(`min_separation`)와 뒤따르는
    #    분석이 `freud`·`numpy`·`scipy` 를 쓴다. 버전이 올라가면 같은 궤적에서
    #    다른 측정값이 나올 수 있고, 기록이 없으면 그 사실을 알 수 없다.
    manifest = {
        "run_hash": cfg.hash(),
        **_provenance(),                    # code_hash·git·env_hash·env — 한 곳에서
        "hoomd_version": __import__("hoomd").version.version,   # 옛 manifest 호환
        "python": platform.python_version(), "seed": cfg.seed,
        "wall_s": round(wall, 3), **(extra_manifest or {}),
    }
    out = {"config": asdict(cfg), "info": info, "guards": guards,
           "manifest": manifest, "wall_s": wall,
           "energy_mean": float(energy[:k + 1].mean()),
           "energy_per_particle": float(energy[:k + 1].mean() / cfg.n_particles)}

    if outdir is not None:
        outdir = Path(outdir)
        outdir.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(outdir / "samples.npz", traj=traj[:k + 1],
                            energy=energy[:k + 1], max_force=max_force[:k + 1],
                            init_pos=init_pos,
                            box=np.array([Lx, Ly]),
                            stride=np.array([stride]))
        (outdir / "manifest.json").write_text(
            json.dumps(out, indent=2, ensure_ascii=False, default=float))
    return out
