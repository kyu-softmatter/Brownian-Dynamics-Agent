"""S5 — the HOOMD BD runner. The harmonic-trap system.

Artefacts:
  - `samples.npz`   : equilibrium position samples + the MSD time series (float32)
  - `manifest.json` : everything reproduction needs (spec/prediction hashes, seed,
                      versions, guard results)

Design principle: **the runner does not judge.** It produces numbers and reports
guard violations only. S7 judges, and a human confirms (CLAUDE.md §verdicts).
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
from bdbot.health import step_displacement_verdict
from .guards import check_finite, configurational_temperature
from .io import provenance as _provenance


@dataclass
class TrapRunConfig:
    """The complete specification of one harmonic-trap run. Reduced units
    (l_trap, kT, tau_trap)."""

    dim: int = 2
    n_particles: int = 1000
    dt_star: float = 5.0e-3
    equil_tau: float = 10.0            # equilibration length [tau_trap]
    prod_tau: float = 40.0             # production length [tau_trap]
    sample_interval_tau: float = 2.0   # independent-sample spacing [tau_trap]
    msd_frames: int = 500              # frames in the MSD time series
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
    # equilibrium samples (independent time points)
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
    # run information
    wall_s: float
    guards: dict
    manifest: dict


def run_trap(cfg: TrapRunConfig, *, outdir: Path | None = None,
             extra_manifest: dict | None = None) -> TrapRunResult:
    """One harmonic-trap BD run."""
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

    # --- equilibration ---
    sim.run(cfg.equil_steps)

    # --- production: collect the MSD series and independent samples together ---
    stride = cfg.msd_stride
    n_frames = cfg.prod_steps // stride
    traj = np.empty((n_frames, cfg.n_particles, cfg.dim), dtype=np.float32)

    indep_var, indep_kT = [], []
    indep_every = max(1, cfg.sample_interval_steps // stride)
    guard_fail: list[str] = []
    max_disp = 0.0
    max_force_star = 0.0
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
        #  ★ The force-driven step displacement, which this runner did not track.
        #    `max_disp` above is the **total** displacement over a stride --
        #    thermal + force, over many steps -- so bdbot's bound does not apply to
        #    it. The comparable quantity is `dt*|F|/gamma`, and in a harmonic trap
        #    `|F| = k*|r|` is computable from the positions we already have
        #    (reduced units: gamma = 1).
        #    Before this, `max_step_displacement_l_trap` went into the manifest and
        #    **nothing ever compared it to a threshold** -- the quiet-box-escape
        #    gap in .claude/rules/overdamped-stability.md. `check_finite` passes on
        #    a box escape.
        max_force_star = max(max_force_star,
                             cfg.k_star * float(np.linalg.norm(p, axis=1).max()))

        if f % indep_every == 0:
            indep_var.append(float(np.mean(p**2)))          # per-component <x*^2>
            indep_kT.append(configurational_temperature(
                cfg.k_star * p, laplacian_U_total=cfg.dim * cfg.k_star))

    wall = time.perf_counter() - t0

    # --- aggregation ---
    def mean_se(a: list[float]) -> tuple[float, float]:
        arr = np.asarray(a, dtype=np.float64)
        if arr.size < 2:
            return float(arr.mean()), float("nan")   # 1 sample claims no error
        return float(arr.mean()), float(arr.std(ddof=1) / np.sqrt(arr.size))

    var_c, var_c_se = mean_se(indep_var)
    kTc, kTc_se = mean_se(indep_kT)
    # <r*^2> = dim * <x*^2>  (isotropic)
    var_r, var_r_se = cfg.dim * var_c, cfg.dim * var_c_se

    # --- MSD (multiple time origins) ---
    max_lag = min(n_frames - 1, int(round(10.0 / (stride * cfg.dt_star))))
    lags = np.unique(np.geomspace(1, max(max_lag, 2), 40).astype(int))
    msd = []
    for lag in lags:
        d = traj[lag:] - traj[:-lag]
        msd.append(float(np.mean(np.sum(d.astype(np.float64) ** 2, axis=2))))

    manifest = {
        "run_hash": cfg.hash(),
        **_provenance(),                    # code_hash·git·env_hash·env — one place
        #  ★ `hoomd_version` and `python` are also inside `env`. They are kept
        #    because **the manifests already on disk use these keys** -- the readers
        #    (`report`, the analysis scripts) have to read old and new runs together.
        "hoomd_version": __import__("hoomd").version.version,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "seed": cfg.seed,
        "wall_s": round(wall, 3),
        "em_bias_expected": euler_maruyama_trap_variance_bias(cfg.dt_star),
        **(extra_manifest or {}),
    }
    #  ⚠ **Reports, does not raise** -- deliberately, and the asymmetry with
    #    `bdbot.run.StepGuard` is by design, not an oversight: this is a batch
    #    runner over seeds, and one diverging seed must not discard the others'
    #    work. What the two halves must share is the **bound**, not the reaction.
    #    `bdbot.health.step_displacement_verdict` is that one bound.
    force_disp_star = cfg.dt_star * max_force_star
    fd_ok, fd_why = step_displacement_verdict(force_disp_star, unit="l_trap")
    if not fd_ok:
        guard_fail.append(f"{fd_why} (|F*|max = {max_force_star:.4g} kT/l_trap)")
    guards = {
        "finite": not any("non-finite" in g for g in guard_fail),
        "failures": guard_fail,
        "max_step_displacement_l_trap": max_disp,
        "max_force_star": max_force_star,
        "force_displacement_star": force_disp_star,
        "force_displacement_ok": fd_ok,
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
# Batch — running independent runs concurrently
# =============================================================================
#  HOOMD is entirely single-threaded (TBB was removed in v3+, and this is not an MPI
#  build). So **running independent runs concurrently is the only parallelization
#  path**, which is also why they are launched as processes.
#  Basis: knowledge/wiki/findings/local-cpu-parallelism.md
def run_trap_batch(configs: list[TrapRunConfig], outroot: Path,
                   concurrency: int = 8, *, on_done=None) -> dict:
    """Run `configs` `concurrency` at a time. A failed run is recorded, not dropped.

    A run that goes missing quietly turns an error bar labelled "4 seeds" into one
    that actually had 3.
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
# The 2D soft-repulsive runner — `A/r^p`
# =============================================================================
#  Card: knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md
#  Reduced units: length `d = n^{-1/2}` · energy `kT` · time `tau_d = d^2/D0`
#  ⇒ kT* = 1, gamma* = 1, D* = 1, n* = 1, tau_d* = 1
@dataclass
class Soft2DRunConfig:
    """The complete specification of one `U = A/r^p` 2D run (reduced units)."""

    amplitude: float = 100.0            # A = beta*U(r=d)
    exponent: float = 3.0
    n_particles: int = 100
    density_star: float = 1.0           # 1 by definition (length unit is n^{-1/2})
    r_cut: float = 0.0                  # 0 means automatic, L/2 * 0.99
    r_min: float = 0.2                  # the Table's lower bound
    nlist_buffer: float = 0.1           # cell buffer (absolute). r_cut+buffer <= L/2
    init: str = "random"                # random | hex
    box_shape: str = "auto"             # auto | square | hex_commensurate
                                        # ★ auto branches on init → it contaminates
                                        #   an initial-condition comparison. State it
                                        #   explicitly when comparing (see below)
    dt_star: float = 2.0e-5
    equil_tau: float = 10.0             # [tau_d]
    prod_tau: float = 20.0
    n_frames: int = 200
    min_sep_init: float = 0.5
    max_tries_init: int = 20000         # rejection-sampling try limit. At n*=1,
                                        # min_sep=0.8 fails within 200 tries
                                        # (measured) — the last few are the hard ones
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
    """Snapshot + pair potential + box info. Shared by the runner and the `max|F|`
    measurement."""
    from .build import (hex_2d_snapshot, hex_tiling_for, random_2d_snapshot,
                        square_box_for)
    from .forces import power_law_table

    # ★★ The box shape is **decoupled** from the initial condition (2026-07-29).
    #  Originally they were tied together: `hex` → a hex-commensurate box (aspect
    #  1.1547), `random` → square (1.0). That made a hex vs random comparison also
    #  change
    #    ① the box aspect ratio  ② r_cut (derived from min(L)/2, so 4.46 vs 4.80)
    #  Measured: the S(k) 6-fold modulation at `A=0.1` came out 0.215–0.461 in the
    #  hex box and 0.004–0.008 in the square one -- **50x apart** -- because in a
    #  hex-commensurate box the Bragg k-vectors land exactly on the allowed k-grid
    #  (the measurement depends on the box shape).
    #  ⇒ To compare initial conditions, `box_shape` has to be **fixed the same**.
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
        raise ValueError(f"box_shape {shape!r} must be 'square' or "
                         f"'hex_commensurate'")

    if cfg.init == "hex":
        if shape != "hex_commensurate":
            raise ValueError(
                "a hexagonal initial placement is only perfect in a "
                "hex-commensurate box. "
                f"With box_shape={shape!r} the periodic boundary cuts the lattice "
                f"and creates artificial defects "
                "(see test_s5_pair.py::test_rotating_a_periodic_crystal...)")
        nx, ny = hex_tiling_for(cfg.n_particles)
        snap, geom = hex_2d_snapshot(hoomd, n_x=nx, n_y=ny,
                                     density_star=cfg.density_star)
    elif cfg.init == "random":
        snap = random_2d_snapshot(hoomd, n=cfg.n_particles, box_x=Lx, box_y=Ly,
                                  min_sep=cfg.min_sep_init, seed=cfg.seed,
                                  max_tries=cfg.max_tries_init)
    else:
        raise ValueError(f"init {cfg.init!r} must be 'random' or 'hex'")
    info = {**geom, "box_shape": shape}

    # ★ HOOMD requires `r_cut + buffer <= L/2` (not r_cut alone).
    #   r_cut is chosen leaving room for the buffer.
    half = min(Lx, Ly) / 2
    r_cut = cfg.r_cut or (0.98 * half - cfg.nlist_buffer)
    if r_cut + cfg.nlist_buffer > half:
        raise ValueError(
            f"r_cut({r_cut:.4f}) + buffer({cfg.nlist_buffer:.4f}) > L/2({half:.4f}) "
            f"— minimum-image violation. Raise N, or reduce r_cut/buffer")
    pair, pinfo = power_law_table(hoomd, amplitude=cfg.amplitude,
                                 exponent=cfg.exponent, r_cut=r_cut,
                                 r_min=cfg.r_min, buffer=cfg.nlist_buffer)
    return snap, pair, {**info, **pinfo}


def measure_max_force_soft2d(cfg: Soft2DRunConfig) -> dict:
    """**Computes the actual force** on the initial placement. Estimating is
    forbidden (master_plan §5.4).

    This value is needed to impose the force-displacement constraint
    `max|F|Δt ≤ δ_F`, and in this system **the force constraint beats the thermal
    displacement one** (the nearest-neighbour force at `A=100` is `~225 kT/d`).
    """
    import hoomd

    snap, pair, info = build_soft2d(hoomd, cfg)
    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=cfg.seed)
    sim.create_state_from_snapshot(snap)
    bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                                  default_gamma=1.0)
    sim.operations.integrator = hoomd.md.Integrator(dt=cfg.dt_star, methods=[bd],
                                                   forces=[pair])
    sim.run(0)                                   # make it compute the force
    f = np.asarray(pair.forces)
    mag = np.linalg.norm(f[:, :2], axis=1)
    return {"max_force_star": float(mag.max()), "mean_force_star": float(mag.mean()),
            "pair_energy_star": float(pair.energy), **info}


def run_soft2d(cfg: Soft2DRunConfig, *, outdir: Path | None = None,
               extra_manifest: dict | None = None) -> dict:
    """One `A/r^p` 2D BD run. **It does not judge** — it returns the trajectory and
    the guards."""
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

    # --- guard: did it leave the Table's lower bound ---
    from .analysis.structure import min_separation
    sep = min_separation(traj[:k + 1], Lx=Lx, Ly=Ly)
    if sep < cfg.r_min:
        guard_fail.append(
            f"minimum separation {sep:.4f} < r_min {cfg.r_min} — it left the Table. "
            f"U(r) may be quietly wrong")

    #  ★ This runner already **measured** `force_displacement_star` and never
    #    compared it to anything. Same bound as `bdbot.run.StepGuard`
    #    (`bdbot.health.step_displacement_verdict`), same quantity `dt*|F|max`,
    #    different reaction -- report, because this is a batch over seeds.
    #    ⚠ Judged on the **production** maximum, not the initial one. The initial
    #      configuration is not where the worst force occurs: measured elsewhere in
    #      this project, the peak was 1062.9 against 244.2 kT/sigma for the last
    #      sample, a factor of 4.4.
    force_disp_prod = float(max_force[:k + 1].max()) * cfg.dt_star
    fd_ok, fd_why = step_displacement_verdict(force_disp_prod)
    if not fd_ok:
        guard_fail.append(fd_why)
    guards = {
        "finite": not any("non-finite" in g for g in guard_fail),
        "failures": guard_fail,
        "min_separation": sep,
        "min_separation_over_r_min": sep / cfg.r_min,
        "max_force_initial": max_force_initial,
        "max_force_production": float(max_force[:k + 1].max()),
        "force_displacement_star": max_force_initial * cfg.dt_star,
        "force_displacement_production_star": force_disp_prod,
        "force_displacement_ok": fd_ok,
        "thermal_displacement_star": float(np.sqrt(2 * cfg.dt_star)),
        "n_frames": int(k + 1),
    }
    #  ★ Recording `hoomd` alone is not enough -- this run's guard
    #    (`min_separation`) and the analysis that follows use `freud`, `numpy` and
    #    `scipy`. A version bump can give a different measurement from the same
    #    trajectory, and with no record there is no way to know that.
    manifest = {
        "run_hash": cfg.hash(),
        **_provenance(),                    # code_hash·git·env_hash·env — one place
        "hoomd_version": __import__("hoomd").version.version,   # old-manifest compat
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
