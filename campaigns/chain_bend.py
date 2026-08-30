"""Colloidal chain bending stiffness -- a sweep over chain length N.

Goal: confirm both the **exponent and the coefficient** of `kappa(N) ~ L^-3`.
The BD counterpart of Pantina & Furst (PRL 94, 138301)
`kappa = kappa_0 (a/s)^(2+d_b)` with `d_b=1`.
**The exponent -3 is convention-independent; the coefficient is convention-dependent**
-- the paper PDF's coefficient is corrupted by two-column extraction (env-log stage 3)
and cannot be trusted, so textbook beam theory is used as the reference instead.

Reduced units: length sigma = 2a = bond length = 1, energy kT = 1,
time tau_D = sigma^2/D0
(gamma* = 1, D0* = 1).

## Prediction (simply-supported beam + the discrete-chain mapping)

    discrete bending energy  U = 1/2 k_th sum (dtheta_i)^2
    continuum counterpart    EI = k_th * b,   b = bond length = 1
    centre point load, simply supported at both ends:  delta = F L^3 / (48 EI)

    ->  kappa* = 48 k_theta* / (N-1)^3

The end particles have their positions fixed but **no angular constraint** (the
vertices of the angle triplets are only 1..N-2), so the boundary is pinned (simply
supported), not clamped. That is why the coefficient is 48.

## Both protocols are run together

- `mode="static"`  : **deterministic.** kT=0; apply a known force F* to the centre
  particle and relax -> kappa = F/delta. With no noise this judges the mapping (the
  coefficient 48) sharply. No statistics needed.
- `mode="thermal"` : fluctuations. Equipartition of the centre particle's deflection
  -> kappa* = 1/<y_c*^2>. This yields error bars and is the route that later extends
  to G'(omega).

If the two agree, the mapping is right. If they disagree, which one is wrong is
separable.

⚠ Distance constraints are NOT used -- they are incompatible with `Brownian` and
diverge silently.
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

LAMBDA = 6.3e6            # throughput constant [particle*steps/s],
                          # config/run_policy.yaml
WALL_BUDGET_S = 600.0     # CLAUDE.md: 10 min/run
MAX_FRAMES = 8000         # cap on the cost of get_snapshot calls

# ★ The dt gate's formulas are owned by `simbot.nondim`, its thresholds by
#   `config/run_policy.yaml`. If this script restated `0.03`, `0.005` or `0.2`, editing
#   the policy would no longer propagate here.
#   Since the reduced units are sigma = kT = gamma = D0 = 1, passing 1.0 to the gate
#   functions returns the reduced bound directly (the gate formulas are
#   unit-independent).
_TS = load_policy().timestep


# =============================================================================
@dataclass
class ChainConfig:
    n_particles: int = 11             # odd only (a centre particle must be defined)
    dim: int = 2
    mode: str = "thermal"             # "static" | "thermal"
    k_theta_star: float = 1.0e4       # angle spring [kT/rad^2]
    k_bond_ratio: float = 100.0       # k_bond* = ratio * kappa(N)*. See findings:
                                      # holding it constant makes the step count grow
                                      # as L^3
    delta_over_span: float = 0.01     # static: target deflection delta/L (stays
                                      # linear)
    dt_star: float = 0.0              # 0 means: set it from the measured force via
                                      # the gate
    relax_tau_bend: float = 30.0      # static: relaxation length [tau_bend]
    equil_tau_bend: float = 30.0      # thermal: equilibration length
    prod_tau_bend: float = 400.0      # thermal: production length
    frames_per_tau_bend: float = 5.0  # thermal: sampled frames per tau_bend
    seed: int = 1
    label: str = ""

    def __post_init__(self):
        if self.n_particles % 2 == 0:
            raise ValueError(f"n={self.n_particles} -- odd only, a centre particle "
                             f"is required")
        if self.n_particles < 5:
            raise ValueError(f"n={self.n_particles} -- at least 3 interior degrees of "
                             f"freedom are needed")
        if self.mode not in ("static", "thermal"):
            raise ValueError(f"mode={self.mode!r}")

    @property
    def span(self) -> float:
        """Support-to-support distance L* = (N-1) * bond length."""
        return float(self.n_particles - 1)

    @property
    def kappa_pred_star(self) -> float:
        return 48.0 * self.k_theta_star / self.span ** 3

    @property
    def tau_bend_star(self) -> float:
        """Estimated bending relaxation time gamma*/kappa*  (gamma* = 1)."""
        return 1.0 / self.kappa_pred_star

    @property
    def k_bond_star(self) -> float:
        """Scale the bond spring with kappa(N), so dt grows as L^3 and the cost is
        nearly independent of N."""
        return self.k_bond_ratio * self.kappa_pred_star

    @property
    def lambda_max_star(self) -> float:
        """Estimated largest eigenvalue of the stiffness matrix.

        4k for a 1D spring chain plus 16k for bending (a fourth-difference operator).
        """
        return 4.0 * self.k_bond_star + 16.0 * self.k_theta_star

    @property
    def dt_stability_star(self) -> float:
        """The explicit-Euler stability limit `2/lambda_max`, times the policy's safety
        factor.

        The measured threshold is 1.22-2.80x this bound -> a safety factor of 0.2
        leaves 6-14x margin.
        Basis: findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
        ★ The displacement gates alone cannot catch this -- a straight chain has
        max|F*| = 0, which disables the force gate.
        """
        return dt_max_stability(_TS["stability_safety_factor"], 1.0,
                                self.lambda_max_star)

    @property
    def load_star(self) -> float:
        """The force applied in static mode. Since delta = F/kappa, F = kappa * delta.
        """
        return self.kappa_pred_star * self.delta_over_span * self.span


# =============================================================================
def _seed_dt(cfg: ChainConfig) -> float:
    """Initial estimate for dt. **The stability term is almost always the binding
    one.**

    Measuring the force on a straight chain gives exactly 0, which disables the force
    gate (measured 2026-07-28: dt came out as 4.5e-4 and the chain blew up even at
    kT=0).
    Basis: findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
    """
    f = [np.sqrt(cfg.k_bond_star), np.sqrt(cfg.k_theta_star), 1.0]
    if cfg.mode == "static":
        f.append(cfg.load_star)
    proxy_force = 3.0 * max(f)          # a temporary proxy before measurement;
                                        # replaced by the measured value after warm-up
    return min(_dt_gates(proxy_force, cfg).values())


def _dt_gates(max_f_star: float, cfg: ChainConfig) -> dict[str, float]:
    """Per-gate dt bounds (reduced units). A disabled gate has its **key omitted.**"""
    gates = {
        "diffusion": dt_max_thermal(_TS["max_thermal_displacement_sigma"], 1.0, 1.0),
        "force": dt_max_force(_TS["max_force_displacement_sigma"], 1.0, 1.0, max_f_star),
        "stability": cfg.dt_stability_star,
    }
    return {k: v for k, v in gates.items() if v is not None}


def _dt_from_force(max_f_star: float, cfg: ChainConfig) -> tuple[float, dict]:
    """The dt gate set -- two displacement gates (accuracy) plus one stability gate.

    The displacement gates are still required but **not sufficient**. A fixed
    `dt/tau_D` gate is not used.
    The formulas and thresholds are owned by `simbot.nondim` + `run_policy.yaml` (do
    not duplicate them).
    """
    gates = _dt_gates(max_f_star, cfg)
    binding = min(gates, key=gates.get)
    return gates[binding], {
        "dt_diffusion_gate": gates["diffusion"],
        # On a straight chain max|F*| = 0, so the force gate **does not exist** (None).
        # Papering over it with something like 1e300 would read as "the gate existed
        # but was loose", which is wrong.
        "dt_force_gate": gates.get("force"),
        "dt_stability_gate": gates["stability"],
        "lambda_max_star": cfg.lambda_max_star,
        "k_bond_star": cfg.k_bond_star,
        "max_force_star_measured": float(max_f_star),
        "binding": binding}


def _integrated_act(x: np.ndarray, c_window: float = 6.0) -> float:
    """Integrated autocorrelation time [frames]. Sokal automatic windowing.

    Bending modes are slow, so independence must not be assumed.
    """
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

    # The two ends (type "E") are not integrated -> perfectly fixed supports.
    # Holding them with traps would mix trap compliance into kappa (the paper's
    # 40 pN/um ceiling problem).
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

    # --- 1) fix dt from the measured force after warm-up ---
    dt_seed = sim.operations.integrator.dt
    warm = max(50, int(round(0.5 * tau_b / dt_seed)))
    sim.run(min(warm, 20000))
    max_f = _max_force(sim)
    # The maximum over several times (looking at one instant underestimates it)
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
        # **Adaptive relaxation.** A fixed length (30 tau_bend) does not work -- the
        # measured true relaxation time is not 1/kappa_point but roughly N times
        # longer, because the whole chain has to be dragged.
        # 2026-07-28: running at fixed length gave a kappa 22 % too high at N=41.
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
            guard_fail.append(f"bond-length violation "
                              f"max_rel_dev={bond_info['max_rel_dev']:.3e}")

        delta = float(p[ic, 1])
        c = np.asarray(curve, dtype=np.float64)
        # Residual rate of change over the final interval
        drift = float(abs(c[-1] - c[-2]) / max(abs(c[-1]), 1e-300)) if c.size > 1 else np.nan
        # Aitken extrapolation (3 equally spaced points) -> delta_inf. Estimates the
        # plateau even if the relaxation has not finished.
        delta_inf, tau_fit = float("nan"), float("nan")
        if c.size >= 9:
            i1, i2, i3 = c.size - 1 - 2 * (c.size // 8), c.size - 1 - (c.size // 8), c.size - 1
            d1, d2, d3 = c[i1], c[i2], c[i3]
            den = d1 + d3 - 2 * d2
            if abs(den) > 1e-300:
                delta_inf = float((d1 * d3 - d2 * d2) / den)
            # Log-linear fit of delta(t) = d_inf (1 - A exp(-t/tau))
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
            "kappa_star": kappa, "kappa_star_se": 0.0,   # deterministic -- no
                                                         # statistical error
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
                    "reason": f"estimated wall {est:.0f}s > budget "
                              f"{WALL_BUDGET_S:.0f}s",
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
                        f"frame {f}: bond-length violation "
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
            assert_statistic_fluctuates(blocks, name="<y_c^2> block means")
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
    """Run independent runs concurrently. HOOMD is single-threaded, so this is the only
    parallelism available."""
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
        seeds = 1 if m == "static" else a.seeds     # seeds are meaningless when
                                                    # deterministic
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
