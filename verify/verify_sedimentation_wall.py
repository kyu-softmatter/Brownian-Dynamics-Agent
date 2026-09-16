"""Can HOOMD do sedimentation at all? Gravity, a wall, and the barometric law.

    PY=./bin/py
    $PY verify/verify_sedimentation_wall.py            # all stages
    $PY verify/verify_sedimentation_wall.py --quick     # skip the dt replicate

★ Why this exists. `docs/hoomd_capabilities.md` lists 15 measured API items and
**a wall is not among them** -- a grep for `hoomd.wall` over the whole repository
returned nothing before this file. `bd-hoomd` documents 20 traps and mentions
walls **zero times**, and its line 421 records that `md.external.field` carries
only `Electric`, `Magnetic` and `Periodic`, so there is no gravity class either.
Rules 4 and 6 point exactly here: do not write it from intuition, settle it by
execution.

Three stages, in the isolation order rule 7 demands; 3 is not interpreted unless
1 and 2 hold:

  1 · `md.force.Constant` under `Brownian`  ->  drift velocity `v = F/gamma`
  2 · a wall + `wall.Plane`                 ->  confinement, and does the wall
      reach through the periodic boundary?
  3 · both, with NO pair potential          ->  `n(z) ~ exp(-z/l_g)`, exact for
      point particles, and therefore the ground truth the interacting case is
      built on

## Two things this file measured the hard way

★ **`wall.Plane` does NOT apply the minimum image.** Its documented signed
distance is `d = n.(r - r_o)` with no wrapping, and the force on a particle at the
far top of a tall box is **exactly 0.0** (stage 2). That is the opposite of
`bd-hoomd` trap 1, where an external *trap* silently needed the minimum image
applied by hand (+1856 % without it). One external force needs it and the other
must not have it, so neither can be assumed from the other.

★ **A Lennard-Jones wall makes the WALL the constraint that binds `dt`, and the
symptom is a quiet pass-through, not NaN.** The first revision of this file used
`sigma_w = 0.5 d`, `epsilon_w = 1 kT`, cut at the LJ minimum. Measured:

      h/d     U [kT]      F [kT/d]    dt_max_force
      0.40    43.9        1517        2.0e-5
      0.47     3.6         141        2.1e-4
      0.50     1.0          48        6.3e-4

  At `dt = 1e-3` a particle at `h = 0.40` moves `F dt/gamma = 1.5 d` in one step,
  so it went **through** the wall and piled up against the periodic lid --
  `min(h) = -1.98 d`, with no NaN and no error. That is
  `.claude/rules/overdamped-stability.md`'s documented failure ("the symptom was
  not NaN but a quiet box escape"), now measured for a wall rather than a pair
  core. The fix is a **bounded** wall: `Gaussian` has `F_max = eps/(sigma sqrt(e))`
  at `h = sigma`, so the escape mode closes by construction rather than by
  choosing `dt` small enough to survive the tail of a divergence.

★ And equilibration is set by the **drift** across the box, not by `l_g^2/D_0`.
The first revision ran `12 l_g^2/D_0` in a box of height `30 l_g` and 56 % of the
particles were still above `10 l_g`: falling the height of the box takes
`L_z/v = L_z l_g/D_0 = 30 tau_sed`, not `tau_sed`. Stage 3 therefore does the
convergence half at a cheap `l_g` and the paper's `l_g` from the analytic start.

Units: `d = 1`, `kT = 1`, `gamma = 1`, so `D_0 = 1` and `tau_d = d^2/D_0 = 1`.
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from bdbot import dt as DT  # noqa: E402

# ── the system, from arXiv:1412.3190 ──────────────────────────────────────
D_UM = 0.8              # 2a, particle diameter [um]          -- paper
LG_UM = 10.7            # gravitational length [um]           -- paper, Table 1
LG_PAPER = LG_UM / D_UM  # = 13.375 d
LG_CHEAP = 4.0          # a second, cheaper l_g for the convergence half

# ── the wall: bounded on purpose (see the module docstring) ───────────────
EPS_W = 20.0            # U(0) = 20 kT  ->  penetration ~ e^-20 = 2e-9
SIG_W = 0.3
RCUT_W = 2.0            # U(2.0) = 4.5e-9 kT, so the fit window starts clean
F_MAX_W = EPS_W / (SIG_W * math.sqrt(math.e))   # Gaussian force peak, at h = sigma
FIT_LO = 2.0            # start the fit above the wall's reach


def _wall_force():
    w = md.external.wall.Gaussian(
        walls=[hoomd.wall.Plane(origin=(0, 0, 0), normal=(0, 0, 1))])
    w.params["A"] = {"epsilon": EPS_W, "sigma": SIG_W, "r_cut": RCUT_W}
    return w


def _make(N, Lxy, Lz, seed, *, lg, dt, start):
    """A box whose bottom is the wall (z_wall = 0, the box spans [0, Lz]).

    Putting the wall at z = 0 and the box at [0, Lz] rather than centring it
    keeps `h = z` and removes one chance to get a sign wrong.
    """
    dev = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=seed)
    rng = np.random.default_rng(seed)

    if start == "uniform":
        z = rng.uniform(FIT_LO, Lz - 1.0, size=N)
    elif start == "analytic":
        z = SIG_W + rng.exponential(lg, size=N)
        z = np.clip(z, 0.05, Lz - 1.0)
    else:
        raise ValueError(start)

    frame = gsd.hoomd.Frame()
    frame.particles.N = N
    frame.particles.position = np.column_stack([
        rng.uniform(-Lxy / 2, Lxy / 2, size=N),
        rng.uniform(-Lxy / 2, Lxy / 2, size=N),
        z - Lz / 2])                 # HOOMD centres the box on the origin
    frame.particles.typeid = np.zeros(N, dtype=int)
    frame.particles.types = ["A"]
    frame.configuration.box = [Lxy, Lxy, Lz, 0, 0, 0]
    sim.create_state_from_snapshot(frame)

    wall = md.external.wall.Gaussian(
        walls=[hoomd.wall.Plane(origin=(0, 0, -Lz / 2), normal=(0, 0, 1))])
    wall.params["A"] = {"epsilon": EPS_W, "sigma": SIG_W, "r_cut": RCUT_W}
    grav = md.force.Constant(filter=hoomd.filter.All())
    grav.constant_force["A"] = (0.0, 0.0, -1.0 / lg)
    grav.constant_torque["A"] = (0.0, 0.0, 0.0)

    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0)
    bd.gamma["A"] = 1.0
    sim.operations.integrator = md.Integrator(dt=dt, methods=[bd],
                                              forces=[wall, grav])
    return sim, Lz


# ── stage 1 ───────────────────────────────────────────────────────────────
def stage1_drift(N=2000, F=1.0, t_run=100.0, dt=1e-3, seed=11) -> dict:
    """A big `F` on purpose: this stage tests the API, not the physics, so the
    drift has to stand out of `sqrt(2 D t / N)`."""
    dev = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=seed)
    L = 1e4                                   # nothing wraps during the run
    frame = gsd.hoomd.Frame()
    frame.particles.N = N
    frame.particles.position = np.zeros((N, 3))
    frame.particles.typeid = np.zeros(N, dtype=int)
    frame.particles.types = ["A"]
    frame.configuration.box = [L, L, L, 0, 0, 0]
    sim.create_state_from_snapshot(frame)
    g = md.force.Constant(filter=hoomd.filter.All())
    g.constant_force["A"] = (0.0, 0.0, -F)
    g.constant_torque["A"] = (0.0, 0.0, 0.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0)
    bd.gamma["A"] = 1.0
    sim.operations.integrator = md.Integrator(dt=dt, methods=[bd], forces=[g])

    z0 = sim.state.get_snapshot().particles.position[:, 2].copy()
    sim.run(int(round(t_run / dt)))
    dz = sim.state.get_snapshot().particles.position[:, 2] - z0
    v = dz.mean() / t_run
    sem = dz.std(ddof=1) / math.sqrt(N) / t_run
    return {"v_meas": v, "sem": sem, "v_want": -F,
            "rel_err_pct": 100 * (v + F) / F, "n_sigma": abs(v + F) / sem}


# ── stage 2 ───────────────────────────────────────────────────────────────
def stage2_static() -> dict:
    """Read `wall.forces` directly, with no integration, at chosen heights.

    The last probe is the one that matters: a particle at the far top is
    `Lz - eps` above the wall and only `eps` below its periodic image.
    """
    Lz = 400.0
    heights = np.array([0.05, SIG_W, 1.0, RCUT_W * 1.01, 5.0, Lz - 2.0])
    sim, _ = _make(len(heights), 20.0, Lz, 7, lg=LG_PAPER, dt=1e-4,
                   start="uniform")
    wall = sim.operations.integrator.forces[0]
    snap = sim.state.get_snapshot()
    pos = np.zeros((len(heights), 3))
    pos[:, 2] = heights - Lz / 2
    snap.particles.position[:] = pos
    sim.state.set_snapshot(snap)
    sim.run(0)
    return {"Lz": Lz, "heights": heights,
            "fz": np.asarray(wall.forces)[:, 2].copy(),
            "U": np.asarray(wall.energies).copy()}


def stage2_confines(dt=5e-4, N=500, seed=13) -> dict:
    """Dynamically, WITH gravity pressing particles into the wall -- which is the
    configuration the first revision failed in."""
    lg = LG_CHEAP
    Lz = 15.0 * lg
    sim, _ = _make(N, 20.0, Lz, seed, lg=lg, dt=dt, start="analytic")
    t_run = 4.0 * lg ** 2
    n_steps = int(round(t_run / dt))
    hmin = np.inf
    for _ in range(40):
        sim.run(n_steps // 40)
        z = sim.state.get_snapshot().particles.position[:, 2] + Lz / 2
        hmin = min(hmin, float(z.min()))
    z = sim.state.get_snapshot().particles.position[:, 2] + Lz / 2
    return {"min_h_ever": hmin, "n_below_zero": int((z < 0).sum()), "N": N,
            "t_run": t_run, "dt": dt}


# ── stage 3 ───────────────────────────────────────────────────────────────
def _fit_decay(h, lo, hi, nbins=40):
    """Weighted least squares on `ln n(h)`. The SLOPE is fitted, so the origin --
    and therefore the wall's own structure layer -- cannot shift the answer."""
    cnt, edges = np.histogram(h, bins=np.linspace(lo, hi, nbins))
    mid = 0.5 * (edges[1:] + edges[:-1])
    ok = cnt >= 20
    if ok.sum() < 5:
        return None
    x, y, w = mid[ok], np.log(cnt[ok]), cnt[ok].astype(float)
    W = w.sum()
    xm, ym = (w * x).sum() / W, (w * y).sum() / W
    sxx = (w * (x - xm) ** 2).sum()
    slope = (w * (x - xm) * (y - ym)).sum() / sxx
    resid = y - (ym + slope * (x - xm))
    dof = max(1, int(ok.sum()) - 2)
    slope_se = math.sqrt((w * resid ** 2).sum() / dof / sxx)
    return {"lg": -1.0 / slope, "se": slope_se / slope ** 2,
            "nbins": int(ok.sum())}


def stage3(lg, dt, start, n_tau_sed, N=800, seed=17, n_samples=200) -> dict:
    Lz = 15.0 * lg                       # lid weight e^-15 = 3e-7
    sim, _ = _make(N, 20.0, Lz, seed, lg=lg, dt=dt, start=start)
    t_run = n_tau_sed * lg ** 2
    n_steps = int(round(t_run / dt))
    every = max(1, n_steps // n_samples)
    zs, t0 = [], time.time()
    for _ in range(n_steps // every):
        sim.run(every)
        zs.append(sim.state.get_snapshot().particles.position[:, 2].copy())
    wall_s = time.time() - t0
    h_all = np.concatenate([z + Lz / 2 for z in zs])
    h = np.concatenate([z + Lz / 2 for z in zs[len(zs) // 2:]])   # second half

    fit = _fit_decay(h, FIT_LO, min(8.0 * lg, Lz - 2.0))
    first = _fit_decay(np.concatenate([z + Lz / 2 for z in zs[:len(zs) // 2]]),
                       FIT_LO, min(8.0 * lg, Lz - 2.0))
    out = {"lg_want": lg, "dt": dt, "start": start, "t_run": t_run,
           "steps": n_steps, "wall_s": wall_s, "N": N,
           "min_h": float(h_all.min()),
           "leak_above_12lg": float((h > 12 * lg).mean())}
    if fit is None:
        out["fit"] = None
        return out
    out.update({"lg_fit": fit["lg"], "lg_se": fit["se"], "nbins": fit["nbins"],
                "rel_err_pct": 100 * (fit["lg"] - lg) / lg,
                "n_sigma": abs(fit["lg"] - lg) / fit["se"],
                "lg_first_half": first["lg"] if first else float("nan")})
    return out


# ── report ────────────────────────────────────────────────────────────────
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="skip the dt replicate")
    a = ap.parse_args()
    fails = []

    print(f"system: 2a = {D_UM} um, l_g = {LG_UM} um  ->  l_g* = {LG_PAPER:.4f} d")
    print("\nwhich constraint binds dt -- all three printed, per "
          "overdamped-stability.md:")
    print(f"  thermal, reference length d       {DT.dt_max_thermal(0.03, 1.0, 1.0):.3e}")
    print(f"  thermal, reference length l_g     {DT.dt_max_thermal(0.03, LG_PAPER, 1.0):.3e}")
    print(f"  gravity  F = kT/l_g = {1 / LG_PAPER:.5f}       "
          f"{DT.dt_max_force(0.03, 1.0, 1.0, 1 / LG_PAPER):.3e}")
    print(f"  ★ the WALL, F_max = {F_MAX_W:.1f} kT/d        "
          f"{DT.dt_max_force(0.03, 1.0, 1.0, F_MAX_W):.3e}   <- binding")

    print("\n" + "=" * 78)
    print("STAGE 1 · md.force.Constant under Brownian -> v = F/gamma")
    print("=" * 78)
    r = stage1_drift()
    print(f"  v_z = {r['v_meas']:.6f} +- {r['sem']:.6f}   analytic {r['v_want']:.6f}"
          f"   -> {r['rel_err_pct']:+.3f} % ({r['n_sigma']:.2f} sigma)")
    ok1 = r["n_sigma"] < 4.0
    print(f"  VERDICT: {'PASS' if ok1 else 'FAIL'}  (< 4 sigma)")
    ok1 or fails.append("stage 1")

    print("\n" + "=" * 78)
    print("STAGE 2 · Gaussian wall + wall.Plane")
    print("=" * 78)
    w = stage2_static()
    print(f"  eps_w = {EPS_W} kT, sigma_w = {SIG_W} d, r_cut = {RCUT_W} d, "
          f"box Lz = {w['Lz']:.0f} d")
    print(f"  {'h [d]':>12} {'U [kT]':>14} {'F_z [kT/d]':>14}")
    for h, u, f in zip(w["heights"], w["U"], w["fz"]):
        print(f"  {h:>12.4f} {u:>14.6e} {f:>14.6e}")
    pushes = w["fz"][0] > 0 and w["fz"][1] > 0
    zero_out = w["fz"][3] == 0.0 and w["fz"][4] == 0.0
    no_image = w["fz"][5] == 0.0
    print(f"  pushes outward inside the cutoff            : {pushes}")
    print(f"  exactly zero beyond the cutoff              : {zero_out}")
    print(f"  ★ exactly zero at the far top of the box    : {no_image}")
    print("    (so wall.Plane does NOT minimum-image -- the opposite of "
          "bd-hoomd trap 1)")
    c = stage2_confines()
    print(f"  dynamic, WITH gravity: {c['n_below_zero']}/{c['N']} below the wall "
          f"after {c['t_run']:.0f} tau_d; lowest h ever seen = {c['min_h_ever']:+.4f} d")
    ok2 = pushes and zero_out and no_image and c["n_below_zero"] == 0 \
        and c["min_h_ever"] > 0.0
    print(f"  VERDICT: {'PASS' if ok2 else 'FAIL'}")
    ok2 or fails.append("stage 2")

    print("\n" + "=" * 78)
    print("STAGE 3 · gravity + wall, no pair potential -> n(z) ~ exp(-z/l_g)")
    print("=" * 78)
    if not (ok1 and ok2):
        print("  SKIPPED -- stage 1 or 2 failed, so this could not be interpreted")
        fails.append("stage 3 not attempted")
    else:
        # ⚠ n_tau_sed has to clear the DRIFT time across the box, not tau_sed.
        #   A uniform start in a box of height 15 l_g needs 15 tau_sed just to
        #   fall, and the first attempt at 8 tau_sed read l_g = 11.5 against an
        #   imposed 4.0 (+187 %) with the two halves 32.5 and 11.5 -- still
        #   relaxing. 60 tau_sed is 4x the fall time.
        runs = [("converges from a uniform start", LG_CHEAP, 5e-4, "uniform", 60.0),
                ("stationary at the paper's l_g", LG_PAPER, 5e-4, "analytic", 5.0)]
        if not a.quick:
            runs.append(("dt replicate (l_g cheap)", LG_CHEAP, 2e-4, "uniform", 60.0))
        rows = []
        for label, lg, dt, start, ntau in runs:
            r3 = stage3(lg=lg, dt=dt, start=start, n_tau_sed=ntau)
            rows.append((label, r3))
            print(f"  {label}")
            print(f"    l_g* = {lg:.4f}, dt = {dt:.1e}, start = {start}, "
                  f"t = {r3['t_run']:.0f} tau_d ({r3['steps']:.3g} steps, "
                  f"{r3['wall_s']:.0f} s)")
            if r3.get("fit", True) is None:
                print("    FIT FAILED -- too few usable bins")
                continue
            print(f"    l_g fitted = {r3['lg_fit']:.4f} +- {r3['lg_se']:.4f} d   "
                  f"imposed {lg:.4f}   -> {r3['rel_err_pct']:+.3f} % "
                  f"({r3['n_sigma']:.2f} sigma, {r3['nbins']} bins)")
            print(f"    first half gave {r3['lg_first_half']:.4f} d "
                  f"(stationarity check)")
            print(f"    lowest h = {r3['min_h']:+.4f} d, leak above 12 l_g = "
                  f"{r3['leak_above_12lg']:.2e}")
        ok3 = all(r.get("n_sigma", 99) < 4.0 and r["min_h"] > 0.0
                  for _, r in rows)
        cheap = [r for lbl, r in rows if r["lg_want"] == LG_CHEAP and "lg_fit" in r]
        if len(cheap) == 2:
            spread = abs(cheap[0]["lg_fit"] - cheap[1]["lg_fit"]) / LG_CHEAP * 100
            print(f"  dt convergence: {spread:.3f} % of l_g between "
                  f"dt = {cheap[0]['dt']:.0e} and {cheap[1]['dt']:.0e}")
            ok3 = ok3 and spread < 3.0
        print(f"  VERDICT: {'PASS' if ok3 else 'FAIL'}  "
              f"(< 4 sigma, nothing below the wall, < 3 % over dt)")
        ok3 or fails.append("stage 3")

    print("\n" + "=" * 78)
    if fails:
        print("  FAILED:", ", ".join(fails))
        print("  Sedimentation is NOT yet a buildable system here.")
        return 1
    print("  All stages pass. Gravity (force.Constant) and a bounded bottom wall")
    print("  (external.wall.Gaussian + wall.Plane) are usable under Brownian, and")
    print("  the non-interacting profile reproduces exp(-z/l_g) -- the analytic")
    print("  ground truth the interacting case can be checked against.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
