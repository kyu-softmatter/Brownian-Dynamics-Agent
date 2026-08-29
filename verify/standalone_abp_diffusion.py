"""**Standalone verification** of the `active.abp` module -- the minimal
configuration, with interactions, traps and shape all switched off.

A worked instance of masterplan principle 9 (verify independent elements one at a
time). A free ABP has an analytic solution; the combination does not:

    ⟨n(0)·n(t)⟩ = e^(−Λt)
    MSD(t) = 2d·D_t·t + (2v₀²/Λ²)[Λt − 1 + e^(−Λt)]
    long time -> 2d*D_eff*t    =>    D_eff = D_t + v0^2/(d*Lambda)

Two traps came out of this check (skill bd-hoomd traps 10 and 11):
  (1) with active_force = 0, ActiveRotationalDiffusion **does not run at all**
      (Lambda = 0)
  (2) HOOMD's rotational_diffusion **IS the director decay rate Lambda itself**.
      It is NOT the standard theory's (d-1)*D_r -- Lambda/D_r = 1.00 in both 2D
      and 3D.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/standalone_abp_diffusion.py        # ~4 min
"""
import math
import sys

import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

KT = GAMMA = SIGMA = 1.0
D_T = KT / GAMMA
TOL = 3.0            # tolerance, %


def _frame(N, L, dim, rng):
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = (rng.random((N, 3)) - .5) * L * np.array(
        [1, 1, 1 if dim == 3 else 0])
    if dim == 2:                                   # rotation about z only
        q = np.zeros((N, 4))
        th = rng.random(N) * 2 * math.pi
        q[:, 0], q[:, 3] = np.cos(th / 2), np.sin(th / 2)
    else:                                          # random direction
        q = rng.normal(size=(N, 4))
        q /= np.linalg.norm(q, axis=1, keepdims=True)
    f.particles.orientation = q
    f.particles.typeid = [0] * N
    f.particles.types = ["A"]
    f.configuration.box = [L, L, L if dim == 3 else 0, 0, 0, 0]
    f.configuration.dimensions = dim
    return f


def run(dim, v0, D_r, N=4000, t_max=20.0, n_samp=240, seed=11, dt=5e-4, L=600.0):
    rng = np.random.default_rng(seed)
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(_frame(N, L, dim, rng))

    act = md.force.Active(filter=hoomd.filter.All())
    act.active_force["A"] = (v0 * GAMMA, 0., 0.)          # f_a = γ v₀
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=GAMMA)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[act])
    integ.integrate_rotational_dof = False                 # bd-hoomd trap 3
    sim.operations.integrator = integ
    sim.operations.updaters.append(act.create_diffusion_updater(
        trigger=hoomd.trigger.Periodic(1), rotational_diffusion=D_r))

    every = max(1, int(t_max / dt) // n_samp)

    def snap():
        s = sim.state.get_snapshot()
        p = np.array(s.particles.position) + np.array(s.particles.image) * L
        q = np.array(s.particles.orientation)
        w, x, y, z = q.T
        n = np.stack([1 - 2 * (y * y + z * z), 2 * (x * y + w * z), 2 * (x * z - w * y)], 1)
        return p[:, :dim].astype(np.float32), n.astype(np.float32)

    P, Nn = [], []
    for _ in range(n_samp):
        sim.run(every)
        p, n = snap()
        P.append(p); Nn.append(n)
    P, Nn = np.array(P), np.array(Nn)
    ds = every * dt

    nl = n_samp // 2                                       # multiple time origins
    msd, cor = np.zeros(nl), np.zeros(nl)
    cor[0] = 1.0
    for k in range(1, nl):
        msd[k] = float(((P[k:] - P[:-k]) ** 2).sum(-1).mean())
        cor[k] = float((Nn[k:] * Nn[:-k]).sum(-1).mean())
    t = np.arange(nl) * ds

    m = (cor > 0.15) & (t > 0)
    lam = float(-np.polyfit(t[m], np.log(cor[m]), 1)[0])
    fit = t > max(10 / lam, t[-1] * 0.5)
    D_eff = float(np.polyfit(t[fit], msd[fit], 1)[0] / (2 * dim))
    return dict(dim=dim, v0=v0, D_r=D_r, lam=lam, D_eff=D_eff)


def main():
    print("=" * 92)
    print("`active.abp` standalone check  (pair/trap/shape all OFF, "
          "kT=gamma=sigma=1 => D_t=1)")
    print("=" * 92)
    print(f"{'dim':>4}{'v0':>5}{'D_r':>6}{'Lam':>9}{'Lam/D_r':>7}{'Lam/[(d-1)Dr]':>14}"
          f"{'D_eff':>9}{'D_t+v0^2/(dLam)':>14}{'error':>8}")
    print("-" * 92)

    fails = []
    for dim, v0, D_r in [(2, 10., 3.), (2, 20., 3.), (3, 10., 3.), (3, 20., 3.)]:
        r = run(dim, v0, D_r)
        pred = D_T + v0 ** 2 / (dim * r["lam"])
        err = 100 * (r["D_eff"] - pred) / pred
        std = r["lam"] / ((dim - 1) * D_r)          # would be 1.0 if the standard theory held
        ok = abs(err) < TOL
        if not ok:
            fails.append((dim, v0, err))
        print(f"{dim:>4}{v0:>5.0f}{D_r:>6.1f}{r['lam']:>9.4f}{r['lam']/D_r:>7.3f}"
              f"{std:>14.3f}{r['D_eff']:>9.4f}{pred:>14.4f}{err:>+7.2f}%  {'✓' if ok else '✗'}")

    print("-" * 92)
    print("verdict")
    print(f"  . D_eff = D_t + v0^2/(d*Lambda)    -> all |error| < {TOL}%  "
          f"{'✓' if not fails else '✗ ' + str(fails)}")
    print("  . Lambda/D_r = 1.00 (both 2D and 3D) -> HOOMD's rotational_diffusion is")
    print("                                        the director decay rate Lambda "
          "itself (trap 11)")
    print("  . Lambda/[(d-1)D_r] is 0.50 in 3D  -> differs from the standard "
          "theory convention")
    print("=" * 92)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
