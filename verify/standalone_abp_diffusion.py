"""`active.abp` 모듈 **단독 검증** — 상호작용·트랩·형상을 전부 끈 최소 구성.

마스터플랜 원칙 9(독립 요소는 하나씩 떼어 검증)의 실례.
자유 ABP에는 해석해가 있다 (조합하면 없다):

    ⟨n(0)·n(t)⟩ = e^(−Λt)
    MSD(t) = 2d·D_t·t + (2v₀²/Λ²)[Λt − 1 + e^(−Λt)]
    장시간 → 2d·D_eff·t    ⟹    D_eff = D_t + v₀²/(d·Λ)

이 검증에서 함정 2건이 나왔다 (skill bd-hoomd 함정 10·11):
  ① active_force = 0 이면 ActiveRotationalDiffusion 이 **아예 동작하지 않는다** (Λ=0)
  ② HOOMD 의 rotational_diffusion 은 **director 감쇠율 Λ 그 자체**다.
     표준 이론의 (d−1)·D_r 이 아니다 — 2D·3D 모두 Λ/D_r = 1.00.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/standalone_abp_diffusion.py        # ~4분
"""
import math
import sys

import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

KT = GAMMA = SIGMA = 1.0
D_T = KT / GAMMA
TOL = 3.0            # 허용 오차 %


def _frame(N, L, dim, rng):
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = (rng.random((N, 3)) - .5) * L * np.array(
        [1, 1, 1 if dim == 3 else 0])
    if dim == 2:                                   # z축 회전만
        q = np.zeros((N, 4))
        th = rng.random(N) * 2 * math.pi
        q[:, 0], q[:, 3] = np.cos(th / 2), np.sin(th / 2)
    else:                                          # 무작위 방향
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
    integ.integrate_rotational_dof = False                 # bd-hoomd 함정 3
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

    nl = n_samp // 2                                       # 다중 시간원점
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
    print("`active.abp` 단독 검증  (pair·trap·shape 전부 OFF, kT=γ=σ=1 ⟹ D_t=1)")
    print("=" * 92)
    print(f"{'dim':>4}{'v₀':>5}{'D_r':>6}{'Λ 측정':>9}{'Λ/D_r':>7}{'Λ/[(d−1)D_r]':>14}"
          f"{'D_eff':>9}{'D_t+v₀²/(dΛ)':>14}{'오차':>8}")
    print("-" * 92)

    fails = []
    for dim, v0, D_r in [(2, 10., 3.), (2, 20., 3.), (3, 10., 3.), (3, 20., 3.)]:
        r = run(dim, v0, D_r)
        pred = D_T + v0 ** 2 / (dim * r["lam"])
        err = 100 * (r["D_eff"] - pred) / pred
        std = r["lam"] / ((dim - 1) * D_r)          # 표준 이론이면 1.0 이어야 함
        ok = abs(err) < TOL
        if not ok:
            fails.append((dim, v0, err))
        print(f"{dim:>4}{v0:>5.0f}{D_r:>6.1f}{r['lam']:>9.4f}{r['lam']/D_r:>7.3f}"
              f"{std:>14.3f}{r['D_eff']:>9.4f}{pred:>14.4f}{err:>+7.2f}%  {'✓' if ok else '✗'}")

    print("-" * 92)
    print("판정")
    print(f"  · D_eff = D_t + v₀²/(d·Λ)          → 전부 |오차| < {TOL}%  "
          f"{'✓' if not fails else '✗ ' + str(fails)}")
    print("  · Λ/D_r = 1.00 (2D·3D 모두)        → HOOMD 의 rotational_diffusion 은")
    print("                                        director 감쇠율 Λ 그 자체 (함정 11)")
    print("  · Λ/[(d−1)D_r] 는 3D 에서 0.50     → 표준 이론 관례와 다름")
    print("=" * 92)
    return 1 if fails else 0


if __name__ == "__main__":
    sys.exit(main())
