#!/usr/bin/env python
"""`network` 케이스 실행 전 능력조사 — 3D BD + `update.BoxResize` (CLAUDE.md 규칙 4·7).

왜 필요한가: 이 프로젝트의 기존 6개 케이스는 **전부 2D** 이고, 박스 압축
(`hoomd.update.BoxResize`)은 **한 번도 써본 적이 없습니다** (skill `bd-hoomd` ·
docs/hoomd_capabilities.md 에 항목 없음 → intake/network/observation.yaml N4·N6).

검사 5종 —
  ① 3D 자유확산 단독:  ⟨r²⟩ = 6·D·t          (격리 검증. N4)
  ② BoxResize 가 입자 좌표를 정말 **아핀 스케일** 하는가 (문서 주장의 실측)
  ③ BoxResize + Brownian + pair.Table + WCA 가 함께 도는가 (셀리스트·유한성)
  ④ ★ **압축이 DLVO 결합을 부수는 문턱** — 이게 겔화 프로토콜의 설계 수치다.
       예측: 트리거당 선형변형 ε_crit = (h_min* − barrier_h*)/ℓ*
       (아핀 스텝이 결합을 장벽 안쪽으로 밀어넣으면 1차극소로 떨어져 비가역)
  ⑤ 비용: 겔화 구성의 steps/s 와 그로부터 나오는 벽시계 (N3 — 판정 전에 재라)

단위: d=1, kT=1, γ=1  ⟹ D_t=1, τ_B = d²/D_t = 1.
DLVO 축약 파라미터는 `cases/chain_bend_dlvo_2d.py` 에서 그대로 가져옵니다
(같은 물리를 두 번 적지 않기 위해 — 그 케이스에서 SI로 검증된 식입니다).

실행:  $PY scratch/verify_3d_boxresize.py [--quick]
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

import gsd.hoomd                                        # noqa: E402
import hoomd                                            # noqa: E402
import hoomd.md as md                                   # noqa: E402

from chain_bend_dlvo_2d import (                         # noqa: E402
    SIGMA_CORE_STAR, CUTOFF_H_STAR,
    build_table_arrays, dlvo_reduced_params, find_well, load_system,
)

R_WCA = 2 ** (1 / 6)
PASS, FAIL = [], []


def check(ok: bool, label: str, detail: str = "") -> None:
    (PASS if ok else FAIL).append(label)
    print(f"  {'✓' if ok else '✗'} {label}" + (f"   {detail}" if detail else ""))


def cpu(seed: int) -> hoomd.Simulation:
    # ★ seed < 65536 (bd-hoomd 함정 12 — 16비트로 잘린다)
    return hoomd.Simulation(device=hoomd.device.CPU(), seed=seed)


def frame_3d(pos, L, types=("A",), typeid=None):
    f = gsd.hoomd.Frame()
    pos = np.asarray(pos, dtype=float)
    f.particles.N = len(pos)
    f.particles.position = pos
    f.particles.typeid = [0] * len(pos) if typeid is None else list(typeid)
    f.particles.types = list(types)
    f.particles.mass = [1.0] * len(pos)
    f.configuration.box = [L, L, L, 0, 0, 0]      # ★ 3D: Lz=L (2D 는 Lz=0 — 함정 9)
    f.configuration.dimensions = 3
    return f


def unwrapped(sim, L=None):
    s = sim.state.get_snapshot()
    pos = np.array(s.particles.position, copy=True)
    img = np.array(s.particles.image, copy=True)
    box = s.configuration.box
    Ls = np.array([box[0], box[1], box[2]], dtype=float) if L is None else np.full(3, L)
    return pos + img * Ls


# ═══════════════════════════════════════════════════════════════════════
# ① 3D 자유확산 — ⟨r²⟩ = 6 D t   (해석해. 격리 검증)
# ═══════════════════════════════════════════════════════════════════════
def check_free_diffusion_3d(n_part=800, dt=1e-4, n_steps=20_000):
    print("\n① 3D 자유확산 (상호작용 없음) — ⟨r²⟩ = 6·D·t")
    L = 40.0
    rng = np.random.default_rng(3)
    pos = rng.uniform(-L / 2, L / 2, size=(n_part, 3))
    sim = cpu(11)
    sim.create_state_from_snapshot(frame_3d(pos, L))
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    r0 = unwrapped(sim)
    sim.run(n_steps)
    r1 = unwrapped(sim)
    t = dt * n_steps
    msd = float(((r1 - r0) ** 2).sum(axis=1).mean())
    pred = 6.0 * 1.0 * t                         # D = kT/γ = 1
    rel = msd / pred - 1.0
    sem = float(((r1 - r0) ** 2).sum(axis=1).std(ddof=1) / math.sqrt(n_part)) / pred
    check(abs(rel) < 4 * sem + 0.02,
          "3D 자유확산이 6Dt 와 일치",
          f"측정 {msd:.5f} / 예측 {pred:.5f} = {1+rel:.4f}  ({rel*100:+.2f}%, SEM {sem*100:.2f}%)")

    # 성분별 등방성 — 3D 로 넘어오며 축 하나를 빠뜨리는 실수를 잡는다
    per_axis = ((r1 - r0) ** 2).mean(axis=0)
    iso = per_axis / (2.0 * t)
    check(np.all(np.abs(iso - 1) < 0.12),
          "세 축이 각각 ⟨Δx²⟩=2Dt (축 누락 없음)",
          f"x/y/z = {iso[0]:.3f} / {iso[1]:.3f} / {iso[2]:.3f}")
    return msd, pred


# ═══════════════════════════════════════════════════════════════════════
# ② BoxResize 가 좌표를 아핀 스케일 하는가
# ═══════════════════════════════════════════════════════════════════════
def check_boxresize_affine():
    print("\n② BoxResize — 입자 좌표를 아핀 스케일 하는가 (문서 주장 실측)")
    L0, s = 10.0, 0.5
    pos = np.array([[1.0, 2.0, 3.0], [-4.0, 0.5, -2.5], [0.0, 0.0, 0.0]])
    sim = cpu(12)
    sim.create_state_from_snapshot(frame_3d(pos, L0))
    # 적분기 없이 updater 만 — 좌표 변화의 원인을 BoxResize 하나로 격리한다
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-12, methods=[bd], forces=[])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    box0 = hoomd.Box(Lx=L0, Ly=L0, Lz=L0)
    box1 = hoomd.Box(Lx=L0 * s, Ly=L0 * s, Lz=L0 * s)
    var = hoomd.variant.box.Interpolate(
        initial_box=box0, final_box=box1,
        variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=0, t_ramp=10))
    sim.operations.updaters.append(
        hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(1), box=var))
    sim.run(11)

    got = np.array(sim.state.get_snapshot().particles.position, copy=True)
    box = sim.state.get_snapshot().configuration.box
    check(abs(box[0] - L0 * s) < 1e-9, "박스가 목표 크기로 줄었다",
          f"Lx {box[0]:.6f} (목표 {L0*s:.6f})")
    err = np.abs(got - pos * s).max()
    check(err < 1e-9, "★ 좌표가 정확히 아핀 스케일된다 (r → s·r)",
          f"최대 오차 {err:.3e} — 즉 **결합길이도 같이 줄어든다**")
    return err


# ═══════════════════════════════════════════════════════════════════════
# DLVO 표 + WCA 코어 (chain-bend-2d-dlvo 와 동일 관례)
# ═══════════════════════════════════════════════════════════════════════
def dlvo_forces(nlist, P, extra_types=()):
    r_min = 1.0 + 1e-4
    r_cut = 1.0 + CUTOFF_H_STAR
    r, U, F = build_table_arrays(P, r_min, r_cut)
    tab = md.pair.Table(nlist=nlist, default_r_cut=r_cut)
    wca = md.pair.LJ(nlist=nlist, default_r_cut=SIGMA_CORE_STAR * R_WCA, mode="shift")
    for a in ("A",) + tuple(extra_types):
        for b in ("A",) + tuple(extra_types):
            tab.params[(a, b)] = dict(r_min=r_min, U=U, F=F)
            wca.params[(a, b)] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
    return [tab, wca], r_cut


# ═══════════════════════════════════════════════════════════════════════
# ③ BoxResize + Brownian + 페어힘이 함께 도는가
# ═══════════════════════════════════════════════════════════════════════
def check_boxresize_with_pair(P):
    print("\n③ BoxResize + Brownian + pair.Table + WCA (셀리스트·유한성)")
    n_side, L0 = 6, 18.0
    a = L0 / n_side
    pos = np.array([[(i + .5) * a - L0 / 2, (j + .5) * a - L0 / 2, (k + .5) * a - L0 / 2]
                    for i in range(n_side) for j in range(n_side) for k in range(n_side)])
    sim = cpu(13)
    sim.create_state_from_snapshot(frame_3d(pos, L0))
    nl = md.nlist.Cell(buffer=0.2)
    forces, r_cut = dlvo_forces(nl, P)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-7, methods=[bd], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)

    L1 = 12.0
    var = hoomd.variant.box.Interpolate(
        initial_box=hoomd.Box(Lx=L0, Ly=L0, Lz=L0),
        final_box=hoomd.Box(Lx=L1, Ly=L1, Lz=L1),
        variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=0, t_ramp=20_000))
    sim.operations.updaters.append(
        hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(100), box=var))
    try:
        sim.run(20_100)
        ran = True
    except Exception as e:                                     # noqa: BLE001
        ran = False
        print(f"      크래시: {type(e).__name__}: {e}")
    check(ran, "압축 중 크래시 없음")
    if not ran:
        return None
    pe = thermo.potential_energy
    box = sim.state.get_snapshot().configuration.box
    check(pe is not None and math.isfinite(pe), "퍼텐셜 에너지가 유한",
          f"PE = {pe:.4f} kT (N={len(pos)})")
    check(abs(box[0] - L1) < 1e-6, "압축이 목표에서 정확히 멈춘다",
          f"Lx {box[0]:.6f} (목표 {L1})")
    check(r_cut < box[0] / 2, "압축 후에도 r_cut < L/2 (함정 6)",
          f"r_cut {r_cut:.4f} < L/2 {box[0]/2:.4f}")
    return pe


# ═══════════════════════════════════════════════════════════════════════
# ④ ★ 압축이 DLVO 결합을 부수는 문턱
# ═══════════════════════════════════════════════════════════════════════
def check_crush_threshold(P, W, quick=False):
    print("\n④ ★ 압축이 DLVO 2차극소 결합을 부수는 문턱")
    h_min, h_bar = W["h_min"], W["barrier_h"]
    ell = 1.0 + h_min
    eps_crit = (h_min - h_bar) / ell
    tau_bond = 1.0 / W["k_bond_star"]
    print(f"   원장: h_min*={h_min:.6f}  barrier*={h_bar:.6f}  ℓ*={ell:.6f}")
    print(f"         k_bond*={W['k_bond_star']:.4g} kT/d²   τ_bond*={tau_bond:.4g} τ_B")
    print(f"   예측 문턱 ε_crit = (h_min*−barrier*)/ℓ* = {eps_crit:.6f}  ({eps_crit*100:.3f}%/트리거)")

    dt = 1e-8                       # dt/τ_bond ≈ 0.0104
    T = 2000                        # 트리거 간격 → 이완시간 T·dt = 21 τ_bond
    L0 = 6.0
    s_tot = 0.95                    # 총 5% 선형 압축
    eps_list = [0.002, 0.004, 0.006, 0.008, 0.012, 0.020]
    if quick:
        eps_list = [0.004, 0.008, 0.020]

    print(f"   dt={dt:g} (dt/τ_bond={dt/tau_bond:.4f}) · 트리거 간격 {T} 스텝"
          f" (= {T*dt/tau_bond:.0f} τ_bond 이완)")
    print("   ε/트리거     최종 h*        판정")
    rows = []
    for eps in eps_list:
        t_ramp = max(int(round((1.0 - s_tot) * T / eps)), T)
        sim = cpu(14)
        sim.create_state_from_snapshot(
            frame_3d([[-ell / 2, 0, 0], [ell / 2, 0, 0]], L0))
        nl = md.nlist.Cell(buffer=0.2)
        forces, _ = dlvo_forces(nl, P)
        # kT=0 결정론적 — 잡음이 문턱을 흐리지 않게 (bd-physics: 적분기 가정 시험은 kT=0)
        bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
        integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        var = hoomd.variant.box.Interpolate(
            initial_box=hoomd.Box(Lx=L0, Ly=L0, Lz=L0),
            final_box=hoomd.Box(Lx=L0 * s_tot, Ly=L0 * s_tot, Lz=L0 * s_tot),
            variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=0, t_ramp=t_ramp))
        sim.operations.updaters.append(
            hoomd.update.BoxResize(trigger=hoomd.trigger.Periodic(T), box=var))
        sim.run(t_ramp + T)
        p = np.array(sim.state.get_snapshot().particles.position, copy=True)
        sep = float(np.linalg.norm(p[1] - p[0]))
        h_end = sep - 1.0
        survived = h_end > h_bar * 3          # 장벽 훨씬 밖 = 2차극소에 남음
        rows.append((eps, h_end, survived))
        print(f"   {eps*100:6.2f}%   {h_end:12.6f}   "
              f"{'2차극소 유지' if survived else '★ 붕괴(1차극소/접촉)'}")

    ok = [e for e, _, s in rows if s]
    bad = [e for e, _, s in rows if not s]
    check(bool(ok) and bool(bad), "문턱이 스윕 범위 안에 있다 (양쪽 다 관측됨)",
          f"유지 ≤{max(ok)*100:.2f}% · 붕괴 ≥{min(bad)*100:.2f}%" if ok and bad else "")
    if ok and bad:
        lo, hi = max(ok), min(bad)
        inside = lo <= eps_crit <= hi or abs(eps_crit - lo) / eps_crit < 0.5
        check(inside, "★ 실측 문턱이 해석 예측 ε_crit 과 정합",
              f"실측 {lo*100:.2f}~{hi*100:.2f}% · 예측 {eps_crit*100:.3f}%")
    return eps_crit, rows


# ═══════════════════════════════════════════════════════════════════════
# ⑤ 비용 — 겔화 구성의 steps/s → 벽시계
# ═══════════════════════════════════════════════════════════════════════
def check_cost(P, W, n_part=1528, phi0=0.02, phi1=0.10, bench_steps=3000):
    print("\n⑤ 비용 — 겔화 구성 steps/s 와 벽시계 (N3: 판정 전에 재라)")
    tau_bond = 1.0 / W["k_bond_star"]
    dt = 1e-2 * tau_bond                          # 설계 관례 dt/τ_fast = 1e-2
    out = {}
    for tag, phi in (("초기(묽음)", phi0), ("최종(압축후)", phi1)):
        L = (n_part * math.pi / (6.0 * phi)) ** (1 / 3)
        n_side = int(math.ceil(n_part ** (1 / 3)))
        a = L / n_side
        pos = np.array([[(i + .5) * a - L / 2, (j + .5) * a - L / 2, (k + .5) * a - L / 2]
                        for i in range(n_side) for j in range(n_side)
                        for k in range(n_side)])[:n_part]
        sim = cpu(15)
        sim.create_state_from_snapshot(frame_3d(pos, L))
        nl = md.nlist.Cell(buffer=0.2)
        forces, _ = dlvo_forces(nl, P)
        bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
        integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        sim.run(200)                                          # 워밍업(셀리스트·오토튠)
        t0 = time.perf_counter()
        sim.run(bench_steps)
        rate = bench_steps / (time.perf_counter() - t0)
        out[tag] = (L, rate)
        print(f"   {tag:12s} φ={phi:.3f}  L/d={L:6.2f}  {rate:9.0f} steps/s")

    rate_min = min(r for _, r in out.values())
    print(f"\n   dt = 1e-2·τ_bond = {dt:.4g} τ_B   (τ_bond*={tau_bond:.4g})")
    for tau_target in (1.0, 5.0, 10.0):
        steps = tau_target / dt
        hours = steps / rate_min / 3600
        print(f"   겔화 {tau_target:4.1f} τ_B  →  {steps:.3g} 스텝  →  "
              f"{hours:8.1f} 시간/시드 (최저 steps/s 기준)")
    check(True, "비용 측정 완료 (판정은 위 표로)", f"최저 {rate_min:.0f} steps/s")
    return out, dt


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="④ 스윕 점 수를 줄인다")
    args = ap.parse_args()

    print("=" * 78)
    print("network 능력조사 — 3D BD + BoxResize   (규칙 4: 감으로 쓰지 않는다)")
    print("=" * 78)
    print(f"hoomd {hoomd.version.version} · CPU · 단위 d=kT=γ=1 (τ_B=1)")

    sys_ = load_system(ROOT / "intake/chain-bend-2d-dlvo/system.yaml")
    P = dlvo_reduced_params(sys_)
    W = find_well(P)
    print(f"DLVO 원장 승계: 장벽 {W['barrier_U']:.2f} kT @ h*={W['barrier_h']:.5f} · "
          f"2차극소 {W['U_min']:.3f} kT @ h*={W['h_min']:.5f}")

    check_free_diffusion_3d()
    check_boxresize_affine()
    check_boxresize_with_pair(P)
    check_crush_threshold(P, W, quick=args.quick)
    check_cost(P, W)

    print("\n" + "=" * 78)
    print(f"{'✓ PASS' if not FAIL else '✗ FAIL'} — {len(PASS)}/{len(PASS)+len(FAIL)} 정상")
    for f in FAIL:
        print(f"   실패: {f}")
    print("=" * 78)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
