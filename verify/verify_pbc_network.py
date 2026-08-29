#!/usr/bin/env python
"""`network` 의 주기경계(PBC)가 맞는지 실측한다 (CLAUDE.md 규칙 6).

왜: 이 프로젝트는 최소이미지 누락으로 **+1856% 오차**를 낸 전례가 있고, 그때도
"강한 조건에서만 테스트해서 통과한 것처럼 보였다". 3D·압축·트리 판정·침투 판정처럼
경계를 여러 곳에서 만지는 코드는 눈으로 읽어서 확신할 수 없다.

검사 4종 —
  ① **브루트포스 27이미지** vs `contacts()` 의 최소이미지 — 쌍거리가 정말 같은가
  ② **병진 불변성** (결정적 시험): 전 입자를 임의 벡터로 옮기고 감싸면
     z·고리·성분·침투·d_f·min_sep·g(r) 피크가 **전부 같아야** 한다.
     경계를 잘못 다루는 코드가 하나라도 있으면 여기서 깨진다.
  ③ **경계 강제 배치**: 입자를 일부러 면·모서리·꼭짓점에 몰아놓고 ①②를 다시
     (약한 조건에서 테스트하라 — bd-hoomd 함정 1의 교훈)
  ④ **HOOMD 자신의 PBC**: 같은 배치를 병진시켜 퍼텐셜 에너지가 불변인가
     (내 코드가 아니라 엔진 쪽을 본다)

실행:  $PY scratch/verify_pbc_network.py [--run <run_id>]
"""
from __future__ import annotations

import argparse
import itertools
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from network_3d import (                                            # noqa: E402
    SIGMA_CORE_STAR, build_table_arrays, contacts, dlvo_reduced_params,
    fractal_dimension, load_system, percolates, rdf, topology,
)

R_WCA = 2 ** (1 / 6)
PASS, FAIL = [], []


def check(ok, label, detail=""):
    (PASS if ok else FAIL).append(label)
    print(f"  {'✓' if ok else '✗'} {label}" + (f"   {detail}" if detail else ""))


def brute_min_dist(pos, L):
    """27 이미지 브루트포스 — 최소이미지 공식을 믿지 않고 직접 센다."""
    shifts = np.array(list(itertools.product((-L, 0.0, L), repeat=3)))
    n = len(pos)
    out = np.full((n, n), np.inf)
    for s in shifts:
        d = pos[:, None, :] - (pos[None, :, :] + s)
        r = np.sqrt((d ** 2).sum(-1))
        np.fill_diagonal(r, np.inf)
        out = np.minimum(out, r)
    return out


def observables(pos, L, r_bond, seed=7):
    pairs, hp, r_all = contacts(pos, L, r_bond)
    t = topology(len(pos), pairs)
    rng = np.random.default_rng(seed)
    mid, g = rdf(r_all, len(pos), L, min(L / 2, 4.0))
    return dict(z=t["z"], loops=t["loops"], comps=t["n_components"],
                dangling=t["dangling"], largest=t["largest_cluster"],
                perc=percolates(pos, L, pairs),
                d_f=fractal_dimension(pos, L, rng),
                min_sep=float(r_all.min()), n_pairs=len(pairs),
                rdf_peak=float(mid[int(np.argmax(g))]),
                g_sum=float(np.nansum(g)))


def compare(a, b, tol=1e-9):
    bad = []
    for k in a:
        x, y = a[k], b[k]
        if x != x and y != y:          # 둘 다 nan
            continue
        if abs(x - y) > tol * max(1.0, abs(x)):
            bad.append(f"{k}: {x!r} vs {y!r}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="network__N512-sprout-mc4-po0.2__8d9baf357248")
    args = ap.parse_args()

    print("=" * 78)
    print("network — 주기경계 실측 (규칙 6: 추론이 아니라 실행으로)")
    print("=" * 78)

    rd = ROOT / "runs" / args.run
    z = np.load(rd / "observables.npz")
    m = json.loads((rd / "metrics.json").read_text())
    pos = np.array(z["final_positions"], dtype=float)
    L = float(m["result"]["L_final"])
    r_bond = float(m["result"]["r_bond_star"])
    n = len(pos)
    print(f"런 {args.run}\n  N={n}  L={L:.6f} d  r_bond={r_bond:.6f}")
    print(f"  좌표 범위 x/y/z: "
          + "  ".join(f"[{pos[:, k].min():+.3f},{pos[:, k].max():+.3f}]" for k in range(3))
          + f"   (박스 [{-L/2:+.3f},{+L/2:+.3f}])")
    inside = np.all(np.abs(pos) <= L / 2 + 1e-9)
    check(inside, "저장된 좌표가 전부 박스 안에 있다 (감싸짐)")

    # ── ① 브루트포스 27이미지 ────────────────────────────────────────
    print("\n① 최소이미지 공식 vs 27이미지 브루트포스")
    sub = pos[:160]                                   # 27×160² 이면 충분히 빠르다
    bf = brute_min_dist(sub, L)
    d = sub[:, None, :] - sub[None, :, :]
    d -= L * np.round(d / L)
    mi = np.sqrt((d ** 2).sum(-1))
    np.fill_diagonal(mi, np.inf)
    fin = np.isfinite(bf) & np.isfinite(mi)
    err = float(np.abs(bf[fin] - mi[fin]).max())
    check(err < 1e-12, "쌍거리가 브루트포스와 일치", f"최대 오차 {err:.3e} d (N={len(sub)} 부분집합)")
    check(bf[np.isfinite(bf)].min() <= L / 2 * np.sqrt(3) + 1e-9,
          "최소이미지 거리가 박스 대각선 반경 안")

    # ── ② 병진 불변성 ────────────────────────────────────────────────
    print("\n② 병진 불변성 — 전 입자를 옮기고 감싸면 관측량이 같아야 한다")
    base = observables(pos, L, r_bond)
    print(f"   기준: z={base['z']:.4f} loops={base['loops']} comps={base['comps']} "
          f"perc={base['perc']:.4f} d_f={base['d_f']:.4f} min_sep={base['min_sep']:.6f} "
          f"pairs={base['n_pairs']}")
    rng = np.random.default_rng(3)
    worst = []
    for trial in range(6):
        sh = rng.uniform(-L, L, 3) if trial else np.array([L / 2, 0.0, 0.0])
        p2 = pos + sh
        p2 -= L * np.round(p2 / L)
        o2 = observables(p2, L, r_bond)
        bad = compare(base, o2)
        crossed = int((np.sign(pos[:, 0]) != np.sign(p2[:, 0])).sum())
        print(f"   이동 ({sh[0]:+.2f},{sh[1]:+.2f},{sh[2]:+.2f})  "
              f"x부호 바뀐 입자 {crossed:3d}개  → {'같음' if not bad else bad}")
        worst += bad
    check(not worst, "★ 6회 병진에서 관측량 전부 불변 (z·고리·성분·침투·d_f·g(r)·min_sep)")

    # ── ③ 경계에 몰아놓고 다시 (약한 조건 테스트) ────────────────────
    print("\n③ 경계 강제 — 입자를 면·모서리·꼭짓점에 붙여놓고")
    for name, sh in (("면(x=±L/2)", np.array([L / 2 - pos[:, 0].max(), 0, 0])),
                     ("꼭짓점", np.full(3, L / 2) - pos.max(axis=0))):
        p2 = pos + sh
        p2 -= L * np.round(p2 / L)
        near = int((np.abs(np.abs(p2) - L / 2) < 1.0).any(axis=1).sum())
        o2 = observables(p2, L, r_bond)
        bad = compare(base, o2)
        print(f"   {name}: 경계 1d 안에 {near:3d}개  → {'같음' if not bad else bad}")
        check(not bad, f"경계 밀착({name})에서도 관측량 불변")

    # ── ④ HOOMD 자신의 PBC (내 코드가 아니라 엔진) ───────────────────
    print("\n④ HOOMD 엔진 — 같은 배치를 병진시켜 퍼텐셜 에너지가 불변인가")
    import gsd.hoomd
    import hoomd
    import hoomd.md as md
    sys_ = load_system(ROOT / "intake/network/system.yaml")
    P = dlvo_reduced_params(sys_)
    r_min, r_cut = 1.0 + 1e-4, 1.06
    _, U_arr, F_arr = build_table_arrays(P, r_min, r_cut)

    def pe_of(p):
        f = gsd.hoomd.Frame()
        f.particles.N = len(p)
        f.particles.position = p
        f.particles.typeid = [0] * len(p)
        f.particles.types = ["A"]
        f.particles.mass = [1.0] * len(p)
        f.configuration.box = [L, L, L, 0, 0, 0]
        f.configuration.dimensions = 3
        sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=2)
        sim.create_state_from_snapshot(f)
        cell = md.nlist.Cell(buffer=0.2)
        tab = md.pair.Table(nlist=cell, default_r_cut=r_cut)
        tab.params[("A", "A")] = dict(r_min=r_min, U=U_arr, F=F_arr)
        wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * R_WCA, mode="shift")
        wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
        bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
        integ = md.Integrator(dt=1e-16, methods=[bd], forces=[tab, wca])
        integ.integrate_rotational_dof = False
        sim.operations.integrator = integ
        sim.run(0)
        return (float(np.array(tab.energies).sum() + np.array(wca.energies).sum()),
                float(np.abs(np.array(tab.forces) + np.array(wca.forces)).max()))

    pe0, f0 = pe_of(pos)
    print(f"   기준        PE = {pe0:.10f} kT   |F|max = {f0:.4f} kT/d")
    rel = []
    for trial in range(3):
        sh = np.array([L / 2, 0, 0]) if trial == 0 else rng.uniform(-L, L, 3)
        p2 = pos + sh
        p2 -= L * np.round(p2 / L)
        pe, fm = pe_of(p2)
        r = abs(pe / pe0 - 1)
        rel.append(r)
        print(f"   이동 ({sh[0]:+.2f},{sh[1]:+.2f},{sh[2]:+.2f})  PE = {pe:.10f}  "
              f"상대차 {r:.3e}   |F|max = {fm:.4f}")
    check(max(rel) < 1e-10, "★ HOOMD 퍼텐셜 에너지가 병진 불변",
          f"최대 상대차 {max(rel):.3e}")
    check(r_cut < L / 2, "r_cut < L/2 (최소이미지 규약 성립)",
          f"{r_cut:.4f} < {L/2:.4f}")

    print("\n" + "=" * 78)
    print(f"{'✓ PASS' if not FAIL else '✗ FAIL'} — {len(PASS)}/{len(PASS)+len(FAIL)}")
    for f in FAIL:
        print(f"   실패: {f}")
    print("=" * 78)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
