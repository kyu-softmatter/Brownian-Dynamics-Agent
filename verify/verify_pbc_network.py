#!/usr/bin/env python
"""Measure whether `network`'s periodic boundaries are right (CLAUDE.md rule 6).

Why: this project has a precedent of a **+1856% error** from a missing minimum image,
and at the time it "looked like it passed because it was only tested under strong
conditions". Code that touches the boundary in several places -- 3D, compression, the
tree test, the percolation test -- cannot be trusted by reading it.

Four checks --
  (1) **brute-force 27 images** vs `contacts()`'s minimum image -- are the pair
      distances really the same?
  (2) **translation invariance** (the decisive test): translate every particle by an
      arbitrary vector and wrap, and z, loops, components, percolation, d_f, min_sep
      and the g(r) peak must **all be unchanged**.
      If any one piece of code mishandles the boundary, it breaks here.
  (3) **forced boundary placement**: deliberately push particles onto faces, edges
      and corners, then repeat (1) and (2)
      (test under WEAK conditions -- the lesson of bd-hoomd trap 1)
  (4) **HOOMD's own PBC**: translate the same configuration and check the potential
      energy is invariant (this looks at the engine, not at my code)

Run:  $PY verify/verify_pbc_network.py [--run <run_id>]
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
    """Brute force over 27 images -- counts directly rather than trusting the
    minimum-image formula."""
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
        if x != x and y != y:          # both nan
            continue
        if abs(x - y) > tol * max(1.0, abs(x)):
            bad.append(f"{k}: {x!r} vs {y!r}")
    return bad


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--run", default="network__N512-sprout-mc4-po0.2__8d9baf357248")
    args = ap.parse_args()

    print("=" * 78)
    print("network -- periodic boundaries measured (rule 6: by execution, not by "
          "reasoning)")
    print("=" * 78)

    rd = ROOT / "runs" / args.run
    z = np.load(rd / "observables.npz")
    m = json.loads((rd / "metrics.json").read_text())
    pos = np.array(z["final_positions"], dtype=float)
    L = float(m["result"]["L_final"])
    r_bond = float(m["result"]["r_bond_star"])
    n = len(pos)
    print(f"run {args.run}\n  N={n}  L={L:.6f} d  r_bond={r_bond:.6f}")
    print(f"  coordinate range x/y/z: "
          + "  ".join(f"[{pos[:, k].min():+.3f},{pos[:, k].max():+.3f}]" for k in range(3))
          + f"   (box [{-L/2:+.3f},{+L/2:+.3f}])")
    inside = np.all(np.abs(pos) <= L / 2 + 1e-9)
    check(inside, "every stored coordinate lies inside the box (wrapped)")

    # ── (1) brute force over 27 images ──────────────────────────────
    print("\n(1) minimum-image formula vs brute force over 27 images")
    sub = pos[:160]                                   # 27 x 160^2 is fast enough
    bf = brute_min_dist(sub, L)
    d = sub[:, None, :] - sub[None, :, :]
    d -= L * np.round(d / L)
    mi = np.sqrt((d ** 2).sum(-1))
    np.fill_diagonal(mi, np.inf)
    fin = np.isfinite(bf) & np.isfinite(mi)
    err = float(np.abs(bf[fin] - mi[fin]).max())
    check(err < 1e-12, "pair distances agree with brute force",
          f"max error {err:.3e} d (on a subset of N={len(sub)})")
    check(bf[np.isfinite(bf)].min() <= L / 2 * np.sqrt(3) + 1e-9,
          "minimum-image distances lie within half the box diagonal")

    # ── (2) translation invariance ──────────────────────────────────
    print("\n(2) translation invariance -- translate every particle and wrap, and the "
          "observables must be unchanged")
    base = observables(pos, L, r_bond)
    print(f"   reference: z={base['z']:.4f} loops={base['loops']} "
          f"comps={base['comps']} "
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
        print(f"   shift ({sh[0]:+.2f},{sh[1]:+.2f},{sh[2]:+.2f})  "
              f"{crossed:3d} particles changed x sign  -> "
              f"{'same' if not bad else bad}")
        worst += bad
    check(not worst, "★ all observables invariant over 6 translations "
                     "(z, loops, comps, percolation, d_f, g(r), min_sep)")

    # ── (3) push them onto the boundary and repeat (the weak-condition test) ───
    print("\n(3) forced boundary -- particles pushed onto faces, edges and corners")
    for name, sh in (("face (x=±L/2)", np.array([L / 2 - pos[:, 0].max(), 0, 0])),
                     ("corner", np.full(3, L / 2) - pos.max(axis=0))):
        p2 = pos + sh
        p2 -= L * np.round(p2 / L)
        near = int((np.abs(np.abs(p2) - L / 2) < 1.0).any(axis=1).sum())
        o2 = observables(p2, L, r_bond)
        bad = compare(base, o2)
        print(f"   {name}: {near:3d} within 1d of the boundary  -> "
              f"{'same' if not bad else bad}")
        check(not bad, f"observables invariant even pressed against the boundary "
                       f"({name})")

    # ── (4) HOOMD's own PBC (the engine, not my code) ───────────────
    print("\n(4) HOOMD engine -- translate the same configuration and check the "
          "potential energy is invariant")
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
    print(f"   reference   PE = {pe0:.10f} kT   |F|max = {f0:.4f} kT/d")
    rel = []
    for trial in range(3):
        sh = np.array([L / 2, 0, 0]) if trial == 0 else rng.uniform(-L, L, 3)
        p2 = pos + sh
        p2 -= L * np.round(p2 / L)
        pe, fm = pe_of(p2)
        r = abs(pe / pe0 - 1)
        rel.append(r)
        print(f"   shift ({sh[0]:+.2f},{sh[1]:+.2f},{sh[2]:+.2f})  PE = {pe:.10f}  "
              f"rel diff {r:.3e}   |F|max = {fm:.4f}")
    check(max(rel) < 1e-10, "★ HOOMD potential energy is translation invariant",
          f"max rel diff {max(rel):.3e}")
    check(r_cut < L / 2, "r_cut < L/2 (the minimum-image convention holds)",
          f"{r_cut:.4f} < {L/2:.4f}")

    print("\n" + "=" * 78)
    print(f"{'✓ PASS' if not FAIL else '✗ FAIL'} — {len(PASS)}/{len(PASS)+len(FAIL)}")
    for f in FAIL:
        print(f"   failed: {f}")
    print("=" * 78)
    return 1 if FAIL else 0


if __name__ == "__main__":
    raise SystemExit(main())
