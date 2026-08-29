"""Particle position correlations -- tangent correlation (persistence length) and the
displacement correlation matrix. Computed after the fact from existing GSD
trajectories.

★ Why these two
  (1) tangent-tangent correlation  <t_i.t_{i+m}> = exp(-m*l/l_p)
      The **standard** diagnostic for a semiflexible polymer, and this system has an
      **analytic prediction**:
        2D discrete chain U = 0.5*kappa_theta*sum(theta^2)
        -> <theta^2> = kT/kappa_theta  ->  l_p = 2*kappa_theta*l/kT
        . JKR  (kappa_theta* = 1.391e6 kT) -> l_p/l = 2.78e6  (effectively rigid)
        . DLVO (no bending term)           -> l_p/l ~ 1  (freely jointed, correlation
                                                          dies immediately)
      So the two systems give predictions **six orders of magnitude apart** -- a far
      bigger contrast than K' (216x).
  (2) displacement correlation matrix
      C_ij = <dy_i dy_j> / sqrt(<dy_i^2><dy_j^2>)
      For a beam, bending ties distant beads together so C_ij is broadly positive;
      for a hinge it dies off locally.

⚠️ Trap correction: beads 0, mid and n-1 are tied to traps, which distorts the
   correlations strongly. So the **free segment between the traps** is measured
   separately, alongside the whole-chain value.

    $PY scratch/analyze_correlations.py
"""
import glob
import json
from pathlib import Path

import gsd.hoomd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent


def load_xy(d, n):
    t = gsd.hoomd.open(str(Path(d) / "traj_A.gsd"), "r")
    return np.array([t[i].particles.position[:n, :2] for i in range(len(t))])


def tangent_corr(xy, lo=None, hi=None):
    """<t_i.t_{i+m}> vs m. lo:hi restricts it to the free segment."""
    p = xy[:, lo:hi, :] if (lo is not None or hi is not None) else xy
    b = np.diff(p, axis=1)                                   # bond vectors (T, nb, 2)
    b /= np.linalg.norm(b, axis=-1, keepdims=True)
    nb = b.shape[1]
    out = []
    for m in range(nb):
        d = np.einsum("tij,tij->ti", b[:, :nb - m], b[:, m:])   # dot product
        out.append(d.mean())
    return np.array(out)


def persistence_length(c, ell=1.0):
    """l_p from the initial slope of log<t.t>.

    If the correlation does not decay, only a lower bound is returned.
    """
    m = np.arange(len(c))
    ok = (c > 0.05) & (m > 0)
    if ok.sum() < 2:
        return np.nan, "no decay (lower bound only)"
    s = np.polyfit(m[ok], np.log(c[ok]), 1)[0]
    if s >= -1e-12:
        return np.inf, "no decay"
    return -ell / s, "fitted"


def disp_corr(xy, lo=None, hi=None):
    y = xy[:, lo:hi, 1] if (lo is not None or hi is not None) else xy[:, :, 1]
    dy = y - y.mean(axis=0, keepdims=True)
    C = np.corrcoef(dy.T)
    return C


def analyze(pat, label, n, free=None):
    ds = sorted(glob.glob(pat))
    if not ds:
        print(f"  [{label}] no runs"); return None
    cs, Cs = [], []
    for d in ds:
        xy = load_xy(d, n)
        cs.append(tangent_corr(xy, *(free or (None, None))))
        Cs.append(disp_corr(xy, *(free or (None, None))))
    c = np.mean(cs, axis=0); C = np.mean(Cs, axis=0)
    lp, how = persistence_length(c)
    seg = f" [free segment {free[0]}:{free[1]}]" if free else " [whole chain]"
    print(f"\n  {label}{seg}   {len(ds)} runs")
    print(f"    ⟨t·t⟩(m) = " + "  ".join(f"{v:+.3f}" for v in c[:min(7, len(c))]))
    print(f"    → ℓ_p/ℓ = {lp:.3g}  ({how})")
    off = C[np.triu_indices_from(C, k=1)]
    print(f"    displacement-correlation off-diagonal mean = {off.mean():+.3f}  "
          f"(min {off.min():+.3f} max {off.max():+.3f})")
    return dict(c=c, C=C, lp=lp, label=label, n=n, free=free)


print("=" * 84)
print("particle position correlations -- DLVO-only vs DLVO+JKR  "
      "(n=9, omega=3000, a=1470nm=1d)")
print("=" * 84)
print("analytic prediction:  JKR l_p/l = 2*kappa_theta*/1 = 2.78e6 (effectively "
      "rigid)   DLVO l_p/l ~ 1 (freely jointed)")
res = {}
res["dlvo_full"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470__*", "DLVO-only", 9)
res["jkr_full"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr__*", "DLVO+JKR", 9)
# For n=9 the traps are at 0, 4 and 8, so the free beads are 1,2,3 and 5,6,7. Taking
# one segment (0..4) alone gives only 4 bonds, which is short
res["dlvo_seg"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470__*", "DLVO-only", 9, free=(0, 5))
res["jkr_seg"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr__*", "DLVO+JKR", 9, free=(0, 5))

print()
print("=" * 84)
print("long chain -- n=25 (traps at 0, 12, 24; free segment of 11 beads), "
      "omega=10 rad/s, DLVO-only")
print("=" * 84)
res["n25_full"] = analyze("runs/chain-bend-2d-dlvo__n25-w10-a632__*", "DLVO-only n=25", 25)
res["n25_seg"] = analyze("runs/chain-bend-2d-dlvo__n25-w10-a632__*", "DLVO-only n=25", 25, free=(0, 13))

np.save("/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/"
        "7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/corr_res.npy",
        res, allow_pickle=True)
print("\nsaved: corr_res.npy")
