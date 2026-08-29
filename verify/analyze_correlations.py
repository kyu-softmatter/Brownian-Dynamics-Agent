"""입자 위치 상관 — 접선상관(지속길이)과 변위 상관행렬. 기존 GSD 궤적에서 사후 계산.

★ 왜 이 두 가지인가
  ① 접선-접선 상관  ⟨t_i·t_{i+m}⟩ = exp(−m·ℓ/ℓ_p)
     반굴곡 고분자의 **표준** 진단량이고, 이 계에는 **해석적 예측이 있다**:
       2D 이산 사슬 U = ½κ_θ Σθ²  →  ⟨θ²⟩ = kT/κ_θ  →  ℓ_p = 2κ_θℓ/kT
       · JKR  (κ_θ* = 1.391e6 kT) → ℓ_p/ℓ = 2.78e6  (사실상 강체)
       · DLVO (굽힘항 없음)        → ℓ_p/ℓ ≈ 1      (자유연결, 즉시 상관 소멸)
     즉 **두 계가 6자릿수 다른 예측**을 내는 양이다. K'(216배)보다 훨씬 큰 대비다.
  ② 변위 상관행렬 C_ij = ⟨δy_i δy_j⟩ / √(⟨δy_i²⟩⟨δy_j²⟩)
     빔이면 굽힘이 먼 비드끼리 묶어 C_ij 가 넓게 양(+)이고, 힌지면 국소적으로 죽는다.

⚠️ 트랩 보정: 비드 0·mid·n−1 이 트랩에 묶여 있어 상관을 강하게 왜곡한다.
   그래서 **트랩 사이의 자유 구간만** 따로 재고, 전체 사슬 값도 같이 낸다.

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
    """⟨t_i·t_{i+m}⟩ vs m. lo:hi 로 자유 구간만 자를 수 있다."""
    p = xy[:, lo:hi, :] if (lo is not None or hi is not None) else xy
    b = np.diff(p, axis=1)                                   # 결합 벡터 (T, nb, 2)
    b /= np.linalg.norm(b, axis=-1, keepdims=True)
    nb = b.shape[1]
    out = []
    for m in range(nb):
        d = np.einsum("tij,tij->ti", b[:, :nb - m], b[:, m:])   # 내적
        out.append(d.mean())
    return np.array(out)


def persistence_length(c, ell=1.0):
    """log⟨t·t⟩ 의 초기 기울기에서 ℓ_p. 상관이 안 죽으면 하한만 준다."""
    m = np.arange(len(c))
    ok = (c > 0.05) & (m > 0)
    if ok.sum() < 2:
        return np.nan, "감쇠 없음(하한만)"
    s = np.polyfit(m[ok], np.log(c[ok]), 1)[0]
    if s >= -1e-12:
        return np.inf, "감쇠 없음"
    return -ell / s, "적합"


def disp_corr(xy, lo=None, hi=None):
    y = xy[:, lo:hi, 1] if (lo is not None or hi is not None) else xy[:, :, 1]
    dy = y - y.mean(axis=0, keepdims=True)
    C = np.corrcoef(dy.T)
    return C


def analyze(pat, label, n, free=None):
    ds = sorted(glob.glob(pat))
    if not ds:
        print(f"  [{label}] 런 없음"); return None
    cs, Cs = [], []
    for d in ds:
        xy = load_xy(d, n)
        cs.append(tangent_corr(xy, *(free or (None, None))))
        Cs.append(disp_corr(xy, *(free or (None, None))))
    c = np.mean(cs, axis=0); C = np.mean(Cs, axis=0)
    lp, how = persistence_length(c)
    seg = f" [자유구간 {free[0]}:{free[1]}]" if free else " [전체 사슬]"
    print(f"\n  {label}{seg}   런 {len(ds)}개")
    print(f"    ⟨t·t⟩(m) = " + "  ".join(f"{v:+.3f}" for v in c[:min(7, len(c))]))
    print(f"    → ℓ_p/ℓ = {lp:.3g}  ({how})")
    off = C[np.triu_indices_from(C, k=1)]
    print(f"    변위상관 비대각 평균 = {off.mean():+.3f}  (최소 {off.min():+.3f} 최대 {off.max():+.3f})")
    return dict(c=c, C=C, lp=lp, label=label, n=n, free=free)


print("=" * 84)
print("입자 위치 상관 — DLVO-only vs DLVO+JKR  (n=9, ω=3000, a=1470nm=1d)")
print("=" * 84)
print("해석적 예측:  JKR ℓ_p/ℓ = 2κ_θ*/1 = 2.78e6 (사실상 강체)   DLVO ℓ_p/ℓ ≈ 1 (자유연결)")
res = {}
res["dlvo_full"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470__*", "DLVO-only", 9)
res["jkr_full"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr__*", "DLVO+JKR", 9)
# n=9 는 트랩이 0·4·8 → 자유 비드는 1,2,3 과 5,6,7. 한쪽 구간(0~4)만 봐도 결합 4개뿐이라 짧다
res["dlvo_seg"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470__*", "DLVO-only", 9, free=(0, 5))
res["jkr_seg"] = analyze("runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr__*", "DLVO+JKR", 9, free=(0, 5))

print()
print("=" * 84)
print("긴 사슬 — n=25 (트랩 0·12·24, 자유 구간 11비드), ω=10 rad/s, DLVO-only")
print("=" * 84)
res["n25_full"] = analyze("runs/chain-bend-2d-dlvo__n25-w10-a632__*", "DLVO-only n=25", 25)
res["n25_seg"] = analyze("runs/chain-bend-2d-dlvo__n25-w10-a632__*", "DLVO-only n=25", 25, free=(0, 13))

np.save("/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/"
        "7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/corr_res.npy",
        res, allow_pickle=True)
print("\n저장: corr_res.npy")
