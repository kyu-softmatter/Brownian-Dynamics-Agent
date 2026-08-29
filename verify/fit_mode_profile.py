"""JKR 사슬의 비드별 A·sin(ωt+φ) 피팅 → 중심거리별 진폭·위상 프로파일.

각 비드의 y(t) 를 구동 주파수에서 락인해 복소 위상자 ŷ_i 를 얻고, A_i=|ŷ_i|,
φ_i=arg(ŷ_i) 를 중심(구동 비드)으로부터의 거리 s=|i−mid| 에 대해 정리한다.

★ **해석적 예측이 있다** — 선형응답:  (iωγI + A_bend + T) ŷ = k_t·a·e_mid
  이 계는 굽힘행렬 A_bend 와 트랩 T 가 전부 알려져 있으므로 ŷ 를 정확히 풀 수 있다.
  측정과 대조하면 구현·물리가 맞는지 **비드마다** 검증된다 (지금까지는 구동 비드
  하나만 봤다 — |ŷ_mid| 가 예측과 +1.3%).

  위상은 구동 비드 기준으로 잡는다 (φ_mid ≡ 0). 절대 위상은 ZOH·표본 지연에
  오염되지만 **차이는 상쇠**된다 (bd-hoomd 함정 17 과 같은 논리).

    $PY scratch/fit_mode_profile.py
"""
import glob
import json
import sys
from pathlib import Path

import gsd.hoomd
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))
from bdbot import lockin as LI            # noqa: E402
from chain_bend_dlvo_2d import bending_matrix, trapped_indices  # noqa: E402

PAT = "runs/chain-bend-2d-dlvo__n9-w3000-a1470-jkr-kt100__*"


def phasors_from_npz(d):
    """observables.npz 의 전 표본에서 비드별 위상자. GSD 보다 표본이 촘촘하다."""
    z = np.load(Path(d) / "observables.npz", allow_pickle=True)
    keys = set(z.files)
    if "shape_y" not in keys:
        return None, None
    t = np.asarray(z["t"], dtype=float)
    ys = np.asarray(z["shape_y"], dtype=float)      # (T, n)
    return t, ys


def phasors_from_gsd(d, n):
    tr = gsd.hoomd.open(str(Path(d) / "traj_A.gsd"), "r")
    ys = np.array([tr[i].particles.position[:n, 1] for i in range(len(tr))])
    s = json.loads((Path(d) / "spec.json").read_text())
    t = np.linspace(0, s["numerics"]["n_prod"] * s["numerics"]["dt_star"], len(ys))
    return t, ys


rows = []
spec0 = None
for d in sorted(glob.glob(PAT)):
    if not (Path(d) / "metrics.json").exists():
        continue
    s = json.loads((Path(d) / "spec.json").read_text())
    spec0 = spec0 or s
    n = int(s["params"]["n_beads"])
    om = float(s["params"]["omega_star"])
    t, ys = phasors_from_npz(d)
    src = "npz"
    if t is None:
        t, ys = phasors_from_gsd(d, n)
        src = "gsd"
    ph = []
    for i in range(n):
        blocks = LI.lockin_blocks(t, ys[:, i] - ys[:, i].mean(), om,
                                  n_blocks=min(10, max(2, len(t) // 20)))
        h, _ = LI.agg(blocks)
        ph.append(h)
    rows.append(np.array(ph))
    print(f"  {Path(d).name[-12:]}  seed={s['numerics']['seed']}  표본원={src} ({len(t)})")

Z = np.array(rows)                      # (n_seeds, n)
n = Z.shape[1]
mid = trapped_indices(n)[1]
# 구동 비드 기준 상대 위상자 (절대 위상의 ZOH·표본 지연을 상쇠)
Zrel = Z / Z[:, [mid]]

A = np.abs(Z)
A_m, A_e = A.mean(0), A.std(0, ddof=1) / np.sqrt(len(A))
phi = np.angle(Zrel)
phi_m = np.angle(Zrel.mean(0))
phi_e = np.abs(Zrel).std(0, ddof=1) / np.sqrt(len(Zrel)) / np.abs(Zrel.mean(0)).clip(1e-30)

# ── 해석적 선형응답 예측 ──
P = spec0["params"]
kth, k_t = float(P["kappa_theta_star"]), float(P["k_t_star"])
ell = 1.0 + float(P["h_min_star"])
amp, om = float(P["amp_star"]), float(P["omega_star"])
Ab = bending_matrix(n, kth, ell)
for i in trapped_indices(n):
    Ab[i, i] += k_t
e = np.zeros(n, dtype=complex); e[mid] = k_t * amp
ypred = np.linalg.solve(1j * om * np.eye(n) + Ab, e)
Apred = np.abs(ypred)
phipred = np.angle(ypred / ypred[mid])

s_idx = np.arange(n) - mid
print()
print("=" * 88)
print(f"JKR 비드별 A·sin(ωt+φ) 피팅  (시드 {len(Z)}개, n={n}, 구동=비드 {mid})")
print("=" * 88)
print(f"{'비드':>4} {'s=i-mid':>8} {'A 측정 [d]':>18} {'A 예측':>10} {'차이%':>8} "
      f"{'φ 측정 [°]':>13} {'φ 예측':>9} {'차이°':>7}")
for i in range(n):
    dA = 100 * (A_m[i] - Apred[i]) / Apred[i]
    dphi = np.degrees(np.angle(np.exp(1j * (phi_m[i] - phipred[i]))))
    print(f"{i:>4} {s_idx[i]:>8} {A_m[i]:>11.5f}±{A_e[i]:<6.4f} {Apred[i]:>10.5f} "
          f"{dA:>+8.2f} {np.degrees(phi_m[i]):>13.2f} {np.degrees(phipred[i]):>9.2f} {dphi:>+7.2f}")

np.savez("/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/"
         "7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/mode_profile.npz",
         A_m=A_m, A_e=A_e, phi_m=phi_m, phi_e=phi_e,
         Apred=Apred, phipred=phipred, s_idx=s_idx, mid=mid, n=n)
print("\n저장: mode_profile.npz")
