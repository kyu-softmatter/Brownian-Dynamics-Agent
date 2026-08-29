"""이동 트랩 골든 테스트 — `bdbot.traps.make_trap(velocity=...)` 이 맞는가.

CLAUDE.md 규칙 7: **독립적인 요소는 하나씩 떼어 검증한다.** `trap-drag-2d-hex300` 은
(이동 트랩) + (소프트 육방 격자) 두 요소의 조합인데, 조합에는 해석해가 없습니다.
**이동 트랩만 켠 최소 구성**에는 있습니다.

등속으로 끌리는 과감쇠 입자:
    γẋ = −k(x − vt) + ξ ,   u ≡ x − vt
    ⟹ γu̇ = −k u − γv + ξ
    ⟹ ⟨u⟩ = −γv/k          (뒤처짐. 부호가 음수인 것이 핵심 — 조합에서 이걸 틀렸다)
       Var(u) = kT/k        ★ **끌어도 분산은 변하지 않는다** (드리프트는 평균만 옮긴다)
       Var(y) = kT/k

무차원(γ=kT=1)으로 ⟨Δx⟩ = −v*/k*, Var = 1/k*.

두 가지를 노립니다:
  ① **부호와 크기** — 조합 코드에서 `F_drag = +k⟨Δx⟩` 로 썼다가 맨 Stokes 대비 −493%
     가 나왔습니다. 여기서 부호 규약을 못박습니다.
  ② **최소 이미지 (함정 1)** — 트랩 중심이 `anchor0 + v·t` 로 **박스를 몇 번이고 넘어갑니다.**
     고정 트랩에서는 입자만 wrap 됐지만 여기서는 **앵커가 무한히 멀어집니다.**
     래핑이 틀리면 약한 트랩에서 조용히 깨집니다 — 그래서 약한 k 로 검증합니다.

    $PY scratch/verify_moving_trap.py
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import sim as SIM, traps as TR  # noqa: E402

L = 40.0
N = 1000                 # 독립 입자 — 통계를 N배로 늘린다 (상호작용 없음)
SEED = 20260804


def run(k_star, v_star, n_tau=200.0, n_samples=200, dt_over_tau=1e-2, seed=SEED):
    """N개의 독립 입자를 각자의 이동 트랩으로 끈다. (측정 dict)

    ★ `seed` 를 조합마다 **다르게** 줘야 합니다. γ=kT=1 이고 dt = h·τ_k 이면 u 의 동역학이
      τ_k 단위로 **k에 무관하게 같아지므로**, 같은 시드를 쓰면 네 행이 전부 같은 궤적이 되어
      독립 검사가 아니라 한 검사를 네 번 인쇄한 것이 됩니다 (처음에 실제로 그랬고,
      네 행의 오차가 +0.15% 로 소수점까지 같아서 알아챘습니다).
      ⚠️ 함정 12: HOOMD 시드는 16비트로 잘리므로 작은 연속 정수를 씁니다.
    """
    tau_k = 1.0 / k_star                      # γ/k, γ=1
    dt = dt_over_tau * tau_k
    n_steps = int(round(n_tau * tau_k / dt))
    every = max(1, n_steps // n_samples)

    rng = np.random.default_rng(seed)
    pos = rng.uniform(-L / 2, L / 2, (N, 2))
    sim = SIM.make_sim(SIM.frame_2d(pos, L), seed=seed)

    vel = np.zeros((N, 3))
    vel[:, 0] = v_star
    trap = TR.make_trap(k_star, pos, L, dt_star=dt, velocity=vel)
    SIM.attach_brownian(sim, dt, [trap])     # ★ 페어 힘 없음 — 트랩만 켠 최소 구성

    sim.run(int(20 * tau_k / dt))             # 평형화 20 τ_k
    dx, dy = [], []
    for _ in range(n_samples):
        sim.run(every)
        d = trap.displacement(sim.state, sim.timestep)
        dx.append(d[:, 0].copy())
        dy.append(d[:, 1].copy())
    dx, dy = np.array(dx), np.array(dy)

    n_ind = N * n_tau / 2.0                   # 상관시간 2τ_k → 독립표본 수
    return {
        "mean_x": float(dx.mean()), "sem_x": float(dx.std() / math.sqrt(n_ind)),
        "var_x": float(dx.var()), "var_y": float(dy.var()),
        "traverse": v_star * n_tau * tau_k / L,      # 트랩이 박스를 몇 번 건넜나
    }


print("=" * 92)
print("이동 트랩 골든 테스트 — ⟨Δx⟩ = −v/k · Var = 1/k (끌어도 분산 불변)")
print("=" * 92)
print(f"  N={N} 독립 입자 · 박스 L={L:g} · 페어 힘 없음 · dt = 10⁻²τ_k")
print("  ★ 판정은 고정 % 가 아니라 **자기 통계오차(SEM) 대비 몇 σ 인가** 로 합니다 —")
print("    SNR=0.0985 인 계에서 '2% 이내'는 애초에 통계가 허용하지 않는 요구입니다.")
print(f"\n{'k*':>7}{'v*':>7}{'SNR':>6}{'횡단':>6} | "
      f"{'⟨Δx⟩ 측정':>13}{'예측 −v/k':>12}{'오차':>8}{'σ':>7} | "
      f"{'Var_x·k':>9}{'Var_y·k':>9}{'이론':>8}")
rows, ok = [], True
# ★ 약한 트랩일수록 최소이미지 결함에 취약하다 (함정 1: k=2 에서 +1856%).
#   SNR = v/√k 를 1 근처로 맞춰 통계가 나오게 한다. 마지막만 실제 케이스 값.
VAR_TH = 1.0 / (1.0 - 1e-2 / 2)          # 이산 OU 정상분산 편향 (bd-physics §1.2)
for i, (k_star, v_star) in enumerate(((2.0, 1.5), (5.0, 2.5), (10.0, 3.5),
                                      (60358.0, 24.205))):
    r = run(k_star, v_star, seed=101 + i)          # ★ 행마다 다른 시드 (위 도크스트링)
    pred = -v_star / k_star
    err = 100 * (r["mean_x"] - pred) / abs(pred)
    nsig = abs(r["mean_x"] - pred) / r["sem_x"]
    vx, vy = r["var_x"] * k_star, r["var_y"] * k_star
    # 분산의 상대 통계오차 ≈ √(2/n_indep)
    v_tol = 4 * math.sqrt(2 / (N * 200.0 / 2))
    good = nsig < 3.0 and abs(vx - VAR_TH) < v_tol and abs(vy - VAR_TH) < v_tol
    ok &= good
    rows.append((k_star, err, nsig, vx, vy, good))
    print(f"{k_star:>7g}{v_star:>7g}{v_star/math.sqrt(k_star):>6.2f}"
          f"{r['traverse']:>6.1f} | {r['mean_x']:>13.6f}{pred:>12.6f}{err:>+7.2f}%"
          f"{nsig:>6.1f}σ | {vx:>9.4f}{vy:>9.4f}{VAR_TH:>8.4f}   {'✓' if good else '✗'}")

# ── dt 수렴 — 남은 계통 오차가 O(dt) 이산화인가, 아니면 버그인가 ────────────
#   위에서 k=2·5·10 이 전부 +0.15% 로 **일정한** 계통 오차를 보였다. 크기는 작지만
#   "일정하다"는 것이 신경 쓰인다 — dt 를 절반으로 줄여 같이 절반이 되면 이산화이고,
#   안 변하면 규약 버그다. 추측하지 않고 돌려서 가른다 (CLAUDE.md 규칙 6).
#   ⭐️ **분산**이 dt 검사의 자리입니다 — 평균이 아니라. 이산 OU 의 평균은 `⟨u⟩ = −vγ/k`
#      로 **dt 에 무관하게 정확**하고(위 유도), 편향은 분산에만 `1/(1−h/2)` 로 들어갑니다.
#      게다가 분산은 상대 통계오차가 √(2/n) ≈ 0.45% 라 평균(0.28% of signal)보다 날카롭습니다.
print("\n  dt 수렴 (k*=5, v*=2.5) — 분산 편향이 이론 1/(1−h/2) 를 따르는가?")
print(f"    {'dt/τ_k':>8}{'Var_x·k':>10}{'이론':>9}{'차이':>9}{'통계오차':>9}")
conv_ok = True
for j, hh in enumerate((4e-2, 2e-2, 1e-2, 5e-3)):
    r = run(5.0, 2.5, dt_over_tau=hh, seed=211 + j)
    th = 1.0 / (1.0 - hh / 2)
    got = r["var_x"] * 5
    tol = 3 * math.sqrt(2 / (N * 200.0 / 2))
    good = abs(got - th) < tol
    conv_ok &= good
    print(f"    {hh:>8.0e}{got:>10.4f}{th:>9.4f}{got-th:>+9.4f}{tol:>9.4f}   "
          f"{'✓' if good else '✗'}")
ok &= conv_ok
print(f"    → {'✓ 분산 편향이 이론과 일치 — 끌어도 분산은 kT/k 그대로다' if conv_ok else '✗ 이론과 어긋난다'}")

print(f"""
  판정 기준: ⟨Δx⟩ 가 예측에서 3σ 이내 · Var 편향이 이론 1/(1−h/2) 와 통계오차 이내
  ⚠️ 처음에 '오차 < 2%' 같은 고정 기준을 쓰고 **모든 행에 같은 시드**를 썼더니
     네 행이 전부 +0.15% 로 소수점까지 같았습니다 — 독립 검사가 아니라 한 검사를
     네 번 인쇄한 것이었고, dt 를 줄여도 안 변해서 '버그'로 오판할 뻔했습니다.
  ★ 마지막 행은 실제 케이스 값(k*=6.04e4, v*=24.2, SNR=0.0985)이다 — 신호가 잡음의
    1/10 이라 N={N}개를 평균해야 겨우 보인다. L3가 예고한 그 통계 문제가 여기서도 그대로다.""")
print("=" * 92)
print(f"{'✓ PASS' if ok else '✗ FAIL'} — 이동 트랩 {sum(1 for r in rows if r[4])}/{len(rows)}")
print("=" * 92)
sys.exit(0 if ok else 1)
