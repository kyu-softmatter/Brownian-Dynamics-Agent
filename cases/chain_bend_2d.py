"""`chain-bend-2d-oscill` — L3 무차원화 + L4 빌더 (실행은 ⑤ 때문에 현재 거부됨).

광집게 3점 굽힘으로 콜로이드 비드 사슬을 y방향으로 진동시켜 G'(ω)·G''(ω)를 재는 계.
양끝 비드는 고정 트랩(힘 센서), 중앙 비드의 트랩을 `y = a sin(ωt)` 로 구동.

물리는 사용자 지정 논문 2편에서 왔습니다 (스케치의 `U_ij` 는 물결선 = 빈칸이었습니다):
  [P1] Pantina & Furst, PRL 94, 138301 (2005)
  [P2] Pantina & Furst, Langmuir 24, 1141-1146 (2008)
★ 결합은 **페어 퍼텐셜이 아닙니다** — 접착 접촉(JKR) + 접선 굽힘 강성 κ_θ = EI/ℓ 입니다.
추출 검증은 `scratch/chain_bend_from_papers.py` (논문 κ₀ 재현 +1.6%, 이산↔연속 매핑은
n=25 에서 −0.35%. system.yaml 이 인용한 −0.08% 는 그 수렴표의 n=51 행입니다).

**앞의 세 케이스와 근본적으로 다른 점**:

  ① `dt` 를 정하는 것이 **관측 대상이 아닙니다.** 강성 행렬의 최대 고유값에서 나오는
     최속 굽힘 모드(τ_fast = 0.279 µs)가 dt를 정하는데, 재려는 것은 τ_chain(1.27 ms)의
     집단 굽힘입니다 — 4570배 떨어져 있습니다. 1-A·1-B·trap-drag 는 전부 dt를 정하는
     척도가 관심 척도였거나 그 근처였습니다.
  ② ★ 최속 모드가 과감쇠가 아닙니다 (τ_p/τ_fast = 0.60, ζ=0.65). BD는 그 모드를 과감쇠로
     다루므로 그 대역의 동역학은 틀리고 **어떤 dt로도 고쳐지지 않습니다**.
     ✔ **해결됨** — `OverdampedViscous` vs `Langevin(kT=0)` 을 7개 ω 전부에서 비교해
     K*(ω) 차이가 **최대 0.159%**. 관측 대역에 영향 없음이 측정됐습니다
     (`scratch/verify_chain_bend_gates.py --gate det --collect`).
  ③ 주기경계가 없습니다 → `box` 역할을 `declare_absent` 로 명시합니다.
  ④ 진폭이 **세 방향에서** 조여 있습니다: ℓ_k ≪ a < δ_max **이면서** min|θ−π| > SMALL.
     위쪽 한계는 M_c (넘으면 결합이 미끄러짐), 아래쪽은 SNR, 그리고 ⑤가 새 제약입니다.

⛔ **⑤ 이 케이스는 지금 `angle.Harmonic` 으로 실행할 수 없습니다** (2026-08-05).
   HOOMD 가 `sin θ` 를 `ANGLE_SIN_SMALL`(=√2×10⁻³)로 클램프해서, 거의 곧은 사슬의
   **힘만** 축소됩니다 (에너지는 정확 → 에너지 검증으로는 안 잡힘). 이 계는 응답
   프로파일의 **23개 각도 전부**가 그 영역입니다 (min|θ−π| = 7.3e-5 = SMALL의 1/19).
   하드 검사 `★angle 힘 유효` 가 막으므로 `--spec` 은 스펙을 쓰지 않고 `--run` 은 거부합니다.
   L4 구성 자체(유령 트랩·구동·락인)는 관문 2종으로 검증돼 있고 **빌더는 완성돼 있습니다** —
   막힌 것은 굽힘 힘 하나입니다. 상세·우회: `assert_angle_force_valid()` · skill bd-hoomd 함정 15.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/chain_bend_2d.py --report              # ω = ω_min (가장 비싼 점)
    $PY cases/chain_bend_2d.py --omega 7853 --report
    $PY cases/chain_bend_2d.py --sweep --spec        # 스윕 전체를 스펙으로
    $PY cases/chain_bend_2d.py --run                 # L4 — 현재 ⑤ 때문에 거부됨
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bdbot import checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import lockin as LI, run as RUN, sim as SIM  # noqa: E402
from bdbot import nondim as ND, scales as SC  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from bdbot.units import Q  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
N_CYCLES = 100.0            # 관측창 = 최저 ω에서 몇 주기
N_SWEEP = 7                 # ω 스윕 점 개수 (로그 등간격)
AC_GATE = 1e-1              # ★제안: 접촉을 점(hinge)으로 볼 수 있는가 — a_c/d 상한
THETA_GATE = 1e-1           # ★제안: 조화 angle 의 소각 선형화 한계 [rad]

# ★★ HOOMD 하드 상수 (실측 확정, `scratch/verify_angle_force_small_theta.py`).
# md.angle.Harmonic 은 토크를 좌표로 옮길 때 1/sin θ 를 쓰고, sin θ 를 이 값으로 **클램프**한다.
# sin θ < SMALL 이면 힘이 sinθ/SMALL 배로 축소된다 → 힘 ∝ κ(θ−π)²/SMALL 로 **2차**가 되어
# 사슬이 실제보다 훨씬 무르고 비선형이 된다. **에너지는 정확하다** (0.000%) — 그래서
# 에너지로 검증하면 통과하고 힘은 틀린 채로 남는다. t0=π 는 평형 자체가 sin θ=0 이라
# 뻣뻣한 사슬에서 항상 이 영역이다. 실측 SMALL = 1.414217e-03 (표준편차 1.4e-7) = √2×10⁻³.
ANGLE_SIN_SMALL = 1.414214e-3


# ════════════════════════════════════════════════════════════════════════
# ① 물리계 (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    P = load_node
    con = raw["interactions"][0]
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "n": int(raw["geometry"]["n_beads"]),
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "E": P(raw["particle"]["youngs_modulus"]),
        "kappa_0": P(con["kappa_0"]),
        "a_c": P(con["contact_radius"]),
        "k_bond": P(con["bond_stiffness"]),
        "M_c": P(con["critical_moment"]),
        "k_t": P(raw["external"]["stiffness"]),
        "amp": P(raw["external"]["amplitude"]),
        "omega_range": [float(x) for x in raw["external"]["omega_range"]["value"]],
        "n_trapped": int(raw["external"]["n_trapped"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② 이 계 고유의 수치 — 이산 사슬의 강성 행렬
#    ★ 감으로 쓰지 않습니다. 논문 값을 재현한 scratch/chain_bend_from_papers.py 와
#      같은 구성이고, 여기서는 **고유값**까지 씁니다 (거기서는 강성만 썼습니다).
# ════════════════════════════════════════════════════════════════════════
def bending_matrix(n: int, kappa_theta: float, ell: float) -> np.ndarray:
    """U = ½ κ_θ Σ θ_i²,  θ_i = (y_{i+1} − 2y_i + y_{i−1})/ℓ  의 2차형식 A (y 방향).

    소변형 근사. `A = κ_θ Bᵀ B`, B 는 2차 미분 행렬 (n−2)×n.
    """
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def bond_matrix(n: int, k_bond: float) -> np.ndarray:
    """결합 신축(x 방향)의 2차형식 — 경로 그래프 라플라시안 × k_b.

    ★ 직선 사슬에서 신축(x)과 굽힘(y)은 **선형 차수에서 분리**됩니다. 그래서 최대
      고유값을 두 블록에서 따로 구해 큰 쪽을 씁니다. 한 행렬로 더하면 (분리된 자유도를
      섞어) λ_max 를 과대평가하고 dt 가 필요 이상으로 작아집니다.
    """
    G = np.zeros((n, n))
    for i in range(n - 1):
        G[i, i] += 1.0
        G[i + 1, i + 1] += 1.0
        G[i, i + 1] -= 1.0
        G[i + 1, i] -= 1.0
    return k_bond * G


def trapped_indices(n: int) -> list[int]:
    """[P2] Fig.1A — 양끝 2개(고정, 힘 센서) + 중앙 1개(구동). 총 3개."""
    return [0, n // 2, n - 1]


def three_point_bending(n: int, kappa_theta: float, ell: float, delta: float):
    """양끝 고정 + 중앙을 delta 만큼 밀 때의 (중앙 강성 F/δ, 최대 결합각 |θ|).

    `scratch/chain_bend_from_papers.discrete_3point_stiffness` 와 같은 구성이고
    (그 수렴표에서 빔 공식과 대조 완료 — **n=25 는 −0.35%**, n=51 이 −0.08%),
    여기서는 각도까지 함께 돌려줍니다 — 조화 angle 퍼텐셜의 소각 선형화가 성립하는지
    보려면 각도가 필요합니다.
    """
    A = bending_matrix(n, kappa_theta, ell)
    fixed = [0, n // 2, n - 1]
    free = [i for i in range(n) if i not in fixed]
    y_fix = np.array([0.0, delta, 0.0])
    y = np.zeros(n)
    y[fixed] = y_fix
    y[free] = np.linalg.solve(A[np.ix_(free, free)], -A[np.ix_(free, fixed)] @ y_fix)
    U = 0.5 * y @ A @ y
    k_center = 2 * U / delta**2                      # U = ½ k δ²
    theta = np.abs(np.diff(y, n=2)) / ell            # |θ_i|
    return k_center, float(theta.max())


def trapped_stiffness_matrix(n: int, kappa_theta: float, ell: float, k_t: float):
    """굽힘 + 트랩 3개의 강성 행렬과 구동 비드 인덱스.

    ★ 양끝은 **강체 고정이 아니라 트랩**이다 (유한 강성 k_t). 그 차이가 30% 넘게 나므로
      경계조건을 강체로 두고 유도한 값을 "트랩이 느끼는 값"이라고 부르면 안 된다.
    """
    A = bending_matrix(n, kappa_theta, ell)
    idx = trapped_indices(n)
    for i in idx:
        A[i, i] += k_t
    return A, idx[1]


def driven_static_stiffness(n: int, kappa_theta: float, ell: float, k_t: float) -> float:
    """구동 트랩이 **실제로** 느끼는 정적 강성 K(ω→0).

    (A_bend + T) y = k_t y_c e_mid 를 풀고 K = k_t(y_c/y_mid − 1). 양끝 트랩이 유한
    강성이라 사슬이 통째로 밀리므로 강체 고정 가정의 48EI/L³ 보다 **작다**.
    실측 대조: 최저 ω(De_true≈1) 에서 HOOMD 가 이 값의 0.95배
    (`scratch/verify_chain_bend_gates.py --gate det`).
    """
    A, mid = trapped_stiffness_matrix(n, kappa_theta, ell, k_t)
    e = np.zeros(n)
    e[mid] = k_t
    return float(k_t * (1.0 / np.linalg.solve(A, e)[mid] - 1.0))


def driven_response(n: int, kappa_theta: float, ell: float, k_t: float,
                    gamma: float, omega: float, amp: float) -> float:
    """구동 비드의 응답 진폭 |ŷ(ω)| — 소각 선형응답 (iωγI + A + T) ŷ = k_t a e_mid.

    ★ SNR 을 재는 데 쓴다. 스펙이 예전에 검사했던 a/ℓ_k 는 **구동 진폭**이고, 관측량은
      **응답**이다 — 사슬과 항력이 저항하므로 |ŷ| ≪ a 이고 고주파에서 ℓ_k 아래로 내려간다.

    ✔ **원인 규명됨 (2026-08-05).** 예전에 HOOMD 실측과 28% 어긋났던 것은 이 예측이 아니라
      **HOOMD 쪽이 틀린 것**이었다 — `md.angle.Harmonic` 이 sin θ 를 `ANGLE_SIN_SMALL`
      로 클램프해서 거의 곧은 사슬의 힘을 축소한다. 이 계는 23개 각도 전부가 그 영역이다.
      이 모델은 정확 비선형 최소화(scipy, 전 2n 좌표)와 **0.32%** 일치하므로 옳다.
      추적: `scratch/diagnose_chain_bend_28pct.py` · `scratch/verify_angle_force_small_theta.py`
    """
    A, mid = trapped_stiffness_matrix(n, kappa_theta, ell, k_t)
    e = np.zeros(n, dtype=complex)
    e[mid] = k_t * amp
    M = 1j * omega * gamma * np.eye(n) + A
    y = np.linalg.solve(M, e)
    return float(abs(y[mid])), np.abs(y)


def response_angles(prof: np.ndarray, ell: float) -> tuple[float, float]:
    """응답 진폭 프로파일에서 결합각 편차 |θ−π| 의 (최대, 최소).

    ★★ `ANGLE_SIN_SMALL` 판정에 쓴다. 힘이 맞으려면 **모든** 각도가 그 위에 있어야 한다 —
    최대만 보면 안 된다. 사슬 양끝의 각도가 중앙보다 한 자릿수 이상 작다.
    """
    th = np.abs(np.diff(prof, n=2)) / ell
    return float(th.max()), float(th.min())


# ════════════════════════════════════════════════════════════════════════
# ③ 스케일 원장 (bd-physics §0 ①②)
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, omega: float, *, dt_scale=1.0, n_cycles=N_CYCLES) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B = b["kT"], b["gamma"], b["tau_B"]
    k_t = sys_["k_t"].value.to("N/m")
    amp = sys_["amp"].value.to("m")
    n = sys_["n"]
    a_rad = d / 2                                    # 입자 반지름 (논문의 a)

    # ── 접촉 → 굽힘 강성 ([P2] eq.5 + 이산 매핑) ──────────────────────
    EI = (sys_["kappa_0"].value.to("N/m") * a_rad**3 / 3).to("N*m^2")
    kappa_theta = (EI / d).to("J")                   # κ_θ = EI/ℓ,  ℓ = d
    L_chain = (n - 1) * d                            # 윤곽길이
    # ★ 힘 정의를 섞으면 논문값과 정확히 2배 어긋난다 (chain_bend_from_papers §②).
    kappa_end = (24 * EI / L_chain**3).to("N/m")     # 논문 정의 (끝 힘 기준)
    # ★ 48EI/L³ 은 **양끝 강체 고정** 3점굽힘 값이다. 이 계의 양끝은 트랩(k_t)이므로
    #   구동 트랩이 느끼는 값이 아니다 — 그건 아래 kappa_drive 다 (32% 차이).
    kappa_center = 2 * kappa_end                     # 빔 공식 48EI/L³ (강체 고정 가정)
    M_c = sys_["M_c"].value.to("N*m")
    delta_max = (M_c * L_chain**2 / (12 * EI)).to("m")    # M<M_c 선형 탄성 한계 진폭

    # ── 시간척도 ──────────────────────────────────────────────────────
    tau_k = C.relaxation_time(gamma, k_t)                       # 트랩
    tau_chain = C.relaxation_time(gamma, kappa_center)          # γ/κ_center (빔 공식 기준)
    k_b = sys_["k_bond"].value.to("N/m")

    # ★ 최속 모드 — 강성 행렬의 최대 고유값. 신축(x)과 굽힘(y)이 분리되므로 따로 구한다.
    idx = trapped_indices(n)
    kth = float(kappa_theta.magnitude)
    ell = float(d.magnitude)
    Ay = bending_matrix(n, kth, ell)
    Ax = bond_matrix(n, float(k_b.magnitude))
    for i in idx:                                    # 트랩은 대각에 k_t 를 더한다
        Ay[i, i] += float(k_t.magnitude)
        Ax[i, i] += float(k_t.magnitude)
    ev_bend = np.linalg.eigvalsh(Ay)
    lam_bend = float(ev_bend[-1])
    lam_bond = float(np.linalg.eigvalsh(Ax)[-1])
    lam_max = max(lam_bend, lam_bond)
    tau_fast = C.relaxation_time(gamma, Q(lam_max, "N/m"))
    tau_bond = C.relaxation_time(gamma, Q(lam_bond, "N/m"))

    # ★★ 최저 고유값 → **최장** 이완시간. 예전에는 이걸 구하지 않고 τ_chain = γ/κ_center
    #    을 지배 척도로 썼는데, τ_max 가 그보다 9.18배 길다. 결과로 (a) De 가 9.18배
    #    과소평가되어 ω 스윕이 준정적 영역에 아예 못 들어갔고 (b) 평형화(20 τ_chain =
    #    2.2 τ_max)가 부족해 K* 가 최대 21% 틀렸다. 실측 확인:
    #    `scratch/verify_chain_bend_gates.py --gate det --eq-steps` (평형화를 10 τ_max
    #    로 늘리면 블록 산포가 1000배 줄고 K′ 이 21% 이동한다).
    lam_min = float(ev_bend[0])
    tau_max = C.relaxation_time(gamma, Q(lam_min, "N/m"))
    # 스윕 **전체**의 저주파 끝이 준정적에 닿는가 — 이 스펙의 ω 가 아니라 범위의 하한으로
    # 판정한다 (점마다 De 를 보면 최저점이 통과해 버려서 범위 문제가 안 드러난다).
    de_lo = float((Q(min(sys_["omega_range"]), "1/s") * tau_max).to("dimensionless").magnitude)
    kappa_drive = driven_static_stiffness(n, kth, ell, float(k_t.magnitude))
    y_resp, y_prof = driven_response(n, kth, ell, float(k_t.magnitude),
                                     float(gamma.to("kg/s").magnitude), omega,
                                     float(amp.magnitude))
    th_hi, th_lo = response_angles(y_prof, ell)      # ★★ HOOMD angle 힘 유효성 판정용

    dt = dt_scale * C.dt_from_gate(tau_fast)         # 최속 모드가 dt를 정한다
    tau_w = Q(1.0 / omega, "s")                      # 구동 (De = τ_chain/τ_w)
    tau_period = Q(2 * math.pi / omega, "s")
    T_obs = (n_cycles * tau_period).to("s")

    l_k = (kT / k_t) ** 0.5
    k_center_disc, theta_max = three_point_bending(n, kth, ell, float(amp.magnitude))

    lg = SC.ScaleLedger()
    lg.add_length("l_k", l_k.to("m"), "√(kT/k_t) 트랩 요동 (잡음 하한)", star=True)
    lg.add_length("y_resp", Q(y_resp, "m"),
                  "★ 구동 비드 응답 진폭 |ŷ(ω)| — 관측량. a 가 아니라 이것이 SNR 의 분자다"
                  " (선형응답 추정, HOOMD 와 고주파에서 28% 차 — 미해명)", star=True)
    lg.add_length("a_c", sys_["a_c"].value.to("m"), "JKR 접촉 반경 (점 hinge 가정)")
    lg.add_length("a", amp, "구동 진폭", star=True)
    lg.add_length("delta_max", delta_max, "M_c 선형 탄성 한계 진폭", star=True)
    lg.add_length("d", d, "비드 지름 = 결합 길이 ℓ (기준)")
    lg.add_length("L_chain", L_chain.to("m"), "사슬 윤곽길이 (n−1)d")
    lg.add_time("tau_p", b["tau_p"], "m/γ 관성 이완", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_fast", tau_fast, "γ/λ_max 최속 굽힘 모드 — dt를 정한다", star=True)
    lg.add_time("tau_bond", tau_bond, "γ/λ_max(신축) 결합 신축")
    lg.add_time("tau_w", tau_w, f"1/ω 구동 (ω = {omega:.0f} rad/s)")
    lg.add_time("tau_k", tau_k, "γ/k_t 트랩")
    lg.add_time("tau_chain", tau_chain, "γ/κ_center (빔 공식 기준 — 지배 척도가 아니다)")
    lg.add_time("tau_max", tau_max,
                "★★ γ/λ_min 최장 이완시간 — **지배 척도**. De·평형화가 이걸 써야 한다",
                star=True)
    lg.add_time("tau_period", tau_period, "2π/ω 구동 주기")
    lg.add_time("tau_B", tau_B, "d²/D_t 확산 (기준)")
    lg.add_time("T_obs", T_obs, f"관측창 ({n_cycles:g}주기)", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.add_energy("k_t_d2", (k_t * d**2).to("J"), "k_t d² 트랩 강성")
    lg.add_energy("kappa_end_d2", (kappa_end * d**2).to("J"), "κ_end d² 사슬 강성 (논문 정의)")
    lg.add_energy("kappa_drive_d2", Q(kappa_drive, "N/m").to("N/m") * d**2,
                  "★ κ_drive d² — 구동 트랩이 **실제로** 느끼는 정적 강성 "
                  "(양끝 트랩이 유한강성. 강체 고정 48EI/L³ 의 0.68배)", star=True)
    lg.add_energy("k_b_d2", (k_b * d**2).to("J"), "k_b d² 결합 신축 강성")
    lg.add_energy("M_c", M_c, "임계 굽힘 모멘트 (모멘트 = 에너지 차원)")
    lg.add_energy("kappa_theta", kappa_theta, "결합각 강성 κ_θ = EI/ℓ", star=True)
    # ★ 주기경계가 없다 — 사슬 하나이고 트랩이 위치를 고정한다. 지어내지 않고 비운다.
    lg.declare_absent(
        "box",
        "주기경계 없음 (geometry.periodic=false). 사슬 하나이고 트랩이 위치를 고정하므로 "
        "최소이미지·유한크기 검사의 분모가 되는 박스가 물리에 없다. L4는 HOOMD 박스를 "
        "사슬 최대 확장 + WCA 컷오프보다 크게만 잡으면 되고, 그 값은 물리를 바꾸지 않는다. "
        "이 계에서 기하 한계 역할을 하는 것은 박스가 아니라 δ_max (M<M_c) 다.")

    lg.derived = dict(gamma=gamma, D_t=b["D_t"], m=b["m"], kT=kT, d=d, tau_B=tau_B,
                      EI=EI, kappa_theta=kappa_theta, kappa_end=kappa_end,
                      kappa_center=kappa_center, delta_max=delta_max, M_c=M_c,
                      L_chain=L_chain.to("m"), l_k=l_k.to("m"), k_t=k_t, k_b=k_b,
                      tau_k=tau_k, tau_chain=tau_chain, tau_fast=tau_fast,
                      tau_max=tau_max, lam_min=lam_min, kappa_drive=kappa_drive,
                      y_resp=y_resp, de_lo=de_lo, th_hi=th_hi, th_lo=th_lo,
                      tau_bond=tau_bond, tau_w=tau_w, tau_period=tau_period,
                      dt=dt, T_obs=T_obs, omega=omega, n=n, amp=amp,
                      lam_bend=lam_bend, lam_bond=lam_bond, lam_max=lam_max,
                      k_center_disc=k_center_disc, theta_max=theta_max,
                      trapped=trapped_indices(n))
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ 이 계는 시간척도가 4자릿수에 걸쳐 있고, dt를 정하는 "
        "척도(τ_fast, 최속 굽힘 모드)가 **관측 대상이 아니다** — 재려는 것은 τ_chain 이다. "
        "그래서 비용이 관측 대역과 무관하게 결정된다.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ④ 무차원수 + ⑤ 분리 검사 (bd-physics §3, §4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg):
    n = sys_["n"]
    D = lg.derived
    r = lg.ratio
    n = sys_["n"]

    groups = [
        ND.Group("kappa_theta/kT", r("energies", "kappa_theta", "kT"),
                 ("energies", "kappa_theta"), ("energies", "kT"), "EI/(ℓ kT)",
                 "★ 사슬이 열적으로 완전히 뻣뻣 (요동이 굽히지 못한다)"),
        ND.Group("M_c/kT", r("energies", "M_c", "kT"),
                 ("energies", "M_c"), ("energies", "kT"), "",
                 "결합이 미끄러지는 모멘트 (조화 angle 의 유효 범위)"),
        ND.Group("k*", r("energies", "k_t_d2", "kT"),
                 ("energies", "k_t_d2"), ("energies", "kT"), "k_t d²/kT",
                 "트랩 vs 열요동"),
        ND.Group("kappa_end/k_t", r("energies", "kappa_end_d2", "k_t_d2"),
                 ("energies", "kappa_end_d2"), ("energies", "k_t_d2"), "",
                 "★ 사슬 vs 트랩 강성 — 비슷해야 측정 창에 들어온다"),
        ND.Group("k_b/k_t", r("energies", "k_b_d2", "k_t_d2"),
                 ("energies", "k_b_d2"), ("energies", "k_t_d2"), "",
                 "결합 신축은 사실상 비신축 (트랩보다 192배 뻣뻣)"),
        ND.Group("a/l_k", r("lengths", "a", "l_k"), ("lengths", "a"), ("lengths", "l_k"),
                 "", "★ SNR — 진폭 vs 트랩 안 열요동"),
        ND.Group("a/delta_max", r("lengths", "a", "delta_max"),
                 ("lengths", "a"), ("lengths", "delta_max"), "",
                 "★ 선형 탄성 여유 (M<M_c). 1을 넘으면 조화 angle 이 무효"),
        ND.Group("a/d", r("lengths", "a", "d"), ("lengths", "a"), ("lengths", "d"), "",
                 "진폭 vs 비드"),
        ND.Group("a_c/d", r("lengths", "a_c", "d"), ("lengths", "a_c"), ("lengths", "d"), "",
                 "접촉이 점(hinge)인가 — JKR 유도의 전제"),
        ND.Group("L_chain/d", r("lengths", "L_chain", "d"),
                 ("lengths", "L_chain"), ("lengths", "d"), "n−1", "사슬 = 비드 몇 개분"),
        ND.Group("n_beads", float(n), None, None, "", "비드 수 (입력, ★제안)"),
        ND.Group("De", r("times", "tau_max", "tau_w"),
                 ("times", "tau_max"), ("times", "tau_w"), "ω τ_max",
                 "★★ Deborah — **최장 이완시간 기준**. 준정적 극한은 De ≪ 1"),
        ND.Group("De_chain_old", r("times", "tau_chain", "tau_w"),
                 ("times", "tau_chain"), ("times", "tau_w"), "ω τ_chain",
                 "★ 예전 정의 (γ/κ_center 기준). De 를 9.18배 과소평가한다 — 보존만"),
        ND.Group("tau_max/tau_chain", r("times", "tau_max", "tau_chain"),
                 ("times", "tau_max"), ("times", "tau_chain"), "κ_center/λ_min",
                 "★★ 지배 척도가 빔 공식보다 이만큼 길다 — 예전 스펙이 놓친 인수"),
        ND.Group("kappa_drive/k_t", r("energies", "kappa_drive_d2", "k_t_d2"),
                 ("energies", "kappa_drive_d2"), ("energies", "k_t_d2"), "",
                 "★ 사슬 vs 트랩 — 트랩 경계를 반영한 값 (κ_end/k_t 가 아니라 이것)"),
        ND.Group("y_resp/l_k", r("lengths", "y_resp", "l_k"),
                 ("lengths", "y_resp"), ("lengths", "l_k"), "",
                 "★★ 실제 SNR — 응답이 열요동보다 큰가. a/ℓ_k 는 이걸 최대 60배 과대평가"),
        ND.Group("De_trap", r("times", "tau_k", "tau_w"),
                 ("times", "tau_k"), ("times", "tau_w"), "ω τ_k", "트랩 기준"),
        ND.Group("tau_fast/tau_chain", r("times", "tau_fast", "tau_chain"),
                 ("times", "tau_fast"), ("times", "tau_chain"), "",
                 "★ 척도 분리 폭 — dt를 정하는 모드가 관심 모드보다 이만큼 빠르다"),
        ND.Group("dt/tau_fast", r("times", "dt", "tau_fast"),
                 ("times", "dt"), ("times", "tau_fast"), "", "적분 해상 (지배 척도)"),
        ND.Group("n_cycles", r("times", "T_obs", "tau_period"),
                 ("times", "T_obs"), ("times", "tau_period"), "", "관측 주기 수"),
        ND.Group("St", r("times", "tau_p", "tau_B"),
                 ("times", "tau_p"), ("times", "tau_B"), "tau_p/tau_B", "관성 vs 확산"),
    ]
    checks = [
        C.Check("model", "inertia negligible   tau_p/tau_max", r("times", "tau_p", "tau_max"),
                C.GATE, "<=",
                "τ_dyn = 관측 대역의 지배 척도 = τ_max (최장 이완). G'(ω)를 재는 대역이 여기다"),
        C.Check("model", "note: tau_p/tau_fast", r("times", "tau_p", "tau_fast"), C.GATE, "<=",
                "★ 최속 굽힘 모드는 과감쇠가 아니다 (ζ = γ/2√(mλ_max) = 0.65 < 1). BD 는 그 "
                "모드를 과감쇠로 다루므로 그 대역의 동역학은 틀리고 **어떤 dt로도 고쳐지지 "
                "않는다**. ✔ **측정 완료** — 같은 파라미터로 OverdampedViscous vs "
                "Langevin(kT=0) 을 7개 ω 전부에서 비교해 K*(ω) 차이가 최대 **0.159%** "
                "(De_old=10). 관측 대역에 영향 없음이 확인됐다 "
                "(`scratch/verify_chain_bend_gates.py --gate det --collect`). "
                "열적 링잉 × 비선형 결합은 이 검정이 덮지 않는다", hard=False),
        C.Check("model", "linear elasticity    a/delta_max", r("lengths", "a", "delta_max"), 1.0, "<=",
                f"★ M < M_c. 넘으면 결합이 미끄러지거나 굴러 ([P2] 결론) 조화 angle "
                f"퍼텐셜이 무효가 된다. δ_max = M_c L²/(12EI) = "
                f"{D['delta_max'].to('nm'):~.0fP}"),
        C.Check("model", "small-angle linear   max|theta| [rad]", D["theta_max"], THETA_GATE, "<=",
                f"조화 angle 은 소각 근사다. 이산 3점 굽힘을 a={D['amp'].to('nm'):~.0fP} 로 "
                f"직접 풀어 얻은 최대 결합각. 상한 {THETA_GATE:g} rad 는 ★제안", hard=False),
        C.Check("model", "point contact        a_c/d", r("lengths", "a_c", "d"), AC_GATE, "<=",
                f"κ₀ = 3πa_c⁴E/(4a³) 는 a_c ≪ a 를 전제한다. 상한 {AC_GATE:g} 는 ★제안",
                hard=False),
        C.Check("model", "*angle force valid   min|theta-pi|", D["th_lo"], ANGLE_SIN_SMALL, ">=",
                f"★★ **HOOMD 하드 제약.** md.angle.Harmonic 은 sin θ 를 "
                f"{ANGLE_SIN_SMALL:.3e} 로 클램프해서, 그 아래에서는 힘이 sinθ/SMALL 배로 "
                f"축소된다 (힘 ∝ κ(θ−π)²  — 선형이 아니라 2차). **에너지는 정확하다**(0.000%) "
                f"→ 에너지로 검증하면 통과하고 힘만 틀린다. 응답 프로파일의 **모든** 각도가 "
                f"이 위에 있어야 한다. 실측·재현: scratch/verify_angle_force_small_theta.py. "
                f"이 계는 max|θ−π|={D['th_hi']:.2e} 인데 최소가 {D['th_lo']:.2e} 라 "
                f"{n-2}개 각도 전부가 깨진 영역이다. "
                + ("→ angle.Harmonic 으로는 실행 불가 (하드)"
                   if BENDING_IMPL == "angle_harmonic" else
                   f"→ ★ 지금은 굽힘을 `{BENDING_IMPL}`(force.Custom 으로 F = −A y 를 "
                   "직접 계산)로 구현하므로 **이 제약이 적용되지 않는다.** 검사는 지우지 "
                   "않고 참고로 남긴다 — 구현을 angle.Harmonic 으로 되돌리면 즉시 하드로 "
                   "돌아온다 (BENDING_IMPL)."),
                # ★ 이건 **구현의 제약**이지 계의 제약이 아니다. 그래서 지우지 않고
                #   구현에 조건부로 건다. 커스텀 힘의 타당성은 아래 '소각 선형화' 가 본다.
                hard=(BENDING_IMPL == "angle_harmonic")),
        C.Check("integration", "fastest mode resolved dt/tau_fast", r("times", "dt", "tau_fast"),
                C.GATE, "<=",
                f"강성 행렬 최대 고유값 λ_max = {D['lam_max']:.4e} N/m (굽힘 "
                f"{D['lam_bend']:.3e} vs 신축 {D['lam_bond']:.3e} — 큰 쪽). "
                "이걸 놓치면 발산한다"),
        C.Check("integration", "stretch resolved     dt/tau_bond", r("times", "dt", "tau_bond"),
                C.GATE, "<=", "결합 신축 모드 (굽힘보다 느려 여유가 크다)"),
        C.Check("integration", "drive resolved       dt/tau_w", r("times", "dt", "tau_w"), C.GATE, "<=",
                f"구동 ω = {D['omega']:.0f} rad/s 를 해상"),
        C.Check("statistics", "SNR   |y_hat(w)|/l_k", r("lengths", "y_resp", "l_k"), 3.0, ">=",
                "★★ **응답** 진폭이 열요동보다 커야 위상 추출이 된다. 예전 검사는 분자에 "
                "구동 진폭 a 를 써서 ω 무관한 9.83 을 돌려주고 통과했지만, 실제 SNR 은 "
                "ω 와 함께 떨어져 고주파에서 1 아래다 (실측 De_old=10 에서 0.165). "
                "trap-drag 의 '하드 검사는 통과하는데 통계가 안 나온다' 와 같은 구멍이었다. "
                "⚠ |ŷ| 은 선형응답 추정이고 HOOMD 와 고주파에서 28% 어긋난다 (미해명)",
                hard=False),
        C.Check("statistics", "drive amplitude      a/l_k", r("lengths", "a", "l_k"), 3.0, ">=",
                "구동이 열요동보다 큰가 — 필요조건이지만 **충분조건이 아니다** "
                "(위 |ŷ|/ℓ_k 가 실제 판정)", hard=False),
        C.Check("statistics", "quasi-static reached De(w_min)", D["de_lo"], 0.1, "<=",
                "★★ De = ω τ_max. 스윕이 준정적 극한(De ≪ 1)을 포함해야 K′ 의 탄성 고원이 "
                "보인다. 예전 정의(τ_chain)로는 0.1~10 을 덮는다고 나왔지만 실제로는 "
                "De ≈ 1~92 라서 고원 영역에 **아예 들어가지 않는다** — 스케치가 요구한 "
                "포화 곡선을 이 스윕으로는 못 낸다. ω 범위는 system.yaml tier 3 "
                "(사용자 승인)이라 여기서 바꾸지 않고 검사로 드러낸다", hard=False),
        C.Check("statistics", "cycles observed      T_obs/(2pi/w)", r("times", "T_obs", "tau_period"),
                N_CYCLES, ">=", "위상 평균에 쓸 주기 수", hard=False),
    ]
    return groups, checks


def report_blocks(sys_, lg, n_eq, n_prod):
    D = lg.derived
    inp = [R.kv(k, f"{sys_[k].value:~.4gP}", sys_[k].tier, sys_[k].source[:44])
           for k in ("d", "T", "eta", "kappa_0", "M_c", "k_t", "amp", "k_bond")]
    inp += [R.kv("n", f"{sys_['n']}", 3, "★제안 [P1] Fig.4 — n=11이면 진폭 창이 닫힌다"),
            R.kv("omega", f"{D['omega']:.0f} rad/s", 3, "★제안 — De 0.1~10 스윕의 한 점")]
    der = [
        f"  EI = {D['EI']:~.4eP}   κ_θ = EI/ℓ = {D['kappa_theta']:~.4eP}"
        f" = {lg.ratio('energies', 'kappa_theta', 'kT'):.3e} kT",
        f"  κ_end = {D['kappa_end'].to('pN/um'):~.3fP} (논문 정의, 끝 힘)"
        f"   κ_center = {D['kappa_center'].to('pN/um'):~.3fP} (빔 48EI/L³, 강체 고정 가정)",
        f"  ★ κ_drive = {Q(D['kappa_drive'], 'N/m').to('pN/um'):~.3fP}"
        f" = κ_center × {D['kappa_drive']/float(D['kappa_center'].magnitude):.3f}"
        f"  ← 구동 트랩이 **실제로** 느끼는 값 (양끝이 트랩이라 유한강성)",
        f"  ★★ λ_min = {D['lam_min']:.4e} N/m → τ_max = {D['tau_max'].to_compact():~.4gP}"
        f" = {float(D['tau_max']/D['tau_chain']):.2f} × τ_chain  ← **지배 척도**",
        f"  ★ 힘 정의를 섞으면 논문값과 정확히 2배 어긋난다"
        f" (scratch/chain_bend_from_papers.py §②에서 확인)",
        f"  이산 3점 굽힘 직접 풀이 = {Q(D['k_center_disc'], 'N/m').to('pN/um'):~.3fP}"
        f"  vs  빔 48EI/L³ = {D['kappa_center'].to('pN/um'):~.3fP}"
        f"  ({100*(D['k_center_disc']/float(D['kappa_center'].magnitude)-1):+.2f}%)",
        f"  λ_max = {D['lam_max']:.4e} N/m  (굽힘 {D['lam_bend']:.3e} · "
        f"신축 {D['lam_bond']:.3e} — 직선 사슬에서 두 블록은 분리된다)",
        f"  δ_max = M_c L²/(12EI) = {D['delta_max'].to('nm'):~.0fP}"
        f"   ℓ_k = {D['l_k'].to('nm'):~.2fP}   → 진폭 창 ℓ_k ≪ a < δ_max",
        f"  트랩 비드 = {D['trapped']}  ([P2] Fig.1A: 양끝 힘센서 + 중앙 구동)",
    ]
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}"
        f"  = {lg.ratio('times', 'dt', 'tau_B'):.3e} τ_B   (최속 굽힘 모드가 정함)",
        f"  ω       = {D['omega']:.0f} rad/s"
        f"  = {D['omega']/(2*math.pi):.1f} Hz   →  De = ω τ_max ="
        f" {lg.ratio('times', 'tau_max', 'tau_w'):.3f}"
        f"   (예전 정의 ω τ_chain = {lg.ratio('times', 'tau_chain', 'tau_w'):.3f})",
        f"  SNR     = |ŷ|/ℓ_k = {lg.ratio('lengths', 'y_resp', 'l_k'):.3f}"
        f"   (구동 a/ℓ_k = {lg.ratio('lengths', 'a', 'l_k'):.2f} — 분자를 응답으로 바꾼 값)",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}  ({N_CYCLES:g}주기)",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}   × n={sys_['n']}",
        f"  ⚠ 비용: prod {n_prod:,} 스텝. κ_θ = "
        f"{lg.ratio('energies', 'kappa_theta', 'kT'):.2e} kT 라 사슬이 뻣뻣해서"
        f" **측정과 무관한** 최속 모드가 dt를 정한다.",
        f"    비용 ∝ 1/ω 이므로 낮은 ω 쪽이 지배한다. 선택지 — (a) 감수"
        f" (b) κ₀를 낮춘다(계면활성제, [P2] Fig.4)",
        f"    (c) 높은 ω만 재고 낮은 ω는 준정적 극한으로 대체.",
    ]
    return inp, der, plan


def sweep_specs(dt_scale: float = 1.0, cycles: float = N_CYCLES, samples: int = 2000):
    """ω 스윕 7점의 스펙 문서를 **메모리에서** 만들어 돌려준다 (파일을 안 씀).

    ★ 하드 검사가 실패하면 `--spec` 은 `specs/` 에 아무것도 쓰지 않는다 (정상 동작, 규칙 2).
    그런데 검증·진단·시각화 스크립트는 파라미터가 필요하다 — 그것들이 스펙 파일에
    의존하면 "검사가 실패해서 그림도 못 그린다" 가 된다. 그래서 창구를 하나 둔다.
    """
    import argparse as _ap
    sys_ = load_system(ROOT / "intake" / "chain-bend-2d-oscill" / "system.yaml")
    args = _ap.Namespace(dt_scale=dt_scale, cycles=cycles, samples=samples)
    lo, hi = sys_["omega_range"]
    out = []
    for om in np.geomspace(lo, hi, N_SWEEP):
        _, spec, _, _, _, _ = build_spec(sys_, float(om), args)
        out.append({"params": spec.params, "numerics": spec.numerics})
    return sorted(out, key=lambda s: s["params"]["omega_star"])


# ════════════════════════════════════════════════════════════════════════
# main — L3 까지만
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, omega, args):
    lg = build_ledger(sys_, omega, dt_scale=args.dt_scale, n_cycles=args.cycles)
    D = lg.derived
    dt = lg.get("times", "dt")
    # 평형화: 구동 **전에** 사슬을 완화시킨다. ★★ 척도는 **τ_max**(최장 이완) 다.
    # 예전에는 20 τ_chain 을 썼는데 그건 2.2 τ_max 밖에 안 되어 잔여 과도가 11% 남았고,
    # K* 가 최대 21% 틀렸다 (실측: 평형화를 10 τ_max 로 늘리니 블록 산포가 1000배 감소).
    n_eq = int(round(20 * float((D["tau_max"] / dt).to(""))))
    n_prod = int(round(float((D["T_obs"] / dt).to(""))))
    sample_every = max(1, n_prod // args.samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(sys_, lg)      # ★ 한 번만 만든다 — 스펙과 리포트가
    tag = f"w{omega:.0f}"                          #    같은 객체를 봐야 갈라지지 않는다
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"
    if args.cycles != N_CYCLES:
        tag += f"-nc{args.cycles:g}"

    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_beads": sys_["n"], "trapped": D["trapped"],
                "L_chain_star": lg.ratio("lengths", "L_chain", "d"),
                # 무차원 강성은 전부 kT/d² 단위 (원장 에너지/kT 와 같은 수)
                "kappa_theta_star": lg.ratio("energies", "kappa_theta", "kT"),
                "k_t_star": lg.ratio("energies", "k_t_d2", "kT"),
                "k_bond_star": lg.ratio("energies", "k_b_d2", "kT"),
                "amp_star": lg.ratio("lengths", "a", "d"),
                "omega_star": float((Q(omega, "1/s") * D["tau_B"]).to("dimensionless").magnitude),
                "De": lg.ratio("times", "tau_max", "tau_w"),
                "De_chain_old": lg.ratio("times", "tau_chain", "tau_w"),
                "kappa_drive_star": lg.ratio("energies", "kappa_drive_d2", "kT"),
                "snr_response": lg.ratio("lengths", "y_resp", "l_k"),
                "M_c_star": lg.ratio("energies", "M_c", "kT"),
                "n_trapped": sys_["n_trapped"],
                "bending_impl": BENDING_IMPL},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "dt_over_tau_fast": args.dt_scale * C.GATE,
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": 20260804},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_eq, n_prod


def emit(sys_, omega, args) -> int:
    lg, spec, groups, checks, n_eq, n_prod = build_spec(sys_, omega, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 무결성 검사")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_eq, n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  ω={omega:.0f} rad/s   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)

    if spec.errors:
        print(f"\n❌ L3 무결성 오류 {len(spec.errors)}건 — 무차원화가 성립하지 않습니다.")
        return 1
    if verdict == "FAIL":
        print("\n❌ 하드 분리 검사 실패 — 스펙을 쓰지 않습니다.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 스펙: {p.relative_to(ROOT)}")
        return 0

    # ── L4 — 디스크의 스펙을 **되읽어서** 실행한다 (해시 검증이 그때 걸린다) ──
    outdir = ROOT / "runs" / run_id
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.txt").write_text(report)
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="L3 리포트만")
    ap.add_argument("--spec", action="store_true", help="L3 스펙 → specs/<run_id>.json")
    ap.add_argument("--omega", type=float, default=None,
                    help="구동 각진동수 [rad/s]. 기본은 스윕 최저값 (가장 비싼 점)")
    ap.add_argument("--sweep", action="store_true", help="ω 스윕 전체")
    ap.add_argument("--cycles", type=float, default=N_CYCLES, help="관측 주기 수")
    ap.add_argument("--force", action="store_true", help="완료된 런을 다시 실행")
    ap.add_argument("--dt-scale", type=float, default=1.0, help="dt 배율 (수렴 확인용)")
    ap.add_argument("--samples", type=int, default=2000, help="주기당 표본이 아니라 총 표본 수")
    ap.add_argument("--run", action="store_true",
                    help="L4 실행. ★ chain-bend 는 angle.Harmonic 힘 버그로 현재 거부됩니다")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/chain-bend-2d-oscill/system.yaml")
    lo, hi = sys_["omega_range"]

    if args.sweep:
        omegas = list(np.geomspace(lo, hi, N_SWEEP))
    else:
        omegas = [args.omega if args.omega is not None else lo]

    if not (args.report or args.spec or args.run):
        print("무엇을 할지 고르세요 — `--report` · `--spec` · `--run`")
        return 3

    rc = 0
    for i, om in enumerate(omegas):
        if i:
            print()
        rc |= emit(sys_, float(om), args)
    return rc




# ════════════════════════════════════════════════════════════════════════
# L4 — 스펙만 읽고 계를 세운다 (bdbot.run 이 돌리고, bdbot.health 가 판정한다)
#
# 구성은 `scratch/verify_chain_bend_gates.py` 의 관문 2종이 검증한 것을 그대로 씁니다:
#   · 트랩 = **유령 입자 + 조화 본드**(r0=0). 커스텀 힘이 필요 없습니다.
#   · 구동 = 중앙 유령을 `y = a sin(ωt)` 로 옮기는 updater
#   · 추정 = `bdbot.lockin` (관문 A 가 해석해와 3σ 이내로 대조)
#   · 페어 힘 없음 — κ_θ=1.4e6 kT 라 사슬이 자기접촉할 만큼 휘지 않습니다
#     (진폭 창 상한 δ_max 가 그것을 보장합니다). 대신 최소 비결합 거리를 감시합니다.
#
# ★★ 추정량에 **공칭 진폭을 쓰지 않습니다** — 유령 y 를 같이 재서 측정 위상자를 씁니다
#    (bdbot/lockin.py 도크스트링: 공칭을 쓰면 De=92 에서 부호까지 틀립니다).
# ════════════════════════════════════════════════════════════════════════
BENDING_IMPL = "custom_linear"   # "custom_linear" | "angle_harmonic"
# ★ 어떤 굽힘 구현을 쓰는가. `angle.Harmonic` 의 sinθ 클램프 제약은 **그 구현의 제약**이지
#   계의 제약이 아니므로, 하드 검사를 구현에 조건부로 겁니다 (지우지 않고 범위를 좁힘).
UPDATE_EVERY = 100      # 유령 이동 주기. ZOH 를 남겨 두고 측정 위상자로 상쇠한다


def make_bending_force(A, n_beads):
    """선형화 굽힘을 `md.force.Custom` 으로 직접 구현. **`angle.Harmonic` 을 못 쓰기 때문**.

    ★★ 왜 내장을 안 쓰는가 — `md.angle.Harmonic` 은 `sin θ` 를 SMALL=1.414e-3 으로
       클램프해서 그 아래에서 힘이 `sinθ/SMALL` 배로 축소된다 (힘 ∝ κ(θ−π)² — 선형이
       아니라 2차). **에너지는 0.000% 정확**해서 에너지 검증으로는 안 잡힌다.
       이 계는 23개 각도 **전부**가 깨진 영역이다 (max|θ−π| = 9.4e-4 < 1.414e-3).
       진폭을 키워 벗어날 수도 없다 — 각도 프로파일이 20배 퍼져 있고 δ_max 가 막는다.
       `angle.Table` 도 같은 문제, `CosineSquared` 는 θ₀=π 에서 4차라 배제.
       재현: `scratch/verify_angle_force_small_theta.py`

    구현하는 것은 L3 원장이 쓰는 **바로 그 선형화 형태**다:
        U = ½ κ_θ Σ_i θ_i²,  θ_i = (y_{i+1} − 2y_i + y_{i−1})/ℓ  →  F_y = −A y,  F_x = 0
    `A = κ_θ BᵀB` 는 `bending_matrix()` 그 자체이므로 **모델과 구현이 정확히 일치**한다
    (λ_max·τ_fast·driven_static_stiffness 가 전부 이 A 에서 나온 값이다).
    B 가 상수와 1차식을 소멸시키므로 **병진·(선형화)회전 불변**이다.

    ⚠️ 큰 변형에서는 진짜 각도 퍼텐셜과 다르다. 여기서는 max|θ| ~ 1e-3 이라 선형 영역
       깊숙이 있고, L3의 '소각 선형화' 검사가 그것을 게이트한다.
    """
    import hoomd.md as md

    class Bending(md.force.Custom):
        def __init__(self):
            super().__init__(aniso=False)
            self.A = np.ascontiguousarray(A, dtype=float)
            self.n = int(n_beads)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)   # ★ tag 인덱싱 필수
                pos = np.array(snap.particles.position, copy=True)
                m = tags < self.n                                # 유령 제외
                y = np.zeros(self.n)
                y[tags[m]] = pos[m, 1]
                fy = -(self.A @ y)
                arr.force[:] = 0.0
                arr.potential_energy[:] = 0.0
                arr.force[m, 1] = fy[tags[m]]
                # 총합이 ½yᵀAy 가 되도록 입자별로 배분 (배분 자체는 규약)
                arr.potential_energy[m] = -0.5 * y[tags[m]] * fy[tags[m]]

    return Bending()


def _move_ghost_action(ghost_tag, amp, omega, dt):
    import hoomd

    class MoveGhost(hoomd.custom.Action):
        def act(self, timestep):
            y = amp * math.sin(omega * timestep * dt)
            with self._state.cpu_local_snapshot as snap:
                tags = np.array(snap.particles.tag, copy=True)
                loc = np.flatnonzero(tags == ghost_tag)
                if len(loc):
                    snap.particles.position[loc[0], 1] = y
    return MoveGhost()


def _chain_frame(n, trapped, ell, L_chain):
    import gsd.hoomd

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:                                  # 유령을 비드 위에 겹쳐 둔다
        pos.append(list(pos[g]))
        typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.particles.mass = [1.0] * len(pos)                # BD 는 질량을 안 쓴다 (함정 5)
    f.configuration.box = [4.0 * L_chain] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(trapped)]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(trapped)
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])
    return f


def assert_angle_force_valid(params: dict) -> tuple[float, float]:
    """★★ `angle.Harmonic` 의 힘이 유효한 영역인지 확인하고, 아니면 거부한다.

    L3 하드 검사와 **같은 판정을 빌더에서 한 번 더** 한다 (방어 이중화). 이유:
    L3 검사는 스펙을 쓰지 않게 막고 `RUN.execute` 는 FAIL 스펙을 거부하지만, 빌더를
    **직접** 호출하는 경로(테스트·수동 호출·예전 스펙)는 아무것도 막지 못한다. 그리고
    이 실패는 크래시가 아니라 **조용히 틀린 물리**다 — 사슬이 실제보다 무르게 나온다.
    """
    n = int(params["n_beads"])
    ell = float(params["L_chain_star"]) / (n - 1)
    _, prof = driven_response(n, float(params["kappa_theta_star"]), ell,
                              float(params["k_t_star"]), 1.0,
                              float(params["omega_star"]), float(params["amp_star"]))
    th_hi, th_lo = response_angles(prof, ell)
    if th_lo < ANGLE_SIN_SMALL:
        raise ValueError(
            f"angle.Harmonic 의 힘이 유효하지 않은 영역입니다 — 실행을 거부합니다.\n"
            f"  min|θ−π| = {th_lo:.3e} < SMALL = {ANGLE_SIN_SMALL:.3e}  "
            f"(max|θ−π| = {th_hi:.3e})\n"
            f"  HOOMD 는 sin θ 를 SMALL 로 클램프하므로 그 아래에서 힘이 sinθ/SMALL 배로\n"
            f"  축소됩니다 (힘 ∝ κ(θ−π)² — 선형이 아니라 2차). **에너지는 정확합니다** →\n"
            f"  에너지로 검증하면 통과하고 힘만 틀립니다. 사슬이 실제보다 무르게 나옵니다.\n"
            f"  재현: scratch/verify_angle_force_small_theta.py (skill bd-hoomd 함정 15)\n"
            f"  길: ① force.Custom 으로 굽힘 직접 구현 (정확하지만 26배 느림)\n"
            f"      ② κ₀ 를 낮춰 사슬을 무르게 ([P2] 계면활성제) → θ 가 커진다\n"
            f"      ③ 이 영역은 선형이므로 MD 없이 해석적으로 푼다 (정확 최소화와 0.32%)")
    return th_hi, th_lo


@RUN.builder("chain-bend-2d-oscill")
def build(spec, outdir=None) -> RUN.Build:
    # ★ 방어 이중화 — 스펙이 선언한 구현과 이 빌더가 실제로 쓰는 구현이 어긋나면 거부.
    #   `angle.Harmonic` 스펙을 커스텀 힘으로 돌리면 조용히 다른 물리가 된다.
    _impl = spec.params.get("bending_impl", "angle_harmonic")
    if _impl != BENDING_IMPL:
        raise ValueError(f"스펙은 굽힘 구현 '{_impl}' 를 선언했는데 빌더는 "
                         f"'{BENDING_IMPL}' 입니다. 스펙을 다시 생성하세요.")
    if _impl == "angle_harmonic":
        assert_angle_force_valid(spec.params)
    # ★ `assert_angle_force_valid` 는 위에서 **구현이 angle_harmonic 일 때만** 부른다.
    #   여기 있던 무조건 호출은 제거했다 — 커스텀 힘에는 sinθ 클램프가 없으므로
    #   그 판정을 그대로 적용하면 유효한 구성을 거부하게 된다.
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_beads"])
    trapped = sorted(int(t) for t in P["trapped"])
    L_chain = float(P["L_chain_star"])
    ell = L_chain / (n - 1)
    k_t, k_b, kth = float(P["k_t_star"]), float(P["k_bond_star"]), float(P["kappa_theta_star"])
    amp, omega = float(P["amp_star"]), float(P["omega_star"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    mid = trapped[len(trapped) // 2]
    ghost_mid = n + trapped.index(mid)

    sim = SIM.make_sim(_chain_frame(n, trapped, ell, L_chain), seed=seed)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=k_b, r0=ell)
    bond.params["trap"] = dict(k=k_t, r0=0.0)          # 트랩 = 유령과의 조화 본드
    # ★ angle.Harmonic 은 이 영역에서 힘이 틀린다 (make_bending_force 도크스트링)
    angle = make_bending_force(bending_matrix(n, kth, ell), n)
    # ★ 유령은 적분하지 않는다 — updater 가 위치를 직접 쓴다. 비드만 BD.
    bd = md.methods.Brownian(filter=hoomd.filter.Type(["A"]), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=[bond, angle])
    integ.integrate_rotational_dof = False             # 함정 3
    sim.operations.integrator = integ
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=_move_ghost_action(ghost_mid, amp, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, int(Nm["n_prod"]) // 200))

    def pe_per_particle():
        return float(np.array(bond.energies).sum() + np.array(angle.energies).sum()) / n

    def sample(timestep, phase):
        snap = sim.state.get_snapshot()
        p = np.array(snap.particles.position, dtype=float)
        # ★ 유령 y 를 **같이** 잰다 — 이게 없으면 ZOH 가 K* 오차로 그대로 넘어간다
        nb = p[:n, :2]
        sep = np.linalg.norm(nb[:, None, :] - nb[None, :, :], axis=-1)
        np.fill_diagonal(sep, np.inf)
        for i in range(n - 1):                          # 결합 이웃은 제외 (겹침 감시용)
            sep[i, i + 1] = sep[i + 1, i] = np.inf
        return {"t": timestep * dt, "y_bead": float(p[mid, 1]),
                "y_ghost": float(p[ghost_mid, 1]),
                "y_end0": float(p[trapped[0], 1]), "y_end1": float(p[trapped[-1], 1]),
                "min_sep_nonbonded": float(sep.min()),
                "max_theta": float(np.abs(np.diff(p[:n, 1], n=2)).max() / ell)}

    def finalize(cols):
        t, yb, yg = cols["t"], cols["y_bead"], cols["y_ghost"]
        blocks_y = LI.lockin_blocks(t, yb, omega, n_blocks=10)
        blocks_g = LI.lockin_blocks(t, yg, omega, n_blocks=10)
        Kb = np.array([LI.k_star(a, b, k_t, omega) for a, b in zip(blocks_y, blocks_g)])
        K, Ksem = LI.agg(Kb)
        yh, _ = LI.agg(blocks_y)
        gh, _ = LI.agg(blocks_g)
        # 공칭 진폭으로도 계산해 **둘의 차이를 남긴다** (ZOH 함정의 사후 증거)
        K_nom, _ = LI.agg(np.array([LI.k_star(a, complex(amp), k_t, omega) for a in blocks_y]))
        # 정적 극한 예측 — **내 모델에서 유도**된 값이므로 implementation_check
        K_static = driven_static_stiffness(n, kth, ell, k_t)
        de = float(P["De"])
        obs = [
            MET.observable("K_prime", float(K.real), K_static if de < 1.5 else None,
                           "kT/d^2", role="implementation_check" if de < 1.5 else "measurement",
                           tol_pct=15.0,
                           note=f"저장 강성 K'(ω). De={de:.3f}. 정적 극한 예측 "
                                f"{K_static:.4g} (driven_static_stiffness — 트랩 유한강성 포함). "
                                f"De<1.5 에서만 대조한다"),
            MET.observable("K_doubleprime", float(K.imag), None, "kT/d^2",
                           role="measurement", note="손실 강성 K''(ω)"),
            MET.observable("K_sem", Ksem, None, "kT/d^2", role="measurement",
                           note="블록 10개 산포 (실·허 중 큰 쪽)"),
            MET.observable("y_response", float(abs(yh)), None, "d", role="measurement",
                           note=f"응답 진폭 |ŷ|. 구동 |ŷ_c|={abs(gh):.5g} (공칭 {amp:g})"),
        ]
        return {"observables": obs,
                "extra": {"K_prime": float(K.real), "K_doubleprime": float(K.imag),
                          "K_sem": Ksem, "K_static_pred": K_static,
                          "K_prime_nominal_amp": float(K_nom.real),
                          "nominal_vs_measured_pct":
                              100 * (K_nom.real - K.real) / abs(K.real) if K.real else None,
                          "y_resp_abs": float(abs(yh)), "drive_abs": float(abs(gh)),
                          "drive_over_nominal": float(abs(gh) / amp),
                          "De": de, "omega_star": omega,
                          "min_sep_nonbonded": float(cols["min_sep_nonbonded"].min()),
                          "max_theta": float(cols["max_theta"].max())}}

    every = max(1, int(Nm["sample_every"]))
    return RUN.Build(
        sim=sim, forces=[bond, angle], n_particles=n,
        sample=sample, pe_per_particle=pe_per_particle, sample_every=every,
        phases=[RUN.Phase("예열", int(Nm["n_eq"]), collect=False,
                          note="구동 ON · 표본 버림 (20 τ_max — 2.2 τ_max 는 11% 과도가 남았다)"),
                RUN.Phase("생산", int(Nm["n_prod"]), every,
                          note=f"락인 수집 · ω*={omega:.4g} · De={P['De']:.3f}")],
        tags=["2D", "chain", "bending", "angle_harmonic", "oscillatory_drive",
              "microrheology", "lockin", "newtonian"],
        physical={"n_beads": n, "De": float(P["De"]), "omega_star": omega,
                  "kappa_theta_star": kth, "k_t_star": k_t, "amp_star": amp,
                  **{k: v for k, v in spec.raw["back_transform"].items()
                     if isinstance(v, float)}},
        finalize=finalize)


if __name__ == "__main__":
    sys.exit(main())
