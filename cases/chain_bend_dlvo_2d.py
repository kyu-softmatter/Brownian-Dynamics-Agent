"""`chain-bend-2d-dlvo` — chain-bend-2d-oscill 의 대안가설 브랜치.

같은 3점 굽힘·오실레이션 기하(양끝 고정 트랩 + 중앙 트랩 진동)를 쓰되, 비드 결합을
[P1][P2]의 JKR 접착접촉+굽힘강성이 아니라 **DLVO 중심 페어 퍼텐셜의 2차극소**로
만든다. **명시적 각도(굽힘) 퍼텐셜은 없다 — 이게 가설이다.**

[P1] Pantina & Furst, PRL 94, 138301 (2005) p.2 좌단이 명시적으로 예측한다:
  "Such assumptions have been made based on DLVO interactions between particles,
  which are centro-symmetric. If particles did undergo free rotations, we would
  expect the aggregates to respond to the bending moment by forming a
  trianglelike structure..."
이 케이스는 그 예측을 실행으로 확인한다 (원칙 6·7 — 불일치는 발견이지 실패가 아니다).

★★ 구조적 사실 (실행 전 유도, `intake/chain-bend-2d-dlvo/observation.yaml` G1):
직선 사슬에서 순수 중심력 결합이 **자연장(natural length, U'=0)에 있으면** 횡변위
y에 대한 결합에너지 변화가 O(y⁴)이다 (O(y²) 항의 계수가 U'(ℓ)이고 힘이 0인 지점에선
0이기 때문) — 즉 **선형 굽힘강성이 정확히 0**이다. chain-bend-2d-oscill 의
bending_matrix() 같은 게 이 계에는 없다. 3점 굽힘에서 매끈한 빔 곡률 대신 [P1]이
말한 "삼각형형" 국소 꺾임이 나올 가능성이 높다 — 지어내지 않고 L4로 확인한다.

물리 파라미터 출처: [P1] d=1.47µm·ψ0=40mV, MgCl2 10mM(이온강도),
Hamaker A=1.05e-20 J(웹검색, tier 3 — 1차 출처 미확인, 불확실성 명시).
DLVO 곡선: HHF(등전위 약한중첩) EDL + Derjaguin 비지연 vdW.
계산: `scratch/dlvo_ledger.py` (장벽 382kT@0.49nm, 2차극소 −11.7kT@11.2nm).

퍼텐셜 구현: **본드 리스트가 아니라 전역 페어 퍼텐셜**이다 (`md.pair.Table`, 전체
이웃쌍에 적용 — 1-B `soft-r3-2d-A-sweep` 과 같은 패턴). "결합"은 인접 비드가 서로의
2차극소에 앉아서 생기는 것이지, 미리 정한 본드 토폴로지가 아니다 — 이게 "중심
페어 퍼텐셜만으로 사슬을 만든다"는 가설 그 자체다. 컷오프(h≲60nm)가 비인접 비드
간격(≈3µm)보다 훨씬 짧아 정상 상태에선 NN만 상호작용하지만, 사슬이 접히면 원거리
비드끼리도 자연스럽게 상호작용한다 — 이것도 실제 물리이므로 막지 않는다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/chain_bend_dlvo_2d.py --n 9 --omega 1000 --amp 100 --report
    $PY cases/chain_bend_dlvo_2d.py --n 9 --omega 1000 --amp 100 --spec
    $PY cases/chain_bend_dlvo_2d.py --n 9 --omega 1000 --amp 100 --run
"""
from __future__ import annotations

import argparse
import json
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
N_CYCLES = 50.0
NA = 6.02214076e23
E_CHARGE = 1.602176634e-19
EPS0 = 8.8541878128e-12
R_WCA = 2 ** (1 / 6)          # WCA 컷오프 (σ=d 단위) — bdbot/pairpot.py 와 동일 관례
CUTOFF_H_STAR = 0.06          # 표 컷오프 — 표면간극 h/d. 2차극소(≈0.0076)의 8배 밖


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
        "psi0": P(raw["particle"]["surface_potential"]),
        "n_list": [int(x) for x in raw["particle"]["count"]["value"]],
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "eps_r": P(raw["medium"]["relative_permittivity"]),
        "ionic_strength": P(raw["medium"]["ionic_strength"]),
        "A_H": P(con["hamaker_constant"]),
        # ★ 대조군 — `--jkr` 일 때만 쓰인다. interactions[1] 은 기본 OFF.
        "kappa_theta": P(raw["interactions"][1]["angle_stiffness"]),
        "k_t": P(raw["external"]["stiffness"]),
        "amp_range": [float(x) for x in raw["external"]["amplitude"]["value"]],
        "omega_range": [float(x) for x in raw["external"]["omega_range"]["value"]],
        "n_trapped": int(raw["external"]["n_trapped"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② DLVO(h) — 축약변수(h*=h/d) 3개(kappa_star, edl_amp, vdw_amp)로 결정된다.
#    ★ 감으로 쓰지 않는다 — scratch/dlvo_ledger.py 에서 SI로 먼저 계산해 검증한
#    바로 그 식(HHF 등전위 약한중첩 + Derjaguin 비지연 vdW)을 축약형으로 옮긴 것.
# ════════════════════════════════════════════════════════════════════════
def dlvo_reduced_params(sys_: dict) -> dict:
    d = sys_["d"].value.to("m").magnitude
    a = d / 2
    T = sys_["T"].value.to("K").magnitude
    kT = 1.380649e-23 * T
    eps_r = sys_["eps_r"].value.to("dimensionless").magnitude
    psi0 = sys_["psi0"].value.to("V").magnitude
    A_H = sys_["A_H"].value.to("J").magnitude
    c_salt = sys_["ionic_strength"].value.to("mol/m^3").magnitude   # MgCl2 몰농도
    I_SI = 0.5 * (c_salt * (2 ** 2) + 2 * c_salt * (1 ** 2)) * NA   # 1/m^3 (이온강도, MgCl2→Mg2++2Cl-)
    kappa = math.sqrt(2 * I_SI * E_CHARGE ** 2 / (EPS0 * eps_r * kT))
    a_star = 0.5
    return {
        "a_star": a_star,
        "kappa_star": kappa * d,
        "edl_amp": 2 * math.pi * EPS0 * eps_r * a * psi0 ** 2 / kT,   # U_edl* 앞 계수
        "vdw_amp": A_H / (12 * kT),                                    # U_vdw* = -vdw_amp*a_star/h*
        "kT": kT, "d": d,
    }


def U_star(h_star, p: dict):
    """U(h)/kT, h*=h/d>0. EDL(HHF 등전위,약한중첩) + vdW(Derjaguin 비지연)."""
    h_star = np.asarray(h_star, dtype=float)
    u_edl = p["edl_amp"] * np.log1p(np.exp(-p["kappa_star"] * h_star))
    u_vdw = -p["vdw_amp"] * p["a_star"] / h_star
    return u_edl + u_vdw


def F_h_star(h_star, p: dict):
    """-dU*/dh* (h* 증가 방향의 힘. 양수=반발). 해석적 도함수."""
    h_star = np.asarray(h_star, dtype=float)
    k = p["kappa_star"]
    dU_edl = p["edl_amp"] * (-k) * np.exp(-k * h_star) / (1 + np.exp(-k * h_star))
    dU_vdw = p["vdw_amp"] * p["a_star"] / h_star ** 2
    return -(dU_edl + dU_vdw)


def find_well(p: dict) -> dict:
    """장벽·2차극소 위치+깊이, 우물 곡률(=결합 방사강성, kT/d² 단위)을 이분법+중심차분으로."""
    hs = np.geomspace(1e-4, CUTOFF_H_STAR * 3, 4000)
    Us = U_star(hs, p)
    ibar = int(np.argmax(Us[: int(len(hs) * 0.5)]))
    barrier_h, barrier_U = float(hs[ibar]), float(Us[ibar])
    tail = Us[ibar:]
    imin = int(np.argmin(tail))
    h_min, U_min = float(hs[ibar + imin]), float(tail[imin])
    dh = h_min * 1e-4
    k_bond_star = (U_star(h_min + dh, p) - 2 * U_star(h_min, p) + U_star(h_min - dh, p)) / dh ** 2
    return {"barrier_h": barrier_h, "barrier_U": float(barrier_U), "h_min": h_min,
            "U_min": float(U_min), "k_bond_star": float(k_bond_star)}


def trapped_indices(n: int) -> list[int]:
    return [0, n // 2, n - 1]


def _bow(y: np.ndarray) -> float:
    """굽음 — 양끝을 잇는 직선 대비 최대 이탈 [d]. 강체 평행이동·기울기를 뺀 순수 굽힘.

    ★ 왜 이게 필요한가: 이 계의 트랩은 유한강성(k_t)이라 사슬이 통째로 밀린다. y 변위를
      그대로 재면 그 강체 운동이 굽힘 신호를 덮는다. 실측 (n=9, ω=3000, a=632nm):
          y 전범위   DLVO 0.303 d  vs JKR 0.090 d  →  3.4배
          굽음      DLVO 0.1175 d vs JKR 0.0060 d →  **19.6배**
      같은 데이터인데 판별력이 5.8배 차이난다.
    """
    n = len(y)
    base = np.linspace(0.0, 1.0, n) * (y[-1] - y[0]) + y[0]
    return float(np.abs(y - base).max())


# ════════════════════════════════════════════════════════════════════════
# ③ 스케일 원장
# ════════════════════════════════════════════════════════════════════════
def bending_matrix(n: int, kappa_theta: float, ell: float) -> np.ndarray:
    """U = ½κ_θ Σθ_i², θ_i=(y_{i+1}−2y_i+y_{i−1})/ℓ 의 2차형식 A = κ_θ BᵀB.

    ★ chain-bend-2d-oscill 과 **같은 식**이다 (거기서 논문 κ₀ 재현·이산↔연속 매핑이
      검증됐다). 여기서는 대조군(`--jkr`)에서만 쓴다.
    """
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def build_ledger(sys_, n: int, omega: float, amp_nm: float, *, dt_scale=1.0,
                  n_cycles=N_CYCLES, jkr: bool = False,
                  kt_scale: float = 1.0) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B = b["kT"], b["gamma"], b["tau_B"]
    # ★ kt_scale — 트랩을 의도적으로 세게 해서 **위치 제어 조건**으로 밀어넣는다.
    #   왜 필요한가: 기본 k_t(스케치 tier 0, 10 pN/µm)로는 구동 비드가 트랩 위치를
    #   전혀 못 따라간다 (실측 추종률 |ŷ|/a — JKR 3.5%, DLVO 29%). 사슬·점성이
    #   트랩보다 강해서 트랩이 늘어날 뿐 비드가 안 끌려온다.
    #   실측 설계표 (n=9, ω=3000, JKR):
    #       배율     1   10   30  100  300  1000
    #       추종률 0.035 0.33 0.59 0.82 0.93 0.98
    #       dt비용  1.00 1.00 1.00 1.01 1.02 1.08 (JKR) / 1.0 1.0 - 2.0 6.0 20 (DLVO)
    #   ★ 100배를 쓴다: 추종 0.82 인데 비용이 거의 없다. 300배부터는 K′ 추정의
    #     차이 소거가 위험해진다 — K′ = k_t(ŷ_c/ŷ − 1) 인데 (ŷ_c/ŷ − 1) 이 0.09까지
    #     줄어 상대 잡음이 증폭된다 (1000배면 0.027).
    #   ⚠️ 물리적 정직성: 100배 = 1000 pN/µm 로, 실제 광집게([P1] ~40 pN/µm)보다
    #     25배 세다. 이건 **위치 제어 극한의 수치 실험**이지 스케치의 계가 아니다.
    #     (trap-drag 가 '뻣뻣한 트랩 = 위치 제어 조건'이라고 적어둔 것과 같은 맥락)
    k_t = sys_["k_t"].value.to("N/m") * float(kt_scale)
    amp = Q(amp_nm, "nm").to("m")

    p = dlvo_reduced_params(sys_)
    w = find_well(p)
    ell = d * (1 + w["h_min"])                       # 결합 자연길이 (중심간, 자연장=U'=0)
    L_chain = (n - 1) * ell
    k_bond = Q(w["k_bond_star"], "dimensionless") * kT / d ** 2   # 방사(신축) 강성
    sigma_bond = (kT / k_bond) ** 0.5

    tau_bond = C.relaxation_time(gamma, k_bond)      # DLVO 결합 신축
    tau_k = C.relaxation_time(gamma, k_t)

    # ★ 대조군(`--jkr`): 굽힘항을 켜면 강성행렬의 최대고유값이 dt 를 정할 수도 있다.
    #   chain-bend-2d-oscill 에서 실제로 그랬다(τ_fast=0.28µs 가 지배) — 여기서도
    #   **더 빠른 쪽**으로 dt 를 잡는다. 안 그러면 조용히 발산한다.
    kappa_theta = None
    tau_bend_fast = None
    if jkr:
        kappa_theta = sys_["kappa_theta"].value.to("J")
        kth_star = float((kappa_theta / kT).to("dimensionless").magnitude)
        A_bend = bending_matrix(n, kth_star, float((ell / d).to("dimensionless").magnitude))
        for i in trapped_indices(n):
            A_bend[i, i] += float((k_t * d ** 2 / kT).to("dimensionless").magnitude)
        lam_max_star = float(np.linalg.eigvalsh(A_bend)[-1])          # kT/d² 단위
        tau_bend_fast = C.relaxation_time(gamma, Q(lam_max_star, "dimensionless") * kT / d ** 2)
    # ★★ dt 는 **가장 빠른 모드**가 정한다 — 후보 셋을 전부 넣는다.
    #   구멍이었다: 예전에는 τ_k(트랩)를 빼고 τ_bond·τ_bend 만 봤다. --kt-scale 로
    #   트랩을 세게 하면 k_t 가 k_bond 를 넘는 순간(배율 200 부근, k_bond*=1.04e6)
    #   트랩이 최속 모드가 되는데 dt 는 그대로여서 **조용히 부족해진다**.
    #   (DLVO 브랜치는 굽힘행렬이 없어 τ_bend 가 None 이라 특히 위험했다)
    cands = [tau_bond, tau_k] + ([tau_bend_fast] if tau_bend_fast is not None else [])
    tau_fast = min(cands, key=lambda q: float(q.to("s").magnitude))
    dt = dt_scale * C.dt_from_gate(tau_fast)

    # ★★ 굽힘(횡) 선형강성은 구조적으로 0 이다 (자연장 평형 + 중심력) — 지어낸 대체
    #    척도를 쓰지 않는다. 대신 사슬 전체의 형태 이완을 "윤곽길이 확산시간"으로
    #    거칠게 어림한다 (Rouse류 정확한 스펙트럼이 아니라 ★제안 상한선).
    tau_chain_diffusion = (L_chain ** 2 / b["D_t"]).to("s")

    tau_w = Q(1.0 / omega, "s")
    tau_period = Q(2 * math.pi / omega, "s")
    T_obs = (n_cycles * tau_period).to("s")

    lg = SC.ScaleLedger()
    lg.add_length("sigma_bond", sigma_bond.to("m"), "결합 방사 열요동 폭 √(kT/k_bond)", star=True)
    lg.add_length("h_min", Q(w["h_min"], "dimensionless") * d, "2차극소 위치(표면간극)", star=True)
    lg.add_length("a", amp, "구동 진폭", star=True)
    lg.add_length("d", d, "비드 지름")
    lg.add_length("ell", ell.to("m"), "결합 자연길이 (중심간, d+h_min)")
    lg.add_length("L_chain", L_chain.to("m"), "사슬 윤곽길이 (n-1)ell")
    lg.add_time("tau_p", b["tau_p"], "m/γ 관성 이완", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_bond", tau_bond, "★ γ/k_bond 결합 신축 — 유일한 선형 강성 모드. dt를 정한다",
                star=True)
    lg.add_time("tau_k", tau_k, "γ/k_t 트랩")
    lg.add_time("tau_w", tau_w, f"1/ω 구동 (ω={omega:.0f} rad/s)")
    lg.add_time("tau_period", tau_period, "2π/ω 구동 주기")
    lg.add_time("tau_chain_diff", tau_chain_diffusion,
                "★제안: L_chain²/D_t — 사슬 형태 이완의 거친 상한 (굽힘강성이 0이라 "
                "정확한 Rouse 스펙트럼 유도는 별도 검증 필요, 지금은 어림)", star=True)
    lg.add_time("tau_B", tau_B, "d²/D_t 확산 (기준)")
    lg.add_time("T_obs", T_obs, f"관측창 ({n_cycles:g}주기)", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.add_energy("k_t_d2", (k_t * d ** 2).to("J"), "k_t d² 트랩 강성")
    lg.add_energy("k_bond_d2", (k_bond * d ** 2).to("J"), "k_bond d² 결합 방사강성", star=True)
    lg.add_energy("well_depth", Q(-w["U_min"], "dimensionless") * kT,
                  "|2차극소 깊이| — 결합이 가역적일 수 있는 열에너지 규모", star=True)
    lg.add_time("tau_fast", tau_fast, "dt를 정하는 최속 모드", role="", star=True)
    lg.declare_absent(
        "box",
        "주기경계 없음 (사슬 하나, 트랩이 위치를 고정). chain-bend-2d-oscill 과 같은 사유.")
    if jkr:
        lg.add_energy("kappa_theta", kappa_theta,
                      "★ 대조군: JKR 접선 굽힘강성 κ_θ = EI/ℓ ([P1][P2])", star=True)
        lg.add_time("tau_bend_fast", tau_bend_fast,
                    "★ 굽힘 강성행렬 최대고유값 기준 최속 모드 — dt를 정할 수 있다", star=True)
    else:
        lg.declare_absent(
            "bending_stiffness",
            "★★ 구조적으로 없다 (지어낸 대체값을 넣지 않는다) — 직선 사슬 + 순수 중심력 + "
            "결합이 자연장(U'=0)에 있으면 횡변위에 대한 선형(O(y²)) 복원력이 정확히 0이다 "
            "(O(y²) 계수가 U'(ℓ)이고 자연장에서 0이기 때문). chain-bend-2d-oscill 의 "
            "bending_matrix/λ_min 에 대응하는 것이 이 계엔 없다 — 이것 자체가 가설(G1)이고 "
            "L4 궤적의 형태(매끈한 곡률 vs 국소 꺾임)로 확인한다. "
            "★ `--jkr` 로 켜면 대조군이 된다.")
    lg.derived = dict(gamma=gamma, D_t=b["D_t"], m=b["m"], kT=kT, d=d, tau_B=tau_B,
                      ell=ell.to("m"), L_chain=L_chain.to("m"), k_bond=k_bond,
                      sigma_bond=sigma_bond.to("m"), tau_bond=tau_bond, tau_k=tau_k,
                      tau_w=tau_w, tau_period=tau_period, tau_chain_diff=tau_chain_diffusion,
                      dt=dt, T_obs=T_obs, omega=omega, n=n, amp=amp,
                      jkr=jkr, kappa_theta=kappa_theta, tau_fast=tau_fast,
                      tau_bend_fast=tau_bend_fast,
                      k_bond_star=w["k_bond_star"], h_min_star=w["h_min"],
                      U_min_star=w["U_min"], barrier_star=w["barrier_U"],
                      reduced=p, trapped=trapped_indices(n))
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ 이 계는 chain-bend-2d-oscill 과 달리 굽힘 선형강성이 "
        "구조적으로 없다 — dt를 정하는 척도(tau_bond, 결합 신축)와 재려는 대상(사슬 형태 "
        "이완, tau_chain_diff)의 관계가 사전에 확정되지 않는다. L4 실행으로 확인한다.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ④ 무차원수 + 검사
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg, n):
    D = lg.derived
    r = lg.ratio

    groups = [
        ND.Group("well_depth/kT", r("energies", "well_depth", "kT"),
                 ("energies", "well_depth"), ("energies", "kT"), "",
                 "★ 결합 깊이 — 열에너지 규모면 가역적(깨지기 쉬움)"),
        ND.Group("k_bond_star", r("energies", "k_bond_d2", "kT"),
                 ("energies", "k_bond_d2"), ("energies", "kT"), "k_bond d²/kT",
                 "결합 방사(신축) 강성 — 굽힘강성 아님"),
        ND.Group("k_t/k_bond", r("energies", "k_t_d2", "k_bond_d2"),
                 ("energies", "k_t_d2"), ("energies", "k_bond_d2"), "",
                 "트랩 vs 결합 강성"),
        ND.Group("sigma_bond/h_min", r("lengths", "sigma_bond", "h_min"),
                 ("lengths", "sigma_bond"), ("lengths", "h_min"), "",
                 "결합 열요동 폭 / 우물 위치 — 1에 가까우면 결합이 불안정"),
        ND.Group("a/sigma_bond", r("lengths", "a", "sigma_bond"),
                 ("lengths", "a"), ("lengths", "sigma_bond"), "",
                 "★ 구동 진폭 vs 결합 열요동 — SNR 대용"),
        ND.Group("a/L_chain", r("lengths", "a", "L_chain"),
                 ("lengths", "a"), ("lengths", "L_chain"), "", "진폭 vs 사슬 길이"),
        ND.Group("n_beads", float(n), None, None, "", "비드 수 (입력, ★제안 스윕)"),
        ND.Group("De_bond", r("times", "tau_bond", "tau_w"),
                 ("times", "tau_bond"), ("times", "tau_w"), "ω τ_bond",
                 "결합 신축 기준 Deborah (매우 작을 것 — τ_bond 가 극히 빠름)"),
        ND.Group("De_chain_diff", r("times", "tau_chain_diff", "tau_w"),
                 ("times", "tau_chain_diff"), ("times", "tau_w"), "ω τ_chain_diff",
                 "★제안: 사슬 형태이완 기준 Deborah (어림)"),
        ND.Group("De_trap", r("times", "tau_k", "tau_w"),
                 ("times", "tau_k"), ("times", "tau_w"), "ω τ_k", "트랩 기준"),
        ND.Group("tau_bond/tau_chain_diff", r("times", "tau_bond", "tau_chain_diff"),
                 ("times", "tau_bond"), ("times", "tau_chain_diff"), "",
                 "★ 척도 분리 폭 (dt를 정하는 모드 vs 형태 이완 어림)"),
        ND.Group("dt/tau_bond", r("times", "dt", "tau_bond"),
                 ("times", "dt"), ("times", "tau_bond"), "", "적분 해상"),
        ND.Group("n_cycles", r("times", "T_obs", "tau_period"),
                 ("times", "T_obs"), ("times", "tau_period"), "", "관측 주기 수"),
        ND.Group("St", r("times", "tau_p", "tau_B"),
                 ("times", "tau_p"), ("times", "tau_B"), "tau_p/tau_B", "관성 vs 확산"),
    ]
    checks = [
        C.Check("model", "참고: τ_p/τ_bond", r("times", "tau_p", "tau_bond"),
                C.GATE, "<=",
                "★ 결합이 깊고 좁은 우물이라 τ_bond 가 τ_p 에 근접한다 (관성 무시 기준을 "
                "~2.8배 넘김, n·ω·amp 와 무관 — 결합 물리 자체의 성질). "
                "chain-bend-2d-oscill 에도 같은 종류의 위반(ζ=0.65, 이쪽보다 훨씬 심함)이 "
                "있었고 OverdampedViscous vs Langevin(kT=0) 대조로 무해함(0.16%)이 검증됐다 "
                "— 이 계는 **아직 그 대조를 하지 않았다**. 미검증 상태로 소프트 처리하고 "
                "구조 스모크테스트 이후 필요하면 검증한다", hard=False),
        C.Check("model", "결합 안정     σ_bond/h_min", r("lengths", "sigma_bond", "h_min"),
                0.5, "<=",
                "★★ 결합 열요동 폭이 우물 위치의 절반을 넘으면 열적으로 자꾸 깨질 "
                "정도로 불안정하다는 뜻 — 그 자체가 결과일 수 있으나 사전에 표시",
                hard=False),
        C.Check("integration", "결합 신축 해상 dt/τ_bond", r("times", "dt", "tau_bond"),
                C.GATE, "<=", "DLVO 결합 신축 모드. 못 맞추면 발산"),
        C.Check("integration", "fastest mode resolved dt/tau_fast", r("times", "dt", "tau_fast"),
                C.GATE, "<=",
                "★ 대조군(--jkr)에서는 굽힘 강성행렬의 최대고유값이 더 빠를 수 있다 — "
                "그쪽으로 dt 를 잡는다 (chain-bend-2d-oscill 에서 실제로 그랬다)"),
        C.Check("integration", "drive resolved       dt/tau_w", r("times", "dt", "tau_w"), C.GATE, "<=",
                f"구동 ω = {lg.derived['omega']:.0f} rad/s 를 해상"),
        C.Check("statistics", "SNR(결합)     a/σ_bond", r("lengths", "a", "sigma_bond"), 3.0, ">=",
                "구동 진폭이 결합 열요동보다 커야 신호가 잡음 위로 나온다", hard=False),
        C.Check("statistics", "cycles observed      T_obs/(2pi/w)", r("times", "T_obs", "tau_period"),
                N_CYCLES, ">=", "위상 평균에 쓸 주기 수", hard=False),
    ]
    return groups, checks


def report_blocks(sys_, lg, n_eq, n_prod, n):
    D = lg.derived
    inp = [R.kv("d", f"{sys_['d'].value:~.4gP}", sys_["d"].tier, sys_["d"].source[:44]),
           R.kv("psi0", f"{sys_['psi0'].value:~.4gP}", sys_["psi0"].tier, sys_["psi0"].source[:44]),
           R.kv("I(MgCl2)", f"{sys_['ionic_strength'].value:~.4gP}",
                sys_["ionic_strength"].tier, sys_["ionic_strength"].source[:44]),
           R.kv("A_H", f"{sys_['A_H'].value:~.4gP}", sys_["A_H"].tier, sys_["A_H"].source[:44]),
           R.kv("n", f"{n}", 3, "★제안 스윕"),
           R.kv("omega", f"{D['omega']:.0f} rad/s", 3, "★제안 스윕"),
           R.kv("amp", f"{D['amp'].to('nm'):~.1fP}", 3, "★제안 스윕")]
    der = [
        f"  λ_D(반지름 무관) 는 build_ledger 밖 scratch/dlvo_ledger.py 참조",
        f"  결합: 장벽 {D['barrier_star']:.2f} kT   2차극소 {D['U_min_star']:.3f} kT"
        f" @ h={D['h_min_star']*float(D['d'].to('nm').magnitude):.2f} nm",
        f"  k_bond = {D['k_bond'].to('pN/um'):~.4fP} = {D['k_bond_star']:.4e} kT/d²"
        f"   σ_bond = {D['sigma_bond'].to('nm'):~.3fP}",
        f"  ell(자연장) = {D['ell'].to('nm'):~.2fP}   L_chain = {D['L_chain'].to('um'):~.3fP}",
        f"  ★★ 굽힘 선형강성 = 0 (구조적, declare_absent) — chain-bend-2d-oscill 의 "
        f"λ_min 에 대응하는 값이 이 계엔 없다",
        f"  τ_bond = {D['tau_bond'].to_compact():~.4gP}  (dt를 정하는 유일한 모드)",
        f"  τ_chain_diff(어림) = {D['tau_chain_diff'].to_compact():~.4gP}"
        f" = {float(D['tau_chain_diff']/D['tau_bond']):.3e} × τ_bond",
    ]
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}  = {lg.ratio('times','dt','tau_B'):.3e} τ_B",
        f"  ω       = {D['omega']:.0f} rad/s  →  De_bond = {lg.ratio('times','tau_bond','tau_w'):.3e}"
        f"   De_chain_diff(어림) = {lg.ratio('times','tau_chain_diff','tau_w'):.3f}",
        f"  SNR     = a/σ_bond = {lg.ratio('lengths','a','sigma_bond'):.3f}",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}  ({N_CYCLES:g}주기)",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}   × n={n}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# main — L3
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, n, omega, amp_nm, args):
    lg = build_ledger(sys_, n, omega, amp_nm, dt_scale=args.dt_scale, n_cycles=args.cycles,
                      jkr=args.jkr, kt_scale=args.kt_scale)
    D = lg.derived
    dt = lg.get("times", "dt")
    # ★★ τ_chain_diff(어림) 기준 평형화(5×)는 n_eq 가 수백억 스텝으로 폭발한다
    #   (τ_chain_diff/τ_bond ~ 1e7~1e8). 굽힘강성이 없어 "완전 평형화"의 의미 자체가
    #   불분명하고, 지금은 **구조 확인이 목표**(smoke)이지 통계 수렴이 아니다 —
    #   그래서 국소(결합) 이완 배수로 저렴하게 잡는다. 전체 형태 이완까지 재려면
    #   --eq-scale 로 늘려 별도 검증한다 (사용자 확인 2026-08-05: 스모크테스트 우선).
    n_eq = int(round(args.eq_scale * float((D["tau_fast"] / dt).to(""))))
    n_prod = int(round(float((D["T_obs"] / dt).to(""))))
    sample_every = max(1, n_prod // args.samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(sys_, lg, n)
    tag = f"n{n}-w{omega:.0f}-a{amp_nm:.0f}"
    if args.jkr:
        tag += "-jkr"
    if args.kt_scale != 1.0:
        tag += f"-kt{args.kt_scale:g}"
    if args.drive_mode != "trap":
        tag += f"-{args.drive_mode}"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"

    p = D["reduced"]
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_beads": n, "trapped": D["trapped"],
                "kappa_star": p["kappa_star"], "edl_amp": p["edl_amp"], "vdw_amp": p["vdw_amp"],
                "a_star": p["a_star"], "cutoff_h_star": CUTOFF_H_STAR,
                "k_t_star": lg.ratio("energies", "k_t_d2", "kT"),
                "amp_star": lg.ratio("lengths", "a", "d"),
                "omega_star": float((Q(omega, "1/s") * D["tau_B"]).to("dimensionless").magnitude),
                "De_bond": lg.ratio("times", "tau_bond", "tau_w"),
                "well_depth_star": D["U_min_star"], "h_min_star": D["h_min_star"],
                "k_bond_star": D["k_bond_star"],
                # ★ 대조군 스위치 — run_id 해시에 들어가야 두 브랜치가 갈린다
                "jkr": bool(args.jkr),
                "drive_mode": args.drive_mode,
                "kappa_theta_star": (lg.ratio("energies", "kappa_theta", "kT")
                                     if args.jkr else 0.0),
                "n_trapped": sys_["n_trapped"]},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": args.seed},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_eq, n_prod


def emit(sys_, n, omega, amp_nm, args) -> int:
    lg, spec, groups, checks, n_eq, n_prod = build_spec(sys_, n, omega, amp_nm, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 무결성 검사")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_eq, n_prod, n)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  n={n} ω={omega:.0f} a={amp_nm:.0f}nm"
              f"   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)

    if spec.errors:
        print(f"\n❌ L3 무결성 오류 {len(spec.errors)}건.")
        return 1
    if verdict == "FAIL":
        print("\n❌ 하드 분리 검사 실패 — 스펙을 쓰지 않습니다.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 스펙: {p.relative_to(ROOT)}")
        return 0

    outdir = ROOT / "runs" / run_id
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.txt").write_text(report)
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    verdict_txt = RUN.render_verdict(v)
    print(verdict_txt)

    # ★★ `result.txt` 를 반드시 쓴다 — 프로젝트의 **완료 마커**다 (다른 케이스 3종이
    #   전부 이걸 쓴다: soft_r3_2d·abp_rod_2d·trap_2d_5um). `bdbot.run.execute` 가
    #   써주는 게 아니라 **케이스 스크립트의 책임**이다. 빠뜨리면 조용히 셋이 깨진다:
    #     ① `bdbot.cli status` 가 런을 0개로 센다 (cli.py:133 이 result.txt 를 본다)
    #     ② `runid.prepare_outdir` 가 완료를 못 알아봐 같은 런을 계속 재실행한다
    #     ③ "미완료 정리" 스크립트가 **완료 런을 지운다** — 이 세션에서 실제로 6개 날렸다
    #   (2026-08-06 에 뒤늦게 추가. 그 전 137런은 scratch/backfill_result_txt.py 로 채웠다)
    if v["status"] != "skipped":
        obs_lines = []
        try:
            mj = json.loads((outdir / "metrics.json").read_text())
            for o in mj.get("observables", []):
                m = o.get("measured")
                p_ = o.get("predicted")
                tail = f"   (예측 {p_:.6g})" if isinstance(p_, (int, float)) else ""
                obs_lines.append(f"  {o['name']:<22} {m:.6g}{tail}" if m is not None
                                 else f"  {o['name']:<22} —")
        except Exception as e:                      # 결과를 못 읽어도 마커는 남긴다
            obs_lines.append(f"  (metrics.json 을 읽지 못함: {e})")
        result = "\n".join(["=" * 84, f"결과 — {run_id}", "=" * 84,
                            *obs_lines, "=" * 84, verdict_txt])
        (outdir / "result.txt").write_text(report + "\n" + result)
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n", type=int, default=None, help="비드 수 (기본: system.yaml 목록 전부)")
    ap.add_argument("--omega", type=float, default=None, help="rad/s (기본: 범위 최저값)")
    ap.add_argument("--amp", type=float, default=None, help="nm (기본: 범위 중앙값)")
    ap.add_argument("--cycles", type=float, default=N_CYCLES)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dt-scale", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eq-scale", type=float, default=200.0,
                    help="평형화 = 이 값 × τ_fast/dt (기본 200 — 국소 이완만. "
                         "★ 전체 형태 이완(τ_chain_diff)까지는 안 덮는다, 구조 스모크용")
    ap.add_argument("--drive-mode", choices=("trap", "position"), default="trap",
                    help="trap=트랩 중심을 움직임(실험적, 컴플라이언스 있음) / "
                         "position=구동 비드 y 를 직접 강제(변형률 제어, 유변학 표준). "
                         "position 은 K_transfer 를 낸다 — trap 의 K′ 과 다른 양")
    ap.add_argument("--kt-scale", type=float, default=1.0,
                    help="트랩 강성 배율. ★ 기본 1 이면 구동 비드가 트랩을 3.5~29%%밖에 "
                         "못 따라간다. 100 이면 82%% (비용 거의 없음). build_ledger 도크스트링에 설계표")
    ap.add_argument("--jkr", action="store_true",
                    help="★ 대조군: JKR 접선 굽힘강성(κ_θ=EI/ℓ)을 DLVO 위에 추가. "
                         "기하·트랩·DLVO·시드를 그대로 두고 굽힘항만 켜서 직접 대조한다")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/chain-bend-2d-dlvo/system.yaml")
    ns = [args.n] if args.n is not None else sys_["n_list"]
    omega = args.omega if args.omega is not None else sys_["omega_range"][0]
    amp = args.amp if args.amp is not None else math.sqrt(
        sys_["amp_range"][0] * sys_["amp_range"][1])

    if not (args.report or args.spec or args.run):
        print("무엇을 할지 고르세요 — `--report` · `--spec` · `--run`")
        return 3

    rc = 0
    for i, n in enumerate(ns):
        if i:
            print()
        rc |= emit(sys_, n, omega, amp, args)
    return rc



# ════════════════════════════════════════════════════════════════════════
# L4 — HOOMD 빌더
#
# ★★ WCA 코어를 particle diameter(σ=d)로 쓰면 안 된다 — WCA 는 r < 2^(1/6)σ 전체에서
#    반발이라, σ=d 로 두면 반발이 r=1.122d(= h=180nm) 까지 침범해 우리 2차극소
#    (h=11.16nm≪180nm)를 완전히 짓밟는다. 대신 WCA 컷오프가 **정확히 r=d 에서 끝나도록**
#    σ_c = d·2^(-1/6) 로 잡는다 — 그러면 WCA 는 r<d(입자 겹침) 에서만 반발하고 r≥d
#    (표면간극 h≥0, DLVO 가 사는 영역) 에서는 정확히 0 이다. 겹침도 이중계산도 없다.
# ════════════════════════════════════════════════════════════════════════
SIGMA_CORE_STAR = 2 ** (-1.0 / 6.0)     # WCA 컷오프가 r*=1(=d, 표면 접촉)에서 끝나도록


def build_table_arrays(P: dict, r_min_star: float, r_cut_star: float, nbins: int = 8000):
    """`md.pair.Table` 용 (U, F) 배열. r* = 1+h*, h*≥0. endpoint=False (bd-hoomd 함정10)."""
    r = np.linspace(r_min_star, r_cut_star, nbins, endpoint=False)
    h = r - 1.0
    h = np.maximum(h, 1e-6)             # r_min_star 가 1+eps 이상이라 실질적으로 무관
    U = U_star(h, P)
    F = F_h_star(h, P)                  # r* 증가 = h* 증가 (같은 방향) → 그대로 쓴다
    U_cut = float(U_star(max(r_cut_star - 1.0, 1e-6), P))
    U = U - U_cut                       # 컷오프에서 0 (bd-hoomd 스니펫 관례)
    return r, U, F


def make_frame(n: int, ell_star: float, trapped: list[int], box_star: float,
               skip_trap: int | None = None):
    """skip_trap: 이 비드에는 트랩 본드를 걸지 않는다 (position 구동 모드의 구동 비드)."""
    import gsd.hoomd
    pos = [[(i - (n - 1) / 2) * ell_star, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:
        pos.append(list(pos[g]))
        typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.particles.mass = [1.0] * len(pos)
    f.configuration.box = [box_star] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    # 트랩 본드만. 백본 본드 없음(가설). ★ position 구동에서는 구동 비드가 적분 대상이
    # 아니라 트랩이 물리적으로 무의미하고, 늘어난 본드가 pe_per_particle 을 오염시킨다.
    grp = [[g, n + j] for j, g in enumerate(trapped) if g != skip_trap]
    f.bonds.N = len(grp)
    f.bonds.types = ["trap"]
    f.bonds.typeid = [0] * len(grp)
    f.bonds.group = np.array(grp)
    return f


def make_bending_force(A, n_beads):
    """선형화 굽힘 U = ½κ_θΣθ_i² 을 `md.force.Custom` 으로 직접 구현 (F_y = −A y).

    ★★ `md.angle.Harmonic` 을 쓰지 않는 이유 — sin θ 를 SMALL=1.414e-3 으로 클램프해서
       거의 곧은 사슬의 **힘만** 축소된다 (에너지는 0.000% 정확 → 에너지 검증으로 안 잡힘).
       chain-bend-2d-oscill 에서 실측 규명됐고 (bd-hoomd 함정 15), 이 계도 같은 영역이다.
       거기서 검증된 구현을 그대로 가져왔다 — 모델(bending_matrix)과 구현이 정확히 일치한다.
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


UPDATE_EVERY = 50


@RUN.builder("chain-bend-2d-dlvo")
def build(spec, outdir=None) -> RUN.Build:
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_beads"])
    trapped = sorted(int(t) for t in P["trapped"])
    k_t = float(P["k_t_star"])
    amp, omega = float(P["amp_star"]), float(P["omega_star"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    mid = trapped[len(trapped) // 2]
    ghost_mid = n + trapped.index(mid)

    reduced = {"kappa_star": P["kappa_star"], "edl_amp": P["edl_amp"],
              "vdw_amp": P["vdw_amp"], "a_star": P["a_star"]}
    h_min_star = float(P["h_min_star"])
    ell_star = 1.0 + h_min_star                      # 결합 자연길이 (중심간, r*=1+h_min)
    r_cut_star = 1.0 + float(P["cutoff_h_star"])
    r_min_star = 1.0 + 1e-6
    box_star = 4.0 * (n - 1) * ell_star

    _pos_drive = str(P.get("drive_mode", "trap")) == "position"
    sim = SIM.make_sim(make_frame(n, ell_star, trapped, box_star,
                                  skip_trap=mid if _pos_drive else None), seed=seed)

    cell = md.nlist.Cell(buffer=0.2)
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_star, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_star, U=U_arr, F=F_arr)
    tab.params[("A", "G")] = dict(r_min=r_min_star, U=U_arr * 0, F=F_arr * 0)
    tab.params[("G", "G")] = dict(r_min=r_min_star, U=U_arr * 0, F=F_arr * 0)
    tab.r_cut[("A", "G")] = tab.r_cut[("G", "G")] = r_min_star     # 사실상 꺼짐

    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
    wca.params[("A", "G")] = wca.params[("G", "G")] = dict(epsilon=0.0, sigma=SIGMA_CORE_STAR)

    bond = md.bond.Harmonic()
    bond.params["trap"] = dict(k=k_t, r0=0.0)          # 트랩 = 유령과의 조화 본드. 백본 없음

    # ★ 대조군 — JKR 굽힘강성. `angle.Harmonic` 을 쓰지 **않는다**: 거의 곧은 사슬에서
    #   sinθ 클램프로 힘이 조용히 틀린다 (bd-hoomd 함정 15, chain-bend-2d-oscill 에서
    #   실측 규명). 그쪽과 **같은** 선형화 굽힘을 force.Custom 으로 직접 건다.
    forces = [tab, wca, bond]
    bend = None
    if bool(P.get("jkr", False)):
        kth_star = float(P["kappa_theta_star"])
        bend = make_bending_force(bending_matrix(n, kth_star, ell_star), n)
        forces.append(bend)

    # ★★ 구동 방식 두 가지 (유변학적으로 의미가 다르다)
    #   "trap"     — 트랩 중심을 y=a·sin(ωt) 로 옮기고, 비드는 트랩에 끌려온다.
    #                실험(광집게)에 가깝지만 **트랩 컴플라이언스** 때문에 비드가 명령
    #                위치를 못 따라간다 (실측 추종률 JKR 3.5% / DLVO 29%).
    #                → K* = k_t(ŷ_c/ŷ − 1)  [구동 트랩이 느끼는 강성]
    #   "position" — 구동 비드의 y 를 **직접** 강제한다 (변형률 제어). 컴플라이언스가
    #                아예 없어 변형이 정확히 부과되고, 유변학의 표준 프로토콜
    #                (변형을 주고 응력을 잰다)에 맞는다. [P1][P2] 도 실험적으로는
    #                중앙을 움직이고 **양끝 비드의 힘**을 센서로 잰다.
    #                → K*_transfer = k_t_sensor·ŷ_end / ŷ_mid  [전달 강성]
    #   ⚠️ 두 K* 는 **같은 양이 아니다** — 전자는 구동점 강성, 후자는 전달 함수다.
    #      값을 직접 비교하지 말고 각각 DLVO vs JKR 대비로 볼 것.
    pos_drive = str(P.get("drive_mode", "trap")) == "position"
    if pos_drive:
        # 구동 비드를 적분에서 제외 → 브라운 운동 없음, updater 가 위치를 직접 씀
        bd_filter = hoomd.filter.SetDifference(hoomd.filter.Type(["A"]),
                                               hoomd.filter.Tags([mid]))
    else:
        bd_filter = hoomd.filter.Type(["A"])
    bd = md.methods.Brownian(filter=bd_filter, kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    # position 모드는 **비드 자신**을, trap 모드는 **유령**을 움직인다
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=_move_ghost_action(mid if pos_drive else ghost_mid, amp, omega, dt),
        trigger=hoomd.trigger.Periodic(UPDATE_EVERY)))
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, int(Nm["n_prod"]) // 200))

    def pe_per_particle():
        return float(sum(np.array(f.energies).sum() for f in forces)) / n

    def sample(timestep, phase):
        snap = sim.state.get_snapshot()
        p = np.array(snap.particles.position, dtype=float)
        nb = p[:n, :2]
        sep = np.linalg.norm(nb[:, None, :] - nb[None, :, :], axis=-1)
        np.fill_diagonal(sep, np.inf)
        nn_sep = np.array([sep[i, i + 1] for i in range(n - 1)])   # NN 간격 (결합 감시)
        return {"t": timestep * dt, "y_bead": float(p[mid, 1]),
                "y_ghost": float(p[ghost_mid, 1]),
                "y_end0": float(p[trapped[0], 1]), "y_end1": float(p[trapped[-1], 1]),
                "min_sep_all": float(sep.min()), "nn_sep_max": float(nn_sep.max()),
                "nn_sep_min": float(nn_sep.min()),
                # ★★ 굽음 — 양끝을 잇는 직선 대비 최대 이탈. **강체 평행이동·회전을 뺀
                #   순수 굽힘 변형**이다. 이게 이 케이스의 핵심 판별량인 이유:
                #   y 전범위로 보면 DLVO/JKR 이 3.4배 차이지만, 굽음으로 보면 19.6배다 —
                #   JKR 사슬이 움직이는 것의 대부분은 (트랩이 유한강성이라) 사슬 통째의
                #   평행이동이지 굽힘이 아니기 때문. shape_localization 은 두 계를 전혀
                #   구별하지 못했다(1.481 vs 1.470) — 열잡음에 지배돼서다. 굽음이 그 대체다.
                "bow": float(_bow(p[:n, 1])),
                "shape_y": p[:n, 1].copy()}       # ★ 전체 y-프로필 — 매끈한 곡률 vs 국소꺾임 판정용

    def finalize(cols):
        t, yb, yg = cols["t"], cols["y_bead"], cols["y_ghost"]
        blocks_y = LI.lockin_blocks(t, yb, omega, n_blocks=min(10, max(2, len(t) // 20)))
        blocks_g = LI.lockin_blocks(t, yg, omega, n_blocks=min(10, max(2, len(t) // 20)))
        Kb = np.array([LI.k_star(a, b, k_t, omega) for a, b in zip(blocks_y, blocks_g)])
        K, Ksem = LI.agg(Kb)
        yh, _ = LI.agg(blocks_y)
        gh, _ = LI.agg(blocks_g)

        # ★★ position 구동이면 **전달 강성**을 따로 낸다.
        #   변형은 정확히 부과됐으므로(컴플라이언스 0) 응력에 해당하는 것은 양끝
        #   센서 비드가 받는 힘 F_end = k_t·y_end 다 ([P1][P2] 의 실험 프로토콜).
        #       K_transfer = k_t·⟨ŷ_end⟩ / ŷ_mid
        #   ⚠️ trap 모드의 K′(구동점 강성)와 **다른 양**이다. 값을 직접 비교하지 말 것.
        pos_mode = str(P.get("drive_mode", "trap")) == "position"
        K_tr = None
        if pos_mode:
            ye = 0.5 * (np.asarray(cols["y_end0"], float) + np.asarray(cols["y_end1"], float))
            be = LI.lockin_blocks(t, ye - ye.mean(), omega,
                                  n_blocks=min(10, max(2, len(t) // 20)))
            eh, esem = LI.agg(be)
            mh, _ = LI.agg(blocks_y)          # 구동 비드 = 부과된 변형 (측정 위상자)
            K_tr = (k_t * eh / mh) if abs(mh) > 0 else complex("nan")

        # ★ 예측은 브랜치에 따라 다르다 — 결과를 보기 전에 고정된다 (원칙 9.2).
        if bool(P.get("jkr", False)):
            # 대조군: 굽힘항이 있으면 구동 트랩이 느끼는 정적 강성이 예측값이다.
            #   (A_bend + T) y = k_t y_c e_mid  →  K = k_t(y_c/y_mid − 1)
            #   chain-bend-2d-oscill.driven_static_stiffness 와 같은 식.
            A_ = bending_matrix(n, float(P["kappa_theta_star"]), ell_star)
            for i in trapped:
                A_[i, i] += k_t
            e_ = np.zeros(n)
            e_[mid] = k_t
            K_pred = float(k_t * (1.0 / np.linalg.solve(A_, e_)[mid] - 1.0))
            K_derivation = ("굽힘항이 있으므로 선형응답 정적극한이 예측을 준다: "
                            "(A_bend+T)y = k_t·y_c·e_mid 를 풀어 K = k_t(y_c/y_mid − 1). "
                            "양끝이 강체가 아니라 유한강성 트랩이라 48EI/L³ 이 아니다. "
                            "저주파(준정적)에서만 성립 — De가 크면 점성이 섞인다.")
            K_note = (f"★ 대조군(JKR 굽힘 ON). 정적극한 예측 {K_pred:.4g} kT/d². "
                      f"DLVO-only 브랜치(예측 0)와 같은 기하·같은 시드에서 대조한다")
        else:
            K_pred = 0.0
            K_derivation = ("중심력 U(r), 결합이 자연장 ℓ(U'(ℓ)=0)에 있을 때 횡변위 "
                            "y_i,y_{i+1}에 대한 결합에너지는 U(ℓ)+U'(ℓ)(y_i-y_{i+1})²/(2ℓ)"
                            "+O(y⁴) — U'(ℓ)=0이라 O(y²) 항이 사라진다. 비드 수·본드 "
                            "토폴로지에 무관한 국소 대칭 논증이라 조합 전체(사슬)에도 "
                            "그대로 성립한다 (극한 조건: 작은 y, 트랩이 사슬을 늘이거나 "
                            "압축하지 않음).")
            K_note = ("★ 예측 0 — G1(직선사슬+순수중심력+자연장 평형이면 선형 굽힘강성"
                      " 정확히 0). 유의미하게 0이 아니면 장력/유한변형 효과가 섞인 것")

        shapes = np.stack(cols["shape_y"])          # (n_samples, n)
        # ★★ G1 검사: 매끈한 곡률이면 프로파일 이차미분이 완만하게 분포, "삼각형형"
        #   좌굴이면 특정 결합에서 각도(이차미분)가 국소적으로 튀어야 한다.
        d2 = np.diff(shapes, n=2, axis=1)            # (n_samples, n-2)
        d2_mean = np.abs(d2).mean(axis=0)
        kurt_like = float(d2_mean.max() / (d2_mean.mean() + 1e-30))   # 1에 가까움=매끈, 크면 국소집중

        # ★★ 굽음 — 이 케이스의 핵심 판별량 (_bow 도크스트링 참조).
        #   두 가지를 따로 낸다:
        #     bow_rms    시간평균 — **열요동 + 구동**이 섞여 있다
        #     bow_drive  구동 주파수 성분만 락인으로 뽑은 것 — **구동에 의한 굽힘만**
        #   후자가 진짜 응답이다. 전자만 보면 열적으로 흐물대는 사슬(DLVO)이 "많이 휜다"고
        #   나오는데 그건 구동에 대한 탄성 응답이 아니다.
        bow = np.asarray(cols["bow"], dtype=float)
        bow_rms = float(np.sqrt((bow ** 2).mean()))
        bb = LI.lockin_blocks(t, bow - bow.mean(), omega,
                              n_blocks=min(10, max(2, len(t) // 20)))
        bow_hat, bow_sem = LI.agg(bb)
        bow_drive = float(abs(bow_hat))

        # ★★ position 구동에서는 K′(구동 트랩 강성)이 **정의되지 않는다** — 트랩이 없다.
        #   그대로 계산하면 유령이 안 움직여 ŷ_c=0 이 되고 K′ = k_t(0/ŷ − 1) = −k_t 라는
        #   **의미 없는 상수**가 나온다 (실측 −5217.11 = −k_t*, K_sem=1.2e-12 로 산포조차
        #   0). 숫자가 그럴듯해서 실제 측정값으로 오독되기 딱 좋다 — 아예 내보내지 않는다.
        #   position 모드의 관측량은 K_transfer_* 다.
        obs = ([] if pos_mode else [
            MET.observable("K_prime", float(K.real), K_pred, "kT/d^2",
                           role="implementation_check", sigma=Ksem, tol_sigma=3.0,
                           derivation=K_derivation,
                           note=K_note),
            MET.observable("K_doubleprime", float(K.imag), None, "kT/d^2", role="measurement"),
            MET.observable("K_sem", Ksem, None, "kT/d^2", role="measurement"),
        ]) + [
            MET.observable("y_response", float(abs(yh)), None, "d", role="measurement"),
            *([MET.observable("K_transfer_prime", float(K_tr.real), None, "kT/d^2",
                              role="measurement",
                              note="★ 위치(변형률) 제어의 전달 강성 실수부 = k_t·ŷ_end/ŷ_mid. "
                                   "변형이 정확히 부과되므로 트랩 컴플라이언스가 없다 "
                                   "([P1][P2] 실험 프로토콜과 같은 구성: 중앙을 움직이고 "
                                   "양끝 센서의 힘을 잰다). ⚠ trap 모드 K′ 과 다른 양"),
               MET.observable("K_transfer_dprime", float(K_tr.imag), None, "kT/d^2",
                              role="measurement", note="전달 강성 허수부")]
              if K_tr is not None else []),
            MET.observable("bow_rms", bow_rms, None, "d", role="measurement",
                           note="★★ 굽음 RMS — 양끝 잇는 직선 대비 최대이탈. 강체 평행이동을 "
                                "뺀 순수 굽힘 변형이다 (트랩이 유한강성이라 사슬이 통째로 "
                                "밀리는데, y 변위로 재면 그게 굽힘 신호를 덮는다). "
                                "실측 판별력: y전범위 3.4배 vs 굽음 19.6배 (DLVO vs JKR). "
                                "★ 단 이 값은 열요동+구동이 섞여 있다 — 구동 응답만 보려면 "
                                "아래 bow_drive"),
            MET.observable("bow_drive", bow_drive, None, "d", role="measurement",
                           sigma=bow_sem,
                           note="★ 구동 주파수 성분만 락인으로 뽑은 굽음 = **구동에 의한 "
                                "굽힘 응답**. bow_rms 만 보면 열적으로 흐물대는 사슬이 '많이 "
                                "휜다'고 나오지만 그건 탄성 응답이 아니다"),
            MET.observable("shape_localization", kurt_like, None, "dimensionless",
                           role="measurement",
                           note="⚠ **판별력 없음이 실측으로 확인됨** (JKR 1.481 vs DLVO 1.470 — "
                                "육안으로는 명백히 다른데 지표는 동일). 양쪽 모두 열잡음이 "
                                "이차미분을 지배해서다. 남겨두되 판정에 쓰지 말 것 — "
                                "대체는 bow_drive. max|θ''|/mean|θ''| 정의"),
        ]
        return {"observables": obs,
                "extra": {**({} if pos_mode else
                              {"K_prime": float(K.real), "K_doubleprime": float(K.imag),
                               "K_sem": Ksem}),
                          "y_resp_abs": float(abs(yh)),
                          "drive_abs": float(abs(gh)), "omega_star": omega,
                          **({"K_transfer_prime": float(K_tr.real),
                              "K_transfer_dprime": float(K_tr.imag)} if K_tr is not None else {}),
                          "drive_mode": str(P.get("drive_mode", "trap")),
                          "bow_rms": bow_rms, "bow_drive": bow_drive,
                          "bow_drive_sem": float(bow_sem),
                          "bow_max": float(bow.max()),
                          "shape_localization": kurt_like,
                          "shape_profile_mean": shapes.mean(axis=0).tolist(),
                          "min_sep_all": float(cols["min_sep_all"].min()),
                          "nn_sep_max": float(cols["nn_sep_max"].max()),
                          "nn_sep_min": float(cols["nn_sep_min"].min())}}

    every = max(1, int(Nm["sample_every"]))
    return RUN.Build(
        sim=sim, forces=forces, n_particles=n,
        sample=sample, pe_per_particle=pe_per_particle, sample_every=every,
        phases=[RUN.Phase("예열", int(Nm["n_eq"]), collect=False,
                          note="구동 ON · 국소(결합) 이완만 — 전체 형태이완 아님(스모크)"),
                RUN.Phase("생산", int(Nm["n_prod"]), every,
                          note=f"ω*={omega:.4g}")],
        tags=["2D", "chain", "dlvo", "pair_table", "oscillatory_drive", "no_bending",
              "hypothesis_test", "structural"],
        physical={"n_beads": n, "omega_star": omega, "amp_star": amp, "k_t_star": k_t,
                  "well_depth_star": float(P["well_depth_star"]),
                  "k_bond_star": float(P["k_bond_star"])},
        finalize=finalize)


if __name__ == "__main__":
    sys.exit(main())
