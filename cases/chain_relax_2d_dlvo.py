"""`chain-relax-2d-dlvo` — chain-bend-2d-dlvo 의 구동 없는(자유 이완) 짝.

같은 물리(PMMA 비드, DLVO 2차극소 중심력 결합, 명시적 굽힘/마찰 없음)를 쓰되,
**트랩·오실레이션을 전부 뗀다** — 사용자가 요청한 "입자간 마찰력은 없고 인력만
있는 체인"을 구동 없이 가장 단순한 형태로 본다 (CLAUDE.md 규칙 8: 정적인 계를
먼저, 움직임은 그 다음). chain-bend-2d-dlvo 는 이 정적 단계를 건너뛰고 곧바로
오실레이션 실험을 했다 — 이 케이스가 그 빠진 단계를 채운다.

두 실험 (--init):
  straight  직선 사슬을 열평형시켜 **결합 신장(radial)** 열요동을 잰다.
            예측: ⟨δℓ*²⟩ = kT*/k_bond* (트랩의 ⟨x²⟩=kT/k 와 같은 종류의 등분배
            골든테스트 — DLVO 표 퍼텐셜을 이 케이스로 이식한 게 맞는지 확인한다).
            implementation_check.
  kink      중앙 결합 하나에 정확한 턴각 Δφ 를 주고(다른 모든 결합은 완전히 곧고
            모든 결합이 정확히 자연장에서 시작 — 신장 신호 없이 순수 굽힘만) 풀어서
            굽음(bow)이 탄성 복원되는지 확산적으로 흩어지는지 본다. G1(횡방향 선형
            굽힘강성이 정확히 0)은 대수적으로만 유도돼 있었고 구동 없이 직접 확인한
            적이 없다 — 사전 정량 예측이 없으므로 measurement.

DLVO 식·표 퍼텐셜 구현은 **재정의하지 않고** chain_bend_dlvo_2d 에서 그대로 임포트한다
(두 번째 케이스 — network 가 이미 같은 방식으로 재사용한 전례를 따름. bdbot/ 로
올릴지는 세 번째 케이스가 나올 때 판단, CLAUDE.md "두 번 나왔는가" 원칙).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/chain_relax_2d_dlvo.py --init straight --report
    $PY cases/chain_relax_2d_dlvo.py --init kink --kink-angle 0.3 --report
    $PY cases/chain_relax_2d_dlvo.py --init straight --smoke --run     # 빠른 정합성 확인
    $PY cases/chain_relax_2d_dlvo.py --init kink --run
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
sys.path.insert(0, str(Path(__file__).resolve().parent))

from bdbot import checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, scales as SC, sim as SIM, stats as ST  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from bdbot.units import Q  # noqa: E402

# ★ DLVO 식은 chain-bend-2d-dlvo 에서 SI로 검증된 것을 **그대로** 쓴다 (두 번 적지 않는다,
#   network 와 같은 전례).
from chain_bend_dlvo_2d import (  # noqa: E402
    CUTOFF_H_STAR, SIGMA_CORE_STAR, build_table_arrays, dlvo_reduced_params, find_well,
    F_h_star, U_star,
)

ROOT = Path(__file__).resolve().parent.parent


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
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② 기하 — 회전/평행이동 불변 형태 서술자
#    ★★ chain-bend-2d-dlvo 의 `_bow(y)` 는 트랩이 방향을 고정해줘서 실험실좌표 y를
#    그대로 썼다. 이 케이스엔 트랩이 없어 사슬 전체가 자유로이 회전·평행이동한다 —
#    반드시 사슬 자신의 몸좌표계(양끝을 잇는 축)로 회전시켜야 한다. 안 그러면
#    "굽음"이 실제로는 그냥 열적 회전 표류를 재는 것이 된다.
# ════════════════════════════════════════════════════════════════════════
def bond_vectors(pos: np.ndarray) -> np.ndarray:
    return pos[1:] - pos[:-1]


def bend_angles(pos: np.ndarray) -> np.ndarray:
    """내부 비드마다 결합방향의 국소 턴각 dtheta_i (n-2,). 0=그 지점이 국소적으로 곧다.

    회전/평행이동 불변 — 연속 결합방향의 **차이**만 본다. chain_bend_dlvo_2d.py 의
    `bending_matrix()` 가 쓰는 이산곡률 θ_i=(y_{i+1}-2y_i+y_{i-1})/ell 과 소각도에서
    동등하다(둘 다 1차 이산곡률) — 여기선 트랩이 없어 y 기반 정의를 못 쓰므로
    턴각으로 재정의했다.
    """
    bv = bond_vectors(pos)
    ang = np.arctan2(bv[:, 1], bv[:, 0])
    dtheta = np.diff(ang)
    return (dtheta + np.pi) % (2 * np.pi) - np.pi


def bow_metrics(pos: np.ndarray) -> tuple[float, float]:
    """양끝을 잇는 축을 x'축으로 정렬한 뒤 그 축에서 벗어난 정도 (max, rms) [d]."""
    d = pos[-1] - pos[0]
    L_ee = float(np.hypot(d[0], d[1]))
    if L_ee < 1e-9:
        return 0.0, 0.0
    u = d / L_ee
    rel = pos - pos[0]
    yp = -rel[:, 0] * u[1] + rel[:, 1] * u[0]        # 몸좌표계 횡성분. yp[0]=yp[-1]=0(구성상)
    return float(np.abs(yp).max()), float(np.sqrt(np.mean(yp ** 2)))


def min_nnn_gap_star(pos: np.ndarray) -> float:
    """|i-j|>=2 인 쌍의 최소 표면간극 h*=r*-1. 조기 비인접 결합("삼각형형") 감지."""
    n = len(pos)
    if n < 4:
        return float("inf")
    diff = pos[:, None, :] - pos[None, :, :]
    dist = np.hypot(diff[..., 0], diff[..., 1])
    idx = np.abs(np.subtract.outer(np.arange(n), np.arange(n)))
    mask = idx >= 2
    return float(dist[mask].min()) - 1.0


def kink_positions(n: int, ell_star: float, kink_angle: float) -> np.ndarray:
    """중앙 결합 하나에만 정확한 턴각 `kink_angle` 을 주고 나머지는 완전히 곧다.

    모든 결합이 정확히 `ell_star` 로 시작한다 — 신장(radial) 신호를 섞지 않고
    순수 굽힘 자유도만 교란한다. `kink_angle=0` 이면 직선(= --init straight).
    """
    mid = n // 2
    pos = np.zeros((n, 2))
    th_l, th_r = -kink_angle / 2.0, kink_angle / 2.0
    for k in range(1, n - mid):
        pos[mid + k] = pos[mid + k - 1] + ell_star * np.array([math.cos(th_r), math.sin(th_r)])
    for k in range(1, mid + 1):
        pos[mid - k] = pos[mid - k + 1] - ell_star * np.array([math.cos(th_l), math.sin(th_l)])
    pos -= pos.mean(axis=0)
    return pos


def bond_variance_boltzmann(p: dict, w: dict, cutoff_h_star: float, nbins: int = 400_000):
    """결합 신장(h-h_min)의 **진짜**(비조화) 열평형 분산 — 2차극소 우물 안(h>barrier_h)
    에서 수치 적분 (basin-restricted Boltzmann average, kT*=1).

    ★★ 조화 근사 `1/k_bond_star` 는 우물 바닥의 국소 곡률만 본다. 이 우물은 안쪽(장벽
    쪽)이 가파르고 바깥쪽(h→∞, U→0⁻)이 훨씬 무른 **비대칭**이라 조화 근사가 진값을
    과소평가한다 — 실측(스모크런, n=9): 조화 예측의 2.9배, 전 우물 수치적분은 4.6배.
    soft-r3 의 "Einstein 케이지 근사가 비조화성으로 어긋난다"(bd-physics §6.2)와
    같은 종류의 함정이다. 여기서는 **체제 판정용으로 낮추지 않고** 적분으로 정확한
    예측값을 만들어 진짜 골든테스트로 쓴다 — 이 계는 조화 근사를 요구하지 않는다
    (자유 이완, 트랩 없음. soft-r3 는 dt 설계에 근사가 필요해서 못 피했다).

    적분 구간을 장벽(barrier_h) 안쪽까지 열지 않는 이유: 장벽이 416 kT 라 시뮬레이션
    시간축에서 **절대 못 넘는다** — 전 구간(h→0의 1차극소 포함)으로 적분하면 vdW
    발산 때문에 분포가 h→0 으로 완전히 무너진다(실측 확인). 시뮬레이션이 실제로
    샘플링하는 것은 **2차극소 우물 안**뿐이므로 그 basin 으로 제한해야 한다.
    """
    h = np.linspace(w["barrier_h"], cutoff_h_star, nbins)
    U = U_star(h, p)
    wt = np.exp(-(U - U.min()))
    Z = np.trapezoid(wt, h)
    mean_h = float(np.trapezoid(h * wt, h) / Z)
    var_h = float(np.trapezoid((h - mean_h) ** 2 * wt, h) / Z)
    return var_h, mean_h


# ════════════════════════════════════════════════════════════════════════
# ③ 스케일 원장
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, n: int, kink_angle: float, *, dt_scale=1.0,
                 T_obs_tau: float, eq_scale: float = 200.0) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B, tau_p = b["kT"], b["gamma"], b["tau_B"], b["tau_p"]

    p = dlvo_reduced_params(sys_)
    w = find_well(p)
    ell = d * (1 + w["h_min"])
    k_bond = Q(w["k_bond_star"], "dimensionless") * kT / d ** 2
    sigma_bond = (kT / k_bond) ** 0.5
    L_chain = (n - 1) * ell

    # ★★ 이 계의 유일한 강성 모드가 결합 신장이다 (굽힘강성은 구조적으로 0, 트랩도
    #   없다) — chain-bend-2d-dlvo 처럼 여러 후보 중 최속을 고를 필요가 없다.
    tau_bond = C.relaxation_time(gamma, k_bond)
    dt = dt_scale * C.dt_from_gate(tau_bond)
    T_obs = Q(T_obs_tau, "dimensionless") * tau_bond

    # 킹크 초기조건의 안전성 — 결합은 전부 정확히 ell 에서 시작하므로 "신장"은 0이지만
    # 킹크 자체가 만드는 NNN(2결합 건너) 간극은 좁아질 수 있다 (기하로 직접 확인).
    ell_star = float((ell / d).to("dimensionless").magnitude)
    nnn_gap0 = min_nnn_gap_star(kink_positions(n, ell_star, kink_angle)) if kink_angle else float("inf")

    # 결합이 h_min 밖으로 얼마나 늘어나야 최대 인장력(F_max)에 닿는지 — 정보용 앵커.
    # (킹크는 신장 신호가 없게 지었으므로 여기서 실제로 쓰이진 않지만, dt/안전 논의의
    # 참조 스케일로 원장에 남긴다.)
    hs = np.linspace(w["h_min"], CUTOFF_H_STAR, 20_000)
    Fs = F_h_star(hs, p)
    F_max = float(-Fs.min())

    # ★★ 골든테스트의 진짜 예측값 — 조화 근사(1/k_bond_star)가 아니라 우물 전체의
    #   비조화 볼츠만 적분(basin-restricted, 위 bond_variance_boltzmann 참조).
    dl_var_boltz, dl_mean_boltz_h = bond_variance_boltzmann(p, w, CUTOFF_H_STAR)

    lg = SC.ScaleLedger()
    lg.add_length("sigma_bond", sigma_bond.to("m"), "결합 신장 열요동 폭 √(kT/k_bond)", star=True)
    lg.add_length("h_min", Q(w["h_min"], "dimensionless") * d, "2차극소 위치(표면간극)")
    lg.add_length("d", d, "비드 지름")
    lg.add_length("ell", ell.to("m"), "결합 자연길이 (중심간, d+h_min)")
    lg.add_length("L_chain", L_chain.to("m"), "사슬 윤곽길이 (n-1)ell")
    lg.add_time("tau_p", b["tau_p"], "m/γ 관성 이완", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_bond", tau_bond,
               "★★ γ/k_bond 결합 신장 — 이 계의 **유일한** 강성 모드. dt를 정한다",
               star=True)
    lg.add_time("tau_B", tau_B, "d²/D_t 확산 (기준)")
    lg.add_time("T_obs", T_obs, "관측창 (τ_bond 배수 — 국소 모드만 필요, 사슬 전체"
               " 형태이완 τ_chain_diff 는 불필요, observation.yaml R2)", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.add_energy("k_bond_d2", (k_bond * d ** 2).to("J"), "k_bond d² 결합 신장강성", star=True)
    lg.add_energy("well_depth", Q(-w["U_min"], "dimensionless") * kT, "|2차극소 깊이|")
    lg.declare_absent(
        "box",
        "주기경계 없음 (사슬 하나, 트랩도 없어 자유 이완·자유 회전). HOOMD 프레임에 "
        "형식적 박스가 필요하므로 사슬 자신의 최대 크기(4×L_chain, chain-bend-2d-dlvo와 "
        "같은 여유)로 크게 잡아 자기 주기이미지와의 상호작용만 피한다 — 물리적으로 "
        "의미 있는 확인 스케일이 아니다.")
    lg.declare_absent(
        "bending_stiffness",
        "★★ 구조적으로 없다 (G1, chain-bend-2d-dlvo 에서 유도) — 이 케이스가 그 "
        "구조적 사실을 구동 없는 최소 구성에서 직접 실행으로 확인하는 것 자체가 "
        "목적이다. 지어낸 대체 척도를 넣지 않는다.")
    lg.derived = dict(gamma=gamma, kT=kT, d=d, tau_B=tau_B, ell=ell.to("m"),
                      L_chain=L_chain.to("m"), ell_star=ell_star, k_bond=k_bond,
                      k_bond_star=w["k_bond_star"], sigma_bond=sigma_bond.to("m"),
                      tau_bond=tau_bond, dt=dt, T_obs=T_obs, reduced=p,
                      h_min_star=w["h_min"], well_star=w["U_min"], barrier_star=w["barrier_U"],
                      F_max_star=F_max, nnn_gap0_star=nnn_gap0, kink_angle=kink_angle, n=n,
                      dl_var_boltz=dl_var_boltz, dl_mean_boltz_h=dl_mean_boltz_h)
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ 이 계엔 굽힘강성도 트랩도 없다 — 유일한 강성 모드가 "
        "결합 신장(tau_bond)뿐이라, dt·평형화 판정이 chain-bend-2d-dlvo 보다 단순하다.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ④ 무차원수 + 분리 검사
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(lg, n, init, kink_angle):
    D = lg.derived
    r = lg.ratio

    groups = [
        ND.Group("k_bond_star", D["k_bond_star"], ("energies", "k_bond_d2"),
                 ("energies", "kT"), "k_bond d²/kT", "결합 신장강성 — 유일한 조화 모드"),
        ND.Group("well_depth/kT", r("energies", "well_depth", "kT"), ("energies", "well_depth"),
                 ("energies", "kT"), "", "결합 깊이 — 열에너지 규모면 가역적"),
        ND.Group("sigma_bond/d", r("lengths", "sigma_bond", "d"), ("lengths", "sigma_bond"),
                 ("lengths", "d"), "", "결합 신장 열요동 (예측 골든테스트)"),
        ND.Group("n_beads", float(n), None, None, "", "사슬 길이 (입력)"),
        ND.Group("kink_angle_rad", float(kink_angle), None, None, "",
                 "초기 킹크각 (입력, straight면 0)"),
        ND.Group("dt/tau_bond", r("times", "dt", "tau_bond"), ("times", "dt"),
                 ("times", "tau_bond"), "", "적분 해상 — 유일한 강성 모드"),
        ND.Group("T_obs/tau_bond", r("times", "T_obs", "tau_bond"), ("times", "T_obs"),
                 ("times", "tau_bond"), "", "관측창 (국소모드 통계용)"),
        ND.Group("St", r("times", "tau_p", "tau_bond"), ("times", "tau_p"),
                 ("times", "tau_bond"), "tau_p/tau_bond", "관성 vs 결합신장"),
    ]
    checks = [
        C.Check("모델", "참고: τ_p/τ_bond", r("times", "tau_p", "tau_bond"), C.GATE, "<=",
              "chain-bend-2d-dlvo 와 동일한 미검증 상태(같은 입자·같은 결합) — 여기서 "
              "재검증하지 않는다. 그 케이스가 OverdampedViscous 대조로 검증하면 이 "
              "케이스도 같이 검증된다", hard=False),
        C.Check("적분", "결합 신장 해상 dt/τ_bond", r("times", "dt", "tau_bond"), C.GATE, "<=",
              "이 계의 유일한 강성 모드. 못 맞추면 발산"),
        C.Check("통계", "관측창 충분     T_obs/τ_bond", r("times", "T_obs", "tau_bond"),
              1000.0, ">=", "국소(결합) 등분배 통계에 필요한 최소 배수 — 사슬 전체 "
              "형태이완(τ_chain_diff)까지는 불필요", hard=False),
    ]
    if kink_angle:
        checks.append(C.Check(
            "기하", "킹크 NNN 간극(초기) vs 컷오프", D["nnn_gap0_star"], CUTOFF_H_STAR, ">=",
            f"2결합 건너 비드가 방출 순간부터 이미 DLVO 컷오프 안(조기 비인접 결합)에 "
            f"들어가 있으면 '순수 굽힘만 교란'했다는 설계 의도가 깨진다", hard=False))
    return groups, checks


def report_blocks(sys_, lg, n, init, kink_angle, n_steps):
    D = lg.derived
    inp = [R.kv("d", f"{sys_['d'].value:~.4gP}", sys_["d"].tier, sys_["d"].source[:44]),
           R.kv("psi0", f"{sys_['psi0'].value:~.4gP}", sys_["psi0"].tier, sys_["psi0"].source[:44]),
           R.kv("I(MgCl2)", f"{sys_['ionic_strength'].value:~.4gP}",
                sys_["ionic_strength"].tier, sys_["ionic_strength"].source[:44]),
           R.kv("n", f"{n}", 3, "chain-bend-2d-dlvo 승계값과 맞춤"),
           R.kv("init", init, 3, "straight=직선 열평형 / kink=방출 이완"),
           R.kv("kink_angle", f"{kink_angle:.3f} rad", 3, "★제안 (observation.yaml R1)")]
    der = [
        f"  결합: 2차극소 {D['well_star']:.3f} kT @ h_min*={D['h_min_star']:.5f}"
        f"   장벽 {D['barrier_star']:.1f} kT",
        f"  k_bond = {D['k_bond'].to('pN/um'):~.4fP} = {D['k_bond_star']:.4e} kT/d²"
        f"   σ_bond(조화) = {D['sigma_bond'].to('nm'):~.4fP}",
        f"  ★ 결합 신장 골든테스트 예측(비조화 볼츠만 적분) = {D['dl_var_boltz']:.4e} d²"
        f"   ({D['dl_var_boltz']*D['k_bond_star']:.2f}× 조화근사 — 우물이 바깥쪽으로 "
        f"무른 비대칭 탓, soft-r3 의 Einstein 케이지와 같은 종류의 비조화 보정)",
        f"  ell(자연장) = {D['ell'].to('nm'):~.2fP}   L_chain = {D['L_chain'].to('um'):~.3fP}",
        f"  F_max(우물 최대 인장력) = {D['F_max_star']:.1f} kT/d — 킹크는 결합신장 "
        f"신호를 안 만들어서(정확히 ell 에서 시작) 참조용",
        f"  ★★ 굽힘 선형강성 = 0 (구조적, declare_absent) — 이 실행이 직접 확인 대상",
    ]
    if kink_angle:
        der.append(f"  킹크 NNN 간극(초기, h*) = {D['nnn_gap0_star']:+.4f}"
                   f"   (컷오프 {CUTOFF_H_STAR:.2f} — 크면 조기 비인접결합 없음)")
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}  = {lg.ratio('times','dt','tau_B'):.3e} τ_B",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}  = {lg.ratio('times','T_obs','tau_bond'):.3g} τ_bond",
        f"  steps   = {n_steps:,}   × n={n}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# main — L3
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, n, init, kink_angle, args):
    T_obs_tau = args.tobs if args.tobs is not None else (3.0e3 if init == "kink" else 2.0e4)
    lg = build_ledger(sys_, n, kink_angle, dt_scale=args.dt_scale, T_obs_tau=T_obs_tau,
                      eq_scale=args.eq_scale)
    D = lg.derived
    dt = lg.get("times", "dt")

    if init == "kink":
        n_eq = 0
        n_release = int(round(float((D["T_obs"] / dt).to(""))))
        sample_every = max(1, n_release // args.samples)
        n_release = (n_release // sample_every) * sample_every
        n_prod = n_release
    else:
        n_eq = int(round(args.eq_scale * float((D["tau_bond"] / dt).to(""))))
        n_prod = int(round(float((D["T_obs"] / dt).to(""))))
        sample_every = max(1, n_prod // args.samples)
        n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(lg, n, init, kink_angle)
    tag = f"n{n}-{init}"
    if init == "kink":
        tag += f"-a{kink_angle:.3f}"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"

    p = D["reduced"]
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_beads": n, "init": init, "kink_angle": kink_angle,
                "kappa_star": p["kappa_star"], "edl_amp": p["edl_amp"], "vdw_amp": p["vdw_amp"],
                "a_star": p["a_star"], "cutoff_h_star": CUTOFF_H_STAR,
                "h_min_star": D["h_min_star"], "k_bond_star": D["k_bond_star"],
                "well_star": D["well_star"],
                # ★ L3에서 한 번 적분한 진짜(비조화) 골든테스트 예측값 — L4는 이 숫자를
                #   그대로 읽는다(스펙이 유일한 계약, 재유도하지 않는다).
                "dl_var_boltz": D["dl_var_boltz"]},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": args.seed},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_eq, n_prod


def emit(sys_, n, init, kink_angle, args) -> int:
    lg, spec, groups, checks, n_eq, n_prod = build_spec(sys_, n, init, kink_angle, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 무결성 검사")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n, init, kink_angle, n_eq + n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  n={n} init={init}"
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

    # ★★ result.txt — 완료 마커 (CLAUDE.md: 안 쓰면 status가 런을 0개로 세고, 미완료
    #   정리 스크립트가 완료 런을 지운다). 케이스 스크립트의 책임.
    if v["status"] != "skipped":
        obs_lines = []
        try:
            mj = json.loads((outdir / "metrics.json").read_text())
            for o in mj.get("observables", []):
                m = o.get("measured")
                p_ = o.get("predicted")
                tail = f"   (예측 {p_:.6g})" if isinstance(p_, (int, float)) else ""
                obs_lines.append(f"  {o['name']:<28} {m:.6g}{tail}" if m is not None
                                 else f"  {o['name']:<28} —")
        except Exception as e:
            obs_lines.append(f"  (metrics.json 을 읽지 못함: {e})")
        result = "\n".join(["=" * 84, f"결과 — {run_id}", "=" * 84,
                            *obs_lines, "=" * 84, verdict_txt])
        (outdir / "result.txt").write_text(report + "\n" + result)
        make_plots(sys_, lg, n, init, kink_angle, outdir)
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


# ════════════════════════════════════════════════════════════════════════
# L4 — HOOMD 빌더. 트랩·구동이 없어 chain-bend-2d-dlvo 보다 훨씬 단순하다 —
# 유령입자·CustomUpdater·force.Custom 전부 불필요. 실제 비드 N개 + Table(DLVO) +
# WCA(코어) + Brownian(전 입자) 뿐이다.
# ════════════════════════════════════════════════════════════════════════
@RUN.builder("chain-relax-2d-dlvo")
def build(spec, outdir=None) -> RUN.Build:
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_beads"])
    init, kink_angle = str(P["init"]), float(P["kink_angle"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])

    reduced = {"kappa_star": P["kappa_star"], "edl_amp": P["edl_amp"],
              "vdw_amp": P["vdw_amp"], "a_star": P["a_star"]}
    h_min_star = float(P["h_min_star"])
    ell_star = 1.0 + h_min_star
    r_cut_star = 1.0 + float(P["cutoff_h_star"])
    r_min_star = 1.0 + 1e-6
    box_star = 4.0 * max(1, n - 1) * ell_star           # 여유 큰 박스 — 자기 주기이미지 회피용

    pos0 = kink_positions(n, ell_star, kink_angle)
    sim = SIM.make_sim(SIM.frame_2d(pos0, box_star), seed=seed)

    cell = md.nlist.Cell(buffer=0.2)
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min_star, r_cut_star)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut_star)
    tab.params[("A", "A")] = dict(r_min=r_min_star, U=U_arr, F=F_arr)

    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * 2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)

    integ, bd = SIM.attach_brownian(sim, dt, [tab, wca])           # 마찰항 없음 — BD 그대로
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, n_prod // 200))
    L = box_star

    def unwrapped_xy():
        # ★ `sim.state.get_snapshot()` 은 (force.Custom 의 cpu_local_snapshot 과 달리)
        #   이미 tag 순서로 모아서 준다 — tag 재색인이 필요 없다(soft_r3_2d.py 의
        #   xy() 와 같은 패턴). tag 인덱싱은 로컬 스냅샷 전용(bd-hoomd 함정 설명).
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, dtype=float)[:, :2]
        img = np.array(snap.particles.image, dtype=float)[:, :2]
        return pos + img * L

    def pe_pp():
        return float(np.array(tab.energies).sum() + np.array(wca.energies).sum()) / n

    def sample(timestep, phase):
        xy = unwrapped_xy()
        bow_max, bow_rms = bow_metrics(xy)
        dtheta = bend_angles(xy)
        bl = np.hypot(*bond_vectors(xy).T)
        return {"bow_max": bow_max, "bow_rms": bow_rms,
                "dtheta_var": float(np.mean(dtheta ** 2)) if len(dtheta) else 0.0,
                "dl_var": float(np.mean((bl - ell_star) ** 2)),
                "min_sep": float(bl.min()) if len(bl) else np.nan,
                "nnn_gap": min_nnn_gap_star(xy),
                "dtheta": dtheta, "bond_len": bl}

    def finalize(cols):
        n_s = len(cols["pe"])
        dtheta_all = np.concatenate(cols["dtheta"]) if n_s else np.zeros(0)
        dl_all = np.concatenate([bl - ell_star for bl in cols["bond_len"]]) if n_s else np.zeros(0)
        dtheta_var = float(np.mean(dtheta_all ** 2)) if dtheta_all.size else float("nan")
        dl_var = float(np.mean(dl_all ** 2)) if dl_all.size else float("nan")
        dl_var_sem = ST.block_sem(np.array([float(np.mean((bl - ell_star) ** 2))
                                            for bl in cols["bond_len"]])) if n_s else float("nan")
        k_bond_star = float(P["k_bond_star"])
        dl_var_pred = float(P["dl_var_boltz"])         # ★ 비조화 볼츠만 적분값 — 아래 참조

        obs = []
        # ── implementation_check — 결합 신장(radial) 등분배 골든테스트 ──────────
        # ★★ 예측은 조화 근사(kT*/k_bond*)가 **아니다** — 스모크런에서 실측해보니
        #   조화 근사가 진값을 4.6배 과소평가했다(우물이 바깥쪽으로 무른 비대칭 —
        #   soft-r3 의 Einstein 케이지 근사와 같은 함정). 대신 build_ledger() 의
        #   `bond_variance_boltzmann()` 이 우물 전체(장벽 안쪽 basin)를 수치적분한
        #   값을 쓴다 — 여전히 **이 케이스가 세우는 같은 모델(DLVO+WCA)에서 그대로
        #   나오는 예측**이라 implementation_check 이 맞다(조화라고 가정하지 않았을
        #   뿐, 이 퍼텐셜 자체의 정확한 결과). DLVO 표 퍼텐셜을 이 케이스로 이식한 게
        #   맞는지 확인하는 역할(chain-bend-2d-dlvo 는 트랩+구동이 섞여 이 자유도를
        #   격리해서 잰 적이 없다).
        obs.append(MET.observable(
            "결합 신장 등분배 ⟨δℓ*²⟩", dl_var, dl_var_pred, "d²",
            "boltzmann_integral", role="implementation_check", scope="module", tol_pct=8.0,
            sigma=dl_var_sem if dl_var_sem > 0 else None,
            note=f"예측 = 우물(basin) 비조화 볼츠만 적분(정확). 참고: 조화근사 "
                 f"kT*/k_bond*={1.0 / k_bond_star:.3e} 는 이 값의 "
                 f"1/{dl_var_pred * k_bond_star:.2f} 뿐 — 조화근사를 예측으로 쓰면 "
                 f"안 된다(체제 판정용으로만)"))
        # ── measurement — 굽힘(각) 자유도. G1 은 '선형강성=0'만 말하고 크기는 "
        #    예측하지 않는다(고차항 지배) — 사전 예측이 없어 measurement.
        obs.append(MET.observable(
            "굽힘각 열요동 ⟨dθ²⟩", dtheta_var, None, "rad²", "none",
            role="measurement",
            note="G1(선형강성=0)의 직접 결과 — 요동 폭 자체는 예측이 없다(2차 이상 "
                 "항·배제부피가 지배). κ_θ,eff=1/⟨dθ²⟩ 로 등분배를 가정하면 "
                 f"{(1.0 / dtheta_var if dtheta_var > 0 else float('nan')):.3g} kT/rad² — "
                 "참고용 숫자일 뿐 조화 모드라는 뜻은 아니다"))

        min_sep = float(np.min(cols["min_sep"]))
        nnn_min = float(np.min(cols["nnn_gap"]))
        post_checks = [
            C.Check("기하", "표 하한 여유 r_table_min/min_sep",
                  (1.0 + 1e-6) / min_sep, 1.0, "<=",
                  f"pair.Table 함정 11: r<r_min 이면 힘이 0. 측정 최소 결합길이 {min_sep:.4f}"),
            C.Check("기하", "NNN 간극 최소(런 전체) vs 컷오프", nnn_min, CUTOFF_H_STAR, ">=",
                  "런 도중 2결합 건너 비드가 DLVO 컷오프 안으로 들어온 적이 있는지 "
                  "(들어오면 '삼각형형' 국소 접힘 후보 — 실패 아니라 관찰 대상)",
                  hard=False),
        ]

        bow_final = float(np.mean(cols["bow_rms"][-max(1, n_s // 10):])) if n_s else float("nan")
        bow_final_sem = ST.block_sem(cols["bow_rms"][-max(1, n_s // 10):]) if n_s >= 8 else float("nan")
        post_dicts = [{**c.as_dict("post_run"), "note": c.note} for c in post_checks]
        return {"observables": obs,
                "extra": {"dtheta_var": dtheta_var, "dl_var": dl_var, "dl_var_sem": dl_var_sem,
                          "dl_var_pred": dl_var_pred, "k_bond_star": k_bond_star,
                          "ell_star": ell_star, "min_sep": min_sep, "nnn_min": nnn_min,
                          "bow_initial": float(cols["bow_rms"][0]) if n_s else float("nan"),
                          "bow_final_rms": bow_final, "bow_final_sem": bow_final_sem,
                          "post_checks": post_dicts,
                          "post_checks_ok": all(c.ok for c in post_checks)},
                "arrays": {**{k: cols[k] for k in
                              ("bow_max", "bow_rms", "min_sep", "nnn_gap", "pe")},
                          "dl_flat": dl_all, "dtheta_flat": dtheta_all}}

    phases = None
    if init == "kink":
        phases = [RUN.Phase("release", n_prod, sample_every=sample_every, collect=True,
                            expect_steady=False,
                            note="방출 순간부터 기록 — 굽음이 변하는 것이 물리(관측 "
                                 "대상), 표류 검사를 걸면 발견을 경고로 부른다 (규칙 7')")]

    return RUN.Build(
        sim=sim, forces=[tab, wca], n_particles=n,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq if init != "kink" else 0, n_prod=n_prod,
        sample_every=sample_every, phases=phases or [],
        tags=["2D", "dlvo_secondary_minimum", "WCA_core", "no_bending", "no_friction",
             "no_drive", "chain", init],
        physical={"n_beads": n, "init": init, "kink_angle": kink_angle},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# 시각화 — 결과는 그래프로 보여준다 (CLAUDE.md)
# ════════════════════════════════════════════════════════════════════════
def make_plots(sys_, lg, n, init, kink_angle, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]     # ★ 라벨은 영어 (CLAUDE.md)
    matplotlib.rcParams["axes.unicode_minus"] = False

    res = np.load(outdir / "observables.npz")
    m = json.loads((outdir / "metrics.json").read_text())
    D = lg.derived
    dt_over_tau_bond = lg.ratio("times", "dt", "tau_bond")     # ★ 이 계의 자연 시간축
    t = np.arange(len(res["bow_rms"])) * (m["numerics"].get("sample_every", 1)) * dt_over_tau_bond

    fig, ax = plt.subplots(2, 3, figsize=(16, 9))

    ax[0, 0].plot(t, res["bow_rms"], "-", lw=1.0, label="bow rms (body frame)")
    ax[0, 0].plot(t, res["bow_max"], "-", lw=0.6, alpha=0.6, label="bow max")
    ax[0, 0].set(xlabel=r"$t / \tau_{bond}$", ylabel="bow [d]",
                title=f"Shape relaxation — n={n}, init={init}"
                      + (f", kink={kink_angle:.2f} rad" if kink_angle else ""))
    ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

    ax[0, 1].plot(t, res["pe"], "-", lw=0.8)
    ax[0, 1].set(xlabel=r"$t / \tau_{bond}$", ylabel=r"$\langle U \rangle / N$ [kT]",
                title="Potential energy per particle")
    ax[0, 1].grid(alpha=.3)

    extra = m.get("result", {})
    dtheta_var = extra.get("dtheta_var", float("nan"))
    dl_var = extra.get("dl_var", float("nan"))
    dl_var_pred = extra.get("dl_var_pred", float("nan"))
    k_bond_star = extra.get("k_bond_star", 1.0)
    if "dl_flat" in res:
        ax[1, 0].hist(res["dl_flat"], bins=60, density=True, alpha=.45,
                      color="tab:blue", label="measured (all bonds x samples)")
    harm_sigma = math.sqrt(1.0 / k_bond_star)
    boltz_sigma = math.sqrt(dl_var_pred) if dl_var_pred > 0 else harm_sigma
    xs2 = np.linspace(-5 * boltz_sigma, 5 * boltz_sigma, 400)
    ax[1, 0].plot(xs2, np.exp(-xs2 ** 2 / (2 * boltz_sigma ** 2)) / math.sqrt(2 * math.pi) / boltz_sigma,
                 "k-", lw=1.4, label=f"predicted (basin Boltzmann integral)  σ={boltz_sigma:.4f} d")
    ax[1, 0].plot(xs2, np.exp(-xs2 ** 2 / (2 * harm_sigma ** 2)) / math.sqrt(2 * math.pi) / harm_sigma,
                 "k:", lw=1.0, alpha=.6, label=f"harmonic approx (NOT the prediction)  σ={harm_sigma:.4f} d")
    ax[1, 0].axvline(0, color="gray", lw=.5)
    ax[1, 0].set(xlabel=r"$\delta\ell^* = \ell^* - \ell^*_{nat}$  [d]", ylabel="density",
                title=f"Bond-length golden test — measured Var={dl_var:.3e} d²"
                      f" vs predicted {dl_var_pred:.3e} d²")
    ax[1, 0].legend(fontsize=7); ax[1, 0].grid(alpha=.3)

    ax[1, 1].plot(t, res["min_sep"], "-", lw=.8, label="min bond length")
    ax[1, 1].axhline(1.0 + D["h_min_star"], color="tab:green", ls=":", label="natural length")
    ax[1, 1].plot(t, res["nnn_gap"] + 1.0, "-", lw=.6, alpha=.7, label="min NNN separation")
    ax[1, 1].axhline(1.0 + float(CUTOFF_H_STAR), color="tab:red", ls=":", label="DLVO cutoff")
    ax[1, 1].set(xlabel=r"$t / \tau_{bond}$", ylabel="r* [d]",
                title="Bond safety / non-adjacent folding")
    ax[1, 1].legend(fontsize=8); ax[1, 1].grid(alpha=.3)

    # ★ 굽힘(각) 자유도 — G1 이 예측이 없다고 말하는 바로 그 관측량. 사전 예측이
    # 없으므로 여기엔 대조선을 긋지 않는다(measurement) — 대신 참고용으로 결합
    # 신장의 golden-test σ(harmonic)와 시각적으로 비교할 수 있게 같은 각도 스케일로
    # 보여준다: 굽힘의 요동 폭이 신장의 요동 폭보다 압도적으로 크다면(래디안 vs
    # 무차원 길이라 직접 비교는 불가하지만) 그 자체가 "선형강성이 없다"는 정성적
    # 신호다 — 하한이 없다는 것은 계에서 정의되지 않은 척도로 잡히지 않는다는 뜻.
    if "dtheta_flat" in res and res["dtheta_flat"].size:
        ax[0, 2].hist(res["dtheta_flat"], bins=60, density=True, color="tab:orange", alpha=.6)
    ax[0, 2].axvline(0, color="gray", lw=.5)
    ax[0, 2].set(xlabel=r"$d\theta_i$ [rad] (local bend, all internal beads x samples)",
                ylabel="density",
                title=f"Bend-angle fluctuation — measured $\\langle d\\theta^2\\rangle$="
                      f"{dtheta_var:.3e} rad$^2$  (G1: no prior prediction)")
    ax[0, 2].grid(alpha=.3)

    ax[1, 2].axis("off")
    ax[1, 2].text(0.02, 0.95,
                  "G1 check (this run)\n"
                  "─────────────────\n"
                  f"bond stretch  Var={dl_var:.3e} d²\n"
                  f"  predicted   ={dl_var_pred:.3e} d²  (basin Boltzmann)\n"
                  f"  err         ={100*(dl_var/dl_var_pred-1):+.1f}%"
                  f"  (implementation_check)\n\n"
                  f"bend angle    Var={dtheta_var:.3e} rad²\n"
                  f"  no prior prediction (measurement)\n"
                  f"  naive equipartition kappa_eff="
                  f"{(1.0/dtheta_var if dtheta_var > 0 else float('nan')):.3g} kT/rad²\n"
                  f"  (reference number only — not a claim of harmonicity)\n\n"
                  + (f"bow: initial={extra.get('bow_initial', float('nan')):.4f} d, "
                     f"final(rms, last decile)={extra.get('bow_final_rms', float('nan')):.4f}"
                     f" +/- {extra.get('bow_final_sem', float('nan')):.4f} d\n"
                     if init == "kink" else ""),
                  transform=ax[1, 2].transAxes, fontsize=9, va="top", family="monospace")

    fig.suptitle(f"{sys_['label']}  n={n}  init={init}", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--init", choices=("straight", "kink"), required=True)
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--n", type=int, default=None, help="비드 수 (기본: system.yaml, n=9)")
    ap.add_argument("--kink-angle", type=float, default=0.30, help="rad (--init kink 전용)")
    ap.add_argument("--tobs", type=float, default=None,
                    help="관측창 (τ_bond 배수). 기본: straight=2e4, kink=3e3")
    ap.add_argument("--cycles", type=float, default=None)   # 미사용 — 인터페이스 대칭용
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--dt-scale", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=2000)
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--eq-scale", type=float, default=200.0,
                    help="(--init straight) 평형화 = 이 값 × τ_bond/dt")
    ap.add_argument("--smoke", action="store_true", help="빠른 정합성 확인용 — 스텝 대폭 축소")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/chain-relax-2d-dlvo/system.yaml")
    n = args.n if args.n is not None else sys_["n_list"][0]
    kink_angle = args.kink_angle if args.init == "kink" else 0.0

    if args.smoke:
        args.tobs = min(args.tobs or 50.0, 50.0)
        args.samples = min(args.samples, 100)
        args.eq_scale = min(args.eq_scale, 20.0)

    if not (args.report or args.spec or args.run):
        print("무엇을 할지 고르세요 — `--report` · `--spec` · `--run`")
        return 3

    return emit(sys_, n, args.init, kink_angle, args)


if __name__ == "__main__":
    sys.exit(main())
