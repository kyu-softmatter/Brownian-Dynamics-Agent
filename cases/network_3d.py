"""`network` — 3D 콜로이드 네트워크. **1단계: 압축 겔화** (CLAUDE.md 규칙 8).

스케치(intake/network/)는 "colloidal network / identical bead / DLVO, & JKR /
move one particle and see the network behaviors / G'(ω), G''(ω) / response time,
stress propagation / x(t) = A sin(ωt)" 라고만 적혀 있고 **수치가 하나도 없습니다.**
물성은 전부 `chain-bend-2d-dlvo` 승계입니다 (system.yaml 의 tier 참조).

이 파일이 하는 것 — **네트워크를 만드는 것까지만** 입니다:
    ① 3D 박스(φ₀=0.02)에 겹침 없이 흩뿌린다
    ② DLVO 2차극소로 응집시킨다 (페어 퍼텐셜만 — 위상 선언 없음)
    ③ 박스를 178단계로 등방 압축해 φ=0.10 에 도달한다
    ④ 압축 후 다시 이완시키고 **구조**를 잰다
        z · 독립 고리 수 · 자유단 비율 · 침투 · d_f · g(r) · 붕괴 결합 비율

구동(x(t)=A sin(ωt))·G'(ω) 는 **2단계**이고 여기 없습니다. 왜 나누는가 —
  · 규칙 8: 정적인 계를 먼저 세우고 움직임을 그 다음에 얹는다
  · N5(observation): HOOMD 의 bond/angle 은 **정적 토폴로지**라 응집 중에 선언할 수
    없다. 겔화가 끝난 뒤 접촉 목록을 뽑아 위상을 동결해야 2단계가 가능하다
  · JKR 굽힘도 그래서 2단계다 (`interactions[1].enabled_stage: 2`)

★★ 압축의 핵심 제약 (실행 전 실측 — scratch/verify_3d_boxresize.py):
    `update.BoxResize` 는 좌표를 **아핀 스케일**한다(오차 8.9e-16) → 이미 결합된 쌍의
    결합길이도 같이 줄어든다. 트리거당 선형변형이 (h_min−h_barrier)/ℓ = 0.703% 를
    넘으면 쌍이 장벽 안쪽으로 밀려 **1차극소(접촉)로 비가역 붕괴**한다.
    실측: 0.40% 유지(h*=0.007591) / 0.80% 붕괴(h*=−0.042257).
    → 0.4%/단계를 쓰고, `Interpolate` 가 L 에 선형이라 최악이 마지막 단계이므로
      단계 수는 134 가 아니라 **178** 이다 (system.yaml geometry.n_stages).

실행:
    $PY cases/network_3d.py --report                 # L3 리포트만
    $PY cases/network_3d.py --spec                   # 스펙 저장 (실행 안 함)
    $PY cases/network_3d.py --n 512 --stage-tau 1e-3 # 파일럿 실행
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
from bdbot import nondim as ND, run as RUN, scales as SC, sim as SIM        # noqa: E402
from bdbot.provenance import load_node                                      # noqa: E402
from bdbot.units import Q                                                   # noqa: E402

# DLVO 식은 chain-bend-2d-dlvo 에서 SI로 검증된 것을 **그대로** 쓴다 (두 번 적지 않는다)
from chain_bend_dlvo_2d import (                                            # noqa: E402
    SIGMA_CORE_STAR, build_table_arrays, dlvo_reduced_params, find_well, U_star,
)

ROOT = Path(__file__).resolve().parent.parent
CUTOFF_H_STAR = 0.06          # 페어 표 컷오프 (표면간극 h/d). 2차극소(0.0076)의 8배 밖
EPS_PER_STAGE = 0.004         # 압축 트리거당 허용 선형변형 (실측 안전선)
R_WCA = 2 ** (1 / 6)


# ════════════════════════════════════════════════════════════════════════
# ① 물리계 (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    P = load_node
    dlvo, jkr = raw["interactions"][0], raw["interactions"][1]
    g = raw["geometry"]
    return {
        "label": raw["label"], "dim": int(raw["dimensions"]),
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "psi0": P(raw["particle"]["surface_potential"]),
        "n_list": [int(x) for x in raw["particle"]["count"]["value"]],
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "eps_r": P(raw["medium"]["relative_permittivity"]),
        "ionic_strength": P(raw["medium"]["ionic_strength"]),
        "A_H": P(dlvo["hamaker_constant"]),
        "kappa_theta": P(jkr["angle_stiffness"]),
        "phi0": float(g["volume_fraction_initial"]["value"]),
        "phi1": float(g["volume_fraction_final"]["value"]),
        "n_stages": int(g["n_stages"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


def bond_edge_h_star(p: dict, w: dict) -> float:
    """결합 판정 경계 — 우물 외곽에서 U 가 깊이의 **절반**으로 돌아오는 h*.

    임의 컷오프를 쓰면 z·고리 수·자유단 비율이 그 선택에 흔들린다. 우물 깊이의
    절반은 원장에서 나오는 값이라 케이스마다 재현 가능하다.
    """
    h = np.geomspace(w["h_min"], CUTOFF_H_STAR, 200_000)
    U = U_star(h, p)
    return float(h[int(np.argmin(np.abs(U - w["U_min"] / 2.0)))])


def box_star(n: int, phi: float) -> float:
    """φ 에서 3D 정육면체 한 변 (d 단위). N(π/6)d³ / L³ = φ."""
    return (n * math.pi / (6.0 * phi)) ** (1.0 / 3.0)


# ════════════════════════════════════════════════════════════════════════
# ② 스케일 원장
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, n: int, *, dt_scale=1.0, stage_tau=1e-3,
                 agg_tau=0.2, post_tau=0.2, init="scatter", max_coord=4,
                 n_seeds=1, loop_bias=1) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B, tau_p = b["kT"], b["gamma"], b["tau_B"], b["tau_p"]

    p = dlvo_reduced_params(sys_)
    w = find_well(p)
    h_edge = bond_edge_h_star(p, w)
    ell = d * (1.0 + w["h_min"])
    k_bond = Q(w["k_bond_star"], "dimensionless") * kT / d ** 2
    sigma_bond = (kT / k_bond) ** 0.5
    lam_D = (Q(1.0, "dimensionless") / (p["kappa_star"] / d)).to("nm")

    tau_bond = C.relaxation_time(gamma, k_bond)
    # ★ 결합 탈출시간 — 우물이 11.7 kT 라 **가역**이다. 압축 준정적성의 기준 시간이
    #   τ_B 가 아니라 이것이다 (Kramers 어림: 전인자 τ_bond, 지수 |U_min|/kT).
    tau_esc = tau_bond * math.exp(abs(w["U_min"]))

    # dt 는 가장 빠른 모드가 정한다. 1단계 후보는 결합 신축뿐(굽힘·트랩은 2단계).
    # ⚠️ 붕괴한 쌍이 앉는 자리(h*<0, WCA↔vdW 균형)는 더 뻣뻣하다 — 붕괴가 생기면
    #    dt 를 다시 봐야 하므로 crushed_bond_fraction 을 관측량으로 낸다.
    tau_fast = tau_bond
    dt = dt_scale * C.dt_from_gate(tau_fast)

    # ★ `init="sprout"` — 목표 φ 에서 직접 망을 만든다 → 응집도 압축도 없다.
    #   L0 = L1 이고 단계 수 0. (사용자 제안 2026-08-06: 압축이 어려우면 중심 입자에서
    #   임의 방향으로 돋아나게 하자. 압축이 어려운 것은 아니었지만, φ₀=0.02 에서
    #   응집이 거의 진행되지 않는 것을 실측했으므로 이 경로가 더 낫다.)
    sprout = (init == "sprout")
    L1 = box_star(n, sys_["phi1"])
    L0 = L1 if sprout else box_star(n, sys_["phi0"])
    n_stages = 0 if sprout else sys_["n_stages"]
    r_cut = 1.0 + CUTOFF_H_STAR
    # 압축 단계당 물리 시간 (stage_tau 는 τ_B 단위) → 관측창
    T_stage = Q(stage_tau, "dimensionless") * tau_B
    T_compress = n_stages * T_stage
    T_obs = (Q(0.0 if sprout else agg_tau, "dimensionless") * tau_B + T_compress
             + Q(post_tau, "dimensionless") * tau_B)

    lg = SC.ScaleLedger(ref=SC.thermal_reference(
        d, kT, tau_B,
        rationale="d·kT·τ_B 열 기준. HOOMD 에서 σ=kT=γ=1 이 되고, 압축은 박스만 바뀐다."))
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("ell", ell, "DLVO 결합길이 (중심간, 자연장)", star=True)
    lg.add_length("h_min", d * w["h_min"], "2차극소 표면간극")
    lg.add_length("h_edge", d * h_edge, "결합 판정 경계 (U=well/2)")
    lg.add_length("sigma_bond", sigma_bond.to("nm"), "결합 열요동 폭")
    lg.add_length("lambda_D", lam_D, "Debye 길이")
    lg.add_length("r_cut", d * r_cut, "페어 표 컷오프")
    lg.add_length("L0", d * L0, f"초기 박스 (φ={sys_['phi0']:.3f})")
    # ★ role="box" 는 **압축 후** 박스로 잡는다 — 기하 검사(r_cut ≤ L/2)의 최악 조건이다
    lg.add_length("L", d * L1, f"박스 한 변 (압축 후, φ={sys_['phi1']:.3f})", role="box")

    lg.add_time("tau_B", tau_B, "브라운 시간 (기준)")
    lg.add_time("tau_p", tau_p, "관성 이완", role="inertia")
    lg.add_time("tau_bond", tau_bond, "★ DLVO 결합 신축 — 최속 모드", star=True)
    lg.add_time("tau_esc", tau_esc, "★ 결합 탈출 (가역 겔의 구조 이완)", star=True)
    lg.add_time("T_stage", T_stage, "압축 한 단계의 물리 시간")
    lg.add_time("T_compress", T_compress if n_stages else dt,
                f"압축 전체 ({n_stages} 단계)" if n_stages else "압축 없음 (sprout)")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("T_obs", T_obs, "관측창 (응집+압축+후이완)", role="observation")

    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("well_depth", Q(-w["U_min"], "dimensionless") * kT, "2차극소 깊이 (부호 반대)")
    lg.add_energy("barrier", Q(w["barrier_U"], "dimensionless") * kT, "DLVO 장벽")
    lg.add_energy("k_bond_d2", k_bond * d ** 2, "결합 신축 강성 × d²")

    total_strain = 1.0 - L1 / L0
    # 압축이 없으면 아핀 변형이 0 이므로 결합을 부술 위험 자체가 없다
    eps_max = 0.0 if n_stages == 0 else (L0 / L1 - 1.0) / n_stages
    eps_crit = (w["h_min"] - w["barrier_h"]) / (1.0 + w["h_min"])

    lg.derived.update(
        reduced=p, well=w, h_edge=h_edge, ell=ell, k_bond=k_bond, tau_fast=tau_fast,
        tau_esc=tau_esc, tau_B=tau_B, T_obs=T_obs, T_stage=T_stage,
        L0_star=L0, L1_star=L1, r_cut_star=r_cut, r_bond_star=1.0 + h_edge,
        total_strain=total_strain, eps_max=eps_max, eps_crit=eps_crit, init=init,
        sprout=sprout, max_coord=max_coord, n_seeds=n_seeds, loop_bias=loop_bias,
        n_stages=n_stages, stage_tau=stage_tau,
        agg_tau=0.0 if sprout else agg_tau, post_tau=post_tau,
        U_min_star=w["U_min"], barrier_star=w["barrier_U"], h_min_star=w["h_min"],
        k_bond_star=w["k_bond_star"])
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ③ 무차원수 + 검사
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg, n):
    D = lg.derived
    r = lg.ratio
    groups = [
        ND.Group("phi_final", sys_["phi1"], None, None,
                 "N(π/6)d³/L³", "부피분율 — 성긴 겔"),
        ND.Group("well_depth/kT", r("energies", "well_depth", "kT"),
                 ("energies", "well_depth"), ("energies", "kT"), "",
                 "★ 11.7 → 결합이 가역 (τ_esc = 0.115 τ_B)"),
        ND.Group("barrier/kT", r("energies", "barrier", "kT"),
                 ("energies", "barrier"), ("energies", "kT"), "",
                 "1차극소 전이 사실상 불가능"),
        ND.Group("k_bond_star", r("energies", "k_bond_d2", "kT"),
                 ("energies", "k_bond_d2"), ("energies", "kT"), "",
                 "결합 신축 강성 — dt 를 정한다"),
        ND.Group("tau_bond/tau_B", r("times", "tau_bond", "tau_B"),
                 ("times", "tau_bond"), ("times", "tau_B"), "",
                 "★ 척도 분리 6자릿수 — dt 비용의 근원"),
        ND.Group("tau_esc/tau_B", r("times", "tau_esc", "tau_B"),
                 ("times", "tau_esc"), ("times", "tau_B"), "",
                 "결합 재배열 / 확산 — 압축 준정적성의 기준"),
        ND.Group("T_stage/tau_esc", r("times", "T_stage", "tau_esc"),
                 ("times", "T_stage"), ("times", "tau_esc"), "",
                 "★ ≥1 이면 준정적 압축, ≪1 이면 아핀 압축 (다른 물리)"),
        ND.Group("ell/d", r("lengths", "ell", "d"),
                 ("lengths", "ell"), ("lengths", "d"), "", "결합길이 / 지름"),
        ND.Group("L/d", r("lengths", "L", "d"),
                 ("lengths", "L"), ("lengths", "d"), "", "압축 후 박스 (전파 유효반경 L/2)"),
        ND.Group("eps_max", D["eps_max"], None, None,
                 "(L0/L1−1)/n_stages", "★ 압축 단계당 최대 선형변형"),
    ]
    checks = [
        C.Check("integration", "dt / tau_fast", lg.ratio("times", "dt", "tau_bond"), C.GATE,
                note="Brownian 은 O(δt) — 최속 모드는 결합 신축"),
        # ★ 소프트 처리는 `chain-bend-2d-dlvo` 의 확립된 관례를 그대로 따른다 —
        #   같은 입자·같은 DLVO 결합이라 이 비는 **글자 그대로 같은 값(0.0282)** 이고,
        #   그 케이스가 이미 같은 판단을 문서화했다. 검사를 없애지 않고 남겨 보고한다.
        C.Check("model", "note: tau_p/tau_bond", lg.ratio("times", "tau_p", "tau_bond"),
                C.GATE, "<=",
                "★ 결합 우물이 깊고 좁아 τ_bond 가 τ_p 에 근접한다 (관성 무시 기준을 "
                "2.82배 넘김 — N·φ·압축과 무관한 **결합 물리 자체의 성질**). "
                "다만 ζ = ½(τ_p/τ_bond)^(−½) = 2.98 로 과감쇠 조건(ζ>1)은 여유 있게 "
                "만족한다. chain-bend-2d-oscill 에서 훨씬 심한 ζ=0.65(τ_p/τ_fast=0.60)를 "
                "OverdampedViscous vs Langevin(kT=0) 로 대조해 관측량 영향이 최대 "
                "0.159% 임을 실측했다 — 여기는 그보다 21배 유리하다. "
                "⚠️ 이 계에서 그 대조를 직접 한 것은 아니다 (not_verified 에 기록)",
                hard=False),
        C.Check("model", "결합 안정  σ_bond/h_min",
                lg.ratio("lengths", "sigma_bond", "h_min"), 0.5, "<=",
                "★ 결합 열요동 폭이 우물 위치의 절반을 넘으면 열적으로 자꾸 깨진다 — "
                "그 자체가 결과일 수 있으나 사전에 표시 (chain-bend 와 같은 검사)",
                hard=False),
        C.Check("geometry", "r_cut / (L/2)",
                lg.ratio("lengths", "r_cut", "L") * 2.0, 1.0,
                note="★ 압축 **후** 박스로 판정 (bd-hoomd 함정 6)"),
        C.Check("geometry", "lambda_D / (L/4)",
                lg.ratio("lengths", "lambda_D", "L") * 4.0, 1.0,
                note="이중층이 박스에 들어가는가"),
        C.Check("geometry", "ell / (L/2)", lg.ratio("lengths", "ell", "L") * 2.0, 1.0,
                note="결합길이가 최소이미지 안"),
        # ★★ 이 케이스의 핵심 검사 — 실측 문턱을 설계에 강제한다
        C.Check("integration", "eps_max / eps_crit", D["eps_max"] / D["eps_crit"], 1.0,
                note="★ 압축이 결합을 부수지 않는 조건. 실측 0.40% 유지/0.80% 붕괴"
                     + ("  (sprout — 압축이 없어 아핀 변형 0)" if D["sprout"] else "")),
        C.Check("statistics", "N", float(n), 512.0, op=">=", hard=False,
                note="구조 통계 — 파일럿 512, 생산 1528"),
        C.Check("statistics", "T_stage / tau_bond",
                lg.ratio("times", "T_stage", "tau_bond") if not D["sprout"] else 1e9,
                10.0, op=">=", hard=False,
                note="단계당 결합이 이완할 시간. 이게 1 미만이면 붕괴 판정이 무의미"
                     + ("  (sprout — 해당 없음)" if D["sprout"] else "")),
        C.Check("finite-size", "L / ell", lg.ratio("lengths", "L", "ell"), 8.0, op=">=",
                hard=False, note="박스가 결합길이의 몇 배 — 메시가 들어가는가"),
    ]
    return groups, checks


def report_blocks(sys_, lg, n, n_agg, n_stage, n_post):
    D = lg.derived
    inp = [R.kv("d", f"{sys_['d'].value:~.4gP}", sys_["d"].tier, sys_["d"].source[:44]),
           R.kv("T", f"{sys_['T'].value:~.4gP}", sys_["T"].tier, "물"),
           R.kv("eta", f"{sys_['eta'].value:~.4gP}", sys_["eta"].tier, "물@300K"),
           R.kv("A_H", f"{sys_['A_H'].value:~.4gP}", sys_["A_H"].tier, sys_["A_H"].source[:44]),
           R.kv("I", f"{sys_['ionic_strength'].value:~.4gP}",
                sys_["ionic_strength"].tier, "MgCl2"),
           R.kv("N", str(n), 3, "★제안 — 그림의 21개는 도식 (A6)"),
           R.kv("phi", f"{sys_['phi0']:.3f} → {sys_['phi1']:.3f}", 1, "압축 겔화")]
    der = [f"  결합: 장벽 {D['barrier_star']:.2f} kT   2차극소 {D['U_min_star']:.3f} kT"
           f" @ h*={D['h_min_star']:.6f}",
           f"  결합 판정 경계 r*<{D['r_bond_star']:.6f} (U=well/2, h*={D['h_edge']:.6f})",
           f"  k_bond* {D['k_bond_star']:.4g} kT/d²  → τ_bond {D['tau_fast']:~.4gP}",
           f"  τ_esc {D['tau_esc']:~.4gP} = {float((D['tau_esc']/D['tau_B']).to('')):.4f} τ_B"
           f"  ← 가역 겔의 구조 이완 시간",
           f"  박스 {D['L0_star']:.3f} d → {D['L1_star']:.3f} d"
           f"   (총 선형변형 {D['total_strain']*100:.2f}%)",
           (f"  위상 생성기 **sprout** — 압축 없음 (아핀 변형 0, 결합 붕괴 위험 0). "
            f"문턱 {D['eps_crit']*100:.4f}% 는 해당 없음"
            if D["eps_max"] == 0 else
            f"  압축 {D['n_stages']} 단계 · 단계당 최대변형 {D['eps_max']*100:.4f}%"
            f"  (문턱 {D['eps_crit']*100:.4f}%, 여유 {D['eps_crit']/D['eps_max']:.2f}×)")]
    if D["sprout"]:
        plan = [f"  ① 돋아나기  (MD 아님)   씨앗→임의방향, ℓ 에 정확히 배치, "
                f"max_coord={D['max_coord']}, seeds={D['n_seeds']}, loop_bias={D['loop_bias']}",
                f"  ② 이완   {n_post:>12,} 스텝  ({D['post_tau']:.4g} τ_B, φ={sys_['phi1']:.3f})",
                f"  합계     {n_post:>12,} 스텝",
                "  ★ 압축 없음 → 아핀 변형 0, 결합 붕괴 위험 0",
                "  ⚠️ 위상은 **부과된 것**이다 (DLVO 동역학이 만든 것이 아님) — "
                "압축 경로와 비교해야 위상 의존성이 보인다 (원칙 7')"]
    else:
        plan = [f"  ① 응집   {n_agg:>12,} 스텝  ({D['agg_tau']:.4g} τ_B, φ={sys_['phi0']:.3f} 고정)",
                f"  ② 압축   {n_stage*D['n_stages']:>12,} 스텝  "
                f"({D['n_stages']} × {n_stage:,} = {D['stage_tau']*D['n_stages']:.4g} τ_B)",
                f"  ③ 후이완 {n_post:>12,} 스텝  ({D['post_tau']:.4g} τ_B, φ={sys_['phi1']:.3f})",
                f"  합계     {n_agg+n_stage*D['n_stages']+n_post:>12,} 스텝",
                f"  ★ T_stage/τ_esc = {D['stage_tau']/float((D['tau_esc']/D['tau_B']).to('')):.4g}"
                f"  ({'준정적' if D['stage_tau'] >= float((D['tau_esc']/D['tau_B']).to('')) else '아핀 압축 — 확산 미개입'})"]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# L3 — 스펙
# ════════════════════════════════════════════════════════════════════════
def build_spec(sys_, n, args):
    lg = build_ledger(sys_, n, dt_scale=args.dt_scale, stage_tau=args.stage_tau,
                      agg_tau=args.agg_tau, post_tau=args.post_tau,
                      init=args.init, max_coord=args.max_coord,
                      n_seeds=args.n_seeds, loop_bias=args.loop_bias)
    D = lg.derived
    dt = lg.get("times", "dt")
    per = lambda tau: max(1, int(round(float((Q(tau, "dimensionless") * D["tau_B"] / dt).to("")))))
    n_agg = 0 if D["sprout"] else per(args.agg_tau)
    n_stage, n_post = per(args.stage_tau), per(args.post_tau)
    n_prod = n_agg + n_stage * D["n_stages"] + n_post
    sample_every = max(1, n_prod // args.samples)

    groups, checks = analyze_scales(sys_, lg, n)
    # ★ tag 는 물리를 갈라야 한다 — init 이 위상 생성기를 바꾸므로 반드시 들어간다
    if D["sprout"]:
        tag = (f"N{n}-sprout-mc{args.max_coord}-sd{args.n_seeds}"
               f"-lb{args.loop_bias}-po{args.post_tau:g}")
    else:
        tag = f"N{n}-st{args.stage_tau:g}-ag{args.agg_tau:g}"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"
    p = D["reduced"]
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"n_particles": n, "phi0": sys_["phi0"], "phi1": sys_["phi1"],
                "L0_star": D["L0_star"], "L1_star": D["L1_star"],
                "n_stages": D["n_stages"], "eps_max": D["eps_max"], "eps_crit": D["eps_crit"],
                "kappa_star": p["kappa_star"], "edl_amp": p["edl_amp"],
                "vdw_amp": p["vdw_amp"], "a_star": p["a_star"],
                "cutoff_h_star": CUTOFF_H_STAR,
                "h_min_star": D["h_min_star"], "h_edge_star": D["h_edge"],
                "r_bond_star": D["r_bond_star"],
                "well_depth_star": D["U_min_star"], "barrier_star": D["barrier_star"],
                "k_bond_star": D["k_bond_star"],
                "tau_esc_star": float((D["tau_esc"] / D["tau_B"]).to("")),
                "stage_tau": args.stage_tau, "agg_tau": D["agg_tau"],
                "post_tau": args.post_tau,
                # ★ 위상 생성기 — 물리를 정하므로 반드시 해시에 들어간다
                "init": D["init"], "max_coord": int(args.max_coord),
                "n_seeds": int(args.n_seeds), "loop_bias": int(args.loop_bias),
                "stage": 1, "jkr": False, "drive": "none"},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "n_eq": 0, "n_prod": n_prod, "n_agg": n_agg,
                  "n_stage": n_stage, "n_post": n_post,
                  "n_samples": args.samples, "sample_every": sample_every,
                  "seed": args.seed},
        tag=tag, nhex=12)
    return lg, spec, groups, checks, n_agg, n_stage, n_post


def emit(sys_, n, args) -> int:
    lg, spec, groups, checks, n_agg, n_stage, n_post = build_spec(sys_, n, args)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n, n_agg, n_stage, n_post)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']} 1단계(겔화)  N={n}"
              f"  T_stage={args.stage_tau:g}τ_B   run_id={run_id}",
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

    # `result.txt` = 프로젝트의 완료 마커. 케이스 스크립트의 책임이다 (안 쓰면
    # status 가 런을 0개로 세고, 정리 스크립트가 완료 런을 지운다)
    if v["status"] != "skipped":
        lines = []
        try:
            mj = json.loads((outdir / "metrics.json").read_text())
            for o in mj.get("observables", []):
                m, pr = o.get("measured"), o.get("predicted")
                tail = f"   (예측 {pr:.6g})" if isinstance(pr, (int, float)) else ""
                lines.append(f"  {o['name']:<26} {m:.6g}{tail}" if m is not None
                             else f"  {o['name']:<26} —")
        except Exception as e:                        # noqa: BLE001
            lines.append(f"  (metrics.json 을 읽지 못함: {e})")
        (outdir / "result.txt").write_text(
            report + "\n" + "\n".join(["=" * 84, f"결과 — {run_id}", "=" * 84,
                                       *lines, "=" * 84, verdict_txt]))
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


# ════════════════════════════════════════════════════════════════════════
# L4 — 구조 분석 (순수 함수. 런과 무관하게 시험 가능)
# ════════════════════════════════════════════════════════════════════════
def contacts(pos: np.ndarray, L: float, r_bond: float):
    """PBC 최소이미지로 결합 쌍 목록 + 각 쌍의 표면간극. O(N²) — N≲2000 이면 충분."""
    n = len(pos)
    d = pos[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d ** 2).sum(-1))
    iu = np.triu_indices(n, 1)
    rr = r[iu]
    m = rr < r_bond
    return np.column_stack([iu[0][m], iu[1][m]]), rr[m], rr


def topology(n: int, pairs: np.ndarray):
    """배위수·연결성분·독립 고리(Betti-1)·자유단."""
    deg = np.zeros(n, dtype=int)
    adj = [[] for _ in range(n)]
    for i, j in pairs:
        deg[i] += 1
        deg[j] += 1
        adj[i].append(j)
        adj[j].append(i)
    seen = np.zeros(n, dtype=bool)
    comps, sizes = 0, []
    for s in range(n):
        if seen[s]:
            continue
        comps += 1
        cnt, stack = 0, [s]
        while stack:
            u = stack.pop()
            if seen[u]:
                continue
            seen[u] = True
            cnt += 1
            stack.extend(v for v in adj[u] if not seen[v])
        sizes.append(cnt)
    e = len(pairs)
    return dict(z=float(2 * e / n), n_links=e, n_components=comps,
                loops=int(e - n + comps), dangling=float((deg == 1).sum() / n),
                isolated=float((deg == 0).sum() / n),
                largest_cluster=float(max(sizes) / n) if sizes else 0.0)


def percolates(pos: np.ndarray, L: float, pairs: np.ndarray) -> float:
    """PBC 스패닝 판정 — 결합을 따라가며 누적 이미지 오프셋이 0 이 아닌 고리를 찾는다.

    최소이미지 변위를 더해가며 BFS 하고, 이미 방문한 노드에 **다른 누적 변위**로
    도달하면 그 성분은 박스를 감싸며 자기 자신에 이어진다 = 침투.
    반환: 침투한 축 수 / 3.
    """
    n = len(pos)
    adj = [[] for _ in range(n)]
    for i, j in pairs:
        dv = pos[j] - pos[i]
        dv -= L * np.round(dv / L)
        adj[i].append((j, dv))
        adj[j].append((i, -dv))
    axes = np.zeros(3, dtype=bool)
    off = np.full((n, 3), np.nan)
    for s in range(n):
        if not np.isnan(off[s, 0]):
            continue
        off[s] = 0.0
        stack = [s]
        while stack:
            u = stack.pop()
            for v, dv in adj[u]:
                cand = off[u] + dv
                if np.isnan(off[v, 0]):
                    off[v] = cand
                    stack.append(v)
                else:
                    w = np.abs(cand - off[v])
                    axes |= w > 0.5 * L          # 감싼 축
    return float(axes.sum() / 3.0)


def fractal_dimension(pos: np.ndarray, L: float, rng) -> float:
    """질량-반지름 스케일링 d_f: N(r) ∝ r^d_f. 중심 여러 개 평균, r ∈ [1.5d, L/4]."""
    rs = np.geomspace(1.5, max(L / 4.0, 2.0), 12)
    cent = pos[rng.choice(len(pos), size=min(40, len(pos)), replace=False)]
    d = cent[:, None, :] - pos[None, :, :]
    d -= L * np.round(d / L)
    r = np.sqrt((d ** 2).sum(-1))
    cnt = np.array([(r < rr).sum(axis=1).mean() for rr in rs])
    ok = cnt > 1
    if ok.sum() < 4:
        return float("nan")
    return float(np.polyfit(np.log(rs[ok]), np.log(cnt[ok]), 1)[0])


def rdf(r_all: np.ndarray, n: int, L: float, r_max: float, nbins: int = 200):
    """g(r) — 전체 쌍거리 배열에서. r_max ≤ L/2."""
    h, edges = np.histogram(r_all, bins=nbins, range=(0.8, r_max))
    mid = 0.5 * (edges[1:] + edges[:-1])
    shell = 4.0 / 3.0 * math.pi * (edges[1:] ** 3 - edges[:-1] ** 3)
    rho = n / L ** 3
    return mid, h / (0.5 * n * rho * shell)


# ════════════════════════════════════════════════════════════════════════
# L4 — HOOMD 빌더
# ════════════════════════════════════════════════════════════════════════
def sprout_network(n: int, L: float, ell: float, rng, *, max_coord: int = 4,
                   n_seeds: int = 1, loop_bias: int = 1,
                   min_center: float | None = None,
                   r_bond: float = 1.0180308, tries: int = 3000):
    """★ 씨앗 입자에서 임의 방향으로 **돋아나게** 해서 망을 직접 만든다 (사용자 제안).

    압축 경로의 두 약점을 없앤다:
      · φ₀=0.02 에서 응집이 거의 진행되지 않는다 (실측: ⟨r_nn⟩ 1.7519→1.7534, **증가**)
      · 아핀 압축이 결합을 부술 위험 (ε_crit 문턱) — 여기선 압축이 아예 없다
    결합을 **정확히 자연장 ℓ 에** 놓으므로 초기 응력이 0 이고, 목표 φ 에서 바로 만든다.

    ⚠️ **인식론적 차이를 명시할 것** (원칙 7'): 압축 경로의 위상은 DLVO 동역학이 만든
    것이고, 이 경로의 위상은 **내가 부과한 것**이다. 그래서 이걸로 얻은 구조는
    `hypothesis` 가 아니라 **입력**이다 — 2단계의 G'(ω) 가 위상에 얼마나 의존하는지를
    보려면 두 생성기를 **둘 다** 돌려 비교해야 한다 (그게 이 방식의 진짜 값이다).

    ★ 고리(loop)에 대하여: 씨앗 하나에서 가지만 뻗으면 위상이 **트리**가 되고, 그건
    스케치 그림과 똑같은 실패 모드다 (고리 0 → G'(ω→0)→0, observation A6).
    그래서 새 입자를 **겹치지만 않으면** 놓는다 (min_center=1.0 = 접촉까지 허용) —
    가지가 서로 다가오면 자동으로 결합이 생기고 **고리가 닫힌다.** 고리 수는
    강제하지 않고 결과로 측정한다.

    `max_coord`: 이 배위수에 도달한 입자는 더 이상 부모로 뽑지 않는다 (Eden 성장류).
                 작으면 성기고 가지가 길다, 크면 조밀하다.
    `min_center`: 허용 최소 중심간 거리. **기본값 = ℓ (결합길이)** 이다.
        ⚠️ 처음에 1.0(접촉)으로 뒀더니 실측 min_sep 이 1.00018 까지 내려갔다 —
        DLVO 장벽(h*=0.000508)의 **안쪽**이라 그 쌍은 이완 첫 스텝에 1차극소(접촉)로
        비가역 붕괴한다. 압축 경로에서 ε_crit 으로 막았던 것과 **정확히 같은 사고**를
        초기 배치로 저지르는 셈이었다. ℓ 로 두면 모든 쌍이 우물 최소 이상이라
        초기 응력이 0 이고 장벽을 넘는 쌍이 없다.
    """
    if min_center is None:
        min_center = ell
    pos = np.zeros((n, 3))
    coord = np.zeros(n, dtype=int)
    # ★ 씨앗을 여러 개 둔다 — 하나면 성장이 한 덩어리라 **주기경계를 감는 하중경로가
    #   안 생긴다** (실측: N=512 에서 고리 14개가 전부 winding 0, percolation 0.00).
    #   씨앗을 흩뿌리면 서로 다른 덩어리의 가지가 경계를 가로질러 만나 감는 고리가 생긴다.
    k = 0
    while k < n_seeds and k < n:
        c = rng.uniform(-L / 2, L / 2, 3) if k else np.zeros(3)
        if k:
            dd = pos[:k] - c
            dd -= L * np.round(dd / L)
            if np.sqrt((dd ** 2).sum(-1)).min() < min_center:
                continue
        pos[k] = c
        k += 1
    stall = 0
    while k < n and stall < tries:
        cand_parents = np.flatnonzero(coord[:k] < max_coord)
        if len(cand_parents) == 0:               # 전부 포화 — 새 씨앗을 띄운다
            c = rng.uniform(-L / 2, L / 2, 3)
            dd = pos[:k] - c
            dd -= L * np.round(dd / L)
            if np.sqrt((dd ** 2).sum(-1)).min() < min_center:
                stall += 1
                continue
            pos[k] = c
            k += 1
            stall = 0
            continue
        par = int(rng.choice(cand_parents))
        # ★ 후보 방향을 loop_bias 개 시험해 **결합이 가장 많이 생기는** 자리를 고른다.
        #   왜 필요한가 (실측): 결합이 성립하는 창은 ℓ=1.00759 ~ r_bond=1.01803 으로
        #   **폭 0.0104 d** 뿐이라, 무작위 한 방향은 부모와의 결합 하나만 만들고
        #   다른 가지와 만날 확률이 사실상 0 이다 → 위상이 트리에 가까워지고
        #   (고리 14/512) **감는 하중경로가 안 생긴다** (percolation 0.00).
        #   MD 이완도 이걸 못 고친다 — 우물이 짧아 0.018 d 밖은 안 끌린다
        #   (실측: 이완 중 ⟨r_nn⟩ 1.0072→1.0069, 위상 거의 불변).
        #   loop_bias=1 이면 예전 동작 그대로다.
        best = None
        for _ in range(max(1, loop_bias)):
            v = rng.normal(size=3)
            v /= np.linalg.norm(v)
            c = pos[par] + ell * v
            c -= L * np.round(c / L)             # PBC
            dd = pos[:k] - c
            dd -= L * np.round(dd / L)
            r = np.sqrt((dd ** 2).sum(-1))
            if r.min() < min_center:             # 겹침 → 이 후보는 버린다
                continue
            nb = int((r < r_bond).sum())
            if best is None or nb > best[0]:
                best = (nb, c, r)
                if nb >= max_coord:              # 더 좋아질 여지 없음
                    break
        if best is None:
            stall += 1
            continue
        _, c, r = best
        pos[k] = c
        # 이 자리에서 실제로 결합이 되는 이웃 전부에 배위수를 준다 → 고리가 세어진다
        coord[k] = int((r < r_bond).sum())
        coord[:k][r < r_bond] += 1
        k += 1
        stall = 0
    if k < n:
        raise RuntimeError(f"돋아나기 실패: {k}/{n} 개만 놓임 (φ 가 너무 높거나 "
                           f"max_coord={max_coord} 가 너무 작다)")
    return pos


def scatter_no_overlap(n: int, L: float, rng, min_sep: float = 1.25, tries: int = 400):
    """겹침·즉시결합 없이 흩뿌린다 (min_sep > r_bond 여야 초기 z=0)."""
    pos = np.empty((n, 3))
    k = 0
    for _ in range(tries * n):
        if k == n:
            break
        c = rng.uniform(-L / 2, L / 2, 3)
        if k:
            d = pos[:k] - c
            d -= L * np.round(d / L)
            if np.sqrt((d ** 2).sum(-1)).min() < min_sep:
                continue
        pos[k] = c
        k += 1
    if k < n:
        raise RuntimeError(f"흩뿌리기 실패: {k}/{n} (φ 가 너무 높거나 min_sep 이 크다)")
    return pos


@RUN.builder("network")
def build(spec, outdir=None) -> RUN.Build:
    import gsd.hoomd
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["n_particles"])
    L0, L1 = float(P["L0_star"]), float(P["L1_star"])
    n_stages = int(P["n_stages"])
    r_bond = float(P["r_bond_star"])
    r_cut = 1.0 + float(P["cutoff_h_star"])
    dt, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_agg, n_stage, n_post = int(Nm["n_agg"]), int(Nm["n_stage"]), int(Nm["n_post"])
    rng = np.random.default_rng(seed)

    init = str(P.get("init", "scatter"))
    if init == "sprout":
        # ★ 목표 φ 에서 바로 만든다 → 응집도 압축도 없다 (L0 == L1 이 스펙에서 강제됨)
        pos0 = sprout_network(n, L1, 1.0 + float(P["h_min_star"]), rng,
                              max_coord=int(P.get("max_coord", 4)),
                              n_seeds=int(P.get("n_seeds", 1)),
                              loop_bias=int(P.get("loop_bias", 1)),
                              r_bond=r_bond)
    else:
        pos0 = scatter_no_overlap(n, L0, rng, min_sep=max(1.25, r_bond * 1.1))

    f = gsd.hoomd.Frame()
    f.particles.N = n
    f.particles.position = pos0
    f.particles.typeid = [0] * n
    f.particles.types = ["A"]
    f.particles.mass = [1.0] * n
    f.configuration.box = [L0, L0, L0, 0, 0, 0]      # ★ 3D: Lz=L (2D 는 0 — 함정 9)
    f.configuration.dimensions = 3
    sim = SIM.make_sim(f, seed=seed)

    cell = md.nlist.Cell(buffer=0.2)
    reduced = {k: P[k] for k in ("kappa_star", "edl_amp", "vdw_amp", "a_star")}
    r_min = 1.0 + 1e-4
    r_arr, U_arr, F_arr = build_table_arrays(reduced, r_min, r_cut)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut)
    tab.params[("A", "A")] = dict(r_min=r_min, U=U_arr, F=F_arr)
    # ★ WCA 컷오프가 정확히 r=d 에서 끝나도록 σ_c = d·2^(−1/6) (chain-bend 와 같은 관례)
    wca = md.pair.LJ(nlist=cell, default_r_cut=SIGMA_CORE_STAR * R_WCA, mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=SIGMA_CORE_STAR)
    forces = [tab, wca]

    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=dt, methods=[bd], forces=forces)
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    # ★★ 압축 — 응집 단계가 끝난 뒤에만 돈다. `Interpolate` 는 L 에 선형이라
    #   상대변형이 뒤로 갈수록 커진다 → 단계 수가 178 (134 아님, system.yaml 참조).
    #   `init="sprout"` 이면 이미 목표 φ 라 압축이 없다 (n_stages=0).
    if n_stages > 0:
        box_var = hoomd.variant.box.Interpolate(
            initial_box=hoomd.Box(Lx=L0, Ly=L0, Lz=L0),
            final_box=hoomd.Box(Lx=L1, Ly=L1, Lz=L1),
            variant=hoomd.variant.Ramp(A=0.0, B=1.0, t_start=n_agg,
                                       t_ramp=n_stage * n_stages))
        sim.operations.updaters.append(hoomd.update.BoxResize(
            trigger=hoomd.trigger.Periodic(period=n_stage, phase=n_agg), box=box_var))
    SIM.add_trajectory_writer(sim, (Path(outdir) / "traj_A.gsd") if outdir else None,
                              max(1, (n_agg + n_stage * n_stages + n_post) // 200))

    # ★★ 평형·FROZEN 판정 지표로 **⟨U⟩/N 을 쓸 수 없다.**
    #   응집 초기에는 어떤 쌍도 컷오프(1.06 d) 안에 없어 PE 가 **정확히 0** 이고,
    #   health 가 그걸 "열잡음이 있는데 변하지 않는다"→FROZEN 으로 판정한다. 오탐인데
    #   대가가 크다: `run.execute` 는 `verdict == OK` 일 때만 `finalize()` 를 부르므로
    #   **관측량이 통째로 안 나온다** (스모크 런에서 실측 — observables 0개, result None).
    #   대신 `abp-rod` 의 전례를 따른다: 계가 살아 있으면 **절대 정확히 상수가 안 되는**
    #   물리량을 쓴다 — 기준 입자 64개의 최근접 이웃거리 평균 ⟨r_nn⟩ [d].
    #     · 응집 → 감소 (결합 형성)   · 압축 → 감소   · 후이완 → 정상상태
    #     · 운동이 멈추면 정확히 고정 → FROZEN 이 **제대로** 잡힌다
    #   ⚠️ health 리포트는 이 값을 "⟨U⟩/N [kT]" 라고 라벨한다 — 단위는 d 다.
    #      PE/N 자체는 표본(`pe_energy`)으로 따로 남긴다.
    n_ref = min(64, n)
    ref = np.arange(n_ref)          # get_snapshot() 은 tag 순서다 (local_snapshot 과 다름)

    def r_nn_mean():
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, copy=True)
        L = float(snap.configuration.box[0])
        dd = pos[ref][:, None, :] - pos[None, :, :]
        dd -= L * np.round(dd / L)
        r = np.sqrt((dd ** 2).sum(-1))
        r[ref, ref] = np.inf                       # 자기 자신 제외
        return float(r.min(axis=1).mean())

    def pe_energy():
        return float(sum(np.array(fo.energies).sum() for fo in forces)) / n

    # ★ 파생량은 클로저에 누적하고 원시배열은 필요한 것만 — 표본마다 전체 N배열을
    #   돌려주면 observables.npz 가 수백배로 부푼다 (chain-bend 에서 448KB→148MB)
    series: dict = {k: [] for k in
                    ("t", "L", "phi", "z", "loops", "dangling", "isolated",
                     "largest_cluster", "n_components", "percolation", "d_f",
                     "crushed", "min_sep", "h_mean", "pe_energy", "r_nn")}
    last: dict = {}

    def sample(timestep, phase):
        snap = sim.state.get_snapshot()
        pos = np.array(snap.particles.position, copy=True)
        L = float(snap.configuration.box[0])
        pairs, h_pairs, r_all = contacts(pos, L, r_bond)
        top = topology(n, pairs)
        phi = n * math.pi / 6.0 / L ** 3
        crushed = float((h_pairs < 1.0).sum() / max(len(h_pairs), 1))   # r<d ⟺ h<0
        perc = percolates(pos, L, pairs) if len(pairs) else 0.0
        df = fractal_dimension(pos, L, rng) if top["largest_cluster"] > 0.1 else float("nan")
        row = dict(t=timestep * dt, L=L, phi=phi, percolation=perc, d_f=df,
                   crushed=crushed, min_sep=float(r_all.min()),
                   h_mean=float(h_pairs.mean() - 1.0) if len(h_pairs) else float("nan"),
                   pe_energy=pe_energy(), r_nn=r_nn_mean(),
                   **{k: top[k] for k in ("z", "loops", "dangling", "isolated",
                                          "largest_cluster", "n_components")})
        for k, v in row.items():
            series[k].append(float(v))
        last.clear()
        last.update(row, pos=pos, pairs=pairs, L=L)
        return {k: row[k] for k in ("z", "loops", "phi", "percolation", "crushed",
                                    "min_sep", "pe_energy", "r_nn")}

    def finalize(_res):
        pos, L = last["pos"], last["L"]
        _, _, r_all = contacts(pos, L, r_bond)
        mid, g = rdf(r_all, n, L, min(L / 2.0, 4.0))
        # ★ 첫 피크는 **따로 촘촘히** 잰다. 위 g(r) 는 빈 폭이 0.016 d 라 2차극소
        #   (h_min=0.0076 d) 를 분해하지 못해 예측 1.00759 대조가 무의미해진다.
        fine_r, fine_g = rdf(r_all, n, L, 1.10, nbins=400)     # 빈 폭 5e-4 d
        peak = float(fine_r[int(np.argmax(fine_g))]) if np.any(fine_g > 0) else float("nan")
        obs = [
            MET.observable(
                "phi_final", last["phi"], predicted=float(P["phi1"]),
                role="implementation_check", tol_pct=1.0,
                note="압축이 목표 부피분율에 정확히 도달했는가",
                derivation="φ = N(π/6)d³/L³ 은 **기하 항등식**이라 조합에서도 성립한다 — "
                           "물리 예측이 아니라 'BoxResize 가 목표 박스 L1 에 도달했는가' "
                           "라는 배선 검사다 (실측: Lx 오차 <1e-6, "
                           "scratch/verify_3d_boxresize.py ②③). 어긋나면 압축 배선 버그."),
            MET.observable(
                "crushed_bond_fraction", last["crushed"], predicted=0.0,
                sigma=1.0 / max(len(last["pairs"]), 1),
                role="implementation_check",
                note="★ 압축이 결합을 1차극소로 부순 비율. 0 이어야 한다",
                derivation="설계가 eps_max=0.399% < eps_crit=0.703% 를 만족하므로 어떤 "
                           "아핀 압축 스텝도 결합쌍을 장벽 안쪽으로 밀지 못한다. 이 문턱은 "
                           "**쌍 하나로 실측**했다 (0.40% 유지 / 0.80% 붕괴, h*=−0.042). "
                           "조합으로 확장되는 근거: 아핀 스케일은 각 쌍에 독립으로 "
                           "같은 상대변형을 걸고, 붕괴 판정은 그 쌍의 h 만으로 결정된다. "
                           "⚠️ 다만 침투한 망이 강직해지면 이완이 방해될 수 있어 "
                           "0 이 아니면 **버그가 아니라 그 가정의 파괴**로 읽고 ε 를 줄인다."),
            MET.observable(
                "rdf_first_peak", peak, predicted=float(1.0 + P["h_min_star"]),
                role="implementation_check", tol_pct=3.0,
                note="g(r) 첫 피크 = DLVO 2차극소 결합거리여야 한다",
                derivation="페어 표를 U_star 에서 직접 만들었으므로 두 입자만의 자유에너지 "
                           "최소는 정의상 h_min 이다. **희박 극한**(φ=0.1, 배위수 z≲3 이라 "
                           "한 결합에 걸리는 다체 압축이 작다)에서 g(r) 최빈값이 그 자리에 "
                           "온다. 다체 효과로 3% 이내 이동은 허용한다 — 그래서 tol_pct=3."),
            MET.observable("coordination_number", last["z"], role="measurement",
                           note="배위수 z — 그림(도식)의 트리는 1.905 였다"),
            MET.observable("independent_loops", last["loops"], role="measurement",
                           note="독립 고리 수 (Betti-1). 그림은 0(트리)"),
            MET.observable("dangling_fraction", last["dangling"], role="measurement",
                           note="자유단 비율. 그림은 6/21 = 0.286"),
            MET.observable("largest_cluster_fraction", last["largest_cluster"],
                           role="measurement"),
            MET.observable("percolation", last["percolation"], role="measurement",
                           note="침투한 축 수 / 3"),
            MET.observable("fractal_dimension", last["d_f"], role="measurement",
                           note="질량-반지름 스케일링. DLCA 는 ~1.8, 조밀 그물은 →3"),
            MET.observable("min_separation", last["min_sep"], role="measurement",
                           note="최소 중심간 거리 [d]. <1 이면 겹침"),
            MET.observable("pe_per_particle", last["pe_energy"], role="measurement",
                           note="⟨U⟩/N [kT]. z·well_depth/2 와 일관되는지 교차확인용"),
            MET.observable("r_nn_mean", last["r_nn"], role="measurement",
                           note="기준 64개의 최근접 이웃거리 평균 [d] — FROZEN 판정 지표"),
        ]
        arrays = {f"series_{k}": np.asarray(v) for k, v in series.items()}
        arrays["rdf_r"], arrays["rdf_g"] = mid, g
        arrays["rdf_fine_r"], arrays["rdf_fine_g"] = fine_r, fine_g
        arrays["final_positions"] = pos
        arrays["final_pairs"] = (last["pairs"] if len(last["pairs"])
                                else np.zeros((0, 2), dtype=int))
        # ★ `finalize` 의 `extra` 는 `run.execute` 가 `extra={"result": extra}` 로 감싸
        #   `metrics["result"]` 에 **그대로** 넣는다 (run.py:393). 즉 여기서 다시
        #   `{"result": …}` 로 싸면 `metrics["result"]["result"]` 로 한 겹 깊어진다 —
        #   실측으로 확인했다. **평평하게** 넣는다.
        return {"observables": obs, "arrays": arrays,
                "extra": {
                    "n_particles": n, "L_final": L, "r_bond_star": r_bond,
                    "n_links": int(len(last["pairs"])),
                    "eps_max": float(P["eps_max"]), "eps_crit": float(P["eps_crit"]),
                    "frozen_indicator": "r_nn_mean [d] — ⟨U⟩/N 은 응집 초기에 정확히 0 이라 "
                                        "FROZEN 오탐이 나고, 그러면 finalize 가 건너뛰어져 "
                                        "관측량이 통째로 사라진다 (abp-rod 전례)",
                }}

    if n_stages > 0:
        phases = [RUN.Phase("응집", n_agg, expect_steady=False, note="φ₀ 고정, 클러스터 형성"),
                  RUN.Phase("압축", n_stage * n_stages, expect_steady=False,
                            note=f"{n_stages} 단계 등방 압축"),
                  RUN.Phase("후이완", n_post, expect_steady=True, note="φ₁ 에서 구조 정착")]
    else:
        phases = [RUN.Phase("relaxation", n_post, expect_steady=True,
                            note="돋아난 망을 φ₁ 에서 이완 (압축 없음)")]
    return RUN.Build(
        sim=sim, forces=forces, n_particles=n, sample=sample,
        pe_per_particle=r_nn_mean, n_eq=0,       # ★ ⟨r_nn⟩ — 위 주석 참조
        n_prod=n_agg + n_stage * n_stages + n_post,
        sample_every=int(Nm["sample_every"]), phases=phases,
        gsd_path=(Path(outdir) / "traj_A.gsd") if outdir else None,
        tags=["3d", "dlvo", "gelation", "stage1",
              "compression" if n_stages > 0 else "sprout"],
        physical={"phi0": float(P["phi0"]), "phi1": float(P["phi1"]),
                  "n_stages": n_stages, "r_bond_star": r_bond},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="리포트만")
    ap.add_argument("--spec", action="store_true", help="스펙 저장, 실행 안 함")
    ap.add_argument("--n", type=int, default=None, help="입자 수 (기본: system.yaml 첫 값)")
    ap.add_argument("--stage-tau", type=float, default=1e-3,
                    help="압축 한 단계의 물리 시간 [τ_B]. τ_esc(0.115)에 가까우면 준정적")
    ap.add_argument("--agg-tau", type=float, default=0.2, help="압축 전 응집 시간 [τ_B]")
    ap.add_argument("--post-tau", type=float, default=0.2, help="압축 후 이완 시간 [τ_B]")
    ap.add_argument("--init", choices=("scatter", "sprout"), default="scatter",
                    help="scatter=흩뿌려 응집→압축 (DLVO 동역학이 위상을 만든다) · "
                         "sprout=씨앗에서 돋아나게 해 목표 φ 에서 직접 (위상을 부과한다)")
    ap.add_argument("--max-coord", type=int, default=4,
                    help="sprout 에서 부모로 뽑히는 최대 배위수 (작으면 성기다)")
    ap.add_argument("--n-seeds", type=int, default=8,
                    help="sprout 씨앗 개수. 흩뿌리면 성장 전선이 경계를 가로질러 만나 "
                         "감는 고리(=침투)가 생긴다. 실측 최적 8 (N=512)")
    ap.add_argument("--loop-bias", type=int, default=100,
                    help="후보 방향 시행 횟수. 결합이 가장 많이 생기는 자리를 고른다. "
                         "1=편향 없음(트리에 가까움). 실측 100 에서 고리 15→190, z 2.05→2.74")
    ap.add_argument("--dt-scale", type=float, default=1.0)
    ap.add_argument("--samples", type=int, default=200)
    ap.add_argument("--seed", type=int, default=1, help="★ <65536 (bd-hoomd 함정 12)")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/network/system.yaml")
    n = args.n or sys_["n_list"][0]
    return emit(sys_, n, args)


if __name__ == "__main__":
    raise SystemExit(main())
