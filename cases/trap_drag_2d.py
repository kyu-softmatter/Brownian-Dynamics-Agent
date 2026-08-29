"""`trap-drag-2d-hex300` — L3 무차원화 (실행은 아직 없음).

육방 콜로이드 격자(소프트 반발 A/r³ + WCA)를 광집게 트랩 하나로 등속 관통시키는 계.
능동 마이크로레올로지 / 끌린 탐침.

**1-A + 1-B의 합집합입니다** — 트랩(1-A)과 소프트 페어 격자(1-B)가 한 계에 같이 있고,
그래서 **강성이 두 개**입니다. 이 케이스가 새로 가져오는 것은 세 가지입니다:

  ① `dt` 를 **두 강성 중 빠른 쪽**이 정한다 — τ_k = γ/k_t 가 τ_int = γ/U''(r_min) 보다
     218배 빠릅니다. 1-A는 트랩만, 1-B는 페어만 있어서 경쟁이 없었습니다.
  ② **이류 시간척도** τ_v = d/v_x 가 처음 등장합니다 (움직이는 경계조건).
  ③ ★ **측정 가능성이 하드 게이트가 아니라 통계 문제로 나타납니다** — 끌림 지연
     Δr_ss = γv/k_t = 2.0 nm 가 트랩 안 열요동 ℓ_k = 20.4 nm 의 1/10 입니다 (SNR = 0.0985).
     분리 검사는 전부 통과하는데 **한 번 횡단으로는 원하는 정밀도가 안 나옵니다.**
     L3가 잡아야 하는 종류의 결함이고, 실행 전에 보입니다.

L4(실행)는 `bdbot.run` 이 돌리고 **수치 건전성만** 판정합니다 (발산·NaN·frozen·표류).
물리 대조는 `metrics.observables` 의 role 체계입니다 — L4가 하지 않습니다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/trap_drag_2d.py --report        # L3 리포트만
    $PY cases/trap_drag_2d.py --spec          # L3 스펙
    $PY cases/trap_drag_2d.py --smoke         # L4 배선 확인 (짧게)
    $PY cases/trap_drag_2d.py                 # L4 본 실행 (~6.8e6 스텝)
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
from bdbot import nondim as ND, run as RUN, scales as SC, sim as SIM  # noqa: E402
from bdbot import stats as ST, traps as TR  # noqa: E402
from bdbot.pairpot import HEX_NN, U2_star, U_star, a_mean_star, approach_distance  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
R_TABLE_MIN = 0.5          # pair.Table 하한 (1-B와 동일). 함정 11 방어의 기준
N_CYCLE_TARGET = 100.0     # 통계 검사 기준 — 격자 주기를 몇 번 지나야 하는가
STAT_TARGET_PCT = 2.0      # ⟨F_drag⟩ 목표 정밀도 (system.yaml numerics.stat_target_pct)
# ★ GSD 저장 주기 [스텝]. None 이면 전체의 1/300. **스펙에 넣지 않습니다** — 출력 설정이지
#   물리가 아니라서 해시에 들어가면 같은 계가 다른 run_id 를 갖게 됩니다.
#   빠른 v 는 끌기 구간이 짧아(v=32 에서 12,544 스텝) 기본 주기로는 프레임이 2개뿐입니다.
GSD_EVERY: int | None = None


# ════════════════════════════════════════════════════════════════════════
# ① 물리계 (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path, v_sel: float | None = None) -> dict:
    """`v_sel` 은 `external.drag_velocity` 리스트에서 고를 값 [µm/s]. None 이면 스케치 값."""
    raw = yaml.safe_load(path.read_text())
    P = load_node
    # ★ 속도 스윕 — 리스트에서 한 점을 고른다. 지어내지 않고 **L2에 있는 값만** 쓴다.
    vnode = dict(raw["external"]["drag_velocity"])
    vlist = vnode["value"] if isinstance(vnode["value"], list) else [vnode["value"]]
    vlist = [float(x) for x in vlist]
    vsel = float(v_sel) if v_sel is not None else (0.5 if 0.5 in vlist else vlist[0])
    if not any(abs(vsel - x) < 1e-12 for x in vlist):
        raise ValueError(f"v={vsel} 는 L2 의 drag_velocity 목록에 없습니다: {vlist}. "
                         "스펙은 물리계에서 유도된 것만 씁니다 (규칙 2) — "
                         "새 속도를 쓰려면 system.yaml 에 먼저 추가하세요.")
    vnode["value"] = vsel
    r3, wca = raw["interactions"][0], raw["interactions"][1]
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "phi": float(raw["geometry"]["area_fraction"]["value"]),
        # ★ 정합 육방의 셀 수. N = n_x·n_y 이므로 N을 따로 읽고 대조한다 (아래 build_ledger).
        "n_x": int(raw["geometry"]["n_x"]["value"]),
        "n_y": int(raw["geometry"]["n_y"]["value"]),
        # L2가 **적어둔** 박스. 아래에서 격자로부터 유도한 값과 대조합니다 — 유도값을
        # 그냥 믿으면 YAML 이 조용히 어긋나도 아무도 모릅니다 (physical.verify 와 같은 태도).
        "box_x": load_node(raw["geometry"]["box_length_x"]),
        "box_y": load_node(raw["geometry"]["box_length_y"]),
        "A": float(r3["amplitude_A"]["value"]),
        "r_c": P(r3["cutoff"]),
        "wca_eps_kT": float(wca["epsilon"]["value"]),
        "k_t": P(raw["external"]["stiffness"]),
        "v_x": P(vnode), "v_list": vlist,
        "n_trapped": int(raw["external"]["n_trapped"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② 스케일 원장 (bd-physics §0 ①②)
#    1-A(길이 3·시간 5)와 1-B(길이 5·시간 5)를 합쳐 길이 7 · 시간 7 이 됩니다.
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, *, dt_scale=1.0, n_traverse=1.0,
                 warm=10.0, equil=20.0, relax=40.0) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B = b["kT"], b["gamma"], b["tau_B"]
    k_t = sys_["k_t"].value.to("N/m")
    v_x = sys_["v_x"].value.to("m/s")
    A, eps, phi, N = sys_["A"], sys_["wca_eps_kT"], sys_["phi"], sys_["N"]

    # ── ★ 정합 육방 격자 (정사각 박스가 아니다) ────────────────────────
    #   L_x = n_x·a_NN , L_y = n_y·(√3/2)·a_NN , N = n_x·n_y , n_y 짝수.
    #   예전에는 L = a_mean√N (정사각)이었는데 그러면 초기 격자가 이음매에서 안 맞고,
    #   이 계의 관측량(격자 변형장·ψ₆)이 바로 그 결함에 민감하다.
    n_x, n_y = sys_["n_x"], sys_["n_y"]
    if n_x * n_y != N:
        raise ValueError(f"정합 육방이 깨졌습니다: n_x·n_y = {n_x}·{n_y} = {n_x*n_y} ≠ N = {N}")
    if n_y % 2:
        raise ValueError(f"n_y = {n_y} 가 홀수입니다 — 엇갈린 행의 주기가 2행이라 "
                         "주기경계에서 이어지지 않습니다")
    a_star = a_mean_star(phi)                       # a_mean/d
    a_nn_star = HEX_NN * a_star                     # a_NN/d — 격자 상수는 이쪽이다
    Lx_star = n_x * a_nn_star
    Ly_star = n_y * (math.sqrt(3) / 2) * a_nn_star
    # φ 는 이 구성에서 **항등적으로** 보존된다 (φ = πd²/(4a_mean²) 이고 L_x·L_y = N a_mean²).
    #   그래서 φ 를 여기서 재계산해 대조해봐야 항상 통과한다 — 돌지 않는 검사다.
    #   실제로 의미 있는 대조는 **L2가 적어둔 박스**와의 대조다: L_x·L_y 를 φ·n_x·n_y 에서
    #   유도하므로, YAML 의 box_length_* 가 조용히 어긋나면 여기서만 잡을 수 있다.
    for sym, got_star, node in (("L_x", Lx_star, sys_["box_x"]), ("L_y", Ly_star, sys_["box_y"])):
        want = float((node.value.to("m") / d).to("dimensionless").magnitude)
        if abs(got_star - want) > 1e-3 * want:
            raise ValueError(
                f"{sym} 가 L2 선언과 어긋납니다: 격자 유도 {got_star*float(d.to('um').magnitude):.4f} µm "
                f"vs system.yaml {node.value.to('um'):~.4gP} "
                f"({100*(got_star-want)/want:+.3f}%). φ·n_x·n_y 와 박스를 함께 고치세요.")
    r_c_star = float((sys_["r_c"].value.to("m") / d).to("dimensionless").magnitude)
    r_min_star, crit, u_rms_rel, state = approach_distance(A, a_star, eps)

    # ★ 두 개의 강성 → 두 개의 이완시간. 둘 다 γ/(국소 강성) 이다 (bdbot.checks).
    tau_k = C.relaxation_time(gamma, k_t)                                  # 트랩
    tau_int = (tau_B / float(U2_star(r_min_star, A, eps))).to("s")          # 페어 (kT/d² 단위)
    tau_v = (d / v_x).to("s")                                              # 이류 (새로 등장)

    # dt 는 **가장 빠른** 척도가 정한다 — 여기서는 트랩이다 (페어보다 218배 빠름).
    tau_dt = min((tau_k, tau_int), key=lambda q: q.to("s").magnitude)
    dt = dt_scale * C.dt_from_gate(tau_dt)
    # ★ 횡단은 **끌기 방향**인 L_x 다 (짧은 변 L_y 가 아니다).
    #   `tau_cross` 는 프로토콜과 **무관한** 계의 성질입니다 — 예전에는 이걸 따로 두지 않고
    #   기하 검사(항적 치유)가 `T_obs` 를 분모로 썼는데, `--traverse` 를 바꾸면 검사의
    #   의미가 조용히 달라졌습니다 (traverse=0.001 로 평형화 연구를 돌리려다 하드 FAIL).
    #   계의 성질을 재는 검사는 계의 척도만 봐야 합니다.
    tau_cross = (Lx_star * d / v_x).to("s")

    # 네 구간 (단위: τ_int, 끌기만 횡단시간). T_obs 는 **전체 프로토콜**이다.
    n_of = lambda x: int(round(x * float((tau_int / dt).to(""))))
    n_warm, n_equil, n_relax = n_of(warm), n_of(equil), n_of(relax)
    n_drag = int(round(float((n_traverse * tau_cross / dt).to(""))))
    T_obs = ((n_warm + n_equil + n_drag + n_relax) * dt).to("s")

    l_k = (kT / k_t) ** 0.5                            # 트랩 안 열요동 폭
    dr_ss = (gamma * v_x / k_t).to("m")                # 끌림의 정상상태 지연

    lg = SC.ScaleLedger()
    lg.add_length("dr_ss", dr_ss, "γv/k_t 끌림 지연 (신호)", star=True)
    lg.add_length("l_k", l_k.to("m"), "√(kT/k_t) 트랩 요동 (잡음)", star=True)
    lg.add_length("d", d, "입자 지름 (기준)")
    lg.add_length("r_min", r_min_star * d, "최근접 접근거리")
    lg.add_length("a_mean", a_star * d, "평균 간격")
    lg.add_length("a_NN", a_nn_star * d, "육방 최근접 = 격자 상수 (a_mean 이 아님)")
    lg.add_length("r_c", r_c_star * d, "컷오프")
    # ★ 직사각 박스라 변이 둘이다. 최소이미지를 정하는 것은 **짧은 변**이므로 거기에
    #   `box` 역할을 준다 (역할은 기호가 아니라 기능이다 — bdbot.scales).
    #   긴 변은 끌기 방향이고 T_obs·격자 주기 수를 정한다.
    short, long_ = (("L_y", Ly_star), ("L_x", Lx_star)) if Ly_star <= Lx_star \
        else (("L_x", Lx_star), ("L_y", Ly_star))
    lg.add_length(short[0], short[1] * d, f"박스 짧은 변 — 최소이미지의 분모 ({n_y}행)"
                  if short[0] == "L_y" else f"박스 짧은 변 — 최소이미지의 분모 ({n_x}열)",
                  role="box")
    lg.add_length(long_[0], long_[1] * d,
                  f"박스 긴 변 = 끌기 방향 ({n_x}열 × a_NN)"
                  if long_[0] == "L_x" else f"박스 긴 변 ({n_y}행)")
    lg.add_time("tau_p", b["tau_p"], "m/γ 관성 이완", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_k", tau_k, "γ/k_t 트랩 — 최속, dt를 정한다", star=True)
    lg.add_time("tau_int", tau_int, f"γ/U''(r_min={r_min_star:.3f}d) 페어")
    lg.add_time("tau_v", tau_v, "d/v_x 이류 (움직이는 트랩)")
    lg.add_time("tau_cell", (HEX_NN * a_star * d / v_x).to("s"), "a_NN/v_x 격자 주기")
    lg.add_time("tau_cross", tau_cross, "L_x/v_x 박스 횡단 (계의 성질 — 프로토콜 무관)")
    lg.add_time("tau_B", tau_B, "d²/D_t 확산 (기준)")
    lg.add_time("T_obs", T_obs, "관측창 = 네 구간 전체", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.add_energy("U_a", (float(U_star(a_star, A, eps)) * kT).to("J"),
                  "U(a_mean) 평균간격 결합 = Γ kT")
    # ★ U(d) = (A+ε) kT — A kT 가 아니다 (1-B에서 L3 검사가 잡은 라벨 오류)
    lg.add_energy("U_d", (float(U_star(1.0, A, eps)) * kT).to("J"),
                  "U(d) 접촉 결합 = (A+ε_WCA) kT")
    lg.add_energy("k_t_d2", (k_t * d**2).to("J"), "k_t d² 트랩 강성")

    lg.derived = dict(gamma=gamma, D_t=b["D_t"], m=b["m"], kT=kT, d=d, tau_B=tau_B,
                      tau_k=tau_k, tau_int=tau_int, tau_v=tau_v, dt=dt, T_obs=T_obs,
                      tau_cross=tau_cross, n_warm=n_warm, n_equil=n_equil,
                      n_drag=n_drag, n_relax=n_relax,
                      tau_int_steps=float((tau_int / dt).to("")),
                      a_star=a_star, a_nn_star=a_nn_star, n_x=n_x, n_y=n_y,
                      Lx_star=Lx_star, Ly_star=Ly_star, r_c_star=r_c_star,
                      r_min_star=r_min_star, crit=crit, u_rms_rel=u_rms_rel, state=state,
                      l_k=l_k.to("m"), dr_ss=dr_ss, k_t=k_t, v_x=v_x,
                      tau_dt_name=("tau_k" if tau_dt is tau_k else "tau_int"))
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ 이 계에는 강성이 **두 개**(트랩 k_t, 페어 U'')이고 "
        "트랩이 218배 빠르다 → dt는 트랩이 정한다. 1-A는 트랩만, 1-B는 페어만 있었다.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ③ 무차원수 + ④ 분리 검사 (bd-physics §3, §4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg):
    D = lg.derived
    f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)
    g = lg.get
    A, eps, phi = sys_["A"], sys_["wca_eps_kT"], sys_["phi"]
    a_star, r_c_star = D["a_star"], D["r_c_star"]
    Lx_star, Ly_star = D["Lx_star"], D["Ly_star"]
    L_min = min(Lx_star, Ly_star)
    Gamma = float(U_star(a_star, A, eps))
    snr = lg.ratio("lengths", "dr_ss", "l_k")

    # 한 번 횡단에서 얻는 ⟨Δr_ss⟩ 의 상대 정밀도.
    #   독립표본 n ≈ T_obs/(2τ_k) (변위 상관시간 τ_k), 표본당 잡음 ℓ_k
    #   ⟹ SEM/Δr_ss = 1/(SNR·√n)
    n_indep = 0.5 * lg.ratio("times", "tau_cross", "tau_k")   # 끌기 한 번 횡단 기준
    prec_pct = 100.0 / (snr * math.sqrt(n_indep))
    n_cell = lg.ratio("times", "tau_cross", "tau_cell")   # 한 번 횡단의 격자 주기 = n_x

    groups = [
        ND.Group("Gamma", Gamma, ("energies", "U_a"), ("energies", "kT"),
                 "A(d/a_mean)³", "페어 결합 → 육방 결정 ★"),
        ND.Group("A", A, None, None, "", "r⁻³ 진폭 (입력, ★제안)"),
        ND.Group("U(d)/kT", float(U_star(1.0, A, eps)),
                 ("energies", "U_d"), ("energies", "kT"), "(A+ε_WCA)", "접촉 결합"),
        ND.Group("k*", lg.ratio("energies", "k_t_d2", "kT"),
                 ("energies", "k_t_d2"), ("energies", "kT"), "k_t d²/kT",
                 "트랩 vs 열요동 (매우 뻣뻣)"),
        ND.Group("phi", phi, None, None, "", "밀집도"),
        ND.Group("Pe_drag", lg.ratio("times", "tau_B", "tau_v"),
                 ("times", "tau_B"), ("times", "tau_v"), "v_x d/D_t", "이류 vs 확산"),
        ND.Group("SNR", snr, ("lengths", "dr_ss"), ("lengths", "l_k"), "Δr_ss/ℓ_k",
                 "★ 끌림 신호가 열요동의 1/10"),
        ND.Group("a_mean/d", a_star, ("lengths", "a_mean"), ("lengths", "d"), "", "평균간격"),
        ND.Group("a_NN/a_mean", HEX_NN, ("lengths", "a_NN"), ("lengths", "a_mean"),
                 "√(2/√3)", "육방 최근접 (파라미터 없는 예측)"),
        ND.Group("L_x/d", Lx_star, ("lengths", "L_x"), ("lengths", "d"), "n_x·a_NN/d",
                 "박스 긴 변 = 끌기 방향"),
        ND.Group("L_y/d", Ly_star, ("lengths", "L_y"), ("lengths", "d"),
                 "n_y·(√3/2)a_NN/d", "박스 짧은 변"),
        ND.Group("L_x/L_y", Lx_star / Ly_star, ("lengths", "L_x"), ("lengths", "L_y"),
                 "", "★ 종횡비 — 정합이 정한다 (정사각이 아니다)"),
        ND.Group("L_x/a_NN", Lx_star / D["a_nn_star"], ("lengths", "L_x"), ("lengths", "a_NN"),
                 "n_x", "★ 정합 — 정수여야 한다"),
        ND.Group("L_y/a_NN", Ly_star / D["a_nn_star"], ("lengths", "L_y"), ("lengths", "a_NN"),
                 "n_y·√3/2", f"★ 정합 — 행 간격의 정수배 (n_y = {D['n_y']}, 짝수)"),
        ND.Group("r_c/a_mean", r_c_star / a_star, ("lengths", "r_c"), ("lengths", "a_mean"),
                 "", "컷오프 (이웃 껍질 수)"),
        ND.Group("tau_k/tau_int", lg.ratio("times", "tau_k", "tau_int"),
                 ("times", "tau_k"), ("times", "tau_int"), "",
                 "★ 두 강성의 비 — 트랩이 218배 빠르다"),
        ND.Group("dt/tau_k", lg.ratio("times", "dt", "tau_k"),
                 ("times", "dt"), ("times", "tau_k"), "", "적분 해상 (지배 척도)"),
        ND.Group("T_obs/tau_B", lg.ratio("times", "T_obs", "tau_B"),
                 ("times", "T_obs"), ("times", "tau_B"), "", "관측창 (네 구간 전체)"),
        ND.Group("tau_cross/tau_B", lg.ratio("times", "tau_cross", "tau_B"),
                 ("times", "tau_cross"), ("times", "tau_B"), "", "박스 횡단"),
        ND.Group("n_cell", n_cell, ("times", "tau_cross"), ("times", "tau_cell"),
                 "L_x/a_NN", "한 번 횡단의 격자 주기 수 = n_x"),
        ND.Group("St", lg.ratio("times", "tau_p", "tau_B"),
                 ("times", "tau_p"), ("times", "tau_B"), "tau_p/tau_B", "관성 vs 확산"),
    ]
    checks = [
        C.Check("model", "관성 무시     τ_p/τ_k", lg.ratio("times", "tau_p", "tau_k"),
                C.GATE, "<=",
                "τ_dyn = 이 계의 최속 관심척도 = τ_k (트랩). dt와 무관 (bd-physics §4)"),
        C.Check("integration", "트랩 해상     dt/τ_k", lg.ratio("times", "dt", "tau_k"),
                C.GATE, "<=",
                f"편향 ≈ {C.bias_from_dt(g('times', 'dt'), g('times', 'tau_k')):.3f}% "
                "(선형계 닫힌 형태 — 트랩은 선형이라 그대로 성립)"),
        C.Check("integration", "페어 해상     dt/τ_int", lg.ratio("times", "dt", "tau_int"),
                C.GATE, "<=",
                f"τ_int = γ/U''(r_min={D['r_min_star']:.3f}d), {D['crit']} 기준. "
                "트랩이 정한 dt라 여유가 크다"),
        C.Check("integration", "이류 해상     dt/τ_v", lg.ratio("times", "dt", "tau_v"),
                C.GATE, "<=",
                "한 스텝에 트랩 중심이 지름의 1% 이상 움직이면 안 된다 (움직이는 경계조건)"),
        C.Check("geometry", "컷오프       r_c/(L_min/2)", r_c_star / (L_min / 2), 1.0, "<=",
                f"최소 이미지 (bd-hoomd 함정 6). 직사각이므로 **짧은 변** "
                f"L_{'y' if Ly_star <= Lx_star else 'x'} 기준. 위반 시 과거 +1856% 사례"),
        # ★ 두 축을 함께 본다. x는 a_NN 의 정수배, y는 행 간격 (√3/2)a_NN 의 **짝수**배.
        #   홀수 행이면 엇갈림이 주기경계에서 어긋난다 (build_ledger 가 먼저 raise 한다).
        C.Check("geometry", "정합 육방    max|정수 이탈|",
                max(abs(Lx_star / D["a_nn_star"] - D["n_x"]),
                    abs(Ly_star / (math.sqrt(3) / 2 * D["a_nn_star"]) - D["n_y"]),
                    float(D["n_y"] % 2)),
                1e-9, "<=",
                f"★ 격자가 주기박스와 정합인가. L_x = {D['n_x']}·a_NN, "
                f"L_y = {D['n_y']}·(√3/2)a_NN, n_y 짝수. 비정합이면 이음매에 결함이 "
                f"주입되고 관측량(격자 변형장·ψ₆)이 바로 거기에 민감하다"),
        C.Check("geometry", "코어 여유    r_table_min/r_min", R_TABLE_MIN / D["r_min_star"],
                1.0, "<=",
                f"pair.Table 함정 11: r<{R_TABLE_MIN}d 면 힘이 0"),
        C.Check("geometry", "항적 치유    v τ_int/L_x", lg.ratio("times", "tau_int", "tau_cross"),
                1.0, "<=",
                "★ 주기박스라 탐침이 자기 항적으로 되돌아온다. 격자가 치유되는 거리 "
                "v·τ_int 가 끌기 방향 박스변보다 짧아야 한다"),
        C.Check("statistics", "관측창       T_obs/τ_k", lg.ratio("times", "T_obs", "tau_k"),
                100.0, ">=", "트랩 상관시간 기준", hard=False),
        C.Check("statistics", "측정가능성   SNR", snr, 1.0, ">=",
                "★ SNR<1 — 끌림 신호가 열요동에 묻힌다. 1회 표본으로는 안 보이고 "
                "평균화로만 볼 수 있다 (아래 정밀도 검사)", hard=False),
        C.Check("statistics", "끌림힘 정밀도 SEM/Δr_ss [%]", prec_pct, STAT_TARGET_PCT, "<=",
                f"★ 한 번 횡단으로 {prec_pct:.2f}% (목표 {STAT_TARGET_PCT:g}%). "
                f"독립표본 T_obs/2τ_k = {n_indep:,.0f}개. "
                f"{math.ceil((prec_pct/STAT_TARGET_PCT)**2):.0f}회 횡단이면 목표에 닿는다 — "
                f"또는 격자 {sys_['N']}개의 변형장을 관측량으로 쓴다", hard=False),
        C.Check("statistics", "격자 주기 수  L_x/a_NN", n_cell, N_CYCLE_TARGET, ">=",
                "★ 끌림힘이 격자 주기로 변조되면(stick-slip) 주기 평균이 필요하다. "
                "한 번 횡단은 √N 방향으로만 늘어난다", hard=False),
    ]
    return groups, checks, Gamma, dict(prec_pct=prec_pct, n_indep=n_indep, n_cell=n_cell,
                                       snr=snr)


def report_blocks(sys_, lg, extra, n_warm, n_equil, n_prod, n_relax,
                  tau_int_steps, args):
    D = lg.derived
    inp = [R.kv(k, f"{sys_[k].value:~.4gP}", sys_[k].tier, sys_[k].source[:44])
           for k in ("d", "T", "eta", "rho_p", "k_t", "v_x")]
    inp += [
        R.kv("N", f"{sys_['N']}", 0,
             f"sketch 'N ~ 300' → 정합 육방 {D['n_x']}×{D['n_y']} (Δ{100*(sys_['N']-300)/300:+.1f}%)"),
        R.kv("A", f"{sys_['A']:g}", 3, "★제안 — 1-B에서 육방 결정 확인된 값"),
        R.kv("phi", f"{sys_['phi']:g}", 3, "★제안 — 1-B와 동일 (스케치 미기재)"),
    ]
    der = [f"  {k:<8} = {D[k].to_compact():~.4gP}" for k in ("gamma", "D_t", "m")]
    der += [
        f"  상태 추정 = {D['state']}  (Lindemann σ_bond/a_NN = {D['u_rms_rel']:.4f}, 기준 0.15)",
        f"  Δr_ss = γv/k_t = {D['dr_ss'].to('nm'):~.3fP}   vs   ℓ_k = {D['l_k'].to('nm'):~.2fP}"
        f"   → SNR = {extra['snr']:.4f}",
        f"  ★ Δr_ss 는 탐침의 **맨 Stokes 지연**이다. 격자가 더하는 힘이 이 위에 얹히는데,",
        f"    그 크기는 예측이 없다 (measurement) — 이 계를 계산하는 이유가 그것이다.",
        f"  ★ 트랩이 페어보다 압도적으로 뻣뻣하다 (k* = {lg.ratio('energies','k_t_d2','kT'):.3g}"
        f"  vs  Γ = {float(U_star(D['a_star'], sys_['A'], sys_['wca_eps_kT'])):.2f})",
        f"    → 탐침은 트랩에 거의 고정된 '단단한' 구동 탐침이다 (일정 속도 경계조건).",
    ]
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}"
        f"  = {lg.ratio('times', 'dt', 'tau_B'):.3e} τ_B"
        f"  = 10⁻²·{D['tau_dt_name']}",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}"
        f"  = {lg.ratio('times', 'T_obs', 'tau_B'):.2f} τ_B  (박스 횡단 L_x/v_x)",
        f"  구간    = 예열 {n_warm:,} + 평형 {n_equil:,} + 끌기 {n_prod:,}"
        f" + 이완 {n_relax:,}  =  {n_warm+n_equil+n_prod+n_relax:,} 스텝  × N={sys_['N']}",
        f"            (예열·평형·이완은 τ_int = {tau_int_steps:,.0f} 스텝의"
        f" {args.warm:g}·{args.equil:g}·{args.relax:g}배)",
        f"  격자    = 정합 육방 {D['n_x']}×{D['n_y']} = {sys_['N']}"
        f"   L_x×L_y = {D['Lx_star']*float(D['d'].to('um').magnitude):.2f}"
        f" × {D['Ly_star']*float(D['d'].to('um').magnitude):.2f} µm"
        f"  (종횡비 {D['Lx_star']/D['Ly_star']:.4f})",
        f"  ⚠ 한 번 횡단의 ⟨F_drag⟩ 정밀도는 {extra['prec_pct']:.2f}% 이고 격자 주기는"
        f" {extra['n_cell']:.0f}회뿐이다.",
        f"    격자 변형장(입자 {sys_['N']}개)을 1차 관측량으로 쓰는 것이 설계 의도다"
        f" (observation.yaml B2).",
        "  ⚠ 정합이 종횡비를 정한다 — 정사각이 아니다. 유한크기 인공효과가 x·y로 비대칭일 수"
        " 있고, 그게 격자 변형장에 보이는지는 미확인이다.",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# ⑤ L4 — 스펙만 읽고 계를 세운다 (bdbot.run 이 돌리고 판정한다)
# ════════════════════════════════════════════════════════════════════════
def hex_lattice(n_x: int, n_y: int, a_nn: float) -> np.ndarray:
    """정합 육방 격자 (n_x·n_y 개). 행은 x 방향, 홀수 행이 a_NN/2 만큼 엇갈린다.

    박스는 `[-L_x/2, L_x/2) × [-L_y/2, L_y/2)`, `L_x = n_x a_NN`,
    `L_y = n_y (√3/2) a_NN`. n_y 가 짝수라야 주기경계에서 엇갈림이 이어진다.
    """
    row = math.sqrt(3) / 2 * a_nn
    Lx, Ly = n_x * a_nn, n_y * row
    j, i = np.divmod(np.arange(n_x * n_y), n_x)
    x = (i + 0.5 * (j % 2)) * a_nn - Lx / 2
    y = j * row - Ly / 2
    return np.c_[x, y]


PH_WARM, PH_EQ, PH_DRAG, PH_RELAX = "예열", "평형", "끌기", "이완"
TAIL_FRAC = 0.2            # 이완 '끝값'으로 볼 뒷부분 비율


def fit_relaxation(t, y, y_eq):
    """이완 `A·exp(−t/τ)+C` 피팅 + **피팅에 의존하지 않는** 회복률.

    ⚠️ 처음에 `curve_fit` 를 경계 없이 돌렸다가 **직선으로 달아났습니다**:
       τ = 6.9e4 τ_int, C = −404 kT (물리적으로 불가능), τ 오차가 값의 400배.
       τ ≫ 관측창이면 `A e^{−t/τ} ≈ A(1 − t/τ)` 라 노이즈 있는 감쇠에 직선이 거의 똑같이
       잘 맞습니다. 그래서 경계를 겁니다:
         C  는 데이터 범위 안        (에너지가 음수로 갈 수 없다)
         τ  는 관측창의 3배 이내      (그보다 길면 애초에 이 창으로 못 잰다)
       그리고 τ 가 상한에 붙으면 "창 안에서 이완이 안 끝났다"로 **보고**합니다.

    회복률은 피팅과 독립입니다 — 끝 20% 평균이 평형선에 얼마나 돌아왔는가.
    피팅이 실패해도 이 숫자는 답을 줍니다.
    """
    t = np.asarray(t, float); y = np.asarray(y, float)
    out = {"n": int(t.size)}
    if t.size < 10:
        return out
    tr = t - t[0]
    win = float(tr[-1])
    n_tail = max(3, int(TAIL_FRAC * t.size))
    y0 = float(y[:n_tail].mean())            # 이완 시작 직후
    y_end = float(y[-n_tail:].mean())        # 이완 끝
    denom = y0 - y_eq
    out.update(y_start=y0, y_end=y_end, window=win,
               recovered_frac=float((y0 - y_end) / denom) if denom else float("nan"),
               residual=float(y_end - y_eq))
    try:
        from scipy.optimize import curve_fit
        f = lambda x, a, tau, c: a * np.exp(-x / tau) + c
        lo = (-10 * abs(denom) - 1, win / 200.0, float(y.min()) - 0.5)
        hi = (10 * abs(denom) + 1, 3.0 * win, float(y.max()) + 0.5)
        p0 = (denom, max(win / 5, win / 100), y_eq)
        p0 = tuple(min(max(v, l), h) for v, l, h in zip(p0, lo, hi))
        pf, pc = curve_fit(f, tr, y, p0=p0, bounds=(lo, hi), maxfev=40000)
        tau = float(pf[1])
        # ★ 판정 기준은 "창 안에 들어왔나"가 아니라 **몇 e-folding 을 봤나** 입니다.
        #   τ = 2.5×창 인데 회복률이 94% 로 나오는 모순이 실제로 있었습니다 — 단일 지수로는
        #   설명이 안 된다는 뜻이고(초반 빠름 + 느린 꼬리), τ 를 인용하면 안 됩니다.
        efold = win / tau if tau > 0 else 0.0
        out.update(amp=float(pf[0]), tau_star=tau, U_inf=float(pf[2]),
                   tau_sem=float(np.sqrt(abs(pc[1, 1]))), tau_over_window=tau / win,
                   e_foldings=efold, resolved=bool(efold >= 1.0))
    except Exception as e:
        out["error"] = str(e)[:120]
    return out


@RUN.builder("trap-drag-2d-hex300")
def build(spec, outdir=None) -> RUN.Build:
    """스펙 → HOOMD 계. **케이스 YAML 을 다시 읽지 않습니다** (L2↔L4 계약은 스펙 하나).

    ★ **네 구간 프로토콜** — 사용자 질문 3종(평형 vs 구동 에너지 · 이완 다이나믹스 ·
      g(r)·결함 변화)이 한 궤적 안에서 답해집니다:

        예열  트랩 고정, 표본 버림   완전격자 IC 에서 열평형으로
        평형  트랩 고정, 표본 수집   ⟨U⟩_eq · g(r)_eq · 결함_eq  ← 기준선
        끌기  트랩 등속 이동         ⟨U⟩_drag · g(r)_drag · 결함_drag + 상승 과도
        이완  트랩 **정지**(끝 위치에 고정)  에너지가 평형으로 돌아가는 다이나믹스

      구동을 켜고 끄는 것을 `traps.make_trap(drive=...)` 의 조각별 함수로 만듭니다 —
      `velocity=` 는 껐다 켤 수 없어서 쓰지 않습니다.
    """
    import freud
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    N, n_x, n_y = int(P["N"]), int(P["n_x"]), int(P["n_y"])
    a_nn, Lx, Ly = P["a_nn_star"], P["Lx_star"], P["Ly_star"]
    A, eps, r_c = P["A"], P["wca_eps"], P["r_c_star"]
    k_star, v_star = P["k_star"], P["v_star"]
    r_tab_min = P["r_table_min_star"]
    dt_star, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_warm, n_equil = int(Nm["n_warm"]), int(Nm["n_equil"])
    n_drag, n_relax = int(Nm["n_prod"]), int(Nm["n_relax"])

    pos0 = hex_lattice(n_x, n_y, a_nn)
    sim = SIM.make_sim(SIM.frame_2d(pos0, (Lx, Ly)), seed=seed)

    cell = md.nlist.Cell(buffer=0.4)
    # r⁻³ 꼬리 — pair.Table. ★ endpoint=False (함정 10), 컷오프에서 시프트
    rr = np.linspace(r_tab_min, r_c, 1000, endpoint=False)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_c)
    tab.params[("A", "A")] = dict(r_min=r_tab_min,
                                  U=A / rr**3 - A / r_c**3, F=3 * A / rr**4)
    wca = SIM.wca(cell, epsilon=eps, sigma=1.0)          # 함정 4·11 (코어를 별도 힘으로)

    # ★ 트랩은 탐침 하나에만 (스케치: X 표시가 원 하나). 나머지는 k=0.
    probe = (n_y // 2) * n_x + n_x // 2
    k_arr = np.zeros(N)
    k_arr[probe] = k_star

    # ── 구동 프로토콜: 정지 → 등속 → 정지 (조각별 선형) ────────────────
    T1 = (n_warm + n_equil) * dt_star            # 끌기 시작 (무차원 시간)
    T2 = T1 + n_drag * dt_star                   # 끌기 종료 → 그 자리에 고정
    _off = np.zeros((N, 3))                      # 매 스텝 할당하지 않으려고 미리 잡음

    def drive(t):
        _off[probe, 0] = v_star * (min(max(t, T1), T2) - T1)
        return _off

    trap = TR.make_trap(k_arr, pos0, (Lx, Ly), dt_star=dt_star, drive=drive)

    SIM.attach_brownian(sim, dt_star, [tab, wca, trap])
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, GSD_EVERY or
                              max(1, (n_warm + n_equil + n_drag + n_relax) // 300))

    box = freud.box.Box(Lx=Lx, Ly=Ly, is2D=True)
    hexatic = freud.order.Hexatic(k=6, weighted=True)
    voro = freud.locality.Voronoi()
    r_max = min(r_c, min(Lx, Ly) / 2 - 1e-6)
    # ★ g(r) 은 구간마다 **따로** 쌓는다. 하나에 몰아 쌓으면 평형과 구동이 섞인다.
    rdf = {ph: freud.density.RDF(bins=300, r_max=r_max, r_min=0.5)
           for ph in (PH_EQ, PH_DRAG, PH_RELAX)}
    # ★ freud 는 compute() 전에 `.rdf` 를 읽으면 AttributeError 를 냅니다 (빈 객체가
    #   아니라 예외). 길이 0 인 구간(`--relax 0`)이 있으면 finalize 가 죽습니다 —
    #   실제로 평형화 연구를 돌리다 런이 끝난 뒤에 죽어서 4.5분을 날렸습니다.
    rdf_used: set = set()
    anchor3 = np.c_[pos0, np.zeros(N)]

    def xy():
        return np.array(sim.state.get_snapshot().particles.position, dtype=float)

    def pe_pair():
        """콜로이드 계의 퍼텐셜 에너지 /N. **트랩 에너지는 제외** — 그건 광집게에
        저장된 것이지 계의 것이 아니다. 트랩 쪽은 따로 잰다."""
        return float(np.array(tab.energies).sum() + np.array(wca.energies).sum()) / N

    def sample(timestep, phase):
        p = xy()
        d_probe = trap.displacement(sim.state, timestep)[probe]
        # 격자 변형장: 초기 격자 자리로부터의 변위 (탐침 제외, 최소이미지).
        # ★ 강체 병진을 뺀다 — 끌린 탐침이 격자 전체를 조금씩 끌고 가는데 그건 변형이
        #   아니다. 빼지 않으면 dev_rms 가 런이 길수록 계속 자란다 (변형이 아니라 표류).
        dev = SIM.minimum_image(p - anchor3, (Lx, Ly))[:, :2]
        dev = np.delete(dev, probe, axis=0)
        dev = dev - dev.mean(axis=0)
        vn = voro.compute((box, p)).nlist
        z = np.asarray(vn.neighbor_counts)               # Voronoi 배위수
        if phase in rdf:
            rdf[phase].compute((box, p), reset=False)
            rdf_used.add(phase)
        u_pair = pe_pair()
        u_trap = float(np.array(trap.energies).sum())
        return {
            "dx_probe": float(d_probe[0]), "dy_probe": float(d_probe[1]),
            "u_pair": u_pair, "u_trap": u_trap,
            "u_total": u_pair + u_trap / N,
            "psi6": float(np.abs(hexatic.compute((box, p), neighbors=vn)
                                 .particle_order).mean()),
            # ★ 결함 = Voronoi 배위수 ≠ 6. 전위는 5–7 쌍으로 나온다.
            "n_def": int((z != 6).sum()), "n5": int((z == 5).sum()),
            "n7": int((z == 7).sum()),
            "dev_rms": float(np.sqrt((dev**2).sum(axis=1).mean())),
            "dev_max": float(np.sqrt((dev**2).sum(axis=1)).max()),
            "min_sep": float(np.asarray(vn.distances).min()),
        }

    def finalize(cols):
        ph = cols["phase"]
        sel = {k: (ph == k) for k in (PH_EQ, PH_DRAG, PH_RELAX)}
        t = cols["_t_step"] * dt_star                     # 무차원 시간
        obs, extra = [], {}

        def stat(key, mask):
            v = cols[key][mask]
            return (float(v.mean()), float(ST.block_sem(v))) if v.size else (float("nan"),) * 2

        # ── ① 평형 vs 구동 에너지 ─────────────────────────────────────
        u_eq, u_eq_e = stat("u_pair", sel[PH_EQ])
        u_dr, u_dr_e = stat("u_pair", sel[PH_DRAG])
        # 끌기 구간의 **후반 절반**만 정상상태로 본다 (앞은 상승 과도).
        dmask = sel[PH_DRAG].copy()
        if dmask.sum() >= 4:
            idx = np.flatnonzero(dmask)
            dmask[idx[: len(idx) // 2]] = False
        u_ss, u_ss_e = stat("u_pair", dmask)
        dU = u_ss - u_eq
        dU_e = math.hypot(u_ss_e, u_eq_e)
        obs += [
            MET.observable("U_pair_equilibrium", u_eq, None, "kT/particle",
                           role="measurement",
                           note=f"트랩 고정 구간의 격자 에너지 (±{u_eq_e:.4g})"),
            MET.observable("U_pair_driven", u_ss, None, "kT/particle", role="measurement",
                           note=f"끌기 정상상태(후반 절반) (±{u_ss_e:.4g})"),
            MET.observable("dU_drive", dU, None, "kT/particle", role="measurement",
                           note=f"구동이 격자에 저장시킨 초과 에너지 = {dU:+.5g} "
                                f"± {dU_e:.4g} kT/입자 ({dU/u_eq*100:+.3f}%). "
                                f"전체로는 {dU*N:+.4g} kT. 예측 없음 — 이게 답이다"),
        ]

        # ── ② 이완 다이나믹스 ─────────────────────────────────────────
        rel = sel[PH_RELAX]
        fit = fit_relaxation(t[rel], cols["u_pair"][rel], u_eq) if rel.sum() else {}
        # ★ 결함도 같이 이완한다 — 에너지와 **다른 시간척도**일 수 있어 따로 잰다.
        fit_def = (fit_relaxation(t[rel], cols["n_def"][rel].astype(float), 0.0)
                   if rel.sum() else {})
        extra["relax_fit"] = fit
        extra["relax_fit_defects"] = fit_def
        tau_int = spec.reduced("times", "tau_int")
        if fit:
            tau_rel = fit.get("tau_star", float("nan"))
            ok_fit = fit.get("resolved", False)
            obs.append(MET.observable(
                "U_recovered_frac", fit.get("recovered_frac", float("nan")), None, "1",
                role="measurement",
                note=f"이완 끝 20% 평균이 평형선까지 돌아온 비율. "
                     f"잔차 {fit.get('residual', float('nan')):+.4g} kT/입자 "
                     f"(평형 요동 ±{u_eq_e:.4g}). ★ 피팅과 독립인 숫자다"))
            obs.append(MET.observable(
                "tau_relax", tau_rel, None, "tau_B", role="measurement",
                note=(f"⟨U⟩ 이완시간 = {tau_rel/tau_int:.2f} τ_int "
                      f"(창의 {fit.get('tau_over_window', float('nan')):.2f}배). "
                      if ok_fit else
                      f"★ 창 안에서 이완이 끝나지 않아 τ 를 정할 수 없다 "
                      f"(τ={tau_rel/tau_int:.3g} τ_int 는 상한에 붙은 값이다). ")
                     + f"국소 케이지 τ_int 는 집단 이완과 일치할 이유가 없다 — 눈금일 뿐"))

        # ── ③ g(r) · 결함 ─────────────────────────────────────────────
        for key, unit in (("n_def", "개"), ("n5", "개"), ("n7", "개"), ("psi6", "1")):
            e_m, _ = stat(key, sel[PH_EQ])
            d_m, _ = stat(key, dmask)
            r_m, _ = stat(key, rel)
            extra[key] = {"equilibrium": e_m, "driven": d_m, "relaxed": r_m}
            obs.append(MET.observable(
                f"{key}_driven_vs_eq", d_m - e_m, None, unit, role="measurement",
                note=f"평형 {e_m:.4g} → 구동 {d_m:.4g} → 이완후 {r_m:.4g} [{unit}]"))

        gr, arrays = {}, {}
        for name in rdf_used:                    # ★ 실제로 compute 된 것만 (위 주석)
            gr[name] = np.asarray(rdf[name].rdf)
        if gr:
            arrays["gr_r"] = np.asarray(rdf[next(iter(rdf_used))].bin_centers)
            key = {PH_EQ: "gr_eq", PH_DRAG: "gr_drag", PH_RELAX: "gr_relax"}
            arrays.update({key[k]: v for k, v in gr.items()})
            # 첫 봉우리 — 구조가 바뀌었는지 한 숫자로 보는 눈금
            for k, v in gr.items():
                i = int(np.argmax(v))
                extra.setdefault("g_r_peak", {})[k] = {
                    "r": float(arrays["gr_r"][i]), "g": float(v[i])}

        # ── 끌림힘 (전 세션에서 이미 있던 것) ──────────────────────────
        dx = cols["dx_probe"][dmask]
        f_drag = -k_star * float(dx.mean()) if dx.size else float("nan")
        sem = k_star * float(ST.block_sem(dx)) if dx.size else float("nan")
        obs.append(MET.observable(
            "F_drag_total", f_drag, None, "kT/d", role="measurement",
            note=f"⟨F_drag⟩ = −k*·⟨Δx⟩ (정상상태). 맨 Stokes γv = {v_star:.4g} 대비 "
                 f"초과 {100*(f_drag/v_star - 1):+.2f}% (±{100*sem/v_star:.2f}%)"))
        extra.update(f_drag_kT_per_d=f_drag, f_drag_sem=sem, f_stokes_bare=v_star,
                     excess_pct=100 * (f_drag / v_star - 1), probe_index=probe,
                     U_pair_eq=u_eq, U_pair_driven=u_ss, dU_drive=dU, dU_sem=dU_e,
                     dU_total_kT=dU * N, tau_relax_star=tau_rel,
                     min_sep_d=float(cols["min_sep"].min()),
                     dev_max_d=float(cols["dev_max"].max()))
        return {"observables": obs, "extra": extra, "arrays": arrays}

    every_dr = max(1, n_drag // 2000)
    phases = [
        RUN.Phase(PH_WARM, n_warm, collect=False,
                  note="트랩 고정 · 표본 버림 (완전격자 IC → 열평형)"),
        RUN.Phase(PH_EQ, n_equil, max(1, n_equil // 200),
                  note="트랩 고정 · 기준선 ⟨U⟩_eq · g(r) · 결함"),
        RUN.Phase(PH_DRAG, n_drag, every_dr,
                  note=f"트랩 등속 이동 v*={v_star:.4g} · 격자 {n_x}주기 횡단"),
        # ★ 이완은 에너지가 **변하는 것이 물리**다 → 표류 검사를 끈다 (거짓 경고 방지)
        RUN.Phase(PH_RELAX, n_relax, max(1, n_relax // 400), expect_steady=False,
                  note="트랩 정지(끝 위치 고정) · 에너지가 평형으로 돌아가는 다이나믹스"),
    ]

    return RUN.Build(
        sim=sim, forces=[tab, wca, trap], n_particles=N,
        sample=sample, pe_per_particle=pe_pair,
        sample_every=every_dr, phases=phases,
        tags=["2D", "soft_repulsion", "r^-3", "WCA_core", "hex_lattice", "moving_trap",
              "microrheology", "newtonian", "driven_relaxation"],
        physical={"N": N, "n_x": n_x, "n_y": n_y, "phi": P["phi"], "A": A,
                  "Gamma": P["Gamma"], "k_star": k_star, "v_star": v_star,
                  **{k: v for k, v in spec.raw["back_transform"].items()
                     if isinstance(v, float)}},
        finalize=finalize)


def make_plots(outdir: Path, spec) -> Path:
    """4패널 — 사용자 질문 3종에 대한 그림.

    ① ⟨U⟩/N 전체 시계열 (구간 음영)   평형 vs 구동 에너지
    ② 이완 구간 확대 + 지수 피팅       이완 다이나믹스
    ③ g(r) 세 구간 겹쳐 그리기         구조 변화
    ④ 결함 수 (전체·5배위·7배위)       전위 생성/소멸
    """
    import json as _json

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    # ★ 로그 축 눈금(10^{-1})은 mathtext 라 한글 폰트에 U+2212 글리프가 없어 깨진다.
    matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
    from matplotlib.ticker import FuncFormatter, NullFormatter

    z = np.load(outdir / "observables.npz", allow_pickle=False)
    m = _json.loads((outdir / "metrics.json").read_text())
    res = m.get("result", {})
    dt_star = float(spec.numerics["dt_star"])
    t = z["_t_step"] * dt_star
    ph = z["phase"]
    tau_int = spec.reduced("times", "tau_int")
    col = {PH_EQ: "tab:green", PH_DRAG: "tab:red", PH_RELAX: "tab:blue"}

    fig, ax = plt.subplots(3, 2, figsize=(13, 13.5))

    # ① 에너지 전체
    a = ax[0, 0]
    for k, c in col.items():
        mk = ph == k
        if mk.any():
            a.plot(t[mk] / tau_int, z["u_pair"][mk], ".", ms=2, color=c, label=k)
    for key, c, lab in (("U_pair_eq", "green", "⟨U⟩ 평형"),
                        ("U_pair_driven", "red", "⟨U⟩ 구동 정상상태")):
        if key in res:
            a.axhline(res[key], color=c, ls="--", lw=1.2, alpha=.8, label=lab)
    a.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$\langle U\rangle/N$  [kT]",
          title=f"① 격자 에너지 — ΔU = {res.get('dU_drive', float('nan')):+.4g} kT/입자"
                f"  (전체 {res.get('dU_total_kT', float('nan')):+.4g} kT)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ② 이완 확대 + 피팅
    # ② 이완 — **log–log**. 단일 지수가 안 맞으므로(0.4 e-folding 인데 회복 94%)
    #    거듭제곱인지 보려면 양쪽 로그가 맞다. 세로축은 U 자체가 아니라 **초과분**
    #    ΔU = U − U_eq 다 — 105 근처의 값은 로그로 그려도 아무것도 안 보인다.
    a = ax[0, 1]
    mk = ph == PH_RELAX
    u_eq = res.get("U_pair_eq", float("nan"))
    fit = {}
    if mk.any():
        tr = (t[mk] - t[mk][0]) / tau_int
        du = z["u_pair"][mk] - u_eq
        fit = fit_relaxation(t[mk], z["u_pair"][mk], u_eq)

        # 로그 구간 평균 — 생표본은 ΔU 가 음수로도 튀어 로그를 못 그린다
        pos = tr > 0
        edges = np.geomspace(max(tr[pos].min(), 1e-3), tr.max(), 26)
        ctr, val = [], []
        for lo_, hi_ in zip(edges[:-1], edges[1:]):
            s = (tr >= lo_) & (tr < hi_)
            if s.sum() >= 2:
                ctr.append(np.sqrt(lo_ * hi_)); val.append(du[s].mean())
        ctr, val = np.array(ctr), np.array(val)
        ok = val > 0
        a.loglog(ctr[ok], val[ok], "o-", ms=5, lw=1.4, color="tab:blue",
                 label=r"$\Delta U = \langle U\rangle/N - U_{eq}$ (로그구간 평균)")
        if (~ok).any():                      # 0 아래로 내려간 구간은 숨기지 않고 표시
            a.loglog(ctr[~ok], np.full((~ok).sum(), val[ok].min() * 0.5), "x",
                     ms=6, color="gray", label="ΔU ≤ 0 (평형 도달·잡음)")
        if "tau_star" in fit:
            tt = np.geomspace(ctr.min(), ctr.max(), 200)
            a.loglog(tt, abs(fit["amp"]) * np.exp(-tt * tau_int / fit["tau_star"]),
                     "-", lw=1.6, color="k",
                     label=(fr"지수 $\tau$={fit['tau_star']/tau_int:.0f}$\tau_{{int}}$ "
                            f"({fit['e_foldings']:.2f} e-folding"
                            f"{'' if fit.get('resolved') else ' — 미해상'})"))
        # 거듭제곱 눈금 — 직선이면 거듭제곱, 휘면 지수
        for p_, c_ in ((-0.5, "tab:orange"), (-1.0, "tab:red")):
            ref = val[ok][0] * (ctr[ok] / ctr[ok][0]) ** p_
            a.loglog(ctr[ok], ref, "--", lw=1, color=c_, alpha=.7, label=fr"$t^{{{p_:g}}}$ 눈금")
        if "n_def" in z:
            a2 = a.twinx()
            nd = np.array([z["n_def"][mk][(tr >= l) & (tr < h)].mean()
                           for l, h in zip(edges[:-1], edges[1:])
                           if ((tr >= l) & (tr < h)).sum() >= 2])
            m2 = nd > 0
            a2.loglog(ctr[m2], nd[m2], "s--", ms=4, lw=1.1, color="tab:purple", alpha=.75)
            a2.set_ylabel("결함 수 (z≠6)", color="tab:purple")
            a2.tick_params(axis="y", colors="tab:purple")
            a2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            a2.yaxis.set_minor_formatter(NullFormatter())
        # 선형 인셋 — 105 근처를 그대로 보고 싶을 때
        ins = a.inset_axes([0.055, 0.06, 0.40, 0.30])
        ins.plot(tr, z["u_pair"][mk], ".", ms=1.5, color="tab:blue", alpha=.45)
        w = max(5, len(tr) // 25)
        ins.plot(tr, np.convolve(z["u_pair"][mk], np.ones(w) / w, mode="same"),
                 "-", lw=1.2, color="tab:cyan")
        ins.axhline(u_eq, color="green", ls="--", lw=1.1)
        lo, hi = np.percentile(z["u_pair"][mk], [0.5, 99.5])
        ins.set_ylim(min(lo, u_eq) - 0.03, hi + 0.03)          # ★ 105 근처로 고정
        ins.set_title("선형 (kT/입자)", fontsize=7)
        ins.tick_params(labelsize=6)
    rec = f"회복 {100*fit['recovered_frac']:.0f}%" if "recovered_frac" in fit else ""
    a.set(xlabel=r"이완 시작 후 $t/\tau_{int}$", ylabel=r"$\Delta U$  [kT/입자]",
          title=f"② 이완 다이나믹스 — log–log  ({rec})")
    # ★ 로그 눈금을 10^{-1} 대신 평범한 숫자로. mathtext 의 U+2212 를 한글 폰트가
    #   못 그려서 "10□1" 로 깨졌습니다 (rcParams 로는 안 고쳐집니다 — \mathdefault 가
    #   본문 폰트를 그대로 씁니다). 값 범위가 0.02~0.4 라 십진 표기가 더 읽기 쉽습니다.
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_formatter(plain)
        axis.set_minor_formatter(NullFormatter())
    a.legend(fontsize=6.5, loc="upper right", framealpha=.92)
    a.grid(alpha=.3, which="both")

    # ③ g(r)
    a = ax[1, 0]
    if "gr_r" in z:
        for key, name in (("gr_eq", PH_EQ), ("gr_drag", PH_DRAG), ("gr_relax", PH_RELAX)):
            if key in z:
                a.plot(z["gr_r"], z[key], "-", lw=1.3, color=col[name], label=name)
        a.axvline(spec.params["a_nn_star"], ls=":", c="gray", label=r"$a_{NN}$")
        a.axhline(1, color="k", lw=.5, alpha=.5)
    a.set(xlabel="r / d", ylabel="g(r)", title="③ 동경분포함수 — 구간별",
          xlim=(0.8, min(6.0, float(z["gr_r"].max()) if "gr_r" in z else 6.0)))
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ④ 결함
    a = ax[1, 1]
    for key, c, lab in (("n_def", "k", "전체 (z≠6)"), ("n5", "tab:orange", "5배위"),
                        ("n7", "tab:purple", "7배위")):
        if key in z:
            a.plot(t / tau_int, z[key], "-", lw=.9, color=c, alpha=.85, label=lab)
    for k, c in col.items():
        mk = ph == k
        if mk.any():
            a.axvspan(t[mk].min() / tau_int, t[mk].max() / tau_int, color=c, alpha=.07)
    nmax = float(z["n_def"].max()) if "n_def" in z else 0.0
    a.set(xlabel=r"$t/\tau_{int}$", ylabel="입자 수", ylim=(0, max(nmax * 1.25, 1.0)),
          title=f"④ 결함 (Voronoi 배위수 ≠ 6) — N={spec.params['N']}")
    if nmax == 0:
        # ★ 0 은 버그가 아니라 결과다. 정합 완전격자에서 시작하고 Γ=29.7 로 강결합이면
        #   열요동이 Voronoi 위상을 못 바꾼다. 1-B(RSA 초기배치)는 결함 1.3% 가 남았다 —
        #   차이를 만드는 것은 결합세기가 아니라 **초기배치와 정합성**이다.
        a.text(0.5, 0.5, "결함 0개 — 정합 완전격자가 끝까지 유지됨\n"
                         "(1-B는 RSA 초기배치라 1.3% 잔류)",
               transform=a.transAxes, ha="center", va="center", fontsize=10,
               bbox=dict(boxstyle="round", fc="lightyellow", alpha=.9))
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ⑤ 탐침에 걸리는 힘 — F = −k*·Δ (정상상태에서 계가 탐침에 주는 항력의 크기)
    #    ★ 한 표본의 잡음이 k*·ℓ_k/d = 246 kT/d 로 **신호(~100)보다 크다** (SNR=0.0985).
    #      생점을 그대로 그리면 아무것도 안 보이므로 이동평균을 겹쳐 그린다.
    k_star, v_star = spec.params["k_star"], spec.params["v_star"]
    a = ax[2, 0]
    fx, fy = -k_star * z["dx_probe"], -k_star * z["dy_probe"]
    a.plot(t / tau_int, fx, ".", ms=1.2, color="0.75", label="$F_x$ 생표본")
    w = max(9, len(fx) // 60)
    kern = np.ones(w) / w
    a.plot(t / tau_int, np.convolve(fx, kern, mode="same"), "-", lw=1.5,
           color="tab:red", label=f"$F_x$ 이동평균 ({w}점)")
    a.plot(t / tau_int, np.convolve(fy, kern, mode="same"), "-", lw=1.1,
           color="tab:blue", alpha=.8, label=f"$F_y$ 이동평균")
    a.axhline(v_star, color="k", ls="--", lw=1.3, label=fr"맨 Stokes $\gamma v$ = {v_star:.1f}")
    a.axhline(0, color="k", lw=.5, alpha=.5)
    if "f_drag_kT_per_d" in res:
        a.axhline(res["f_drag_kT_per_d"], color="darkred", ls=":", lw=1.4,
                  label=f"⟨$F_x$⟩ 정상상태 = {res['f_drag_kT_per_d']:.0f}")
    for k, c in col.items():
        mk2 = ph == k
        if mk2.any():
            a.axvspan(t[mk2].min() / tau_int, t[mk2].max() / tau_int, color=c, alpha=.07)
    lo, hi = np.percentile(np.convolve(fx, kern, mode="same"), [1, 99])
    a.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$F = -k^*\Delta$  [kT/d]",
          ylim=(min(lo, -50) - 100, max(hi, v_star) + 150),
          title="⑤ 탐침에 걸리는 힘 (트랩이 가하는 힘 = 계가 주는 항력)")
    a.legend(fontsize=7, ncol=2); a.grid(alpha=.3)

    # ⑥ 격자 주기로 접기 — L3가 "stick-slip 이면 주기 평균이 필요하다"고 경고한 자리.
    #    끌기 구간만, 탐침이 지나온 거리를 a_NN 으로 나눈 나머지에 대해 F_x 를 평균한다.
    a = ax[2, 1]
    dmk = ph == PH_DRAG
    a_nn = spec.params["a_nn_star"]
    if dmk.sum() > 40:
        t_d = t[dmk] - t[dmk][0]
        x_travel = v_star * t_d                      # 트랩 중심이 간 거리 [d]
        phase_x = (x_travel % a_nn) / a_nn           # 0~1 (한 격자 주기)
        fxd = fx[dmk]
        nb = 12
        edges = np.linspace(0, 1, nb + 1)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        mu = np.array([fxd[(phase_x >= l) & (phase_x < h)].mean()
                       for l, h in zip(edges[:-1], edges[1:])])
        se = np.array([fxd[(phase_x >= l) & (phase_x < h)].std()
                       / max(1, np.sqrt(((phase_x >= l) & (phase_x < h)).sum()))
                       for l, h in zip(edges[:-1], edges[1:])])
        a.errorbar(np.r_[ctr, ctr + 1], np.r_[mu, mu], yerr=np.r_[se, se],
                   fmt="o-", ms=5, lw=1.4, capsize=3, color="tab:red",
                   label=f"⟨$F_x$⟩ (주기 {x_travel[-1]/a_nn:.0f}회 평균, 2주기 반복 표시)")
        a.axhline(np.nanmean(mu), color="darkred", ls=":", lw=1.3,
                  label=f"전체 평균 {np.nanmean(mu):.0f}")
        a.axhline(v_star, color="k", ls="--", lw=1.2, label="맨 Stokes")
        a.axvline(1, color="0.7", lw=.8)
        amp = (np.nanmax(mu) - np.nanmin(mu)) / 2
        a.set_title(f"⑥ 격자 주기로 접은 힘 — 변조 진폭 ±{amp:.0f} kT/d "
                    f"(평균 대비 {100*amp/max(abs(np.nanmean(mu)),1e-9):.0f}%)")
    else:
        a.set_title("⑥ 격자 주기 접기 — 표본 부족")
    a.set(xlabel="격자 주기 내 위상  (x mod $a_{NN}$)/$a_{NN}$", ylabel=r"$F_x$  [kT/d]")
    a.legend(fontsize=7); a.grid(alpha=.3)

    fig.suptitle(f"{spec.label}   Γ={spec.params['Gamma']:.2f}  "
                 f"φ={spec.params['phi']}  {spec.params['n_x']}×{spec.params['n_y']}"
                 f"={spec.params['N']}   run_id={spec.run_id}", fontsize=11)
    fig.tight_layout()
    p = outdir / "observables.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# main — L3 리포트/스펙 + L4 실행
# ════════════════════════════════════════════════════════════════════════
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="L3 리포트만")
    ap.add_argument("--spec", action="store_true", help="L3 스펙 → specs/<run_id>.json")
    ap.add_argument("--dt-scale", type=float, default=1.0, help="dt 배율 (수렴 확인용)")
    ap.add_argument("--traverse", type=float, default=1.0, help="박스 횡단 횟수 (T_obs 배율)")
    ap.add_argument("--samples", type=int, default=2000, help="끌기 구간 표본 수")
    ap.add_argument("--warm", type=float, default=10.0, help="예열 구간 [τ_int]")
    ap.add_argument("--equil", type=float, default=20.0, help="평형 측정 구간 [τ_int]")
    ap.add_argument("--relax", type=float, default=40.0, help="이완 구간 [τ_int]")
    ap.add_argument("--gsd-every", type=int, default=None,
                    help="궤적 저장 주기 [스텝]. 영상용 — 스펙/해시에 들어가지 않는다")
    ap.add_argument("--v", type=float, default=None,
                    help="끌기 속도 [µm/s]. L2 의 drag_velocity 목록 중 하나여야 한다")
    ap.add_argument("--seed", type=int, default=20260804,
                    help="시드. ★ 함정 12: HOOMD 는 16비트로 자르므로 앙상블은 "
                         "작은 연속 정수(1,2,3…)를 쓸 것 — 65536 차이면 같은 궤적이 된다")
    ap.add_argument("--force", action="store_true", help="완료된 런을 다시 실행")
    ap.add_argument("--smoke", action="store_true", help="짧게 (L4 배선 확인용)")
    args = ap.parse_args()

    global GSD_EVERY
    GSD_EVERY = args.gsd_every
    sys_ = load_system(ROOT / "intake/trap-drag-2d-hex300/system.yaml", args.v)
    if args.smoke:
        # ★ 배선 확인용. 물리 결과가 아니다 — 횡단의 1/200 만 돈다.
        args.traverse, args.samples = min(args.traverse, 0.01), min(args.samples, 40)
        args.warm, args.equil, args.relax = 0.5, 1.0, 2.0
    lg = build_ledger(sys_, dt_scale=args.dt_scale, n_traverse=args.traverse,
                      warm=args.warm, equil=args.equil, relax=args.relax)
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")

    # ── 네 구간 프로토콜 (사용자 질문: 평형 vs 구동 에너지 · 이완 · g(r)·결함) ──
    #   시간 단위는 전부 **τ_int**(격자의 국소 이완시간)입니다 — 끌기만 L_x/v_x 로 정해집니다.
    #   ★ 이완 구간을 τ_int 의 40배로 잡은 것은 ★제안입니다. 집단 모드가 τ_int 보다
    #     느리면 꼬리가 잘리므로, 런 후에 피팅된 τ 와 구간 길이를 반드시 대조해야 합니다.
    tau_int_steps = D["tau_int_steps"]
    n_warm, n_equil, n_relax = D["n_warm"], D["n_equil"], D["n_relax"]
    n_prod = D["n_drag"]
    sample_every = max(1, n_prod // args.samples) if n_prod else 1
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, Gamma, extra = analyze_scales(sys_, lg)
    tag = None
    if args.dt_scale != 1.0:
        tag = f"dt{args.dt_scale:g}"
    if args.traverse != 1.0:
        tag = (tag + "-" if tag else "") + f"tr{args.traverse:g}"
    if (args.warm, args.equil, args.relax) != (10.0, 20.0, 40.0):
        tag = (tag + "-" if tag else "") + f"w{args.warm:g}e{args.equil:g}r{args.relax:g}"
    if args.v is not None and abs(args.v - 0.5) > 1e-12:
        tag = (tag + "-" if tag else "") + f"v{args.v:g}"
    if args.seed != 20260804:
        tag = (tag + "-" if tag else "") + f"s{args.seed}"
    if args.smoke:
        tag = (tag + "-" if tag else "") + "smoke"

    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"A": sys_["A"], "phi": sys_["phi"], "N": sys_["N"],
                # ★ 직사각 박스 — L4는 정사각을 가정하면 안 된다.
                "Lx_star": D["Lx_star"], "Ly_star": D["Ly_star"],
                "n_x": D["n_x"], "n_y": D["n_y"], "a_nn_star": D["a_nn_star"],
                "lattice": "hex_commensurate", "drag_axis": "x",
                "r_c_star": D["r_c_star"],
                "wca_eps": sys_["wca_eps_kT"], "Gamma": Gamma,
                # 트랩·끌기를 무차원으로. L4는 이것만 보고 돌린다.
                "k_star": lg.ratio("energies", "k_t_d2", "kT"),
                "v_star": lg.ratio("times", "tau_B", "tau_v"),
                "n_trapped": sys_["n_trapped"],
                "r_table_min_star": R_TABLE_MIN},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "dt_over_tau_k": args.dt_scale * C.GATE,
                  # ★ 구간 길이는 **물리를 정하는 값**이라 해시에 들어갑니다 —
                  #   프로토콜이 다르면 다른 런입니다 (같은 계여도).
                  "n_warm": n_warm, "n_equil": n_equil,
                  "n_prod": n_prod, "n_relax": n_relax,
                  "n_samples": args.samples,
                  "sample_every": sample_every, "seed": args.seed},
        tag=tag, nhex=12)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 무결성 검사")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, extra, n_warm, n_equil, n_prod, n_relax,
                               tau_int_steps, args)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}   run_id={run_id}",
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

    # ── L4 — 스펙을 **되읽어서** 실행한다 (메모리의 객체가 아니라) ──────
    #   되읽는 것이 핵심입니다: 디스크의 스펙만으로 돌 수 있어야 L2↔L4 계약이 성립하고,
    #   해시 검증도 그때 걸립니다 (`execute` 가 verify_hash 를 먼저 봅니다).
    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.txt").write_text(report)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    if v["status"] == RUN.OK:
        print(f"플롯: {make_plots(outdir, loaded).relative_to(ROOT)}")
    return 0 if v["status"] in (RUN.OK, "skipped") else 1


if __name__ == "__main__":
    sys.exit(main())
