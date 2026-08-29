"""Phase 1-A — trap-2d-5um 관통.

물리계(SI) → 스케일 표 → 분리 검사 → 무차원화 → 실행 → 역변환 → 해석해 대조.

**Phase 1-C에서 공통 부분을 `bdbot/`로 옮겼습니다.** 여기 남은 것은 이 계 고유의 것뿐입니다:
조화 트랩 힘 · 앵커 변위 관측량 4종 · 해석해 · 평형 지표(앵커 변위) · 플롯.
무엇을 올리고 무엇을 남길지는 skill `bd-physics` §6.3 대조표로 판정했습니다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/trap_2d_5um.py              # 전체 실행 (~3분)
    $PY cases/trap_2d_5um.py --smoke      # 짧게 (~20초)
    $PY cases/trap_2d_5um.py --report     # 리포트만, 실행 안 함
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

from bdbot import Q, checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC  # noqa: E402
from bdbot import sim as SIM, stats as ST, traps as TR  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent


# ════════════════════════════════════════════════════════════════════════
# ① 물리계 (SI) — YAML에서 읽고 단위를 붙인다
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": load_node(raw["particle"]["diameter"]),
        "rho_p": load_node(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": load_node(raw["medium"]["temperature"]),
        "eta": load_node(raw["medium"]["viscosity"]),
        "k_t": load_node(raw["external"]["stiffness"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② 스케일 표 (bd-physics §0 ①)
#    ★ 케이스 고유: 트랩이 만드는 ℓ_k, τ_k. 지배 시간척도는 τ_B가 아니라 τ_k다.
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_: dict, L_box, num: dict) -> SC.ScaleLedger:
    """원장. ★ `dt`·`T_obs` 도 여기서 정합니다 — 둘 다 τ_k 에서 나오므로 원장의 일부입니다.

    예전에는 main() 이 dt/T_obs 를 따로 계산해 원장 밖에 두었습니다. 그러면 시간척도
    정렬 표에 둘이 안 나타나 "분리 위반이 눈에 보인다"는 원장의 목적이 반감됩니다
    (bdbot.scales.MANDATORY_ROLES).
    """
    d = sys_["d"].value.to("m")
    k = sys_["k_t"].value.to("N/m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma = b["kT"], b["gamma"]
    tau_k = C.relaxation_time(gamma, k)
    dt = C.dt_from_bias(tau_k, num["target_bias"])       # 편향에서 역산 (bd-physics §1.2)
    T_obs = num["production_tau"] * tau_k

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "입자 지름 (기준)")
    lg.add_length("l_k", ((kT / k) ** 0.5).to("m"), "√(kT/k) 트랩 요동폭")
    lg.add_length("L", L_box.to("m"), "박스", role="box")
    lg.add_time("tau_p", b["tau_p"], "m/γ 관성 이완", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_k", tau_k, "γ/k 트랩 이완", star=True)
    lg.add_time("tau_B", b["tau_B"], "d²/D_t 확산 (기준)")
    lg.add_time("T_obs", T_obs, "관측창", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.add_energy("k_d2", (k * d**2).to("J"), "k·d² 트랩 강성")
    lg.derived = {"gamma": gamma, "D_t": b["D_t"], "m": b["m"], "kT": kT, "k": k, "d": d,
                  "tau_k": tau_k, "dt": dt, "T_obs": T_obs}
    lg.ref = SC.thermal_reference(
        d, kT, b["tau_B"],
        SC.THERMAL_RATIONALE + " 이 계의 지배 시간척도는 τ_k다 — τ_B=242s는 트랩이 4ms에 "
        "붙잡아 실현되지 않는다. 결과 보고는 τ_k 단위로 재척도한다.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ③ 무차원수 + ④ 분리 검사 (bd-physics §3, §4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg):
    """무차원수와 분리 검사. dt·T_obs 는 원장에서 읽습니다 (더 이상 인자로 받지 않음)."""
    d, kT, k = lg.derived["d"], lg.derived["kT"], lg.derived["k"]
    tau_k, tau_p = lg.get("times", "tau_k"), lg.get("times", "tau_p")
    tau_B, dt, T_obs = lg.get("times", "tau_B"), lg.get("times", "dt"), lg.get("times", "T_obs")
    l_k, L = lg.get("lengths", "l_k"), lg.get("lengths", "L")

    f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)
    k_star = f(k * d**2 / kT)
    # num/den 을 붙인 것은 `NondimSpec.validate()` 가 원장에서 비를 재계산해 대조합니다.
    # k* 는 원장 두 항목의 비가 아니라 에너지비(k·d²/kT)라서 원장 기호로 표현됩니다.
    groups = [
        ND.Group("k*", k_star, ("energies", "k_d2"), ("energies", "kT"),
                 "k d²/kT", "트랩 vs 열요동"),
        ND.Group("l_k/d", f(l_k / d), ("lengths", "l_k"), ("lengths", "d"),
                 "1/√k*", "요동폭 vs 입자"),
        ND.Group("tau_k/tau_B", f(tau_k / tau_B), ("times", "tau_k"), ("times", "tau_B"),
                 "1/k*", "트랩 이완 vs 확산"),
        ND.Group("dt/tau_k", f(dt / tau_k), ("times", "dt"), ("times", "tau_k"),
                 "", "적분 해상"),
        ND.Group("T_obs/tau_k", f(T_obs / tau_k), ("times", "T_obs"), ("times", "tau_k"),
                 "", "관측창"),
    ]
    checks = [
        C.Check("model", "관성 무시   τ_p/τ_k", f(tau_p / tau_k), C.GATE, "<=",
                "BD(과감쇠)가 타당한가. dt와 무관 (bd-physics §4)"),
        C.Check("integration", "트랩 해상   dt/τ_k", f(dt / tau_k), C.GATE, "<=",
                f"편향 ≈ (dt/τ)/2 = {C.bias_from_dt(dt, tau_k):.3f}% (bd-physics §1.2)"),
        C.Check("geometry", "요동 vs 박스 2ℓ_k/L", f(2 * l_k / L), 0.5, "<=",
                "최소 이미지 안전 여유 (bd-hoomd 함정 1·6)"),
        C.Check("statistics", "관측창     T_obs/τ_k", f(T_obs / tau_k), 100.0, ">=",
                "정상상태 통계 충분성", hard=False),
    ]
    return groups, checks


def report_blocks(sys_, lg, n_steps):
    """리포트의 케이스별 블록 (공통 골격은 bdbot.report.render)."""
    tau_k, dt = lg.get("times", "tau_k"), lg.get("times", "dt")
    T_obs = lg.get("times", "T_obs")
    inp = [R.kv(key, f"{sys_[key].value:~.4gP}", sys_[key].tier, sys_[key].source[:44], val_w=20)
           for key in ("d", "T", "eta", "k_t", "rho_p")]
    inp.append(f"  N      = {sys_['N']}")
    der = [f"  {k_:<6} = {lg.derived[k_].to_compact():~.4gP}"
           for k_ in ("gamma", "D_t", "m", "kT")]
    plan = [
        f"  dt      = {dt.to_compact():~.4gP}   (= {float((dt / tau_k).to('')):.2e} τ_k)",
        f"  T_obs   = {T_obs.to_compact():~.4gP}",
        f"  steps   = {n_steps:,}   × N={sys_['N']}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# ⑤ 해석해 (골든 테스트 기준) — ★ 케이스 고유. 1-B에는 해석해가 없다.
# ════════════════════════════════════════════════════════════════════════
def analytic(lg):
    kT, k, gamma = lg.derived["kT"], lg.derived["k"], lg.derived["gamma"]
    tau_k = C.relaxation_time(gamma, k)
    x2 = (kT / k).to("um^2")
    return {
        "x2": x2,                                      # ⟨x²⟩ 자유도당
        "sigma": (x2 ** 0.5).to("nm"),                 # P(x) 폭
        "tau_k": tau_k,
        "f_c": (1 / (2 * math.pi * tau_k)).to("Hz"),   # PSD corner
        "S0": (4 * x2 * tau_k).to("um^2/Hz"),          # PSD f→0 극한
    }


# ════════════════════════════════════════════════════════════════════════
# ⑥ L4 — 스펙만 읽고 계를 세운다 (bdbot.run 이 돌리고 판정한다)
#    ★ 케이스 고유: 조화 트랩(고정) + 앵커 변위로부터의 관측량 4종
# ════════════════════════════════════════════════════════════════════════
def analytic_star(k_star: float) -> dict:
    """무차원(σ=d, E=kT, τ=τ_B) 닫힌 형태. `analytic(lg)` 와 같은 식 — k*만 있으면 된다.

    두 곳(플롯용 `analytic(lg)`, 관측량 대조용 `finalize`)이 각자 값을 만들면 갈라질 수
    있어 여기 하나로 모은다.
    """
    return {"x2_star": 1.0 / k_star, "tau_k_star": 1.0 / k_star,
            "fc_star": k_star / (2 * math.pi), "S0_star": 4.0 / k_star ** 2}


@RUN.builder("trap-2d-5um")
def build(spec, outdir=None) -> RUN.Build:
    """스펙 → HOOMD 계. 트랩은 `bdbot.traps.make_trap` (velocity·drive 없음 = 고정).

    ★ 이전에는 이 케이스가 자기 안에 `HarmonicTrap` 을 따로 갖고 있었습니다 — `traps.py`가
      "세 번 나와서 올렸다"고 적어둔 것이 바로 이 클래스입니다. 최소 이미지·tag 인덱싱
      규약이 정확히 같아 (bdbot/sim.py의 `minimum_image` 와 동일한 `period` 계산) 여기서
      교체해도 결과가 바뀌지 않습니다 (아래에서 재실행으로 대조).
    """
    P, Nm = spec.params, spec.numerics
    N, L_star, k_star = int(P["N"]), float(P["L_star"]), float(P["k_star"])
    dt_star = float(Nm["dt_star"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])
    seed = int(Nm["seed"])

    n_side = int(math.ceil(math.sqrt(N)))
    a = L_star / n_side
    pos0 = np.array([[(i % n_side + .5) * a - L_star / 2,
                      (i // n_side + .5) * a - L_star / 2, 0.0] for i in range(N)])
    sim = SIM.make_sim(SIM.frame_2d(pos0, L_star), seed=seed)

    trap = TR.make_trap(k_star, pos0, L_star, dt_star=dt_star)
    SIM.attach_brownian(sim, dt_star, [trap])
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, max(1, n_prod // 200))

    def pe_pp():
        return float(np.array(trap.energies).sum()) / N

    # ── 표본 누적기 — 원래(1-C 이전) 구현과 같은 메모리 사용량을 유지한다 ──
    #   `sample()`이 매 표본의 전체 N-입자 배열을 그대로 `cols`에 돌려주면 RUN.execute가
    #   그걸 통째로 observables.npz에 저장한다(스택돼서 (n_samp,N,2)). 1000입자×2e4표본
    #   이면 148MB — 원래는 448KB였다(파생량만 저장). 그래서 P(x)·⟨x²⟩는 여기 클로저에서
    #   즉시 누적하고, C(t)·PSD용 부분집합(n_trace개)만 미리 잡은 배열에 채운다 — 둘 다
    #   `cols`를 거치지 않으므로 저장되지 않는다.
    n_samp = n_prod // sample_every if sample_every else 0
    n_trace = min(250, N)
    trace = np.empty((max(n_samp, 1), n_trace, 2), dtype=np.float32)
    sum_x2 = np.zeros(2)
    per_sample_x2 = np.empty(max(n_samp, 1))         # 블록SEM 용 (⟨d²⟩ 전체 N 기준)
    hist_edges = np.linspace(-6 / math.sqrt(k_star), 6 / math.sqrt(k_star), 121)
    hist = np.zeros(len(hist_edges) - 1)
    i_sample = [0]

    def sample(timestep, phase):
        dxy = trap.displacement(sim.state, timestep)[:, :2]
        i = i_sample[0]
        if i < n_samp:
            trace[i] = dxy[:n_trace]
        sum_x2[:] += (dxy ** 2).mean(axis=0)
        if i < n_samp:
            per_sample_x2[i] = float((dxy ** 2).mean())
        hist[:] += np.histogram(dxy.ravel(), bins=hist_edges)[0]
        i_sample[0] = i + 1
        return {}

    def finalize(cols):
        from scipy import optimize, signal

        tr = trace[:i_sample[0]]                    # C(t)·PSD 용 (원래도 이 부분집합만 썼다)
        dt_sample_star = dt_star * sample_every

        x2 = float(sum_x2.sum() / (2 * i_sample[0]))  # ⟨x²⟩ 자유도당 — 전체 N 입자 기준
        x2_sem = ST.block_sem(per_sample_x2[:i_sample[0]], 20)

        # P(x) — 전체 N·전체 표본의 히스토그램 (누적)
        centers = 0.5 * (hist_edges[1:] + hist_edges[:-1])
        px = hist / (hist.sum() * (hist_edges[1] - hist_edges[0]))

        # C(t) — 표본평균을 빼지 않는다 (bd-physics §5.1)
        ac = ST.autocorr_unbiased(tr)
        t = np.arange(len(ac)) * dt_sample_star
        tau_guess = 1.0 / k_star
        fit_n = min(max(10, int(3 * tau_guess / dt_sample_star)), len(ac))
        popt, _ = optimize.curve_fit(lambda tt, A, tau: A * np.exp(-tt / tau),
                                     t[:fit_n], ac[:fit_n],
                                     p0=[ac[0], tau_guess], maxfev=20000)
        tau_fit = float(abs(popt[1]))

        # PSD (one-sided density; ∫S df = variance)
        x = tr[:, :, 0].astype(np.float64).T
        fs = 1.0 / dt_sample_star
        nper = min(len(tr), 4096)
        f, S = signal.welch(x, fs=fs, nperseg=nper, axis=-1, detrend="constant")
        S = S.mean(axis=0)
        psd_f, psd = f[1:], S[1:]

        def lorentz(ff, S0, fc):
            return S0 / (1 + (ff / fc) ** 2)

        try:
            popt, _ = optimize.curve_fit(lorentz, psd_f, psd,
                                         p0=[psd[0], 1.0 / (2 * math.pi / k_star)],
                                         maxfev=20000)
            psd_S0, psd_fc = float(popt[0]), float(abs(popt[1]))
        except Exception as e:
            psd_S0, psd_fc = float("nan"), float("nan")
            print(f"    (PSD 피팅 실패: {e})")

        # ── 해석해 대조 (역할: implementation_check·module — bd-physics §7.5) ──
        #   단일 트랩+BD는 시뮬레이션이 그대로 푸는 모델의 닫힌 형태다. 조합 가정이 없다.
        ana = analytic_star(k_star)
        # (표시 이름, 측정*, 예측*, L, T, 단위) — 표시 이름·단위는 원래(1-C 이전) 구현과
        # 그대로 맞춘다 (postmortem·리포트가 이 문자열로 사람이 읽는다).
        rows = [
            ("⟨x²⟩", x2, ana["x2_star"], 2, 0, "µm²"),
            ("σ = √⟨x²⟩", math.sqrt(x2), math.sqrt(ana["x2_star"]), 1, 0, "nm"),
            ("τ (C(t) 피팅)", tau_fit, ana["tau_k_star"], 0, 1, "ms"),
            ("f_c (PSD 피팅)", psd_fc, ana["fc_star"], 0, -1, "Hz"),
            ("S(0) (PSD 피팅)", psd_S0, ana["S0_star"], 2, 1, "µm²/Hz"),
        ]
        obs = []
        for name, meas_star, pred_star, L, T, unit in rows:
            meas = float(spec.physical(meas_star, L=L, T=T).to(unit).magnitude)
            pred = float(spec.physical(pred_star, L=L, T=T).to(unit).magnitude)
            obs.append(MET.observable(
                name, meas, pred, unit, "analytic", role="implementation_check",
                scope="module", tol_pct=5.0,
                note="단일 조화 트랩 + BD — OU 과정의 닫힌 형태 (조합 가정 없음)"))

        bias_pct = 50.0 * dt_star * k_star     # bias ≈ (dt/τ_k)/2, τ_k* = 1/k*
        return {"observables": obs,
                "extra": {"x2_star": x2, "x2_sem_pct": 100 * x2_sem / x2 if x2 else None,
                          "tau_fit_star": tau_fit, "psd_fc_star": psd_fc,
                          "psd_S0_star": psd_S0, "bias_predicted_pct": bias_pct,
                          "stat_target_pct": 0.5},
                "arrays": {"px_centers": centers, "px": px, "ac_t": t, "ac": ac,
                          "psd_f": psd_f, "psd": psd}}

    return RUN.Build(
        sim=sim, forces=[trap], n_particles=N,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq, n_prod=n_prod, sample_every=sample_every,
        tags=["2D", "harmonic_trap", "newtonian", "single_particle", "no_pair_interaction"],
        physical={"N": N, "k_star": k_star, "L_star": L_star},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="L3 스펙을 specs/<run_id>.json 으로 쓰고 종료 (실행 안 함)")
    ap.add_argument("--force", action="store_true", help="같은 run_id 결과가 있어도 재실행")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/trap-2d-5um/system.yaml")
    num = sys_["numerics"]
    if args.smoke:
        sys_["N"], num = 200, dict(num, production_tau=120, equilibration_tau=10)

    # 박스: 격자 간격 = 1d, 요동폭(~0.004d)보다 압도적으로 큼
    n_side = int(math.ceil(math.sqrt(sys_["N"])))
    L_star = float(n_side)
    lg = build_ledger(sys_, Q(L_star, "dimensionless") * sys_["d"].value, num)

    tau_k, dt = lg.get("times", "tau_k"), lg.get("times", "dt")
    tau_B, T_obs = lg.ref["time"][1], lg.get("times", "T_obs")
    n_eq = int(round(float((num["equilibration_tau"] * tau_k / dt).to(""))))
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, int(round(float((tau_k / num["samples_per_tau"] / dt).to("")))))
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks = analyze_scales(sys_, lg)
    # ★ L3 산출물. `run_id` 해시는 {system, params, numerics} 뿐이고 physics_only 가
    #   적용됩니다 (주석·출처를 고쳐서 런이 무효화되면 콘텐츠 주소가 쓸모없다).
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"N": sys_["N"], "L_star": L_star,
                "k_star": float((lg.derived["k"] * lg.derived["d"] ** 2
                                 / lg.derived["kT"]).to(""))},
        numerics={"dt_star": float((dt / tau_B).to("")),
                  "dt_over_tau_k": float((dt / tau_k).to("")),
                  "n_eq": n_eq, "n_prod": n_prod, "sample_every": sample_every,
                  "seed": 20260803},
        nhex=12)
    run_id = spec.run_id()

    # L3 무결성 (원장 완전성 · 무차원수가 정말 그 비인가) — 물리 검사와 다른 층입니다.
    l3 = spec.validate()
    if l3:
        print("L3 무결성 검사")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)
    if spec.errors:
        print(f"\n❌ L3 무결성 오류 {len(spec.errors)}건 — 무차원화가 성립하지 않습니다.")
        return 1
    if verdict == "FAIL":
        print("\n❌ 하드 분리 검사 실패 — 실행하지 않습니다.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 스펙: {p.relative_to(ROOT)}")
        return 0

    # ── L4 — 디스크의 스펙을 되읽어서 실행한다 (bdbot.run, 해시 검증이 그때 걸린다) ──
    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    if v["status"] == "skipped":
        return 0
    if v["status"] != RUN.OK:
        return 1
    # ★ report.txt는 execute() 이후에 쓴다 — execute()가 시작할 때 outdir를 청소하므로
    #   먼저 쓰면 지워진다 (prepare_outdir가 result.txt/record.json 외 전부 지운다).
    (outdir / "report.txt").write_text(report)

    # ── 해석해 대조 — finalize()가 이미 계산해 metrics.json에 넣어둔 것을 표로 찍는다 ──
    #    (이 케이스의 핵심 산출물이라 콘솔·result.txt에도 남긴다)
    m = json.loads((outdir / "metrics.json").read_text())
    obs_out = m["observables"]
    res_extra = m.get("result", {})
    lines = ["", "=" * R.W, "결과 — 물리 단위로 역변환 후 해석해 대조", "=" * R.W,
             f"{'관측량':<18}{'측정':>18}{'해석해':>18}{'오차':>10}   판정"]
    all_ok = True
    for o in obs_out:
        ok = o["err_pct"] is not None and abs(o["err_pct"]) < (o["tol_pct"] or 5.0)
        all_ok &= ok
        lines.append(f"{o['name']:<18}{o['measured']:>18.6g}{o['predicted']:>18.6g}"
                     f"{o['err_pct']:>+9.2f}%   {'✓' if ok else '✗'}")
    lines += ["",
              f"  (⟨x²⟩ 예상 계통 편향 = (dt/τ_k)/2 = "
              f"+{res_extra.get('bias_predicted_pct', float('nan')):.3f}%  — bd-physics §1.2)",
              f"  (⟨x²⟩ 통계 오차 ±{res_extra.get('x2_sem_pct', float('nan')):.3f}%)",
              "=" * R.W,
              f"VERDICT: {'✓ PASS — 4개 관측량 모두 해석해와 일치' if all_ok else '✗ FAIL'}",
              "=" * R.W]
    result = "\n".join(lines)
    print(result)
    (outdir / "result.txt").write_text(report + "\n" + result)

    ana = analytic(lg)
    make_plots(ana, lg, outdir)
    print("\n".join(RID.list_artifacts(outdir, ROOT)))
    return 0 if all_ok else 1


def make_plots(ana, lg, outdir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    matplotlib.rcParams["mathtext.fontset"] = "dejavusans"

    obs = np.load(outdir / "observables.npz")
    d, tau_B = lg.derived["d"], lg.get("times", "tau_B")
    tau_k = ana["tau_k"]
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))

    xc = obs["px_centers"] * float(d.to("nm").magnitude)
    sig = float(ana["sigma"].to("nm").magnitude)
    ax[0, 0].plot(xc, obs["px"] / float(d.to("nm").magnitude), "o", ms=3, label="측정")
    xx = np.linspace(xc.min(), xc.max(), 400)
    ax[0, 0].plot(xx, np.exp(-xx**2 / (2 * sig**2)) / (sig * math.sqrt(2 * math.pi)),
                  "-", lw=2, label=f"Gaussian σ={sig:.2f} nm")
    ax[0, 0].set(xlabel="x [nm]", ylabel="P(x) [1/nm]", title="① 위치 분포")
    ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

    ax[0, 1].semilogy(xc, obs["px"] / float(d.to("nm").magnitude), "o", ms=3)
    ax[0, 1].semilogy(xx, np.exp(-xx**2 / (2 * sig**2)) / (sig * math.sqrt(2 * math.pi)),
                      "-", lw=2)
    ax[0, 1].set(xlabel="x [nm]", ylabel="P(x)",
                 title="② 위치 분포 (log) — 꼬리까지 Gaussian?")
    ax[0, 1].grid(alpha=.3)

    t_ms = obs["ac_t"] * float(tau_B.to("ms").magnitude)
    c_um2 = obs["ac"] * float((d**2).to("um^2").magnitude)
    n = min(len(t_ms), int(6 * float(tau_k.to("ms").magnitude) / max(t_ms[1], 1e-12)))
    ax[1, 0].semilogy(t_ms[1:n], c_um2[1:n], "o", ms=3, label="측정")
    tt = np.linspace(0, t_ms[n - 1], 300)
    ax[1, 0].semilogy(tt, float(ana["x2"].to("um^2").magnitude)
                      * np.exp(-tt / float(tau_k.to("ms").magnitude)),
                      "-", lw=2, label=f"exp(−t/τ), τ={float(tau_k.to('ms').magnitude):.3f} ms")
    ax[1, 0].set(xlabel="t [ms]", ylabel="⟨x(0)x(t)⟩ [µm²]", title="③ 위치 자기상관")
    ax[1, 0].legend(); ax[1, 0].grid(alpha=.3)

    f_hz = obs["psd_f"] / float(tau_B.to("s").magnitude)
    s_phys = obs["psd"] * float((d**2 * tau_B).to("um^2*s").magnitude)
    ax[1, 1].loglog(f_hz, s_phys, "-", lw=1, alpha=.7, label="측정")
    S0, fc = float(ana["S0"].to("um^2/Hz").magnitude), float(ana["f_c"].to("Hz").magnitude)
    ax[1, 1].loglog(f_hz, S0 / (1 + (f_hz / fc) ** 2), "-", lw=2,
                    label=f"Lorentzian, $f_c$={fc:.1f} Hz")
    ax[1, 1].axvline(fc, ls="--", c="k", alpha=.5)
    ax[1, 1].set(xlabel="f [Hz]", ylabel="S(f) [µm²/Hz]", title="④ 파워 스펙트럼")
    ax[1, 1].legend(); ax[1, 1].grid(alpha=.3, which="both")

    for a in ax.ravel():
        a.title.set_fontsize(10)
    fig.suptitle("trap-2d-5um — 측정 vs 해석해 (물리 단위)", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
