"""Phase 1-B — soft-r3-2d-A-sweep 관통.

물리계(SI) → 스케일 표 → 분리 검사 → 무차원화 → 실행 → 역변환 → 검증.

**Phase 1-C에서 공통 부분을 `bdbot/`로 옮겼습니다.** 여기 남은 것은 이 계 고유의 것뿐입니다:
r⁻³ 표 퍼텐셜 + WCA 코어 · 최근접 접근거리에서 dt 역산 · RSA 초기배치 ·
구조 관측량(g(r)·ψ₆·Voronoi) · 검증 3종 · 플롯.
판정 근거는 skill `bd-physics` §6.3 대조표입니다.

1-A와 다른 점: 해석해가 없습니다. 대신 세 가지로 검증합니다.
  ① 희박극한  g(r) → exp(−βU)      퍼텐셜 구현이 맞는가 (--dilute)
  ② 에너지 일관성  ⟨U⟩/N  vs  (ρ/2)∫U g(r) 2πr dr
  ③ 최소 이웃거리 감시 — pair.Table 함정 11 (r<r_min에서 힘이 0) 방어

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/soft_r3_2d.py --A 100            # 한 개 A (독립 런 → 병렬 가능)
    $PY cases/soft_r3_2d.py --A 1 --smoke      # 짧게
    $PY cases/soft_r3_2d.py --dilute           # 희박극한 검증 런
    $PY cases/soft_r3_2d.py --A 100 --report   # 리포트만
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
from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC  # noqa: E402
from bdbot import sim as SIM, stats as ST  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
# ★ 퍼텐셜 수치(U·U''·r_min)는 `trap-drag` 가 같은 것을 쓰면서 bdbot 으로 올렸습니다.
#   여기서 재정의하면 두 케이스의 dt 가 갈라집니다 (bdbot/pairpot.py 참조).
from bdbot.pairpot import HEX_NN, R_WCA, U2_star, U_star, approach_distance  # noqa: E402,F401

ROOT = Path(__file__).resolve().parent.parent
R_TABLE_MIN = 0.5          # pair.Table 하한 (trap-drag와 동일). 함정 11 방어의 기준


# ════════════════════════════════════════════════════════════════════════
# ① 물리계 (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    P = load_node
    r3 = raw["interactions"][0]
    wca = raw["interactions"][1]
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "phi": float(raw["geometry"]["area_fraction"]["value"]),
        "A_list": list(r3["amplitude_A"]["value"]),
        "r_c": P(r3["cutoff"]),
        "wca_eps_kT": float(wca["epsilon"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② 스케일 표 (bd-physics §0 ①)
#    1-A와 비교: 길이가 3개→5개로 늘고, τ_int가 상수가 아니라 r의 함수다.
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, A, N, phi, r_c_star, dt_scale=1.0, T_obs_tau=None) -> SC.ScaleLedger:
    """원장. ★ `dt`·`T_obs` 도 원장에 올립니다 — 시간척도 정렬에서 보이게 하려고.

    `dt` 는 τ_int(r_min) 에서, `T_obs` 는 τ_B 배수로 나옵니다. 둘 다 이 원장의 값에서
    유도되므로 원장 밖에 두면 "분리 위반이 눈에 보인다"는 목적이 반감됩니다.
    """
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, D_t, tau_B, m = b["kT"], b["gamma"], b["D_t"], b["tau_B"], b["m"]

    a_star = math.sqrt(math.pi / (4 * phi))             # a_mean/d
    L_star = a_star * math.sqrt(N)
    r_min_star, crit, u_rms_rel, state = approach_distance(A, a_star, sys_["wca_eps_kT"])
    # τ_int = γ/U''(r_min) — 트랩의 τ_k = γ/k 와 같은 구조 (bdbot.checks.relaxation_time).
    # 강성이 kT/d² 단위의 무차원값이라 여기서는 τ_B로 환산해 쓴다.
    tau_int = (tau_B / float(U2_star(r_min_star, A, sys_["wca_eps_kT"]))).to("s")

    dt = dt_scale * 1e-2 * tau_int                       # 기본은 하드 게이트에 딱 맞춤
    T_obs = (T_obs_tau if T_obs_tau is not None
             else float(sys_["numerics"]["production_tau_B"])) * tau_B

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "입자 지름 (기준)")
    lg.add_length("r_min", r_min_star * d, "최근접 접근거리", star=True)
    lg.add_length("a_mean", a_star * d, "평균 간격")
    lg.add_length("r_c", r_c_star * d, "컷오프")
    lg.add_length("L", L_star * d, "박스", role="box")
    lg.add_time("tau_p", b["tau_p"], "m/γ 관성 이완", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_int", tau_int, "γ/U''(r_min) 상호작용", star=True)
    lg.add_time("tau_B", tau_B, "d²/D_t 확산 (기준)")
    lg.add_time("T_obs", T_obs, "관측창", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.add_energy("U_a", (float(U_star(a_star, A, sys_["wca_eps_kT"])) * kT).to("J"),
                  "U(a_mean) 평균간격 결합 = Γ kT")
    # ★ U(d) = A + ε_WCA 다 — **A kT 가 아니다.** r=d 에서 WCA 코어가 정확히 ε 을 냅니다
    #   (4ε(1−1)+ε). 예전 원장은 이 항목을 "= A kT" 로 적어뒀고 A=100 에서 1% 어긋났습니다.
    #   L3 무결성 검사(groups.A ≠ U_d/kT = 101)가 이 라벨 오류를 잡았습니다.
    lg.add_energy("U_d", (float(U_star(1.0, A, sys_["wca_eps_kT"])) * kT).to("J"),
                  "U(d) 접촉 결합 = (A+ε_WCA) kT")
    lg.derived = dict(gamma=gamma, D_t=D_t, m=m, kT=kT, d=d, tau_B=tau_B,
                      a_star=a_star, L_star=L_star, r_min_star=r_min_star,
                      crit=crit, u_rms_rel=u_rms_rel, state=state, tau_int=tau_int,
                      dt=dt, T_obs=T_obs)
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " 1-A와 달리 여기서는 τ_B가 실제 지배 시간척도이고, "
        "dt는 τ_int(r_min)으로 정한다.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ③ 무차원수 + ④ 분리 검사 (bd-physics §3, §4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg, A, phi, r_c_star):
    """무차원수와 분리 검사. dt·T_obs 는 원장에서 읽습니다 (더 이상 인자로 받지 않음)."""
    D = lg.derived
    f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)
    a_star, L_star = D["a_star"], D["L_star"]
    Gamma = float(U_star(a_star, A, sys_["wca_eps_kT"]))
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    tau_p = lg.get("times", "tau_p")

    # ★ Γ 는 진짜 제어 파라미터입니다 (A 단독이 아님 — bd-physics §6.2).
    #   원장의 U_a/kT 와 대조되므로 A(d/a)³ 계산이 어긋나면 validate() 가 잡습니다.
    groups = [
        ND.Group("Gamma", Gamma, ("energies", "U_a"), ("energies", "kT"),
                 "U(a_mean)/kT", "결합 vs 열요동 ★"),
        # A 는 r⁻³ 항의 진폭(입력)이고 원장 두 항목의 비가 **아니다** — U(d)/kT = A+ε 이다.
        ND.Group("A", A, None, None, "", "r⁻³ 진폭 (입력, 스케치값)"),
        ND.Group("U(d)/kT", float(U_star(1.0, A, sys_["wca_eps_kT"])),
                 ("energies", "U_d"), ("energies", "kT"), "(A+ε_WCA)", "접촉 결합"),
        ND.Group("phi", phi, None, None, "", "밀집도"),
        ND.Group("a_mean/d", a_star, ("lengths", "a_mean"), ("lengths", "d"), "", "평균간격"),
        ND.Group("L/d", L_star, ("lengths", "L"), ("lengths", "d"), "", "박스 크기"),
        ND.Group("r_c/d", r_c_star, ("lengths", "r_c"), ("lengths", "d"), "", "컷오프"),
        ND.Group("r_c/a_mean", r_c_star / a_star, ("lengths", "r_c"), ("lengths", "a_mean"),
                 "", "컷오프(이웃 껍질 수)"),
        ND.Group("dt/tau_int", f(dt / D["tau_int"]), ("times", "dt"), ("times", "tau_int"),
                 "", "적분 해상"),
        ND.Group("T_obs/tau_B", f(T_obs / D["tau_B"]), ("times", "T_obs"), ("times", "tau_B"),
                 "", "관측창"),
        ND.Group("St", f(tau_p / D["tau_B"]), ("times", "tau_p"), ("times", "tau_B"),
                 "tau_p/tau_B", "관성 vs 확산"),
    ]
    checks = [
        C.Check("모델", "관성 무시    τ_p/τ_int",
              f(tau_p / D["tau_int"]), C.GATE, "<=",
              "BD(과감쇠)가 타당한가. dt와 무관 (bd-physics §4)"),
        C.Check("적분", "상호작용 해상 dt/τ_int", f(dt / D["tau_int"]), C.GATE, "<=",
              f"τ_int = γ/U''(r_min={D['r_min_star']:.3f}d), {D['crit']} 기준. "
              f"선형계 기준 편향 ≈ {C.bias_from_dt(dt, D['tau_int']):.3f}% — 비선형이라 수렴확인 별도"),
        C.Check("기하", "컷오프       r_c/(L/2)", r_c_star / (L_star / 2), 1.0, "<=",
              "최소 이미지 (bd-hoomd 함정 6). 위반 시 과거 +1856% 사례"),
        C.Check("기하", "코어 여유    r_table_min/r_min",
              R_TABLE_MIN / D["r_min_star"], 1.0, "<=",
              "pair.Table 함정 11: r<r_min이면 힘이 0. 접근거리가 표 하한 위인가"),
        C.Check("통계", "관측창       T_obs/τ_B", f(T_obs / D["tau_B"]), 100.0, ">=",
              "구조 완화 통계. τ_B 기준 — 강결합 결함 어닐링은 더 느릴 수 있음",
              hard=False),
    ]
    return groups, checks, Gamma


def report_blocks(sys_, lg, A, phi, N, n_eq, n_prod):
    """리포트의 케이스별 블록 (공통 골격은 bdbot.report.render)."""
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    inp = [R.kv(key, f"{sys_[key].value:~.4gP}", sys_[key].tier, sys_[key].source[:46])
           for key in ("d", "T", "eta", "rho_p")]
    inp += [
        R.kv("A", f"{A}", 0, "sketch 'A = 0.1, 1, 10, 100' (무차원 해석 — 확인 필요)"),
        R.kv("phi", f"{phi}", 1, "사용자 확인 2026-08-03 (스케치 미기재)"),
        R.kv("N", f"{N}", 1, "스케치 100 → 최소이미지 위해 상향"),
    ]
    der = [f"  {k_:<8} = {D[k_].to_compact():~.4gP}" for k_ in ("gamma", "D_t", "m")]
    der += [
        f"  상태 추정 = {D['state']}  (Lindemann σ_bond/a_NN = {D['u_rms_rel']:.4f}, 기준 0.15)",
        f"  a_NN(육방 예측) = {HEX_NN*D['a_star']:.4f} d   ← a_mean 이 아님",
    ]
    plan = [
        f"  dt      = {dt.to_compact():~.4gP}  = {float((dt/D['tau_B']).to('')):.3e} τ_B",
        f"  T_obs   = {T_obs.to_compact():~.4gP}  = {float((T_obs/D['tau_B']).to('')):.1f} τ_B",
        f"  steps   = eq {n_eq:,} + prod {n_prod:,}   × N={N}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# ⑥ 실행 (무차원 단위)
# ════════════════════════════════════════════════════════════════════════
def rsa_positions(N, L, min_sep, rng):
    """랜덤 순차 배치 — 격자에서 시작하면 구조 결과가 편향된다."""
    pos = np.empty((N, 2))
    n = 0
    tries = 0
    while n < N:
        tries += 1
        if tries > 400 * N:
            raise RuntimeError(f"RSA 실패: {n}/{N} (min_sep={min_sep} 너무 큼)")
        p = rng.uniform(-L / 2, L / 2, 2)
        if n:
            dr = pos[:n] - p
            dr -= L * np.round(dr / L)
            if (dr**2).sum(axis=1).min() < min_sep**2:
                continue
        pos[n] = p
        n += 1
    return pos


RDF_BINS = 300             # freud g(r) 빈 수 — 출력 설정이지 물리가 아니라 스펙에 안 넣는다


@RUN.builder("soft-r3-2d-A-sweep")
def build(spec, outdir=None) -> RUN.Build:
    """스펙 → HOOMD 계. r⁻³ 소프트 반발 + WCA 코어, RSA 초기배치.

    ★ 검증 3종(에너지 일관성 · 육방 NN거리 · 코어 여유)은 여기 `finalize`가 계산합니다.
      해석해가 없는 케이스라 `checks`(하드/소프트 분리 검사)와는 다른 층입니다 — 원래
      `post_checks`로 부르던 사후 가드는 `bdbot.run`의 checks 스키마에 아직 자리가 없어
      (두 번째로 필요해지면 그때 올립니다) `extra.post_checks`에 `Check.as_dict()`로 담습니다.
    """
    import freud
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    N, A = int(P["N"]), float(P["A"])
    phi, r_c_star, eps = float(P["phi"]), float(P["r_c_star"]), float(P["wca_eps"])
    a_star = math.sqrt(math.pi / (4 * phi))
    L_star = a_star * math.sqrt(N)
    r_min_star = approach_distance(A, a_star, eps)[0]
    dt_star, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_eq, n_prod = int(Nm["n_eq"]), int(Nm["n_prod"])
    sample_every = int(Nm["sample_every"])

    np_seed, _ = SIM.resolve_seed(seed)
    rng = np.random.default_rng(np_seed)
    pos = rsa_positions(N, L_star, 1.0, rng)          # ★ 격자에서 시작하면 구조가 편향된다
    sim = SIM.make_sim(SIM.frame_2d(pos, L_star), seed=seed)

    cell = md.nlist.Cell(buffer=0.4)
    # r⁻³ 꼬리 — pair.Table. ★ endpoint=False (함정 10), 컷오프에서 시프트
    nbins = 1000
    rr = np.linspace(R_TABLE_MIN, r_c_star, nbins, endpoint=False)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_c_star)
    tab.params[("A", "A")] = dict(r_min=R_TABLE_MIN, U=A / rr**3 - A / r_c_star**3,
                                  F=3 * A / rr**4)
    # 배제부피 코어 — 별도 힘. 함정 11(r<r_min 힘 0) 방어의 1차선
    wca = SIM.wca(cell, epsilon=eps, sigma=1.0)                        # 함정 4
    SIM.attach_brownian(sim, dt_star, [tab, wca])
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, max(1, n_prod // 200))

    box = freud.box.Box.square(L_star)
    rdf = freud.density.RDF(bins=RDF_BINS, r_max=min(r_c_star, L_star / 2 - 1e-6), r_min=0.3)
    hexatic = freud.order.Hexatic(k=6, weighted=True)
    voro = freud.locality.Voronoi()
    # ★ 원래(1-C 이전) 구현과 같은 메모리 사용량을 유지한다 — 배위수 히스토그램은
    #   여기 클로저에서 즉시 누적한다 (`sample()`이 매번 13칸을 돌려주면 될 일이지만
    #   원래 코드처럼 한 번에 정규화하는 편이 finalize와 대칭이라 그대로 둔다).
    coord_hist = np.zeros(13)

    def xy():
        return np.array(sim.state.get_snapshot().particles.position, dtype=float)

    def pe_pp():
        return float(np.array(tab.energies).sum() + np.array(wca.energies).sum()) / N

    def sample(timestep, phase):
        p = xy()
        rdf.compute((box, p), reset=False)
        vn = voro.compute((box, p)).nlist
        counts = np.asarray(vn.neighbor_counts)
        coord_hist[:] += np.bincount(np.clip(counts, 0, 12), minlength=13)
        dists = np.asarray(vn.distances)
        psi6 = float(np.abs(hexatic.compute((box, p), neighbors=vn).particle_order).mean())
        return {"psi6": psi6, "min_sep": float(dists.min()),
                "bond_mean": float(dists.mean()), "bond_std": float(dists.std())}

    def finalize(cols):
        psi6 = float(cols["psi6"].mean())
        psi6_sem = ST.block_sem(cols["psi6"])
        min_sep = float(cols["min_sep"].min())
        bond_mean = float(cols["bond_mean"].mean())
        bond_std = float(cols["bond_std"].mean())
        pe_mean = float(cols["pe"].mean())
        pe_sem = ST.block_sem(cols["pe"])
        rdf_r, rdf_g = np.array(rdf.bin_centers), np.array(rdf.rdf)
        pe_rdf = energy_from_rdf(rdf_r, rdf_g, A, eps, phi, r_c_star)

        # ── 검증 — 해석해가 없으므로 극한·일관성으로 (bd-physics §7.6) ──────
        obs = []
        # 에너지 항등식은 g(r) 이 무엇이든 성립하는 **정확한 식**이다(근사 없음) —
        # 그래서 결정/유체 어느 쪽이든 implementation_check 다. 불일치 = g(r) 빈닝/퍼텐셜
        # 정의 버그. scope=composite 인데도 유도가 그대로 옮겨오는 이유가 이것이다(원칙 9.1).
        obs.append(MET.observable(
            "에너지 일관성 ⟨U⟩/N", pe_mean, pe_rdf, "kT", "consistency",
            role="implementation_check", scope="composite", tol_pct=2.0,
            note="consistency: HOOMD 힘 합 vs (ρ/2)∫U g(r) 2πr dr",
            derivation="⟨U⟩/N = (ρ/2)∫U(r)g(r)2πr dr 는 g(r)이 무엇이든 성립하는 항등식"
                       "(평균장·희박 근사가 아니다) — 결정/유체 어느 구조에도 그대로 옮겨온다"))
        # 반면 완전 육방 결정을 가정하는 것은 시뮬레이션이 강제하지 않은 가정이다 —
        # r⁻³ 퍼텐셜이 실제로 그 구조로 응결하는지가 이 케이스의 관심사라 hypothesis 다.
        a_nn_pred = HEX_NN * a_star
        if psi6 > 0.6:      # 결정일 때만 의미 있는 예측
            obs.append(MET.observable(
                "육방 NN 거리", bond_mean, a_nn_pred, "d", "lattice",
                role="hypothesis", tol_pct=2.0,
                note="lattice: a_NN = √(2/√3)·a_mean (완전 육방). ψ₆>0.6 일 때만 적용"))

        # ── 사후 가드 — 런이 설계 가정 안에 머물렀는가 ──────────────────
        post_checks = [
            C.Check("기하", "표 하한 여유  r_tab_min/min_sep", R_TABLE_MIN / min_sep, 1.0, "<=",
                  f"pair.Table 함정 11: r<{R_TABLE_MIN}d 면 힘이 0. 측정 최소 {min_sep:.3f}d"),
            C.Check("적분", "설계 r_min 준수 r_min/min_sep",
                  r_min_star / min_sep, 1.15, "<=",
                  f"dt는 r_min={r_min_star:.3f}d 의 국소 강성에서 정했다. "
                  f"측정 최소 {min_sep:.3f}d — 크게 안쪽이면 dt 재검토"),
        ]

        _, _, u_rms_rel, state = approach_distance(A, a_star, eps)
        final_xy = xy()[:, :2]
        post_dicts = [{**c.as_dict("post_run"), "note": c.note} for c in post_checks]
        return {"observables": obs,
                "extra": {"psi6": psi6, "psi6_sem": psi6_sem,
                          "nn_distance_d": bond_mean, "nn_std_rel": bond_std / bond_mean,
                          "min_sep_d": min_sep, "Gamma": float(U_star(a_star, A, eps)),
                          "coord_hist": list(map(float, coord_hist / coord_hist.sum())),
                          "state_predicted": state, "u_rms_rel_einstein": float(u_rms_rel),
                          "pe_mean": pe_mean, "pe_sem": pe_sem, "pe_rdf": pe_rdf,
                          "post_checks": post_dicts,
                          "post_checks_ok": all(c.ok for c in post_checks)},
                "arrays": {"rdf_r": rdf_r, "rdf_g": rdf_g, "final_xy": final_xy,
                          "coord_hist": coord_hist / coord_hist.sum()}}

    return RUN.Build(
        sim=sim, forces=[tab, wca], n_particles=N,
        sample=sample, pe_per_particle=pe_pp,
        n_eq=n_eq, n_prod=n_prod, sample_every=sample_every,
        tags=["2D", "soft_repulsion", "r^-3", "WCA_core", "newtonian",
             "pair_interaction", "structure"],
        physical={"N": N, "A": A, "phi": phi, "r_c_star": r_c_star, "L_star": L_star},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# ⑦ 검증 — 해석해가 없으므로 극한·일관성으로 (bd-physics §7.6)
# ════════════════════════════════════════════════════════════════════════
def energy_from_rdf(r, g, A, eps, phi, r_c_star):
    """⟨U⟩/N = (ρ/2) ∫ U(r) g(r) 2πr dr   (2D).  ρ = 4φ/(π d²) → 무차원 ρ* = 4φ/π"""
    rho_star = 4 * phi / math.pi
    m = r < r_c_star
    Ur = U_star(r[m], A, eps) - A / r_c_star**3          # 시프트된 실제 퍼텐셜
    integ = np.trapezoid(Ur * g[m] * 2 * math.pi * r[m], r[m])
    return 0.5 * rho_star * integ


# ════════════════════════════════════════════════════════════════════════
# main
# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--A", type=float, default=None, help="r⁻³ 진폭 (무차원)")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="L3 스펙을 specs/<run_id>.json 으로 쓰고 종료 (실행 안 함)")
    ap.add_argument("--dilute", action="store_true", help="희박극한 검증 런")
    ap.add_argument("--tobs", type=float, default=None, help="T_obs/τ_B 덮어쓰기")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--N", type=int, default=None, help="입자 수 덮어쓰기")
    ap.add_argument("--samples", type=int, default=400, help="생산 구간 표본 수")
    ap.add_argument("--dt-scale", type=float, default=1.0,
                    help="dt 배율. CV1(수렴 확인)용 — 0.5 로 절반")
    ap.add_argument("--rc-shells", type=float, default=5.0,
                    help="r_c = (이 값) × a_mean. CV2(컷오프 수렴)용")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake/soft-r3-2d-A-sweep/system.yaml")
    num = sys_["numerics"]

    if args.dilute:
        # 희박극한 검증: g(r) → exp(−βU). bin별 통계가 관건이라 N·T_obs를 크게 잡는다.
        # 첫 시도(N=200, 60τ_B)는 r<1.5d 에서 기대 쌍이 6개 미만이라 대조가 불가능했다.
        A, phi, N = 10.0, 0.01, 800
        r_c_star = 8.0
        T_obs_tau = args.tobs or 200.0
        tag = "dilute"
    else:
        if args.A is None:
            ap.error("--A 또는 --dilute 가 필요합니다")
        A, phi, N = args.A, sys_["phi"], sys_["N"]
        a_star = math.sqrt(math.pi / (4 * phi))
        r_c_star = args.rc_shells * a_star
        T_obs_tau = args.tobs or float(num["production_tau_B"])
        tag = f"A{A:g}"
    if args.N:
        N = args.N
    if args.smoke:
        N, T_obs_tau = min(N, 144), min(T_obs_tau, 3.0)
        tag += "-smoke"
    if args.dt_scale != 1.0:
        tag += f"-dt{args.dt_scale:g}"
    if args.rc_shells != 5.0:
        tag += f"-rc{args.rc_shells:g}"

    lg = build_ledger(sys_, A, N, phi, r_c_star, args.dt_scale, T_obs_tau)
    D = lg.derived
    tau_B, tau_int = D["tau_B"], D["tau_int"]
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")

    n_eq = int(round(float((0.2 * T_obs / dt).to(""))))   # 관측창의 20%를 평형화에
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, n_prod // args.samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, Gamma = analyze_scales(sys_, lg, A, phi, r_c_star)
    # ★ `system` 이 스펙에 들어갑니다. 예전 스펙에는 물리계가 없어서 d·η·ρ_p 를 바꿔도
    #   run_id 가 같았고(τ_B 16배 차이에도), 완료된 다른 계의 런으로 오인됐습니다
    #   (`scratch/verify_l3_spec_gaps.py` 결함 ①).
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"A": A, "phi": phi, "N": N, "r_c_star": r_c_star,
                "wca_eps": sys_["wca_eps_kT"], "Gamma": Gamma},
        numerics={"dt_star": float((dt / tau_B).to("")),
                  "dt_over_tau_int": args.dt_scale * 1e-2,
                  "n_eq": n_eq, "n_prod": n_prod, "n_samples": args.samples,
                  "sample_every": sample_every, "seed": 20260803},
        tag=tag, nhex=10)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 무결성 검사")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, A, phi, N, n_eq, n_prod)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}  A={A}   run_id={run_id}",
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
    # ★ report.txt는 execute() 이후에 쓴다 — 먼저 쓰면 prepare_outdir가 지운다.
    (outdir / "report.txt").write_text(report)

    # ── 검증 — finalize()가 이미 계산해 metrics.json에 넣어둔 것을 표로 찍는다 ──
    m = json.loads((outdir / "metrics.json").read_text())
    obs_out = m["observables"]
    res_extra = m.get("result", {})
    post_checks = res_extra.get("post_checks", [])
    psi6, psi6_sem = res_extra["psi6"], res_extra["psi6_sem"]
    bond_mean = res_extra["nn_distance_d"]
    bond_std = bond_mean * res_extra["nn_std_rel"]
    min_sep = res_extra["min_sep_d"]
    pe_mean, pe_sem = res_extra["pe_mean"], res_extra["pe_sem"]
    coord_hist = res_extra["coord_hist"]

    lines = ["", "=" * R.W, f"결과 — {sys_['label']} A={A} (Γ={Gamma:.4f})", "=" * R.W]
    lines.append(f"{'검증 (측정 vs 예측)':<26}{'측정':>14}{'예측':>14}{'차이':>10}   판정")
    for o in obs_out:
        ok = o["err_pct"] is not None and abs(o["err_pct"]) < (o["tol_pct"] or 2.0)
        lines.append(f"{o['name']:<26}{o['measured']:>14.6g}{o['predicted']:>14.6g}"
                     f"{o['err_pct']:>+9.2f}%   {'✓' if ok else '✗'}   [{o['unit']}]")
        lines.append(f"    {o['note']}")
    lines.append("")
    lines.append(f"{'사후 가드':<26}{'value':>14}{'limit':>14}{'margin':>10}   판정")
    for c in post_checks:
        lines.append(f"{c['name']:<26}{c['value']:>14.4g}{c['limit']:>14.4g}"
                     f"{c['margin']:>9.2f}×   {'✓' if c['ok'] else '✗'}")
        lines.append(f"    {c['note']}")
    all_ok = (all(abs(o["err_pct"]) < (o["tol_pct"] or 2.0) for o in obs_out
                 if o["err_pct"] is not None)
             and res_extra.get("post_checks_ok", True))

    lines += ["", "관측량 (구조)",
              f"  ⟨U⟩/N        = {pe_mean:.5f} ± {pe_sem:.5f} kT",
              f"  ψ₆ (Voronoi 가중) = {psi6:.4f} ± {psi6_sem:.4f}",
              f"  NN 거리      = {bond_mean:.4f} d   (표준편차 {bond_std:.4f} d"
              f" = {100*bond_std/bond_mean:.2f}%)",
              f"  최소 이웃거리 = {min_sep:.4f} d   (전 표본 최소값)",
              f"  Voronoi 배위수 분포: " +
              "  ".join(f"{i}:{coord_hist[i]:.3f}" for i in range(4, 10)
                        if coord_hist[i] > 0.005),
              "",
              f"  Einstein 케이지 예측 σ_bond/a_NN = {D['u_rms_rel']:.4f}"
              f"  vs 측정 {bond_std/bond_mean:.4f}"
              f"  ({100*(bond_std/bond_mean - D['u_rms_rel'])/D['u_rms_rel']:+.1f}%)",
              "    ※ 2D 결정의 절대 u_rms는 Mermin-Wagner로 로그 발산한다. 유한한 것은",
              "      상대(NN) 요동이다. Einstein 근사는 조화 + 이웃 변위 무상관을 가정하므로",
              "      비조화성(r⁻³는 바깥쪽이 무르다)만큼 측정이 더 클 수 있다. 정량 예측이 아니라",
              "      체제 판정용 지표로만 쓴다.",
              ]

    # 역변환 (물리 단위)
    d_um = float(D["d"].to("um").magnitude)
    lines += ["",
              "역변환 (물리 단위)",
              f"  a_mean = {D['a_star']*d_um:.3f} µm      L = {D['L_star']*d_um:.1f} µm",
              f"  NN 거리 = {bond_mean*d_um:.3f} µm      τ_B = {float(tau_B.to('s').magnitude):.2f} s",
              f"  dt = {float(dt.to('ms').magnitude):.4f} ms   "
              f"T_obs = {float(T_obs.to('s').magnitude):.0f} s = {T_obs_tau:.0f} τ_B",
              "=" * R.W,
              f"VERDICT: {'✓ PASS' if all_ok else '✗ FAIL'}",
              "=" * R.W]
    result = "\n".join(lines)
    print(result)
    (outdir / "result.txt").write_text(report + "\n" + result)

    make_plots(sys_, lg, A, phi, r_c_star, Gamma, outdir, args.dilute)
    print("\n".join(RID.list_artifacts(outdir, ROOT)))
    return 0 if all_ok else 1


def make_plots(sys_, lg, A, phi, r_c_star, Gamma, outdir, dilute):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    res = np.load(outdir / "observables.npz")
    fig, ax = plt.subplots(2, 2, figsize=(12, 9))
    eps = sys_["wca_eps_kT"]

    # ① g(r)  (+ 희박극한이면 exp(-βU) 대조)
    r, g = res["rdf_r"], res["rdf_g"]
    ax[0, 0].plot(r, g, "-", lw=1.4, label="측정 g(r)")
    if dilute:
        Ur = U_star(r, A, eps) - A / r_c_star**3
        ax[0, 0].plot(r, np.exp(-Ur), "--", lw=1.8, label=r"희박극한 $e^{-\beta U}$")
    ax[0, 0].axhline(1, color="k", lw=.5, alpha=.5)
    ax[0, 0].axvline(lg.derived["a_star"], ls=":", c="gray", label=r"$a_{mean}$")
    ax[0, 0].set(xlabel="r / d", ylabel="g(r)",
                 title=f"① 동경분포함수  (A={A:g}, Γ={Gamma:.3f})",
                 xlim=(0.5, min(r.max(), 6 * lg.derived["a_star"])))
    ax[0, 0].legend(); ax[0, 0].grid(alpha=.3)

    # ② 최종 배치 + Voronoi 배위수
    xy = res["final_xy"]
    L = lg.derived["L_star"]
    ax[0, 1].plot(xy[:, 0], xy[:, 1], "o", ms=max(1.5, 260 / math.sqrt(len(xy)) / 4))
    ax[0, 1].set(xlim=(-L / 2, L / 2), ylim=(-L / 2, L / 2), aspect="equal",
                 xlabel="x / d", ylabel="y / d", title="② 최종 배치")
    ax[0, 1].grid(alpha=.2)

    # ③ 평형화 + 에너지 시계열
    eq = res["eq_trace"]
    n_eq_pts = len(eq)
    ax[1, 0].plot(np.arange(n_eq_pts), eq[:, 1], "-", label="평형화")
    ax[1, 0].plot(np.linspace(n_eq_pts, n_eq_pts + 20, len(res["pe"])), res["pe"],
                  "-", lw=.8, alpha=.8, label="생산")
    ax[1, 0].set(xlabel="구간 (평형화 20구간 + 생산)", ylabel="⟨U⟩/N [kT]",
                 title="③ 퍼텐셜 에너지 — 평형 도달?")
    ax[1, 0].legend(); ax[1, 0].grid(alpha=.3)

    # ④ ψ₆ 시계열 + 배위수 분포
    ax[1, 1].plot(res["psi6"], "-", lw=.9)
    ax[1, 1].set(xlabel="표본", ylabel=r"$|\psi_6|$", ylim=(0, 1),
                 title=f"④ 육방 질서  ⟨ψ₆⟩={res['psi6'].mean():.3f}")
    ax[1, 1].grid(alpha=.3)
    a2 = ax[1, 1].twinx()
    ch = res["coord_hist"]
    a2.bar(np.arange(13), ch, alpha=.25, color="tab:orange", width=.6)
    a2.set_ylabel("Voronoi 배위수 분포", color="tab:orange")
    a2.set_xlim(-0.5, 12.5)

    for a in ax.ravel():
        a.title.set_fontsize(10)
    fig.suptitle(f"{sys_['label']}  A={A:g}  φ={phi}  N={len(xy)}", fontsize=12)
    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=140)
    plt.close(fig)


if __name__ == "__main__":
    sys.exit(main())
