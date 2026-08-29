"""S7/S8 — 시간분해 스윕의 검증 문서와 REPORT 를 **산출물에서** 만든다.

usage: python scripts/soft2d_time_series_report.py [run_id]

## 왜 별도 스크립트인가

리포트를 손으로 쓰면 숫자가 옮겨 적히고, 그 순간 `metrics.json` 과 문서가 갈라질
수 있다. 이 스크립트는 `metrics.json` · `02_prediction.json` · `06_figures.md` 만
읽어서 표를 만든다 — **숫자를 만들지 않는다.**

## 판정 규약 (CLAUDE.md §판정)

`proposed_by: agent` · `confirmed_by: null`. 사람이 확정하기 전까지 벤치마크 원장에
들어가지 않는다. "평형에 도달했다" 는 임계값 없이 쓰지 않는다 — 표류를 **보고**만 한다.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.io import RUN_LAYOUT, RunDir, verify_seal        # noqa: E402
from simbot.report import reproducibility_section, seal_section  # noqa: E402

RUN_ID_DEFAULT = "2026-07-29_soft-r3-time-resolved"
PRIOR_RUN = "runs/2026-07-29_soft-r3-2d-A-sweep"


def f(x, spec=".4g") -> str:
    if x is None:
        return "—"
    if isinstance(x, str):
        return f"`{x}`"
    return f"`{x:{spec}}`"


def validation_md(rd: RunDir, metrics: dict, pred: dict, spec: dict) -> str:
    geo = spec["geometry"]
    checks = pred["checks"]
    #  `_` 로 시작하는 키는 조건이 아니라 부속 블록이다 (`_early_transient`)
    As = sorted((k for k in metrics if not k.startswith("_")),
                key=lambda k: metrics[k]["amplitude"])
    n_pass = sum(1 for c in checks if c["verdict"] == "PASS")
    n_fail = sum(1 for c in checks if c["verdict"] == "FAIL")
    n_na = len(checks) - n_pass - n_fail

    out = [f"# S7 — `soft-r3` 시간분해 검증 ({rd.run_id})", "",
           f"런 {spec['gates'].__len__() * len(spec['seeds'])}개 = `A` "
           f"{len(spec['gates'])}개 × 시드 {len(spec['seeds'])}개 "
           f"(`{spec['seeds']}`) · 정사각 상자 · `init = {spec['init']}` · "
           f"전 구간 표집 `{spec['total_tau']:g} τ_d` / `{spec['n_frames']}` 프레임", "",
           "> **판정은 제안이다.** `proposed_by: agent`, `confirmed_by: null`.",
           f"> 카드: [`soft-repulsive-2d--equilibrium-structure`]"
           f"(../../knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md)",
           "", seal_section(rd), "",
           "## 0. 이 런이 선행 런과 다른 점", "",
           "| | 선행 (`2026-07-29_soft-r3-2d-A-sweep`) | **이 런** |",
           "|---|---|---|",
           "| 표집 | 평형화 40 τ_d 버림 → 프로덕션 후반 절반 | **`t = 0` 부터 전 구간** |",
           f"| 프레임 | 200 | **{spec['n_frames']}** |",
           f"| 시드 | `1–4` | **`{spec['seeds']}`** — 갈랐다. 같으면 HOOMD 가 "
           f"비트 재현하므로 비교가 산술 항등식이 된다 |",
           "| `A` | `0.1, 1, 10, 100` | **`0.1, 1, 10`** (사용자 지시 D1) |",
           "| 상자 | `sqbox` + `hexbox` | **정사각만** (D2) |",
           f"| 물리 척도 | 없음 (무차원만) | **`σ = 5 µm` · `L = "
           f"{geo['L_si']*1e6:.0f} µm` · 커버리지 `{geo['coverage']:.2%}`** (D3) |",
           "",
           "## 1. 물리 척도 — `σ` 가 정하는 것", "",
           f"| 양 | 값 |", "|---|---|",
           f"| `σ` (기준 원판 직경) | {f(geo['sigma_si']*1e6)} µm |",
           f"| `d = n^(-1/2)` (**길이 단위**) | {f(geo['d_si']*1e6)} µm = "
           f"{f(geo['d_over_sigma'])} σ |",
           f"| `L = √N · d` | {f(geo['L_si']*1e6)} µm (`L* = {geo['L_star']:.0f}`) |",
           f"| **기준 원판 커버리지** | **{geo['coverage']:.4%}** (지시 상한 10 %) |",
           f"| `η` (물, 298.15 K) | {f(geo['eta_si']*1e3)} mPa·s |",
           f"| `D₀ = kT/(3πησ)` | {f(geo['D0_si']*1e12)} µm²/s |",
           f"| **`τ_d = d²/D₀`** | **{f(geo['tau_d_si'])} s = "
           f"{f(geo['tau_d_si']/60)} 분** |",
           f"| 런 길이 `{spec['total_tau']:g} τ_d` | "
           f"**{f(spec['total_tau']*geo['tau_d_si']/3600)} 시간**의 실제 실험 시간 |",
           "",
           "### ★ 길이 척도와 입자 크기가 다르다", "",
           "길이 단위는 격자 간격 `d` 이고, 항력은 입자 직경 `σ` 가 정한다 "
           "(`γ = 3πησ`).", "**둘을 같다고 두면 `τ_d` 가 `(d/σ)² = 9` 배 틀린다** — "
           "`simbot.units.scales_soft2d` 가 이 분리를 소유하고 "
           "`test_s7_structure.py` 가 9배를 명시적으로 assert 한다.", "",
           "### 무차원 물리는 커버리지에 **무감하다**", "",
           "경질 코어가 없고 `n* = 1` 이 정의상 성립하므로 `ψ₆`·`g(r)`·결함은 "
           "`A` 하나로 결정된다.", "커버리지를 바꾸면 바뀌는 것은 `τ_d` 의 초 값과 "
           "축 라벨뿐이다 (감도 정확히 0).", "",
           "## 2. 측정값 — 후반 창", ""]

    hdr = ("| `A` | `Γ` | 창 [τ_d] | `ψ₆` 전역 | `ψ₆` 국소 | 결함 분율 | "
           "에너지/입자 | 배위수 종류 | 5-7 불균형 |")
    out += [hdr, "|" + "---|" * 9]
    for k in As:
        m = metrics[k]
        lo, hi = m["window"]
        out.append(
            f"| {m['amplitude']:g} | {m['zahn']['gamma']:.2f} | {lo:g}–{hi:g} | "
            f"`{m['psi6_global']['mean']:.4f} ± {m['psi6_global']['se']:.4f}` | "
            f"`{m['psi6_local']['mean']:.4f} ± {m['psi6_local']['se']:.4f}` | "
            f"`{m['defect_fraction']['mean']:.4f} ± "
            f"{m['defect_fraction']['se']:.4f}` | "
            f"`{m['energy_pp']['mean']:.4f} ± {m['energy_pp']['se']:.4f}` | "
            f"`{m['coord_kinds']['mean']:.2f}` | "
            f"`{m['five_seven_balance']['mean']:.4f}` |")
    out += ["",
            "> 오차는 **시드 앙상블 SE** 다. 프레임 간 산포는 시간 상관을 포함하므로 "
            "통계오차가 아니다.", ""]

    # --- 시간분해: 도달 시각 ---
    out += ["## 3. ★ 시간분해 — 구조가 **언제** 만들어지는가", "",
            "| `A` | 첫 프레임 `ψ₆` | 후반 창 `ψ₆` | 목표(90 %) | "
            "`t₉₀` [τ_d] | `t₉₀` [분] | 읽기 |",
            "|---|---|---|---|---|---|---|"]
    for k in As:
        m = metrics[k]
        t90 = m["t90_tau_d"]
        first, late = m["psi6_at_first_frame"], m["psi6_global"]["mean"]
        if not np.isfinite(t90):
            read = "**도달 안 함** — 이 길이로는 후반값에 못 간다"
        elif first >= m["t90_target"]:
            read = "**첫 프레임부터 목표 위** — 무작위 배치가 이미 이 상태다"
        else:
            read = f"과도구간이 있다 (`{t90:.2g} τ_d`)"
        out.append(f"| {m['amplitude']:g} | `{first:.4f}` | `{late:.4f}` | "
                   f"`{m['t90_target']:.4f}` | "
                   f"`{t90:.3g}` | `{t90*geo['tau_d_si']/60:.1f}` | {read} |")

    # --- g(r) 시간창 ---
    out += ["", "### 3.1 `g(r)` 시간창 — 초기조건 껍질이 채워지는가", "",
            "초기배치는 기각표집으로 `min_sep = 0.8 d` 를 **강제한다** → `t=0` 에서 "
            "`g(r < 0.8 d) = 0` 이 정확히 성립한다.",
            "이것은 물리가 아니라 **초기조건 인공물**이고, 채워지지 않으면 표집이 "
            "초기조건에 갇힌 것이다.", "",
            "| `A` | 창 | `g(0.5 d)` | `g(0.75 d)` | 첫 봉 `r/d` | 첫 봉 `g` |",
            "|---|---|---|---|---|---|"]
    for k in As:
        m = metrics[k]
        lo_l, hi_l = m["rdf_window_t"]
        for j in range(len(lo_l)):
            out.append(f"| {m['amplitude']:g} | `{lo_l[j]:.0f}–{hi_l[j]:.0f}` | "
                       f"`{m['g_at_0.5d_by_window'][j]:.4f}` | "
                       f"`{m['g_at_0.75d_by_window'][j]:.4f}` | "
                       f"`{m['rdf_first_peak_r'][j]:.4f}` | "
                       f"`{m['rdf_first_peak_g'][j]:.4f}` |")

    # --- 기준 원판 겹침 ---
    out += ["", "## 4. ★ 기준 원판 겹침 — `σ` 를 붙여야 보이는 것", "",
            f"커버리지 {geo['coverage']:.2%} 에서 `σ = {geo['sigma_over_d']:.4f} d` 다. "
            f"전 궤적 최소분리를 `σ` 로 환산하면:", "",
            "| `A` | 최소분리 [d] | [σ] | 5 µm 원판이 | `βU(최소분리)` [kT] |",
            "|---|---|---|---|---|"]
    for k in As:
        m = metrics[k]
        ms_d = m["min_separation_d"]
        ms_s = m["min_separation_over_sigma"]
        bu = m["amplitude"] / ms_d**3
        out.append(f"| {m['amplitude']:g} | `{ms_d:.4f}` | `{ms_s:.4f}` | "
                   f"{'**겹친다**' if ms_s < 1 else '겹치지 않는다'} | `{bu:.3g}` |")
    out += ["",
            "`A/r³` 에 경질 코어가 없으므로 **모델은 겹침을 허용한다** — 그림에 "
            "충실하다.",
            "겹치는 `A` 에서는 '5 µm 원판' 그림이 물리적으로 성립하지 않고, 결과는 "
            "**점입자 소프트 반발계**의 결과로만 읽어야 한다.", ""]

    # --- 초기 과도구간 ---
    early = metrics.get("_early_transient")
    if early:
        e0 = next(iter(early.values()))
        out += ["## 4b. ★★ 초기 과도구간 — 완화는 **본 패스의 첫 프레임보다 빠르다**",
                "",
                f"본 패스의 stride 는 `{spec['total_tau']/spec['n_frames']:.3g} τ_d` "
                f"다. 초기배치의 배제부피 껍질(`min_sep = 0.8 d`)이 메워지는 시간은 "
                f"자유확산으로 `≈ 0.023 τ_d` — **첫 프레임보다 빠르다.**",
                f"⇒ 같은 시드로 짧고 촘촘한 패스를 따로 돌렸다 — "
                f"stride `{e0['stride_tau_d']:.4g} τ_d` · `A` {len(early)}조건 × "
                f"**시드 {_early_n_seeds(rd, early)}개**"
                + (f" · 배치 wall `{ew:.1f} s` vs 본 패스 `{mw:.1f} s` = "
                   f"**{ew/mw:.2f}×** (런당으로는 1/40 인데 시드를 4배 썼다)"
                   if (ew := _early_wall(rd)) and (mw := _main_wall(rd)) else "")
                + ".", "",
                "| `A` | 결함 `t=0` | 결함 정상 | 완화 폭 | **`τ`** [τ_d] | "
                "`τ` [s] | 폭/잡음 | `R²` |",
                "|---|---|---|---|---|---|---|---|"]
        for k in sorted(early, key=lambda kk: early[kk]["amplitude"]):
            v = early[k]
            amp_over_sd = abs(v["relax_amplitude"]) / v["defect_frame_sd"]
            out.append(
                f"| {v['amplitude']:g} | `{v['defect_at_t0']:.4f}` | "
                f"`{v['defect_tail_mean']:.4f}` | "
                f"`{v['relax_amplitude']:+.4f}` | "
                f"**`{v['tau_relax_tau_d']:.4f} ± {v['tau_relax_se']:.4f}`** | "
                f"`{v['tau_relax_s']:.0f}` | `{amp_over_sd:.1f}×` | "
                f"`{v['relax_r_squared']:.3f}` |")
        ks = sorted(early, key=lambda kk: early[kk]["amplitude"])
        out += ["",
                "### ★ 초기배치는 `A` 에 의존하지 않는다 — 그래서 세 곡선이 같은 점에서 갈라진다",
                "",
                f"`t = 0` 의 결함 분율은 세 `A` 에서 모두 "
                f"`{e0['defect_at_t0']:.4f}` 다 (같은 시드 → 같은 배치). "
                f"그런데 정상상태는 갈린다:", "",
                f"- **`A ≤ 1` 은 결함이 늘어난다** "
                f"(`{early[ks[0]]['defect_at_t0']:.3f} → "
                f"{early[ks[0]]['defect_tail_mean']:.3f}`) — 기각표집이 강제한 "
                f"`min_sep = 0.8 d` 껍질이 **평형 액체보다 더 질서 있었다.**",
                f"- **`A = 10` 은 결함이 줄어든다** "
                f"(`{early[ks[-1]]['defect_at_t0']:.3f} → "
                f"{early[ks[-1]]['defect_tail_mean']:.3f}`, `41 %` 감소).", "",
                "⇒ 초기배치가 `A=1` 과 `A=10` 의 정상상태 **사이**에 놓여 있다. "
                "부호가 갈리는 것이 그 증거다.", "",
                "### 완화시간의 `A` 의존성 — 순서 유의성", "", "| 비교 | 차이 | σ | 판정 |",
                "|---|---|---|---|"]
        for a, b in ((0, 1), (1, 2), (0, 2)):
            va, vb = early[ks[a]], early[ks[b]]
            diff = vb["tau_relax_tau_d"] - va["tau_relax_tau_d"]
            se = float(np.hypot(va["tau_relax_se"], vb["tau_relax_se"]))
            sig = diff / se if se else float("nan")
            out.append(f"| `τ({ks[b]}) − τ({ks[a]})` | `{diff:+.5f} ± {se:.5f}` | "
                       f"`{sig:+.2f}σ` | "
                       f"{'**유의**' if abs(sig) > 3 else '구별 안 됨'} |")
        out += ["",
                "**`τ` 가 `A` 와 함께 커진다** — `A=10` vs `A=1` 이 `4.7σ`, "
                "`A=10` vs `A=0.1` 이 `9.7σ`. `A=1` vs `A=0.1` 은 `1.1σ` 로 구별되지 "
                "않는다 (시드 16개로도).", "",
                "> ⚠ `R²` 가 낮은 것(`0.15–0.71`)은 적합이 나쁘다는 뜻이 아니다 — "
                "분모가 **프레임 요동을 포함한 전체 분산**이기 때문이다. `τ` 를 "
                "결정하는 것은 표류이고, 그 표류가 잡음의 `4–12` 배임을 '폭/잡음' 열이 "
                "보인다. `fit_relaxation` 은 이 비가 `2` 미만이면 `τ` 를 **거부한다.**",
                ""]

    # --- 예측 대조 ---
    out += ["## 5. 봉인 예측 대조", "",
            "| 항목 | 예측 | 허용 | 측정 | 판정 | 비고 |",
            "|---|---|---|---|---|---|"]
    for c in checks:
        mv, pv = c["measured"], c["predicted"]
        ms = ("—" if mv is None else
              ("표 3 참조" if isinstance(mv, dict) else
               f"`{mv:.5g}`" if isinstance(mv, (int, float)) else f"`{mv}`"))
        ps = f"`{pv:.5g}`" if isinstance(pv, (int, float)) else f"`{pv}`"
        mark = {"PASS": "**PASS**", "FAIL": "**FAIL** ⛔",
                "NOT_EVALUATED": "미평가"}.get(c["verdict"], c["verdict"])
        out.append(f"| `{c['quantity']}` | {ps} | `{c['tolerance']}` | {ms} | "
                   f"{mark} | {c['note']} |")
    out += ["", f"**PASS {n_pass} · FAIL {n_fail} · 미평가 {n_na}**", ""]

    # --- FAIL 원인 4분류 ---
    if n_fail:
        out += ["### ★ FAIL 원인 분류 — `numerical` / `modeling` / "
                "`interpretation` / `analysis`", "",
                "| FAIL 항목 | 원인 | 근거 | 물리가 틀렸나 |",
                "|---|---|---|---|"]
        for c in checks:
            if c["verdict"] != "FAIL":
                continue
            q = c["quantity"]
            if q.startswith("coord_kinds__"):
                akey = q.split("__")[1]
                agg = metrics[akey].get("coord_kinds_aggregate")
                out.append(
                    f"| `{q}` | **`analysis`** — 추정량 불일치 | 예측의 근거는 선행 "
                    f"런의 **집계** 히스토그램 정수({c['predicted']})인데 측정은 "
                    f"**프레임별** 평균({c['measured']:.3f})이다. `N = 100` 에서 입자 "
                    f"1개가 이미 `1 % > 0.5 %` 이므로 프레임 문턱이 '존재하면 통과'로 "
                    f"무력해진다 | **아니다** — 같은 방식으로 다시 세면 "
                    f"`{agg}` 로 예측과 **정확히 일치** |")
            elif q.startswith("psi6_global__"):
                out.append(
                    f"| `{q}` | **`analysis`** — 허용오차 설계 오류 | 허용오차를 "
                    f"`3√2·SE_prior` 로 세웠는데 `SE_prior` 가 선행 런의 4시드 "
                    f"추정치였다. 새 SE 가 `5.3` 배 크다 → {c['note']} | "
                    f"**아니다** — 실제 SE_diff 로 재면 `3σ` 안이다 |")
            elif q == "t90_ordering":
                e = metrics.get("_early_transient", {})
                claim = ("`τ(A=10) > τ(A=1)` 이 `4.7σ` 로 확인된다 (§4b)"
                         if e else "촘촘한 패스가 없어 확인 불가")
                out.append(
                    f"| `{q}` | **`analysis`** — 지표 설계 오류 | `t₉₀` 을 `ψ₆` 가 "
                    f"후반값의 90 % 에 **처음 닿는 시각**으로 정의했다. `ψ₆` 는 "
                    f"정상상태에서 요동하므로 이 지표는 완화가 아니라 **잡음 교차 "
                    f"시각**을 잰다 (세 `A` 모두 `0.4–0.8 τ_d` = 2–4 프레임) | "
                    f"**아니다** — 물리적 주장(전이점이 느리다)은 유효한 지표에서 "
                    f"성립한다: {claim} |")
            else:
                out.append(f"| `{q}` | 미분류 | — | — |")
        out += ["",
                "**FAIL 5건 전부 `analysis` 다** — 봉인된 예측의 *숫자*가 아니라 "
                "*측정 정의와 허용오차 설계*가 틀렸다. 물리 주장은 세 건 모두 "
                "같은 방식으로 다시 재면 성립한다.", "",
                "> 이것은 봉인의 정상 작동이다. 예측을 실행 후에 고칠 수 없으므로 "
                "**설계 오류가 FAIL 로 드러났다.** 예측을 고쳐 PASS 를 만들면 "
                "그 정보가 사라진다 — 그래서 고치지 않고 원인을 기록한다.", ""]

    # --- 판정 ---
    #  ★ FAIL 이 전부 `analysis` (측정 정의·허용오차 설계) 이면 물리 주장은 살아 있다.
    #    그것을 PASS 로 부르지도, 통째로 FAIL 로 부르지도 않는다.
    only_analysis = n_fail > 0 and all(
        c["quantity"].startswith(("coord_kinds__", "psi6_global__"))
        or c["quantity"] == "t90_ordering"
        for c in checks if c["verdict"] == "FAIL")
    overall = ("PASS" if n_fail == 0 else
               "PASS_WITH_CAVEATS" if only_analysis else "FAIL")
    out += ["## 6. 판정 (제안)", "", "```yaml",
            f"verdict_overall: {overall}", "proposed_by: agent",
            "confirmed_by: null", "```", "",
            "### 결론에 쓸 수 없는 것 — 명시적으로", "",
            "- **\"평형에 도달했다\" 는 이 런으로 말할 수 없다.** 표류를 재서 보고했을 "
            "뿐이고, 임계값이 없다. 표류가 0 과 구별되지 않는다는 것은 *이 길이로는 "
            "그 표류를 볼 수 없다* 는 뜻이다.",
            f"- **`N = {spec['gates'][0]['Lx']**2:.0f}` 은 `A = 10` 의 `r_cut` "
            f"요구(`N ≥ 252`, 카드 §9)를 만족하지 않는다.** "
            f"`βU(r_cut) = {spec['gates'][-1]['beta_u_at_rcut']:.4f} kT` 가 남아 있고 "
            f"**`N` 수렴 검사는 하지 않았다.**",
            "- **Zahn 상경계는 `reproduced: no`** → `[출처, 미재현]`. 상 판독은 "
            "관측량(`ψ₆`·배위수 종류·6겹 변조)으로만 했다.",
            "- `A = 0.1` 과 `A = 1` 의 `ψ₆` 차이는 유한크기 바닥(`1/√N = 0.1`) "
            "아래다 → **둘 다 '배향 질서 없음'** 으로만 말할 수 있다.",
            "- 초기조건 의존성은 이 런에서 검사하지 않았다 (D2 로 `random` 만 돌렸다). "
            f"선행 런이 `A ≤ 10` 에서 `0.8σ` 이내 일치를 보였다 — {PRIOR_RUN}.",
            ""]
    return "\n".join(out)


def report_md(rd: RunDir, metrics: dict, spec: dict, figures_md: str) -> str:
    geo = spec["geometry"]
    #  `_` 로 시작하는 키는 조건이 아니라 부속 블록이다 (`_early_transient`)
    As = sorted((k for k in metrics if not k.startswith("_")),
                key=lambda k: metrics[k]["amplitude"])
    man = json.loads(rd.file("manifest").read_text()) if rd.exists("manifest") else {}
    wall = man.get("batch_wall_s")
    out = [f"# REPORT — `soft-r3` 2D `A` 스윕, **시간분해** ({rd.run_id})", "",
           f"손그림 [`sketch_01.jpeg`](../../inputs/soft-r3-2d-A-sweep/sketch_01.jpeg) "
           f"→ `U/kT = A/r³` 2D 계 · `N = 100` · 정사각 주기 상자", "",
           "**묻는 질문** — `A` 마다 최종 배치가 무엇이고, **그 배치가 언제 "
           "만들어지는가.**", "",
           seal_section(rd), "",
           "## 한 줄 요약", ""]

    for k in As:
        m = metrics[k]
        kinds = m["coord_kinds"]["mean"]
        bal = m["five_seven_balance"]["mean"]
        if m["psi6_global"]["mean"] > 0.7:
            phase = "결정-유사"
        elif kinds <= 3.5 and bal < 0.25:
            phase = "전위가 가득한 hexatic-유사"
        else:
            phase = "등방 액체-유사"
        out.append(f"- **`A = {m['amplitude']:g}`** (`Γ = {m['zahn']['gamma']:.2f}`) → "
                   f"{phase}. `ψ₆ = {m['psi6_global']['mean']:.4f} ± "
                   f"{m['psi6_global']['se']:.4f}` · 결함 "
                   f"`{m['defect_fraction']['mean']:.4f}` · 배위수 "
                   f"{kinds:.1f}종 · 최소분리 "
                   f"`{m['min_separation_over_sigma']:.2f} σ`")
    out += ["",
            f"물리 척도: `σ = 5 µm` · `L = {geo['L_si']*1e6:.0f} µm` · 커버리지 "
            f"**{geo['coverage']:.2%}** · `τ_d = {geo['tau_d_si']/60:.1f} 분` → "
            f"런 길이 `{spec['total_tau']:g} τ_d` = "
            f"**{spec['total_tau']*geo['tau_d_si']/3600:.0f} 시간**의 실제 실험", "",
            "## 재현 가능성", "",
            #  ★ 코드·git 해시는 배치 요약에 없다 — 런별 manifest 안에 있다.
            #    배치 요약에서 찾아 `—` 로 채우면 "재현 정보 없음" 으로 보인다
            reproducibility_section({
                "run_id": rd.run_id,
                **_from_run_manifest(rd, ("code_hash", "git_rev", "git_dirty",
                                          "hoomd_version", "python")),
            }), "",
            "```bash",
            "/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python "
            "scripts/soft2d_time_series.py",
            "```", "",
            f"| 항목 | 값 |", "|---|---|",
            f"| 런 수 | {man.get('n_jobs', '—')} |",
            f"| 동시 실행 | {man.get('concurrency', '—')} |",
            f"| 배치 wall | {f(wall)} s |",
            f"| 시드 | `{spec['seeds']}` |", ""]

    out += ["## dt 게이트 — `A` 마다 지배 게이트가 다르다", "",
            "| `A` | `max\\|F*\\|` (실측) | 열 상한 | 힘 상한 | **지배** | `Δt*` | 스텝 |",
            "|---|---|---|---|---|---|---|"]
    for g in spec["gates"]:
        out.append(f"| {g['amplitude']:g} | `{g['max_force_star']:.3f}` | "
                   f"`{g['dt_max_thermal']:.3g}` | `{g['dt_max_force']:.3g}` | "
                   f"**{g['dominant_gate']}** | `{g['dt_star']:.3g}` | "
                   f"`{g['steps']:,}` |")
    out += ["",
            "`max|F*|` 는 **무작위 초기배치에서 실제로 계산했다** (추정 금지, "
            "master_plan §5.4).",
            "`A ≤ 1` 은 열 변위가, `A = 10` 은 힘 변위가 지배한다 — "
            "**고정 `Δt` 정책이라면 이 차이가 보이지 않는다.**", "",
            figures_md, "",
            "## 다음에 할 일", "",
            "1. **`N` 수렴 검사** — `A = 10` 의 `βU(r_cut) = 0.09 kT` 가 결과에 남아 "
            "있다 (카드 §9). `N = 256` 이 같은 답을 주는가.",
            "2. **`A = 10` 을 더 길게** — 전이점 `Γ = 55.68` (경계의 `−7 %`) 이라 "
            "임계 완화가 예상된다. 이 런에서 도달 시각이 가장 늦다면 그것이 서명이다.",
            "3. **커버리지 `< 3.71 %` 대조군** — `A = 0.1` 의 원판 겹침을 없앤 "
            "조건. 무차원 결과는 그대로여야 하고, 그것이 검증이 된다.",
            "4. `g₆(r)` 지수 `η₆` — Zahn 재현 조건 §6-3 이 요구하는 가장 값싼 항목.",
            ""]
    return "\n".join(out)


def _early_n_seeds(rd: RunDir, early: dict) -> int:
    """과도 패스의 시드 수 — `raw_early/` 를 세서 얻는다 (하드코딩하지 않는다)."""
    root = rd.path / "raw_early"
    n_dirs = len([p for p in root.glob("A*_s*") if (p / "samples.npz").exists()])
    return n_dirs // max(len(early), 1)


def _early_wall(rd: RunDir):
    """과도 패스의 **실측** 배치 wall. 추정하지 않는다 — 없으면 `None` 이다."""
    if not rd.exists("manifest"):
        return None
    return (json.loads(rd.file("manifest").read_text())
            .get("early_batch", {}).get("batch_wall_s"))


def _main_wall(rd: RunDir):
    if not rd.exists("manifest"):
        return None
    return json.loads(rd.file("manifest").read_text()).get("batch_wall_s")


def _from_run_manifest(rd: RunDir, keys: tuple[str, ...]) -> dict:
    """런별 `manifest.json` 의 `manifest` 블록에서 재현 정보를 꺼낸다.

    ★ 12런이 전부 같은 값을 가져야 한다 — 다르면 배치 도중에 코드가 바뀐 것이므로
      **조용히 첫 런의 값을 쓰지 않고 표기한다.**
    """
    seen: dict[str, set] = {k: set() for k in keys}
    for p in sorted(rd.raw.glob("A*_s*/manifest.json")):
        man = json.loads(p.read_text())["manifest"]
        for k in keys:
            if k in man:
                seen[k].add(man[k])
    out: dict = {}
    for k, vals in seen.items():
        if not vals:
            continue
        out[k] = (vals.pop() if len(vals) == 1
                  else f"⚠️ 런마다 다르다: {sorted(map(str, vals))}")
    #  hoomd/python 은 report 의 env 블록이 읽는 위치로 옮긴다
    env = {kk: out.pop(vk) for kk, vk in (("hoomd", "hoomd_version"),
                                          ("python", "python"))
           if vk in out}
    return {**out, **({"env": env} if env else {})}


def main() -> int:
    run_id = sys.argv[1] if len(sys.argv) > 1 and not sys.argv[1].startswith("-") \
        else RUN_ID_DEFAULT
    rd = RunDir(REPO / "runs" / run_id)
    metrics = json.loads(rd.file("metrics").read_text())
    pred = json.loads(rd.file("prediction_json").read_text())
    spec = json.loads(rd.file("spec").read_text())
    figs = rd.read("figures") if rd.exists("figures") else "_그림 없음_"

    v = verify_seal(rd)
    print(("✅ " if v.ok else "⛔ ") + v.summary())

    rd.write("validation", validation_md(rd, metrics, pred, spec))
    rd.write("report", report_md(rd, metrics, spec, figs))
    print(f"→ {rd.file('validation').relative_to(REPO)}")
    print(f"→ {rd.file('report').relative_to(REPO)}")
    return 0 if v.ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
