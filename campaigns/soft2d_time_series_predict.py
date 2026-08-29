"""S2 — 시간분해 스윕의 예측을 **선행 런의 측정값에서 유도해서** 봉인한다.

usage: python scripts/soft2d_time_series_predict.py [--out examples/.../prediction.yaml]

## 왜 생성 스크립트인가 (손으로 쓰지 않고)

예측의 `tolerance` 는 **이론 통계오차에서 뽑아야** 한다 (`tests/conftest.py` 규칙 ①).
선행 런의 시드 SE 를 눈으로 읽어 옮기면 그 순간 손계산이 되고, 값이 바뀌었을 때
예측이 조용히 낡는다. 이 스크립트는 `runs/<선행>/metrics.json` 을 읽어
`3√2 · SE` 를 계산한다:

  · 시드 4개 평균끼리의 차이는 `SE_diff = √(SE₁² + SE₂²) ≈ √2 · SE`
  · 3σ 를 허용오차로 쓴다 → `tolerance = 3√2 · SE`

## 이 스크립트가 조심하는 것

**선행 런과 새 런은 시드가 다르다** (선행 `1–4`, 새 `5–8`). 같은 시드를 쓰면
HOOMD `Brownian` 이 비트 단위로 재현되므로 후반 창 비교가 **산술 항등식**이 된다
(`conftest` 규칙 ③: 요동하지 않는 측정값은 항등식이다). 시드를 갈라야 진짜 검정이다.

**비교 창은 `A` 마다 다르다.** 선행 런의 총 길이가 달랐다 (`A ≤ 1` 은 `30 τ_d`,
`A = 10` 은 `80 τ_d`). 후반 절반이 각각 `[20,30]` 과 `[60,80]` 이므로 그 창으로
비교한다 — 같은 창 이름을 쓰면 다른 시각을 비교하게 된다.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.analysis.structure import zahn_phase                    # noqa: E402
from simbot.build import (box_si_for_coverage,                      # noqa: E402
                          coverage_from_sigma_over_d)
from simbot.units import scales_soft2d, water_viscosity_si          # noqa: E402

PRIOR = REPO / "runs" / "2026-07-29_soft-r3-2d-A-sweep" / "metrics.json"
OUT_DEFAULT = REPO / "examples" / "soft-r3-time-resolved" / "prediction.yaml"

SIGMA_SI = 5.0e-6
D_OVER_SIGMA = 3.0
COVERAGE_MAX = 0.10
T_SI = 298.15
N_PARTICLES = 100
AMPLITUDES = (0.1, 1.0, 10.0)
NEW_SEEDS = (5, 6, 7, 8)

#  선행 런에서 후반 절반이 덮은 시각 [τ_d] — 총 길이가 A 마다 달랐다
PRIOR_WINDOW = {0.1: (20.0, 30.0), 1.0: (20.0, 30.0), 10.0: (60.0, 80.0)}
#  선행 런의 궤적 전체 최소분리 (단위 d) — 프레임 수가 늘면 더 작아진다
PRIOR_MIN_SEP = {0.1: 0.2172, 1.0: 0.4586, 10.0: 0.7164}

N_SIGMA = 3.0


def tol_from_se(se: float, n_sigma: float = N_SIGMA) -> float:
    """시드 4개 평균끼리의 차이에 대한 `n_sigma` 허용오차."""
    return float(n_sigma * np.sqrt(2.0) * se)


def main() -> int:
    out_path = Path(sys.argv[sys.argv.index("--out") + 1]) if "--out" in sys.argv \
        else OUT_DEFAULT
    prior = json.loads(PRIOR.read_text())

    geom = box_si_for_coverage(n_particles=N_PARTICLES, sigma_si=SIGMA_SI,
                               coverage_max=COVERAGE_MAX,
                               d_over_sigma_round=D_OVER_SIGMA)
    scales = scales_soft2d(d_si=geom["d_si"], sigma_si=SIGMA_SI, T_si=T_SI)
    eta_si, extrapolated = water_viscosity_si(T_SI)
    assert not extrapolated, "물 점도 표 밖이면 provenance 를 낮춰야 한다"

    items: list[dict] = []

    # --- P1/P2/P7: 후반 창이 선행 런(독립 시드)과 일치하는가 -------------------
    for A in AMPLITUDES:
        key = f"A{A:g}_random_sqbox"
        p = prior[key]
        lo, hi = PRIOR_WINDOW[A]
        for quantity, field, unit in (
                ("psi6_global", "psi6_global", ""),
                ("defect_fraction", "defect_fraction", ""),
                ("energy_per_particle", "energy_pp", "kT")):
            m, se = p[field]["mean"], p[field]["se"]
            items.append({
                "quantity": f"{quantity}__A{A:g}__t{lo:g}_{hi:g}",
                "value": float(m),
                "tolerance": f"±{tol_from_se(se):.6g}",
                "basis": (f"선행 런 runs/2026-07-29_soft-r3-2d-A-sweep 의 "
                          f"{key} 후반 절반 시드 4개 평균 {m:.6g} ± {se:.4g} (SE). "
                          f"허용오차 = 3√2·SE (시드 4 평균끼리 차이의 3σ)"),
                "discriminates": ("t=0 부터 연속 표집하는 새 방식이 '평형화 버리고 "
                                  "후반만' 방식과 같은 답을 주는가. 다르면 과도구간이 "
                                  "후반 창을 오염시켰거나 분석 경로가 갈라진 것이다"),
                "competing_value": None,
                **({"unit": unit} if unit else {}),
            })

    # --- P3: 배위수 종류 수 — 결함의 성격 ------------------------------------
    for A in AMPLITUDES:
        kinds = int(prior[f"A{A:g}_random_sqbox"]["coord_kinds"])
        items.append({
            "quantity": f"coord_kinds__A{A:g}",
            "value": kinds,
            "tolerance": "±0",
            "basis": (f"선행 런 A{A:g}_random_sqbox 의 배위수 히스토그램에서 "
                      f"분율 > 0.5 % 인 종류가 {kinds}개. 정수량이므로 오차 0 을 요구한다"),
            "discriminates": ("결함의 **성격**. 분율은 같아도 종류 수가 다르면 다른 "
                              "물리다 — 액체(6종) vs 전위만(3종). 카드 §8.2"),
            "competing_value": None,
        })

    # --- P4: 초기조건의 min_sep=0.8 d 껍질이 채워진다 -------------------------
    #  `random_2d_snapshot` 은 기각표집으로 `min_sep = 0.8 d` 를 강제한다 →
    #  t=0 에서 g(r < 0.8 d) = 0 이 **정확히** 성립한다. 평형 액체는 그렇지 않다.
    #  자유확산으로 간극 0.8→0.5 d 를 메우는 시간: MSD* = 4t* ⇒ t* ≈ 0.3²/4 = 0.0225.
    fill_t = float(0.3**2 / 4.0)
    items.append({
        "quantity": "rdf_at_0.5d__A0.1__first_window",
        "value": ">0",
        "tolerance": ">0",
        "basis": (f"초기배치는 기각표집으로 min_sep = 0.8 d 를 강제하므로 t=0 에서 "
                  f"g(r<0.8d) = 0 이 정확히 성립한다. 이것은 **초기조건 인공물**이고 "
                  f"자유확산이 메운다: MSD* = 4t* ⇒ 간극 0.3 d 는 t* ≈ {fill_t:.4g} τ_d. "
                  f"첫 창(0–20 τ_d)은 그보다 {20/fill_t:.0f}배 길다"),
        "discriminates": ("초기조건의 배제부피 껍질이 실제로 완화되는가. 남아 있으면 "
                          "표집이 초기조건에 갇힌 것이다"),
        "competing_value": None,
    })

    # --- P5: 기준 원판 겹침 — σ 를 붙인 결과로만 말할 수 있는 예측 -------------
    sigma_over_d = geom["sigma_over_d"]
    for A in AMPLITUDES:
        prior_ms = PRIOR_MIN_SEP[A]
        overlaps = prior_ms < sigma_over_d
        items.append({
            "quantity": f"min_separation_over_sigma__A{A:g}",
            "value": f"<{1.0:g}" if overlaps else f">{1.0:g}",
            "tolerance": f"<{1.0:g}" if overlaps else f">{1.0:g}",
            "basis": (f"선행 런 최소분리 {prior_ms:.4g} d = "
                      f"{prior_ms / sigma_over_d:.4g} σ (σ/d = {sigma_over_d:.6g}). "
                      f"새 런은 프레임이 2배·길이가 "
                      f"{80.0 / (PRIOR_WINDOW[A][1]):.3g}배이므로 극값이 더 작아진다 "
                      f"→ 겹침 판정은 유지되거나 강화된다"),
            "discriminates": ("'5 µm 원판 100개' 라는 그림이 이 A 에서 물리적으로 "
                              "성립하는가. A/r³ 에는 경질 코어가 없으므로 모델은 "
                              "겹침을 허용한다 — 성립하지 않는다면 그것이 결과다"),
            "competing_value": None,
            "unit": "sigma",
        })

    # --- P6: 도달 시각의 순서 — 전이점이 느리다 ------------------------------
    items.append({
        "quantity": "t90_ordering",
        "value": "t90(A=10) > t90(A=1) 이고 t90(A=10) > 1 tau_d",
        "tolerance": "정성 (순서)",
        "basis": ("A=10 은 Γ = π^{3/2}·10 = 55.68 로 Zahn 전이점 59.88 의 −7.0 % 다 "
                  "[출처, 미재현] → 임계 완화가 예상된다. A ≤ 1 은 깊은 액체이고 "
                  "무작위 초기배치가 이미 액체 배치에 가까우므로 ψ₆ 가 거의 평탄할 것이다"),
        "discriminates": "전이점 근처에서 구조 완화가 느려지는가 (임계 완화의 서명)",
        "competing_value": "세 A 의 도달 시각이 모두 O(1) tau_d 로 구별되지 않는다",
    })

    doc = {
        "items": items,
        "regimes": {
            "amplitudes": list(AMPLITUDES),
            "gamma_zahn": {f"A{A:g}": zahn_phase(A)["gamma"] for A in AMPLITUDES},
            "phase_zahn": {f"A{A:g}": zahn_phase(A)["phase_zahn"]
                           for A in AMPLITUDES},
            "sigma_si": SIGMA_SI,
            "d_si": geom["d_si"],
            "L_si": geom["L_si"],
            "L_star": geom["L_star"],
            "coverage": geom["coverage"],
            "coverage_max": COVERAGE_MAX,
            "sigma_over_d": geom["sigma_over_d"],
            "tau_d_si": scales.time_si,
            "D0_si": scales.diffusivity_si,
            "eta_si": eta_si,
            "T_si": T_SI,
            "seeds": list(NEW_SEEDS),
            "note": ("무차원 물리는 A 하나로 결정된다 (n* = 1, 경질 코어 없음). "
                     "커버리지는 τ_d 의 초 값과 점입자 이상화의 타당성만 바꾼다 — "
                     "ψ₆·g(r)·결함에 대한 감도는 정확히 0 이다"),
        },
        "alternatives": [
            "선행 런과 시드를 갈랐다 (1–4 → 5–8). 같은 시드면 HOOMD Brownian 이 "
            "비트 단위로 재현되므로 후반 창 비교가 산술 항등식이 된다.",
            "INCONCLUSIVE 예상 — '평형에 도달했다'는 이 런으로 판정할 수 없다. "
            "선행 런에서 A=10 의 ψ₆ 표류가 −0.013 ± 0.008 (1.6σ) 로 0 과 구별되지 "
            "않았고, 구별 안 됨은 '평형'이 아니라 '이 길이로는 안 보인다'는 뜻이다.",
            "INCONCLUSIVE 예상 — A=0.1 과 A=1 의 ψ₆ 차이(0.049 vs 0.063)는 유한크기 "
            "바닥(1/√N = 0.1) 아래다. 둘의 배향 질서는 '둘 다 없음' 으로만 말할 수 있다.",
            "N=100 은 A=10 의 r_cut 요구(N ≥ 252)를 만족하지 않는다 (카드 §9). "
            "βU(r_cut) = 0.09 kT 가 A=10 결과에 남아 있다 — N 수렴 검사는 이 런의 범위 밖이다.",
            "A=0.1 에서 기준 원판이 겹치는 것은 코드 오류가 아니라 모델의 성질이다. "
            "겹침을 없애려면 커버리지를 3.71 % 아래로 내려야 하고, 그러면 무차원 결과는 "
            "그대로이고 τ_d 만 바뀐다.",
        ],
    }
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "# S2 PREDICTION — 시뮬레이션 전에 답을 적는다  (자동 생성: "
        "scripts/soft2d_time_series_predict.py)\n"
        "#\n"
        "# ⚠ 이 파일은 S5 실행 전에 봉인된다 (SEALED.sha256). 실행 후 수정 금지.\n"
        "# 허용오차는 선행 런의 시드 SE 에서 유도했다 (3√2·SE) — 관측값에 맞춰 재단하지 않았다.\n"
        "#\n"
        + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8")
    print(f"→ {out_path.relative_to(REPO)}  ({len(items)} 항목)")
    for it in items:
        print(f"   {it['quantity']:<46} {str(it['value']):>12}  {it['tolerance']}")
    print(f"\n  d = {geom['d_si']*1e6:.3f} µm · L = {geom['L_si']*1e6:.2f} µm · "
          f"coverage = {geom['coverage']:.4%} · τ_d = {scales.time_si:.1f} s "
          f"= {scales.time_si/60:.2f} min")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
