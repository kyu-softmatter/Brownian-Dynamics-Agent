"""S2 — `N` 수렴 런의 예측을 **`N=100` 실측에서 유도해서** 봉인한다.

usage: python scripts/soft2d_nconv_predict.py

## 이 런의 1급 질문은 절단오차가 아니다

`N=100 → 256` 은 `βU(r_cut)` 을 `0.0904 → 0.0216 kT` 로 내린다 (카드 §9 의 미결).
그런데 **더 값있는 것은 `ψ₆` 의 유한크기 지수다.**

`|⟨ψ₆⟩| ~ N^{-p}`, `η₆ = 4p`:
  · 액체 (`g₆` 지수 감쇠)      → `p = 1/2`
  · hexatic (`g₆ ~ r^{-η₆}`)   → `p = η₆/4 ≤ 1/16` (KTHNY 경계 `η₆ = 1/4`)
  · 결정                        → `p → 0`

`A=10` 은 `ψ₆ = 0.248` 로 "결정은 아니고 전위가 가득한 hexatic-유사" 로 판독됐는데
(카드 §8.3), **그 판독을 지수가 검정한다.** `p ≈ 0.5` 면 그 `0.248` 은 유한크기
바닥일 뿐이고, `p` 가 작으면 진짜 준장거리 배향 질서다.

⇒ 이것이 Zahn 재현 조건 §6-3 이 요구하는 `η₆` 를 **`g₆(r)` 적합 없이** 얻는 경로다.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import yaml

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

from simbot.analysis.structure import (KTHNY_ETA6_HEXATIC_LIQUID,  # noqa: E402
                                       LIQUID_EXPONENT_P, zahn_phase)
from simbot.build import box_si_for_coverage                        # noqa: E402
from simbot.units import scales_soft2d                              # noqa: E402

PARENT = REPO / "runs" / "2026-07-29_soft-r3-time-resolved" / "metrics.json"
OUT = REPO / "examples" / "soft-r3-nconv" / "prediction.yaml"

N_REF, N_NEW = 100, 256
AMPLITUDES = (0.1, 1.0, 10.0)
SEEDS = (21, 22, 23, 24)          # ★ 부모 런(5–8)·완화 스윕(5–1540) 과 겹치지 않게
SIGMA_SI, D_OVER_SIGMA, T_SI = 5.0e-6, 3.0, 298.15

#  N=256 에서 실측한 게이트 (scripts/soft2d_nconv.py 가 같은 값을 재계산해 대조한다)
BETA_U_RCUT = {0.1: 0.0002, 1.0: 0.0022, 10.0: 0.0216}
N_SIGMA = 3.0


def main() -> int:
    parent = json.loads(PARENT.read_text())
    ratio = N_NEW / N_REF
    ln_ratio = float(np.log(ratio))
    items: list[dict] = []

    for A in AMPLITUDES:
        ref = parent[f"A{A:g}"]
        p6, p6se = ref["psi6_global"]["mean"], ref["psi6_global"]["se"]

        # --- ① ψ₆ 유한크기 예측: 두 극단을 모두 적는다 ---
        liquid_pred = p6 * ratio ** (-LIQUID_EXPONENT_P)
        hexatic_pred = p6 * ratio ** (-KTHNY_ETA6_HEXATIC_LIQUID / 4.0)
        #  허용오차: 부모 SE 를 스케일링한 뒤 3√2 (새 SE 가 비슷하다고 가정)
        tol = 3.0 * np.sqrt(2.0) * p6se * ratio ** (-LIQUID_EXPONENT_P)
        items.append({
            "quantity": f"psi6_global__A{A:g}__N256",
            "value": float(liquid_pred),
            "tolerance": f"±{tol:.6g}",
            "basis": (
                f"액체 가설 `|⟨ψ₆⟩| ~ N^-1/2`: {p6:.6g} × ({ratio:g})^-0.5 = "
                f"{liquid_pred:.6g}. **경쟁 가설(hexatic, η₆=1/4 → p=1/16)은 "
                f"{hexatic_pred:.6g}** 이고 둘의 간격이 "
                f"{abs(hexatic_pred - liquid_pred):.6g} = 허용오차의 "
                f"{abs(hexatic_pred - liquid_pred) / tol:.1f}배다 — 판별 가능하다"),
            "discriminates": (
                "A=%g 의 ψ₆ 가 유한크기 바닥인가 진짜 준장거리 질서인가" % A),
            "competing_value": float(hexatic_pred),
        })

        # --- ② 유한크기 지수 자체 ---
        items.append({
            "quantity": f"psi6_exponent_p__A{A:g}",
            "value": LIQUID_EXPONENT_P,
            "tolerance": "±0.15",
            "basis": (
                f"`p = -dln|⟨ψ₆⟩|/dlnN`. 액체면 0.5, KTHNY hexatic 경계면 "
                f"{KTHNY_ETA6_HEXATIC_LIQUID / 4.0:.4f}. 허용오차 ±0.15 는 두 점 "
                f"추정의 오차 전파(부모 SE {p6se:.4g} 상대 "
                f"{p6se / p6:.1%}, /ln{ratio:g}={ln_ratio:.4f}) 규모다"),
            "discriminates": f"A={A:g} 의 상 판독 (η₆ = 4p)",
            "competing_value": KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
        })

        # --- ③ 국소량은 N 에 무관해야 한다 ---
        for q, field in (("defect_fraction", "defect_fraction"),
                         ("psi6_local", "psi6_local")):
            m, se = ref[field]["mean"], ref[field]["se"]
            items.append({
                "quantity": f"{q}__A{A:g}__N256",
                "value": float(m),
                "tolerance": f"±{N_SIGMA * np.sqrt(2.0) * se:.6g}",
                "basis": (
                    f"**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 "
                    f"{m:.6g} ± {se:.4g}. 허용오차 3√2·SE. 어긋나면 `N=100` 에 "
                    f"유한크기 효과가 있었다는 뜻이다"),
                "discriminates": f"{q} 가 N=100 에서 이미 수렴했는가",
                "competing_value": None,
            })

    # --- ④ 절단오차 (카드 §9 의 원래 미결) ---
    for A in AMPLITUDES:
        items.append({
            "quantity": f"beta_u_at_rcut__A{A:g}__N256",
            "value": BETA_U_RCUT[A],
            "tolerance": "±0.0005",
            "basis": (
                f"`r_cut = 0.98·L/2 − buffer`, `L* = √256 = 16` → `r_cut = 7.740`. "
                f"`βU(r_cut) = {A:g}/7.740³`. 결정론적이므로 오차는 반올림뿐이다. "
                f"`N=100` 에서는 `r_cut = 4.800` 이었다"),
            "discriminates": "카드 §9 의 절단오차가 실제로 줄었는가",
            "competing_value": None,
        })

    # --- ⑤ 5-7 불균형·배위수 종류 (결함의 성격) ---
    items.append({
        "quantity": "coord_kinds_aggregate__A10__N256",
        "value": int(parent["A10"]["coord_kinds_aggregate"]),
        "tolerance": "±0",
        "basis": (
            f"`N=100` 집계에서 {parent['A10']['coord_kinds_aggregate']}종. "
            f"**집계 추정량**을 쓴다 — 프레임별 문턱은 `N` 에 의존한다 "
            f"(입자 1개가 `N=100` 에서 1 %, `N=256` 에서 0.39 %). "
            f"근거: findings/fraction-threshold-flips-meaning-*.md"),
        "discriminates": "결함의 성격이 N 에 무관한가",
        "competing_value": None,
    })

    geo_ref = box_si_for_coverage(n_particles=N_REF, sigma_si=SIGMA_SI,
                                  coverage_max=0.10,
                                  d_over_sigma_round=D_OVER_SIGMA)
    geo_new = box_si_for_coverage(n_particles=N_NEW, sigma_si=SIGMA_SI,
                                  coverage_max=0.10,
                                  d_over_sigma_round=D_OVER_SIGMA)
    sc = scales_soft2d(d_si=geo_new["d_si"], sigma_si=SIGMA_SI, T_si=T_SI)

    doc = {
        "question": (
            "N=100 → 256 에서 ① ψ₆ 가 1/√N 으로 줄어드는가 (= 유한크기 바닥인가) "
            "② 국소량(결함 분율·ψ₆ 국소)이 변하지 않는가 ③ 절단오차가 줄어드는가"),
        "parent_run": "runs/2026-07-29_soft-r3-time-resolved",
        "card_open_item": "카드 §9 (N ≥ 252 요구) · §10 (g₆ 지수 η₆ 미구현)",
        "items": items,
        "regimes": {
            "n_ref": N_REF, "n_new": N_NEW, "n_ratio": ratio,
            "seeds": list(SEEDS),
            "amplitudes": list(AMPLITUDES),
            "gamma_zahn": {f"A{A:g}": zahn_phase(A)["gamma"] for A in AMPLITUDES},
            "L_star_ref": geo_ref["L_star"], "L_star_new": geo_new["L_star"],
            "coverage": geo_new["coverage"],
            "note_coverage": (
                "★ 커버리지는 N 에 무관하다 (n* = 1 이라 입자당 면적이 정확히 d²). "
                "상자만 √N 으로 커진다 — L = 150 → 240 µm"),
            "L_si_ref": geo_ref["L_si"], "L_si_new": geo_new["L_si"],
            "tau_d_si": sc.time_si,
            "exponent_liquid": LIQUID_EXPONENT_P,
            "exponent_hexatic_boundary": KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
            "eta6_hexatic_boundary": KTHNY_ETA6_HEXATIC_LIQUID,
        },
        "alternatives": [
            "★ 가장 중요한 경쟁 가설: A=10 의 ψ₆ = 0.248 이 **유한크기 바닥**일 "
            "가능성. p ≈ 0.5 가 나오면 카드 §8.3 의 'hexatic-유사' 판독을 약화시켜야 "
            "한다 — 0.248 은 그저 N=100 의 바닥이었다는 뜻이 된다.",
            "두 점으로는 멱함수 **형태**를 검증할 수 없다. 지수만 뽑는 것이고, "
            "그 한계를 psi6_finite_size_exponent 의 n_points 가 나른다. "
            "형태를 논하려면 셋 이상(예: N=64·144·400)이 필요하다.",
            "A=0.1 과 A=1 은 p ≈ 0.5 가 나와야 정상이다 (깊은 액체). 만약 그것도 "
            "작게 나오면 지수 추정 자체가 틀린 것이므로 **A ≤ 1 이 대조군 역할**을 한다.",
            "에너지/입자는 N 사이에서 직접 비교할 수 없다 — power_law_table 이 "
            "U(r_cut) = 0 으로 이동시키므로 r_cut 이 바뀌면 이동량도 바뀐다. "
            "기록은 하되 예측 항목으로 세우지 않았다.",
            "커버리지 대조군(3.71 %)은 **돌리지 않는다.** 축약 config 가 비트 단위로 "
            "같아서 산술 항등식이다 — tests/test_s7_structure.py::"
            "test_coverage_does_not_touch_the_reduced_config_at_all 로 고정했다.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# S2 PREDICTION — N 수렴 + ψ₆ 유한크기 지수  (자동 생성: "
        "scripts/soft2d_nconv_predict.py)\n"
        "#\n# ⚠ 이 파일은 S5 실행 전에 봉인된다 (SEALED.sha256). 실행 후 수정 금지.\n"
        "# 허용오차는 부모 런(N=100)의 시드 SE 를 스케일링해서 유도했다.\n#\n"
        + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8")
    print(f"→ {OUT.relative_to(REPO)}  ({len(items)} 항목)\n")
    print(f"{'항목':<40} {'예측':>12} {'경쟁':>12} {'허용':>12}")
    print("-" * 80)
    for it in items:
        cv = it.get("competing_value")
        cvs = f"{cv:.5g}" if isinstance(cv, (int, float)) else "—"
        v = it["value"]
        vs = f"{v:.5g}" if isinstance(v, (int, float)) else str(v)
        print(f"{it['quantity']:<40} {vs:>12} {cvs:>12} {it['tolerance']:>12}")
    print(f"\n  N: {N_REF} → {N_NEW} · L*: {geo_ref['L_star']:.0f} → "
          f"{geo_new['L_star']:.0f} · L: {geo_ref['L_si']*1e6:.0f} → "
          f"{geo_new['L_si']*1e6:.0f} µm · 커버리지 {geo_new['coverage']:.4%} (불변)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
