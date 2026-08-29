"""S2 — 유한크기 사다리(S17) 예측을 앞선 두 런에서 유도해서 봉인한다.

usage: python scripts/soft2d_fss_predict.py

## 이 런이 검정하는 것은 **방법**이다

`runs/2026-07-29_soft-r3-nconv` 가 두 점(`N=100, 256`)으로 `p` 를 얻어
"`A=10` 은 hexatic 이 아니다"(`3.5σ`) 를 주장했다. 두 점은 멱함수를 **가정**한다.
이 런은 4점으로 그 가정을 검증한다 — 형태가 틀렸으면 `η₆ = 1.46` 도 못 쓴다.

동시에 세 가지가 달라진다 (**독립 측정**이 되도록):
  · `r_cut` 3.80 고정 (전 N 동일) vs 앞 런의 자연 4.80/7.74
  · `prod_tau` 30 τ_d vs 80 · 창 [20,30] vs [60,80]
  · 시드 32–35 vs 21–24
⇒ 같은 `p` 가 나오면 **방법이 견고하다**는 뜻이고, 다르면 어느 축이 원인인지
  좁혀야 한다.
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
                                       LIQUID_EXPONENT_P)

NCONV = REPO / "runs" / "2026-07-29_soft-r3-nconv" / "metrics.json"
OUT = REPO / "examples" / "soft-r3-fss" / "prediction.yaml"
N_LADDER = (64, 144, 256, 400)
AMPLITUDES = (0.1, 1.0, 10.0)
CHI2_MAX = 3.0


def main() -> int:
    nc = json.loads(NCONV.read_text())
    items: list[dict] = []

    # --- ① 형태 검증: 이 런의 존재 이유 -------------------------------------
    for A in AMPLITUDES:
        items.append({
            "quantity": f"chi2_reduced__A{A:g}",
            "value": 1.0,
            "tolerance": f"<{CHI2_MAX:g}",
            "basis": (
                f"로그-로그 직선 적합의 `χ²/dof` (dof = 4−2 = 2). 멱함수가 맞고 "
                f"오차막대가 정직하면 기대값 1 이다. `< {CHI2_MAX:g}` 를 통과 기준으로 "
                f"사전등록한다 — 넘으면 단일 지수로 요약하면 안 된다"),
            "discriminates": "ψ₆(N) 이 실제로 멱함수인가 (두 점 추정의 전제)",
            "competing_value": None,
        })

    # --- ② 지수: 액체 대조군 + A=10 재측정 ----------------------------------
    for A in AMPLITUDES:
        f = nc[f"A{A:g}"]["finite_size"]
        if A <= 1.0:
            items.append({
                "quantity": f"psi6_exponent_p__A{A:g}",
                "value": LIQUID_EXPONENT_P,
                "tolerance": "±0.12",
                "basis": (
                    f"깊은 액체이므로 `p = 1/2`. 2점 런에서 {f['p']:.3f} ± "
                    f"{f['p_se']:.3f} 였고 4점이면 SE 가 줄어야 한다. "
                    f"허용오차 ±0.12 는 2점 SE({f['p_se']:.3f})의 약 1.5배 — "
                    f"**대조군이므로 좁게 잡는다**"),
                "discriminates": "지수 추정이 건강한가 (대조군)",
                "competing_value": KTHNY_ETA6_HEXATIC_LIQUID / 4.0,
            })
        else:
            items.append({
                "quantity": f"psi6_exponent_p__A{A:g}",
                "value": float(f["p"]),
                "tolerance": f"±{3.0 * np.sqrt(2.0) * f['p_se']:.4g}",
                "basis": (
                    f"2점 런의 {f['p']:.3f} ± {f['p_se']:.3f} 를 재측정한다. "
                    f"`r_cut`·런 길이·시드가 모두 다르므로 **독립 측정**이다. "
                    f"허용오차 3√2·SE"),
                "discriminates": (
                    "2점 추정이 r_cut·런길이·시드에 견고한가"),
                "competing_value": LIQUID_EXPONENT_P,
            })

    # --- ③ A=10 이 여전히 hexatic 을 기각하는가 ------------------------------
    items.append({
        "quantity": "eta6_minus_3sigma__A10",
        "value": 0.5,
        "tolerance": f">{KTHNY_ETA6_HEXATIC_LIQUID:g}",
        "basis": (
            f"`η₆ − 3·SE` 가 hexatic 상한 {KTHNY_ETA6_HEXATIC_LIQUID:g} 보다 크면 "
            f"hexatic 을 `3σ` 로 기각한다. 2점 런에서 "
            f"{nc['A10']['finite_size']['eta6']:.2f} − 3×"
            f"{nc['A10']['finite_size']['eta6_se']:.2f} = "
            f"{nc['A10']['finite_size']['eta6'] - 3*nc['A10']['finite_size']['eta6_se']:.2f} "
            f"였다. 예측값 0.5 는 '넉넉히 넘는다'는 뜻의 자리표시자다"),
        "discriminates": "A=10 의 hexatic 기각이 4점에서도 유지되는가",
        "competing_value": None,
    })

    # --- ④ r_cut 8.4배 여행에도 ψ₆ 가 같은가 (S16 확장) ----------------------
    for A in AMPLITUDES:
        ref = nc[f"A{A:g}"]
        m, se = ref["psi6_global"]["mean"], ref["psi6_global"]["se"]
        items.append({
            "quantity": f"psi6_global__A{A:g}__N256",
            "value": float(m),
            "tolerance": f"±{3.0 * np.sqrt(2.0) * se:.4g}",
            "basis": (
                f"같은 `N=256` 을 `r_cut = 3.80` 으로 다시 잰다 (앞 런은 7.740 — "
                f"**절단오차 8.4배 차이**, βU: 0.0216 → 0.182). 앞 런 "
                f"{m:.5g} ± {se:.4g}. 일치하면 S16(절단오차가 관측량을 편향시키지 "
                f"않는다)이 더 큰 여행에서도 성립한다. 런 길이·시드도 다르다"),
            "discriminates": "절단오차 8.4배가 ψ₆ 를 편향시키는가",
            "competing_value": None,
        })

    doc = {
        "question": (
            "ψ₆(N) 이 실제로 멱함수인가 (4점) · 그 지수가 r_cut·런길이·시드에 "
            "견고한가 · A=10 의 hexatic 기각이 유지되는가"),
        "parent_runs": ["runs/2026-07-29_soft-r3-nconv",
                        "runs/2026-07-29_soft-r3-time-resolved"],
        "card_open_item": "카드 §8.5 S17 (멱함수 형태 검증)",
        "design": {
            "n_ladder": list(N_LADDER),
            "lever_arm_ln": float(np.log(N_LADDER[-1] / N_LADDER[0])),
            "amplitudes": list(AMPLITUDES),
            "seeds": [32, 33, 34, 35],
            "seed_screening": (
                "min_sep = 0.8 d 기각표집이 시드마다 실패한다 (성공률 N=64 98.3 % · "
                "N=144 95.0 % · N=256/400 100 %, 시드 31–90 실측). 전 N 에서 "
                "성공하는 시드를 골라 **짝지은 설계**로 만들었다 — 시드가 N 마다 "
                "다르면 ψ₆(N) 비교에 다른 초기배치 앙상블이 섞인다"),
            "r_cut_fixed": 3.80,
            "r_cut_rationale": (
                "자연 r_cut 은 N 과 함께 커져 절단오차가 8배 변한다 "
                "(N=64: βU=0.179 → N=400: 0.011). 그러면 ψ₆(N) 기울기에 절단오차 "
                "추세가 섞인다. 최소 상자(N=64, L/2=4.0)가 허용하는 3.80 으로 "
                "고정하면 절단오차가 0.182 로 커지지만 **전 N 동일**하므로 거짓 "
                "기울기를 만들 수 없다"),
            "prod_tau": 30.0, "window": [20.0, 30.0],
            "prod_tau_rationale": (
                "τ_relax ≈ 0.098 τ_d (§8.4) 이므로 20 τ_d 는 완화시간의 200배다. "
                "80 τ_d 를 쓰면 A=10·N=400 이 예산(600 s/런)을 넘는다"),
            "chi2_max": CHI2_MAX,
        },
        "items": items,
        "alternatives": [
            "★ 형태가 깨지는 경우(χ²/dof > 3)가 가장 정보량이 크다. 그러면 "
            "'A=10 은 hexatic 이 아니다' 를 단일 지수로 주장한 것이 무효가 되고, "
            "카드 §8.5 를 다시 약화시켜야 한다. 그 가능성을 열어 둔다.",
            "N=64 는 L* = 8 이라 r_cut = 3.80 이 L/2 의 95 % 다 — 최소이미지 여유가 "
            "거의 없다. 그 점이 잔차를 지배하면 3점(144·256·400)으로 다시 봐야 한다.",
            "A=10 의 p 가 2점 런과 다르게 나올 수 있다. r_cut·런길이·시드 세 축이 "
            "동시에 바뀌므로 어긋나면 **어느 축이 원인인지 이 런만으로는 못 가른다** — "
            "그때는 한 축씩 되돌리는 런이 필요하다.",
            "이 런은 hexatic 창(A = 10.03–10.75)을 훑지 않는다. 방법 검증이 먼저이고, "
            "형태가 확인된 뒤에야 그 창의 η₆ 를 믿을 수 있다.",
        ],
    }
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# S2 PREDICTION — 유한크기 사다리 (S17, 형태 검증)  (자동 생성: "
        "scripts/soft2d_fss_predict.py)\n#\n"
        "# ⚠ 이 파일은 S5 실행 전에 봉인된다 (SEALED.sha256). 실행 후 수정 금지.\n#\n"
        + yaml.safe_dump(doc, allow_unicode=True, sort_keys=False, width=100),
        encoding="utf-8")
    print(f"→ {OUT.relative_to(REPO)}  ({len(items)} 항목)\n")
    print(f"{'항목':<34} {'예측':>10} {'허용':>12} {'경쟁':>10}")
    print("-" * 70)
    for it in items:
        cv = it.get("competing_value")
        print(f"{it['quantity']:<34} {it['value']:>10.5g} "
              f"{it['tolerance']:>12} "
              f"{(f'{cv:.5g}' if isinstance(cv, (int, float)) else '—'):>10}")
    print(f"\n  사다리 {list(N_LADDER)} · 레버암 ln(400/64) = "
          f"{np.log(400/64):.3f} · dof = {len(N_LADDER)-2}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
