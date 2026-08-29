---
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
status: proposed
cause_class: analysis
stage: S2
question: "tolerance 대역 안에 들어온 PASS 는 예측이 맞았다는 뜻인가?"
answer: "아니다. 넓은 tolerance 는 통계적으로 유의한 어긋남을 가린다. 실제로 3.54σ 어긋난 항목이 ±2% 대역 안에서 PASS 로 나왔고, 원인은 알려진 편향을 예측에 넣지 않은 것이었다."
cites:
  - simbot/validate.py
  - examples/trap-2d-5um/prediction.yaml
  - knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
runs: [2026-07-28_cli-e2e-test]
affects_docs: [master_plan.md]
---

# `PASS` 가 3.54σ 어긋날 수 있다 — 넓은 tolerance 가 가리는 것

## 증상

`cli.py` 로 파이프라인을 처음 관통했을 때 대조표가 이렇게 나왔다:

| 양 | 예측 | 측정 | tolerance | 편차 | verdict |
|---|---|---|---|---|---|
| `msd_plateau_star` | `4` | `4.02149 ± 0.00607` | `±2%` | **+0.54 % (3.54σ)** | **PASS** |

편차가 대역(`±2 %` = `±0.08`) 안이므로 `PASS` 다. 그런데 통계오차 기준으로는
**3.54σ** — 우연이라면 4000분의 1 확률이다. **PASS 와 "예측이 맞았다"가 갈렸다.**

## 진단 경로

1. **측정이 틀렸나** — 아니다. `plateau_over_2d_var = 0.99961 ± 0.00151` 이 독립
   경로(시계열 vs 스냅샷)의 일관성을 0.04 % 로 확인한다.
2. **통계오차가 과소평가됐나** — 아니다. 시드 4개 앙상블이고 `⟨x*²⟩` 의 SE 는
   0.45 % 로 크다. `plateau` 의 SE(0.15 %)가 더 작은 것은 MSD 피팅이 500 프레임을
   쓰기 때문이고, 자기일관성 검사가 이를 뒷받침한다.
3. **예측이 틀렸다** — 맞다. 예측을 `2d = 4.0` **정확히** 로 썼다.

## 근본 원인

`cause_class: analysis` (예측 작성 단계). 조화 트랩 MSD plateau 는

```
plateau = 2d · ⟨x*²⟩
```

이고, `⟨x*²⟩` 는 **Euler–Maruyama 편향을 갖는다** (`1/(1−Δt*/2) = 1.0025063`).
따라서 plateau 의 옳은 예측은

```
2d(1 + bias) = 4 × 1.0025063 = 4.010025      (2D, Δt* = 5e-3)
```

`4.0` 을 쓴 것은 **알려진 편향을 예측에 넣지 않은 것**이다. 편향 0.25 % 가
plateau 의 통계오차 0.15 % 보다 크기 때문에 유의하게 드러났다.

> ★ 같은 편향이 `⟨x*²⟩` 에서는 드러나지 않았다 — 거기서는 SE(0.45 %)가 편향(0.25 %)보다
> 커서 `INCONCLUSIVE` 였다. **같은 물리를 두 관측량으로 보면 검정력이 다르다.**
> plateau 가 더 정밀한 창이었고, 그래서 예측의 허술함이 거기서 먼저 터졌다.

## 처방

### ① 예측 수정 (근본 원인)

```diff
  - quantity: msd_plateau_star
-   value: 4.0
-   basis: MSD(t->inf) = 2d <x*^2> (축약 단위)
+   value: 4.010025062656641
+   basis: '... ★ <x*^2> 가 EM 편향을 갖고 있으므로 plateau 도 그만큼 높다:
+     2d(1+bias) = 4 x 1.0025063. 정확히 2d = 4 를 예측하면 3.54 sigma 어긋난다'
```

수정 후: 편차 `+0.29 % (1.89σ)` — 유의하지 않다.

### ② 구조적 방지 — 판정기가 이 상태를 스스로 표시한다

`simbot.validate.compare` 가 `PASS` 이면서 `|편차|/SE > 3` 인 행에 플래그를 붙이고,
`ValidationReport.problems` 에 올린다:

```python
if row.verdict == PASS and row.sigma > significance_sigma:
    row.flags.append("significant_deviation_within_tolerance")
```

리포트에는 `PASS ⚑` 로 나오고, 이유가 적힌다:

> `msd_plateau_star`: **PASS 이지만 편차가 3.54σ** 다 — tolerance 대역 안이지만
> 통계적으로 유의한 어긋남이다. tolerance 가 너무 넓어 검증이 무력화됐는지, 아니면
> 예측에 넣지 않은 알려진 편향이 있는지 확인할 것

## 왜 중요한가

`master_plan.md` §S2 실패모드에 이미 있던 항목이다:

> tolerance 를 너무 넓게 잡아 어떤 결과든 PASS (검증 무력화) — **금지, 리뷰 대상**

그런데 그것을 **기계적으로 검사할 방법**이 없었다. tolerance 가 넓은지는 그 자체로는
알 수 없고, **측정의 통계오차와 비교해야** 알 수 있다. 이 플래그가 그 비교다.

두 가지 원인이 같은 증상을 만들고, 둘 다 알아야 한다:

| 원인 | 처방 |
|---|---|
| tolerance 가 통계정밀도보다 훨씬 넓다 | tolerance 를 좁힌다 (또는 왜 넓은지 적는다) |
| 예측에 알려진 편향을 넣지 않았다 | 예측을 고친다 ← **이번 경우** |

## 재발 방지

| 장치 | 위치 |
|---|---|
| PASS + 3σ 를 표시하는지 | `simbot/validate.py::compare` (`significance_sigma=3.0`) |
| 리포트에 드러나는지 | `tests/test_s8_report.py::test_validation_problems_are_surfaced` |
| 예시 예측이 편향을 포함하는지 | `examples/trap-2d-5um/prediction.yaml` `msd_plateau_star` |

## 적용 범위 / 한계

- `stat_err` 가 없으면 이 검사는 작동하지 않는다 — 그 경우는 애초에 `INCONCLUSIVE` 다.
- 단측 경계(`R² > 0.99`)에서는 σ 가 크게 나오는 것이 정상이다 (경계에서 멀수록 좋으므로).
  현재 구현은 단측 경계에도 플래그를 걸 수 있으므로, 그런 항목이 늘면 예외 처리가 필요하다.
  → 미결로 남긴다.

## 참고

- [[hoomd-brownian-scheme-and-noise]] — EM 편향 `dt*/2` 의 출처
- [[ks-test-needs-independent-samples]] — 같은 부류(분석 단계에서 나온 잘못된 판정)의 앞선 사례
- `master_plan.md` §S2 실패모드
