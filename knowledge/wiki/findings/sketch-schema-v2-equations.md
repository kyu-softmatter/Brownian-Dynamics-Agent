---
type: finding
author: agent
drafted: 2026-07-27
confirmed_by:
question: "손으로 그린 스케치에서 무엇을 읽을 수 있는가? S1 판독 스키마는 무엇을 담아야 하는가?"
cites:
  - bdkit/reading/sketch.py
  - bdkit/reading/fidelity.py
  - tests/test_image.jpeg
  - tests/test_image2.jpeg
related_systems: []
---

# 스케치는 "개념도"가 아니라 "손으로 쓴 문제 정의서"다

**질문:** S1 INTAKE의 이미지 판독 스키마는 무엇을 담아야 하는가?

**답:** 기하만으로는 턱없이 부족하다. **수식·파라미터·목적**을 담아야 한다.
v1 스키마를 실제 그림 2장에 대보고 알았다.

---

## v1이 세운 가정과 그것이 틀린 방식

v1은 그림을 **"계의 위상을 그린 개념도"** 로 가정했다. 프롬프트에까지 못박았다.

> 다음은 그림의 모양으로 알 수 없다. 라벨에 적혀 있지 않으면 무조건 `unknown`이다.
> 온도 · 점도 · 용매 · **상호작용 퍼텐셜** · 전하 · 시간간격 …

규칙 자체는 맞다. 문제는 **담을 필드를 아예 안 만든 것**이다. 그래서 라벨에 적혀 있어도
버릴 수밖에 없었다.

실제 그림 (`tests/test_image.jpeg`, `test_image2.jpeg`) 에 적혀 있던 것:

```
T = 300 K (~27°C)          U(r) = ½ k_t r²        k_t = 10 pN/µm
U_ij = A/r_ij³ ≫ k_BT      r_ij = |r_i − r_j|     Lx = Ly
100 particles (원은 13개만 그려짐)                  no z-direction
free diffusion             find final configuration by minimizing U_tot
```

**절반이 수식이다.** 이건 개념도가 아니라 문제 정의서다.

## `F8` 누락 검사가 스스로 잡았다

가장 값진 대목. v1 스키마로 그림 1을 판독하니:

```
[DOUBT] F8 annotations: 라벨의 숫자 300 이 어느 필드에도 들어가지 않았다
    원인: 판독 중 흘렸거나, 스키마에 해당 필드가 없다. 원문: 'T = 300 K. (~27 C)'
```

**검사가 설계 결함을 검출했다.** "원문 라벨을 그대로 전사하고, 그 숫자가 구조화된 필드에
들어갔는지 대조한다"는 규칙은 판독기의 실수뿐 아니라 **스키마의 부족**도 잡는다. 이건
의도한 것보다 강한 성질이었다.

## v2가 세우는 구분 — "없다" vs "눈으로는 못 본다"

| | v1 | v2 |
|---|---|---|
| 온도 | 필드 없음 | `medium.temperature`. 적혀 있으면 `sketch:label`, 없으면 `unknown` |
| 퍼텐셜 | 필드 없음 | `interactions[]` · `external_potential[]` — `raw`/`form`/`params` 분리 |
| `≫ k_BT` | 개념 없음 | `Potential.regime` — 값이 아니라 부등식 |
| `Lx = Ly` | 개념 없음 | `relations[]` — 값이 아니라 제약 |
| 목적 | **필드 없음** | `objective{kind, target, criterion, raw}` |
| 100 vs 13 | 충돌 처리 없음 | `count`(정본) + `drawn_count`(진단) + `F11` |

`F2`는 그대로 살아 있다 — 이 값들을 `sketch:visual`로 채우면 여전히 FAIL이다.
**적혀 있으면 읽고, 없으면 비운다.**

## 왜 `objective`가 가장 큰 구멍이었나

`"find final configuration by rearrange & minimize total energy"` 는 계 기술이 아니라
**무엇을 계산할지**다. 그리고 이 과제는 **BD 동역학이 아니라 에너지 최소화**다.

스펙이 목적을 못 담으면 S5 PLAN이 엉뚱한 것을 계획한다 — 최소화 문제에 시간적분을 붙이고,
평형화 판정을 하고, MSD를 재고 있을 것이다. 전부 통과하면서 전부 무의미하다. 정확히
`master_plan.md` §6-A가 말하는 **"과정이 조용히 틀리는"** 실패다.

## 수식은 원문과 정규화를 둘 다 남긴다

`Potential.raw`(그림 그대로) + `form`(정규화) + `params`(기호별 값·단위)로 갈랐다.
자유 텍스트 한 덩어리로 두면 S3 무차원화가 `k_t`에 단위를 붙일 수 없다.

`F12`가 이걸 지킨다 — `form`만 있고 `raw`가 비면 DOUBT다. **손글씨 수식은 오독이 잦아
원문 없이는 되짚을 수 없다.** `k_t`인지 `k_B`인지, `r³`인지 `r⁵`인지.

## 남은 한계

- **`27` (섭씨)이 F8에 걸린다.** 300 K의 재기술이라 실질 누락은 아니다. 단위 변환을 아는
  누락 검사는 아직 없다. DOUBT로 두는 게 정직하다 — 사람이 1초에 판단한다.
- **`free diffusion` + `U_ij ≫ k_BT` 의 상충**은 코드가 못 잡는다. 판독기가
  `reader_notes`에 적었고, 게이트1에서 사람이 판단할 몫이다. 이건 `V4`(단계 간 일관성)의
  후보 검사이기도 하다.
- **판독 커버리지가 낮다** — 그림 1은 19칸 중 6칸, 그림 2는 10칸. 정상이다. 스케치는
  원래 부분 정보이고, 나머지는 S2 ELICIT가 물어본다. 커버리지를 높이려 빈칸을 채우는
  판독기가 나쁜 판독기다.
