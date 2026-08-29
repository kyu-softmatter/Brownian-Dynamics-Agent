---
type: source
kind: paper
lab_authored: false
title: Direct measurements of the effects of salt and surfactant on interaction forces between colloidal particles at water-oil interfaces
authors:
  - "Park BJ"
  - "Pantina JP"
  - "Furst EM"
  - "Oettel M"
  - "Reynaert S"
  - "Vermant J"
year: 2008
journal: Langmuir 24(5), 1686-1694
doi: 10.1021/la7008804
source_url: "https://doi.org/10.1021/la7008804"
raw_file: knowledge/raw/lab/2008-park-salt-surfactant-interface-forces.pdf
si_available: false
access: paywall
engine: "실험 (광집게)"
reproduced: no
parameters_extracted: yes
affiliation: University of Delaware (Furst EM group) · Mainz (Oettel) · KU Leuven (Vermant)
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "method"
  - "interfacial-colloid"
  - "drag-coefficient"
---
# Direct measurements of the effects of salt and surfactant on interaction forces between colloidal particles at water−oil interfaces

## 왜 이 위키에 있는가

★ **계면 유효 점도 `η_eff` 의 원출처.** BD 항력계수가 여기서 정해진다.

**Park BJ 와 Furst EM 의 공저** — Park 교수가 Delaware Furst 랩에 있던 시기.
[[goal-autonomous-paper-to-sim-verification]] 이 지목한 두 학자가 여기서 만난다.

`Received March 26, 2007. In Final Form: October 29, 2007`

---

## 1. ★ `η_eff` — 확정

본문 p.1689, *Materials and Methods*. 계면에 포획된 입자를 잡고 스테이지를 등속 `U`로 움직여
**Stokes 항력과 광집게 힘의 균형**으로 트랩을 보정한다.

$$F_S = 6\pi a\,\eta_{eff}\,U, \qquad
\boxed{\;\eta_{eff} = \frac{\eta_{oil}(1-\cos\theta) + \eta_{water}(1+\cos\theta)}{2}\;}$$

`a` = 입자 **반지름**, `θ` = 삼상 접촉각.

### 이게 무슨 뜻인가 — 표면적 가중 평균이다

계면에 걸친 구에서 각 상에 잠긴 **표면적 분율**이 정확히 저 계수다:

| | 분율 | |
|---|---|---|
| 오일 상 | `(1 − cos θ)/2` | |
| 물 상 | `(1 + cos θ)/2` | |

따라서 `η_eff = η_oil·(오일 분율) + η_water·(물 분율)` — **두 벌크 점도의 표면적 가중 평균.**

**극한 검산**

| `θ` | `cos θ` | `η_eff` | 해석 |
|---|---|---|---|
| 0° | +1 | `η_water` | 완전히 물 속 |
| 90° | 0 | `(η_oil + η_water)/2` | 단순 평균 |
| 180° | −1 | `η_oil` | 완전히 오일 속 |

> **`θ` 는 물 상 기준으로 잰다.** 본문이 *"The contact angle increases as SDS is added; that is,
> particles are pushed into the oil phase"* 라고 명시한다 — `θ` 증가 = 오일 쪽으로 이동.
> 규약을 뒤집어 넣으면 `η_eff` 가 두 상 사이에서 반대로 간다.

### BD 항력계수

$$\gamma = 6\pi a\,\eta_{eff} = 3\pi d\,\eta_{eff}$$

`2020-choi` SI 의 `κ_t = 3πdη_eff u/Δx` 와 **같은 식이다** (`a` = 반지름, `d` = 지름).
두 논문의 표기가 달랐을 뿐 내용은 동일 — 이제 `η_eff` 의 정의까지 채워졌다.

## 2. 성립 조건과 한계 — 반드시 함께 기록할 것

| 항목 | 상태 |
|---|---|
| 유체쌍 | **n-decane / 탈이온수.** decane 은 알루미나 흡착으로 극성 성분 제거 |
| 표면점도 항 | ❌ **없다.** 논문 전체에 surface viscosity · interfacial rheology · Boussinesq 언급 0건 |
| 계면활성제(SDS) 존재 시 | ⚠️ **같은 식을 그대로 쓴다.** SDS 는 `θ` 를 바꾸는 경로로만 반영되고, SDS 가 만드는 **표면 점탄성은 `η_eff` 에 들어가지 않는다** |
| `η_oil`·`η_water` 수치 | ❌ 논문에 없음 (표준 물성으로 대입한 것으로 보임) |
| `θ` 수치 | ❌ 본문에 표 없음 — Figure 2A 그래프, 출처는 ref 14 (Reynaert et al.) |

> **가장 중요한 한계:** 이 `η_eff` 는 **벌크 두 상의 기하 평균일 뿐 표면점도를 무시한다.**
> 계면활성제가 있는 계에서 `Bq = η_s/((η₁+η₂)a)` 가 1을 넘으면 이 식은 항력을 과소평가한다.
> 이 논문은 SDS 를 CMC 훨씬 아래로만 썼으므로 그 영역에서 검증된 것으로 본다.
> **SDS 농도가 높거나 다른 계면활성제를 쓰는 계로 옮길 때는 재검토 대상.**

### decane/water 에서는 `θ` 의존성이 사실상 사라진다

`η_water ≈ 0.89`, `η_decane ≈ 0.84 mPa·s` (25 °C, **문헌 상온값 — 이 논문에 없음, 카드 등재 전
출처 확보 필요**). 두 값이 6% 이내라 가중치가 어떻게 바뀌어도 `η_eff` 가 거의 안 움직인다:

| `θ` | `η_eff` (추정) |
|---|---|
| 99.4° (SPS) | ≈ 0.861 mPa·s |
| 142.6° | ≈ 0.845 mPa·s |
| 145.4° | ≈ 0.844 mPa·s |

99°~145° 전 구간에서 변동 ~2%, `142.6°` vs `145.4°` 는 **0.1% 미만**.

> 이것이 [`2020-choi`](2020-choi-electrostatic-self-potential-heterogeneity.md) 의 APS/CPS
> 접촉각 불일치를 "아무거나 써도 된다"로 처리해도 되는 **정량적 근거**다.
> 단, 이 논거는 **decane/water 에만** 유효하다. 물/실리콘오일처럼 점도차가 큰 계면에서는
> `θ` 가 `η_eff` 를 크게 흔든다.

## 3. 보조 보정법 — 진동 트랩

메니스커스 곡률이 커서 등속 견인이 안 되는 경우, 광집게를 정현파로 진동시켜
**"phase lock → phase slip" 전이 주파수**로 보정한다 (Faucheux, Stolovitzky, Libchaber,
*Phys. Rev. E* 1995, 51). 계면이 평평하지 않은 셋업에서 쓸 수 있는 대안.

레이저 세기를 바꿔도 상호작용력 측정이 불변임을 확인했다 — 광집게가 입자에 알짜힘을
주지 않는다는 근거.

## 4. 계

| 항목 | 값 |
|---|---|
| 오일 | n-decane (Acros, 99+%), 극성 성분 제거 |
| 물 | 탈이온수 |
| 입자 | 설페이트 안정화 폴리스티렌, **지름 3.1 ± 0.2 μm**, 표면전하 **9.1 μC/cm²** (IDC) |
| 첨가물 | SDS (CMC 이하) · NaCl |
| 접촉각 측정 | ① PS 캐스트 필름 + 고니오미터(정·역방향으로 이력 검사) ② gel trapping (Paunov) |

## 5. 이 논문의 물리 주장 (참고)

계면 콜로이드의 장거리 반발이 **오일 상의 소수 해리 전하**에 의한 것이라는 가설을 염·계면활성제
효과로 검증. 거리 의존 지수 실측 `0.43 ± 0.04` 관련 논의, 쌍마다 편차가 커서 **다수 쌍 평균이
필수**라는 방법론적 결론. 입자 쌍별 이질성의 원인 후보로 표면전하 불균질·나노 거칠기·젖음
불균질을 든다 — [`2020-choi`](2020-choi-electrostatic-self-potential-heterogeneity.md) 의
`Ω` 불균질 논의로 이어지는 대목.

## 6. 우리 쪽 쓸모

**① 계면 콜로이드 카드 §3 기준 단위의 항력 행이 이걸로 채워진다.** `γ = 6πaη_eff` 에
`η_eff` 정의식과 성립 조건까지 `[출처]` 로 붙는다.

**② 인용 사슬의 종점.** `2020-choi`(SI ref 2) 와 `2023-lee`(ref 22) 가 독립적으로 여기를
가리켰고, 실제로 정의가 여기 있다. 추적 종료.

**③ 한계 기록이 곧 S0 필터 입력.** 표면점도를 무시한다는 사실이 명시됐으므로, 계면활성제가
중요한 계에서 이 `η_eff` 로 BD 를 돌린 결과는 **그 근사를 리포트에 표기**해야 한다.

## 7. 남은 것

| 항목 | 상태 |
|---|---|
| `η_oil`·`η_water` 수치 출처 | ⚠️ 표준 물성표에서 확보 필요 (이 논문에 없음) |
| `θ` 수치 | ⚠️ Figure 2A 그래프. 필요하면 이미지에서 판독하거나 ref 14 확보 |
| 온도 | ❌ 명시 없음 (ambient 추정) |
