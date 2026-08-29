---
type: source
kind: arxiv
lab_authored: true
title: Motility-induced buckling and glassy dynamics regulate three-dimensional transitions of bacterial monolayers
authors:
  - "Takatori SC"
  - "Mandadapu KK"
year: 2020
journal: arXiv preprint (2020)
arxiv_id: 2003.05618
source_url: "https://arxiv.org/abs/2003.05618"
si_available: true
raw_file: knowledge/raw/lab/2020-takatori-motility-induced-buckling.pdf
access: open (arXiv)
engine: "BD (ABP, 2D 이분산)"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---
# Motility-induced buckling and glassy dynamics regulate three-dimensional transitions of bacterial monolayers

## 왜 이 위키에 있는가

**MIPS·유리동역학**. 조밀 능동계의 BD — 상거동 판별 관측량 참조.

**관련도:** `direct`

## 초록 (원문)

Many mature bacterial colonies and biofilms are complex three-dimensional (3D) structures. One
key step in their developmental program is a transition from a two-dimensional (2D) monolayer
into a 3D architecture. Despite the importance of controlling the growth of microbial colonies
and biofilms in a variety of medical and industrial settings, the underlying physical mechanisms
behind single-cell dynamics, collective behaviors of densely-packed cells, and 3D complex colony
expansion remain largely unknown. In this work, we explore the mechanisms behind the 2D-to-3D
transition of motile Pseudomonas aeruginosa colonies; we provide a new motility-induced, rate-
dependent buckling mechanism for their out-of-plane growth. We find that swarming of motile
bacterial colonies generate sustained in-plane flows. We show that the viscous shear stresses
and dynamic pressures arising from these flows allow cells to overcome cell-substrate adhesion,
leading to buckling of bacterial monolayers and growth into the third dimension. Modeling
bacterial monolayers as 2D fluid films, we identify universal relationships that elucidate the
competition between in-plane viscous stresses, pressure and cell-substrate adhesion.
Furthermore, we show that bacterial monolayers can exhibit crossover from swarming to
kinetically-arrested, glassy-like states above an onset density, resulting in distinct 2D-to-3D
transition mechanisms. Combining experimental observations of P. aeruginosa colonies at single-
cell resolution, molecular dynamics simulations of active systems, and theories of glassy
dynamics and 2D fluid films, we develop a dynamical state diagram that predicts the state of the
colony, and the mechanisms governing their 2D-to-3D transitions.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2003.05618
- 저널판: —
- arXiv comment: *See Supplementary Information in Ancillary files*

## 요약

박테리아 단층이 조밀해지면서 **좌굴(buckling)** 로 3차원 전이를 일으키는 과정을, 능동
콜로이드의 **MIPS·유리 동역학** 관점에서 본다. 실험(박테리아 콜로니)과 ABP 시뮬레이션,
그리고 dynamical facilitation(DF) 이론을 결합했다.

세 가지 영역: 저밀도에서 집단 스웜 → 중간에서 동역학적 이질성 → **고밀도에서 유리 상태로
운동학적 정지(kinetically arrested)**. 고밀도 계에서 개별 입자가 활성화된 hopping을 보이는
것이 수동 유리계와 같다.

## 시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 모델 | **ABP** — 자기추진력 `F = ζU₀q(θ)`, 상호작용력 `Fᴾᵢⱼ` |
| 차원 | **2D**, `q = (cosθ, sinθ)` |
| **병진 브라운 운동** | **무시함** — 자기추진과 상호작용만 |
| 재배향 | `D_R = 1/τ_R`, 연속 무작위 과정 (가우시안 백색잡음 `Λ_R`) |
| **활성도** | **`Pe = U₀τ_R/σ`** = 런 길이 / 입자 크기 |
| 면적분율 | **`φ = nπσ²`** (`n` = 수밀도) |
| **결정화 억제** | **이분산 혼합** — 지름 `σ`, `σ₁`, 비율 **1.4**, 몰분율 2/3 : 1/3 |
| 조사 범위 | `φ` 최대 **0.83** 까지 (유리 영역), `Pe = 2` 등 |
| 압력 분해 | 상호작용 기여 vs 능동력 기여 — **조밀계에서는 상호작용이 압도적** |

> ⚠️ `Δt`, 입자 수, 총 스텝 미확보.

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### ★ `D10`(다분산도) 결정에 직접적인 근거

**"결정화 방해를 피하려고 이분산 혼합을 쓴다"** 를 명시하고, 비율까지 준다 — **1.4, 2/3 : 1/3**.
이건 소프트매터에서 표준 관행(Kob–Andersen 계열)이고, 우리 `D10`이 "v1 단분산"으로 되어 있는데
**조밀계로 가는 순간 필요해진다**는 것을 확인해준다.

> `D10` 갱신 제안: v1 단분산 유지하되, **`φ > 0.6` 정도부터는 이분산을 강제**하는 게이트를
> `05_validation_gates.md`에 넣을 것. 안 그러면 의도치 않게 결정화된 계를 "유리"라고 부르게 된다.

### `Pe` 정의 — 세 가지 해석

논문이 같은 `Pe`를 세 가지로 읽는다. 우리 무차원수 원장에 그대로 옮길 가치가 있다.

| 표현 | 읽는 법 |
|---|---|
| `Pe = U₀τ_R/σ` | 런 길이 / 입자 크기 |
| `Pe = U₀/(D_Rσ)` | 병진 이류 / 회전 확산 |
| `Pe = τ_R/(σ/U₀)` | 회전확산 시간 / 지름 통과 시간 |

### 주의 — 병진 브라운 운동을 껐다

`(2.8)`에서 **병진 브라운 운동을 무시**한다. 우리 BD 기본 모델과 다르다. 재현할 때
`kT = 0`으로 두고 자기추진만 남기는 것인지 확인 필요. **`τ_R`이 열적 기원일 필요가 없다**는
점도 명시되어 있다 (박테리아는 편모로 재배향).

### 벤치마크 후보 (조건부)

- `takatori_active_glass_phi` — `Pe` 고정 시 운동학적 정지가 일어나는 `φ`. **`Δt`·N 미확보로 보류**
- 압력-면적분율 곡선 `P(φ)` — 조밀 능동계에서 상호작용 기여가 능동 기여를 압도. 정성적 검사로 유용

---

> **재현 상태: `reproduced: no`.** 위 파라미터는 논문에 적힌 값을 옮긴 것이고,
> 우리 코드에서 돌려본 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
