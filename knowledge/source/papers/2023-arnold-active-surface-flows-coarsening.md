---
type: source
kind: arxiv
lab_authored: true
title: Active surface flows accelerate the coarsening of lipid membrane domains
authors:
  - "Arnold DP"
  - "Gubbala A"
  - "Takatori SC"
year: 2023
journal: Phys. Rev. Lett. 131, 128402
doi: 10.1103/PhysRevLett.131.128402
arxiv_id: 2306.00218
source_url: "https://arxiv.org/abs/2306.00218"
raw_file: knowledge/raw/lab/2023-arnold-active-surface-flows-coarsening.pdf
si_available: true
access: open (arXiv)
engine: "실험 + 해석 (Smoluchowski 응집 + phase field)"
reproduced: no
parameters_extracted: yes
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "adjacent"
---
# Active surface flows accelerate the coarsening of lipid membrane domains

## 왜 이 위키에 있는가

능동 흐름이 조대화를 가속. **조대화 지수(coarsening exponent)** 는 우리 클러스터 분석 관측량과 겹침.

**관련도:** `adjacent`

## 초록 (원문)

Phase separation of multicomponent lipid membranes is characterized by the nucleation and
coarsening of circular membrane domains that grow slowly in time as $\sim t^{1/3}$, following
classical theories of coalescence and Ostwald ripening. In this work, we study the coarsening
kinetics of phase-separating lipid membranes subjected to nonequilibrium forces and flows
transmitted by motor-driven gliding actin filaments. We experimentally observe that the
activity-induced surface flows trigger rapid coarsening of non-circular membrane domains that
grow as $\sim t^{2/3}$, a 2$\times$ acceleration in the growth exponent compared to passive
coalescence and Ostwald ripening. We analyze these results by developing analytical theories
based on the Smoluchowski coagulation model and the phase field model to predict the domain
growth in the presence of active flows. Our work demonstrates that active matter forces may be
used to control the growth and morphology of membrane domains driven out of equilibrium.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2306.00218
- 저널판: https://doi.org/10.1103/PhysRevLett.131.128402
- arXiv comment: *Main text is 5 pages with 3 figures. Supplemental materials attached include a supplemental appendix (includes supplemental movie legends, detailed methods, figures, and derivations) and five supplemental movies (S1-S5)*

## 요약

다성분 지질막의 상분리는 원형 도메인이 핵생성 후 **`~t^{1/3}`으로 느리게 자라는** 고전적
조대화(coarsening) 이론을 따른다. 이 논문은 거기에 **모터 구동 액틴 흐름**이라는 비평형 힘을
가한다.

결과: 활성 유도 표면 흐름이 **비원형 도메인의 빠른 조대화**를 유발하고, 성장 지수가
**`~t^{2/3}`으로 2배 가속**된다.

Smoluchowski 응집 모형과 phase field 모형 두 갈래로 해석 이론을 세워 이 성장을 예측했다.

## 파라미터

| 항목 | 값 |
|---|---|
| 관측량 | 도메인 크기 `R(t)` |
| **수동 조대화 지수** | **`R(t) ~ t^{1/3}`** — 고전 이론 |
| **능동 조대화 지수** | **`R(t) ~ t^{2/3}`** — 2배 가속 |
| 도메인 형태 | 수동 = 원형 / 능동 = 비원형 |
| 이론 | Smoluchowski 응집 모형 + phase field |
| 구동 | 모터 구동 액토미오신 수축 흐름 |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### ★ 조대화 지수는 전용 가능한 정량 벤치마크다

`R(t) ~ t^{1/3}`은 **해석적으로 알려진 지수**(Lifshitz–Slyozov / 확산 제한 응집)다. 계가
지질막이든 콜로이드 겔이든, **상분리 후 조대화의 보편 지수**라는 점에서 우리 계에도 적용된다.

> `master_plan.md` §6의 **증거층 ②(해석적 극한)** 에 바로 쓸 수 있다.
> 우리가 인력 콜로이드(depletion gel)의 스피노달 분해를 시뮬레이션하면,
> **조대화 지수가 `1/3`으로 나오는지**가 강력한 자기검사다.

### 벤치마크 후보 ★

- `coarsening_exponent_diffusive_1_3` — 수동 상분리에서 `R(t) ~ t^{1/3}`
  - 계 무관하게 성립하는 보편 지수 → **우리 콜로이드 계에 직접 적용 가능**
  - 측정법: 도메인 크기를 시간의 함수로 → log–log 기울기
  - **주의:** 유체역학이 있으면 `t^{1}`(관성/점성 지배)로 바뀐다. 우리는 HI를 무시하므로(`D11`)
    확산 지배 `1/3`이 맞는 예측이다 — 이것 자체가 `D11`의 자기일관성 검사가 된다

### 관측량 구현 참고

같은 그룹의 후속 논문(`2024-arnold-lipid-domains-actin-viscoelasticity`)에서
**도메인 크기를 정적 구조인자 `S(k)`의 평균 파수로 정의**한다. 우리 `07_observables.md`에
`S(k)`가 이미 있으므로 **`⟨k⟩ → R = 2π/⟨k⟩` 변환만 추가하면 조대화 분석이 된다.**

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
