---
type: source
kind: arxiv
lab_authored: true
title: Learning continuum-level closures for control of interacting active particles
authors:
  - "Quah T"
  - "Takatori SC"
  - "Rawlings JB"
year: 2025
journal: J. Chem. Phys. 164, 044126
doi: 10.1063/5.0300697
arxiv_id: 2501.18809
source_url: "https://arxiv.org/abs/2501.18809"
raw_file: knowledge/raw/lab/2025-quah-continuum-closures-active-control.pdf
si_available: false
access: open (arXiv)
engine: "ABP BD (agent-based) + 신경 연산자 UDE"
reproduced: no
parameters_extracted: yes
ingested_at: 2026-07-27
ingested_by: agent
code_repo: knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md
tags:
  - "direct"
---
# Learning continuum-level closures for control of interacting active particles

## 왜 이 위키에 있는가

**상호작용하는 ABP의 BD**에서 연속체 폐쇄를 학습. 입자↔연속체 대응은 우리 관측량 검증에 직결.

**관련도:** `direct`

## 초록 (원문)

Active matter swarms -- collectives of self-propelled particles that could self-assemble, ferry
microscopic cargo, or endow materials with dynamic properties -- remain hard to steer. In
crowded systems, tracking or controlling individual agents becomes challenging, so strategies
should operate on macroscopic fields like particle density. Yet predicting how density evolves
is difficult due to inter-agent interactions. For model-based feedback control methods -- like
Model Predictive Control (MPC) -- fast, accurate, and differentiable models are crucial.
Detailed agent-based simulations are too slow, necessitating coarse-grained continuum models.
However, constructing accurate closures -- approximations expressing the effect of unresolved
microscopic states (e.g., agent positions) on continuum dynamics -- is challenging for active
matter swarms. We present a learning-for-control framework that learns continuum closures from
agent simulations, demonstrated with active Brownian particles under a controllable external
field. Our Universal Differential Equation (UDE) framework represents the continuum as an
advection-diffusion equation. A neural operator learns the advection term, providing closure
relations for microscopic effects like self-propulsion, interactions, and external field
responses. This UDE approach, embedding universal function approximators in differential
equations, ensures adherence to physical laws (e.g., conservation) while learning complex
dynamics directly from data. We embed this learned continuum model into MPC for precise agent
simulation control. We demonstrate this framework's capabilities by dynamically exchanging
particle densities between two groups, and simultaneously controlling particle density and mean
flux to follow a prescribed sinusoidal profile. These results highlight the framework's
potential to control complex active matter dynamics.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2501.18809
- 저널판: https://doi.org/10.1063/5.0300697

## 요약

**입자 시뮬레이션에서 연속체 폐쇄(closure)를 학습**하는 틀. 이류-확산 방정식 형태를
고정하고, 이류 항만 신경 연산자로 학습한다(Universal Differential Equation). 그래야 질량 보존
같은 물리 법칙은 구조적으로 지켜지면서 복잡한 미시 효과(자기추진·상호작용·외부장 반응)만
데이터에서 배운다.

학습된 연속체 모델을 MPC에 넣어 ABP 무리의 **밀도장을 실시간 제어**한다.

## 시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 모델 | **ABP**, 조종 가능한 정렬 외부장 |
| 배제부피 | **WCA** `U₂(r) = 4ε[(a/r)¹² − (a/r)⁶] + ε`, `r < 2^{1/6}a` |
| 입자 수 | **N = 10⁴** |
| 무차원수 1 | **`(ℓ/δ)² = U₀²τ_R/D_T = 100`** — 런 길이 / 미시 확산길이 |
| 무차원수 2 | `δ = √(D_T/τ_R)` |
| 면적분율 | **`φ = n̄πa²/4`** |
| 박스 | 폭 `W` = 외부장 주기와 일치, 높이는 φ 맞춰 조정 |
| 연속체 해법 | 스펙트럴 + method of lines, 측정 간격 `h` 로 이산화 |
| 비선형 항 | 수정 후진 오일러 — 비선형 항은 현재 스텝 값만 사용 |

## 공개 코드 ★

[`titusswsquah/graybox_abp_mpc`](https://github.com/titusswsquah/graybox_abp_mpc) (MIT).
**논문에 없던 `Δt`·초기배치·평형화 프로토콜이 전부 코드에 있다.**
BD 부분 전체 증류는 → [`2024-quah-graybox-abp-mpc-repo.md`](2024-quah-graybox-abp-mpc-repo.md)

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 참고사항

- **`(ℓ/δ)² = U₀²τ_R/D_T`** 는 우리 원장에 없던 무차원수다. `Pe`와 다른 축 —
  **런 길이를 열적 확산 길이로 잰 것**. 능동계를 다룰 때 `03_units_nondim.md`에 추가할 것
- **`φ = n̄πa²/4`** 에서 `a`가 지름인지 반지름인지 표기가 미묘하다. WCA에서 `2^{1/6}a`를
  차단거리로 쓰므로 `a` = 지름 규약. **φ 정의 시 항상 지름/반지름을 명시해야 한다는 예**
- 입자 시뮬을 **연속체 모델의 학습 데이터**로 쓰는 구조 — 우리 S10 ANALYZE가 내는 관측량이
  이런 용도로도 쓰일 수 있다는 참고
- `Δt`·총 스텝 미확보

### 우리 에이전트와의 관계

이 논문은 **제어**가 목적이라 우리 범위(주어진 계를 올바르게 시뮬레이션) 밖이다. 다만
"입자 → 연속체 관측량" 변환이 우리 S10과 겹치고, `(ℓ/δ)²` 같은 무차원수는 바로 쓸 수 있다.

---

> **재현 상태: `reproduced: no`.** 위 파라미터는 논문에 적힌 값을 옮긴 것이고,
> 우리 코드에서 돌려본 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
