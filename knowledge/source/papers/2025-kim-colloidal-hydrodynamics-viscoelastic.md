---
type: source
kind: arxiv
lab_authored: true
title: Colloidal hydrodynamic interactions in viscoelastic fluids
authors:
  - "Kim DY"
  - "Nagella SG"
  - "Malik S"
  - "Park N"
  - "Nam J"
  - "Shaqfeh ESG"
  - "Takatori SC"
year: 2025
journal: Soft Matter
doi: 10.1039/D5SM00874C
arxiv_id: 2508.11948
source_url: "https://arxiv.org/abs/2508.11948"
raw_file: knowledge/raw/lab/2025-kim-colloidal-hydrodynamics-viscoelastic.pdf
si_available: true
access: open (arXiv)
engine: "실험 (광집게) + 점탄성 이론"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "method"
---
# Colloidal hydrodynamic interactions in viscoelastic fluids

## 왜 이 위키에 있는가

점탄성 유체 내 콜로이드 HI. **우리 v1은 HI를 무시한다(D11)** — 그 근사가 언제 깨지는지의 경계 자료.

**관련도:** `method`

## 초록 (원문)

The motion of suspended colloidal particles generates fluid disturbances in the surrounding
medium that create interparticle interactions. While such colloidal hydrodynamic interactions
(HIs) have been extensively studied in viscous Newtonian media, comprehensive understanding of
HIs in viscoelastic fluids is lacking. We develop a framework to quantify HIs in viscoelastic
fluids with high spatiotemporal precision by trapping colloids and inducing translation-rotation
hydrodynamic coupling. Using solutions of wormlike micelles (WLMs) as a case study, we discover
that HIs are strongly time-dependent and depend on the structural memory generated in the
viscoelastic fluid, in contrast to "instantaneous" HIs in viscous Newtonian fluids. We directly
measure time-dependent HIs between a stationary probe and a driven particle during transient
start-up, developing on the WLM relaxation timescale. Following the sudden cessation of the
driven particle, we observe an intriguing flow reversal in the opposing direction, lasting for a
time about ten times larger than the WLM relaxation time. We corroborate our observations with
analytical microhydrodynamic theory, direct numerical solutions of a continuum model, and
particle-based Stokesian dynamics simulations. We find that the structural recovery of the WLMs
from a nonlinear strain can generate anisotropic and heterogeneous stresses that produce flow
reversals and hydrodynamic attraction among colloids. Measured heterogeneities indicate a
breakdown of standard continuum models for constitutive relations when the size of colloids is
comparable to the length scales of the polymeric constituents and their entanglement lengths.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2508.11948
- 저널판: https://doi.org/10.1039/D5SM00874C
- arXiv comment: *22 pages, 7 figures. Supplementary Information and videos available as ancillary files*

## 요약

위 논문(광집게 다체 HI)의 **점탄성 유체 확장판**. Newtonian 유체에서 HI가 "순간적
(instantaneous)"인 것과 달리, 점탄성 유체에서는 **HI가 시간에 의존**한다 — 유체에 남은
구조적 기억(structural memory) 때문이다.

worm-like micelle(WLM) 용액을 모델계로 써서, 구동 시작(transient start-up) 동안 정지 탐침과
구동 입자 사이의 **시간의존 HI**를 직접 측정했다.

## 실험·시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 매질 | worm-like micelle (WLM) 용액 |
| 그물 크기 | **`ξ_M ≈ 120 – 200 nm`** |
| 윤곽 길이 | `L_c` (본문 값 참조) |
| 입자 | 지름 5 µm |
| 구동 속도 | **`U₁ = 60 µm/s`** |
| 구성방정식 | Oldroyd-B 계열 — `λ ∂τ/∂t + τ + O(Wi) = 2η₀γ̇` |
| 무차원수 | **Weissenberg 수 `Wi`** |
| 주파수 영역 | 복소 점도를 쓴 Stokes 흐름

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 참고사항

- **`Wi`(Weissenberg 수)** 가 등장한다. 우리 무차원 원장에 없는 축 — 점탄성을 다루게 되면
  추가해야 한다. `Wi = λγ̇` (완화시간 × 전단율)
- **시간의존 HI** 는 우리 BD 모델이 구조적으로 표현할 수 없는 물리다.
  마찰계수 `γ`가 상수인 순간 이미 기억이 없는 유체를 가정한 것
- `ξ_M ≈ 120–200 nm` vs 입자 5 µm — **그물 크기 ≪ 입자 크기**라 연속체 취급이 성립.
  우리가 복잡 유체를 다루게 되면 이 스케일 분리를 게이트로 걸어야 한다

### 우리 범위 밖임을 명시하는 용도

`master_plan.md` §12 비목표에 *"HI를 하지 않는다"* 가 있는데, 이 논문은 **점탄성까지 가면
HI가 시간의존까지 된다**는 것을 보여준다. v1 범위를 방어하는 근거로 인용할 수 있다.

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
