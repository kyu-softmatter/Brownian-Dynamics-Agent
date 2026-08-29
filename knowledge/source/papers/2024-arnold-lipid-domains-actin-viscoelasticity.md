---
type: source
kind: arxiv
lab_authored: true
title: Lipid membrane domains control actin network viscoelasticity
authors:
  - "Arnold DP"
  - "Takatori SC"
year: 2024
journal: Langmuir
doi: 10.1021/acs.langmuir.4c03463
arxiv_id: 2406.13218
source_url: "https://arxiv.org/abs/2406.13218"
raw_file: knowledge/raw/lab/2024-arnold-lipid-domains-actin-viscoelasticity.pdf
si_available: true
access: open (arXiv)
engine: "실험 (미세유변학 + 형광)"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "adjacent"
---
# Lipid membrane domains control actin network viscoelasticity

## 왜 이 위키에 있는가

지질막 도메인과 액틴망 점탄성. 랩 산출물.

**관련도:** `adjacent`

## 초록 (원문)

The mammalian cell membrane is embedded with biomolecular condensates of protein and lipid
clusters, which interact with an underlying viscoelastic cytoskeleton network to organize the
cell surface and mechanically interact with the extracellular environment. However, the
mechanical and thermodynamic interplay between the viscoelastic network and liquid-liquid phase
separation of 2-dimensional (2D) lipid condensates remains poorly understood. Here, we engineer
materials composed of 2D lipid membrane condensates embedded within a thin viscoelastic actin
network. The network generates localized anisotropic stresses that deform lipid condensates into
triangular morphologies with sharp edges and corners, shapes unseen in 3D composite gels.
Kinetic coarsening of phase-separating lipid condensates accelerates the viscoelastic relaxation
of the network, leading to an effectively softer composite material over intermediate
timescales. We dynamically manipulate the membrane composition to control the elastic-to-viscous
crossover of the network. Such viscoelastic composite membranes may enable the development of
coatings, catalytic surfaces, separation membranes, and other interfaces with tunable spatial
organization and plasticity mechanisms.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2406.13218
- 저널판: https://doi.org/10.1021/acs.langmuir.4c03463
- arXiv comment: *There are 29 pages with four main figures and one table of contents (TOC) figure in the main text. Ancillary files include a supplemental appendix and three supplemental videos. The supplemental appendix contains three video legends, a supplemental note, and four supplemental figures*

## 요약

지질막 도메인과 액틴 네트워크가 **양방향으로 결합**되어 있음을 보인다. 상분리하는 지질
응축물이 조대화하면 네트워크의 점탄성 완화가 가속되고, 반대로 액틴 네트워크는 도메인 성장을
제약한다.

액틴 그물 크기(1–10 µm)가 재구성 막의 지질 도메인 크기와 비슷한 스케일이라 두 계가 강하게 결합한다.

## 파라미터

| 항목 | 값 |
|---|---|
| **도메인 크기 정의** | **정적 구조인자 `S(k)`의 평균 파수** |
| 액틴 그물 크기 | **1 – 10 µm** — 지질 도메인 크기와 동급 |
| 선행 결과 | 능동 미오신 흐름이 조대화를 2× 가속 |
| 측정 | 미세유변학 (점탄성 완화) + 형광 이미징 |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### ★ 관측량 정의를 그대로 가져올 수 있다

> **"We find the size of lipid domains and actin using the mean wave vector of the static
> structure factor `S(k)`."**

우리 `07_observables.md`에 `S(k)`가 이미 있다. 여기에 **`⟨k⟩ = ∫k·S(k)dk / ∫S(k)dk`,
`R = 2π/⟨k⟩`** 를 추가하면 **클러스터·도메인 크기 관측량**이 생긴다.

산란 실험과 직접 비교 가능한 정의라는 점도 이점 — 우리가 겔·응집 계를 다루면 필요하다.

### 스케일 분리 게이트 참고

그물 크기 ≈ 도메인 크기 → **두 계가 결합**한다. 반대로 스케일이 분리되면 하나를 연속체로
취급할 수 있다. 이런 **스케일 비교를 무차원수로 명시**하는 태도는 `03_units_nondim.md`에
반영할 가치가 있다.

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
