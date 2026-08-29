---
type: source
kind: arxiv
lab_authored: true
title: Anisotropic diffusion in lyotropic chromonic liquid crystal using fluorescence recovery after photobleaching
authors:
  - "Choi KH"
  - "Cheon J"
  - "Jeong J"
  - "Takatori SC"
year: 2026
journal: J. Colloid Interface Sci. 721, 140705
doi: 10.1016/j.jcis.2026.140705
arxiv_id: 2603.02147
source_url: "https://arxiv.org/abs/2603.02147"
raw_file: knowledge/raw/lab/2026-choi-anisotropic-diffusion-lclc-frap.pdf
si_available: true
access: open (arXiv)
engine: "실험 (FRAP)"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "method"
---
# Anisotropic diffusion in lyotropic chromonic liquid crystal using fluorescence recovery after photobleaching

## 왜 이 위키에 있는가

이방성 확산계수를 실험으로 측정. **비등방 D 텐서 검증 데이터** — 등방 BD의 확장 검증에 쓸 수 있음.

**관련도:** `method`

## 초록 (원문)

Anisotropic diffusion governs transport in a wide range of soft and biological materials, where
microstructure and molecular interactions jointly shape how matter moves. Here, we
quantitatively investigate anisotropic molecular transport in lyotropic chromonic liquid
crystals (LCLCs) using fluorescence recovery after photobleaching (FRAP). Disodium cromoglycate
(DSCG) serves as a model LCLC system, and diffusion is measured across isotropic, nematic, and
columnar phases as concentration and temperature are varied. To disentangle the roles of
microstructure and molecular interactions, we employ two fluorescent tracers with distinct
affinities for the LCLC aggregates: Acridine Orange (AO), which intercalates into DSCG
aggregates, and Bodipy, which interacts weakly and remains largely in the aqueous phase.
Fourier-space FRAP analysis independently resolves the parallel and perpendicular diffusion
coefficients for both dyes relative to the liquid-crystal alignment. In the nematic phase,
diffusion becomes anisotropic, with faster transport along the liquid-crystal director. As the
DSCG concentration increases, AO dye molecules that are strongly coupled to the aggregates
exhibit a slowdown in all directions, reflecting enhanced packing and steric confinement of the
LC microstructure. In contrast, weakly interacting Bodipy dye molecules display enhanced
transport along the alignment direction as the DSCG concentration increases in the nematic
regime, suggesting the emergence of microscopic channels that guide motion, analogous to
transport in oriented porous media. These results reveal how the evolving microstructure of
LCLCs controls effective diffusion and provide a quantitative framework for understanding and
designing anisotropic transport in aligned soft materials.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2603.02147
- 저널판: https://doi.org/10.1016/j.jcis.2026.140705
- arXiv comment: *11 pages, 3 figures, 9 videos, 1 supporting information*

## 요약

liotropic chromonic liquid crystal(LCLC)에서 **비등방 분자 확산**을 FRAP으로 정량 측정.
disodium cromoglycate(DSCG)를 모델계로, 농도와 온도를 바꿔 **등방(isotropic) · 네마틱(nematic) ·
컬럼나(columnar) 상**에 걸쳐 확산계수를 쟀다.

미세구조 효과와 분자 상호작용 효과를 분리하기 위해 **친화도가 다른 형광 추적자 두 종**을 썼다:
- **Acridine Orange (AO)** — DSCG 응집체에 삽입(intercalate)됨
- **Bodipy** — 약하게만 상호작용, 대부분 수상에 잔류

**Fourier 공간 FRAP 해석**으로 평행/수직 확산계수 `D∥`, `D⊥`를 독립적으로 분해한다.

## 실험·시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 계 | DSCG 수용액 (LCLC) |
| 상 | isotropic / nematic / columnar |
| 추적자 1 | Acridine Orange, 최종 농도 **3 µM** — 응집체 삽입형 |
| 추적자 2 | Bodipy-NHS, **50 nM** — 약한 상호작용 |
| 여기 파장 | 488 nm |
| 광표백 스팟 | 지름 **70–100 µm** (원형 조리개) |
| 시료 두께 | 조절됨 (µm 스케일) |
| 해석 | **Fourier-space FRAP** → `D∥`, `D⊥` 분리 |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 비등방 확산 — 우리 모델의 확장 경계

우리 v1은 **등방 확산**(스칼라 `D`)을 가정한다. `bdkit/units/`의 `UnitMap`도 `D = kT/γ` 하나만
다룬다. 이 논문은 **`D`가 텐서가 되는 계**를 다루므로, 확장 시 참조점이 된다.

- HOOMD `md.methods.Brownian`은 `gamma`를 입자 타입별 스칼라로 받는다.
  **비등방 항력은 `gamma_r`(회전)과는 다른 문제** — 방향 의존 병진 항력은 기본 제공되지 않는다
- 확장하려면 `md.methods.Langevin`의 이방성 입자 경로나 커스텀 힘이 필요

### 실험 대조 오라클 후보

`D∥/D⊥` 비는 **무차원량**이라 단위 환산 없이 바로 비교 가능하다. 다만 우리가 LCLC를
시뮬레이션할 계획은 없으므로, **지금은 벤치마크가 아니라 "우리가 다루지 않는 물리"의 기록**이다.

### 참고 — 사용자 본인 논문

`Kyu Hwan Choi` 제1저자. 실험 프로토콜과 원 데이터에 접근 가능할 것이므로,
나중에 **실험–시뮬 대조 기능**(`D23+` 미항목화 목록)을 만들 때 첫 시험대상으로 적합.

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
