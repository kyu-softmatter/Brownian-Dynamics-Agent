---
type: source
kind: arxiv
lab_authored: true
title: Antibody binding reports spatial heterogeneities in cell membrane organization
authors:
  - "Arnold DP"
  - "Xu Y"
  - "Takatori SC"
year: 2022
journal: Nature Commun. 14, 2884
doi: 10.1038/s41467-023-38525-2
arxiv_id: 2211.12022
source_url: "https://arxiv.org/abs/2211.12022"
si_available: true
raw_file: knowledge/raw/lab/2022-arnold-antibody-binding-heterogeneities.pdf
access: open (arXiv)
engine: "실험 + 모델"
reproduced: no
parameters_extracted: no
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "adjacent"
---
# Antibody binding reports spatial heterogeneities in cell membrane organization

## 왜 이 위키에 있는가

막 조직의 공간 이질성. 랩 산출물.

**관련도:** `adjacent`

## 초록 (원문)

The spatial organization of cell membrane glycoproteins and glycolipids is critical for
mediating the binding of ligands, receptors, and macromolecules on the plasma membrane. However,
we currently do not have the methods to quantify the spatial heterogeneities of macromolecular
crowding on live cell surfaces. In this work, we combine experiment and simulation to report
crowding heterogeneities on reconstituted membranes and live cell membranes with nanometer
spatial resolution. By quantifying the effective binding affinity of IgG monoclonal antibodies
to engineered antigen sensors, we discovered sharp gradients in crowding within a few nanometers
of the crowded membrane surface. Our measurements on human cancer cells support the hypothesis
that raft-like membrane domains exclude bulky membrane proteins and glycoproteins. Our facile
and high-throughput method to quantify spatial crowding heterogeneities on live cell membranes
may facilitate monoclonal antibody design and provide a mechanistic understanding of plasma
membrane biophysical organization.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2211.12022
- 저널판: https://doi.org/10.1038/s41467-023-38525-2
- arXiv comment: *31 pages, 4 figures. See Supporting Information in Ancillary files (1 appendix, 3 movies)*

## 요약

항체 결합이 **세포막 조직의 공간적 이질성**을 보고하는 프로브가 된다는 것을 보인다.
세포외 단백질의 copy number와 세포외 도메인 크기로 막 표면의 붐빔(crowding)을 특성화한다.

## 파라미터

본 스캔에서 시뮬레이션 파라미터를 확인하지 못했다. **미확보.**

| 항목 | 값 |
|---|---|
| 접근 | 실험 (항체 결합) + 모델 |
| 특성화 축 | 세포외 단백질 copy number, 세포외 도메인 크기 |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 우리 범위 밖

세포막 생물물리. 콜로이드 BD/MC와 직접 관련이 없다. **랩 산출물 기록으로만 보관.**

굳이 연결점을 찾자면 "표면 붐빔(crowding)"이 유효 배제부피 개념과 닿아 있지만,
현재 v1 범위(`D3`: BD + 구형 콜로이드 3D)에서 활용할 지점이 없다.

> 이 항목은 **`parameters_extracted: no`로 남긴다.** 억지로 채우면
> 위키가 "관련 있는 척하는 페이지"로 오염된다.

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
