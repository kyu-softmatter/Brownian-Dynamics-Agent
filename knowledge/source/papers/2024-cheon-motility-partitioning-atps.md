---
type: source
kind: arxiv
lab_authored: true
title: Motility modulates the partitioning of bacteria in aqueous two-phase systems
authors:
  - "Cheon J"
  - "Choi KH"
  - "Modica KJ"
  - "Mitchell RJ"
  - "Takatori SC"
  - "Jeong J"
year: 2024
journal: Phys. Rev. Lett. (2025)
doi: 10.1103/6gm5-cnv1
arxiv_id: 2405.08995
source_url: "https://arxiv.org/abs/2405.08995"
raw_file: knowledge/raw/lab/2024-cheon-motility-partitioning-atps.pdf
si_available: true
access: open (arXiv)
engine: "실험 (광집게) + BD"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---
# Motility modulates the partitioning of bacteria in aqueous two-phase systems

## 왜 이 위키에 있는가

**운동성 입자의 상분배**. ABP BD와 실험 대조가 같은 논문에 있음 — 증거층 ③④ 후보.

**관련도:** `direct`

## 초록 (원문)

We study the partitioning of motile bacteria in an aqueous two-phase mixture of dextran (DEX)
and polyethylene glycol (PEG), which can phase separate into DEX-rich and PEG-rich phases. While
non-motile bacteria partition exclusively into the DEX-rich phase in all conditions tested, we
observed that motile bacteria penetrate the soft DEX/PEG interface and partition variably among
the two phases. For our model organism \textit{Bacillus subtilis}, the fraction of motile
bacteria in the DEX-rich phase increased from 0.58 to 1 as we increased DEX composition within
the two-phase region. We hypothesized that the chemical affinity between DEX and the bacteria
cell wall acts to weakly confine the bacteria within the DEX-rich phase; however, motility can
generate sufficient mechanical forces to overcome the soft confinement and propel the bacteria
into the PEG-rich phase. Using optical tweezers to drag a bacterium across the DEX/PEG
interface, we demonstrate that the overall bacteria partitioning is determined by a competition
between the interfacial forces and bacterial propulsive forces. Our measurements are supported
by a theoretical model of dilute active rods embedded within a periodic soft confinement
potential.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2405.08995
- 저널판: https://doi.org/10.1103/6gm5-cnv1
- arXiv comment: *5 pages, 3 figures, See Supporting Information in ancillary files (1 appendix, 3 movies)*

## 요약

덱스트란(DEX)/PEG **수성 이상계(ATPS)** 에서 운동성 박테리아가 어떻게 분배되는지 본다.

- **비운동성** 박테리아는 모든 조건에서 **DEX-rich 상에만** 들어간다
- **운동성** 박테리아는 부드러운 DEX/PEG 계면을 **뚫고 지나가** 두 상에 나뉘어 분포한다

*B. subtilis* 기준, DEX 조성을 이상영역 안에서 늘리면 DEX-rich 상의 운동성 박테리아 분율이
**0.58 → 1** 로 증가한다.

가설: DEX와 세포벽 사이 화학적 친화도가 박테리아를 DEX-rich 상에 **약하게 구속**하는데,
운동성이 그 구속을 이길 만큼의 기계적 힘을 만든다는 것.

## 실험·시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 계 | DEX / PEG 수성 이상계 |
| 생물 | *Bacillus subtilis* (평균 몸길이 ≈ 5 µm) |
| 챔버 두께 | **5 µm** — 박테리아가 수직으로 서지 못하게 |
| 계면장력 | **`~10⁻⁵ N/m`** (매우 낮음) |
| 측정 분율 | DEX-rich 상 운동성 박테리아 분율 **0.58 → 1** |
| 대조 | 비운동성 박테리아 = 1.0 (전량 DEX-rich) |
| 보조 측정 | 광집게로 계면 통과 힘 |

> ⚠️ BD 쪽 `Δt`·`N`·퍼텐셜은 SI에서 확인 필요 — 이번 스캔에서는 잡히지 않았다. **미확보.**

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### ★ "운동성 vs 구속" 경쟁 — 우리 무차원 원장에 넣을 축

`kT` 단독으로는 마이크론 입자를 계면에 가둔다 (`계면 에너지 ≫ kT`). 그런데 **자기추진력이
그 장벽을 넘는다.** 이건 다음 무차원수로 표현된다:

```
자기추진력 / 구속력  =  ζU₀ / (계면장력 × 특성길이)
```

선행 프로젝트의 `_derived_quantities()`에 이미 유사한 항목이 있다 —
`active_escape_ratio_propulsion_over_trap = γv₀/(k·cutoff)`, *">1이면 직진하는 ABP가 트랩을
이기고 탈출"*. **같은 구조의 무차원수**다.

> `03_units_nondim.md`에 **"능동력 대 구속력 비"** 를 일반 항목으로 넣을 것. 트랩이든
> 계면이든 같은 형태다.

### 벤치마크 후보 (조건부)

- `cheon_motile_partition_fraction` — DEX 조성 대비 운동성 박테리아의 DEX-rich 분율 (0.58 → 1)
- **단, 우리 v1로 재현하려면** 두 상 사이 화학 친화도 + 계면을 모델링해야 하고, 박테리아는
  구형 ABP가 아니다. **현재 범위 밖.**

### 참고 — 사용자 공동제1저자

`Kyu Hwan Choi` 공동제1저자, `Kevin J. Modica` 공저. SI에 BD 세부가 있을 가능성이 높으니
`raw/lab/si/2024-cheon-motility-partitioning-atps__SI.pdf` 를 직접 확인할 것.

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
