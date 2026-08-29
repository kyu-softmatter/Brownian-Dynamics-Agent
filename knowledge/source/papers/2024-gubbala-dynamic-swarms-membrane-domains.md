---
type: source
kind: arxiv
lab_authored: true
title: Dynamic swarms regulate the morphology and distribution of soft membrane domains
authors:
  - "Gubbala A"
  - "Arnold DP"
  - "Jena A"
  - "Anujarerat S"
  - "Takatori SC"
year: 2024
journal: Phys. Rev. E 110, 014410
doi: 10.1103/PhysRevE.110.014410
arxiv_id: 2402.06518
source_url: "https://arxiv.org/abs/2402.06518"
raw_file: knowledge/raw/lab/2024-gubbala-dynamic-swarms-membrane-domains.pdf
si_available: true
access: open (arXiv)
engine: "연속체 (Toner–Tu + Cahn–Hilliard)"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "adjacent"
---
# Dynamic swarms regulate the morphology and distribution of soft membrane domains

## 왜 이 위키에 있는가

능동 스웜이 막 도메인 형태를 조절. 조대화·클러스터 관측량 참조.

**관련도:** `adjacent`

## 초록 (원문)

We study the dynamic structure of lipid domain inclusions embedded within a phase-separated
reconstituted lipid bilayer in contact with a swarming flow of gliding filamentous actin.
Passive circular domains transition into highly-deformed morphologies that continuously
elongate, rotate, and pinch off into smaller fragments, leading to a dynamic steady state with
approximately 23x speed up in the relaxation of the intermediate scattering function compared to
passive membrane domains driven by purely thermal forces. To corroborate experimental results,
we develop a phase-field model of the lipid domains with two-way coupling to the Toner-Tu
equations. We report phase domains that become entrained in the chaotic eddy patterns, with
oscillating waves of domains that correlate with the dominant wavelengths of the Toner-Tu flow
fields.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2402.06518
- 저널판: https://doi.org/10.1103/PhysRevE.110.014410
- arXiv comment: *Main text is 6 pages with 4 figures. Supplemental Material includes a supplemental appendix with supplemental methods and theory, 9 supplemental figures, and 8 supplemental video captions. There are 8 supplemental videos as well*

## 요약

능동 네마틱 스웜이 상분리된 지질 도메인의 **형태와 분포**를 어떻게 조절하는지, **Toner–Tu
방정식과 Cahn–Hilliard phase field를 결합**해 본다.

앞선 연구(수축성 액토미오신)에서는 도메인 성장이 2× 가속됐는데, **혼돈적(chaotic) 능동 흐름**
아래에서는 다른 결과가 나온다 — 가속된 조대화 대신 **도메인이 섞이는 동적 정상상태**가 된다.

## 파라미터

| 항목 | 값 |
|---|---|
| 모형 | **Toner–Tu + Cahn–Hilliard** 결합 |
| 비교 대상 | 수축성 흐름 → 2× 가속 / 혼돈 흐름 → 동적 정상상태 |
| 관측량 | 도메인 형태, 분포 |
| 선행 이론 | nematohydrodynamic + Cahn–Hilliard |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 입자 기반이 아니다 — 직접 비교는 불가

연속체 모형이라 우리 BD와 1:1 대응이 안 된다. 다만:

- **능동 흐름의 성격(수축 vs 혼돈)에 따라 조대화 결과가 질적으로 달라진다**는 것은
  중요한 경고다. "능동이면 빨라진다"가 아니라 **흐름 구조에 의존**한다
- 우리가 능동 콜로이드의 상분리를 다룰 때, **조대화 지수 하나로 판정하면 안 된다**는 신호.
  동적 정상상태일 수도 있으므로 **정상상태 도달 여부를 먼저 확인**해야 한다

> `05_validation_gates.md`의 평형화 판정(`pymbar.detect_equilibration`)이
> **능동계에서는 "평형"이 아니라 "정상상태"** 를 봐야 한다는 점을 여기서 확인.
> 능동계는 열평형에 도달하지 않는다 — 게이트 문구를 정밀화할 것.

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
