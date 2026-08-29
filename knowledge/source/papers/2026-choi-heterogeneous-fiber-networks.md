---
type: source
kind: arxiv
lab_authored: true
title: Mechanics of heterogeneous fiber networks
authors:
  - "Choi KH"
  - "Ray S"
  - "Sweeney R"
  - "Dogic Z"
  - "Takatori SC"
year: 2026
journal: submitted (arXiv preprint)
arxiv_id: 2605.11342
source_url: "https://arxiv.org/abs/2605.11342"
raw_file: knowledge/raw/lab/2026-choi-heterogeneous-fiber-networks.pdf
si_available: false
access: open (arXiv)
engine: "실험 (능동 미세유변학 + 형광 이미징)"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "adjacent"
---
# Mechanics of heterogeneous fiber networks

## 왜 이 위키에 있는가

이질적 섬유망의 역학. 네트워크 계 모델링 시 참조.

**관련도:** `adjacent`

## 초록 (원문)

Internally generated active stresses drive soft materials into architectures inaccessible to
thermal self-assembly. We use a microtubule-based active fluid to assemble and irreversibly
restructure actin-fascin networks. Subsequently, we probe the mesoscale mechanics of such
networks by combining active microrheology with fluorescence imaging of the strain field around
the probe. Increasing motor concentration broadens the pore-size distribution and thickens load-
bearing bundles, raising the mean local elastic modulus and its spatial variability.
Displacement fields of actively-processed networks propagate over longer range when compared to
unprocessed networks. At large strains, both networks strain soften and plastically restructure.
The combined microrheology and strain-imaging approach show that tunable active stresses
reprogram the structure and viscoelastic response of fiber networks at the scale of their
structural heterogeneity.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2605.11342
- 저널판: —
- arXiv comment: *7 pages, 4 figures, 4 videos*

## 요약

미세소관 기반 **능동 유체**로 actin–fascin 네트워크를 조립하고 비가역적으로 재구조화한 뒤,
그 **중간스케일(mesoscale) 역학**을 능동 미세유변학 + 탐침 주변 변형장 이미징으로 측정한다.

모터 농도를 올리면 **공극 크기 분포가 넓어지고 하중지지 다발이 굵어져서**, 평균 국소 탄성률과
그 공간적 편차가 함께 커진다. 능동 처리된 네트워크는 변위장이 더 멀리 전파된다. 큰 변형에서는
두 네트워크 모두 **strain softening과 소성 재구조화**를 보인다.

## 실험·시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 능동 유체 | 형광 미세소관 (25 µm, 3% Alexa 647 표지) |
| 모터 | K401-streptavidin 클러스터 (KSA), 25–250 nm |
| **KSA 농도 시리즈** | **190, 140, 90, 50 nM** |
| 네트워크 | actin + fascin |
| Fascin 농도 | 2–4 µM |
| Actin 표지 | phalloidin-Alexa 488, 670 nM (고정) |
| 탐침 비드 | 폴리스티렌 **2R = 5 µm**, gelsolin 부착 |
| ATP 유사체 | AMP-PNP (비가수분해성) |
| 측정 | 능동 미세유변학 + 수동 2점 미세유변학 |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 우리 v1 범위 밖 — 하지만 방법론 참고

섬유 네트워크는 우리 v1(구형 콜로이드 3D, `D3`)에 없다. 다만 두 가지가 참고된다.

- **공극 크기 분포**를 관측량으로 쓴다 — 우리 `07_observables.md`의 구조 관측량 목록에
  Voronoi·클러스터는 있지만 **공극(pore) 분석은 없다.** 다공성 계를 다루면 필요
- **국소 탄성률의 공간 편차**를 보고한다 — 평균만이 아니라 분산을 관측량으로 삼는 태도.
  우리 "오차막대 없는 숫자 금지"(`A2`)와 같은 정신

### 참고 — 사용자 본인 공동제1저자

Dogic 랩과의 공동연구. 능동 물질 + 네트워크 역학. 미발표(arXiv 프리프린트, 2026-05).

---

> **재현 상태: `reproduced: no`.** 위 값은 논문에 적힌 것을 옮긴 것이고,
> 우리 코드에서 확인한 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
