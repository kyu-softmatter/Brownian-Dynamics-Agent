---
type: source
kind: paper
lab_authored: true
title: Mapping anisotropic and heterogeneous colloidal interactions via optical laser tweezers
authors:
  - "Kang DW"
  - "Choi KH"
  - "Lee SJ"
  - "Park BJ"
year: 2019
journal: J. Phys. Chem. Lett. 10(8), 1691-1697
doi: 10.1021/acs.jpclett.9b00232
source_url: "https://doi.org/10.1021/acs.jpclett.9b00232"
access: paywall
engine: 실험 (광집게)
reproduced: no
parameters_extracted: no
affiliation: Kyung Hee University (Park BJ group) — Takatori 랩 이전
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "method"
---
# Mapping anisotropic and heterogeneous colloidal interactions via optical laser tweezers

## 왜 이 위키에 있는가

**이방성·불균질** 쌍상호작용 맵핑. 표면 전하가 균일하지 않을 때 유효 퍼텐셜이 어떻게 갈리는지.

**관련도:** `method` · **소속:** Park BJ 그룹 (경희대) — Takatori 랩 합류 이전 연구

## 원문 접근

- DOI: https://doi.org/10.1021/acs.jpclett.9b00232
- ⚠️ **paywall.** 기관 접근권으로 받아 `knowledge/raw/lab/2019-kang-mapping-anisotropic-colloidal-interactions.pdf` 에 두면 증류 가능

## 계면 항력 — 광집게 셋업 출처 (η_eff 원출처는 아님)

`2023-lee-colloidal-debye-force`가 드래그 보정법의 출처로 **ref 24 = 이 논문**을 ref 23
([`2020-choi-electrostatic-self-potential-heterogeneity`](2020-choi-electrostatic-self-potential-heterogeneity.md))
과 **함께** 인용한다.

**추적 종료 (2026-07-27):** ref 23의 SI에 `κ_t = 3πdη_eff u/Δx` 형태까지는 있으나
`η_eff` 정의의 원출처는 그 SI의 ref 2 —
[Park, Pantina, Furst, Oettel, Reynaert, Vermant, *Langmuir* 2008, 24(5),
1686−1694](2008-park-salt-surfactant-interface-forces.md) 이고, **확보해서 확인 완료**다
(`η_eff = [η_oil(1−cosθ) + η_water(1+cosθ)]/2`).

이 논문은 광집게 셋업·시분할 트랩 쪽 출처이지 `η_eff` 원출처가 아니다. **확보 우선순위 하향** —
쌍상호작용 이방성/불균질 맵핑이 필요해질 때 다시 본다.

## 추출 대기

> ⚠️ **수집 기록**이다. 추출 전까지 근거로 인용하지 않는다.

측정 기반 논문에서 뽑을 것:

| 항목 | 값 | 위치 |
|---|---|---|
| 계 (입자 재질·크기·매질) | | |
| 물성 (η, T, 이온세기, ζ-전위) | | |
| **측정된 쌍퍼텐셜 `U(r)` 형태·크기** | | |
| Debye 길이 `κ⁻¹` / `κσ` | | |
| 불확실도 | | |
| **BD 퍼텐셜과 대조 가능한 수치** | | |

→ 대조 가능한 수치가 있으면 `knowledge/wiki/benchmarks/`에 등재하고
`systems/charged-colloid--pair-potential` 카드에 연결한다.
