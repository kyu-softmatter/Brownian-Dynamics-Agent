---
type: source
kind: paper
lab_authored: true
title: "Porous media microstructure determines the diffusion of active matter: experiments and simulations"
authors:
  - "Modica KJ"
  - "Xi Y"
  - "Takatori SC"
year: 2022
journal: Front. Phys. 10, 869175
doi: 10.3389/fphy.2022.869175
source_url: "https://www.frontiersin.org/articles/10.3389/fphy.2022.869175"
raw_file: knowledge/raw/lab/2022-modica-porous-media-active-diffusion.pdf
access: open
engine: "BD (ABP) + Janus 입자 실험"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---
# Porous media microstructure determines the diffusion of active matter: experiments and simulations

## 왜 이 위키에 있는가

**다공성 매질 내 능동물질 확산**. 실험과 BD를 같은 논문에서 대조 — 증거층 ③④ 동시 제공.

**관련도:** `direct`

## 원문 접근

- 오픈액세스 전문: https://www.frontiersin.org/articles/10.3389/fphy.2022.869175
- DOI: https://doi.org/10.3389/fphy.2022.869175

## 요약

고정 장애물 배열 속 **능동 브라운 입자(ABP)** 의 유효 확산계수를 다룬다. 능동 입자는
경계에 축적되는 성질 때문에, 같은 장애물 밀도에서도 **수동 입자보다 훨씬 크게 느려진다.**

핵심 수치: **장애물 면적분율 `φ ≈ 12%` 라는 희박한 조건에서도 능동 입자의 유효 확산계수가
25% 감소** — 같은 조건 수동 브라운 입자 감소폭의 약 2배.

Janus 입자 실험 + BD 시뮬레이션 + Smoluchowski 이론 세 갈래를 같은 계에 적용했다.

## 시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 모델 | **ABP** — 과감쇠 Langevin |
| 차원 | **2D** (실리카 비드가 챔버 바닥에 침강) |
| 능동 입자 상호작용 | **장애물과는 순수 척력, 입자끼리는 상호작용 없음** ("ideal gas") |
| 자기추진 속도 | **`U₀ = 0.84 ± 0.01 µm/s`** (실험에서 측정, 시뮬에 입력) |
| 재배향 시간 | **`τ_R = 14 ± 2 s`** (MSD 피팅으로 결정) |
| 장애물 | 5 µm 실리카, 지질이중층 코팅 |
| 능동 입자 | 4 µm 실리카 Janus (한쪽 Pt, 한쪽 지질이중층) |
| 장애물 면적분율 | **`φ ≈ 12%`** |
| 초기 배치 | **실험에서 측정한 장애물 좌표를 그대로 시뮬 입력으로 사용** |
| 연료 | 2% H₂O₂ |
| 경계조건 | 장애물 표면 무유속 `n̂·[U₀qP − D_T∇P] = 0`, 단위셀 주기 |

> ⚠️ `Δt`, 총 스텝, 입자 수는 본문에 없다 (Materials and Methods에 있을 것). **미확보.**

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### ★ 방법론이 우리 검증 철학과 같다

이 논문은 **"실험 좌표를 시뮬 입력으로 그대로 넣어서" 계를 1:1 재현**한다. 그리고 MSD를 비교해
차이가 나면 그것이 미고려 메커니즘(유체역학, 입자간 인력)의 증거라고 판정한다.

> 우리 `master_plan.md` §6의 **증거층 ④(독립 방법)** 을 실제로 수행한 사례다. 특히
> *"MSD 비교로 유체역학·인력이 원인이 아님을 배제했다"* 는 논리는 `D11`(HI 무시)의
> 정당화 논거로 쓸 수 있다.

### 벤치마크 후보 ★

- `modica_active_diffusion_reduction_phi12` — `φ = 12%` 고정 장애물 배열에서
  **ABP 유효 `D` 25% 감소**. 수동 입자 대비 약 2배.
  - 값이 명시적이고 (25%), 조건이 좁고 (`φ=12%`, `U₀`, `τ_R` 주어짐), 2D라 **싸다**
  - `U₀ = 0.84 µm/s`, `τ_R = 14 s` → `Pe = U₀τ_R/σ`를 우리 무차원 원장으로 환산해서 등재
- `modica_passive_vs_active_ratio` — 같은 `φ`에서 능동/수동 감소폭 비 ≈ 2

### 주의

- **2D다.** 우리 `D9` 기본값은 3D — 이 벤치마크를 쓰려면 2D 경로가 필요
- 장애물 좌표가 실험 유래라 **정확한 재현에는 원 데이터가 필요**. 무작위 배열로 대체하면
  통계적으로만 비교 가능

---

> **재현 상태: `reproduced: no`.** 위 파라미터는 논문에 적힌 값을 옮긴 것이고,
> 우리 코드에서 돌려본 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
