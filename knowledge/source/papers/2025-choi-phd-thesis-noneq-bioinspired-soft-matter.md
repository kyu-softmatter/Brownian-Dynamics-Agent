---
type: source
kind: book
subtype: dissertation
lab_authored: true
title: "Non-equilibrium Dynamics in Bio-inspired Soft Matter (박사학위논문)"
authors:
  - "Choi KH"
year: 2025
journal: "UC Santa Barbara — Electronic Theses and Dissertations"
advisor: "Takatori SC"
source_url: "https://escholarship.org/uc/item/10f9g565"
raw_file: knowledge/raw/lab/2025-choi-phd-thesis-noneq-bioinspired-soft-matter.pdf
pages: 148
access: open (eScholarship)
engine: "HOOMD-blue (BD, GPU)"
reproduced: no
parameters_extracted: yes
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---

# Non-equilibrium Dynamics in Bio-inspired Soft Matter — 박사학위논문

**148쪽, 5개 챕터.** 저자의 UCSB 학위논문(지도교수 Takatori). eScholarship 오픈액세스.

> ★ **이 문서 하나가 논문 4편을 덮는다.** 각 챕터가 발표 논문에 대응하면서,
> **논문에는 없는 "Theoretical and Experimental Framework"와 "Experimental Design" 절**을 갖고 있다.
> 우리가 원하던 BD 프로토콜이 여기 있다.

## 챕터 구성

| Ch | 제목 | 대응 논문 |
|---|---|---|
| 1 | Introduction | — |
| 2 | **Dynamic Interfaces for Contact-time Control of Colloidal Interactions** | [Xu 2023](2023-xu-dynamic-interfaces-contact-time.md) |
| 3 | Anisotropic Diffusion Mechanism in Lyotropic Chromonic Liquid Crystals | [Choi 2026 FRAP](2026-choi-anisotropic-diffusion-lclc-frap.md) |
| 4 | Mechanics of Fibrous Networks Assembled by an Active Fluid | [Choi 2026 fiber](2026-choi-heterogeneous-fiber-networks.md) |
| 5 | Motility Modulates the Partitioning of Bacteria in ATPS | [Cheon 2024](2024-cheon-motility-partitioning-atps.md) |

---

## ★ Ch.2 §2.2.2 — 완전한 BD 프로토콜

논문 본문에 없던 것이 전부 있다.

### 지배 방정식

과감쇠 Langevin. 상호작용력 `F^P`와 열적 힘 `F^R`(요동-소산 만족)을 항력 `ζ`로 나눈 것:

```
Δx_i / Δt = ( F^P_i + F^R_i ) / ζ
```

### ★ 항력 규약 — 지름 기준이다

```
ζ = 3πη·d_ρ            ← d_ρ 는 비드 "지름"
```

> **이것이 [`passive-sphere--equilibrium-structure`](../../wiki/systems/passive-sphere--equilibrium-structure.md)
> 카드에서 경고한 규약 함정의 실제 사례다.** 교과서 관례 `6πηa`(반지름)와 **같은 값**이지만,
> 코드에서 `a` 자리에 지름을 넣으면 **2배 틀린다.**
> `UnitMap`은 어느 규약인지 필드로 기록해야 한다.

### 수치 파라미터 — `Δt` 독립 확인 ★

| 항목 | 값 |
|---|---|
| **시간 간격** | **`Δt = 2×10⁻⁵ · (σ²/D_ρ)`** |
| 샘플링 | **10³ 스텝마다** |
| 입자 수 | **600–2000** 조대화 폴리머 비드 |
| 독립 실현 | **20–30회** |
| 엔진 | HOOMD-blue, GPU 가속 |

> **`Δt = 2e-5 τ_D`가 학위논문에서 독립적으로 확인됐다.** Xu 2023 SI에서 읽은 값과 일치한다.
> 스텝당 확산 변위로 환산하면 **0.0063σ** — [[dt-gate-should-be-displacement-based]]의 표본 3건 중
> 가장 보수적인 값.

### 폴리머 모델

- **Kremer–Grest bead-spring** + 반유연성(semi-flexibility)
- 그래프팅 지점이 표면을 따라 **확산 이동** 가능 — 이게 "dynamic interface"의 핵심
- 유효 콜로이드 상호작용: `F = −n_ρ M ∂_H V`, 중심선 방향 힘을 평균

### ★ 접근-이완 프로토콜 (재사용 가치 높음)

시간의존 상호작용을 재려면 **비평형 접근 과정 자체를 제어**해야 한다.

| 단계 | 설정 |
|---|---|
| 초기 분리 | `H = 30σ` (브러시 비중첩) |
| 접근 | `ΔH = 0.005σ – 0.08σ` 씩, `Δt_incr = 10²–10³` 스텝마다 |
| 접근 속도 | `ΔH`와 `Δt_incr`을 조절해 `v` 고정 |
| 이완 | `H` 고정 후 **`10⁶–10⁷` 스텝** |
| 힘 측정 | `10²–10³` 스텝마다 |

> **속도를 "설정"하는 게 아니라 증분과 간격의 비로 만든다.** 우리 에이전트가 시간의존
> 프로토콜을 다룰 때 그대로 쓸 수 있는 패턴.

### 계 파라미터

| 항목 | 값 |
|---|---|
| F-actin 면적분율 | `φ = 0.11 – 0.46` |
| F-actin 단량체 크기 | `d = 7 nm` |
| 브러시 평균 높이 | `h₀ ≈ 9 µm` (시뮬레이션으로 설정) |
| 광집게 강성 | `κ_t = 0.5 – 0.7 pN/µm` |
| 접근 속도 (실험) | `v = 0.5 – 10 µm/s` |

### 이론 골격

- **Smoluchowski 방정식** — 단량체 확률분포 `P_N(h₁..h_N, t; r)`
- 유체역학 상호작용 **없음** → `D_ij = I_ij kT/ζ_i` (Stokes–Einstein–Sutherland)
- BBGKY 유사 닫힘. **반희박(semidilute, 부피분율 30% 미만)에서 유효**하다고 명시
- 평형 단량체 밀도를 가우시안으로 근사: `n_eq(z) ≈ n_ρ exp(−(z²−h₀²)/h₀²)`
- 평균장 쌍퍼텐셜: `V_single(z) = −kT ln(n_eq(z)/n_eq(0))`

> **닫힘의 유효 범위가 명시돼 있다** (`φ < 30%`). 우리 `S4 LIT-GROUND`가 문헌 범위를 대조할 때
> 이런 명시적 유효구간이 정확히 필요한 것이다.

---

## 우리 에이전트에 반영할 것

| 반영처 | 내용 |
|---|---|
| `03_units_nondim.md` | **`ζ = 3πηd` vs `6πηa` 규약을 `UnitMap` 필드로 강제** |
| [[dt-gate-should-be-displacement-based]] | `Δt = 2e-5 τ_D` 표본 독립 확인 (변위 0.0063σ) |
| `06_repair_policy.md` | 접근-이완 프로토콜 — 시간의존 계의 표준 패턴 |
| `systems/` 신규 카드 | `brush-colloid--nonequilibrium-contact` — 이 챕터가 근거 |
| `07_observables.md` | **20–30 독립 실현** — 랩 관행. block averaging과 병행 (`D16`) |

## 남은 것

- **Ch.3–5는 아직 안 읽었다.** LCLC 이방성 확산(D 텐서), 섬유망 역학, ATPS 분배
- Ch.5 ATPS는 BD를 쓰므로 파라미터가 더 있을 가능성
- `reproduced: no` — 우리 환경에서 돌린 것이 아니다
