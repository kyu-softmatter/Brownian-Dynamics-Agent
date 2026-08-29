---
type: source
kind: arxiv
lab_authored: true
title: Dynamic surfactants drive anisotropic colloidal assembly
authors:
  - "Xu Y"
  - "Jandhyala P"
  - "Takatori SC"
year: 2024
journal: J. Chem. Phys. 161, 064901
doi: 10.1063/5.0220112
arxiv_id: 2405.08936
source_url: "https://arxiv.org/abs/2405.08936"
raw_file: knowledge/raw/lab/2024-xu-dynamic-surfactants-anisotropic-assembly.pdf
si_available: false
access: open (arXiv)
engine: "HOOMD-blue (BD, GPU)"
reproduced: no
parameters_extracted: partial
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---
# Dynamic surfactants drive anisotropic colloidal assembly

## 왜 이 위키에 있는가

**비등방 콜로이드 조립의 BD**. 시간의존 상호작용을 다룸 — 퍼텐셜 설계 참조.

**관련도:** `direct`

## 초록 (원문)

Colloidal building blocks with re-configurable shapes and dynamic interactions can exhibit
unusual self-assembly behaviors and pathways. In this work, we consider the phase behavior of
colloids coated with surface-mobile polymer brushes that behave as "dynamic surfactants." Unlike
traditional polymer-grafted colloids, we show that colloids coated with dynamic surfactants can
acquire anisotropic macroscopic assemblies, even for spherical colloids with isotropic
attractive interactions. We use Brownian Dynamics simulations and dynamic density functional
theory (DDFT) to demonstrate that time-dependent reorganization of the dynamic surfactants leads
to phase diagrams with anisotropic assemblies. We observed that the microscopic polymer
distributions impose unique geometric constraints between colloids that control their packing
into lamellar, string, and vesicle phases. Our work may help discover versatile building blocks
and provide extensive design freedom for assembly out of thermodynamic equilibrium.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2405.08936
- 저널판: https://doi.org/10.1063/5.0220112
- arXiv comment: *8 pages, 4 figures*

## 요약

표면 이동성 폴리머(“동적 계면활성제”)로 덮인 콜로이드가 **이방성 조립**을 하는 과정.
폴리머가 접촉면에서 밀려나면서 유효 인력에 방향성이 생기고, 결과적으로 **string phase → lamellar
(이중층 시트가 쌓인) phase** 로 진행한다.

BD 시뮬레이션 + Smoluchowski 장이론 + DDFT를 조합해 **스피노달 곡선**(환산온도 × 유효 패킹분율 ×
표면 피복률의 3차원 상도)을 냈다.

## 시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 엔진 | **HOOMD-blue** (GPU 가속) |
| 콜로이드 수 | **N_c = 1,000–1,200** |
| **총 스텝** | **8×10⁸** |
| 박스 | 3D `L³`, **`L = 20σ – 60σ`** (φ 맞추려고 변화) |
| 특성 시간 | **`t* = 4R_c²/D_c`** — 콜로이드가 지름만큼 확산하는 브라운 시간 |
| 배제부피 | **WCA** ("hard sphere-like") |
| 패킹분율 | `φ = N_c·(π/3)·4R_c³/L³`.  대표값 **`φ = 0.25`** |
| 유효 패킹분율 | `φ_eff = φ + N_c·m·n·πd_p³/(6L³)` — 폴리머 비드 부피 보정 |
| 표면 피복률 | `f = 0.05` |
| 사슬 길이 | `m = 2` |
| 환산온도 | `k_BT·R_c/a = 0.24` |
| 관측 시점 | `t/t* ≤ 500` 에서 string phase, 이후 lamellar |

> ⚠️ **`Δt` 값이 본문에 없다.** 총 스텝(8×10⁸)과 `t*`만 주어져 있어 `Δt/t*`를 역산할 수 없다.

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### 참고사항

- **`t* = 4R_c²/D_c`** 를 특성 시간으로 쓴다. 반지름 `R_c` 기준이라 `4R_c² = σ²` — 즉
  **우리의 `τ_D = σ²/D`와 같은 양**이다. 표기만 다르다. `UnitMap`에서 이 환산을 명시해둘 것
- **`φ_eff` 보정** — 폴리머 비드가 차지하는 부피를 더해 유효 패킹분율을 다시 정의한다.
  복합 입자를 다룰 때 `03_units_nondim.md`의 φ 정의가 모호해질 수 있다는 신호
- **8×10⁸ 스텝은 CUDA 없는 M4에서 비현실적이다.** `master_plan.md` §8의 비용 게이트가
  이런 런을 걸러내야 한다. 재현하려면 클러스터 또는 축소 스케일이 필요
- 상도(스피노달)가 있으므로 **장기적으로는 벤치마크 후보** — 다만 `Δt` 미확보라 지금은 불가

### 벤치마크 후보 (조건부)

- `xu_lamellar_spinodal` — 주어진 `(f, m, k_BTR_c/a)`에서 스피노달 φ. **`Δt` 확보 후에만.**

---

> **재현 상태: `reproduced: no`.** 위 파라미터는 논문에 적힌 값을 옮긴 것이고,
> 우리 코드에서 돌려본 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
