---
type: source
kind: arxiv
lab_authored: true
title: Dynamic interfaces for contact-time control of colloidal interactions
authors:
  - "Xu Y"
  - "Choi KH"
  - "Nagella SG"
  - "Takatori SC"
year: 2023
journal: Soft Matter
doi: 10.1039/D3SM00673E
arxiv_id: 2303.08880
source_url: "https://arxiv.org/abs/2303.08880"
raw_file: knowledge/raw/lab/2023-xu-dynamic-interfaces-contact-time.pdf
si_available: true
access: open (arXiv)
engine: "HOOMD-blue (BD, GPU)"
reproduced: no
parameters_extracted: yes
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---
# Dynamic interfaces for contact-time control of colloidal interactions

## 왜 이 위키에 있는가

**접촉시간 의존 콜로이드 상호작용**. 시간의존 퍼텐셜을 BD로 다룬 사례.

**관련도:** `direct`

## 초록 (원문)

Understanding multibody interactions between colloidal particles out of equilibrium has a
profound impact on dynamical processes such as colloidal self assembly. However, traditional
colloidal interactions are effectively quasi-static on colloidal timescales and cannot be
modulated out of equilibrium. A mechanism to dynamically tune the interactions during colloidal
contacts can provide new avenues for self assembly and material design. In this work, we develop
a framework based on polymer-coated colloids and demonstrate that in-plane surface mobility and
mechanical relaxation of polymers at colloidal contact interfaces enable an effective, dynamic
interaction. Combining analytical theory, simulations, and optical tweezer experiments, we
demonstrate precise control of dynamic pair interactions over a range of pico-Newton forces and
seconds timescales. Our model may be used to engineer colloids with exquisite control over the
kinetics and thermodynamics of colloidal self-assembly dynamics via interface modulation and
nonequilibrium processing.

## 원문 접근

- arXiv 전문: https://arxiv.org/abs/2303.08880
- 저널판: https://doi.org/10.1039/D3SM00673E
- arXiv comment: *6 pages, 4 figures. See Supporting Information in ancillary files (1 appendix, 8 movies)*

## 요약

표면에 **이동 가능한 폴리머(브러시)** 를 붙인 콜로이드 두 개를 서로 접근시켰다가 고정하고,
접촉 시간에 따라 상호작용이 어떻게 달라지는지 본다. 폴리머가 접촉면에서 빠져나가는 데 시간이
걸리므로 **상호작용이 이력(history)에 의존**한다.

콜로이드 코어는 결정론적으로 움직이거나(접근 단계) 고정되고(완화 단계), **폴리머 비드만 BD로
움직인다.** 비평형 접근 → 완화 프로토콜이 핵심.

## 시뮬레이션 파라미터

| 항목 | 값 |
|---|---|
| 엔진 | **HOOMD-blue** (GPU 가속) |
| 적분 | 과감쇠 Langevin, `Δxᵢ/Δt = (Fᴾᵢ + Fᴿᵢ)/ζ` |
| 항력계수 | `ζ = 3πη d_ρ` → `D_ρ = kT/ζ` (Stokes–Einstein) |
| **시간간격** | **`Δt = 2×10⁻⁵ · (σ²/D_ρ)`** |
| 샘플링 | 10³ 스텝마다 |
| 앙상블 | **20 realization 평균** |
| 길이 단위 | `σ = 1 µm` |
| 박스 | `L_y = L_z = 33σ`, `L_x = 65σ`, 주기 경계 |
| 폴리머 | Kremer–Grest bead-spring, `d_ρ = 0.8 µm`, **17 beads/chain** |
| 비결합 상호작용 | **WCA**, `ε = α = k_BT` |
| 결합 | **FENE** `r₀ = 1.5σ`, `k_FENE = 30 k_BT` |
| 각도 | `V_ang = (1−cos(θ−π))·l_ρ/σ`, 지속길이 `l_ρ = 13σ` |
| 표면 피복률 | `φ = n_ρ(d_ρ²/4d_c²) = 0.15` (비교 시 0.15–0.43) |
| 그래프팅 | 두 동심 구면 벽 사이에 비드 구속 — 표면 위 이동은 가능, 이탈 불가 |
| 평형화 | `H = 30σ` (브러시 비중첩)에서 사전 평형화 |
| 접근 | `ΔH = 0.005σ–0.08σ` / `Δt_incr = 10²–10³` 스텝 |
| 완화 | `10⁶–10⁷` 스텝, 콜로이드 고정 |

## 참고사항 — 우리 에이전트에 어떻게 쓰는가

### ★ `Δt` 규약이 우리 게이트와 일치한다

**`Δt = 2×10⁻⁵ σ²/D`** 는 정확히 우리가 `03_units_nondim.md`에서 쓰려는 형태 — **브라운 확산시간
`τ_D = σ²/D`로 무차원화한 시간간격**이다. `master_plan.md` §5의 S3 게이트가 `Δt/τ_D ≤ 1e-4`인데
이 논문은 `2×10⁻⁵`로 **게이트 안쪽**이다.

> 우리 게이트가 실제 랩 관행과 어긋나지 않는다는 **첫 외부 확인**이다. `D7`(기준 시간 단위 =
> `τ_D`)의 근거로 쓸 수 있다.

### 그 외 참고

- **WCA를 `ε = k_BT`로 썼다.** 선행 프로젝트는 WCA의 `r⁻¹³` 코어가 오버댐프에서 위험하다고
  경고했는데(`forces.py:117`), 여기서는 `Δt`가 충분히 작아서 문제가 안 된 것으로 보인다.
  **`D8`(배제부피 기본 퍼텐셜) 판단에 중요한 자료** — "WCA가 무조건 위험한 게 아니라 `Δt`와
  묶인 문제"라는 방향
- **20 realization 평균** — 우리 `D16`(오차막대) 결정과 비교할 것. 독립 시드 앙상블 방식이다
- `ζ = 3πη d_ρ` — 반지름이 아니라 **지름** 기준. 항력계수 관례가 그룹마다 다르니 `UnitMap`에서
  주의할 지점

### 벤치마크로 쓰기는 어렵다

비평형 프로토콜(접근 속도 의존)이라 정상상태 관측량이 없다. **파라미터 사전으로서 값지고
벤치마크로는 부적합.**

---

> **재현 상태: `reproduced: no`.** 위 파라미터는 논문에 적힌 값을 옮긴 것이고,
> 우리 코드에서 돌려본 것이 아니다. 검증 근거로 인용하지 않는다 — 계약: [`knowledge/wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
