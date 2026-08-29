---
type: finding
author: agent
drafted: 2026-07-27
confirmed_by:
question: "Kyu Hwan Choi의 전체 저작이 이 BD 에이전트에 무엇을 줄 수 있는가?"
answer: "두 가지. ① 학위논문이 완전한 BD 프로토콜을 준다. ② 랩 이전 광집게 연구가 쌍퍼텐셜의 실측 대조군을 준다 — 증거층 ④."
source: "Google Scholar (FqLTe9AAAAAJ), 2026-07-27 기준 30편 + Show more"
cites:
  - knowledge/source/papers/2025-choi-phd-thesis-noneq-bioinspired-soft-matter.md
  - knowledge/source/papers/2023-lee-colloidal-debye-force.md
  - knowledge/source/papers/2019-choi-electrostatic-pmma-optical-tweezers.md
related_systems:
  - passive-sphere--equilibrium-structure
---

# Choi 저작 전체 — 이 에이전트에 무엇을 주는가

**Google Scholar 기준 30편** (2026-07-27, "Show more" 남아 있음 · 인용 256 · h-index 11).
경력이 두 시기로 뚜렷이 갈리고, **각 시기가 우리에게 주는 것이 다르다.**

| 시기 | 소속 | 주제 | 우리에게 주는 것 |
|---|---|---|---|
| **2017–2023** | 경희대 (Park BJ 그룹) | **광집게 콜로이드 쌍상호작용 측정** · 계면 포집 입자 | ★★ **증거층 ④ — 쌍퍼텐셜 실측 대조군** |
| **2023–2026** | UCSB → Stanford (Takatori) | 능동물질 · 브러시 콜로이드 · 생체모사 | ★ BD 프로토콜 · 계-동역학 카드 |
| **2025** | UCSB 학위논문 | 위 4편을 묶은 148쪽 | ★★★ **완전한 BD 프로토콜** |

---

## 1. ★★★ 학위논문이 가장 값지다

[`2025-choi-phd-thesis`](../../source/papers/2025-choi-phd-thesis-noneq-bioinspired-soft-matter.md) — 148쪽, eScholarship 오픈액세스.

**논문 4편을 덮으면서, 논문에는 없는 "Theoretical Framework"와 "Experimental Design" 절이 있다.**
우리가 paywall 때문에 못 구하던 것이 대부분 여기 있었다.

Ch.2에서 건진 것:

| 항목 | 값 | 왜 중요한가 |
|---|---|---|
| **`Δt = 2×10⁻⁵ (σ²/D)`** | 변위 **0.0063σ** | `Δt` 표본의 **독립 확인** |
| **`ζ = 3πη·d`** (지름 기준!) | | **규약 함정의 실제 사례** — `6πηa`와 같은 값이지만 코드에서 2배 틀림 |
| 입자 수 · 실현 | 600–2000 · **20–30회** | `D16` 오차 산출 — 랩은 독립 시드 앙상블 |
| 접근-이완 프로토콜 | `ΔH`/`Δt_incr` 비로 속도 고정 | 시간의존 계의 재사용 패턴 |
| 닫힘 유효구간 | **`φ < 30%`** 명시 | `S4 LIT-GROUND`가 대조할 수 있는 명시적 범위 |

> **교훈:** 랩 논문이 paywall이면 **학위논문을 먼저 찾아라.** 오픈액세스이고, 방법이 더 자세하고,
> 여러 논문을 한 번에 덮는다. 다음 단계에서 다른 랩 구성원 학위논문도 찾아볼 가치가 있다.

---

## 2. ★★ 랩 이전 광집게 연구 — 우리가 몰랐던 자산

**이 발견이 이번 조사의 핵심이다.** Takatori 랩 이전 연구가 **콜로이드 쌍퍼텐셜의 직접 측정**이다.
그건 우리 에이전트가 `md.pair.Yukawa` / `md.pair.DLVO`로 **넣을 바로 그 상호작용**의 실측값이다.

| 논문 | 무엇을 쟀나 | 원문 |
|---|---|---|
| [Lee 2023 Nat Commun](../../source/papers/2023-lee-colloidal-debye-force.md) | **콜로이드 Debye 힘 직접 측정** | ✅ 확보 (CC BY) |
| [Choi 2019 Soft Matter](../../source/papers/2019-choi-electrostatic-pmma-optical-tweezers.md) | PMMA 미소구 **정전기 상호작용** (제1저자) | 🔒 |
| [Kang 2019 JPCL](../../source/papers/2019-kang-mapping-anisotropic-colloidal-interactions.md) | **이방성·불균질** 상호작용 맵핑 | 🔒 |
| [Choi 2020 ACS APM](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md) | 표면 불균질과 자기퍼텐셜 (제1저자) | 🔒 |

### 왜 이게 증거층 ④인가

`master_plan.md` §6의 4층 증거 체계에서 **④ 독립 방법**은 "다른 경로로 같은 답이 나오는가"다.
쌍퍼텐셜에 대해 이보다 독립적인 경로는 없다 — **시뮬레이션이 가정하는 것을 실험이 직접 잰다.**

```
우리 BD:   U(r) = DLVO/Yukawa 를 가정하고 넣는다
이 논문들: U(r) 을 광집게로 직접 잰다
```

**그래서 새 계-동역학 카드가 필요하다:** `charged-colloid--pair-potential`.
목적 동역학이 "구조"도 "수송"도 아니라 **"쌍퍼텐셜 자체의 검증"** 인 쌍이다.
`κ⁻¹`, `κσ`, `ζ`-전위, 이온세기가 주요 파라미터가 된다.

> ⚠️ **단, 4편 중 3편이 paywall이다.** Debye 힘 논문(오픈)만으로 시작하고,
> 나머지는 기관 접근권으로 `raw/`에 넣어야 증류할 수 있다.

### 곁가지 — 계면 포집 입자 (`adjacent`)

모세관력·타원체·렌즈형·흡착 확률 6편. v1(구형 3D 벌크)과는 거리가 있지만,
**비구형 입자**나 **계면 계**로 범위가 넓어지면 그때 볼 것.

---

## 3. 관련 없는 것 — 명시해둔다

Park 그룹 시절의 재료·촉매 논문 6편(두부 유래 탄소 ORR, 메탄-메탄올 전환, 바이오연료,
Pickering 안정제, 미세먼지 에멀전, 폴리모폼 입자)은 **이 에이전트와 무관하다.**
수집 목록에는 남기되 `source/`에 항목을 만들지 않았다.

> **안 만든 것을 적어두는 이유:** 나중에 "왜 30편 중 24편만 있나"를 다시 묻지 않기 위해서다.

학회 초록 4건(AIChE 2025 ×2, APS 2025 ×2)도 항목을 만들지 않았다 — 초록은 방법이 없다.
다만 **APS 2025 "DNA Liquid Droplets Embedded in Cytoskeletal Networks"** (Wilken, Choi, Hopkins,
Marchetti, Dogic, Takatori)는 향후 논문이 될 가능성이 있어 여기 기록해둔다.

---

## 4. 수집 현황

| | |
|---|---|
| Scholar 등재 | **30편** (+ "Show more" 미확인분) |
| `source/` 항목 생성 | **12편** (학위논문 1 + 광집게·계면 11) |
| 이미 있던 Takatori 논문 | 6편 |
| 원문 확보 | 학위논문 ✅ · Debye 힘 ✅ |
| paywall | 9편 |
| 의도적 제외 | 재료·촉매 6편 · 학회초록 4건 |

**전체 `source/papers/` 38편**이 됐다.

### 확보 실패와 그 이유

| 대상 | 결과 |
|---|---|
| Google Scholar 31편 이후 | 페이지네이션 차단 (`cstart=30`이 빈 페이지 반환) |
| RSC (ChemComm) | 스크립트 다운로드 차단 |
| Elsevier (JIEC) | 리다이렉트 차단 |
| PMC (Barakat 2024) | 봇 차단 — 오픈액세스인데도 |

> Crossref가 ACS 논문 5편을 **보충자료 DOI(`.sNNN`)** 에 잘못 매칭했다. 접미사를 벗겨 교정했다.
> 자동 수집 스크립트를 짤 때 같은 함정이 있다.

---

## 5. 이 조사가 바꾸는 것

| 문서 | 내용 |
|---|---|
| `03_units_nondim.md` | ★ **`ζ = 3πηd` vs `6πηa` 규약을 `UnitMap` 필드로 강제** — 학위논문이 실제 사례 |
| [[dt-gate-should-be-displacement-based]] | `Δt = 2e-5 τ_D` 독립 확인. 표본 3건 유지 |
| `systems/` | **신규 카드 필요:** `charged-colloid--pair-potential`, `brush-colloid--nonequilibrium-contact` |
| `07_observables.md` (`D16`) | 랩 관행 = **20–30 독립 실현**. block averaging과 병행 |
| `benchmarks/` | Debye 힘 논문에서 **측정 `U(r)`** 을 추출하면 DLVO 구현 검증 오라클이 된다 |

## 남은 불확실성

- **Scholar 31편 이후를 못 봤다.** h-index 11에 30편이면 대체로 다 나온 것으로 보이지만 확언 못 함
- 학위논문 **Ch.3–5를 아직 안 읽었다** — ATPS 챕터에 BD 파라미터가 더 있을 가능성
- paywall 9편 중 **광집게 3편이 가장 아깝다** — 증거층 ④의 핵심
- 광집게 측정은 **2입자 계**다. 다체 BD의 검증 오라클로 쓸 때 그 한계를 명시해야 한다
