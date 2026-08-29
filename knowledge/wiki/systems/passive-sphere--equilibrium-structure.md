---
type: system
author: agent
drafted: 2026-07-27
confirmed_by:
system: passive-sphere
dynamics: equilibrium-structure
status: draft
cites:
  - knowledge/wiki/benchmarks/candidates.md
  - knowledge/wiki/findings/wca-reproduces-carnahan-starling.md
  - knowledge/wiki/findings/dt-gate-should-be-displacement-based.md
---

# `수동 구형 콜로이드` × `평형 구조`

**묻는 질문:** 주어진 `φ`와 상호작용에서 평형 미세구조는 무엇인가 — `g(r)`, `S(k)`, 상거동.

> ★ **이 카드가 v1의 주 표적이다** (`D3`). M1–M2가 여기서 완주된다.
> 랩 논문에 직접 대응하는 것이 없어 **대부분 표준 문헌 기반**이다 — 그래서 `status: draft`.

---

## 1. 계

| 항목 | 값 |
|---|---|
| 입자 | 구형, 단분산 (`D10`) |
| 차원 | **3D** 기본 (`D9`) |
| 상호작용 | 배제부피 **WCA** (`D8` 닫힘 2026-07-28) + 선택적 Yukawa / Morse / DLVO |
| 경계 | 주기 |
| 용매 | 암묵적, **HI 없음** (`D11`) |

## 2. 목적 동역학

| | |
|---|---|
| 질문 | `g(r)` · `S(k)` · 압축인자 `Z(φ)` · 상거동 (액체/결정/유리) |
| 평형/비평형 | **평형** — 이 카드의 핵심 특징 |
| 정상상태 | 평형 상태. **평형화 판정이 필수적으로 유효하다** |

## 3. 기준 단위 — 표준

| 차원 | 선택 | 왜 |
|---|---|---|
| 길이 | **입자 지름 `σ`** | 계의 유일한 자연 길이 |
| 에너지 | **`kT`** — **입력값** | 평형계이므로 온도가 제어 변수 |
| 시간 | **`τ_D = σ²/D₀ = σ²γ/kT`** | 한 지름을 확산해 지나가는 브라운 시간 |

**유도와 검산**

```
D₀ = kT/γ                      ← Stokes–Einstein
γ  = 3πη·σ  (지름 기준)  또는  6πη·a  (반지름 기준)   ← ★ 규약 명시 필수
τ_D = σ²/D₀ = σ²γ/kT
τ_B = m/γ                      ← 관성 완화

검산: kT, γ, D₀ 중 둘이 주어지면 셋째는 유도된다. 셋 다 주어지면 상대오차 1e-3 안에서 일치할 것
      (선행 프로젝트 config.py:105 규칙 계승)
```

> **규약 함정 ★** Xu 2023은 `ζ = 3πη·d_ρ` (**지름** 기준)를 쓴다. 교과서 관례인 `6πηa`
> (반지름)와 **같은 값**이지만, 코드에서 `a`에 지름을 넣으면 2배 틀린다.
> `UnitMap`은 **어느 규약인지 필드로 기록**해야 한다.

## 4. 무차원수 원장

| 기호 | 정의 | 의미 | 문헌 범위 |
|---|---|---|---|
| **`φ`** | `Nπσ³/(6V)` (3D) | 패킹분율 | 액체 0–0.494 · **동결 0.494** · **융해 0.545** · RCP 0.64 |
| `ε/kT` | 인력 세기 | | 겔화는 대략 `> 2–3` |
| **`B₂*`** | `B₂/B₂^{HS}` | 환산 2차 비리얼 — **서로 다른 퍼텐셜을 한 축으로** | Noro–Frenkel: `≈ −1.5`에서 임계 |
| `λ/σ` | 인력 사거리 | 짧으면 `B₂*`로 환원 | |
| `κσ` | Debye 스크리닝 | 하전 콜로이드 | |
| `L/σ`, `L/(2r_cut)` | 박스 | **minimum image 위반 검사** | `r_cut ≤ L/2` |
| `τ_B/τ_D` | 오버댐프 타당성 | 콜로이드는 `~1e-6` | `< 1e-3` 요구 |
| `Δt/τ_D` | 기록용 | **게이트 아님** ([[dt-gate-should-be-displacement-based]]) | 2e-5 – 1e-3 (랩 실측) |

> **`B₂*`가 이 카드의 열쇠다.** 인력의 형태(Morse / square-well / depletion)가 달라도
> `B₂*`가 같으면 임계 거동이 거의 같다 (확장 대응상태 원리). **퍼텐셜 선택(`D8`)의
> 자유도를 하나로 줄여준다.**

## 5. 주요 파라미터 — 아직 랩 실측이 없다

| 파라미터 | 값 | 출처 |
|---|---|---|
| `dt` | **미정** — 변위 `√(2D₀Δt)/σ ≤ 0.03`에서 역산 | [[dt-gate-should-be-displacement-based]] |
| `N` | 미정. 3D `g(r)`는 보통 `1e3–1e4` | — |
| `L` | `r_cut ≤ L/2` 만족하도록 | — |
| 배제부피 | **WCA, `ε`=1 k_BT** (`φ_eff`≲0.32) / **`ε`=10** (조밀계). `D8` 닫힘 | [[wca-reproduces-carnahan-starling]] |
| nlist buffer | 0.4 (랩 관행 2건 일치) | Quah 코드 · 선행 slit |
| 초기배치 | **격자 → 부분추출 → 박스 압축** | Quah 코드에서 차용 |

> ⚠️ **이 표가 비어 있는 것이 현재 상태를 정확히 반영한다.** 랩은 능동물질·브러시·막을
> 연구하지 왔지 **수동 하드스피어 평형 구조를 하지 않았다.** 그래서 이 카드는
> 랩 지식이 아니라 **표준 문헌**으로 채워야 한다 (Frenkel–Smit, Allen–Tildesley, Hansen–McDonald).

## 6. 관측량

| 관측량 | 산출 | 오차 |
|---|---|---|
| `g(r)` | `freud.density.RDF` | block avg over frames |
| `S(k)` | `freud.diffraction.StaticStructureFactor` | block avg |
| `Z = P/(ρkT)` | 비리얼 | block avg + `τ_ac` |
| 도메인 크기 `R` | `2π/⟨k⟩` from `S(k)` | 조대화 분석용 |
| Steinhardt `q₆` | `freud.order` | 결정 판별 |

## 7. 적용 게이트 — 이 카드는 전부 켠다

| 게이트 | 적용 | 이유 |
|---|---|---|
| 스텝당 변위 `√(2D₀Δt)/σ ≤ 0.03` | ✅ | |
| 오버댐프 `τ_B/τ_D < 1e-3` | ✅ | |
| **평형화 (`pymbar`)** | ✅ **여기서는 유효** | 평형계이므로. [[abp--dense-collective]]와 정반대 |
| 자기일관성 `D_msd = kT/γ` | ✅ 희박 tracer | 능동 기여가 없어 그대로 성립 |
| 유한크기 `L` vs `1.5L` | ✅ | 임계점 근처에서 특히 |
| `r_cut ≤ L/2` | ✅ | minimum image |
| 이분산 강제 | ⚠️ `φ > 0.58` 근처 | 의도치 않은 결정화 방지. **문턱 미확정** |

## 8. 벤치마크 ★ — 이 카드가 벤치마크 밀도가 가장 높다

| ID | 검사 | 기대값 | 증거층 | 비용 |
|---|---|---|---|---|
| **`free_bd_stokes_einstein`** | 자유 BD `MSD = 6D₀t`, `D₀ = kT/γ` | 1% | ② | **극저** |
| **`carnahan_starling_hs_eos`** | `Z(φ) = (1+φ+φ²−φ³)/(1−φ)³` | 2% | ②③ | **저** |
| `hoover_ree_hs_freezing` | `φ_freeze = 0.494`, `φ_melt = 0.545` | 0.01 | ③ | 중 |
| `lj_eos_jzg` | LJ 상태방정식 표 상태점 | 2% | ③ | 중 |
| `noro_frenkel_b2star` | 해석 `B₂` vs 수치적분 | 1% | ② | 극저 |
| `coarsening_exponent_1_3` | 인력계 `R(t) ~ t^{1/3}` | ±0.05 | ② | 중 |

> ✅ **실현됐다** (2026-07-28) — [[wca-reproduces-carnahan-starling]]. 아래 예고대로 쟀고 WCA 가 `φ_eff`≲0.32 에서 2% 이내였다. 남은 것: harmonic core · Wang–Frenkel 은 재지 않았다.
>
> **`carnahan_starling_hs_eos`가 `D8`을 결정한다.** 배제부피 퍼텐셜 후보(harmonic / WCA /
> Wang–Frenkel)를 각각 넣고 `Z(φ)`를 재서, **CS와 가장 잘 맞는 것**을 고른다.
> 추측하지 않고 실측으로 정하는 방식.

## 9. 출처

랩 논문에 직접 대응하는 것이 없다. 표준 문헌으로 채운다 — 단계 F에서.

| | |
|---|---|
| 벤치마크 후보 | [`benchmarks/candidates.md`](../benchmarks/candidates.md) |
| `dt` 게이트 | [[dt-gate-should-be-displacement-based]] |
| 채워야 할 문헌 | Carnahan–Starling 1969 · Hoover–Ree 1968 · Noro–Frenkel 2000 · Frenkel–Smit |
