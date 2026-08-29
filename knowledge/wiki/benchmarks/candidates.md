---
type: benchmark
subtype: candidate-list
author: agent
drafted: 2026-07-27
confirmed_by:
status: draft
cites:
  - knowledge/source/papers/2022-barakat-enhanced-dispersion-harmonic-traps.md
  - knowledge/source/papers/2022-modica-porous-media-active-diffusion.md
  - knowledge/source/papers/2023-arnold-active-surface-flows-coarsening.md
  - knowledge/source/papers/2025-gubbala-phase-field-viscous-inclusions.md
---

# 벤치마크 후보 — 랩 논문에서 추출

랩 논문 16편을 읽고 **`benchmarks.yaml`에 등재할 수 있는 정량 검사**를 추렸다.

> **선별 기준**
> ① 값이 명시적일 것 ② 조건이 좁게 특정될 것 ③ 10코어 M4에서 감당 가능할 것
> ④ **우리 v1 모델(BD, 구형, HI 없음)으로 재현 가능할 것**
>
> ④를 통과 못 하는 것도 기록한다 — **"우리가 재현할 수 없는 물리"의 목록**도 자산이다.

---

## A. 등재 권장 — v1으로 재현 가능

### A1 · `coarsening_exponent_diffusive` ★ 최우선

| | |
|---|---|
| **검사** | 확산 지배 상분리에서 도메인 크기 `R(t) ~ t^{1/3}` |
| **기대값** | 지수 `1/3` |
| **허용오차** | log–log 기울기 `0.33 ± 0.05` |
| **증거층** | **②** 해석적 극한 (Lifshitz–Slyozov / 확산 제한 응집) |
| **비용** | 중 — 상분리가 일어날 만큼은 돌려야 함 |
| **출처** | Arnold, Gubbala, Takatori, *PRL* **131**, 128402 (2023) |

**왜 좋은가:** 계에 무관한 **보편 지수**다. 지질막에서 나온 값이지만 콜로이드 겔의 스피노달
분해에도 그대로 적용된다. 그리고 **우리가 HI를 무시한다는 사실 자체를 검사한다** —
유체역학이 있으면 지수가 `1`(점성 지배)로 바뀌므로, `1/3`이 나오는 것이
`D11`(free-draining)의 자기일관성 확인이 된다.

**구현:** `S(k)` 평균 파수 → `R = 2π/⟨k⟩` → `log R` vs `log t` 기울기.
`07_observables.md`에 `S(k)`가 이미 있으므로 **변환만 추가**하면 된다.

### A2 · `barakat_trap_D_limit`

| | |
|---|---|
| **검사** | 조화 트랩 배열에서 `κ → 0` 극한에 `D → kT/γ` |
| **기대값** | `D(κ→0) = kT/γ` |
| **허용오차** | 1% |
| **증거층** | **②** 해석적 극한 |
| **비용** | **저** — 비상호작용 BD |
| **출처** | Barakat & Takatori, *Phys. Rev. E* **107**, 014601 (2023) |

**왜 좋은가:** `free_bd_stokes_einstein`의 자연스러운 확장. 트랩을 켜고 강성을 0으로 보내면
자유 확산으로 수렴해야 한다 — **힘 구현이 맞는지를 극한으로 검사**한다.
선행 프로젝트에 `HarmonicTrapForce`가 이미 있어 구현 비용이 낮다.

### A3 · `barakat_D_monotonic_in_kappa`

| | |
|---|---|
| **검사** | 정지 트랩에서 `D(κ)`가 `κ`에 대해 **단조 감소** |
| **기대값** | 단조성 (정성적) |
| **증거층** | ③ 문헌 |
| **비용** | 저 |

정성적 검사지만 싸고, 부호 오류·단위 오류를 잘 잡는다.

---

## B. 조건부 — 확보·확장 필요

### B1 · `modica_active_diffusion_reduction`

| | |
|---|---|
| **검사** | 고정 장애물 `φ = 12%`에서 ABP 유효 `D` **25% 감소** |
| **부가** | 수동 입자 대비 감소폭 약 **2배** |
| **증거층** | ③ 문헌 + ④ 실험 |
| **출처** | Modica, Xi, Takatori, *Front. Phys.* **10**, 869175 (2022) |
| **막는 것** | **2D다** (우리 `D9` 기본값은 3D) · 장애물 좌표가 실험 유래 |

값이 명시적(25%)이고 조건이 좁아 **매력적이다.** 2D 경로가 생기면 1순위로 등재.
`U₀ = 0.84 µm/s`, `τ_R = 14 s`, 장애물 5 µm → `Pe`로 환산해 등재할 것.

### B2 · `coarsening_elastic_suppression`

| | |
|---|---|
| **검사** | 탄성 매질 속 상분리에서 `R(t) ~ t^α`, **`α < 1/4`** |
| **증거층** | ② 이론 (Cahn–Hilliard + Landau–de Gennes) |
| **출처** | Gubbala, Jena, Arnold, Takatori, *Soft Matter* **21**, 6243 (2025) |
| **막는 것** | 탄성 네트워크 모델링 필요 — v1 범위 밖 |

### B3 · `xu_lamellar_spinodal`

| | |
|---|---|
| **검사** | 주어진 `(f, m, k_BTR_c/a)`에서 스피노달 `φ` |
| **출처** | Xu, Jandhyala, Takatori, *J. Chem. Phys.* **161**, 064901 (2024) |
| **막는 것** | **`Δt` 미확보** · 8×10⁸ 스텝은 M4에서 비현실적 |

---

## C. 재현 불가 — "우리가 다루지 않는 물리"의 기록

이 항목들은 `benchmarks.yaml`에 **등재하지 않는다.** 대신 결과에 `unverified` 배지를 붙일 때
**"왜 검증할 수 없는가"의 근거**로 인용한다.

| 항목 | 왜 불가 | 출처 |
|---|---|---|
| 다체 유체역학 상호작용 | **HOOMD BD/Langevin에 HI가 없다.** RPY 이동도 텐서 필요 | Kim et al., *PRFluids* **10**, 064301 (2025) |
| 시간의존 HI (점탄성) | 마찰계수가 상수인 순간 유체 기억을 표현 못 함 | Kim et al., *Soft Matter* (2025) |
| 비등방 확산 텐서 `D∥/D⊥` | `md.methods.Brownian`의 `gamma`는 타입별 스칼라 | Choi et al., *JCIS* **721**, 140705 (2026) |
| 운동성 박테리아의 상 분배 | 두 상 화학친화도 + 계면 모델링 필요, 비구형 | Cheon et al., *PRL* (2025) |

> **이 표가 `master_plan.md` §12(명시적 비목표)에 정량적 근거를 준다.**
> "HI를 하지 않는다"를 그냥 선언하는 대신, **어떤 실측 데이터를 재현할 수 없는지**를
> 구체적으로 댈 수 있다.

---

## 다음 할 일

1. **A1·A2를 `benchmarks.yaml`에 실제로 등재** — 문헌 코퍼스 단계(F)에서
2. Modica 3편 paywall 해소 → B1 승격 가능성 재평가
3. `05_validation_gates.md`에 조대화 지수를 **S8 진단 항목**으로 추가 검토
   (상분리 계에서만 적용되는 조건부 게이트)
