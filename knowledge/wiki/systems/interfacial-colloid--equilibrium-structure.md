---
type: system
author: agent
drafted: 2026-07-27
confirmed_by:
system: interfacial-colloid
dynamics: equilibrium-structure
status: draft
cites:
  - knowledge/source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md
  - knowledge/source/papers/2008-park-salt-surfactant-interface-forces.md
  - knowledge/source/papers/2023-lee-colloidal-debye-force.md
---

# `유체계면 콜로이드` × `평형 구조`

**묻는 질문:** 오일-물 계면에 포획된 하전 콜로이드 단층막의 평형 미세구조는 무엇인가
— `g(r)`, Voronoi 결함, 육방 질서.

> ★ **랩 대응 논문이 있는 첫 계면 카드다.** [`2020-choi`](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md)
> 가 이 계를 MC로 실제로 돌렸고 실험 `g(r)`과 대조까지 했다. 파라미터가 실측으로 채워진다.
> 짝 카드: [`interfacial-colloid--transport`](interfacial-colloid--transport.md)

---

## 1. 계

| 항목 | 값 |
|---|---|
| 입자 | 구형 폴리스티렌, 설페이트/아민/카복실 표면. 단분산 (`d ≈ 3 μm`) |
| **차원** | **2D** — 아래 참조 |
| 상호작용 | **쌍극자 정전 반발** `U_ij/kT = Ω_i Ω_j (d/r_ij)³` · `F ∼ r⁻⁴` |
| 다분산 | 크기는 단분산이나 **`Ω`가 입자마다 다르다** (자기퍼텐셜 불균질) ← 이 계의 특징 |
| 경계 | 실험은 열린 단층막. 시뮬은 주기 또는 **실험 스냅샷 고정배치** |
| 용매 | 암묵적, **HI 없음** (`D11`) — ⚠️ 준2D에서 더 위험, §7 참조 |

### 왜 2D로 다뤄도 되는가

탈착 에너지가 `ΔE = πR²γ_int(1−|cos θ_c|)²`이고, `R ~ 1.5 μm`·`γ_int ~ 50 mN/m`이면
**10⁶ kT 규모**다. z 자유도는 열적으로 완전히 얼어 있다 — 슬릿 포어보다 깨끗한 2D다.
따라서 `dimension: 2`는 근사가 아니라 사실상 정확하다.

## 2. 목적 동역학

| | |
|---|---|
| 질문 | `g(r)` · Voronoi 결함(5/6/7 이웃) · 육방 질서 `ψ₆` · `Ω` 불균질이 구조에 미치는 영향 |
| 평형/비평형 | **평형** |
| 정상상태 | 평형 상태. **평형화 판정이 유효하다** |
| 엔진 | **MC로 충분** — 원 논문이 MC를 썼고 브라운 확산을 무시했다 (§5) |

## 3. 기준 단위 ★

| 차원 | 선택 | 왜 |
|---|---|---|
| 길이 | **입자 지름 `d`** | 문헌 규약이 `d`다 — `(d/r)³`, `r/d ≈ 8.3`. **반지름 아님** |
| 에너지 | **`kT`** — 입력값 | 평형계. `Ω`가 이미 `kT` 단위로 무차원화되어 있다 |
| 시간 | **없음 (MC)** | MC 사이클은 물리 시간이 아니다. BD로 돌릴 경우 → [수송 카드](interfacial-colloid--transport.md) §3 |

**유도와 검산**

```
U_ij/kT = Ω_i Ω_j (d/r_ij)³        ← Ω 는 이미 무차원. 추가 환산 없음
a_ij    = Ω_i Ω_j                   ← 쌍상호작용 크기
F_ij/kT = 3 Ω_i Ω_j d³ / r_ij⁴      ← F = −dU/dr

검산: a_ij 실측값과 Ω_i·Ω_j 곱이 일치할 것
      CPS: Ω = 571.7 → Ω² = 326,841  vs  실측 a_ij = 326,985  → 0.04% 일치 ✅
```

> **규약 함정 ★ — `d`는 지름이다.** `U ∝ (d/r)³` 에 반지름을 넣으면 **8배** 틀린다.
> [수송 카드](interfacial-colloid--transport.md)의 항력 규약 함정(`3πd` vs `6πa`)과 짝이다.

## 4. 무차원수 원장 ★

| 기호 | 정의 | 의미 | 문헌 범위 |
|---|---|---|---|
| **`Ω_i`** | 입자 자기퍼텐셜 | 이 계의 주 파라미터 | **84.7 – 571.7** (실측, §5) |
| **`a_ij`** | `Ω_i Ω_j` | 쌍상호작용 크기 | `9.4e3 – 3.3e5` (실측) |
| `r/d` | 평균 입자간 거리 | 희박도 | **≈ 8.3** (실측) |
| **`U_ij/kT` @ `r̄`** | `a_ij (d/r̄)³` | 반발이 열에너지를 얼마나 압도하는가 | **≈ 577** (실측) |
| `φ_2D` | `Nπd²/(4A)` | 면적분율 | `≈ 0.0091` (`r/d`=8.3에서 유도) |
| `θ_c` | 삼상 접촉각 (물 상 기준) | 각 상 노출 면적 결정 | **97 – 145°** (실측) |
| `Bo` | `Δρ g R²/γ_int` | 중력 vs 모세관 (flotation) | `~1e-6` — **무시 가능** (`R ~ 1.5 μm`) |
| **`r_cut/d`** | 절단 거리 | §7 참조 — **이 계의 최대 함정** | **≥ 21 – 69** 필요 |

> **`r⁻³`는 2D에서 열역학적 장거리는 아니다** (`∫ r⁻³·r dr` 수렴, 조건은 `r⁻²`보다 빠른 감쇠).
> 그러나 **계수 `a_ij`가 `1e5` 규모**라 실질 사거리가 수십 지름에 이른다. §7의 절단 게이트 참조.

## 5. 주요 파라미터 — 실측값 ★

출처: [`2020-choi`](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md) Table 1 · §4.
**전부 `reproduced: no`** — 사실 기록이지 검증된 근거가 아니다 (`[출처, 미재현]`으로 인용할 것).

| 입자 | `Ω` | `a_ij` | `θ_c` (deg) | `ψ` (mV) | `d` (μm) |
|---|---|---|---|---|---|
| NPS | 인력 (반발 없음) | — | 97.0 ± 0.9 | −51.2 ± 2.4 | 2.93 ± 0.03 |
| SPS | 84.71 ± 19.45 | 9422 ± 6978 | 99.4 ± 2.1 | −57.5 ± 2.2 | 2.96 ± 0.05 |
| APS | 136.18 ± 40.48 | 18100 ± 9108 | 145.4 ± 0.9 | −69.5 ± 2.4 | 2.79 ± 0.11 |
| CPS | 571.7 ± 94.16 | 326985 ± 79935 | 142.6 ± 1.2 | −65.6 ± 3.2 | 3.16 ± 0.07 |

> ⚠️ `θ_c(APS)`·`θ_c(CPS)`는 원문 Table 1과 SI Figure S8이 **서로 뒤바뀌어 있다.**
> 정본은 Table 1로 고정했다 (제1저자 판단, 2026-07-27). 상세는 증류본 §3.

**MC 설정 (원 논문이 실제로 쓴 값)**

| 파라미터 | 값 | 비고 |
|---|---|---|
| `N` | **247** | 실험 스냅샷의 입자 수 |
| 초기배치 | **실험 현미경 이미지에서 취득** | Figure S4A,B — 무작위 아님 |
| MC 이동 스텝 | `Δ(r/d) = 0.05` (radial) | |
| 사이클 | `10⁶` | |
| 평형 판정 | 평균 입자간 거리 **plateau** | `pymbar` 아님 |
| 샘플링 | 4000+ 사이클, 사이클마다 `g(r)` 후 평균 | |
| 브라운 확산 | **무시** | `U_ij ≈ 577 kT` 라 정당화 |
| 비교 조건 | `Ω_hetero` · `Ω_max = 804.33` · `Ω_min = 202.18` · `a_hetero` | 네 가지 |

## 6. 관측량

| 관측량 | 산출 | 오차 |
|---|---|---|
| **`g(r)` (2D)** | `freud.density.RDF` (`box`를 2D로) | block avg over cycles |
| **Voronoi 결함** | `freud.locality.Voronoi` → 최근접이웃 수 5/6/7 분포 | block avg |
| **육방 질서 `ψ₆`** | `freud.order.Hexatic(k=6)` | 2D 전용. 원 논문엔 없으나 자연스러운 확장 |
| `Ω` 분포 vs 국소구조 상관 | 자체 구현 | 이 계 고유 — 불균질의 구조적 귀결 |

> **Voronoi 결함이 `g(r)`보다 나은 검사일 수 있다.** 원 논문 Figure S4C–F가 실험과 MC 네 조건의
> 결함 패턴을 나란히 놓는다 — `g(r)`보다 조건 구분력이 크다.

## 7. 적용 게이트 ★

| 게이트 | 적용 | 이유 |
|---|---|---|
| **평형화 (`pymbar`)** | ✅ **유효** | 평형계다. 단 원 논문은 plateau 육안 판정을 썼으므로 우리가 더 엄격해진다 |
| 스텝당 변위 `√(2D₀Δt)/d` | ❌ **해당 없음 (MC)** | MC에는 `Δt`가 없다. MC 스텝 `Δ(r/d)=0.05`가 대응물 — 수용률로 관리 |
| 오버댐프 `τ_B/τ_D` | ❌ **해당 없음 (MC)** | 시간축이 없다 |
| 자기일관성 `D_msd = kT/γ` | ❌ **해당 없음 (MC)** | 〃. BD로 돌리면 → [수송 카드](interfacial-colloid--transport.md) |
| **`r_cut` 절단** | ⚠️ **이 카드의 최대 함정 — 아래 참조** | |
| 유한크기 `L` vs `1.5L` | ⚠️ **불가** | 초기배치가 실험 스냅샷이라 `N=247` 고정. 재현 시 이 한계를 명시할 것 |
| MC 수용률 | ✅ 20–50% 목표 | `Δ(r/d)=0.05`가 이 범위를 주는지 확인 |

### ★ 절단 거리 — 박스가 상호작용 사거리를 못 담는다

`U_ij/kT = a_ij (d/r)³ = 1` 이 되는 거리:

```
r/d = a_ij^(1/3)
  SPS  a=9,422    → r/d ≈ 21.1
  APS  a=18,100   → r/d ≈ 26.3
  CPS  a=326,985  → r/d ≈ 68.9      ← 최악
```

그런데 `N=247`, `r̄/d = 8.3`이면 2D 박스는 `L/d ≈ √(247 × 8.3²) ≈ 130`.
**minimum image 요구 `r_cut ≤ L/2 = 65 d` < CPS가 필요로 하는 68.9 d.**

> **CPS 조건에서는 박스가 `U = kT` 사거리조차 담지 못한다.** 원 논문은 절단을 언급하지 않는다
> (전 쌍 합산으로 보이나 미확인). 재현할 때 **반드시 절단 처리를 명시**하고, 절단 반경 민감도를
> 돌려야 한다. 이 항목이 재현 실패의 1순위 후보다.

## 8. 벤치마크

| ID | 검사 | 기대값 | 증거층 | 비용 |
|---|---|---|---|---|
| **[`choi2020_interfacial_rdf`](../benchmarks/choi2020-interfacial-rdf.md)** | `g(r)` 봉우리 높이·위치 | **1차 피크 `3.88 ±15%` @ `r/d 8.0 ±0.3`** 외 | ③ 코드-코드 + ④ 실험 | 중 |
| ↳ **판별 검사** ★ | `Ω_hetero` 통과 · `Ω_max`(+135%)·`Ω_min`(+83%) 실패 | 세 조건의 **순서** | ④ | 중 |
| `choi2020_voronoi_defects` | 5/6/7 이웃 분포가 실험과 일치 | ⚠️ 정성 — 정량화 필요 | ④ | 중 |
| `omega_squared_consistency` | `a_ij = Ω_i Ω_j` 항등식 | 0.1% | ② 내부 | **극저** |

> ✅ **`g(r)` 기대값은 곡선 판독으로 정량화했다** (2026-07-27).
> 원 논문의 "excellent agreement"는 게이트가 아니었으므로(원칙 7) Figure 4B를 색분리로
> 디지타이즈했다. 축 보정 검증은 `exp` 꼬리가 `1.008 ± 0.114`로 `g(r)→1`에 수렴한 것.
> 판독 스크립트 [`docs/tools/digitize_fig4b.py`](../../../docs/tools/digitize_fig4b.py)로 재현 가능.
>
> **판별력이 크다.** 1차 피크 하나로 세 조건이 갈린다 — `Ω_hetero` +0%, `Ω_min` +83%,
> `Ω_max` +135%. 불균질을 안 넣으면 피크가 2배 이상 솟는다.
>
> ⚠️ **1차 최소는 게이트로 쓰지 않는다** — `Ω_hetero`조차 실험 대비 −57%다.
> 원 논문의 일치 주장은 1차 피크와 전체 형상에 대한 것이지 모든 지점이 아니다.
>
> `omega_squared_consistency`는 지금 바로 쓸 수 있다 — CPS에서 0.044% 일치 확인됨 (§3).

## 9. 출처

| | |
|---|---|
| 주 출처 | [`2020-choi-electrostatic-self-potential-heterogeneity`](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md) — 파라미터·MC 설정·벤치마크 |
| 항력 (수송 카드용) | [`2008-park-salt-surfactant-interface-forces`](../../source/papers/2008-park-salt-surfactant-interface-forces.md) |
| 같은 계 추가 | [`2023-lee-colloidal-debye-force`](../../source/papers/2023-lee-colloidal-debye-force.md) · [`2019-kang-mapping-anisotropic-colloidal-interactions`](../../source/papers/2019-kang-mapping-anisotropic-colloidal-interactions.md) |
| 짝 카드 | [`interfacial-colloid--transport`](interfacial-colloid--transport.md) |
