---
type: system
author: agent
drafted: 2026-07-27
confirmed_by:
system: interfacial-colloid
dynamics: transport
status: draft
cites:
  - knowledge/source/papers/2008-park-salt-surfactant-interface-forces.md
  - knowledge/source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md
  - knowledge/source/papers/2025-kim-manybody-hydrodynamics-optical-tweezers.md
---

# `유체계면 콜로이드` × `수송`

**묻는 질문:** 계면에 포획된 콜로이드는 어떻게 확산하는가 — `MSD`, `D`, `Ω` 불균질이
확산에 미치는 영향, 구조 완화 시간.

> ★ **이 카드에는 벤치마크가 없다. 그게 요점이다.**
> 확보한 계면 논문들은 **구조는 검증했지만 동역학은 시뮬레이션한 적이 없다.**
> [`2020-choi`](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md)는
> MC로 브라운 확산을 **의도적으로 뺐다**. 즉 이 카드가 다루는 영역은 **미검증 영역**이고,
> [[goal-autonomous-paper-to-sim-verification]]의 첫 구체적 표적이다.
> 짝 카드: [`interfacial-colloid--equilibrium-structure`](interfacial-colloid--equilibrium-structure.md)

---

## 1. 계

구조 카드 §1과 동일하다 — 2D, 쌍극자 반발 `U/kT = Ω_iΩ_j(d/r)³`, 암묵적 용매.
**차이는 시간축이 살아 있다는 것뿐이고, 그래서 항력이 필요하다.**

## 2. 목적 동역학

| | |
|---|---|
| 질문 | `MSD(t)` · `D` · `Ω` 불균질 → 확산 이질성 · 구조 완화 시간 |
| 평형/비평형 | **평형** (평형 상태에서의 수송) |
| 정상상태 | 평형. 단 **MSD는 과도 구간을 봐야 한다** |
| 엔진 | **BD 필수** — MC로는 답할 수 없다 |

## 3. 기준 단위 ★ — 여기가 구조 카드와 갈린다

| 차원 | 선택 | 왜 |
|---|---|---|
| 길이 | **입자 지름 `d`** | 구조 카드와 동일 (`Ω` 정의가 `d` 기준) |
| 에너지 | **`kT`** — 입력값 | 평형계 |
| 시간 | **`τ_D = d²/D₀`**, `D₀ = kT/γ` | 한 지름을 확산해 지나가는 시간 |

### ★ 항력계수 — `η`가 아니라 `η_eff`다

$$\boxed{\;\gamma = 3\pi\,\eta_{eff}\,d = 6\pi\,\eta_{eff}\,a\;}$$

$$\eta_{eff} = \frac{\eta_{oil}(1-\cos\theta_c) + \eta_{water}(1+\cos\theta_c)}{2}$$

**두 벌크 점도의 표면적 가중 평균이다.** `(1−cos θ_c)/2`와 `(1+cos θ_c)/2`가 각각
오일·물에 잠긴 구의 표면적 분율이다. `θ_c`는 **물 상 기준**.

출처: [`2008-park-salt-surfactant-interface-forces`](../../source/papers/2008-park-salt-surfactant-interface-forces.md)
(Park·Pantina·**Furst**·Oettel·Reynaert·Vermant, Langmuir 2008, 24, 1686) — `F_S = 6πaη_eff U`.
[`2020-choi`](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md) SI의
`κ_t = 3πdη_eff u/Δx`와 같은 식이다.

**극한 검산**

| `θ_c` | `η_eff` |
|---|---|
| 0° | `η_water` (완전히 물 속) |
| 90° | `(η_oil + η_water)/2` (단순 평균) |
| 180° | `η_oil` (완전히 오일 속) |

**유도와 검산**

```
η_eff = [η_oil(1−cos θ_c) + η_water(1+cos θ_c)]/2
γ     = 3π η_eff d                    ← 지름 기준. 반지름이면 6π η_eff a
D₀    = kT/γ
τ_D   = d²/D₀ = d²γ/kT
τ_B   = m/γ

검산 1: θ_c = 90° 를 넣으면 η_eff = (η_oil+η_water)/2 로 떨어질 것
검산 2: η_oil = η_water 를 넣으면 θ_c 와 무관하게 η_eff = η 가 될 것
검산 3: kT, γ, D₀ 중 둘이 주어지면 셋째가 유도된다 (상대오차 1e-3)
```

> **규약 함정 ★ 두 개가 동시에 걸린다.**
> ① `γ = 3πηd`(지름) vs `6πηa`(반지름) — 같은 값이지만 코드에서 뒤바꾸면 **2배**
> ② `θ_c`를 오일 상 기준으로 넣으면 `η_eff`가 두 상 사이에서 **반대로** 간다
> `UnitMap`은 **길이 규약과 접촉각 기준을 둘 다 필드로 기록**해야 한다.

## 4. 무차원수 원장 ★

구조 카드 §4의 항목(`Ω`, `a_ij`, `r/d`, `φ_2D`, `Bo`, `r_cut/d`)을 전부 승계하고, 아래를 추가한다.

| 기호 | 정의 | 의미 | 범위 |
|---|---|---|---|
| **`η_eff/η_water`** | 위 식 | 계면 항력 보정 | decane/water: **≈ 0.95** (거의 1) |
| `τ_B/τ_D` | `m kT/(γ²d²)` | 오버댐프 타당성 | 콜로이드 `~1e-6`. **`< 1e-3` 요구** |
| `Δt` 변위 | `√(4D₀Δt)/d` | **2D는 4** (3D는 6) | `≤ 0.03` ([[dt-gate-should-be-displacement-based]]) |
| **`Bq`** | `η_s/((η_oil+η_water)a)` | **Boussinesq — 표면점도 대 벌크** | ⚠️ **`η_eff`가 무시하는 것.** §7 |
| **`ℓ_SD/d`** | `η_s/(2η d)` | Saffman–Delbrück 스크리닝 길이 | ⚠️ HI 관련. §7 |

### decane/water에서는 `θ_c` 의존성이 사실상 사라진다

`η_water ≈ 0.89`, `η_decane ≈ 0.84 mPa·s` (25 °C, **문헌 상온값 — 원 논문에 없음, 출처 확보 필요**):

| `θ_c` | `η_eff` (mPa·s) |
|---|---|
| 99.4° | ≈ 0.861 |
| 142.6° | ≈ 0.845 |
| 145.4° | ≈ 0.844 |

99–145° 전 구간 변동 **~2%**, `142.6°` vs `145.4°`는 **0.1% 미만**.
따라서 구조 카드 §5의 APS/CPS 접촉각 불일치는 이 카드에 **실질 영향이 없다.**

> ⚠️ **이 논거는 decane/water 한정이다.** 물/실리콘오일처럼 점도차가 10배 이상인 계면에서는
> `θ_c`가 `η_eff`를 크게 흔든다. 유체쌍이 바뀌면 이 절을 다시 계산할 것.

## 5. 주요 파라미터

구조 카드 §5의 실측표(`Ω`·`a_ij`·`θ_c`·`ψ`·`d`)를 그대로 쓴다. 추가로 필요한 것:

| 파라미터 | 값 | 출처 |
|---|---|---|
| `η_eff` | `≈ 0.85 mPa·s` (decane/water, `θ_c ≈ 100–145°`) | [Park 2008] 식 + 문헌 점도 **[미재현]** |
| `η_oil`·`η_water` | **미확보** — 표준 물성표 필요 | ❌ 원 논문에 없음 |
| 온도 `T` | **미확보** (ambient 추정) | ❌ 두 논문 모두 명시 없음 |
| `dt` | **미정** — 변위 `√(4D₀Δt)/d ≤ 0.03`에서 역산 | [[dt-gate-should-be-displacement-based]] |
| `N` | 247 (원 논문 배치 재사용) 또는 자체 설정 | [2020-choi] |
| 초기배치 | 실험 스냅샷 또는 평형화된 MC 배치 | |

> ⚠️ **`T`와 `η_oil`·`η_water`가 없으면 `τ_D`를 SI로 못 박는다.** 무차원 결과는 낼 수 있으나
> 실험과 시간축을 맞추려면 이 둘이 필요하다. **S2 ELICIT에서 물어야 할 항목.**

## 6. 관측량

| 관측량 | 산출 | 오차 |
|---|---|---|
| **`MSD(t)`** | `⟨\|r(t)−r(0)\|²⟩`, **2D → `MSD = 4Dt`** | block avg + `τ_ac` |
| `D` | MSD 장시간 기울기 / 4 | block-averaged 오차막대 필수 |
| `D_i` 개별 확산계수 | 입자별 MSD | **`Ω_i`와의 상관이 이 카드 고유 질문** |
| VACF | 속도 자기상관 | 관성 무시 확인용 |
| 구조 완화 시간 | `g(r)` 또는 `ψ₆`의 시간 상관 | |

> **2D 계수 함정:** `MSD = 2·dim·D·t` 이므로 **2D는 `4Dt`**다. 3D 습관대로 `6Dt`를 쓰면
> `D`가 1.5배 틀린다.

## 7. 적용 게이트 ★

| 게이트 | 적용 | 이유 |
|---|---|---|
| **자기일관성 `D_msd = kT/γ`** | ✅ **필수 — 이 카드의 핵심 검사** | 희박 극한에서 성립. **`γ`에 `η_eff`를 넣었는지 검증하는 유일한 장치** |
| 스텝당 변위 `√(4D₀Δt)/d ≤ 0.03` | ✅ | **2D 계수 4** 확인할 것 |
| 오버댐프 `τ_B/τ_D < 1e-3` | ✅ | |
| **평형화 (`pymbar`)** | ✅ 유효 | 평형계 |
| 유한크기 `L` vs `1.5L` | ✅ | 구조 카드와 달리 배치를 자유롭게 정할 수 있으므로 **여기서는 가능** |
| `r_cut` 절단 | ⚠️ 구조 카드 §7과 동일 — `r_cut/d ≥ 21–69` | |
| **HI 무시** | ⚠️ **3D보다 위험. 아래 참조** | |
| **표면점도 무시** | ⚠️ **`η_eff`의 내재 한계. 아래 참조** | |

### ⚠️ HI — 준2D는 3D보다 까다롭다

`master_plan.md` §12가 HI를 명시적 비목표로 두었고 그 결정은 유효하다. 다만 **준2D에서
free-draining 근사의 실패 방식이 3D와 다르다**는 것을 리포트에 적어야 한다.

순수 2D Stokes 유동에는 유계해가 없고(Stokes 역설), 운동량이 벌크 상으로 새면서
**Saffman–Delbrück 스크리닝 길이 `ℓ_SD = η_s/2η`** 가 생긴다. 3D처럼 `1/r` 원거리 감쇠
뒤에 숨을 수가 없다.

> 랩 자산 [`2025-kim-manybody-hydrodynamics-optical-tweezers`](../../source/papers/2025-kim-manybody-hydrodynamics-optical-tweezers.md)
> 가 다체 HI를 **직접 측정**한 논문이다 — 이 축이 살아 있는 문제라는 증거.
>
> **S0 후보 선별 규칙:** 검증하려는 주장이 HI에 의존하면 no-HI BD는 검증자가 될 수 없다.
> 후보에서 제외한다. 이 필터가 없으면 "논문이 틀림"과 "우리 엔진에 그 물리가 없음"이
> 구분되지 않는 불일치가 쏟아진다.

### ⚠️ 표면점도 — `η_eff`가 안고 있는 근사

[Park 2008] 논문 전체에 surface viscosity · interfacial rheology · Boussinesq 언급이
**0건**이다. 이 `η_eff`는 벌크 두 상의 표면적 가중 평균일 뿐이다.

- 계면활성제가 없는 순수 decane/water: **문제 없음** — 이 카드의 기본 조건
- SDS 등이 있으면: `Bq > 1` 에서 **항력을 과소평가**한다. [Park 2008]은 CMC 훨씬 아래에서만
  검증했다. **계면활성제 농도가 올라가면 이 카드의 `γ`는 무효**

## 8. 벤치마크 — **없다**

| ID | 검사 | 상태 |
|---|---|---|
| `free_bd_2d_stokes_einstein` | 상호작용 없는 2D 자유 BD `MSD = 4D₀t`, `D₀ = kT/γ` | ✅ **쓸 수 있음** — 증거층 ②, 비용 극저 |
| `eta_eff_limits` | `θ_c=90°→` 단순평균, `η_oil=η_water→η` | ✅ **쓸 수 있음** — 해석 항등식, 비용 극저 |
| 계면 콜로이드 `D` 실측 대조 | — | ❌ **없음.** 확보한 논문 어디에도 계면 `MSD`/`D` 실측이 없다 |
| `Ω` 불균질 → 확산 이질성 | — | ❌ **없음.** 아무도 안 했다 |

> **위 두 개는 이 계의 검증이 아니라 엔진 위생 검사다.** 계면 특유의 동역학을 검증할
> 오라클은 현재 **하나도 없다.** 따라서 이 카드로 낸 결과에는 반드시 `unverified` 배지가
> 붙는다 (`master_plan.md` §12).

### 오라클을 만드는 두 가지 길

| 길 | 내용 | 비용 |
|---|---|---|
| **A. 실험 데이터 확보** | Park 랩의 계면 입자 트래킹 궤적. `D24`(실험을 5번째 증거층으로)가 여기서 발동된다 | 사람 협의 필요 |
| **B. 구조 카드로 교차검증** | 같은 퍼텐셜·같은 배치로 BD를 돌려 **평형 `g(r)`이 MC 결과를 재현**하는지 먼저 확인. 구조가 맞으면 동역학을 신뢰할 근거가 하나 생긴다 | **낮음 — 먼저 할 것** |

> **B를 먼저 한다.** [구조 카드](interfacial-colloid--equilibrium-structure.md) §8의
> `choi2020_rdf_hetero`를 BD로 재현하면, 이 카드의 `γ`·`Δt`·절단이 최소한 구조를 망치지는
> 않는다는 것이 확인된다. 그 위에서 MSD를 주장한다.

## 9. 출처

| | |
|---|---|
| **항력 `η_eff`** | [`2008-park-salt-surfactant-interface-forces`](../../source/papers/2008-park-salt-surfactant-interface-forces.md) — Park·Furst |
| 퍼텐셜·파라미터 | [`2020-choi-electrostatic-self-potential-heterogeneity`](../../source/papers/2020-choi-electrostatic-self-potential-heterogeneity.md) |
| HI 경고 근거 | [`2025-kim-manybody-hydrodynamics-optical-tweezers`](../../source/papers/2025-kim-manybody-hydrodynamics-optical-tweezers.md) |
| 짝 카드 | [`interfacial-colloid--equilibrium-structure`](interfacial-colloid--equilibrium-structure.md) |
| `Δt` 게이트 | [[dt-gate-should-be-displacement-based]] |
