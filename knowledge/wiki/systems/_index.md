---
type: system
subtype: index
author: agent
drafted: 2026-07-27
confirmed_by:
---

# 계-동역학 카드 — 색인

> **왜 계 하나로는 부족한가**
>
> 무차원화 방법과 주요 파라미터는 **계만으로 정해지지 않는다.** 같은 콜로이드라도
> 평형 구조를 보느냐 수송을 보느냐에 따라 **기준 시간이 달라지고**, 따라서 `Δt` 게이트도
> 관측량도 벤치마크도 달라진다.
>
> 실제 증거: 같은 랩 안에서 **기준 단위가 셋으로 갈린다.**
>
> | 계 · 목적 동역학 | 기준 길이 | 기준 시간 | `kT` |
> |---|---|---|---|
> | ABP × 제어 (Quah) | **런 길이 `ℓ`** | **`τ_r = 1/D_r`** | **유도량** |
> | 브러시 콜로이드 × 비평형 접촉 (Xu) | `σ` | `τ_D = σ²/D` | 입력값 |
> | 수동 tracer × 수송 (선행 slit) | `σ` | `τ_D` | 입력값 |
>
> 하나의 무차원화 규약을 강요하면 **셋 중 둘이 부자연스러워진다.**

## 에이전트는 이 카드를 어떻게 쓰는가

```
S1 INTAKE    → (계, 목적 동역학) 쌍을 분류한다.  분류 실패 시 unknowns[] 로
S3 NONDIM    → 해당 카드의 §3 기준 단위 · §4 무차원수 원장을 적용한다
S5 PLAN      → §5 주요 파라미터를 사전(prior)으로 쓴다
S6/S8 게이트  → §7의 켜고 끄기를 따른다  ← ★ 게이트는 쌍마다 다르다
S10 ANALYZE  → §6 관측량 목록
검증          → §8 벤치마크
```

> 카드가 없는 쌍을 만나면 **새 카드를 `_TEMPLATE.md`로 만들고 `status: draft`로 시작한다.**
> 카드 없이 진행하면 무차원화를 즉흥으로 하게 되고, 그게 `master_plan.md` §2-a가 말한
> "파라미터 선택이 암묵지다"의 재발이다.

---

## 매트릭스

행 = 계 · 열 = 목적 동역학. **✅ 카드 있음 · ○ 필요하나 없음 · — 해당 없음**

| 계 \ 목적 동역학 | 평형 구조 | 수송 (`D`, MSD) | 상거동·조대화 | 비평형 과도 | 조밀 집단 |
|---|---|---|---|---|---|
| **수동 구형** | **✅ [카드](passive-sphere--equilibrium-structure.md)** | ○ | ○ | **✅ [트랩 카드](passive-sphere--harmonic-trap.md)** | ○ |
| **인력 콜로이드** (depletion/Morse) | ○ | — | ○ | — | ○ |
| **하전 콜로이드** (Yukawa/DLVO) | ○ | — | ○ | — | — |
| **유체계면 콜로이드** | **✅ [카드](interfacial-colloid--equilibrium-structure.md)** | **✅ [카드](interfacial-colloid--transport.md)** | ○ | ○ | — |
| **브러시 콜로이드** | — | — | — | ○ | — |
| **ABP** | — | ○ | — | — | **✅ [카드](abp--dense-collective.md)** |
| **다공성 매질 속 tracer** | — | ○ | — | — | — |
| **2D 콜로이드 결정** ★신규 | **✅ [소프트 반발 `A/r³`](soft-repulsive-2d--equilibrium-structure.md)** | — | ○ | **✅ [끌기](colloidal-crystal-2d--driven-probe.md)** · **✅ [진동](colloidal-crystal-2d--oscillatory-microrheology.md)** | — |
| **콜로이드 사슬** ★신규 | — | — | — | **✅ [굽힘 강성](colloidal-chain--bending-rigidity.md)** | — |

**현재 9/20.** v1(`D3`: BD + 구형 콜로이드 3D)에 필요한 것은 왼쪽 위 블록이다.

> **계면 콜로이드 행은 v1 범위 밖이지만 카드가 먼저 생겼다.** 랩 논문 11편이 이 계에 몰려 있어
> 파라미터가 실측으로 채워졌기 때문이다 — 카드는 "다룰 때 만든다"가 원칙이나(`D26`),
> **근거가 이미 손에 있는 경우는 예외로 둔다.** 쓰이지 않은 증류가 검증되지 않는다는 우려는
> 여기선 약하다: 값이 논문 표에서 그대로 왔고 출처가 붙어 있다.

---

## 카드 목록

| 카드 | 상태 | 기준 시간 | 근거 |
|---|---|---|---|
| [`passive-sphere--equilibrium-structure`](passive-sphere--equilibrium-structure.md) | `draft` | `τ_D = σ²/D₀` | **표준 문헌** — 랩 대응 논문 없음 |
| [`abp--dense-collective`](abp--dense-collective.md) | `usable` | **`τ_r = 1/D_r`** | 랩 공개 코드 + 논문 2편 |
| [`interfacial-colloid--equilibrium-structure`](interfacial-colloid--equilibrium-structure.md) | `draft` | **없음 (MC)** | Choi 2020 — 파라미터·MC 설정 실측 |
| [`colloidal-crystal-2d--driven-probe`](colloidal-crystal-2d--driven-probe.md) | `draft` | `τ_d = d²/D₀` | 사용자 손그림 `trap-drag-2d-hex300` + Zahn 1999 상도 |
| [`colloidal-crystal-2d--oscillatory-microrheology`](colloidal-crystal-2d--oscillatory-microrheology.md) | `draft` | `τ_d = d²/D₀` (동일) | 사용자 지시 2026-07-28 |
| [`interfacial-colloid--transport`](interfacial-colloid--transport.md) | `draft` | **`τ_D = d²/D₀`**, `γ = 3πη_eff d` | Park & Furst 2008 (`η_eff`) + Choi 2020 |
| [`passive-sphere--harmonic-trap`](passive-sphere--harmonic-trap.md) | **`usable`** | **`τ_trap = γ/k`**, 길이 `ℓ_trap = √(kT/k)` | **사용자 손그림 + 11개 벤치마크 실측 (run 2026-07-28)** |
| [`soft-repulsive-2d--equilibrium-structure`](soft-repulsive-2d--equilibrium-structure.md) | `draft` | `τ_d = d²/D₀`, 길이 **`d = n^{-1/2}`** (최근접거리 아님) | **사용자 손그림 `soft-r3-2d-A-sweep` + 52런 실측 (2026-07-28/29) + Zahn 1999 상도 `[미재현]`** |
| [`colloidal-chain--bending-rigidity`](colloidal-chain--bending-rigidity.md) | `draft` | — 관측량이 `κ` 다 | 사용자 손그림 `chain-bend-2d-oscill` + Furst 계수 |

---

## ★ 카드가 일반 게이트를 뒤집은 사례

카드를 나눠보니 **`05_validation_gates.md`의 게이트가 계-동역학에 무관하다고 가정한 것이
틀렸다**는 게 드러났다.

| 게이트 | 수동 구형 × 평형 구조 | ABP × 조밀 집단 | 계면 × 평형 구조 | 계면 × 수송 | 콜로이드 사슬 × 굽힘 ★ |
|---|---|---|---|---|---|
| 평형화 (`pymbar.detect_equilibration`) | ✅ 유효 | ❌ **무의미** — 능동계는 열평형에 안 간다 | ✅ 유효 | ✅ 유효 | ✅ 유효 (`thermal` 모드) · — (`static`) |
| 자기일관성 `D_msd = kT/γ` | ✅ 성립 | ⚠️ **성립 안 함** — `D_eff = D_t + U₀²τ_r/2` | ❌ **해당 없음 (MC)** | ✅ **필수** — `γ`에 `η_eff` | — 관측량이 `κ` 다 |
| 스텝당 변위 | ✅ `√(2D₀Δt)/σ` | ✅ | ❌ **MC엔 `Δt`가 없다** | ✅ **`√(4D₀Δt)/d` — 2D는 4** | ✅ **결합 상대로 켜진다** (쌍 없음) |
| 힘 변위 `max\|F\|Δt/γσ` | ✅ | ✅ | ❌ | ✅ | ⚠️ **무력** — 곧은 사슬은 `max\|F\| = 0` 정류점 |
| **강성 안정성** `Δt ≤ 0.2·2γ/λ_max` ★ | — 결합 없음 | — | ❌ **`Δt` 없음** | — | ✅ **1급. 거의 항상 지배** |
| 결합 길이 `\|b−b₀\|/b₀ ≤ 5 %` ★ | — | — | — | — | ✅ **필수** — `check_finite` 로는 못 잡는다 |
| 오버댐프 `τ_B/τ_D` | ✅ | ✅ | ❌ **시간축 없음** | ✅ | ✅ |
| 이류 변위 `u₀Δt/σ` | — 해당 없음 | ✅ **필수** | — | — | — |
| 이분산 강제 | `φ > 0.58` 근처 | `φ` 고임계 (1.4, 2/3:1/3) | — (`Ω` 불균질이 대신) | — | — |
| 유한크기 `L` vs `1.5L` | ✅ | ✅ | ⚠️ **불가** — 실험 스냅샷 배치라 `N` 고정 | ✅ | ⚠️ `N` 스윕이 관측량이다 (`κ ∝ L⁻³`) |
| `r_cut` 절단 | `r_cut ≤ L/2` | — | ⚠️ **`r_cut/d ≥ 21–69` 필요. 박스가 못 담을 수 있다** | ⚠️ 동일 | — 쌍 포텐셜 없음 |

> ★ **사슬 열은 게이트 두 개를 새로 요구했다** (2026-07-28). 정확도 게이트 셋만으로는
> 결합 계의 `Δt` 를 정할 수 없다 — 곧은 사슬은 힘이 정확히 0인 정류점이라 힘 게이트가
> 무력해지고, `kT = 0` 에서도 터진다. 그리고 그 폭발이 **조용하다** (`1.4×10⁷` 도 유한하다).
> 근거: [`dt-gate-needs-a-stability-term-for-stiff-bonds`](../findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md)
> · [`dead-end-distance-constraint-with-brownian`](../findings/dead-end-distance-constraint-with-brownian.md)
>
> ⚠ **이 열에는 아직 카드가 없다.** `scripts/chain_bend.py` 가 카드 없이 돌고 있다 —
> 즉흥 무차원화 금지 규약 위반이다. 아래 우선순위 표 1번.

> **계면 두 카드가 이 표의 논지를 가장 세게 보여준다.** 같은 계·같은 퍼텐셜인데
> **엔진이 MC냐 BD냐로 게이트 네 개가 켜지고 꺼진다.** MC에는 `Δt`도 시간축도 없으므로
> 변위·오버댐프·자기일관성 게이트가 통째로 무의미하다 — 그런데 이걸 모르고 일반 게이트를
> 걸면 **통과할 수 없는 검사**를 요구하게 된다.

> **`05_validation_gates.md`는 게이트를 쌍별로 켜고 끄는 구조로 써야 한다.**
> 지금 계획은 5범주를 모든 런에 똑같이 적용하게 되어 있는데, 그러면
> 능동계에서 **통과할 수 없는 게이트를 걸거나, 성립하지 않는 검사를 통과시킨다.**

## 다음에 만들 카드

| 우선 | 카드 | 왜 |
|---|---|---|
| **1** | **`colloidal-chain--bending`** | ★ **이미 돌고 있다.** `scripts/chain_bend.py` 가 카드 없이 자체 축약단위(`σ = b = 1`, `τ_D`)를 쓰고, 게이트 2개(강성 안정성·결합 길이)를 새로 요구했다. `κ ∝ L⁻³` 벤치마크와 `k_bond* = C·κ(N)*` 비용 규칙이 여기 붙는다. 큐의 `chain-bend-2d-oscill` 이 이 카드를 거쳐야 `G'(ω)` 로 간다 |
| 2 | `attractive-colloid--coarsening` | 조대화 지수 `t^{1/3}` 벤치마크가 여기 붙는다 |
| 3 | `passive-sphere--transport` | 자유 BD `D = kT/γ` — M1의 표적. `CARD_SCALE_RULES` 에 등록만 되어 있고 카드가 없다 |
| 4 | `abp--dilute-transport` | Modica 다공성 매질 · Barakat 트랩 |
| 5 | `brush-colloid--nonequilibrium-contact` | Xu 2023 — `Δt` 실측 표본 |


---

## ★ 2026-07-28 — 카드의 단위는 "척도가 다른 쌍"이 아니다

새로 만든 두 카드([`끌기`](colloidal-crystal-2d--driven-probe.md) ·
[`진동`](colloidal-crystal-2d--oscillatory-microrheology.md))는 **기준 단위가 같다**
(`d`, `kT`, `τ_d = d²/D₀`). 이 색인 서두의 정당화("목적동역학마다 기준 시간이 달라진다")가
여기서는 성립하지 않는다.

**그래도 카드를 나누는 것이 옳다.** 갈리는 것은 척도가 아니라:

| | 등속 끌기 | 진동 |
|---|---|---|
| 지배 무차원수 | `Pe_d = vd/D₀` | `ωτ_trap`, `a/d` |
| 정상상태 | 비평형 정상상태 | **주기적** 정상상태 |
| 관측량 | 항력·일·결함 축적 | `G'`, `G''` |
| 1급 게이트 | 힘 변위 · 결함 판정 | **선형성 (`a` 스윕)** · `ωΔt ≪ 1` |

⇒ **카드가 소유하는 것 중 실제로 갈리는 것은 §7 게이트와 §6 관측량이다.**
기준 단위는 갈릴 수도 있고 안 갈릴 수도 있다. 트랩 카드에서 `τ_trap` vs `τ_D` 가
24만 배 갈린 것은 **강한 사례였을 뿐 정의는 아니다.**

`master_plan.md` §5·§10.3 의 서술이 척도 차이만 강조하고 있어 이 사례로 넓혀야 한다.
