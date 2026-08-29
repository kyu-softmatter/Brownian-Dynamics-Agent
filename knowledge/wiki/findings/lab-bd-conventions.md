---
type: finding
author: agent
drafted: 2026-07-27
confirmed_by:
question: "우리 랩은 BD 시뮬레이션을 실제로 어떤 파라미터로 돌려왔는가?"
cites:
  - knowledge/source/papers/2023-xu-dynamic-interfaces-contact-time.md
  - knowledge/source/papers/2024-xu-dynamic-surfactants-anisotropic-assembly.md
  - knowledge/source/papers/2022-barakat-enhanced-dispersion-harmonic-traps.md
  - knowledge/source/papers/2020-takatori-motility-induced-buckling.md
  - knowledge/source/papers/2025-quah-continuum-closures-active-control.md
  - knowledge/source/papers/2022-modica-porous-media-active-diffusion.md
  - knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md
related_systems:
  - abp--dense-collective
  - passive-sphere--equilibrium-structure
cross_cut: true   # 계-동역학을 가로지르는 감사 문서. 1차 조직은 systems/ 카드
---

# 우리 랩의 BD 관행 — 발표 논문에서 추출한 파라미터 사전

**질문:** 새 BD 시뮬레이션을 설계할 때, 우리 랩이 실제로 써온 값은 무엇인가?

**답:** 발표 논문 6편 + 공개 코드 1건에서 추출했다.

> ### ⚠️ 이 페이지는 **가로지르기(cross-cut) 감사**지 1차 조직이 아니다
>
> 계-동역학 쌍을 무시하고 랩 전체를 한 축으로 훑은 문서다. 그 덕에 **`dt` 게이트 결함을
> 찾아낼 수 있었지만**, 무차원화 규약은 쌍마다 정당하게 다르므로 이 페이지를 사전(prior)으로
> 직접 쓰면 안 된다.
>
> **파라미터를 찾으려면 → [`../systems/_index.md`](../systems/_index.md) 의 해당 카드를 볼 것.**
> 이 페이지는 "랩 전체에 걸친 관행과 그 편차"를 보는 용도다.

> ⚠️ **전부 `reproduced: no`다.** 논문에 적힌 값을 옮긴 것이지 우리 코드에서 돌려본 것이 아니다.
> 검증 근거가 아니라 **출발점**으로만 쓴다. 계약: [`../CLAUDE.md`](../CLAUDE.md)

---

## 1. 엔진 — HOOMD-blue가 표준이다

확인된 6편 중 **BD를 돌린 논문은 전부 HOOMD-blue**를 썼다. LAMMPS 사용 사례는 없다.

| 논문 | 엔진 | 비고 |
|---|---|---|
| **Quah 2025 (공개 코드)** | **HOOMD 3.8.1** | `environment.yml`에 버전 명시. GPU 우선, CPU 폴백 |
| Xu 2023 (dynamic interfaces) | HOOMD-blue | **GPU 가속** 명시 |
| Xu 2024 (dynamic surfactants) | HOOMD-blue | **GPU 가속** 명시 |
| Barakat 2022 (harmonic traps) | HOOMD-blue | N = 10,000 비상호작용 |

> **함의:** 우리 `D2`(로컬 CPU 우선)와 랩 관행 사이에 간극이 있다. 랩은 GPU를 쓴다.
> Xu 2024의 **8×10⁸ 스텝**은 CUDA 없는 M4에서 재현 불가다.
> `master_plan.md` §8의 비용 게이트가 이런 런을 반드시 걸러야 한다.

## 2. 시간간격 — `τ_D`로 무차원화한다 ★

**확인된 값 2건.** 하나는 논문 SI, 하나는 **공개 코드**에서 나왔다.

```
Δt = 2×10⁻⁵ · (σ²/D)
```

`σ²/D`는 입자가 자기 지름만큼 확산하는 브라운 시간 — **우리가 `D7`에서 기준 시간으로 채택한
`τ_D`와 같은 양**이다.

두 번째 값은 Quah의 공개 코드 `graybox_abp_mpc` (`md_dt = 5e-4`, `τ_r` 단위):

| 출처 | `Δt/τ_D` | 스텝당 확산변위 |
|---|---|---|
| 선행 slit 프로젝트 | `1.0e-3` | 0.045 σ |
| Xu 2023 | `2.0e-5` | 0.0063 σ |
| **Quah 코드** | **`1.67e-4`** | 0.018 σ |
| 우리 게이트 | `≤ 1e-4` | — |

> ### ★ 우리 게이트가 3건 중 2건을 기각한다
>
> 선행 slit 프로젝트와 Quah 코드가 모두 `1e-4`를 넘는다. **실제로 돌아가고 논문까지 나온
> 시뮬레이션을 거부하는 게이트는 잘못 맞춰진 것이다.**
>
> 스텝당 변위로 재면 셋 다 0.006–0.045σ 안에 든다 (7배 폭 vs `dt/τ_D`의 50배 폭).
> → **[[dt-gate-should-be-displacement-based]]** 에서 게이트 교체를 제안한다.
>
> `D7`(기준 시간 = `τ_D`) 자체는 유지된다 — 뒤집히는 것은 **게이트지 단위가 아니다.**

**표기 차이 주의** — Xu 2024는 같은 양을 `t* = 4R_c²/D_c`로 쓴다. `2R_c = σ`이므로
`4R_c² = σ²` — **동일하다.** `UnitMap`에 이 환산을 명시해둘 것.

## 3. 배제부피 — WCA를 쓴다 (`D8` 관련)

확인된 4편이 **WCA**를 썼다. Xu 2023은 `ε = α = k_BT`로 명시.

```
U₂(r) = 4ε[(a/r)¹² − (a/r)⁶] + ε,   r < 2^{1/6}a        (Quah 2025)
```

> **`D8` 판단에 중요하다.** 선행 프로젝트 `brownian_slit_sim/src/forces.py:117`은
> *"WCA의 `r⁻¹³` 코어는 오버댐프에서 위험하다 — 작은 겹침도 입자를 박스 밖으로 날린다"*
> 고 경고하며 bounded harmonic을 권했다.
>
> 그런데 랩 논문들은 WCA를 문제없이 쓴다. **모순이 아니다** — Xu 2023의 `Δt = 2e-5 τ_D`가
> 선행 프로젝트보다 충분히 작기 때문일 가능성이 크다.

**Quah 코드가 결정적 증거를 준다:** `ε/kT = 500`이라는 매우 강한 WCA를 쓰고도 안정적으로 돈다.
선행 프로젝트보다 `dt/τ_D`는 6배 작지만 **스텝당 변위가 0.018σ 대 0.045σ로 2.5배 작다.**

> **가설 (수정):** WCA의 위험성은 퍼텐셜 자체도, `Δt`도 아니라 **스텝당 변위**와의 결합이다.
> 변위를 0.02σ 이하로 두면 `ε/kT = 500`도 안정적이다.
> → 재현 실험으로 확정할 것. [[dt-gate-should-be-displacement-based]] 참조.

## 4. 결정화 억제 — 이분산 (`D10` 관련) ★

Takatori 2020 (조밀 능동계, `φ` 최대 0.83):

| 항목 | 값 |
|---|---|
| 지름 비 | **1.4** |
| 몰분율 | **2/3 : 1/3** |
| 목적 | *"To avoid interference by crystallization"* |

> **`D10`(v1 단분산) 갱신 제안:** 단분산을 기본으로 유지하되,
> **`φ`가 높아지면 이분산을 강제하는 게이트**를 `05_validation_gates.md`에 넣을 것.
> 안 그러면 의도치 않게 결정화된 계를 "유리"라고 부르게 된다.
> 문턱값은 미정 — 랩 논문은 `φ = 0.83`에서 썼다.

## 5. 통계 — 독립 시드 앙상블 (`D16` 관련)

| 논문 | 방식 |
|---|---|
| Xu 2023 | **20 realization 평균**, 10³ 스텝마다 샘플링 |
| Barakat 2022 | N = 10,000 입자 (비상호작용 → 입자 자체가 앙상블) |

> 우리 `D16` 기본값은 **block averaging**이다. 랩 관행은 **독립 시드 앙상블**에 가깝다.
> **둘 다 필요하다** — 시드 앙상블은 시드 간 분산을, block averaging은 시계열 내
> 자기상관을 잡는다. `07_observables.md`에서 둘을 모두 요구할 것.

## 6. 무차원수 — 랩이 실제로 쓰는 것

| 무차원수 | 정의 | 출처 |
|---|---|---|
| `Pe` | `U₀τ_R/σ` = `U₀/(D_Rσ)` = `τ_R/(σ/U₀)` — **세 가지로 읽는다** | Takatori 2020 |
| **`(ℓ/δ)²`** | `U₀²τ_R/D_T` — 런 길이 / 열확산 길이 | Quah 2025 |
| `φ` (2D) | `nπσ²` (Takatori) · `n̄πa²/4` (Quah) — **정의가 다르다** | 여러 편 |
| `φ_eff` | 복합 입자의 폴리머 부피 보정 | Xu 2024 |
| 능동력/구속력 | `ζU₀ / (구속력)` — 계면·트랩 탈출 판정 | Cheon 2024 |

> ⚠️ **`φ` 정의가 논문마다 다르다.** `nπσ²`(지름 기준)와 `n̄πa²/4`(다른 규약)가 섞여 있다.
> `bdkit/units/`는 **φ를 계산할 때 지름/반지름 규약을 항상 기록**해야 한다.
> 이건 조용히 4배 틀릴 수 있는 지점이다.

---

## 이 발견이 바꾸는 것

| 문서 | 반영할 내용 |
|---|---|
| `master_plan.md` §5 S3 | ★ **게이트 교체** — `Δt/τ_D ≤ 1e-4` → 스텝당 변위 기준. [[dt-gate-should-be-displacement-based]] |
| `docs/00_decision_log.md` `D7` | 기준 시간 `τ_D`는 **유지**. 뒤집히는 것은 게이트지 단위가 아니다 |
| `docs/00_decision_log.md` `D8` | WCA 위험성 = **스텝당 변위**와의 결합. `ε/kT=500`도 변위 0.02σ면 안정 |
| `docs/00_decision_log.md` `D10` | 조밀계에서 이분산 강제 게이트 필요 (1.4, 2/3:1/3) |
| `docs/00_decision_log.md` `D16` | block averaging **과** 시드 앙상블 둘 다 |
| `docs/03_units_nondim.md` | `(ℓ/δ)²` 추가 · `t*↔τ_D` 환산 · **φ 규약 명시 강제** |
| `docs/05_validation_gates.md` | `φ` 고임계에서 이분산 게이트 |

## 남은 빈칸

- **`Δt` 명시 표본이 2개뿐이다** (Xu 2023 SI, Quah 코드) → paywall 논문 확보 후 재조사
- **HOOMD 버전이 3.8.1 대 우리 7.1.0.** API 이식 필요 — `gamma.default` 문법 등
- Barakat 2022의 BD `Δt`·스텝 수 — SI에 있다는데 arXiv ancillary에 미포함
- Modica 3편(BD 직결)이 전부 paywall → **가장 아쉬운 공백**
