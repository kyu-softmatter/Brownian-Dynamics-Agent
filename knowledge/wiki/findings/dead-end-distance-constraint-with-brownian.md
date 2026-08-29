---
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
question: "hoomd.md.constrain.Distance 를 md.methods.Brownian 과 함께 쓸 수 있는가?"
verdict: NO
reproduced: yes
---

# 막힌 길 — `constrain.Distance` × `methods.Brownian`

> **결론: 못 쓴다. 그리고 조용히 틀린다** — 예외를 던지지 않고 결합 길이가 발산하는데
> `check_finite` 가드는 **통과한다.**

## 왜 시도했나

콜로이드 사슬 굽힘(Pantina & Furst 재현)에서 접촉 결합의 강성이 열적 스케일보다
`κ₀/(kT/σ²) ≈ 3.4×10⁷` 배 크다. 조화 결합으로 이 강성을 표현하면 변위 기준 `dt` 게이트가
`Δt* ~ 10⁻⁸ τ_D` 를 요구해 비용이 폭발한다. **신축 모드를 아예 제거**하면 `dt` 가
굽힘 모드에만 묶이므로 거리 구속이 자연스러운 해법으로 보였다.

## 실측 (2026-07-28, hoomd 7.1.0, 2D, N=5 직선 사슬, `dt*=1e-5`, 3000 step)

| 결합 처리 | 초기 `d` | 3000 step 후 `d` | `max|d−1|` | `check_finite` |
|---|---|---|---|---|
| `bond.Harmonic` `k*=1e4` | `1,1,1,1` | `0.997 … 1.008` | `1.4e-2` | pass |
| **`constrain.Distance`** | `1,1,1,1` | **`3.3e7, 5.8e7, 4.3e7, 8.5e6`** | **`5.8e7`** | **pass** |
| `bond.Harmonic` `k*=1e2` | `1,1,1,1` | `0.888 … 1.029` | `1.1e-1` | pass |

HOOMD 는 죽지 않고 경고만 반복한다:

```
*Warning*: Constraint 3 between particles 3 and 4 violated!
(distance 18.5696 exceeds 1 within relative tolerance 0.001)
```

## 왜 이렇게 되는가

`Distance` 는 Yoneya 행렬법(비반복 구속 MD)이다 — **가속도에서 Lagrange 승수를 푼다.**
`Brownian` 은 과감쇠 적분기라 속도·가속도를 쓰지 않는다(질량은 무시된다). 따라서
구속력이 물리적으로 무의미한 값으로 계산되어 **매 스텝 에너지를 주입한다.**
`Langevin`(관성 있음)에서는 다를 수 있으나 이 프로젝트의 계는 전부 과감쇠다.

## 그래서 무엇을 하는가

**조화 결합을 쓰고, 강성을 물리값이 아니라 "굽힘을 오염시키지 않는 최소값"으로 고른다.**
필요 조건은 `k_bond* ≫ κ(N)*` 이고, 그 여유가 얼마여야 하는지는
`cli.py converge` 로 `k_bond*` 를 흔들어 답이 변하지 않음을 확인해서 정한다.
물리 강성 `κ₀* ≈ 3.4×10⁷` 을 그대로 넣는 것은 **목표가 아니다** — 지수 `−3` 은
`κ₀` 에 무관하므로 접근 가능한 강성에서 검증하고 `κ₀` 는 계수로 환산한다.

## ★ 가드 교훈 — 이게 이 기록의 본론이다

`check_finite` 만으로는 이 실패를 못 잡는다. `5.8e7` 은 유한하다.
**결합 계에는 "결합 길이가 목표값 근처인가" 가드가 따로 필요하다.**

```
assert max|d_ij − d_target| / d_target  <=  tol       # 결합 계 전용 게이트
```

없으면 폭발한 사슬에서 `κ` 를 측정해 그럴듯한 `s⁻³` 을 보고할 수 있다 —
[[hoomd-brownian-scheme-and-noise]] 의 "요동하지 않는 측정값" 과 같은 종류의 함정이다.

관련: [[dt-gate-should-be-displacement-based]]
