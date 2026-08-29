---
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
status: proposed
question: "변위 기준 `dt` 게이트를 모든 계에 걸면 되는가?"
answer: "안 된다. 트랩 계에서는 변위 게이트 상한이 완화시간 상한의 1086배다 — 게이트가 아무것도 막지 않는다. 변위 게이트는 쌍 상호작용 계의 게이트이고, 구속 계의 게이트는 완화시간이다."
cites:
  - knowledge/wiki/findings/dt-gate-should-be-displacement-based.md
  - knowledge/wiki/systems/passive-sphere--harmonic-trap.md
  - simbot/nondim.py
runs: [2026-07-28_trap-2d-5um_2dfb9d]
affects_docs: [CLAUDE.md, config/run_policy.yaml, knowledge/wiki/systems/_index.md]
---

# 변위 게이트는 트랩 계에서 1086배 느슨하다 — 게이트의 적용 범위

## 질문

[[dt-gate-should-be-displacement-based]] 는 `dt/τ_D` 고정 게이트를 기각하고
**스텝당 변위** 기준을 채택했다. 그 결론은 유효하다. 그런데 **변위 게이트를 모든 계에
걸어도 되는가?**

## 답 — 안 된다. 계마다 진짜 제약이 다르다.

첫 손그림 계(`a = 5 μm`, `k = 10 pN/μm`, `T = 300 K`, 물, 2D)에서 네 제약을 모두
SI 로 계산해 같은 축에 올렸다 (`simbot.nondim.choose_dt`):

| 제약 | `Δt` 상한 [s] | `Δt*` 상한 (`τ_trap` 단위) | 이 계에서 |
|---|---|---|---|
| 열 변위 `√(2D₀Δt) ≤ 0.03 σ` | `0.876` | **`108.6`** | 무의미 |
| 힘 변위 `max\|F\|Δt/γ ≤ 0.005 σ` | — | — | 쌍 힘 없음 |
| **완화시간 `Δt ≤ 0.01 τ_trap`** | **`8.06e-5`** | **`0.01`** | **← 지배** |
| 활성 변위 | — | — | 능동 없음 |

```
변위 상한 / 완화 상한 = 108.6 / 0.01 = 1086배
```

**변위 게이트만 켰다면 `Δt*` 를 108 까지 허용한다.** 안정성 한계가 `Δt* < 2` 이므로
그 값은 발산하지만, `Δt* = 1` 같은 값은 **발산하지 않고 통과한다** — 완화 과정을
통째로 건너뛰면서.

## 왜 이렇게 벌어지는가

두 게이트가 서로 다른 길이를 기준으로 삼기 때문이다.

```
변위 게이트의 기준 길이 = σ        (= 10 μm,  겹침을 막는 척도)
구속 계의 실제 기준 길이 = ℓ_trap  (= 20.35 nm, 입자가 탐색하는 척도)

σ / ℓ_trap = 491
```

입자는 **자기 지름의 0.2 %** 만 움직인다. "지름의 3 % 이상 움직이지 마라"는 제약은
탐색 영역의 15배를 허용하는 셈이다. 변위 게이트는 **겹침을 막는 게이트**이고,
겹칠 상대가 없으면 막을 것이 없다.

## 두 게이트는 경쟁하지 않는다 — 보완한다

| | 진짜 제약 | 왜 |
|---|---|---|
| **쌍 상호작용 계** (WCA·LJ·Yukawa) | **변위** | 한 스텝에 코어를 관통하면 폭발한다 |
| **구속 계** (트랩·조화) | **완화시간** | `Δt > τ_relax` 면 동역학을 건너뛴다 |
| **능동계** | 변위 + `1/D_r` | 이류와 회전 완화가 각각 걸린다 |

원 finding 의 실측 3건(선행 slit / Xu 2023 / Quah 코드)은 **전부 쌍 상호작용 계**였다.
그래서 그 결론의 적용 범위는 "쌍 상호작용 계"이고, 트랩 계로 외삽하면 안 된다.

⇒ **`(계 × 목적동역학)` 카드가 게이트 표를 소유하는 이유가 정량적으로 확인됐다.**
게이트를 전역 상수로 두면, 어느 계에서든 하나는 반드시 헐거워진다.

## 처방

`simbot.nondim.choose_dt` 는 제약마다 `active` 플래그와 **끈 이유**를 갖는다.
꺼진 제약도 상한을 계산해 표에 남긴다 — 그래야 "얼마나 헐거웠는지"가 보인다.

```python
DtConstraint(
    name="thermal_displacement",
    dt_si_max=(0.03 * sigma) ** 2 / (2 * D0),
    active=has_pair,                     # ← 쌍 상호작용이 있을 때만
    off_reason="쌍 상호작용이 없어 겹칠 상대가 없다 — 변위 상한이 무의미")
```

활성 제약이 하나도 없으면 **기본값을 쓰지 않고 예외를 던진다.** `dt` 를 고를 근거가
없는 상태로 런이 나가는 것이 가장 나쁘다.

## 부수 확인 — 정확도 목표는 별도 제약이다

트랩 계에서는 "정확도를 얼마나 원하는가"가 `Δt` 를 정할 수 있다. Euler–Maruyama
편향이 해석적으로 알려져 있으므로 (`1/(1−Δt*/2) − 1`), 목표 편향 `0.25 %` 를 넣으면
`Δt* = 0.004988` 이 나온다.

**첫 런에서 사람이 고른 `Δt* = 5e-3` 과 0.24 % 이내로 일치한다** — 사람의 선택과
코드의 유도가 독립적으로 같은 값에 도달했다.

## 재발 방지

| 장치 | 위치 |
|---|---|
| 1086배 비를 고정 | `tests/test_s4_nondim.py::test_displacement_gate_is_1000x_too_loose_for_a_trap` (`benchmark` 마커) |
| 쌍 상호작용 추가 시 게이트가 켜지는지 | `::test_displacement_gate_turns_on_with_pair_interactions` |
| 근거 없이 dt 를 고르지 않는지 | `::test_choose_dt_refuses_when_no_constraint_applies` |
| `max\|F\|` 추정 금지 | `::test_force_constraint_is_na_without_measured_force` |
| τ_D 규약의 파국 재현 | `::test_universal_tau_D_convention_would_be_catastrophic` |

## 적용 범위 / 한계

- 이 1086배는 **`k*_σ = kσ²/kT = 2.41e5`** 인 계의 값이다. 약한 트랩(`k*_σ ~ 1`)에서는
  `ℓ_trap ~ σ` 가 되어 두 게이트가 만난다. 비는 대략 `√(k*_σ)` 로 커진다.
- 힘 변위 제약은 `max|F|` 를 **실제로 계산**해야 켜진다 (§5.4: 추정 금지). 아직 초기배치
  힘 계산 경로가 없어 `n/a` 로 남는다 — 트랩+WCA 계에서 처음 필요해진다.

## 참고

- [[dt-gate-should-be-displacement-based]] — 이 finding 이 적용 범위를 좁힌다 (기각 아님)
- [[hoomd-brownian-scheme-and-noise]] — EM 편향의 해석식
- `runs/2026-07-28_trap-2d-5um_2dfb9d/08_conclusion.md` §3 — `τ_D` 규약이 24만 배 틀리는 실측
