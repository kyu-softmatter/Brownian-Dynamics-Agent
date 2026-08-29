# S3 · S4 · S5 — 명세 · 무차원화 · 실행

> 이 세 단계는 **거의 전부 코드**다. 너의 일은 `spec.yaml` 을 쓰는 것과
> 게이트 결과를 읽는 것이다. `cli.py run` 이 S3→S8 을 한 번에 돌린다.

## S3 — `spec.yaml`

### 가장 빠른 길: 예시를 복사한다

```bash
cp examples/trap-2d-5um/spec.yaml runs/<topic>/spec.yaml
```

[`examples/trap-2d-5um/spec.yaml`](../../../../examples/trap-2d-5um/spec.yaml) 은
provenance 18필드가 채워진 검증된 명세다.

### 모든 값에 `provenance` + `basis`

```yaml
eta_si:
  value: 8.5566e-4
  unit: Pa*s
  provenance: assumed        # 값은 정확하나 "물이다"라는 선택이 가정
  basis: IAPWS 표 보간 (외삽 아님). knowledge/wiki/concepts/water-298k.md
  confidence: medium
  affects: [tau_trap, D0]    # <x^2> 에는 영향 없음
```

| `provenance` | 언제 | 누가 쓸 수 있나 |
|---|---|---|
| `observation` / `from_drawing` | 자료에서 직접 읽음 | 누구나 (Haiku 포함) |
| `derived` | 다른 필드에서 계산 | 누구나 |
| `rule` | 정책에서 유도 (`run_policy.yaml`) | 누구나 |
| `from_knowledge` / `from_paper` | wiki·문헌 증류 | 누구나 |
| **`inference`** | 자료 + 물리지식으로 유도 | **Opus 만** |
| **`assumed`** | 자료에 없어 채움 | **Opus 만** |
| `user` | 사람이 이 런에서 지정 (`session set`) | 코드가 자동 표시 |
| `measured` | 실험값 | — |

`inference`·`assumed` 는 `confidence` 가 **필수**다. 없으면 `validate` 가 잡는다.

### ⚠ YAML 1.1 함정 — 지수 표기

```yaml
value: 5e-3      # ❌ 문자열이다! 만티사에 소수점이 없다
value: 5.0e-3    # ✅ float
value: 6.3e6     # ❌ 문자열. 지수부 부호가 없다
value: 6.3e+6    # ✅ float
```

**두 번 겪었다** (처리량 상수, `session set`). 전문:
[`yaml-scientific-notation-parsed-as-string`](../../../../knowledge/wiki/findings/yaml-scientific-notation-parsed-as-string.md)

`simbot.spec.dump_yaml` 로 쓰면 안전하다 (파이썬 `float` 이 항상 소수점을 붙인다).
**손으로 쓴 YAML 이 위험하다.**

### 게이트 선언 — 카드가 켜고 끈다

```yaml
gates:
  equipartition: {status: required, reason: 이 카드의 1급 게이트}
  step_displacement_vs_sigma:
    status: off
    reason: 쌍 상호작용이 없어 겹침이 발생할 수 없다      # ★ off 에는 이유 필수
```

- 게이트 이름은 **등록된 것만** (`simbot.spec.KNOWN_GATES`). 오타는
  **한 번도 실행되지 않는 검사**가 된다 → `validate` 가 거부한다
- 어느 게이트가 켜지는지는 카드 §7 이 정한다
- `required` 는 **결과가 아니라 선언**이다. S3 리포트에 `⏳ S7 판정` 으로 나온다

### 파생값은 저장하지 않는다

`kT_si`, `tau_trap_si` 같은 파생값을 `spec.yaml` 에 적지 않는다 —
`simbot.spec.derive()` 가 계산한다. 적혀 있으면 `validate(stored_derived=...)` 가
재계산과 대조해서 불일치를 잡는다 (손으로 고친 파생값을 잡는 유일한 방법).

### 검사

```bash
<PY> -c "
from simbot.spec import SystemSpec, validate
r = validate(SystemSpec.load('<spec.yaml>'))
print(r.table())
print()
print('규약 위반:', r.problems or '없음')
print('S7 판정 대기:', [c.name for c in r.deferred()])"
```

---

## S4 — 무차원화

**전부 자동이다.** 너가 할 일은 결과를 읽는 것.

```bash
<PY> -c "
from simbot.spec import SystemSpec
from simbot.nondim import reduce_spec, roundtrip_errors, nondim_table
sp = SystemSpec.load('<spec.yaml>')
r = reduce_spec(sp)
print(nondim_table(sp, r))
print('왕복오차:', max(roundtrip_errors(sp, r).values()))
print('dt* =', r.dt_star, '지배 제약:', r.dt_dominant)"
```

### 확인할 것 셋

1. **척도 출처** — `scales_harmonic_trap: (l_trap, kT, tau_trap)` 처럼 카드가 나와야 한다
2. **왕복오차 `< 1e-12`** — 크면 계산 실수가 아니라 **규약 위반**이다
   (예: `τ_D` 로 나누고 `τ_trap` 으로 되돌리기)
3. **지배 제약** — 어느 `dt` 제약이 이겼는지. `spec(명시값)` 이면 사람이 정한 것

### `dt` 제약은 계마다 다르게 켜진다

| 제약 | 켜지는 조건 |
|---|---|
| 열 변위 `√(2D₀Δt) ≤ 0.03σ` | **쌍 상호작용 있을 때만** |
| 힘 변위 | 쌍 상호작용 + `max|F|` **실측** (추정 금지) |
| 완화시간 `Δt ≤ 0.01 τ` | 구속·활성 있을 때 |
| 활성 변위 | 능동 구동 있을 때 |
| 정확도 목표 | 조화 트랩 + 목표 편향 명시 |

★ **변위 게이트는 만능이 아니다.** 트랩 계에서 변위 상한이 완화시간 상한의 **1086배**다 —
게이트가 아무것도 막지 못한다. 전문:
[`displacement-gate-is-1000x-loose-for-traps`](../../../../knowledge/wiki/findings/displacement-gate-is-1000x-loose-for-traps.md)

`dt/τ_D` 는 **기록만** 한다. 게이트로 쓰면 논문까지 나온 런을 기각한다.

### 카드가 없으면 예외가 난다

```
KeyError: 카드 'colloid--new-thing' 의 척도 규칙이 등록되지 않았다.
즉흥 무차원화는 금지다 — _TEMPLATE.md 로 draft 카드를 먼저 만들고
CARD_SCALE_RULES 에 등록할 것
```

**우회하지 말고 카드를 만든다.** 카드 없이 돌린 결과는 나중에 재현할 수 없다.

---

## S5 — 실행

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>
```

### CLI 가 실행 **전에** 멈추는 조건

| 조건 | 메시지 |
|---|---|
| S3 게이트 실패 | `S3 게이트 미통과` |
| 왕복오차 `≥ 1e-12` | `S4 왕복 게이트 위반 — 척도 규약이 어긋났다` |
| 시드 `< 4` | `오차 막대 없는 프로덕션 런은 금지다` |
| 예산 초과 예상 | `예산 초과 예상 — 실행하지 않고 보고한다` |

**`--force` 는 사용자가 명시적으로 요구할 때만.** 게이트를 우회한 결과에는
그 사실을 리포트에 적는다.

### 티어 사다리

처음 보는 카드면 CLI 가 경고한다:

> ⚠️ 이 카드의 첫 런이다. 정책상 사다리 `[smoke, pilot, explore]` 를 건너뛸 수 없다 —
> 이 런이 사다리의 첫 칸이 된다. **결과를 프로덕션으로 인용하지 말 것**

### 오차 막대는 공짜다

`k ≤ 4` 효율이 93 % 이므로 **시드 4개 비용 ≈ 1개**. `run_policy.yaml` 이 최소 4개를
강제한다. **긴 런 1개보다 짧은 런 4개.**

### 실패한 런은 버리지 않는다

배치에서 일부가 죽으면 `05_run_manifest.json` 의 `batch.failed` 에 남고 CLI 가 보고한다.
조용히 빠지면 "시드 4개"라고 적힌 오차막대가 실제로는 3개짜리가 된다.

### 파라미터를 흔들어 보기 — 실행 없이

```bash
<PY> -m simbot.session new <spec.yaml>
<PY> -m simbot.session set numerics.dt_star=2.5e-3 species.0.n_simulated=4000
<PY> -m simbot.session show
```

`set` 은 **비용 추정만** 한다. 실행은 `cli.py run` 이 한다.
바뀐 값은 `provenance: user` 로 표시되고 **이전 값과 원래 근거가 basis 에 남는다.**

### 수렴 확인

```bash
<PY> cli.py converge <spec.yaml>
```

`dt/2`, `dt×2`, `N×2`, 시드 이동을 돌려 **통계오차 대비**로 판정한다.
`3σ` 이내면 "구별 안 됨" — **같다는 증명이 아니라 이 오차로는 차이를 볼 수 없다는 뜻이다.**
