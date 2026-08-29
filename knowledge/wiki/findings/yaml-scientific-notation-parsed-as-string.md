---
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
status: proposed
cause_class: environment
stage: S3
question: "YAML 로 읽은 `6.3e6` 이 왜 숫자가 아닌가?"
answer: "YAML 1.1 의 float 정규식은 만티사에 소수점을, 또는 지수부에 부호를 요구한다. `6.3e6`·`5e-3` 은 문자열이 된다. PyYAML 은 YAML 1.1 이므로 `safe_load` 를 숫자 파싱에 그대로 쓸 수 없다."
cites:
  - config/run_policy.yaml
  - simbot/policy.py
  - simbot/session.py
affects_docs: [config/run_policy.yaml, CLAUDE.md]
---

# `6.3e6` 은 숫자가 아니다 — YAML 1.1 지수 표기 함정

## 증상

**하루에 서로 다른 두 곳에서 같은 증상이 났다.**

### ① `config/run_policy.yaml` 의 처리량 상수

```python
>>> yaml.safe_load(open('config/run_policy.yaml'))['hardware']['throughput_particle_steps_per_s']
'6.3e6'          # ← str. float 이 아니다
```

비용 추정 `wall ≈ N·steps/(Λ·η_k)` 가 `Λ` 로 문자열을 받는다.
`estimators.estimate_wall_time_s` 는 모듈 상수를 쓰기 때문에 살아 있었고,
**정책 파일 경로로만 조회하면 조용히 깨지는 상태**로 M0 내내 있었다.

### ② `session set numerics.dt_star=5e-3`

```python
>>> yaml.safe_load('5e-3')
'5e-3'           # ← str
>>> yaml.safe_load('2.5e-3')
0.0025           # ← float. 소수점이 있으면 된다
```

`dt_star` 에 문자열이 들어가 `reduce_spec` 이 `TypeError` 로 죽었다.
증상이 "비용 추정 실패"로 나타나서 원인이 두 단계 떨어져 있었다.

## 진단 경로

1. **값이 틀렸나** — 아니다. 문자열의 내용은 정확히 `'6.3e6'` 이다.
2. **YAML 파일이 잘못됐나** — YAML 1.2 라면 유효한 float 이다. 파일이 아니라 **파서**가 원인이다.
3. **PyYAML 버전 문제인가** — 아니다. PyYAML 은 설계상 YAML 1.1 이고 이것은 명세대로의 동작이다.
4. **어느 표기가 통과하는가** — 실측:

| 표기 | `yaml.safe_load` | 이유 |
|---|---|---|
| `6.3e6` | `'6.3e6'` (str) | 지수부에 부호 없음 |
| `6.3e+6` | `6300000.0` | ✅ 부호 있음 |
| `5e-3` | `'5e-3'` (str) | **만티사에 소수점 없음** |
| `2.5e-3` | `0.0025` | ✅ 소수점 있음 |
| `5.0e-3` | `0.005` | ✅ 소수점 있음 |

## 근본 원인

`cause_class: environment`. PyYAML 의 float resolver 정규식이 요구하는 것:

```
[-+]? [0-9][0-9_]* \. [0-9_]* ([eE][-+][0-9]+)?      # 만티사에 소수점 필수
| \. [0-9_]+ ([eE][-+][0-9]+)?                        # 그리고 지수부 부호 필수
```

**소수점과 지수부 부호 중 어느 하나라도 빠지면 스칼라가 문자열로 해석된다.**
`run_policy.yaml` 의 다른 값들(`5.0e-5`, `1.0e-7`, `4.5e-4`)이 전부 통과한 것은
소수점과 음수 지수를 둘 다 갖고 있었기 때문이다 — **우연히 안전했다.**

## 왜 위험한가

문자열이 되면 두 갈래로 갈린다:

| 경로 | 결과 |
|---|---|
| 산술에 바로 쓰임 | `TypeError` — 즉시 발견된다. **운이 좋은 경우** |
| 비교·직렬화·로깅만 됨 | `'6.3e6' > 1e6` 은 파이썬 3 에서 `TypeError` 지만, `if value:` 는 True, `str(value)` 는 정상, JSON 왕복도 정상 → **리포트에 그대로 실린다** |

두 번째가 진짜 위험하다. 예산 게이트가 문자열을 받고도 통과할 경로가 있으면
"10분 예산"이 검사되지 않은 상태로 런이 나간다.

## 처방

### ① 데이터 — 지수부 부호를 명시한다

```diff
- throughput_particle_steps_per_s: 6.3e6
+ throughput_particle_steps_per_s: 6.3e+6
```

### ② 로더 — 숫자처럼 보이는 문자열을 거부한다

`simbot/policy.py::_find_numeric_strings` 가 정책 트리 전체를 훑고, 하나라도 있으면
`load_policy` 가 **예외를 던진다.** 다음에 누가 `1.5e8` 을 추가하면 로드 시점에 걸린다.

### ③ 파서 — 숫자를 먼저 시도하고, 아닐 때만 YAML 에 넘긴다

`simbot/session.py::_parse_scalar`:

```python
try:    return int(s)
except ValueError: pass
try:    return float(s)      # ← 파이썬 float() 은 5e-3 을 읽는다
except ValueError: pass
return yaml.safe_load(s)     # bool·str·list 만 여기로
```

`float()` 은 YAML 보다 관대하므로 순서를 이렇게 두면 함정이 사라진다.

## 재발 방지

| 장치 | 위치 |
|---|---|
| 정책 파일 전체 감시 | `tests/test_s4_nondim.py::test_real_policy_has_no_string_numbers` |
| 로더가 거부하는지 | `tests/test_s4_nondim.py::test_policy_rejects_number_parsed_as_string` |
| `set` 파싱 (7가지 표기) | `tests/test_cli_session.py::test_scientific_notation_never_becomes_a_string` |
| 하류에서 실제로 쓰이는지 | `tests/test_cli_session.py::test_set_scientific_notation_is_usable_downstream` |
| 파일 안 주석 경고 | `config/run_policy.yaml` `hardware:` 블록 |

## 적용 범위

**YAML 로 숫자를 읽는 모든 곳.** 이 프로젝트에서는 `config/run_policy.yaml`,
`examples/*/spec.yaml`, `examples/*/prediction.yaml`, `sessions/*/session.yaml`.

`spec.yaml` 은 `simbot.spec.dump_yaml` 이 파이썬 `float` 을 써서 쓰기 때문에 안전하다
(`5.0e-06` 처럼 소수점을 항상 포함해 직렬화된다) — **손으로 쓴 YAML 이 위험하다.**

## 참고

- [[dt-gate-should-be-displacement-based]] — `dt_star` 를 다루는 게이트
- `config/run_policy.yaml` §1 hardware
