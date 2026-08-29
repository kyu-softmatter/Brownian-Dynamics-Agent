---
name: bd-diagnose
description: 터진·이상한 Brownian Dynamics 런을 진단한다. 시뮬레이션이 NaN 으로 죽었을 때, 결과가 문헌·해석해와 자릿수로 다를 때, 예측이 FAIL 났을 때, MSD·RDF 모양이 이상할 때 쓴다. knowledge/wiki/findings 의 진단 경로를 순서대로 따라간다 — 실패가 쌓일수록 강해진다. 새 런을 설계·실행하려면 bd-pipeline, 진단 결과를 지식으로 남기려면 bd-knowledge 를 쓴다.
---

# bd-diagnose — 터진 런의 원인 찾기

> **틀린 `FAIL` 은 틀린 `PASS` 만큼 나쁘다.** 있지도 않은 물리 문제를 좇게 만든다.
> 그래서 **분석 코드를 먼저 의심한다.**

## 0. 원인 5분류 — 어느 쪽인지 먼저 정한다

| 분류 | 뜻 | 빈도 (첫 완주 실측) |
|---|---|---|
| **`analysis`** | 측정·통계 코드가 틀렸다 | **1/4** ← 가장 위험 |
| `numerical` | `dt` 과대, 적분기, 수렴 부족 | 1/4 |
| `interpretation` | S1 판독이 틀렸다 (차원·단위·화살표) | 0/4 |
| `modeling` | 포텐셜·근사가 부적절 | 0/4 |
| `environment` | 패키지·파서·플랫폼 | 2/4 |

**물리를 마지막에 의심한다.** 첫 완주의 4건 중 물리 문제는 0건이었다.

---

## 1. 배제 순서 — 이대로 따라간다

### ① 통계량이 요동하는가

```python
from simbot.guards import assert_statistic_fluctuates
assert_statistic_fluctuates(samples, "이름")
```

**요동하지 않는 "측정값"은 산술 항등식이다.** 2026-07-28 에 실제로 겪었다:
변위에서 평균을 뺀 뒤 교차상관을 재면 `cross/auto = −1/(n−1)` 이 **항등적으로** 나오고,
200회 반복의 표준편차가 `6.7e-20` 이었다. **결과가 그럴듯해서 통과할 뻔했다.**

### ② 독립 경로 둘이 일치하는가 (자기일관성)

같은 양을 두 방법으로 재서 비교한다:

| 양 | 경로 A | 경로 B | 일치해야 |
|---|---|---|---|
| `⟨x²⟩` | 스냅샷 분산 | MSD plateau / 2d | `1.0` |
| `D` | MSD 기울기 | `kT/γ` | 단시간에만 |

어긋나면 **분석 코드**의 문제다. 물리가 아니다.

### ③ 표본이 독립인가

**상관된 표본에 KS·χ² 를 걸면 항상 기각한다.** `n` 이 커질수록 더 확실히 기각하므로
"표본이 많아서 신뢰할 수 있다"는 직관이 **정확히 반대로** 작동한다.

전문: [`ks-test-needs-independent-samples`](../../../knowledge/wiki/findings/ks-test-needs-independent-samples.md)

검사: 프레임 간격이 상관시간(`~2τ`)보다 긴가.

```python
frames_per_tau = 1.0 / (frame_interval_steps * dt_star)
step = ceil(2.0 * frames_per_tau)      # 이만큼 띄운 프레임만 쓴다
```

### ④ 단위·차원

| 의심 | 검사 |
|---|---|
| `γ = 6πηa` 에 직경을 넣었다 | 시간척도가 정확히 2배 틀리다 |
| 시간 척도를 `τ_D` 로 썼다 | `roundtrip_errors` 가 `1e-3` 이상 |
| `kT` 와 `ε` 를 둘 다 1로 잡았다 | `T* = 1` 로 고정돼 버렸다 |
| 2D 인데 3D Stokes 를 썼다 | 의도면 명시, 아니면 버그 |

```bash
<PY> -c "
from simbot.spec import SystemSpec
from simbot.nondim import roundtrip_errors
print(roundtrip_errors(SystemSpec.load('<spec.yaml>')))"
```

### ⑤ 수치 — `dt` 와 가드

```bash
<PY> cli.py converge <spec.yaml>
```

`dt/2` 에서 답이 바뀌면 `dt` 가 과대다. 바뀌지 않으면 **그 오차로는 차이를 볼 수 없다는
뜻**이고 `dt` 문제가 아니라는 증명은 아니다.

가드 결과: `05_run_manifest.json` 의 `batch` + 각 `raw/*/manifest.json` 의 `guards`.

### ⑥ 그 다음에 물리를 의심한다

- HI 무시가 깨지는 조건인가 (`concepts/no-hydrodynamics.md`)
- 상수-`γ` Stokes 가 부족한가 (`concepts/stokes-drag-corrections.md`) —
  `a ≳ 5 μm` 면 Basset `4.3 %`, Faxén `2.3배`
- 과감쇠가 유효한가 (`τ_i/τ_process ≪ 1`)

---

## 2. 증상별 첫 의심 — knowledge 조회

```bash
ls knowledge/wiki/findings/
grep -rl "<증상 키워드>" knowledge/wiki/
```

| 증상 | 먼저 볼 것 |
|---|---|
| NaN / 폭발 | `dt` 과대, 초기 겹침. `guards.check_finite` 가 어느 프레임에서 잡았나 |
| 스텝당 변위가 `√3σ_step` 을 넘는다 | **불가능하다** — HOOMD 노이즈는 균일분포다. 측정 코드를 의심 |
| 분포 꼬리가 Gaussian 과 다르다 | **정상이다.** `findings/hoomd-brownian-scheme-and-noise.md` |
| 첨도가 정확히 3.000 | **오히려 이상하다.** `3 − 1.2 dt*` 가 나와야 한다 |
| `⟨x²⟩` 가 예측보다 `dt*/2` 만큼 높다 | **정상이다.** EM 편향 |
| KS p = 0.0000 | 상관 표본. ③ 으로 |
| `dt` 게이트가 통과했는데 이상하다 | 변위 게이트가 트랩 계에서 1086배 느슨하다 |
| 비용 추정이 이상하다 / TypeError | YAML 1.1 지수 표기 (`5e-3` 이 문자열) |
| 카드 없는 계를 돌리려다 KeyError | **정상 동작이다.** 카드를 먼저 만든다 |

## 3. 진단 도구

```bash
# 봉인 상태
<PY> -c "
from simbot.io import RunDir, verify_seal
print(verify_seal(RunDir('runs/<id>')).summary())"

# 게이트 재검사
<PY> -c "
from simbot.spec import SystemSpec, validate
r = validate(SystemSpec.load('runs/<id>/03_spec.yaml'))
print(r.table()); print(r.problems)"

# 궤적을 직접 본다
<PY> -c "
from simbot.analysis.trap import load_run
d = load_run('runs/<id>/raw/<label>')
print({k: (v.shape, v.dtype) for k, v in d.items()})
print('var range', d['indep_var'].min(), d['indep_var'].max())"

# 그림만 다시 만든다 (런 재실행 없이)
<PY> cli.py resume runs/<id>
```

## 4. 진단이 끝나면 — 기록한다

**실패도 산출물이다.** 지우면 "이 봇이 가능한가"에 대한 증거의 절반이 사라진다.

```
knowledge/wiki/findings/<slug>.md              원인을 찾았다
knowledge/wiki/findings/dead-end-<slug>.md     막힌 길이다
```

`dead-end` 는 `why_it_failed` 에 **원인**을 적는다. "발산했다"는 증상이고,
**"WCA 의 `r⁻¹³` 코어 때문에 오버댐프에서 `F·dt/γ` 가 폭발했다"** 가 원인이다.

**재발 방지에 "없음"을 적는 것이 반복되면 가드를 만들 신호다.**

## 5. 하지 않는 것

- 원인을 못 찾았는데 "고쳤다"고 하지 않는다. `INCONCLUSIVE` 를 사실로 고정한다
- 증상이 사라졌다고 원인이 사라진 게 아니다 — 무엇이 왜 고쳐졌는지 적을 수 없으면
  고친 것이 아니다
- 게이트를 끄는 것으로 문제를 "해결"하지 않는다. 끌 때는 **이유가 카드에 있어야 한다**
