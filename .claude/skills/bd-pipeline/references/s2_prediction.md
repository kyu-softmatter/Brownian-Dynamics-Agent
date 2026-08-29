# S2 — 예측 프로토콜

> **시뮬레이션 전에 답을 적는다.** 이 단계가 프로젝트의 과학적 정직성을 담보한다.
> 예측을 나중에 쓰면 그것은 예측이 아니라 설명이다.

## 0. 수치는 호출해서 얻는다

```python
from simbot.spec import SystemSpec, derive
from simbot.estimators import harmonic_trap, euler_maruyama_trap_variance_bias
from simbot.analysis.trap import em_uniform_noise_excess_kurtosis
```

**손계산 금지.** `4.14e-21 × 2` 도 하지 않는다 — 코드가 하면 오차가 없고 기록이 남는다.

예측 파일을 코드로 생성하는 것이 가장 안전하다. 예시가 그렇게 만들어졌다:
[`examples/trap-2d-5um/prediction.yaml`](../../../../examples/trap-2d-5um/prediction.yaml)

## 1. 예측 1개의 4요소

```yaml
- quantity: var_x_star            # 측정 이름과 **정확히** 일치해야 한다
  value: 1.0025062656641603       # simbot 함수 출력. 절단하지 않는다
  tolerance: ±1%                  # 이 밖이면 FAIL
  basis: 'Euler-Maruyama 정상분산 1/(1-dt*/2). estimators.euler_maruyama_...'
  discriminates: '적분기 스킴이 EM 인가 exact 인가'
  competing_value: 1.0            # ★ 경쟁 가설 — 검정력 계산에 쓰인다
```

`quantity` 가 측정 이름과 다르면 `validate_run` 이 **"대응하는 측정이 없다"** 로 보고한다.
측정 이름은 `cli.py::measure_trap_batch` 가 정한다.

### `value` 를 절단하지 않는다

`1.00251` 로 적으면 검정력이 `0.5628σ` 로 나오고, 전정밀도 `1.0025063` 이면 `0.5618σ` 다.
**해석식의 값을 그대로 쓴다.**

## 2. `tolerance` — 형식과 함정

| 형식 | 의미 | 예 |
|---|---|---|
| `±X%` | 예측값의 상대 대역 | `±1.5%` |
| `±X` | 절대 대역 | `±0.03` |
| `>X` / `p>X` / `R^2>X` | 단측 하한 | `>0.99`, `p>0.05` |
| `<X` | 단측 상한 | `<1e-12` |

파서: `simbot.validate.parse_tolerance`. **읽지 못하면 예외를 던진다** — 조용히 넘어가면
그 항목이 판정 없이 통과한다.

### ❌ tolerance 를 넓게 잡아 어떤 결과든 PASS

**금지, 리뷰 대상** (master_plan §S2 실패모드). 판정기가 이제 이것을 잡는다:
`PASS` 인데 편차가 `3σ` 를 넘으면 `PASS ⚑` 로 표시하고 문제로 올린다.

★ 실제로 걸렸다: MSD plateau 를 `2d = 4.0` 정확히로 예측했더니 `±2%` 대역 안에서
`3.54σ` 어긋났다. 원인은 **알려진 편향을 예측에 넣지 않은 것** —
`plateau = 2d⟨x*²⟩` 이므로 EM 편향이 곱해진다.
전문: [`wide-tolerance-hides-significant-deviation`](../../../../knowledge/wiki/findings/wide-tolerance-hides-significant-deviation.md)

⇒ **알려진 계통 편향은 예측에 포함시킨다.** "이상적인 값"이 아니라 "이 스킴에서
나와야 하는 값"을 적는다.

## 3. `competing_value` — 검정력을 미리 계산한다

경쟁 가설을 적으면 판정기가 **설계 검정력**을 계산한다:

```
검정력 = |예측 − 경쟁| / SE
```

`< 1σ` 면 그 측정은 두 가설을 구별하지 못하므로 `INCONCLUSIVE` 다.

### ★ 이것을 예측 단계에서 계산한다

시드 4개의 SE 를 미리 알 수 있으면 (`estimators.samples_for_variance_precision`),
**어떤 예측이 판별 불가일지 사전에 안다.** 첫 예시에서 `P8`·`P9` 가 그랬고,
예측 문서에 "INCONCLUSIVE 예상"을 적어뒀다 — 그래서 결과가 나왔을 때
실패가 아니라 예견된 한계로 읽혔다.

**판별 불가를 미리 알면 두 선택지가 있다:**
1. 그대로 진행하고 `INCONCLUSIVE` 를 사실로 고정한다 (결론이 그에 의존하지 않으면)
2. 검정력이 생기는 조건으로 바꾼다 — 적분기 검증은 **일부러 큰 `dt*`** 에서 한다
   (`dt*=2e-2` 에서 `3.8σ`, `5e-3` 에서 `0.56σ`)

> **프로덕션 `dt*` 에서는 적분기를 검증할 수 없다.** `dt*` 를 줄이면 편향도 줄지만
> 검증 가능성도 함께 줄어든다.

## 4. `alternatives` — 예측이 틀릴 방식과 그때의 신호

```yaml
alternatives:
  - 'dt* 가 크면 <x*^2> 가 1/(1-dt*/2) 만큼 높게 나온다 — 이것이 정상이다.'
  - '첨도는 정확히 3 이 아니라 3 - 1.2 dt* 다 (균일 노이즈). 정확히 3.000 이면 이상하다.'
  - 'tau_trap 이 어긋나면 eta 가정 또는 a 해석이 틀렸다는 신호다.
     <x^2> 는 영향받지 않으므로 둘을 분리해서 판별할 수 있다.'
```

**"틀리면 어떻게 보일지"를 미리 적어두면 S7 에서 원인 추론이 훨씬 빠르다.**

## 5. `regimes` — 무차원수와 경계 근접도

```python
from simbot.nondim import groups
groups(spec)     # 계산 가능한 것만 나온다. 불가능한 것은 넣지 않는다
```

각 무차원수가 **레짐 경계에서 얼마나 떨어져 있는지** 적는다.
`k*_σ = 2.4e5` 는 경계(`~1`)에서 5.4 decade — "극단적으로 멀다"까지 적어야
"이 값이 조금 틀려도 결론이 안 바뀐다"를 말할 수 있다.

경계 근처(`Pe = 45`, MIPS `Pe_c ≈ 40–60`)면 **사용자에게 묻는다** — 질문 예산의
정당한 사용이다.

## 6. 봉인

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>
```

CLI 가 `SEALED.sha256` 을 **실행 전에** 쓴다. 표준 `sha256sum` 형식이라
우리 코드 없이도 검증된다:

```bash
shasum -a 256 -c SEALED.sha256
```

**봉인 후 예측을 고치면 S7 이 대조표를 만들지 않는다.** 우회 경로를 찾지 말 것 —
그 경로가 있으면 봉인이 아무것도 보증하지 않는다.

`--prediction` 없이 돌리면 CLI 가 경고한다:
> ⚠️ 예측 파일이 없다 — S7 대조 없이 진행한다. 봉인할 것이 없으므로 사후합리화를
> 막을 장치가 없다.

**탐색 단계에서는 괜찮다. 결론을 낼 런에서는 안 된다.**

## 7. 게이트

- 정량 예측 ≥ 1개, 각각 `tolerance` + `basis` (`Prediction.problems()` 가 검사)
- 그림에 사용자가 그린 예상 곡선이 있으면 그것과의 정합/불일치를 명시
- 알려진 계통 편향이 예측값에 반영됐는가 ← **2026-07-28 에 추가된 항목**
