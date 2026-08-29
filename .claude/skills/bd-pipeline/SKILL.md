---
name: bd-pipeline
description: 자료(손그림·텍스트·데이터) 한 장에서 검증된 Brownian Dynamics 결론까지 S1→S8 파이프라인을 수행한다. 사용자가 그림/사진을 주며 "이거 시뮬레이션 해줘", "이 계를 돌려봐", "MSD 를 재줘" 라고 할 때, 또는 BD/콜로이드/광집게/ABP 계를 설계·실행·검증할 때 쓴다. 이미 돌아간 런을 진단하려면 bd-diagnose, knowledge 검색·추가는 bd-knowledge 를 쓴다.
---

# bd-pipeline — S1 → S8

> **너는 판단하고, `simbot` 이 계산한다.** 숫자를 머리로 만들지 않는다.
> 모든 수치는 함수 호출 결과여야 한다 (CLAUDE.md §결정론 코어).

## 먼저 읽을 것 (매번)

1. [`CLAUDE.md`](../../../CLAUDE.md) — 제1원칙, 판정 규약, 실행 환경
2. [`knowledge/wiki/systems/_index.md`](../../../knowledge/wiki/systems/_index.md) — 이 계의 카드가 있는가
3. 카드가 있으면 그 카드 — **무차원화와 게이트는 카드가 소유한다**

읽지 말 것: `docs/history/` 의 예전 master_plan 전문. 필요한 절만 링크로 간다.
현행 설계는 [`docs/01-architecture.md`](../../../docs/01-architecture.md).

### 도메인 스킬 — 단계마다 **코드를 쓰기 전에** 읽는다

이 세 스킬은 파이프라인과 경쟁하지 않는다. 파이프라인이 *언제* 무엇을 할지 정하고,
이들이 *어떻게* 해야 조용히 틀리지 않는지를 정한다. 2026-08-28 병합으로 들어왔다.

| 단계 | 읽을 스킬 | 왜 |
|---|---|---|
| **S1** 판독 | **`bd-intake`** | 전사 우선 · 모호성 명시 · 결측은 `null`. 지어내기 억제 8규칙 |
| **S3·S4** 명세·무차원화 | **`bd-physics`** | 스케일 표 작성법 · 기준 척도 선택 · 스케일 분리 검사 · 역변환 |
| **S5** 실행 | **`bd-hoomd`** | HOOMD 함정 20개. **여럿이 에러 없이 조용히 틀린 결과를 낸다** (+1856% 오차 / `angle.Harmonic` 의 힘만 96% 틀리는데 에너지는 정확) |

⚠️ **S5 에서 HOOMD 코드를 한 줄이라도 새로 쓰면 `bd-hoomd` 를 먼저 읽는다.**
이것은 권고가 아니다 — 함정 목록은 전부 실제로 당한 것이고, 에너지 검증으로도
안 잡히는 것이 그 안에 있다.

## 인터프리터

```
/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```

`conda activate` 는 non-interactive shell 에서 불안정하다. **절대경로를 쓴다.**

---

## 단계와 게이트

각 단계는 **게이트를 통과하지 못하면 다음으로 가지 않는다.** 되돌아가거나 사용자에게 보고한다.

| 단계 | 하는 일 | 도구 | 게이트 |
|---|---|---|---|
| **S1** 판독 | 자료 → 관찰/추론/가정 3단 분리 + 모호성 후보 | 너 (Opus) | 차원·경계·구동 확정, `question` 이 반증 가능한가 |
| **S2** 예측 | 시뮬레이션 **전에** 답을 적고 봉인 | `simbot.estimators` | 정량 예측 ≥1, 각각 `tolerance` + `basis` |
| **S3** 명세 | provenance 붙은 `spec.yaml` | `simbot.spec` | 빈 필드 없음, 타당성 검사 통과 |
| **S4** 무차원화 | 카드 척도 + `dt` | `simbot.nondim` | 왕복오차 `< 1e-12`, `dt` 제약 기록 |
| **S5** 실행 | 시드 앙상블 배치 | `simbot.run` | 완주, 가드 무위반 |
| **S6** 그림 | 필수 진단 세트 | `simbot.viz` | 캡션·이중축 전부 존재 |
| **S7** 검증 | 봉인 검증 → 대조표 → 판정 **제안** | `simbot.validate` | 봉인 무결, FAIL 에 원인 분류 |
| **S8** 결론 | 질문에 답 + knowledge 커밋 | 너 (Opus) | knowledge 최소 1항목 |

**S3–S8 은 한 명령이다:**

```bash
/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python cli.py run <spec.yaml> --prediction <prediction.yaml>
```

너가 손으로 하는 것은 **S1·S2 작성**과 **S7 해석·S8 서술**이다.

---

## 절차

### 1. S1 — 자료 판독 → `01_intake.md`

프로토콜: [`references/s1_intake_drawing.md`](references/s1_intake_drawing.md) ← **반드시 읽는다**

핵심만: 손그림의 **절대 크기를 신뢰하지 않는다.** 믿는 것은 토폴로지·비율·개수·대칭성·
**명시된 숫자**. 모호한 요소는 임의 해석 대신 **후보 2~3개를 명시**하고 각 후보의
결과 차이를 예측한다.

### 2. 카드 확인

```bash
ls knowledge/wiki/systems/
```

- 카드가 있다 → §3 기준단위, §7 게이트 표를 그대로 따른다
- **카드가 없다 → 즉흥 무차원화 금지.** `_TEMPLATE.md` 로 `status: draft` 카드를 먼저 만들고
  `simbot/nondim.py::CARD_SCALE_RULES` 에 등록한다

### 3. S2 — 예측 → `prediction.yaml`

프로토콜: [`references/s2_prediction.md`](references/s2_prediction.md)

수치는 `simbot.estimators` 를 **호출해서** 얻는다. 예시:
[`examples/trap-2d-5um/prediction.yaml`](../../../examples/trap-2d-5um/prediction.yaml)

### 4. S3 — 명세 → `spec.yaml`

프로토콜: [`references/s3_s5_execute.md`](references/s3_s5_execute.md)

예시를 복사해서 고치는 것이 가장 빠르다:
[`examples/trap-2d-5um/spec.yaml`](../../../examples/trap-2d-5um/spec.yaml)

**모든 값에 `provenance` 와 `basis`.** 검사:

```bash
<PY> -c "
from simbot.spec import SystemSpec, validate
r = validate(SystemSpec.load('<spec.yaml>'))
print(r.table()); print(r.problems or '규약 위반 없음')"
```

### 5. 비용 확인 → 실행

```bash
<PY> cli.py run <spec.yaml> --prediction <prediction.yaml>
```

CLI 가 예산 초과·게이트 위반이면 **실행 전에 멈춘다.** `--force` 는 사용자가
명시적으로 요구할 때만.

파라미터를 흔들어 보려면 **실행 없이** 비용만 본다:

```bash
<PY> -m simbot.session new <spec.yaml>
<PY> -m simbot.session set numerics.dt_star=2.5e-3 species.0.n_simulated=4000
```

### 6. S7 해석

프로토콜: [`references/s6_s7_validate.md`](references/s6_s7_validate.md)

CLI 가 대조표와 판정을 **제안**한다. 너의 일은:
- `INCONCLUSIVE` 가 **예견된 것인가** 아니면 설계 실수인가
- `FAIL` 의 원인 분류: `numerical` / `modeling` / `interpretation` / `analysis`
- `PASS ⚑` (유의한 편차인데 tolerance 안) 가 나오면 **예측이 허술했는지** 확인

### 7. S8 결론 → `08_conclusion.md` + knowledge

프로토콜: [`references/s8_knowledge.md`](references/s8_knowledge.md)

`REPORT.md` 는 CLI 가 만든다. **너가 쓰는 것은 질문에 대한 답·원인 가설·다음 실험이다.**

---

## 절대 하지 않는 것

| ❌ | 왜 |
|---|---|
| 숫자를 머리로 계산 | 물리가 틀렸는지 너가 틀렸는지 영영 알 수 없다 |
| `confirmed_by` 채우기 | 사람만 쓴다. 아무도 안 본 합격 도장이 찍힌다 |
| 오차막대 없는 수를 결론에 쓰기 | 시드 1개 프로덕션 런은 금지 |
| "검증됐다"·"평형에 도달했다" | 임계값이 있어야 할 수 있는 말이다 |
| 문헌값을 기억으로 인용 | `knowledge/source/papers/` 의 증류를 인용한다 |
| `reproduced: no` 를 근거로 쓰기 | 사실 기록이지 근거가 아니다. `[출처, 미재현]` 으로 표기 |
| 카드 없는 쌍에 즉흥 무차원화 | `nondim.py` 가 예외를 던진다 — 우회하지 말 것 |
| API 키 요구 | 판독하는 LLM 이 이미 이 세션 안에 있다 |

## 질문 예산 — 라운드당 3개

**대화가 스무고개가 되면 도구가 실패한 것이다.** 다음 셋 중 하나일 때만 묻는다:

1. 결론이 뒤집히는 값 (감도 `|S| > 1` 예상)
2. 레짐 경계 근처 (해석이 갈리는 지점)
3. 자릿수조차 모름 (knowledge·문헌에 근거 없음)

그 외에는 `provenance: assumed` 로 채우고 **제안표를 보여주고** 진행한다.
제안에는 항상 3요소: **값 · 근거 · 신뢰도**. 그리고 **바꿀 위치**를 함께 알려준다.

## 실패했을 때

`bd-diagnose` 스킬로 넘어간다. 그리고 **실패를 기록한다** —
지우면 "이 봇이 가능한가"에 대한 증거의 절반이 사라진다.
