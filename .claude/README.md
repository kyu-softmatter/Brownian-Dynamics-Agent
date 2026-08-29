# `.claude/` — L1 에이전트 층

`simbot/` 이 계산하고, 여기가 판단한다. 이 층에는 **숫자를 만드는 코드가 없다.**

## 구성

```
.claude/
├── settings.json                  권한 — 인터프리터 허용 + 봉인 문서 편집 금지
├── skills/
│   ├── bd-pipeline/               [메인] S1→S8 오케스트레이터
│   │   ├── SKILL.md               단계·게이트·금지사항 체크리스트
│   │   └── references/
│   │       ├── s1_intake_drawing.md   ★ 손그림 판독 — 스킬 층의 유일한 고유 내용
│   │       ├── s2_prediction.md       예측·봉인·검정력
│   │       ├── s3_s5_execute.md       명세·무차원화·실행 (대부분 cli.py 호출)
│   │       ├── s6_s7_validate.md      그림·판정
│   │       └── s8_knowledge.md        결론·지식 커밋
│   ├── bd-diagnose/SKILL.md       터진 런 진단 (배제 순서)
│   └── bd-knowledge/SKILL.md      knowledge 검색·추가·정리
└── agents/                        9개. model: frontmatter 로 티어링
```

## 참조문서를 5개로 쪼갠 이유 (master_plan Q6 결정, 2026-07-28)

설계(§12.3)는 단계별 8개였다. 그러나 **결정론 코어가 완성된 뒤 S3·S4·S5 는
`cli.py run` 한 줄**이 되었고, 각각에 별도 문서를 두면 "이 함수를 호출한다"만 적힌
얇은 파일 3개가 생긴다.

대신 **내용이 있는 곳에 문서를 둔다:**

| 문서 | 왜 독립인가 |
|---|---|
| `s1_intake_drawing.md` | **코드로 표현할 수 없는 유일한 단계.** 가장 비싼 오류 지점 |
| `s2_prediction.md` | 봉인·tolerance·검정력 규율. 여기서 허술하면 검증이 무력화된다 |
| `s3_s5_execute.md` | 셋 다 `cli.py` 호출 + 게이트 읽기. 합쳐야 흐름이 보인다 |
| `s6_s7_validate.md` | 그림과 판정은 같은 판단(무엇이 이상한가)에 쓰인다 |
| `s8_knowledge.md` | 결론 서술 + knowledge 계약 |

## 모델 티어링 (master_plan §12)

> 원칙: **추출은 저가, 해석은 고가. 그리고 계산은 LLM 이 아니라 코드.**

| 에이전트 | model | 담당 |
|---|---|---|
| `bd-intake-extract` | haiku | S1 텍스트·숫자·파일 추출 |
| **`bd-intake-interpret`** | **opus** | S1 물리 해석 — 여기서 틀리면 전부 틀린다 |
| **`bd-predict`** | **opus** | S2 봉인되는 주장 |
| `bd-spec` | sonnet | S3 규칙 적용 |
| **`bd-validate`** | **opus** | S7 판정·원인 추론 |
| **`bd-conclude`** | **opus** | S8 최종 주장 |
| `bd-lit-distill` | sonnet | 문헌 증류 (식 변환은 Opus 검토) |
| `bd-lit-scan` | haiku | 서지·INDEX 대량 처리 |
| **`bd-diagnose`** | **opus** | 실패 진단 |

### 안전장치 — 기계적으로 검사 가능하다

> **`provenance` 가 `inference` 또는 `assumed` 인 필드는 Opus 만 쓸 수 있다.**

`observation` · `derived` · `rule` · `from_knowledge` 는 저가 모델이 채워도 된다.
`bd-spec`(sonnet)·`bd-intake-extract`(haiku) 의 지시문에 이 경계가 명시돼 있고,
`simbot.spec.Quantity.problems()` 가 `written_by` 를 검사한다.

**비용 절감이 목적이 아니라 속도와 품질의 배분이 목적이다.**

## `settings.json` 의 두 결정

**① 인터프리터 절대경로를 허용한다.** `conda activate` 는 non-interactive shell 에서
불안정하므로 **거부 목록에 넣었다** — 실수로 쓰는 것을 막는다.

**② 봉인 문서의 `Edit` 을 거부한다.**

```json
"deny": ["Edit(./runs/**/02_prediction.md)", "Edit(./runs/**/01_intake.md)", ...]
```

예측은 `cli.py` 가 파이썬으로 **쓰고**(생성), 그 뒤 에이전트가 텍스트 편집으로
**고치는 것**을 막는다. 사후합리화의 가장 쉬운 경로를 구조적으로 닫는다.
봉인 검증(`SEALED.sha256`)이 사후에 잡지만, **애초에 못 하게 하는 것이 낫다.**

## 이 층이 하지 않는 것

- 숫자를 만들지 않는다 — 전부 `simbot` 호출 결과
- `confirmed_by` 를 채우지 않는다 — 사람만
- `master_plan.md` 전문을 읽지 않는다 (1300줄) — 필요한 절만 링크로
- 물리를 다시 적지 않는다 — `knowledge/wiki/` 와 카드를 인용한다
