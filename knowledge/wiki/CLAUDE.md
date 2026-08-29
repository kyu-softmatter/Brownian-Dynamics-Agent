---
wiki: bd-agent-wiki
status: pilot
contract_version: 0.1

layers:
  raw: ../raw/
  source: ../source/
  wiki: ./

structure: by_publication  # source/ 는 *발표 여부*로 가른다 — 저작자가 아니라.
                           #   papers/ = 발표된 것 (우리 랩 논문 포함, lab_authored: true)
                           #   lab/    = 미발표 랩 자산 (코드·노트·미발표 파라미터)
                           # 공개 경계가 폴더 단위로 걸리기 때문이다 (master_plan.md §4).
                           # 주제 라우팅은 wiki/ 에서 한다.

source_kinds: [paper, arxiv, lab, repo, book, dataset]

filename_conventions:
  paper: "<year>-<firstAuthor>-<slug>"
  arxiv: "<year>-<firstAuthor>-<slug>"
  lab:   "<year>-<firstAuthor>-<slug>"
  repo:  "<org>-<repo>"
  book:  "<year>-<firstAuthor>-<slug>"

source_frontmatter_required: [type, kind, title, authors, year, source_url, ingested_at, ingested_by]
source_frontmatter_optional: [doi, journal, arxiv_id, raw_file, si_available, si_file, tags,
                              cites, related_systems, summary, access, lab_authored]

# lab_authored: true 이거나 kind: lab 이면 추가로 요구되는 필드
lab_frontmatter_required: [engine, reproduced]
lab_frontmatter_optional: [code_location, reproduction_finding, si_file, parameters_extracted]

wiki_frontmatter_required: [type, author, drafted]
wiki_types: [concept, technique, system, benchmark, finding, question]
author_values: [agent, human, hybrid]

promotion:
  finding_to_concept: human_gated

ingest:
  human_gated: true
  agent_writes_raw: never_without_approval
  agent_writes_source: with_approval
  agent_writes_wiki: with_attribution_either_author

precedence:
  L0: ../../.claude/rules/
  L1: ./            # wiki — 해석·합성·결정
  L2: ../source/    # 원문이 실제로 무엇을 말하는가
  L3: ../raw/       # 원본 자체
  rule: "낮은 L이 '무엇이 사실인가'를 이긴다. 높은 L이 '그래서 무엇을 할 것인가'를 이긴다"
---

# bd-agent-wiki — 콜로이드 시뮬레이션 지식 계약

BD/MC 시뮬레이션 에이전트의 지식 층. **위 YAML frontmatter가 기계가 읽는 계약**이고,
아래 산문은 사람이 읽는 근거다. 코드는 경로를 하드코딩하지 말고 frontmatter를 파싱한다.

## 이 위키의 목적

일반 문헌 저장소가 아니다. **두 가지 일**을 한다.

1. **검증 오라클 공급** — `benchmarks/`가 `pytest`로 실행되는 회귀 테스트의 근거가 된다.
   우리 도메인에는 채점자가 없으므로 문헌이 그 역할을 대신한다 (`master_plan.md` §6).
2. **파라미터 사전 공급** — "이 계에서 `dt`를 얼마로 했는가"에 **출처가 붙은 답**을 준다.

읽기용 컨텍스트로만 쓰이는 페이지는 이 위키에 있을 이유가 약하다.

## 계층

```
../raw/        # GITIGNORED. PDF·코드 원본 로컬 캐시
../source/     # 원본 1개당 증류 .md 1개. 출처는 frontmatter에
  ├── papers/  #   발표된 문헌 — 우리 랩 논문 포함 (lab_authored: true).  in-git, 공개 나감
  └── lab/     #   미발표 랩 자산 — 코드·노트·미발표 파라미터.  확보 시 gitignore, 공개 안 나감
./             # 합성 층
  concepts/    # WHAT-IS  — 무차원수, 상거동, 퍼텐셜 종류, MIPS 같은 개념
  techniques/  # HOW-TO   — 평형화 판정, 오차막대, 초기배치, 렌더링
  systems/     # ★ 계-동역학 카드 — 무차원화 규약이 여기서 정해진다
  benchmarks/  # 검증 오라클. benchmarks.yaml + 항목별 근거 페이지
  findings/    # Q→A + 인용.  dead-end-<slug>.md 포함
  questions/   # 아직 답 없는 것. 삭제하지 않고 status로 닫는다
```

## `systems/` — 계-동역학 카드 ★

> **무차원화 방법과 주요 파라미터는 계 하나만으로 정해지지 않는다.**
> 같은 콜로이드라도 평형 구조를 보느냐 수송을 보느냐에 따라 기준 시간이 달라지고,
> 따라서 `Δt` 게이트도 관측량도 벤치마크도 달라진다.

실측 근거 — 같은 랩 안에서 기준 단위가 셋으로 갈린다:

| 계 · 목적 동역학 | 기준 길이 | 기준 시간 | `kT` |
|---|---|---|---|
| ABP × 제어 | **런 길이 `ℓ`** | **`τ_r = 1/D_r`** | **유도량** |
| 브러시 콜로이드 × 비평형 접촉 | `σ` | `τ_D = σ²/D` | 입력값 |
| 수동 tracer × 수송 | `σ` | `τ_D` | 입력값 |

따라서 지식의 **1차 조직 단위는 `(계, 목적 동역학)` 쌍**이다. 파일명 `<계>--<동역학>.md`.

```yaml
type: system
system: <계 슬러그>          # passive-sphere · abp · attractive-colloid · brush-colloid ...
dynamics: <동역학 슬러그>     # equilibrium-structure · transport · coarsening · dense-collective ...
status: draft | usable | validated
```

**카드가 소유하는 것** (템플릿: [`systems/_TEMPLATE.md`](systems/_TEMPLATE.md))

| § | 내용 |
|---|---|
| 3 | **기준 단위** — 길이/에너지/시간 선택과 그 이유 |
| 4 | **무차원수 원장** — 이 쌍에서 의미 있는 것만 |
| 5 | 주요 파라미터 실측값 |
| 6 | 관측량 |
| 7 | **적용 게이트 — 켜고 끄기** |
| 8 | 벤치마크 |

### 게이트는 쌍마다 다르다 — 이게 카드의 핵심 효용

| 게이트 | 수동 구형 × 평형 구조 | ABP × 조밀 집단 |
|---|---|---|
| 평형화 (`pymbar`) | ✅ 유효 | ❌ **무의미** — 능동계는 열평형에 안 간다 |
| `D_msd = kT/γ` | ✅ 성립 | ⚠️ **성립 안 함** — `D_eff = D_t + U₀²τ_r/2` |
| 이류 변위 `u₀Δt/σ` | 해당 없음 | ✅ 필수 |

**카드가 없는 쌍을 만나면 즉흥으로 무차원화하지 않는다.** `_TEMPLATE.md`로 새 카드를 만들고
`status: draft`로 시작한다. 즉흥 무차원화는 `master_plan.md` §2-a가 지적한
"파라미터 선택이 암묵지다"의 재발이다.

색인과 매트릭스: [`systems/_index.md`](systems/_index.md)

`raw/`는 gitignore한다 — 저작권 있는 PDF를 저장소에 넣지 않기 위해서다. `source/`의 증류가
in-git 정본이고, 원본이 필요하면 frontmatter의 `source_url`로 각자 받는다.

### 왜 저작자가 아니라 발표 여부로 가르는가

공개 경계가 **폴더 단위**로 걸리기 때문이다 (`master_plan.md` §4). 우리 랩 논문을 미발표
자산과 같은 폴더에 두면, 저장소를 공개할 때 **이미 DOI가 붙어 공개된 논문의 증류까지 함께
빠진다.** 보호해야 할 것은 저작자가 아니라 **아직 발표되지 않았다는 사실**이다.

우리 랩 논문은 `papers/`에 두고 `lab_authored: true`로 표시한다 — 그러면 폴더는 공개 가능한
상태로 유지되면서도 "이건 우리 것"을 질의할 수 있다.

## 연구실 산출물은 다르게 다룬다 (`lab_authored: true` 또는 `kind: lab`)

> **논문은 "발표된 결과"고, 연구실 코드·노트는 "실제로 돌아간 파라미터 세트"다.**

논문에는 `dt`나 평형화 스텝 수가 안 적혀 있는 경우가 많다. 그래서 파라미터 사전으로서는
연구실 산출물이 훨씬 값지다. 대신 **검증되지 않은 관행도 함께 딸려온다.**

```yaml
lab_authored: true                    # 우리 랩 저작 — 발표 여부와 무관하게 표시
engine: hoomd | lammps | 자체구현 | 해석해 | 실험 | 없음
reproduced: yes | no | partial        # 우리가 직접 재현해봤는가
reproduction_finding: findings/<slug>.md   # 재현 시도 결과
code_location: <경로 또는 URL>
parameters_extracted: yes | no        # 파라미터를 구조화해서 뽑아냈는가
```

**`reproduced` 추적은 발표 여부와 무관하게 적용한다.** 발표된 랩 논문의 파라미터도
우리가 재현하기 전까지는 근거가 아니다 — 논문에 실렸다는 것이 우리 코드에서 그 값이
동작한다는 뜻은 아니기 때문이다.

### 규율 — 이걸 어기면 위키가 소문 저장소가 된다

**`reproduced: no`인 파라미터를 문헌 근거처럼 인용하지 않는다.**
재현 전까지는 "이렇게 했었다"는 **사실 기록**이지 "이게 맞다"는 **근거**가 아니다.

리포트와 `grounding.md`에서 두 가지를 다르게 표기한다.

| 표기 | 의미 |
|---|---|
| `[출처]` | 검증된 근거 — 문헌 벤치마크 또는 `reproduced: yes` |
| `[출처, 미재현]` | 사실 기록 — 참고는 하되 검증 주장에 쓰지 않음 |

## 저작 대칭성

사람도 에이전트도 위키를 쓸 수 있다. 정직성 장치는 **frontmatter 표기**다.

```yaml
type: concept | technique | system | benchmark | finding | question
author: agent | human | hybrid
drafted: YYYY-MM-DD
confirmed_by: human        # 선택. 사람이 검토한 뒤
cites: [경로들]
```

- **`author: agent` 비율 자체가 자기개선 지표다** (`master_plan.md` §10).
- **승격(finding → concept)은 항상 사람 승인.** 에이전트가 스스로 개념을 인플레이션시키는 것을 막는다.
- 저자와 무관하게 품질 기준은 같다: **근거가 `source/`나 확립된 문헌으로 추적되고, 인용이 구체적
  파일·URL이며, 주장이 확인 가능할 것.** "학습 데이터에서 봤다"는 근거가 아니다.

## 페이지가 생기는 계기

| 계기 | 결과 |
|---|---|
| 질문이 생겼는데 위키에 답이 없음 | `questions/<날짜>-<slug>.md` (`status: open`) |
| 질문에 답함 | `findings/<slug>.md` + 원래 질문을 `status: answered`로 |
| **접근이 실패함** | `findings/dead-end-<slug>.md` — **왜 안 됐는지**를 적는다 |
| 문헌에서 검증값을 추출 | `benchmarks/<id>.md` + `benchmarks.yaml` 항목 |
| 새 계를 다루기 시작 | `systems/<slug>.md` |
| finding이 반복 인용됨 | concept으로 승격 (**사람 승인**) |

### `dead-end` 최소 형식

```yaml
type: finding
subtype: dead-end
author: agent
drafted: YYYY-MM-DD
system: <어떤 계에서>
what_was_tried: <무엇을>
why_it_failed: <왜 — 증상이 아니라 원인>
evidence: <측정값·로그 경로>
what_to_try_instead: <다음 후보>
```

`why_it_failed`가 "발산했다"면 그건 증상이다. **"WCA의 `r⁻¹³` 코어 때문에 오버댐프에서
`F·dt/γ`가 폭발했다"**가 원인이다.

## 우선순위(precedence)

frontmatter의 `precedence` 참조.

- **L0 `.claude/rules/`** — 범위·게이트·행동을 결정한다
- **L1 `wiki/`** — 해석·합성·결정을 결정한다
- **L2 `source/`** — "그 문헌이 실제로 뭐라고 했는가"를 결정한다
- **L3 `raw/`** — 원본 자체에 대한 최종 근거

## 이 위키가 아닌 것

- 런별 실행 상태가 아니다 — `outputs/<run_id>/`
- 코드가 아니다 — `bdkit/`
- 결과 파일이 아니다 — `outputs/`
- PDF 보관소가 아니다 — `raw/`는 gitignore된 로컬 캐시일 뿐이다
