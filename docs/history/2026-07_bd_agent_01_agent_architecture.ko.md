# 01 · 에이전트 아키텍처

> 이 문서는 **에이전트가 무엇으로 이루어져 있고, 어떤 순서로 무엇을 하며, 상태를 어디에 저장하는가**를
> 정의한다. 개별 스테이지의 물리·수치 내용은 `02`~`09`에서 다룬다.

---

## 0. 먼저: 에이전트 구축 일반 원칙

에이전트를 처음 만들 때 대부분의 실패는 물리가 아니라 **구조**에서 온다.
아래 10가지는 도메인과 무관하게 적용되는 원칙이고, 이 프로젝트의 모든 설계 결정이 여기서 파생된다.

| # | 원칙 | 이 프로젝트에서의 구현 |
|---|---|---|
| 1 | **상태는 파일에 둔다. 대화에 두지 않는다.** 대화 컨텍스트는 사라지고 검사할 수 없다 | `run_state.yaml` |
| 2 | **LLM을 수치 루프 안에 넣지 않는다.** 비결정론적이고 테스트 불가하며 조용히 틀린다 | `bdkit/`은 LLM 의존 0 |
| 3 | **경계마다 구조화 출력을 강제한다.** 자유 텍스트를 파싱하는 순간 버그가 시작된다 | 모든 LLM 출력은 JSON Schema 검증 통과 |
| 4 | **모든 결정을 저널에 남긴다. 누가 내렸는지(actor)까지.** 에이전트 디버깅의 유일한 수단 | `decision_journal.jsonl` |
| 5 | **예산을 먼저 정한다.** 무한 루프는 에이전트의 기본 실패 모드다 | `run_state.budget` (D18) |
| 6 | **작은 티어부터 실행한다.** 첫 실행이 비싼 실행이면 안 된다 | smoke → pilot → production |
| 7 | **판정 기준은 숫자여야 한다.** "그럴듯해 보인다"는 게이트가 아니다 | `\|D_msd/(kT/γ) − 1\| < 0.02` |
| 8 | **모르는 것은 모른다고 하게 만든다.** 과학 작업에서 LLM의 최대 실패는 그럴듯한 파라미터를 조용히 지어내는 것 | `unknowns[]` + `confidence` + `assumed: true` |
| 9 | **되돌아갈 수 있게 만든다.** 체크포인트, 멱등한 스테이지, 재개 | 스테이지별 아티팩트 + `resume` |
| 10 | **사람 게이트는 비용이 뛰기 직전에 둔다.** 그 외에는 두지 않는다 | 게이트 2곳 (D4) |

> **핵심 한 줄: LLM은 제안하고, 결정론적 코드가 판정한다.**

---

## 1. 레이어 구조

```
┌─────────────────────────────────────────────────────────┐
│ 사용자      텍스트 · 스케치 · 사진 · 녹음 · 논문 PDF      │
└───────────────────────────┬─────────────────────────────┘
                            ↓
┌─────────────────────────────────────────────────────────┐
│ agent/     LLM 레이어                                    │
│            추출 · 문헌대조 · 실패 triage · 리포트 서술    │
│            → 반드시 구조화 출력(schema)만 반환            │
└───────────────────────────┬─────────────────────────────┘
                            ↓  JSON Schema 검증 (경계)
┌─────────────────────────────────────────────────────────┐
│ bdkit/     결정론 코어  ── LLM 의존 0, pytest 단독 검증   │
│            계산 · 검증 · 실행 · 진단 · 분석 · 시각화      │
└───────────────────────────┬─────────────────────────────┘
                            ↓
             HOOMD-blue 7.1.0 · freud · OVITO/fresnel
```

### `bdkit/` — 결정론 코어

| 모듈 | 책임 | 관련 문서 |
|---|---|---|
| `spec/` | SystemSpec 스키마·검증. **에러 전체 수집** | `02` |
| `units/` | `UnitMap`(SI ↔ reduced), 무차원수 원장. **계-동역학 카드를 읽어 기준 단위를 정한다** | `03` |
| `plan/` | `RunPlan` 생성·검증, 비용 추정 | `05` |
| `build/` | HOOMD 스크립트 생성 + preflight 정적 검사 | `04`, `05` |
| `run/` | `Runner` 인터페이스 · `LocalRunner`(v1) · `SlurmRunner`(v2) | D2 |
| `diagnose/` | 안정성·열역학·자기일관성·평형화·유한크기 지표. **카드 §7이 켜고 끈다** | `05` |
| `repair/` | 실패 → 조치 규칙 테이블 | `06` |
| `analyze/` | freud 관측량 + block-averaged 오차막대 | `07` |
| `viz/` | 입자 렌더 + matplotlib 플롯 | `08` |
| `report/` | 단일 HTML 조립 | `08` |

**불변 조건: `bdkit/` 어디에서도 LLM을 호출하지 않는다.** `grep -r "anthropic\|claude" bdkit/`가 비어야 한다.

### `agent/` — LLM 레이어

스테이지별로 `프롬프트 + 출력 스키마 + 검증기` 세 쌍.
LLM이 개입하는 스테이지는 **S1, S2, S4, S5(제안만), S9(규칙 미적중 시), S12** 여섯 곳뿐이다.

**호출 방식은 `D23` — `DECIDED` (2026-07-27): Anthropic SDK 직접**, `agent/llm.py` 한 곳에서만.
tool-use로 출력 스키마를 **구조적으로 강제**한다 (원칙 3). 대화 왕복이 필요한 S2 ELICIT만
나중에 Claude Code로 뺄 여지를 남긴다. 비교표는 [`master_plan.md`](../../docs/history/2026-08_simulation_auto_master_plan.ko.md) §4.5,
확정 경위는 [`findings/d23-sdk-backend.md`](../../knowledge/wiki/findings/d23-sdk-backend.md).

세 쌍 중 **프롬프트·출력 스키마·검증기는 호출 방식과 무관하게 재사용된다** — 방식이 바뀌어도
`agent/llm.py` 하나만 교체된다. 이 이음새는 계속 유지한다.

---

## 2. 상태기계

**진입점은 사람이다** (`D25`). 바깥 자율 루프·큐·스케줄러는 v1에 없다.

```
사람:  bd-agent new "<자연어 기술>" [--image ...] [--pdf ...]     → outputs/<run_id>/ 생성
사람:  bd-agent resume outputs/<run_id> [--from S5_PLAN]         → 게이트 승인 후 재개
```

```
S1 INTAKE → S2 ELICIT →🚦게이트1→ S2.5 PREREGISTER(v0)
    → S3 NONDIM → S4 LIT-GROUND → S2.5′ PREREGISTER(v1) → S5 PLAN → S6 PREFLIGHT
    → S7 EXECUTE(smoke → pilot →🚦게이트2→ production)
    → S7.5 EYEBALL → S8 DIAGNOSE ⇄ S9 REPAIR
    → S10 ANALYZE → S11 VISUALIZE → S12 REPORT → DONE
                                                        ↑ S8=REDESIGN이면 S5로 복귀
```

| 상태 | 주체 | 성공 시 → | 실패 시 → |
|---|---|---|---|
| `S1 INTAKE` | LLM | `S2` | `BLOCKED_INPUT` (입력 모순) |
| ↳ | | **`(계, 목적 동역학)` 쌍 분류** → 카드 조회. 분류 실패 시 `unknowns[]` | |
| `S2 ELICIT` | LLM + 사람 | `GATE_SPEC` | `S1` |
| `🚦GATE_SPEC` | 사람 | `S2.5` | `S2` (수정 요청) |
| `S2.5 PREREGISTER` | 사람 + LLM | `S3` | `S2` (예상이 스펙과 모순) |
| `S3 NONDIM` | 코드 | `S4` | `S2` (오버댐프 부적합 → 엔진 변경 협의) |
| ↳ | | **카드 기반** — `wiki/systems/<계>--<동역학>.md` §3·§4·§7 적용. 카드 없으면 신규 생성 | |
| `S4 LIT-GROUND` | LLM | `S2.5′` (경고는 통과) | — (차단 없음) |
| `S2.5′ PREREGISTER` | 사람 + LLM | `S5` | — (차단 없음. v0만으로도 진행 가능) |
| `S5 PLAN` | LLM 제안 → 코드 검증 | `S6` | `S5` (재제안, 예산 차감) |
| `S6 PREFLIGHT` | 코드 | `S7` | `S5` |
| `S7 EXECUTE` | Runner | `S7.5` | `S7.5` (크래시도 진단 대상) |
| `🚦GATE_PROD` | 사람 | `S7:production` | `S5` \| `DONE`(pilot에서 종료) |
| `S7.5 EYEBALL` | 코드 → 사람 | `S8` | `S8` (육안 이상도 진단 입력) |
| `S8 DIAGNOSE` | 코드 + LLM triage | `S10` | `S9`(REPAIR) \| `S5`(REDESIGN) |
| `S9 REPAIR` | 규칙 → LLM → 사람 | `S6` (재검사 후 재실행) | `ESCALATED` (예산 소진) |
| `S10 ANALYZE` | 코드 | `S11` | `S9` (`N_eff` 부족 → 런 연장) |
| `S11 VISUALIZE` | 코드 | `S12` | `S11` (렌더러 폴백) |
| `S12 REPORT` | LLM + 코드 | `DONE` | — |

**S2.5·S7.5 추가 근거 (`D28`, 2026-07-27):** 사전등록은 사후해석을 막는 유일한 장치라 부속물이
아니라 스테이지여야 하고, 육안 검사는 *분석 앞에* 있어야 한다 — 결정화·클러스터·겹침처럼 **눈으로
1초, 숫자로는 어려운** 실패가 실재하기 때문이다. 상세는 [`master_plan.md`](../../docs/history/2026-08_simulation_auto_master_plan.ko.md) §5.
`S2.5′`는 새 상태가 아니라 **S2.5의 두 번째 통과**다 (정성 v0 → 정량 v1).

**종료 상태:** `DONE` · `ESCALATED`(사람 개입 필요) · `BLOCKED_INPUT` · `ABORTED`(사용자 중단)

**멱등성 규칙:** 각 스테이지는 입력 아티팩트가 같으면 같은 출력을 낸다(LLM 스테이지는 온도 0 + 시드 고정).
따라서 어느 스테이지에서든 재개할 수 있다.

---

## 3. `run_state.yaml`

런 하나의 단일 진실 공급원. **모든 스테이지 전이 후 원자적으로 갱신**한다(임시 파일 → `os.replace`).

```yaml
schema_version: 1
run_id: "2026-07-27T11-17-03Z__silica-depletion-gel"
state: S8_DIAGNOSE
tier: pilot                       # smoke | pilot | production
created_at: "2026-07-27T11:17:03Z"
updated_at: "2026-07-27T11:42:18Z"

stage_history:
  - {stage: S1_INTAKE,  entered: "...", exited: "...", result: PASS}
  - {stage: S2_ELICIT,  entered: "...", exited: "...", result: PASS}
  - {stage: GATE_SPEC,  entered: "...", exited: "...", result: APPROVED, by: human}
  - {stage: S2_5_PREREGISTER, entered: "...", exited: "...", result: PASS, version: v0}
  - {stage: S3_NONDIM,  entered: "...", exited: "...", result: PASS}

gates:
  spec:       {status: approved, at: "2026-07-27T11:25:00Z", by: human}
  production: {status: pending}

artifacts:                        # 상대 경로. 없으면 아직 생성 안 됨
  spec_draft:   artifacts/spec_draft.yaml
  spec:         artifacts/spec.yaml
  hypothesis_v0: artifacts/hypothesis-v0.yaml   # S2.5  정성. 근거=직관·경험
  hypothesis_v1: artifacts/hypothesis-v1.yaml   # S2.5′ 정량. 근거=무차원수·문헌
  reduced_spec: artifacts/reduced_spec.yaml
  unit_map:     artifacts/unit_map.yaml
  dimensionless: artifacts/dimensionless.yaml
  grounding:    artifacts/grounding.md
  run_plan:     artifacts/run_plan.yaml
  cost_estimate: artifacts/cost_estimate.yaml
  preflight:    artifacts/preflight_report.md
  sim_script:   artifacts/simulate.py
  eyeball:      figures/eyeball/                # S7.5 스냅샷 3장 (저해상도)
  diagnosis:    artifacts/diagnosis.yaml

budget:                           # D18
  max_total_walltime_s: 21600
  max_repair_iterations: 8
  max_disk_gb: 20
  max_llm_calls: 100
  spent_walltime_s: 412
  repair_iterations_used: 1
  disk_used_gb: 0.8
  llm_calls_used: 14

provenance:                       # 재현성의 핵심. D6 참조
  bdkit_version: "0.1.0"
  git_commit: "TBD"               # D6 미결정 시 "no-vcs"
  hoomd_version: "7.1.0"
  hoomd_gpu_enabled: false
  freud_version: "TBD"
  python: "3.12.13"
  platform: "macOS-15-arm64"
  master_seed: 12345
```

---

## 4. `decision_journal.jsonl`

**에이전트를 디버깅하는 유일한 수단.** 한 줄 = 하나의 결정. append-only.

`actor` 필드로 `rule`(결정론 규칙) / `llm`(모델 판단) / `human`(사람)을 구분한다.
이게 있어야 사후에 **"이 결과에 LLM 판단이 몇 번 개입했나"** 를 셀 수 있다.

```jsonl
{"ts":"...","stage":"S1_INTAKE","actor":"llm","action":"extract_field","field":"particle_radius_m","value":5.0e-7,"confidence":0.6,"provenance":"사진 속 스케일바 '1 µm'에서 추정","note":"confidence<0.8 → unknowns[]에도 등재"}
{"ts":"...","stage":"GATE_SPEC","actor":"human","action":"approve","note":"반지름 500nm 확인함"}
{"ts":"...","stage":"S3_NONDIM","actor":"rule","rule_id":"U03_overdamped_validity","observation":{"tau_B_over_tau_D":2.1e-7},"verdict":"PASS","threshold":1e-3}
{"ts":"...","stage":"S8_DIAGNOSE","actor":"rule","rule_id":"G02_max_displacement","observation":{"max_disp_per_step_sigma":0.34,"threshold":0.10},"verdict":"FAIL"}
{"ts":"...","stage":"S9_REPAIR","actor":"rule","rule_id":"R02_reduce_dt","action":{"param":"run_plan.dt","from":1.0e-4,"to":5.0e-5},"iteration":1,"rationale":"규칙 테이블 적중"}
{"ts":"...","stage":"S9_REPAIR","actor":"llm","action":"triage","observation":"규칙 미적중: nlist 재구성 빈도 이상","proposal":{"param":"run_plan.nlist_buffer","from":0.4,"to":0.8},"iteration":3,"rationale":"..."}
{"ts":"...","stage":"S9_REPAIR","actor":"human","action":"escalate_ack","note":"예산 소진, dt 대신 퍼텐셜을 harmonic으로 교체하기로"}
```

**필수 필드:** `ts` `stage` `actor` `action`
**권장 필드:** `rule_id`(actor=rule) · `observation`(측정값+임계값) · `rationale` · `iteration`

---

## 5. 런 디렉터리 레이아웃

```
outputs/<run_id>/
├── run_state.yaml                  ← 단일 진실 공급원
├── decision_journal.jsonl          ← append-only
├── artifacts/
│   ├── spec_draft.yaml  spec.yaml  hypothesis-v0.yaml  hypothesis-v1.yaml
│   ├── reduced_spec.yaml  unit_map.yaml  dimensionless.yaml
│   ├── grounding.md
│   ├── run_plan.yaml  cost_estimate.yaml
│   ├── preflight_report.md
│   ├── simulate.py                 ← 생성된 HOOMD 스크립트 (사람이 읽고 손으로 돌릴 수 있어야 함)
│   ├── diagnosis.yaml
│   └── observables.parquet  observables_summary.csv
├── raw/
│   ├── smoke/       trajectory.gsd  log.h5
│   ├── pilot/       trajectory.gsd  log.h5  checkpoint.gsd
│   └── production/  trajectory.gsd  log.h5  checkpoint.gsd
├── figures/                        ← PNG (리포트에 base64로 임베드)
│   └── eyeball/                    ← S7.5 저해상도 스냅샷 3장. 분석 전 육안 검사용
└── report.html                     ← 최종 산출물, 자체완결 단일 파일
```

**`artifacts/simulate.py`는 반드시 사람이 직접 실행 가능해야 한다.** 에이전트를 못 믿겠을 때
손으로 돌려볼 수 있는 탈출구이자, 에이전트가 실제로 무엇을 했는지 확인하는 가장 빠른 방법이다.

`run_id` 규약: `<UTC ISO8601, 콜론→하이픈>__<slug>` — 사전순 정렬 = 시간순 정렬.

---

## 6. 재개(resume) 규약

```
bd-agent resume outputs/<run_id> [--from S5_PLAN]
```

1. `run_state.yaml`을 읽어 `state`를 복원
2. `--from`이 없으면 마지막 미완료 스테이지부터
3. 이미 존재하는 아티팩트는 재계산하지 않음 (`--force`로 무시)
4. `S7 EXECUTE`는 `checkpoint.gsd`가 있으면 거기서 이어감
5. 재개 사실 자체를 저널에 기록

---

## 7. 실패 처리 원칙

| 실패 유형 | 처리 |
|---|---|
| **검증 실패** (예상된 것) | 규칙 테이블(`06`) → 예산 내 자동 수정 |
| **규칙 미적중** | LLM triage 1회 → 제안을 결정론 검증기에 통과시킨 뒤 적용 |
| **예산 소진** | `ESCALATED`. **증상 · 시도한 것 전부 · 다음 후보**를 정리해 사람에게 제시 |
| **크래시** (예상 못 한 것) | 스택트레이스 + `run_state` + 저널 마지막 20줄을 묶어 보고. 상태는 보존 |
| **입력 모순** | `BLOCKED_INPUT`. 모순 지점을 짚어 반환 (예: "φ=0.7인데 hard sphere — 물리적으로 불가능") |

**절대 하지 않을 것:** 실패를 조용히 삼키고 다음 스테이지로 넘어가기.
검증 없이 통과한 스테이지는 리포트에 `UNVERIFIED` 배지로 표시한다.

---

## 8. 이 아키텍처의 검증 방법

| 대상 | 방법 |
|---|---|
| `bdkit/` 순수 로직 | `pytest tests/` — HOOMD 없이 동작 (선행 프로젝트의 33개 테스트 관례 계승) |
| 물리 정합성 | `pytest tests/test_benchmarks.py` — 짧은 시뮬 실행 → 문헌값 대조 (`09`) |
| 상태기계 | 각 스테이지를 목(mock) 아티팩트로 전이시키는 테스트. 전이표(§2)와 1:1 |
| LLM 추출 정확도 | 자연어 프롬프트 ~10개 → 기대 `SystemSpec` 골든 파일 대조 |
| 저널 완전성 | 런 종료 후 `actor` 집계. LLM 판단 개입 횟수가 예상 범위인지 |
