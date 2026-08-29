# Brownian Dynamics Simulation Bot — Master Plan

> 작성일: 2026-08-03 · 최종 갱신 **2026-08-04** · **버전 v0.5**
> 대상 엔진: [HOOMD-blue v7.1.1](https://hoomd-blue.readthedocs.io/en/v7.1.1/)
> 실행 환경: macOS (Apple Silicon, 10 cores / 16 GB), CPU 전용, 로컬
> v0.1 → v0.2: 멀티모달 인테이크 · 문헌 지식베이스 · 명시적 무차원화 레이어 ·
> 로우데이터 저장/재분석 · 사람 승인 게이트 · 사후분석 학습 루프
> v0.2 → v0.3: 차원 우선을 하드 불변식으로 격상 · Claude Code를 에이전트 런타임으로 ·
> 스코프 확장 필요 발견
> v0.3 → v0.4: **스코프 개방 + 물리 모듈 레지스트리**(§5.6) · **Phase 0 완료**(능력 실측,
> 함정 9건) · **Phase 8-min 완료**(CLAUDE.md + 스킬 2종) · **Phase 1-A 완료**(관측량 4종
> 해석해 일치) · **분리 검사를 모델/적분/기하/통계로 재분류**(§6.4 — `τ_p/dt` 오류 수정) ·
> **dt를 편향에서 역산**하는 정량 규칙 · **KB는 `record.json`부터**(§7.0)
> v0.4 → v0.5: **원칙 8**(검증≠일치) · **원칙 9**(독립 요소 단독 검증) +
> **9.1 조합에 기존 이론 금지** + **9.2 예측 사전 등록** · 원칙 순서 1→9 재배치 ·
> `bdbot/` 16모듈 + CLI + L3 `NondimSpec` 완료 · 함정 13·14(ABP) · KB 37건

---

## 0. 한 줄 정의

**손으로 그린 스케치·메모·문헌으로부터 물리계를 해석하고, 축적된 지식을 근거로 파라미터를 제안·무차원화하여 Brownian dynamics 시뮬레이션을 실행하고, 그 성공/실패 경험까지 다시 지식으로 환류시키는 LLM 에이전트.**

핵심은 **닫힌 루프**입니다. 봇은 쓸수록 똑똑해져야 합니다.

```
       ┌─────────────────── 지식베이스 (KB) ───────────────────┐
       │  문헌 증류 · 무차원수 좌표 · 우리 런의 성공/실패 경험    │
       └──▲──────────────────────────────────────────────┬────┘
          │ (참조: 파라미터 제안 근거)      (환류: 사후분석) │
          │                                              ▼
  스케치/메모/그림 → 물리계 해석 → [사람 승인] → 무차원화 → 시뮬레이션 → 분석
                     (SI, 차원 있음)              (스케일 원장 기반)
```

**차원이 먼저입니다.** 모든 계는 SI 단위의 물리계로 먼저 확정하고, 계에 존재하는
주요 길이·시간·에너지 스케일을 원장에 열거한 뒤, 그중 기준을 골라 무차원화합니다.
무차원 값으로 시작하는 우회 경로는 없습니다 ([원칙 3](#원칙-3--차원이-먼저다-무차원화는-그-다음이다--불변식)).

---

## 0.1 현재 상태 (2026-08-04) — 이 절이 상태의 단일 진실 공급원

| 단계 | 상태 | 실체 |
|---|---|---|
| **0** 범용 환경 + HOOMD 능력 실측 | ✅ | `environment.yml` · [`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md) · 15 API 실동작 |
| **8-min** 지식 캡처 | ✅ | [`CLAUDE.md`](../../CLAUDE.md) 절대 규칙 9개 (7 이 둘 — `7`·`7'` 관례) · skill **3종** (`bd-hoomd` `bd-physics` `bd-intake`) |
| **1-A** `trap-2d-5um` 관통 | ✅ | 관측량 4종 해석해 일치 · 2런 |
| **1-B** `soft-r3-2d-A-sweep` 관통 | ✅ | 8런 · 검증 5종 + 수렴 2종 |
| **1-C** 추상화 → `bdbot/` 패키지 | ✅ | **16모듈** + CLI(`status/intake/interactions/system/nondim/run`) |
| **앞단** 인테이크 스키마 + 검사 | ✅ | 적대적 검사 27/27 |
| **L3** `NondimSpec` (L2↔L4 유일 계약) | ✅ | `specs/` 3건 · 적대적 검사 33/33 |
| **KB** `record.json` + `kb/entries/` | ✅ 가동 | **37건** (런 8 · 런 없는 지식 29) |
| **L4** 수치 건전성 판정 | ✅ | [`bdbot/health.py`](../../bdbot/health.py) · [`tools/health.py`](../../tools/health.py) · `cli health` · 적대적 31/31 · 런 33/33 HEALTHY |
| 1-D 두 케이스 스크립트 (`chain-bend` · `trap-drag`) | ⬅️ **다음** | — |
| 1-D 세 케이스 · 2 · 8-rest · 9 | 대기/차단 | §0.2 |

**케이스 실측** (`$PY -m bdbot.cli status`)

| 케이스 | L0 | L2 | 스크립트 | 런 |
|---|---|---|---|---|
| `trap-2d-5um` | READY | READY | ✅ | 2 |
| `soft-r3-2d-A-sweep` | READY | READY | ✅ | 8 |
| `abp-rod-2d-run-flip` | READY | READY | ✅ | 3 |
| `chain-bend-2d-oscill` | READY | READY | ❌ | 0 |
| `trap-drag-2d-hex300` | READY | READY | ❌ | 0 |

**검증 스크립트 7종 전부 PASS**: `verify_{1c_equivalence, bdbot, intake_guards,
l3_spec_gaps, nondim_guards, pair_table, skill_snippets, health}`

**아직 없는 것** (계획엔 있으나 미구현 — 의도적):
`PhysicsModule` 레지스트리(§5.6은 설계만) · L4 실행 레이어 · 스윕 오케스트레이션 ·
훅/서브에이전트/슬래시명령(Phase 8-rest) · SQLite KB(§7.1) · 로우데이터 Tier C

---

## 0.2 남은 것 — "스케치가 시뮬레이션으로 잘 넘어가는가"를 기준으로 (2026-08-04)

**스케치 5장 중 3장**이 전 경로를 통과했습니다. §0.2-A 의 L2 차단은 **전부 해소됐습니다**
(`bdbot.cli status`: 5개 전부 L0·L2 READY).

| 스케치 | 인테이크 | 물리계 | 무차원화 | 실행 | 런 |
|---|---|---|---|---|---|
| `trap-2d-5um` | ✅ | ✅ | ✅ | ✅ | 2 |
| `soft-r3-2d-A-sweep` | ✅ | ✅ | ✅ | ✅ | 8 |
| `abp-rod-2d-run-flip` | ✅ | ✅ | ✅ | ✅ | 3 |
| `chain-bend-2d-oscill` | ✅ | ✅ | ⏸ | ⏸ | 0 — **스크립트 없음** |
| `trap-drag-2d-hex300` | ✅ | ✅ | ⏸ | ⏸ | 0 — **스크립트 없음** |

### A. ~~지금 막혀 있는 것~~ → **해소됨** (2026-08-04)

세 케이스의 물리계 결측(`U_ij` · 페어 퍼텐셜 · `R`)이 사람 확인으로 채워져
전부 L2 READY 가 됐습니다. `abp-rod` 는 실행까지 갔습니다(3런).

**다만 남은 판단 하나** — `abp-rod` 의 **이방성 병진 마찰은 BD로 불가**로 결론났습니다
(HI 부재. `γ⊥/γ∥ = 1.000000` 실측). 등방 평균 `γ̄` 를 쓰고 있으며, 그 한계가
단시간 MSD 에 남습니다 (§20 A · [`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md) §5).

### B. 파이프라인 자체에 남은 것

> **방향 (사용자 지시 2026-08-04)**: 엄밀한 시뮬레이션 실험·검증보다 **앞단(L0→L3) 구축**에
> 먼저 집중한다. 즉 "스케치 한 장을 던지면 물리계·무차원화까지 절차대로 나오는가"를
> 먼저 완성하고, 골든 테스트 정식화(Phase 2)와 로우데이터·스윕은 뒤로 미룬다.
>
> **방향 보강 (사용자 지시 2026-08-04, L3 작업 중)** ⭐️: 런 **후**의 과학적·물리적 검증은
> 비중을 낮춘다 — 계마다 다르고, 보고된(문헌에 있는) 계가 아닐 수 있다. 런 후에는
> **수치해석적 오류만** 본다: 발산 · NaN/Inf · 이상한 값으로의 수렴. 그 대신 **스케치에서
> 파라미터를 뽑아 무차원화를 제대로 해서 시뮬레이션을 잘 세우는 것**이 훨씬 중요하다.
>
> → 그래서 L4는 "물리 검증기"가 아니라 **수치 건전성 판정기**로 만든다. 해석해 대조·문헌
> 대조 기계를 더 만들지 않는다 (`bdbot.metrics` 의 `role` 체계로 이미 충분하다).
> → 엄밀성 투자는 L0→L3 쪽으로 몰아준다. `bdbot.nondim` 의 `validate()` 층(원장 완전성 ·
> 무차원수가 정말 그 비인가 · 역변환 가능성)이 이 지시의 직접적 산물이다 (§6.4).

**✅ 앞단 (완료 2026-08-04)**

| # | 무엇 | 산출물 | 결과 |
|---|---|---|---|
| 1 | `Observation` 스키마 + `intake check` | [`bdbot/intake.py`](../../bdbot/intake.py) | 스키마를 **5개 파일 실사용 빈도에서 도출**. `resolution` 키가 47/47 전부에 있고 값은 자주 `null` — 그게 §8.3의 장치. 5개 전부 오류 0 |
| 2 | `PhysicalSystem` 로더 공통화 | [`bdbot/physical.py`](../../bdbot/physical.py) | `derived_from` 불변식 · tier 게이트(§12.4) · **유도값 재계산 검증** · "L0이 BLOCKED면 L2는 존재할 수 없다" 교차검사 |
| 3 | skill `bd-intake` | [`.claude/skills/bd-intake/SKILL.md`](../../.claude/skills/bd-intake/SKILL.md) | 지어내기 억제 8규칙. `physical`/`choice` 구분 · 자기일관성 역산 · 차원 모호성 처리 |
| 4 | `bdbot` CLI | [`bdbot/cli.py`](../../bdbot/cli.py) | `status`·`intake`·`system`·`nondim`·`run`. 종료코드 0/1/2/3 으로 판정을 노출 (훅이 읽을 표면) |
| 5 | **L3 `NondimSpec` + `specs/`** | [`bdbot/nondim.py`](../../bdbot/nondim.py) | §6의 ③④⑤를 검사 가능하게 만든 것. L2↔L4 **유일한 계약**. 케이스 3종 이관 · 적대적 검사 30/30 · 세 케이스 공통 스키마 16키 |

**L3에서 도구가 잡아낸 것** (§6.4 / `scratch/verify_l3_spec_gaps.py`)
- ⭐️ **`run_id`가 물리계를 덮지 않았다** — 1-B 스펙에 물리계가 없어서 `d` 5µm→0.5µm·
  `η` 62배로 바꿔도(τ_B 16.1배 차이) run_id가 동일. `prepare_outdir`가 "이미 완료된 런"으로
  건너뛰어 **예전 계의 결과를 새 계의 결과로 보고**한다. → `system`을 해시에 넣음
- **무차원수가 정말 그 비인지 검사할 수 없었다** → `Group(num,den)` 이 원장을 가리키고
  `validate()`가 재계산 대조. **첫 실사용에서 라벨 오류 1건 즉시 검출**:
  `U(d)`를 "A kT"로 적어뒀으나 실제는 `(A+ε_WCA)kT` (A=100에서 101, 1% 어긋남)
- **기준 스케일이 원장에 없었다** — 리포트는 기준 길이를 `d`라고 적는데 타원체 케이스
  원장에는 `d_eq`만 있었다. 기호 불일치를 검사가 잡음
- **원장에서 빠진 스케일은 돌지 않는 검사다** — 1-B가 `dt`·`T_obs`를 원장 밖에 둬서
  시간척도 정렬표에 안 보였다. → 필수 **역할** 4종 강제

**앞단에서 도구가 잡아낸 것** (전부 도구를 실제 데이터에 돌려서 나왔습니다)
- **판정이 틀렸다** — 이미 완주한 2개 케이스가 BLOCKED로 나왔고, 원인은 `missing_required`에
  **물리적 미지값과 시뮬레이션 선택이 섞여** 있던 것. → `kind: physical|choice` 추가
- **검사기가 크래시했다** — 단위를 `furlong^2`로 바꾸니 오류 보고 대신 pint 예외.
  → `verify_intake_guards.py` (적대적 테스트 27종)를 만들어 잡음
- **`run_id`가 문서 수정에 반응했다** — `derived_from` 추가만으로 1-A run_id 변경.
  → `runid.physics_only()` 로 물리 필드만 해시. 재실행해 96필드 동일 확인
- **`derived_from`이 주석에만 있었다** — 기계가 못 읽는 출처. 필드로 승격

**남은 것**

| # | 무엇 | 왜 미루나 | 예상 |
|---|---|---|---|
| ~~4~~ | ✅ **L4 수치 건전성 판정기 — 완료** | `Guard`(런타임 즉시중단) · `judge_series`(발산·NaN·정지·붕괴) · **`step_health`(L3 원장 되먹임)** · `gate`(실행 전, 손댄 스펙 거부). 물리 검증기가 아님 | — |
| **4b** | `chain-bend` · `trap-drag` 케이스 스크립트 ⬅️ **다음** | L2 는 READY 인데 스크립트가 없다. **원칙 9.1대로 쪼갠다** — `chain-bend` 는 ① 평형 사슬(구동 OFF) → ② 구동 추가 | 각 반나절 |
| 5 | Phase 2 골든 테스트 정식화 | 검증이 `scratch/verify_*.py` 4종으로 이미 돌고 있다. `pytest` 정식화는 앞단이 안정된 뒤 | 1일 |
| 6 | Phase 8-rest 훅 + 슬래시명령 | CLI(#4)가 있어야 지킬 대상이 생긴다 | 1일 |
| 7 | Phase 3 로우데이터 Tier B~D | 지금 Tier A(GSD 위치)만으로 충분했다 | 1.5일 |
| 8 | Phase 4 나머지 | LLM 서술 · 규칙 승격(표본 부족) · 승인 원장 | 0.5일 |
| 9 | Phase 5 SQLite KB + 문헌 증류 | 문헌 **0편**. `chain` 스케치의 "Eric Furst 논문"이 첫 시드 후보 | 2.5일 |
| 10 | Phase 6 병렬 스윕 | 1-B에서 7런을 손으로 띄웠다 | 1.5일 |

### C. 확인하지 않은 채 남은 것 (정직하게 열거)

- `soft-r3`의 **`A`가 무차원이라는 해석**은 스케치 목표에서 역추론한 것 (`system.yaml.not_verified`)
- **유한크기 CV3**: `N=400` 하나만 돌렸다. `N=900` 대조 미실시
- **초기조건 의존성**: RSA 랜덤 배치 하나만. 결정에서 시작한 런과 대조 안 함
- **문헌 대조 0건**: 2D `r⁻³` 녹음 전이의 문헌 `Γ` 값을 확인하지 않았다
- **§20 A 이방성 병진 마찰**: `abp-rod`의 전제. 옵션 A(등방 평균)로 진행 가능해 보이나 미확정
- ⭐️ **가설 검증이 프로젝트 전체에 0건**: 관측량 역할을 세면 `implementation_check` 10 ·
  `measurement` 9 · **`hypothesis` 0** (+구 런 미표기 34). 지금까지 한 것은 전부
  "파이프라인이 맞다"는 증명이고 **물리적 발견은 0건**입니다.
  파이프라인 구축 단계에서는 정상이고, **원칙 9.1과 위 사용자 지시("런 후엔 수치 건전성만")가
  같은 결론**을 가리킵니다 — 조합 결과에 기존 이론을 갖다 대지 않는다. 다만 언젠가
  발견을 원한다면 케이스 설계 때 "이론이 추가하는 가정"을 비워두지 않아야 합니다 (원칙 8).
- **`abp-rod` 의 구현검사 5종은 여전히 순환적**: 원칙 9.1 게이트를 통과시키려 `derivation` 을
  채웠으나, 그 내용이 전부 "내가 넣은 값이 그대로 나오는가" 입니다. 게이트는 순환을
  **드러내지만 막지는 못합니다** — 그 케이스에 `hypothesis` 관측량 설계가 별도로 필요합니다.
- **`PhysicsModule` 레지스트리는 설계만**(§5.6). 원칙 9의 `standalone_check()` 도
  따라서 미구현입니다. 케이스 3개로는 아직 모듈 경계가 확정되지 않았습니다.
- **함정 14(3D `rotational_diffusion` 규약)는 2D 케이스에만 무해**합니다. 3D 액티브를
  하는 순간 물립니다.
- `bdbot`에 **일부러 올리지 않은 것**: 평형 지표 · 관측량 · 검증 전략 · 지배 시간척도 선택 ·
  초기 배치 · 표본 수집 루프. 세 번째 케이스에서 또 나오면 그때 올린다

---

## 1. 확정된 스코프

| 항목 | 결정 |
|---|---|
| 자율성 | **Claude Code가 에이전트 런타임** (자연어 + 이미지 → 시뮬레이션), **초기에는 전 단계 사람 승인** |
| 실행 환경 | 로컬 맥북, CPU 전용 (10 physical cores, 16 GB RAM) |
| 대상 시스템 | ⭐️ **제한하지 않음.** 특정 물리계를 하드코딩하지 않고, **물리 모듈 레지스트리**로 조합·확장한다. 새 물리가 필요하면 모듈 파일 하나를 추가한다 (코어 수정 없음) → **[§5.6](#56-physicsmodule--확장의-단위-)** |
| 산출물 | GSD 궤적(로우데이터 포함) + 관측량 테이블 + 자동 분석 플롯 |
| **지식** | **문헌 증류 KB + 자체 런 경험 KB, 파라미터 제안의 근거로 사용** |
| **단위** | **차원 우선 (하드 불변식). 모든 계는 SI 물리계로 먼저 확정하고, 주요 길이·시간·에너지 스케일로 무차원화. 양방향 변환 필수.** |
| **로우데이터** | **위치·방향·힘·토크를 계층적으로 저장, 온디맨드 재분석** |

### 비스코프 (지금은 안 함)
HPC/SLURM · 멀티 GPU · 웹 UI · 자율 파라미터 최적화(능동학습) · 고분자/rigid body/patchy

---

## 2. 핵심 설계 원칙

### 원칙 1 — LLM은 코드를 쓰지 않는다. **스펙과 지식**을 쓴다 ⭐️
```
❌ 자연어/그림 → LLM이 HOOMD 스크립트 생성 → exec()
✅ 자연어/그림 → LLM이 PhysicalSystem(JSON) 생성 → 결정론적 무차원화 → SimSpec → 빌더
```
재현성(`run_id = sha256(SimSpec)`), 검증 가능성, 디버깅 가능성, 안전성, 비용 전부가 여기서 나옵니다.

### 원칙 2 — **모든 숫자는 출처를 갖는다** ⭐️ (신규)
KB에서 온 파라미터든 LLM이 추정한 값이든, **어디서 왔는지 추적**합니다.

```python
particle_diameter = Provenanced(
    value=1.2 * ureg.micrometer,
    source="user_sketch:2026-08-03/note1.png#annotation_3",
    confidence=1,              # 사람이 확인함
)
solvent_viscosity = Provenanced(
    value=1.0e-3 * ureg.Pa * ureg.s,
    source="kb:paper/10.1103-PhysRevLett.110.238301#table1",
    confidence=2,              # 논문 추출, 미검증
)
```
- LLM이 논문을 증류하는 순간 **환각 위험**이 생깁니다. 출처 추적이 유일한 방어선입니다.
- 스펙 리포트에 "이 스펙의 각 값이 어디서 왔는가" 계보 테이블을 항상 첨부합니다.

### 원칙 3 — **차원이 먼저다. 무차원화는 그 다음이다.** ⭐️⭐️ (불변식)

이 프로젝트의 **하드 불변식(hard invariant)** 입니다. 우회 경로를 두지 않습니다.

```
        모든 계는 반드시 이 순서를 통과한다 — 예외 없음

PhysicalSystem (SI, pint Quantity, 모든 필드에 차원 있음)
      │
      │  ① 스케일 원장 작성 — 계에 존재하는 모든 길이·시간·에너지 스케일을 열거
      │  ② 기준 스케일 선택 — 주요 길이 σ*, 주요 시간 τ*, 주요 에너지 E* 를 명시적으로 지정
      │  ③ 무차원수 유도 — 모든 무차원수를 "두 스케일의 비"로 유도 (해석이 따라옴)
      │  ④ 스케일 분리 검사 — 무시하려는 스케일이 정말 분리되어 있는지 검증
      ▼
SimSpec (reduced units) + DimensionlessReport + ScaleLedger
      │
      ▼  시뮬레이션 → 무차원 결과
      │  ⑤ redimensionalize() — 역변환
      ▼
물리 단위 결과 (µm²/s, Pa, s, ...)
```

**세 가지 강제 사항:**

1. **`SimSpec`은 단독으로 존재할 수 없다.** 모든 `SimSpec`은 `derived_from: PhysicalSystemRef`를 가지며,
   그것을 만들어낸 `PhysicalSystem`과 `ScaleLedger`가 함께 저장됩니다.
   손으로 무차원 스펙을 바로 쓰는 경로는 **없습니다** — 무차원 값으로 시작하고 싶어도
   먼저 그 값이 함의하는 물리계를 명시해야 합니다 (§6.6 역구성 경로).

2. **기준 스케일은 자동으로 정해지지 않는다.** 어떤 길이를, 어떤 시간을 기준으로 삼았는지
   `ScaleLedger.reference`에 **선택 근거와 함께** 기록되고, 사람 확인 #3의 대상이 됩니다.

3. **역변환이 없으면 KB에 못 넣는다.** 무차원 결과만으로는 문헌·실험과 대조할 수 없습니다.
   `observables.parquet`은 무차원 값과 물리 단위 값을 **항상 쌍으로** 저장합니다.

> 이 원칙을 지키면 얻는 것: 무차원수가 "어디서 온 숫자"가 아니라 "두 물리 스케일의 비"가 되어
> 물리적 해석이 자동으로 따라오고, 스케일 분리 검사(§6.4)라는 강력한 자동 검증이 가능해집니다.

### 원칙 4 — 사람이 승인하고, 승인 이력이 자율성의 근거가 된다 (신규)
초기에는 `Intake → PhysicalSystem → SimSpec` 각 단계마다 사람이 확인합니다.
승인/거부/수정 이력을 **승인 원장(approval ledger)** 에 남깁니다.
나중에 "이 유형은 30번 중 30번 무수정 승인됨 → 자동화 후보"를 데이터로 판단합니다.

### 원칙 5 — 실패는 빠르고 크게, 그리고 **기록된다**
실행 전 검증 / 실행 중 감시 / 실행 후 **구조화된 사후분석**. 실패 경험이 다음 검증의 규칙이 됩니다.

### 원칙 6 — 맥북 CPU 현실 직시
HOOMD CPU 실행은 **단일 런 = 단일 코어**. 병렬성은 "여러 런 동시 실행"에서 나옵니다.
워커 기본 8개, 현실적 규모 N = 10³~2×10⁴, 10⁶~10⁸ steps.

---

### 원칙 7 — **물리 주장은 검증하고 말한다** ⭐️ (신규, v0.4)

이 프로젝트에서 추론으로 단정했다가 실측에서 뒤집힌 사례가 **두 번** 있었습니다:

| 주장 | 실측 |
|---|---|
| "트랩 검증 통과" (k=10만 확인) | k=2에서 **+1856%** — 최소 이미지 누락 (§11 함정 7) |
| "강체로 묶으면 마찰이 이방성" | `γ⊥/γ∥ = 1.000000` — BD엔 HI가 없음 |

둘 다 **그럴듯했고 둘 다 틀렸습니다.** HOOMD 동작이나 물리 결과에 대한 주장은
추론이 아니라 **실행으로** 확인합니다. 확인 못 했으면 "확인 안 함"이라고 씁니다
— §10의 `record.json.not_verified` 필드가 이걸 강제합니다.

### 원칙 8 — **"검증"은 기존 가설과 일치한다는 뜻이 아니다** ⭐️⭐️ (신규, 2026-08-04)

`intake/` 의 계들을 계산하는 이유는 **기존 가설과 다를 수 있기 때문**입니다.
그런데 판정 로직을 "예측과 다르면 FAIL"로 짜면 **발견을 실패로 부릅니다.**

모든 대조에 역할을 붙입니다 (`bdbot.metrics.ROLES`, `judge()`):

| 역할 | 예측의 출처 | 불일치 |
|---|---|---|
| `implementation_check` | **구현한 모델에서 해석적으로 유도** | **버그** → FAIL |
| `hypothesis` | 시뮬레이션이 **부과하지 않은** 가정 (연속체·희박극한·유효매질·문헌 모델) | **결과** — FAIL 아님 |
| `measurement` | 없음 | 시뮬레이션이 답 |

**실제로 물렸습니다.** `abp-rod` 의 예측 5종(`τ_eff`, `D_eff`, `D̄`, `τ`, 텀블 빈도)이
0.66% 이내로 일치했지만, 다섯 개가 **전부** 제가 구현한 모델(활성힘 + 회전확산 updater +
포아송 텀블 + BD)에서 유도되는 것이었습니다. 일치는 코드 검증이고 물리 발견이 아닙니다
— 거의 순환입니다. **그 케이스의 가설 검증은 0건이었습니다.**

따라서 케이스 설계 단계에서 두 목록을 따로 적습니다:

```
시뮬레이션이 부과하는 가정   과감쇠 BD · 등방 병진마찰 · 독립 회전확산 · 포아송 텀블 · 희박
이론이 추가하는 가정         ← 여기가 비면 발견이 없다
```

**역할은 결과를 보기 전에 정합니다.** 결과를 본 뒤 역할을 바꾸면 (`hypothesis` 불일치를
`measurement` 로 강등하거나, 우연한 일치를 `implementation_check` 성공으로 올리거나)
이 원칙이 무력해집니다. 시점과 구성 수준 규칙은 **원칙 9.1·9.2** 참조.

이 원칙은 §12(검증 레이어)와 §16 Phase 2(골든 테스트)의 의미도 바꿉니다 —
**골든 테스트는 `implementation_check` 만 담당하고, 물리적 발견은 별도 설계가 필요합니다.**

### 원칙 9 — **독립적인 요소는 하나씩 떼어 검증한다** ⭐️⭐️ (신규, 2026-08-04)

계에 서로 독립적인 요소가 여럿이면, **각 요소만 켠 최소 구성**을 먼저 돌려 단독으로
맞는지 확인하고 나서 조합합니다.

**왜 — 세 가지, 그중 첫째가 결정적**

1. **단독 구성에는 해석해가 있는 경우가 많습니다. 합치면 없습니다.**
   격리는 디버깅 편의가 아니라 **ground truth를 만드는 유일한 방법**일 때가 있습니다.
2. **틀렸을 때 범인이 특정됩니다.** 조합에서 틀리면 N개를 동시에 의심해야 하고,
   단독 버그인지 상호작용 탓인지도 구분 못 합니다.
3. **조합에서만 나타나는 것을 진짜 상호작용으로 식별**할 수 있습니다.
   단독 검증이 없으면 `cross_check`(§5.6)가 의미를 갖지 못합니다.

**증거 — 이 프로젝트에서 실제로**

| 단독 검증 | 해석해 | 결과 |
|---|---|---|
| `external.harmonic_trap` (상호작용 OFF, 독립 입자 N개) | `⟨x²⟩=kT/k` · `τ=γ/k` · Lorentzian | 관측량 4종 전부 일치 |
| `shape.rigid_rod` 마찰 (등속력 → 종단속도) | `v=F/γ` | `γ⊥/γ∥ = 1.000000` — **예상 뒤집음** |
| 적분기 편향 (선형계) | `편향=(dt/τ)/2` | 이론과 0.07% |
| `active.abp` (pair·trap·shape OFF) | `D_eff = D_t + v₀²/(d·Λ)` | **함정 2건 발견** ↓ |

`active.abp` 단독 검증(`scratch/standalone_abp_diffusion.py`)에서 나온 것:

- **활성력이 0이면 회전확산이 아예 동작하지 않습니다** — Λ=0 (4개 조합 전부).
- **HOOMD의 `rotational_diffusion`은 director 감쇠율 Λ 그 자체입니다.**
  2D·3D 모두 `Λ/D_r = 1.00`. 표준 이론의 `(d−1)D_r`이 아니라 **3D에서 2배 어긋납니다.**
- 따라서 §17에 있던 `D_eff = D_t + v₀²/[2(d−1)D_r]`는 **2D에서만 우연히 맞았습니다**
  (3D 실측 +29~31% 불일치). 올바른 형태는 `D_t + v₀²/(d·Λ)` — 실측 오차 1.5% 이내.

**둘 다 조합 시뮬레이션 안에서는 못 잡았을 것입니다.** 활성력·상호작용·트랩이 섞이면
MSD 한 곡선에 전부 묻힙니다. 실제로 처음엔 `D_eff`만 재서 "식이 틀렸나" 했는데,
**Λ와 v₀를 따로 떼어 재고 나서야** 원인이 회전확산 쪽임이 드러났습니다.

**절차**

1. 모듈 M을 추가할 때 **M만 켠 최소 케이스**를 먼저 돌립니다
2. 해석해나 알려진 극한과 대조합니다. 대조가 안 맞으면 **더 쪼갭니다**
   (ABP → Λ 따로, v₀ 따로)
3. 통과해야 조합에 투입합니다
4. 조합이 실패하면 **관련 모듈의 단독 검증이 최신인지 먼저** 봅니다 (회귀)

`PhysicsModule`에 `standalone_check()`를 둡니다 (§5.6) — 최소 구성 · 예측 · 허용오차.

### 9.1 단계마다 **역할이 다릅니다** ⭐️⭐️ — 조합에 기존 이론을 갖다 대지 마세요

원칙 9의 두 단계는 인식론적 지위가 다릅니다. 섞으면 원칙 8을 어기게 됩니다.

| 단계 | 역할 | 왜 | 불일치의 뜻 |
|---|---|---|---|
| **모듈 단독** | `implementation_check` | 최소 구성에서는 **기존 이론이 곧 내 모델**이다. 자유 ABP의 `D_eff` 는 내가 구현한 방정식에서 유도된다 | **버그** |
| **조합 전체** | `measurement` 또는 `hypothesis` | 조합의 해석해는 대개 **없다** — 그래서 시뮬레이션하는 것이다 | **결과** |

**조합 결과에 기존 이론을 갖다 대는 것이 위험한 이유:**

- 맞으면 "검증됐다"고 말하게 되는데, 사실은 그 이론이 부과되지 않았음을 확인한 게 아니다
- 안 맞으면 "이 영역엔 그 이론이 안 맞는다"고 말하게 되는데, 이건 반증 불가능하다
- **어느 쪽이든 배우는 게 없고, 나쁘면 발견을 지워버린다**

조합 단계에서 기존 이론을 쓰려면 `hypothesis` 로 등록하고 **불일치를 결과로 보고**합니다.
`composite` + `implementation_check` 를 쓰려면 **왜 그 해석식이 조합에도 유도되는지**
(`derivation`) 를 반드시 적어야 합니다 — 대개 극한(희박·선형응답·단시간)입니다.

### 9.2 예측은 **결과를 보기 전에** 고정합니다 ⭐️

결과를 본 뒤 해석을 맞추면, 맞춘 것인지 맞은 것인지 구분할 수 없습니다.

**구조적 보장** (자기신고가 아니라): 예측 함수는 **시뮬레이션 결과를 인자로 받지 않습니다.**
`cases/*.py` 의 `analytic(lg)` 가 이미 그 형태입니다 — 스케일 원장만 받으므로
결과에 의존할 방법이 없습니다. 실행 전에 `predictions.json` 으로 떨궈 감사 가능하게 합니다.

**결과를 본 뒤 해석을 바꿔야 한다면** — 있을 수 있는 일입니다. 다만 조건이 붙습니다:

> 새 해석은 **실패한 그 관측량이 아닌 다른 관측량**으로 독립 확인되어야 합니다.
> 그리고 원래 예측이 실패했다는 사실을 기록에 남깁니다.

**실제 사례 (2026-08-04, `active.abp` 단독 검증)** — 이 규칙이 왜 필요한지:

| 단계 | 무슨 일 |
|---|---|
| 사전 등록 | `D_eff = D_t + v₀²/[d(d−1)D_r]` (실행 전 유도) |
| 결과 | 3D에서 **+90% 불일치 — 등록한 예측 실패** |
| 원인 추적 | `Λ`(director 감쇠율)를 **MSD가 아닌 자기상관**에서 독립 측정 → `Λ = D_r` |
| 재예측 | `D_eff = D_t + v₀²/(d·Λ)` → 4개 조건 −1.5~−0.6% |

곡선맞춤이 아닌 근거: ① `Λ` 가 **다른 관측량**에서 나왔다 ② 넣은 뒤 자유 파라미터 0
③ 4개 조건에서 동시에 맞았다 (맞췄다면 한 점만 0%이고 나머지는 나빠야 한다).

**진짜 교훈은 따로 있습니다.** 틀린 것은 물리식이 아니라 **등록되지 않은 파라미터 매핑**이었습니다.
`D_eff = f(D_r)` 로 등록했지만 실제 구조는 `D_eff = f(Λ)` 이고, `Λ = (d−1)D_r` 이라는
**보이지 않는 가정**이 하나 끼어 있었습니다. HOOMD는 `Λ = D_r` 이었습니다.

> **예측을 등록할 때 파라미터 매핑 가정도 함께 등록하세요.**
> "이 식의 `D_r` 은 도구의 `rotational_diffusion` 과 같다" 도 하나의 가정입니다.

**한계 — 정직하게**

- 단독으로 맞다고 조합이 맞지는 않습니다. `cross_check`가 별도로 필요합니다.
- **원리적으로 분리 불가능한 것이 있습니다.** `shape`는 `medium` 없이 마찰이 정의되지
  않습니다. §5.6의 의존성 DAG를 따라 "최소 구성"을 정의해야 합니다.
- **분리하면 물리적으로 무의미해지는 조합이 있습니다.** 위 함정 ①이 정확히 그 경우입니다 —
  활성력을 끄고 회전확산만 보려 했더니 HOOMD가 회전확산을 꺼버렸습니다.

**로드맵 영향** — 남은 케이스 중 둘이 모듈 2개를 동시에 넣습니다:

| 케이스 | 추가 모듈 | 단독성 | 제안 |
|---|---|---|---|
| `trap-2d-5um` ✅ | external | 완전 단독 | — |
| `soft-r3` ✅ | pair | 페어만 | — |
| `trap-drag` | external(검증됨) + driving | 사실상 driving 단독 | 그대로 |
| `chain-bend` | bonded + driving | **2개 동시** | ① 평형 사슬(구동 OFF) → ② 구동 추가 |
| `abp-rod` | shape + active | **2개 동시** | ① 구형 ABP ✅ **완료** → ② 형상 추가 |

## 3. 아키텍처 — 닫힌 루프

```
┌── L0 INTAKE (멀티모달) ────────────────────────────────────────┐
│  스케치 이미지 / 손메모 / 화이트보드 사진 / 논문 PDF / 자연어      │
│  → Claude Code Read 툴 → Observation(구조화) ─ [사람 확인 #1] ►│
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L1 KNOWLEDGE ────────────────────────────────────────────────┐
│  KB 검색: 유사 시스템 · 무차원수 좌표 근방 · 과거 우리 런         │
│  결과를 근거로 결측 파라미터 제안 (출처+신뢰도 명시)              │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L2 PHYSICAL SYSTEM (차원 있음, SI) ──────────────────────────┐
│  σ=1.2µm, T=298K, η=1mPa·s, φ=0.6, v₀=8µm/s, ...              │
│  각 필드에 Provenance 부착   ─── [사람 확인 #2, 필드별 승인] ──► │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L3 NONDIMENSIONALIZATION ────────────────────────────────────┐
│  ① 스케일 원장  ② 기준 선택(d, τ_B, kT)  ③ 비로 유도  ④ 분리 검사 │
│  → SimSpec + DimensionlessReport + ScaleLedger                 │
│  Pe=38.4, φ=0.60, D_r*=3.0, dt*=5e-5 ── [사람 확인 #3] ───────► │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L4 VALIDATION ── 물리·수치·비용 + KB 기반 경고 ───────────────┐
│  "과거 Pe>80, dt=1e-4 조합에서 3회 발산했음"                     │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L5 BUILD → L6 RUN (프로세스 풀 8) → L7 RAW DATA STORE ───────┐
│  GSD 궤적(위치/방향/힘/토크, 계층적 빈도) + HDF5 로그            │
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L8 ANALYZE (온디맨드) → L9 PLOT / 문헌 대조 ─────────────────┐
└──────────────────────┬─────────────────────────────────────────┘
                       ▼
┌── L10 POST-MORTEM ── 성공/실패 요인 구조화 ────────────────────┐
│  → KB에 환류 (다음 제안과 검증 규칙을 개선)  ──────► L1 로 순환  │
└────────────────────────────────────────────────────────────────┘
```

**중요**: L3~L9는 LLM 없이 CLI/파이썬 API로 단독 동작해야 합니다. LLM은 L0~L2와 L10의 "해석·판단" 부분만 담당합니다. 이렇게 해야 파이프라인을 LLM 없이 검증할 수 있습니다.

---

## 4. 디렉토리 구조

### 4.1 현재 (2026-08-04) — 케이스 주도로 실제로 존재하는 것

```
simulation_auto/
├── CLAUDE.md                      불변식 (항상 로드)
├── mater_plan.md · environment.yml
├── .claude/skills/{bd-hoomd,bd-physics}/SKILL.md
├── docs/hoomd_capabilities.md     HOOMD 능력 실측 매트릭스
├── intake/<case>/                 스케치 + observation.yaml + system.yaml
├── bdbot/                         ⭐️ Phase 1-C 산출물 — **두 번 나온 것만**
│   ├── units.py                   단일 pint 레지스트리
│   ├── provenance.py              Provenanced (값+출처+tier)
│   ├── materials.py               γ=3πηd · D_t=kT/γ · τ_B · m · τ_p
│   ├── scales.py                  ScaleLedger + thermal 기준
│   ├── checks.py                  Check(모델/적분/기하/통계) + dt = 10⁻²·γ/(국소 강성)
│   ├── report.py                  DimensionlessReport 렌더러
│   ├── runid.py                   콘텐츠 주소 지정 + 재실행 방지
│   ├── metrics.py                 metrics.json 스키마
│   ├── stats.py                   블록평균 · 자기상관 보정 · 불편 자기상관
│   └── sim.py                     2D 프레임 · BD 적분기 · GSD · 시드 · 최소이미지
├── cases/{trap_2d_5um,soft_r3_2d}.py   관통 스크립트 (케이스 고유 물리만 남음)
├── tools/{postmortem,kb}.py       사후분석 + KB 질의
├── scratch/*.py                   검증 스크립트 11종 (재현 가능)
└── runs/<run_id>/                 report·result·metrics·record·traj·plots
```

`bdbot` **CLI는 아직 없습니다** — 케이스가 직접 `bdbot`을 임포트합니다. CLI는 지킬 대상
(`bdbot run`)이 필요한 Phase 8-rest의 훅과 함께 만드는 것이 맞습니다.

**의도적으로 넣지 않은 것**: 평형 판정 지표 · 관측량 · 검증 전략 · 지배 시간척도 선택 ·
초기 배치 · 표본 수집 루프. 케이스마다 달랐습니다 (skill `bd-physics` §6.3).

### 4.2 목표 (Phase 1-C 이후 수렴할 형태)

아래는 **설계 의도**이며 현재 구현 상태가 아닙니다. Phase 1-C에서 1-A/1-B의
공통분모를 뽑을 때 참고하되, **실제로 두 번 나온 것만** 만듭니다.

```
simulation_auto/
├── mater_plan.md
├── environment.yml · pyproject.toml
│
├── bdbot/
│   ├── intake/                    # ⭐️ L0
│   │   ├── vision.py              # 이미지 → Observation (Claude vision)
│   │   ├── pdf.py                 # 논문 PDF → 증류 후보
│   │   ├── observation.py         # Observation 스키마
│   │   └── review.py              # 사람 확인 UI (터미널 폼)
│   │
│   ├── knowledge/                 # ⭐️ L1, L10
│   │   ├── store.py               # SQLite + FTS5
│   │   ├── schema.py              # KnowledgeEntry, Claim, Source, Provenance
│   │   ├── distill.py             # 논문 → KnowledgeEntry (LLM)
│   │   ├── search.py              # 키워드 + 무차원수 좌표 근방 검색
│   │   └── feedback.py            # 런 사후분석 → KB 환류
│   │
│   ├── physical/                  # ⭐️ L2
│   │   ├── system.py              # PhysicalSystem (pint 단위 부착)
│   │   ├── materials.py           # 물·글리세롤 점도, 실리카/PS 밀도 등 상수 라이브러리
│   │   └── provenance.py
│   │
│   ├── nondim/                    # ⭐️ L3
│   │   ├── scales.py              # 스케일 선택 전략
│   │   ├── forward.py             # PhysicalSystem → SimSpec
│   │   ├── inverse.py             # 무차원 결과 → 물리 단위
│   │   ├── groups.py              # Pe, φ, Re, Pé_r, λ_D/σ ... 무차원수 계산
│   │   └── report.py              # DimensionlessReport
│   │
│   ├── spec/                      # SimSpec (무차원)
│   ├── validate/                  # L4 (+ KB 연계 경고)
│   ├── build/                     # L5 (HOOMD 조립)
│   ├── run/                       # L6 (풀, 큐, 가드, 체크포인트)
│   ├── rawdata/                   # ⭐️ L7
│   │   ├── policy.py              # 저장 계층 정책 (tier A/B/C)
│   │   ├── writers.py             # GSD/Burst/per-particle 로거
│   │   ├── index.py               # 프레임 인덱스 (빠른 랜덤 액세스)
│   │   └── loader.py              # 지연 로딩, 메모리맵, 슬라이싱
│   │
│   ├── analyze/                   # L8 (레지스트리 + 캐시)
│   ├── plot/                      # L9 (+ 문헌 대조 오버레이)
│   ├── postmortem/                # ⭐️ L10
│   │   ├── taxonomy.py            # 실패 분류 체계
│   │   ├── diagnose.py            # 자동 진단 (수치 지표 기반)
│   │   └── narrate.py             # LLM 요약 → KnowledgeEntry
│   │
│   ├── approval/                  # ⭐️ 승인 원장
│   │   ├── ledger.py
│   │   └── gates.py
│   │
│   ├── distill_batch.py           # (선택) anthropic Batch API 대량 증류 — §15.9
│   └── cli.py                     # ⭐️ 유일한 진입점. Claude Code는 이것만 호출.
│
├── CLAUDE.md                      # ⭐️ 항상 로드되는 불변식 (§15.4)
├── .claude/                       # ⭐️ Claude Code 통합 표면 (§15.3)
│   ├── settings.json              #   권한 + 훅
│   ├── skills/                    #   bd-physics, bd-hoomd, bd-intake, bd-distill
│   ├── agents/                    #   bd-distiller, bd-analyst, bd-reviewer
│   ├── commands/                  #   /bd-intake, /bd-spec, /bd-run, ...
│   └── hooks/                     #   guard_invariant, guard_separation, guard_cost
│
├── kb/
│   ├── knowledge.db               # SQLite (FTS5)
│   ├── papers/                    # 원본 PDF
│   └── figures/                   # 추출한 그림
│
├── intake/<case>/                 # 스케치/메모 (사용자가 던져 넣음)
│   ├── sketch_01.jpeg
│   ├── observation.yaml           # L0 산출물 (사람 확인 후 확정)
│   └── system.yaml                # L2 PhysicalSystem (사람 확인 후 확정)
├── specs/                         # L3 산출물 (bdbot nondim이 생성, 손으로 안 씀)
├── runs/<run_id>/                 # 로우데이터 + 분석 결과
└── tests/
```

> **`cli.py`가 유일한 진입점**인 것이 중요합니다. Claude Code 세션 밖에서도(cron, 스크립트,
> 다른 사람) 똑같이 동작해야 하며, 이것이 §15.10에서 지적한 유일한 위험에 대한 방어책입니다.

---

## 5. 데이터 모델 (5개 핵심 객체)

### 5.1 `Observation` — L0 산출물
```python
class Observation(BaseModel):
    """스케치/메모/그림에서 읽어낸 것. LLM이 채우고 사람이 확인."""
    source_files: list[str]
    raw_transcription: str              # 이미지에서 읽은 글자 그대로
    system_guess: str                   # "2D active colloid monolayer"
    entities: list[Entity]              # 입자, 벽, 장, 화살표 …
    stated_quantities: list[StatedQuantity]  # "d ≈ 1 µm", "Pe 20~100"
    stated_goals: list[str]             # "MIPS 상경계 찾기"
    ambiguities: list[str]              # ⭐️ LLM이 확신 못 한 것 명시
    unread_regions: list[str]           # ⭐️ 읽지 못한 부분 명시
```
`ambiguities` / `unread_regions`를 **강제 필드**로 둡니다. 모르는 걸 모른다고 말하게 하는 장치입니다.

### 5.2 `KnowledgeEntry` — KB 단위 레코드
```python
class Source(BaseModel):
    kind: Literal["paper", "book", "our_run", "user_input", "handbook"]
    doi: str | None; title: str | None; authors: list[str]; year: int | None
    locator: str | None                 # "Fig.3b", "Table 1", "p.4 eq.(7)"
    local_path: str | None

class Claim(BaseModel):
    statement: str                      # "MIPS onset near Pe≈35 at φ=0.6 (2D ABP)"
    dimensionless_coords: dict[str, float]   # {"Pe": 35, "phi": 0.6}
    kind: Literal["parameter", "phase_boundary", "scaling", "method_note", "pitfall"]

class KnowledgeEntry(BaseModel):
    id: str
    system_tags: list[str]              # ["2D", "ABP", "WCA", "monodisperse"]
    source: Source
    claims: list[Claim]
    physical_params: dict               # 차원 있는 값 (단위 문자열 포함)
    dimensionless_params: dict          # Pe, φ, D_r*, ...
    numerics: dict                      # 그 논문이 쓴 dt, N, 적분기, 박스
    confidence: Literal[0, 1, 2, 3]     # 아래 표
    extracted_by: str                   # "claude-opus-5 / distill_prompt_v3"
    verified_by: str | None; verified_at: datetime | None
    notes: str
```

**신뢰 등급 (confidence tier)** ⭐️
| tier | 의미 | 프로덕션 런에 단독 사용 |
|---|---|---|
| 0 | 사람이 직접 입력 / 물성 핸드북 | ✅ |
| 1 | 문헌 추출 + 사람 검증 완료 | ✅ |
| 2 | 문헌 추출, 미검증 (LLM 단독) | ⚠️ 승인 필요 |
| 3 | 우리 시뮬레이션에서 귀납 | ⚠️ 승인 필요 |

검증기가 "tier 2 이하 값만으로 구성된 스펙"에는 반드시 사람 승인을 요구합니다.

### 5.3 `PhysicalSystem` — 차원 있는 물리계
```python
class PhysicalSystem(BaseModel):
    label: str
    dimensions: Literal[2, 3]

    # 입자
    particle_diameter: Provenanced[Quantity]        # µm
    particle_density: Provenanced[Quantity] | None  # kg/m³
    n_particles: Provenanced[int] | None
    area_or_volume_fraction: Provenanced[float]

    # 매질
    temperature: Provenanced[Quantity]              # K
    solvent_viscosity: Provenanced[Quantity]        # Pa·s
    # → γ = 3πηd (Stokes) 로 유도, 또는 직접 지정
    drag_coefficient: Provenanced[Quantity] | None

    # 상호작용
    interaction: Provenanced[InteractionSpec]       # WCA/LJ/Yukawa + ε(J or kT), κ(1/m)

    # 액티브
    self_propulsion_speed: Provenanced[Quantity] | None   # µm/s
    rotational_diffusion: Provenanced[Quantity] | None    # 1/s

    # 관측 목표
    target_observables: list[str]
    target_physical_time: Provenanced[Quantity]     # s (얼마나 오래 볼 것인가)
```
모든 필드가 `Provenanced` → 값 + 출처 + 신뢰도.
`pint`로 단위를 강제해 단위 실수를 차단합니다.

### 5.4 `SimSpec` — 무차원 스펙 (HOOMD 직행)
```python
class SimSpec(BaseModel):
    schema_version: str = "0.2"
    label: str

    # ⭐️ 불변식: 무차원 스펙은 반드시 물리계에서 유도된다
    derived_from: PhysicalSystemRef          # 해시 + 저장 경로
    scale_ledger: ScaleLedger                # 기준 스케일 + 전체 스케일 목록
    dimensionless: dict[str, float]          # {Pe, phi, D_r_star, dt_star, Re, St, ...}

    box: BoxSpec
    types: list[ParticleType]
    modules: list[ModuleSpec]                # ⭐️ 물리 모듈의 조합 (§5.6)
    integrator: IntegratorSpec               # 전부 reduced units
    run: RunSpec
    raw_data: RawDataPolicy
    observables: list[ObservableSpec]
```
`derived_from`이 없는 `SimSpec`은 **빌더가 거부합니다**. 무차원 값으로 시작하고 싶다면
§6.6의 역구성 경로를 통해 먼저 물리계를 명시해야 합니다.

### 5.5 `RunRecord` + `PostMortem`
```python
class PostMortem(BaseModel):
    run_id: str
    outcome: Literal["success", "partial", "failure"]
    failure_modes: list[FailureMode]    # 아래 taxonomy
    diagnostics: dict                   # 에너지 드리프트, 평형 판정, 유한크기 지표 …
    narrative: str                      # LLM 요약
    lessons: list[Claim]                # ⭐️ KB로 환류될 항목
    dimensionless_coords: dict          # 어느 파라미터 영역이었는지
```

---

### 5.6 `PhysicsModule` — 확장의 단위 ⭐️⭐️

**스코프를 제한하지 않기로 했으므로, 물리계를 하드코딩하지 않고 모듈로 조합합니다.**
새 물리가 필요하면 **모듈 파일 하나를 추가**하고 레지스트리에 등록합니다. 코어는 건드리지 않습니다.

핵심은 하나의 모듈이 **7개 레이어 전부에 스스로 기여**한다는 것입니다. 그래야 새 물리를 붙여도
스케일 원장·무차원화·분리 검사·검증이 자동으로 따라옵니다 (원칙 3이 확장에도 유지됨).

```python
class ModuleContribution(BaseModel):
    scales:  dict[str, Quantity]          # ① 원장에 추가할 특성 스케일 (§6.1)
    groups:  list[DimensionlessGroup]     # ② 이 모듈이 만드는 무차원수 (§6.3)
    checks:  list[SeparationCheck]        # ③ 이 모듈이 요구하는 분리 검사 (§6.4)
    reduced: dict[str, float]             # ④ 무차원 파라미터

class PhysicsModule(ABC):
    kind:     ClassVar[str]               # "external.harmonic_trap"
    requires: ClassVar[set[str]] = set()  # 의존 모듈 (예: {"shape.*"})
    PhysicalParams: ClassVar[type[BaseModel]]   # 차원 있는 파라미터 스키마
    ReducedParams:  ClassVar[type[BaseModel]]   # 무차원 파라미터 스키마

    @abstractmethod
    def contribute(self, phys, ctx: ScaleContext) -> ModuleContribution: ...  # L3
    @abstractmethod
    def build(self, sim, spec, ctx: BuildContext) -> None: ...                # L5

    def standalone_check(self) -> StandaloneCheck | None: ...  # ⭐️ 원칙 9 단독 검증
    def cross_check(self, others) -> list[Issue]: return []   # L4 조합 검사
    def default_observables(self) -> list[str]: return []     # L8
    def periodic_safe(self) -> bool: return True              # ⭐️ 함정 7 (§11)

class StandaloneCheck(BaseModel):
    """이 모듈만 켠 최소 구성 + 해석해 + 허용오차 (원칙 9)."""
    minimal_spec: dict          # 다른 모듈 전부 OFF
    predictions: dict           # 관측량 → 해석식
    tolerance_pct: float
    script: str | None = None   # 재현 스크립트 경로
```

#### 해결 순서 (의존성 DAG)

```
① medium.*   →  η, T, kT, ρ_f
② shape.*    →  마찰 텐서 γ∥, γ⊥, γ_r  →  D_t, D_r     ← 구/타원체 분기점
③ 나머지 모듈 →  각자 스케일·무차원수·검사 기여
④ 기준 스케일 선택 (§6.2)
⑤ 무차원화 + 전체 분리 검사 (§6.4)
```

`shape`가 먼저 풀려야 하는 이유: **`D_r = 3D_t/d²`는 구에만 성립**합니다. 이건 코어 공식이 아니라
`shape.sphere` 모듈의 기여입니다. `shape.ellipsoid`는 Perrin 마찰 인자를 대신 기여합니다.
이렇게 두면 `abp-rod` 같은 케이스가 코어 수정 없이 들어옵니다.

> ⚠️ **실측 결과 (2026-08-03)**: 강체(`constrain.Rigid`)로 구를 묶어도, 본드로 비드사슬을 만들어도
> **병진 마찰은 등방입니다** (`γ⊥/γ∥ = 1.000000`, 두 경로 모두). 막대의 `γ⊥/γ∥ → 2`는
> **유체역학적 상호작용(HI)의 효과**이고 BD에는 HI가 없기 때문입니다 — HOOMD의 한계가 아니라
> **BD 모델 자체의 성질**입니다. 회전 마찰은 `gamma_r` 텐서로 정상 동작합니다(비 0.2000 실측).
> → `shape.*` 모듈은 `translational_friction: "isotropic_average" | "anisotropic"` 을 선언하고,
> 검증기가 "단시간 MSD 이방성이 목표인데 등방 근사를 씀" 경고를 냅니다.
> 근거·선택지: [`docs/hoomd_capabilities.md` §5.1–5.4](../../docs/hoomd_capabilities.md)

#### 모듈이 기여하는 것 — 예시

| 모듈 | 추가 스케일 | 무차원수 | 분리 검사 |
|---|---|---|---|
| `shape.sphere` | `d` | — | — |
| `shape.ellipsoid` | `a`(장반축) `b`(단반축) | `p = a/b` 종횡비 | ⚠️ 이방성 병진확산 **BD 미지원** (아래) |
| `pair.wca` | `r_c` | `ε*` | `r_c < L/2` |
| `pair.table` | `r_c`, 퍼텐셜 스케일 | `A*` | `dt/τ_int ≤ 1e-2` |
| `external.harmonic_trap` | `ℓ_k=√(kT/k)`, `τ_k=γ/k` | `k* = k d²/kT` | `dt/τ_k ≤ 1e-2`, `ℓ_k < L/2` |
| `bonded.bond_harmonic` | `ℓ_b=√(kT/k_b)`, `τ_b=γ/k_b` | `k_b*` | `dt/τ_b ≤ 1e-2` |
| `bonded.angle_harmonic` | — | `κ* = κ_bend/kT` (지속길이) | `ℓ_p^chain ≤ L/4` |
| `active.abp` | `ℓ_p=v₀/D_r`, `τ_v`, `τ_r` | `Pe`, `D_r*` | `dt·D_r ≤ 1e-2`, `ℓ_p ≤ L/4` |
| `active.run_and_flip` | `τ_flip` | `Pe`, `τ_flip/τ_B` | `dt/τ_flip ≤ 1e-2` |
| `driving.oscillate` | `τ_ω = 1/ω` | `De = τ_relax·ω` (Deborah) | `dt/τ_ω ≤ 1e-2`, `T_obs ≫ τ_ω` |

> 각 행이 "새 물리를 붙이면 스케일 원장·무차원수·분리 검사가 **자동으로** 늘어난다"는 뜻입니다.
> 원칙 3의 이점이 확장에도 그대로 유지됩니다.

#### `ModuleSpec` — 스펙 안에서의 표현

```python
class ModuleSpec(BaseModel):
    kind: str                    # 레지스트리 키. "external.harmonic_trap"
    params: dict                 # 해당 모듈의 ReducedParams
    targets: str = "all"         # 적용 대상 필터 (hoomd.filter 표현)
```

레지스트리는 `bdbot/modules/__init__.py`의 데코레이터 등록:
```python
@register("external.harmonic_trap")
class HarmonicTrap(PhysicsModule): ...
```

#### 현재 상태 — 실측 기반

어떤 모듈을 만들 수 있는지는 추측이 아니라 조사 결과에 근거합니다:
**[`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md)** — hoomd 7.1.0 실측,
15개 API 실동작 검증 완료, `intake/` 5개 케이스 전부 구현 가능 확인.

---

## 6. 무차원화 엔진 (L3) — 상세

원칙 3의 5단계(① 스케일 원장 → ② 기준 선택 → ③ 비로 유도 → ④ 분리 검사 → ⑤ 역변환)를 구현합니다.

### 6.1 ① 스케일 원장 (Scale Ledger) ⭐️

무차원화의 첫 단계는 **계에 존재하는 모든 특성 스케일을 SI 단위로 빠짐없이 열거**하는 것입니다.
무차원수를 먼저 계산하고 스케일을 나중에 유추하는 순서를 금지합니다.

```python
class ScaleLedger(BaseModel):
    lengths: dict[str, Quantity]     # 이름 → SI 값
    times:   dict[str, Quantity]
    energies: dict[str, Quantity]
    reference: ReferenceScales       # 기준으로 택한 셋 + 선택 근거
    separations: list[SeparationCheck]
```

**반드시 계산해 원장에 올릴 스케일** (해당 없으면 `None`으로 명시적 기록):

| 종류 | 이름 | 정의 | 의미 |
|---|---|---|---|
| **길이** | `d` | 입자 지름 | 기본 길이 |
| | `a_mean` | `ρ^(-1/dim)` = 평균 입자 간격 | 조밀/희박 판정 |
| | `r_c` | 상호작용 컷오프 | 힘의 도달 거리 |
| | `λ_D` | Debye 길이 (하전계) | 정전 스크리닝 |
| | `ℓ_p` | `v₀/D_r` (액티브) | 방향 지속 길이 |
| | `L` | 박스 한 변 | 시스템 크기 |
| | `ξ` | 상관 길이 (사후, 측정값) | 유한크기 판정 |
| **시간** | `τ_p` | `m/γ` | 운동량 이완 (관성) |
| | `dt` | 적분 시간 스텝 | 수치 해상도 |
| | `τ_v` | `d/v₀` (액티브) | 이류/탄도 시간 |
| | `τ_int` | `d²γ/ε` | 상호작용 응답 시간 |
| | `τ_r` | `1/D_r` | 방향 상관 시간 |
| | `τ_B` | `d²/D_t` | 확산(Brownian) 시간 |
| | `T_obs` | 프로덕션 물리 시간 | 관측 창 |
| **에너지** | `k_BT` | 열에너지 | 요동 척도 |
| | `ε` | 상호작용 깊이/세기 | 결합 세기 |
| | `f_a·d` | 자기추진 일 | 액티브 구동 |

기본 물성 유도식 (`physical/materials.py`가 제공, 전부 차원 있음):
```
γ    = 3πηd                      (Stokes 항력, 구)
D_t  = k_BT/γ                    (Stokes–Einstein)
D_r  = k_BT/(πηd³) = 3D_t/d²     (Stokes–Einstein–Debye, 구)
m    = ρ_p (π/6) d³              (입자 질량)
v₀   = f_a/γ                     (자기추진 속도)
```

### 6.2 ② 기준 스케일 선택 — 명시적, 근거 기록

```python
class ReferenceScales(BaseModel):
    length: tuple[str, Quantity]     # ("d", 1.2 µm)
    time:   tuple[str, Quantity]     # ("tau_B", 3.96 s)
    energy: tuple[str, Quantity]     # ("kT", 4.11e-21 J)
    rationale: str                   # 왜 이걸 골랐는지 (사람 확인 대상)
```

| 전략 | 길이 | 시간 | 에너지 | 언제 |
|---|---|---|---|---|
| `thermal` **(기본)** | `d` | `τ_B = d²/D_t` | `k_BT` | 열요동 지배 — 콜로이드, ABP |
| `interaction` | `d` | `τ_int = d²γ/ε` | `ε` | 상호작용 지배 — 강한 결합, 저온 |
| `active` | `ℓ_p` | `τ_r` | `f_a·ℓ_p` | 액티브 지배 — 매우 큰 Pe |
| `custom` | 사용자 지정 | | | 문헌 관례를 따를 때 |

`thermal` 전략에서 HOOMD 스펙은 `σ=1, kT=1, γ=1` → `τ_B = 1`, `D_t = 1`이 됩니다.

⚠️ **기준을 바꾸면 같은 물리계가 전혀 다른 무차원수를 갖습니다.** 문헌과 무차원수를 비교할 때
상대 논문이 어떤 기준을 썼는지 확인해야 하며, `KnowledgeEntry`에도 기준을 기록합니다.

### 6.3 ③ 무차원수는 "두 스케일의 비"로 유도한다 ⭐️

숫자를 던지지 않고 **어떤 두 스케일의 비인지** 함께 기록합니다. 물리적 해석이 따라옵니다.

| 무차원수 | 스케일 비 | 표현 | 물리적 의미 |
|---|---|---|---|
| `φ` | 점유부피/전체부피 | `N v_p / V` | 밀집도 |
| `Pe` | `τ_B / τ_v` | `v₀d/D_t = f_a d/k_BT` | 이류 vs 확산 |
| `D_r*` | `τ_B / τ_r` | `D_r τ_B` (구면 Stokes → 3) | 회전 vs 병진 확산 |
| `ℓ_p/d` | `ℓ_p / d` | `Pe / D_r*` | 지속 길이가 입자 몇 개분인가 |
| `ε*` | `ε / k_BT` | | 결합 세기 vs 열요동 |
| `κ*` | `d / λ_D` | `κd` | 스크리닝 vs 입자 크기 |
| `dt*` | `dt / τ_B` | | 수치 해상도 |
| `St` | `τ_p / τ_B` | `m/(γτ_B)` | 관성 vs 확산 |
| `Re` | 관성 vs 점성 (유체) | `ρ_f v₀ d / η` | 유체 관성 |
| `N*` | `L / d` | `L/d` | 시스템이 입자 몇 개분인가 |
| `T*` | `T_obs / τ_B` | | 관측 창이 확산시간 몇 배인가 |

`groups.py`의 각 무차원수는 `DimensionlessGroup(name, value, numerator_scale, denominator_scale, expression, interpretation)`으로 저장되어, 리포트가 자동으로 "무엇 대 무엇"인지 설명합니다.

### 6.4 ④ 스케일 분리 검사 ⭐️⭐️ — 하드 게이트

무차원화의 진짜 가치는 여기 있습니다. **무시하기로 한 스케일이 정말 분리되어 있는가**를 자동 검증합니다.
원장의 시간 스케일을 크기순 정렬해 놓으면 위반이 눈에 보입니다.

**시간 스케일 기준 규칙 (하나로 통일)**: 적분 스텝은 **가장 빠른 물리 시간척도의 1% 이하**여야 하고,
관성 완화는 **적분 스텝보다 최소 100배 빨라야** 합니다. 모든 임계가 `10⁻²`로 통일되어 기억하기 쉽습니다.

| 검사 | 조건 | 등가 표현 | 위반 시 |
|---|---|---|---|
| **관성 무시 (모델)** | `τ_p / τ_dyn ≤ 10⁻²` | `τ_dyn` = 관심 최속 시간척도 | ❌ 하드 — 과감쇠 BD 무효, Langevin 필요. **`dt`와 비교하지 않음** — BD엔 관성이 없어서 `dt≫τ_p`가 불필요 (실측: `τ_p/dt=4000`에서도 0.38% 정확) |
| **이류 해상** | `dt / τ_v ≤ 10⁻²` | `v₀ dt ≤ 0.01 d` | ❌ 하드 — 한 스텝에 지름 1% 초과 이동 |
| **회전 해상** | `dt / τ_r ≤ 10⁻²` | `dt·D_r ≤ 0.01` | ❌ 하드 — ABP 방향 동역학 붕괴 |
| **상호작용 해상** | `dt / τ_int ≤ 10⁻²` | `dt·ε/(d²γ) ≤ 0.01` | ❌ 하드 — 힘 적분 부정확 |
| **저 Reynolds** | `Re ≤ 10⁻³` | | ❌ 하드 — 유체 관성 무시 불가 |
| **컷오프 vs 박스** | `r_c < L/2` | | ❌ 하드 — minimum image 위반 |
| **유한크기 (지속길이)** | `ℓ_p ≤ L/4` | | ⚠️ 경고 — 액티브 인공효과 |
| **유한크기 (상관길이)** | `ξ ≤ L/4` (사후 측정) | | ⚠️ 경고 — 유한크기 아티팩트 |
| **스크리닝** | `λ_D ≤ L/4` | | ⚠️ 경고 — 스크리닝 불완전 |
| **관측 충분** | `T_obs ≥ 10² · max(τ_B, τ_r)` | | ⚠️ 경고 — 통계 부족 |
| **여유 부족** | 하드 검사 중 한계의 1/5 이내 | | ⚠️ 경고 — 파라미터 상향 여지 없음 |
| **조밀/희박** | `a_mean / d` 보고 | | ℹ️ 정보 — 체제 판정 |

> Brownian 적분기는 O(δt) 오차이므로(§11 함정 2) 위 상한을 통과해도 **정확도**가 보장되진 않습니다.
> 논문급 결과에는 `dt`를 절반으로 줄인 수렴 확인(§17 수렴 연구)이 별도로 필요합니다.

`SeparationCheck(name, ratio, threshold, verdict, message)` 리스트가 `ScaleLedger`에 저장되고,
하나라도 ❌면 **검증기가 실행을 거부**합니다.

### 6.4b `dt`는 편향에서 역산한다 ⭐️ (v0.4, 실측 확인)

"충분히 작게"는 애매합니다. **선형계에서는 계통 편향을 정확히 계산할 수 있습니다.**
조화 트랩의 Euler–Maruyama 이산 정상분산:

```
⟨x²⟩_discrete = (kT/k) / (1 − h/2),   h ≡ dt/τ    ⟹   상대 편향 ≈ h/2
```

**`dt/τ`의 절반이 곧 계통 편향입니다.** HOOMD Brownian이 이 법칙을 따름을 실측
(`scratch/dt_convergence.py`, N=2000):

| `dt/τ` | 편향 측정 | 이론 | 차이 |
|---|---|---|---|
| 0.10 | 5.262% | 5.263% | −0.002% |
| 0.05 | 2.539% | 2.564% | −0.025% |
| 0.02 | 1.079% | 1.010% | +0.069% |
| 0.01 | 0.569% | 0.503% | +0.067% |

| 목표 정확도 | `dt/τ_fast` | 상대 비용 |
|---|---|---|
| 1% | 2e-2 | 5× |
| **0.5%** | **1e-2** ← §6.4 하드 게이트 | 10× |
| 0.1% | 2e-3 | 50× |
| 0.025% | 5e-4 | 200× |

→ **§6.4의 하드 게이트 `1e-2`는 "0.5% 편향"을 뜻합니다.**
비선형계는 닫힌 형태가 없으니 수렴 연구가 필요하나, O(dt) 스케일링은 유지되므로
두 점 Richardson 외삽이 잘 듣습니다. 상세: skill `bd-physics` §1.2.

### 6.5 `DimensionlessReport` — 사람 확인 #3의 대상

```
System: 2D active colloid  (label: abp_silica_1p2um)
Reference scales: length=d, time=τ_B, energy=kT   [strategy: thermal]
  rationale: 열요동 지배 영역(Pe~40), 문헌 관례가 τ_B 기준
══════════════════════════════════════════════════════════════════════

INPUT (dimensional, SI)
  d   = 1.20 µm      [sketch:note1.png#ann3, tier 1]
  T   = 298 K        [user, tier 0]
  η   = 1.00 mPa·s   [handbook:water@25C, tier 0]
  ρ_p = 2000 kg/m³   [kb:silica, tier 0]
  φ   = 0.600        [user, tier 0]
  N   = 4000         [user, tier 0]
  v₀  = 11.6 µm/s    [kb:10.1103/PhysRevLett.110.238301#fig2, tier 2]  ⚠ 미검증

DERIVED (dimensional)
  γ   = 1.131e-8 kg/s        D_t = 0.364 µm²/s       D_r = 0.758 1/s
  m   = 1.81e-15 kg          L   = 86.8 µm           ℓ_p = 15.4 µm

SCALE LEDGER
  lengths   d=1.20µm  <  a_mean=1.37µm  <  r_c=1.35µm  <  ℓ_p=15.4µm  <  L=86.8µm
  times     τ_p=1.6e-7s  <  dt=1.98e-4s  <  τ_v=0.103s  <  τ_r=1.32s
            <  τ_B=3.96s  <  T_obs=3960s
  energies  kT=4.11e-21 J  =  ε=4.11e-21 J  <  f_a·d=1.58e-19 J

DIMENSIONLESS GROUPS
  φ      = 0.600      점유부피비
  Pe     = 38.4       τ_B/τ_v      이류 vs 확산
  D_r*   = 3.00       τ_B/τ_r      회전 vs 병진  (Stokes 예측 3.00 ✓)
  ℓ_p/d  = 12.8       Pe/D_r*      지속길이 = 입자 12.8개분
  ε*     = 1.00       ε/kT         WCA
  dt*    = 5.0e-5     dt/τ_B
  L/d    = 72.4       박스 = 입자 72개분
  T*     = 1000       관측창 = 1000 τ_B
  St     = 4.0e-8     τ_p/τ_B
  Re     = 1.4e-5     유체 관성

SCALE SEPARATION CHECKS                        value      limit    margin
  ✓ 관성 무시        τ_p/dt        =  8.1e-4   ≤ 1e-2     12.4×
  ✓ 이류 해상        dt/τ_v        =  1.9e-3   ≤ 1e-2      5.2×   ← 가장 타이트
  ✓ 회전 해상        dt/τ_r        =  1.5e-4   ≤ 1e-2     66×
  ✓ 상호작용 해상    dt/τ_int      =  5.0e-5   ≤ 1e-2    200×
  ✓ 저 Reynolds      Re            =  1.4e-5   ≤ 1e-3     71×
  ✓ 컷오프           r_c/L         =  0.016    <  0.5      31×
  ⚠ 유한크기         ℓ_p/L         =  0.178    ≤ 0.25      1.4×   ← 여유 부족
  ✓ 관측 충분        T_obs/τ_B     =  1000     ≥ 100       10×
  ℹ 조밀도           a_mean/d      =  1.14              → 조밀계 (접촉 지배)

VERDICT: PASS (경고 2건)
  ⚠ ℓ_p/L = 0.178 — 지속길이가 박스의 1/5.6. 한계(1/4)까지 1.4배 여유뿐.
     → N을 16000으로 늘리면 L=174µm, ℓ_p/L=0.089 (2.8배 여유). 벽시계 시간은 4배.
  ⚠ v₀는 tier 2 (미검증 문헌값) — 사람 확인 필요.
  ℹ Pe를 상향할 계획이면 이류 해상(5.2배 여유)이 먼저 걸립니다. Pe=200이면 dt*를 1e-5로.

RESOURCE
  production 2.0e7 steps = 1000 τ_B = 66 min (물리) ≈ 47 min (벽시계)
  raw data: A 0.5GB + B 0.2GB + D 0.05GB = 0.75 GB
```

이 리포트가 **사람 확인 #3**의 대상입니다. 무차원수 숫자만 보여주는 게 아니라
스케일 원장과 분리 검사를 함께 제시해, 사람이 "물리적으로 말이 되는가"를 판단할 수 있게 합니다.

### 6.6 ⑤ 역변환과 역구성

**역변환 (결과 → 물리 단위)** — 항상 수행:
```python
D_eff = D_eff_star * (sigma**2 / tau_B)      # 1.83 → 0.666 µm²/s
P     = P_star * (kT / sigma**dim)           # 무차원 압력 → Pa
t     = t_star * tau_B                       # 스텝 → 초
```
`observables.parquet`은 `y`(무차원)와 `y_physical`+`y_unit_si`를 항상 쌍으로 저장합니다.

**역구성 (무차원 → 물리계)** — 무차원 값으로 시작하고 싶을 때의 **유일한** 경로:
```python
# "Pe=40, φ=0.6인 2D ABP를 돌려줘"  ← 물리계가 없음
system = PhysicalSystem.from_dimensionless(
    groups={"Pe": 40, "phi": 0.6, "D_r_star": 3.0},
    anchors={                                    # ⭐️ 앵커를 반드시 지정해야 함
        "d": 1.0 * ureg.micrometer,              #    (아니면 물리계가 결정 안 됨)
        "T": 298 * ureg.kelvin,
        "eta": 1.0 * ureg.mPa * ureg.s,
    },
    note="문헌 관례 좌표. 앵커는 전형적 실리카 콜로이드 기준.",
)
```
앵커 없이는 물리계가 결정되지 않으므로 **에러**입니다. 앵커 기본값은 `materials.py`가 제공하되,
"임의로 정한 앵커"임을 `Provenanced.confidence = 3`으로 표시하고 리포트에 명기합니다.
이렇게 하면 무차원 작업을 하면서도 스케일 분리 검사(§6.4)를 그대로 받을 수 있습니다.

---

## 7. 지식베이스 (L1 / L10) — 상세

### 7.0 현재 상태 — `record.json` 부터 (2026-08-03) ⭐️

**SQLite KB는 아직 만들지 않았습니다.** 런 1건·문헌 0편에서는 조기 추상화입니다.
대신 **데이터 형식만 먼저 고정**해 흘려보냅니다:

```
cases/*.py            → runs/<id>/metrics.json   (기계가 읽는 결과)
tools/postmortem.py   → runs/<id>/record.json    (자동 진단 + tier3 KB 엔트리)
tools/kb.py           → list | query | lessons   (glob + 필터)
```

이렇게 하면 Phase 5가 도착했을 때 **환류할 히스토리가 이미 존재**합니다.
런이 100건을 넘거나 문헌이 들어오면 그때 §7.1의 SQLite로 옮깁니다 (형식은 그대로).

`record.json` 은 §5.2 `KnowledgeEntry` 의 축소판입니다 — `system_tags`,
`dimensionless`, `observables`(측정 vs 해석해), `outcome`, `failure_modes`,
`not_verified`, `lessons`(tier 3).

**`not_verified` 필드가 중요합니다.** "직접 확인하지 않은 것"을 명시적으로 남깁니다
(예: `dt_convergence_direct`). 성공 기록이 과잉 신뢰로 이어지지 않게 하는 장치입니다.

### 7.1 저장소 선택 (Phase 5에서)
**SQLite + FTS5**로 시작합니다. 이유:
- 로컬 단일 머신, 항목 수 수백~수천 규모 → 벡터 DB는 과잉
- 무차원수 좌표 검색(`Pe BETWEEN 30 AND 50 AND phi BETWEEN 0.55 AND 0.65`)이 핵심인데, 이건 **구조화된 SQL 쿼리**가 임베딩보다 정확
- FTS5로 전문 검색 병행
- 항목이 수천 개를 넘으면 그때 임베딩 컬럼 추가 (스키마는 미리 열어둠)

### 7.2 세 가지 검색 모드
| 모드 | 쿼리 | 용도 |
|---|---|---|
| **좌표 검색** | 무차원수 근방 | "Pe≈40, φ≈0.6 근처 선행 연구" |
| **태그 검색** | `system_tags` 매칭 | "2D ABP WCA" |
| **전문 검색** | FTS5 | "MIPS onset" |

세 결과를 합치고 confidence tier로 정렬해 LLM 컨텍스트에 주입합니다.

### 7.3 논문 증류 파이프라인
```
PDF → (a) 텍스트/표 추출  → LLM 증류 → KnowledgeEntry 초안
   → (b) 그림 추출        → vision  ↗       ↓
                                     [사람 검증] → confidence 2 → 1
```
- 증류 프롬프트는 **버전 관리**합니다 (`distill_prompt_v3`). 프롬프트를 바꾸면 재증류 가능.
- LLM에게 **"논문에 없으면 null"** 을 강하게 지시하고, `locator`(그림/표 번호)를 반드시 채우게 합니다.
- 검증 UI: 항목별로 원문 스니펫을 옆에 띄우고 ✓/✗/수정.

### 7.4 자체 런 환류
`PostMortem.lessons`가 `Claim`으로 변환되어 tier 3 엔트리로 들어갑니다.
```
"dt*=1e-4 at Pe>80 diverged (3/3 runs). Safe: dt* ≤ 2e-5."
   kind = "pitfall",  coords = {"Pe": 80, "dt_star": 1e-4}
```
→ 다음 검증 시 **자동 경고 규칙**으로 승격 (동일 실패 3회 이상이면 하드 에러로 승격 제안).

---

## 8. 인테이크 (L0) — 스케치·메모·그림 해석

### 8.1 입력 형태
| 형태 | 처리 |
|---|---|
| 손그림 스케치 (사진/스캔) | Claude vision → `Observation` |
| 손메모 텍스트 | vision (필기) 또는 plain text |
| 화이트보드 사진 | vision, 여러 장 지원 |
| 논문 그림 캡처 | vision + 문맥 질의 |
| 논문 PDF | 텍스트 + 그림 분리 후 증류 (7.3) |
| 자연어 | 그대로 |

기술: `client.messages.create()`에 `{"type": "image", "source": {"type": "base64", ...}}` 블록.
여러 장은 한 메시지에 여러 이미지 블록으로. PDF는 `{"type": "document"}` 블록(base64 PDF) 또는 Files API.

### 8.2 해석 프로토콜 (환각 억제)
LLM에게 요구하는 것:
1. **먼저 그대로 옮겨 적기** (`raw_transcription`) — 해석 전에 전사
2. 그 다음 구조화 (`entities`, `stated_quantities`)
3. **`ambiguities`와 `unread_regions`를 반드시 채우기** — 빈 리스트면 의심
4. 스케치에 **없는** 값은 절대 지어내지 않기 → 결측은 `null`로 두고, L1(KB)에서 근거와 함께 채우기

### 8.3 사람 확인 #1 — 필드별 승인
자유서술 요약을 보여주고 "맞나요?"라고 묻지 않습니다. **항목별 폼**으로 제시:
```
[1] system_guess: "2D active colloidal monolayer"        [✓ 승인] [✗] [수정: ___]
[2] particle diameter ≈ 1 µm    (출처: 그림 좌상단 "d~1um")  [✓] [✗] [수정: ___]
[3] Pe range 20–100             (출처: 그림 우측 화살표 라벨)  [✓] [✗] [수정: ___]
[!] ambiguity: 화살표가 자기추진 방향인지 전단 방향인지 불명    [해석: ___]
[!] unread: 우하단 손글씨 판독 불가                          [입력: ___]
```
승인 결과는 승인 원장에 기록됩니다.

---

## 9. 로우데이터 저장 전략 (L7) ⭐️

"입자 움직임과 힘을 저장하고 나중에 꺼내 쓴다" — 그대로 하면 디스크가 터집니다.
N=10⁴, 10⁸ steps에서 매 스텝 위치+힘 저장 = **수십 TB**.
→ **계층적 저장 정책**으로 해결합니다.

### 9.1 저장 계층

| Tier | 내용 | 빈도 | 크기 (N=10⁴, 10⁷ steps 기준) | 기본 | **구현** |
|---|---|---|---|---|---|
| **A** | 위치 + 방향 (전체 입자) | `10⁴` steps마다 | ~0.5 GB | 항상 ON | ✅ **1-A에서 사용** |
| **B** | 위치 + 방향 + **속도/힘/토크** | `10⁵` steps마다 | ~0.2 GB | 기본 ON | 실동작만 검증 |
| **C** | 고빈도 버스트 (짧은 구간) | `10` steps × 1000 프레임 | ~2.5 GB | 요청 시 | 실동작만 검증 |
| **D** | 추적 입자 (subset, 고빈도) | 100개, `10` steps마다 | ~0.05 GB | 기본 ON | 실동작만 검증 |
| **L** | 전역 스칼라 로그 | `10³` steps마다 | < 10 MB | 항상 ON | 실동작만 검증 |

> 5계층 **전부 실동작 확인**됨 (Phase 0, `scratch/smoke.py` 15/15). 1-A는 입자 간
> 상호작용이 없어 Tier A만 필요했습니다. 힘이 의미를 갖는 첫 케이스(1-B soft-r3)에서
> Tier B를 켭니다.

**Tier C가 핵심 아이디어**: `hoomd.write.Burst`가 슬라이딩 윈도우로 최근 N 프레임을 메모리에 들고 있다가 `.dump()` 호출 시에만 디스크에 씁니다. → "관심 사건 발생 시점 주변만 고해상도로" 저장할 수 있습니다.
예: 클러스터가 형성되는 순간, 에너지가 급변하는 순간.

**Tier D**: 소수 입자를 끝까지 고빈도로 추적 → MSD, 속도 자기상관, 단일입자 궤적 통계에 충분. 전체 입자를 고빈도로 저장할 필요가 거의 없습니다.

### 9.2 HOOMD 구현 매핑
| Tier | 구현 |
|---|---|
| A | `hoomd.write.GSD(trigger=Periodic(1e4), dynamic=['property'])` |
| B | `GSD(trigger=Periodic(1e5), dynamic=['property','momentum'])` + `logger`에 per-particle 힘 추가 → GSD `log/` 네임스페이스 |
| C | `hoomd.write.Burst(max_burst_size=1000, trigger=Periodic(10))` + 커스텀 Action이 조건 만족 시 `.dump()` |
| D | `GSD(filter=hoomd.filter.Tags([...]), trigger=Periodic(10))` |
| L | `hoomd.logging.Logger` + `hoomd.write.HDF5Log` |

per-particle 힘 로깅:
```python
logger = hoomd.logging.Logger(categories=['particle'])
logger.add(lj, quantities=['forces', 'energies'])
logger.add(active, quantities=['forces', 'torques'])
gsd_b = hoomd.write.GSD(filename='traj_forces.gsd', trigger=Periodic(int(1e5)),
                        mode='xb', logger=logger, dynamic=['property','momentum'])
```

### 9.3 정책은 스펙에 명시
```python
class RawDataPolicy(BaseModel):
    tier_a_every: int = 10_000
    tier_b_every: int | None = 100_000       # None이면 끔
    tier_c: BurstPolicy | None = None        # 조건부 고빈도
    tier_d_n_tracers: int = 100
    tier_d_every: int = 10
    log_every: int = 1_000
    estimated_bytes: int                     # 검증기가 계산, 사람에게 보여줌
```
검증기가 예상 용량을 계산해 **5 GB 초과 시 사람 승인**을 요구합니다.

### 9.4 온디맨드 재분석
- **프레임 인덱스**: 런 종료 시 `index.parquet` 생성 (프레임 번호 ↔ step ↔ 파일 오프셋 ↔ 시간)
- **지연 로더**: `loader.frames(run_id, tier='A', steps=slice(1e6, 2e6))` → 필요한 프레임만 읽음
- **관측량 캐시**: 계산한 관측량은 `observables.parquet`에 저장, 같은 요청 재계산 안 함
- **에이전트 툴**: `query_raw_data(run_id, what, when, who)` — 예: "3e6~4e6 스텝 구간에서 클러스터 안쪽 입자들의 힘 분포"

### 9.5 보관 정책
- 완료 런: 90일 후 Tier C/D 삭제 (A/B/L만 보존)
- 실패 런: 30일 후 궤적 삭제 (`PostMortem`과 로그는 영구 보존)
- `bdbot gc --dry-run`으로 확인 후 정리

---

## 10. 사후분석 & 학습 루프 (L10) ⭐️

> **구현 상태 (2026-08-03)**: [`tools/postmortem.py`](../../tools/postmortem.py)가 §10.1 분류체계와
> §10.2 자동 진단을 구현하고 `record.json`(tier 3 KB 엔트리)을 방출합니다.
> [`tools/kb.py`](../../tools/kb.py)로 질의합니다. **아직 없는 것**: LLM 서술(§10.3),
> 규칙 승격(§10.4). 런 1건이라 승격을 논할 표본이 없습니다.
>
> 실제 산출 예시는 `runs/trap-2d-5um__70b9394e7310/record.json`.

### 10.1 실패 분류 체계 (taxonomy)
자유서술로 두면 축적이 안 됩니다. **미리 정의된 분류**에 매핑합니다.

| 코드 | 의미 | 자동 감지 지표 |
|---|---|---|
| `NUM_DIVERGE` | 수치 발산 | NaN/Inf, PE가 초기 대비 >100× |
| `NUM_DRIFT` | 에너지/온도 드리프트 | 정상상태 후 PE 추세 유의 |
| `EQ_INSUFFICIENT` | 평형화 부족 | 전반부/후반부 블록 평균 불일치 |
| `STAT_INSUFFICIENT` | 통계 부족 | 관측량 오차막대 > 임계 |
| `FINITE_SIZE` | 유한크기 효과 | 상관길이 > L/4 |
| `WRONG_REGIME` | 목표 현상 미발현 | 목표 관측량이 기대 범위 밖 |
| `RESOURCE` | 시간/디스크 초과 | 러너 기록 |
| `SPEC_ERROR` | 스펙 자체가 틀림 | 사람 판정 |
| `SUCCESS` | 성공 | 위 전부 통과 |

### 10.2 자동 진단 (LLM 없이)
런 종료 시 항상 실행:
- **평형 판정**: 궤적을 5블록으로 나눠 PE·압력 블록 평균의 정상성 검정
- **에너지 드리프트**: 선형 회귀 기울기 유의성
- **유한크기**: g(r) 상관길이 vs L/4
- **통계 충분성**: 블록 평균법으로 오차막대 산출, 상대오차 확인
- **목표 달성**: `target_observables`가 기대 범위 안인지

### 10.3 LLM 서술 + 교훈 추출
자동 진단 결과 + 스펙 + 무차원 좌표를 LLM에 주고:
1. 자연어 서술 (`narrative`)
2. **재사용 가능한 교훈** (`lessons: list[Claim]`) — 무차원 좌표를 반드시 포함
3. 다음 시도 제안 (파라미터 수정안)

### 10.4 규칙 승격
같은 좌표 영역에서 같은 실패가 반복되면:
```
3회 반복 → 검증기 경고 규칙 자동 생성 제안 (사람 승인)
5회 반복 → 하드 에러 승격 제안
```
`validate/learned_rules.py`에 규칙이 코드로 생성되고, 사람이 리뷰 후 머지합니다.
**LLM이 검증 규칙을 직접 코드로 쓰지 않습니다.** 템플릿에 파라미터만 채웁니다.

### 10.5 성공 요인도 기록
실패만 기록하면 편향됩니다. 성공 런도:
```
"φ=0.6, Pe=40, dt*=5e-5, N=4000, prod=1000τ_B → MIPS 관측, 통계 충분, 47분"
   kind = "parameter", tier = 3
```
→ 유사 요청 시 **검증된 출발점**으로 제안됩니다.

---

## 11. HOOMD API 매핑 (문서 v7.1.1 / **설치본 7.1.0** — 조사한 API는 동일)

| SimSpec | HOOMD v7.1.1 |
|---|---|
| device | `hoomd.device.CPU()` |
| seed | `hoomd.Simulation(device=dev, seed=spec.run.seed)` |
| box + init | `gsd.hoomd.Frame` 구성 → `sim.create_state_from_snapshot(...)` |
| pair `wca` | `md.pair.LJ(nlist=cell, default_r_cut=2**(1/6)*σ, mode='shift')` |
| pair `lj` | `md.pair.LJ(nlist=cell, default_r_cut=r_cut, mode='shift')` |
| pair `yukawa` | `md.pair.Yukawa(nlist=cell, default_r_cut=r_cut)` |
| nlist | `md.nlist.Cell(buffer=0.4)` |
| BD 적분 | `md.methods.Brownian(filter=All(), kT=kT, default_gamma=γ)` |
| 적분기 | `md.Integrator(dt=dt, methods=[bd], forces=[...])` |
| ABP 힘 | `md.force.Active(filter=All())`; `active.active_force['A'] = (f_a, 0, 0)` |
| ABP 회전확산 | `active.create_diffusion_updater(trigger, rotational_diffusion=D_r)` |
| 궤적 | `hoomd.write.GSD(...)` |
| 고빈도 버스트 | `hoomd.write.Burst(...)` + `.dump()` |
| 로그 | `hoomd.logging.Logger` + `hoomd.write.HDF5Log` |
| 열역학량 | `md.compute.ThermodynamicQuantities(filter=All())` |
| 재시작 | `GSD(mode='wb', truncate=True, dynamic=[...])` |
| 겹침 제거 | `md.minimize.FIRE(...)` 짧게 → 본 적분기로 교체 |

### ⚠️ 반드시 지킬 함정 목록
1. **WCA 전용 클래스 없음** → `md.pair.LJ` + `r_cut=2^(1/6)σ` + `mode='shift'`. (`ForceShiftedLJ`는 WCA가 아님)
2. **Brownian은 O(δt) 오차** — 문서가 명시. Langevin보다 훨씬 작은 dt 필요. 기본 `dt* = 1e-4`, 힘 강하면 더 축소.
3. **ABP 회전은 적분기가 아니라 updater로**:
   ```python
   integrator.integrate_rotational_dof = False   # ← 필수
   sim.operations.updaters.append(
       active.create_diffusion_updater(trigger, rotational_diffusion=D_r))
   ```
   `True`로 두면 관성 회전이 섞여 ABP가 아니게 됨.
4. **2D**: `Lz=0`, `dimensions=2`. orientation quaternion이 z축 회전만 하도록 초기화.
5. **BD는 과감쇠** — 속도에 물리적 의미 없음. `thermalize_particle_momenta` 불필요. 속도 기반 MSD 금지.
6. **`r_cut < L/2`** — 아니면 minimum image 위반.
7. ⭐️ **외부 힘 + 주기경계: 최소 이미지를 적용하지 않으면 조용히 틀린다.** 실측으로 확인한 함정.
   고정 앵커를 향한 트랩에서 `d = pos - anchor`를 그대로 쓰면, 입자가 박스를 넘어 wrap되는 순간
   거리가 L만큼 점프해 거대한 **반대 방향** 복원력을 받습니다. **터지지 않고 조용히 틀립니다** —
   트랩이 강하면 정확한 값이 나오고 약할수록 오차가 커집니다 (k=10 +0.2% ✓ / k=2 **+1856%** ✗).
   ```python
   d = pos - anchors[tags]
   d -= L * np.round(d / L)      # ← 이 한 줄. external.* 모듈 전체에 해당
   ```
   → `PhysicsModule.periodic_safe()`를 선언 필수 필드로 두고, 검증기가 확인합니다.
8. **2D에서 최소 이미지의 z성분 NaN** — 비주기 축의 박스 길이를 `inf`로 두면
   `inf * round(0/inf) = nan`. 힘 배열에 NaN이 들어가 런타임 가드가 오탐합니다. 주기 축만 마스크.
9. **`write.Burst`는 새 파일에 `write_at_start=True` 필요** — 없으면
   `RuntimeError: Must set write_at_start to write to a new file.`

> 함정 7~9는 실측에서 발견했습니다. 전체 조사·검증 결과: **[`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md)**

---

## 12. 검증 레이어 (L4)

### 12.0 구조 검사 (불변식) ⭐️
- `SimSpec.derived_from`이 비어 있으면 **거부** — 물리계 없는 무차원 스펙은 실행 불가 (원칙 3)
- `ScaleLedger`가 없거나 기준 스케일에 `rationale`이 비어 있으면 **거부**
- `PhysicalSystem`의 모든 필드가 `pint` 차원 검사를 통과하는지 확인 (단위 불일치 = 즉시 거부)

### 12.1 스케일 분리 하드 게이트 ⭐️
§6.4의 ❌ 항목이 검증기의 하드 에러입니다. **전부 `10⁻²` 기준으로 통일**:
`τ_p/τ_dyn ≤ 10⁻²`(모델) · `dt/τ_v ≤ 10⁻²` · `dt/τ_r ≤ 10⁻²` · `dt/τ_int ≤ 10⁻²` · `Re ≤ 10⁻³` · `r_c < L/2`

에러 메시지는 항상 **한계까지의 여유(margin)와 구체적 수정안**을 포함합니다:
```
❌ 이류 해상 위반: dt/τ_v = 3.2e-2 (한계 1e-2, 3.2배 초과)
   원인: Pe=200에서 τ_v=0.0198s인데 dt=6.3e-4s
   수정: dt* 를 5e-5 → 1.5e-5 로 낮추세요 (벽시계 시간 3.3배 증가)
```

### 12.1b 여유 경고
하드 검사 중 한계의 1/5 이내면 경고 — 파라미터를 조금만 올려도 위반하게 되므로,
스윕을 계획 중이라면 스윕 전 구간에 대해 미리 검사합니다.

### 12.2 스펙 범위 검사
`φ` 상한(2D 0.9, 3D 0.74) · `100 ≤ N ≤ 10⁵` · Pe/force 중복 지정 금지 · 쓰기 빈도 하한 · 디스크 상한

### 12.2b 물리 일관성 경고
`D_r*`가 Stokes 예측(구=3)에서 크게 벗어나면 경고 — 의도적이면 사유 기록 요구.
`a_mean/d`가 1에 가까우면 조밀계, ≫1이면 희박계로 분류해 이후 해석에 태그.

### 12.3 KB 기반 학습 경고 ⭐️
```
⚠ 과거 경고: Pe=90, dt*=1e-4 영역에서 3회 발산 (run_ids: a3f9…, b71c…, e024…)
             권장: dt* ≤ 2e-5
⚠ 문헌 대조: Redner et al. 2013은 이 영역에서 N≥10⁴를 사용 (현재 N=2000)
             → 유한크기 효과 가능
```

### 12.4 출처 신뢰도 검사 ⭐️
tier 2 이하 값만으로 구성된 스펙 → 사람 승인 강제.

### 12.5 런타임 가드
`10⁴` 스텝마다 NaN/Inf, PE 폭발 검사 → 중단 + `status: "diverged"` + 마지막 스냅샷 보존.

---

## 13. 분석 & 관측량 (L8)

### 13.1 관측량 목록
**구조**: `rdf` g(r) · `sq` S(q) · `psi6` 육방오더(2D) · `voronoi_density`
**동역학**: `msd` · `diffusion_coefficient` · `fskt` · `vacf`(Tier D 활용)
**액티브**: `cluster_size` · `local_density_hist`(→ MIPS 이봉성) · `polar_order` · `giant_number_fluctuations`
**힘 통계** ⭐️ (Tier B/C 활용): `force_distribution` · `stress_per_particle` · `virial_pressure_decomposition`
**열역학**: `potential_energy` · `pressure` · `pressure_tensor`

### 13.1b 산출물 형식 (2026-08-03 확정)

| 파일 | 내용 | 상태 |
|---|---|---|
| `report.txt` | DimensionlessReport (스케일 표 + 분리 검사 + 여유) | ✅ |
| `result.txt` | + 물리 단위 역변환 후 해석해 대조 | ✅ |
| `metrics.json` | **기계가 읽는** 결과 — 좌표·관측량·검사 | ✅ |
| `record.json` | 자동 진단 + `outcome` + `lessons` + `not_verified` | ✅ |
| `observables.npz` | P(x)·C(t)·PSD 원자료 | ✅ |
| `observables.png` | 4패널 플롯 (측정 vs 해석해) | ✅ |
| `traj_A.gsd` | 원시 궤적 (Tier A) | ✅ |
| `observables.parquet` | long-format 테이블 (§13.2) | ❌ 미구현 — npz로 대체 중 |

### 13.2 출력 스키마 (목표)
```
run_id | observable | x | y | y_err | x_unit_reduced | y_unit_reduced
       | x_physical | y_physical | x_unit_si | y_unit_si | metadata(json)
```
⭐️ **무차원 값과 물리 단위 값을 동시에** 저장 (원칙 3).

### 13.3 표준 플롯
1. 에너지/압력 vs 시간 (**항상**, 평형 판정 근거)
2. g(r) / MSD(log-log, 기울기 1 가이드) / S(q)
3. ABP: 국소밀도 히스토그램, 클러스터 분포, 최종 스냅샷
4. ⭐️ **문헌 대조 플롯**: 같은 무차원 좌표에서 KB의 문헌 값을 점으로 오버레이

### 13.3b 결과 확인은 **그래프 + 애니메이션이 기본값** ⭐️⭐️ (2026-08-05)

> 결과를 말과 표로만 보고하지 않습니다. 사람이 결과를 확인하는 통상적인 방법이
> **그림과 영상**이므로, 리포트에는 둘 다 만들어 첨부하는 것이 기본 동작입니다.

**그래프 — "무엇이 틀렸는지"가 눈에 보여야 합니다.** 측정과 예측을 같은 축에 놓고,
검사 문턱·해석해·문헌값을 선으로 겹칩니다. 본보기: [`scratch/viz_chain_bend.py`](../../verify/viz_chain_bend.py)
— 6패널이 `chain-bend` 스펙 오류 5건에 하나씩 대응하고, 패널마다 "틀린 값"과 "맞는 값"이
같은 그림 안에 있습니다 (예: SNR 검사가 그은 ω-무관 직선 vs 실제로 떨어지는 곡선).

**애니메이션 — 계가 실제로 어떻게 움직이는지.** 특히 `kT=0`(결정론, 모드 형태)과
`kT>0`(열요동)을 **나란히** 두면 SNR 문제가 숫자보다 빨리 보입니다. `chain-bend` 에서는
위아래 두 패널로 "구동 응답이 ±ℓ_k 띠 안에 묻힌다"가 한눈에 드러났습니다.
애니메이션은 **큰 `dt`로 싸게** 만들어도 됩니다 (안정 한계의 10% 정도). 단 **생산 측정이
아님을 반드시 명시**하고, 정적 평형량(ℓ_k 등)이 그 `dt`에서도 옳은지 확인합니다.

**⚠️ 그래프 라벨은 한글이 아니라 영어로 쓴다** (실측). matplotlib 기본 `DejaVu Sans`에
한글 글리프가 없어 라벨이 전부 `□`가 되는데, 한글을 갖춘 폰트로 바꾸면 이번엔 기호가
빠집니다 — `AppleGothic` · `Apple SD Gothic Neo` · `NanumGothic` 셋 다 `−`(U+2212)와
`ŷ`(U+0177)가 **없습니다**. 한글과 기호를 모두 갖춘 것은 `Arial Unicode MS` 뿐이라
폰트 전환으로는 완전히 못 고칩니다 (`axes.unicode_minus = False`도 눈금 라벨만 고치고
**문자열 리터럴의 `−`는 그대로 깨집니다**). → **폰트를 맞추려 하지 말고 축·범례·타이틀·
주석을 처음부터 전부 영어로 씁니다.** 생성 후 `missing from font` 경고가 **0건**인지
반드시 확인하세요 — 경고를 무시하면 사람이 읽을 수 없는 그림을 결과라고 내놓게 됩니다.
폰트가 문제없더라도 이 규칙은 유지합니다 — 이 분야 물리·통계 용어(`de-correlation
time`·`shear thinning`·`yield force` 등)는 한국어 표현이 어색하거나 자리잡지 않은
경우가 많아, 영어 쪽이 뜻을 더 정확하고 간결하게 전달합니다.

ffmpeg는 이 환경에 **없습니다** → `PillowWriter`로 GIF를 씁니다 (`FFMpegWriter.isAvailable()` = False).

---

## 14. 실행 / 잡 관리 (L6)

- **프로세스 풀 8개** (10 cores − 2). 워커는 완전 독립 프로세스.
- **SQLite 잡 큐** (`runs.db`): run_id, status, pid, progress, error
- **체크포인트**: `restart.gsd` truncate 덮어쓰기, `bdbot resume <run_id>`
- **상태 머신**: `queued → running → {completed | failed | diverged | interrupted}`
- **스윕**: `SweepSpec.axes` 전개 → run_id 해시로 자동 중복 제거

---

## 15. Claude Code 통합 레이어 ⭐️ (v0.3에서 전면 변경)

**Claude Code가 에이전트 런타임입니다.** 우리는 에이전트를 만들지 않고, 그 아래 놓일
**결정론적 엔진(`bdbot`)** 과 **통합 표면(`.claude/`)** 만 만듭니다.

### 15.1 역할 분담 — 무엇이 사라지는가

| 계획 요소 | v0.2 (자체 구현) | v0.3 (Claude Code) |
|---|---|---|
| 에이전트 루프 | `tool_runner()` + `@beta_tool` | **내장** |
| 모델 호출·캐싱·스트리밍 | anthropic SDK 직접 | **내장** |
| 스케치 읽기 (vision) | base64 image 블록 | **Read 툴 — 무료** |
| 논문 PDF 읽기 | document 블록 / Files API | **Read 툴 (`pages` 인자) — 무료** |
| 대화 이력 관리 | 직접 구현 | **내장** |
| 툴 10종 | `@beta_tool` 함수 | **`bdbot` CLI (Bash)** |
| 도메인 지식 주입 | 시스템 프롬프트 + `cache_control` | **Skill + `CLAUDE.md`** |
| 승인 게이트 | 툴 내부 분기 | **Hooks + AskUserQuestion + 권한** |
| 워크플로 규범 | 프롬프트 지시 | **Slash commands** |
| 전문 역할 분리 | (없음) | **Subagents (컨텍스트 격리)** |
| 장기 실행 폴링 | 직접 구현 | **Bash `run_in_background`** |

**삭제**: `bdbot/agent/` 대부분 · `bdbot chat` · `anthropic` 필수 의존성
**남는 anthropic 용도**: 논문 **배치** 증류만 (§15.9)

### 15.2 통합 방식: **CLI + Skill** (MCP 아님)

| 방식 | 장점 | 단점 | 판정 |
|---|---|---|---|
| **CLI + Skill** | 서버 불필요 · 파일 기반이라 검사·diff·버전관리 가능 · `Bash(bdbot:*)` 허용목록으로 권한 제어 · 사람도 똑같이 씀 | 스키마가 모델에 직접 노출 안 됨 (→ Skill이 보완) | ✅ **채택** |
| MCP 서버 | 타입 스키마 노출 · 구조화 반환 | 서버 프로세스 유지 · 로컬 파이썬 패키지엔 과잉 · 디버깅 어려움 | 나중에 필요하면 CLI 위에 얇게 |
| 파일만 + 프롬프트 | 최소 | 검증·게이트를 강제할 수 없음 | ❌ |

**핵심 근거**: 이 파이프라인의 산출물은 전부 파일입니다(`system.yaml` → `spec.json` → `traj.gsd` → `observables.parquet`). Claude Code는 파일 작업이 본업이므로,
`Write(system.yaml)` → `Bash(bdbot nondim system.yaml)` → `Read(report.txt)` 가 가장 자연스럽습니다.
게다가 중간 산출물이 전부 사람이 열어볼 수 있는 파일로 남습니다 — 원칙 2(출처 추적)와 잘 맞습니다.

### 15.3 `.claude/` 구조

```
.claude/
├── settings.json                    # 권한 허용목록 + 훅 등록
├── skills/
│   ├── bd-physics/SKILL.md          # 단위계·스케일 원장·무차원화 규약 (§2, §6)
│   ├── bd-hoomd/SKILL.md            # HOOMD v7.1.1 API 매핑 + 함정 6종 (§11)
│   ├── bd-intake/SKILL.md           # 스케치 해석 프로토콜 (§8.2)
│   └── bd-distill/SKILL.md          # 논문 증류 프로토콜 (§7.3)
├── agents/
│   ├── bd-distiller.md              # 논문 1편 → KnowledgeEntry
│   ├── bd-analyst.md                # 런 결과 → PostMortem
│   └── bd-reviewer.md               # 스펙 적대적 물리 검토
├── commands/
│   ├── bd-intake.md                 # /bd-intake <folder>
│   ├── bd-spec.md                   # /bd-spec <observation>
│   ├── bd-run.md · bd-sweep.md
│   ├── bd-analyze.md · bd-postmortem.md
│   └── bd-distill.md
└── hooks/
    ├── guard_invariant.py           # ⭐️ derived_from 없는 스펙 실행 차단
    ├── guard_separation.py          # ⭐️ 분리 검사 미통과 차단
    ├── guard_cost.py                # 비용 초과 시 확인 요구
    └── log_approval.py              # 승인 원장 기록

CLAUDE.md                            # 항상 로드되는 불변식 (짧게)
```

### 15.4 `CLAUDE.md` — 항상 로드되는 불변식만

Skill은 필요할 때 로드되지만 `CLAUDE.md`는 **항상** 컨텍스트에 있습니다. 짧게 유지하고
**절대 규칙**만 넣습니다. 상세는 Skill로 미룹니다.

```markdown
# bdbot — Brownian dynamics simulation pipeline

## 절대 규칙
1. 차원이 먼저다. 모든 계는 SI 단위 PhysicalSystem으로 먼저 확정하고,
   스케일 원장을 거쳐 무차원화한다. 무차원 값으로 시작하는 경로는 없다.
2. SimSpec을 손으로 쓰지 않는다. `bdbot nondim`이 생성한 것만 실행한다.
3. 모든 파라미터에 출처(tier 0–3)를 붙인다. KB에 없으면 없다고 말한다.
   추정값은 추정이라고 표시한다. 절대 지어내지 않는다.
4. HOOMD 스크립트를 직접 작성하지 않는다. `bdbot` CLI만 사용한다.
5. 스케치를 읽을 때는 먼저 그대로 전사하고, 읽지 못한 부분과
   모호한 부분을 반드시 명시한다.

## 상세 지식
- 무차원화·스케일 원장 → skill `bd-physics`
- HOOMD API·함정 → skill `bd-hoomd`
- 워크플로 → `/bd-intake`, `/bd-spec`, `/bd-run`, `/bd-postmortem`
```

### 15.5 Hooks — 프롬프트가 아니라 **하네스**가 강제한다 ⭐️⭐️

이번 변경에서 가장 큰 소득입니다. 원칙 3 같은 하드 불변식을 프롬프트로 "부탁"하면
모델이 어길 수 있습니다. **PreToolUse 훅은 하네스가 실행**하므로 모델이 우회할 수 없습니다.

| 훅 | 이벤트 | 동작 |
|---|---|---|
| `guard_invariant` | PreToolUse `Bash(bdbot run*)` | 스펙에 `derived_from`/`scale_ledger` 없으면 **거부** |
| `guard_separation` | PreToolUse `Bash(bdbot run*)` | §6.4 하드 게이트 미통과면 **거부**, 수정안 반환 |
| `guard_cost` | PreToolUse `Bash(bdbot run*\|sweep*)` | 추정 비용 초과 시 사용자 확인 요구 |
| `log_approval` | PostToolUse | 승인/거부/수정을 원장에 기록 |
| `capture_postmortem` | Stop | 완료된 런에 사후분석 미실행분이 있으면 알림 |

훅이 거부할 때는 **구체적 수정안을 stderr로 반환**해 Claude Code가 스스로 고치게 합니다:
```
❌ blocked: spec has no `derived_from`.
   Create a PhysicalSystem first: bdbot init-system --template abp_2d > system.yaml
```

> 훅 설정은 `.claude/settings.json`에 들어갑니다. 이 파일은 `update-config` 스킬로 작성하는 게 안전합니다.

### 15.6 Subagents — 컨텍스트 격리가 목적

| 서브에이전트 | 역할 | 격리 이유 |
|---|---|---|
| `bd-distiller` | 논문 1편 → `KnowledgeEntry` | **논문 전문이 메인 대화를 오염시키지 않음.** 5편을 5개 동시 실행 가능 |
| `bd-analyst` | 런 결과 → `PostMortem` | 대용량 수치 출력이 메인 컨텍스트에 안 들어옴 |
| `bd-reviewer` | 스펙 적대적 검토 | **독립 관점** — 스펙을 만든 컨텍스트를 모른 채 물리적 타당성만 봄 |

`bd-reviewer`는 §17의 "에이전트 평가"를 자동화합니다. 스펙 제출 전에 "이 스펙의 물리적 문제를
반박해봐"를 독립 컨텍스트에서 돌리면, 자기가 만든 걸 자기가 검토하는 편향을 피할 수 있습니다.

### 15.7 Slash commands — 워크플로 고정

| 명령 | 하는 일 |
|---|---|
| `/bd-intake <folder>` | 폴더의 스케치/메모 전부 읽기 → `Observation` YAML → 사람 확인 |
| `/bd-spec <observation>` | KB 검색 → `PhysicalSystem` 제안(출처 부착) → 사람 확인 → `bdbot nondim` → 리포트 |
| `/bd-run <spec>` | 검증 → 비용 추정 → 승인 → 백그라운드 실행 |
| `/bd-sweep <sweep>` | 스윕 전개 → 전 구간 분리 검사 → 승인 → 큐 등록 |
| `/bd-analyze <run_id>` | 관측량 계산 → 플롯 → 문헌 대조 |
| `/bd-postmortem <run_id>` | `bd-analyst` 서브에이전트 → KB 환류 |
| `/bd-distill <pdf...>` | `bd-distiller` 서브에이전트 병렬 → 검증 대기열 |

워크플로를 프롬프트 규범이 아니라 **명령으로 고정**하면 순서가 지켜집니다.

### 15.8 승인 게이트 재매핑

| 게이트 | v0.2 구현 | v0.3 구현 |
|---|---|---|
| #1 Observation | 툴 내부 폼 | AskUserQuestion + YAML 파일 리뷰 |
| #2 PhysicalSystem | 필드별 폼 | YAML 파일 리뷰 (diff로 수정 확인) + AskUserQuestion |
| #3 SimSpec/무차원화 | 리포트 확인 | 리포트 파일 Read + AskUserQuestion |
| #4 비용 | 툴 내부 분기 | **`guard_cost` 훅** (하네스 강제) |
| #5 지식 tier 승격 | 툴 | AskUserQuestion |
| #6 규칙 승격 | 툴 | 사람이 PR처럼 리뷰 (파일 diff) |

초기 단계에서는 **plan mode**를 활용해, 실행 전에 전체 계획을 사람이 승인하는 흐름도 유효합니다.

### 15.9 `anthropic` SDK가 남는 유일한 자리 — 배치 증류

논문 50편을 대화형으로 증류하면 느리고 토큰이 비쌉니다. 대량 증류는 별도 스크립트로:
- **Message Batches API** (`client.messages.batches.create`) — 표준가 대비 **50% 할인**, 최대 10만 요청
- `client.messages.parse(output_format=KnowledgeEntry)` 로 Pydantic 검증된 결과 수령
- 모델 `claude-opus-5`

`pyproject.toml`에서 **선택 의존성**(`bdbot[distill]`)으로 둡니다. 소량 증류는 Claude Code에서
`/bd-distill`로 처리하므로 평소엔 필요 없습니다.

### 15.10 이 변경의 효과

| 항목 | 변화 |
|---|---|
| 개발 일정 | Phase 7+8이 4.5일 → **1.5일** (총 18~20일 → **15~17일**) |
| 코드량 | `bdbot/agent/` ~800줄 삭제, `.claude/` ~300줄 추가 |
| 불변식 강제력 | 프롬프트 부탁 → **하네스 강제** (더 강함) |
| vision/PDF | 직접 구현 → 무료 |
| 사람 개입 | 별도 UI 필요 → Claude Code 대화가 곧 UI |
| 위험 | Claude Code 세션 밖에서는 자동화가 안 됨 → **CLI가 단독 동작해야 하는 이유가 더 커짐** |

---

## 16. 개발 로드맵

각 Phase에 **완료 조건(DoD)** 명시. 통과해야 다음으로.

### ✅ Phase 0 — 범용 환경 & API 실증 **(완료 2026-08-03)**
```bash
conda env update -f environment.yml -n simulation_bot --prune
```
- [x] 환경 구축 — 기존 `simulation_bot` 확장 (hoomd **7.1.0** CPU, gsd 5.0.1, freud 3.5.0,
      numpy 2.5.1, + pyarrow/pydantic/pint/typer/rich 추가). `environment.yml`로 고정.
- [x] **HOOMD 능력 전수 조사** (`scratch/survey.py`) — 등방 페어 28종, 이방 페어 17종,
      본드/각도/이면각, 마찰 접촉, 다체, 메시, 장거리, manifold, 강체, HPMC, MPCD 확인
- [x] **15개 API 실동작 검증** (`scratch/smoke.py`) — **15/15 PASS**.
      로우데이터 5계층(§9) 전부 실동작 확인
- [x] **조화 트랩 골든 물리 검증** (`scratch/golden_trap.py`) — `⟨x²⟩=kT/k`를
      k 10배 범위에서 오차 0.6% 이내 재현, `⟨x²⟩·k` 변동계수 **0.28%**
- [x] **함정 3건 신규 발견** (§11 함정 7~9) — 특히 함정 7은 조용히 틀리는 유형
- [x] 결과 문서화 → **[`docs/hoomd_capabilities.md`](../../docs/hoomd_capabilities.md)**
- **DoD 달성**: `intake/`의 5개 케이스 전부 구현 가능함을 실측으로 확인.
  §11 매핑 표에 틀린 항목 없음(함정 3건 추가). BD 적분기 정확성 확인.

> **부수 소득**: `trap-2d-5um`이 골든 물리 테스트로 성립함을 미리 증명했습니다 → Phase 2 단축.

### Phase 1 — **케이스 주도 증분 구축** ⭐️ (v0.4에서 방식 변경)

> **왜 바꿨나**: 이전 계획은 `ScaleLedger` 전체·기준 전략 4종·무차원수 11종·분리 검사 8종을
> **먼저** 만드는 2일짜리 프레임워크 선행 설계였습니다. 그런데 무차원화는 계마다 다르고,
> 어떤 스케일이 실제로 필요한지는 케이스를 관통해봐야 압니다. 추측으로 일반화하면
> 안 쓰는 코드를 만들고 정작 필요한 건 빠집니다.
>
> **바꾼 방식**: 케이스 하나를 끝까지 관통 → 다음 케이스에서 공통점 발견 → **그때** 추상화.
> 각 단계는 반나절이고, 매 단계 끝에 **사람과 함께 물리를 검증**합니다.

#### 진행 규칙
1. 케이스 하나를 고른다 (쉬운 것부터)
2. 스케치 읽기 → 물리계 초안 → **함께 확인**
3. 그 계의 스케일 표를 손으로 작성 → **함께 확인** ← 여기서 무차원화 지식이 쌓임
4. 무차원화 → 실행 → 역변환 → 골든/문헌 대조 → **함께 결과 검토**
5. 배운 것을 기록 (스케일 표, 검사 임계, 함정) → 다음 케이스로
6. **케이스 2개를 마친 뒤에만** 공통 구조를 뽑아 추상화한다

> 3번이 핵심입니다. "이 계에 어떤 시간척도가 있고 무엇을 기준으로 삼을 것인가"는
> 코드가 아니라 물리 판단이고, 지식이 쌓이기 전에는 사람이 같이 봐야 합니다.

#### ✅ Phase 1-A · `trap-2d-5um` 관통 **(완료 2026-08-03)**
- [x] [`intake/trap-2d-5um/observation.yaml`](../../intake/trap-2d-5um/observation.yaml) — 스케치 전사·모호점 3건 해소
- [x] [`intake/trap-2d-5um/system.yaml`](../../intake/trap-2d-5um/system.yaml) — 물리계(SI), 필드별 tier
- [x] [`cases/trap_2d_5um.py`](../../cases/trap_2d_5um.py) — 스케일 표 → 분리 검사 → 무차원화 →
      실행 → 역변환 → 해석해 대조. **이 케이스 전용** (범용 프레임워크 아님)
- [x] 분리 검사 4종을 **모델/적분/기하/통계**로 분류, 여유(margin) 표시
- [x] `run_id` 콘텐츠 주소 지정 + 재실행 방지 (§14)
- **DoD 달성**: 5 µm·물·300K YAML → 물리 단위 결과 → **관측량 4종 전부 해석해와 일치**
  (`⟨x²⟩` +0.02%, `σ` +0.01%, `τ` −0.25%, `f_c` +1.17%, `S(0)` −0.33%)
- **발견**: 변위 자기상관에서 표본평균을 빼면 `τ`가 −7.75% 어긋남 → 안 빼면 −0.26%
  (bd-physics §5.1에 기록)
- **산출**: `runs/trap-2d-5um__70b9394e7310/` (리포트·궤적·관측량·플롯)

#### ✅ Phase 1-B · `soft-r3-2d-A-sweep` 관통 **(완료 2026-08-04)**
- [x] `intake/soft-r3-2d-A-sweep/{observation,system}.yaml` — 밀도·코어·N·앵커를 사람이 확정
- [x] [`cases/soft_r3_2d.py`](../../cases/soft_r3_2d.py) — `pair.Table`(r⁻³) + WCA 코어
- [x] **런 7건**: 진폭 스윕 4 (A=0.1·1·10·100) + 희박극한 검증 1 + 수렴 확인 2
- [x] 검증 5종 — 2입자 직접 대조 **0.000%** · 에너지 일관성 **+0.00~0.67%** ·
      육방 NN 거리 **+0.45%** · 희박극한(O(ρ) 포함) RMS **2.43%** · dt 절반 **−0.004%**
- [x] 1-A와의 대조 판정 완료 → skill `bd-physics` **§6.3**
- **DoD 달성**: 무엇이 두 번 나왔고(기준 스케일·물성 유도·`τ_p` 용도·`dt=10⁻²·γ/국소강성`·
  검사 분류·리포트·run_id·metrics 스키마) 무엇이 케이스마다 다른지(평형 지표·관측량·
  검증 전략·지배 시간척도) 표로 확정.
- **결정**: `ScaleLedger` 추상화를 1-C에서 도입한다. 단 §6.3의 "공통화하지 말 것"은 제외.

**1-B에서 새로 알아낸 것**
- **함정 3건** (skill `bd-hoomd` 10·11·12): `pair.Table` 격자가 `endpoint=False`
  (아니면 힘 −1.65%) · `r<r_min`에서 힘·에너지가 **0**(발산 퍼텐셜에서 위험) ·
  `seed`가 16비트로 잘려 **다른 seed가 같은 런**이 됨
- **`A` 단독은 제어 파라미터가 아니다.** `Γ = A(d/a_mean)³`가 구조를 정한다.
  밀도 없이는 진폭 스윕이 정의되지 않는다 (스케치에 밀도가 없었다)
- **`Γ_max = N^{3/2}u_c/8`** — 절대 kT 컷오프 기준과 최소이미지를 함께 풀면
  달성 가능한 결합세기에 `A`와 무관한 상한이 생긴다. 스케치의 N=100은 Γ≲1.25
- **2D `r⁻³`는 구조만 수렴한다.** `r_c` 5a→7a에서 ψ₆·NN·배위수는 0.15% 이내 불변,
  절대 `⟨U⟩/N`은 +7.5% (꼬리가 `1/r_c`로만 줄어서)
- **약결합이 수치적으로 더 비싸다.** Γ=0.03이 Γ=30보다 12.8배 (WCA 코어가 dt를 정함)
- **사후분석 도구의 두 가정이 1-A 전용이었다** — 아래 Phase 4 항목 참조

#### ✅ Phase 1-C · 추상화 **(완료 2026-08-04)**
1-A/1-B에서 **실제로 두 번 나온 것만** 공통화했습니다 → `bdbot/` 10개 모듈 (§4.1).
- [x] `ScaleLedger` + `thermal_reference` — 실제 쓴 스케일만
- [x] `Check` — 모델/적분/기하/통계 4분류 + 하드/소프트 + 여유
- [x] `dt = 10⁻²·γ/(국소 강성)` — 트랩 `γ/k`와 페어 `γ/U''(r_min)`이 같은 공식이었다
- [x] `report.render` · `runid` · `metrics` · `stats` · `sim` · `materials` · `provenance`
- [x] 두 케이스를 리팩터해 같은 경로로 돌림. 케이스에는 고유 물리만 남음
      (`trap_2d_5um.py` 610→436줄, `soft_r3_2d.py` 755→703줄)
- **DoD 달성 — 측정으로 확인**: [`scratch/verify_1c_equivalence.py`](../../verify/verify_1c_equivalence.py)
  - `run_id` **8개 전부 보존** (스펙 해시가 바이트 단위로 동일)
  - 1-A 재실행: metrics **77필드 동일**, 위반 0. 관측량 4종 6자리까지 일치
  - 1-B A=100 재실행: metrics **124필드 동일**, 위반 0
  - 허용된 변경 3종뿐:
    ① 벽시계 ② 블록 SEM이 float32→float64 업캐스트로 1.2e-6 달라짐(**더 정확한 쪽**.
    원본 trace가 float32였음) ③ 에너지 일관성 `err_pct`의 **부호 규약 통일 — 버그 수정**
    (1-B 원본은 이 한 행만 `(예측−측정)/|측정|` 이라 같은 파일의 다른 행들과도, 1-A와도
    부호가 반대였다. 크기·판정은 그대로)
- 부수 검증: [`scratch/verify_bdbot.py`](../../verify/verify_bdbot.py) — 공통 모듈이 케이스
  실측값을 재현하는지 (물성 5종 · dt 규약 · 하드/소프트 판정 · 자기상관 보정 · 함정 방어)

- **`SimSpec` + `derived_from` + Builder/Runner는 만들지 않았습니다.** 두 케이스가
  공유한 것은 "스펙 dict → sha256 → run_id"라는 **규약**이었고, 그건 `runid.py`로
  충분합니다. pydantic 스키마와 빌더는 아직 한 번도 두 번 나오지 않았습니다.

#### Phase 1-D 이후 · 케이스 추가마다 반나절
`chain-bend` → `trap-drag` → `abp-rod` 순. 각 케이스가 모듈 하나씩을 추가하고,
그때마다 스케일 원장이 자연히 자랍니다 (§5.6 모듈이 기여하는 구조).

**총 예상: 2~2.5일** (이전 3.5일). 더 중요한 건 **매 반나절마다 동작하는 결과가 나온다**는 점입니다.

### Phase 2 — 물리 검증 (골든 테스트) (1.5일) ⭐️
**LLM보다 먼저 물리가 맞는지 증명합니다.**

> ⚠️ **원칙 8**: 아래는 전부 `implementation_check` 입니다 — 예측이 구현한 모델에서
> 유도되므로 불일치는 버그입니다. **물리적 발견은 여기서 나오지 않습니다.**
> 발견을 원하면 시뮬레이션이 부과하지 않는 가정을 시험하는 케이스를 따로 설계해야 합니다
> (예: `soft-r3` 의 녹음 Γ, `trap-drag` 의 유효점도 성립 여부).
- [ ] `test_free_particle_msd`: `⟨Δr²⟩ = 2dD_t t`, D를 5% 이내 복원
- [ ] `test_diffusion_from_gamma`: γ 변화 → `D = kT/γ` 확인
- [ ] `test_wca_rdf_dilute`: 희박 극한 g(r) → `exp(-βU(r))`
- [x] `test_abp_effective_diffusion`: **`D_eff = D_t + v₀²/(d·Λ)`**, `Λ` = director 감쇠율.
      ~~`v₀²/(2(d−1)D_r)`~~ 은 2D에서만 우연히 맞음 (3D 실측 +29~31% 오차).
      HOOMD는 `Λ = rotational_diffusion` (2D·3D 모두). → `scratch/standalone_abp_diffusion.py`
- [x] `test_nondim_roundtrip`: 물리 → 무차원 → 물리 항등성 — **왕복 오차 ≤ 1.2e-16**
      (길이²·시간·확산계수·에너지 4종). `NondimSpec.physical(v, L=,T=,E=)` 을 차원 지수로
      호출하고 기준 스케일로 되나눠 항등성 확인. 저장→로드한 스펙도 동일 값
      (`scratch/verify_nondim_guards.py` ⑥⑦). pytest 로 옮길 때 그대로 쓸 수 있다.
- [ ] `test_scale_invariance` ⭐️: **같은 물리계를 `thermal` 기준과 `interaction` 기준으로 각각
      무차원화해 실행 → 물리 단위로 역변환한 결과가 일치**. 무차원화가 올바름을 증명하는 핵심 테스트.
- [ ] `test_separation_gates`: 각 분리 검사가 인위적 위반 케이스를 정확히 잡는지 (8종)
      — L3 층(원장 완전성·비 정합성·역변환)은 `verify_nondim_guards.py` 30종으로 완료.
      물리 분리 검사(§6.4의 게이트) 쪽은 아직 남았다.
- [ ] `test_reproducibility`: 같은 seed → 동일 궤적
- **DoD**: `pytest tests/test_golden_physics.py` 전부 통과.

### Phase 3 — 로우데이터 계층 + 온디맨드 분석 (2일)
- [ ] Tier A/B/C/D/L 라이터 구현
- [ ] 프레임 인덱스 + 지연 로더 + 슬라이싱
- [ ] 관측량 레지스트리 + 캐시 (13.1)
- [ ] 표준 플롯 세트
- [ ] CLI: `bdbot analyze <run_id> --obs msd,rdf`, `bdbot raw <run_id> --steps 1e6:2e6 --what forces`
- **DoD**: 완료된 런에서 임의 구간의 힘 분포를 30초 안에 뽑을 수 있다.

### 🟡 Phase 4 — 자동 진단 + 사후분석 (**자동 진단·record.json 완료**, 나머지 0.5일)
완료: 실패 taxonomy · 자동 진단(평형·드리프트·통계·편향) · `record.json` 방출 · KB 질의
→ [`tools/postmortem.py`](../../tools/postmortem.py), [`tools/kb.py`](../../tools/kb.py)

**1-B에서 고친 것 — 진단 3건이 1-A 전용 가정이었습니다** (전부 *틀린 판정*을 냈습니다):
1. **평형 지표가 '앵커로부터의 변위'**였습니다 → 속박계(트랩) 전용. 확산계는 변위가
   무한히 자라 항상 `EQ_INSUFFICIENT`. → 케이스가 `metrics.equilibration`으로 지표를 선언
2. **드리프트 t 검정에 자기상관 보정이 없었습니다** → 전 구간 변화가 평균의 **−0.026%**인
   런이 `t=−3.3`으로 `NUM_DRIFT` 판정. Sokal 자동창으로 `n_eff`를 구해 SE를 팽창시키고,
   **유의성과 크기를 함께** 봅니다(0.5% 미만은 실패로 보지 않음).
   → 1-A의 "오차막대는 블록 평균으로" 교훈과 **같은 종류의 실수**였습니다
3. **분리 검사 실패를 전부 하드로 취급**했습니다 → bd-physics §4는 통계·유한크기를
   ⚠ 경고로 규정. 1-A에서는 전부 통과해서 구분이 드러나지 않았습니다
4. 통계 목표 0.5%가 하드코딩(1-A의 `⟨x²⟩` 기준) → `numerics.stat_target_pct`로 케이스가 선언
남음: LLM 서술(§10.3) · 규칙 승격(§10.4, 표본 부족) · 승인 원장
- [ ] 실패 taxonomy + 자동 진단기 (10.2)
- [ ] `PostMortem` 생성 (LLM 서술 포함)
- [ ] 승인 원장
- **DoD**: 일부러 dt를 키운 런이 `NUM_DIVERGE`로 자동 분류되고, 평형화 짧은 런이 `EQ_INSUFFICIENT`로 잡힌다.

### Phase 5 — 지식베이스 (2.5일)
- [ ] SQLite + FTS5 스키마, 세 가지 검색 모드
- [ ] 논문 증류 파이프라인 (텍스트 + 그림)
- [ ] 검증 UI
- [ ] 사후분석 → KB 환류
- [ ] KB 기반 검증 경고 (12.3)
- [ ] 문헌 대조 플롯
- **DoD**: ABP 논문 3편을 증류해 넣고, "Pe≈40, φ≈0.6" 검색으로 관련 항목이 나온다. 우리 런의 실패가 tier 3 엔트리로 쌓인다.

### Phase 6 — 병렬 스윕 (1.5일)
- [ ] SQLite 잡 큐 + 프로세스 풀 8
- [ ] `SweepSpec` 전개 + 중복 제거 + 스윕 플롯 (상평형도 heatmap)
- **DoD**: Pe×φ 30런 스윕이 8워커로 돌고 heatmap이 나온다.

### Phase 7 — 인테이크 스키마 + 승인 절차 (0.5일) ⬇️ 축소
vision은 Claude Code가 무료로 제공하므로, 우리가 만드는 건 **스키마와 확정 절차**뿐입니다.
- [ ] `Observation` 스키마 + YAML 직렬화
- [ ] `bdbot intake init <folder>` — 빈 `observation.yaml` 템플릿 생성
- [ ] `bdbot intake check <folder>` — `ambiguities`/`unread_regions` 미기재 시 거부
- [ ] skill `bd-intake` — 해석 프로토콜 (§8.2: 전사 우선, 모호성 명시, 결측은 null)
- **DoD**: `/bd-intake intake/abp-rod-2d-run-flip` → 스케치를 읽고 `observation.yaml`이 생성되며,
  모호한 항목이 명시적으로 나열되고 사람이 확인·수정한다.

### Phase 8 — Claude Code 통합 표면

Claude Code가 에이전트 런타임이므로(§15), 우리가 만드는 건 **지식 캡처 + 강제 장치**입니다.
두 덩이로 나눕니다 — 앞의 것은 지금 당장 가치가 있고, 뒤의 것은 CLI가 생긴 뒤라야 의미가 있습니다.

#### ✅ Phase 8-min · 지식 캡처 **(완료 2026-08-03)**
세션이 끝나면 이번에 알아낸 것이 전부 사라집니다. 그걸 막는 게 목적이었습니다.
- [x] [`CLAUDE.md`](../../CLAUDE.md) — 불변식 6개 (54줄). 6번은 "물리 주장은 검증하고 말한다"
- [x] skill [`bd-hoomd`](../../.claude/skills/bd-hoomd/SKILL.md) (407줄) — 함정 9개(★ 5개는
      조용히 틀리는 유형) + 하드 제약(병진 이방성 불가, 실측 근거) + 검증된 스니펫 14개 + API 참조
- [x] skill [`bd-physics`](../../.claude/skills/bd-physics/SKILL.md) (209줄) — 5단계 절차 ·
      단위 규약 · 무차원수=스케일비 · 분리검사(전부 10⁻²) · **케이스별 스케일 표**(트랩 검증 완료)
- [x] `scratch/verify_skill_snippets.py` — **문서의 코드를 추출해 실행하는 회귀 테스트**
- **DoD 달성**: 스킬 문서에서 추출한 트랩 코드가 그대로 돌고 오차 +0.46%, NaN 0개.
  문법 검사 14/14. 문서가 깨지면 테스트가 깨집니다.
- **왜 했나**: 이번 세션에서 두 번 틀렸습니다 — 최소 이미지 누락(k=2에서 +1856%),
  강체 이방성 오판. 기록하지 않으면 다음 세션에서 같은 실수를 합니다.

#### Phase 8-rest · 강제 장치 + 워크플로 (1일) ← **CLI가 생긴 뒤**
- [ ] hooks 4종: `guard_invariant`, `guard_separation`, `guard_cost`, `log_approval`
- [ ] `.claude/settings.json` — 권한 허용목록 (`Bash(bdbot:*)`),
      `specs/`·`runs/` 쓰기 금지(§부록 B.6 두 번째 방어선)
- [ ] subagents 3종: `bd-distiller`, `bd-analyst`, `bd-reviewer`
- [ ] slash commands 7종
- [ ] skill `bd-intake`, `bd-distill`
- **DoD**: `derived_from` 없는 스펙으로 `bdbot run`을 시도하면 **훅이 차단**하고 수정안을 제시한다.
  스케치 폴더 하나에서 `/bd-intake` → `/bd-spec` → `/bd-run` → `/bd-postmortem`이 끝까지 흐른다.
- **왜 나중**: 훅은 `bdbot run`을 가로채는 장치입니다. `bdbot run`이 없으면 지킬 대상이 없습니다.

### Phase 9 — (선택) 확장
수렴 연구 자동화 · 자동 리포트 생성 · GPU 지원 · 능동학습 · **MCP 서버화(필요 시)** ·
이방성 병진 마찰 모듈(§20 질문 10 옵션 B)

**총 예상: 13~15일** (풀타임 기준). **Phase 1(케이스 관통) → 2(골든 검증)가 임계 경로**입니다.

### 16.5 권장 실행 순서

```
✅ Phase 0       범용 환경 + HOOMD 능력 조사 + 15개 API 검증        완료
✅ Phase 8-min   CLAUDE.md + bd-hoomd/bd-physics skill              완료
✅ Phase 1-A     trap-2d-5um 관통 — 관측량 4종 해석해 일치           완료
   Phase 1-B     soft-r3 관통 → 공통점 대조 (함께 판단)            반나절   ⭐️ ← 다음
   Phase 1-C     실제로 두 번 나온 것만 추상화 + Builder/Runner     반나절
   Phase 2       골든 테스트 정식화 (trap은 이미 통과)              1일
   Phase 1-D~    chain-bend → trap-drag → abp-rod                각 반나절
   Phase 3       로우데이터 계층 Tier B~D + 온디맨드 분석            1.5일
   Phase 4       사후분석 나머지 (LLM 서술 · 규칙 승격)              0.5일  ← 자동 진단은 완료
   Phase 7       인테이크 스키마                                   0.5일
   Phase 8-rest  훅 + 워크플로 + 서브에이전트                       1일
   Phase 5       지식베이스                                        2.5일
   Phase 6       병렬 스윕                                         1.5일
```

**Phase 1은 매 반나절마다 동작하는 결과가 나옵니다.** 프레임워크를 먼저 만들지 않고
케이스를 관통하므로, 중간에 방향을 바꿔도 버리는 게 적습니다.

---

## 17. 테스트 전략

| 층위 | 내용 | 현황 (2026-08-03) |
|---|---|---|
| **골든 물리** | 해석해로 적분기 정확성 증명 | ✅ 트랩 `⟨x²⟩=kT/k` k 4종 오차<0.6%, `⟨x²⟩·k` 변동 0.28% · 1-A 관측량 4종 |
| **문서 코드 검증** ⭐️ | 스킬 스니펫을 추출해 실행 | ✅ `scratch/verify_skill_snippets.py` (문법 14/14 + 트랩 물리) |
| **dt 편향 법칙** ⭐️ | `편향 = (dt/τ)/2` 확인 | ✅ `scratch/dt_convergence.py` 4점, 이론과 0.07% 이내 |
| **자동 진단** | 평형·드리프트·통계·편향 일관성 | ✅ `tools/postmortem.py` |
| **재현성** | 고정 seed → 동일 결과 | ✅ 1-A 재실행 자릿수 동일 |
| **무차원화 왕복** | 물리 → 무차원 → 물리 항등성 | 부분 — 1-A가 왕복하나 `test_scale_invariance`(기준 2종 비교) 미작성 |
| **수렴 연구** ⭐️ | `dt` 절반, `N` 2배로 불변 | ❌ **미실시** — 1-A `record.json.not_verified`에 명시 |
| 회귀 | 고정 seed 짧은 런 → 스냅샷 diff | ❌ |
| 스키마 | SimSpec 라운드트립, run_id 안정성 | 부분 — run_id 콘텐츠 주소 지정 동작 |
| 검증 로직 | 잘못된 스펙 → 정확한 에러 | 부분 — 분리 검사 게이트가 실제로 스모크를 차단함 |
| **증류 정확도** ⭐️ | 논문 5편 정답셋 | ❌ (문헌 0편) |
| 에이전트 평가 | 요청 10종 사람 채점 | ❌ |

**빠른 반복**: 케이스 스크립트에 `--smoke`(N=200, 짧게) · `--report`(리포트만) · `--force`(재실행).
`--smoke`도 분리 검사를 통과해야 합니다 — 실제로 초기 설정이 통계 검사에 걸려 차단됐습니다.

---

## 18. 리스크 & 완화

| 리스크 | 완화 |
|---|---|
| HOOMD API가 예상과 다름 | Phase 0에서 전부 실증 후 진행 |
| **LLM이 논문 값을 환각** ⭐️ | 출처(locator) 필수 · tier 시스템 · 사람 검증 · 증류 정확도 테스트 |
| **스케치 오독** ⭐️ | 전사 우선 · `ambiguities` 강제 · 필드별 승인 |
| LLM이 그럴듯하지만 틀린 물리 스펙 | 검증 레이어 + 골든 테스트 + 관계식 하드코딩 |
| dt가 커서 조용히 부정확 | 하드 상한 + 런타임 감시 + 수렴 연구 유틸 |
| **로우데이터 디스크 폭발** ⭐️ | 계층 저장 + 용량 사전 추정 + 승인 게이트 + GC 정책 |
| **KB가 쓰레기로 오염** ⭐️ | tier 시스템 · 검증 워크플로 · 반증되면 엔트리에 `retracted` 표시 (삭제 아님, 이력 보존) |
| BD 가정 위반 계에 BD 적용 | `Re`/`St` 자동 검사 |
| 맥북 CPU 한계 | 비용 추정 + 승인 게이트 + smoke 모드 |
| 평형화 부족 모르고 분석 | 에너지 시계열 **항상** 플롯 + 자동 평형 판정 |
| API 비용 | 프롬프트 캐싱 · 짧은 구조화 출력 · 로컬 파이프라인은 LLM 없이 동작 |

---

## 19. 기술 스택

| 용도 | 선택 | 이유 |
|---|---|---|
| 시뮬레이션 | HOOMD-blue 7.1.1 (`*cpu*`) | 지정 |
| 스키마 | Pydantic v2 | CLI 검증 + YAML 직렬화 + (배치 증류 시) `messages.parse()` 직결 |
| **단위** | **pint** | 차원 있는 계산의 단위 실수 차단 |
| 궤적 | GSD (+ `Burst`) | HOOMD 네이티브 |
| 로그 | HDF5 (`hoomd.write.HDF5Log`) | 대용량 시계열 |
| 분석 | freud + numpy/scipy | HOOMD와 궁합 |
| 테이블 | pandas + pyarrow (Parquet) | 스윕 비교 |
| **KB** | **SQLite + FTS5** | 로컬 규모엔 벡터DB 과잉, 좌표 검색은 SQL이 정확 |
| 잡 큐 | SQLite | 서버 불필요 |
| 플롯 | matplotlib | 의존성 최소 |
| **에이전트 런타임** | **Claude Code** | 지정 — vision·PDF·대화·승인이 전부 내장 |
| **통합 방식** | **CLI + Skill + Hooks** (MCP 아님) | §15.2 |
| LLM 직접 호출 | `anthropic` — **선택 의존성** `bdbot[distill]` | 논문 배치 증류에만 (§15.9) |
| CLI | typer | 서브커맨드 구조 + `--format json` |
| 테스트 | pytest | — |

---

## 20. 열린 질문

### 해결됨 (v0.4 기준)
| # | 질문 | 결론 |
|---|---|---|
| 3 | 스케치 형식 | 종이+펜 사진, Claude Code Read로 판독 가능. 정보 밀도 낮음(파라미터 2~3개) |
| 4 | 기본 기준 스케일 | `thermal` 고정 (σ=d, E=kT, τ=τ_B). **단위계는 계마다 바꾸지 않음** |
| 4b | 분리 검사 임계 | 전부 `10⁻²`로 통일. `dt/τ=1e-2` ⟺ 편향 0.5% (실측 확인) |
| 6 | 인터페이스 | Claude Code 런타임 + CLI/Skill/Hooks (MCP 아님) |
| 9 | 스코프 | 제한하지 않음. 물리 모듈 레지스트리(§5.6)로 확장 |
| — | 매질 | 물(뉴턴). 점탄성은 별도 케이스로 (`medium.*` 모듈) |
| — | KB 시점 | `record.json`부터, SQLite는 런 100건/문헌 유입 시 (§7.0) |

### 남은 질문

**A. 이방성 병진 마찰** — 아래 §10에 상세. `abp-rod` 케이스의 전제 조건. **결정 필요.**

**B. 수렴 연구를 언제 할 것인가** — 1-A는 `dt` 절반 재실행을 하지 않고 편향 법칙 예측과의
대조로 갈음했습니다(`record.json.not_verified`에 명시). 선형계라 정당하지만, 비선형계
(1-B soft-r3부터)에는 닫힌 형태가 없어 **직접 수렴 연구가 필요**합니다. 1-B에서 표준
절차로 넣을까요, 아니면 `tools/converge.py` 유틸로 분리할까요?

**C. 관측량 우선순위** — Phase 3에서 먼저 구현할 5개.
추천: `msd`, `rdf`, `local_density_hist`, `cluster_size`, `force_distribution`.
(1-A에서 이미 `P(x)`, `⟨x²⟩`, `C(t)`, `PSD`가 케이스 코드 안에 구현됨 — 1-C에서 공통화 대상)

**D. 논문 소스** — PDF 직접 투입 vs arXiv/DOI 취득. 후자면 네트워크 접근 정책 필요.
KB 시드로 넣을 논문 5~10편을 정해주시면 증류 프롬프트를 그 스타일에 맞춰 튜닝합니다.
(현재 문헌 0편 — Phase 5까지는 급하지 않음)

**E. 노트북 지원** — 분석 결과를 Jupyter에서 만지실 거면 Phase 1-C 추상화 때
`bdbot.api`를 1급으로 설계해야 합니다.

**F. `mater_plan.md` 파일명** — `master_plan.md` 오타로 보이나 요청하신 이름 그대로
유지 중입니다. 바꿀까요?

### A. 이방성 병진 마찰 — 어떻게 할까요? (**결정 필요**)

강체로 묶으면 이방성이 나온다는 앞선 판단은 **실측 결과 틀렸습니다** (§5.6 주석,
[`docs/hoomd_capabilities.md` §5.1](../../docs/hoomd_capabilities.md)). 강체·비드사슬 모두 `γ⊥/γ∥ = 1.000000`.
이방성은 HI의 효과이고 BD에는 HI가 없습니다.

`abp-rod-2d-run-flip`에 대한 실제 영향:

| 관측량 | 정확도 |
|---|---|
| **MSAD** | ✓ 정확 (`γ_r,z`만 의존, 2D라 회전축 하나뿐) |
| **MSD 장시간** (t ≫ τ_r) | ✓ 정확 (등방 평균 `γ̄`만 의존) |
| **MSD 단시간** (t < τ_r) | ✗ 이방성 손실 |

run-and-flip은 flip 사이 방향이 오래 유지되므로 단시간 이방성이 실제로 보일 수 있습니다.
스케치의 "measure MSD, MSAD"가 그 신호를 겨냥한 것인지가 판단 기준입니다.

| 옵션 | 방법 | 비용 | 결과 |
|---|---|---|---|
| **A (권장 기본)** | Perrin/slender-body로 종횡비에서 등방 평균 `γ̄` 계산. 회전은 `γ_r` 텐서로 정확히. 모듈이 한계 선언 | 낮음 | MSAD·장시간 MSD 정확 |
| **B** | 커스텀 적분기로 이방성 부과 (`md.half_step_hook`). 결정론 항은 보정력으로 되나 **잡음항이 요동-소산을 깸** — 신중한 설계 필요 | 높음 | 정확 |
| **C** | `hoomd.mpcd`로 HI 도입 | 매우 높음 | 정확, 사실상 다른 프로젝트 |

**A로 시작하고 B는 별도 모듈로 남겨두는 것**을 권합니다 — 모듈 구조라 나중에 붙여도
코어를 건드리지 않습니다. 다만 단시간 MSD 이방성이 이번 연구의 **핵심 관측량**이라면
B를 Phase 1b에 넣어야 하니, 스케치 의도를 알려주세요.


`intake/`에 들어온 5개 케이스를 보니 §1에서 잡은 스코프(구형 콜로이드 + 구형 ABP)로는 **하나도 다룰 수 없습니다.**

| 케이스 | 필요한 물리 | 현재 스코프 | HOOMD 대응 |
|---|---|---|---|
| `trap-2d-5um` | **외부 조화 트랩** | ❌ | `md.force.Custom` 또는 `md.external.field` |
| `trap-drag-2d-hex300` | 트랩 + **이동 구동**, 육방 초기배치 | ❌ | 위 + `hoomd.variant`로 트랩 중심 이동 |
| `chain-bend-2d-oscill` | **본드 + 굽힘 강성 + 시간의존 구동** | ❌ (비스코프) | `md.bond.FENE/Harmonic` + `md.angle` + `variant` |
| `abp-rod-2d-run-flip` | **타원체(비구형)** + **run-and-flip**(이산 사건) | ❌ | rigid body 또는 이방성 퍼텐셜 + **커스텀 updater** |
| `soft-r3-2d-A-sweep` | **소프트 `r⁻³` 퍼텐셜**, 진폭 스윕 | ❌ | `md.pair.Table` 또는 커스텀 |

특히 `abp-rod`는 계획의 물리 공식에 직접 영향을 줍니다:
- **§6.1의 `D_r = 3D_t/d²`는 구에만 성립**합니다. 타원체는 Perrin friction factor가 필요하고,
  병진 확산이 **이방성**(장축/단축)이며 회전 확산도 축마다 다릅니다.
- **run-and-flip은 ABP가 아닙니다.** ABP는 연속 회전확산(`create_diffusion_updater`)이지만,
  run-and-flip은 **이산 사건**(포아송 과정으로 180° 반전)입니다. HOOMD 내장 기능으로 안 되고
  커스텀 updater가 필요합니다.
- 스케치의 `τ_R = 0.5 s`가 회전확산 시간인지 flip 간격인지 **모호**합니다 (전형적인 `ambiguities` 항목).

**제안**: §1 스코프를 아래 6개 **물리 모듈**로 재정의하고, 각각을 플러그인으로 만듭니다.
`SimSpec`을 모듈 조합으로 설계하면 5개 케이스를 전부 커버하고 확장도 열립니다.

| 모듈 | 내용 | 우선순위 |
|---|---|---|
| `M1 pair` | WCA/LJ/Yukawa/**Table(임의 r⁻ⁿ)** | 필수 |
| `M2 external` | **조화 트랩**, 이동 트랩, 벽 | 높음 (5개 중 2개) |
| `M3 bonded` | 본드 + **굽힘 각도** (사슬) | 높음 (1개) |
| `M4 active` | 연속 ABP + **이산 run-and-tumble/flip** | 높음 (1개) |
| `M5 anisotropic` | **타원체/막대** (Perrin 마찰, 이방성 확산) | 중간 (1개, 가장 어려움) |
| `M6 driving` | **시간 의존 구동** (`hoomd.variant`: 진동, 램프) | 높음 (2개) |

이건 **일정에 영향**을 줍니다: Phase 1b(Builder)가 1.5일 → 3~4일로 늘고, M5(비구형)는 별도 Phase가 필요할 수 있습니다.

세 가지 선택지:
- **(a)** 5개 케이스 전부 커버하도록 스코프 확장 — 정직하지만 일정 +4~6일
- **(b)** 쉬운 것부터: M1+M2+M6 먼저(트랩 2개 케이스), M3/M4/M5는 나중 — 2주 안에 첫 결과
- **(c)** 대표 1개만 골라 끝까지 관통(수직 슬라이스) 후 확장 — 파이프라인 검증이 가장 빠름

**제 추천은 (c)입니다.** `trap-2d-5um`이 가장 단순하면서(입자 1개 + 조화 트랩) 해석해가 있어서
**골든 물리 테스트로 바로 쓸 수 있습니다** — 트랩 안 BD 입자의 위치 분포는 정확히
`P(x) ∝ exp(-kx²/2k_BT)`이고 `⟨x²⟩ = k_BT/k`, 이완시간은 `τ = γ/k`. 무차원화·검증·로우데이터·
분석·사후분석·KB 환류 전 경로를 최소 비용으로 관통 검증한 뒤 M3/M4/M5를 붙이는 게 안전합니다.

---

## 부록 A — Phase 0 스모크 테스트 스크립트

```python
# scratch/hello_bd.py — API 실증용 최소 예제
import itertools, math
import numpy as np, gsd.hoomd, hoomd

N_SIDE, PHI, KT, GAMMA, DT = 40, 0.5, 1.0, 1.0, 1e-4
N = N_SIDE ** 2
L = math.sqrt(N * math.pi / (4 * PHI))            # 2D 면적분율 → 박스 길이 (σ=1)

a = L / N_SIDE
pos = np.array([[(i + .5) * a - L / 2, (j + .5) * a - L / 2, 0.]
                for i, j in itertools.product(range(N_SIDE), repeat=2)])

frame = gsd.hoomd.Frame()
frame.particles.N = N
frame.particles.position = pos
frame.particles.orientation = [(1, 0, 0, 0)] * N
frame.particles.typeid = [0] * N
frame.particles.types = ['A']
frame.configuration.box = [L, L, 0, 0, 0, 0]      # Lz=0 → 2D
frame.configuration.dimensions = 2

sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_snapshot(frame)

cell = hoomd.md.nlist.Cell(buffer=0.4)
lj = hoomd.md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode='shift')  # = WCA
lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)

bd = hoomd.md.methods.Brownian(filter=hoomd.filter.All(), kT=KT, default_gamma=GAMMA)
integrator = hoomd.md.Integrator(dt=DT, methods=[bd], forces=[lj])
integrator.integrate_rotational_dof = False       # ABP 확장 대비 명시
sim.operations.integrator = integrator

# --- Tier A: 위치/방향 ---
sim.operations.writers.append(
    hoomd.write.GSD(filename='traj_A.gsd', trigger=hoomd.trigger.Periodic(10_000),
                    mode='xb', dynamic=['property']))

# --- Tier B: per-particle 힘 (Phase 0에서 반드시 실증할 것) ---
plog = hoomd.logging.Logger(categories=['particle'])
plog.add(lj, quantities=['forces', 'energies'])
sim.operations.writers.append(
    hoomd.write.GSD(filename='traj_B.gsd', trigger=hoomd.trigger.Periodic(50_000),
                    mode='xb', logger=plog, dynamic=['property', 'momentum']))

# --- Tier L: 전역 스칼라 ---
thermo = hoomd.md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)
glog = hoomd.logging.Logger()
glog.add(thermo, quantities=['potential_energy', 'pressure'])
sim.operations.writers.append(
    hoomd.write.Table(trigger=hoomd.trigger.Periodic(10_000), logger=glog))

sim.run(100_000)
```

> Phase 0의 목적은 **이 스크립트가 v7.1.1에서 그대로 도는지** 확인하는 것입니다.
> API 시그니처가 하나라도 다르면 11절 매핑 표와 빌더 설계를 그때 고칩니다.
> 특히 `logger=` 인자를 통한 per-particle 힘 저장은 Tier B의 전제이므로 반드시 검증하세요.

---

## 부록 B — Claude Code 통합 산출물 스케치

### B.1 `bdbot` CLI 표면 (Claude Code가 호출하는 전부)

```
bdbot init-system --template <name> [--from-observation FILE]  → system.yaml 골격
bdbot nondim system.yaml [--strategy thermal] [-o specs/x.json] → SimSpec + 리포트
bdbot validate specs/x.json [--format json]                     → 검증 리포트 (종료코드로 통과 여부)
bdbot estimate specs/x.json                                     → wall-time / 디스크 추정
bdbot run specs/x.json [--background]                           → run_id
bdbot sweep sweeps/y.yaml                                       → sweep_id (전 구간 사전 검사)
bdbot status [run_id|sweep_id] [--format json]                  → 진행률/상태
bdbot resume <run_id>                                           → 체크포인트에서 재시작
bdbot analyze <run_id> --obs msd,rdf [--format json]            → observables.parquet
bdbot raw <run_id> --what forces --steps 3e6:4e6 [--who tracers] → Parquet 추출 + 요약통계
bdbot plot <run_id|sweep_id> --kind standard|phase|literature   → PNG 경로들
bdbot redimensionalize <run_id> --obs msd                       → 물리 단위 값
bdbot kb search --tags 2D,ABP --coords Pe=40,phi=0.6 [--format json]
bdbot kb add entry.yaml   |   bdbot kb verify <id>   |   bdbot kb conflicts
bdbot postmortem <run_id> [--format json]                       → 자동 진단 (LLM 없이)
bdbot intake init <folder>  |  bdbot intake check <folder>
bdbot gc --dry-run
```

**설계 규칙**
- 모든 명령에 `--format json` — Claude Code가 파싱하기 쉽게
- 종료 코드로 통과/실패 — 훅이 판정에 사용
- 에러 메시지에 **항상 구체적 수정안** 포함 (§12.1 예시)
- 사람이 직접 써도 똑같이 동작 (Claude Code 의존성 없음)

### B.2 훅 — 하드 불변식 강제 (`.claude/hooks/guard_invariant.py`)

```python
#!/usr/bin/env python3
"""PreToolUse hook: block `bdbot run` on specs that bypass the dimensional layer.

원칙 3(차원 우선)은 프롬프트가 아니라 하네스가 강제한다.
stdin으로 훅 입력(JSON)을 받고, 차단 시 stderr에 사유+수정안을 쓰고 exit 2.
"""
import json, re, sys
from pathlib import Path

payload = json.load(sys.stdin)
cmd = payload.get("tool_input", {}).get("command", "")

m = re.search(r"\bbdbot\s+run\s+(\S+)", cmd)
if not m:
    sys.exit(0)                                   # 관련 없는 명령 — 통과

spec_path = Path(m.group(1))
if not spec_path.exists():
    sys.exit(0)                                   # 존재 여부는 CLI가 판단

spec = json.loads(spec_path.read_text())
missing = [k for k in ("derived_from", "scale_ledger") if not spec.get(k)]
if missing:
    print(
        f"BLOCKED: {spec_path} is missing {missing}.\n"
        "원칙 3(차원 우선) 위반 — 무차원 스펙은 반드시 PhysicalSystem에서 유도되어야 합니다.\n"
        "수정: bdbot init-system --template <name> > system.yaml\n"
        f"      bdbot nondim system.yaml -o {spec_path}",
        file=sys.stderr,
    )
    sys.exit(2)                                   # 2 = 차단, stderr가 Claude에게 전달됨

ledger = spec["scale_ledger"]
failed = [c for c in ledger.get("separations", [])
          if c["verdict"] == "fail"]
if failed:
    lines = "\n".join(
        f"  - {c['name']}: {c['ratio']:.2e} (한계 {c['threshold']:.0e}) — {c['message']}"
        for c in failed
    )
    print(f"BLOCKED: 스케일 분리 검사 실패 {len(failed)}건\n{lines}", file=sys.stderr)
    sys.exit(2)

sys.exit(0)
```

> 훅 등록은 `.claude/settings.json`에 들어갑니다. 이 파일은 손으로 편집하기보다
> `update-config` 스킬로 작성하는 편이 안전합니다.

### B.3 슬래시 명령 (`.claude/commands/bd-intake.md`)

```markdown
---
description: 스케치 폴더를 읽어 Observation을 작성하고 사람 확인을 받는다
---

인자로 받은 폴더: $ARGUMENTS

1. `bdbot intake init $ARGUMENTS` 로 `observation.yaml` 골격을 만든다.
2. 폴더 안의 모든 이미지를 Read 툴로 읽는다.
3. skill `bd-intake` 프로토콜을 따른다:
   - 먼저 `raw_transcription`에 **보이는 그대로** 전사한다. 해석하지 않는다.
   - 그 다음 `entities`, `stated_quantities`를 채운다.
   - 스케치에 **없는** 값은 절대 지어내지 않는다. `null`로 둔다.
   - `ambiguities`와 `unread_regions`를 반드시 채운다. 비어 있으면 다시 본다.
4. `observation.yaml`을 Write 한다.
5. `bdbot intake check $ARGUMENTS` 로 필수 필드를 검사한다.
6. 사용자에게 **항목별로** 확인을 받는다 (AskUserQuestion). 특히 모호 항목은
   해석 후보를 선택지로 제시한다.
7. 승인된 내용으로 `observation.yaml`을 갱신한다.

다음 단계는 `/bd-spec $ARGUMENTS` 이지만, 지금 자동으로 넘어가지 않는다.
```

### B.4 서브에이전트 (`.claude/agents/bd-reviewer.md`)

```markdown
---
name: bd-reviewer
description: 시뮬레이션 스펙을 적대적으로 검토한다. 스펙 제출 전에 호출.
tools: Read, Bash, Grep
---

너는 회의적인 시뮬레이션 물리 리뷰어다. 이 스펙을 **만든 사람이 아니며**,
그들의 의도를 모른다. 스펙 파일과 무차원 리포트만 보고 판단한다.

너의 임무는 승인이 아니라 **반박**이다. 다음을 순서대로 점검한다:

1. 무차원수가 목표 현상을 관측할 수 있는 영역에 있는가? 문헌과 대조하라
   (`bdbot kb search`).
2. 스케일 분리 검사에서 여유가 5배 미만인 항목이 있는가? 그것이 결과를
   왜곡할 수 있는가?
3. 박스가 관심 길이척도(상관길이, 지속길이)에 비해 충분한가?
4. 관측 시간이 가장 느린 시간척도의 100배 이상인가?
5. 물성값 중 tier 2 이하(미검증)가 결론을 좌우하는가?
6. 이 계에 Brownian dynamics를 쓰는 것이 타당한가? (Re, St, 과감쇠 가정)

각 항목에 대해 **문제 없음 / 우려 / 치명적** 중 하나로 판정하고, 우려·치명적에는
구체적 수정안을 제시한다. 확신이 없으면 "확신 없음"이라고 말한다.
찾지 못한 문제를 지어내지 않는다.
```

### B.5 스킬 (`.claude/skills/bd-physics/SKILL.md` — 발췌)

```markdown
---
name: bd-physics
description: |
  Brownian dynamics 시뮬레이션의 단위계, 스케일 원장, 무차원화 규약.
  물리계를 정의하거나, 파라미터를 제안하거나, 무차원수를 해석하거나,
  스케일 분리를 판단할 때 읽어라.
---

# 차원 우선 워크플로

... (마스터플랜 §6.1 스케일 목록 표 전체) ...
... (§6.2 기준 선택 전략 표) ...
... (§6.3 무차원수 = 두 스케일의 비 표) ...
... (§6.4 분리 검사 임계 표) ...

## 절대 금지
- 무차원 값을 먼저 정하고 물리계를 나중에 유추하지 않는다.
- 구형 Stokes 관계(D_r = 3D_t/d²)를 비구형 입자에 적용하지 않는다.
- 스케치나 문헌에 없는 물성값을 지어내지 않는다.
```

### B.6 권한 설정 (`.claude/settings.json` 골자)

```json
{
  "permissions": {
    "allow": [
      "Bash(bdbot:*)",
      "Read(//Users/kyuhwan/Desktop/simulation_auto/**)",
      "Write(//Users/kyuhwan/Desktop/simulation_auto/intake/**)",
      "Write(//Users/kyuhwan/Desktop/simulation_auto/kb/**)"
    ]
  },
  "hooks": {
    "PreToolUse": [
      {
        "matcher": "Bash",
        "hooks": [
          {"type": "command", "command": ".claude/hooks/guard_invariant.py"},
          {"type": "command", "command": ".claude/hooks/guard_cost.py"}
        ]
      }
    ]
  }
}
```

읽기는 프로젝트 전체를 허용하되, **쓰기는 `intake/`와 `kb/`로 제한**합니다.
`specs/`와 `runs/`는 `bdbot` CLI만 쓰게 해서, 손으로 만든 스펙이 섞이는 것을 구조적으로 막습니다
(원칙 3 강제의 두 번째 방어선).
