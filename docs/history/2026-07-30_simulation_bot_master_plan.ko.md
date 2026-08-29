# Simulation Bot — Master Plan

> Brownian Dynamics 시뮬레이션 에이전트. 손그림 한 장에서 검증된 결론까지.
>
> 상태: **설계 확정 v0.1 / 결정론 코어 완성 / 첫 물리 캠페인 완주** · 최종 수정 2026-07-30
>
> **규모** 커밋 28 · `simbot` 19 모듈 · `scripts` 18 · 테스트 **562** ·
> systems 카드 11 · findings 23 · 프로덕션 런 **3,856** (배치 wall 85 분)

---

## 진행상황 대시보드

**표기 규약** — `[O]` 완료·검증됨 · `[~]` 진행중 · `[X]` 미착수 · `[-]` 이번 버전 비범위
> 이 대시보드가 진행상황의 단일 진실 원천이다. 작업을 마치면 **코드보다 먼저 여기를 갱신**한다.

### 마일스톤
| | 마일스톤 | 상태 |
|---|---|---|
| `[O]` | **M0** 골격 (env·디렉터리·문서·knowledge·테스트) | 완료 |
| `[O]` | **M1+M2** 수직 슬라이스 — **손그림으로 S1→S8 전체 관통** | **완료.** 자유확산 대신 조화 트랩으로. 예측 9개 중 7 PASS / 2 INCONCLUSIVE / **0 FAIL** · 16런 wall **10.6 s** |
| `[O]` | **M2.5** 결정론 코어 — **관통을 코드로 재현 가능하게** | **완료 (2026-07-28).** `spec` `nondim` `io` `validate` `report` `policy` `session` `viz` + `cli.py`. 손그림 관통이 `cli.py run` 한 줄로 재현된다 (**2.3 s**, 첫 런과 비트 일치) |
| `[O]` | **M2.7** L1 에이전트 층 (`.claude/`) | **완료 (2026-07-28).** 스킬 3 + 참조문서 5 + 서브에이전트 9 + `settings.json`. 구조를 `test_agent_layer.py` (64개) 가 감시한다 |
| `[~]` | **M3** 물리 도메인 A→D 확장 | 트랩 분기 + **2D 소프트 반발계 완주** (§4.1). 도메인 B·C·D 는 카드만 |
| `[O]` | **M4** 검증 강화 | **완료 (2026-07-29).** dt래더·시드앙상블·자기일관성·판별력·`converge` + **부트스트랩·유한크기 스캔·설계검정력 사전계산·순차 사전등록** |
| `[O]` | **M4.5** 첫 물리 캠페인 — **소프트 반발 2D** | **완료 (2026-07-29~30).** 6 런 · 3,856 프로덕션 런. 카드가 `draft` → 벤치마크 31개. §4.1 |
| `[X]` | **M5** 인풋 모달리티 확장 | 손그림 완료 → 다음은 실험 데이터 대조 (A2 를 닫으려면 필요) |

### 첫 완주 결과 (run `2026-07-28_trap-2d-5um_2dfb9d`)
| | |
|---|---|
| 입력 | 손그림 1장 (2D 광집게, `R=5 μm`, `k=10 pN/μm`, `T=300 K`) |
| 확정 답 | `⟨x²⟩ = 416.58 ± 1.85 nm²` · `τ_trap = 8.0567 ± 0.0300 ms` · `f_c = 19.7 Hz` |
| 최고 정밀도 | `τ_fit/τ = 0.9998 ± 0.0008` (**0.08 %**) · MSD `R² = 0.99998` |
| 총 계산 | **10.6 s** (16런, 동시 8) |
| 잡힌 오류 | 4건 — S0 지식 수치 · S3 절단 solver · S5 GC detach · **S7 상관표본 KS 거짓기각** |

### M0 세부
| | 항목 | 비고 |
|---|---|---|
| `[O]` | `simulation_bot` conda env 생성 | python 3.12 · hoomd 7.1.0 · gsd 5.0.1 · numpy · scipy |
| `[O]` | 로컬 CPU 병렬성 검토 + 최적 코어 수 실측 | §7.3 · `wiki/findings/local-cpu-parallelism.md` |
| `[O]` | `master_plan.md` | 이 문서 |
| `[O]` | 진행상황 대시보드 | 이 절 |
| `[O]` | `config/run_policy.yaml` — 최적 조건 제안 + 사람 override | §7.4 |
| `[O]` | **`knowledge/` — BD_agent에서 이식** | §10. `source/papers` 42편 · `wiki/{systems 11, findings 23, benchmarks 5, concepts 2}` (2026-07-30 현재) + `wiki/CLAUDE.md` 계약 |
| `[O]` | `CLAUDE.md` 프로젝트 규약 | 제1원칙(제안하되 모르면 묻는다) 중심 |
| `[O]` | 감도 분석 설계 (S7b) | §11 |
| `[O]` | 모델 티어링 설계 | §12 |
| `[O]` | 첫 예시 인풋 저장 | `inputs/trap-2d-5um/sketch_01.jpeg` + sha256 |
| `[O]` | `environment.yml` + `wiki/techniques/env-log.md` | §7.1. HOOMD 7.1 호환 최신 11개 패키지 |
| `[O]` | `passive-sphere--harmonic-trap` 카드 | 첫 예시 그림의 (계×동역학) 쌍 |
| `[O]` | CLI·세션 층 범위 결정 | **전장 채택** (세션 + converge + params) |
| `[O]` | git 초기화 + 초기 커밋 | `4ac2a53`, 78파일. BD_agent와 독립. 현재 커밋 28개 |
| `[O]` | **S1 판독 + S2 예측 (첫 손그림)** | `runs/2026-07-28_trap-2d-5um_2dfb9d/` |
| `[O]` | `simbot/` 4개 모듈 | `units` `estimators` `forces` `guards` — 전부 테스트됨 |
| `[O]` | **테스트 스위트 (단계별)** | **562 통과 / 1 skip / 94 s.** `pytest -m "not slow"` 는 14 s |
| `[O]` | HOOMD 스킴 실측 판정 | EM 확정 · 노이즈 균일분포 · 입자 독립 · 배위온도 |
| `[O]` | **`simbot/` 결정론 코어** | `units` `estimators` `forces` `guards` `build` `run` `cutoff` `analysis/trap` `spec` `nondim` `io` `validate` `report` `policy` `session` `viz` — **16 모듈, 전부 테스트됨** |
| `[O]` | **`cli.py`** | `run` `resume` `converge` `params` `calibrate` 전부 동작 확인 |
| `[O]` | **`examples/trap-2d-5um/`** | 기계가 읽는 `spec.yaml` + `prediction.yaml`. 손으로 쓴 첫 런의 파생값 10개를 재현 |
| `[O]` | **`.claude/` 에이전트 층** | 스킬 3 · 참조문서 5 · 서브에이전트 9 · `settings.json`. §12.3 |
| `[O]` | **`simbot/viz.py`** | 그림 5장 자동 생성. **캡션·`shows` 를 생성 시점에 강제**, 건너뛴 그림에 이유 필수 |
| `[~]` | `simbot/analysis/` | `trap` `structure` **완료** (`structure` 는 캠페인이 열었다 — RDF·`ψ₆`·Voronoi·`S(k)`·시간분해·유한크기·부트스트랩, 테스트 37). 남음: `msd` `microrheo` `active` `equilibration` — **계가 생길 때.** 이 원칙이 옳았다: `structure` 는 검증 런이 생기고 나서 만들었고 그래서 버그 3개를 테스트가 잡았다 |

### 테스트 현황 — **562 통과 / 1 skip / 94 s** (`-m "not slow"` 는 14 s)
| 파일 | 대상 | 개수 |
|---|---|---|
| `test_s0_units.py` | 단위·상수·척도 왕복 | 28 |
| `test_s2_estimators.py` | 해석해 항등식·극한·비용모델 | 32 |
| `test_s3_cutoff.py` | `r_cut` 제안 (WCA/LJ/Yukawa/Morse) | 36 |
| `test_s3_spec.py` | provenance·게이트·파생값 회귀·YAML 왕복 | 47 |
| `test_s4_nondim.py` | 카드 척도·왕복<1e-12·dt 제약·정책 | 55 |
| `test_s5_forces.py` | `HarmonicTrap` 힘/에너지/수치미분 | 9 (1 skip) |
| `test_s5_guards.py` | 배위온도 + 가드 발동 | 15 |
| `test_s5_scheme.py` | B1·B2·B5·B7·B9 + 재현성 | 8 |
| **`test_s6_viz.py`** | **캡션 강제**·이중축·건너뛴 이유·독립 프레임 | **34** |
| **`test_s5_pair.py`** | **쌍 상호작용 러너**·상자 모양·Table 포텐셜 | **42** |
| **`test_s7_structure.py`** | **시간분해·유한크기 지수·부트스트랩·상 판독** | **37** |
| `test_s7_validate.py` | tolerance 파싱·판정·검정력·봉인 | 53 |
| `test_s8_io.py` | 해시·run 디렉터리·**봉인**·provenance | 31 |
| `test_s8_report.py` | REPORT.md — 나쁜 소식이 빠지지 않는가 | 33 |
| `test_cli_session.py` | 세션 append-only·예산 게이트·`session run`·전체 관통 | 54 |
| **`test_agent_layer.py`** | **`.claude/` 구조** — frontmatter·링크·티어링·권한 | **64** |

### 파이프라인 단계 구현 (§2)
| | 단계 | 코어 모듈 | 상태 |
|---|---|---|---|
| `[O]` | S1 Intake | — (LLM) | 프로토콜: `.claude/skills/bd-pipeline/references/s1_intake_drawing.md` |
| `[O]` | S2 Predict | `estimators.py` · `spec.Prediction` | 예측 YAML + 봉인. 9항목 예시 |
| `[O]` | S3 Specify | `spec.py` `units.py` `cutoff.py` | provenance 강제 · 게이트 선언 · 파생값 재계산 대조 |
| `[O]` | S4 Nondim | `nondim.py` `policy.py` | 카드별 척도 · 왕복 `1.6e-16` · dt 4제약 |
| `[O]` | S5 Run | `build.py` `forces.py` `run.py` `guards.py` | 배치 동시 실행 + 실패 런 기록 |
| `[O]` | S6 Visualize | `viz.py` | 그림 5장 + `06_figures.md`. 캡션 없는 그림은 **만들 수 없다** |
| `[O]` | S7 Validate | `analysis/trap.py` `validate.py` | PASS/FAIL/**INCONCLUSIVE** + 설계 검정력 |
| `[O]` | S8 Conclude | `report.py` `io.py` | `REPORT.md` 자동 생성. 결론 **서술**은 에이전트 |

### 물리 도메인 (§4)
| | 도메인 | 참조 케이스 | 회귀 기준 | 상태 |
|---|---|---|---|---|
| `[~]` | **A** 자유·구속 BD | `examples/trap-2d-5um/` | 등분배 `⟨x²⟩=kT/k` | **구속 분기만.** ★ **자유확산 `D*=1.00±0.03` 이 아직 없다** (Q9) |
| `[X]` | **B** 미세유변학 | — | 뉴턴유체 극한 `G''=ηω` | 카드 없음 |
| `[X]` | **C** 활성물질 (ABP) | — | ABP MSD 해석식 | 카드 `draft`, 척도 규칙 미구현 |
| `[X]` | **D** 응집 콜로이드 | — | `B₂` vs `βU_min` | 카드 없음 |
| `[O]` | **E** **2D 소프트 반발 `A/r³`** ★신규 | `runs/2026-07-29_soft-r3-*` (6런) | 완벽격자 `ψ₆=1`·액체 지수 `p=1/2` | **가장 발전된 도메인.** 벤치마크 31개 · §4.1 |

> ★ **도메인 E 가 이 저장소 최대 지식체다** — 6런 · 3,856 프로덕션 런 ·
> 벤치마크 31개 · 방법론 finding 6개. 전문은 **§4.1**.

> ★ **도메인 A 의 기본 검증이 비어 있다.** 트랩 분기는 등분배·완화시간·EM 편향까지
> 검증됐지만, **자유확산 `D* = 1`** 은 한 번도 측정하지 않았다. M1 의 원래 DoD 였고
> 손그림이 트랩이어서 건너뛰었다 (§8). 그 결과:
> - `CARD_SCALE_RULES` 의 `brownian` 경로(`σ`, `τ_D`)가 **end-to-end 로 실행된 적 없다**
> - `dt` **변위 게이트가 실제로 구속한 런이 0건**이다 (트랩 계에서는 꺼진다)
> - `analysis/msd.py` 를 만들 검증 대상이 없다
>
> **정정 (2026-07-28, 위 3항 재검증):** 둘째 항은 **파이프라인에 한해서만** 참이다.
> `scripts/chain_bend.py:113` 이 같은 임계값(`0.03`/`0.005`)을 재구현한 경로에서
> 변위 게이트가 **이미 구속했다** — `runs/chain-bend/smoke/batch.json`:
> `binding: "force"`, `dt_force = 4.82e-6` vs `dt_diffusion = 4.5e-4`,
> 실측 `max|F*| = 1037.7` (힘 제약이 `dt` 를 100배 줄였다).
> ⇒ `simbot.nondim.choose_dt` 는 여전히 구속 0건이고, **게이트 로직이 두 곳에 복제됐다.**
> 그리고 `choose_dt` 의 변위 제약은 `active=has_pair`, `has_pair = bool(spec.pair)` 인데
> `SystemSpec` 에 bond·angle 필드가 없다 → **결합만 있는 계를 파이프라인에 올리면
> 게이트가 조용히 꺼진다.** 위 실측이 바로 그 게이트가 필요한 계에서 나왔다.
> (첫째·셋째 항은 재검증에서 그대로 참: `scales_brownian` 호출자는 `nondim.py:68` 과
> `test_s0_units.py` 뿐 · `simbot/analysis/` 에는 `trap.py` 하나뿐)

### 인풋 모달리티 (§8 M5)
| | 모달리티 | 상태 |
|---|---|---|
| `[O]` | **손그림 사진** ← v1 목표 | ✅ **완주.** 판독 프로토콜 + 실제 사용자 그림 1장 |
| `[~]` | 텍스트 설명 | spec YAML 을 직접 쓰면 동작. 자연어 → spec 경로는 없음 |
| `[X]` | 실험 화면/영상 | M5. **A2(`a` vs `R`)를 측정으로 닫으려면 필요** |
| `[X]` | 논문 PDF | M5 (`bd-lit-distill` 이 증류는 한다) |
| `[X]` | 음성 | M5 |

### 추가 제안 채택 여부 (§9) — 상세는 §9
| | 제안 | 결정 |
|---|---|---|
| `[O]` | 1. git 초기화 | ✅ 커밋 28개 |
| `[O]` | 2. **예측 봉인** | ✅ `io.py`. `shasum -c` 호환 · 깨지면 대조표를 만들지 않는다 |
| `[X]` | 3. 파일럿 런 | 미구현. 정책에 `mandatory: true` 로 선언됐으나 `cli.py` 가 실행하지 않는다 |
| `[O]` | 4. 단위 접미사 강제 (`_si`/`_star`) | ✅ `Quantity.si` 가 문자열·bool 거부 + 왕복 테스트 |
| `[~]` | 5. 질문 예산 | 규약은 `CLAUDE.md`·스킬에. **코드 강제 장치는 없다** |
| `[O]` | 6. **`bd-diagnose` 스킬** | ✅ `.claude/skills/bd-diagnose/` |
| `[~]` | 7. 파라미터 스윕 (`sweep: [...]`) | **`spec.yaml` 에는 아직 없다.** 그러나 `scripts/soft2d_*.py` 6개가 `A`·`N`·시드 스윕을 실제로 돌렸다 (3,856 런) — **패턴이 확립됐으므로 이제 일반화할 재료가 있다** |
| `[~]` | 8. run 캐시 | 재료는 있다 (`spec.hash()`, `completed_stages()`). **`spec_hash` 조회는 없다** |
| `[X]` | 9. HTML 리포트 | 미구현. 그림이 생겼으므로 이제 이득이 있다 |
| `[X]` | 10. **손그림 작성 가이드** | 미작성. 첫 그림 모호성 2건 + **`soft-r3` 의 `r` 단위 공백** (Zahn 규약 `d = n^{-1/2}` 로 해석해야 `Γ = π^{3/2}A` 가 성립 — 다르게 읽으면 물리가 달라진다). 근거가 늘었다 |
| `[-]` | 11. Langevin 폴백 | 검토 후 결정. `overdamped` 게이트가 위반을 잡고 권고만 한다 |
| `[-]` | 12. HI 근사 (RPY) | v1 비범위 |
| `[-]` | 13. 실험 데이터 직접 대조 | 검토 후 결정. A2 를 닫으려면 필요 |

---

## 0. 목표와 범위

### 0.1 한 줄 정의
사용자가 제공한 자료(v1: **손그림**)를 해석해 Brownian Dynamics 시뮬레이션을
**설계 → 예측 → 실행 → 검증 → 결론**까지 자율 수행하고, 그 과정에서 얻은 판단 근거를
지식 베이스에 축적하는 Claude Code 네이티브 에이전트.

### 0.2 확정된 설계 결정
| 항목 | 결정 | 비고 |
|---|---|---|
| 챗봇 런타임 | **Claude Code 네이티브** | 이 대화창이 곧 챗봇. API 키 불필요 |
| 물리 엔진 | **HOOMD-Blue 7.1.0** | `md.methods.Brownian` (overdamped Langevin) |
| 실행 환경 | **로컬 CPU 단독** | Apple Silicon, GPU 없음. N ≲ 10⁴ 스케일 |
| conda env | **`simulation_bot`** (신규) | 패키지 단계적 적립, `knowledge/env_log.md`에 이력 |
| v1 인풋 | **손그림 사진** | 이후 실험영상 → 논문 PDF → 음성 → 텍스트로 확장 |
| 물리 도메인 | 4개 전부 (§4) | 자유/구속 BD · 미세유변학 · 활성물질 · 응집콜로이드 |

### 0.3 명시적 비범위 (v1에서 안 함)
- 유체동역학 상호작용(HI): Oseen/RPY/Stokesian dynamics **없음**. 자유배수(free-draining) 근사.
  → 이 근사가 깨지는 조건은 `knowledge/models/no_hydrodynamics.md`에 명시.
- GPU / MPI / 클러스터 제출
- 반응(화학), 유동장 결합(CFD), 전자기 완전결합

### 0.4 설계 원칙 (전 단계 관통)
1. **에이전트는 판단하고, 코어는 계산한다.**
   LLM이 숫자를 "머리로" 계산하는 일은 금지. 모든 수치는 `simbot/` 함수를 호출해서 얻는다.
   LLM의 역할은 *어떤 모델을 쓸지, 어떤 값을 가정할지, 왜 그런지*를 결정하고 기록하는 것.
2. **예측을 먼저 봉인한다.**
   S2 예측 문서는 S5 실행 전에 해시로 봉인. 사후합리화(post-hoc rationalization)를 구조적으로 차단.
3. **모든 숫자에 출처(provenance)가 있다.**
   `from_drawing` / `from_paper` / `from_knowledge` / `assumed` / `derived` 중 하나.
4. **단위를 타입으로 취급한다.**
   변수명 접미사 `_si`(물리 단위) vs `_star`(무차원)를 강제. 혼용은 테스트로 잡는다.
5. **실패도 산출물이다.**
   폭발한 시뮬레이션, 틀린 예측, 잘못된 그림 해석 → 전부 `knowledge/failures/`에 기록.
6. **재현성이 기본값이다.**
   run 하나가 `spec + seed + env + 코드 해시`만으로 완전 복원 가능해야 한다.

---

## 1. 시스템 아키텍처

### 1.1 4개 레이어

```
┌──────────────────────────────────────────────────────────────┐
│ L1  에이전트 레이어 —  .claude/skills/ , CLAUDE.md            │
│     Claude Code가 직접 수행. 해석·판단·리즈닝·질문·기록.       │
│     결정론적 계산은 절대 하지 않고 L2를 호출한다.              │
└───────────────┬──────────────────────────────────────────────┘
                │ 호출 (python -m simbot.… / import)
┌───────────────▼──────────────────────────────────────────────┐
│ L2  결정적 코어 —  simbot/                                    │
│     순수 Python. LLM 없음. pytest로 전수 검증.                │
│     단위변환·무차원화·해석해·HOOMD 실행·분석·플롯·리포트.       │
└───────────────┬──────────────────────────────────────────────┘
                │ 읽기/쓰기
┌───────────────▼──────────────┐  ┌───────────────────────────┐
│ L3  지식 레이어 — knowledge/  │  │ L4  산출물 — runs/<run_id>/│
│     시스템 아키타입, 파라미터  │  │     단계별 아티팩트 전체,   │
│     근거, 모델링 근거, 실패    │  │     궤적, 그림, 리포트.     │
│     사례, 검증 벤치마크.       │  │     자기완결·재현가능.      │
│     ← 누적, 버전관리 대상      │  │     ← .gitignore 대상      │
└──────────────────────────────┘  └───────────────────────────┘
```

**L3가 이 프로젝트의 진짜 자산이다.** L2 코드는 다시 쓸 수 있지만,
"왜 이 물질에 η=1.2 mPa·s를 썼는가", "왜 dt=5e-5에서 터졌는가"는 축적하지 않으면 사라진다.

### 1.2 디렉터리 구조

> 표기: `[O]` 존재하고 테스트됨 · `[X]` 계획만 · `[~]` 부분 구현
> (2026-07-28 실측. 이 절은 실제 트리와 일치해야 한다 — 어긋나면 설계 문서가 아니라 소설이다)

```
Simulation_bot/
├── master_plan.md            [O] 이 문서. 전체 설계의 단일 진실 원천
├── CLAUDE.md                 [O] 에이전트가 매 세션 읽는 프로젝트 규약
├── environment.yml           [O] simulation_bot env 재현
├── pyproject.toml            [O] pytest 설정 (마커·경로)
├── cli.py                    [O] run · resume · converge · params · calibrate
├── README.md                 [X] 사람용 사용법 — 아직 없다
│
├── .claude/                  [O] L1 에이전트 층. §12.3–12.4
│   ├── README.md                 [O] 구성 + 티어링 표 + Q6 결정 근거
│   ├── skills/
│   │   ├── bd-pipeline/          [O] [메인] S1→S8 오케스트레이터
│   │   │   ├── SKILL.md              단계·게이트·금지사항 체크리스트
│   │   │   └── references/       [O] s1_intake_drawing ★ · s2_prediction
│   │   │                             s3_s5_execute · s6_s7_validate · s8_knowledge
│   │   ├── bd-diagnose/          [O] 터진 런 진단 (배제 순서)
│   │   └── bd-knowledge/         [O] knowledge/ 검색·추가·정리
│   ├── agents/                   [O] 서브에이전트 9개. model: frontmatter 로 티어링
│   └── settings.json             [O] 인터프리터 허용 + **봉인 문서 편집 거부**
│
├── simbot/                   ← L2 결정적 코어. LLM 0줄
│   ├── units.py              [O] 물리상수, Scales, 카드별 척도 팩토리
│   ├── spec.py               [O] Quantity/SystemSpec/Prediction + 게이트 검사
│   ├── nondim.py             [O] 카드별 무차원화, 무차원수, dt 4제약
│   ├── policy.py             [O] run_policy.yaml 로더 + overrides 깊은 병합
│   ├── estimators.py         [O] 해석해·스케일링 (S2 예측 엔진)
│   ├── cutoff.py             [O] r_cut 제안 (WCA/LJ/Yukawa/Morse)
│   ├── build.py              [~] 트랩 스냅샷만. 겹침 제거는 아직 없음
│   ├── forces.py             [~] HarmonicTrap (md.force.Custom) 만
│   ├── run.py                [O] 트랩 BD 러너 + 배치 동시 실행
│   ├── guards.py             [O] NaN/변위/배위온도/요동 검사
│   ├── analysis/
│   │   ├── trap.py           [O] MSD 피팅, 시드 집계, 분포 검정
│   │   ├── msd.py            [X] 일반 MSD, D, 블록평균 오차
│   │   ├── structure.py      [X] RDF, S(q), 클러스터, 밀도프로파일
│   │   ├── microrheo.py      [X] GSER → G'(ω), G''(ω)
│   │   ├── active.py         [X] ABP MSD crossover, MIPS 판정
│   │   └── equilibration.py  [X] 평형 도달 판정
│   ├── validate.py           [O] 예측 vs 측정, PASS/FAIL/INCONCLUSIVE, 검정력
│   ├── report.py             [O] REPORT.md 생성
│   ├── session.py            [O] 세션 상태 (턴 append-only), set = 추정만
│   ├── io.py                 [O] run 디렉터리, 해시, **봉인**
│   └── viz.py                [X] ★ 없음. S6 그림은 scripts/ 에 일회성으로
│
├── scripts/                  [~] 일회성 — simbot 으로 옮겨야 할 코드
│   ├── trap_batch.py             → run.run_trap_batch 로 흡수됨 (중복)
│   └── trap_analyze.py           → viz.py 가 없어서 아직 여기 산다
│
├── examples/                 [O] 기계가 읽는 참조 케이스
│   └── trap-2d-5um/
│       ├── spec.yaml             S3 명세 (provenance 18필드)
│       └── prediction.yaml       S2 예측 9항목 (봉인 대상)
│
├── knowledge/                ← L3 지식 베이스. 스키마는 §10 (BD_agent 정본)
│   ├── source/papers/        [O] 문헌 증류 42편 + INDEX
│   └── wiki/                 [O] CLAUDE.md 계약 + systems 5 · findings 13
│       │                         concepts 3 · benchmarks 2 · questions 2
│       ├── systems/              ★ (계 × 목적동역학) 카드 — 무차원화·게이트 소유
│       ├── findings/             Q→A + dead-end
│       ├── concepts/  techniques/  benchmarks/  questions/
│
├── inputs/                   [O] 사용자 제공 원자료 (gitignore, .sha256 만 추적)
│   └── <topic>/…
│
├── sessions/                 [O] 세션 상태 (gitignore)
│   └── <session_id>/session.yaml + spec_turnNN.yaml
│
├── config/run_policy.yaml    [O] 자원·티어·dt 정책. 사람의 overrides: 우선
│
├── tests/                    [O] pytest — 373 통과 / 1 skip / 28 s
│   └── test_s0_units · test_s2_estimators · test_s3_cutoff · test_s3_spec
│       test_s4_nondim · test_s5_{forces,guards,scheme} · test_s7_validate
│       test_s8_{io,report} · test_cli_session
│
└── runs/                     [O] 산출물 (gitignore, 단 .md/.json/.yaml 추적)
    └── 2026-07-28_trap-2d-5um_2dfb9d/     첫 손그림 완주 (사람+에이전트)
        2026-07-28_cli-e2e-test/           같은 계를 cli.py 로 재현 (2.3 s)
```

### 1.3 run 디렉터리 규약 (자기완결성)

```
runs/<run_id>/
├── 00_input/            원자료 사본 (손그림 사진 등) + sha256
├── 01_intake.md         관찰/추론/가정 분리 기록
├── 01_intake.json
├── 02_prediction.md     ⚠ 봉인됨. S5 이후 수정 금지
├── 02_prediction.json
├── 03_spec.yaml         물리 단위 완전 명세
├── 03_spec_rationale.md 각 값의 출처와 근거
├── 04_reduced.yaml      무차원 명세 + 역변환 계수
├── 04_nondim.md         변환표
├── 05_run_manifest.json 코드해시·env해시·seed·HOOMD버전·wall time
├── traj.gsd             궤적
├── thermo.h5            열역학 로그
├── 06_figures.md        그림 목록 + 캡션
├── figs/*.png
├── 07_validation.md     예측 vs 측정 대조표
├── metrics.json         측정값 + 오차
├── 08_conclusion.md     결론
└── REPORT.md            전체 요약 (사람이 읽는 최종 산출물)
```

`run_id = <ISO시각>_<슬러그>_<spec해시 앞6자리>`

---

## 2. 8단계 파이프라인 — 상세

각 단계는 **입력 → 처리 → 산출물 → 게이트(통과조건) → 실패모드**로 정의된다.
게이트를 통과하지 못하면 다음 단계로 넘어가지 않고, 사용자에게 보고하거나 이전 단계로 되돌아간다.

---

### S1. Intake — 자료 해석

**입력** `inputs/<topic>/` 의 손그림 사진 (+ 사용자의 구두 설명)

**처리**
1. 이미지를 읽고 다음을 목록화한다:
   - 입자: 개수(정확/대략), 크기 차이, 종류(색·해칭·라벨로 구분), 특별 표시된 개체(프로브 등)
   - 경계: 상자 테두리, 벽(실선), 주기경계(점선/화살표), 슬릿/원통/구 형상, 차원(2D/3D)
   - 화살표: 위치·방향·길이 → **힘 / 속도 / 흐름 / 시간진행 중 무엇인지 후보 나열**
   - 텍스트: 숫자, 단위, 기호(η, T, k, φ, v₀…), 축 라벨, 캡션
   - 그래프: 손으로 그린 예상 곡선(있으면 S2 예측과 대조할 근거)
2. **관찰/추론/가정 3단 분리** — 이것이 S1의 핵심 산출물이다.

   | 등급 | 정의 | 예 |
   |---|---|---|
   | `observation` | 그림에서 직접 읽음 | "입자 약 30개", "왼쪽·오른쪽에 실선 벽" |
   | `inference` | 그림 + 물리지식으로 유도 | "실선 벽 + 상하 점선 → 슬릿 기하, y·z 주기" |
   | `assumption` | 그림에 없어 내가 채움 | "매질은 물, η=1.0 mPa·s, T=298 K" |

   각 항목에 `confidence: high/medium/low`와 한 줄 근거를 붙인다.
3. **gaps** — 시뮬레이션에 필수인데 없는 정보를 나열하고 처리방침 결정:
   `ask_user` (질문 예산 내) / `fill_from_knowledge` / `assume_and_flag` / `sweep` (파라미터 스캔)
4. `knowledge/systems/` 검색 → 유사 아키타입이 있으면 재사용.

**손그림 특화 규칙** (`references/s1_intake_drawing.md`)
- 손그림의 **절대 크기는 신뢰하지 않는다.** 신뢰하는 것은 ① 토폴로지(무엇이 무엇 안에/옆에)
  ② 비율(입자:상자 ≈ 1:20) ③ 개수 ④ 대칭성 ⑤ 명시된 숫자·단위.
- 화살표 굵기/길이의 절대값은 무의미. 상대 비교만 사용.
- 모호한 요소는 임의 해석 대신 **후보 2~3개를 명시**하고 S2에서 각 후보의 결과 차이를 예측.
  → 어느 해석이 맞는지 사용자가 즉시 판별할 수 있게 된다.
- 그림에 스케일 정보가 전무하면 φ(부피분율)를 자유 파라미터로 두고 sweep 후보로 표시.

**산출물** `01_intake.md`, `01_intake.json`

**게이트**
- 필수 필드 확정: 공간차원 `d`, 입자 종류 수, 경계조건, 구동/활성 유무, 질문(question)이 무엇인지
- `question` 필드가 반증 가능한 형태인가 (예: "확산이 얼마나 느려지는가" ✅ / "어떻게 되는가" ❌)

**실패모드** → `knowledge/failures/intake_*.md`
- 화살표를 힘으로 읽었는데 실제론 속도장이었음
- 2D 그림을 2D 시뮬레이션으로 읽었는데 실제론 3D의 단면
- 손그림 입자 개수를 실제 N으로 착각 (그림은 스케치, N은 통계적으로 필요한 수)

---

### S2. Predict — 예상 결과 리즈닝

> **시뮬레이션 전에 답을 적는다.** 이 단계가 이 프로젝트의 과학적 정직성을 담보한다.

**입력** `01_intake.json`, `knowledge/`

**처리**
1. **지배 물리 식별** — 어떤 힘/시간척도가 경쟁하는가. 무엇이 결과를 결정하는가.
2. **무차원수 자릿수 추정** — φ, Pe, κσ, T*=kT/ε, k σ²/kT, D_r τ_B … (`simbot.nondim`)
   각 무차원수가 어느 레짐에 있는지, 레짐 경계에서 얼마나 떨어져 있는지.
3. **정량 예측** — `simbot.estimators` 의 해석해/스케일링 호출. 예:
   - Stokes–Einstein `D₀ = k_BT / 6πηa`
   - 자유확산 `⟨Δr²⟩ = 2d D t`
   - 광집게 `⟨Δr²⟩(t) = (2d k_BT/k)(1−e^{−t/τ_k})`, `τ_k = γ/k`, 등분배 `⟨x²⟩ = k_BT/k`
   - 침강 `ρ(z) ∝ e^{−z/ℓ_g}`, `ℓ_g = k_BT/(Δρ V g)`
   - ABP `⟨Δr²⟩ = 2dD_t t + (2v₀²/λ²)(λt − 1 + e^{−λt})`, `λ = (d−1)D_r`
     `D_eff = D_t + v₀²/[d(d−1)D_r]`
   - 농후계 장시간 확산 감소, MIPS 상경계 등은 문헌 상관식 (출처 명기)
4. **반증 가능한 형태로 봉인** — 각 예측은 다음 4요소를 반드시 가진다:

   ```yaml
   - quantity: D_long / D_0
     value: 0.42
     tolerance: "±25%"          # 이 밖이면 FAIL
     basis: "Batchelor + φ=0.35 준희박 보정, knowledge/validation/dense_diffusion.md"
     discriminates: "HI 무시가 타당한지 여부를 가른다"
   ```
5. **대안 시나리오** — 예측이 틀릴 수 있는 방식과, 그때 나타날 신호를 미리 적는다.
   (예: "dt가 너무 크면 D가 과대평가된다 → dt 절반 재실행에서 D가 바뀌면 수치 문제")

**산출물** `02_prediction.md`, `02_prediction.json`
→ 두 파일의 sha256을 `05_run_manifest.json`에 기록해 **봉인**. S7은 이 해시를 검증한다.

**게이트**
- 정량 예측 ≥ 1개, 각각에 `tolerance`와 판정기준 존재
- 그림에 사용자가 그린 예상 곡선이 있으면 그것과의 정합/불일치를 명시

**실패모드** → `knowledge/failures/prediction_*.md`
- 무차원수 자릿수를 틀려 레짐을 오판
- 문헌 상관식을 적용범위 밖에 적용
- tolerance를 너무 넓게 잡아 어떤 결과든 PASS (검증 무력화) — **금지, 리뷰 대상**

---

### S3. Specify — 시스템 구체화 (물리 단위)

**입력** `01_intake.json`, `knowledge/parameters/`, `knowledge/systems/`

**처리** 완전한 `SystemSpec`(SI 단위)을 채운다. 빈 필드가 남아선 안 된다.

| 그룹 | 필드 |
|---|---|
| 기하 | `dim`, `box_si` (Lx,Ly,Lz), `boundary` (pbc/wall/slit/cylinder/sphere) |
| 입자 | `species[]`: `name, N, radius_si, mass_si, density_si, charge, active` |
| 매질 | `T_si` (K), `eta_si` (Pa·s), `rho_fluid_si`, `epsilon_r` |
| 마찰 | `gamma_si` = 6πηa (Stokes) 또는 명시값; 벽 근접 보정 여부 |
| 상호작용 | `pair[]`: 타입쌍별 포텐셜(WCA/LJ/Yukawa/Morse/DLVO) + 파라미터 + `r_cut` |
| 외장 | `traps[]`, `gravity`, `shear_rate`, `field` |
| 활성 | `v0_si`, `D_r_si` 또는 `tau_r_si` |
| 시간 | `t_total_si`, `t_equil_si`, `dump_interval_si`, `thermo_interval_si` |
| 수치 | `seed`, `dt_policy` |

**모든 필드에 provenance 필수:**
```yaml
eta_si:
  value: 0.890e-3          # 298.15 K. 20 °C 값(1.002e-3)과 혼동 금지
  unit: "Pa*s"
  provenance: from_knowledge
  source: "knowledge/parameters/water_298k.md"
  note: "298.15 K 순수 물, IAPWS"
```

**물리적 타당성 자동 체크** (`simbot.spec.validate`)
- φ < 0.64 (RCP), 2D면 φ_A < 0.9
- Reynolds `Re = ρ v a/η ≪ 1` (BD 전제)
- 관성 시간척도 `τ_i = m/γ ≪ dt` (overdamped 전제) — 위반 시 Langevin 권고
- Debye 길이 vs 입자간 거리 정합
- `r_cut` < L/2 (최소이미지)
- 활성계: `v₀ τ_r` (persistence length) vs box 크기

**산출물** `03_spec.yaml`, `03_spec_rationale.md`

**게이트** 모든 필드 채워짐 + provenance 존재 + 타당성 체크 전항 통과(또는 명시적 예외 승인)

**실패모드** → `knowledge/failures/spec_*.md`
- γ를 6πηa 대신 6πηd로 계산 (반지름/직경 혼동) — **가장 흔한 실수**
- 2D 시뮬레이션에서 3D Stokes 마찰 사용 (의도적이면 명시)
- box가 persistence length보다 작아 유한크기 인공물

---

### S4. Nondimensionalize — 무차원화

> 상세 규약은 §5. 여기서는 파이프라인 관점만.

**입력** `03_spec.yaml`

**처리**
1. 기준 척도 3개 선택 (기본: `L*=σ`, `E*=k_BT`, `T*=τ_B=σ²/D₀`)
2. 모든 SI 값 → 무차원 값 변환, 역변환 계수 저장
3. 무차원수 전량 계산 및 표로 정리
4. **dt 선택** — §5.4의 4개 제약 중 최소값 채택, 근거 기록
5. 왕복 변환 테스트: `to_reduced → to_si` 상대오차 < 1e-12

**산출물** `04_reduced.yaml`, `04_nondim.md`(변환표: 물리량 | SI | 무차원 | 역변환계수)

**게이트** 왕복오차 통과 + dt 제약 전항 만족 + 무차원수 표 완성

**실패모드** → `knowledge/failures/nondim_*.md`
- HOOMD 시간 단위(τ_LJ = σ√(m/ε))와 Brownian 시간(τ_B = σ²/D₀) 혼동 → **치명적**
- 2D에서 τ_B 정의의 차원 계수(2d) 누락
- kT와 ε를 둘 다 1로 잡아 T*=1로 고정해버림 (의도 아니면 버그)

---

### S5. Run — HOOMD-Blue 실행

**입력** `04_reduced.yaml`

**처리**
1. **초기배치** (`simbot.build`) — 격자/랜덤/비중첩 삽입. 겹침 있으면 soft pushoff(점진적 σ 증가)로 제거.
2. **파일럿 런** (기본 활성) — 본실행의 0.5% 스텝만 돌려서
   ① 가드 위반 없는지 ② 예상 wall time ③ dt 타당성 확인. 문제 있으면 여기서 중단.
3. **포스/적분기** (`simbot.forces`, `simbot.run`)
   - `md.methods.Brownian(filter, kT, default_gamma)` — D = kT/γ
   - 활성: `md.force.Active` + `md.update.ActiveRotationalDiffusion`
   - 벽: `md.external.wall.{LJ,Morse,Yukawa,ForceShiftedLJ}` + `hoomd.wall.{Plane,Sphere,Cylinder}`
   - 광집게/조화구속: `md.force.Custom` 서브클래스 (내장 없음)
   - 중력/전기장: `md.force.Constant` 또는 `md.external.field.Electric`
4. **평형화 → 프로덕션** 분리. 평형 판정은 `simbot.analysis.equilibration`.
5. **런타임 가드** (`simbot.guards`) — 매 `thermo_interval`마다 검사, 위반 시 즉시 중단 + 진단 저장:
   - NaN/Inf 위치·힘
   - 스텝당 최대 변위 > 0.1σ (폭발 징후).
     ⚠ **변위 분포를 Gaussian으로 가정하지 말 것** — HOOMD 노이즈는 균일분포다 (`max/σ = √3`)
   - 벽 관통 (경계 밖 입자)
   - **배위 온도** `kT_conf = ⟨|∇U|²⟩/⟨∇²U⟩` 가 입력 `kT`와 일치하는가
     — 실측 `1.00382 ± 0.00480` (조화 트랩). 힘이 이미 계산되어 있어 추가 비용 거의 없음
   - 압력/포텐셜 에너지 발산

   > ❌ **삭제됨: "운동에너지 온도가 목표 kT에서 이탈"** — 작동할 수 없는 가드였다.
   > HOOMD `Brownian`은 속도를 적분하지 않고 **매 스텝 목표 분포에서 뽑는다.**
   > `kinetic_temperature`는 입력값을 되풀어 말할 뿐 계통적으로 이탈하지 못한다.
   > 근거: [`findings/hoomd-brownian-scheme-and-noise.md`](../../knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md)
   >
   > **`kT` 자체는 여전히 1급 입력이다** — `U/kT`가 Boltzmann 가중치를 정하고 모든
   > 무차원수(`k* = kσ²/kT`, `T* = kT/ε`, `Pe = Fσ/kT`)에 들어간다. 쓸 수 없는 것은
   > 운동에너지 되읽기뿐이고, 그 자리를 배위온도가 대신한다.
6. **로깅** — GSD 궤적(위치·방향·이미지 플래그), `ThermodynamicQuantities` → HDF5.
   재현성: seed, HOOMD 버전, `simbot` 코드 해시, env 해시, spec 해시, wall time.

**산출물** `traj.gsd`, `thermo.h5`, `05_run_manifest.json`

**게이트** 완주 + 가드 무위반 + 평형 도달 + 프로덕션 구간 길이 ≥ 요구 상관시간의 10배

**실패모드** → `knowledge/failures/run_*.md`
- ABP에서 `moment_inertia`를 0으로 안 두고 `Brownian`을 써서 회전 자유도가 이중 적분됨
- `r_cut` > L/2 로 최소이미지 위반
- dt가 커서 WCA 코어를 관통 (입자 겹침 후 폭발)
- 초기 겹침을 제거하지 않고 시작 → 1스텝에 발산

---

### S6. Visualize — 시각화

**입력** `traj.gsd`, `thermo.h5`

**처리**
1. **필수 진단 세트** (도메인 무관, 항상 생성)
   - MSD log–log + 기울기 1 참조선
   - RDF g(r)
   - thermo 시계열 (PE, P, 온도추정)
   - 스텝당 변위 분포 (가드 사후 확인)
   - 초기/중간/최종 스냅샷 3연
2. **도메인별 플롯** — §4 표 참조
3. **애니메이션** — 2D는 matplotlib, 3D는 fresnel 레이트레이싱. GIF/MP4.
   (ffmpeg 미설치 → GIF 우선, 필요 시 ffmpeg 설치 후 MP4)
4. **모든 축에 무차원 값과 물리 단위 이중 표기** (`t/τ_B` 위, `t [ms]` 아래)

**산출물** `figs/*.png`, `figs/*.gif`, `06_figures.md`(각 그림 캡션 = 무엇을 보이려는 그림인가)

**게이트** 축 라벨·단위·캡션 전부 존재. 캡션 없는 그림은 산출물로 인정하지 않는다.

---

### S7. Analyze & Validate — 분석 및 검증

**입력** `traj.gsd`, `02_prediction.json`(봉인 해시 검증), `04_reduced.yaml`

**처리**
1. **봉인 검증** — `02_prediction.*`의 sha256이 manifest와 일치하는지 확인. 불일치면 **중단하고 보고**.
2. **정량 측정** (`simbot.analysis`) — 각 측정값에 통계오차 동반
   - MSD 피팅 → D (블록평균 + bootstrap 신뢰구간)
   - 구조: g(r), S(q), 클러스터 크기 분포, 밀도 프로파일
   - 미세유변학: GSER(Mason) → G'(ω), G''(ω)
   - 활성: MSD crossover 시각, D_eff, MIPS 판정지표
3. **예측 vs 측정 대조표** — 각 항목 `PASS` / `FAIL` / `INCONCLUSIVE`
   `INCONCLUSIVE`는 통계오차가 tolerance보다 커서 판정 불가인 경우 → 더 긴 런 필요
4. **수치 수렴 체크**
   - dt 래더: dt/2로 재실행해 측정값 변화 < 통계오차인지 (기본: 핵심 지표 1개에 대해 수행)
   - 유한크기: N 또는 L을 √2배 해서 재실행 (요청 시)
   - 유한시간: 궤적 전반부/후반부 측정값 일치 여부
5. **알려진 해석해 대조 (sanity)** — 항상 실행
   - 무차원 자유확산이면 `D* = 1.00 ± stat`
   - 조화구속이면 등분배 `⟨x*²⟩ = 1/k*`
   - ABP면 해석 MSD 곡선과 중첩
6. **FAIL 항목마다 원인 가설** — 4분류: `numerical` / `modeling` / `interpretation` / `analysis`

**산출물** `07_validation.md`, `metrics.json`

**게이트** sanity 체크 통과 + 모든 FAIL에 원인 가설 존재

---

### S8. Conclude — 결론 요약 + 지식 커밋

**입력** 전 단계 산출물

**처리**
1. **질문 → 답** — S1의 `question`에 직접 답한다. 한 문단.
2. **근거 3줄** — 어떤 측정이 그 답을 지지하는가.
3. **신뢰도와 한계** — HI 무시, 유한 N/t, dt 수렴, 가정의 영향.
4. **다음 실험 제안** — 이 결론을 더 굳히거나 반증할 최소 비용 실험 1~2개.
5. **knowledge/ 업데이트 (필수, 생략 불가)**
   - `systems/`: 이 시스템 아키타입 요약 (또는 기존 항목 갱신)
   - `parameters/`: 새로 확정한 파라미터 근거
   - `models/`: 모델링 결정 근거 (왜 이 포텐셜/근사인가)
   - `failures/`: 발생한 모든 실패 (사소한 것도)
   - `validation/`: 이번에 확인된 벤치마크 수치 → pytest 회귀 테스트로 승격 후보

**산출물** `08_conclusion.md`, `REPORT.md`, knowledge 항목 N개

**게이트** knowledge에 최소 1개 항목 추가/갱신됨. `REPORT.md`가 그림/수치 링크 포함해 자기완결적.

---

## 3. 데이터 모델

**구현 완료** (2026-07-28). dataclass 기반, pydantic 없음. 설계와 달라진 곳은 ★ 로 표시.

```python
# simbot/spec.py
Quantity:        value, unit, provenance, basis, confidence, ambiguity,
                 sensitivity, affects[], written_by
                 ★ Provenance 를 별도 dataclass 로 두지 않았다 — 필드가 5개뿐이라
                   래퍼를 한 겹 더 두면 YAML 이 두 단 깊어지고 손으로 못 쓴다
Species:         name, n_simulated, radius_si, density_si, n_physical, charge, active
Geometry:        dim, boundary, box_si | box_over_ref
                 ★ box_over_ref 추가 — 트랩 계는 박스를 ℓ_trap 배수로 준다
Medium:          T_si, eta_si, rho_fluid_si, species
Friction:        model, gamma_si, wall_correction, note
PairInteraction: type_a, type_b, potential, params, r_cut_si
ExternalField:   kind, params{Quantity}, implementation, note
                 ★ traps/gravity/shear/active 를 개별 필드로 두지 않고 하나로 합쳤다 —
                   외장 종류마다 필드를 늘리면 스키마가 계마다 달라진다
Timing:          equil_in_tau, prod_in_tau, sample_interval_in_tau, target_precision
Numerics:        dt_star, seed_base, n_seeds, integrator, scheme, noise_distribution
Gate:            status(required|pass|fail|off|applicable|unknown), reason
                 ★ 신규 — 게이트를 카드가 켜고 끈다. off 에는 이유가 필수
SystemSpec:      card, question, geometry, species[], medium, friction, pair[],
                 external[], timing, numerics, gates{}, tier, notes[]
PredictionItem:  quantity, value, tolerance, basis, discriminates, unit,
                 competing_value, note        ★ competing_value 추가 — 검정력 계산
Prediction:      items[], regimes, alternatives[]
                 ★ sealed_hash 를 여기 두지 않았다 — 봉인은 파일 해시이므로
                   문서 안에 자기 해시를 넣을 수 없다. io.write_seal 이 소유

# simbot/nondim.py
Scales:          length_si, energy_si, time_si, origin     (units.py)
DtConstraint:    name, dt_si_max, active, basis, off_reason     ★ 신규
DtChoice:        dt_si, dt_star, dominant, constraints[], logged{}   ★ 신규
ReducedSpec:     card, scales, dim, n_particles, box_star, kT_star, gamma_star,
                 D_star, sigma_star, dt_star, dt_dominant, k_star,
                 equil_steps, prod_steps, sample_interval_steps, groups{}, logged{}
                 ★ sigma_star 는 1 이 아니다 — 트랩 카드에서 491.358 이다

# simbot/validate.py
Tolerance:       kind(relative|absolute|lower_bound|upper_bound), magnitude, text
Measurement:     quantity, value, stat_err, method, n_samples, spread, unit
ValidationRow:   quantity, predicted, measured, tolerance, verdict, stat_err,
                 deviation, deviation_rel, sigma, design_power,
                 samples_needed_for_3sigma, cause_class, reason, note, flags[]
                 ★ design_power·samples_needed_for_3sigma·flags 추가

# simbot/io.py
RunDir:          path + RUN_LAYOUT 키로 접근
SealVerdict:     ok, changed[], missing[], unsealed[], entries{}   ★ 신규
RunManifest:     ★ dataclass 로 두지 않고 build_manifest() → dict.
                 필드가 env 버전표를 포함해 가변이라 고정 스키마가 방해된다
```

**불변식** (테스트로 강제)
- `*_si` 필드는 항상 `Quantity`(단위 있음). `*_star` 필드는 항상 순수 float.
  → `Quantity.si` 가 문자열·bool 을 거부한다 (`test_s3_spec.py`)
- `SystemSpec` → YAML → `SystemSpec` 왕복 오차 **0** (파생값 비트 일치)
- `SystemSpec` → `ReducedSpec` → SI 왕복 오차 **< 1e-12** (실측 `1.6e-16`)
- 모든 `Quantity`에 `provenance` 존재. **기본값이어도 파일에 적힌다** —
  `assumed` 를 생략하면 "가정했다"와 "적기를 잊었다"가 구별되지 않는다
- `provenance ∈ {inference, assumed}` 인 필드는 저가 모델이 쓸 수 없다 (§12.2)

---

## 4. 물리 도메인 커버리지

| # | 도메인 | 모델 | HOOMD 구성 | 핵심 무차원수 | 검증 기준 (해석해/문헌) |
|---|---|---|---|---|---|
| **A** | 자유·구속 BD | 점입자/구형, WCA 배제부피 | `Brownian` + `pair.LJ`(WCA 모드) + `external.wall` | φ, `kσ²/kT`, `σ/ℓ_g` | `D*=1` 자유확산 · 등분배 `⟨x²⟩=kT/k` · 침강 지수분포 `ℓ_g` · 슬릿 plateau MSD |
| **B** | 미세유변학 | 프로브 + 배경 매질/네트워크; 수동(열요동) & 능동(집게 견인) | `Brownian` + `force.Custom`(트랩) + 배경 pair | `Pe = Fσ/kT`, `kσ²/kT`, `ωτ_B` | 뉴턴유체 극한: `G'=0, G''=ηω` · GSER 왕복 일관성 · 트랩 코너주파수 `f_c = k/2πγ` |
| **C** | 활성물질 | ABP: 자기추진 + 회전확산 | `Brownian`(`moment_inertia=0`) + `force.Active` + `update.ActiveRotationalDiffusion` | `Pe = v₀/(σD_r)` 또는 `v₀τ_r/σ`, φ | 단일입자 ABP MSD 해석식 · `D_eff = D_t + v₀²/[d(d−1)D_r]` · MIPS 상경계 (2D: Pe≳40–60, φ≳0.4) |
| **D** | 응집 콜로이드 | Yukawa(DLVO) / Morse 인력 + WCA | `pair.{Yukawa,Morse,DLVO,ExpandedLJ}` | `T*=kT/ε`, `κσ`, `βU_min`, φ | 2체 결합확률 vs `βU_min` · 낮은 φ에서 2차 비리얼 계수 `B₂` · Smoluchowski 응집속도 초기기울기 · RDF 접촉피크 |
| **E** ★ | **2D 소프트 반발** | 점입자 `U/kT = A/r³` (경질 코어 **없음**) | `Brownian` + `pair.Table` (멱함수 표) | **`Γ = π^{3/2}A`** · `n* ≡ 1` · `βU(r_cut)` · **`η₆ = 4p`** | 완벽격자 `ψ₆ = 1`·결함 `0` · 액체 지수 `p = 1/2` · Zahn 상도 `[출처, 미재현]` |

각 도메인은 `examples/`에 **검증된 참조 케이스 1개**를 두고, 그 수치를 `tests/`에 회귀 테스트로 고정한다.

### 4.1 ★ 첫 물리 캠페인 — 2D 소프트 반발계 (2026-07-29~30)

손그림 `soft-r3-2d-A-sweep` 한 장에서 시작해 **6개 런 · 3,856 프로덕션 런**.
카드 [`soft-repulsive-2d--equilibrium-structure`](../../knowledge/wiki/systems/soft-repulsive-2d--equilibrium-structure.md)
가 `draft` 에서 벤치마크 31개(S1–S31)를 가진 이 저장소 최대 지식체가 됐다.

| 런 | 무엇을 물었나 | 런 수 | 결과 |
|---|---|---|---|
| `soft-r3-2d-A-sweep` | `A` 스윕, 시드 4 | 40 | 상자 모양이 초기조건 비교를 교란함을 발견 |
| `soft-r3-time-resolved` | **언제** 구조가 생기나 | 60 | `τ_relax = 0.03–0.10 τ_d` · 물리 척도(`σ=5 µm`) 최초 부착 |
| `soft-r3-relax-seeds` (2단계) | `τ(A=1)` vs `τ(A=0.1)` | 3,584 | **INCONCLUSIVE** (사전등록 `3σ` 미달) |
| `soft-r3-nconv` | `N` 수렴 + `ψ₆` 유한크기 | 12 | **`A=10` 은 hexatic 이 아니다** (16/16 PASS) |
| `soft-r3-fss` | 멱함수 **형태** 검증 | 48 | `A=0.1`·`A=10` 형태 확인 (`χ²/dof` `0.55`·`0.58`) |
| `soft-r3-hexwin` | Zahn hexatic 창 | 72 | 창 안 3점 전부 **등방 액체** · 절단오차가 지수를 `2.9σ` 편향 |

**물리 결론 (현재 최선값)**

| `A` | `Γ` | `η₆` | 상 | 근거 |
|---|---|---|---|---|
| 0.1 | 0.56 | `2.03 ± 0.08` | 등방 액체 | `p = 0.508 ± 0.020` (액체값 `0.5` 에서 `0.4σ`) |
| 1 | 5.57 | `1.86 ± 0.05` | 등방 액체 | 형태 검증은 SE-제한으로 미결 |
| 10 | 55.68 | **`2.06 ± 0.23`** | 등방 액체 | `r_cut=7.80` 최선값. hexatic 기각 `7.9σ` |
| 10.2–10.8 | 56.8–60.1 | `1.5–2.2` | 등방 액체 | **Zahn 창 안인데 hexatic 이 없다** |
| 31.6 | 176.0 | `0.27 ± 0.25` | 결정 | 결함 `0.020` |

⇒ **액체→결정 브래킷 `A = 10.8–31.6`. hexatic 은 관측되지 않았다.**
⚠ 단, **무작위 시작으로는 경계를 찾을 수 없다** (과냉각). 위에서만 묶은 값이다.

### ★ 이 캠페인이 만든 **방법론** finding 6개 — 다른 카드에도 적용된다

| finding | 규칙 |
|---|---|
| [[coarse-sampling-hides-the-whole-transient]] | `stride ≲ τ_relax/5`. 완화시간을 **먼저** 추정한다 |
| [[fraction-threshold-flips-meaning-between-per-frame-and-aggregate]] | 분율 문턱이 `1/N` 보다 작으면 문턱이 없는 것이다 |
| [[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]] | 시드 4개 SE 로 `3σ` 를 세우려면 `t(3)=5.84` |
| [[low-seed-pilots-give-optimistic-design-power]] | 표본 수를 **미리** 계산한다. 편향 지표는 **폭/잡음 ≳ 20** |
| [[order-parameter-magnitude-cannot-identify-a-phase]] | 크기도 기울기도 각각으로는 부족하다 — 둘 다 본다 |
| [[provenance-must-have-one-definition-and-three-capture-points]] | 봉인·궤적·분석 세 시점에서 잡는다 |

**이 여섯이 `master_plan` 의 §11(감도)·§S2(예측)를 실질적으로 개정한다** — §8.1 참조.

---

## 5. 무차원화 규약 (§S4 상세)

### 5.1 기준 척도
| 기호 | 선택 | 이유 |
|---|---|---|
| 길이 `L*` | `σ` = 대표 입자 **직경** | 상호작용 거리와 배제부피의 자연 척도. (반지름 아님 — 혼동이 가장 흔한 버그) |
| 에너지 `E*` | `k_BT` | BD는 열요동이 주역. `T*=kT/ε` 를 자유 무차원수로 남김 |
| 시간 `T*` | `τ_B = σ²/D₀` | 자기 직경만큼 확산하는 시간. `D₀ = k_BT/γ₀`, `γ₀ = 6πη(σ/2)` |

이 선택 하에 HOOMD 입력은 **`σ*=1, kT*=1, γ*=1 ⟹ D₀*=1, τ_B*=1`** 로 고정된다.
따라서 **HOOMD 시간 1단위 = 1 τ_B**. 물리 시간 환산은 `t_si = t_star · τ_B,si`.

> ⚠ HOOMD 문서의 기본 시간단위 `τ_LJ = σ√(m/ε)`와 혼동 금지.
> BD는 과감쇠라 질량이 동역학에 안 들어가므로, `m*=1`로 두고 시간척도는 `γ`가 정한다.
> 이 규약을 어기면 모든 시간이 조용히 틀린다. `tests/test_nondim.py`가 이를 감시한다.

### 5.2 주요 변환
```
σ_si   = 2 a_si                          대표 직경
γ₀_si  = 6 π η_si a_si                   Stokes 항력 (구, 무한매질)
D₀_si  = k_B T_si / γ₀_si
τ_B_si = σ_si² / D₀_si
F* = F_si σ_si / k_BT_si                 힘
k* = k_si σ_si² / k_BT_si                스프링 상수
v* = v_si σ_si / D₀_si                   속도 (= Pe)
D_r* = D_r,si · τ_B_si                   회전확산
ω* = ω_si · τ_B_si                       각주파수
G* = G_si σ_si³ / k_BT_si                탄성률
```

### 5.3 무차원수 목록 (`simbot.nondim.groups`)
| 기호 | 정의 | 물리적 의미 | 레짐 경계 |
|---|---|---|---|
| φ | `N v_p / V` | 부피(면적)분율 | 희박 <0.05, 농후 >0.3, RCP 0.64 |
| `Pe_F` | `F σ / k_BT` | 구동력 vs 열운동 | ~1에서 전이 |
| `Pe_a` | `v₀ / (σ D_r)` | 활성 지속성 | MIPS ≳ 40 (2D) |
| `T*` | `k_BT / ε` | 열 vs 인력 | 응집 ≲ 0.3 |
| `κσ` | 스크리닝 | 전기이중층 두께 | 장거리 <1, 단거리 >5 |
| `k*` | `k σ²/k_BT` | 구속 강도 | 강구속 ≫1 |
| `σ/ℓ_g` | 중력 Pe | 침강 vs 확산 | ~1에서 전이 |
| `Re` | `ρ v a / η` | 관성 (BD 전제 검증) | 반드시 ≪1 |
| `τ_i/τ_B` | `m/(γ τ_B)` | 과감쇠 타당성 | 반드시 ≪1 |

### 5.4 dt 선택 규칙 (`simbot.nondim.choose_dt`) — 구현됨

BD 업데이트: `Δr = (F/γ) Δt + √(2 D₀ Δt) ξ`

**★ 제약은 SI 로 계산하고 마지막에 카드의 시간 척도로 나눈다.** 무차원 단위로 바로
쓰면 "`Δt*` 의 `*` 가 어느 시간인가"를 매번 따져야 하고, `τ_D` 와 `τ_trap` 이 24만 배
차이 나는 계에서는 그 실수가 조용히 통과한다.

제약들의 **최소값**을 채택하고 어느 제약이 지배했는지 기록한다:

| # | 제약 | 식 (SI) | 기본 목표값 | 켜지는 조건 |
|---|---|---|---|---|
| 1 | 열 변위 | `√(2 D₀ Δt) ≤ δ_th σ` | `δ_th = 0.03` (성분별) | **겹칠 상대 ∪ 결합 상대** 있을 때 |
| 2 | 힘 변위 | `max\|F\| Δt/γ ≤ δ_F σ` | `δ_F = 0.005` | 위와 같음 + `max\|F\|` 실측 |
| 3 | **강성 안정성** | `Δt ≤ s · 2γ/λ_max` | `s = 0.2` | **결합·각 있을 때** |
| 4 | 최단 완화시간 | `Δt ≤ ζ · min(τ_trap, 1/D_r, …)` | `ζ = 0.01` | 구속·활성 있을 때 |
| 5 | 활성 변위 | `v₀ Δt ≤ δ_a σ` | `δ_a = 0.01` | 능동 구동 있을 때 |
| 6 | **정확도 목표** | `Δt ≤ 2b/(1+b) · τ_trap` | 목표 편향 `b` | 조화 트랩 + `b` 명시 |

`λ_max = 4 k_bond + 16 k_angle/b²` (`derive()` 가 파생). `1·2·5` 는 **정확도**, `3` 은
**안정성**, `4·6` 은 **관측 가능성**이다 — 세 종류를 하나로 합칠 수 없다.

**활성 제약이 하나도 없으면 예외를 던진다.** 근거 없이 기본값으로 런이 나가는 것이
가장 나쁘다.

**게이트 식은 `nondim.dt_max_{thermal,force,active,stability}` 가 소유한다** — 단위 무관
함수라 `choose_dt` 는 SI 로, 축약 단위로 도는 스크립트는 `σ=γ=D₀=1` 로 같은 함수를 부른다.
문턱은 전부 `config/run_policy.yaml` §timestep 에서 온다. 스크립트가 숫자를 다시 쓰면
정책을 고쳐도 따라오지 않는다 — 2026-07-28 에 `scripts/chain_bend.py` 가 실제로 그랬다.

#### ★ 변위 게이트는 만능이 아니다 (2026-07-28 실측)

첫 손그림 계에서 네 제약을 같은 축에 올린 결과:

| 제약 | `Δt*` 상한 (`τ_trap` 단위) |
|---|---|
| 열 변위 | **`108.6`** |
| 완화시간 | **`0.01`** ← 지배 |

**비 1086배.** 변위 게이트만 켜면 아무것도 막지 못한다 — 기준 길이가 `σ` 인데
입자가 탐색하는 거리는 `ℓ_trap = σ/491` 이기 때문이다.
[`dt-gate-should-be-displacement-based`](../../knowledge/wiki/findings/dt-gate-should-be-displacement-based.md)
의 실측 3건은 **전부 쌍 상호작용 계**였고, 그 결론의 적용 범위가 그것이다.
전문: [`displacement-gate-is-1000x-loose-for-traps`](../../knowledge/wiki/findings/displacement-gate-is-1000x-loose-for-traps.md)

⇒ 두 게이트는 경쟁하지 않고 보완한다. **어느 것이 켜지는지는 카드가 정한다.**

#### ★ 변위 게이트는 정확도 게이트일 뿐이다 (2026-07-28 실측, 결합 계)

곧은 사슬에 변위 게이트를 그대로 걸었더니 `Δt* = 4.5×10⁻⁴` 가 선택되고 사슬이 터졌다.
게이트가 **무력화**된 것이다 — 완전히 곧은 사슬은 `max|F*| = 0` 인 **정류점**이므로
힘 게이트가 상한을 주지 못한다. `kT = 0` 에서도 터지므로 확률적 현상도 아니다.

| `k_bond*` | `k_angle*` | `Δt_crit` 실측 (이분법) | `2/λ_max*` 하한 | 실측/하한 |
|---|---|---|---|---|
| `10⁶` | `10⁴` | `1.00e-6` (N=5) · `5.87e-7` (N=9) | `4.81e-7` | `2.08` · **`1.22`** |
| `10³` | `10³` | `1.84e-4` (N=5) · `1.48e-4` (N=9) | `1.00e-4` | `1.84` · `1.48` |

10조합 전부에서 **하한이 실측보다 작다** (비율 `1.22–2.80`) → 게이트로 쓸 수 있다.
`s = 0.2` 는 최악 `1.22` 대비 `6–14` 배 여유다.
전문: [`dt-gate-needs-a-stability-term-for-stiff-bonds`](../../knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md)

**`hard_floor` 는 정확도 제약에만 적용한다.** `k_bond* = 10⁶` 의 안정성 상한 `9.6×10⁻⁸` 은
floor `10⁻⁷` 아래인데, 실측은 그 계가 `Δt* = 10⁻⁶` 까지 안정임을 보였다 — floor 쪽이 틀렸다.
안정성 상한은 협상 대상이 아니므로(낮추라는 요구가 아니라 `k_bond` 를 낮추라는 요구다)
기각하지 않고 `logged["dt_star_below_hard_floor"]` 에 남긴다.
⇒ 미결: `hard_floor` 는 **카드 기준시간 단위**라 카드마다 뜻이 달라진다. `dt/τ_D` 게이트에서
이미 폐기한 "보편 규약" 함정이고, 비용 검사라면 `Λ` 예산 추정이 할 일이다.

`max|F|`는 초기배치에서 실제 힘을 계산해 얻는다(**추정 금지**). **"재봤더니 0"과 "아직 안
재봤다"를 같은 문장으로 적지 않는다** — 앞은 물리(정류점), 뒤는 절차 위반이고,
`DtChoice.table()` 이 둘을 다르게 표시한다.
`Δt/τ_D` 는 `logged` 에 **기록만** 한다. 게이트로 쓰지 않는다.

---

## 6. `knowledge/` 지식 축적 스키마

> ## ⚠️ 이 절은 폐기됐다 — §10 이 정본이다
>
> 아래 §6.1~6.4 는 2026-07-28 오전에 설계한 스키마다. 같은 날 `BD_agent/knowledge/`
> 를 이관하면서 **그쪽 계약이 더 낫다고 판정하고 정본으로 채택했다** (§10.2).
> 실제 디렉터리는 `source/papers` + `wiki/{systems,findings,concepts,techniques,
> benchmarks,questions}` 이고, 계약은
> [`knowledge/wiki/CLAUDE.md`](../../knowledge/wiki/CLAUDE.md) 가 소유한다.
>
> **지우지 않고 남기는 이유:** §10.4 가 "이관이 오늘 작성한 것을 반박한 사례"를
> 인용하고 있어서, 반박당한 원문이 사라지면 그 기록이 읽히지 않는다.
> 아래를 **구현 지침으로 읽지 말 것.**

### 6.1 구조 (폐기 — §10 참조)
```
knowledge/
├── INDEX.md              전체 목록 + 한 줄 요약 (에이전트가 매번 먼저 읽음)
├── systems/              시스템 아키타입: 무엇을 시뮬레이션했는가
├── parameters/           파라미터 값의 출처와 근거 (η, T, ε, dt, r_cut…)
├── models/               모델링 결정 근거 (왜 WCA인가, 왜 HI 무시인가)
├── failures/             실패 사례: 증상 → 원인 → 처방
├── validation/           검증 벤치마크 수치 (pytest 회귀 후보)
└── env_log.md            패키지 적립 이력
```

### 6.2 항목 포맷 (모든 파일 공통)
```markdown
---
id: dense-diffusion-hardsphere
kind: validation            # systems | parameters | models | failures | validation
tags: [diffusion, hard-sphere, dense, phi]
created: 2026-07-28
updated: 2026-07-28
runs: [2026-07-28T14-30_free-diff_a1b2c3]   # 이 지식을 만든 run
confidence: medium          # high | medium | low
supersedes: []              # 대체한 이전 항목 id
---

## 요약
한 문단. 이 항목이 주장하는 것.

## 근거
데이터/문헌/식. 수치는 오차와 함께.

## 적용 범위 / 한계
언제 이 지식을 믿어도 되고 언제 안 되는가.

## 참고
문헌, run 링크, 관련 항목 [[다른-id]]
```

### 6.3 `failures/` 특화 포맷 (가장 중요)
```markdown
## 증상
관측된 것. 에러 메시지, 그래프 모양, 이상한 숫자.

## 진단 경로
무엇을 의심하고 어떻게 배제했는가. (다음에 같은 증상 만나면 이 순서로 따라간다)

## 근본 원인
분류: numerical | modeling | interpretation | analysis | environment

## 처방
구체적 수정. 코드/설정 diff.

## 재발 방지
추가한 테스트 / 가드 / 문서. 없으면 "없음"이라고 적을 것.
```

### 6.4 지식 축적 규칙
- **S8에서 knowledge 업데이트 없이 파이프라인을 종료할 수 없다.** (게이트)
- 기존 항목과 모순되는 결과가 나오면 새 항목을 만들고 `supersedes`로 연결. **기존 항목을 조용히 덮어쓰지 않는다.**
- `validation/` 항목 중 재현 가능한 수치는 `tests/test_knowledge_regression.py`로 승격 → 코드 변경 시 지식이 깨지는지 자동 감시.
- `INDEX.md`는 항목 추가 시 자동 갱신 (`simbot.io.reindex_knowledge`).

---

## 7. 환경 및 패키지 적립

### 7.1 env: `simulation_bot` (신규 생성)
기존 `hoomd_slit`은 건드리지 않는다. 새 env로 시작해 필요한 것만 단계적으로 추가하고,
**추가할 때마다 `knowledge/env_log.md`에 "왜 필요했는가"를 기록**한다.

| 단계 | 패키지 | 목적 | 시점 |
|---|---|---|---|
| **1. 코어** | `python=3.12`, `hoomd`, `gsd`, `numpy`, `scipy` | BD 실행, 궤적 I/O, 수치 | 지금 |
| **2. 분석·플롯** | `matplotlib`, `freud`, `pandas`, `h5py` | RDF/S(q)/클러스터, 플롯, 로그 | S6/S7 착수 시 |
| **3. 인풋 처리** | `pillow` | 손그림 이미지 메타/전처리 | S1 착수 시 |
| **4. 개발** | `pytest`, `pyyaml` | 테스트, spec 직렬화 | 코어 구현 시 |
| **5. 3D 렌더** | `fresnel` | 3D 스냅샷 레이트레이싱 | 3D 케이스 등장 시 |
| **6. 동영상** | `ffmpeg` | MP4 애니메이션 | GIF로 부족할 때 |

`environment.yml`은 단계가 올라갈 때마다 갱신하고, `conda env export --from-history`로 고정한다.

### 7.2 실행 규약
- 모든 Python 실행은 `simulation_bot` env의 인터프리터 절대경로 사용
  (`conda activate`는 non-interactive shell에서 불안정)
- `CLAUDE.md`에 인터프리터 경로를 기록해 에이전트가 매번 참조

---

### 7.3 로컬 CPU 병렬 컴퓨팅 — 실측 검토

> 측정일 2026-07-28 · 기기 Apple M4 (`Mac16,12`), **4 P-core + 6 E-core = 10코어**, 16 GB
> 벤치마크 커널: 3D WCA + `md.methods.Brownian`, φ=0.30, `dt=1e-4`, Cell nlist(buffer 0.3)
> 원자료: `knowledge/parameters/local_cpu_parallelism.md`

#### 7.3.1 결론 요약 (먼저 읽을 것)

| 발견 | 근거 |
|---|---|
| **HOOMD 7.1은 이 환경에서 완전 단일스레드다.** | 측정 `CPU util = 1.00x` (N=500~32000 전부) |
| **CPU 병렬화 경로가 아예 없다.** HOOMD v3+는 TBB를 제거했고, 남은 경로는 MPI 도메인 분할뿐인데 — | `hoomd.version.mpi_enabled = False` |
| **conda-forge에 MPI 빌드가 osx-arm64·linux-64 **양쪽 모두** 존재하지 않는다.** | 전체 빌드 문자열 검색 결과 `mpi` 매치 0건. 모두 `cpu*`/`gpu*` |
| ⇒ **강한 스케일링(한 시뮬레이션을 여러 코어로)은 불가능.** MPI가 필요하면 소스 빌드 필수. | |
| ⇒ **처리량 병렬화(독립 런 동시 실행)가 유일하고, 이 프로젝트에는 오히려 더 적합하다.** | §7.4 |

#### 7.3.2 단일 프로세스 처리량 — N에 거의 무관

| N | TPS (steps/s) | 입자·스텝/s | CPU util |
|---|---|---|---|
| 500 | 13 313 | 6.66e6 | 1.00x |
| 2 000 | 3 210 | 6.42e6 | 1.00x |
| 4 000 | 1 629 | 6.52e6 | 1.00x |
| 8 000 | 785 | 6.28e6 | 1.00x |
| 32 000 | 189 | 6.05e6 | 1.00x |

**처리량 상수 `Λ ≈ 6.3 × 10⁶ 입자·스텝/s` (P-코어 1개).**
N을 64배 늘려도 9 %만 저하 → Cell 리스트가 잘 동작하고 작업집합이 캐시에 들어간다.
이 상수 하나로 모든 런의 wall time을 예측할 수 있다: `wall ≈ N × steps / Λ`.

#### 7.3.3 동시 실행 스케일링 (N=4000, 12 000스텝, 동기화 시작)

| k | 총 TPS | 스피드업 | 코어당 효율 | 프로세스당 TPS | 추가 1개당 이득 |
|---:|---:|---:|---:|---:|---:|
| 1 | 1 629 | 1.00x | 100 % | 1 629 | — |
| 2 | 3 090 | 1.90x | 95 % | 1 545 | +1 461 |
| 3 | 4 524 | 2.78x | 93 % | 1 508 | +1 434 |
| **4** | **6 027** | **3.70x** | **93 %** | 1 507 | +1 503 |
| 5 | 6 304 | 3.87x | 77 % | 1 261 | +277 ← **절벽** |
| 6 | 6 804 | 4.18x | 70 % | 1 134 | +500 |
| **8** | **8 041** | **4.94x** | 62 % | 1 005 | +309/개 |
| **10** | **8 906** | **5.47x** | 55 % | 891 | +432/개 |
| 12 | 8 663 | 5.32x | 44 % | 722 | **−243 (회귀)** |

**읽는 법**
1. **k ≤ 4는 거의 선형(93 %).** 4개의 P-코어에 하나씩 얹히는 구간.
2. **k=5에서 절벽** — 5번째부터 P-코어를 공유하기 시작한다. 프로세스별 TPS 산포가
   `1.02x`에 불과한 것으로 보아, macOS는 프로세스를 P/E에 **고정하지 않고 시분할 이동**시킨다.
   따라서 "느린 런 하나가 배치를 지연시키는" 문제는 없다 — 배치 스케줄링에 유리한 성질.
3. **E-코어의 실효 성능은 P-코어의 약 1/3이다.**
   4 P × 1 507 = 6 027, k=10 총 8 906 ⇒ 6 E가 2 879 기여 ⇒ E당 480 TPS ⇒ `480/1507 = 32 %`.
4. **k=12는 k=10보다 느리다(−2.7 %).** 오버서브스크립션은 순손실. **상한은 코어 수 10.**

#### 7.3.4 권장 동시 실행 수

| 시나리오 | **k** | 총 TPS | 최대치 대비 | 코어당 효율 | 채택 이유 |
|---|:---:|---:|---:|---:|---|
| 대화형 (사용자가 동시에 작업 중) | **4** | 6 027 | 68 % | 93 % | P-코어만 사용. 머신 반응성 유지, 코어당 효율 최고 |
| **기본값 (권장)** | **8** | 8 041 | **90 %** | 62 % | 최대 처리량의 90 %를 얻으면서 **2코어를 여유로 남긴다** — 에이전트 자신, freud/numpy 분석(Accelerate BLAS는 멀티스레드), 플로팅, OS |
| 배치 (무인 실행) | **10** | 8 906 | 100 % | 55 % | 절대 최대. 다른 작업이 없을 때만 |
| — | 12+ | 8 663 | 97 % | 44 % | **측정상 회귀. 금지** |

**기본값을 10이 아니라 8로 두는 이유**: 10을 써서 얻는 추가 이득은 +11 %인데, 그 대가로
분석·플로팅·에이전트 자신이 실행될 코어가 사라진다. 파이프라인은 시뮬레이션만 돌리는 게
아니라 S6·S7 분석을 계속 병행하므로, 8이 실효 처리량이 더 높다.

---

### 7.4 최적 시뮬레이션 조건 (제안) — 사람이 덮어쓸 수 있음

> 아래는 §7.3 실측에서 **유도된 제안값**이다. 확정값이 아니다.
> 사람이 바꿀 곳: **`config/run_policy.yaml`의 `overrides:` 블록 한 곳.**
> 에이전트는 매 런 시작 시 조건표를 제시하고 승인을 받되, `overrides`가 있으면 그것을 우선한다.

#### 7.4.1 핵심 정책: **오차 막대는 공짜다**

HOOMD가 단일스레드이고 코어가 10개이므로, **독립 시드 4개를 동시에 돌리는 비용은 1개와 같다**
(k≤4는 93 % 효율). 따라서:

> **기본 정책: 모든 프로덕션 런은 시드 ≥ 4개를 동시 실행한다.**
> 단일 런 + 시드 1개는 금지 (오차 막대 없는 결과가 되기 때문).

이것이 이 하드웨어에서 가장 중요한 설계 결론이다. 긴 런 1개보다 짧은 런 4개가 거의 항상 낫다.

#### 7.4.2 실행 티어

`Λ = 6.3e6 입자·스텝/s`, `η_k` = 프로세스당 효율(§7.3.3), `dt* = 5e-5 τ_B` 기준.
`wall ≈ N × steps / (Λ · η_k)` , `t_total = steps × dt*`

| 티어 | 목적 | N | steps | `t_total` | 동시 k | 프로세스당 wall | 배치 총 wall |
|---|---|---:|---:|---:|:---:|---:|---:|
| **T0** `smoke` | 코드가 도는지만 확인 | 200 | 2e3 | 0.1 τ_B | 1 | **0.06 s** | 즉시 |
| **T1** `pilot` | 가드·dt 검증 + wall 예측 (S5-2 **필수**) | 프로덕션과 동일 | 프로덕션의 0.5 % | — | 1 | **≤30 s** | ≤30 s |
| **T2** `explore` | 레짐 탐색 / 파라미터 스윕 | 1 000 | 4e5 | 20 τ_B | **8** | 103 s | **8점 스윕 1.7분** |
| **T3** `production` | 본 측정 (기본 티어) | 4 000 | 1e6 | 50 τ_B | **4** | 11.4 분 | **4시드 11.4분** |
| **T4** `long` | GSER · MIPS · 농후 장시간 | 4 000 | 1e7 | 500 τ_B | 2–4 | 1.9 h | **사용자 승인 필수** |

**티어 선택 규칙 (에이전트가 자동 적용)**
- 처음 보는 시스템 → 반드시 `T0 → T1 → T2` 순서. T2 결과로 레짐을 확인한 뒤 T3.
- `t_total` 요구는 측정하려는 양이 정한다:
  - 확산계수 `D` → `t_total ≳ 10 τ_B` (T2도 충분)
  - 농후계 장시간 `D_L` → `t_total ≳ 100 τ_B` (T3 하한)
  - **GSER `G'(ω), G''(ω)` → MSD가 3~4 decade 필요 → `t_total ≳ 100 τ_B` (T3는 marginal, T4 권장)**
  - **MIPS 조대화 → `t_total ~ 10³ τ_B` + `N ≳ 10⁴` → T4 전용**
- `wall_time_budget_s = 600` (기본). 초과 예상 시 **실행하지 않고 사용자에게 보고**한다.

#### 7.4.3 N 선택 논리 (제안)
- 통계 오차는 `1/√(N × 독립 시간원점 수)`. **희박계에서 입자는 서로 독립**이므로,
  N을 키우는 것과 시간을 늘리는 것이 통계적으로 등가 → **wall time이 `N × steps`에 비례하니 무차별.**
  → 이 경우 **시드 수를 늘리는 것이 유일한 진짜 이득** (계통오차까지 드러남).
- 구조·상거동은 다르다: `L > 2 × 상관길이`가 강제 조건.
  `L = (N v_p/φ)^{1/3}`이므로 φ가 낮으면 같은 L에 N이 적게 든다.
- 활성계는 `L > 2 × v₀τ_r` (persistence length) 필수.

#### 7.4.4 사람 개입 지점 (3곳뿐)

| 언제 | 어디서 | 무엇을 |
|---|---|---|
| 영구 변경 | `config/run_policy.yaml` → `overrides:` | 티어 기본값, k, wall 예산. 이유를 함께 적으면 S8이 `knowledge/`에 반영 |
| 런 1회 한정 | 에이전트가 S4 종료 시 제시하는 조건표에 답 | "T3 대신 T2로", "N=8000으로" |
| 실행 중 | 예산 초과 보고를 받았을 때 | 승인 / 축소 / 중단 |

그 외에는 에이전트가 위 규칙으로 **스스로 결정하고 근거를 기록한다** (질문 예산, §9-5).

---

## 8. 개발 로드맵

| | 마일스톤 | 내용 | 완료 기준 (DoD) | 실제 |
|---|---|---|---|---|
| `[O]` | **M0** 골격 | env, 디렉터리, 문서 3종, knowledge 스키마 | `pytest` 초기 통과, 문서 존재 | ✅ |
| `[O]` | **M1+M2** 수직 슬라이스 | S1→S8 **전체** 관통 | `REPORT.md` 생성 + 해석해 검증 | ✅ **손그림으로 직행.** 자유확산 대신 조화 트랩 — 그림이 트랩이었다. 7 PASS / 2 INCONCLUSIVE / 0 FAIL |
| `[O]` | **M2.5** 결정론 코어 | 관통을 **코드로** 재현 가능하게 | `cli.py run <spec>` 한 줄로 완주 | ✅ **2.3 s, 첫 런과 비트 일치.** 373 테스트 |
| `[O]` | **M2.7** L1 에이전트 층 | `.claude/skills/` + `agents/` (§12.3–12.4) | 새 손그림 → 스킬이 S1→S8 위임 | ✅ 스킬 3 · 참조문서 5 · 서브에이전트 9. **위임은 아직 시험되지 않았다** (두 번째 그림 필요) |
| `[~]` | **M3** 도메인 확장 | A→D 순차, 각각 `examples/` 참조 케이스 | 도메인별 검증 기준 PASS | A(트랩)만. **트랩+WCA 는 사용자 보류** |
| `[~]` | **M4** 검증 강화 | dt 래더, 유한크기, bootstrap, 회귀 승격 | `INCONCLUSIVE` 판정 로직 동작 | ✅ INCONCLUSIVE·검정력·`converge` 동작. 남음: bootstrap, 유한크기 |
| `[X]` | **M5** 인풋 확장 | 실험영상 → 논문 PDF → 음성 → 텍스트 | 각 모달리티로 M1 재현 | 손그림만 |

**지금 진행: M2.7 완료 → 다음은 §8.1.**

### M1 이 설계와 다르게 진행된 이유 (기록)

DoD 는 "자유확산으로 먼저 관통, `D* = 1.00 ± 0.03`" 이었다. 실제로는 **사용자 그림이
광집게였고, 그 계에서는 자유확산이 관측되지 않는다** (`τ_D/τ_trap = 2.4e5` 이므로
자유확산 구간이 `10⁻⁵ τ_trap` 아래에 있다).

자유확산 케이스를 먼저 만들고 손그림으로 넘어가면 M1 이 **버려지는 작업**이 됐을 것이다.
대신 트랩 카드로 직행했고, 그 결과 `(계 × 목적동역학)` 카드 체계가 실측으로
정당화됐다 — 보편 `τ_D` 규약이 이 계에서 `Δt = 12 τ_trap` 을 만든다는 것을 확인했다.

⇒ **자유확산 회귀 케이스(`D* = 1.00 ± 0.03`)는 아직 없다.** §8.1-①.

---

## 8.1 다음 할 일 — **2026-07-30 갱신**

> 판단 기준: **무엇이 다른 것을 막고 있는가**, 그리고 **무엇이 조용히 틀릴 수 있는가.**

### ① 결정에서 녹여 상 경계 묶기 (S30) — **물리 쪽 최우선** (~1시간)

`soft-r3-hexwin` 이 Zahn 창 안 세 점을 전부 등방 액체로 읽었지만, **무작위 시작으로는
경계를 찾을 수 없다** — 1차 전이 근처의 핵생성 장벽 때문에 임의로 오래 과냉각으로
남는다. 정상상태처럼 보이는 것도 구별에 도움이 안 된다.

⇒ **육방 결정에서 출발해 녹는지** 본다. 무작위는 위에서, 결정은 아래에서 묶어
그 사이가 참 경계다. 육방정합 상자를 써야 하고, `ψ₆` 는 상자 모양에 무관하므로
비교가 유효하다 ([[box-shape-confounds-initial-condition-comparison]]).
**이것 없이는 "hexatic 이 없다" 를 주장할 수 없다.**

### ② `A = 13.3` 재측정 (S31) — 브래킷 하한 (~15분)

§8.7 의 `A=13.3` 결정 판정은 `βU(r_cut) = 0.24 kT` 로 쟀다. 절단오차가 지수를
`2.9σ` 움직였으므로(§8.8) **결정화 여부도 영향받을 수 있다.** `r_cut = 7.80` 으로
다시 재면 브래킷이 `10.8–31.6` (193 %) 에서 크게 좁아진다.

### ③ 비용 모델 수정 — **정책이 거짓을 말한다** (30분)

`hexwin` 이 예상 22분 vs **실측 54분**이었다. 오버헤드 계수 `3.4` 를 `N=400` 에서
얻었는데 프레임 추출 비용이 `N` 과 함께 커진다. `estimate_wall_time_s` 가
큰 `N` 배치에서 **2.5배 낙관적**이다 — 예산 게이트가 통과시키면 안 될 런을 통과시킨다.

### ④ §11 감도 분석을 캠페인 결과로 개정 (1시간)

§4.1 의 방법론 finding 6개가 §11 을 실질적으로 개정한다. 특히:
- §11.5 판정 규칙에 **설계 검정력 사전계산**(`seeds_for_target_sigma`)을 넣는다
- 허용오차 유도에 **`t(ν)` 보정**을 강제한다 (현재는 정규 분위수)
- **순차 설계 사전등록** 규약을 §S2 에 추가한다 (`no_stage_N` 을 미리 박는다)

### ⑤ 자유확산 회귀 케이스 (Q9) — **아직도 비어 있다** (반나절)

`brownian` 척도 경로가 여전히 end-to-end 로 실행된 적 없다. 소프트 반발계는
`Soft2DRunConfig` 라는 별도 경로를 썼으므로 **이 부채를 갚지 않았다.**

### ⑥ 보류 — 사용자 판단 대기

| | 왜 보류인가 |
|---|---|
| `τ_relax` 3단계 (`k = 2145`) | 사전등록에 `no_stage_3: true` 를 박았다. 돌리려면 **사람이 결정**해야 한다 |
| `A = 100` 유한크기 사다리 | **이 기계에서 불가능** — `N=1024` 런 하나가 예산의 3.5배 |
| 트랩 + WCA (M3) | 사용자 보류 (2026-07-28). 큐의 `trap-drag-2d-hex300` 이 요구한다 |
| 큐의 남은 손그림 2장 | `abp-rod-2d-run-flip` · `trap-drag-2d-hex300` |

---

## 8.2 ⚠️ 폐기 — 2026-07-28 시점의 큐 판독 (기록으로 남긴다)

> 아래는 손그림 4장이 큐에 들어왔을 때의 triage 다. **`soft-r3-2d-A-sweep` 은 완주했고
> `chain-bend-2d-oscill` 도 닫혔다** (커밋 `fefd5c9` 이전). 남은 것은 2장이다.
> 원문을 지우지 않는 이유는 "무엇이 병목이라고 판단했고 그것이 맞았는가" 가
> 읽혀야 하기 때문이다 — **쌍 러너가 병목이라는 판단은 맞았다.**

### (구) 다음 할 일 — 우선순위와 근거 (2026-07-28 기준)

> 판단 기준: **무엇이 다른 것을 막고 있는가**, 그리고 **무엇이 조용히 틀릴 수 있는가.**

### ★ 큐에 들어온 손그림 4장 (2026-07-28, 사용자 추가)

`inputs/` 에 새 손그림 4장과 논문 2편이 들어왔다. **네 장 중 어느 것도 현재 코드로
돌아가지 않는다.** 각각 무엇을 요구하는지 (판독 triage — 정식 S1 은 착수 시 수행):

| 그림 | 계 | 도메인 | 명시된 값 | 필요한 것 | 막힌 이유 |
|---|---|---|---|---|---|
| **`soft-r3-2d-A-sweep`** | `U/kT = A/r³` 소프트 반발, 2D 정사각 | **A/D** 구조 | `A = 0.1, 1, 10, 100` · `N=100` · `Lx=Ly` | `r⁻³` 쌍 포텐셜 · **sweep** · RDF · Voronoi | 쌍 러너 없음 · `sweep:` 미지원 · `analysis/structure.py` 없음 |
| **`trap-drag-2d-hex300`** | 2D 육방 평형에서 프로브 1개를 끌기 | **B** 능동 미세유변학 | `N≈300` · `R=5 μm` · `k_t=10 pN/μm` · `v=0.5 μm/s` | 쌍 상호작용 · **움직이는 트랩** `r_trap(t)=r₀+vt` · 항력 측정 | 위 + `HarmonicTrap` 중심이 고정 |
| **`chain-bend-2d-oscill`** | 광집게로 잡은 사슬을 `y=a sin(ωt)` 로 굽힘 | **B** 능동 미세유변학 | `k_t=10 pN/μm` · `R=5 μm` | **결합 포텐셜** · 진동 트랩 · GSER → `G'(ω)`, `G''(ω)` | ★ **`U_ij` 가 그림에 공백** ("Eric Furst 논문 참고") · `analysis/microrheo.py` 없음 |
| **`abp-rod-2d-run-flip`** | 단일 활성 타원체, run-and-flip | **C** 활성 | `τ_R=0.5 s` · `v ≤ 5 μm/s` | 이방성 입자 · 방향 자유도 · 180° 뒤집기 · **MSAD** | `abp` 카드 `draft` · `active_run_length` 척도 `NotImplementedError` · 이방성 항력 |

**공통분모가 답을 정한다:**
- **4장 중 3장이 쌍 상호작용을 요구한다** → 쌍 러너가 최우선 병목
- **4장 중 2장이 시간의존 트랩을 요구한다** (끌기·진동) → `HarmonicTrap` 에 `center(t)`
- 논문 2편(`PhysRevLett.94.138301`, `la7023617`)은 `chain-bend` 의 공백 `U_ij` 를
  메우려고 추가된 것으로 보인다. `pypdf` 설치도 같은 이유다 (env-log 3단계)

⇒ 사용자가 보류했던 **트랩+WCA 가 이제 큐의 임계 경로에 있다** (`trap-drag-2d-hex300`).

### ⓪ 쌍 상호작용 러너 + `soft-r3-2d-A-sweep` — **큐의 병목** (1일)

★ **큐를 읽고 나서 최우선으로 올렸다.** 4장 중 3장이 막혀 있는 지점이고,
`soft-r3` 이 그중 **가장 깨끗한 입구**다 — 트랩도 구동도 없는 순수 쌍 상호작용이다.

이 하나가 동시에 해결하는 것:

| | 어떻게 |
|---|---|
| 쌍 러너 (3장 공통 병목) | `run.py` 에 `RUNNERS["...equilibrium-structure"]` 추가 → **디스패치 설계가 처음 시험된다** |
| `dt` **변위 게이트** | 쌍 상호작용이 생기면 게이트가 켜진다 → **처음으로 `dt` 를 구속한다** (①의 관심사 절반이 여기서 해결) |
| `brownian` 척도 경로 (`σ`, `τ_D`) | 이 계의 카드가 그것을 쓴다 → **end-to-end 로 처음 실행** (①의 나머지 절반) |
| `sweep:` 지원 (§9-7) | `A = 0.1, 1, 10, 100` 이 그림에 **명시**돼 있다. 넷 다 돌려야 답이 나온다 |
| `analysis/structure.py` | RDF · Voronoi · ψ₆. **검증 대상이 처음 생긴다** (`freud` 설치돼 있음) |
| `max\|F\|` 실측 경로 | 힘 변위 제약이 `n/a` 에서 벗어난다 (§5.4 "추정 금지") |

**검증 가능하다:** `A=100` 이면 육방 결정화, `A=0.1` 이면 거의 이상기체.
2D 융해는 `knowledge/source/papers/1999-zahn-two-stage-melting-2d.md` 에 증류가 있다.
`N=100` 이라 비용도 낮다.

⚠ `r⁻³` 는 HOOMD 내장 포텐셜이 아니다 → `md.pair.Mie(n=3, m=0)` 또는 `pair.Table`
중 어느 것이 맞는지 확인이 먼저다.

### ① 자유확산 회귀 케이스 — **검증 부채** (반나절)

> ⓪ 이 이 항목의 상당 부분을 흡수한다 (`brownian` 척도 · 변위 게이트 · 러너 분기).
> 그래도 `D* = 1.00 ± 0.03` **정확해 대조 자체**는 따로 남는다 — 쌍 상호작용이 있으면
> `D` 가 1 이 아니게 되므로 순수 자유확산 케이스가 별도로 필요하다. ⓪ 의 smoke 로 붙인다.

`D* = 1.00 ± 0.03`. M1 의 원래 DoD 이고 아직 없다. 이것이 없어서 **세 가지가 미검증
상태로 남아 있다:**

| 미검증 | 왜 위험한가 |
|---|---|
| `CARD_SCALE_RULES` 의 `brownian` 경로 (`σ`, `τ_D`) | **end-to-end 로 한 번도 실행된 적 없다.** 카드 체계의 중심 주장("계마다 척도가 다르다")이 절반만 시험됐다 |
| `dt` **변위 게이트** | 트랩 계에서는 꺼진다 → **실제로 `dt` 를 구속한 런이 0건.** 단위 테스트만 있다 |
| `RUNNERS` 분기 | 항목이 1개뿐이라 러너 디스패치 설계가 시험되지 않았다 |

필요한 것: `run.py` 에 자유 BD 러너 (트랩 포스 없음) + `passive-sphere--transport`
카드 등록 + `analysis/msd.py` (여기서 처음 **검증 대상**이 생긴다).

**해석해가 정확하고(`D* = 1`) 비용이 낮다.** 틀리면 즉시 드러난다.

### ② 손그림 작성 가이드 — **가장 싼 레버리지** (30분)

`docs/drawing_guide.md`. 첫 그림에서 모호성 2건이 나왔고 **둘 다 이 가이드로 없앨 수
있었다:**

| 모호성 | 가이드가 요구할 것 |
|---|---|
| A1 (2D vs 3D 단면) | "축과 차원을 적어주세요" |
| A2 (`a` vs `R` — 반지름인가 지름인가) | "숫자에 단위와 **무엇의 크기인지**를" |

A2 는 **실험 `f_c` 가 없어서 측정으로 닫지 못했다** — 간결성으로 닫고 반증 조건만
기록했다. 그림에 한 글자만 더 있었으면 애초에 생기지 않았을 모호성이다.

③ 의 성공률을 직접 올린다.

### ③ 큐의 나머지 3장 — 순서와 이유

⓪ 이 끝나면 쌍 러너가 있으므로 다음 순서가 자연스럽다:

| 순 | 그림 | ⓪ 에서 재사용 | 새로 필요한 것 |
|---|---|---|---|
| 1 | **`trap-drag-2d-hex300`** | 쌍 러너 · RDF | 시간의존 트랩 중심 · 항력 측정. **사용자가 보류했던 트랩+WCA 가 여기서 목적을 갖는다** |
| 2 | **`chain-bend-2d-oscill`** | 트랩 (⇧ 에서) | 결합 포텐셜 · GSER. **논문 2편을 먼저 읽어 `U_ij` 공백을 메워야 한다** (`bd-lit-distill`) |
| 3 | **`abp-rod-2d-run-flip`** | — | 이방성 입자 · 방향 자유도 · MSAD. **가장 큰 새 물리** |

각 그림이 시험하는 것 (물리와 별도로):
- **S1 판독 프로토콜** — 첫 그림에서 귀납한 규칙이 새 그림에도 통하는가
- **스킬 위임** — `.claude/` 가 실제로 작동하는가 (Q10, 지금은 **시험되지 않았다**)
- **카드 없는 쌍의 처리** — `nondim` 이 예외를 던지고 카드를 먼저 만들게 하는가
  (4장 전부 카드가 없다 → 이 경로가 4번 발동한다)

### ④ 파일럿 런 — **정책과 코드의 불일치** (1시간)

`run_policy.yaml` 이 `pilot: {mandatory: true}` 로 선언하는데 `cli.py` 가 실행하지
않는다. **지금 정책 파일이 거짓을 말한다.** 트랩 계는 wall 이 2 s 라 필요가 없었고
그래서 미루어졌지만, 런이 분 단위가 되는 순간(트랩+WCA, `N=8000`) 의미를 갖는다.

최소 조치: 구현하거나, 정책에서 `mandatory: false` + 이유를 적는다.

### 보류 중 — 사용자 판단 대기

| | 왜 보류인가 |
|---|---|
| 트랩 + WCA (M3) | **사용자 보류** (2026-07-28). `kT_conf` 가 독립 검사가 되고 `r_cut` 이 의미를 갖는 지점 |
| `analysis/` 5개 | 호출자도 검증 런도 없다. ① 이 `msd.py` 를 열어준다 |
| 실험 데이터 대조 | A2 를 측정으로 닫으려면 필요. 측정값이 없다 |

---

## 9. 추가 제안 (기타 사항)

사용자 요청에 따라 제안하는 항목. 채택 여부는 별도 확인.

### 강력 권장
1. `[O]` **git 초기화** — ✅ 완료. 커밋 28개, `4ac2a53` 이후. BD_agent 와 독립 repo.
2. `[O]` **예측 봉인 (prediction sealing)** — ✅ `simbot/io.py`.
   `SEALED.sha256` 은 표준 `sha256sum` 형식이라 **우리 코드 없이 `shasum -c` 로 검증된다** —
   봉인의 신뢰성이 우리 코드에 의존하면 안 되기 때문이다.
   봉인이 깨지면 `validate_run` 이 **대조표를 만들지 않는다.**
   ★ 설계에 없던 것: 봉인 **후에** 만든 문서를 `unsealed` 로 따로 보고한다.
   "봉인 파일이 통과했다"가 "예측이 봉인됐다"를 뜻하지 않기 때문이다.
3. `[X]` **파일럿 런 (pilot run)** — 미구현. `run_policy.yaml` 에 `mandatory: true` 로
   선언돼 있지만 `cli.py` 가 아직 실행하지 않는다. **트랩 계에서는 wall 이 2 s 라
   필요가 없었고, 그래서 미루어졌다.** 트랩+WCA 에서 처음 필요해진다.
4. `[O]` **단위 접미사 강제** (`_si` / `_star`) — ✅ `Quantity.si` 가 문자열·bool 을
   거부하고, 왕복 테스트가 척도 혼동을 잡는다 (`test_s4_nondim.py`).
5. `[~]` **질문 예산** — 규약은 `CLAUDE.md` 에 있고 **사람이 지킨다.** 코드에 강제 장치는
   없다 (에이전트 층이 없으므로 강제할 지점도 아직 없다). M2.7 에서 스킬이 소유한다.

### 권장
6. `[O]` **`bd-diagnose` 스킬** — ✅ `.claude/skills/bd-diagnose/SKILL.md`.
   배제 순서(통계량 요동 → 자기일관성 → 표본 독립성 → 단위 → 수치 → **그 다음 물리**)를
   담았다. 첫 완주 4건 중 물리 문제가 0건이었던 실측이 근거다.
7. `[X]` **파라미터 스윕 지원** — `spec` 에 `sweep: [...]` 미지원.
   `cli.py converge` 가 dt·N·시드를 흔드는 것은 **수렴 확인**이고 레짐 지도가 아니다.
8. `[~]` **캐시** — 재료는 있다 (`spec.hash()`, `RunDir.completed_stages()`,
   `cli.py resume` 가 완주한 런을 재사용). **`spec_hash` 로 기존 run 을 찾아오는
   조회는 없다** — `resume` 은 디렉터리를 직접 받아야 한다.
9. `[X]` **`REPORT.md` → HTML** — 미구현. **그림 5장이 생겼으므로 이제 이득이 있다**
   (현재는 `figs/*.png` 상대 링크라 리포트만 옮기면 그림이 깨진다).
10. `[X]` **손그림 작성 가이드** (`docs/drawing_guide.md`) — 미작성.
    **첫 그림 판독에서 모호성 2개(A1 차원, A2 `a` vs `R`)가 나왔고 둘 다 이 가이드가
    있었으면 없었을 것이다.** 비용 대비 이득이 가장 큰 미결 항목.

### 검토 후 결정
11. **Langevin 폴백** — `τ_i/τ_B ≪ 1`이 깨지는 계(작은 입자, 저점도)에서는 BD가 부적절.
    자동 감지 후 `md.methods.Langevin`으로 전환할지, 경고만 하고 중단할지.
12. **HI 근사 도입** — v1 비범위지만, 농후계 결과가 문헌과 체계적으로 안 맞으면 필요해진다.
    최소 도입안: Rotne–Prager 이동도 행렬 + Cholesky (N ≲ 10³ 한정, HOOMD 외부 구현).
13. **실험 데이터 직접 대조** — 사용자가 실험 MSD/영상을 주면 시뮬레이션과 같은 파이프라인으로
    분석해 나란히 비교. S7의 확장.

---

## 10. 지식 계층 이관 — `BD_agent/knowledge/` 흡수

`/Users/kyuhwan/Desktop/BD_agent/knowledge/`에 이미 상당히 성숙한 지식 계약이 있다.
**§6에서 내가 새로 설계한 스키마보다 이쪽이 낫다. 이쪽을 정본으로 채택하고 §6을 폐기한다.**

### 10.1 거기에 있는 것
```
BD_agent/knowledge/
├── raw/lab/          원본 PDF + LaTeX 소스 tarball (~20편, gitignored)
├── source/
│   ├── papers/       ★ 논문별 증류 .md 42편 + INDEX.md  ← 즉시 자산
│   └── lab/          미발표 랩 자산
└── wiki/
    ├── CLAUDE.md     ★ 지식 계약 (frontmatter가 기계 계약, 산문이 근거)
    ├── systems/      ★ (계 × 목적동역학) 카드 + _TEMPLATE.md + _index.md
    ├── concepts/  techniques/  benchmarks/  findings/  questions/
```
42편 증류본은 Choi·Gubbala·Arnold·Kim·Xu·Cheon·Takatori·Barakat·Quah·Modica 등
**사용자 랩의 실제 논문 코퍼스**다. 새로 만들 수 없는 자산이다.

### 10.2 채택하는 설계 (내 §6보다 우수한 점)

| 그들의 설계 | 왜 우수한가 |
|---|---|
| **3계층 `raw / source / wiki`** | 저작권(원본)·정본(증류)·해석(합성)이 분리된다. 내 flat 구조는 이걸 못 한다 |
| **`(계, 목적동역학)` 쌍 카드** ★ | **§5를 반박한다** — 아래 10.3 |
| **`reproduced: yes/no/partial`** | 내 `verified`보다 강하다. "논문에 실렸다 ≠ 우리 코드에서 돈다" |
| **`[출처]` vs `[출처, 미재현]` 표기** | 검증 근거와 사실 기록을 리포트에서 구분 강제 |
| **발표여부로 폴더 분할** (`papers/` vs `lab/`) | 공개 경계가 폴더 단위이므로. 랩 논문을 미발표물과 섞으면 공개 시 함께 빠진다 |
| **`precedence L0–L3`** | "낮은 L이 무엇이 사실인가를 이기고, 높은 L이 그래서 무엇을 할까를 이긴다" |
| **`questions/`** 를 삭제하지 않고 `status`로 닫음 | 미해결 문제의 목록이 남는다 |
| `dead-end-<slug>.md` | 막힌 길도 자산 |

### 10.3 ★ 이 이관이 **§5를 수정한다** — 무차원화는 계 하나로 정해지지 않는다

그들의 `wiki/CLAUDE.md`가 실측으로 보여주는 것:

| 계 × 목적동역학 | 기준 길이 | 기준 시간 | `kT` |
|---|---|---|---|
| ABP × 제어 | **런 길이 `ℓ`** | **`τ_r = 1/D_r`** | **유도량** |
| 브러시 콜로이드 × 비평형 접촉 | `σ` | `τ_D = σ²/D` | 입력값 |
| 수동 tracer × 수송 | `σ` | `τ_D` | 입력값 |

> **§5는 `(σ, k_BT, τ_B)`를 보편 규약으로 못박았다. 이는 도메인 C(활성물질)에 틀리다.**
> ABP에서 자연 단위는 런 길이와 회전완화시간이고, `kT`는 유도량이 된다.
>
> **수정:** §5는 **도메인 A·B·D의 기본 규약**으로 격하한다.
> 무차원화 규약의 소유자는 `wiki/systems/<계>--<동역학>.md` 카드다.
> **카드가 없는 쌍을 만나면 즉흥 무차원화 금지** — `_TEMPLATE.md`로 `status: draft` 카드를 먼저 만든다.

게이트도 쌍마다 다르다 (그들의 표):

| 게이트 | 수동 구형 × 평형구조 | ABP × 조밀집단 |
|---|---|---|
| 평형화 판정 | ✅ 유효 | ❌ **무의미** — 능동계는 열평형에 안 간다 |
| `D_msd = kT/γ` | ✅ 성립 | ⚠️ 성립 안 함 — `D_eff = D_t + v₀²τ_r/2` |
| 이류 변위 `v₀Δt/σ` | 해당 없음 | ✅ 필수 |

### 10.4 이관이 **오늘 작성한 것을 반박한 사례** (이미 수정 완료)

`wiki/findings/dt-gate-should-be-displacement-based.md`:

| | 내가 오늘 쓴 값 | 그들의 실측 근거 | 조치 |
|---|---|---|---|
| `dt` 상한 | `hard_ceiling: 1e-4` | 실제 돌아간 런 3건 중 **2건을 기각**한다 (선행 slit 1e-3, Quah 코드 1.67e-4) | ✅ `4.5e-4`로 수정 (`config/run_policy.yaml`) |
| 변위 규약 | `√(2d·Δt)` (d차원 총변위) | 랩 관행은 `√(2Δt)` (성분별) | ✅ `per_component`로 명시 + 환산 `×√d` 기록 |
| 변위 문턱 | `0.02σ` | 실측 0.006 / 0.018 / 0.045σ → 권장 `0.03σ` | ✅ `0.03σ`로 수정, 기본 dt는 0.010σ 유지 |

**교훈:** 게이트를 문헌·실측 없이 정하면 실제로 작동하는 설정을 기각한다.
이관해야 할 이유가 이것 하나로 충분하다.

### 10.5 추가로 얻는 실측 사전 (`findings/lab-bd-conventions.md`)

| 항목 | 랩 실측값 | 우리 정책에 대한 함의 |
|---|---|---|
| 엔진 | BD 논문 **전부 HOOMD-blue** (Quah: 3.8.1) | 우리는 7.1.0 → **API 이식 필요** |
| 실행 하드웨어 | **GPU 가속** 다수. Xu 2024는 `8×10⁸` 스텝 | M4 CPU로 재현 불가 (`N=1000`이면 ~35시간). §7.4 예산 게이트가 반드시 걸러야 함 |
| 배제부피 | WCA 표준. Quah는 **`ε/kT = 500`** 으로도 안정 | 강한 WCA 자체는 위험하지 않다. 위험은 **변위와의 결합** |
| 결정화 억제 | Takatori: 지름비 **1.4**, 몰분율 **2/3:1/3** (`φ`≤0.83) | 고 `φ`에서 **이분산 강제 게이트** 필요. 안 하면 결정을 "유리"로 오판 |
| 통계 | Xu 2023: **20 realization 평균** | 우리 T3 기본 4시드는 랩 관행보다 약하다 → §10.6 미해결 |
| `φ` 정의 | `nπσ²`(지름) vs `n̄πa²/4` — **논문마다 다름** | **조용히 4배 틀릴 수 있다.** `simbot`은 φ 계산 시 지름/반지름 규약을 항상 기록 |

### 10.6 이관 계획 및 미해결

| | 단계 | 비고 |
|---|---|---|
| `[O]` | `wiki/CLAUDE.md` 계약 이식 | ✅ `knowledge/wiki/CLAUDE.md` 로 그대로 |
| `[O]` | `source/papers/` 42편 + `INDEX.md` | ✅ 복사 완료. 즉시 자산 |
| `[X]` | `wiki/systems/` 카드 + `_TEMPLATE.md` + `_index.md` 복사 | §5 수정의 근거 |
| `[X]` | `wiki/findings/`, `benchmarks/` 복사 | `benchmarks.yaml`은 pytest 회귀의 근거 |
| `[X]` | 오늘 만든 항목 3개를 새 스키마로 재작성 | `water_298k` → `wiki/concepts/`, `no_hydrodynamics` → `wiki/concepts/`, `local_cpu_parallelism` → `wiki/techniques/` |
| `[X]` | `raw/` 는 복사 여부 결정 필요 | 2.3 MB. gitignore 대상이므로 원본 위치 참조만으로도 충분할 수 있음 |
| `[X]` | **미해결:** 기본 시드 수 4 vs 랩 관행 20 | k=8이면 8시드가 1.7분(T2). T3에서 몇 개로? → 사용자 판단 필요 |
| `[X]` | **미해결:** HOOMD 3.8.1 → 7.1.0 API 이식 표 작성 | `gamma.default` 문법 등 |

> ⚠️ 이관은 **복사가 아니라 채택**이다. 그들의 `master_plan.md`·`docs/00_decision_log.md`도
> 읽고 결정 이력(D1–D16 등)을 승계해야 한다. 아직 읽지 않았다.

---

## 11. 감도 분석 (Sensitivity Analysis) — S7b

> 파이프라인의 어느 단계에도 없었던 것. **S7 검증 직후, S8 결론 직전에 삽입한다.**

### 11.1 왜 필요한가 — provenance와 직결된다
S1·S3은 그림에 없는 값을 `provenance: assumed`로 채운다. 결론이 그 가정에 의존하면
**결론이 아니라 추측이다.** 감도 분석은 이 질문에 답한다:

> "내가 임의로 채운 값이 틀렸다면 결론이 바뀌는가?"

**자동 연결 규칙:** `provenance: assumed` 인 모든 필드가 감도 분석 후보다.
사람이 후보를 고르지 않는다 — 스펙이 스스로 후보를 지정한다.

### 11.2 ★ 감도는 **무차원수 공간에서** 계산한다
원시 SI 파라미터로 감도를 재면 낭비다. `η, T, a`는 오직 `D₀`와 `τ_B`를 통해서만 들어가므로
셋을 따로 흔드는 것은 **같은 방향을 세 번 흔드는 것**이다.

> **규칙:** 감도는 §5(또는 해당 systems 카드)의 무차원수 원장에 대해 계산한다.
> `m`개 SI 파라미터가 `n`개 무차원수로 줄면 런 수가 `2m → 2n`으로 줄고, 보통 `n ≪ m`이다.

### 11.3 4단계 (싼 것부터)

| 단계 | 방법 | 비용 | 언제 |
|---|---|---|---|
| **A. 레짐 근접도** | 각 무차원수가 레짐 경계에서 얼마나 떨어졌나. `d = \|log(X/X_c)\|` | **런 0회** | **S4에서** 항상 |
| **B. 국소 1차 (OAT)** | 무차원수 하나씩 `×2, ÷2` → 무차원 감도지수 `S_i = ∂lnQ/∂lnX_i` | `2n` 런 (T2) | 기본. 항상 |
| **C. 2차 상호작용** | 강한 `S_i` 상위 2~3개만 격자로 조합 | `~9` 런 (T2) | B에서 비선형 징후 시 |
| **D. 글로벌 (Sobol/LHS)** | 가정 상자 전체 표본 | `≥64` 런 | C가 강한 상호작용을 보일 때만. 사용자 승인 |

**A는 런이 필요 없는데 가장 중요할 수 있다.** `Pe = 45`이고 MIPS 경계가 `Pe_c ≈ 40–60`이면,
그 계는 정의상 감도가 극대다 — 돌려보기 전에 알 수 있다.

### 11.4 우리 하드웨어에서의 비용 — 사실상 공짜
§7.3 실측 기준, T2(`N=1000`, 4e5스텝) 1런 = 103 s, 동시 8개 가능:

| 무차원수 개수 `n` | OAT 런 수 | 배치 수 (k=8) | **총 wall** |
|---|---|---|---|
| 2 | 4 | 1 | **1.7분** |
| 4 | 8 | 1 | **1.7분** |
| 8 | 16 | 2 | **3.4분** |

> **결론: 감도 분석을 생략할 이유가 없다.** 무차원수 4개까지는 프로덕션 런 1회의
> 1/7 비용으로 끝난다. **기본 활성화한다.**

### 11.5 판정 규칙

| `\|S_i\|` | 해석 | 조치 |
|---|---|---|
| `> 1` | 결론이 이 가정에 **강하게 의존** | 리포트 상단에 경고. 가정을 좁혀야 함 (문헌 조회 / 사용자 질문 / 조건부 결론) |
| `0.2 – 1` | 보통 의존 | 리포트에 명시 |
| `< 0.2` | **무관** | "이 가정은 결론을 바꾸지 않는다"를 **명시적으로 보고** — 이것도 결과다 |

`S_i`의 부호와 크기가 **예측(S2)과 일치하는지도 검사한다.** 불일치는 모델 이해의 결함 신호.

### 11.6 산출물
- `07b_sensitivity.md` — 감도표 + 레짐 근접도 + 판정
- `figs/tornado.png` — 토네이도 플롯 (`|S_i|` 내림차순 수평 막대)
- `08_conclusion.md`의 "신뢰도와 한계"가 이 결과를 **인용해야 한다** (게이트)

### 11.7 실패모드
- 감도를 SI 파라미터로 계산해 중복 방향에 런을 낭비
- `×2, ÷2` 섭동이 레짐 경계를 넘어가서 `S_i`가 의미를 잃음 (→ 섭동 폭을 `±20%`로 축소 후 재시도)
- T2의 통계오차가 `S_i`보다 커서 전부 `INCONCLUSIVE` (→ 시드 늘리거나 T3로 승격)

---

## 12. 모델 티어링 — 어디에 비싼 리즈닝을 쓸까

> 원칙: **추출은 저가, 해석은 고가.** 그리고 **계산은 LLM이 아니라 코드.**

### 12.1 배분

| 단계 / 작업 | 모델 | 근거 |
|---|---|---|
| **S1 손그림 해석** — 기하·경계·차원 판정, 화살표 의미, 모호성 후보 생성 | **Opus 5** | 멀티모달 + 물리 리즈닝. **여기서 틀리면 뒤의 전부가 틀린다.** 가장 비싼 오류 지점 |
| S1 추출 — 텍스트·숫자·라벨 판독, EXIF·해상도·파일 인덱싱 | **Haiku 4.5** | 정형 추출. 리즈닝 불필요 |
| **S2 예측 리즈닝** | **Opus 5** | 봉인되는 과학적 주장 |
| S3 스펙 채우기 (knowledge 조회 + 규칙 적용) | **Sonnet 5** | 판단 여지 적음 |
| S3 YAML 직렬화·provenance 정리 | **Haiku 4.5** | 정형 |
| S4 무차원화 | **코드**(`simbot.nondim`) + Sonnet 5 검토 | 숫자는 LLM이 만들지 않는다 |
| S5 실행 | **코드만** | LLM 없음 |
| S6 그림 생성 / 캡션 | 코드 / **Sonnet 5** | |
| **S7 판정 + 원인 가설** | **Opus 5** | 인과 추론. FAIL 원인 오판이 가장 비싸다 |
| S7b 감도 해석 | **Sonnet 5** | 수치는 코드, 해석만 |
| **S8 결론** | **Opus 5** | 최종 과학적 주장 |
| S8 knowledge 항목 초안 | **Sonnet 5** | 템플릿 채우기 |
| **실패 진단** (`bd-diagnose`) | **Opus 5** | 가설 생성·배제 추론 |
| 문헌 증류 (`source/` 항목 작성) | **Sonnet 5** | 대량. 단 **식 변환은 Opus 검토** |
| 문헌 대량 스캔·서지 추출·`INDEX` 갱신 | **Haiku 4.5** | 정형 대량 |

### 12.2 안전장치 — 저가 모델이 물리 판단을 하지 못하게

> **규칙: `provenance`가 `inference` 또는 `assumed`인 필드는 Opus만 쓸 수 있다.**
> `observation`·`derived`는 저가 모델이 채워도 된다.

이 규칙이 좋은 이유: 기존 provenance 스키마에 그대로 얹히고, **기계적으로 검사 가능**하다.
`simbot.spec.validate`가 각 필드에 `written_by` 필드를 요구하고 위반을 잡는다.

### 12.3 구현 — ✅ 완료 (2026-07-28)
`.claude/agents/*.md` 의 `model:` frontmatter 로 정의했고, `bd-pipeline` 스킬이 단계별로 위임한다.
**배분표와 실제 `model:` 이 일치하는지 `tests/test_agent_layer.py` 가 검사한다.**

| 서브에이전트 | model | 담당 |
|---|---|---|
| `bd-intake-extract` | haiku | S1 추출 |
| `bd-intake-interpret` | opus | S1 해석 |
| `bd-predict` | opus | S2 |
| `bd-spec` | sonnet | S3 |
| `bd-validate` | opus | S7 판정 |
| `bd-conclude` | opus | S8 |
| `bd-lit-distill` | sonnet | 문헌 증류 |
| `bd-lit-scan` | haiku | 문헌 스캔·인덱싱 |
| `bd-diagnose` | opus | 실패 진단 |

> 비용 절감이 목적이 아니라 **속도와 품질의 배분**이 목적이다.
> 인풋 스캔을 Haiku로 돌리면 손그림 10장 목록화가 몇 초에 끝나고, 그 절약분을
> S1 해석과 S7 판정에 쓸 수 있다.

### 12.4 구현 결과 (2026-07-28) — 참조문서는 8개가 아니라 5개다

**Q6 결정.** 설계는 단계별 참조문서 8개였다. 그런데 결정론 코어가 완성된 뒤
**S3·S4·S5 는 `cli.py run` 한 줄**이 되었고, 각각에 별도 문서를 두면 "이 함수를
호출한다"만 적힌 얇은 파일 3개가 생긴다.

**내용이 있는 곳에 문서를 둔다:**

| 문서 | 왜 독립인가 | 크기 |
|---|---|---|
| `s1_intake_drawing.md` | **코드로 표현할 수 없는 유일한 단계.** 가장 비싼 오류 지점 | 가장 두껍다 |
| `s2_prediction.md` | 봉인·tolerance·검정력 규율. 허술하면 검증이 무력화된다 | |
| `s3_s5_execute.md` | 셋 다 `cli.py` 호출 + 게이트 읽기 → **합쳐야 흐름이 보인다** | |
| `s6_s7_validate.md` | 그림과 판정은 같은 판단(무엇이 이상한가)에 쓰인다 | |
| `s8_knowledge.md` | 결론 서술 + knowledge 계약 | |

**스킬이 물리를 다시 적지 않는다.** `simbot` 함수와 `knowledge/wiki/` 카드를 인용한다.
스킬 층의 유일한 고유 내용은 **S1 손그림 판독 프로토콜**이다.

구조를 [`tests/test_agent_layer.py`](../../tests/test_agent_layer.py) (64개)가 감시한다 —
frontmatter 유효성, 링크 무결성, **§12.1 배분표와 `model:` 일치**, 저가 모델의
권한 경계 명시, `settings.json` 의 봉인 문서 편집 거부.

★ **`settings.json` 이 봉인 문서의 `Edit` 을 거부한다:**

```json
"deny": ["Edit(./runs/**/02_prediction.md)", "Edit(./runs/**/01_intake.md)",
         "Edit(./runs/**/SEALED.sha256)", "Bash(conda activate:*)"]
```

예측은 `cli.py` 가 파이썬으로 **생성**하고, 그 뒤 에이전트가 텍스트 편집으로
**고치는 것**을 막는다. 봉인 검증이 사후에 잡지만 **애초에 못 하게 하는 것이 낫다.**

상세: [`.claude/README.md`](../../.claude/README.md)

---

## 13. 열린 질문

### 닫힌 것

| # | 질문 | 답 (2026-07-28) |
|---|---|---|
| Q1 | git 초기화할까? | ✅ **했다.** 커밋 28개 (2026-07-30) |
| Q2 | 첫 손그림의 주제는? | ✅ **2D 광집게** (`R=5 μm`, `k=10 pN/μm`, `T=300 K`). 텍스트 자유확산 단계를 건너뛰었다 (§8) |
| Q3 | 리포트 언어? | ✅ **한국어 본문 + 영어 기술용어.** 단 **matplotlib 그림의 글자는 영문** (기본 폰트에 한글 글리프 없음) |
| Q4 | 런 1회 허용 wall time 상한? | ✅ **10분.** `budget.wall_time_per_run_s: 600`. `cli.py` 가 초과 예상 시 **실행하지 않고 중단**한다 (테스트로 고정) |
| Q5 | 2D를 1급으로 지원? | ✅ **예.** 첫 계가 2D 였고 `Lz=0` 경로가 검증됐다 (`plateau = 2d`, `⟨r²⟩` 비 `3/2` 를 0.3 % 로 확인) |

### 새로 열린 것

| # | 질문 | 상태 |
|---|---|---|
| Q6 | **`.claude/` 스킬을 어느 입도로 쪼갤까** | ✅ **닫힘 — 스킬 3 + 참조문서 5** (§12.4) |
| Q7 | 트랩 커널 처리량이 기준선의 1.2–1.4배인 이유 | 열림 → [`questions/trap-kernel-throughput-vs-wca-baseline.md`](../../knowledge/wiki/questions/trap-kernel-throughput-vs-wca-baseline.md) |
| Q8 | **`scripts/trap_batch.py` 를 지울까** — `run.run_trap_batch` 가 같은 일을 한다 | 미결. 첫 완주의 재현 경로라 남겨둠 |
| **Q9** | **자유확산 회귀 케이스를 언제 만들까** | **아직 열림** → §8.1-⑤. 소프트 반발계는 `Soft2DRunConfig` 별도 경로를 썼으므로 이 부채를 **갚지 않았다** |
| Q10 | 스킬 위임이 실제로 작동하는가 | ✅ **닫힘 (2026-07-29).** `bd-pipeline` 스킬로 `soft-r3` 캠페인 6런을 완주했다 |
| Q11 | `pilot: {mandatory: true}` 를 구현할까, 정책을 고칠까 | 미결. **지금 정책 파일이 거짓을 말한다** |

### 2026-07-30 에 새로 열린 것

| # | 질문 | 상태 |
|---|---|---|
| **Q12** | **무작위 시작으로 상 경계를 찾을 수 있는가** | ❌ **아니다 — 닫힘.** 과냉각 때문에 상한만 준다. 결정에서 녹여야 한다 (§8.1-①) |
| **Q13** | 절단오차 허용치를 **값 기준**으로 세워도 되는가 | ❌ **아니다 — 닫힘.** `βU(r_cut)` 이 값을 `3σ` 안에 두어도 **지수를 `2.9σ` 편향**시킨다. 지수 기준 허용치가 따로 필요하다 |
| **Q14** | `estimate_wall_time_s` 의 오버헤드 계수 | 열림 → §8.1-③. 큰 `N` 에서 **2.5배 낙관적**이다 (`hexwin` 22분 예상 vs 54분 실측) |
| **Q15** | `χ²` 형태 검정을 4시드로 할 수 있는가 | ❌ **아니다 — 닫힘.** `χ² ∝ 1/SE²` 라 4시드 SE 의 41 % 불확실성이 4배 증폭된다 |
| **Q16** | `A = 13.3` 이 제대로 된 절단에서도 결정인가 | 열림 → §8.1-② (S31) |
| **Q17** | `master_plan` §11(감도)·§S2(예측)를 캠페인 결과로 개정할까 | 열림 → §8.1-④. finding 6개가 실질적 개정을 요구한다 |

---

*이 문서가 설계의 단일 진실 원천이다. 설계 변경 시 코드보다 먼저 이 문서를 고친다.*
*폐기된 절은 지우지 않고 `⚠️ 폐기` 로 표시한다 — 반박당한 원문이 사라지면 그 반박이 읽히지 않는다.*
