---
type: source
kind: paper
lab_authored: true
title: Interpretation of electrostatic self-potential measurements using interface-trapped microspheres with surface heterogeneity
authors:
  - "Choi KH"
  - "Kang DW"
  - "Yoo S"
  - "Lee S"
  - "Park BJ"
year: 2020
journal: ACS Applied Polymer Materials 2(3), 1304-1311
doi: 10.1021/acsapm.9b01189
source_url: "https://doi.org/10.1021/acsapm.9b01189"
raw_file: knowledge/raw/lab/2020-choi-electrostatic-self-potential-heterogeneity.pdf
si_file: knowledge/raw/lab/si/2020-choi-electrostatic-self-potential-heterogeneity__SI.pdf
si_available: true
access: paywall
engine: "실험 (광집게) + MC 시뮬레이션"
reproduced: no
parameters_extracted: yes
affiliation: Kyung Hee University (Park BJ group) — Takatori 랩 이전
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "method"
  - "interfacial-colloid"
  - "pair-potential"
  - "benchmark-candidate"
---
# Interpretation of electrostatic self-potential measurements using interface-trapped microspheres with surface heterogeneity

## 왜 이 위키에 있는가

★ **제1저자.** 표면 불균질이 정전 자기퍼텐셜 해석에 미치는 영향.
계면 포획 입자의 **쌍퍼텐셜 실측값 + MC 시뮬레이션 파라미터 + 항력 처리** 셋을 모두 준다.

**관련도:** `method` · **소속:** Park BJ 그룹 (경희대) — Takatori 랩 합류 이전 연구

## 원문 접근

- DOI: https://doi.org/10.1021/acsapm.9b01189
- ✅ 본문·SI 확보 (2026-07-27, Stanford 기관 접근권). `raw_file`·`si_file` 참조

---

## 1. 계 — 무엇을 측정했나

n-decane / 초순수 평면 계면에 포획된 폴리스티렌 미구. 4종(NPS·SPS·APS·CPS)이 표면
작용기만 다르다. 시분할 광집게로 3입자를 잡아 **쌍상호작용력을 직접 측정**한다.

## 2. 쌍퍼텐셜 — 쌍극자 정전 반발 `r⁻³`

$$\frac{U_{ij}}{k_BT} = a_{ij}\left(\frac{d}{r_{ij}}\right)^{3}, \qquad a_{ij} = \Omega_i \Omega_j$$

- `d` = 입자 **지름**, `r_ij` = 중심간 거리 → **힘은 `F ∼ r⁻⁴`**
- 계면 건너편 비대칭 전하해리로 생긴 쌍극자에서 옴 (apolar 상의 강한 전기장)
- `Ω_i` = 입자 개별 **자기퍼텐셜**(dimensionless). 3입자 쌍측정 3개로 `Ω₁·Ω₂·Ω₃`를 역산
- 쌍대 가산성(pairwise additivity) 가정은 장거리 정전 상호작용에서 타당하다고 본문이 인용(ref 53)

> **규약 함정:** `d`가 **지름**이다. 반지름으로 읽으면 `U`가 8배 틀린다.

## 3. Table 1 — 실측값 (본문 p.1306)

| | NPS | SPS | APS | CPS |
|---|---|---|---|---|
| `Ω` | 인력 (반발 없음) | 84.71 ± 19.45 | 136.18 ± 40.48 | 571.7 ± 94.16 |
| RSD | — | 22.9% | 29.72% | 16.47% |
| `a_ij` | — | 9422.18 ± 6977.75 | 18099.63 ± 9107.79 | 326985.20 ± 79935.24 |
| `θ_c` (deg) | 97.0 ± 0.9 | **99.4 ± 2.1** | 145.4 ± 0.9 | 142.6 ± 1.2 |
| `ψ` (mV) | −51.2 ± 2.4 | −57.5 ± 2.2 | −69.5 ± 2.4 | −65.6 ± 3.2 |
| `d` (μm) | 2.93 ± 0.03 | 2.96 ± 0.05 | 2.79 ± 0.11 | 3.16 ± 0.07 |

접촉각은 gel trapping method(SI ref 3, Paunov)로 SEM 측정. `θ_c`가 90°를 넘으므로
입자 표면의 과반이 **오일 상**에 있다.

> SPS의 `θ_c = 99.4 ± 2.1°` 가 `2023-lee-colloidal-debye-force` 가 "~99.4°" 로 인용한 값이다
> — 인용 사슬이 수치로 확인된다.

### ⚠️ 원문 내부 불일치 — APS/CPS 접촉각이 뒤바뀌어 있다

| 출처 | APS | CPS |
|---|---|---|
| **본문 Table 1** (p.1306) | **145.4 ± 0.9°** | **142.6 ± 1.2°** |
| **SI Figure S8 캡션** (S-9) | **142.6°** | **145.4°** |

Table 1 쪽 열 매핑은 PDF 단어 좌표로 확인했다 — 헤더 `APS@x=383` · `CPS@x=501` 에 대해
`145.4@358` · `142.6@472` 이고, 같은 표의 `ψ` 행(`−69.5@358` · `−65.6@472`)이 본문 문장
*"ψ_APS ≈ −69.5 mV … ψ_CPS ≈ −65.6 mV"* 와 일치하므로 **열 매핑 자체는 독립 검증됨.**
따라서 추출 오류가 아니라 **원문 두 곳이 서로 어긋난다.**

**해소 — 제1저자 판단 (Choi KH, 2026-07-27):** *"값에 따른 차이가 크지 않으니 아무거나 사용."*
두 값 차이가 2% 미만이고 `η_eff` 의 접촉각 의존도 그 범위에서는 완만하다.

**따라서 정본은 본문 Table 1로 고정한다** — `θ_c(APS) = 145.4 ± 0.9°`, `θ_c(CPS) = 142.6 ± 1.2°`.

> "아무거나"가 **조회할 때마다 다른 값**을 뜻하지는 않는다. 파라미터 사전은 결정론적이어야
> 하므로 한쪽을 정본으로 못박고, 어긋난 사실은 여기 남긴다. 정본을 Table 1로 고른 이유는
> 열 매핑이 `ψ` 행을 통해 본문 문장과 교차검증되기 때문이지 SI가 틀렸다는 판정이 아니다.

인용 시 표기: `[Choi 2020 Table 1]` — 2% 이내 정밀도를 요구하는 주장에는 쓰지 않는다.

## 4. MC 시뮬레이션 — 이 논문이 실제로 돌린 것

| 항목 | 값 |
|---|---|
| 차원 | **2D** (계면) |
| `N` | 247 |
| 초기배치 | **실험 스냅샷에서 취득** (Figure S4A,B) — 무작위 아님 |
| 평균 간격 | `r/d ≈ 8.3` |
| MC 이동 스텝 | `Δ(r/d) = 0.05` (radial) |
| 사이클 | `10⁶` |
| 평형 판정 | 평균 입자간 거리의 **plateau** |
| 샘플링 | 4000+ 사이클, 사이클마다 rdf 계산 후 평균 |
| 상호작용 | `U_tot,i/k_BT = Σ_{j≠i} Ω_i Ω_j (d/r_ij)³` |
| **브라운 확산** | **무시함** — `r/d ≈ 8.3` 에서 `U_ij ≈ 577 k_BT` 로 충분히 강한 반발 |
| 비교 조건 | `Ω_hetero`(분포) vs `Ω_max = 804.33` vs `Ω_min = 202.18` |

**결과:** `Ω_hetero`를 넣은 rdf가 실험과 "excellent agreement", 균일값(`Ω_max`/`Ω_min`)은 덜 맞음.

### SI 그림이 주는 추가 정보 (S-4 ~ S-12)

| 그림 | 내용 | 우리 쪽 쓸모 |
|---|---|---|
| **S4** | 실험 미세구조(CPS 247입자)의 원본·이진화 이미지 + **Voronoi 다이어그램**. 최근접이웃 5/6/7개를 적/회/청으로 칠해 **결함**을 표시 | 초기배치 출처. Voronoi·이웃수는 `freud.locality.Voronoi` 로 재현 가능 → **벤치마크 관측량 후보** |
| **S5** | `Ω_max` 평형 미세구조의 정전 에너지 지형 (log scale) | — |
| **S6** | rdf 비교 — **`Ω_hetero`(자기퍼텐셜 불균질) vs `a_hetero`(쌍상호작용 불균질)** 두 모델 | 불균질을 어느 층에 넣느냐의 비교. 재현 시 두 조건 다 필요 |
| **S7** | 4종 입자 ζ-전위 측정 | Table 1 `ψ` 의 원자료 |
| **S8** | gel trapping 접촉각 측정 이미지 | 위 §불일치 참조 |
| **S9** | **입자 크기에 따른 `a_ij` 의존성** (SPS·CPS) | `d` 를 바꿀 때 `a_ij` 를 어떻게 옮길지의 단서 |
| **S10–S12** | SPS·APS·CPS 각 4개 입자의 **AFM 표면 형상** (1×1 μm, non-contact phase) | 표면 거칠기 = 접촉선 요동의 물리적 근거. 모세관 사중극자 논의와 직결 |

> **S2 보충:** 크기 효과 검증에 쓴 입자는 `d ≈ 1.04, 5.05, 1.00 μm` (본 측정의 ~3 μm 와 별개).

**SI에 `η_eff` 추가 정보는 없다.** 산문은 S-1~S-3 뿐이고 전부 읽었다 — S-4 이후는 그림,
S-13은 참고문헌이다. `η_eff` 는 §5 의 원출처(Park & Furst 2008)로 가야만 얻을 수 있다.

## 5. ★ 계면 유효 점도 `η_eff` — 형태는 확정, 값은 미확보

SI S-3, *Optical trapping and drag calibration*. 계면에 포획된 입자를 잡고 스테이지를 등속
`u`로 움직여 크리핑 유동을 걸고, **Stokes 항력 ↔ 광집게 힘의 균형**에서 trap stiffness를 얻는다.

$$\kappa_t = \frac{3\pi d\,\eta_{eff}\,u}{\Delta x}$$

| 항목 | 값 | 근거 |
|---|---|---|
| **형태** | `F_drag = 3π d η_eff u` = **`6π R η_eff u`** — 표준 Stokes 항력에 `η → η_eff` 치환 | SI S-3 |
| ⇒ BD 항력계수 | **`γ = 3π η_eff d = 6π η_eff R`** | 위에서 유도 |
| `η_eff` 정의식 | ✅ **해소** — `[η_oil(1−cosθ) + η_water(1+cosθ)]/2` | [SI ref 2 확보 완료](2008-park-salt-surfactant-interface-forces.md) |
| `η_eff` 수치 | ⚠️ 원출처에도 수치는 없음 (표준 물성 대입) | 〃 |
| 접촉각 의존성 | ✅ **표면적 가중** — `θ` 는 물 상 기준 | 〃 |
| 계면활성제/`Bq` 보정 | ❌ **원출처도 표면점도를 무시한다** | 〃 |

> **`d`는 지름이다.** `3πd`가 낯설어 보이지만 `3π(2R) = 6πR` 로 표준형과 같다.
> 반지름으로 오독하면 항력이 2배 틀리고, `D = kT/γ` 를 통해 **시간축 전체가 2배 틀어진다.**

### ✅ 값의 원출처 — 확보 완료 (2026-07-27)

> **SI ref 2 · 본문 ref 48/52:**
> [**Park, B. J.; Pantina, J. P.; Furst, E. M.**; Oettel, M.; Reynaert, S.; Vermant, J.
> *Direct Measurements of the Effects of Salt and Surfactant on Interaction Forces between
> Colloidal Particles at Water−Oil Interfaces.* **Langmuir 2008, 24(5), 1686−1694.**
> ](2008-park-salt-surfactant-interface-forces.md)

$$\eta_{eff} = \frac{\eta_{oil}(1-\cos\theta) + \eta_{water}(1+\cos\theta)}{2}$$

**두 벌크 점도의 표면적 가중 평균**이고, `θ`는 물 상 기준이다. 상세·한계는 위 증류본 §1–2.
핵심 한계: **표면점도를 무시한다** — 계면활성제가 중요한 계에서는 재검토 대상.

`2023-lee-colloidal-debye-force`의 ref 22와 **같은 논문**이다 — 두 논문이 독립적으로 같은 곳을
가리켰고 실제로 정의가 거기 있었다. **추적 종료.** Park·Furst 공저이므로
[[goal-autonomous-paper-to-sim-verification]]의 두 축이 여기서 만난다.

관련(미확보): SI ref 5 = Park, B. J.; Furst, E. M. *Optical Trapping Forces for Colloids at the
Oil−Water Interface.* Langmuir 2008, 24(23), 13383−13392.

---

## 6. 우리 쪽 쓸모

**① 벤치마크 후보 — rdf.** §4의 MC 설정이 완전히 명세되어 있고(`N`·초기배치·스텝·사이클·
퍼텐셜·파라미터 전부) 실험 rdf와 대조까지 되어 있다. 재현 가능성이 높은 편.
`knowledge/wiki/benchmarks/`에 등재 검토 대상.

**② 파라미터 사전.** `Ω`·`a_ij`·`θ_c`·`ψ`·`d`가 오차와 함께 있다 (§3). 계면 콜로이드 계를
다룰 때 `U(r)` 크기의 실측 근거로 쓸 수 있다.

**③ [[goal-autonomous-paper-to-sim-verification]] 관점 — 첫 구체적 S0 후보.**
이 논문은 **구조(rdf)는 MC로 검증했지만 동역학은 시뮬레이션한 적이 없다.**
저자들이 브라운 확산을 의도적으로 무시했기 때문에(`U_ij ≈ 577 k_BT`) 애초에 `η_eff`가
시뮬레이션에 들어가지 않았다 — 실험 보정에만 쓰였다. 따라서:

> **BD로 이 계의 동역학(MSD·구조 완화·`Ω` 불균질이 확산에 미치는 영향)을 돌리는 것은
> 이 논문이 하지 않은 검증이다.** S0가 찾아야 할 후보의 첫 실례.

단, S0 필터 주의 — 이 계는 준2D라 HI가 살아 있다. no-HI BD로 정량 주장을 하려면
`master_plan.md` §12(HI 비목표)를 리포트에 명시해야 한다.

## 7. 증류 범위와 남은 것

**읽은 것:** 본문 8쪽 전체 + SI 13쪽 전체(S-1~S-3 산문, S-4~S-12 그림 캡션, S-13 참고문헌).
SI는 더 볼 것이 없다 — `η_eff` 추가 정보 없음을 확인했다.

| 항목 | 상태 |
|---|---|
| `η_eff` 정의식 | ✅ **해소** — [Park & Furst 2008](2008-park-salt-surfactant-interface-forces.md) 확보 완료 |
| Debye 길이 `κ⁻¹` | ❌ 순수 물 계라 이온세기 명시 없음 — 원문에 부재 |
| 온도 `T` | ❌ 명시 없음 (ambient 추정) — 벤치마크 등재 전 확인 |
| rdf 정규화 상수 | ⚠️ 본문 수식이 PDF 텍스트 추출에서 깨짐. 벤치마크 등재 시 그림/원문 재확인 |
| `θ_c(APS)` vs `θ_c(CPS)` | ⚠️ 원문 내부 불일치. **정본 = Table 1** 로 고정 (§3 참조) |
| 표면 거칠기 정량값 | ⚠️ AFM 이미지(S10–S12)만 있고 RMS 등 수치는 미기재 — 필요하면 이미지에서 재추출 |
