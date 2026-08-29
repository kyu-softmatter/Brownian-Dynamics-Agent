---
type: technique
author: agent
drafted: 2026-07-28
confirmed_by:
question: "simulation_bot env에 어떤 패키지를 언제 왜 넣었는가?"
---

# 패키지 적립 이력 — `simulation_bot`

> **append-only.** 과거 행을 고치지 않고 `정정:` 으로 새 행을 추가한다.
> 재현: `conda env create -f environment.yml`

env 위치: `/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot`
인터프리터: `.../envs/simulation_bot/bin/python` — **절대경로로 호출한다**
(`conda activate`는 non-interactive shell에서 불안정)

---

## 2026-07-28 · 1단계 코어

```
python=3.12  hoomd  gsd  numpy  scipy
```

| 패키지 | 버전 | 왜 필요했나 |
|---|---|---|
| `python` | 3.12 | hoomd 7.1의 conda-forge 빌드가 있는 버전 중 안정 |
| `hoomd` | **7.1.0** (`cpu_py312_hf901771_0`) | BD 엔진. `md.methods.Brownian` |
| `gsd` | 5.0.1 | 궤적 I/O (HOOMD 네이티브 포맷) |
| `numpy` | 2.5.1 | 수치 |
| `scipy` | 1.18.0 | 피팅·통계·특수함수(`Γ`) |

### ★ 이 단계에서 확인된 제약 — 기록해둘 가치가 있다
```
hoomd.version.mpi_enabled = False
CPU util (프로세스당)      = 1.00x   → 완전 단일스레드
conda-forge MPI 빌드       = osx-arm64 0건, linux-64 0건
```
> **HOOMD의 MPI 빌드는 conda-forge에 아예 없다.** 플랫폼 문제가 아니다.
> v3+에서 TBB도 제거되었으므로 **CPU 병렬화 경로가 없다.**
> ⇒ 처리량 병렬화(독립 런 동시 실행)만 가능. 자세히 [[local-cpu-parallelism]]

기존 `hoomd_slit` env는 건드리지 않았다 (거기에도 hoomd·freud·fresnel이 있지만
버전 고정과 이력 관리를 새로 시작하기 위해 별도 env로).

## 2026-07-28 · 2단계 분석·플롯·개발

```
pyyaml  matplotlib  freud  pandas  h5py  pillow  pytest
```

사용자 지시: **"HOOMD7 충돌이 없는 선에서 최신 버전으로."** 전부 conda-forge 최신으로 설치했고
설치 후 `hoomd.Simulation` 생성이 정상 동작함을 확인했다.

| 패키지 | 버전 | 왜 필요했나 |
|---|---|---|
| `pyyaml` | 6.0.3 | `config/run_policy.yaml`, `spec.yaml` 직렬화. **1단계에서 빠뜨려서 즉시 막혔다** |
| `matplotlib` | 3.11.1 | S6 플롯 전부. **그림 글자는 영문** (한글 글리프 없음) |
| `freud` | 3.5.0 | RDF·S(q)·클러스터·보로노이. 직접 구현보다 빠르고 검증됨 |
| `pandas` | 3.0.5 | 측정값 표·`params` 원장 |
| `h5py` | 3.16.0 | `thermo.h5` — HOOMD `write.HDF5Log` |
| `pillow` | 12.3.0 | 손그림 이미지 메타·EXIF. **첫 예시에서 orientation=6을 잡아냈다** |
| `pytest` | 9.1.1 | 회귀 테스트. 벤치마크의 실행 주체 |

### 아직 넣지 않은 것과 이유
| 패키지 | 언제 넣나 |
|---|---|
| `fresnel` | 3D 레이트레이싱 스냅샷이 필요해질 때. v1 첫 케이스는 2D라 불필요 |
| `ffmpeg` | GIF로 부족할 때(MP4). 용량이 크므로 미룬다 |
| `pymbar` | 평형 판정을 자동화할 때. 단 능동계엔 무의미하므로 카드별로 켜고 끈다 |

---

## 2026-07-28 · 3단계 논문 원본 판독

```
pypdf==6.14.2      # pip (conda 아님)
```

| 패키지 | 버전 | 왜 필요했나 |
|---|---|---|
| `pypdf` | 6.14.2 | **논문 PDF에서 텍스트를 뽑을 경로가 하나도 없었다.** `knowledge/source/papers/` 증류를 만들려면 원본을 읽어야 하는데, env에 PDF 파서가 없고 시스템에도 `pdftotext`·`mutool`·`qpdf` 전부 없었다 |

### 왜 이 선택이었나 — 기각한 대안들
| 대안 | 기각 이유 |
|---|---|
| Claude Code 내장 PDF 리더 | `pdftoppm`(poppler) 의존. 시스템에 없어서 실패 |
| `brew install poppler` | conda env 밖의 시스템 변경. `pdftotext`+`pdftoppm`을 동시에 얻어 **2단 조판 추출 품질이 pypdf보다 낫고 내장 리더도 살아난다**는 장점은 있다 → 그림·수식을 봐야 할 때 재검토 |
| `pymupdf` | 텍스트+렌더링 둘 다 되지만 바이너리 의존이 크다. 텍스트만 필요한 지금은 과하다 |

`pypdf`는 순수 파이썬·무의존이라 HOOMD와 충돌 가능성이 없다.
설치 후 `import hoomd` → `7.1.0` 정상 확인했다.

### ★ 한계 — 알고 쓸 것
2단 조판 저널 PDF(PRL 등)에서 **컬럼이 섞이고 특수문자가 `/.0020` 같은 글리프 코드로 깨진다**
(실측: `κ` → `/.0020`, `μ` → `/.0022`, `=` → `/.0136`). 숫자와 영문은 살아 있어서
값 추출은 되지만 **수식을 그대로 신뢰하면 안 된다.** 계수·지수는 사람이 눈으로 확인할 것.
그림은 아예 못 읽는다 → Fig.에서 값을 읽어야 하면 poppler 또는 이미지 렌더링이 필요하다.

---

## 2026-07-29 · 4단계 논문 **그림·수식** 판독 (poppler)

```
brew install poppler        # 26.07.0 — conda 아님, 시스템 설치
```

| 도구 | 왜 필요했나 |
|---|---|
| `pdftotext -layout` | 2단 조판에서 컬럼이 안 섞인다. `pypdf` 는 섞는다 |
| **`pdftoppm`** | ★ **본론.** 페이지를 PNG 로 렌더링 → **수식과 그림을 눈으로 읽는다.** Claude Code 내장 PDF 리더도 이걸 요구한다 |
| `pdfimages` | (아직 안 씀) 그림만 뽑을 때 |

### ★ 왜 3단계(`pypdf`)로 부족했나 — 실제로 막혔다
`pypdf` 도 `pdftotext` 도 **수학 기호 글리프를 잃는다.** 저널 PDF 의 심볼 폰트에
`ToUnicode` 가 없어서다. 실측:

```
pypdf     :  "0  3a4c E=4a3"        ← κ₀ = 3π a_c⁴E/4a³  (κ·π·= 소실)
pdftotext :  "0  3a4c E=4a3 ,"      ← 같은 손실
렌더링 이미지: "κ₀ = 3π a_c⁴E/4a³"    ← 읽힌다
```

**계수를 텍스트 층에서 읽으면 안 된다.** `Eq. (1)` 의 `Lx²/4 − |x³|/6` 에서 `4` 와 `6` 은
살아남았지만 `π`·`κ` 는 사라졌고, 어느 것이 사라졌는지 텍스트만 봐서는 알 수 없다.

그리고 **그림은 텍스트 층에 아예 없다.** `Fig. 3(a)` 의 자료점과 inset 이
`κ₀` 의 규약(`(a/s)³` vs `(2a/s)³`)을 확정한 결정적 증거였다 —
이건 렌더링 없이는 접근 불가였다.
근거: [[2005-pantina-furst-bending-coefficient]]

### 절차 (다음에 논문 정독할 때)
```bash
pdftotext -layout paper.pdf out.txt          # 산문·수치는 여기서
pdftoppm -r 300 -f <p> -l <p> -png paper.pdf out   # 수식·그림은 여기서
# 그 뒤 PIL 로 crop 해서 확대해 읽는다 (전체 페이지는 수식이 작다)
```

`conda` env 밖의 시스템 설치이므로 `environment.yml` 에는 주석으로만 남긴다.

---

## 교훈
- **1단계에 `pyyaml`을 빼먹었다.** "코어"에 직렬화를 안 넣어서 설정 파일을 쓰는 순간 막혔다.
  → 다음 env를 만들 때: 직렬화(`pyyaml`)와 테스트(`pytest`)는 코어에 포함시킬 것.
- conda 단발 solve가 여러 번 solve보다 훨씬 빠르다. "하나씩 적립"은
  **설치 순서**가 아니라 **기록의 단위**로 해석하는 것이 맞다.
