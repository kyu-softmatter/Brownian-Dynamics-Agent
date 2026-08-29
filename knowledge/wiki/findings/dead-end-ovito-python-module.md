---
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
question: "OVITO Python 모듈을 `hoomd_slit` env 에 넣을 수 있는가? (`D14`)"
cites:
  - docs/00_decision_log.md
  - simbot/environment.yaml
related_systems: []
---

# dead-end — OVITO Python 모듈은 `hoomd` 와 같은 env 에 넣지 않는다

**질문:** 입자 렌더러로 OVITO 를 쓰려면 Python 모듈이 필요하다. `hoomd_slit` env 에 넣을 수 있는가?

**답: 넣을 수 있지만 넣으면 안 된다.** `tbb` 를 내리고, `hoomd` 는 그 아래 남는다.

---

## 함정 ① — conda-forge 의 `ovito` 는 GUI 앱이다

```bash
conda install -n hoomd_slit -c conda-forge ovito      # 77개 패키지, 성공. 그런데
python -c "import ovito"                              # ModuleNotFoundError
```

`bin/ovito` (실행파일)만 들어오고 `site-packages` 에는 아무것도 없다. ffmpeg · hdf5 ·
libopenvino · aws-c-* 로 77개를 끌고 오는데 **헤드리스 워크플로에서 쓸 수 없는 GUI 앱**이다.

**신호는 빌드 문자열에 있었다.**

```
conda-forge          ovito 3.15.5  h999aba6_1      ← py 없음 = Python 모듈 아님
conda.ovito.org      ovito 3.15.5  py312hf032ba3_0 ← 이게 모듈
```

`--dry-run --json` 을 이미 돌려놓고도 `h999aba6_1` 을 그냥 지나쳤다. **패키지 이름이 같아도
빌드 문자열이 다른 물건이다.** `py3XX` 접두가 없으면 Python 모듈이 아니다.

되돌리기는 깨끗했다 — `conda remove -n hoomd_slit ovito` 가 정확히 그 77개를 가져가고
`hoomd`·`gsd`·`numpy`·`freud`·`fresnel`·`tbb` 는 건드리지 않았다 (`LINK: 0`).

## 함정 ② — 진짜 모듈은 `tbb` 를 내린다

OVITO 자체 채널의 `py312` 빌드를 `--freeze-installed --dry-run` 으로 재보면:

```
- tbb       2023.0.0  →  + tbb       2022.3.0      다운그레이드
- fresnel   0.13.8    →  + fresnel   0.13.8        떼었다 붙임 (embree·tbb 재링크)
- qt6-main  6.11.1 → 6.10.2 · libjxl 0.12.0 → 0.11.2 · lame 4.0 → 3.100 · svt-av1 4.2.0 → 4.0.1
UNLINK 25 · LINK 35
```

**`hoomd 7.1.0` 은 UNLINK 목록에 없다.** 즉 `tbb` 가 바뀌어도 hoomd 는 재빌드·재링크되지 않고
바뀐 `tbb` 위에 그대로 남는다. 이 저장소의 전제가 `hoomd 7.1.0` + `numpy 2.5.1` 기준선의
재현성(`D34` · `machine_profile`)이므로, 렌더러 하나를 얻으려고 엔진의 스레딩 라이브러리를
내리는 것은 교환이 성립하지 않는다.

**`--freeze-installed` 는 이걸 막지 못했다.** 그 플래그는 명시한 spec 을 고정할 뿐,
전이 의존성의 다운그레이드는 막지 않는다. 그래서 `--dry-run --json` 의 `UNLINK` 를
직접 읽어야 한다 — 요약 출력만 보면 "성공" 으로 보인다.

## 그래서 무엇을 잃었나 — 사실상 없다

OVITO 를 원한 이유는 렌더와 분석 모디파이어 둘이었고, 둘 다 대체가 있다.

| 원한 것 | 대체 | 확인 |
|---|---|---|
| 3D 입자 렌더 | `fresnel 0.13.8` | macOS 헤드리스 `pathtrace` 성공 (320×320, 구 60개, embree CPU) |
| 2D 단층 렌더 | matplotlib | 이미 쓰고 있다. 2D 에는 음영 줄 깊이가 없어 레이트레이싱이 주는 것이 없다 |
| Voronoi · ψ₆ · CNA | `freud 3.5.0` | `Hexatic` · `Voronoi` · `StaticStructureFactorDirect` · `DiffractionPattern` |

## 다시 필요해지면

**별 env 에 두고 GSD 만 넘긴다.** `hoomd_slit` 을 건드리지 않는 것이 조건이다.

```bash
conda create -n ovito_viz -c https://conda.ovito.org -c conda-forge python=3.12 "ovito=*=py312*"
```

렌더는 궤적 파일만 읽으므로 env 를 공유할 이유가 애초에 없다. 같은 env 에 넣으려 한 것이
실수였고, 비용은 77개 패키지의 왕복이었다.

## 규칙 후보

`silent-default-is-a-lie` 계열의 인프라 판본이 하나 더 생겼다 —
**"conda 가 `success: true` 라고 해도 `UNLINK` 를 읽지 않았으면 모르는 것이다."**
`SD5`(GSD `dynamic`) · `SD13`(`write_at_start`) · `SD14`(로그 간격)와 같은 부류다:
기본값·요약이 조용히 틀린 답을 준다.

See also: [d23-sdk-backend](d23-sdk-backend.md) · `docs/00_decision_log.md` `D14`
