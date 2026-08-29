---
type: finding
author: agent
drafted: 2026-07-29
confirmed_by:
system: any
dynamics: any
cites:
  - simbot/io.py
  - simbot/run.py
  - runs/2026-07-29_soft-r3-time-resolved/metrics.json
  - tests/test_s8_io.py
---

# provenance 는 **정의가 하나**여야 하고 **포착 지점은 셋**이다

> `report.py` 는 `env_hash` 를 읽는데 두 러너가 자기 manifest 를 손으로 만들면서
> 그 키를 빼먹었다. **리포트의 재현 정보가 빈칸으로 렌더되고 있었고, 아무도 몰랐다.**

## 무엇이 일어났나

`simbot.io.build_manifest()` 가 정본 provenance 빌더였고 `code_hash`·`git_rev`·
`git_dirty`·`env_hash`·`env` 를 전부 담고 있었다. `cli.py` 는 그것을 썼다.
그런데 **`run_trap` 과 `run_soft2d` 는 각자 manifest dict 를 손으로 만들었고**
`env_hash`·`env` 를 넣지 않았다. `report.reproducibility_section` 은 그 키를 읽으므로
러너로 만든 런의 리포트에서는 그 줄이 **조용히 사라졌다.**

고치려다 같은 실수를 반복했다: `metrics.json` 용 provenance 를 새로 만들면서
키 이름을 `code_hash_simbot` 으로 지었다 — **읽는 곳이 0개인 어휘**였다.

## 규칙 ① — 정의는 하나

```python
# simbot/io.py — provenance() 가 유일한 정의
def provenance(driver=None) -> dict:
    return {"code_hash": ..., "git_rev": ..., "git_dirty": ...,
            "env_hash": ..., "env": ..., **driver_keys}

# build_manifest · run_trap · run_soft2d · 분석 스크립트가 전부 이것을 쓴다
manifest = {"run_hash": cfg.hash(), **provenance(), ...}
```

호출처마다 손으로 나열하면 **키 이름이 갈라지고, 어긋남이 예외가 아니라 빈칸으로
나타난다.** 빈칸은 "정보 없음" 과 구별되지 않으므로 아무도 알아채지 못한다.

감시: `tests/test_s8_io.py::test_provenance_supplies_every_key_the_report_renders`
와 `::test_build_manifest_gets_provenance_from_the_single_definition`.

## 규칙 ② — 포착 지점은 셋이고, 서로를 대체하지 못한다

| 지점 | 어디에 | 무엇을 답하는가 | 없으면 |
|---|---|---|---|
| **봉인 시점** (실행 전) | `03_spec.yaml` (`SEALED.sha256` 대상) | "무엇을 하겠다고 **약속**했나" | 사후합리화를 막을 수 없다 |
| **궤적 생성 시점** | `raw/<run>/manifest.json` | 이 **궤적**을 무엇이 만들었나 | 궤적을 특정할 수 없다 |
| **분석 시점** | `metrics.json` 의 `_provenance_at_analysis` | 이 **측정값**을 무엇이 만들었나 | 아래 참조 |

### 분석 시점이 왜 따로 필요한가 — 실측

분석은 궤적과 **따로 돌 수 있다** (`--analyze-only`). 2026-07-29 세션에서
`metrics.json` 을 4번 재생성했는데 그때마다 `simbot` 의 `code_hash` 가 달랐다:

```
60b7adc8932e → 911d3274ba33 → 0417e17e46fa → 361096e5fddc → 16cf434faf67
```

그동안 `raw/*/manifest.json` 은 계속 `60b7adc8932e` 를 가리켰다 —
**궤적 시점 해시로는 측정값을 특정할 수 없다.** 가설이 아니라 한 세션에서 5번 일어났다.

### `freud` 가 빠지면 안 되는 이유

`hoomd` 만 기록하는 것으로는 부족하다. **`freud` 가 Voronoi·`ψ₆` 를, `scipy` 가
완화시간 적합을 계산한다.** 같은 궤적에서 `freud` 버전이 바뀌면 결함 분율이 바뀔 수
있고, 기록이 없으면 그 사실을 사후에 알 수 없다.
`ENV_PACKAGES = ("hoomd", "numpy", "scipy", "gsd", "freud")`.

## 규칙 ③ — `code_hash` 는 `simbot/` 만 덮는다

런의 `A` 목록·시드·런 길이·**분석 창**을 정하는 것은 `scripts/` 의 드라이버다.
`code_hash` 에 그것이 없으므로 `driver_hash` 를 따로 잡는다.

```
code_hash    = sha256(simbot/**/*.py)      ← 엔진
driver_hash  = sha256(scripts/<driver>.py) ← 이 런의 파라미터와 분석 창
```

⚠ **`driver_hash` 는 파일 1개의 해시다.** 드라이버가 `scripts/` 안의 다른 모듈을
import 하면 덮이지 않는다. 현재 드라이버들은 `simbot` 만 import 하므로 둘이 합쳐
전부를 덮는다 — **그 전제가 깨지면 이 문장이 틀린다.**

감시: `::test_driver_hash_changes_with_content_but_code_hash_does_not`

## 규칙 ④ — 봉인 파일이 `.gitignore` 에 걸려 있지 않은지 확인한다

`.gitignore` 가 `runs/**` 를 제외하고 `.md`/`.json`/`.yaml` 만 재포함하고 있었다.
`SEALED.sha256` 은 확장자가 달라서 **전 런에서 한 번도 커밋되지 않았다** (실측: 0개,
`CLAUDE.md` 가 첫 완주로 인용하는 `2026-07-28_trap-2d-5um_2dfb9d` 까지 포함).

결과가 조용하다:
- `REPORT.md` 는 독자에게 `shasum -a 256 -c SEALED.sha256` 을 안내하는데
  **클론에는 그 파일이 없다.**
- 예측 `.md` 는 커밋되므로 해시를 다시 계산할 수는 있다. 그러나 봉인의 요점은
  **실행 전에 기록된 해시**이고, 그것이 히스토리에 없으면 "나중에 고치지 않았다" 를
  히스토리로 주장할 수 없다.

⇒ `!runs/**/SEALED.sha256` 을 재포함한다. **텍스트 산출물을 확장자로 필터링하면
봉인처럼 확장자가 다른 1급 산출물이 빠진다.**

## 그래도 provenance 는 재현을 **증명하지 않는다**

해시는 "무엇이 만들었나" 를 기록할 뿐 "다시 돌리면 같은가" 를 증명하지 않는다.
직접 증거는 재실행 비교뿐이다:
`tests/test_s5_pair.py::test_same_config_reproduces_the_trajectory_bit_for_bit`
(같은 `Soft2DRunConfig` 두 번 → `traj`·`energy`·`max_force`·`init_pos` 바이트 동일).

그리고 **커밋되지 않은 상태(`git_dirty: true`)면 `git_rev` 는 아무것도 특정하지 않는다.**
2026-07-29 에 작업 트리가 되돌아가는 사고가 있었고, `runs/`(gitignored) 만 살아남아
**산출물은 있고 그것을 만든 코드는 없는** 상태가 실제로 발생했다.
→ 런을 결론에 쓰려면 **먼저 커밋한다.**

## 함께 볼 것

- [[coarse-sampling-hides-the-whole-transient]] · [[tolerance-from-a-4-seed-se-is-not-a-3-sigma-test]]
  — 같은 런에서 나온 다른 실패들
