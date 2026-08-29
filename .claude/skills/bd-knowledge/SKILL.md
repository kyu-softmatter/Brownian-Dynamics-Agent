---
name: bd-knowledge
description: knowledge/ 지식 베이스를 검색·추가·정리한다. "이 파라미터 근거가 뭐야", "논문에서 이 값 찾아줘", "이걸 기록해둬", "카드 만들어줘" 라고 할 때, 또는 새 (계×동역학) 카드·finding·benchmark 를 만들 때 쓴다. 시뮬레이션 실행은 bd-pipeline, 실패 진단은 bd-diagnose.
---

# bd-knowledge — 지식 검색·추가·정리

계약: [`knowledge/wiki/CLAUDE.md`](../../../knowledge/wiki/CLAUDE.md) ← **먼저 읽는다**

> **이 위키는 일반 문헌 저장소가 아니다.** 두 가지 일만 한다:
> ① `pytest` 회귀의 **검증 오라클** 공급 ② 출처 붙은 **파라미터 사전** 공급.
> 읽기용 컨텍스트로만 쓰이는 페이지는 여기 있을 이유가 약하다.

## 구조

```
knowledge/
├── source/papers/   논문별 증류 42편 + INDEX.md      (원본 PDF 는 raw/, gitignored)
└── wiki/
    ├── systems/     ★ (계 × 목적동역학) 카드 — 무차원화·게이트를 소유한다
    ├── findings/    Q→A + dead-end-<slug>.md
    ├── concepts/    WHAT-IS (무차원수, 상거동, 퍼텐셜)
    ├── techniques/  HOW-TO (평형화 판정, 오차막대, env 이력)
    ├── benchmarks/  검증 오라클
    └── questions/   아직 답 없는 것. 삭제하지 않고 status 로 닫는다
```

## 검색

```bash
ls knowledge/wiki/systems/                      # 이 계의 카드가 있는가
grep -rl "<키워드>" knowledge/wiki/
grep -rl "<키워드>" knowledge/source/papers/
sed -n '1,20p' knowledge/source/papers/INDEX.md
```

**frontmatter 를 먼저 읽는다** — `status`, `reproduced`, `confirmed_by` 가 그 페이지를
근거로 쓸 수 있는지 결정한다.

## 인용 규율 — 어기면 위키가 소문 저장소가 된다

| 표기 | 의미 | 언제 |
|---|---|---|
| `[출처]` | **검증된 근거** | 문헌 벤치마크 또는 `reproduced: yes` |
| `[출처, 미재현]` | **사실 기록** | `reproduced: no` — 참고는 하되 검증 주장에 쓰지 않음 |

- **문헌값을 기억으로 인용하지 않는다.** `source/papers/` 의 증류를 인용한다.
  "학습 데이터에서 봤다"는 근거가 아니다
- **`reproduced: no` 를 근거처럼 쓰지 않는다.** 논문에 실렸다는 것이 우리 코드에서
  그 값이 동작한다는 뜻은 아니다
- 예: `N=1000` 을 앙상블로 쓰는 관행은 Barakat 2022 인데 `reproduced: no` 다 →
  `[출처, 미재현]`

## 새 (계 × 목적동역학) 카드 만들기

**카드 없는 쌍을 만나면 즉흥 무차원화 금지.** `simbot.nondim` 이 예외를 던진다.

```bash
cp knowledge/wiki/systems/_TEMPLATE.md \
   knowledge/wiki/systems/<계>--<동역학>.md
```

frontmatter:

```yaml
type: system
system: passive-sphere | abp | attractive-colloid | brush-colloid | interfacial-colloid
dynamics: equilibrium-structure | transport | harmonic-trap | dense-collective | coarsening
status: draft            # draft → usable → validated. 승격은 실측 근거가 있을 때만
```

**카드가 소유하는 것** (다른 곳에 적지 않는다):

| § | 내용 | 왜 카드가 소유하나 |
|---|---|---|
| 3 | **기준 단위** (길이·에너지·시간) | 같은 계라도 목적동역학에 따라 다르다 |
| 4 | 무차원수 원장 | 이 쌍에서 의미 있는 것만 |
| 6 | 관측량 | |
| 7 | **적용 게이트 — 켜고 끄기** | 능동계에 평형화 판정을 걸면 통과 불가 |
| 8 | 벤치마크 | `pytest` 승격 후보 |
| 10 | 남은 빈칸 | 닫힌 항목은 지우지 않고 표시만 바꾼다 |

카드를 만들면 **코드에도 등록한다**:

```python
# simbot/nondim.py
CARD_SCALE_RULES = {
    "<계>--<동역학>": "harmonic_trap" | "brownian" | "active_run_length",
}
```

등록하지 않으면 `scales_for` 가 예외를 던진다 — **의도된 동작이다.**

## 새 finding

`findings/_TEMPLATE.md` 를 쓴다. **진단 경로가 이 문서의 핵심 가치다** —
다음에 같은 증상을 만나면 그 순서를 그대로 따라간다.

`재발 방지`에 **"없음"을 적는 것이 반복되는 유형은 가드를 만들 신호다.**

## 새 benchmark → pytest 승격

```markdown
| # | 벤치마크 | 예상값 | 상태 |
|---|---|---|---|
| B1 | 무차원 등분배 · EM 편향 | `⟨x*²⟩ = 1/(1−dt*/2)` | `[O]` 측정 `1.01041±0.00264` (`0.1σ`) |
```

`tests/` 로 올릴 때:

```python
@pytest.mark.benchmark
def test_B1(...):
    """허용오차는 **이론 통계오차**에서 뽑는다.
    관측값에 맞춰 재단하면 검증이 아니라 사후합리화다."""
```

★ **문서 값과 대조할 때는 문서에 적힌 유효숫자에서 허용오차를 뽑는다.**
전역 상수 하나면 느슨한 쪽이 엄격한 값을 봐준다
(`tests/test_s3_spec.py::rounding_halfwidth`).

★ **경쟁 가설도 기각한다.** "예측과 맞다"만으로는 약하다. 단,
설계 검정력이 `3σ` 를 못 만드는 곳에서 `3σ` 기각을 요구하지 않는다 —
그 구간은 `INCONCLUSIVE` 를 사실로 고정한다.

## 저작 대칭성과 승격

```yaml
author: agent | human | hybrid
confirmed_by:              # 비워둔다. 사람이 검토한 뒤 채운다
```

- **`author: agent` 비율 자체가 자기개선 지표다**
- **`finding` → `concept` 승격은 항상 사람 승인.** 에이전트가 스스로 개념을
  인플레이션시키는 것을 막는다
- 저자와 무관하게 품질 기준은 같다: 근거가 `source/` 나 확립된 문헌으로 추적되고,
  인용이 구체적 파일·URL 이며, 주장이 확인 가능할 것

## 모순이 생겼을 때

**기존 항목을 조용히 덮어쓰지 않는다.** 새 항목을 만들고 `supersedes: [<이전 id>]` 로
연결한다. 적용 범위만 좁히는 것이면 `supersedes` 가 아니라 상호 링크다 —
[`displacement-gate-is-1000x-loose-for-traps`](../../../knowledge/wiki/findings/displacement-gate-is-1000x-loose-for-traps.md)
가 [`dt-gate-should-be-displacement-based`](../../../knowledge/wiki/findings/dt-gate-should-be-displacement-based.md)
를 **기각하지 않고 범위를 좁힌** 예다.

## 환경이 바뀌면

패키지를 추가하면 **왜 필요했는가**를
[`techniques/env-log.md`](../../../knowledge/wiki/techniques/env-log.md) 에 적고
`environment.yml` 도 갱신한다.
