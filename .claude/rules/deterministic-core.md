# deterministic-core — 코어는 LLM을 부르지 않는다

`simbot/` 과 `bdkit/` 에 LLM 의존을 넣지 않는다. LLM 호출은 `agent/llm.py` **한 곳**에서만
일어나고, `cli.py` 가 둘을 엮는다. 의존 방향은 한쪽이다 — `agent → bdkit` 은 되고
`bdkit → agent` 는 안 된다.

```
pytest tests/test_invariants.py     # 이것이 검사다
```

**`grep` 으로 검사하지 않는다.** 원래 `A4` 는 `grep -rE "anthropic|claude" simbot/ bdkit/`
로 적혀 있었는데, 그 검사는 **규칙을 설명하는 산문 자체에 걸린다** — 이 파일을 인용하는
docstring 이 늘어날 때마다 걸린다 (2026-07-28 실측 7건). 문자열 검사는 코드와 산문을
구분하지 못한다. `tests/test_invariants.py` 는 AST 를 파싱해 **실제 임포트만** 본다.

**Why (the triggering incident):** 사고가 아니라 아키텍처 결정이다 (`A4` · `D1` · `SD2`,
2026-07-27). 이유는 하나뿐이고 그게 전부다 — **결과가 틀렸을 때 물리가 틀린 건지 LLM이
틀린 건지 구분할 수 있어야 한다.** 시뮬레이션에서 이건 치명적이다: 그럴듯하지만 틀린
`g(r)` 은 눈으로 구분이 안 되고, 채점자가 없는 도메인이라 스스로 틀렸다는 걸 알 방법도 없다.

이 분리가 실제로 값을 한 사례는 있다. `SQ6` — 트랩 자기상관 `τ` 가 이론보다 70% 높게
나왔을 때 후보가 둘이었다: 시뮬레이션이 틀렸나, 추정기가 틀렸나. 코어가 결정론이라
**정답을 아는 합성 OU 데이터를 넣어보는 것**으로 갈랐다 (추정기였다, `SD12`).
코어에 LLM 이 섞여 있었으면 후보가 셋이 되고, 셋째는 재현조차 안 된다.

**How to apply:**
- `simbot/` `bdkit/` 에 무엇을 더할 때 `import anthropic` 이 필요해지면 **위치가 틀린 것**이다.
  그 함수는 `agent/` 로 간다
- `agent/` 의 산출물은 항상 **스키마 고정 + 결정론 검증기 통과** 뒤에 코어로 들어간다.
  LLM 이 낸 값이 검증 없이 config 에 도달하는 경로를 만들지 않는다
- 새 분석 추정기는 **정답이 있는 합성 데이터**로 시험한다 (`tests/test_simbot_estimators.py`,
  `tests/test_structure.py`). 이건 물리 검증이 아니라 "코드가 답을 아는 입력에 답을 내는가"다
- `tests/test_invariants.py` 가 이 규칙을 강제한다. 끄지 않는다

**Anti-patterns explicitly forbidden:**
- **편의 임포트** — "여기서 한 번만 부르면 되는데" 로 코어에 LLM 을 들이는 것
- **검증기 우회** — LLM 제안을 "이번엔 확실하니까" 로 검증기 없이 적용하는 것
- **세 번째 후보 만들기** — 물리·추정기 외에 재현 불가능한 원인을 디버깅 경로에 추가하는 것

See also: [axioms](axioms.md) · `docs/00_decision_log.md` `D23`
