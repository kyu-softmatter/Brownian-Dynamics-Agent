---
type: finding
author: agent
drafted: 2026-07-27
confirmed_by:
question: "`agent/` 층은 LLM을 어떤 경로로 불러야 하는가? (`D23`)"
cites:
  - master_plan.md
  - agent/llm.py
related_systems: []
---

# `D23` — SDK 직접 호출로 확정. 구독 경로는 조직 정책으로 막혀 있다

**질문:** `agent/`가 모델을 부르는 경로를 (a) Anthropic SDK 직접 / (b) Claude Code `claude -p` /
(c) Agent SDK 중 무엇으로 할 것인가?

**답: (a) SDK 직접.** 문서의 원래 기본값과 같지만, **한 번 뒤집었다가 되돌아왔다.** 그 과정이
이 finding의 내용이다.

---

## 뒤집으려 했던 이유 (틀렸다)

환경을 조사하니 이랬다.

| | |
|---|---|
| `anthropic` 패키지 | 미설치 |
| `ANTHROPIC_API_KEY` | 없음 |
| `claude` CLI | **설치됨** (`~/.local/bin/claude`) |

여기서 *"(b)는 오늘 당장 돈다. 의존성도 API 키도 필요 없다"* 고 판단해 `cli`를 기본 백엔드로
두었다. 이 저장소의 "표준 라이브러리만 쓴다" 관례와도 맞아 보였다.

## 실제로 돌려보니 (2026-07-27)

```
$ claude -p --output-format json -- "hi"
{"is_error": true,
 "result": "Your organization has disabled Claude subscription access for Claude Code
            · Use an Anthropic API key instead, or ask your admin to enable access"}
```

**어차피 API 키가 필요하다.** 조직이 Claude Code의 구독 접근을 막아두었기 때문이다.

## 그래서 (a)로 되돌아온다

`cli`의 유일한 장점이 *"키 없이 구독으로 돈다"* 였는데 그게 사라졌다. 남는 것은 단점뿐이다.

| | (a) SDK | (b) CLI |
|---|---|---|
| API 키 | 필요 | **필요** ← 여기가 뒤집혔다 |
| 스키마 강제 | **tool-use로 구조적 강제** | 불가. 프롬프트의 예시가 전부 |
| 테스트 | 함수 1개, mock 주입 | 서브프로세스, 비결정론 |
| 이미지 | base64 직접 첨부 | Read 도구 경유 |

스키마 강제가 특히 크다. `SketchReading`은 필드마다 `source`·`where`·`confidence`를 요구하는
중첩 구조라 프롬프트만으로는 형식 이탈이 잦다. tool-use는 **파싱 실패라는 실패 모드 자체를
없앤다.**

## 무엇을 배웠는가

**환경 조사만으로 결정하면 안 된다. 실제로 한 번 불러봐야 한다.** "CLI가 설치되어 있다"와
"CLI가 동작한다"는 다른 명제였고, 둘을 구분하는 비용은 호출 한 번이었다.

`D23`을 *"M1에서 S1 INTAKE를 실제로 구현할 때"* 결정하기로 미뤄둔 것이 맞았다. 설계 시점에
정했다면 이 사실을 모른 채 정했을 것이다.

## 부수 발견 — 오류 꼬리만 자르면 원인이 날아간다

처음 `_call_cli`는 실패 시 출력의 **마지막 600자**를 보여줬다. 그랬더니 이렇게 나왔다.

```
claude -p exit 1: ur admin to enable access","stop_reason":"stop_sequence", ...
```

`claude -p`는 진짜 이유를 **JSON 봉투의 `result` 필드 맨 앞**에 넣는데, 뒤에 긴 usage 블록이
붙어서 꼬리를 자르면 정확히 그 부분이 날아간다. 종료코드와 무관하게 **봉투를 먼저 파싱**하도록
고쳤다 (`agent/llm.py`).

일반화하면 — **오류를 자를 때는 뒤가 아니라 앞을 남긴다.** 구조화된 오류는 앞이 원인이고
뒤가 부속물이다.

## 되돌리는 조건

조직이 Claude Code 구독 접근을 열어주면 (b)가 다시 후보가 된다. 다만 그때도 스키마 강제와
테스트 용이성 때문에 (a)가 유리하고, **S2 ELICIT(대화 왕복)만** (b)로 빼는 원래 설계가
그대로 성립한다. 백엔드 이음새는 `agent/llm.py` 하나에 있으므로 교체 비용은 여전히 낮다.
