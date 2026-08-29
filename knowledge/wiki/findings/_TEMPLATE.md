---
id: <kebab-case-id>
kind: failures
tags: []
created: YYYY-MM-DD
updated: YYYY-MM-DD
runs: []                    # 이 실패가 발생한 run_id
confidence: high            # 증상-원인 인과가 확인되면 high
supersedes: []
cause_class: numerical      # numerical | modeling | interpretation | analysis | environment
stage: S5                   # S1..S8
---

## 증상
<관측된 것을 그대로. 에러 메시지 전문, 그래프 모양, 이상한 숫자값.
 "이상했다"가 아니라 "MSD가 t^2로 자랐다", "300스텝에서 NaN" 처럼 구체적으로.>

## 진단 경로
<무엇을 의심하고 어떻게 배제했는가. 순서대로.
 다음에 같은 증상을 만나면 이 순서를 그대로 따라간다. 이게 이 문서의 핵심 가치다.>

1.
2.
3.

## 근본 원인
<한 문단. `cause_class`와 일치해야 한다.>

## 처방
<구체적 수정. 코드/설정 diff 또는 변경한 값. "조심한다" 같은 서술 금지.>

```diff
```

## 재발 방지
<추가한 테스트 / 런타임 가드 / 문서. 없으면 "없음"이라고 명시할 것.
 "없음"이 반복되는 실패 유형은 우선적으로 가드를 만들 신호다.>

## 참고
<관련 항목 [[other-id]], 문헌, run 링크>
