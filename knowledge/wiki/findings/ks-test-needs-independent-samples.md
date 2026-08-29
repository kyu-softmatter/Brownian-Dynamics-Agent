---
id: ks-test-needs-independent-samples
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
question: "분포 적합도 검정(KS)을 시뮬레이션 궤적에 걸 때 무엇을 주의해야 하는가?"
answer: "상관된 표본에 KS 를 걸면 항상 기각한다. 표본이 많을수록 더 확실히 기각한다 — 직관이 정확히 반대로 작동한다."
status: measured
cause_class: analysis
stage: S7
runs: [2026-07-28_trap-2d-5um_2dfb9d]
cites:
  - knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
affects_docs: [master_plan.md, simbot/analysis/trap.py]
reproduced: yes
---

# 상관표본 KS 검정 — 거짓 기각

## 증상
조화 트랩 런 16개 전부에서 위치 분포의 KS p값이 **`0.0000`** 이었다.
예측 `P7`("위치는 Gaussian, `p > 0.05`")이 **FAIL** 로 보였다.

그런데 같은 데이터의 **첨도는 `2.94 – 3.00`** 으로 Gaussian(`3.00`)과 잘 맞았다.
두 진단이 정반대를 말한 것이 단서였다.

## 진단 경로
1. **물리를 의심** — 균일 노이즈(HOOMD)가 정말 Gaussian 이 안 되나?
   → 기각. 완화시간이 `200` 스텝이라 CLT 가 충분히 작동한다. 첨도가 그것을 뒷받침한다.
2. **폭 차이를 의심** — EM 편향으로 `⟨x²⟩ = 1.0058` 인데 `N(0,1)` 과 대조했다.
   → 부분 원인이지만 크기가 안 맞는다. 폭 0.29 % 차이가 만드는 KS 통계량은
   `φ(1) × 0.0029 ≈ 0.0007` 인데 관측값은 `0.0106` 이었다.
3. **표본 독립성을 의심** ← **적중.**
   마지막 20 프레임을 썼는데 프레임 간격이 `16` 스텝 = `0.08 τ_trap` 이었다.
   **20 프레임 전체가 `1.6 τ_trap` 구간 안**에 있었다 — 사실상 한 개의 배위(configuration).
   50,000 점을 iid 로 취급했으나 실효 독립표본은 `~1000`(입자 수) 수준이었다.

## 근본 원인
`분류: analysis`

KS 검정의 귀무분포는 **iid 표본**을 전제한다. 임계값이 `1.36/√n` 로 `n` 에 의존하므로,
상관된 표본을 `n` 개로 세면 임계값이 실제보다 `√(n/n_eff)` 배 작아진다.

> **`n` 을 키우면 기각이 더 확실해진다.** "표본이 많아서 믿을 수 있다"는 직관이
> **정확히 반대로** 작동한다. 이것이 이 실패가 위험한 이유다.

## 처방
```python
# 이전 — 마지막 20 프레임을 무조건 사용
frames = traj[-20:]
x = frames.reshape(-1)
ks = kstest(x, norm(0, 1).cdf)          # ← 상관 무시 + 폭까지 검정

# 이후 — 독립 프레임만, 형태만
frames_per_tau = 1.0 / (frame_interval_steps * dt_star)
step = ceil(2.0 * frames_per_tau)        # 2 tau 이상 떨어진 프레임만
x = traj[::step].reshape(-1)
x = (x - x.mean()) / x.std(ddof=1)       # 폭은 P1 의 일. 여기선 형태만
ks = kstest(x, norm(0, 1).cdf)
```
결과: **`p = 0.0000 → 0.29 – 0.53`**

두 번째 수정도 중요하다. 규격화하지 않으면 KS 가 **폭 차이(EM 편향)** 까지 잡아
"형태가 Gaussian 인가"가 아니라 "폭이 정확히 1 인가"를 검정하게 된다.
그건 `P1`(등분배)의 일이고, 한 검정에 두 질문을 섞으면 어느 쪽이 실패했는지 알 수 없다.

## 재발 방지
`simbot.analysis.trap.check_position_distribution` 이 **`dt_star` 와
`frame_interval_steps` 를 필수 키워드 인자로 요구한다.**
독립성을 계산하지 않고는 호출 자체가 불가능하다 — 기본값을 주면 다시 실수한다.

반환값에 `n_independent_frames` 와 `n_effective_samples` 를 실어 보고서에 노출한다.

## 일반화 — 이 함정이 적용되는 다른 곳
| 검정 | 상관표본에서 무슨 일이 생기나 |
|---|---|
| KS · Anderson–Darling · χ² 적합도 | **거짓 기각.** `n` 이 클수록 심함 |
| `t` 검정 · 평균 비교 | 오차막대 과소평가 → 거짓 유의 |
| 분산 추정 `sqrt(2/n)` | 과소평가 → **거짓 정밀도** |
| 상관계수 유의성 | 거짓 유의 |

> **규칙: 시계열에서 뽑은 표본에 통계검정을 걸기 전, 상관시간으로 나눈 실효표본 수를
> 먼저 계산한다.** 그 수를 보고서에 적는다.

관련: 실효표본 수 자체를 실측한 사례는 [[hoomd-brownian-scheme-and-noise]] §3
(입자간 독립성 `996/1000`). 거기서는 **평균을 빼서 항등식을 측정**하는 반대 방향의
실수를 했다 — 통계 진단은 양쪽으로 틀릴 수 있다.
