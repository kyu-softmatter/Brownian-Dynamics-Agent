---
type: finding
author: agent
drafted: 2026-07-28
confirmed_by:
question: "HOOMD 7.1 md.methods.Brownian 은 정확히 어떤 스킴인가? 노이즈는? N개 비상호작용 입자는 독립표본인가? BD 에서 온도를 어떻게 검증하나?"
answer: "Euler-Maruyama(1차, 편향 dt*/2) 확인. 노이즈는 **균일분포**(첨도 1.81). 입자간 독립 확인(실효표본 996/1000). 운동온도는 진단 불가, 배위온도를 써야 한다."
status: measured
cause_class: numerical
cites:
  - knowledge/wiki/systems/passive-sphere--harmonic-trap.md
  - runs/2026-07-28_trap-2d-5um_2dfb9d/02_prediction.md
affects_docs: [master_plan.md, CLAUDE.md, knowledge/wiki/systems/passive-sphere--harmonic-trap.md]
engine: hoomd 7.1.0 (cpu_py312)
reproduced: yes
---

# HOOMD `Brownian` 의 적분 스킴·노이즈·온도 진단

**측정** 2026-07-28 · 조화 트랩 축약 단위(`ℓ_trap=1, kT=1, τ_trap=1`) · 2D · `simbot.forces.HarmonicTrap`

---

## 1. 스킴은 Euler–Maruyama다 — 확인

docstring: *"uses the numerical integration method from I. Snook 2007 §6.2.5 …
This numerical method has errors in `O(δt)`"*

EM의 정상분산은 **정확히** `1/(1−dt*/2)`다. dt 래더로 판정:

| `dt*` | EM 예측 | exact 예측 | 측정 | ±stat | 판정 |
|---|---|---|---|---|---|
| `2.0e-2` | **1.01010** | 1.00000 | **1.01041** | 0.00264 | **EM `0.1σ` · exact `3.9σ`** ✅ 결정적 |
| `1.0e-2` | **1.00503** | 1.00000 | **1.00512** | 0.00284 | **EM `0.0σ`** ✅ |
| `5.0e-3` | 1.00251 | 1.00000 | 1.00100 | 0.00310 | ⚠️ **INCONCLUSIVE** |
| `2.5e-3` | 1.00125 | 1.00000 | 0.99722 | 0.00301 | ⚠️ **INCONCLUSIVE** |

> **결론: EM 확정.** 판별력이 있는 두 지점(`dt* ≥ 1e-2`)에서 예측과 `0.1σ` 이내로 일치하고,
> exact 스킴 가설은 `3.9σ`로 기각된다.
>
> **작은 `dt*`에서 "불명확"인 것은 EM에 대한 반증이 아니다** — 편향(0.25 %)이
> 통계오차(0.31 %)보다 작아서 판별이 불가능한 것이다. 이것이 `master_plan` §S7의
> `INCONCLUSIVE` 판정이 필요한 정확한 이유다. 판별하려면 표본을 **36배** 늘려야 한다
> (`0.25 %`를 `3σ`로 보려면 SE `≤ 0.08 %` 필요).

**함의:** `dt*`를 줄이면 편향은 줄지만 **검증 가능성도 함께 줄어든다.**
적분기 검증은 **일부러 큰 `dt*`에서** 해야 한다. 프로덕션 `dt*`에서는 검증할 수 없다.

## 2. ★ 노이즈가 균일분포다 — Gaussian이 아니다

docstring: *"with the exception that `F_R` is drawn from a **uniform** random number distribution"*

자유 BD(힘=0) 한 스텝 변위, 표본 20 000개:

| 측정 | 값 | Gaussian | 균일분포 |
|---|---|---|---|
| 표준편차 | `4.4622e-2` | `√(2D·dt) = 4.4721e-2` | (분산은 요동-소산으로 맞춰짐) |
| **첨도** | **`1.8125`** | 3.00 | **1.80** ✅ |
| `max/σ` | **`1.7358`** | 무한 | **`√3 = 1.7321`** ✅ |

**균일분포 확정.**

### 무엇이 영향받지 않는가
`⟨x²⟩`, MSD, `D`, 모든 2차 모멘트 — **선형 SDE는 노이즈의 2차 모멘트만 관여**하므로
분포 형태와 무관하다. 위 표에서 표준편차가 예측과 0.2 % 이내로 맞는 것이 그 확인이다.

### 무엇이 영향받는가 — ⚠ 주의
| 항목 | 영향 |
|---|---|
| **위치 분포의 Gaussian성** | 근거가 바뀐다. "Boltzmann이라서"가 아니라 **"균일 증분 `τ/dt`개가 CLT로 합쳐져서"** 다. `dt* = 5e-3`이면 200스텝이 합쳐지므로 통과하지만, `dt*`가 크면 깨진다 |
| **스텝당 변위 분포** | **Gaussian이 아니다.** 런타임 가드가 변위 분포를 Gaussian으로 가정하면 안 된다. `max/σ = √3`이 상한이다 |
| **꼬리 사건 (rare events)** | **꼬리가 유한하다.** 큰 킥이 아예 발생하지 않는다 |
| **장벽 넘기 · 탈출률** | ⚠️ **속도가 틀린다.** Kramers 탈출은 꼬리 사건이 지배하는데 그 꼬리가 잘려 있다. 트랩 탈출·응집 해리 문제에서는 이 스킴을 신뢰할 수 없다 |

> **가장 중요한 실무 결론: 활성화 장벽 넘기(`βΔU ≳ 5`)가 관여하는 문제에
> `md.methods.Brownian`을 쓰면 안 된다.** `dt`를 아주 작게 하면 CLT로 회복되지만
> (스텝 수가 많아져 합이 Gaussian화) 필요한 `dt`가 얼마나 작아야 하는지는 미측정이다.
> → `questions/` 항목으로 남길 것.

## 3. N개 비상호작용 입자는 독립표본이다 — 확인

정확한 검사: **표본평균의 분산이 개별분산/N 인가.**
```
독립이면   Var_t(mean_i d_i) = <Var_i(d_i)> / N
양의 상관   좌변이 커진다
```
스텝 3000개 × 입자 1000개, 자유 BD의 한 스텝 x변위:

| 양 | 값 |
|---|---|
| `⟨Var_i(d_i)⟩` | `9.9961e-3` |
| `Var_t(mean_i d_i)` | `1.00370e-5` |
| 독립 예측 | `9.9961e-6` |
| **비 `R`** | **`1.0041 ± 0.0259`** → `1`에서 `0.16σ` ✅ |
| **실효 표본 수** | **`N/R = 996` / 1000** |

실용 검증 (시드 40개의 `⟨x²⟩` 산포 vs `√(2/(N·d))` 예측):
`실측/예측 = 0.872 ± 0.099` → `1.30σ` ✅ 통과

> **결론: HOOMD의 counter-based RNG는 입자별로 독립 스트림을 준다.**
> `N`을 앙상블로 쓰는 것이 정당하다 → S2의 "1 % 정밀도 1.3초" 근거가 확보됐다.

### ⚠️ 진단 경로 — 첫 시도는 무의미한 테스트였다
처음엔 "입자간 변위의 교차상관"을 재려 했고, 변위에서 평균을 뺀 뒤 계산했다:
```python
d = d - d.mean()
cross = (d.sum()**2 - (d**2).sum()) / (n*(n-1))
```
**평균을 빼면 `sum(d) = 0`이 항등적으로 성립**하므로
`cross = -Σd²/(n(n-1))`, 따라서 `cross/auto = -1/(n-1)` 이 **언제나 정확히** 나온다.

측정이 아니라 **산술 항등식**이었다. 결과가 `-1.001e-3 = -1/999`로 나왔고,
**200스텝 반복의 표준편차가 `6.7e-20`**(= 부동소수점 잡음)이었던 것이 결정적 증거다.

> **재발 방지 규칙: 상관을 측정할 때 표준편차가 비정상적으로 작으면 (통계 요동이
> 보이지 않으면) 측정이 아니라 항등식을 계산한 것이다.**
> 통계량은 **반드시 요동해야 한다.** 요동하지 않는 "측정값"은 버그다.
> 시드 8개로 잰 `실측/예측 = 0.64`가 40개에서 `0.872`로 움직인 것이 정상적인 통계 거동이다.

## 4. 온도 진단 — 운동온도는 못 쓴다. 배위온도를 쓴다

docstring: *"At each time step, `Brownian` **draws** a new velocity distribution
consistent with the current set temperature so that `ThermodynamicQuantities`
will report appropriate temperatures"*

> ⚠️ **속도는 적분되지 않고 매 스텝 목표 분포에서 뽑힌다.**
> `kinetic_temperature`는 입력값을 되풀어 말할 뿐이고, **계통적으로 이탈할 수 없다.**
> 측정값 `0.991818`은 유한표본 요동일 뿐이다.
> ⇒ `master_plan` §S5의 가드 "온도가 목표 `kT`에서 이탈"은 **작동할 수 없다. 삭제.**

### 대체: 배위 온도 (configurational temperature)
```
kT_conf = <|grad U|^2> / <laplacian U>          (위치와 힘만 사용)
조화 트랩: = k^2<r^2>/(d k) = k<r^2>/d = kT
```
| 측정 | 값 |
|---|---|
| `⟨\|∇U\|²⟩` | `2.00763` |
| `⟨∇²U⟩ = d·k` | `2.0` |
| **`kT_conf`** | **`1.00382 ± 0.00480`** (입력 `1.0`) → `0.80σ` ✅ |

> **주의: 순수 조화 트랩에서 `kT_conf`는 `⟨x²⟩ = kT/k` 검사와 대수적으로 동일하므로
> 새 정보가 아니다.** 독립적인 검사가 되는 것은 **쌍 상호작용(WCA 등)이 있을 때**다 —
> 그때 `∇U`가 트랩만이 아니라 이웃 입자들로부터도 오기 때문이다.
>
> ⇒ **`kT_conf`를 `simbot.guards`의 표준 가드로 채택한다.** 힘은 이미 계산되어 있으므로
> 추가 비용이 거의 없다.

**`kT` 자체는 여전히 1급 입력이다.** `U/kT`가 Boltzmann 가중치를 정하고 모든 무차원수에
들어간다 (`k* = kσ²/kT`, `T* = kT/ε`, `Pe = Fσ/kT`, `ℓ_trap = √(kT/k)`).
쓸 수 없는 것은 **운동에너지 되읽기**뿐이고, 그 자리를 배위온도가 대신한다.

## 5. 부수 확인
- `moment_inertia = 0` 이면 회전 자유도가 적분되지 않는다 (docstring: *"About axes where `I^i > 0`"*).
  → ABP에서 `force.Active` + `update.ActiveRotationalDiffusion` 과 함께 쓸 때 필수. 함정 회피 확인.
- `snap.configuration.dimensions` 는 **setter가 없다.** `box`의 `Lz = 0` 으로 2D를 지정한다.
- `md.force.Custom` 은 `self._state.cpu_local_snapshot` + `self.cpu_local_force_arrays` 로
  구현한다. `virial` 순서는 `xx, xy, xz, yy, yz, zz`.

## 6. 남은 질문
- **균일 노이즈로 장벽 넘기를 신뢰하려면 `dt`가 얼마나 작아야 하나** (`βΔU` 의존)
- HOOMD가 균일분포를 쓰는 이유 (속도? GPU 친화성?) — 성능 이득이 정확도 손실을 정당화하는가
- `dt* ≤ 5e-3` 에서 EM 편향을 실제로 확인하려면 36배 표본 — 할 가치가 있나
  (조화 트랩은 6초면 되므로 사실 저렴하다)

## 참고
관련: [[dt-gate-should-be-displacement-based]], [[local-cpu-parallelism]],
[[passive-sphere--harmonic-trap]]
