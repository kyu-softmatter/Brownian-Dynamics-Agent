---
name: bd-hoomd
description: |
  HOOMD-blue v7.1 API 매핑, 실행 검증된 코드 스니펫, 실측으로 발견한 함정 18개.
  Brownian dynamics 시뮬레이션 코드를 쓰기 전에 반드시 읽어라 — 함정 중 여럿은
  에러 없이 조용히 틀린 결과를 낸다(+1856% 오차 / 사슬이 거의 곧으면 angle.Harmonic 의
  **힘만** 최대 96% 틀리는데 에너지는 정확해서 에너지 검증으로는 안 잡힌다 /
  BoxResize 압축이 좁은 우물의 결합을 비가역 붕괴시킨다). HOOMD/hoomd.md 코드 작성,
  트랩·액티브·사슬·비구형 입자 구현, **3D·박스 압축(겔화)**, 궤적/힘 저장, 재시작,
  런타임 감시에 해당.
---

# HOOMD-blue v7.1 — 검증된 사용법

> 이 문서의 모든 스니펫과 수치는 **설치본에서 실제로 실행해 확인한 것**입니다.
> 근거: `scratch/` 스크립트, [docs/hoomd_capabilities.md](../../../docs/hoomd_capabilities.md)
> 추측을 추가하지 마세요. 새 사실은 실행으로 확인한 뒤 여기에 기록합니다.

## 환경

```bash
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
```
hoomd **7.1.0** · `gpu=False` · **`mpi=False`** · macOS arm64 · double precision

**MPI가 없으므로 단일 런 = 단일 코어입니다.** 병렬성은 "여러 런 동시 실행"에서만 나옵니다.

---

# ⚠️ 함정 — 코드 쓰기 전에 확인

**★ 표시는 에러 없이 조용히 틀리는 것.** 가장 위험합니다.

### ★ 1. 외부 힘 + 주기경계 → 최소 이미지 필수

고정 앵커를 향한 힘에서 `d = pos - anchor`를 그대로 쓰면, 입자가 박스를 넘어 wrap되는 순간
거리가 L만큼 점프해 **반대 방향의 거대한 힘**을 받습니다. **터지지 않습니다.**
트랩이 강하면 맞는 값이 나오고 약할수록 틀립니다:

| k | 최소이미지 없음 | 있음 |
|---|---|---|
| 2 | **+1856%** ✗ | +0.38% ✓ |
| 5 | +344% ✗ | +0.56% ✓ |
| 10 | +0.16% ✓ | −0.02% ✓ |

강한 조건만 테스트하면 통과했다고 착각합니다. **약한 조건으로 검증하세요.**

### ★ 2. Brownian은 O(δt) 오차

문서가 명시합니다. Langevin/ConstantVolume보다 훨씬 작은 `dt`가 필요합니다.
σ=kT=γ=1 기준 `dt ≈ 1e-4`부터 시작하고, 힘이 강하면(Pe 큼, ε 큼, k 큼) 더 줄입니다.
**상한을 통과해도 정확도는 보장되지 않습니다** — 논문급 결과엔 `dt` 절반 수렴 확인이 별도로 필요.

### ★ 3. ABP 회전은 적분기가 아니라 updater로

```python
integrator.integrate_rotational_dof = False        # ← 필수
sim.operations.updaters.append(
    active.create_diffusion_updater(trigger, rotational_diffusion=D_r))
```
`True`로 두면 관성 회전이 섞여 ABP가 아니게 됩니다. 조용히 다른 물리가 됩니다.

### ★ 4. WCA 전용 클래스는 없다

```python
md.pair.LJ(nlist=cell, default_r_cut=2**(1/6)*sigma, mode='shift')   # = WCA
```
`ForceShiftedLJ`는 WCA가 **아닙니다** (다른 퍼텐셜). 이름에 속지 마세요.

### ★ 5. BD는 과감쇠 — 속도에 물리적 의미가 없다

`thermalize_particle_momenta()` 불필요. **속도 기반 MSD/VACF 금지.** 위치 차분만 씁니다.

### 6. `r_cut < L/2`

아니면 minimum image convention 위반. 작은 박스에서 자주 걸립니다.

### 7. 2D 최소이미지의 z성분 NaN

비주기 축의 주기를 `np.inf`로 두면 `inf * round(0/inf) = inf*0 = nan`.
힘 배열에 NaN이 들어가 런타임 가드가 오탐합니다. **주기 축만 마스크**하세요 (아래 스니펫).

### 8. `write.Burst`는 새 파일에 `write_at_start=True` 필요

없으면 `RuntimeError: Must set write_at_start to write to a new file.`
결과 파일에 초기 프레임 1개가 추가로 들어갑니다 (버퍼 10 → 파일 11).

### 9. 2D는 `Lz=0` + `dimensions=2`

박스 z를 0으로 두고 `configuration.dimensions = 2`. 둘 다 필요합니다.

### ★ 10. `pair.Table`의 격자는 `endpoint=False`

문서: *"implicitly defined r values are those returned by
`numpy.linspace(r_min, r_cut, len(U), endpoint=False)`"*.
`endpoint=True`로 만들면 표가 통째로 어긋나 **힘이 조용히 틀립니다.**
실측 (`scratch/verify_pair_table.py`, `U=A/r³`, r_min=0.5, r_cut=3.0, 200 bins):

| sep | `endpoint=False` | `endpoint=True` |
|---|---|---|
| 0.70 | +0.000% | −0.572% |
| 1.50 | +0.000% | −1.329% |
| 2.90 | +0.000% | **−1.646%** |

`endpoint=False`는 **0.000%로 정확**합니다. 컷오프에 가까울수록 오차가 커집니다.

### ★ 11. `pair.Table`은 `r < r_min`에서 힘·에너지가 **0**

문서에 `F=0, U=0 for r < r_min`이라고 적혀 있는데, **발산 퍼텐셜에서는 이게 함정**입니다.
반발이 사라지므로 입자가 `r_min`을 한 번 뚫으면 그대로 겹쳐 있습니다. 터지지 않습니다.

실측 (`A/r³`, r_min=0.5): sep=0.51 → F=443.9 (정상) / sep=0.49 → **F=0.0, U=0.0**.

→ 방어 두 겹: **배제부피 코어를 별도 힘으로** 두어 `r_min` 도달을 막고
(WCA + Table은 정확히 합해집니다 — 실측 오차 0.000%), 런 중 **최소 이웃거리를 감시**합니다.

### ★ 12. `seed`는 16비트로 잘린다 — 다른 seed가 같은 런이 될 수 있다

`Simulation(seed=...)`에 65535를 넘는 값을 주면 경고와 함께 `seed % 65536`으로 잘립니다.
경고는 뜨지만 **런은 그대로 진행**되므로, 반복 런의 seed를 65536 간격으로 주면
전부 **완전히 동일한 궤적**이 나옵니다 (독립 표본인 줄 알고 평균을 냄).

실측 (자유 BD, N=40, 2000 스텝, 좌표 최대차):

| seed 쌍 | 차이 |
|---|---|
| `20260803` vs `10179` (= mod 65536) | **0.000e+00** |
| `20260803` vs `20260803+65536` | **0.000e+00** |
| `20260803` vs `10180` | 9.97 (정상적으로 다름) |

→ 반복 런은 **작은 연속 정수**(1,2,3,…)를 쓰세요.

---

### ★ 13. `active_force = 0` 이면 `ActiveRotationalDiffusion` 이 **동작하지 않는다**

방향이 전혀 안 풀립니다 — `⟨n(0)·n(t)⟩` 감쇠율 **Λ = 0.0000** (4개 조합 전부 실측).

```python
act.active_force["A"] = (0.0, 0.0, 0.0)     # ← 회전확산도 같이 꺼진다
```

**활성력을 끄고 회전 통계만 따로 보려는 시도가 막힙니다.** 방향 확산만 검증하려면
작은 활성력을 켜 두고 director 를 재세요 (`scratch/standalone_abp_diffusion.py`).

### ★ 14. `rotational_diffusion` 은 **director 감쇠율 Λ 그 자체** — 3D에서 2배 어긋난다

HOOMD 에 `rotational_diffusion=D_r` 를 주면 실측 결과:

| dim | `Λ/D_r` | `Λ/[(d−1)D_r]` (표준 이론) |
|---|---|---|
| 2 | **1.00** | 1.00 ✓ |
| 3 | **1.00** | **0.50** ✗ |

즉 `⟨n(0)·n(t)⟩ = exp(−D_r·t)` 가 **2D·3D 모두** 성립합니다.
표준 구형 회전확산은 `exp(−(d−1)D_r^phys·t)` 이므로, **3D에서 물리적 `D_r^phys` 를
재현하려면 `rotational_diffusion = 2·D_r^phys` 를 넣어야 합니다.** 2D는 그대로.

따라오는 결과 — 자유 ABP 유효확산:
```
D_eff = D_t + v₀² / (d · Λ)          Λ = HOOMD 의 rotational_diffusion
```
실측 오차 1.5% 이내 (2D·3D, v₀ 2종). 흔히 쓰는 `v₀²/[2(d−1)D_r]` 는
**2D에서만 우연히 일치**하고 3D에서 +29~31% 틀립니다.

> 우리 케이스는 전부 2D라 지금은 영향이 없습니다. **3D 액티브를 하는 순간 물립니다.**

---

### ★★ 15. `angle.Harmonic` 은 **거의 곧은** 사슬에서 에너지는 맞고 **힘이 틀린다**

`md.angle.Harmonic` 은 토크를 좌표로 옮길 때 `1/sin θ` 를 쓰고, **`sin θ` 를 아래로 클램프**한다.
실측 확정: **SMALL = 1.414217e-03 = √2×10⁻³** (표준편차 1.4e-7, 6개 진폭에서).

```
sin θ < SMALL  →  힘이 정확히  sinθ/SMALL  배로 축소된다
```

`t0 = π`(곧은 사슬이 평형)이면 **평형 자체가 sin θ = 0 인 특이점**이고, 뻣뻣한 사슬은
항상 거의 곧다. 즉 이 함정은 **강체 필라멘트에서 가장 심하게** 작동한다.
`sin θ ≈ |θ−π|` 이므로 힘이 `∝ κ(θ−π)·(θ−π)/SMALL = κ(θ−π)²/SMALL` — **선형이 아니라 2차**다.
사슬이 실제보다 훨씬 **무르고 비선형**이 되고, 곧을수록 심해진다.

실측 (κ_θ = 1.39e6 kT, n=5, 매끄러운 굽힘 모드):

| max\|θ−π\| | 에너지 오차 | 힘비 (정답 대비) |
|---|---|---|
| 1.757e-02 | 0.0000% | **1.000000** ✓ |
| 1.757e-03 | 0.0000% | **1.000000** ✓ |
| 8.787e-04 | 0.0000% | 0.621320 ✗ |
| 2.929e-04 | 0.0000% | 0.207107 ✗ |
| 5.858e-05 | 0.0000% | **0.041421** ✗ |

**★ 에너지는 전 구간에서 정확히 0.0000% 다.** 그래서 **에너지로 검증하면 통과하고 힘은
틀린 채로 남는다.** 실제로 이 프로젝트에서 `scratch/verify_angle_matrix.py` 가
`k = 2U/δ²` 로 재서 0.55% 로 **통과했고**, 힘 버그를 놓쳤다 (입자 순서 문제를 피하려고
일부러 에너지를 골랐던 것이 함정을 가렸다).

> **교훈: 퍼텐셜을 에너지로만 검증하면 힘을 검증한 것이 아니다.**
> 새 퍼텐셜은 **힘을 에너지의 수치 기울기와 직접 대조**해야 한다.

**우회책이 마땅치 않다** (둘 다 실측/해석으로 배제):
- `angle.Table` 에 같은 `U`·`tau` 를 줘도 어긋난다 (힘비 3.27 — 표 해상도 Δθ=3.1e-3 가
  `|θ−π|~3e-4` 를 못 담는다). 폭을 키워도 같은 `1/sin θ` 경로를 탄다.
- `angle.CosineSquared`: `U = ½k(cosθ − cosθ₀)²`. `θ₀=π` 에서 `cosθ+1 ≈ (π−θ)²/2` 이므로
  `U ≈ k(π−θ)⁴/8` — **4차**다. 조화 굽힘 강성을 낼 수 없다.

→ **대응**: ① 계를 `min|θ−π| > 1.41e-3` 안에 두거나(전 각도가! 최대만 보면 안 된다 —
사슬 양끝 각도가 중앙보다 한 자릿수 작다) ② `force.Custom` 으로 굽힘을 직접 구현하거나
(정확하지만 26배 느림 — 함정 16) ③ **선형 영역이면 MD 를 쓰지 말고 해석적으로 푼다.**

재현: `scratch/verify_angle_force_small_theta.py` · 추적 경로:
`scratch/diagnose_chain_bend_28pct.py` (dt·비선형·x자유도·에너지를 차례로 배제)

### ★ 16. `force.Custom` 은 매 스텝 파이썬 — 작은 계에서 26배 느리다

`set_forces` 가 매 스텝 호출된다. 실측 (n=25 사슬, 2D, CPU):

| 구성 | steps/s |
|---|---|
| `bond` + `angle` 만 (컴파일) | 61,264 |
| **+ `force.Custom` 트랩** | **2,339** (26배 느림) |
| **+ 유령입자 `bond.Harmonic(r0=0)` 트랩** | **55,551** |

`bond.Harmonic(r0=0)` 은 `U = ½k r²` 로 **조화 트랩과 함수가 정확히 같다.** 유령입자를
적분기 `filter` 에서 빼면 움직이지 않으므로 고정 트랩이 된다 → 전부 컴파일 경로.
시간의존 구동은 유령을 `CustomUpdater` 로 옮긴다 (매 스텝일 필요가 없다 — 아래 함정 17).

→ **"비싸서 못 돌린다" 판정 전에 구현 비용인지 물리 비용인지 재라.**
재현: `scratch/bench_chain_bend.py`

### ★ 17. 앵커를 U 스텝마다 옮기면 구동이 **영차 유지(ZOH)** 가 된다

`CustomUpdater` 로 트랩 중심을 `U` 스텝마다 갱신하면 구동 사인파가 계단이 되어
기본파가 `sinc(ωΔt/2)` 배로 줄고 위상이 `ωΔt/2` 늦는다 (`Δt = U·dt`).
실측 `|ŷ_c|/a = 0.98999`, 위상 `−0.2522 rad` / ZOH 예측 `0.99040`, `−0.2404 rad`.

→ **응답함수 추정량에 공칭 진폭을 쓰지 마라.** 앵커(유령) 위치를 **같이 재서** 측정된
위상자를 쓰면 ZOH 감쇠가 분자·분모에서 정확히 상쇠된다. 공칭을 쓴 `K′` 은 De=10 에서
**−6559 (부호까지 틀림, 오차 236%)**, 측정 위상자로는 5863 (21%).
재현: `scratch/verify_chain_bend_gates.py --gate lockin`

### ★ 18. `update.BoxResize` 는 좌표를 **아핀 스케일**한다 — 좁은 우물의 결합이 부서진다

문서대로 `r → s·r` 로 스케일합니다 (**실측 오차 8.9e-16**). 즉 압축은
**이미 결합된 쌍의 결합길이도 같이 줄입니다.** DLVO 2차극소처럼 우물이 좁으면
한 트리거의 변형이 우물 폭을 넘는 순간 쌍이 장벽 안쪽으로 밀려
**1차극소(접촉)로 비가역 붕괴**합니다. 터지지 않습니다 — 결합 종류만 조용히 바뀝니다.

```
ε_crit = (h_min − h_barrier) / ℓ          # 트리거당 허용 선형 변형
```

실측 (DLVO: h_min*=0.007593, 장벽*=0.000508, ℓ*=1.007593 → 예측 **0.703%**):

| ε/트리거 | 최종 h* | 결과 |
|---|---|---|
| 0.40% | 0.007591 | 유지 ✓ |
| 0.80% | **−0.042257** | 붕괴 (겹침 4.2%) |
| 2.00% | −0.042519 | 붕괴 |

→ 압축 전에 **ε_crit 을 원장에서 계산하고 단계 수를 그것으로 나눠 정하라.**
총 변형 41.5%(φ 0.02→0.10)를 0.4% 로 쪼개면 **134단계**다.
⚠️ **단계 수는 물리가 정하고, 단계당 이완시간은 목적이 정한다** — 준정적(≫τ_B)으로
하면 비용이 1800시간이고, τ_bond 규모만 이완하면 42초인데 **구조가 다르다**
(확산이 개입한 겔 vs 아핀 압축된 초기 배치). 어느 쪽인지 반드시 명시하라.
또 압축이 `r_cut < L/2` 를 깰 수 있으니 **최종 박스에서 다시 확인**하라 (함정 6).

재현: `scratch/verify_3d_boxresize.py` · 상세: docs/hoomd_capabilities.md §5.5

### ★ 19. `md.pair.friction` 의 지렛대 `R` 은 `sigma` 가 아니라 `particles.diameter/2`

접선력·토크를 내는 `FrictionLJ*` 3종은 `params` 에 `sigma` 를 받지만, 접촉점 지렛대는
**입자 속성 `diameter`** 에서 가져옵니다. `diameter` 기본값은 **1.0** 이라 `sigma≠1` 이면
조용히 어긋납니다. 실측 (`|τ|/F_tan` 로 R 을 역산):

| `particles.diameter` | 측정 `R` | |
|---|---|---|
| 기본(1.0), `sigma=0.8909` | **0.5000000000** | σ/2 보다 **1.1225배 과대** ✗ |
| `= sigma` 로 명시 | 0.4454493591 | `= σ/2` 정확 ✓ |

에러가 없고 힘의 크기는 그럴듯하게 나옵니다. **토크와 무슬립 구름 조건만 틀립니다.**
→ 이 힘을 쓸 때는 `frame.particles.diameter` 를 반드시 같이 설정하세요.

### ★★ 20. `Brownian` 의 `velocity` 는 0 이 아니라 **무상관 열잡음**이다

함정 5("BD 는 과감쇠 — 속도에 물리적 의미가 없다")의 구체적 귀결인데, 위험한 방향이
반대입니다. `methods.Brownian` 은 속도를 0 으로 두지 않고 **매 스텝 Maxwell–Boltzmann
에서 새로 뽑아 채웁니다** (실측 `⟨v²⟩ = 2.771` vs `3kT/m = 3.000`, 비 0.924, N=80).
그런데 그 값은 실제 변위율과 **무관**합니다 — 자유 BD 1스텝 실측:

```
저장된 v      = [-0.634, -0.306, -0.625]
실제 Δx/dt    = [14.36, -95.99, -21.53]      비 = [-0.044, +0.003, +0.029]
```
과감쇠 변위는 `√(2Dδt) ∝ δt^½` 라 `Δx/dt` 는 `δt→0` 에서 발산하고, `v` 는 그것과 별개로
추첨된 값입니다.

⟹ **속도 의존 페어힘(`md.pair.friction`, DPD 계열)을 `methods.Brownian` 과 결합하면
마찰이 보는 것은 실제 상대운동이 아니라 열잡음입니다.** 힘이 0 이 아니라 **쓰레기값**이
나오므로 "결과가 나왔다"로는 절대 못 잡습니다.

재현: `scratch/verify_pair_friction.py` (23/23) · 3종의 구조적 성질(산일성·구름 면제·
컷오프 밖 0)은 [docs/hoomd_capabilities.md §마찰 접촉](../../../docs/hoomd_capabilities.md)

### 3D 로 넘어갈 때 (실측 확인) ✓

`network` 가 이 프로젝트 최초의 3D 케이스입니다. 확인된 것:

```python
f.configuration.box = [L, L, L, 0, 0, 0]      # 3D 는 Lz=L (2D 는 Lz=0 — 함정 9)
f.configuration.dimensions = 3
```
자유 BD 3D 가 `⟨r²⟩ = 6·D·t` 를 재현합니다 (**−1.40%, SEM 2.81%**), 세 축이 각각
`⟨Δx²⟩=2Dt` (x/y/z = 1.023/0.938/0.997 — 축 누락 없음).
⚠️ 언랩 좌표는 `snapshot.particles.image × L` 을 더해서 만드세요.
⚠️ 3D 액티브는 함정 14(`rotational_diffusion` 2배)에 물립니다 — 이 케이스는 액티브가 없습니다.

---

# 하드 제약 (우회 불가 — 실측 확인)

### 병진 마찰은 스칼라뿐. 이방성을 만들 수 없다

`Brownian(default_gamma=float)` — 타입별 **스칼라**. 회전만 `default_gamma_r=(x,y,z)` 텐서.

두 경로 모두 실측했고 **둘 다 등방**입니다:

| 구성 | γ⊥/γ∥ |
|---|---|
| 강체 `constrain.Rigid` (구 5개) | **1.000000** |
| 자유배수 비드 막대 (본드+각도, 구 5개) | **1.000000** |
| 이론 (slender body, HI 있음) | → 2 |

강체 구성 입자에 `gamma['A']=10`을 줘도 병진 속도가 전혀 안 변합니다.
**강체의 병진 항력 = 중심 입자의 스칼라 gamma 하나.**

**왜**: 자유배수에서 각 비드는 독립적으로 Stokes 항력을 받아 총 항력이 `N·γ_bead`로 방향 무관.
실제 막대의 `γ⊥/γ∥→2`는 **유체역학적 상호작용(HI)** 효과이고 BD엔 HI가 없습니다.
**HOOMD 한계가 아니라 BD 모델 자체의 성질**입니다. 기하로 우회하려 하지 마세요.

→ 대응: 등방 평균 `γ̄`를 쓰고 모듈이 한계를 선언. 상세는 `docs/hoomd_capabilities.md` §5.1–5.4.

### 회전 마찰 텐서는 정상 동작 ✓

| `gamma_r` | τ∥x 회전각 | τ∥z 회전각 | 비 |
|---|---|---|---|
| (1,1,1) | 2.00000 | 2.00000 | 1.0000 |
| (1,1,5) | 2.00000 | 0.40000 | **0.2000** |

MSAD·회전 동역학은 정확히 재현됩니다.

### 조화 트랩은 내장에 없다

`md.external.field`에는 `Electric`/`Magnetic`/`Periodic`뿐.
`hoomd.hpmc.external.Harmonic`은 몬테카를로 전용이라 MD/BD에서 못 씁니다.
→ `md.force.Custom`으로 직접 구현 (아래 스니펫).

---

# 검증된 스니펫

## 2D 프레임 + BD + WCA (기본형)

```python
import itertools, math
import numpy as np, gsd.hoomd, hoomd, hoomd.md as md

def frame_2d(n_side=40, phi=0.5, sigma=1.0):
    """면적분율 phi인 2D 정사각 격자."""
    N = n_side ** 2
    L = math.sqrt(N * math.pi * sigma**2 / (4 * phi))
    a = L / n_side
    pos = np.array([[(i + .5)*a - L/2, (j + .5)*a - L/2, 0.]
                    for i, j in itertools.product(range(n_side), repeat=2)])
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = pos
    f.particles.orientation = [(1, 0, 0, 0)] * N        # 액티브/비구형이면 필요
    f.particles.typeid = [0] * N
    f.particles.types = ['A']
    f.configuration.box = [L, L, 0, 0, 0, 0]            # Lz=0 → 2D (함정 9)
    f.configuration.dimensions = 2
    return f, N, L

f, N, L = frame_2d()
sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=1)
sim.create_state_from_snapshot(f)

cell = md.nlist.Cell(buffer=0.4)
lj = md.pair.LJ(nlist=cell, default_r_cut=2**(1/6), mode='shift')   # WCA (함정 4)
lj.params[('A', 'A')] = dict(epsilon=1.0, sigma=1.0)

bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
integrator.integrate_rotational_dof = False
sim.operations.integrator = integrator
sim.run(100_000)
```

## 조화 트랩 (최소 이미지 적용 — 함정 1·7 해결)

```python
class HarmonicTrap(md.force.Custom):
    """입자별 앵커로 끌어당기는 조화 트랩. 주기 축에만 최소 이미지 적용."""

    def __init__(self, k, anchors, box_L, dimensions=2):
        super().__init__(aniso=False)
        self.k = float(k)
        self.anchors = np.asarray(anchors, dtype=float)
        # 주기 축만 래핑. 2D의 z는 0으로 둬서 제외 (inf 쓰면 inf*0=nan → 함정 7)
        self.period = np.array([box_L, box_L, box_L if dimensions == 3 else 0.0])

    def set_forces(self, timestep):
        with self._state.cpu_local_snapshot as snap, \
             self.cpu_local_force_arrays as arr:
            tags = np.array(snap.particles.tag, copy=True)      # ← tag 인덱싱 필수
            d = np.array(snap.particles.position, copy=True) - self.anchors[tags]
            m = self.period > 0
            d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
            arr.force[:] = -self.k * d
            arr.potential_energy[:] = 0.5 * self.k * (d ** 2).sum(axis=1)

    def displacements(self, state):
        d = np.array(state.get_snapshot().particles.position) - self.anchors
        m = self.period > 0
        d[:, m] -= self.period[m] * np.round(d[:, m] / self.period[m])
        return d
```

**`tags` 인덱싱이 필수**입니다 — HOOMD의 `ParticleSorter`가 메모리 순서를 캐시 효율에 맞게
재배열하므로, 로컬 스냅샷의 순서는 tag 순서가 아닙니다.

**검증**: k=2(가장 취약)에서 `⟨x²⟩` 오차 −0.64%, 힘 배열 NaN 0개, z 힘 정확히 0.
전체 검증은 `scratch/golden_trap.py` (k=2·5·10·20에서 전부 5% 이내, `⟨x²⟩·k` 변동 0.28%).

## 이동 트랩 / 시간 의존 구동

```python
ramp = hoomd.variant.Ramp(A=0.0, B=5.0, t_start=0, t_ramp=1000)
cyc  = hoomd.variant.Cycle(A=-1.0, B=1.0, t_start=0, t_A=10, t_AB=100, t_B=10, t_BA=100)
# set_forces 안에서 center = f(ramp(timestep)) 로 앵커를 움직인다
```
확인: `ramp(0,500,1000) = 0.00, 2.50, 5.00` / `cycle(0,60,160) = -1.00, 0.00, 0.20`

## ABP (연속 회전확산)

```python
active = md.force.Active(filter=hoomd.filter.All())
active.active_force['A'] = (f_a, 0.0, 0.0)      # 입자 로컬 프레임
active.active_torque['A'] = (0.0, 0.0, 0.0)

integrator = md.Integrator(dt=dt, methods=[bd], forces=[lj, active])
integrator.integrate_rotational_dof = False      # ← 함정 3
sim.operations.integrator = integrator
sim.operations.updaters.append(
    active.create_diffusion_updater(trigger=hoomd.trigger.Periodic(1),
                                    rotational_diffusion=D_r))
```

## run-and-flip / run-and-tumble (이산 사건)

ABP와 **다릅니다.** 연속 회전확산이 아니라 포아송 과정입니다. 내장 기능이 없어 커스텀 Action:

```python
class RunAndFlip(hoomd.custom.Action):
    """포아송 과정으로 방향을 180° 반전."""
    def __init__(self, rate, dt, seed=7):
        self.p = rate * dt
        self.rng = np.random.default_rng(seed)
        self.n_flips = 0

    def act(self, timestep):
        with self._state.cpu_local_snapshot as snap:
            q = np.array(snap.particles.orientation, copy=True)
            flip = self.rng.random(len(q)) < self.p
            self.n_flips += int(flip.sum())
            # z축 180° 회전
            q[flip] = np.column_stack([-q[flip,3], -q[flip,2], q[flip,1], q[flip,0]])
            snap.particles.orientation[:] = q

sim.operations.updaters.append(hoomd.update.CustomUpdater(
    action=RunAndFlip(rate=2.0, dt=dt), trigger=hoomd.trigger.Periodic(1)))
```
`rate·dt ≪ 1`이어야 포아송 근사가 성립합니다 (`dt/τ_flip ≤ 1e-2`).

## 사슬 (본드 + 굽힘)

```python
f.bonds.N = M - 1
f.bonds.types = ['backbone']
f.bonds.typeid = [0] * (M - 1)
f.bonds.group = [[i, i+1] for i in range(M-1)]
f.angles.N = M - 2
f.angles.types = ['bend']
f.angles.typeid = [0] * (M - 2)
f.angles.group = [[i, i+1, i+2] for i in range(M-2)]

bond = md.bond.Harmonic();  bond.params['backbone'] = dict(k=100.0, r0=1.0)
angle = md.angle.Harmonic(); angle.params['bend'] = dict(k=10.0, t0=math.pi)
```

> ⚠️ **`t0=math.pi` 를 쓰기 전에 함정 15 를 읽으세요.** 거의 곧은 사슬(`min|θ−π| < 1.41e-3`)
> 에서는 **에너지가 맞는데 힘이 틀립니다.** 뻣뻣한 사슬은 항상 그 영역입니다.
`md.bond.FENEWCA`, `md.angle.CosineSquared`, `md.bond.Table`도 있습니다.

## 임의 퍼텐셜 (r⁻ⁿ 등)

```python
r_min, r_cut, nbins = 0.5, 3.0, 200
r = np.linspace(r_min, r_cut, nbins, endpoint=False)   # ★ endpoint=False (함정 10)
U = A / r**3;  F = 3*A / r**4
U = U - A / r_cut**3                             # 컷오프에서 0. r[-1]이 아니라 r_cut 값으로
tab = md.pair.Table(nlist=cell, default_r_cut=r_cut)
tab.params[('A','A')] = dict(r_min=r_min, U=U, F=F)

# ★ r < r_min 이면 힘이 0이 된다 (함정 11) → 코어를 따로 걸어 도달을 막는다
wca = md.pair.LJ(nlist=cell, default_r_cut=2**(1/6), mode='shift')
wca.params[('A','A')] = dict(epsilon=1.0, sigma=1.0)
integrator = md.Integrator(dt=dt, methods=[bd], forces=[tab, wca])   # 두 힘은 합해진다
```
`F = -dU/dr` 부호에 주의. `U`는 시프트하되 `F`는 시프트하지 않습니다.
**검증** (`scratch/verify_pair_table.py`): 위 조합의 힘·에너지가 해석해와 **0.000%** 일치
(sep = 0.90 ~ 2.00, WCA 영역 포함).

## 비구형 입자

```python
gb = md.pair.aniso.GayBerne(nlist=cell, default_r_cut=4.0)
gb.params[('A','A')] = dict(epsilon=1.0, lperp=0.5, lpar=1.5)   # aspect = 3
bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                         default_gamma=1.0,               # 병진은 스칼라뿐 (위 제약)
                         default_gamma_r=(0.5, 0.5, 2.0)) # 회전은 축별 가능
integrator.integrate_rotational_dof = True
```
`md.pair.aniso.ALJ`, `md.constrain.Rigid`도 가능. **단, 병진 이방성은 어느 쪽도 안 됩니다.**

## 데이터 저장 5계층

```python
# Tier A — 위치/방향 (저빈도, 전체)
hoomd.write.GSD(filename='traj_A.gsd', trigger=hoomd.trigger.Periodic(10_000),
                mode='xb', dynamic=['property'])

# Tier B — per-particle 힘/에너지 (중빈도)
plog = hoomd.logging.Logger(categories=['particle'])
plog.add(lj, quantities=['forces', 'energies'])
hoomd.write.GSD(filename='traj_B.gsd', trigger=hoomd.trigger.Periodic(100_000),
                mode='xb', logger=plog, dynamic=['property', 'momentum'])
# 읽기: frame.log['particles/md/pair/LJ/forces']  → shape (N, 3)

# Tier C — 슬라이딩 윈도우, 사건 발생 시에만 디스크로
burst = hoomd.write.Burst(filename='burst.gsd', trigger=hoomd.trigger.Periodic(10),
                          mode='xb', max_burst_size=1000,
                          write_at_start=True)      # ← 함정 8
# ... 조건 만족 시:  burst.dump()

# Tier D — 추적 입자 소수를 고빈도로
hoomd.write.GSD(filename='tracers.gsd', trigger=hoomd.trigger.Periodic(10),
                mode='xb', filter=hoomd.filter.Tags(list(range(100))),
                dynamic=['property'])

# Tier L — 전역 스칼라
thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
sim.operations.computes.append(thermo)
glog = hoomd.logging.Logger(categories=['scalar', 'sequence'])
glog.add(thermo, quantities=['potential_energy', 'pressure', 'kinetic_temperature'])
hoomd.write.HDF5Log(trigger=hoomd.trigger.Periodic(1000), filename='log.h5',
                    logger=glog, mode='x')
```
`writer.flush()`를 호출해야 파일에 확실히 씁니다 (테스트에서 특히).

## 재시작

```python
hoomd.write.GSD(filename='restart.gsd', trigger=hoomd.trigger.Periodic(100_000),
                mode='wb', truncate=True)          # 항상 1프레임만 유지
# 복원:
sim2.create_state_from_gsd(filename='restart.gsd')
```

## 런타임 가드 (NaN / 에너지 폭발)

```python
class Guard(hoomd.custom.Action):
    def __init__(self, thermo):
        self.thermo = thermo
    def act(self, timestep):
        pe = self.thermo.potential_energy
        if pe is not None and not math.isfinite(pe):
            raise RuntimeError(f"non-finite PE at step {timestep}")

sim.operations.writers.append(hoomd.write.CustomWriter(
    action=Guard(thermo), trigger=hoomd.trigger.Periodic(10_000)))
```

## 초기 겹침 제거

```python
fire = md.minimize.FIRE(dt=1e-4, force_tol=1e-2, angmom_tol=1e-2, energy_tol=1e-7,
                        methods=[...], forces=[lj])
sim.operations.integrator = fire
while not fire.converged:
    sim.run(1000)
sim.operations.integrator = integrator          # 본 적분기로 교체
```

---

# API 빠른 참조 (실측 시그니처)

```
md.methods.Brownian(filter, kT, default_gamma=1.0, default_gamma_r=(1.0,1.0,1.0))
md.methods.Langevin(filter, kT, tally_reservoir_energy=False, default_gamma=..., default_gamma_r=...)
md.methods.OverdampedViscous(filter, default_gamma=1.0, default_gamma_r=(1.0,1.0,1.0))
md.force.Active(filter)                  # .active_force[type], .active_torque[type]
md.force.Custom(aniso=False)             # set_forces(self, timestep) 구현
md.force.Constant(filter)                # .constant_force[type], .constant_torque[type]
md.update.ActiveRotationalDiffusion(trigger, active_force, rotational_diffusion)
md.pair.LJ(nlist, default_r_cut=None, default_r_on=0.0, mode='none', tail_correction=False)
md.pair.Table(nlist, default_r_cut=None)
md.nlist.Cell(buffer=...)                # Stencil, Tree도 있음
md.minimize.FIRE(dt, force_tol, angmom_tol, energy_tol, ...)
hoomd.Simulation(device, seed=None)
hoomd.write.GSD(trigger, filename, filter=All(), mode='ab', truncate=False,
                dynamic=None, logger=None, precision='single')
hoomd.write.Burst(trigger, filename, ..., max_burst_size=-1, write_at_start=False,
                  clear_whole_buffer_after_dump=True)
hoomd.variant.Ramp(A, B, t_start, t_ramp)
hoomd.variant.Cycle(A, B, t_start, t_A, t_AB, t_B, t_BA)
hoomd.filter.Rigid(flags=('center',))    # 'center' | 'constituent' | 'free'
```

`OverdampedViscous`는 **열잡음 없는 과감쇠** — 결정론적 검증(`v=F/γ`, 이완 `τ=γ/k`)에 씁니다.

전체 목록(등방 페어 28종·이방 17종·본드/각도/이면각·마찰·다체·메시·장거리·manifold·HPMC·MPCD):
[docs/hoomd_capabilities.md](../../../docs/hoomd_capabilities.md)

---

# 새 사실을 발견하면

1. 재현 스크립트를 `scratch/`에 남긴다
2. 이 문서에 추가한다 (함정이면 ★ 여부를 판단 — 에러 없이 틀리는가?)
3. `docs/hoomd_capabilities.md`에 실측 수치를 기록한다

**추측을 이 문서에 쓰지 마세요.** 여기 있는 것은 전부 실행으로 확인된 것이어야 합니다.
