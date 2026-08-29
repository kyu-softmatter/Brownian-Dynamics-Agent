---
type: source
kind: repo
lab_authored: true
title: "graybox_abp_mpc — Learning continuum-level closures for control of interacting active particles (코드)"
authors:
  - "Quah T"
  - "Takatori SC"
  - "Rawlings JB"
year: 2024
source_url: "https://github.com/titusswsquah/graybox_abp_mpc"
homepage: "https://titusswsquah.github.io/graybox_abp_mpc/"
license: MIT
commit: b6df95712be7
commit_date: 2025-07-10
raw_file: knowledge/raw/lab/graybox_abp_mpc/
access: open (MIT)
engine: "HOOMD-blue 3.8.1 (BD, ABP)"
reproduced: no
parameters_extracted: yes
paper: knowledge/source/papers/2025-quah-continuum-closures-active-control.md
ingested_at: 2026-07-27
ingested_by: agent
tags:
  - "direct"
---

# graybox_abp_mpc — BD 시뮬레이션 부분

Quah 2025 논문([arXiv:2501.18809](https://arxiv.org/abs/2501.18809))의 공개 코드.
**MIT 라이선스.** 논문에 없던 `Δt`·초기배치·평형화 프로토콜이 전부 여기 있다.

> **범위:** 이 문서는 **Brownian dynamics 시뮬레이션을 만드는 부분까지만** 다룬다.
> MPC 제어기(`lib/hd_controller.py`), 신경 ODE(`lib/hd_neural_ode.py`), 학습 파이프라인은
> 우리 v1 범위(`D3`) 밖이라 제외한다.

**핵심 파일**

| 파일 | 크기 | 역할 |
|---|---|---|
| `lib/hd_simulator.py` | 38 KB | **BD 설정 전체** — 단위계, 스냅샷, 적분기, 힘, 박스 압축 |
| `lib/hoomd_utils.py` | 3.6 KB | 쿼터니언 변환, `XFilter` (위치 기반 CustomFilter) |
| `environment.yml` | — | **HOOMD 3.8.1**, gsd 2.5.2, Python 3.9, numba, cupy 경로 |

---

## 1. 단위계 — 능동물질 고유 방식 ★

`hd_simulator.py:300-315`. **우리 설계와 방향이 반대**라 특히 볼 가치가 있다.

```python
ell    = 1                              # 런 길이 ℓ  ← 길이 단위
D_r    = 1                              # 회전확산   ← 시간 단위 τ_r = 1/D_r
zeta_t = 1                              # 병진 항력  ← 질량 단위
delta  = ell / sqrt(ld2)                # δ: 미시 확산길이,  ld2 = (ℓ/δ)²
pcl_rad = sqrt(3/4) * delta             # a
eta    = zeta_t / (6*pi*pcl_rad)        # 용매 점도 ← Stokes 항력에서 역산
zeta_r = 8*pi*eta*pcl_rad**3            # 회전 항력 ← Stokes
kbt    = D_r * zeta_r                   # ★ kT를 회전 Stokes–Einstein–Debye에서 유도
D_t    = 4/3 * pcl_rad**2 * D_r         # 병진 확산
u0     = sqrt(ld2 * D_t / tau)          # 자기추진 속도
fp     = zeta_t * u0                    # 능동력
plate_width = n_pcls*pi*pcl_rad**2 / (plate_gap * vol_frac)   # 목표 면적분율에서 박스 폭 역산
```

### 우리와 반대 방향이다

| | 이 코드 | 우리 설계 (`D7`) |
|---|---|---|
| 기준 길이 | **런 길이 `ℓ`** | 입자 지름 `σ` |
| 기준 시간 | **`τ_r = 1/D_r`** (회전) | `τ_D = σ²/D_t` (병진 확산) |
| `kT` | **유도량** — `kT = D_r·ζ_r` | 입력값 |
| `D_t` | 유도량 — `(4/3)a²D_r` | `kT/γ` |

**능동 입자에서는 회전이 주인공**이므로 이 선택이 자연스럽다. 하지만 수동 콜로이드에는
맞지 않는다. → `03_units_nondim.md`에 **"능동계 단위계"** 를 별도 절로 둘 것.

### 내부 일관성 검산 — 통과 ✅

```
D_t = (4/3)a²D_r  vs  kT/ζ_t = D_r·ζ_r/ζ_t = D_r·8πηa³/(6πηa) = (4/3)a²D_r   ✓ 항등
```

수치 확인 (`ld2 = 100`): `D_t = 0.01`, `kT/ζ_t = 0.01` — **정확히 일치.**

> 이건 우리 `bdkit/units/`가 해야 할 검산의 실제 예다. 선행 프로젝트 `config.py:105`의
> *"`kT`·`γ`·`D_t` 셋이 다 주어지면 서로 맞아야 한다"* 규칙과 같은 정신.

### `ld2 = 100`에서의 실제 수치

| 양 | 값 |
|---|---|
| `δ` | 0.1 |
| `a` (반지름) | 0.0866 · `σ = 2a` = **0.173** |
| `η` | 0.613 |
| `ζ_t` = 1 · `ζ_r` = 0.01 |
| **`kT`** | **0.01** |
| `D_t` | 0.01 |
| `u0` | **1** · `f_p` = 1 |
| `τ_r` = 1 · **`τ_D = σ²/D_t` = 3** (= 3 τ_r) |
| **`Pe = u0τ_r/σ`** | **5.77** |

---

## 2. HOOMD 설정 (3.8.1)

`hd_simulator.py:413-437`

```python
device = hoomd.device.GPU() if hoomd.version.gpu_enabled else hoomd.device.CPU()
sim = hoomd.Simulation(device=device, seed=self.md_seed)
sim.create_state_from_snapshot(snapshot)

integrator = hoomd.md.Integrator(dt=self.md_dt, integrate_rotational_dof=True)
brownian = hoomd.md.methods.Brownian(kT=self.kbt, filter=hoomd.filter.All())
brownian.gamma.default   = self.zeta_t
brownian.gamma_r.default = [self.zeta_r, self.zeta_r, self.zeta_r]
integrator.methods.append(brownian)
sim.operations.integrator = integrator

active1 = hoomd.md.force.Active(filter=hoomd.filter.All())
active1.active_force["A"] = (self.fp, 0, 0)
integrator.forces.append(active1)

if self.interact:
    cell = hoomd.md.nlist.Cell(buffer=0.4)
    lj_pcl = hoomd.md.pair.LJ(nlist=cell)
    lj_pcl.params[("A","A")] = dict(epsilon=self.epsilon,
                                    sigma=(2*self.pcl_rad / (2**(1/6))))
    lj_pcl.r_cut[("A","A")] = 2*self.pcl_rad
    integrator.forces.append(lj_pcl)
```

### ★ 회전확산을 얻는 두 가지 방법 — 이 코드는 다른 쪽을 골랐다

| 방법 | 이 코드 | 선행 slit 프로젝트 |
|---|---|---|
| 구현 | `Brownian(gamma_r=ζ_r)` + `integrate_rotational_dof=True` | `Active.create_diffusion_updater(rotational_diffusion=D_r)` |
| `D_r`의 성격 | **열적** — `D_r = kT/ζ_r`로 묶임 | **독립** — `kT`와 무관하게 지정 |
| `create_diffusion_updater` 호출 | **없음** (grep 0건) | 있음 |

**물리적으로 다르다.** Takatori 2020은 *"회전확산이 반드시 열적 기원일 필요는 없다 —
박테리아는 편모로 재배향한다"* 고 명시한다. 합성 Janus 입자는 열적, 박테리아는 비열적.

> → `04_hoomd_mapping.md`에 **"ABP 회전확산: 두 경로와 선택 기준"** 항목을 넣을 것.
> 조용히 잘못 고르면 `D_r`이 `kT`에 묶여버려 원하는 `Pe`를 못 만든다.

### ★ WCA를 LJ 절단으로 만든다 — 에너지가 shift되지 않는다

```python
sigma_LJ = d / 2**(1/6)        # d = 2·pcl_rad = 입자 지름
r_cut    = d                   # = 2^(1/6)·sigma_LJ = LJ 최솟점
```

절단점이 정확히 LJ 최솟점이라 **힘은 연속(그 점에서 0)** 이다. 하지만 `mode="shift"`가 없어서
**퍼텐셜 에너지는 `-ε`만큼 불연속**이다.

- **BD에는 문제없다** — 과감쇠 동역학은 힘만 쓴다
- **에너지를 보고하면 틀린다** — 쌍마다 `ε`씩 오프셋이 쌓인다
- 우리 S8 진단의 *"퍼텐셜에너지 드리프트"* 항목은 이런 계에서 **절대값이 아니라 변화율**만
  봐야 한다

> `04_hoomd_mapping.md`에 기록: **WCA = `md.pair.LJ` + `σ→σ/2^{1/6}` + `r_cut=σ` + `mode="shift"`.**
> `mode="shift"`를 빼면 힘은 맞고 에너지는 틀린다.

### 2D 설정

```python
snapshot.configuration.box = [width, plate_gap, 0, 0, 0, 0]   # Lz=0 → 2D
snapshot.particles.moment_inertia = ones((N,3)) * 2/5          # 균질 구, m=1
snapshot.particles.orientation = 쿼터니언(θ ~ U[0,2π])
```

---

## 3. 초기배치와 평형화 — 그대로 채택할 만한 프로토콜 ★

`hd_simulator.py:379-505`

```python
positions = utils.get_fcc_pts(width_scale*plate_width, plate_gap, 2*pcl_rad, xp)
positions = positions[random.choice(positions.shape[0], n_pcls, replace=False), :]
```

1. **FCC 격자점을 만들고 무작위 부분추출** → **겹침이 구조적으로 불가능**
2. 박스를 `width_scale = 3.0` 배 넓게 시작 (= 목표 밀도의 1/3)
3. **`BoxResize` + `variant.Ramp`로 목표 밀도까지 서서히 압축**

```python
ramp_time = init_gsd_time / 5                       # = 20 τ_r
t_ramp    = int(ramp_time / md_dt)
ramp      = hoomd.variant.Ramp(A=0, B=1, t_start=0, t_ramp=t_ramp)
box_resize = hoomd.update.BoxResize(box1=initial_box, box2=final_box,
                                    variant=ramp, trigger=hoomd.trigger.Periodic(10))
n_steps = init_gsd_time / md_dt + 1                 # = 200,001 스텝
```

> **우리 S6 PREFLIGHT의 "초기배치 겹침" 검사와 S9 REPAIR의 "겹침 완화" 규칙에 직접 대응한다.**
> 현재 우리 규칙표는 `DisplacementCapped`를 쓰라고 되어 있는데, **격자→부분추출→박스압축**이
> 더 깔끔하다. 겹침을 만들었다가 푸는 게 아니라 애초에 안 만든다.
>
> → `06_repair_policy.md` 갱신 후보.

---

## 4. 파라미터 기본값

`hd_simulator.py:130-190`

| 파라미터 | 값 | 무차원 환산 |
|---|---|---|
| `n_pcls` | **1e4** | 논문과 일치 |
| `vol_frac` | **0.2** | 면적분율 (2D) |
| `md_dt` | **5e-4** | `dt/τ_r = 5e-4` · **`dt/τ_D = 1.67e-4`** |
| `md_samp_dt` | 0.01 | 20 스텝마다 샘플 |
| `md_ctrl_dt` | 1e-2 | 제어 주기 |
| `md_seed` | **10** | |
| `epsilon` | **5** | **`ε/kT = 500`** ← 매우 강함 |
| `plate_gap` | 10 | `y` 방향 폭 |
| `width_scale` | 3.0 | 초기 박스 확대 배율 |
| `init_gsd_time` | 100 | = 200,000 스텝 평형화 |
| `interact` | True | `False`면 이상기체 ABP |
| nlist buffer | **0.4** | 선행 slit 프로젝트와 동일 |

---

## 5. ★ 가장 중요한 발견 — 우리 `dt` 게이트가 잘못 맞춰져 있다

이 코드의 `dt`를 우리 게이트에 넣어보면 **기각된다.** 그래서 실제로 돌아가는 시뮬레이션
세 개를 나란히 놓고 봤다.

| 출처 | `dt/τ_D` | `√(2D·dt)/σ` | `ε/kT` |
|---|---|---|---|
| 선행 slit 프로젝트 | `1.0e-3` | 0.045 | — (harmonic) |
| Xu 2023 (dynamic interfaces) | `2.0e-5` | 0.0063 | 1 |
| **이 코드 (Quah)** | **`1.67e-4`** | **0.018** | **500** |
| | **50배 폭** | **7배 폭** | |

**우리 게이트 `dt/τ_D ≤ 1e-4`로 판정하면:**

| | 판정 |
|---|---|
| 선행 slit 프로젝트 | ★ **기각** (1.0e-3) |
| Xu 2023 | 통과 |
| 이 코드 | ★ **기각** (1.67e-4) |

**실제로 돌아가고 논문까지 나온 시뮬레이션 3개 중 2개를 기각한다.**
사용자 본인의 선행 프로젝트까지 포함해서.

**스텝당 변위 기준 `√(2D·dt) ≤ 0.1σ`로 판정하면 셋 다 통과한다** (0.006–0.045σ).

> ### 제안: 게이트를 `dt/τ_D` 고정값이 아니라 **스텝당 변위**로 바꿀 것
>
> `dt/τ_D`는 계마다 50배씩 벌어지는데 변위는 7배 안에 든다. 그리고 변위는
> **물리적으로 의미가 있다** — 한 스텝에 입자가 지름의 몇 %를 움직이는가.
>
> 실제로 우리 S8 진단에는 이미 *"스텝당 최대 변위 < 0.1σ"* 가 있다.
> **S3의 사전 게이트와 S8의 사후 진단이 같은 양을 봐야 일관된다.**
>
> 능동계에서는 이류 변위 `u₀·dt/σ`도 함께 봐야 한다 (이 코드에서 0.0029σ).

이 발견은 [[dt-gate-should-be-displacement-based]] 로 별도 기록한다.

---

## 6. `ε/kT = 500`인데 WCA가 안정적이다 (`D8` 관련)

선행 프로젝트 `forces.py:117`은 *"WCA의 `r⁻¹³` 코어는 오버댐프에서 위험하다"* 고 경고하며
bounded harmonic을 권했다. 그런데 이 코드는 **`ε/kT = 500`이라는 매우 강한 WCA**를 쓰고도
문제없이 돈다.

**차이는 `dt`가 아니다** — 오히려 이 코드가 선행 프로젝트보다 `dt/τ_D`는 6배 작지만
스텝당 변위는 0.018σ 대 0.045σ로 **2.5배 작다.**

> **`D8` 잠정 결론:** WCA의 위험성은 퍼텐셜 자체가 아니라 **스텝당 변위와의 결합** 문제다.
> 변위를 0.02σ 이하로 유지하면 `ε/kT = 500`도 안정적이다.
> 확정하려면 재현 실험이 필요하다 — `reproduced: no` 상태.

---

## 재현 계획

| 단계 | 내용 |
|---|---|
| 1 | **HOOMD 3.8.1 → 7.1.0 API 이식.** `gamma.default` 문법 등 확인 필요 |
| 2 | `interact=False`(이상기체 ABP)로 `D_t = kT/ζ_t` 검산 — 가장 싼 검사 |
| 3 | `vol_frac=0.2`, `ld2=100`으로 `dt` 스윕 → **변위 기준 게이트 실측 보정** |
| 4 | `ε/kT` 스윕 → `D8` 확정 |

> **재현 상태: `reproduced: no`.** 위 값은 코드에서 읽은 것이고 우리 환경에서 돌린 것이 아니다.
> 검증 근거로 인용하지 않는다 — 계약: [`../../wiki/CLAUDE.md`](../../wiki/CLAUDE.md)
