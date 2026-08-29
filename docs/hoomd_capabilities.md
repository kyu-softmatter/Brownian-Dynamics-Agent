# HOOMD 능력 매트릭스 — 우리가 만들 수 있는 물리 모듈

> 조사일 2026-08-03 · **hoomd 7.1.0** (설치본; 문서는 v7.1.1 — 패치 차이, 조사한 API는 동일)
> 환경: macOS 26.5 arm64, `simulation_bot` conda env
> 빌드: `gpu=False`, `mpi=False`, `compile flags: DOUBLE[SINGLE]`
> 방법: 설치본 introspection(`scratch/survey.py`) + 실동작 스모크 테스트(`scratch/smoke.py`)
>
> **이 문서는 추측이 아니라 실측입니다.** 마스터플랜 §5.6 모듈 레지스트리 설계의 근거.

---

## 0. 요약

- **스코프를 열어둔 것이 옳았습니다.** HOOMD는 예상보다 훨씬 넓습니다:
  등방 페어 28종, 이방(비구형) 페어 17종, 본드/각도/이면각, 마찰 접촉, 다체 퍼텐셜,
  메시(막·vesicle), 장거리 정전기, 곡면 구속(manifold), 강체, HPMC, MPCD까지.
- **`intake/`의 5개 케이스 전부 구현 가능**합니다 (§3).
- **구멍은 하나**: 광집게(조화 트랩)가 내장에 없습니다 → `md.force.Custom`으로 직접 구현.
  실동작 확인 완료, 골든 물리 검증까지 통과(§4).
- **가장 중요한 제약**: **BD에서 막대의 이방성 병진 마찰은 기하만으로 만들 수 없습니다** — 강체로
  묶든 비드사슬로 묶든 `γ⊥/γ∥ = 1.000000`. 이방성은 유체역학적 상호작용(HI)의 효과이고 BD에는
  HI가 없기 때문입니다. HOOMD의 한계가 아니라 모델 자체의 성질입니다 (§5.1). 회전 마찰은
  텐서로 정상 동작합니다 (§5.2). 그 외: MPI 없음(단일 코어), GPU 없음.

---

## 1. 실동작 검증 결과 (15개 API)

| # | 항목 | 결과 | 대응 설계 |
|---|---|---|---|
| 1 | BD + WCA 2D | ✓ | §11, 부록 A |
| 2 | `GSD(logger=)` per-particle 힘 저장 | ✓ `(64,3)` 배열 확인 | §9.2 Tier B |
| 3 | `write.Burst` 슬라이딩 윈도우 + `dump()` | ✓ (단, `write_at_start=True` 필요) | §9.2 Tier C |
| 4 | `md.force.Custom` 조화 트랩 | ✓ | `external.harmonic_trap` |
| 5 | `variant.Ramp` / `Cycle` 시간의존 구동 | ✓ | `driving.*` |
| 6 | `bond.Harmonic` + `angle.Harmonic` 사슬 | ✓ | `bonded.*` |
| 7 | `force.Active` + `create_diffusion_updater` | ✓ | `active.abp` |
| 8 | 커스텀 Action updater (run-and-flip) | ✓ 5 flips (기대 3.6) | `active.run_and_flip` |
| 9 | `pair.Table` 임의 r⁻ⁿ | ✓ | `pair.table` |
| 10 | `pair.aniso.GayBerne` + 축별 `gamma_r` | ✓ aspect 3.0 | `shape.ellipsoid` |
| 11 | `methods.OverdampedViscous` (결정론적) | ✓ | 골든 테스트용 |
| 12 | `write.HDF5Log` 전역 스칼라 | ✓ | §9.2 Tier L |
| 13 | restart GSD (`truncate=True`) | ✓ 1프레임 유지 + 재로드 | §14 체크포인트 |
| 14 | 커스텀 Action 런타임 감시 | ✓ | §12.5 런타임 가드 |
| 15 | `GSD(filter=Tags)` 추적 입자 서브셋 | ✓ 8/64 | §9.2 Tier D |

**로우데이터 5계층(§9)이 전부 실동작 확인되었습니다.**

### 조화 트랩 골든 물리 검증 (`scratch/golden_trap.py`)

2D 조화 트랩 안 BD 입자의 평형 분포는 `⟨x²⟩ = k_BT/k`, 이완시간 `τ = γ/k`.
N=400, 각 340,000 스텝:

| k | τ=γ/k | dt | ⟨x²⟩ 측정 | 예측 kT/k | 오차 |
|---|---|---|---|---|---|
| 2.0 | 0.500 | 2.5e-4 | 0.50188 ± 0.0016 | 0.50000 | +0.38% |
| 5.0 | 0.200 | 1.0e-4 | 0.20112 ± 0.0006 | 0.20000 | +0.56% |
| 10.0 | 0.100 | 5.0e-5 | 0.09998 ± 0.0003 | 0.10000 | −0.02% |
| 20.0 | 0.050 | 2.5e-5 | 0.04993 ± 0.0001 | 0.05000 | −0.14% |

`⟨x²⟩·k = [1.0038, 1.0056, 0.9998, 0.9986]` — k를 10배 바꿔도 kT로 일정 (변동계수 **0.28%**).

→ **BD 적분기가 정확하고, `trap-2d-5um`을 골든 물리 테스트로 채택할 수 있습니다.**

---

## 2. ⚠️ 새로 발견한 함정 (마스터플랜 §11 함정 목록에 추가)

### 함정 7 — 외부 힘 + 주기경계: 최소 이미지를 적용하지 않으면 조용히 틀린다 ⭐️

고정 앵커를 향한 트랩을 짤 때 `d = pos - anchor`를 그대로 쓰면, 입자가 박스를 넘어
wrap되는 순간 거리가 L만큼 점프하고 거대한 **잘못된 방향의** 복원력을 받습니다.

증상이 고약합니다 — **터지지 않고 조용히 틀립니다.** 트랩이 강하면(경계 도달 확률이 낮으면)
정확한 값이 나오고, 약할수록 오차가 커집니다:

| k | 최소 이미지 없음 | 최소 이미지 적용 |
|---|---|---|
| 2.0 | +1856% ✗ | +0.38% ✓ |
| 5.0 | +344% ✗ | +0.56% ✓ |
| 10.0 | +0.16% ✓ | −0.02% ✓ |

k=10만 테스트했다면 통과했다고 착각했을 것입니다. **`external.*` 모듈 전체에 해당**하며,
검증기에 "외부 힘 모듈은 최소 이미지 적용 여부를 선언해야 함" 규칙이 필요합니다.

```python
d = pos - anchors[tags]
d -= L * np.round(d / L)      # ← 이 한 줄
```

### 함정 8 — 2D에서 최소 이미지의 z성분 NaN

위 코드에서 z 방향 박스 길이를 `np.inf`로 두면 `inf * round(0/inf) = inf*0 = nan`이 되어
z 힘이 NaN이 됩니다. 2D라 HOOMD가 z를 무시해 결과는 맞았지만, **NaN이 힘 배열에 들어가는 것
자체가 위험**합니다(§12.5 런타임 가드가 오탐할 수 있음). 주기 차원에만 마스크로 적용하세요.

### 함정 9 — `write.Burst`는 새 파일에 `write_at_start=True` 필요

없으면 `RuntimeError: Must set write_at_start to write to a new file.`
결과적으로 파일에 초기 프레임 1개가 추가로 들어갑니다(버퍼 10 + 시작 1 = 11).

---

## 3. `intake/` 5개 케이스 매핑 (전부 구현 가능)

| 케이스 | 필요 물리 | HOOMD 대응 | 난이도 |
|---|---|---|---|
| `trap-2d-5um` | 조화 트랩 | `md.force.Custom` (내장 없음) ✓검증 | 낮음 |
| `trap-drag-2d-hex300` | 트랩 + 이동 구동 | `force.Custom` + `variant.Ramp/Cycle` ✓검증 | 낮음 |
| `chain-bend-2d-oscill` | 본드 + 굽힘 + 진동 | `bond.Harmonic` + `angle.Harmonic` + `variant.Cycle` ✓검증 | 낮음 |
| `soft-r3-2d-A-sweep` | 소프트 r⁻³ | `pair.Table` (또는 `pair.Mie`) ✓검증 | 낮음 |
| `abp-rod-2d-run-flip` | 타원체 + run-and-flip | `pair.aniso.GayBerne`/`ALJ` + 커스텀 updater ✓검증 | **중간** |

`abp-rod`가 유일하게 까다롭습니다 — 이방성 병진 마찰 문제 때문입니다 (§5.1~5.4).
다만 2D라 회전은 z축 하나뿐이고, MSAD와 장시간 MSD는 정확히 재현됩니다 (§5.3).

---

## 4. 전체 능력 목록 (실측)

### 등방 페어 퍼텐셜 — `md.pair` (28종)
`Buckingham` `DLVO` `DPD` `DPDConservative` `DPDLJ` `Ewald` `ExpandedGaussian`
`ExpandedLJ` `ExpandedMie` `ForceShiftedLJ` `Fourier` `Gaussian` **`LJ`** `LJ0804`
`LJ1208` `LJGauss` `Mie` `Moliere` `Morse` `OPP` `ReactionField` `TWF` **`Table`**
`WangFrenkel` **`Yukawa`** `ZBL` `Zetterling`

> WCA 전용 클래스는 여전히 없음 → `LJ(r_cut=2^(1/6)σ, mode='shift')` (§11 함정 1)

### 이방 페어 (비구형 입자) — `md.pair.aniso` (17종)
**`GayBerne`** (타원체) · **`ALJ`** (다면체/타원체 이방 LJ) · `Dipole` · `YLZ` ·
`Patchy` `PatchyLJ` `PatchyGaussian` `PatchyMie` `PatchyYukawa`
`PatchyExpandedLJ` `PatchyExpandedGaussian` `PatchyExpandedMie`

### 마찰 접촉 — `md.pair.friction` ✓ 실측 (2026-08-06)
`FrictionLJCoulomb` `FrictionLJCoulombNewton` `FrictionLJLinear`
(Hofmann et al. 2025, [arXiv:2507.16388](https://doi.org/10.48550/arXiv.2507.16388))

`params = {epsilon, sigma, kT, gamma_f?, kappa_f?}` · 법선은 WCA 고정 ·
접선력 = `w(r)·f(u)`, `w(r) = −dU_WCA/dr`, `u` = 접촉점 상대 표면속도 ·
힘과 **토크**를 같이 낸다. 실측 (`scratch/verify_pair_friction.py`, 23/23, kT=0, V=1, r*=0.95):

| 모델 | `f(u)` | F_tan | 잔류 슬립 u/V=1e−4 |
|---|---|---|---|
| `FrictionLJLinear` | `γ_f·u` | 6.1926 | 6.19e−4 (∝u) |
| `FrictionLJCoulomb` | `κ_f` | 18.578 | **18.578 (안 줄어듦)** |
| `FrictionLJCoulombNewton` | `min[γ_f u, κ_f]` | 1.0 | 1e−4 (∝u) |

**세 가지 구조적 성질** (전부 실측):
- **컷오프 밖에서 정확히 0** — `w(r)=0`. DLVO 2차극소(`r*=1.00759`)는 밖이라 `F_tan=|τ|=0.000e+00`.
- **산일성** — 접선 변위가 있어도 `u=0`이면 힘 0 (법선 WCA는 6.1926으로 살아있음).
  `K″`에만 기여하고 `K′`에는 0 → 준정적 극한에서 소멸. 정지마찰(stick) 상태가 없다.
- **무슬립 구름에 면제** — `ω=V/2R`이면 3종 전부 정확히 0. 굽힘(=구름)을 막지 못한다.

> ⛔ **JKR 굽힘강성의 대체재가 아니다.** 필요한 것은 이력 있는 접선 스프링 + 구름 저항
> (Cundall–Strack / Dominik–Tielens)이고 **HOOMD에 없다** → `force.Custom(aniso=True)`.
> ★ 함정 19·20 (지렛대 `R = diameter/2` · `Brownian` 속도)을 skill `bd-hoomd`에서 먼저 읽으세요.

### 결합 상호작용
| 종류 | 클래스 |
|---|---|
| `md.bond` | **`Harmonic`** **`FENEWCA`** `Table` `Tether` |
| `md.angle` | **`Harmonic`** `CosineSquared` `Table` |
| `md.dihedral` | `OPLS` `Periodic` `Table` |
| `md.improper` | `Harmonic` `Periodic` |
| `md.special_pair` | `Coulomb` `LJ` |

### 외부 힘
| 종류 | 클래스 | 비고 |
|---|---|---|
| `md.external.field` | `Electric` `Magnetic` `Periodic` | **조화 트랩 없음** |
| `md.external.wall` | `LJ` `Gaussian` `Morse` `Mie` `Yukawa` `ForceShiftedLJ` | + `hoomd.wall`: `Plane` `Sphere` `Cylinder` |
| `md.force` | **`Custom`** **`Constant`** **`Active`** `ActiveOnManifold` | `Custom`이 확장의 주 통로 |

### 적분 방법 — `md.methods`
**`Brownian`** `Langevin` **`OverdampedViscous`** `ConstantVolume` `ConstantPressure`
`DisplacementCapped`
`md.methods.rattle`: 위 방법들 + **Manifold 구속** (곡면 위 BD)
`md.methods.thermostats`: `Berendsen` `Bussi` `MTTK`

### 업데이터 · 구속
`md.update`: **`ActiveRotationalDiffusion`** `ReversePerturbationFlow` `ZeroMomentum` `Mesh*`
`md.constrain`: **`Rigid`** `Distance`
`hoomd.update`: **`BoxResize`** `CustomUpdater` `FilterUpdater` `RemoveDrift`

#### `BoxResize` 실측 (`scratch/verify_3d_boxresize.py`, 2026-08-06) ⭐️

```
hoomd.update.BoxResize(trigger, box: variant.box.BoxVariant, filter=All())
hoomd.variant.box.Interpolate(initial_box, final_box, variant)      # variant: 0→1 스칼라
hoomd.variant.box.InverseVolumeRamp(initial_box, final_volume, t_start, t_ramp)
```

| 확인한 것 | 결과 |
|---|---|
| 좌표를 아핀 스케일하는가 (`r → s·r`) | **예 — 최대 오차 8.9e-16** (문서 주장 그대로) |
| 목표 크기에서 정확히 멈추는가 | 예 (`Lx` 오차 < 1e-6) |
| `Brownian` + `pair.Table` + WCA 와 공존 | 예. 크래시 없음, PE 유한, 셀리스트 자동 갱신 |
| 압축 후 `r_cut < L/2` | 별도로 확인해야 함 (압축이 이 조건을 깰 수 있다 — 함정 6) |

★ **아핀 스케일이 곧 함정입니다** — 압축은 **이미 결합된 쌍의 결합길이도 같이 줄입니다.**
좁은 우물(DLVO 2차극소)에서는 한 트리거의 변형이 우물 폭을 넘으면 쌍이 장벽 안쪽으로
밀려 **1차극소(접촉)로 비가역 붕괴**합니다. 실측 문턱 → 아래 §5.5.

### 다체 · 장거리 · 메시
`md.many_body`: `Tersoff` `SquareDensity` `RevCross`
`md.long_range.pppm`: `Coulomb`
`md.mesh`: `bending.Helfrich` `bending.BendingRigidity` `bond.*` `conservation.{Area,Volume,TriangleArea}`
> 막·vesicle 시뮬레이션 가능

### 시간 의존 구동 — `hoomd.variant`
`Constant` **`Ramp`** **`Cycle`** `Power`
`variant.box`: `Interpolate` `InverseVolumeRamp`

### 출력 — `hoomd.write`
**`GSD`** **`Burst`** **`HDF5Log`** `Table` `DCD` `CustomWriter`

### 필터 — `hoomd.filter`
`All` **`Tags`** `Type` `Rigid` `Union` `Intersection` `SetDifference` `Null` `CustomFilter`

### 기타 축의 확장
`hoomd.hpmc` — 하드입자 몬테카를로 (`Ellipsoid` `ConvexPolyhedron` `Sphinx` … + `external.Harmonic`)
`hoomd.mpcd` — 다입자 충돌 동역학 (유체역학 결합)
`md.alchemy` — 알케미컬 변환
`md.tune` / `hoomd.tune` — 자동 튜닝 (`NeighborListBuffer` `ParticleSorter` `GridOptimizer`)

---

## 5. 제약 사항 (모듈 설계에 반영 필수)

| 제약 | 내용 | 영향 |
|---|---|---|
| **병진 마찰이 스칼라 — 우회 불가** ⭐️ | 아래 §5.1에서 실측 확인. 강체·비드사슬 어느 쪽도 이방성을 만들지 못함 | 이방성 병진확산은 **명시적으로 부과**해야 함. `shape.*` 모듈이 한계를 선언 |
| **회전 마찰은 텐서 — 정상 동작** ✓ | `default_gamma_r=(x,y,z)` 실측 확인 (§5.1) | MSAD·회전 동역학은 정확히 재현 가능 |
| **MPI 없음** | `mpi_enabled=False` | 단일 런 = 단일 코어 확정. 병렬성은 "여러 런 동시 실행"에서만 (§2 원칙 6 유효) |
| **GPU 없음** | `gpu_enabled=False` | osx-arm64에 CUDA 빌드 없음. 규모 한계 유지 |
| **조화 트랩 미내장** | `md.external.field`에 없음 | `force.Custom` 필수 + 함정 7(최소 이미지) 준수 |
| **`hpmc.external.Harmonic`은 못 씀** | HPMC(몬테카를로) 전용 | MD/BD 경로와 호환 안 됨 |

### 5.1 이방성 병진 마찰 — 실측 결과 ⭐️

**결론: BD에서 막대의 이방성 병진 마찰은 기하만으로 나오지 않는다. 이건 HOOMD의 한계가 아니라
BD 모델 자체의 성질이다.**

두 경로를 모두 실측했습니다 (`scratch/rigid_rod_friction.py`, `scratch/anisotropy_probe.py`).
방향을 고정한 막대에 일정한 힘을 걸고 종단속도로 `γ = F/v`를 직접 측정 (`OverdampedViscous`, 잡음 없음):

| 구성 | γ∥ | γ⊥ | γ⊥/γ∥ | 판정 |
|---|---|---|---|---|
| **강체** (`constrain.Rigid`, 구 5개) | 1.00000 | 1.00000 | **1.000000** | 등방 |
| **자유배수 비드 막대** (본드+각도로 묶은 구 5개) | 5.00000 | 5.00000 | **1.000000** | 등방 |
| 이론 기대 (slender body, HI 있음) | — | — | → 2 | — |

추가 확인: 강체의 구성 입자 타입에 `gamma['A'] = 10`을 줘도 병진 속도가 전혀 변하지 않습니다
(v∥ = 1.000000 동일). **강체의 병진 항력은 중심 입자의 스칼라 `gamma` 하나로 완전히 결정**됩니다.

**왜 그런가**: 자유배수(free-draining) 근사에서 각 비드는 독립적으로 Stokes 항력을 받으므로
총 항력은 `N·γ_bead`로 방향에 무관합니다. 실제 막대의 `γ⊥/γ∥ → 2`는 **비드끼리 서로의 유동장을
가려주는 유체역학적 상호작용(HI)** 때문이고, 축방향 정렬일 때 가림이 더 강해 `γ∥ < γ⊥`가 됩니다.
**HI가 없는 BD에서는 기하를 아무리 정교하게 만들어도 이방성이 생기지 않습니다.**

측정값 `γ = 5.000 = N·γ_bead`가 이 설명과 정확히 일치합니다.

### 5.2 회전 마찰 텐서 — 정상 동작 ✓

`default_gamma_r=(x,y,z)` 축별 지정은 완벽히 동작합니다. 일정 토크에 대한 회전각 측정:

| `gamma_r` | τ∥x 회전각 | τ∥z 회전각 | 비 z/x | 기대 |
|---|---|---|---|---|
| (1, 1, 1) | 2.00000 rad | 2.00000 rad | 1.0000 | 1.0 ✓ |
| (1, 1, 5) | 2.00000 rad | 0.40000 rad | **0.2000** | 0.2 ✓ |

→ **MSAD와 회전 동역학은 정확히 재현할 수 있습니다.**

### 5.3 `abp-rod-2d-run-flip`에 대한 실제 영향

2D이므로 회전은 z축 하나뿐 — `gamma_r` 이방성이 애초에 필요 없습니다.

| 관측량 | 의존 | 정확도 |
|---|---|---|
| **MSAD** | `γ_r,z` 만 | ✓ 정확 |
| **MSD (장시간, t ≫ τ_r)** | 등방 평균 `γ̄` 만 (방향이 평균됨) | ✓ 정확 |
| **MSD (단시간, t < τ_r)** | `γ∥`, `γ⊥` 개별 | ✗ **이방성 손실** |

run-and-flip은 flip 사이에 방향이 오래 유지되므로 단시간 이방성이 실제로 관측될 수 있습니다.
그 신호가 목표라면 아래 선택지 중 하나가 필요합니다.

### 5.4 선택지

| 옵션 | 방법 | 비용 | 정확도 |
|---|---|---|---|
| **A. 등방 평균 γ̄ (권장 기본)** | Perrin/slender-body 이론으로 종횡비에서 `γ̄ = 3/(1/γ∥ + 2/γ⊥)` 계산. 회전은 `γ_r` 텐서로 정확히. 모듈이 한계를 선언 | 낮음 | 장시간 MSD·MSAD 정확, 단시간 이방성 없음 |
| **B. 커스텀 적분기로 이방성 부과** | 결정론적 항은 보정력으로 가능하나, **잡음항이 문제** — HOOMD가 등방 잡음을 생성하므로 요동-소산 정리가 깨짐. `md.half_step_hook` 경로 필요 | 높음 | 정확 |
| **C. MPCD로 HI 도입** | `hoomd.mpcd` (다입자 충돌 동역학) — HI가 자연히 생김 | 매우 높음 | 정확, 다른 축의 프로젝트 |

→ **A를 기본값으로 두고, `shape.ellipsoid` 모듈이 `translational_friction: "isotropic_average"`를
선언**하게 합니다. 검증기가 "단시간 MSD 이방성을 목표로 하는데 등방 근사를 쓰고 있음" 경고를
낼 수 있습니다. B가 필요해지면 별도 모듈로 추가합니다 (모듈 구조의 이점).

### 5.5 박스 압축이 좁은 우물의 결합을 부순다 — 실측 문턱 ⭐️

`scratch/verify_3d_boxresize.py` (2026-08-06). `network` 케이스가 "응집 후 스퀴즈"
프로토콜이라 실행 전에 재봤습니다. 단위 d=kT=γ=1.

DLVO 원장 (`chain-bend-2d-dlvo` 승계): 2차극소 `h*=0.007593`, 장벽 `h*=0.000508`,
결합길이 `ℓ*=1.007593`, `k_bond*=1.042e6 kT/d²`, `τ_bond*=9.594e-7 τ_B`.

**해석 예측**: 아핀 스텝이 쌍을 장벽 안쪽으로 밀면 되돌아오지 못하므로
```
ε_crit = (h_min* − barrier_h*) / ℓ*  =  0.703 %  (트리거당 선형 변형)
```

**실측** (kT=0 결정론적, 트리거 간격 = 21 τ_bond 이완):

| ε/트리거 | 최종 `h*` | 결과 |
|---|---|---|
| 0.40 % | 0.007591 | 2차극소 유지 ✓ |
| 0.80 % | **−0.042257** | 붕괴 — 접촉/겹침 |
| 2.00 % | −0.042519 | 붕괴 |

**실측 문턱 0.40~0.80% 가 예측 0.703% 를 포함합니다.** 붕괴 후 `h*≈−0.042`(4.2% 겹침)는
vdW 발산과 WCA 코어가 균형하는 자리이고 **비가역**입니다 — 가역적 2차극소 결합이
영구 접촉으로 바뀝니다.

**따라오는 비용**: φ 0.02→0.10 은 총 선형변형 41.5% 이므로 ε=0.4% 로 나누면
**최소 134 단계**입니다. 각 단계를 준정적(≫τ_B)으로 하면 ~400 τ_B 가 되고,
`dt = 1e-2·τ_bond = 9.6e-9 τ_B` 에서 4.2e10 스텝 — N=1528 에서 측정한
**6384 steps/s** 기준 **1800 시간**입니다. → 준정적 압축은 이 dt 에서 불가능.
빠른 압축(단계 간 이완을 τ_bond 규모만)은 134×2000 = 2.7e5 스텝 ≈ **42초**로 가능하지만,
그때 구조는 확산이 개입하지 않은 **아핀 압축된 초기 배치**입니다 — 다른 물리입니다.

---

## 6. 재현 방법

```bash
conda env update -f environment.yml -n simulation_bot --prune
CONDA=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
$CONDA scratch/survey.py     # 능력 조사 (설치본 introspection)
$CONDA scratch/smoke.py      # 15개 API 실동작
$CONDA scratch/golden_trap.py         # 조화 트랩 골든 물리 검증 (~2분)
$CONDA scratch/rigid_rod_friction.py  # 강체 막대 이방성 (등방으로 나옴)
$CONDA scratch/anisotropy_probe.py    # gamma_r 이방성 + 자유배수 비드막대
```
