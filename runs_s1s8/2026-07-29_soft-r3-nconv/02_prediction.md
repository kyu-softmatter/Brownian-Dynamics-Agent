# S2 PREDICTION — N 수렴 + ψ₆ 유한크기 지수  (자동 생성: scripts/soft2d_nconv_predict.py)
#
# ⚠ 이 파일은 S5 실행 전에 봉인된다 (SEALED.sha256). 실행 후 수정 금지.
# 허용오차는 부모 런(N=100)의 시드 SE 를 스케일링해서 유도했다.
#
question: N=100 → 256 에서 ① ψ₆ 가 1/√N 으로 줄어드는가 (= 유한크기 바닥인가) ② 국소량(결함 분율·ψ₆ 국소)이 변하지 않는가 ③ 절단오차가 줄어드는가
parent_run: runs/2026-07-29_soft-r3-time-resolved
card_open_item: 카드 §9 (N ≥ 252 요구) · §10 (g₆ 지수 η₆ 미구현)
items:
- quantity: psi6_global__A0.1__N256
  value: 0.029821626374541662
  tolerance: ±0.00508785
  basis: '액체 가설 `|⟨ψ₆⟩| ~ N^-1/2`: 0.0477146 × (2.56)^-0.5 = 0.0298216. **경쟁 가설(hexatic, η₆=1/4 → p=1/16)은
    0.0449921** 이고 둘의 간격이 0.0151705 = 허용오차의 3.0배다 — 판별 가능하다'
  discriminates: A=0.1 의 ψ₆ 가 유한크기 바닥인가 진짜 준장거리 질서인가
  competing_value: 0.0449921046953466
- quantity: psi6_exponent_p__A0.1
  value: 0.5
  tolerance: ±0.15
  basis: '`p = -dln|⟨ψ₆⟩|/dlnN`. 액체면 0.5, KTHNY hexatic 경계면 0.0625. 허용오차 ±0.15 는 두 점 추정의 오차 전파(부모 SE 0.001919
    상대 4.0%, /ln2.56=0.9400) 규모다'
  discriminates: A=0.1 의 상 판독 (η₆ = 4p)
  competing_value: 0.0625
- quantity: defect_fraction__A0.1__N256
  value: 0.6444999999999999
  tolerance: ±0.0120275
  basis: '**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 0.6445 ± 0.002835. 허용오차 3√2·SE. 어긋나면 `N=100` 에 유한크기 효과가
    있었다는 뜻이다'
  discriminates: defect_fraction 가 N=100 에서 이미 수렴했는가
  competing_value: null
- quantity: psi6_local__A0.1__N256
  value: 0.3611549796164035
  tolerance: ±0.00517919
  basis: '**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 0.361155 ± 0.001221. 허용오차 3√2·SE. 어긋나면 `N=100` 에 유한크기 효과가
    있었다는 뜻이다'
  discriminates: psi6_local 가 N=100 에서 이미 수렴했는가
  competing_value: null
- quantity: psi6_global__A1__N256
  value: 0.03732159153005341
  tolerance: ±0.00656739
  basis: '액체 가설 `|⟨ψ₆⟩| ~ N^-1/2`: 0.0597145 × (2.56)^-0.5 = 0.0373216. **경쟁 가설(hexatic, η₆=1/4 → p=1/16)은
    0.0563074** 이고 둘의 간격이 0.0189858 = 허용오차의 2.9배다 — 판별 가능하다'
  discriminates: A=1 의 ψ₆ 가 유한크기 바닥인가 진짜 준장거리 질서인가
  competing_value: 0.05630735669569704
- quantity: psi6_exponent_p__A1
  value: 0.5
  tolerance: ±0.15
  basis: '`p = -dln|⟨ψ₆⟩|/dlnN`. 액체면 0.5, KTHNY hexatic 경계면 0.0625. 허용오차 ±0.15 는 두 점 추정의 오차 전파(부모 SE 0.002477
    상대 4.1%, /ln2.56=0.9400) 규모다'
  discriminates: A=1 의 상 판독 (η₆ = 4p)
  competing_value: 0.0625
- quantity: defect_fraction__A1__N256
  value: 0.5488999999999999
  tolerance: ±0.0129684
  basis: '**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 0.5489 ± 0.003057. 허용오차 3√2·SE. 어긋나면 `N=100` 에 유한크기 효과가
    있었다는 뜻이다'
  discriminates: defect_fraction 가 N=100 에서 이미 수렴했는가
  competing_value: null
- quantity: psi6_local__A1__N256
  value: 0.3916287848353386
  tolerance: ±0.00358632
  basis: '**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 0.391629 ± 0.0008453. 허용오차 3√2·SE. 어긋나면 `N=100` 에 유한크기
    효과가 있었다는 뜻이다'
  discriminates: psi6_local 가 N=100 에서 이미 수렴했는가
  competing_value: null
- quantity: psi6_global__A10__N256
  value: 0.15493777482770385
  tolerance: ±0.0210782
  basis: '액체 가설 `|⟨ψ₆⟩| ~ N^-1/2`: 0.2479 × (2.56)^-0.5 = 0.154938. **경쟁 가설(hexatic, η₆=1/4 → p=1/16)은
    0.233756** 이고 둘의 간격이 0.078818 = 허용오차의 3.7배다 — 판별 가능하다'
  discriminates: A=10 의 ψ₆ 가 유한크기 바닥인가 진짜 준장거리 질서인가
  competing_value: 0.23375574821979261
- quantity: psi6_exponent_p__A10
  value: 0.5
  tolerance: ±0.15
  basis: '`p = -dln|⟨ψ₆⟩|/dlnN`. 액체면 0.5, KTHNY hexatic 경계면 0.0625. 허용오차 ±0.15 는 두 점 추정의 오차 전파(부모 SE 0.007949
    상대 3.2%, /ln2.56=0.9400) 규모다'
  discriminates: A=10 의 상 판독 (η₆ = 4p)
  competing_value: 0.0625
- quantity: defect_fraction__A10__N256
  value: 0.294975
  tolerance: ±0.0287428
  basis: '**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 0.294975 ± 0.006775. 허용오차 3√2·SE. 어긋나면 `N=100` 에 유한크기 효과가
    있었다는 뜻이다'
  discriminates: defect_fraction 가 N=100 에서 이미 수렴했는가
  competing_value: null
- quantity: psi6_local__A10__N256
  value: 0.5687030825018882
  tolerance: ±0.0175345
  basis: '**국소량이므로 `N` 에 무관해야 한다.** `N=100` 관측 0.568703 ± 0.004133. 허용오차 3√2·SE. 어긋나면 `N=100` 에 유한크기 효과가
    있었다는 뜻이다'
  discriminates: psi6_local 가 N=100 에서 이미 수렴했는가
  competing_value: null
- quantity: beta_u_at_rcut__A0.1__N256
  value: 0.0002
  tolerance: ±0.0005
  basis: '`r_cut = 0.98·L/2 − buffer`, `L* = √256 = 16` → `r_cut = 7.740`. `βU(r_cut) = 0.1/7.740³`. 결정론적이므로
    오차는 반올림뿐이다. `N=100` 에서는 `r_cut = 4.800` 이었다'
  discriminates: 카드 §9 의 절단오차가 실제로 줄었는가
  competing_value: null
- quantity: beta_u_at_rcut__A1__N256
  value: 0.0022
  tolerance: ±0.0005
  basis: '`r_cut = 0.98·L/2 − buffer`, `L* = √256 = 16` → `r_cut = 7.740`. `βU(r_cut) = 1/7.740³`. 결정론적이므로
    오차는 반올림뿐이다. `N=100` 에서는 `r_cut = 4.800` 이었다'
  discriminates: 카드 §9 의 절단오차가 실제로 줄었는가
  competing_value: null
- quantity: beta_u_at_rcut__A10__N256
  value: 0.0216
  tolerance: ±0.0005
  basis: '`r_cut = 0.98·L/2 − buffer`, `L* = √256 = 16` → `r_cut = 7.740`. `βU(r_cut) = 10/7.740³`. 결정론적이므로
    오차는 반올림뿐이다. `N=100` 에서는 `r_cut = 4.800` 이었다'
  discriminates: 카드 §9 의 절단오차가 실제로 줄었는가
  competing_value: null
- quantity: coord_kinds_aggregate__A10__N256
  value: 3
  tolerance: ±0
  basis: '`N=100` 집계에서 3종. **집계 추정량**을 쓴다 — 프레임별 문턱은 `N` 에 의존한다 (입자 1개가 `N=100` 에서 1 %, `N=256` 에서 0.39
    %). 근거: findings/fraction-threshold-flips-meaning-*.md'
  discriminates: 결함의 성격이 N 에 무관한가
  competing_value: null
regimes:
  n_ref: 100
  n_new: 256
  n_ratio: 2.56
  seeds:
  - 21
  - 22
  - 23
  - 24
  amplitudes:
  - 0.1
  - 1.0
  - 10.0
  gamma_zahn:
    A0.1: 0.5568327996831708
    A1: 5.568327996831708
    A10: 55.68327996831708
  L_star_ref: 10.0
  L_star_new: 16.0
  coverage: 0.08726646259971647
  note_coverage: ★ 커버리지는 N 에 무관하다 (n* = 1 이라 입자당 면적이 정확히 d²). 상자만 √N 으로 커진다 — L = 150 → 240 µm
  L_si_ref: 0.00015000000000000001
  L_si_new: 0.00024000000000000003
  tau_d_si: 2292.427238845818
  exponent_liquid: 0.5
  exponent_hexatic_boundary: 0.0625
  eta6_hexatic_boundary: 0.25
alternatives:
- '★ 가장 중요한 경쟁 가설: A=10 의 ψ₆ = 0.248 이 **유한크기 바닥**일 가능성. p ≈ 0.5 가 나오면 카드 §8.3 의 ''hexatic-유사'' 판독을 약화시켜야
  한다 — 0.248 은 그저 N=100 의 바닥이었다는 뜻이 된다.'
- '두 점으로는 멱함수 **형태**를 검증할 수 없다. 지수만 뽑는 것이고, 그 한계를 psi6_finite_size_exponent 의 n_points 가 나른다. 형태를 논하려면
  셋 이상(예: N=64·144·400)이 필요하다.'
- A=0.1 과 A=1 은 p ≈ 0.5 가 나와야 정상이다 (깊은 액체). 만약 그것도 작게 나오면 지수 추정 자체가 틀린 것이므로 **A ≤ 1 이 대조군 역할**을 한다.
- 에너지/입자는 N 사이에서 직접 비교할 수 없다 — power_law_table 이 U(r_cut) = 0 으로 이동시키므로 r_cut 이 바뀌면 이동량도 바뀐다. 기록은 하되 예측
  항목으로 세우지 않았다.
- 커버리지 대조군(3.71 %)은 **돌리지 않는다.** 축약 config 가 비트 단위로 같아서 산술 항등식이다 — tests/test_s7_structure.py::test_coverage_does_not_touch_the_reduced_config_at_all
  로 고정했다.
