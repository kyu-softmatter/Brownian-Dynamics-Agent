# S2 PREDICTION — 시뮬레이션 전에 답을 적는다
#
# ⚠ 이 파일은 S5 실행 전에 봉인된다 (SEALED.sha256). 실행 후 수정 금지.
# 모든 수치는 simbot 함수 출력이다 — 손계산 없음:
#   simbot.spec.derive · simbot.estimators.euler_maruyama_trap_variance_bias
#   simbot.analysis.trap.em_uniform_noise_excess_kurtosis
#
# tolerance 형식: ±X% (예측값 상대) · ±X (절대) · >X / <X (단측 경계)
# competing_value: 경쟁 가설. 설계 검정력 |예측-경쟁|/SE 계산에 쓰인다.

items:
  - quantity: var_x_star
    value: 1.0025062656641603
    tolerance: ±1%
    basis: Euler-Maruyama 정상분산 1/(1-dt*/2). estimators.euler_maruyama_trap_variance_bias
    discriminates: 적분기 스킴이 EM 인가 exact 인가
    competing_value: 1.0
  - quantity: var_x_nm2
    value: 414.19469999999995
    tolerance: ±1.5%
    basis: 등분배 <x^2> = kT/k (정확해, 근사 없음). a·eta 에 무관
    discriminates: k 와 T 만으로 결정 — 가정에 독립
    unit: nm^2
  - quantity: var_r_nm2
    value: 828.3893999999999
    tolerance: ±1.5%
    basis: <r^2> = d kT/k
    discriminates: 차원 (2D vs 3D 는 3/2 배 차이)
    unit: nm^2
  - quantity: msd_plateau_star
    value: 4.010025062656641
    tolerance: ±2%
    basis: 'MSD(t->inf) = 2d <x*^2> (축약 단위).  ★ <x*^2> 가 EM 편향을 갖고 있으므로
      plateau 도 그만큼 높다: 2d(1+bias) = 4 x 1.0025063 = 4.010025.
      정확히 2d = 4 를 예측하면 3.54 sigma 어긋난다 — 알려진 편향을 예측에 넣지 않은
      것이고, tolerance(±2%) 가 넓어서 PASS 로 가려진다 (2026-07-28 발견)'
    note: 'plateau 와 <x*^2> 는 독립 경로(시계열 vs 스냅샷)이지만 **같은 편향을 공유한다**.
      경로 독립성은 plateau_over_2d_var 이 따로 검사한다'
  - quantity: tau_star
    value: 1.0
    tolerance: ±5%
    basis: tau_fit/tau_trap = 1 (카드 단위에서 정의상)
    discriminates: MSD 피팅이 gamma/k 를 회수하는가
  - quantity: tau_trap_ms
    value: 8.064405509911927
    tolerance: ±5%
    basis: tau_trap = gamma/k = 6 pi eta a / k. eta·a 에 선형
    discriminates: eta 와 a 를 검증하는 유일한 창
    unit: ms
  - quantity: msd_r_squared
    value: '>0.99'
    tolerance: '>0.99'
    basis: MSD 가 단일지수 형태 plateau(1-exp(-t/tau))
  - quantity: kT_conf_star
    value: 1.0
    tolerance: ±2%
    basis: 배위 온도 <|grad U|^2>/<lap U> = kT. guards.configurational_temperature
    note: 순수 조화 트랩에서는 <x^2> 검사와 대수적으로 동일 — 독립 검사가 아니다
  - quantity: plateau_over_2d_var
    value: 1.0
    tolerance: ±1%
    basis: '자기일관성: MSD plateau(시계열) / 2d<x^2>(스냅샷) = 1'
    discriminates: 분석 코드의 오류 — 독립 경로 둘이 어긋나면 코드 문제
regimes:
  k_star_sigma: 241432.3505346641
  l_trap_over_sigma: 0.002035177387846081
  tau_D_over_tau_trap: 241432.35053466403
  note: k*_sigma 가 경계(~1)에서 5.4 decade 떨어져 있다 — 극도로 강한 구속. 배제부피가 동역학에 무관하다.
alternatives:
  - dt* 가 크면 <x*^2> 가 1/(1-dt*/2) = 1.00251 만큼 높게 나온다 — 이것이 정상이다.
  - 위치분포 첨도는 정확히 3 이 아니라 2.9940 다 (HOOMD 노이즈가 균일분포이므로). 정확히 3.000 이면 오히려 이상하다.
  - tau_trap 이 예측과 어긋나면 eta 가정 또는 a 해석(A2) 이 틀렸다는 신호다. <x^2> 는 영향받지 않으므로 둘을 분리해서 판별할 수 있다.
  - 'INCONCLUSIVE 예상: var_x_star 는 EM 편향 0.25 % 가 통계오차(시드 4개, ~0.45 %)보다 작아 exact 스킴과 구별되지 않는다. 사전에 예견된
    한계다.'
