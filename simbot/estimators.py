"""S2 예측 엔진 — 해석해와 스케일링. LLM 0줄.

여기 있는 모든 함수는 **시뮬레이션 전에** 답을 내놓는다.
S7이 이 값들과 측정값을 대조한다.
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math

from .units import K_B, kT_si, stokes_drag_si, stokes_einstein_D_si


# =============================================================================
# 조화 트랩 (광집게) — 도메인 A/B
# =============================================================================
@dataclass
class HarmonicTrapPrediction:
    dim: int
    # 입력 (SI)
    T_si: float
    eta_si: float
    radius_si: float
    k_si: float
    # 파생 (SI)
    kT_si: float
    gamma_si: float
    D0_si: float
    sigma_si: float
    tau_trap_si: float
    tau_D_si: float
    l_trap_si: float
    corner_freq_si: float
    var_per_component_si: float
    msd_plateau_si: float
    rms_radial_si: float
    # 무차원
    k_star_sigma: float          # k sigma^2 / kT  — 트랩 강성 (sigma 기준)
    l_trap_over_sigma: float
    tau_sep: float               # tau_D / tau_trap = k_star_sigma
    notes: list[str] = field(default_factory=list)

    def msd_si(self, t_si: float) -> float:
        """<|r(t)-r(0)|^2> = 2 d (kT/k) (1 - exp(-t/tau_trap))."""
        return self.msd_plateau_si * (1.0 - math.exp(-t_si / self.tau_trap_si))


def harmonic_trap(
    *, T_si: float, eta_si: float, radius_si: float, k_si: float, dim: int = 2
) -> HarmonicTrapPrediction:
    """조화 트랩에 갇힌 단일 구형 입자의 과감쇠 BD 예측.

    U(r) = 1/2 k r^2,  r^2 = sum over `dim` 성분

    해석해 (정확, 근사 없음):
      성분별 등분배   <x^2> = kT/k                       (dim 무관)
      반경 제곱       <r^2> = dim * kT/k
      MSD             <|dr(t)|^2> = 2 dim (kT/k)(1 - e^{-t/tau})
      완화시간        tau_trap = gamma/k
      코너주파수      f_c = k/(2 pi gamma)
      구속길이        l_trap = sqrt(kT/k)
    """
    kT = kT_si(T_si)
    gamma = stokes_drag_si(eta_si, radius_si)
    D0 = stokes_einstein_D_si(T_si, gamma)
    sigma = 2.0 * radius_si

    tau_trap = gamma / k_si
    tau_D = sigma**2 / D0
    l_trap = math.sqrt(kT / k_si)
    var_1c = kT / k_si                       # 성분별 <x^2>
    msd_plateau = 2.0 * dim * var_1c         # MSD(t->inf) = 2 <r^2>
    rms_radial = math.sqrt(dim * var_1c)

    k_star_sigma = k_si * sigma**2 / kT

    notes: list[str] = []
    if k_star_sigma > 1e3:
        notes.append(
            f"k*_sigma = {k_star_sigma:.3g} >> 1 — 매우 강한 트랩. "
            f"입자는 자기 지름의 {l_trap/sigma:.2e} 배만 움직인다. 배제부피 무관."
        )
    if tau_D / tau_trap > 10:
        notes.append(
            f"tau_D/tau_trap = {tau_D/tau_trap:.3g}. tau_D 를 기준 시간으로 쓰면 "
            f"dt 가 완화시간을 넘는다 — tau_trap 을 써야 한다."
        )
    return HarmonicTrapPrediction(
        dim=dim, T_si=T_si, eta_si=eta_si, radius_si=radius_si, k_si=k_si,
        kT_si=kT, gamma_si=gamma, D0_si=D0, sigma_si=sigma,
        tau_trap_si=tau_trap, tau_D_si=tau_D, l_trap_si=l_trap,
        corner_freq_si=k_si / (2.0 * math.pi * gamma),
        var_per_component_si=var_1c, msd_plateau_si=msd_plateau,
        rms_radial_si=rms_radial,
        k_star_sigma=k_star_sigma, l_trap_over_sigma=l_trap / sigma,
        tau_sep=tau_D / tau_trap, notes=notes,
    )


# =============================================================================
# 수치 스킴의 계통 오차 — 예측 가능하므로 예측한다
# =============================================================================
def euler_maruyama_trap_variance_bias(dt_star: float) -> float:
    """Euler-Maruyama로 조화 트랩을 적분할 때의 <x^2> 계통 편향 (상대값).

    스킴:      x_{n+1} = x_n (1 - dt*) + sqrt(2 dt*) xi     (tau_trap 단위, D*=1)
    정상분산:  Var* = 2 dt* / (2 dt* - dt*^2) = 1/(1 - dt*/2)
    정확값:    Var* = 1
    ⇒ 상대편향 = dt*/2 / (1 - dt*/2)  ≈ dt*/2

    이 값은 **알려진 오차**다. S7에서 측정 <x^2>가 이만큼 높게 나오는 것이 정상이고,
    그렇지 않으면 다른 문제가 있다는 신호다.
    """
    if not 0 < dt_star < 2:
        raise ValueError(f"dt_star must be in (0,2) for stability; got {dt_star}")
    return 1.0 / (1.0 - dt_star / 2.0) - 1.0


def dt_star_for_trap_bias(target_rel_bias: float) -> float:
    """목표 계통 편향을 만족하는 dt* (tau_trap 단위). 위 식의 역함수."""
    b = target_rel_bias
    return 2.0 * b / (1.0 + b)


# =============================================================================
# 과감쇠·연속체 근사 타당성 (S3 게이트)
# =============================================================================
@dataclass
class ValidityChecks:
    tau_inertial_si: float
    inertial_ratio: float        # tau_inertial / tau_process — << 1 이어야 과감쇠
    reynolds: float              # << 1 이어야 Stokes
    passed: bool
    failures: list[str] = field(default_factory=list)


def overdamped_validity(
    *, gamma_si: float, mass_si: float, tau_process_si: float,
    velocity_scale_si: float, radius_si: float, eta_si: float, rho_fluid_si: float,
    inertial_tol: float = 1e-2, reynolds_tol: float = 1e-2,
) -> ValidityChecks:
    """BD(과감쇠 + Stokes)의 전제를 검사한다."""
    tau_i = mass_si / gamma_si
    ratio = tau_i / tau_process_si
    Re = rho_fluid_si * velocity_scale_si * radius_si / eta_si
    fails: list[str] = []
    if ratio > inertial_tol:
        fails.append(f"관성 시간척도 비 {ratio:.3g} > {inertial_tol:g} — Langevin 검토 필요")
    if Re > reynolds_tol:
        fails.append(f"Reynolds {Re:.3g} > {reynolds_tol:g} — Stokes 항력 의심")
    return ValidityChecks(tau_i, ratio, Re, not fails, fails)


# =============================================================================
# 통계 정밀도 — 몇 개 표본이 필요한가
# =============================================================================
def samples_for_variance_precision(rel_err: float) -> float:
    """분산 추정의 상대 표준오차를 rel_err 이하로 만드는 독립표본 수.

    가우시안 표본의 분산 추정량: SE(s^2)/s^2 = sqrt(2/(n-1)) ≈ sqrt(2/n)
    """
    return 2.0 / rel_err**2


def seeds_for_target_sigma(*, diff: float, se_diff: float, k_current: int,
                           n_sigma: float = 3.0, t_correction: bool = True) -> dict:
    """관측된 차이를 `n_sigma` 로 분해하는 데 필요한 시드 수. **돌리기 전에 계산한다.**

    ★ 이것을 먼저 계산하지 않으면 달성 불가능한 검정을 요구하게 된다
      (CLAUDE.md §테스트 4규칙 ④: "설계 검정력이 3σ를 못 만드는 곳에서 3σ 기각을
      요구하지 않는다"). 그리고 유의해질 때까지 시드를 늘리는 것은 **optional
      stopping** 이라 p 값을 망친다 — 목표 `k` 를 **미리 고정**하는 데 쓴다.

    `SE ∝ 1/√k` 를 가정한다. 시드 평균의 표준오차가 그렇게 줄기 때문이다.

    Args:
        t_correction: `True` 면 `n_sigma` 를 `t(ν)` 분위수로 환산한다.
            `k` 가 작으면 `SE` 자체가 크게 흔들리므로 정규 분위수는 과소추정이다
            (근거: findings/tolerance-from-a-4-seed-se-is-not-a-3-sigma-test.md).
    """
    if se_diff <= 0.0:
        raise ValueError(f"se_diff = {se_diff} — 양수여야 한다")
    if k_current < 2:
        raise ValueError(f"k_current = {k_current} — 최소 2개")
    from scipy.stats import t as _t

    sigma_now = abs(diff) / se_diff

    def needed(target_quantile: float) -> float:
        #  SE ∝ 1/√k  ⇒  k_needed = k_current · (target/현재σ)²
        return k_current * (target_quantile / sigma_now) ** 2

    #  정규 기준
    k_normal = needed(n_sigma)
    out = {"sigma_now": sigma_now, "k_current": k_current,
           "n_sigma": n_sigma, "k_needed_normal": k_normal}

    if t_correction:
        #  t 분위수는 ν = k−1 에 의존하고 ν 는 k 에 의존한다 → 자기일관 반복
        two_sided_p = 2.0 * (1.0 - _norm_cdf(n_sigma))
        k = max(k_normal, 3.0)
        for _ in range(60):
            q = float(_t.ppf(1.0 - two_sided_p / 2.0, df=max(k - 1.0, 1.0)))
            k_new = needed(q)
            if abs(k_new - k) < 1e-9:
                k = k_new
                break
            k = k_new
        out.update({"t_quantile": q, "k_needed": k,
                    "two_sided_p": two_sided_p})
    else:
        out["k_needed"] = k_normal
    out["k_needed_int"] = int(math.ceil(out["k_needed"]))
    out["se_diff_target"] = abs(diff) / out.get("t_quantile", n_sigma)
    return out


def _norm_cdf(z: float) -> float:
    return 0.5 * (1.0 + math.erf(z / math.sqrt(2.0)))


def trap_run_length(
    *, n_particles: int, rel_err_target: float, tau_trap_si: float,
    decorrelation_in_tau: float = 2.0,
) -> dict[str, float]:
    """비상호작용 입자 N개를 같은 트랩에 넣었을 때 필요한 런 길이.

    입자들이 상호작용하지 않으므로 **스냅샷 하나가 독립표본 N개**다
    (Barakat 2022 방식: knowledge/source/papers/2022-barakat-enhanced-dispersion-harmonic-traps.md).
    시간 방향 독립성은 ~2 tau_trap 마다 확보된다.
    """
    n_needed = samples_for_variance_precision(rel_err_target)
    n_timepoints = max(1.0, n_needed / n_particles)
    t_total_tau = n_timepoints * decorrelation_in_tau
    return {
        "independent_samples_needed": n_needed,
        "independent_timepoints_needed": n_timepoints,
        "t_total_in_tau_trap": t_total_tau,
        "t_total_si": t_total_tau * tau_trap_si,
    }


# =============================================================================
# 비용 추정 — knowledge/wiki/findings/local-cpu-parallelism.md 실측 상수
# =============================================================================
THROUGHPUT_PARTICLE_STEPS_PER_S = 6.3e6
EFFICIENCY_BY_K = {1: 1.0, 2: 0.948, 3: 0.926, 4: 0.925, 5: 0.774,
                   6: 0.696, 8: 0.617, 10: 0.547, 12: 0.443}


def estimate_wall_time_s(n_particles: int, n_steps: int, concurrency: int = 1) -> float:
    """단일 프로세스 wall time [s]. 동시 실행 k개일 때 프로세스당 소요시간.

    ⚠ 이 모델은 N >= 500 에서 실측되었다. N 이 작으면 스텝당 고정 오버헤드가
      지배하므로 과소추정한다. 작은 N 은 pilot 런으로 실측할 것.
    """
    eta = EFFICIENCY_BY_K.get(concurrency)
    if eta is None:
        ks = sorted(EFFICIENCY_BY_K)
        lo = max([k for k in ks if k <= concurrency], default=ks[0])
        eta = EFFICIENCY_BY_K[lo]
    return n_particles * n_steps / (THROUGHPUT_PARTICLE_STEPS_PER_S * eta)
