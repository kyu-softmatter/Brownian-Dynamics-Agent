"""절단 반경(`r_cut`) 제안 — 관습이 아니라 **명시된 오차 허용치**에서 유도한다.

## 왜 별도 모듈인가

`r_cut = 2.5 σ` 는 LJ 액체 문헌의 관습이고, 우리 계에는 근거가 없다.
지식 베이스가 이 문제의 실제 비용을 보여준다:

- `knowledge/wiki/benchmarks/choi2020-interfacial-rdf.md` — *"절단 거리가 미확정이다.
  원 논문이 언급하지 않는다. **재현 실패 시 1순위 용의자.**"*
- `knowledge/wiki/systems/interfacial-colloid--equilibrium-structure.md` §7 —
  `r_cut/d ≥ 21–69` 이 필요한데 minimum image 는 `r_cut ≤ L/2 = 65 d` 를 강제한다.
  **박스가 상호작용 사거리를 담지 못하는 계가 실재한다.**

그래서 이 모듈은 세 가지를 항상 함께 돌려준다:
  ① `r_cut` 값  ② **어느 제약이 지배했는가**  ③ 절단 지점에서 남은 에너지·힘

## 규약

| 포텐셜 | `r_cut` | 성격 |
|---|---|---|
| WCA | `2^{1/6} σ_LJ` | **선택이 아니라 정의.** LJ 최솟점 |
| 절단 LJ · Morse · Yukawa · DLVO | 오차 허용치에서 수치 해 | 선택 — 근거를 기록해야 한다 |

절대 제약: `r_cut ≤ L/2` (minimum image)
"""
from __future__ import annotations

from dataclasses import dataclass, field
import math
from typing import Callable

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

WCA_RCUT_OVER_SIGMA_LJ = 2.0 ** (1.0 / 6.0)   # 1.1224620483...

# =============================================================================
# 허용치 프리셋 — **관측량이 허용치를 정한다**
# =============================================================================
# 문헌 관습 `r_cut = 2.5 sigma` 는 LJ 액체의 구조 계산에 맞춰진 값이다.
# beta*eps=1 에서 그 지점의 잔여량을 실측하면 beta U = 1.63e-2, beta F sigma = 3.90e-2 다.
# 즉 관습값은 "beta U ~ 2e-2 까지는 봐준다"는 선택이고, 그것이 모든 관측량에
# 충분한 것은 아니다.
#
# 그래서 허용치를 관측량으로 키(key) 한다. 값을 고르면 **근거가 함께 기록된다.**
TOLERANCE_PRESETS: dict[str, tuple[float, float]] = {
    # 이름          : (beta_u_tol, beta_f_sigma_tol)
    "convention":    (2.0e-2, 4.0e-2),   # 문헌 관습 2.5 sigma 에 상당. 구조 계산용
    "structure":     (1.0e-2, 2.0e-2),   # g(r), 배위수, S(q)
    "thermodynamics": (1.0e-3, 1.0e-2),  # B2, 상경계, 삼투압 — 기본값
    "precision":     (1.0e-4, 1.0e-3),   # 임계점 근처, 정밀 자유에너지
}
DEFAULT_TOLERANCE_PRESET = "thermodynamics"
DEFAULT_BETA_U_TOL, DEFAULT_BETA_F_SIGMA_TOL = TOLERANCE_PRESETS[DEFAULT_TOLERANCE_PRESET]

# 비용 비교 기준. 셀 리스트에서 이웃 수는 r_cut^dim 에 비례한다.
CONVENTIONAL_LJ_RCUT = 2.5


@dataclass
class CutoffProposal:
    """`r_cut` 제안 — 값 + 지배 제약 + 잔여 오차. 셋을 분리하지 않는다."""

    potential: str
    r_cut_over_sigma: float
    criterion: str                     # 어느 제약이 값을 정했는가
    beta_U_at_cut: float               # 절단점에 남은 포텐셜 [kT]
    beta_F_sigma_at_cut: float         # 절단점에 남은 힘 [kT/sigma]
    min_box_over_sigma: float          # minimum image 가 요구하는 최소 L
    exact: bool                        # 포텐셜 정의상 고정된 값인가
    tolerance_preset: str = ""         # 어느 프리셋을 썼는가
    rationale: str = ""
    warnings: list[str] = field(default_factory=list)
    hoomd_note: str = ""

    def cost_vs_convention(self, dim: int = 3) -> float:
        """관습값 `2.5 sigma` 대비 이웃 쌍 수 배수 (= 쌍힘 계산 비용 배수).

        셀 리스트에서 이웃 수는 `r_cut^dim` 에 비례한다. 이 값이 `run_policy.yaml`
        의 wall time 예산과 직결되므로 항상 함께 보고한다.
        """
        return (self.r_cut_over_sigma / CONVENTIONAL_LJ_RCUT) ** dim

    def box_ok(self, box_over_sigma: float) -> tuple[bool, str]:
        """minimum image 검사. `r_cut <= L/2` 여야 한다."""
        if box_over_sigma >= self.min_box_over_sigma:
            return True, ""
        return False, (
            f"minimum image 위반: r_cut/sigma = {self.r_cut_over_sigma:.4g} 이므로 "
            f"L/sigma >= {self.min_box_over_sigma:.4g} 가 필요한데 {box_over_sigma:.4g} 다. "
            f"박스를 키우거나 r_cut 을 줄이고 그 대가를 기록할 것."
        )


# =============================================================================
# WCA — 선택이 아니다
# =============================================================================
def wca_cutoff() -> CutoffProposal:
    """WCA 의 `r_cut` 은 LJ 최솟점 `2^{1/6} σ_LJ` 로 **정의상 고정**이다.

    그 점에서 힘이 정확히 0 이므로 **힘이 연속**이다 — 절단 인공물이 없다.
    포텐셜은 `-ε` 만큼 어긋나므로 `mode="shift"` 로 올려야 순수 척력이 된다.

    ⚠ 랩 코드(Quah `graybox_abp_mpc`)는 `mode="shift"` 가 **없다** —
      힘은 맞지만 에너지에 `-ε` 오프셋이 남는다. BD 동역학에는 영향이 없으나
      포텐셜 에너지를 로깅하거나 상태방정식을 계산하면 틀린다.
      근거: knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md §"WCA를 LJ 절단으로"
    """
    return CutoffProposal(
        potential="WCA",
        r_cut_over_sigma=WCA_RCUT_OVER_SIGMA_LJ,
        criterion="potential_definition",
        beta_U_at_cut=0.0,
        beta_F_sigma_at_cut=0.0,
        min_box_over_sigma=2.0 * WCA_RCUT_OVER_SIGMA_LJ,
        exact=True,
        rationale=(
            "WCA = LJ 를 최솟점에서 잘라 올린 것. r_cut = 2^{1/6} sigma_LJ 는 정의이고, "
            "그 점에서 F = 0 이라 힘이 연속이다. 오차 허용치가 개입할 여지가 없다."
        ),
        hoomd_note=(
            'md.pair.LJ(nlist, default_r_cut=2**(1/6), mode="shift"); '
            'params = dict(epsilon=eps, sigma=sigma_LJ).  '
            '입자 지름 d 를 쓰려면 sigma_LJ = d / 2**(1/6) 로 두고 r_cut = d.'
        ),
    )


def barker_henderson_diameter(beta_epsilon: float, sigma_lj: float = 1.0) -> float:
    """WCA 의 유효 경구 지름  d_BH = int_0^{r_min} [1 - exp(-beta u(r))] dr.

    `phi` 를 Carnahan-Starling 같은 경구 이론에 넣으려면 이 보정이 필요하다.
    근거: knowledge/wiki/findings/wca-reproduces-carnahan-starling.md
      — *"phi 를 CS 에 바로 넣으면 안 된다. beta*eps=1 에서 d_eff = 1.017 sigma"*
    """
    if beta_epsilon <= 0:
        raise ValueError("beta_epsilon must be > 0")
    r_min = WCA_RCUT_OVER_SIGMA_LJ * sigma_lj

    def integrand(r: float) -> float:
        s6 = (sigma_lj / r) ** 6
        bu = beta_epsilon * (4.0 * (s6 * s6 - s6) + 1.0)
        return 1.0 - math.exp(-bu) if bu < 700 else 1.0

    # r -> 0 에서 integrand -> 1 이고 매끄럽다. 하한을 0 으로 두면 exp 오버플로가 나므로
    # 해석적으로 1 인 구간을 분리한다.
    r_lo = 0.3 * sigma_lj
    val, _ = quad(integrand, r_lo, r_min, limit=200)
    return r_lo + val


# =============================================================================
# 오차 허용치에서 r_cut 을 푸는 일반 루틴
# =============================================================================
def cutoff_from_tolerance(
    *,
    potential: str,
    beta_u: Callable[[float], float],
    beta_f_sigma: Callable[[float], float],
    r_start: float,
    r_max: float = 100.0,
    preset: str = DEFAULT_TOLERANCE_PRESET,
    beta_u_tol: float | None = None,
    beta_f_sigma_tol: float | None = None,
    hoomd_note: str = "",
) -> CutoffProposal:
    """`|βU(r_cut)| ≤ beta_u_tol` **그리고** `|βFσ(r_cut)| ≤ beta_f_sigma_tol` 을
    만족하는 가장 작은 `r_cut` (sigma 단위).

    두 제약 중 어느 것이 지배했는지 기록한다 — 이것이 없으면 `r_cut` 을 왜 그 값으로
    잡았는지 나중에 재구성할 수 없다.

    Args:
        preset: `TOLERANCE_PRESETS` 의 키. 관측량이 허용치를 정한다.
        beta_u_tol / beta_f_sigma_tol: 명시하면 preset 을 덮어쓴다.
    """
    if preset not in TOLERANCE_PRESETS:
        raise KeyError(f"unknown preset {preset!r}; known: {sorted(TOLERANCE_PRESETS)}")
    p_u, p_f = TOLERANCE_PRESETS[preset]
    beta_u_tol = p_u if beta_u_tol is None else beta_u_tol
    beta_f_sigma_tol = p_f if beta_f_sigma_tol is None else beta_f_sigma_tol
    def solve(f: Callable[[float], float], tol: float) -> float | None:
        """`r > 해` 전 구간에서 `|f| <= tol` 인 **마지막 교차점**을 찾는다.

        ⚠ 첫 교차를 찾으면 안 된다. 인력이 있는 포텐셜의 **힘은 비단조**다 —
          LJ 는 최솟점 `2^{1/6} sigma` 에서 힘이 정확히 0 이고 그 바깥에서 다시 커진다.
          첫 교차를 찾으면 "시작점에서 이미 무시할 수준"이라 판단하고 `r_start` 를
          돌려준다. 2026-07-28 에 실제로 그렇게 구현했고 테스트가 잡았다.
          (Morse 의 힘도 같은 이유로 비단조다.)
        """
        g = lambda r: abs(f(r)) - tol
        grid = np.geomspace(r_start, r_max, 3000)
        over = np.nonzero(np.array([g(r) for r in grid]) > 0.0)[0]
        if over.size == 0:
            return r_start                      # 전 구간에서 무시할 수준
        last = int(over[-1])
        if last == grid.size - 1:
            return None                         # r_max 까지 가도 허용치에 못 든다
        return float(brentq(g, grid[last], grid[last + 1], xtol=1e-13, rtol=1e-14))

    r_u = solve(beta_u, beta_u_tol)
    r_f = solve(beta_f_sigma, beta_f_sigma_tol)

    warnings: list[str] = []
    if r_u is None or r_f is None:
        warnings.append(
            f"r/sigma = {r_max:g} 까지 가도 허용치를 만족하지 못한다 — "
            f"장거리 상호작용이다. Ewald 또는 명시적 절단 타협이 필요하다."
        )
        r_cut = r_max
        criterion = "r_max_reached"
    elif r_u >= r_f:
        r_cut, criterion = r_u, "potential_tolerance"
    else:
        r_cut, criterion = r_f, "force_tolerance"

    return CutoffProposal(
        potential=potential,
        r_cut_over_sigma=r_cut,
        criterion=criterion,
        beta_U_at_cut=abs(beta_u(r_cut)),
        beta_F_sigma_at_cut=abs(beta_f_sigma(r_cut)),
        min_box_over_sigma=2.0 * r_cut,
        exact=False,
        rationale=(
            f"|beta U| <= {beta_u_tol:g} (해 r/sigma = "
            f"{'n/a' if r_u is None else format(r_u, '.4g')}) 와 "
            f"|beta F sigma| <= {beta_f_sigma_tol:g} (해 "
            f"{'n/a' if r_f is None else format(r_f, '.4g')}) 중 "
            f"{'포텐셜' if criterion == 'potential_tolerance' else '힘'} 제약이 지배."
        ),
        warnings=warnings,
        hoomd_note=hoomd_note,
    )


# =============================================================================
# 포텐셜별 래퍼
# =============================================================================
def lj_cutoff(beta_epsilon: float, **kw) -> CutoffProposal:
    """인력 있는 절단 LJ.  U = 4eps[(s/r)^12 - (s/r)^6],  거리는 sigma_LJ 단위.

    관습값 `2.5 sigma` 와 비교해 보고한다 — 관습을 쓸지 말지는 근거를 보고 정한다.
    """
    bu = lambda r: beta_epsilon * 4.0 * (r ** -12 - r ** -6)
    bf = lambda r: beta_epsilon * 24.0 * (2.0 * r ** -13 - r ** -7)
    p = cutoff_from_tolerance(
        potential="LJ (truncated)", beta_u=bu, beta_f_sigma=bf,
        r_start=WCA_RCUT_OVER_SIGMA_LJ, r_max=50.0,
        hoomd_note='md.pair.LJ(..., mode="shift") 또는 md.pair.ForceShiftedLJ '
                   '(힘 불연속까지 없애려면 후자)',
        **kw)
    p.rationale += (
        f"  [관습 r_cut = 2.5 sigma 에서는 beta U = {abs(bu(2.5)):.3e}, "
        f"beta F sigma = {abs(bf(2.5)):.3e}]"
    )
    return p


def yukawa_cutoff(beta_epsilon: float, kappa_sigma: float, **kw) -> CutoffProposal:
    """스크리닝 Coulomb.  U = eps * exp(-kappa r)/r  (HOOMD md.pair.Yukawa 규약).

    거리는 sigma 단위, `kappa_sigma` = kappa*sigma.
    ★ `r_cut` 이 스크리닝 길이에 지배되므로 `kappa_sigma` 가 작으면 (약한 스크리닝)
      사거리가 폭발한다 — 계면 콜로이드 카드가 겪은 문제가 이것이다.
    """
    ks = kappa_sigma
    bu = lambda r: beta_epsilon * math.exp(-ks * r) / r
    bf = lambda r: beta_epsilon * math.exp(-ks * r) * (1.0 / r ** 2 + ks / r)
    return cutoff_from_tolerance(
        potential=f"Yukawa (kappa*sigma={ks:g})", beta_u=bu, beta_f_sigma=bf,
        r_start=1.0, r_max=200.0,
        hoomd_note="md.pair.Yukawa(nlist, default_r_cut=r_cut); "
                   'params = dict(epsilon=eps, kappa=kappa).  mode="shift" 권장',
        **kw)


def morse_cutoff(beta_D0: float, alpha_sigma: float, r0_over_sigma: float = 1.0,
                 **kw) -> CutoffProposal:
    """Morse 인력.  U = D0[exp(-2a(r-r0)) - 2exp(-a(r-r0))]  (HOOMD 규약).

    응집 콜로이드(도메인 D)의 기본 인력.
    """
    a, r0 = alpha_sigma, r0_over_sigma
    def bu(r: float) -> float:
        e = math.exp(-a * (r - r0))
        return beta_D0 * (e * e - 2.0 * e)
    def bf(r: float) -> float:
        e = math.exp(-a * (r - r0))
        return beta_D0 * 2.0 * a * (e * e - e)
    return cutoff_from_tolerance(
        potential=f"Morse (alpha*sigma={a:g}, r0/sigma={r0:g})",
        beta_u=bu, beta_f_sigma=bf, r_start=r0, r_max=100.0,
        hoomd_note="md.pair.Morse(nlist, default_r_cut=r_cut); "
                   'params = dict(D0=D0, alpha=alpha, r0=r0).  mode="shift" 권장',
        **kw)


# =============================================================================
# ★ 최종 제안 — 관습을 먼저 보되, 감당 가능하면 늘린다
# =============================================================================
def neighbors_per_particle(r_cut_over_sigma: float, phi: float, dim: int = 3) -> float:
    """`r_cut` 안에 들어오는 이웃 수 (입자당).

    3D:  n = 6 phi/(pi sigma^3)  =>  N_nb = (4pi/3) n r_cut^3 = 8 phi (r_cut/sigma)^3
    2D:  n = 4 phi/(pi sigma^2)  =>  N_nb = pi n r_cut^2      = 4 phi (r_cut/sigma)^2

    ★ **이것이 `r_cut` 비용의 전부다.** 희박계(`phi` 작음)에서는 `r_cut` 을 크게 잡아도
      이웃이 거의 없으므로 사실상 공짜다. 단일 입자면 0 이다.
    """
    if not 0.0 <= phi:
        raise ValueError("phi must be >= 0")
    if dim == 3:
        return 8.0 * phi * r_cut_over_sigma ** 3
    if dim == 2:
        return 4.0 * phi * r_cut_over_sigma ** 2
    raise ValueError(f"dim must be 2 or 3, got {dim}")


# Lambda(처리량 상수) 를 측정한 조건. 비용 환산의 기준선.
#   knowledge/wiki/findings/local-cpu-parallelism.md — WCA, phi=0.30, 3D
LAMBDA_BASELINE_RCUT = WCA_RCUT_OVER_SIGMA_LJ
LAMBDA_BASELINE_PHI = 0.30
DEFAULT_NEIGHBOR_BUDGET = 150.0     # 입자당 이웃 수 상한. 관습(phi=0.3, 2.5s)의 약 4배


def pair_cost_vs_lambda_baseline(r_cut_over_sigma: float, phi: float,
                                 dim: int = 3) -> float:
    """`estimators.THROUGHPUT_PARTICLE_STEPS_PER_S` 기준선 대비 쌍힘 비용 배수.

    ⚠ **처리량 상수는 WCA(`r_cut = 1.122 sigma`, `phi = 0.30`)에서 측정되었다.**
      그 조건의 이웃 수는 3.4 개뿐이다 — 매우 싼 구간이다.
      `r_cut` 을 늘리면 wall time 예측이 그만큼 낙관적으로 틀린다.
      이 배수를 `run_policy.yaml` 의 예산 검사에 곱해야 한다.
    """
    base = neighbors_per_particle(LAMBDA_BASELINE_RCUT, LAMBDA_BASELINE_PHI, dim)
    here = neighbors_per_particle(r_cut_over_sigma, phi, dim)
    return max(here, 1e-12) / max(base, 1e-12)


def propose_cutoff(
    tolerance_proposal: CutoffProposal,
    *,
    phi: float,
    n_particles: int,
    dim: int = 3,
    box_over_sigma: float | None = None,
    conventional: float = CONVENTIONAL_LJ_RCUT,
    neighbor_budget: float = DEFAULT_NEIGHBOR_BUDGET,
) -> CutoffProposal:
    """★ 최종 `r_cut` 결정. 사용자 규칙(2026-07-28)을 구현한다:

    > *"관습 값을 먼저 고려하되, 시스템이 다일루트(싱글파티클 또는 적은 입자)하거나
    >   롱 레인지 인터랙션의 경우는 더 길게 고려해도 될 것 같네."*

    구현 논리:

    1. **정의상 고정이면 그대로 쓴다** (WCA). 논의 종료.
    2. 허용치 기반 `r_cut` 이 관습값보다 **작으면** 그것을 쓴다 — 공짜 이득.
    3. 크면 **이웃 수로 감당 가능한지** 본다. `N_nb = 8 phi (r_cut/sigma)^3`.
       - 감당 가능 (`<= neighbor_budget`) → **허용치 기반을 쓴다.** 정확도 이득이 공짜에 가깝다.
         희박계·단일입자·장거리 상호작용이 전부 여기 해당한다.
       - 감당 불가 → 예산에 맞는 최대 `r_cut` 으로 줄이고, **남은 오차를
         알려진 한계로 기록한다** (조용히 넘어가지 않는다).
    4. `r_cut <= L/2` 는 절대 제약이므로 마지막에 적용한다. 위반하면 경고를 남긴다 —
       이 경우 자동으로 해결할 방법이 없다 (계면 콜로이드 카드 §7).
    """
    p = tolerance_proposal
    if p.exact:
        p.rationale += "  [정의상 고정 — 관습·비용 논의가 개입하지 않는다]"
        return p

    r_tol = p.r_cut_over_sigma
    nb_tol = neighbors_per_particle(r_tol, phi, dim)
    nb_conv = neighbors_per_particle(conventional, phi, dim)
    notes: list[str] = list(p.warnings)

    single = n_particles <= 1
    if single:
        chosen, why = r_tol, "single_particle_free"
        notes.append("입자가 1개다 — 쌍 상호작용이 없으므로 r_cut 이 동역학에 무관하다.")
    elif r_tol <= conventional:
        chosen, why = r_tol, "tolerance_below_convention"
    elif nb_tol <= neighbor_budget:
        chosen, why = r_tol, "dilute_or_longrange_affordable"
        notes.append(
            f"관습 {conventional:g} sigma 보다 크지만 이웃이 {nb_tol:.1f}개뿐이라 "
            f"(예산 {neighbor_budget:g}) 감당 가능 — 정확도를 취한다. phi={phi:g}."
        )
    else:
        # 예산에 맞는 최대 r_cut 을 역산
        if dim == 3:
            chosen = (neighbor_budget / (8.0 * phi)) ** (1.0 / 3.0)
        else:
            chosen = (neighbor_budget / (4.0 * phi)) ** 0.5
        chosen = max(chosen, conventional)
        why = "neighbor_budget_limited"
        notes.append(
            f"허용치는 r_cut = {r_tol:.3g} sigma 를 요구하지만 이웃이 {nb_tol:.0f}개로 "
            f"예산 {neighbor_budget:g}을 넘는다. {chosen:.3g} sigma 로 줄였다 — "
            f"★ 남은 오차: beta U = {abs(p.beta_U_at_cut):.2e} 목표였으나 실제로는 더 크다. "
            f"이것은 **알려진 한계**이며 S7 에서 절단 민감도를 확인해야 한다."
        )

    out = CutoffProposal(
        potential=p.potential,
        r_cut_over_sigma=chosen,
        criterion=why,
        beta_U_at_cut=p.beta_U_at_cut if chosen >= r_tol else float("nan"),
        beta_F_sigma_at_cut=p.beta_F_sigma_at_cut if chosen >= r_tol else float("nan"),
        min_box_over_sigma=2.0 * chosen,
        exact=False,
        tolerance_preset=p.tolerance_preset,
        rationale=(
            f"{p.rationale}  || 최종 결정: {why}. "
            f"이웃 수 {neighbors_per_particle(chosen, phi, dim):.1f}개 "
            f"(관습 {conventional:g} sigma 에서는 {nb_conv:.1f}개). "
            f"Lambda 기준선 대비 쌍힘 비용 "
            f"{pair_cost_vs_lambda_baseline(chosen, phi, dim):.2f}x."
        ),
        warnings=notes,
        hoomd_note=p.hoomd_note,
    )

    if box_over_sigma is not None:
        ok, msg = out.box_ok(box_over_sigma)
        if not ok:
            out.warnings.append(msg)
    return out


# =============================================================================
# 이웃 리스트 버퍼 — 정확도가 아니라 성능
# =============================================================================
def neighbor_list_buffer(r_cut_over_sigma: float,
                         max_step_displacement_over_sigma: float,
                         rebuild_every: int = 20) -> float:
    """Cell 리스트 버퍼 제안.  `r_buff >= rebuild_every * (스텝당 최대 변위) * 2`.

    두 입자가 서로 마주 다가올 수 있으므로 2배.
    ⚠ 버퍼는 **정확도가 아니라 성능** 문제다 — 다만 너무 작으면 이웃을 놓쳐
      힘이 조용히 틀린다. HOOMD 는 자동 재구축하지만 버퍼가 성능을 정한다.
    """
    if max_step_displacement_over_sigma <= 0:
        raise ValueError("변위는 > 0 이어야 한다")
    need = 2.0 * rebuild_every * max_step_displacement_over_sigma
    # 실무 범위로 클램프: 너무 크면 이웃 수가 (1+buff/r_cut)^3 로 늘어난다
    return float(min(max(need, 0.1 * r_cut_over_sigma), 0.5 * r_cut_over_sigma))
