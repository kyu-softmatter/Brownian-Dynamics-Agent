"""S4 — 무차원화 · 무차원수 · `dt` 선택. LLM 0줄.

## 이 모듈이 지키는 두 규약

**① 기준 척도는 (계 × 목적동역학) 카드가 소유한다.**
보편 규약을 강요하면 트랩 계에서 `Δt` 가 완화시간의 24만 배가 된다 — 발산하지
않으므로 더 위험하다. 완화 과정을 통째로 건너뛰고도 그럴듯한 숫자가 나온다.
실측: `runs/2026-07-28_trap-2d-5um_2dfb9d/08_conclusion.md` §3.

**② `dt` 게이트는 변위 기준이다. `dt/τ_D` 는 기록만 한다.**
`dt/τ_D` 고정 게이트는 실제로 논문까지 나온 런 3건 중 2건을 기각한다.
근거: `knowledge/wiki/findings/dt-gate-should-be-displacement-based.md`.

## 왜 SI 에서 계산하는가

제약을 무차원 단위로 바로 쓰면 "`Δt*` 의 `*` 가 어느 시간인가"를 매번 따져야 하고,
`τ_D` 와 `τ_trap` 이 24만 배 차이 나는 계에서는 그 실수가 조용히 통과한다.
**전부 SI 로 계산하고 마지막에 카드의 시간 척도로 나눈다.** 그러면 단위 혼동이
구조적으로 불가능해진다.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .policy import Policy, load_policy
from .spec import SystemSpec, derive, reference_length_si
from .units import Scales, scales_brownian, scales_harmonic_trap

# =============================================================================
# 카드 → 기준 척도
# =============================================================================
#  카드 이름의 접미사가 목적동역학이다. 같은 계라도 무엇을 보느냐에 따라 다르다.
#  카드: knowledge/wiki/systems/_index.md
CARD_SCALE_RULES: dict[str, str] = {
    "passive-sphere--harmonic-trap": "harmonic_trap",
    "passive-sphere--transport": "brownian",
    "passive-sphere--equilibrium-structure": "brownian",
    "interfacial-colloid--transport": "brownian",
    "interfacial-colloid--equilibrium-structure": "brownian",
    "abp--dense-collective": "active_run_length",
}


def scale_rule_for(card: str, spec: SystemSpec | None = None) -> str:
    """카드가 요구하는 척도 규칙. 등록되지 않은 카드면 **즉흥 무차원화를 거부한다.**"""
    if card in CARD_SCALE_RULES:
        return CARD_SCALE_RULES[card]
    raise KeyError(
        f"카드 {card!r} 의 척도 규칙이 등록되지 않았다. 즉흥 무차원화는 금지다 —\n"
        f"  knowledge/wiki/systems/_TEMPLATE.md 로 `status: draft` 카드를 먼저 만들고\n"
        f"  CARD_SCALE_RULES 에 등록할 것 (CLAUDE.md §무차원화 규약).\n"
        f"  등록된 카드: {sorted(CARD_SCALE_RULES)}")


def scales_for(spec: SystemSpec, derived: dict | None = None) -> Scales:
    """카드가 정한 (길이, 에너지, 시간) 척도."""
    d = derived if derived is not None else derive(spec)
    rule = scale_rule_for(spec.card, spec)
    T = spec.medium.T_si.si
    gamma = d["gamma_si"]

    if rule == "harmonic_trap":
        if "k_si" not in d:
            raise ValueError("트랩 카드인데 조화 트랩의 k_si 가 없다")
        return scales_harmonic_trap(k_si=d["k_si"], T_si=T, gamma_si=gamma)
    if rule == "brownian":
        return scales_brownian(sigma_si=d["sigma_si"], T_si=T, gamma_si=gamma)
    if rule == "active_run_length":
        raise NotImplementedError(
            "ABP 카드의 척도 (런 길이 ℓ, τ_r = 1/D_r) 는 아직 구현되지 않았다 — "
            "카드가 usable 로 승격될 때 함께 구현할 것")
    raise KeyError(f"알 수 없는 척도 규칙 {rule!r}")


# =============================================================================
# dt 선택 — 변위 기준
# =============================================================================
@dataclass
class DtConstraint:
    """제약 1개. `dt_si_max=None` 은 '이 계에 해당 없음'이다."""

    name: str
    dt_si_max: float | None
    active: bool
    basis: str
    off_reason: str = ""

    @property
    def applies(self) -> bool:
        return self.active and self.dt_si_max is not None


@dataclass
class DtChoice:
    dt_si: float
    dt_star: float
    dominant: str
    constraints: list[DtConstraint]
    logged: dict[str, float] = field(default_factory=dict)

    def table(self) -> str:
        rows = ["| 제약 | 활성 | `Δt` 상한 [s] | `Δt*` 상한 | 근거 |",
                "|---|---|---|---|---|"]
        for c in self.constraints:
            if c.dt_si_max is None:
                lim_si, lim_star = "—", "—"
            else:
                lim_si = f"`{c.dt_si_max:.4g}`"
                lim_star = f"`{c.dt_si_max / self.dt_si * self.dt_star:.4g}`"
            mark = "✅" if c.applies else ("— off" if not c.active else "— n/a")
            note = c.basis if c.active else (c.off_reason or c.basis)
            star = " **←지배**" if c.name == self.dominant else ""
            rows.append(f"| `{c.name}`{star} | {mark} | {lim_si} | {lim_star} | {note} |")
        return "\n".join(rows)


# =============================================================================
# 게이트 식 — **단위 무관.** 정본은 `bdbot.dt` 다
# =============================================================================
#  ★ 스크립트가 임계값과 식을 다시 쓰면 두 곳이 갈라진다. 2026-07-28 에 실제로
#    `scripts/chain_bend.py` 가 `0.03`·`0.005` 를 손으로 박아 두고 있었고, 그 사이
#    정책 파일의 문턱이 바뀌어도 스크립트는 따라오지 않는 상태였다.
#    SI 로 부르면 SI 상한, 축약단위로 부르면 축약 상한이 나온다 — 인자만 맞추면 된다.
#
#  ★★ 2026-08-29: "this is the only definition" **was a false sentence.**
#     `bdbot.checks` held its own version of the same job, and worse, **on a
#     different criterion** -- displacement here, timescale ratio (`dt = 1e-2*tau`)
#     there. `.claude/rules/overdamped-stability.md` forbids the latter explicitly
#     ("displacement goes as the **square root** of Δt, so when the dimensionality
#     changes the same Δt/τ_D gives a different displacement"). The equations
#     therefore moved into `bdbot/dt.py` and this module only re-exports them.
#     No value changed -- the function bodies moved verbatim.
#     To see both criteria side by side: `bdbot.dt.compare_criteria`.
#
#  정확도 게이트(thermal·force·active)와 **안정성 게이트**(stability)는 다른 것을 지킨다:
#  정확도는 위반해도 답이 나오고, 안정성은 위반하면 답이 없다. 합칠 수 없다.
from bdbot.dt import (dt_max_active, dt_max_force, dt_max_stability,  # noqa: E402
                      dt_max_thermal)


def choose_dt(spec: SystemSpec, *, derived: dict | None = None,
              scales: Scales | None = None, policy: Policy | None = None,
              max_force_si: float | None = None,
              target_em_bias: float | None = None) -> DtChoice:
    """제약들의 **최소값**을 채택하고 어느 제약이 지배했는지 기록한다.

    Args:
        max_force_si: 초기배치에서 **실제로 계산한** 최대 힘 [N]. 추정 금지 —
            없으면 힘 제약은 `n/a` 로 남고 그 사실이 표에 드러난다.
        target_em_bias: 조화 트랩의 목표 Euler–Maruyama 편향 (상대값).
            주면 `Δt ≤ 2b/(1+b) · τ_trap` 제약이 켜진다.

    ⚠ 변위 게이트는 **겹칠 상대나 결합 상대가 있을 때만** 의미가 있다. 상대가 없으면
      "겹침을 막는 변위 상한"은 존재하지 않는 문제를 검사하는 것이다.
      ★ 그 조건을 `bool(spec.pair)` 로 쓰면 결합만 있는 계에서 게이트가 조용히 꺼진다.
    """
    d = derived if derived is not None else derive(spec)
    sc = scales if scales is not None else scales_for(spec, d)
    pol = policy if policy is not None else load_policy()
    ts = pol.timestep

    sigma = d["sigma_si"]
    D0 = d["D0_si"]
    gamma = d["gamma_si"]
    has_partner = spec.has_neighbor_interaction
    cs: list[DtConstraint] = []

    # --- 1. 열 변위:  sqrt(2 D0 dt) <= delta_th * sigma  (성분별) ---
    delta_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    cs.append(DtConstraint(
        name="thermal_displacement",
        dt_si_max=dt_max_thermal(delta_th, sigma, D0),
        active=has_partner,
        basis=f"√(2 D₀ Δt) ≤ {delta_th:g} σ (성분별, 랩 관행)",
        off_reason="겹칠 상대도 결합 상대도 없다 — 변위 상한이 무의미"))

    # --- 2. 힘 변위:  (F_max/gamma) dt <= delta_F * sigma ---
    delta_F = float(ts.get("max_force_displacement_sigma", 0.005))
    cs.append(DtConstraint(
        name="force_displacement",
        dt_si_max=dt_max_force(delta_F, sigma, gamma, max_force_si),
        active=has_partner,
        # ★ "재봤더니 0" 과 "아직 안 재봤다" 를 같은 문장으로 적지 않는다. 앞은 물리
        #   (정류점이라 이 게이트가 무력하다), 뒤는 절차 위반이다.
        basis=("실측 max|F| = 0 — 정류점이라 이 게이트는 **무력하다** "
               "(안정성 게이트가 있어야 한다)" if max_force_si == 0.0 else
               f"max|F|Δt/γ ≤ {delta_F:g} σ, 실측 max|F| = {max_force_si:.4g} N"
               if max_force_si else
               "max|F| 를 아직 계산하지 않았다 — **추정 금지** (§5.4)"),
        off_reason="쌍·결합 상호작용 없음"))

    # --- 3. 강성 안정성:  dt <= s * 2 gamma / lambda_max ---
    #  ★ 정확도 게이트가 아니다. 위반하면 답이 틀리는 게 아니라 **적분이 발산한다** —
    #    그런데 `check_finite` 는 통과한다 (결합길이 1 → 1.4e7 도 유한하다).
    #    직선 사슬은 max|F| = 0 인 정류점이라 힘 게이트가 이 발산을 못 막는다.
    #    실측 표와 안전계수 근거:
    #    knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
    safety = float(ts.get("stability_safety_factor", 0.2))
    lam = d.get("lambda_max_si")
    cs.append(DtConstraint(
        name="stiff_stability",
        dt_si_max=dt_max_stability(safety, gamma, lam),
        active=lam is not None,
        basis=(f"Δt ≤ {safety:g} · 2γ/λ_max, λ_max = 4k_bond + 16k_angle/b² "
               f"= {lam:.4g} N/m" if lam else "결합·각 강성 없음"),
        off_reason="결합·각이 없어 강성 불안정 모드가 없다"))

    # --- 4. 최단 완화시간 ---
    #  ⚠ 결합·각의 완화시간(`tau_bond_si`·`tau_angle_si`)은 **여기 넣지 않는다.**
    #    `zeta = 0.01` 은 "완화 과정을 관측하려면 그만큼 잘게 쪼개야 한다"에서 온 계수이고,
    #    강성 결합은 관측 대상이 아니라 구속이다 (요동을 일부러 무시하는 자유도).
    #    `k_bond* = 1e6` 에서 `zeta·γ/k_bond = 1e-8` 이라 실측 교정된 안정성 게이트
    #    `9.6e-8` 보다 10배 엄격해지고, 스텝 수가 10배로 늘어난다 — 근거 없는 비용이다.
    #    결합의 발산은 `stiff_stability` 가 막고, 결합길이는 `guards.check_bond_lengths`
    #    가 사후에 확인한다.
    zeta = float(ts.get("relaxation_safety_factor", 0.01))
    relax_times = {k: d[k] for k in ("tau_trap_si",) if k in d}
    tau_min = min(relax_times.values()) if relax_times else None
    cs.append(DtConstraint(
        name="relaxation_time",
        dt_si_max=(zeta * tau_min if tau_min else None),
        active=tau_min is not None,
        basis=(f"Δt ≤ {zeta:g} · min({', '.join(relax_times)}) = "
               f"{zeta:g} × {tau_min:.4g} s" if tau_min else "완화시간 척도 없음"),
        off_reason="구속·활성 완화시간이 없다"))

    # --- 5. 활성 변위 ---
    v0 = next((e.params["v0_si"].si for e in spec.external
               if e.kind == "active" and "v0_si" in e.params), None)
    delta_a = float(ts.get("max_active_displacement_sigma", 0.01))
    cs.append(DtConstraint(
        name="active_displacement",
        dt_si_max=dt_max_active(delta_a, sigma, v0),
        active=v0 is not None,
        basis=(f"v₀Δt ≤ {delta_a:g} σ" if v0 else "능동 구동 없음"),
        off_reason="능동·흐름 없음"))

    # --- 6. 조화 트랩의 적분기 편향 목표 ---
    #  ★ 이 제약은 정확도 목표에서 나온다. EM 편향은 알려진 해석식이므로
    #    "얼마나 틀려도 되는가"를 정하면 Δt 가 정해진다.
    if target_em_bias is not None and "tau_trap_si" in d:
        from .estimators import dt_star_for_trap_bias
        cs.append(DtConstraint(
            name="em_bias_target",
            dt_si_max=dt_star_for_trap_bias(target_em_bias) * d["tau_trap_si"],
            active=True,
            basis=f"⟨x*²⟩ 편향 ≤ {target_em_bias:.3%} "
                  f"(1/(1−Δt*/2)−1, estimators)"))

    applicable = [c for c in cs if c.applies]
    if not applicable:
        raise ValueError(
            "활성 제약이 하나도 없다 — dt 를 고를 근거가 없다. "
            "쌍·결합 상호작용도 강성도 완화시간도 없는 계인지 확인하고, 그렇다면 "
            "target_em_bias 를 명시할 것")

    winner = min(applicable, key=lambda c: c.dt_si_max)
    dt_si = winner.dt_si_max
    #  `hard_floor` 는 **정확도 제약에만** 적용한다.
    #  ★ 안정성이 지배할 때 floor 로 기각하면 실제로 도는 런을 막는다: `k_bond* = 1e6` 의
    #    게이트는 `9.6e-8` 인데 floor 는 `1e-7` 이다. finding 의 이분법 실측은 그 계가
    #    `Δt* = 1e-6` 까지 안정이고 `9.6e-8` 에서 4000 step 을 완주함을 보였다 —
    #    즉 floor 가 틀렸다. 그리고 안정성 상한은 **협상 대상이 아니다.** 낮추라는 요구가
    #    아니라 `k_bond` 를 낮추라는 요구다 (`k_bond* = C·κ(N)*`, finding §비용 함의).
    #    대신 floor 아래라는 사실을 기록해서 리포트에서 보이게 한다.
    floor = float(ts.get("hard_floor", 1e-7)) * sc.time_si
    below_floor = dt_si < floor
    if below_floor and winner.name != "stiff_stability":
        raise ValueError(
            f"Δt* = {dt_si / sc.time_si:.3g} 가 hard_floor 아래다 — "
            f"모델링 자체를 재검토할 것 (지배 제약: {winner.name})")

    return DtChoice(
        dt_si=dt_si, dt_star=dt_si / sc.time_si, dominant=winner.name,
        constraints=cs,
        # 기록용 — 게이트 아님. 다른 논문과 비교할 때만 쓴다.
        logged={
            "dt_over_tau_D": dt_si / d["tau_D_si"],
            "thermal_displacement_over_sigma": math.sqrt(2 * D0 * dt_si) / sigma,
            **({"dt_over_tau_trap": dt_si / d["tau_trap_si"]}
               if "tau_trap_si" in d else {}),
            # 안정성 여유. 선형 한계가 `dt/τ_stiff = 2` 이므로 이 값이 2 에 가까우면
            # 발산 직전이다. `stiff_stability` 가 지배하면 정확히 2·safety 가 된다.
            **({"dt_over_tau_stiff": dt_si / d["tau_stiff_si"]}
               if "tau_stiff_si" in d else {}),
            # 안정성이 정확도 floor 아래를 요구했다 — 비용 지레는 k_bond 다
            **({"dt_star_below_hard_floor": dt_si / sc.time_si}
               if below_floor else {}),
        })


# =============================================================================
# 무차원수
# =============================================================================
def groups(spec: SystemSpec, derived: dict | None = None,
           scales: Scales | None = None) -> dict[str, float]:
    """계산 가능한 무차원수 전량 (master_plan §5.3).

    계산 **불가능한** 것은 넣지 않는다 — `None` 이나 `0` 으로 채우면 리포트에서
    "이 무차원수는 0 이다"로 읽힌다.
    """
    d = derived if derived is not None else derive(spec)
    sc = scales if scales is not None else scales_for(spec, d)
    out: dict[str, float] = {}

    sigma, kT = d["sigma_si"], d["kT_si"]
    out["sigma_over_ref_length"] = sigma / sc.length_si
    out["tau_D_over_ref_time"] = d["tau_D_si"] / sc.time_si

    if "k_si" in d:
        out["k_star_sigma"] = d["k_si"] * sigma**2 / kT          # kσ²/kT
        out["k_star"] = d["k_si"] * sc.length_si**2 / kT         # 카드 단위 (=1)
        out["l_trap_over_sigma"] = d["l_trap_si"] / sigma
        out["tau_D_over_tau_trap"] = d["tau_D_si"] / d["tau_trap_si"]

    if "tau_inertial_si" in d:
        out["tau_inertial_over_ref_time"] = d["tau_inertial_si"] / sc.time_si

    if spec.medium.rho_fluid_si is not None:
        v = sc.velocity_si
        out["reynolds"] = (spec.medium.rho_fluid_si.si * v
                           * spec.primary.radius_si.si / spec.medium.eta_si.si)

    from .spec import packing_fraction
    box = spec.box_lengths_si(reference_length_si(spec, d))
    phi = packing_fraction(spec, box)
    if phi is not None and spec.pair:
        out["phi"] = phi                # 쌍 상호작용 없으면 의미 없음 (spec.py 참조)

    return out


# =============================================================================
# ReducedSpec
# =============================================================================
@dataclass
class ReducedSpec:
    """무차원 명세 + 역변환 계수. HOOMD 에 그대로 들어가는 값들."""

    card: str
    scales: Scales
    dim: int
    n_particles: int
    box_star: list[float]
    kT_star: float
    gamma_star: float
    D_star: float
    sigma_star: float
    dt_star: float
    dt_dominant: str
    k_star: float | None
    equil_steps: int
    prod_steps: int
    sample_interval_steps: int
    groups: dict[str, float]
    logged: dict[str, float]

    @property
    def inverse(self) -> dict[str, float]:
        """무차원 → SI 역변환 계수. 리포트의 변환표가 이걸 쓴다."""
        sc = self.scales
        return {"length": sc.length_si, "energy": sc.energy_si, "time": sc.time_si,
                "force": sc.force_si, "stiffness": sc.stiffness_si,
                "velocity": sc.velocity_si, "diffusivity": sc.diffusivity_si,
                "rate": sc.rate_si, "area": sc.length_si**2}

    def to_si(self, value_star: float, kind: str) -> float:
        return self.scales.to_si(value_star, kind)


def reduce_spec(spec: SystemSpec, *, dt_star: float | None = None,
                policy: Policy | None = None,
                max_force_si: float | None = None,
                target_em_bias: float | None = None) -> ReducedSpec:
    """SI 명세 → 무차원 명세. `dt_star` 를 주지 않으면 `choose_dt` 가 고른다."""
    d = derive(spec)
    sc = scales_for(spec, d)

    if dt_star is None:
        if spec.numerics.dt_star is not None:
            dt_star, dominant = spec.numerics.dt_star.si, "spec(명시값)"
            logged = {"dt_over_tau_D": dt_star * sc.time_si / d["tau_D_si"]}
        else:
            ch = choose_dt(spec, derived=d, scales=sc, policy=policy,
                           max_force_si=max_force_si,
                           target_em_bias=target_em_bias)
            dt_star, dominant, logged = ch.dt_star, ch.dominant, ch.logged
    else:
        dominant, logged = "caller(명시값)", {
            "dt_over_tau_D": dt_star * sc.time_si / d["tau_D_si"]}

    box_si = spec.box_lengths_si(sc.length_si)
    if box_si is None:
        raise ValueError("박스를 알 수 없어 무차원화할 수 없다")

    def steps(q, default_tau: float) -> int:
        tau = q.si if q is not None else default_tau
        return int(round(tau / dt_star))

    t = spec.timing
    return ReducedSpec(
        card=spec.card, scales=sc, dim=spec.geometry.d,
        n_particles=int(spec.primary.n_simulated.si),
        box_star=[x / sc.length_si for x in box_si],
        kT_star=d["kT_si"] / sc.energy_si,                      # 정의상 1
        gamma_star=d["gamma_si"] * sc.length_si**2 / (sc.energy_si * sc.time_si),
        D_star=d["D0_si"] / sc.diffusivity_si,
        sigma_star=d["sigma_si"] / sc.length_si,
        dt_star=dt_star, dt_dominant=dominant,
        k_star=(d["k_si"] / sc.stiffness_si if "k_si" in d else None),
        equil_steps=steps(t.equil_in_tau, 10.0),
        prod_steps=steps(t.prod_in_tau, 40.0),
        sample_interval_steps=max(1, steps(t.sample_interval_in_tau, 2.0)),
        groups=groups(spec, d, sc), logged=logged)


# =============================================================================
# 왕복 검증 — master_plan §S4 게이트
# =============================================================================
def roundtrip_errors(spec: SystemSpec, reduced: ReducedSpec | None = None
                     ) -> dict[str, float]:
    """`to_reduced → to_si` 상대오차. 게이트는 `< 1e-12`.

    무차원화는 나눗셈 한 번이므로 오차가 클 수 없다 — **큰 오차는 규약이 어긋났다는
    신호다** (예: 시간 척도를 `τ_D` 로 나누고 `τ_trap` 으로 되돌리기).
    """
    d = derive(spec)
    r = reduced if reduced is not None else reduce_spec(spec)
    sc = r.scales

    pairs: dict[str, tuple[float, float]] = {
        "kT": (d["kT_si"], sc.to_si(r.kT_star, "energy")),
        "D0": (d["D0_si"], sc.to_si(r.D_star, "diffusivity")),
        "sigma": (d["sigma_si"], sc.to_si(r.sigma_star, "length")),
        "gamma": (d["gamma_si"],
                  r.gamma_star * sc.energy_si * sc.time_si / sc.length_si**2),
    }
    if r.k_star is not None:
        pairs["k"] = (d["k_si"], sc.to_si(r.k_star, "stiffness"))
    box_si = spec.box_lengths_si(sc.length_si)
    if box_si:
        pairs["box_x"] = (box_si[0], sc.to_si(r.box_star[0], "length"))

    return {k: abs(back - orig) / abs(orig) if orig else abs(back)
            for k, (orig, back) in pairs.items()}


def nondim_table(spec: SystemSpec, reduced: ReducedSpec | None = None) -> str:
    """`04_nondim.md` 의 변환표: 물리량 | SI | 무차원 | 역변환계수."""
    d = derive(spec)
    r = reduced if reduced is not None else reduce_spec(spec)
    sc = r.scales
    rows = ["| 물리량 | SI | 무차원 | 역변환 계수 |", "|---|---|---|---|"]

    def row(name, si, star, coeff, unit=""):
        rows.append(f"| {name} | `{si:.6g}`{f' {unit}' if unit else ''} "
                    f"| `{star:.6g}` | `{coeff:.6g}` |")

    row("길이 척도", sc.length_si, 1.0, sc.length_si, "m")
    row("에너지 척도 `kT`", sc.energy_si, 1.0, sc.energy_si, "J")
    row("시간 척도", sc.time_si, 1.0, sc.time_si, "s")
    row("`σ`", d["sigma_si"], r.sigma_star, sc.length_si, "m")
    row("`γ`", d["gamma_si"], r.gamma_star,
        sc.energy_si * sc.time_si / sc.length_si**2, "kg/s")
    row("`D₀`", d["D0_si"], r.D_star, sc.diffusivity_si, "m²/s")
    if r.k_star is not None:
        row("`k`", d["k_si"], r.k_star, sc.stiffness_si, "N/m")
    row("`Δt`", r.dt_star * sc.time_si, r.dt_star, sc.time_si, "s")
    rows.append(f"| **척도 출처** | {sc.origin} | | |")
    return "\n".join(rows)
