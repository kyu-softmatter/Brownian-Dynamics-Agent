"""S4 — non-dimensionalization · dimensionless groups · choosing `dt`. 0 LLM lines.

## The two conventions this module enforces

**① The reference scales are owned by the (system × target dynamics) card.**
Impose a universal convention and in a trap system `Δt` becomes 240,000x the
relaxation time -- which is more dangerous because it does not diverge. It skips
the entire relaxation process and still produces plausible-looking numbers.
Measured: `runs/2026-07-28_trap-2d-5um_2dfb9d/08_conclusion.md` §3.

**② The `dt` gate is displacement-based. `dt/τ_D` is only recorded.**
A fixed `dt/τ_D` gate actually rejects 2 of the 3 runs that made it into papers.
Basis: `knowledge/wiki/findings/dt-gate-should-be-displacement-based.md`.

## Why the computation is done in SI

Writing the constraints directly in reduced units means asking "which time is the
`*` in `Δt*`" every single time, and in a system where `τ_D` and `τ_trap` differ by
240,000x that mistake passes quietly.
**Everything is computed in SI and divided by the card's time scale at the end.**
That makes a unit mix-up structurally impossible.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

from .policy import Policy, load_policy
from .spec import SystemSpec, derive, reference_length_si
from .units import Scales, scales_brownian, scales_harmonic_trap

# =============================================================================
# card → reference scales
# =============================================================================
#  The suffix of a card's name is the target dynamics. The same system gets
#  different scales depending on what is being looked at.
#  Cards: knowledge/wiki/systems/_index.md
CARD_SCALE_RULES: dict[str, str] = {
    "passive-sphere--harmonic-trap": "harmonic_trap",
    "passive-sphere--transport": "brownian",
    "passive-sphere--equilibrium-structure": "brownian",
    "interfacial-colloid--transport": "brownian",
    "interfacial-colloid--equilibrium-structure": "brownian",
    "abp--dense-collective": "active_run_length",
}


def scale_rule_for(card: str, spec: SystemSpec | None = None) -> str:
    """The scale rule the card requires. For an unregistered card it **refuses to
    improvise a non-dimensionalization.**"""
    if card in CARD_SCALE_RULES:
        return CARD_SCALE_RULES[card]
    raise KeyError(
        f"no scale rule is registered for card {card!r}. Improvised "
        f"non-dimensionalization is forbidden —\n"
        f"  create a `status: draft` card from "
        f"knowledge/wiki/systems/_TEMPLATE.md first\n"
        f"  and register it in CARD_SCALE_RULES "
        f"(CLAUDE.md §non-dimensionalization convention).\n"
        f"  registered cards: {sorted(CARD_SCALE_RULES)}")


def scales_for(spec: SystemSpec, derived: dict | None = None) -> Scales:
    """The (length, energy, time) scales the card fixes."""
    d = derived if derived is not None else derive(spec)
    rule = scale_rule_for(spec.card, spec)
    T = spec.medium.T_si.si
    gamma = d["gamma_si"]

    if rule == "harmonic_trap":
        if "k_si" not in d:
            raise ValueError("a trap card with no harmonic-trap k_si")
        return scales_harmonic_trap(k_si=d["k_si"], T_si=T, gamma_si=gamma)
    if rule == "brownian":
        return scales_brownian(sigma_si=d["sigma_si"], T_si=T, gamma_si=gamma)
    if rule == "active_run_length":
        raise NotImplementedError(
            "the ABP card's scales (run length ℓ, τ_r = 1/D_r) are not implemented "
            "yet — implement them when the card is promoted to usable")
    raise KeyError(f"unknown scale rule {rule!r}")


# =============================================================================
# choosing dt — displacement-based
# =============================================================================
@dataclass
class DtConstraint:
    """One constraint. `dt_si_max=None` means 'not applicable to this system'."""

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
        rows = ["| constraint | active | `Δt` cap [s] | `Δt*` cap | basis |",
                "|---|---|---|---|---|"]
        for c in self.constraints:
            if c.dt_si_max is None:
                lim_si, lim_star = "—", "—"
            else:
                lim_si = f"`{c.dt_si_max:.4g}`"
                lim_star = f"`{c.dt_si_max / self.dt_si * self.dt_star:.4g}`"
            mark = "✅" if c.applies else ("— off" if not c.active else "— n/a")
            note = c.basis if c.active else (c.off_reason or c.basis)
            star = " **←dominant**" if c.name == self.dominant else ""
            rows.append(f"| `{c.name}`{star} | {mark} | {lim_si} | {lim_star} | {note} |")
        return "\n".join(rows)


# =============================================================================
# The gate formulas — **unit-agnostic.** The canonical copy lives in `bdbot.dt`
# =============================================================================
#  ★ When a script re-writes the thresholds and formulas, the two copies drift
#    apart. On 2026-07-28 `scripts/chain_bend.py` actually had `0.03` and `0.005`
#    nailed in by hand, and the policy file's thresholds could change without the
#    script following.
#    Call it in SI and you get an SI cap; call it in reduced units and you get a
#    reduced cap -- only the arguments have to match.
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
#  The accuracy gates (thermal, force, active) and the **stability gate** protect
#  different things: violate accuracy and you still get an answer; violate
#  stability and there is no answer. They cannot be merged.
from bdbot.dt import (dt_max_active, dt_max_force, dt_max_stability,  # noqa: E402
                      dt_max_thermal)


def choose_dt(spec: SystemSpec, *, derived: dict | None = None,
              scales: Scales | None = None, policy: Policy | None = None,
              max_force_si: float | None = None,
              target_em_bias: float | None = None) -> DtChoice:
    """Take the **minimum** of the constraints and record which one dominated.

    Args:
        max_force_si: the maximum force [N] **actually computed** on the initial
            placement. Estimating it is forbidden -- without it the force
            constraint stays `n/a` and that fact shows up in the table.
        target_em_bias: the target Euler–Maruyama bias (relative) for a harmonic
            trap. Given it, the `Δt ≤ 2b/(1+b) · τ_trap` constraint turns on.

    ⚠ A displacement gate only means something **when there is something to overlap
      with or be bonded to.** With no partner, "the displacement cap that prevents
      overlap" is checking a problem that does not exist.
      ★ Write that condition as `bool(spec.pair)` and the gate turns quietly off in
      a bond-only system.
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

    # --- 1. thermal displacement: sqrt(2 D0 dt) <= delta_th * sigma (per component)
    delta_th = float(ts.get("max_thermal_displacement_sigma", 0.03))
    cs.append(DtConstraint(
        name="thermal_displacement",
        dt_si_max=dt_max_thermal(delta_th, sigma, D0),
        active=has_partner,
        basis=f"√(2 D₀ Δt) ≤ {delta_th:g} σ (per component, lab convention)",
        off_reason="nothing to overlap with and nothing bonded — a displacement "
                   "cap is meaningless"))

    # --- 2. force displacement: (F_max/gamma) dt <= delta_F * sigma ---
    delta_F = float(ts.get("max_force_displacement_sigma", 0.005))
    cs.append(DtConstraint(
        name="force_displacement",
        dt_si_max=dt_max_force(delta_F, sigma, gamma, max_force_si),
        active=has_partner,
        # ★ "I measured it and got 0" and "I have not measured it yet" are never
        #   written as the same sentence. The first is physics (a stationary point,
        #   so this gate is toothless); the second is a procedure violation.
        basis=("measured max|F| = 0 — a stationary point, so this gate is "
               "**toothless** (a stability gate is required)"
               if max_force_si == 0.0 else
               f"max|F|Δt/γ ≤ {delta_F:g} σ, measured max|F| = "
               f"{max_force_si:.4g} N"
               if max_force_si else
               "max|F| has not been computed yet — **estimating is forbidden** "
               "(§5.4)"),
        off_reason="no pair or bond interaction"))

    # --- 3. stiffness stability: dt <= s * 2 gamma / lambda_max ---
    #  ★ Not an accuracy gate. Violating it does not make the answer wrong -- **the
    #    integration diverges** -- and yet `check_finite` passes (a bond length of
    #    1 → 1.4e7 is still finite).
    #    A straight chain is a stationary point with max|F| = 0, so the force gate
    #    cannot stop this divergence.
    #    Measured table and the basis for the safety factor:
    #    knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
    safety = float(ts.get("stability_safety_factor", 0.2))
    lam = d.get("lambda_max_si")
    cs.append(DtConstraint(
        name="stiff_stability",
        dt_si_max=dt_max_stability(safety, gamma, lam),
        active=lam is not None,
        basis=(f"Δt ≤ {safety:g} · 2γ/λ_max, λ_max = 4k_bond + 16k_angle/b² "
               f"= {lam:.4g} N/m" if lam else "no bond or angle stiffness"),
        off_reason="no bonds or angles, so no stiff unstable mode"))

    # --- 4. shortest relaxation time ---
    #  ⚠ The bond and angle relaxation times (`tau_bond_si`, `tau_angle_si`) are
    #    **deliberately not put here.** `zeta = 0.01` comes from "to observe a
    #    relaxation you have to slice it that finely", and a stiff bond is not
    #    something being observed -- it is a constraint (a degree of freedom whose
    #    fluctuation is being ignored on purpose).
    #    At `k_bond* = 1e6`, `zeta·γ/k_bond = 1e-8`, which is 10x stricter than the
    #    measured-calibrated stability gate `9.6e-8` and multiplies the step count
    #    by 10 -- a cost with no basis.
    #    A bond's divergence is stopped by `stiff_stability`, and the bond lengths
    #    are confirmed afterwards by `guards.check_bond_lengths`.
    zeta = float(ts.get("relaxation_safety_factor", 0.01))
    relax_times = {k: d[k] for k in ("tau_trap_si",) if k in d}
    tau_min = min(relax_times.values()) if relax_times else None
    cs.append(DtConstraint(
        name="relaxation_time",
        dt_si_max=(zeta * tau_min if tau_min else None),
        active=tau_min is not None,
        basis=(f"Δt ≤ {zeta:g} · min({', '.join(relax_times)}) = "
               f"{zeta:g} × {tau_min:.4g} s" if tau_min
               else "no relaxation-time scale"),
        off_reason="no confinement or active relaxation time"))

    # --- 5. active displacement ---
    v0 = next((e.params["v0_si"].si for e in spec.external
               if e.kind == "active" and "v0_si" in e.params), None)
    delta_a = float(ts.get("max_active_displacement_sigma", 0.01))
    cs.append(DtConstraint(
        name="active_displacement",
        dt_si_max=dt_max_active(delta_a, sigma, v0),
        active=v0 is not None,
        basis=(f"v₀Δt ≤ {delta_a:g} σ" if v0 else "no active drive"),
        off_reason="no active drive or flow"))

    # --- 6. integrator-bias target for a harmonic trap ---
    #  ★ This constraint comes from an accuracy target. The EM bias has a known
    #    analytic form, so fixing "how wrong is acceptable" fixes Δt.
    if target_em_bias is not None and "tau_trap_si" in d:
        from .estimators import dt_star_for_trap_bias
        cs.append(DtConstraint(
            name="em_bias_target",
            dt_si_max=dt_star_for_trap_bias(target_em_bias) * d["tau_trap_si"],
            active=True,
            basis=f"⟨x*²⟩ bias ≤ {target_em_bias:.3%} "
                  f"(1/(1−Δt*/2)−1, estimators)"))

    applicable = [c for c in cs if c.applies]
    if not applicable:
        raise ValueError(
            "not a single active constraint — there is no basis for choosing dt. "
            "Confirm this really is a system with no pair or bond interaction, no "
            "stiffness and no relaxation time, and if so state target_em_bias")

    winner = min(applicable, key=lambda c: c.dt_si_max)
    dt_si = winner.dt_si_max
    #  `hard_floor` applies to the **accuracy constraints only.**
    #  ★ Rejecting on the floor when stability dominates blocks a run that actually
    #    works: at `k_bond* = 1e6` the gate is `9.6e-8` while the floor is `1e-7`.
    #    The finding's bisection measurement showed that system stable out to
    #    `Δt* = 1e-6` and completing 4000 steps at `9.6e-8` -- so the floor is what
    #    is wrong. And a stability cap is **not negotiable**: it is not a demand to
    #    lower it but a demand to lower `k_bond` (`k_bond* = C·κ(N)*`,
    #    finding §cost implications).
    #    Instead the fact that it is below the floor is recorded so the report
    #    shows it.
    floor = float(ts.get("hard_floor", 1e-7)) * sc.time_si
    below_floor = dt_si < floor
    if below_floor and winner.name != "stiff_stability":
        raise ValueError(
            f"Δt* = {dt_si / sc.time_si:.3g} is below hard_floor — "
            f"reconsider the modelling itself (dominant constraint: {winner.name})")

    return DtChoice(
        dt_si=dt_si, dt_star=dt_si / sc.time_si, dominant=winner.name,
        constraints=cs,
        # For the record — not a gate. Used only when comparing to another paper.
        logged={
            "dt_over_tau_D": dt_si / d["tau_D_si"],
            "thermal_displacement_over_sigma": math.sqrt(2 * D0 * dt_si) / sigma,
            **({"dt_over_tau_trap": dt_si / d["tau_trap_si"]}
               if "tau_trap_si" in d else {}),
            # Stability margin. The linear limit is `dt/τ_stiff = 2`, so a value
            # near 2 is on the edge of divergence. When `stiff_stability`
            # dominates this comes out exactly 2·safety.
            **({"dt_over_tau_stiff": dt_si / d["tau_stiff_si"]}
               if "tau_stiff_si" in d else {}),
            # Stability demanded below the accuracy floor — the cost lever is k_bond
            **({"dt_star_below_hard_floor": dt_si / sc.time_si}
               if below_floor else {}),
        })


# =============================================================================
# dimensionless groups
# =============================================================================
def groups(spec: SystemSpec, derived: dict | None = None,
           scales: Scales | None = None) -> dict[str, float]:
    """Every dimensionless group that can be computed (master_plan §5.3).

    Anything that **cannot** be computed is left out -- filling it with `None` or
    `0` reads in the report as "this group is 0".
    """
    d = derived if derived is not None else derive(spec)
    sc = scales if scales is not None else scales_for(spec, d)
    out: dict[str, float] = {}

    sigma, kT = d["sigma_si"], d["kT_si"]
    out["sigma_over_ref_length"] = sigma / sc.length_si
    out["tau_D_over_ref_time"] = d["tau_D_si"] / sc.time_si

    if "k_si" in d:
        out["k_star_sigma"] = d["k_si"] * sigma**2 / kT          # kσ²/kT
        out["k_star"] = d["k_si"] * sc.length_si**2 / kT         # card units (=1)
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
        out["phi"] = phi                # meaningless with no pair interaction (spec.py)

    return out


# =============================================================================
# ReducedSpec
# =============================================================================
@dataclass
class ReducedSpec:
    """The reduced specification + the back-conversion factors. These values go
    straight into HOOMD."""

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
        """Reduced → SI back-conversion factors. The report's conversion table
        uses these."""
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
    """SI spec → reduced spec. Without `dt_star`, `choose_dt` picks it."""
    d = derive(spec)
    sc = scales_for(spec, d)

    if dt_star is None:
        if spec.numerics.dt_star is not None:
            dt_star, dominant = spec.numerics.dt_star.si, "spec(stated)"
            logged = {"dt_over_tau_D": dt_star * sc.time_si / d["tau_D_si"]}
        else:
            ch = choose_dt(spec, derived=d, scales=sc, policy=policy,
                           max_force_si=max_force_si,
                           target_em_bias=target_em_bias)
            dt_star, dominant, logged = ch.dt_star, ch.dominant, ch.logged
    else:
        dominant, logged = "caller(stated)", {
            "dt_over_tau_D": dt_star * sc.time_si / d["tau_D_si"]}

    box_si = spec.box_lengths_si(sc.length_si)
    if box_si is None:
        raise ValueError("the box is unknown, so it cannot be non-dimensionalized")

    def steps(q, default_tau: float) -> int:
        tau = q.si if q is not None else default_tau
        return int(round(tau / dt_star))

    t = spec.timing
    return ReducedSpec(
        card=spec.card, scales=sc, dim=spec.geometry.d,
        n_particles=int(spec.primary.n_simulated.si),
        box_star=[x / sc.length_si for x in box_si],
        kT_star=d["kT_si"] / sc.energy_si,                      # 1 by definition
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
# round-trip verification — the master_plan §S4 gate
# =============================================================================
def roundtrip_errors(spec: SystemSpec, reduced: ReducedSpec | None = None
                     ) -> dict[str, float]:
    """The `to_reduced → to_si` relative error. The gate is `< 1e-12`.

    Non-dimensionalization is one division, so the error cannot be large -- **a
    large error is a signal that a convention was broken** (dividing the time scale
    by `τ_D` and converting it back with `τ_trap`, say).
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
    """The conversion table for `04_nondim.md`:
    quantity | SI | reduced | back-conversion factor."""
    d = derive(spec)
    r = reduced if reduced is not None else reduce_spec(spec)
    sc = r.scales
    rows = ["| quantity | SI | reduced | back-conversion factor |",
            "|---|---|---|---|"]

    def row(name, si, star, coeff, unit=""):
        rows.append(f"| {name} | `{si:.6g}`{f' {unit}' if unit else ''} "
                    f"| `{star:.6g}` | `{coeff:.6g}` |")

    row("length scale", sc.length_si, 1.0, sc.length_si, "m")
    row("energy scale `kT`", sc.energy_si, 1.0, sc.energy_si, "J")
    row("time scale", sc.time_si, 1.0, sc.time_si, "s")
    row("`σ`", d["sigma_si"], r.sigma_star, sc.length_si, "m")
    row("`γ`", d["gamma_si"], r.gamma_star,
        sc.energy_si * sc.time_si / sc.length_si**2, "kg/s")
    row("`D₀`", d["D0_si"], r.D_star, sc.diffusivity_si, "m²/s")
    if r.k_star is not None:
        row("`k`", d["k_si"], r.k_star, sc.stiffness_si, "N/m")
    row("`Δt`", r.dt_star * sc.time_si, r.dt_star, sc.time_si, "s")
    rows.append(f"| **scale origin** | {sc.origin} | | |")
    return "\n".join(rows)
