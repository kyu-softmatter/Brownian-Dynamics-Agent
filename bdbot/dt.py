"""`dt` selection -- **the single source of truth for the gate equations.**

Pure floats, `math` only, **unit-agnostic**: call it with SI and you get an SI
bound, call it with reduced units and you get a reduced bound. Only the arguments
have to be consistent.

    bdbot.checks     the pint-carrying wrappers        -> imports this
    simbot.nondim    `choose_dt` policy orchestration  -> imports this

**Why this module exists (measured 2026-08-29).** The `dt` rule existed in three
places with **two different criteria**:

    bdbot.checks.dt_from_gate       dt = 1e-2 * tau            timescale ratio
    simbot.nondim.dt_max_*          displacement per step      displacement
    campaigns/chain_bend.py:113     hand-inlined thresholds    displacement

`.claude/rules/overdamped-stability.md` settles which is right, and it is not the
first one: *"Do not use `Δt/τ_D` as the criterion: displacement goes as the
**square root** of `Δt`, so when the dimensionality changes the same `Δt/τ_D`
gives a different displacement."* The displacement gates are therefore the
primary criterion here and the ratio gate is kept as a **secondary, reported**
one -- the same rule file asks for exactly that: *"Also record what the rule you
did not adopt would have given. It is a cheap way to make a silent error
visible."* -> `compare_criteria()`.

Accuracy gates (thermal / force / active) and the **stability** gate protect
different things and cannot be merged: violate an accuracy gate and you still
get an answer, violate stability and there is no answer.
(`simbot/nondim.py` carried this note first; it is the reason the two kinds are
separate functions rather than one `min`.)
"""
from __future__ import annotations

# The historical `bdbot` gate: dt <= GATE * tau_fast. See the module docstring
# for why it is secondary. 1e-2 corresponds to 0.5 % variance bias in a linear
# system (exactly 0.5025 %, `em_variance_bias(1e-2)`).
GATE = 1e-2


# ════════════════════════════════════════════════════════════════════════════
# 1 · Displacement gates -- the primary criterion
# ════════════════════════════════════════════════════════════════════════════
def dt_max_thermal(delta: float, length: float, D0: float) -> float:
    """Diffusive displacement: `sqrt(2 D0 dt) <= delta * length` (per component).

    `length` is **that system's reference length**, not always `sigma`: for a
    harmonic trap it is `l_trap = sqrt(kT/k_t)`, for a repulsive pair system the
    mean spacing `d` (`.claude/rules/overdamped-stability.md`).
    """
    return (delta * length) ** 2 / (2.0 * D0)


def dt_max_force(delta: float, length: float, gamma: float,
                 max_force: float | None) -> float | None:
    """Force displacement: `(max|F|/gamma) dt <= delta * length`.

    `None` when `max_force` is missing or exactly `0` -- a straight chain sits at
    a stationary point where this gate is **powerless**, and blurring that into
    `inf` hides it. Distinguish "measured it, got 0" (physics) from "have not
    measured it" (a procedure violation); the caller writes that sentence.
    """
    if not max_force or max_force <= 0.0:
        return None
    return delta * length * gamma / max_force


def dt_max_active(delta: float, length: float, v0: float | None) -> float | None:
    """Advective displacement: `v0 dt <= delta * length`."""
    if not v0 or v0 <= 0.0:
        return None
    return delta * length / v0


def dt_max_stability(safety: float, gamma: float,
                     lambda_max: float | None) -> float | None:
    """Stiffness stability: `dt <= safety * 2 gamma / lambda_max`.

    The linear stability limit of explicit overdamped Euler. **Not an accuracy
    gate** -- past it the integration diverges, and `check_finite` does not catch
    it (a bond length of 1.4e7 is still finite). The measured threshold is
    `1.22-2.80x` this bound, so `safety = 0.2` leaves 6-14x headroom.
    Basis: `knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md`
    """
    if not lambda_max or lambda_max <= 0.0:
        return None
    return safety * 2.0 * gamma / lambda_max


# ════════════════════════════════════════════════════════════════════════════
# 2 · Euler-Maruyama variance bias in a linear system (harmonic trap)
# ════════════════════════════════════════════════════════════════════════════
#  Scheme:   x_{n+1} = x_n (1 - dt*) + sqrt(2 dt*) xi   (units of tau_trap, D* = 1)
#  Steady:   Var* = 2 dt* / (2 dt* - dt*^2) = 1/(1 - dt*/2)
#  Exact:    Var* = 1
#  =>        relative bias = (dt*/2)/(1 - dt*/2)
#
#  ★ The exact form and its first-order expansion were **both** in the repository,
#    in different packages, under different names, and nothing tied them together:
#      simbot.estimators.euler_maruyama_trap_variance_bias   exact
#      bdbot.checks.bias_from_dt                             (dt/tau)/2, linearized
#    They agree to the precision each one states (0.5025 % vs 0.5 % at dt* = 1e-2),
#    so this was never going to show up as a wrong answer -- only as two numbers
#    that drift apart the day one of them is "improved".
def em_variance_bias(dt_star: float) -> float:
    """Relative `<x^2>` bias of Euler-Maruyama on a harmonic trap. **Exact.**

    This is a **known error**, not a defect: in S7 the measured `<x^2>` is
    *supposed* to come out this much high, and its absence is the signal that
    something else is wrong.
    """
    if not 0 < dt_star < 2:
        raise ValueError(f"dt_star must be in (0,2) for stability; got {dt_star}")
    return 1.0 / (1.0 - dt_star / 2.0) - 1.0


def dt_star_for_em_bias(target_rel_bias: float) -> float:
    """`dt*` achieving a target bias -- the exact inverse of `em_variance_bias`."""
    b = target_rel_bias
    return 2.0 * b / (1.0 + b)


def em_variance_bias_linearized(dt_star: float) -> float:
    """First-order form `dt*/2`. Kept because runs depend on it, not because it
    is better.

    `bdbot.checks.bias_from_dt` and `cases/trap_2d_5um.py` were built on this,
    and `run_id` is the hash of the spec content -- so silently swapping in the
    exact form would invalidate existing runs for a `dt` change of order
    `target_bias` itself. Use `em_variance_bias` for anything new.
    """
    return dt_star / 2.0


def dt_star_for_em_bias_linearized(target_rel_bias: float) -> float:
    """`dt* = 2 b`. The inverse of `em_variance_bias_linearized`."""
    return 2.0 * target_rel_bias


def em_bias_form_gap(target_rel_bias: float) -> float:
    """Relative gap between the linearized and exact `dt*` at a given target.

    **Equals `b` exactly** (measured 2026-08-29, and it falls out algebraically:
    `|2b - 2b/(1+b)| / (2b/(1+b)) = b`). So at a 0.1 % bias target the two forms
    differ by 0.1 % in `dt` -- the gap is the target itself.
    ⚠ `b/(1+b)` is *not* the answer; it is off by `b^2`. Cheap to get wrong and
      impossible to notice, which is why this is a function and not a comment.
    """
    exact = dt_star_for_em_bias(target_rel_bias)
    return abs(dt_star_for_em_bias_linearized(target_rel_bias) - exact) / exact


# ════════════════════════════════════════════════════════════════════════════
# 3 · The timescale-ratio gate, and what it costs relative to displacement
# ════════════════════════════════════════════════════════════════════════════
def relaxation_time(gamma: float, stiffness: float) -> float:
    """`tau = gamma / k`. Trap stiffness, `U''(r_min)`, bond stiffness -- one
    structure that showed up in several cases under several names."""
    return gamma / stiffness


def dt_from_gate(tau: float, gate: float = GATE) -> float:
    """`dt = gate * tau`. The secondary criterion -- see the module docstring.

    Safe when `tau` really is the fastest mode of the system. The failure mode is
    not divergence: with `tau_D/tau_trap = 2.4e5`, a `dt` chosen from `tau_D`
    **cannot see the trap at all** and still runs to completion.
    """
    return gate * tau


def compare_criteria(*, tau_fast: float, gate: float = GATE,
                     delta: float, length: float, D0: float,
                     gamma: float | None = None,
                     max_force: float | None = None) -> dict:
    """Both criteria side by side, plus which one binds.

    `.claude/rules/overdamped-stability.md` asks for the rule you *did not* adopt
    to be recorded, because that is what makes a silent `dt` error visible after
    the fact. Returns the two bounds, their ratio, and the resulting thermal
    displacement of the ratio-gate `dt` -- the last one is the number that
    answers "would this `dt` have thrown a particle out of the box".
    """
    dt_ratio = dt_from_gate(tau_fast, gate)
    dt_disp = dt_max_thermal(delta, length, D0)
    dt_f = (None if gamma is None
            else dt_max_force(delta, length, gamma, max_force))
    candidates = {"ratio_gate": dt_ratio, "thermal_displacement": dt_disp}
    if dt_f is not None:
        candidates["force_displacement"] = dt_f
    binding = min(candidates, key=lambda k: candidates[k])
    return {
        "candidates": candidates,
        "binding": binding,
        "dt": candidates[binding],
        "ratio_over_displacement": dt_ratio / dt_disp,
        # what the ratio gate would actually have moved a particle, in units of
        # `length`. Compare against `delta`.
        "thermal_disp_of_ratio_gate": (2.0 * D0 * dt_ratio) ** 0.5 / length,
    }


__all__ = ["GATE", "dt_max_thermal", "dt_max_force", "dt_max_active",
           "dt_max_stability", "em_variance_bias", "dt_star_for_em_bias",
           "em_variance_bias_linearized", "dt_star_for_em_bias_linearized",
           "em_bias_form_gap", "relaxation_time", "dt_from_gate",
           "compare_criteria"]
