"""Block C — the canonical dimensionless-group registry. LLM 0 lines.

## Why a registry exists

Measured across the 8 case scripts at `4b4503a`: **73 distinct `Group` names, 56
of them used exactly once.** Every case invented its own. Worse, the two names
with enough reuse to support a cross-case comparison were both broken:

    St        6 uses, TWO definitions -- tau_p/tau_B in abp_rod_2d:169,
              chain_bend_2d:470, chain_bend_dlvo_2d:397, soft_r3_2d:174 and
              trap_drag_2d:312; but tau_p/tau_bond in chain_relax_2d_dlvo:312
    Pe        ONE group under two names -- "Pe" (abp_rod_2d:155) and "Pe_drag"
              (trap_drag_2d:282), both tau_B/tau_v

So aggregating either today gives a wrong answer, quietly. `Group.mismatch()`
cannot catch it: each is a true ratio of its own declared num/den.

## What this module does and does not do

It **does not** replace `nondim.Group`. A case still builds `Group` objects from
its own ledger, and `NondimSpec.validate()` still recomputes every ratio. This
module adds one thing: if a case uses a **registered name**, the roles of its
numerator and denominator must match the registry, and any ambiguous group must
name its free choice. A case keeps the right to declare a bespoke group; it loses
the right to call a bespoke group `St`.

## Admission rule

A group enters `REGISTRY` at its **second independent use** — the repository's
existing promotion rule, unchanged (`docs/01-architecture.md` §"has it appeared
twice?"). Everything in `REGISTRY` below is either already used twice or is a
textbook group the roadmap names (`Re`, `Pe_shear`, `Wi`). Those three are marked
`uses=0` honestly rather than being presented as established here.

## Reynolds is not new

`simbot/spec.py:598-609` already computes and gates `stokes_reynolds` at the same
`1e-2` threshold. `Re` below records that, and a case should call that check
rather than re-implement it.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

SCHEMA = "bdbot.groups/0.1"

#: Comparison operators a threshold may use.
OPS = ("<=", ">=")


@dataclass(frozen=True)
class GroupDef:
    """One canonical dimensionless group.

    `num_role` / `den_role` are **ledger roles**, not symbols: `St` is
    `tau_p / tau_gov`, and which time plays `tau_gov` is the case's choice. What
    the registry fixes is that the numerator is the inertial time and the
    denominator is *the governing time the case named*.
    """
    symbol: str
    expr: str
    num_role: str
    den_role: str
    meaning: str
    requires: tuple            # the SI inputs needed to form it
    limit: float | None = None
    op: str = "<="
    mismatch_means: str = ""
    free_choice: str = ""      # non-empty => the case MUST name this
    undefined_when: str = ""   # when the group does not exist for a system
    aliases: tuple = ()        # historical names that must resolve here
    uses: int = 0              # independent case uses at 4b4503a
    note: str = ""

    @property
    def ambiguous(self) -> bool:
        """Does using this name require naming a free choice (e.g. `tau_gov`)."""
        return bool(self.free_choice)


REGISTRY: dict = {d.symbol: d for d in (
    GroupDef(
        symbol="St", expr="tau_p / tau_gov",
        num_role="inertia", den_role="governing",
        meaning="particle inertia vs the governing dynamics",
        requires=("m", "gamma", "tau_gov"),
        limit=1e-2, op="<=",
        mismatch_means="inertia is not negligible: overdamped BD is invalid, "
                       "and HOOMD's Brownian integrator cannot represent it",
        free_choice="tau_gov",
        uses=6,
        note="★ THE collision. 5 cases mean tau_p/tau_B, 1 means tau_p/tau_bond. "
             "Both are defensible Stokes numbers and they are different "
             "quantities. Naming tau_gov is now mandatory.",
    ),
    GroupDef(
        symbol="Pe_adv", expr="v d / D_t  =  tau_B / tau_v",
        num_role="diffusion", den_role="advection",
        meaning="advection vs diffusion",
        requires=("v", "d", "D_t"),
        mismatch_means="no threshold: this is a regime label, not a validity check",
        aliases=("Pe", "Pe_drag"),
        uses=2,
        note="★ THE split. `Pe` (abp_rod_2d:155) and `Pe_drag` (trap_drag_2d:282) "
             "are the same group. Both resolve here through `aliases`.",
    ),
    GroupDef(
        symbol="Pe_shear", expr="gammadot d^2 / (4 D_t)",
        num_role="shear", den_role="diffusion",
        meaning="shear vs diffusion",
        requires=("gammadot", "d", "D_t"),
        undefined_when="no imposed deformation rate (any quiescent or "
                       "passive-probe system) -- declare_absent it",
        uses=0,
        note="Distinct symbol from Pe_adv on purpose. An unsubscripted `Pe` in a "
             "sheared system is the second-most common way to mean two things.",
    ),
    GroupDef(
        symbol="Re", expr="rho_f v d / eta",
        num_role="inertia_fluid", den_role="viscous",
        meaning="fluid inertia vs viscous stress",
        requires=("rho_f", "v", "d", "eta"),
        limit=1e-2, op="<=",
        mismatch_means="Stokes drag is suspect; gamma = 3 pi eta d no longer holds",
        free_choice="v",
        uses=0,
        note="Already implemented as `stokes_reynolds` in simbot/spec.py:598-609 "
             "at this same 1e-2 limit. Call that, do not re-derive it.",
    ),
    GroupDef(
        symbol="Wi", expr="lambda_p * gammadot",
        num_role="polymer_relaxation", den_role="deformation",
        meaning="elastic stretch vs relaxation",
        requires=("lambda_p", "gammadot"),
        undefined_when="a Newtonian solvent has no polymer relaxation time, and a "
                       "passive (unforced) measurement has no deformation rate -- "
                       "Wi does not exist for either. declare_absent it",
        uses=0,
        note="Wi uses the deformation RATE; De uses a ratio of times. They "
             "coincide in steady shear and do not in general.",
    ),
    GroupDef(
        symbol="De", expr="lambda_relax / t_obs   (lambda * omega if oscillatory)",
        num_role="relaxation", den_role="observation",
        meaning="relaxation vs the observation or driving timescale",
        requires=("lambda_relax", "t_obs or omega"),
        free_choice="lambda_relax",
        uses=2,
        note="chain_bend_2d declares De, De_chain_old and De_trap; "
             "chain_bend_dlvo_2d declares De_bond, De_chain_diff, De_trap. Three "
             "relaxation times in one case is legitimate -- which is exactly why "
             "lambda_relax must be named.",
    ),
    GroupDef(
        symbol="phi", expr="packing fraction",
        num_role="occupied", den_role="box",
        meaning="crowding",
        requires=("N", "d", "box"),
        undefined_when="no pair interaction -- the value is computable and "
                       "meaningless (simbot/nondim.py already guards this with "
                       "`if spec.pair`)",
        uses=2,
    ),
    GroupDef(
        symbol="k_star", expr="k d^2 / kT",
        num_role="trap_energy", den_role="thermal",
        meaning="trap stiffness vs thermal fluctuation",
        requires=("k", "d", "kT"),
        aliases=("k*",),
        uses=3,
        note="l_k/d = 1/sqrt(k_star) and tau_k/tau_B = 1/k_star are derived from "
             "it; both appear as separate groups in trap_2d_5um and are not "
             "registered separately.",
    ),
)}

#: alias -> canonical symbol. Built once so a lookup cannot drift.
ALIASES: dict = {a: d.symbol for d in REGISTRY.values() for a in d.aliases}


def canonical(symbol: str) -> str | None:
    """Canonical registry symbol for `symbol`, or None if it is bespoke.

    Bespoke is allowed. `None` means "the registry has no opinion", not "wrong".
    """
    s = str(symbol).strip()
    if s in REGISTRY:
        return s
    return ALIASES.get(s)


def get(symbol: str) -> GroupDef | None:
    c = canonical(symbol)
    return REGISTRY[c] if c else None


@dataclass
class Usage:
    """How a case actually used a registered name. What `validate_usage` checks."""
    symbol: str
    value: float
    num_symbol: str = ""       # the ledger symbol used as numerator
    den_symbol: str = ""       # ... and denominator
    free_choice: str = ""      # what the case named for `free_choice`
    justification: str = ""    # required when using a group declared undefined


def validate_usage(usages, *, strict_unregistered=False) -> list:
    """Check a case's group usage against the registry. Returns a list of strings;
    empty means clean.

    Two rules, and both come from a measured failure:

      1. An **ambiguous** group must name its free choice. `St` without a named
         `tau_gov` is the collision above, and it is an error rather than a
         warning because the two meanings differ by orders of magnitude
         (`tau_B / tau_bond` in `chain-relax` is not O(1)).
      2. An **alias** must be reported, so `Pe` and `Pe_drag` aggregate together
         rather than as two groups of one.

    `strict_unregistered=True` additionally flags bespoke names. Off by default:
    bespoke groups are legitimate and 56 of them exist.
    """
    out = []
    seen_canon: dict = {}
    for u in usages:
        c = canonical(u.symbol)
        if c is None:
            if strict_unregistered:
                out.append(f"[info] '{u.symbol}' is not in the registry "
                           f"(bespoke groups are allowed; flagged because "
                           f"strict_unregistered=True)")
            continue
        d = REGISTRY[c]
        if c != u.symbol:
            out.append(f"[info] '{u.symbol}' is an alias of '{c}' -- aggregate "
                       f"under '{c}'")
        if d.ambiguous and not str(u.free_choice).strip():
            out.append(f"[error] '{u.symbol}' must name its {d.free_choice}: "
                       f"{d.note.splitlines()[0] if d.note else d.meaning}")
        if d.undefined_when and not str(u.justification).strip():
            out.append(f"[warn] '{u.symbol}' may be undefined for this system "
                       f"({d.undefined_when}). Give a justification or "
                       f"declare it absent.")
        # the same canonical group used twice with different free choices is
        # legitimate (three De's in one case) but must be distinguishable
        key = (c, str(u.free_choice).strip())
        if key in seen_canon:
            what = d.free_choice or "no free choice"
            out.append(f"[error] '{u.symbol}' with {what}="
                       f"'{u.free_choice}' is declared twice")
        seen_canon[key] = u
    return out


def errors(messages) -> list:
    return [m for m in messages if m.startswith("[error]")]


# -- the formulas, so a case computes rather than retypes ---------------------

def st(tau_p: float, tau_gov: float) -> float:
    """`St = tau_p / tau_gov`. The denominator is the case's named governing time."""
    if tau_gov <= 0:
        raise ValueError(f"tau_gov must be > 0, got {tau_gov}")
    return float(tau_p) / float(tau_gov)


def pe_adv(v: float, d: float, D_t: float) -> float:
    """`Pe_adv = v d / D_t`. Equals `tau_B/tau_v` with `tau_B = d^2/D_t`."""
    if D_t <= 0:
        raise ValueError(f"D_t must be > 0, got {D_t}")
    return float(v) * float(d) / float(D_t)


def pe_shear(gammadot: float, d: float, D_t: float) -> float:
    """`Pe_shear = gammadot d^2 / (4 D_t)`. The 4 is the a^2 = (d/2)^2 convention."""
    if D_t <= 0:
        raise ValueError(f"D_t must be > 0, got {D_t}")
    return float(gammadot) * float(d) ** 2 / (4.0 * float(D_t))


def reynolds(rho_f: float, v: float, d: float, eta: float) -> float:
    """`Re = rho_f v d / eta`. Same definition and limit as
    `simbot/spec.py`'s `stokes_reynolds`, which uses the radius; this uses the
    diameter, so the two differ by 2. **Call simbot's check, not this, when a
    verdict is wanted** -- this exists so a report can print the number."""
    if eta <= 0:
        raise ValueError(f"eta must be > 0, got {eta}")
    return float(rho_f) * float(v) * float(d) / float(eta)


def wi(lambda_p: float, gammadot: float) -> float:
    """`Wi = lambda_p gammadot`. Undefined for a Newtonian solvent or an unforced
    measurement -- the caller must not reach here in those cases."""
    return float(lambda_p) * float(gammadot)


def de(lambda_relax: float, t_obs: float = None, omega: float = None) -> float:
    """`De = lambda/t_obs`, or `lambda*omega` when oscillatory. Exactly one of
    `t_obs` and `omega`."""
    if (t_obs is None) == (omega is None):
        raise ValueError("give exactly one of t_obs and omega")
    if omega is not None:
        return float(lambda_relax) * float(omega)
    if t_obs <= 0:
        raise ValueError(f"t_obs must be > 0, got {t_obs}")
    return float(lambda_relax) / float(t_obs)


def k_star(k: float, d: float, kT: float) -> float:
    """`k* = k d^2 / kT`."""
    if kT <= 0:
        raise ValueError(f"kT must be > 0, got {kT}")
    return float(k) * float(d) ** 2 / float(kT)


def overdamped_verdict(tau_p: float, tau_gov: float) -> tuple:
    """`(ok, St, limit, margin)` for the overdamped/underdamped decision.

    ★ This is **deterministic**. The model's job is to choose `tau_gov`; whether
    `St` clears `1e-2` is arithmetic, and `trap-2d-5um` already carries it as a
    check (`8.139e-04` against `1e-02`, margin 12.286).
    """
    d = REGISTRY["St"]
    s = st(tau_p, tau_gov)
    ok = s <= d.limit
    margin = (d.limit / s) if s > 0 else math.inf
    return ok, s, d.limit, margin


def render_registry() -> str:
    L = ["=" * 78, "canonical dimensionless groups", "=" * 78]
    for s, d in REGISTRY.items():
        lim = f"{d.op} {d.limit:g}" if d.limit is not None else "—"
        L.append(f"\n{s:<10} {d.expr}")
        L.append(f"           {d.meaning}   [limit {lim}]   uses={d.uses}")
        if d.aliases:
            L.append(f"           aliases: {', '.join(d.aliases)}")
        if d.free_choice:
            L.append(f"           ★ must name: {d.free_choice}")
        if d.undefined_when:
            L.append(f"           undefined when: {d.undefined_when}")
    L.append("=" * 78)
    return "\n".join(L)


__all__ = ["SCHEMA", "OPS", "GroupDef", "REGISTRY", "ALIASES", "Usage",
           "canonical", "get", "validate_usage", "errors", "st", "pe_adv",
           "pe_shear", "reynolds", "wi", "de", "k_star", "overdamped_verdict",
           "render_registry"]
