"""The crossing point between the two halves of the core.

**This module is the "two engines" seam, closed** (docs/00-merge-decisions.md
section 5). Before it, the two halves could not reach each other:

    a case run through `bdbot`     got a health verdict but **no sealed prediction**
    a spec run through `simbot`    got sealing but **could not reach bdbot's cases**

Both halves model the same thing -- an SI physical system whose every number
carries a provenance (CLAUDE.md rule 3) -- in two vocabularies:

    bdbot.physical.PhysicalSystem   a `system.yaml` tree of `Provenanced`
                                    leaves: pint value + free-text `source` +
                                    ordinal `tier` 0-3
    simbot.spec.SystemSpec          typed dataclasses of `Quantity` leaves: raw
                                    float + `unit` string + categorical
                                    `provenance` + `basis`

**Why the models were NOT collapsed into one.** They are the input side of two
content-addressed archives. `PhysicalSystem` feeds `bdbot.nondim` ->
`specs/<run_id>.json`, and `run_id` is the hash of the spec content: 278 spec
files and 263 run directories are *named* by it. `SystemSpec` feeds
`simbot.nondim.reduce_spec` and the sealed S2 prediction documents, which
`SEALED.sha256` pins. Unifying the serialisation would rename one archive or break
the other's seals. So instead: **one definition of the physics, two
serialisations, and a verified crossing.** The physics was already single-sourced
by the earlier merge (`bdbot.constants`, `bdbot.dt`,
`tests/test_cross_package_equivalence.py`); what was still duplicated is the
*meaning of a provenance leaf*, and that is what lives here.

⚠ **The tier <-> provenance mapping is not a bijection, and pretending otherwise
  would be an invention path.** Ten `provenance` kinds collapse onto four `tier`
  values, so `provenance -> tier -> provenance` is lossy by construction. Carrying
  a value across the seam and *re-guessing* its tier on the far side is exactly the
  "never invent a value" failure (rule 3), which is why the mapping is written down
  once here rather than inferred at each call site.
  `tests/test_bridge_equivalence.py` pins both the forward map and the lossiness.
"""
from __future__ import annotations

import math

from bdbot import materials as _mat
from bdbot.provenance import Provenanced
from bdbot.units import Q as _pintQ

from .spec import (Geometry, Medium, Quantity, Species, SystemSpec,
                   PROVENANCE_KINDS, derive)

# ════════════════════════════════════════════════════════════════════════════
# 1 · the provenance vocabulary, mapped once
# ════════════════════════════════════════════════════════════════════════════
#  tier 0 = directly given or handbook
#       1 = literature + verification, or a confirmed convention
#       2 = literature, unverified
#       3 = arbitrary assumption
#  (bdbot/provenance.py; `TIER_MEANING` in bdbot/physical.py says the same)
PROVENANCE_TO_TIER: dict[str, int] = {
    # tier 0 -- the value was given to us, not reasoned to
    "from_drawing": 0,     # read off the sketch
    "observation": 0,      # read off the source material, any modality
    "measured": 0,         # an experimental value
    "user": 0,             # a person fixed it for this run. Deliberately NOT an
                           # assumption -- S7b excludes it from the sensitivity sweep
    "derived": 0,          # computed from other fields by a simbot function.
                           # bdbot's equivalent is DERIVED_SECTIONS, which are
                           # exempt from needing a source because `verify()`
                           # recomputes them
    # tier 1 -- literature or convention, with something backing it
    "from_paper": 1,       # a distillation in knowledge/source/papers
    "from_knowledge": 1,   # a knowledge/wiki entry
    "rule": 1,             # a confirmed convention from config/run_policy.yaml
    # tier 2 -- reasoned, unverified
    "inference": 2,
    # tier 3 -- filled because the source did not say
    "assumed": 3,
}

#  The reverse is a **choice of representative**, because the forward map is
#  many-to-one. The representative is the most conservative kind at that tier:
#  it must never overstate where a number came from.
TIER_TO_PROVENANCE: dict[int, str] = {
    0: "observation",      # not "measured" -- claiming measurement is a stronger
                           # statement than tier 0 supports
    1: "from_knowledge",   # not "from_paper" -- naming a paper we cannot cite is worse
    2: "inference",
    3: "assumed",
}

#  Kinds that vanish under a round trip. Listed rather than derived so the loss is
#  visible in the source, not just computable from it.
LOSSY_UNDER_ROUNDTRIP = tuple(
    sorted(k for k, t in PROVENANCE_TO_TIER.items() if TIER_TO_PROVENANCE[t] != k))


def tier_of(provenance: str) -> int:
    """`provenance` kind -> `tier`. Raises on an unregistered kind.

    Refusing an unknown kind is deliberate: defaulting it to a tier would silently
    launder an unrecognised provenance into a confident one.
    """
    if provenance not in PROVENANCE_KINDS:
        raise ValueError(f"{provenance!r} is not a registered provenance kind; "
                         f"known: {sorted(PROVENANCE_KINDS)}")
    if provenance not in PROVENANCE_TO_TIER:
        raise ValueError(
            f"provenance {provenance!r} is registered in simbot.spec but has no "
            f"tier here. Add it to PROVENANCE_TO_TIER with a reason -- do not let "
            f"it fall through to a default.")
    return PROVENANCE_TO_TIER[provenance]


def provenance_of(tier: int) -> str:
    """`tier` -> the conservative representative `provenance` kind."""
    if tier not in TIER_TO_PROVENANCE:
        raise ValueError(f"tier {tier!r} is outside 0-3")
    return TIER_TO_PROVENANCE[tier]


# ════════════════════════════════════════════════════════════════════════════
# 2 · leaf conversion
# ════════════════════════════════════════════════════════════════════════════
def quantity_to_provenanced(q: Quantity) -> Provenanced:
    """`simbot` leaf -> `bdbot` leaf. The value becomes a pint Quantity.

    `basis` carries into `source`, because that is the field a human reads to
    decide whether the number is defensible.
    """
    unit = q.unit or "dimensionless"
    return Provenanced(value=_pintQ(q.value, unit),
                       source=q.basis or f"(no basis recorded; provenance={q.provenance})",
                       tier=tier_of(q.provenance))


def provenanced_to_quantity(p: Provenanced, **kw) -> Quantity:
    """`bdbot` leaf -> `simbot` leaf. **Converts to SI base units** before stripping.

    ⚠⚠ The conversion is the whole point, and getting it wrong is silent. ⚠⚠
      The two models have different unit contracts:
        `bdbot`   a pint Quantity -- the unit travels WITH the value, so
                  `0.851 mPa*s` is a complete statement
        `simbot`  fields named `*_si`, and `Quantity.si` returns `float(value)`
                  with **no conversion** -- so the value must ALREADY be in SI and
                  `unit` is documentation only
      Handing over `(0.851, "mPa*s")` therefore reads back as `0.851 Pa*s`:
      **gamma comes out 1000x too large, D 1000x too small, every timescale
      1000x off.** Measured, not reasoned -- `derived_agreement()` caught exactly
      this on 5 of the 8 real `system.yaml` files when this function was first
      written, which is why that check exists.

    ⚠ Lossy in the provenance field -- see the module docstring. The `basis`
      records both the tier and the original unit, so nothing a human would need is
      dropped even though the machine-readable kind is a representative.
    """
    val = p.value
    if hasattr(val, "magnitude"):
        as_read = f"{val:~}"
        si = val.to_base_units()
        magnitude, unit = float(si.magnitude), f"{si.units:~}"
        unit_note = f"; SI base from {as_read}" if str(si.units) != str(val.units) else ""
    else:
        magnitude, unit, unit_note = val, "", ""
    kind = provenance_of(p.tier)
    return Quantity(value=magnitude, unit=unit, provenance=kind,
                    basis=kw.pop("basis", None) or (
                        f"{p.source} [via bdbot tier {p.tier}{unit_note}; the "
                        f"categorical provenance is a representative, not a source "
                        f"claim]"),
                    **kw)


# ════════════════════════════════════════════════════════════════════════════
# 3 · the crossing: PhysicalSystem -> SystemSpec
# ════════════════════════════════════════════════════════════════════════════
def _scalar_or_refuse(p: Provenanced, name: str, sweep_index: int | None):
    """A core leaf that holds a sweep is N systems, not one. Refuse or index it.

    ★ Found by measurement, 2026-08-29: `bdbot`'s L2 model allows a **parameter
      sweep inside a leaf** -- `chain-bend-2d-dlvo` has `particle.count =
      [5, 9, 15, 25]` and `network` has `[512, 1528]`. `SystemSpec` describes
      exactly one system. That is a **structural asymmetry between the two
      models**, not a formatting difference, and it is the reason this crossing
      cannot be a pure field-by-field map.
      Taking element 0 silently would be inventing which system was meant, so the
      caller has to say. Three of the eight cases hit this.
    """
    val = p.value
    mag = getattr(val, "magnitude", val)
    if getattr(mag, "shape", ()) == ():
        return p
    try:
        n = len(mag)
    except TypeError:
        return p
    if sweep_index is None:
        raise ValueError(
            f"{name} holds a sweep of {n} values ({list(mag)!r}), and a SystemSpec "
            f"describes ONE system. Pass sweep_index=0..{n - 1} to choose which. "
            f"Do not assume 0 -- that is picking a physical system by accident.")
    if not 0 <= sweep_index < n:
        raise ValueError(f"sweep_index {sweep_index} is outside 0..{n - 1} for {name}")
    picked = mag[sweep_index]
    unit = getattr(val, "units", None)
    return Provenanced(
        value=(_pintQ(picked, str(unit)) if unit is not None else picked),
        source=f"{p.source} [sweep element {sweep_index} of {n}: {list(mag)!r}]",
        tier=p.tier)


def physical_to_systemspec(ps, *, card: str = "", question: str = "",
                           sweep_index: int | None = None) -> SystemSpec:
    """A `bdbot` L2 `PhysicalSystem` -> a `simbot` `SystemSpec`.

    This is the direction that closes the seam: once a `bdbot` case's
    `system.yaml` is a `SystemSpec`, that case can use `simbot.io.write_seal`
    (prediction sealing), `simbot.validate` (PASS / FAIL / **INCONCLUSIVE**) and
    `simbot.report` -- none of which `bdbot` has.

    **Only the five core leaves cross** (`d`, `rho_p`, `N`, `T`, `eta`), because
    those are the ones `bdbot.physical.CORE_PROVENANCED` guarantees are present and
    parsed in every case. Interactions do not: `bdbot` stores them as free-form
    YAML per case while `SystemSpec` needs typed `PairInteraction` /
    `BondInteraction` / `ExternalField`, and inventing that structure per case is
    the invention path rule 3 forbids. So the result is a **minimal** spec, good
    for the derived scales and the overdamped/Stokes gates, and it says so rather
    than pretending to be complete.

    Raises `ValueError` if a core leaf is missing -- a half-populated spec that
    then silently derives the wrong `gamma` is worse than a refusal.
    """
    missing = [k for k in ("d", "T", "eta", "N") if k not in ps.core]
    if missing:
        raise ValueError(
            f"cannot cross the seam: {ps.path} is missing core leaf/leaves "
            f"{missing}. bdbot.physical.CORE_PROVENANCED expects "
            f"particle.diameter, particle.count, medium.temperature, "
            f"medium.viscosity (and particle.density for the inertial gate).")

    d = ps.core["d"]
    radius = Provenanced(value=(d.value / 2).to("m"),
                         source=f"{d.source} [halved: SystemSpec takes a radius, "
                                f"bdbot stores a diameter]",
                         tier=d.tier)
    species = Species(
        name=str(ps.raw.get("particle", {}).get("material", "primary")),
        n_simulated=provenanced_to_quantity(
            _scalar_or_refuse(ps.core["N"], "particle.count", sweep_index)),
        radius_si=provenanced_to_quantity(radius),
        density_si=(provenanced_to_quantity(ps.core["rho_p"])
                    if "rho_p" in ps.core else None),
    )
    medium = Medium(T_si=provenanced_to_quantity(ps.core["T"]),
                    eta_si=provenanced_to_quantity(ps.core["eta"]),
                    species=Quantity(value=str(ps.raw.get("medium", {}).get("name", "")),
                                     provenance="observation",
                                     basis="medium.name from system.yaml"))
    dim = int(ps.raw.get("dimensions", 2))
    geometry = Geometry(
        dim=Quantity(value=dim, provenance="observation",
                     basis="dimensions from system.yaml"),
        boundary=Quantity(value=str(ps.raw.get("boundary", "periodic")),
                          provenance="assumed",
                          basis="bdbot's system.yaml does not record a boundary "
                                "condition; periodic is this project's convention"),
    )
    return SystemSpec(
        card=card or str(ps.raw.get("label", "")),
        question=question or str(ps.raw.get("description", "")),
        geometry=geometry, species=[species], medium=medium,
        notes=[f"crossed from bdbot.physical via simbot.bridge: {ps.path}",
               "MINIMAL -- interactions did not cross (see physical_to_systemspec)"],
    )


# ════════════════════════════════════════════════════════════════════════════
# 4 · the behavioural pin
# ════════════════════════════════════════════════════════════════════════════
#  ★ The check with power here is **not** syntactic. A peer session's
#    token-stream comparator proves "only comments and strings changed", which is
#    the wrong claim for a data-model crossing: it is not invariant to moving code
#    between files, and it would accept a silently changed provenance string.
#    What matters is that both halves compute the *same physics* from the *same
#    system.yaml*. So: derive the SI bundle down both paths and diff it.
def derived_agreement(ps, *, rtol: float = 0.0,
                      sweep_index: int | None = None) -> dict:
    """Derive the SI scales down both paths from one `PhysicalSystem` and compare.

    `rtol = 0.0` means **bit-identical**, which is the bar: the two paths already
    share `bdbot.constants` and the Stokes/Einstein relations, so any difference is
    a defect in this bridge rather than a rounding question.

    Returns per-quantity `{bdbot, simbot, rel_gap}` plus `all_agree`.
    """
    bulk = ps.bulk()
    if bulk is None:
        return {"all_agree": None,
                "why": "bdbot.physical.bulk() returned None -- needs d, T, eta, rho_p"}
    spec = physical_to_systemspec(ps, sweep_index=sweep_index)
    d = derive(spec)

    pairs = {
        "kT":    (float(bulk["kT"].to("J").magnitude),        d["kT_si"]),
        "gamma": (float(bulk["gamma"].to("kg/s").magnitude),  d["gamma_si"]),
        "D_t":   (float(bulk["D_t"].to("m^2/s").magnitude),   d["D0_si"]),
        "tau_B": (float(bulk["tau_B"].to("s").magnitude),     d["tau_D_si"]),
    }
    if "m" in bulk and "mass_si" in d:
        pairs["mass"] = (float(bulk["m"].to("kg").magnitude), d["mass_si"])
        pairs["tau_p"] = (float(bulk["tau_p"].to("s").magnitude), d["tau_inertial_si"])

    out, agree = {}, True
    for name, (a, b) in pairs.items():
        gap = 0.0 if a == b else abs(a - b) / abs(b) if b else math.inf
        ok = (a == b) if rtol == 0.0 else gap <= rtol
        agree &= ok
        out[name] = {"bdbot": a, "simbot": b, "rel_gap": gap, "identical": a == b}
    return {"all_agree": bool(agree), "quantities": out,
            "label": spec.card, "path": str(ps.path)}


__all__ = ["PROVENANCE_TO_TIER", "TIER_TO_PROVENANCE", "LOSSY_UNDER_ROUNDTRIP",
           "tier_of", "provenance_of", "quantity_to_provenanced",
           "provenanced_to_quantity", "physical_to_systemspec", "derived_agreement"]
