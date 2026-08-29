"""Provenance tracking -- "every number carries a source" (CLAUDE.md rule 3).

The first two cases had each grown a `Provenanced` plus a YAML node parser.

tier: 0 = directly given or handbook . 1 = literature plus verification, or a
      confirmed convention . 2 = literature, unverified . 3 = arbitrary assumption

WARNING: tier 1 by inheritance is the dangerous case. `T = 300 K` is recorded as
tier 1 across every case here and is in fact a *choice* inherited from a sketch
with no temperature -- worth -4% to -14% on every timescale, because water's
viscosity is 2.06 %/K sensitive. Inheriting is legitimate; recording it as if it
were measured is not.
"""
from __future__ import annotations

from dataclasses import dataclass

from .units import Q


@dataclass
class Provenanced:
    """Value + source + confidence tier. A bare value is never left floating."""

    value: object
    source: str
    tier: int

    def __repr__(self):
        try:
            return f"{self.value:~.4gP} [tier{self.tier}]"
        except (TypeError, ValueError):
            return f"{self.value!r} [tier{self.tier}]"


def load_node(node: dict) -> Provenanced:
    """A `{value, unit, source, tier}` YAML node -> `Provenanced[Quantity]`."""
    return Provenanced(Q(node["value"], node["unit"]), node["source"], int(node["tier"]))


def tier_summary(items) -> dict:
    """Counts per tier. Used by the validator to catch a spec built only from
    tier 2 and below.
    """
    out: dict[int, int] = {}
    for p in items:
        out[p.tier] = out.get(p.tier, 0) + 1
    return dict(sorted(out.items()))


__all__ = ["Provenanced", "load_node", "tier_summary"]
