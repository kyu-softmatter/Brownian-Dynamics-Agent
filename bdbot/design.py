"""Block D — which dimensionless groups the design preserves, and which it gave
up. LLM 0 lines.

## Why abandonment has to be written down

Matching every group at once is generally **over-determined**. BD is routinely run
at an unphysical Reynolds number on purpose; a bead-spring chain does not match a
real polymer's `N`; a 2D simulation of a 3D experiment abandons the whole
dimensionality. None of that is wrong. What is wrong is doing it silently, because
the failure looks exactly like success: the agent matches the subset it can hit
and reports "matched", and the groups it could not hit are simply absent from the
report.

So this module makes the union a **hard requirement**: every group the registry
says the system supports must appear in `matched` or in `abandoned`, and an
abandoned group must carry both a `reason` and a `safe_because`. A group in
neither is an error, not a warning.

That is `scales.declare_absent(role, reason)` — which already raises on an empty
reason — lifted from ledger *roles* to *groups*.

## `safe_because` is a claim, not a formality

If a case says a group is safe to abandon *because* some other check clears, and
that check fails, the two disagree and the design is wrong. `verify_safety` takes
the case's real checks and looks for that contradiction. `A1`'s rule applies: if
two pieces of evidence disagree, stop — do not average them.

## What this does NOT do

It does not decide *whether* abandoning is safe. That is a judgment plus a
threshold, and the model may make the judgment while the code owns the threshold
(CLAUDE.md, the division of labour). This module records the judgment, checks the
bookkeeping is complete, and reports contradictions.
"""
from __future__ import annotations

from dataclasses import dataclass, field

from . import groups as _groups

SCHEMA = "bdbot.design/0.1"

#: Relative tolerance for calling an experiment/simulation group pair "matched".
#: Loose on purpose: matching a group to 1 % is already a strong statement, and a
#: tighter default would push cases into abandoning groups they had in fact hit.
MATCH_RTOL = 1e-2


@dataclass
class Matched:
    symbol: str
    experiment: float
    simulation: float

    @property
    def rel_err(self) -> float:
        e = float(self.experiment)
        if e == 0.0:
            return 0.0 if float(self.simulation) == 0.0 else float("inf")
        return abs(float(self.simulation) - e) / abs(e)

    @property
    def ok(self) -> bool:
        return self.rel_err <= MATCH_RTOL

    def to_json(self) -> dict:
        return {"symbol": self.symbol, "experiment": float(self.experiment),
                "simulation": float(self.simulation),
                "rel_err": round(self.rel_err, 9), "ok": self.ok}


@dataclass
class Abandoned:
    symbol: str
    reason: str
    safe_because: str
    experiment: float | None = None
    simulation: float | None = None
    #: names of checks that must pass for `safe_because` to hold. Optional, but
    #: it is what turns the claim into something `verify_safety` can test.
    depends_on_checks: tuple = ()

    def __post_init__(self):
        if not str(self.reason).strip():
            raise ValueError(
                f"abandoning '{self.symbol}' requires a reason -- same rule as "
                f"scales.declare_absent (rule 3)")
        if not str(self.safe_because).strip():
            raise ValueError(
                f"abandoning '{self.symbol}' requires `safe_because`: why is it "
                f"safe to give this group up in THIS system? A reason says what "
                f"you did; safe_because says why the result survives it.")

    def to_json(self) -> dict:
        return {"symbol": self.symbol, "reason": self.reason,
                "safe_because": self.safe_because,
                "experiment": self.experiment, "simulation": self.simulation,
                "depends_on_checks": list(self.depends_on_checks)}


@dataclass
class DesignLedger:
    """The two-column record. `supported` is what the registry says this system
    can form -- the coverage requirement is measured against it."""
    case: str
    supported: tuple = ()
    matched: list = field(default_factory=list)
    abandoned: list = field(default_factory=list)

    # -- construction -------------------------------------------------------
    def match(self, symbol: str, experiment: float, simulation: float):
        self.matched.append(Matched(symbol, experiment, simulation))
        return self

    def abandon(self, symbol: str, *, reason: str, safe_because: str,
                experiment=None, simulation=None, depends_on_checks=()):
        self.abandoned.append(Abandoned(
            symbol, reason, safe_because, experiment, simulation,
            tuple(depends_on_checks)))
        return self

    # -- reading ------------------------------------------------------------
    def covered(self) -> set:
        return ({m.symbol for m in self.matched}
                | {a.symbol for a in self.abandoned})

    def uncovered(self) -> list:
        """Registry-supported groups in neither column. **The silent failure.**"""
        return sorted(set(self.supported) - self.covered())

    def double_counted(self) -> list:
        return sorted({m.symbol for m in self.matched}
                      & {a.symbol for a in self.abandoned})

    def validate(self) -> list:
        """Bookkeeping errors. Returns strings; empty means clean."""
        out = []
        for s in self.uncovered():
            out.append(f"[error] '{s}' is supported by this system but appears in "
                       f"neither matched nor abandoned. A group in neither is how "
                       f"a partial match gets reported as a match.")
        for s in self.double_counted():
            out.append(f"[error] '{s}' is both matched and abandoned")
        for m in self.matched:
            if not m.ok:
                out.append(f"[error] '{m.symbol}' is listed as matched but differs "
                           f"by {m.rel_err * 100:.2f} % (limit "
                           f"{MATCH_RTOL * 100:g} %) -- move it to abandoned with "
                           f"a reason, or fix the design")
        for a in self.abandoned:
            c = _groups.canonical(a.symbol)
            if c and _groups.REGISTRY[c].undefined_when and not a.experiment:
                out.append(f"[info] '{a.symbol}' may simply not exist for this "
                           f"system ({_groups.REGISTRY[c].undefined_when}) -- that "
                           f"is `undefined`, not `abandoned`. Consider "
                           f"ScaleLedger.declare_absent instead.")
        return out

    def verify_safety(self, checks) -> list:
        """Cross-check each `safe_because` against the case's real checks.

        `checks` is any iterable of objects or dicts with `name` and `ok`. If an
        abandonment depends on a check that FAILED, the design's own justification
        contradicts its own measurement -- report it rather than averaging (`A1`).
        """
        by_name = {}
        for c in checks or ():
            n = c.get("name") if isinstance(c, dict) else getattr(c, "name", None)
            o = c.get("ok") if isinstance(c, dict) else getattr(c, "ok", None)
            if n is not None:
                by_name[str(n)] = bool(o)
        out = []
        for a in self.abandoned:
            for want in a.depends_on_checks:
                if want not in by_name:
                    out.append(f"[warn] '{a.symbol}' says it is safe because of "
                               f"check '{want}', which this case does not declare")
                elif not by_name[want]:
                    out.append(f"[error] '{a.symbol}' is abandoned on the grounds "
                               f"that '{want}' holds, and '{want}' FAILED. The "
                               f"justification contradicts the measurement -- do "
                               f"not average, decide (A1).")
        return out

    def to_json(self) -> dict:
        return {"schema": SCHEMA, "case": self.case,
                "supported": list(self.supported),
                "groups_matched": [m.to_json() for m in self.matched],
                "groups_abandoned": [a.to_json() for a in self.abandoned],
                "uncovered": self.uncovered()}

    def render(self) -> str:
        L = ["=" * 78, f"design ledger — {self.case}", "=" * 78]
        L.append("\nMATCHED")
        if not self.matched:
            L.append("  (none)")
        for m in self.matched:
            mark = "✓" if m.ok else "✗"
            L.append(f"  {mark} {m.symbol:<12} exp {m.experiment:<12.6g} "
                     f"sim {m.simulation:<12.6g} err {m.rel_err * 100:.3f} %")
        L.append("\nABANDONED (deliberately given up -- with why that is safe)")
        if not self.abandoned:
            L.append("  (none)")
        for a in self.abandoned:
            L.append(f"  · {a.symbol}")
            L.append(f"      reason        {a.reason[:58]}")
            L.append(f"      safe because  {a.safe_because[:58]}")
            if a.depends_on_checks:
                L.append(f"      depends on    {', '.join(a.depends_on_checks)}")
        msgs = self.validate()
        L.append("")
        L.append("=" * 78)
        if any(m.startswith("[error]") for m in msgs):
            L.append("VERDICT: FAIL")
        else:
            L.append("VERDICT: OK -- every supported group is accounted for")
        for m in msgs:
            L.append(f"  {m}")
        L.append("=" * 78)
        return "\n".join(L)


def errors(messages) -> list:
    return [m for m in messages if m.startswith("[error]")]


__all__ = ["SCHEMA", "MATCH_RTOL", "Matched", "Abandoned", "DesignLedger", "errors"]
