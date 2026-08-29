"""Scale-separation checks -- skill `bd-physics` section 4.

What both of the first two cases used identically:
`Check(kind, name, value, limit, op)` plus `ok`/`margin`, and the convention of
classifying every check as one of **model / integration / geometry / statistics**.

The `hard` flag became necessary in the second case -- in the first, every check
passed, so the hard/soft distinction never surfaced. bd-physics section 4 defines
model, integration and geometry as hard failures, and statistics and finite-size
as warnings.

`dt` is set to **at most 1% of the fastest physical timescale.** In both cases
that timescale had the form `gamma/(local stiffness)` -- `gamma/k` for the trap,
`gamma/U''(r_min)` for the soft pair.

WARNING: the candidate list is not static. Whenever a knob can reorder the
timescales, re-derive it -- with a trap stiffness scaled by ~200 the trap becomes
the *fastest* mode, and `dt` had not been recomputed.
"""
from __future__ import annotations

from dataclasses import dataclass

GATE = 1e-2          # threshold for every separation check. dt/tau=1e-2 <=> 0.5% bias
MARGIN_WARN = 5.0    # below this margin, warn "no room to raise a parameter"
SOFT_KINDS = frozenset({"statistics", "finite-size"})
# * Floating-point boundary tolerance. Prevents a check whose limit is set
#   **equal to the construction** from failing by accident.
#   This actually bit: chain-bend defines T_obs == 100*(2*pi/omega), and
#   `(100x)/x = 99.99999999999999` made "cycles observed >= 100" emit a false
#   warning. The same spec's `dt/tau_fast <= 1e-2` happened to divide exactly and
#   passed -- do not rely on that coincidence. Same intent as
#   `bdbot.nondim.RATIO_RTOL`.
CMP_RTOL = 1e-12


@dataclass
class Check:
    kind: str                    # model / integration / geometry / statistics
    name: str
    value: float
    limit: float
    op: str = "<="
    note: str = ""
    hard: bool = True

    @property
    def ok(self) -> bool:
        tol = CMP_RTOL * max(abs(self.value), abs(self.limit))
        return (self.value <= self.limit + tol if self.op == "<="
                else self.value >= self.limit - tol)

    @property
    def margin(self) -> float:
        if self.op == "<=":
            return self.limit / self.value if self.value else float("inf")
        return self.value / self.limit if self.limit else float("inf")

    def as_dict(self, phase: str = "design") -> dict:
        # `op` is needed when reading an L3 spec back -- without it a ">=" check
        # (an observation window, say) is restored as "<=" and the verdict flips
        # (bdbot.nondim.load).
        return {"kind": self.kind, "name": self.name.strip(), "value": self.value,
                "limit": self.limit, "op": self.op, "ok": bool(self.ok),
                "margin": self.margin, "hard": bool(self.hard), "phase": phase}


def soft(kind: str) -> bool:
    """Is this kind a warning? (the hard/soft definition in bd-physics section 4)"""
    return kind in SOFT_KINDS


def verdict(checks) -> tuple[str, list, list, list]:
    """(verdict, hard failures, soft failures, thin margins).

    A single broken hard check refuses the run.
    """
    hard_fail = [c for c in checks if c.hard and not c.ok]
    soft_fail = [c for c in checks if not c.hard and not c.ok]
    tight = [c for c in checks if c.ok and c.margin < MARGIN_WARN]
    if hard_fail:
        v = "FAIL"
    elif soft_fail:
        v = f"PASS ({len(soft_fail)} warnings)"
    else:
        v = "PASS"
    return v, hard_fail, soft_fail, tight


# -- dt selection: two cases independently arrived at "1% of gamma/(local stiffness)" --
def relaxation_time(gamma, stiffness):
    """tau = gamma/k. For a trap k is the trap stiffness; for a soft pair k = U''(r_min).

    The same formula appeared in two cases under different names (tau_k, tau_int) —
    it is one structure.
    """
    return (gamma / stiffness).to("s")


def dt_from_gate(tau, gate: float = GATE):
    """dt = gate*tau, matched to the hard gate. gate=1e-2 is 0.5% bias in a
    linear system.
    """
    return gate * tau


def dt_from_bias(tau, bias: float):
    """Invert from a target bias (bd-physics section 1.2):
    bias ~ (dt/tau)/2  =>  dt = 2*bias*tau.

    WARNING: the closed form holds **only for a linear system.** A nonlinear
    system needs a dt-halved convergence check.
    """
    return 2 * bias * tau


def bias_from_dt(dt, tau) -> float:
    """Expected systematic bias [%] from dt (linear system)."""
    return 50.0 * float((dt / tau).to("dimensionless").magnitude)


__all__ = ["Check", "GATE", "MARGIN_WARN", "SOFT_KINDS", "soft", "verdict",
           "relaxation_time", "dt_from_gate", "dt_from_bias", "bias_from_dt"]
