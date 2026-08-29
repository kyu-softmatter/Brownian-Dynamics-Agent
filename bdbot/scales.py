"""The scale ledger -- skill `bd-physics` section 0.

Enumerate every length, time and energy in the system in SI, then choose the
references **with a stated basis.** Fixing the dimensionless groups first and
inferring the scales afterwards is forbidden.

Both of the first two cases had the same structure
(lengths/times/energies + ref + rationale). **What differs per case is which
scales go in**, and the case fills those.

## What the L3 work changed

A ledger entry's key used to be a **single blob with the symbol and the
description glued together**, like `"d        particle diameter (reference)"`.
That was enough to render a human-readable table but a machine could not use it:

  * using `ratio("times", a, b)` meant matching the decoration exactly.
  * there was **no way to check** that a dimensionless group really was the ratio
    of those two scales -- you could write `dt/tau_int` and put a different value
    in and nobody would catch it. This project actually evaluated the cage
    stiffness at `a_mean` when it should have been `a_NN`, a 41% error, and there
    was no place that would catch that class of mistake.
  * nobody could catch **a scale missing from the ledger.** A missing scale is
    **a check that never runs** -- and "silently passing" is different from "not
    checking."

So it was split into `Scale(symbol, value, note, role)`. `role` names "the scale
that plays the box / time-step / observation-window / inertia role in this
system", independently of what the symbol is called -- an ellipsoid case names its
reference length `d_eq`, so required entries cannot be enforced by symbol.

A scale that does not exist is emptied **explicitly** with
`declare_absent(role, reason)` -- the same attitude as rule 3: do not invent what
you do not know; record that it is absent.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

THERMAL_RATIONALE = (
    "thermal convention, fixed (sigma=d, E=kT, tau=tau_B). The unit system is not "
    "changed per system, so that literature comparison stays possible "
    "(bd-physics 1.1). The governing timescale may differ per system, and dt is "
    "set from that one."
)

# Roles that must be present in the ledger. If absent, `declare_absent` must
# record the reason.
#   box         the denominator of the geometry checks (minimum image, finite size)
#   dt          the numerator of the integration-resolution check -- if it is not
#               in the ledger it is invisible in the timescale ordering
#   observation the numerator of the statistics check (T_obs/tau)
#   inertia     the model-validity check (tau_p/tau_dyn) -- is BD admissible here
MANDATORY_ROLES = ("box", "dt", "observation", "inertia")
ROLE_MEANING = {
    "box": "box edge L -- the denominator of the geometry checks",
    "dt": "integration time step",
    "observation": "observation window T_obs",
    "inertia": "momentum relaxation tau_p -- for the model-validity verdict",
    "ref_length": "reference length", "ref_time": "reference time",
    "ref_energy": "reference energy",
}


@dataclass
class Scale:
    """One ledger entry. The symbol and the description are **separated** so a
    machine can address it by symbol.
    """

    symbol: str
    value: Any                   # pint Quantity (SI)
    note: str = ""
    role: str = ""               # one of MANDATORY_ROLES, or the empty string
    star: bool = False           # highlight in the report (this system's governing scale)

    @property
    def display(self) -> str:
        """The name for one report line -- reproduces the old key shape.

        WARNING: the padding is a **minimum width.** When `chain-bend` used a
           12-character symbol like `kappa_end_d2`, a fixed 9 characters glued the
           symbol and the description together. Display only, so it affects
           neither the ledger values nor the run_id.
        """
        s = f"{self.symbol:<12} {self.note}"
        return (s + " ★") if self.star else s


def _qty(v):
    """Extract the value whether it is a Scale or a Quantity (the ledger accepts
    both).
    """
    return v.value if isinstance(v, Scale) else v


@dataclass
class ScaleLedger:
    """The length / time / energy ledger, plus the reference scales and the basis
    for choosing them.
    """

    lengths: dict = field(default_factory=dict)     # symbol → Scale
    times: dict = field(default_factory=dict)
    energies: dict = field(default_factory=dict)
    derived: dict = field(default_factory=dict)     # intermediates the case uses (not the ledger)
    ref: dict = field(default_factory=dict)
    rationale: str = ""
    absent: dict = field(default_factory=dict)      # role -> the reason it is absent

    # -- writing -------------------------------------------------------------
    def add(self, cat: str, symbol: str, value, note: str = "", role: str = "",
            star: bool = False) -> "ScaleLedger":
        """Add one ledger entry. `cat` is "lengths" / "times" / "energies"."""
        getattr(self, cat)[symbol] = Scale(symbol, value, note, role, star)
        return self

    def add_length(self, symbol, value, note="", role="", star=False):
        return self.add("lengths", symbol, value, note, role, star)

    def add_time(self, symbol, value, note="", role="", star=False):
        return self.add("times", symbol, value, note, role, star)

    def add_energy(self, symbol, value, note="", role="", star=False):
        return self.add("energies", symbol, value, note, role, star)

    def declare_absent(self, role: str, reason: str) -> "ScaleLedger":
        """Declare explicitly that this system has **no** scale in that role.

        It cannot be emptied without a reason.
        """
        if not reason:
            raise ValueError(f"emptying role '{role}' requires a reason (rule 3)")
        self.absent[role] = reason
        return self

    # -- reading -------------------------------------------------------------
    def categories(self):
        return (("lengths", self.lengths), ("times", self.times), ("energies", self.energies))

    def sorted_items(self, cat: dict):
        """(display name, value), smallest first. Sorted, a separation violation
        becomes visible.
        """
        items = sorted(cat.items(), key=lambda kv: _qty(kv[1]).to_base_units().magnitude)
        return [((v.display if isinstance(v, Scale) else k), _qty(v)) for k, v in items]

    def get(self, cat: str, symbol: str):
        """One value by symbol. KeyError if absent -- never a silent None."""
        return _qty(getattr(self, cat)[symbol])

    def has(self, cat: str, symbol: str) -> bool:
        return symbol in getattr(self, cat)

    def by_role(self, role: str):
        """Look up by role. (symbol, value) or None."""
        for _, cat in self.categories():
            for sym, sc in cat.items():
                if isinstance(sc, Scale) and sc.role == role:
                    return sym, sc.value
        return None

    def ratio(self, cat: str, a: str, b: str) -> float:
        """The ratio of two scales -- every dimensionless group must have this
        form (bd-physics section 3).
        """
        dd = getattr(self, cat)
        return float((_qty(dd[a]) / _qty(dd[b])).to("dimensionless").magnitude)

    def span(self, cat: str) -> float:
        """Largest scale / smallest scale. How many decades the system spans."""
        dd = getattr(self, cat)
        if len(dd) < 2:
            return 1.0
        vs = [_qty(v).to_base_units().magnitude for v in dd.values()]
        return max(vs) / min(vs)

    # -- completeness --------------------------------------------------------
    def missing_roles(self) -> list[str]:
        """Required roles that are neither in the ledger nor `declare_absent`ed.

        **A missing scale is a check that never runs.** One case had been passing
        `dt` and `T_obs` separately instead of putting them in the ledger -- and
        neither was visible in the timescale ordering table.
        """
        return [r for r in MANDATORY_ROLES
                if self.by_role(r) is None and r not in self.absent]


def thermal_reference(d, kT, tau_B, rationale: str | None = None, *,
                      length_symbol: str = "d", energy_symbol: str = "kT",
                      time_symbol: str = "tau_B") -> dict:
    """The `thermal` reference triple. In the HOOMD spec this becomes
    sigma=1, kT=1, gamma=1, tau_B=1.

    Why the symbol names are overridable: an ellipsoid case's reference length is
    the equivalent-volume sphere diameter and the ledger calls it `d_eq`. Pinning
    the symbol to "d" would make the report point at a reference that is not in the
    ledger, and `NondimSpec.validate()` raises that as an error (it actually did).
    """
    return {"length": (length_symbol, d), "energy": (energy_symbol, kT),
            "time": (time_symbol, tau_B),
            "strategy": "thermal", "rationale": rationale or THERMAL_RATIONALE}


__all__ = ["ScaleLedger", "Scale", "thermal_reference", "THERMAL_RATIONALE",
           "MANDATORY_ROLES", "ROLE_MEANING"]
