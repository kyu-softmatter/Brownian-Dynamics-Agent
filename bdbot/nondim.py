"""The L3 `NondimSpec` -- the artefact of non-dimensionalization.
skill `bd-physics` section 0, steps 3-5.

The **only contract** between L2 (`system.yaml`, SI) and L4 (execution). A run
must be possible from this spec alone, and a result must be convertible back into
physical units from this spec alone.

## Why it exists -- four measured defects (`verify/verify_l3_spec_gaps.py`)

Three cases had each been hand-building a `spec = {...}` dictionary with mutually
different schemas (only three keys in common). So:

  1. * **`run_id` did not cover the physical system.** One spec contained no
     physical system at all, so changing `d` 5um -> 0.5um, `eta` water -> glycerol
     (62x) and `rho_p` silica -> polystyrene left the `run_id` **exactly the
     same.** That id belongs to an already-completed run, so `prepare_outdir`
     skips with "this run is already complete" and **reports the old system's
     result as the new system's result** -- for a system whose tau_B differs by
     16.1x.
     -> The spec **must** include `system` (physics_only).
  2. Inversion was impossible from the spec alone (no SI values for sigma, tau,
     kT) -- violating bd-physics section 5.
     -> `back_transform` pins the three anchors as SI floats.
  3. There was no way to check that a dimensionless group really was the ratio of
     those two scales.
     -> `Group.num`/`den` point at ledger symbols and `validate()` recomputes and
     compares.
  4. Nobody could catch a scale missing from the ledger.
     -> `ScaleLedger.missing_roles()` is used as a hard gate.

## What was not built (it has not appeared even once -- the abstraction rule)

  * reference strategies other than `thermal` (`interaction`, `active`, `custom`)
  * the inverse construction `from_dimensionless(groups, anchors)` -- no case
    starts from dimensionless values
  * an engine that **computes the ledger for you** -- which scales exist differs
    per system, and that is the physics judgment deliberately left in the case
    scripts (bd-physics section 6.3)

Only `{system, params, numerics}` enter the `run_id` hash -- `schema`, the ledger,
the checks, the bases and the derived values are excluded. If bumping a schema
version or editing a comment invalidates a run, content addressing is useless
(the lesson in `bdbot.runid`).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import intake as _intake
from . import runid as _runid
from . import scales as _scales
from .checks import Check
from .checks import verdict as _verdict
from .units import Q

SCHEMA = "bdbot.nondim/0.1"
RATIO_RTOL = 1e-9          # a group is the ratio of two ledger entries. Float error only.


# ══════════════════════════════════════════════════════════════════════
# 2. reference scales
# ══════════════════════════════════════════════════════════════════════
@dataclass
class Reference:
    """The length, time and energy chosen as references, plus the basis.

    WARNING: a different reference gives the same system entirely different
       dimensionless groups. Comparing against the literature means checking which
       one they used, so `strategy` and `rationale` are stored alongside.
    """

    length: tuple           # (symbol, Quantity)
    time: tuple
    energy: tuple
    strategy: str = "thermal"
    rationale: str = ""

    @classmethod
    def from_dict(cls, d: dict) -> "Reference":
        """Accept `scales.thermal_reference()`'s return value directly (case code unchanged)."""
        return cls(length=tuple(d["length"]), time=tuple(d["time"]),
                   energy=tuple(d["energy"]), strategy=d.get("strategy", "thermal"),
                   rationale=d.get("rationale", ""))

    def as_dict(self) -> dict:
        """Convert back to the dict shape `bdbot.report.render` expects."""
        return {"length": self.length, "time": self.time, "energy": self.energy,
                "strategy": self.strategy, "rationale": self.rationale}

    def si(self, kind: str):
        return {"length": self.length, "time": self.time, "energy": self.energy}[kind][1]

    def to_json(self) -> dict:
        out = {"strategy": self.strategy, "rationale": self.rationale}
        for kind in ("length", "time", "energy"):
            sym, q = getattr(self, kind)
            qb = q.to_base_units()
            out[kind] = {"symbol": sym, "value": float(qb.magnitude), "unit": str(qb.units)}
        return out

    @classmethod
    def from_json(cls, d: dict) -> "Reference":
        g = lambda k: (d[k]["symbol"], Q(d[k]["value"], d[k]["unit"]))
        return cls(length=g("length"), time=g("time"), energy=g("energy"),
                   strategy=d.get("strategy", "thermal"), rationale=d.get("rationale", ""))


# ══════════════════════════════════════════════════════════════════════
# 3. a dimensionless group is a ratio of two scales
# ══════════════════════════════════════════════════════════════════════
@dataclass
class Group:
    """One dimensionless group.

    Given `num`/`den`, it is **recomputed from the ledger and compared.** That is
    the reason this class exists -- a dict holding only numbers cannot catch "the
    name is dt/tau_int but the value is something else."

        Group("dt/tau_int", 0.01, num=("times", "dt"), den=("times", "tau_int"),
              meaning="integration resolution")

    Leave `num`/`den` off for anything that is not a ratio (phi, A and other
    inputs that are dimensionless by definition).
    """

    name: str
    value: float
    num: tuple | None = None        # (category, symbol)
    den: tuple | None = None
    expr: str = ""
    meaning: str = ""

    @property
    def label(self) -> str:
        """The name for one report line -- shaped `name  = expression   meaning`."""
        head = f"{self.name:<12}" + (f"= {self.expr:<14}" if self.expr else " " * 16)
        return (head + self.meaning).rstrip()

    def recompute(self, ledger) -> float | None:
        """Recompute the ratio from the ledger. None when num/den are absent."""
        if not (self.num and self.den):
            return None
        nc, ns = self.num
        dc, ds = self.den
        return float((ledger.get(nc, ns) / ledger.get(dc, ds)).to("dimensionless").magnitude)

    def mismatch(self, ledger) -> str | None:
        """A mismatch message, or None. A KeyError (a symbol absent from the ledger)
        counts as a defect too.
        """
        try:
            want = self.recompute(ledger)
        except KeyError as e:
            return f"points at a symbol absent from the ledger: {e}"
        except Exception as e:                      # a dimensionally inconsistent ratio -> not dimensionless
            return f"cannot compute the ratio ({e})"
        if want is None:
            return None
        if abs(self.value - want) > RATIO_RTOL * max(abs(want), 1e-300):
            return (f"written {self.value:.12g} != ledger recomputation "
                    f"{self.num[1]}/{self.den[1]} = {want:.12g}")
        return None

    def to_json(self) -> dict:
        return {"name": self.name, "value": self.value, "expr": self.expr,
                "meaning": self.meaning,
                "num": list(self.num) if self.num else None,
                "den": list(self.den) if self.den else None}

    @classmethod
    def from_json(cls, d: dict) -> "Group":
        return cls(name=d["name"], value=d["value"], expr=d.get("expr", ""),
                   meaning=d.get("meaning", ""),
                   num=tuple(d["num"]) if d.get("num") else None,
                   den=tuple(d["den"]) if d.get("den") else None)


def groups_dict(groups) -> dict:
    """As the {display name: value} that `bdbot.report.render` takes -- **for the
    human-readable report.**
    """
    return {g.label: g.value for g in groups}


def metrics_dict(groups) -> dict:
    """For `metrics.json`'s `dimensionless` -- **machine-readable keys.**

    These used to be the report's display names verbatim
    (`'k*     = k d^2/kT   trap vs thermal'`). metrics.json is the post-mortem's
    only input, and such a key cannot be queried. Symbols only (`'k*'`,
    `'dt/tau_k'`).
    """
    return {g.name: g.value for g in groups}


# ══════════════════════════════════════════════════════════════════════
# the L3 artefact
# ══════════════════════════════════════════════════════════════════════
@dataclass
class NondimSpec:
    """The complete non-dimensionalization result. Stored as
    `specs/<run_id>.json`, and L4 reads only this.
    """

    case: str
    system: dict                                  # the L2 original (before physics_only)
    reference: Reference
    ledger: Any                                   # ScaleLedger
    groups: list = field(default_factory=list)
    checks: list = field(default_factory=list)
    numerics: dict = field(default_factory=dict)  # dimensionless run parameters (dt_star, n_prod, seed, ...)
    params: dict = field(default_factory=dict)    # per-case dimensionless knobs (A, phi, r_c_star, ...)
    tag: str | None = None
    nhex: int = 12                                # 12 hex chars in one case, 10 in the other (kept for reproducibility)
    label: str | None = None

    def __post_init__(self):
        if self.label is None:
            self.label = self.system.get("label", self.case)
        if isinstance(self.reference, dict):       # accept thermal_reference() directly
            self.reference = Reference.from_dict(self.reference)

    # -- content addressing ---------------------------------------------------
    def hash_payload(self) -> dict:
        """The target of the `run_id` hash -- **only what fixes the physics.**

        `schema`, the ledger, the checks, the bases and the derived values are left
        out. If bumping a schema version or editing a comment invalidates a run,
        content addressing is useless (the `bdbot.runid` DOC_KEYS lesson).
        """
        return {"system": _runid.physics_only(self.system),
                "params": _runid.physics_only(self.params),
                "numerics": _runid.physics_only(self.numerics)}

    def run_id(self) -> str:
        return _runid.content_run_id(self.label, self.hash_payload(),
                                     tag=self.tag, nhex=self.nhex)

    # -- 4. verdict -----------------------------------------------------------
    def verdict(self) -> str:
        return _verdict(self.checks)[0]

    def validate(self) -> list:
        """L3's own integrity check -- **a different layer** from the physics
        checks (`checks`).

        What it looks at: is the ledger complete; is each group really that ratio;
        is the reference actually in the ledger; do the inversion anchors hold. That
        is "was the non-dimensionalization done correctly", whereas `checks` asks
        "is BD valid for this system and integrated finely enough".
        """
        I = _intake.Issue
        out = []

        # 1. ledger completeness -- a missing scale is a check that never runs
        for role in self.ledger.missing_roles():
            out.append(I("error", f"ledger.{role}",
                         f"required role '{role}' ({_scales.ROLE_MEANING.get(role, '')}) "
                         f"is not in the ledger. Add it, or state it explicitly with "
                         f"`declare_absent('{role}', reason)`."))

        # 2. is each group the ratio of two ledger entries
        n_checked = 0
        for g in self.groups:
            msg = g.mismatch(self.ledger)
            if msg:
                out.append(I("error", f"groups.{g.name}", msg))
            elif g.num and g.den:
                n_checked += 1
        if self.groups and n_checked == 0:
            out.append(I("warn", "groups",
                         "not one dimensionless group was compared against the ledger "
                         "-- every num/den is empty. Attach ledger symbols to anything "
                         "that is a ratio."))

        # 3. is the reference scale actually in the ledger (a reference that is not
        #    in the ledger cannot be compared against anything)
        for kind, cat in (("length", "lengths"), ("time", "times"), ("energy", "energies")):
            sym = getattr(self.reference, kind)[0]
            if not self.ledger.has(cat, sym):
                out.append(I("error", f"reference.{kind}",
                             f"reference '{sym}' is not in ledger {cat}."))

        # 4. are the inversion anchors positive and finite (0 cannot be inverted)
        for kind in ("length", "time", "energy"):
            q = self.reference.si(kind)
            try:
                v = float(q.to_base_units().magnitude)
            except Exception as e:
                out.append(I("error", f"reference.{kind}", f"SI conversion failed: {e}"))
                continue
            if not (v > 0) or v != v or v in (float("inf"), float("-inf")):
                out.append(I("error", f"reference.{kind}",
                             f"the reference scale is not positive and finite ({v!r}) -- cannot invert."))

        # 5. required run parameters (L4 runs from these alone)
        for k in ("dt_star", "n_prod"):
            if k not in self.numerics:
                out.append(I("error", f"numerics.{k}",
                             "L4 cannot run from the spec alone (a dimensionless run parameter is missing)."))
        dt_star = self.numerics.get("dt_star")
        if isinstance(dt_star, (int, float)) and not dt_star > 0:
            out.append(I("error", "numerics.dt_star", f"dt* is not positive ({dt_star!r})."))

        # 6. `dt_star` **must equal** the ledger's dt / reference time. The two are
        #    computed separately, so they can diverge, and then HOOMD runs at a
        #    different step from the ledger -- the separation checks pass on the
        #    ledger's dt while the actual integration uses another value. Silently
        #    wrong.
        if isinstance(dt_star, (int, float)) and self.ledger.has("times", "dt"):
            try:
                want = float((self.ledger.get("times", "dt")
                              / self.reference.si("time")).to("dimensionless").magnitude)
            except Exception as e:
                out.append(I("error", "numerics.dt_star", f"comparison against the ledger dt failed: {e}"))
            else:
                if abs(dt_star - want) > RATIO_RTOL * max(abs(want), 1e-300):
                    out.append(I("error", "numerics.dt_star",
                                 f"disagrees with the ledger: numerics.dt_star = {dt_star:.12g} != "
                                 f"dt/{self.reference.time[0]} = {want:.12g}. "
                                 "HOOMD would run at a different step from the ledger."))
        return out

    @property
    def errors(self) -> list:
        return [i for i in self.validate() if i.level == "error"]

    # -- 5. inversion (bd-physics section 5 -- always performed) ---------------
    def physical(self, value, L: int = 0, T: int = 0, E: int = 0):
        """A dimensionless value -> physical units. Specified by dimension exponents.

            spec.physical(D_star, L=2, T=-1)     # D_eff → m²/s
            spec.physical(x2_star, L=2)          # ⟨x²⟩  → m²
            spec.physical(step * dt_star, T=1)   # time      -> s
            spec.physical(P_star, E=1, L=-2)     # 2D pressure -> N/m
        """
        sigma = self.reference.si("length")
        tau = self.reference.si("time")
        en = self.reference.si("energy")
        return value * sigma**L * tau**T * en**E

    def back_transform(self) -> dict:
        """The three inversion anchors as SI floats. With the spec alone, inversion
        works without pint.
        """
        out = {}
        for kind, key in (("length", "sigma_SI"), ("time", "tau_SI"), ("energy", "energy_SI")):
            qb = self.reference.si(kind).to_base_units()
            out[key] = float(qb.magnitude)
            out[key + "_unit"] = str(qb.units)
        return out

    # -- serialization --------------------------------------------------------
    def _ledger_json(self) -> dict:
        """Store the ledger both in SI and dimensionless (divided by the
        reference). L4 uses the starred side.
        """
        out = {}
        ref_of = {"lengths": "length", "times": "time", "energies": "energy"}
        for cat_name, cat in self.ledger.categories():
            rows = []
            base = self.reference.si(ref_of[cat_name])
            for sym, sc in cat.items():
                q = sc.value if isinstance(sc, _scales.Scale) else sc
                qb = q.to_base_units()
                row = {"symbol": sym, "value": float(qb.magnitude), "unit": str(qb.units)}
                if isinstance(sc, _scales.Scale):
                    row.update(note=sc.note, role=sc.role, star=sc.star)
                try:
                    row["reduced"] = float((q / base).to("dimensionless").magnitude)
                except Exception:
                    row["reduced"] = None       # an entry of a different kind (should not happen)
                rows.append(row)
            out[cat_name] = rows
        return out

    def to_json(self) -> dict:
        return {
            "schema": SCHEMA,
            "case": self.case,
            "label": self.label,
            "tag": self.tag,
            "run_id": self.run_id(),
            "system": self.system,
            "reference": self.reference.to_json(),
            "back_transform": self.back_transform(),
            "ledger": self._ledger_json(),
            "ledger_absent": dict(self.ledger.absent),
            "groups": [g.to_json() for g in self.groups],
            "checks": [c.as_dict() for c in self.checks],
            "verdict": self.verdict(),
            "params": self.params,
            "numerics": self.numerics,
            "l3_issues": [{"level": i.level, "where": i.where, "msg": i.msg}
                          for i in self.validate()],
        }

    def write(self, path) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.to_json(), indent=2, ensure_ascii=False, default=str))
        return p


# ══════════════════════════════════════════════════════════════════════
# reading -- L4 uses only this function (it does not import case scripts)
# ══════════════════════════════════════════════════════════════════════
@dataclass
class LoadedSpec:
    """A `specs/*.json` read back. The ledger is restored as Quantities so
    inversion works.

    Why a `NondimSpec` is not restored verbatim: a stored spec has no
    `ScaleLedger.derived` (the case's intermediates) and does not need one. All L4
    needs is the dimensionless parameters and the inversion anchors.
    """

    raw: dict
    reference: Reference
    groups: list
    checks: list

    @property
    def case(self) -> str:
        return self.raw["case"]

    @property
    def label(self) -> str:
        return self.raw.get("label", self.raw["case"])

    @property
    def run_id(self) -> str:
        return self.raw["run_id"]

    @property
    def params(self) -> dict:
        return self.raw.get("params", {})

    @property
    def numerics(self) -> dict:
        return self.raw.get("numerics", {})

    @property
    def verdict(self) -> str:
        return self.raw.get("verdict", "?")

    def reduced(self, cat: str, symbol: str) -> float:
        """A dimensionless value (divided by the reference). L4 gets L*, r_c*, dt*
        from here.
        """
        for row in self.raw["ledger"][cat]:
            if row["symbol"] == symbol:
                return row["reduced"]
        raise KeyError(f"{cat}.{symbol} is not in the spec ledger")

    def si(self, cat: str, symbol: str):
        for row in self.raw["ledger"][cat]:
            if row["symbol"] == symbol:
                return Q(row["value"], row["unit"])
        raise KeyError(f"{cat}.{symbol} is not in the spec ledger")

    def group(self, name: str) -> float:
        for g in self.groups:
            if g.name == name:
                return g.value
        raise KeyError(f"dimensionless group '{name}' is not in the spec")

    def physical(self, value, L: int = 0, T: int = 0, E: int = 0):
        bt = self.raw["back_transform"]
        sigma = Q(bt["sigma_SI"], bt["sigma_SI_unit"])
        tau = Q(bt["tau_SI"], bt["tau_SI_unit"])
        en = Q(bt["energy_SI"], bt["energy_SI_unit"])
        return value * sigma**L * tau**T * en**E

    def render(self) -> str:
        """Redraw the report from **the spec alone** -- with no case script.

        If this works the spec is self-sufficient and L4 need not import case code.
        If it does not, something is missing from the spec.
        """
        r = self.raw
        L: list[str] = []
        w = L.append
        W = 88
        w("=" * W)
        w(f"NondimSpec — {self.label}   run_id={self.run_id}")
        w("=" * W)
        ref = r["reference"]
        w(f"reference scales: length={ref['length']['symbol']}  energy={ref['energy']['symbol']}"
          f"  time={ref['time']['symbol']}   [strategy: {ref['strategy']}]")
        w(f"  basis: {ref.get('rationale', '')}")

        ok, want = self.verify_hash()
        w("")
        w(f"hash self-check: {'ok -- run_id matches the spec contents' if ok else f'MISMATCH -- expected {want}'}")

        w("")
        w("SCALE LEDGER  (SI . reduced against the reference)")
        for cat in ("lengths", "times", "energies"):
            rows = sorted(r["ledger"].get(cat, []), key=lambda x: x["value"])
            if not rows:
                continue
            w(f"  {cat}")
            for row in rows:
                star = " ★" if row.get("star") else ""
                role = f"  [{row['role']}]" if row.get("role") else ""
                red = row.get("reduced")
                red_s = f"{red:.6g}" if isinstance(red, (int, float)) else "—"
                w(f"    {row['symbol']:<11}{row['value']:>13.5e} {row['unit']:<22}"
                  f"{red_s:>13}{role}{star}")
        if r.get("ledger_absent"):
            for role, why in r["ledger_absent"].items():
                w(f"    (absent) {role}: {why}")

        w("")
        w("DIMENSIONLESS GROUPS")
        for g in self.groups:
            ratio = f"{g.num[1]}/{g.den[1]}" if (g.num and g.den) else "-- (not a ratio)"
            w(f"  {g.name:<16}{g.value:>13.6g}   {ratio:<22}{g.meaning}")

        w("")
        w(f"{'SEPARATION CHECKS':<58}{'value':>10}{'limit':>10}{'margin':>9}")
        for c in self.checks:
            mark = "✓" if c.ok else ("✗" if c.hard else "⚠")
            w(f"  {mark} [{c.kind}] {c.name:<42}{c.value:>10.3e}{c.limit:>10.0e}"
              f"{c.margin:>8.1f}×")

        w("")
        w("BACK TRANSFORM  (result -> physical units)")
        bt = r["back_transform"]
        w(f"  σ = {bt['sigma_SI']:.6e} {bt['sigma_SI_unit']}    "
          f"τ = {bt['tau_SI']:.6e} {bt['tau_SI_unit']}")
        w(f"  E = {bt['energy_SI']:.6e} {bt['energy_SI_unit']}")

        w("")
        w("RUN PARAMETERS  (what L4 reads)")
        for k, v in r.get("params", {}).items():
            w(f"  params.{k:<20} {v}")
        for k, v in r.get("numerics", {}).items():
            w(f"  numerics.{k:<18} {v}")

        issues = r.get("l3_issues", [])
        if issues:
            w("")
            w("L3 INTEGRITY")
            for i in issues:
                mark = {"error": "✗", "warn": "⚠", "info": "ℹ"}.get(i["level"], "·")
                w(f"  {mark} [{i['where']}] {i['msg']}")

        w("")
        w(f"VERDICT: {self.verdict}")
        w("=" * W)
        return "\n".join(L)

    def verify_hash(self) -> tuple[bool, str]:
        """Does the stored `run_id` match the spec contents -- catches a hand-edited
        spec.

        This is where a machine confirms rule 2 ("never hand-write a dimensionless
        spec").
        """
        payload = {"system": _runid.physics_only(self.raw.get("system", {})),
                   "params": _runid.physics_only(self.params),
                   "numerics": _runid.physics_only(self.numerics)}
        stored = self.raw["run_id"]
        nhex = len(stored.rsplit("__", 1)[-1])
        want = _runid.content_run_id(self.label, payload,
                                     tag=self.raw.get("tag"), nhex=nhex)
        return want == stored, want


def load(path) -> LoadedSpec:
    raw = json.loads(Path(path).read_text())
    got = raw.get("schema")
    if got != SCHEMA:
        raise ValueError(f"schema differs: {got!r} (expected {SCHEMA!r})")
    return LoadedSpec(
        raw=raw,
        reference=Reference.from_json(raw["reference"]),
        groups=[Group.from_json(g) for g in raw.get("groups", [])],
        checks=[Check(kind=c["kind"], name=c["name"], value=c["value"], limit=c["limit"],
                      op=c.get("op", "<="), note=c.get("note", ""), hard=c.get("hard", True))
                for c in raw.get("checks", [])],
    )


__all__ = ["SCHEMA", "Reference", "Group", "NondimSpec", "LoadedSpec", "load",
           "groups_dict", "RATIO_RTOL"]
