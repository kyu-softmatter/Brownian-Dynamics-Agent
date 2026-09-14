"""Rule 10 — the parameter manifest. **Every number, before the run.** LLM 0 lines.

> *"Every time you formulate a problem, lay out all the numbers you are using —
> particle radius, box size, kT, trap stiffness if applicable, timestep and total
> steps, and so on — before you start the run, and ask for permission."*
> — user instruction, 2026-09-02, now CLAUDE.md rule 10

## Why a manifest rather than a habit

This repository already had the numbers. `trap-2d-5um`'s report prints `d`, `T`,
`eta`, `k_t`, `rho_p`, `N`, `gamma`, `D_t`, `tau_B`, `dt`, `T_obs` and the box.
What it did **not** have is anything that fails when one of them is missing, so
"lay out the numbers" was a practice rather than a gate — and practices in this
repository have a measured record of not holding:

  - `bd-intake` §2.1's empty-goal blocker: written twice, enforced zero times,
    walked past by 2 of 8 cases which then produced 85 runs
  - `A4`'s grep: 7 false hits, never a real check, unenforced for a month
  - `health.gate()`: reachable only from `tools/health.py:138`, absent from
    `run.py` entirely, so no run has ever gated itself

Three for three. So this is a dataclass with a completeness check, not a
docstring.

## What it enforces

**1 · Every number carries `value`, `unit` and `provenance`.** A bare float is
rejected. `unit="1"` is how a dimensionless quantity says so, because blank must
mean *forgotten* rather than *dimensionless*.

**2 · Every required category is present or explicitly not applicable.**
`not_applicable(category, reason)` mirrors `scales.declare_absent(role, reason)`
and `design.Abandoned`: it raises on an empty reason. "If applicable" in the
instruction is exactly this — a trap stiffness is genuinely absent for a free
probe, and *saying so* is different from omitting it.

**3 · A human approves before the run.** `approved_by` is set by a person and by
nothing else. `blockers()` reports it unset, `runcard` refuses to seal, and
`run.execute(require_approval=True)` refuses to run.

## The categories, and where each number comes from

    geometry     d or radius, box L, dimensionality, N
    energy       T, kT
    medium       eta, rho_fluid
    interaction  trap k, pair potential parameters, bond/angle  (may be N/A)
    numerics     dt, n_eq, n_prod, sample_every, seed
    derived      gamma, D_t, tau_B, and the GOVERNING timescale

`derived` is listed separately and checked against a recomputation, because
`physical.verify()` already catches a hand-edited derived value and this is the
same invariant one layer up.

## ★ Four for four (2026-09-14)

The paragraph above was **false for twelve days**. `derived` was listed
separately and nothing recomputed it: `blockers()` checked that each
`REQUIRED_KEYS` name was *present*, never that its value was *real*. Measured by
using it -- the first campaign to build a manifest wrote

    m.add("derived", "gamma",   0.0, "1", "6 pi eta a -- in the ledger, recomputed by L3")
    m.add("derived", "D_t",     0.0, "1", "kT/gamma -- in the ledger, recomputed by L3")
    m.add("derived", "tau_gov", 0.0, "1", "tau_B = d^2/D_t -- in the ledger, recomputed by L3")

and `blockers()` returned `[]`. Three physically impossible zeros -- no friction,
no diffusion, no timescale -- cleared the gate whose entire purpose is to refuse
a number that is not there, in the run written to demonstrate that gate. The
provenance string was itself the confession ("in the ledger") and nothing read
it.

So the list of three at the top of this docstring was really a list of four, and
the fourth was in this file. `check_derived()` below is the check the paragraph
had been claiming. `unknown()` is what "it is in the ledger" should have used.
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path

SCHEMA = "bdbot.params/0.1"

#: Required categories. A run card is not sealable until each is present or
#: explicitly declared not applicable with a reason.
CATEGORIES = ("geometry", "energy", "medium", "interaction", "numerics", "derived")

#: The minimum a category must name, by category. Not exhaustive -- a case may add
#: more -- but these cannot be missing, because each has been the cause of a
#: measured error in this repository.
REQUIRED_KEYS = {
    "geometry": ("d", "L", "dim", "N"),
    "energy": ("T", "kT"),
    "medium": ("eta",),
    "interaction": (),                 # per system; may be entirely N/A
    "numerics": ("dt", "n_prod", "seed"),
    "derived": ("gamma", "D_t", "tau_gov"),
}

#: Unit strings that mean "expressed in this system's own reference scales".
#: In that convention the scales are chosen so that kT = D = gamma = 1, which
#: makes every `derived` entry exactly 1 -- checkable without an SI
#: recomputation. `"1"` is included because rule 10 requires a dimensionless
#: quantity to say so rather than leave the field blank.
REDUCED_UNITS = ("1", "tau_B", "tau_trap", "tau_int", "tau_bond", "tau_v",
                 "d", "sigma", "kT")

#: Provenance tiers, same scale as `provenance.py`.
#:   0 given/handbook · 1 literature+verified or a confirmed convention
#:   2 literature unverified · 3 arbitrary assumption
TIERS = (0, 1, 2, 3)


@dataclass(frozen=True)
class Param:
    """One number. `value`, `unit` and `provenance` are all mandatory."""
    name: str
    value: float
    unit: str
    provenance: str
    tier: int = 3
    note: str = ""

    def __post_init__(self):
        if not str(self.name).strip():
            raise ValueError("a parameter needs a name")
        if self.value is None:
            raise ValueError(
                f"{self.name!r} has no value. If it is unknown, say so with "
                f"Manifest.unknown() -- a null that looks like a number is how an "
                f"invented value gets in (rule 3).")
        if not str(self.unit).strip():
            raise ValueError(
                f"{self.name!r} has no unit. Use '1' for a dimensionless "
                f"quantity: blank must mean forgotten, not dimensionless.")
        if not str(self.provenance).strip():
            raise ValueError(
                f"{self.name!r} has no provenance. Sketch, literature, handbook "
                f"or estimate -- say which (rule 3).")
        if self.tier not in TIERS:
            raise ValueError(f"{self.name!r}: tier must be one of {TIERS}, "
                             f"got {self.tier}")

    def line(self, w: int = 14) -> str:
        u = "" if self.unit == "1" else f" {self.unit}"
        v = (f"{self.value:.6g}" if isinstance(self.value, (int, float))
             else str(self.value))
        return (f"{self.name:<{w}} = {v}{u}".ljust(38)
                + f"[tier {self.tier}] {self.provenance[:34]}")


@dataclass
class Manifest:
    """Every number for one run, by category, plus the human approval."""
    case: str
    params: dict = field(default_factory=dict)      # category -> {name: Param}
    absent: dict = field(default_factory=dict)      # category -> reason
    unknowns: dict = field(default_factory=dict)    # name -> what is needed
    approved_by: str | None = None

    # -- construction -------------------------------------------------------
    def add(self, category: str, name: str, value, unit: str, provenance: str,
            tier: int = 3, note: str = ""):
        if category not in CATEGORIES:
            raise ValueError(f"category must be one of {CATEGORIES}, got {category!r}")
        if category in self.absent:
            raise ValueError(
                f"category {category!r} was declared not applicable "
                f"({self.absent[category]!r}); it cannot also carry {name!r}")
        self.params.setdefault(category, {})[name] = Param(
            name, value, unit, provenance, tier, note)
        return self

    def not_applicable(self, category: str, reason: str):
        """Declare a whole category absent. **Requires a reason.**

        Same rule as `scales.declare_absent` and `design.Abandoned`: "if
        applicable" means the absence is stated, not that the line is omitted.
        """
        if category not in CATEGORIES:
            raise ValueError(f"category must be one of {CATEGORIES}, got {category!r}")
        if not str(reason).strip():
            raise ValueError(
                f"declaring {category!r} not applicable requires a reason -- "
                f"otherwise it is indistinguishable from forgetting it (rule 3)")
        if self.params.get(category):
            raise ValueError(
                f"category {category!r} already carries "
                f"{sorted(self.params[category])}; it cannot also be N/A")
        self.absent[category] = str(reason)
        return self

    def unknown(self, name: str, needed_from: str):
        """Record a number that is **not known**. Blocks, and names the supplier.

        `BLOCKED` is a success in this repository -- it narrows the missing input
        to one thing. Inventing a value to avoid it is the failure (rule 3).
        """
        if not str(needed_from).strip():
            raise ValueError(f"{name!r}: say who or what would supply it")
        self.unknowns[str(name)] = str(needed_from)
        return self

    # -- reading ------------------------------------------------------------
    def get(self, category: str, name: str):
        return (self.params.get(category) or {}).get(name)

    def all_params(self) -> list:
        return [p for c in CATEGORIES for p in (self.params.get(c) or {}).values()]

    def tier_counts(self) -> dict:
        out = {}
        for p in self.all_params():
            out[p.tier] = out.get(p.tier, 0) + 1
        return dict(sorted(out.items()))

    #: How close a declared `derived` value must sit to its recomputation.
    #:
    #: ⚠ Was 1e-6 for about an hour, and 1e-6 REFUSED A CORRECT MANIFEST.
    #:   `trap-2d-5um`'s card writes `gamma = 4.0102e-8 kg/s`, which is the right
    #:   number to five significant figures and sits 1.07e-5 from the
    #:   recomputation. The tolerance of a comparison comes from how the value
    #:   was WRITTEN, not from how precisely it could be computed -- the same
    #:   lesson this repository already recorded for a 5-s.f. finite-size number.
    #:   0.1 % accepts anything written to four figures and still catches what
    #:   this check is for: a placeholder, a wrong unit, an order-of-magnitude
    #:   slip, a stale value. The tight comparison belongs one layer down, in
    #:   `physical.verify()`, where both sides are full-precision.
    DERIVED_RTOL = 1e-3

    def check_derived(self) -> list:
        """Recompute `derived` from `geometry.d`, `energy.T` and `medium.eta`
        and report every disagreement.

        This is `physical.verify()`'s invariant one layer up: a derived number
        that does not follow from the inputs beside it was either typed or
        stale, and in the measured case it was a placeholder zero with an
        excuse in its provenance field.

        Returns a list of strings, empty when it agrees. A missing input is
        REPORTED, never silently skipped -- an unrunnable check that returns
        `[]` is indistinguishable from a check that ran and passed, which is
        this repository's most-repeated defect.
        """
        from . import materials as MAT
        from .units import Q

        have = {c: (self.params.get(c) or {}) for c in CATEGORIES}
        need = {"d": ("geometry", "d"), "T": ("energy", "T"),
                "eta": ("medium", "eta")}
        missing = [f"{c}.{k}" for k, (c, k) in need.items() if k not in have[c]]
        if missing:
            return [f"`derived` cannot be recomputed: {missing} absent. "
                    f"State them, or the derived block is unchecked."]
        try:
            q = {k: Q(float(have[c][name].value), have[c][name].unit)
                 for k, (c, name) in need.items()}
            b = MAT.sphere_bulk(q["d"], q["T"], q["eta"])
        except Exception as exc:                       # a bad unit string
            return [f"`derived` cannot be recomputed ({type(exc).__name__}: "
                    f"{exc}). A check that cannot run is not a check that passed."]

        # ★ `tau_gov` is NOT tau_B, and assuming it was refused a correct
        #   manifest too. CLAUDE.md rule 10 asks for "gamma, D_t, tau_B, and
        #   WHICH TIMESCALE GOVERNS" -- four things, and the fourth is a
        #   case-dependent choice. `trap-2d-5um` declares
        #   `tau_gov = 4.010e-3 s` with the provenance "tau_k = gamma/k governs,
        #   NOT tau_B = 242 s", and it is right: the trap relaxation is five
        #   orders faster than diffusion across the particle. So `tau_gov` is
        #   checked for being a positive time; `tau_B` is checked against the
        #   recomputation when the manifest carries it.
        want = {"gamma": b["gamma"], "D_t": b["D_t"], "tau_B": b["tau_B"]}
        out = []
        gov = have["derived"].get("tau_gov")
        if gov is not None:
            if str(gov.unit).strip() in REDUCED_UNITS:
                if float(gov.value) <= 0:
                    out.append(
                        f"derived.tau_gov = {gov.value} {gov.unit!r}: a timescale "
                        f"cannot be zero or negative -- provenance says "
                        f"{gov.provenance[:40]!r}")
            else:
                try:
                    secs = float(Q(float(gov.value), gov.unit).to("s").magnitude)
                except Exception as exc:
                    out.append(f"derived.tau_gov = {gov.value} {gov.unit!r} is not "
                               f"a time ({type(exc).__name__})")
                else:
                    if secs <= 0:
                        out.append(
                            f"derived.tau_gov = {secs:.6g} s: a timescale cannot be "
                            f"zero or negative -- provenance says "
                            f"{gov.provenance[:40]!r}")
        for name, expect in want.items():
            p = have["derived"].get(name)
            if p is None:
                continue                # `blockers` already reports it missing
            # ★ Two conventions, both checked, neither skipped.
            #   A manifest written in REDUCED units says `gamma = 1` because the
            #   reference scales are chosen to make kT = D = gamma = 1 -- there
            #   is no SI recomputation to compare against, but the value is
            #   fixed at 1 BY CONSTRUCTION, so it is still checkable.
            #   `bead-water-3d`'s manifest is of this kind and is correct.
            #   The placeholder that started all this was `0.0` with unit "1",
            #   which this branch refuses.
            if str(p.unit).strip() in REDUCED_UNITS:
                if float(p.value) != 1.0:
                    out.append(
                        f"derived.{name} = {p.value} {p.unit!r}: in the reduced "
                        f"convention the reference scales make kT = D = gamma = 1, "
                        f"so this is 1 by construction. {p.value} is not a value, "
                        f"it is a placeholder -- provenance says "
                        f"{p.provenance[:40]!r}")
                continue
            try:
                got = Q(float(p.value), p.unit).to(expect.units)
            except Exception as exc:
                out.append(f"derived.{name} = {p.value} {p.unit!r} is not "
                           f"convertible to {expect.units:~P} "
                           f"({type(exc).__name__})")
                continue
            e = float(expect.magnitude)
            g = float(got.magnitude)
            if abs(g - e) > self.DERIVED_RTOL * abs(e):
                out.append(
                    f"derived.{name} = {g:.6g} {expect.units:~P} but "
                    f"d, T and eta give {e:.6g} "
                    f"({'placeholder zero' if g == 0 else f'{100*(g-e)/e:+.3g} %'})"
                    f" -- provenance says {p.provenance[:44]!r}")
        return out

    def blockers(self) -> list:
        """Every reason this manifest is not ready to run."""
        out = []
        for c in CATEGORIES:
            if c in self.absent:
                continue
            have = set(self.params.get(c) or {})
            if not have:
                out.append(f"category {c!r} is empty and not declared N/A -- "
                           f"state the numbers or say why there are none")
                continue
            missing = [k for k in REQUIRED_KEYS[c] if k not in have]
            if missing:
                out.append(f"category {c!r} is missing {missing}")
        #  ★ presence is not validity. Added 2026-09-14 after three placeholder
        #    zeros cleared this function -- see the docstring's "four for four".
        if "derived" not in self.absent:
            out += self.check_derived()
        for name, who in self.unknowns.items():
            out.append(f"{name} is UNKNOWN -- needs {who}. BLOCKED is correct here; "
                       f"inventing it is not (rule 3)")
        if not self.approved_by:
            out.append("approved_by is unset. Rule 10: lay out every number and "
                       "ASK before running. No code sets this field.")
        return out

    @property
    def ready(self) -> bool:
        return not self.blockers()

    # -- output -------------------------------------------------------------
    def render(self) -> str:
        W = 78
        L = ["=" * W, f"PARAMETER MANIFEST — {self.case}",
             "rule 10: every number, before the run", "=" * W]
        for c in CATEGORIES:
            if c in self.absent:
                L.append(f"\n{c.upper()}   — NOT APPLICABLE")
                L.append(f"  reason: {self.absent[c]}")
                continue
            ps = self.params.get(c) or {}
            L.append(f"\n{c.upper()}")
            if not ps:
                L.append("  ** EMPTY and not declared N/A **")
            for p in ps.values():
                L.append(f"  {p.line()}")
                if p.note:
                    L.append(f"      {p.note[:66]}")
        if self.unknowns:
            L.append("\nNOT KNOWN (blocks — rule 3: do not invent)")
            for n, who in self.unknowns.items():
                L.append(f"  {n:<14} needs {who}")
        tc = self.tier_counts()
        if tc:
            L.append("\nTIER DISTRIBUTION")
            mean = {0: "given / handbook", 1: "literature, verified",
                    2: "literature, unverified", 3: "arbitrary assumption"}
            for t, n in tc.items():
                L.append(f"  tier {t}  {n:>3}   {mean.get(t, '')}")
            if tc.get(3):
                L.append(f"  ⚠ {tc[3]} value(s) at tier 3 are arbitrary choices, "
                         f"not measurements")
        L.append("")
        L.append("=" * W)
        b = self.blockers()
        if b:
            L.append("VERDICT: NOT READY")
            for x in b:
                L.append(f"  ✗ {x}")
        else:
            L.append(f"VERDICT: READY — approved by {self.approved_by}")
        L.append("=" * W)
        return "\n".join(L)

    def to_json(self) -> dict:
        return {
            "schema": SCHEMA, "case": self.case,
            "params": {c: {n: {"value": p.value, "unit": p.unit,
                               "provenance": p.provenance, "tier": p.tier,
                               "note": p.note}
                           for n, p in (self.params.get(c) or {}).items()}
                       for c in CATEGORIES if self.params.get(c)},
            "not_applicable": dict(self.absent),
            "unknown": dict(self.unknowns),
            "approved_by": self.approved_by,
            "tier_counts": {str(k): v for k, v in self.tier_counts().items()},
            "blockers": self.blockers(),
        }

    def write(self, outdir) -> Path:
        p = Path(outdir) / "params.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.to_json(), indent=1, ensure_ascii=False) + "\n")
        return p


def load(path) -> Manifest:
    """Read back a written manifest. Used to check a run against what was approved."""
    p = Path(path)
    if p.is_dir():
        p = p / "params.json"
    d = json.loads(p.read_text())
    m = Manifest(case=str(d.get("case", p.parent.name)),
                 approved_by=d.get("approved_by"))
    for c, items in (d.get("params") or {}).items():
        for n, v in items.items():
            m.add(c, n, v["value"], v["unit"], v["provenance"],
                  int(v.get("tier", 3)), v.get("note", ""))
    for c, r in (d.get("not_applicable") or {}).items():
        m.not_applicable(c, r)
    for n, who in (d.get("unknown") or {}).items():
        m.unknown(n, who)
    return m


def diff_against_spec(manifest: Manifest, spec_numerics: dict,
                      rtol: float = 1e-9) -> list:
    """Do the approved numbers match what the run is about to execute?

    Approving a manifest and then running something else is the failure this
    guards. Compares `numerics` names that appear in both.
    """
    out = []
    have = manifest.params.get("numerics") or {}
    for name, p in have.items():
        for key in (name, f"{name}_star", f"n_{name}"):
            if key in spec_numerics:
                a, b = float(p.value), float(spec_numerics[key])
                if a == 0.0 and b == 0.0:
                    break
                if abs(a - b) > rtol * max(abs(a), abs(b), 1e-300):
                    out.append(f"{name}: manifest says {a:.6g}, the spec says "
                               f"{b:.6g} ({key}) -- the approved numbers are not "
                               f"the numbers about to run")
                break
    return out


__all__ = ["SCHEMA", "CATEGORIES", "REQUIRED_KEYS", "TIERS", "Param", "Manifest",
           "load", "diff_against_spec"]
