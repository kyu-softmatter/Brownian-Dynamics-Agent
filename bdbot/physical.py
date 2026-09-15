"""The L2 `PhysicalSystem` -- the dimensional physical system (SI). The starting
point of rule 1 ("dimensions come first").

The first two cases had each grown their own `load_system()` -- **it appeared
twice**, so it is shared. The schema was split out of how the two real
`system.yaml` files were actually used:

    in both      label . description . dimensions . particle . medium .
                interactions · external · targets · numerics
    in one       geometry . derived_scales . dimensionless .
                required_convergence_checks · not_verified

The **Provenanced leaves** (`value` + `unit` + `source` + `tier`) were identical
in both: `particle.{diameter,density,count}`, `medium.{temperature,viscosity}`,
plus the per-case interactions.

`derived_scales` is **a deliberate exception** -- derived values like gamma, D_t
and tau_B are verified not by a source but **by recomputation** (`verify()` does
that).

Invariants enforced here:
  1. `derived_from` -- the `observation.yaml` this physical system came from.
     Rejected if absent. (The two original files recorded this only in a
     **comment**, which a machine cannot read.)
  2. If L0 is BLOCKED, L2 cannot exist -- settling a physical system while an
     unresolved physical gap remains means a value was invented somewhere
     (rule 3).
  3. A value composed only of tier >= 2 (unverified) needs human approval.
"""
from __future__ import annotations

import pathlib

import math
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import intake as _intake
from . import materials as _mat
from .provenance import Provenanced, load_node
from .units import Q

SCHEMA = "bdbot.system/0.1"

REQUIRED_TOP = ("label", "dimensions", "particle", "medium", "targets", "numerics")
OPTIONAL_TOP = ("description", "geometry", "interactions", "external", "derived_scales",
                "dimensionless", "required_convergence_checks", "not_verified",
                "derived_from")
# Sections allowed to carry a value with no source -- they are derived, and are
# verified by recomputation
DERIVED_SECTIONS = ("derived_scales", "dimensionless", "friction")
# Provenanced paths that were identical in both cases.
# * A generalization forced by the third case (abp-rod, an ellipsoid): some shapes
#   do not obey the sphere formulas.
#   gamma_bar (Perrin, 2D harmonic mean) = 7.21e-9 vs 3*pi*eta*d_eq = 6.37e-9 --
#   a 13% difference. Running the sphere recomputation check unchanged would
#   **flag a correct spec as an error.**
SPHERICAL_SHAPES = ("sphere", None)

CORE_PROVENANCED = {
    "d": ("particle", "diameter"),
    "rho_p": ("particle", "density"),
    "N": ("particle", "count"),
    "T": ("medium", "temperature"),
    "eta": ("medium", "viscosity"),
}
TIER_MEANING = {0: "given / handbook", 1: "literature + verified, or a confirmed convention",
                2: "literature, unverified", 3: "arbitrary assumption"}

# ════════════════════════════════════════════════════════════════════════════
# `structure.dim` -- teeth for D9
# ════════════════════════════════════════════════════════════════════════════
#  `knowledge/wiki/concepts/dimensionality-has-no-default.md` closes
#  `D9 · 차원 (2D / 3D)`, which had stood `OPEN` with a 3D default that **7 of
#  the 8 cases contradict**. The decision: there is no default; every case
#  declares WHY its dimension is what it is.
#
#  This block is what stops that decision from being the fourth entry in rule
#  10's list of practices that were written and never enforced (`bd-intake`
#  §2.1's empty-goal blocker, `A4`'s grep, `health.gate()` -- three for three).
#  Measured 2026-09-10, before this existed: `dimensions:` was a bare scalar in
#  **8 of 8** `system.yaml`, with no `source` and no `tier` -- less provenance
#  than `particle.density`, which only feeds the `tau_p` sanity check.
#
#  ⚠ Deliberately NOT hashed into `run_id` (`runid.DOC_KEYS`). The physics value
#  stays where it was -- top-level `dimensions` -- so 2 -> 3 still re-ids every
#  run, while writing down *why* it is 2 does not. Both directions of
#  `runid.py`'s warning are preserved, and `_check_dim` cross-checks the two so
#  they cannot drift.
STRUCTURE_SECTION = "structure"

#: `given` beats the other three. If the source -- sketch, paper, or a human
#: asked directly -- states the dimension, it is an input and is not re-derived
#: (rule 3: never invent; rule 5: transcribe before you interpret).
DIM_BASES = ("given", "required", "inherited", "sufficient")

DIM_BASIS_MEANING = {
    "given":      "the input states it -- not a modelling choice",
    "required":   "the question does not exist in another dimension",
    "inherited":  "matched to an existing case so the two stay comparable",
    "sufficient": "the DOFs separate and the observable is a separated component",
}

#: Each basis owes a different piece of evidence, and two of them **expire**:
#: `inherited` when the comparison target changes, `sufficient` when an
#: observable is added. The field is what makes the expiry detectable.
DIM_BASIS_REQUIRES = {
    "given":      ("source",),
    "required":   (),
    "inherited":  ("compared_with",),
    "sufficient": ("checked_observables",),
}

#: `basis: given` asserts "the input states it". But WHICH input? A `source`
#: field satisfies the check below by existing, and it cannot tell a file from a
#: sentence someone is recalling -- the same "presence is not validity" shape
#: that `params.blockers()` had for `derived` (see bdbot/params.py, "four for
#: four"). Measured 2026-09-15: of the three `given` cases, two cite an artefact
#: and one rested entirely on a conversation, and the gate passed all three
#: identically. Its own prose said "The sketch itself is silent."
#:
#: So the author DECLARES which kind, and each kind owes a different thing. Not
#: inferred from the text: a keyword scan for "sketch"/"observation.yaml" in the
#: source string mislabelled exactly the one conversational case, because its
#: source says the sketch is silent. Parse, do not grep.
DIM_GIVEN_KINDS = ("artefact", "human")
DIM_GIVEN_MEANING = {
    "artefact": "a file in the record states it -- and the file must exist",
    "human":    "a person stated it, in conversation -- so it needs confirming",
}
DIM_GIVEN_REQUIRES = {
    "artefact": ("source_file",),
    "human":    ("confirmed_by",),
}

#: Always required, whatever the basis. `what_would_change` is the load-bearing
#: one: empty means "I do not know whether this choice is safe", which is a
#: different state from "it is safe", and the whole point is that the two stop
#: looking alike.
DIM_ALWAYS = ("value", "basis", "alternatives", "what_would_change")

# ── `structure.size` ───────────────────────────────────────────────────────
#  Read out of the 8 cases (2026-09-10). `particle.count` was carrying three
#  different quantities under one name:
#
#    object       N is the size of the THING -- a 25-bead chain. Raising it makes
#                 a different object, so an N-sweep is physics, not convergence.
#                 `chain-bend-2d-dlvo` reports its n-dependence as evidence for
#                 the conclusion; calling that a finite-size artefact would file
#                 a finding as a failure (rule 7').
#    system       N is how much of an infinite system is simulated. phi links N
#                 and L but fixes neither -- a SECOND constraint fixes one and
#                 phi gives the other, and that second constraint is what has to
#                 be recorded. `soft-r3` went 100 -> 400 because of `r_c < L/2`,
#                 not because phi moved.
#    replication  the physical system is ONE particle and N is a numerical
#                 multiplier. `trap-2d-5um` records exactly this in prose
#                 ("simulation choice ... the sketch has 1") while filing 1000
#                 under `particle.count`, next to diameter and density, where it
#                 reads as a property of the physical system. It is not one; it
#                 belongs with `T_obs` and the seed count as sampling.
#
#  The three carry opposite obligations, which is why one field cannot serve all
#  three: `converge --only N_double` passing is reassurance for `replication`, a
#  category error for `object`, and possibly a false negative for `system` (the
#  two sizes may simply have been too close -- findings/
#  order-parameter-magnitude-cannot-identify-a-phase).
SIZE_ROLES = ("object", "system", "replication")

SIZE_ROLE_MEANING = {
    "object":      "N is the size of the thing being studied",
    "system":      "N sets how much of an infinite system is simulated",
    "replication": "the physical system is smaller; N is a sampling multiplier",
}

#: What actually fixed the size. Valid values depend on the role, because N is
#: never chosen directly -- something else is chosen and N follows.
SIZE_FIXED_BY = {
    "object":      ("given", "inherited", "design_power", "swept"),
    "system":      ("minimum_image", "measurement_radius", "commensurability",
                    "correlation_length", "cost"),
    # `confirmed_by_run` is weaker than `target_precision` and is named
    # separately because it is what `trap-2d-5um` actually did: a round number,
    # then `confirmed_by: run` once the four observables matched the analytic
    # solution. Recording it as a derived target would be tidying after the fact.
    "replication": ("target_precision", "confirmed_by_run", "inherited", "cost"),
}

SIZE_ALWAYS = ("role", "n", "fixed_by", "what_would_change", "stated_in_source")

#: REQUIRED, and `null` is the answer when the source is silent -- because
#: omitting the key is how the override obligation below quietly disappears.
#: `_require` treats an explicit `null` as present-but-blank, so it is listed in
#: `SIZE_NULLABLE` and checked for presence only.
#:
#: What the INPUT said, if anything.
#: Unlike `dim`, a stated N does NOT automatically win: `soft-r3`'s sketch says
#: N=100 and that value cannot satisfy `r_c < L/2`, so physics legitimately
#: overrode it. What may not happen is overriding it QUIETLY -- so when
#: `stated_in_source` differs from `n`, `approved_by` becomes required.
SIZE_STATED = "stated_in_source"

#: Fields where an explicit `null` is a real answer, not a blank. "The source is
#: silent on N" has to be sayable, and has to be distinguishable from "nobody
#: filled this in" -- which is exactly why the key is mandatory.
SIZE_NULLABLE = ("stated_in_source",)
SIZE_APPROVAL = "approved_by"

#: `system` owes three more.
#:
#: ⚠ `pins` was called `follows` and held the PINNED quantity -- i.e. the exact
#:   opposite of what the name reads, and nothing could tell the two conventions
#:   apart, so a file written to the name rather than the docstring would be
#:   silently backwards. Renamed after adversarial review. `pins: L` = the
#:   second constraint fixes L and N follows from phi.
#:   2 of the 3 cases pin L; `trap-drag` pins N (an integer pair n_x*n_y).
#:
#: `finite_size` may say "not done", but it may not be blank: measured
#: 2026-09-10, only 1 of the 3 system-role cases has an N-sweep, while
#: `trap-drag` reports collective observables from 81 runs at a single size.
SIZE_SYSTEM_REQUIRES = ("phi", "pins", "finite_size", "finite_size_status")

#: Where the run's actual packing fraction lives, by dimensionality. `size.phi`
#: is cross-checked against it the way `n` is against `particle.count` -- before
#: this, `phi: banana` and a 2.6x-wrong phi both validated clean.
GEOMETRY_PHI_KEYS = ("area_fraction", "volume_fraction_final", "volume_fraction")

#: A status word, not prose. The prose still has to be there (`finite_size`), but
#: what the gate branches on is this.
SIZE_FS_STATUS = ("done", "partial", "not_done")
SIZE_PINS = ("n", "L")

#: `replication` owes the split it is asserting, and `n_physical * n_replicas`
#: is cross-checked against `particle.count` -- otherwise the claim "the system
#: is really one particle" is unfalsifiable.
SIZE_REPLICATION_REQUIRES = ("n_physical", "n_replicas")


@dataclass
class PhysicalSystem:
    path: Path
    raw: dict
    core: dict = field(default_factory=dict)      # name -> Provenanced
    issues: list = field(default_factory=list)

    @property
    def label(self) -> str:
        return self.raw.get("label", "?")

    @property
    def dim(self) -> int:
        """`dimensions` as an int, or 0 when it is not one.

        ⚠ CLASS A, instance #5. This was `int(self.raw.get("dimensions", 0))`,
        which raises on `dimensions: 2D` -- and `render_check` reads it, so
        `bdbot system check` tracebacked on exactly the value the new gate was
        built to REPORT. `_int_or_none` was added one round earlier for the same
        expression inside `_phi_closure` and was not applied here, which is the
        shape this class keeps taking: the helper lands on the site that was
        found, not on the class.
        """
        return _int_or_none(self.raw.get("dimensions")) or 0

    @property
    def shape(self):
        """`particle.shape` -- treated as a sphere when absent (as in the first two cases)."""
        return _mapping(self.raw.get("particle")).get("shape")

    @property
    def is_spherical(self) -> bool:
        return self.shape in SPHERICAL_SHAPES

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.level == "error"]

    def node(self, *path, required: bool = True):
        """Access a per-case Provenanced node. `sys_.node("external","stiffness")`."""
        cur = self.raw
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                if required:
                    raise KeyError(f"{'.'.join(map(str, path))} not found ({self.path})")
                return None
            cur = cur[k]
        return load_node(cur)

    def tiers(self) -> dict:
        out: dict[int, list[str]] = {}
        for name, p in self.core.items():
            out.setdefault(p.tier, []).append(name)
        for name, p in self._extra_provenanced().items():
            out.setdefault(p.tier, []).append(name)
        return dict(sorted(out.items()))

    def _extra_provenanced(self) -> dict:
        """Provenanced leaves outside core (per-case interactions, geometry, ...)."""
        out: dict[str, Provenanced] = {}
        for path, node in _walk_provenanced(self.raw):
            top = path.split(".")[0].split("[")[0]
            if top in DERIVED_SECTIONS or top == STRUCTURE_SECTION:
                continue
            if path in {".".join(p) for p in CORE_PROVENANCED.values()}:
                continue
            try:
                out[path] = load_node(node)
            except Exception:
                pass
        return out

    def bulk(self) -> dict | None:
        """Basic properties of a sphere in a Newtonian fluid. Both of the first two
        cases used exactly this bundle.

        * Raises nothing. If the checker crashes on a spec with broken units, the
          user gets a traceback instead of "what is wrong" -- which is exactly what
          an adversarial test produced (feeding `furlong^2` killed it with a
          DimensionalityError).
        """
        if not all(k in self.core for k in ("d", "T", "eta", "rho_p")):
            return None
        try:
            return _mat.sphere_bulk(self.core["d"].value, self.core["T"].value,
                                    self.core["eta"].value, self.core["rho_p"].value)
        except Exception:
            return None


def _walk_provenanced(d, pre=""):
    """Enumerate leaves carrying `value` plus (`unit`|`source`|`tier`) as (path, node)."""
    out = []
    if isinstance(d, dict):
        if "value" in d and any(k in d for k in ("unit", "source", "tier")):
            out.append((pre, d))
        else:
            for k, v in d.items():
                out += _walk_provenanced(v, f"{pre}.{k}" if pre else str(k))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out += _walk_provenanced(v, f"{pre}[{i}]")
    return out


def load(path) -> PhysicalSystem:
    p = Path(path)
    if p.is_dir():
        p = p / "system.yaml"
    if not p.exists():
        s = PhysicalSystem(p, {})
        s.issues.append(_intake.Issue("error", str(p), "system.yaml is missing."))
        return s
    raw = yaml.safe_load(p.read_text()) or {}
    s = PhysicalSystem(p, raw)
    for name, keys in CORE_PROVENANCED.items():
        cur = raw
        ok = True
        for k in keys:
            if not isinstance(cur, dict) or k not in cur:
                ok = False
                break
            cur = cur[k]
        if ok:
            try:
                s.core[name] = load_node(cur)
            except Exception as e:
                s.issues.append(_intake.Issue("error", ".".join(keys), f"parse failed: {e}"))
    s.issues += validate(s)
    return s


def validate(s: PhysicalSystem) -> list:
    I = _intake.Issue
    raw = s.raw
    out: list = []

    # 1. required sections
    for k in REQUIRED_TOP:
        if k not in raw:
            out.append(I("error", k, "required section missing"))
    for name, keys in CORE_PROVENANCED.items():
        if name not in s.core:
            out.append(I("error", ".".join(keys), "required Provenanced node missing (present in both reference cases)"))

    # 2. * invariant: derived_from
    df = raw.get("derived_from")
    if not df:
        out.append(I("error", "derived_from",
                     "there is no record of which observation.yaml this physical "
                     "system came from. It has to be a field, not a comment, for a "
                     "machine to verify it."))
    else:
        ref = (s.path.parent / Path(str(df)).name)
        if not ref.exists():
            out.append(I("error", "derived_from", f"referenced file does not exist: {df}"))
        else:
            # 3. * if L0 is BLOCKED, L2 cannot exist
            obs = _intake.load(ref)
            if obs.errors:
                out.append(I("error", "derived_from",
                             f"the source observation.yaml has {len(obs.errors)} schema error(s)."))
            if obs.open_missing:
                names = ", ".join(m.get("symbol", "?") for m in obs.open_missing)
                out.append(I("error", "derived_from",
                             f"L0 has unresolved physical gaps but the physical system is "
                             f"settled: {names}. A value may have been invented somewhere "
                             f"(rule 3)."))

    # 4. completeness of the Provenanced leaves, plus unit parsing
    for path, node in _walk_provenanced(raw):
        top = path.split(".")[0].split("[")[0]
        # `structure` holds Choices, not Provenanced numbers -- different schema,
        # validated by `check_dim`.
        if top in DERIVED_SECTIONS or top == STRUCTURE_SECTION:
            continue
        for need in ("source", "tier"):
            if need not in node:
                out.append(I("error", path, f"'{need}' missing (rule 3: every number carries a source)"))
        if "unit" in node and node["unit"] is not None:
            try:
                Q(1.0, str(node["unit"]))
            except Exception as e:
                out.append(I("error", path, f"cannot parse the unit: {node['unit']} ({e})"))
        t = node.get("tier")
        if t is not None and t not in TIER_MEANING:
            out.append(I("error", path, f"tier must be 0-3 (got {t})"))

    # 5. recomputation check on the derived values (when present)
    out += verify(s)

    # 5b. the structural choices -- dimensionality (D9) and system size
    out += check_dim(s)
    out += check_size(s)

    # 6. the tier approval gate
    low = [n for n, p in {**s.core, **s._extra_provenanced()}.items() if p.tier >= 2]
    if low:
        out.append(I("warn", "tier", f"{len(low)} value(s) at tier >= 2 (unverified): {', '.join(low[:6])}"
                                     f"{' ...' if len(low) > 6 else ''} -- needs human approval"))
    return out


def _int_or_none(v):
    """`int(v)` or `None` -- never raises.

    ⚠ THE FOURTH instance of one regression class. `_mapping` ended it for
    *containers* (`particle`, `geometry`, `structure`), and the comment there
    says "one helper, so the next container cannot be missed" -- which was true
    and beside the point, because `dimensions` is a top-level SCALAR. An
    unguarded `int(s.raw.get("dimensions", 0) or 0)` in `_phi_closure` then
    raised out of `load()` on `dimensions: 2D` (a plausible hand-authoring typo;
    `intake/` has no generator and rule 2 protects only `specs/` and `runs/`),
    taking the whole `bdbot status` table down -- `cli.py` loops every case with
    no guard and prints only after the loop. `check_dim` already REPORTS a
    non-numeric `dimensions` by comparing without converting, so nothing here
    needs to convert unsafely.
    """
    try:
        return int(v)
    except (TypeError, ValueError):
        return None


def _mapping(v) -> dict:
    """`v` if it is a mapping, else `{}`.

    Replaces `x = raw.get(k) or {}`, which rescues only FALSY values -- so a
    hand-authored `geometry: periodic_square_2d` or `particle: sphere` survived
    it and was then dotted, and `validate()` raised instead of returning its
    issues. The round-1 fix guarded `particle.count` and left both enclosing
    lookups on the old idiom; round 3 found them. One helper, so the next
    container cannot be missed.
    """
    return v if isinstance(v, dict) else {}


def _seq(v):
    """`list(v)` for a genuine sequence, else `None`.

    ⚠ `isinstance(v, (list, tuple))` was the test, and pint hands back a NUMPY
    ARRAY for a list-valued quantity (`core["N"].value.magnitude` is
    `array([5, 9, 15, 25])`). So the sequence branch below was DEAD for exactly
    the 3 cases that carry list counts, and `_same_number` returned an ARRAY of
    booleans instead of a bool -- `not <array>` then raises. Strings are
    excluded on purpose: "400" is a scalar here, not a 3-element sequence.
    """
    if isinstance(v, (str, bytes)) or isinstance(v, dict):
        return None
    try:
        return list(v)
    except TypeError:
        return None


def _same_number(a, b, rtol: float = 1e-9) -> bool:
    """Scalar-or-sequence equality, ALWAYS a bool. `particle.count` is a list on
    the sweep cases (`chain-bend-2d-dlvo`: [5, 9, 15, 25]), and pint turns it
    into an ndarray, so the cross-check has to survive both.
    """
    sa, sb = _seq(a), _seq(b)
    if sa is not None or sb is not None:
        sa = sa if sa is not None else [a]
        sb = sb if sb is not None else [b]
        return len(sa) == len(sb) and all(bool(_same_number(x, y, rtol))
                                          for x, y in zip(sa, sb))
    try:
        return bool(abs(float(a) - float(b)) <= rtol * max(1.0, abs(float(b))))
    except (TypeError, ValueError):
        try:
            return bool(a == b)
        except Exception:
            return False


def _blank(v) -> bool:
    """Present but empty. Blank must mean *forgotten*, never *nothing to say*."""
    return (v is None
            or (isinstance(v, str) and not v.strip())
            or (isinstance(v, (list, tuple, dict)) and not v))


def _require(e: dict, where: str, fields, why: str = "", nullable=()) -> list:
    I = _intake.Issue
    out = []
    for f in fields:
        if f in nullable:
            if f not in e:
                out.append(I("error", f"{where}.{f}",
                             "missing -- write `null` if the source is silent. Omitting "
                             "the key is how an obligation disappears quietly."))
            continue
        if f not in e:
            out.append(I("error", f"{where}.{f}", "missing" + (f" -- {why}" if why else "")))
        elif _blank(e[f]):
            out.append(I("error", f"{where}.{f}",
                         "present but empty. Blank must mean *forgotten*, never "
                         "*nothing to say*"
                         + (" -- an empty `what_would_change` means you do not know "
                            "whether this choice is safe, which is not the same as "
                            "it being safe."
                            if f == "what_would_change" else ".")))
    return out


def _source_stated_n(s: PhysicalSystem):
    """What the L0 transcription says about `N` -- `(value, where)` or `(None, where)`.

    `structure.size.stated_in_source` was self-reported and checked against
    nothing, so the override gate could be bypassed by writing the run's own `n`
    into it. The transcription is right there: `validate()` already loads
    `derived_from`, and `stated_quantities` carries `symbol: N` on exactly the
    two cases whose sketches state it.
    """
    ref = s.path.parent / Path(str(s.raw.get("derived_from") or "observation.yaml")).name
    if not ref.exists():
        return None, None
    try:
        obs = _intake.load(ref)
    except Exception:
        return None, None
    for q in (obs.raw.get("stated_quantities") or []):
        if isinstance(q, dict) and str(q.get("symbol")).strip() == "N":
            return q.get("value"), ref.name
    return None, ref.name


def check_size(s: PhysicalSystem) -> list:
    """`structure.size` -- what fixed N, and which kind of N it is.

    `particle.count` was carrying three different quantities under one name, and
    the three owe different evidence. See `SIZE_ROLES` for the reading of the 8
    cases that produced this.

    The load-bearing asymmetry: **N is never chosen directly.** For `object`
    something outside the box fixes it (a paper, a design-power inequality, a
    comparison, or it is the swept variable); for `system` phi links N and L
    while a *second* constraint pins one of them; for `replication` the target
    precision buys it. `fixed_by` records which, and its valid values depend on
    the role.
    """
    I = _intake.Issue
    raw = s.raw
    if not raw:
        return []

    st = _mapping(raw.get(STRUCTURE_SECTION))
    e = st.get("size")
    #  ⚠ `or {}` only rescues FALSY values, so the shorthand `count: 400` (or a
    #     list, or a string) survived it and was then dotted -- turning a
    #     reported schema error into a traceback, and taking the whole
    #     `bdbot status` table down with it (cli.py:170 loops with no guard).
    #     A validator that crashes on malformed input is the failure its own
    #     `bulk()` docstring already names. Found by adversarial review.
    #  ⚠ ROUND 3: the guard below was added for `count` while the ENCLOSING
    #     lookups kept the condemned idiom, so a non-mapping `particle:` or
    #     `geometry:` still raised -- a strict regression for the 3 `role: system`
    #     cases, with the same `bdbot status` blast radius. Normalize containers.
    cnt = _mapping(raw.get("particle")).get("count")
    cnt = _mapping(cnt)
    if not isinstance(e, dict):
        return [I("error", f"{STRUCTURE_SECTION}.size",
                  f"missing. `particle.count: {cnt.get('value')}` sits beside "
                  f"diameter and density, so it reads as a property of the physical "
                  f"system -- and for a replicated single-particle case it is not "
                  f"one. Declare a role: {', '.join(SIZE_ROLES)}. "
                  f"See knowledge/wiki/concepts/dimensionality-has-no-default.md")]

    out: list = []
    where = f"{STRUCTURE_SECTION}.size"
    role = e.get("role")
    if not isinstance(role, str):        # a list/dict role is unhashable -> `in` raises
        role = None
    if role not in SIZE_ROLES:
        out.append(I("error", f"{where}.role",
                     f"must be one of {SIZE_ROLES} (got {role!r}). "
                     + " . ".join(f"{k}={v}" for k, v in SIZE_ROLE_MEANING.items())))

    need = SIZE_ALWAYS
    if role == "system":
        need += SIZE_SYSTEM_REQUIRES
    elif role == "replication":
        need += SIZE_REPLICATION_REQUIRES
    out += _require(e, where, need,
                    f"required when role is {role!r}" if role in SIZE_ROLES else "",
                    nullable=SIZE_NULLABLE)

    #  the value itself has to be checkable, not prose
    stated = e.get(SIZE_STATED)
    if stated is not None and not isinstance(stated, (int, float, list, tuple)):
        out.append(I("error", f"{where}.{SIZE_STATED}",
                     f"must be a number, a list, or `null` (got {stated!r}). Prose "
                     f"here cannot be compared with `n`, so the override gate "
                     f"below would never fire."))
        stated = None

    #  ...and it has to match the transcription, or the gate is bypassed by
    #  simply writing the run's own n into the field
    src_n, src_file = _source_stated_n(s)
    if src_file is not None and not _same_number(stated, src_n):
        out.append(I("error", f"{where}.{SIZE_STATED}",
                     f"= {stated!r} but {src_file} `stated_quantities` says "
                     f"N = {src_n!r}. This field records what the INPUT said; it "
                     f"is not a free slot, and writing `n` into it is how the "
                     f"approval gate below gets bypassed."))

    #  a stated N may be overridden -- but not silently (rule 10's `approved_by`)
    if SIZE_STATED in e and not _blank(e[SIZE_STATED]) \
            and not _same_number(e[SIZE_STATED], e.get("n")):
        if _blank(e.get(SIZE_APPROVAL)):
            out.append(I("error", f"{where}.{SIZE_APPROVAL}",
                         f"the source states N = {e[SIZE_STATED]} and this run uses "
                         f"{e.get('n')}. Overriding a stated value is allowed -- "
                         f"`soft-r3` had to -- but it needs a human: record who "
                         f"approved it and when."))

    ok = e.get("fixed_by")
    if role in SIZE_FIXED_BY and ok is not None and ok not in SIZE_FIXED_BY[role]:
        out.append(I("error", f"{where}.fixed_by",
                     f"for role {role!r} must be one of {SIZE_FIXED_BY[role]} "
                     f"(got {ok!r}). N is never chosen directly -- name what was."))

    # the choice and the physics field must be one number
    if "n" in e and isinstance(cnt, dict) and "value" in cnt \
            and not _same_number(e["n"], cnt["value"]):
        out.append(I("error", f"{where}.n",
                     f"disagrees with `particle.count.value` = {cnt['value']} "
                     f"(structure says {e['n']})."))

    if role == "system":
        #  phi is the equation linking N and L, so a phi that is not a number --
        #  or not the one the run uses -- makes the whole `fixed_by` story
        #  unfalsifiable. Same invariant `n` carries against `particle.count`.
        phi = e.get("phi")
        if phi is not None:
            try:
                phi_f = float(phi)
            except (TypeError, ValueError):
                out.append(I("error", f"{where}.phi",
                             f"must be a number (got {phi!r}). phi is the equation "
                             f"that links N and L; prose here cannot be checked."))
            else:
                geo = _mapping(raw.get("geometry"))
                node = next((geo[k] for k in GEOMETRY_PHI_KEYS
                             if isinstance(geo.get(k), dict) and "value" in geo[k]), None)
                if node is None:
                    out.append(I("warn", f"{where}.phi",
                                 f"no {' / '.join(GEOMETRY_PHI_KEYS)} in `geometry` to "
                                 f"check {phi_f:g} against. Recorded, NOT verified."))
                elif not _same_number(phi_f, node["value"], rtol=1e-6):
                    out.append(I("error", f"{where}.phi",
                                 f"= {phi_f:g} disagrees with the packing fraction the "
                                 f"run uses, `geometry` = {node['value']}."))
                #  ...and agreeing with geometry is not enough: phi = 3.5 (350%
                #  packing) passed both, because nothing closed N, L and phi
                #  against each other. Recompute from the box.
                out += _phi_closure(s, geo, phi_f, where)

        f = e.get("pins")

        if f is not None and (not isinstance(f, str) or f not in SIZE_PINS):
            out.append(I("error", f"{where}.pins",
                         f"must be one of {SIZE_PINS} -- phi links N and L but "
                         f"fixes neither, so name the one the second constraint "
                         f"PINS; the other follows from phi (got {f!r})."))
        #  ⚠ This was a substring search for "not done" in free prose, which fails
        #     in BOTH directions -- "not done" appearing anywhere in a paragraph
        #     that describes a study that WAS done, and any other phrasing of
        #     "no sweep" going unwarned. A status word decides; the prose explains.
        st_ = e.get("finite_size_status")
        if st_ is not None and (not isinstance(st_, str) or st_ not in SIZE_FS_STATUS):
            out.append(I("error", f"{where}.finite_size_status",
                         f"must be one of {SIZE_FS_STATUS} (got {st_!r})."))
        elif st_ == "partial":
            out.append(I("warn", f"{where}.finite_size",
                         "N is a system parameter and the N-sweep is PARTIAL -- a "
                         "ladder exists but does not cover the production point, so a "
                         "collective observable here still carries no N-scaling at the "
                         "point it is reported for "
                         "(findings/order-parameter-magnitude-cannot-identify-a-phase)"))
        elif st_ != "done":
            out.append(I("warn", f"{where}.finite_size",
                         f"finite_size_status={st_!r}: N is a system parameter and NO "
                         f"N-sweep is recorded, so a collective observable here carries "
                         f"no N-scaling "
                         f"(findings/order-parameter-magnitude-cannot-identify-a-phase)"))

    if role == "replication":
        #  "the physical system is really one particle" has to be FALSIFIABLE, so a
        #  non-numeric operand must be reported rather than skipped. The original
        #  `except: pass` swallowed exactly the case the repo writes elsewhere --
        #  a Provenanced `{value: 1, unit: ..., tier: ...}` -- and the check then
        #  passed on a document it had never actually checked.
        if all(f in e and not _blank(e[f]) for f in SIZE_REPLICATION_REQUIRES):
            nums, bad = [], []
            for f in SIZE_REPLICATION_REQUIRES:
                v = e[f]
                if isinstance(v, dict):       # a Provenanced node, not a bare number
                    bad.append(f"{f} is a mapping -- give a bare number here")
                    continue
                try:
                    nums.append(float(v))
                except (TypeError, ValueError):
                    bad.append(f"{f}={v!r} is not a number")
            if bad:
                out.append(I("error", f"{where}.n_replicas",
                             "the falsifiability check could not run: " + "; ".join(bad)))
            elif "value" in cnt and not _same_number(nums[0] * nums[1], cnt["value"]):
                out.append(I("error", f"{where}.n_replicas",
                             f"n_physical x n_replicas = {nums[0] * nums[1]:g} does not "
                             f"equal `particle.count.value` = {cnt['value']}."))
    return out


def _phi_closure(s: PhysicalSystem, geo: dict, phi: float, where: str) -> list:
    """phi, N, d and L are one equation -- close it when the box is recorded.

    `phi = N pi d^2 / (4 L^2)` in 2D, `N pi d^3 / (6 L^3)` in 3D. Agreeing with
    `geometry.area_fraction` is not enough: both can be 3.5 (350 % packing) and
    nothing noticed, because no check ever recomputed phi from the box.
    """
    I = _intake.Issue
    #  ⚠ The box is not always a cube. `trap-drag` is deliberately RECTANGULAR
    #     -- commensurability fixes the aspect ratio at 1.0906 -- so the first
    #     version of this check, which squared `box_length_x`, reported the one
    #     case whose phi is exactly right (0.3209 vs 0.35). Take the PRODUCT of
    #     the per-axis lengths, and fall back to L^dim only for a stated cube.
    #  ⚠ This returned [] in SILENCE when the box was not recorded, and
    #     `network` is exactly that case: 3D, compressed (volume_fraction_initial
    #     -> _final), so phi is the recorded quantity and no box_length_* exists.
    #     The 350 %-packing regression this closure was written to kill therefore
    #     still validated clean on 1 of the 3 system cases. An unrunnable check
    #     has to say so -- "Recorded, NOT verified" -- the same way the
    #     missing-geometry-node path above already does.
    per_axis = [geo[k] for k in ("box_length_x", "box_length_y", "box_length_z")
                if isinstance(geo.get(k), dict) and "value" in geo[k]]
    cube = geo.get("box_length")
    dim0 = _int_or_none(s.raw.get("dimensions"))
    have_box = bool(per_axis) or (isinstance(cube, dict) and "value" in cube)
    if dim0 not in (2, 3):
        return [I("warn", f"{where}.phi",
                  f"= {phi:g} could not be closed: `dimensions` is "
                  f"{s.raw.get('dimensions')!r}, not 2 or 3. Recorded, NOT "
                  f"verified. (`structure.dim` reports the bad value itself.)")]
    if not have_box:
        return [I("warn", f"{where}.phi",
                  f"= {phi:g} could not be closed against N and the box: no "
                  f"box_length / box_length_x.. in `geometry`. Recorded, NOT "
                  f"verified. Agreement with the geometry fraction alone does "
                  f"not exclude an impossible phi.")]
    if "d" not in s.core or "N" not in s.core:
        #  ⚠ was `info`. The reason first recorded here -- "`checks.verdict` does
        #    not surface `info`" -- is WRONG: `checks.verdict` takes `Check`
        #    objects and never receives an `Issue`. The conclusion stands for a
        #    different reason: `run.execute`'s loop filters `i.level == "warn"`,
        #    so an `info` is invisible on exactly the path that needed it. Kept
        #    with the correction rather than quietly reworded (round 6).
        #    So on the `run.execute` path (where `core` is empty by
        #    construction) an impossible phi produced NOTHING while the gate
        #    printed "structure OK -- dim/size validated". A skipped check is a
        #    warning, not an info.
        return [I("warn", f"{where}.phi",
                  f"= {phi:g} NOT closed against N and the box: d / N are not "
                  f"loaded on this path (run.execute passes the spec's system "
                  f"document, whose Provenanced leaves are not parsed). "
                  f"`bdbot.cli system check` runs the closure in full.")]
    if per_axis and len(per_axis) < dim0:
        return [I("warn", f"{where}.phi",
                  f"= {phi:g} could not be closed: {len(per_axis)} of {dim0} box "
                  f"axes recorded. Recorded, NOT verified.")]
    try:
        dim = dim0
        if len(per_axis) >= dim:
            vol = 1.0
            for b in per_axis[:dim]:
                vol *= Q(b["value"], b["unit"]).to("m").magnitude
        else:
            vol = Q(cube["value"], cube["unit"]).to("m").magnitude ** dim
        d = s.core["d"].value.to("m").magnitude
        N = s.core["N"].value.magnitude
        _n = _seq(N)                       # ndarray as well as list/tuple
        N = max(float(x) for x in _n) if _n else float(N)
        phi_box = (N * math.pi * d ** 2 / 4 / vol) if dim == 2 \
            else (N * math.pi * d ** 3 / 6 / vol)
    except Exception as exc:
        #  ⚠ this was a silent `return []`, so a phi that could not be closed for
        #    ANY reason -- a bad unit, a missing key -- read exactly like a phi
        #    that had been closed and agreed. Two such paths survived round 4.
        return [I("warn", f"{where}.phi",
                  f"= {phi:g} could not be closed against N and the box "
                  f"({type(exc).__name__}: {exc}). Recorded, NOT verified.")]
    if not _same_number(phi, phi_box, rtol=2e-2):
        return [I("error", f"{where}.phi",
                  f"= {phi:g} is not what the recorded box gives: N={N:g}, "
                  f"d={d:.4g} m, box volume={vol:.4g} m^{dim} -> "
                  f"phi = {phi_box:.4g}. "
                  f"phi, N and L are one equation; two of them fix the third.")]
    return []


def check_dim(s: PhysicalSystem) -> list:
    """`structure.dim` -- the basis for this case's dimensionality (D9).

    Three kinds of failure, all errors:

      absent        the historical state, 8 of 8 on 2026-09-10. The gate exists
                    to turn that into 8 written answers
      inconsistent  `structure.dim.value` disagrees with the hashed `dimensions`.
                    Two numbers where there should be one -- the same invariant
                    `params.json` carries against the spec's numerics
      empty         a required field present but blank, above all
                    `what_would_change`. "Not stated" and "no consequence" must
                    not render the same

    Two paths reach this. `bdbot.cli system check` / `status` run it in full; and
    `run.execute()` runs it on the spec's own `system` document when that
    document carries a `structure` block (wired 2026-09-11).
    ⚠ This docstring used to say "`run.execute()` never calls this at all" --
    true when written, made false by the wiring, and caught by adversarial
    review 2026-09-13. On the `execute` path `core` is empty, so `_phi_closure`
    and the `stated_in_source` cross-check announce that they did not run rather
    than passing silently.
    The only case script this can stop is `cases/abp_rod_2d.py:518`, the one that
    gates on `s.errors`; every other case reads `system.yaml` directly. With all
    8 bases filled in none is blocked -- measured, `abp_rod_2d.py --smoke` runs.
    """
    I = _intake.Issue
    raw = s.raw
    if not raw:
        return []

    st = _mapping(raw.get(STRUCTURE_SECTION))
    e = st.get("dim")
    if not isinstance(e, dict):
        return [I("error", f"{STRUCTURE_SECTION}.dim",
                  f"missing. `dimensions: {raw.get('dimensions')}` is a bare scalar "
                  f"with no source and no tier -- less provenance than "
                  f"`particle.density`. Declare a basis: "
                  f"{', '.join(DIM_BASES)}. "
                  f"See knowledge/wiki/concepts/dimensionality-has-no-default.md")]

    out: list = []
    where = f"{STRUCTURE_SECTION}.dim"
    basis = e.get("basis")
    if not isinstance(basis, str):   # a list/dict basis is unhashable -> `in` raises
        basis = None
    if basis not in DIM_BASES:
        out.append(I("error", f"{where}.basis",
                     f"must be one of {DIM_BASES} (got {basis!r}). "
                     + " . ".join(f"{k}={v}" for k, v in DIM_BASIS_MEANING.items())))

    #  ★ `given` owes one more thing: WHICH kind of input, and the obligation
    #    that kind carries. See DIM_GIVEN_KINDS.
    if basis == "given":
        kind = e.get("source_kind")
        if kind is None:
            out.append(I("error", f"{where}.source_kind",
                         "missing -- `basis: given` says the input states the "
                         "dimension, and this says WHICH input. One of "
                         + " . ".join(f"{k}={v}" for k, v in
                                      DIM_GIVEN_MEANING.items())))
        elif kind not in DIM_GIVEN_KINDS:
            out.append(I("error", f"{where}.source_kind",
                         f"{kind!r} is not one of {DIM_GIVEN_KINDS}"))
        else:
            for f in DIM_GIVEN_REQUIRES[kind]:
                if _blank(e.get(f)):
                    out.append(I("error", f"{where}.{f}",
                                 f"missing -- required when source_kind is "
                                 f"{kind!r} ({DIM_GIVEN_MEANING[kind]})"))
            if kind == "artefact":
                rel = e.get("source_file")
                if not _blank(rel):
                    #  ⚠ `ROOT` does not exist in this module. The first version
                    #    of this line read `s.path.parent if s.path else ROOT`,
                    #    a latent NameError that the eight real cases never
                    #    reach because `s.path` is always set for them.
                    repo = pathlib.Path(__file__).resolve().parent.parent
                    base = s.path.parent if s.path else repo
                    if not base.exists():
                        #  ★ `run.execute` rebuilds a PhysicalSystem around the
                        #    SPEC's system document, and its `path` does not
                        #    point at the case directory -- so existence cannot
                        #    be checked there. Say so. A check that cannot run
                        #    is not a check that passed, which is why the phi
                        #    closure and the `stated_in_source` cross-check
                        #    already warn on that path instead of going silent.
                        out.append(I(
                            "warn", f"{where}.source_file",
                            f"{rel!r} recorded, NOT verified: the case directory "
                            f"{str(base)!r} is not present here, so the file's "
                            f"existence cannot be checked. `bdbot.cli system "
                            f"check` is the full gate"))
                    #  resolve against the case directory first, then the repo
                    #  root -- `knowledge/source/papers/x.md` is a legitimate
                    #  artefact and lives at the root.
                    elif not [q for q in (base / str(rel), repo / str(rel))
                              if q.exists()]:
                        out.append(I(
                            "error", f"{where}.source_file",
                            f"{rel!r} does not exist, relative to the case "
                            f"directory or the repo root. `source_kind: "
                            f"artefact` claims a file states the dimension; a "
                            f"path that resolves to nothing is the one thing "
                            f"this check exists to catch"))

    need = DIM_ALWAYS + DIM_BASIS_REQUIRES.get(basis, ())
    for f in need:
        if f not in e:
            out.append(I("error", f"{where}.{f}",
                         "missing" + (f" -- required when basis is {basis!r} "
                                      f"({DIM_BASIS_MEANING[basis]})"
                                      if basis in DIM_BASIS_REQUIRES
                                      and f in DIM_BASIS_REQUIRES[basis] else "")))
        elif _blank(e[f]):      # was inline, and did not treat `{}` as blank

            out.append(I("error", f"{where}.{f}",
                         "present but empty. Blank must mean *forgotten*, never "
                         "*nothing to say*"
                         + (" -- an empty `what_would_change` means you do not know "
                            "whether this choice is safe, which is not the same as "
                            "it being safe."
                            if f == "what_would_change" else ".")))

    # the choice and the hashed physics field must be one number
    if "value" in e and "dimensions" in raw and e["value"] != raw["dimensions"]:
        out.append(I("error", f"{where}.value",
                     f"disagrees with `dimensions: {raw['dimensions']}` (structure "
                     f"says {e['value']}). `run_id` is hashed on `dimensions`, so "
                     f"this records one choice and runs another."))

    # the two bases that expire say so, every time they are read
    if basis == "inherited":
        tgt = e.get("compared_with") or []
        out.append(I("warn", where,
                     f"basis `inherited` EXPIRES if the comparison is dropped "
                     f"(target: {', '.join(map(str, tgt)) if tgt else '?'})"))
    elif basis == "sufficient":
        obs = e.get("checked_observables") or []
        out.append(I("warn", where,
                     f"basis `sufficient` EXPIRES when an observable is added -- "
                     f"checked against {len(obs)}: "
                     f"{', '.join(map(str, obs))[:60]}"))
    return out


def verify(s: PhysicalSystem, rtol: float = 1e-3) -> list:
    """Recompute the values written in `derived_scales` from the material formulas
    and compare (verified by reproduction rather than by provenance).
    """
    I = _intake.Issue
    ds = s.raw.get("derived_scales")
    if not ds:
        return []
    if not s.is_spherical:
        # The sphere formulas (3*pi*eta*d, kT/gamma) do not hold -> skip the
        # recomputation check and instead require **where it was derived from.**
        # The Perrin factors were not promoted into bdbot: they have appeared in
        # only one case so far (the "twice" rule in CLAUDE.md).
        src = str(ds.get("source", ""))
        if not src:
            return [I("error", "derived_scales",
                      f"shape is '{s.shape}', so the sphere formulas cannot recompute it. "
                      "Name the derivation script in `derived_scales.source` so it "
                      "stays reproducible.")]
        return [I("info", "derived_scales",
                  f"shape '{s.shape}' -- skipping the sphere recomputation check. "
                  f"basis: {src[:60]}")]
    b = s.bulk()
    if b is None:
        return [I("error", "derived_scales",
                  "could not recompute the material properties, so the derived values "
                  "were not compared (fix the unit/node errors above first).")]
    want = {"gamma": b["gamma"], "D_t": b["D_t"], "tau_B": b["tau_B"], "tau_p": b["tau_p"],
            "kT": b["kT"]}
    out = []
    for k, expect in want.items():
        if k not in ds or not isinstance(ds[k], dict):
            continue
        try:
            got = Q(ds[k]["value"], ds[k]["unit"])
            rel = abs(float((got - expect).to(expect.units).magnitude)
                      / float(expect.magnitude))
        except Exception as e:
            out.append(I("error", f"derived_scales.{k}", f"comparison failed: {e}"))
            continue
        if not math.isfinite(rel) or rel > rtol:
            out.append(I("error", f"derived_scales.{k}",
                         f"disagrees with the recomputation: written {got:~.5gP} vs computed {expect.to(got.units):~.5gP} "
                         f"({100*rel:.3f}%)"))
    return out


def render_check(s: PhysicalSystem) -> str:
    L: list[str] = []
    w = L.append
    w("=" * 78)
    w(f"system check — {s.path.parent.name}")
    w("=" * 78)
    if not s.raw:
        w("\n".join(str(i) for i in s.issues))
        return "\n".join(L)

    n_err = len(s.errors)
    n_warn = len([i for i in s.issues if i.level == "warn"])
    w(f"{s.label}   {s.dim}D   schema: {n_err} error(s) . {n_warn} warning(s)")
    if s.raw.get("derived_from"):
        w(f"  source (L0): {s.raw['derived_from']}")
    if s.issues:
        w("")
        for i in s.issues:
            w(str(i))

    w("")
    w("PHYSICAL SYSTEM (SI)")
    for name, p in s.core.items():
        w(f"  {name:<8} = {str(f'{p.value:~.4gP}'):<18} [tier {p.tier}] {p.source[:40]}")
    extra = s._extra_provenanced()
    for path, p in extra.items():
        w(f"  {path:<28} = {str(f'{p.value:~.4gP}')[:14]:<14} [tier {p.tier}]")

    b = s.bulk()
    if b is not None:
        w("")
        w("DERIVED PROPERTIES (recomputed)")
        w(f"  γ = {b['gamma']:~.4eP}   D_t = {b['D_t'].to('um^2/s'):~.4fP}   "
          f"τ_B = {b['tau_B']:~.4gP}   τ_p = {b['tau_p'].to('us'):~.3fP}")

    w("")
    w("TIER DISTRIBUTION")
    for t, names in s.tiers().items():
        w(f"  tier {t} ({TIER_MEANING[t]:<44}) {len(names):>2}  {', '.join(names[:5])}"
          f"{' …' if len(names) > 5 else ''}")

    nv = s.raw.get("not_verified")
    if nv:
        w("")
        w("NOT CONFIRMED (rule 3)")
        for x in nv:
            w(f"  · {str(x).splitlines()[0][:70]}")

    w("")
    w("=" * 78)
    if n_err:
        w(f"VERDICT: FAIL -- {n_err} error(s). Not advancing to L3 (non-dimensionalization).")
    else:
        w("VERDICT: READY -- L3 (non-dimensionalization) can proceed.")
        if n_warn:
            w(f"         ({n_warn} warning(s) -- check the tier approvals)")
    w("=" * 78)
    return "\n".join(L)


__all__ = [
           "DIM_GIVEN_KINDS", "DIM_GIVEN_MEANING", "DIM_GIVEN_REQUIRES",
           "SCHEMA", "PhysicalSystem", "load", "validate", "check_dim", "check_size",
           "verify", "render_check",
           "REQUIRED_TOP", "OPTIONAL_TOP", "CORE_PROVENANCED", "DERIVED_SECTIONS",
           "TIER_MEANING", "STRUCTURE_SECTION", "DIM_BASES",
           "DIM_BASIS_MEANING", "DIM_BASIS_REQUIRES", "DIM_ALWAYS",
           "SIZE_ROLES", "SIZE_ROLE_MEANING", "SIZE_FIXED_BY",
           "SIZE_ALWAYS", "SIZE_SYSTEM_REQUIRES", "SIZE_PINS", "SIZE_FS_STATUS",
           "SIZE_STATED", "SIZE_APPROVAL", "SIZE_NULLABLE", "GEOMETRY_PHI_KEYS",
           "check_size", "check_dim",
           "SIZE_REPLICATION_REQUIRES"]
