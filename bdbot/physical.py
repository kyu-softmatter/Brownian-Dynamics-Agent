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
        return int(self.raw.get("dimensions", 0))

    @property
    def shape(self):
        """`particle.shape` -- treated as a sphere when absent (as in the first two cases)."""
        return (self.raw.get("particle") or {}).get("shape")

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
            if top in DERIVED_SECTIONS:
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
        if top in DERIVED_SECTIONS:
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

    # 6. the tier approval gate
    low = [n for n, p in {**s.core, **s._extra_provenanced()}.items() if p.tier >= 2]
    if low:
        out.append(I("warn", "tier", f"{len(low)} value(s) at tier >= 2 (unverified): {', '.join(low[:6])}"
                                     f"{' ...' if len(low) > 6 else ''} -- needs human approval"))
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


__all__ = ["SCHEMA", "PhysicalSystem", "load", "validate", "verify", "render_check",
           "REQUIRED_TOP", "OPTIONAL_TOP", "CORE_PROVENANCED", "DERIVED_SECTIONS",
           "TIER_MEANING"]
