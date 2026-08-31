"""S3 — the system-specification data model + validity checks. 0 lines of LLM.

Three things this module enforces:

1. **Every physical quantity has a provenance.** A bare float cannot be entered
   without a `Quantity`. A number with no basis cannot later answer "where did
   this value come from", and then neither sensitivity analysis nor reproduction
   is possible.

2. **Derived values are recomputed, not stored.** If a `derived:` block is in the
   file it is compared against the recomputation to catch a mismatch. A
   hand-edited derived value is quietly wrong -- on 2026-07-28 a 4th-digit error
   in `kT(293.15 K)` actually happened.

3. **Turning a gate off requires writing down why.** And only registered gate
   names may be used -- a mistyped gate name becomes **a check that never runs
   once**. Which gates are on is decided by the (system × target dynamics) card:
   `knowledge/wiki/systems/_index.md`
"""
from __future__ import annotations

import math
from dataclasses import MISSING as _MISSING
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path

import yaml

from .io import sha256_payload
from .units import kT_si, stokes_drag_si, stokes_einstein_D_si
from bdbot.constants import sphere_mass_si   # one definition, shared with bdbot.materials

# =============================================================================
# provenance
# =============================================================================
PROVENANCE_KINDS: frozenset[str] = frozenset({
    "from_drawing",    # read directly off the drawing
    "observation",     # read directly off the drawing/source (from_drawing's
                       # modality-general form)
    "inference",       # derived from the source + physics knowledge
    "assumed",         # filled in because the source lacks it → an S7b sensitivity
                       # target
    "derived",         # computed from other fields (the return of a simbot call)
    "rule",            # derived from policy (config/run_policy.yaml)
    "from_knowledge",  # a knowledge/wiki entry
    "from_paper",      # a knowledge/source/papers distillation
    "measured",        # an experimental value
    "user",            # a human set this directly in this run (session `set`). Not
                       # an assumption -- excluded from what S7b shakes
})

# master_plan §12.2 — these provenances cannot be filled in by a cheap model
LLM_RESTRICTED: frozenset[str] = frozenset({"inference", "assumed"})
CHEAP_MODELS: frozenset[str] = frozenset({"haiku", "sonnet"})

CONFIDENCE_LEVELS: frozenset[str] = frozenset({"high", "medium", "low", ""})


@dataclass
class Quantity:
    """A physical-quantity wrapper. **One value, one basis.**

    `value` is a float by default but an int, a string or a list is allowed too
    (`dim=2`, `boundary="periodic"`, `center=[0,0,0]`).
    """

    value: float | int | str | bool | list
    unit: str = ""
    provenance: str = "assumed"
    basis: str = ""
    confidence: str = ""
    ambiguity: str = ""          # ambiguity id from 01_intake.md (A1, A2 …)
    sensitivity: str = ""        # none | low | high — S7b's result written back
    affects: list[str] = field(default_factory=list)
    written_by: str = ""         # for the model-tiering check (§12.2)

    def __post_init__(self) -> None:
        if self.provenance not in PROVENANCE_KINDS:
            raise ValueError(
                f"provenance {self.provenance!r} is not registered. "
                f"allowed: {sorted(PROVENANCE_KINDS)}")
        if self.confidence not in CONFIDENCE_LEVELS:
            raise ValueError(f"confidence {self.confidence!r} must be one of "
                             f"{sorted(CONFIDENCE_LEVELS)}")

    @property
    def si(self) -> float:
        """The numeric value. A string or list raises — that blocks it from being
        used in unit arithmetic."""
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)):
            raise TypeError(f"{self.value!r} is not numeric (unit={self.unit!r})")
        return float(self.value)

    def problems(self) -> list[str]:
        """The convention violations of this one value."""
        out = []
        if not str(self.basis).strip():
            out.append("basis is empty — a number with no basis allows neither "
                       "sensitivity analysis nor reproduction")
        if self.provenance in LLM_RESTRICTED and not self.confidence:
            out.append(f"provenance={self.provenance} with no confidence "
                       f"(an inference or assumption has to state its confidence)")
        if (self.provenance in LLM_RESTRICTED
                and self.written_by.lower() in CHEAP_MODELS):
            out.append(f"provenance={self.provenance} was filled in by "
                       f"{self.written_by} — master_plan §12.2 violation "
                       f"(Opus only)")
        return out


def Q(value, unit: str = "", provenance: str = "assumed", basis: str = "", **kw):
    """A short `Quantity` constructor. For when the YAML is not hand-written."""
    return Quantity(value=value, unit=unit, provenance=provenance, basis=basis, **kw)


# =============================================================================
# gates — the card turns them on and off
# =============================================================================
GATE_STATUSES: frozenset[str] = frozenset({"required", "pass", "fail", "off",
                                           "applicable", "unknown"})

# Registered gate names. A mistyped name becomes "a check that never runs", so it
# is rejected.
KNOWN_GATES: frozenset[str] = frozenset({
    # premises
    "overdamped", "stokes_reynolds", "hydrodynamics_neglected",
    # equilibrium and structure
    "equilibration_detection", "equipartition", "configurational_temperature",
    "self_consistency_D", "polydispersity",
    # numerics
    "em_bias_reproduced", "dt_over_tau_trap", "thermal_displacement",
    "force_displacement", "active_displacement", "step_displacement_vs_sigma",
    # geometry
    "box_much_larger_than_l_trap", "r_cut_le_half_box", "finite_size_L",
    "persistence_length_vs_box", "packing_fraction", "debye_length_consistency",
})


@dataclass
class Gate:
    status: str = "unknown"
    reason: str = ""

    def __post_init__(self) -> None:
        if self.status not in GATE_STATUSES:
            raise ValueError(f"gate status {self.status!r} must be one of "
                             f"{sorted(GATE_STATUSES)}")


# =============================================================================
# components
# =============================================================================
@dataclass
class Species:
    name: str
    n_simulated: Quantity
    radius_si: Quantity
    density_si: Quantity | None = None
    n_physical: Quantity | None = None
    charge: Quantity | None = None
    active: bool = False

    @property
    def sigma_si(self) -> float:
        """The representative diameter. ⚠ `σ = 2a`. Confuse it with the radius and
        every timescale is wrong by 2x."""
        return 2.0 * self.radius_si.si

    def mass_si(self) -> float | None:
        """Sphere-assumed mass. `None` without a density (no overdamped check then).

        ★ The expression is `bdbot.constants.sphere_mass_si`, shared with
          `bdbot.materials.sphere_mass`. It used to be written here as
          `rho*(4/3)*pi*a**3`, which is the same number on paper and **1 ULP
          different in floating point** from bdbot's `rho*(pi/6)*d**3` -- so the
          two halves could never compare equal. Merged 2026-08-29.
          ⚠ The kernel takes a **diameter**; this class stores a radius.
        """
        if self.density_si is None:
            return None
        return sphere_mass_si(self.density_si.si, 2.0 * self.radius_si.si)


@dataclass
class Medium:
    T_si: Quantity
    eta_si: Quantity
    rho_fluid_si: Quantity | None = None
    species: Quantity | None = None      # "water" and the like


@dataclass
class Geometry:
    dim: Quantity
    boundary: Quantity
    box_si: Quantity | None = None       # [Lx, Ly, Lz] stated explicitly
    box_over_ref: Quantity | None = None  # multiples of the card's reference length
                                          # (the convention for trap systems)

    @property
    def d(self) -> int:
        return int(self.dim.value)


@dataclass
class PairInteraction:
    type_a: str
    type_b: str
    potential: str                       # wca | lj | yukawa | morse | dlvo
    params: dict = field(default_factory=dict)
    r_cut_si: Quantity | None = None


@dataclass
class BondInteraction:
    """One bond kind. Maps to HOOMD `md.bond.Harmonic`.

    `params`: `k_si` [N/m = J/m²] · `r0_si` [m].

    ⚠ **Do not merge this into one class with the angle (`AngleInteraction`).** A
      bond's `k` is N/m and an angle's `k` is J/rad² -- put two different
      dimensions in a same-named field and building `λ_max` goes quietly wrong.
      The stability gate uses that value.

    ⚠ Do not use a distance constraint (`constrain.Distance`) together with
      `Brownian`. Basis:
      `knowledge/wiki/findings/dead-end-distance-constraint-with-brownian.md`
    """

    name: str = "backbone"               # the HOOMD bond type name
    potential: str = "harmonic"
    params: dict = field(default_factory=dict)


@dataclass
class AngleInteraction:
    """One angle kind. Maps to HOOMD `md.angle.Harmonic`.

    `params`: `k_si` [J/rad²] · `t0_rad` [rad].
    """

    name: str = "backbone"               # the HOOMD angle type name
    potential: str = "harmonic"
    params: dict = field(default_factory=dict)


@dataclass
class ExternalField:
    kind: str                            # harmonic_trap | gravity | shear | electric
    params: dict[str, Quantity] = field(default_factory=dict)
    implementation: str = ""
    note: str = ""


@dataclass
class Friction:
    model: str = "stokes_infinite_medium"
    gamma_si: Quantity | None = None     # if absent, derived as 6πηa
    wall_correction: str = "none"
    note: str = ""


@dataclass
class Timing:
    equil_in_tau: Quantity | None = None
    prod_in_tau: Quantity | None = None
    sample_interval_in_tau: Quantity | None = None
    target_precision: Quantity | None = None


@dataclass
class Numerics:
    dt_star: Quantity | None = None
    seed_base: int = 1
    n_seeds: Quantity | None = None
    integrator: str = "hoomd.md.methods.Brownian"
    scheme: str = "euler_maruyama"
    noise_distribution: str = "uniform"  # ★ NOT Gaussian. findings §2


# =============================================================================
# SystemSpec
# =============================================================================
@dataclass
class SystemSpec:
    """The complete specification in physical units. **Inputs only** — derived
    values are produced by `derive()`."""

    card: str
    question: str
    geometry: Geometry
    species: list[Species]
    medium: Medium
    friction: Friction = field(default_factory=Friction)
    pair: list[PairInteraction] = field(default_factory=list)
    bonds: list[BondInteraction] = field(default_factory=list)
    angles: list[AngleInteraction] = field(default_factory=list)
    external: list[ExternalField] = field(default_factory=list)
    timing: Timing = field(default_factory=Timing)
    numerics: Numerics = field(default_factory=Numerics)
    gates: dict[str, Gate] = field(default_factory=dict)
    tier: str = ""
    notes: list[str] = field(default_factory=list)

    # --- convenience accessors ---
    @property
    def primary(self) -> Species:
        return self.species[0]

    def trap(self) -> ExternalField | None:
        for e in self.external:
            if e.kind == "harmonic_trap":
                return e
        return None

    @property
    def has_neighbor_interaction(self) -> bool:
        """Is there something to overlap with (a pair) or **to be bonded to** (a
        bond or angle)?

        ★ This is the displacement gate's activation condition. Deciding it with
          `bool(spec.pair)` turns the gate quietly off in a bond-only system (a
          colloidal chain) -- and it actually was off that way. Basis:
          `knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md`
        """
        return bool(self.pair or self.bonds or self.angles)

    def bond_stiffness_si(self) -> float | None:
        """The stiffest bond spring [N/m]. Must be the **maximum** — stability is
        set by the worst mode."""
        ks = [b.params["k_si"].si for b in self.bonds if "k_si" in b.params]
        return max(ks) if ks else None

    def angle_stiffness_si(self) -> float | None:
        """The stiffest angle spring [J/rad²]."""
        ks = [a.params["k_si"].si for a in self.angles if "k_si" in a.params]
        return max(ks) if ks else None

    def bond_length_si(self) -> float | None:
        """The shortest bond equilibrium length [m]. It is the lever arm when
        converting an angular stiffness into a transverse one, and shorter means
        stiffer, so the minimum is used."""
        rs = [b.params["r0_si"].si for b in self.bonds if "r0_si" in b.params]
        return min(rs) if rs else None

    def gamma_si(self) -> float:
        """Drag coefficient. The stated value if there is one, else Stokes `6πηa`."""
        if self.friction.gamma_si is not None:
            return self.friction.gamma_si.si
        return stokes_drag_si(self.medium.eta_si.si, self.primary.radius_si.si)

    def box_lengths_si(self, ref_length_si: float | None = None) -> list[float] | None:
        """Box edge lengths [m]. With only `box_over_ref`, `ref_length_si` is
        required."""
        if self.geometry.box_si is not None:
            return [float(x) for x in self.geometry.box_si.value]
        if self.geometry.box_over_ref is not None and ref_length_si is not None:
            L = self.geometry.box_over_ref.si * ref_length_si
            d = self.geometry.d
            return [L, L, (L if d == 3 else 0.0)]
        return None

    def hash(self) -> str:
        """The spec hash — `run_id` and the cache key. Derived values excluded."""
        return sha256_payload(to_dict(self))

    # --- serialization ---
    def to_yaml(self) -> str:
        return dump_yaml(self)

    @classmethod
    def from_yaml(cls, text: str) -> SystemSpec:
        return from_dict(cls, yaml.safe_load(text))

    @classmethod
    def load(cls, path: str | Path) -> SystemSpec:
        return cls.from_yaml(Path(path).read_text(encoding="utf-8"))

    def save(self, path: str | Path) -> Path:
        p = Path(path)
        p.write_text(self.to_yaml(), encoding="utf-8")
        return p


# =============================================================================
# predictions (S2) — the scientific claim that gets sealed
# =============================================================================
@dataclass
class PredictionItem:
    """One prediction in falsifiable form. All 4 elements are required
    (master_plan §S2-4)."""

    quantity: str
    value: float | str
    tolerance: str                  # "±1.5%" | "±0.03" | "p>0.05" | ">0.99"
    basis: str
    discriminates: str = ""
    unit: str = ""
    competing_value: float | None = None   # the competing hypothesis (for power)
    note: str = ""

    def problems(self) -> list[str]:
        out = []
        for name in ("tolerance", "basis"):
            if not str(getattr(self, name)).strip():
                out.append(f"{self.quantity}: {name} is empty")
        return out


@dataclass
class Prediction:
    items: list[PredictionItem]
    regimes: dict = field(default_factory=dict)
    alternatives: list[str] = field(default_factory=list)

    def problems(self) -> list[str]:
        if not self.items:
            return ["0 quantitative predictions — the S2 gate requires at least 1"]
        return [p for it in self.items for p in it.problems()]


# =============================================================================
# derived values — computed, never stored
# =============================================================================
def derive(spec: SystemSpec) -> dict[str, float]:
    """Every SI scale derived from the spec. **Every value is a function's return.**

    A derived value that is not here does not go in the report -- that removes the
    place a hand calculation could slip in.
    """
    sp = spec.primary
    T = spec.medium.T_si.si
    eta = spec.medium.eta_si.si
    a = sp.radius_si.si
    sigma = sp.sigma_si
    gamma = spec.gamma_si()
    kT = kT_si(T)
    D0 = stokes_einstein_D_si(T, gamma)

    out = {
        "kT_si": kT,
        "sigma_si": sigma,
        "gamma_si": gamma,
        "D0_si": D0,
        "tau_D_si": sigma**2 / D0,
    }
    mass = sp.mass_si()
    if mass is not None:
        out["mass_si"] = mass
        out["tau_inertial_si"] = mass / gamma

    trap = spec.trap()
    if trap is not None and "k_si" in trap.params:
        k = trap.params["k_si"].si
        out["k_si"] = k
        out["tau_trap_si"] = gamma / k
        out["l_trap_si"] = math.sqrt(kT / k)
        out["corner_freq_si"] = k / (2.0 * math.pi * gamma)
        out["var_per_component_si"] = kT / k
        out["msd_plateau_si"] = 2.0 * spec.geometry.d * kT / k
        out["k_star_sigma"] = k * sigma**2 / kT

    # --- bonds and angles → λ_max, the largest eigenvalue of the stiffness matrix -
    #  Approximation for a chain's stiffness matrix: 4k for a 1D spring chain, and
    #  16k for bending because it is a 4th-order difference. An angular spring is
    #  in J/rad², so it is divided by the lever arm b² to become a transverse
    #  stiffness [N/m].
    #  This value sets **explicit Euler's stability limit** Δt ≤ 2γ/λ_max. It is
    #  the quantity that decides divergence, not accuracy, and the displacement
    #  gate cannot catch it (a straight chain has |F| = 0).
    #  Measured calibration (ratio 1.22–2.80):
    #  knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
    k_bond = spec.bond_stiffness_si()
    k_angle = spec.angle_stiffness_si()
    if k_bond is not None or k_angle is not None:
        b = spec.bond_length_si() or sigma
        out["bond_length_si"] = b
        lam = 0.0
        if k_bond is not None:
            out["k_bond_si"] = k_bond
            out["k_bond_star_sigma"] = k_bond * sigma**2 / kT
            out["tau_bond_si"] = gamma / k_bond
            lam += 4.0 * k_bond
        if k_angle is not None:
            out["k_angle_si"] = k_angle
            out["k_angle_star"] = k_angle / kT        # J/rad² uses no length
            out["tau_angle_si"] = gamma * b**2 / k_angle
            lam += 16.0 * k_angle / b**2
        out["lambda_max_si"] = lam
        out["tau_stiff_si"] = gamma / lam
    return out


def reference_length_si(spec: SystemSpec, derived: dict | None = None) -> float:
    """The card's reference length [m]. `ℓ_trap` for a trap system, `σ` otherwise.

    The choice is owned by the (system × target dynamics) card —
    CLAUDE.md §non-dimensionalization convention.
    """
    d = derived if derived is not None else derive(spec)
    return d.get("l_trap_si", d["sigma_si"])


# =============================================================================
# validity checks
# =============================================================================
@dataclass
class Check:
    """One check.

    `declared` is **not a result** -- it means the card declared this gate on but
    it is a quantity S3 cannot compute (equipartition, EM bias and the like are
    decided by S7). Write `declared` as `pass` and a pass stamp gets applied that
    no human ever looked at.
    """

    name: str
    status: str            # pass | fail | off | na | warn | declared
    detail: str = ""
    value: float | None = None
    threshold: float | None = None


@dataclass
class SpecReport:
    """The check results. **It does not judge** — `ok` only means whether there is
    a convention violation."""

    checks: list[Check]
    problems: list[str]
    derived: dict[str, float]

    @property
    def ok(self) -> bool:
        return not self.problems and not any(c.status == "fail" for c in self.checks)

    def failed(self) -> list[Check]:
        return [c for c in self.checks if c.status == "fail"]

    def deferred(self) -> list[Check]:
        """Gates S7 has to decide — S3 cannot compute them."""
        return [c for c in self.checks if c.status == "declared"]

    def table(self) -> str:
        """A markdown table. Used by `03_spec_rationale.md` and `REPORT.md`."""
        rows = ["| check | status | value | threshold | note |",
                "|---|---|---|---|---|"]
        mark = {"pass": "✅ pass", "fail": "❌ **fail**", "off": "— off",
                "na": "— n/a", "warn": "⚠️ warn", "declared": "⏳ S7 decides"}
        for c in self.checks:
            v = "" if c.value is None else f"`{c.value:.4g}`"
            t = "" if c.threshold is None else f"`{c.threshold:.4g}`"
            rows.append(f"| `{c.name}` | {mark.get(c.status, c.status)} | {v} | {t} "
                        f"| {c.detail} |")
        return "\n".join(rows)


def _iter_quantities(obj, path: str = "") -> list[tuple[str, Quantity]]:
    """Pull every `Quantity` out of a nested structure as (path, value)."""
    out: list[tuple[str, Quantity]] = []
    if isinstance(obj, Quantity):
        out.append((path, obj))
    elif is_dataclass(obj):
        for f in fields(obj):
            out += _iter_quantities(getattr(obj, f.name), f"{path}.{f.name}".lstrip("."))
    elif isinstance(obj, dict):
        for k, v in obj.items():
            out += _iter_quantities(v, f"{path}.{k}".lstrip("."))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            out += _iter_quantities(v, f"{path}[{i}]")
    return out


def packing_fraction(spec: SystemSpec, box_si: list[float] | None) -> float | None:
    """φ (3D volume fraction / 2D area fraction). `None` if the box is unknown."""
    if not box_si:
        return None
    d = spec.geometry.d
    if d == 2:
        area = box_si[0] * box_si[1]
        tot = sum(s.n_simulated.si * math.pi * s.radius_si.si**2 for s in spec.species)
        return tot / area if area > 0 else None
    vol = box_si[0] * box_si[1] * box_si[2]
    tot = sum(s.n_simulated.si * (4.0 / 3.0) * math.pi * s.radius_si.si**3
              for s in spec.species)
    return tot / vol if vol > 0 else None


def validate(spec: SystemSpec, *, stored_derived: dict | None = None,
             rel_tol: float = 1e-3) -> SpecReport:
    """Convention violations + physical validity checks.

    Args:
        stored_derived: the derived values as written in the file. Given them, they
            are compared against the recomputation (the only way to catch a
            hand-edited derived value).
        rel_tol: relative tolerance for that comparison. Documents record 4 to 5
            significant figures, so the default is `1e-3`.
    """
    problems: list[str] = []
    d = derive(spec)
    computed: dict[str, Check] = {}

    def put(c: Check) -> None:
        computed[c.name] = c

    # --- 1. provenance completeness -------------------------------------------
    for path, q in _iter_quantities(spec):
        problems += [f"{path}: {p}" for p in q.problems()]

    # --- 2. the form of a gate declaration ------------------------------------
    unknown = sorted(set(spec.gates) - KNOWN_GATES)
    if unknown:
        problems.append(
            f"unregistered gate name {unknown} — a typo becomes 'a check that "
            f"never runs'. For a new gate, add it to KNOWN_GATES first")
    for name, g in spec.gates.items():
        if g.status == "off" and not g.reason.strip():
            problems.append(f"gate {name} was turned off with no reason — "
                            f"the basis for turning it off belongs in the card")

    # --- 3. overdamped --------------------------------------------------------
    tau_proc = d.get("tau_trap_si", d["tau_D_si"])
    if "tau_inertial_si" in d:
        ratio = d["tau_inertial_si"] / tau_proc
        put(Check("overdamped", "pass" if ratio < 1e-2 else "fail", value=ratio,
                  threshold=1e-2,
                  detail=("τ_i/τ_process ≪ 1 — the BD premise holds"
                          if ratio < 1e-2
                          else "inertia is not negligible → consider Langevin")))
    else:
        put(Check("overdamped", "na",
                  detail="no density, so the mass is unknown — cannot check"))

    # --- 4. Reynolds --------------------------------------------------------
    rho_f = spec.medium.rho_fluid_si
    if rho_f is not None:
        # characteristic velocity: the reference length per reference time
        v = reference_length_si(spec, d) / tau_proc
        Re = rho_f.si * v * spec.primary.radius_si.si / spec.medium.eta_si.si
        put(Check("stokes_reynolds", "pass" if Re < 1e-2 else "fail",
                  value=Re, threshold=1e-2,
                  detail="Re ≪ 1 — Stokes drag holds" if Re < 1e-2
                  else "inertial fluid effects suspected"))
    else:
        put(Check("stokes_reynolds", "na", detail="no ρ_fluid"))

    # --- 5. box · φ · r_cut ---------------------------------------------------
    ref = reference_length_si(spec, d)
    box = spec.box_lengths_si(ref)
    if box is None:
        problems.append("box size unknown — box_si or box_over_ref is required")
    else:
        n_ref = min(x for x in box[:spec.geometry.d]) / ref
        put(Check("box_much_larger_than_l_trap",
                  "pass" if n_ref >= 10 else "fail", value=n_ref, threshold=10.0,
                  detail=f"box is {n_ref:.3g}x the reference length"))
        phi = packing_fraction(spec, box)
        if phi is not None:
            cap = 0.9 if spec.geometry.d == 2 else 0.64
            kind = "area" if spec.geometry.d == 2 else "volume"
            if not spec.pair:
                # ★ With no pair interaction there is no excluded volume → φ has
                #   no physical meaning. N particles in the same trap are
                #   **independent replicas**, not a suspension.
                #   Gating on φ here gives φ=4741 and a can-never-pass verdict --
                #   a problem that does not exist. Show the value, do not judge it.
                put(Check("packing_fraction", "off", value=phi, threshold=cap,
                          detail=f"{kind} fraction is meaningless — no pair "
                                 f"interaction, so no excluded volume (the "
                                 f"particles are independent replicas)"))
            else:
                put(Check("packing_fraction", "pass" if phi < cap else "fail",
                          value=phi, threshold=cap,
                          detail=f"{kind} fraction (cap = "
                                 f"{'2D limit' if spec.geometry.d == 2 else 'RCP'})"))
        if spec.pair:
            half = min(box[:spec.geometry.d]) / 2.0
            worst = max((p.r_cut_si.si for p in spec.pair
                         if p.r_cut_si is not None), default=None)
            if worst is None:
                problems.append("there is a pair interaction but no r_cut")
            else:
                put(Check("r_cut_le_half_box", "pass" if worst <= half else "fail",
                          value=worst, threshold=half,
                          detail="minimum-image convention"))
        else:
            put(Check("r_cut_le_half_box", "off", detail="no pair potential"))

    # --- 6. comparing the derived values --------------------------------------
    if stored_derived:
        bad, n_compared, skipped = [], 0, []
        for k, stored in stored_derived.items():
            if k not in d or not isinstance(stored, (int, float)):
                skipped.append(k)
                continue
            n_compared += 1
            rel = abs(float(stored) - d[k]) / max(abs(d[k]), 1e-300)
            if rel > rel_tol:
                bad.append(f"{k}: file {stored:.6g} vs recomputed {d[k]:.6g} "
                           f"(relative diff {rel:.2e})")
        problems += [f"derived-value mismatch — {b} — possibly a hand-edited value"
                     for b in bad]
        # ★ Report how many were compared. Reporting 0 comparisons as pass makes
        #   "it was checked" a false statement.
        detail = f"{n_compared} compared"
        if bad:
            detail += f" · {len(bad)} mismatched"
        if skipped:
            detail += f" · {len(skipped)} uncomparable ({', '.join(skipped[:3])}…)" \
                if len(skipped) > 3 else f" · uncomparable {skipped}"
        put(Check("derived_consistency",
                  "fail" if bad else ("na" if n_compared == 0 else "pass"),
                  value=float(n_compared), detail=detail))

    # --- 7. gates ∪ computed results -------------------------------------------
    # Computed results fill in the declarations. A gate that cannot be computed
    # stays `declared` and passes to S7 -- stamping pass here would be a pass
    # nobody looked at.
    checks: list[Check] = []
    for name, g in spec.gates.items():
        if name in computed:
            c = computed.pop(name)
            parts = [x for x in (c.detail, g.reason) if x]
            if len(parts) == 2 and parts[0] == parts[1]:
                parts = parts[:1]
            checks.append(Check(name, c.status, value=c.value,
                                threshold=c.threshold, detail=" · ".join(parts)))
        elif g.status == "off":
            checks.append(Check(name, "off", detail=g.reason))
        else:
            checks.append(Check(name, "declared",
                                detail=g.reason or f"card declares: {g.status}"))
    # computed but never declared by the card (this surfaces a missing declaration)
    for c in computed.values():
        checks.append(c)

    return SpecReport(checks=checks, problems=problems, derived=d)


# =============================================================================
# YAML serialization — the round-trip error must be 0
# =============================================================================
_QUANTITY_DEFAULTS = {f.name: f.default for f in fields(Quantity)
                      if f.name not in ("value",)}


def to_dict(obj):
    """A dataclass tree → a plain dict. A `Quantity` omits default-valued fields to
    stay short."""
    if isinstance(obj, Quantity):
        # provenance is **always written**, even at its default. Omitting it
        # because `assumed` is the default makes "I assumed it" and "I forgot to
        # write it" indistinguishable in the file -- and `assumed` is the field
        # that decides S7b's sensitivity target list.
        out = {"value": obj.value, "provenance": obj.provenance}
        if obj.unit:
            out = {"value": obj.value, "unit": obj.unit,
                   "provenance": obj.provenance}
        for name, default in _QUANTITY_DEFAULTS.items():
            if name in ("unit", "provenance"):
                continue
            v = getattr(obj, name)
            if name == "affects":
                if v:
                    out[name] = list(v)
                continue
            if v != default:
                out[name] = v
        return out
    if isinstance(obj, Gate):
        return {"status": obj.status, **({"reason": obj.reason} if obj.reason else {})}
    if is_dataclass(obj):
        out = {}
        for f in fields(obj):
            v = getattr(obj, f.name)
            if v is None:
                continue
            default = f.default if f.default is not _MISSING else None
            if isinstance(v, (list, dict)) and not v:
                continue
            if default is not None and v == default and not isinstance(v, Quantity):
                continue
            out[f.name] = to_dict(v)
        return out
    if isinstance(obj, dict):
        return {k: to_dict(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_dict(v) for v in obj]
    return obj


def _build_quantity(raw) -> Quantity:
    if isinstance(raw, Quantity):
        return raw
    if isinstance(raw, dict) and "value" in raw:
        return Quantity(**raw)
    # a bare value arrived — it has no provenance, so let it through as a violation
    return Quantity(value=raw, provenance="assumed", basis="")


def from_dict(cls, raw):
    """dict → dataclass. A `Quantity` field takes the `{value, provenance, …}` form."""
    if raw is None:
        return None
    if cls is Quantity:
        return _build_quantity(raw)
    if cls is Gate:
        return Gate(**raw) if isinstance(raw, dict) else Gate(status=str(raw))
    if not is_dataclass(cls):
        return raw

    kwargs = {}
    for f in fields(cls):
        if f.name not in raw:
            continue
        v = raw[f.name]
        ann = f.type if not isinstance(f.type, str) else f.type
        kwargs[f.name] = _coerce(f.name, ann, v, cls)
    return cls(**kwargs)


# field name → element type. An explicit table rather than parsing string annotations.
_LIST_FIELD_TYPES = {
    ("SystemSpec", "species"): Species,
    ("SystemSpec", "pair"): PairInteraction,
    ("SystemSpec", "bonds"): BondInteraction,
    ("SystemSpec", "angles"): AngleInteraction,
    ("SystemSpec", "external"): ExternalField,
    ("Prediction", "items"): PredictionItem,
}
_NESTED_FIELD_TYPES = {
    ("SystemSpec", "geometry"): Geometry,
    ("SystemSpec", "medium"): Medium,
    ("SystemSpec", "friction"): Friction,
    ("SystemSpec", "timing"): Timing,
    ("SystemSpec", "numerics"): Numerics,
}
# Fields that must be wrapped in a Quantity (by name). Values inside `params` are
# all Quantity too.
_QUANTITY_FIELDS = {
    "n_simulated", "n_physical", "radius_si", "density_si", "charge",
    "T_si", "eta_si", "rho_fluid_si", "species", "dim", "boundary", "box_si",
    "box_over_ref", "gamma_si", "r_cut_si", "dt_star", "n_seeds",
    "equil_in_tau", "prod_in_tau", "sample_interval_in_tau", "target_precision",
}


def _coerce(name: str, ann, v, owner):
    key = (owner.__name__, name)
    if key in _LIST_FIELD_TYPES:
        return [from_dict(_LIST_FIELD_TYPES[key], x) for x in v]
    if key in _NESTED_FIELD_TYPES:
        return from_dict(_NESTED_FIELD_TYPES[key], v)
    if name == "gates":
        return {k: from_dict(Gate, g) for k, g in v.items()}
    if name == "params" and isinstance(v, dict):
        # a trap's k_si and so on are Quantity; pure config like active_axes is not
        return {k: (_build_quantity(x) if isinstance(x, dict) and "value" in x else x)
                for k, x in v.items()}
    if name in _QUANTITY_FIELDS:
        return _build_quantity(v)
    return v


class _Dumper(yaml.SafeDumper):
    """Indentation a person can read. The same shape as a hand-written 03_spec.yaml."""

    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)


def dump_yaml(obj) -> str:
    return yaml.dump(to_dict(obj), Dumper=_Dumper, sort_keys=False,
                     allow_unicode=True, default_flow_style=False, width=100)


def load_prediction(path: str | Path) -> Prediction:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return from_dict(Prediction, raw)
