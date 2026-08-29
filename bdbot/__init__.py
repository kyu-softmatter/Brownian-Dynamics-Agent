"""bdbot -- the shared half of the Brownian dynamics pipeline (23 modules).

Everything here **actually appeared twice** across at least two cases. The basis
for each promotion is the comparison table in skill `bd-physics` section 6.3.
(`pairpot`, `traps`, `lockin` and `run` were promoted later, as a second or third
case appeared -- each module docstring records which.)

    cli          * the front end -- status . intake . system . nondim . run
    constants    ** k_B and the water eta(T)/rho(T) tables -- the SINGLE source
                 of truth, shared with `simbot.units`. Pure floats, no pint
    dt           ** the `dt` gate equations -- displacement gates, exact
                 Euler-Maruyama bias, and the timescale-ratio gate side by
                 side. SINGLE source of truth, shared with `simbot.nondim`
    intake       L0 Observation schema plus checks (derived from reading 5 real sketches)
    interactions catalogue of colloidal interactions plus a recommender (for when U_ij is blank)
    physical     L2 PhysicalSystem loader plus tier / derived_from / recomputation checks
    units        one pint registry (mixing registries makes pint refuse)
    provenance   Provenanced -- every number carries a source and a tier
    materials    gamma=3*pi*eta*d, D_t=kT/gamma, tau_B=d^2/D_t, m, tau_p (sphere in a Newtonian fluid)
    pairpot      numerics of a soft repulsive pair -- U, U'', closest approach r_min
                 (the physics that sets dt)
    scales       ScaleLedger -- the length/time/energy ledger, references, basis,
                 and required-role completeness
    nondim       * the L3 NondimSpec -- the only contract between L2 and L4.
                 Ledger, dimensionless groups, checks, inversion
    checks       Check(model/integration/geometry/statistics) + hard/soft verdict
                 + dt = 10^-2 * gamma/(local stiffness)
    traps        one harmonic trap expressing fixed, constant-velocity and
                 oscillatory driving (promoted after three appearances)
    sim          2D frame . Brownian integrator . GSD writer . seed handling . minimum image
    run          * the L5 assembly contract (`@builder`/`Build`) + the L6 execution
                 loop + L7 artefact storage. A case supplies only
                 `build(spec) -> Build`; the equilibration and production loops,
                 the guard calls and saving metrics.json are common. **All 8 cases**
                 use this contract. Three of them (`trap-2d-5um`,
                 `soft-r3-2d-A-sweep`, `abp-rod-2d-run-flip`) were re-run during
                 migration and compared against the old results -- the observables
                 matched to display precision (15 decimal places), meaning the
                 refactor did not change the physics.
    health       the L4 numerical-health verdict -- Guard, judge, step_health
                 (which feeds back into L3)
    lockin       complex stiffness K*(omega) of an oscillatory-driven system
                 (promoted after two appearances)
    report       the DimensionlessReport renderer (the case supplies its own blocks)
    runid        content-addressed run_id plus re-run prevention
    metrics      the metrics.json schema (the post-mortem's only input)
    stats        block averaging . autocorrelation correction . unbiased autocorrelation

`cli` exists so the same command gives the same result outside a Claude Code
session, and so a hook has something to intercept.
`bdbot.cli` does not import the heavy dependencies (hoomd, freud) -- the front end
stays fast when only reading and specifying. `run` defers its hoomd import into
`execute()` for the same reason.

** marks a module shared with `simbot/`. The dependency is one-way,
`simbot -> bdbot`, and it has to stay that way: `simbot.viz` pulls matplotlib.

**Submodules load lazily** (PEP 562). `from bdbot import checks as C` and
`bdbot.sim.frame_2d` behave exactly as before, but `import bdbot.constants`
no longer drags in pint and numpy: 190 ms -> 1 ms, measured 2026-08-29. That
is what lets `simbot.units` -- which every simbot module imports -- share this
package's constants without paying for pint.

**Deliberately not here** (appeared only once, or differs per system):
    the equilibrium indicator . observables . verification strategy . choice of the
    governing timescale . initial placement . the sampling loop
    -> those stay in the case scripts. If a third case needs one, it gets promoted
    then.

The absolute rules are in [CLAUDE.md](../CLAUDE.md), the physics procedure is in
skill `bd-physics`, and the HOOMD traps are in skill `bd-hoomd`.
"""
# ── lazy submodule / name resolution (PEP 562) ──────────────────────────────
#  Why lazy: `simbot.units` imports `bdbot.constants`, and `simbot` is imported by
#  the S1-S8 front end where a 190 ms package import is felt. Eagerly importing
#  every submodule here also contradicted this package's own stated goal that
#  "`bdbot.cli` does not import the heavy dependencies".
#  ⚠ Adding a module means adding it to `_SUBMODULES`, and a re-exported name
#    means adding it to `_NAMES`. `tests/test_bdbot_lazy_api.py` fails if the two
#    lists drift from what is actually importable -- an unwired list is worse than
#    no list (docs/02-verification.md section 6).
import importlib as _importlib

_SUBMODULES = (
    "checks", "cli", "constants", "dt", "health", "intake", "interactions",
    "lockin", "materials", "metrics", "nondim", "pairpot", "physical",
    "provenance", "report", "run", "runid", "scales", "sim", "stats", "traps",
    "units",
)

#  name -> the submodule that owns it
_NAMES = {
    "u": "units", "Q": "units", "kB": "units",
    "Provenanced": "provenance", "load_node": "provenance",
    "ScaleLedger": "scales", "Scale": "scales", "thermal_reference": "scales",
    "NondimSpec": "nondim", "Reference": "nondim", "Group": "nondim",
    "Check": "checks", "GATE": "checks", "verdict": "checks",
    "relaxation_time": "checks", "dt_from_gate": "checks",
    "dt_from_bias": "checks", "bias_from_dt": "checks",
}

__version__ = "0.1.0"

__all__ = [
    "u", "Q", "kB",
    "Provenanced", "load_node",
    "ScaleLedger", "Scale", "thermal_reference",
    "NondimSpec", "Reference", "Group",
    "Check", "GATE", "verdict", "relaxation_time", "dt_from_gate", "dt_from_bias",
    "bias_from_dt",
    *_SUBMODULES,
]


#  ⚠ `__name__` and not a literal, so renaming the package keeps working -- but
#    strip a trailing `.__init__`: `import bdbot.__init__` is legal and creates a
#    second module object whose `__name__` is `bdbot.__init__`, and the relative
#    import then fails with the misleading `'bdbot.__init__' is not a package`.
#    Nothing in this repo imports it that way; the two lines are here so that if
#    something ever does, the failure is not a red herring. (Found while auditing
#    this very refactor -- the audit script derived module names from file paths.)
_PKG = __name__.removesuffix(".__init__")


def __getattr__(name: str):
    if name in _SUBMODULES:
        mod = _importlib.import_module(f".{name}", _PKG)
        globals()[name] = mod                     # cache: __getattr__ runs once
        return mod
    if name in _NAMES:
        obj = getattr(_importlib.import_module(f".{_NAMES[name]}", _PKG), name)
        globals()[name] = obj
        return obj
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return sorted(set(__all__) | set(globals()))
