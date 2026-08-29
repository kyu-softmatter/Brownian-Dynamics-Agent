"""bdbot -- the shared half of the Brownian dynamics pipeline (21 modules).

Everything here **actually appeared twice** across at least two cases. The basis
for each promotion is the comparison table in skill `bd-physics` section 6.3.
(`pairpot`, `traps`, `lockin` and `run` were promoted later, as a second or third
case appeared -- each module docstring records which.)

    cli          * the front end -- status . intake . system . nondim . run
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

**Deliberately not here** (appeared only once, or differs per system):
    the equilibrium indicator . observables . verification strategy . choice of the
    governing timescale . initial placement . the sampling loop
    -> those stay in the case scripts. If a third case needs one, it gets promoted
    then.

The absolute rules are in [CLAUDE.md](../CLAUDE.md), the physics procedure is in
skill `bd-physics`, and the HOOMD traps are in skill `bd-hoomd`.
"""
from . import (checks, intake, interactions, materials, metrics, nondim, pairpot, physical,
               provenance, report, runid, scales, sim, stats)
from .checks import GATE, Check, bias_from_dt, dt_from_bias, dt_from_gate, relaxation_time, verdict
from .nondim import Group, NondimSpec, Reference
from .provenance import Provenanced, load_node
from .scales import Scale, ScaleLedger, thermal_reference
from .units import Q, kB, u

__version__ = "0.1.0"

__all__ = [
    "u", "Q", "kB",
    "Provenanced", "load_node",
    "ScaleLedger", "Scale", "thermal_reference",
    "NondimSpec", "Reference", "Group",
    "Check", "GATE", "verdict", "relaxation_time", "dt_from_gate", "dt_from_bias",
    "bias_from_dt",
    "units", "provenance", "materials", "scales", "nondim", "checks", "report", "runid",
    "metrics", "stats", "sim", "intake", "physical", "interactions", "pairpot",
]
