"""`metrics.json` -- the machine-readable run result. The **only** input to the
post-mortem and the knowledge feedback loop.

Both of the first two cases used the same top-level keys. `tools/postmortem.py`
depends on this schema, so changing a field means looking there too (backward
compatibility is handled by postmortem using `.get`).

Three fields added by the second case, all from the same lesson -- "what differs
per case is declared by the case":
  * `equilibration`  the time series used for the equilibrium verdict. The first
                     case's 'anchor displacement' cannot be used in a diffusive system
  * `checks[].hard`  the hard/soft distinction. In the first case everything passed,
                     so the distinction never surfaced
  * `numerics.stat_target_pct`  the statistical target. The first case's 0.5% must
                     not be carried to another system

WARNING: this file is written by `bdbot.run`, but `result.txt` -- which
`bdbot.cli status` counts runs by -- is written by the *case script*. A case that
never added that line has its runs reported as zero while its metrics.json files
sit on disk.
"""
from __future__ import annotations

import json
from pathlib import Path

SCHEMA = "bdbot.metrics/0.3"

# * The **role** of an observable. Without this distinction the verdict is wrong.
#
#   implementation_check  the prediction **follows analytically from the model I
#                         implemented.** Agreement means "the code is right" and is
#                         not a physical discovery (very nearly circular).
#                         A mismatch = **a bug** -> FAIL.
#   hypothesis            the prediction rests on **an assumption the simulation
#                         does not impose** (continuum approximation, dilute limit,
#                         effective medium, a literature model).
#                         A mismatch = **a result.** Not a FAIL -- it is the thing
#                         we wanted to know.
#   measurement           there is no prediction. The simulation is the answer.
#
# Why this is needed: the reason to compute these systems is that they **may
# differ** from the standard picture. If the verdict logic is "differs from the
# prediction, therefore FAIL", it calls discoveries failures.
ROLES = ("implementation_check", "hypothesis", "measurement")


# Composition level (rule 7 -- isolation). A single module and a full combination
# have different epistemic status.
#   module     the minimal configuration with only this module on -- the standard
#              theory *is* my model here, so implementation_check is legitimate
#   composite  several modules combined -- usually no analytic solution exists.
#              Applying a standard theory here is dangerous
SCOPES = ("module", "composite")


def observable(name, measured, predicted=None, unit="1", source="none",
               role="measurement", tol_pct=None, sigma=None, tol_sigma=None, note="",
               scope="composite", derivation="") -> dict:
    """One observable.

    Always think about `role` and set it (see ROLES above). The default is the most
    conservative, `measurement` (no verdict) -- so a discovery is never accidentally
    called a failure.
    `tol_pct` only means anything for `implementation_check`.

    * If the prediction is exactly 0, a percentage error is undefined (division by
      zero). Pass `sigma` (the measurement's standard error, e.g. the block-mean
      SEM) alongside and the verdict uses the z-score
      `err_sigma = (measured-predicted)/sigma` (default tolerance `tol_sigma=3.0`).
      Using `predicted=0` without `sigma` makes it undecidable (treated as
      measurement) -- `judge()` does that so it cannot be silently wrong.

    `scope` says whether this observable came from **a single module** or from **the
    full combination** (rule 7). The default `composite` is the conservative one.

    * `composite` + `implementation_check` requires a `derivation`.
      A combination usually has no analytic solution -- that is why we are
      simulating it. If you nevertheless want an implementation check, you have to
      write down **why that expression is also derivable for the combination**
      (usually a limit: dilute, linear response, short time). Without it, you are
      applying a standard theory to a combination, and then agreement reads as
      "verified" and disagreement as "the theory does not hold here" -- so
      **either way nothing is learned.**

    A prediction must be fixed **before the result is seen.** The structural
    guarantee is that the prediction function does not take simulation results as
    arguments -- the `analytic(ledger)` pattern in `cases/*.py`.
    """
    if role not in ROLES:
        raise ValueError(f"role must be one of {ROLES} (got: {role!r})")
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES} (got: {scope!r})")
    if scope == "composite" and role == "implementation_check" and not derivation:
        raise ValueError(
            f"[{name}] composite + implementation_check requires a derivation "
            f"(rule 7). State why that analytic expression is also derivable for "
            f"the combination -- usually a limit. If you cannot state one, "
            f"role='hypothesis' is the correct choice: a mismatch is then a result, "
            f"not a bug.")
    err = None
    if predicted not in (None, 0) and measured is not None:
        err = 100.0 * (float(measured) - float(predicted)) / abs(float(predicted))
    err_sigma = None
    if sigma and predicted is not None and measured is not None:
        err_sigma = (float(measured) - float(predicted)) / float(sigma)
    return {"name": name, "measured": None if measured is None else float(measured),
            "predicted": None if predicted is None else float(predicted),
            "unit": unit, "err_pct": err, "err_sigma": err_sigma,
            "sigma": None if sigma is None else float(sigma),
            "prediction_source": source,
            "role": role, "scope": scope, "derivation": derivation,
            "tol_pct": tol_pct, "tol_sigma": tol_sigma, "note": note}


def judge(observables) -> tuple:
    """(verdict, failed implementation checks, hypothesis deviations, measurements
    only) -- handled differently per role.

    * A `hypothesis` deviation is not a FAIL. It is **a result** to report.
    """
    bad_impl, dev_hypo, meas = [], [], []
    for o in observables:
        role, err = o.get("role", "measurement"), o.get("err_pct")
        err_sigma = o.get("err_sigma")
        if role == "implementation_check":
            if err_sigma is not None:
                tol_s = o.get("tol_sigma") or 3.0
                if abs(err_sigma) > tol_s:
                    bad_impl.append(o)
            else:
                tol = o.get("tol_pct") or 5.0
                if err is None or abs(err) > tol:
                    bad_impl.append(o)
        elif role == "hypothesis":
            if err is not None and abs(err) > (o.get("tol_pct") or 10.0):
                dev_hypo.append(o)
        else:
            meas.append(o)
    # rule 7: an item that implementation-checks the full combination carries a
    # different weight in the verdict -- flag it
    comp_impl = [o for o in observables
                 if o.get("role") == "implementation_check"
                 and o.get("scope", "composite") == "composite"]
    if bad_impl:
        v = f"FAIL -- {len(bad_impl)} implementation check(s) mismatched (a bug)"
    elif dev_hypo:
        v = f"PASS (implementation sound) . {len(dev_hypo)} hypothesis deviation(s) <- results"
    else:
        v = "PASS"
    if comp_impl:
        v += f"  [{len(comp_impl)} composite implementation check(s) -- rule 7: verify the derivation]"
    return v, bad_impl, dev_hypo, meas


def build(*, run_id, case, system_tags, reference_scales, physical, dimensionless,
          checks, observables, numerics, equilibration, wall_seconds,
          steps_per_second=None, extra=None) -> dict:
    """`checks` is a list of `(Check, phase)` tuples. phase = design | post_run."""
    m = {
        "schema": SCHEMA,
        "run_id": run_id,
        "case": case,
        "system_tags": list(system_tags),
        "reference_scales": dict(reference_scales),
        "physical": dict(physical),
        "dimensionless": {k: float(v) for k, v in dimensionless.items()},
        "checks": [c.as_dict(phase) for c, phase in checks],
        "observables": list(observables),
        "equilibration": dict(equilibration),
        "numerics": dict(numerics),
        "wall_seconds": float(wall_seconds),
    }
    if steps_per_second is not None:
        m["steps_per_second"] = float(steps_per_second)
    if extra:
        m.update(extra)
    return m


def write(outdir: Path, metrics: dict) -> Path:
    p = Path(outdir) / "metrics.json"
    p.write_text(json.dumps(metrics, indent=2))
    return p


def equilibration_series(series_key: str, label: str, source: str = "observables.npz") -> dict:
    """Declare the equilibrium indicator. The case names the time series that suits
    its own system.

    Bound system (a trap): displacement from the anchor.
    Diffusive or structural system: potential energy <U>/N.
    """
    return {"source": source, "series_key": series_key, "label": label}


__all__ = ["SCHEMA", "ROLES", "build", "write", "observable", "judge",
           "equilibration_series"]
