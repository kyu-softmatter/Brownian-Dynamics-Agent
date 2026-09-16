"""Does the sealed analysis plan name anything nobody produces?

    ./bin/py verify/verify_plan_is_implemented.py

★ Why this exists. Reading `campaigns/sediment_preregistration/analysis_plan.yaml`
against the code that is supposed to satisfy it has found **two real gaps**, both
of which would have surfaced only at the end of a 1.7-hour campaign:

| what the plan promised | what existed |
|---|---|
| step 3: the first half of the production window compared against the second | nothing. `Build.gsd_path` is declared in `bdbot/run.py` and never read, so no trajectory is written and the split cannot be recovered afterwards |
| figure F5: every sealed prediction against its measurement | `campaigns/s31_analyze.py` produced F1-F4, F6 and F7 |

Both were found by eye. A pre-registered plan that names an artefact nobody
produces is decoration, and the seal makes it *look* rigorous — which is the
defect class this repository keeps recording: a check whose success is
indistinguishable from its absence.

⚠ **This gate is about EXISTENCE, not correctness.** It asserts that every
figure id in the plan is emitted by the analyzer and that every observable the
plan's computations name is emitted by the case. It cannot tell whether the
figure shows the right thing. The half-window bug is the example: after the
accumulation existed, the statistic still split 6 frames against 30, and only a
balance assertion inside the run caught that.
"""
from __future__ import annotations

import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
PLAN = ROOT / "campaigns" / "sediment_preregistration" / "analysis_plan.yaml"
PRED = ROOT / "campaigns" / "sediment_preregistration" / "prediction.yaml"
ANALYZER = ROOT / "campaigns" / "s31_analyze.py"
CASE = ROOT / "cases" / "sediment_3d.py"


def figures_promised() -> list[str]:
    plan = yaml.safe_load(PLAN.read_text())
    return [f["id"] for f in plan["figures"]]


def figures_emitted() -> set[str]:
    """Figure ids the analyzer actually writes, by the filenames it constructs.

    ⚠ By PARSE of the string literals, not by running it -- but the strings are
    the filenames, so a promised id with no file is caught. `docs/05-pitfalls.md`
    records a grep that matched the prose explaining a rule, which is why this
    looks only at `FIGDIR / "..."` expressions.
    """
    src = ANALYZER.read_text()
    names = re.findall(r'FIGDIR\s*/\s*"([^"]+)"', src)
    out = set()
    for n in names:
        for tok in re.findall(r"F\d+", n):
            out.add(tok)
    return out


def observables_promised() -> set[str]:
    """Observables the two sealed documents name as things to read."""
    want = set()
    pred = yaml.safe_load(PRED.read_text())
    for p in pred["predictions"]:
        want.add(p["quantity"])
    want.add(pred["primary_statistic"]["name"])
    want.add(pred["secondary_statistic"]["name"])
    #  and anything the plan's `computations` block names as a key
    plan = yaml.safe_load(PLAN.read_text())
    for key in plan.get("computations", {}):
        want.add(key)
    return want


def observables_emitted() -> set[str]:
    """Observable names the case constructs. By parse of `MET.observable("...")`."""
    src = CASE.read_text()
    return set(re.findall(r'MET\.observable\(\s*\n?\s*"([A-Za-z0-9_]+)"', src))


#: Names the plan's `computations` block uses as headings for a QUANTITY rather
#: than for an observable. Listed explicitly so the gate does not pass by having
#: a loose matcher -- each one says why it is not an observable.
NOT_OBSERVABLES = {
    "profile": "the histogram itself, carried in observables.npz as arrays",
    "l_g_fitted": "the heading for the two windowed observables below it",
    "Z": "the estimator, not an observable -- `eos_Z` is the array",
    "phi_eff": "an axis, carried as the `eos_phi_eff` array",
    "Z_dev_wmean": "the heading for the two mapped versions",
}


def main() -> int:
    bad = []

    want_f, got_f = figures_promised(), figures_emitted()
    missing_f = [f for f in want_f if f not in got_f]
    print(f"figures: plan promises {want_f}")
    print(f"         analyzer emits {sorted(got_f)}")
    for f in missing_f:
        bad.append(f"figure {f} is named in the sealed plan and the analyzer "
                   f"emits no file for it")

    want_o, got_o = observables_promised(), observables_emitted()
    print(f"\nobservables: {len(want_o)} named across the two sealed documents, "
          f"{len(got_o)} emitted by the case")
    for name in sorted(want_o):
        if name in got_o or name in NOT_OBSERVABLES:
            continue
        bad.append(f"observable {name!r} is named in a sealed document and the "
                   f"case emits nothing by that name")
    for name in sorted(NOT_OBSERVABLES):
        if name in got_o:
            bad.append(f"{name!r} is on the NOT_OBSERVABLES list but the case "
                       f"now emits it -- remove it from the list")

    #  ★ and the gate on the gate: the exemption list must not be able to grow
    #    silently into a way of passing. Every entry carries a reason and every
    #    entry must still be named by the plan.
    plan_keys = set(yaml.safe_load(PLAN.read_text()).get("computations", {}))
    for name, why in NOT_OBSERVABLES.items():
        if not why.strip():
            bad.append(f"NOT_OBSERVABLES[{name!r}] has no reason")
        if name not in plan_keys:
            bad.append(f"NOT_OBSERVABLES names {name!r}, which the plan's "
                       f"computations block no longer mentions -- a stale "
                       f"exemption is a hole")

    print("\n" + "=" * 70)
    for b in bad:
        print(f"  x {b}")
    print(f"\n{len(bad)} promise(s) with nothing behind them" if bad else
          "\nevery figure and every observable the sealed documents name is "
          "produced")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
