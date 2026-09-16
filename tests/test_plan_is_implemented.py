"""The sealed analysis plan must not name anything nobody produces.

★ `verify/verify_plan_is_implemented.py` found **three** gaps in this campaign,
and the third one it found by itself:

| what the plan promised | found by |
|---|---|
| step 3's first-half/second-half profile comparison | reading the plan against the code, 12 min into the first production run |
| figure F5, every sealed prediction against its measurement | the same, by eye |
| figure F4, the two cross-sections | **this gate, on its first run** |

A pre-registered plan that names an artefact nobody produces is decoration, and
the seal makes it look rigorous — which is this repository's most-repeated defect
class. `docs/06-roadmap.md` records three prior instances of a practice written
down and enforced nowhere, so the review is a gate rather than a habit.
"""
from __future__ import annotations

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify_plan_is_implemented as VP  # noqa: E402


def test_every_promised_figure_is_emitted():
    want, got = VP.figures_promised(), VP.figures_emitted()
    assert want, "the plan names no figures -- the check would be vacuous"
    missing = [f for f in want if f not in got]
    assert not missing, (
        f"the sealed plan promises {missing} and the analyzer emits no file for "
        f"them; it emits {sorted(got)}")


def test_every_promised_observable_is_emitted():
    want, got = VP.observables_promised(), VP.observables_emitted()
    assert len(want) >= 10 and len(got) >= 10, (len(want), len(got))
    missing = [n for n in sorted(want)
               if n not in got and n not in VP.NOT_OBSERVABLES]
    assert not missing, missing


def test_the_exemption_list_cannot_grow_into_a_hole():
    """★ A gate with an exemption list is only as good as the list's discipline.
    Every entry must carry a reason AND still be named by the plan, so an
    exemption cannot outlive the thing it exempts."""
    plan_keys = set(VP.yaml.safe_load(VP.PLAN.read_text()).get("computations", {}))
    for name, why in VP.NOT_OBSERVABLES.items():
        assert why.strip(), f"{name!r} is exempted with no reason"
        assert name in plan_keys, (
            f"{name!r} is exempted but the plan's computations block no longer "
            f"names it -- a stale exemption is a hole")
        assert name not in VP.observables_emitted(), (
            f"{name!r} is exempted as 'not an observable' but the case now emits "
            f"it; remove the exemption")


def test_the_gate_catches_a_missing_figure():
    """Mutation, in memory: drop a figure from what the analyzer emits."""
    real = VP.figures_emitted
    try:
        VP.figures_emitted = lambda: real() - {"F7"}
        want, got = VP.figures_promised(), VP.figures_emitted()
        assert [f for f in want if f not in got] == ["F7"]
    finally:
        VP.figures_emitted = real
    assert "F7" in VP.figures_emitted()


def test_the_script_exits_nonzero_when_it_finds_a_gap():
    """The gate has to FAIL, not warn. Run it as CI would."""
    r = subprocess.run([sys.executable, str(VP.__file__)], cwd=ROOT,
                       capture_output=True, text=True)
    assert r.returncode == 0, r.stdout[-2000:]
    assert "every figure and every observable" in r.stdout
