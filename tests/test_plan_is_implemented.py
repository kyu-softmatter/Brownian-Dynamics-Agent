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

import re
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


def test_the_analyzer_derives_its_gates_from_the_declared_roles():
    """★ The analyzer must not keep its own copy of which observables gate.

    Measured: revision 4 demoted `Z_dilute_tail` from implementation_check to
    measurement -- a rule 7' correction, because Z -> 1 in the dilute tail holds
    only AT equilibrium and equilibrium is what the campaign measures. The
    analyzer's hardcoded `GATES` tuple went on gating it and reported
    IMPLEMENTATION FAILURE on a run behaving exactly as the physics requires.
    A second copy of a role declaration is a second thing to forget.
    """
    src = (ROOT / "campaigns" / "s31_analyze.py").read_text()
    assert "GATES = (" not in src, "a hardcoded gate list is back"
    assert 'ob.get("role") != "implementation_check"' in src
    #  and a run with no implementation_check at all must not silently pass
    assert "there is nothing to gate on, which is not" in src


def test_the_case_and_the_sealed_prediction_agree_on_every_role():
    """A role in the sealed document that the case does not implement is a
    promise with nothing behind it -- the same defect as a missing figure, one
    level down."""
    pred = VP.yaml.safe_load(VP.PRED.read_text())
    sealed = {p["quantity"]: p["role"] for p in pred["predictions"]}
    src = (ROOT / "cases" / "sediment_3d.py").read_text()
    #  ⚠ Split on `MET.observable(` FIRST. The first version searched for the
    #  name and then scanned forward, which found the name inside a comment or a
    #  neighbouring observable's note and read that block's role -- three false
    #  hits. `docs/05-pitfalls.md` records the same shape: a string search that
    #  matched the prose explaining the thing it was checking.
    blocks = src.split("MET.observable(")[1:]
    roles = {}
    for b in blocks:
        m = re.match(r'\s*\n?\s*"([A-Za-z0-9_]+)"', b)
        if not m:
            continue
        mr = re.search(r'role="([a-z_]+)"', b)
        roles[m.group(1)] = mr.group(1) if mr else "measurement"   # the default
    bad = []
    for name, role in sealed.items():
        if name not in roles:
            bad.append(f"{name}: the sealed document names it and the case emits "
                       f"no observable by that name")
        elif roles[name] != role:
            bad.append(f"{name}: the sealed document says {role!r}, the case "
                       f"sets {roles[name]!r}")
    assert not bad, bad


def test_a_sealed_copy_and_the_source_never_share_a_revision_number_while_differing():
    """★★ The failure this pins actually happened, mid-campaign.

    Two vestigial `role:` fields were corrected in the SOURCE
    `campaigns/sediment_preregistration/prediction.yaml` while eight run
    directories were already sealed against it. The seals stayed valid — they
    cover the copies — but both files then called themselves *revision 4* with
    different content, which makes "what was pre-registered?" unanswerable from
    the revision number alone.

    So: for every sealed copy, either it is byte-identical to the source, or its
    `revision` differs. Never the same number and different bytes.
    """
    import hashlib

    src_path = VP.PRED
    src_bytes = src_path.read_bytes()
    src_rev = VP.yaml.safe_load(src_bytes)["revision"]
    sealed = sorted((ROOT / "runs").glob("*/prediction.yaml"))
    assert sealed, "no sealed copy exists -- this check would be vacuous"

    bad = []
    for q in sealed:
        b = q.read_bytes()
        if b == src_bytes:
            continue
        rev = VP.yaml.safe_load(b).get("revision")
        if rev == src_rev:
            bad.append(
                f"{q.parent.name}: sealed copy and source both say revision "
                f"{rev!r} but differ ("
                f"{hashlib.sha256(b).hexdigest()[:12]} vs "
                f"{hashlib.sha256(src_bytes).hexdigest()[:12]})")
    assert not bad, "\n".join(bad)


def test_the_divergence_from_the_running_campaign_is_declared_at_the_top():
    """A source that has moved past the sealed copies must SAY so where a reader
    will see it, not only in a revision-history entry near the bottom."""
    head = VP.PRED.read_text()[:3000]
    sealed = sorted((ROOT / "runs").glob("*/prediction.yaml"))
    if not sealed:
        return
    src_rev = VP.yaml.safe_load(VP.PRED.read_text())["revision"]
    revs = {VP.yaml.safe_load(q.read_bytes()).get("revision") for q in sealed}
    if revs == {src_rev}:
        return                     # nothing has diverged; nothing to declare
    assert "SEALED AGAINST REVISION" in head.upper(), (
        f"the source is revision {src_rev} and the sealed copies are {revs}, and "
        f"the top of the file does not say so")
    #  and it must name what differs, rather than just that something does
    assert "role:" in head and "measurement" in head
