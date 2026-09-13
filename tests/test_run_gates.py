"""L4 — the pre-run gates inside `run.execute()`, and whether each can fire.

**Why this file exists.** `verify/verify_gates_bite.py` deletes each enforcement
point in the repository in turn and counts how many tests fail. Measured
2026-09-13, five gates in `bdbot/run.py` had **zero** tests fail:

    execute() refuses a hand-edited spec (rule 2)                 0
    execute() refuses a spec whose L3 verdict is FAIL             0
    rule 10: require_approval=True with no params.json            0
    rule 10: an unapproved params.json does not run               0
    rule 10: approved one set of numbers, running another         0

★ Three of those are rule 10's entire enforcement, and rule 10 was introduced
  *because* "this repository's record on practices-without-gates is three for
  three against." It was four for four.

**And two of the five were not merely untested — they were unreachable.**
`RID.prepare_outdir` (run.py:399) deletes every file in an existing output
directory except those in `runid.PRESERVE`, and it runs BEFORE the seal check
(416) and the `params.json` branch (424). Measured on a temp directory before
the fix:

    a VALID seal        -> "no SEALED.sha256. Refusing to run"
    a TAMPERED seal     -> "[warn] ... is unsealed", and the run PROCEEDED
    an APPROVED params  -> "no params.json. Rule 10 requires every number"
    a DRIFTING params   -> the drift check never ran

So `require_seal=True` and `require_approval=True` could only ever fail, never
be satisfied. The fix is on `runid.PRESERVE`, which now keeps the human-authored
inputs and whatever the seal covers.

`TestPreserve` below is the regression test for that root cause, and it is the
one to keep if any other is dropped: with `PRESERVE` reverted to
`{"record.json"}`, every gate test in the two rule-10 classes still passes or
fails for the *wrong reason*, because the file under test is gone before the
gate looks.

**Each gate is tested twice, deliberately.** Once that it FIRES, with a
`build_fn` that raises a sentinel so the refusal has to arrive before the build;
and once that it can be SATISFIED, reaching `build_fn`. A gate that refuses
everything is worse than no gate — this repository's own pre-run gate once
rejected 80 of 83 specs with zero real failures among them.
"""
from __future__ import annotations

import hashlib
import json
import pathlib

import pytest

from bdbot import nondim as ND, params as PRM, run as RUN, runid as RID
from bdbot import runcard as RC

ROOT = pathlib.Path(__file__).resolve().parent.parent
SENTINEL = "build_fn was reached"


def _spec_json() -> dict:
    """An archived spec, hash-valid and verdict PASS. Archived specs carry no
    `system.structure`, so the three structure gates above these ones are
    skipped and cannot mask what is under test."""
    src = sorted((ROOT / "specs").glob("trap-2d-5um__*.json"))[0]
    return json.loads(src.read_text())


@pytest.fixture
def spec_factory(tmp_path):
    def make(mutate=None, name="spec.json"):
        d = _spec_json()
        if mutate is not None:
            mutate(d)
        p = tmp_path / name
        p.write_text(json.dumps(d))
        return ND.load(p)
    return make


def _build_fn(_spec, _outdir):
    raise AssertionError(SENTINEL)


def _run(spec, outdir, **kw):
    return RUN.execute(spec, _build_fn, outdir, progress=False, **kw)


def _reaches_build(spec, outdir, **kw):
    """The gate let it through. AssertionError is not ValueError, so a gate
    firing cannot be mistaken for this."""
    with pytest.raises(AssertionError, match=SENTINEL):
        _run(spec, outdir, **kw)
    return True


# ── the root cause ─────────────────────────────────────────────────────────

class TestPreserve:
    """`prepare_outdir` clears the OUTPUTS of a partial run. It must not delete
    the human-authored INPUTS that the gates downstream read."""

    def test_the_gate_inputs_survive_prepare_outdir(self, tmp_path):
        d = tmp_path / "run"
        d.mkdir()
        for n in ("params.json", "SEALED.sha256", "record.json"):
            (d / n).write_text("{}" if n.endswith("json") else "")
        (d / "metrics.json").write_text("{}")       # a real output
        (d / "traj.gsd").write_bytes(b"")           # a real output
        go, _ = RID.prepare_outdir(d, force=False)
        assert go
        left = {p.name for p in d.iterdir()}
        assert left == {"params.json", "SEALED.sha256", "record.json"}, left

    def test_the_documents_a_seal_covers_also_survive(self, tmp_path):
        """Keeping `SEALED.sha256` while deleting what it hashes leaves a
        dangling seal, which `verify_seal` then reports as a missing document —
        trading a silent failure for a loud one that blocks every re-run."""
        d = tmp_path / "run"
        d.mkdir()
        doc = d / "02_prediction.md"
        doc.write_text("K' = 0\n")
        (d / "SEALED.sha256").write_text(
            f"{hashlib.sha256(doc.read_bytes()).hexdigest()}  "
            f"runs/whatever/02_prediction.md\n")
        (d / "metrics.json").write_text("{}")
        RID.prepare_outdir(d, force=False)
        assert doc.exists(), "the sealed document was deleted; the seal dangles"
        assert not (d / "metrics.json").exists(), "an output was not cleared"

    def test_outputs_are_still_cleared(self, tmp_path):
        """...and PRESERVE has not grown into 'delete nothing'. Without this,
        the test above passes for a `prepare_outdir` that does no work at all."""
        d = tmp_path / "run"
        d.mkdir()
        for n in ("metrics.json", "traj.gsd", "result.partial", "log.txt"):
            (d / n).write_text("x")
        RID.prepare_outdir(d, force=False)
        assert [p.name for p in d.iterdir()] == [], "outputs were not cleared"

    def test_preserve_is_what_this_file_thinks_it_is(self):
        """Escapes self-reference: the tests above are written against these
        three names, so the set itself must be pinned rather than read."""
        assert RID.PRESERVE == {"record.json", "params.json", "SEALED.sha256"}


# ── rule 2 · a hand-edited spec ────────────────────────────────────────────

class TestSpecHash:
    def test_a_hand_edited_spec_is_refused(self, spec_factory, tmp_path):
        spec = spec_factory(lambda d: d["system"].__setitem__("dimensions", 3))
        with pytest.raises(ValueError, match="hand-edited"):
            _run(spec, tmp_path / "r")

    def test_an_untouched_spec_passes_the_hash_gate(self, spec_factory, tmp_path):
        """Non-vacuity: the same call path with the document unmodified must get
        past this gate. Otherwise the test above would pass for a gate that
        refuses every spec."""
        assert _reaches_build(spec_factory(), tmp_path / "r")

    def test_editing_a_doc_key_does_not_trip_it(self, spec_factory, tmp_path):
        """`runid.DOC_KEYS` are excluded from the hash on purpose — writing down
        *why* a value is what it is must not rename the run."""
        spec = spec_factory(lambda d: d["system"].__setitem__(
            "description", "edited prose, not physics"))
        assert _reaches_build(spec, tmp_path / "r")


# ── the L3 verdict ─────────────────────────────────────────────────────────

class TestVerdict:
    def test_a_failed_spec_is_not_run(self, spec_factory, tmp_path):
        """`verdict` is outside the hash payload, so this reaches the verdict
        gate rather than tripping the hash gate one line above it."""
        spec = spec_factory(lambda d: d.__setitem__("verdict", "FAIL (2 hard)"))
        with pytest.raises(ValueError, match="L3 verdict is FAIL"):
            _run(spec, tmp_path / "r")

    def test_setting_the_verdict_did_not_break_the_hash(self, spec_factory):
        """The premise of the test above, asserted rather than assumed. If
        `verdict` ever enters the hash, that test would still pass — on the
        wrong gate, with the wrong message."""
        spec = spec_factory(lambda d: d.__setitem__("verdict", "FAIL (2 hard)"))
        ok, _ = spec.verify_hash()
        assert ok, "verdict is now hashed; the verdict test is testing the hash"

    def test_a_passing_spec_runs(self, spec_factory, tmp_path):
        assert _reaches_build(spec_factory(), tmp_path / "r")

    @pytest.mark.parametrize("verdict", ["PASS", "PASS (3 warnings)"])
    def test_warnings_do_not_block(self, spec_factory, tmp_path, verdict):
        """`health.gate()` once tested `verdict != "PASS"` and rejected
        `"PASS (3 warnings)"` — 80 false rejections out of 83 specs. `execute()`
        reads `startswith("FAIL")`, which is the correct half; pin it."""
        spec = spec_factory(lambda d: d.__setitem__("verdict", verdict))
        assert _reaches_build(spec, tmp_path / "r")


# ── rule 10 · the approved numbers ─────────────────────────────────────────

def _manifest(dt=3.3135576e-08, approved=None) -> PRM.Manifest:
    m = PRM.Manifest(case="trap-2d-5um")
    m.add("numerics", "dt", dt, "1", "the spec", tier=0)
    m.approved_by = approved
    return m


class TestRule10:
    def test_no_params_json_blocks(self, spec_factory, tmp_path):
        with pytest.raises(ValueError, match="no params.json"):
            _run(spec_factory(), tmp_path / "r", require_approval=True)

    def test_an_unapproved_manifest_blocks(self, spec_factory, tmp_path):
        d = tmp_path / "r"
        d.mkdir()
        _manifest().write(d)
        with pytest.raises(ValueError, match="not approved"):
            _run(spec_factory(), d, require_approval=True)

    def test_an_approved_matching_manifest_runs(self, spec_factory, tmp_path):
        """★ The one that was impossible before `PRESERVE` was fixed: the gate
        has to be satisfiable, not only failable. Until 2026-09-13 this raised
        "no params.json" because `prepare_outdir` had deleted it."""
        d = tmp_path / "r"
        d.mkdir()
        _manifest(approved="a-human").write(d)
        assert _reaches_build(spec_factory(), d, require_approval=True)

    def test_drift_blocks(self, spec_factory, tmp_path):
        d = tmp_path / "r"
        d.mkdir()
        _manifest(dt=9.9e-09, approved="a-human").write(d)
        with pytest.raises(ValueError, match="does not match the spec"):
            _run(spec_factory(), d, require_approval=True)

    def test_drift_blocks_even_without_require_approval(self, spec_factory,
                                                        tmp_path):
        """★ CLAUDE.md rule 10: *"if a `params.json` exists it is **always**
        checked against the spec's `numerics` — approving one set of numbers and
        running another is worse than never writing them down."* That sentence
        was false: the check sat behind a branch the file never survived to
        reach."""
        d = tmp_path / "r"
        d.mkdir()
        _manifest(dt=9.9e-09, approved="a-human").write(d)
        with pytest.raises(ValueError, match="does not match the spec"):
            _run(spec_factory(), d, require_approval=False)

    def test_a_matching_manifest_does_not_trip_the_drift_check(
            self, spec_factory, tmp_path):
        """...and the drift check compares, rather than always firing."""
        d = tmp_path / "r"
        d.mkdir()
        _manifest(approved="a-human").write(d)
        assert _reaches_build(spec_factory(), d, require_approval=False)

    def test_the_drift_check_reads_the_value_this_file_uses(self, spec_factory):
        """Pins the premise both drift tests stand on: the manifest's `dt` is
        matched against the spec's `dt_star`, and the default above really is
        that value. If the spec changes, these tests must not quietly become
        no-ops that pass because nothing was compared."""
        spec = spec_factory()
        assert PRM.diff_against_spec(_manifest(), dict(spec.numerics)) == []
        drift = PRM.diff_against_spec(_manifest(dt=9.9e-09),
                                      dict(spec.numerics))
        assert drift and "dt" in drift[0], drift


# ── the seal ───────────────────────────────────────────────────────────────

def _seal(d: pathlib.Path, body: str = "K' = 0\n", sealed: str | None = None):
    doc = d / "prediction.md"
    doc.write_text(sealed if sealed is not None else body)
    digest = hashlib.sha256(body.encode()).hexdigest()
    (d / RC.SEAL_NAME).write_text(f"{digest}  prediction.md\n")


class TestSeal:
    def test_require_seal_blocks_an_unsealed_run(self, spec_factory, tmp_path):
        with pytest.raises(RC.SealBroken, match="Refusing to run"):
            _run(spec_factory(), tmp_path / "r", require_seal=True)

    def test_a_valid_seal_satisfies_it(self, spec_factory, tmp_path):
        """★ Also impossible before the `PRESERVE` fix."""
        d = tmp_path / "r"
        d.mkdir()
        _seal(d)
        assert _reaches_build(spec_factory(), d, require_seal=True)

    def test_a_tampered_seal_is_a_hard_stop_even_unrequired(self, spec_factory,
                                                            tmp_path):
        """★ `runcard.verify_or_raise`: *"unsealed is a known historical state,
        broken is tampering."* Before the fix the seal file was deleted first,
        so tampering was downgraded to "[warn] ... is unsealed" and the run
        proceeded. `require_seal` is deliberately left False here — that is the
        whole asymmetry."""
        d = tmp_path / "r"
        d.mkdir()
        _seal(d, sealed="K' = 999 TAMPERED\n")
        with pytest.raises(RC.SealBroken, match="seal broken"):
            _run(spec_factory(), d)

    def test_sealbroken_is_not_a_valueerror(self):
        """Every other gate here raises ValueError. If `SealBroken` ever became
        one, the two tests above would still pass while no longer distinguishing
        a seal failure from any other refusal."""
        assert issubclass(RC.SealBroken, RuntimeError)
        assert not issubclass(RC.SealBroken, ValueError)
