"""L2 — `structure.dim`, the gate that gives `D9` teeth.

`knowledge/wiki/concepts/dimensionality-has-no-default.md` closes
`D9 · 차원 (2D / 3D)`, which had stood `OPEN` with a **3D default that 7 of the 8
cases contradict**, none of them recording why. The decision: no default; every
case declares a `basis`.

**Why this file exists rather than only `verify/verify_dim_gate.py`.** CI runs
`pytest -q -rs` and nothing else (`.github/workflows/ci.yml:127`), so a verify
script is not enforcement — it is a claim someone has to remember to run. That
distinction is this repository's most-repeated failure: `health.gate()` is
reachable only from `tools/health.py` and `run.py` never mentions it, so **no run
has ever gated itself**. `test_s1_intake_goal.py` was added for the same reason
and says so in its own docstring.

Measured before this file existed (2026-09-10): commenting out the single line
`out += check_dim(s)` in `physical.validate` left **1185 passed, 2 skipped** —
the gate could be silently unwired and the suite would not notice. `test_wiring`
below is the test that closes that, and it is the one to keep if any other is
ever dropped.
"""
from __future__ import annotations

import copy
import hashlib
import json
import pathlib

import pytest
import yaml

from bdbot import physical as P
from bdbot import runid as RID

ROOT = pathlib.Path(__file__).resolve().parent.parent
DONOR = ROOT / "intake/soft-r3-2d-A-sweep"
CASES = sorted(p.parent for p in ROOT.glob("intake/*/system.yaml"))


@pytest.fixture(scope="module")
def donor_sys() -> dict:
    return yaml.safe_load((DONOR / "system.yaml").read_text())


@pytest.fixture
def case_dir(tmp_path, donor_sys):
    """A writable copy of the donor case. Returns a writer that installs a
    mutated `system.yaml` and hands back the loaded `PhysicalSystem`.
    """
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    (tmp_path / "observation.yaml").write_text((DONOR / "observation.yaml").read_text())

    def write(mutate=None):
        d = copy.deepcopy(donor_sys)
        if mutate is not None:
            mutate(d)
        (tmp_path / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
        return P.load(tmp_path)

    return write


def _dim_errors(s) -> list:
    return [i for i in s.errors if i.where.startswith("structure.dim")]


# ── the one that matters ───────────────────────────────────────────────────

def test_wiring(case_dir):
    """★ The gate is reached through `physical.load` -> `validate`, not only by
    calling `check_dim` directly.

    This is the test that fails if someone removes `out += check_dim(s)`. Every
    other test in this file would still pass with the gate unwired, because they
    could all be satisfied by `check_dim` in isolation.
    """
    s = case_dir(lambda d: d.pop("structure"))
    assert _dim_errors(s), "structure.dim missing produced no error via load()"


def test_the_donor_passes_unmutated(case_dir):
    """Guard on the guard. Every break test below starts from this document, so
    if it does not pass clean, a 'caught' result proves nothing.
    """
    s = case_dir()
    assert not _dim_errors(s), [str(i) for i in _dim_errors(s)]


# ── each rule fires ────────────────────────────────────────────────────────

@pytest.mark.parametrize("mutate,where", [
    (lambda d: d.pop("structure"), "structure.dim"),
    (lambda d: d["structure"].pop("dim"), "structure.dim"),
    (lambda d: d["structure"]["dim"].pop("basis"), "structure.dim.basis"),
    (lambda d: d["structure"]["dim"].update(basis="because"), "structure.dim.basis"),
    (lambda d: d["structure"]["dim"].pop("what_would_change"),
     "structure.dim.what_would_change"),
    (lambda d: d["structure"]["dim"].update(what_would_change="   "),
     "structure.dim.what_would_change"),
    (lambda d: d["structure"]["dim"].pop("alternatives"), "structure.dim.alternatives"),
    (lambda d: d["structure"]["dim"].update(alternatives=[]), "structure.dim.alternatives"),
    (lambda d: d["structure"]["dim"].pop("value"), "structure.dim.value"),
], ids=["no-section", "no-dim", "no-basis", "bad-basis", "no-wwc", "blank-wwc",
        "no-alternatives", "empty-alternatives", "no-value"])
def test_a_broken_field_is_caught(case_dir, mutate, where):
    s = case_dir(mutate)
    assert any(i.where == where for i in s.errors), \
        f"{where} not caught: {[str(i) for i in s.errors]}"


def test_blank_what_would_change_says_why_it_matters(case_dir):
    """"Not stated" and "no consequence" must not render the same, and the
    message is what carries that. A bare "missing" would let a reader fill it
    with "nothing" and move on.
    """
    s = case_dir(lambda d: d["structure"]["dim"].update(what_would_change=" "))
    msg = next(i.msg for i in s.errors if i.where == "structure.dim.what_would_change")
    assert "not the same as it being safe" in msg


# ── the cross-check: one choice, one number ────────────────────────────────

def test_value_must_agree_with_the_hashed_dimensions(case_dir):
    """`run_id` is hashed on top-level `dimensions`; `structure.dim.value` is
    not. If they drift, the file records one choice and runs another — the same
    invariant `params.json` carries against the spec's numerics.
    """
    s = case_dir(lambda d: d["structure"]["dim"].update(value=3))
    assert any(i.where == "structure.dim.value" for i in s.errors)


def test_agreement_is_not_asserted_vacuously(case_dir):
    """...and the check passes when they do agree, so the test above is testing
    the comparison rather than the presence of the key."""
    s = case_dir(lambda d: (d["structure"]["dim"].update(value=3),
                            d.update(dimensions=3)))
    assert not [i for i in s.errors if i.where == "structure.dim.value"]


# ── each basis owes different evidence ─────────────────────────────────────

@pytest.mark.parametrize("basis,owed", sorted(
    (b, f) for b, fs in P.DIM_BASIS_REQUIRES.items() for f in fs))
def test_each_basis_requires_its_own_evidence(case_dir, basis, owed):
    """`given` owes a source, `inherited` a named comparison, `sufficient` the
    observable list it was checked against. Without the field, `inherited` is
    indistinguishable from inertia and `sufficient` from a guess.
    """
    def mutate(d):
        e = d["structure"]["dim"]
        e["basis"] = basis
        e.pop(owed, None)
    s = case_dir(mutate)
    assert any(i.where == f"structure.dim.{owed}" for i in s.errors), \
        f"basis={basis} without {owed} was accepted"


def test_required_owes_nothing_extra(case_dir):
    """A basis with no extra obligation must not be blocked by one. A gate that
    refuses everything is worse than no gate — this repository's own pre-run
    gate once rejected 80 of 83 specs with zero real failures.
    """
    s = case_dir(lambda d: d["structure"]["dim"].update(basis="required"))
    assert not _dim_errors(s), [str(i) for i in _dim_errors(s)]


@pytest.mark.parametrize("basis", ["inherited", "sufficient"])
def test_the_expiring_bases_announce_themselves(case_dir, basis):
    """An expired basis looks exactly like a valid one, so the warning is the
    only thing that makes the expiry visible. It is a warn, not an error:
    whether the comparison still holds is not knowable from this file.
    """
    def mutate(d):
        e = d["structure"]["dim"]
        e["basis"] = basis
        e["compared_with"] = ["some-other-case"]
        e["checked_observables"] = ["x2"]
    s = mutate and case_dir(mutate)
    warns = [i for i in s.issues
             if i.level == "warn" and i.where.startswith("structure.dim")]
    assert any("EXPIRES" in i.msg for i in warns), [str(i) for i in s.issues]
    assert not _dim_errors(s), "an expiring basis must warn, not block"


# ── the real tree ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_every_real_case_declares_a_basis(case):
    """The eight worked cases were not written to satisfy this checker — they
    are the strongest available evidence that it lets a correct answer through.
    """
    s = P.load(case)
    assert not _dim_errors(s), [str(i) for i in _dim_errors(s)]
    basis = ((s.raw.get("structure") or {}).get("dim") or {}).get("basis")
    assert basis in P.DIM_BASES, f"{case.name}: basis={basis!r}"


def test_the_case_count_matches_the_parametrisation():
    """If a case appears, the parametrisation above must have picked it up.

    Was `== 8` and named `..._eight_cases` until 2026-09-16, when
    `sediment-pmma-3d` made it nine. Renamed rather than re-numbered in place: a
    test called `eight` that asserts nine is a lie in the one place a reader
    looks first.
    """
    assert len(CASES) == 9, [c.name for c in CASES]


# ── run_id must not move ───────────────────────────────────────────────────

def _strip(node, keys):
    if isinstance(node, dict):
        return {k: _strip(v, keys) for k, v in node.items() if k not in keys}
    if isinstance(node, list):
        return [_strip(v, keys) for v in node]
    return node


def _h(d) -> str:
    return hashlib.sha256(
        json.dumps(d, sort_keys=True, default=str).encode()).hexdigest()[:12]


def test_structure_is_excluded_from_the_run_id():
    """Writing down *why* the dimension is 2 must not rename a run."""
    a = {"dimensions": 2}
    b = {"dimensions": 2, "structure": {"dim": {"value": 2, "basis": "given"}}}
    assert _h(_strip(a, RID.DOC_KEYS)) == _h(_strip(b, RID.DOC_KEYS))


def test_dimensions_is_still_hashed():
    """...and the reverse direction, which `runid.py` warns about just as loudly:
    changing the physics must still re-id the run.
    """
    assert _h(_strip({"dimensions": 2}, RID.DOC_KEYS)) \
        != _h(_strip({"dimensions": 3}, RID.DOC_KEYS))


def test_a_spec_carrying_the_block_hashes_the_same_as_one_without():
    """The invariant, stated on the object it is about.

    ⚠️ This started as *"no archived spec carries a `structure` key, so adding it
    to DOC_KEYS is a no-op"* — and that premise died within the hour. Re-running
    `cases/trap_2d_5um.py` regenerates its spec **from the system.yaml**, so the
    block appears in `specs/` the moment any case is re-run. The measured proof
    that the exclusion works is better than the premise was: that spec now
    carries the block and is **still named `trap-2d-5um__ed80885e01b3`**, the
    run_id it had before the block existed.
    """
    without = {"case": "x", "system": {"dimensions": 2, "label": "x"}}
    with_ = json.loads(json.dumps(without))
    with_["system"]["structure"] = {
        "dim": {"value": 2, "basis": "given", "source": "s",
                "alternatives": [3], "what_would_change": "w"}}
    assert _h(RID.physics_only(without)) == _h(RID.physics_only(with_))


def test_every_archived_spec_is_named_by_its_own_run_id():
    """256 run directories are named by the `run_id` their spec records. If a
    field ever silently enters or leaves the hash, the two drift apart — this is
    the general form of the check, and it survives specs being regenerated.
    """
    specs = sorted((ROOT / "specs").glob("*.json"))
    assert len(specs) > 250, f"only {len(specs)} specs found"
    bad = [f.name for f in specs
           if json.load(open(f)).get("run_id") != f.stem]
    assert not bad, bad[:5]


# ── `given` owes more than a non-empty string ──────────────────────────────
#
# ★ `basis: given` asserts "the input states it", and `DIM_BASIS_REQUIRES` used
#   to discharge that with the mere PRESENCE of a `source` field. A file and a
#   half-remembered conversation satisfied it identically.
#
#   Measured 2026-09-15, auditing which `structure` claims an artefact could
#   settle: of the three `given` cases, `network` quotes `observation.yaml` A4
#   verbatim and `trap-2d-5um` cites the sketch's own formula, while
#   `chain-bend-2d-oscill` rested entirely on an exchange dated 2026-09-10 that
#   appears nowhere in that case's intake -- its own `source` prose says "The
#   sketch itself is silent." The gate passed all three the same way. (The claim
#   was then put to the user and confirmed, so the `given` stands; what did not
#   stand was the check.)
#
#   ⚠ And the first attempt to tell them apart was a keyword scan of the source
#     string for "sketch"/"observation.yaml"/".jpeg". It labelled the one
#     conversational case "artefact", because its text contains the word
#     "sketch" while saying the sketch says nothing. Hence `source_kind`: the
#     author declares, the gate checks. Parse, do not grep.

@pytest.fixture
def given_case(tmp_path, donor_sys):
    """A `basis: given` document whose artefacts exist inside the case dir."""
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    (tmp_path / "observation.yaml").write_text(
        (DONOR / "observation.yaml").read_text())

    def write(mutate=None):
        d = copy.deepcopy(donor_sys)
        d["structure"]["dim"].update(basis="given", source="the sketch says so",
                                     source_kind="artefact",
                                     source_file="sketch_01.jpeg")
        if mutate is not None:
            mutate(d)
        (tmp_path / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
        return P.load(tmp_path)

    return write


def test_a_declared_artefact_that_exists_is_accepted(given_case):
    """Guard on the guard. Every break below starts here."""
    assert not _dim_errors(given_case()), [str(i) for i in _dim_errors(given_case())]


@pytest.mark.parametrize("mutate,where", [
    (lambda d: d["structure"]["dim"].pop("source_kind"), "structure.dim.source_kind"),
    (lambda d: d["structure"]["dim"].update(source_kind="vibes"),
     "structure.dim.source_kind"),
    (lambda d: d["structure"]["dim"].pop("source_file"), "structure.dim.source_file"),
    (lambda d: d["structure"]["dim"].update(source_file="no_such_file.jpeg"),
     "structure.dim.source_file"),
    (lambda d: (d["structure"]["dim"].update(source_kind="human"),
                d["structure"]["dim"].pop("source_file")),
     "structure.dim.confirmed_by"),
    (lambda d: (d["structure"]["dim"].update(source_kind="human", confirmed_by="  "),
                d["structure"]["dim"].pop("source_file")),
     "structure.dim.confirmed_by"),
], ids=["no-kind", "bad-kind", "artefact-no-file", "artefact-file-missing",
        "human-no-confirmation", "human-blank-confirmation"])
def test_each_given_obligation_fires(given_case, mutate, where):
    s = given_case(mutate)
    assert any(i.where == where for i in s.errors), \
        f"{where} not caught: {[str(i) for i in s.errors]}"


def test_a_missing_artefact_is_the_point_of_the_check(given_case):
    """The message has to say the path resolved to nothing, because "the field
    is present" and "the file is there" are exactly what this distinguishes."""
    s = given_case(lambda d: d["structure"]["dim"].update(source_file="ghost.jpeg"))
    msg = next(i.msg for i in s.errors if i.where == "structure.dim.source_file")
    assert "does not exist" in msg and "ghost.jpeg" in msg, msg


def test_a_human_source_is_accepted_when_it_is_confirmed(given_case):
    """A conversation is a legitimate source -- `chain-bend-2d-oscill`'s
    dimension really was specified by a person. It just has to say who and be
    confirmed, rather than being indistinguishable from a file."""
    s = given_case(lambda d: (
        d["structure"]["dim"].update(source_kind="human",
                                     confirmed_by="user, 2026-09-15"),
        d["structure"]["dim"].pop("source_file")))
    assert not _dim_errors(s), [str(i) for i in _dim_errors(s)]


def test_an_artefact_at_the_repo_root_resolves(given_case):
    """A paper distillation is a legitimate artefact and lives at the root, not
    in the case directory. Both roots are tried."""
    s = given_case(lambda d: d["structure"]["dim"].update(
        source_file="knowledge/source/papers/2024-quah-graybox-abp-mpc-repo.md"))
    assert not _dim_errors(s), [str(i) for i in _dim_errors(s)]


def test_the_other_bases_are_not_asked_for_a_source_kind(given_case):
    """`required` and `inherited` do not claim the input states anything, so the
    obligation must not leak onto them -- a gate that refuses everything is
    worse than none."""
    for basis, extra in (("required", {}),
                         ("inherited", {"compared_with": ["some-case"]})):
        s = given_case(lambda d, b=basis, x=extra: (
            d["structure"]["dim"].update(basis=b, **x),
            [d["structure"]["dim"].pop(k, None)
             for k in ("source_kind", "source_file")]))
        assert not [i for i in s.errors
                    if i.where.endswith(("source_kind", "source_file"))], \
            f"basis={basis} was asked for a source_kind"


@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_every_real_given_case_declares_its_source_kind(case):
    """The eight worked cases, which were not written for this checker."""
    s = P.load(case)
    dim = ((s.raw.get("structure") or {}).get("dim") or {})
    if dim.get("basis") != "given":
        pytest.skip(f"{case.name} is basis={dim.get('basis')}")
    assert dim.get("source_kind") in P.DIM_GIVEN_KINDS, dim.get("source_kind")
    assert not _dim_errors(s), [str(i) for i in _dim_errors(s)]


def test_the_run_path_says_it_could_not_check_the_artefact():
    """★ `run.execute` rebuilds a PhysicalSystem around the SPEC's system
    document, and its `path` does not point at the case directory -- so the
    artefact's existence cannot be checked there.

    That must be a WARNING, not silence and not an error. Silence would let the
    run path print "structure OK" over an unverifiable claim, which is the shape
    `_phi_closure` already had to be fixed for; an error would refuse every real
    run, and this repository's pre-run gate once rejected 80 of 83 specs with
    zero real failures among them.
    """
    raw = yaml.safe_load((ROOT / "intake/trap-2d-5um/system.yaml").read_text())
    s = P.PhysicalSystem(path=pathlib.Path("nowhere/system.yaml"), raw=raw)
    issues = P.check_dim(s)
    w = [i for i in issues if i.where == "structure.dim.source_file"]
    assert w, [str(i) for i in issues]
    assert w[0].level == "warn", f"{w[0].level}: an invisible level is not saying so"
    assert "NOT verified" in w[0].msg, w[0].msg
    assert not [i for i in issues if i.level == "error"], \
        [str(i) for i in issues if i.level == "error"]


def test_the_full_gate_does_check_it():
    """...and the same document loaded from its real case directory DOES get
    the existence check, or the warning above would be the only behaviour and
    the check would exist nowhere."""
    s = P.load(ROOT / "intake/trap-2d-5um")
    assert not [i for i in s.issues if i.where == "structure.dim.source_file"], \
        "the real case emitted a source_file issue"
    #  and a broken path there is an error, not a warning
    raw = yaml.safe_load((ROOT / "intake/trap-2d-5um/system.yaml").read_text())
    raw["structure"]["dim"]["source_file"] = "ghost.jpeg"
    s2 = P.PhysicalSystem(path=ROOT / "intake/trap-2d-5um/system.yaml", raw=raw)
    e = [i for i in P.check_dim(s2) if i.where == "structure.dim.source_file"]
    assert e and e[0].level == "error", [str(i) for i in P.check_dim(s2)]
