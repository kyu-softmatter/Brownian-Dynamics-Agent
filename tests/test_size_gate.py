"""L2 — `structure.size`, what fixed N and which kind of N it is.

`particle.count` was carrying three different quantities under one name, read
out of the 8 cases on 2026-09-10:

    object       N is the size of the thing — a 25-bead chain
    system       N is how much of an infinite system is simulated
    replication  the physical system is one particle; N is a sampling multiplier

They carry opposite obligations, which is the whole reason for the field: the
same `converge --only N_double` is reassurance for `replication`, a category
error for `object` (n=50 is not the converged limit of n=25, it is a different
chain), and possibly a false negative for `system`.

**Sibling of `test_dim_gate.py`, and here for the same reason:** CI runs
`pytest` and nothing else, so `verify/` is a claim and `tests/` is enforcement.
`test_wiring` is the one to keep if any other is dropped.

⚠️ This file exists partly because the first check of the eight cases lied. The
loop counted `grep -c "✗ \\[structure"` and reported `errors=0` for all eight —
while every one of them was raising `NameError` and printing a traceback, in
which that pattern does not appear. "A check whose success is indistinguishable
from its own failure" (docs/05 §2), reproduced by hand. The tests below assert
on parsed `Issue` objects, not on rendered text.
"""
from __future__ import annotations

import copy
import pathlib

import pytest
import yaml

from bdbot import physical as P

ROOT = pathlib.Path(__file__).resolve().parent.parent
DONOR = ROOT / "intake/soft-r3-2d-A-sweep"          # role: system
CASES = sorted(p.parent for p in ROOT.glob("intake/*/system.yaml"))


@pytest.fixture(scope="module")
def donor_sys() -> dict:
    return yaml.safe_load((DONOR / "system.yaml").read_text())


@pytest.fixture
def case_dir(tmp_path, donor_sys):
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    (tmp_path / "observation.yaml").write_text((DONOR / "observation.yaml").read_text())

    def write(mutate=None):
        d = copy.deepcopy(donor_sys)
        if mutate is not None:
            mutate(d)
        (tmp_path / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
        return P.load(tmp_path)

    return write


def _size_errors(s) -> list:
    return [i for i in s.errors if i.where.startswith("structure.size")]


# ── the one that matters ───────────────────────────────────────────────────

def test_wiring(case_dir):
    """★ Reached through `physical.load` -> `validate`, not only by calling
    `check_size` directly. Fails if `out += check_size(s)` is removed."""
    s = case_dir(lambda d: d["structure"].pop("size"))
    assert _size_errors(s), "structure.size missing produced no error via load()"


def test_the_donor_passes_unmutated(case_dir):
    """Guard on the guard: every break test starts from this document."""
    s = case_dir()
    assert not _size_errors(s), [str(i) for i in _size_errors(s)]


def test_load_never_raises_on_a_broken_size_block(case_dir):
    """The checker must report, not crash. Feeding `furlong^2` as a unit once
    killed `bulk()` with a DimensionalityError instead of producing an error,
    and this file's own first run raised `NameError` on all eight cases while a
    grep-based check read the traceback as a pass.
    """
    for mutate in (lambda d: d["structure"].__setitem__("size", "not a mapping"),
                   lambda d: d["structure"]["size"].update(n="banana"),
                   lambda d: d["structure"]["size"].update(role=None),
                   lambda d: d["structure"]["size"].update(phi=[1, 2, 3])):
        s = case_dir(mutate)          # must not raise
        assert isinstance(s.errors, list)


# ── each rule fires ────────────────────────────────────────────────────────

@pytest.mark.parametrize("mutate,where", [
    (lambda d: d["structure"].pop("size"), "structure.size"),
    (lambda d: d["structure"]["size"].pop("role"), "structure.size.role"),
    (lambda d: d["structure"]["size"].update(role="biggish"), "structure.size.role"),
    (lambda d: d["structure"]["size"].pop("n"), "structure.size.n"),
    (lambda d: d["structure"]["size"].pop("fixed_by"), "structure.size.fixed_by"),
    (lambda d: d["structure"]["size"].pop("what_would_change"),
     "structure.size.what_would_change"),
    (lambda d: d["structure"]["size"].update(what_would_change="  "),
     "structure.size.what_would_change"),
    (lambda d: d["structure"]["size"].update(n=99), "structure.size.n"),
], ids=["no-size", "no-role", "bad-role", "no-n", "no-fixed_by", "no-wwc",
        "blank-wwc", "n-disagrees"])
def test_a_broken_field_is_caught(case_dir, mutate, where):
    s = case_dir(mutate)
    assert any(i.where == where for i in s.errors), \
        f"{where} not caught: {[str(i) for i in s.errors]}"


def test_fixed_by_is_validated_against_the_role(case_dir):
    """`minimum_image` is meaningful for a `system` and meaningless for an
    `object`, so the vocabulary is per-role. N is never chosen directly — the
    field names what was.
    """
    s = case_dir(lambda d: d["structure"]["size"].update(fixed_by="design_power"))
    assert any(i.where == "structure.size.fixed_by" for i in s.errors), \
        "a `object` reason was accepted on a `system` role"


# ── the coverage gap adversarial review found ──────────────────────────────

@pytest.mark.parametrize("count", [400, [5, 9, 15, 25], "400", 0.0],
                         ids=["scalar", "list", "string", "zero"])
def test_a_non_mapping_particle_count_is_reported_not_raised(case_dir, count):
    """★ The gap that let a real crash ship.

    `cnt = ((raw.get("particle") or {}).get("count") or {})` — `or {}` rescues
    only FALSY values, so the shorthand `count: 400` survived and was then
    dotted. `bdbot status` loops over every case with no guard (cli.py:170), so
    one hand-authored file took down the whole table, and HEAD had *reported*
    the same file cleanly. `test_load_never_raises_on_a_broken_size_block`
    above did not catch it because it only ever mutates `structure.size`.
    """
    #  ⚠ This assertion used to be `any(where.startswith("particle.count"))` —
    #     which `load()` emits with or without `check_size`, so the test passed
    #     with the gate stubbed out AND with the `or {}` regression restored.
    #     Round-2 review caught it. It now asserts the gate itself ran (it must
    #     report the missing/!= size block) and, above all, that nothing raised.
    for pop_size in (False, True):          # both call sites of `cnt`
        def mutate(d, pop=pop_size, c=count):
            d["particle"]["count"] = c
            if pop:
                d["structure"].pop("size")
        s = case_dir(mutate)                # ← the invariant: MUST NOT RAISE
        assert isinstance(s.errors, list)
        assert any(i.where.startswith("particle.count") for i in s.errors), \
            "the unreadable count node itself must still be reported"
        #  With `size` present the gate correctly says nothing more: `cnt` is {},
        #  so the n-cross-check has nothing to compare and the count node is
        #  already reported above. With `size` popped it must speak — that branch
        #  is the one that proves this test is not vacuous, and it is the branch
        #  where the original AttributeError was raised from the message itself.
        if pop_size:
            assert any(i.where.startswith("structure.size") for i in s.errors), \
                f"the gate did not run: {[str(i) for i in s.errors]}"


@pytest.mark.parametrize("container", ["particle", "geometry", "structure"])
def test_a_non_mapping_CONTAINER_is_reported_not_raised(case_dir, container):
    """★★ The scoping gap that let the same bug ship twice.

    Round 1 guarded `particle.count`; the ENCLOSING `(raw.get("particle") or {})`
    and `geo = raw.get("geometry") or {}` kept the idiom the fix comment itself
    condemns, and round 3 found a non-mapping `geometry:` raising from inside the
    new phi block — a strict regression for the 3 `role: system` cases, with the
    whole `bdbot status` table as blast radius.

    The test that was supposed to prevent this mutated `d["particle"]["count"]`,
    one level too deep — the same scoping mistake its own docstring diagnoses in
    `test_load_never_raises_on_a_broken_size_block`. This one mutates the
    CONTAINER, and is parametrised so the next one cannot be missed.
    """
    for bad in ("a string", 3, ["a"], 0.0):
        s = case_dir(lambda d, c=container, b=bad: d.__setitem__(c, b))
        assert isinstance(s.errors, list), f"{container}={bad!r} raised"


def test_a_non_mapping_count_does_not_crash_the_replication_branch(case_dir):
    """The second crash site: `except (KeyError, TypeError, ValueError)` did not
    list AttributeError, so a compliant `role: replication` block plus a scalar
    count still died."""
    def mutate(d):
        _as_replication(d)
        d["particle"]["count"] = 400
    s = case_dir(mutate)
    assert isinstance(s.errors, list)


@pytest.mark.parametrize("bad", [{"value": 1, "unit": "1"}, "one", None],
                         ids=["provenanced-dict", "word", "none"])
def test_a_non_numeric_replication_operand_is_reported(case_dir, bad):
    """The falsifiability check silently no-opped on anything float() refused —
    including the Provenanced `{value: ...}` form this repo writes everywhere
    else. A check that skips itself on malformed input has not checked."""
    s = case_dir(lambda d: _as_replication(d, n_physical=bad))
    assert _size_errors(s), "a non-numeric operand passed silently"


def test_pins_names_the_quantity_that_is_pinned(case_dir):
    """The field was called `follows` and held the PINNED quantity — the exact
    opposite of what the name reads, with nothing able to tell the two
    conventions apart. Renamed; this pins the meaning."""
    s = case_dir(lambda d: d["structure"]["size"].update(pins="phi"))
    assert any(i.where == "structure.size.pins" for i in s.errors)
    assert "PINS" in next(i.msg for i in s.errors if i.where == "structure.size.pins")


@pytest.mark.parametrize("status", ["done", "partial", "not_done"])
def test_finite_size_branches_on_a_status_word_not_prose(case_dir, status):
    """It was a substring search for "not done" in free prose, which fails in
    both directions. Only `done` may pass without a warning."""
    s = case_dir(lambda d: d["structure"]["size"].update(
        finite_size_status=status, finite_size="whatever the prose says"))
    warned = any("N-scaling" in i.msg for i in s.issues if i.level == "warn")
    assert warned == (status != "done"), f"{status}: warned={warned}"


@pytest.mark.parametrize("field,where", [("role", "structure.size.role"),
                                         ("basis", "structure.dim.basis")],
                         ids=["size-role", "dim-basis"])
def test_a_non_scalar_enum_field_is_reported_on_both_gates(case_dir, field, where):
    """Round 1 fixed the unhashable-`role` crash in `check_size` and left the
    identical `basis` line in `check_dim` alone; round 2 found it still raised.
    Both gates, one test, so the next enum field cannot be fixed on one side."""
    key = "size" if field == "role" else "dim"
    s = case_dir(lambda d: d["structure"][key].update(**{field: ["x"]}))
    assert any(i.where == where for i in s.errors), [str(i) for i in s.errors]


@pytest.mark.parametrize("phi,why", [("banana", "must be a number"), (0.9, "disagrees")],
                         ids=["not-a-number", "wrong-value"])
def test_phi_is_validated_and_cross_checked(case_dir, phi, why):
    """`phi` was required-present and never checked — `phi: banana` and a
    2.6x-wrong phi both validated clean. phi is the equation linking N and L, so
    an unchecked phi makes the whole `fixed_by` story unfalsifiable."""
    s = case_dir(lambda d: d["structure"]["size"].update(phi=phi))
    e = [i for i in s.errors if i.where == "structure.size.phi"]
    assert e, [str(i) for i in s.errors]
    assert why in e[0].msg


def test_phi_agreeing_with_geometry_passes(case_dir):
    """...and the check above is testing the comparison, not the key."""
    geo = case_dir().raw["geometry"]["area_fraction"]["value"]
    s = case_dir(lambda d: d["structure"]["size"].update(phi=geo))
    assert not [i for i in s.errors if i.where == "structure.size.phi"]


def test_stated_in_source_must_match_the_transcription(case_dir):
    """★ It was self-reported and checked against nothing, so the override gate
    was bypassed by writing the run's own `n` into it. The transcription is right
    there — `validate()` already loads `derived_from`, and `stated_quantities`
    carries `symbol: N` on the two cases whose sketches state it."""
    n = case_dir().raw["structure"]["size"]["n"]
    def mutate(d):
        d["structure"]["size"]["stated_in_source"] = n     # = the run's own N
        d["structure"]["size"].pop("approved_by", None)
    s = case_dir(mutate)
    e = [i for i in s.errors if i.where == "structure.size.stated_in_source"]
    assert e, [str(i) for i in s.errors]
    assert "stated_quantities" in e[0].msg


@pytest.mark.parametrize("bad", ["about a hundred", {"value": 100}],
                         ids=["prose", "mapping"])
def test_stated_in_source_must_be_comparable(case_dir, bad):
    """Prose in this field cannot be compared with `n`, so the override gate
    would never fire — the only mandatory size field whose value went unchecked."""
    s = case_dir(lambda d: d["structure"]["size"].update(stated_in_source=bad))
    e = [i for i in s.errors if i.where == "structure.size.stated_in_source"]
    assert e, [str(i) for i in s.errors]
    #  ⚠ THE FOURTH same-`where` collision in this file. Asserting only on the
    #     `where` passes with the type check deleted, because the transcription
    #     cross-check fires at the identical `where` (the donor's sketch says
    #     N=100; prose never equals 100). Assert the message this check owns.
    assert "must be a number" in e[0].msg, e[0].msg


def test_phi_must_close_against_N_and_the_box(case_dir):
    """★ Agreeing with `geometry.area_fraction` was not enough: phi = 3.5
    (350 % packing) passed both, because nothing recomputed phi from the box.
    phi, N and L are one equation; two of them fix the third."""
    def mutate(d):
        d["geometry"]["area_fraction"]["value"] = 3.5
        d["structure"]["size"]["phi"] = 3.5
    s = case_dir(mutate)
    e = [i for i in s.errors if i.where == "structure.size.phi"]
    assert e, [str(i) for i in s.errors]
    assert "box" in e[0].msg


def test_the_closure_survives_a_rectangular_box():
    """...and it must not fire on the one case that is deliberately NOT square.
    `trap-drag`'s aspect ratio is fixed at 1.0906 by commensurability, and the
    first version of the closure squared `box_length_x` and reported it.

    ⚠ This filtered `s.errors` only, so EVERY degradation path satisfied it —
    the closure not running at all looked identical to the closure running and
    agreeing. It now requires the closure to have actually RUN: no error AND no
    "NOT verified" / "NOT closed" warning on this case.
    """
    s = P.load(ROOT / "intake/trap-drag-2d-hex300")
    phi_iss = [i for i in s.issues if i.where == "structure.size.phi"]
    assert not phi_iss, f"the closure did not run cleanly: {[str(i) for i in phi_iss]}"


def test_the_partial_and_not_done_warnings_differ(case_dir):
    """The warning told the one case that records two FSS ladders that "no
    N-sweep is recorded". `partial` and `not_done` are different states."""
    msgs = {}
    for st in ("partial", "not_done"):
        s = case_dir(lambda d, v=st: d["structure"]["size"].update(finite_size_status=v))
        msgs[st] = next(i.msg for i in s.issues
                        if i.level == "warn" and "N-scaling" in i.msg)
    assert msgs["partial"] != msgs["not_done"]
    assert "PARTIAL" in msgs["partial"] and "NO" in msgs["not_done"]


def test_omitting_stated_in_source_is_itself_an_error(case_dir):
    """★ The override gate was OPT-IN: dropping the key dropped the
    `approved_by` obligation with it, so the one case that had to override a
    stated N could have hidden that by deleting one line. `null` is the answer
    when the source is silent; absence is not."""
    def mutate(d):
        d["structure"]["size"].pop("stated_in_source")
        d["structure"]["size"].pop("approved_by", None)
    s = case_dir(mutate)
    e = [i for i in s.errors if i.where == "structure.size.stated_in_source"]
    assert e, [str(i) for i in s.errors]
    #  ⚠ Asserting only on the `where` was VACUOUS: with the mandate removed from
    #     SIZE_ALWAYS the cross-check against the transcription fires at the SAME
    #     `where` (the donor's sketch says N=100, a missing field reads as None,
    #     mismatch) -- so the test passed for a different reason. Measured.
    #     Assert on the message the mandate itself produces.
    assert "write `null`" in e[0].msg, e[0].msg


def test_an_explicit_null_stated_in_source_is_accepted(tmp_path, donor_sys):
    """The six cases whose sketches are silent on N say so with `null`, and that
    must not be confused with the blank `test_omitting_...` rejects.

    This one cannot use the shared donor: soft-r3's transcription DOES state
    N = 100, and the cross-check would correctly reject a `null`. So the
    observation is stripped of its `N` here — which makes the test assert the
    pairing (silent source <-> null) rather than just the null.
    """
    import copy as _c
    obs = yaml.safe_load((DONOR / "observation.yaml").read_text())
    obs["stated_quantities"] = [q for q in obs["stated_quantities"]
                                if str(q.get("symbol")).strip() != "N"]
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    (tmp_path / "observation.yaml").write_text(yaml.safe_dump(obs, allow_unicode=True))
    d = _c.deepcopy(donor_sys)
    d["structure"]["size"]["stated_in_source"] = None
    d["structure"]["size"].pop("approved_by", None)
    (tmp_path / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    s = P.load(tmp_path)
    assert not [i for i in s.errors if i.where.startswith("structure.size")], \
        [str(i) for i in s.errors]


def test_a_silent_source_may_not_claim_a_stated_N(tmp_path, donor_sys):
    """The other half of the pairing: if the transcription says nothing about N,
    a non-null `stated_in_source` is a claim about the input that the input does
    not make (rule 3)."""
    import copy as _c
    obs = yaml.safe_load((DONOR / "observation.yaml").read_text())
    obs["stated_quantities"] = [q for q in obs["stated_quantities"]
                                if str(q.get("symbol")).strip() != "N"]
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    (tmp_path / "observation.yaml").write_text(yaml.safe_dump(obs, allow_unicode=True))
    d = _c.deepcopy(donor_sys)
    d["structure"]["size"]["stated_in_source"] = 100
    (tmp_path / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    s = P.load(tmp_path)
    assert any(i.where == "structure.size.stated_in_source" for i in s.errors)


@pytest.mark.parametrize("count,n,caught", [
    ([5, 9, 15, 25], [5, 9, 15, 50], True),    # same length, one value differs
    ([5, 9, 15, 25], [5, 9, 15],     True),    # ⚠ genuinely different LENGTHS
    ([5, 9, 15, 25], [5, 9, 15, 25], False),   # identical -> must pass
], ids=["different-values", "different-length", "identical"])
def test_a_list_valued_n_is_cross_checked(case_dir, donor_sys, count, n, caught):
    """3 of the 8 real cases carry list counts.

    ⚠ The "different-length" case used to pass `[5, 9]` against the donor's
    SCALAR count of 400 — caught as a type mismatch, never exercising the length
    branch of `_same_number`. Both sides are lists here.
    """
    import copy as _c
    def mutate(d):
        d["particle"]["count"] = _c.deepcopy(donor_sys["particle"]["count"])
        d["particle"]["count"]["value"] = count
        d["structure"]["size"]["n"] = n
    s = case_dir(mutate)
    assert bool([i for i in s.errors if i.where == "structure.size.n"]) == caught


def test_require_structure_refuses_a_spec_without_the_block(tmp_path, hoomd_mod):
    """★ The opt-in branch had ZERO behavioural coverage — deleting it left the
    whole suite green. It defaults to False (the 254 archived specs predate the
    block and must stay re-runnable), so this is the only thing that pins it.
    """
    import json
    from bdbot import nondim as ND, run as RUN
    src = sorted((ROOT / "specs").glob("trap-2d-5um__*.json"))[0]
    raw = json.loads(src.read_text())
    (raw.get("system") or {}).pop(P.STRUCTURE_SECTION, None)
    sp = tmp_path / "spec.json"; sp.write_text(json.dumps(raw))
    spec = ND.load(sp)

    def build_fn(_s, _o):
        raise AssertionError("build_fn reached — require_structure did not refuse")

    with pytest.raises(ValueError, match="no `structure` block"):
        RUN.execute(spec, build_fn, tmp_path / "run", progress=False,
                    require_structure=True)

    #  ...and the default stays permissive: the same spec must run
    reached = []
    def build_ok(_s, _o):
        reached.append(True); raise AssertionError("sentinel")
    with pytest.raises(AssertionError, match="sentinel"):
        RUN.execute(spec, build_ok, tmp_path / "run2", progress=False)
    assert reached, "require_structure=False must not block an archived spec"


@pytest.mark.parametrize("mutate,phrase", [
    (lambda d: [d["geometry"].pop(k, None) for k in
                ("box_length", "box_length_x", "box_length_y", "box_length_z")],
     "no box_length"),
    (lambda d: d.__setitem__("dimensions", "2D"), "`dimensions` is"),
    (lambda d: d.__setitem__("dimensions", [2]), "`dimensions` is"),
    (lambda d: d["geometry"]["box_length"].__setitem__("unit", "furlong^2"),
     "could not be closed"),
], ids=["no-box", "bad-dimensions-str", "bad-dimensions-list", "bad-unit"])
def test_an_unrunnable_closure_says_so_visibly(case_dir, mutate, phrase):
    """`_phi_closure` had FOUR paths that returned `[]` or an invisible `info`,
    so a phi that could not be closed read exactly like one that had been closed
    and agreed. `network` is the real instance of the first (3D, compressed, no
    box_length), and a non-numeric `dimensions` used to RAISE.

    ⚠ My earlier version of this test had a compound assertion containing
    `"not" in i.msg.lower()` — which matches "NOT verified" and equally matches
    a message saying the closure ran and found nothing wrong. It accepted a
    message claiming the opposite. Assert the level and the exact phrase.
    """
    s = case_dir(mutate)
    e = [i for i in s.issues if i.where == "structure.size.phi"]
    assert e, [str(i) for i in s.issues]
    assert e[0].level == "warn", f"{e[0].level}: an invisible level is not saying so"
    assert "NOT verified" in e[0].msg, e[0].msg
    #  ⚠ "NOT verified" alone was not enough: removing the `dim0 not in (2, 3)`
    #     guard leaves `dim0 = None`, which then raises inside the try and is
    #     absorbed by the generic `except` — producing the same phrase from a
    #     different path, so the guard was untested. Pin the phrase each path
    #     owns, so the specific reason has to survive.
    assert phrase in e[0].msg, f"expected {phrase!r} in: {e[0].msg}"


def test_network_phi_is_reported_as_unverified():
    """The real case, not a mutation: `network` is 3D and compressed, so no
    box_length exists and its phi cannot be closed. That must be a WARNING —
    an `info` is not surfaced by `checks.verdict` and no reader sees it."""
    s = P.load(ROOT / "intake/network")
    e = [i for i in s.issues if i.where == "structure.size.phi"]
    assert e, [str(i) for i in s.issues]
    assert e[0].level == "warn" and "NOT verified" in e[0].msg, str(e[0])


def test_the_run_path_warns_rather_than_claiming_validation():
    """★ `run.execute` builds a PhysicalSystem with an EMPTY `core`, so the phi
    closure cannot run there. It used to emit an invisible `info` while the gate
    printed "structure OK -- dim/size validated" — a 350 %-packing phi passed the
    run path in silence. The skipped check must be visible."""
    import copy as _c, pathlib as _pl
    raw = _c.deepcopy(yaml.safe_load((DONOR / "system.yaml").read_text()))
    raw["structure"]["size"]["phi"] = 3.5
    raw["geometry"]["area_fraction"]["value"] = 3.5
    ps = P.PhysicalSystem(path=_pl.Path("nonexistent/system.yaml"), raw=raw)
    e = [i for i in P.check_size(ps) if i.where == "structure.size.phi"]
    assert e and e[0].level == "warn" and "NOT closed" in e[0].msg, \
        [str(i) for i in e]


def test_a_matching_list_valued_n_passes(case_dir, donor_sys):
    """...and the list comparison is not simply rejecting every list."""
    import copy as _c
    lst = [5, 9, 15, 25]
    def mutate(d):
        d["particle"]["count"] = _c.deepcopy(donor_sys["particle"]["count"])
        d["particle"]["count"]["value"] = lst
        d["structure"]["size"]["n"] = lst
    s = case_dir(mutate)
    assert not [i for i in s.errors if i.where == "structure.size.n"]


def test_an_unknown_finite_size_status_is_an_error(case_dir):
    s = case_dir(lambda d: d["structure"]["size"].update(finite_size_status="mostly"))
    assert any(i.where == "structure.size.finite_size_status" for i in s.errors)


def test_overriding_a_stated_N_needs_a_human(case_dir):
    """A stated N may be overridden — `soft-r3`'s sketch says 100 and that value
    cannot satisfy `r_c < L/2` — but not silently."""
    s = case_dir(lambda d: d["structure"]["size"].update(
        stated_in_source=100, approved_by=None))
    assert any(i.where == "structure.size.approved_by" for i in s.errors)


def test_agreeing_with_the_stated_N_needs_no_approval(case_dir):
    """...and the check is testing the disagreement, not the presence of a key."""
    n = case_dir().raw["structure"]["size"]["n"]
    s = case_dir(lambda d: d["structure"]["size"].update(stated_in_source=n))
    assert not [i for i in s.errors if i.where == "structure.size.approved_by"]


def test_the_blank_what_would_change_message_says_why_it_matters(case_dir):
    """The sibling dim gate tests its copy of this message; size did not, so
    degrading it to the generic text failed nothing."""
    s = case_dir(lambda d: d["structure"]["size"].update(what_would_change="  "))
    msg = next(i.msg for i in s.errors if i.where == "structure.size.what_would_change")
    assert "not the same as it being safe" in msg


#: ⚠ SELF-REFERENTIAL PARAMETRISATION, the 6th instance of the vacuity class and
#: a new shape of it: the tests below used `parametrize(..., P.SIZE_*_REQUIRES)`,
#: so the case list WAS the subject. Deleting an element from the constant
#: deleted its own test case and the suite stayed green — measured: dropping
#: `phi` gave 1298 passed, and emptying the tuple gave 1295 passed / 3 skipped,
#: exit 0, while `check_size` on a `role: system` block with all four keys
#: removed reported NO errors. `phi` is the worst one, because every phi check
#: (the geometry cross-check AND `_phi_closure`) nests under `if phi is not
#: None`, so dropping it from the tuple kills them too.
#: Literal lists, plus a membership assertion, so the constant is pinned
#: independently of the tests that consume it.
OWED_BY_SYSTEM = ("phi", "pins", "finite_size", "finite_size_status")
OWED_BY_REPLICATION = ("n_physical", "n_replicas")


def test_the_owed_field_sets_are_what_this_file_thinks_they_are():
    """The pin. If a future edit narrows either constant, this fails by name
    instead of silently deleting the test cases below."""
    assert tuple(P.SIZE_SYSTEM_REQUIRES) == OWED_BY_SYSTEM, P.SIZE_SYSTEM_REQUIRES
    assert tuple(P.SIZE_REPLICATION_REQUIRES) == OWED_BY_REPLICATION, \
        P.SIZE_REPLICATION_REQUIRES
    assert tuple(P.SIZE_ROLES) == ("object", "system", "replication")
    assert tuple(P.SIZE_FS_STATUS) == ("done", "partial", "not_done")
    assert tuple(P.SIZE_PINS) == ("n", "L")


@pytest.mark.parametrize("field", OWED_BY_SYSTEM)
def test_system_owes_phi_follows_and_finite_size(case_dir, field):
    """phi links N and L but fixes neither, so `follows` records which one the
    second constraint pinned — 2 of the 3 real cases fix L, `trap-drag` fixes N.
    `finite_size` may say "not done"; it may not be blank.
    """
    s = case_dir(lambda d: d["structure"]["size"].pop(field))
    assert any(i.where == f"structure.size.{field}" for i in s.errors)


def test_a_missing_finite_size_study_warns_but_does_not_block(case_dir):
    """Whether the N-sweep was run is a property of the campaign, not of this
    file, so it warns. Measured 2026-09-10: only 1 of the 3 system-role cases
    has one, while `trap-drag` reports collective observables from 81 runs at a
    single size.
    """
    #  ⚠ The mutation used to change only the PROSE, which is a no-op on a donor
    #     that already warns — so the test never constructed the state it names.
    #     It now sets the status word (the thing that decides) and keeps the prose
    #     change as the thing that must NOT matter.
    s = case_dir(lambda d: d["structure"]["size"].update(
        finite_size_status="not_done", finite_size="single size only, no sweep"))
    assert not _size_errors(s), "a missing N-sweep must not block"
    assert any("N-scaling" in i.msg for i in s.issues
               if i.level == "warn"), [str(i) for i in s.issues]


# ── replication: the claim has to be falsifiable ───────────────────────────

def _as_replication(d, **over):
    e = {"role": "replication", "n": 400, "n_physical": 1, "n_replicas": 400,
         "fixed_by": "confirmed_by_run", "what_would_change": "nothing physical",
         #  the donor's observation.yaml states N = 100, and `stated_in_source`
         #  is cross-checked against it, so a synthetic block has to agree — and
         #  then needs the approval its disagreement with `n` demands.
         "stated_in_source": 100, "approved_by": "test fixture"}
    e.update(over)
    d["structure"]["size"] = e


def test_replication_product_must_match_particle_count(case_dir):
    """"The system is really one particle" is unfalsifiable unless
    `n_physical x n_replicas` is checked against the count the run uses.
    """
    s = case_dir(lambda d: _as_replication(d, n_replicas=7))
    assert any(i.where == "structure.size.n_replicas" for i in s.errors), \
        [str(i) for i in s.errors]


def test_a_correct_replication_block_passes(case_dir):
    """...and the check above is testing the product, not the presence of a key."""
    s = case_dir(_as_replication)
    assert not _size_errors(s), [str(i) for i in _size_errors(s)]


@pytest.mark.parametrize("field", OWED_BY_REPLICATION)
def test_replication_owes_the_split_it_asserts(case_dir, field):
    s = case_dir(lambda d: _as_replication(d, **{field: None}))
    assert any(i.where == f"structure.size.{field}" for i in s.errors)


def test_object_owes_neither_set(case_dir):
    """A role with no extra obligation must not be blocked by one — a gate that
    refuses everything is worse than no gate.
    """
    s = case_dir(lambda d: d["structure"].__setitem__("size", {
        "role": "object", "n": 400, "fixed_by": "swept",
        "stated_in_source": 100, "approved_by": "test fixture",
        "what_would_change": "n is the independent variable"}))
    assert not _size_errors(s), [str(i) for i in _size_errors(s)]


# ── the real tree ──────────────────────────────────────────────────────────

@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_every_real_case_declares_a_role(case):
    s = P.load(case)
    assert not _size_errors(s), [str(i) for i in _size_errors(s)]
    role = ((s.raw.get("structure") or {}).get("size") or {}).get("role")
    assert role in P.SIZE_ROLES, f"{case.name}: role={role!r}"


def test_the_roles_are_distributed_as_measured():
    """object=3 · system=3 · replication=2 on 2026-09-10. Not a constraint on
    the physics — a tripwire, so that a role silently flipping is noticed.
    """
    dist: dict = {}
    for c in CASES:
        r = ((P.load(c).raw.get("structure") or {}).get("size") or {}).get("role")
        dist[r] = dist.get(r, 0) + 1
    assert dist == {"object": 3, "system": 3, "replication": 2}, dist


def test_every_system_case_without_a_complete_n_sweep_says_so():
    """All three `system` cases warn. `soft-r3` was written `done` because the
    FSS ladders exist (campaigns/soft2d_fss.py, soft2d_nconv.py) — but A = 100,
    the production point its n = 400 was chosen for, has no N-sweep and the
    file's own CV3 never ran. Adversarial review caught that as reading the
    existence of a ladder as coverage of this point; it is now `partial`.
    An unrecorded gap and a recorded one must not look alike, and neither must
    a partial one and a complete one.
    """
    warned = sorted(c.name for c in CASES
                    if any("N-scaling" in i.msg for i in P.load(c).issues
                           if i.level == "warn"))
    assert warned == ["network", "soft-r3-2d-A-sweep", "trap-drag-2d-hex300"], warned


def test_no_case_currently_claims_a_complete_finite_size_study():
    """The honest state on 2026-09-10: 0 of 3. If one ever reaches `done`, this
    test fails and the claim gets looked at — which is the point.
    """
    done = [c.name for c in CASES
            if ((P.load(c).raw.get("structure") or {}).get("size") or {})
            .get("finite_size_status") == "done"]
    assert done == [], done


# ── the gate has to be able to stop a RUN, not only an inspection ──────────

def test_execute_carries_a_structure_gate():
    """★ Adversarial review found the gate blocked no run: `validate()` is
    reachable only from `bdbot.cli system check` and `status`, so a case script
    ran to PASS with the whole block missing — the `health.gate()` pattern the
    seal check was moved into `execute()` to escape.

    This asserts the wiring exists and carries the same shape as its two
    siblings: default False (the 254 archived specs predate the block and must
    stay re-runnable), but a block that IS present is always validated.
    """
    import inspect
    from bdbot import run as RUN
    sig = inspect.signature(RUN.execute)
    assert "require_structure" in sig.parameters
    assert sig.parameters["require_structure"].default is False
    src = inspect.getsource(RUN.execute)
    assert "check_dim" in src and "check_size" in src, \
        "execute() must call the gate itself, not delegate to a sibling tool"
    #  ⚠ This assertion alone is VACUOUS — it survives `if _struct:` being
    #     changed to `if False:` (measured: 71/71 still green). The test that
    #     actually catches the disablement is the next one, which calls
    #     execute() and requires the gate to fire before build_fn is reached.


def test_a_broken_structure_block_in_a_spec_is_an_error(donor_sys):
    """The always-check half: `execute()` rebuilds a PhysicalSystem around the
    spec's own `system` document and refuses on errors. Reproduced here at the
    level `execute()` uses, so it holds without hoomd or a build_fn.

    `core` is empty and `derived_from` does not resolve from a run directory —
    both of those checks degrade to SILENCE by design, so this also pins that
    they degrade rather than producing a false pass.
    """
    import copy as _c
    bad = _c.deepcopy(donor_sys)
    bad["structure"]["size"]["role"] = "biggish"
    ps = P.PhysicalSystem(path=ROOT / "nonexistent/system.yaml", raw=bad)
    errs = [i for i in (P.check_dim(ps) + P.check_size(ps)) if i.level == "error"]
    assert any(i.where == "structure.size.role" for i in errs), [str(i) for i in errs]

    good = P.PhysicalSystem(path=ROOT / "nonexistent/system.yaml", raw=donor_sys)
    clean = [i for i in (P.check_dim(good) + P.check_size(good)) if i.level == "error"]
    assert not clean, [str(i) for i in clean]


def test_execute_refuses_a_broken_structure_block_before_building(tmp_path, hoomd_mod):
    """★★ The non-vacuous half. `test_execute_carries_a_structure_gate` above
    only greps the source and survives the gate being turned off; this one calls
    `execute()` with a spec whose `structure` block is broken and requires the
    ValueError to arrive BEFORE `build_fn` runs.

    The gate sits before `build_fn(spec, outdir)`, so a build_fn that raises a
    sentinel proves the ordering: if the sentinel escapes, the gate did not fire.
    """
    import json
    from bdbot import nondim as ND, run as RUN

    src = sorted((ROOT / "specs").glob("trap-2d-5um__*.json"))[0]
    raw = json.loads(src.read_text())
    own = yaml.safe_load((ROOT / "intake/trap-2d-5um/system.yaml").read_text())
    raw.setdefault("system", {})["structure"] = copy.deepcopy(own["structure"])
    raw["system"]["structure"]["size"]["role"] = "biggish"      # the one break
    spec_path = tmp_path / "spec.json"
    spec_path.write_text(json.dumps(raw))
    spec = ND.load(spec_path)

    def build_fn(_spec, _outdir):
        raise AssertionError("build_fn was reached — the structure gate did not fire")

    with pytest.raises(ValueError) as ei:
        RUN.execute(spec, build_fn, tmp_path / "run", progress=False)
    assert "structure" in str(ei.value)
    assert "structure.size.role" in str(ei.value)


def test_execute_lets_a_valid_structure_block_through(tmp_path, hoomd_mod):
    """...and the check above is testing the block, not the presence of the key:
    the same spec with an intact block must reach build_fn."""
    import json
    from bdbot import nondim as ND, run as RUN

    src = sorted((ROOT / "specs").glob("trap-2d-5um__*.json"))[0]
    raw = json.loads(src.read_text())
    #  the case's OWN block — grafting another case's (n=400 vs count=1000)
    #  makes the `n` cross-check fire, which is the gate working, not a bug
    own = yaml.safe_load((ROOT / "intake/trap-2d-5um/system.yaml").read_text())
    raw.setdefault("system", {})["structure"] = copy.deepcopy(own["structure"])
    spec_path = tmp_path / "spec.json"
    spec_path.write_text(json.dumps(raw))
    spec = ND.load(spec_path)

    reached = []

    def build_fn(_spec, _outdir):
        reached.append(True)
        raise AssertionError("sentinel: reached build_fn")

    with pytest.raises(AssertionError, match="sentinel"):
        RUN.execute(spec, build_fn, tmp_path / "run", progress=False)
    assert reached, "a valid structure block must not be blocked"


# ── the two round-6 fixes my own break-test found untested ─────────────────

@pytest.mark.parametrize("bad", ["2D", [2], {"value": 2}, "two"],
                         ids=["str", "list", "dict", "word"])
def test_render_check_does_not_crash_on_a_non_numeric_dimensions(tmp_path,
                                                                 donor_sys, bad):
    """★ CLASS A instance #5, and the test for it.

    `PhysicalSystem.dim` was `int(self.raw.get("dimensions", 0))`. `load()` was
    fixed to report, but `render_check` READS `s.dim`, so
    `bdbot.cli system check` still tracebacked on exactly the value the gate was
    built to report. `_int_or_none` had existed for a round by then and had been
    applied only to the site that was found — which is the shape this class
    keeps taking.
    """
    import copy as _c
    d = _c.deepcopy(donor_sys)
    d["dimensions"] = bad
    (tmp_path / "sketch_01.jpeg").write_bytes(b"")
    (tmp_path / "observation.yaml").write_text((DONOR / "observation.yaml").read_text())
    (tmp_path / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    s_ = P.load(tmp_path)                       # must not raise
    txt = P.render_check(s_)                    # ← the path that still crashed
    assert "VERDICT: FAIL" in txt
    assert any(i.where.startswith("structure.dim") for i in s_.errors), \
        [str(i) for i in s_.errors]


def test_the_run_gate_reads_n_beads_too(tmp_path, hoomd_mod):
    """★ The gate's N lookup was `params.get("N", params.get("n_particles"))`
    and missed `n_beads` — which **183 of the 278** archived specs use, i.e. the
    chain cases, whose `n` IS the chain length. The comparison silently did not
    happen for the majority of the archive.
    """
    import json
    from bdbot import nondim as ND, run as RUN
    assert "n_beads" in RUN._N_PARAM_KEYS

    src = sorted((ROOT / "specs").glob("chain-bend-2d-dlvo__n9-*.json"))[0]
    raw = json.loads(src.read_text())
    assert "n_beads" in (raw.get("params") or {}), sorted(raw.get("params") or {})
    own = yaml.safe_load((ROOT / "intake/chain-bend-2d-dlvo/system.yaml").read_text())
    raw.setdefault("system", {})["structure"] = copy.deepcopy(own["structure"])
    sp = tmp_path / "spec.json"; sp.write_text(json.dumps(raw))
    spec = ND.load(sp)

    seen = []
    def build_fn(_s, _o):
        seen.append(True); raise AssertionError("sentinel")
    import io, contextlib
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf), pytest.raises(AssertionError, match="sentinel"):
        RUN.execute(spec, build_fn, tmp_path / "run", progress=False)
    out = buf.getvalue()
    assert seen, "a valid block must not be blocked"
    #  n_beads = 9 is a member of the declared [5, 9, 15, 25], so it must report
    #  a COMPARISON, not "carries no N"
    assert "carries no N" not in out, out
    assert "N compared against spec.params (9)" in out, out


# ── the round-6 mutation survivors, and the amplifier ──────────────────────

def test_a_partial_box_says_so(case_dir):
    """Mutation survivor: the `len(per_axis) < dim0` path. The parametrised
    closure test covers no-box / bad-dimensions x2 / bad-unit but not a box with
    SOME axes recorded, so that warn had no test."""
    def mutate(d):
        g = d["geometry"]
        g.pop("box_length", None)
        g["box_length_x"] = {"value": 30.0, "unit": "um", "source": "t", "tier": 3}
        g.pop("box_length_y", None)          # 1 of 2 axes in 2D
    s_ = case_dir(mutate)
    e = [i for i in s_.issues if i.where == "structure.size.phi"]
    assert e and e[0].level == "warn", [str(i) for i in s_.issues]
    assert "of 2 box axes" in e[0].msg, e[0].msg


def test_a_missing_geometry_fraction_says_so(case_dir):
    """Mutation survivor, and a 5th instance of the CLASS C shape: the
    `node is None` warn ("no area_fraction ... recorded, not verified") had no
    test, so demoting it to `info` — the exact mistake the round closed
    elsewhere — was green."""
    def mutate(d):
        for k in ("area_fraction", "volume_fraction_final", "volume_fraction"):
            d["geometry"].pop(k, None)
    s_ = case_dir(mutate)
    e = [i for i in s_.issues if i.where == "structure.size.phi"]
    assert e and e[0].level == "warn", [str(i) for i in s_.issues]
    #  one phrase across every degradation path, so the assertion does not have
    #  to know which path produced it (the casing used to differ per path)
    assert "NOT verified" in e[0].msg, e[0].msg


def test_finite_size_status_done_validates_cleanly(case_dir):
    """Mutation survivor `FSSTATUS-del-done`: no test asserted that the
    ACCEPTING state validates. `test_finite_size_branches_on_a_status_word...`
    only checks `warned == (status != "done")`, which stays true if `done`
    becomes an outright error, and `test_no_case_currently_claims_a_complete_
    finite_size_study` actively asserts no case uses it — so the first case to
    finish an N-sweep would have been the first reader of that path.
    """
    assert "done" in P.SIZE_FS_STATUS
    s_ = case_dir(lambda d: d["structure"]["size"].update(finite_size_status="done"))
    assert not _size_errors(s_), [str(i) for i in _size_errors(s_)]
    assert not [i for i in s_.issues if i.level == "warn" and "N-scaling" in i.msg]


def test_the_status_table_survives_one_unreadable_case(tmp_path, monkeypatch, capsys):
    """★★ THE AMPLIFIER. `cmd_status` looped every case with no guard and printed
    only after the loop, so ONE unreadable `system.yaml` took down all 8 rows —
    and six rounds of review closed coercion sites one at a time without
    touching the line that converts any single site into total loss. ~145
    provenance sites in `load_node`/`render_check` can still raise and they
    predate this branch, so closing them one by one is not a strategy.
    """
    import shutil
    from bdbot import cli
    (tmp_path / "intake").mkdir()
    for d in sorted((ROOT / "intake").iterdir()):
        if d.is_dir():
            shutil.copytree(d, tmp_path / "intake" / d.name)
    bad = tmp_path / "intake" / "soft-r3-2d-A-sweep" / "system.yaml"
    raw = yaml.safe_load(bad.read_text())
    raw["particle"]["diameter"] = {"value": [1, 2], "unit": "um",
                                   "source": "t", "tier": 0}   # pint refuses
    bad.write_text(yaml.safe_dump(raw, allow_unicode=True))

    monkeypatch.setattr(cli, "_cases",
                        lambda: sorted(p for p in (tmp_path / "intake").iterdir()
                                       if p.is_dir()))
    cli.cmd_status(type("A", (), {})())          # must not raise
    out = capsys.readouterr().out
    assert "soft-r3-2d-A-sweep" in out and "ERROR" in out
    #  ...and the other seven rows survive, which is the whole point
    for name in ("abp-rod-2d-run-flip", "chain-bend-2d-dlvo", "network",
                 "trap-2d-5um", "trap-drag-2d-hex300"):
        assert name in out, f"{name} lost from the table"
    assert "READY" in out


def test_the_run_gate_refuses_rather_than_tracebacks_on_an_unreadable_block(tmp_path,
                                                                           hoomd_mod):
    """The gate's own refusal path was unguarded: a block it could not even
    validate produced a traceback out of `execute()` instead of the ValueError.
    `verify_hash()` does not cover `structure` (it is a DOC_KEY), so nothing
    upstream guarantees the block is well-formed."""
    import json
    from bdbot import nondim as ND, run as RUN
    src = sorted((ROOT / "specs").glob("trap-2d-5um__*.json"))[0]
    raw = json.loads(src.read_text())
    own = yaml.safe_load((ROOT / "intake/trap-2d-5um/system.yaml").read_text())
    raw.setdefault("system", {})["structure"] = copy.deepcopy(own["structure"])
    #  a shape that makes the CHECKER itself fall over, not one it reports
    raw["system"]["structure"]["size"] = {"role": "system", "n": 1000,
                                          "fixed_by": "minimum_image",
                                          "what_would_change": "x",
                                          "stated_in_source": None,
                                          "phi": 0.35, "pins": "L",
                                          "finite_size": "x",
                                          "finite_size_status": "not_done"}
    raw["system"]["geometry"] = {"box_length": {"value": "wide", "unit": "um"}}
    sp = tmp_path / "spec.json"; sp.write_text(json.dumps(raw))
    spec = ND.load(sp)

    def build_fn(_s, _o):
        raise AssertionError("build_fn reached")
    #  either it refuses (ValueError) or it validates and reaches build_fn --
    #  what it must NEVER do is raise something else out of execute()
    with pytest.raises((ValueError, AssertionError)) as ei:
        RUN.execute(spec, build_fn, tmp_path / "run", progress=False)
    assert isinstance(ei.value, (ValueError, AssertionError)), type(ei.value)
