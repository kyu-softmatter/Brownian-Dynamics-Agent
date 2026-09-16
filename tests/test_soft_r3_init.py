"""`soft-r3`'s crystal-start arm — the two knobs, and what must not move.

The melting bracket needs a run started from a perfect hexagonal crystal. That
needs a commensurate rectangular box, which means the initial condition and the
box shape are two knobs and not one: switching only `init` would change the box
as well, and the comparison would confound them. `simbot.run` reached the same
conclusion independently and refuses `init="hex"` with `box_shape="square"`.

The invariant these tests exist to protect is that **no archived run is
renamed**. `params` is hashed, so writing `"init": "rsa"` on the default path
would re-id all nine soft-r3 specs and their run directories.
"""
from __future__ import annotations

import copy
import json
import math
import pathlib
import struct
import subprocess
import sys

import pytest

from bdbot import runid as RID
from bdbot.pairpot import a_mean_star

ROOT = pathlib.Path(__file__).resolve().parent.parent
CASE = ROOT / "cases/soft_r3_2d.py"
ARCHIVED = "soft-r3-2d-A-sweep__A100__30caa5c9e0"


def _run(*args, **kw):
    return subprocess.run([sys.executable, str(CASE), *args],
                          capture_output=True, text=True, cwd=ROOT, **kw)


def _cleanup(rundir) -> None:
    """Remove a run directory and its spec. A test must not leave the working
    tree dirty -- `git status` after `pytest` has to stay clean, or a real
    change is invisible among the residue."""
    import shutil
    if rundir is None:
        return
    spec = ROOT / "specs" / f"{rundir.name}.json"
    shutil.rmtree(rundir, ignore_errors=True)
    spec.unlink(missing_ok=True)


# ── the forbidden combination ──────────────────────────────────────────────

def test_a_hex_lattice_in_a_square_box_is_refused():
    """It does not tile, so the periodic images do not match at the seam — and a
    seam is a line of defects that melts on its own. In a melting experiment that
    would answer the question in the wrong direction."""
    r = _run("--A", "34.938", "--init", "hex", "--box", "square", "--report")
    assert r.returncode != 0
    assert "--init hex needs --box hex" in (r.stderr + r.stdout)


def test_the_matched_random_arm_is_allowed():
    """★ The reverse is not merely permitted, it is required: the random arm has
    to sit in the SAME box as the crystal arm. A gate that refused this would
    force the confound it exists to prevent."""
    r = _run("--A", "34.938", "--box", "hex", "--init", "rsa", "--report")
    assert r.returncode == 0, r.stderr[-800:]
    assert "box20x20" in r.stdout


# ── run_id stability ───────────────────────────────────────────────────────
#
# ⚠ **The literal digest is not portable, and that is measured rather than
#   suspected.** `params.Gamma = A / a_mean^3` with `a_mean = sqrt(pi/(4 phi))`.
#   `sqrt` is correctly rounded by IEEE-754 and is therefore identical
#   everywhere; `pow` is **not required to be**, and the two libms disagree:
#
#       a_mean          1.4979969134027407     identical on both
#       a_mean**3       3.361497213033026      Apple libm,  osx-arm64
#       a_mean*a*a      3.3614972130330254     glibc + exact mults, linux-64
#       Gamma(A=100)   29.748648790272707  vs 29.74864879027271   (1 ULP)
#
#   `runid.spec_hash` serialises floats at full `repr` precision, so one ULP
#   renames the run: `__A100__30caa5c9e0` on osx-arm64 against
#   `__A100__079a25f073` on linux-64. Measured on CI run 35026488825, and
#   localised by perturbing each of the payload's 30 float leaves by one ULP in
#   turn -- `params.Gamma` is the only leaf that reproduces the linux digest.
#
#   **The archive cannot be re-identified to fix this.** The sealed
#   pre-registration `campaigns/s30_preregistration/prediction.yaml` cites
#   `runs/soft-r3-2d-A-sweep__A100__30caa5c9e0` by name, and its sha256 is locked
#   into 12 `SEALED.sha256` files; renaming means editing a sealed document after
#   the fact, which is the one thing pre-registration exists to prevent. 104 of
#   296 specs carry a hashed `params.Gamma` and all 104 move under a 1-ULP shift.
#   So the platform dependence is pinned here instead of papered over.
#
#   `LoadedSpec.verify_hash()` is **unaffected**: it re-hashes STORED content, so
#   it is portable, and it passed on linux-64. Only re-deriving a spec from the
#   physics is platform-dependent.

def _archived_payload() -> dict:
    """The hashed payload of the archived spec, read off the artefact."""
    spec = json.loads((ROOT / "specs" / f"{ARCHIVED}.json").read_text())
    return {"system": RID.physics_only(spec.get("system", {})),
            "params": RID.physics_only(spec["params"]),
            "numerics": RID.physics_only(spec["numerics"])}


def _float_leaves(node, path=()):
    if isinstance(node, dict):
        for k, v in node.items():
            yield from _float_leaves(v, path + (k,))
    elif isinstance(node, list):
        for i, v in enumerate(node):
            yield from _float_leaves(v, path + (i,))
    elif isinstance(node, float):
        yield path


def _ulp_neighbourhood(payload: dict, nhex: int, radius: int = 1) -> set:
    """Every digest the archived payload yields when ONE float leaf lands on an
    adjacent double -- i.e. the set of run_ids meaning *the same physics,
    computed against a different libm*.

    One leaf at a time, because that is what was measured. Two leaves diverging
    together would fail the test that uses this, and that is the intended
    behaviour: it would be a new fact and should be measured, not absorbed.

    ⚠ `radius` is 1, and it was 2 for one commit. At 2 the accepted set is 121
      digests rather than 61 while every name, docstring and failure message in
      this file said "one ULP" -- undocumented slack, and slack in exactly the
      wrong place: `pow`'s permitted error is +-1 ULP around the correctly
      rounded result and this archive sits at one edge of that band (Apple libm
      = glibc + 1 ULP), so a third libm at the *opposite* edge is two ULP from
      the archive. That is precisely the divergence this file promises will
      "fail, which is news rather than noise", and radius=2 absorbed it
      silently. Both measured digests lie at distance 0 and 1, so 1 is enough.
    """
    payload = copy.deepcopy(payload)
    leaves = list(_float_leaves(payload))
    assert leaves, "no float leaves -- this is not the payload the test means"

    def at(path, value=None):
        d = payload
        for k in path[:-1]:
            d = d[k]
        if value is None:
            return d[path[-1]]
        d[path[-1]] = value

    out = {RID.spec_hash(payload, nhex)}
    for path in leaves:
        orig = at(path)
        for direction in (math.inf, -math.inf):
            y = orig
            for _ in range(radius):
                y = math.nextafter(y, direction)
                at(path, y)
                out.add(RID.spec_hash(payload, nhex))
        at(path, orig)
    return out


def _ulp_distance(x: float, y: float) -> int:
    ix, iy = (struct.unpack("<q", struct.pack("<d", v))[0] for v in (x, y))
    return abs(ix - iy)


def test_the_archived_run_id_does_not_move():
    """The default path must still produce the archived run_id -- to within the
    last bit of one float, which is all the portability `pow` offers.

    If `init` were written as "rsa" rather than left absent, this is what catches
    it: adding a key to `params` moves the digest somewhere no single-ULP
    perturbation of the archived payload can reach, and nine specs and their run
    directories would be renamed.

    ★ This asserts strictly MORE than the literal string it replaces. The string
    said "the digest is 30caa5c9e0" -- true only on the machine the archive was
    built on. This says "the payload IS the archived payload, to within one ULP
    of one float, and in no other respect", which holds on every platform.
    """
    r = _run("--A", "100", "--report")
    assert r.returncode == 0, r.stderr[-800:]
    got = next(ln for ln in r.stdout.splitlines()
               if "run_id=" in ln).split("run_id=")[1].strip()

    assert got.startswith("soft-r3-2d-A-sweep__A100__"), got
    digest = got.rsplit("__", 1)[1]
    payload = _archived_payload()
    ok = _ulp_neighbourhood(payload, len(digest))
    assert digest in ok, (
        f"{got} is not the archived payload. No one-ULP perturbation of any of "
        f"{len(list(_float_leaves(payload)))} float leaves in {ARCHIVED} reaches "
        f"{digest}, so something changed in `params` or `numerics` -- which "
        f"renames nine specs and their run directories.")


def test_both_measured_platform_digests_are_the_same_physics():
    """Pin both measured values: `30caa5c9e0` on osx-arm64 (the archive) and
    `079a25f073` on linux-64 (CI run 35026488825). If a libm update moves the
    digest further than one ULP this fails, which is news rather than noise."""
    ok = _ulp_neighbourhood(_archived_payload(), 10)
    assert "30caa5c9e0" in ok, "the archive's own digest is not reproducible"
    assert "079a25f073" in ok, "the linux-64 digest is no longer one ULP away"


def test_the_ulp_tolerance_still_refuses_what_it_exists_to_catch():
    """★ Guard on the guard. A tolerance that accepts everything is not a
    tolerance, and this repository's record on unexercised checks is bad enough
    that the discriminating power is measured rather than assumed.

    Each mutation below is a real change to the hashed payload -- the first two
    are exactly the two that survived this file's first version -- and each must
    land OUTSIDE the one-ULP neighbourhood.

    ⚠ The last two were measured by hand when this fix was written, quoted as
      CAUGHT in `docs/05-pitfalls.md` and in the wiki finding, and then not
      wired here: the mutation list held four entries while three documents
      cited six, and one document named *this test* as where the `Gamma`
      rounding was caught. The word `Gamma` did not appear in this function.
      That is presence-not-validity committed in the prose rather than the code,
      which is worse than the original defect because the prose is the part a
      reader trusts. A mutation takes a callable now, because "drop the key"
      cannot be expressed by `dict.update`.
    """
    payload = _archived_payload()
    ok = _ulp_neighbourhood(payload, 10)

    def _set(**kw):
        return lambda params: params.update(kw)

    mutations = [
        ("init written on the default path", _set(init="rsa")),
        ("lattice written on the default path", _set(n_x=20, n_y=20)),
        ("phi moved in the 10th digit", _set(phi=0.3500000001)),
        ("A moved in the 12th digit", _set(A=100.00000000001)),
        ("Gamma dropped from params", lambda params: params.pop("Gamma")),
        ("Gamma rounded to 12 significant figures",
         lambda params: params.__setitem__("Gamma", float(f"{params['Gamma']:.12g}"))),
    ]
    for why, mutate in mutations:
        m = copy.deepcopy(payload)
        mutate(m["params"])
        assert RID.spec_hash(m, 10) not in ok, f"accepted: {why}"

    # ...and the boundary, which must be a measurement and not decoration.
    # `len(ok) <= 2*n_leaves + 1` is the construction's arithmetic maximum and
    # can therefore never fail; assert the exact size instead, and that the
    # first shift OUTSIDE the radius is refused.
    n_leaves = len(list(_float_leaves(payload)))
    assert len(ok) == 2 * n_leaves + 1, (len(ok), n_leaves)
    beyond = copy.deepcopy(payload)
    g = beyond["params"]["Gamma"]
    for _ in range(2):
        g = math.nextafter(g, math.inf)
    beyond["params"]["Gamma"] = g
    assert RID.spec_hash(beyond, 10) not in ok, "a two-ULP shift is accepted"


def test_pow_is_what_moves_and_it_moves_by_one_ulp():
    """The cause, asserted portably on both platforms.

    `sqrt` is correctly rounded by IEEE-754, so `a_mean` is bit-identical
    everywhere. `pow` carries no such guarantee, so `a**3` and `a*a*a` are
    *allowed* to differ -- they do on osx-arm64 and do not on linux-64. What is
    true on both is that they agree to 15 significant figures and differ by at
    most one ULP.

    ⚠ **Do not generalise the 15 from here.** This file's first version said a
      payload serialised at 15 significant figures "would have been portable",
      which is true for `params.Gamma` (0 of 104 specs disagree at 15 digits)
      and false for `params.k_bond_star`, hashed in 186 specs:
      `cases/chain_bend_dlvo_2d.py`'s `find_well` takes a central second
      difference with `dh = h_min*1e-4`, so one ULP in a single `U_star`
      evaluation becomes 2.96e-9 *relative* in `k_bond_star` -- measured
      1042362.8817700658 against 1042362.8848514813, which agree at 9
      significant figures and disagree at 12 and at 15. The portable digit count
      is set by each field's conditioning, not by the payload.
    """
    a = a_mean_star(0.35)
    assert a == math.sqrt(math.pi / (4 * 0.35))       # exact, everywhere
    p, m = 100.0 / a**3, 100.0 / (a * a * a)
    assert p == pytest.approx(m, rel=1e-15)
    assert _ulp_distance(p, m) <= 1, (repr(p), repr(m))

    # and one ULP is enough to rename the run -- the reason any of this matters
    payload = _archived_payload()
    shifted = copy.deepcopy(payload)
    shifted["params"]["Gamma"] = math.nextafter(payload["params"]["Gamma"],
                                                -math.inf)
    assert RID.spec_hash(shifted, 10) != RID.spec_hash(payload, 10)


def test_the_archived_spec_is_still_named_by_its_own_run_id():
    """...asserted on the artefact rather than on stdout, because regenerating
    the spec rewrites the file. The run_id must survive that rewrite."""
    p = ROOT / "specs/soft-r3-2d-A-sweep__A100__30caa5c9e0.json"
    assert json.loads(p.read_text())["run_id"] == p.stem


def _run_id(*args) -> str:
    r = _run("--A", "34.938", "--N", "400", *args, "--report")
    assert r.returncode == 0, r.stderr[-800:]
    line = next(ln for ln in r.stdout.splitlines() if "run_id=" in ln)
    return line.split("run_id=")[1].strip()


@pytest.mark.parametrize("args,expect", [
    (("--box", "hex", "--init", "hex"), "hex20x20"),
    (("--box", "hex",), "box20x20"),
    (("--seed", "7"), "-s7"),
])
def test_each_knob_reaches_the_tag(args, expect):
    """The directory name says which arm it is, without opening the spec."""
    assert expect in _run_id(*args)


@pytest.mark.parametrize("args", [
    ("--box", "hex", "--init", "hex"),
    ("--box", "hex"),
    ("--seed", "7"),
])
def test_each_knob_reaches_the_HASH_not_only_the_tag(args):
    """★★ The tag is cosmetic; the hash is the content address.

    ⚠ The first version of this file asserted on the tag alone, and two
    mutations survived it: deleting `init` from `params`, and deleting
    `n_x`/`n_y` from `params`. Both leave the tag intact, so the directory names
    still differ while two physically different runs share a hash -- which
    breaks "same spec -> same run_id -> do not re-run" and means `verify_hash`
    would not notice a swapped initial condition. Measured: 0 tests failed.
    Assert on the part that carries the physics.
    """
    base = _run_id().rsplit("__", 1)[1]
    got = _run_id(*args).rsplit("__", 1)[1]
    assert got != base, (
        f"hash {got} is unchanged by {args} -- the knob is not in `params`, so "
        f"the spec does not record which run this is")


def test_the_two_arms_differ_only_in_the_initial_condition_and_still_re_id():
    """The sharpest case: same box, same N, same seed, same everything except
    the positions. If `init` is not hashed these two collide."""
    hexes = _run_id("--box", "hex", "--init", "hex").rsplit("__", 1)[1]
    rsa = _run_id("--box", "hex", "--init", "rsa").rsplit("__", 1)[1]
    assert hexes != rsa, "the crystal arm and the matched random arm share a hash"


# ── the minimum image is read off the short side ───────────────────────────

def test_the_min_image_check_uses_the_short_side_of_the_commensurate_box():
    """★ The commensurate box has the same AREA as the square one but is 7.5 %
    shorter in y, and the minimum image is set by the short side. Checking
    against the area-equivalent square passes a cutoff that is actually too
    large. Both numbers are printed so the optimism is visible, not just fixed.
    """
    r = _run("--A", "34.938", "--N", "400", "--rc-shells", "7.8",
             "--box", "hex", "--init", "hex", "--report")
    assert r.returncode == 0, r.stderr[-800:]
    line = next(ln for ln in r.stdout.splitlines() if "r_c/(L/2)" in ln)
    assert "8.382e-01" in line, line          # 11.684 / (27.881/2)
    assert "0.7800" in r.stdout               # what the square box would have said
    assert "7.5 % optimistic" in r.stdout


def test_the_square_box_reading_is_unchanged():
    """Guard on the guard: the default path must still report the square-box
    ratio, or the test above would pass for a check that simply moved."""
    r = _run("--A", "34.938", "--N", "400", "--rc-shells", "7.8", "--report")
    assert r.returncode == 0, r.stderr[-800:]
    line = next(ln for ln in r.stdout.splitlines() if "r_c/(L/2)" in ln)
    assert "7.800e-01" in line, line
    assert "optimistic" not in line


# ── the artefact path, which is where a melting run actually broke ─────────

@pytest.mark.slow
def test_a_run_with_no_equilibration_phase_writes_its_artefacts():
    """★ `--eq-frac 0` crashed AFTER the verdict had printed PASS.

    The equilibration phase is built with `collect=False`, so a melting run
    needs `--eq-frac 0` to record from t = 0 -- and then `observables.npz` has
    no `eq_trace` array, which `make_plots` panel 3 read unconditionally. The
    physics was complete and correct; the `KeyError` arrived during artefact
    writing, so `result.txt` was never written, the exit code was 1, and the
    campaign driver read that as a failed run and stopped. It cost a 12-run
    sealed campaign its first run, measured 2026-09-14.

    So the assertion is on `result.txt`, not on the exit code alone: that file
    is what `RID.prepare_outdir` treats as the completion marker, and its
    absence is what made a passing run look failed.
    """
    #  ⚠ its own seed, and cleaned up. Without that the test runs `--force` over
    #     a COMMITTED smoke directory, so every `pytest` left four tracked files
    #     modified and `git status` was never clean after a test run.
    rundir = None
    try:
        r = _run("--A", "34.938", "--N", "144", "--rc-shells", "5", "--smoke",
                 "--box", "hex", "--init", "hex", "--eq-frac", "0",
                 "--seed", "424242")
        assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-1500:]
        line = next(ln for ln in r.stdout.splitlines() if "run_id=" in ln)
        rundir = ROOT / "runs" / line.split("run_id=")[1].strip()
        for name in ("result.txt", "metrics.json", "observables.npz",
                     "observables.png"):
            assert (rundir / name).exists(), f"{name} was not written to {rundir}"
        import numpy as np
        res = np.load(rundir / "observables.npz")
        assert "eq_trace" not in res.files, \
            "an eq_frac=0 run grew an eq_trace -- this no longer covers the bug"
        assert "psi6_global" in res.files, "the decision statistic is not on disk"
    finally:
        _cleanup(rundir)


@pytest.mark.slow
def test_the_default_run_still_has_an_equilibration_trace():
    """Guard on the guard: the branch must be a branch. If `eq_trace` vanished
    from every run, the test above would pass for the wrong reason."""
    rundir = None
    try:
        r = _run("--A", "34.938", "--N", "144", "--rc-shells", "5", "--smoke",
                 "--box", "hex", "--init", "hex", "--seed", "424243")
        assert r.returncode == 0, r.stdout[-1500:]
        line = next(ln for ln in r.stdout.splitlines() if "run_id=" in ln)
        rundir = ROOT / "runs" / line.split("run_id=")[1].strip()
        import numpy as np
        res = np.load(rundir / "observables.npz")
        assert "eq_trace" in res.files, \
            "the default path lost its equilibration trace"
    finally:
        _cleanup(rundir)


def test_report_is_read_only():
    """★ `--report` must not write to `specs/`.

    It used to. Every invocation left a spec file behind, and once this file
    grew 13 tests that call `--report` to read a run_id off stdout, running the
    suite added one spec per test. Measured: 17 of the specs committed in a
    single session had no run behind them, most of them test residue, in a
    directory the README presents as the artefact ledger.

    `--spec` is the path that writes. This test is what keeps the two apart.
    """
    before = {p.name for p in (ROOT / "specs").glob("*.json")}
    r = _run("--A", "34.938", "--N", "400", "--box", "hex", "--init", "hex",
             "--seed", "999", "--report")
    assert r.returncode == 0, r.stderr[-800:]
    after = {p.name for p in (ROOT / "specs").glob("*.json")}
    assert after == before, f"--report wrote {sorted(after - before)}"


def test_spec_still_writes():
    """Guard on the guard: `--spec` must still write, or the test above passes
    for a case script that can no longer produce a spec at all -- and the
    campaign driver reads run_ids from exactly that path."""
    r = _run("--A", "34.938", "--N", "400", "--box", "hex", "--init", "hex",
             "--seed", "998", "--spec")
    assert r.returncode == 0, r.stderr[-800:]
    line = next(ln for ln in r.stdout.splitlines() if "L3 spec:" in ln)
    p = ROOT / line.split("L3 spec:")[1].strip()
    assert p.exists(), p
    p.unlink()                      # a probe, not an artefact
