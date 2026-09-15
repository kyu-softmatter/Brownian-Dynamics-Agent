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

import json
import pathlib
import subprocess
import sys

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
CASE = ROOT / "cases/soft_r3_2d.py"


def _run(*args, **kw):
    return subprocess.run([sys.executable, str(CASE), *args],
                          capture_output=True, text=True, cwd=ROOT, **kw)


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

def test_the_archived_run_id_does_not_move():
    """The default path must produce the run_id the archive already carries. If
    `init` were written as "rsa" rather than left absent, this is what would
    catch it — nine specs and their run directories would be renamed."""
    r = _run("--A", "100", "--report")
    assert r.returncode == 0, r.stderr[-800:]
    assert "soft-r3-2d-A-sweep__A100__30caa5c9e0" in r.stdout, r.stdout[:300]


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
    r = _run("--A", "34.938", "--N", "144", "--rc-shells", "5", "--smoke",
             "--box", "hex", "--init", "hex", "--eq-frac", "0", "--force")
    assert r.returncode == 0, r.stdout[-1500:] + r.stderr[-1500:]
    line = next(ln for ln in r.stdout.splitlines() if "run_id=" in ln)
    rundir = ROOT / "runs" / line.split("run_id=")[1].strip()
    for name in ("result.txt", "metrics.json", "observables.npz",
                 "observables.png"):
        assert (rundir / name).exists(), f"{name} was not written to {rundir}"
    import numpy as np
    res = np.load(rundir / "observables.npz")
    assert "eq_trace" not in res.files, \
        "an eq_frac=0 run grew an eq_trace -- this test no longer covers the bug"
    assert "psi6_global" in res.files, "the decision statistic is not on disk"


@pytest.mark.slow
def test_the_default_run_still_has_an_equilibration_trace():
    """Guard on the guard: the branch must be a branch. If `eq_trace` vanished
    from every run, the test above would pass for the wrong reason."""
    r = _run("--A", "34.938", "--N", "144", "--rc-shells", "5", "--smoke",
             "--box", "hex", "--init", "hex", "--force")
    assert r.returncode == 0, r.stdout[-1500:]
    line = next(ln for ln in r.stdout.splitlines() if "run_id=" in ln)
    import numpy as np
    res = np.load(ROOT / "runs" / line.split("run_id=")[1].strip() / "observables.npz")
    assert "eq_trace" in res.files, "the default path lost its equilibration trace"
