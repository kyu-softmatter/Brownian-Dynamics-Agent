"""The S8 skeleton — run directories, hashes, sealing.

What this file protects: **does the seal actually work.** If the seal passes
silently, the pipeline runs with nothing at all preventing post-hoc
rationalisation.
"""
from __future__ import annotations

import datetime as dt
import json
from pathlib import Path

import pytest

from simbot import io


# =============================================================================
# hashes
# =============================================================================
def test_sha256_text_matches_known_value():
    # standard sha256 — must equal what an external tool gives
    assert io.sha256_text("") == (
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")


def test_sha256_file_matches_text(tmp_path):
    p = tmp_path / "a.txt"
    p.write_text("hello ηλ\n", encoding="utf-8")
    assert io.sha256_file(p) == io.sha256_text("hello ηλ\n")


def test_payload_hash_is_key_order_independent():
    """Dicts differing only in key order are the same system. Otherwise the cache
    is invalidated every time."""
    assert io.sha256_payload({"a": 1, "b": 2}) == io.sha256_payload({"b": 2, "a": 1})


def test_payload_hash_distinguishes_values():
    assert io.sha256_payload({"a": 1}) != io.sha256_payload({"a": 1.0000001})


def test_code_hash_covers_analysis_subpackage(tmp_path):
    """Leave `analysis/` out and the hash survives an analysis-code change -- a
    hole in the seal."""
    pkg = tmp_path / "pkg"
    (pkg / "analysis").mkdir(parents=True)
    (pkg / "a.py").write_text("x = 1")
    (pkg / "analysis" / "b.py").write_text("y = 1")
    before = io.code_hash(pkg)
    (pkg / "analysis" / "b.py").write_text("y = 2")
    assert io.code_hash(pkg) != before


def test_code_hash_ignores_pycache(tmp_path):
    pkg = tmp_path / "pkg"
    (pkg / "__pycache__").mkdir(parents=True)
    (pkg / "a.py").write_text("x = 1")
    before = io.code_hash(pkg)
    (pkg / "__pycache__" / "junk.py").write_text("noise")
    assert io.code_hash(pkg) == before


def test_code_hash_detects_file_rename(tmp_path):
    """Hash only the contents and a rename is missed -- check the path is in it too."""
    pkg = tmp_path / "pkg"
    pkg.mkdir()
    (pkg / "a.py").write_text("x = 1")
    before = io.code_hash(pkg)
    (pkg / "a.py").rename(pkg / "b.py")
    assert io.code_hash(pkg) != before


def test_env_versions_reports_absent_without_raising():
    v = io.env_versions()
    assert "python" in v
    for name in io.ENV_PACKAGES:
        assert isinstance(v[name], str) and v[name]


# =============================================================================
# provenance — **built in exactly one place**
# =============================================================================
#  ★ 2026-07-29: `report.reproducibility_section` reads `env_hash`, and both runners
#    built their own manifest by hand and left that key out → the report's
#    reproducibility information rendered as a silent blank. The tests below pin
#    that mismatch.
def test_provenance_supplies_every_key_the_report_renders():
    """`provenance()` gives every key the report reads -- a drifted vocabulary
    renders as a blank."""
    p = io.provenance()
    for key in ("code_hash", "git_rev", "git_dirty", "env_hash", "env"):
        assert key in p, f"{key} is not in provenance"
    for name in ("hoomd", "numpy", "scipy", "freud"):
        assert name in p["env"], f"{name} is not in env"


def test_build_manifest_gets_provenance_from_the_single_definition():
    """`build_manifest` uses `provenance()` -- listing it by hand makes it drift."""
    man = io.build_manifest(run_id="r", spec_hash="abc", seed=1)
    prov = io.provenance()
    for key in ("code_hash", "env_hash"):
        assert man[key] == prov[key]
    assert man["env"] == prov["env"]


def test_provenance_records_the_driver_as_a_repo_relative_path():
    """The driver is recorded as a **repo-relative path** -- an absolute one bakes
    in the home directory."""
    p = io.provenance(driver=io.REPO_ROOT / "simbot" / "io.py")
    assert p["driver"] == "simbot/io.py"
    assert p["driver_hash"] == io.file_hash(io.REPO_ROOT / "simbot" / "io.py")
    assert len(p["driver_hash"]) == 12


def test_provenance_hashes_every_driver_when_several_are_given():
    """★ Several driver files means hashing all of them -- catch one and the claim
    becomes false."""
    a = io.REPO_ROOT / "simbot" / "io.py"
    b = io.REPO_ROOT / "simbot" / "run.py"
    p = io.provenance(driver=[a, b])

    assert p["driver"] == ["simbot/io.py", "simbot/run.py"]
    assert set(p["drivers"]) == {"simbot/io.py", "simbot/run.py"}
    assert p["drivers"]["simbot/io.py"] == io.file_hash(a)
    #  the composite hash must change if any single file changes
    other = dict(p["drivers"]); other["simbot/io.py"] = "0" * 12
    assert io.sha256_payload(other)[:12] != p["driver_hash"]
    #  a single path keeps the old shape (reader compatibility)
    single = io.provenance(driver=a)
    assert single["driver"] == "simbot/io.py"
    assert single["driver_hash"] == io.file_hash(a)


def test_provenance_marks_a_missing_driver_instead_of_raising():
    """A missing driver is marked `"?"` -- it neither kills the run with an
    exception nor hides the fact."""
    p = io.provenance(driver=io.REPO_ROOT / "scripts" / "does_not_exist.py")
    assert p["driver_hash"] == "?"


def test_driver_hash_changes_with_content_but_code_hash_does_not(tmp_path):
    """★ `code_hash` covers only `simbot/` -- it cannot catch a driver change.

    That is why `driver_hash` is needed separately. What sets a run's `A` list,
    seeds and analysis windows is the driver in `scripts/`, and if that is not in
    the hash then the artefacts alone cannot answer "what made this run".
    """
    drv = tmp_path / "driver.py"
    drv.write_text("AMPLITUDES = (0.1, 1.0, 10.0)\n")
    before_driver = io.file_hash(drv)
    before_code = io.code_hash()

    drv.write_text("AMPLITUDES = (0.1, 1.0, 100.0)\n")     # a change to the physics
    assert io.file_hash(drv) != before_driver
    assert io.code_hash() == before_code, \
        "simbot did not change but code_hash did"


# =============================================================================
# run_id
# =============================================================================
def test_run_id_is_deterministic_given_date():
    rid = io.new_run_id("Trap 2D 5um!", "c59e93fd24a2", dt.date(2026, 7, 28))
    assert rid == "2026-07-28_trap-2d-5um_c59e93"


def test_run_id_uses_six_hash_chars():
    rid = io.new_run_id("x", "abcdef0123456789", dt.date(2026, 1, 2))
    assert rid.endswith("_abcdef")


def test_slugify_never_returns_empty():
    assert io.slugify("!!!") == "run"


# =============================================================================
# RunDir
# =============================================================================
def test_rundir_rejects_unknown_stage(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    with pytest.raises(KeyError):
        rd.file("no_such_stage")


def test_rundir_roundtrips_json(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write_json("metrics", {"a": 1.5})
    assert rd.read_json("metrics") == {"a": 1.5}


def test_completed_stages_lists_only_existing(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    assert rd.completed_stages() == []
    rd.write("prediction", "# p")
    assert rd.completed_stages() == ["prediction"]


# =============================================================================
# sealing — this section is why the file exists
# =============================================================================
@pytest.fixture
def sealed_run(tmp_path):
    rd = io.RunDir.create(tmp_path, "2026-07-28_x_abc123")
    rd.write("intake", "# S1\n")
    rd.write("prediction", "# S2\nD = 1.00 ± 0.03\n")
    rd.write("spec", "card: x\n")
    io.write_seal(rd)
    return rd


def test_seal_passes_when_untouched(sealed_run):
    v = io.verify_seal(sealed_run)
    assert v.ok and len(v.entries) == 3


def test_seal_catches_edited_prediction(sealed_run):
    """★ If this test does not pass, nothing prevents post-hoc rationalisation."""
    sealed_run.write("prediction",
                     "# S2\nD = 0.42 ± 0.25   (edited after seeing the result)\n")
    v = io.verify_seal(sealed_run)
    assert not v.ok
    assert io.RUN_LAYOUT["prediction"] in v.changed


def test_seal_catches_whitespace_only_edit(sealed_run):
    """It must catch a single space -- grant an exception for a 'trivial edit' and
    it is not a seal."""
    p = sealed_run.file("prediction")
    p.write_text(p.read_text() + " ", encoding="utf-8")
    assert not io.verify_seal(sealed_run).ok


def test_seal_catches_deleted_document(sealed_run):
    sealed_run.file("spec").unlink()
    v = io.verify_seal(sealed_run)
    assert not v.ok and io.RUN_LAYOUT["spec"] in v.missing


def test_seal_catches_missing_seal_file(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("prediction", "# p")
    assert not io.verify_seal(rd).ok


def test_seal_reports_document_added_after_sealing(tmp_path):
    """A prediction written after sealing is not sealed -- it must not report pass."""
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("intake", "# S1\n")
    io.write_seal(rd)                        # no prediction exists at this point
    rd.write("prediction", "# S2 (written after the run)\n")
    v = io.verify_seal(rd)
    assert not v.ok
    assert io.RUN_LAYOUT["prediction"] in v.unsealed


def test_write_seal_refuses_when_nothing_to_seal(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    with pytest.raises(FileNotFoundError):
        io.write_seal(rd)


def test_seal_file_is_shasum_compatible(sealed_run):
    """It must verify under `shasum -a 256 -c`, without our code."""
    lines = sealed_run.read("seal").strip().splitlines()
    for line in lines:
        digest, sep, rel = line.partition("  ")
        assert sep == "  " and len(digest) == 64
        assert int(digest, 16) >= 0            # is it hexadecimal
        assert rel and not rel.startswith(" ")


def test_verify_seal_on_real_first_run():
    """Is the seal of the actual first completed run still valid (regression)."""
    p = io.REPO_ROOT / "runs" / "2026-07-28_trap-2d-5um_2dfb9d"
    if not p.exists():
        pytest.skip("runs/ is gitignored — not present in this checkout")
    v = io.verify_seal(io.RunDir(p))
    assert v.ok, v.summary()


# =============================================================================
# manifest
# =============================================================================
def test_manifest_records_reproducibility_fields(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("prediction", "# p")
    io.write_seal(rd)
    man = io.build_manifest(run_id="r1", spec_hash="abc", seed=[1, 2, 3, 4], rundir=rd)
    for k in ("run_id", "spec_hash", "seed", "code_hash", "git_rev", "env_hash",
              "env", "sealed"):
        assert k in man, k
    json.dumps(man)                            # must be serializable


def test_manifest_dirty_flag_is_tristate():
    """Reporting False when dirty is undecidable falsely claims 'reproducible'."""
    assert io.git_dirty(Path("/")) in (True, False, None)
