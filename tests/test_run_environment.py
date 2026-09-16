"""The run environment and the permission rules, checked by executing them.

★ Why this file exists. Two things this repository documents as load-bearing were
measured on 2026-09-15 to be pointing at nothing:

  * `CLAUDE.md`'s `PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python`
    **does not exist** -- not on the development machine (miniforge, not a
    Homebrew miniconda), not on CI (micromamba under `/home/runner`), not
    anywhere. It was the documented instruction in 34 files, and
    `.claude/settings.json` allow-listed the same dead path, so the allow rules
    matched nothing either and every command prompted.
  * three of the four `Edit(...)` deny rules protecting sealed predictions matched
    **zero files**: they named `./runs/**/02_prediction.md` and
    `./runs/**/01_intake.md`, the simbot-era layout, which lives under
    `runs_s1s8/`. Meanwhile the 24 documents actually under seal today --
    `prediction.yaml` and `analysis_plan.yaml` in the 12 `runs/` campaign
    directories -- were **not protected by any rule**.

Both are the same failure: a rule written once, never executed, and therefore
never wrong out loud. So this file executes them.

⚠ A deny rule that matches zero files is NOT itself a defect -- a rule for a
document that does not exist yet is forward protection. The property that matters
is **coverage**: every document named in a seal must be matched by some rule.
That is asserted from the seals themselves, not from a hand-written list, because
a hand-written list is what drifted.
"""
from __future__ import annotations

import json
import pathlib
import re
import subprocess

import pytest

ROOT = pathlib.Path(__file__).resolve().parent.parent
PY_WRAPPER = ROOT / "bin" / "py"
SETTINGS = ROOT / ".claude" / "settings.json"

#: Files that quote the dead path on purpose and must keep it verbatim.
#: `docs/history/` holds copies of the predecessor repositories' own documents,
#: where the path was correct; `runs_s1s8/` holds run artefacts from that era.
#: Rewriting either would falsify a historical record. `CLAUDE.md` and `bin/py`
#: quote it while explaining that it is dead.
QUOTES_THE_DEAD_PATH_ON_PURPOSE = (
    "docs/history/",
    "runs_s1s8/",
    "CLAUDE.md",
    "bin/py",
    "tests/test_run_environment.py",
)


def _tracked_text_files() -> list[pathlib.Path]:
    out = subprocess.run(["git", "ls-files", "-z"], cwd=ROOT,
                         capture_output=True, check=True).stdout
    names = [f.decode("utf-8", "surrogateescape") for f in out.split(b"\0") if f]
    skip_suffix = {".jpeg", ".jpg", ".png", ".gsd", ".npz", ".pdf", ".h5"}
    return [ROOT / n for n in names if pathlib.PurePath(n).suffix not in skip_suffix]


def _rel(p: pathlib.Path) -> str:
    return p.relative_to(ROOT).as_posix()


# ── the interpreter the documentation tells you to use ────────────────────
def test_the_wrapper_exists_and_is_executable():
    assert PY_WRAPPER.exists(), "bin/py is gone; CLAUDE.md tells every reader to use it"
    assert PY_WRAPPER.stat().st_mode & 0o111, "bin/py is not executable"


def test_the_wrapper_resolves_to_this_environment():
    """★ The whole point. `bin/py` must land on `simulation_bot`, not on a system
    python -- a system python has no `hoomd`, so the failure would look like an
    import error rather than an environment error."""
    r = subprocess.run([str(PY_WRAPPER), "-c", "import sys; print(sys.prefix)"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    assert r.stdout.strip().endswith("/simulation_bot"), r.stdout


def test_the_wrapper_reaches_the_engine_version_the_docs_claim():
    r = subprocess.run([str(PY_WRAPPER), "-c",
                        "import hoomd; print(hoomd.version.version)"],
                       cwd=ROOT, capture_output=True, text=True)
    assert r.returncode == 0, r.stderr
    got = r.stdout.strip()
    assert got in ROOT.joinpath("CLAUDE.md").read_text(), (
        f"the wrapper reaches hoomd {got}, which CLAUDE.md does not mention")


def test_the_wrapper_refuses_a_wrong_override_instead_of_using_it():
    """Guard on the guard: `BD_PY` pointing somewhere else must fail loudly. A
    wrapper that silently accepts the wrong interpreter is worse than no
    wrapper, because the error surfaces later and elsewhere."""
    import sys
    r = subprocess.run([str(PY_WRAPPER), "-c", "print('should not run')"],
                       cwd=ROOT, capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "HOME": "/nonexistent",
                            "BD_PY": "/usr/bin/python3"})
    assert r.returncode != 0, r.stdout
    assert "not the simulation_bot environment" in r.stderr, r.stderr
    assert "should not run" not in r.stdout


def test_the_wrapper_fails_loudly_when_there_is_no_environment():
    r = subprocess.run([str(PY_WRAPPER), "-c", "print(1)"], cwd=ROOT,
                       capture_output=True, text=True,
                       env={"PATH": "/usr/bin:/bin", "HOME": "/nonexistent"})
    assert r.returncode != 0
    assert "could not find" in r.stderr and "Looked at:" in r.stderr, r.stderr


def test_no_tracked_file_names_an_interpreter_path_that_does_not_exist():
    """★ The general form of the defect: an absolute interpreter path written into
    documentation ages into a lie, silently, because nothing runs it.

    Anything that still quotes the dead path must be a historical record or the
    explanation itself -- and that exemption list is explicit, so adding a file to
    it is a visible decision.
    """
    ABS_PY = re.compile(r"(/(?:opt|usr|home|Users)/[\w./+-]*/bin/(?:python[\d.]*|pytest))")
    offenders = []
    for p in _tracked_text_files():
        rel = _rel(p)
        if rel.startswith(QUOTES_THE_DEAD_PATH_ON_PURPOSE) or rel in QUOTES_THE_DEAD_PATH_ON_PURPOSE:
            continue
        try:
            text = p.read_text()
        except (UnicodeDecodeError, IsADirectoryError):
            continue
        for hit in set(ABS_PY.findall(text)):
            if not pathlib.Path(hit).exists():
                offenders.append(f"{rel}: {hit}")
    assert not offenders, (
        "these name an interpreter that does not exist:\n  " + "\n  ".join(offenders))


def test_the_allow_list_points_at_the_wrapper_and_not_at_an_absolute_path():
    allow = json.loads(SETTINGS.read_text())["permissions"]["allow"]
    py_rules = [a for a in allow if "py" in a and a.startswith("Bash(")]
    assert any("bin/py" in a for a in py_rules), py_rules
    for a in allow:
        m = re.fullmatch(r"Bash\((/[^:*)]+)[:*].*\)", a)
        if m:
            assert pathlib.Path(m.group(1)).exists(), (
                f"allow rule names a path that does not exist: {a}")


# ── the deny rules that protect a sealed prediction ───────────────────────
def _deny_path_rules() -> list[str]:
    deny = json.loads(SETTINGS.read_text())["permissions"]["deny"]
    out = []
    for r in deny:
        m = re.fullmatch(r"(?:Edit|Write|MultiEdit)\(\./(.*)\)", r)
        if m:
            out.append(m.group(1))
    return out


def _glob_to_regex(pattern: str) -> re.Pattern:
    """`a/**/b` -> matches `a/x/b` and `a/x/y/b`. Only the constructs the deny
    rules actually use are supported; anything else should fail here rather than
    match by accident."""
    assert not set(pattern) & set("[]{}?"), f"unsupported glob: {pattern}"
    parts = pattern.split("/")
    out = []
    for part in parts:
        if part == "**":
            out.append(r"(?:[^/]+/)*")
        else:
            out.append(re.escape(part).replace(r"\*", r"[^/]*") + "/")
    return re.compile("^" + "".join(out)[:-1] + "$")


def test_every_sealed_document_is_covered_by_a_deny_rule():
    """★★ The property, derived from the seals rather than from a list.

    The previous rules named the simbot-era layout under `./runs/**`, where it
    does not live, so three of four matched zero files while the 24 documents
    actually under seal -- `prediction.yaml` and `analysis_plan.yaml` in the 12
    campaign directories -- had no rule at all.
    """
    rules = [_glob_to_regex(p) for p in _deny_path_rules()]
    assert rules, "there are no path-shaped deny rules left"

    sealed: set[str] = set()
    seals = sorted(ROOT.glob("runs*/*/SEALED.sha256"))
    assert len(seals) >= 21, f"only {len(seals)} seals found"
    for seal in seals:
        sealed.add(_rel(seal))
        for line in seal.read_text().splitlines():
            line = line.strip()
            if not line:
                continue
            _, _, rel = line.partition("  ")
            # the recorded path may predate a directory rename; the document
            # verified is the seal's sibling (simbot.io.verify_seal)
            sealed.add(_rel(seal.parent / pathlib.PurePath(rel.strip()).name))

    uncovered = sorted(d for d in sealed if not any(r.match(d) for r in rules))
    assert not uncovered, (
        f"{len(uncovered)} of {len(sealed)} sealed documents are editable:\n  "
        + "\n  ".join(uncovered[:15]))


def test_the_deny_rules_cover_both_sealing_vocabularies():
    """There are two, and the reason is recorded in `bdbot/runcard.py`'s module
    docstring (`bdbot` cannot import `simbot`). A rule set that knows only one of
    them is how this broke."""
    rules = " ".join(_deny_path_rules())
    for doc in ("prediction.yaml", "analysis_plan.yaml",          # bdbot.runcard
                "02_prediction.md", "01_intake.md", "03_spec.yaml",  # simbot.io
                "SEALED.sha256"):
        assert doc in rules, f"no deny rule mentions {doc}"
    for root in ("runs/", "runs_s1s8/"):
        assert root in rules, f"no deny rule covers {root}"


@pytest.mark.parametrize("pattern,path,expect", [
    ("runs/**/prediction.yaml", "runs/x/prediction.yaml", True),
    ("runs/**/prediction.yaml", "runs/a/b/prediction.yaml", True),
    ("runs/**/prediction.yaml", "runs_s1s8/x/prediction.yaml", False),
    ("runs/**/prediction.yaml", "runs/x/prediction.yaml.bak", False),
    ("runs_s1s8/**/SEALED.sha256", "runs_s1s8/y/SEALED.sha256", True),
])
def test_the_glob_translation_used_above_is_itself_right(pattern, path, expect):
    """The coverage test is only as good as this translation, so it is tested
    directly -- otherwise a too-permissive regex would report full coverage."""
    assert bool(_glob_to_regex(pattern).match(path)) is expect
