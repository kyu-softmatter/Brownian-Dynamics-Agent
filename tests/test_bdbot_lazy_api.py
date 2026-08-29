"""`bdbot`'s lazy package API (PEP 562) exposes exactly what it used to.

Why this file exists: `bdbot/__init__.py` was changed on 2026-08-29 from eagerly
importing every submodule to resolving them through `__getattr__`, so that
`simbot.units` could import `bdbot.constants` without paying for pint
(measured: `import bdbot` 0.19 s -> 0.01 s). A lazy `__init__` fails **silently** --
a name missing from `_SUBMODULES` or `_NAMES` does not error at import time, it
errors the first time somebody reaches for it, possibly in a case script an hour
into a run.

`bdbot/__init__.py` names this file in a comment. That reference is only worth
anything if the file exists and actually fails when the lists drift, which is what
`test_submodule_list_matches_disk` checks.
"""
from __future__ import annotations

import importlib
from pathlib import Path

import pytest

import bdbot

ROOT = Path(bdbot.__file__).parent

# The import styles that appear in cases/, verify/, tools/ and tests/, collected
# by grepping the tree. If a case adds a new style, add it here.
HISTORICAL_IMPORTS = [
    "from bdbot import Q, checks as C, materials as M, metrics as MET, report as R",
    "from bdbot import Q, checks as C, materials as M, scales as SC, sim as SIM, stats as ST",
    "from bdbot import Q, checks as C, metrics as MET, physical as PH, report as R",
    "from bdbot import Q, checks as C, physical as P, report as R, scales as SC",
    "from bdbot import Q, nondim as ND, scales as SC",
    "from bdbot import Q, physical as PH",
    "from bdbot import health as H",
    "from bdbot import intake as I",
    "from bdbot import interactions as X",
    "from bdbot import lockin as LI, run as RUN, sim as SIM",
    "from bdbot import nondim as ND, run as RUN, runid as RID, scales as SC",
    "from bdbot import sim as SIM, stats as ST, traps as TR",
    "from bdbot.checks import Check",
    "from bdbot.lockin import agg, k_star, lockin_blocks",
    "from bdbot.pairpot import HEX_NN, R_WCA, U2_star, U_star, approach_distance",
    "from bdbot.provenance import load_node",
    "from bdbot.units import Q",
]


def test_submodule_list_matches_disk():
    """`_SUBMODULES` must list every module in the package, and nothing else.

    Both directions matter. A module missing from the list is unreachable as
    `bdbot.<name>`; a name in the list with no file behind it raises only when
    touched.
    """
    on_disk = {p.stem for p in ROOT.glob("*.py") if p.stem != "__init__"}
    listed = set(bdbot._SUBMODULES)
    assert listed == on_disk, (
        f"bdbot._SUBMODULES has drifted from the package contents.\n"
        f"  on disk but not listed: {sorted(on_disk - listed)}\n"
        f"  listed but not on disk: {sorted(listed - on_disk)}")


@pytest.mark.parametrize("name", sorted(bdbot._SUBMODULES))
def test_every_listed_submodule_imports(name):
    mod = getattr(bdbot, name)
    assert mod is importlib.import_module(f"bdbot.{name}")


@pytest.mark.parametrize("name", sorted(bdbot._NAMES))
def test_every_reexported_name_resolves(name):
    """A re-exported name has to come from the module `_NAMES` claims it does."""
    obj = getattr(bdbot, name)
    home = importlib.import_module(f"bdbot.{bdbot._NAMES[name]}")
    assert obj is getattr(home, name)


def test_all_is_fully_resolvable():
    missing = [n for n in bdbot.__all__ if not hasattr(bdbot, n)]
    assert not missing, f"__all__ advertises unresolvable names: {missing}"


@pytest.mark.parametrize("stmt", HISTORICAL_IMPORTS)
def test_historical_import_style_still_works(stmt):
    """Every import form the tree actually uses. The lazy switch must be invisible."""
    exec(compile(stmt, "<historical>", "exec"), {})


def test_unknown_attribute_still_raises_attribute_error():
    """`__getattr__` must not turn a typo into something truthy."""
    with pytest.raises(AttributeError):
        bdbot.no_such_module


def test_the_check_can_fail(monkeypatch):
    """Deliberately break it and see (CLAUDE.md working practice).

    A lazy-API test that cannot fail is the unwired-checker failure this project
    keeps documenting, so prove the drift detector actually fires.
    """
    monkeypatch.setattr(bdbot, "_SUBMODULES",
                        tuple(n for n in bdbot._SUBMODULES if n != "constants"))
    with pytest.raises(AssertionError, match="constants"):
        test_submodule_list_matches_disk()


def test_lazy_means_constants_does_not_pull_pint():
    """The whole point of the change, as a test rather than a claim.

    `import bdbot.constants` in a fresh interpreter must not import pint. Measured
    in-process by checking what a subprocess ends up with in `sys.modules`.
    """
    import subprocess
    import sys

    out = subprocess.run(
        [sys.executable, "-c",
         "import bdbot.constants, sys; "
         "print('pint' in sys.modules, 'numpy' in sys.modules)"],
        capture_output=True, text=True, cwd=ROOT.parent, check=True).stdout.strip()
    assert out == "False False", (
        f"importing bdbot.constants pulled in heavy dependencies ({out}) -- the "
        f"lazy __init__ has regressed, and simbot.units pays for it on every import")
