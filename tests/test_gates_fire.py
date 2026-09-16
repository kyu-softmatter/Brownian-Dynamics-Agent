"""Every load-bearing gate must refuse when handed the thing it exists to refuse.

`verify/verify_gates_actually_fire.py` is the enumeration; this wires it into the
suite so it cannot rot. Without this file the procedure that found 88 unfired
error paths would itself be a script nobody runs -- an unwired checker for
unwired checkers, which is where this started.

Clause 3 of the general form: assert the failure's IDENTITY, not that a failure
occurred. Each case pins the exception type AND a distinctive message fragment,
because the first version of the enumeration counted any exception as a pass and
7 of 12 "successes" were TypeError from wrong call signatures.
"""
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT / "verify"))
import verify_gates_actually_fire as G          # noqa: E402


@pytest.mark.parametrize("what,why,want,frag,fn", G.CASES,
                         ids=[c[0] for c in G.CASES])
def test_gate_refuses(what, why, want, frag, fn):
    try:
        fn()
    except want as e:
        assert frag in str(e), (
            f"{what}: right exception type, wrong message. Wanted {frag!r}, got {e!r}. "
            f"A negative test that accepts any message tests nothing.")
    except Exception as e:                       # noqa: BLE001
        pytest.fail(f"{what}: {type(e).__name__} is not {want.__name__} -- this is the "
                    f"harness calling the gate wrongly, not the gate failing. {e}")
    else:
        pytest.fail(f"{what}: returned without raising. {why}")


def test_every_case_pins_a_message():
    """Clause 3, enforced on the enumeration itself."""
    for what, _why, want, frag, _fn in G.CASES:
        assert isinstance(want, type) and issubclass(want, BaseException), what
        assert frag and len(frag) >= 8, f"{what}: message fragment {frag!r} is too weak"


# ---------------------------------------------------------------------------
# Guards on this file itself. Both holes below were live before they were added.
# ---------------------------------------------------------------------------
sys.path.insert(0, str(Path(__file__).resolve().parent))
import _gate_inventory as INV                                   # noqa: E402


def test_the_enumeration_is_not_empty():
    """With `CASES` empty, the parametrised test is SKIPPED and the suite is green.

    Measured 2026-09-16: clearing `CASES` gave "1 passed, 1 skipped", exit 0 --
    the same vacuous pass the bridge found in its own `--selftest`, which printed
    "0 fixtures, all pinned" with every fixture deleted. A check that finds
    nothing to check must fail, not congratulate itself.
    """
    assert len(G.CASES) >= 12, (
        f"the gate enumeration has {len(G.CASES)} cases. If it has shrunk, say why; "
        "if it is empty this whole file tests nothing while reporting green.")
    names = " | ".join(c[0] for c in G.CASES)
    for required in ("health", "metrics", "params", "scales", "nondim"):
        assert required in names, f"no case covers {required}, which gates a run"


def test_no_new_gate_arrives_unobserved():
    """A new error site in a load-bearing module must be classified by a human.

    Not "every site must have a case" -- 52 do not, and writing 52 cases to go
    green would be the cosmetic version. The requirement is that a NEW gate cannot
    appear in neither list, because then nobody has decided whether anyone has
    ever watched it run. That decision is what was missing when `health.Guard`'s
    divergence aborts sat unfired through 1440 tests.
    """
    inv = INV.load()
    known = set(inv["covered"]) | set(inv["known_unobserved"])
    found = set(INV.error_sites())
    new = sorted(found - known)
    gone = sorted(known - found)
    assert not new, (
        f"{len(new)} error site(s) in a load-bearing module are in neither list of "
        f"tests/_gate_inventory.json. Add a case in verify/verify_gates_actually_fire.py "
        f"and put it in `covered`, or put it in `known_unobserved` and say so:\n  "
        + "\n  ".join(new[:6]))
    assert not gone, (
        f"{len(gone)} inventoried error site(s) no longer exist. If a gate was "
        f"removed that is a decision, not a refactor -- drop it from the inventory "
        f"in the same commit:\n  " + "\n  ".join(gone[:6]))
    assert len(inv["covered"]) >= 12, "the covered list shrank"
