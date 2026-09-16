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
