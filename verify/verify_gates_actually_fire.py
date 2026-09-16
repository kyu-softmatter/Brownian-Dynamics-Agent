"""Deliberately trip every load-bearing gate whose error path has never executed.

`verify_error_paths_fire.py` found 88 of bdbot's 154 error sites never executed
under the 1440-test suite, 38 of them in gate modules. This one trips the ones
whose failure is load-bearing -- a run aborted, a rule refused, a spec rejected --
and reports whether each fires, with what message.

CLAUDE.md: "When you build a checker, deliberately break it and see. 'Silently
passing' and 'not checking' are different things." That is the practice; this is
the procedure form, which does not need anyone to remember.

    $PY verify/verify_gates_actually_fire.py
"""
from __future__ import annotations
import math, sys, traceback
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

CASES = []
def case(what, why, want, frag):
    """`want` is the exception type the gate must raise; `frag` a distinctive
    fragment of its message.

    ⚠ The first version of this harness counted ANY exception as the gate firing.
    Seven of twelve "FIRES" were `TypeError: unexpected keyword argument` from my
    own wrong call signatures -- a checker passing for the wrong reason, in the
    script written to find checkers that pass for the wrong reason. Requiring the
    type and the message is what makes the result mean anything.
    """
    def deco(fn):
        CASES.append((what, why, want, frag, fn)); return fn
    return deco


# ---- health: the runtime guards that abort a diverging run ------------------
@case("health NUM_NONFINITE (positions)", "a run with a NaN position must abort, not finish",
      RuntimeError, "[NUM_NONFINITE]")
def _():
    from bdbot import health as H
    g = H.Guard(box_L=10.0)
    pos = np.zeros((4, 3)); pos[2, 0] = float("nan")
    g.check(timestep=100, positions=pos, pe=0.0)

@case("health NUM_DIVERGE (positions)", "a particle far outside the box must abort",
      RuntimeError, "[NUM_DIVERGE]")
def _():
    from bdbot import health as H
    g = H.Guard(box_L=10.0)
    pos = np.zeros((4, 3)); pos[1, 0] = 1e30
    g.check(timestep=100, positions=pos, pe=0.0)

@case("health NUM_NONFINITE (PE)", "a non-finite potential energy must abort",
      RuntimeError, "step 100: PE=nan")
def _():
    from bdbot import health as H
    g = H.Guard(box_L=10.0)
    g.check(timestep=100, positions=np.zeros((4, 3)), pe=float("nan"))

@case("health fdt_residual: too few points", "4 points minimum, or the slope is noise",
      ValueError, "need >= 4 points")
def _():
    from bdbot import health as H
    H.fdt_residual(np.array([1.0, 2.0]), np.array([1.0, 2.0]), D_expected=1.0, dim=2)

@case("health fdt_residual: D_expected <= 0", "a non-positive D is not a diffusivity",
      ValueError, "D_expected must be > 0")
def _():
    from bdbot import health as H
    H.fdt_residual(np.arange(1.0, 9.0), np.arange(1.0, 9.0), D_expected=-1.0, dim=2)


# ---- metrics: rule 7' -------------------------------------------------------
@case("metrics: unknown role", "rule 7' -- only the three roles have meaning",
      ValueError, "role must be one of")
def _():
    from bdbot import metrics as M
    M.observable("x", 1.0, predicted=1.0, tol_pct=5.0, role="looks_right")

@case("metrics: composite + implementation_check without a derivation",
      "rule 7' -- a composite check derived from nothing is a hypothesis wearing a bug's clothes",
      ValueError, "requires a derivation")
def _():
    from bdbot import metrics as M
    M.observable("x", 1.0, predicted=1.0, tol_pct=5.0,
                 role="implementation_check", scope="composite")


# ---- scales / params: rule 3 ------------------------------------------------
@case("scales.declare_absent without a reason", "rule 3 -- absence stated, not omitted",
      ValueError, "requires a reason")
def _():
    from bdbot import scales as SC
    lg = SC.ScaleLedger()
    lg.declare_absent("times", "")

@case("params.not_applicable without a reason", "rule 10 -- 'if applicable' means stated",
      ValueError, "requires a reason")
def _():
    from bdbot import params as PM
    PM.Manifest(case="x").not_applicable("interaction", "  ")

@case("params.add with an unknown category", "a typo'd category must not vanish silently",
      ValueError, "category must be one of")
def _():
    from bdbot import params as PM
    PM.Manifest(case="x").add("geomerty", "d", 1.0, "um", "typo", 0)

@case("params.unknown without a supplier", "BLOCKED must name who would supply it",
      ValueError, "say who or what would supply it")
def _():
    from bdbot import params as PM
    PM.Manifest(case="x").unknown("eta", "")


# ---- nondim: rule 2 ---------------------------------------------------------
@case("nondim: schema mismatch on a stored spec", "rule 2 -- a spec from another schema is not runnable",
      ValueError, "schema differs")
def _():
    from bdbot import nondim as ND
    import json, tempfile
    f = Path(tempfile.mkstemp(suffix=".json")[1])
    f.write_text(json.dumps({"schema": "not.the.schema/0.0"}))
    ND.load(f)


def main():
    print("=" * 96)
    print("Deliberately tripping the load-bearing gates whose error paths never executed")
    print("=" * 96)
    fired = dead = 0
    wrong = 0
    for what, why, wantT, frag, fn in CASES:
        try:
            fn()
        except wantT as e:
            msg = str(e).replace("\n", " ")
            if frag in msg:
                fired += 1
                print(f"  FIRES   {what}\n          {wantT.__name__}: {msg[:76]}")
            else:
                wrong += 1
                print(f"  ⚠ WRONG {what}\n          right type, wrong message -- wanted {frag!r}, got: {msg[:60]}")
        except Exception as e:
            wrong += 1
            print(f"  ⚠ WRONG {what}\n          {type(e).__name__} is not {wantT.__name__}: {str(e)[:66]}\n"
                  f"          The harness, not the gate. Fix the call, then re-run.")
        else:
            dead += 1
            print(f"  DEAD    {what}\n          returned without raising. {why}")
    print("=" * 96)
    print(f"  {fired} fired correctly, {wrong} harness errors, {dead} DEAD")
    dead += wrong
    if dead:
        print("  A gate that does not refuse when handed the thing it exists to refuse is the")
        print("  unwired-checker shape. Fix before anything else.")
    else:
        print("  Every gate tripped above is now something someone has watched run.")
    return 1 if dead else 0


if __name__ == "__main__":
    sys.exit(main())
