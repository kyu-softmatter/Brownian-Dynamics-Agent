"""Test the L3 checker by **deliberately breaking it** (CLAUDE.md working practice).

"Silently passing" and "not checking" are different things. `verify_intake_guards.py`
caught a real crash bug this way. What is tested here:

  (1)  does a healthy spec pass (are there no false positives)
  (2)  is a missing required role in the ledger caught
  (3)  is a dimensionless value inconsistent with the ledger caught
       <- the kind that could not be caught before
  (4)  is a group pointing at a nonexistent symbol caught (as an error, not a crash)
  (5)  does changing the physical system change the run_id
       <- regression guard for defect (1)
  (6)  does the inverse transform round-trip (reduced -> physical -> reduced)
  (7)  do save->load preserve the run_id and the groups (L4 uses only this path)
  (8)  does `verify_hash()` catch a hand-edited spec
  (9)  is a reference scale missing from the ledger caught
  (10) is a zero or negative dt* caught
  (11) is a dimensionally inconsistent ratio (length/time) caught when declared
       dimensionless

    $PY scratch/verify_nondim_guards.py
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import Q, nondim as ND, scales as SC  # noqa: E402
from bdbot.checks import Check  # noqa: E402

PASS, FAIL = [], []


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASS if ok else FAIL).append(name)
    print(f"  {'✓' if ok else '✗'} {name}" + (f"\n        {detail}" if detail else ""))


def errs(spec) -> list[str]:
    return [f"{i.where}: {i.msg}" for i in spec.validate() if i.level == "error"]


# ── Assemble one healthy spec by hand (deliberately not depending on any case
#    script) ────────
def good_spec() -> ND.NondimSpec:
    d = Q(5.0, "um").to("m")
    kT = Q(4.141947e-21, "J")
    tau_B = Q(242.1, "s")
    tau_k = Q(4.01, "ms").to("s")
    dt = Q(8.02, "us").to("s")

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("L", Q(160.0, "um").to("m"), "box", role="box")
    lg.add_time("tau_p", Q(3.264, "us").to("s"), "inertia", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_k", tau_k, "trap relaxation", star=True)
    lg.add_time("tau_B", tau_B, "diffusion (reference)")
    lg.add_time("T_obs", Q(8.02, "s"), "observation window", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.ref = SC.thermal_reference(d, kT, tau_B)
    lg.rationale = lg.ref["rationale"]

    groups = [
        ND.Group("dt/tau_k", float((dt / tau_k).to("")), ("times", "dt"), ("times", "tau_k"),
                 "", "integration resolution"),
        ND.Group("L/d", float((lg.get("lengths", "L") / d).to("")),
                 ("lengths", "L"), ("lengths", "d"), "", "box size"),
        ND.Group("phi", 0.35, None, None, "", "packing fraction (an input, not a ratio)"),
    ]
    checks = [Check("integration", "trap resolved dt/tau_k",
                    float((dt / tau_k).to("")), 1e-2, "<=")]
    system = {"label": "guard-test", "dimensions": 2,
              "particle": {"diameter": {"value": 5.0, "unit": "um", "source": "test", "tier": 0}},
              "medium": {"temperature": {"value": 300, "unit": "K", "source": "t", "tier": 0}}}
    # ★ dt* is DERIVED from the ledger. Writing 3.3e-8 by hand disagreed with the
    #   ledger's dt/tau_B (3.3127e-8) and check (10b) caught it -- the checker caught
    #   its own fixture.
    dt_star = float((dt / tau_B).to("dimensionless").magnitude)
    return ND.NondimSpec(
        case="guard-test", system=system, reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"N": 1000}, numerics={"dt_star": dt_star, "n_prod": 1000000})


def main() -> int:
    print("=" * 80)
    print("L3 (bdbot.nondim) adversarial checks")
    print("=" * 80)

    # (1) no false positives
    print("\n[1] a healthy spec")
    s = good_spec()
    e = errs(s)
    check("a healthy spec reports 0 errors", not e, "; ".join(e))

    # (2) a missing required role
    print("\n[2] removing a required role from the ledger")
    for role, sym, cat in (("dt", "dt", "times"), ("observation", "T_obs", "times"),
                           ("box", "L", "lengths"), ("inertia", "tau_p", "times")):
        s = good_spec()
        # Removing a symbol from the ledger must also remove the group that used it,
        # or the test is not a fair one
        getattr(s.ledger, cat).pop(sym)
        s.groups = [g for g in s.groups
                    if sym not in ((g.num or ("", ""))[1], (g.den or ("", ""))[1])]
        e = errs(s)
        check(f"a missing '{role}' role is caught",
              any(f"ledger.{role}" in x for x in e),
              "; ".join(e) or "no error -- NOT caught")

    print("\n[2b] emptying an absent role WITH a reason passes")
    s = good_spec()
    s.ledger.times.pop("tau_p")
    s.ledger.declare_absent("inertia", "test: declared as a system with no inertia")
    check("declare_absent lets it pass", not errs(s), "; ".join(errs(s)))
    try:
        good_spec().ledger.declare_absent("box", "")
        check("emptying it WITHOUT a reason is rejected", False,
              "an empty reason was accepted")
    except ValueError:
        check("emptying it WITHOUT a reason is rejected", True)

    # (3) ⭐️ a group inconsistent with the ledger
    print("\n[3] a dimensionless value inconsistent with the ledger "
          "(the kind that could not be caught before)")
    for factor, label in ((1.41, "41% off -- an a_mean vs a_NN class of mistake"),
                          (1 + 1e-6, "1e-6 off -- a tiny arithmetic slip")):
        s = good_spec()
        s.groups[0].value *= factor
        e = errs(s)
        check(f"scaling dt/tau_k by x{factor} is caught ({label})",
              any("groups.dt/tau_k" in x for x in e),
              "; ".join(e) or "no error -- NOT caught")

    print("\n[3b] a floating-point-level discrepancy (1e-12) must pass "
          "(no false positives)")
    s = good_spec()
    s.groups[0].value *= (1 + 1e-12)
    check("a 1e-12 discrepancy passes", not errs(s), "; ".join(errs(s)))

    # (4) pointing at a nonexistent symbol -> an error, not a crash
    print("\n[4] a group pointing at a symbol absent from the ledger")
    s = good_spec()
    s.groups[0].den = ("times", "tau_nonexistent")
    try:
        e = errs(s)
        check("a nonexistent symbol is reported as an error (not a crash)",
              any("groups.dt/tau_k" in x for x in e), "; ".join(e))
    except Exception as ex:
        check("a nonexistent symbol is reported as an error (not a crash)", False,
              f"died on an exception: {type(ex).__name__}: {ex}")

    # (5) ⭐️ regression guard for defect (1) -- is the physical system in the run_id?
    print("\n[5] does changing the physical system change the run_id "
          "(regression guard for defect (1))")
    a = good_spec()
    b = good_spec()
    b.system = copy.deepcopy(b.system)
    b.system["particle"]["diameter"]["value"] = 0.5
    check("changing d changes the run_id", a.run_id() != b.run_id(),
          f"{a.run_id()} vs {b.run_id()}")

    c = good_spec()
    c.system = copy.deepcopy(c.system)
    c.system["particle"]["diameter"]["source"] = "provenance only, edited"
    check("editing only the source (a documentation field) preserves the run_id",
          a.run_id() == c.run_id(),
          f"{a.run_id()} vs {c.run_id()}")

    dd = good_spec()
    dd.ledger.rationale = "rationale wording only, edited"
    dd.checks = []
    check("editing the ledger rationale or the checks preserves the run_id",
          a.run_id() == dd.run_id(),
          f"{a.run_id()} vs {dd.run_id()}")

    # (6) inverse-transform round trip
    print("\n[6] inverse-transform round trip (reduced -> physical -> reduced)")
    s = good_spec()
    for val, kw, unit in ((0.5, dict(L=2), "um^2"), (3.0, dict(T=1), "s"),
                          (1.83, dict(L=2, T=-1), "um^2/s"), (2.0, dict(E=1), "J")):
        phys = s.physical(val, **kw)
        back = phys
        for kind, expo in (("length", kw.get("L", 0)), ("time", kw.get("T", 0)),
                           ("energy", kw.get("E", 0))):
            back = back / s.reference.si(kind) ** expo
        rel = abs(float(back.to("dimensionless").magnitude) - val) / val
        check(f"{kw} round-trip error < 1e-12  ({phys.to(unit):~.5gP})",
              rel < 1e-12, f"rel={rel:.2e}")

    # (7) save -> load preservation
    print("\n[7] save -> load (the only path L4 uses)")
    s = good_spec()
    tmp = ROOT / "verify" / "_tmp_spec.json"
    s.write(tmp)
    ls = ND.load(tmp)
    check("the run_id is preserved", ls.run_id == s.run_id(),
          f"{ls.run_id} vs {s.run_id()}")
    check("the groups are preserved",
          abs(ls.group("dt/tau_k") - s.groups[0].value) < 1e-15)
    check("the reduced values are computed -- L4 reads L* from here",
          abs(ls.reduced("lengths", "L") - 32.0) < 1e-9,
          f"L/d = {ls.reduced('lengths', 'L')}")
    r1 = float(s.physical(1.83, L=2, T=-1).to("um^2/s").magnitude)
    r2 = float(ls.physical(1.83, L=2, T=-1).to("um^2/s").magnitude)
    check("the loaded spec's inverse transform matches the original",
          abs(r1 - r2) / r1 < 1e-12, f"{r1} vs {r2}")
    ok, want = ls.verify_hash()
    check("hash self-verification passes", ok, f"expected {want}")

    # (8) a hand-edited spec
    print("\n[8] hand-editing a spec (rule 2 -- specs are never written by hand)")
    raw = json.loads(tmp.read_text())
    raw["params"]["N"] = 4000
    tmp.write_text(json.dumps(raw))
    ok, want = ND.load(tmp).verify_hash()
    check("hand-editing params is caught as a hash mismatch", not ok,
          f"expected {want}")

    raw = json.loads(tmp.read_text())
    raw["params"]["N"] = 1000
    raw["ledger_absent"] = {"note": "documentation field only, edited"}
    tmp.write_text(json.dumps(raw))
    ok, _ = ND.load(tmp).verify_hash()
    check("editing a field outside the hash payload still passes", ok)
    tmp.unlink()

    # (9) a reference missing from the ledger
    print("\n[9] a reference scale absent from the ledger")
    s = good_spec()
    s.ledger.times.pop("tau_B")
    s.groups = [g for g in s.groups if "tau_B" not in str(g.num) + str(g.den)]
    e = errs(s)
    check("a reference time absent from the ledger is caught",
          any("reference.time" in x for x in e), "; ".join(e) or "no error")

    # (10) out-of-range run parameters
    print("\n[10] out-of-range run parameters")
    for bad, why in ((0.0, "0"), (-1e-8, "negative")):
        s = good_spec()
        s.numerics["dt_star"] = bad
        e = errs(s)
        check(f"dt* = {why} is caught", any("numerics.dt_star" in x for x in e),
              "; ".join(e))
    s = good_spec()
    s.numerics.pop("n_prod")
    check("a missing n_prod is caught", any("numerics.n_prod" in x for x in errs(s)))

    print("\n[10b] dt* inconsistent with the ledger's dt/tau_B "
          "(HOOMD would run with a different step)")
    s = good_spec()
    true_dt_star = float((s.ledger.get("times", "dt")
                          / s.reference.si("time")).to("dimensionless").magnitude)
    check("a dt* derived from the ledger passes", not errs(s),
          f"ledger dt/tau_B = {true_dt_star:.8e}")
    for factor, why in ((1.01, "1%"), (1 + 1e-6, "1e-6")):
        s2 = good_spec()
        s2.numerics["dt_star"] = true_dt_star * factor
        e = errs(s2)
        check(f"a {why} discrepancy is caught",
              any("numerics.dt_star" in x for x in e), "; ".join(e))

    # (11) a dimensionally inconsistent ratio
    print("\n[11] a ratio that is not dimensionless (length/time)")
    s = good_spec()
    s.groups[0].den = ("lengths", "d")          # dt / d -> not dimensionless
    e = errs(s)
    check("a dimensionally inconsistent ratio is caught",
          any("groups.dt/tau_k" in x for x in e),
          "; ".join(e) or "no error -- NOT caught")

    print("\n" + "=" * 80)
    print(f"passed {len(PASS)} . failed {len(FAIL)}")
    for f in FAIL:
        print(f"   ✗ {f}")
    print("✓ PASS -- the L3 checker really does catch these" if not FAIL
          else "✗ FAIL -- the items above are NOT caught")
    print("=" * 80)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
