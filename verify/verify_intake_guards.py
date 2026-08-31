"""Adversarial test of the front-end checkers -- do they really catch a
deliberately broken spec?

"Passed silently" and "did not check" are different things. This script feeds broken
input to the L0 (observation) and L2 (system) checkers to confirm **each rule
actually fires**.

This test caught a real bug: changing a unit to `furlong^2` made the checker crash
with a pint `DimensionalityError` instead of reporting an error
(bdbot.physical.bulk).

    $PY scratch/verify_intake_guards.py
"""
import copy
import pathlib
import sys
import tempfile

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import intake as I  # noqa: E402
from bdbot import interactions as X  # noqa: E402
from bdbot import physical as P  # noqa: E402
CASE_OK = ROOT / "intake/soft-r3-2d-A-sweep"

results = []


def report(ok, label, detail=""):
    results.append(ok)
    print(f"  {'✓' if ok else '✗'} {label:<44}{detail}")


# ══════════════════════════════════════════════════════════════════════
print("=" * 78)
print("(1) L0 observation checker")
print("=" * 78)
obs_base = yaml.safe_load((CASE_OK / "observation.yaml").read_text())
tmp = pathlib.Path(tempfile.mkdtemp())
(tmp / "sketch_01.jpeg").write_bytes(b"")


def obs_check(mutate, label, want_error=True):
    d = copy.deepcopy(obs_base)
    mutate(d)
    (tmp / "observation.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    try:
        o = I.load(tmp)
    except Exception as e:
        report(False, label, f"CRASH! {type(e).__name__}")
        return
    errs = o.errors
    ok = (len(errs) > 0) == want_error
    report(ok, label, f"{len(errs)} error(s)"
                      + (f"  -> {errs[0].msg[:44]}" if errs else ""))


obs_check(lambda d: None, "unmodified (0 errors is the right answer)", want_error=False)
obs_check(lambda d: d.pop("ambiguities"), "delete the ambiguities key (§8.3)")
obs_check(lambda d: d.pop("unread_regions"), "delete the unread_regions key (§8.3)")
obs_check(lambda d: d.pop("raw_transcription"), "delete the transcription (rule 5)")
obs_check(lambda d: d.__setitem__("raw_transcription", "too short"),
          "transcription too short")
obs_check(lambda d: d["ambiguities"][0].pop("resolution"),
          "delete an ambiguity's resolution key")
obs_check(lambda d: d["missing_required"][0].pop("confidence"),
          "delete tier from an assumed value (rule 3) ★")
obs_check(lambda d: d["stated_quantities"][0].__setitem__("source", ""),
          "blank the source of an explicit number (principle 2)")
obs_check(lambda d: d["ambiguities"].__setitem__(1, dict(d["ambiguities"][1], id="A1")),
          "duplicate ambiguity id")
obs_check(lambda d: d["missing_required"][2].__setitem__("kind", "maybe"),
          "set kind outside the allowed values")

# Does the choice/physical distinction change the verdict?
d = copy.deepcopy(obs_base)
(tmp / "observation.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
o = I.load(tmp)
report(not o.open_missing and len(o.open_choices) == 1,
       "kind:choice does not block",
       f"blocking {len(o.open_missing)} . open choices {len(o.open_choices)}")
for m in d["missing_required"]:
    if m.get("kind") == "choice":
        m.pop("kind")
(tmp / "observation.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
o2 = I.load(tmp)
report(len(o2.open_missing) == 1,
       "with no kind, it falls back conservatively to physical",
       f"{len(o2.open_missing)} blocking")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("(2) L2 system checker")
print("=" * 78)
sys_base = yaml.safe_load((CASE_OK / "system.yaml").read_text())
t2 = pathlib.Path(tempfile.mkdtemp())
(t2 / "observation.yaml").write_text((CASE_OK / "observation.yaml").read_text())
(t2 / "sketch_01.jpeg").write_bytes(b"")


def sys_check(mutate, label, want_error=True):
    d = copy.deepcopy(sys_base)
    mutate(d)
    (t2 / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    try:
        s = P.load(t2)
    except Exception as e:
        report(False, label, f"CRASH! {type(e).__name__}")
        return
    errs = s.errors
    ok = (len(errs) > 0) == want_error
    report(ok, label, f"{len(errs)} error(s)"
                      + (f"  -> {errs[0].msg[:44]}" if errs else ""))


sys_check(lambda d: None, "unmodified (0 errors is the right answer)",
          want_error=False)
sys_check(lambda d: d["derived_scales"].__setitem__("tau_B", {"value": 300.0, "unit": "s"}),
          "change tau_B 242->300 (recomputation mismatch) ★")
sys_check(lambda d: d["particle"]["diameter"].pop("source"),
          "delete d's source (principle 2)")
sys_check(lambda d: d["particle"]["diameter"].__setitem__("unit", "furlong^2"),
          "make the unit dimensionally inconsistent (used to crash) ★")
sys_check(lambda d: d["particle"]["diameter"].__setitem__("unit", "nonsense_unit"),
          "make the unit unparseable")
sys_check(lambda d: d["particle"]["diameter"].__setitem__("tier", 7), "set tier to 7")
sys_check(lambda d: d.pop("derived_from"),
          "delete derived_from (§5.4 invariant) ★")
sys_check(lambda d: d.__setitem__("derived_from", "nope.yaml"),
          "reference a file that does not exist")
sys_check(lambda d: d.pop("medium"), "delete the medium section")

# ★ L0 BLOCKED while L2 is settled.
# Built from **synthetic input** so it does not depend on the real state of any case
# -- once a case is resolved the test would be silently neutered, which is exactly
# what happened when abp-rod was resolved.
obs_blocked = copy.deepcopy(obs_base)
obs_blocked["missing_required"].append({
    "symbol": "made_up_param", "kind": "physical",
    "what": "synthetic unresolved physical gap for this test",
    "assumed_value": None, "resolution": None})
(t2 / "observation.yaml").write_text(yaml.safe_dump(obs_blocked, allow_unicode=True))
(t2 / "system.yaml").write_text(yaml.safe_dump(sys_base, allow_unicode=True))
s = P.load(t2)
# ⚠ This filter must match what bdbot/physical.py ACTUALLY emits, not the `what`
#   text injected above -- the message quotes the SYMBOL, not the description.
#   It read the Korean for "unresolved physical gaps" until 2026-08-29, which
#   stopped matching the moment
#   physical.py was translated to English, so this guard -- the one enforcing rule 3,
#   "never invent a value" -- was failing and nothing surfaced it (this script is not
#   part of the pytest suite). Anchored on the stable English phrase instead.
blk = [i for i in s.errors if "unresolved physical gaps" in i.msg]
report(bool(blk), "L0 BLOCKED while L2 is settled (rule 3) ★",
       blk[0].msg[:44] if blk else "NOT CAUGHT!")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("(3) verdicts on the five real cases")
print("=" * 78)
# 2026-08-04: the three blocked cases were resolved -- user confirmation (abp-rod
# shape and tumbling), paper distillation (chain-bend U_ij), and user confirmation
# plus ★proposals (trap-drag pair and density).
# Since tier-3 proposals are mixed in, READY means "may proceed to L3", NOT approved.
expect = {"trap-2d-5um": "READY", "soft-r3-2d-A-sweep": "READY",
          "abp-rod-2d-run-flip": "READY", "chain-bend-2d-oscill": "READY",
          "trap-drag-2d-hex300": "READY"}
for name, want in expect.items():
    o = I.load(ROOT / "intake" / name)
    got = "FAIL" if o.errors else ("BLOCKED" if o.open_missing else "READY")
    report(got == want, f"{name} -> {got}", f"(expected {want})")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("(4) the interaction recommender -- what does it suggest for an unspecified "
      "U_ij?")
print("=" * 78)
EXPECT_TOP = {
    # case -> (top-ranked key, why it should be that)
    "chain-bend-2d-oscill": ("contact.adhesive_bending",
                             "bead chain + modulus measurement -> tangential contact "
                             "force (the Furst papers)"),
    "trap-drag-2d-hex300": ("pair.soft_power",
                            "hexagonal lattice + structure -> soft repulsion "
                            "(matches the user's confirmation)"),
    "abp-rod-2d-run-flip": ("pair.none",
                            "MSD/MSAD single-particle observables -> no interaction "
                            "needed"),
}
for name, (want, why) in EXPECT_TOP.items():
    o = I.load(ROOT / "intake" / name)
    recs, tags = X.recommend(o)
    got = recs[0][0].key if recs else "(none)"
    report(got == want, f"{name[:26]} → {got}", f"({why})")

# False-positive regression: the word 'rheology' alone must not promote a contact
# model
o = I.load(ROOT / "intake/trap-drag-2d-hex300")
tags = X.infer_tags(o)
report("tangential" not in tags and "gel" not in tags,
       "the word 'microrheology' alone does not attach a contact tag",
       f"(tags: {', '.join(sorted(tags))})")

# Catalogue integrity
bad = [k for k, it in X.CATALOG.items()
       if not it.form or not it.use_when or not it.avoid_when or not it.hoomd]
report(not bad,
       "every catalogue entry has a form, a use, a caveat and a HOOMD mapping",
       f"(empty entries: {bad})" if bad else f"({len(X.CATALOG)} entries)")
needs_ok = [k for k, it in X.CATALOG.items() if it.key != "pair.none" and not it.needs]
report(not needs_ok,
       "every entry except pair.none states what values must be supplied",
       f"(missing: {needs_ok})" if needs_ok else "")
# Is the verification status stated honestly -- nothing marked verified that was
# never actually used?
unused_but_verified = [k for k in ("pair.yukawa", "pair.dlvo", "pair.ao_depletion")
                       if X.CATALOG[k].verified]
report(not unused_but_verified,
       "no never-used interaction is marked as 'verified'",
       f"(false claims: {unused_but_verified})" if unused_but_verified else "")

print()
print("=" * 78)
print(f"{'✓ PASS' if all(results) else '✗ FAIL'} -- "
      f"{sum(results)}/{len(results)} OK")
print("=" * 78)
raise SystemExit(0 if all(results) else 1)
