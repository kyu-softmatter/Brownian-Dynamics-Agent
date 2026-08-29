"""The three pipeline-level merges, pinned by behaviour rather than by syntax.

`docs/00-merge-decisions.md` section 5 called these the "two engines" seam:

    ① the SI data model     bdbot.physical.PhysicalSystem vs simbot.spec.SystemSpec
    ② the runners           bdbot.run.StepGuard          vs simbot.run.run_*
    ③ the CLIs              bdbot/cli.py                 vs cli.py

They could not be collapsed onto one implementation, because each half is the
input side of a **content-addressed archive**: `run_id` names 278 spec files and
263 run directories, and `SEALED.sha256` pins the S2 prediction documents.
Unifying the serialisation would rename one archive or break the other's seals.

So the merge is: **one definition of the physics, two serialisations, and a
verified crossing.** That claim is only worth anything if the crossing is checked,
which is what this file does.

⚠ The check here is deliberately **behavioural, not syntactic.** A peer session
  offered a token-stream comparator that proves "only comments and strings
  changed"; that is the wrong claim for a data-model crossing -- it is not
  invariant to moving code between files, and it would accept a silently changed
  provenance string. What has power is: build the same physical system down both
  paths and demand the derived SI scales come out **bit-identical.**
"""
from __future__ import annotations

import glob
import math

import pytest

from bdbot import constants as BC
from bdbot import health as BH
from bdbot import materials as BM
from bdbot import physical as PH
from bdbot.provenance import Provenanced
from bdbot.units import Q as pintQ
from simbot import bridge as BR
from simbot.spec import PROVENANCE_KINDS, Quantity

SYSTEM_YAMLS = sorted(glob.glob("intake/*/system.yaml"))


def test_there_are_system_yamls_to_check():
    """Guards against a silently empty suite -- the failure family in docs/05."""
    assert len(SYSTEM_YAMLS) >= 8, (
        f"found only {SYSTEM_YAMLS}; this file's central check would be measuring "
        f"nothing")


# ═══════════════════════════════════════════════════════════════════════════
# ① the data model -- the crossing
# ═══════════════════════════════════════════════════════════════════════════
def test_every_provenance_kind_has_a_tier():
    """A kind with no tier must **raise**, never fall through to a default.

    Defaulting would launder an unrecognised provenance into a confident one,
    which is the "never invent a value" failure (CLAUDE.md rule 3).
    """
    for kind in PROVENANCE_KINDS:
        tier = BR.tier_of(kind)
        assert tier in (0, 1, 2, 3), f"{kind} -> {tier}"
    with pytest.raises(ValueError, match="not a registered provenance kind"):
        BR.tier_of("wishful_thinking")


def test_the_reverse_map_is_a_valid_kind_and_conservative():
    for tier in (0, 1, 2, 3):
        kind = BR.provenance_of(tier)
        assert kind in PROVENANCE_KINDS
        assert BR.tier_of(kind) == tier, "the representative must sit at its own tier"
    #  the representative must never overstate the source
    assert BR.provenance_of(0) != "measured", "tier 0 does not imply a measurement"
    assert BR.provenance_of(1) != "from_paper", "tier 1 does not name a citable paper"


def test_the_roundtrip_is_lossy_and_says_so():
    """⚠ Ten kinds collapse onto four tiers, so the round trip **cannot** be
    faithful. Pinning the loss explicitly is the point: a future reader must not
    assume a value carried across the seam kept its categorical provenance.
    """
    assert BR.LOSSY_UNDER_ROUNDTRIP, "if this is empty the map became a bijection"
    for kind in BR.LOSSY_UNDER_ROUNDTRIP:
        assert BR.provenance_of(BR.tier_of(kind)) != kind
    #  and the survivors really do survive
    survivors = set(BR.PROVENANCE_TO_TIER) - set(BR.LOSSY_UNDER_ROUNDTRIP)
    for kind in survivors:
        assert BR.provenance_of(BR.tier_of(kind)) == kind


def test_leaf_conversion_lands_in_SI_BASE_units():
    """★ The bug this test exists for, pinned.

    `bdbot` carries the unit with the value (`0.851 mPa*s` is complete); `simbot`
    names its fields `*_si` and `Quantity.si` does **no conversion**. Handing over
    `(0.851, "mPa*s")` therefore reads back as `0.851 Pa*s` -- gamma 1000x too
    large, D 1000x too small, every timescale 1000x off. It broke 5 of the 8 real
    `system.yaml` files when `simbot.bridge` was first written, and
    `derived_agreement()` is what caught it.
    """
    eta = Provenanced(value=pintQ(0.851, "mPa*s"), source="handbook", tier=0)
    q = BR.provenanced_to_quantity(eta)
    assert q.si == pytest.approx(8.51e-4, rel=1e-15)
    assert q.si != 0.851
    #  a micrometre diameter must arrive in metres
    d = Provenanced(value=pintQ(5.0, "um"), source="sketch", tier=0)
    assert BR.provenanced_to_quantity(d).si == pytest.approx(5e-6, rel=1e-15)
    #  and the original unit is not lost -- it moves into the basis
    assert "mPa" in BR.provenanced_to_quantity(eta).basis


def test_a_sweep_leaf_is_refused_not_silently_indexed():
    """★ A structural asymmetry between the models, found by measurement.

    `bdbot`'s L2 allows a **parameter sweep inside a leaf** --
    `chain-bend-2d-dlvo` has `particle.count = [5, 9, 15, 25]`. A `SystemSpec`
    describes exactly one system. Taking element 0 silently would be choosing a
    physical system by accident. 3 of the 8 cases hit this.
    """
    sweeps = []
    for f in SYSTEM_YAMLS:
        ps = PH.load(f)
        try:
            BR.physical_to_systemspec(ps)
        except ValueError as e:
            if "sweep" in str(e):
                sweeps.append(f)
                assert "sweep_index=" in str(e), "the refusal must say how to resolve it"
    assert len(sweeps) >= 3, f"expected >=3 sweep cases, got {sweeps}"
    #  and with an explicit index it goes through
    ps = PH.load(sweeps[0])
    spec = BR.physical_to_systemspec(ps, sweep_index=0)
    assert spec.primary.n_simulated.si > 0
    with pytest.raises(ValueError, match="outside"):
        BR.physical_to_systemspec(ps, sweep_index=99)


@pytest.mark.parametrize("path", SYSTEM_YAMLS,
                         ids=[p.split("/")[1] for p in SYSTEM_YAMLS])
def test_both_paths_derive_bit_identical_SI_scales(path):
    """★ **The central check of the whole merge.**

    One `system.yaml`, derived down both halves' code paths, compared with
    `rtol = 0` -- bit-identical, not `approx`. Anything else means the two halves
    disagree about the physics, which is the thing the merge claims they no longer
    do.
    """
    ps = PH.load(path)
    try:
        report = BR.derived_agreement(ps)
    except ValueError as e:
        if "sweep" not in str(e):
            raise
        report = BR.derived_agreement(ps, sweep_index=0)
    if report["all_agree"] is None:
        pytest.skip(report["why"])
    bad = {k: v for k, v in report["quantities"].items() if not v["identical"]}
    assert report["all_agree"], f"not bit-identical: {bad}"
    assert len(report["quantities"]) >= 4, "kT, gamma, D_t and tau_B at minimum"


# ═══════════════════════════════════════════════════════════════════════════
# the shared kernels the crossing depends on
# ═══════════════════════════════════════════════════════════════════════════
def test_sphere_mass_is_one_expression_in_both_halves():
    """They were **exactly 1 ULP** apart: `rho*(pi/6)*d^3` vs `rho*(4/3)*pi*a^3`.

    Equal on paper, unequal in floating point, so the two halves could never
    compare equal. Safe to unify because mass reaches no hash -- it lands at
    `.system.derived_scales.tau_p` and `derived_scales` is in
    `bdbot.runid.DOC_KEYS`, so `physics_only()` strips it.
    """
    from simbot.spec import Species
    rho, d = 1180.0, 1.47e-6
    kernel = BC.sphere_mass_si(rho, d)
    bdbot_side = float(BM.sphere_mass(pintQ(rho, "kg/m^3"),
                                      pintQ(d, "m")).to("kg").magnitude)
    simbot_side = Species(
        name="x",
        n_simulated=Quantity(value=1, provenance="assumed", basis="test"),
        radius_si=Quantity(value=d / 2, unit="m", provenance="assumed", basis="test"),
        density_si=Quantity(value=rho, unit="kg/m^3", provenance="assumed",
                            basis="test")).mass_si()
    assert bdbot_side == kernel, "bdbot must call the kernel, not its own expression"
    assert simbot_side == kernel, "simbot must call the kernel, not its own expression"
    #  the old simbot expression, kept here to show the 1 ULP is real and not folklore
    old_simbot = rho * (4.0 / 3.0) * math.pi * (d / 2) ** 3
    assert old_simbot != kernel
    assert abs(old_simbot - kernel) == pytest.approx(math.ulp(kernel), rel=0.5)


# ═══════════════════════════════════════════════════════════════════════════
# ② the runners -- one bound, two reactions
# ═══════════════════════════════════════════════════════════════════════════
def test_the_step_displacement_bound_has_one_definition():
    from bdbot.run import StepGuard
    assert StepGuard(dt_star=1e-4, n_particles=1).step_disp_max == BH.STEP_DISP_MAX
    ok, why = BH.step_displacement_verdict(BH.STEP_DISP_MAX)
    assert ok and why == ""
    ok, why = BH.step_displacement_verdict(BH.STEP_DISP_MAX * 1.001)
    assert not ok and "Reduce dt" in why


def test_the_two_runners_react_differently_on_purpose():
    """⚠ `bdbot.run` **raises**, `simbot.run` **reports**. That is not an
    inconsistency to fix: `bdbot.run` drives one run and should abort it, while
    `simbot.run` drives a batch over seeds and one diverging seed must not discard
    the other seeds' work. What they must share is the **bound**, not the reaction.
    """
    import inspect

    from bdbot import run as BR_run
    from simbot import run as SR

    assert "step_displacement_verdict" in inspect.getsource(BR_run.StepGuard.check)
    assert "raise Diverged" in inspect.getsource(BR_run.StepGuard.check)
    trap_src = inspect.getsource(SR.run_trap)
    assert "step_displacement_verdict" in trap_src
    assert "raise" not in trap_src.split("step_displacement_verdict")[1][:400]
    assert "guard_fail.append" in trap_src


def test_simbot_runners_now_judge_the_quantity_they_measured():
    """The gap this closed: both simbot runners **measured** the force displacement
    and never compared it to a threshold. `check_finite` passes on a box escape --
    `.claude/rules/overdamped-stability.md`: *"the symptom was not NaN but a quiet
    box escape, so the log looked normal."*
    """
    import inspect

    from simbot import run as SR
    for fn in (SR.run_trap, SR.run_soft2d):
        src = inspect.getsource(fn)
        assert "force_displacement_ok" in src, f"{fn.__name__} still does not judge it"


# ═══════════════════════════════════════════════════════════════════════════
# ③ the CLIs -- one entry point
# ═══════════════════════════════════════════════════════════════════════════
def test_one_entry_point_reaches_both_halves():
    from bdbot.cli import build_parser
    actions = [a for a in build_parser()._actions if hasattr(a, "choices") and a.choices]
    cmds = set()
    for a in actions:
        cmds |= set(a.choices or {})
    for expected in ("status", "intake", "nondim", "health", "run", "pipeline"):
        assert expected in cmds, f"`{expected}` missing from the merged surface: {sorted(cmds)}"


def test_the_root_shim_reexports_the_whole_surface():
    """A shim covering only today's call sites breaks on the next one."""
    import cli

    from simbot import cli as impl
    for name in dir(impl):
        if name.startswith("__"):
            continue
        assert hasattr(cli, name), f"cli.{name} does not resolve through the shim"
    with pytest.raises(AttributeError, match="simbot.cli"):
        cli.definitely_not_a_command


def test_bdbot_cli_does_not_pay_for_the_other_half():
    """`bdbot.cli` must not import matplotlib. The delegation is inside the handler.

    This property is stated in `bdbot/__init__.py` and it is why `pipeline` uses a
    function-local import.
    """
    import subprocess
    import sys
    out = subprocess.run(
        [sys.executable, "-c",
         "import sys; import bdbot.cli; "
         "print('matplotlib' in sys.modules, 'simbot.viz' in sys.modules)"],
        capture_output=True, text=True, check=True)
    assert out.stdout.strip() == "False False", out.stdout
