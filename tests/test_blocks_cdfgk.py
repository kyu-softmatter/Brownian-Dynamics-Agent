"""Blocks C, D, F, G, K — the five modules built 2026-09-02.

Structure per block: does it fire · is it wired · does turning it off change the
verdict · is its scope narrow on purpose. Same shape as `test_a_goal.py`.

Every measured number cited here was taken from the tree at `4b4503a`.
"""
from __future__ import annotations

import copy
import math
from pathlib import Path

import pytest
import yaml

from bdbot import analysisplan as AP
from bdbot import cost as CO
from bdbot import design as DE
from bdbot import goal as GO
from bdbot import groups as GR
from bdbot import runcard as RC

ROOT = Path(__file__).resolve().parent.parent


# ══════════════════════════════════════════════════════════════════════════
# C — the group registry
# ══════════════════════════════════════════════════════════════════════════

def test_the_two_measured_collisions_resolve():
    """★ The reason the registry exists.

    `Pe` and `Pe_drag` were the same group under two names (abp_rod_2d:155 and
    trap_drag_2d:282, both tau_B/tau_v). `k*` and `k_star` likewise.
    """
    assert GR.canonical("Pe") == "Pe_adv"
    assert GR.canonical("Pe_drag") == "Pe_adv"
    assert GR.canonical("k*") == "k_star"
    assert GR.canonical("St") == "St"
    assert GR.canonical("tau_bond/tau_chain_diff") is None   # bespoke is allowed


def test_st_must_name_its_governing_time():
    """★ The other collision. `St` is tau_p/tau_B in 5 cases and tau_p/tau_bond in
    chain_relax_2d_dlvo:312. Both are real Stokes numbers, so the name alone is
    not enough and an unnamed tau_gov is an ERROR, not a warning."""
    bad = [GR.Usage(symbol="St", value=1e-3)]
    msgs = GR.validate_usage(bad)
    assert GR.errors(msgs), msgs
    assert "tau_gov" in msgs[0]

    good = [GR.Usage(symbol="St", value=1e-3, free_choice="tau_B")]
    assert not GR.errors(GR.validate_usage(good))


def test_the_same_group_twice_with_different_free_choices_is_allowed():
    """chain_bend_2d legitimately declares De, De_chain_old and De_trap — three
    relaxation times in one case. That must not be an error; the same free choice
    twice must be."""
    ok = [GR.Usage("De", 1.0, free_choice="tau_max"),
          GR.Usage("De", 2.0, free_choice="tau_trap")]
    assert not GR.errors(GR.validate_usage(ok))
    dup = [GR.Usage("De", 1.0, free_choice="tau_max"),
           GR.Usage("De", 2.0, free_choice="tau_max")]
    assert GR.errors(GR.validate_usage(dup))


def test_an_alias_is_reported_so_aggregation_works():
    msgs = GR.validate_usage([GR.Usage("Pe_drag", 40.0)])
    assert any("alias of 'Pe_adv'" in m for m in msgs), msgs


def test_undefined_groups_ask_for_a_justification():
    """Wi does not exist for a Newtonian solvent or an unforced measurement — the
    passive-microrheology-in-water case is exactly both."""
    msgs = GR.validate_usage([GR.Usage("Wi", 0.0)])
    assert any(m.startswith("[warn]") and "Wi" in m for m in msgs), msgs
    msgs2 = GR.validate_usage([GR.Usage("Wi", 0.0, justification="Boger fluid, "
                                        "lambda_p measured")])
    assert not any(m.startswith("[warn]") for m in msgs2)


def test_overdamped_verdict_reproduces_the_archived_check():
    """`trap-2d-5um`'s spec records
    `{"name": "inertia negligible tau_p/tau_k", "value": 8.139e-04,
      "limit": 0.01, "margin": 12.286}`. The registry must agree with it."""
    tau_p, tau_k = 3.2637e-6, 4.010243022307369e-3
    ok, st, limit, margin = GR.overdamped_verdict(tau_p, tau_k)
    assert ok
    assert math.isclose(st, 8.1385e-4, rel_tol=2e-3), st
    assert limit == 1e-2
    assert math.isclose(margin, 12.286, rel_tol=2e-3), margin


def test_formulas_refuse_bad_input():
    for fn, args in [(GR.st, (1.0, 0.0)), (GR.pe_adv, (1.0, 1.0, 0.0)),
                     (GR.k_star, (1.0, 1.0, 0.0)), (GR.reynolds, (1, 1, 1, 0))]:
        with pytest.raises(ValueError):
            fn(*args)
    with pytest.raises(ValueError):
        GR.de(1.0)                                # neither t_obs nor omega
    with pytest.raises(ValueError):
        GR.de(1.0, t_obs=1.0, omega=1.0)          # both


# ══════════════════════════════════════════════════════════════════════════
# D — the matched / abandoned ledger
# ══════════════════════════════════════════════════════════════════════════

def test_a_group_in_neither_column_is_an_error():
    """★ The silent failure this block exists for: match the subset you can hit,
    report 'matched', and never mention the rest."""
    led = DE.DesignLedger("synthetic", supported=("St", "Re", "Pe_adv"))
    led.match("St", 1e-3, 1e-3)
    msgs = led.validate()
    assert DE.errors(msgs)
    assert led.uncovered() == ["Pe_adv", "Re"]
    assert all(s in " ".join(msgs) for s in ("Re", "Pe_adv"))


def test_covering_every_group_clears_it():
    led = DE.DesignLedger("synthetic", supported=("St", "Re"))
    led.match("St", 1e-3, 1.0005e-3)
    led.abandon("Re", reason="BD has no fluid inertia",
                safe_because="both experiment and simulation are << 1e-2, so "
                             "Stokes drag holds in each")
    assert not DE.errors(led.validate()), led.validate()
    assert led.uncovered() == []


def test_abandoning_without_a_reason_raises():
    led = DE.DesignLedger("synthetic")
    with pytest.raises(ValueError, match="requires a reason"):
        led.abandon("Re", reason="", safe_because="x")
    with pytest.raises(ValueError, match="safe_because"):
        led.abandon("Re", reason="x", safe_because="   ")


def test_a_matched_group_that_does_not_match_is_caught():
    led = DE.DesignLedger("synthetic", supported=("St",))
    led.match("St", 1.0, 1.5)                     # 50 % off, listed as matched
    msgs = led.validate()
    assert DE.errors(msgs) and "50.00 %" in " ".join(msgs)


def test_safe_because_is_cross_checked_against_real_checks():
    """★ `safe_because` is a claim. If it depends on a check that FAILED, the
    justification contradicts the measurement — report it, do not average (A1)."""
    led = DE.DesignLedger("synthetic", supported=("Re",))
    led.abandon("Re", reason="no fluid inertia in BD",
                safe_because="the overdamped check clears",
                depends_on_checks=("inertia negligible",))
    passing = [{"name": "inertia negligible", "ok": True}]
    assert not DE.errors(led.verify_safety(passing))

    failing = [{"name": "inertia negligible", "ok": False}]
    msgs = led.verify_safety(failing)
    assert DE.errors(msgs) and "contradicts the measurement" in msgs[0]

    absent = [{"name": "something else", "ok": True}]
    assert any(m.startswith("[warn]") for m in led.verify_safety(absent))


def test_double_counting_is_caught():
    led = DE.DesignLedger("s", supported=("St",))
    led.match("St", 1.0, 1.0)
    led.abandon("St", reason="r", safe_because="s")
    assert any("both matched and abandoned" in m for m in led.validate())


# ══════════════════════════════════════════════════════════════════════════
# F — the cost gate
# ══════════════════════════════════════════════════════════════════════════

def test_no_budget_is_not_a_pass():
    """A3: 'do not start without an estimate' — and an estimate with nothing to
    compare against is not a gate. An unset budget is INFEASIBLE."""
    est = CO.from_measured("trap-2d-5um", n_steps=1_000_000, n_particles=1000)
    assert est.user_budget_s is None
    assert est.verdict == "INFEASIBLE"
    assert any("A3" in b for b in est.blockers)
    assert CO.blocks(est)


def test_the_archived_run_reproduces_its_own_wall_time():
    """The only measured throughput in the repo: 332.507375 steps/s and 3037.53 s
    wall. ★ `n_steps` must be TOTAL steps: 3037.53 x 332.507375 = 1,010,000 =
    n_eq + n_prod, not n_prod. This test found that ambiguity."""
    n = CO.total_steps(n_eq=10_000, n_prod=1_000_000)
    assert n == 1_010_000
    est = CO.from_measured("trap-2d-5um", n_steps=n, n_particles=1000,
                           user_budget_s=3600, urgency="same-day")
    assert math.isclose(est.estimate_s, 3037.53, rel_tol=1e-4), est.estimate_s
    assert est.verdict == "TIGHT", est.margin      # 3600/3037.5 = 1.185 < 1.5
    assert not est.blockers


def test_using_n_prod_alone_under_predicts_by_one_percent():
    """The mistake, quantified, so the guard is known to be real. A relaxation
    case with n_eq ~ n_prod would be out by ~50 % rather than 1 %."""
    full = CO.from_measured("trap-2d-5um", n_steps=1_010_000, n_particles=1000)
    prod_only = CO.from_measured("trap-2d-5um", n_steps=1_000_000, n_particles=1000)
    assert math.isclose(1 - prod_only.estimate_s / full.estimate_s, 0.0099,
                        abs_tol=1e-4)
    with pytest.raises(ValueError):
        CO.total_steps(-1, 10)


def test_borrowing_another_cases_throughput_is_refused():
    """cli.py:628 already refuses to overwrite the global constant because the
    kernel differs. Same refusal here rather than a silent guess."""
    with pytest.raises(KeyError, match="no measured throughput"):
        CO.from_measured("bead-water-3d", n_steps=1000, n_particles=1)


def test_infeasible_must_carry_what_to_give_up():
    est = CO.from_measured("trap-2d-5um", n_steps=10_000_000, n_particles=1000,
                           user_budget_s=600)
    assert est.verdict == "INFEASIBLE"
    assert any("give_up" in b for b in est.blockers)

    est.give_up = CO.propose_give_up(est)
    assert est.give_up
    assert not any("give_up" in b for b in est.blockers)
    # ★ priced in statistical error, not step count
    assert any("statistical error x" in g for g in est.give_up), est.give_up


def test_error_scaling_is_one_over_sqrt_n():
    assert math.isclose(CO.error_scaling(1000, 250), 2.0)
    assert math.isclose(CO.error_scaling(1000, 1000), 1.0)
    with pytest.raises(ValueError):
        CO.error_scaling(1000, 0)


def test_smoke_route_is_preferred_and_works(tmp_path):
    est = CO.from_smoke("bead-water-3d", n_steps=1_000_000, n_particles=1,
                        smoke_steps=20_000, smoke_wall_s=4.0,
                        user_budget_s=1000, urgency="today")
    assert math.isclose(est.steps_per_second, 5000.0)
    assert math.isclose(est.estimate_s, 200.0)
    assert est.verdict == "FEASIBLE"
    p = est.write(tmp_path)
    assert p.name == "cost.json"
    import json
    d = json.loads(p.read_text())
    assert d["schema"] == CO.SCHEMA and d["verdict"] == "FEASIBLE"


# ══════════════════════════════════════════════════════════════════════════
# G — the run card and the seal
# ══════════════════════════════════════════════════════════════════════════

def _seal_fixture(tmp_path):
    (tmp_path / "prediction.yaml").write_text("k_prime: 0.0\n")
    (tmp_path / "analysis_plan.yaml").write_text("planned: [P1]\n")
    return RC.write_seal(tmp_path, root=tmp_path)


def test_seal_round_trips_and_is_plain_sha256sum_format(tmp_path):
    p = _seal_fixture(tmp_path)
    lines = p.read_text().strip().splitlines()
    assert len(lines) == 2
    for line in lines:
        h, _, rel = line.partition("  ")
        assert len(h) == 64 and all(c in "0123456789abcdef" for c in h)
    ok, problems = RC.verify_seal(tmp_path, root=tmp_path)
    assert ok, problems


def test_a_broken_seal_raises_and_is_never_a_warning(tmp_path):
    _seal_fixture(tmp_path)
    (tmp_path / "prediction.yaml").write_text("k_prime: 999.0\n")   # tampered
    ok, problems = RC.verify_seal(tmp_path, root=tmp_path)
    assert not ok and any("seal broken" in p for p in problems)
    with pytest.raises(RC.SealBroken, match="hard stop"):
        RC.verify_or_raise(tmp_path, root=tmp_path)


def test_both_documents_are_sealed_together(tmp_path):
    """Sealing the prediction alone leaves the choice of figure open."""
    (tmp_path / "prediction.yaml").write_text("k: 0\n")
    with pytest.raises(FileNotFoundError, match="analysis_plan.yaml"):
        RC.write_seal(tmp_path, root=tmp_path)


def test_an_empty_seal_is_refused(tmp_path):
    with pytest.raises((FileNotFoundError, ValueError)):
        RC.write_seal(tmp_path, docs=(), root=tmp_path)


def test_unsealed_is_permitted_but_reported_and_broken_never_is(tmp_path):
    """The asymmetry: 254 archived runs are unsealed and must stay re-runnable,
    but a broken seal is tampering."""
    notes = RC.verify_or_raise(tmp_path, require=False)
    assert any("unsealed" in n for n in notes)
    with pytest.raises(RC.SealBroken, match="Refusing to run"):
        RC.verify_or_raise(tmp_path, require=True)


def test_execute_calls_the_seal_check_from_inside_itself():
    """★ The lesson from `health.gate()`: a gate reachable only from a sibling tool
    is not enforcement. Assert the call site is inside `run.execute`."""
    src = (ROOT / "bdbot" / "run.py").read_text()
    body = src.split("def execute(")[1].split("\ndef ")[0]
    assert "verify_or_raise" in body, "the seal check is not inside execute()"
    assert "runcard" in body


def test_run_card_refuses_to_seal_without_a_prediction_or_a_plan():
    card = RC.RunCard(case="c", run_id="c__abc", question="Does X hold?",
                      cost_line="200 s vs 1000 s budget")
    b = card.blockers()
    assert any("no prediction" in x for x in b)
    assert any("no analysis plan" in x for x in b)
    assert any("cannot_decide" in x for x in b)
    assert "NOT SEALABLE" in card.render()


def test_a_prediction_without_a_role_is_a_blocker():
    card = RC.RunCard(case="c", run_id="r", question="Does X hold?",
                      cost_line="x", analysis=(("P1", "f", "d", "measurement"),),
                      cannot_decide=("nothing about hydrodynamics",),
                      predictions=(("K'", 0.0, 5.0, ""),))
    assert any("has no role" in x for x in card.blockers())


def test_a_complete_card_is_sealable_but_not_approved_by_code():
    """★ Updated 2026-09-02 for rule 10. This test previously asserted a card with
    NO parameter manifest was sealable -- i.e. it encoded the contract as it stood
    before "lay out every number, then ask" became a gate. A tightened gate
    breaking a test that asserted the looser rule is the gate working; the test is
    what was wrong, so it now supplies a manifest.
    """
    from bdbot import params as PRM
    man = PRM.Manifest("c", approved_by="a human, in conversation")
    man.add("geometry", "d", 5.0, "um", "sketch", 0)
    man.add("geometry", "L", 160.0, "um", "n_side * d", 0)
    man.add("geometry", "dim", 2, "1", "sketch: a 2D trap", 0)
    man.add("geometry", "N", 1000, "1", "choice: ensemble for statistics", 3)
    man.add("energy", "T", 300.0, "K", "sketch", 0)
    man.add("energy", "kT", 4.142e-21, "J", "k_B T at 300 K", 0)
    man.add("medium", "eta", 0.851e-3, "Pa*s", "Welty, water at 300 K", 0)
    man.add("interaction", "k_t", 10.0, "pN/um", "sketch, stated explicitly", 0)
    man.add("numerics", "dt", 8.02e-6, "s", "inverted from a 0.1 % EM bias", 0)
    man.add("numerics", "n_prod", 1_000_000, "1", "T_obs / dt", 0)
    man.add("numerics", "seed", 20260803, "1", "choice", 3)
    man.add("derived", "gamma", 4.0102e-8, "kg/s", "3 pi eta d", 0)
    man.add("derived", "D_t", 0.1033, "um**2/s", "kT/gamma", 0)
    man.add("derived", "tau_gov", 4.010e-3, "s",
            "tau_k = gamma/k governs, NOT tau_B = 242 s", 0)
    assert man.ready, man.blockers()

    card = RC.RunCard(
        case="c", run_id="r", question="Does X hold?",
        answering_quantity="f_c [Hz]", cost_line="200 s vs 1000 s",
        predictions=(("f_c", 39.69, 5.0, "implementation_check", "k/(2 pi gamma)"),),
        analysis=(("P1", "PSD with Lorentzian fit", "is f_c recovered",
                   "implementation_check"),),
        cannot_decide=("anything requiring hydrodynamic coupling",),
        manifest=man)
    assert card.blockers() == [], card.blockers()
    text = card.render()
    assert "SEALABLE" in text and "awaiting human approval" in text
    assert card.approved_by is None          # the CARD's approval is separate
    assert "PARAMETER MANIFEST" in text


# ══════════════════════════════════════════════════════════════════════════
# K — the analysis plan
# ══════════════════════════════════════════════════════════════════════════

VALID_PLAN = {
    "schema": AP.SCHEMA,
    "case": "synthetic",
    "answers_goal": "intake/synthetic/goal.yaml",
    "planned": [{
        "id": "P1",
        "figure": "PSD of x(t) with a Lorentzian fit, f_c marked",
        "estimator": "bdbot.microrheo.msd_to_gstar",
        "decides": "is f_c recovered from the trajectory",
        "role": "implementation_check",
        "answers": "f_c",
    }],
    "not_planned": [{"figure": "g(r)", "reason": "no pair interaction here"}],
    "post_hoc": [],
}


@pytest.fixture
def write_plan(tmp_path):
    def go(raw):
        if raw is None:
            return AP.load(tmp_path)
        (tmp_path / "analysis_plan.yaml").write_text(
            yaml.safe_dump(raw, allow_unicode=True))
        return AP.load(tmp_path)
    go.dir = tmp_path
    return go


def test_a_valid_plan_passes(write_plan):
    pl = write_plan(VALID_PLAN)
    assert pl.errors == [], [str(i) for i in pl.errors]
    assert pl.known


def test_an_unresolvable_estimator_is_caught_at_seal_time(write_plan):
    """★ The self-inflicted example: `docs/07`'s first draft named
    `simbot.analysis.trap.lockin_k_star`, which does not exist. The real one is
    `bdbot.lockin.k_star`."""
    raw = copy.deepcopy(VALID_PLAN)
    raw["planned"][0]["estimator"] = "simbot.analysis.trap.lockin_k_star"
    pl = write_plan(raw)
    assert pl.errors, "an unresolvable estimator must fail the plan"
    assert any("lockin_k_star" in i.msg or "AttributeError" in i.msg
               for i in pl.errors), [str(i) for i in pl.errors]

    raw["planned"][0]["estimator"] = "bdbot.lockin.k_star"
    assert write_plan(raw).errors == []


def test_resolve_estimator_messages_are_usable():
    assert AP.resolve_estimator("bdbot.microrheo.msd_to_gstar") is not None
    with pytest.raises(ValueError):
        AP.resolve_estimator("notdotted")
    with pytest.raises(ImportError):
        AP.resolve_estimator("bdbot.nosuchmodule.f")
    with pytest.raises(AttributeError):
        AP.resolve_estimator("bdbot.microrheo.nosuchfunction")


@pytest.mark.parametrize("field", ["id", "figure", "estimator", "decides", "role"])
def test_every_planned_item_field_is_required(write_plan, field):
    raw = copy.deepcopy(VALID_PLAN)
    raw["planned"][0][field] = ""
    pl = write_plan(raw)
    assert any(field in i.msg for i in pl.errors), [str(i) for i in pl.errors]


def test_an_invalid_role_is_caught(write_plan):
    raw = copy.deepcopy(VALID_PLAN)
    raw["planned"][0]["role"] = "probably_fine"
    assert write_plan(raw).errors


def test_an_omission_needs_a_reason(write_plan):
    """`FigureSet.skip()`'s discipline raised to the plan level: an absent figure
    with no reason is indistinguishable from forgetting."""
    raw = copy.deepcopy(VALID_PLAN)
    raw["not_planned"] = [{"figure": "g(r)"}]
    pl = write_plan(raw)
    assert any("reason" in i.msg for i in pl.errors)


def test_an_empty_plan_is_an_error(write_plan):
    raw = copy.deepcopy(VALID_PLAN)
    raw["planned"] = []
    assert write_plan(raw).errors


def test_the_plan_must_mention_the_goals_answering_quantity(write_plan, tmp_path):
    """★ The A -> K spine, enforced. A plan that never mentions the quantity the
    goal says answers the question cannot answer it."""
    goal_raw = {
        "schema": GO.SCHEMA, "question": "Does the PSD recover the stiffness?",
        "answering_quantity": {"symbol": "f_c", "unit": "Hz", "why": "w"},
        "decisive_precision": {"value": 5.0, "basis": "b",
                               "what_would_change_my_mind": "c"},
        "physics_that_matters": [{"symbol": "k_t", "why": "sets f_c"}],
        "analysis_implied": ["PSD with a Lorentzian fit"],
        "provenance": {"source": "s", "confirmed_by": None},
    }
    gdir = tmp_path / "g"
    gdir.mkdir()
    (gdir / "goal.yaml").write_text(yaml.safe_dump(goal_raw))
    goal = GO.load(gdir)
    assert goal.errors == []

    good = write_plan(VALID_PLAN)
    assert not AP.validate(good, goal=goal, check_imports=False) \
        or not [i for i in AP.validate(good, goal=goal, check_imports=False)
                if i.level == "error"]

    raw = copy.deepcopy(VALID_PLAN)
    raw["planned"][0]["answers"] = ""
    raw["planned"][0]["figure"] = "some unrelated histogram"
    bad = write_plan(raw)
    issues = AP.validate(bad, goal=goal, check_imports=False)
    assert any(i.level == "error" and "answering_quantity" in i.msg
               for i in issues), [str(i) for i in issues]


def test_post_hoc_may_not_be_promoted_into_planned(write_plan, tmp_path):
    """★ What makes the labelling load-bearing rather than cosmetic."""
    sealed_raw = copy.deepcopy(VALID_PLAN)
    sealed_raw["post_hoc"] = [{"id": "X1", "figure": "an interesting extra"}]
    sdir = tmp_path / "sealed"
    sdir.mkdir()
    (sdir / "analysis_plan.yaml").write_text(yaml.safe_dump(sealed_raw))
    sealed = AP.load(sdir)

    later = copy.deepcopy(VALID_PLAN)
    later["planned"].append({
        "id": "X1", "figure": "an interesting extra",
        "estimator": "bdbot.stats.block_sem", "decides": "d",
        "role": "measurement"})
    current = write_plan(later)

    viol = AP.post_hoc_is_clean(sealed, current)
    assert viol and any("moved from post_hoc" in v for v in viol), viol


def test_an_absent_plan_says_how_to_make_one(write_plan):
    pl = write_plan(None)
    assert not pl.exists and pl.errors
    assert "goal.yaml" in pl.errors[0].msg
    assert "ABSENT" in AP.render_check(pl)


def test_the_template_has_every_required_key(tmp_path):
    made, _ = AP.init_template(tmp_path, case="synthetic")
    assert made
    pl = AP.load(tmp_path)
    assert [i for i in pl.errors if i.msg == "required field missing"] == []
    assert pl.errors, "an untouched template must not pass as a real plan"
