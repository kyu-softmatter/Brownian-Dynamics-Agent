"""Block K — the analysis plan, fixed before the data exists. LLM 0 lines.

Derived from `goal.yaml` (block A), sealed at the run card (block G), executed
after the run. Not the same object as the goal's `analysis_implied`, which is a
hint written before anything is known about cost.

## Why the plan is sealed rather than chosen afterwards

*"Decide the best analysis for this system from the goal"* is principled.
*"Decide it after seeing the trajectory"* is a garden of forking paths, and this
repository has the branches for it: 145 runs in `chain-bend-2d-dlvo`, **3,856** in
the `soft-r3` campaign. With that many, something always looks significant.

`post_hoc` is not forbidden — it is **labelled**, which is the whole difference.
An item may move from `planned` to nowhere, and from nowhere to `post_hoc`. It may
never move from `post_hoc` to `planned`; `validate` checks that against the seal.

## Why every estimator must be importable

The first draft of `docs/07` named `simbot.analysis.trap.lockin_k_star` as an
estimator. No such function exists — the real one is `bdbot/lockin.py:52
k_star()`. A plan naming an unresolvable estimator is a plan nobody can execute,
and it looks exactly like one that can. So `resolve_estimator` imports it at
**seal** time, not at plot time.

## Roles, and why `not_planned` needs a reason

Each item carries a `role` from `metrics.ROLES`, because the role decides what a
mismatch *means* (rule 7'). And an omission needs a reason for the same purpose
`simbot/viz.py:99 FigureSet.skip()` serves: *"if not drawing the RDF was 'no pair
interaction' rather than 'forgot', an absent figure is a silent omission."*
"""
from __future__ import annotations

import importlib
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import metrics as _metrics

SCHEMA = "bdbot.analysisplan/0.1"

REQUIRED_TOP = ("answers_goal", "planned", "not_planned")
OPTIONAL_TOP = ("schema", "case", "post_hoc", "notes")

ITEM_REQUIRED = ("id", "figure", "estimator", "decides", "role")
ITEM_OPTIONAL = ("overlays", "note", "unit", "answers")


@dataclass
class Issue:
    level: str
    where: str
    msg: str

    def __str__(self):
        mark = {"error": "✗", "warn": "⚠", "info": "ℹ"}[self.level]
        return f"  {mark} [{self.where}] {self.msg}"


@dataclass
class Plan:
    path: Path
    raw: dict
    issues: list = field(default_factory=list)
    exists: bool = True

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.level == "error"]

    @property
    def planned(self) -> list:
        return list(self.raw.get("planned") or [])

    @property
    def post_hoc(self) -> list:
        return list(self.raw.get("post_hoc") or [])

    @property
    def known(self) -> bool:
        return self.exists and not self.errors

    def ids(self) -> set:
        return {str(i.get("id")) for i in self.planned if isinstance(i, dict)}


def resolve_estimator(dotted: str):
    """Import `pkg.mod.func` and return it. Raises with a usable message.

    Called at seal time so an unresolvable name is caught before the run, not
    after — see the module docstring.
    """
    s = str(dotted).strip()
    if "." not in s:
        raise ValueError(f"estimator must be a dotted path, got {s!r}")
    mod_name, _, attr = s.rpartition(".")
    try:
        mod = importlib.import_module(mod_name)
    except ImportError as e:
        raise ImportError(
            f"estimator {s!r}: cannot import module {mod_name!r} ({e}). "
            f"A plan naming an unresolvable estimator looks exactly like one "
            f"that can be executed.") from e
    if not hasattr(mod, attr):
        raise AttributeError(
            f"estimator {s!r}: {mod_name!r} has no attribute {attr!r}. "
            f"Available: {', '.join(sorted(n for n in dir(mod) if not n.startswith('_'))[:12])}…")
    return getattr(mod, attr)


def load(path) -> Plan:
    p = Path(path)
    if p.is_dir():
        p = p / "analysis_plan.yaml"
    if not p.exists():
        return Plan(p, {}, issues=[Issue(
            "error", p.name,
            "analysis_plan.yaml is missing. It is derived from goal.yaml and "
            "sealed with the prediction at the run card (block G).")], exists=False)
    raw = yaml.safe_load(p.read_text()) or {}
    pl = Plan(p, raw)
    pl.issues = validate(pl)
    return pl


def validate(pl: Plan, *, goal=None, check_imports=True) -> list:
    """Schema + content. Pass `goal` (a `bdbot.goal.Goal`) to also check coverage."""
    raw = pl.raw
    out: list[Issue] = []

    for k in REQUIRED_TOP:
        if k not in raw:
            out.append(Issue("error", k, "required field missing"))

    sc = raw.get("schema")
    if sc is not None and sc != SCHEMA:
        out.append(Issue("warn", "schema", f"expected {SCHEMA}, got {sc!r}"))

    planned = raw.get("planned")
    if planned is not None:
        if not isinstance(planned, list):
            out.append(Issue("error", "planned", "must be a list"))
        elif not planned:
            out.append(Issue("error", "planned",
                             "empty. A run with no planned analysis cannot answer "
                             "the goal's question."))
        else:
            seen = set()
            for i, it in enumerate(planned):
                at = f"planned[{i}]"
                if not isinstance(it, dict):
                    out.append(Issue("error", at, "must be a mapping"))
                    continue
                at = f"planned[{i}] {it.get('id', '?')}"
                for need in ITEM_REQUIRED:
                    if not str(it.get(need) or "").strip():
                        out.append(Issue("error", at, f"'{need}' is empty"))
                r = it.get("role")
                if r and r not in _metrics.ROLES:
                    out.append(Issue("error", at,
                                     f"role must be one of {_metrics.ROLES}, got {r!r}"))
                iid = str(it.get("id") or "")
                if iid and iid in seen:
                    out.append(Issue("error", at, f"duplicate id {iid!r}"))
                seen.add(iid)
                est = it.get("estimator")
                if est and check_imports:
                    try:
                        resolve_estimator(est)
                    except Exception as e:                    # noqa: BLE001
                        out.append(Issue("error", at, f"{type(e).__name__}: {e}"))

    npl = raw.get("not_planned")
    if npl is not None:
        if not isinstance(npl, list):
            out.append(Issue("error", "not_planned", "must be a list"))
        else:
            for i, it in enumerate(npl):
                at = f"not_planned[{i}]"
                if not isinstance(it, dict):
                    out.append(Issue("error", at, "must be a mapping {figure, reason}"))
                    continue
                if not str(it.get("reason") or "").strip():
                    out.append(Issue("error", f"{at} {it.get('figure', '?')}",
                                     "'reason' is empty. An omission with no reason "
                                     "is indistinguishable from forgetting "
                                     "(simbot/viz.py FigureSet.skip)."))

    ph = raw.get("post_hoc")
    if ph is not None and not isinstance(ph, list):
        out.append(Issue("error", "post_hoc", "must be a list"))

    # coverage against the goal
    if goal is not None:
        implied = [str(x) for x in (goal.raw.get("analysis_implied") or [])]
        covered = " \n".join(
            [str(i.get("figure", "")) for i in pl.planned if isinstance(i, dict)]
            + [str(i.get("figure", "")) for i in (raw.get("not_planned") or [])
               if isinstance(i, dict)])
        for h in implied:
            head = h.split(",")[0].split(" with ")[0].strip().lower()
            if head and head[:18] not in covered.lower():
                out.append(Issue("warn", "answers_goal",
                                 f"goal.analysis_implied item not obviously covered "
                                 f"in planned or not_planned: {h[:56]!r}"))
        sym = goal.answering_symbol
        if sym:
            answers = " ".join(str(i.get("answers", "")) + " " + str(i.get("figure", ""))
                               for i in pl.planned if isinstance(i, dict))
            if sym.lower() not in answers.lower():
                out.append(Issue("error", "answers_goal",
                                 f"no planned item mentions the goal's "
                                 f"answering_quantity {sym!r}. The plan cannot "
                                 f"answer the question it is derived from."))
    return out


def post_hoc_is_clean(sealed_plan: Plan, current: Plan) -> list:
    """An item may not migrate from `post_hoc` into `planned` after sealing.

    Returns violation strings; empty means clean. This is the check that makes
    the labelling load-bearing rather than cosmetic.
    """
    sealed_ids = sealed_plan.ids()
    sealed_ph = {str(i.get("id")) for i in sealed_plan.post_hoc
                 if isinstance(i, dict)}
    out = []
    for iid in current.ids():
        if iid not in sealed_ids:
            out.append(f"planned item {iid!r} was not in the sealed plan "
                       f"(new planned items after sealing are post_hoc)")
        if iid in sealed_ph:
            out.append(f"item {iid!r} moved from post_hoc into planned after "
                       f"sealing -- forbidden")
    return out


def render_check(pl: Plan) -> str:
    L = ["=" * 78, f"analysis plan — {pl.path.parent.name}", "=" * 78]
    if not pl.exists:
        L += [""] + [str(i) for i in pl.issues]
        L += ["", "=" * 78, "VERDICT: ABSENT -- block K has no plan for this run.", "=" * 78]
        return "\n".join(L)
    n_err = len(pl.errors)
    n_warn = len([i for i in pl.issues if i.level == "warn"])
    L.append(f"schema: {n_err} error(s) . {n_warn} warning(s)")
    L.append(f"answers goal: {pl.raw.get('answers_goal', '(unset)')}")
    if pl.issues:
        L.append("")
        L += [str(i) for i in pl.issues]
    L.append("")
    L.append("PLANNED (fixed before the data exists)")
    for it in pl.planned:
        if isinstance(it, dict):
            L.append(f"  [{it.get('id','?'):<4}] {str(it.get('figure',''))[:52]}")
            L.append(f"         estimator {it.get('estimator','?')}")
            L.append(f"         decides   {str(it.get('decides',''))[:52]}"
                     f"   role={it.get('role','?')}")
    npl = pl.raw.get("not_planned") or []
    if npl:
        L.append("")
        L.append("NOT PLANNED (with a reason -- an omission needs one)")
        for it in npl:
            if isinstance(it, dict):
                L.append(f"  - {str(it.get('figure',''))[:40]}: "
                         f"{str(it.get('reason',''))[:44]}")
    if pl.post_hoc:
        L.append("")
        L.append("POST HOC (added after the run -- rendered as post-hoc)")
        for it in pl.post_hoc:
            L.append(f"  - {str(it.get('figure', it))[:66]}")
    L.append("")
    L.append("=" * 78)
    L.append("VERDICT: FAIL -- fix the errors" if n_err else "VERDICT: OK -- sealable")
    L.append("=" * 78)
    return "\n".join(L)


TEMPLATE = """\
# Block K -- the analysis plan. Derived from goal.yaml, SEALED at the run card.
#
# `planned` is fixed before the data exists. Adding to it after the run is
# forbidden; add to `post_hoc` instead, which is rendered as post-hoc.
# Every `estimator` is imported when the plan is sealed, so it must resolve.

schema: {schema}
case: {case}
answers_goal: intake/{case}/goal.yaml

planned:
  - id: P1
    figure: ""                 # axes, scale, overlays -- what a reader will see
    estimator: ""              # dotted path, e.g. bdbot.microrheo.msd_to_gstar
    decides: ""                # what question this figure settles
    role: measurement          # implementation_check | hypothesis | measurement
    answers: ""                # the goal's answering_quantity, if this is it

not_planned: []                # - figure: "g(r)"
                               #   reason: "no pair interaction -- uninformative here"

post_hoc: []                   # appended AFTER the run. Never promoted to planned.
"""


def init_template(folder, case: str | None = None, force: bool = False) -> tuple:
    p = Path(folder)
    p.mkdir(parents=True, exist_ok=True)
    target = p / "analysis_plan.yaml"
    if target.exists() and not force:
        return False, f"already exists: {target}  (--force to overwrite)"
    target.write_text(TEMPLATE.format(case=case or p.name, schema=SCHEMA))
    return True, f"template created: {target}"


__all__ = ["SCHEMA", "REQUIRED_TOP", "OPTIONAL_TOP", "ITEM_REQUIRED", "Plan",
           "Issue", "load", "validate", "resolve_estimator", "post_hoc_is_clean",
           "render_check", "init_template"]
