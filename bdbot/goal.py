"""Block A -- the goal. `intake/<case>/goal.yaml`, schema `bdbot.goal/0.1`.

**The first artifact of the pipeline, and the only one written before the sketch
is read.** It answers one question -- *what does the experimentalist want to know?*
-- and it is read again at the very end, by the analysis stage, to decide what to
plot. Nothing else in this repository spans the pipeline that way.

## Why it is a separate file from `observation.yaml`

They record different things and only one of them may be empty.

    stated_goals (observation.yaml)   what the SKETCH says. Rule 5: anything
                                      absent from the sketch stays null, so an
                                      empty list is a CORRECT transcription
    goal.yaml                         what the EXPERIMENTALIST wants. Never
                                      legitimately empty -- if it is unknown,
                                      ask; if it cannot be asked, stop

Measured: `trap-2d-5um`'s sketch carries no goal at all, and the goal it actually
ran for (fit the PSD) was settled in conversation. Under the first version of this
check (2026-09-02, blocking on `stated_goals` directly) that case was BLOCKED for
transcribing the sketch honestly -- punishing rule 5. This file is the fix: the
sketch may say nothing, and `goal.yaml` must then say something.

## What it must contain, and why each field is not optional

  question                      the thing to answer. Without it, `T_obs`, the
                                sample interval and the success criterion are all
                                unconstrained -- trap-2d-5um needed
                                `T_obs = 2000 tau_k` and `tau_k/10` sampling, and
                                neither follows from the sketch
  answering_quantity            which single measured number answers it. This is
                                what the analysis stage resolves against
                                `metrics.observable` names
  what_would_change_my_mind     the falsifiability gate S1 has always stated in
                                prose. A question that nothing could refute is not
                                a question
  physics_that_matters          the parameters the goal makes load-bearing. This
                                is the "extract physics" half -- the goal is what
                                tells you which numbers you actually care about
  analysis_implied              a HINT to the analysis stage, not the plan. The
                                plan is fixed later, once cost is known

`decisive_precision` may carry a null `value` -- often you do not know it until
the design exists -- but then `basis` has to say so. That field is deliberately
written here, before anything is known about what is affordable, so that it cannot
quietly become "whatever this design can reach".

`confirmed_by` is written by a human only. This module never sets it.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

SCHEMA = "bdbot.goal/0.1"

REQUIRED_TOP = ("question", "answering_quantity", "decisive_precision",
                "physics_that_matters", "analysis_implied", "provenance")
OPTIONAL_TOP = ("schema", "case", "asked_by", "out_of_scope", "notes")

ASKED_BY = ("experimentalist", "self", "paper")

# A question shorter than this is a label, not a question. 15 chars is not a
# principled number -- it is the shortest real one in the six non-empty
# `stated_goals` lists ("measure MSD, MSAD" is 17).
MIN_QUESTION_CHARS = 15


@dataclass
class Issue:
    level: str          # error | warn | info
    where: str
    msg: str

    def __str__(self):
        mark = {"error": "✗", "warn": "⚠", "info": "ℹ"}[self.level]
        return f"  {mark} [{self.where}] {self.msg}"


@dataclass
class Goal:
    path: Path
    raw: dict
    issues: list = field(default_factory=list)
    exists: bool = True

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.level == "error"]

    @property
    def question(self) -> str:
        return str(self.raw.get("question") or "").strip()

    @property
    def answering_symbol(self) -> str:
        aq = self.raw.get("answering_quantity") or {}
        return str(aq.get("symbol") or "").strip() if isinstance(aq, dict) else ""

    @property
    def confirmed(self) -> bool:
        """Has a human signed this off. `simbot.validate`'s convention: a judgment
        with `confirmed_by: null` does not count as confirmed, and no code fills it.
        """
        prov = self.raw.get("provenance") or {}
        return bool(isinstance(prov, dict) and prov.get("confirmed_by"))

    @property
    def known(self) -> bool:
        """Is the goal usable downstream -- the file exists and its schema holds.

        Deliberately does NOT require `confirmed`. An unconfirmed goal still tells
        the analysis stage what to aim at; it just is not signed off. Requiring
        confirmation here would block every case the moment it was drafted.
        """
        return self.exists and not self.errors


def load(path) -> Goal:
    p = Path(path)
    if p.is_dir():
        p = p / "goal.yaml"
    if not p.exists():
        return Goal(p, {}, issues=[Issue(
            "error", str(p.name),
            "goal.yaml is missing. Create it with `bdbot goal init <folder>`. "
            "The sketch may legitimately state no goal (rule 5) -- this file is "
            "where the answer to 'what do you want to know?' goes.")], exists=False)
    raw = yaml.safe_load(p.read_text()) or {}
    g = Goal(p, raw)
    g.issues = validate(g)
    return g


def validate(g: Goal) -> list:
    """Schema + content check. A `null` is allowed only where stated below."""
    raw = g.raw
    out: list[Issue] = []

    for k in REQUIRED_TOP:
        if k not in raw:
            out.append(Issue("error", k, "required field missing"))

    # 1. the question has to be a question
    q = g.question
    if not q:
        out.append(Issue("error", "question",
                         "empty. This is the field the whole pipeline exists to "
                         "answer; if it is unknown, ask rather than guess (rule 3)."))
    elif len(q) < MIN_QUESTION_CHARS:
        out.append(Issue("error", "question",
                         f"too short ({len(q)} chars) to be a question -- that is a "
                         f"label. Say what would be learned."))

    # 2. one named quantity answers it, and the analysis stage resolves that name
    aq = raw.get("answering_quantity")
    if aq is not None:
        if not isinstance(aq, dict):
            out.append(Issue("error", "answering_quantity", "must be a mapping"))
        else:
            if not g.answering_symbol:
                out.append(Issue("error", "answering_quantity.symbol",
                                 "empty. The analysis stage matches this against a "
                                 "measured observable name; without it the goal "
                                 "cannot be checked against any result."))
            if not str(aq.get("unit") or "").strip():
                out.append(Issue("warn", "answering_quantity.unit",
                                 "empty. Use '1' for a dimensionless quantity, so "
                                 "that blank means forgotten rather than dimensionless."))
            if not str(aq.get("why") or "").strip():
                out.append(Issue("warn", "answering_quantity.why",
                                 "empty -- why does this quantity answer the question?"))

    # 3. ★ the falsifiability gate, in code for the first time
    dp = raw.get("decisive_precision")
    if dp is not None:
        if not isinstance(dp, dict):
            out.append(Issue("error", "decisive_precision", "must be a mapping"))
        else:
            change = str(dp.get("what_would_change_my_mind") or "").strip()
            if not change:
                out.append(Issue("error", "decisive_precision.what_would_change_my_mind",
                                 "empty. A question nothing could refute is not a "
                                 "question -- this is S1's falsifiability gate."))
            basis = str(dp.get("basis") or "").strip()
            if dp.get("value") is None and not basis:
                out.append(Issue("error", "decisive_precision",
                                 "value is null and basis is empty. A null value is "
                                 "allowed -- often it is unknown until the design "
                                 "exists -- but then basis must say that."))

    # 4. the goal names the parameters it makes load-bearing
    pm = raw.get("physics_that_matters")
    if pm is not None:
        if not isinstance(pm, list):
            out.append(Issue("error", "physics_that_matters", "must be a list"))
        elif not pm:
            out.append(Issue("error", "physics_that_matters",
                             "empty. The point of stating a goal first is that it "
                             "tells you which parameters matter."))
        else:
            for i, it in enumerate(pm):
                at = f"physics_that_matters[{i}]"
                if not isinstance(it, dict):
                    out.append(Issue("error", at, "must be a mapping {symbol, why}"))
                    continue
                for need in ("symbol", "why"):
                    if not str(it.get(need) or "").strip():
                        out.append(Issue("error", f"{at}  {it.get('symbol', '?')}",
                                         f"'{need}' is empty"))

    # 5. a hint for the analysis stage -- not the plan, but not nothing either
    ai = raw.get("analysis_implied")
    if ai is not None:
        if not isinstance(ai, list):
            out.append(Issue("error", "analysis_implied", "must be a list"))
        elif not ai:
            out.append(Issue("error", "analysis_implied",
                             "empty. This is the field the analysis stage reads to "
                             "decide what to plot; with nothing here the goal cannot "
                             "reach the end of the pipeline."))

    # 6. provenance, and the field no code may write
    prov = raw.get("provenance")
    if prov is not None:
        if not isinstance(prov, dict):
            out.append(Issue("error", "provenance", "must be a mapping"))
        else:
            if not str(prov.get("source") or "").strip():
                out.append(Issue("error", "provenance.source",
                                 "empty. Where did this goal come from -- a "
                                 "conversation, a paper, the sketch? (rule 3)"))
            if "confirmed_by" not in prov:
                out.append(Issue("warn", "provenance.confirmed_by",
                                 "key absent. State it as null explicitly -- that is "
                                 "the slot recording that nobody has signed off yet."))

    ab = raw.get("asked_by")
    if ab is not None and ab not in ASKED_BY:
        out.append(Issue("warn", "asked_by", f"expected one of {ASKED_BY}, got {ab!r}"))

    sc = raw.get("schema")
    if sc is not None and sc != SCHEMA:
        out.append(Issue("warn", "schema", f"expected {SCHEMA}, got {sc!r}"))

    return out


def status(case_dir, obs=None) -> tuple[str, str]:
    """Is the goal known for this case -- `(verdict, reason)`.

    `verdict` is one of `CONFIRMED` · `DRAFT` · `FAIL` · `ABSENT`.
    Only `ABSENT` and `FAIL` block. This is the single place that answers the
    question, so a caller cannot re-derive it and drift (which is exactly what
    happened to `intake.ready_for_system`).

    `obs` is accepted but not required: a goal written on the sketch is evidence
    that the goal is known, and is reported, but it does not substitute for this
    file -- a transcription is not a decision.
    """
    g = load(case_dir)
    if not g.exists:
        sketch = bool(obs is not None and (obs.raw.get("stated_goals") or []))
        hint = (" (the sketch does state a goal -- transcribe it into goal.yaml "
                "and say what would refute it)" if sketch else "")
        return "ABSENT", f"goal.yaml does not exist{hint}"
    if g.errors:
        return "FAIL", f"{len(g.errors)} schema error(s): {g.errors[0].msg[:60]}"
    if not g.confirmed:
        return "DRAFT", "schema holds; provenance.confirmed_by is null (awaiting a human)"
    return "CONFIRMED", f"confirmed by {g.raw['provenance']['confirmed_by']}"


def blocks(case_dir, obs=None) -> bool:
    """Does the goal block L2. `ABSENT` and `FAIL` block; `DRAFT` does not."""
    return status(case_dir, obs)[0] in ("ABSENT", "FAIL")


def render_check(g: Goal) -> str:
    L: list[str] = []
    w = L.append
    w("=" * 78)
    w(f"goal check — {g.path.parent.name}")
    w("=" * 78)
    if not g.exists:
        w("")
        for i in g.issues:
            w(str(i))
        w("")
        w("=" * 78)
        w("VERDICT: ABSENT -- block A has not been written for this case.")
        w("=" * 78)
        return "\n".join(L)

    n_err = len(g.errors)
    n_warn = len([i for i in g.issues if i.level == "warn"])
    w(f"schema: {n_err} error(s) . {n_warn} warning(s)")
    if g.issues:
        w("")
        for i in g.issues:
            w(str(i))

    w("")
    w("THE QUESTION")
    for line in _wrap(g.question or "(empty)", 74):
        w(f"  {line}")

    aq = g.raw.get("answering_quantity") or {}
    if isinstance(aq, dict) and aq:
        w("")
        w("ANSWERED BY")
        unit = f" [{aq.get('unit')}]" if aq.get("unit") else ""
        w(f"  {aq.get('symbol', '?')}{unit}")
        if aq.get("why"):
            for line in _wrap(str(aq["why"]), 70):
                w(f"      {line}")

    dp = g.raw.get("decisive_precision") or {}
    if isinstance(dp, dict) and dp:
        w("")
        w("DECISIVE PRECISION")
        w(f"  value  {dp.get('value') if dp.get('value') is not None else 'not yet known'}")
        if dp.get("basis"):
            for line in _wrap(f"basis  {dp['basis']}", 70):
                w(f"  {line}")
        if dp.get("what_would_change_my_mind"):
            w("  * WOULD CHANGE MY MIND")
            for line in _wrap(str(dp["what_would_change_my_mind"]), 70):
                w(f"      {line}")

    pm = g.raw.get("physics_that_matters") or []
    if isinstance(pm, list) and pm:
        w("")
        w("PHYSICS THE GOAL MAKES LOAD-BEARING")
        for it in pm:
            if isinstance(it, dict):
                w(f"  {str(it.get('symbol', '?')):<14} {str(it.get('why', ''))[:56]}")

    ai = g.raw.get("analysis_implied") or []
    if isinstance(ai, list) and ai:
        w("")
        w("IMPLIED ANALYSIS (a hint for the analysis stage -- not the plan)")
        for a in ai:
            w(f"  - {str(a)[:70]}")

    oos = g.raw.get("out_of_scope") or []
    if isinstance(oos, list) and oos:
        w("")
        w("OUT OF SCOPE")
        for a in oos:
            w(f"  - {str(a)[:70]}")

    verdict, reason = status(g.path.parent)
    w("")
    w("=" * 78)
    if verdict == "FAIL":
        w(f"VERDICT: FAIL -- {n_err} schema error(s). Not advancing until fixed.")
    elif verdict == "DRAFT":
        w("VERDICT: DRAFT -- the schema holds and downstream stages may use it,")
        w("         but `provenance.confirmed_by` is null. A human has not signed")
        w("         this off, and no code will do it (CLAUDE.md, the judgment split).")
    else:
        w(f"VERDICT: CONFIRMED -- {reason}")
    w("=" * 78)
    return "\n".join(L)


def _wrap(text: str, width: int) -> list:
    out, line = [], ""
    for word in str(text).split():
        if line and len(line) + 1 + len(word) > width:
            out.append(line)
            line = word
        else:
            line = f"{line} {word}".strip()
    if line:
        out.append(line)
    return out


TEMPLATE = """\
# Block A -- the goal. **DRAFT: awaiting human confirmation**
#
# Written BEFORE the sketch is read, and read again at the very end by the
# analysis stage. The sketch may legitimately state no goal (rule 5); this file
# is where the answer to "what do you want to know?" goes instead.
#
# Check:  $PY -m bdbot.cli goal check intake/{case}
#
# `confirmed_by` is filled by a human only. No code in this repository writes it.

schema: {schema}
case: {case}
asked_by: experimentalist        # experimentalist | self | paper

# -- what do you want to know? ------------------------------------------------
# One sentence. If you cannot write it, ask rather than guess (rule 3).
question: ""

# -- which single measured number answers it? ---------------------------------
# `symbol` is matched against a measured observable name at the end of the
# pipeline, so it has to be the name the analysis will actually produce.
answering_quantity:
  symbol: ""
  unit: ""                       # use "1" for dimensionless -- blank means forgotten
  why: ""

# -- how well would you have to measure it to decide? -------------------------
# `value` may be null (often unknown until the design exists) but then `basis`
# must say so. Written here, before cost is known, so it cannot quietly become
# "whatever this design can reach".
decisive_precision:
  value: null
  basis: ""
  what_would_change_my_mind: ""  # * the falsifiability gate. Required

# -- which parameters does this goal make load-bearing? -----------------------
# The "extract physics" half: the goal is what tells you which numbers matter.
physics_that_matters: []         # - symbol: k_t
                                 #   why: "sets the corner frequency the PSD fit needs"

# -- what analysis does this imply? -------------------------------------------
# A HINT for the analysis stage, not the plan. The plan is fixed later, once the
# cost of the run is known.
analysis_implied: []             # - "PSD with a Lorentzian fit, f_c marked"

out_of_scope: []                 # - "anything requiring hydrodynamic coupling"

provenance:
  source: ""                     # a conversation (with a date), a paper, the sketch
  confirmed_by: null             # <- a human only
"""


def init_template(folder, case: str | None = None, force: bool = False) -> tuple[bool, str]:
    p = Path(folder)
    p.mkdir(parents=True, exist_ok=True)
    target = p / "goal.yaml"
    if target.exists() and not force:
        return False, f"already exists: {target}  (--force to overwrite)"
    target.write_text(TEMPLATE.format(case=case or p.name, schema=SCHEMA))
    return True, f"template created: {target}"


__all__ = ["SCHEMA", "Goal", "Issue", "load", "validate", "render_check", "status",
           "blocks", "init_template", "REQUIRED_TOP", "OPTIONAL_TOP", "ASKED_BY",
           "MIN_QUESTION_CHARS"]
