"""L0 intake -- the `Observation` schema and its checks.

**The schema was derived from actual use, not guessed from a plan.** The basis is
the field-usage frequency across the five `intake/*/observation.yaml` files written
by hand from five real sketches:

    5/5 (required)  source_files . raw_transcription . system_guess . entities .
                stated_quantities · stated_goals · ambiguities · unread_regions ·
                missing_required
    2/5 (optional)  prerun_findings . references
    1/5 (optional)  hard_constraints . model_notes . contrast_with_*

The sub-keys were split the same way. The most important finding:

  * The `resolution` key is present in **all** of them -- 23/23 ambiguities,
    24/24 missing_required -- and its value is frequently `null`. That is the
    device that makes the reader say what they do not know.
    -> The check **enforces the key's presence** and allows a `null` value. An
    omission is rejected.

  * An `assumed_value` with no `confidence` is **an invented value with no tier.**
    That violates rule 3 in CLAUDE.md, so it is rejected.

pydantic was deliberately not used -- the value of this check is telling a person
exactly what is missing and why, and for a schema mixing optionality with
meaningful `null`s, a hand-written diagnosis is better.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

SCHEMA = "bdbot.observation/0.1"

# -- the schema (derived from actual usage frequency) --------------------------
REQUIRED_TOP = (
    "source_files", "raw_transcription", "system_guess", "entities",
    "stated_quantities", "stated_goals", "ambiguities", "unread_regions",
    "missing_required",
)
OPTIONAL_TOP = ("prerun_findings", "references", "hard_constraints", "model_notes")

# per list item: (required keys, optional keys). Required = present in 100% of the
# five files, and nothing else.
ITEM_KEYS = {
    "entities": (("kind",), ("note", "form", "expression", "shape", "center")),
    "stated_quantities": (("symbol", "value", "unit", "source"), ("confidence", "note")),
    "ambiguities": (("id", "issue", "resolution"),
                    ("impact", "note", "lean", "detail", "confirmed_by", "options",
                     "consequence", "evidence")),
    "missing_required": (("symbol", "resolution"),
                         ("assumed_value", "assumed_unit", "confidence", "what", "note",
                          "confirmed_by", "proposed", "followup", "lean", "kind")),
    "prerun_findings": (("id", "finding"), ("severity", "consequence", "verified", "options")),
    "references": (("kind",), ("authors", "doi", "title", "note", "resolved", "locator")),
}
# These keys must be **stated explicitly** even when empty -- they are the ones you
# must not omit.
MUST_BE_EXPLICIT = ("ambiguities", "unread_regions", "missing_required")

# The kinds of `missing_required`. * A distinction that became necessary after
# running the tool over the five real files -- the verdict was wrong: trap-2d-5um
# and soft-r3, both already completed end to end, came back BLOCKED. The cause was
# that things like `L` (box size) and `T_obs` (observation window) -- **simulation
# choices, not unknown physics** -- sat in the same list. The hand-written files had
# recorded that distinction only in prose, inside `note`.
#   physical  a property of the system. Without a human or the KB supplying it, L2
#             cannot be written -> **blocks**
#   choice    a simulation choice (box, observation window, sample count). Settled
#             in the numerics section -> does not block
MISSING_KINDS = ("physical", "choice")


@dataclass
class Issue:
    level: str          # error | warn | info
    where: str
    msg: str

    def __str__(self):
        mark = {"error": "✗", "warn": "⚠", "info": "ℹ"}[self.level]
        return f"  {mark} [{self.where}] {self.msg}"


@dataclass
class Observation:
    path: Path
    raw: dict
    issues: list = field(default_factory=list)

    # -- readiness ------------------------------------------------------------
    @property
    def open_ambiguities(self) -> list:
        return [a for a in (self.raw.get("ambiguities") or [])
                if isinstance(a, dict) and a.get("resolution") in (None, "")]

    @property
    def open_missing(self) -> list:
        """Unresolved gaps that block L2.

        An `assumed_value` counts as provisionally proceedable.
        `kind: choice` (box, observation window and other simulation choices) does
        not block -- see MISSING_KINDS above.
        The default is `physical` (conservative: an unstated kind blocks).
        """
        out = []
        for m in (self.raw.get("missing_required") or []):
            if not isinstance(m, dict):
                continue
            if m.get("kind", "physical") == "choice":
                continue
            if m.get("resolution") in (None, "") and m.get("assumed_value") in (None, ""):
                out.append(m)
        return out

    @property
    def open_choices(self) -> list:
        """Undecided simulation choices -- they do not block, but L3 must settle them."""
        return [m for m in (self.raw.get("missing_required") or [])
                if isinstance(m, dict) and m.get("kind") == "choice"
                and m.get("resolution") in (None, "") and m.get("assumed_value") in (None, "")]

    @property
    def assumed(self) -> list:
        """Values filled by assumption -- they must appear in the report with their
        tier (rule 3).
        """
        return [m for m in (self.raw.get("missing_required") or [])
                if isinstance(m, dict) and m.get("assumed_value") not in (None, "")]

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.level == "error"]

    @property
    def ready_for_system(self) -> bool:
        """Can `system.yaml` (L2) be written -- the schema must be intact and no
        unresolved gap may remain.
        """
        return not self.errors and not self.open_missing


def load(path) -> Observation:
    p = Path(path)
    if p.is_dir():
        p = p / "observation.yaml"
    if not p.exists():
        obs = Observation(p, {})
        obs.issues.append(Issue("error", str(p), "observation.yaml is missing. "
                                "Create a template with `bdbot intake init <folder>`."))
        return obs
    raw = yaml.safe_load(p.read_text()) or {}
    obs = Observation(p, raw)
    obs.issues = validate(obs)
    return obs


def validate(obs: Observation) -> list:
    """The schema check. A `null` value is allowed; **a missing key is rejected.**"""
    raw, root = obs.raw, obs.path.parent.parent.parent
    out: list[Issue] = []

    # 1. required top-level fields
    for k in REQUIRED_TOP:
        if k not in raw:
            hint = ("state it explicitly even if it is an empty list -- this is the "
                    "device that makes you say what you do not know"
                    if k in MUST_BE_EXPLICIT else "")
            out.append(Issue("error", k, f"required field missing. {hint}".strip()))

    # 2. was the transcription actually filled in
    tr = raw.get("raw_transcription")
    if isinstance(tr, str) and len(tr.strip()) < 20:
        out.append(Issue("error", "raw_transcription",
                         "the transcription is empty or too short. Rule 5: transcribe "
                         "before you interpret."))

    # 3. do the source files actually exist
    for s in (raw.get("source_files") or []):
        if not (obs.path.parent / Path(s).name).exists():
            out.append(Issue("warn", "source_files", f"file not found: {s}"))

    # 4. the keys of each list item
    for key, (req, opt) in ITEM_KEYS.items():
        items = raw.get(key)
        if items is None:
            continue
        if not isinstance(items, list):
            out.append(Issue("error", key, f"must be a list (got {type(items).__name__})"))
            continue
        known = set(req) | set(opt)
        seen_ids: set = set()
        for i, it in enumerate(items):
            at = f"{key}[{i}]"
            if not isinstance(it, dict):
                if key in ("unread_regions",):
                    continue
                out.append(Issue("error", at, "must be a mapping"))
                continue
            for r in req:
                if r not in it:
                    out.append(Issue("error", at, f"required key '{r}' missing "
                                                  f"(present in 100% of the five real files)"))
            for k2 in it:
                if k2 not in known:
                    out.append(Issue("info", at, f"key '{k2}' is not in the schema -- ignore if intended"))
            if "id" in it:
                if it["id"] in seen_ids:
                    out.append(Issue("error", at, f"duplicate id: {it['id']}"))
                seen_ids.add(it["id"])

    # 5. * the core defence -- reject an invented value with no confidence tier
    #    (CLAUDE.md rule 3)
    for i, m in enumerate(raw.get("missing_required") or []):
        if not isinstance(m, dict):
            continue
        at = f"missing_required[{i}]  {m.get('symbol', '?')}"
        has_val = m.get("assumed_value") not in (None, "")
        if has_val and m.get("confidence") is None:
            out.append(Issue("error", at,
                             "there is an assumed_value but no confidence (tier). "
                             "A value with no source does not go in (rule 3)."))
        if has_val and m.get("note") in (None, "") and m.get("resolution") in (None, ""):
            out.append(Issue("warn", at, "the assumed value's basis (note) is empty."))
        if m.get("kind") not in (None,) + MISSING_KINDS:
            out.append(Issue("error", at, f"kind must be one of {MISSING_KINDS}."))
        blocking = (m.get("kind", "physical") != "choice" and not has_val
                    and m.get("resolution") in (None, ""))
        if blocking and not (m.get("what") or m.get("note")):
            out.append(Issue("warn", at, "this item blocks L2 but what/note is empty -- "
                                         "there is no way to tell a person what is needed."))

    # 6. a value written on the sketch must carry a source (rule 3)
    for i, q in enumerate(raw.get("stated_quantities") or []):
        if not isinstance(q, dict):
            continue
        at = f"stated_quantities[{i}]  {q.get('symbol', '?')}"
        if not (q.get("source") or "").strip():
            out.append(Issue("error", at, "source is empty (rule 3: every number carries a source)."))
        if q.get("value") is None and q.get("unit") is not None:
            out.append(Issue("info", at, "value is null -- read as the sketch giving only a symbol."))

    # 7. be suspicious of a claim that there are no ambiguities at all
    if isinstance(raw.get("ambiguities"), list) and not raw["ambiguities"]:
        out.append(Issue("warn", "ambiguities",
                         "zero ambiguities. That is rare for a hand-drawn sketch -- look again."))
    return out


def render_check(obs: Observation) -> str:
    """The human-readable check report."""
    L: list[str] = []
    w = L.append
    w("=" * 78)
    w(f"intake check — {obs.path.parent.name}")
    w("=" * 78)
    if not obs.raw:
        w("\n".join(str(i) for i in obs.issues))
        return "\n".join(L)

    n_err = len(obs.errors)
    n_warn = len([i for i in obs.issues if i.level == "warn"])
    n_info = len([i for i in obs.issues if i.level == "info"])
    w(f"schema: {n_err} error(s) . {n_warn} warning(s) . {n_info} info")
    if obs.issues:
        w("")
        for i in obs.issues:
            w(str(i))

    w("")
    w("CONTENT SUMMARY")
    w(f"  transcription {len((obs.raw.get('raw_transcription') or '').splitlines())} line(s) . "
      f"entities {len(obs.raw.get('entities') or [])} . "
      f"stated quantities {len(obs.raw.get('stated_quantities') or [])} . "
      f"goals {len(obs.raw.get('stated_goals') or [])}")
    w(f"  ambiguities {len(obs.raw.get('ambiguities') or [])} "
      f"({len(obs.open_ambiguities)} unresolved) . "
      f"unread {len(obs.raw.get('unread_regions') or [])} . "
      f"gaps {len(obs.raw.get('missing_required') or [])} "
      f"({len(obs.open_missing)} unresolved, {len(obs.assumed)} filled by assumption)")

    if obs.assumed:
        w("")
        w("FILLED BY ASSUMPTION (tier needs checking -- rule 3)")
        for m in obs.assumed:
            unit = f" {m.get('assumed_unit')}" if m.get("assumed_unit") else ""
            w(f"  {m.get('symbol', '?'):<16} = {m.get('assumed_value')}{unit}"
              f"   [tier {m.get('confidence')}]  {(m.get('note') or '')[:44]}")

    if obs.open_missing:
        w("")
        w("* WHY THE PHYSICAL SYSTEM CANNOT BE SETTLED (stopped rather than inventing -- rule 3)")
        for m in obs.open_missing:
            desc = m.get("what") or m.get("note") or "(no description -- fill in `what`)"
            w(f"  {m.get('symbol', '?'):<16} {str(desc)[:56]}")
    if obs.open_choices:
        w("")
        w("UNDECIDED SIMULATION CHOICES (do not block -- settled at L3)")
        for m in obs.open_choices:
            desc = m.get("what") or m.get("proposed") or m.get("note") or ""
            w(f"  {m.get('symbol', '?'):<16} {str(desc)[:56]}")

    if obs.open_ambiguities:
        w("")
        w("UNRESOLVED AMBIGUITIES (for the first human check)")
        for a in obs.open_ambiguities:
            lean = f"  -> leaning: {a['lean']}" if a.get("lean") else ""
            w(f"  [{a.get('id', '?')}] {str(a.get('issue', ''))[:58]}{lean}")

    w("")
    w("=" * 78)
    if n_err:
        w(f"VERDICT: FAIL -- {n_err} schema error(s). Not advancing until they are fixed.")
    elif obs.open_missing:
        w(f"VERDICT: BLOCKED -- the schema is intact but {len(obs.open_missing)} gap(s) "
          f"are unresolved.")
        w("         L2 (system.yaml) cannot be written. A human must supply the value, "
          "or it must be found in the KB.")
    else:
        w("VERDICT: READY -- L2 (system.yaml) can be written.")
        if obs.open_ambiguities:
            w(f"         ({len(obs.open_ambiguities)} unresolved ambiguity/ies affect "
              f"how the result is interpreted)")
    w("=" * 78)
    return "\n".join(L)


TEMPLATE = """\
# L0 Observation -- what was read out of the sketch. **DRAFT: awaiting human confirmation**
# Protocol (skill `bd-intake`): transcribe first -> structure -> state ambiguities and
# unread regions -> anything absent stays null
#
# Check:  $PY -m bdbot.cli intake check intake/{case}
#   * ambiguities / unread_regions / missing_required must **keep their keys even when empty**
#   * using assumed_value makes confidence (tier) mandatory -- a value with no source does not go in
#   * leave resolution as null. That is the slot for saying explicitly that you do not know

source_files:
  - intake/{case}/sketch_01.jpeg

# -- 1. transcription (exactly as visible. No interpretation) -------------------
raw_transcription: |
  [Transcribe the characters you read in the image, exactly as written.
   Mark anything you could not read as [illegible] and add it to unread_regions.]

# -- 2. structuring ------------------------------------------------------------
system_guess: ""

entities: []          # - kind: particle / pair_interaction / external_potential / box ...
                      #   note: "..."

stated_quantities: []  # - symbol: k_t
                       #   value: 10
                       #   unit: pN/um
                       #   source: "right side of the sketch, written explicitly"
                       #   confidence: 1

stated_goals: []      # the measurement goals written on the sketch. If there are none, leave
                      # the list empty and record that in ambiguities -- an empty
                      # stated_goals is itself a blocker: without knowing what will be
                      # measured you cannot set T_obs, the sample interval, or a success
                      # criterion. Ask the user.

# -- 3. what is ambiguous (keep the key even when empty) -----------------------
ambiguities: []       # - id: A1
                      #   issue: "what is ambiguous"
                      #   impact: "how this choice changes the result -- write it as a number,
                      #            and split it per observable"
                      #   lean: "which way I lean (with the basis)"
                      #   resolution: null        # <- filled by a human
                      #
                      # Zero ambiguities is rare for a hand-drawn sketch, and the checker
                      # warns. Also: do not gain confidence from listing several bases in
                      # `lean` -- an exclusion argument can be strong while the remaining
                      # options stay undistinguished.

# -- 4. what could not be read -------------------------------------------------
unread_regions: []    # - "two characters at the bottom left, illegible"

# -- 5. values absent from the sketch that must be supplied --------------------
#       (do not invent -- record the source)
missing_required: []  # - symbol: eta
                      #   kind: physical          # physical blocks L2; choice does not
                      #   what: "solvent viscosity"
                      #   assumed_value: 0.851
                      #   assumed_unit: mPa*s
                      #   confidence: 1           # mandatory whenever assumed_value is present
                      #   note: "water@300K handbook. The medium is not stated in this sketch."
                      #   resolution: null
                      #
                      # tier: 0 given/handbook . 1 literature+verified, or a CONFIRMED
                      # convention . 2 literature unverified . 3 arbitrary assumption.
                      # Tier 1 by inheritance is the dangerous case -- T = 300 K is recorded
                      # as tier 1 across every case here and is actually a *choice* inherited
                      # from a sketch with no temperature, worth -4% to -14% on every
                      # timescale. Inheriting is fine; recording it as measured is not.
"""


def init_template(folder, case: str | None = None, force: bool = False) -> tuple[bool, str]:
    p = Path(folder)
    p.mkdir(parents=True, exist_ok=True)
    target = p / "observation.yaml"
    if target.exists() and not force:
        return False, f"already exists: {target}  (--force to overwrite)"
    target.write_text(TEMPLATE.format(case=case or p.name))
    return True, f"template created: {target}"


__all__ = ["SCHEMA", "Observation", "Issue", "load", "validate", "render_check",
           "init_template", "REQUIRED_TOP", "OPTIONAL_TOP", "ITEM_KEYS", "MUST_BE_EXPLICIT",
           "MISSING_KINDS"]
