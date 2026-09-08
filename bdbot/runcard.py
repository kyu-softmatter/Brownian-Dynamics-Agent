"""Block G — the run card, and the seal. LLM 0 lines.

The gate the pipeline did not have, and the cheapest place to close
[roadmap item 1](../docs/06-roadmap.md) — *"wire sealing into `bdbot`"*.

## What is sealed, and why here

At this point the goal is fixed (A), the physics and design are fixed (B–E), the
cost is known (F), and **nothing has executed**. So this is the last moment at
which a prediction can be written down without having seen a result, which is the
only moment at which writing one down means anything.

Two documents are sealed together: `prediction.yaml` and `analysis_plan.yaml`.
Sealing the prediction without the plan leaves the forking path open — you would
be committed to a number and free to choose which figure to show.

## Why the seal is verified by reading, not by importing

`SEALED.sha256` is plain `sha256sum` format. `.github/workflows/ci.yml:53` says
why: *"this check does not import a line of this repository's code -- which is the
point. A seal verified by the same codebase that wrote it is not verified."*
`simbot/io.py:341 write_seal` exists and is wired to `cli.py`, but `bdbot` cannot
import `simbot` (that direction is a cycle: `simbot` imports `bdbot` in four
modules). So this module re-implements *reading* the format — 20 lines of
`hashlib` — rather than reaching across the seam.

## Why the check goes INSIDE execute()

`health.gate()` is called from `tools/health.py:138` and `bdbot/run.py` contains
no reference to `gate` at all, so a run never gates itself. Its own docstring
records the consequence: *"an unwired checker cannot be wrong out loud."* A gate on
an optional path is not enforcement. `verify_or_raise` is therefore called from
`run.execute`, beside its existing `verify_hash()`.
"""
from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path

SCHEMA = "bdbot.runcard/0.1"

SEAL_NAME = "SEALED.sha256"

#: The documents sealed together at the run card. Order is the file order.
SEALED_DOCS = ("prediction.yaml", "analysis_plan.yaml")


class SealBroken(RuntimeError):
    """A sealed document changed after it was sealed. A hard stop, never a warning."""


def sha256_file(path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(65536), b""):
            h.update(chunk)
    return h.hexdigest()


def write_seal(rundir, docs=SEALED_DOCS, *, root=None) -> Path:
    """Write `SEALED.sha256` in plain `sha256sum` format.

    Paths are recorded relative to `root` (the repo root) so that
    `shasum -a 256 -c SEALED.sha256` verifies from there without this code.
    """
    rundir = Path(rundir)
    root = Path(root) if root else Path(__file__).resolve().parent.parent
    lines, missing = [], []
    for d in docs:
        p = rundir / d
        if not p.exists():
            missing.append(d)
            continue
        try:
            rel = p.resolve().relative_to(root)
        except ValueError:
            rel = p
        lines.append(f"{sha256_file(p)}  {rel}")
    if missing:
        raise FileNotFoundError(
            f"cannot seal {rundir.name}: missing {missing}. Both the prediction "
            f"and the analysis plan are sealed together -- sealing the prediction "
            f"alone leaves the choice of figure open (see the module docstring).")
    if not lines:
        raise ValueError("refusing to write an empty seal -- it would pass on nothing")
    out = rundir / SEAL_NAME
    out.write_text("\n".join(lines) + "\n")
    return out


def read_seal(rundir) -> list:
    """`[(want_hash, recorded_path)]` from `SEALED.sha256`. No repo code needed."""
    p = Path(rundir) / SEAL_NAME
    if not p.exists():
        return []
    out = []
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        want, _, rel = line.partition("  ")
        out.append((want.strip(), rel.strip()))
    return out


def verify_seal(rundir, *, root=None) -> tuple:
    """`(ok, problems)`. `ok` is True when a seal exists and every entry matches.

    A run with **no** seal is `(False, ["no seal"])` — absence is reported, not
    silently accepted, because "we never sealed it" and "the seal holds" are the
    two states this whole mechanism exists to distinguish.
    """
    rundir = Path(rundir)
    root = Path(root) if root else Path(__file__).resolve().parent.parent
    entries = read_seal(rundir)
    if not entries:
        return False, [f"{rundir.name}: no {SEAL_NAME} -- this run carries no "
                       f"sealed prediction"]
    problems, drifted = [], 0
    for want, rel in entries:
        target = root / rel
        if not target.exists():                 # renamed run dir: fall back, warn
            target = rundir / Path(rel).name
            drifted += 1
        if not target.exists():
            problems.append(f"sealed document missing: {rel}")
            continue
        got = sha256_file(target)
        if got != want:
            problems.append(f"seal broken for {rel}: recorded {want[:12]}…, "
                            f"computed {got[:12]}…")
    if drifted:
        problems.append(f"[warn] {drifted} recorded path(s) did not resolve from "
                        f"the repo root; verified next to the seal instead")
    hard = [p for p in problems if not p.startswith("[warn]")]
    return (not hard), problems


def verify_or_raise(rundir, *, root=None, require=True) -> list:
    """Raise `SealBroken` if a seal exists and does not hold.

    `require=False` permits an unsealed run (the 254 archived runs predate this
    module), but a **broken** seal always raises. That asymmetry is deliberate:
    unsealed is a known historical state, broken is tampering.
    """
    entries = read_seal(rundir)
    if not entries:
        if require:
            raise SealBroken(
                f"{Path(rundir).name}: no {SEAL_NAME}. Refusing to run without a "
                f"sealed prediction. Pass require_seal=False only for a "
                f"deliberately exploratory run, and record that choice.")
        return [f"[warn] {Path(rundir).name} is unsealed"]
    ok, problems = verify_seal(rundir, root=root)
    if not ok:
        raise SealBroken(
            f"{Path(rundir).name}: " + "; ".join(problems)
            + ". A broken seal is a hard stop: the comparison table is not built.")
    return problems


# -- the card itself ---------------------------------------------------------

@dataclass
class RunCard:
    """One screen a human reads before approving. Assembly only, no judgment."""
    case: str
    run_id: str
    question: str = ""
    answering_quantity: str = ""
    physical_lines: tuple = ()
    design_lines: tuple = ()
    dt_line: str = ""
    dt_bound_by: str = ""
    cost_line: str = ""
    predictions: tuple = ()        # (name, value, tol, role, basis)
    analysis: tuple = ()           # (id, figure, decides, role)
    cannot_decide: tuple = ()
    approved_by: str | None = None
    #: ★ rule 10. A `params.Manifest`. Without it the card cannot be sealed:
    #: "lay out all the numbers before you start the run, and ask."
    manifest: object = None

    def blockers(self) -> list:
        out = []
        if self.manifest is None:
            out.append("no parameter manifest -- rule 10 requires every number "
                       "(radius, box, kT, trap stiffness if applicable, dt, total "
                       "steps, ...) laid out before the run. See bdbot/params.py")
        else:
            out += [f"manifest: {b}" for b in self.manifest.blockers()]
        if not self.question.strip():
            out.append("no question -- block A is missing or empty")
        if not self.predictions:
            out.append("no prediction -- there is nothing to seal, so running now "
                       "makes the result un-preregisterable")
        if not self.analysis:
            out.append("no analysis plan -- the prediction would be sealed while "
                       "the choice of figure stays open")
        if not self.cost_line:
            out.append("no cost estimate (A3)")
        for p in self.predictions:
            if len(p) < 4 or not str(p[3]).strip():
                out.append(f"prediction {p[0] if p else '?'!r} has no role -- the "
                           f"role decides whether a mismatch is a bug or a result "
                           f"(rule 7')")
        if not self.cannot_decide:
            out.append("`cannot_decide` is empty -- state in advance which items "
                       "this design CANNOT settle, or the design-power check is "
                       "not being made")
        return out

    def render(self) -> str:
        W = 78
        L = ["=" * W, f"RUN CARD — {self.case}", f"run_id: {self.run_id}", "=" * W]
        L.append("\nTHE QUESTION")
        L.append(f"  {self.question or '(unset)'}")
        if self.answering_quantity:
            L.append(f"  answered by: {self.answering_quantity}")
        # ★ rule 10: the numbers come FIRST, before the design and the cost,
        #   because they are what the reader is being asked to approve.
        if self.manifest is not None:
            L.append("")
            L.append(self.manifest.render())
        if self.physical_lines:
            L.append("\nPHYSICAL SYSTEM (SI, with tiers)")
            L += [f"  {s}" for s in self.physical_lines]
        if self.design_lines:
            L.append("\nDESIGN — groups matched and abandoned")
            L += [f"  {s}" for s in self.design_lines]
        if self.dt_line:
            L.append("\nNUMERICS")
            L.append(f"  {self.dt_line}")
            if self.dt_bound_by:
                L.append(f"  bound by: {self.dt_bound_by}")
        if self.cost_line:
            L.append("\nCOST")
            L.append(f"  {self.cost_line}")
        L.append("\n★ PREDICTION (sealed on approval — this is the preregistration)")
        for p in self.predictions:
            name, val, tol, role = (list(p) + ["", "", "", ""])[:4]
            basis = p[4] if len(p) > 4 else ""
            L.append(f"  {str(name):<22} {str(val):<14} ±{str(tol):<8} [{role}]")
            if basis:
                L.append(f"      basis: {str(basis)[:60]}")
        L.append("\n★ ANALYSIS PLAN (sealed with it — fixed before the data exists)")
        for a in self.analysis:
            aid, fig, decides, role = (list(a) + ["", "", "", ""])[:4]
            L.append(f"  [{aid}] {str(fig)[:52]}")
            L.append(f"        decides {str(decides)[:48]}  [{role}]")
        L.append("\nWHAT THIS DESIGN CANNOT DECIDE")
        for c in self.cannot_decide or ("(unstated — this is a blocker)",):
            L.append(f"  - {c}")
        L.append("")
        L.append("=" * W)
        b = self.blockers()
        if b:
            L.append("VERDICT: NOT SEALABLE")
            for x in b:
                L.append(f"  ✗ {x}")
        elif self.approved_by:
            L.append(f"VERDICT: APPROVED by {self.approved_by} — sealed")
        else:
            L.append("VERDICT: SEALABLE — awaiting human approval.")
            L.append("         No code sets `approved_by` (CLAUDE.md, the "
                     "judgment split).")
        L.append("=" * W)
        return "\n".join(L)

    def write(self, rundir) -> Path:
        p = Path(rundir) / "RUNCARD.md"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text("```text\n" + self.render() + "\n```\n")
        return p


__all__ = ["SCHEMA", "SEAL_NAME", "SEALED_DOCS", "SealBroken", "sha256_file",
           "write_seal", "read_seal", "verify_seal", "verify_or_raise", "RunCard"]
