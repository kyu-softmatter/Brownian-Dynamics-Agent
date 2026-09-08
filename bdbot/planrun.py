"""Block K, second half — **executing** the sealed analysis plan. LLM 0 lines.

`analysisplan.py` says what will be measured. This module runs it. Together they
close the row the crosswalk left open: *"K analysis — exists, but plots hand-coded
per case."*

## What changes, and what deliberately does not

Before: each case script hand-wrote its own `make_plots()`. That is why the
crosswalk called K "exists but hand-coded" — the *figures* existed and the
*choice* of figure was reinvented per case, after the data was in hand.

After: the plan names an estimator and a figure; this module resolves the
estimator, calls it, and emits the figure with a mandatory caption. What stays
per case is the **data adapter** — how a run directory's arrays map onto an
estimator's arguments — because that genuinely differs per system and has not
appeared twice in the same shape.

So the promotion is narrow on purpose (CLAUDE.md: abstract only what appeared
twice). What is promoted is *plan -> figure*, not *run -> data*.

## Two rules carried over rather than reinvented

**An uncaptioned figure cannot be created.** `simbot/viz.py:84` rejects an empty
caption at `save()` time rather than auditing afterwards, *"because by then the
figure already exists."* `emit()` does the same: no caption, no file.

**A `post_hoc` figure is rendered as post-hoc.** Not hidden, not omitted --
labelled, in the caption itself, because that is the only difference between an
exploratory figure and a preregistered one.

## English labels, and why this is not a style preference

CLAUDE.md, measured: matplotlib's default `DejaVu Sans` has no Hangul so labels
render as `□`, and the fonts that do have Hangul are missing `−` (U+2212) and
`ŷ` (U+0177). `assert_ascii_labels()` therefore refuses a non-ASCII label instead
of producing a figure full of boxes.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

from . import analysisplan as _ap

SCHEMA = "bdbot.planrun/0.1"


class CaptionRequired(ValueError):
    """A figure was emitted without a caption. Never a warning."""


def assert_ascii_labels(*texts):
    """Refuse non-ASCII figure text. See the module docstring for the measurement."""
    for t in texts:
        s = str(t)
        bad = [c for c in s if ord(c) > 127]
        if bad:
            raise ValueError(
                f"non-ASCII characters {sorted(set(bad))!r} in figure text {s[:40]!r}. "
                f"matplotlib's DejaVu Sans renders Hangul as boxes, and the fonts "
                f"that have Hangul lack the minus sign (U+2212). Write axes, "
                f"legends and titles in English (CLAUDE.md).")
    return True


@dataclass
class Figure:
    """One emitted figure and its provenance."""
    id: str
    path: Path
    caption: str
    estimator: str
    role: str
    decides: str
    post_hoc: bool = False
    values: dict = field(default_factory=dict)

    def caption_line(self) -> str:
        head = "[POST HOC, not preregistered] " if self.post_hoc else ""
        return (f"{head}{self.caption}  "
                f"[{self.id} · role={self.role} · estimator={self.estimator}]")


@dataclass
class PlanRun:
    """The result of executing a plan: figures, skips and what could not run."""
    case: str
    figures: list = field(default_factory=list)
    skipped: list = field(default_factory=list)
    failed: list = field(default_factory=list)

    @property
    def ok(self) -> bool:
        return not self.failed

    def render(self) -> str:
        L = ["=" * 78, f"analysis plan executed — {self.case}", "=" * 78]
        L.append(f"\nEMITTED ({len(self.figures)})")
        for f in self.figures:
            L.append(f"  [{f.id}] {f.path.name}")
            L.append(f"        {f.caption_line()[:70]}")
            for k, v in list(f.values.items())[:4]:
                L.append(f"          {k} = {v}")
        if self.skipped:
            L.append(f"\nNOT PLANNED, with reasons ({len(self.skipped)})")
            for s in self.skipped:
                L.append(f"  - {s}")
        if self.failed:
            L.append(f"\n★ FAILED ({len(self.failed)}) -- a planned item that could "
                     f"not run is NOT a skip")
            for f in self.failed:
                L.append(f"  ✗ {f}")
        L.append("")
        L.append("=" * 78)
        L.append("VERDICT: " + ("OK" if self.ok else "INCOMPLETE -- see FAILED"))
        L.append("=" * 78)
        return "\n".join(L)


def emit(outdir, fig_id: str, caption: str, *, estimator: str, role: str,
         decides: str, values=None, post_hoc: bool = False, plot_fn=None,
         data=None) -> Figure:
    """Emit one figure. **Refuses without a caption.**

    `plot_fn(path, data, values)` does the drawing if given; with `plot_fn=None`
    only the record is produced, which is what lets a plan be executed and
    audited on a machine with no matplotlib.
    """
    if not str(caption).strip():
        raise CaptionRequired(
            f"figure {fig_id!r} has no caption. simbot/viz.py rejects this at "
            f"save() time rather than auditing afterwards, because by then the "
            f"figure already exists.")
    assert_ascii_labels(caption)
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    path = outdir / f"{fig_id}.png"
    f = Figure(id=str(fig_id), path=path, caption=str(caption),
               estimator=str(estimator), role=str(role), decides=str(decides),
               post_hoc=bool(post_hoc), values=dict(values or {}))
    if plot_fn is not None:
        plot_fn(path, data, f.values)
    return f


def execute(plan: _ap.Plan, outdir, adapter, *, plot_fn=None,
            require_clean: bool = True) -> PlanRun:
    """Run every `planned` item, then every `post_hoc` item.

    `adapter(item) -> (kwargs, caption_extra)` is the per-case data mapping: it
    turns a plan item into the keyword arguments its estimator needs. That stays
    per case deliberately -- see the module docstring.

    A planned item whose estimator raises lands in `failed`, **not** in `skipped`.
    The distinction matters: `not_planned` means "decided in advance not to", and
    `failed` means "said we would and did not". Collapsing them is how a missing
    result becomes invisible.
    """
    if require_clean and plan.errors:
        raise ValueError(
            f"refusing to execute a plan with {len(plan.errors)} error(s): "
            f"{plan.errors[0].msg[:80]}. Fix the plan before the run, not after.")

    run = PlanRun(case=str(plan.raw.get("case") or plan.path.parent.name))

    for it in plan.planned:
        if not isinstance(it, dict):
            continue
        fid = str(it.get("id", "?"))
        try:
            fn = _ap.resolve_estimator(it["estimator"])
            kwargs, extra = adapter(it)
            values = fn(**kwargs) if kwargs is not None else {}
            if not isinstance(values, dict):
                values = {"result": values}
            run.figures.append(emit(
                outdir, fid, str(it.get("figure", "")) + (extra or ""),
                estimator=str(it["estimator"]), role=str(it.get("role", "")),
                decides=str(it.get("decides", "")), values=values,
                post_hoc=False, plot_fn=plot_fn, data=values))
        except Exception as e:                       # noqa: BLE001
            run.failed.append(f"{fid}: {type(e).__name__}: {e}")

    for it in (plan.raw.get("not_planned") or []):
        if isinstance(it, dict):
            run.skipped.append(f"{it.get('figure', '?')}: {it.get('reason', '?')}")

    for it in plan.post_hoc:
        if not isinstance(it, dict):
            continue
        fid = str(it.get("id", "PH"))
        try:
            kwargs, extra = adapter(it)
            values = {}
            if it.get("estimator"):
                values = _ap.resolve_estimator(it["estimator"])(**kwargs) or {}
            run.figures.append(emit(
                outdir, fid, str(it.get("figure", "")) + (extra or ""),
                estimator=str(it.get("estimator", "(none)")),
                role=str(it.get("role", "measurement")),
                decides=str(it.get("decides", "exploratory")),
                values=values if isinstance(values, dict) else {"result": values},
                post_hoc=True, plot_fn=plot_fn, data=values))
        except Exception as e:                       # noqa: BLE001
            run.failed.append(f"{fid} (post_hoc): {type(e).__name__}: {e}")

    return run


def captions_md(run: PlanRun) -> str:
    """The caption block for `REPORT.md`. Post-hoc items are labelled inline."""
    L = [f"## Figures — {run.case}", ""]
    for f in run.figures:
        L.append(f"**{f.id}** `{f.path.name}`")
        L.append("")
        L.append(f.caption_line())
        L.append("")
    if run.skipped:
        L.append("### Not plotted, with reasons")
        L.append("")
        for s in run.skipped:
            L.append(f"- {s}")
        L.append("")
    if run.failed:
        L.append("### ★ Planned but not produced")
        L.append("")
        for f in run.failed:
            L.append(f"- {f}")
        L.append("")
    return "\n".join(L)


__all__ = ["SCHEMA", "CaptionRequired", "assert_ascii_labels", "Figure",
           "PlanRun", "emit", "execute", "captions_md"]
