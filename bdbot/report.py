"""The `DimensionlessReport` renderer -- the object of the third human check.

Both of the first two cases had grown the same skeleton:
  header -> reference scales + basis -> INPUT -> DERIVED -> SCALE LEDGER ->
  DIMENSIONLESS GROUPS → SEPARATION CHECKS → RUN PLAN → VERDICT

**What is common**: the frame, the ledger, the dimensionless groups, the checks,
the verdict (implemented below).
**What differs per case**: what goes into INPUT / DERIVED / RUN PLAN -> the case
passes those as line lists.

The width was unified here (one case used 84, the other 88). It is rendering
rather than a numerical result, so unifying it changes no result.
"""
from __future__ import annotations

from .checks import verdict as _verdict

W = 88            # total width
LEDGER_W = 38     # ledger entry name
GROUP_W = 44      # dimensionless group name
CHECK_W = 42      # check name


def render(*, title, ref, ledger, groups, checks,
           input_lines=(), derived_lines=(), run_plan_lines=()) -> tuple[str, str]:
    """Return the report text and the verdict.

    If the verdict is FAIL, the caller refuses to run.
    """
    out: list[str] = []
    w = out.append

    w("=" * W)
    w(title)
    w("=" * W)
    w(f"reference scales: length={ref['length'][0]}  energy={ref['energy'][0]}  time={ref['time'][0]}"
      + (f"   [strategy: {ref['strategy']}]" if ref.get("strategy") else ""))
    w(f"  basis: {ref.get('rationale', ledger.rationale)}")

    if input_lines:
        w("")
        w("INPUT (dimensional, SI)")
        for ln in input_lines:
            w(ln)

    if derived_lines:
        w("")
        w("DERIVED")
        for ln in derived_lines:
            w(ln)

    w("")
    w("SCALE LEDGER  (smallest first)")
    for cat_name, cat in ledger.categories():
        if not cat:
            continue
        w(f"  {cat_name}")
        for nm, v in ledger.sorted_items(cat):
            w(f"    {nm:<{LEDGER_W}} {v.to_compact():~.4gP}")

    w("")
    w("DIMENSIONLESS GROUPS  (each a ratio of two scales)")
    for nm, v in groups.items():
        w(f"  {nm:<{GROUP_W}} {v:.4g}")

    w("")
    w(f"{'SEPARATION CHECKS':<{CHECK_W + 14}}{'value':>10}{'limit':>10}{'margin':>9}")
    for c in checks:
        mark = "✓" if c.ok else ("✗" if c.hard else "⚠")
        w(f"  {mark} [{c.kind}] {c.name:<{CHECK_W}}{c.value:>10.3e}{c.limit:>10.0e}"
          f"{c.margin:>8.1f}×")
        if c.note:
            w(f"        {c.note}")

    if run_plan_lines:
        w("")
        w("RUN PLAN")
        for ln in run_plan_lines:
            w(ln)

    v, hard_fail, soft_fail, tight = _verdict(checks)
    w("")
    w(f"VERDICT: {v}")
    for c in hard_fail:
        w(f"  x HARD VIOLATION: {c.name.strip()} = {c.value:.3g} (limit {c.limit:g}) -- {c.note}")
    for c in soft_fail:
        w(f"  ! warning: {c.name.strip()} = {c.value:.3g} (limit {c.limit:g}) -- {c.note}")
    for c in tight:
        w(f"  ! thin margin: {c.name.strip()} ({c.margin:.1f}x)")
    w("=" * W)
    return "\n".join(out), v


def kv(key: str, value: str, tier=None, source: str = "", key_w: int = 6, val_w: int = 18) -> str:
    """One INPUT line. The tier and the source follow the value (rule 3)."""
    s = f"  {key:<{key_w}} = {value:<{val_w}}"
    if tier is not None:
        s += f" [tier {tier}]"
    if source:
        s += f" {source}"
    return s


__all__ = ["render", "kv", "W"]
