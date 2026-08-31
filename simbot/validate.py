"""S7 — the sealed prediction vs the measurement. **It proposes a verdict; it never
confirms one.**

## What this module never does

It does not fill in `confirmed_by`. Only a human writes that field
(CLAUDE.md §verdicts). A verdict with `confirmed_by: null` does not enter the
benchmark-ledger tally -- count it and a pass stamp gets applied that no human ever
looked at.

## `INCONCLUSIVE` is not a failure

On the first completed run, 2 of 9 predictions came out `INCONCLUSIVE`, and **for
both the lack of discriminating power was already foreseen when the prediction
document was written.** If the statistical error exceeds the tolerance then that
measurement cannot decide that tolerance -- calling it `PASS` is luck, not
verification.

Two things are computed separately:

| | question | on failure |
|---|---|---|
| **decidability** | is `SE` smaller than the tolerance band | `INCONCLUSIVE` |
| **design power** | the gap to the competing hypothesis, in `SE` | `INCONCLUSIVE` + sample multiple |

★ Do not demand a `3σ` rejection where the power cannot produce `3σ` -- that is an
unachievable assert. In that regime `INCONCLUSIVE` is **pinned as a fact**
(CLAUDE.md, the 4 statistics rules).
"""
from __future__ import annotations

import math
import re
from dataclasses import dataclass, field

from .io import RunDir, SealVerdict, verify_seal
from .spec import Prediction, PredictionItem

PASS, FAIL, INCONCLUSIVE = "PASS", "FAIL", "INCONCLUSIVE"

# The 4 FAIL cause classes (master_plan §S7-6). A FAIL must always carry one.
CAUSE_CLASSES: frozenset[str] = frozenset({
    "numerical", "modeling", "interpretation", "analysis", "environment"})


# =============================================================================
# parsing a tolerance
# =============================================================================
@dataclass(frozen=True)
class Tolerance:
    """The machine form of a string like `"±1.5%"`.

    kind:
      relative     — ±X % of the predicted value
      absolute     — ±X (in the measurement's unit)
      lower_bound  — the measurement must exceed this (`R² > 0.99`, `p > 0.05`)
      upper_bound  — the measurement must be below this
    """

    kind: str
    magnitude: float
    text: str

    def half_width(self, predicted: float | None) -> float | None:
        """Half-width of the decision band. `None` for a one-sided bound, which
        has no band."""
        if self.kind == "absolute":
            return self.magnitude
        if self.kind == "relative":
            if predicted is None:
                return None
            return abs(predicted) * self.magnitude
        return None


_TOL_RE = re.compile(r"""
    ^\s*
    (?:(?P<name>[A-Za-z_][\w^]*)\s*)?           # p, R2, R^2 …
    (?P<op>±|\+/-|>=|<=|>|<)\s*
    (?P<num>[-+]?[0-9]*\.?[0-9]+(?:[eE][-+]?[0-9]+)?)
    \s*(?P<pct>%)?\s*$
""", re.VERBOSE)


def parse_tolerance(text: str) -> Tolerance:
    """A tolerance string → `Tolerance`. Unreadable means **raise** (never skip
    quietly).

    Skipping an unreadable tolerance lets that item pass with no verdict at all --
    the easiest way there is to neutralise the verification
    (master_plan §S2 failure modes).
    """
    m = _TOL_RE.match(str(text))
    if not m:
        raise ValueError(
            f"tolerance {text!r} cannot be read. Accepted forms: "
            f"'±1.5%', '±0.03', '>0.99', 'p>0.05', '<1e-3'")
    num = float(m.group("num"))
    op = m.group("op")
    if op in ("±", "+/-"):
        if num <= 0:
            raise ValueError(f"tolerance {text!r}: the half-width is 0 or less — "
                             f"every result would FAIL")
        return Tolerance("relative" if m.group("pct") else "absolute",
                         num / 100.0 if m.group("pct") else num, str(text))
    kind = "lower_bound" if op in (">", ">=") else "upper_bound"
    return Tolerance(kind, num / 100.0 if m.group("pct") else num, str(text))


# =============================================================================
# measurements
# =============================================================================
@dataclass
class Measurement:
    """One measurement. **A value with no `stat_err` cannot go in a conclusion**
    (CLAUDE.md).

    `stat_err=None` means "the error is unknown", and in that state no verdict is
    issued. A single-seed production run is stopped right here.
    """

    quantity: str
    value: float
    stat_err: float | None = None
    method: str = ""
    n_samples: int | None = None
    spread: float | None = None          # seed spread — a systematic-error indicator
    unit: str = ""

    def problems(self) -> list[str]:
        out = []
        if self.stat_err is None:
            out.append(f"{self.quantity}: no statistical error — either a "
                       f"single-seed run or the error was never computed. Cannot "
                       f"go in a conclusion")
        elif self.stat_err <= 0:
            out.append(f"{self.quantity}: stat_err = {self.stat_err} — a "
                       f"'measurement' that does not fluctuate is most likely an "
                       f"arithmetic identity "
                       f"(see guards.assert_statistic_fluctuates)")
        if not math.isfinite(self.value):
            out.append(f"{self.quantity}: the measured value is non-finite")
        return out


# =============================================================================
# one verdict row
# =============================================================================
@dataclass
class ValidationRow:
    quantity: str
    predicted: float | str
    measured: float
    tolerance: str
    verdict: str
    stat_err: float | None = None
    deviation: float | None = None        # absolute deviation
    deviation_rel: float | None = None    # relative deviation
    sigma: float | None = None            # deviation / SE
    design_power: float | None = None     # |predicted − competing| / SE
    samples_needed_for_3sigma: float | None = None
    cause_class: str = ""
    reason: str = ""
    note: str = ""
    flags: list[str] = field(default_factory=list)

    def problems(self) -> list[str]:
        out = []
        if self.verdict == FAIL and self.cause_class not in CAUSE_CLASSES:
            out.append(f"{self.quantity}: a FAIL with no cause class — "
                       f"one of {sorted(CAUSE_CLASSES)} is required (S7 gate)")
        if "significant_deviation_within_tolerance" in self.flags:
            out.append(
                f"{self.quantity}: **a PASS, but the deviation is "
                f"{self.sigma:.2f}σ** — inside the tolerance band yet a "
                f"statistically significant mismatch. Check whether the tolerance "
                f"is so wide that it neutralised the verification, or whether "
                f"there is a known bias that was left out of the prediction "
                f"(master_plan §S2 failure modes)")
        return out


def compare(item: PredictionItem, meas: Measurement, *,
            power_threshold: float = 1.0,
            significance_sigma: float = 3.0) -> ValidationRow:
    """One prediction vs one measurement → one **proposed** verdict row.

    Args:
        power_threshold: design power below this value gives `INCONCLUSIVE`.
            Default `1.0` — if the gap to the competing hypothesis is under `1σ`,
            that measurement cannot tell the two hypotheses apart. Raising it to
            `3.0` is more conservative.
        significance_sigma: flag a PASS whose deviation is this significant.
            ★ A wide tolerance turns any result into a PASS (§S2 failure modes).
            Landing inside the band and the prediction being right are two
            different statements — a PASS off by 3σ is a signal that "there is a
            bias the prediction left out".
    """
    tol = parse_tolerance(item.tolerance)
    se = meas.stat_err
    row = ValidationRow(quantity=item.quantity, predicted=item.value,
                        measured=meas.value, tolerance=item.tolerance,
                        verdict=INCONCLUSIVE, stat_err=se, note=item.note)

    # --- with an unknown error, no verdict is issued ---
    if se is None:
        row.reason = "no statistical error — an unbarred number is not usable"
        return row

    # --- one-sided bounds (R² > 0.99, p > 0.05) ---
    if tol.kind in ("lower_bound", "upper_bound"):
        bound = tol.magnitude
        row.deviation = meas.value - bound
        row.sigma = abs(row.deviation) / se if se > 0 else None
        ok = meas.value > bound if tol.kind == "lower_bound" else meas.value < bound
        if se > 0 and abs(meas.value - bound) < se:
            row.verdict = INCONCLUSIVE
            row.reason = (f"bound {bound:g}, gap {abs(row.deviation):.3g}, "
                          f"SE = {se:.3g} — indistinguishable from the bound")
        else:
            row.verdict = PASS if ok else FAIL
            row.reason = (f"measured {meas.value:.6g} "
                          f"{'>' if tol.kind == 'lower_bound' else '<'} {bound:g}"
                          + ("" if ok else " violated"))
        return row

    # --- two-sided band ---
    if not isinstance(item.value, (int, float)):
        row.reason = f"predicted value is not numeric ({item.value!r}) — no band"
        return row

    pred = float(item.value)
    half = tol.half_width(pred)
    row.deviation = meas.value - pred
    row.deviation_rel = row.deviation / pred if pred else None
    row.sigma = abs(row.deviation) / se if se > 0 else None

    # ① decidability — if SE exceeds the band, this measurement cannot decide
    #    this tolerance
    if half is not None and se > half:
        row.verdict = INCONCLUSIVE
        row.reason = (f"SE = {se:.4g} > tolerance half-width {half:.4g} — "
                      f"this measurement cannot decide this band. "
                      f"sample multiple needed ≈ {(se / half) ** 2:.3g}×")
    else:
        inside = half is not None and abs(row.deviation) <= half
        row.verdict = PASS if inside else FAIL
        row.reason = (f"deviation {row.deviation:+.4g}"
                      + (f" ({row.deviation_rel:+.3%})" if row.deviation_rel else "")
                      + f", band ±{half:.4g}")

    # ② design power — is it distinguishable from the competing hypothesis
    if item.competing_value is not None and se > 0:
        gap = abs(pred - float(item.competing_value))
        row.design_power = gap / se
        row.samples_needed_for_3sigma = (3.0 / row.design_power) ** 2 \
            if row.design_power > 0 else float("inf")
        if row.design_power < power_threshold:
            row.verdict = INCONCLUSIVE
            row.reason = (
                f"design power {row.design_power:.2f}σ < {power_threshold:g}σ — "
                f"indistinguishable from the competing hypothesis "
                f"({item.competing_value:g}). "
                f"{row.samples_needed_for_3sigma:.3g}x the samples needed to "
                f"resolve at 3σ. "
                f"★ a foreseen limit, not a failure")

    # ③ a significant mismatch hidden by a wide tolerance
    if (row.verdict == PASS and row.sigma is not None
            and row.sigma > significance_sigma):
        row.flags.append("significant_deviation_within_tolerance")
    return row


# =============================================================================
# report
# =============================================================================
@dataclass
class ValidationReport:
    rows: list[ValidationRow]
    seal: SealVerdict | None = None
    problems: list[str] = field(default_factory=list)
    sanity: list[ValidationRow] = field(default_factory=list)

    # --- tallies ---
    def count(self, verdict: str) -> int:
        return sum(1 for r in self.all_rows() if r.verdict == verdict)

    def all_rows(self) -> list[ValidationRow]:
        return [*self.rows, *self.sanity]

    @property
    def verdict_overall(self) -> str:
        if self.seal is not None and not self.seal.ok:
            return "SEAL_BROKEN"
        if self.count(FAIL):
            return "FAIL"
        if self.count(INCONCLUSIVE):
            return "PASS_WITH_INCONCLUSIVE"
        return PASS

    def yaml_block(self) -> str:
        """The verdict block. ★ `confirmed_by` is **always** `null` — code never
        fills it in."""
        return (
            "```yaml\n"
            f"verdict_overall: {self.verdict_overall}\n"
            "proposed_by: agent\n"
            "confirmed_by: null            # ← awaiting human confirmation\n"
            f"pass: {self.count(PASS)}\n"
            f"inconclusive: {self.count(INCONCLUSIVE)}\n"
            f"fail: {self.count(FAIL)}\n"
            "```")

    def table(self) -> str:
        rows = ["| quantity | predicted (sealed) | measured | tolerance | "
                "deviation | power | verdict |",
                "|---|---|---|---|---|---|---|"]
        mark = {PASS: "**PASS**", FAIL: "❌ **FAIL**",
                INCONCLUSIVE: "⚠ **INCONCLUSIVE**"}
        for r in self.all_rows():
            suffix = " ⚑" if r.flags else ""
            pred = f"`{r.predicted:.6g}`" if isinstance(r.predicted, (int, float)) \
                else f"`{r.predicted}`"
            meas = f"`{r.measured:.6g}`"
            if r.stat_err is not None:
                meas += f" ± `{r.stat_err:.3g}`"
            dev = ""
            if r.deviation_rel is not None:
                dev = f"{r.deviation_rel:+.2%}"
            elif r.deviation is not None:
                dev = f"{r.deviation:+.3g}"
            if r.sigma is not None:
                dev += f" ({r.sigma:.2f}σ)"
            power = f"`{r.design_power:.2f}σ`" if r.design_power is not None else "—"
            rows.append(f"| `{r.quantity}` | {pred} | {meas} | `{r.tolerance}` "
                        f"| {dev} | {power} | {mark.get(r.verdict, r.verdict)}{suffix} |")
        if any(r.flags for r in self.all_rows()):
            rows.append("")
            rows.append("⚑ = a PASS whose deviation is statistically significant "
                        "— see below.")
        return "\n".join(rows)

    def reasons(self) -> str:
        """Spell out the reasons for `INCONCLUSIVE` and `FAIL`. PASS is omitted."""
        out = []
        for r in self.all_rows():
            if r.verdict == PASS:
                continue
            head = f"- **`{r.quantity}`** — {r.verdict}"
            if r.cause_class:
                head += f" (`{r.cause_class}`)"
            out.append(f"{head}: {r.reason}")
        return "\n".join(out) if out else "_none — every item PASS._"

    def notes(self) -> str:
        """Notes attached to the prediction document. **Kept even on a PASS** --
        this is where limits like "that is not an independent check" live, and
        losing them overstates the conclusion."""
        out = [f"- **`{r.quantity}`**: {r.note}" for r in self.all_rows() if r.note]
        return "\n".join(out) if out else ""


def validate_run(prediction: Prediction, measurements: dict[str, Measurement], *,
                 rundir: RunDir | None = None, power_threshold: float = 1.0,
                 causes: dict[str, str] | None = None) -> ValidationReport:
    """Verify the seal → compare against the prediction → propose a verdict.

    ⚠ **With a broken seal, no verdict is produced.** Drawing a comparison table
      while the prediction may have been edited after seeing the result gives a
      table that looks like a verification and is not one.

    Args:
        measurements: `quantity` → the measured value. A prediction with no
            matching measurement is reported as a problem.
        causes: cause classes for FAIL items (supplied by the agent). Without
            them, the gate catches it.
    """
    problems: list[str] = list(prediction.problems())
    seal = verify_seal(rundir) if rundir is not None else None

    if seal is not None and not seal.ok:
        return ValidationReport(
            rows=[], seal=seal,
            problems=[f"★ seal violation — {seal.summary()}. "
                      f"The prediction may have been edited after the run, so "
                      f"**no comparison table is built** "
                      f"(master_plan §S7-1)."] + problems)

    rows: list[ValidationRow] = []
    for item in prediction.items:
        meas = measurements.get(item.quantity)
        if meas is None:
            problems.append(f"no measurement matches the prediction "
                            f"{item.quantity!r} — a sealed prediction cannot be "
                            f"quietly dropped")
            continue
        problems += meas.problems()
        row = compare(item, meas, power_threshold=power_threshold)
        if causes and item.quantity in causes:
            row.cause_class = causes[item.quantity]
        rows.append(row)
        problems += row.problems()

    extra = sorted(set(measurements) - {i.quantity for i in prediction.items})
    if extra:
        problems.append(f"unsealed measurements {extra} — an item added after the "
                        f"fact must be reported separately as an 'observation', "
                        f"not as a 'prediction comparison'")

    return ValidationReport(rows=rows, seal=seal, problems=problems)
