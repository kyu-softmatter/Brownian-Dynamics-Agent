"""S8 — generating `REPORT.md`. 0 lines of LLM.

**The report is an audit record, not a summary.** A person must be able to read
this one file and answer "can I believe this conclusion". So it does not carry only
the good news:

- the seal status (**at the very top** if it is broken)
- the `INCONCLUSIVE` items and why
- reproducibility (`git_dirty` included)
- gates not yet decided
- `confirmed_by: null`

What `simbot` cannot produce — the answer to the question, a cause hypothesis, the
next experiment — is **quoted** from the agent-written `08_conclusion.md`. It is
never invented.
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

from .io import RUN_LAYOUT, RunDir, verify_seal
from .nondim import ReducedSpec, nondim_table, reduce_spec, roundtrip_errors
from .spec import SpecReport, SystemSpec, validate as validate_spec
from .validate import ValidationReport

ROUNDTRIP_GATE = 1e-12


@dataclass
class ReportInputs:
    """What the report needs. Leave a missing piece `None` and that section renders
    as 'none'."""

    spec: SystemSpec | None = None
    spec_report: SpecReport | None = None
    reduced: ReducedSpec | None = None
    validation: ValidationReport | None = None
    manifest: dict | None = None
    figures: dict[str, str] | None = None      # filename → caption
    wall_s: float | None = None
    n_runs: int | None = None


def _fmt(x: float | None, spec: str = ".6g") -> str:
    return "—" if x is None else f"`{x:{spec}}`"


# =============================================================================
# per-section renderers — each is tested independently
# =============================================================================
def seal_section(rundir: RunDir) -> str:
    v = verify_seal(rundir)
    if v.ok:
        return (f"**seal check** ✅ {v.summary()}\n\n"
                f"Verification command (works without this code):\n\n"
                f"```bash\nshasum -a 256 -c {RUN_LAYOUT['seal']}\n```")
    return ("> ## ⛔ seal violation\n>\n"
            f"> {v.summary()}\n>\n"
            "> **The prediction may have been edited after the run.** The comparison "
            "table below must not be read as a verification.\n"
            "> Basis: master_plan §S7-1.")


def reproducibility_section(manifest: dict | None) -> str:
    if not manifest:
        return "_no manifest — reproducibility cannot be claimed._"
    rows = ["| item | value |", "|---|---|"]
    for key, label in (("run_id", "run_id"), ("spec_hash", "spec hash"),
                       ("code_hash", "code hash"), ("git_rev", "git commit"),
                       ("env_hash", "env hash"), ("seed", "seed")):
        if key in manifest:
            rows.append(f"| {label} | `{manifest[key]}` |")
    env = manifest.get("env", {})
    for k in ("hoomd", "numpy", "scipy", "python"):
        if k in env:
            rows.append(f"| {k} | `{env[k]}` |")

    dirty = manifest.get("git_dirty")
    if dirty is True:
        rows.append("| **git state** | ⚠️ **uncommitted changes present** — "
                    "this run does not reproduce from `git_rev` alone |")
    elif dirty is False:
        rows.append("| git state | ✅ clean — `git_rev` pins the code |")
    else:
        rows.append("| git state | ❔ undecidable |")
    return "\n".join(rows)


def gates_section(spec_report: SpecReport | None) -> str:
    if spec_report is None:
        return "_no S3 check results._"
    out = [spec_report.table()]
    deferred = spec_report.deferred()
    if deferred:
        out.append("")
        out.append(f"⏳ **{len(deferred)} gates S7 has to decide** — "
                   f"quantities S3 cannot compute: "
                   + ", ".join(f"`{c.name}`" for c in deferred))
    if spec_report.problems:
        out.append("")
        out.append("### ⚠️ convention violation")
        out += [f"- {p}" for p in spec_report.problems]
    return "\n".join(out)


def nondim_section(spec: SystemSpec | None, reduced: ReducedSpec | None) -> str:
    if spec is None or reduced is None:
        return "_no S4 non-dimensionalization results._"
    errs = roundtrip_errors(spec, reduced)
    worst = max(errs.values()) if errs else 0.0
    gate = "✅ pass" if worst < ROUNDTRIP_GATE else "❌ **violation**"
    out = [f"reference scales: **{reduced.scales.origin}**", "",
           nondim_table(spec, reduced), "",
           f"**round-trip error** max `{worst:.2e}` "
           f"(gate `< {ROUNDTRIP_GATE:g}`) — {gate}", "",
           f"`Δt*` = `{reduced.dt_star:.6g}` · dominant constraint "
           f"**{reduced.dt_dominant}**"]
    if reduced.logged:
        out.append("")
        out.append("for the record (not a gate — used only when comparing to "
                   "another paper):")
        out += [f"- `{k}` = `{v:.4g}`" for k, v in reduced.logged.items()]
    if reduced.groups:
        out.append("")
        out.append("| dimensionless group | value |")
        out.append("|---|---|")
        out += [f"| `{k}` | `{v:.6g}` |" for k, v in reduced.groups.items()]
    return "\n".join(out)


def validation_section(validation: ValidationReport | None) -> str:
    if validation is None:
        return "_no S7 verdict._"
    out = [validation.table(), ""]
    n_pass = validation.count("PASS")
    n_inc = validation.count("INCONCLUSIVE")
    n_fail = validation.count("FAIL")
    out.append(f"**{len(validation.all_rows())} items: {n_pass} PASS · "
               f"{n_inc} INCONCLUSIVE · {n_fail} FAIL.**")
    if n_inc or n_fail:
        out += ["", "### items that are not PASS", validation.reasons()]
    if n_inc:
        out += ["", "> `INCONCLUSIVE` is not a failure — it is the fact that the "
                    "statistical error exceeds the tolerance, so **no verdict is "
                    "possible**. The sample multiple needed is above."]
    notes = validation.notes()
    if notes:
        out += ["", "### notes (PASS items included)", notes,
                "", "> These are the limits written into the prediction document. "
                    "A PASS does not make them go away — this is where the 'that "
                    "is not an independent check' kind lives, and dropping it "
                    "overstates the conclusion."]
    if validation.problems:
        out += ["", "### ⚠️ problems with the verification procedure",
                *[f"- {p}" for p in validation.problems]]
    return "\n".join(out)


def figures_section(rundir: RunDir, figures: dict[str, str] | None) -> str:
    """An uncaptioned figure is not accepted as an artefact (master_plan §S6 gate)."""
    files = sorted(p.name for p in rundir.figs.glob("*.png")) if rundir.figs.exists() \
        else []
    if not files:
        return "_no figures._"
    caps = figures or {}
    out = []
    for f in files:
        cap = caps.get(f)
        if cap:
            # the alt text is cut to one line — a multi-line caption inside
            # `![...]` breaks the markdown image syntax. The full caption goes
            # below the figure.
            alt = " ".join(cap.split())
            if len(alt) > 90:
                alt = alt[:89].rstrip() + "…"
            out.append(f"### {f}\n\n![{alt}](figs/{f})\n\n{cap}")
        else:
            out.append(f"### {f}\n\n![{f}](figs/{f})\n\n"
                       f"⚠️ **no caption** — what the figure is meant to show has "
                       f"to be written down (§S6 gate).")
    return "\n\n".join(out)


def cost_section(wall_s: float | None, n_runs: int | None) -> str:
    if wall_s is None:
        return "_no compute-cost record._"
    parts = [f"**total compute time `{wall_s:.1f} s`**"]
    if n_runs:
        parts.append(f"{n_runs} runs")
        parts.append(f"mean per run `{wall_s / n_runs:.2f} s`")
    return " · ".join(parts)


def _excerpt(rundir: RunDir, stage: str, title: str) -> str:
    """**Quote** the agent-written document. If absent, say so — never invent it."""
    if not rundir.exists(stage):
        return f"_`{RUN_LAYOUT[stage]}` is missing — the {title} has to be " \
               f"written by the agent._"
    return f"→ [`{RUN_LAYOUT[stage]}`]({RUN_LAYOUT[stage]})"


# =============================================================================
# assembly
# =============================================================================
def render(rundir: RunDir, inputs: ReportInputs) -> str:
    spec = inputs.spec
    sr = inputs.spec_report
    if sr is None and spec is not None:
        sr = validate_spec(spec)
    reduced = inputs.reduced
    if reduced is None and spec is not None:
        try:
            reduced = reduce_spec(spec)
        except Exception:                      # card with no registered scales etc.
            reduced = None

    v = inputs.validation
    headline = v.verdict_overall if v is not None else "no verdict"

    parts = [
        f"# REPORT — `{rundir.run_id}`",
        "",
        f"**verdict (proposed)** `{headline}` · **awaiting confirmation** "
        f"(`confirmed_by: null`)",
    ]
    if spec is not None:
        parts += ["", f"**question** {spec.question.strip()}",
                  f"**card** [`{spec.card}`]"
                  f"(../../knowledge/wiki/systems/{spec.card}.md)"]
    parts += ["", cost_section(inputs.wall_s, inputs.n_runs), "",
              seal_section(rundir), "", "---", ""]

    # if the seal is broken the comparison table does not go in the report
    seal_ok = verify_seal(rundir).ok

    parts += ["## 1. Verdict summary", ""]
    if not seal_ok:
        parts += ["The comparison table was **not generated** — because the seal "
                  "is broken.", ""]
    else:
        parts += [validation_section(v), ""]
    if v is not None:
        parts += [v.yaml_block(), "",
                  "> The verdict is a **proposal**. Until a human fills in "
                  "`confirmed_by` it does not enter the benchmark-ledger tally "
                  "(CLAUDE.md §verdicts).", ""]

    parts += ["---", "", "## 2. System specification and gates (S3)", "",
              gates_section(sr), "",
              "---", "", "## 3. Non-dimensionalization (S4)", "",
              nondim_section(spec, reduced), "",
              "---", "", "## 4. Figures (S6)", "",
              figures_section(rundir, inputs.figures), "",
              "---", "", "## 5. Reproducibility", "",
              reproducibility_section(inputs.manifest), "",
              "---", "", "## 6. Agent-written documents", "",
              f"- intake (S1) {_excerpt(rundir, 'intake', 'intake')}",
              f"- prediction (S2, sealed) "
              f"{_excerpt(rundir, 'prediction', 'prediction')}",
              f"- verification narrative (S7) "
              f"{_excerpt(rundir, 'validation', 'verification narrative')}",
              f"- conclusion (S8) {_excerpt(rundir, 'conclusion', 'conclusion')}",
              "",
              "> `simbot` produces the numbers; the answer to the question, the "
              "cause hypothesis and the next experiment are written by the agent. "
              "This report **quotes** those documents and does not write them "
              "instead.",
              ""]
    return "\n".join(parts)


def write_report(rundir: RunDir, inputs: ReportInputs) -> Path:
    return rundir.write("report", render(rundir, inputs))
