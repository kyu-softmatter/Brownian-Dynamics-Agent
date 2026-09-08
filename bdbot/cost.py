"""Block F — the cost gate, with a budget the user actually stated. LLM 0 lines.

`A3` says: *"If the estimate exceeds the budget, do not run. In a chat session,
**ask** rather than refuse. Either way, do not start without an estimate."*

An estimator existed (`cli.py cmd_run`, `cmd_calibrate`). What did not exist is a
place where the **budget** is recorded, so "exceeds the budget" had no left-hand
side. Measured: `trap-2d-5um`'s production run took **3037.5 s** wall at
**332.5 steps/s**, and nothing anywhere records whether anyone had 51 minutes.

## The two rules this module enforces

**1 · No estimate, no run.** `A3`, unchanged.

**2 · An INFEASIBLE verdict must carry what to give up, priced in the currency
that matters.** "Reduce N" is not a trade-off; *"N 1000 -> 400, statistical error
x1.6"* is. Cutting `N` to fit a clock changes the error bar, and by `A2` a number
without an error bar is not a result — so a silent reduction converts a run into
a non-result. This module therefore **proposes** and never applies.

## Why throughput is stored per case and never globally

`cli.py:628` already refuses to overwrite the global constant, because the kernel
differs per system: *"you must NOT change the global Lambda to this value — WCA
estimates would be optimistically wrong by {ratio}x."* Same rule here. A
throughput measured on a trap kernel does not predict a WCA kernel, and a
free-diffusion kernel (no forces, no neighbour list) is faster than both.
"""
from __future__ import annotations

import json
import math
from dataclasses import dataclass, field, asdict
from pathlib import Path

SCHEMA = "bdbot.cost/0.1"

VERDICTS = ("FEASIBLE", "TIGHT", "INFEASIBLE")

#: A run inside this factor of the budget is FEASIBLE but flagged. Chosen to match
#: `checks.MARGIN_WARN = 5.0`'s intent: report the ratio, not just the verdict.
TIGHT_MARGIN = 1.5

#: Measured throughputs, per case, with the kernel named. **Never a global.**
#: source: runs/<run_id>/metrics.json `steps_per_second`.
MEASURED_THROUGHPUT: dict = {
    "trap-2d-5um": {
        "steps_per_second": 332.5073748089591,
        "n_particles": 1000,
        "kernel": "md.force.Custom harmonic trap, no pair interaction, 2D",
        "source": "runs/trap-2d-5um__a5ef4f45d589/metrics.json",
        "wall_seconds": 3037.526612997055,
    },
}


@dataclass
class Estimate:
    """A cost estimate and the budget it is judged against.

    ⚠️ **`n_steps` is TOTAL steps, `n_eq + n_prod`.** Measured: the archived
    `trap-2d-5um` run recorded `wall_seconds = 3037.53` and
    `steps_per_second = 332.507375`, and `3037.53 x 332.507375 = 1,010,000`,
    which is `steps_done`, **not** `n_prod = 1,000,000`. So the throughput
    denominator includes equilibration. Estimating from `n_prod` alone
    under-predicts by 0.99 % here, and by more whenever equilibration is a larger
    fraction -- a relaxation case with `n_eq ~ n_prod` would be out by ~50 %.
    `total_steps()` is provided so a caller does not have to remember.
    """
    case: str
    n_steps: int
    n_particles: int
    steps_per_second: float
    basis: str                      # where steps_per_second came from
    user_budget_s: float | None = None
    urgency: str = ""               # free text from the user, e.g. "same-day"
    give_up: list = field(default_factory=list)

    @property
    def estimate_s(self) -> float:
        if self.steps_per_second <= 0:
            raise ValueError("steps_per_second must be > 0")
        scale = (self.n_particles / MEASURED_THROUGHPUT.get(
            self.case, {}).get("n_particles", self.n_particles))
        return self.n_steps / self.steps_per_second * max(scale, 1e-12)

    @property
    def margin(self) -> float:
        """budget / estimate. > 1 means it fits."""
        if self.user_budget_s is None:
            return math.nan
        e = self.estimate_s
        return math.inf if e <= 0 else self.user_budget_s / e

    @property
    def verdict(self) -> str:
        if self.user_budget_s is None:
            return "INFEASIBLE"          # A3: no budget is not a pass
        m = self.margin
        if m < 1.0:
            return "INFEASIBLE"
        return "TIGHT" if m < TIGHT_MARGIN else "FEASIBLE"

    @property
    def blockers(self) -> list:
        out = []
        if self.user_budget_s is None:
            out.append("user_budget_s is unset -- A3 requires a budget from the "
                       "user, not a default. Ask.")
        elif self.margin < 1.0:
            out.append(f"estimate {self.estimate_s:.0f} s exceeds the budget "
                       f"{self.user_budget_s:.0f} s by "
                       f"{self.estimate_s / self.user_budget_s:.2f}x")
        if self.verdict == "INFEASIBLE" and not self.give_up:
            out.append("INFEASIBLE with no `give_up` list -- a refusal must say "
                       "what could be traded away, priced in statistical error "
                       "rather than step count")
        return out

    def to_json(self) -> dict:
        d = {"schema": SCHEMA, "case": self.case,
             "n_steps": int(self.n_steps), "n_particles": int(self.n_particles),
             "steps_per_second": float(self.steps_per_second),
             "estimate_basis": self.basis,
             "estimate_s": round(self.estimate_s, 3),
             "user_budget_s": self.user_budget_s,
             "urgency": self.urgency,
             "margin": (None if math.isnan(self.margin) else round(self.margin, 4)),
             "verdict": self.verdict,
             "if_infeasible_give_up": list(self.give_up),
             "blockers": self.blockers}
        return d

    def write(self, outdir) -> Path:
        p = Path(outdir) / "cost.json"
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.to_json(), indent=1, ensure_ascii=False) + "\n")
        return p


def from_measured(case: str, n_steps: int, n_particles: int, *,
                  user_budget_s=None, urgency="", give_up=()) -> Estimate:
    """Build an estimate from a throughput measured on **this case**.

    Raises if the case has no measurement. That refusal is the point: guessing a
    throughput from a different kernel is how a 51-minute job gets started
    casually, which is the incident `A3` came from.
    """
    m = MEASURED_THROUGHPUT.get(case)
    if m is None:
        raise KeyError(
            f"no measured throughput for case '{case}'. Known: "
            f"{sorted(MEASURED_THROUGHPUT)}. Run a smoke test first and record "
            f"steps_per_second -- do NOT borrow another case's number, the "
            f"kernel differs (cli.py:628).")
    return Estimate(case=case, n_steps=int(n_steps), n_particles=int(n_particles),
                    steps_per_second=m["steps_per_second"],
                    basis=f"measured on {m['kernel']} at N={m['n_particles']}, "
                          f"{m['source']}",
                    user_budget_s=user_budget_s, urgency=urgency,
                    give_up=list(give_up))


def from_smoke(case: str, n_steps: int, n_particles: int, *,
               smoke_steps: int, smoke_wall_s: float,
               user_budget_s=None, urgency="", give_up=()) -> Estimate:
    """Build an estimate from a smoke run of the same kernel — the preferred route.

    Block H exists to produce exactly these two numbers.
    """
    if smoke_wall_s <= 0 or smoke_steps <= 0:
        raise ValueError("smoke_steps and smoke_wall_s must both be > 0")
    sps = smoke_steps / smoke_wall_s
    return Estimate(case=case, n_steps=int(n_steps), n_particles=int(n_particles),
                    steps_per_second=sps,
                    basis=f"smoke run of this kernel: {smoke_steps} steps in "
                          f"{smoke_wall_s:.3f} s = {sps:.1f} steps/s",
                    user_budget_s=user_budget_s, urgency=urgency,
                    give_up=list(give_up))


def total_steps(n_eq: int, n_prod: int) -> int:
    """`n_eq + n_prod` — what `steps_per_second` is measured over.

    Exists so the 0.99 % (and, for a relaxation case, ~50 %) error of using
    `n_prod` alone cannot be made by forgetting. See `Estimate`'s docstring.
    """
    if n_eq < 0 or n_prod < 0:
        raise ValueError("step counts must be >= 0")
    return int(n_eq) + int(n_prod)


def error_scaling(n_from: int, n_to: int) -> float:
    """Factor the statistical error changes by when N goes `n_from -> n_to`.

    `SE ~ 1/sqrt(N)`, so shrinking N to fit a clock **inflates** the error bar by
    this factor. The number a `give_up` entry must quote.
    """
    if n_from <= 0 or n_to <= 0:
        raise ValueError("N must be > 0")
    return math.sqrt(n_from / n_to)


def propose_give_up(est: Estimate, *, n_floor=100) -> list:
    """Suggest trades that would bring `est` inside its budget, each priced.

    Proposes only. Applying one is the user's decision — see the module docstring.
    """
    if est.user_budget_s is None or est.margin >= 1.0:
        return []
    need = est.estimate_s / est.user_budget_s          # factor too slow
    out = []
    n_to = max(n_floor, int(round(est.n_particles / need)))
    if n_to < est.n_particles:
        out.append(f"N {est.n_particles} -> {n_to} "
                   f"(statistical error x{error_scaling(est.n_particles, n_to):.2f})")
    steps_to = int(est.n_steps / need)
    out.append(f"n_steps {est.n_steps} -> {steps_to} "
               f"(loses {math.log10(est.n_steps / max(steps_to, 1)):.2f} decades "
               f"of the observable's time range)")
    out.append(f"or raise the budget to {est.estimate_s:.0f} s and say why this "
               f"system costs that much (A3 anti-pattern: raising it without a reason)")
    return out


def render(est: Estimate) -> str:
    L = ["=" * 78, f"cost gate — {est.case}", "=" * 78]
    L.append(f"  steps        {est.n_steps:,}   x   N = {est.n_particles}")
    L.append(f"  throughput   {est.steps_per_second:.1f} steps/s")
    L.append(f"               {est.basis}")
    e = est.estimate_s
    L.append(f"  estimate     {e:,.0f} s  =  {e / 60:.1f} min  =  {e / 3600:.2f} h")
    if est.user_budget_s is None:
        L.append("  budget       ** UNSET ** -- A3 requires one from the user")
    else:
        L.append(f"  budget       {est.user_budget_s:,.0f} s"
                 f"   (urgency: {est.urgency or 'unstated'})")
        L.append(f"  margin       {est.margin:.2f}x")
    if est.give_up:
        L.append("")
        L.append("  IF INFEASIBLE, WHAT TO GIVE UP (proposed, not applied)")
        for g in est.give_up:
            L.append(f"    - {g}")
    L.append("")
    L.append("=" * 78)
    L.append(f"VERDICT: {est.verdict}")
    for b in est.blockers:
        L.append(f"  ✗ {b}")
    L.append("=" * 78)
    return "\n".join(L)


def blocks(est: Estimate) -> bool:
    return est.verdict == "INFEASIBLE"


__all__ = ["SCHEMA", "VERDICTS", "TIGHT_MARGIN", "MEASURED_THROUGHPUT",
           "Estimate", "from_measured", "from_smoke", "error_scaling",
           "propose_give_up", "render", "blocks", "total_steps"]
