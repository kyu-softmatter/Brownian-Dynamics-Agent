"""Block H — the smoke profile. One convention, not eight bespoke flags. LLM 0 lines.

## Why this is shared and not per case

Measured 2026-09-02: `--smoke` existed in **5 of 8** case scripts
(`abp_rod_2d`, `chain_relax_2d_dlvo`, `soft_r3_2d`, `trap_2d_5um`,
`trap_drag_2d`) and was absent from `chain_bend_2d`, `chain_bend_dlvo_2d` and
`network_3d`. Each of the five implemented it inline, with its own choice of what
to shrink.

Adding three more inline versions would have made eight conventions. The knobs
differ per case — `--cycles` and `--samples` for the oscillatory ones, `--n` and
the stage times for `network` — so what is shared is not the *values* but the
**contract**:

  1. A smoke run must finish in **seconds**, not minutes. It answers *"does this
     path execute on this machine, and at what rate"* and nothing else.
  2. It must exercise the **same kernel** as production. Shrinking `N` and the
     step count is fine; switching off a force is not, because then the measured
     throughput does not predict the production run's — which is exactly the
     error `cli.py:628` refuses to make with the global constant.
  3. Its output is **two numbers**: steps executed and wall seconds. Block F
     (`cost.from_smoke`) consumes precisely those.
  4. It is **not** the production measurement and must never be reported as one
     (CLAUDE.md: *"cheap large-dt animations are fine but must be labelled as not
     being the production measurement"*).

## ⚠ Verification status

The five pre-existing profiles are transcribed from what those scripts already
did, so they are as verified as those scripts. **The three new ones
(`chain-bend-2d-oscill`, `chain-bend-2d-dlvo`, `network`) have never been
executed** — they need HOOMD, which is not available where they were written.
They are scaled-down versions of each case's own defaults, so the risk is low,
but `verified: False` says so rather than implying otherwise. Run
`--smoke` once on each and flip the flag.
"""
from __future__ import annotations

from dataclasses import dataclass, field

SCHEMA = "bdbot.smoke/0.1"

#: A smoke run above this wall time is not a smoke run. Advisory: the profile
#: cannot enforce it without running, but `check_wall` reports it afterwards.
SMOKE_WALL_BUDGET_S = 60.0


@dataclass
class Profile:
    """Per-case argument overrides that turn a production run into a smoke run.

    `overrides` maps `argparse` destination names to values. Unknown destinations
    raise in `apply()` rather than being silently ignored -- a typo'd override is
    a smoke run that quietly ran at production size.
    """
    case: str
    overrides: dict
    why: str
    verified: bool = False
    expect_wall_s: float | None = None


PROFILES: dict = {p.case: p for p in (
    # ── the five that already had --smoke, transcribed from those scripts ──
    Profile("trap-2d-5um", {"N": 200, "production_tau": 120, "equilibration_tau": 10},
            why="cases/trap_2d_5um.py:358 -- N 1000->200, production 2000->120 tau_k",
            verified=True, expect_wall_s=20.0),
    Profile("soft-r3-2d-A-sweep", {"smoke": True},
            why="the case script owns its own scaling (cases/soft_r3_2d.py:409)",
            verified=True),
    Profile("abp-rod-2d-run-flip", {"smoke": True},
            why="the case script owns its own scaling (cases/abp_rod_2d.py:509)",
            verified=True),
    Profile("chain-relax-2d-dlvo", {"smoke": True},
            why="the case script owns its own scaling (cases/chain_relax_2d_dlvo.py:762)",
            verified=True),
    Profile("trap-drag-2d-hex300", {"smoke": True},
            why="the case script owns its own scaling (cases/trap_drag_2d.py:1041)",
            verified=True),

    # ── the three that had none. NOT YET EXECUTED -- see the module docstring ──
    Profile("chain-bend-2d-oscill",
            {"cycles": 1.0, "samples": 100},
            why="1 drive cycle and 100 samples instead of 10 cycles / 2000 samples. "
                "The bond, angle and trap forces all stay ON, so the throughput "
                "still predicts production (contract 2). One cycle cannot measure "
                "K'(omega) -- it is not meant to.",
            verified=False),
    Profile("chain-bend-2d-dlvo",
            {"n": 5, "cycles": 1.0, "samples": 100, "eq_scale": 20.0},
            why="the shortest chain in system.yaml (n=5), 1 cycle, 100 samples, and "
                "equilibration 200->20 x tau_fast/dt. DLVO pair force, WCA core and "
                "traps all stay on.",
            verified=False),
    Profile("network",
            {"n": 64, "stage_tau": 1e-4, "agg_tau": 0.01, "post_tau": 0.01,
             "n_seeds": 1},
            why="N 512->64, each stage time /20, one seed. The 3D DLVO+WCA kernel "
                "and the BoxResize compression both stay on -- compression is the "
                "expensive part and dropping it would break contract 2.",
            verified=False),
)}


def apply(case: str, args) -> list:
    """Apply `case`'s smoke overrides to an `argparse.Namespace`, in place.

    Returns the list of `"dest: old -> new"` strings, so the caller can print
    exactly what was shrunk. Raises `KeyError` for an unknown case and
    `AttributeError` for an override whose destination the parser does not have --
    both loudly, because a silently ignored override means a "smoke" run at
    production size, which is worse than no smoke run at all.
    """
    p = PROFILES.get(case)
    if p is None:
        raise KeyError(f"no smoke profile for case {case!r}. Known: "
                       f"{sorted(PROFILES)}. Add one to bdbot/smoke.py rather "
                       f"than inlining a fourth convention.")
    changed = []
    for dest, val in p.overrides.items():
        if not hasattr(args, dest):
            raise AttributeError(
                f"smoke profile for {case!r} sets {dest!r}, which this parser does "
                f"not define. A silently dropped override gives a smoke run at "
                f"production size.")
        old = getattr(args, dest)
        setattr(args, dest, val)
        if old != val:
            changed.append(f"{dest}: {old} -> {val}")
    return changed


def banner(case: str, changed) -> str:
    p = PROFILES[case]
    L = ["── SMOKE RUN " + "─" * 62,
         "  NOT the production measurement. Answers only: does this path execute,",
         "  and at what rate (block F consumes steps + wall seconds).",
         f"  why these values: {p.why}"]
    if not p.verified:
        L.append("  ⚠ this profile has never been executed -- verified: False")
    for c in changed:
        L.append(f"    {c}")
    L.append("─" * 74)
    return "\n".join(L)


def check_wall(case: str, wall_s: float, budget_s: float = SMOKE_WALL_BUDGET_S) -> list:
    """Advisory notes after a smoke run. Empty means it behaved like a smoke run."""
    out = []
    p = PROFILES.get(case)
    if wall_s > budget_s:
        out.append(f"[warn] smoke took {wall_s:.1f} s, above the {budget_s:g} s "
                   f"budget -- shrink the profile further, or it is not answering "
                   f"'does this execute' quickly enough to be worth having")
    if p is not None and p.expect_wall_s and wall_s > 5.0 * p.expect_wall_s:
        out.append(f"[warn] {wall_s:.1f} s is more than 5x the recorded "
                   f"{p.expect_wall_s:g} s for this profile -- the machine or the "
                   f"kernel changed")
    return out


def unverified() -> list:
    """Profiles that have never been run. The honest coverage number."""
    return sorted(c for c, p in PROFILES.items() if not p.verified)


__all__ = ["SCHEMA", "SMOKE_WALL_BUDGET_S", "Profile", "PROFILES", "apply",
           "banner", "check_wall", "unverified"]
