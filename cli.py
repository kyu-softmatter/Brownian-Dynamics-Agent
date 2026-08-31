#!/usr/bin/env python
"""Root shim for the S5-S8 command surface. The implementation is
[`simbot/cli.py`](simbot/cli.py).

    python cli.py run <spec.yaml>            one spec -> REPORT.md
    python cli.py resume runs/<id>           take over a dead run
    python cli.py converge <spec.yaml>       is the answer the same when dt, N and the initial condition are shaken
    python cli.py params [--path runs]       the parameters of several runs, side by side
    python cli.py calibrate                  measure this machine's throughput

**Why this file is now a shim (2026-08-29).** This repository had **two CLIs** —
this one (S5-S8, prediction sealing, `INCONCLUSIVE`) and `bdbot/cli.py` (L0-L4,
intake, non-dimensionalization, health). That was the "two engines" seam in
[docs/00-merge-decisions.md](docs/00-merge-decisions.md) section 5: a case run
through `bdbot` got a health verdict but **no sealed prediction**, and a spec run
through here got sealing but **could not use bdbot's cases**.

The merge was mostly additive, because the two command sets turned out to be
almost disjoint — the only collision was `run`, which means different things in
the two halves. So the implementation moved into the package as `simbot.cli`, and
`bdbot.cli` exposes it under `pipeline`:

    python -m bdbot.cli pipeline run <spec.yaml>     one entry point
    python cli.py run <spec.yaml>                    this shim, byte-identical behaviour
    python -m simbot.cli run <spec.yaml>             the module directly

**This file is kept rather than deleted, and not for politeness.** Three things
depend on the name: `tests/test_cli_session.py` does `import cli` and reaches for
`cli.main`, `cli._runner_for`, `cli.build_parser`, `cli.BASELINE_KERNEL` and
`cli.CALIBRATE_KERNEL`; `tests/test_agent_layer.py` asserts the string
`"cli.py run"` appears in the pipeline skill; and four `.claude/skills/` documents
tell the agent to invoke `cli.py run <spec.yaml>`. A shim keeps all three true
while there is still exactly one implementation.
    ⚠ Re-export everything public, not just what today's tests touch. A shim that
      covers only the current call sites is a shim that breaks on the next one.
"""
from __future__ import annotations

import sys
from pathlib import Path

REPO = Path(__file__).resolve().parent
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

from simbot import cli as _impl                                        # noqa: E402
from simbot.cli import (BASELINE_KERNEL, CALIBRATE_KERNEL,              # noqa: E402,F401
                        RUNNERS, build_parser, cmd_calibrate, cmd_converge,
                        cmd_params, cmd_resume, cmd_run, main,
                        measure_trap_batch)

#  Anything public in `simbot.cli` that is not named above is still reachable as
#  `cli.<name>` -- so the shim cannot go stale when a command is added there.
#  (`_runner_for`, `_trap_configs`, `_print_cost`, ... are private-by-underscore but
#  `tests/test_cli_session.py` uses `_runner_for`, so they must resolve too.)
def __getattr__(name: str):
    try:
        return getattr(_impl, name)
    except AttributeError:
        raise AttributeError(
            f"module 'cli' has no attribute {name!r} — the implementation is "
            f"simbot.cli, and it has no such name either") from None


def __dir__():
    return sorted(set(dir(_impl)) | set(globals()))


if __name__ == "__main__":
    raise SystemExit(main())
