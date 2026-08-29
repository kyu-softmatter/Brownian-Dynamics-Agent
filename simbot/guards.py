"""Runtime guards -- catches a simulation being silently wrong.

**The implementations moved to [`bdbot/health.py`](../bdbot/health.py) section 0
on 2026-08-29.** This module is now a re-export so existing imports keep working.

Why they moved rather than one side being deleted: both packages had runtime
guards and **neither was a superset.**

    bdbot/health.py   minimum image in the displacement measure; the force-driven
                      vs thermal displacement split; the "worst value over the
                      whole run" convention (measured: peak force 1062.9 against
                      244.2 kT/sigma for the last sample, a factor of 4.4)
    simbot/guards.py  the configurational thermometer; the bond-length check; the
                      does-it-actually-fluctuate assertion; the uniform-noise
                      sqrt(3) bound

So the merge is a union, hosted in `bdbot.health` because that is the L4 layer the
engine calls and `simbot` may depend on `bdbot` but not the reverse.

The organising principle is unchanged and is the reason this module existed:
**a guard has to watch a quantity that CAN drift systematically.** HOOMD
`Brownian`'s kinetic temperature is redrawn from the target distribution every
step, so it cannot drift and cannot be a guard. The configurational temperature
takes that seat.
Basis: knowledge/wiki/findings/hoomd-brownian-scheme-and-noise.md
"""
from __future__ import annotations

from bdbot.health import (DisplacementReport, assert_statistic_fluctuates,
                          check_bond_lengths, check_finite, check_inside_box,
                          check_step_displacements, configurational_temperature)

__all__ = ["configurational_temperature", "DisplacementReport",
           "check_step_displacements", "check_finite", "check_inside_box",
           "check_bond_lengths", "assert_statistic_fluctuates"]
