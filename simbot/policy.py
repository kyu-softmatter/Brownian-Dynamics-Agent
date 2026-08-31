"""The `config/run_policy.yaml` loader. 0 lines of LLM.

**There is exactly one place a human writes: `overrides:`** (run_policy.yaml §8).
This module's only job is to merge that block deeply over the rest. Merge it
shallowly and an attempt to change just `tiers.production.N` wipes out all of
`tiers.production` — quietly running a different system.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

import yaml

from .io import REPO_ROOT

DEFAULT_POLICY_PATH = REPO_ROOT / "config" / "run_policy.yaml"


def deep_merge(base: dict, over: dict) -> dict:
    """Merge `over` recursively over `base`. Descends into dicts only; a list is
    replaced wholesale.

    Not merging lists is deliberate -- someone writing
    `required_tier_ladder: [pilot]` to **shorten** the ladder must not end up with
    `[smoke, pilot, explore, pilot]`.
    """
    out = copy.deepcopy(base)
    for k, v in over.items():
        if isinstance(v, dict) and isinstance(out.get(k), dict):
            out[k] = deep_merge(out[k], v)
        else:
            out[k] = copy.deepcopy(v)
    return out


@dataclass
class Policy:
    """The policy after merging. `raw` is the whole merged result."""

    raw: dict
    overridden_paths: list[str]

    # --- lookup ---
    def get(self, dotted: str, default=None):
        cur = self.raw
        for part in dotted.split("."):
            if not isinstance(cur, dict) or part not in cur:
                return default
            cur = cur[part]
        return cur

    @property
    def timestep(self) -> dict:
        return self.raw.get("timestep", {})

    @property
    def seeds_default(self) -> int:
        return int(self.get("seeds.default", 4))

    @property
    def wall_budget_s(self) -> float:
        return float(self.get("budget.wall_time_per_run_s", 600))

    def concurrency(self, mode: str = "default") -> int:
        """Concurrency. Cannot exceed `hard_max` (a regression was measured at k=12)."""
        k = int(self.get(f"concurrency.{mode}", self.get("concurrency.default", 8)))
        return min(k, int(self.get("concurrency.hard_max", 10)))

    def tier(self, name: str) -> dict:
        t = self.get(f"tiers.{name}")
        if t is None:
            raise KeyError(f"tier {name!r} is not in the policy. "
                           f"available: {sorted(self.raw.get('tiers', {}))}")
        return t

    def tier_ladder(self) -> list[str]:
        return list(self.get("budget.required_tier_ladder", ["smoke", "pilot", "explore"]))

    def efficiency(self, k: int) -> float:
        """Per-process efficiency at concurrency k. A step function over the
        measured table."""
        table = {int(kk): float(vv)
                 for kk, vv in self.get("hardware.efficiency_by_k", {}).items()}
        if not table:
            return 1.0
        if k in table:
            return table[k]
        below = [kk for kk in table if kk <= k]
        return table[max(below)] if below else table[min(table)]


def _find_numeric_strings(node, prefix: str = "") -> list[str]:
    """Paths of values that should read as numbers but parsed as strings.

    ★ YAML 1.1 does not read an exponent as a float unless it carries a sign --
      `6.3e6` is the **string** `'6.3e6'`. That actually happened to
      `throughput_particle_steps_per_s`, and the cost estimate was quietly being
      handed a string (found 2026-07-28). The sign has to be written:
      `6.3e+6`.
    """
    bad: list[str] = []
    items = node.items() if isinstance(node, dict) else enumerate(node)
    for k, v in items:
        key = f"{prefix}.{k}" if prefix else str(k)
        if isinstance(v, (dict, list)):
            bad += _find_numeric_strings(v, key)
        elif isinstance(v, str):
            try:
                float(v)
            except ValueError:
                continue
            bad.append(f"{key} = {v!r}")
    return bad


def load_policy(path: str | Path | None = None) -> Policy:
    """Read the policy and apply `overrides:`. Returns the overridden paths too.

    Why it returns them: the report has to mark "a human chose this value" for the
    `params` command's ⚠ (a default nobody picked) to mean anything.
    """
    p = Path(path) if path else DEFAULT_POLICY_PATH
    raw = yaml.safe_load(p.read_text(encoding="utf-8")) or {}

    bad = _find_numeric_strings(raw)
    if bad:
        raise ValueError(
            f"{p}: a number parsed as a string — {bad}.\n"
            f"  YAML 1.1 does not read an exponent as a float without a sign: "
            f"`6.3e6` → string, `6.3e+6` → float.\n"
            f"  Left as is, the cost estimate is quietly wrong.")

    over = raw.pop("overrides", None) or {}
    over = {k: v for k, v in over.items() if not k.startswith("_")}  # `_why` is a note
    paths: list[str] = []

    def walk(d: dict, prefix: str = "") -> None:
        for k, v in d.items():
            key = f"{prefix}.{k}" if prefix else str(k)
            if isinstance(v, dict):
                walk(v, key)
            else:
                paths.append(key)

    walk(over)
    return Policy(raw=deep_merge(raw, over), overridden_paths=sorted(paths))
