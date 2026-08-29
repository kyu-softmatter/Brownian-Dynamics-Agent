"""`config/run_policy.yaml` 로더. LLM 0줄.

**사람이 쓰는 곳은 `overrides:` 한 곳이다** (run_policy.yaml §8). 이 모듈의 일은
그 블록을 나머지 위에 깊게 덮어쓰는 것뿐이다. 얕게 덮어쓰면 `tiers.production.N`
하나만 바꾸려다 `tiers.production` 전체가 사라진다 — 조용히 다른 계를 돌리게 된다.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass
from pathlib import Path

import yaml

from .io import REPO_ROOT

DEFAULT_POLICY_PATH = REPO_ROOT / "config" / "run_policy.yaml"


def deep_merge(base: dict, over: dict) -> dict:
    """`over` 를 `base` 위에 재귀 병합. dict 만 파고들고 리스트는 통째로 교체한다.

    리스트를 병합하지 않는 것은 의도다 — `required_tier_ladder: [pilot]` 로
    사다리를 **줄이려는** 사람이 `[smoke, pilot, explore, pilot]` 을 얻으면 안 된다.
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
    """병합이 끝난 정책. `raw` 는 병합 결과 전체."""

    raw: dict
    overridden_paths: list[str]

    # --- 조회 ---
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
        """동시 실행 수. `hard_max` 를 넘지 못한다 (k=12 에서 회귀 측정됨)."""
        k = int(self.get(f"concurrency.{mode}", self.get("concurrency.default", 8)))
        return min(k, int(self.get("concurrency.hard_max", 10)))

    def tier(self, name: str) -> dict:
        t = self.get(f"tiers.{name}")
        if t is None:
            raise KeyError(f"티어 {name!r} 가 정책에 없다. "
                           f"있는 것: {sorted(self.raw.get('tiers', {}))}")
        return t

    def tier_ladder(self) -> list[str]:
        return list(self.get("budget.required_tier_ladder", ["smoke", "pilot", "explore"]))

    def efficiency(self, k: int) -> float:
        """동시 실행 k 개일 때 프로세스당 효율. 측정된 표의 계단 함수."""
        table = {int(kk): float(vv)
                 for kk, vv in self.get("hardware.efficiency_by_k", {}).items()}
        if not table:
            return 1.0
        if k in table:
            return table[k]
        below = [kk for kk in table if kk <= k]
        return table[max(below)] if below else table[min(table)]


def _find_numeric_strings(node, prefix: str = "") -> list[str]:
    """숫자로 읽혀야 하는데 문자열로 파싱된 값들의 경로.

    ★ YAML 1.1 은 지수부에 부호가 없으면 float 으로 읽지 않는다 — `6.3e6` 은
      **문자열** `'6.3e6'` 이다. 실제로 `throughput_particle_steps_per_s` 가 그랬고,
      비용 추정이 조용히 문자열을 받고 있었다 (2026-07-28 발견).
      `6.3e+6` 처럼 부호를 적어야 한다.
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
    """정책을 읽고 `overrides:` 를 적용한다. 덮어쓴 경로를 함께 돌려준다.

    덮어쓴 경로를 반환하는 이유: 리포트에 "이 값은 사람이 정했다"를 표시해야
    `params` 명령의 ⚠(아무도 안 고른 기본값) 표시가 의미를 가진다.
    """
    p = Path(path) if path else DEFAULT_POLICY_PATH
    raw = yaml.safe_load(p.read_text(encoding="utf-8")) or {}

    bad = _find_numeric_strings(raw)
    if bad:
        raise ValueError(
            f"{p}: 숫자가 문자열로 파싱됐다 — {bad}.\n"
            f"  YAML 1.1 은 지수부에 부호가 없으면 float 으로 읽지 않는다: "
            f"`6.3e6` → 문자열, `6.3e+6` → float.\n"
            f"  이대로 두면 비용 추정이 조용히 틀린다.")

    over = raw.pop("overrides", None) or {}
    over = {k: v for k, v in over.items() if not k.startswith("_")}  # `_why` 는 주석
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
