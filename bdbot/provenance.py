"""출처 추적 — 마스터플랜 원칙 2 "모든 숫자는 출처를 갖는다".

1-A·1-B가 똑같이 `Provenanced` + YAML 노드 파서를 갖고 있었습니다.

tier: 0=직접입력/핸드북 · 1=문헌+검증 또는 확인된 관례 · 2=문헌 미검증 · 3=임의 가정
"""
from __future__ import annotations

from dataclasses import dataclass

from .units import Q


@dataclass
class Provenanced:
    """값 + 출처 + 신뢰등급. 값만 떠다니게 두지 않습니다."""

    value: object
    source: str
    tier: int

    def __repr__(self):
        try:
            return f"{self.value:~.4gP} [tier{self.tier}]"
        except (TypeError, ValueError):
            return f"{self.value!r} [tier{self.tier}]"


def load_node(node: dict) -> Provenanced:
    """`{value, unit, source, tier}` YAML 노드 → `Provenanced[Quantity]`."""
    return Provenanced(Q(node["value"], node["unit"]), node["source"], int(node["tier"]))


def tier_summary(items) -> dict:
    """tier별 개수. 검증기가 "tier 2 이하만으로 구성된 스펙"을 잡을 때 씁니다."""
    out: dict[int, int] = {}
    for p in items:
        out[p.tier] = out.get(p.tier, 0) + 1
    return dict(sorted(out.items()))


__all__ = ["Provenanced", "load_node", "tier_summary"]
