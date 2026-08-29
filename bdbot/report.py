"""`DimensionlessReport` 렌더러 — 사람 확인 #3의 대상 (마스터플랜 §6.5).

1-A·1-B가 똑같은 골격을 갖고 있었습니다:
  헤더 → 기준 스케일 + 근거 → INPUT → DERIVED → SCALE LEDGER →
  DIMENSIONLESS GROUPS → SEPARATION CHECKS → RUN PLAN → VERDICT

**공통인 것**: 프레임, 원장, 무차원수, 검사, 판정 (아래 구현).
**케이스마다 다른 것**: INPUT/DERIVED/RUN PLAN에 무엇을 적을지 → 케이스가 줄 목록으로 넘김.

폭은 여기서 하나로 통일했습니다 (1-A는 84, 1-B는 88을 쓰고 있었음).
숫자 결과가 아니라 렌더링이라 통일해도 결과는 바뀌지 않습니다.
"""
from __future__ import annotations

from .checks import verdict as _verdict

W = 88            # 전체 폭
LEDGER_W = 38     # 원장 항목 이름
GROUP_W = 44      # 무차원수 이름
CHECK_W = 42      # 검사 이름


def render(*, title, ref, ledger, groups, checks,
           input_lines=(), derived_lines=(), run_plan_lines=()) -> tuple[str, str]:
    """리포트 텍스트와 판정을 반환. 판정이 FAIL이면 호출부가 실행을 거부합니다."""
    out: list[str] = []
    w = out.append

    w("=" * W)
    w(title)
    w("=" * W)
    w(f"기준 스케일: length={ref['length'][0]}  energy={ref['energy'][0]}  time={ref['time'][0]}"
      + (f"   [strategy: {ref['strategy']}]" if ref.get("strategy") else ""))
    w(f"  근거: {ref.get('rationale', ledger.rationale)}")

    if input_lines:
        w("")
        w("INPUT (차원 있음, SI)")
        for ln in input_lines:
            w(ln)

    if derived_lines:
        w("")
        w("DERIVED")
        for ln in derived_lines:
            w(ln)

    w("")
    w("SCALE LEDGER  (작은 것부터)")
    for cat_name, cat in ledger.categories():
        if not cat:
            continue
        w(f"  {cat_name}")
        for nm, v in ledger.sorted_items(cat):
            w(f"    {nm:<{LEDGER_W}} {v.to_compact():~.4gP}")

    w("")
    w("DIMENSIONLESS GROUPS  (두 스케일의 비)")
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
        w(f"  ✗ 하드 위반: {c.name.strip()} = {c.value:.3g} (기준 {c.limit:g}) — {c.note}")
    for c in soft_fail:
        w(f"  ⚠ 경고: {c.name.strip()} = {c.value:.3g} (기준 {c.limit:g}) — {c.note}")
    for c in tight:
        w(f"  ⚠ 여유 부족: {c.name.strip()} ({c.margin:.1f}×)")
    w("=" * W)
    return "\n".join(out), v


def kv(key: str, value: str, tier=None, source: str = "", key_w: int = 6, val_w: int = 18) -> str:
    """INPUT 한 줄. 값 뒤에 tier와 출처를 붙입니다 (원칙 2)."""
    s = f"  {key:<{key_w}} = {value:<{val_w}}"
    if tier is not None:
        s += f" [tier {tier}]"
    if source:
        s += f" {source}"
    return s


__all__ = ["render", "kv", "W"]
