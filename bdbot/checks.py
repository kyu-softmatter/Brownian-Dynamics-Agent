"""분리 검사 — 마스터플랜 원칙 3 ④, skill `bd-physics` §4.

1-A·1-B가 똑같이 쓴 것: `Check(kind, name, value, limit, op)` + `ok`/`margin`,
그리고 검사를 **모델 / 적분 / 기하 / 통계** 네 종류로 분류하는 규약.

`hard` 플래그는 1-B에서 필요해졌습니다 — 1-A는 모든 검사를 통과해서 하드/소프트
구분이 드러나지 않았습니다. bd-physics §4는 모델·적분·기하를 ❌ 하드,
통계·유한크기를 ⚠ 경고로 규정합니다.

`dt`는 **가장 빠른 물리 시간척도의 1% 이하**로 잡습니다. 두 케이스 모두 그 시간척도가
`γ/(국소 강성)` 형태였습니다 — 트랩은 `γ/k`, 소프트 페어는 `γ/U''(r_min)`.
"""
from __future__ import annotations

from dataclasses import dataclass

GATE = 1e-2                      # 모든 분리 검사 임계 (bd-physics §4). dt/τ=1e-2 ⟺ 편향 0.5%
MARGIN_WARN = 5.0                # 여유가 이보다 작으면 "파라미터 상향 여지 없음" 경고
SOFT_KINDS = frozenset({"statistics", "finite-size"})
# ★ 부동소수 경계 허용. **구성상 한계와 같게 잡은** 검사가 우연히 실패하는 것을 막습니다.
#   실제로 물렸습니다: chain-bend 는 T_obs ≡ 100·(2π/ω) 로 정의하는데
#   `(100x)/x = 99.99999999999999` 가 나와 "관측 주기 수 ≥ 100" 이 거짓 경고를 냈습니다.
#   같은 스펙의 `dt/τ_fast ≤ 1e-2` 는 우연히 정확히 떨어져 통과했습니다 — 그 우연에
#   기대면 안 됩니다. `bdbot.nondim.RATIO_RTOL` 과 같은 취지입니다.
CMP_RTOL = 1e-12


@dataclass
class Check:
    kind: str                    # 모델 / 적분 / 기하 / 통계
    name: str
    value: float
    limit: float
    op: str = "<="
    note: str = ""
    hard: bool = True

    @property
    def ok(self) -> bool:
        tol = CMP_RTOL * max(abs(self.value), abs(self.limit))
        return (self.value <= self.limit + tol if self.op == "<="
                else self.value >= self.limit - tol)

    @property
    def margin(self) -> float:
        if self.op == "<=":
            return self.limit / self.value if self.value else float("inf")
        return self.value / self.limit if self.limit else float("inf")

    def as_dict(self, phase: str = "design") -> dict:
        # `op` 는 L3 스펙을 되읽을 때 필요합니다 — 없으면 ">=" 검사(관측창 등)가
        # "<=" 로 복원되어 판정이 뒤집힙니다 (bdbot.nondim.load).
        return {"kind": self.kind, "name": self.name.strip(), "value": self.value,
                "limit": self.limit, "op": self.op, "ok": bool(self.ok),
                "margin": self.margin, "hard": bool(self.hard), "phase": phase}


def soft(kind: str) -> bool:
    """이 종류가 경고인가? (bd-physics §4의 하드/소프트 규정)"""
    return kind in SOFT_KINDS


def verdict(checks) -> tuple[str, list, list, list]:
    """(verdict, hard failures, soft failures, thin margins).

    A single broken hard check refuses the run.
    """
    hard_fail = [c for c in checks if c.hard and not c.ok]
    soft_fail = [c for c in checks if not c.hard and not c.ok]
    tight = [c for c in checks if c.ok and c.margin < MARGIN_WARN]
    if hard_fail:
        v = "FAIL"
    elif soft_fail:
        v = f"PASS ({len(soft_fail)} warnings)"
    else:
        v = "PASS"
    return v, hard_fail, soft_fail, tight


# -- dt selection: two cases independently arrived at "1% of gamma/(local stiffness)" --
def relaxation_time(gamma, stiffness):
    """tau = gamma/k. For a trap k is the trap stiffness; for a soft pair k = U''(r_min).

    The same formula appeared in two cases under different names (tau_k, tau_int) —
    it is one structure.
    """
    return (gamma / stiffness).to("s")


def dt_from_gate(tau, gate: float = GATE):
    """하드 게이트에 맞춘 dt = gate·τ. gate=1e-2 는 선형계 기준 편향 0.5%."""
    return gate * tau


def dt_from_bias(tau, bias: float):
    """목표 편향에서 역산 (bd-physics §1.2): 편향 ≈ (dt/τ)/2 ⟹ dt = 2·bias·τ.

    ⚠️ 닫힌 형태는 **선형계에만** 성립합니다. 비선형계는 dt 절반 수렴 확인이 필요합니다.
    """
    return 2 * bias * tau


def bias_from_dt(dt, tau) -> float:
    """dt에서 예상 계통 편향 [%] (선형계)."""
    return 50.0 * float((dt / tau).to("dimensionless").magnitude)


__all__ = ["Check", "GATE", "MARGIN_WARN", "SOFT_KINDS", "soft", "verdict",
           "relaxation_time", "dt_from_gate", "dt_from_bias", "bias_from_dt"]
