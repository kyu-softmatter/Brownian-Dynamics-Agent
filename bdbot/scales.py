"""스케일 원장 — 마스터플랜 원칙 3 ①②, skill `bd-physics` §0.

계에 있는 모든 길이·시간·에너지를 SI로 열거하고, 그중 기준을 **근거와 함께** 고릅니다.
무차원수를 먼저 정하고 스케일을 나중에 유추하는 순서는 금지입니다.

1-A·1-B가 똑같은 구조(lengths/times/energies + ref + rationale)를 갖고 있었습니다.
**케이스마다 다른 것은 어떤 스케일이 들어가는지**이고, 그건 케이스가 채웁니다.

## L3 작업에서 바꾼 것 (2026-08-04)

원장 항목의 키가 `"d        입자 지름 (기준)"` 처럼 **기호와 설명이 붙은 한 덩어리**였습니다.
사람이 읽는 표를 만들기엔 충분했지만 기계가 쓸 수 없었습니다:

  · `ratio("times", a, b)` 를 쓰려면 장식까지 정확히 맞춰야 했습니다.
  · 무차원수가 "정말 그 두 스케일의 비인가"를 **검사할 수 없었습니다** —
    `dt/τ_int` 이라 적어놓고 다른 값을 넣어도 아무도 잡지 못합니다.
    실제로 이 프로젝트는 케이지 강성을 `a_mean` 에서 평가하는 (a_NN 이어야 함) 실수를
    했고 41% 어긋났습니다. 그 종류의 실수를 잡는 자리가 없었습니다.
  · 원장에서 **빠뜨린 스케일**을 아무도 잡지 못했습니다. 빠진 스케일은 곧
    **돌지 않는 검사**입니다 — "조용히 통과"와 "검사를 안 함"은 다릅니다.

그래서 `Scale(symbol, value, note, role)` 로 쪼갰습니다. `role` 은 기호 이름과 무관하게
"이 계에서 박스/적분스텝/관측창/관성 역할을 하는 스케일"을 가리킵니다 — 타원체 케이스가
기준 길이를 `d_eq` 로 부르므로 기호로 필수 항목을 강제할 수 없습니다.

없는 스케일은 `declare_absent(role, 이유)` 로 **명시적으로** 비웁니다 (규칙 3과 같은 태도 —
모르는 것을 지어내지 않고, 없다는 것을 적는다).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

THERMAL_RATIONALE = (
    "thermal 규약 고정 (σ=d, E=kT, τ=τ_B). 단위계는 문헌 비교를 위해 계마다 바꾸지 "
    "않는다 (bd-physics §1.1). 지배 시간척도는 계마다 다를 수 있고, dt는 그쪽으로 정한다."
)

# 원장에 반드시 있어야 하는 역할. 없으면 `declare_absent` 로 이유를 적어야 합니다.
# 근거: 마스터플랜 §6.1 "반드시 계산해 원장에 올릴 스케일 (해당 없으면 None으로 명시적 기록)".
#   box         기하 검사(최소 이미지·유한크기)의 분모
#   dt          적분 해상 검사의 분자 — 원장에 없으면 시간척도 정렬에서 안 보인다
#   observation 통계 검사(T_obs/τ)의 분자
#   inertia     모델 타당성 검사(τ_p/τ_dyn) — BD를 써도 되는 계인가
MANDATORY_ROLES = ("box", "dt", "observation", "inertia")
ROLE_MEANING = {
    "box": "박스 한 변 L — 기하 검사의 분모",
    "dt": "적분 시간 스텝",
    "observation": "관측창 T_obs",
    "inertia": "관성 이완 τ_p — 모델 타당성 판정용",
    "ref_length": "기준 길이", "ref_time": "기준 시간", "ref_energy": "기준 에너지",
}


@dataclass
class Scale:
    """원장의 한 항목. 기호와 설명을 **분리**해서 기계가 기호로 접근할 수 있게 합니다."""

    symbol: str
    value: Any                   # pint Quantity (SI)
    note: str = ""
    role: str = ""               # MANDATORY_ROLES 중 하나이거나 빈 문자열
    star: bool = False           # 리포트에서 ★ 강조 (이 계의 지배 스케일)

    @property
    def display(self) -> str:
        """리포트 한 줄의 이름 — 예전 키 모양(`d        입자 지름 (기준)`)을 그대로 만듭니다.

        ⚠️ 패딩은 **최소 폭**입니다. `chain-bend` 가 `kappa_end_d2`(12자) 같은 기호를 쓰면서
           고정 9자로는 기호와 설명이 붙어버렸습니다 (`delta_maxM_c 선형 탄성 한계`).
           표시 전용이라 원장 값·run_id 에는 영향이 없습니다.
        """
        s = f"{self.symbol:<12} {self.note}"
        return (s + " ★") if self.star else s


def _qty(v):
    """Scale 이든 Quantity 든 값을 꺼냅니다 (원장은 둘 다 받습니다)."""
    return v.value if isinstance(v, Scale) else v


@dataclass
class ScaleLedger:
    """길이·시간·에너지 원장 + 기준 스케일 + 선택 근거."""

    lengths: dict = field(default_factory=dict)     # symbol → Scale
    times: dict = field(default_factory=dict)
    energies: dict = field(default_factory=dict)
    derived: dict = field(default_factory=dict)     # 케이스가 쓰는 중간값 (원장 아님)
    ref: dict = field(default_factory=dict)
    rationale: str = ""
    absent: dict = field(default_factory=dict)      # role → 없는 이유

    # ── 쓰기 ──────────────────────────────────────────────────────────
    def add(self, cat: str, symbol: str, value, note: str = "", role: str = "",
            star: bool = False) -> "ScaleLedger":
        """원장에 한 항목. `cat` 은 "lengths"/"times"/"energies"."""
        getattr(self, cat)[symbol] = Scale(symbol, value, note, role, star)
        return self

    def add_length(self, symbol, value, note="", role="", star=False):
        return self.add("lengths", symbol, value, note, role, star)

    def add_time(self, symbol, value, note="", role="", star=False):
        return self.add("times", symbol, value, note, role, star)

    def add_energy(self, symbol, value, note="", role="", star=False):
        return self.add("energies", symbol, value, note, role, star)

    def declare_absent(self, role: str, reason: str) -> "ScaleLedger":
        """이 계에 해당 역할의 스케일이 **없다**고 명시. 이유 없이는 비울 수 없습니다."""
        if not reason:
            raise ValueError(f"role '{role}' 을 비우려면 이유가 필요합니다 (규칙 3)")
        self.absent[role] = reason
        return self

    # ── 읽기 ──────────────────────────────────────────────────────────
    def categories(self):
        return (("lengths", self.lengths), ("times", self.times), ("energies", self.energies))

    def sorted_items(self, cat: dict):
        """작은 것부터 (표시이름, 값). 정렬해 놓으면 분리 위반이 눈에 보입니다."""
        items = sorted(cat.items(), key=lambda kv: _qty(kv[1]).to_base_units().magnitude)
        return [((v.display if isinstance(v, Scale) else k), _qty(v)) for k, v in items]

    def get(self, cat: str, symbol: str):
        """기호로 값 하나. 없으면 KeyError — 조용히 None을 주지 않습니다."""
        return _qty(getattr(self, cat)[symbol])

    def has(self, cat: str, symbol: str) -> bool:
        return symbol in getattr(self, cat)

    def by_role(self, role: str):
        """역할로 찾기. (기호, 값) 또는 None."""
        for _, cat in self.categories():
            for sym, sc in cat.items():
                if isinstance(sc, Scale) and sc.role == role:
                    return sym, sc.value
        return None

    def ratio(self, cat: str, a: str, b: str) -> float:
        """두 스케일의 비 — 무차원수는 전부 이 형태여야 합니다 (bd-physics §3)."""
        dd = getattr(self, cat)
        return float((_qty(dd[a]) / _qty(dd[b])).to("dimensionless").magnitude)

    def span(self, cat: str) -> float:
        """가장 큰 스케일 / 가장 작은 스케일. 계가 몇 자릿수에 걸쳐 있는가."""
        dd = getattr(self, cat)
        if len(dd) < 2:
            return 1.0
        vs = [_qty(v).to_base_units().magnitude for v in dd.values()]
        return max(vs) / min(vs)

    # ── 완전성 ────────────────────────────────────────────────────────
    def missing_roles(self) -> list[str]:
        """필수 역할 중 원장에도 없고 `declare_absent` 도 안 된 것.

        **빠진 스케일 = 돌지 않는 검사**입니다. 1-B가 `dt`·`T_obs` 를 원장에 올리지 않고
        따로 넘기고 있었습니다 — 시간척도 정렬 표에 둘이 안 보였습니다.
        """
        return [r for r in MANDATORY_ROLES
                if self.by_role(r) is None and r not in self.absent]


def thermal_reference(d, kT, tau_B, rationale: str | None = None, *,
                      length_symbol: str = "d", energy_symbol: str = "kT",
                      time_symbol: str = "tau_B") -> dict:
    """`thermal` 기준 셋. HOOMD 스펙에서 σ=1, kT=1, γ=1, τ_B=1 이 됩니다.

    기호 이름을 바꿀 수 있게 둔 이유: 타원체 케이스의 기준 길이는 등가부피 구 지름이고
    원장에서 `d_eq` 로 부릅니다. 기호를 "d" 로 고정해두면 리포트가 원장에 없는 기준을
    가리키게 되고, `NondimSpec.validate()` 가 그걸 오류로 잡습니다 (실제로 잡았습니다).
    """
    return {"length": (length_symbol, d), "energy": (energy_symbol, kT),
            "time": (time_symbol, tau_B),
            "strategy": "thermal", "rationale": rationale or THERMAL_RATIONALE}


__all__ = ["ScaleLedger", "Scale", "thermal_reference", "THERMAL_RATIONALE",
           "MANDATORY_ROLES", "ROLE_MEANING"]
