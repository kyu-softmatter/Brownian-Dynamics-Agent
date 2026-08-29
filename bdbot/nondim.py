"""L3 `NondimSpec` — 무차원화의 산출물. 마스터플랜 §6, skill `bd-physics` §0 ③④⑤.

L2(`system.yaml`, SI) 와 L4(실행) 사이의 **유일한 계약**입니다. 실행은 이 스펙만 보고
할 수 있어야 하고, 결과는 이 스펙만 보고 물리 단위로 되돌릴 수 있어야 합니다.

## 왜 만들었나 — 실측한 결함 4건 (`scratch/verify_l3_spec_gaps.py`)

세 케이스가 각자 `spec = {...}` 딕셔너리를 손으로 만들고 있었고, 서로 스키마가 달랐습니다
(공통 키가 `N`·`n_eq`·`n_prod` 세 개뿐). 그래서:

  ① ⭐️ **`run_id` 가 물리계를 덮지 않았습니다.** 1-B 스펙에는 물리계가 아예 없어서
     `d` 5µm→0.5µm, `η` 물→글리세롤(62배), `ρ_p` 실리카→폴리스티렌으로 바꿔도
     `run_id` 가 `soft-r3-2d-A-sweep__A100__27f70deab9` 로 **완전히 같았습니다**.
     그건 이미 완료된 런의 이름이므로 `prepare_outdir` 가 "이미 완료된 런입니다" 하고
     건너뛰고, **예전 계의 결과를 새 계의 결과로 보고합니다.** τ_B 가 16.1배 다른 계인데도.
     → 스펙은 `system`(physics_only)을 **반드시** 포함합니다.
  ② 스펙만으로 역변환이 불가능했습니다 (σ·τ·kT의 SI 값이 없음) — bd-physics §5 위반.
     → `back_transform` 에 세 앵커를 SI 부동소수로 박아둡니다.
  ③ 무차원수가 "정말 그 두 스케일의 비인가"를 검사할 수 없었습니다.
     → `Group.num`/`den` 이 원장 기호를 가리키고, `validate()` 가 재계산해 대조합니다.
  ④ 원장에서 빠진 스케일을 아무도 잡지 못했습니다.
     → `ScaleLedger.missing_roles()` 를 하드 게이트로 씁니다.

## 하지 않은 것 (아직 한 번도 안 나왔으므로 — CLAUDE.md 추상화 규칙)

  · `thermal` 이외의 기준 전략 (`interaction`·`active`·`custom`) — bd-physics §1 참조
  · 역구성 `from_dimensionless(groups, anchors)` — 무차원에서 출발하는 케이스가 없습니다
  · 원장을 **대신 계산해주는** 엔진 — 어떤 스케일이 있는지는 계마다 다르고,
    그게 케이스 스크립트에 남겨둔 물리 판단입니다 (bd-physics §6.3)

`run_id` 해시에 들어가는 것은 `{system, params, numerics}` 뿐입니다 — `schema`·원장·검사·
근거·유도값은 제외합니다. 스키마 버전을 올리거나 주석을 고쳐서 런이 무효화되면
콘텐츠 주소는 쓸모가 없습니다 (`bdbot.runid` 의 교훈).
"""
from __future__ import annotations

import json
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from . import intake as _intake
from . import runid as _runid
from . import scales as _scales
from .checks import Check
from .checks import verdict as _verdict
from .units import Q

SCHEMA = "bdbot.nondim/0.1"
RATIO_RTOL = 1e-9          # 무차원수 = 원장 두 항목의 비. 부동소수 오차만 허용.


# ══════════════════════════════════════════════════════════════════════
# ② 기준 스케일
# ══════════════════════════════════════════════════════════════════════
@dataclass
class Reference:
    """기준으로 택한 길이·시간·에너지 + 선택 근거 (마스터플랜 §6.2).

    ⚠️ 기준이 다르면 같은 계가 전혀 다른 무차원수를 갖습니다. 문헌과 대조할 때
       상대가 무엇을 썼는지 확인해야 하므로 `strategy` 와 `rationale` 을 함께 저장합니다.
    """

    length: tuple           # (기호, Quantity)
    time: tuple
    energy: tuple
    strategy: str = "thermal"
    rationale: str = ""

    @classmethod
    def from_dict(cls, d: dict) -> "Reference":
        """`scales.thermal_reference()` 의 반환값을 그대로 받습니다 (케이스 코드 무변경)."""
        return cls(length=tuple(d["length"]), time=tuple(d["time"]),
                   energy=tuple(d["energy"]), strategy=d.get("strategy", "thermal"),
                   rationale=d.get("rationale", ""))

    def as_dict(self) -> dict:
        """리포트 렌더러(`bdbot.report.render`)가 기대하는 dict 모양으로 되돌립니다."""
        return {"length": self.length, "time": self.time, "energy": self.energy,
                "strategy": self.strategy, "rationale": self.rationale}

    def si(self, kind: str):
        return {"length": self.length, "time": self.time, "energy": self.energy}[kind][1]

    def to_json(self) -> dict:
        out = {"strategy": self.strategy, "rationale": self.rationale}
        for kind in ("length", "time", "energy"):
            sym, q = getattr(self, kind)
            qb = q.to_base_units()
            out[kind] = {"symbol": sym, "value": float(qb.magnitude), "unit": str(qb.units)}
        return out

    @classmethod
    def from_json(cls, d: dict) -> "Reference":
        g = lambda k: (d[k]["symbol"], Q(d[k]["value"], d[k]["unit"]))
        return cls(length=g("length"), time=g("time"), energy=g("energy"),
                   strategy=d.get("strategy", "thermal"), rationale=d.get("rationale", ""))


# ══════════════════════════════════════════════════════════════════════
# ③ 무차원수 = 두 스케일의 비
# ══════════════════════════════════════════════════════════════════════
@dataclass
class Group:
    """무차원수 하나 (마스터플랜 §6.3).

    `num`/`den` 을 주면 **원장에서 재계산해 대조합니다.** 이게 이 클래스의 존재 이유입니다 —
    숫자만 담은 dict 는 "이름은 dt/τ_int 인데 값은 다른 것"을 잡을 수 없었습니다.

        Group("dt/tau_int", 0.01, num=("times", "dt"), den=("times", "tau_int"),
              meaning="적분 해상")

    비가 아닌 것(φ, A 처럼 정의상 무차원인 입력)은 `num`/`den` 없이 둡니다.
    """

    name: str
    value: float
    num: tuple | None = None        # (category, symbol)
    den: tuple | None = None
    expr: str = ""
    meaning: str = ""

    @property
    def label(self) -> str:
        """리포트 한 줄 이름 — `이름  = 식   의미` 모양."""
        head = f"{self.name:<12}" + (f"= {self.expr:<14}" if self.expr else " " * 16)
        return (head + self.meaning).rstrip()

    def recompute(self, ledger) -> float | None:
        """원장에서 비를 다시 계산. num/den 이 없으면 None."""
        if not (self.num and self.den):
            return None
        nc, ns = self.num
        dc, ds = self.den
        return float((ledger.get(nc, ns) / ledger.get(dc, ds)).to("dimensionless").magnitude)

    def mismatch(self, ledger) -> str | None:
        """불일치 메시지 또는 None. KeyError(원장에 없는 기호)도 결함으로 봅니다."""
        try:
            want = self.recompute(ledger)
        except KeyError as e:
            return f"원장에 없는 기호를 가리킵니다: {e}"
        except Exception as e:                      # 단위가 안 맞는 비 → 무차원이 아님
            return f"비를 계산할 수 없습니다 ({e})"
        if want is None:
            return None
        if abs(self.value - want) > RATIO_RTOL * max(abs(want), 1e-300):
            return (f"적힌 값 {self.value:.12g} ≠ 원장 재계산 "
                    f"{self.num[1]}/{self.den[1]} = {want:.12g}")
        return None

    def to_json(self) -> dict:
        return {"name": self.name, "value": self.value, "expr": self.expr,
                "meaning": self.meaning,
                "num": list(self.num) if self.num else None,
                "den": list(self.den) if self.den else None}

    @classmethod
    def from_json(cls, d: dict) -> "Group":
        return cls(name=d["name"], value=d["value"], expr=d.get("expr", ""),
                   meaning=d.get("meaning", ""),
                   num=tuple(d["num"]) if d.get("num") else None,
                   den=tuple(d["den"]) if d.get("den") else None)


def groups_dict(groups) -> dict:
    """`bdbot.report.render` 이 받는 {표시이름: 값} 으로 — **사람이 읽는 리포트용**."""
    return {g.label: g.value for g in groups}


def metrics_dict(groups) -> dict:
    """`metrics.json` 의 `dimensionless` 용 — **기계가 읽는 키**.

    예전에는 리포트 표시이름을 그대로 키로 썼습니다
    (`'k*     = k d²/kT   트랩 vs 열요동'`). metrics.json 은 postmortem 의 유일한
    입력인데 이런 키로는 질의를 할 수 없습니다. 기호만 씁니다 (`'k*'`, `'dt/tau_k'`).
    """
    return {g.name: g.value for g in groups}


# ══════════════════════════════════════════════════════════════════════
# L3 산출물
# ══════════════════════════════════════════════════════════════════════
@dataclass
class NondimSpec:
    """무차원화 결과 전체. `specs/<run_id>.json` 으로 저장되고 L4가 이것만 읽습니다."""

    case: str
    system: dict                                  # L2 원본 (physics_only 적용 전)
    reference: Reference
    ledger: Any                                   # ScaleLedger
    groups: list = field(default_factory=list)
    checks: list = field(default_factory=list)
    numerics: dict = field(default_factory=dict)  # 무차원 실행 파라미터 (dt_star, n_prod, seed…)
    params: dict = field(default_factory=dict)    # 케이스별 무차원 손잡이 (A, phi, r_c_star…)
    tag: str | None = None
    nhex: int = 12                                # 1-A는 12자, 1-B는 10자 (재현성 유지)
    label: str | None = None

    def __post_init__(self):
        if self.label is None:
            self.label = self.system.get("label", self.case)
        if isinstance(self.reference, dict):       # thermal_reference() 를 바로 받기
            self.reference = Reference.from_dict(self.reference)

    # ── 콘텐츠 주소 ────────────────────────────────────────────────────
    def hash_payload(self) -> dict:
        """`run_id` 해시의 대상 — **물리를 정하는 것만**.

        `schema`·원장·검사·근거·유도값은 넣지 않습니다. 스키마 버전을 올리거나 주석을
        고쳐서 런이 무효화되면 콘텐츠 주소가 쓸모없어집니다 (`bdbot.runid` DOC_KEYS 교훈).
        """
        return {"system": _runid.physics_only(self.system),
                "params": _runid.physics_only(self.params),
                "numerics": _runid.physics_only(self.numerics)}

    def run_id(self) -> str:
        return _runid.content_run_id(self.label, self.hash_payload(),
                                     tag=self.tag, nhex=self.nhex)

    # ── ④ 판정 ────────────────────────────────────────────────────────
    def verdict(self) -> str:
        return _verdict(self.checks)[0]

    def validate(self) -> list:
        """L3 자체의 무결성 검사 — 물리 검사(`checks`)와 **다른 층**입니다.

        여기서 보는 것: 원장이 완전한가 · 무차원수가 정말 그 비인가 · 기준이 원장에 있는가 ·
        역변환 앵커가 성립하는가. 즉 "무차원화를 제대로 했는가"이고,
        `checks` 는 "이 계에 BD가 타당하고 충분히 잘게 적분하는가" 입니다.
        """
        I = _intake.Issue
        out = []

        # ① 원장 완전성 — 빠진 스케일은 돌지 않는 검사다
        for role in self.ledger.missing_roles():
            out.append(I("error", f"ledger.{role}",
                         f"필수 역할 '{role}' ({_scales.ROLE_MEANING.get(role, '')}) 이 "
                         f"원장에 없습니다. 올리거나 "
                         f"`declare_absent('{role}', 이유)` 로 명시하세요."))

        # ② 무차원수 = 원장 두 항목의 비인가
        n_checked = 0
        for g in self.groups:
            msg = g.mismatch(self.ledger)
            if msg:
                out.append(I("error", f"groups.{g.name}", msg))
            elif g.num and g.den:
                n_checked += 1
        if self.groups and n_checked == 0:
            out.append(I("warn", "groups",
                         "원장과 대조된 무차원수가 하나도 없습니다 — 전부 num/den 이 "
                         "비어 있습니다. 비인 것에는 원장 기호를 붙이세요 (§6.3)."))

        # ③ 기준 스케일이 원장에 실제로 있는가 (기준만 있고 원장에 없으면 대조가 불가능)
        for kind, cat in (("length", "lengths"), ("time", "times"), ("energy", "energies")):
            sym = getattr(self.reference, kind)[0]
            if not self.ledger.has(cat, sym):
                out.append(I("error", f"reference.{kind}",
                             f"기준 '{sym}' 이 원장 {cat} 에 없습니다."))

        # ④ 역변환 앵커가 양수·유한인가 (0이면 되돌릴 수 없다)
        for kind in ("length", "time", "energy"):
            q = self.reference.si(kind)
            try:
                v = float(q.to_base_units().magnitude)
            except Exception as e:
                out.append(I("error", f"reference.{kind}", f"SI 변환 실패: {e}"))
                continue
            if not (v > 0) or v != v or v in (float("inf"), float("-inf")):
                out.append(I("error", f"reference.{kind}",
                             f"기준 스케일이 양수·유한이 아닙니다 ({v!r}) — 역변환 불가."))

        # ⑤ 필수 실행 파라미터 (L4가 이것만 보고 돌린다)
        for k in ("dt_star", "n_prod"):
            if k not in self.numerics:
                out.append(I("error", f"numerics.{k}",
                             "L4가 스펙만으로 실행할 수 없습니다 (무차원 실행 파라미터 누락)."))
        dt_star = self.numerics.get("dt_star")
        if isinstance(dt_star, (int, float)) and not dt_star > 0:
            out.append(I("error", "numerics.dt_star", f"dt* 가 양수가 아닙니다 ({dt_star!r})."))

        # ⑥ `dt_star` 는 원장의 dt/기준시간과 **같아야** 합니다. 둘을 따로 계산하고 있으므로
        #    어긋날 수 있고, 어긋나면 HOOMD가 원장과 다른 스텝으로 돕니다 — 분리 검사는
        #    원장의 dt로 통과하는데 실제 적분은 다른 값이라, 조용히 틀리는 유형입니다.
        if isinstance(dt_star, (int, float)) and self.ledger.has("times", "dt"):
            try:
                want = float((self.ledger.get("times", "dt")
                              / self.reference.si("time")).to("dimensionless").magnitude)
            except Exception as e:
                out.append(I("error", "numerics.dt_star", f"원장 dt와 대조 실패: {e}"))
            else:
                if abs(dt_star - want) > RATIO_RTOL * max(abs(want), 1e-300):
                    out.append(I("error", "numerics.dt_star",
                                 f"원장과 어긋납니다: numerics.dt_star = {dt_star:.12g} ≠ "
                                 f"dt/{self.reference.time[0]} = {want:.12g}. "
                                 "HOOMD가 원장과 다른 스텝으로 돌게 됩니다."))
        return out

    @property
    def errors(self) -> list:
        return [i for i in self.validate() if i.level == "error"]

    # ── ⑤ 역변환 (bd-physics §5 — 항상 수행) ──────────────────────────
    def physical(self, value, L: int = 0, T: int = 0, E: int = 0):
        """무차원 값 → 물리 단위. 차원 지수로 지정합니다.

            spec.physical(D_star, L=2, T=-1)     # D_eff → m²/s
            spec.physical(x2_star, L=2)          # ⟨x²⟩  → m²
            spec.physical(step * dt_star, T=1)   # 시간  → s
            spec.physical(P_star, E=1, L=-2)     # 2D 압력 → N/m
        """
        sigma = self.reference.si("length")
        tau = self.reference.si("time")
        en = self.reference.si("energy")
        return value * sigma**L * tau**T * en**E

    def back_transform(self) -> dict:
        """역변환 3앵커를 SI 부동소수로. 스펙만 있으면 pint 없이도 되돌릴 수 있습니다."""
        out = {}
        for kind, key in (("length", "sigma_SI"), ("time", "tau_SI"), ("energy", "energy_SI")):
            qb = self.reference.si(kind).to_base_units()
            out[key] = float(qb.magnitude)
            out[key + "_unit"] = str(qb.units)
        return out

    # ── 직렬화 ────────────────────────────────────────────────────────
    def _ledger_json(self) -> dict:
        """원장을 SI + 무차원(기준으로 나눈 값) 둘 다 저장. L4는 star 쪽을 씁니다."""
        out = {}
        ref_of = {"lengths": "length", "times": "time", "energies": "energy"}
        for cat_name, cat in self.ledger.categories():
            rows = []
            base = self.reference.si(ref_of[cat_name])
            for sym, sc in cat.items():
                q = sc.value if isinstance(sc, _scales.Scale) else sc
                qb = q.to_base_units()
                row = {"symbol": sym, "value": float(qb.magnitude), "unit": str(qb.units)}
                if isinstance(sc, _scales.Scale):
                    row.update(note=sc.note, role=sc.role, star=sc.star)
                try:
                    row["reduced"] = float((q / base).to("dimensionless").magnitude)
                except Exception:
                    row["reduced"] = None       # 같은 종류가 아닌 항목 (있으면 안 되지만)
                rows.append(row)
            out[cat_name] = rows
        return out

    def to_json(self) -> dict:
        return {
            "schema": SCHEMA,
            "case": self.case,
            "label": self.label,
            "tag": self.tag,
            "run_id": self.run_id(),
            "system": self.system,
            "reference": self.reference.to_json(),
            "back_transform": self.back_transform(),
            "ledger": self._ledger_json(),
            "ledger_absent": dict(self.ledger.absent),
            "groups": [g.to_json() for g in self.groups],
            "checks": [c.as_dict() for c in self.checks],
            "verdict": self.verdict(),
            "params": self.params,
            "numerics": self.numerics,
            "l3_issues": [{"level": i.level, "where": i.where, "msg": i.msg}
                          for i in self.validate()],
        }

    def write(self, path) -> Path:
        p = Path(path)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(json.dumps(self.to_json(), indent=2, ensure_ascii=False, default=str))
        return p


# ══════════════════════════════════════════════════════════════════════
# 읽기 — L4는 이 함수만 씁니다 (케이스 스크립트를 임포트하지 않습니다)
# ══════════════════════════════════════════════════════════════════════
@dataclass
class LoadedSpec:
    """`specs/*.json` 을 되읽은 것. 원장이 Quantity 로 복원되어 역변환이 됩니다.

    `NondimSpec` 을 그대로 복원하지 않는 이유: 저장된 스펙에는 `ScaleLedger.derived`
    (케이스 중간값)가 없고 있을 필요도 없습니다. L4가 필요한 건 무차원 파라미터와
    역변환 앵커뿐입니다.
    """

    raw: dict
    reference: Reference
    groups: list
    checks: list

    @property
    def case(self) -> str:
        return self.raw["case"]

    @property
    def label(self) -> str:
        return self.raw.get("label", self.raw["case"])

    @property
    def run_id(self) -> str:
        return self.raw["run_id"]

    @property
    def params(self) -> dict:
        return self.raw.get("params", {})

    @property
    def numerics(self) -> dict:
        return self.raw.get("numerics", {})

    @property
    def verdict(self) -> str:
        return self.raw.get("verdict", "?")

    def reduced(self, cat: str, symbol: str) -> float:
        """무차원 값 (기준으로 나눈 것). L4가 L*·r_c*·dt* 를 여기서 가져옵니다."""
        for row in self.raw["ledger"][cat]:
            if row["symbol"] == symbol:
                return row["reduced"]
        raise KeyError(f"{cat}.{symbol} 이 스펙 원장에 없습니다")

    def si(self, cat: str, symbol: str):
        for row in self.raw["ledger"][cat]:
            if row["symbol"] == symbol:
                return Q(row["value"], row["unit"])
        raise KeyError(f"{cat}.{symbol} 이 스펙 원장에 없습니다")

    def group(self, name: str) -> float:
        for g in self.groups:
            if g.name == name:
                return g.value
        raise KeyError(f"무차원수 '{name}' 이 스펙에 없습니다")

    def physical(self, value, L: int = 0, T: int = 0, E: int = 0):
        bt = self.raw["back_transform"]
        sigma = Q(bt["sigma_SI"], bt["sigma_SI_unit"])
        tau = Q(bt["tau_SI"], bt["tau_SI_unit"])
        en = Q(bt["energy_SI"], bt["energy_SI_unit"])
        return value * sigma**L * tau**T * en**E

    def render(self) -> str:
        """스펙**만** 보고 리포트를 다시 그립니다 — 케이스 스크립트 없이.

        이게 되면 스펙이 자족적이라는 뜻이고, L4가 케이스 코드를 임포트하지 않아도
        됩니다. 안 되면 스펙에 무언가 빠진 것입니다.
        """
        r = self.raw
        L: list[str] = []
        w = L.append
        W = 88
        w("=" * W)
        w(f"NondimSpec — {self.label}   run_id={self.run_id}")
        w("=" * W)
        ref = r["reference"]
        w(f"기준 스케일: length={ref['length']['symbol']}  energy={ref['energy']['symbol']}"
          f"  time={ref['time']['symbol']}   [strategy: {ref['strategy']}]")
        w(f"  근거: {ref.get('rationale', '')}")

        ok, want = self.verify_hash()
        w("")
        w(f"해시 자기검증: {'✓ 스펙 내용과 run_id 가 일치' if ok else f'✗ 불일치 — 기대 {want}'}")

        w("")
        w("SCALE LEDGER  (SI · 기준 대비 환산)")
        for cat in ("lengths", "times", "energies"):
            rows = sorted(r["ledger"].get(cat, []), key=lambda x: x["value"])
            if not rows:
                continue
            w(f"  {cat}")
            for row in rows:
                star = " ★" if row.get("star") else ""
                role = f"  [{row['role']}]" if row.get("role") else ""
                red = row.get("reduced")
                red_s = f"{red:.6g}" if isinstance(red, (int, float)) else "—"
                w(f"    {row['symbol']:<11}{row['value']:>13.5e} {row['unit']:<22}"
                  f"{red_s:>13}{role}{star}")
        if r.get("ledger_absent"):
            for role, why in r["ledger_absent"].items():
                w(f"    (없음) {role}: {why}")

        w("")
        w("DIMENSIONLESS GROUPS")
        for g in self.groups:
            ratio = f"{g.num[1]}/{g.den[1]}" if (g.num and g.den) else "— (비 아님)"
            w(f"  {g.name:<16}{g.value:>13.6g}   {ratio:<22}{g.meaning}")

        w("")
        w(f"{'SEPARATION CHECKS':<58}{'value':>10}{'limit':>10}{'margin':>9}")
        for c in self.checks:
            mark = "✓" if c.ok else ("✗" if c.hard else "⚠")
            w(f"  {mark} [{c.kind}] {c.name:<42}{c.value:>10.3e}{c.limit:>10.0e}"
              f"{c.margin:>8.1f}×")

        w("")
        w("BACK TRANSFORM  (결과 → 물리 단위)")
        bt = r["back_transform"]
        w(f"  σ = {bt['sigma_SI']:.6e} {bt['sigma_SI_unit']}    "
          f"τ = {bt['tau_SI']:.6e} {bt['tau_SI_unit']}")
        w(f"  E = {bt['energy_SI']:.6e} {bt['energy_SI_unit']}")

        w("")
        w("RUN PARAMETERS  (L4가 읽는 것)")
        for k, v in r.get("params", {}).items():
            w(f"  params.{k:<20} {v}")
        for k, v in r.get("numerics", {}).items():
            w(f"  numerics.{k:<18} {v}")

        issues = r.get("l3_issues", [])
        if issues:
            w("")
            w("L3 무결성")
            for i in issues:
                mark = {"error": "✗", "warn": "⚠", "info": "ℹ"}.get(i["level"], "·")
                w(f"  {mark} [{i['where']}] {i['msg']}")

        w("")
        w(f"VERDICT: {self.verdict}")
        w("=" * W)
        return "\n".join(L)

    def verify_hash(self) -> tuple[bool, str]:
        """저장된 `run_id` 가 스펙 내용과 일치하는가 — 손으로 고친 스펙을 잡습니다.

        마스터플랜 §16 규칙 2("스펙을 손으로 쓰지 않는다")를 기계가 확인하는 자리입니다.
        """
        payload = {"system": _runid.physics_only(self.raw.get("system", {})),
                   "params": _runid.physics_only(self.params),
                   "numerics": _runid.physics_only(self.numerics)}
        stored = self.raw["run_id"]
        nhex = len(stored.rsplit("__", 1)[-1])
        want = _runid.content_run_id(self.label, payload,
                                     tag=self.raw.get("tag"), nhex=nhex)
        return want == stored, want


def load(path) -> LoadedSpec:
    raw = json.loads(Path(path).read_text())
    got = raw.get("schema")
    if got != SCHEMA:
        raise ValueError(f"스키마가 다릅니다: {got!r} (기대 {SCHEMA!r})")
    return LoadedSpec(
        raw=raw,
        reference=Reference.from_json(raw["reference"]),
        groups=[Group.from_json(g) for g in raw.get("groups", [])],
        checks=[Check(kind=c["kind"], name=c["name"], value=c["value"], limit=c["limit"],
                      op=c.get("op", "<="), note=c.get("note", ""), hard=c.get("hard", True))
                for c in raw.get("checks", [])],
    )


__all__ = ["SCHEMA", "Reference", "Group", "NondimSpec", "LoadedSpec", "load",
           "groups_dict", "RATIO_RTOL"]
