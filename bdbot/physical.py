"""L2 `PhysicalSystem` — 차원 있는 물리계 (SI). 마스터플랜 원칙 3의 출발점.

1-A·1-B가 각자 `load_system()`을 갖고 있었습니다 — **두 번 나왔으므로** 공통화합니다.
스키마는 두 `system.yaml`의 실사용에서 갈랐습니다:

    2/2 (공통)  label · description · dimensions · particle · medium ·
                interactions · external · targets · numerics
    1/2 (선택)  geometry · derived_scales · dimensionless ·
                required_convergence_checks · not_verified

**Provenanced 잎** (`value`+`unit`+`source`+`tier`)은 2/2에서 동일했습니다:
`particle.{diameter,density,count}` · `medium.{temperature,viscosity}` + 케이스별 상호작용.

`derived_scales`는 **일부러 예외**입니다 — γ·D_t·τ_B 같은 유도값은 출처가 아니라
**재계산으로** 검증됩니다 (`verify()`가 그걸 합니다).

강제하는 불변식 (마스터플랜 §5.4의 L2 버전):
  ① `derived_from` — 이 물리계가 나온 `observation.yaml`. 없으면 거부.
     (기존 두 파일은 이걸 **주석**에만 적어뒀습니다 — 기계가 못 읽습니다)
  ② L0이 BLOCKED면 L2는 존재할 수 없다 — 미해소 물리 결측이 있는데 물리계를
     확정했다면 어딘가에서 값을 지어낸 것입니다 (규칙 3).
  ③ tier ≥ 2 (미검증)만으로 구성된 값이 있으면 사람 승인 필요 (§12.4).
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path

import yaml

from . import intake as _intake
from . import materials as _mat
from .provenance import Provenanced, load_node
from .units import Q

SCHEMA = "bdbot.system/0.1"

REQUIRED_TOP = ("label", "dimensions", "particle", "medium", "targets", "numerics")
OPTIONAL_TOP = ("description", "geometry", "interactions", "external", "derived_scales",
                "dimensionless", "required_convergence_checks", "not_verified",
                "derived_from")
# 출처 없이 값만 있어도 되는 섹션 — 유도값이라 재계산으로 검증한다
DERIVED_SECTIONS = ("derived_scales", "dimensionless", "friction")
# 2/2 에서 같았던 Provenanced 경로
# ★ 세 번째 케이스(abp-rod, 타원체)에서 필요해진 일반화: 구 공식이 안 맞는 형상이 있다.
#   γ̄(Perrin, 2D 조화평균) = 7.21e-9 vs 3πηd_eq = 6.37e-9 — 13% 차이.
#   구 재계산 검증을 그대로 돌리면 **정상인 스펙을 오류로 잡는다.**
SPHERICAL_SHAPES = ("sphere", None)

CORE_PROVENANCED = {
    "d": ("particle", "diameter"),
    "rho_p": ("particle", "density"),
    "N": ("particle", "count"),
    "T": ("medium", "temperature"),
    "eta": ("medium", "viscosity"),
}
TIER_MEANING = {0: "직접입력/핸드북", 1: "문헌+검증 또는 확인된 관례",
                2: "문헌 미검증", 3: "임의 가정"}


@dataclass
class PhysicalSystem:
    path: Path
    raw: dict
    core: dict = field(default_factory=dict)      # 이름 → Provenanced
    issues: list = field(default_factory=list)

    @property
    def label(self) -> str:
        return self.raw.get("label", "?")

    @property
    def dim(self) -> int:
        return int(self.raw.get("dimensions", 0))

    @property
    def shape(self):
        """`particle.shape` — 없으면 구로 본다 (1-A·1-B가 그랬다)."""
        return (self.raw.get("particle") or {}).get("shape")

    @property
    def is_spherical(self) -> bool:
        return self.shape in SPHERICAL_SHAPES

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.level == "error"]

    def node(self, *path, required: bool = True):
        """케이스별 Provenanced 노드 접근. `sys_.node("external","stiffness")`."""
        cur = self.raw
        for k in path:
            if not isinstance(cur, dict) or k not in cur:
                if required:
                    raise KeyError(f"{'.'.join(map(str, path))} 없음 ({self.path})")
                return None
            cur = cur[k]
        return load_node(cur)

    def tiers(self) -> dict:
        out: dict[int, list[str]] = {}
        for name, p in self.core.items():
            out.setdefault(p.tier, []).append(name)
        for name, p in self._extra_provenanced().items():
            out.setdefault(p.tier, []).append(name)
        return dict(sorted(out.items()))

    def _extra_provenanced(self) -> dict:
        """core 밖의 Provenanced 잎 (케이스별 상호작용·기하 등)."""
        out: dict[str, Provenanced] = {}
        for path, node in _walk_provenanced(self.raw):
            top = path.split(".")[0].split("[")[0]
            if top in DERIVED_SECTIONS:
                continue
            if path in {".".join(p) for p in CORE_PROVENANCED.values()}:
                continue
            try:
                out[path] = load_node(node)
            except Exception:
                pass
        return out

    def bulk(self) -> dict | None:
        """구 + 뉴턴 유체 기본 물성. 두 케이스가 똑같이 이 묶음을 씁니다.

        ★ 예외를 던지지 않습니다. 단위가 깨진 스펙에서 검사기가 크래시하면
          "무엇이 틀렸는지" 대신 트레이스백이 나옵니다 — 적대적 테스트에서 실제로
          그랬습니다 (`furlong^2` 를 넣었더니 DimensionalityError 로 죽음).
        """
        if not all(k in self.core for k in ("d", "T", "eta", "rho_p")):
            return None
        try:
            return _mat.sphere_bulk(self.core["d"].value, self.core["T"].value,
                                    self.core["eta"].value, self.core["rho_p"].value)
        except Exception:
            return None


def _walk_provenanced(d, pre=""):
    """`value` + (`unit`|`source`|`tier`) 를 가진 잎을 (경로, 노드)로 열거."""
    out = []
    if isinstance(d, dict):
        if "value" in d and any(k in d for k in ("unit", "source", "tier")):
            out.append((pre, d))
        else:
            for k, v in d.items():
                out += _walk_provenanced(v, f"{pre}.{k}" if pre else str(k))
    elif isinstance(d, list):
        for i, v in enumerate(d):
            out += _walk_provenanced(v, f"{pre}[{i}]")
    return out


def load(path) -> PhysicalSystem:
    p = Path(path)
    if p.is_dir():
        p = p / "system.yaml"
    if not p.exists():
        s = PhysicalSystem(p, {})
        s.issues.append(_intake.Issue("error", str(p), "system.yaml 이 없습니다."))
        return s
    raw = yaml.safe_load(p.read_text()) or {}
    s = PhysicalSystem(p, raw)
    for name, keys in CORE_PROVENANCED.items():
        cur = raw
        ok = True
        for k in keys:
            if not isinstance(cur, dict) or k not in cur:
                ok = False
                break
            cur = cur[k]
        if ok:
            try:
                s.core[name] = load_node(cur)
            except Exception as e:
                s.issues.append(_intake.Issue("error", ".".join(keys), f"파싱 실패: {e}"))
    s.issues += validate(s)
    return s


def validate(s: PhysicalSystem) -> list:
    I = _intake.Issue
    raw = s.raw
    out: list = []

    # ① 필수 섹션
    for k in REQUIRED_TOP:
        if k not in raw:
            out.append(I("error", k, "필수 섹션 누락"))
    for name, keys in CORE_PROVENANCED.items():
        if name not in s.core:
            out.append(I("error", ".".join(keys), f"필수 Provenanced 노드 누락 (2/2 공통)"))

    # ② ★ 불변식: derived_from (§5.4의 L2 버전)
    df = raw.get("derived_from")
    if not df:
        out.append(I("error", "derived_from",
                     "이 물리계가 어느 observation.yaml에서 나왔는지 없습니다. "
                     "주석이 아니라 필드로 적어야 기계가 검증할 수 있습니다 (§5.4)."))
    else:
        ref = (s.path.parent / Path(str(df)).name)
        if not ref.exists():
            out.append(I("error", "derived_from", f"참조 파일이 없습니다: {df}"))
        else:
            # ③ ★ L0이 BLOCKED면 L2는 존재할 수 없다
            obs = _intake.load(ref)
            if obs.errors:
                out.append(I("error", "derived_from",
                             f"근거 observation.yaml에 스키마 오류 {len(obs.errors)}건이 있습니다."))
            if obs.open_missing:
                names = ", ".join(m.get("symbol", "?") for m in obs.open_missing)
                out.append(I("error", "derived_from",
                             f"L0에 미해소 물리 결측이 있는데 물리계가 확정돼 있습니다: {names}. "
                             "어딘가에서 값을 지어냈을 수 있습니다 (규칙 3)."))

    # ④ Provenanced 잎의 완전성 + 단위 파싱
    for path, node in _walk_provenanced(raw):
        top = path.split(".")[0].split("[")[0]
        if top in DERIVED_SECTIONS:
            continue
        for need in ("source", "tier"):
            if need not in node:
                out.append(I("error", path, f"'{need}' 누락 (원칙 2: 모든 숫자는 출처를 갖는다)"))
        if "unit" in node and node["unit"] is not None:
            try:
                Q(1.0, str(node["unit"]))
            except Exception as e:
                out.append(I("error", path, f"단위를 해석할 수 없습니다: {node['unit']} ({e})"))
        t = node.get("tier")
        if t is not None and t not in TIER_MEANING:
            out.append(I("error", path, f"tier 는 0~3 이어야 합니다 (지금 {t})"))

    # ⑤ 유도값 재계산 검증 (있으면)
    out += verify(s)

    # ⑥ tier 승인 게이트 (§12.4)
    low = [n for n, p in {**s.core, **s._extra_provenanced()}.items() if p.tier >= 2]
    if low:
        out.append(I("warn", "tier", f"tier ≥ 2 (미검증) 값 {len(low)}건: {', '.join(low[:6])}"
                                     f"{' …' if len(low) > 6 else ''} — 사람 승인 대상 (§12.4)"))
    return out


def verify(s: PhysicalSystem, rtol: float = 1e-3) -> list:
    """`derived_scales`에 적힌 값을 물성식으로 재계산해 대조 (출처 대신 재현으로 검증)."""
    I = _intake.Issue
    ds = s.raw.get("derived_scales")
    if not ds:
        return []
    if not s.is_spherical:
        # 구 공식(3πηd, kT/γ)이 성립하지 않는다 → 재계산 검증을 하지 않고, 대신
        # **어디서 유도했는지**를 요구한다. Perrin 인자를 bdbot에 올리지는 않았다:
        # 아직 한 케이스에서만 나왔다 (두 번 나오면 그때 — CLAUDE.md 추상화 규칙).
        src = str(ds.get("source", ""))
        if not src:
            return [I("error", "derived_scales",
                      f"형상이 '{s.shape}' 라서 구 공식으로 재계산할 수 없습니다. "
                      "`derived_scales.source` 에 유도 스크립트를 적어 재현 가능하게 하세요.")]
        return [I("info", "derived_scales",
                  f"형상 '{s.shape}' — 구 공식 재계산 검증을 건너뜁니다. "
                  f"근거: {src[:60]}")]
    b = s.bulk()
    if b is None:
        return [I("error", "derived_scales",
                  "물성을 재계산할 수 없어 유도값을 대조하지 못했습니다 "
                  "(위의 단위/노드 오류를 먼저 고치세요).")]
    want = {"gamma": b["gamma"], "D_t": b["D_t"], "tau_B": b["tau_B"], "tau_p": b["tau_p"],
            "kT": b["kT"]}
    out = []
    for k, expect in want.items():
        if k not in ds or not isinstance(ds[k], dict):
            continue
        try:
            got = Q(ds[k]["value"], ds[k]["unit"])
            rel = abs(float((got - expect).to(expect.units).magnitude)
                      / float(expect.magnitude))
        except Exception as e:
            out.append(I("error", f"derived_scales.{k}", f"대조 실패: {e}"))
            continue
        if not math.isfinite(rel) or rel > rtol:
            out.append(I("error", f"derived_scales.{k}",
                         f"재계산과 불일치: 적힌 값 {got:~.5gP} vs 계산 {expect.to(got.units):~.5gP} "
                         f"({100*rel:.3f}%)"))
    return out


def render_check(s: PhysicalSystem) -> str:
    L: list[str] = []
    w = L.append
    w("=" * 78)
    w(f"system check — {s.path.parent.name}")
    w("=" * 78)
    if not s.raw:
        w("\n".join(str(i) for i in s.issues))
        return "\n".join(L)

    n_err = len(s.errors)
    n_warn = len([i for i in s.issues if i.level == "warn"])
    w(f"{s.label}   {s.dim}D   스키마: 오류 {n_err} · 경고 {n_warn}")
    if s.raw.get("derived_from"):
        w(f"  근거(L0): {s.raw['derived_from']}")
    if s.issues:
        w("")
        for i in s.issues:
            w(str(i))

    w("")
    w("물리계 (SI)")
    for name, p in s.core.items():
        w(f"  {name:<8} = {str(f'{p.value:~.4gP}'):<18} [tier {p.tier}] {p.source[:40]}")
    extra = s._extra_provenanced()
    for path, p in extra.items():
        w(f"  {path:<28} = {str(f'{p.value:~.4gP}')[:14]:<14} [tier {p.tier}]")

    b = s.bulk()
    if b is not None:
        w("")
        w("유도 물성 (재계산)")
        w(f"  γ = {b['gamma']:~.4eP}   D_t = {b['D_t'].to('um^2/s'):~.4fP}   "
          f"τ_B = {b['tau_B']:~.4gP}   τ_p = {b['tau_p'].to('us'):~.3fP}")

    w("")
    w("신뢰등급 분포 (§12.4)")
    for t, names in s.tiers().items():
        w(f"  tier {t} ({TIER_MEANING[t]:<22}) {len(names):>2}건  {', '.join(names[:5])}"
          f"{' …' if len(names) > 5 else ''}")

    nv = s.raw.get("not_verified")
    if nv:
        w("")
        w("확인하지 않은 것 (원칙 7)")
        for x in nv:
            w(f"  · {str(x).splitlines()[0][:70]}")

    w("")
    w("=" * 78)
    if n_err:
        w(f"VERDICT: FAIL — 오류 {n_err}건. L3(무차원화)로 넘어가지 않습니다.")
    else:
        w("VERDICT: READY — L3(무차원화) 가능.")
        if n_warn:
            w(f"         (경고 {n_warn}건 — tier 승인 대상 확인)")
    w("=" * 78)
    return "\n".join(L)


__all__ = ["SCHEMA", "PhysicalSystem", "load", "validate", "verify", "render_check",
           "REQUIRED_TOP", "OPTIONAL_TOP", "CORE_PROVENANCED", "DERIVED_SECTIONS",
           "TIER_MEANING"]
