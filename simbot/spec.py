"""S3 — 시스템 명세 데이터 모델 + 타당성 검사. LLM 0줄.

이 모듈이 지키는 것 셋:

1. **모든 물리량에 provenance 가 있다.** `Quantity` 없이 맨 float 을 넣을 수 없다.
   근거 없는 숫자는 나중에 "이 값이 어디서 왔는가"에 답할 수 없고, 그러면 감도
   분석도 재현도 불가능해진다.

2. **파생값은 저장하지 않고 다시 계산한다.** `derived:` 블록이 파일에 있으면
   재계산과 대조해서 불일치를 잡는다. 손으로 고친 파생값은 조용히 틀린다 —
   2026-07-28 에 `kT(293.15 K)` 4번째 자리 오류를 실제로 겪었다.

3. **게이트를 끄려면 이유를 적어야 한다.** 그리고 게이트 이름은 등록된 것만
   쓸 수 있다 — 오타 난 게이트 이름은 **한 번도 실행되지 않는 검사**가 된다.
   어느 게이트가 켜지는지는 (계 × 목적동역학) 카드가 정한다:
   `knowledge/wiki/systems/_index.md`
"""
from __future__ import annotations

import math
from dataclasses import MISSING as _MISSING
from dataclasses import dataclass, field, fields, is_dataclass
from pathlib import Path

import yaml

from .io import sha256_payload
from .units import kT_si, stokes_drag_si, stokes_einstein_D_si
from bdbot.constants import sphere_mass_si   # one definition, shared with bdbot.materials

# =============================================================================
# provenance
# =============================================================================
PROVENANCE_KINDS: frozenset[str] = frozenset({
    "from_drawing",    # 그림에서 직접 읽음
    "observation",     # 그림/자료에서 직접 읽음 (from_drawing 의 모달리티 일반형)
    "inference",       # 자료 + 물리지식으로 유도
    "assumed",         # 자료에 없어 채움 → S7b 감도 분석 대상
    "derived",         # 다른 필드에서 계산 (simbot 함수 호출 결과)
    "rule",            # 정책에서 유도 (config/run_policy.yaml)
    "from_knowledge",  # knowledge/wiki 항목
    "from_paper",      # knowledge/source/papers 증류
    "measured",        # 실험값
    "user",            # 사람이 이 런에서 직접 지정 (세션 `set`). 가정이 아니다 —
                       # S7b 감도 분석이 흔들 대상에서 제외된다
})

# master_plan §12.2 — 이 provenance 는 값싼 모델이 채울 수 없다
LLM_RESTRICTED: frozenset[str] = frozenset({"inference", "assumed"})
CHEAP_MODELS: frozenset[str] = frozenset({"haiku", "sonnet"})

CONFIDENCE_LEVELS: frozenset[str] = frozenset({"high", "medium", "low", ""})


@dataclass
class Quantity:
    """물리량 래퍼. **값 하나에 근거 하나.**

    `value` 는 float 이 기본이지만 정수·문자열·리스트도 허용한다
    (`dim=2`, `boundary="periodic"`, `center=[0,0,0]`).
    """

    value: float | int | str | bool | list
    unit: str = ""
    provenance: str = "assumed"
    basis: str = ""
    confidence: str = ""
    ambiguity: str = ""          # 01_intake.md 의 모호성 id (A1, A2 …)
    sensitivity: str = ""        # none | low | high — S7b 결과를 되기록
    affects: list[str] = field(default_factory=list)
    written_by: str = ""         # 모델 티어링 검사용 (§12.2)

    def __post_init__(self) -> None:
        if self.provenance not in PROVENANCE_KINDS:
            raise ValueError(
                f"provenance {self.provenance!r} 는 등록되지 않았다. "
                f"허용: {sorted(PROVENANCE_KINDS)}")
        if self.confidence not in CONFIDENCE_LEVELS:
            raise ValueError(f"confidence {self.confidence!r} 는 "
                             f"{sorted(CONFIDENCE_LEVELS)} 중 하나여야 한다")

    @property
    def si(self) -> float:
        """수치 값. 문자열/리스트면 예외 — 단위 산술에 쓰지 못하게 막는다."""
        if isinstance(self.value, bool) or not isinstance(self.value, (int, float)):
            raise TypeError(f"{self.value!r} 은 수치가 아니다 (unit={self.unit!r})")
        return float(self.value)

    def problems(self) -> list[str]:
        """이 값 하나에 대한 규약 위반 목록."""
        out = []
        if not str(self.basis).strip():
            out.append("basis 가 비어 있다 — 근거 없는 숫자는 감도 분석도 재현도 불가")
        if self.provenance in LLM_RESTRICTED and not self.confidence:
            out.append(f"provenance={self.provenance} 인데 confidence 가 없다 "
                       f"(추론·가정은 신뢰도를 밝혀야 한다)")
        if (self.provenance in LLM_RESTRICTED
                and self.written_by.lower() in CHEAP_MODELS):
            out.append(f"provenance={self.provenance} 를 {self.written_by} 가 채웠다 — "
                       f"master_plan §12.2 위반 (Opus 만 허용)")
        return out


def Q(value, unit: str = "", provenance: str = "assumed", basis: str = "", **kw):
    """`Quantity` 짧은 생성자. YAML 을 손으로 쓰지 않을 때 쓴다."""
    return Quantity(value=value, unit=unit, provenance=provenance, basis=basis, **kw)


# =============================================================================
# 게이트 — 카드가 켜고 끈다
# =============================================================================
GATE_STATUSES: frozenset[str] = frozenset({"required", "pass", "fail", "off",
                                           "applicable", "unknown"})

# 등록된 게이트 이름. 오타 난 이름은 "한 번도 실행되지 않는 검사"가 되므로 거부한다.
KNOWN_GATES: frozenset[str] = frozenset({
    # 전제
    "overdamped", "stokes_reynolds", "hydrodynamics_neglected",
    # 평형·구조
    "equilibration_detection", "equipartition", "configurational_temperature",
    "self_consistency_D", "polydispersity",
    # 수치
    "em_bias_reproduced", "dt_over_tau_trap", "thermal_displacement",
    "force_displacement", "active_displacement", "step_displacement_vs_sigma",
    # 기하
    "box_much_larger_than_l_trap", "r_cut_le_half_box", "finite_size_L",
    "persistence_length_vs_box", "packing_fraction", "debye_length_consistency",
})


@dataclass
class Gate:
    status: str = "unknown"
    reason: str = ""

    def __post_init__(self) -> None:
        if self.status not in GATE_STATUSES:
            raise ValueError(f"게이트 status {self.status!r} 는 "
                             f"{sorted(GATE_STATUSES)} 중 하나여야 한다")


# =============================================================================
# 구성 요소
# =============================================================================
@dataclass
class Species:
    name: str
    n_simulated: Quantity
    radius_si: Quantity
    density_si: Quantity | None = None
    n_physical: Quantity | None = None
    charge: Quantity | None = None
    active: bool = False

    @property
    def sigma_si(self) -> float:
        """대표 직경. ⚠ `σ = 2a`. 반지름과 혼동하면 모든 시간척도가 2배 틀린다."""
        return 2.0 * self.radius_si.si

    def mass_si(self) -> float | None:
        """Sphere-assumed mass. `None` without a density (no overdamped check then).

        ★ The expression is `bdbot.constants.sphere_mass_si`, shared with
          `bdbot.materials.sphere_mass`. It used to be written here as
          `rho*(4/3)*pi*a**3`, which is the same number on paper and **1 ULP
          different in floating point** from bdbot's `rho*(pi/6)*d**3` -- so the
          two halves could never compare equal. Merged 2026-08-29.
          ⚠ The kernel takes a **diameter**; this class stores a radius.
        """
        if self.density_si is None:
            return None
        return sphere_mass_si(self.density_si.si, 2.0 * self.radius_si.si)


@dataclass
class Medium:
    T_si: Quantity
    eta_si: Quantity
    rho_fluid_si: Quantity | None = None
    species: Quantity | None = None      # "water" 등


@dataclass
class Geometry:
    dim: Quantity
    boundary: Quantity
    box_si: Quantity | None = None       # [Lx, Ly, Lz] 명시
    box_over_ref: Quantity | None = None  # 카드 기준길이의 배수 (트랩 계 관행)

    @property
    def d(self) -> int:
        return int(self.dim.value)


@dataclass
class PairInteraction:
    type_a: str
    type_b: str
    potential: str                       # wca | lj | yukawa | morse | dlvo
    params: dict = field(default_factory=dict)
    r_cut_si: Quantity | None = None


@dataclass
class BondInteraction:
    """결합 1종. HOOMD `md.bond.Harmonic` 대응.

    `params`: `k_si` [N/m = J/m²] · `r0_si` [m].

    ⚠ **각(`AngleInteraction`)과 한 클래스로 합치지 않는다.** 결합의 `k` 는 N/m 이고
      각의 `k` 는 J/rad² 다 — 같은 이름의 필드에 다른 차원을 담으면 `λ_max` 를 만들 때
      조용히 틀린다. 안정성 게이트가 그 값을 쓴다.

    ⚠ 거리 구속(`constrain.Distance`)은 `Brownian` 과 함께 쓰지 않는다.
      근거: `knowledge/wiki/findings/dead-end-distance-constraint-with-brownian.md`
    """

    name: str = "backbone"               # HOOMD bond type 이름
    potential: str = "harmonic"
    params: dict = field(default_factory=dict)


@dataclass
class AngleInteraction:
    """각 1종. HOOMD `md.angle.Harmonic` 대응.

    `params`: `k_si` [J/rad²] · `t0_rad` [rad].
    """

    name: str = "backbone"               # HOOMD angle type 이름
    potential: str = "harmonic"
    params: dict = field(default_factory=dict)


@dataclass
class ExternalField:
    kind: str                            # harmonic_trap | gravity | shear | electric
    params: dict[str, Quantity] = field(default_factory=dict)
    implementation: str = ""
    note: str = ""


@dataclass
class Friction:
    model: str = "stokes_infinite_medium"
    gamma_si: Quantity | None = None     # 없으면 6πηa 로 파생
    wall_correction: str = "none"
    note: str = ""


@dataclass
class Timing:
    equil_in_tau: Quantity | None = None
    prod_in_tau: Quantity | None = None
    sample_interval_in_tau: Quantity | None = None
    target_precision: Quantity | None = None


@dataclass
class Numerics:
    dt_star: Quantity | None = None
    seed_base: int = 1
    n_seeds: Quantity | None = None
    integrator: str = "hoomd.md.methods.Brownian"
    scheme: str = "euler_maruyama"
    noise_distribution: str = "uniform"  # ★ Gaussian 아님. findings §2


# =============================================================================
# SystemSpec
# =============================================================================
@dataclass
class SystemSpec:
    """물리 단위 완전 명세. **입력만 담는다** — 파생값은 `derive()` 가 만든다."""

    card: str
    question: str
    geometry: Geometry
    species: list[Species]
    medium: Medium
    friction: Friction = field(default_factory=Friction)
    pair: list[PairInteraction] = field(default_factory=list)
    bonds: list[BondInteraction] = field(default_factory=list)
    angles: list[AngleInteraction] = field(default_factory=list)
    external: list[ExternalField] = field(default_factory=list)
    timing: Timing = field(default_factory=Timing)
    numerics: Numerics = field(default_factory=Numerics)
    gates: dict[str, Gate] = field(default_factory=dict)
    tier: str = ""
    notes: list[str] = field(default_factory=list)

    # --- 편의 접근 ---
    @property
    def primary(self) -> Species:
        return self.species[0]

    def trap(self) -> ExternalField | None:
        for e in self.external:
            if e.kind == "harmonic_trap":
                return e
        return None

    @property
    def has_neighbor_interaction(self) -> bool:
        """겹칠 상대(쌍) 또는 **결합 상대**(결합·각)가 있는가.

        ★ 변위 게이트의 활성 조건이다. `bool(spec.pair)` 로 판정하면 결합만 있는 계
          (콜로이드 사슬)에서 게이트가 조용히 꺼진다 — 실제로 그렇게 꺼져 있었다.
          근거: `knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md`
        """
        return bool(self.pair or self.bonds or self.angles)

    def bond_stiffness_si(self) -> float | None:
        """가장 강한 결합 스프링 [N/m]. **최댓값**이어야 한다 — 안정성은 최악 모드가 정한다."""
        ks = [b.params["k_si"].si for b in self.bonds if "k_si" in b.params]
        return max(ks) if ks else None

    def angle_stiffness_si(self) -> float | None:
        """가장 강한 각 스프링 [J/rad²]."""
        ks = [a.params["k_si"].si for a in self.angles if "k_si" in a.params]
        return max(ks) if ks else None

    def bond_length_si(self) -> float | None:
        """가장 짧은 결합 평형길이 [m]. 각 강성을 횡방향 강성으로 환산할 때의 지레팔이고,
        짧을수록 강성이 커지므로 최솟값을 쓴다."""
        rs = [b.params["r0_si"].si for b in self.bonds if "r0_si" in b.params]
        return min(rs) if rs else None

    def gamma_si(self) -> float:
        """항력계수. 명시값이 있으면 그것, 없으면 Stokes `6πηa`."""
        if self.friction.gamma_si is not None:
            return self.friction.gamma_si.si
        return stokes_drag_si(self.medium.eta_si.si, self.primary.radius_si.si)

    def box_lengths_si(self, ref_length_si: float | None = None) -> list[float] | None:
        """박스 변 길이 [m]. `box_over_ref` 만 있으면 `ref_length_si` 가 필요하다."""
        if self.geometry.box_si is not None:
            return [float(x) for x in self.geometry.box_si.value]
        if self.geometry.box_over_ref is not None and ref_length_si is not None:
            L = self.geometry.box_over_ref.si * ref_length_si
            d = self.geometry.d
            return [L, L, (L if d == 3 else 0.0)]
        return None

    def hash(self) -> str:
        """spec 해시 — `run_id` 와 캐시 키. 파생값은 포함하지 않는다."""
        return sha256_payload(to_dict(self))

    # --- 직렬화 ---
    def to_yaml(self) -> str:
        return dump_yaml(self)

    @classmethod
    def from_yaml(cls, text: str) -> SystemSpec:
        return from_dict(cls, yaml.safe_load(text))

    @classmethod
    def load(cls, path: str | Path) -> SystemSpec:
        return cls.from_yaml(Path(path).read_text(encoding="utf-8"))

    def save(self, path: str | Path) -> Path:
        p = Path(path)
        p.write_text(self.to_yaml(), encoding="utf-8")
        return p


# =============================================================================
# 예측 (S2) — 봉인되는 과학적 주장
# =============================================================================
@dataclass
class PredictionItem:
    """반증 가능한 형태의 예측 1개. 4요소가 전부 있어야 한다 (master_plan §S2-4)."""

    quantity: str
    value: float | str
    tolerance: str                  # "±1.5%" | "±0.03" | "p>0.05" | ">0.99"
    basis: str
    discriminates: str = ""
    unit: str = ""
    competing_value: float | None = None   # 경쟁 가설 (검정력 계산용)
    note: str = ""

    def problems(self) -> list[str]:
        out = []
        for name in ("tolerance", "basis"):
            if not str(getattr(self, name)).strip():
                out.append(f"{self.quantity}: {name} 가 비어 있다")
        return out


@dataclass
class Prediction:
    items: list[PredictionItem]
    regimes: dict = field(default_factory=dict)
    alternatives: list[str] = field(default_factory=list)

    def problems(self) -> list[str]:
        if not self.items:
            return ["정량 예측이 0개다 — S2 게이트는 최소 1개를 요구한다"]
        return [p for it in self.items for p in it.problems()]


# =============================================================================
# 파생값 — 저장하지 않고 계산한다
# =============================================================================
def derive(spec: SystemSpec) -> dict[str, float]:
    """spec 에서 파생되는 SI 스케일 전량. **모든 값이 함수 호출 결과다.**

    여기 없는 파생값을 리포트에 쓰지 않는다 — 손계산이 끼어들 자리를 없앤다.
    """
    sp = spec.primary
    T = spec.medium.T_si.si
    eta = spec.medium.eta_si.si
    a = sp.radius_si.si
    sigma = sp.sigma_si
    gamma = spec.gamma_si()
    kT = kT_si(T)
    D0 = stokes_einstein_D_si(T, gamma)

    out = {
        "kT_si": kT,
        "sigma_si": sigma,
        "gamma_si": gamma,
        "D0_si": D0,
        "tau_D_si": sigma**2 / D0,
    }
    mass = sp.mass_si()
    if mass is not None:
        out["mass_si"] = mass
        out["tau_inertial_si"] = mass / gamma

    trap = spec.trap()
    if trap is not None and "k_si" in trap.params:
        k = trap.params["k_si"].si
        out["k_si"] = k
        out["tau_trap_si"] = gamma / k
        out["l_trap_si"] = math.sqrt(kT / k)
        out["corner_freq_si"] = k / (2.0 * math.pi * gamma)
        out["var_per_component_si"] = kT / k
        out["msd_plateau_si"] = 2.0 * spec.geometry.d * kT / k
        out["k_star_sigma"] = k * sigma**2 / kT

    # --- 결합·각 → 강성행렬 최대고유값 λ_max ---
    #  사슬 강성행렬의 최대고유값 근사: 1D 스프링 사슬이 4k, 굽힘은 4계 차분이라 16k.
    #  각 스프링은 J/rad² 이므로 지레팔 b² 로 나눠 횡방향 강성 [N/m] 으로 환산한다.
    #  이 값이 **명시적 오일러의 안정 한계** Δt ≤ 2γ/λ_max 를 정한다. 정확도가 아니라
    #  발산 여부를 정하는 양이고, 변위 게이트로는 잡히지 않는다 (직선 사슬은 |F| = 0).
    #  실측 교정 (비율 1.22–2.80):
    #  knowledge/wiki/findings/dt-gate-needs-a-stability-term-for-stiff-bonds.md
    k_bond = spec.bond_stiffness_si()
    k_angle = spec.angle_stiffness_si()
    if k_bond is not None or k_angle is not None:
        b = spec.bond_length_si() or sigma
        out["bond_length_si"] = b
        lam = 0.0
        if k_bond is not None:
            out["k_bond_si"] = k_bond
            out["k_bond_star_sigma"] = k_bond * sigma**2 / kT
            out["tau_bond_si"] = gamma / k_bond
            lam += 4.0 * k_bond
        if k_angle is not None:
            out["k_angle_si"] = k_angle
            out["k_angle_star"] = k_angle / kT        # J/rad² 는 길이를 안 쓴다
            out["tau_angle_si"] = gamma * b**2 / k_angle
            lam += 16.0 * k_angle / b**2
        out["lambda_max_si"] = lam
        out["tau_stiff_si"] = gamma / lam
    return out


def reference_length_si(spec: SystemSpec, derived: dict | None = None) -> float:
    """카드의 기준 길이 [m]. 트랩 계면 `ℓ_trap`, 그 외 `σ`.

    (계 × 목적동역학) 카드가 소유하는 선택이다 — CLAUDE.md §무차원화 규약.
    """
    d = derived if derived is not None else derive(spec)
    return d.get("l_trap_si", d["sigma_si"])


# =============================================================================
# 타당성 검사
# =============================================================================
@dataclass
class Check:
    """검사 1건.

    `declared` 는 **결과가 아니다** — 카드가 이 게이트를 켜라고 선언했지만
    S3 이 계산할 수 없는 양이라는 뜻이다 (등분배·EM 편향 등은 S7 이 판정한다).
    `declared` 를 `pass` 로 적으면 사람이 한 번도 보지 않은 합격 도장이 찍힌다.
    """

    name: str
    status: str            # pass | fail | off | na | warn | declared
    detail: str = ""
    value: float | None = None
    threshold: float | None = None


@dataclass
class SpecReport:
    """검사 결과. **판정하지 않는다** — `ok` 는 규약 위반 유무일 뿐이다."""

    checks: list[Check]
    problems: list[str]
    derived: dict[str, float]

    @property
    def ok(self) -> bool:
        return not self.problems and not any(c.status == "fail" for c in self.checks)

    def failed(self) -> list[Check]:
        return [c for c in self.checks if c.status == "fail"]

    def deferred(self) -> list[Check]:
        """S7 이 판정해야 하는 게이트 — S3 에서는 계산할 수 없다."""
        return [c for c in self.checks if c.status == "declared"]

    def table(self) -> str:
        """마크다운 표. `03_spec_rationale.md` 와 `REPORT.md` 가 쓴다."""
        rows = ["| 검사 | 상태 | 값 | 문턱 | 비고 |", "|---|---|---|---|---|"]
        mark = {"pass": "✅ pass", "fail": "❌ **fail**", "off": "— off",
                "na": "— n/a", "warn": "⚠️ warn", "declared": "⏳ S7 판정"}
        for c in self.checks:
            v = "" if c.value is None else f"`{c.value:.4g}`"
            t = "" if c.threshold is None else f"`{c.threshold:.4g}`"
            rows.append(f"| `{c.name}` | {mark.get(c.status, c.status)} | {v} | {t} "
                        f"| {c.detail} |")
        return "\n".join(rows)


def _iter_quantities(obj, path: str = "") -> list[tuple[str, Quantity]]:
    """중첩 구조에서 모든 `Quantity` 를 (경로, 값) 으로 뽑는다."""
    out: list[tuple[str, Quantity]] = []
    if isinstance(obj, Quantity):
        out.append((path, obj))
    elif is_dataclass(obj):
        for f in fields(obj):
            out += _iter_quantities(getattr(obj, f.name), f"{path}.{f.name}".lstrip("."))
    elif isinstance(obj, dict):
        for k, v in obj.items():
            out += _iter_quantities(v, f"{path}.{k}".lstrip("."))
    elif isinstance(obj, (list, tuple)):
        for i, v in enumerate(obj):
            out += _iter_quantities(v, f"{path}[{i}]")
    return out


def packing_fraction(spec: SystemSpec, box_si: list[float] | None) -> float | None:
    """φ (3D 부피분율 / 2D 면적분율). 박스를 모르면 `None`."""
    if not box_si:
        return None
    d = spec.geometry.d
    if d == 2:
        area = box_si[0] * box_si[1]
        tot = sum(s.n_simulated.si * math.pi * s.radius_si.si**2 for s in spec.species)
        return tot / area if area > 0 else None
    vol = box_si[0] * box_si[1] * box_si[2]
    tot = sum(s.n_simulated.si * (4.0 / 3.0) * math.pi * s.radius_si.si**3
              for s in spec.species)
    return tot / vol if vol > 0 else None


def validate(spec: SystemSpec, *, stored_derived: dict | None = None,
             rel_tol: float = 1e-3) -> SpecReport:
    """규약 위반 + 물리적 타당성 검사.

    Args:
        stored_derived: 파일에 적혀 있던 파생값. 주면 재계산과 대조한다
            (손으로 고친 파생값을 잡는 유일한 방법).
        rel_tol: 파생값 대조 허용 상대오차. 문서에 유효숫자 4~5자리로 적히므로
            기본 `1e-3`.
    """
    problems: list[str] = []
    d = derive(spec)
    computed: dict[str, Check] = {}

    def put(c: Check) -> None:
        computed[c.name] = c

    # --- 1. provenance 완전성 -------------------------------------------------
    for path, q in _iter_quantities(spec):
        problems += [f"{path}: {p}" for p in q.problems()]

    # --- 2. 게이트 선언 형식 --------------------------------------------------
    unknown = sorted(set(spec.gates) - KNOWN_GATES)
    if unknown:
        problems.append(
            f"등록되지 않은 게이트 이름 {unknown} — 오타는 '한 번도 실행되지 않는 "
            f"검사'가 된다. 새 게이트면 KNOWN_GATES 에 먼저 추가할 것")
    for name, g in spec.gates.items():
        if g.status == "off" and not g.reason.strip():
            problems.append(f"게이트 {name} 를 껐는데 이유가 없다 — "
                            f"끄는 근거는 카드에 있어야 한다")

    # --- 3. 과감쇠 -----------------------------------------------------------
    tau_proc = d.get("tau_trap_si", d["tau_D_si"])
    if "tau_inertial_si" in d:
        ratio = d["tau_inertial_si"] / tau_proc
        put(Check("overdamped", "pass" if ratio < 1e-2 else "fail", value=ratio,
                  threshold=1e-2,
                  detail=("τ_i/τ_process ≪ 1 — BD 전제 유효" if ratio < 1e-2
                          else "관성이 무시되지 않는다 → Langevin 검토")))
    else:
        put(Check("overdamped", "na", detail="밀도가 없어 질량을 모른다 — 검사 불가"))

    # --- 4. Reynolds --------------------------------------------------------
    rho_f = spec.medium.rho_fluid_si
    if rho_f is not None:
        # 특성 속도: 기준 길이를 기준 시간에 지나는 속도
        v = reference_length_si(spec, d) / tau_proc
        Re = rho_f.si * v * spec.primary.radius_si.si / spec.medium.eta_si.si
        put(Check("stokes_reynolds", "pass" if Re < 1e-2 else "fail",
                  value=Re, threshold=1e-2,
                  detail="Re ≪ 1 — Stokes 항력 유효" if Re < 1e-2
                  else "관성 유체효과 의심"))
    else:
        put(Check("stokes_reynolds", "na", detail="ρ_fluid 없음"))

    # --- 5. 박스 · φ · r_cut ------------------------------------------------
    ref = reference_length_si(spec, d)
    box = spec.box_lengths_si(ref)
    if box is None:
        problems.append("박스 크기를 알 수 없다 — box_si 또는 box_over_ref 필요")
    else:
        n_ref = min(x for x in box[:spec.geometry.d]) / ref
        put(Check("box_much_larger_than_l_trap",
                  "pass" if n_ref >= 10 else "fail", value=n_ref, threshold=10.0,
                  detail=f"박스가 기준길이의 {n_ref:.3g}배"))
        phi = packing_fraction(spec, box)
        if phi is not None:
            cap = 0.9 if spec.geometry.d == 2 else 0.64
            kind = "면적" if spec.geometry.d == 2 else "부피"
            if not spec.pair:
                # ★ 쌍 상호작용이 없으면 배제부피가 없다 → φ 는 물리적 의미가 없다.
                #   N 개 입자는 같은 트랩 안의 **독립 복제**이지 서스펜션이 아니다.
                #   여기서 φ 를 게이트로 걸면 φ=4741 로 통과 불가 판정이 나온다 —
                #   존재하지 않는 문제다. 값은 보여주되 판정하지 않는다.
                put(Check("packing_fraction", "off", value=phi, threshold=cap,
                          detail=f"{kind}분율은 의미 없음 — 쌍 상호작용이 없어 "
                                 f"배제부피가 없다 (입자는 독립 복제)"))
            else:
                put(Check("packing_fraction", "pass" if phi < cap else "fail",
                          value=phi, threshold=cap,
                          detail=f"{kind}분율 (상한 = "
                                 f"{'2D 한계' if spec.geometry.d == 2 else 'RCP'})"))
        if spec.pair:
            half = min(box[:spec.geometry.d]) / 2.0
            worst = max((p.r_cut_si.si for p in spec.pair
                         if p.r_cut_si is not None), default=None)
            if worst is None:
                problems.append("쌍 상호작용이 있는데 r_cut 이 없다")
            else:
                put(Check("r_cut_le_half_box", "pass" if worst <= half else "fail",
                          value=worst, threshold=half, detail="최소이미지 규약"))
        else:
            put(Check("r_cut_le_half_box", "off", detail="쌍 포텐셜 없음"))

    # --- 6. 파생값 대조 ------------------------------------------------------
    if stored_derived:
        bad, n_compared, skipped = [], 0, []
        for k, stored in stored_derived.items():
            if k not in d or not isinstance(stored, (int, float)):
                skipped.append(k)
                continue
            n_compared += 1
            rel = abs(float(stored) - d[k]) / max(abs(d[k]), 1e-300)
            if rel > rel_tol:
                bad.append(f"{k}: 파일 {stored:.6g} vs 재계산 {d[k]:.6g} "
                           f"(상대차 {rel:.2e})")
        problems += [f"파생값 불일치 — {b} — 손으로 고친 값일 가능성" for b in bad]
        # ★ 대조한 개수를 보고한다. 0건 대조를 pass 로 보고하면 "검사했다"가 거짓이 된다.
        detail = f"{n_compared}개 대조"
        if bad:
            detail += f" · {len(bad)}개 불일치"
        if skipped:
            detail += f" · 대조 불가 {len(skipped)}개 ({', '.join(skipped[:3])}…)" \
                if len(skipped) > 3 else f" · 대조 불가 {skipped}"
        put(Check("derived_consistency",
                  "fail" if bad else ("na" if n_compared == 0 else "pass"),
                  value=float(n_compared), detail=detail))

    # --- 7. 게이트 ∪ 계산결과 --------------------------------------------------
    # 계산된 결과가 선언을 채운다. 계산할 수 없는 게이트는 `declared` 로 남아
    # S7 로 넘어간다 — 여기서 pass 를 찍으면 아무도 보지 않은 합격이 된다.
    checks: list[Check] = []
    for name, g in spec.gates.items():
        if name in computed:
            c = computed.pop(name)
            parts = [x for x in (c.detail, g.reason) if x]
            if len(parts) == 2 and parts[0] == parts[1]:
                parts = parts[:1]
            checks.append(Check(name, c.status, value=c.value,
                                threshold=c.threshold, detail=" · ".join(parts)))
        elif g.status == "off":
            checks.append(Check(name, "off", detail=g.reason))
        else:
            checks.append(Check(name, "declared",
                                detail=g.reason or f"카드 선언: {g.status}"))
    # 카드가 선언하지 않았는데 계산된 것들 (선언 누락을 드러낸다)
    for c in computed.values():
        checks.append(c)

    return SpecReport(checks=checks, problems=problems, derived=d)


# =============================================================================
# YAML 직렬화 — 왕복 오차 0 이어야 한다
# =============================================================================
_QUANTITY_DEFAULTS = {f.name: f.default for f in fields(Quantity)
                      if f.name not in ("value",)}


def to_dict(obj):
    """dataclass 트리 → 순수 dict. `Quantity` 는 기본값 필드를 생략해 짧게 쓴다."""
    if isinstance(obj, Quantity):
        # provenance 는 기본값이어도 **항상 적는다.** `assumed` 가 기본값이라는 이유로
        # 생략하면 "가정했다"와 "적기를 잊었다"가 파일에서 구별되지 않는다 —
        # 그런데 `assumed` 는 S7b 감도 분석의 대상 목록을 정하는 필드다.
        out = {"value": obj.value, "provenance": obj.provenance}
        if obj.unit:
            out = {"value": obj.value, "unit": obj.unit,
                   "provenance": obj.provenance}
        for name, default in _QUANTITY_DEFAULTS.items():
            if name in ("unit", "provenance"):
                continue
            v = getattr(obj, name)
            if name == "affects":
                if v:
                    out[name] = list(v)
                continue
            if v != default:
                out[name] = v
        return out
    if isinstance(obj, Gate):
        return {"status": obj.status, **({"reason": obj.reason} if obj.reason else {})}
    if is_dataclass(obj):
        out = {}
        for f in fields(obj):
            v = getattr(obj, f.name)
            if v is None:
                continue
            default = f.default if f.default is not _MISSING else None
            if isinstance(v, (list, dict)) and not v:
                continue
            if default is not None and v == default and not isinstance(v, Quantity):
                continue
            out[f.name] = to_dict(v)
        return out
    if isinstance(obj, dict):
        return {k: to_dict(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [to_dict(v) for v in obj]
    return obj


def _build_quantity(raw) -> Quantity:
    if isinstance(raw, Quantity):
        return raw
    if isinstance(raw, dict) and "value" in raw:
        return Quantity(**raw)
    # 맨 값이 들어온 경우 — provenance 가 없으므로 규약 위반으로 흘려보낸다
    return Quantity(value=raw, provenance="assumed", basis="")


def from_dict(cls, raw):
    """dict → dataclass. `Quantity` 필드는 `{value, provenance, …}` 형태를 받는다."""
    if raw is None:
        return None
    if cls is Quantity:
        return _build_quantity(raw)
    if cls is Gate:
        return Gate(**raw) if isinstance(raw, dict) else Gate(status=str(raw))
    if not is_dataclass(cls):
        return raw

    kwargs = {}
    for f in fields(cls):
        if f.name not in raw:
            continue
        v = raw[f.name]
        ann = f.type if not isinstance(f.type, str) else f.type
        kwargs[f.name] = _coerce(f.name, ann, v, cls)
    return cls(**kwargs)


# 필드 이름 → 원소 타입. 문자열 annotation 을 파싱하지 않고 명시 표로 둔다.
_LIST_FIELD_TYPES = {
    ("SystemSpec", "species"): Species,
    ("SystemSpec", "pair"): PairInteraction,
    ("SystemSpec", "bonds"): BondInteraction,
    ("SystemSpec", "angles"): AngleInteraction,
    ("SystemSpec", "external"): ExternalField,
    ("Prediction", "items"): PredictionItem,
}
_NESTED_FIELD_TYPES = {
    ("SystemSpec", "geometry"): Geometry,
    ("SystemSpec", "medium"): Medium,
    ("SystemSpec", "friction"): Friction,
    ("SystemSpec", "timing"): Timing,
    ("SystemSpec", "numerics"): Numerics,
}
# Quantity 로 감싸야 하는 필드 (이름 기준). `params` 안의 값도 전부 Quantity 다.
_QUANTITY_FIELDS = {
    "n_simulated", "n_physical", "radius_si", "density_si", "charge",
    "T_si", "eta_si", "rho_fluid_si", "species", "dim", "boundary", "box_si",
    "box_over_ref", "gamma_si", "r_cut_si", "dt_star", "n_seeds",
    "equil_in_tau", "prod_in_tau", "sample_interval_in_tau", "target_precision",
}


def _coerce(name: str, ann, v, owner):
    key = (owner.__name__, name)
    if key in _LIST_FIELD_TYPES:
        return [from_dict(_LIST_FIELD_TYPES[key], x) for x in v]
    if key in _NESTED_FIELD_TYPES:
        return from_dict(_NESTED_FIELD_TYPES[key], v)
    if name == "gates":
        return {k: from_dict(Gate, g) for k, g in v.items()}
    if name == "params" and isinstance(v, dict):
        # 트랩의 k_si 등은 Quantity, active_axes 같은 순수 설정은 그대로
        return {k: (_build_quantity(x) if isinstance(x, dict) and "value" in x else x)
                for k, x in v.items()}
    if name in _QUANTITY_FIELDS:
        return _build_quantity(v)
    return v


class _Dumper(yaml.SafeDumper):
    """들여쓰기를 사람이 읽기 좋게. 손으로 쓴 03_spec.yaml 과 같은 모양."""

    def increase_indent(self, flow=False, indentless=False):
        return super().increase_indent(flow, False)


def dump_yaml(obj) -> str:
    return yaml.dump(to_dict(obj), Dumper=_Dumper, sort_keys=False,
                     allow_unicode=True, default_flow_style=False, width=100)


def load_prediction(path: str | Path) -> Prediction:
    raw = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    return from_dict(Prediction, raw)
