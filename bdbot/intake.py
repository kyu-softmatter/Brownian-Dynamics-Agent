"""L0 인테이크 — `Observation` 스키마와 검사 (마스터플랜 §5.1, §8).

**스키마를 계획에서 추측하지 않고 실제 사용에서 도출했습니다.** 스케치 5장을 손으로
읽어 쓴 `intake/*/observation.yaml` 5개의 필드 사용 빈도가 근거입니다:

    5/5 (필수)  source_files · raw_transcription · system_guess · entities ·
                stated_quantities · stated_goals · ambiguities · unread_regions ·
                missing_required
    2/5 (선택)  prerun_findings · references
    1/5 (선택)  hard_constraints · model_notes · contrast_with_*

하위 키도 같은 방식으로 갈랐습니다. 가장 중요한 발견:

  ★ `resolution` 키가 ambiguities 23/23, missing_required 24/24 **전부**에 있고
    값은 자주 `null` 입니다. 이게 §8.3의 "모르는 걸 모른다고 말하게 하는 장치"입니다.
    → 검사는 **키의 존재를 강제**하고 `null` 값은 허용합니다. 누락은 거부합니다.

  ★ `assumed_value`가 있는데 `confidence`가 없으면 **지어낸 값에 tier가 없는 것**입니다.
    CLAUDE.md 규칙 3 위반이므로 거부합니다.

pydantic을 쓰지 않았습니다 — 이 검사의 가치는 "무엇이 왜 빠졌는지 사람에게 정확히
말해주는 것"이고, 선택/`null` 의미가 섞인 스키마에서는 직접 쓴 진단이 더 낫습니다.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import yaml

SCHEMA = "bdbot.observation/0.1"

# ── 스키마 (실사용 빈도에서 도출) ──────────────────────────────────────
REQUIRED_TOP = (
    "source_files", "raw_transcription", "system_guess", "entities",
    "stated_quantities", "stated_goals", "ambiguities", "unread_regions",
    "missing_required",
)
OPTIONAL_TOP = ("prerun_findings", "references", "hard_constraints", "model_notes")

# 리스트 항목별 (필수 키, 선택 키). 필수는 5개 파일에서 100% 등장한 것만.
ITEM_KEYS = {
    "entities": (("kind",), ("note", "form", "expression", "shape", "center")),
    "stated_quantities": (("symbol", "value", "unit", "source"), ("confidence", "note")),
    "ambiguities": (("id", "issue", "resolution"),
                    ("impact", "note", "lean", "detail", "confirmed_by", "options",
                     "consequence", "evidence")),
    "missing_required": (("symbol", "resolution"),
                         ("assumed_value", "assumed_unit", "confidence", "what", "note",
                          "confirmed_by", "proposed", "followup", "lean", "kind")),
    "prerun_findings": (("id", "finding"), ("severity", "consequence", "verified", "options")),
    "references": (("kind",), ("authors", "doi", "title", "note", "resolved", "locator")),
}
# 이 키들은 "빠뜨리면 안 되는 것"이라 빈 리스트라도 **명시**해야 합니다 (§8.3).
MUST_BE_EXPLICIT = ("ambiguities", "unread_regions", "missing_required")

# `missing_required` 의 종류. ★ 도구를 5개 실제 파일에 돌려보고 필요해진 구분입니다 —
# 판정이 틀렸습니다: 이미 완주한 trap-2d-5um·soft-r3 가 BLOCKED 로 나왔고, 원인은
# `L`(박스 크기)·`T_obs`(관측 창)처럼 **물리적 미지값이 아니라 시뮬레이션 선택 사항**이
# 같은 목록에 섞여 있던 것이었습니다. 손으로 쓴 파일은 이 구분을 note 산문에만 적어뒀습니다.
#   physical  계의 성질. 사람이 주거나 KB에서 찾아야 L2를 쓸 수 없다 → **차단**
#   choice    시뮬레이션 선택 (박스·관측창·표본수). 수치 섹션에서 정한다 → 차단하지 않음
MISSING_KINDS = ("physical", "choice")


@dataclass
class Issue:
    level: str          # error | warn | info
    where: str
    msg: str

    def __str__(self):
        mark = {"error": "✗", "warn": "⚠", "info": "ℹ"}[self.level]
        return f"  {mark} [{self.where}] {self.msg}"


@dataclass
class Observation:
    path: Path
    raw: dict
    issues: list = field(default_factory=list)

    # ── 준비도 ────────────────────────────────────────────────────────
    @property
    def open_ambiguities(self) -> list:
        return [a for a in (self.raw.get("ambiguities") or [])
                if isinstance(a, dict) and a.get("resolution") in (None, "")]

    @property
    def open_missing(self) -> list:
        """L2를 막는 미해소 결측.

        `assumed_value`가 있으면 잠정 진행 가능으로 봅니다.
        `kind: choice`(박스·관측창 등 시뮬레이션 선택)는 차단하지 않습니다 — 위 MISSING_KINDS.
        기본값은 `physical`(보수적: 종류를 안 적었으면 막습니다).
        """
        out = []
        for m in (self.raw.get("missing_required") or []):
            if not isinstance(m, dict):
                continue
            if m.get("kind", "physical") == "choice":
                continue
            if m.get("resolution") in (None, "") and m.get("assumed_value") in (None, ""):
                out.append(m)
        return out

    @property
    def open_choices(self) -> list:
        """미정 시뮬레이션 선택 — 차단하지 않지만 L3에서 정해야 합니다."""
        return [m for m in (self.raw.get("missing_required") or [])
                if isinstance(m, dict) and m.get("kind") == "choice"
                and m.get("resolution") in (None, "") and m.get("assumed_value") in (None, "")]

    @property
    def assumed(self) -> list:
        """가정으로 채운 값 — 리포트에 tier와 함께 반드시 드러나야 합니다 (원칙 2)."""
        return [m for m in (self.raw.get("missing_required") or [])
                if isinstance(m, dict) and m.get("assumed_value") not in (None, "")]

    @property
    def errors(self) -> list:
        return [i for i in self.issues if i.level == "error"]

    @property
    def ready_for_system(self) -> bool:
        """`system.yaml`(L2)을 쓸 수 있는가 — 스키마가 온전하고 미해소 결측이 없어야."""
        return not self.errors and not self.open_missing


def load(path) -> Observation:
    p = Path(path)
    if p.is_dir():
        p = p / "observation.yaml"
    if not p.exists():
        obs = Observation(p, {})
        obs.issues.append(Issue("error", str(p), "observation.yaml 이 없습니다. "
                                "`bdbot intake init <folder>` 로 템플릿을 만드세요."))
        return obs
    raw = yaml.safe_load(p.read_text()) or {}
    obs = Observation(p, raw)
    obs.issues = validate(obs)
    return obs


def validate(obs: Observation) -> list:
    """스키마 검사. `null` 값은 허용하고 **키 누락은 거부**합니다."""
    raw, root = obs.raw, obs.path.parent.parent.parent
    out: list[Issue] = []

    # ① 최상위 필수 필드
    for k in REQUIRED_TOP:
        if k not in raw:
            hint = ("빈 리스트라도 명시해야 합니다 (§8.3: 모르는 걸 모른다고 말하는 장치)"
                    if k in MUST_BE_EXPLICIT else "")
            out.append(Issue("error", k, f"필수 필드 누락. {hint}".strip()))

    # ② 전사가 실제로 채워졌는가
    tr = raw.get("raw_transcription")
    if isinstance(tr, str) and len(tr.strip()) < 20:
        out.append(Issue("error", "raw_transcription",
                         "전사가 비어 있거나 너무 짧습니다. 규칙 5: 해석보다 전사가 먼저입니다."))

    # ③ 소스 파일이 실제로 있는가
    for s in (raw.get("source_files") or []):
        if not (obs.path.parent / Path(s).name).exists():
            out.append(Issue("warn", "source_files", f"파일을 찾을 수 없습니다: {s}"))

    # ④ 리스트 항목의 키
    for key, (req, opt) in ITEM_KEYS.items():
        items = raw.get(key)
        if items is None:
            continue
        if not isinstance(items, list):
            out.append(Issue("error", key, f"리스트여야 합니다 (지금 {type(items).__name__})"))
            continue
        known = set(req) | set(opt)
        seen_ids: set = set()
        for i, it in enumerate(items):
            at = f"{key}[{i}]"
            if not isinstance(it, dict):
                if key in ("unread_regions",):
                    continue
                out.append(Issue("error", at, "매핑이어야 합니다"))
                continue
            for r in req:
                if r not in it:
                    out.append(Issue("error", at, f"필수 키 '{r}' 누락 "
                                                  f"(실사용 5파일에서 100% 등장)"))
            for k2 in it:
                if k2 not in known:
                    out.append(Issue("info", at, f"스키마에 없는 키 '{k2}' — 의도한 것이면 무시"))
            if "id" in it:
                if it["id"] in seen_ids:
                    out.append(Issue("error", at, f"id 중복: {it['id']}"))
                seen_ids.add(it["id"])

    # ⑤ ★ 핵심 방어 — 지어낸 값에 신뢰등급이 없으면 거부 (CLAUDE.md 규칙 3)
    for i, m in enumerate(raw.get("missing_required") or []):
        if not isinstance(m, dict):
            continue
        at = f"missing_required[{i}]  {m.get('symbol', '?')}"
        has_val = m.get("assumed_value") not in (None, "")
        if has_val and m.get("confidence") is None:
            out.append(Issue("error", at,
                             "assumed_value 가 있는데 confidence(tier)가 없습니다. "
                             "출처 없는 값은 넣지 않습니다 (규칙 3)."))
        if has_val and m.get("note") in (None, "") and m.get("resolution") in (None, ""):
            out.append(Issue("warn", at, "가정값의 근거(note)가 비어 있습니다."))
        if m.get("kind") not in (None,) + MISSING_KINDS:
            out.append(Issue("error", at, f"kind 는 {MISSING_KINDS} 중 하나여야 합니다."))
        blocking = (m.get("kind", "physical") != "choice" and not has_val
                    and m.get("resolution") in (None, ""))
        if blocking and not (m.get("what") or m.get("note")):
            out.append(Issue("warn", at, "L2를 막는 항목인데 what/note 가 비어 있습니다 — "
                                         "사람에게 무엇이 필요한지 알려줄 수 없습니다."))

    # ⑥ 스케치에 적힌 값에는 출처가 있어야 한다 (원칙 2)
    for i, q in enumerate(raw.get("stated_quantities") or []):
        if not isinstance(q, dict):
            continue
        at = f"stated_quantities[{i}]  {q.get('symbol', '?')}"
        if not (q.get("source") or "").strip():
            out.append(Issue("error", at, "source 가 비어 있습니다 (원칙 2: 모든 숫자는 출처를 갖는다)."))
        if q.get("value") is None and q.get("unit") is not None:
            out.append(Issue("info", at, "값이 null — 스케치에 기호만 있는 경우로 읽었습니다."))

    # ⑦ 모호점이 하나도 없다고 주장하면 의심한다
    if isinstance(raw.get("ambiguities"), list) and not raw["ambiguities"]:
        out.append(Issue("warn", "ambiguities",
                         "모호점이 0건입니다. 손으로 그린 스케치에서 이런 경우는 드뭅니다 — 다시 보세요."))
    return out


def render_check(obs: Observation) -> str:
    """사람이 읽는 검사 리포트."""
    L: list[str] = []
    w = L.append
    w("=" * 78)
    w(f"intake check — {obs.path.parent.name}")
    w("=" * 78)
    if not obs.raw:
        w("\n".join(str(i) for i in obs.issues))
        return "\n".join(L)

    n_err = len(obs.errors)
    n_warn = len([i for i in obs.issues if i.level == "warn"])
    n_info = len([i for i in obs.issues if i.level == "info"])
    w(f"스키마: 오류 {n_err} · 경고 {n_warn} · 정보 {n_info}")
    if obs.issues:
        w("")
        for i in obs.issues:
            w(str(i))

    w("")
    w("내용 요약")
    w(f"  전사 {len((obs.raw.get('raw_transcription') or '').splitlines())}줄 · "
      f"엔티티 {len(obs.raw.get('entities') or [])} · "
      f"명시 수치 {len(obs.raw.get('stated_quantities') or [])} · "
      f"목표 {len(obs.raw.get('stated_goals') or [])}")
    w(f"  모호점 {len(obs.raw.get('ambiguities') or [])}건 "
      f"(미해소 {len(obs.open_ambiguities)}) · "
      f"미판독 {len(obs.raw.get('unread_regions') or [])}건 · "
      f"결측 {len(obs.raw.get('missing_required') or [])}건 "
      f"(미해소 {len(obs.open_missing)}, 가정으로 채움 {len(obs.assumed)})")

    if obs.assumed:
        w("")
        w("가정으로 채운 값 (tier 확인 대상 — 원칙 2)")
        for m in obs.assumed:
            unit = f" {m.get('assumed_unit')}" if m.get("assumed_unit") else ""
            w(f"  {m.get('symbol', '?'):<16} = {m.get('assumed_value')}{unit}"
              f"   [tier {m.get('confidence')}]  {(m.get('note') or '')[:44]}")

    if obs.open_missing:
        w("")
        w("★ 물리계를 확정할 수 없는 이유 (지어내지 않고 멈춤 — 규칙 3)")
        for m in obs.open_missing:
            desc = m.get("what") or m.get("note") or "(설명 없음 — what 을 채우세요)"
            w(f"  {m.get('symbol', '?'):<16} {str(desc)[:56]}")
    if obs.open_choices:
        w("")
        w("미정 시뮬레이션 선택 (차단하지 않음 — L3에서 정함)")
        for m in obs.open_choices:
            desc = m.get("what") or m.get("proposed") or m.get("note") or ""
            w(f"  {m.get('symbol', '?'):<16} {str(desc)[:56]}")

    if obs.open_ambiguities:
        w("")
        w("미해소 모호점 (사람 확인 #1 대상)")
        for a in obs.open_ambiguities:
            lean = f"  → 기울어짐: {a['lean']}" if a.get("lean") else ""
            w(f"  [{a.get('id', '?')}] {str(a.get('issue', ''))[:58]}{lean}")

    w("")
    w("=" * 78)
    if n_err:
        w(f"VERDICT: FAIL — 스키마 오류 {n_err}건. 고치기 전에는 다음 단계로 넘어가지 않습니다.")
    elif obs.open_missing:
        w(f"VERDICT: BLOCKED — 스키마는 온전하나 결측 {len(obs.open_missing)}건이 미해소입니다.")
        w("         L2(system.yaml)를 쓸 수 없습니다. 사람이 값을 주거나 KB에서 찾아야 합니다.")
    else:
        w("VERDICT: READY — L2(system.yaml) 작성 가능.")
        if obs.open_ambiguities:
            w(f"         (미해소 모호점 {len(obs.open_ambiguities)}건은 결과 해석에 영향)")
    w("=" * 78)
    return "\n".join(L)


TEMPLATE = """\
# L0 Observation — 스케치에서 읽어낸 것. **초안: 사람 확인 대기**
# 프로토콜 (skill `bd-intake`): 전사 먼저 → 구조화 → 모호/미판독 명시 → 없는 값은 null
#
# 검사:  $PY -m bdbot.cli intake check intake/{case}
#   · ambiguities / unread_regions / missing_required 는 **비어 있어도 키를 남겨야** 합니다
#   · assumed_value 를 쓰면 confidence(tier)가 필수입니다 — 출처 없는 값은 넣지 않습니다
#   · resolution 은 null 로 두세요. 모른다는 것을 명시하는 자리입니다

source_files:
  - intake/{case}/sketch_01.jpeg

# ── ① 전사 (보이는 그대로. 해석 안 함) ─────────────────────────────────
raw_transcription: |
  [여기에 이미지에서 읽은 글자를 그대로 옮겨 적으세요.
   읽지 못한 부분은 [판독 불가] 로 표시하고 unread_regions 에도 적으세요.]

# ── ② 구조화 ──────────────────────────────────────────────────────────
system_guess: ""

entities: []          # - kind: particle / pair_interaction / external_potential / box ...
                      #   note: "..."

stated_quantities: []  # - symbol: k_t
                       #   value: 10
                       #   unit: pN/um
                       #   source: "스케치 우측, 명시"
                       #   confidence: 1

stated_goals: []      # 스케치에 적힌 측정 목표. 없으면 빈 리스트로 두고 ambiguities에 적으세요

# ── ③ 모호한 것 (비어 있어도 키는 남길 것) ─────────────────────────────
ambiguities: []       # - id: A1
                      #   issue: "무엇이 모호한가"
                      #   impact: "이 선택이 결과를 어떻게 바꾸는가"
                      #   lean: "내가 기울어진 쪽 (근거와 함께)"
                      #   resolution: null        # ← 사람이 채움

# ── ④ 판독 못 한 부분 ─────────────────────────────────────────────────
unread_regions: []    # - "좌측 하단 두 글자 판독 불가"

# ── ⑤ 스케치에 없어서 채워야 하는 값 (지어내지 않고 출처 표기) ──────────
missing_required: []  # - symbol: eta
                      #   what: "용매 점도"
                      #   assumed_value: 0.851
                      #   assumed_unit: mPa*s
                      #   confidence: 1           # assumed_value 가 있으면 필수
                      #   note: "물@300K 핸드북. 매질은 스케치 미기재"
                      #   resolution: null
"""


def init_template(folder, case: str | None = None, force: bool = False) -> tuple[bool, str]:
    p = Path(folder)
    p.mkdir(parents=True, exist_ok=True)
    target = p / "observation.yaml"
    if target.exists() and not force:
        return False, f"이미 있습니다: {target}  (--force 로 덮어쓰기)"
    target.write_text(TEMPLATE.format(case=case or p.name))
    return True, f"템플릿 생성: {target}"


__all__ = ["SCHEMA", "Observation", "Issue", "load", "validate", "render_check",
           "init_template", "REQUIRED_TOP", "OPTIONAL_TOP", "ITEM_KEYS", "MUST_BE_EXPLICIT",
           "MISSING_KINDS"]
