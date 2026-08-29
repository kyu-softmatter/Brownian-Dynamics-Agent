"""콜로이드 상호작용 카탈로그 + 추천기 — 스케치에 `U_ij` 가 없을 때.

**왜 필요한가**: 스케치 5장 중 2장이 페어 퍼텐셜을 비워뒀습니다 (`chain-bend` 의 물결선,
`trap-drag` 는 아예 없음). CLAUDE.md 규칙 3 때문에 지어낼 수 없고, 그렇다고 "없습니다"로
멈추면 사람이 뭘 골라야 하는지 알 수 없습니다.

→ **추천하고 묻습니다.** 표준 후보를 근거·필요 파라미터·무차원 결합세기와 함께 제시하고,
   `resolution` 에 붙여넣을 YAML까지 만들어 줍니다. 다른 걸 쓰겠다면 사람이 직접 적습니다.
   결정은 여전히 사람이 하고, 도구는 **선택지와 그 결과를 보여주는 역할**만 합니다.

단일 매질(뉴턴 유체) 안의 콜로이드에 대한 표준 골격은 **DLVO** 입니다:
    U = U_전기이중층(스크린된 Coulomb) + U_van der Waals   (+ 배제부피)
여기에 계에 따라 고갈 인력·장 유도 쌍극자·접촉 접선력이 붙습니다.

⚠️ 함수 형태는 표준이지만 **구체적 수치(Hamaker 상수·유효전하·λ_D)는 핸드북/문헌을 봐야
   합니다.** 카탈로그는 값을 주지 않고 "무엇을 알아야 하는지"를 줍니다 (원칙 2).
"""
from __future__ import annotations

from dataclasses import dataclass, field

SCHEMA = "bdbot.interactions/0.1"


@dataclass
class Interaction:
    key: str
    name: str
    form: str
    needs: tuple                 # 지정해야 하는 차원 있는 파라미터
    coupling: str                # 이 상호작용이 만드는 무차원 결합세기
    use_when: str
    avoid_when: str
    hoomd: str
    verified: str = ""           # 이 프로젝트에서 실행 검증했는가 (근거)
    not_verified: str = ""       # ⛔ 확인되지 않은 것 · 실측으로 드러난 하드 제약
    tags: frozenset = frozenset()
    note: str = ""


# ══════════════════════════════════════════════════════════════════════
# 카탈로그 — 단일 매질 콜로이드의 표준 상호작용
# ══════════════════════════════════════════════════════════════════════
CATALOG: dict[str, Interaction] = {
    "pair.none": Interaction(
        key="pair.none", name="상호작용 없음 (독립 입자)",
        form="U_ij = 0",
        needs=(),
        coupling="—",
        use_when="관측량이 단일입자 성질(MSD·MSAD·트랩 안 분포)이고 밀도가 낮을 때. "
                 "독립 입자 앙상블로 통계만 얻는 용도.",
        avoid_when="구조(g(r)·voronoi)나 집단 거동이 목표일 때.",
        hoomd="forces 목록에서 페어 퍼텐셜을 빼면 된다.",
        verified="1-A trap-2d-5um — 관측량 4종 해석해 일치",
        tags=frozenset({"single_particle", "dilute"})),

    "pair.wca": Interaction(
        key="pair.wca", name="WCA (순수 배제부피)",
        form="U = 4ε[(σ/r)¹² − (σ/r)⁶] + ε,  r < 2^(1/6)σ",
        needs=("σ (= 입자 지름)", "ε (관례상 kT)"),
        coupling="없음 — 사실상 하드코어. 구조는 φ 하나가 정한다.",
        use_when="전하가 잘 스크리닝되어 반발이 짧고, 겹침만 막으면 될 때. "
                 "다른 상호작용의 코어로도 거의 항상 함께 쓴다.",
        avoid_when="장거리 반발/인력이 구조를 지배할 때 (그것만으로는 부족).",
        hoomd="md.pair.LJ(r_cut=2^(1/6)σ, mode='shift')  ← 전용 클래스 없음 (함정 4)",
        verified="1-B soft-r3 — 코어로 사용, 2입자 힘 대조 0.000%",
        tags=frozenset({"excluded_volume", "screened", "structure"})),

    "pair.yukawa": Interaction(
        key="pair.yukawa", name="Yukawa / 스크린된 Coulomb (DLVO 반발 분기)",
        form="U(r)/kT = ε_Y · exp(−κ(r−σ)) / (r/σ),    κ = 1/λ_D",
        needs=("유효 표면전하 또는 ζ 전위", "이온세기 → Debye 길이 λ_D",
               "매질 유전율 ε_r"),
        coupling="ε_Y = U(σ)/kT  (접촉 결합)  ·  κσ = σ/λ_D  (스크리닝 vs 입자)",
        use_when="하전 콜로이드 + 전해질. **단일 매질 수용액계의 가장 표준적 선택**이고, "
                 "λ_D 로 반발 거리를 실험적으로 조절할 수 있다.",
        avoid_when="염 농도가 매우 높아 λ_D ≪ σ 이면 WCA와 구별되지 않는다. "
                   "반대로 매우 낮으면 λ_D 가 박스를 넘어 최소이미지가 깨진다.",
        hoomd="md.pair.Yukawa(nlist, default_r_cut=...)  — params: kappa, epsilon",
        verified="",   # 이 프로젝트에서 아직 안 씀
        tags=frozenset({"charged", "screened", "repulsive", "structure", "dlvo"}),
        note="λ_D 가 새 길이척도로 원장에 들어오고, 분리 검사 `λ_D ≤ L/4` 가 추가된다."),

    "pair.soft_power": Interaction(
        key="pair.soft_power", name="소프트 멱함수 반발 r⁻ⁿ",
        form="U(r)/kT = A (σ/r)^n     (n=3 이면 2D 쌍극자계의 표준 모델)",
        needs=("진폭 A (무차원 또는 length^n)", "지수 n", "컷오프 r_c"),
        coupling="Γ = U(a_mean)/kT = A (σ/a_mean)^n   ★ A 단독이 아니라 밀도와의 조합",
        use_when="외부장으로 유도된 쌍극자(상자성/유전 콜로이드), 또는 소프트 코로나. "
                 "n=3·2D 는 계면 상자성 콜로이드의 고전 모델.",
        avoid_when="n ≤ dim 이면 에너지 적분이 발산해 컷오프로 처리할 수 없다. "
                   "n 이 작으면 꼬리가 길어 최소이미지가 먼저 걸린다.",
        hoomd="md.pair.Table (endpoint=False! 함정 10) + WCA 코어 (함정 11)",
        verified="1-B soft-r3 전체 — 2입자 0.000%, 에너지 일관성 +0.00~0.67%, "
                 "Γ 3~30 사이에서 육방 결정 전이 확인",
        tags=frozenset({"repulsive", "soft", "field_induced", "structure", "dipolar"}),
        note="2D r⁻³ 는 꼬리 에너지가 1/r_c 로만 줄어 **절대 kT 컷오프 기준이 성립하지 "
             "않는다.** 컷오프를 a_mean 배수로 잡고 수렴을 확인해야 한다 (1-B 교훈)."),

    "pair.dlvo": Interaction(
        key="pair.dlvo", name="DLVO 전체 (스크린된 반발 + van der Waals 인력)",
        form="U = U_Yukawa(r) − (A_H/12)·(σ/(r−σ)) 류 + 배제부피",
        needs=("λ_D·유효전하 (위 Yukawa와 동일)", "Hamaker 상수 A_H",
               "접촉 차단 거리 (vdW 발산 방지)"),
        coupling="ε_Y, κσ, 그리고 A_H/kT  (1차 최소 깊이)",
        use_when="응집·젤화의 시작을 보려 할 때. 2차 최소·에너지 장벽이 물리의 핵심일 때.",
        avoid_when="vdW가 접촉에서 발산하므로 반드시 코어/차단이 필요하고, "
                   "1차 최소가 깊으면 사실상 비가역 접착이 되어 페어 퍼텐셜 그림이 무너진다 "
                   "(→ contact.adhesive_bending 을 보라).",
        hoomd="md.pair.Table 로 합성 (Yukawa+vdW+코어를 하나의 표로)",
        verified="",
        tags=frozenset({"charged", "attractive", "aggregation", "dlvo"}),
        note="A_H 는 재질·매질 조합마다 다르다 (~1e-21 ~ 1e-19 J 범위). "
             "핸드북/문헌에서 가져오고 tier 를 붙여야 한다 — 지어내지 말 것."),

    "pair.ao_depletion": Interaction(
        key="pair.ao_depletion", name="Asakura–Oosawa 고갈 인력",
        form="U_dep(r) < 0, 유효 거리 ≈ 2R_g (폴리머 크기). 깊이 ∝ 폴리머 삼투압",
        needs=("폴리머 회전반경 R_g", "폴리머 농도(삼투압)"),
        coupling="깊이/kT  ·  거리비 2R_g/σ  (단거리성)",
        use_when="비흡착 폴리머를 넣어 인력을 **거리와 깊이를 따로** 조절하고 싶을 때. "
                 "콜로이드-폴리머 상분리·젤의 표준 도구.",
        avoid_when="폴리머가 흡착하면(입체 안정화) 부호가 뒤집힌다. 모델이 안 맞는다.",
        hoomd="md.pair.Table 로 구현",
        verified="",
        tags=frozenset({"attractive", "short_range", "aggregation", "tunable"})),

    "contact.adhesive_bending": Interaction(
        key="contact.adhesive_bending",
        name="접착 접촉 + 접선(굽힘) 강성  ★ 페어 퍼텐셜이 아니다",
        form="반경 방향: 접착 접촉(JKR) — 사실상 비신축 결합\n"
             "        접선 방향: 결합각 강성 κ_θ = EI/ℓ,  EI = πE a_c⁴/4\n"
             "        임계 모멘트 M_c 위에서는 미끄러짐/구름 (소성)",
        needs=("입자 탄성률 E, 포아송비 ν", "접착 에너지 W_SL (또는 접촉반경 a_c)",
               "임계 굽힘 모멘트 M_c"),
        coupling="κ_θ/kT (열적으로 뻣뻣한가) · M/M_c (선형 탄성 범위 안인가)",
        use_when="입자가 **접촉해 응집체/사슬/젤**을 이룰 때. 중심 페어 퍼텐셜만으로는 "
                 "단일 결합이 토크를 지탱하는 현상을 재현할 수 없다.",
        avoid_when="입자가 분산되어 접촉하지 않는 계. 그땐 보통의 페어 퍼텐셜을 쓴다.",
        hoomd="md.bond.Harmonic(뻣뻣, r0=σ) + md.angle.Harmonic(k=κ_θ, t0=π) "
              "(+ WCA). M_c 초과 소성은 내장에 없음 → 커스텀 또는 M<M_c 로 제한. "
              "⛔ **angle.Harmonic 은 거의 곧은 사슬에서 힘이 틀린다** — 아래 not_verified",
        verified="κ_θ = EI/ℓ 매핑을 이산 사슬 굽힘으로 검증 "
                 "(scratch/chain_bend_from_papers.py: n=51에서 빔 공식과 −0.08%; "
                 "n=25 는 −0.35%). 이산 행렬 ↔ HOOMD angle.Harmonic 도 정적으로 0.55% "
                 "일치 (scratch/verify_angle_matrix.py) — 단 **에너지 기준**이다",
        not_verified=(
            "⛔ **HOOMD 하드 제약**: md.angle.Harmonic 은 sin θ 를 SMALL = √2×10⁻³ "
            "(실측 1.414217e-3) 로 클램프한다. sin θ < SMALL 이면 힘이 sinθ/SMALL 배로 "
            "축소되어 ∝ κ(θ−π)² — **선형이 아니라 2차**가 되고, 사슬이 실제보다 무르고 "
            "비선형이 된다. **에너지는 전 구간 정확**하므로 에너지로 검증하면 통과한다. "
            "t0=π 는 평형 자체가 sin θ=0 이라 **κ_θ/kT 가 클수록(뻣뻣할수록) 심하다** — "
            "이 상호작용이 바로 그 영역이다. 설계 시 응답 진폭에서 **모든** 결합각의 "
            "min|θ−π| > SMALL 인지 확인할 것 (최대만 보면 안 된다 — 사슬 끝 각도가 중앙보다 "
            "한 자릿수 작다). chain-bend-2d-oscill 은 23개 각도 전부가 깨진 영역이라 "
            "실행이 거부됐다. 우회: force.Custom 직접 구현(26배 느림) · κ_θ 를 낮춤 · "
            "선형 영역이면 해석적으로 풀기. angle.Table 과 CosineSquared 는 둘 다 대체 불가. "
            "재현: scratch/verify_angle_force_small_theta.py (skill bd-hoomd 함정 15)"),
        tags=frozenset({"bonded", "aggregate", "contact", "gel", "tangential"}),
        note="Pantina & Furst (PRL 94 138301; Langmuir 24 1141) 의 핵심 결과. "
             "단일 결합이 토크를 지탱한다는 것이 실측이고, DLVO 같은 중심력 모델은 "
             "이걸 못 낸다 — 논문이 명시적으로 그 점을 지적한다."),
}


# ══════════════════════════════════════════════════════════════════════
# 추천기
# ══════════════════════════════════════════════════════════════════════
INTERACTION_SYMBOLS = frozenset({
    "U_ij", "pair_potential", "interaction", "interactions", "pair", "U", "potential",
})


def looks_like_interaction(symbol: str) -> bool:
    s = (symbol or "").strip()
    return s in INTERACTION_SYMBOLS or "potential" in s.lower() or s.lower().startswith("u_")


def infer_tags(obs) -> set:
    """observation.yaml 에서 문맥 태그를 뽑는다. 추측이 아니라 적혀 있는 것만 본다."""
    raw = obs.raw if hasattr(obs, "raw") else dict(obs)
    text = " ".join([
        str(raw.get("system_guess", "")),
        str(raw.get("raw_transcription", "")),
        " ".join(str(e) for e in (raw.get("entities") or [])),
        " ".join(str(g) for g in (raw.get("stated_goals") or [])),
    ]).lower()
    tags: set[str] = set()
    kinds = {str((e or {}).get("kind", "")).lower() for e in (raw.get("entities") or [])
             if isinstance(e, dict)}

    if "bonded_interaction" in kinds or any(w in text for w in
                                           ("chain", "사슬", "aggregate", "응집", "비드", "bead")):
        tags |= {"bonded", "aggregate", "contact"}
    if any(w in text for w in ("rdf", "voronoi", "structure", "구조", "결정", "configuration")):
        tags.add("structure")
    if any(w in text for w in ("msd", "msad", "확산", "diffusion")):
        tags.add("single_particle")
    if any(w in text for w in ("charge", "전하", "zeta", "salt", "mgcl", "nacl", "전해질")):
        tags.add("charged")
    if any(w in text for w in ("magnetic", "자성", "dipole", "쌍극자", "field", "외부장")):
        tags.add("field_induced")
    # ★ '레올로지'라는 단어만으로 접촉/젤 모델을 추천하면 안 된다.
    #   trap-drag 는 "능동 마이크로레올로지"라 적혀 있지만 끌린 탐침이고 접촉계가 아니다
    #   (추천기를 실제 케이스에 돌려보고 잡은 오탐).
    #   탄성률 측정 + **접촉/응집 증거**가 함께 있을 때만 접선 상호작용을 후보로 올린다.
    rheo = any(w in text for w in ("g'", "g''", "modulus", "탄성률", "rheolog", "레올로지"))
    if rheo and ({"bonded", "aggregate", "contact"} & tags):
        tags |= {"gel", "tangential"}
    elif rheo:
        tags.add("rheology")
    if any(w in text for w in ("hexagonal", "육방", "lattice", "격자")):
        tags |= {"structure", "repulsive"}
    if any(w in text for w in ("trap", "트랩", "tweezer", "집게")):
        tags.add("trapped")
    return tags


def recommend(obs, top: int = 3) -> list:
    """(Interaction, 점수, 근거) 목록. 점수는 태그 일치 수 + 프로젝트 검증 보너스."""
    tags = infer_tags(obs)
    scored = []
    for it in CATALOG.values():
        hit = tags & it.tags
        if not hit:
            continue
        score = len(hit) + (0.5 if it.verified else 0.0)
        why = f"문맥 태그 일치: {', '.join(sorted(hit))}"
        if it.verified:
            why += f" · 이 프로젝트에서 검증됨"
        scored.append((it, score, why))
    scored.sort(key=lambda x: -x[1])
    return scored[:top], tags


def yaml_snippet(it: Interaction, symbol: str = "U_ij") -> str:
    """`observation.yaml` 의 resolution 에 붙여넣을 초안."""
    needs = "\n".join(f"    #   {n}" for n in it.needs) or "    #   (추가 파라미터 없음)"
    return (f'  - symbol: {symbol}\n'
            f'    resolution: |\n'
            f'      {it.key} — {it.name}\n'
            f'      {it.form.splitlines()[0]}\n'
            f'      결합세기: {it.coupling}\n'
            f'      ★ 아래 값을 채워야 물리계가 확정된다 (각각 출처+tier 필요):\n'
            f'{needs}\n'
            f'    confirmed_by: user        # ← 확인하면 채우세요\n')


def render_suggestion(obs, case_name: str = "") -> str:
    """추천 + 질문. 결정은 사람이 한다."""
    recs, tags = recommend(obs)
    open_int = [m for m in obs.open_missing
                if looks_like_interaction(str(m.get("symbol", "")))]

    L: list[str] = []
    w = L.append
    w("=" * 80)
    w(f"intake suggest — {case_name or obs.path.parent.name}")
    w("=" * 80)
    w("단일 매질 콜로이드의 표준 골격은 DLVO 입니다:")
    w("  U = U_전기이중층(스크린된 Coulomb) + U_van der Waals  (+ 배제부피)")
    w("계에 따라 고갈 인력 · 장 유도 쌍극자 · 접촉 접선력이 붙습니다.")
    w("")
    w(f"이 케이스에서 읽어낸 문맥 태그: {', '.join(sorted(tags)) or '(없음)'}")

    if not open_int:
        w("")
        w("상호작용이 미해소 결측으로 잡혀 있지 않습니다.")
        if obs.raw.get("missing_required"):
            w("(이미 해소되었거나, 애초에 상호작용이 없는 계입니다)")
    else:
        w("")
        w(f"미지정 상호작용 {len(open_int)}건: "
          + ", ".join(str(m.get("symbol")) for m in open_int))

    w("")
    w("─" * 80)
    w("추천 (근거 순)")
    w("─" * 80)
    for i, (it, score, why) in enumerate(recs, 1):
        mark = "★ 1순위" if i == 1 else f"  {i}순위"
        w("")
        w(f"{mark}  {it.key}  —  {it.name}")
        w(f"    형태     {it.form.splitlines()[0]}")
        for extra in it.form.splitlines()[1:]:
            w(f"             {extra.strip()}")
        w(f"    결합세기 {it.coupling}")
        w(f"    쓸 때     {it.use_when}")
        w(f"    피할 때   {it.avoid_when}")
        w(f"    HOOMD    {it.hoomd}")
        w(f"    검증     {it.verified or '이 프로젝트에서 아직 안 씀 (미검증)'}")
        if it.needs:
            w(f"    ★ 채워야 하는 값: {', '.join(it.needs)}")
        if it.note:
            for ln in _wrap(it.note, 72):
                w(f"    비고     {ln}" if ln == _wrap(it.note, 72)[0] else f"             {ln}")
        w(f"    근거     {why}")

    w("")
    w("=" * 80)
    w("질문 — 어떻게 하시겠습니까?")
    w("=" * 80)
    if recs:
        top = recs[0][0]
        w(f"  (a) 1순위 `{top.key}` 를 쓴다")
        w(f"      → 위의 '채워야 하는 값'을 주시면 system.yaml 을 씁니다.")
        w(f"  (b) 다른 후보를 쓴다 — 위 목록에서 고르거나 카탈로그에 없는 것을 직접 지정")
        w(f"  (c) 상호작용 없음 (`pair.none`) — 단일입자 관측량만 볼 때")
        w("")
        w("  붙여넣기용 초안 (observation.yaml 의 missing_required 항목):")
        w("")
        for ln in yaml_snippet(top, str(open_int[0].get("symbol")) if open_int else "U_ij").splitlines():
            w(f"  {ln}")
    else:
        w("  문맥 태그로 후보를 좁히지 못했습니다. 카탈로그 전체를 보세요:")
        w("    $PY -m bdbot.cli interactions list")
    w("")
    w("  ⚠️ 어느 쪽이든 **구체적 수치는 지어내지 않습니다.** 핸드북/문헌 값을 주시거나,")
    w("     제가 제안하면 tier 3(임의 가정)으로 표시해 system.yaml 에 남깁니다.")
    w("=" * 80)
    return "\n".join(L)


def render_catalog() -> str:
    L = ["=" * 80, "콜로이드 상호작용 카탈로그 (단일 매질)", "=" * 80,
         "DLVO = 스크린된 Coulomb + van der Waals 가 표준 골격.", ""]
    for it in CATALOG.values():
        v = "검증됨" if it.verified else "미검증"
        if it.not_verified:
            v += " · ⛔제약있음"
        L.append(f"  {it.key:<28} {it.name}")
        L.append(f"  {'':<28} {it.form.splitlines()[0]}")
        L.append(f"  {'':<28} 결합세기: {it.coupling}   [{v}]")
        if it.not_verified:                    # 하드 제약은 접지 않고 그대로 보여준다
            for ln in _wrap(it.not_verified, 74):
                L.append(f"  {'':<28} {ln}")
        L.append("")
    L.append("상세: $PY -m bdbot.cli intake suggest <case>")
    L.append("=" * 80)
    return "\n".join(L)


def _wrap(s: str, width: int) -> list:
    words, lines, cur = s.split(), [], ""
    for wd in words:
        if len(cur) + len(wd) + 1 > width:
            lines.append(cur)
            cur = wd
        else:
            cur = f"{cur} {wd}".strip()
    if cur:
        lines.append(cur)
    return lines or [""]


__all__ = ["SCHEMA", "Interaction", "CATALOG", "recommend", "infer_tags",
           "render_suggestion", "render_catalog", "yaml_snippet",
           "looks_like_interaction", "INTERACTION_SYMBOLS"]
