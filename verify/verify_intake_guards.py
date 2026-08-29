"""앞단 검사기 적대적 테스트 — 일부러 망가뜨린 스펙을 정말 잡는가.

"조용히 통과한 것"과 "검사를 안 한 것"은 다릅니다. 이 스크립트는 L0(observation)·
L2(system) 검사기에 고장난 입력을 넣어 **각 규칙이 실제로 발동하는지** 확인합니다.

이 테스트로 실제 버그를 하나 잡았습니다: 단위를 `furlong^2`로 바꿨을 때 검사기가
오류를 보고하는 대신 pint `DimensionalityError`로 크래시했습니다 (bdbot.physical.bulk).

    $PY scratch/verify_intake_guards.py
"""
import copy
import pathlib
import sys
import tempfile

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import intake as I  # noqa: E402
from bdbot import interactions as X  # noqa: E402
from bdbot import physical as P  # noqa: E402
CASE_OK = ROOT / "intake/soft-r3-2d-A-sweep"

results = []


def report(ok, label, detail=""):
    results.append(ok)
    print(f"  {'✓' if ok else '✗'} {label:<44}{detail}")


# ══════════════════════════════════════════════════════════════════════
print("=" * 78)
print("① L0 observation 검사기")
print("=" * 78)
obs_base = yaml.safe_load((CASE_OK / "observation.yaml").read_text())
tmp = pathlib.Path(tempfile.mkdtemp())
(tmp / "sketch_01.jpeg").write_bytes(b"")


def obs_check(mutate, label, want_error=True):
    d = copy.deepcopy(obs_base)
    mutate(d)
    (tmp / "observation.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    try:
        o = I.load(tmp)
    except Exception as e:
        report(False, label, f"크래시! {type(e).__name__}")
        return
    errs = o.errors
    ok = (len(errs) > 0) == want_error
    report(ok, label, f"오류 {len(errs)}건" + (f"  → {errs[0].msg[:44]}" if errs else ""))


obs_check(lambda d: None, "원본 (오류 0이 정답)", want_error=False)
obs_check(lambda d: d.pop("ambiguities"), "ambiguities 키 삭제 (§8.3)")
obs_check(lambda d: d.pop("unread_regions"), "unread_regions 키 삭제 (§8.3)")
obs_check(lambda d: d.pop("raw_transcription"), "전사 삭제 (규칙 5)")
obs_check(lambda d: d.__setitem__("raw_transcription", "짧음"), "전사가 너무 짧음")
obs_check(lambda d: d["ambiguities"][0].pop("resolution"),
          "ambiguity 의 resolution 키 삭제")
obs_check(lambda d: d["missing_required"][0].pop("confidence"),
          "가정값에서 tier 삭제 (규칙 3) ★")
obs_check(lambda d: d["stated_quantities"][0].__setitem__("source", ""),
          "명시 수치의 source 를 비움 (원칙 2)")
obs_check(lambda d: d["ambiguities"].__setitem__(1, dict(d["ambiguities"][1], id="A1")),
          "ambiguity id 중복")
obs_check(lambda d: d["missing_required"][2].__setitem__("kind", "maybe"),
          "kind 를 허용값 밖으로")

# choice/physical 구분이 판정을 바꾸는가
d = copy.deepcopy(obs_base)
(tmp / "observation.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
o = I.load(tmp)
report(not o.open_missing and len(o.open_choices) == 1,
       "kind:choice 는 차단하지 않음",
       f"차단 {len(o.open_missing)} · 미정선택 {len(o.open_choices)}")
for m in d["missing_required"]:
    if m.get("kind") == "choice":
        m.pop("kind")
(tmp / "observation.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
o2 = I.load(tmp)
report(len(o2.open_missing) == 1, "kind 없으면 physical 로 보수적 판정",
       f"차단 {len(o2.open_missing)}건")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("② L2 system 검사기")
print("=" * 78)
sys_base = yaml.safe_load((CASE_OK / "system.yaml").read_text())
t2 = pathlib.Path(tempfile.mkdtemp())
(t2 / "observation.yaml").write_text((CASE_OK / "observation.yaml").read_text())
(t2 / "sketch_01.jpeg").write_bytes(b"")


def sys_check(mutate, label, want_error=True):
    d = copy.deepcopy(sys_base)
    mutate(d)
    (t2 / "system.yaml").write_text(yaml.safe_dump(d, allow_unicode=True))
    try:
        s = P.load(t2)
    except Exception as e:
        report(False, label, f"크래시! {type(e).__name__}")
        return
    errs = s.errors
    ok = (len(errs) > 0) == want_error
    report(ok, label, f"오류 {len(errs)}건" + (f"  → {errs[0].msg[:44]}" if errs else ""))


sys_check(lambda d: None, "원본 (오류 0이 정답)", want_error=False)
sys_check(lambda d: d["derived_scales"].__setitem__("tau_B", {"value": 300.0, "unit": "s"}),
          "τ_B 를 242→300 (재계산 불일치) ★")
sys_check(lambda d: d["particle"]["diameter"].pop("source"), "d 의 source 삭제 (원칙 2)")
sys_check(lambda d: d["particle"]["diameter"].__setitem__("unit", "furlong^2"),
          "단위를 차원 불일치로 (예전엔 크래시) ★")
sys_check(lambda d: d["particle"]["diameter"].__setitem__("unit", "nonsense_unit"),
          "단위를 해석 불가로")
sys_check(lambda d: d["particle"]["diameter"].__setitem__("tier", 7), "tier 를 7 로")
sys_check(lambda d: d.pop("derived_from"), "derived_from 삭제 (§5.4 불변식) ★")
sys_check(lambda d: d.__setitem__("derived_from", "nope.yaml"), "없는 파일 참조")
sys_check(lambda d: d.pop("medium"), "medium 섹션 삭제")

# ★ L0이 BLOCKED인데 L2가 확정돼 있으면
# 실제 케이스 상태에 의존하지 않게 **합성 입력**으로 만든다 — 케이스가 해소되면
# 테스트가 조용히 무력화되기 때문이다 (abp-rod 가 해소되자 실제로 그랬다).
obs_blocked = copy.deepcopy(obs_base)
obs_blocked["missing_required"].append({
    "symbol": "made_up_param", "kind": "physical",
    "what": "합성 테스트용 미해소 물리 결측", "assumed_value": None, "resolution": None})
(t2 / "observation.yaml").write_text(yaml.safe_dump(obs_blocked, allow_unicode=True))
(t2 / "system.yaml").write_text(yaml.safe_dump(sys_base, allow_unicode=True))
s = P.load(t2)
blk = [i for i in s.errors if "미해소 물리 결측" in i.msg]
report(bool(blk), "L0 BLOCKED인데 L2 확정 (규칙 3) ★",
       blk[0].msg[:44] if blk else "못 잡음!")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("③ 실제 5개 케이스 판정")
print("=" * 78)
# 2026-08-04: 막혀 있던 3케이스가 해소됐다 — 사용자 확정(abp-rod 형상·텀블) ·
# 논문 증류(chain-bend U_ij) · 사용자 확정 + ★제안(trap-drag 페어·밀도).
# tier 3 제안이 섞여 있으므로 READY 는 "L3로 넘어갈 수 있다"는 뜻이고 승인 완료가 아니다.
expect = {"trap-2d-5um": "READY", "soft-r3-2d-A-sweep": "READY",
          "abp-rod-2d-run-flip": "READY", "chain-bend-2d-oscill": "READY",
          "trap-drag-2d-hex300": "READY"}
for name, want in expect.items():
    o = I.load(ROOT / "intake" / name)
    got = "FAIL" if o.errors else ("BLOCKED" if o.open_missing else "READY")
    report(got == want, f"{name} → {got}", f"(기대 {want})")

# ══════════════════════════════════════════════════════════════════════
print()
print("=" * 78)
print("④ 상호작용 추천기 — 미지정 U_ij 에 무엇을 권하는가")
print("=" * 78)
EXPECT_TOP = {
    # 케이스 → (1순위 키, 왜 그래야 하는가)
    "chain-bend-2d-oscill": ("contact.adhesive_bending",
                             "비드 사슬 + 탄성률 측정 → 접촉 접선력 (Furst 논문)"),
    "trap-drag-2d-hex300": ("pair.soft_power",
                            "육방 격자 + 구조 → 소프트 반발 (사용자 확정과 일치)"),
    "abp-rod-2d-run-flip": ("pair.none",
                            "MSD·MSAD 단일입자 관측량 → 상호작용 불필요"),
}
for name, (want, why) in EXPECT_TOP.items():
    o = I.load(ROOT / "intake" / name)
    recs, tags = X.recommend(o)
    got = recs[0][0].key if recs else "(없음)"
    report(got == want, f"{name[:26]} → {got}", f"({why})")

# 오탐 회귀: '레올로지'라는 단어만으로 접촉 모델을 올리지 않는다
o = I.load(ROOT / "intake/trap-drag-2d-hex300")
tags = X.infer_tags(o)
report("tangential" not in tags and "gel" not in tags,
       "'마이크로레올로지' 단어만으로 접촉 태그 안 붙음",
       f"(태그: {', '.join(sorted(tags))})")

# 카탈로그 무결성
bad = [k for k, it in X.CATALOG.items()
       if not it.form or not it.use_when or not it.avoid_when or not it.hoomd]
report(not bad, "카탈로그 항목이 전부 형태·용도·회피·HOOMD 매핑을 가짐",
       f"(빈 항목: {bad})" if bad else f"({len(X.CATALOG)}개)")
needs_ok = [k for k, it in X.CATALOG.items() if it.key != "pair.none" and not it.needs]
report(not needs_ok, "pair.none 을 뺀 모든 항목이 '채워야 하는 값'을 명시",
       f"(누락: {needs_ok})" if needs_ok else "")
# 검증 여부가 정직하게 표시되는가 (안 쓴 것을 검증됨으로 적지 않았는가)
unused_but_verified = [k for k in ("pair.yukawa", "pair.dlvo", "pair.ao_depletion")
                       if X.CATALOG[k].verified]
report(not unused_but_verified, "안 써본 상호작용을 '검증됨'으로 적지 않음",
       f"(거짓 주장: {unused_but_verified})" if unused_but_verified else "")

print()
print("=" * 78)
print(f"{'✓ PASS' if all(results) else '✗ FAIL'} — {sum(results)}/{len(results)} 정상")
print("=" * 78)
raise SystemExit(0 if all(results) else 1)
