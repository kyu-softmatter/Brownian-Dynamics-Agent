"""L3 검사기를 **일부러 망가뜨려** 시험합니다 (CLAUDE.md 작업 관행).

"조용히 통과"와 "검사를 안 함"은 다릅니다. `verify_intake_guards.py` 가 이 방식으로
실제 크래시 버그 1건을 잡았습니다. 여기서 시험하는 것:

  ① 정상 스펙은 통과하는가 (거짓 양성이 없는가)
  ② 원장에서 필수 역할을 빼면 잡는가
  ③ 무차원수 값을 원장과 어긋나게 바꾸면 잡는가  ← 예전엔 잡을 수 없었던 종류
  ④ 무차원수가 없는 기호를 가리키면 잡는가 (크래시가 아니라 오류로)
  ⑤ 물리계를 바꾸면 run_id 가 바뀌는가          ← 결함 ① 회귀 방지
  ⑥ 역변환이 왕복하는가 (무차원 → 물리 → 무차원)
  ⑦ 스펙을 저장→로드하면 run_id·무차원수가 보존되는가 (L4가 이 경로만 씁니다)
  ⑧ 손으로 고친 스펙을 `verify_hash()` 가 잡는가
  ⑨ 기준 스케일이 원장에 없으면 잡는가
  ⑩ dt* 가 0/음수면 잡는가
  ⑪ 단위가 안 맞는 비(길이/시간)를 무차원수라고 하면 잡는가

    $PY scratch/verify_nondim_guards.py
"""
from __future__ import annotations

import copy
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import Q, nondim as ND, scales as SC  # noqa: E402
from bdbot.checks import Check  # noqa: E402

PASS, FAIL = [], []


def check(name: str, ok: bool, detail: str = "") -> None:
    (PASS if ok else FAIL).append(name)
    print(f"  {'✓' if ok else '✗'} {name}" + (f"\n        {detail}" if detail else ""))


def errs(spec) -> list[str]:
    return [f"{i.where}: {i.msg}" for i in spec.validate() if i.level == "error"]


# ── 정상 스펙 하나를 손으로 조립 (케이스 스크립트에 의존하지 않습니다) ────────
def good_spec() -> ND.NondimSpec:
    d = Q(5.0, "um").to("m")
    kT = Q(4.141947e-21, "J")
    tau_B = Q(242.1, "s")
    tau_k = Q(4.01, "ms").to("s")
    dt = Q(8.02, "us").to("s")

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "입자 지름 (기준)")
    lg.add_length("L", Q(160.0, "um").to("m"), "박스", role="box")
    lg.add_time("tau_p", Q(3.264, "us").to("s"), "관성", role="inertia")
    lg.add_time("dt", dt, "적분 스텝", role="dt")
    lg.add_time("tau_k", tau_k, "트랩 이완", star=True)
    lg.add_time("tau_B", tau_B, "확산 (기준)")
    lg.add_time("T_obs", Q(8.02, "s"), "관측창", role="observation")
    lg.add_energy("kT", kT, "열에너지 (기준)")
    lg.ref = SC.thermal_reference(d, kT, tau_B)
    lg.rationale = lg.ref["rationale"]

    groups = [
        ND.Group("dt/tau_k", float((dt / tau_k).to("")), ("times", "dt"), ("times", "tau_k"),
                 "", "적분 해상"),
        ND.Group("L/d", float((lg.get("lengths", "L") / d).to("")),
                 ("lengths", "L"), ("lengths", "d"), "", "박스 크기"),
        ND.Group("phi", 0.35, None, None, "", "밀집도 (비가 아닌 입력)"),
    ]
    checks = [Check("적분", "트랩 해상 dt/τ_k", float((dt / tau_k).to("")), 1e-2, "<=")]
    system = {"label": "guard-test", "dimensions": 2,
              "particle": {"diameter": {"value": 5.0, "unit": "um", "source": "test", "tier": 0}},
              "medium": {"temperature": {"value": 300, "unit": "K", "source": "t", "tier": 0}}}
    # ★ dt* 는 원장에서 유도합니다 — 손으로 3.3e-8 이라고 적었더니 원장의 dt/τ_B
    #   (3.3127e-8)와 어긋나 검사 ⑩b 가 잡았습니다. 검사기가 자기 픽스처를 잡은 셈입니다.
    dt_star = float((dt / tau_B).to("dimensionless").magnitude)
    return ND.NondimSpec(
        case="guard-test", system=system, reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"N": 1000}, numerics={"dt_star": dt_star, "n_prod": 1000000})


def main() -> int:
    print("=" * 80)
    print("L3 (bdbot.nondim) 적대적 검사")
    print("=" * 80)

    # ① 거짓 양성 없음
    print("\n[①] 정상 스펙")
    s = good_spec()
    e = errs(s)
    check("정상 스펙은 오류 0건", not e, "; ".join(e))

    # ② 필수 역할 누락
    print("\n[②] 원장에서 필수 역할 제거")
    for role, sym, cat in (("dt", "dt", "times"), ("observation", "T_obs", "times"),
                           ("box", "L", "lengths"), ("inertia", "tau_p", "times")):
        s = good_spec()
        # 원장에서 빼면 그 기호를 쓰던 무차원수도 같이 빠져야 공정한 시험이 됩니다
        getattr(s.ledger, cat).pop(sym)
        s.groups = [g for g in s.groups
                    if sym not in ((g.num or ("", ""))[1], (g.den or ("", ""))[1])]
        e = errs(s)
        check(f"'{role}' 역할 누락을 잡는다",
              any(f"ledger.{role}" in x for x in e), "; ".join(e) or "오류 없음 — 못 잡았다")

    print("\n[②b] 없는 역할을 이유와 함께 비우면 통과")
    s = good_spec()
    s.ledger.times.pop("tau_p")
    s.ledger.declare_absent("inertia", "테스트: 관성이 없는 계라고 선언")
    check("declare_absent 로 비우면 통과", not errs(s), "; ".join(errs(s)))
    try:
        good_spec().ledger.declare_absent("box", "")
        check("이유 없이 비우는 것은 거부된다", False, "빈 이유가 통과했다")
    except ValueError:
        check("이유 없이 비우는 것은 거부된다", True)

    # ③ ⭐️ 무차원수를 원장과 어긋나게
    print("\n[③] 무차원수 값을 원장과 어긋나게 (예전엔 잡을 수 없던 종류)")
    for factor, label in ((1.41, "41% 어긋남 — a_mean vs a_NN 급 실수"),
                          (1 + 1e-6, "1e-6 어긋남 — 미세한 계산 실수")):
        s = good_spec()
        s.groups[0].value *= factor
        e = errs(s)
        check(f"dt/tau_k 를 ×{factor} 하면 잡는다 ({label})",
              any("groups.dt/tau_k" in x for x in e), "; ".join(e) or "오류 없음 — 못 잡았다")

    print("\n[③b] 부동소수 오차 수준(1e-12)은 통과해야 한다 (거짓 양성 방지)")
    s = good_spec()
    s.groups[0].value *= (1 + 1e-12)
    check("1e-12 어긋남은 통과", not errs(s), "; ".join(errs(s)))

    # ④ 없는 기호를 가리킴 → 크래시가 아니라 오류
    print("\n[④] 무차원수가 원장에 없는 기호를 가리킴")
    s = good_spec()
    s.groups[0].den = ("times", "tau_없음")
    try:
        e = errs(s)
        check("없는 기호를 오류로 보고한다 (크래시 아님)",
              any("groups.dt/tau_k" in x for x in e), "; ".join(e))
    except Exception as ex:
        check("없는 기호를 오류로 보고한다 (크래시 아님)", False,
              f"예외로 죽었다: {type(ex).__name__}: {ex}")

    # ⑤ ⭐️ 결함 ① 회귀 — 물리계가 run_id 에 들어가는가
    print("\n[⑤] 물리계를 바꾸면 run_id 가 바뀌는가 (결함 ① 회귀 방지)")
    a = good_spec()
    b = good_spec()
    b.system = copy.deepcopy(b.system)
    b.system["particle"]["diameter"]["value"] = 0.5
    check("d 를 바꾸면 run_id 가 바뀐다", a.run_id() != b.run_id(),
          f"{a.run_id()} vs {b.run_id()}")

    c = good_spec()
    c.system = copy.deepcopy(c.system)
    c.system["particle"]["diameter"]["source"] = "출처만 고침"
    check("출처(문서 필드)만 고치면 run_id 가 유지된다", a.run_id() == c.run_id(),
          f"{a.run_id()} vs {c.run_id()}")

    dd = good_spec()
    dd.ledger.rationale = "근거 문구만 고침"
    dd.checks = []
    check("원장 근거·검사를 고치면 run_id 가 유지된다", a.run_id() == dd.run_id(),
          f"{a.run_id()} vs {dd.run_id()}")

    # ⑥ 역변환 왕복
    print("\n[⑥] 역변환 왕복 (무차원 → 물리 → 무차원)")
    s = good_spec()
    for val, kw, unit in ((0.5, dict(L=2), "um^2"), (3.0, dict(T=1), "s"),
                          (1.83, dict(L=2, T=-1), "um^2/s"), (2.0, dict(E=1), "J")):
        phys = s.physical(val, **kw)
        back = phys
        for kind, expo in (("length", kw.get("L", 0)), ("time", kw.get("T", 0)),
                           ("energy", kw.get("E", 0))):
            back = back / s.reference.si(kind) ** expo
        rel = abs(float(back.to("dimensionless").magnitude) - val) / val
        check(f"{kw} 왕복 오차 < 1e-12  ({phys.to(unit):~.5gP})", rel < 1e-12, f"rel={rel:.2e}")

    # ⑦ 저장 → 로드 보존
    print("\n[⑦] 저장 → 로드 (L4가 쓰는 경로)")
    s = good_spec()
    tmp = ROOT / "scratch" / "_tmp_spec.json"
    s.write(tmp)
    ls = ND.load(tmp)
    check("run_id 보존", ls.run_id == s.run_id(), f"{ls.run_id} vs {s.run_id()}")
    check("무차원수 보존", abs(ls.group("dt/tau_k") - s.groups[0].value) < 1e-15)
    check("환산값(reduced) 계산됨 — L4가 L* 를 여기서 읽는다",
          abs(ls.reduced("lengths", "L") - 32.0) < 1e-9,
          f"L/d = {ls.reduced('lengths', 'L')}")
    r1 = float(s.physical(1.83, L=2, T=-1).to("um^2/s").magnitude)
    r2 = float(ls.physical(1.83, L=2, T=-1).to("um^2/s").magnitude)
    check("로드한 스펙의 역변환이 원본과 같다", abs(r1 - r2) / r1 < 1e-12, f"{r1} vs {r2}")
    ok, want = ls.verify_hash()
    check("해시 자기검증 통과", ok, f"기대 {want}")

    # ⑧ 손으로 고친 스펙
    print("\n[⑧] 스펙을 손으로 고침 (§16 규칙 2 — 스펙을 손으로 쓰지 않는다)")
    raw = json.loads(tmp.read_text())
    raw["params"]["N"] = 4000
    tmp.write_text(json.dumps(raw))
    ok, want = ND.load(tmp).verify_hash()
    check("params 를 손으로 고치면 해시 불일치를 잡는다", not ok, f"기대 {want}")

    raw = json.loads(tmp.read_text())
    raw["params"]["N"] = 1000
    raw["ledger_absent"] = {"note": "문서 필드만 고침"}
    tmp.write_text(json.dumps(raw))
    ok, _ = ND.load(tmp).verify_hash()
    check("해시 대상이 아닌 필드를 고치면 통과한다", ok)
    tmp.unlink()

    # ⑨ 기준이 원장에 없음
    print("\n[⑨] 기준 스케일이 원장에 없음")
    s = good_spec()
    s.ledger.times.pop("tau_B")
    s.groups = [g for g in s.groups if "tau_B" not in str(g.num) + str(g.den)]
    e = errs(s)
    check("기준 시간이 원장에 없으면 잡는다", any("reference.time" in x for x in e),
          "; ".join(e) or "오류 없음")

    # ⑩ dt* 이상값
    print("\n[⑩] 실행 파라미터 이상값")
    for bad, why in ((0.0, "0"), (-1e-8, "음수")):
        s = good_spec()
        s.numerics["dt_star"] = bad
        e = errs(s)
        check(f"dt* = {why} 를 잡는다", any("numerics.dt_star" in x for x in e), "; ".join(e))
    s = good_spec()
    s.numerics.pop("n_prod")
    check("n_prod 누락을 잡는다", any("numerics.n_prod" in x for x in errs(s)))

    print("\n[⑩b] dt* 가 원장의 dt/τ_B 와 어긋남 (HOOMD가 다른 스텝으로 돈다)")
    s = good_spec()
    true_dt_star = float((s.ledger.get("times", "dt")
                          / s.reference.si("time")).to("dimensionless").magnitude)
    check("원장에서 유도한 dt* 는 통과", not errs(s),
          f"원장 dt/τ_B = {true_dt_star:.8e}")
    for factor, why in ((1.01, "1%"), (1 + 1e-6, "1e-6")):
        s2 = good_spec()
        s2.numerics["dt_star"] = true_dt_star * factor
        e = errs(s2)
        check(f"{why} 어긋나면 잡는다", any("numerics.dt_star" in x for x in e), "; ".join(e))

    # ⑪ 단위가 안 맞는 비
    print("\n[⑪] 무차원이 아닌 비 (길이/시간)")
    s = good_spec()
    s.groups[0].den = ("lengths", "d")          # dt / d → 무차원 아님
    e = errs(s)
    check("차원이 안 맞는 비를 잡는다", any("groups.dt/tau_k" in x for x in e),
          "; ".join(e) or "오류 없음 — 못 잡았다")

    print("\n" + "=" * 80)
    print(f"통과 {len(PASS)} · 실패 {len(FAIL)}")
    for f in FAIL:
        print(f"   ✗ {f}")
    print("✓ PASS — L3 검사기가 실제로 잡는다" if not FAIL else "✗ FAIL — 위 항목이 안 잡힌다")
    print("=" * 80)
    return 0 if not FAIL else 1


if __name__ == "__main__":
    sys.exit(main())
