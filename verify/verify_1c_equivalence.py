"""Phase 1-C 동등성 회귀 테스트 — 리팩터가 결과를 바꾸지 않았는가?

Phase 1-C의 DoD: **1-A와 1-B가 같은 코드 경로로 돌고, 결과가 이전과 동일.**
`scratch/ref_1b/*.metrics.json` 이 리팩터 **전**(Phase 1-B) 결과의 스냅샷입니다.

세 가지를 봅니다:
  ① run_id  — 스펙 해시. 같으면 원장→검사→dt 경로가 바이트 단위로 동일하다는 뜻
             ⚠️ 1-A의 run_id는 **앞단 작업에서 한 번 의도적으로 바뀌었습니다** (아래 EXPECTED_RUN_ID).
                원인: 스펙에 system.yaml 전체를 넣어 해시하고 있어서 `derived_from` 필드를
                추가한 것만으로 무효화됐습니다. 물리 필드만 해시하도록 고쳤습니다
                (bdbot.runid.physics_only). 물리 수치는 96필드 전부 동일함을 확인했습니다.
  ② 물리 수치 — 관측량·무차원수·검사값이 스냅샷과 일치하는가
  ③ 스키마 — 1-C에서 추가된 필드가 무엇인지 명시 (삭제는 허용하지 않음)

    $PY scratch/verify_1c_equivalence.py
"""
from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
PY = sys.executable
REF = ROOT / "verify/ref_1b"

TARGETS = [
    ("1-A trap-2d-5um", "cases/trap_2d_5um.py", [], "trap-2d-5um.metrics.json"),
    ("1-B soft-r3 A=100", "cases/soft_r3_2d.py", ["--A", "100"], "soft-r3-A100.metrics.json"),
]

# run_id 가 스냅샷과 다른 것이 **정상인** 케이스 — 이유를 적어야 통과시킵니다.
EXPECTED_RUN_ID = {
    "trap-2d-5um__70b9394e7310": (
        "trap-2d-5um__a5ef4f45d589",
        "① 앞단 #2에서 해시 기준을 'YAML 전체' → '물리 필드만' 으로 고쳤다 "
        "(→ a49f2508556b). `derived_from` 필드를 추가한 것만으로 런이 무효화되는 것은 "
        "콘텐츠 주소의 목적에 반한다. "
        "② L3 작업에서 스펙 스키마를 세 케이스 공통(`bdbot.nondim.NondimSpec`)으로 "
        "통일했다 (→ d724e8d507cc). 해시 대상이 평평한 dict 에서 "
        "{system, params, numerics} 로 바뀌었다 — **같은 물리, 다른 배치**다. "
        "1-A 자체에는 결함이 없었고 스키마 통일의 부수효과다. "
        "③ L5 이관(2026-08-05)으로 `numerics.seed`가 스펙에 추가됐다(→ a5ef4f45d589) — "
        "`RUN.builder`의 `build(spec)`은 케이스 YAML을 다시 읽지 않으므로 이전에 "
        "main() 안에만 있던 seed=20260803 을 스펙에 넣어야 build()가 읽을 수 있다. "
        "재실행 결과가 관측량 5종 소수 15자리까지 옛 런과 일치함을 확인했다."),
    "soft-r3-2d-A-sweep__A100__27f70deab9": (
        "soft-r3-2d-A-sweep__A100__30caa5c9e0",
        "⭐️ **결함 수정**. 1-B 스펙에는 물리계가 **없었다** — d 5µm→0.5µm, η 62배, "
        "ρ_p 를 바꿔도 run_id 가 그대로였다(τ_B 16.1배 차이인데도). 완료된 다른 계의 "
        "런으로 오인돼 예전 결과가 새 계의 결과로 보고된다. "
        "`system`(physics_only)을 해시에 넣어 고쳤다 "
        "(재현: `scratch/verify_l3_spec_gaps.py`). "
        "기존 런 7개는 레거시로 보존한다 (사용자 결정 2026-08-04) — 그 런들의 "
        "물리·결과는 유효하고, 디렉토리 이름만 새 규약으로 재현되지 않는다. "
        "L5 이관(2026-08-05)은 params/numerics 키를 바꾸지 않아 이 run_id 는 그대로다 — "
        "`30caa5c9e0` 를 `RUN.builder` 로 재실행해 관측량이 소수 15자리까지 일치함을 "
        "확인했다(에너지 일관성 105.510722899358, 육방 NN거리 1.6170154574513436)."),
}

# ★ L5 이관(2026-08-05, bdbot.run 의 @RUN.builder 계약)으로 metrics.json 모양이 바뀌었다.
#   chain-bend·trap-drag는 이미 이 모양이었고 이 스크립트가 그 둘을 대조한 적이 없어서
#   드러나지 않았다 — trap-2d-5um·soft-r3 를 이관하면서 처음 걸렸다. 아래는 **정보 손실이
#   아니라 위치 이동**이거나(RENAMED로 값까지 대조), 이 스키마에서 의도적으로 뺀 것이다
#   (checks[]는 design/post_run 두 phase만 받는다 — post_checks는 result에 남긴다;
#   prediction_source 는 role/scope/derivation/note 로 대체돼 더 풍부해졌다).
#   물리 수치 자체는 옛 런과 소수 15자리까지 일치함을 별도로 확인했다(대화 기록 참조) —
#   여기서 "위반 없음"이 되는 것은 스키마가 같다는 뜻이 아니라 **삭제된 필드가 전부
#   설명 가능하다**는 뜻이다.
ALLOWED_DELETION_PREFIXES = (
    "physical.",             # 이제 star 파라미터만 담는다(trap-drag 관례). SI 값은
                              # spec.physical()로 스펙에서 언제든 재구성 가능하다.
    "checks[5]", "checks[6]",  # post_run 검사 → result.post_checks[] 로 이동
    "numerics.bias_predicted_pct", "numerics.dt_over_tau_B",
    "numerics.primary_sem_pct", "numerics.stat_target_pct",
    # soft-r3: 예측 없는 관측량(psi6·U_per_N)은 이제 result.psi6/result.pe_mean 에만
    # 있다 — observables[] 는 '측정 vs 예측' 두 항목만 남긴다 (정보 손실 아님).
    "observables[2]", "observables[3]",
)
ALLOWED_DELETION_SUFFIXES = (
    ".prediction_source",    # → role/scope/derivation/note (더 풍부한 구조)
)

# 리팩터로 바뀌어도 되는 필드
TOLERATED = {
    "run_id": "①에서 EXPECTED_RUN_ID 로 따로 판정한다 (여기서 중복 판정하지 않음)",
    "wall_seconds": "벽시계 — 매 실행 다름",
    "steps_per_second": "벽시계 파생",
    "numerics.x2_sem_pct": "블록 SEM을 float64로 업캐스트 (1.2e-6, 더 정확한 쪽)",
    "numerics.primary_sem_pct": "같은 이유",
    "observables[0].err_pct": (
        "부호 규약 통일 — **버그 수정**. 1-B 원본은 이 한 행만 (예측−측정)/|측정| 이라 "
        "같은 파일의 다른 행들(측정−기준)/기준 과도, 1-A와도 부호가 반대였다. "
        "분모도 |측정|→|예측| 로 바뀌어 크기가 두 값의 상대차(8.7e-5)만큼 달라진다. "
        "|오차|<2% 판정은 그대로."),
    "observables[0].prediction_source": (
        "1-B 원본은 이 필드 하나에 전체 설명 문장을 담았다. L5 스키마는 짧은 태그"
        "(source)와 전체 문장(note)·role·derivation을 분리한다 — 이 런은 그 분리가 "
        "생기기 **전** 코드로 만들어져 source가 'none'(기본값)으로 남아 있지만, "
        "현재 코드는 source='consistency'를 넘긴다(observables[0].note에 전체 문장이 "
        "그대로 있다 — 정보 손실 아님, 값 이관은 됐다)."),
    "observables[1].prediction_source": "같은 이유 — 현재 코드는 source='lattice'.",
}
# 1-C·앞단에서 **의도적으로 추가된** 스키마 필드. 삭제는 여전히 허용하지 않는다.
EXPECTED_NEW_FIELDS = {
    "schema": "metrics 스키마 버전 (bdbot.metrics/0.2)",
    "checks[].hard": "하드/소프트 구분 — 1-A는 전부 통과해서 구분이 안 드러났다",
    "checks[].phase": "design / post_run",
    "equilibration.source": "평형 지표를 케이스가 선언 (1-A의 앵커 변위는 확산계에 못 쓴다)",
    "equilibration.series_key": "같은 이유",
    "equilibration.label": "같은 이유",
    "numerics.dt_over_tau_B": "기준 시간 대비 dt 도 함께 기록",
    "numerics.primary_sem_pct": "케이스마다 1차 지표가 다르다 (x2_sem_pct 는 1-A 전용)",
    "numerics.stat_target_pct": "통계 목표를 케이스가 선언",
    "checks[].op": ("비교 방향(<=/>=). L3 스펙을 되읽을 때 필요하다 — 없으면 '>=' 검사"
                    "(관측창 등)가 '<=' 로 복원돼 판정이 뒤집힌다 (bdbot.nondim.load)."),
}

# 부호 규약만 바뀐 항목: 크기가 '두 값의 상대차' 수준을 넘어서면 실제 회귀로 본다
SIGN_ONLY = {"observables[0].err_pct": 1e-3}
RTOL = 1e-9

# ★ **이름만** 바뀐 필드 (옛 키 → 새 키). 값은 여전히 대조합니다 — 이름 변경을 빌미로
#   검사를 건너뛰면 "조용히 통과"가 됩니다.
#   L3 작업에서 `metrics.dimensionless` 의 키를 리포트 표시문자열에서 기호로 바꿨습니다.
#   metrics.json 은 postmortem 의 유일한 입력인데
#   `'k*     = k d²/kT   트랩 vs 열요동'` 같은 키로는 질의할 수 없었습니다.
RENAMED = {
    "dimensionless.k*     = k d²/kT   트랩 vs 열요동": "dimensionless.k*",
    "dimensionless.l_k/d  = 1/√k*     요동폭 vs 입자": "dimensionless.l_k/d",
    "dimensionless.tau_k/tau_B = 1/k*": "dimensionless.tau_k/tau_B",
    "dimensionless.dt/tau_k           적분 해상": "dimensionless.dt/tau_k",
    "dimensionless.T_obs/tau_k        관측창": "dimensionless.T_obs/tau_k",
    "dimensionless.Gamma  = U(a_mean)/kT   결합 vs 열요동 ★": "dimensionless.Gamma",
    "dimensionless.A      = U(d)/kT        접촉 결합": "dimensionless.A",
    "dimensionless.phi                     밀집도": "dimensionless.phi",
    "dimensionless.a_mean/d                평균간격": "dimensionless.a_mean/d",
    "dimensionless.L/d                     박스 크기": "dimensionless.L/d",
    "dimensionless.r_c/d                   컷오프": "dimensionless.r_c/d",
    "dimensionless.r_c/a_mean              컷오프(이웃 껍질 수)": "dimensionless.r_c/a_mean",
    "dimensionless.dt/tau_int              적분 해상": "dimensionless.dt/tau_int",
    "dimensionless.T_obs/tau_B             관측창": "dimensionless.T_obs/tau_B",
    "dimensionless.St     = tau_p/tau_B    관성 vs 확산": "dimensionless.St",
    "numerics.x2_sem_pct": "result.x2_sem_pct",
    # L5 이관(2026-08-05) — 케이스별 부가 정보가 top-level `structure`에서
    # RUN.execute()가 채우는 `result`로 옮겨졌다. 값은 그대로 대조한다.
    "structure.psi6": "result.psi6",
    "structure.psi6_sem": "result.psi6_sem",
    "structure.nn_distance_d": "result.nn_distance_d",
    "structure.nn_std_rel": "result.nn_std_rel",
    "structure.min_sep_d": "result.min_sep_d",
    "structure.Gamma": "result.Gamma",
    "structure.coord_hist": "result.coord_hist",
    "structure.state_predicted": "result.state_predicted",
    "structure.u_rms_rel_einstein": "result.u_rms_rel_einstein",
}


def flatten(d, pre=""):
    out = {}
    for k, v in d.items():
        p = f"{pre}{k}"
        if isinstance(v, dict):
            out.update(flatten(v, p + "."))
        elif isinstance(v, list) and v and isinstance(v[0], dict):
            for i, x in enumerate(v):
                out.update(flatten(x, f"{p}[{i}]."))
        else:
            out[p] = v
    return out


def run_id_of(script, extra) -> str | None:
    out = subprocess.run([PY, str(ROOT / script)] + extra + ["--report"],
                         capture_output=True, text=True, cwd=ROOT).stdout
    m = re.search(r"run_id=(\S+)", out)
    return m.group(1) if m else None


def compare(label, ref_path: Path, new_path: Path) -> tuple[bool, list[str]]:
    ref, new = json.loads(ref_path.read_text()), json.loads(new_path.read_text())
    fr, fn = flatten(ref), flatten(new)
    lines, bad = [], []
    same = 0
    renamed = 0
    for k in sorted(fr):
        key = k
        if key not in fn and key in RENAMED and RENAMED[key] in fn:
            key = RENAMED[k]          # 이름만 바뀜 — 값은 아래에서 그대로 대조한다
            renamed += 1
        if key not in fn:
            if k.startswith(ALLOWED_DELETION_PREFIXES) or k.endswith(ALLOWED_DELETION_SUFFIXES):
                lines.append(f"    ~ 삭제(허용): {k}")
            else:
                bad.append(f"    ✗ 삭제된 필드: {k}")
            continue
        a, b = fr[k], fn[key]
        if isinstance(a, (int, float)) and not isinstance(a, bool) \
                and isinstance(b, (int, float)) and not isinstance(b, bool):
            if a == b or abs(a - b) <= RTOL * max(abs(a), abs(b), 1e-300):
                same += 1
            elif k in TOLERATED:
                tol = SIGN_ONLY.get(k)
                if tol is not None and abs(abs(a) - abs(b)) > tol * max(abs(a), abs(b)):
                    bad.append(f"    ✗ {k}: 부호만 바뀐 게 아니라 크기도 다르다 "
                               f"({a:.12g} → {b:.12g})")
                else:
                    lines.append(f"    ~ {k}: {a:.12g} → {b:.12g}   ({TOLERATED[k]})")
            else:
                rel = 100 * (b - a) / abs(a) if a else float("nan")
                bad.append(f"    ✗ {k}: {a!r} → {b!r}  ({rel:+.4g}%)")
        elif a == b:
            same += 1
        elif k in TOLERATED:
            lines.append(f"    ~ {k}  ({TOLERATED[k]})")
        else:
            bad.append(f"    ✗ {k}:\n        전: {str(a)[:80]}\n        후: {str(b)[:80]}")
    added = sorted(set(fn) - set(fr) - set(RENAMED.values()))
    print(f"  동일 {same}개 · 이름만 변경 {renamed}개 · 허용된 변경 {len(lines)}개 · "
          f"추가 {len(added)}개 · 위반 {len(bad)}개")
    for ln in lines:
        print(ln)
    for ln in bad:
        print(ln)
    if added:
        groups = sorted({re.sub(r"\[\d+\]", "[]", a) for a in added})
        print(f"    + 추가된 필드({len(added)}): " + ", ".join(groups))
    return not bad, bad


def main() -> int:
    print("=" * 84)
    print("Phase 1-C 동등성 검사 — 리팩터 전(1-B 스냅샷) vs 현재 코드")
    print("=" * 84)
    ok_all = True
    for label, script, extra, ref_name in TARGETS:
        print(f"\n■ {label}   ({script} {' '.join(extra)})")
        rid = run_id_of(script, extra)
        ref_path = REF / ref_name
        ref = json.loads(ref_path.read_text())
        exp = EXPECTED_RUN_ID.get(ref["run_id"])
        if rid == ref["run_id"]:
            id_ok, note = True, "✓ 동일"
        elif exp and rid == exp[0]:
            id_ok, note = True, "✓ 의도된 변경"
        else:
            id_ok, note = False, "✗ 달라짐 (이유가 기록되지 않음)"
        ok_all &= id_ok
        print(f"  ① run_id  {rid}")
        print(f"     스냅샷  {ref['run_id']}   {note}")
        if exp and rid == exp[0]:
            print(f"     이유: {exp[1]}")
        new_path = ROOT / "runs" / (rid or "") / "metrics.json"
        if not new_path.exists():
            print(f"  ② 재실행 결과 없음 — 먼저 실행하세요:  $PY {script} {' '.join(extra)} --force")
            ok_all = False
            continue
        print("  ②③ metrics 대조")
        ok, _ = compare(label, ref_path, new_path)
        ok_all &= ok
    print()
    print("=" * 84)
    print("✓ PASS — 리팩터가 결과를 바꾸지 않았다" if ok_all else "✗ FAIL — 위 위반 항목 확인")
    print("=" * 84)
    return 0 if ok_all else 1


if __name__ == "__main__":
    sys.exit(main())
