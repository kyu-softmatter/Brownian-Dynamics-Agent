"""`metrics.json` — 기계가 읽는 런 결과. 사후분석·KB 환류의 **유일한** 입력.

두 케이스가 똑같은 최상위 키를 썼습니다. `tools/postmortem.py`가 이 스키마에 의존하므로
필드를 바꿀 때는 거기도 함께 봐야 합니다 (하위호환은 postmortem이 `.get`으로 처리).

1-B에서 추가된 세 필드 — 전부 "케이스마다 다른 것을 케이스가 선언한다"는 교훈에서 나왔습니다:
  · `equilibration`  평형 판정에 쓸 시계열. 1-A의 '앵커 변위'는 확산계에 못 씀
  · `checks[].hard`  하드/소프트 구분. 1-A는 전부 통과해서 구분이 안 드러났음
  · `numerics.stat_target_pct`  통계 목표. 1-A의 0.5%를 다른 계에 쓰면 안 됨
"""
from __future__ import annotations

import json
from pathlib import Path

SCHEMA = "bdbot.metrics/0.3"

# ★ 관측량의 **역할** — 2026-08-04 추가. 이걸 안 가르면 판정이 틀린다.
#
#   implementation_check  예측이 **내가 구현한 모델에서 해석적으로 따라 나온다.**
#                         일치는 "코드가 맞다"는 뜻이고 물리 발견이 아니다 (거의 순환).
#                         불일치 = **버그** → FAIL.
#   hypothesis            예측이 시뮬레이션이 **부과하지 않은 가정** 위에 서 있다
#                         (연속체 근사·희박 극한·유효매질·문헌 모델 등).
#                         불일치 = **결과**. FAIL 이 아니다 — 그게 알고 싶던 것이다.
#   measurement           예측이 없다. 시뮬레이션이 답이다.
#
# 왜 필요한가: 이 계들을 계산하는 이유는 기존 가설과 **다를 수 있기** 때문이다.
# 그런데 판정 로직이 "예측과 다르면 FAIL" 이면, 발견을 실패로 부르게 된다.
ROLES = ("implementation_check", "hypothesis", "measurement")


# 구성 수준 (마스터플랜 원칙 9.1). 단독 모듈과 조합 전체는 인식론적 지위가 다르다.
#   module     이 모듈만 켠 최소 구성 — 기존 이론이 곧 내 모델이므로 implementation_check 정당
#   composite  여러 모듈 조합 — 해석해가 대개 없다. 기존 이론을 갖다 대면 위험하다
SCOPES = ("module", "composite")


def observable(name, measured, predicted=None, unit="1", source="none",
               role="measurement", tol_pct=None, sigma=None, tol_sigma=None, note="",
               scope="composite", derivation="") -> dict:
    """관측량 한 개.

    `role` 을 반드시 생각해서 넣으세요 (위 ROLES). 기본값은 가장 보수적인 `measurement`
    (판정하지 않음) 입니다 — 실수로 발견을 실패로 부르지 않도록.
    `tol_pct` 는 `implementation_check` 에서만 의미가 있습니다.

    ★ `predicted=0` 처럼 예측이 정확히 0이면 퍼센트 오차가 정의되지 않습니다(0으로 나눔).
      이때는 `sigma`(측정의 표준오차, 예: 블록평균 SEM)를 같이 넘기면 z-점수
      `err_sigma = (measured-predicted)/sigma` 로 판정합니다 (기본 허용 `tol_sigma=3.0`).
      `sigma` 없이 `predicted=0`을 쓰면 판정 불가(measurement 취급) — 조용히 틀리지 않도록
      `judge()`가 그렇게 처리합니다.

    `scope` 는 이 관측량이 **단독 모듈**에서 나온 것인지 **조합 전체**에서 나온 것인지입니다
    (원칙 9.1). 기본값 `composite` 가 보수적입니다.

    ★ `composite` + `implementation_check` 는 `derivation` 을 요구합니다.
      조합에는 해석해가 대개 없습니다 — 그래서 시뮬레이션하는 것입니다. 그런데도 구현 검사를
      걸겠다면 **왜 그 식이 조합에도 유도되는지**(보통 희박·선형응답·단시간 같은 극한) 를
      적어야 합니다. 안 적으면 기존 이론을 조합에 갖다 대는 것이고, 맞으면 "검증됐다",
      안 맞으면 "이 영역엔 안 맞는다" 가 되어 **어느 쪽이든 배우는 게 없습니다.**

    예측은 **결과를 보기 전에** 고정되어야 합니다 (원칙 9.2). 구조적 보장은
    예측 함수가 시뮬레이션 결과를 인자로 받지 않는 것입니다 — `cases/*.py` 의
    `analytic(ledger)` 패턴.
    """
    if role not in ROLES:
        raise ValueError(f"role 은 {ROLES} 중 하나여야 합니다 (받은 값: {role!r})")
    if scope not in SCOPES:
        raise ValueError(f"scope 는 {SCOPES} 중 하나여야 합니다 (받은 값: {scope!r})")
    if scope == "composite" and role == "implementation_check" and not derivation:
        raise ValueError(
            f"[{name}] composite + implementation_check 에는 derivation 이 필요합니다 "
            f"(원칙 9.1). 조합에도 그 해석식이 유도되는 근거 — 보통 극한 — 를 적으세요. "
            f"근거를 못 대겠으면 role='hypothesis' 가 맞습니다: 불일치가 버그가 아니라 결과입니다.")
    err = None
    if predicted not in (None, 0) and measured is not None:
        err = 100.0 * (float(measured) - float(predicted)) / abs(float(predicted))
    err_sigma = None
    if sigma and predicted is not None and measured is not None:
        err_sigma = (float(measured) - float(predicted)) / float(sigma)
    return {"name": name, "measured": None if measured is None else float(measured),
            "predicted": None if predicted is None else float(predicted),
            "unit": unit, "err_pct": err, "err_sigma": err_sigma,
            "sigma": None if sigma is None else float(sigma),
            "prediction_source": source,
            "role": role, "scope": scope, "derivation": derivation,
            "tol_pct": tol_pct, "tol_sigma": tol_sigma, "note": note}


def judge(observables) -> tuple:
    """(판정, 실패한 구현검사, 가설과 어긋난 것, 측정만) — 역할별로 다르게 다룹니다.

    ★ `hypothesis` 가 어긋나도 FAIL 이 아닙니다. 보고할 **결과**입니다.
    """
    bad_impl, dev_hypo, meas = [], [], []
    for o in observables:
        role, err = o.get("role", "measurement"), o.get("err_pct")
        err_sigma = o.get("err_sigma")
        if role == "implementation_check":
            if err_sigma is not None:
                tol_s = o.get("tol_sigma") or 3.0
                if abs(err_sigma) > tol_s:
                    bad_impl.append(o)
            else:
                tol = o.get("tol_pct") or 5.0
                if err is None or abs(err) > tol:
                    bad_impl.append(o)
        elif role == "hypothesis":
            if err is not None and abs(err) > (o.get("tol_pct") or 10.0):
                dev_hypo.append(o)
        else:
            meas.append(o)
    # 원칙 9.1: 조합 전체를 구현 검사한 항목은 판정의 무게가 다르다 — 표시해 둔다
    comp_impl = [o for o in observables
                 if o.get("role") == "implementation_check"
                 and o.get("scope", "composite") == "composite"]
    if bad_impl:
        v = f"FAIL — 구현 검사 {len(bad_impl)}건 불일치 (버그)"
    elif dev_hypo:
        v = f"PASS (구현 정상) · 가설과 다름 {len(dev_hypo)}건 ← 결과"
    else:
        v = "PASS"
    if comp_impl:
        v += f"  [조합 구현검사 {len(comp_impl)}건 — 원칙 9.1: derivation 확인]"
    return v, bad_impl, dev_hypo, meas


def build(*, run_id, case, system_tags, reference_scales, physical, dimensionless,
          checks, observables, numerics, equilibration, wall_seconds,
          steps_per_second=None, extra=None) -> dict:
    """`checks`는 `(Check, phase)` 튜플 목록. phase = design | post_run."""
    m = {
        "schema": SCHEMA,
        "run_id": run_id,
        "case": case,
        "system_tags": list(system_tags),
        "reference_scales": dict(reference_scales),
        "physical": dict(physical),
        "dimensionless": {k: float(v) for k, v in dimensionless.items()},
        "checks": [c.as_dict(phase) for c, phase in checks],
        "observables": list(observables),
        "equilibration": dict(equilibration),
        "numerics": dict(numerics),
        "wall_seconds": float(wall_seconds),
    }
    if steps_per_second is not None:
        m["steps_per_second"] = float(steps_per_second)
    if extra:
        m.update(extra)
    return m


def write(outdir: Path, metrics: dict) -> Path:
    p = Path(outdir) / "metrics.json"
    p.write_text(json.dumps(metrics, indent=2))
    return p


def equilibration_series(series_key: str, label: str, source: str = "observables.npz") -> dict:
    """평형 판정 지표 선언. 케이스가 자기 계에 맞는 시계열을 지정합니다.

    속박계(트랩): 앵커로부터의 변위  ·  확산/구조계: 퍼텐셜 에너지 ⟨U⟩/N
    """
    return {"source": source, "series_key": series_key, "label": label}


__all__ = ["SCHEMA", "ROLES", "build", "write", "observable", "judge",
           "equilibration_series"]
