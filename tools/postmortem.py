"""사후분석 — 완료된 런을 자동 진단하고 `record.json`(KB 엔트리)을 만든다.

마스터플랜 §10. 선언이 아니라 **측정**으로 성공/실패를 판정한다.
LLM 없이 동작한다 — 판정은 수치 지표로만 한다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY tools/postmortem.py runs/<run_id>
    $PY tools/postmortem.py runs/<run_id> --lesson "교훈" --kind pitfall

`record.json` 은 나중에 SQLite KB(§7)가 들어올 때의 입력이다. 지금은 파일 하나면 충분하고,
런이 100개를 넘거나 문헌이 들어오면 그때 DB로 옮긴다.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
SCHEMA = "bdbot.record/0.1"

# 실패 분류 체계 (마스터플랜 §10.1)
TAXONOMY = ["NUM_DIVERGE", "NUM_DRIFT", "EQ_INSUFFICIENT", "STAT_INSUFFICIENT",
            "FINITE_SIZE", "WRONG_REGIME", "RESOURCE", "SPEC_ERROR"]


def _tau_int(y: np.ndarray, c: float = 5.0) -> float:
    """적분 자기상관 시간 (스텝 단위, 자동 창 절단). n_eff = n/(2τ+1)."""
    y = np.asarray(y, dtype=float) - y.mean()
    n = len(y)
    if n < 8 or y.std() == 0:
        return 0.0
    nfft = 1 << (2 * n - 1).bit_length()
    F = np.fft.rfft(y, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    ac /= ac[0]
    tau, k = 0.0, 1
    while k < n:                                  # Sokal 자동 창
        tau += ac[k]
        if k >= c * (2 * tau + 1):
            break
        k += 1
    return max(0.0, float(tau))


def _stationarity(series: np.ndarray, steps: np.ndarray) -> dict:
    """정상성 두 지표: 전반/후반 z 검정 + 선형 추세 t 검정.

    ★ 자기상관 보정 필수. 1-B에서 발견: 보정 없이 t 검정을 하면 전 구간 변화가
      평균의 −0.026% 인 런이 t=−3.3 으로 '드리프트'로 잡혔다. 시계열 표본은
      상관되어 있어 naive SE가 유효 표본 수를 과대평가한다
      (skill bd-physics §5.1의 '오차막대는 블록 평균으로'와 같은 실수).
    """
    n = len(series)
    tau = _tau_int(series)
    infl = math.sqrt(2 * tau + 1)                 # SE 팽창 인자
    half = n // 2
    a, b = series[:half], series[half:]
    pooled = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) * infl
    z = float((b.mean() - a.mean()) / pooled) if pooled > 0 else 0.0
    x = (steps - steps.mean()) / (steps.std() or 1.0)
    slope, icept = np.polyfit(x, series, 1)
    resid = series - (slope * x + icept)
    se = resid.std(ddof=2) / math.sqrt(n) * infl
    span = float(slope * (x.max() - x.min()))     # 전 구간 변화량
    mean = float(series.mean())
    return {"equilibrium_z": z, "trend_t": float(slope / se) if se > 0 else 0.0,
            "first_half": float(a.mean()), "second_half": float(b.mean()),
            "tau_int_samples": tau, "n_eff": float(n / (2 * tau + 1)),
            "drift_span": span,
            "drift_span_rel_pct": 100 * span / abs(mean) if mean else None}


def series_diagnostics(run: Path, m: dict) -> dict | None:
    """케이스가 지정한 평형 지표 시계열로 정상성을 본다 (metrics.equilibration).

    ★ 1-B에서 필요해졌다. 아래 궤적 기반 진단은 '앵커로부터의 변위'를 쓰는데
      이건 트랩 계 전용이다. 유체는 확산하므로 변위가 무한히 자라 항상 EQ 실패가 된다.
      구조계의 표준 평형 지표는 퍼텐셜 에너지 시계열이다.
    """
    spec = m.get("equilibration")
    if not spec:
        # 케이스가 명시하지 않았으면 퍼텐셜 에너지 시계열을 찾아본다.
        # 구조계의 표준 평형 지표이고, 없으면 궤적 변위 폴백으로 넘어간다.
        npz0 = run / "observables.npz"
        if not npz0.exists():
            return None
        with np.load(npz0) as z0:
            if "pe" not in z0:
                return None
        spec = {"source": "observables.npz", "series_key": "pe",
                "label": "⟨U⟩/N [kT] (기본 추정)"}
    npz = run / spec.get("source", "observables.npz")
    if not npz.exists():
        return None
    with np.load(npz) as z:
        key = spec["series_key"]
        if key not in z:
            return None
        y = np.asarray(z[key], dtype=float)
    if len(y) < 16:
        return None
    d = _stationarity(y, np.arange(len(y), dtype=float))
    d.update({"available": True, "source": f"{spec['source']}:{spec['series_key']}",
              "label": spec.get("label", spec["series_key"]),
              "n_frames_production": len(y)})
    return d


def equilibrium_diagnostics(run: Path, n_eq_steps: int) -> dict:
    """궤적에서 평형·드리프트를 측정한다. traj_A.gsd 가 없으면 None 들로 채운다.

    ⚠️ 앵커(초기 위치)로부터의 변위를 쓰므로 **속박계(트랩) 전용**이다.
       확산하는 계에는 series_diagnostics 를 쓴다.
    """
    traj = run / "traj_A.gsd"
    if not traj.exists():
        return {"available": False}
    import gsd.hoomd

    with gsd.hoomd.open(str(traj), mode="r") as t:
        steps = np.array([f.configuration.step for f in t])
        anchors = np.array(t[0].particles.position)     # frame 0 = 초기 배치 = 앵커
        L = float(t[0].configuration.box[0])
        dim = int(t[0].configuration.dimensions)
        r2 = np.empty(len(t))
        for i, f in enumerate(t):
            d = np.array(f.particles.position) - anchors
            d[:, :2] -= L * np.round(d[:, :2] / L)      # 주기 축만 (bd-hoomd 함정 7·8)
            r2[i] = (d[:, :dim] ** 2).sum(axis=1).mean()

    prod = steps >= n_eq_steps
    if prod.sum() < 8:
        return {"available": False, "reason": "프로덕션 프레임 부족"}
    r2p, sp = r2[prod], steps[prod].astype(float)

    half = len(r2p) // 2
    a, b = r2p[:half], r2p[half:]
    pooled = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b))
    z = float((b.mean() - a.mean()) / pooled) if pooled > 0 else 0.0

    x = (sp - sp.mean()) / sp.std()
    slope, icept = np.polyfit(x, r2p, 1)
    resid = r2p - (slope * x + icept)
    se = resid.std(ddof=2) / math.sqrt(len(x))
    t_stat = float(slope / se) if se > 0 else 0.0

    return {"available": True, "n_frames_production": int(prod.sum()),
            "initial_rms_displacement": float(math.sqrt(r2[0])),
            "source": "traj_A.gsd:displacement_from_anchor", "label": "⟨Δr²⟩ (앵커 기준)",
            "equilibrium_z": z, "trend_t": t_stat,
            "r2_first_half": float(a.mean()), "r2_second_half": float(b.mean())}


def diagnose(run: Path) -> dict:
    m = json.loads((run / "metrics.json").read_text())
    num = m["numerics"]
    # 케이스가 평형 지표를 지정했으면 그걸 쓰고, 없으면 궤적 변위(속박계 전용)로 폴백
    diag = series_diagnostics(run, m) or equilibrium_diagnostics(run, num["n_eq"])

    findings, failure_modes, not_verified = [], [], []

    # ① 평형
    if diag.get("available"):
        src = diag.get("label", "?")
        if abs(diag["equilibrium_z"]) < 3:
            findings.append(f"평형 ✓ (전반/후반 z={diag['equilibrium_z']:+.2f}, 지표={src})")
        else:
            failure_modes.append("EQ_INSUFFICIENT")
            findings.append(f"평형 ✗ z={diag['equilibrium_z']:+.2f} (지표={src})")
        rel = diag.get("drift_span_rel_pct")
        neff = diag.get("n_eff")
        extra = (f", 전구간 {rel:+.3f}%" if rel is not None else "")
        extra += (f", n_eff={neff:.0f}/{diag.get('n_frames_production', 0)}"
                  if neff is not None else "")
        # 유의성만으로 판정하지 않는다: 크기가 0.5% 미만이면 물리적으로 무해하다고 본다
        significant = abs(diag["trend_t"]) >= 3
        material = rel is None or abs(rel) >= 0.5
        if not significant:
            findings.append(f"드리프트 없음 ✓ (t={diag['trend_t']:+.2f}{extra})")
        elif not material:
            findings.append(f"드리프트 유의하나 미미 ⚠ (t={diag['trend_t']:+.2f}{extra})"
                            f" — 크기 기준 0.5% 미만이므로 실패로 보지 않음")
        else:
            failure_modes.append("NUM_DRIFT")
            findings.append(f"드리프트 ✗ t={diag['trend_t']:+.2f}{extra}")
    else:
        not_verified.append("equilibrium_from_trajectory")

    # ② 분리 검사 — 하드/소프트를 구분한다 (bd-physics §4).
    #    1-B에서 필요해졌다: 통계·유한크기 검사는 ⚠ 경고이지 실행 거부 사유가 아니다.
    #    hard 키가 없는 옛 런은 전부 하드로 읽는다 (1-A 호환).
    bad_hard = [c for c in m["checks"] if not c["ok"] and c.get("hard", True)]
    bad_soft = [c for c in m["checks"] if not c["ok"] and not c.get("hard", True)]
    if bad_hard:
        failure_modes.append("SPEC_ERROR")
        findings += [f"분리 검사 실패(하드): {c['name']}" for c in bad_hard]
    if bad_soft:
        findings += [f"분리 검사 경고: {c['name']} = {c['value']:.3g} (기준 {c['limit']:g})"
                     for c in bad_soft]
        not_verified.append("소프트 검사 미충족: " + ", ".join(c["name"] for c in bad_soft))
    if not bad_hard:
        tight = [c for c in m["checks"] if c["ok"] and c["margin"] < 5]
        findings.append(f"분리 검사 {len(m['checks'])}종 중 하드 전부 통과 ✓"
                        + (f" (여유 부족 {len(tight)}건)" if tight else ""))

    # ③ 목표 달성 — 예측이 있는 관측량만 판정한다 (err_pct=None 은 예측 없음)
    predicted = [o for o in m["observables"] if o.get("err_pct") is not None]
    # predicted=0 인 관측량은 퍼센트 오차가 정의되지 않는다 — err_sigma(z-점수)로 따로 판정한다
    # (bdbot/metrics.py `observable(sigma=...)`). 둘 다 없으면 진짜 "판정 불가".
    sigma_judged = [o for o in m["observables"]
                    if o.get("err_sigma") is not None and o.get("err_pct") is None]
    n_nopred = len(m["observables"]) - len(predicted) - len(sigma_judged)
    if predicted:
        worst = max(predicted, key=lambda o: abs(o["err_pct"]))
        if all(abs(o["err_pct"]) < 5 for o in predicted):
            findings.append(f"관측량 {len(predicted)}종 예측과 일치 ✓ "
                            f"(최대 오차 {worst['err_pct']:+.2f}% @ {worst['name']})")
        else:
            failure_modes.append("WRONG_REGIME")
            findings.append(f"관측량 불일치 ✗ 최대 {worst['err_pct']:+.2f}% @ {worst['name']}")
    if sigma_judged:
        worst_s = max(sigma_judged, key=lambda o: abs(o["err_sigma"]))
        tol_s = worst_s.get("tol_sigma") or 3.0
        if all(abs(o["err_sigma"]) < (o.get("tol_sigma") or 3.0) for o in sigma_judged):
            findings.append(f"관측량 {len(sigma_judged)}종 0-예측과 통계적으로 일치 ✓ "
                            f"(최대 {worst_s['err_sigma']:+.2f}σ @ {worst_s['name']}, "
                            f"기준 {tol_s:g}σ)")
        else:
            failure_modes.append("WRONG_REGIME")
            findings.append(f"관측량 0-예측과 불일치 ✗ 최대 {worst_s['err_sigma']:+.2f}σ "
                            f"@ {worst_s['name']}")
    if n_nopred:
        not_verified.append(f"예측값 없는 관측량 {n_nopred}종 (측정만 기록, 판정 안 함)")

    # ④ 통계 + 편향 일관성
    sem = num.get("x2_sem_pct") or num.get("primary_sem_pct")
    # 통계 목표는 계마다 다르다. 케이스가 선언하면 그걸 쓰고, 없으면 0.5%(1-A 기본).
    target = num.get("stat_target_pct", 0.5)
    if sem is not None:
        if sem < target:
            findings.append(f"통계 충분 ✓ (±{sem:.3f}% < 목표 {target:g}%)")
        else:
            failure_modes.append("STAT_INSUFFICIENT")
            findings.append(f"통계 부족 ✗ (±{sem:.3f}% ≥ 목표 {target:g}%)")
        bp = num.get("bias_predicted_pct")
        x2 = next((o for o in m["observables"] if "x²" in o["name"] or "x2" in o["name"]), None)
        if bp is not None and x2 is not None:
            d_ = abs(x2["err_pct"] - bp)
            ok = d_ < 3 * sem
            findings.append(f"편향 법칙 {'✓' if ok else '✗'} 예측 {bp:+.3f}% vs "
                            f"측정 {x2['err_pct']:+.3f}% ({d_/sem:.1f} SEM)")
            if not ok:
                failure_modes.append("NUM_DRIFT")

    # ⑤ 직접 확인하지 않은 것 — 정직하게 남긴다
    conv = m.get("convergence_checked") or []
    if "dt" in conv:
        findings.append("dt 수렴 직접 확인 ✓ " + str(m.get("convergence_notes", {}).get("dt", "")))
    else:
        not_verified.append("dt_convergence_direct  (dt/2 재실행 안 함)")
    if "r_cut" in conv:
        findings.append("r_c 수렴 직접 확인 ✓ " + str(m.get("convergence_notes", {}).get("r_cut", "")))

    outcome = "success" if not failure_modes else (
        "partial" if len(failure_modes) == 1 else "failure")
    return {"metrics": m, "trajectory": diag, "findings": findings,
            "failure_modes": sorted(set(failure_modes)),
            "not_verified": not_verified, "outcome": outcome}


def build_record(run: Path, d: dict, lessons: list) -> dict:
    m = d["metrics"]
    return {
        "schema": SCHEMA,
        "run_id": m["run_id"],
        "case": m["case"],
        "kind": "our_run",
        "tier": 3,                       # 우리 시뮬레이션에서 귀납 (§5.2)
        "system_tags": m["system_tags"],
        "reference_scales": m["reference_scales"],
        "physical": m["physical"],
        "dimensionless": m["dimensionless"],
        "observables": m["observables"],
        "numerics": m["numerics"],
        "diagnostics": {k: v for k, v in d["trajectory"].items() if k != "available"},
        "findings": d["findings"],
        "outcome": d["outcome"],
        "failure_modes": d["failure_modes"],
        "not_verified": d["not_verified"],
        "lessons": lessons,
        "artifacts": sorted(f.name for f in run.iterdir()),
    }


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("run", type=Path)
    ap.add_argument("--lesson", action="append", default=[],
                    help='"주장::종류::좌표키=값,..." 형식 또는 그냥 주장')
    args = ap.parse_args()

    run = (ROOT / args.run) if not args.run.is_absolute() else args.run
    if not (run / "metrics.json").exists():
        print(f"✗ {run}/metrics.json 이 없습니다. 케이스 스크립트를 최신본으로 재실행하세요.")
        return 1

    d = diagnose(run)

    lessons = []
    for raw in args.lesson:
        parts = raw.split("::")
        claim = parts[0]
        kind = parts[1] if len(parts) > 1 else "method_note"
        coords = {}
        if len(parts) > 2 and parts[2]:
            for kv in parts[2].split(","):
                k_, v_ = kv.split("=")
                coords[k_.strip()] = float(v_)
        lessons.append({"claim": claim, "kind": kind, "coords": coords, "tier": 3})

    # ★ 교훈은 누적한다. 덮어쓰면 안 된다 — 1-C에서 사후분석을 재실행하며 6건을 날렸다.
    #   run_id 가 콘텐츠 주소라 같은 record 는 같은 스펙의 런이고, 이전 교훈이 그대로 유효하다.
    prev_path = run / "record.json"
    if prev_path.exists():
        try:
            prev = json.loads(prev_path.read_text()).get("lessons", [])
        except (json.JSONDecodeError, OSError):
            prev = []
        seen = {l.get("claim") for l in lessons}
        kept = [l for l in prev if l.get("claim") not in seen]
        if kept:
            print(f"  (이전 교훈 {len(kept)}건 유지)")
        lessons = kept + lessons

    rec = build_record(run, d, lessons)
    (run / "record.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False))

    print("=" * 78)
    print(f"사후분석 — {rec['run_id']}")
    print("=" * 78)
    for f_ in d["findings"]:
        print(f"  {f_}")
    if d["not_verified"]:
        print("\n  직접 확인하지 않은 것:")
        for nv in d["not_verified"]:
            print(f"    · {nv}")
    if lessons:
        print("\n  교훈 (tier 3):")
        for l_ in lessons:
            c = f"  {l_['coords']}" if l_["coords"] else ""
            print(f"    [{l_['kind']}] {l_['claim']}{c}")
    print("\n" + "=" * 78)
    print(f"OUTCOME: {d['outcome'].upper()}"
          + (f"   failure_modes={d['failure_modes']}" if d["failure_modes"] else ""))
    print(f"→ {(run / 'record.json').relative_to(ROOT)}")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
