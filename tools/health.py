"""L4 판정기를 완료된 런에 적용한다 (사후) — `bdbot.health` 의 CLI 어댑터.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY tools/health.py runs/<run_id>        # 한 건
    $PY tools/health.py --all                # 전부 훑기
    $PY tools/health.py --gate specs/x.json  # 실행 **전** 게이트만

`metrics.json` 의 `equilibration.series_key` 를 1순위로 보고, 없으면 `observables.npz`
안의 1차원 시계열을 자동으로 고릅니다. **물리적 옳고 그름은 판정하지 않습니다** —
발산·NaN·정지·붕괴, 그리고 L3 원장과의 대조뿐입니다.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bdbot import health as H          # noqa: E402
from bdbot import nondim as ND         # noqa: E402

ROOT = Path(__file__).resolve().parent.parent

# 평형·건전성 지표로 쓸 만한 1차원 시계열. 값이 항상 양수여야 하는 것은 positive=True.
# (positive, cumulative, 표시이름).  cumulative = 자라는 게 정상인 누적량
SERIES_HINTS = {
    "pe":               (False, False, "⟨U⟩/N"),
    "psi6":             (True,  False, "ψ₆"),
    "min_sep":          (True,  False, "최소 간격"),
    "x2":               (True,  False, "⟨x²⟩"),
    "cos_theta_series": (False, False, "⟨cos θ⟩"),
    "msd":              (True,  True,  "MSD"),
    "msad_folded":      (True,  True,  "MSAD"),
}
SKIP_PREFIX = ("rdf_", "psd_", "ac_", "px", "t", "disp_corr_", "final_")


def pick_series(npz_path: Path, metrics: dict) -> list[tuple[str, np.ndarray, bool]]:
    out = []
    if not npz_path.exists():
        return out
    with np.load(npz_path) as z:
        keys = list(z.files)
        eq = (metrics.get("equilibration") or {}).get("series_key")
        ordered = ([eq] if eq in keys else []) + [k for k in keys if k != eq]
        # 누적량의 성장 지수는 **실제 시간축**으로 재야 한다 (lag 가 로그 간격일 수 있음)
        tax = np.asarray(z["t"], float) if "t" in keys else None
        for k in ordered:
            if k in SKIP_PREFIX or any(k.startswith(p) for p in SKIP_PREFIX):
                continue
            y = z[k]
            if y.ndim != 1 or y.size < 20 or not np.issubdtype(y.dtype, np.number):
                continue
            pos, cum, label = SERIES_HINTS.get(k, (False, False, k))
            tt = tax if (cum and tax is not None and tax.size == y.size) else None
            out.append((label, np.asarray(y, float), pos, cum, tt))
            if len(out) >= 3:
                break
    return out


def judge_run(run: Path, verbose=True) -> H.HealthReport:
    rep = H.HealthReport()
    mfile = run / "metrics.json"
    metrics = json.loads(mfile.read_text()) if mfile.exists() else {}

    series = pick_series(run / "observables.npz", metrics)
    if not series:
        rep.add(True, None, "시계열", "판정 가능한 시계열 없음 — 생략")
    for label, y, pos, cum, tt in series:
        H.judge_series(label, y, rep, positive=pos, cumulative=cum, t=tt)

    # L3 대조 — 스펙이 있고 스텝 변위가 기록돼 있을 때만
    num = metrics.get("numerics", {})
    # ⓐ 힘 기반(`run.Guard` 가 런타임에 측정, 런 전체 최악값) ⓑ 위치 차분(하위호환).
    # ⓐ 를 우선합니다 — 열잡음을 뺄 필요가 없어 정확합니다 (health.step_health 참조).
    drift = num.get("step_drift_max_sigma")
    step_rms = num.get("step_rms_sigma")
    spec_path = _find_spec(metrics.get("run_id", run.name))
    pred = None
    if spec_path:
        try:
            pred = H.predicted_dt_over_tau(ND.load(spec_path))
        except Exception:
            pred = None
    dim = metrics.get("dimensions", 2)
    if drift is not None:
        H.step_health(None, num.get("dt_star", 0.0), dim, pred, rep, drift_direct=drift)
    elif step_rms is not None and num.get("dt_star"):
        H.step_health(step_rms, num["dt_star"], dim, pred, rep)
    else:
        # ⚠️ 여기에 들어오면 **이 모듈의 핵심 검사가 돌지 않은 것**입니다.
        #    HEALTHY 로 보이지만 스텝 해상은 검사되지 않았습니다 — 그걸 명시합니다.
        #    ★ 소급 측정은 불가능합니다: ① run_id 가 콘텐츠 주소라 현재 코드로 재실행하면
        #      **다른 id 의 새 런**이 생기고 이 런은 그대로 미측정으로 남습니다.
        #      ② GSD 재생으로 힘만 다시 계산하는 우회는 **시간 의존 구동에서 무효**입니다
        #      (트랩 앵커가 t=0 에 고정되어 실측 16배 과대 — scratch/probe_gsd_replay.py).
        rep.add(True, None, "스텝 변위 (미측정)",
                "이 런은 step_drift 측정이 배선된 2026-08-05 이전에 실행됐습니다 — "
                "스텝 해상이 **검사되지 않았습니다**. 소급 측정은 불가능합니다(새 런만 가능)")
        rep.measured["step_method"] = "none"
        if pred is not None:
            rep.measured["dt_over_tau_fast_predicted"] = pred
    return rep


def _find_spec(run_id: str) -> Path | None:
    p = ROOT / "specs" / f"{run_id}.json"
    if p.exists():
        return p
    hits = sorted((ROOT / "specs").glob(f"{run_id.split('__')[0]}*.json"))
    return hits[0] if hits else None


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("run", nargs="?", type=Path)
    ap.add_argument("--all", action="store_true")
    ap.add_argument("--gate", type=Path, help="실행 전 게이트만 (스펙 경로)")
    a = ap.parse_args()

    if a.gate:
        spec = ND.load(a.gate)
        probs = H.gate(spec)
        notes = H.gate_notes(spec)
        print(f"gate — {spec.run_id}   (L3 verdict: {spec.verdict})")
        for p in probs:
            print(f"  ✗ {p}")
        # 막지 않는 것도 반드시 보여줍니다 — 조용히 통과시키면 게이트가 무의미합니다.
        for n in notes:
            print(f"  ⚠ {n}")
        print("  OK — cleared to run" if not probs else "  -> RUN REFUSED")
        return 1 if probs else 0

    runs = sorted(p.parent for p in (ROOT / "runs").glob("*/metrics.json")) if a.all \
        else [(ROOT / a.run) if not a.run.is_absolute() else a.run]
    bad = 0
    unmeasured = 0
    for r in runs:
        rep = judge_run(r)
        m = rep.measured
        if m.get("step_method", "none") == "none":
            unmeasured += 1
        if a.all:
            extra = ""
            if "ledger_ratio" in m:
                extra = (f"  dt/τ 측정 {m['dt_over_tau_fast_measured']:.1e}"
                         f" / 예측 {m['dt_over_tau_fast_predicted']:.1e}"
                         f" = {m['ledger_ratio']:.2f}×")
            elif "dt_over_tau_fast_measured" in m:
                extra = f"  dt/τ 측정 {m['dt_over_tau_fast_measured']:.1e} (L3 예측 없음)"
            elif "dt_over_tau_fast_predicted" in m:
                extra = (f"  L3예측 dt/τ={m['dt_over_tau_fast_predicted']:.1e}"
                         f"  ⚠ 스텝 미검사")
            else:
                extra = "  ⚠ 스텝 미검사"
            print(f"  {'✓' if rep.verdict == 'HEALTHY' else '✗'} {r.name[:44]:<46}"
                  f"{rep.verdict:<11}{','.join(rep.failure_modes)}{extra}")
        else:
            print(rep.render())
        bad += rep.verdict != "HEALTHY"
    if a.all:
        n = len(runs)
        print(f"\n{n - bad}/{n} HEALTHY")
        # ★ 커버리지를 **따로** 보고합니다. "81/81 HEALTHY" 만 찍으면 핵심 검사가 한 번도
        #   돌지 않았는데 전부 통과한 것처럼 읽힙니다 — 침묵은 성공이 아닙니다.
        print(f"스텝 해상 측정: {n - unmeasured}/{n} 런")
        if unmeasured:
            print(f"  ⚠ 미측정 {unmeasured}런은 step_drift 측정이 배선된 2026-08-05 이전에 "
                  f"실행된 **레거시**입니다.")
            # ★ 2026-08-06 정정: "재실행하면 무조건 다른 id 가 된다"고 적어뒀는데 **틀렸습니다.**
            #   trap-drag 를 원래 인자(--traverse 0.117647)로 재실행하니 run_id 가 레거시와
            #   **정확히 일치**해 제자리에서 채워졌습니다(그리고 옛 데이터를 덮어썼습니다).
            #   같은 인자를 주면 채워지고, 코드가 그 케이스의 스펙을 바꿨으면 새 id 가 됩니다.
            print("    채우려면 **원래 CLI 인자 그대로** 재실행해야 합니다 — 그러면 run_id 가 "
                  "일치해 제자리에서 채워집니다.")
            print("    ⚠ 그때 옛 데이터(metrics·observables·GSD)는 덮여 사라집니다 "
                  "(record.json 만 보존). 필요하면 먼저 백업하세요.")
            print("    스펙이 바뀐 케이스는 새 id 가 생기고 레거시는 미측정으로 남습니다. "
                  "GSD 재생 우회는")
            print("    시간 의존 구동(이동 트랩·진동 구동)에서 무효입니다 (실측 16배 과대).")
            print("    이 런들의 HEALTHY 는 '발산·정지·붕괴 없음'이고 'dt 가 충분히 작다'가 "
                  "아닙니다.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
