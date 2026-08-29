"""Post-mortem -- diagnose a completed run automatically and write
`record.json` (a KB entry).

Success and failure are decided by **measurement**, not by declaration.
Runs with no LLM -- every verdict comes from a numerical indicator.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY tools/postmortem.py runs/<run_id>
    $PY tools/postmortem.py runs/<run_id> --lesson "the lesson" --kind pitfall

`record.json` is the input for a future SQLite KB. One file per run is enough for
now; move to a DB once runs exceed a hundred or the literature comes in.
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

# The failure classification
TAXONOMY = ["NUM_DIVERGE", "NUM_DRIFT", "EQ_INSUFFICIENT", "STAT_INSUFFICIENT",
            "FINITE_SIZE", "WRONG_REGIME", "RESOURCE", "SPEC_ERROR"]


def _tau_int(y: np.ndarray, c: float = 5.0) -> float:
    """Integrated autocorrelation time (in steps, automatic window truncation).

    n_eff = n/(2*tau+1).
    """
    y = np.asarray(y, dtype=float) - y.mean()
    n = len(y)
    if n < 8 or y.std() == 0:
        return 0.0
    nfft = 1 << (2 * n - 1).bit_length()
    F = np.fft.rfft(y, n=nfft)
    ac = np.fft.irfft(F * np.conj(F), n=nfft)[:n]
    ac /= ac[0]
    tau, k = 0.0, 1
    while k < n:                                  # Sokal automatic window
        tau += ac[k]
        if k >= c * (2 * tau + 1):
            break
        k += 1
    return max(0.0, float(tau))


def _stationarity(series: np.ndarray, steps: np.ndarray) -> dict:
    """Two stationarity indicators: a first-half/second-half z test plus a
    linear-trend t test.

    * Autocorrelation correction is mandatory. Found in the second case: without
      it, a run whose total change was -0.026% of the mean was flagged as 'drift'
      at t=-3.3. Time-series samples are correlated, so a naive SE overestimates
      the effective sample count (the same mistake as 'error bars come from block
      averaging' in skill bd-physics section 5.1).
    """
    n = len(series)
    tau = _tau_int(series)
    infl = math.sqrt(2 * tau + 1)                 # SE inflation factor
    half = n // 2
    a, b = series[:half], series[half:]
    pooled = math.sqrt(a.var(ddof=1) / len(a) + b.var(ddof=1) / len(b)) * infl
    z = float((b.mean() - a.mean()) / pooled) if pooled > 0 else 0.0
    x = (steps - steps.mean()) / (steps.std() or 1.0)
    slope, icept = np.polyfit(x, series, 1)
    resid = series - (slope * x + icept)
    se = resid.std(ddof=2) / math.sqrt(n) * infl
    span = float(slope * (x.max() - x.min()))     # total change across the window
    mean = float(series.mean())
    return {"equilibrium_z": z, "trend_t": float(slope / se) if se > 0 else 0.0,
            "first_half": float(a.mean()), "second_half": float(b.mean()),
            "tau_int_samples": tau, "n_eff": float(n / (2 * tau + 1)),
            "drift_span": span,
            "drift_span_rel_pct": 100 * span / abs(mean) if mean else None}


def series_diagnostics(run: Path, m: dict) -> dict | None:
    """Check stationarity on the equilibrium series the case named
    (metrics.equilibration).

    * Became necessary in the second case. The trajectory-based diagnosis below
      uses 'displacement from the anchor', which is trap-specific: a fluid
      diffuses, so the displacement grows without bound and it always fails the
      equilibrium check. The standard equilibrium indicator for a structural system
      is the potential-energy series.
    """
    spec = m.get("equilibration")
    if not spec:
        # If the case did not name one, look for a potential-energy series.
        # That is the standard indicator for a structural system; failing that, fall
        # through to the trajectory-displacement fallback.
        npz0 = run / "observables.npz"
        if not npz0.exists():
            return None
        with np.load(npz0) as z0:
            if "pe" not in z0:
                return None
        spec = {"source": "observables.npz", "series_key": "pe",
                "label": "<U>/N [kT] (default guess)"}
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
    """Measure equilibrium and drift from the trajectory. Fills with Nones when
    traj_A.gsd is absent.

    WARNING: this uses displacement from the anchor (the initial position), so it
       is **for bound (trap) systems only.** Use series_diagnostics for a diffusive
       system.
    """
    traj = run / "traj_A.gsd"
    if not traj.exists():
        return {"available": False}
    import gsd.hoomd

    with gsd.hoomd.open(str(traj), mode="r") as t:
        steps = np.array([f.configuration.step for f in t])
        anchors = np.array(t[0].particles.position)     # frame 0 = initial placement = the anchors
        L = float(t[0].configuration.box[0])
        dim = int(t[0].configuration.dimensions)
        r2 = np.empty(len(t))
        for i, f in enumerate(t):
            d = np.array(f.particles.position) - anchors
            d[:, :2] -= L * np.round(d[:, :2] / L)      # periodic axes only (bd-hoomd traps 7, 8)
            r2[i] = (d[:, :dim] ** 2).sum(axis=1).mean()

    prod = steps >= n_eq_steps
    if prod.sum() < 8:
        return {"available": False, "reason": "too few production frames"}
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
            "source": "traj_A.gsd:displacement_from_anchor", "label": "<dr^2> (from the anchor)",
            "equilibrium_z": z, "trend_t": t_stat,
            "r2_first_half": float(a.mean()), "r2_second_half": float(b.mean())}


def diagnose(run: Path) -> dict:
    m = json.loads((run / "metrics.json").read_text())
    num = m["numerics"]
    # Use the equilibrium indicator the case named; failing that, fall back to the
    # trajectory displacement (bound systems only)
    diag = series_diagnostics(run, m) or equilibrium_diagnostics(run, num["n_eq"])

    findings, failure_modes, not_verified = [], [], []

    # 1. equilibrium
    if diag.get("available"):
        src = diag.get("label", "?")
        if abs(diag["equilibrium_z"]) < 3:
            findings.append(f"equilibrated OK (half/half z={diag['equilibrium_z']:+.2f}, indicator={src})")
        else:
            failure_modes.append("EQ_INSUFFICIENT")
            findings.append(f"NOT equilibrated z={diag['equilibrium_z']:+.2f} (indicator={src})")
        rel = diag.get("drift_span_rel_pct")
        neff = diag.get("n_eff")
        extra = (f", total {rel:+.3f}%" if rel is not None else "")
        extra += (f", n_eff={neff:.0f}/{diag.get('n_frames_production', 0)}"
                  if neff is not None else "")
        # Do not judge on significance alone: below 0.5% in magnitude it is treated
        # as physically harmless
        significant = abs(diag["trend_t"]) >= 3
        material = rel is None or abs(rel) >= 0.5
        if not significant:
            findings.append(f"no drift OK (t={diag['trend_t']:+.2f}{extra})")
        elif not material:
            findings.append(f"drift significant but negligible ! (t={diag['trend_t']:+.2f}{extra})"
                            f" -- below 0.5% in magnitude, so not counted as a failure")
        else:
            failure_modes.append("NUM_DRIFT")
            findings.append(f"DRIFTING t={diag['trend_t']:+.2f}{extra}")
    else:
        not_verified.append("equilibrium_from_trajectory")

    # 2. separation checks -- distinguish hard from soft (bd-physics section 4).
    #    Became necessary in the second case: statistics and finite-size checks are
    #    warnings, not grounds for refusing to run.
    #    An older run with no `hard` key is read as all-hard (backward compatible).
    bad_hard = [c for c in m["checks"] if not c["ok"] and c.get("hard", True)]
    bad_soft = [c for c in m["checks"] if not c["ok"] and not c.get("hard", True)]
    if bad_hard:
        failure_modes.append("SPEC_ERROR")
        findings += [f"separation check FAILED (hard): {c['name']}" for c in bad_hard]
    if bad_soft:
        findings += [f"separation check warning: {c['name']} = {c['value']:.3g} (limit {c['limit']:g})"
                     for c in bad_soft]
        not_verified.append("soft checks unmet: " + ", ".join(c["name"] for c in bad_soft))
    if not bad_hard:
        tight = [c for c in m["checks"] if c["ok"] and c["margin"] < 5]
        findings.append(f"all hard checks passed OK of {len(m['checks'])} separation checks"
                        + (f" ({len(tight)} thin margin(s))" if tight else ""))

    # 3. targets met -- judge only observables that have a prediction
    #    (err_pct=None means no prediction)
    predicted = [o for o in m["observables"] if o.get("err_pct") is not None]
    # An observable with predicted=0 has no defined percentage error -- judge it
    # separately by err_sigma (a z-score), see bdbot/metrics.py
    # `observable(sigma=...)`. With neither, it is genuinely undecidable.
    sigma_judged = [o for o in m["observables"]
                    if o.get("err_sigma") is not None and o.get("err_pct") is None]
    n_nopred = len(m["observables"]) - len(predicted) - len(sigma_judged)
    if predicted:
        worst = max(predicted, key=lambda o: abs(o["err_pct"]))
        if all(abs(o["err_pct"]) < 5 for o in predicted):
            findings.append(f"{len(predicted)} observable(s) agree with prediction OK "
                            f"(worst error {worst['err_pct']:+.2f}% @ {worst['name']})")
        else:
            failure_modes.append("WRONG_REGIME")
            findings.append(f"observable MISMATCH worst {worst['err_pct']:+.2f}% @ {worst['name']}")
    if sigma_judged:
        worst_s = max(sigma_judged, key=lambda o: abs(o["err_sigma"]))
        tol_s = worst_s.get("tol_sigma") or 3.0
        if all(abs(o["err_sigma"]) < (o.get("tol_sigma") or 3.0) for o in sigma_judged):
            findings.append(f"{len(sigma_judged)} observable(s) statistically agree with a "
                            f"zero prediction OK (worst {worst_s['err_sigma']:+.2f} sigma @ "
                            f"{worst_s['name']}, limit {tol_s:g} sigma)")
        else:
            failure_modes.append("WRONG_REGIME")
            findings.append(f"observable MISMATCH against a zero prediction, worst "
                            f"{worst_s['err_sigma']:+.2f} sigma "
                            f"@ {worst_s['name']}")
    if n_nopred:
        not_verified.append(f"{n_nopred} observable(s) with no prediction (recorded as "
                            f"measurements, not judged)")

    # 4. statistics plus bias consistency
    sem = num.get("x2_sem_pct") or num.get("primary_sem_pct")
    # The statistical target differs per system. Use what the case declares; the
    # default 0.5% comes from the first case and must not be assumed elsewhere.
    target = num.get("stat_target_pct", 0.5)
    if sem is not None:
        if sem < target:
            findings.append(f"statistics sufficient OK (+/-{sem:.3f}% < target {target:g}%)")
        else:
            failure_modes.append("STAT_INSUFFICIENT")
            findings.append(f"statistics INSUFFICIENT (+/-{sem:.3f}% >= target {target:g}%)")
        bp = num.get("bias_predicted_pct")
        x2 = next((o for o in m["observables"] if "x²" in o["name"] or "x2" in o["name"]), None)
        if bp is not None and x2 is not None:
            d_ = abs(x2["err_pct"] - bp)
            ok = d_ < 3 * sem
            findings.append(f"bias law {'OK' if ok else 'MISMATCH'} predicted {bp:+.3f}% vs "
                            f"measured {x2['err_pct']:+.3f}% ({d_/sem:.1f} SEM)")
            if not ok:
                failure_modes.append("NUM_DRIFT")

    # 5. what was not directly confirmed -- recorded honestly
    conv = m.get("convergence_checked") or []
    if "dt" in conv:
        findings.append("dt convergence directly confirmed OK " + str(m.get("convergence_notes", {}).get("dt", "")))
    else:
        not_verified.append("dt_convergence_direct  (no dt/2 re-run)")
    if "r_cut" in conv:
        findings.append("r_c convergence directly confirmed OK " + str(m.get("convergence_notes", {}).get("r_cut", "")))

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
        "tier": 3,                       # induced from our own simulation
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
                    help='"claim::kind::coordkey=value,..." or just the claim')
    args = ap.parse_args()

    run = (ROOT / args.run) if not args.run.is_absolute() else args.run
    if not (run / "metrics.json").exists():
        print(f"x {run}/metrics.json is missing. Re-run the case script from the current version.")
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

    # * Lessons accumulate. Never overwrite -- re-running the post-mortem once
    #   destroyed 6 of them. run_id is content-addressed, so the same record means a
    #   run of the same spec, and the previous lesson is still valid.
    prev_path = run / "record.json"
    if prev_path.exists():
        try:
            prev = json.loads(prev_path.read_text()).get("lessons", [])
        except (json.JSONDecodeError, OSError):
            prev = []
        seen = {l.get("claim") for l in lessons}
        kept = [l for l in prev if l.get("claim") not in seen]
        if kept:
            print(f"  ({len(kept)} previous lesson(s) preserved)")
        lessons = kept + lessons

    rec = build_record(run, d, lessons)
    (run / "record.json").write_text(json.dumps(rec, indent=2, ensure_ascii=False))

    print("=" * 78)
    print(f"post-mortem -- {rec['run_id']}")
    print("=" * 78)
    for f_ in d["findings"]:
        print(f"  {f_}")
    if d["not_verified"]:
        print("\n  not directly confirmed:")
        for nv in d["not_verified"]:
            print(f"    · {nv}")
    if lessons:
        print("\n  lessons (tier 3):")
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
