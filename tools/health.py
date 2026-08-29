"""Apply the L4 judge to a completed run (post hoc) -- the CLI adapter for
`bdbot.health`.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY tools/health.py runs/<run_id>        # one run
    $PY tools/health.py --all                # sweep everything
    $PY tools/health.py --gate specs/x.json  # the **pre-run** gate only

Looks first at `metrics.json`'s `equilibration.series_key`, and failing that picks
a 1-D time series out of `observables.npz` automatically. **It does not judge
physical correctness** -- divergence, NaN, stalling, collapse, and the comparison
against the L3 ledger, and nothing else.
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

# 1-D series usable as an equilibrium or health indicator. positive=True when the
# value must always be positive.
# (positive, cumulative, display name).  cumulative = an accumulated quantity for
# which growing is normal
SERIES_HINTS = {
    "pe":               (False, False, "⟨U⟩/N"),
    "psi6":             (True,  False, "ψ₆"),
    "min_sep":          (True,  False, "min separation"),
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
        # the growth exponent of a cumulative quantity must be measured against the
        # **real time axis** (the lags may be logarithmically spaced)
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
        rep.add(True, None, "time series", "no judgeable series -- skipped")
    for label, y, pos, cum, tt in series:
        H.judge_series(label, y, rep, positive=pos, cumulative=cum, t=tt)

    # the L3 comparison -- only when a spec exists and a step displacement was recorded
    num = metrics.get("numerics", {})
    # (a) force-based (measured at runtime by `run.Guard`, the worst value over the
    #     whole run)  (b) position difference (backward compatible).
    # Prefer (a) -- there is no thermal noise to subtract, so it is exact
    # (see health.step_health).
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
        # WARNING: reaching here means **this module's core check did not run.**
        #    It will look HEALTHY, but the step resolution was not checked -- say so
        #    explicitly.
        #    * Retroactive measurement is impossible: (1) run_id is content-addressed,
        #      so re-running under the current code produces **a new run with a
        #      different id** and this one stays unmeasured. (2) The workaround of
        #      replaying the GSD to recompute forces is **invalid under
        #      time-dependent driving** (the trap anchor is pinned at t=0, measured
        #      16x overestimate -- verify/probe_gsd_replay.py).
        rep.add(True, None, "step displacement (NOT MEASURED)",
                "this run predates the step_drift measurement being wired (2026-08-05) "
                "-- the step resolution was **not checked**. Retroactive measurement "
                "is impossible (new runs only)")
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
    ap.add_argument("--gate", type=Path, help="the pre-run gate only (a spec path)")
    a = ap.parse_args()

    if a.gate:
        spec = ND.load(a.gate)
        probs = H.gate(spec)
        notes = H.gate_notes(spec)
        print(f"gate — {spec.run_id}   (L3 verdict: {spec.verdict})")
        for p in probs:
            print(f"  ✗ {p}")
        # Always show what does not block -- a gate that passes silently is not a gate.
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
                extra = (f"  dt/tau measured {m['dt_over_tau_fast_measured']:.1e}"
                         f" / predicted {m['dt_over_tau_fast_predicted']:.1e}"
                         f" = {m['ledger_ratio']:.2f}×")
            elif "dt_over_tau_fast_measured" in m:
                extra = f"  dt/tau measured {m['dt_over_tau_fast_measured']:.1e} (no L3 prediction)"
            elif "dt_over_tau_fast_predicted" in m:
                extra = (f"  L3 predicted dt/tau={m['dt_over_tau_fast_predicted']:.1e}"
                         f"  ! step NOT checked")
            else:
                extra = "  ! step NOT checked"
            print(f"  {'✓' if rep.verdict == 'HEALTHY' else '✗'} {r.name[:44]:<46}"
                  f"{rep.verdict:<11}{','.join(rep.failure_modes)}{extra}")
        else:
            print(rep.render())
        bad += rep.verdict != "HEALTHY"
    if a.all:
        n = len(runs)
        print(f"\n{n - bad}/{n} HEALTHY")
        # * Report coverage **separately.** Printing only "81/81 HEALTHY" reads as
        #   everything passing when the core check never ran once -- silence is not
        #   success.
        print(f"step resolution measured: {n - unmeasured}/{n} runs")
        if unmeasured:
            print(f"  ! the {unmeasured} unmeasured run(s) are **legacy** -- executed "
                  f"before the step_drift measurement was wired (2026-08-05).")
            # * Correction (2026-08-06): this used to say "re-running always gives a
            #   different id", which was **wrong.** Re-running trap-drag with its
            #   original arguments (--traverse 0.117647) produced a run_id that
            #   **matched the legacy one exactly** and filled it in place (overwriting
            #   the old data). Same arguments -> filled in place; if the code changed
            #   that case's spec -> a new id.
            print("    to fill them, re-run with **the original CLI arguments** -- then "
                  "the run_id matches and it fills in place.")
            print("    ! the old data (metrics, observables, GSD) is overwritten and "
                  "lost when that happens (only record.json is preserved). Back it up "
                  "first if needed.")
            print("    a case whose spec changed gets a new id, and the legacy run stays "
                  "unmeasured. The GSD-replay workaround is")
            print("    invalid under time-dependent driving (a moving trap, oscillatory "
                  "driving) -- measured 16x overestimate.")
            print("    HEALTHY for these runs means 'no divergence, no stall, no "
                  "collapse' and NOT 'dt is small enough'.")
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())
