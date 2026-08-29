"""Adversarial checks on the L4 numerical-health judge.

CLAUDE.md working practice: "when you build a checker, deliberately break it and
see -- silently passing and not checking are different things." Each failure mode
is manufactured, and the test is whether it is caught **as exactly that mode**.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_health.py
"""
from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bdbot import health as H          # noqa: E402
from bdbot import nondim as ND         # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
rng = np.random.default_rng(4)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'✓' if cond else '✗'} {name}" + (f"   {detail}" if detail else ""))


def modes_of(series, **kw):
    r = H.HealthReport()
    H.judge_series("y", series, r, **kw)
    return r.failure_modes, r


print("=" * 78)
print("(1) time-series verdict -- each failure mode is manufactured")
print("=" * 78)

n = 400
healthy = 1.0 + 0.05 * rng.normal(size=n)
m, _ = modes_of(healthy)
check("a healthy series passes", m == [], f"modes={m}")

bad = healthy.copy(); bad[137] = np.nan
m, _ = modes_of(bad)
check("NaN → NUM_NONFINITE", m == ["NUM_NONFINITE"], f"modes={m}")

bad = healthy.copy(); bad[200] = np.inf
m, _ = modes_of(bad)
check("Inf → NUM_NONFINITE", m == ["NUM_NONFINITE"], f"modes={m}")

explode = np.exp(np.linspace(0, 25, n)) * (1 + 0.01 * rng.normal(size=n))
m, _ = modes_of(explode)
check("exponential blow-up -> NUM_DIVERGE", "NUM_DIVERGE" in m, f"modes={m}")

const = np.full(n, 3.14159)
m, _ = modes_of(const)
check("constant -> NUM_FROZEN", "NUM_FROZEN" in m, f"modes={m}")

collapse = np.abs(healthy) * np.exp(-np.linspace(0, 60, n))
m, _ = modes_of(collapse, positive=True)
check("collapse to 0 -> NUM_COLLAPSE", "NUM_COLLAPSE" in m, f"modes={m}")

m, _ = modes_of(collapse, positive=False)
check("with positive=False, collapse is not judged", "NUM_COLLAPSE" not in m,
      f"modes={m}")

m, _ = modes_of(healthy[:5])
check("too few samples -> judgment skipped (passes)", m == [], f"modes={m}")

# ★ A false positive that real data caught -- a regression test is written from the
#   data that actually broke.
#   A synthetic LINEAR MSD has a back/front ratio of only 7x, which cannot reproduce
#   the false positive. The real abp-rod MSD samples lag LOGARITHMICALLY, giving
#   1010x, and that is what was misjudged as NUM_DIVERGE.
# ★ The fixture was made **self-contained and synthetic** (2026-08-05). It used to
#   pick a real npz out of `runs/abp-rod*`, and when that run disappeared (run
#   directories change with re-runs and cleanups) the back/front ratio fell from
#   1010x to 343x, so **the false positive could no longer be reproduced and the
#   test broke.**
#   A regression test that depends on a particular run existing cannot guard the
#   regression.
#   The cause was never "real data" -- it was that **the lag is log-spaced** -- so
#   synthesising only that is sufficient: log lag + a diffusive (alpha=1) MSD.
t_log = np.logspace(0, 4.5, 40)                     # log-spaced lag (the same
                                                    # structure as the real npz)
msd_log = 4.0 * t_log * (1 + 0.02 * rng.normal(size=t_log.size))
q_ = max(2, msd_log.size // 4)
ratio_log = np.abs(msd_log[-q_:]).mean() / np.abs(msd_log[:q_]).mean()
check("a log-lag MSD exceeds a back/front ratio of 1e3 (reproduces the false-positive "
      "condition)", ratio_log > 1e3,
      f"back/front = {ratio_log:.3g}x  (a linear lag gives only 7x, which does not "
      f"reproduce it)")

m, _ = modes_of(msd_log, positive=True)
check("treating an MSD as a steady-state series gives a false positive (regression)",
      "NUM_DIVERGE" in m,
      f"back/front = {ratio_log:.3g}x -> {m}")

m, _ = modes_of(msd_log, positive=True, cumulative=True)
check("cumulative alone, with no time axis, is not enough (log lag)",
      "NUM_DIVERGE" in m,
      f"alpha measured against the index exceeds the limit -> {m}")

m, r_ = modes_of(msd_log, positive=True, cumulative=True, t=t_log)
# ★ Match what health.py actually names this finding. It filtered for the Korean
#   "성장 지수" while health.py has emitted "<name> growth exponent" since it was
#   translated, so this list was always empty. It only feeds the message, not the
#   assertion, so it degraded the diagnostic rather than breaking the check --
#   but the same coupling broke two real assertions elsewhere in this file.
alpha = [f.detail for f in r_.findings if "growth exponent" in f.name]
check("supplying the real time axis makes it pass", "NUM_DIVERGE" not in m,
      f"{alpha[0] if alpha else m}")

# If a real run happens to be present, compare against it **as well** (nice to have;
# the synthetic check above guards the regression on its own)
real = sorted((ROOT / "runs").glob("abp-rod*/observables.npz"))
real = [f for f in real if "msd" in np.load(f).files]
if real:
    with np.load(real[0]) as z:
        msd_real, t_real = np.asarray(z["msd"], float), np.asarray(z["t"], float)
    m, r_ = modes_of(msd_real, positive=True, cumulative=True, t=t_real)
    check("a real abp-rod MSD also passes once given the time axis",
          "NUM_DIVERGE" not in m,
          f"{real[0].parent.name[:40]} -> {m or 'no modes'}")
else:
    print("      (no real abp-rod run present -- the synthetic check alone still "
          "guards the regression)")

t = np.arange(1, n + 1, dtype=float)
msd_diffusive = 4.0 * t * (1 + 0.02 * rng.normal(size=n))       # α=1
m, _ = modes_of(msd_diffusive, positive=True, cumulative=True)
check("a diffusive MSD (alpha ~ 1) passes", m == [], f"modes={m}")

msd_ballistic = 9.0 * t ** 2 * (1 + 0.02 * rng.normal(size=n))  # α=2
m, _ = modes_of(msd_ballistic, positive=True, cumulative=True)
check("a ballistic MSD (alpha = 2) passes too", m == [], f"modes={m}")

msd_blowup = t ** 3.4                              # alpha=3.4 -- physically impossible
m, _ = modes_of(msd_blowup, positive=True, cumulative=True)
check("super-ballistic alpha=3.4 -> NUM_DIVERGE", "NUM_DIVERGE" in m, f"modes={m}")

print()
print("=" * 78)
print("(2) step displacement -> measuring dt/tau_fast + comparing against the L3 "
      "ledger")
print("=" * 78)

dt_star, dim = 1e-4, 2
thermal = math.sqrt(2 * dim * dt_star)

r = H.HealthReport()
H.step_health(thermal, dt_star, dim, predicted_dt_over_tau=2e-3, rep=r)
check("a purely thermal step -> drift 0, passes", r.verdict == "HEALTHY",
      f"drift={r.measured['dt_over_tau_fast_measured']:.2e}")

drift_true = 5e-3
r = H.HealthReport()
H.step_health(math.hypot(thermal, drift_true), dt_star, dim, 4e-3, r)
got = r.measured["dt_over_tau_fast_measured"]
check("drift is recovered separately from the thermal part",
      abs(got - drift_true) / drift_true < 0.02,
      f"injected {drift_true:.1e} -> recovered {got:.3e}")

r = H.HealthReport()
H.step_health(math.hypot(thermal, 5e-2), dt_star, dim, 4e-2, r)
check("step too large -> NUM_STEP_TOO_COARSE",
      "NUM_STEP_TOO_COARSE" in r.failure_modes,
      f"modes={r.failure_modes}")

r = H.HealthReport()
H.step_health(math.hypot(thermal, 8e-3), dt_star, dim, predicted_dt_over_tau=1e-4, rep=r)
check("measured >> L3 prediction -> LEDGER_INCOMPLETE",
      "LEDGER_INCOMPLETE" in r.failure_modes,
      f"measured/predicted = {r.measured['ledger_ratio']:.0f}x  <- a timescale "
      f"missing from the ledger")

r = H.HealthReport()
H.step_health(math.hypot(thermal, 4e-3), dt_star, dim, predicted_dt_over_tau=4e-3, rep=r)
check("measured ~ L3 prediction -> the ledger is complete",
      "LEDGER_INCOMPLETE" not in r.failure_modes,
      f"measured/predicted = {r.measured['ledger_ratio']:.2f}x")

# Is the minimum image applied to the step displacement? (traps 1 and 7)
L = 20.0
p0 = np.array([[9.9, 0.0, 0.0]])
p1 = np.array([[-9.9, 0.0, 0.0]])       # wraps across the boundary -- the real
                                        # displacement is 0.2
d = H.measure_step_displacement(p0, p1, L, 2)
check("the minimum image is applied to the step displacement", abs(d - 0.2) < 1e-9,
      f"{d:.4f} (true 0.2; 19.8 if the wrap is ignored)")

print()
print("=" * 78)
print("(2b) the force-based measurement path (drift_direct) -- the hole that left "
      "all 81 runs with 'no measurement'")
print("=" * 78)
# `run.Guard` measures dt*|F|max while health was looking for `step_rms_sigma`, so
# the names did not line up and **this module's central check never ran once**.
r = H.HealthReport()
H.step_health(None, dt_star, dim, predicted_dt_over_tau=2e-3, rep=r, drift_direct=1.8e-3)
check("the force-based path does NOT subtract the thermal part",
      abs(r.measured["dt_over_tau_fast_measured"] - 1.8e-3) < 1e-15,
      f"injected 1.8e-03 -> unchanged at "
      f"{r.measured['dt_over_tau_fast_measured']:.3e} "
      f"(a position difference would have to subtract the thermal part "
      f"{thermal:.2e}, giving 0)")
check("the measurement method is recorded",
      r.measured.get("step_method") == "force",
      f"step_method={r.measured.get('step_method')}")
check("a healthy force-based value gives HEALTHY", r.verdict == "HEALTHY",
      f"modes={r.failure_modes}")

r = H.HealthReport()
H.step_health(None, dt_star, dim, 2e-3, r, drift_direct=5e-2)
check("a too-large force-based step -> NUM_STEP_TOO_COARSE",
      "NUM_STEP_TOO_COARSE" in r.failure_modes, f"modes={r.failure_modes}")

r = H.HealthReport()
H.step_health(None, dt_star, dim, 1e-4, r, drift_direct=8e-3)
check("ledger incompleteness is caught on the force-based path too",
      "LEDGER_INCOMPLETE" in r.failure_modes,
      f"measured/predicted = {r.measured['ledger_ratio']:.0f}x")

# ★ Why the force-based path beats a position difference: **statistical error**, not
#   a floating-point problem. In exact arithmetic a position difference recovers the
#   drift too. The issue is that the rms of a finite sample itself fluctuates by
#   ~rms/sqrt(2N), and when drift << thermal noise that fluctuation swallows the
#   drift. Tested on a real particle ensemble.
rng = np.random.default_rng(20260805)
N_p = 1000
drift_small = 1e-4                                  # 0.5% of the thermal noise (2e-2)
sig = math.sqrt(2.0 * dt_star)                      # thermal std dev per component
noise = rng.normal(0.0, sig, size=(N_p, dim))
step = noise + np.array([drift_small] + [0.0] * (dim - 1))   # deterministic drift
                                                             # along x
rms_sampled = float(np.sqrt((step ** 2).sum(axis=1).mean()))

r_pos = H.HealthReport()
H.step_health(rms_sampled, dt_star, dim, drift_small, r_pos)
r_frc = H.HealthReport()
H.step_health(None, dt_star, dim, drift_small, r_frc, drift_direct=drift_small)
err_pos = abs(r_pos.measured["dt_over_tau_fast_measured"] - drift_small) / drift_small
err_frc = abs(r_frc.measured["dt_over_tau_fast_measured"] - drift_small) / drift_small
# If the drift the position difference recovers is under 10% of the true value, count
# it as "the drift was lost".
# Measured: the sample rms comes out **below** the thermal expectation, so the
# subtraction goes negative and clips to 0.
pos_recovered = r_pos.measured["dt_over_tau_fast_measured"]
check("thermally dominated regime: a finite-sample position difference loses the "
      "drift (force-based is exact)",
      err_frc < 1e-12 and pos_recovered < 0.1 * drift_small,
      f"drift {drift_small:.0e} vs thermal noise {thermal:.1e} "
      f"= {100*drift_small/thermal:.1f}%, N={N_p}: "
      f"position difference recovered {pos_recovered:.2e} "
      f"({100*pos_recovered/drift_small:.0f}% of true -- the sample rms fell below "
      f"the expectation so it clipped to 0) vs force-based error {100*err_frc:.0e}%")

print()
print("=" * 78)
print("(3) the runtime Guard -- does it abort immediately?")
print("=" * 78)

g = H.Guard(box_L=32.0)
try:
    g.check(0, np.zeros((10, 3)), 1.0); ok = True
except RuntimeError:
    ok = False
check("a healthy state passes", ok)

pos = np.zeros((10, 3)); pos[3, 1] = np.nan
try:
    H.Guard(32.0).check(100, pos, 1.0); ok = False
except RuntimeError as e:
    ok = "NUM_NONFINITE" in str(e)
check("NaN in a position -> aborts immediately", ok)

pos = np.zeros((10, 3)); pos[0, 0] = 1e6
try:
    H.Guard(32.0).check(100, pos, 1.0); ok = False
except RuntimeError as e:
    ok = "NUM_DIVERGE" in str(e)
check("runaway coordinates -> aborts immediately", ok)

g = H.Guard(32.0, pe_blowup=100)
g.check(0, np.zeros((5, 3)), 1.0)
try:
    g.check(10, np.zeros((5, 3)), 1e5); ok = False
except RuntimeError as e:
    ok = "NUM_DIVERGE" in str(e)
check("energy blow-up -> aborts immediately", ok)

try:
    H.Guard(32.0).check(0, np.zeros((5, 3)), float("nan")); ok = False
except RuntimeError as e:
    ok = "NUM_NONFINITE" in str(e)
check("NaN in the potential energy -> aborts immediately", ok)

g = H.Guard(32.0, pe_blowup=100)
g.check(0, np.zeros((5, 3)), 1.0)
try:
    g.check(10, np.zeros((5, 3)), 50.0); ok = True     # 50x -- within the limit
except RuntimeError:
    ok = False
check("an increase within the limit passes (not over-sensitive)", ok)

print()
print("=" * 78)
print("(4) the pre-run gate -- against real specs")
print("=" * 78)

specs = sorted((ROOT / "specs").glob("*.json"))
check("specs/ contains specs", len(specs) > 0, f"{len(specs)} found")

spec = ND.load(specs[0])
probs = H.gate(spec)
check(f"a genuine spec passes the gate ({spec.run_id[:34]})", probs == [], f"{probs}")

# A hand-edited spec -- the run_id hash must disagree (rule 2)
raw = json.loads(specs[0].read_text())
key = next(iter(raw["params"]))
raw["params"][key] = raw["params"][key] * 1.5 if isinstance(raw["params"][key], (int, float)) else "X"
tmp = Path(tempfile.mkdtemp()) / "tampered.json"
tmp.write_text(json.dumps(raw))
probs = H.gate(ND.load(tmp))
# ★ Assert on what health.py ACTUALLY emits. This filtered for the Korean
#   "run_id 불일치" while health.py has emitted "run_id mismatch" since it was
#   translated, so the check FAILED while the gate it tests worked correctly --
#   the failure output literally contained the rejection it was looking for.
#   Now matched on the stable `run_id` + `mismatch` tokens rather than a
#   sentence, so rewording the message cannot break it again.
check("hand-edited spec -> the gate rejects it",
      any("run_id" in p and "mismatch" in p for p in probs),
      f"{probs[:1]}")

# Extracting the L3 prediction
p = H.predicted_dt_over_tau(spec)
# ⚠ CURRENTLY FAILING, and correctly so -- it reports a real defect in
#   health.predicted_dt_over_tau, which filters `kind == "integration"` while 274 of
#   278 archived specs still carry the Korean '적분'. Route that call through
#   bdbot.checks.canon_kind and this passes. Not fixed here because
#   bdbot/health.py has uncommitted work from another session.
check("the L3 dt/tau prediction is extracted from the spec",
      p is not None and p > 0, f"dt/tau_fast(L3) = {p}")

n_ok = 0
for sp in specs[:12]:
    s = ND.load(sp)
    if H.predicted_dt_over_tau(s) is not None:
        n_ok += 1
check("the prediction can be extracted for several cases", n_ok >= 8,
      f"{n_ok}/12")

print()
print("=" * 78)
print("(4b) does the gate avoid blocking soft warnings? -- the bug that falsely "
      "rejected 80 of 83")
print("=" * 78)
# ★ A verdict string reads in TWO vocabularies, and both must be accepted.
#   `checks.verdict()` now emits "PASS (N warnings)", but specs/ is a frozen
#   archive written by the older code: measured 2026-08-29, 271 of 278 archived
#   specs carry the Korean "PASS (경고 N건)" and only 1 carries the English form.
#   Testing either spelling alone makes this regression check silently sample
#   almost none of the archive -- which is the failure it exists to prevent.
#   Same read-both-write-one rule as bdbot.checks.canon_kind and runid.py's dual
#   result markers.
def _is_soft_pass(verdict: str) -> bool:
    """Is this a PASS carrying soft warnings, in either vocabulary?"""
    return verdict.startswith("PASS (") and verdict != "PASS"


# The old gate tested `verdict != "PASS"`. `checks.verdict()` returns a
# "PASS (N warnings)" form whenever a soft check trips, so **every spec with only
# warnings was rejected**. bd-physics section 4 defines statistics and finite-size
# as warnings, not failures, and `run.execute` read startswith("FAIL") correctly --
# the reason the two disagreed unnoticed is that `execute` never called `gate()`.
blocked_wrongly = []
n_soft = n_pass = 0
for sp in specs:
    s = ND.load(sp)
    probs = H.gate(s)
    hard = [c for c in s.checks if getattr(c, "hard", True) and not H._ok(c)]
    errs = [i for i in s.raw.get("l3_issues", []) if i.get("level") == "error"]
    if _is_soft_pass(s.verdict):
        n_soft += 1
    if not probs:
        n_pass += 1
    if probs and not hard and not errs and not s.verdict.startswith("FAIL"):
        blocked_wrongly.append((s.run_id, s.verdict))
check("a spec with only soft warnings is not blocked", not blocked_wrongly,
      f"of {len(specs)} specs: {n_soft} warnings-only, {n_pass} passed the gate, "
      f"{len(blocked_wrongly)} falsely rejected"
      + (f"  e.g. {blocked_wrongly[0]}" if blocked_wrongly else ""))

soft_specs = [sp for sp in specs if _is_soft_pass(ND.load(sp).verdict)]
if soft_specs:
    s = ND.load(soft_specs[0])
    check("warnings do not block but MUST be surfaced via gate_notes",
          len(H.gate_notes(s)) > 0,
          f"{s.run_id[:34]}: {len(H.gate_notes(s))} note(s) -- "
          f"{H.gate_notes(s)[0][:60]}")
else:
    check("a spec with warnings must exist for this regression to be testable",
          False, "none found")

# A hard failure must still block (guards the opposite malfunction)
raw = json.loads(specs[0].read_text())
raw["verdict"] = "FAIL"
tmp2 = Path(tempfile.mkdtemp()) / "failverdict.json"
tmp2.write_text(json.dumps(raw))
probs = H.gate(ND.load(tmp2))
check("verdict=FAIL still blocks", any("verdict=FAIL" in p for p in probs),
      f"{probs}")

# An L3 integrity **error** blocks; warn and info do not
for lvl, should_block in (("error", True), ("warn", False), ("info", False)):
    raw = json.loads(specs[0].read_text())
    raw["l3_issues"] = [{"level": lvl, "where": "test", "msg": "injected"}]
    t3 = Path(tempfile.mkdtemp()) / f"iss_{lvl}.json"
    t3.write_text(json.dumps(raw))
    # ★ Filter on what health.py ACTUALLY emits. This looked for the Korean
    #   "무결성" while health.py has emitted "L3 integrity errors (...)" since
    #   it was translated, so `probs` was always empty and the error case
    #   reported "does not block" when it does. Matched on "integrity", the
    #   stable token, rather than the sentence.
    probs = [p for p in H.gate(ND.load(t3)) if "integrity" in p.lower()]
    check(f"l3_issues level={lvl} -> "
          f"{'blocks' if should_block else 'does not block'}",
          bool(probs) == should_block, f"{probs}")

print()
print("=" * 78)
print(f"{len(PASS)}/{len(PASS) + len(FAIL)} PASS")
if FAIL:
    for f in FAIL:
        print(f"  ✗ {f}")
print("=" * 78)
sys.exit(1 if FAIL else 0)
