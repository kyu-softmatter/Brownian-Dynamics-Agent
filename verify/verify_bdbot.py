"""Self-verification of `bdbot/` -- do the shared modules reproduce the case numbers?

Checks that what was factored out in Phase 1-C does not disagree with the measured
values from 1-A and 1-B.
Every expected value here is **a number already confirmed by execution**
(skill `bd-physics` §6.1, §6.2).

    $PY scratch/verify_bdbot.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from bdbot import Q, checks as C, materials as M, scales as SC, sim as SIM, stats as ST

fails = []


def check(label, got, want, rtol=1e-4, unit=""):
    ok = abs(got - want) <= rtol * abs(want)
    print(f"  {'✓' if ok else '✗'} {label:<42} {got:>14.6g} {unit:<8} (expected {want:.6g})")
    if not ok:
        fails.append(label)


print("=" * 84)
print("(1) material property derivation -- reproducing the 1-A/1-B measured values  "
      "(d=5 µm silica, water@300 K)")
print("=" * 84)
b = M.sphere_bulk(Q(5.0, "um"), Q(300, "K"), Q(0.851, "mPa*s"), Q(2000, "kg/m^3"))
check("kT", float(b["kT"].to("J").m), 4.1419e-21, 1e-4, "J")
check("γ = 3πηd", float(b["gamma"].to("kg/s").m), 4.0102e-8, 1e-4, "kg/s")
check("D_t = kT/γ", float(b["D_t"].to("um^2/s").m), 0.1033, 1e-3, "µm²/s")
check("τ_B = d²/D_t", float(b["tau_B"].to("s").m), 242.05, 1e-4, "s")
check("τ_p = m/γ", float(b["tau_p"].to("us").m), 3.264, 1e-3, "µs")

print("\n(2) the dt convention -- do the trap and the soft pair share one formula? "
      "(the key finding of 1-C)")
k_trap = Q(10, "pN/um").to("N/m")
tau_k = C.relaxation_time(b["gamma"], k_trap)
check("τ_k = γ/k  (1-A)", float(tau_k.to("ms").m), 4.010, 1e-3, "ms")
check("dt(bias 0.1%) = 2*bias*tau_k", float(C.dt_from_bias(tau_k, 1e-3).to("us").m),
      8.020, 1e-3, "µs")
check("inverting the bias, bias_from_dt",
      C.bias_from_dt(C.dt_from_bias(tau_k, 1e-3), tau_k),
      0.1, 1e-6, "%")
# The soft pair uses the same formula: give the stiffness in kT/d^2 units and
# tau_int = tau_B/U''*
# ★ Per-case comparison against measured values is done by snapshot in
#   scratch/verify_1c_equivalence.py. Here only **the identity of the formula** is
#   checked (using the rounded values printed in a report disagrees by 0.9% from
#   3-digit rounding alone, because U'' goes as r^-5 -- this actually happened once).
Upp_star = 276.4                      # the value the code used for 1-B at A=100
tau_int = b["tau_B"] / Upp_star
check("dt = 10⁻²·τ_int → dt/τ_B", float((C.dt_from_gate(tau_int) / b["tau_B"]).to("").m),
      1e-2 / Upp_star, 1e-12)
check("do the two cases share one formula: gamma/k == gamma/U''",
      float((C.relaxation_time(b["gamma"], Q(1, "N/m")) ).to("s").m),
      float((b["gamma"] / Q(1, "N/m")).to("s").m), 1e-12, "s")

print("\n(3) checks -- hard/soft classification and the verdict "
      "(a distinction 1-B made necessary)")
ck = [C.Check("model", "inertia", 8.1e-4, C.GATE),
      C.Check("integration", "resolution", C.GATE, C.GATE),
      C.Check("geometry", "cutoff", 0.5, 1.0),
      C.Check("statistics", "observation window", 50.0, 100.0, ">=", hard=False)]
v, hf, sf, tight = C.verdict(ck)
# ★ Assert the PROPOSITION, not the wording. This line used to string-match the
#   rendered verdict, and that broke silently on 2026-08-28 (c0074a2) when
#   bdbot/checks.py was translated and began emitting "PASS (1 warnings)" instead
#   of the Korean it was compared against. The script then FAILed for a day
#   without anyone noticing, because it is not part of the pytest suite. What this
#   check actually means is "no hard failure, exactly one soft failure, and the
#   verdict still passes" -- so that is what it now asserts, and rewording the
#   message can no longer break it.
verdict_ok = v.startswith("PASS") and not hf and len(sf) == 1
print(f"  {'✓' if verdict_ok else '✗'} verdict = {v!r}  "
      f"(hard failures {len(hf)} / soft failures {len(sf)} / thin margins {len(tight)})")
if not verdict_ok:
    fails.append("verdict")
ck2 = [C.Check("geometry", "cutoff", 1.5, 1.0)]
v2, _, _, _ = C.verdict(ck2)
print(f"  {'✓' if v2 == 'FAIL' else '✗'} with a hard violation present: {v2!r}")
if v2 != "FAIL":
    fails.append("verdict-hard")

print("\n(4) statistics -- does the autocorrelation correction actually work?")
rng = np.random.default_rng(0)
white = rng.normal(0, 1, 4000)
walk = np.cumsum(rng.normal(0, 1, 4000))
nw, nk = ST.n_eff(white), ST.n_eff(walk)
print(f"  {'✓' if nw > 3000 else '✗'} white noise   n_eff = {nw:8.1f} / 4000  "
      f"(almost all independent)")
print(f"  {'✓' if nk < 500 else '✗'} random walk   n_eff = {nk:8.1f} / 4000  "
      f"(strongly correlated)")
if not (nw > 3000 and nk < 500):
    fails.append("n_eff")
# Drift: does it report significance and magnitude separately?
trend = np.linspace(0, 0.001, 4000) + rng.normal(0, 1e-5, 4000) + 105.5
st = ST.stationarity(trend)
print(f"  {'✓' if abs(st['drift_span_rel_pct']) < 0.01 else '✗'} negligible drift: "
      f"t={st['trend_t']:+.1f} yet only {st['drift_span_rel_pct']:+.5f}% across the "
      f"whole span -> you have to look at the magnitude too, or you get a false alarm")
if abs(st["drift_span_rel_pct"]) >= 0.01:
    fails.append("drift")
# Unbiased autocorrelation: does it recover the tau of an OU process?
tau_true = 40.0
n = 200000
x = np.empty(n)
x[0] = 0.0
a = math.exp(-1 / tau_true)
s = math.sqrt(1 - a * a)
g = rng.normal(0, 1, n)
for i in range(1, n):
    x[i] = a * x[i - 1] + s * g[i]
ac = ST.autocorr_unbiased(x[:, None])
tau_fit = -1.0 / math.log(ac[1] / ac[0])
check("recovering tau from the unbiased autocorrelation (OU)",
      tau_fit, tau_true, 2e-2)

print("\n(5) sim -- are the trap guards actually baked into the code?")
ns, hs = SIM.resolve_seed(20260803)
print(f"  {'✓' if hs == 10179 else '✗'} trap 12 seed: {ns} -> HOOMD {hs} "
      f"(the warning value 10179)")
if hs != 10179:
    fails.append("seed")
d = SIM.minimum_image(np.array([[9.0, -9.0, 0.0]]), 10.0, dims=2)
ok = np.allclose(d, [[-1.0, 1.0, 0.0]]) and np.isfinite(d).all()
print(f"  {'✓' if ok else '✗'} traps 1 and 7, minimum image: "
      f"[9,-9,0] -> {d[0].tolist()} (z is not NaN)")
if not ok:
    fails.append("minimum_image")

print("\n(6) the ledger -- can it derive by ratio and show the span in decades?")
lg = SC.ScaleLedger()
lg.lengths = {"d": Q(5, "um").to("m"), "l_k": Q(20.35, "nm").to("m"), "L": Q(150, "um").to("m")}
lg.times = {"tau_p": b["tau_p"], "tau_k": tau_k, "tau_B": b["tau_B"]}
lg.ref = SC.thermal_reference(b["d"], b["kT"], b["tau_B"])
order = [k for k, _ in lg.sorted_items(lg.times)]
print(f"  {'✓' if order == ['tau_p', 'tau_k', 'tau_B'] else '✗'} "
      f"timescales, sorted: {order}")
if order != ["tau_p", "tau_k", "tau_B"]:
    fails.append("sorted")
check("ratio(times, tau_p, tau_k)", lg.ratio("times", "tau_p", "tau_k"), 8.139e-4, 1e-3)
check("span(times) in decades", math.log10(lg.span("times")), 7.87, 1e-2, "decades")

print()
print("=" * 84)
print("✓ PASS -- the bdbot shared modules reproduce the cases' measured values"
      if not fails else f"✗ FAIL -- {len(fails)}: {fails}")
print("=" * 84)
sys.exit(0 if not fails else 1)
