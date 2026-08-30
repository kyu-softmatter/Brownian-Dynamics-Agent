"""System-level rheology -- from **the whole system's energy budget** rather than a
per-bead fit.

★ Why this viewpoint
  Until now K*(omega) came from the phasor of a single driven bead
  (`lockin.k_star`).
  Here the whole system is treated as one viscoelastic body, looking only at **the
  work done per cycle and the energy stored**. It needs to know neither how many beads
  there are nor which bead has which phase.

  (1) total dissipation (model-independent) -- in steady state, all the work the
      driving trap does leaves into the fluid
        W_cycle = ∮ F·dy,   F = k_t(y_ghost − y_bead)
        K″_total ≡ W_cycle / (π|ŷ|²)
  (2) decomposition -- that dissipation arises in two places
        K″_total = K″_chain + ωγ
        (chain dissipation)  (**the driven bead's own solvent drag**, present even
                              with no chain)
     `lockin.k_star` subtracts -i*omega*gamma by definition and returns K''_chain
     alone (see its docstring).
     -> the two routes must differ by **exactly omega*gamma**. That is a
     cross-check with no free parameters.
     Measured: JKR -0.71%, DLVO -0.11% -- they agree.
  (3) total storage -- the **2*omega** component of the system's total potential
      energy U(t) (U ~ y^2, hence the second harmonic)
        amplitude of U_osc = 0.25*K'_sys*|y_hat|^2
        ->  K'_sys = 4|U_hat_2omega|/|y_hat|^2
     ⚠ K'_sys is NOT K'_lockin. The lock-in value is the stiffness felt **at the
       drive point**, while K'_sys is the energy stored across the whole chain
       (bending, bond stretching, traps). The difference IS the "fraction stored away
       from the drive point", so it is itself information.

★★ The system-level conclusion in one line:
   the DLVO chain is **rheologically invisible** -- the whole system's dissipation
   equals that of a single driven bead being dragged freely. For JKR the chain carries
   76% of the dissipation.

    $PY scratch/system_moduli.py
"""
import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
kT = 1.380649e-23 * 300.0
d = 1.47e-6
KTD2 = kT / d ** 2                      # kT/d² → N/m


def analyze(pat, label):
    rows = []
    for dd in sorted(glob.glob(str(ROOT / pat))):
        if not (Path(dd) / "metrics.json").exists():
            continue
        z = np.load(Path(dd) / "observables.npz", allow_pickle=True)
        s = json.loads((Path(dd) / "spec.json").read_text())
        om = float(s["params"]["omega_star"])
        n = int(s["params"]["n_beads"])
        k_t = float(s["params"]["k_t_star"])
        pos = str(s["params"].get("drive_mode", "trap")) == "position"
        t, pe, yb, yg = (np.asarray(z[k], float) for k in ("t", "pe", "y_bead", "y_ghost"))
        U = pe * n                                       # pe is per particle -> the
                                                         # whole system
        m = {x["name"]: x["measured"] for x in
             json.loads((Path(dd) / "metrics.json").read_text())["observables"]}

        Ay = abs(2 * np.mean(yb * np.exp(-1j * om * t)))          # driven-bead
                                                                  # amplitude |y_hat|

        # (1) total dissipation -- the work the trap does on the driven bead. In
        #     position mode there is no trap, so it is undefined.
        if pos:
            K2_tot = np.nan
        else:
            F = k_t * (yg - yb)
            Pmean = np.trapezoid(F * np.gradient(yb, t), t) / (t[-1] - t[0])   # mean power
            K2_tot = Pmean * (2 * np.pi / om) / (np.pi * Ay ** 2)

        # (3) total storage -- the 2*omega component of U.
        #     U_osc = -0.25*K'|y_hat|^2 cos(2*omega*t)
        #     -> |U_hat_2omega| = 0.25*K'|y_hat|^2
        Uhat2 = 2 * np.mean((U - U.mean()) * np.exp(-2j * om * t))
        Kp_sys = 4 * abs(Uhat2) / Ay ** 2

        rows.append(dict(seed=s["numerics"]["seed"], Ay=Ay, om=om, n=n, pos=pos,
                         Kp_lock=m.get("K_prime"), K2_lock=m.get("K_doubleprime"),
                         K2_tot=K2_tot, Kp_sys=Kp_sys))
    if not rows:
        print(f"  [{label}] no runs")
        return None
    g = lambda k: (np.nanmean([r[k] for r in rows]),
                   np.nanstd([r[k] for r in rows], ddof=1) / np.sqrt(len(rows)))
    om = rows[0]["om"]
    print(f"\n  {label}   {len(rows)} seeds   |y_hat| = {g('Ay')[0]:.5f} d")
    k2t, k2te = g("K2_tot")
    if np.isfinite(k2t):
        k2l, k2le = g("K2_lock")
        drag = om                                        # omega*gamma, with gamma*=1
                                                         # in reduced units
        print(f"    (1) total dissipation  K''_total = {k2t:>10.4g} ± {k2te:<8.3g}")
        print(f"    (2) decomposition   chain {k2l:>10.4g}  +  driven-bead drag "
              f"omega*gamma = {drag:.4g}"
              f"  =  {k2l + drag:>10.4g}   "
              f"({100*(k2t-k2l-drag)/k2t:+.2f}% vs the total)")
        print(f"       -> fraction of dissipation carried by the chain = "
              f"{100*k2l/k2t:>6.1f}%")
    kps, kpse = g("Kp_sys")
    kpl, _ = g("Kp_lock")
    print(f"    (3) total storage   K'_sys = {kps:>10.4g} ± {kpse:<8.3g}"
          + (f"   (drive-point K'_lockin {kpl:.4g}, ratio {kps/kpl:.2f}x)"
             if kpl and abs(kpl) > 1 else
             f"   (drive-point K'_lockin = {kpl:.3g} ~ 0, yet storage exists "
             f"-- bond stretching)"))
    return rows


print("=" * 98)
print("system-level rheology -- the whole system's energy budget "
      "(no per-bead fitting)")
print("=" * 98)
P = "runs/chain-bend-2d-dlvo__n9-w3000-a1470"
rj = analyze(f"{P}-jkr-kt100__*", "DLVO+JKR  (k_t×100)")
rd = analyze(f"{P}-kt100__*", "DLVO-only (k_t×100)")

if rj and rd:
    om = rj[0]["om"]
    tj = np.nanmean([r["K2_tot"] for r in rj])
    td = np.nanmean([r["K2_tot"] for r in rd])
    print()
    print("=" * 98)
    print("★★ the system-level conclusion")
    print("=" * 98)
    print(f"  dissipation of the driven bead **alone** (no chain) : "
          f"omega*gamma = {om:.4g}")
    print(f"  whole-system dissipation with the DLVO chain attached : {td:.4g}"
          f"   -> {td/om:.3f}x  (chain contributes {100*(td-om)/om:+.1f}%)")
    print(f"  whole-system dissipation with the JKR  chain attached : {tj:.4g}"
          f"   -> {tj/om:.3f}x  (chain contributes {100*(tj-om)/om:+.1f}%)")
    print()
    print("  => the DLVO chain is **rheologically invisible**: the dissipation of")
    print("     the whole system equals shaking one bead in water, so attaching the")
    print("     chain changes nothing. For JKR the chain dominates the response.")
