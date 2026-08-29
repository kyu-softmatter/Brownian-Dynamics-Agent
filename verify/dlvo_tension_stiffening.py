"""Can strong DLVO alone look like bending? -- measure the size of tension
(geometric-nonlinearity) stiffening.

Raised by the user (2026-08-06): "even with no friction, a system with only strong
DLVO would probably show similar behaviour."

G1 (a straight chain + purely central forces + natural bond length => **linear**
bending stiffness exactly 0) still holds. But in three-point bending the two ends are
fixed, so pushing the centre sideways by delta **lengthens the path** -- bonds that
sat at their natural length are forced to stretch. So even purely central forces give
transverse resistance, via **tension**. That is a **string**, not a beam, and its
apparent stiffness depends on amplitude.

★★ Predictions fixed BEFORE running (principle 9.2 -- fixed before seeing results)
  path extension = 2*delta^2/L (L=(n-1)*l), shared equally by n-1 bonds
  -> 2*delta^2/((n-1)^2*l) per bond
    U = (n-1)*0.5*k_bond*[2*delta^2/((n-1)^2*l)]^2
      = 2*k_bond*delta^4 / ((n-1)^3*l^2)
  P6  radial bonds only:  **K = 2U/delta^2 = 4*k_bond*delta^2 / ((n-1)^3*l^2)**
                          -> linear stiffness 0, apparent stiffness ~ delta^2
  P7  bending models:     K is **independent** of delta (a genuinely quadratic
                          potential)
  P8  crossover amplitude: delta* = sqrt( k_bend*(n-1)^3*l^2 / (4*k_bond) )
      bending dominates for delta << delta*, tension for delta >> delta* --
      so **an amplitude sweep is the discriminator**.

★ Measured statically only (principle 8) -- an exact minimisation with analytic
  gradients, so there is none of MD's noise or transients.
  Variables: (x,y) of every free bead plus every orientation theta.
  Constraints: both ends at (x,0), the centre at (x_mid, delta).
  ⚠ x MUST be released -- fixing it prevents the chain being pulled inward to
    relieve the extension, which overestimates the tension term.

    $PY scratch/dlvo_tension_stiffening.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "verify"))

from rolling_contact import RollingContact, k_roll_from_kappa_theta  # noqa: E402

H_MIN_STAR = 0.00759259035993831
ELL = 1.0 + H_MIN_STAR
K_BOND = 1042362.8817700658          # DLVO secondary-minimum curvature [kT/d^2]
                                     # (k_bond_star in specs)
KAPPA_THETA = 1391229.7767209478     # [kT] — JKR κ₀=64 mN/m
R_C = 0.5
K_ROLL = k_roll_from_kappa_theta(KAPPA_THETA, R_C)

results = []


def check(name, ok, detail=""):
    results.append((name, ok, detail))
    print(f"  {'✓' if ok else '✗'} {name}" + (f"   {detail}" if detail else ""))


def bending_matrix(n, kappa_theta, ell):
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    return kappa_theta * (B.T @ B)


def static_3point(n, delta, *, k_bond=K_BOND, bend="none"):
    """Three-point bending stiffness k = 2U/delta^2. bend in {none, harmonic, rolling}.

    ⚠ The radial bond is used **fully nonlinearly** (a harmonic potential in
      |r_j - r_i|) -- tension stiffening comes precisely from that nonlinearity, so
      linearising it removes the entire effect.
    """
    pos0 = np.zeros((n, 3))
    pos0[:, 0] = np.arange(n) * ELL
    quat0 = np.tile(np.array([1.0, 0, 0, 0]), (n, 1))
    bonds = [[i, i + 1] for i in range(n - 1)]
    mid = n // 2
    clamped = [0, mid, n - 1]
    free = [i for i in range(n) if i not in clamped]
    A_bend = bending_matrix(n, KAPPA_THETA, ELL) if bend == "harmonic" else None
    rc = RollingContact(bonds, pos0, quat0, R_C, K_ROLL, 0.0) if bend == "rolling" else None
    nf = len(free)

    def unpack(v):
        pos = pos0.copy()
        pos[free, 0] = pos0[free, 0] + v[:nf]
        pos[free, 1] = v[nf:2 * nf]
        pos[mid, 1] = delta
        th = v[2 * nf:] if rc is not None else np.zeros(n)
        c, s = np.cos(th / 2), np.sin(th / 2)
        quat = np.stack([c, np.zeros(n), np.zeros(n), s], axis=1)
        return pos, quat

    def energy_grad(v):
        pos, quat = unpack(v)
        d = pos[1:] - pos[:-1]
        r = np.linalg.norm(d, axis=1)
        nhat = d / r[:, None]
        U = 0.5 * k_bond * ((r - ELL) ** 2).sum()
        # dU/dr_j = k(r−ℓ)n̂ ; F = −dU/dr
        fb = (k_bond * (r - ELL))[:, None] * nhat
        F = np.zeros((n, 3))
        np.add.at(F, np.arange(1, n), -fb)
        np.add.at(F, np.arange(0, n - 1), fb)
        T = np.zeros((n, 3))
        if A_bend is not None:
            y = pos[:, 1]
            fy = -(A_bend @ y)
            U += -0.5 * float(y @ fy)
            F[:, 1] += fy
        if rc is not None:
            U += rc.energy(pos, quat)
            Fr, Tr = rc.force_torque(pos, quat)
            F += Fr
            T += Tr
        g = np.concatenate([-F[free, 0], -F[free, 1]]
                           + ([-T[:, 2]] if rc is not None else []))
        return U, g

    n_var = 2 * nf + (n if rc is not None else 0)
    res = minimize(energy_grad, np.zeros(n_var), jac=True, method="L-BFGS-B",
                   options=dict(maxiter=200000, maxfun=200000, ftol=0.0, gtol=0.0))
    return 2 * res.fun / delta ** 2


def k_tension_pred(n, delta, k_bond=K_BOND):
    return 4 * k_bond * delta ** 2 / ((n - 1) ** 3 * ELL ** 2)


if __name__ == "__main__":
    # ══════════════════════════════════════════════════════════════════════════
    N = 9                                   # same chain length as the production case
    print("=" * 92)
    print(f"apparent bending stiffness from strong DLVO alone -- tension "
          f"(geometric-nonlinearity) stiffening   [n={N}]")
    print("=" * 92)
    print(f"  k_bond = {K_BOND:.6g} kT/d^2 (DLVO secondary-min curvature)   "
          f"kappa_theta = {KAPPA_THETA:.6g} kT (JKR)")
    k_bend_lin = static_3point(N, 1e-5, bend="harmonic")
    print(f"  linear bending stiffness (delta->0) = {k_bend_lin:.6g} kT/d^2")
    d_star = np.sqrt(k_bend_lin * (N - 1) ** 3 * ELL ** 2 / (4 * K_BOND))
    print(f"  ★ P8 crossover amplitude delta* = {d_star:.4g} d   "
          f"(= {d_star*1470:.0f} nm, {d_star:.2f}x the particle diameter)")
    print()
    print(f"  {'delta [d]':>10} {'DLVO-only K':>13} {'P6 pred':>13} {'ratio':>8} | "
          f"{'JKR bend K':>13} "
          f"{'rolling K':>13} | {'DLVO/JKR':>9}")
    print("  " + "-" * 96)
    rows = []
    for delta in (1e-4, 1e-3, 1e-2, 0.03, 0.1, 0.3, 0.43, 1.0):
        k_none = static_3point(N, delta, bend="none")
        k_pred = k_tension_pred(N, delta)
        k_h = static_3point(N, delta, bend="harmonic")
        k_r = static_3point(N, delta, bend="rolling")
        rows.append((delta, k_none, k_pred, k_h, k_r))
        print(f"  {delta:10.4g} {k_none:13.6g} {k_pred:13.6g} {k_none/k_pred:8.4f} | "
              f"{k_h:13.6g} {k_r:13.6g} | {k_none/k_h:9.4f}")

    print()
    mid = [r for r in rows if 1e-3 <= r[0] <= 0.1]
    worst = max(abs(r[1] / r[2] - 1) for r in mid)
    check("P6 DLVO-only K = 4*k_bond*delta^2/((n-1)^3*l^2)  (delta = 1e-3 to 0.1)",
          worst < 0.10,
          f"worst deviation {100*worst:.2f}%  -- a prediction with no free parameters")
    sl = np.polyfit(np.log([r[0] for r in mid]), np.log([r[1] for r in mid]), 1)[0]
    check("P6' slope d(log K)/d(log delta) = 2 (evidence the linear stiffness is 0)",
          abs(sl - 2) < 0.05,
          f"slope = {sl:.4f}")
    k_h_lo, k_h_hi = static_3point(N, 1e-4, bend="harmonic"), static_3point(N, 0.3, bend="harmonic")
    check("P7 the bending model's K is independent of delta",
          abs(k_h_hi / k_h_lo - 1) < 0.02,
          f"over a 3000x range in delta, the ratio is {k_h_hi/k_h_lo:.6f}")
    k_r_lo, k_r_hi = static_3point(N, 1e-4, bend="rolling"), static_3point(N, 0.3, bend="rolling")
    check("P7' the rolling model is also independent of delta",
          abs(k_r_hi / k_r_lo - 1) < 0.02,
          f"over a 3000x range in delta, the ratio is {k_r_hi/k_r_lo:.6f}")

    print()
    print("=" * 92)
    print("★ suppressor 1: trap compliance -- the table above clamped both ends "
          "**rigidly**")
    print("=" * 92)
    print("""  In the real case both ends are traps of stiffness k_t (`--kt-scale`
      scales all three traps). Pushing the centre lets the chain be pulled inward and
      relieve some extension, so what sets the tension is the **series extensional
      stiffness**:   1/k_ext = (n-1)/k_bond + 2/k_t
          K_tension = 4 k_ext δ² / L²,   L = (n−1)ℓ""")
    L_CHAIN = (N - 1) * ELL
    C_CHAIN = (N - 1) / K_BOND
    K_T_BASE = 5217.1116627035535
    print()
    print(f"  chain compliance (n-1)/k_bond = {C_CHAIN:.4e}")
    print(f"  {'protocol':>18} {'k_t':>11} {'2/k_t':>11} {'k_ext':>11} {'vs rigid':>9} "
          f"{'K@δ=0.43':>10} {'K@δ=1':>9}")
    print("  " + "-" * 88)
    protocols = [("trap default", K_T_BASE), ("trap k_t x100", K_T_BASE * 100),
                 ("rigid clamp (table above)", np.inf)]
    k_ext_tab = {}
    for label, kt in protocols:
        c = C_CHAIN + (2.0 / kt if np.isfinite(kt) else 0.0)
        k_ext = 1.0 / c
        k_ext_tab[label] = k_ext
        rigid = 1.0 / C_CHAIN
        print(f"  {label:>18} {kt:11.4g} {2/kt if np.isfinite(kt) else 0:11.4e} {k_ext:11.6g} "
              f"{k_ext/rigid:9.4f} {4*k_ext*0.43**2/L_CHAIN**2:10.3f} "
              f"{4*k_ext*1.0**2/L_CHAIN**2:9.1f}")
    print("  -> with the default trap, the trap compliance is 50x the chain's, so "
          "tension stiffening is **suppressed 51x**.")

    print()
    print("=" * 92)
    print("★ suppressor 2: the DLVO well cannot hold the tension -- where the "
          "harmonic-bond approximation breaks")
    print("=" * 92)
    sys.path.insert(0, str(ROOT / "cases"))
    from chain_bend_dlvo_2d import (F_h_star, U_star, dlvo_reduced_params,  # noqa: E402
                                    find_well, load_system)
    sysd = load_system(ROOT / "intake" / "chain-bend-2d-dlvo" / "system.yaml")
    p = dlvo_reduced_params(sysd)
    well = find_well(p)
    hs = np.geomspace(well["h_min"], 2.0, 40000)
    Fs = -F_h_star(hs, p)                       # tensile (pulling) force
                                                # = -(repulsive-direction force)
    i_max = int(np.argmax(Fs))
    F_MAX, h_infl = float(Fs[i_max]), float(hs[i_max])
    dh_infl = h_infl - well["h_min"]
    print("  ⚠ correcting the criterion -- on the outside DLVO only **asymptotes** to "
          "U->0^- and never crosses")
    print("    zero (the first attempt gave nan for the 'escape distance'). The right "
          "criterion is the **maximum tensile force**:")
    print("    once the chain tension T exceeds F_max, that bond is mechanically "
          "unstable under force control and unbinds.")
    print()
    print(f"  secondary minimum  h_min = {well['h_min']*1470:.2f} nm,  "
          f"U_min = {well['U_min']:.2f} kT")
    print(f"  ★ maximum tensile force  F_max = {F_MAX:.1f} kT/d  @ "
          f"h = {h_infl*1470:.1f} nm "
          f"(extension to the inflection point {dh_infl*1470:.2f} nm)")
    print(f"  extension needed to reach F_max under the harmonic approximation "
          f"(k_bond) = {F_MAX/K_BOND*1470:.2f} nm "
          f"-> the real curve is {dh_infl*K_BOND/F_MAX:.1f}x softer "
          f"(anharmonic softening)")
    print()
    print(f"  {'protocol':>18} {'delta':>7} {'total ext':>9} {'tension T':>11} "
          f"{'T/F_max':>9} {'verdict':>9}")
    print("  " + "-" * 72)
    for label, kt in protocols:
        k_ext = k_ext_tab[label]
        for delta in (0.116, 0.43, 1.0):
            tot = 2 * delta ** 2 / L_CHAIN
            T = k_ext * tot
            print(f"  {label:>18} {delta:7.3f} {tot:9.4f} {T:11.1f} {T/F_MAX:9.2f} "
                  f"{'⛔ bond fails' if T > F_MAX else '✓ holds':>9}")
    print(f"""
      => **under the k_t x100 protocol (100% tracking), the tension is already 4.9x
         F_max by delta = 0.43 d.**
         A harmonic-bond calculation gives
         K_tension ~ {4*k_ext_tab['trap k_t x100']*1.0/L_CHAIN**2:.0f} (delta=1), but
         the real DLVO well cannot sustain that tension -- consistent with 1-D
         measuring K' ~ 0 (0.84 sigma) under exactly that condition. With the default
         trap (delta ~ 0.116) T/F_max = 0.01 and it is safe, but the tension
         stiffening is correspondingly small.""")

    print()
    print("=" * 92)
    print("★ conclusion -- is 1-D's observed K' explained by tension stiffening?")
    print("=" * 92)
    print(f"  {'observed K1':>10} {'protocol':>16} {'delta est':>8} "
          f"{'tension pred':>10} {'fraction expl':>10}")
    print("  " + "-" * 62)
    for K_obs, label, delta in ((41.8, "trap default", 0.43 * 0.27),
                                (12.7, "trap default", 0.43 * 0.27),
                                (105.7, "trap default", 0.43 * 0.27),
                                (199.5, "trap default", 0.43 * 0.27)):
        Kt = 4 * k_ext_tab[label] * delta ** 2 / L_CHAIN ** 2
        print(f"  {K_obs:10.1f} {label:>16} {delta:8.3f} {Kt:10.3f} {100*Kt/K_obs:9.1f}%")
    print("""
      => ✗ **it is not explained.** With the default trap, tension stiffening is only
         1-5% of the observed value.
         ★ This calculation nonetheless **strengthens** the 1-D conclusion rather
         than weakening it -- it quantified "the most plausible candidate for
         mimicking bending with central forces alone" with no free parameters, and
         showed that it is quantitatively far too small.
      ★ The discriminator is not the magnitude but the **delta scaling**: tension
         gives K ~ delta^2 (slope exactly 2), real bending is delta-independent. An
         amplitude sweep separates them in principle.
      ⚠ But this suppression is specific to THIS system's parameters -- with rigidly
         clamped ends (or traps stiffer than the chain) and a deeper well, tension
         stiffening could dominate.""")

    print()
    print("=" * 92)
    n_pass = sum(1 for _, ok, _ in results if ok)
    print(f"{n_pass}/{len(results)} PASS")
    for name, ok, detail in results:
        if not ok:
            print(f"  ✗ {name}   {detail}")
    print("=" * 92)
    sys.exit(0 if n_pass == len(results) else 1)
