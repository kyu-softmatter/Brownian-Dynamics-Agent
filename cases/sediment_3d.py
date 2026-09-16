"""`sediment-pmma-3d` — sedimentation equilibrium of a dilute hard-sphere colloid.

    PY=./bin/py
    $PY cases/sediment_3d.py --report                 # the dimensionless report, no run
    $PY cases/sediment_3d.py --spec                   # write specs/<run_id>.json, no run
    $PY cases/sediment_3d.py --seed 1 --require-seal  # run one seed, gated on the seal
    $PY cases/sediment_3d.py --lxy 24 --seed 1        # the finite-size arm

Source: [Newman & Yethiraj 2014](../knowledge/source/papers/2014-newman-yethiraj-sedimentation-equilibrium.md)
(`arXiv:1412.3190v3`). Card:
[`passive-sphere--sedimentation`](../knowledge/wiki/systems/passive-sphere--sedimentation.md).

## What this case is for

The profile `Phi(z)` is a Boltzmann distribution in a gravitational potential, so
**one run sweeps the equation of state over a whole range of local volume
fraction** — `Phi` runs from about 0.159 at the wall down to 1e-4 at the top of
the cell, where Carnahan-Starling goes from 98 % above the ideal gas down to 1.
That is the measurement; the fitted `l_g` is the correctness check that licenses
it.

★ **The EOS is read from the PROFILE, never from HOOMD's pressure.**
`knowledge/wiki/findings/wca-reproduces-carnahan-starling.md` trap 2: HOOMD's
pressure carries a kinetic term built from velocities, and overdamped BD
velocities are not physical (measured `kT_kin` = 1.0009-1.015). The hydrostatic
route `dP/dz = -rho kT/l_g` never touches them.

## Conventions that silently bite, both measured on 2026-09-16

★ **`sigma_LJ = 2^(-1/6) d`, `r_cut = d`.** `md.pair.LJ(sigma=1.0,
r_cut=2**(1/6))` puts the WCA minimum at 1.1225 d, so the effective particle
diameter is 1.1225 and `phi` comes out `1.1225^3 = 1.414x` too large — **+41.4 %,
with no symptom**, because every other number still looks reasonable.
`simbot/cutoff.py`'s `hoomd_note` states the right convention.

★ **`phi` must be put into Carnahan-Starling through the Barker-Henderson
diameter.** At `eps = 1 kT` and this `sigma_LJ`, `d_BH = 0.9048 d`, so
`phi_eff/phi = 0.7407` and `CS(phi_eff)` differs from `CS(phi_nominal)` by
**-17.05 %** at the wall. Both are reported, because which one the measurement
agrees with is the answer rather than a detail.

Three more gates the card carries, all from
`verify/verify_sedimentation_wall.py`: the **wall** binds `dt` (an LJ wall reaches
1517 kT/d and a particle passes straight through with no NaN), the **cell height**
sets how empty the top is, and equilibration is the **drift** time
`L_z l_g / D_0`, not `l_g^2/D_0` — 1672 vs 179 `tau_d` here.
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from bdbot import Q, checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import dt as DT, nondim as ND, placement as PL, run as RUN  # noqa: E402
from bdbot import scales as SC  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402
from simbot.cutoff import barker_henderson_diameter  # noqa: E402

CASE = "sediment-pmma-3d"
SIGMA_LJ = 2.0 ** (-1.0 / 6.0)      # so the WCA minimum sits at r = d  ★
R_CUT = 1.0
DISPLACEMENT_GATE = 0.03            # overdamped-stability.md
FIT_LO = 2.0                        # the wall's Gaussian tail is 4.5e-9 kT here
BIN_D = 0.5                         # profile bin width, in diameters
#  ★ the TAIL window, where the fluid is ideal enough that l_g carries an
#    implementation check. Over [55, 100] d the CS-equilibrium profile has
#    Z_CS(phi_eff) <= 1.019 and an apparent decay length of 13.495 d, i.e. +0.90 %
#    against the imposed 13.375 -- so a 5 % tolerance there is a real check.
#    ⚠ Over the DENSE window [2, 6 l_g] the same profile gives 15.44 d (+15.5 %),
#      so a check against 13.375 THERE would call a correct result a bug. That is
#      what the first version of this case did.
TAIL_LO, TAIL_HI = 55.0, 100.0


# ════════════════════════════════════════════════════════════════════════
# 1. the physical system (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path) -> dict:
    raw = yaml.safe_load(path.read_text())
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": load_node(raw["particle"]["diameter"]),
        "rho_p": load_node(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": load_node(raw["medium"]["temperature"]),
        "eta": load_node(raw["medium"]["viscosity"]),
        "rho_f": load_node(raw["medium"]["density"]),
        "l_g": load_node(raw["external"]["gravitational_length"]),
        "g": load_node(raw["external"]["g"]),
        "phi": float(raw["geometry"]["volume_fraction"]["value"]),
        "L_z": load_node(raw["geometry"]["cell_height"]),
        "L_xy": load_node(raw["geometry"]["cross_section"]),
        "eps_w": float(raw["external"]["boundary"]["epsilon_kT"]),
        "sig_w": float(raw["external"]["boundary"]["sigma_d"]),
        "eps_wca": float(raw["interactions"]["epsilon_kT"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# 2. the scale table
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_: dict, lxy_d: float, num: dict) -> SC.ScaleLedger:
    """The ledger. `dt` and `T_obs` are in it because both derive from the
    timescales it orders, and a `dt` outside the ledger cannot be seen to violate
    a separation.

    ⚠ The arrow runs `l_g -> T`, not the other way: the paper states `l_g` and
      not `T`, and `T` was solved from it. So recomputing `l_g` here from `T`
      would be circular and is NOT done. Same for `D_0` and `eta`. The four
      independent measured inputs are `d`, `l_g`, `D_0`, `phi`
      (`intake/sediment-pmma-3d/system.yaml`, the `independent_inputs` note).
    """
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, D0, tau_d = b["kT"], b["gamma"], b["D_t"], b["tau_B"]
    l_g = sys_["l_g"].value.to("m")
    L_z = sys_["L_z"].value.to("m")
    L_xy = (Q(lxy_d, "dimensionless") * d).to("m")

    tau_sed = (l_g ** 2 / D0).to("s")          # the slowest mode AT equilibrium
    tau_fall = (L_z * l_g / D0).to("s")        # ★ the DRIFT time across the cell
    dt = Q(float(num["dt_star"]), "dimensionless") * tau_d
    T_obs = float(num["production_tau_sed"]) * tau_sed

    lg = SC.ScaleLedger()
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("l_g", l_g, "kT/(buoyant weight) -- MEASURED, Table 1", star=True)
    lg.add_length("L_xy", L_xy, "lateral box (periodic)")
    lg.add_length("L_z", L_z, "cell height (walls at both ends)", role="box")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_d", tau_d, "d^2/D_0 diffusion (reference)")
    lg.add_time("tau_sed", tau_sed, "l_g^2/D_0 -- the slowest equilibrium mode", star=True)
    lg.add_time("tau_fall", tau_fall, "L_z l_g/D_0 -- drift across the cell ★")
    lg.add_time("T_obs", T_obs, "observation window", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("eps_wall", (sys_["eps_w"] * kT).to("J"), "wall depth at contact")
    lg.add_energy("eps_wca", (sys_["eps_wca"] * kT).to("J"), "WCA amplitude")
    lg.derived = {"gamma": gamma, "D_0": D0, "m": b["m"], "kT": kT, "d": d,
                  "l_g": l_g, "L_z": L_z, "L_xy": L_xy, "tau_d": tau_d,
                  "tau_sed": tau_sed, "tau_fall": tau_fall, "dt": dt,
                  "T_obs": T_obs}
    lg.ref = SC.thermal_reference(
        d, kT, tau_d, time_symbol="tau_d",
        rationale=SC.THERMAL_RATIONALE
        + " The length unit is the DIAMETER and the clock is "
        "tau_d = d^2/D_0, shared with passive-sphere--equilibrium-structure so the "
        "two cards are comparable. The gravitational length is the OTHER length "
        "and becomes a dimensionless group (l_g/d), not a unit -- it spans orders "
        "of magnitude across the source's three solvent mixtures while d does not.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# 3. dimensionless groups + 4. separation checks
# ════════════════════════════════════════════════════════════════════════
def f_wca(r: float, eps: float) -> float:
    """|F| of the WCA core at separation `r` [d], in kT/d. Used only for the `dt`
    gate, so it is the force the run can actually sample that matters."""
    x = SIGMA_LJ / r
    return 4.0 * eps * (12.0 * x ** 12 / r - 6.0 * x ** 6 / r)


def u_wca(r: float, eps: float) -> float:
    x = SIGMA_LJ / r
    return 4.0 * eps * (x ** 12 - x ** 6) + eps


def closest_sampled_separation(eps: float, u_kT: float = 8.0) -> float:
    """The separation at which `U_WCA = u_kT`. Below that the Boltzmann weight is
    `e^-8 = 3.4e-4`, so it is the closest approach the run realistically samples
    and therefore the force the `dt` gate should be set from. Solved, not guessed."""
    lo, hi = 0.5, R_CUT
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        if u_wca(mid, eps) > u_kT:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def analyze_scales(sys_, lg, lxy_d: float):
    d, kT = lg.derived["d"], lg.derived["kT"]
    l_g, L_z = lg.derived["l_g"], lg.derived["L_z"]
    dt, T_obs = lg.derived["dt"], lg.derived["T_obs"]
    lg_star = float((l_g / d).to(""))
    lz_star = float((L_z / d).to(""))
    phi = sys_["phi"]
    #  the local phi the profile reaches at the wall: mean over [0, L_z] of
    #  phi0 exp(-z/l_g) equals phi, solved for phi0
    phi_wall = phi * lz_star / (lg_star * (1.0 - math.exp(-lz_star / lg_star)))
    d_bh = barker_henderson_diameter(sys_["eps_wca"], sigma_lj=SIGMA_LJ)

    groups = [
        ND.Group("phi", phi, None, None, "N(pi/6)/(L_xy^2 L_z)",
                 "mean volume fraction (the paper's 0.017)"),
        ND.Group("l_g/d", lg_star, ("lengths", "l_g"), ("lengths", "d"),
                 "", "gravitational length in diameters *"),
        ND.Group("L_z/l_g", lz_star / lg_star, ("lengths", "L_z"), ("lengths", "l_g"),
                 "", "how many gravitational lengths the cell is"),
        ND.Group("L_xy/d", lxy_d, ("lengths", "L_xy"), ("lengths", "d"),
                 "", "lateral box in diameters (ours, not the paper's 96)"),
        ND.Group("phi_wall", phi_wall, None, None,
                 "phi(L_z/l_g)/(1-e^-Lz/lg)",
                 "predicted local phi at the wall -- how far up CS this reaches"),
        ND.Group("d_BH/d", d_bh, None, None, "Barker-Hend.",
                 "effective hard-sphere diameter at this eps_WCA ★"),
        ND.Group("eps_wall/kT", sys_["eps_w"], ("energies", "eps_wall"),
                 ("energies", "kT"), "", "wall depth -- penetration ~ e^-eps"),
    ]

    #  the dt gates, all three, so the binding one is visible (the card's rule)
    dt_star = float((dt / lg.derived["tau_d"]).to(""))
    r_close = closest_sampled_separation(sys_["eps_wca"])
    f_wall = sys_["eps_w"] / (sys_["sig_w"] * math.sqrt(math.e))
    gates = {
        "thermal": DT.dt_max_thermal(DISPLACEMENT_GATE, 1.0, 1.0),
        "wall": DT.dt_max_force(DISPLACEMENT_GATE, 1.0, 1.0, f_wall),
        "wca": DT.dt_max_force(DISPLACEMENT_GATE, 1.0, 1.0,
                               f_wca(r_close, sys_["eps_wca"])),
        "gravity": DT.dt_max_force(DISPLACEMENT_GATE, 1.0, 1.0, 1.0 / lg_star),
    }
    binding = min(gates, key=lambda k: gates[k])

    checks = [
        C.Check("geometry", "r_c/(L_xy/2)", R_CUT / (lxy_d / 2.0), 1.0,
                note="minimum image -- the WCA cutoff must fit in half the lateral box"),
        C.Check("geometry", "exp(-L_z/l_g)", math.exp(-lz_star / lg_star), 1e-3,
                note="how empty the top of the cell is. With a wall at each end "
                     "there is no periodic lid, so this bounds how much the TOP "
                     "wall can matter rather than a wrap-around artefact"),
        C.Check("geometry", "l_g/d", 1.0 / lg_star, 1.0,
                note="l_g must exceed d, or the profile is not resolved by the "
                     "particle size"),
        C.Check("integration", f"dt/dt_max[{binding}]", dt_star / gates[binding], 1.0,
                note=f"the BINDING gate is {binding}. thermal="
                     f"{gates['thermal']:.2e} wall={gates['wall']:.2e} "
                     f"wca(U=8kT at r={r_close:.3f}d)={gates['wca']:.2e} "
                     f"gravity={gates['gravity']:.2e}"),
        C.Check("model", "tau_p/dt", float((lg.derived["m"] / lg.derived["gamma"]
                                            / dt).to("")), 1.0,
                note="overdamped: momentum relaxation must be unresolved"),
        C.Check("statistics", "tau_sed/T_obs",
                float((lg.derived["tau_sed"] / T_obs).to("")), 1.0,
                note="the window must cover the slowest equilibrium mode. ⚠ the "
                     "DRIFT time is L_z l_g/D_0 = "
                     f"{float((lg.derived['tau_fall'] / lg.derived['tau_d']).to('')):.0f}"
                     " tau_d, 9x longer -- which is why the run starts FROM the "
                     "analytic profile rather than from uniform",
                hard=False),
    ]
    return groups, checks, phi_wall, d_bh, gates, binding


def report_blocks(sys_, lg, n_prod, lxy_d, gates, binding):
    inp = [R.kv(k, f"{sys_[k].value:~.4gP}", sys_[k].tier, sys_[k].source[:44], val_w=20)
           for k in ("d", "l_g", "T", "eta", "rho_p", "rho_f")]
    inp.append(f"  N      = {sys_['N']}   phi = {sys_['phi']}")
    der = [f"  {k:<8} = {lg.derived[k].to_compact():~.4gP}"
           for k in ("gamma", "D_0", "kT", "tau_d", "tau_sed", "tau_fall")]
    der.append("  ★ only d, l_g, D_0, phi are independent; T and eta are derived "
               "from them, so recomputing either is circular")
    plan = [
        f"  dt       = {lg.derived['dt'].to_compact():~.4gP}"
        f"   (= {float((lg.derived['dt'] / lg.derived['tau_d']).to('')):.2e} tau_d)",
        f"  binding dt gate = {binding}  ({gates[binding]:.2e} tau_d)",
        f"  T_obs    = {lg.derived['T_obs'].to_compact():~.4gP}"
        f"   (= {float((lg.derived['T_obs'] / lg.derived['tau_sed']).to('')):.1f} tau_sed)",
        f"  box      = {lxy_d:g} x {lxy_d:g} x "
        f"{float((lg.derived['L_z'] / lg.derived['d']).to('')):g} d",
        f"  steps    = {n_prod:,}   x N={sys_['N']}",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# 5. the analytic solution -- the ground truth this case is built on
# ════════════════════════════════════════════════════════════════════════
def cs_hydrostatic_profile(l_g: float, lz: float, phi_mean: float, n: int = 4001):
    """The equilibrium profile a Carnahan-Starling fluid would have in this cell.

    Hydrostatic balance `dP/dz = -rho kT/l_g` with `P = rho kT Z(phi)` gives

        (Z + phi dZ/dphi) dphi/dz = -phi/l_g

    integrated down from `phi(0)`, with `phi(0)` solved so the mean over the cell
    is `phi_mean`. Returns `(z, phi)`.

    ★ Why this exists, and why using it as an INITIAL CONDITION is not circular.
    The ideal-gas profile is a long way from this one -- `phi(0) = 0.1589` against
    `0.1031` -- and the relaxation between them is the slowest thing in the
    problem. Measured: from an ideal-gas start, 8 `tau_sed` moved `phi(0)` only
    77 % of the way and the drift check still read 6.3 sigma. Extrapolating, ~25
    `tau_sed` per seed would be needed, which is 4.7 h for the primary arm alone.

    So the primary arm starts HERE and the measurement is whether the profile
    **stays**. That is falsifiable: if the WCA fluid's equilibrium is not the CS
    one, the profile moves and the L4 drift check reports it, and `Z` is read from
    the profile the run actually produces rather than from this one. A separate
    arm still starts from the ideal gas, so the convergence evidence is not
    assumed either.
    """
    from scipy.integrate import solve_ivp
    from scipy.optimize import brentq

    def dz_dphi(p):
        p = max(p, 1e-14)
        num = 1.0 + p + p * p - p ** 3
        den = (1.0 - p) ** 3
        dzdp = ((1.0 + 2.0 * p - 3.0 * p * p) * den + num * 3.0 * (1.0 - p) ** 2) / den ** 2
        return -(p / l_g) / (z_carnahan_starling(np.array([p]))[0] + p * dzdp)

    zz = np.linspace(0.0, lz, n)

    def solve(phi0):
        sol = solve_ivp(lambda z, y: [dz_dphi(y[0])], [0.0, lz], [phi0],
                        dense_output=True, rtol=1e-10, atol=1e-16)
        return np.clip(sol.sol(zz)[0], 0.0, None)

    phi0 = brentq(lambda p: float(np.trapezoid(solve(p), zz)) / lz - phi_mean,
                  1e-3, 0.40, xtol=1e-12)
    return zz, solve(phi0)


def cs_height_sampler(l_g: float, lz: float, phi_mean: float, *, margin: float):
    """A `placement.rsa_slab` height sampler drawing from `cs_hydrostatic_profile`
    by inverse CDF, truncated to the margins."""
    zz, phi = cs_hydrostatic_profile(l_g, lz, phi_mean)
    keep = (zz >= margin) & (zz <= lz - margin)
    z, w = zz[keep], np.maximum(phi[keep], 0.0)
    #  ⚠ TRAPEZOID, not cumsum. A `cumsum` CDF left the sampled mean 0.44 % below
    #    the target and KS sat on the 95 % critical value (D = 0.00291 against
    #    0.00304 at N = 200,000) -- marginal rather than correct. The trapezoid
    #    rule removes it: the bias was the half-interval a left-Riemann sum drops
    #    at each step, which is the same error class as the half-bin in
    #    `eos_from_profile`.
    cdf = np.concatenate([[0.0], np.cumsum(0.5 * (w[1:] + w[:-1]) * np.diff(z))])
    cdf = cdf / cdf[-1]
    #  strictly increasing, so np.interp inverts it
    uniq = np.concatenate([[True], np.diff(cdf) > 0])

    def sample(rng, k):
        return np.interp(rng.random(k), cdf[uniq], z[uniq])

    return sample


def z_carnahan_starling(phi):
    """`Z = P/(rho kT)` for hard spheres. Diverges at phi = 1; the caller must
    keep phi below freezing (0.494) for it to mean anything."""
    p = np.asarray(phi, dtype=float)
    return (1.0 + p + p * p - p ** 3) / (1.0 - p) ** 3


def eos_from_profile(phi_z, bin_d: float, l_g: float):
    """`Z = P/(rho kT)` read from the profile by the hydrostatic relation.

    In equilibrium `dP/dz = -rho(z) kT/l_g`, so integrating downward from the top
    of the cell gives `P(z) = (kT/l_g) * integral_z^Lz rho dz'` and therefore

        Z(z) = (1/(l_g rho(z))) * integral_z^Lz rho dz'

    ★ The `6/pi` converting `phi` to a number density cancels between numerator
    and denominator, so `phi` goes in directly.

    ★★ **This never touches HOOMD's pressure**, which is the point:
    `knowledge/wiki/findings/wca-reproduces-carnahan-starling.md` trap 2 records
    that HOOMD's pressure carries a kinetic term built from velocities, and
    overdamped BD velocities are not physical (measured `kT_kin` 1.0009-1.015),
    so `P/(rho kT)` taken from the engine is wrong by construction here.

    For an ideal gas the profile is exactly `exp(-z/l_g)`, the integral above `z`
    is exactly `l_g rho(z)`, and `Z == 1` identically -- which is what
    `tests/test_sediment_eos.py` checks before any run is believed.

    ★★ **The half-bin matters, and by exactly the size of the threshold.**
    A plain `cumsum` of the bins at or above `z` integrates from the bin's LOWER
    EDGE, not from its midpoint, so it overestimates the column by `bin/2` and
    `Z` comes out high by `bin/(2 l_g)`. Measured on an exact exponential before
    this correction existed:

        l_g = 13.375, bin = 0.25 -> +0.71 %   (bin/2l_g = 0.93 %)
        l_g = 13.375, bin = 0.50 -> +1.65 %   (bin/2l_g = 1.87 %)   <- production
        l_g =  4.000, bin = 1.00 -> +12.79 %  (bin/2l_g = 12.5 %)

    **+1.65 % at the production settings, against a pre-registered 2 % band** --
    it would have consumed the whole threshold and been reported as a deviation
    from Carnahan-Starling. The fix is to take only HALF of the bin containing
    `z` plus the whole bins above it, which is the midpoint rule on
    `[z_mid, z_hi]` and leaves an O((bin/l_g)^2) residual: 0.14 % here.
    """
    phi_z = np.asarray(phi_z, dtype=float)
    above = np.concatenate([np.cumsum(phi_z[::-1])[::-1][1:], [0.0]])   # bins ABOVE
    tail = (0.5 * phi_z + above) * bin_d
    with np.errstate(divide="ignore", invalid="ignore"):
        return np.where(phi_z > 0, tail / (l_g * phi_z), np.nan)


def fit_decay_length(h_mid, counts, lo, hi):
    """Weighted least squares on `ln n(h)`. The SLOPE is fitted, so neither the
    origin nor the wall's own structure layer can shift the answer.

    ⚠ This estimator is biased at small counts -- measured +1.93 % with 11 usable
      bins and +9.17 % with 5, on data drawn from an exact exponential
      (`tests/test_placement.py`). The window below is chosen to keep the bin
      count high, and the fit is reported with its standard error so the
      comparison is not made on the point value alone.
    """
    m = (h_mid >= lo) & (h_mid <= hi) & (counts >= 20)
    if m.sum() < 5:
        return None
    x, y, w = h_mid[m], np.log(counts[m]), counts[m].astype(float)
    W = w.sum()
    xm, ym = (w * x).sum() / W, (w * y).sum() / W
    sxx = (w * (x - xm) ** 2).sum()
    slope = (w * (x - xm) * (y - ym)).sum() / sxx
    resid = y - (ym + slope * (x - xm))
    se = math.sqrt((w * resid ** 2).sum() / max(1, int(m.sum()) - 2) / sxx)
    return {"l_g": -1.0 / slope, "se": se / slope ** 2, "nbins": int(m.sum())}


# ════════════════════════════════════════════════════════════════════════
# 6. the build -- L5
# ════════════════════════════════════════════════════════════════════════
@RUN.builder(CASE)
def build(spec, outdir=None) -> RUN.Build:
    import gsd.hoomd
    import hoomd
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    n = int(P["N"])
    lxy, lz, l_g = float(P["L_xy_star"]), float(P["L_z_star"]), float(P["l_g_star"])
    eps_w, sig_w = float(P["eps_wall_kT"]), float(P["sig_wall_d"])
    rcut_w, eps_wca = float(P["rcut_wall_d"]), float(P["eps_wca_kT"])
    margin, min_sep = float(P["wall_margin_d"]), float(P["min_sep_d"])
    dt = float(Nm["dt_star"])
    seed = int(Nm["seed"])

    init = str(P.get("init", "ideal"))
    if init == "ideal":
        sampler = PL.exponential_height(l_g, lz, margin=margin)
    elif init == "cs":
        sampler = cs_height_sampler(l_g, lz, float(P["phi"]), margin=margin)
    else:
        raise ValueError(f"init must be 'ideal' or 'cs', got {init!r}")
    rng = np.random.default_rng(seed)
    pos = PL.rsa_slab(n, lxy=lxy, lz=lz, min_sep=min_sep, rng=rng, margin=margin,
                      height_sampler=sampler)
    placed_min_sep = PL.min_separation_slab(pos, lxy=lxy)

    dev = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=seed)
    frame = gsd.hoomd.Frame()
    frame.particles.N = n
    frame.particles.position = pos
    frame.particles.typeid = np.zeros(n, dtype=int)
    frame.particles.types = ["A"]
    frame.particles.diameter = np.ones(n)
    frame.configuration.box = [lxy, lxy, lz, 0, 0, 0]
    sim.create_state_from_snapshot(frame)

    #  WCA: LJ cut at its minimum and shifted up. sigma_LJ = 2^(-1/6) d puts that
    #  minimum at r = d, so the nominal diameter IS d.  ★ see the module docstring
    lj = md.pair.LJ(nlist=md.nlist.Cell(buffer=0.3), default_r_cut=R_CUT,
                    mode="shift")
    lj.params[("A", "A")] = dict(epsilon=eps_wca, sigma=SIGMA_LJ)
    #  a BOUNDED wall. An LJ wall reaches 1517 kT/d at h = 0.4 d and a particle
    #  passes through with no NaN (bd-hoomd trap 22).
    wall = md.external.wall.Gaussian(walls=[
        hoomd.wall.Plane(origin=(0, 0, -lz / 2), normal=(0, 0, 1)),
        hoomd.wall.Plane(origin=(0, 0, lz / 2), normal=(0, 0, -1))])
    wall.params["A"] = {"epsilon": eps_w, "sigma": sig_w, "r_cut": rcut_w}
    grav = md.force.Constant(filter=hoomd.filter.All())
    grav.constant_force["A"] = (0.0, 0.0, -1.0 / l_g)
    grav.constant_torque["A"] = (0.0, 0.0, 0.0)

    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0)
    bd.gamma["A"] = 1.0
    forces = [lj, wall, grav]
    sim.operations.integrator = md.Integrator(dt=dt, methods=[bd], forces=forces)

    #  ── accumulators. `sample()` must NOT return the N-particle array every
    #     frame: observables.npz went 448 KB -> 148 MB that way
    #     (docs/05-pitfalls.md).
    nbins = int(round(lz / BIN_D))
    edges = np.linspace(0.0, lz, nbins + 1)
    acc = {"hist": np.zeros(nbins), "frames": 0, "ref_xy": None, "ref_t": None}

    def _snap():
        s = sim.state.get_snapshot()
        p = np.array(s.particles.position, copy=True)
        img = np.array(s.particles.image, copy=True)
        return p, img

    def sample(timestep, phase):
        p, img = _snap()
        h = p[:, 2] + lz / 2
        t = timestep * dt
        if phase == "production":
            acc["hist"] += np.histogram(h, bins=edges)[0]
            acc["frames"] += 1
            unwrapped = p[:, :2] + img[:, :2] * lxy
            if acc["ref_xy"] is None:
                acc["ref_xy"], acc["ref_t"] = unwrapped, t
                msd_xy = 0.0
            else:
                dxy = unwrapped - acc["ref_xy"]
                msd_xy = float((dxy ** 2).sum(axis=1).mean())
        else:
            msd_xy = float("nan")
        dr = p[:, None, :] - p[None, :, :]
        dr[:, :, :2] -= lxy * np.round(dr[:, :, :2] / lxy)
        d2 = (dr ** 2).sum(-1)
        np.fill_diagonal(d2, np.inf)
        return dict(
            t=t, h_mean=float(h.mean()), h_min=float(h.min()), h_max=float(h.max()),
            n_below_wall=float((h < 0).sum()), n_above_cell=float((h > lz).sum()),
            min_sep=float(math.sqrt(d2.min())), msd_xy=msd_xy,
            u_pair=float(np.asarray(lj.energies).sum()) / n,
            u_wall=float(np.asarray(wall.energies).sum()) / n)

    def pe_per_particle() -> float:
        """Total potential energy per particle in kT, gravity included. The
        gravitational term is what actually equilibrates here, so leaving it out
        would make the drift check blind to the one mode that is slow."""
        p, _ = _snap()
        h = p[:, 2] + lz / 2
        u = (float(np.asarray(lj.energies).sum())
             + float(np.asarray(wall.energies).sum())
             + float(h.sum()) / l_g)
        return u / n

    def finalize(cols) -> dict:
        frames = max(acc["frames"], 1)
        mid = 0.5 * (edges[1:] + edges[:-1])
        counts = acc["hist"] / frames                      # per frame
        slab_vol = lxy * lxy * BIN_D
        phi_z = counts * (math.pi / 6.0) / slab_vol
        #  Poisson error on the accumulated counts, propagated to phi
        with np.errstate(divide="ignore", invalid="ignore"):
            phi_se = np.where(acc["hist"] > 0,
                              phi_z / np.sqrt(np.maximum(acc["hist"], 1.0)), 0.0)

        hi = min(6.0 * l_g, lz - 2.0)
        fit = fit_decay_length(mid, acc["hist"], FIT_LO, hi)
        fit_tail = fit_decay_length(mid, acc["hist"], TAIL_LO, min(TAIL_HI, lz - 2.0))

        #  ── the EOS, from the hydrostatic relation. Z = 1 in the dilute tail is
        #     a built-in check that this read is right, and it needs no input.
        #     P*(z) = (1/l_g) integral_z^Lz rho* dz'  ->  Z = P*/rho*
        #     the (6/pi) in rho* cancels, so phi can be used directly.
        z_eos = eos_from_profile(phi_z, BIN_D, l_g)

        d_bh = barker_henderson_diameter(eps_wca, sigma_lj=SIGMA_LJ)
        phi_eff = phi_z * d_bh ** 3
        cs_eff = z_carnahan_starling(np.clip(phi_eff, 0, 0.45))
        cs_nom = z_carnahan_starling(np.clip(phi_z, 0, 0.45))

        #  the comparison window: above the wall layer, and only where a slab has
        #  enough counts for Z to mean anything
        ok = (mid >= FIT_LO) & (acc["hist"] >= 200) & np.isfinite(z_eos)
        dev_eff = (100.0 * (z_eos[ok] - cs_eff[ok]) / cs_eff[ok]
                   if ok.any() else np.array([np.nan]))
        dev_nom = (100.0 * (z_eos[ok] - cs_nom[ok]) / cs_nom[ok]
                   if ok.any() else np.array([np.nan]))

        #  ★★ The PRIMARY statistic is a Poisson-WEIGHTED MEAN, not the maximum.
        #
        #  ⚠ The max over slabs is unusable at this sampling, and that was
        #    measured rather than suspected. A 0.5 d slab at phi = 0.05 collects
        #    6.9 particles per frame, so 200 frames give ~1370 counts = 2.7 %
        #    Poisson on that slab's Z. The maximum over ~40 such slabs is then
        #    about 2.5 sigma * 2.7 % = 7 % from NOISE ALONE -- it would break a
        #    2 % band whatever the physics did. The weighted mean over the same
        #    window averages that down to a few tenths of a percent.
        #
        #    `goal.yaml`'s answering quantity was changed from the max to this
        #    BEFORE sealing, which is what the pre-registration window is for.
        #    The max is still reported, as a `measurement`, so the noise level
        #    stays visible next to the statistic that carries the verdict.
        w = acc["hist"][ok].astype(float) if ok.any() else np.array([1.0])
        def _wmean(dev):
            m = np.isfinite(dev)
            return float(np.average(dev[m], weights=w[m])) if m.any() else float("nan")
        wdev_eff, wdev_nom = _wmean(dev_eff), _wmean(dev_nom)

        #  and the single most discriminating point: the densest usable slab
        if ok.any():
            i_dense = int(np.argmax(np.where(ok, phi_z, -1.0)))
            dense = {"phi": float(phi_z[i_dense]), "phi_eff": float(phi_eff[i_dense]),
                     "Z": float(z_eos[i_dense]), "cs_eff": float(cs_eff[i_dense]),
                     "cs_nom": float(cs_nom[i_dense]),
                     "counts": float(acc["hist"][i_dense])}
        else:
            dense = {k: float("nan") for k in
                     ("phi", "phi_eff", "Z", "cs_eff", "cs_nom", "counts")}

        #  ── the dilute tail. ⚠ A separate, LOWER count threshold on purpose:
        #     the tail is exponentially sparse by construction, so requiring the
        #     EOS window's 200 counts there is self-defeating -- it produced
        #     `nan` on the first smoke run. At phi = 0.005 a slab collects
        #     phi L_xy^2 dz/(pi/6) = 0.69 particles per frame, so 200 frames give
        #     137 and the two conditions barely intersect. 50 counts is 14 %
        #     Poisson on one slab, and Z is a ratio whose numerator is a
        #     cumulative sum over everything above, so the denominator dominates;
        #     averaging several slabs brings it down again.
        dilute = (mid >= FIT_LO) & (acc["hist"] >= 50) & np.isfinite(z_eos) \
            & (phi_eff < 0.01)
        z_tail = float(np.nanmean(z_eos[dilute])) if dilute.any() else float("nan")
        n_tail = int(dilute.sum())

        #  lateral diffusion, from the linear part of the MSD
        t = np.asarray(cols.get("t", []), dtype=float)
        msd = np.asarray(cols.get("msd_xy", []), dtype=float)
        good = np.isfinite(msd) & (msd > 0)
        d_xy = float("nan")
        if good.sum() >= 10:
            tt, mm = t[good], msd[good]
            tt = tt - tt[0]
            half = tt > 0.5 * tt.max()          # the late, linear part
            if half.sum() >= 5:
                d_xy = float(np.polyfit(tt[half], mm[half], 1)[0] / 4.0)  # 2D: 4Dt

        obs = [
            MET.observable(
                "l_g_fitted_tail", fit_tail["l_g"] if fit_tail else float("nan"),
                predicted=l_g, unit="d", role="implementation_check", tol_pct=5.0,
                sigma=fit_tail["se"] if fit_tail else None,
                source="paper Table 1, l_g = 10.7 um / d = 0.8 um",
                note=f"the decay length over the DILUTE TAIL, h in "
                     f"[{TAIL_LO:g}, {TAIL_HI:g}] d. ★ This is the "
                     f"implementation check, not the dense-window fit",
                derivation="Where phi -> 0 the fluid is ideal and the profile is "
                           "exactly n0 exp(-z/l_g), so the apparent decay length "
                           "IS l_g. Quantified on the CS-equilibrium profile: "
                           "over [55, 100] d it has Z_CS(phi_eff) <= 1.019 and an "
                           "apparent l_g of 13.495 d, +0.90 % against the imposed "
                           "13.375 -- inside the 5 % tolerance, so excluded volume "
                           "cannot fake a failure here. ⚠ The DENSE window gives "
                           "15.44 d (+15.5 %) on the same profile, so a check "
                           "against 13.375 there would call a correct result a "
                           "bug; that fit is reported separately as a hypothesis."),
            MET.observable(
                "l_g_fitted_dense", fit["l_g"] if fit else float("nan"),
                predicted=15.442, unit="d", role="hypothesis", tol_pct=8.0,
                sigma=fit["se"] if fit else None,
                source="cs_hydrostatic_profile, integrated for this cell",
                note=f"the apparent decay length over the dense window h in "
                     f"[{FIT_LO:g}, {hi:.1f}] d, where excluded volume flattens "
                     f"the profile",
                derivation="HYPOTHESIS: the prediction 15.442 d comes from "
                           "integrating the Carnahan-Starling hydrostatic ODE for "
                           "this cell and fitting the same window, so a mismatch "
                           "says the WCA fluid's profile differs from the CS one "
                           "-- which is the same question P3 asks, read off the "
                           "profile instead of the EOS. The ideal-gas value 13.375 "
                           "is 13.4 % away, so the two are separable."),
            MET.observable(
                "Z_dilute_tail", z_tail, predicted=1.0, role="implementation_check",
                tol_pct=5.0, source="ideal gas",
                note="Z -> 1 as phi -> 0. Costs nothing and catches a wrong "
                     "hydrostatic read",
                derivation="With n = n0 exp(-z/l_g), the integral of rho above z "
                           "is exactly l_g rho(z), so Z = 1 identically. Any "
                           "error in the cumulative sum, the bin width or the "
                           "l_g used shows up here before it reaches the EOS."),
            MET.observable(
                "Z_dev_wmean_vs_CS_eff", wdev_eff, predicted=0.0, unit="percent",
                role="hypothesis", sigma=2.0, tol_sigma=1.0,
                source="knowledge/wiki/findings/wca-reproduces-carnahan-starling.md",
                note="Poisson-weighted mean of (Z - CS(phi_eff))/CS over the "
                     "window, in percent. ★ THE answering quantity of "
                     "intake/sediment-pmma-3d/goal.yaml",
                derivation="HYPOTHESIS, not a check: the simulation imposes a WCA "
                           "core and imposes no equation of state at all, so a "
                           "mismatch is a result about how well a soft core "
                           "stands in for a hard sphere. The 2 % band is "
                           "INHERITED from this repository's own prior "
                           "measurement of the same substitution (-0.32 % at "
                           "phi ~ 0.10, -0.51 % at phi ~ 0.21 at eps = 1 kT), not "
                           f"invented. phi_eff = phi (d_BH/d)^3 with d_BH/d = "
                           f"{d_bh:.4f}."),
            MET.observable(
                "Z_dev_wmean_vs_CS_nominal", wdev_nom, unit="percent",
                role="measurement",
                note="the same weighted mean WITHOUT the Barker-Henderson "
                     "mapping. The GAP between this and the line above is what "
                     "says whether the mapping is required",
                derivation="CS(phi_eff) and CS(phi_nominal) differ by -17.05 % at "
                           "the wall of this cell, so the two are separable by "
                           "this measurement. No prediction is attached: this "
                           "number IS the discriminator."),
            MET.observable(
                "Z_dev_max_vs_CS_eff", float(np.nanmax(np.abs(dev_eff))),
                unit="percent", role="measurement",
                note="the LARGEST single-slab deviation. Reported without a "
                     "verdict on purpose: at 200 frames a slab carries ~2.7 % "
                     "Poisson error, so the max over ~40 slabs is ~7 % from "
                     "noise alone. It is here to show the noise floor beside the "
                     "weighted mean, not to be compared with 2 %"),
            MET.observable(
                "Z_at_densest_slab", dense["Z"], predicted=dense["cs_eff"],
                role="hypothesis", tol_pct=8.0,
                note=f"the single most discriminating point: phi = "
                     f"{dense['phi']:.4f}, phi_eff = {dense['phi_eff']:.4f}, "
                     f"{dense['counts']:.0f} counts. CS(phi_nominal) there is "
                     f"{dense['cs_nom']:.4f}, which is the alternative this "
                     f"point separates from",
                derivation="One slab, so it carries the full Poisson error "
                           "(~1/sqrt(counts)); the 8 % tolerance is that error "
                           "plus room, not a physics claim. Its value is that "
                           "CS(phi_eff) and CS(phi_nominal) are furthest apart "
                           "here, so it is where the mapping question is "
                           "sharpest."),
            MET.observable(
                "D_xy", d_xy, predicted=1.0, unit="d^2/tau_d",
                role="implementation_check", tol_pct=10.0,
                source="D = kT/gamma, the input",
                note="lateral self-diffusion in reduced units -- the clock, "
                     "checked against itself",
                derivation="x and y are unconfined and periodic, so the lateral "
                           "MSD is 4 D t with D = D_0 = 1 in these units. ⚠ the "
                           "VERTICAL MSD saturates at ~l_g^2 and must not be fit "
                           "(the card's gate table). Excluded volume at "
                           "phi <~ 0.16 lowers the long-time D by a few percent, "
                           "hence 10 % rather than 5 %."),
            MET.observable("phi_wall_measured", float(phi_z[mid < 1.5].max())
                           if (mid < 1.5).any() else float("nan"),
                           role="measurement",
                           note="the highest local volume fraction the profile "
                                "reached -- how far up the CS curve this run got"),
            MET.observable("min_sep_placed", placed_min_sep, unit="d",
                           role="measurement",
                           note="closest pair at t=0. Below ~0.8 d the WCA core "
                                "starts to set dt"),
        ]
        return {
            "observables": obs,
            "extra": {
                "l_g_star": l_g, "d_BH_over_d": d_bh, "frames": acc["frames"],
                "n_bins_fitted": fit["nbins"] if fit else 0,
                "n_bins_fitted_tail": fit_tail["nbins"] if fit_tail else 0,
                "fit_window_d": [FIT_LO, hi],
                "tail_window_d": [TAIL_LO, min(TAIL_HI, lz - 2.0)],
                "init": str(P.get("init", "ideal")),
                "eos_window_bins": int(ok.sum()),
                "dilute_tail_bins": n_tail,
                "densest_slab": dense,
                "phi_mean_measured": float(phi_z.mean()),
            },
            "arrays": {
                "profile_h": mid, "profile_phi": phi_z, "profile_phi_se": phi_se,
                "profile_counts": acc["hist"], "eos_Z": z_eos,
                "eos_phi_eff": phi_eff, "eos_CS_eff": cs_eff,
                "eos_CS_nominal": cs_nom, "eos_window": ok.astype(float),
            },
        }

    n_prod = int(Nm["n_prod"])
    return RUN.Build(
        sim=sim, forces=forces, n_particles=n, sample=sample,
        pe_per_particle=pe_per_particle,
        n_eq=int(Nm["n_eq"]), n_prod=n_prod,
        sample_every=int(Nm["sample_every"]),
        gsd_path=(Path(outdir) / "traj_A.gsd") if outdir else None,
        tags=["3d", "sedimentation", "gravity", "wall", "wca", "hard-sphere-eos"],
        physical={"l_g_star": l_g, "L_z_star": lz, "L_xy_star": lxy,
                  "phi": float(P["phi"]), "eps_wca_kT": eps_wca},
        finalize=finalize)


# ════════════════════════════════════════════════════════════════════════
# 7. main
# ════════════════════════════════════════════════════════════════════════
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--lxy", type=float, default=None,
                    help="lateral box in diameters (default: the intake's 12)")
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--tau-sed", type=float, default=None,
                    help="observation window in multiples of tau_sed")
    ap.add_argument("--init", choices=("ideal", "cs"), default="ideal",
                    help="initial profile: the ideal-gas exponential, or the "
                         "Carnahan-Starling hydrostatic equilibrium for this cell")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--spec", action="store_true",
                    help="write specs/<run_id>.json and exit (does not run)")
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--require-seal", action="store_true")
    ap.add_argument("--require-approval", action="store_true")
    args = ap.parse_args()

    sys_ = load_system(ROOT / "intake" / CASE / "system.yaml")
    num = dict(sys_["numerics"])
    d_um = sys_["d"].value.to("um").magnitude
    lxy_d = args.lxy if args.lxy is not None else \
        sys_["L_xy"].value.to("um").magnitude / d_um
    lz_d = sys_["L_z"].value.to("um").magnitude / d_um
    if args.tau_sed is not None:
        num["production_tau_sed"] = args.tau_sed
    if args.smoke:
        num["production_tau_sed"] = 0.3

    #  N is NOT free: phi and the box fix it (structure.size.fixed_by = cost)
    n = int(round(sys_["phi"] * lxy_d * lxy_d * lz_d / (math.pi / 6.0)))
    sys_["N"] = n

    lg = build_ledger(sys_, lxy_d, num)
    dt, T_obs = lg.derived["dt"], lg.derived["T_obs"]
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, n_prod // int(round(
        float(num["samples_per_tau_sed"]) * float(num["production_tau_sed"]))))
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, phi_wall, d_bh, gates, binding = analyze_scales(sys_, lg, lxy_d)

    tag = f"lxy{lxy_d:g}-{args.init}"
    if args.smoke:
        tag += "-smoke"
    if args.tau_sed is not None:
        tag += f"-t{args.tau_sed:g}"
    seed = args.seed if args.seed is not None else int(num["seed_base"])
    tag += f"-s{seed}"

    spec = ND.NondimSpec(
        case=CASE, system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={
            "N": n, "phi": sys_["phi"],
            "L_xy_star": lxy_d, "L_z_star": lz_d,
            "l_g_star": float((lg.derived["l_g"] / lg.derived["d"]).to("")),
            "eps_wca_kT": sys_["eps_wca"],
            "eps_wall_kT": sys_["eps_w"], "sig_wall_d": sys_["sig_w"],
            "rcut_wall_d": float(sys_["_raw"]["external"]["boundary"].get(
                "rcut_d", 2.0)),
            "wall_margin_d": float(num["wall_margin_d"]),
            "min_sep_d": float(num["min_sep_d"]),
            #  ★ `init` is in `params`, so it IS hashed into the run_id: two runs
            #    that differ only in where they started are two different runs.
            #    (`soft_r3_2d` records what happens when an initial condition is
            #    left out of the hash -- two physically different runs collide.)
            "init": args.init,
        },
        numerics={"dt_star": float((dt / lg.derived["tau_d"]).to("")),
                  "n_eq": 0, "n_prod": n_prod, "sample_every": sample_every,
                  "seed": seed},
        tag=tag, nhex=10)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, n_prod, lxy_d, gates, binding)
    plan.append(f"  phi(wall) predicted = {phi_wall:.5f}   "
                f"CS(phi_eff) there = {float(z_carnahan_starling(phi_wall * d_bh**3)):.4f}"
                f"   CS(phi_nom) = {float(z_carnahan_starling(phi_wall)):.4f}")
    report, verdict = R.render(
        title=f"DimensionlessReport — {CASE}   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)
    if spec.errors:
        print(f"\nx {len(spec.errors)} L3 integrity error(s).")
        return 1
    if verdict == "FAIL":
        print("\nx a hard separation check failed -- not running.")
        return 1

    #  ★ ONE guard for both read-only flags. Splitting it let `--report` fall
    #    through and RUN the simulation (soft_r3_2d records the incident).
    if args.spec or args.report:
        if args.spec:
            p = spec.write(ROOT / "specs" / f"{run_id}.json")
            print(f"\nL3 spec: {p.relative_to(ROOT)}")
        return 0
    p = spec.write(ROOT / "specs" / f"{run_id}.json")

    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True,
                    require_seal=args.require_seal,
                    require_approval=args.require_approval)
    print(RUN.render_verdict(v))
    if v["status"] == "skipped":
        return 0
    if v["status"] != RUN.OK:
        return 1

    metrics = json.loads((outdir / "metrics.json").read_text())
    res = metrics.get("result", {})
    lines = ["", "=" * 78, f"RESULT -- {CASE}  {tag}", "=" * 78]
    for o in metrics.get("observables", []):
        pred = "" if o.get("predicted") is None else f"  vs {o['predicted']:.5g}"
        lines.append(f"  {o['name']:<26} {o['measured']:>12.5g}{pred:<16} "
                     f"[{o.get('role','')}] {o.get('verdict','')}")
    (outdir / "result.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    try:
        make_plots(outdir, res)
    except Exception as exc:                      # a plot must not fail a run
        print(f"  (plots skipped: {exc})")
    return 0


def make_plots(outdir: Path, res: dict) -> None:
    """The profile and the EOS. Labels in English — `DejaVu Sans` has no Hangul
    and the fonts that do are missing the minus sign (CLAUDE.md)."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    z = np.load(outdir / "observables.npz")
    h, phi = z["profile_h"], z["profile_phi"]
    se, ok = z["profile_phi_se"], z["eos_window"].astype(bool)
    l_g = float(res.get("l_g_star", np.nan))

    fig, ax = plt.subplots(1, 3, figsize=(16, 4.6))

    m = phi > 0
    ax[0].errorbar(h[m], phi[m], yerr=se[m], fmt="o", ms=3, lw=0.8,
                   label="measured")
    if np.isfinite(l_g) and m.any():
        h0 = h[m][0]
        ax[0].plot(h[m], phi[m][0] * np.exp(-(h[m] - h0) / l_g), "--",
                   label=f"exp(-z/l_g), l_g={l_g:.3f} d (paper)")
    ax[0].set_yscale("log")
    ax[0].set_xlabel("height above the wall  h / d")
    ax[0].set_ylabel(r"local volume fraction  $\phi(h)$")
    ax[0].set_title("1. density profile")
    ax[0].legend(fontsize=8)

    Z, pe = z["eos_Z"], z["eos_phi_eff"]
    cse, csn = z["eos_CS_eff"], z["eos_CS_nominal"]
    ax[1].plot(pe[ok], Z[ok], "o", ms=4, label="measured  $Z$")
    o = np.argsort(pe[ok])
    ax[1].plot(pe[ok][o], cse[ok][o], "-", label=r"Carnahan-Starling$(\phi_{eff})$")
    ax[1].plot(phi[ok][np.argsort(phi[ok])], csn[ok][np.argsort(phi[ok])], ":",
               label=r"Carnahan-Starling$(\phi_{nominal})$")
    ax[1].axhline(1.0, color="k", lw=0.5)
    ax[1].set_xlabel(r"$\phi_{eff} = \phi\,(d_{BH}/d)^3$")
    ax[1].set_ylabel(r"$Z = P/(\rho k_BT)$")
    ax[1].set_title("2. equation of state, both mappings")
    ax[1].legend(fontsize=8)

    with np.errstate(divide="ignore", invalid="ignore"):
        de = 100 * (Z - cse) / cse
        dn = 100 * (Z - csn) / csn
    ax[2].plot(pe[ok], de[ok], "o", ms=4, label=r"vs CS$(\phi_{eff})$")
    ax[2].plot(pe[ok], dn[ok], "s", ms=4, mfc="none",
               label=r"vs CS$(\phi_{nominal})$")
    for y in (-2, 2):
        ax[2].axhline(y, color="r", lw=0.6, ls="--")
    ax[2].axhline(0, color="k", lw=0.5)
    ax[2].set_xlabel(r"$\phi_{eff}$")
    ax[2].set_ylabel("deviation from CS  [%]")
    ax[2].set_title("3. the discriminator (red: the pre-registered 2 % band)")
    ax[2].legend(fontsize=8)

    fig.tight_layout()
    fig.savefig(outdir / "observables.png", dpi=130)
    plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
