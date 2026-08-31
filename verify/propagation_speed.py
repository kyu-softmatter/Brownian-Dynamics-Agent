"""Stress propagation speed -- from the spatial gradient of the per-bead phase lag.

In steady-state oscillation, if the bead s bonds from the driven one lags by phi(s),
the wavenumber is k = -dphi/ds [rad/bond] and the phase velocity is

    v_phase = omega / k = omega*l / |dphi/ds|   [m/s]   (l = bond length)

★ What is actually being measured -- this system is **overdamped** (BD, no inertia),
  so this is not an elastic wave like sound. It is the speed at which bending
  deformation **diffuses** along the chain against viscous resistance.
  The dispersion relation for overdamped bending is omega ~ (kappa_theta*l/gamma)k^4,
  so v_phase depends on omega (it is NOT dispersionless).
  -> this is **a value at one omega**, not a material constant. Always quote omega
  with it.

★ The amplitude decay length lambda is reported too, but **must not be read as a
  "penetration depth".** In three-point bending, most of the amplitude falling off
  away from the centre is the **static beam shape** (the ends are held by traps, so of
  course it is smaller), not decay. The measurement shows this -- at k_t x100,
  lambda=2.6 << 1/k=42.5, so the amplitude dies quickly while the phase barely turns.
  A genuinely diffusive wave would have lambda ~ 1/k.
  -> **the dynamical information is in the phase gradient alone.** lambda is read only
  as a boundary-condition diagnostic.

★★ Coherence is checked FIRST -- if the phase does not reproduce between seeds, v is
  meaningless. For DLVO, coh ~ 0.46-0.59: the interior beads' phase is thermal noise
  unrelated to the drive. Fitting still returns a number, but with an error larger
  than the value and a negative lambda -- do not quote it.

⚠️ The end beads (0 and n-1) are trapped force sensors dominated by the boundary
condition -- they are excluded from the fit.

    $PY scratch/propagation_speed.py
"""
import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
import sys                                                    # noqa: E402
sys.path.insert(0, str(ROOT))
from bdbot import lockin as LI                                # noqa: E402

d_SI = 1.47e-6
ELL = 1.0076                       # bond length [d]
TAU_B = 6.1510                     # s -- for converting omega* -> omega[rad/s]
P = "runs/chain-bend-2d-dlvo__n9-w3000-a1470"

CONDS = [("trap, default k_t",  "",        "-jkr"),
         ("trap, k_t×100",   "-kt100",     "-jkr-kt100"),
         ("position forcing", "-position",  "-jkr-position")]


def phasors(dd):
    """Per-bead complex phasor y_hat_i (lock-in over all samples)."""
    z = np.load(Path(dd) / "observables.npz", allow_pickle=True)
    if "shape_y" not in set(z.files):
        return None, None, None
    s = json.loads((Path(dd) / "spec.json").read_text())
    if s["numerics"]["n_prod"] < 600_000:              # exclude smoke runs
        return None, None, None
    om = float(s["params"]["omega_star"])
    t = np.asarray(z["t"], float)
    ys = np.asarray(z["shape_y"], float)
    out = []
    for i in range(ys.shape[1]):
        b = LI.lockin_blocks(t, ys[:, i] - ys[:, i].mean(), om,
                             n_blocks=min(10, max(2, len(t) // 20)))
        h, _ = LI.agg(b)
        out.append(h)
    return np.array(out), om, int(s["params"]["n_beads"])


def fit(pat):
    Zs, om, n = [], None, None
    for dd in sorted(glob.glob(str(ROOT / pat))):
        if not (Path(dd) / "metrics.json").exists():
            continue
        z, o, nn = phasors(dd)
        if z is None:
            continue
        Zs.append(z); om, n = o, nn
    if not Zs:
        return None
    Z = np.array(Zs)
    mid = n // 2
    Zrel = Z / Z[:, [mid]]                       # referenced to the driven bead
                                                 # (cancels absolute-phase pollution)

    # Interior beads only (the trapped ends excluded), one arm at a time, folded onto
    # |s| and averaged
    # ★ Coherence check -- does the phase profile reproduce between seeds?
    #   coh = |seed-averaged phasor| / mean|phasor|. 1 is fully coherent, 0 is random.
    #   When the chain cannot transmit force (as with DLVO) the interior beads' phase
    #   is thermal noise, coh -> 0, and fitting v then yields **a meaningless number**
    #   (measured: error > value, lambda < 0).
    coh = np.abs(Zrel.mean(0)) / np.abs(Zrel).mean(0)

    res = []
    for zr in Zrel:
        s = np.arange(n) - mid
        keep = (np.abs(s) >= 1) & (np.abs(s) <= mid - 1)       # s=±1..±3 (n=9)
        ss = np.abs(s[keep]).astype(float)
        ph = np.unwrap(np.angle(zr[keep]))
        amp = np.abs(zr[keep])
        # gradient of phi(s) [rad/bond] and of ln|A| -> the decay length
        kk = -np.polyfit(ss, ph, 1)[0]
        lam = -1.0 / np.polyfit(ss, np.log(amp), 1)[0]
        res.append((kk, lam))
    res = np.array(res)
    k_m, k_e = res[:, 0].mean(), res[:, 0].std(ddof=1) / np.sqrt(len(res))
    l_m, l_e = res[:, 1].mean(), res[:, 1].std(ddof=1) / np.sqrt(len(res))
    om_SI = om / TAU_B                                          # rad/s
    v = om_SI * ELL * d_SI / k_m if k_m else np.inf             # m/s
    v_e = abs(v) * (k_e / abs(k_m)) if k_m else np.nan
    ss_all = np.arange(n) - mid
    keep_all = (np.abs(ss_all) >= 1) & (np.abs(ss_all) <= mid - 1)
    coh_in = float(coh[keep_all].mean())
    return dict(k=k_m, k_e=k_e, lam=l_m, lam_e=l_e, v=v, v_e=v_e,
                om_SI=om_SI, n=len(res), coh=coh_in)


print("=" * 104)
print("stress propagation speed -- from the phase gradient, "
      "v = omega*l/|dphi/ds|   (omega = 3000 rad/s, n=9, a=1d)")
print("=" * 104)
COH_MIN = 0.5     # ★ below this the phase is incoherent and v cannot be defined
print(f"{'condition':<16} {'branch':<6} {'coh':>7} {'dphi/ds [rad/bond]':>19} "
      f"{'v [µm/s]':>20} "
      f"{'lambda [bond]':>15} {'1/k':>7}")
rows = []
for lab, sd, sj in CONDS:
    # ★ DLVO is excluded (user instruction, 2026-08-06). Coherence 0.46-0.59 already
    #   established that v is undefined there -- with no stress to propagate there is
    #   no speed either.
    for br, suf in (("JKR", sj),):
        r = fit(f"{P}{suf}__*")
        if r is None:
            print(f"{lab:<16} {br:<6} {'—':>19}")
            continue
        if r["coh"] < COH_MIN:
            print(f"{lab:<16} {br:<6} {r['coh']:>7.3f}   ★ phase incoherent -- "
                  f"v undefined (interior beads move independently of the drive)")
        else:
            print(f"{lab:<16} {br:<6} {r['coh']:>7.3f} {r['k']:>12.5f}±{r['k_e']:<6.4f} "
                  f"{r['v']*1e6:>13.0f}±{r['v_e']*1e6:<6.0f} "
                  f"{r['lam']:>9.3f}±{r['lam_e']:<5.3f} {1/r['k']:>7.2f}")
        rows.append(dict(lab=lab, br=br, **r))
    print()

import pickle                                                  # noqa: E402
pickle.dump(rows, open("/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/"
                       "7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/prop_speed.pkl", "wb"))
print("=" * 104)
print("How to read this -- (1) below coherence 0.5, v is not defined at all "
      "(DLVO is such a case).")
print("           (2) lambda is **NOT** a penetration depth -- the static beam shape "
      "of three-point bending dominates.")
print("              At k_t x100, lambda=2.6 << 1/k=42.5 is the evidence "
      "(a diffusive wave would give lambda ~ 1/k).")
print("           (3) default k_t and position forcing agree at v ~ 30,000 µm/s -- "
      "because in both the end traps")
print("              have the default stiffness. At k_t x100 the ends are effectively "
      "fixed and it is 6x faster.")
print("★ This is **a value at omega=3000 rad/s**. Overdamped bending goes as "
      "omega ~ k^4, so it is dispersive and")
print("  v_phase depends on omega -- it is not a material constant.")
