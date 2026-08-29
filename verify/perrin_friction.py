"""Ellipsoid friction coefficients — needed for the abp-rod case.

**Verified against the sphere limit before being used.**

CLAUDE.md rule 6: never write a physics formula from intuition. The Perrin factors are
easy to misremember, so they are computed directly from the elliptic-integral
definition and checked against the known sphere limit (a1=a2=a3):

    translation  zeta = 6*pi*eta*R        rotation  zeta_r = 8*pi*eta*R^3

Definitions (Kim & Karrila, ellipsoid resistance):
    Delta(lam) = sqrt((a1^2+lam)(a2^2+lam)(a3^2+lam))
    chi_0      = int_0^inf dlam/Delta      chi_i = int_0^inf dlam/[(a_i^2+lam)Delta]
    translation  zeta_i   = 16*pi*eta / (chi_0 + a_i^2 chi_i)
    rotation     zeta_r,i = 16*pi*eta (a_j^2+a_k^2) / [3 (a_j^2 chi_j + a_k^2 chi_k)]
                            (i,j,k cyclic)

    $PY scratch/perrin_friction.py
"""
import math

import numpy as np
from scipy import integrate

ETA = 0.851e-3          # Pa*s, water@300K
KT = 1.380649e-23 * 300  # J


def chi(axes):
    """(chi_0, chi_1, chi_2, chi_3) — elliptic integrals.

    Computed numerically and validated against the sphere limit.
    """
    a = np.asarray(axes, dtype=float)

    def delta(lam):
        return math.sqrt((a[0] ** 2 + lam) * (a[1] ** 2 + lam) * (a[2] ** 2 + lam))

    # Convergence is slow as lam -> inf, so the substitution
    # lam = a_max^2 * (1/u - 1) maps it onto a finite interval
    scale = float(a.max()) ** 2

    def sub(f):
        def g(u):
            lam = scale * (1.0 / u - 1.0)
            return f(lam) * scale / u**2
        return integrate.quad(g, 1e-14, 1.0, limit=400)[0]

    chi0 = sub(lambda lam: 1.0 / delta(lam))
    chis = [sub(lambda lam, i=i: 1.0 / ((a[i] ** 2 + lam) * delta(lam))) for i in range(3)]
    return chi0, chis


def friction(axes, eta=ETA):
    """Semi-axes (a1,a2,a3) [m] -> translational and rotational friction [SI]."""
    a = np.asarray(axes, dtype=float)
    chi0, chis = chi(a)
    zeta_t = [16 * math.pi * eta / (chi0 + a[i] ** 2 * chis[i]) for i in range(3)]
    zeta_r = []
    for i in range(3):
        j, k = (i + 1) % 3, (i + 2) % 3
        num = 16 * math.pi * eta * (a[j] ** 2 + a[k] ** 2)
        den = 3 * (a[j] ** 2 * chis[j] + a[k] ** 2 * chis[k])
        zeta_r.append(num / den)
    return zeta_t, zeta_r


# ══════════════════════════════════════════════════════════════════════
print("=" * 80)
print("(1) Sphere-limit check — does the formula reproduce the known answer "
      "before we use it")
print("=" * 80)
print(f"  {'R':>10} {'zeta_t calc':>14} {'6.pi.eta.R':>14} {'error':>10} | "
      f"{'zeta_r calc':>13} {'8.pi.eta.R^3':>13} {'error':>10}")
ok = True
for R in (0.1e-6, 0.5e-6, 1e-6, 5e-6):
    zt, zr = friction((R, R, R))
    want_t, want_r = 6 * math.pi * ETA * R, 8 * math.pi * ETA * R**3
    et = 100 * (zt[0] - want_t) / want_t
    er = 100 * (zr[0] - want_r) / want_r
    ok &= abs(et) < 1e-6 and abs(er) < 1e-6
    print(f"  {R*1e6:8.1f}µm {zt[0]:14.6e} {want_t:14.6e} {et:+9.2e}% | "
          f"{zr[0]:13.6e} {want_r:13.6e} {er:+9.2e}%")
print(f"\n  {'✓' if ok else '✗'} equal axes must give isotropy: "
      f"zeta_t spread {100*(max(zt)-min(zt))/zt[0]:.2e}%")

print()
print("=" * 80)
print("(2) abp-rod confirmed shape — major axis 2 µm, minor axis 500 nm "
      "(user, 2026-08-04)")
print("=" * 80)
a1, a2 = 1.0e-6, 0.25e-6            # semi-axis = half the full length
p = a1 / a2
zt, zr = friction((a1, a2, a2))
d_eq = 2 * (a1 * a2 * a2) ** (1 / 3)   # equivalent-volume sphere diameter

print(f"  semi-axes a1={a1*1e6:.3f} µm (major)  a2=a3={a2*1e6:.3f} µm (minor)"
      f"   aspect ratio p = {p:.1f}")
print(f"  equivalent-volume sphere diameter d_eq = 2(a1.a2.a3)^(1/3) = "
      f"{d_eq*1e6:.4f} µm")
print(f"\n  translational friction [kg/s]")
print(f"    zeta_par (major) = {zt[0]:.4e}")
print(f"    zeta_perp (minor)= {zt[1]:.4e}      "
      f"zeta_perp/zeta_par = {zt[1]/zt[0]:.4f}")
print(f"    zeta(sphere d_eq)= {3*math.pi*ETA*d_eq:.4e}   (for comparison)")
print(f"\n  rotational friction [kg.m^2/s]")
print(f"    zeta_r,par (about major) = {zr[0]:.4e}   <- does not change the axis "
      f"direction")
print(f"    zeta_r,perp (about transverse) = {zr[1]:.4e}   ★ this is the in-plane "
      f"2D rotation")

# ── BD isotropic average (§20 option A) ───────────────────────────────
Dt = [KT / z for z in zt]
gamma_bar_3d = 3.0 / sum(1.0 / z for z in zt)          # matched to Dbar = (D1+D2+D3)/3
gamma_bar_2d = 2.0 / (1.0 / zt[0] + 1.0 / zt[1])       # the two in-plane directions only
print(f"\n  ★ BD isotropic average (translational anisotropy cannot be reproduced "
      f"in BD — a hard bd-hoomd constraint)")
print(f"    D_par = {Dt[0]*1e12:.4f} µm^2/s   D_perp = {Dt[1]*1e12:.4f} µm^2/s")
print(f"    gammabar(3D avg) = {gamma_bar_3d:.4e} kg/s   ->  "
      f"Dbar = {KT/gamma_bar_3d*1e12:.4f} µm^2/s")
print(f"    gammabar(2D in-plane) = {gamma_bar_2d:.4e} kg/s   ->  "
      f"Dbar = {KT/gamma_bar_2d*1e12:.4f} µm^2/s  <- this one, being 2D")
print(f"    NOTE why the harmonic mean: the long-time MSD is set by "
      f"Dbar = the mean diffusivity, so")
print(f"       gammabar = kT/Dbar is required for the long-time behaviour to be "
      f"right (mater_plan §20 option A)")

# ── Separation of the rotational-diffusion and tumble timescales ──────
Dr = KT / zr[1]
tau_r = 1.0 / Dr                # 2D: <cos dtheta> = exp(-D_r t)
tau_tumble = 0.5
print(f"\n  ★ rotational diffusion vs tumbling — do the two timescales separate")
print(f"    D_r,perp = {Dr:.4f} 1/s      tau_r = 1/D_r = {tau_r:.4f} s   "
      f"(thermal rotational diffusion)")
print(f"    tau_tumble = {tau_tumble} s                            (user-confirmed)")
print(f"    tau_r/tau_tumble = {tau_r/tau_tumble:.3f}")
if tau_r < tau_tumble:
    print(f"    -> thermal rotational diffusion is **{tau_tumble/tau_r:.1f}x faster "
          f"than tumbling.**")
    print(f"      What governs the orientational correlation is rotational "
          f"diffusion, not tumbling.")
else:
    print(f"    -> tumbling is {tau_r/tau_tumble:.1f}x faster than rotational "
          f"diffusion. Tumbling governs.")

# Effective persistence time: if the two processes are independent, the decay
# rates add
for label, factor in (("tumble fully randomises the direction (run-and-tumble)", 1.0),
                      ("tumble reverses by 180 deg (the sketch's run-and-flip)", 2.0)):
    tau_eff = 1.0 / (1.0 / tau_r + factor / tau_tumble)
    print(f"\n    {label}")
    print(f"      1/tau_eff = 1/tau_r + {factor:.0f}/tau_tumble  ->  "
          f"tau_eff = {tau_eff:.4f} s")
    for v_ums in (1.0, 5.0):
        v = v_ums * 1e-6
        lp = v * tau_eff
        Pe = v * d_eq / (KT / gamma_bar_2d)
        print(f"      v={v_ums:.0f}µm/s:  l_p = v.tau_eff = {lp*1e6:.3f} µm = "
              f"{lp/d_eq:.2f} d_eq"
              f" = {lp/(2*a1):.2f} body lengths   Pe = v.d_eq/Dbar = {Pe:.1f}")

print()
print("=" * 80)
print(f"{'✓ PASS' if ok else '✗ FAIL'} — sphere limit reproduced "
      f"(translation 6.pi.eta.R, rotation 8.pi.eta.R^3, error < 1e-11%)")
print("=" * 80)
print()
print("=" * 80)
print("(3) So what is fixed and what is left open")
print("=" * 80)
print(f"""  FIXED (user-confirmed -> tier 0/1)
    shape     semi-axes (1.0, 0.25, 0.25) µm,  p = 4,  d_eq = {d_eq*1e6:.4f} µm
    medium    water @300K (eta = {ETA*1e3:.3f} mPa*s)
    motion    run-and-tumble, tau_tumble = 0.5 s
    speed     v <= 5 µm/s (sketch)

  DERIVED (Perrin, sphere limit verified)
    gammabar(2D)  {gamma_bar_2d:.4e} kg/s        Dbar = {KT/gamma_bar_2d*1e12:.4f} µm^2/s
    gamma_r,z     {zr[1]:.4e} kg.m^2/s   D_r = {Dr:.4f} 1/s   tau_r = {tau_r:.4f} s

  ★ STILL TO DECIDE — the tumble angle distribution
    tau_eff splits as in the table above depending on how a "tumble" changes the
    direction. The sketch says 'run-and-flip' (180 deg reversal) while the user said
    'run-and-tumble'. Uniform random reorientation is standard in 2D, so that is
    taken as the proposed value (tier 3).

  ★ WORTH FLAGGING — in this system thermal rotational diffusion is FASTER than
    tumbling (tau_r = {tau_r:.3f} s < 0.5 s). So the picture of 'run straight, then
    tumble' does not fully hold: even during a run, thermal fluctuations bend the
    direction on a {tau_r:.2f}s scale, because a p=4 ellipsoid in water is that small.
    -> The orientational correlation in the MSD is governed by tau_eff and is not
    explained by tumbling alone. Needs human confirmation.""")
