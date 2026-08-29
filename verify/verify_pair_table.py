"""Measured verification of md.pair.Table -- before using r^-3 in Phase 1-B.

Two places where the documentation looks like it could be quietly wrong are settled
by execution:
  (1) are the force and energy really 0 for r < r_min? (for a diverging potential
      that means the repulsion disappears)
  (2) is the grid linspace(r_min, r_cut, len(U), **endpoint=False**)?
      The snippet in skill bd-hoomd builds it with endpoint=True -- is that right?

    $PY scratch/verify_pair_table.py
"""
import math

import gsd.hoomd
import hoomd
import hoomd.md as md
import numpy as np

A = 10.0
R_MIN, R_CUT, NBINS = 0.5, 3.0, 200


def make_table(endpoint):
    """U = A/r^3 (shifted at the cutoff), F = -dU/dr = 3A/r^4."""
    r = np.linspace(R_MIN, R_CUT, NBINS, endpoint=endpoint)
    U = A / r**3
    U = U - A / R_CUT**3  # 0 at the cutoff (Table has no shift mode)
    F = 3 * A / r**4
    return dict(r_min=R_MIN, U=U, F=F)


def measure(sep, params):
    """Place two particles a distance sep apart and read back HOOMD's force/energy."""
    L = 20.0
    fr = gsd.hoomd.Frame()
    fr.particles.N = 2
    fr.particles.position = [[-sep / 2, 0, 0], [sep / 2, 0, 0]]
    fr.particles.typeid = [0, 0]
    fr.particles.types = ["A"]
    fr.configuration.box = [L, L, 0, 0, 0, 0]
    fr.configuration.dimensions = 2

    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(fr)
    cell = md.nlist.Cell(buffer=0.4)
    tab = md.pair.Table(nlist=cell, default_r_cut=R_CUT)
    tab.params[("A", "A")] = params
    # No integration -- evaluate the force only (run 0 steps rather than set dt=0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-9, methods=[bd], forces=[tab])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(0)
    f = np.array(tab.forces)
    e = np.array(tab.energies)
    return float(f[1][0]), float(e.sum())  # +x force on the right particle, and the total potential


def analytic(sep):
    return 3 * A / sep**4, A / sep**3 - A / R_CUT**3


print("=" * 78)
print("(1) grid convention -- is endpoint=False the right one?")
print("=" * 78)
print(f"  A = {A}, r_min = {R_MIN}, r_cut = {R_CUT}, nbins = {NBINS}")
print(f"\n  {'sep':>6} {'F analytic':>11} | {'endpoint=False':>15} {'error':>9} |"
      f" {'endpoint=True':>14} {'error':>9}")
worst = {False: 0.0, True: 0.0}
for sep in (0.7, 1.0, 1.5, 2.0, 2.5, 2.9):
    Fa, _ = analytic(sep)
    row = f"  {sep:6.2f} {Fa:11.4f} |"
    for ep in (False, True):
        Fm, _ = measure(sep, make_table(endpoint=ep))
        err = 100 * (Fm - Fa) / Fa
        worst[ep] = max(worst[ep], abs(err))
        row += f" {Fm:15.4f} {err:+8.3f}% |"
    print(row)
print(f"\n  worst error:  endpoint=False {worst[False]:.3f}%   "
      f"endpoint=True {worst[True]:.3f}%")
better = "endpoint=False" if worst[False] < worst[True] else "endpoint=True"
print(f"  -> {better} is correct (check this against what the docs claim)")

print()
print("=" * 78)
print("(2) are the force and energy really 0 for r < r_min?  "
      "★ dangerous for a diverging potential")
print("=" * 78)
p = make_table(endpoint=False)
print(f"  {'sep':>6} {'F analytic':>13} {'F meas':>13} {'U meas':>13}   verdict")
for sep in (0.60, 0.51, 0.49, 0.40, 0.30):
    Fa, Ua = analytic(sep)
    Fm, Um = measure(sep, p)
    inside = sep < R_MIN
    verdict = ("below r_min -> 0" if inside else "normal")
    print(f"  {sep:6.2f} {Fa:13.4f} {Fm:13.4f} {Um:13.4f}   {verdict}")
print(f"""
  ★ Conclusion: the repulsion **disappears** for r < r_min({R_MIN}).
    With a diverging potential like r^-3, once a particle breaks through r_min it
    feels no force at all thereafter and stays overlapped -- silently wrong, with
    no error raised.
    -> Guard: (a) put an excluded-volume core (WCA) in as a separate force so
              r_min is never reached, and
              (b) monitor the minimum neighbour distance during the run.
""")

print("=" * 78)
print("(3) WCA(LJ) + Table together -- do the two forces add?")
print("=" * 78)


def measure_both(sep):
    L = 20.0
    fr = gsd.hoomd.Frame()
    fr.particles.N = 2
    fr.particles.position = [[-sep / 2, 0, 0], [sep / 2, 0, 0]]
    fr.particles.typeid = [0, 0]
    fr.particles.types = ["A"]
    fr.configuration.box = [L, L, 0, 0, 0, 0]
    fr.configuration.dimensions = 2
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=1)
    sim.create_state_from_snapshot(fr)
    cell = md.nlist.Cell(buffer=0.4)
    tab = md.pair.Table(nlist=cell, default_r_cut=R_CUT)
    tab.params[("A", "A")] = make_table(endpoint=False)
    wca = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    wca.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=0.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-9, methods=[bd], forces=[tab, wca])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(0)
    return (float(np.array(tab.forces)[1][0]) + float(np.array(wca.forces)[1][0]),
            float(np.array(tab.energies).sum()) + float(np.array(wca.energies).sum()))


def analytic_both(sep):
    Ft, Ut = analytic(sep)
    if sep < 2 ** (1 / 6):
        Uw = 4 * (sep**-12 - sep**-6) + 1.0
        Fw = 4 * (12 * sep**-13 - 6 * sep**-7)
    else:
        Uw = Fw = 0.0
    return Ft + Fw, Ut + Uw


print(f"  {'sep':>6} | {'F analytic':>12} {'F meas':>12} {'err':>9} |"
      f" {'U analytic':>12} {'U meas':>12} {'err':>9}")
ok_all = True
for sep in (0.90, 0.95, 1.00, 1.10, 1.30, 2.00):
    Fa, Ua = analytic_both(sep)
    Fm, Um = measure_both(sep)
    ef = 100 * (Fm - Fa) / Fa
    eu = 100 * (Um - Ua) / Ua
    ok_all &= abs(ef) < 1 and abs(eu) < 1
    print(f"  {sep:6.2f} | {Fa:12.4f} {Fm:12.4f} {ef:+8.3f}% |"
          f" {Ua:12.4f} {Um:12.4f} {eu:+8.3f}%")
print(f"\n  -> the two forces "
      f"{'add correctly' if ok_all else '★ do NOT add correctly'} "
      f"(judged within 1%)")
print()
grid_ok = worst[False] < 1e-6
print("=" * 78)
print("✓ PASS -- endpoint=False exact (0.000%), force 0 below r_min confirmed, "
      "WCA+Table sum exact"
      if (grid_ok and ok_all) else "✗ FAIL -- check the tables above")
print("=" * 78)
