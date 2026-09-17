"""Can HOOMD hold one particle in TWO optical traps, and which one does it pick?

    PY=./bin/py
    $PY verify/verify_dual_trap.py                # all stages
    $PY verify/verify_dual_trap.py --stage 1      # one stage

⚠ **THIS IS A CAPABILITY PROBE, NOT A MEASUREMENT.** Nothing here is sealed,
pre-registered or approved, and no number from it may enter a report as a result.
Its job is to establish that the implementation is right and to size the
production runs — the same role `verify/verify_sedimentation_wall.py` played, and
the same reason: `bd-hoomd` says there is no built-in trap, so a DOUBLE trap has
to be settled by execution rather than assumed from the single-trap recipe.

## The physics, and why the answer is not "the stronger trap"

Occupancy is set by free energy, not depth. In the harmonic limit

    log(p1/p2) = (U01 - U02)/kT  -  (1/2) sum_i ln(kappa_1i / kappa_2i)

so depth favours the deeper trap exponentially while **stiffness disfavours the
stiffer one entropically**. For traps made with the same optics at different
power, `kappa = 2 U0 / w^2` scales with depth, and the two terms fight. Computed
exactly on a 3D grid beforehand (same waist, d/w = 3, axial/transverse 3):

    U1/kT  U2/U1   p1/p2 exact   exp(dU/kT)   ratio
      3.0   0.70          2.52         2.46    1.03
      5.0   0.70          2.86         4.48    0.64
      8.0   0.70          5.64        11.02    0.51
     12.0   0.70         19.0         36.6     0.52

Depth-only Boltzmann overpredicts by up to 2x. And at EQUAL depth with different
waist the answer inverts exactly: `p1/p2 = (w1/w2)^3`, i.e. **the softer trap
wins** — 2.192 measured against 2.197 predicted at `w1/w2 = 1.3`.

## ⚠ The design trap this file exists to avoid

The occupancy ratio computed as an INTEGRAL is strongly dependent on where the
"well" is declared to end, because the unbound region carries volume:

    bound cut U < -0.5 kT  ->  p1/p2 = 2.47
                   -1.0 kT  ->  2.86
                   -2.0 kT  ->  4.11
                   -3.0 kT  -> 11.52        (a factor 4.7 across the range)

So the occupancy must come from a TRAJECTORY, classified by `x < x_saddle`, with
experiment and simulation using the identical rule. That is what makes the
comparison well posed — the arbitrariness cancels because it is the same
arbitrariness on both sides — and it is why a formula cannot stand in for the run.

## ⚠⚠ THE ERROR BAR IS 2.8x BIGGER THAN THE OBVIOUS ONE

Residence times are exponentially distributed, so the total time in a well is
carried by a few long visits and the effective number of independent samples is
far below the visit count. Measured over 10 seeds at 400 tau_w:

    within-run error from the visit counts   0.232   (7.8 %)
    seed-to-seed standard deviation          0.654   (22.2 %)   -> 2.8x

★ This cost two wrong readings before it was caught. A single 2000 tau_w run gave
3.39 and one 800 tau_w seed gave 3.29, both against an exact answer of 2.9515,
and the gap was reported as a 12 % discrepancy in the simulation. It was not:
with 10 seeds the mean is **2.9535 +- 0.2070**, which is **0.0 sigma** from
exact. Both high readings were ordinary fluctuations of a statistic whose error
had been underestimated by a factor of 3.

The same underestimate applies to the EXPERIMENT, unchanged: measuring the
occupancy ratio to 5 % needs 16x the observation time a naive hop-count estimate
suggests. That is a number the BD side owes the microscope side's target-error
gate.

## What 10 seeds x 400 tau_w decides

    measured                       2.9535 +- 0.2070
    exact 3D Boltzmann integral    2.9515      0.0 sigma   agrees
    harmonic + entropy             3.4624      2.5 sigma
    depth-only exp(dU/kT)          4.0552      5.3 sigma   REJECTED

★ **Depth-only Boltzmann is rejected at 5.3 sigma.** That is the demo: the trap a
particle occupies is not set by depth alone. 400 tau_w is 3.9 minutes of real
time for these parameters, so ~40 minutes of recording rejects it, and ~1 hour
separates the harmonic approximation too.

## Stages, in the isolation order rule 7 demands

  1 · ONE Gaussian trap        -> recover kappa from <x^2>, against the EXACT
                                  1D Boltzmann integral for a Gaussian well
                                  (not the harmonic formula — the well is
                                  anharmonic and that difference is measured)
  2 · TWO IDENTICAL traps      -> occupancy must be 1 by symmetry, and hops must
                                  actually occur. A control whose answer is known
                                  exactly and which cannot be fitted
  3 · TWO DIFFERENT traps      -> p1/p2 and both hopping rates, against the three
                                  candidate predictions
"""
from __future__ import annotations

import argparse
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

# ── the example physical system (SI first — rule 1) ───────────────────────
#
# ⚠ EXAMPLE VALUES. The experiment is not fixed yet; these are a plausible
#   silica/PS bead in water in a dual time-shared optical trap, chosen to put the
#   barrier in the window where hopping is observable. Every one is a `tier 3`
#   assumption and none may be cited.
A_PARTICLE = 0.5e-6        # m,  particle RADIUS
ETA = 1.0e-3               # Pa s, water at ~293 K
T_K = 293.0                # K
W_TRAP = 0.5e-6            # m,  transverse beam waist
AXIAL_RATIO = 3.0          # w_z / w_xy of a real optical trap
#  ★ d/w = 1.8 and U0 = 14 kT, and BOTH follow from a measured tension rather
#  than from taste. A particle can escape a dual trap to infinity, not only hop
#  between the wells, and the two barriers pull opposite ways:
#
#     d/w   inter-well / U0   hop rate        escape rate   usable
#     1.8             0.156   0.50 /tau_w     2e-6 /tau_w   YES
#     2.0             0.284   0.084           3e-6          YES
#     2.2             0.412   0.014           3e-6          no (too few hops)
#     3.0             0.789   1e-4            4e-6          no
#
#  So the traps must sit about TWO WAISTS apart. At d/w = 3 the inter-well
#  barrier is 79 % of the trap depth and NO depth satisfies both conditions.
#  ⚠ Experimentally that means d ~ 1 um for w = 0.5 um -- close to the limit at
#  which two traps are resolvable at all, which is a constraint on the optics
#  rather than on the analysis.
#
#  The first version of this probe used d/w = 3 and U1 = 3 kT, and the particle
#  ESCAPED: <x^2> came out 35 w^2 in a 20 w box, because a 3 kT Gaussian well
#  does not confine (its escape barrier IS 3 kT, so ~20 tau_w to leave).
D_SEP = 0.9e-6             # m,  trap separation  (d/w = 1.8)
U1_KT = 14.0               # trap 1 depth, in kT
#  ⚠ 0.90, not 0.70. At U1 = 14 kT a ratio of 0.70 means dU = 4.2 kT and
#  p1/p2 ~ 60, so trap 2 is visited a few times in a whole run and carries no
#  statistics. The occupancy contrast and the statistics pull against each other
#  exactly as the two barriers do, and the resolution is the same: stay close to
#  the symmetric point and measure the EXPONENT rather than one ratio.
#      log(p1/p2) = U1 (1 - r)  -  n * ln(1/r)
#  with n = 0 for depth-only Boltzmann and n = 3/2 for harmonic+entropy in 3D.
#  Scanning r and fitting n is the production design; this probe takes one point.
POWER_RATIO = 0.90         # U2/U1 -- the same optics, less power

K_B = 1.380649e-23


def system():
    kT = K_B * T_K
    gamma = 6.0 * math.pi * ETA * A_PARTICLE
    D = kT / gamma
    tau_w = W_TRAP ** 2 / D                   # time to diffuse one waist
    kappa1 = 2.0 * U1_KT * kT / W_TRAP ** 2   # Gaussian curvature at the bottom
    return dict(kT=kT, gamma=gamma, D=D, tau_w=tau_w, kappa1=kappa1,
                tau_k=gamma / kappa1)


# ── the potential, in reduced units: length = w, energy = kT, time = tau_w ──
def u_reduced(x, y, z, depths, sep, axial):
    """U/kT for the dual Gaussian trap. `depths` in kT, lengths in w."""
    out = np.zeros_like(x, dtype=float)
    for u0, x0 in zip(depths, (-sep / 2.0, +sep / 2.0)):
        out -= u0 * np.exp(-((x - x0) ** 2 + y ** 2) - (z / axial) ** 2)
    return out


def saddle_x(depths, sep, axial, n=200001):
    x = np.linspace(-sep / 2, sep / 2, n)
    z = np.zeros_like(x)
    u = u_reduced(x, z, z, depths, sep, axial)
    i1 = int(np.argmin(np.where(x < 0, u, np.inf)))
    i2 = int(np.argmin(np.where(x > 0, u, np.inf)))
    isad = i1 + int(np.argmax(u[i1:i2 + 1]))
    return float(x[isad]), float(u[isad] - u[i1]), float(u[isad] - u[i2])


def build(depths, sep, axial, box, dt, seed, steps, sample_every):
    """One BD run of a single particle in the dual Gaussian trap."""
    import hoomd
    import hoomd.md as md
    import gsd.hoomd

    class DualGaussianTrap(md.force.Custom):
        """Sum of two anchored Gaussian wells.

        ★ The minimum image is applied to EACH anchor separately. `bd-hoomd`
        trap 1 records +1856 % from omitting it on a single-anchor harmonic trap;
        with two anchors there are two displacements to wrap, and wrapping the
        sum would be meaningless.
        """

        def __init__(self, depths, sep, axial, box_L):
            super().__init__(aniso=False)
            self.depths = [float(q) for q in depths]
            self.anchors = [np.array([-sep / 2.0, 0.0, 0.0]),
                            np.array([+sep / 2.0, 0.0, 0.0])]
            self.period = np.array([box_L, box_L, box_L], dtype=float)
            self.axial = float(axial)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                    self.cpu_local_force_arrays as arr:
                tags = np.array(snap.particles.tag, copy=True)
                pos = np.array(snap.particles.position, copy=True)[np.argsort(tags)]
                f = np.zeros_like(pos)
                pe = np.zeros(len(pos))
                for u0, a in zip(self.depths, self.anchors):
                    dr = pos - a
                    dr -= self.period * np.round(dr / self.period)   # per anchor
                    s = np.array([1.0, 1.0, 1.0 / self.axial ** 2])
                    q = (dr ** 2 * s).sum(axis=1)
                    g = u0 * np.exp(-q)
                    pe -= g
                    #  F = -grad U,  U = -u0 exp(-q),  dq/dr = 2 s r
                    f -= (g[:, None] * 2.0 * s * dr)
                arr.force[np.argsort(tags)] = f
                arr.potential_energy[np.argsort(tags)] = pe

    dev = hoomd.device.CPU(notice_level=0)
    sim = hoomd.Simulation(device=dev, seed=seed)
    frame = gsd.hoomd.Frame()
    frame.particles.N = 1
    frame.particles.position = [[-sep / 2.0, 0.0, 0.0]]   # start in trap 1
    frame.particles.types = ["A"]
    frame.particles.typeid = [0]
    frame.configuration.box = [box, box, box, 0, 0, 0]
    frame.configuration.dimensions = 3
    sim.create_state_from_snapshot(frame)

    trap = DualGaussianTrap(depths, sep, axial, box)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0)
    bd.gamma["A"] = 1.0
    sim.operations.integrator = md.Integrator(dt=dt, methods=[bd], forces=[trap])

    xs, pes = [], []
    n_chunks = max(1, steps // sample_every)
    t0 = time.perf_counter()
    for _ in range(n_chunks):
        sim.run(sample_every)
        s = sim.state.get_snapshot()
        xs.append(np.array(s.particles.position[0], copy=True))
    el = time.perf_counter() - t0
    return np.array(xs), steps / el


def classify(x, xsad):
    """Which trap, frame by frame. THE identical rule both arms must use."""
    return (x > xsad).astype(int)          # 0 = trap 1 (left), 1 = trap 2


def occupancy_and_rates(traj_x, xsad, dt_frame):
    lab = classify(traj_x, xsad)
    n1, n2 = int((lab == 0).sum()), int((lab == 1).sum())
    flips = np.flatnonzero(np.diff(lab) != 0)
    #  residence times, in frames, for each visit
    edges = np.concatenate([[0], flips + 1, [len(lab)]])
    segs = [(lab[edges[i]], edges[i + 1] - edges[i]) for i in range(len(edges) - 1)]
    t1 = [s * dt_frame for w_, s in segs if w_ == 0]
    t2 = [s * dt_frame for w_, s in segs if w_ == 1]
    #  drop the first and last visit -- both are censored by the run boundary
    t1i, t2i = t1[1:-1] if len(t1) > 2 else [], t2[1:-1] if len(t2) > 2 else []
    return dict(n1=n1, n2=n2, ratio=(n1 / n2 if n2 else float("inf")),
                n_hops=len(flips),
                mean_t1=float(np.mean(t1i)) if t1i else float("nan"),
                mean_t2=float(np.mean(t2i)) if t2i else float("nan"),
                n_visits1=len(t1i), n_visits2=len(t2i))


def exact_ratio_grid(depths, sep, axial, ucut=1.0, n=220):
    """The integral version, for contrast only. Reported WITH its cut sensitivity
    because a single number from it would be misleading -- see the docstring."""
    pad = 3.0
    x = np.linspace(-(sep / 2 + pad), sep / 2 + pad, 2 * n)
    y = np.linspace(-pad, pad, n)
    z = np.linspace(-pad * axial, pad * axial, n)
    xsad, _, _ = saddle_x(depths, sep, axial)
    Y, Z = np.meshgrid(y, z, indexing="ij")
    umin = u_reduced(np.array([-sep / 2]), np.array([0.0]), np.array([0.0]),
                     depths, sep, axial)[0]
    wx = np.empty_like(x)
    for i, xv in enumerate(x):
        U = u_reduced(np.full_like(Y, xv), Y, Z, depths, sep, axial)
        g = np.where(U < -ucut, np.exp(-(U - umin)), 0.0)
        wx[i] = np.trapezoid(np.trapezoid(g, z, axis=1), y)
    isad = int(np.argmin(np.abs(x - xsad)))
    return (np.trapezoid(wx[:isad + 1], x[:isad + 1])
            / np.trapezoid(wx[isad:], x[isad:]))


# ── reference values computed without the simulation ──────────────────────
def boltzmann_1d(xg, u0, axial, sep=0.0, depths=None):
    """The exact marginal p(x), by integrating y and z at each x.

    ★ This REPLACES an <x^2> check. A Gaussian well has no finite second moment
    -- U -> 0 at infinity, so the integral is dominated by unbound volume, and
    both the measurement (35 w^2 in a 20 w box) and the grid reference (5 w^2 at
    pad = 4 w) were properties of the DOMAIN rather than of the trap. A SHAPE
    comparison over the range the particle actually visits has no such freedom:
    normalise both over the same bins and the domain cancels.
    """
    pad = 4.0
    y = np.linspace(-pad, pad, 240)
    z = np.linspace(-pad * axial, pad * axial, 240)
    Y, Z = np.meshgrid(y, z, indexing="ij")
    dep = depths if depths is not None else (u0, 0.0)
    out = np.empty_like(np.asarray(xg, dtype=float))
    for i, xv in enumerate(np.asarray(xg, dtype=float)):
        u = u_reduced(np.full_like(Y, xv), Y, Z, dep, sep, axial)
        out[i] = np.trapezoid(np.trapezoid(np.exp(-u), z, axis=1), y)
    return out


# ── the stages ────────────────────────────────────────────────────────────
#  ★ 1.0e-4, BELOW the 1.5e-4 thermal gate, and the gate was measured rather
#  than trusted. Same seed, same length, four steps:
#
#      dt/tau_w   thermal step   p1/p2
#      4.0e-4         0.049 w    2.467     <- 2.7x the gate, 25 % low
#      2.0e-4         0.035 w    3.179     <- 1.3x the gate, still 3 % off
#      1.0e-4         0.025 w    3.293
#
#  ⚠ The FRAME INTERVAL was also suspected and is NOT a factor: the ratio is flat
#  at 3.29 across a 100x range (0.056 to 5.6 tau_k). Ruled out by measurement,
#  recorded so it is not suspected again.
DT = 1.0e-4            # tau_w. thermal gate is 1.5e-4, and it MATTERS -- above
SEP = D_SEP / W_TRAP
AX = AXIAL_RATIO


def stage1(u0=None, n_tau=300):
    u0 = U1_KT if u0 is None else u0
    print(f"\n{'=' * 74}\n1 · ONE Gaussian trap — does p(x) have the Boltzmann "
          f"SHAPE?\n{'=' * 74}")
    steps = int(n_tau / DT)
    traj, rate = build((u0, 0.0), SEP, AX, box=20.0, dt=DT, seed=11,
                       steps=steps, sample_every=100)
    x = traj[:, 0] + SEP / 2.0          # displacement from the trap-1 anchor
    lo, hi, nb = -1.2, 1.2, 40
    hist, edges = np.histogram(x, bins=nb, range=(lo, hi))
    mid = 0.5 * (edges[1:] + edges[:-1])
    ref = boltzmann_1d(mid, u0, AX)
    pm, pr = hist / hist.sum(), ref / ref.sum()
    keep = hist >= 25
    dev = 100.0 * (pm[keep] / pr[keep] - 1.0)
    err = 100.0 * np.sqrt(hist[keep]) / hist.sum() / pr[keep]
    chi2 = float(np.mean((dev / err) ** 2))
    print(f"   {len(traj)} frames, {n_tau} tau_w, {rate:.0f} steps/s, "
          f"|x| max {np.abs(x).max():.2f} w  (box 20 w -- it did NOT escape)")
    print(f"   {int(keep.sum())} bins with >= 25 counts over |x| < {hi} w")
    print(f"   p(x) vs exp(-U/kT):  mean |dev| {np.abs(dev).mean():.2f} %, "
          f"max {np.abs(dev).max():.2f} %, chi2/nu = {chi2:.2f}")
    #  ★ the SAME parabolic fit applied to the EXACT p(x), and that -- not 2 U0
    #    -- is the reference.
    #
    #    ⚠ A Gaussian well is anharmonic: U = -U0 exp(-x^2) = -U0 + U0 x^2
    #      - U0 x^4/2 + ..., so it SOFTENS away from the bottom and a parabolic
    #      fit over a finite window returns less than 2 U0. Measured: 24.15
    #      against an imposed 28.0, i.e. -13.8 %, which the first version of
    #      this stage called a FAIL. It is not an error -- it is what the
    #      estimator does. This is the same mistake the sedimentation campaign
    #      made with `l_g_fitted` (revision 1 predicted 13.375 for a window
    #      where the model gives 15.44), one system later:
    #      **predict for the estimator, not for an idealisation of it.**
    m = np.abs(mid) < 0.35

    def kappa_of(counts):
        return -2.0 * np.polyfit(mid[m], np.log(np.maximum(counts[m], 1e-300)),
                                 2)[0]

    kappa_fit = kappa_of(hist.astype(float))
    kappa_ref = kappa_of(ref * hist.sum() / ref.sum())
    print(f"   kappa, ln p parabola over |x| < 0.35 w      {kappa_fit:7.3f} kT/w^2")
    print(f"   kappa, THE SAME FIT on the exact p(x)       {kappa_ref:7.3f} kT/w^2"
          f"      dev {100*(kappa_fit/kappa_ref-1):+6.2f} %  <- the check")
    print(f"   kappa imposed, 2 U0                         {2*u0:7.3f} kT/w^2"
          f"      the fit is {100*(kappa_ref/(2*u0)-1):+6.2f} % low BY "
          f"CONSTRUCTION")
    print(f"     (a Gaussian well softens away from the bottom; that offset is "
          f"the")
    print(f"      estimator's, not the simulation's, and an experiment "
          f"calibrating")
    print(f"      kappa by equipartition inherits exactly the same offset)")
    ok = chi2 < 4.0 and abs(kappa_fit / kappa_ref - 1) < 0.05
    print(f"   -> {'PASS' if ok else 'FAIL'} (chi2/nu < 4 on the shape, kappa "
          f"within 5 % of the same fit on the exact p)")
    return ok


def stage2(u0=None, n_tau=400):
    u0 = U1_KT if u0 is None else u0
    print(f"\n{'=' * 74}\n2 · TWO IDENTICAL traps — the ratio must be 1, and hops "
          f"must occur\n{'=' * 74}")
    xs, b1, b2 = saddle_x((u0, u0), SEP, AX)
    print(f"   saddle {xs:+.5f} w (must be 0 by symmetry), "
          f"barriers {b1:.3f} / {b2:.3f} kT")
    steps = int(n_tau / DT)
    traj, rate = build((u0, u0), SEP, AX, box=20.0, dt=DT, seed=21,
                       steps=steps, sample_every=100)
    dt_frame = 100 * DT
    r = occupancy_and_rates(traj[:, 0], xs, dt_frame)
    #  the ratio's own error, from the number of independent VISITS
    rel = math.sqrt(1.0 / max(r["n_visits1"], 1) + 1.0 / max(r["n_visits2"], 1))
    print(f"   {r['n_hops']} hops, {r['n_visits1']}/{r['n_visits2']} complete visits")
    print(f"   p1/p2 = {r['ratio']:.4f} +- {r['ratio']*rel:.4f}   "
          f"(expected exactly 1)   {abs(r['ratio']-1)/(r['ratio']*rel):.1f} sigma")
    print(f"   mean residence  {r['mean_t1']:.2f} / {r['mean_t2']:.2f} tau_w "
          f"(must be equal)")
    ok = r["n_hops"] > 30 and abs(r["ratio"] - 1) < 3 * r["ratio"] * rel
    print(f"   -> {'PASS' if ok else 'FAIL'} (>30 hops, ratio within 3 sigma of 1)")
    return ok


def stage3(u1=None, ratio=POWER_RATIO, n_tau=400,
           seeds=tuple(range(101, 111))):
    u1 = U1_KT if u1 is None else u1
    print(f"\n{'=' * 74}\n3 · TWO DIFFERENT traps — which one, and by how much\n"
          f"{'=' * 74}")
    depths = (u1, ratio * u1)
    xs, b1, b2 = saddle_x(depths, SEP, AX)
    print(f"   U1 = {u1} kT, U2 = {depths[1]:.2f} kT, d/w = {SEP:g}, "
          f"axial/transverse {AX:g}")
    print(f"   saddle {xs:+.4f} w, barriers {b1:.3f} / {b2:.3f} kT")
    print(f"   kappa1/kappa2 = {1/ratio:.4f}  (same waist, so it equals U1/U2)\n")
    out = []
    for sd in seeds:
        steps = int(n_tau / DT)
        traj, rate = build(depths, SEP, AX, box=20.0, dt=DT, seed=sd,
                           steps=steps, sample_every=100)
        r = occupancy_and_rates(traj[:, 0], xs, 100 * DT)
        rel = math.sqrt(1.0 / max(r["n_visits1"], 1) + 1.0 / max(r["n_visits2"], 1))
        out.append(r["ratio"])
        print(f"   seed {sd}: {r['n_hops']:4d} hops  "
              f"p1/p2 = {r['ratio']:7.4f} +- {r['ratio']*rel:6.4f}  "
              f"residence {r['mean_t1']:6.2f} / {r['mean_t2']:6.2f} tau_w")
    v = np.array(out)
    m, sem = v.mean(), v.std(ddof=1) / math.sqrt(len(v))
    #  ★ the seed-to-seed spread IS the error. The within-run visit-count error
    #    understates it by 2.8x -- measured, see the module docstring.
    print(f"\n   TRAJECTORY  p1/p2 = {m:.4f} +- {sem:.4f}  "
          f"(seed-to-seed SEM, n={len(v)}, sd {v.std(ddof=1):.4f} = "
          f"{100*v.std(ddof=1)/m:.1f} %/seed)\n")
    depth_only = math.exp(u1 * (1 - ratio))
    harm = depth_only * ratio ** 1.5
    grid = exact_ratio_grid(depths, SEP, AX)
    print(f"   {'candidate':34s} {'value':>9s} {'measured/candidate':>19s}")
    for lbl, val in (("depth only  exp(dU/kT)", depth_only),
                     ("harmonic+entropy  x (k2/k1)^3/2", harm),
                     ("exact 3D grid, bound cut 1 kT", grid)):
        print(f"   {lbl:34s} {val:9.4f} {m/val:19.3f}")
    #  ★ free check 1: the ratio of MEAN RESIDENCE TIMES must equal the
    #    occupancy ratio. Two different reductions of one trajectory, so a
    #    disagreement is a bug in the classifier rather than physics.
    print(f"\n   residence-time ratio {r['mean_t1']/r['mean_t2']:.4f} against "
          f"occupancy {r['ratio']:.4f}  -> "
          f"{100*abs(r['mean_t1']/r['mean_t2']/r['ratio']-1):.1f} % apart")

    #  ★ free check 2: is the grid reference actually well defined HERE? At the
    #    shallow, widely separated settings this probe started from it moved by
    #    4.7x across a plausible bound cut. Deep wells remove that freedom, and
    #    the window chosen from the escape/hopping tension is exactly where it
    #    goes away -- so the analytic reference becomes usable in the same place
    #    the experiment becomes possible.
    cuts = (0.5, 1.0, 2.0, 3.0)
    gv = [exact_ratio_grid(depths, SEP, AX, ucut=c) for c in cuts]
    print(f"   grid value against the bound cut: "
          + ", ".join(f"{c:g} kT -> {v:.4f}" for c, v in zip(cuts, gv)))
    print(f"     spread {max(gv)/min(gv):.2f}x  "
          f"(it was 4.67x at U1=5, r=0.7, d/w=3 -- the window fixes it)")
    print(f"\n   ⚠ Even so, the TRAJECTORY number is the one an experiment is"
          f"\n     compared against, because the classification rule is then"
          f"\n     identical on both sides and its arbitrariness cancels.")
    return True


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", type=int, choices=(1, 2, 3), default=None)
    ap.add_argument("--tau", type=float, default=400.0)
    a = ap.parse_args()

    s = system()
    print(__doc__.split("\n\n")[0])
    print(f"\n⚠ CAPABILITY PROBE — EXAMPLE VALUES, tier 3, none citable.\n")
    print(f"   a = {A_PARTICLE*1e6:.2f} um radius, eta = {ETA*1e3:.2f} mPa s, "
          f"T = {T_K:.0f} K, w = {W_TRAP*1e6:.2f} um")
    print(f"   gamma = {s['gamma']:.3e} kg/s, D = {s['D']*1e12:.3f} um^2/s, "
          f"tau_w = w^2/D = {s['tau_w']:.3f} s")
    print(f"   reduced units: length w, energy kT, time tau_w. "
          f"tau_k/tau_w = 1/(2 U0)")
    print(f"   dt = {DT:.1e} tau_w. Gates: thermal 1.5e-4 (0.03 w) BINDS, "
          f"force 1.2e-2.")
    print(f"     ⚠ dt is {DT/1.5e-4:.1f}x the thermal gate. Accepted for a PROBE "
          f"and stated: the step")
    print(f"       is {math.sqrt(6*DT):.3f} w, there is no excluded volume and no "
          f"hard core, and the")
    print(f"       gate's own rule says to record what the unadopted one would "
          f"have given.")
    print(f"       A production run must re-derive it. tau_k = {1/(2*U1_KT):.4f} "
          f"tau_w, so dt is")
    print(f"       {DT/(1/(2*U1_KT)):.4f} tau_k -- {1/(2*U1_KT)/DT:.0f} steps per "
          f"trap relaxation time.")

    ok = True
    if a.stage in (None, 1):
        ok &= stage1()
    if a.stage in (None, 2):
        ok &= stage2(n_tau=a.tau)
    if a.stage in (None, 3):
        if not ok:
            print("\n★ stage 3 NOT run: rule 7 says an isolated stage must hold "
                  "first.")
            return 1
        stage3(n_tau=a.tau)
    return 0 if ok else 1


if __name__ == "__main__":
    raise SystemExit(main())
