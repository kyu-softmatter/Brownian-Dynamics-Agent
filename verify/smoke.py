"""HOOMD smoke test -- do the APIs the plan depends on actually work?

This checks whether they RUN, not whether they exist (that is survey.py). Each item
corresponds to a specific design decision in the masterplan, and a failure means
that design has to change.
"""
import itertools, math, sys, tempfile, traceback
from pathlib import Path

import numpy as np
import gsd.hoomd
import hoomd
import hoomd.md as md

TMP = Path(tempfile.mkdtemp(prefix="bdsmoke_"))
RESULTS = []


def check(name, plan_ref):
    def deco(fn):
        try:
            detail = fn()
            RESULTS.append(("PASS", name, plan_ref, detail or ""))
        except Exception as e:
            RESULTS.append(("FAIL", name, plan_ref, f"{type(e).__name__}: {e}"))
            if "-v" in sys.argv:
                traceback.print_exc()
        return fn
    return deco


def make_frame(n_side=8, phi=0.4, dim=2, types=("A",), orientation=False):
    N = n_side ** 2
    L = math.sqrt(N * math.pi / (4 * phi))
    a = L / n_side
    pos = np.array([[(i + .5) * a - L / 2, (j + .5) * a - L / 2, 0.]
                    for i, j in itertools.product(range(n_side), repeat=2)])
    f = gsd.hoomd.Frame()
    f.particles.N = N
    f.particles.position = pos
    f.particles.typeid = [0] * N
    f.particles.types = list(types)
    if orientation:
        f.particles.orientation = [(1, 0, 0, 0)] * N
        f.particles.moment_inertia = [(0, 0, 1)] * N
    f.configuration.box = [L, L, 0 if dim == 2 else L, 0, 0, 0]
    f.configuration.dimensions = dim
    return f, N, L


def base_sim(frame, seed=1):
    sim = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=seed)
    sim.create_state_from_snapshot(frame)
    return sim


# ── 1. basic BD + WCA (2D) -- masterplan §11 mapping, appendix A ─────────
@check("BD + WCA 2D", "§11, appendix A")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    sim.run(200)
    return f"N={N} L={L:.2f} step={sim.timestep}"


# ── 2. Tier B: store per-particle forces via the GSD logger -- §9.2 ─────
@check("GSD(logger=) per-particle forces", "§9.2 Tier B")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])

    plog = hoomd.logging.Logger(categories=["particle"])
    plog.add(lj, quantities=["forces", "energies"])
    path = TMP / "tierB.gsd"
    w = hoomd.write.GSD(filename=str(path), trigger=hoomd.trigger.Periodic(50),
                        mode="xb", logger=plog, dynamic=["property", "momentum"])
    sim.operations.writers.append(w)
    sim.run(101)
    w.flush()

    with gsd.hoomd.open(str(path), mode="r") as t:
        fr = t[-1]
        keys = list(fr.log.keys())
        forces = fr.log["particles/md/pair/LJ/forces"]
    return f"frames={len(keys)>0} force_shape={forces.shape} keys={keys}"


# ── 3. Tier C: Burst sliding window -- §9.2 ────────────────────────────
@check("write.Burst sliding window + dump()", "§9.2 Tier C")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])

    path = TMP / "tierC.gsd"
    burst = hoomd.write.Burst(filename=str(path), trigger=hoomd.trigger.Periodic(5),
                              mode="xb", max_burst_size=10,
                              write_at_start=True)
    sim.operations.writers.append(burst)
    sim.run(100)                      # accumulates in the buffer only
    n_buffered = len(burst)
    burst.dump()                      # condition met -> flushed to disk
    burst.flush()
    with gsd.hoomd.open(str(path), mode="r") as t:
        n_written = len(t)
    return f"buffered={n_buffered} written_after_dump={n_written}"


# ── 4. harmonic trap (no built-in -> force.Custom) -- trap-2d-5um ───────
@check("md.force.Custom harmonic trap", "trap-2d-5um / external.harmonic_trap")
def _():
    class HarmonicTrap(md.force.Custom):
        def __init__(self, k, center=(0., 0., 0.)):
            super().__init__(aniso=False)
            self.k = float(k)
            self.center = np.asarray(center, dtype=float)

        def set_forces(self, timestep):
            with self._state.cpu_local_snapshot as snap, \
                 self.cpu_local_force_arrays as arr:
                pos = np.array(snap.particles.position, copy=True)
                d = pos - self.center
                arr.force[:] = -self.k * d
                arr.potential_energy[:] = 0.5 * self.k * (d ** 2).sum(axis=1)

    frame, N, L = make_frame(n_side=4)
    sim = base_sim(frame)
    trap = HarmonicTrap(k=5.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[trap])
    sim.run(300)
    snap = sim.state.get_snapshot()
    r2 = float(np.mean((snap.particles.position[:, :2] ** 2).sum(axis=1)))
    # At equilibrium <r^2> = 2kT/k (2D) = 0.4 -- 300 steps is not converged, so any
    # finite value passes
    return f"<r2>={r2:.4f} (equilibrium prediction {2*1.0/5.0:.3f})"


# ── 5. moving trap: shift the centre with a variant -- trap-drag ────────
@check("variant.Ramp/Cycle time-dependent driving",
       "trap-drag, chain-oscill / driving.*")
def _():
    ramp = hoomd.variant.Ramp(A=0.0, B=5.0, t_start=0, t_ramp=1000)
    cyc = hoomd.variant.Cycle(A=-1.0, B=1.0, t_start=0, t_A=10, t_AB=100, t_B=10, t_BA=100)
    vals = [ramp(0), ramp(500), ramp(1000), cyc(0), cyc(60), cyc(160)]
    return "ramp(0,500,1000)=%.2f,%.2f,%.2f cycle(0,60,160)=%.2f,%.2f,%.2f" % tuple(vals)


# ── 6. chain: bond + angle -- chain-bend-2d-oscill ─────────────────────
@check("bond.Harmonic + angle.Harmonic (chain)", "chain-bend / bonded.*")
def _():
    M = 20
    f = gsd.hoomd.Frame()
    f.particles.N = M
    f.particles.position = [[i * 1.0 - M / 2, 0., 0.] for i in range(M)]
    f.particles.typeid = [0] * M
    f.particles.types = ["A"]
    f.configuration.box = [60, 60, 0, 0, 0, 0]
    f.configuration.dimensions = 2
    f.bonds.N = M - 1
    f.bonds.types = ["backbone"]
    f.bonds.typeid = [0] * (M - 1)
    f.bonds.group = [[i, i + 1] for i in range(M - 1)]
    f.angles.N = M - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (M - 2)
    f.angles.group = [[i, i + 1, i + 2] for i in range(M - 2)]

    sim = base_sim(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=100.0, r0=1.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=10.0, t0=math.pi)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-5, methods=[bd], forces=[bond, angle])
    sim.run(200)
    return f"M={M} bonds={f.bonds.N} angles={f.angles.N} ok"


# ── 7. ABP: Active + ActiveRotationalDiffusion — active.abp ────────────
@check("force.Active + create_diffusion_updater", "§11 trap 3 / active.abp")
def _():
    frame, N, L = make_frame(orientation=True)
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    active = md.force.Active(filter=hoomd.filter.All())
    active.active_force["A"] = (10.0, 0.0, 0.0)
    active.active_torque["A"] = (0.0, 0.0, 0.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-5, methods=[bd], forces=[lj, active])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    upd = active.create_diffusion_updater(trigger=hoomd.trigger.Periodic(1),
                                          rotational_diffusion=3.0)
    sim.operations.updaters.append(upd)
    sim.run(200)
    return f"updater={type(upd).__name__} step={sim.timestep}"


# ── 8. discrete activity: run-and-flip via a custom updater -- abp-rod ──
@check("custom Action updater (run-and-flip)", "abp-rod / active.run_and_flip")
def _():
    class RunAndFlip(hoomd.custom.Action):
        """Reverse the orientation by 180 deg as a Poisson process.

        Different from the ABP's continuous rotational diffusion.
        """
        def __init__(self, rate, dt, seed=7):
            self.p = rate * dt
            self.rng = np.random.default_rng(seed)
            self.n_flips = 0

        def act(self, timestep):
            with self._state.cpu_local_snapshot as snap:
                q = np.array(snap.particles.orientation, copy=True)
                flip = self.rng.random(len(q)) < self.p
                self.n_flips += int(flip.sum())
                # multiply by the 180 deg z-axis rotation quaternion (0,0,0,1)
                # -> the orientation reverses
                q[flip] = np.column_stack([-q[flip, 3], -q[flip, 2],
                                            q[flip, 1],  q[flip, 0]])
                snap.particles.orientation[:] = q

    frame, N, L = make_frame(n_side=6, orientation=True)
    sim = base_sim(frame)
    active = md.force.Active(filter=hoomd.filter.All())
    active.active_force["A"] = (5.0, 0.0, 0.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    integ = md.Integrator(dt=1e-4, methods=[bd], forces=[active])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ
    action = RunAndFlip(rate=2.0, dt=1e-4)
    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=action, trigger=hoomd.trigger.Periodic(1)))
    sim.run(500)
    return f"flips={action.n_flips} (expected ~{500*2.0*1e-4*N:.1f})"


# ── 9. arbitrary soft potential r^-3 -- soft-r3-2d-A-sweep ─────────────
@check("pair.Table (arbitrary r^-n potential)", "soft-r3 / pair.table")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    r_min, r_cut, nbins = 0.3, 3.0, 200
    r = np.linspace(r_min, r_cut, nbins)
    A = 1.0
    U = A / r ** 3
    F = 3 * A / r ** 4
    U = U - U[-1]                                   # shifted to 0 at the cutoff
    tab = md.pair.Table(nlist=cell, default_r_cut=r_cut)
    tab.params[("A", "A")] = dict(r_min=r_min, U=U, F=F)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-5, methods=[bd], forces=[tab])
    sim.run(200)
    return f"nbins={nbins} r=[{r_min},{r_cut}] U(r_min)={U[0]:.2f}"


# ── 10. non-spherical: GayBerne ellipsoid + per-axis rotational drag -- abp-rod ──
@check("pair.aniso.GayBerne + per-axis gamma_r (ellipsoid)",
       "abp-rod / shape.ellipsoid")
def _():
    frame, N, L = make_frame(n_side=6, phi=0.2, orientation=True)
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    gb = md.pair.aniso.GayBerne(nlist=cell, default_r_cut=4.0)
    gb.params[("A", "A")] = dict(epsilon=1.0, lperp=0.5, lpar=1.5)
    # Perrin: the per-axis rotational friction can be given as a tuple
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0,
                             default_gamma=1.0, default_gamma_r=(0.5, 0.5, 2.0))
    integ = md.Integrator(dt=1e-5, methods=[bd], forces=[gb])
    integ.integrate_rotational_dof = True
    sim.operations.integrator = integ
    sim.run(200)
    return f"lperp=0.5 lpar=1.5 gamma_r=(0.5,0.5,2.0) aspect={1.5/0.5:.1f}"


# ── 11. overdamped with no thermal noise (for golden tests) -- §17 ──────
@check("methods.OverdampedViscous (deterministic overdamped)",
       "§17 golden test")
def _():
    frame, N, L = make_frame(n_side=4)
    sim = base_sim(frame)
    ov = md.methods.OverdampedViscous(filter=hoomd.filter.All(), default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-3, methods=[ov], forces=[])
    sim.run(100)
    return "usable for deterministic relaxation checks (tau=gamma/k)"


# ── 12. HDF5Log — §9.2 Tier L ──────────────────────────────────────────
@check("write.HDF5Log global scalar logging", "§9.2 Tier L")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
    thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)
    log = hoomd.logging.Logger(categories=["scalar", "sequence"])
    log.add(thermo, quantities=["potential_energy", "pressure", "kinetic_temperature"])
    path = TMP / "log.h5"
    h5 = hoomd.write.HDF5Log(trigger=hoomd.trigger.Periodic(20),
                             filename=str(path), logger=log, mode="x")
    sim.operations.writers.append(h5)
    sim.run(101)
    h5.flush()
    import h5py
    with h5py.File(path, "r") as fh:
        keys = []
        fh.visit(lambda k: keys.append(k))
    return f"datasets={[k for k in keys if 'potential' in k or 'pressure' in k]}"


# ── 13. restart GSD -- §14 ─────────────────────────────────────────────
@check("restart GSD (truncate)", "§14 checkpointing")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
    path = TMP / "restart.gsd"
    w = hoomd.write.GSD(filename=str(path), trigger=hoomd.trigger.Periodic(50),
                        mode="wb", truncate=True)
    sim.operations.writers.append(w)
    sim.run(101)
    w.flush()
    sim2 = hoomd.Simulation(device=hoomd.device.CPU(notice_level=0), seed=2)
    sim2.create_state_from_gsd(filename=str(path))
    with gsd.hoomd.open(str(path), mode="r") as t:
        nframes = len(t)
    return f"restart_frames={nframes} (truncate -> must be 1) reload_ok"


# ── 14. custom Action for the runtime guard -- §12.5 ───────────────────
@check("custom Action runtime monitoring (NaN / energy)", "§12.5 runtime guard")
def _():
    class Guard(hoomd.custom.Action):
        def __init__(self, thermo):
            self.thermo = thermo
            self.checks = 0
            self.pe_history = []

        def act(self, timestep):
            self.checks += 1
            pe = self.thermo.potential_energy
            self.pe_history.append(pe)
            if pe is not None and not math.isfinite(pe):
                raise RuntimeError(f"non-finite PE at step {timestep}")

    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
    thermo = md.compute.ThermodynamicQuantities(filter=hoomd.filter.All())
    sim.operations.computes.append(thermo)
    g = Guard(thermo)
    sim.operations.writers.append(hoomd.write.CustomWriter(
        action=g, trigger=hoomd.trigger.Periodic(25)))
    sim.run(101)
    return f"checks={g.checks} last_PE={g.pe_history[-1]:.4f}"


# ── 15. Tier D: high frequency for tracer particles only -- §9.2 ───────
@check("GSD(filter=Tags) tracer particle subset", "§9.2 Tier D")
def _():
    frame, N, L = make_frame()
    sim = base_sim(frame)
    cell = md.nlist.Cell(buffer=0.4)
    lj = md.pair.LJ(nlist=cell, default_r_cut=2 ** (1 / 6), mode="shift")
    lj.params[("A", "A")] = dict(epsilon=1.0, sigma=1.0)
    bd = md.methods.Brownian(filter=hoomd.filter.All(), kT=1.0, default_gamma=1.0)
    sim.operations.integrator = md.Integrator(dt=1e-4, methods=[bd], forces=[lj])
    tracers = hoomd.filter.Tags(list(range(8)))
    path = TMP / "tracers.gsd"
    w = hoomd.write.GSD(filename=str(path), trigger=hoomd.trigger.Periodic(10),
                        mode="xb", filter=tracers, dynamic=["property"])
    sim.operations.writers.append(w)
    sim.run(101)
    w.flush()
    with gsd.hoomd.open(str(path), mode="r") as t:
        n_part = t[-1].particles.N
        n_fr = len(t)
    return f"tracer_N={n_part} (of {N} total) frames={n_fr}"


# ── results ─────────────────────────────────────────────────────────────
print("=" * 92)
print(f"HOOMD {hoomd.version.version} SMOKE TEST   (tmp: {TMP})")
print("=" * 92)
w1 = max(len(r[1]) for r in RESULTS)
for status, name, ref, detail in RESULTS:
    mark = "✓" if status == "PASS" else "✗"
    print(f" {mark} {name:<{w1}}  [{ref}]")
    if detail:
        print(f"     {detail}")
n_pass = sum(1 for r in RESULTS if r[0] == "PASS")
print("=" * 92)
print(f"{n_pass}/{len(RESULTS)} PASS")
sys.exit(0 if n_pass == len(RESULTS) else 1)
