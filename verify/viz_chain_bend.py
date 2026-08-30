"""`chain-bend-2d-oscill` results check -- figures plus animation.

**Do not report results in prose alone** (working practice). Two artefacts:

  scratch/_viz/chain_bend_results.png   6 panels -- see the gate results and the
                                        spec errors with your own eyes
  scratch/_viz/chain_bend_motion.gif    chain motion -- kT=0 (the mode shape) vs
                                        kT=1 (why SNR is the problem)

The animation is made with **kT=0 determinism and a large dt** (dt*lambda_max = 0.22,
11% of the stability limit). Its purpose is to show the SHAPE of the motion; it is
**not a production measurement** -- the production dt is 4.53e-10.
The thermal amplitude l_k = sqrt(kT/k_t) is a static equilibrium quantity, so it
still comes out right at this dt (dt*k_t/gamma = 5e-5).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/viz_chain_bend.py
"""
from __future__ import annotations

import glob
import json
import math
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt          # noqa: E402
import numpy as np                       # noqa: E402
from matplotlib.animation import FuncAnimation, PillowWriter   # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "verify" / "_viz"
GATES = ROOT / "verify" / "_gates"
sys.path.insert(0, str(ROOT))

C_MEAS, C_THEORY, C_BAD, C_GOOD, C_GREY = "#1f77b4", "#d62728", "#ff7f0e", "#2ca02c", "#888888"

plt.rcParams["font.family"] = ["DejaVu Sans"]     # * labels in English (CLAUDE.md)
plt.rcParams["axes.unicode_minus"] = False


def load_specs():
    """Read the spec if it exists; otherwise build it from system.yaml directly.

    ★ When an L3 hard check fails, `nondim spec` **writes no spec** (correct
    behaviour). The results figure must not become undrawable because of that, so
    this falls back to the case module.
    """
    sp = [json.loads(Path(q).read_text())
          for q in glob.glob(str(ROOT / "specs" / "chain-bend-2d-oscill__*.json"))]
    if sp:
        return sorted(sp, key=lambda s: s["params"]["omega_star"])
    import argparse as _ap
    from cases import chain_bend_2d as CB
    sys_ = CB.load_system(ROOT / "intake" / "chain-bend-2d-oscill" / "system.yaml")
    args = _ap.Namespace(dt_scale=1.0, cycles=CB.N_CYCLES, samples=2000)
    lo, hi = sys_["omega_range"]
    out = []
    for om in np.geomspace(lo, hi, CB.N_SWEEP):
        _, spec, _, _, _, _ = CB.build_spec(sys_, float(om), args)
        out.append(json.loads(json.dumps(spec.to_doc(), default=str))
                   if hasattr(spec, "to_doc") else
                   dict(params=spec.params, numerics=spec.numerics))
    return sorted(out, key=lambda s: s["params"]["omega_star"])


def chain_matrices(p):
    n = int(p["n_beads"]); trapped = sorted(int(t) for t in p["trapped"])
    kth = float(p["kappa_theta_star"]); kt = float(p["k_t_star"])
    ell = float(p["L_chain_star"]) / (n - 1)
    B = np.zeros((n - 2, n))
    for i in range(n - 2):
        B[i, i], B[i, i + 1], B[i, i + 2] = 1.0, -2.0, 1.0
    B /= ell
    A = kth * (B.T @ B)
    T = np.zeros((n, n))
    for e in trapped:
        T[e, e] += kt
    return n, trapped, kth, kt, ell, A, T


# ════════════════════════════════════════════════════════════════════════
# Figures
# ════════════════════════════════════════════════════════════════════════
def make_figure():
    specs = load_specs()
    p = specs[0]["params"]
    n, trapped, kth, kt, ell, A, T = chain_matrices(p)
    amp = float(p["amp_star"]); mid = trapped[len(trapped) // 2]
    lk = math.sqrt(1.0 / kt)
    em = np.eye(n)[mid]
    ev = np.linalg.eigvalsh(A + T)
    kappa_center = 48 * (kth * ell) / float(p["L_chain_star"]) ** 3
    Kstat = kt * (1.0 / np.linalg.solve(A + T, kt * em)[mid] - 1.0)

    det, deq = {}, {}
    for f in glob.glob(str(GATES / "det_*.json")):
        r = json.load(open(f))
        if r["method"] == "ov":
            det[round(r["de"], 3)] = r
    for f in glob.glob(str(GATES / "deq_*.json")):
        r = json.load(open(f))
        deq[round(r["de"], 3)] = r

    om_all = np.array([s["params"]["omega_star"] for s in specs])
    de_all = np.array([s["params"]["De"] for s in specs])
    snr_all = np.array([s["params"]["snr_response"] for s in specs])

    fig, ax = plt.subplots(2, 3, figsize=(16.5, 9.2))
    fig.suptitle("chain-bend-2d-oscill -- L3 gate results (2026-08-05)  .  "
                 "a record of 5 spec-derivation errors caught by execution",
                 fontsize=13, y=0.98)

    # (1) K*(omega) -- measured vs linear response
    a = ax[0, 0]
    om_f = np.geomspace(om_all.min() * 0.7, om_all.max() * 1.4, 200)
    KA = np.array([kt * amp / np.linalg.solve(1j * w * np.eye(n) + A + T, kt * amp * em)[mid]
                   - kt - 1j * w for w in om_f])
    a.loglog(om_f, KA.real, "-", c=C_THEORY, lw=1.6, label="K' linear response (0.32% vs exact minimisation)")
    a.loglog(om_f, KA.imag, "--", c=C_THEORY, lw=1.6, label="K'' linear response")
    if deq:
        k = sorted(deq)
        a.loglog([deq[d]["omega"] for d in k], [deq[d]["K_re"] for d in k], "o",
                 c=C_MEAS, ms=8, label="K' HOOMD (equilibrated 20 tau_max)")
        a.loglog([deq[d]["omega"] for d in k], [deq[d]["K_im"] for d in k], "s",
                 c=C_MEAS, ms=7, mfc="none", label="K″ HOOMD")
    a.axhline(kappa_center, c=C_GREY, ls=":", lw=1.2)
    a.text(om_f[-1], kappa_center * 1.15, "48EI/L^3 (rigid clamp)", fontsize=8,
           c=C_GREY, ha="right")
    a.axhline(Kstat, c=C_GOOD, ls=":", lw=1.4)
    a.text(om_f[-1], Kstat * 1.15, f"kappa_drive={Kstat:.0f} (trap boundary)",
           fontsize=8,
           c=C_GOOD, ha="right")
    a.set_ylim(Kstat * 0.55, None)
    a.set_xlabel("ω*  [1/τ_B]"); a.set_ylabel("K*  [kT/d²]")
    a.set_title("(1) K*(omega) -- the measurement is 28% low. The cause is HOOMD\n"
                "angle.Harmonic's sin(theta) clamp (force wrong, energy exact)",
                fontsize=10)
    a.legend(fontsize=7.5, loc="upper left"); a.grid(alpha=0.25, which="both")

    # (2) SNR -- what the spec checked vs the real thing
    a = ax[0, 1]
    a.loglog(de_all, snr_all, "o-", c=C_MEAS, ms=7,
             label=r"real  $|\hat{y}(\omega)|/\ell_k$")
    a.axhline(amp / lk, c=C_BAD, ls="-", lw=2,
              label=f"what the spec checked, a/l_k = {amp/lk:.2f} "
                    f"(omega-independent)")
    a.axhline(3, c=C_GOOD, ls="--", lw=1.3, label="check threshold 3")
    a.axhline(1, c=C_GREY, ls=":", lw=1.3, label="equal to the thermal fluctuation")
    for d, s in zip(de_all, snr_all):
        if s < 1:
            a.plot(d, s, "x", c=C_BAD, ms=11, mew=2.2)
    a.set_xlabel("De = ω τ_max"); a.set_ylabel("SNR")
    a.set_title("(2) the SNR check had the wrong numerator -- overestimating by up "
                "to 60x\n"
                "(x = points where the response is below the thermal fluctuation, "
                "4 of 7)", fontsize=10)
    a.legend(fontsize=7.5, loc="lower left"); a.grid(alpha=0.25, which="both")

    # (3) relaxation spectrum + the sweep range
    a = ax[0, 2]
    a.semilogy(range(1, n + 1), ev, "o", c=C_MEAS, ms=5)
    a.axhline(ev[0], c=C_GOOD, ls="--", lw=1.4)
    a.text(1.4, ev[0] * 1.5,
           f"lambda_min={ev[0]:.0f} -> tau_max  ★governing scale", fontsize=8,
           c=C_GOOD)
    a.axhline(kappa_center, c=C_BAD, ls="--", lw=1.4)
    a.text(1.4, kappa_center * 1.5,
           f"kappa_center={kappa_center:.0f} -> tau_chain (what the spec used)",
           fontsize=8, c=C_BAD)
    a.axhline(ev[-1], c=C_GREY, ls=":", lw=1.2)
    a.text(1.4, ev[-1] * 0.35, f"lambda_max={ev[-1]:.2e} -> sets dt", fontsize=8,
           c=C_GREY)
    a.fill_between([1, n], om_all.min(), om_all.max(), color=C_MEAS, alpha=0.13)
    a.text(n * 0.97, om_all.max() * 1.6, "omega sweep range", fontsize=8.5,
           c=C_MEAS, ha="right")
    a.set_xlabel("mode index"); a.set_ylabel(r"eigenvalue $\lambda$  [kT/d$^2$]")
    a.set_title(f"(3) the sweep never reaches the quasi-static region\n"
                f"tau_max/tau_chain = {kappa_center/ev[0]:.2f}x", fontsize=10)
    a.grid(alpha=0.25, which="both")

    # (4) gate A -- the nominal amplitude gets even the sign wrong
    a = ax[1, 0]
    de_A = [0.11, 0.23, 0.49, 1.04, 2.21, 4.70, 10.0]
    Kn = [4805.4, 4840.1, 4823.4, 4777.6, 4076.3, 2106.5, -6559.1]
    Km = [4807.8, 4846.7, 4863.7, 4913.2, 4765.0, 4712.8, 5863.2]
    Ks = [57.8, 62.6, 61.0, 77.9, 108.3, 493.4, 1715.2]
    ks_true = 4830.66
    a.semilogx(de_A, Kn, "o-", c=C_BAD, ms=7,
               label="using the nominal amplitude a -> collapses")
    a.errorbar(de_A, Km, yerr=Ks, fmt="s-", c=C_MEAS, ms=6, capsize=3,
               label=r"using the measured phasor $\hat{y}_c$ -> flat")
    a.axhline(ks_true, c=C_GOOD, ls="--", lw=1.5,
              label=f"analytic k_s = {ks_true:.0f}")
    a.axhline(0, c="k", lw=0.8)
    a.annotate("even the sign is wrong\n(236% error)", xy=(10.0, -6559),
               xytext=(1.5, -5200),
               fontsize=8.5, c=C_BAD, arrowprops=dict(arrowstyle="->", color=C_BAD, lw=1.2))
    a.set_xlabel("De (gate A, single bead)"); a.set_ylabel(r"K'  [kT/d$^2$]")
    a.set_title("(4) zero-order hold on the drive -- using the nominal amplitude\n"
                "in the estimator is silently wrong", fontsize=10)
    a.legend(fontsize=7.5, loc="lower left"); a.grid(alpha=0.25)

    # (5) equilibration -- sigma fell by 1000x
    a = ax[1, 1]
    if det and deq:
        k = sorted(deq)
        x = np.arange(len(k)); w = 0.36
        a.bar(x - w/2, [det[d]["K_sem"] for d in k], w, color=C_BAD,
              label="equilibrated 5 tau_chain = 0.54 tau_max")
        a.bar(x + w/2, [max(deq[d]["K_sem"], 1e-2) for d in k], w, color=C_GOOD,
              label="equilibrated 20 tau_max")
        a.set_yscale("log")
        a.set_xticks(x); a.set_xticklabels([f"{d:g}" for d in k])
        for i, d in enumerate(k):
            r = det[d]["K_sem"] / max(deq[d]["K_sem"], 1e-9)
            a.text(i - w/2, det[d]["K_sem"] * 1.7, f"×{r:,.0f}", ha="center", fontsize=8)
        a.set_ylim(None, max(det[d]["K_sem"] for d in k) * 12)
    a.set_xlabel("De (old definition)"); a.set_ylabel(r"block scatter $\sigma$(K*)  [kT/d$^2$]")
    a.set_title("(5) equilibration at 2.2 tau_max was insufficient\n"
                "switching to tau_max cut the scatter ~1000x", fontsize=10)
    a.legend(fontsize=7.5); a.grid(alpha=0.25, axis="y")

    # (6) gate B' -- the effect of inertia
    a = ax[1, 2]
    rows = []
    for f in glob.glob(str(GATES / "det_*.json")):
        rows.append(json.load(open(f)))
    pair = {}
    for r in rows:
        pair.setdefault(round(r["de"], 3), {})[r["method"]] = r
    dd = sorted([d for d, v in pair.items() if len(v) >= 2])
    if dd:
        dre = [100 * abs(pair[d]["ov"]["K_re"] - pair[d]["lang0"]["K_re"])
               / abs(pair[d]["ov"]["K_re"]) for d in dd]
        dim = [100 * abs(pair[d]["ov"]["K_im"] - pair[d]["lang0"]["K_im"])
               / abs(pair[d]["ov"]["K_im"]) for d in dd]
        a.loglog(dd, dre, "o-", c=C_MEAS, ms=7, label="K' difference")
        a.loglog(dd, dim, "s--", c=C_GOOD, ms=6, label="K'' difference")
    a.axhline(47, c=C_BAD, ls="-", lw=2,
              label="power limit of the thermal comparison (47%)")
    a.axhline(1, c=C_GREY, ls=":", lw=1.2, label="1% reference")
    a.set_xlabel("De (old definition)"); a.set_ylabel("|overdamped - inertial| / K  [%]")
    a.set_title("(6) tau_p/tau_fast=0.60 is harmless -- at most 0.159%\n"
                "the kT=0 deterministic difference is 300x more sensitive than the "
                "thermal comparison", fontsize=10)
    a.legend(fontsize=7.5, loc="upper left"); a.grid(alpha=0.25, which="both")

    fig.tight_layout(rect=[0, 0, 1, 0.955])
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "chain_bend_results.png"
    fig.savefig(out, dpi=125)
    plt.close(fig)
    return out


# ════════════════════════════════════════════════════════════════════════
# Animation -- kT=0 (the mode shape) vs kT=1 (seeing the SNR problem)
# ════════════════════════════════════════════════════════════════════════
def simulate_frames(kT, de_target, n_cycles=3, n_frames=90):
    """Run the ghost-trap chain at a large dt and collect frames.

    For the animation only.
    """
    import gsd.hoomd, hoomd, hoomd.md as md

    specs = load_specs()
    sp = min(specs, key=lambda s: abs(s["params"]["De"] - de_target))
    p = sp["params"]
    n, trapped, kth, kt, ell, A, T = chain_matrices(p)
    amp = float(p["amp_star"]); omega = float(p["omega_star"])
    lam_max = float(np.linalg.eigvalsh(A + T)[-1])
    lam_min = float(np.linalg.eigvalsh(A + T)[0])
    dt = 0.22 / lam_max                      # 11% of the stability limit 2/lambda

    pos = [[(i - (n - 1) / 2) * ell, 0.0, 0.0] for i in range(n)]
    typeid = [0] * n
    for g in trapped:
        pos.append(list(pos[g])); typeid.append(1)
    f = gsd.hoomd.Frame()
    f.particles.N = len(pos)
    f.particles.position = np.array(pos)
    f.particles.typeid = typeid
    f.particles.types = ["A", "G"]
    f.configuration.box = [4 * float(p["L_chain_star"])] * 2 + [0, 0, 0, 0]
    f.configuration.dimensions = 2
    grp = [[i, i + 1] for i in range(n - 1)] + [[g, n + j] for j, g in enumerate(trapped)]
    f.bonds.N = len(grp)
    f.bonds.types = ["backbone", "trap"]
    f.bonds.typeid = [0] * (n - 1) + [1] * len(trapped)
    f.bonds.group = np.array(grp)
    f.angles.N = n - 2
    f.angles.types = ["bend"]
    f.angles.typeid = [0] * (n - 2)
    f.angles.group = np.array([[i, i + 1, i + 2] for i in range(n - 2)])

    sim = hoomd.Simulation(device=hoomd.device.CPU(), seed=7)
    sim.create_state_from_snapshot(f)
    bond = md.bond.Harmonic()
    bond.params["backbone"] = dict(k=float(p["k_bond_star"]), r0=ell)
    bond.params["trap"] = dict(k=kt, r0=0.0)
    angle = md.angle.Harmonic()
    angle.params["bend"] = dict(k=kth, t0=math.pi)
    filt = hoomd.filter.Type(["A"])
    meth = (md.methods.OverdampedViscous(filter=filt, default_gamma=1.0) if kT == 0
            else md.methods.Brownian(filter=filt, kT=kT, default_gamma=1.0))
    integ = md.Integrator(dt=dt, methods=[meth], forces=[bond, angle])
    integ.integrate_rotational_dof = False
    sim.operations.integrator = integ

    mid = trapped[len(trapped) // 2]
    ghost_mid = n + trapped.index(mid)

    class Move(hoomd.custom.Action):
        def act(self, timestep):
            y = amp * math.sin(omega * timestep * dt)
            with self._state.cpu_local_snapshot as s:
                tg = np.array(s.particles.tag, copy=True)
                loc = np.flatnonzero(tg == ghost_mid)
                if len(loc):
                    s.particles.position[loc[0], 1] = y

    sim.operations.updaters.append(hoomd.update.CustomUpdater(
        action=Move(), trigger=hoomd.trigger.Periodic(10)))

    n_eq = int(round(10.0 / (lam_min * dt)))          # 10 τ_max
    sim.run(n_eq)
    period = 2 * math.pi / omega
    total = int(round(n_cycles * period / dt))
    every = max(1, total // n_frames)
    frames, drive = [], []
    for _ in range(n_frames):
        sim.run(every)
        snap = sim.state.get_snapshot()               # a global snapshot is in tag order
        q = np.array(snap.particles.position)
        frames.append(q[:n, :2].copy())
        drive.append(q[ghost_mid, 1])
    return dict(frames=np.array(frames), drive=np.array(drive), n=n, trapped=trapped,
                amp=amp, ell=ell, lk=math.sqrt(1.0 / kt),
                de=float(p["De"]), dt=dt, n_eq=n_eq)


def make_animation():
    hot = simulate_frames(1.0, 9.5)
    cold = simulate_frames(0.0, 9.5)
    nf = len(cold["frames"])
    ylim = max(np.abs(hot["frames"][:, :, 1]).max(), cold["amp"]) * 1.35

    fig, axes = plt.subplots(2, 1, figsize=(9.6, 6.4), sharex=True)
    fig.suptitle(f"chain-bend chain motion  .  De = {cold['de']:.1f}  .  "
                 f"drive amplitude a = {cold['amp']:.3f} d  .  "
                 f"l_k = {cold['lk']:.4f} d",
                 fontsize=11.5)
    arts = []
    for a, dat, ttl, col in ((axes[0], cold,
                              "kT = 0 (deterministic) -- the three-point bending "
                              "mode shape", C_GOOD),
                             (axes[1], hot,
                              "kT = 1 (thermal) -- the response is buried in the "
                              "thermal fluctuation", C_MEAS)):
        a.set_xlim(-cold["ell"] * cold["n"] / 2 * 1.06, cold["ell"] * cold["n"] / 2 * 1.06)
        a.set_ylim(-ylim, ylim)
        a.axhline(0, c=C_GREY, lw=0.7, ls=":")
        a.axhspan(-dat["lk"], dat["lk"], color=C_GREY, alpha=0.18)
        ln, = a.plot([], [], "-", c=col, lw=1.4, zorder=2)
        pt, = a.plot([], [], "o", c=col, ms=6.5, zorder=3)
        tr, = a.plot([], [], "v", c=C_BAD, ms=11, zorder=4)
        a.set_ylabel("y  [d]"); a.set_title(ttl, fontsize=9.5, loc="left")
        a.grid(alpha=0.2)
        arts.append((ln, pt, tr, dat))
    axes[0].text(0.995, 0.05, r"grey band = $\pm\ell_k$ (thermal fluctuation)",
                 transform=axes[0].transAxes,
                 ha="right", fontsize=8, c=C_GREY)
    axes[1].set_xlabel("x  [d]   (▼ = driven trap centre)")

    def upd(i):
        out = []
        for ln, pt, tr, dat in arts:
            q = dat["frames"][i]
            ln.set_data(q[:, 0], q[:, 1]); pt.set_data(q[:, 0], q[:, 1])
            trp = dat["trapped"]
            tr.set_data([q[trp[len(trp) // 2], 0]], [dat["drive"][i]])
            out += [ln, pt, tr]
        return out

    fig.tight_layout(rect=[0, 0, 1, 0.94])
    anim = FuncAnimation(fig, upd, frames=nf, blit=True, interval=55)
    OUT.mkdir(parents=True, exist_ok=True)
    out = OUT / "chain_bend_motion.gif"
    anim.save(out, writer=PillowWriter(fps=18))
    plt.close(fig)
    return out, cold


def main() -> int:
    png = make_figure()
    print(f"figures   -> {png.relative_to(ROOT)}")
    if "--fig-only" in sys.argv:           # do not re-run the simulation when only fixing labels
        return 0
    gif, meta = make_animation()
    print(f"animation -> {gif.relative_to(ROOT)}   "
          f"(De={meta['de']:.1f}, dt={meta['dt']:.2e}, "
          f"equilibrated {meta['n_eq']:,} steps)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
