"""`trap-drag-2d-hex300` — L3 non-dimensionalization + L4 builder.

(The Korean original said "no run yet". That was stale as of 2026-08-29: measured
84 specs carry this label and 80 of their run directories hold a completed
`metrics.json`. Corrected rather than translated forward. Note those runs predate
`observables.npz` being written, so `make_plots` cannot replot them.)

A hexagonal colloidal lattice (soft repulsion A/r³ + WCA) driven through at constant
speed by a single optical trap. Active microrheology / dragged probe.

**This is the union of 1-A and 1-B** — the trap (1-A) and the soft pair lattice (1-B)
coexist in one system, so there are **two stiffnesses**. Three things are new here:

  ① `dt` is set by **whichever of the two stiffnesses is faster** — τ_k = γ/k_t is
     218× faster than τ_int = γ/U''(r_min). 1-A had only the trap and 1-B only the
     pair, so there was no competition.
  ② An **advective timescale** τ_v = d/v_x appears for the first time (a moving
     boundary condition).
  ③ ★ **Measurability shows up as a statistics problem, not as a hard gate** — the
     drag lag Δr_ss = γv/k_t = 2.0 nm is 1/10 of the in-trap thermal fluctuation
     ℓ_k = 20.4 nm (SNR = 0.0985).
     Every separation check passes, yet **one traverse does not reach the precision
     we want.**
     That is exactly the kind of flaw L3 has to catch, and it is visible before running.

L4 (the run) is driven by `bdbot.run` and judges **numerical health only**
(blow-up, NaN, frozen, drift).
The physics comparison is the role system in `metrics.observables` — not L4's job.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY cases/trap_drag_2d.py --report        # the L3 report only
    $PY cases/trap_drag_2d.py --spec          # the L3 spec
    $PY cases/trap_drag_2d.py --smoke         # L4 wiring check (short)
    $PY cases/trap_drag_2d.py                 # the real L4 run (~6.8e6 steps)
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np
import yaml

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from bdbot import checks as C, materials as M, metrics as MET, report as R  # noqa: E402
from bdbot import nondim as ND, run as RUN, scales as SC, sim as SIM  # noqa: E402
from bdbot import stats as ST, traps as TR  # noqa: E402
from bdbot.pairpot import HEX_NN, U2_star, U_star, a_mean_star, approach_distance  # noqa: E402
from bdbot.provenance import load_node  # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
R_TABLE_MIN = 0.5          # pair.Table lower bound (same as 1-B). The basis of the trap-11 guard
N_CYCLE_TARGET = 100.0     # statistics-check basis -- how many lattice periods to cross
STAT_TARGET_PCT = 2.0      # ⟨F_drag⟩ target precision (system.yaml numerics.stat_target_pct)
# ★ GSD write period [steps]. None means 1/300 of the total. **Not put in the spec** —
#   it is an output setting, not physics, so putting it in the hash would give the same
#   system a different run_id.
#   A fast v has a short drag phase (12,544 steps at v=32), so the default period would
#   give only 2 frames.
GSD_EVERY: int | None = None


# ════════════════════════════════════════════════════════════════════════
# ① The physical system (SI)
# ════════════════════════════════════════════════════════════════════════
def load_system(path: Path, v_sel: float | None = None) -> dict:
    """`v_sel` is the value to pick from the `external.drag_velocity` list [µm/s].

    None means the sketch value.
    """
    raw = yaml.safe_load(path.read_text())
    P = load_node
    # ★ Velocity sweep — pick one point from the list. Invent nothing; use **only values
    #   that are in L2**.
    vnode = dict(raw["external"]["drag_velocity"])
    vlist = vnode["value"] if isinstance(vnode["value"], list) else [vnode["value"]]
    vlist = [float(x) for x in vlist]
    vsel = float(v_sel) if v_sel is not None else (0.5 if 0.5 in vlist else vlist[0])
    if not any(abs(vsel - x) < 1e-12 for x in vlist):
        raise ValueError(f"v={vsel}  is not in the L2 drag_velocity list: {vlist}. "
                         "A spec only ever uses what was derived from the physical "
                         "system (rule 2) — "
                         "to use a new velocity, add it to system.yaml first.")
    vnode["value"] = vsel
    r3, wca = raw["interactions"][0], raw["interactions"][1]
    return {
        "label": raw["label"],
        "dim": raw["dimensions"],
        "d": P(raw["particle"]["diameter"]),
        "rho_p": P(raw["particle"]["density"]),
        "N": int(raw["particle"]["count"]["value"]),
        "T": P(raw["medium"]["temperature"]),
        "eta": P(raw["medium"]["viscosity"]),
        "phi": float(raw["geometry"]["area_fraction"]["value"]),
        # ★ Cell counts of the commensurate hexagon. N = n_x·n_y, so N is read
        #   separately and cross-checked (build_ledger below).
        "n_x": int(raw["geometry"]["n_x"]["value"]),
        "n_y": int(raw["geometry"]["n_y"]["value"]),
        # The box as **written down** by L2. It is compared below against the value
        # derived from the lattice — just trusting the derived value means nobody
        # notices when the YAML silently drifts (same attitude as physical.verify).
        "box_x": load_node(raw["geometry"]["box_length_x"]),
        "box_y": load_node(raw["geometry"]["box_length_y"]),
        "A": float(r3["amplitude_A"]["value"]),
        "r_c": P(r3["cutoff"]),
        "wca_eps_kT": float(wca["epsilon"]["value"]),
        "k_t": P(raw["external"]["stiffness"]),
        "v_x": P(vnode), "v_list": vlist,
        "n_trapped": int(raw["external"]["n_trapped"]["value"]),
        "numerics": raw["numerics"],
        "targets": [t["name"] for t in raw["targets"]],
        "_raw": raw,
    }


# ════════════════════════════════════════════════════════════════════════
# ② The scale ledger (bd-physics §0 ①②)
#    1-A (3 lengths, 5 times) and 1-B (5 lengths, 5 times) combine into 7 lengths,
#    7 times.
# ════════════════════════════════════════════════════════════════════════
def build_ledger(sys_, *, dt_scale=1.0, n_traverse=1.0,
                 warm=10.0, equil=20.0, relax=40.0) -> SC.ScaleLedger:
    d = sys_["d"].value.to("m")
    b = M.sphere_bulk(d, sys_["T"].value, sys_["eta"].value, sys_["rho_p"].value)
    kT, gamma, tau_B = b["kT"], b["gamma"], b["tau_B"]
    k_t = sys_["k_t"].value.to("N/m")
    v_x = sys_["v_x"].value.to("m/s")
    A, eps, phi, N = sys_["A"], sys_["wca_eps_kT"], sys_["phi"], sys_["N"]

    # ── ★ Commensurate hexagonal lattice (the box is NOT square) ──────────
    #   L_x = n_x·a_NN , L_y = n_y·(√3/2)·a_NN , N = n_x·n_y , n_y even.
    #   It used to be L = a_mean√N (square), but then the initial lattice does not match
    #   at the seam — and this system's observables (lattice deformation field, ψ₆) are
    #   sensitive to exactly that defect.
    n_x, n_y = sys_["n_x"], sys_["n_y"]
    if n_x * n_y != N:
        raise ValueError(f"commensurate hexagon broken: n_x·n_y = {n_x}·{n_y} = {n_x*n_y} ≠ N = {N}")
    if n_y % 2:
        raise ValueError(f"n_y = {n_y} is odd — the staggered rows have a period of 2 rows, so "
                         "they do not join across the periodic boundary")
    a_star = a_mean_star(phi)                       # a_mean/d
    a_nn_star = HEX_NN * a_star                     # a_NN/d — this is the lattice constant
    Lx_star = n_x * a_nn_star
    Ly_star = n_y * (math.sqrt(3) / 2) * a_nn_star
    # φ is conserved **identically** in this construction (φ = πd²/(4a_mean²) and
    #   L_x·L_y = N a_mean²). So recomputing φ here and comparing always passes — it is
    #   a check that never fires.
    #   The comparison that actually means something is against **the box L2 wrote
    #   down**: L_x·L_y is derived from φ·n_x·n_y, so a silent drift in the YAML's
    #   box_length_* can only be caught here.
    for sym, got_star, node in (("L_x", Lx_star, sys_["box_x"]), ("L_y", Ly_star, sys_["box_y"])):
        want = float((node.value.to("m") / d).to("dimensionless").magnitude)
        if abs(got_star - want) > 1e-3 * want:
            raise ValueError(
                f"{sym} disagrees with the L2 declaration: lattice-derived {got_star*float(d.to('um').magnitude):.4f} µm "
                f"vs system.yaml {node.value.to('um'):~.4gP} "
                f"({100*(got_star-want)/want:+.3f}%). Fix φ·n_x·n_y and the box together.")
    r_c_star = float((sys_["r_c"].value.to("m") / d).to("dimensionless").magnitude)
    r_min_star, crit, u_rms_rel, state = approach_distance(A, a_star, eps)

    # ★ Two stiffnesses → two relaxation times. Both are γ/(local stiffness)
    #   (bdbot.checks).
    tau_k = C.relaxation_time(gamma, k_t)                                  # trap
    tau_int = (tau_B / float(U2_star(r_min_star, A, eps))).to("s")          # pair (in kT/d² units)
    tau_v = (d / v_x).to("s")                                              # advection (new here)

    # dt is set by the **fastest** scale — here that is the trap (218× faster than the
    # pair).
    tau_dt = min((tau_k, tau_int), key=lambda q: q.to("s").magnitude)
    dt = dt_scale * C.dt_from_gate(tau_dt)
    # ★ The traverse is along L_x, the **drag direction** (not the short side L_y).
    #   `tau_cross` is a property of the system, **independent of the protocol**. It
    #   used not to be separated out, and the geometry check (wake healing) used `T_obs`
    #   as the denominator — so changing `--traverse` silently changed what the check
    #   meant (a hard FAIL while trying to run an equilibration study at traverse=0.001).
    #   A check that measures a property of the system must look only at the system's
    #   own scales.
    tau_cross = (Lx_star * d / v_x).to("s")

    # Four phases (unit: τ_int; the drag phase alone uses the traverse time).
    # T_obs is the **whole protocol**.
    n_of = lambda x: int(round(x * float((tau_int / dt).to(""))))
    n_warm, n_equil, n_relax = n_of(warm), n_of(equil), n_of(relax)
    n_drag = int(round(float((n_traverse * tau_cross / dt).to(""))))
    T_obs = ((n_warm + n_equil + n_drag + n_relax) * dt).to("s")

    l_k = (kT / k_t) ** 0.5                            # in-trap thermal fluctuation width
    dr_ss = (gamma * v_x / k_t).to("m")                # steady-state drag lag

    lg = SC.ScaleLedger()
    lg.add_length("dr_ss", dr_ss, "gamma*v/k_t drag lag (signal)", star=True)
    lg.add_length("l_k", l_k.to("m"), "sqrt(kT/k_t) trap fluctuation (noise)", star=True)
    lg.add_length("d", d, "particle diameter (reference)")
    lg.add_length("r_min", r_min_star * d, "closest approach distance")
    lg.add_length("a_mean", a_star * d, "mean spacing")
    lg.add_length("a_NN", a_nn_star * d, "hex nearest-neighbour = lattice constant (NOT a_mean)")
    lg.add_length("r_c", r_c_star * d, "cutoff")
    # ★ The box is rectangular, so there are two sides. Minimum image is set by the
    #   **short side**, so the `box` role goes there (a role is a function, not a
    #   symbol — bdbot.scales).
    #   The long side is the drag direction and sets T_obs and the lattice-period count.
    short, long_ = (("L_y", Ly_star), ("L_x", Lx_star)) if Ly_star <= Lx_star \
        else (("L_x", Lx_star), ("L_y", Ly_star))
    lg.add_length(short[0], short[1] * d, f"box short side — denominator of minimum image ({n_y} rows)"
                  if short[0] == "L_y" else f"box short side — denominator of minimum image ({n_x} cols)",
                  role="box")
    lg.add_length(long_[0], long_[1] * d,
                  f"box long side = drag direction ({n_x} cols x a_NN)"
                  if long_[0] == "L_x" else f"box long side ({n_y} rows)")
    lg.add_time("tau_p", b["tau_p"], "m/gamma momentum relaxation", role="inertia")
    lg.add_time("dt", dt, "integration step", role="dt")
    lg.add_time("tau_k", tau_k, "gamma/k_t trap — fastest, sets dt", star=True)
    lg.add_time("tau_int", tau_int, f"γ/U''(r_min={r_min_star:.3f}d) pair")
    lg.add_time("tau_v", tau_v, "d/v_x advection (moving trap)")
    lg.add_time("tau_cell", (HEX_NN * a_star * d / v_x).to("s"), "a_NN/v_x lattice period")
    lg.add_time("tau_cross", tau_cross, "L_x/v_x box traverse (a system property — protocol-independent)")
    lg.add_time("tau_B", tau_B, "d^2/D_t diffusion (reference)")
    lg.add_time("T_obs", T_obs, "observation window = all four phases", role="observation")
    lg.add_energy("kT", kT, "thermal energy (reference)")
    lg.add_energy("U_a", (float(U_star(a_star, A, eps)) * kT).to("J"),
                  "U(a_mean) mean-spacing coupling = Gamma*kT")
    # ★ U(d) = (A+ε) kT — NOT A kT (a label error the L3 check caught in 1-B)
    lg.add_energy("U_d", (float(U_star(1.0, A, eps)) * kT).to("J"),
                  "U(d) contact coupling = (A+eps_WCA)*kT")
    lg.add_energy("k_t_d2", (k_t * d**2).to("J"), "k_t*d^2 trap stiffness")

    lg.derived = dict(gamma=gamma, D_t=b["D_t"], m=b["m"], kT=kT, d=d, tau_B=tau_B,
                      tau_k=tau_k, tau_int=tau_int, tau_v=tau_v, dt=dt, T_obs=T_obs,
                      tau_cross=tau_cross, n_warm=n_warm, n_equil=n_equil,
                      n_drag=n_drag, n_relax=n_relax,
                      tau_int_steps=float((tau_int / dt).to("")),
                      a_star=a_star, a_nn_star=a_nn_star, n_x=n_x, n_y=n_y,
                      Lx_star=Lx_star, Ly_star=Ly_star, r_c_star=r_c_star,
                      r_min_star=r_min_star, crit=crit, u_rms_rel=u_rms_rel, state=state,
                      l_k=l_k.to("m"), dr_ss=dr_ss, k_t=k_t, v_x=v_x,
                      tau_dt_name=("tau_k" if tau_dt is tau_k else "tau_int"))
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ This system has **two** stiffnesses (trap k_t, "
        "pair U'') and "
        "the trap is 218x faster -> the trap sets dt. 1-A had only the trap, "
        "1-B only the pair.")
    lg.rationale = lg.ref["rationale"]
    return lg


# ════════════════════════════════════════════════════════════════════════
# ③ Dimensionless groups + ④ separation checks (bd-physics §3, §4)
# ════════════════════════════════════════════════════════════════════════
def analyze_scales(sys_, lg):
    D = lg.derived
    f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)
    g = lg.get
    A, eps, phi = sys_["A"], sys_["wca_eps_kT"], sys_["phi"]
    a_star, r_c_star = D["a_star"], D["r_c_star"]
    Lx_star, Ly_star = D["Lx_star"], D["Ly_star"]
    L_min = min(Lx_star, Ly_star)
    Gamma = float(U_star(a_star, A, eps))
    snr = lg.ratio("lengths", "dr_ss", "l_k")

    # Relative precision of ⟨Δr_ss⟩ obtained from one traverse.
    #   independent samples n ≈ T_obs/(2τ_k) (displacement correlation time τ_k),
    #   per-sample noise ℓ_k
    #   ⟹ SEM/Δr_ss = 1/(SNR·√n)
    n_indep = 0.5 * lg.ratio("times", "tau_cross", "tau_k")   # per one drag traverse
    prec_pct = 100.0 / (snr * math.sqrt(n_indep))
    n_cell = lg.ratio("times", "tau_cross", "tau_cell")   # lattice periods in one traverse = n_x

    groups = [
        ND.Group("Gamma", Gamma, ("energies", "U_a"), ("energies", "kT"),
                 "A(d/a_mean)³", "pair coupling -> hexagonal crystal ★"),
        ND.Group("A", A, None, None, "", "r^-3 amplitude (input, ★proposed)"),
        ND.Group("U(d)/kT", float(U_star(1.0, A, eps)),
                 ("energies", "U_d"), ("energies", "kT"), "(A+ε_WCA)", "contact coupling"),
        ND.Group("k*", lg.ratio("energies", "k_t_d2", "kT"),
                 ("energies", "k_t_d2"), ("energies", "kT"), "k_t d²/kT",
                 "trap vs thermal fluctuation (very stiff)"),
        ND.Group("phi", phi, None, None, "", "packing"),
        ND.Group("Pe_drag", lg.ratio("times", "tau_B", "tau_v"),
                 ("times", "tau_B"), ("times", "tau_v"), "v_x d/D_t", "advection vs diffusion"),
        ND.Group("SNR", snr, ("lengths", "dr_ss"), ("lengths", "l_k"), "Δr_ss/ℓ_k",
                 "★ the drag signal is 1/10 of the thermal fluctuation"),
        ND.Group("a_mean/d", a_star, ("lengths", "a_mean"), ("lengths", "d"), "", "mean spacing"),
        ND.Group("a_NN/a_mean", HEX_NN, ("lengths", "a_NN"), ("lengths", "a_mean"),
                 "√(2/√3)", "hex nearest-neighbour (a parameter-free prediction)"),
        ND.Group("L_x/d", Lx_star, ("lengths", "L_x"), ("lengths", "d"), "n_x·a_NN/d",
                 "box long side = drag direction"),
        ND.Group("L_y/d", Ly_star, ("lengths", "L_y"), ("lengths", "d"),
                 "n_y·(√3/2)a_NN/d", "box short side"),
        ND.Group("L_x/L_y", Lx_star / Ly_star, ("lengths", "L_x"), ("lengths", "L_y"),
                 "", "★ aspect ratio — set by commensurability (NOT square)"),
        ND.Group("L_x/a_NN", Lx_star / D["a_nn_star"], ("lengths", "L_x"), ("lengths", "a_NN"),
                 "n_x", "★ commensurate — must be an integer"),
        ND.Group("L_y/a_NN", Ly_star / D["a_nn_star"], ("lengths", "L_y"), ("lengths", "a_NN"),
                 "n_y·√3/2", f"★ commensurate — integer multiple of the row spacing (n_y = {D['n_y']}, even)"),
        ND.Group("r_c/a_mean", r_c_star / a_star, ("lengths", "r_c"), ("lengths", "a_mean"),
                 "", "cutoff (number of neighbour shells)"),
        ND.Group("tau_k/tau_int", lg.ratio("times", "tau_k", "tau_int"),
                 ("times", "tau_k"), ("times", "tau_int"), "",
                 "★ ratio of the two stiffnesses — the trap is 218x faster"),
        ND.Group("dt/tau_k", lg.ratio("times", "dt", "tau_k"),
                 ("times", "dt"), ("times", "tau_k"), "", "integration resolution (governing scale)"),
        ND.Group("T_obs/tau_B", lg.ratio("times", "T_obs", "tau_B"),
                 ("times", "T_obs"), ("times", "tau_B"), "", "observation window (all four phases)"),
        ND.Group("tau_cross/tau_B", lg.ratio("times", "tau_cross", "tau_B"),
                 ("times", "tau_cross"), ("times", "tau_B"), "", "box traverse"),
        ND.Group("n_cell", n_cell, ("times", "tau_cross"), ("times", "tau_cell"),
                 "L_x/a_NN", "lattice periods per traverse = n_x"),
        ND.Group("St", lg.ratio("times", "tau_p", "tau_B"),
                 ("times", "tau_p"), ("times", "tau_B"), "tau_p/tau_B", "inertia vs diffusion"),
    ]
    checks = [
        C.Check("model", "inertia negligible  tau_p/tau_k", lg.ratio("times", "tau_p", "tau_k"),
                C.GATE, "<=",
                "tau_dyn = the fastest scale of interest in this system = tau_k (trap). "
                "Independent of dt (bd-physics §4)"),
        C.Check("integration", "trap resolved       dt/tau_k", lg.ratio("times", "dt", "tau_k"),
                C.GATE, "<=",
                f"bias ≈ {C.bias_from_dt(g('times', 'dt'), g('times', 'tau_k')):.3f}% "
                "(closed form for a linear system — the trap is linear, so it holds "
                "exactly)"),
        C.Check("integration", "pair resolved       dt/tau_int", lg.ratio("times", "dt", "tau_int"),
                C.GATE, "<=",
                f"τ_int = γ/U''(r_min={D['r_min_star']:.3f}d), criterion {D['crit']}. "
                "dt was set by the trap, so there is a lot of margin"),
        C.Check("integration", "advection resolved   dt/tau_v", lg.ratio("times", "dt", "tau_v"),
                C.GATE, "<=",
                "the trap centre must not move more than 1% of a diameter in one step "
                "(a moving boundary condition)"),
        C.Check("geometry", "cutoff              r_c/(L_min/2)", r_c_star / (L_min / 2), 1.0, "<=",
                f"minimum image (bd-hoomd trap 6). Rectangular, so the **short side** "
                f"L_{'y' if Ly_star <= Lx_star else 'x'} is the basis. Violating it gave the historical +1856% case"),
        # ★ Both axes together. x must be an integer multiple of a_NN; y an **even**
        #   multiple of the row spacing (√3/2)a_NN. With an odd row count the stagger
        #   fails to join across the periodic boundary (build_ledger raises first).
        C.Check("geometry", "commensurate hex    max|integer deviation|",
                max(abs(Lx_star / D["a_nn_star"] - D["n_x"]),
                    abs(Ly_star / (math.sqrt(3) / 2 * D["a_nn_star"]) - D["n_y"]),
                    float(D["n_y"] % 2)),
                1e-9, "<=",
                f"★ Is the lattice commensurate with the periodic box? "
                f"L_x = {D['n_x']}·a_NN, "
                f"L_y = {D['n_y']}·(√3/2)a_NN, n_y even. If incommensurate, defects are "
                f"injected at the seam — and the observables (lattice deformation "
                f"field, ψ₆) are sensitive to exactly that"),
        C.Check("geometry", "core margin        r_table_min/r_min", R_TABLE_MIN / D["r_min_star"],
                1.0, "<=",
                f"pair.Table trap 11: the force is 0 for r<{R_TABLE_MIN}d"),
        C.Check("geometry", "wake healing        v tau_int/L_x", lg.ratio("times", "tau_int", "tau_cross"),
                1.0, "<=",
                "★ The box is periodic, so the probe returns into its own wake. The "
                "healing distance v·tau_int must be shorter than the box side along "
                "the drag direction"),
        C.Check("statistics", "observation window  T_obs/tau_k", lg.ratio("times", "T_obs", "tau_k"),
                100.0, ">=", "in units of the trap correlation time", hard=False),
        C.Check("statistics", "measurability       SNR", snr, 1.0, ">=",
                "★ SNR<1 — the drag signal is buried in the thermal fluctuation. It is "
                "invisible in a single sample and only recoverable by averaging (see "
                "the precision check below)", hard=False),
        C.Check("statistics", "drag-force precision SEM/dr_ss [%]", prec_pct, STAT_TARGET_PCT, "<=",
                f"★ One traverse gives {prec_pct:.2f}% (target {STAT_TARGET_PCT:g}%). "
                f"Independent samples T_obs/2tau_k = {n_indep:,.0f}. "
                f"{math.ceil((prec_pct/STAT_TARGET_PCT)**2):.0f} traverses would reach "
                f"the target — "
                f"or use the deformation field of all {sys_['N']} lattice particles as "
                f"the observable", hard=False),
        C.Check("statistics", "lattice periods     L_x/a_NN", n_cell, N_CYCLE_TARGET, ">=",
                "★ If the drag force is modulated at the lattice period (stick-slip), a "
                "period average is needed. "
                "One traverse only grows as sqrt(N)", hard=False),
    ]
    return groups, checks, Gamma, dict(prec_pct=prec_pct, n_indep=n_indep, n_cell=n_cell,
                                       snr=snr)


def report_blocks(sys_, lg, extra, n_warm, n_equil, n_prod, n_relax,
                  tau_int_steps, args):
    D = lg.derived
    inp = [R.kv(k, f"{sys_[k].value:~.4gP}", sys_[k].tier, sys_[k].source[:44])
           for k in ("d", "T", "eta", "rho_p", "k_t", "v_x")]
    inp += [
        R.kv("N", f"{sys_['N']}", 0,
             f"sketch 'N ~ 300' -> commensurate hex {D['n_x']}×{D['n_y']} (Δ{100*(sys_['N']-300)/300:+.1f}%)"),
        R.kv("A", f"{sys_['A']:g}", 3, "★proposed — the value for which 1-B confirmed a hexagonal crystal"),
        R.kv("phi", f"{sys_['phi']:g}", 3, "★proposed — same as 1-B (not written on the sketch)"),
    ]
    der = [f"  {k:<8} = {D[k].to_compact():~.4gP}" for k in ("gamma", "D_t", "m")]
    der += [
        f"  state estimate = {D['state']}  (Lindemann sigma_bond/a_NN = "
        f"{D['u_rms_rel']:.4f}, criterion 0.15)",
        f"  Δr_ss = γv/k_t = {D['dr_ss'].to('nm'):~.3fP}   vs   ℓ_k = {D['l_k'].to('nm'):~.2fP}"
        f"   → SNR = {extra['snr']:.4f}",
        f"  ★ Δr_ss is the probe's **bare Stokes lag**. What the lattice adds sits on",
        f"    top of it, and its size has no prediction (measurement) — which is the",
        f"    reason to compute this system at all.",
        f"  ★ The trap is overwhelmingly stiffer than the pair "
        f"(k* = {lg.ratio('energies','k_t_d2','kT'):.3g}"
        f"  vs  Γ = {float(U_star(D['a_star'], sys_['A'], sys_['wca_eps_kT'])):.2f})",
        f"    -> the probe is a 'hard' driven probe, nearly pinned to the trap "
        f"(constant-velocity boundary condition).",
    ]
    plan = [
        f"  dt      = {D['dt'].to_compact():~.4gP}"
        f"  = {lg.ratio('times', 'dt', 'tau_B'):.3e} τ_B"
        f"  = 10⁻²·{D['tau_dt_name']}",
        f"  T_obs   = {D['T_obs'].to_compact():~.4gP}"
        f"  = {lg.ratio('times', 'T_obs', 'tau_B'):.2f} τ_B  (box traverse L_x/v_x)",
        f"  phases  = warm {n_warm:,} + equil {n_equil:,} + drag {n_prod:,}"
        f" + relax {n_relax:,}  =  {n_warm+n_equil+n_prod+n_relax:,} steps  x N={sys_['N']}",
        f"            (warm/equil/relax are"
        f" {args.warm:g}/{args.equil:g}/{args.relax:g} x tau_int ="
        f" {tau_int_steps:,.0f} steps)",
        f"  lattice = commensurate hex {D['n_x']}x{D['n_y']} = {sys_['N']}"
        f"   L_x×L_y = {D['Lx_star']*float(D['d'].to('um').magnitude):.2f}"
        f" × {D['Ly_star']*float(D['d'].to('um').magnitude):.2f} µm"
        f"  (aspect ratio {D['Lx_star']/D['Ly_star']:.4f})",
        f"  ⚠ One traverse gives ⟨F_drag⟩ to {extra['prec_pct']:.2f}% and only"
        f" {extra['n_cell']:.0f} lattice periods.",
        f"    Using the lattice deformation field ({sys_['N']} particles) as the primary"
        f" observable is the design intent"
        f" (observation.yaml B2).",
        "  ⚠ Commensurability sets the aspect ratio — it is not square. Finite-size"
        " artefacts may therefore be asymmetric in x and y, and whether that shows up"
        " in the deformation field is unverified.",
    ]
    return inp, der, plan


# ════════════════════════════════════════════════════════════════════════
# ⑤ L4 — build the system from the spec alone (bdbot.run runs it and judges)
# ════════════════════════════════════════════════════════════════════════
def hex_lattice(n_x: int, n_y: int, a_nn: float) -> np.ndarray:
    """Commensurate hexagonal lattice (n_x·n_y particles).

    Rows run along x; odd rows are staggered by a_NN/2.

    The box is `[-L_x/2, L_x/2) x [-L_y/2, L_y/2)`, `L_x = n_x a_NN`,
    `L_y = n_y (√3/2) a_NN`. n_y must be even for the stagger to join across the
    periodic boundary.
    """
    row = math.sqrt(3) / 2 * a_nn
    Lx, Ly = n_x * a_nn, n_y * row
    j, i = np.divmod(np.arange(n_x * n_y), n_x)
    x = (i + 0.5 * (j % 2)) * a_nn - Lx / 2
    y = j * row - Ly / 2
    return np.c_[x, y]


PH_WARM, PH_EQ, PH_DRAG, PH_RELAX = "warm-up", "equilibrium", "drag", "relaxation"
TAIL_FRAC = 0.2            # tail fraction taken as the relaxation 'end value'


def fit_relaxation(t, y, y_eq):
    """Relaxation fit `A·exp(-t/tau)+C` plus a recovery fraction that does **not**
    depend on the fit.

    ⚠️ Running `curve_fit` unbounded the first time **ran away to a straight line**:
       tau = 6.9e4 tau_int, C = -404 kT (physically impossible), and the error on tau
       was 400x the value itself.
       When tau >> the observation window, `A e^{-t/tau} ≈ A(1 - t/tau)`, so a straight
       line fits a noisy decay almost equally well. Hence the bounds:
         C   within the data range        (the energy cannot go negative)
         tau within 3x the window        (longer than that and this window cannot
                                          measure it in the first place)
       And if tau pins to the upper bound, that is **reported** as "the relaxation did
       not finish inside the window".

    The recovery fraction is independent of the fit — how far the mean of the last 20%
    has returned to the equilibrium line.
    Even when the fit fails, that number still answers the question.
    """
    t = np.asarray(t, float); y = np.asarray(y, float)
    out = {"n": int(t.size)}
    if t.size < 10:
        return out
    tr = t - t[0]
    win = float(tr[-1])
    n_tail = max(3, int(TAIL_FRAC * t.size))
    y0 = float(y[:n_tail].mean())            # just after relaxation starts
    y_end = float(y[-n_tail:].mean())        # end of relaxation
    denom = y0 - y_eq
    out.update(y_start=y0, y_end=y_end, window=win,
               recovered_frac=float((y0 - y_end) / denom) if denom else float("nan"),
               residual=float(y_end - y_eq))
    try:
        from scipy.optimize import curve_fit
        f = lambda x, a, tau, c: a * np.exp(-x / tau) + c
        lo = (-10 * abs(denom) - 1, win / 200.0, float(y.min()) - 0.5)
        hi = (10 * abs(denom) + 1, 3.0 * win, float(y.max()) + 0.5)
        p0 = (denom, max(win / 5, win / 100), y_eq)
        p0 = tuple(min(max(v, l), h) for v, l, h in zip(p0, lo, hi))
        pf, pc = curve_fit(f, tr, y, p0=p0, bounds=(lo, hi), maxfev=40000)
        tau = float(pf[1])
        # ★ The criterion is not "did it fit inside the window" but **how many
        #   e-foldings did we see**. A real contradiction occurred: tau = 2.5x the
        #   window while the recovery fraction came out 94%. That means a single
        #   exponential does not explain it (fast early + slow tail), so tau must not
        #   be quoted.
        efold = win / tau if tau > 0 else 0.0
        out.update(amp=float(pf[0]), tau_star=tau, U_inf=float(pf[2]),
                   tau_sem=float(np.sqrt(abs(pc[1, 1]))), tau_over_window=tau / win,
                   e_foldings=efold, resolved=bool(efold >= 1.0))
    except Exception as e:
        out["error"] = str(e)[:120]
    return out


@RUN.builder("trap-drag-2d-hex300")
def build(spec, outdir=None) -> RUN.Build:
    """Spec -> HOOMD system. **The case YAML is never re-read** (the L2<->L4 contract is
    the spec alone).

    ★ **The four-phase protocol** — all three user questions (equilibrium vs driven
      energy, relaxation dynamics, change in g(r) and defects) are answered inside one
      trajectory:

        warm   trap fixed, samples discarded   perfect-lattice IC -> thermal equilibrium
        equil  trap fixed, samples kept        ⟨U⟩_eq · g(r)_eq · defects_eq  <- baseline
        drag   trap moves at constant v        ⟨U⟩_drag · g(r)_drag · defects_drag
                                               + the rising transient
        relax  trap **stopped** (pinned at its final position)
                                               the dynamics of the energy returning to
                                               equilibrium

      Turning the drive on and off is done with a piecewise function passed to
      `traps.make_trap(drive=...)` — `velocity=` is not used because it cannot be
      switched off and on.
    """
    import freud
    import hoomd.md as md

    P, Nm = spec.params, spec.numerics
    N, n_x, n_y = int(P["N"]), int(P["n_x"]), int(P["n_y"])
    a_nn, Lx, Ly = P["a_nn_star"], P["Lx_star"], P["Ly_star"]
    A, eps, r_c = P["A"], P["wca_eps"], P["r_c_star"]
    k_star, v_star = P["k_star"], P["v_star"]
    r_tab_min = P["r_table_min_star"]
    dt_star, seed = float(Nm["dt_star"]), int(Nm["seed"])
    n_warm, n_equil = int(Nm["n_warm"]), int(Nm["n_equil"])
    n_drag, n_relax = int(Nm["n_prod"]), int(Nm["n_relax"])

    pos0 = hex_lattice(n_x, n_y, a_nn)
    sim = SIM.make_sim(SIM.frame_2d(pos0, (Lx, Ly)), seed=seed)

    cell = md.nlist.Cell(buffer=0.4)
    # r^-3 tail — pair.Table. ★ endpoint=False (trap 10), shifted at the cutoff
    rr = np.linspace(r_tab_min, r_c, 1000, endpoint=False)
    tab = md.pair.Table(nlist=cell, default_r_cut=r_c)
    tab.params[("A", "A")] = dict(r_min=r_tab_min,
                                  U=A / rr**3 - A / r_c**3, F=3 * A / rr**4)
    wca = SIM.wca(cell, epsilon=eps, sigma=1.0)          # traps 4 and 11 (core as a separate force)

    # ★ The trap acts on one probe only (sketch: the X mark is on a single circle).
    #   All the rest get k=0.
    probe = (n_y // 2) * n_x + n_x // 2
    k_arr = np.zeros(N)
    k_arr[probe] = k_star

    # ── Drive protocol: stopped -> constant v -> stopped (piecewise linear) ─────
    T1 = (n_warm + n_equil) * dt_star            # drag starts (reduced time)
    T2 = T1 + n_drag * dt_star                   # drag ends -> pinned there
    _off = np.zeros((N, 3))                      # preallocated to avoid allocating every step

    def drive(t):
        _off[probe, 0] = v_star * (min(max(t, T1), T2) - T1)
        return _off

    trap = TR.make_trap(k_arr, pos0, (Lx, Ly), dt_star=dt_star, drive=drive)

    SIM.attach_brownian(sim, dt_star, [tab, wca, trap])
    gsd = (Path(outdir) / "traj_A.gsd") if outdir else None
    SIM.add_trajectory_writer(sim, gsd, GSD_EVERY or
                              max(1, (n_warm + n_equil + n_drag + n_relax) // 300))

    box = freud.box.Box(Lx=Lx, Ly=Ly, is2D=True)
    hexatic = freud.order.Hexatic(k=6, weighted=True)
    voro = freud.locality.Voronoi()
    r_max = min(r_c, min(Lx, Ly) / 2 - 1e-6)
    # ★ g(r) accumulates **separately** per phase. Piling it into one mixes equilibrium
    #   with driven.
    rdf = {ph: freud.density.RDF(bins=300, r_max=r_max, r_min=0.5)
           for ph in (PH_EQ, PH_DRAG, PH_RELAX)}
    # ★ freud raises AttributeError if `.rdf` is read before compute() (an exception,
    #   not an empty object). A zero-length phase (`--relax 0`) therefore kills
    #   finalize — this actually happened during an equilibration study, dying after
    #   the run had finished and wasting 4.5 minutes.
    rdf_used: set = set()
    anchor3 = np.c_[pos0, np.zeros(N)]

    def xy():
        return np.array(sim.state.get_snapshot().particles.position, dtype=float)

    def pe_pair():
        """Potential energy of the colloidal system per particle.

        **The trap energy is excluded** — that is stored in the optical tweezer, not in
        the system. The trap side is measured separately.
        """
        return float(np.array(tab.energies).sum() + np.array(wca.energies).sum()) / N

    def sample(timestep, phase):
        p = xy()
        d_probe = trap.displacement(sim.state, timestep)[probe]
        # Lattice deformation field: displacement from the initial lattice site
        # (probe excluded, minimum image).
        # ★ Subtract the rigid translation — the dragged probe pulls the whole lattice
        #   along slightly, and that is not deformation. Without subtracting it,
        #   dev_rms keeps growing with run length (drift, not deformation).
        dev = SIM.minimum_image(p - anchor3, (Lx, Ly))[:, :2]
        dev = np.delete(dev, probe, axis=0)
        dev = dev - dev.mean(axis=0)
        vn = voro.compute((box, p)).nlist
        z = np.asarray(vn.neighbor_counts)               # Voronoi coordination number
        if phase in rdf:
            rdf[phase].compute((box, p), reset=False)
            rdf_used.add(phase)
        u_pair = pe_pair()
        u_trap = float(np.array(trap.energies).sum())
        return {
            "dx_probe": float(d_probe[0]), "dy_probe": float(d_probe[1]),
            "u_pair": u_pair, "u_trap": u_trap,
            "u_total": u_pair + u_trap / N,
            "psi6": float(np.abs(hexatic.compute((box, p), neighbors=vn)
                                 .particle_order).mean()),
            # ★ Defect = Voronoi coordination != 6. Dislocations appear as 5-7 pairs.
            "n_def": int((z != 6).sum()), "n5": int((z == 5).sum()),
            "n7": int((z == 7).sum()),
            "dev_rms": float(np.sqrt((dev**2).sum(axis=1).mean())),
            "dev_max": float(np.sqrt((dev**2).sum(axis=1)).max()),
            "min_sep": float(np.asarray(vn.distances).min()),
        }

    def finalize(cols):
        ph = cols["phase"]
        sel = {k: (ph == k) for k in (PH_EQ, PH_DRAG, PH_RELAX)}
        t = cols["_t_step"] * dt_star                     # reduced time
        obs, extra = [], {}

        def stat(key, mask):
            v = cols[key][mask]
            return (float(v.mean()), float(ST.block_sem(v))) if v.size else (float("nan"),) * 2

        # ── ① Equilibrium vs driven energy ────────────────────────────────
        u_eq, u_eq_e = stat("u_pair", sel[PH_EQ])
        u_dr, u_dr_e = stat("u_pair", sel[PH_DRAG])
        # Only the **second half** of the drag phase counts as steady state (the first
        # half is the rising transient).
        dmask = sel[PH_DRAG].copy()
        if dmask.sum() >= 4:
            idx = np.flatnonzero(dmask)
            dmask[idx[: len(idx) // 2]] = False
        u_ss, u_ss_e = stat("u_pair", dmask)
        dU = u_ss - u_eq
        dU_e = math.hypot(u_ss_e, u_eq_e)
        obs += [
            MET.observable("U_pair_equilibrium", u_eq, None, "kT/particle",
                           role="measurement",
                           note=f"lattice energy while the trap is fixed "
                                f"(±{u_eq_e:.4g})"),
            MET.observable("U_pair_driven", u_ss, None, "kT/particle", role="measurement",
                           note=f"drag steady state (second half) (±{u_ss_e:.4g})"),
            MET.observable("dU_drive", dU, None, "kT/particle", role="measurement",
                           note=f"excess energy the drive stored in the lattice = "
                                f"{dU:+.5g} "
                                f"± {dU_e:.4g} kT/particle ({dU/u_eq*100:+.3f}%). "
                                f"{dU*N:+.4g} kT in total. No prediction — this IS the "
                                f"answer"),
        ]

        # ── ② Relaxation dynamics ─────────────────────────────────────────
        rel = sel[PH_RELAX]
        fit = fit_relaxation(t[rel], cols["u_pair"][rel], u_eq) if rel.sum() else {}
        # ★ Defects relax too — possibly on a **different timescale** from the energy,
        #   so they are measured separately.
        fit_def = (fit_relaxation(t[rel], cols["n_def"][rel].astype(float), 0.0)
                   if rel.sum() else {})
        extra["relax_fit"] = fit
        extra["relax_fit_defects"] = fit_def
        tau_int = spec.reduced("times", "tau_int")
        if fit:
            tau_rel = fit.get("tau_star", float("nan"))
            ok_fit = fit.get("resolved", False)
            obs.append(MET.observable(
                "U_recovered_frac", fit.get("recovered_frac", float("nan")), None, "1",
                role="measurement",
                note=f"fraction by which the mean of the last 20% of the relaxation has "
                     f"returned to the equilibrium line. "
                     f"Residual {fit.get('residual', float('nan')):+.4g} kT/particle "
                     f"(equilibrium fluctuation ±{u_eq_e:.4g}). "
                     f"★ This number is independent of the fit"))
            obs.append(MET.observable(
                "tau_relax", tau_rel, None, "tau_B", role="measurement",
                note=(f"⟨U⟩ relaxation time = {tau_rel/tau_int:.2f} tau_int "
                      f"({fit.get('tau_over_window', float('nan')):.2f}x the window). "
                      if ok_fit else
                      f"★ The relaxation did not finish inside the window, so tau "
                      f"cannot be determined "
                      f"(tau={tau_rel/tau_int:.3g} tau_int is pinned to the upper "
                      f"bound). ")
                     + f"The local cage tau_int has no reason to match the collective "
                       f"relaxation — it is only a yardstick"))

        # ── ③ g(r) and defects ────────────────────────────────────────────
        for key, unit in (("n_def", ""), ("n5", ""), ("n7", ""), ("psi6", "1")):
            e_m, _ = stat(key, sel[PH_EQ])
            d_m, _ = stat(key, dmask)
            r_m, _ = stat(key, rel)
            extra[key] = {"equilibrium": e_m, "driven": d_m, "relaxed": r_m}
            obs.append(MET.observable(
                f"{key}_driven_vs_eq", d_m - e_m, None, unit, role="measurement",
                note=f"equilibrium {e_m:.4g} -> driven {d_m:.4g} -> after relaxation "
                     f"{r_m:.4g} [{unit}]"))

        gr, arrays = {}, {}
        for name in rdf_used:                    # ★ only what was actually computed (comment above)
            gr[name] = np.asarray(rdf[name].rdf)
        if gr:
            arrays["gr_r"] = np.asarray(rdf[next(iter(rdf_used))].bin_centers)
            key = {PH_EQ: "gr_eq", PH_DRAG: "gr_drag", PH_RELAX: "gr_relax"}
            arrays.update({key[k]: v for k, v in gr.items()})
            # First peak — a one-number yardstick for whether the structure changed
            for k, v in gr.items():
                i = int(np.argmax(v))
                extra.setdefault("g_r_peak", {})[k] = {
                    "r": float(arrays["gr_r"][i]), "g": float(v[i])}

        # ── Drag force (already present from the previous session) ────────
        dx = cols["dx_probe"][dmask]
        f_drag = -k_star * float(dx.mean()) if dx.size else float("nan")
        sem = k_star * float(ST.block_sem(dx)) if dx.size else float("nan")
        obs.append(MET.observable(
            "F_drag_total", f_drag, None, "kT/d", role="measurement",
            note=f"⟨F_drag⟩ = -k*·⟨Δx⟩ (steady state). Against bare Stokes "
                 f"gamma*v = {v_star:.4g}, "
                 f"excess {100*(f_drag/v_star - 1):+.2f}% (±{100*sem/v_star:.2f}%)"))
        extra.update(f_drag_kT_per_d=f_drag, f_drag_sem=sem, f_stokes_bare=v_star,
                     excess_pct=100 * (f_drag / v_star - 1), probe_index=probe,
                     U_pair_eq=u_eq, U_pair_driven=u_ss, dU_drive=dU, dU_sem=dU_e,
                     dU_total_kT=dU * N, tau_relax_star=tau_rel,
                     min_sep_d=float(cols["min_sep"].min()),
                     dev_max_d=float(cols["dev_max"].max()))
        return {"observables": obs, "extra": extra, "arrays": arrays}

    every_dr = max(1, n_drag // 2000)
    phases = [
        RUN.Phase(PH_WARM, n_warm, collect=False,
                  note="trap fixed, samples discarded "
                       "(perfect-lattice IC -> thermal equilibrium)"),
        RUN.Phase(PH_EQ, n_equil, max(1, n_equil // 200),
                  note="trap fixed, baseline ⟨U⟩_eq, g(r), defects"),
        RUN.Phase(PH_DRAG, n_drag, every_dr,
                  note=f"trap moves at constant v*={v_star:.4g}, "
                       f"crossing {n_x} lattice periods"),
        # ★ In the relaxation phase a **changing** energy IS the physics -> turn the
        #   drift check off (avoids a false warning)
        RUN.Phase(PH_RELAX, n_relax, max(1, n_relax // 400), expect_steady=False,
                  note="trap stopped (pinned at its final position), the dynamics of "
                       "the energy returning to equilibrium"),
    ]

    return RUN.Build(
        sim=sim, forces=[tab, wca, trap], n_particles=N,
        sample=sample, pe_per_particle=pe_pair,
        sample_every=every_dr, phases=phases,
        tags=["2D", "soft_repulsion", "r^-3", "WCA_core", "hex_lattice", "moving_trap",
              "microrheology", "newtonian", "driven_relaxation"],
        physical={"N": N, "n_x": n_x, "n_y": n_y, "phi": P["phi"], "A": A,
                  "Gamma": P["Gamma"], "k_star": k_star, "v_star": v_star,
                  **{k: v for k, v in spec.raw["back_transform"].items()
                     if isinstance(v, float)}},
        finalize=finalize)


def make_plots(outdir: Path, spec) -> Path:
    """Six panels — the figures answering the three user questions.

    ① ⟨U⟩/N full time series (phases shaded)   equilibrium vs driven energy
    ② relaxation phase zoom + exponential fit   relaxation dynamics
    ③ g(r) for the three phases overlaid   structural change
    ④ defect counts (total, 5-fold, 7-fold)   dislocation creation/annihilation
    """
    import json as _json

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams["font.family"] = ["DejaVu Sans"]     # * labels in English (CLAUDE.md)
    matplotlib.rcParams["axes.unicode_minus"] = False
    # ★ Log-axis ticks (10^{-1}) are mathtext, so keep the mathtext font set explicitly.
    matplotlib.rcParams["mathtext.fontset"] = "dejavusans"
    from matplotlib.ticker import FuncFormatter, NullFormatter

    z = np.load(outdir / "observables.npz", allow_pickle=False)
    m = _json.loads((outdir / "metrics.json").read_text())
    res = m.get("result", {})
    dt_star = float(spec.numerics["dt_star"])
    t = z["_t_step"] * dt_star
    ph = z["phase"]
    tau_int = spec.reduced("times", "tau_int")
    col = {PH_EQ: "tab:green", PH_DRAG: "tab:red", PH_RELAX: "tab:blue"}

    fig, ax = plt.subplots(3, 2, figsize=(13, 13.5))

    # ① Energy, full series
    a = ax[0, 0]
    for k, c in col.items():
        mk = ph == k
        if mk.any():
            a.plot(t[mk] / tau_int, z["u_pair"][mk], ".", ms=2, color=c, label=k)
    for key, c, lab in (("U_pair_eq", "green", "⟨U⟩ equilibrium"),
                        ("U_pair_driven", "red", "⟨U⟩ driven steady state")):
        if key in res:
            a.axhline(res[key], color=c, ls="--", lw=1.2, alpha=.8, label=lab)
    a.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$\langle U\rangle/N$  [kT]",
          title=f"① Lattice energy — ΔU = "
                f"{res.get('dU_drive', float('nan')):+.4g} kT/particle"
                f"  (total {res.get('dU_total_kT', float('nan')):+.4g} kT)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ② Relaxation zoom + fit
    # ② Relaxation — **log-log**. A single exponential does not fit (0.4 e-foldings yet
    #    94% recovery), so log on both axes is the right way to see whether it is a
    #    power law. The y axis is the **excess** ΔU = U - U_eq, not U itself — values
    #    near 105 show nothing on a log axis.
    a = ax[0, 1]
    mk = ph == PH_RELAX
    u_eq = res.get("U_pair_eq", float("nan"))
    fit = {}
    if mk.any():
        tr = (t[mk] - t[mk][0]) / tau_int
        du = z["u_pair"][mk] - u_eq
        fit = fit_relaxation(t[mk], z["u_pair"][mk], u_eq)

        # Log-bin averages — raw samples let ΔU go negative, which cannot be plotted
        # on a log axis
        pos = tr > 0
        edges = np.geomspace(max(tr[pos].min(), 1e-3), tr.max(), 26)
        ctr, val = [], []
        for lo_, hi_ in zip(edges[:-1], edges[1:]):
            s = (tr >= lo_) & (tr < hi_)
            if s.sum() >= 2:
                ctr.append(np.sqrt(lo_ * hi_)); val.append(du[s].mean())
        ctr, val = np.array(ctr), np.array(val)
        ok = val > 0
        a.loglog(ctr[ok], val[ok], "o-", ms=5, lw=1.4, color="tab:blue",
                 label=r"$\Delta U = \langle U\rangle/N - U_{eq}$ (log-bin mean)")
        if (~ok).any():                      # bins that went below 0 are shown, not hidden
            a.loglog(ctr[~ok], np.full((~ok).sum(), val[ok].min() * 0.5), "x",
                     ms=6, color="gray", label="ΔU ≤ 0 (equilibrium reached / noise)")
        if "tau_star" in fit:
            tt = np.geomspace(ctr.min(), ctr.max(), 200)
            a.loglog(tt, abs(fit["amp"]) * np.exp(-tt * tau_int / fit["tau_star"]),
                     "-", lw=1.6, color="k",
                     label=(fr"exponential $\tau$={fit['tau_star']/tau_int:.0f}"
                            fr"$\tau_{{int}}$ "
                            f"({fit['e_foldings']:.2f} e-folding"
                            f"{'' if fit.get('resolved') else ' — unresolved'})"))
        # Power-law guides — straight means a power law, curved means exponential
        for p_, c_ in ((-0.5, "tab:orange"), (-1.0, "tab:red")):
            ref = val[ok][0] * (ctr[ok] / ctr[ok][0]) ** p_
            a.loglog(ctr[ok], ref, "--", lw=1, color=c_, alpha=.7, label=fr"$t^{{{p_:g}}}$ guide")
        if "n_def" in z:
            a2 = a.twinx()
            nd = np.array([z["n_def"][mk][(tr >= l) & (tr < h)].mean()
                           for l, h in zip(edges[:-1], edges[1:])
                           if ((tr >= l) & (tr < h)).sum() >= 2])
            m2 = nd > 0
            a2.loglog(ctr[m2], nd[m2], "s--", ms=4, lw=1.1, color="tab:purple", alpha=.75)
            a2.set_ylabel("defect count (z!=6)", color="tab:purple")
            a2.tick_params(axis="y", colors="tab:purple")
            a2.yaxis.set_major_formatter(FuncFormatter(lambda v, _: f"{v:g}"))
            a2.yaxis.set_minor_formatter(NullFormatter())
        # Linear inset — for looking at the values near 105 directly
        ins = a.inset_axes([0.055, 0.06, 0.40, 0.30])
        ins.plot(tr, z["u_pair"][mk], ".", ms=1.5, color="tab:blue", alpha=.45)
        w = max(5, len(tr) // 25)
        ins.plot(tr, np.convolve(z["u_pair"][mk], np.ones(w) / w, mode="same"),
                 "-", lw=1.2, color="tab:cyan")
        ins.axhline(u_eq, color="green", ls="--", lw=1.1)
        lo, hi = np.percentile(z["u_pair"][mk], [0.5, 99.5])
        ins.set_ylim(min(lo, u_eq) - 0.03, hi + 0.03)          # ★ pinned near 105
        ins.set_title("linear (kT/particle)", fontsize=7)
        ins.tick_params(labelsize=6)
    rec = f"recovered {100*fit['recovered_frac']:.0f}%" if "recovered_frac" in fit else ""
    a.set(xlabel=r"$t/\tau_{int}$ after relaxation starts",
          ylabel=r"$\Delta U$  [kT/particle]",
          title=f"② Relaxation dynamics — log-log  ({rec})")
    # ★ Plain numbers on the log ticks instead of 10^{-1}. With a Hangul font the
    #   mathtext U+2212 could not be drawn and rendered as "10□1" (rcParams does not
    #   fix it — \mathdefault uses the body font as-is). The labels are now English so
    #   that specific failure is gone, but the range is 0.02-0.4, where decimal
    #   notation simply reads better. Kept for that reason.
    plain = FuncFormatter(lambda v, _: f"{v:g}")
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_formatter(plain)
        axis.set_minor_formatter(NullFormatter())
    a.legend(fontsize=6.5, loc="upper right", framealpha=.92)
    a.grid(alpha=.3, which="both")

    # ③ g(r)
    a = ax[1, 0]
    if "gr_r" in z:
        for key, name in (("gr_eq", PH_EQ), ("gr_drag", PH_DRAG), ("gr_relax", PH_RELAX)):
            if key in z:
                a.plot(z["gr_r"], z[key], "-", lw=1.3, color=col[name], label=name)
        a.axvline(spec.params["a_nn_star"], ls=":", c="gray", label=r"$a_{NN}$")
        a.axhline(1, color="k", lw=.5, alpha=.5)
    a.set(xlabel="r / d", ylabel="g(r)",
          title="③ Radial distribution function — by phase",
          xlim=(0.8, min(6.0, float(z["gr_r"].max()) if "gr_r" in z else 6.0)))
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ④ Defects
    a = ax[1, 1]
    for key, c, lab in (("n_def", "k", "total (z!=6)"), ("n5", "tab:orange", "5-fold"),
                        ("n7", "tab:purple", "7-fold")):
        if key in z:
            a.plot(t / tau_int, z[key], "-", lw=.9, color=c, alpha=.85, label=lab)
    for k, c in col.items():
        mk = ph == k
        if mk.any():
            a.axvspan(t[mk].min() / tau_int, t[mk].max() / tau_int, color=c, alpha=.07)
    nmax = float(z["n_def"].max()) if "n_def" in z else 0.0
    a.set(xlabel=r"$t/\tau_{int}$", ylabel="particle count", ylim=(0, max(nmax * 1.25, 1.0)),
          title=f"④ Defects (Voronoi coordination != 6) — N={spec.params['N']}")
    if nmax == 0:
        # ★ Zero is a result, not a bug. Starting from a commensurate perfect lattice
        #   with strong coupling (Γ=29.7), thermal fluctuations cannot change the
        #   Voronoi topology. 1-B (RSA initial placement) retained 1.3% defects — what
        #   makes the difference is **the initial placement and commensurability**, not
        #   the coupling strength.
        a.text(0.5, 0.5, "0 defects — the commensurate perfect lattice survives\n"
                         "(1-B used RSA placement and retained 1.3%)",
               transform=a.transAxes, ha="center", va="center", fontsize=10,
               bbox=dict(boxstyle="round", fc="lightyellow", alpha=.9))
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ⑤ Force on the probe — F = -k*·Δ (in steady state, the drag the system exerts
    #    on the probe)
    #    ★ Single-sample noise is k*·ℓ_k/d = 246 kT/d, **larger than the signal (~100)**
    #      (SNR=0.0985). Raw points show nothing, so a moving average is overlaid.
    k_star, v_star = spec.params["k_star"], spec.params["v_star"]
    a = ax[2, 0]
    fx, fy = -k_star * z["dx_probe"], -k_star * z["dy_probe"]
    a.plot(t / tau_int, fx, ".", ms=1.2, color="0.75", label="$F_x$ raw samples")
    w = max(9, len(fx) // 60)
    kern = np.ones(w) / w
    a.plot(t / tau_int, np.convolve(fx, kern, mode="same"), "-", lw=1.5,
           color="tab:red", label=f"$F_x$ moving average ({w} pts)")
    a.plot(t / tau_int, np.convolve(fy, kern, mode="same"), "-", lw=1.1,
           color="tab:blue", alpha=.8, label=f"$F_y$ moving average")
    a.axhline(v_star, color="k", ls="--", lw=1.3, label=fr"bare Stokes $\gamma v$ = {v_star:.1f}")
    a.axhline(0, color="k", lw=.5, alpha=.5)
    if "f_drag_kT_per_d" in res:
        a.axhline(res["f_drag_kT_per_d"], color="darkred", ls=":", lw=1.4,
                  label=f"⟨$F_x$⟩ steady state = {res['f_drag_kT_per_d']:.0f}")
    for k, c in col.items():
        mk2 = ph == k
        if mk2.any():
            a.axvspan(t[mk2].min() / tau_int, t[mk2].max() / tau_int, color=c, alpha=.07)
    lo, hi = np.percentile(np.convolve(fx, kern, mode="same"), [1, 99])
    a.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$F = -k^*\Delta$  [kT/d]",
          ylim=(min(lo, -50) - 100, max(hi, v_star) + 150),
          title="⑤ Force on the probe (trap force = drag from the system)")
    a.legend(fontsize=7, ncol=2); a.grid(alpha=.3)

    # ⑥ Folded onto the lattice period — the place L3 warned "if stick-slip, a period
    #    average is needed". Drag phase only: average F_x against the distance the
    #    probe has travelled, modulo a_NN.
    a = ax[2, 1]
    dmk = ph == PH_DRAG
    a_nn = spec.params["a_nn_star"]
    if dmk.sum() > 40:
        t_d = t[dmk] - t[dmk][0]
        x_travel = v_star * t_d                      # distance the trap centre travelled [d]
        phase_x = (x_travel % a_nn) / a_nn           # 0..1 (one lattice period)
        fxd = fx[dmk]
        nb = 12
        edges = np.linspace(0, 1, nb + 1)
        ctr = 0.5 * (edges[:-1] + edges[1:])
        mu = np.array([fxd[(phase_x >= l) & (phase_x < h)].mean()
                       for l, h in zip(edges[:-1], edges[1:])])
        se = np.array([fxd[(phase_x >= l) & (phase_x < h)].std()
                       / max(1, np.sqrt(((phase_x >= l) & (phase_x < h)).sum()))
                       for l, h in zip(edges[:-1], edges[1:])])
        a.errorbar(np.r_[ctr, ctr + 1], np.r_[mu, mu], yerr=np.r_[se, se],
                   fmt="o-", ms=5, lw=1.4, capsize=3, color="tab:red",
                   label=f"⟨$F_x$⟩ (mean over {x_travel[-1]/a_nn:.0f} periods, "
                         f"2 periods shown)")
        a.axhline(np.nanmean(mu), color="darkred", ls=":", lw=1.3,
                  label=f"overall mean {np.nanmean(mu):.0f}")
        a.axhline(v_star, color="k", ls="--", lw=1.2, label="bare Stokes")
        a.axvline(1, color="0.7", lw=.8)
        amp = (np.nanmax(mu) - np.nanmin(mu)) / 2
        a.set_title(f"⑥ Force folded onto the lattice period — modulation amplitude "
                    f"±{amp:.0f} kT/d "
                    f"({100*amp/max(abs(np.nanmean(mu)),1e-9):.0f}% of the mean)")
    else:
        a.set_title("⑥ Lattice-period folding — too few samples")
    a.set(xlabel="phase within the lattice period  (x mod $a_{NN}$)/$a_{NN}$",
          ylabel=r"$F_x$  [kT/d]")
    a.legend(fontsize=7); a.grid(alpha=.3)

    fig.suptitle(f"{spec.label}   Γ={spec.params['Gamma']:.2f}  "
                 f"φ={spec.params['phi']}  {spec.params['n_x']}×{spec.params['n_y']}"
                 f"={spec.params['N']}   run_id={spec.run_id}", fontsize=11)
    fig.tight_layout()
    p = outdir / "observables.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    return p


# ════════════════════════════════════════════════════════════════════════
# main — L3 report/spec + the L4 run
# ════════════════════════════════════════════════════════════════════════
def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--report", action="store_true", help="the L3 report only")
    ap.add_argument("--spec", action="store_true", help="the L3 spec -> specs/<run_id>.json")
    ap.add_argument("--dt-scale", type=float, default=1.0, help="dt multiplier (for a convergence check)")
    ap.add_argument("--traverse", type=float, default=1.0, help="number of box traverses (T_obs multiplier)")
    ap.add_argument("--samples", type=int, default=2000, help="number of samples in the drag phase")
    ap.add_argument("--warm", type=float, default=10.0, help="warm-up phase [tau_int]")
    ap.add_argument("--equil", type=float, default=20.0, help="equilibrium measurement phase [tau_int]")
    ap.add_argument("--relax", type=float, default=40.0, help="relaxation phase [tau_int]")
    ap.add_argument("--gsd-every", type=int, default=None,
                    help="trajectory write period [steps]. For animation — does NOT "
                         "enter the spec or the hash")
    ap.add_argument("--v", type=float, default=None,
                    help="drag velocity [µm/s]. Must be one of the L2 drag_velocity "
                         "list")
    ap.add_argument("--seed", type=int, default=20260804,
                    help="seed. ★ trap 12: HOOMD truncates to 16 bits, so an ensemble "
                         "must use small consecutive integers (1,2,3...) — seeds "
                         "differing by 65536 give the same trajectory")
    ap.add_argument("--force", action="store_true", help="re-run a completed run")
    ap.add_argument("--smoke", action="store_true", help="short run (to check the L4 wiring)")
    args = ap.parse_args()

    global GSD_EVERY
    GSD_EVERY = args.gsd_every
    sys_ = load_system(ROOT / "intake/trap-drag-2d-hex300/system.yaml", args.v)
    if args.smoke:
        # ★ Wiring check only. NOT a physics result — it runs 1/200 of a traverse.
        args.traverse, args.samples = min(args.traverse, 0.01), min(args.samples, 40)
        args.warm, args.equil, args.relax = 0.5, 1.0, 2.0
    lg = build_ledger(sys_, dt_scale=args.dt_scale, n_traverse=args.traverse,
                      warm=args.warm, equil=args.equil, relax=args.relax)
    D = lg.derived
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")

    # ── The four-phase protocol (user questions: equilibrium vs driven energy,
    #    relaxation, g(r) and defects) ──
    #   All time units are **tau_int** (the lattice's local relaxation time) — only the
    #   drag phase is set by L_x/v_x.
    #   ★ Taking the relaxation phase as 40x tau_int is a ★proposal. If the collective
    #     mode is slower than tau_int the tail is truncated, so after the run the fitted
    #     tau MUST be compared against the phase length.
    tau_int_steps = D["tau_int_steps"]
    n_warm, n_equil, n_relax = D["n_warm"], D["n_equil"], D["n_relax"]
    n_prod = D["n_drag"]
    sample_every = max(1, n_prod // args.samples) if n_prod else 1
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, Gamma, extra = analyze_scales(sys_, lg)
    tag = None
    if args.dt_scale != 1.0:
        tag = f"dt{args.dt_scale:g}"
    if args.traverse != 1.0:
        tag = (tag + "-" if tag else "") + f"tr{args.traverse:g}"
    if (args.warm, args.equil, args.relax) != (10.0, 20.0, 40.0):
        tag = (tag + "-" if tag else "") + f"w{args.warm:g}e{args.equil:g}r{args.relax:g}"
    if args.v is not None and abs(args.v - 0.5) > 1e-12:
        tag = (tag + "-" if tag else "") + f"v{args.v:g}"
    if args.seed != 20260804:
        tag = (tag + "-" if tag else "") + f"s{args.seed}"
    if args.smoke:
        tag = (tag + "-" if tag else "") + "smoke"

    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"A": sys_["A"], "phi": sys_["phi"], "N": sys_["N"],
                # ★ Rectangular box — L4 must not assume a square.
                "Lx_star": D["Lx_star"], "Ly_star": D["Ly_star"],
                "n_x": D["n_x"], "n_y": D["n_y"], "a_nn_star": D["a_nn_star"],
                "lattice": "hex_commensurate", "drag_axis": "x",
                "r_c_star": D["r_c_star"],
                "wca_eps": sys_["wca_eps_kT"], "Gamma": Gamma,
                # Trap and drag in reduced units. L4 runs from these alone.
                "k_star": lg.ratio("energies", "k_t_d2", "kT"),
                "v_star": lg.ratio("times", "tau_B", "tau_v"),
                "n_trapped": sys_["n_trapped"],
                "r_table_min_star": R_TABLE_MIN},
        numerics={"dt_star": lg.ratio("times", "dt", "tau_B"),
                  "dt_over_tau_k": args.dt_scale * C.GATE,
                  # ★ The phase lengths **determine the physics**, so they enter the
                  #   hash — a different protocol is a different run (even for the
                  #   same system).
                  "n_warm": n_warm, "n_equil": n_equil,
                  "n_prod": n_prod, "n_relax": n_relax,
                  "n_samples": args.samples,
                  "sample_every": sample_every, "seed": args.seed},
        tag=tag, nhex=12)
    run_id = spec.run_id()

    l3 = spec.validate()
    if l3:
        print("L3 INTEGRITY CHECK")
        for i in l3:
            print(str(i))
        print()

    inp, der, plan = report_blocks(sys_, lg, extra, n_warm, n_equil, n_prod, n_relax,
                               tau_int_steps, args)
    report, verdict = R.render(
        title=f"DimensionlessReport — {sys_['label']}   run_id={run_id}",
        ref=lg.ref, ledger=lg, groups=ND.groups_dict(groups), checks=checks,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(report)

    if spec.errors:
        print(f"\n❌ {len(spec.errors)} L3 integrity error(s) — the "
              f"non-dimensionalization does not hold.")
        return 1
    if verdict == "FAIL":
        print("\n❌ A hard separation check failed — no spec is written.")
        return 1
    p = spec.write(ROOT / "specs" / f"{run_id}.json")
    if args.spec or args.report:
        if args.spec:
            print(f"\nL3 spec: {p.relative_to(ROOT)}")
        return 0

    # ── L4 — run by **re-reading** the spec from disk (not the in-memory object) ────
    #   Re-reading is the point: the L2<->L4 contract only holds if the run works from
    #   the on-disk spec alone, and that is also when the hash check fires (`execute`
    #   looks at verify_hash first).
    outdir = ROOT / "runs" / run_id
    loaded = ND.load(p)
    outdir.mkdir(parents=True, exist_ok=True)
    (outdir / "report.txt").write_text(report)
    v = RUN.execute(loaded, RUN.get_builder(loaded.case), outdir,
                    force=args.force, progress=True)
    print(RUN.render_verdict(v))
    if v["status"] == RUN.OK:
        print(f"plot: {make_plots(outdir, loaded).relative_to(ROOT)}")

    return 0 if v["status"] in (RUN.OK, "skipped") else 1


if __name__ == "__main__":
    sys.exit(main())
