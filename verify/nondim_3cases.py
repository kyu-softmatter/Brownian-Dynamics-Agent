"""L3 non-dimensionalization report -- the 3 newly confirmed cases (nothing is run).

User instruction (2026-08-04): "fill in the proposed values and show me the report
first".

Uses `bdbot`'s shared parts (ScaleLedger, Check, report.render) as-is -- which also
makes this a test bench for whether the abstraction pulled out in 1-C holds for
**three more cases**.
What differs per case (which scales enter, which checks are needed) is filled in
here.

    $PY scratch/nondim_3cases.py [abp|trap-drag|chain]
"""
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import Q, checks as C, physical as P, report as R, scales as SC  # noqa: E402

f = lambda q: float(q.to("dimensionless").magnitude) if hasattr(q, "to") else float(q)


def node(s, *path):
    """A (value, unit) node from system.yaml -> a pint Quantity."""
    cur = s.raw
    for k in path:
        cur = cur[k]
    u = cur.get("unit", "dimensionless")
    return Q(cur["value"], u if u else "dimensionless")


def ds(s, key):
    return node(s, "derived_scales", key)


# ══════════════════════════════════════════════════════════════════════
def abp_rod():
    s = P.load(ROOT / "intake/abp-rod-2d-run-flip")
    d = node(s, "particle", "diameter").to("m")
    kT = Q(1.380649e-23, "J/K") * node(s, "medium", "temperature")
    gbar = node(s, "friction", "gamma_bar_2d").to("kg/s")
    v = node(s, "active", "speed").to("m/s")
    tau_tumble = node(s, "active", "tumble_interval").to("s")
    tau_r, tau_B, tau_v = ds(s, "tau_r").to("s"), ds(s, "tau_B").to("s"), ds(s, "tau_v").to("s")
    tau_eff, l_p = ds(s, "tau_eff").to("s"), ds(s, "l_p").to("m")
    tau_p, D_r = ds(s, "tau_p").to("s"), ds(s, "D_r").to("1/s")
    L = node(s, "geometry", "box_length").to("m")
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    N = int(node(s, "particle", "count").magnitude)
    semi = node(s, "particle", "semi_axes").to("m")

    lg = SC.ScaleLedger()
    lg.lengths = {
        "d_eq     equivalent-volume sphere diameter (reference)": d,
        "2a_minor minor axis length": 2 * semi[1],
        "2a_major major axis length": 2 * semi[0],
        "l_p      v*tau_eff persistence length ★": l_p,
        "L        box": L,
    }
    lg.times = {
        "tau_p    m/gammabar inertial relaxation": tau_p,
        "dt       integration step": dt,
        "tau_v    d_eq/v advection ★": tau_v,
        "tau_eff  orientational correlation (tumble + rotational diffusion)": tau_eff,
        "tau_tumble tumble interval": tau_tumble,
        "tau_r    1/D_r thermal rotational diffusion": tau_r,
        "tau_B    d_eq^2/Dbar diffusion (reference)": tau_B,
        "T_obs    observation window": T_obs,
    }
    lg.energies = {
        "kT       thermal energy (reference)": kT.to("J"),
        "f_a*d_eq self-propulsion work": (gbar * v * d).to("J"),
    }
    lg.derived = {"gamma": gbar, "d": d, "kT": kT.to("J")}
    lg.ref = SC.thermal_reference(
        d, kT.to("J"), tau_B,
        SC.THERMAL_RATIONALE + " The reference length is the equivalent-volume sphere "
        "diameter d_eq -- an ellipsoid has separate major and minor axes, but the "
        "non-dimensionalization needs exactly one reference. The translational "
        "friction is the isotropic mean gammabar, forced by the BD constraint.")
    lg.rationale = lg.ref["rationale"]

    groups = {
        "Pe     = v d_eq/Dbar   advection vs diffusion": f(v * d * gbar / kT.to("J")),
        "D_r*   = D_r tau_B     rotation vs translation": f(D_r * tau_B),
        "l_p/d_eq               persistence length": f(l_p / d),
        "p      = a_major/a_minor aspect ratio": f(semi[0] / semi[1]),
        "tau_tumble/tau_r       ★ tumbling vs rotational diffusion":
            f(tau_tumble / tau_r),
        "zeta_perp/zeta_par     the real friction anisotropy (BD cannot produce it)":
            1.287,
        "L/d_eq                 box": f(L / d),
        "St     = tau_p/tau_B   inertia vs diffusion": f(tau_p / tau_B),
    }
    ck = [
        C.Check("model", "inertia negligible  tau_p/tau_v", f(tau_p / tau_v),
                C.GATE, "<=",
                "tau_dyn = the fastest scale of interest = tau_v (advection). "
                "BD validity; independent of dt"),
        C.Check("integration", "advection resolved  dt/tau_v", f(dt / tau_v),
                C.GATE, "<=",
                f"{100*f(dt/tau_v):.1f}% of d_eq travelled per step"),
        C.Check("integration", "rotation resolved   dt*D_r", f(dt * D_r),
                C.GATE, "<=",
                "resolves the orientational dynamics"),
        C.Check("integration", "tumble resolved     dt/tau_tumble",
                f(dt / tau_tumble), C.GATE, "<=",
                "required for the Poisson tumble approximation to hold "
                "(the bd-hoomd run-and-flip snippet)"),
        C.Check("geometry", "finite size         l_p/(L/4)", f(l_p / (L / 4)),
                1.0, "<=",
                "the persistence length must stay within 1/4 of the box or active "
                "artefacts appear", hard=False),
        C.Check("statistics", "observation window  T_obs/tau_eff",
                f(T_obs / tau_eff), 100.0, ">=",
                f"in units of the orientational correlation time. {N} independent "
                f"particles supply the statistical multiplier", hard=False),
    ]
    inp = [
        R.kv("semi_axes", "(1.0, 0.25, 0.25) µm", 1, "user-confirmed 2026-08-04",
             val_w=22),
        R.kv("medium", "water @300K", 1,
             "user-confirmed (medium) + inherited from 1-A (temperature)", val_w=22),
        R.kv("v", f"{v.to('um/s'):~.3gP}", 0, "sketch 'v <= 5 µm/s', upper bound",
             val_w=22),
        R.kv("tau_tumble", f"{tau_tumble:~.3gP}", 1,
             "user-confirmed 'tumbles every 0.5 s'", val_w=22),
        R.kv("N", str(N), 3, "★proposed independent ensemble", val_w=22),
        R.kv("tumble_angle", "uniform random (2D)", 3,
             "★proposed -- the sketch says 'flip' (180 deg)", val_w=22),
    ]
    der = [
        f"  gammabar(2D) = {gbar:~.4eP}   (Perrin, sphere limit verified)",
        f"  D̄       = {(kT.to('J')/gbar).to('um^2/s'):~.4fP}",
        f"  γ_r,z    = {node(s,'friction','gamma_rot_z'):~.4eP}   D_r = {D_r:~.4fP}",
        f"  ★ the real friction is anisotropic (zeta_perp/zeta_par=1.287) but BD can "
        f"only do isotropic -> the short-time MSD anisotropy is lost",
    ]
    plan = [
        f"  dt      = {dt.to('ms'):~.4gP}  = {f(dt/tau_B):.3e} τ_B",
        f"  T_obs   = {T_obs:~.4gP} = {f(T_obs/tau_eff):.0f} τ_eff",
        f"  steps   = {f(T_obs/dt):,.0f}   × N={N}",
    ]
    return "abp-rod (run-and-tumble ellipsoid)", s, lg, groups, ck, inp, der, plan


# ══════════════════════════════════════════════════════════════════════
def trap_drag():
    s = P.load(ROOT / "intake/trap-drag-2d-hex300")
    d = node(s, "particle", "diameter").to("m")
    kT = ds(s, "kT").to("J")
    gam = ds(s, "gamma").to("kg/s")
    k_t = node(s, "external", "stiffness").to("N/m")
    # ★ Also stale from the same rework: drag_velocity became a 7-point sweep
    #   list, so this returns an array and every downstream f() raised
    #   "only 0-dimensional arrays can be converted to Python scalars".
    #   Same selection rule as cases/trap_drag_2d.py: prefer 0.5 µm/s (the
    #   sketch value), else the first entry. Never invent a velocity.
    _v = node(s, "external", "drag_velocity").to("m/s")
    _vm = getattr(_v.magnitude, "tolist", lambda: _v.magnitude)()
    if isinstance(_vm, list):
        _pick = 0.5e-6 if any(abs(x - 0.5e-6) < 1e-15 for x in _vm) else _vm[0]
        v = Q(_pick, "m/s")
    else:
        v = _v
    tau_k, tau_int = ds(s, "tau_k").to("s"), ds(s, "tau_int").to("s")
    tau_v, tau_B, tau_p = ds(s, "tau_v").to("s"), ds(s, "tau_B").to("s"), ds(s, "tau_p").to("s")
    l_k, dr_ss = ds(s, "l_k").to("m"), ds(s, "dr_ss").to("m")
    a_mean, a_nn = node(s, "geometry", "a_mean").to("m"), node(s, "geometry", "a_nn").to("m")
    # ★ This box is RECTANGULAR since the commensurate-hexagon rework, so the
    #   single `geometry.box_length` this draft read no longer exists -- it became
    #   box_length_x / box_length_y and this script was left behind, crashing with
    #   KeyError. Which side goes where is NOT decided here:
    #   cases/trap_drag_2d.py establishes it and this follows.
    #     minimum image -> the SHORT side (it is the binding denominator)
    #     box traverse  -> the DRAG DIRECTION, L_x
    #   (abp_rod() above still reads a plain `box_length`; that box IS square.)
    L_x = node(s, "geometry", "box_length_x").to("m")   # drag direction
    L_y = node(s, "geometry", "box_length_y").to("m")
    L_min = min(L_x, L_y)                                # sets the minimum image
    r_c = node(s, "interactions", 0, "cutoff").to("m")
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    N = int(node(s, "particle", "count").magnitude)
    A = f(node(s, "interactions", 0, "amplitude_A"))
    phi = f(node(s, "geometry", "area_fraction"))
    Gamma = A / f(a_mean / d) ** 3

    lg = SC.ScaleLedger()
    lg.lengths = {
        "dr_ss    gamma*v/k drag lag ★": dr_ss,
        "l_k      sqrt(kT/k) trap fluctuation": l_k,
        "d        particle diameter (reference)": d,
        "a_mean   mean spacing": a_mean,
        "a_NN     hexagonal nearest neighbour": a_nn,
        "r_c      cutoff": r_c,
        "L_x      box, drag direction": L_x,
        "L_y      box, short side (minimum image)": L_y,
    }
    lg.times = {
        "tau_p    m/gamma inertia": tau_p,
        "dt       integration step": dt,
        "tau_k    gamma/k_t trap relaxation ★": tau_k,
        "tau_int  gamma/U''(r_min) pair": tau_int,
        "tau_v    d/v_x dragging": tau_v,
        "tau_B    d^2/D_t diffusion (reference)": tau_B,
        "T_obs    observation window (one box crossing)": T_obs,
    }
    lg.energies = {
        "kT       thermal energy (reference)": kT,
        "Gamma*kT pair coupling at the mean spacing": (Gamma * kT).to("J"),
        "k_t d^2  trap stiffness": (k_t * d**2).to("J"),
    }
    lg.derived = {"gamma": gam, "d": d, "kT": kT}
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " This system has **two stiffnesses** (the trap k_t and "
        "the pair U''), and the trap is 218x faster -> the trap sets dt.")
    lg.rationale = lg.ref["rationale"]

    groups = {
        "Gamma  = A(d/a_mean)^3 pair coupling ★": Gamma,
        "k*     = k_t d^2/kT    trap vs thermal fluctuation": f(k_t * d**2 / kT),
        "phi                    packing fraction": phi,
        "Pe     = v d/D_t       advection vs diffusion": f(v * d * gam / kT),
        "SNR    = dr_ss/l_k     ★ signal vs noise": f(dr_ss / l_k),
        "r_c/a_mean             number of neighbour shells": f(r_c / a_mean),
        "tau_k/tau_int          ★ ratio of the two stiffnesses": f(tau_k / tau_int),
        "L_x/d                  box, drag direction": f(L_x / d),
        "St     = tau_p/tau_B": f(tau_p / tau_B),
    }
    ck = [
        C.Check("model", "inertia negligible  tau_p/tau_k", f(tau_p / tau_k),
                C.GATE, "<=",
                "tau_dyn = the fastest scale of interest = tau_k (the trap)"),
        C.Check("integration", "trap resolved       dt/tau_k", f(dt / tau_k),
                C.GATE, "<=",
                f"bias ~ {C.bias_from_dt(dt, tau_k):.3f}% (for a linear system)"),
        C.Check("integration", "pair resolved       dt/tau_int", f(dt / tau_int),
                C.GATE, "<=",
                "the pair is 218x slower than the trap, so there is plenty of margin"),
        C.Check("geometry", "cutoff              r_c/(L_min/2)",
                f(r_c / (L_min / 2)),
                1.0, "<=",
                "minimum image (bd-hoomd trap 6). In 1-B this was the real gate"),
        C.Check("statistics", "observation window  T_obs/tau_k", f(T_obs / tau_k),
                100.0, ">=",
                "in units of the trap correlation time", hard=False),
        C.Check("statistics", "measurability       SNR", f(dr_ss / l_k), 1.0, ">=",
                "★ the drag signal must exceed the thermal fluctuation to be visible "
                "in a single sample. Below SNR=1 it is only recoverable by averaging",
                hard=False),
    ]
    inp = [
        R.kv("d", f"{d.to('um'):~.3gP}", 1,
             "sketch 'R=5µm' + the 1-A diameter convention", val_w=22),
        R.kv("k_t", f"{k_t.to('pN/um'):~.3gP}", 0, "sketch", val_w=22),
        R.kv("v_x", f"{v.to('um/s'):~.3gP}", 0, "sketch", val_w=22),
        R.kv("N", str(N), 0, "sketch 'N ~ 300'", val_w=22),
        R.kv("A", f"{A:g}", 3,
             "★proposed -- the value for which 1-B confirmed a hexagonal crystal",
             val_w=22),
        R.kv("phi", f"{phi:g}", 3, "★proposed -- same as 1-B", val_w=22),
    ]
    der = [
        f"  γ = {gam:~.4eP}   D_t = {(kT/gam).to('um^2/s'):~.4fP}",
        f"  Gamma = {Gamma:.2f} -> in the hexagonal-crystal regime, per the 1-B "
        f"measurement (psi6=0.885, 6-fold 98.7%)",
        f"  Δr_ss = γv/k_t = {dr_ss.to('nm'):~.3fP}  vs  ℓ_k = {l_k.to('nm'):~.2fP}",
    ]
    plan = [
        f"  dt      = {dt.to('us'):~.4gP}  = {f(dt/tau_B):.3e} τ_B   (set by the trap)",
        f"  T_obs   = {T_obs:~.4gP} (one box crossing, L_x/v_x)",
        f"  steps   = {f(T_obs/dt):,.0f}   × N={N}",
    ]
    return "trap-drag (hexagonal lattice + moving trap)", s, lg, groups, ck, inp, der, plan


# ══════════════════════════════════════════════════════════════════════
def chain_bend():
    s = P.load(ROOT / "intake/chain-bend-2d-oscill")
    d = node(s, "particle", "diameter").to("m")
    kT = ds(s, "kT").to("J")
    gam = ds(s, "gamma").to("kg/s")
    k_t = node(s, "external", "stiffness").to("N/m")
    amp = node(s, "external", "amplitude").to("m")
    tau_k, tau_chain = ds(s, "tau_k").to("s"), ds(s, "tau_chain").to("s")
    tau_fast, tau_B, tau_p = (ds(s, "tau_fast").to("s"), ds(s, "tau_B").to("s"),
                              ds(s, "tau_p").to("s"))
    l_k, dmax = ds(s, "l_k").to("m"), ds(s, "delta_max").to("m")
    k_end = ds(s, "kappa_end").to("N/m")
    L = node(s, "geometry", "contour_length").to("m")
    dt = node(s, "numerics", "dt").to("s")
    T_obs = Q(s.raw["numerics"]["production_s"], "s")
    n = int(s.raw["geometry"]["n_beads"])
    kth = node(s, "interactions", 0, "angle_stiffness").to("J")
    Mc = node(s, "interactions", 0, "critical_moment").to("N*m")
    om_lo, om_hi = s.raw["external"]["omega_range"]["value"]
    tau_w_lo = Q(1.0 / om_lo, "s")

    lg = SC.ScaleLedger()
    lg.lengths = {
        "l_k      sqrt(kT/k_t) trap fluctuation": l_k,
        "a        drive amplitude ★": amp,
        "delta_max M_c linear limit ★": dmax,
        "d        bead diameter (reference)": d,
        "L        chain contour length": L,
    }
    lg.times = {
        "tau_p    m/gamma inertia": tau_p,
        "dt       integration step": dt,
        "tau_fast gamma/lambda_max fastest bending mode ★": tau_fast,
        "tau_k    gamma/k_t trap": tau_k,
        "tau_chain gamma/kappa_center collective bending ★": tau_chain,
        "tau_w    1/omega_min lowest drive": tau_w_lo,
        "tau_B    d^2/D_t diffusion (reference)": tau_B,
        "T_obs    observation window": T_obs,
    }
    lg.energies = {
        "kT       thermal energy (reference)": kT,
        "kappa_th bond-angle stiffness ★": kth,
    }
    lg.derived = {"gamma": gam, "d": d, "kT": kT}
    lg.ref = SC.thermal_reference(
        d, kT, tau_B,
        SC.THERMAL_RATIONALE + " ★ This system's timescales span four decades -- the "
        "fastest bending mode (0.28 µs) sets dt while the mode of interest is "
        "tau_chain (1.27 ms).")
    lg.rationale = lg.ref["rationale"]

    groups = {
        "kappa_th/kT            ★ the chain is thermally stiff": f(kth / kT),
        "a/l_k                  amplitude vs thermal fluctuation = SNR": f(amp / l_k),
        "a/delta_max            ★ linear-elastic margin": f(amp / dmax),
        "kappa_end/k_t          chain vs trap stiffness": f(k_end / k_t),
        "De_max = w_max*tau_k   highest Deborah": om_hi * f(tau_k / Q(1, "s")),
        "tau_fast/tau_chain     ★ scale-separation width": f(tau_fast / tau_chain),
        "n_beads": float(n),
        "St     = tau_p/tau_B": f(tau_p / tau_B),
    }
    ck = [
        C.Check("model", "inertia negligible  tau_p/tau_chain", f(tau_p / tau_chain),
                C.GATE, "<=",
                "tau_dyn = the fastest scale **of interest** = tau_chain (collective "
                "bending). That is the band in which G'(omega) is measured"),
        C.Check("model", "note: tau_p/tau_fast", f(tau_p / tau_fast), C.GATE, "<=",
                "★ the fastest bending mode is in fact **not overdamped** "
                "(ratio 0.6). BD treats it as overdamped -- it sits 4570x away from "
                "the observed band (tau_chain), so it probably does not affect "
                "G'(omega), but that has **not been checked**", hard=False),
        C.Check("integration", "fastest mode resolved dt/tau_fast", f(dt / tau_fast),
                C.GATE, "<=",
                "against the largest eigenvalue of the stiffness matrix. Miss this "
                "and it blows up"),
        C.Check("integration", "drive resolved      dt/tau_w", f(dt / tau_w_lo),
                C.GATE, "<=",
                "resolves the lowest drive period"),
        C.Check("geometry", "linear elasticity   a/delta_max", f(amp / dmax),
                1.0, "<=",
                "★ M < M_c. Beyond it the bond slips and the harmonic angle "
                "potential is invalid (the conclusion of paper [P2])"),
        C.Check("statistics", "SNR          a/ℓ_k", f(amp / l_k), 3.0, ">=",
                "the amplitude must exceed the in-trap thermal fluctuation by enough "
                "for phase extraction to work", hard=False),
        C.Check("statistics", "observation window  T_obs*w_min/2pi",
                f(T_obs / (2 * math.pi * tau_w_lo)),
                100.0, ">=", "number of cycles at the lowest frequency", hard=False),
    ]
    inp = [
        R.kv("d", f"{d.to('um'):~.3gP}", 2,
             "[P1] '2a = 1.47 µm' -- disagrees 3.4x with the sketch's R=5µm",
             val_w=22),
        R.kv("kappa_0", "64 mN/m", 2, "[P1] measured at 10 mM MgCl2", val_w=22),
        R.kv("M_c", f"{Mc.to('pN*um'):~.3gP}", 2, "[P1] the high-salt plateau",
             val_w=22),
        R.kv("k_t", f"{k_t.to('pN/um'):~.3gP}", 0, "sketch (the paper says ~40)",
             val_w=22),
        R.kv("n_beads", str(n), 3,
             "★proposed [P1] Fig.4 -- at n=11 the amplitude window closes", val_w=22),
        R.kv("amplitude", f"{amp.to('nm'):~.3gP}", 3,
             "★proposed -- the window is 20nm << a < 429nm", val_w=22),
    ]
    der = [
        f"  EI = {node(s,'interactions',0,'flexural_rigidity'):~.4eP}"
        f"   κ_θ = EI/ℓ = {kth:~.4eP} = {f(kth/kT):.3e} kT",
        f"  kappa_end = {k_end.to('pN/um'):~.3fP} (the papers' definition, end force)"
        f"   kappa_center = {(2*k_end).to('pN/um'):~.3fP} (what the driving trap feels)",
        f"  ★ mixing the two force definitions disagrees with the papers by exactly "
        f"2x (confirmed by verification)",
        f"  δ_max = M_c L²/(12EI) = {dmax.to('nm'):~.0fP}   ℓ_k = {l_k.to('nm'):~.2fP}",
    ]
    plan = [
        f"  dt      = {dt.to('ns'):~.4gP}  = {f(dt/tau_B):.3e} τ_B   (set by the fastest bending mode)",
        f"  omega sweep = {om_lo:.0f} to {om_hi:.0f} rad/s  "
        f"(De = omega*tau_k, 0.1 to 10)",
        f"  T_obs   = {T_obs:~.4gP} (100 cycles at the lowest omega)",
        f"  steps   = {f(T_obs/dt):,.0f}   x n={n}   ★ expensive -- see the VERDICT "
        f"below",
    ]
    return "chain-bend (three-point bending microrheology)", s, lg, groups, ck, inp, der, plan


# ══════════════════════════════════════════════════════════════════════
BUILDERS = {"abp": abp_rod, "trap-drag": trap_drag, "chain": chain_bend}
want = sys.argv[1:] or list(BUILDERS)
for key in want:
    title, s, lg, groups, ck, inp, der, plan = BUILDERS[key]()
    txt, verdict = R.render(
        title=f"DimensionlessReport — {s.label}   [{title}]",
        ref=lg.ref, ledger=lg, groups=groups, checks=ck,
        input_lines=inp, derived_lines=der, run_plan_lines=plan)
    print(txt)
    print()
