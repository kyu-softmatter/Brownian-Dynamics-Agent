"""Draft scale table for trap-2d-5um -- computed with pint-attached units."""
import pint

u = pint.UnitRegistry()
Q = u.Quantity

kB = Q(1.380649e-23, "J/K")

# ── values read off the sketch (dimensional) ──────────────────────────
k_t = Q(10, "pN/um").to("N/m")
T = Q(300, "K")
kT = (kB * T).to("J")

# ── values absent from the sketch that have to be filled in ───────────
eta = Q(0.851, "mPa*s").to("Pa*s")     # water @300K (handbook). Medium not stated -- assumed
rho_p = Q(2000, "kg/m^3")              # silica assumed. Used only to compute tau_p


def table(d, label):
    d = d.to("m")
    gamma = (3 * 3.141592653589793 * eta * d).to("kg/s")
    D_t = (kT / gamma).to("um^2/s")
    tau_B = (d**2 / D_t).to("s")
    tau_k = (gamma / k_t).to("s")
    l_k = ((kT / k_t) ** 0.5).to("nm")
    x2 = (kT / k_t).to("um^2")
    m = (rho_p * (3.141592653589793 / 6) * d**3).to("kg")
    tau_p = (m / gamma).to("s")
    k_star = float((k_t * d**2 / kT).to("dimensionless"))

    print(f"\n{'='*74}\n{label}   d = {d.to('um'):~.3fP}\n{'='*74}")
    print(f"  γ      = {gamma:~.4eP}")
    print(f"  D_t    = {D_t:~.4fP}")
    print(f"  m      = {m:~.3eP}   (rho_p = {rho_p:~P} assumed)")
    print()
    print(f"  timescales (smallest first)")
    print(f"    tau_p = m/gamma = {tau_p.to('us'):~.3fP}   inertial relaxation")
    print(f"    tau_k = gamma/k = {tau_k.to('ms'):~.3fP}   "
          f"★ trap relaxation -- the governing timescale")
    print(f"    tau_B = d^2/D_t = {tau_B:~.1fP}  diffusion (never realised, because "
          f"the trap is there)")
    print()
    print(f"  length scales")
    print(f"    l_k = sqrt(kT/k) = {l_k:~.2fP}   in-trap fluctuation width "
          f"(independent of d!)")
    print(f"    d              = {d.to('um'):~.1fP}")
    print(f"    ℓ_k/d          = {float((l_k/d).to('dimensionless')):.2e}")
    print()
    print(f"  dimensionless groups")
    print(f"    k* = k d^2/kT  = {k_star:.3e}   trap vs thermal fluctuation")
    print(f"    tau_p/tau_k    = {float((tau_p/tau_k).to('dimensionless')):.3e}   "
          f"inertia vs trap")
    print(f"    τ_k/τ_B        = {float((tau_k/tau_B).to('dimensionless')):.3e}")
    print()
    print(f"  analytic solution (the golden test)")
    print(f"    <x^2> per degree of freedom = kT/k = {x2:~.5eP}  -> rms {l_k:~.2fP}")
    print(f"    2D ⟨r²⟩       = 2kT/k = {(2*x2):~.5eP}")
    print(f"    relaxation time tau = gamma/k = {tau_k.to('ms'):~.3fP}")

    # the dt window
    dt_max = 0.01 * tau_k
    print()
    print(f"  choosing dt")
    print(f"    trap resolved: dt <= 0.01 tau_k = {dt_max.to('us'):~.2fP}")
    print(f"    recommended dt = tau_k/2000 = {(tau_k/2000).to('us'):~.3fP}   "
          f"(the value used in the golden check)")
    print(f"    compare with tau_p: tau_p = {tau_p.to('us'):~.3fP}")
    ratio = float((tau_p / (tau_k / 2000)).to("dimensionless"))
    print(f"    → τ_p/dt = {ratio:.2f}")
    if ratio > 0.01:
        print(f"       ⚠️ trips bd-physics' 'inertia negligible: tau_p/dt <= 1e-2' "
              f"check (discussed below)")
    return dict(gamma=gamma, tau_k=tau_k, tau_p=tau_p, l_k=l_k, k_star=k_star)


print(f"read off the sketch:  k_t = {k_t:~.3eP}   T = {T:~P}   kT = {kT:~.4eP}")
print(f"filled in:          eta = {eta:~P} (water@300K, assumed)   "
      f"rho_p = {rho_p:~P} (silica, assumed)")

A = table(Q(10, "um"), "reading A -- R=5µm is the particle RADIUS -> d = 10 µm")
B = table(Q(5, "um"), "reading B -- the particle DIAMETER is 5 µm -> d = 5 µm")

print(f"\n{'='*74}")
print("difference between the two readings")
print(f"{'='*74}")
print(f"  <x^2> (analytic) : identical -- kT/k does not depend on d")
print(f"  tau_k (relaxation) : {A['tau_k'].to('ms'):~.3fP}  vs  "
      f"{B['tau_k'].to('ms'):~.3fP}   (exactly 2x)")
print(f"  k*             : {A['k_star']:.3e}  vs  {B['k_star']:.3e}")
print("  -> the equilibrium result is the same; only the timescale differs by 2x. "
      "No effect on the golden check.")
