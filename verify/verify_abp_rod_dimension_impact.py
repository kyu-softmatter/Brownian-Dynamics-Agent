"""What `abp-rod-2d-run-flip` would cost if it were 3D — the impact analysis
behind its `dim.basis`.

`D9` is closed by `knowledge/wiki/concepts/dimensionality-has-no-default.md`,
which requires every case to declare *why* its dimension is what it is and
*what would change* if it flipped. This script computes the second half for
`abp-rod`, the only case where the dimension is provably load-bearing.

Everything is read from `intake/abp-rod-2d-run-flip/system.yaml`; nothing is
hard-coded. **The self-check is section ①**: the 2D formula is re-derived and
compared against the `gamma_bar_2d` the case already records. If that does not
reproduce, the 3D number below is meaningless and the script fails.

Two things are cited rather than derived here:
  · `<n(0)·n(t)> = exp(-(d-1) D_r t)`   skill `bd-hoomd` trap 14 (measured)
  · `D_eff = D_t + v0^2/(d*Lambda)`     skill `bd-hoomd` trap 14 (measured to 1.5%)

⚠️ The `tau_eff` composition `1/(1/tau_r + 1/tau_tumble)` is the case's own, and
`intake/.../observation.yaml` model_notes marks it **not_verified** (rule 7). So
section ③ reports a *range* over the one assumption that is not settled, rather
than a single 3D number.

    $PY verify/verify_abp_rod_dimension_impact.py
"""
import pathlib
import sys

import yaml

ROOT = pathlib.Path(__file__).resolve().parent.parent
CASE = ROOT / "intake/abp-rod-2d-run-flip/system.yaml"
K_B = 1.380649e-23        # J/K -- bdbot.constants is the single source of truth;
                          # duplicated here only so this script stays standalone

results = []


def report(ok, label, detail=""):
    results.append(bool(ok))
    print(f"  {'OK ' if ok else 'XX '} {label:<48}{detail}")


sys_ = yaml.safe_load(CASE.read_text())
fr, ds = sys_["friction"], sys_["derived_scales"]
z_par = fr["zeta_parallel"]["value"]          # long axis
z_perp = fr["zeta_perp"]["value"]             # short axes: 1 of them in 2D, 2 in 3D
kT = K_B * sys_["medium"]["temperature"]["value"]
v0 = sys_["active"]["speed"]["value"]                  # um/s
tau_tumble = sys_["active"]["tumble_interval"]["value"]  # s
tau_r_2d = ds["tau_r"]["value"]               # s, thermal director decay in 2D

print("=" * 78)
print("(1) self-check -- re-derive the 2D friction the case already records")
print("=" * 78)
#  A prolate ellipsoid has one long axis and two short ones. In the plane the
#  particle samples 1 long + 1 short; in 3D, 1 long + 2 short.
g_2d = 2.0 / (1 / z_par + 2 / z_perp - 1 / z_perp)     # == 2/(1/zp + 1/zt)
g_3d = 3.0 / (1 / z_par + 2 / z_perp)
recorded = fr["gamma_bar_2d"]["value"]
rel = abs(g_2d - recorded) / recorded
#  Tolerance is set by how the value is WRITTEN, not by float precision:
#  `gamma_bar_2d: 7.2122e-9` carries 5 significant figures, so agreement can
#  never be better than ~1e-5. A tighter bar fails on the rounding and says
#  nothing about the physics -- it did, at 1e-9, on this script's first run.
#  Still 100x stricter than `physical.verify`'s rtol=1e-3 for the same job.
report(rel < 1e-5, "gamma_bar_2d reproduced from zeta_par, zeta_perp   ",
       f"{g_2d:.6e} vs {recorded:.6e}  ({100*rel:.5f}%)")
if not all(results):
    print("\nthe 2D formula does not reproduce -- the 3D numbers below mean nothing")
    sys.exit(1)

D_t_2d, D_t_3d = kT / g_2d * 1e12, kT / g_3d * 1e12    # um^2/s

print()
print("=" * 78)
print("(2) translational -- small")
print("=" * 78)
print(f"  gamma_bar   2D {g_2d:.4e}   3D {g_3d:.4e} kg/s      "
      f"{g_3d/g_2d:+.1%} " .replace("+1", "+"))
print(f"  D_t         2D {D_t_2d:.4f}       3D {D_t_3d:.4f} um^2/s   "
      f"{(D_t_3d/D_t_2d - 1):+.1%}")
report(abs(g_3d / g_2d - 1) < 0.10,
       "friction moves by less than 10%", f"{(g_3d/g_2d - 1):+.2%}")

print()
print("=" * 78)
print("(3) rotational -- this is what actually moves")
print("=" * 78)
print("  trap 14: <n(0).n(t)> = exp(-(d-1) D_r t), so the director decays")
print("           TWICE as fast in 3D at the same physical D_r.")
print("  trap 14: D_eff = D_t + v0^2/(d*Lambda),  Lambda = 1/tau_eff")
print()


def d_eff(D_t, tau_r, dim):
    tau_eff = 1.0 / (1 / tau_r + 1 / tau_tumble)
    return D_t + v0 * v0 * tau_eff / dim, tau_eff


scenarios = [
    ("2D  (as built)", 2, D_t_2d, tau_r_2d),
    ("3D  Lambda_r = 2 D_r", 3, D_t_3d, tau_r_2d / 2),   # trap 14
    ("3D  Lambda_r unchanged", 3, D_t_3d, tau_r_2d),     # sensitivity bound
]
print(f"  {'scenario':26}{'tau_eff':>9}{'D_eff':>10}{'2*d*D_eff':>11}"
      f"{'D_eff rel':>11}")
base = None
out = {}
for name, dim, D_t, tau_r in scenarios:
    De, te = d_eff(D_t, tau_r, dim)
    pre = 2 * dim * De
    base = base or De
    out[name] = (De, pre)
    print(f"  {name:26}{te:>9.4f}{De:>10.3f}{pre:>11.2f}{De/base:>10.2f}x")

lo = min(out[n][0] for n in out if n.startswith("3D")) / base
hi = max(out[n][0] for n in out if n.startswith("3D")) / base
report(hi < 0.75, "D_eff moves by more than 25% under EVERY assumption",
       f"3D/2D in [{lo:.2f}, {hi:.2f}]")
report(abs(g_3d / g_2d - 1) < abs(1 - hi),
       "and rotation dominates friction", "the (d-1) factor, not Perrin")

print()
print("=" * 78)
print("(4) what a 3D switch would additionally require")
print("=" * 78)
for line in [
    "friction     gamma_bar_3d = 3/(1/z_par + 2/z_perp) replaces gamma_bar_2d",
    "constraint   H1 'MSAD exact -- 2D has one rotation axis' is rewritten",
    "HOOMD        rotational_diffusion must be passed as 2*D_r_phys (trap 14)",
    "archive      the 4 existing runs no longer describe the system",
    "first-ever   trap 14 says 'every case here is 2D, so this has no effect",
    "             today. It bites the moment anyone does 3D active matter.'",
]:
    print(f"  {line}")

print()
print("=" * 78)
n_ok = sum(results)
print(f"{n_ok}/{len(results)} passed")
print("=" * 78)
sys.exit(0 if n_ok == len(results) else 1)
