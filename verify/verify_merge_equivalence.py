#!/usr/bin/env python
"""Proves the 2026-08-29 bdbot/simbot de-duplication changed **no number**.

    $PY verify/verify_merge_equivalence.py

Every merged function is checked against the implementation it replaced, held
inline here as `_old_*`. Bit-identity, not `approx` -- `run_id` is the hash of the
spec content, so a last-bit change in `gamma` or `dt` renames runs, and 263 run
directories plus 279 spec files are named that way.

Part 2 does the opposite: it **deliberately breaks** each merge and confirms the
checker fires. CLAUDE.md working practice -- "silently passing" and "not checking"
are different things, and this project has shipped an unwired checker before
(`A4` was true but unchecked for a month; L4 `step_health` never ran across 81
runs and printed `82/82 HEALTHY`).

Exit code 0 = every merge is numerically inert and every guard against re-forking
works.
"""
from __future__ import annotations

import ast
import hashlib
import importlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import bdbot.checks as BCK          # noqa: E402
import bdbot.dt as BDT              # noqa: E402
import bdbot.health as BH           # noqa: E402
import simbot.estimators as SE      # noqa: E402
import simbot.io as SIO             # noqa: E402
import simbot.nondim as SN          # noqa: E402
import simbot.units as SU           # noqa: E402
from bdbot import constants as BC   # noqa: E402
from bdbot import materials as BM   # noqa: E402
from bdbot.runid import spec_hash   # noqa: E402
from bdbot.sim import minimum_image, period_array  # noqa: E402
from bdbot.units import Q           # noqa: E402

fails: list[str] = []


def check(label: str, got, want, *, exact: bool = True, rel: float = 0.0) -> None:
    if exact:
        ok = (got == want) if not isinstance(got, np.ndarray) else np.array_equal(got, want)
    else:
        ok = abs(got - want) <= rel * abs(want)
    print(f"  {'✓' if ok else '✗'} {label:<58} {got!r}")
    if not ok:
        fails.append(label)
        print(f"      expected {want!r}")


# ════════════════════════════════════════════════════════════════════════════
# PART 1 · the merges are numerically inert
# ════════════════════════════════════════════════════════════════════════════
print("=" * 84)
print("PART 1 — every merged function reproduces what it replaced, bit for bit")
print("=" * 84)

print("\n① water eta(T)/rho(T): the table moved from simbot.units to bdbot.constants")
#  ⚠⚠ THESE TWO DICTS ARE A DELIBERATE COPY AND MUST STAY ONE. ⚠⚠
#  They are a frozen snapshot of `simbot/units.py` **as it was before the merge**,
#  and the entire job of this file is to prove the merged code reproduces them. You
#  cannot verify "the move changed nothing" by importing the thing you moved.
#   ★ Saying so loudly because a peer session, while fixing a bug whose root cause
#     was copies of this table diverging, added a **fourth** copy of it to a
#     verifier -- and a test in `tests/test_cross_package_equivalence.py` copied the
#     Welty rows into its own body. A de-duplication pass that spawns copies while
#     it runs is the failure it exists to prevent. The difference here is that these
#     are labelled `_OLD_*`, are never read by anything but the comparison below,
#     and are wrong-by-construction the moment someone "updates" them.
#   Do NOT change these to track `bdbot.constants`. If a physical value legitimately
#   changes, this file should FAIL and the change should be argued, not absorbed.
_OLD_ETA = {293.15: 1.0016e-3, 298.15: 0.8900e-3, 303.15: 0.7972e-3, 308.15: 0.7191e-3}
_OLD_RHO = {293.15: 998.21, 298.15: 997.05, 303.15: 995.65, 308.15: 994.03}


def _old_interp(table, T_si):
    ts = sorted(table)
    if T_si in table:
        return table[T_si], False
    if T_si < ts[0] or T_si > ts[-1]:
        lo, hi = (ts[0], ts[1]) if T_si < ts[0] else (ts[-2], ts[-1])
        f = (T_si - lo) / (hi - lo)
        return table[lo] + f * (table[hi] - table[lo]), True
    for lo, hi in zip(ts, ts[1:]):
        if lo <= T_si <= hi:
            f = (T_si - lo) / (hi - lo)
            return table[lo] + f * (table[hi] - table[lo]), False
    raise AssertionError


for T in (280.0, 293.15, 295.0, 298.15, 300.0, 303.15, 306.0, 308.15, 330.0):
    check(f"eta({T} K) == pre-merge simbot", SU.water_viscosity_si(T), _old_interp(_OLD_ETA, T))
    check(f"rho({T} K) == pre-merge simbot", SU.water_density_si(T), _old_interp(_OLD_RHO, T))
check("k_B unchanged", SU.K_B, 1.380649e-23)
#  ⚠ Compare against the interpolation, NOT a rounded literal. A first pass wrote
#    `8.556640000000000e-04` and this line "failed" -- the actual float is
#    `8.556639999999996e-04`. Pinning a physical value to a hand-typed decimal is
#    how you manufacture a false failure; the seal itself is stated at rel 1e-4.
check("the sealed S2 value eta(300 K)", SU.water_viscosity_si(300.0)[0],
      _old_interp(_OLD_ETA, 300.0)[0])
check("...and it still matches the sealed 8.5566e-4 to the sealed precision",
      SU.water_viscosity_si(300.0)[0], 8.5566e-4, exact=False, rel=1e-4)

print("\n② dt gate equations: moved from simbot.nondim to bdbot.dt")
for d, l, D in ((0.03, 1.0, 1.0), (0.005, 3.0, 0.25), (0.1, 1e-6, 4.9e-13)):
    check(f"dt_max_thermal({d},{l},{D})", SN.dt_max_thermal(d, l, D), (d * l) ** 2 / (2.0 * D))
for d, l, g, F in ((0.005, 1.0, 1.0, 1037.7), (0.03, 1.0, 1.0, 0.0), (0.01, 2.0, 3.0, 5.0)):
    want = None if not F or F <= 0 else d * l * g / F
    check(f"dt_max_force({d},{l},{g},{F})", SN.dt_max_force(d, l, g, F), want)
for s, g, lam in ((0.2, 1.0, 4e6), (0.2, 1.0, None), (0.5, 2.0, 1e3)):
    want = None if not lam or lam <= 0 else s * 2.0 * g / lam
    check(f"dt_max_stability({s},{g},{lam})", SN.dt_max_stability(s, g, lam), want)
for d, l, v in ((0.01, 1.0, 2.5), (0.01, 1.0, 0.0)):
    check(f"dt_max_active({d},{l},{v})", SN.dt_max_active(d, l, v),
          None if not v or v <= 0 else d * l / v)

print("\n③ Euler-Maruyama bias: moved from simbot.estimators to bdbot.dt")
for x in (1e-4, 1e-3, 1e-2, 2e-2, 0.5):
    check(f"em_variance_bias({x})", SE.euler_maruyama_trap_variance_bias(x),
          1.0 / (1.0 - x / 2.0) - 1.0)
for b in (1e-4, 1e-3, 2.5e-3, 1e-2):
    check(f"dt_star_for_trap_bias({b})", SE.dt_star_for_trap_bias(b), 2.0 * b / (1.0 + b))
#  the linearized twin that bdbot.checks deliberately keeps
tau = Q(1.234567e-3, "s")
for b in (1e-4, 1e-3, 5e-3, 1e-2):
    check(f"checks.dt_from_bias(tau,{b}) == 2*b*tau",
          BCK.dt_from_bias(tau, b).to("s").magnitude, (2 * b * tau).to("s").magnitude)
for r in (1e-4, 1e-3, 1e-2, 3.7e-3):
    dt = r * tau
    rr = float((dt / tau).to("dimensionless").magnitude)
    check(f"checks.bias_from_dt(dt,tau) == 50*r  (r={r})", BCK.bias_from_dt(dt, tau), 50.0 * rr)
check("checks.GATE is bdbot.dt.GATE", BCK.GATE is BDT.GATE, True)
#  and the gap between the two forms is exactly b -- NOT b/(1+b)
for b in (1e-4, 1e-3, 1e-2):
    check(f"em_bias_form_gap({b}) == b", BDT.em_bias_form_gap(b), b, exact=False, rel=1e-12)

print("\n④ minimum image: three inline copies collapsed into bdbot.sim")


def _old_mi(delta, L, dims=2):
    d = np.asarray(delta, float).copy()
    Lx, Ly = (L, L) if np.isscalar(L) else (L[0], L[1])
    period = np.array([float(Lx), float(Ly), float(Lx) if dims == 3 else 0.0])
    m = period > 0
    d[:, m] -= period[m] * np.round(d[:, m] / period[m])
    return d


def _old_step_disp(p0, p1, L, dim):
    d = np.asarray(p1, float) - np.asarray(p0, float)
    d[:, :2] -= L * np.round(d[:, :2] / L)
    if dim == 3:
        d[:, 2] -= L * np.round(d[:, 2] / L)
    return float(np.sqrt((d[:, :dim] ** 2).sum(axis=1).mean()))


rng = np.random.default_rng(20260829)
for L in (10.0, 32.0, (20.0, 17.320508)):
    for dims in (2, 3):
        d = rng.uniform(-40, 40, (300, 3))
        check(f"minimum_image(L={L}, dims={dims})", minimum_image(d, L, dims), _old_mi(d, L, dims))
for L in (10.0, 32.0):
    for dim in (2, 3):
        p0 = rng.uniform(-L / 2, L / 2, (200, 3))
        p1 = p0 + rng.uniform(-L, L, (200, 3))
        check(f"measure_step_displacement(L={L}, dim={dim})",
              BH.measure_step_displacement(p0.copy(), p1.copy(), L, dim),
              _old_step_disp(p0.copy(), p1.copy(), L, dim))
#  ★ the one place the three copies did NOT agree, now surfaced
check("period_array(L,dims=2)[2] == 0 (2D: z never wrapped)", period_array(32.0, 2)[2], 0.0)
check("period_array(L,dims=3)[2] == L (3D: z IS wrapped -- traps.py used to say 0)",
      period_array(32.0, 3)[2], 32.0)

print("\n⑤ guard primitives: moved from simbot.guards to bdbot.health")
import simbot.guards as SG                                          # noqa: E402
for name in ("configurational_temperature", "check_finite", "check_inside_box",
             "check_bond_lengths", "check_step_displacements",
             "assert_statistic_fluctuates", "DisplacementReport"):
    check(f"simbot.guards.{name} is bdbot.health.{name}",
          getattr(SG, name) is getattr(BH, name), True)
f = rng.normal(0, 1, (500, 3))
check("configurational_temperature reproduces the direct formula",
      SG.configurational_temperature(f, 3.0),
      float(np.mean(np.sum(f ** 2, axis=1)) / 3.0))

print("\n⑥ drag / diffusion: NOT merged (pint vs float), so pinned instead")
for T, d_um in ((293.15, 1.0), (298.15, 0.5), (300.0, 5.0)):
    eta, _ = SU.water_viscosity_si(T)
    d_si = d_um * 1e-6
    g_bd = BM.sphere_drag(Q(eta, "Pa*s"), Q(d_si, "m")).to("kg/s").magnitude
    check(f"3*pi*eta*d == 6*pi*eta*(d/2)  (T={T}, d={d_um} um)",
          g_bd, SU.stokes_drag_si(eta, d_si / 2.0), exact=False, rel=1e-15)

print("\n⑦ spec hash: NOT merged. Pin where the two agree and where they do not")
for payload, same in (({"d": 5e-6, "T": 300}, True),
                      ({"N": 512, "k_star": 1.0}, True),
                      ({"src": "물@300K 핸드북"}, False),
                      ({"점도": 8.51e-4}, False)):
    check(f"spec_hash == sha256_payload for {list(payload)[:1]}... ",
          spec_hash(payload, 12) == SIO.sha256_payload(payload)[:12], same)

print("\n⑧ the provenance gap is reported, not resolved")
gap = BC.water_viscosity_provenance_gap(300.0)
check("rel_gap(0.851 vs table) at 300 K", round(gap["rel_gap"], 5), 0.00548)
check("1 K of T costs more than 4x that gap",
      BC.water_viscosity_sensitivity_per_K(300.0) / gap["rel_gap"] > 4.0, True)


# ════════════════════════════════════════════════════════════════════════════
# PART 2 · break each merge, confirm the checker fires
# ════════════════════════════════════════════════════════════════════════════
print()
print("=" * 84)
print("PART 2 — deliberately re-fork each merge; the checker must notice")
print("=" * 84)
import tests.test_cross_package_equivalence as T2                   # noqa: E402

asleep: list[str] = []


def expect_fail(label, fn, *a) -> None:
    try:
        fn(*a)
    except AssertionError:
        print(f"  ✓ {label:<58} fires")
        return
    print(f"  ✗ {label:<58} DID NOT FIRE")
    asleep.append(label)


_o = SU.water_viscosity_si
SU.water_viscosity_si = lambda T_si: (0.851e-3, False)      # recreate the 0.545 % split
expect_fail("water table re-forked", T2.test_water_viscosity_agrees_across_packages, 300.0)
expect_fail("sealed 300 K value moved", T2.test_the_sealed_300K_value_is_unchanged)
SU.water_viscosity_si = _o

_o = SN.dt_max_thermal
SN.dt_max_thermal = lambda d, l, D: (d * l) ** 2 / (2.0 * D)        # a copy, not the object
expect_fail("dt gate re-forked", T2.test_dt_gates_are_the_same_function, "dt_max_thermal")
SN.dt_max_thermal = _o

_o = SE.euler_maruyama_trap_variance_bias
SE.euler_maruyama_trap_variance_bias = lambda d: d / 2.0
importlib.reload(T2)
expect_fail("EM bias re-forked", T2.test_em_bias_is_the_same_function)
SE.euler_maruyama_trap_variance_bias = _o
importlib.reload(T2)

_o = SIO.sha256_payload
SIO.sha256_payload = lambda ob: hashlib.sha256(
    json.dumps(ob, sort_keys=True, default=str).encode()).hexdigest()
importlib.reload(T2)
expect_fail("spec-hash divergence silently 'fixed'",
            T2.test_the_two_spec_hashes_agree_on_ascii_and_differ_otherwise,
            {"점도": 8.51e-4}, False)
SIO.sha256_payload = _o
importlib.reload(T2)

_o = SU.stokes_drag_si
SU.stokes_drag_si = lambda eta, a: 3.0 * math.pi * eta * a          # diameter formula, radius arg
importlib.reload(T2)
expect_fail("radius/diameter convention broken",
            T2.test_drag_and_diffusion_agree_across_conventions, 300.0, 5.0)
SU.stokes_drag_si = _o

# ════════════════════════════════════════════════════════════════════════════
# PART 3 · did any mechanical edit silently drop a name?
# ════════════════════════════════════════════════════════════════════════════
#  ★ Added after a peer session lost 41 lines to a search-delimited replacement
#    whose *end* marker overran (recovered via `ast.parse` + `git checkout`). The
#    rule it produced: an index- or search-delimited edit needs an asserted END
#    marker, not just an asserted start. This pass is the asserted end marker for
#    a whole refactor -- it re-resolves, at runtime, every module-level name that
#    existed before the edits.
#  ⚠ A first version compared AST *definitions* and flagged 5 files -- but 4 were
#    `def`s that became re-export `import`s, i.e. the merge working as designed. An
#    audit that cannot tell "moved" from "deleted" is the same defect as a
#    doc-scraper that cannot tell 0 matches from a pass. Hence: runtime resolution.
print()
print("=" * 84)
print("PART 3 — every pre-edit module-level name still resolves")
print("=" * 84)

BASELINE_REV = "34bc0fe"          # the tree as it was before this refactor began
#  Names that were only ever a module's own imports, never part of its API. They
#  legitimately disappear when a module becomes a re-export shim.
INCIDENTAL = {"dataclass", "field", "math", "np", "annotations", "hashlib", "json",
              "Path", "platform", "re", "subprocess", "_date", "importlib"}


def _public_names(src: str, label: str) -> set:
    tree = ast.parse(src, filename=label)
    names = set()
    for n in tree.body:
        if isinstance(n, (ast.FunctionDef, ast.AsyncFunctionDef, ast.ClassDef)):
            names.add(n.name)
        elif isinstance(n, ast.Assign):
            names |= {t.id for t in n.targets if isinstance(t, ast.Name)}
        elif isinstance(n, (ast.Import, ast.ImportFrom)):
            names |= {a.asname or a.name.split(".")[0] for a in n.names}
    return {x for x in names if not x.startswith("__")}


#  ⚠⚠ `git diff --name-only HEAD`, **not** `git diff --name-only`. ⚠⚠
#   The bare form lists UNSTAGED changes only. This script was written with it, and
#   the moment the 25 paths of this refactor got staged, Part 3 silently audited
#   **zero files** and the script still printed PASS -- measured 2026-08-29.
#   That is precisely the failure family this session wrote into
#   `docs/05-pitfalls.md`: *the check cannot distinguish the case it is testing from
#   the case where it did not run.* Found only because another session staged the
#   work and the empty output was noticed by eye, which is not a mechanism.
#   Hence the explicit emptiness check below -- an audit with nothing to audit is a
#   FAIL, never a pass.
edited = [f for f in subprocess.run(
    ["git", "diff", "--name-only", "HEAD"], capture_output=True, text=True,
    cwd=REPO).stdout.split() if f.endswith(".py") and not f.startswith("tests/")]

if not edited:
    print("  ✗ no .py changes found against HEAD — this audit covered NOTHING.")
    print("    Either the work is already committed (then run it against the")
    print("    parent commit instead), or the diff command is wrong.")
    fails.append("part 3 audited zero files (silently-empty check)")
for f in edited:
    old_src = subprocess.run(["git", "show", f"{BASELINE_REV}:{f}"],
                             capture_output=True, text=True, cwd=REPO)
    if old_src.returncode != 0:
        print(f"  · {f:<34} new file, no baseline")
        continue
    mod_name = (f[:-len("/__init__.py")] if f.endswith("/__init__.py")
                else f[:-3]).replace("/", ".")
    try:
        want = _public_names(old_src.stdout, f)
        mod = importlib.import_module(mod_name)
    except Exception as e:                                   # report, never crash
        print(f"  ✗ {f:<34} {type(e).__name__}: {e}")
        fails.append(f"{f} unresolvable")
        continue
    gone = sorted(n for n in want if not hasattr(mod, n))
    real = [n for n in gone if n not in INCIDENTAL]
    if real:
        print(f"  ✗ {f:<34} LOST {real}")
        fails.append(f"{f} lost {real}")
    else:
        note = f"  (shim dropped its own imports: {gone})" if gone else ""
        print(f"  ✓ {f:<34} all {len(want)} pre-edit names resolve{note}")

print()
print("=" * 84)
if fails or asleep:
    if fails:
        print(f"✗ FAIL — {len(fails)} merge(s) changed a number: {fails[:6]}")
    if asleep:
        print(f"✗ FAIL — {len(asleep)} checker(s) do not fire when broken: {asleep}")
    print("=" * 84)
    sys.exit(1)
print("✓ PASS — every merge is numerically inert, and every guard fires when broken")
print("=" * 84)
