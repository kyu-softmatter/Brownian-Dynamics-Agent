"""Adversarial L3 checks -- `trap-drag-2d-hex300` and `chain-bend-2d-oscill`
(2026-08-04).

Where `verify_nondim_guards.py` tests `NondimSpec` **itself**, this deliberately
breaks **the ledger, the dimensionless groups and the hash of the two new cases**.
"Silently passing" and "not checking" are different things (CLAUDE.md working
practice).

What is specific to these two cases:
  . chain-bend has no periodic boundary, so it empties the `box` role via
    `declare_absent` -- **forgetting** to empty it must be caught as a missing
    required role.
  . chain-bend's `lambda_max` is the **larger** of the bending and stretching
    blocks. Adding them (mixing decoupled degrees of freedom) underestimates dt by
    18% and only raises the cost.
  . does trap-drag's physical system cover the run_id -- this is a hole that was
    actually breached in 1-B.

    $PY scratch/verify_l3_two_cases.py
"""
import sys, json, copy, argparse, tempfile
from pathlib import Path
ROOT = Path(__file__).resolve().parent.parent; sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))
import chain_bend_2d as CB, trap_drag_2d as TD
from bdbot import nondim as ND

A = argparse.Namespace(dt_scale=1.0, cycles=CB.N_CYCLES, samples=2000, spec=False, report=True)
ok = lambda b: "✓" if b else "✗ NOT caught"
res = []

# ── (1) chain-bend: without declare_absent on box, is it caught as a missing
#        required role?
s = CB.load_system(ROOT/"intake/chain-bend-2d-oscill/system.yaml")
lg, spec, g, c, *_ = CB.build_spec(s, 85.0, A)
lg.absent.pop("box")
errs = [i for i in spec.validate() if i.level == "error"]
res.append(("chain: box undeclared -> missing required role",
            any("ledger.box" in i.where for i in errs)))
lg.declare_absent("box", "restored")

# ── (2) chain-bend: is a dt_star inconsistent with the ledger caught? (otherwise
#        HOOMD runs with a different step than the ledger claims)
spec.numerics["dt_star"] *= 1.0000001
res.append(("chain: dt_star perturbed by 1e-7", any("numerics.dt_star" in i.where
            for i in spec.validate() if i.level == "error")))
spec.numerics["dt_star"] /= 1.0000001
res.append(("chain: passes once reverted", not spec.errors))

# ── (3) chain-bend: does recomputing from the ledger catch a slightly wrong group?
De = next(x for x in spec.groups if x.name == "De"); orig = De.value
De.value = orig * 1.000001
res.append(("chain: De perturbed by 1e-6", any("groups.De" in i.where
            for i in spec.validate() if i.level == "error")))
De.value = orig

# ── (4) trap-drag: does changing the physical system change the run_id?
#        (a hole actually breached in 1-B)
t = TD.load_system(ROOT/"intake/trap-drag-2d-hex300/system.yaml")
Ta = argparse.Namespace(dt_scale=1.0, traverse=1.0, samples=2000)
lg2 = TD.build_ledger(t, dt_scale=1.0, n_traverse=1.0)
g2, c2, Gam, ex = TD.analyze_scales(t, lg2)
mk = lambda raw: ND.NondimSpec(case=t["label"], system=raw, reference=lg2.ref, ledger=lg2,
        groups=g2, checks=c2, params={"A": t["A"]},
        numerics={"dt_star": lg2.ratio("times","dt","tau_B"), "n_prod": 1}, nhex=12).run_id()
base = mk(t["_raw"])
alt = copy.deepcopy(t["_raw"]); alt["particle"]["diameter"]["value"] = 0.5   # 5µm → 0.5µm
res.append(("trap-drag: run_id changes when d is scaled 10x", base != mk(alt)))
doc = copy.deepcopy(t["_raw"]); doc["description"] = "comment only, edited"
doc["particle"]["diameter"]["source"] = "provenance wording only, edited"
res.append(("trap-drag: run_id unchanged when only comment/source is edited",
            base == mk(doc)))

# ── (5) does hash verification catch a hand-edited stored spec?
# ★ The run_id is NOT hardcoded -- it changes when the physical system changes
#   (N 300 -> 306 for the commensurate hexagon).
p = max(ROOT.glob("specs/trap-drag-2d-hex300__*.json"), key=lambda f: f.stat().st_mtime)
raw = json.loads(p.read_text()); raw["params"]["A"] = 42.0
tmp = Path(tempfile.mkdtemp()) / "tampered.json"
tmp.write_text(json.dumps(raw))
res.append(("trap-drag: hand-edited spec -> hash mismatch",
            not ND.load(tmp).verify_hash()[0]))
res.append(("trap-drag: the original matches its hash", ND.load(p).verify_hash()[0]))

# ── (6) trap-drag: the commensurate-hexagon guard -- does it stop when the lattice
#        does not fit the periodic box?
#     ⚠️ The first attempt recomputed phi from the lattice and compared, and that
#        **passed identically** (because L_x*L_y is derived FROM phi). It was a check
#        that never fired. The comparison that means something is against **the box
#        L2 wrote down**.
def broken(mut):
    tt = copy.deepcopy(t); mut(tt)
    try:
        TD.build_ledger(tt); return False
    except ValueError:
        return True

from bdbot.units import Q as _Q
res.append(("trap-drag: odd n_y -> rejected",
            broken(lambda x: x.update(n_y=17, N=17*17))))
res.append(("trap-drag: N != n_x*n_y -> rejected", broken(lambda x: x.update(N=300))))
res.append(("trap-drag: phi changed so it disagrees with the box -> rejected",
            broken(lambda x: x.update(phi=0.30))))
res.append(("trap-drag: only the YAML box_length_x touched -> rejected",
            broken(lambda x: setattr(x["box_x"], "value", _Q(129.73, "um")))))
_lg = TD.build_ledger(t); _D = _lg.derived
res.append(("trap-drag: L_x/a_NN is the integer n_x",
            abs(_D["Lx_star"]/_D["a_nn_star"] - _D["n_x"]) < 1e-12))
res.append(("trap-drag: L_y/(sqrt(3)/2 a_NN) is the even n_y",
            abs(_D["Ly_star"]/(3**0.5/2*_D["a_nn_star"]) - _D["n_y"]) < 1e-12 and _D["n_y"] % 2 == 0))

# ── (7) chain-bend: is lambda_max the larger of the two blocks? (adding them makes
#        dt needlessly small)
D = lg.derived
res.append(("chain: lambda_max = max(bending, stretching)",
            abs(D["lam_max"] - max(D["lam_bend"], D["lam_bond"])) < 1e-18))
res.append(("chain: bending is faster than stretching", D["lam_bend"] > D["lam_bond"]))

print("="*72); print("adversarial checks -- L3 guards of the two new cases"); print("="*72)
for name, good in res: print(f"  {ok(good):<10} {name}")
n = sum(1 for _, g_ in res if g_)
print("="*72); print(f"{'✓ PASS' if n==len(res) else '✗ FAIL'}  {n}/{len(res)}"); print("="*72)
sys.exit(0 if n == len(res) else 1)
