"""S31 — the sedimentation EOS campaign. Seal first, then run.

    PY=./bin/py
    $PY campaigns/s31_sedimentation.py --check     # do the docs determine the runs?
    $PY campaigns/s31_sedimentation.py --prepare   # write + seal the pre-registration
    $PY campaigns/s31_sedimentation.py --run       # run, REFUSED without the seals
    $PY campaigns/s31_sedimentation.py --status

The pre-registration is `campaigns/sediment_preregistration/{prediction,analysis_plan}.yaml`.
`--prepare` copies both into every run directory, writes the rule-10 manifest and
seals all three; `run.execute(require_seal=True, require_approval=True)` then
refuses to start without them.

## Three things this driver does that the obvious version would not

★ **Everything is validated before anything is written.** `--prepare` builds and
checks all eight manifests first. The predecessor campaign crashed mid-loop once
and left one run directory holding UNSEALED copies of the two documents — an
unsealed pre-registration sitting exactly where a sealed one belongs, which is
worse than no file at all.

★ **`--check` parses the sealed documents and compares them against the specs the
runs will actually use.** A pre-registration that does not determine the runs is
decoration. The check is by PARSE, never by substring: the predecessor recorded
that a grep for `sealable: false` matched the prose explaining why a correction
had been made.

★ **The primary arm runs first.** It is also the arm the analysis plan reads
first, so an implementation failure is known before the 4x-cost finite-size arm
has started.
"""
from __future__ import annotations

import argparse
import json
import math
import shutil
import subprocess
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
from bdbot import params as PRM, runcard as RC  # noqa: E402

CASE = ROOT / "cases" / "sediment_3d.py"
DOCS = ROOT / "campaigns" / "sediment_preregistration"
PY = sys.executable


def _prereg() -> dict:
    return yaml.safe_load((DOCS / "prediction.yaml").read_text())


#: `arm` -> the letter the run ids use. Keyed so a new arm has to be named here
#: rather than silently falling into the else-branch of a conditional.
ARM_TAG = {"primary": "P", "convergence": "C", "finite_size": "F"}

#: HOOMD truncates `Simulation(seed=...)` to 16 bits — measured 2026-09-16, see
#: `prediction.yaml: runs.⚠_the_seeds_must_stay_distinct_INSIDE_hoomd`.
SEED_BITS = 0xFFFF


def runs() -> list[dict]:
    """The run table, read OUT OF the sealed document. Typing it twice is how a
    campaign and its pre-registration drift apart.

    `production_tau_sed` and `init` live under `common` as DEFAULTS and may be
    overridden per entry: the arms deliberately differ in both.
    """
    d = _prereg()
    common = d["runs"]["common"]
    out = []
    for i, e in enumerate(d["runs"]["entries"]):
        arm = e["arm"]
        if arm not in ARM_TAG:
            raise SystemExit(f"entry {i}: unknown arm {arm!r}; add it to ARM_TAG")
        out.append({
            "id": f"{ARM_TAG[arm]}{i + 1}",
            "arm": arm, "lxy": float(e["lxy"]), "seed": int(e["seed"]),
            "N": int(e["N"]),
            "tau_sed": float(e.get("tau_sed", common["production_tau_sed"])),
            "init": str(e.get("init", common["init"])),
        })
    #  primary first -- the analysis plan reads it first, so a failure is cheap
    order = {"primary": 0, "convergence": 1, "finite_size": 2}
    return sorted(out, key=lambda r: (order[r["arm"]], r["id"]))


def seed_collisions() -> list[str]:
    """★ Two DISTINCT seeds that HOOMD truncates to the SAME 16-bit value give
    two different run_ids driven by identical thermal noise, and a seed-to-seed
    sd computed over them is silently too small. Checked per (arm, lxy, init),
    which is the group the sd is actually taken over."""
    groups: dict[tuple, dict[int, list[int]]] = {}
    for r in runs():
        key = (r["arm"], r["lxy"], r["init"], r["tau_sed"])
        groups.setdefault(key, {}).setdefault(r["seed"] & SEED_BITS, []).append(r["seed"])
    bad = []
    for key, byt in groups.items():
        for trunc, seeds in byt.items():
            if len(set(seeds)) > 1:
                bad.append(f"{key}: seeds {sorted(set(seeds))} all become "
                           f"{trunc} inside HOOMD")
    return bad


def case_args(r: dict) -> list[str]:
    return ["--lxy", f"{r['lxy']:g}", "--seed", str(r["seed"]),
            "--tau-sed", f"{r['tau_sed']:g}", "--init", r["init"]]


def run_id_of(r: dict) -> str:
    """Ask the case, do not reconstruct. The run_id is the content hash of the
    spec, so any reconstruction here is a second implementation that can drift."""
    out = subprocess.run([PY, str(CASE), *case_args(r), "--spec"],
                         cwd=ROOT, capture_output=True, text=True)
    if out.returncode != 0:
        raise SystemExit(f"{r['id']}: --spec failed\n{out.stdout[-2000:]}{out.stderr[-2000:]}")
    for line in out.stdout.splitlines():
        if "run_id=" in line:
            return line.split("run_id=")[1].strip()
    raise SystemExit(f"{r['id']}: no run_id on stdout\n{out.stdout[-2000:]}")


# ── rule 10, built FROM the spec so no number is typed twice ───────────────
def manifest_for(r: dict, run_id: str) -> PRM.Manifest:
    spec = json.loads((ROOT / "specs" / f"{run_id}.json").read_text())
    p, nm, sysd = spec["params"], spec["numerics"], spec["system"]
    led = {g["symbol"]: g for g in spec["ledger"]["lengths"]}
    ledt = ledt_pre = {g["symbol"]: g for g in spec["ledger"]["times"]}
    m = PRM.Manifest(case=spec["case"])

    def si(node):
        return node["value"], node["unit"], node["source"][:70], node.get("tier", 3)

    v, u, s, t = si(sysd["particle"]["diameter"])
    m.add("geometry", "d", v, u, s, tier=t)
    m.add("geometry", "dim", 3, "1",
          "3D: the profile IS the third dimension, and Carnahan-Starling is the "
          "3D hard-sphere EOS", tier=0)
    m.add("geometry", "N", int(p["N"]), "1", "fixed by phi and the box, not chosen", tier=1)
    m.add("geometry", "L_z", float(led["L_z"]["value"]), "m",
          "the paper's stated cell thickness", tier=0)
    #  rule 10's REQUIRED_KEYS asks for `L`, and this cell has two box lengths
    #  that are not interchangeable: L_z is the paper's and L_xy is ours. `L` is
    #  bound to the one the confinement is along, which is the one every gate in
    #  the card's table is written against.
    m.add("geometry", "L", float(led["L_z"]["value"]), "m",
          "= L_z. The confined axis, which is the box length every gate uses: "
          "exp(-L_z/l_g) is the lid weight and L_z l_g/D_0 is the settling time",
          tier=0)
    m.add("geometry", "L_xy", float(led["L_xy"]["value"]), "m",
          "OURS, tier 3, approved by the user 2026-09-16", tier=3)
    m.add("geometry", "phi", float(p["phi"]), "1", "the paper's mean volume fraction", tier=0)
    m.add("geometry", "l_g", float(led["l_g"]["value"]), "m",
          "the paper's measured gravitational length, Table 1", tier=0)

    v, u, s, t = si(sysd["medium"]["temperature"])
    m.add("energy", "T", v, u, s, tier=t)
    m.add("energy", "kT", 1.0, "1", "reduced units: kT is the energy scale", tier=0)
    m.add("energy", "eps_wca", float(p["eps_wca_kT"]), "kT",
          "WCA amplitude. eps=10 is 3x more accurate at phi=0.3 for 1.5x the "
          "cost; 1 is inside 2 % over THIS phi range", tier=1)
    m.add("energy", "eps_wall", float(p["eps_wall_kT"]), "kT",
          "bounded wall. U(0)=20 kT -> penetration ~e^-20", tier=1)

    v, u, s, t = si(sysd["medium"]["viscosity"])
    m.add("medium", "eta", v, u, s, tier=t)
    v, u, s, t = si(sysd["medium"]["density"])
    m.add("medium", "rho_fluid", v, u, s, tier=t)
    v, u, s, t = si(sysd["particle"]["density"])
    m.add("medium", "rho_p", v, u, s, tier=t)

    m.add("interaction", "sigma_wall", float(p["sig_wall_d"]), "d",
          "Gaussian wall width. F_max = eps/(sigma sqrt(e)) is BOUNDED, which is "
          "the point -- an LJ wall reaches 1517 kT/d (bd-hoomd trap 22)", tier=1)
    m.add("interaction", "r_cut_wall", float(p["rcut_wall_d"]), "d",
          "U(2 d) = 4.5e-9 kT, so the analysis windows start clean", tier=1)
    m.add("interaction", "sigma_lj", 2.0 ** (-1.0 / 6.0), "d",
          "so the WCA minimum sits at r = d. sigma_LJ = 1 would make the "
          "effective diameter 1.1225 d and phi 1.414x too large", tier=0)

    #  ⚠ `dt` is recorded in TAU_D, not seconds. `run.execute` compares the
    #     approved manifest against the spec, and `diff_against_spec` matches a
    #     manifest `dt` to the spec's `dt_star` -- which is reduced. The first
    #     version put the SI step here and the gate refused the run, correctly
    #     but with a message about values rather than units. The SI step is kept
    #     beside it under a name that does not collide.
    m.add("numerics", "dt", float(nm["dt_star"]), "tau_d",
          "the step in the run's own reference time, which is the number the "
          "integrator actually receives", tier=1)
    m.add("numerics", "dt_seconds", float(ledt_pre["dt"]["value"]), "s",
          "the same step dimensionally, so it can be read against gamma and "
          "D_t. ⚠ NOT compared against the spec -- see the unit note in "
          "bdbot.params.diff_against_spec", tier=1)
    m.add("numerics", "dt_star", float(nm["dt_star"]), "tau_d",
          "set BY the WCA gate at U = 8 kT, not fitted to it: 2.0e-4 was refused "
          "by the hard check at 1.08x", tier=1)
    m.add("numerics", "n_eq", int(nm["n_eq"]), "1",
          "ZERO by design: the approach to equilibrium IS part of the "
          "measurement (stationarity in the primary arm, convergence in the "
          "ideal arm), so it is sampled rather than discarded blind", tier=1)
    m.add("numerics", "init", 0 if str(p["init"]) == "ideal" else 1, "1",
          f"init={p['init']!r} encoded 0=ideal-gas exponential, 1=CS hydrostatic "
          f"profile. It is in `params`, so it is hashed into the run_id", tier=1)
    m.add("numerics", "seed_hoomd", int(nm["seed"]) & SEED_BITS, "1",
          "★ what HOOMD ACTUALLY uses: Simulation(seed=) truncates to 16 bits. "
          "Recorded because the spec hashes the untruncated value, so two runs "
          "could differ in name and share a noise stream", tier=1)
    m.add("numerics", "n_prod", int(nm["n_prod"]), "1",
          "production_tau_sed x tau_sed / dt, from the sealed document", tier=1)
    m.add("numerics", "sample_every", int(nm["sample_every"]), "1",
          "samples_per_tau_sed from the sealed document", tier=1)
    m.add("numerics", "seed", int(nm["seed"]), "1", "one of the sealed seed list", tier=0)
    m.add("numerics", "min_sep", float(p["min_sep_d"]), "d",
          "RSA rejection radius. Uniform random placement threw a particle 39 d "
          "out of the box on the first thousand steps", tier=1)

    #  ⚠ The first version called `m.not_applicable("interaction", ...)` here and
    #     `bdbot.params` refused it, correctly: the category ALREADY carries
    #     sigma_lj, sigma_wall and r_cut_wall, so declaring it absent would have
    #     been a contradiction. What is absent is not the category -- it is the
    #     TRAP and the BONDS, which are items inside it. Rule 10 wants the
    #     absence STATED, and zero stiffness is exactly what "no trap" means in
    #     this model, so it is stated as a number rather than as a category-level
    #     N/A that would also erase the wall parameters.
    m.add("interaction", "k_trap", 0.0, "kT/d^2",
          "NO TRAP. Nothing is held: the only external force is uniform gravity "
          "and the only confinement is the two walls", tier=0)
    m.add("interaction", "n_bonds", 0, "1",
          "NO BONDS and no angles. 584 independent spheres; the only pair "
          "interaction is the WCA core above", tier=0)

    #  ★ gamma, D_t and tau_B come from `bdbot.materials.sphere_bulk`, the same
    #    function `Manifest.check_derived` recomputes them with -- so the check
    #    is a round-trip through the SI inputs, not a comparison of one typed
    #    number against another.
    from bdbot import materials as MAT
    from bdbot.units import Q
    bulk = MAT.sphere_bulk(
        Q(float(sysd["particle"]["diameter"]["value"]),
          sysd["particle"]["diameter"]["unit"]),
        Q(float(sysd["medium"]["temperature"]["value"]),
          sysd["medium"]["temperature"]["unit"]),
        Q(float(sysd["medium"]["viscosity"]["value"]),
          sysd["medium"]["viscosity"]["unit"]))
    m.add("derived", "gamma", float(bulk["gamma"].magnitude), "kg/s",
          "3 pi eta d, bdbot.materials.sphere_drag", tier=1)
    m.add("derived", "D_t", float(bulk["D_t"].to("m^2/s").magnitude), "m^2/s",
          "kT/gamma, bdbot.materials.diffusion. The paper's measured D_0 is "
          "0.44 um^2/s and eta was DERIVED from it, so this is a round trip and "
          "not independent evidence -- see the disclosure in prediction.yaml",
          tier=1)
    m.add("derived", "tau_B", float(bulk["tau_B"].to("s").magnitude), "s",
          "d^2/D_t. ⚠ NOT the governing timescale here -- tau_sed is", tier=1)

    g = float(ledt["tau_d"]["value"])
    m.add("derived", "tau_d", g, "s", "d^2/D_0, the clock", tier=1)
    m.add("derived", "tau_sed", float(ledt["tau_sed"]["value"]), "s",
          "l_g^2/D_0, the slowest equilibrium mode", tier=1)
    m.add("derived", "tau_fall", float(ledt["tau_fall"]["value"]), "s",
          "L_z l_g/D_0, the DRIFT time across the cell -- 9x tau_sed, and the "
          "reason the run starts from the analytic profile", tier=1)
    m.add("derived", "tau_gov", float(ledt["tau_sed"]["value"]), "s",
          "the governing timescale is tau_sed: it is what T_obs is measured in", tier=1)
    m.approved_by = ("user, 2026-09-16 — the cross-section and the two-arm "
                     "finite-size scope were approved against a cost table")
    return m


def check_docs_determine_the_runs() -> list[str]:
    """Does the sealed pre-registration actually fix what will run? Compared by
    PARSE against the spec each run will use -- never by substring."""
    bad = []
    d = _prereg()
    common = d["runs"]["common"]
    if not d.get("sealable"):
        bad.append("prediction.yaml has `sealable: false` -- refusing to seal it")
    plan = yaml.safe_load((DOCS / "analysis_plan.yaml").read_text())
    if plan.get("sealed_with") != "prediction.yaml":
        bad.append("analysis_plan.yaml does not name prediction.yaml in `sealed_with`")
    if plan.get("campaign") != d.get("campaign"):
        bad.append("the two documents name different campaigns")

    for r in runs():
        rid = run_id_of(r)
        spec = json.loads((ROOT / "specs" / f"{rid}.json").read_text())
        p, nm = spec["params"], spec["numerics"]
        for key, got, want in (
            ("phi", float(p["phi"]), float(common["phi"])),
            ("l_g_star", float(p["l_g_star"]), float(common["l_g_star"])),
            ("L_z_star", float(p["L_z_star"]), float(common["L_z_star"])),
            ("eps_wca_kT", float(p["eps_wca_kT"]), float(common["eps_wca_kT"])),
            ("eps_wall_kT", float(p["eps_wall_kT"]), float(common["eps_wall_kT"])),
            ("sig_wall_d", float(p["sig_wall_d"]), float(common["sig_wall_d"])),
            ("min_sep_d", float(p["min_sep_d"]), float(common["min_sep_d"])),
            ("dt_star", float(nm["dt_star"]), float(common["dt_star"])),
            ("n_eq", int(nm["n_eq"]), int(common["n_eq"])),
            ("seed", int(nm["seed"]), int(r["seed"])),
            ("L_xy_star", float(p["L_xy_star"]), float(r["lxy"])),
            ("N", int(p["N"]), int(r["N"])),
        ):
            if not math.isclose(float(got), float(want), rel_tol=1e-9):
                bad.append(f"{r['id']}: {key} = {got!r} in the spec but {want!r} "
                           f"in the sealed document")
        #  `init` is a string, so it does not go through isclose -- and it is the
        #  one field that distinguishes the two arms, so it is checked explicitly
        if str(p.get("init")) != r["init"]:
            bad.append(f"{r['id']}: init = {p.get('init')!r} in the spec but "
                       f"{r['init']!r} in the sealed document")
        #  ⚠ T_obs is the most cost-determining number in the document, and the
        #     first version of this check was `n_prod * dt_star > 0` -- which is
        #     arithmetically always true and therefore not a check at all. Caught
        #     by dropping `--tau-sed` from `case_args` and watching it pass. The
        #     real comparison is against the LEDGER's tau_sed, in tau_d.
        ledt = {g["symbol"]: g for g in spec["ledger"]["times"]}
        tau_sed_in_tau_d = (float(ledt["tau_sed"]["value"])
                            / float(ledt["tau_d"]["value"]))
        t_obs = int(nm["n_prod"]) * float(nm["dt_star"]) / tau_sed_in_tau_d
        #  n_prod is floored to a whole number of sampling intervals, so this is
        #  not exact; 1e-3 is far tighter than any arm separation (3, 4, 8).
        if not math.isclose(t_obs, r["tau_sed"], rel_tol=1e-3):
            bad.append(f"{r['id']}: T_obs = {t_obs:.6g} tau_sed in the spec but "
                       f"{r['tau_sed']:g} in the sealed document")
    bad += seed_collisions()
    return bad


# ── the three commands ─────────────────────────────────────────────────────
def prepare() -> int:
    for name in RC.SEALED_DOCS:
        if not (DOCS / name).exists():
            raise SystemExit(f"{name} is missing from {DOCS}")
    rs = runs()
    print(f"{len(rs)} runs\n")

    #  ★ validate everything, write nothing
    plan = []
    for r in rs:
        rid = run_id_of(r)
        man = manifest_for(r, rid)
        blockers = man.blockers()
        if blockers:
            raise SystemExit(f"{r['id']}: rule 10 blocks:\n  " + "\n  ".join(blockers))
        plan.append((r, rid, man))
    bad = check_docs_determine_the_runs()
    if bad:
        raise SystemExit("the sealed document does not determine the runs:\n  "
                         + "\n  ".join(bad))
    print(f"  all {len(plan)} manifests clear rule 10, and every number in "
          f"prediction.yaml matches the spec it will run\n")

    for r, rid, man in plan:
        d = ROOT / "runs" / rid
        d.mkdir(parents=True, exist_ok=True)
        for name in RC.SEALED_DOCS:
            shutil.copy2(DOCS / name, d / name)
        man.write(d)
        seal = RC.write_seal(d, root=ROOT)
        ok, problems = RC.verify_seal(d, root=ROOT)
        print(f"  {r['id']:3s} {rid}")
        print(f"      sealed {len(RC.read_seal(d))} doc(s) -> {seal.name}   "
              f"verify: {'OK' if ok else problems}")
        if not ok:
            raise SystemExit(f"{r['id']}: the seal does not verify immediately "
                             f"after writing")
    print("\nprepared. `--run` now REFUSES to start without these seals.")
    return 0


def run() -> int:
    for r in runs():
        print(f"\n{'=' * 70}\n{r['id']}  arm={r['arm']}  L_xy={r['lxy']:g} d  "
              f"N={r['N']}  seed={r['seed']} (hoomd {r['seed'] & SEED_BITS})  "
              f"init={r['init']}  T_obs={r['tau_sed']:g} tau_sed\n"
              f"{'=' * 70}", flush=True)
        cmd = [PY, str(CASE), *case_args(r), "--require-seal", "--require-approval"]
        rc = subprocess.run(cmd, cwd=ROOT).returncode
        if rc:
            print(f"\n{r['id']} exited {rc} -- stopping. Nothing after it has run.")
            return rc
    return 0


def status() -> int:
    print(f"{'id':4s} {'arm':12s} {'init':6s} {'run_id':56s} {'seal':5s} "
          f"{'params':7s} {'metrics':8s}")
    for r in runs():
        rid = run_id_of(r)
        d = ROOT / "runs" / rid
        ok = RC.verify_seal(d, root=ROOT)[0] if (d / RC.SEAL_NAME).exists() else False
        print(f"{r['id']:4s} {r['arm']:12s} {r['init']:6s} {rid:56s} "
              f"{'OK' if ok else '-':5s} "
              f"{'yes' if (d / 'params.json').exists() else '-':7s} "
              f"{'yes' if (d / 'metrics.json').exists() else '-':8s}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--check", action="store_true")
    g.add_argument("--prepare", action="store_true")
    g.add_argument("--run", action="store_true")
    g.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.check:
        bad = check_docs_determine_the_runs()
        for b in bad:
            print(f"  x {b}")
        print(f"\n{len(bad)} problem(s)" if bad else
              "\nthe sealed documents determine every run")
        return 1 if bad else 0
    if a.prepare:
        return prepare()
    if a.run:
        return run()
    return status()


if __name__ == "__main__":
    raise SystemExit(main())
