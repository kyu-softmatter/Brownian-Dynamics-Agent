"""soft-r3 S30 — the melting hysteresis bracket, and this repository's first
pre-registered run.

    $PY campaigns/soft2d_hysteresis.py --prepare   # specs, params.json, SEAL
    $PY campaigns/soft2d_hysteresis.py --run       # execute, seal REQUIRED
    $PY campaigns/soft2d_hysteresis.py --status

## Why it is in two steps

`--prepare` writes `prediction.yaml` and `analysis_plan.yaml` into each run
directory and seals them; `--run` then refuses to start without that seal
(`require_seal=True`) and without an approved `params.json` (`require_approval=
True`, rule 10). Splitting them is the point: the prediction is on disk and
hashed before any trajectory exists, and the hash is checked by CI thereafter.

Until 2026-09-13 this was impossible. `RID.prepare_outdir` deleted every file
except `record.json` and it ran BEFORE the seal check, so a valid seal was
removed and then reported missing — `require_seal=True` could only refuse,
never be satisfied, and a TAMPERED seal was downgraded to "[warn] unsealed" and
the run proceeded. 0 of 260 archived run directories carry a seal. Fixed in
63fde87; this campaign is the first use.

## The experiment

`hexwin` measured Gamma_Zahn = 57.9 — the centre of Zahn's hexatic window —
from RANDOM starts and read isotropic liquid. Correction 3 of
`findings/order-parameter-magnitude-cannot-identify-a-phase` argues that is an
upper bound, not the boundary: a random start can sit behind a nucleation
barrier arbitrarily long, and a supercooled liquid is also steady. This runs the
other arm — a perfect hexagonal crystal at the same coupling — plus the matched
random arm in the SAME box with the SAME seeds, so the only difference is the
initial positions.

⚠ THREE conventions for the coupling, and this repository has confused them four
times. Every amplitude below is `A_case` — the case script's own, where
`U(r)/kT = A (d/r)^3` with `d` the particle DIAMETER:

    Gamma_case = A_case * 0.2974865      (what metrics.json calls Gamma)
    A_campaign = Gamma_case              (what simbot/campaigns call "A")
    Gamma_Zahn = A_case * 1.6565023      (pi^1.5 * beta U(a_mean))
    Zahn window  Gamma_Z 55.87-59.88  <->  A_case 33.728-36.148
"""
from __future__ import annotations

import argparse
import json
import pathlib
import shutil
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import params as PRM          # noqa: E402
from bdbot import runcard as RC          # noqa: E402

CASE = ROOT / "cases/soft_r3_2d.py"
DOCS = pathlib.Path(__file__).resolve().parent / "s30_preregistration"
SEEDS = (20260914, 20260915, 20260916)

#: `A_case`. See the module docstring for the other two conventions.
A_PRIMARY = 34.938        # Gamma_Zahn 57.87 -- the centre of Zahn's window
A_CRYSTAL = 106.224       # Gamma_Zahn 175.96 -- deep crystal control
A_LIQUID = 10.000         # Gamma_Zahn 16.57  -- deep liquid control

RUNS = (
    [dict(id=f"H{i+1}", A=A_PRIMARY, tobs=100.0, init="hex", seed=s)
     for i, s in enumerate(SEEDS)]
    + [dict(id=f"R{i+1}", A=A_PRIMARY, tobs=100.0, init="rsa", seed=s)
       for i, s in enumerate(SEEDS)]
    + [dict(id="C1", A=A_CRYSTAL, tobs=10.0, init="hex", seed=SEEDS[0]),
       dict(id="C2", A=A_LIQUID, tobs=10.0, init="hex", seed=SEEDS[0])]
)

COMMON = ["--N", "400", "--rc-shells", "7.8", "--box", "hex"]
APPROVED_BY = "Takuya Kobayashi, 2026-09-14 (chat: 'run it')"


def case_args(r) -> list[str]:
    return [*COMMON, "--A", f"{r['A']:g}", "--tobs", f"{r['tobs']:g}",
            "--init", r["init"], "--seed", str(r["seed"])]


def run_id_of(r) -> str:
    out = subprocess.run([sys.executable, str(CASE), *case_args(r), "--spec"],
                         capture_output=True, text=True, cwd=ROOT)
    if out.returncode:
        raise SystemExit(f"{r['id']}: --spec failed\n{out.stdout[-2000:]}\n{out.stderr[-2000:]}")
    line = next(ln for ln in out.stdout.splitlines() if "run_id=" in ln)
    return line.split("run_id=")[1].strip()


def manifest_for(r, run_id: str) -> PRM.Manifest:
    """Rule 10, built FROM THE SPEC. No number here is typed twice: every value
    is read back out of `specs/<run_id>.json`, so an approved manifest that
    disagrees with the spec is impossible to write by hand -- and
    `diff_against_spec` still checks it at run time."""
    spec = json.loads((ROOT / "specs" / f"{run_id}.json").read_text())
    p, nm, sysd = spec["params"], spec["numerics"], spec["system"]
    # ★ ledger entries are keyed `symbol`, not `name`. Keying on "name"
    #   raised KeyError here and `prepare()` died AFTER copying the two
    #   documents into the run directory -- leaving an UNSEALED
    #   pre-registration sitting where a sealed one belongs.
    led = {g["symbol"]: g for g in spec["ledger"]["lengths"]}
    m = PRM.Manifest(case=spec["case"])

    def si(node):
        return node["value"], node["unit"], node["source"][:70], node.get("tier", 3)

    v, u, s, t = si(sysd["particle"]["diameter"])
    m.add("geometry", "d", v, u, s, tier=t)
    m.add("geometry", "N", int(p["N"]), "1", "the spec", tier=1)
    m.add("geometry", "dim", 2, "1", "2D: psi_6 and the eta_6 = 1/4 boundary are 2D constructs", tier=0)
    m.add("geometry", "L", float(led["L"]["value"]), "m",
          "a_mean*sqrt(N); the commensurate box has the same area", tier=1)
    m.add("geometry", "n_x", int(p["n_x"]), "1", "commensurate hexagonal tiling", tier=0)
    m.add("geometry", "n_y", int(p["n_y"]), "1", "commensurate hexagonal tiling; must be even", tier=0)
    m.add("geometry", "phi", float(p["phi"]), "1", "confirmed by the user", tier=1)

    v, u, s, t = si(sysd["medium"]["temperature"])
    m.add("energy", "T", v, u, s, tier=t)
    m.add("energy", "kT", 1.0, "1", "reduced units: kT is the energy scale", tier=0)
    v, u, s, t = si(sysd["medium"]["viscosity"])
    m.add("medium", "eta", v, u, s, tier=t)
    v, u, s, t = si(sysd["particle"]["density"])
    m.add("medium", "rho_p", v, u, s, tier=t)

    m.add("interaction", "A", float(p["A"]), "1",
          "A_case: U/kT = A (d/r)^3, d = particle DIAMETER. Gamma_Zahn = 1.6565023*A", tier=0)
    m.add("interaction", "r_c", float(p["r_c_star"]), "d",
          "7.8 a_mean; beta U(r_c) = 0.0219 kT", tier=1)
    m.add("interaction", "wca_eps", float(p["wca_eps"]), "kT", "excluded-volume core", tier=1)

    m.add("numerics", "dt", float(nm["dt_star"]), "1",
          "dt/tau_B; set from tau_int(r_min), dt/tau_int = 0.01", tier=0)
    m.add("numerics", "n_eq", int(nm["n_eq"]), "1", "20% of the window", tier=1)
    m.add("numerics", "n_prod", int(nm["n_prod"]), "1", "T_obs/dt", tier=1)
    m.add("numerics", "sample_every", int(nm["sample_every"]), "1", "400 samples", tier=1)
    m.add("numerics", "seed", int(nm["seed"]), "1", "fixed per run; H and R share seeds", tier=0)

    # ★ Real numbers, and one of them is an independent cross-check.
    #   `tau_gov` comes from the L3 LEDGER (`spec.ledger.times.tau_B`), which was
    #   computed by `build_ledger`, so `Manifest.check_derived` comparing it
    #   against its own recomputation is a genuine two-source agreement.
    #   `gamma` and `D_t` are not in the ledger and are recomputed here from the
    #   same closed form the check uses, so for THOSE two the check protects
    #   against a placeholder, a typo and a wrong unit, but not against a bug in
    #   `materials.sphere_bulk` itself. That is covered by tests/test_s0_units.py
    #   and by a hand cross-check recorded in tests/test_params.py (15 s.f.).
    #
    #   ⚠ This block used to be three literal 0.0 values with the provenance
    #     "in the ledger, recomputed by L3" -- and `blockers()` returned []. The
    #     rule-10 gate checked that the key was present, never that the value
    #     was real, in the manifest written to demonstrate the rule-10 gate.
    #     See bdbot/params.py "four for four".
    from bdbot import materials as MAT
    from bdbot.units import Q
    d_q = Q(float(sysd["particle"]["diameter"]["value"]), sysd["particle"]["diameter"]["unit"])
    T_q = Q(float(sysd["medium"]["temperature"]["value"]), sysd["medium"]["temperature"]["unit"])
    eta_q = Q(float(sysd["medium"]["viscosity"]["value"]), sysd["medium"]["viscosity"]["unit"])
    bulk = MAT.sphere_bulk(d_q, T_q, eta_q)
    tau_B_ledger = next(e for e in spec["ledger"]["times"] if e["symbol"] == "tau_B")

    m.add("derived", "gamma", float(bulk["gamma"].to("kg/s").magnitude), "kg/s",
          "3 pi eta d (Stokes); recomputed from d, T, eta above", tier=1)
    m.add("derived", "D_t", float(bulk["D_t"].to("m^2/s").magnitude), "m^2/s",
          "kT/gamma; recomputed from d, T, eta above", tier=1)
    m.add("derived", "tau_gov", float(tau_B_ledger["value"]), "s",
          "tau_B governs here -- unlike the trap case, dt is set from "
          "tau_int(r_min) but tau_B is the governing scale (case docstring)", tier=1)
    #  ...and tau_B separately, because THAT one is recomputable and therefore
    #  checked against d, T and eta. tau_gov is only checked for being a
    #  positive time -- which timescale governs is a case-dependent choice.
    m.add("derived", "tau_B", float(tau_B_ledger["value"]), "s",
          "d^2/D_t, taken from the L3 ledger (cross-checks the recomputation)", tier=1)

    m.approved_by = APPROVED_BY
    return m


def prepare() -> int:
    if not (DOCS / "prediction.yaml").exists():
        raise SystemExit(f"missing {DOCS}/prediction.yaml")
    # ⛔ The document says whether it may be sealed, and this refuses to seal one
    #    that says no. Four structural objections from the pre-seal review are
    #    still open (see the header of prediction.yaml). A pre-registration is
    #    permanent once hashed, so the refusal lives here rather than in a note.
    head = (DOCS / "prediction.yaml").read_text()
    if "sealable: false" in head:
        raise SystemExit(
            "prediction.yaml says `sealable: false` -- refusing to seal it.\n"
            "Open objections are listed in its header (K1-K4). Once sealed the "
            "document cannot be edited, so this is the last gate there is.")
    print(f"{len(RUNS)} runs\n")
    # ★ Build and validate every manifest FIRST, writing nothing. A crash in
    #   the middle of the loop used to leave one run directory holding UNSEALED
    #   copies of prediction.yaml and analysis_plan.yaml -- an unsealed
    #   pre-registration exactly where a sealed one is supposed to be, which is
    #   worse than no file at all.
    plan = []
    for r in RUNS:
        rid = run_id_of(r)
        man = manifest_for(r, rid)
        blockers = man.blockers()
        if blockers:
            raise SystemExit(f"{r['id']}: rule 10 blocks:\n  " + "\n  ".join(blockers))
        plan.append((r, rid, man))
    print("  all 8 manifests build and clear rule 10; writing\n")

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
            raise SystemExit(f"{r['id']}: the seal does not verify immediately after writing")
    print("\nprepared. `--run` now REFUSES to start without these seals.")
    return 0


def run() -> int:
    for r in RUNS:
        print(f"\n{'=' * 70}\n{r['id']}  A_case={r['A']}  init={r['init']}  "
              f"seed={r['seed']}  T_obs={r['tobs']} tau_B\n{'=' * 70}")
        cmd = [sys.executable, str(CASE), *case_args(r),
               "--require-seal", "--require-approval"]
        rc = subprocess.run(cmd, cwd=ROOT).returncode
        if rc:
            print(f"\n{r['id']} exited {rc} -- stopping. Nothing after it has run.")
            return rc
    return 0


def status() -> int:
    print(f"{'id':4s} {'run_id':52s} {'seal':6s} {'params':7s} {'metrics':8s}")
    for r in RUNS:
        rid = run_id_of(r)
        d = ROOT / "runs" / rid
        ok = RC.verify_seal(d, root=ROOT)[0] if (d / RC.SEAL_NAME).exists() else False
        print(f"{r['id']:4s} {rid:52s} "
              f"{'OK' if ok else '-':6s} "
              f"{'yes' if (d / 'params.json').exists() else '-':7s} "
              f"{'yes' if (d / 'metrics.json').exists() else '-':8s}")
    return 0


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--prepare", action="store_true")
    ap.add_argument("--run", action="store_true")
    ap.add_argument("--status", action="store_true")
    a = ap.parse_args()
    if a.prepare:
        return prepare()
    if a.run:
        return run()
    return status()


if __name__ == "__main__":
    sys.exit(main())
