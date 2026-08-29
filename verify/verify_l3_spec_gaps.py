"""**Regression test** for four L3 spec defects.

CLAUDE.md rule 6: claims about physics or tool behaviour are settled by execution,
not by reasoning. This script was originally written to **reproduce** the defects,
and all four reproduced (2026-08-04):

  (1) the run_id did not cover the physical system -- the 1-B spec carried no
      physical system, so changing d from 5µm to 0.5µm, eta from water to glycerol
      (62x) and rho_p from silica to polystyrene left the run_id **exactly the
      same**, `soft-r3-2d-A-sweep__A100__27f70deab9`, despite tau_B differing by
      16.1x. That is the name of an already-completed run, so `prepare_outdir`
      skips it and **reports the old system's results as the new system's**.
  (2) inverting back to physical units from the spec alone was impossible (the SI
      values of sigma, tau and kT were absent).
  (3) the three cases' spec.json schemas differed (only 3 keys in common:
      N, n_eq, n_prod).
  (4) there was no way to check that a dimensionless group really is the ratio of
      two ledger entries -> `verify_nondim_guards.py`

This confirms (1), (2) and (3) disappeared once `bdbot.nondim.NondimSpec` was
introduced. **Passing (= no defect) is now the correct outcome.**

    $PY scratch/verify_l3_spec_gaps.py
"""
from __future__ import annotations

import copy
import math
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import nondim as ND  # noqa: E402

sys.path.insert(0, str(ROOT / "cases"))
import soft_r3_2d as S3  # noqa: E402

FAILS = []
# ★ The marker the case scripts print. It lives in ONE place because it is a
#   cross-file coupling: cases/*.py print it, this script parses it. It used to be
#   the Korean "L3 스펙:" inline in two places, and when cases/ was translated in
#   3d68b57 this script started reporting "the case does not write a spec" for all
#   three cases -- a misleading diagnosis, since they do write one.
SPEC_MARKER = "L3 spec:"
LEGACY_RUN_ID = "soft-r3-2d-A-sweep__A100__27f70deab9"   # the value at the time of defect (1)


def check(name: str, ok: bool, detail: str) -> None:
    print(f"  {'✓' if ok else '✗'} {name}\n      {detail}")
    if not ok:
        FAILS.append(name)


def soft_spec(raw: dict, A=100.0, N=400, rc_shells=5.0, dt_scale=1.0, samples=400):
    """Reproduce the **current** spec-construction path of `cases/soft_r3_2d.py`
    main(), verbatim."""
    tmp = ROOT / "verify" / "_tmp_sys.yaml"
    tmp.write_text(yaml.safe_dump(raw))
    sys_ = S3.load_system(tmp)
    tmp.unlink()

    phi = sys_["phi"]
    a_star = math.sqrt(math.pi / (4 * phi))
    r_c_star = rc_shells * a_star
    T_obs_tau = float(sys_["numerics"]["production_tau_B"])

    lg = S3.build_ledger(sys_, A, N, phi, r_c_star, dt_scale, T_obs_tau)
    tau_B = lg.derived["tau_B"]
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    n_eq = int(round(float((0.2 * T_obs / dt).to(""))))
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, n_prod // samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, Gamma = S3.analyze_scales(sys_, lg, A, phi, r_c_star)
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"A": A, "phi": phi, "N": N, "r_c_star": r_c_star,
                "wca_eps": sys_["wca_eps_kT"], "Gamma": Gamma},
        numerics={"dt_star": float((dt / tau_B).to("")),
                  "dt_over_tau_int": dt_scale * 1e-2, "n_eq": n_eq, "n_prod": n_prod,
                  "n_samples": samples, "sample_every": sample_every, "seed": 20260803},
        tag=f"A{A:g}", nhex=10)
    return spec, dt.to("s"), tau_B.to("s")


def main() -> int:
    print("=" * 78)
    print("L3 spec defect regression test -- soft-r3-2d-A-sweep")
    print("=" * 78)

    raw = yaml.safe_load((ROOT / "intake" / "soft-r3-2d-A-sweep" / "system.yaml").read_text())
    spec_a, dt_a, tauB_a = soft_spec(raw)
    id_a = spec_a.run_id()

    other = copy.deepcopy(raw)
    other["particle"]["diameter"]["value"] = 0.5
    other["medium"]["viscosity"]["value"] = 53.0
    other["medium"]["name"] = "glycerol-water 60%"
    other["particle"]["density"]["value"] = 1050
    spec_b, dt_b, tauB_b = soft_spec(other)
    id_b = spec_b.run_id()

    # ── (1) does the run_id cover the physical system? ───────────────────
    print("\n[1] does the run_id cover the physical system (d, eta, rho_p)?")
    print(f"    system A: d=5.0µm  eta=0.851mPa*s  ->  tau_B={tauB_a:~.4gP}  "
          f"dt={dt_a.to('ms'):~.4gP}")
    print(f"    system B: d=0.5µm  eta=53mPa*s   ->  tau_B={tauB_b:~.4gP}  "
          f"dt={dt_b.to('ms'):~.4gP}")
    print(f"    physical time factor: tau_B differs by {float(tauB_a/tauB_b):.1f}x")
    check("a different physical system gives a different run_id", id_a != id_b,
          f"A={id_a}  B={id_b}")
    check("the spec records the physical system", bool(spec_a.to_json().get("system")),
          f"system.label = {spec_a.to_json()['system'].get('label')}")
    check("the run_id from the time of the defect is no longer produced",
          LEGACY_RUN_ID not in (id_a, id_b),
          f"legacy {LEGACY_RUN_ID} -> current {id_a} "
          f"(the 7 existing runs are preserved under the legacy id)")

    print("\n[1b] documentation fields do not change the run_id "
          "(the condition for content addressing)")
    doc = copy.deepcopy(raw)
    doc["particle"]["diameter"]["source"] = "provenance wording only, edited"
    doc["description"] = "description only, edited"
    check("run_id survives editing source and description",
          soft_spec(doc)[0].run_id() == id_a,
          f"{soft_spec(doc)[0].run_id()} vs {id_a}")

    # ── (2) the three inversion anchors ─────────────────────────────────
    print("\n[2] can a result be returned to physical units from the spec alone?")
    bt = spec_a.to_json()["back_transform"]
    need = ("sigma_SI", "tau_SI", "energy_SI")
    check("the three inversion anchors are in the spec", all(k in bt for k in need),
          "  ".join(f"{k}={bt[k]:.4g} {bt[k+'_unit']}" for k in need if k in bt))
    d_si = float(raw["particle"]["diameter"]["value"]) * 1e-6
    check("sigma_SI matches the actual particle diameter",
          abs(bt["sigma_SI"] - d_si) / d_si < 1e-12,
          f"{bt['sigma_SI']:.6e} m vs {d_si:.6e} m")

    # ── (3) a schema shared by the three cases ──────────────────────────
    print("\n[3] do the three cases use the same schema?")
    keysets = {}
    for name, script, extra in (("trap-2d-5um", "trap_2d_5um", []),
                                ("soft-r3", "soft_r3_2d", ["--A", "100"]),
                                ("abp-rod", "abp_rod_2d", [])):
        import subprocess
        out = subprocess.run([sys.executable, str(ROOT / "cases" / f"{script}.py"),
                              *extra, "--spec"], capture_output=True, text=True, cwd=ROOT)
        line = [x for x in out.stdout.splitlines() if SPEC_MARKER in x]
        if not line:
            check(f"{name} writes a spec", False,
                  out.stdout[-300:] + out.stderr[-300:])
            keysets[name] = set()
            continue
        p = ROOT / line[0].split(SPEC_MARKER)[1].strip()
        ls = ND.load(p)
        keysets[name] = set(ls.raw)
        ok, want = ls.verify_hash()
        print(f"    {name:<12} {p.name}   hash check {'✓' if ok else '✗ ' + want}")
        if not ok:
            FAILS.append(f"{name} hash check")
    common = set.intersection(*keysets.values()) if all(keysets.values()) else set()
    required = {"schema", "reference", "back_transform", "ledger", "groups", "checks",
                "system", "params", "numerics", "run_id", "verdict"}
    check("the three cases share a common schema", required <= common,
          f"{len(common)} common keys, missing {sorted(required - common) or 'none'}")

    print("\n" + "=" * 78)
    if FAILS:
        print(f"✗ FAIL -- {len(FAILS)} defect(s):")
        for f in FAILS:
            print(f"   ✗ {f}")
        print("=" * 78)
        return 1
    print("✓ PASS -- defects (1), (2) and (3) are all resolved.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
