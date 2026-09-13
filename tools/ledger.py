"""The corpus ledger — what the pipeline has actually *decided*, from artifacts.

    $PY tools/ledger.py                    # the tables
    $PY tools/ledger.py --miss             # every recorded implementation_check failure
    $PY tools/ledger.py --json out.json

Reads `runs/*/record.json` (schema `bdbot.record/0.1`) and renders a verdict for
every observable. Nothing here is typed by hand; every number is a return value.

**Two threshold paths, and no invented third.** Some observables carry
`(err_pct, tol_pct)`, others `(err_sigma, tol_sigma)`. The percent path is
*undefined* when the prediction is exactly zero, which is the case for this
repository's central result — `chain-bend-2d-dlvo`'s `K_prime` is predicted at
**0** from `U'(l)=0`, so it reports in sigma and carries `err_pct: null`. A
percent-only reader therefore sees nothing at all for 132 of the most important
records. Both paths are read here.

An observable with a prediction and **no** threshold in either path is reported
as `no-threshold`, not quietly passed and not assigned a default. Rule 3: a
number that is not known blocks; inventing the tolerance would be the failure.

**Why `role` decides the consequence, not the verdict.** Rule 7'
(CLAUDE.md): a mismatch on an `implementation_check` is a **bug**, on a
`hypothesis` it is a **result**. So `--miss` lists the failures by role -- an
`implementation_check` MISS is an open item whether or not anyone reopened it.

**The pathology this surfaces.** A 3-sigma rule applied to a measurement precise
to a fraction of a percent turns a small model imperfection into a large-sigma
FAIL. Rows where |err_pct| is small and |err_sigma| is large are flagged
`tol-vs-precision`: the disagreement is with the tolerance convention, not
necessarily with the physics. docs/05 section 3 is the general form -- *"a
metric's discriminating power depends on the protocol."*
"""
from __future__ import annotations

import argparse
import collections
import glob
import json
import pathlib
import sys

ROOT = pathlib.Path(__file__).resolve().parent.parent

#: |err_pct| below this with a sigma-path MISS is a tolerance/precision clash,
#: not a large physical disagreement. Reporting threshold only -- it changes no
#: verdict.
SMALL_PCT = 5.0


def load():
    out = []
    for f in sorted(glob.glob(str(ROOT / "runs/*/record.json"))):
        try:
            r = json.load(open(f))
        except Exception as exc:                      # a corrupt record is a row
            out.append((dict(case="?", run_id=pathlib.Path(f).parent.name,
                             broken=f"{type(exc).__name__}: {exc}"), None))
            continue
        obs = r.get("observables") or []
        if not obs:
            out.append((r, None))
        for o in obs:
            out.append((r, o))
    return out


def verdict(o) -> tuple[str, str]:
    """-> (verdict, path). Never invents a threshold."""
    if o is None:
        return "no-observable", "-"
    if o.get("predicted") is None:
        return "no-prediction", "-"
    ep, tp = o.get("err_pct"), o.get("tol_pct")
    es, ts = o.get("err_sigma"), o.get("tol_sigma")
    if ep is not None and tp is not None:
        return ("PASS" if abs(ep) <= tp else "MISS"), "pct"
    if es is not None and ts is not None:
        return ("PASS" if abs(es) <= ts else "MISS"), "sigma"
    return "no-threshold", "-"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--miss", action="store_true")
    ap.add_argument("--json", default=None)
    a = ap.parse_args()

    rows = load()
    recs = {r["run_id"] for r, _ in rows if r.get("run_id")}
    broken = [r for r, _ in rows if r.get("broken")]
    obs = [(r, o) for r, o in rows if o is not None]

    V = collections.Counter()
    RV = collections.Counter()
    PATH = collections.Counter()
    CASE = collections.Counter()
    misses = []
    for r, o in obs:
        v, path = verdict(o)
        V[v] += 1
        RV[(o.get("role"), v)] += 1
        if v in ("PASS", "MISS"):
            PATH[path] += 1
        CASE[(r.get("case"), v)] += 1
        if v == "MISS":
            misses.append((r, o, path))

    print(f"records {len(recs)}   observables {len(obs)}"
          + (f"   UNREADABLE {len(broken)}" if broken else ""))
    print()
    print("| verdict | n | what it means |")
    print("|---|---:|---|")
    MEAN = {
        "no-prediction": "the simulation is the answer; no verdict exists to give",
        "PASS": "prediction and measurement agree inside a stated tolerance",
        "MISS": "they do not -- role decides whether that is a bug or a result",
        "no-threshold": "a prediction was written and no tolerance was; undecidable",
    }
    for k in ("no-prediction", "PASS", "MISS", "no-threshold"):
        if V[k]:
            print(f"| `{k}` | {V[k]} | {MEAN[k]} |")
    dec = V["PASS"] + V["MISS"]
    print()
    print(f"decidable: **{dec} of {len(obs)}** "
          f"({100*dec/max(len(obs),1):.1f} %)   "
          f"via percent {PATH['pct']} · via sigma {PATH['sigma']}")
    print()

    print("| role | no-prediction | PASS | MISS | no-threshold | rule 7' consequence "
          "of a MISS |")
    print("|---|---:|---:|---:|---:|---|")
    CONS = {"implementation_check": "**a bug** -> fix it",
            "hypothesis": "**a result** -> report it",
            "measurement": "n/a -- no prediction",
            None: "⚠ no role recorded"}
    roles = sorted({k[0] for k in RV}, key=lambda x: (x is None, str(x)))
    for role in roles:
        print(f"| `{role}` | " + " | ".join(
            str(RV[(role, v)] or "·") for v in
            ("no-prediction", "PASS", "MISS", "no-threshold"))
            + f" | {CONS.get(role, '?')} |")
    print()

    hyp = sum(n for (role, _), n in RV.items() if role == "hypothesis")
    print(f"⚠ `hypothesis` observables in the whole corpus: **{hyp}**. Rule 7' — "
          f"*\"if the second list is empty the case can only validate, never "
          f"discover.\"*")
    print()

    print("| case | observables | decidable | MISS |")
    print("|---|---:|---:|---:|")
    for case in sorted({k[0] for k in CASE}):
        tot = sum(n for (c, _), n in CASE.items() if c == case)
        d = CASE[(case, "PASS")] + CASE[(case, "MISS")]
        print(f"| `{case}` | {tot} | {d} | {CASE[(case,'MISS')] or '·'} |")

    if misses:
        print()
        print(f"### {len(misses)} recorded MISS")
        print()
        print("| case | observable | n | role | measured | predicted | err | tol | flag |")
        print("|---|---|---:|---|---:|---:|---:|---:|---|")
        # ★ group on the *prediction*, not only the name. `K_prime` MISSes span
        #   predicted=0 (DLVO, percent undefined), 9628 (JKR soft trap, +248 to
        #   +422 %) and 1.113e5 (JKR stiff trap, +1.3 % at 15-29σ). One row for
        #   all three read as "3.0-28.6σ (+1.18-+421.57 %)", which is three
        #   opposite conclusions averaged into noise.
        def gkey(r, o, path):
            return (r.get("case"), o.get("name"), path,
                    f"{o.get('predicted'):.4g}")
        seen = collections.Counter()
        for r, o, path in misses:
            seen[gkey(r, o, path)] += 1
        for k, n in sorted(seen.items(), key=lambda x: (x[0][0], -x[1])):
            case, name, path, _pred = k
            grp = [o for r, o, p in misses if gkey(r, o, p) == k]
            ex = grp[0]
            # ★ a group spans a range -- one example hid 3.2σ and 28.6σ in the
            #   same row. Report the span, not a representative.
            key = "err_sigma" if path == "sigma" else "err_pct"
            mags = sorted(abs(g[key]) for g in grp if g.get(key) is not None)
            pcs = [g["err_pct"] for g in grp if g.get("err_pct") is not None]
            unit = "σ" if path == "sigma" else " %"
            err = (f"{mags[0]:.1f}–{mags[-1]:.1f}{unit}"
                   if len(mags) > 1 and mags[0] != mags[-1]
                   else f"{mags[0]:.1f}{unit}" if mags else "?")
            tol = (f"{ex.get('tol_sigma')}σ" if path == "sigma"
                   else f"{ex.get('tol_pct')} %")
            pc = ""
            if path == "sigma" and pcs:
                lo, hi = min(pcs, key=abs), max(pcs, key=abs)
                pc = (f" ({lo:+.2f}–{hi:+.2f} %)" if lo != hi
                      else f" ({lo:+.2f} %)")
            # flag on the *smallest* percent disagreement in the group: if even
            # one row is a large-sigma MISS at a sub-SMALL_PCT percent error,
            # the convention is what is being refuted there.
            flag = ("`tol-vs-precision`"
                    if path == "sigma" and pcs
                    and abs(min(pcs, key=abs)) < SMALL_PCT else "")
            preds = {g.get("predicted") for g in grp}
            pred = (f"{ex.get('predicted'):.4g}" if len(preds) == 1
                    else f"{min(preds):.4g}–{max(preds):.4g}")
            meas = sorted(g["measured"] for g in grp if g.get("measured") is not None)
            ms = (f"{meas[0]:.4g}–{meas[-1]:.4g}"
                  if len(meas) > 1 and meas[0] != meas[-1]
                  else f"{meas[0]:.4g}" if meas else "?")
            print(f"| `{case}` | {name} | {n} | `{ex.get('role')}` | "
                  f"{ms} | {pred} | {err}{pc} | {tol} | {flag} |")
        print()
        print("`tol-vs-precision` — the percent disagreement is small and the sigma "
              "disagreement is large, so what is being refuted is the 3σ convention "
              "applied to a high-precision measurement, not the model. That "
              "distinction is not recorded anywhere in the artifact; it has to be "
              "made by a human, and nothing currently requires it.")

    if a.miss:
        print()
        for r, o, path in sorted(misses, key=lambda m: -abs(
                m[1].get("err_sigma") or m[1].get("err_pct") or 0)):
            print(f"{r['run_id']}")
            print(f"    {o.get('name')}  role={o.get('role')}  path={path}")
            print(f"    measured={o.get('measured'):.6g}  "
                  f"predicted={o.get('predicted'):.6g}  "
                  f"err_pct={o.get('err_pct')}  err_sigma={o.get('err_sigma')}")
            print(f"    source={o.get('prediction_source')}")

    if broken:
        print()
        print(f"⚠ {len(broken)} unreadable record(s): "
              + ", ".join(f"{b['run_id']} ({b['broken']})" for b in broken[:5]))

    if a.json:
        pathlib.Path(a.json).write_text(json.dumps(dict(
            records=len(recs), observables=len(obs),
            verdicts=dict(V), by_role={f"{k[0]}|{k[1]}": v for k, v in RV.items()},
            paths=dict(PATH), hypothesis=hyp,
            miss=[dict(case=r.get("case"), run_id=r.get("run_id"),
                       name=o.get("name"), role=o.get("role"), path=p,
                       measured=o.get("measured"), predicted=o.get("predicted"),
                       err_pct=o.get("err_pct"), err_sigma=o.get("err_sigma"),
                       tol_pct=o.get("tol_pct"), tol_sigma=o.get("tol_sigma"))
                  for r, o, p in misses],
        ), indent=1, ensure_ascii=False))
        print(f"\nwrote {a.json}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
