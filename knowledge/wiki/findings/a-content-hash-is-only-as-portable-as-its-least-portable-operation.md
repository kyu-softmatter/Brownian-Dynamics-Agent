---
type: finding
author: agent
drafted: 2026-09-15
confirmed_by:
system: tooling
dynamics: none
id: a-content-hash-is-only-as-portable-as-its-least-portable-operation
kind: findings
tags: [run-id, content-addressing, reproducibility, floating-point, libm, ci, platform]
created: 2026-09-15
updated: 2026-09-15
confidence: high
supersedes: []
cause_class: tooling
stage: S3
question: "Does the same physical system get the same run_id on a different machine?"
answer: "No. `run_id` hashes floats at full `repr` precision and `params.Gamma = A/a_mean**3` goes through `pow`, which IEEE-754 does not require to be correctly rounded. Apple libm and glibc differ by one ULP, so `soft-r3-2d-A-sweep__A100` is `30caa5c9e0` on osx-arm64 and `079a25f073` on linux-64. `verify_hash()` is unaffected; only re-deriving a spec is."
reproduced: yes
runs:
  - soft-r3-2d-A-sweep__A100__30caa5c9e0
cites:
  - bdbot/runid.py
  - bdbot/pairpot.py
  - tests/test_soft_r3_init.py
  - docs/05-pitfalls.md
  - campaigns/s30_preregistration/prediction.yaml
---

# A content hash is only as portable as its least portable operation

> **The same physical system has two `run_id`s, and which one you get depends on
> the machine.** Rule 2 says `run_id` is the hash of the content, so editing the
> content breaks the hash. What it did not say is that *computing* the content
> twice can break it too.

## What was asked

CI run 35026488825, the first run on `main` after 13 commits merged, failed with
**one** test out of 1428: `tests/test_soft_r3_init.py::test_the_archived_run_id_does_not_move`.

```
expected   soft-r3-2d-A-sweep__A100__30caa5c9e0     (asserted as a literal)
linux-64   soft-r3-2d-A-sweep__A100__079a25f073
```

A fresh clone at the same commit on `osx-arm64` gave `1421 passed, 7 skipped` —
identical collection, 1428 tests, no failure. So the divergence was the
platform, and the question was which byte moved.

## How it was localised

`runid.spec_hash` is `sha256(json.dumps(payload, sort_keys=True))[:nhex]`, so any
single changed byte moves the whole digest and the digest says nothing about
*what* changed. The payload was reconstructed from the archived artefact
`specs/soft-r3-2d-A-sweep__A100__30caa5c9e0.json` (which re-hashes to its own
name — confirming the reconstruction), giving **30 float leaves**. Each was then
perturbed by ±1…6 ULP, one at a time, and the digest recomputed:

```
HIT  params.Gamma: 29.748648790272707 -> 29.74864879027271  (+1 ulp)
```

Exactly one leaf, exactly one ULP, and no other single perturbation within ±6
ULP of any leaf reproduces the linux digest. That is the whole difference.

## The cause

`Gamma = A / a_mean³` with `a_mean = √(π/4φ)`.

| | value | |
|---|---|---|
| `a_mean` | `1.4979969134027407` | **identical on both** — `sqrt` is required by IEEE-754 to be correctly rounded |
| `a_mean**3` | `3.361497213033026` | Apple libm, `osx-arm64` |
| `a_mean*a_mean*a_mean` | `3.3614972130330254` | exact multiplications — **and glibc's `pow` agrees with this** |

★ **`pow` carries no correct-rounding guarantee.** Every `pow`-based route on
`osx-arm64` — `a**3`, `math.pow(a,3)`, `np.asarray(a)**3`, `a**-3.0` — returns
the same value, and plain repeated multiplication returns the value `linux-64`
reports. So the divergence is one library's `pow`, not a numpy build, not a BLAS
difference, and not a version drift: `environment.yml` was untouched in all 13
commits, so micromamba's cache key was unchanged and CI ran the *same* resolved
environment as the previous green run. python, numpy, freud and hoomd versions
are identical; only `scipy` differs (1.18.0 vs 1.18.1) and it is not on this
path.

It was invisible for six weeks because every run in this repository was computed
on one laptop. CI has existed since 2026-08-31 and had run green eight times —
it had simply never been asked to *re-derive* a spec and compare the name.

## What this does and does not break

**Does not.** `LoadedSpec.verify_hash()` re-hashes the **stored** JSON, so it
reproduces the stored name on any platform. Rule 2's hand-edit detector is
intact, and it passed on `linux-64` across all 296 specs.

**Does.** Re-deriving a spec on another machine produces a different `run_id`,
so the archive cannot be reproduced or extended elsewhere under its own names.
104 of 296 specs carry a hashed float `params.Gamma`, and **all 104** move under
a 1-ULP shift (measured, not estimated).

## ★ Why it is pinned and not fixed

Normalising the payload — rounding `Gamma` to 12 significant figures — makes the
derivation portable. It also renames the archive, and that rename is not
available:

- the **sealed** `campaigns/s30_preregistration/prediction.yaml` cites
  `runs/soft-r3-2d-A-sweep__A100__30caa5c9e0` by name;
- its sha256 is locked into **12** `SEALED.sha256` files.

Renaming therefore means editing a pre-registered document after the results are
known, which is the single thing pre-registration exists to prevent. The
trade-off is measured rather than asserted: the mutation "`Gamma` rounded to 12
sig figs" is **CAUGHT** by the replacement test, i.e. it does move the identity.

This follows the precedent already recorded in `spec_hash`'s own docstring,
where the `ensure_ascii` divergence against `simbot.io.sha256_payload` was
deliberately left unified-never because 263 run directories are named by it, and
pinned by a test so it could not be "fixed" by accident.

## How the test was rewritten

The failing assertion was a literal digest, which made a payload test into a
platform test. It now asserts the payload:

- the run_id still starts `soft-r3-2d-A-sweep__A100__`;
- its digest lies in the **one-ULP neighbourhood** of the archived payload —
  every digest reachable by moving one float leaf to an adjacent double;
- both measured digests are pinned as members of that set: `30caa5c9e0`
  (`osx-arm64`) and `079a25f073` (`linux-64`). A libm change that moved it
  *further* fails, which is news rather than noise;
- `pow`'s freedom is asserted portably: `a**3` and `a*a*a` must agree to 15
  significant figures and differ by at most one ULP — true on both platforms.

★ **Guard on the guard.** A tolerance that accepts everything is not a
tolerance, so its discriminating power is measured. Four mutations must land
outside the neighbourhood and all four do:

```
CAUGHT   init written on the default path        (survived this file's v1)
CAUGHT   n_x/n_y written on the default path     (survived this file's v1)
CAUGHT   phi moved in the 10th digit
CAUGHT   A moved in the 12th digit
CAUGHT   Gamma dropped from params
CAUGHT   Gamma rounded to 12 significant figures
```

This is strictly more checking than the literal string did: the string said
"the digest is `30caa5c9e0`", true on one machine; the neighbourhood says "the
payload *is* the archived payload, to within one ULP of one float, and in no
other respect", which holds everywhere.

## The general form

**A content address is only as portable as the least portable operation upstream
of it.** Hash rounded values, or accept that the address is machine-local — but
say which. A content hash advertises "same content, same name"; ours advertises
"same content, same machine, same name", and that had not been written down.

If the archive is ever re-identified — which can only happen at a campaign
boundary, where re-sealing is legitimate — normalise the hashed payload to 15
significant figures at the same time. Both platforms' `Gamma` collapse to the
same 15-digit string, so 15 digits would have been portable from the start.

## See also

`docs/05-pitfalls.md#a-content-hash-over-a-pow-result-is-not-portable` ·
[[a-crystal-start-does-not-survive-at-the-zahn-window-centre]] ·
[[provenance-must-have-one-definition-and-three-capture-points]] ·
[[yaml-scientific-notation-parsed-as-string]]
