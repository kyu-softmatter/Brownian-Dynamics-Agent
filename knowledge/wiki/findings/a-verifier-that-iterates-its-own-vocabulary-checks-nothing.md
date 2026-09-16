---
type: finding
author: agent
drafted: 2026-09-15
confirmed_by:
system: tooling
dynamics: none
id: a-verifier-that-iterates-its-own-vocabulary-checks-nothing
kind: findings
tags: [seal, pre-registration, unwired-checker, silent-skip, tamper, verification, s8]
created: 2026-09-15
updated: 2026-09-15
confidence: high
supersedes: []
cause_class: tooling
stage: S8
question: "Has the seal ever verified an archived run?"
answer: "No. `simbot.io.verify_seal` hashed zero documents on all 21 archived seal directories — 12 of them returning ok=True while claiming two documents were verified unchanged, 9 returning 'never sealed' on intact content. It iterated SEALED_STAGES instead of the seal's own entries. Tampering produced a bit-identical verdict. `bdbot.runcard.verify_seal` was correct throughout; the broken one was the one wired to cli.py, report.py and validate.py."
reproduced: yes
runs:
  - runs_s1s8/2026-07-29_soft-r3-hexwin
  - runs/soft-r3-2d-A-sweep__A34.938-eq0-rc7.8-hex20x20-s20260914__5e60542041
cites:
  - simbot/io.py
  - simbot/report.py
  - bdbot/runcard.py
  - tests/test_s8_io.py
  - tests/test_s8_report.py
  - verify/verify_gates_bite.py
  - .github/workflows/ci.yml
  - docs/05-pitfalls.md
---

# A verifier that iterates its own vocabulary checks nothing

> **The seal is the centrepiece of this repository's claim to pre-registration,
> and it had never verified an archived run.** Not "verified weakly" — hashed
> **zero bytes**, on all 21 sealed directories, while 12 of them printed
> *"봉인 검증 통과 — 2개 문서 실행 후 미변경"*.

## How it was measured

Not by reading the code. By instrumenting `simbot.io.sha256_file` and counting
calls per directory:

```
dir                                    simbot.ok  hashed  entries   runcard.ok  hashed
runs/soft-r3-…-s20260914__5e60542041   True       0       2         True        2   <-- ok=True, hashed NOTHING
…  (12 such)
runs_s1s8/2026-07-29_soft-r3-hexwin    False      0       1         True        1
…  (9 such)
                                       21 of 21 hashed 0.  runcard: 21 of 21 correct.
```

`entries` is the number of **lines in the seal file**, and `summary()` was
reporting it as the number of documents verified. So the pass was not merely
empty, it was *specific about a quantity it had not computed*.

## The cause is one loop

```python
for stage in stages:              # ("prediction", "intake", "spec")
    p = rundir.file(stage)        # 02_prediction.md · 01_intake.md · 03_spec.yaml
    rel = _seal_relpath(p)        # where the document is NOW
    if rel not in sealed:         # the seal records where it WAS
        if p.exists(): unsealed.append(...)
        continue                  # ← content never hashed
```

It iterated **the filename vocabulary this module knows** rather than **the
seal's own entries**, and that is wrong in two independent ways at once:

| | what happens | verdict |
|---|---|---|
| 12 `runs/` (bdbot-era) | seals `prediction.yaml` + `analysis_plan.yaml`; no `SEALED_STAGES` file exists and no stage path is a key, so every branch `continue`s and all three lists stay empty | `ok = not(…)` → **True on an empty check** |
| 9 `runs_s1s8/` (renamed at the 2026-08-28 merge) | the seal records `runs/<id>/…`; the current path is `runs_s1s8/<id>/…`, not a key | **`unsealed`** — "never sealed", on intact content |

★ The second row is the sharper one. **"We never sealed it" and "we sealed it
and then renamed the directory" became the same output** — the two states this
whole mechanism exists to tell apart.

## It could not see tampering

Because the branch taken was `unsealed` and not `changed`, the bytes were never
read. Appending a falsified prediction to one of those 18 documents gave a
**bit-identical** verdict:

```
pristine: ok=False unsealed=['03_spec.yaml'] changed=[]
tampered: ok=False unsealed=['03_spec.yaml'] changed=[]
VERDICTS IDENTICAL: True
```

Meanwhile CI's bash check — which does have a rename fallback — reported
`42 sealed documents checked, 18 with a stale recorded path` and was green. Two
checkers, opposite verdicts on the same 18 documents, both shipping.

## ★ Why nothing caught it

Every unit test builds its run in `tmp_path` and seals it in place — the one
configuration where the stage path and the seal key coincide. The archive is the
only place they diverge, and the single test that touched the archive
**never ran**:

```python
p = io.REPO_ROOT / "runs" / "2026-07-28_trap-2d-5um_2dfb9d"
if not p.exists():
    pytest.skip("runs/ 는 gitignore 대상 — 이 체크아웃에 없다")
```

That path is absent from **every commit in the history** — the run has been under
`runs_s1s8/` since the initial commit. And the skip reason is false twice over:
`git check-ignore runs` exits 1, and 1311 files under `runs/` are tracked.
Pointed at the real path, the assertion failed. A silent skip, on a false
premise, guarding the seal.

[[provenance-must-have-one-definition-and-three-capture-points]] is the same
family: a mechanism that looks present at every capture point and is checked at
none.

## Two implementations, and the broken one was wired

`bdbot.runcard.verify_seal` was correct all along — it iterates the seal's
entries, falls back to the sibling on a rename, and separates `[warn]` from hard
problems. It is called from `run.execute`, which is why the S30 campaign's seals
really were enforced. The broken one is what `cli.py:365`,
`simbot/report.py:50,222` and `simbot/validate.py:355` use.

The two are **deliberately** not unified, and the reason is recorded in
`bdbot/runcard.py`'s module docstring: `bdbot` cannot import `simbot` (that
direction is a cycle). Fine — but then divergence has to be a test failure, so
there is now one asserting the two agree on every archived seal.

## The recorded path is provenance, not a resolution strategy

CI preferred the recorded path and fell back to the sibling. That ordering is
unsafe in a way the rename concealed: a **copied** run directory carries the
original's paths, so its seal would verify the *original's* documents and pass
while its own went unchecked.

A seal can only ever cover documents in its own directory — `write_seal` cannot
reach outside it. So the rule is now: **verify the sibling; report a recorded
path that no longer matches as drift, and let content alone gate.** Both
checkers follow it, and both now report 42 documents and 18 drifts.

## The advertised external check did not work either

Reports carry, under *"verified without this code"*:

```bash
shasum -a 256 -c SEALED.sha256
```

This resolves in **no** working directory. From the run directory the
repo-relative entries become `runs/<id>/runs/<id>/…`; from the repo root the
seal file is not there. Two committed reports print it directly beside
*"✅ 봉인 검증 통과"*. The seal's only independent verification path was a
command that fails.

It is generated by `report.seal_check_command` now, and **two tests execute its
output** over all 21 seals: one asserting every command passes, one asserting it
prints `FAILED` on an edited document — because a command that only ever passes
is indistinguishable from `true`.

## After

```
21 / 21   verify, with content hashed
42        documents hashed  (= what CI's independent bash check counts)
18        drifts reported, on exactly the 9 renamed directories
21 / 21   the two implementations agree
5 / 5     gate mutants CAUGHT   (verify_gates_bite.py, SEAL-*)
1434 passed, 6 skipped          (was 1424 / 7 — the false skip is gone)
```

`ok` now requires a non-empty `verified` list, and the summary states the count.

## Prevention

- **A pass must say how much it checked.** `SealVerdict.verified` exists because
  `len(entries)` was being reported as a verification count
- **Iterate the artefact, not the vocabulary.** The seal defines its own scope
- **Never name an archive path in a test.** `_archived_seal_dirs()` globs, and a
  separate test fails if it finds fewer than 21 — the same guard CI has against
  `checked -eq 0`
- **Execute an advertised command; never assert its text.** The previous test
  asserted the exact broken string, and so pinned the defect in place
- If two verifiers may not be unified, **assert that they agree**

## See also

[[a-content-hash-is-only-as-portable-as-its-least-portable-operation]] ·
[[a-crystal-start-does-not-survive-at-the-zahn-window-centre]] ·
[[provenance-must-have-one-definition-and-three-capture-points]] ·
`docs/05-pitfalls.md#a-verifier-that-iterates-its-own-vocabulary-instead-of-the-artefact`
