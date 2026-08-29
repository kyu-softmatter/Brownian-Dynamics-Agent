# NOTICE — what is not published here

This repository is public. It was assembled on 2026-08-28 by merging three
private predecessor repositories (see [`docs/00-merge-decisions.md`](docs/00-merge-decisions.md)),
and four classes of file were held back at assembly time rather than removed
later. **They are in no commit of this repository** — the history starts from
the merged tree, so there is nothing to scrub.

Nothing below is load-bearing for a conclusion. Every number this project
relies on is either recomputed by code in this repository or carried inline
with its citation in a distillation under `knowledge/source/`.

---

## 1 · Copyrighted papers and book chapters — obtain from the publisher

This repository has no licence to redistribute them.

| Held back | What it was | Where the content lives now |
|---|---|---|
| `kb/cambridge-core_…/*.pdf` (21 chapter PDFs) | Leal, *Microstructural Rheology of Complex Fluids*, Cambridge 2026 | [`knowledge/source/books/leal_microstructural_rheology.md`](knowledge/source/books/leal_microstructural_rheology.md) — our distillation |
| `kb/file_…pdf` | Welty et al., *Momentum, Heat and Mass Transfer*, 5th ed. | [`knowledge/source/books/welty_transport.md`](knowledge/source/books/welty_transport.md) — our distillation |
| `intake/chain-bend-*/PhysRevLett.94.138301.pdf` | Pantina & Furst, *PRL* **94**, 138301 (2005) | [`knowledge/source/papers/2005-pantina-furst-bending-coefficient.md`](knowledge/source/papers/2005-pantina-furst-bending-coefficient.md) |
| `intake/chain-bend-*/la7023617.pdf` | Pantina & Furst, *Langmuir* **24**, 1141 (2008) | distilled in the same file |
| `knowledge/raw/lab/**` (32 PDFs + SI) | Group and reference papers collected for literature scan | 42 distillations in [`knowledge/source/papers/`](knowledge/source/papers/) |

The distillations are our own summaries: what the paper gives us in one line,
its equations **converted into this project's non-dimensional conventions**, the
values we reproduced, and what we could not reproduce. Thirty-five of the 42
carry a DOI in their frontmatter, so a reader can fetch the original; seven
carry no DOI field and have to be found by title. Fifty-six claims taken out
of the two books are re-derived numerically by
[`verify/verify_book_claims.py`](verify/verify_book_claims.py) (56/56 pass), and
that script labels each check `[BOOK]` (the book's own reported number, i.e. did
we read it right), `[DERIV]` (the book's formulas checked against each other),
or `[OURS]` (checked against our measurements). None of those checks needs the
PDF.

**Bibliographic data in these distillations is largely unchecked.** One
carries an explicit `verified: false` — meaning its citation was written from
memory — and the other 41 carry no `verified` field at all, which operationally
means the same thing. The physics in a distillation was checked (that is what
the `[BOOK]`/`[DERIV]`/`[OURS]` labels are for); the volume-and-page line
mostly was not. Confirm the citation before it goes into a manuscript.

Thirty-eight of the 42 are marked `lab_authored: true` — the literature base is
mostly the group's own published work, which is why it grew fast and why it is
narrow. That is a real limitation on any claim of the form *"the literature
says…"* in this repository.

## 2 · Trajectories and raw arrays — regenerate, do not download

Measured on the 254 completed runs carried here: the text artifacts are
**5.7 MB**, the `.gsd` trajectories and `.npz` observable arrays are **542 MB** —
95× larger, and reproducible from `spec + seed` because a run is content-
addressed (see [`bdbot/runid.py`](bdbot/runid.py)).

So `runs/` here holds `metrics.json`, `record.json`, `result.txt`, `l4.json`
and `spec.json` — the numbers, the verdicts and the lessons — and not the
trajectory that produced them. `runs_s1s8/` additionally keeps `SEALED.sha256`,
because that file is the evidence that a prediction was fixed before the run
and `shasum -a 256 -c SEALED.sha256` has to be runnable from a fresh clone.

⚠️ **A retained `spec.json` does not by itself prove a run is reproducible.**
`run_id` hashes the physics fields, so when the code evolves the same case
re-runs under a *different* id. Measured 2026-08-05: of 79 legacy runs missing
step-resolution telemetry, **zero** reproduced under the current code. Stale
specs stay on disk. Treat `specs/` as a record of what was asked for, not a
guarantee of what re-runs today.

## 3 · Unpublished lab assets — never present

`knowledge/source/lab/` was reserved in the predecessor design for
senior-students' code, notes and unpublished parameters, and was gitignored
pre-emptively so it could never need scrubbing. **It was still empty when this
repository was assembled**, and the `.gitignore` rule is carried forward.

The dividing line is publication status, not authorship: a group-authored paper
that already has a DOI is distilled like any other and marked
`lab_authored: true`, rather than being hidden alongside unpublished material.
Hiding by authorship would have dropped already-public work from the public
repository.

## 4 · Hand sketches — published, downscaled

The six input sketches **are** here, at reduced resolution. Originals were
8.6–10 MB JPEGs (~55 MB total, phone camera); they are carried at a 2000 px
long edge, ~0.6–1.1 MB each, which stays legible for every annotation the
intake reading cites. `intake/network/transcription_annotated.png` (2.9 MB)
became a JPEG for the same reason.

If a reading is ever disputed on the grounds of image resolution, the originals
are on the author's machine — the sketches are the one input class that cannot
be regenerated, so they are backed up outside this repository rather than
trusted to it.

---

## Licence

**No licence is granted.** This repository is published for reading, and all
rights are reserved. If you want to use any part of it, ask.
