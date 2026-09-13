"""Content-addressed run_id plus run-directory preparation.

Same spec -> same `run_id` -> **do not re-run.** Both of the first two cases had
adopted this convention independently.

Why `nhex` is an argument: one case used 12 hex characters and the other 10. The
existing run directories and their `run_id`s have to stay valid (reproducibility),
so rather than unify it, the case declares it.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

#: Files `prepare_outdir` must NOT delete. The distinction is input vs output:
#: `prepare_outdir` exists to clear the *outputs* of a partial run, and every
#: name here is a human-authored *input* that is written before the run and
#: that a gate downstream has to be able to see.
#:
#: ★ Measured 2026-09-13, and the reason this set grew. `execute()` checks the
#:   seal at run.py:416 and `params.json` at run.py:424 -- both AFTER
#:   `prepare_outdir` at run.py:399 had already deleted them. Verified by
#:   running `execute()` on a temp directory:
#:
#:     a VALID seal          -> "no SEALED.sha256. Refusing to run"
#:     a TAMPERED seal       -> "[warn] ... is unsealed", and the run PROCEEDS
#:     an APPROVED params    -> "no params.json. Rule 10 requires every number"
#:     a DRIFTING params     -> the drift check never runs at all
#:
#:   So `require_seal=True` and `require_approval=True` could only ever fail,
#:   never be satisfied; the two branches behind them were unreachable; and
#:   `runcard.verify_or_raise`'s documented asymmetry -- *"unsealed is a known
#:   historical state, broken is tampering"* -- was inverted, with tampering
#:   silently downgraded to "unsealed". It also explains why 0 of 260 run
#:   directories carry a seal: sealing before a run was undone by the run.
#:
#:   The structure gate does not appear in that list only because it was hoisted
#:   above `prepare_outdir` in a previous pass (run.py:324 says so). The two
#:   gates below it had the identical defect and were walked past.
SEAL_FILE = "SEALED.sha256"    # = runcard.SEAL_NAME; runcard imports us, not vice versa
PRESERVE = {
    "record.json",             # a lesson must outlive the run artefacts
    "params.json",             # rule 10: the approved numbers
    SEAL_FILE,                 # the pre-registration
}

# * Keys excluded from the run_id hash -- documentation, provenance, derived values.
#   This actually bit: adding a `derived_from` field to a system file changed one
#   case's run_id from 70b9394e7310 to dc67e4e2b825, because the whole YAML had
#   been put into the spec and hashed.
#   **If editing one comment line invalidates a run, content addressing is
#   useless.** A run_id must track only "what was simulated".
#   (The other case had listed the physics fields explicitly and did not have this
#   problem -- that was the correct design.)
#
#   WARNING: the reverse direction is equally dangerous and also happened. One
#   spec contained no physical system at all, so changing d from 5um to 0.5um and
#   eta by 62x -- a 16.1x change in tau_B -- left the run_id *identical*, and an
#   old result got reported as the new system's result. The hash must cover
#   everything that fixes the physics and nothing that documents it. Both
#   directions are guarded by verify/verify_l3_spec_gaps.py.
DOC_KEYS = frozenset({
    "description", "derived_from", "not_verified", "required_convergence_checks",
    "derived_scales", "dimensionless", "source", "note", "source_note",
    "interpretation", "deviates_from_sketch", "expr", "role", "meaning", "tier",
    "what", "proposed", "followup", "lean", "confirmed_by",
    # `structure` (physical.STRUCTURE_SECTION) documents WHY the dimensionality
    # is what it is. The choice itself stays hashed where it already was --
    # top-level `dimensions` -- so 2 -> 3 still re-ids the run while writing
    # down why it is 2 does not. Both directions of the warning above are
    # therefore preserved, and `physical.check_dim` cross-checks the two so they
    # cannot drift. Existing specs have no `structure` key, so adding this
    # renames nothing: verified against all 278 specs in
    # `verify/verify_dim_gate.py` section (1).
    "structure",
})


def physics_only(node):
    """Recursively strip documentation, provenance and derived fields.

    What remains is only what fixes the physics.
    """
    if isinstance(node, dict):
        return {k: physics_only(v) for k, v in node.items() if k not in DOC_KEYS}
    if isinstance(node, list):
        return [physics_only(v) for v in node]
    return node


def spec_hash(spec: dict, nhex: int = 12) -> str:
    """sha256 of the spec's sorted JSON -> the first nhex characters.

    Deterministic regardless of key order or whitespace.

    ⚠ **NOT interchangeable with `simbot.io.sha256_payload`,** which does the same
      job for the S2 sealing path. Measured 2026-08-29:

          payload            spec_hash(12)   sha256_payload()[:12]   same
          ASCII only         93a7e2f22fd8    93a7e2f22fd8            yes
          non-ASCII value    5ef45ee9ed7c    062a55bdd263            NO
          non-ASCII key      6b178a2bd90d    df013e47cf6b            NO

      The single difference is `ensure_ascii`: this function takes json's default
      (`True`, so non-ASCII becomes `\\uXXXX`), `sha256_payload` passes `False`
      (raw UTF-8). They agree on every ASCII payload, which is why the divergence
      never surfaced -- and this repository is full of Korean strings, so it was
      one non-stripped field away from surfacing.
      **Deliberately not unified**: 263 `runs/` directories and 279 `specs/` files
      are named by this function's output, so changing the serialization renames
      all of them. `tests/test_cross_package_equivalence.py` pins the difference so
      it cannot be "fixed" by accident.
      What makes the ASCII agreement hold at all is `physics_only`, which strips
      `source`, `note`, `description` and the rest of `DOC_KEYS` -- i.e. exactly
      the fields that carry prose -- before this ever sees the spec.
    """
    blob = json.dumps(spec, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:nhex]


def content_run_id(label: str, spec: dict, tag: str | None = None, nhex: int = 12) -> str:
    mid = f"{tag}__" if tag else ""
    return f"{label}__{mid}{spec_hash(spec, nhex)}"


def prepare_outdir(outdir: Path, force: bool = False) -> tuple[bool, str]:
    """(should we run, message). False when a completed run already exists.

    The presence of `result.txt` counts as complete (a directory holding only
    partial artefacts is cleared and rewritten).

    WARNING: `result.txt` is written by the *case script*, not by `bdbot.run`. A
    case that never added that line reports zero runs here and in
    `bdbot.cli status` while its `metrics.json` files sit on disk -- and the
    convention once caused a "clean up incomplete runs" pass to delete 6
    completed runs.
    """
    if (outdir / "result.txt").exists() and not force:
        prev = (outdir / "result.txt").read_text()
        # Accept both markers: 167 result.txt files in the archive were written
        # before the case scripts were translated (2026-08-29), and runs/ is not
        # rewritten because it is the content-addressed evidence ledger.
        marker = next((m for m in ("result —", "결과 —") if m in prev), None)
        tail = prev.split(marker)[-1] if marker else prev[-1200:]
        return False, (f"\nthis run is already complete: runs/{outdir.name}/  "
                       f"(--force to re-run)\n{tail}")
    if outdir.exists():
        # * Keep PRESERVE. A lesson (a KB entry) must outlive the run artefacts;
        #   run_id is content-addressed, so the same directory means the same
        #   spec and the previous lesson is still valid. (Added after a --force
        #   re-run destroyed 6 lessons.) The approved numbers and the seal are
        #   kept for the reason recorded on PRESERVE above.
        # * ...and whatever the seal COVERS. Keeping `SEALED.sha256` while
        #   deleting the document it hashes leaves a dangling seal, which
        #   `verify_seal` then reports as a missing document -- trading a silent
        #   failure for a loud one that blocks every re-run of a sealed case.
        keep = set(PRESERVE)
        seal = outdir / SEAL_FILE
        if seal.exists():
            for line in seal.read_text().splitlines():
                part = line.split(None, 1)
                if len(part) == 2:
                    keep.add(Path(part[1].strip()).name)
        for f in outdir.iterdir():
            if f.is_file() and f.name not in keep:
                f.unlink()
    outdir.mkdir(parents=True, exist_ok=True)
    return True, ""


def write_spec(outdir: Path, spec: dict) -> None:
    (outdir / "spec.json").write_text(json.dumps(spec, indent=2, default=str))


def list_artifacts(outdir: Path, root: Path) -> list[str]:
    lines = [f"\nartefacts: {outdir.relative_to(root)}/"]
    for f in sorted(outdir.iterdir()):
        lines.append(f"   {f.name:<22} {f.stat().st_size / 1024:8.1f} KB")
    return lines


__all__ = ["spec_hash", "content_run_id", "prepare_outdir", "write_spec", "list_artifacts"]
