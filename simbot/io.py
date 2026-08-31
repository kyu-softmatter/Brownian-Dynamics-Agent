"""Run directories · hashes · sealing. 0 lines of LLM.

**Sealing is why this module exists.** Without pinning the prediction document by
hash before the run, there is no structural way to prevent post-hoc
rationalisation -- edit the prediction after seeing the result and no record is
left. S7 compares against the pre-run hash and **stops** on a mismatch.

`SEALED.sha256` is in the standard `sha256sum` format
(`<hash>  <repo-relative path>`), so `shasum -a 256 -c SEALED.sha256` verifies it
without this code -- the seal's trustworthiness must not depend on our own code.
"""
from __future__ import annotations

import hashlib
import json
import platform
import re
import subprocess
from dataclasses import dataclass, field
from datetime import date as _date
from pathlib import Path

REPO_ROOT = Path(__file__).resolve().parent.parent

# run directory layout — master_plan.md §1.3
RUN_LAYOUT: dict[str, str] = {
    "input": "00_input",
    "intake": "01_intake.md",
    "intake_json": "01_intake.json",
    "prediction": "02_prediction.md",
    "prediction_json": "02_prediction.json",
    "spec": "03_spec.yaml",
    "spec_rationale": "03_spec_rationale.md",
    "reduced": "04_reduced.yaml",
    "nondim": "04_nondim.md",
    "manifest": "05_run_manifest.json",
    "figures": "06_figures.md",
    "validation": "07_validation.md",
    "sensitivity": "07b_sensitivity.md",
    "conclusion": "08_conclusion.md",
    "metrics": "metrics.json",
    "report": "REPORT.md",
    "seal": "SEALED.sha256",
}

# what gets sealed — documents that must exist before S5 and never change after
SEALED_STAGES: tuple[str, ...] = ("prediction", "intake", "spec")


# =============================================================================
# hashes
# =============================================================================
def sha256_bytes(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def sha256_text(text: str) -> str:
    return sha256_bytes(text.encode("utf-8"))


def sha256_file(path: str | Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_payload(obj) -> str:
    """Hash a dict/list by sorted serialization. For structure hashes such as spec_hash.

    `sort_keys=True` is mandatory -- if a mere change of key order makes it look
    like a different system, the cache is invalidated and "the same run" becomes
    undecidable.

    ⚠ **NOT interchangeable with `bdbot.runid.spec_hash`,** which does the same job
      for the `runs/` content-addressing path. The single difference is
      `ensure_ascii`: this one passes `False` (raw UTF-8), `spec_hash` takes json's
      default `True` (`\\uXXXX` escapes). Measured 2026-08-29: identical on every
      ASCII payload, **different on any payload with a non-ASCII key or value.**
      Deliberately not unified -- 263 run directories are named by the other one.
      `tests/test_cross_package_equivalence.py` pins it.
    """
    return sha256_text(json.dumps(obj, sort_keys=True, ensure_ascii=False,
                                  default=str))


def code_hash(root: Path | None = None) -> str:
    """Hash of the whole `simbot` source (12 chars).

    If the code changes, comparing results is meaningless. `analysis/` is included
    -- a change to the analysis code changes the measurement, so excluding it
    would leave a hole in the seal.
    """
    base = Path(root) if root else Path(__file__).parent
    h = hashlib.sha256()
    for p in sorted(base.rglob("*.py")):
        if "__pycache__" in p.parts:
            continue
        h.update(p.relative_to(base).as_posix().encode())
        h.update(p.read_bytes())
    return h.hexdigest()[:12]


def file_hash(path: str | Path) -> str:
    """Hash of a single file (12 chars). **Used to pin down the driver script.**

    ★ `code_hash` covers only `simbot/`. But what sets a run's parameters (the `A`
      list, seeds, run length, geometry) is the driver in `scripts/` -- and if that
      is not in the hash, **the artefacts alone cannot say what made this run.**
      It was an actual hole in `soft-r3-time-resolved`, 2026-07-29.
    """
    return sha256_file(path)[:12]




def git_rev(cwd: Path | None = None) -> str:
    """The current commit (short form). `"?"` if git is absent or outside a repo."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5,
                             cwd=cwd or REPO_ROOT)
        return out.stdout.strip() or "?"
    except Exception:
        return "?"


def git_dirty(cwd: Path | None = None) -> bool | None:
    """Are there uncommitted changes to tracked files? `None` if undecidable.

    A run made in a dirty tree does not reproduce from `git_rev` alone -- it is
    recorded in the manifest so that "can this run be reproduced" can be answered
    honestly later.
    """
    try:
        out = subprocess.run(["git", "status", "--porcelain"],
                             capture_output=True, text=True, timeout=5,
                             cwd=cwd or REPO_ROOT)
        if out.returncode != 0:
            return None
        return bool(out.stdout.strip())
    except Exception:
        return None


# Only packages that actually affect reproducibility. Include everything and the
# hash changes on harmless updates too.
ENV_PACKAGES: tuple[str, ...] = ("hoomd", "numpy", "scipy", "gsd", "freud")


def env_versions() -> dict[str, str]:
    """Versions of the packages that affect the numbers. A failed import is
    `"absent"`."""
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ENV_PACKAGES:
        try:
            mod = __import__(name)
        except Exception:
            out[name] = "absent"
            continue
        v = getattr(mod, "__version__", None)
        if v is None:                                   # hoomd uses version.version
            v = getattr(getattr(mod, "version", None), "version", "unknown")
        out[name] = str(v)
    return out


def env_hash() -> str:
    return sha256_payload(env_versions())[:12]


def provenance(driver: str | Path | None = None) -> dict:
    """"What made this artefact" — **the single definition of the provenance block.**

    `build_manifest`, both runners and the analysis scripts all use this. Built by
    hand at each call site the key names drift apart, and then the names
    `report.reproducibility_section` reads stop matching the ones the runner writes
    -- so **the reproducibility information renders as a silent blank.**
    That actually happened on 2026-07-29: `report.py` reads `env_hash`, and both
    runners built their own manifest by hand without that key.

    **It is called at analysis time too.** A trajectory's manifest holds the hashes
    from *trajectory-generation* time, and the analysis can be run separately later
    (`--analyze-only`) -- and then there is no telling which analysis code and which
    `freud` produced `metrics.json`.

    Args:
        driver: the script that defined the run. `code_hash` covers only `simbot/`,
            so the driver in `scripts/` -- which sets the `A` list, seeds, run
            length, analysis windows and so on -- is pinned by this argument alone.
            ⚠ **It is the hash of a single file.** If the driver imports another
            module from `scripts/`, that one is not covered -- the current drivers
            import only `simbot`, so `code_hash` + `driver_hash` together cover
            everything.
    """
    out = {
        "code_hash": code_hash(),
        "git_rev": git_rev(),
        "git_dirty": git_dirty(),
        "env_hash": env_hash(),
        "env": env_versions(),
    }
    if driver is not None:
        #  ★ It takes several files. If the driver imports another module from
        #    `scripts/`, that one also sets run parameters and must be hashed too --
        #    catch only one and the claim "code_hash + driver_hash cover everything"
        #    becomes false.
        paths = ([driver] if isinstance(driver, (str, Path))
                 else list(driver))
        pairs = {}
        for d in paths:
            p = Path(d).resolve()
            pairs[_seal_relpath(p)] = file_hash(p) if p.exists() else "?"
        if len(pairs) == 1:
            (rel, h), = pairs.items()
            out["driver"] = rel
            out["driver_hash"] = h
        else:
            out["driver"] = sorted(pairs)
            out["drivers"] = pairs
            #  a composite hash — it changes if any single file changes
            out["driver_hash"] = sha256_payload(pairs)[:12]
    return out


# =============================================================================
# run_id
# =============================================================================
_SLUG_OK = re.compile(r"[^a-z0-9]+")


def slugify(text: str) -> str:
    return _SLUG_OK.sub("-", text.strip().lower()).strip("-") or "run"


def new_run_id(slug: str, spec_hash: str, when: _date | None = None) -> str:
    """`run_id = <date>_<slug>_<first 6 of the spec hash>`.

    Deterministic when `when` is given (for tests and reproduction). Today if
    omitted.
    """
    day = (when or _date.today()).isoformat()
    return f"{day}_{slugify(slug)}_{spec_hash[:6]}"


# =============================================================================
# run directory
# =============================================================================
@dataclass(frozen=True)
class RunDir:
    """A path calculator for a run directory. It creates no files, only knows paths.

    Access is by `RUN_LAYOUT` key -- scatter the filename strings through the code
    and the name S6 writes quietly diverges from the name S7 reads.
    """

    path: Path

    @classmethod
    def create(cls, runs_root: str | Path, run_id: str) -> RunDir:
        p = Path(runs_root) / run_id
        (p / RUN_LAYOUT["input"]).mkdir(parents=True, exist_ok=True)
        (p / "figs").mkdir(exist_ok=True)
        (p / "raw").mkdir(exist_ok=True)
        return cls(p)

    @property
    def run_id(self) -> str:
        return self.path.name

    @property
    def figs(self) -> Path:
        return self.path / "figs"

    @property
    def raw(self) -> Path:
        return self.path / "raw"

    def file(self, stage: str) -> Path:
        if stage not in RUN_LAYOUT:
            raise KeyError(f"unknown stage {stage!r}; known: {sorted(RUN_LAYOUT)}")
        return self.path / RUN_LAYOUT[stage]

    def exists(self, stage: str) -> bool:
        return self.file(stage).exists()

    def write(self, stage: str, text: str) -> Path:
        p = self.file(stage)
        p.parent.mkdir(parents=True, exist_ok=True)
        p.write_text(text, encoding="utf-8")
        return p

    def write_json(self, stage: str, obj) -> Path:
        return self.write(stage, json.dumps(obj, indent=2, ensure_ascii=False,
                                            default=str) + "\n")

    def read(self, stage: str) -> str:
        return self.file(stage).read_text(encoding="utf-8")

    def read_json(self, stage: str):
        return json.loads(self.read(stage))

    def completed_stages(self) -> list[str]:
        """Stages that already have artefacts — the basis for `resume` skipping
        recomputation."""
        return [s for s in RUN_LAYOUT if s != "input" and self.exists(s)]


# =============================================================================
# sealing
# =============================================================================
@dataclass(frozen=True)
class SealEntry:
    digest: str
    relpath: str


@dataclass
class SealVerdict:
    """Seal-verification result. **On `ok=False`, S7 must stop.**"""

    ok: bool
    changed: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    unsealed: list[str] = field(default_factory=list)   # should be sealed, not listed
    entries: dict[str, str] = field(default_factory=dict)

    def summary(self) -> str:
        if self.ok:
            n = len(self.entries)
            return f"seal verified — nothing changed after the run ({n} sealed)"
        parts = []
        if self.changed:
            parts.append(f"changed {self.changed}")
        if self.missing:
            parts.append(f"vanished {self.missing}")
        if self.unsealed:
            parts.append(f"not sealed {self.unsealed}")
        return "seal violation — " + " · ".join(parts)


def _seal_relpath(path: Path) -> str:
    """Repo-relative path. Absolute if outside the repo (a test tmpdir, say)."""
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def write_seal(rundir: RunDir, stages: tuple[str, ...] = SEALED_STAGES) -> Path:
    """Write the hashes of the existing sealable documents into `SEALED.sha256`.

    ⚠ **Must be called before S5 runs.** Seal after the run and the seal guarantees
      nothing -- this function cannot detect that, so the pipeline has to keep the
      order.
    """
    lines = []
    for stage in stages:
        p = rundir.file(stage)
        if p.exists():
            lines.append(f"{sha256_file(p)}  {_seal_relpath(p)}")
    if not lines:
        raise FileNotFoundError(
            f"nothing to seal — not one of {list(stages)} exists. "
            f"Write the S2 prediction first.")
    out = rundir.file("seal")
    out.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return out


def read_seal(rundir: RunDir) -> dict[str, str]:
    """`SEALED.sha256` → `{relpath: digest}`."""
    text = rundir.file("seal").read_text(encoding="utf-8")
    out = {}
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        digest, _, rel = line.partition("  ")
        if not rel:                                    # a single-space separator too
            digest, _, rel = line.partition(" ")
        out[rel.strip()] = digest.strip()
    return out


def verify_seal(rundir: RunDir, stages: tuple[str, ...] = SEALED_STAGES) -> SealVerdict:
    """Check that the sealed documents did not change after the run.

    A document in `stages` that is not in the seal list at all is reported as
    `unsealed` -- needed so that "the seal file passed" does not come to mean "the
    prediction was sealed".
    """
    if not rundir.exists("seal"):
        return SealVerdict(ok=False, missing=[RUN_LAYOUT["seal"]])

    sealed = read_seal(rundir)
    changed, missing, unsealed = [], [], []
    for stage in stages:
        p = rundir.file(stage)
        rel = _seal_relpath(p)
        if rel not in sealed:
            if p.exists():
                unsealed.append(RUN_LAYOUT[stage])
            continue
        if not p.exists():
            missing.append(RUN_LAYOUT[stage])
            continue
        if sha256_file(p) != sealed[rel]:
            changed.append(RUN_LAYOUT[stage])

    return SealVerdict(ok=not (changed or missing or unsealed),
                       changed=changed, missing=missing, unsealed=unsealed,
                       entries=sealed)


# =============================================================================
# manifest
# =============================================================================
def build_manifest(*, run_id: str, spec_hash: str, seed, extra: dict | None = None,
                   rundir: RunDir | None = None) -> dict:
    """The contents of `05_run_manifest.json`. Everything reproduction needs.

    Pinning the prediction hash here keeps the seal alive even if the
    `SEALED.sha256` file is deleted (deleting the manifest too loses it, but that
    creates the friction of having to edit two places at once).
    """
    man = {
        "run_id": run_id,
        "spec_hash": spec_hash,
        "seed": seed,
        #  ★ provenance is built in exactly one place — list it by hand here and
        #    the key names drift from the runners and the analysis (see the
        #    `provenance()` docstring)
        **provenance(),
    }
    if rundir is not None and rundir.exists("seal"):
        man["sealed"] = read_seal(rundir)
    man.update(extra or {})
    return man
