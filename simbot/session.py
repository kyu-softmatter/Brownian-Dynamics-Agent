"""Conversational session state. 0 lines of LLM.

> **The state lives in a file, not in the conversation.** A conversation disappears
> and cannot be inspected.

`sessions/<id>/session.yaml` holds the per-turn history (`runs/<id>/` is per run).
It has to be possible to pick up with `show` after the conversational context is
gone.

## `set` does not run anything

`set` changes a parameter and **only estimates the cost**. Running is `run`'s job.
Without that separation, "let's try N = 8000" quietly starts an 11-minute run.

## The history is append-only

A turn is never edited; a new turn is stacked on top. Overwriting a past turn loses
"what was tried and what did not work" -- which is half of this project's output.
"""
from __future__ import annotations

import json
from dataclasses import asdict, dataclass, field, fields, is_dataclass
from datetime import datetime
from pathlib import Path

import yaml

from .estimators import estimate_wall_time_s
from .io import REPO_ROOT, slugify
from .policy import Policy, load_policy
from .spec import Quantity, SystemSpec, dump_yaml, validate as validate_spec

SESSIONS_ROOT = REPO_ROOT / "sessions"


# =============================================================================
# dotted-key assignment
# =============================================================================
def _navigate(root, path: list[str]):
    """Follows `['species', '0', 'radius_si']` and returns (parent, last key)."""
    cur = root
    for part in path[:-1]:
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict):
            cur = cur[part]
        elif is_dataclass(cur) and part in {f.name for f in fields(cur)}:
            cur = getattr(cur, part)
        else:
            raise KeyError(f"path {'.'.join(path)}: {part!r} not found")
        if cur is None:
            raise KeyError(f"path {'.'.join(path)}: {part!r} is None")
    return cur, path[-1]


def _read(parent, key: str):
    if isinstance(parent, list):
        return parent[int(key)]
    if isinstance(parent, dict):
        return parent[key]
    return getattr(parent, key)


def _parse_scalar(text: str):
    """string → int/float/bool/str/list.

    ★ **`yaml.safe_load` must not be used as is.** YAML 1.1's float regex requires a
      decimal point in the mantissa, so `5e-3` and `2e-5` read as **strings**
      (`2.5e-3` is a float). Then a string lands in `dt_star` and the
      non-dimensionalization dies with a TypeError -- or worse, a comparison passes
      quietly. The same trap was hit with `6.3e6` in `config/run_policy.yaml` --
      see `simbot/policy.py::_find_numeric_strings`.

    So **a number is tried first**, and only otherwise is it handed to YAML.
    """
    s = text.strip()
    try:
        return int(s)
    except ValueError:
        pass
    try:
        return float(s)
    except ValueError:
        pass
    return yaml.safe_load(s)


def set_by_path(spec: SystemSpec, dotted: str, raw: str, *, turn: int) -> tuple:
    """Applies `numerics.dt_star=2.5e-3` and returns `(old value, new value)`.

    For a `Quantity` field, only `value` changes and **the provenance becomes
    `user`** -- a value a human set directly in this run is neither an assumption,
    nor an inference, nor policy. That distinction decides S7b's sensitivity target
    list (a human-chosen value is not something to shake).
    """
    parent, key = _navigate(spec, dotted.split("."))
    old = _read(parent, key)
    new_value = _parse_scalar(raw)

    if isinstance(old, Quantity):
        prev = old.value
        old.value = new_value
        old.provenance = "user"
        old.basis = (f"set by a human on turn {turn} (previously {prev!r}"
                     + (f", original basis: {old.basis}" if old.basis else "") + ")")
        old.confidence = old.confidence or "high"
        return prev, new_value

    if isinstance(parent, list):
        parent[int(key)] = new_value
    elif isinstance(parent, dict):
        parent[key] = new_value
    else:
        setattr(parent, key, new_value)
    return old, new_value


# =============================================================================
# Turn
# =============================================================================
@dataclass
class Turn:
    index: int
    timestamp: str
    kind: str                                  # new | set | run
    changes: dict[str, list] = field(default_factory=dict)
    spec_hash: str = ""
    cost: dict = field(default_factory=dict)
    run_id: str = ""
    metrics: dict = field(default_factory=dict)
    note: str = ""
    problems: list[str] = field(default_factory=list)


# =============================================================================
# Session
# =============================================================================
@dataclass
class Session:
    session_id: str
    spec_path: str
    spec: SystemSpec
    turns: list[Turn] = field(default_factory=list)
    root: Path = SESSIONS_ROOT

    # --- paths ---
    @property
    def dir(self) -> Path:
        return Path(self.root) / self.session_id

    @property
    def state_file(self) -> Path:
        return self.dir / "session.yaml"

    def derived_spec_file(self, turn: int) -> Path:
        return self.dir / f"spec_turn{turn:02d}.yaml"

    # --- create / save / load ---
    @classmethod
    def create(cls, spec_path: str | Path, *, root: Path | None = None,
               when: datetime | None = None) -> Session:
        spec = SystemSpec.load(spec_path)
        stamp = (when or datetime.now()).strftime("%Y-%m-%dT%H-%M")
        sid = f"{stamp}_{slugify(spec.card)}"
        s = cls(session_id=sid, spec_path=str(spec_path), spec=spec,
                root=Path(root) if root else SESSIONS_ROOT)
        s.dir.mkdir(parents=True, exist_ok=True)
        rep = validate_spec(spec)
        s._append(Turn(index=0, timestamp=_now(when), kind="new",
                       spec_hash=spec.hash(), note=f"spec: {spec_path}",
                       problems=rep.problems))
        return s

    @classmethod
    def load(cls, session_id: str, *, root: Path | None = None) -> Session:
        base = Path(root) if root else SESSIONS_ROOT
        f = base / session_id / "session.yaml"
        if not f.exists():
            raise FileNotFoundError(f"session {session_id} does not exist ({f})")
        raw = yaml.safe_load(f.read_text(encoding="utf-8"))
        turns = [Turn(**t) for t in raw.get("turns", [])]
        # the latest turn's derived spec is the canonical one — the original spec is
        # never modified
        spec_file = base / session_id / raw["current_spec"]
        return cls(session_id=session_id, spec_path=raw["spec_path"],
                   spec=SystemSpec.load(spec_file), turns=turns, root=base)

    @classmethod
    def latest(cls, *, root: Path | None = None) -> Session:
        base = Path(root) if root else SESSIONS_ROOT
        cands = sorted(p.name for p in base.glob("*") if (p / "session.yaml").exists())
        if not cands:
            raise FileNotFoundError(f"no session under {base}. Run "
                                    f"`session new <spec>` first.")
        return cls.load(cands[-1], root=base)

    def _append(self, turn: Turn) -> Turn:
        self.turns.append(turn)
        self.derived_spec_file(turn.index).write_text(dump_yaml(self.spec),
                                                     encoding="utf-8")
        self.save()
        return turn

    def save(self) -> Path:
        self.dir.mkdir(parents=True, exist_ok=True)
        payload = {
            "session_id": self.session_id,
            "spec_path": self.spec_path,
            "card": self.spec.card,
            "question": self.spec.question,
            "current_spec": self.derived_spec_file(self.turns[-1].index).name,
            "turns": [asdict(t) for t in self.turns],
        }
        self.state_file.write_text(
            yaml.dump(payload, sort_keys=False, allow_unicode=True, width=100),
            encoding="utf-8")
        return self.state_file

    # --- commands ---
    def set(self, assignments: list[str], *, policy: Policy | None = None,
            when: datetime | None = None) -> Turn:
        """Change a parameter → a derived spec + **a cost estimate only. It does not
        run.**"""
        changes: dict[str, list] = {}
        for a in assignments:
            if "=" not in a:
                raise ValueError(f"{a!r} is not in `key=value` form")
            key, _, raw = a.partition("=")
            old, new = set_by_path(self.spec, key.strip(), raw.strip(),
                                   turn=len(self.turns))
            changes[key.strip()] = [old, new]

        rep = validate_spec(self.spec)
        return self._append(Turn(
            index=len(self.turns), timestamp=_now(when), kind="set",
            changes=changes, spec_hash=self.spec.hash(),
            cost=self.estimate_cost(policy=policy), problems=rep.problems,
            note="cost estimate only — to run it, `session run`"))

    def estimate_cost(self, *, policy: Policy | None = None) -> dict:
        """How long running this spec now would take. **If it is over budget it says
        so.**"""
        pol = policy or load_policy()
        from .nondim import reduce_spec
        try:
            r = reduce_spec(self.spec)
        except Exception as e:                        # a card with no registered scales
            return {"error": f"non-dimensionalization failed: {e}"}

        n_seeds = int(self.spec.numerics.n_seeds.si) if self.spec.numerics.n_seeds \
            else pol.seeds_default
        k = min(n_seeds, pol.concurrency())
        steps = r.equil_steps + r.prod_steps
        wall = estimate_wall_time_s(r.n_particles, steps, k)
        budget = pol.wall_budget_s
        out = {
            "n_particles": r.n_particles,
            "steps_per_seed": steps,
            "n_seeds": n_seeds,
            "concurrency": k,
            "efficiency": pol.efficiency(k),
            "wall_s_per_seed": round(wall, 3),
            "wall_s_batch": round(wall, 3),           # run concurrently, so batch ≈ 1
            "budget_s": budget,
            "over_budget": wall > budget,
        }
        if out["over_budget"]:
            out["action"] = (f"over the {budget:g} s budget — policy says **report "
                             f"without running** (run_policy §5)")
        if r.n_particles < 500:
            out["warning"] = ("N < 500 — the throughput model was measured at "
                              "N ≥ 500. This estimate underestimates (measure it "
                              "with a pilot)")
        return out

    def record_run(self, *, run_id: str, metrics: dict, wall_s: float | None = None,
                   problems: list[str] | None = None,
                   when: datetime | None = None) -> Turn:
        """Records a run's result in a turn. The running itself is `cli.py`'s job."""
        return self._append(Turn(
            index=len(self.turns), timestamp=_now(when), kind="run",
            spec_hash=self.spec.hash(), run_id=run_id, metrics=metrics,
            cost=({"wall_s_measured": round(wall_s, 3)} if wall_s else {}),
            problems=list(problems or [])))

    def run(self, *, prediction: str | Path | None = None,
            runs_root: str | Path | None = None, extra_args: list[str] | None = None,
            when: datetime | None = None) -> Turn:
        """Runs `cli.py run` on the current turn's derived spec and records the
        result in a turn.

        ★ **Separated from `set`.** `set` only estimates; `run` runs. Without that
          boundary, "let's try N = 8000" quietly starts an 11-minute run.

        On failure the stderr goes into `problems` and **the run turn is recorded
        anyway** -- a failed attempt vanishing from the history loses half of "what
        was tried".
        """
        import subprocess
        import sys

        turn_index = self.turns[-1].index
        spec_file = self.derived_spec_file(turn_index)
        # ★ The session sets the run_id. `cli.py`'s default slug is the spec's parent
        #   directory name, and a session spec lives in `sessions/<date>_<card>/`, so
        #   the date would appear twice.
        run_id = f"{self.session_id}_t{turn_index:02d}_{self.spec.hash()[:6]}"
        cmd = [sys.executable, str(REPO_ROOT / "cli.py"), "run", str(spec_file),
               "--run-id", run_id]
        if prediction:
            cmd += ["--prediction", str(prediction)]
        if runs_root:
            cmd += ["--runs-root", str(runs_root)]
        cmd += list(extra_args or [])

        proc = subprocess.run(cmd, capture_output=True, text=True, cwd=str(REPO_ROOT))
        print(proc.stdout)
        if proc.returncode != 0:
            return self._append(Turn(
                index=len(self.turns), timestamp=_now(when), kind="run",
                spec_hash=self.spec.hash(),
                problems=[f"cli.py run failed (exit {proc.returncode})",
                          *(proc.stderr.strip().splitlines()[-5:] or ["(no stderr)"])],
                note="a failed attempt stays in the history too"))

        metrics, wall = {}, None
        for line in proc.stdout.splitlines():
            if "measured wall" in line:
                try:
                    wall = float(line.split("measured wall")[1].split()[0])
                except (IndexError, ValueError):
                    pass
        root = Path(runs_root) if runs_root else REPO_ROOT / "runs"
        mfile = root / run_id / "metrics.json"
        if mfile.exists():
            raw = json.loads(mfile.read_text(encoding="utf-8"))
            metrics = {k: v["value"] for k, v in raw.items()
                       if isinstance(v, dict) and "value" in v}
        return self.record_run(run_id=run_id, metrics=metrics, wall_s=wall, when=when)

    # --- views ---
    def show(self) -> str:
        """A summary for picking a session back up. Readable with no conversational
        context."""
        out = [f"# session `{self.session_id}`", "",
               f"**card** `{self.spec.card}`",
               f"**question** {self.spec.question.strip()}",
               f"**original spec** `{self.spec_path}`",
               f"**{len(self.turns)} turns** · current spec_hash "
               f"`{self.spec.hash()[:12]}`",
               "", "| turn | time | kind | what | problems |", "|---|---|---|---|---|"]
        for t in self.turns:
            if t.kind == "set":
                what = " · ".join(f"`{k}`: {v[0]!r} → {v[1]!r}"
                                  for k, v in t.changes.items())
                if t.cost.get("wall_s_batch") is not None:
                    what += f" (estimate `{t.cost['wall_s_batch']:g} s`)"
                    if t.cost.get("over_budget"):
                        what += " ⚠️ **over budget**"
            elif t.kind == "run":
                what = f"run `{t.run_id}`"
                if t.metrics:
                    what += " · " + ", ".join(f"{k}=`{_short(v)}`"
                                              for k, v in list(t.metrics.items())[:3])
            else:
                what = t.note
            out.append(f"| {t.index} | {t.timestamp} | {t.kind} | {what} "
                       f"| {len(t.problems) or ''} |")
        probs = [(t.index, p) for t in self.turns for p in t.problems]
        if probs:
            out += ["", "### unresolved problems",
                    *[f"- (turn {i}) {p}" for i, p in probs]]
        out += ["", "next: `session set k.v=…` (estimate only) or `session run` "
                    "(execute)"]
        return "\n".join(out)

    def compare(self, i: int, j: int) -> str:
        """A comparison table of two turns' parameters and measurements."""
        a, b = self._turn(i), self._turn(j)
        rows = [f"# turn {i} vs turn {j}", "",
                f"| | turn {i} | turn {j} |", "|---|---|---|",
                f"| time | {a.timestamp} | {b.timestamp} |",
                f"| kind | {a.kind} | {b.kind} |",
                f"| spec_hash | `{a.spec_hash[:12]}` | `{b.spec_hash[:12]}` |",
                f"| run_id | `{a.run_id or '—'}` | `{b.run_id or '—'}` |"]
        if a.spec_hash == b.spec_hash:
            rows.append("| **difference** | the same spec — the parameters are "
                        "identical | |")

        # parameter differences — accumulate every set between the two turns
        lo, hi = min(i, j), max(i, j)
        acc: dict[str, list] = {}
        for t in self.turns[lo + 1:hi + 1]:
            for k, v in t.changes.items():
                acc[k] = [acc.get(k, v)[0], v[1]]
        if acc:
            rows += ["", "## parameter changes", "", "| key | before | after |",
                     "|---|---|---|"]
            rows += [f"| `{k}` | `{v[0]}` | `{v[1]}` |" for k, v in acc.items()]

        keys = sorted(set(a.metrics) | set(b.metrics))
        if keys:
            rows += ["", "## measurements", "",
                     f"| quantity | turn {i} | turn {j} | change |",
                     "|---|---|---|---|"]
            for k in keys:
                va, vb = a.metrics.get(k), b.metrics.get(k)
                delta = ""
                if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and va:
                    delta = f"{(vb - va) / abs(va):+.2%}"
                rows.append(f"| `{k}` | `{_short(va)}` | `{_short(vb)}` | {delta} |")
        else:
            rows += ["", "_neither turn has measurements — `session run` has to be "
                     "executed before there is anything to compare._"]
        return "\n".join(rows)

    def _turn(self, i: int) -> Turn:
        for t in self.turns:
            if t.index == i:
                return t
        raise KeyError(f"turn {i} does not exist (available: "
                       f"{[t.index for t in self.turns]})")


def _now(when: datetime | None = None) -> str:
    return (when or datetime.now()).strftime("%Y-%m-%dT%H:%M:%S")


def _short(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    return "—" if v is None else str(v)


# =============================================================================
# CLI entry point — `python -m simbot.session …`
# =============================================================================
def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(prog="python -m simbot.session",
                                description="Session state (the state is in a file)")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("new", help="start a session")
    q.add_argument("spec")

    q = sub.add_parser("set", help="change a parameter → cost estimate only (no run)")
    q.add_argument("assignments", nargs="+", metavar="key=value")
    q.add_argument("--session", default=None)

    q = sub.add_parser("show", help="pick a session back up")
    q.add_argument("--session", default=None)

    q = sub.add_parser("run", help="run the current spec and record the result")
    q.add_argument("--prediction", default=None)
    q.add_argument("--runs-root", default=None)
    q.add_argument("--session", default=None)

    q = sub.add_parser("compare", help="compare two turns")
    q.add_argument("i", type=int)
    q.add_argument("j", type=int)
    q.add_argument("--session", default=None)

    q = sub.add_parser("list", help="list the sessions")

    a = p.parse_args(argv)

    if a.cmd == "list":
        for d in sorted(SESSIONS_ROOT.glob("*")):
            if (d / "session.yaml").exists():
                print(d.name)
        return 0

    if a.cmd == "new":
        s = Session.create(a.spec)
        print(f"session created: {s.session_id}")
        print(s.show())
        return 0

    s = Session.load(a.session) if a.session else Session.latest()

    if a.cmd == "set":
        t = s.set(a.assignments)
        print(f"turn {t.index} recorded — **nothing was run.**\n")
        for k, (old, new) in t.changes.items():
            print(f"  {k}: {old!r} → {new!r}")
        print()
        if "error" in t.cost:
            print(f"  cost cannot be estimated: {t.cost['error']}")
        else:
            print(f"  N={t.cost['n_particles']} · steps={t.cost['steps_per_seed']} "
                  f"· seeds {t.cost['n_seeds']} "
                  f"· concurrency {t.cost['concurrency']}")
            print(f"  estimated wall ≈ {t.cost['wall_s_batch']:g} s "
                  f"(budget {t.cost['budget_s']:g} s)")
            for key in ("action", "warning"):
                if key in t.cost:
                    print(f"  ⚠️  {t.cost[key]}")
        if t.problems:
            print("\n  convention violations:")
            for x in t.problems:
                print(f"    - {x}")
        print("\nto run it: python cli.py run "
              f"{s.derived_spec_file(t.index)}")
        return 0

    if a.cmd == "run":
        t = s.run(prediction=a.prediction, runs_root=a.runs_root)
        print(f"\nturn {t.index} recorded — run `{t.run_id or 'failed'}`")
        if t.problems:
            print("  problems:")
            for x in t.problems:
                print(f"    - {x}")
            return 1
        for k, v in list(t.metrics.items())[:6]:
            print(f"    {k:<22} {_short(v)}")
        prev = next((x for x in reversed(s.turns[:-1]) if x.kind == "run"), None)
        if prev:
            print(f"\ncompare against the previous run: python -m simbot.session "
                  f"compare "
                  f"{prev.index} {t.index}")
        return 0

    if a.cmd == "show":
        print(s.show())
        return 0

    if a.cmd == "compare":
        print(s.compare(a.i, a.j))
        return 0
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
