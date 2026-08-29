"""대화 세션 상태. LLM 0줄.

> **상태는 대화가 아니라 파일에 있다.** 대화는 사라지고 검사할 수 없다.

`sessions/<id>/session.yaml` 이 턴 단위 이력을 갖는다 (`runs/<id>/` 는 런 단위).
대화 컨텍스트가 날아가도 `show` 로 이어받을 수 있어야 한다.

## `set` 은 실행하지 않는다

`set` 은 파라미터를 바꾸고 **비용 추정만** 한다. 실행은 `run` 이 한다.
이 분리가 없으면 "N 을 8000 으로 해보자"가 11분짜리 런을 조용히 시작한다.

## 이력은 append-only

턴을 고치지 않고 새 턴을 쌓는다. 과거 턴을 덮어쓰면 "무엇을 시도했고 무엇이
안 됐는가"가 사라진다 — 그것이 이 프로젝트 산출물의 절반이다.
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
# dotted-key 설정
# =============================================================================
def _navigate(root, path: list[str]):
    """`['species', '0', 'radius_si']` 를 따라가 (부모, 마지막키) 를 돌려준다."""
    cur = root
    for part in path[:-1]:
        if isinstance(cur, list):
            cur = cur[int(part)]
        elif isinstance(cur, dict):
            cur = cur[part]
        elif is_dataclass(cur) and part in {f.name for f in fields(cur)}:
            cur = getattr(cur, part)
        else:
            raise KeyError(f"경로 {'.'.join(path)} 에서 {part!r} 를 찾을 수 없다")
        if cur is None:
            raise KeyError(f"경로 {'.'.join(path)} 의 {part!r} 가 None 이다")
    return cur, path[-1]


def _read(parent, key: str):
    if isinstance(parent, list):
        return parent[int(key)]
    if isinstance(parent, dict):
        return parent[key]
    return getattr(parent, key)


def _parse_scalar(text: str):
    """문자열 → int/float/bool/str/list.

    ★ **`yaml.safe_load` 를 그대로 쓰면 안 된다.** YAML 1.1 의 float 정규식은
      만티사에 소수점을 요구하므로 `5e-3`·`2e-5` 가 **문자열**로 읽힌다
      (`2.5e-3` 은 float). 그러면 `dt_star` 에 문자열이 들어가고 무차원화가
      TypeError 로 죽거나, 더 나쁘게는 비교 연산이 조용히 통과한다.
      같은 함정을 `config/run_policy.yaml` 의 `6.3e6` 에서도 겪었다 —
      `simbot/policy.py::_find_numeric_strings` 참조.

    따라서 **숫자를 먼저 시도**하고, 아닐 때만 YAML 에 넘긴다.
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
    """`numerics.dt_star=2.5e-3` 을 적용하고 `(이전값, 새값)` 을 돌려준다.

    `Quantity` 필드면 `value` 만 바꾸고 **provenance 를 `user` 로 바꾼다** —
    사람이 이 런에서 직접 지정한 값은 가정도 추론도 정책도 아니다. 이 구별이
    S7b 감도 분석의 대상 목록을 정한다 (사람이 정한 값은 흔들 대상이 아니다).
    """
    parent, key = _navigate(spec, dotted.split("."))
    old = _read(parent, key)
    new_value = _parse_scalar(raw)

    if isinstance(old, Quantity):
        prev = old.value
        old.value = new_value
        old.provenance = "user"
        old.basis = (f"턴 {turn} 에서 사람이 지정 (이전 {prev!r}"
                     + (f", 원래 근거: {old.basis}" if old.basis else "") + ")")
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
# 턴
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
# 세션
# =============================================================================
@dataclass
class Session:
    session_id: str
    spec_path: str
    spec: SystemSpec
    turns: list[Turn] = field(default_factory=list)
    root: Path = SESSIONS_ROOT

    # --- 경로 ---
    @property
    def dir(self) -> Path:
        return Path(self.root) / self.session_id

    @property
    def state_file(self) -> Path:
        return self.dir / "session.yaml"

    def derived_spec_file(self, turn: int) -> Path:
        return self.dir / f"spec_turn{turn:02d}.yaml"

    # --- 생성 / 저장 / 로드 ---
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
            raise FileNotFoundError(f"세션 {session_id} 가 없다 ({f})")
        raw = yaml.safe_load(f.read_text(encoding="utf-8"))
        turns = [Turn(**t) for t in raw.get("turns", [])]
        # 최신 턴의 파생 spec 을 정본으로 읽는다 — 원본 spec 은 바뀌지 않는다
        spec_file = base / session_id / raw["current_spec"]
        return cls(session_id=session_id, spec_path=raw["spec_path"],
                   spec=SystemSpec.load(spec_file), turns=turns, root=base)

    @classmethod
    def latest(cls, *, root: Path | None = None) -> Session:
        base = Path(root) if root else SESSIONS_ROOT
        cands = sorted(p.name for p in base.glob("*") if (p / "session.yaml").exists())
        if not cands:
            raise FileNotFoundError(f"{base} 에 세션이 없다. `session new <spec>` 먼저.")
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

    # --- 명령 ---
    def set(self, assignments: list[str], *, policy: Policy | None = None,
            when: datetime | None = None) -> Turn:
        """파라미터 변경 → 파생 spec + **비용 추정만. 실행하지 않는다.**"""
        changes: dict[str, list] = {}
        for a in assignments:
            if "=" not in a:
                raise ValueError(f"{a!r} 은 `키=값` 형식이 아니다")
            key, _, raw = a.partition("=")
            old, new = set_by_path(self.spec, key.strip(), raw.strip(),
                                   turn=len(self.turns))
            changes[key.strip()] = [old, new]

        rep = validate_spec(self.spec)
        return self._append(Turn(
            index=len(self.turns), timestamp=_now(when), kind="set",
            changes=changes, spec_hash=self.spec.hash(),
            cost=self.estimate_cost(policy=policy), problems=rep.problems,
            note="비용 추정만 — 실행하려면 `session run`"))

    def estimate_cost(self, *, policy: Policy | None = None) -> dict:
        """이 spec 을 지금 돌리면 얼마나 걸리는가. **예산 초과면 그렇게 적는다.**"""
        pol = policy or load_policy()
        from .nondim import reduce_spec
        try:
            r = reduce_spec(self.spec)
        except Exception as e:                        # 척도 미등록 카드 등
            return {"error": f"무차원화 실패: {e}"}

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
            "wall_s_batch": round(wall, 3),           # 동시 실행이므로 배치 ≈ 1개
            "budget_s": budget,
            "over_budget": wall > budget,
        }
        if out["over_budget"]:
            out["action"] = (f"예산 {budget:g} s 초과 — 정책상 **실행하지 않고 "
                             f"보고한다** (run_policy §5)")
        if r.n_particles < 500:
            out["warning"] = ("N < 500 — 처리량 모델은 N ≥ 500 에서 실측됐다. "
                              "이 추정은 과소추정이다 (pilot 으로 실측할 것)")
        return out

    def record_run(self, *, run_id: str, metrics: dict, wall_s: float | None = None,
                   problems: list[str] | None = None,
                   when: datetime | None = None) -> Turn:
        """실행 결과를 턴에 기록한다. 실행 자체는 `cli.py` 가 한다."""
        return self._append(Turn(
            index=len(self.turns), timestamp=_now(when), kind="run",
            spec_hash=self.spec.hash(), run_id=run_id, metrics=metrics,
            cost=({"wall_s_measured": round(wall_s, 3)} if wall_s else {}),
            problems=list(problems or [])))

    def run(self, *, prediction: str | Path | None = None,
            runs_root: str | Path | None = None, extra_args: list[str] | None = None,
            when: datetime | None = None) -> Turn:
        """현재 턴의 파생 spec 으로 `cli.py run` 을 실행하고 결과를 턴에 기록한다.

        ★ **`set` 과 분리되어 있다.** `set` 은 추정만, `run` 은 실행. 이 경계가 없으면
          "N 을 8000 으로 해보자"가 11분짜리 런을 조용히 시작한다.

        실패하면 `problems` 에 stderr 를 담아 **run 턴을 그대로 기록한다** —
        실패한 시도가 이력에서 사라지면 "무엇을 시도했는가"의 절반이 없어진다.
        """
        import subprocess
        import sys

        turn_index = self.turns[-1].index
        spec_file = self.derived_spec_file(turn_index)
        # ★ run_id 를 세션이 정한다. `cli.py` 의 기본 슬러그는 spec 의 부모 디렉터리명인데,
        #   세션 spec 은 `sessions/<날짜>_<카드>/` 에 있어서 날짜가 두 번 들어간다.
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
                problems=[f"cli.py run 실패 (exit {proc.returncode})",
                          *(proc.stderr.strip().splitlines()[-5:] or ["(stderr 없음)"])],
                note="실패한 시도도 이력에 남는다"))

        metrics, wall = {}, None
        for line in proc.stdout.splitlines():
            if "실측 wall" in line:
                try:
                    wall = float(line.split("실측 wall")[1].split()[0])
                except (IndexError, ValueError):
                    pass
        root = Path(runs_root) if runs_root else REPO_ROOT / "runs"
        mfile = root / run_id / "metrics.json"
        if mfile.exists():
            raw = json.loads(mfile.read_text(encoding="utf-8"))
            metrics = {k: v["value"] for k, v in raw.items()
                       if isinstance(v, dict) and "value" in v}
        return self.record_run(run_id=run_id, metrics=metrics, wall_s=wall, when=when)

    # --- 보기 ---
    def show(self) -> str:
        """세션 이어받기용 요약. 대화 컨텍스트가 없어도 읽힌다."""
        out = [f"# 세션 `{self.session_id}`", "",
               f"**카드** `{self.spec.card}`",
               f"**질문** {self.spec.question.strip()}",
               f"**원본 spec** `{self.spec_path}`",
               f"**턴 {len(self.turns)}개** · 현재 spec_hash `{self.spec.hash()[:12]}`",
               "", "| 턴 | 시각 | 종류 | 내용 | 문제 |", "|---|---|---|---|---|"]
        for t in self.turns:
            if t.kind == "set":
                what = " · ".join(f"`{k}`: {v[0]!r} → {v[1]!r}"
                                  for k, v in t.changes.items())
                if t.cost.get("wall_s_batch") is not None:
                    what += f" (추정 `{t.cost['wall_s_batch']:g} s`)"
                    if t.cost.get("over_budget"):
                        what += " ⚠️ **예산 초과**"
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
            out += ["", "### 미해결 문제", *[f"- (턴 {i}) {p}" for i, p in probs]]
        out += ["", "다음: `session set k.v=…` (추정만) 또는 `session run` (실행)"]
        return "\n".join(out)

    def compare(self, i: int, j: int) -> str:
        """턴 두 개의 파라미터·측정값 대조표."""
        a, b = self._turn(i), self._turn(j)
        rows = [f"# 턴 {i} vs 턴 {j}", "",
                f"| | 턴 {i} | 턴 {j} |", "|---|---|---|",
                f"| 시각 | {a.timestamp} | {b.timestamp} |",
                f"| 종류 | {a.kind} | {b.kind} |",
                f"| spec_hash | `{a.spec_hash[:12]}` | `{b.spec_hash[:12]}` |",
                f"| run_id | `{a.run_id or '—'}` | `{b.run_id or '—'}` |"]
        if a.spec_hash == b.spec_hash:
            rows.append("| **차이** | 같은 spec — 파라미터가 동일하다 | |")

        # 파라미터 차이 — 두 턴 사이의 모든 set 을 누적한다
        lo, hi = min(i, j), max(i, j)
        acc: dict[str, list] = {}
        for t in self.turns[lo + 1:hi + 1]:
            for k, v in t.changes.items():
                acc[k] = [acc.get(k, v)[0], v[1]]
        if acc:
            rows += ["", "## 파라미터 변경", "", "| 키 | 이전 | 이후 |", "|---|---|---|"]
            rows += [f"| `{k}` | `{v[0]}` | `{v[1]}` |" for k, v in acc.items()]

        keys = sorted(set(a.metrics) | set(b.metrics))
        if keys:
            rows += ["", "## 측정값", "", f"| 양 | 턴 {i} | 턴 {j} | 변화 |",
                     "|---|---|---|---|"]
            for k in keys:
                va, vb = a.metrics.get(k), b.metrics.get(k)
                delta = ""
                if isinstance(va, (int, float)) and isinstance(vb, (int, float)) and va:
                    delta = f"{(vb - va) / abs(va):+.2%}"
                rows.append(f"| `{k}` | `{_short(va)}` | `{_short(vb)}` | {delta} |")
        else:
            rows += ["", "_두 턴에 측정값이 없다 — `session run` 으로 실행해야 "
                     "대조할 것이 생긴다._"]
        return "\n".join(rows)

    def _turn(self, i: int) -> Turn:
        for t in self.turns:
            if t.index == i:
                return t
        raise KeyError(f"턴 {i} 가 없다 (있는 것: "
                       f"{[t.index for t in self.turns]})")


def _now(when: datetime | None = None) -> str:
    return (when or datetime.now()).strftime("%Y-%m-%dT%H:%M:%S")


def _short(v) -> str:
    if isinstance(v, float):
        return f"{v:.6g}"
    return "—" if v is None else str(v)


# =============================================================================
# CLI 진입점 — `python -m simbot.session …`
# =============================================================================
def main(argv: list[str] | None = None) -> int:
    import argparse
    p = argparse.ArgumentParser(prog="python -m simbot.session",
                                description="세션 상태 관리 (상태는 파일에 있다)")
    sub = p.add_subparsers(dest="cmd", required=True)

    q = sub.add_parser("new", help="세션 시작")
    q.add_argument("spec")

    q = sub.add_parser("set", help="파라미터 변경 → 비용 추정만 (실행 안 함)")
    q.add_argument("assignments", nargs="+", metavar="키=값")
    q.add_argument("--session", default=None)

    q = sub.add_parser("show", help="세션 이어받기")
    q.add_argument("--session", default=None)

    q = sub.add_parser("run", help="현재 spec 으로 실행하고 결과를 턴에 기록")
    q.add_argument("--prediction", default=None)
    q.add_argument("--runs-root", default=None)
    q.add_argument("--session", default=None)

    q = sub.add_parser("compare", help="턴 대조")
    q.add_argument("i", type=int)
    q.add_argument("j", type=int)
    q.add_argument("--session", default=None)

    q = sub.add_parser("list", help="세션 목록")

    a = p.parse_args(argv)

    if a.cmd == "list":
        for d in sorted(SESSIONS_ROOT.glob("*")):
            if (d / "session.yaml").exists():
                print(d.name)
        return 0

    if a.cmd == "new":
        s = Session.create(a.spec)
        print(f"세션 생성: {s.session_id}")
        print(s.show())
        return 0

    s = Session.load(a.session) if a.session else Session.latest()

    if a.cmd == "set":
        t = s.set(a.assignments)
        print(f"턴 {t.index} 기록 — **실행하지 않았다.**\n")
        for k, (old, new) in t.changes.items():
            print(f"  {k}: {old!r} → {new!r}")
        print()
        if "error" in t.cost:
            print(f"  비용 추정 불가: {t.cost['error']}")
        else:
            print(f"  N={t.cost['n_particles']} · steps={t.cost['steps_per_seed']} "
                  f"· 시드 {t.cost['n_seeds']} · 동시 {t.cost['concurrency']}")
            print(f"  추정 wall ≈ {t.cost['wall_s_batch']:g} s "
                  f"(예산 {t.cost['budget_s']:g} s)")
            for key in ("action", "warning"):
                if key in t.cost:
                    print(f"  ⚠️  {t.cost[key]}")
        if t.problems:
            print("\n  규약 위반:")
            for x in t.problems:
                print(f"    - {x}")
        print("\n실행하려면: python cli.py run "
              f"{s.derived_spec_file(t.index)}")
        return 0

    if a.cmd == "run":
        t = s.run(prediction=a.prediction, runs_root=a.runs_root)
        print(f"\n턴 {t.index} 기록 — run `{t.run_id or '실패'}`")
        if t.problems:
            print("  문제:")
            for x in t.problems:
                print(f"    - {x}")
            return 1
        for k, v in list(t.metrics.items())[:6]:
            print(f"    {k:<22} {_short(v)}")
        prev = next((x for x in reversed(s.turns[:-1]) if x.kind == "run"), None)
        if prev:
            print(f"\n이전 런과 대조: python -m simbot.session compare "
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
