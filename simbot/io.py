"""run 디렉터리 · 해시 · 봉인. LLM 0줄.

**봉인(sealing)이 이 모듈의 존재 이유다.** 예측 문서를 실행 전에 해시로 고정하지
않으면 사후합리화를 구조적으로 막을 방법이 없다 — 결과를 보고 나서 예측을 손대도
아무 기록이 남지 않는다. S7 은 실행 전 해시와 대조하고, 불일치하면 **중단한다**.

`SEALED.sha256` 은 표준 `sha256sum` 형식이다 (`<hash>  <repo 상대경로>`).
`shasum -a 256 -c SEALED.sha256` 으로 이 코드 없이도 검증된다 — 봉인의 신뢰성이
우리 코드에 의존하지 않아야 한다.
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

# run 디렉터리 레이아웃 — master_plan.md §1.3
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

# 봉인 대상 — S5 실행 전에 존재해야 하고 이후 바뀌면 안 되는 문서
SEALED_STAGES: tuple[str, ...] = ("prediction", "intake", "spec")


# =============================================================================
# 해시
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
    """dict/list 를 정렬 직렬화해서 해시. spec_hash 등 구조체 해시용.

    `sort_keys=True` 가 필수다 — 키 순서가 바뀌었을 뿐인데 다른 계로 보이면
    캐시가 무효화되고 "같은 런"을 판정할 수 없다.
    """
    return sha256_text(json.dumps(obj, sort_keys=True, ensure_ascii=False,
                                  default=str))


def code_hash(root: Path | None = None) -> str:
    """`simbot` 소스 전체의 해시 (12자).

    코드가 바뀌면 결과 비교가 무의미해진다. `analysis/` 하위까지 포함한다 —
    분석 코드가 바뀌면 측정값이 바뀌므로 제외하면 봉인에 구멍이 난다.
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
    """파일 1개의 해시 (12자). **드라이버 스크립트를 특정하는 데 쓴다.**

    ★ `code_hash` 는 `simbot/` 만 덮는다. 그런데 런의 파라미터(`A` 목록·시드·런 길이·
      기하)를 정하는 것은 `scripts/` 의 드라이버다 — 그것이 해시에 없으면
      **무엇이 이 런을 만들었는지 산출물만으로 특정할 수 없다.**
      2026-07-29 `soft-r3-time-resolved` 에서 실제로 구멍이었다.
    """
    return sha256_file(path)[:12]




def git_rev(cwd: Path | None = None) -> str:
    """현재 커밋 (짧은 형식). git 이 없거나 repo 밖이면 `"?"`."""
    try:
        out = subprocess.run(["git", "rev-parse", "--short", "HEAD"],
                             capture_output=True, text=True, timeout=5,
                             cwd=cwd or REPO_ROOT)
        return out.stdout.strip() or "?"
    except Exception:
        return "?"


def git_dirty(cwd: Path | None = None) -> bool | None:
    """추적 중인 파일에 커밋 안 된 변경이 있는가. 판정 불가면 `None`.

    dirty 인 상태의 런은 `git_rev` 만으로 재현되지 않는다 — manifest 에 기록해서
    나중에 "이 런을 재현할 수 있는가"에 정직하게 답할 수 있게 한다.
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


# 재현성에 실제로 영향을 주는 패키지만. 전부 넣으면 무해한 변경에도 해시가 바뀐다.
ENV_PACKAGES: tuple[str, ...] = ("hoomd", "numpy", "scipy", "gsd", "freud")


def env_versions() -> dict[str, str]:
    """수치에 영향을 주는 패키지의 버전. import 실패는 `"absent"`."""
    out = {"python": platform.python_version(), "platform": platform.platform()}
    for name in ENV_PACKAGES:
        try:
            mod = __import__(name)
        except Exception:
            out[name] = "absent"
            continue
        v = getattr(mod, "__version__", None)
        if v is None:                                   # hoomd 는 version.version
            v = getattr(getattr(mod, "version", None), "version", "unknown")
        out[name] = str(v)
    return out


def env_hash() -> str:
    return sha256_payload(env_versions())[:12]


def provenance(driver: str | Path | None = None) -> dict:
    """"이 산출물을 무엇이 만들었는가" — **provenance 블록의 유일한 정의.**

    `build_manifest` 와 두 러너와 분석 스크립트가 모두 이것을 쓴다. 호출처마다
    손으로 만들면 키 이름이 갈라지고, 그러면 `report.reproducibility_section` 이
    읽는 이름과 러너가 쓰는 이름이 어긋나 **재현 정보가 조용히 빈칸으로 렌더된다.**
    2026-07-29 에 실제로 그랬다: `report.py` 는 `env_hash` 를 읽는데 두 러너가
    자기 manifest 를 손으로 만들면서 그 키를 넣지 않았다.

    **분석 단계에서도 호출한다.** 궤적의 manifest 는 *궤적 생성 시점*의 해시를 담고
    분석은 나중에 따로 돌 수 있다 (`--analyze-only`) — 그러면 `metrics.json` 이
    어느 분석 코드·어느 `freud` 에서 나왔는지 알 수 없다.

    Args:
        driver: 런을 정의한 스크립트. `code_hash` 는 `simbot/` 만 덮으므로
            `A` 목록·시드·런 길이·분석 창 같은 것을 정하는 `scripts/` 의 드라이버는
            이 인자로만 특정된다.
            ⚠ **파일 1개의 해시다.** 드라이버가 `scripts/` 안의 다른 모듈을
            import 하면 그것은 덮이지 않는다 — 현재 드라이버들은 `simbot` 만
            import 하므로 `code_hash` + `driver_hash` 가 함께 전부를 덮는다.
    """
    out = {
        "code_hash": code_hash(),
        "git_rev": git_rev(),
        "git_dirty": git_dirty(),
        "env_hash": env_hash(),
        "env": env_versions(),
    }
    if driver is not None:
        #  ★ 여러 파일을 받는다. 드라이버가 `scripts/` 안의 다른 모듈을 import 하면
        #    그것도 런의 파라미터를 정하므로 함께 해싱해야 한다 — 하나만 잡으면
        #    "code_hash + driver_hash 가 전부를 덮는다"는 주장이 거짓이 된다.
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
            #  합성 해시 — 파일 하나만 바뀌어도 달라진다
            out["driver_hash"] = sha256_payload(pairs)[:12]
    return out


# =============================================================================
# run_id
# =============================================================================
_SLUG_OK = re.compile(r"[^a-z0-9]+")


def slugify(text: str) -> str:
    return _SLUG_OK.sub("-", text.strip().lower()).strip("-") or "run"


def new_run_id(slug: str, spec_hash: str, when: _date | None = None) -> str:
    """`run_id = <날짜>_<슬러그>_<spec해시 앞6자리>`.

    `when` 을 명시하면 결정적이다 (테스트·재현용). 생략하면 오늘.
    """
    day = (when or _date.today()).isoformat()
    return f"{day}_{slugify(slug)}_{spec_hash[:6]}"


# =============================================================================
# run 디렉터리
# =============================================================================
@dataclass(frozen=True)
class RunDir:
    """run 디렉터리의 경로 계산기. 파일을 만들지 않고 경로만 안다.

    `RUN_LAYOUT` 의 키로 접근한다 — 파일명 문자열이 코드 곳곳에 흩어지면
    S6 이 쓰는 이름과 S7 이 읽는 이름이 조용히 달라진다.
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
        """이미 산출물이 있는 단계들 — `resume` 이 재계산을 건너뛰는 근거."""
        return [s for s in RUN_LAYOUT if s != "input" and self.exists(s)]


# =============================================================================
# 봉인
# =============================================================================
@dataclass(frozen=True)
class SealEntry:
    digest: str
    relpath: str


@dataclass
class SealVerdict:
    """봉인 검증 결과. **`ok=False` 면 S7 은 중단해야 한다.**"""

    ok: bool
    changed: list[str] = field(default_factory=list)
    missing: list[str] = field(default_factory=list)
    unsealed: list[str] = field(default_factory=list)   # 봉인 대상인데 목록에 없음
    entries: dict[str, str] = field(default_factory=dict)

    def summary(self) -> str:
        if self.ok:
            n = len(self.entries)
            return f"봉인 검증 통과 — {n}개 문서 실행 후 미변경"
        parts = []
        if self.changed:
            parts.append(f"변경됨 {self.changed}")
        if self.missing:
            parts.append(f"사라짐 {self.missing}")
        if self.unsealed:
            parts.append(f"봉인 안 됨 {self.unsealed}")
        return "봉인 위반 — " + " · ".join(parts)


def _seal_relpath(path: Path) -> str:
    """repo 상대경로. repo 밖이면 절대경로 (테스트용 tmpdir 등)."""
    p = Path(path).resolve()
    try:
        return p.relative_to(REPO_ROOT).as_posix()
    except ValueError:
        return p.as_posix()


def write_seal(rundir: RunDir, stages: tuple[str, ...] = SEALED_STAGES) -> Path:
    """존재하는 봉인 대상 문서의 해시를 `SEALED.sha256` 에 쓴다.

    ⚠ **S5 실행 전에 호출해야 한다.** 실행 후에 봉인하면 봉인이 아무것도 보증하지
      않는다 — 이 함수는 그것을 감지할 수 없으므로 파이프라인이 순서를 지켜야 한다.
    """
    lines = []
    for stage in stages:
        p = rundir.file(stage)
        if p.exists():
            lines.append(f"{sha256_file(p)}  {_seal_relpath(p)}")
    if not lines:
        raise FileNotFoundError(
            f"봉인할 문서가 없다 — {list(stages)} 중 하나도 존재하지 않는다. "
            f"S2 예측을 먼저 쓸 것.")
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
        if not rel:                                    # 한 칸 구분자도 허용
            digest, _, rel = line.partition(" ")
        out[rel.strip()] = digest.strip()
    return out


def verify_seal(rundir: RunDir, stages: tuple[str, ...] = SEALED_STAGES) -> SealVerdict:
    """봉인된 문서가 실행 후 바뀌지 않았는지 확인한다.

    `stages` 중 봉인 목록에 아예 없는 문서는 `unsealed` 로 보고한다 —
    "봉인 파일이 통과했다"가 "예측이 봉인됐다"를 뜻하지 않게 하려면 이것이 필요하다.
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
    """`05_run_manifest.json` 의 내용. 재현에 필요한 전부.

    예측 해시를 여기에 박아두면 `SEALED.sha256` 파일이 지워져도 봉인이 남는다
    (manifest 도 지우면 사라지지만, 두 곳을 동시에 고쳐야 한다는 마찰이 생긴다).
    """
    man = {
        "run_id": run_id,
        "spec_hash": spec_hash,
        "seed": seed,
        #  ★ provenance 는 한 곳에서만 만든다 — 여기서 손으로 나열하면 러너·분석과
        #    키 이름이 갈라진다 (`provenance()` docstring 참조)
        **provenance(),
    }
    if rundir is not None and rundir.exists("seal"):
        man["sealed"] = read_seal(rundir)
    man.update(extra or {})
    return man
