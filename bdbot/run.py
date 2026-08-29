r"""L4 — 스펙만 읽고 실행하는 층 + **수치 건전성 판정** (마스터플랜 §12.5 · §16 Phase 4).

⭐️ **L4는 물리 검증기가 아닙니다** (사용자 지시 2026-08-04, mater_plan §16 방향 보강).
   런 후에 보는 것은 **발산 · NaN/Inf · 이상수렴** 뿐입니다. 해석해 대조·문헌 대조는
   `bdbot.metrics` 의 `role` 체계가 이미 담당하고, 엄밀성 투자는 앞단(L0→L3)에 몰아줍니다.

## 무엇이 여기 있고 무엇이 케이스에 남는가

여기(L4, 계와 무관):
  · `Guard`     런타임 감시 — 비유한 PE/힘 · PE 폭발 · 한 스텝 과대 변위 (§12.5)
  · `judge`     사후 판정 — nan / diverged / frozen / drifting / ok
  · `execute`   런 디렉토리 생명주기 · 평형화+생산 루프 · 표본 수집 · metrics 방출

케이스(고유 물리):
  · `build(spec) -> Build`  스펙의 무차원 파라미터로 HOOMD 계를 세우고 표본 함수를 준다

**데이터 계약은 스펙 하나입니다** — `build` 는 `LoadedSpec` 만 받고 케이스 YAML 을
다시 읽지 않습니다. 물리를 세우는 **코드**는 계마다 다를 수밖에 없어 레지스트리로
붙입니다 (`@builder("케이스이름")`). `nondim show` 가 자족성을 보장하는 것은
**데이터** 쪽이고, 여기서 깨지지 않습니다.

## 판정 기준 — 왜 이 넷인가

| 상태 | 무엇을 봤나 | 왜 수치 문제인가 |
|---|---|---|
| `nan` | PE 또는 힘이 비유한 | 적분이 깨졌다 |
| `diverged` | \|PE\|/N 폭발, 또는 한 스텝 변위 > `step_disp_max` | dt 가 너무 크다 |
| `frozen` | 생산 구간 PE/N 의 분산이 0 | 열잡음이 있는데 안 움직인다 — 힘이 0이 됐다는 뜻 (pair.Table 함정 11이 정확히 이 모양) |
| `drifting` | 전반/후반 평균차 > 3·블록SEM | 평형화 부족 (경고 — 실패 아님) |

`frozen` 은 이 프로젝트의 실제 함정에서 나왔습니다: `pair.Table` 은 `r < r_min` 에서
힘·에너지가 **0** 이라 입자가 겹친 채로 조용히 멈춥니다. 터지지 않으므로 NaN 검사로는
안 잡히고, PE 도 폭발하지 않습니다. **변하지 않는다는 것**이 유일한 신호입니다.
"""
from __future__ import annotations

import json
import math
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Callable

import numpy as np

from . import metrics as MET
from . import runid as RID
from . import stats as ST

SCHEMA = "bdbot.run/0.1"

OK, DIVERGED, NAN, FROZEN = "ok", "diverged", "nan", "frozen"
DRIFT_SIGMA = 3.0            # 전반/후반 평균차가 블록SEM의 이 배수를 넘으면 평형화 부족
GUARD_EVERY = 10_000         # §12.5 — 10⁴ 스텝마다


# ══════════════════════════════════════════════════════════════════════
# 케이스가 돌려주는 것
# ══════════════════════════════════════════════════════════════════════
@dataclass
class Phase:
    """런의 한 구간. **구동 프로토콜이 있는 계**를 위해 도입했습니다.

    `trap-drag` 의 질문이 "평형 에너지 vs 끌 때 에너지 vs 멈춘 뒤 이완"이라
    한 궤적 안에 성격이 다른 구간이 셋 필요해졌습니다. eq→prod 두 단계로는
    표현할 수 없습니다.

    ⭐️ `expect_steady=False` 가 중요합니다. 이완 구간은 **에너지가 변하는 것이 물리**라서
       표류 검사를 돌리면 거짓 경고가 납니다. "정상상태를 기대하는 구간"에서만 판정합니다.
    """

    name: str
    n_steps: int
    sample_every: int = 0            # 0 이면 Build.sample_every
    collect: bool = True             # False 면 표본을 버린다 (초기 완화)
    expect_steady: bool = True       # False 면 표류 검사를 건너뛴다
    note: str = ""


@dataclass
class Build:
    """케이스가 스펙에서 만들어 주는 것. L4는 이것만 받아 돌립니다."""

    sim: Any                                   # hoomd.Simulation
    forces: list                               # 진단용 (힘 배열을 읽는다)
    n_particles: int
    sample: Callable[..., dict]                # (timestep, phase) -> {이름: 값}
    pe_per_particle: Callable[[], float]       # ⟨U⟩/N [kT] — 평형·frozen 판정의 지표
    n_eq: int = 0
    n_prod: int = 0
    sample_every: int = 1
    phases: list = field(default_factory=list)   # 비면 n_eq/n_prod 로 만든다
    gsd_path: Any = None
    tags: list = field(default_factory=list)   # metrics 의 system_tags
    physical: dict = field(default_factory=dict)
    finalize: Callable[[dict], dict] | None = None    # 표본 → 관측량·extra

    def plan(self) -> list:
        """실행할 구간 목록. 구간을 안 준 케이스는 예전처럼 평형화→생산 두 개."""
        if self.phases:
            return list(self.phases)
        return [Phase("평형화", self.n_eq, GUARD_EVERY, collect=False),
                Phase("생산", self.n_prod, self.sample_every)]


# ══════════════════════════════════════════════════════════════════════
# 레지스트리 — 케이스가 자기 빌더를 등록한다
# ══════════════════════════════════════════════════════════════════════
BUILDERS: dict = {}


def builder(case: str):
    def deco(fn):
        BUILDERS[case] = fn
        return fn
    return deco


def get_builder(case: str):
    if case not in BUILDERS:
        raise KeyError(f"'{case}' 의 L4 빌더가 등록되지 않았습니다. "
                       f"등록된 것: {', '.join(sorted(BUILDERS)) or '(없음)'}")
    return BUILDERS[case]


# ══════════════════════════════════════════════════════════════════════
# §12.5 런타임 가드
# ══════════════════════════════════════════════════════════════════════
class Diverged(RuntimeError):
    """가드가 걸렸을 때. `status` 로 어떤 종류인지 구분합니다."""

    def __init__(self, status: str, step: int, msg: str):
        super().__init__(f"[{status}] step {step}: {msg}")
        self.status, self.step, self.msg = status, step, msg


@dataclass
class Guard:
    """`10⁴` 스텝마다 비유한값·폭발·과대변위를 본다 (마스터플랜 §12.5).

    `step_disp_max` 는 **한 스텝의 힘 구동 변위** `dt·|F|/γ` 의 상한입니다 (무차원 γ=1).
    BD 는 O(dt) 라 한 스텝에 지름의 10% 를 넘게 밀리면 적분이 이미 신뢰할 수 없습니다.
    ⚠️ 열잡음 변위 `√(2 dt)` 는 여기 포함하지 않습니다 — 그건 물리이고, 판정 대상은
       **힘 항** 입니다.
    """

    dt_star: float
    n_particles: int
    pe_abs_max: float = 1e8          # |PE|/N [kT] — 이 위는 물리가 아니라 폭발
    step_disp_max: float = 0.1       # [σ]
    trips: list = field(default_factory=list)
    # ★ 런 전체의 **최악값**. 마지막 표본만 남기면 중간의 최악 순간을 놓칩니다 —
    #   `dt/τ_fast` 판정은 평균이 아니라 최악값으로 해야 합니다 (안정성 기준).
    f_max_seen: float = 0.0
    step_disp_seen: float = 0.0

    def check(self, step: int, pe_per_n: float, forces) -> None:
        if not math.isfinite(pe_per_n):
            raise Diverged(NAN, step, f"PE/N 이 비유한값입니다 ({pe_per_n!r})")
        if abs(pe_per_n) > self.pe_abs_max:
            raise Diverged(DIVERGED, step,
                           f"PE/N = {pe_per_n:.4e} kT 가 상한 {self.pe_abs_max:.0e} 초과")
        fmax = 0.0
        for f in forces:
            arr = np.asarray(f.forces, dtype=float)
            if arr.size == 0:
                continue
            if not np.all(np.isfinite(arr)):
                raise Diverged(NAN, step, f"{type(f).__name__} 의 힘 배열에 비유한값")
            fmax = max(fmax, float(np.abs(arr).max()))
        disp = self.dt_star * fmax
        self.f_max_seen = max(self.f_max_seen, fmax)
        self.step_disp_seen = max(self.step_disp_seen, disp)
        if disp > self.step_disp_max:
            raise Diverged(DIVERGED, step,
                           f"한 스텝 힘 변위 dt·|F|max = {disp:.4g} σ 가 상한 "
                           f"{self.step_disp_max} σ 초과 (|F|max = {fmax:.4g} kT/σ). "
                           f"dt 를 줄이세요")
        self.trips.append({"step": step, "pe_per_n": pe_per_n, "f_max": fmax,
                           "step_disp": disp})


# ══════════════════════════════════════════════════════════════════════
# 사후 판정 — 수치만
# ══════════════════════════════════════════════════════════════════════
def judge(pe_series, *, status: str = OK, expect_steady: bool = True,
          label: str = "생산") -> dict:
    """`⟨U⟩/N` 시계열 하나에서 수치 건전성을 판정.

    `frozen` 이 핵심입니다 — 열잡음이 있는 계에서 PE/N 이 **정확히 상수**면 힘이 0이
    됐다는 뜻이고, 그건 터지지 않으므로 NaN·폭발 검사로는 잡히지 않습니다.

    ⭐️ `expect_steady=False` 면 표류 검사를 건너뜁니다 — 이완 구간처럼 **에너지가 변하는
       것이 물리인** 구간에 표류 경고를 내면 발견을 경고로 부르게 됩니다 (규칙 7').
    """
    s = np.asarray(pe_series, dtype=float)
    out = {"status": status, "n_samples": int(s.size), "warnings": []}
    if status != OK:
        return out
    if s.size == 0:
        out["status"] = NAN
        out["warnings"].append(f"{label} 표본이 0개입니다")
        return out
    if not np.all(np.isfinite(s)):
        out["status"] = NAN
        return out

    out["pe_mean"] = float(s.mean())
    out["pe_std"] = float(s.std())
    if s.size >= 2 and out["pe_std"] == 0.0:
        out["status"] = FROZEN
        out["warnings"].append(
            f"[{label}] ⟨U⟩/N 이 {s.size}개 표본에서 정확히 상수({s[0]:.6g})입니다 — "
            "열잡음이 있는데 변하지 않습니다. 힘이 0이 됐을 가능성 "
            "(pair.Table 함정 11: r<r_min 에서 F=0)")
        return out

    # 평형화 판정 — 전반/후반 평균차 vs 블록 SEM. 실패가 아니라 경고입니다.
    if expect_steady and s.size >= 8:
        h = s.size // 2
        d = float(s[h:].mean() - s[:h].mean())
        sem = float(ST.block_sem(s))
        out["drift"] = d
        out["drift_sem"] = sem
        out["drift_sigma"] = abs(d) / sem if sem > 0 else float("inf")
        if sem > 0 and abs(d) > DRIFT_SIGMA * sem:
            out["warnings"].append(
                f"[{label}] ⟨U⟩/N 이 표류합니다: 후반−전반 = {d:+.4g} kT "
                f"= {abs(d)/sem:.1f}·블록SEM (기준 {DRIFT_SIGMA:g}). 평형화가 부족합니다")
    return out


# ══════════════════════════════════════════════════════════════════════
# 실행 — 스펙 하나 → 런 디렉토리 하나
# ══════════════════════════════════════════════════════════════════════
def execute(spec, build_fn, outdir, *, force: bool = False, progress: bool = True,
            guard_every: int = GUARD_EVERY, extra_metrics=None) -> dict:
    """스펙 하나를 실행하고 `metrics.json` 을 남깁니다. 반환은 판정 dict.

    `spec` 은 `bdbot.nondim.LoadedSpec` — 즉 `specs/<run_id>.json` 을 되읽은 것입니다.
    """
    from . import sim as SIM       # hoomd 를 여기서만 끌어온다

    outdir = Path(outdir)
    go, msg = RID.prepare_outdir(outdir, force)
    if not go:
        print(msg)
        return {"status": "skipped", "run_id": spec.run_id}

    ok_hash, want = spec.verify_hash()
    if not ok_hash:
        raise ValueError(f"손으로 고친 스펙입니다 — run_id 가 내용과 맞지 않습니다 "
                         f"(기대 {want}). §16 규칙 2.")
    if spec.verdict.startswith("FAIL"):
        raise ValueError(f"L3 판정이 FAIL 인 스펙은 실행하지 않습니다: {spec.verdict}")

    b = build_fn(spec, outdir)          # ★ outdir 을 넘긴다 — GSD 경로를 스펙에 넣지 않으려고
    dt_star = float(spec.numerics["dt_star"])
    guard = Guard(dt_star=dt_star, n_particles=b.n_particles)
    t0 = time.time()
    status, trip = OK, None

    def _loop(n_steps, label, collect, every):
        """평형화·생산 공통 루프. `every` 스텝마다 멈춰 가드를 돌리고 표본을 받는다.

        ★ 가드 주기와 표본 주기를 **같게** 맞춥니다 (둘 중 작은 쪽). 따로 돌리면 표본
          사이에서 발산한 것을 표본에 남기지 못합니다.
        """
        nonlocal status, trip
        done = 0
        every = max(1, min(every, n_steps)) if n_steps else 1
        next_print = 0
        while done < n_steps:
            chunk = min(every, n_steps - done)
            b.sim.run(chunk)
            done += chunk
            pe = b.pe_per_particle()
            try:
                guard.check(b.sim.timestep, pe, b.forces)
            except Diverged as e:
                status, trip = e.status, {"step": e.step, "msg": e.msg, "phase": label}
                print(f"\n  ✗ 런타임 가드: {e}", flush=True)
                return False, done
            collect(b.sim.timestep, pe)
            if progress and done >= next_print:
                print(SIM.progress(done, n_steps, time.time() - t0,
                                   f"{label}  ⟨U⟩/N={pe:9.4f}"), flush=True)
                next_print += max(1, n_steps // 10)
        return True, done

    # ── 구간을 순서대로 ────────────────────────────────────────────────
    plan = b.plan()
    samples: list = []
    pe_series: list = []          # 판정용 — 정상상태를 기대하는 마지막 구간
    phase_of: list = []
    eq_trace: list = []
    per_phase: dict = {}
    n_done = 0
    alive = True
    for ph in plan:
        if not alive or ph.n_steps <= 0:
            continue
        pe_ph: list = []
        print(f"\n  ▸ {ph.name}  {ph.n_steps:,} 스텝"
              + (f"   {ph.note}" if ph.note else ""), flush=True)

        def collect(ts, pe, _ph=ph, _pe=pe_ph):
            _pe.append(pe)
            if _ph.collect:
                s = b.sample(ts, _ph.name)
                s["_t_step"] = ts
                samples.append(s)
                phase_of.append(_ph.name)
                pe_series.append(pe)
            else:
                eq_trace.append((ts, pe))

        every = ph.sample_every or b.sample_every
        alive, done = _loop(ph.n_steps, ph.name, collect,
                            min(every, guard_every) if ph.collect else guard_every)
        n_done += done
        per_phase[ph.name] = judge(pe_ph, status=OK if alive else status,
                                   expect_steady=ph.expect_steady, label=ph.name)
        per_phase[ph.name]["steps"] = done
        per_phase[ph.name]["expect_steady"] = ph.expect_steady

    SIM.flush_writers(b.sim)
    wall = time.time() - t0

    # 전체 판정 = 구간 판정의 최악값. 표류는 구간별로 이미 봤으므로 여기선 넘긴다.
    verdict = judge(pe_series, status=status, expect_steady=False)
    for name, pv in per_phase.items():
        if pv["status"] != OK and verdict["status"] == OK:
            verdict["status"] = pv["status"]
        verdict["warnings"].extend(pv.get("warnings", []))
    verdict["phases"] = per_phase
    verdict.update(run_id=spec.run_id, wall_seconds=wall,
                   steps_done=n_done, steps_planned=sum(p.n_steps for p in plan),
                   trip=trip, guard_samples=len(guard.trips))
    if guard.trips:
        last = guard.trips[-1]
        verdict["f_max_last"] = last["f_max"]
        verdict["step_disp_last"] = last["step_disp"]
        verdict["f_max_seen"] = guard.f_max_seen
        verdict["step_disp_seen"] = guard.step_disp_seen

    # ── 산출물 ────────────────────────────────────────────────────────
    cols = {}
    if samples:
        for k in samples[0]:
            try:
                cols[k] = np.asarray([s[k] for s in samples])
            except Exception:
                pass
    cols["pe"] = np.asarray(pe_series, dtype=float)
    cols["phase"] = np.asarray(phase_of)            # 구간별로 잘라 볼 수 있게
    if eq_trace:
        cols["eq_trace"] = np.asarray(eq_trace, dtype=float)

    obs, extra = [], {}
    if b.finalize is not None and verdict["status"] == OK:
        res = b.finalize(cols) or {}
        obs, extra = res.get("observables", []), res.get("extra", {})
        # 케이스가 파생 배열(g(r) 등)을 npz 로 보낼 수 있게. JSON 에 넣으면 metrics.json 이
        # 수천 개 숫자로 부풀고 postmortem 이 읽기 어려워집니다.
        cols.update({k: np.asarray(v) for k, v in res.get("arrays", {}).items()})

    np.savez_compressed(outdir / "observables.npz",
                        **{k: v for k, v in cols.items() if isinstance(v, np.ndarray)})

    m = MET.build(
        run_id=spec.run_id, case=spec.case,
        system_tags=list(b.tags),
        reference_scales={k: spec.raw["reference"][k]["symbol"]
                          for k in ("length", "energy", "time")},
        physical=dict(b.physical),
        dimensionless={g.name: g.value for g in spec.groups},
        checks=[(c, "design") for c in spec.checks],
        observables=obs,
        equilibration=MET.equilibration_series("pe", "⟨U⟩/N [kT]"),
        numerics={**{k: v for k, v in spec.numerics.items()
                     if isinstance(v, (int, float))},
                  "steps_done": n_done,
                  # ★ L4→L3 되먹임의 **입력**. `health.step_health()` 가 이 키를 읽어
                  #   L3 가 예측한 dt/τ_fast 와 대조합니다. 예전에는 가드가 이 값을
                  #   계산해 놓고 `l4` 안에만 넣어서, health 가 찾는 키
                  #   (`step_rms_sigma`)와 이름이 달라 **81런 전부 "측정 없음"** 이었습니다
                  #   — 모듈의 핵심 검사가 한 번도 돌지 않았습니다.
                  **({"step_drift_max_sigma": guard.step_disp_seen,
                      "f_max_kT_per_sigma": guard.f_max_seen} if guard.trips else {})},
        wall_seconds=wall,
        steps_per_second=(n_done / wall) if wall > 0 else None,
        extra={"l4": verdict, **({"result": extra} if extra else {})})
    MET.write(outdir, m)
    (outdir / "spec.json").write_text(json.dumps(spec.raw, indent=2, ensure_ascii=False))
    (outdir / "l4.json").write_text(json.dumps(verdict, indent=2, ensure_ascii=False,
                                               default=str))
    return verdict


def render_verdict(v: dict) -> str:
    """사람이 읽는 한 문단. `cli run` 과 케이스 스크립트가 같이 씁니다."""
    mark = {OK: "✓", DIVERGED: "✗", NAN: "✗", FROZEN: "✗"}.get(v["status"], "?")
    L = ["", "=" * 84,
         f"L4 수치 건전성 판정: {mark} {v['status'].upper()}   run_id={v.get('run_id','?')}",
         "=" * 84,
         f"  스텝     {v.get('steps_done', 0):,} / {v.get('steps_planned', 0):,}"
         f"   벽시계 {v.get('wall_seconds', 0):.1f}s",
         f"  가드     {v.get('guard_samples', 0)}회 통과"
         + (f"   최대 힘 {v['f_max_last']:.4g} kT/σ"
            f"   한 스텝 변위 {v['step_disp_last']:.3e} σ" if "f_max_last" in v else "")]
    if v.get("trip"):
        L.append(f"  ✗ 중단   [{v['trip']['phase']}] step {v['trip']['step']}: "
                 f"{v['trip']['msg']}")
    if v.get("phases"):
        L.append(f"  {'구간':<12}{'스텝':>12}{'⟨U⟩/N [kT]':>16}{'표준편차':>11}"
                 f"{'표류':>8}   정상상태 기대")
        for name, p in v["phases"].items():
            drift = f"{p['drift_sigma']:.1f}σ" if "drift_sigma" in p else "—"
            L.append(f"  {name:<12}{p.get('steps', 0):>12,}"
                     f"{p.get('pe_mean', float('nan')):>16.6g}"
                     f"{p.get('pe_std', float('nan')):>11.4g}{drift:>8}"
                     f"   {'O' if p.get('expect_steady') else '— (변하는 게 물리)'}")
    elif "pe_mean" in v:
        L.append(f"  ⟨U⟩/N   {v['pe_mean']:.6g} ± {v.get('pe_std', 0):.4g} kT"
                 + (f"   표류 {v['drift_sigma']:.1f}σ" if "drift_sigma" in v else ""))
    for w in v.get("warnings", []):
        L.append(f"  ⚠ {w}")
    L.append("  ※ L4는 수치만 봅니다 — 물리 대조는 metrics.observables 의 role 체계입니다")
    L.append("=" * 84)
    return "\n".join(L)


__all__ = ["SCHEMA", "OK", "DIVERGED", "NAN", "FROZEN", "Build", "Guard", "Diverged",
           "builder", "get_builder", "BUILDERS", "judge", "execute", "render_verdict"]
