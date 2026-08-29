"""L4 — 수치 건전성 판정기.

**물리 검증기가 아닙니다.** 사용자 지시(2026-08-04, 마스터플랜 §0.2-B):

    런 후의 과학적·물리적 검증은 비중을 낮춘다 — 계마다 다르고, 보고된 계가
    아닐 수 있다. 런 후에는 **수치해석적 오류만** 본다: 발산 · NaN/Inf ·
    이상한 값으로의 수렴.

원칙 9.1도 같은 결론입니다 — 조합 결과에 기존 이론을 갖다 대지 않습니다.
그래서 여기에는 해석해도 문헌값도 없습니다. **시계열의 수치적 성질만** 봅니다.

세 부분:

  ① `Guard`          실행 중 감시. NaN/Inf/폭주를 만나면 **즉시 중단**한다.
  ② `judge(...)`     실행 후 시계열 판정. 발산·정지·붕괴.
  ③ `step_health()`  ⭐️ **L3 로 되먹임**. L3가 예측한 `dt/τ_fast` 를 L4가 측정한다.
                     어긋나면 스케일 원장에 **빠진 시간척도**가 있다는 뜻이다.

③ 이 이 모듈의 핵심입니다. 무차원 규약(σ=kT=γ=1)에서 한 스텝의 결정론적 변위는

    drift_per_step / σ  =  F* dt* / (γ* σ)  =  dt / τ_fast

이므로, **측정한 스텝 변위가 곧 `dt/τ_fast` 입니다.** L3의 예측값(원장에 있는
시간척도로 계산한 것)과 비교하면 원장의 완전성을 사후에 검사할 수 있습니다.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

# 실패 분류 (마스터플랜 §10.1 중 수치 관련만)
NUMERIC_MODES = ("NUM_NONFINITE", "NUM_DIVERGE", "NUM_FROZEN", "NUM_COLLAPSE",
                 "NUM_STEP_TOO_COARSE", "LEDGER_INCOMPLETE")

STEP_HARD = 1e-2        # dt/τ_fast 하드 한계 (bd-physics §4 와 같은 값)
LEDGER_TOL = 3.0        # L3 예측 대비 측정이 이 배수를 넘으면 원장 의심


# ════════════════════════════════════════════════════════════════════════
# ① 실행 중 감시
# ════════════════════════════════════════════════════════════════════════
class Guard:
    """`hoomd.custom.Action` 로 감싸 쓰는 런타임 감시자.

    HOOMD 를 임포트하지 않습니다 — 테스트 가능하도록 순수 함수로 두고,
    `as_action()` 에서만 hoomd 를 끌어옵니다.

        g = Guard(box_L=32.0)
        sim.operations.writers.append(g.as_action(period=10_000, state=sim.state,
                                                  thermo=thermo))
    """

    def __init__(self, box_L: float, pe_blowup: float = 1e3):
        self.box_L = float(box_L)
        self.pe_blowup = float(pe_blowup)
        self.pe0: float | None = None
        self.n_checks = 0
        self.history: list[tuple[int, float]] = []

    def check(self, timestep: int, positions: np.ndarray, pe) -> None:
        """위반이면 RuntimeError. 조용히 넘어가지 않는다 (§2 원칙 5)."""
        self.n_checks += 1
        if not np.all(np.isfinite(positions)):
            n = int((~np.isfinite(positions)).sum())
            raise RuntimeError(f"[NUM_NONFINITE] step {timestep}: 위치에 non-finite {n}개")
        far = np.abs(positions).max()
        if far > 50 * self.box_L:
            raise RuntimeError(
                f"[NUM_DIVERGE] step {timestep}: |r|max={far:.3g} > 50·L={50*self.box_L:.3g}")
        if pe is not None:
            pe = float(pe)
            if not math.isfinite(pe):
                raise RuntimeError(f"[NUM_NONFINITE] step {timestep}: PE={pe}")
            self.history.append((int(timestep), pe))
            if self.pe0 is None and pe != 0:
                self.pe0 = abs(pe)
            elif self.pe0 and abs(pe) > self.pe_blowup * self.pe0:
                raise RuntimeError(
                    f"[NUM_DIVERGE] step {timestep}: PE {pe:.3g} 가 초기 {self.pe0:.3g} 의 "
                    f"{abs(pe)/self.pe0:.0f}배 (한계 {self.pe_blowup:.0f}배)")

    def as_action(self, period: int, state, thermo=None):
        import hoomd

        outer = self

        class _A(hoomd.custom.Action):
            def act(self, timestep):
                snap = state.get_snapshot()
                pe = None if thermo is None else thermo.potential_energy
                outer.check(timestep, np.asarray(snap.particles.position), pe)

        return hoomd.write.CustomWriter(action=_A(), trigger=hoomd.trigger.Periodic(period))


# ════════════════════════════════════════════════════════════════════════
# ② 실행 후 시계열 판정
# ════════════════════════════════════════════════════════════════════════
@dataclass
class Finding:
    ok: bool
    mode: str | None
    name: str
    detail: str


@dataclass
class HealthReport:
    findings: list = field(default_factory=list)
    failure_modes: list = field(default_factory=list)
    measured: dict = field(default_factory=dict)

    @property
    def verdict(self) -> str:
        return "HEALTHY" if not self.failure_modes else "UNHEALTHY"

    def add(self, ok, mode, name, detail):
        self.findings.append(Finding(ok, None if ok else mode, name, detail))
        if not ok and mode and mode not in self.failure_modes:
            self.failure_modes.append(mode)

    def render(self) -> str:
        L = ["=" * 78, "L4 수치 건전성 판정  (물리 검증 아님 — 발산·NaN·이상수렴만)", "=" * 78]
        for f in self.findings:
            L.append(f"  {'✓' if f.ok else '✗'} {f.name:<26} {f.detail}")
        L += ["", f"VERDICT: {self.verdict}"
              + (f"   modes={self.failure_modes}" if self.failure_modes else ""), "=" * 78]
        return "\n".join(L)


def judge_series(name: str, y, rep: HealthReport, *, positive: bool = False,
                 cumulative: bool = False, t=None) -> None:
    """시계열 하나의 수치 건전성. **물리적 옳고 그름은 판정하지 않는다.**

    `cumulative=True` 는 **자라는 것이 정상인 양** (MSD·MSAD 처럼 누적된 것) 입니다.
    ★ 이걸 구분하지 않으면 오탐합니다. 실제로 물렸습니다 — `abp-rod` 두 런의 MSD 가
      "뒤/앞 1010배"로 NUM_DIVERGE 판정됐는데, 그건 그냥 확산이었습니다.
      합성 정상시계열로만 적대적 시험을 해서 못 잡았고 **실제 데이터가 잡았습니다.**

    누적량에는 대신 **초탄도 검사**를 겁니다: log-log 기울기 α (y ~ t^α).
    어떤 과감쇠/관성 동역학도 α ≤ 2 (탄도)를 넘지 못하므로, α > 2.5 는 물리 가정과
    무관하게 수치적 발산입니다.

    ★ `t` (실제 시간축)를 반드시 주세요. 인덱스를 시간으로 쓰면 **lag 가 로그 간격일 때
      α 가 엉뚱하게 나옵니다.** 이것도 실제 데이터가 잡았습니다 — abp-rod 의 MSD 는
      lag 가 로그 간격이라 인덱스 기준 α 가 2.5 를 넘어 발산으로 오판됐습니다.
      `t=None` 이면 등간격으로 가정합니다.
    """
    y = np.asarray(y, dtype=float)
    if y.size < 8:
        rep.add(True, None, f"{name}", f"표본 {y.size}개 — 판정 생략")
        return

    if not np.all(np.isfinite(y)):
        rep.add(False, "NUM_NONFINITE", f"{name} 유한성",
                f"non-finite {int((~np.isfinite(y)).sum())}/{y.size}개")
        return
    rep.add(True, None, f"{name} 유한성", f"{y.size}개 전부 유한")

    m = float(np.mean(y))
    scale = abs(m) if m else float(np.max(np.abs(y)) or 1.0)

    # 정지 — 적분기가 돌지 않았거나 완전히 얼어붙음
    rel_std = float(np.std(y)) / scale
    if rel_std < 1e-12:
        rep.add(False, "NUM_FROZEN", f"{name} 변동", f"상대 표준편차 {rel_std:.1e} — 상수")
    else:
        rep.add(True, None, f"{name} 변동", f"상대 표준편차 {rel_std:.3g}")

    q = max(2, y.size // 4)
    if cumulative:
        # 누적량 — 자라는 게 정상. 대신 탄도 상한 α ≤ 2 를 넘는지 본다.
        tt = (np.arange(1, y.size + 1, dtype=float) if t is None
              else np.asarray(t, dtype=float))
        axis = "인덱스(등간격 가정)" if t is None else "t"
        m_ = (np.abs(y) > 0) & (tt > 0) & np.isfinite(tt)
        if m_.sum() >= 8:
            alpha = float(np.polyfit(np.log(tt[m_]), np.log(np.abs(y[m_])), 1)[0])
            rep.add(alpha <= 2.5, "NUM_DIVERGE", f"{name} 성장 지수",
                    f"α = {alpha:.3f}  (y ~ t^α, 축={axis}; 탄도 상한 2, 한계 2.5)")
        else:
            rep.add(True, None, f"{name} 성장 지수", "양수 표본 부족 — 생략")
    else:
        # 정상상태량 — 뒤 1/4 이 앞 1/4 대비 폭증하면 발산
        a, b = np.abs(y[:q]).mean(), np.abs(y[-q:]).mean()
        growth = b / a if a else float("inf")
        rep.add(growth <= 1e3, "NUM_DIVERGE", f"{name} 폭증",
                f"뒤/앞 = {growth:.4g}배" + ("" if growth <= 1e3 else " (한계 1e3)"))

    # 붕괴 — 양수여야 하는 양이 0으로 수렴 (누적량은 처음이 0이라 해당 없음)
    if positive and not cumulative:
        tail = float(np.abs(y[-q:]).mean())
        if tail < 1e-12 * scale:
            rep.add(False, "NUM_COLLAPSE", f"{name} 붕괴", f"뒤 1/4 평균 {tail:.1e} ≈ 0")
        else:
            rep.add(True, None, f"{name} 붕괴", f"뒤 1/4 평균 {tail:.4g}")


# ════════════════════════════════════════════════════════════════════════
# ③ ⭐️ L3 되먹임 — 예측한 dt/τ_fast 를 측정한다
# ════════════════════════════════════════════════════════════════════════
def step_health(step_disp_rms: float | None, dt_star: float, dim: int,
                predicted_dt_over_tau: float | None, rep: HealthReport,
                *, drift_direct: float | None = None) -> None:
    """한 스텝 변위에서 `dt/τ_fast` 를 **측정**하고 L3 예측과 대조한다.

    측정 경로가 **두 개**입니다. 있으면 ⓐ를 씁니다 — ⓑ보다 엄격하게 낫습니다.

    ⓐ **힘 기반** (`drift_direct`, `bdbot.run.Guard` 가 런타임에 측정):

        drift = dt* · |F*|max / γ*  =  dt/τ_fast        (γ*=σ=1)

      열잡음이 섞이지 않으므로 **뺄 것이 없습니다.** 그리고 런 전체의 **최악값**이라
      안정성 판정에 맞습니다.

    ⓑ **위치 차분** (`step_disp_rms`, 스냅샷 두 장에서):

        Δr = (F*/γ*)·dt*  +  √(2 D* dt*)·ξ ,    D* = 1
        drift = √(max(0, ⟨Δr²⟩ − 2·dim·dt*))

      열적 성분을 제곱에서 빼야 하는데, 표류 ≪ 열잡음이면 **비슷한 두 수의 차**라
      신뢰할 수 없습니다 (그래서 아래에서 "열적 지배 — 대조 무의미" 로 빠집니다).
      `run.execute` 를 쓰지 않는 예전 런의 하위호환 경로입니다.

    `predicted_dt_over_tau` 는 L3 가 **원장의 시간척도로 계산한** 값이다.
    측정이 예측보다 크게 크면 **원장에 없는 더 빠른 시간척도가 있다** — 스케일 표가
    불완전하다는 뜻이고, 이건 L4 가 앞단(L3)에 돌려줄 수 있는 유일한 신호다.
    """
    if drift_direct is not None:
        drift = float(drift_direct)
        rep.measured["dt_over_tau_fast_measured"] = drift
        rep.measured["step_method"] = "force"
        detail = (f"dt/τ_fast = {drift:.3e}  (한계 {STEP_HARD:.0e}, "
                  f"힘 기반 dt·|F|max — 런 전체 최악값)")
    else:
        thermal2 = 2.0 * dim * dt_star
        meas2 = float(step_disp_rms) ** 2
        drift = math.sqrt(max(0.0, meas2 - thermal2))
        rep.measured["step_rms_sigma"] = float(step_disp_rms)
        rep.measured["thermal_rms_sigma"] = math.sqrt(thermal2)
        rep.measured["dt_over_tau_fast_measured"] = drift
        rep.measured["step_method"] = "position"
        detail = (f"dt/τ_fast = {drift:.3e}  (한계 {STEP_HARD:.0e}, "
                  f"열적분 {math.sqrt(thermal2):.3e} 제외)")

    rep.add(drift <= STEP_HARD, "NUM_STEP_TOO_COARSE", "스텝 변위 (측정)", detail)

    if predicted_dt_over_tau is None:
        rep.add(True, None, "L3 원장 대조", "L3 예측값 없음 — 대조 생략")
        return
    rep.measured["dt_over_tau_fast_predicted"] = float(predicted_dt_over_tau)
    if drift <= 0 or predicted_dt_over_tau <= 0:
        rep.add(True, None, "L3 원장 대조", "표류가 0 — 열적 지배, 대조 무의미")
        return
    ratio = drift / predicted_dt_over_tau
    rep.measured["ledger_ratio"] = ratio
    rep.add(ratio <= LEDGER_TOL, "LEDGER_INCOMPLETE", "L3 원장 완전성",
            f"측정/예측 = {ratio:.2f}× "
            + (f"— 원장에 없는 더 빠른 시간척도 의심 (한계 {LEDGER_TOL:.0f}×)"
               if ratio > LEDGER_TOL else "— 원장이 최속 척도를 담고 있음"))

    # ── 반대 방향: 설계가 과보수적이면 **비용**이다 (실패가 아님, 정보) ──────
    # 원장 검사는 ratio ≫ 1 (빠진 척도)만 봅니다. ratio ≪ 1 은 건전하지만
    # "dt 를 필요 이상으로 작게 잡았다" = 벽시계를 그만큼 더 쓴 것입니다.
    # 이 프로젝트에서 비용은 반복해서 문제였습니다 (스윕 25일 → 1.16일).
    # ⚠️ 단정하지 않습니다: ① 가드는 GUARD_EVERY 스텝마다만 보므로 표본 사이의 최댓값을
    #    놓칠 수 있고 ② L3 의 r_min 은 **설계 최악 접근거리**여서 이 런이 거기까지
    #    가지 않았을 수 있습니다. 그래서 "확인해 볼 여지"로만 보고합니다.
    if ratio < 1.0 / LEDGER_TOL and ratio > 0:
        rep.measured["dt_headroom"] = 1.0 / ratio
        rep.add(True, None, "dt 여유 (비용)",
                f"측정이 예측의 {ratio:.3f}× — 실제 만난 최대 힘 기준으로는 dt 를 "
                f"약 {1.0/ratio:.0f}배까지 키울 여지가 있습니다 (벽시계 ÷{1.0/ratio:.0f}). "
                f"단정 금지: 가드 표본 사이의 최댓값을 놓쳤거나, L3 의 설계 최악 "
                f"접근거리까지 이 런이 가지 않았을 수 있습니다 — 수렴 확인으로 검증하세요")


def measure_step_displacement(positions_t0, positions_t1, L: float, dim: int) -> float:
    """연속한 두 스냅샷에서 한 스텝 변위의 rms (σ 단위). 최소 이미지 적용 (함정 1·7)."""
    d = np.asarray(positions_t1, float) - np.asarray(positions_t0, float)
    d[:, :2] -= L * np.round(d[:, :2] / L)
    if dim == 3:
        d[:, 2] -= L * np.round(d[:, 2] / L)
    return float(np.sqrt((d[:, :dim] ** 2).sum(axis=1).mean()))


# ════════════════════════════════════════════════════════════════════════
# 스펙에서 L3 예측 뽑기
# ════════════════════════════════════════════════════════════════════════
def predicted_dt_over_tau(spec) -> float | None:
    """`LoadedSpec` 의 적분 해상 검사에서 L3가 예측한 `dt/τ_fast` 를 꺼낸다.

    검사 이름은 케이스마다 다르므로 `kind == '적분'` 중 **가장 큰 값**을 씁니다
    (가장 빠른 시간척도가 가장 큰 비를 만든다).
    """
    vals = [c.value for c in getattr(spec, "checks", [])
            if getattr(c, "kind", "") == "적분" and isinstance(c.value, (int, float))]
    return max(vals) if vals else None


def gate(spec) -> list[str]:
    """실행 **전** 게이트. **막아야 할** 이유들을 돌려준다 (빈 리스트 = 통과).

    L4 는 스펙만 읽습니다 — 케이스 코드를 임포트하지 않습니다 (마스터플랜 L2↔L4 계약).

    ⚠️ **소프트 경고는 막지 않습니다.** 예전에는 `verdict != "PASS"` 로 판정해서
       `"PASS (경고 3건)"` 을 거부했습니다 — 실측 83개 스펙 중 **80개가 거짓 거부**였고
       그중 진짜 하드 실패는 **0개**였습니다. bd-physics §4 는 통계·유한크기를 ⚠ 경고로
       규정합니다(❌ 아님). `run.execute` 는 `startswith("FAIL")` 로 옳게 보고 있었고,
       둘이 어긋난 채로 아무도 눈치채지 못한 이유는 **`execute` 가 `gate()` 를 부르지
       않아서** 입니다 — 배선되지 않은 검사기는 틀려도 드러나지 않습니다.
       막는 것은 ① 해시 불일치 ② 하드 검사 실패(FAIL) ③ L3 무결성 **오류**뿐입니다.
    """
    problems = []
    ok, want = spec.verify_hash()
    if not ok:
        problems.append(f"run_id 불일치: 저장 {spec.run_id} vs 계산 {want} "
                        f"— 스펙을 손으로 고쳤습니까? (규칙 2)")
    if spec.verdict.startswith("FAIL"):
        bad = [c.name for c in spec.checks if getattr(c, "hard", True) and not _ok(c)]
        problems.append(f"L3 verdict={spec.verdict}"
                        + (f" — 하드 검사 실패: {bad}" if bad else ""))
    errs = [i for i in spec.raw.get("l3_issues", []) if i.get("level") == "error"]
    if errs:
        problems.append(f"L3 무결성 오류 {len(errs)}건: {errs[:2]}")
    return problems


def gate_notes(spec) -> list[str]:
    """게이트를 **막지는 않지만** 사람이 봐야 하는 것들. 조용히 넘기지 않기 위해서."""
    notes = []
    soft = [c for c in spec.checks if not getattr(c, "hard", True) and not _ok(c)]
    for c in soft:
        notes.append(f"소프트 경고 [{c.kind}] {c.name.strip()} = {c.value:.3g} "
                     f"(기준 {c.limit:g}) — 막지 않지만 통계·유한크기 한계입니다")
    tight = [c for c in spec.checks
             if _ok(c) and getattr(c, "hard", True) and _margin(c) < 5.0]
    for c in tight:
        notes.append(f"여유 부족 [{c.kind}] {c.name.strip()} — 한계까지 {_margin(c):.1f}배뿐")
    warn = [i for i in spec.raw.get("l3_issues", []) if i.get("level") != "error"]
    for i in warn:
        notes.append(f"L3 {i.get('level')} [{i.get('where')}] {i.get('msg')}")
    return notes


def _margin(c) -> float:
    if getattr(c, "op", "<=") == "<=":
        return c.limit / c.value if c.value else float("inf")
    return c.value / c.limit if c.limit else float("inf")


def _ok(c) -> bool:
    return c.value <= c.limit if getattr(c, "op", "<=") == "<=" else c.value >= c.limit


__all__ = ["Guard", "HealthReport", "judge_series", "step_health",
           "measure_step_displacement", "predicted_dt_over_tau", "gate", "gate_notes",
           "NUMERIC_MODES", "STEP_HARD", "LEDGER_TOL"]
