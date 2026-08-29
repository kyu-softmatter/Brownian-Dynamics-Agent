"""L4 수치 건전성 판정기 적대적 검사.

CLAUDE.md 관행: "검사기를 만들면 일부러 망가뜨려 시험한다 — 조용히 통과와 검사를
안 함은 다르다." 각 실패 모드를 인위적으로 만들어 **정확히 그 모드로** 잡히는지 본다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_health.py
"""
from __future__ import annotations

import json
import math
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
from bdbot import health as H          # noqa: E402
from bdbot import nondim as ND         # noqa: E402

ROOT = Path(__file__).resolve().parent.parent
rng = np.random.default_rng(4)
PASS, FAIL = [], []


def check(name, cond, detail=""):
    (PASS if cond else FAIL).append(name)
    print(f"  {'✓' if cond else '✗'} {name}" + (f"   {detail}" if detail else ""))


def modes_of(series, **kw):
    r = H.HealthReport()
    H.judge_series("y", series, r, **kw)
    return r.failure_modes, r


print("=" * 78)
print("① 시계열 판정 — 각 실패 모드를 인위적으로 만든다")
print("=" * 78)

n = 400
healthy = 1.0 + 0.05 * rng.normal(size=n)
m, _ = modes_of(healthy)
check("정상 시계열은 통과", m == [], f"modes={m}")

bad = healthy.copy(); bad[137] = np.nan
m, _ = modes_of(bad)
check("NaN → NUM_NONFINITE", m == ["NUM_NONFINITE"], f"modes={m}")

bad = healthy.copy(); bad[200] = np.inf
m, _ = modes_of(bad)
check("Inf → NUM_NONFINITE", m == ["NUM_NONFINITE"], f"modes={m}")

explode = np.exp(np.linspace(0, 25, n)) * (1 + 0.01 * rng.normal(size=n))
m, _ = modes_of(explode)
check("지수 폭증 → NUM_DIVERGE", "NUM_DIVERGE" in m, f"modes={m}")

const = np.full(n, 3.14159)
m, _ = modes_of(const)
check("상수 → NUM_FROZEN", "NUM_FROZEN" in m, f"modes={m}")

collapse = np.abs(healthy) * np.exp(-np.linspace(0, 60, n))
m, _ = modes_of(collapse, positive=True)
check("0으로 붕괴 → NUM_COLLAPSE", "NUM_COLLAPSE" in m, f"modes={m}")

m, _ = modes_of(collapse, positive=False)
check("positive=False 면 붕괴 판정 안 함", "NUM_COLLAPSE" not in m, f"modes={m}")

m, _ = modes_of(healthy[:5])
check("표본 부족은 판정 생략(통과)", m == [], f"modes={m}")

# ★ 실제 데이터가 잡아낸 오탐 — 회귀 테스트는 **깨졌던 진짜 데이터**로 한다.
#   합성 선형 MSD 는 뒤/앞이 7배뿐이라 오탐을 재현하지 못했다. 실제 abp-rod 의 MSD 는
#   lag 가 로그 간격이라 1010배가 나왔고, 그래서 NUM_DIVERGE 로 오판됐다.
# ★ 픽스처를 **합성으로 자족화**했습니다 (2026-08-05). 예전에는 `runs/abp-rod*` 의 실제
#   npz 를 골라 썼는데, 그 런이 사라지자(런 디렉토리는 재실행·정리로 바뀝니다) 뒤/앞 비가
#   1010배 → 343배로 떨어져 **오탐을 재현할 수 없어 테스트가 깨졌습니다.**
#   회귀 테스트가 특정 런의 존재에 의존하면 회귀를 못 지킵니다.
#   원인은 "실제 데이터"가 아니라 **lag 가 로그 간격**인 것이었으므로, 그것만 합성하면
#   충분합니다 — 로그 lag + 확산(α=1) MSD.
t_log = np.logspace(0, 4.5, 40)                     # 로그 간격 lag (실제 npz 와 같은 구조)
msd_log = 4.0 * t_log * (1 + 0.02 * rng.normal(size=t_log.size))
q_ = max(2, msd_log.size // 4)
ratio_log = np.abs(msd_log[-q_:]).mean() / np.abs(msd_log[:q_]).mean()
check("로그 lag MSD 는 뒤/앞 비가 1e3 을 넘는다 (오탐 조건 재현)", ratio_log > 1e3,
      f"뒤/앞 = {ratio_log:.3g}배  (선형 lag 이면 7배뿐이라 재현 안 됨)")

m, _ = modes_of(msd_log, positive=True)
check("MSD 를 정상상태로 보면 오탐한다 (회귀)", "NUM_DIVERGE" in m,
      f"뒤/앞 = {ratio_log:.3g}배 → {m}")

m, _ = modes_of(msd_log, positive=True, cumulative=True)
check("시간축 없이 cumulative 만으로는 부족 (로그 lag)", "NUM_DIVERGE" in m,
      f"인덱스 기준 α 가 상한 초과 → {m}")

m, r_ = modes_of(msd_log, positive=True, cumulative=True, t=t_log)
alpha = [f.detail for f in r_.findings if "성장 지수" in f.name]
check("실제 시간축을 주면 통과", "NUM_DIVERGE" not in m, f"{alpha[0] if alpha else m}")

# 실제 런이 있으면 **추가로** 대조한다 (있으면 좋고, 없어도 위 합성 검사가 회귀를 지킨다)
real = sorted((ROOT / "runs").glob("abp-rod*/observables.npz"))
real = [f for f in real if "msd" in np.load(f).files]
if real:
    with np.load(real[0]) as z:
        msd_real, t_real = np.asarray(z["msd"], float), np.asarray(z["t"], float)
    m, r_ = modes_of(msd_real, positive=True, cumulative=True, t=t_real)
    check("실제 abp-rod MSD 도 시간축을 주면 통과", "NUM_DIVERGE" not in m,
          f"{real[0].parent.name[:40]} → {m or 'modes 없음'}")
else:
    print("      (실제 abp-rod 런 없음 — 합성 검사만으로 회귀는 지켜집니다)")

t = np.arange(1, n + 1, dtype=float)
msd_diffusive = 4.0 * t * (1 + 0.02 * rng.normal(size=n))       # α=1
m, _ = modes_of(msd_diffusive, positive=True, cumulative=True)
check("확산 MSD (α≈1) 통과", m == [], f"modes={m}")

msd_ballistic = 9.0 * t ** 2 * (1 + 0.02 * rng.normal(size=n))  # α=2
m, _ = modes_of(msd_ballistic, positive=True, cumulative=True)
check("탄도 MSD (α=2) 도 통과", m == [], f"modes={m}")

msd_blowup = t ** 3.4                                            # α=3.4 — 물리 불가
m, _ = modes_of(msd_blowup, positive=True, cumulative=True)
check("초탄도 α=3.4 → NUM_DIVERGE", "NUM_DIVERGE" in m, f"modes={m}")

print()
print("=" * 78)
print("② 스텝 변위 → dt/τ_fast 측정 + L3 원장 대조")
print("=" * 78)

dt_star, dim = 1e-4, 2
thermal = math.sqrt(2 * dim * dt_star)

r = H.HealthReport()
H.step_health(thermal, dt_star, dim, predicted_dt_over_tau=2e-3, rep=r)
check("순수 열적 스텝 → 표류 0, 통과", r.verdict == "HEALTHY",
      f"drift={r.measured['dt_over_tau_fast_measured']:.2e}")

drift_true = 5e-3
r = H.HealthReport()
H.step_health(math.hypot(thermal, drift_true), dt_star, dim, 4e-3, r)
got = r.measured["dt_over_tau_fast_measured"]
check("표류를 열적분과 분리해 복원", abs(got - drift_true) / drift_true < 0.02,
      f"넣은 값 {drift_true:.1e} → 복원 {got:.3e}")

r = H.HealthReport()
H.step_health(math.hypot(thermal, 5e-2), dt_star, dim, 4e-2, r)
check("스텝 과대 → NUM_STEP_TOO_COARSE", "NUM_STEP_TOO_COARSE" in r.failure_modes,
      f"modes={r.failure_modes}")

r = H.HealthReport()
H.step_health(math.hypot(thermal, 8e-3), dt_star, dim, predicted_dt_over_tau=1e-4, rep=r)
check("측정 ≫ L3 예측 → LEDGER_INCOMPLETE", "LEDGER_INCOMPLETE" in r.failure_modes,
      f"측정/예측 = {r.measured['ledger_ratio']:.0f}×  ← 원장에 빠진 시간척도")

r = H.HealthReport()
H.step_health(math.hypot(thermal, 4e-3), dt_star, dim, predicted_dt_over_tau=4e-3, rep=r)
check("측정 ≈ L3 예측 → 원장 완전", "LEDGER_INCOMPLETE" not in r.failure_modes,
      f"측정/예측 = {r.measured['ledger_ratio']:.2f}×")

# 최소 이미지가 스텝 변위에 적용되는가 (함정 1·7)
L = 20.0
p0 = np.array([[9.9, 0.0, 0.0]])
p1 = np.array([[-9.9, 0.0, 0.0]])           # 경계를 넘어 wrap — 실제 이동은 0.2
d = H.measure_step_displacement(p0, p1, L, 2)
check("스텝 변위에 최소 이미지 적용", abs(d - 0.2) < 1e-9, f"{d:.4f} (참값 0.2, wrap 무시 시 19.8)")

print()
print("=" * 78)
print("②b 힘 기반 측정 경로 (drift_direct) — 81런 전부 '측정 없음' 이던 구멍")
print("=" * 78)
# `run.Guard` 가 dt·|F|max 를 측정하는데 health 는 `step_rms_sigma` 를 찾고 있어서
# 이름이 어긋났고, 그 결과 이 모듈의 **핵심 검사가 한 번도 돌지 않았습니다**.
r = H.HealthReport()
H.step_health(None, dt_star, dim, predicted_dt_over_tau=2e-3, rep=r, drift_direct=1.8e-3)
check("힘 기반 경로가 열적분을 빼지 않는다",
      abs(r.measured["dt_over_tau_fast_measured"] - 1.8e-3) < 1e-15,
      f"넣은 값 1.8e-03 → 그대로 {r.measured['dt_over_tau_fast_measured']:.3e} "
      f"(위치 차분이면 열적분 {thermal:.2e} 를 빼야 해서 0 이 됨)")
check("측정 방법이 기록된다", r.measured.get("step_method") == "force",
      f"step_method={r.measured.get('step_method')}")
check("힘 기반 정상값은 HEALTHY", r.verdict == "HEALTHY", f"modes={r.failure_modes}")

r = H.HealthReport()
H.step_health(None, dt_star, dim, 2e-3, r, drift_direct=5e-2)
check("힘 기반 과대 스텝 → NUM_STEP_TOO_COARSE",
      "NUM_STEP_TOO_COARSE" in r.failure_modes, f"modes={r.failure_modes}")

r = H.HealthReport()
H.step_health(None, dt_star, dim, 1e-4, r, drift_direct=8e-3)
check("힘 기반에서도 원장 불완전을 잡는다", "LEDGER_INCOMPLETE" in r.failure_modes,
      f"측정/예측 = {r.measured['ledger_ratio']:.0f}×")

# ★ 힘 기반이 위치 차분보다 나은 이유 — **통계 오차** 때문입니다 (부동소수 문제가 아님).
#   정확한 산술에서는 위치 차분도 표류를 복원합니다. 문제는 유한 표본의 rms 자체가
#   ~rms/√(2N) 만큼 흔들리고, 표류 ≪ 열잡음이면 그 흔들림이 표류를 삼킨다는 것입니다.
#   실제 입자 앙상블로 시험합니다.
rng = np.random.default_rng(20260805)
N_p = 1000
drift_small = 1e-4                                  # 열잡음(2e-2)의 0.5%
sig = math.sqrt(2.0 * dt_star)                      # 성분당 열적 표준편차
noise = rng.normal(0.0, sig, size=(N_p, dim))
step = noise + np.array([drift_small] + [0.0] * (dim - 1))   # x 방향 결정론적 표류
rms_sampled = float(np.sqrt((step ** 2).sum(axis=1).mean()))

r_pos = H.HealthReport()
H.step_health(rms_sampled, dt_star, dim, drift_small, r_pos)
r_frc = H.HealthReport()
H.step_health(None, dt_star, dim, drift_small, r_frc, drift_direct=drift_small)
err_pos = abs(r_pos.measured["dt_over_tau_fast_measured"] - drift_small) / drift_small
err_frc = abs(r_frc.measured["dt_over_tau_fast_measured"] - drift_small) / drift_small
# 위치차분이 복원한 표류가 참값의 10% 미만이면 "표류를 잃었다"고 봅니다.
# 실측에서는 표본 rms 가 열적 기대값보다 **작게** 나와 뺄셈이 음수→0 으로 잘립니다.
pos_recovered = r_pos.measured["dt_over_tau_fast_measured"]
check("열적 지배 영역: 유한 표본 위치차분은 표류를 잃는다 (힘 기반은 정확)",
      err_frc < 1e-12 and pos_recovered < 0.1 * drift_small,
      f"표류 {drift_small:.0e} = 열잡음 {thermal:.1e} 의 {100*drift_small/thermal:.1f}%, N={N_p}: "
      f"위치차분 복원 {pos_recovered:.2e} (참값의 {100*pos_recovered/drift_small:.0f}% — "
      f"표본 rms 가 기대값보다 작아 0 으로 잘림) vs 힘기반 오차 {100*err_frc:.0e}%")

print()
print("=" * 78)
print("③ 런타임 Guard — 즉시 중단하는가")
print("=" * 78)

g = H.Guard(box_L=32.0)
try:
    g.check(0, np.zeros((10, 3)), 1.0); ok = True
except RuntimeError:
    ok = False
check("정상 상태는 통과", ok)

pos = np.zeros((10, 3)); pos[3, 1] = np.nan
try:
    H.Guard(32.0).check(100, pos, 1.0); ok = False
except RuntimeError as e:
    ok = "NUM_NONFINITE" in str(e)
check("위치 NaN → 즉시 중단", ok)

pos = np.zeros((10, 3)); pos[0, 0] = 1e6
try:
    H.Guard(32.0).check(100, pos, 1.0); ok = False
except RuntimeError as e:
    ok = "NUM_DIVERGE" in str(e)
check("좌표 폭주 → 즉시 중단", ok)

g = H.Guard(32.0, pe_blowup=100)
g.check(0, np.zeros((5, 3)), 1.0)
try:
    g.check(10, np.zeros((5, 3)), 1e5); ok = False
except RuntimeError as e:
    ok = "NUM_DIVERGE" in str(e)
check("에너지 폭증 → 즉시 중단", ok)

try:
    H.Guard(32.0).check(0, np.zeros((5, 3)), float("nan")); ok = False
except RuntimeError as e:
    ok = "NUM_NONFINITE" in str(e)
check("PE NaN → 즉시 중단", ok)

g = H.Guard(32.0, pe_blowup=100)
g.check(0, np.zeros((5, 3)), 1.0)
try:
    g.check(10, np.zeros((5, 3)), 50.0); ok = True     # 50배 — 한계 이내
except RuntimeError:
    ok = False
check("한계 이내 증가는 통과 (과민하지 않음)", ok)

print()
print("=" * 78)
print("④ 실행 전 게이트 — 실제 스펙으로")
print("=" * 78)

specs = sorted((ROOT / "specs").glob("*.json"))
check("specs/ 에 스펙이 있다", len(specs) > 0, f"{len(specs)}건")

spec = ND.load(specs[0])
probs = H.gate(spec)
check(f"진짜 스펙은 게이트 통과 ({spec.run_id[:34]})", probs == [], f"{probs}")

# 손으로 고친 스펙 — run_id 해시가 어긋나야 한다 (규칙 2)
raw = json.loads(specs[0].read_text())
key = next(iter(raw["params"]))
raw["params"][key] = raw["params"][key] * 1.5 if isinstance(raw["params"][key], (int, float)) else "X"
tmp = Path(tempfile.mkdtemp()) / "tampered.json"
tmp.write_text(json.dumps(raw))
probs = H.gate(ND.load(tmp))
check("손으로 고친 스펙 → 게이트 거부", any("run_id 불일치" in p for p in probs),
      f"{probs[:1]}")

# L3 예측 추출
p = H.predicted_dt_over_tau(spec)
check("L3 예측 dt/τ 를 스펙에서 뽑는다", p is not None and p > 0, f"dt/τ_fast(L3) = {p}")

n_ok = 0
for sp in specs[:12]:
    s = ND.load(sp)
    if H.predicted_dt_over_tau(s) is not None:
        n_ok += 1
check("여러 케이스에서 예측 추출 가능", n_ok >= 8, f"{n_ok}/12 건")

print()
print("=" * 78)
print("④b 게이트가 소프트 경고를 막지 않는가 — 83개 중 80개를 거짓 거부하던 버그")
print("=" * 78)
# 예전 게이트는 `verdict != "PASS"` 였습니다. `checks.verdict()` 는 소프트 검사가 걸리면
# "PASS (경고 N건)" 을 돌려주므로 **경고만 있는 스펙이 전부 거부**됐습니다.
# bd-physics §4 는 통계·유한크기를 ⚠(❌ 아님)로 규정하고, `run.execute` 는
# startswith("FAIL") 로 옳게 봤습니다 — 둘이 어긋난 채 아무도 몰랐던 이유는
# `execute` 가 `gate()` 를 부르지 않아서입니다.
blocked_wrongly = []
n_soft = n_pass = 0
for sp in specs:
    s = ND.load(sp)
    probs = H.gate(s)
    hard = [c for c in s.checks if getattr(c, "hard", True) and not H._ok(c)]
    errs = [i for i in s.raw.get("l3_issues", []) if i.get("level") == "error"]
    if s.verdict.startswith("PASS (경고"):
        n_soft += 1
    if not probs:
        n_pass += 1
    if probs and not hard and not errs and not s.verdict.startswith("FAIL"):
        blocked_wrongly.append((s.run_id, s.verdict))
check("소프트 경고만인 스펙을 막지 않는다", not blocked_wrongly,
      f"스펙 {len(specs)}개 중 경고만 {n_soft}개 · 게이트 통과 {n_pass}개 · "
      f"거짓 거부 {len(blocked_wrongly)}개"
      + (f"  예: {blocked_wrongly[0]}" if blocked_wrongly else ""))

soft_specs = [sp for sp in specs if ND.load(sp).verdict.startswith("PASS (경고")]
if soft_specs:
    s = ND.load(soft_specs[0])
    check("경고는 막지 않되 gate_notes 로 반드시 드러낸다",
          len(H.gate_notes(s)) > 0,
          f"{s.run_id[:34]}: {len(H.gate_notes(s))}건 — {H.gate_notes(s)[0][:60]}")
else:
    check("경고 있는 스펙이 있어야 이 회귀를 시험할 수 있다", False, "없음")

# 하드 실패는 여전히 막아야 한다 (반대 방향 오작동 방지)
raw = json.loads(specs[0].read_text())
raw["verdict"] = "FAIL"
tmp2 = Path(tempfile.mkdtemp()) / "failverdict.json"
tmp2.write_text(json.dumps(raw))
probs = H.gate(ND.load(tmp2))
check("verdict=FAIL 은 여전히 막는다", any("verdict=FAIL" in p for p in probs), f"{probs}")

# L3 무결성 **오류**는 막고, warn/info 는 막지 않는다
for lvl, should_block in (("error", True), ("warn", False), ("info", False)):
    raw = json.loads(specs[0].read_text())
    raw["l3_issues"] = [{"level": lvl, "where": "test", "msg": "주입"}]
    t3 = Path(tempfile.mkdtemp()) / f"iss_{lvl}.json"
    t3.write_text(json.dumps(raw))
    probs = [p for p in H.gate(ND.load(t3)) if "무결성" in p]
    check(f"l3_issues level={lvl} → {'막는다' if should_block else '막지 않는다'}",
          bool(probs) == should_block, f"{probs}")

print()
print("=" * 78)
print(f"{len(PASS)}/{len(PASS) + len(FAIL)} PASS")
if FAIL:
    for f in FAIL:
        print(f"  ✗ {f}")
print("=" * 78)
sys.exit(1 if FAIL else 0)
