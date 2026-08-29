"""`bdbot/` 자체 검증 — 공통 모듈이 케이스 수치를 재현하는가.

Phase 1-C에서 공통화한 것이 1-A/1-B의 실측값과 어긋나지 않는지 봅니다.
기대값은 전부 **이미 실행으로 확인된 수치**입니다 (skill `bd-physics` §6.1, §6.2).

    $PY scratch/verify_bdbot.py
"""
from __future__ import annotations

import math
import sys

import numpy as np

sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))

from bdbot import Q, checks as C, materials as M, scales as SC, sim as SIM, stats as ST

fails = []


def check(label, got, want, rtol=1e-4, unit=""):
    ok = abs(got - want) <= rtol * abs(want)
    print(f"  {'✓' if ok else '✗'} {label:<42} {got:>14.6g} {unit:<8} (기대 {want:.6g})")
    if not ok:
        fails.append(label)


print("=" * 84)
print("① 물성 유도 — 1-A/1-B 실측값 재현  (d=5 µm 실리카, 물@300 K)")
print("=" * 84)
b = M.sphere_bulk(Q(5.0, "um"), Q(300, "K"), Q(0.851, "mPa*s"), Q(2000, "kg/m^3"))
check("kT", float(b["kT"].to("J").m), 4.1419e-21, 1e-4, "J")
check("γ = 3πηd", float(b["gamma"].to("kg/s").m), 4.0102e-8, 1e-4, "kg/s")
check("D_t = kT/γ", float(b["D_t"].to("um^2/s").m), 0.1033, 1e-3, "µm²/s")
check("τ_B = d²/D_t", float(b["tau_B"].to("s").m), 242.05, 1e-4, "s")
check("τ_p = m/γ", float(b["tau_p"].to("us").m), 3.264, 1e-3, "µs")

print("\n② dt 규약 — 트랩과 소프트 페어가 같은 공식인가 (1-C의 핵심 발견)")
k_trap = Q(10, "pN/um").to("N/m")
tau_k = C.relaxation_time(b["gamma"], k_trap)
check("τ_k = γ/k  (1-A)", float(tau_k.to("ms").m), 4.010, 1e-3, "ms")
check("dt(편향 0.1%) = 2·bias·τ_k", float(C.dt_from_bias(tau_k, 1e-3).to("us").m),
      8.020, 1e-3, "µs")
check("편향 역산 bias_from_dt", C.bias_from_dt(C.dt_from_bias(tau_k, 1e-3), tau_k),
      0.1, 1e-6, "%")
# 소프트 페어도 같은 공식: 강성을 kT/d² 단위로 주면 τ_int = τ_B/U''*
# ★ 케이스별 실측값 대조는 scratch/verify_1c_equivalence.py 가 스냅샷으로 한다.
#   여기서는 **공식의 항등성**만 본다 (리포트에 표시된 반올림값을 쓰면 U''∝r⁻⁵ 라
#   3자리 반올림만으로 0.9% 어긋난다 — 실제로 한 번 걸렸다).
Upp_star = 276.4                      # 1-B A=100 에서 코드가 쓴 값
tau_int = b["tau_B"] / Upp_star
check("dt = 10⁻²·τ_int → dt/τ_B", float((C.dt_from_gate(tau_int) / b["tau_B"]).to("").m),
      1e-2 / Upp_star, 1e-12)
check("두 케이스가 같은 공식인가: γ/k ≡ γ/U''",
      float((C.relaxation_time(b["gamma"], Q(1, "N/m")) ).to("s").m),
      float((b["gamma"] / Q(1, "N/m")).to("s").m), 1e-12, "s")

print("\n③ 검사 — 하드/소프트 분류와 판정 (1-B에서 필요해진 구분)")
ck = [C.Check("모델", "관성", 8.1e-4, C.GATE),
      C.Check("적분", "해상", C.GATE, C.GATE),
      C.Check("기하", "컷오프", 0.5, 1.0),
      C.Check("통계", "관측창", 50.0, 100.0, ">=", hard=False)]
v, hf, sf, tight = C.verdict(ck)
print(f"  {'✓' if v == 'PASS (경고 1건)' else '✗'} 판정 = {v!r}  "
      f"(하드실패 {len(hf)} / 소프트실패 {len(sf)} / 여유부족 {len(tight)})")
if v != "PASS (경고 1건)":
    fails.append("verdict")
ck2 = [C.Check("기하", "컷오프", 1.5, 1.0)]
v2, _, _, _ = C.verdict(ck2)
print(f"  {'✓' if v2 == 'FAIL' else '✗'} 하드 위반이 있으면 {v2!r}")
if v2 != "FAIL":
    fails.append("verdict-hard")

print("\n④ 통계 — 자기상관 보정이 실제로 작동하는가")
rng = np.random.default_rng(0)
white = rng.normal(0, 1, 4000)
walk = np.cumsum(rng.normal(0, 1, 4000))
nw, nk = ST.n_eff(white), ST.n_eff(walk)
print(f"  {'✓' if nw > 3000 else '✗'} 백색잡음  n_eff = {nw:8.1f} / 4000  (거의 전부 독립)")
print(f"  {'✓' if nk < 500 else '✗'} 랜덤워크  n_eff = {nk:8.1f} / 4000  (강한 상관)")
if not (nw > 3000 and nk < 500):
    fails.append("n_eff")
# 드리프트: 유의성과 크기를 분리해 보고하는가
trend = np.linspace(0, 0.001, 4000) + rng.normal(0, 1e-5, 4000) + 105.5
st = ST.stationarity(trend)
print(f"  {'✓' if abs(st['drift_span_rel_pct']) < 0.01 else '✗'} 미미한 드리프트: "
      f"t={st['trend_t']:+.1f} 인데 전구간 {st['drift_span_rel_pct']:+.5f}% "
      f"→ 크기를 함께 봐야 오탐을 막는다")
if abs(st["drift_span_rel_pct"]) >= 0.01:
    fails.append("drift")
# 불편 자기상관: OU 과정의 τ 를 복원하는가
tau_true = 40.0
n = 200000
x = np.empty(n)
x[0] = 0.0
a = math.exp(-1 / tau_true)
s = math.sqrt(1 - a * a)
g = rng.normal(0, 1, n)
for i in range(1, n):
    x[i] = a * x[i - 1] + s * g[i]
ac = ST.autocorr_unbiased(x[:, None])
tau_fit = -1.0 / math.log(ac[1] / ac[0])
check("불편 자기상관에서 τ 복원 (OU)", tau_fit, tau_true, 2e-2)

print("\n⑤ sim — 함정 방어가 코드에 박혀 있는가")
ns, hs = SIM.resolve_seed(20260803)
print(f"  {'✓' if hs == 10179 else '✗'} 함정 12 시드: {ns} → HOOMD {hs} (경고값 10179)")
if hs != 10179:
    fails.append("seed")
d = SIM.minimum_image(np.array([[9.0, -9.0, 0.0]]), 10.0, dims=2)
ok = np.allclose(d, [[-1.0, 1.0, 0.0]]) and np.isfinite(d).all()
print(f"  {'✓' if ok else '✗'} 함정 1·7 최소이미지: [9,-9,0] → {d[0].tolist()} (z는 NaN 아님)")
if not ok:
    fails.append("minimum_image")

print("\n⑥ 원장 — 비로 유도하고 자릿수 폭을 볼 수 있는가")
lg = SC.ScaleLedger()
lg.lengths = {"d": Q(5, "um").to("m"), "l_k": Q(20.35, "nm").to("m"), "L": Q(150, "um").to("m")}
lg.times = {"tau_p": b["tau_p"], "tau_k": tau_k, "tau_B": b["tau_B"]}
lg.ref = SC.thermal_reference(b["d"], b["kT"], b["tau_B"])
order = [k for k, _ in lg.sorted_items(lg.times)]
print(f"  {'✓' if order == ['tau_p', 'tau_k', 'tau_B'] else '✗'} 시간척도 정렬: {order}")
if order != ["tau_p", "tau_k", "tau_B"]:
    fails.append("sorted")
check("ratio(times, tau_p, tau_k)", lg.ratio("times", "tau_p", "tau_k"), 8.139e-4, 1e-3)
check("span(times) 자릿수 폭", math.log10(lg.span("times")), 7.87, 1e-2, "자릿수")

print()
print("=" * 84)
print("✓ PASS — bdbot 공통 모듈이 케이스 실측값을 재현한다" if not fails
      else f"✗ FAIL — {len(fails)}건: {fails}")
print("=" * 84)
sys.exit(0 if not fails else 1)
