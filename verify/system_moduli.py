"""시스템 수준 유변학 — 개별 비드 피팅이 아니라 **계 전체의 에너지 수지**로.

★ 왜 이 관점인가
  지금까지 K*(ω) 는 구동 비드 하나의 위상자에서 나왔다 (`lockin.k_star`).
  여기서는 계 전체를 하나의 점탄성체로 보고, **한 주기당 들어간 일과 저장된 에너지**만
  본다. 비드가 몇 개인지, 어느 비드가 어떤 위상인지 몰라도 된다.

  ① 총 소산 (모델 무관) — 정상상태에서 구동 트랩이 한 일이 전부 유체로 나간다
        W_cycle = ∮ F·dy,   F = k_t(y_ghost − y_bead)
        K″_total ≡ W_cycle / (π|ŷ|²)
  ② 분해 — 이 소산은 두 곳에서 난다
        K″_total = K″_chain + ωγ
        └ 사슬이 하는 소산  └ **구동 비드 자신의 용매 항력** (사슬이 없어도 나는 몫)
     `lockin.k_star` 는 정의상 −iωγ 를 빼서 K″_chain 만 낸다 (그 도크스트링 참조).
     → 두 경로가 **정확히 ωγ 만큼 달라야 한다.** 이게 자유 파라미터 없는 교차검증이다.
     실측: JKR −0.71%, DLVO −0.11% 로 일치.
  ③ 총 저장 — 계 전체 퍼텐셜 에너지 U(t) 의 **2ω** 성분 (U ∝ y² 이므로 배주파)
        U_osc 진폭 = ¼ K′_sys |ŷ|²   →   K′_sys = 4|Û_2ω|/|ŷ|²
     ⚠ K′_sys ≠ K′_락인 이다. 락인은 **구동점**에서 느끼는 강성이고, K′_sys 는
       사슬 전체(굽힘·결합신축·트랩)에 저장된 에너지다. 차이가 곧 "구동점 밖에
       저장된 몫"이라 그 자체가 정보다.

★★ 시스템 관점의 결론 한 줄:
   DLVO 사슬은 **유변학적으로 보이지 않는다** — 계 전체의 소산이 구동 비드 하나가
   자유롭게 끌려다닐 때와 같다. JKR 은 사슬이 소산의 76% 를 담당한다.

    $PY scratch/system_moduli.py
"""
import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
kT = 1.380649e-23 * 300.0
d = 1.47e-6
KTD2 = kT / d ** 2                      # kT/d² → N/m


def analyze(pat, label):
    rows = []
    for dd in sorted(glob.glob(str(ROOT / pat))):
        if not (Path(dd) / "metrics.json").exists():
            continue
        z = np.load(Path(dd) / "observables.npz", allow_pickle=True)
        s = json.loads((Path(dd) / "spec.json").read_text())
        om = float(s["params"]["omega_star"])
        n = int(s["params"]["n_beads"])
        k_t = float(s["params"]["k_t_star"])
        pos = str(s["params"].get("drive_mode", "trap")) == "position"
        t, pe, yb, yg = (np.asarray(z[k], float) for k in ("t", "pe", "y_bead", "y_ghost"))
        U = pe * n                                       # pe 는 입자당 → 계 전체
        m = {x["name"]: x["measured"] for x in
             json.loads((Path(dd) / "metrics.json").read_text())["observables"]}

        Ay = abs(2 * np.mean(yb * np.exp(-1j * om * t)))          # 구동 비드 진폭 |ŷ|

        # ① 총 소산 — 트랩이 구동 비드에 한 일. position 모드는 트랩이 없어 정의 불가.
        if pos:
            K2_tot = np.nan
        else:
            F = k_t * (yg - yb)
            Pmean = np.trapezoid(F * np.gradient(yb, t), t) / (t[-1] - t[0])   # 평균 일률
            K2_tot = Pmean * (2 * np.pi / om) / (np.pi * Ay ** 2)

        # ③ 총 저장 — U 의 2ω 성분.  U_osc = −¼K′|ŷ|²cos2ωt  →  |Û_2ω| = ¼K′|ŷ|²
        Uhat2 = 2 * np.mean((U - U.mean()) * np.exp(-2j * om * t))
        Kp_sys = 4 * abs(Uhat2) / Ay ** 2

        rows.append(dict(seed=s["numerics"]["seed"], Ay=Ay, om=om, n=n, pos=pos,
                         Kp_lock=m.get("K_prime"), K2_lock=m.get("K_doubleprime"),
                         K2_tot=K2_tot, Kp_sys=Kp_sys))
    if not rows:
        print(f"  [{label}] 런 없음")
        return None
    g = lambda k: (np.nanmean([r[k] for r in rows]),
                   np.nanstd([r[k] for r in rows], ddof=1) / np.sqrt(len(rows)))
    om = rows[0]["om"]
    print(f"\n  {label}   시드 {len(rows)}개   |ŷ| = {g('Ay')[0]:.5f} d")
    k2t, k2te = g("K2_tot")
    if np.isfinite(k2t):
        k2l, k2le = g("K2_lock")
        drag = om                                        # ωγ, 무차원계에서 γ*=1
        print(f"    ① 총 소산      K″_total = {k2t:>10.4g} ± {k2te:<8.3g}")
        print(f"    ② 분해         사슬 {k2l:>10.4g}  +  구동비드 항력 ωγ = {drag:.4g}"
              f"  =  {k2l + drag:>10.4g}   (총계와 {100*(k2t-k2l-drag)/k2t:+.2f}%)")
        print(f"       → 사슬이 담당하는 소산 비율 = {100*k2l/k2t:>6.1f}%")
    kps, kpse = g("Kp_sys")
    kpl, _ = g("Kp_lock")
    print(f"    ③ 총 저장      K′_sys = {kps:>10.4g} ± {kpse:<8.3g}"
          + (f"   (구동점 K′_락인 {kpl:.4g} 의 {kps/kpl:.2f}배)" if kpl and abs(kpl) > 1 else
             f"   (구동점 K′_락인 = {kpl:.3g} ≈ 0 인데도 저장은 있다 — 결합 신축)"))
    return rows


print("=" * 98)
print("시스템 수준 유변학 — 계 전체의 에너지 수지 (개별 비드 피팅 없음)")
print("=" * 98)
P = "runs/chain-bend-2d-dlvo__n9-w3000-a1470"
rj = analyze(f"{P}-jkr-kt100__*", "DLVO+JKR  (k_t×100)")
rd = analyze(f"{P}-kt100__*", "DLVO-only (k_t×100)")

if rj and rd:
    om = rj[0]["om"]
    tj = np.nanmean([r["K2_tot"] for r in rj])
    td = np.nanmean([r["K2_tot"] for r in rd])
    print()
    print("=" * 98)
    print("★★ 시스템 관점의 결론")
    print("=" * 98)
    print(f"  구동 비드가 **혼자** 있을 때의 소산 (사슬 없음) : ωγ = {om:.4g}")
    print(f"  DLVO 사슬을 붙였을 때 계 전체 소산              : {td:.4g}"
          f"   → {td/om:.3f}배  (사슬 기여 {100*(td-om)/om:+.1f}%)")
    print(f"  JKR  사슬을 붙였을 때 계 전체 소산              : {tj:.4g}"
          f"   → {tj/om:.3f}배  (사슬 기여 {100*(tj-om)/om:+.1f}%)")
    print()
    print("  ⟹ DLVO 사슬은 **유변학적으로 보이지 않는다**. 계를 통째로 보면 소산이")
    print("     비드 하나를 그냥 물속에서 흔드는 것과 같다. 사슬이 붙어 있으나 마나다.")
    print("     JKR 은 사슬이 소산의 대부분을 담당해 계의 응답을 지배한다.")
