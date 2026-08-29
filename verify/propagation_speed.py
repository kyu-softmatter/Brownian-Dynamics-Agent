"""응력 전파 속도 — 비드별 위상 지연의 공간 기울기에서.

정상상태 진동에서 구동 비드로부터 s 결합 떨어진 비드의 위상이 φ(s) 만큼 뒤처지면
파수 k = −dφ/ds [rad/결합] 이고, 위상속도는

    v_phase = ω / k = ω·ℓ / |dφ/ds|        [m/s]      (ℓ = 결합길이)

★ 무엇을 재는 것인가 — 이 계는 **과감쇠**(BD, 관성 없음)라 음파 같은 탄성파가 아니다.
  굽힘 변형이 점성 저항을 뚫고 사슬을 타고 **확산적으로 스며드는** 속도다.
  과감쇠 굽힘의 분산관계는 ω ~ (κ_θℓ/γ)k⁴ 이라 v_phase 가 ω 에 의존한다 (비분산 아님).
  → **한 ω 에서의 값**이며, 물질 상수가 아니다. 반드시 ω 를 같이 적을 것.

★ 진폭 감쇠길이 λ 도 같이 내지만 **"침투깊이"로 읽으면 안 된다.**
  3점 굽힘에서 진폭이 중심에서 멀어지며 줄어드는 것은 대부분 **정적 빔 형태**
  (양끝이 트랩에 잡혀 있으니 당연히 작아진다)이지 감쇠가 아니다. 실측이 그것을
  보여준다 — k_t×100 은 λ=2.6 ≪ 1/k=42.5 로, 진폭은 빨리 죽는데 위상은 거의
  안 돈다. 진짜 확산형 파동이면 λ ≈ 1/k 여야 한다.
  → **동역학 정보는 위상 기울기에만 있다.** λ 는 경계조건 진단용으로만 본다.

★★ 결맞음(coherence) 을 먼저 본다 — 시드 간 위상이 재현되지 않으면 v 는 무의미하다.
  DLVO 는 coh ≈ 0.46~0.59 로, 내부 비드 위상이 구동과 무관한 열잡음이다.
  적합하면 숫자는 나오지만 오차가 값보다 크고 λ 가 음수로 나온다 — 인용 금지.

⚠️ 양끝 비드(0·n−1)는 트랩된 힘센서라 경계조건이 지배한다 — 적합에서 뺀다.

    $PY scratch/propagation_speed.py
"""
import glob
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
import sys                                                    # noqa: E402
sys.path.insert(0, str(ROOT))
from bdbot import lockin as LI                                # noqa: E402

d_SI = 1.47e-6
ELL = 1.0076                       # 결합길이 [d]
TAU_B = 6.1510                     # s — ω* → ω[rad/s] 환산용
P = "runs/chain-bend-2d-dlvo__n9-w3000-a1470"

CONDS = [("trap, k_t 기본",  "",           "-jkr"),
         ("trap, k_t×100",   "-kt100",     "-jkr-kt100"),
         ("position 강제",    "-position",  "-jkr-position")]


def phasors(dd):
    """비드별 복소 위상자 ŷ_i (전 표본 락인)."""
    z = np.load(Path(dd) / "observables.npz", allow_pickle=True)
    if "shape_y" not in set(z.files):
        return None, None, None
    s = json.loads((Path(dd) / "spec.json").read_text())
    if s["numerics"]["n_prod"] < 600_000:              # 스모크 런 제외
        return None, None, None
    om = float(s["params"]["omega_star"])
    t = np.asarray(z["t"], float)
    ys = np.asarray(z["shape_y"], float)
    out = []
    for i in range(ys.shape[1]):
        b = LI.lockin_blocks(t, ys[:, i] - ys[:, i].mean(), om,
                             n_blocks=min(10, max(2, len(t) // 20)))
        h, _ = LI.agg(b)
        out.append(h)
    return np.array(out), om, int(s["params"]["n_beads"])


def fit(pat):
    Zs, om, n = [], None, None
    for dd in sorted(glob.glob(str(ROOT / pat))):
        if not (Path(dd) / "metrics.json").exists():
            continue
        z, o, nn = phasors(dd)
        if z is None:
            continue
        Zs.append(z); om, n = o, nn
    if not Zs:
        return None
    Z = np.array(Zs)
    mid = n // 2
    Zrel = Z / Z[:, [mid]]                       # 구동 비드 기준 (절대 위상 오염 상쇠)

    # 내부 비드만 (양끝 트랩 제외), 한쪽 팔씩 → |s| 로 접어서 평균
    # ★ 결맞음 검사 — 위상 프로파일이 시드 간 재현되는가.
    #   coh = |시드평균 위상자| / 평균|위상자|.  1 이면 완전 결맞음, 0 이면 무작위.
    #   DLVO 처럼 사슬이 힘을 전달 못 하면 내부 비드 위상이 열잡음이라 coh → 0 이고,
    #   그때 v 를 적합해봐야 **의미 없는 숫자**가 나온다 (실측: 오차 > 값, λ < 0).
    coh = np.abs(Zrel.mean(0)) / np.abs(Zrel).mean(0)

    res = []
    for zr in Zrel:
        s = np.arange(n) - mid
        keep = (np.abs(s) >= 1) & (np.abs(s) <= mid - 1)       # s=±1..±3 (n=9)
        ss = np.abs(s[keep]).astype(float)
        ph = np.unwrap(np.angle(zr[keep]))
        amp = np.abs(zr[keep])
        # φ(s) 기울기 [rad/결합]  ·  ln|A| 기울기 → 침투깊이
        kk = -np.polyfit(ss, ph, 1)[0]
        lam = -1.0 / np.polyfit(ss, np.log(amp), 1)[0]
        res.append((kk, lam))
    res = np.array(res)
    k_m, k_e = res[:, 0].mean(), res[:, 0].std(ddof=1) / np.sqrt(len(res))
    l_m, l_e = res[:, 1].mean(), res[:, 1].std(ddof=1) / np.sqrt(len(res))
    om_SI = om / TAU_B                                          # rad/s
    v = om_SI * ELL * d_SI / k_m if k_m else np.inf             # m/s
    v_e = abs(v) * (k_e / abs(k_m)) if k_m else np.nan
    ss_all = np.arange(n) - mid
    keep_all = (np.abs(ss_all) >= 1) & (np.abs(ss_all) <= mid - 1)
    coh_in = float(coh[keep_all].mean())
    return dict(k=k_m, k_e=k_e, lam=l_m, lam_e=l_e, v=v, v_e=v_e,
                om_SI=om_SI, n=len(res), coh=coh_in)


print("=" * 104)
print("응력 전파 속도 — 위상 기울기  v = ωℓ/|dφ/ds|   (ω = 3000 rad/s, n=9, a=1d)")
print("=" * 104)
COH_MIN = 0.5     # ★ 이 아래면 위상이 결맞지 않아 v 를 정의할 수 없다
print(f"{'조건':<16} {'분기':<6} {'결맞음':>7} {'dφ/ds [rad/결합]':>19} {'v [µm/s]':>20} "
      f"{'λ [결합]':>15} {'1/k':>7}")
rows = []
for lab, sd, sj in CONDS:
    # ★ DLVO 는 제외한다 (사용자 지정 2026-08-06). 이미 결맞음 0.46~0.59 로
    #   v 가 정의되지 않음을 확인했다 — 전파할 응력 자체가 없으니 속도도 없다.
    for br, suf in (("JKR", sj),):
        r = fit(f"{P}{suf}__*")
        if r is None:
            print(f"{lab:<16} {br:<6} {'—':>19}")
            continue
        if r["coh"] < COH_MIN:
            print(f"{lab:<16} {br:<6} {r['coh']:>7.3f}   ★ 위상이 결맞지 않음 — "
                  f"v 정의 불가 (내부 비드가 구동과 무관하게 움직인다)")
        else:
            print(f"{lab:<16} {br:<6} {r['coh']:>7.3f} {r['k']:>12.5f}±{r['k_e']:<6.4f} "
                  f"{r['v']*1e6:>13.0f}±{r['v_e']*1e6:<6.0f} "
                  f"{r['lam']:>9.3f}±{r['lam_e']:<5.3f} {1/r['k']:>7.2f}")
        rows.append(dict(lab=lab, br=br, **r))
    print()

import pickle                                                  # noqa: E402
pickle.dump(rows, open("/private/tmp/claude-501/-Users-kyuhwan-Desktop-simulation-auto/"
                       "7f5025d1-46c3-455c-a0ba-f595731412cc/scratchpad/prop_speed.pkl", "wb"))
print("=" * 104)
print("읽는 법 —  ① 결맞음 < 0.5 면 v 자체가 정의되지 않는다 (DLVO 가 그렇다).")
print("           ② λ 는 침투깊이가 **아니다** — 3점굽힘의 정적 빔 형태가 지배한다.")
print("              k_t×100 의 λ=2.6 ≪ 1/k=42.5 가 그 증거 (확산형이면 λ≈1/k).")
print("           ③ k_t 기본과 position 이 v≈30,000 µm/s 로 같다 — 둘 다 양끝 트랩이")
print("              기본 강성이라서다. k_t×100 은 양끝이 사실상 고정돼 6배 빨라진다.")
print("★ 이 값은 **ω=3000 rad/s 에서의 값**이다. 과감쇠 굽힘은 ω~k⁴ 이라 분산성이라서")
print("  v_phase 가 ω 에 의존한다 — 물질 상수가 아니다.")
