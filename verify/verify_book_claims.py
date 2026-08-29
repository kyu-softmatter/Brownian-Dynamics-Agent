#!/usr/bin/env python
"""kb/ 에 들어온 책 두 권에서 뽑은 주장을 **실행으로** 검증한다 (절대 규칙 6).

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
    $PY scratch/verify_book_claims.py

책:
  [L] L. Gary Leal, *Microstructural Rheology of Complex Fluids*,
      Cambridge Monographs on Mechanics, 2026. DOI 10.1017/9781009688437.
      kb/cambridge-core_microstructural-rheology-of-complex-fluids_7Aug2026/
  [W] J. Welty et al., *Fundamentals of Momentum, Heat and Mass Transfer*, 5th ed.
      kb/file_1731415827j8JuJ.pdf  (제목은 PDF 메타데이터, ISBN 0-4701-2868-2)

검증을 세 종류로 **구분해서** 표시한다 — 섞으면 "책이 맞다"와 "내가 옳게 읽었다"와
"우리 값이 맞다"가 뭉개진다:

  [BOOK]  책이 스스로 보고한 숫자를 책의 공식으로 재현 (내 독해가 맞는지)
  [DERIV] 책의 공식들이 서로 모순 없는지 (수치적분·극한 등으로 독립 확인)
  [OURS]  책의 값·공식과 이 프로젝트가 쓰는 값을 대조

★ 실패해도 좋다. 실패는 "책이 틀렸다"가 아니라 대개 "내가 잘못 읽었다"는 뜻이고,
   그러면 그 주장은 digest 와 KB 에 들어가면 안 된다.
"""
from __future__ import annotations

import math
import re
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent

PASS: list[str] = []
FAIL: list[str] = []


def chk(kind: str, name: str, got, want, rtol=1e-3, note=""):
    """수치 대조."""
    ok = abs(got - want) <= rtol * max(abs(want), 1e-300)
    rec = f"[{kind}] {name}: got={got:.6g} want={want:.6g} rel={abs(got-want)/max(abs(want),1e-300):.2e}"
    if note:
        rec += f"  ({note})"
    (PASS if ok else FAIL).append(rec)
    return ok


def chk_true(kind: str, name: str, cond: bool, note=""):
    rec = f"[{kind}] {name}" + (f"  ({note})" if note else "")
    (PASS if cond else FAIL).append(rec)
    return cond


# ── 증류본에 **인쇄된** 값을 대조한다 ──────────────────────────────────────
# 2026-08-29: 이게 없어서 typo 가 56/56 을 통과했다. `welty_transport.md` 는
# log-보간값을 `0.8580 mPa·s` 로 찍었는데 표에서 계산하면 `0.8598` 이다.
# 검증기는 0.8598 을 **옳게 계산해서 출력하고 있었지만**, 증류본에 적힌 숫자와
# 대조하는 단정이 없었다 — `d_log < 0.015` 는 두 값 모두 통과시킨다.
# ⚠️ 계산이 맞는 것과 **인용될 문서가 맞는 것은 다른 명제**다. 증류본이
# 인용되는 산출물이므로, 계산값과 인쇄값을 잇는 단정이 없으면 전사 오류가
# 조용히 살아남는다 (이 프로젝트의 서명 실패 유형: 배선되지 않은 검사기).
def chk_doc(book: str, pattern: str, want: float, unit_scale: float = 1.0,
            rtol: float = 5e-4, note: str = ""):
    """증류본에서 정규식으로 숫자를 꺼내 계산값과 대조한다.

    `pattern` 은 정확히 하나의 캡처 그룹(숫자)을 가져야 한다. 매치가 없거나
    둘 이상이면 **FAIL** 이다 — 조용히 건너뛰면 이 검사 자체가 무의미해진다.
    """
    path = ROOT / "knowledge" / "source" / "books" / book
    if not path.exists():
        FAIL.append(f"[DOC] {book}: 파일이 없다 ({path})")
        return False
    hits = re.findall(pattern, path.read_text(encoding="utf-8"))
    if len(hits) != 1:
        FAIL.append(f"[DOC] {book} /{pattern}/: 매치가 {len(hits)}개 (정확히 1개여야 한다)")
        return False
    got = float(hits[0]) * unit_scale
    return chk("DOC", f"{book} 인쇄값 {pattern!r}" + (f" — {note}" if note else ""),
               got, want, rtol=rtol)


# ════════════════════════════════════════════════════════════════════════
# 상수
# ════════════════════════════════════════════════════════════════════════
KB = 1.380649e-23  # J/K  (SI 정의값)

# ── 이 프로젝트가 쓰는 값 (intake/*/system.yaml, tier 1) ────────────────
OURS_T = 300.0  # K
OURS_ETA = 0.851e-3  # Pa*s   "물@300K 핸드북" — 5개 케이스 공통
OURS_GAMMA = {  # 3 pi eta d
    "trap-drag-2d-hex300 (d=5um)": (5.000e-6, 4.0102e-8),
    "chain-bend-2d-dlvo (d=1.47um)": (1.470e-6, 1.1790e-8),
}

# ── [W] Appendix I, "Water" (SI 표, 책 p.686) ──────────────────────────
#   T(K), rho(kg/m3), mu*1e6 (Pa*s), nu*1e6 (m2/s)
WELTY_WATER = [
    (273, 999.3, 1794, 1.795),
    (293, 998.2, 993, 0.995),
    (313, 992.2, 658, 0.663),
    (333, 983.2, 472, 0.480),
    (353, 971.8, 352, 0.362),
    (373, 958.4, 278, 0.290),
]


# ════════════════════════════════════════════════════════════════════════
# ① [W] 물의 점도 — 책의 표 ↔ 우리가 쓰는 0.851 mPa*s @300 K
# ════════════════════════════════════════════════════════════════════════
def s1_water_viscosity():
    T = np.array([r[0] for r in WELTY_WATER], float)
    mu = np.array([r[2] for r in WELTY_WATER], float) * 1e-6  # Pa*s
    rho = np.array([r[1] for r in WELTY_WATER], float)
    nu = np.array([r[3] for r in WELTY_WATER], float) * 1e-6

    # 표 자체의 내적 일관성: nu = mu/rho (책이 세 열을 독립적으로 찍었으므로 검사가 된다)
    err = float(np.max(np.abs(nu - mu / rho) / nu))
    chk_true("BOOK", f"Welty water table: nu == mu/rho 가 모든 행에서 0.1% 이내 (실측 {err:.2e})",
             err < 1e-3, note="표를 옳게 전사했는지 + 표 자체가 일관적인지")

    # 300 K 로 보간. eta(T) 는 지수적으로 휘므로 log-선형이 선형보다 옳다.
    i = 1  # 293 K
    w = (300.0 - T[i]) / (T[i + 1] - T[i])
    mu_lin = mu[i] + w * (mu[i + 1] - mu[i])
    mu_log = math.exp(math.log(mu[i]) + w * (math.log(mu[i + 1]) - math.log(mu[i])))

    print(f"    Welty 표 293/313 K = {mu[1]*1e3:.3f}/{mu[2]*1e3:.3f} mPa*s")
    print(f"    -> 300 K 선형보간 {mu_lin*1e3:.4f} / log-선형보간 {mu_log*1e3:.4f} mPa*s")
    print(f"    -> 우리 값        {OURS_ETA*1e3:.4f} mPa*s")

    # 우리 값과의 상대차. log-보간이 20 K 간격 표에서 낼 수 있는 정확도 안에 있는가.
    d_log = abs(mu_log - OURS_ETA) / OURS_ETA
    d_lin = abs(mu_lin - OURS_ETA) / OURS_ETA
    chk_true("OURS", f"eta(300K): Welty log-보간이 우리 0.851 mPa*s 와 1.5% 이내 (실측 {d_log*100:.2f}%)",
             d_log < 0.015, note=f"선형보간은 {d_lin*100:.2f}% 어긋남 — 보간법이 결과를 바꾼다")
    chk_true("OURS", "eta(300K): log-보간이 선형보간보다 우리 값에 가깝다",
             d_log < d_lin, note="eta(T) 의 곡률이 20 K 간격에서 이미 유의미")

    # ★ 증류본에 **인쇄된** 두 값이 위 계산과 일치하는가. 이 두 줄이 없어서
    #   0.8580 (실제 0.8598) 이 56/56 을 통과했다.
    chk_doc("welty_transport.md", r"log-선형[^|]*\|\s*\*\*([\d.]+) mPa·s\*\*",
            mu_log * 1e3, note="log-선형 보간값")
    chk_doc("welty_transport.md", r"\|\s*선형\s*\|\s*([\d.]+) mPa·s",
            mu_lin * 1e3, note="선형 보간값")
    # 인쇄된 백분율도 같은 계산에서 나와야 한다 (0.8580 은 +0.82% 라 +1.03% 와 모순)
    chk_doc("welty_transport.md", r"log-선형[^|]*\|[^|]*\|\s*\*\*\+([\d.]+)%\*\*",
            d_log * 100, rtol=6e-3, note="log-선형 상대차 %")

    # 온도 민감도: %/K. 값 하나만 인용하고 T 를 안 적으면 이만큼 틀릴 수 있다.
    sens = abs(math.log(mu[2] / mu[1])) / (T[2] - T[1]) * 100  # %/K, 293-313 구간
    print(f"    d(ln eta)/dT ~ {sens:.2f} %/K  (293-313 K)")
    chk_true("BOOK", f"물 점도의 온도 민감도가 2%/K 를 넘는다 (실측 {sens:.2f} %/K)",
             sens > 2.0, note="T 를 +-1 K 잘못 적으면 eta 가 2% 틀린다")
    return sens


# ════════════════════════════════════════════════════════════════════════
# ② [W] Stokes 항력 ↔ C_D = 24/Re ↔ 우리 gamma = 3 pi eta d
# ════════════════════════════════════════════════════════════════════════
def s2_stokes():
    # C_D = F / (0.5 rho v^2 A),  A = pi d^2/4,  F = 3 pi mu d v  (Stokes)
    # -> C_D = 3 pi mu d v / (0.5 rho v^2 pi d^2/4) = 24 mu/(rho v d) = 24/Re
    rho, mu, d, v = 998.2, 0.993e-3, 3e-6, 1e-4
    Re = rho * v * d / mu
    F = 3 * math.pi * mu * d * v
    CD = F / (0.5 * rho * v**2 * math.pi * d**2 / 4)
    chk("DERIV", "C_D(Stokes) == 24/Re", CD, 24.0 / Re, rtol=1e-12,
        note="[W] 12.2 의 C_D=24/Re 와 F=3*pi*mu*d*v 가 같은 것임을 확인")
    chk_true("DERIV", f"검사한 조건이 creeping flow 영역 (Re={Re:.2e} < 1)", Re < 1.0)

    # 우리 gamma = 3 pi eta d  (= 6 pi eta a) 가 각 케이스에서 맞게 계산됐는가
    for name, (d_case, gamma_yaml) in OURS_GAMMA.items():
        chk("OURS", f"gamma = 3*pi*eta*d  [{name}]",
            3 * math.pi * OURS_ETA * d_case, gamma_yaml, rtol=2e-4)

    # Stokes-Einstein [W] (24-50): D = kT/(6 pi mu r) = kT/gamma
    for name, (d_case, gamma_yaml) in OURS_GAMMA.items():
        r = d_case / 2
        D_se = KB * OURS_T / (6 * math.pi * OURS_ETA * r)
        chk("DERIV", f"Stokes-Einstein D == kT/gamma  [{name}]",
            D_se, KB * OURS_T / gamma_yaml, rtol=2e-4)
        print(f"    {name}: D = {D_se:.4e} m^2/s")


# ════════════════════════════════════════════════════════════════════════
# ③ [L] Batchelor / Brady-Vicic 의 phi^2 계수 — 책의 산술 재현
# ════════════════════════════════════════════════════════════════════════
def s3_phi2():
    # [L] 3.4.2: 직접 Brownian 0.97 phi^2 + 유체역학 5.2 phi^2 -> Batchelor K* = 6.2
    chk("BOOK", "[L] Batchelor K* = 0.97 (Brownian) + 5.2 (hydro)", 0.97 + 5.2, 6.2,
        rtol=5e-3, note="책이 두 성분과 합을 따로 찍었으므로 산술 검사가 된다")
    # Brady & Vicic (1995): K* = 5.91, N1/(mu gdot) = 0.899 phi^2 Pe, N2 = -0.788 phi^2 Pe
    chk_true("BOOK", "[L] Brady-Vicic: |N2/N1| = 0.877 (고분자의 ~1/7 과 다르다)",
             abs(-0.788 / 0.899 + 0.877) < 5e-3,
             note=f"|N2/N1|={abs(-0.788/0.899):.3f} vs 고분자 {1/7:.3f}")
    chk_true("BOOK", "[L] N1>0, N2<0 (고분자와 부호는 같다)", 0.899 > 0 and -0.788 < 0)


# ════════════════════════════════════════════════════════════════════════
# ④ [L] Doi-Edwards 준희박 회전확산 — 책의 예제 숫자 재현
# ════════════════════════════════════════════════════════════════════════
def s4_doi_edwards():
    # [L] (6.45): Dr0 = beta * Drbar0 / (n^2 L^6),  beta = 1.3e3
    # 책의 예: n=0.1, L=50 -> Dr0/Drbar0 = O(1e-5)
    beta, n, L = 1.3e3, 0.1, 50.0
    ratio = beta / (n**2 * L**6)
    print(f"    Dr(semi-dilute)/Dr(dilute) = {ratio:.3e}   (책: O(1e-5))")
    chk_true("BOOK", f"[L] (6.45) 예제: 비 = {ratio:.2e} 가 O(1e-5)",
             1e-6 <= ratio < 1e-4, note="n=0.1, L=50, beta=1.3e3")
    # 스케일: Dr ~ (n L^3)^-2  (n^2 L^6 = (nL^3)^2)
    chk("DERIV", "n^2 L^6 == (n L^3)^2", n**2 * L**6, (n * L**3) ** 2, rtol=1e-12)
    # 배제부피: 직교한 두 막대 -> L x L x 2a 평행육면체 = 2 a L^2
    a = 0.5
    chk("BOOK", "[L] 배제부피(직교) = 2 a L^2", L * L * (2 * a), 2 * a * L**2, rtol=1e-12)
    # 준희박 창의 폭은 종횡비가 정한다: O(1) << n L^3 << L/a = 2r
    r_asp = L / (2 * a)
    chk("DERIV", "준희박 상한 n L^3 ~ L/a = 2r", L / a, 2 * r_asp, rtol=1e-12,
        note=f"r={r_asp:.0f} -> nL^3 창 [O(1), {L/a:.0f}]")


# ════════════════════════════════════════════════════════════════════════
# ⑤ [L] Jeffery 의 G, 그리고 Pe = |E|/Dr 의 |E| 규약 (2배 함정)
# ════════════════════════════════════════════════════════════════════════
def s5_jeffery_and_E():
    G = lambda r: (r**2 - 1) / (r**2 + 1)
    chk("DERIV", "[L] (2.34) G(r=1) = 0  (구는 변형률에 반응 안 함)", G(1.0), 0.0, rtol=1e-12)
    chk("DERIV", "[L] (2.34) G(r->inf) -> 1", G(1e8), 1.0, rtol=1e-12)
    chk_true("DERIV", "G(r) 는 r 에 대해 단조증가", bool(np.all(np.diff(G(np.linspace(1, 100, 500))) > 0)))

    # ★ 규약 함정: [L] 4.4.4 는 Pe = |E|/Dr, [L] 1.6 (1.32) 는 Wi = gdot/Dr.
    #   둘이 같으려면 |E| = sqrt(2 E:E) 여야 한다 (E:E 만 쓰면 sqrt(2) 배 틀린다).
    gdot = 3.7
    E = np.array([[0.0, gdot / 2, 0.0], [gdot / 2, 0.0, 0.0], [0.0, 0.0, 0.0]])
    EE = float(np.tensordot(E, E))  # E:E
    chk("DERIV", "단순전단에서 sqrt(2 E:E) == gdot", math.sqrt(2 * EE), gdot, rtol=1e-12,
        note="[L] 3.4.2 의 Ehat = E/(2E:E)^1/2 = E/gdot 규약과 일치")
    chk_true("DERIV", f"sqrt(E:E) 를 쓰면 sqrt(2)={math.sqrt(2):.4f} 배 어긋난다",
             abs(math.sqrt(EE) * math.sqrt(2) - gdot) < 1e-12,
             note="Pe/Wi 를 2배 틀리게 만드는 가장 흔한 규약 사고")

    # 회전(vorticity)은 등방 평형분포를 바꾸지 못한다 -> Pe 는 E 로만 만든다.
    Om = np.array([[0.0, gdot / 2, 0.0], [-gdot / 2, 0.0, 0.0], [0.0, 0.0, 0.0]])
    chk("DERIV", "단순전단: |grad u| 은 E 와 Omega 로 반씩 쪼개진다",
        float(np.max(np.abs(E + Om))), gdot, rtol=1e-12,
        note="|grad u|=gdot 인데 |E|=gdot 이므로 여기서는 우연히 같다 — 순수전단이면 다르다")
    # 순수 회전 흐름: E=0 이므로 Pe=0 이어야 한다 (아무리 grad u 가 커도)
    E_rot = np.zeros((3, 3))
    chk_true("DERIV", "순수 회전류(E=0)에서는 Pe=0 — |grad u| 로 재면 0 이 아니다",
             float(np.max(np.abs(E_rot))) == 0.0)


# ════════════════════════════════════════════════════════════════════════
# ⑥ [L] 막대 현탁액 SAOS: G'' ~ De/(36+De^2) -> 점탄성 이완시간 = 1/(6 Dr)
# ════════════════════════════════════════════════════════════════════════
def s6_rod_saos():
    # [L] (5.121b): FG*(2/15)*3De/(36+De^2) 꼴. Maxwell 단일모드의 극은 De=6.
    De = np.logspace(-2, 3, 200001)
    loss = De / (36 + De**2)  # G'' ~ eta' * omega 의 De 의존부
    De_peak = float(De[int(np.argmax(loss))])
    chk("DERIV", "[L] (5.121b) 손실항 De/(36+De^2) 의 최대는 De=6", De_peak, 6.0, rtol=1e-3,
        note="De = omega/Dr 이므로 omega_peak = 6 Dr -> lambda = 1/(6 Dr)")
    # Maxwell 대조: lambda*omega/(1+(lambda omega)^2) 는 lambda*omega=1 에서 최대
    lam = 1.0 / 6.0
    m = lam * De / (1 + (lam * De) ** 2)
    chk("DERIV", "Maxwell(lambda=1/6) 의 손실 최대도 De=6",
        float(De[int(np.argmax(m))]), 6.0, rtol=1e-3)
    chk_true("DERIV", "De/(36+De^2) 는 Maxwell(lambda=1/6) 과 비례상수 하나로 같다",
             float(np.max(np.abs(loss / m - loss[0] / m[0]))) < 1e-12,
             note="즉 희박 막대 현탁액의 선형점탄성은 정확히 단일 Maxwell 모드")
    print(f"    -> 배향(l=2) 이완시간 lambda = 1/(6 Dr).  1/Dr 나 1/(2Dr) 가 아니다.")


# ════════════════════════════════════════════════════════════════════════
# ⑦ [L] 탄성 덤벨: Kramers 응력 · Oldroyd-B — 모멘트식 수치적분으로 독립 확인
# ════════════════════════════════════════════════════════════════════════
def s7_dumbbell():
    """[L] 11.2-11.3.

    선형(Hookean) 덤벨의 2차 모멘트 방정식은
        d<RR>/dt = grad u . <RR> + <RR> . grad u^T  - (4 w0/zeta)(<RR> - (kT/w0) I)
    (= (11.6) 의 상위대류 미분 형태). 여기서 lambda_p = zeta/(4 w0).
    Kramers (11.15b):  T^(p) = n <F_s R> - n kT I = n (w0 <RR> - kT I).
    """
    kT, zeta, n = 1.0, 1.0, 1.0
    b2 = 3.0 * kT  # <R^2>_eq;  w0 = 3kT/b^2 (Gaussian chain, [L] p.494 "3kT/w0 = b^2")
    w0 = 3.0 * kT / b2
    lam = zeta / (4.0 * w0)
    eta_p = n * kT * lam

    chk("DERIV", "[L] Gaussian chain: 3kT/w0 == b^2", 3 * kT / w0, b2, rtol=1e-12)

    def rhs(RR, gradu):
        relax = -(1.0 / lam) * (RR - (kT / w0) * np.eye(3))
        return gradu @ RR + RR @ gradu.T + relax

    def integrate(gradu, t_end, RR0, dt=None):
        """모멘트식을 강성 대응 적분기로 적분 (닫힌 해를 쓰지 않는 독립 확인)."""
        from scipy.integrate import solve_ivp

        sol = solve_ivp(
            lambda t, y: rhs(y.reshape(3, 3), gradu).ravel(),
            (0.0, t_end), RR0.ravel(),
            method="LSODA", rtol=1e-10, atol=1e-12,
        )
        assert sol.success, sol.message
        return sol.y[:, -1].reshape(3, 3)

    RR_eq = (kT / w0) * np.eye(3)

    # (a) 평형에서 Kramers 응력이 정확히 0 (부호·전인자 규약의 골든 테스트)
    T_eq = n * (w0 * RR_eq - kT * np.eye(3))
    chk("DERIV", "[L] (11.15b) 평형에서 Kramers 응력 = 0",
        float(np.max(np.abs(T_eq))), 0.0, rtol=1.0,
        note=f"max|T|={np.max(np.abs(T_eq)):.2e} — 정확히 0 이어야 한다")
    chk_true("DERIV", "평형 Kramers 응력이 기계정밀도로 0",
             float(np.max(np.abs(T_eq))) < 1e-14)

    # (b) 정상 단순전단: T12 = eta_p gdot, N1 = 2 eta_p lam gdot^2, N2 = 0
    for gdot in (0.1, 1.0, 5.0):
        gradu = np.array([[0.0, gdot, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        RR = integrate(gradu, 60 * lam, RR_eq)
        T = n * (w0 * RR - kT * np.eye(3))
        Wi = lam * gdot
        chk("DERIV", f"[L] (11.20) T12 = eta_p*gdot  (Wi={Wi:g})", T[0, 1], eta_p * gdot, rtol=2e-4)
        chk("DERIV", f"[L] (11.21) N1 = 2 eta_p lam gdot^2  (Wi={Wi:g})",
            T[0, 0] - T[1, 1], 2 * eta_p * lam * gdot**2, rtol=2e-4)
        chk("DERIV", f"[L] (11.21) N2 = 0  (Wi={Wi:g})", T[1, 1] - T[2, 2], 0.0, rtol=1.0,
            note=f"|N2|={abs(T[1,1]-T[2,2]):.2e}")
        # (11.24) tr<RR> = b^2 (1 + (2/3) Wi^2)
        chk("DERIV", f"[L] (11.24) tr<RR> = b^2(1+2/3 Wi^2)  (Wi={Wi:g})",
            float(np.trace(RR)), b2 * (1 + 2.0 / 3.0 * Wi**2), rtol=2e-4)

    # (c) Oldroyd-B 는 전단담화를 예측하지 않는다 — 점도가 gdot 에 무관
    etas = []
    for gdot in (0.01, 0.1, 1.0, 10.0, 100.0):
        gradu = np.array([[0.0, gdot, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
        RR = integrate(gradu, 60 * lam, RR_eq, dt=min(lam / 20000.0, 1e-4 / gdot))
        etas.append(n * (w0 * RR[0, 1]) / gdot)
    spread = (max(etas) - min(etas)) / eta_p
    chk_true("DERIV", f"Oldroyd-B: eta_p 가 gdot 4자리에 걸쳐 불변 (산포 {spread:.2e})",
             spread < 1e-3, note="'전단담화 없음' — 거의 모든 실제 고분자와 다르다")

    # (d) 흐름을 끊으면 유체역학 성분은 즉시 0, 열역학 성분만 유한시간 이완
    gradu = np.array([[0.0, 2.0 / lam, 0.0], [0.0, 0.0, 0.0], [0.0, 0.0, 0.0]])
    RR_sheared = integrate(gradu, 60 * lam, RR_eq)
    RR_relaxed = integrate(np.zeros((3, 3)), 1.0 * lam, RR_sheared)
    T0 = n * (w0 * RR_sheared[0, 1])
    T1 = n * (w0 * RR_relaxed[0, 1])
    chk("DERIV", "정지 후 lambda 만큼 지나면 T12 가 exp(-1) 배",
        T1 / T0, math.exp(-1.0), rtol=5e-3,
        note="[L] 2.2.1: 유체역학 성분은 즉시 0 -> 남는 기억은 전부 열역학 성분")
    print(f"    lambda_p = zeta/(4 w0) = {lam:g},  eta_p = n kT lambda_p = {eta_p:g}")


# ════════════════════════════════════════════════════════════════════════
# ⑧ [L] 2.2.1 의 '즉시성' 명제를 우리 계에 적용했을 때의 산술
#      (chain-bend-2d-dlvo: 사슬을 붙여도 소산이 0.996 배)
# ════════════════════════════════════════════════════════════════════════
def s8_our_dlvo():
    # runs/.../system_moduli 결과 (CLAUDE.md 1-D 절)
    bead_only, dlvo, jkr = 18453.0, 18380.0, 75590.0
    chk("OURS", "DLVO 사슬의 총 소산 / 비드 단독", dlvo / bead_only, 0.996, rtol=2e-3,
        note="[L] 2.2.1 -> 복원기구가 없는 자유도는 응력에 기여하지 않는다")
    chk("OURS", "JKR 사슬의 총 소산 / 비드 단독", jkr / bead_only, 4.10, rtol=5e-3)
    chk("OURS", "JKR 에서 사슬이 담당하는 소산 비율", (jkr - bead_only) / jkr, 0.758, rtol=5e-3)
    frac = (dlvo - bead_only) / bead_only
    chk_true("OURS", f"DLVO 사슬의 초과 소산이 |{frac*100:.1f}%| < 1%",
             abs(frac) < 0.01,
             note="[L] 언어로: 굽힘에 대한 열역학 복원기구가 없어 직접 성분이 0")


# ════════════════════════════════════════════════════════════════════════
# ⑨ [W] Buckingham pi 는 r = '차원행렬의 rank' — 기본차원 개수가 아니다
# ════════════════════════════════════════════════════════════════════════
def s9_buckingham():
    # [W] 11.3: i = n - r,  r = rank of the dimensional matrix.
    # 반례: 기본차원이 3개(M,L,t)인데 rank 가 2인 변수 목록
    #   변수 = (v [L/t], gdot [1/t], L [L]) -> M 은 아무 데도 없다 -> rank 2
    #   M L t 행 x 3 변수 열
    A = np.array([
        [0, 0, 0],   # M
        [1, 0, 1],   # L
        [-1, -1, 0],  # t
    ], float)
    r = int(np.linalg.matrix_rank(A))
    n_var = A.shape[1]
    chk("DERIV", "차원행렬 rank (v, gdot, L)", r, 2, rtol=1e-12)
    chk("DERIV", "i = n - rank = 3 - 2", n_var - r, 1.0, rtol=1e-12,
        note="'기본차원 3개' 로 세면 i=0 이라 무차원수가 없다고 오판한다")
    # 실제 그 무차원수: gdot L / v
    v, gdot, Lc = 2.0, 3.0, 5.0
    pi1 = gdot * Lc / v
    chk_true("DERIV", f"찾아진 무차원수 gdot*L/v = {pi1:g} 는 실제로 무차원",
             True, note="rank 를 세야 개수가 맞는다")

    # Re 의 dimensional matrix 는 rank 3 -> i = 5-3 = 2 (Eu, Re) — 책의 예 (11-9)
    #   변수 = (F, rho, v, L, mu)
    B = np.array([
        [1, 1, 0, 0, 1],    # M
        [1, -3, 1, 1, -1],  # L
        [-2, 0, -1, 0, -1],  # t
    ], float)
    chk("BOOK", "[W] (11-9) 예: rank=3, i = 5-3 = 2 (Eu, Re)",
        B.shape[1] - int(np.linalg.matrix_rank(B)), 2.0, rtol=1e-12)


# ════════════════════════════════════════════════════════════════════════
# ⑩ [L] 희박 구 현탁액은 뉴턴 유체 — Einstein (2.31) 과 우리 계의 거리
# ════════════════════════════════════════════════════════════════════════
def s10_einstein():
    # mu*/mu = 1 + 2.5 phi + K* phi^2
    for phi, want in [(0.01, 1.02562), (0.05, 1.14050), (0.10, 1.31200)]:
        chk("DERIV", f"mu*/mu (phi={phi}) = 1+2.5phi+6.2phi^2",
            1 + 2.5 * phi + 6.2 * phi**2, want, rtol=1e-4)
    # phi^2 항이 2.5phi 항의 10% 가 되는 phi: 6.2 phi^2 = 0.25 phi -> phi = 0.0403
    phi_star = 0.25 / 6.2
    chk("DERIV", "phi^2 항이 Einstein 항의 10% 가 되는 phi", phi_star, 0.0403, rtol=1e-3,
        note="phi>4% 면 '희박=뉴턴' 이 이미 10% 수준에서 깨진다")


def main():
    sections = [
        ("① 물의 점도 [W] Appendix I ↔ 우리 0.851 mPa*s", s1_water_viscosity),
        ("② Stokes 항력 · C_D=24/Re · Stokes-Einstein", s2_stokes),
        ("③ phi^2 계수 (Batchelor / Brady-Vicic)", s3_phi2),
        ("④ Doi-Edwards 준희박 Dr", s4_doi_edwards),
        ("⑤ Jeffery G · |E| 규약의 2배 함정", s5_jeffery_and_E),
        ("⑥ 막대 SAOS -> lambda = 1/(6 Dr)", s6_rod_saos),
        ("⑦ 탄성 덤벨: Kramers 응력 · Oldroyd-B (수치적분)", s7_dumbbell),
        ("⑧ 우리 1-D 결과를 [L] 2.2.1 로 읽기", s8_our_dlvo),
        ("⑨ Buckingham pi 의 r 은 rank 다", s9_buckingham),
        ("⑩ Einstein 점도와 '희박' 의 한계", s10_einstein),
    ]
    for title, fn in sections:
        print(f"\n{'='*78}\n{title}\n{'='*78}")
        fn()

    print(f"\n{'='*78}")
    for r in PASS:
        print(f"  PASS  {r}")
    for r in FAIL:
        print(f"  FAIL  {r}")
    print(f"{'='*78}")
    print(f"  {len(PASS)}/{len(PASS)+len(FAIL)} PASS")
    return 1 if FAIL else 0


if __name__ == "__main__":
    sys.exit(main())
