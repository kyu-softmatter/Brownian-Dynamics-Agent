"""`trap-drag` 앙상블 평균 — 시드만 다른 런 여러 개를 모아 실현 간 산포를 본다.

**왜**: 단일 런에서 얻은 "ΔU = 4.7·n_def", "결함이 24.7 τ_int 에 사라진다", "이완이
평탄역 + 급락" 이 전부 **한 실현**의 결과였습니다. 전위 생성은 확률적이라 급락 시각이
시드마다 흔들릴 것이고, 그러면 **앙상블 평균에서는 급락이 뭉개져** 매끄러운 곡선이
나올 수 있습니다. 그 둘을 가르는 것이 이 스크립트의 목적입니다.

같이 보는 것: 단일 런의 블록 SEM 이 **실현 간 표준편차**와 맞는가.
블록 SEM 은 한 궤적 안의 상관만 보정하므로, 실현마다 다른 전위 배치 같은
'느린 자유도'를 놓치면 오차를 과소평가합니다. 앙상블이 있어야 확인됩니다.

    $PY scratch/trap_drag_ensemble.py [--glob 'trap-drag-2d-hex300__tr0.117647-s*']
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "cases"))

from bdbot import nondim as ND, stats as ST  # noqa: E402

PH_EQ, PH_DRAG, PH_RELAX = "평형", "끌기", "이완"


def load(d: Path):
    spec = ND.load(d / "spec.json")
    z = np.load(d / "observables.npz")
    r = json.loads((d / "metrics.json").read_text())["result"]
    ti = spec.reduced("times", "tau_int")
    t = z["_t_step"] * float(spec.numerics["dt_star"]) / ti
    return spec, z, r, t, z["phase"]


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--glob", default="trap-drag-2d-hex300__tr0.117647-s*")
    ap.add_argument("--out", default="runs/trap-drag-2d-hex300__ENSEMBLE")
    ap.add_argument("--include-legacy", action="store_true",
                    help="스텝 해상 미측정(2026-08-05 이전) 런도 섞는다 — 기본은 제외")
    args = ap.parse_args()

    dirs = sorted(p for p in (ROOT / "runs").glob(args.glob)
                  if (p / "observables.npz").exists() and (p / "metrics.json").exists())

    # ★ 태그가 같고 해시만 다른 **레거시 런이 같은 글롭에 걸립니다.**
    #   `tr0.117647-s1__<옛해시>` 와 `tr0.117647-s1__<새해시>` 가 둘 다 매칭되어,
    #   그냥 두면 옛 코드 런과 새 런을 **같은 평균에 섞고 시드를 두 번 셉니다.**
    #   판별자는 `step_drift_max_sigma` 유무입니다 — 그게 곧 "현재 코드로 돈 런"입니다.
    #   ⚠️ 조용히 버리지 않고 무엇을 뺐는지 반드시 찍습니다.
    def is_current(p: Path) -> bool:
        try:
            m = json.loads((p / "metrics.json").read_text())
            return "step_drift_max_sigma" in m.get("numerics", {})
        except Exception:
            return False

    if not args.include_legacy:
        legacy = [p for p in dirs if not is_current(p)]
        dirs = [p for p in dirs if is_current(p)]
        if legacy:
            print(f"⚠ 레거시 {len(legacy)}런 제외 (스텝 해상 미측정 = 옛 코드). "
                  f"섞으려면 --include-legacy:")
            for p in legacy:
                print(f"    - {p.name}")
    if len(dirs) < 2:
        print(f"런이 부족합니다 ({len(dirs)}개)", file=sys.stderr)
        return 1
    print(f"앙상블 대상 {len(dirs)}런: {', '.join(p.name.split('__')[1] for p in dirs)}")

    rows, rel_u, rel_n, drag_f = [], [], [], []
    t_rel = None
    for d in dirs:
        spec, z, r, t, ph = load(d)
        N = spec.params["N"]
        rel = ph == PH_RELAX
        tr = t[rel] - t[rel][0]
        if t_rel is None or len(tr) < len(t_rel):
            t_rel = tr
        rel_u.append((z["u_pair"][rel] - r["U_pair_eq"]) * N)     # ΔU_total(t)
        rel_n.append(z["n_def"][rel].astype(float))
        nd = z["n_def"][rel]
        t0 = tr[np.argmax(nd == 0)] if (nd == 0).any() else np.nan
        dmask = ph == PH_DRAG
        idx = np.flatnonzero(dmask)
        ss = idx[len(idx) // 2:]                                   # 정상상태 = 후반 절반
        fx = -spec.params["k_star"] * z["dx_probe"][ss]
        drag_f.append(fx)
        # ΔU vs n_def 기울기 (끌기+이완)
        m = dmask | rel
        x = z["n_def"][m].astype(float)
        y = (z["u_pair"][m] - r["U_pair_eq"]) * N
        slope = np.linalg.lstsq(np.vstack([x, np.ones_like(x)]).T, y, rcond=None)[0][0]
        rows.append(dict(seed=int(spec.numerics["seed"]), dU=r["dU_drive"],
                         dU_sem=r["dU_sem"], F=float(fx.mean()),
                         F_sem=float(ST.block_sem(fx)),
                         nd_drag=float(z["n_def"][ss].mean()),
                         nd_max=float(z["n_def"][dmask].max()),
                         t_zero=float(t0), slope=float(slope),
                         psi6=r["psi6"]["driven"]))

    n = min(len(a) for a in rel_u)
    U = np.array([a[:n] for a in rel_u])         # (M, n)
    D = np.array([a[:n] for a in rel_n])
    t_rel = t_rel[:n]
    M = len(dirs)

    print("=" * 96)
    print(f"앙상블 — 시드 {M}개 · 2격자 (3.22 d) · N={rows[0].get('N', 306)}")
    print("=" * 96)
    print(f"{'시드':>5}{'ΔU/입자':>10}{'(런 SEM)':>10}{'F_drag':>9}{'(런 SEM)':>10}"
          f"{'결함 구동':>10}{'최대':>6}{'결함0 시각':>11}{'기울기':>8}{'ψ6':>7}")
    for r_ in rows:
        tz = f"{r_['t_zero']:.1f}" if np.isfinite(r_["t_zero"]) else "안 감"
        print(f"{r_['seed']:>5}{r_['dU']:>10.4f}{r_['dU_sem']:>10.4f}{r_['F']:>9.1f}"
              f"{r_['F_sem']:>10.1f}{r_['nd_drag']:>10.2f}{r_['nd_max']:>6.0f}"
              f"{tz:>11}{r_['slope']:>8.2f}{r_['psi6']:>7.3f}")

    def summ(key):
        v = np.array([r_[key] for r_ in rows], dtype=float)
        v = v[np.isfinite(v)]
        return v.mean(), v.std(ddof=1), v.std(ddof=1) / np.sqrt(len(v)), len(v)

    print("-" * 96)
    print(f"{'':>5}{'앙상블 평균':>14}{'실현간 std':>12}{'앙상블 SEM':>12}"
          f"{'런 SEM 평균':>13}{'비(std/런SEM)':>14}")
    for key, lab, semkey in (("dU", "ΔU/입자", "dU_sem"), ("F", "F_drag", "F_sem")):
        m_, s_, e_, k_ = summ(key)
        run_sem = np.mean([r_[semkey] for r_ in rows])
        print(f"{lab:>5}{m_:>14.4f}{s_:>12.4f}{e_:>12.4f}{run_sem:>13.4f}"
              f"{s_/run_sem:>14.2f}")
    for key, lab in (("nd_drag", "결함(구동)"), ("t_zero", "결함0 시각"), ("slope", "기울기")):
        m_, s_, e_, k_ = summ(key)
        print(f"{lab:>5}{m_:>14.3f}{s_:>12.3f}{e_:>12.3f}"
              + (f"      (유효 {k_}/{M}개)" if k_ < M else ""))

    # ── 그림 ──────────────────────────────────────────────────────────
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.ticker import FuncFormatter, NullFormatter
    matplotlib.rcParams["font.family"] = ["Apple SD Gothic Neo", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False
    plain = FuncFormatter(lambda v, _: f"{v:g}")

    out = ROOT / args.out
    out.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))

    # ① ΔU(t) 개별 + 앙상블
    a = ax[0, 0]
    for u in U:
        a.plot(t_rel, u, "-", lw=.6, color="0.75")
    mu, se = U.mean(0), U.std(0, ddof=1) / np.sqrt(M)
    a.plot(t_rel, mu, "-", lw=2.2, color="tab:red", label=f"앙상블 평균 (M={M})")
    a.fill_between(t_rel, mu - se, mu + se, color="tab:red", alpha=.25, label="±SEM")
    a.axhline(0, color="green", ls="--", lw=1.2, label="평형")
    a.set(xlabel=r"이완 시작 후 $t/\tau_{int}$", ylabel=r"$\Delta U_{total}$ [kT]",
          title="① 이완 — 개별 실현(회색) vs 앙상블 평균")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ② log–log — 급락이 평균에서도 남는가
    a = ax[0, 1]
    pos = t_rel > 0
    edges = np.geomspace(t_rel[pos][0], t_rel[-1], 26)
    ctr, mm, ss = [], [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        s = (t_rel >= lo) & (t_rel < hi)
        if s.sum() >= 1:
            ctr.append(np.sqrt(lo * hi)); mm.append(U[:, s].mean()); ss.append(U[:, s].mean(1).std(ddof=1) / np.sqrt(M))
    ctr, mm, ss = map(np.array, (ctr, mm, ss))
    ok = mm > 0
    a.errorbar(ctr[ok], mm[ok], yerr=ss[ok], fmt="o-", ms=4, lw=1.5, capsize=2,
               color="tab:red", label="앙상블 ΔU")
    for u in U:
        ub = np.array([u[(t_rel >= lo) & (t_rel < hi)].mean()
                       for lo, hi in zip(edges[:-1], edges[1:])
                       if ((t_rel >= lo) & (t_rel < hi)).sum() >= 1])
        m2 = ub > 0
        a.loglog(ctr[m2], ub[m2], "-", lw=.6, color="0.8", zorder=0)
    a.set_xscale("log"); a.set_yscale("log")
    for p_, c_ in ((-0.5, "tab:orange"), (-1.0, "tab:blue")):
        a.loglog(ctr[ok], mm[ok][0] * (ctr[ok] / ctr[ok][0]) ** p_, "--", lw=1,
                 color=c_, alpha=.7, label=fr"$t^{{{p_:g}}}$")
    for axis in (a.xaxis, a.yaxis):
        axis.set_major_formatter(plain); axis.set_minor_formatter(NullFormatter())
    a.set(xlabel=r"$t/\tau_{int}$", ylabel=r"$\Delta U_{total}$ [kT]",
          title="② log–log — 개별의 급락이 평균에도 남는가")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    # ③ 결함
    a = ax[1, 0]
    for dd in D:
        a.plot(t_rel, dd, "-", lw=.6, color="0.8")
    mud, sed = D.mean(0), D.std(0, ddof=1) / np.sqrt(M)
    a.plot(t_rel, mud, "-", lw=2.2, color="tab:purple", label=f"앙상블 평균 (M={M})")
    a.fill_between(t_rel, mud - sed, mud + sed, color="tab:purple", alpha=.25)
    tz = [r_["t_zero"] for r_ in rows if np.isfinite(r_["t_zero"])]
    for x in tz:
        a.axvline(x, color="tab:green", lw=.8, alpha=.6)
    a.set(xlabel=r"이완 시작 후 $t/\tau_{int}$", ylabel="결함 수 (z≠6)",
          title=f"③ 결함 이완 — 초록선 = 각 실현이 0 도달 ({len(tz)}/{M}개)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    # ④ 실현 간 산포
    a = ax[1, 1]
    dU = np.array([r_["dU"] for r_ in rows])
    F = np.array([r_["F"] for r_ in rows])
    a.errorbar(dU, F, xerr=[r_["dU_sem"] for r_ in rows],
               yerr=[r_["F_sem"] for r_ in rows], fmt="o", ms=7, capsize=3,
               color="tab:blue", label="개별 런 (막대 = 런 자체 SEM)")
    a.errorbar([dU.mean()], [F.mean()], xerr=[dU.std(ddof=1) / np.sqrt(M)],
               yerr=[F.std(ddof=1) / np.sqrt(M)], fmt="*", ms=20, capsize=4,
               color="crimson", label="앙상블 평균 ±SEM")
    a.axhline(F.mean(), color="crimson", ls=":", lw=1)
    a.axvline(dU.mean(), color="crimson", ls=":", lw=1)
    a.set(xlabel=r"$\Delta U$ [kT/입자]", ylabel=r"$F_{drag}$ [kT/d]",
          title="④ 실현 간 산포 — 런 오차막대가 정직한가")
    a.legend(fontsize=8); a.grid(alpha=.3)

    fig.suptitle(f"trap-drag-2d-hex300 앙상블 (시드 {M}개, 2격자=3.22 d)", fontsize=12)
    fig.tight_layout()
    p = out / "ensemble.png"
    fig.savefig(p, dpi=140)
    plt.close(fig)
    (out / "ensemble.json").write_text(json.dumps(
        {"n_runs": M, "runs": rows,
         "t_rel": t_rel.tolist(), "dU_mean": U.mean(0).tolist(),
         "dU_sem": (U.std(0, ddof=1) / np.sqrt(M)).tolist(),
         "ndef_mean": D.mean(0).tolist()}, indent=2, ensure_ascii=False))
    # ══ 힘 앙상블 ══════════════════════════════════════════════════════
    #   단일 런의 F 한 표본 잡음은 k*·ℓ_k = 246 kT/d 로 신호(~98)보다 크다.
    #   8런을 겹치면 √8 = 2.8배 줄고, **격자 주기로 접으면** 2주기×8런 = 16주기가 되어
    #   17격자 단일런과 같은 통계가 된다.
    Fx, Fy, PHz, Feq = [], [], [], []
    for d in dirs:
        spec, z, r, t, ph = load(d)
        k = spec.params["k_star"]; v = spec.params["v_star"]; a_nn = spec.params["a_nn_star"]
        dm = ph == PH_DRAG
        Fx.append(-k * z["dx_probe"][dm]); Fy.append(-k * z["dy_probe"][dm])
        td = t[dm] - t[dm][0]
        tiv = spec.reduced("times", "tau_int")
        PHz.append(((v * td * tiv) % a_nn) / a_nn)      # 격자 주기 내 위상
        Feq.append(-k * z["dx_probe"][ph == PH_EQ])     # ★ 대조군: 트랩 정지 → ⟨F⟩=0 이어야
    nf = min(len(a) for a in Fx)
    FX = np.array([a[:nf] for a in Fx]); FY = np.array([a[:nf] for a in Fy])
    tdrag = t_rel[:0]  # placeholder
    spec0, z0, r0, t0, ph0 = load(dirs[0])
    td0 = (t0[ph0 == PH_DRAG] - t0[ph0 == PH_DRAG][0])[:nf]
    v_star = spec0.params["v_star"]

    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    a = ax[0, 0]
    w = max(9, nf // 40); kern = np.ones(w) / w
    for f in FX:
        a.plot(td0, np.convolve(f, kern, "same"), "-", lw=.6, color="0.8")
    mu, se = FX.mean(0), FX.std(0, ddof=1) / np.sqrt(M)
    a.plot(td0, np.convolve(mu, kern, "same"), "-", lw=2.2, color="tab:red",
           label=f"앙상블 $F_x$ (M={M}, {w}점 이동평균)")
    a.fill_between(td0, np.convolve(mu - se, kern, "same"),
                   np.convolve(mu + se, kern, "same"), color="tab:red", alpha=.25)
    a.axhline(v_star, color="k", ls="--", lw=1.3, label=f"맨 Stokes γv = {v_star:.1f}")
    a.axhline(FX.mean(), color="darkred", ls=":", lw=1.4,
              label=f"⟨$F_x$⟩ = {FX.mean():.1f} ± {FX.mean(1).std(ddof=1)/np.sqrt(M):.1f}")
    a.set(xlabel=r"끌기 시작 후 $t/\tau_{int}$", ylabel=r"$F_x$ [kT/d]",
          title="① 탐침에 걸리는 힘 — 개별(회색) vs 앙상블")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[0, 1]
    allp = np.concatenate([p_[:nf] for p_ in PHz]); allf = FX.ravel()
    nb = 16; edges = np.linspace(0, 1, nb + 1); ctr = .5 * (edges[:-1] + edges[1:])
    mu2 = np.array([allf[(allp >= l) & (allp < h)].mean() for l, h in zip(edges[:-1], edges[1:])])
    se2 = np.array([allf[(allp >= l) & (allp < h)].std() /
                    max(1, np.sqrt(((allp >= l) & (allp < h)).sum()))
                    for l, h in zip(edges[:-1], edges[1:])])
    a.errorbar(np.r_[ctr, ctr + 1], np.r_[mu2, mu2], yerr=np.r_[se2, se2], fmt="o-",
               ms=5, lw=1.5, capsize=3, color="tab:red",
               label=f"앙상블 ({M}런 × 2주기 = {2*M}주기)")
    a.axhline(np.nanmean(mu2), color="darkred", ls=":", lw=1.3, label=f"평균 {np.nanmean(mu2):.0f}")
    a.axhline(v_star, color="k", ls="--", lw=1.2, label="맨 Stokes")
    a.axvline(1, color="0.7", lw=.8)
    amp2 = (np.nanmax(mu2) - np.nanmin(mu2)) / 2
    a.set(xlabel="격자 주기 내 위상 (x mod $a_{NN}$)/$a_{NN}$", ylabel=r"$F_x$ [kT/d]",
          title=f"② 접은 힘 — 변조 ±{amp2:.0f} kT/d ({100*amp2/abs(np.nanmean(mu2)):.0f}%)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[1, 0]
    muy = FY.mean(0); sey = FY.std(0, ddof=1) / np.sqrt(M)
    a.plot(td0, np.convolve(muy, kern, "same"), "-", lw=2, color="tab:blue", label=r"앙상블 $F_y$")
    a.fill_between(td0, np.convolve(muy - sey, kern, "same"),
                   np.convolve(muy + sey, kern, "same"), color="tab:blue", alpha=.25)
    a.axhline(0, color="k", ls="--", lw=1.2)
    fy_m = FY.mean(); fy_e = FY.mean(1).std(ddof=1) / np.sqrt(M)
    a.set(xlabel=r"끌기 시작 후 $t/\tau_{int}$", ylabel=r"$F_y$ [kT/d]",
          title=f"③ 횡방향 힘 — 0 이어야 한다 (온전성): "
                f"{fy_m:+.1f} ± {fy_e:.1f} = {abs(fy_m)/fy_e:.1f}σ")
    a.legend(fontsize=8); a.grid(alpha=.3)

    a = ax[1, 1]
    FEQ = np.concatenate(Feq)
    for lab, arr, c in (("평형 (트랩 정지)", FEQ, "tab:green"), ("끌기", allf, "tab:red")):
        a.hist(arr, bins=70, histtype="step", lw=1.6, color=c, density=True,
               label=f"{lab}: ⟨F⟩={arr.mean():+.1f}")
    a.axvline(0, color="k", ls="--", lw=1)
    a.axvline(v_star, color="0.4", ls=":", lw=1.4, label=f"맨 Stokes {v_star:.0f}")
    eq_e = FEQ.std() / np.sqrt(len(FEQ))
    a.set(xlabel=r"$F_x$ [kT/d]", ylabel="확률밀도",
          title=f"④ 대조군 — 평형 ⟨F⟩ = {FEQ.mean():+.1f} ± {eq_e:.1f} "
                f"({abs(FEQ.mean())/eq_e:.1f}σ, 0이어야 함)")
    a.legend(fontsize=8); a.grid(alpha=.3)

    fig.suptitle(f"trap-drag 앙상블 — 탐침에 걸리는 힘 (시드 {M}개)", fontsize=12)
    fig.tight_layout()
    pf = out / "ensemble_force.png"; fig.savefig(pf, dpi=140); plt.close(fig)
    print(f"\n힘 앙상블: ⟨F_x⟩ = {FX.mean():.2f} ± {FX.mean(1).std(ddof=1)/np.sqrt(M):.2f} kT/d"
          f"  (맨 Stokes {v_star:.2f} 대비 +{100*(FX.mean()/v_star-1):.0f}%)")
    print(f"  F_y = {fy_m:+.2f} ± {fy_e:.2f} ({abs(fy_m)/fy_e:.1f}σ)   "
          f"평형 대조군 F_x = {FEQ.mean():+.2f} ± {eq_e:.2f} ({abs(FEQ.mean())/eq_e:.1f}σ)")
    print(f"  격자 변조 ±{amp2:.1f} kT/d = 평균의 {100*amp2/abs(np.nanmean(mu2)):.0f}%")
    print(f"\n{p}\n{pf}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
