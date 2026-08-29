"""`trap-drag` 속도 스윕 재분석 — **시드 앙상블 오차막대로**.

## 왜 다시 내는가

기존에 발표된 v-스윕 결과(CLAUDE.md: 항복력 `F(v→0) = 35~45 kT/d`, "0에서 11~20σ",
전단 담화 `F/γv` 16.9→1.87, 결함 수 v-무관)는 **속도마다 단일 런**이고 오차로
**블록 SEM**을 썼습니다. 그런데 시드 앙상블에서 블록 SEM 이 실현 간 산포를
과소평가한다는 것이 측정됐습니다 — γv=24.21 에서 **ΔU/입자 2.35× · F_drag 1.37×**
(`trap_drag_ensemble.py`). ⚠️ 배수는 **관측량마다 다릅니다**. 널리 인용된 2.35 는
ΔU/입자 값이고, 여기서 다루는 끌림힘은 1.37 이었습니다 — 섞어 쓰면 안 됩니다.
어느 쪽이든 기존 유의도(11~20σ)는 과신입니다.

이제 7속도 × 9시드가 모두 있으므로 **실현 간 산포**로 다시 냅니다. 두 오차를
나란히 찍어 과소평가 배수가 속도에 따라 어떻게 변하는지도 봅니다.

⚠️ 이 스크립트는 **수치를 다시 내는 것**이고 물리 해석을 바꾸지 않습니다.
   `f_drag_kT_per_d` 등은 각 런의 `metrics.json["result"]` 에 이미 있는 값을 씁니다 —
   여기서 재계산하면 케이스 코드와 갈라질 수 있습니다.

## 항복력에 대한 태도

`F(v→0)` 은 **외삽**이라 함수형에 의존합니다. 하나만 고르지 않고 셋을 다 냅니다:
최저속 실측값 · 로그선형 외삽 · Herschel–Bulkley 꼴 적합. 셋이 갈리면 그게 결론입니다
(원칙: 단정하지 않고 무엇에 의존하는지 보여준다).

    $PY scratch/trap_drag_vsweep.py
    $PY scratch/trap_drag_vsweep.py --include-legacy   # 옛 코드 런도 섞기 (비권장)
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import nondim as ND  # noqa: E402

OUT = ROOT / "runs" / "trap-drag-2d-hex300__ENSEMBLE"


def is_current(p: Path) -> bool:
    """현재 코드로 돈 런인가 — 판별자는 `step_drift_max_sigma` (L4 측정 배선 이후)."""
    try:
        m = json.loads((p / "metrics.json").read_text())
        return "step_drift_max_sigma" in m.get("numerics", {})
    except Exception:
        return False


def collect(include_legacy: bool) -> dict:
    """속도 → 런 목록. `v_star`(무차원 γv)로 묶습니다 — 태그 문자열보다 신뢰할 수 있습니다."""
    by_v: dict[float, list] = {}
    skipped = []
    for p in sorted((ROOT / "runs").glob("trap-drag-2d-hex300__tr0.117647*")):
        if not (p / "metrics.json").exists() or not (p / "spec.json").exists():
            continue
        if "smoke" in p.name:
            continue
        if not include_legacy and not is_current(p):
            skipped.append(p.name)
            continue
        try:
            spec = ND.load(p / "spec.json")
            m = json.loads((p / "metrics.json").read_text())
            r = m["result"]
        except Exception as e:
            skipped.append(f"{p.name} (읽기 실패: {e})")
            continue
        # 예열/평형/이완 구간 설정이 기본과 다른 변형 런은 스윕에서 제외 (w2e3r10 등)
        tag = p.name.split("__")[1]
        if "-w" in tag:
            skipped.append(f"{p.name} (구간 설정 변형)")
            continue
        by_v.setdefault(float(spec.params["v_star"]), []).append({
            "name": p.name,
            "seed": int(spec.numerics["seed"]),
            "v_star": float(spec.params["v_star"]),
            "F": float(r["f_drag_kT_per_d"]),
            "F_blocksem": float(r["f_drag_sem"]),
            "stokes": float(r["f_stokes_bare"]),
            "n_def": float(r["n_def"]["driven"]),
            "psi6": float(r["psi6"]["driven"]),
            "rec": float(r["relax_fit_defects"]["recovered_frac"]),
            "dev": float(r["dev_max_d"]),
        })
    return by_v, skipped


def stats(xs):
    a = np.asarray([x for x in xs if np.isfinite(x)], float)
    if a.size == 0:
        return dict(mean=np.nan, std=np.nan, sem=np.nan, n=0)
    return dict(mean=a.mean(), std=a.std(ddof=1) if a.size > 1 else np.nan,
                sem=(a.std(ddof=1) / np.sqrt(a.size)) if a.size > 1 else np.nan, n=a.size)


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--include-legacy", action="store_true")
    a = ap.parse_args()

    by_v, skipped = collect(a.include_legacy)
    if skipped:
        print(f"⚠ 제외 {len(skipped)}런 (옛 코드 또는 변형 설정) — 조용히 버리지 않습니다:")
        for s in skipped[:10]:
            print(f"    - {s}")
        if len(skipped) > 10:
            print(f"    … 외 {len(skipped)-10}건")
    if len(by_v) < 3:
        print(f"속도 점이 부족합니다 ({len(by_v)}개)", file=sys.stderr)
        return 1

    vs = sorted(by_v)
    print("=" * 104)
    print(f"trap-drag v-스윕 재분석 — 속도 {len(vs)}점 · 시드 앙상블")
    print("=" * 104)
    print(f"{'γv':>8}{'시드':>5}{'⟨F⟩':>9}{'앙상블SEM':>11}{'블록SEM':>10}"
          f"{'과소평가':>9}{'F/γv':>8}{'결함':>7}{'ψ6':>7}{'회복%':>7}")
    rows = []
    for v in vs:
        g = by_v[v]
        F = stats([x["F"] for x in g])
        blk = float(np.mean([x["F_blocksem"] for x in g]))
        ratio = F["std"] / blk if blk else np.nan
        nd = stats([x["n_def"] for x in g])
        ps = stats([x["psi6"] for x in g])
        rc = stats([x["rec"] for x in g])
        rows.append(dict(v=v, F=F, blk=blk, ratio=ratio, nd=nd, psi6=ps, rec=rc,
                         thin=F["mean"] / v))
        print(f"{v:>8.2f}{F['n']:>5}{F['mean']:>9.1f}{F['sem']:>11.1f}{blk:>10.1f}"
              f"{ratio:>9.2f}{F['mean']/v:>8.2f}{nd['mean']:>7.2f}{ps['mean']:>7.3f}"
              f"{100*rc['mean']:>7.1f}")

    # ── 블록 SEM 과소평가 배수 ─────────────────────────────────────────
    rr = np.array([r["ratio"] for r in rows], float)
    rr = rr[np.isfinite(rr)]
    print("-" * 104)
    # ⚠️ 비교 기준을 관측량별로 정확히: 단일-v 앙상블(γv=24.21)의 과소평가 배수는
    #    **ΔU/입자 2.35× · F_drag 1.37×** 였습니다. 여기서 재는 것은 F_drag 이므로
    #    1.37 과 비교해야 합니다. 2.35 와 비교하면 다른 관측량을 섞는 것입니다.
    print(f"블록 SEM 과소평가 배수 (F_drag): 중앙 {np.median(rr):.2f}× · "
          f"범위 {rr.min():.2f}~{rr.max():.2f}×")
    print(f"  대조 — 단일-v 앙상블(γv=24.21): F_drag 1.37× · ΔU/입자 2.35× "
          f"(관측량마다 다르다)")

    # ── 항복력 — 외삽 3종 ─────────────────────────────────────────────
    V = np.array([r["v"] for r in rows])
    FM = np.array([r["F"]["mean"] for r in rows])
    FE = np.array([r["F"]["sem"] for r in rows])
    print()
    print("항복력 F(v→0) — **외삽이라 함수형에 의존합니다. 셋을 다 냅니다**")
    lo = int(np.argmin(V))
    print(f"  ⓐ 최저속 실측      γv={V[lo]:.2f} 에서 F = {FM[lo]:.1f} ± {FE[lo]:.1f} kT/d"
          f"   → 0 에서 {FM[lo]/FE[lo]:.1f}σ   (외삽 없음, 가장 방어 가능)")

    # ⓑ 로그선형: F = a + b·ln(γv)  — 저속 3점
    k = min(3, len(V))
    idx = np.argsort(V)[:k]
    A = np.vstack([np.ones(k), np.log(V[idx])]).T
    w = 1.0 / np.where(FE[idx] > 0, FE[idx], np.nan)
    good = np.isfinite(w)
    if good.sum() >= 2:
        coef, *_ = np.linalg.lstsq(A[good] * w[good, None], FM[idx][good] * w[good], rcond=None)
        # v→0 은 ln→-inf 로 발산 → 이 꼴은 항복력을 정의하지 못한다는 것이 결론
        print(f"  ⓑ 로그선형 적합    F = {coef[0]:.1f} + {coef[1]:.1f}·ln(γv)  "
              f"→ v→0 에서 **발산**(정의 안 됨). 로그 꼴은 항복력을 주지 못합니다")

    # ⓒ Herschel–Bulkley: F = F_y + c·(γv)^n  — 전 구간, n 을 격자로 탐색
    best = None
    for n_ in np.linspace(0.1, 1.5, 141):
        X = np.vstack([np.ones(len(V)), V ** n_]).T
        ww = 1.0 / np.where(FE > 0, FE, np.nan)
        m_ = np.isfinite(ww)
        if m_.sum() < 3:
            continue
        c, res, *_ = np.linalg.lstsq(X[m_] * ww[m_, None], FM[m_] * ww[m_], rcond=None)
        pred = X @ c
        chi2 = float(np.nansum(((FM - pred) * ww) ** 2))
        if best is None or chi2 < best[0]:
            best = (chi2, n_, c)
    if best:
        chi2, n_, c = best
        dof = max(1, int(np.isfinite(FE).sum()) - 3)
        # F_y 불확실도: 가중 최소제곱 공분산
        X = np.vstack([np.ones(len(V)), V ** n_]).T
        ww = 1.0 / np.where(FE > 0, FE, np.nan)
        m_ = np.isfinite(ww)
        XtX = (X[m_] * ww[m_, None]).T @ (X[m_] * ww[m_, None])
        cov = np.linalg.inv(XtX) * max(1.0, chi2 / dof)
        sF = float(np.sqrt(cov[0, 0]))
        print(f"  ⓒ Herschel–Bulkley  F = F_y + c·(γv)^n,  n={n_:.2f} 최적 "
              f"→ F_y = {c[0]:.1f} ± {sF:.1f} kT/d   → 0 에서 {abs(c[0])/sF:.1f}σ"
              f"   (χ²/dof = {chi2/dof:.2f})")

    print()
    print("전단 담화 F/γv")
    print(f"  γv={V[lo]:.2f}: {FM[lo]/V[lo]:.2f}  →  γv={V[-1]:.2f}: {FM[-1]/V[-1]:.2f}"
          f"   ({(FM[lo]/V[lo])/(FM[-1]/V[-1]):.1f}배 감소)")

    nd_m = np.array([r["nd"]["mean"] for r in rows])
    nd_e = np.array([r["nd"]["sem"] for r in rows])
    print()
    print("결함 수 — v 에 의존하는가")
    print(f"  범위 {nd_m.min():.2f}~{nd_m.max():.2f} · 앙상블 SEM 평균 {np.nanmean(nd_e):.2f}")
    spread = (nd_m.max() - nd_m.min()) / np.nanmean(nd_e)
    print(f"  최대-최소 = {spread:.1f}·SEM → "
          + ("**v 의존성 있음**" if spread > 3 else "v 무관과 구별 안 됨"))

    _plot(rows, V, FM, FE, best)
    _save(rows, best)
    return 0


def _plot(rows, V, FM, FE, best) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    # ★ 한글 + '−'(U+2212) 를 모두 갖춘 것은 Arial Unicode MS 뿐 (CLAUDE.md)
    matplotlib.rcParams["font.family"] = ["Arial Unicode MS", "AppleGothic", "DejaVu Sans"]
    matplotlib.rcParams["axes.unicode_minus"] = False

    OUT.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(2, 2, figsize=(13, 9))
    blk = np.array([r["blk"] for r in rows])

    a = ax[0, 0]
    a.errorbar(V, FM, yerr=FE, fmt="o-", capsize=4, lw=1.8, color="tab:red",
               label="앙상블 SEM (시드 간)")
    a.errorbar(V, FM, yerr=blk, fmt="none", capsize=8, lw=1.0, color="tab:blue",
               alpha=.8, label="블록 SEM (단일 런) — 과소평가")
    a.plot(V, V, "k--", lw=1.2, label="맨 Stokes γv")
    if best:
        _, n_, c = best
        vv = np.geomspace(V.min() * .5, V.max() * 1.2, 100)
        a.plot(vv, c[0] + c[1] * vv ** n_, ":", color="darkgreen", lw=1.6,
               label=f"HB 적합 (n={n_:.2f}, $F_y$={c[0]:.0f})")
        a.axhline(c[0], color="darkgreen", ls="-.", lw=1.0, alpha=.6)
    a.set(xscale="log", xlabel=r"$\gamma v$ [kT/d]", ylabel=r"$\langle F_x\rangle$ [kT/d]",
          title="① 끌림힘 vs 속도 — 두 오차 비교")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[0, 1]
    a.plot(V, FM / V, "o-", lw=1.8, color="tab:purple")
    a.axhline(1, color="k", ls="--", lw=1.2, label="맨 Stokes (=1)")
    a.set(xscale="log", yscale="log", xlabel=r"$\gamma v$ [kT/d]",
          ylabel=r"$F/\gamma v$", title="② 전단 담화 — 유효 항력 배수")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    a = ax[1, 0]
    nd = np.array([r["nd"]["mean"] for r in rows])
    nde = np.array([r["nd"]["sem"] for r in rows])
    a.errorbar(V, nd, yerr=nde, fmt="s-", capsize=4, lw=1.8, color="tab:orange")
    a.set(xscale="log", xlabel=r"$\gamma v$ [kT/d]", ylabel="구동 중 결함 수",
          title="③ 결함 수 — v 의존성")
    a.grid(alpha=.3, which="both")

    a = ax[1, 1]
    ratio = np.array([r["ratio"] for r in rows])
    a.plot(V, ratio, "D-", lw=1.8, color="tab:green")
    a.axhline(1, color="k", ls="--", lw=1.2, label="일치 (=1)")
    a.axhline(1.37, color="tab:red", ls=":", lw=1.4,
              label=r"단일-v 측정 $F_{drag}$ 1.37×")
    a.set(xscale="log", xlabel=r"$\gamma v$ [kT/d]",
          ylabel="실현간 std / 블록 SEM", title="④ 블록 SEM 과소평가 배수")
    a.legend(fontsize=8); a.grid(alpha=.3, which="both")

    fig.tight_layout()
    p = OUT / "vsweep_ensemble.png"
    fig.savefig(p, dpi=130)
    print(f"\n그림: {p}")


def _save(rows, best) -> None:
    out = {"schema": "bdbot.trap_drag_vsweep/0.2",
           "note": "시드 앙상블 오차막대. 오차는 실현 간 산포(앙상블 SEM)이며 "
                   "블록 SEM 이 아니다.",
           "points": [{"v_star": r["v"], "n_seeds": r["F"]["n"],
                       "F_mean": r["F"]["mean"], "F_ensemble_sem": r["F"]["sem"],
                       "F_realization_std": r["F"]["std"], "F_block_sem_mean": r["blk"],
                       "block_sem_underestimate": r["ratio"],
                       "F_over_gamma_v": r["thin"],
                       "n_def_mean": r["nd"]["mean"], "n_def_sem": r["nd"]["sem"],
                       "psi6_mean": r["psi6"]["mean"],
                       "recovered_frac_mean": r["rec"]["mean"]} for r in rows]}
    if best:
        chi2, n_, c = best
        out["herschel_bulkley"] = {"n": n_, "F_yield": c[0], "c": c[1]}
    p = OUT / "vsweep_ensemble.json"
    p.write_text(json.dumps(out, indent=2, ensure_ascii=False))
    print(f"데이터: {p}")


if __name__ == "__main__":
    sys.exit(main())
