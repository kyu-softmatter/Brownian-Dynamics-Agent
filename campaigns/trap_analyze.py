"""S6 + S7 — 트랩 배치의 시각화와 분석.

usage: python scripts/trap_analyze.py <run_dir>
"""
from __future__ import annotations

import json
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(REPO))

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

from simbot.analysis.trap import (
    aggregate_seeds, check_position_distribution, em_uniform_noise_excess_kurtosis,
    fit_msd, load_run, msd_model,
)
from simbot.estimators import euler_maruyama_trap_variance_bias, harmonic_trap

SKETCH = dict(T_si=300.0, eta_si=8.5566e-4, radius_si=5e-6, k_si=1e-5)
DTS = (5.0e-3, 2.5e-3)
DIMS = (2, 3)
SEEDS = (1, 2, 3, 4)

run_dir = Path(sys.argv[1])
figs = run_dir / "figs"
figs.mkdir(parents=True, exist_ok=True)

# --- 로드 ---------------------------------------------------------------
runs: dict[tuple[int, float, int], dict] = {}
for d in DIMS:
    for dt in DTS:
        for s in SEEDS:
            p = run_dir / "raw" / f"d{d}_dt{dt:g}_s{s}"
            if p.exists():
                runs[(d, dt, s)] = load_run(p)

print(f"# 로드 {len(runs)} 런\n")

# =========================================================================
# 1. <x*^2> — 시드 앙상블 + EM 편향 대조
# =========================================================================
print("=" * 76)
print("P1/P2/P8  <x*^2> — 등분배와 Euler-Maruyama 편향")
print("=" * 76)
print(f"{'dim':>4} {'dt*':>8} {'측정':>10} {'±SE':>8} {'시드산포':>9} | "
      f"{'EM예측':>9} {'exact':>7} | {'EM에서':>8} {'판별력':>7}")
print("-" * 76)

var_summary = {}
for d in DIMS:
    for dt in DTS:
        vals = [float(runs[(d, dt, s)]["indep_var"].mean()) for s in SEEDS
                if (d, dt, s) in runs]
        agg = aggregate_seeds(vals)
        bias = euler_maruyama_trap_variance_bias(dt)
        em = 1.0 + bias
        dev_em = abs(agg.mean - em) / agg.se
        power = bias / agg.se
        var_summary[(d, dt)] = (agg, em, dev_em, power)
        print(f"{d:>4} {dt:>8.1e} {agg.mean:>10.5f} {agg.se:>8.5f} {agg.spread:>9.5f} | "
              f"{em:>9.5f} {1.0:>7.3f} | {dev_em:>7.2f}σ {power:>6.2f}σ")

# =========================================================================
# 2. MSD 피팅
# =========================================================================
print()
print("=" * 76)
print("P4/P5/P6  MSD = 2d(1 - exp(-t/tau))  — plateau 와 완화시간")
print("=" * 76)
print(f"{'dim':>4} {'dt*':>8} | {'plateau':>9} {'±':>7} {'예측':>6} | "
      f"{'tau':>8} {'±':>7} {'예측':>5} | {'R^2':>8}")
print("-" * 76)

msd_summary = {}
for d in DIMS:
    for dt in DTS:
        pls, taus, r2s, curves = [], [], [], []
        for s in SEEDS:
            if (d, dt, s) not in runs:
                continue
            r = runs[(d, dt, s)]
            lags_tau = r["lags_steps"] * dt
            f = fit_msd(lags_tau, r["msd"], d)
            pls.append(f.plateau); taus.append(f.tau); r2s.append(f.r_squared)
            curves.append((lags_tau, r["msd"]))
        ap, at = aggregate_seeds(pls), aggregate_seeds(taus)
        msd_summary[(d, dt)] = (ap, at, float(np.mean(r2s)), curves)
        print(f"{d:>4} {dt:>8.1e} | {ap.mean:>9.4f} {ap.se:>7.4f} {2*d:>6d} | "
              f"{at.mean:>8.4f} {at.se:>7.4f} {1.0:>5.1f} | {np.mean(r2s):>8.6f}")

# =========================================================================
# 3. 위치 분포
# =========================================================================
print()
print("=" * 76)
print("P7  위치 분포 — 균일 노이즈가 CLT 로 Gaussian 이 되는가")
print("=" * 76)
print(f"{'dim':>4} {'dt*':>8} | {'독립프레임':>9} {'n_eff':>7} {'KS':>8} {'p값':>7} | "
      f"{'첨도':>7} {'예측':>7} {'±SE':>6} {'편차':>6}")
print("-" * 76)
dist_summary = {}
for d in DIMS:
    for dt in DTS:
        ks, ps, ku, nf, ne = [], [], [], 0, 0
        for s in SEEDS:
            if (d, dt, s) not in runs:
                continue
            r = runs[(d, dt, s)]
            stride = int(r["lags_steps"][0])
            c = check_position_distribution(r["traj"], dt_star=dt,
                                            frame_interval_steps=stride)
            ks.append(c.ks_stat); ps.append(c.ks_p); ku.append(c.kurtosis)
            nf, ne, kpred, kse = c.n_independent_frames, c.n_effective_samples, \
                c.kurtosis_predicted, c.kurtosis_se
        kmean = float(np.mean(ku))
        # 시드 4개 평균이므로 SE 는 sqrt(4) 만큼 줄어든다
        kse_agg = kse / np.sqrt(len(ku))
        dev = abs(kmean - kpred) / kse_agg
        dist_summary[(d, dt)] = (np.mean(ks), np.min(ps), kmean, kpred, kse_agg, dev)
        print(f"{d:>4} {dt:>8.1e} | {nf:>9d} {ne:>7d} {np.mean(ks):>8.5f} "
              f"{np.min(ps):>7.4f} | {kmean:>7.4f} {kpred:>7.4f} {kse_agg:>6.4f} "
              f"{dev:>5.2f}σ")

# =========================================================================
# 4. 물리 단위 환산
# =========================================================================
print()
print("=" * 76)
print("물리 단위 환산 (봉인된 02_prediction.md 대조)")
print("=" * 76)
phys = {}
for d in DIMS:
    p = harmonic_trap(dim=d, **SKETCH)
    agg = var_summary[(d, 5e-3)][0]
    var_x_nm2 = agg.mean * p.var_per_component_si * 1e18
    se_x_nm2 = agg.se * p.var_per_component_si * 1e18
    var_r_nm2 = d * var_x_nm2
    plateau_nm2 = msd_summary[(d, 5e-3)][0].mean * p.var_per_component_si * 1e18
    tau_ms = msd_summary[(d, 5e-3)][1].mean * p.tau_trap_si * 1e3
    tau_se_ms = msd_summary[(d, 5e-3)][1].se * p.tau_trap_si * 1e3
    phys[d] = dict(var_x=var_x_nm2, se_x=se_x_nm2, var_r=var_r_nm2,
                   plateau=plateau_nm2, tau=tau_ms, tau_se=tau_se_ms)
    print(f"  dim={d}:  <x^2> = {var_x_nm2:7.2f} ± {se_x_nm2:.2f} nm^2 "
          f"(예측 {p.var_per_component_si*1e18:.2f})")
    print(f"          <r^2> = {var_r_nm2:7.2f} nm^2 "
          f"(예측 {d*p.var_per_component_si*1e18:.2f})")
    print(f"          MSD plateau = {plateau_nm2:8.2f} nm^2 "
          f"(예측 {p.msd_plateau_si*1e18:.2f})")
    print(f"          tau_trap = {tau_ms:.4f} ± {tau_se_ms:.4f} ms "
          f"(예측 {p.tau_trap_si*1e3:.4f})")

# =========================================================================
# 그림  (글자는 영문 — matplotlib 한글 글리프 없음)
# =========================================================================
plt.rcParams.update({"font.size": 9, "figure.dpi": 130,
                     "axes.grid": True, "grid.alpha": 0.25})

# --- Fig 1: MSD ---
fig, axes = plt.subplots(1, 2, figsize=(9, 3.6))
for ax, d in zip(axes, DIMS):
    _, _, _, curves = msd_summary[(d, 5e-3)]
    for i, (t, m) in enumerate(curves):
        ax.loglog(t, m, "o", ms=2.5, alpha=0.5,
                  label="simulation (4 seeds)" if i == 0 else None)
    tt = np.geomspace(curves[0][0][0], curves[0][0][-1], 200)
    ax.loglog(tt, msd_model(tt, 2 * d, 1.0), "k-", lw=1.6,
              label=r"analytic $2d(1-e^{-t/\tau})$")
    ax.loglog(tt, 2 * d * tt, "r--", lw=1.0, alpha=0.7,
              label=r"free diffusion $2dD_0t$")
    ax.axhline(2 * d, color="gray", ls=":", lw=1)
    ax.set_xlabel(r"lag  $t/\tau_{\rm trap}$")
    ax.set_ylabel(r"MSD  $\langle\Delta r^{*2}\rangle$")
    ax.set_title(f"{d}D   plateau $=2d={2*d}$")
    ax.legend(fontsize=7, loc="lower right")
fig.suptitle("Harmonic trap MSD — reduced units ($\\ell_{\\rm trap}$, $kT$, $\\tau_{\\rm trap}$)",
             fontsize=10)
fig.tight_layout()
fig.savefig(figs / "01_msd.png", bbox_inches="tight")
plt.close(fig)

# --- Fig 2: <x^2> vs dt* with EM line ---
fig, ax = plt.subplots(figsize=(4.6, 3.4))
dts = np.geomspace(1e-3, 3e-2, 100)
ax.plot(dts, 1 + np.array([euler_maruyama_trap_variance_bias(x) for x in dts]),
        "k-", lw=1.5, label=r"Euler-Maruyama  $1/(1-\Delta t^*/2)$")
ax.axhline(1.0, color="r", ls="--", lw=1.2, label="exact scheme (no bias)")
for d, mk, col in zip(DIMS, ("o", "s"), ("C0", "C1")):
    xs = [dt for dt in DTS]
    ys = [var_summary[(d, dt)][0].mean for dt in DTS]
    es = [var_summary[(d, dt)][0].se for dt in DTS]
    ax.errorbar(xs, ys, yerr=es, fmt=mk, color=col, capsize=3, ms=5,
                label=f"{d}D measured (4 seeds)")
ax.set_xscale("log")
ax.set_xlabel(r"$\Delta t^* = \Delta t/\tau_{\rm trap}$")
ax.set_ylabel(r"$\langle x^{*2}\rangle$")
ax.set_title("Equipartition and integrator bias")
ax.legend(fontsize=7)
fig.tight_layout()
fig.savefig(figs / "02_equipartition_dt.png", bbox_inches="tight")
plt.close(fig)

# --- Fig 3: position distribution ---
fig, axes = plt.subplots(1, 2, figsize=(9, 3.4))
for ax, d in zip(axes, DIMS):
    tr = runs[(d, 5e-3, 1)]["traj"][-20:].reshape(-1)
    ax.hist(tr, bins=120, density=True, alpha=0.6, color=f"C{d-2}",
            label="simulation")
    xs = np.linspace(-4.5, 4.5, 300)
    ax.plot(xs, np.exp(-xs**2 / 2) / np.sqrt(2 * np.pi), "k-", lw=1.5,
            label=r"$\mathcal{N}(0,1)$ (Boltzmann)")
    ax.set_xlabel(r"$x/\ell_{\rm trap}$")
    ax.set_ylabel("PDF")
    ax.set_yscale("log")
    ax.set_ylim(1e-5, 1)
    ax.set_title(f"{d}D   component positions")
    ax.legend(fontsize=7)
fig.suptitle("Position distribution — uniform noise becomes Gaussian by CLT",
             fontsize=10)
fig.tight_layout()
fig.savefig(figs / "03_distribution.png", bbox_inches="tight")
plt.close(fig)

# --- Fig 4: 2D vs 3D radial — the A1 discriminator ---
fig, ax = plt.subplots(figsize=(4.6, 3.4))
p2 = harmonic_trap(dim=2, **SKETCH)
for d, col in zip(DIMS, ("C0", "C1")):
    tr = runs[(d, 5e-3, 1)]["traj"][-20:]
    r = np.sqrt(np.sum(tr.astype(np.float64) ** 2, axis=2)).reshape(-1)
    r_nm = r * p2.l_trap_si * 1e9
    ax.hist(r_nm, bins=100, density=True, alpha=0.55, color=col,
            label=f"{d}D  " + r"$\langle r^2\rangle$ = "
                  f"{phys[d]['var_r']:.0f} nm$^2$")
ax.set_xlabel("radial displacement  $r$  [nm]")
ax.set_ylabel("PDF")
ax.set_title("A1 discriminator: 2D vs 3D radial distribution")
ax.legend(fontsize=8)
fig.tight_layout()
fig.savefig(figs / "04_A1_discriminator.png", bbox_inches="tight")
plt.close(fig)

print()
print(f"# 그림 4개 저장: {figs}")

# --- 요약 JSON ---
out = {
    "var_x_star": {f"d{d}_dt{dt:g}": dict(
        mean=var_summary[(d, dt)][0].mean, se=var_summary[(d, dt)][0].se,
        spread=var_summary[(d, dt)][0].spread,
        em_prediction=var_summary[(d, dt)][1],
        sigma_from_em=var_summary[(d, dt)][2],
        discrimination_power=var_summary[(d, dt)][3])
        for d in DIMS for dt in DTS},
    "msd": {f"d{d}_dt{dt:g}": dict(
        plateau=msd_summary[(d, dt)][0].mean, plateau_se=msd_summary[(d, dt)][0].se,
        plateau_predicted=2 * d,
        tau=msd_summary[(d, dt)][1].mean, tau_se=msd_summary[(d, dt)][1].se,
        r_squared=msd_summary[(d, dt)][2])
        for d in DIMS for dt in DTS},
    "distribution": {f"d{d}_dt{dt:g}": dict(
        ks_stat=dist_summary[(d, dt)][0], ks_p_min=dist_summary[(d, dt)][1],
        kurtosis=dist_summary[(d, dt)][2],
        kurtosis_predicted=dist_summary[(d, dt)][3],
        kurtosis_se=dist_summary[(d, dt)][4],
        kurtosis_sigma=dist_summary[(d, dt)][5])
        for d in DIMS for dt in DTS},
    "physical_units": phys,
}
(run_dir / "metrics.json").write_text(json.dumps(out, indent=2, default=float))
print(f"# metrics.json 저장")
