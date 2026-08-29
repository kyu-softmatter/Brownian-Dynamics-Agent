"""L3 스펙 결함 4건의 **회귀 테스트**.

CLAUDE.md 규칙 6: 물리·도구 동작 주장은 추론이 아니라 실행으로 확인한다.
이 스크립트는 처음에 결함을 **재현**하려고 썼고 4건 모두 재현됐습니다 (2026-08-04):

  ① run_id 가 물리계를 덮지 않았다 — 1-B 스펙에 물리계가 없어서 d 5µm→0.5µm,
     η 물→글리세롤(62배), ρ_p 실리카→폴리스티렌 으로 바꿔도 run_id 가
     `soft-r3-2d-A-sweep__A100__27f70deab9` 로 **완전히 같았다** (τ_B 16.1배 차이).
     그건 이미 완료된 런의 이름이라 `prepare_outdir` 가 건너뛰고 **예전 계의 결과를
     새 계의 결과로 보고**한다.
  ② 스펙만으로 역변환이 불가능했다 (σ·τ·kT 의 SI 값이 없음).
  ③ 세 케이스의 spec.json 스키마가 달랐다 (공통 키가 N·n_eq·n_prod 3개뿐).
  ④ 무차원수가 원장 두 항목의 비인지 검사할 수 없었다 → `verify_nondim_guards.py`

`bdbot.nondim.NondimSpec` 도입 후 ①②③ 이 사라졌는지 확인합니다.
**이제는 통과(=결함 없음)가 정상입니다.**

    $PY scratch/verify_l3_spec_gaps.py
"""
from __future__ import annotations

import copy
import math
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from bdbot import nondim as ND  # noqa: E402

sys.path.insert(0, str(ROOT / "cases"))
import soft_r3_2d as S3  # noqa: E402

FAILS = []
LEGACY_RUN_ID = "soft-r3-2d-A-sweep__A100__27f70deab9"   # 결함 ① 당시의 값


def check(name: str, ok: bool, detail: str) -> None:
    print(f"  {'✓' if ok else '✗'} {name}\n      {detail}")
    if not ok:
        FAILS.append(name)


def soft_spec(raw: dict, A=100.0, N=400, rc_shells=5.0, dt_scale=1.0, samples=400):
    """`cases/soft_r3_2d.py` main() 의 **현재** 스펙 구성 경로를 그대로 재현."""
    tmp = ROOT / "scratch" / "_tmp_sys.yaml"
    tmp.write_text(yaml.safe_dump(raw))
    sys_ = S3.load_system(tmp)
    tmp.unlink()

    phi = sys_["phi"]
    a_star = math.sqrt(math.pi / (4 * phi))
    r_c_star = rc_shells * a_star
    T_obs_tau = float(sys_["numerics"]["production_tau_B"])

    lg = S3.build_ledger(sys_, A, N, phi, r_c_star, dt_scale, T_obs_tau)
    tau_B = lg.derived["tau_B"]
    dt, T_obs = lg.get("times", "dt"), lg.get("times", "T_obs")
    n_eq = int(round(float((0.2 * T_obs / dt).to(""))))
    n_prod = int(round(float((T_obs / dt).to(""))))
    sample_every = max(1, n_prod // samples)
    n_prod = (n_prod // sample_every) * sample_every

    groups, checks, Gamma = S3.analyze_scales(sys_, lg, A, phi, r_c_star)
    spec = ND.NondimSpec(
        case=sys_["label"], system=sys_["_raw"], reference=lg.ref, ledger=lg,
        groups=groups, checks=checks,
        params={"A": A, "phi": phi, "N": N, "r_c_star": r_c_star,
                "wca_eps": sys_["wca_eps_kT"], "Gamma": Gamma},
        numerics={"dt_star": float((dt / tau_B).to("")),
                  "dt_over_tau_int": dt_scale * 1e-2, "n_eq": n_eq, "n_prod": n_prod,
                  "n_samples": samples, "sample_every": sample_every, "seed": 20260803},
        tag=f"A{A:g}", nhex=10)
    return spec, dt.to("s"), tau_B.to("s")


def main() -> int:
    print("=" * 78)
    print("L3 스펙 결함 회귀 테스트 — soft-r3-2d-A-sweep")
    print("=" * 78)

    raw = yaml.safe_load((ROOT / "intake" / "soft-r3-2d-A-sweep" / "system.yaml").read_text())
    spec_a, dt_a, tauB_a = soft_spec(raw)
    id_a = spec_a.run_id()

    other = copy.deepcopy(raw)
    other["particle"]["diameter"]["value"] = 0.5
    other["medium"]["viscosity"]["value"] = 53.0
    other["medium"]["name"] = "glycerol-water 60%"
    other["particle"]["density"]["value"] = 1050
    spec_b, dt_b, tauB_b = soft_spec(other)
    id_b = spec_b.run_id()

    # ── ① run_id 가 물리계를 덮는가 ──────────────────────────────────────
    print("\n[①] run_id 가 물리계(d·η·ρ_p)를 덮는가")
    print(f"    계 A: d=5.0µm  η=0.851mPa·s  →  τ_B={tauB_a:~.4gP}  dt={dt_a.to('ms'):~.4gP}")
    print(f"    계 B: d=0.5µm  η=53mPa·s     →  τ_B={tauB_b:~.4gP}  dt={dt_b.to('ms'):~.4gP}")
    print(f"    물리 시간 배율: τ_B 가 {float(tauB_a/tauB_b):.1f}배 다름")
    check("물리계가 다르면 run_id 가 다르다", id_a != id_b, f"A={id_a}  B={id_b}")
    check("스펙이 물리계를 기록한다", bool(spec_a.to_json().get("system")),
          f"system.label = {spec_a.to_json()['system'].get('label')}")
    check("결함 당시의 run_id 를 더 이상 만들지 않는다",
          LEGACY_RUN_ID not in (id_a, id_b),
          f"레거시 {LEGACY_RUN_ID} → 현재 {id_a} (기존 런 7개는 레거시로 보존)")

    print("\n[①b] 문서 필드는 run_id 를 바꾸지 않는다 (콘텐츠 주소의 조건)")
    doc = copy.deepcopy(raw)
    doc["particle"]["diameter"]["source"] = "출처 문구만 고침"
    doc["description"] = "설명만 고침"
    check("출처·설명을 고쳐도 run_id 유지", soft_spec(doc)[0].run_id() == id_a,
          f"{soft_spec(doc)[0].run_id()} vs {id_a}")

    # ── ② 역변환 3앵커 ──────────────────────────────────────────────────
    print("\n[②] 스펙만으로 결과를 물리 단위로 되돌릴 수 있는가")
    bt = spec_a.to_json()["back_transform"]
    need = ("sigma_SI", "tau_SI", "energy_SI")
    check("역변환 3앵커가 스펙에 있다", all(k in bt for k in need),
          "  ".join(f"{k}={bt[k]:.4g} {bt[k+'_unit']}" for k in need if k in bt))
    d_si = float(raw["particle"]["diameter"]["value"]) * 1e-6
    check("sigma_SI 가 실제 입자 지름과 일치", abs(bt["sigma_SI"] - d_si) / d_si < 1e-12,
          f"{bt['sigma_SI']:.6e} m vs {d_si:.6e} m")

    # ── ③ 세 케이스 공통 스키마 ─────────────────────────────────────────
    print("\n[③] 세 케이스가 같은 스키마를 쓰는가")
    keysets = {}
    for name, script, extra in (("trap-2d-5um", "trap_2d_5um", []),
                                ("soft-r3", "soft_r3_2d", ["--A", "100"]),
                                ("abp-rod", "abp_rod_2d", [])):
        import subprocess
        out = subprocess.run([sys.executable, str(ROOT / "cases" / f"{script}.py"),
                              *extra, "--spec"], capture_output=True, text=True, cwd=ROOT)
        line = [x for x in out.stdout.splitlines() if "L3 스펙:" in x]
        if not line:
            check(f"{name} 이 스펙을 쓴다", False, out.stdout[-300:] + out.stderr[-300:])
            keysets[name] = set()
            continue
        p = ROOT / line[0].split("L3 스펙:")[1].strip()
        ls = ND.load(p)
        keysets[name] = set(ls.raw)
        ok, want = ls.verify_hash()
        print(f"    {name:<12} {p.name}   해시검증 {'✓' if ok else '✗ ' + want}")
        if not ok:
            FAILS.append(f"{name} 해시검증")
    common = set.intersection(*keysets.values()) if all(keysets.values()) else set()
    required = {"schema", "reference", "back_transform", "ledger", "groups", "checks",
                "system", "params", "numerics", "run_id", "verdict"}
    check("세 케이스가 공통 스키마를 공유한다", required <= common,
          f"공통 키 {len(common)}개 · 누락 {sorted(required - common) or '없음'}")

    print("\n" + "=" * 78)
    if FAILS:
        print(f"✗ FAIL — 결함 {len(FAILS)}건:")
        for f in FAILS:
            print(f"   ✗ {f}")
        print("=" * 78)
        return 1
    print("✓ PASS — 결함 ①②③ 이 모두 해소됐습니다.")
    print("=" * 78)
    return 0


if __name__ == "__main__":
    sys.exit(main())
