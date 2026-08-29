"""`bdbot` CLI — 앞단(L0→L3)의 단일 진입점.

마스터플랜 §15.10: Claude Code 세션 **밖에서도** 똑같이 동작해야 합니다. cron·스크립트·
다른 사람이 같은 명령으로 같은 결과를 얻어야 하고, Phase 8-rest의 훅이 가로챌 대상도
이 명령들입니다.

    PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python

    $PY -m bdbot.cli status                     모든 케이스의 파이프라인 진행도
    $PY -m bdbot.cli intake init  <folder>      observation.yaml 템플릿
    $PY -m bdbot.cli intake check <folder>      L0 스키마 + 준비도 판정
    $PY -m bdbot.cli system check <folder>      L2 스키마 + tier + 유도값 재계산
    $PY -m bdbot.cli nondim report <case>       L3 무차원화 리포트 (실행 안 함)
    $PY -m bdbot.cli nondim spec   <case>       L3 스펙 → specs/<run_id>.json
    $PY -m bdbot.cli nondim show   <run_id>     스펙만 보고 리포트 재현 + 해시 검증
    $PY -m bdbot.cli nondim list                specs/ 목록
    $PY -m bdbot.cli run <case> [-- ...]        실행 (케이스 스크립트 위임)

종료 코드: 0 정상 · 1 FAIL(스키마 오류) · 2 BLOCKED(결측 미해소) · 3 사용법 오류
→ 스크립트·훅이 판정을 코드로 읽을 수 있습니다.
"""
from __future__ import annotations

import argparse
import subprocess
import sys
from pathlib import Path

from . import intake as _intake
from . import interactions as _inter
from . import physical as _physical

ROOT = Path(__file__).resolve().parent.parent
INTAKE_DIR = ROOT / "intake"

EXIT_OK, EXIT_FAIL, EXIT_BLOCKED, EXIT_USAGE = 0, 1, 2, 3

# 케이스 → 스크립트. `L3_ONLY` 는 무차원화(--report/--spec)까지만 있고 L4(실행)가 없습니다 —
# `status` 가 "O" 대신 "L3" 로 표시해 관통 스크립트와 구분합니다. 있는 것보다 **없는 것**을
# 정확히 보여주는 것이 이 표의 목적입니다.
CASE_SCRIPTS = {
    "trap-2d-5um": "cases/trap_2d_5um.py",
    "soft-r3-2d-A-sweep": "cases/soft_r3_2d.py",
    "abp-rod-2d-run-flip": "cases/abp_rod_2d.py",
    "trap-drag-2d-hex300": "cases/trap_drag_2d.py",
    "chain-bend-2d-oscill": "cases/chain_bend_2d.py",
    "chain-bend-2d-dlvo": "cases/chain_bend_dlvo_2d.py",
    "chain-relax-2d-dlvo": "cases/chain_relax_2d_dlvo.py",
    # 2026-08-28 (merge): `cases/network_3d.py` 는 있고 런도 냈는데 이 레지스트리에만
    # 빠져 있어서 `bdbot.cli run network` 가 "관통 스크립트가 없습니다" 로 거부했다.
    # status 표의 스크립트 칸이 `—` 였던 것도 스크립트가 없다는 뜻이 아니었다.
    "network": "cases/network_3d.py",
}
L3_ONLY = frozenset({"trap-drag-2d-hex300", "chain-bend-2d-oscill"})


def _resolve(folder: str) -> Path:
    p = Path(folder)
    if not p.exists():
        alt = INTAKE_DIR / folder
        if alt.exists():
            return alt
        cand = [d.name for d in sorted(INTAKE_DIR.iterdir()) if d.is_dir()]
        print(f"경로를 찾을 수 없습니다: {folder}\n케이스: {', '.join(cand)}", file=sys.stderr)
        raise SystemExit(EXIT_USAGE)
    return p


def _cases() -> list[Path]:
    if not INTAKE_DIR.exists():
        return []
    return [d for d in sorted(INTAKE_DIR.iterdir())
            if d.is_dir() and (d / "observation.yaml").exists()]


# ══════════════════════════════════════════════════════════════════════
def cmd_intake_init(args) -> int:
    p = Path(args.folder) if "/" in args.folder else INTAKE_DIR / args.folder
    ok, msg = _intake.init_template(p, force=args.force)
    print(msg)
    if ok:
        print("\n다음: 이미지를 읽고 전사부터 채우세요 (skill `bd-intake` §0 ①).")
        print(f"      확인:  $PY -m bdbot.cli intake check {p}")
    return EXIT_OK if ok else EXIT_USAGE


def cmd_intake_check(args) -> int:
    obs = _intake.load(_resolve(args.folder))
    print(_intake.render_check(obs))
    if obs.errors:
        return EXIT_FAIL
    unspec = [m for m in obs.open_missing
              if _inter.looks_like_interaction(str(m.get("symbol", "")))]
    if unspec:
        syms = ", ".join(str(m.get("symbol")) for m in unspec)
        print(f"\n상호작용이 미지정입니다 ({syms}).")
        print(f"  → 표준 후보 추천을 보려면:  "
              f"$PY -m bdbot.cli intake suggest {args.folder}")
    return EXIT_BLOCKED if obs.open_missing else EXIT_OK


def cmd_intake_suggest(args) -> int:
    """스케치에 상호작용이 없을 때 추천하고 **묻는다**. 결정은 사람이 한다."""
    obs = _intake.load(_resolve(args.folder))
    if obs.errors:
        print(_intake.render_check(obs))
        return EXIT_FAIL
    print(_inter.render_suggestion(obs))
    return EXIT_BLOCKED if obs.open_missing else EXIT_OK


def cmd_interactions_list(args) -> int:
    print(_inter.render_catalog())
    return EXIT_OK


def cmd_system_check(args) -> int:
    s = _physical.load(_resolve(args.folder))
    print(_physical.render_check(s))
    return EXIT_FAIL if s.errors else EXIT_OK


def cmd_status(args) -> int:
    """L0 → L2 → L3/런 진행도를 한 표로. 무엇이 어디서 막혀 있는지 보는 화면."""
    rows = []
    for d in _cases():
        obs = _intake.load(d)
        l0 = "FAIL" if obs.errors else ("BLOCKED" if obs.open_missing else "READY")
        has_sys = (d / "system.yaml").exists()
        if has_sys:
            s = _physical.load(d)
            l2 = "FAIL" if s.errors else "READY"
        else:
            l2 = "—"
        script = CASE_SCRIPTS.get(d.name)
        mark = ("L3" if d.name in L3_ONLY else "O") if script else "—"
        runs = sorted(ROOT.glob(f"runs/{d.name}__*"))
        runs = [r for r in runs if (r / "result.txt").exists()]
        # L3 산출물 — `nondim spec` 이 쓴 것. 런보다 앞선 단계라 따로 셉니다.
        specs = sorted((ROOT / "specs").glob(f"{d.name}__*.json")) \
            if (ROOT / "specs").exists() else []
        blockers = ", ".join(m.get("symbol", "?") for m in obs.open_missing) or "—"
        rows.append((d.name, l0, l2, mark, len(specs), len(runs), blockers))

    w = max([len(r[0]) for r in rows] + [8])
    print("=" * (w + 62))
    print("bdbot status — 스케치에서 런까지")
    print("=" * (w + 62))
    print(f"{'케이스':<{w}}  {'L0':<8}{'L2':<8}{'스크립트':<9}{'스펙':>5}{'런':>5}   차단 결측")
    print("-" * (w + 62))
    for name, l0, l2, sc, ns, nr, blk in rows:
        print(f"{name:<{w}}  {l0:<8}{l2:<8}{sc:<9}{ns:>5}{nr:>5}   {blk[:30]}")
    print("-" * (w + 62))
    n_ready = sum(1 for r in rows if r[1] == "READY")
    n_spec = sum(1 for r in rows if r[4] > 0)
    n_run = sum(1 for r in rows if r[5] > 0)
    print(f"케이스 {len(rows)}개 · L0 READY {n_ready} · L3 스펙 있는 케이스 {n_spec} · "
          f"런 있는 케이스 {n_run}")
    print("스크립트: O = 관통(L4 실행까지) · L3 = 무차원화까지만 (--report/--spec)")
    if any(r[1] == "FAIL" for r in rows):
        print("\nFAIL 케이스: `intake check <case>` 로 스키마 오류를 보세요.")
    blocked = [r[0] for r in rows if r[1] == "BLOCKED"]
    if blocked:
        print(f"\nBLOCKED {len(blocked)}개 — 스케치에 없는 값이라 지어내지 않고 멈춘 상태입니다.")
        print("  사람이 값을 주거나 KB에서 찾아야 합니다. `intake check <case>` 에 목록이 있습니다.")
    print("=" * (w + 58))
    return EXIT_OK


def _dispatch(case: str, extra: list[str]) -> int:
    script = CASE_SCRIPTS.get(case)
    if not script:
        print(f"'{case}' 에는 아직 관통 스크립트가 없습니다.\n"
              f"있는 것: {', '.join(CASE_SCRIPTS)}\n"
              f"(나머지는 물리계가 확정되지 않았습니다 — `bdbot.cli status` 참조)",
              file=sys.stderr)
        return EXIT_USAGE
    cmd = [sys.executable, str(ROOT / script), *extra]
    print(f"→ {' '.join(cmd[1:])}\n", flush=True)
    return subprocess.call(cmd, cwd=ROOT)


def cmd_nondim_report(args) -> int:
    """L3 리포트. **케이스 스크립트에 위임합니다** — 스케일 원장은 계마다 다릅니다.

    공통화된 부분(원장 구조·검사 분류·렌더러)은 `bdbot.scales/checks/report` 에 있고,
    어떤 스케일이 들어가는지는 케이스가 정합니다 (skill `bd-physics` §6.3).
    """
    return _dispatch(args.case, ["--report", *args.extra])


def cmd_nondim_spec(args) -> int:
    """L3 스펙을 `specs/<run_id>.json` 으로 씁니다 (실행 안 함) — 사람 확인 #3의 대상."""
    return _dispatch(args.case, ["--spec", *args.extra])


def cmd_nondim_show(args) -> int:
    """저장된 스펙**만** 보고 리포트를 다시 그립니다 — 케이스 스크립트 없이.

    L4가 쓸 경로를 사람이 먼저 확인하는 자리입니다. 여기서 다 보이면 스펙이 자족적이고,
    안 보이면 스펙에 무언가 빠진 것입니다. 손으로 고친 스펙은 해시 검증이 잡습니다.
    """
    from . import nondim as _nd
    p = Path(args.spec)
    if not p.exists():
        alt = ROOT / "specs" / args.spec
        alt = alt if alt.exists() else ROOT / "specs" / f"{args.spec}.json"
        if not alt.exists():
            have = sorted(x.name for x in (ROOT / "specs").glob("*.json")) \
                if (ROOT / "specs").exists() else []
            print(f"스펙을 찾을 수 없습니다: {args.spec}\n"
                  + ("있는 것:\n  " + "\n  ".join(have) if have else
                     "specs/ 가 비어 있습니다 — `nondim spec <case>` 로 먼저 만드세요."),
                  file=sys.stderr)
            return EXIT_USAGE
        p = alt
    try:
        spec = _nd.load(p)
    except Exception as e:
        print(f"스펙을 읽을 수 없습니다: {e}", file=sys.stderr)
        return EXIT_FAIL
    print(spec.render())
    ok, _ = spec.verify_hash()
    n_err = sum(1 for i in spec.raw.get("l3_issues", []) if i["level"] == "error")
    if not ok:
        print("\n손으로 고친 스펙입니다 — run_id 가 내용과 맞지 않습니다 (§16 규칙 2).")
        return EXIT_FAIL
    return EXIT_FAIL if (n_err or spec.verdict.startswith("FAIL")) else EXIT_OK


def cmd_nondim_list(args) -> int:
    from . import nondim as _nd
    d = ROOT / "specs"
    files = sorted(d.glob("*.json")) if d.exists() else []
    if not files:
        print("specs/ 가 비어 있습니다 — `nondim spec <case>` 로 만드세요.")
        return EXIT_OK
    print(f"{'run_id':<44}{'판정':<22}해시")
    print("-" * 78)
    for f in files:
        try:
            s = _nd.load(f)
            ok, _ = s.verify_hash()
            print(f"{s.run_id:<44}{s.verdict:<22}{'✓' if ok else '✗ 손으로 고침'}")
        except Exception as e:
            print(f"{f.name:<44}{'읽기 실패':<22}{e}")
    return EXIT_OK


def cmd_health(args) -> int:
    """L4 — 수치 건전성 판정 (마스터플랜 §0.2-B, `bdbot.health`).

    물리 검증기가 아닙니다. 발산·NaN/Inf·정지·붕괴, 그리고 L3 원장과의 대조뿐입니다.
    """
    import subprocess
    cmd = [sys.executable, str(ROOT / "tools" / "health.py")]
    if args.all:
        cmd.append("--all")
    elif args.gate:
        cmd += ["--gate", args.target]
    elif args.target:
        cmd.append(args.target)
    else:
        cmd.append("--all")
    return subprocess.call(cmd)


def cmd_run(args) -> int:
    return _dispatch(args.case, list(args.extra))


# ══════════════════════════════════════════════════════════════════════
def build_parser() -> argparse.ArgumentParser:
    ap = argparse.ArgumentParser(
        prog="bdbot", description="Brownian dynamics 파이프라인 — 앞단 진입점",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="종료 코드: 0 정상 · 1 FAIL · 2 BLOCKED · 3 사용법 오류")
    sub = ap.add_subparsers(dest="cmd", required=True)

    sub.add_parser("status", help="모든 케이스의 진행도").set_defaults(fn=cmd_status)

    p_in = sub.add_parser("intake", help="L0 인테이크")
    s_in = p_in.add_subparsers(dest="sub", required=True)
    q = s_in.add_parser("init", help="observation.yaml 템플릿 생성")
    q.add_argument("folder")
    q.add_argument("--force", action="store_true")
    q.set_defaults(fn=cmd_intake_init)
    q = s_in.add_parser("check", help="스키마 + 준비도 판정")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_intake_check)
    q = s_in.add_parser("suggest", help="미지정 상호작용에 대한 추천 + 질문")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_intake_suggest)
    q = s_in.add_parser("list", help="status 의 별칭")
    q.set_defaults(fn=cmd_status)

    q = sub.add_parser("interactions", help="콜로이드 상호작용 카탈로그")
    q.set_defaults(fn=cmd_interactions_list)

    p_sy = sub.add_parser("system", help="L2 물리계")
    s_sy = p_sy.add_subparsers(dest="sub", required=True)
    q = s_sy.add_parser("check", help="스키마 + tier + 유도값 재계산")
    q.add_argument("folder")
    q.set_defaults(fn=cmd_system_check)

    p_nd = sub.add_parser("nondim", help="L3 무차원화")
    s_nd = p_nd.add_subparsers(dest="sub", required=True)
    q = s_nd.add_parser("report", help="무차원화 리포트 (실행 안 함)")
    q.add_argument("case")
    q.add_argument("extra", nargs=argparse.REMAINDER)
    q.set_defaults(fn=cmd_nondim_report)
    q = s_nd.add_parser("spec", help="L3 스펙을 specs/<run_id>.json 으로 저장")
    q.add_argument("case")
    q.add_argument("extra", nargs=argparse.REMAINDER)
    q.set_defaults(fn=cmd_nondim_spec)
    q = s_nd.add_parser("show", help="저장된 스펙만 보고 리포트 재현 (케이스 코드 없이)")
    q.add_argument("spec", help="specs/ 아래 파일명 또는 run_id")
    q.set_defaults(fn=cmd_nondim_show)
    q = s_nd.add_parser("list", help="specs/ 목록 + 해시 검증")
    q.set_defaults(fn=cmd_nondim_list)

    q = sub.add_parser("health", help="L4 수치 건전성 판정 (물리 검증 아님)")
    q.add_argument("target", nargs="?", help="runs/<run_id> 또는 --gate 용 스펙 경로")
    q.add_argument("--all", action="store_true", help="완료된 런 전부")
    q.add_argument("--gate", action="store_true", help="실행 전 게이트만 (스펙 경로)")
    q.set_defaults(fn=cmd_health)

    q = sub.add_parser("run", help="케이스 실행 (스크립트 위임)")
    q.add_argument("case")
    q.add_argument("extra", nargs=argparse.REMAINDER)
    q.set_defaults(fn=cmd_run)
    return ap


def main(argv=None) -> int:
    args = build_parser().parse_args(argv)
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main())
