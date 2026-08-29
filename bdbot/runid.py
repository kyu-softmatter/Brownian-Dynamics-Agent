"""콘텐츠 주소 지정 run_id + 런 디렉토리 준비 (마스터플랜 §14).

같은 스펙 → 같은 `run_id` → **재실행하지 않습니다.** 두 케이스가 똑같이 이 규약을 썼습니다.

`nhex`를 인자로 둔 이유: 1-A는 12자, 1-B는 10자를 쓰고 있었습니다. 기존 런 디렉토리와
`run_id`를 그대로 유지해야 하므로(재현성) 통일하지 않고 케이스가 선언합니다.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

PRESERVE = {"record.json"}     # 재실행해도 지우지 않는 파일

# ★ run_id 해시에서 빼는 키 — 문서·출처·유도값.
#   실제로 물렸습니다: `derived_from` 필드를 system.yaml 에 추가했더니 1-A의 run_id가
#   70b9394e7310 → dc67e4e2b825 로 바뀌었습니다. 스펙에 YAML 전체를 넣어 해시했기 때문입니다.
#   **주석 한 줄을 고쳐서 런이 무효화되면 콘텐츠 주소는 쓸모가 없습니다.**
#   run_id 는 "무엇을 시뮬레이션했는가"만 추적해야 합니다.
#   (1-B는 물리 필드를 직접 나열해 이 문제가 없었습니다 — 그쪽이 옳은 설계였습니다)
DOC_KEYS = frozenset({
    "description", "derived_from", "not_verified", "required_convergence_checks",
    "derived_scales", "dimensionless", "source", "note", "source_note",
    "interpretation", "deviates_from_sketch", "expr", "role", "meaning", "tier",
    "what", "proposed", "followup", "lean", "confirmed_by",
})


def physics_only(node):
    """문서·출처·유도 필드를 재귀적으로 제거. 남는 것은 물리를 정하는 값뿐입니다."""
    if isinstance(node, dict):
        return {k: physics_only(v) for k, v in node.items() if k not in DOC_KEYS}
    if isinstance(node, list):
        return [physics_only(v) for v in node]
    return node


def spec_hash(spec: dict, nhex: int = 12) -> str:
    """스펙의 정렬된 JSON을 sha256 → 앞 nhex자. 순서·공백에 무관하게 결정적."""
    blob = json.dumps(spec, sort_keys=True, default=str).encode()
    return hashlib.sha256(blob).hexdigest()[:nhex]


def content_run_id(label: str, spec: dict, tag: str | None = None, nhex: int = 12) -> str:
    mid = f"{tag}__" if tag else ""
    return f"{label}__{mid}{spec_hash(spec, nhex)}"


def prepare_outdir(outdir: Path, force: bool = False) -> tuple[bool, str]:
    """(실행할 것인가, 안내 메시지). 완료된 런이 있으면 False.

    `result.txt` 존재 = 완료로 봅니다 (부분 산출물이 남은 디렉토리는 지우고 다시 씁니다).
    """
    if (outdir / "result.txt").exists() and not force:
        prev = (outdir / "result.txt").read_text()
        tail = prev.split("결과 —")[-1] if "결과 —" in prev else prev[-1200:]
        return False, f"\n이미 완료된 런입니다: runs/{outdir.name}/  (--force 로 재실행)\n{tail}"
    if outdir.exists():
        for f in outdir.iterdir():
            # ★ record.json 은 남긴다. 교훈(KB 엔트리)은 런 산출물보다 오래 살아야 한다.
            #   run_id 가 콘텐츠 주소이므로 같은 디렉토리면 같은 스펙이고, 이전 교훈이 그대로 유효하다.
            #   (1-C에서 --force 재실행으로 교훈 6건을 날린 뒤 넣은 방어)
            if f.is_file() and f.name not in PRESERVE:
                f.unlink()
    outdir.mkdir(parents=True, exist_ok=True)
    return True, ""


def write_spec(outdir: Path, spec: dict) -> None:
    (outdir / "spec.json").write_text(json.dumps(spec, indent=2, default=str))


def list_artifacts(outdir: Path, root: Path) -> list[str]:
    lines = [f"\n산출물: {outdir.relative_to(root)}/"]
    for f in sorted(outdir.iterdir()):
        lines.append(f"   {f.name:<22} {f.stat().st_size / 1024:8.1f} KB")
    return lines


__all__ = ["spec_hash", "content_run_id", "prepare_outdir", "write_spec", "list_artifacts"]
