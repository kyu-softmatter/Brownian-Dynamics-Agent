"""S8 뼈대 — run 디렉터리·해시·봉인.

이 파일이 지키는 것: **봉인이 실제로 작동하는가.** 봉인이 조용히 통과하면
사후합리화를 막는 장치가 하나도 없는 상태로 파이프라인이 돌아간다.
"""
from __future__ import annotations

import datetime as dt
import json
import sys
from pathlib import Path

import pytest

from simbot import io


# =============================================================================
# 해시
# =============================================================================
def test_sha256_text_matches_known_value():
    # 표준 sha256 — 외부 도구와 같은 값이어야 한다
    assert io.sha256_text("") == (
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855")


def test_sha256_file_matches_text(tmp_path):
    p = tmp_path / "a.txt"
    p.write_text("hello 한글\n", encoding="utf-8")
    assert io.sha256_file(p) == io.sha256_text("hello 한글\n")


def test_payload_hash_is_key_order_independent():
    """키 순서만 다른 dict 은 같은 계다. 다르면 캐시가 매번 무효화된다."""
    assert io.sha256_payload({"a": 1, "b": 2}) == io.sha256_payload({"b": 2, "a": 1})


def test_payload_hash_distinguishes_values():
    assert io.sha256_payload({"a": 1}) != io.sha256_payload({"a": 1.0000001})


def test_code_hash_covers_analysis_subpackage(tmp_path):
    """`analysis/` 를 빼면 분석 코드가 바뀌어도 해시가 그대로다 — 봉인에 구멍."""
    pkg = tmp_path / "pkg"
    (pkg / "analysis").mkdir(parents=True)
    (pkg / "a.py").write_text("x = 1")
    (pkg / "analysis" / "b.py").write_text("y = 1")
    before = io.code_hash(pkg)
    (pkg / "analysis" / "b.py").write_text("y = 2")
    assert io.code_hash(pkg) != before


def test_code_hash_ignores_pycache(tmp_path):
    pkg = tmp_path / "pkg"
    (pkg / "__pycache__").mkdir(parents=True)
    (pkg / "a.py").write_text("x = 1")
    before = io.code_hash(pkg)
    (pkg / "__pycache__" / "junk.py").write_text("noise")
    assert io.code_hash(pkg) == before


def test_code_hash_detects_file_rename(tmp_path):
    """내용 합만 해시하면 파일 이름 변경을 놓친다 — 경로도 해시에 넣었는지 확인."""
    pkg = tmp_path / "pkg"
    pkg.mkdir()
    (pkg / "a.py").write_text("x = 1")
    before = io.code_hash(pkg)
    (pkg / "a.py").rename(pkg / "b.py")
    assert io.code_hash(pkg) != before


def test_env_versions_reports_absent_without_raising():
    v = io.env_versions()
    assert "python" in v
    for name in io.ENV_PACKAGES:
        assert isinstance(v[name], str) and v[name]


# =============================================================================
# provenance — **한 곳에서만 만든다**
# =============================================================================
#  ★ 2026-07-29: `report.reproducibility_section` 은 `env_hash` 를 읽는데 두 러너가
#    자기 manifest 를 손으로 만들면서 그 키를 빼먹었다 → 리포트의 재현 정보가
#    조용히 빈칸으로 렌더됐다. 아래 테스트들이 그 어긋남을 고정한다.
def test_provenance_supplies_every_key_the_report_renders():
    """`provenance()` 가 리포트가 읽는 키를 전부 준다 — 어휘가 갈라지면 빈칸이 된다."""
    p = io.provenance()
    for key in ("code_hash", "git_rev", "git_dirty", "env_hash", "env"):
        assert key in p, f"{key} 가 provenance 에 없다"
    for name in ("hoomd", "numpy", "scipy", "freud"):
        assert name in p["env"], f"{name} 이 env 에 없다"


def test_build_manifest_gets_provenance_from_the_single_definition():
    """`build_manifest` 가 `provenance()` 를 쓴다 — 손으로 나열하면 갈라진다."""
    man = io.build_manifest(run_id="r", spec_hash="abc", seed=1)
    prov = io.provenance()
    for key in ("code_hash", "env_hash"):
        assert man[key] == prov[key]
    assert man["env"] == prov["env"]


def test_provenance_records_the_driver_as_a_repo_relative_path():
    """드라이버는 **repo 상대경로**로 적는다 — 절대경로면 홈 디렉터리가 박힌다."""
    p = io.provenance(driver=io.REPO_ROOT / "simbot" / "io.py")
    assert p["driver"] == "simbot/io.py"
    assert p["driver_hash"] == io.file_hash(io.REPO_ROOT / "simbot" / "io.py")
    assert len(p["driver_hash"]) == 12


def test_provenance_hashes_every_driver_when_several_are_given():
    """★ 드라이버가 여러 파일이면 전부 해싱한다 — 하나만 잡으면 주장이 거짓이 된다."""
    a = io.REPO_ROOT / "simbot" / "io.py"
    b = io.REPO_ROOT / "simbot" / "run.py"
    p = io.provenance(driver=[a, b])

    assert p["driver"] == ["simbot/io.py", "simbot/run.py"]
    assert set(p["drivers"]) == {"simbot/io.py", "simbot/run.py"}
    assert p["drivers"]["simbot/io.py"] == io.file_hash(a)
    #  합성 해시는 파일 하나만 바뀌어도 달라져야 한다
    other = dict(p["drivers"]); other["simbot/io.py"] = "0" * 12
    assert io.sha256_payload(other)[:12] != p["driver_hash"]
    #  단일 경로는 예전 형식을 유지한다 (리더 호환)
    single = io.provenance(driver=a)
    assert single["driver"] == "simbot/io.py"
    assert single["driver_hash"] == io.file_hash(a)


def test_provenance_marks_a_missing_driver_instead_of_raising():
    """없는 드라이버는 `"?"` 로 표기한다 — 예외로 런을 죽이지 않되 숨기지도 않는다."""
    p = io.provenance(driver=io.REPO_ROOT / "scripts" / "does_not_exist.py")
    assert p["driver_hash"] == "?"


def test_driver_hash_changes_with_content_but_code_hash_does_not(tmp_path):
    """★ `code_hash` 는 `simbot/` 만 덮는다 — 드라이버 변경을 못 잡는다.

    이것이 `driver_hash` 가 따로 필요한 이유다. 런의 `A` 목록·시드·분석 창을
    정하는 것은 `scripts/` 의 드라이버이고, 그것이 해시에 없으면 산출물만으로
    "무엇이 이 런을 만들었는가" 에 답할 수 없다.
    """
    drv = tmp_path / "driver.py"
    drv.write_text("AMPLITUDES = (0.1, 1.0, 10.0)\n")
    before_driver = io.file_hash(drv)
    before_code = io.code_hash()

    drv.write_text("AMPLITUDES = (0.1, 1.0, 100.0)\n")     # 물리가 바뀌는 변경
    assert io.file_hash(drv) != before_driver
    assert io.code_hash() == before_code, "simbot 은 안 바뀌었는데 code_hash 가 바뀌었다"


# =============================================================================
# run_id
# =============================================================================
def test_run_id_is_deterministic_given_date():
    rid = io.new_run_id("Trap 2D 5um!", "c59e93fd24a2", dt.date(2026, 7, 28))
    assert rid == "2026-07-28_trap-2d-5um_c59e93"


def test_run_id_uses_six_hash_chars():
    rid = io.new_run_id("x", "abcdef0123456789", dt.date(2026, 1, 2))
    assert rid.endswith("_abcdef")


def test_slugify_never_returns_empty():
    assert io.slugify("!!!") == "run"


# =============================================================================
# RunDir
# =============================================================================
def test_rundir_rejects_unknown_stage(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    with pytest.raises(KeyError):
        rd.file("no_such_stage")


def test_rundir_roundtrips_json(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write_json("metrics", {"a": 1.5})
    assert rd.read_json("metrics") == {"a": 1.5}


def test_completed_stages_lists_only_existing(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    assert rd.completed_stages() == []
    rd.write("prediction", "# p")
    assert rd.completed_stages() == ["prediction"]


# =============================================================================
# 봉인 — 이 절이 이 파일의 존재 이유다
# =============================================================================
@pytest.fixture
def sealed_run(tmp_path):
    rd = io.RunDir.create(tmp_path, "2026-07-28_x_abc123")
    rd.write("intake", "# S1\n")
    rd.write("prediction", "# S2\nD = 1.00 ± 0.03\n")
    rd.write("spec", "card: x\n")
    io.write_seal(rd)
    return rd


def test_seal_passes_when_untouched(sealed_run):
    v = io.verify_seal(sealed_run)
    assert v.ok and len(v.entries) == 3


def test_seal_catches_edited_prediction(sealed_run):
    """★ 이 테스트가 통과하지 않으면 사후합리화를 막는 장치가 없다."""
    sealed_run.write("prediction", "# S2\nD = 0.42 ± 0.25   (결과 보고 나서 고침)\n")
    v = io.verify_seal(sealed_run)
    assert not v.ok
    assert io.RUN_LAYOUT["prediction"] in v.changed


def test_seal_catches_whitespace_only_edit(sealed_run):
    """공백 한 칸도 잡아야 한다 — '사소한 수정'에 예외를 두면 봉인이 아니다."""
    p = sealed_run.file("prediction")
    p.write_text(p.read_text() + " ", encoding="utf-8")
    assert not io.verify_seal(sealed_run).ok


def test_seal_catches_deleted_document(sealed_run):
    sealed_run.file("spec").unlink()
    v = io.verify_seal(sealed_run)
    assert not v.ok and io.RUN_LAYOUT["spec"] in v.missing


def test_seal_catches_missing_seal_file(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("prediction", "# p")
    assert not io.verify_seal(rd).ok


def test_seal_reports_document_added_after_sealing(tmp_path):
    """봉인 후에 만든 예측은 봉인되지 않았다 — 통과로 보고하면 안 된다."""
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("intake", "# S1\n")
    io.write_seal(rd)                        # 이 시점에 prediction 이 없다
    rd.write("prediction", "# S2 (실행 후 작성)\n")
    v = io.verify_seal(rd)
    assert not v.ok
    assert io.RUN_LAYOUT["prediction"] in v.unsealed


def test_write_seal_refuses_when_nothing_to_seal(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    with pytest.raises(FileNotFoundError):
        io.write_seal(rd)


def test_seal_file_is_shasum_compatible(sealed_run):
    """`shasum -a 256 -c` 로 우리 코드 없이 검증되어야 한다."""
    lines = sealed_run.read("seal").strip().splitlines()
    for line in lines:
        digest, sep, rel = line.partition("  ")
        assert sep == "  " and len(digest) == 64
        assert int(digest, 16) >= 0            # 16진수인지
        assert rel and not rel.startswith(" ")


# ── the archive: this section never ran ───────────────────────────────────
#
# ★ The previous revision was a single test, and that test pointed at a path
#   **absent from every commit in the history**: `runs/2026-07-28_trap-2d-5um_2dfb9d`
#   has been under `runs_s1s8/` since the initial commit. So it always skipped,
#   and its skip reason ("runs/ is gitignored") was **false** too --
#   `git check-ignore runs` exits 1 and 1311 files under `runs/` are tracked.
#   Pointed at the real path, the assertion failed. A silent skip, on a false
#   premise, sitting at the centre of the seal machinery -- and the consequence
#   was that `verify_seal` went 18 days without hashing the archive once.
#   So this section **names no individual path**: it globs every tracked seal,
#   and finding zero of them is itself a failure.

def _archived_seal_dirs() -> list[io.RunDir]:
    """Every run directory holding a tracked `SEALED.sha256`."""
    return [io.RunDir(p.parent)
            for p in sorted(io.REPO_ROOT.glob("runs*/*/SEALED.sha256"))]


def test_there_is_at_least_one_archived_seal_to_check():
    """Guard on the gate. The tests below iterate a list, so an empty list makes
    all of them pass silently -- the same reason CI's bash check treats
    `checked -eq 0` as an error."""
    dirs = _archived_seal_dirs()
    assert len(dirs) >= 21, [d.run_id for d in dirs]


def test_every_archived_seal_verifies_and_none_passes_on_nothing(monkeypatch):
    """★ A pass has to say **how many documents it hashed**.

    The previous revision's `ok=True` hashed nothing on all 21 archived
    directories, and 12 of them printed "2 documents unchanged after the run"
    alongside it. `entries` is the number of lines in the seal file, not the
    number of documents hashed.
    """
    hashed = []
    real = io.sha256_file
    monkeypatch.setattr(io, "sha256_file",
                        lambda q: (hashed.append(str(q)), real(q))[1])
    total = 0
    for rd in _archived_seal_dirs():
        hashed.clear()
        v = io.verify_seal(rd)
        assert v.ok, f"{rd.run_id}: {v.summary()}"
        assert len(hashed) == len(v.verified) == len(v.entries), (
            f"{rd.run_id}: the seal has {len(v.entries)} lines but "
            f"{len(hashed)} files were hashed and {len(v.verified)} verified")
        assert len(hashed) > 0, f"{rd.run_id}: passed having hashed 0 documents"
        total += len(hashed)
    #  ⚠ This was the literal `42` until 2026-09-16, and it broke the moment the
    #     sedimentation campaign sealed 8 more directories -- a hand-maintained
    #     count inside the very suite whose counts gate exists because 48 of 102
    #     documented counts had drifted. It now reads the SAME counter the gate
    #     and the README read, so the three cannot disagree.
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "verify"))
    import verify_counts
    expected = verify_counts.COUNTERS["sealed_documents"]()
    assert total == expected, (total, expected)
    #  CI re-counts the same thing with `shasum -a 256 -c`, importing none of
    #  this repository's code -- that independence is the point of the seal
    #  format, and it is why this assertion is allowed to import the counter.


def test_a_relocated_run_is_reported_drifted_and_not_unsealed():
    """The 9 runs whose recorded paths stopped resolving at the `runs/` →
    `runs_s1s8/` rename.

    The content is intact, so they must pass, and the fact that the path moved
    must be reported. The previous revision called these 9 `unsealed` -- "never
    sealed" -- which makes "we never sealed it" and "we sealed it and then
    renamed the directory" the same output. Those are exactly the two states this
    mechanism exists to distinguish.
    """
    drifted = [rd for rd in _archived_seal_dirs() if io.verify_seal(rd).drifted]
    assert len(drifted) == 9, [rd.run_id for rd in drifted]
    for rd in drifted:
        v = io.verify_seal(rd)
        assert v.ok and not v.unsealed and not v.changed, v.summary()
        assert len(v.drifted) == len(v.entries)
        assert "기록된 경로" in v.summary()      # summary() is Korean by design


def test_a_tampered_archived_document_is_caught(tmp_path):
    """★★ The absence of this test is why nobody noticed for 18 days.

    Appending a line to a document of a renamed run left the previous revision's
    verdict **bit-identical**: the recorded path was not a seal key, so the
    content was never hashed at all, and the result was
    `ok=False unsealed=[...]` both before and after the tamper.
    """
    import shutil
    src = next(rd for rd in _archived_seal_dirs() if io.verify_seal(rd).drifted)
    dst = tmp_path / src.path.name
    shutil.copytree(src.path, dst)
    rd = io.RunDir(dst)
    assert io.verify_seal(rd).ok, "the copy must pass first, or the comparison is void"

    target = dst / Path(next(iter(io.read_seal(rd)))).name
    target.write_bytes(target.read_bytes()
                       + "\n# edited after seeing the result\n".encode("utf-8"))
    v = io.verify_seal(rd)
    assert not v.ok, v.summary()
    assert target.name in v.changed, v
    assert target.name not in v.verified


def test_the_two_seal_implementations_agree_on_every_archived_seal():
    """★ The original defect was that two checkers gave different verdicts on the
    same documents.

    `bdbot.runcard.verify_seal` iterates the seal's lines and carries a fallback
    for a renamed directory, so it was right on all 21. `simbot.io.verify_seal`
    iterated the stage list and was wrong on all 21 -- 12 false passes and 9
    false violations. The two are deliberately not unified, and the reason is in
    `bdbot.runcard`'s module docstring (`bdbot` cannot import `simbot`). If they
    may not be unified, they must at least fail when they diverge.
    """
    from bdbot import runcard as RC
    for rd in _archived_seal_dirs():
        mine = io.verify_seal(rd)
        theirs_ok, problems = RC.verify_seal(rd.path, root=io.REPO_ROOT)
        assert mine.ok == theirs_ok, (rd.run_id, mine.summary(), problems)
        hard = [q for q in problems if not q.startswith("[warn]")]
        assert bool(mine.changed or mine.missing or mine.unsealed) == bool(hard), (
            rd.run_id, mine, problems)


def test_a_seal_naming_documents_that_are_not_there_cannot_pass(tmp_path):
    """A seal that exists while none of its documents do is not a pass -- hashing
    zero documents and passing was the core of the original defect."""
    rd = io.RunDir.create(tmp_path, "r1")
    rd.file("seal").write_text(
        "0" * 64 + "  runs/gone/02_prediction.md\n", encoding="utf-8")
    v = io.verify_seal(rd)
    assert not v.ok
    assert v.missing == ["02_prediction.md"] and not v.verified


def test_an_empty_seal_file_cannot_pass(tmp_path):
    """`write_seal` refuses to write an empty seal, but a truncated file must not
    pass either.

    With a document present, `unsealed` is the accurate word -- a prediction
    sitting beside an empty seal really is "not sealed".
    """
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("prediction", "# p\n")
    rd.file("seal").write_text("", encoding="utf-8")
    v = io.verify_seal(rd)
    assert not v.ok and not v.verified
    assert v.unsealed == ["02_prediction.md"], v


def test_a_seal_that_covers_nothing_at_all_says_so(tmp_path):
    """★ Having nothing to say is where the original defect lived: an empty seal
    and no sealable document either. The previous revision left all three problem
    lists empty and returned `ok=True`."""
    rd = io.RunDir.create(tmp_path, "r1")
    rd.file("seal").write_text("\n# a seal with only a comment\n", encoding="utf-8")
    v = io.verify_seal(rd)
    assert not v.ok, v
    assert not (v.changed or v.missing or v.unsealed or v.verified), v
    assert "0개" in v.summary(), v.summary()


# =============================================================================
# manifest
# =============================================================================
def test_manifest_records_reproducibility_fields(tmp_path):
    rd = io.RunDir.create(tmp_path, "r1")
    rd.write("prediction", "# p")
    io.write_seal(rd)
    man = io.build_manifest(run_id="r1", spec_hash="abc", seed=[1, 2, 3, 4], rundir=rd)
    for k in ("run_id", "spec_hash", "seed", "code_hash", "git_rev", "env_hash",
              "env", "sealed"):
        assert k in man, k
    json.dumps(man)                            # 직렬화 가능해야 한다


def test_manifest_dirty_flag_is_tristate():
    """dirty 를 판정할 수 없을 때 False 로 보고하면 '재현 가능'을 거짓 주장한다."""
    assert io.git_dirty(Path("/")) in (True, False, None)
