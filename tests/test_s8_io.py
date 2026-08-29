"""S8 뼈대 — run 디렉터리·해시·봉인.

이 파일이 지키는 것: **봉인이 실제로 작동하는가.** 봉인이 조용히 통과하면
사후합리화를 막는 장치가 하나도 없는 상태로 파이프라인이 돌아간다.
"""
from __future__ import annotations

import datetime as dt
import json
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


def test_verify_seal_on_real_first_run():
    """실제 첫 완주 런의 봉인이 지금도 유효한가 (회귀)."""
    p = io.REPO_ROOT / "runs" / "2026-07-28_trap-2d-5um_2dfb9d"
    if not p.exists():
        pytest.skip("runs/ 는 gitignore 대상 — 이 체크아웃에 없다")
    v = io.verify_seal(io.RunDir(p))
    assert v.ok, v.summary()


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
