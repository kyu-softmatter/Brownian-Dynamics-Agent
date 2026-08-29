#!/bin/zsh
# chain-bend-2d-dlvo — 미완 런 재개 (2026-08-06 중단분)
#
# 중단 시점: kt100 10/12 · position 6/12 완료. 아래가 남은 6런이다.
# 전부 **JKR 분기**라 느리다 — force.Custom 이 매 스텝 파이썬을 호출해서
# ~5,900 steps/s (DLVO 의 1/17). 런당 13.9M 스텝 = **약 40분**.
#
#   $ zsh scratch/resume_dlvo_runs.sh          # 4병렬, 총 ~80분
#   $ zsh scratch/resume_dlvo_runs.sh 2        # 2병렬 (CPU 여유 두기)
#
# ★ 이미 끝난 런은 run_id 가 콘텐츠 주소라 자동으로 건너뛴다 — 중복 실행 걱정 없다.
#   단, `--force` 를 붙이면 덮어쓰므로 쓰지 말 것.
#
# 끝난 뒤 할 것:
#   $PY scratch/backfill_bow.py                      # 굽음 소급 채움
#   $PY scratch/analyze_correlations.py              # 위치 상관
#   → 3조건(trap kt1 / trap kt100 / position) × 2분기(DLVO/JKR) 비교 그래프

set -e
cd "$(dirname "$0")/.."
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
NPAR=${1:-4}

JOBS=$(mktemp)
cat > "$JOBS" << 'EOF'
kt100 5
kt100 6
position 1
position 2
position 3
position 4
position 5
position 6
EOF

RUNNER=$(mktemp)
cat > "$RUNNER" << 'EOF'
#!/bin/zsh
cd "$RESUME_ROOT"
PY=/opt/homebrew/Caskroom/miniconda/base/envs/simulation_bot/bin/python
mode=$1; s=$2
# ★ zsh 은 bash 와 달리 `$extra` 를 단어분할하지 않는다 — 배열로 두거나 `${=extra}` 를
#   써야 한다. 처음에 문자열로 뒀다가 `--drive-mode position` 이 통짜 한 인자로 넘어가
#   8런 전부 argparse 에서 죽었다 (unrecognized arguments).
case $mode in
  kt100)    extra=(--kt-scale 100) ;;
  position) extra=(--drive-mode position) ;;
  *)        echo "FAILED $mode seed=$s (알 수 없는 모드)"; exit 1 ;;
esac
log=/tmp/resume_${mode}_s${s}.log
$PY cases/chain_bend_dlvo_2d.py --n 9 --omega 3000 --amp 1470 --cycles 20 \
    "${extra[@]}" --jkr --seed $s --run > "$log" 2>&1
rc=$?
# ★ 종료코드를 **반드시** 본다. 예전엔 무조건 "done" 을 찍어서 크래시가 성공으로 보였다
#   (같은 함정으로 배치 하나가 48/48 done 인데 metrics 는 37개였다 — KB tooling 참조).
if [ $rc -ne 0 ]; then
  echo "FAILED $mode seed=$s (rc=$rc) — $log"
else
  echo "done $mode seed=$s"
fi
EOF
chmod +x "$RUNNER"
export RESUME_ROOT="$PWD"

echo "재개: 남은 JKR 런 8개 (kt100 2 + position 6), ${NPAR}병렬"
echo "로그: /tmp/resume_<mode>_s<seed>.log"
xargs -P "$NPAR" -L 1 "$RUNNER" < "$JOBS"
echo "RESUME_DONE"
rm -f "$JOBS" "$RUNNER"

# ★ '완료 개수'가 아니라 **metrics.json 개수**로 검증한다 (KB tooling: xargs 는 자식이
#   크래시해도 done 을 찍는다). 기대치와 다르면 여기서 드러난다.
echo ""
echo "검증 — 실제 산출물 개수:"
for pat in "jkr-kt100" "jkr-position"; do
  cnt=$(ls -d runs/chain-bend-2d-dlvo__n9-w3000-a1470-${pat}__*/metrics.json 2>/dev/null | wc -l | tr -d ' ')
  echo "  ${pat}: ${cnt} 런"
done
echo "(기대: jkr-kt100 6 · jkr-position 6)"
