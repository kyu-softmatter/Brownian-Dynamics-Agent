#!/bin/zsh
# chain-bend-2d-dlvo — resume the unfinished runs (the 2026-08-06 interruption)
#
# Where it stopped: kt100 10/12 · position 6/12 done. The 6 below are what is left.
# All of them are the **JKR branch**, so they are slow — force.Custom calls into
# python every step, giving ~5,900 steps/s (1/17 of DLVO). 13.9M steps per run =
# **about 40 minutes**.
#
#   $ zsh scratch/resume_dlvo_runs.sh          # 4-way parallel, ~80 min total
#   $ zsh scratch/resume_dlvo_runs.sh 2        # 2-way (leave the CPU some room)
#
# ★ Runs that already finished are skipped automatically, because run_id is a
#   content address — no need to worry about running them twice. But `--force`
#   overwrites, so do not pass it.
#
# When it finishes:
#   $PY scratch/backfill_bow.py                      # backfill the bow
#   $PY scratch/analyze_correlations.py              # position correlations
#   → comparison plots over 3 conditions (trap kt1 / trap kt100 / position)
#     × 2 branches (DLVO/JKR)

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
# ★ Unlike bash, zsh does **not** word-split `$extra` — it has to be an array, or
#   written `${=extra}`. Holding it as a string the first time made
#   `--drive-mode position` arrive as one single argument, and all 8 runs died in
#   argparse (unrecognized arguments).
case $mode in
  kt100)    extra=(--kt-scale 100) ;;
  position) extra=(--drive-mode position) ;;
  *)        echo "FAILED $mode seed=$s (unknown mode)"; exit 1 ;;
esac
log=/tmp/resume_${mode}_s${s}.log
$PY cases/chain_bend_dlvo_2d.py --n 9 --omega 3000 --amp 1470 --cycles 20 \
    "${extra[@]}" --jkr --seed $s --run > "$log" 2>&1
rc=$?
# ★ Look at the exit code, **always**. An earlier version printed "done"
#   unconditionally, so a crash looked like a success (the same trap once left a
#   batch reading 48/48 done with only 37 metrics — see KB tooling).
if [ $rc -ne 0 ]; then
  echo "FAILED $mode seed=$s (rc=$rc) — $log"
else
  echo "done $mode seed=$s"
fi
EOF
chmod +x "$RUNNER"
export RESUME_ROOT="$PWD"

echo "resume: 8 JKR runs left (kt100 2 + position 6), ${NPAR}-way parallel"
echo "logs: /tmp/resume_<mode>_s<seed>.log"
xargs -P "$NPAR" -L 1 "$RUNNER" < "$JOBS"
echo "RESUME_DONE"
rm -f "$JOBS" "$RUNNER"

# ★ Verify by the **number of metrics.json**, not by the 'done' count (KB tooling:
#   xargs prints done even when the child crashed). A mismatch shows up here.
echo ""
echo "verification — the actual number of artefacts:"
for pat in "jkr-kt100" "jkr-position"; do
  cnt=$(ls -d runs/chain-bend-2d-dlvo__n9-w3000-a1470-${pat}__*/metrics.json 2>/dev/null | wc -l | tr -d ' ')
  echo "  ${pat}: ${cnt} runs"
done
echo "(expected: jkr-kt100 6 · jkr-position 6)"
