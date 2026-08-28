#!/usr/bin/env bash
# Idempotent resume — run this ANY time (after a shutdown, next session, etc.).
# It (re)launches 11_run_extra_replicas.sh detached (survives closing this
# terminal/VSCode) and continues from wherever it left off (every stage checks
# for its own output file / .cpt before redoing work).

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

cd "$CSIMONI_BASE/md" || exit 1
if pgrep -f "11_run_extra_replicas.sh" >/dev/null || pgrep -f "gmx mdrun" >/dev/null; then
  echo "already running"; exit 0
fi
echo "launching (detached)"
setsid nohup bash 11_run_extra_replicas.sh > run_extra_replicas.log 2>&1 < /dev/null &
echo "monitor: tail -f md/run_extra_replicas.log   |   stop: pkill -f 11_run_extra_replicas.sh"
