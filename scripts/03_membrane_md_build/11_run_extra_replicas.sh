#!/usr/bin/env bash
# Extend all 4 complexes from 3 -> 8 replicas (5 new 200ns runs each, 20 total)
# for statistical power (see 03_power_analysis.py). Equilibration already done for
# all 4 - this calls 09_prod.sh directly (idempotent: reps 1-3 skip immediately).

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -uo pipefail
cd "$CSIMONI_BASE/md"

bash 09_prod.sh cstx2_kv42
bash 09_prod.sh cstx3_kv42
bash 09_prod.sh cstx2_kv43
bash 09_prod.sh cstx3_kv43

echo "=========================================================="
echo ">>> ALL 8-REPLICA PRODUCTION DONE  $(date)"
echo "=========================================================="
