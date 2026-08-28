#!/usr/bin/env bash
# Chain the 3 remaining complexes through the full pipeline, SEQUENTIALLY
# (single GPU). Each 01_pipeline.sh call is itself idempotent (skips finished
# stages), so this whole script is safe to re-launch after an interruption.
#
# Order: CsTx3-Kv4.3 (cleanest interface in the Boltz-2 co-folding stage) -> CsTx2-Kv4.3 -> CsTx3-Kv4.2

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -uo pipefail
cd "$CSIMONI_BASE/md"

bash 01_pipeline.sh cstx3_kv43 33 Kv43
bash 01_pipeline.sh cstx2_kv43 32 Kv43
bash 01_pipeline.sh cstx3_kv42 33 Kv42

echo "=========================================================="
echo ">>> ALL REMAINING COMPLEXES DONE  $(date)"
echo "=========================================================="
