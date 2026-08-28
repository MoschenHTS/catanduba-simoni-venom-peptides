#!/usr/bin/env bash

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -uo pipefail
cd "$CSIMONI_BASE/md"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=$CSIMONI_BASE/.mamba
ENV=$CSIMONI_BASE/boltz_env
declare -A NTOX=( [cstx2_kv42]=32 [cstx3_kv42]=33 [cstx2_kv43]=32 [cstx3_kv43]=33 )
for sys in cstx2_kv42 cstx3_kv42 cstx2_kv43 cstx3_kv43; do
  for rep in prod_rep1 prod_rep2 prod_rep3 prod_rep4 prod_rep5 prod_rep6 prod_rep7 prod_rep8; do
    [ -f "$sys/${rep}_membrane.dat" ] && { echo "[$sys] $rep membrane done, skip"; continue; }
    echo "[$sys] $rep membrane analysis  $(date +%T)"
    $MM run -p "$ENV" python 04_analyze_membrane.py "$sys" "$rep" "${NTOX[$sys]}" 2>&1 | grep -vE "Warning|deprecat"
  done
done
echo ">>> MEMBRANE ANALYSIS ALL DONE  $(date)"
