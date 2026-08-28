#!/usr/bin/env bash
# Batch driver for 02_analyze_md.py (pose RMSD / paddle contact / R3 salt bridge),
# 4 systems x 8 replicas, idempotent (skips replicas whose _rmsd.dat exists).
#
# RECONSTRUCTED 2026-08. The original loop that produced the committed
# prod_repN_rmsd.dat files was run ad hoc and never saved as a script (only its
# output survived, in rmsd_extra_run.log). This file restores it, mirroring
# 03_run_membrane_batch.sh exactly.
#
# WHY THIS MATTERS: 02_analyze_md.py takes <replica> <n_tox_res> <vsd_type> POSITIONALLY
# and silently defaults to 32/Kv42. Running it on a CsTx3 system (33 residues) without
# arguments therefore mis-assigns the toxin/VSD chain boundary and produces wrong RMSD
# with no error. The per-system values below are the authoritative ones - CsTx2 = 32,
# CsTx3 = 33 - matching 03_run_membrane_batch.sh and the residue counts recorded in the
# _ss.dat headers. Do not run 02_analyze_md.py by hand without them.
#
# Unlike the other analyses, 02_analyze_md.py must be run from INSIDE the system dir.

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -uo pipefail
cd "$(dirname "$0")"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=${MAMBA_ROOT_PREFIX:-$CSIMONI_BASE/.mamba}
ENV=${CSIMONI_BOLTZ_ENV:-$CSIMONI_BASE/boltz_env}

declare -A NTOX=( [cstx2_kv42]=32 [cstx3_kv42]=33 [cstx2_kv43]=32 [cstx3_kv43]=33 )
declare -A VSD=(  [cstx2_kv42]=Kv42 [cstx3_kv42]=Kv42 [cstx2_kv43]=Kv43 [cstx3_kv43]=Kv43 )

for sys in cstx2_kv42 cstx3_kv42 cstx2_kv43 cstx3_kv43; do
  for rep in prod_rep1 prod_rep2 prod_rep3 prod_rep4 prod_rep5 prod_rep6 prod_rep7 prod_rep8; do
    [ -f "$sys/${rep}_rmsd.dat" ] && { echo "[$sys] $rep RMSD done, skip"; continue; }
    echo "[$sys] $rep RMSD analysis  $(date +%T)"
    ( cd "$sys" && $MM run -p "$ENV" python ../02_analyze_md.py "$rep" "${NTOX[$sys]}" "${VSD[$sys]}" ) \
      2>&1 | grep -vE "Warning|deprecat"
  done
done
echo ">>> RMSD ANALYSIS ALL DONE  $(date)"
