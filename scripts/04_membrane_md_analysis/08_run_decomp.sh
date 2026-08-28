#!/usr/bin/env bash
# Per-residue MM-GBSA decomposition, all 4 complexes x 8 replicas (extended
# 2026-08 from n=3, matching the RMSD/MM-GBSA totals extension). Reuses the
# tv_index_full.ndx built by 07_run_mmpbsa.sh. 100 frames/replica (interval=10).

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -uo pipefail
BASE=$CSIMONI_BASE/md
cd "$BASE"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=$CSIMONI_BASE/.mamba
GENV=$CSIMONI_BASE/gmx_env
MMENV=$CSIMONI_BASE/mmpbsa_env
export PATH="$GENV/bin:$PATH"

run_one () {
  local SYS="$1"
  cd "$BASE/$SYS/mmpbsa"
  for rep in 1 2 3 4 5 6 7 8; do
    local OUT="rep${rep}_FINAL_DECOMP_MMPBSA.dat"
    [ -f "$OUT" ] && { echo "[$SYS] rep$rep decomp done, skip"; continue; }
    echo "[$SYS] rep$rep decomposition  $(date +%T)"
    rm -rf COM_* RECEPTOR_* LIGAND_* _GMXMMPBSA_* gmx_MMPBSA.log leap.log 2>/dev/null
    $MM run -p "$MMENV" gmx_MMPBSA -O -i ../../09_mmpbsa_decomp.in -cs "../prod_rep${rep}.tpr" \
        -ci tv_index_full.ndx -cg 1 0 -ct "../prod_rep${rep}.xtc" -nogui \
        > "rep${rep}_decomp_run.log" 2>&1
    if [ -f FINAL_DECOMP_MMPBSA.dat ]; then
      mv FINAL_DECOMP_MMPBSA.dat "$OUT"
      echo "[$SYS] rep$rep decomp OK"
    else
      echo "[$SYS] rep$rep decomp FAILED - see mmpbsa/rep${rep}_decomp_run.log"
    fi
  done
}

run_one cstx2_kv42
run_one cstx3_kv42
run_one cstx2_kv43
run_one cstx3_kv43
echo ">>> DECOMP ALL DONE  $(date)"
