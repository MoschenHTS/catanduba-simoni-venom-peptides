#!/usr/bin/env bash
# Continuous production: 8 replicas x 200 ns (extended from 3 -> 8 for statistical
# power, see 03_power_analysis.py), checkpointed (survives shutdown), nice'd (PC
# usable), idempotent/resumable. Reps 1-3 already done are skipped automatically.
# Usage: 09_prod.sh <system_dir>

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -uo pipefail
SYS="$1"
cd "$CSIMONI_BASE/md/$SYS"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=$CSIMONI_BASE/.mamba
GENV=$CSIMONI_BASE/gmx_env
G="$MM run -p $GENV gmx"
MD="nice -n 10 $MM run -p $GENV gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -cpt 5"
NSTEPS=100000000     # 200 ns @ 2 fs
for rep in 1 2 3 4 5 6 7 8; do
  d="prod_rep${rep}"
  [ -f "$d.gro" ] && { echo ">> [$SYS] $d complete, skip"; continue; }
  [ -f "$d.tpr" ] || $G grompp -f ../mdp/prod.mdp -c npt2.gro -p system.top -n index.ndx -o "$d.tpr" -maxwarn 5
  echo ">> [$SYS] $d : 200 ns  $(date)"
  $MD -deffnm "$d" -nsteps "$NSTEPS" -cpi "$d.cpt"
done
echo ">> [$SYS] ALL PRODUCTION DONE  $(date)"
