#!/usr/bin/env bash
# Build the solvated, ionized, minimized membrane system for one complex.
# Usage: 07_build.sh <system_dir> <n_prot_atoms>

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -e
SYS="$1"; NPROT="$2"
cd "$CSIMONI_BASE/md/$SYS"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=$CSIMONI_BASE/.mamba
GENV=$CSIMONI_BASE/gmx_env
PENV=$CSIMONI_BASE/boltz_env
G="$MM run -p $GENV gmx"
PY="$MM run -p $PENV python"

echo ">>> [$SYS][1] solvate  $(date +%T)"
$G solvate -cp embedded.gro -cs spc216.gro -p system.top -o solv.gro
$PY ../05_fix_top.py solv.gro system.top

echo ">>> [$SYS][2] clean core waters"
$PY ../04_clean_core_water.py solv.gro cleaned.gro system.top "$NPROT"
$PY ../05_fix_top.py cleaned.gro system.top

echo ">>> [$SYS][3] genion — 0.15 M KCl, neutralize"
$G grompp -f ../mdp/ions.mdp -c cleaned.gro -p system.top -o ions.tpr -maxwarn 5
printf "SOL\n" | $G genion -s ions.tpr -o ionized.gro -p system.top -pname POT -nname CLA -conc 0.15 -neutral

echo ">>> [$SYS][4] index groups (MEMBRANE / SOLV)"
$PY ../06_make_index.py ionized.gro index.ndx "$NPROT"

echo ">>> [$SYS][5] energy minimization  $(date +%T)"
$G grompp -f ../mdp/min.mdp -c ionized.gro -p system.top -n index.ndx -o min.tpr -maxwarn 5
nice -n 10 $MM run -p $GENV gmx mdrun -deffnm min -ntmpi 1 -ntomp 8 -nb gpu
echo ">>> [$SYS] MINIMIZATION DONE  $(date +%T)"
grep -E "Potential Energy|Maximum force|Norm of force|steps" min.log | tail -4
