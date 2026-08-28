#!/usr/bin/env bash
# Equilibrate (NVT -> restrained NPT -> unrestrained NPT) and validate the membrane.
# Usage: 08_equil.sh <system_dir>

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -e
SYS="$1"
cd "$CSIMONI_BASE/md/$SYS"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=$CSIMONI_BASE/.mamba
GENV=$CSIMONI_BASE/gmx_env
G="$MM run -p $GENV gmx"
MD="nice -n 10 $MM run -p $GENV gmx mdrun -ntmpi 1 -ntomp 8 -nb gpu -pme gpu -cpt 5"
stage(){ local d=$1 c=$2 r=$3 t=$4 ex=""
  [ -f "$d.gro" ] && { echo ">> [$SYS] $d done"; return; }
  [ -n "$r" ] && ex="-r $r"; [ -n "$t" ] && ex="$ex -t $t"
  echo ">> [$SYS] $d  $(date +%T)"
  $G grompp -f ../mdp/$d.mdp -c "$c" -p system.top -n index.ndx -o "$d.tpr" $ex -maxwarn 5
  $MD -deffnm "$d" -cpi "$d.cpt"; }
stage nvt  min.gro  min.gro  ""
stage npt1 nvt.gro  min.gro  nvt.cpt
stage npt2 npt1.gro ""       npt1.cpt
echo ">> [$SYS] EQUILIBRATION DONE  $(date +%T)"
python3 - <<PY
L=open("npt2.gro").read().splitlines(); bx=L[-1].split()
nP=sum('POPC' in l[5:10] for l in L[2:-1])//134
print(f"[$SYS] final box: {float(bx[0]):.2f} x {float(bx[1]):.2f} x {float(bx[2]):.2f} nm")
print(f"[$SYS] area/lipid ~ {float(bx[0])*float(bx[1])*2/nP:.3f} nm^2 (POPC target 0.63-0.65; minus protein area)")
PY
