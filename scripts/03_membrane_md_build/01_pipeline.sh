#!/usr/bin/env bash
# Full membrane-MD pipeline for ONE complex: pdb2gmx -> orient -> embed -> build
# (solvate/ionize/minimize) -> equilibrate -> 3x200ns production.
# Usage: 01_pipeline.sh <system_dir> <n_tox_res> <vsd_type: Kv42|Kv43>

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -euo pipefail
SYS="$1"; NTOX="$2"; VSD="$3"
BASE=$CSIMONI_BASE/md
cd "$BASE"
MM="$MICROMAMBA"
export MAMBA_ROOT_PREFIX=$CSIMONI_BASE/.mamba
GENV=$CSIMONI_BASE/gmx_env
PENV=$CSIMONI_BASE/boltz_env

echo "=========================================================="
echo ">>> PIPELINE START: $SYS  (toxin=$NTOX res, VSD=$VSD)  $(date)"
echo "=========================================================="

if [ ! -f "$SYS/prot.gro" ]; then
  echo ">>> [$SYS] pdb2gmx (CHARMM36m)"
  ( cd "$SYS" && GMXLIB="$BASE" $MM run -p "$GENV" gmx pdb2gmx -f start.pdb -o prot.gro \
      -p topol.top -i posre.itp -ff charmm36-jul2022 -water tip3p -ignh -merge all )
fi

if [ ! -f "$SYS/prot_oriented.gro" ]; then
  echo ">>> [$SYS] orient in membrane frame"
  $MM run -p "$PENV" python 02_orient_place.py "$SYS" "$NTOX" "$VSD"
fi

NPROT=$(sed -n '2p' "$SYS/prot_oriented.gro" | tr -d ' ')
echo ">>> [$SYS] NPROT=$NPROT"

if [ ! -f "$SYS/embedded.gro" ]; then
  echo ">>> [$SYS] embed into POPC bilayer"
  $MM run -p "$PENV" python 03_embed.py "$SYS" "$NPROT"
fi

bash 07_build.sh "$SYS" "$NPROT"
bash 08_equil.sh "$SYS"
bash 09_prod.sh  "$SYS"

echo ">>> [$SYS] PIPELINE COMPLETE  $(date)"
