#!/usr/bin/env bash
# MM-GBSA (single-trajectory) on the last 50ns of each production replica, for
# all 4 complexes. 200 frames/replica (interval=5 of the 1000 available frames).
# Uses the patched gmx_MMPBSA (OT1/OT2 CHARMM C-terminal naming fix applied to
# GMXMMPBSA/make_top.py fixparm2amber - see notes) and ff14SB (not ff99SB, which
# fails on these termini regardless of the OT1/OT2 fix - ff14SB's templates are
# more complete). No GAFF needed (pure protein-protein, no small molecule).

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
G="$MM run -p $GENV gmx"

cat > mmpbsa_prod.in <<'EOF'
Production GB-MMPBSA, last 50ns, 200 frames
&general
startframe=1, endframe=1000, interval=5,
verbose=1,
forcefields="leaprc.protein.ff14SB"
/
&gb
igb=5, saltcon=0.150,
/
EOF

run_one () {
  local SYS="$1" NTOX="$2"
  cd "$BASE/$SYS"
  local NPROT; NPROT=$(sed -n '2p' prot_oriented.gro | tr -d ' ')
  echo "=========================================================="
  echo ">>> [$SYS] NPROT=$NPROT NTOX=$NTOX  $(date +%T)"
  echo "=========================================================="
  mkdir -p mmpbsa
  if [ ! -f mmpbsa/tv_index_full.ndx ]; then
    python3 - "$NPROT" "$NTOX" <<PYEOF
NPROT, NTOX = $NPROT, $NTOX
L = open("npt2.gro").read().splitlines()
n = int(L[1]); atoms = L[2:2+n]
resid = [int(a[0:5]) for a in atoms[:NPROT]]
boundary = next(i for i in range(1, NPROT) if resid[i] == 1 and resid[i-1] != 1)
tox_idx = list(range(1, boundary+1))
vsd_idx = list(range(boundary+1, NPROT+1))
with open("mmpbsa/tv_index_full.ndx", "w") as f:
    f.write("[ TOXIN ]\n")
    for j in range(0, len(tox_idx), 15): f.write(" ".join(map(str, tox_idx[j:j+15])) + "\n")
    f.write("[ VSD ]\n")
    for j in range(0, len(vsd_idx), 15): f.write(" ".join(map(str, vsd_idx[j:j+15])) + "\n")
print(f"[$SYS] TOXIN={len(tox_idx)} VSD={len(vsd_idx)}")
PYEOF
  fi
  for rep in 1 2 3 4 5 6 7 8; do
    local OUT="mmpbsa/rep${rep}_FINAL_RESULTS_MMPBSA.dat"
    [ -f "$OUT" ] && { echo "[$SYS] rep$rep done, skip"; continue; }
    echo "[$SYS] rep$rep MM-GBSA  $(date +%T)"
    cd mmpbsa
    rm -rf COM_* RECEPTOR_* LIGAND_* _GMXMMPBSA_* gmx_MMPBSA.log leap.log 2>/dev/null
    $MM run -p "$MMENV" gmx_MMPBSA -O -i ../../mmpbsa_prod.in -cs "../prod_rep${rep}.tpr" \
        -ci tv_index_full.ndx -cg 1 0 -ct "../prod_rep${rep}.xtc" -nogui \
        > "rep${rep}_run.log" 2>&1
    if [ -f FINAL_RESULTS_MMPBSA.dat ]; then
      mv FINAL_RESULTS_MMPBSA.dat "rep${rep}_FINAL_RESULTS_MMPBSA.dat"
      echo "[$SYS] rep$rep OK"
    else
      echo "[$SYS] rep$rep FAILED - see mmpbsa/rep${rep}_run.log"
    fi
    cd ..
  done
}

run_one cstx2_kv42 32
run_one cstx3_kv42 33
run_one cstx2_kv43 32
run_one cstx3_kv43 33

echo ">>> MM-GBSA ALL DONE  $(date)"
