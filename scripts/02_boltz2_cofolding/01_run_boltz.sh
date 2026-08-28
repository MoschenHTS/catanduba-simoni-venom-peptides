#!/usr/bin/env bash
# =============================================================================
# Boltz-2 run — C. simoni peptides (GPU: RTX 4080 in the original run).
#   * verified S1-S4 Kv4.2/Kv4.3 VSD receptors
#   * 50 diffusion samples, fixed seed (reproducible)
#   * positive benchmark (ProTx2/6N4I, DockQ-scored) + negative controls
#   * ensemble + interface-PAE evaluation, not best-of-N cherry-picking
# Run:  bash 01_run_boltz.sh
# =============================================================================

# Project paths. Set CSIMONI_BASE to a checked-out copy of this repository, and
# MICROMAMBA to your micromamba binary.
CSIMONI_BASE="${CSIMONI_BASE:-/path/to/this/repo}"
MICROMAMBA="${MICROMAMBA:-/path/to/micromamba}"

set -euo pipefail
cd "$(dirname "$0")"
BASE="$(pwd)"

# --- paths (all under $BASE, not the default home-directory locations) -------
export MAMBA_ROOT_PREFIX="$BASE/.mamba"
export PIP_CACHE_DIR="$BASE/.pipcache"
MM="$MICROMAMBA"
ENVDIR="$BASE/boltz_env"
CACHE="$BASE/.boltz_cache"        # Boltz weights/CCD (~7 GB) -> $BASE, NOT ~/.boltz
PY="$MM run -p $ENVDIR python"

# --- run params --------------------------------------------------------------
ACCEL=gpu
SAMPLES=50           # diffusion samples per input
RECYCLE=3
MAXPAR=5             # cap simultaneous diffusion samples to bound 16 GB VRAM
SEED=15              # reproducibility

# --- 1. environment ----------------------------------------------------------
# torch MUST be a CUDA-12 build (driver 535 / CUDA 12.2); cu130 wheels run on CPU.
if [ ! -x "$ENVDIR/bin/boltz" ]; then
  echo ">> creating boltz env at $ENVDIR"
  $MM create -y -p "$ENVDIR" -c conda-forge python=3.11
  $MM run -p "$ENVDIR" python -m pip install -U pip boltz
  $MM run -p "$ENVDIR" python -m pip install \
    --index-url https://download.pytorch.org/whl/cu121 "torch==2.5.1"
fi
# validation deps (biopython ships with boltz; tmtools optional for TM-align)
$MM run -p "$ENVDIR" python -c "import tmtools" 2>/dev/null || \
  $MM run -p "$ENVDIR" python -m pip install -q tmtools || echo ">> tmtools install skipped"

RUN="$MM run -p $ENVDIR boltz"
command -v nvidia-smi >/dev/null 2>&1 || { echo ">> no GPU -> CPU"; ACCEL=cpu; }
$PY -c "import torch; assert torch.cuda.is_available(), 'CUDA NOT available'; \
  print('>> GPU:', torch.cuda.get_device_name(0))" || \
  { echo '!! CUDA unavailable — fix torch before running (would be CPU-slow)'; exit 1; }

# --- 2. regenerate inputs (idempotent) --------------------------------------
$PY 02_make_inputs.py

# --- 3. predict --------------------------------------------------------------
# --use_msa_server fetches MSAs for BOTH chains; --write_full_pae for interface PAE.
mkdir -p out "$CACHE"
$RUN predict ./inputs \
  --accelerator "$ACCEL" \
  --use_msa_server \
  --output_format pdb \
  --out_dir ./out \
  --cache "$CACHE" \
  --seed "$SEED" \
  --diffusion_samples "$SAMPLES" \
  --max_parallel_samples "$MAXPAR" \
  --recycling_steps "$RECYCLE" \
  --no_kernels \
  --write_full_pae

# --- 4. evaluate + validate --------------------------------------------------
{
  $PY 03_evaluate_boltz.py
  echo; echo "########## POSITIVE BENCHMARK (ground truth: PDB 6N4I) ##########"
  $PY 04_validate_benchmark.py
  echo; echo "########## FOLD / DISULFIDE VALIDATION ##########"
  $PY 05_validate_fold.py
} | tee out/SUMMARY.txt

echo ""
echo ">> DONE. To hand off just the structures + scores for ChimeraX rendering"
echo "   elsewhere (skipping the multi-GB MSA/cache files):"
echo "   cd out && zip -r ../boltz_results.zip \\"
echo "       boltz_results_inputs/predictions/*/*.pdb \\"
echo "       boltz_results_inputs/predictions/*/confidence_*.json SUMMARY.txt"
