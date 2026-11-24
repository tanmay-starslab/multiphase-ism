#!/bin/bash
#SBATCH --job-name=ism-6regimes
#SBATCH --partition=htc,public
#SBATCH --qos=public
#SBATCH --array=0-5
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/tsingh65/ism-absorber-explorer/slurm_ism_build_%A_%a.log

set -euo pipefail
cd /scratch/tsingh65/ism-absorber-explorer

export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

PYTHON="/scratch/tsingh65/envs/gal3d/bin/python"
if [ ! -x "$PYTHON" ]; then
  PYTHON="python3"
fi

VOLUMES=("RTML_50pc" "CloudCrush_200pc" "GMC_STARFORGE_20pc" "Superbubble_500pc" "ShatteringMist_10pc" "ISM_1kpc")
IDX="${SLURM_ARRAY_TASK_ID:-0}"
VOL="${VOLUMES[$IDX]}"

echo "=================================================================="
echo "[SLURM] Node: $(hostname) | CPUs: ${SLURM_CPUS_PER_TASK:-16} | Task: ${IDX} (${VOL}) | Start: $(date)"
echo "=================================================================="

"$PYTHON" scripts/build_ism_showcase.py \
  --only "$VOL" \
  --levels L0,L1,L2,L3,L4 \
  --n 256 \
  --nrays 120 \
  --workers "${SLURM_CPUS_PER_TASK:-16}" \
  --out-root /scratch/tsingh65/ism-absorber-explorer/build_5lvl

echo "=================================================================="
echo "[SLURM] Finished ${VOL} at $(date)"
echo "=================================================================="
