#!/bin/bash
#SBATCH --job-name=ism_viz_build
#SBATCH --partition=htc
#SBATCH --qos=public
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem=48G
#SBATCH --time=01:30:00
#SBATCH --output=/scratch/tsingh65/ism-absorber-explorer/logs/ism_build_%j.out
#SBATCH --error=/scratch/tsingh65/ism-absorber-explorer/logs/ism_build_%j.err

set -euo pipefail
cd /scratch/tsingh65/ism-absorber-explorer
mkdir -p logs build_L012/rtml_50pc build_L012/ism_1kpc

PY="/scratch/tsingh65/envs/gal3d/bin/python"
if [ ! -x "$PY" ]; then
  PY="python3"
fi

echo "[$(date -Is)] Starting 3-level (L0,L1,L2) 256^3 ISM & Radiative Turbulent Mixing Layer deep-zoom build..."
"$PY" -m tngsight.ingest --preset mixing_layer --out build_L012/rtml_50pc --half-pc 25 --levels L0,L1,L2 --n 256 --nrays 120 --workers 16
"$PY" -m tngsight.vizpack validate build_L012/rtml_50pc/viz
"$PY" -m tngsight.ingest --preset ism_kpc --out build_L012/ism_1kpc --half-pc 500 --levels L0,L1,L2 --n 256 --nrays 120 --workers 16
"$PY" -m tngsight.vizpack validate build_L012/ism_1kpc/viz
echo "[$(date -Is)] Completed and validated 3-level L0,L1,L2 deep-zoom cubes!"
