#!/usr/bin/env bash
#SBATCH --job-name=vmh-vp
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=02:00:00
#SBATCH --output=logs/vmh-vp-%j.out
#SBATCH --error=logs/vmh-vp-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh
echo "host=$(hostname)  started=$(date)"
"$PY" 21_variance_partition.py
echo "finished=$(date)"
