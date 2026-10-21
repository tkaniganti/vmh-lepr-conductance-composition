#!/usr/bin/env bash
#SBATCH --job-name=vmh-corr
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=08:00:00
#SBATCH --output=logs/vmh-corr-%j.out
#SBATCH --error=logs/vmh-corr-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh
echo "host=$(hostname)  started=$(date)"
"$PY" 27_correctness_checks.py
echo "finished=$(date)"
