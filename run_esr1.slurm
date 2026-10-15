#!/usr/bin/env bash
#SBATCH --job-name=vmh-e1
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=03:00:00
#SBATCH --output=logs/vmh-e1-%j.out
#SBATCH --error=logs/vmh-e1-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh
echo "host=$(hostname)  started=$(date)"
"$PY" 28_esr1_1_markers.py
echo "finished=$(date)"
