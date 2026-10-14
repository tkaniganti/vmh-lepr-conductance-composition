#!/usr/bin/env bash
#SBATCH --job-name=vmh-st
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=00:20:00
#SBATCH --output=logs/vmh-st-%j.out
#SBATCH --error=logs/vmh-st-%j.err
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh
"$PY" -c "import openpyxl" 2>/dev/null || "$ENVPATH/bin/pip" install -q openpyxl
"$PY" 25_supplementary_tables.py
