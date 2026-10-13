#!/usr/bin/env bash
#SBATCH --job-name=vmh-lepr
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=06:00:00
#SBATCH --output=logs/vmh-lepr-%j.out
#SBATCH --error=logs/vmh-lepr-%j.err

# RM-shared caps memory at 2000M per core, and --mem is rejected on this account.
# Memory is requested by raising --cpus-per-task: 8 cores gives about 16 GB, which is
# well above the 1 to 2 GB these streaming passes need.
# -A and -q are stated explicitly because two accounts are active and the default
# account does not carry every QOS.
#
# Run 01_download.sh on a LOGIN NODE first. Compute nodes have no outbound network.

set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh

echo "host=$(hostname)  started=$(date)"
for f in data/GSE172204_matrix.mtx.gz data/GSE172460_matrix.mtx.gz \
         data/Nr5a1_TRAP.csv.gz data/Lepr_TRAP.csv.gz; do
  [[ -s "$f" ]] || { echo "MISSING $f - run 01_download.sh on a login node"; exit 1; }
done

"$PY" 02_snrna_core.py
"$PY" 03_trpc5_depth.py
"$PY" 04_discreteness.py
"$PY" 05_trap.py
"$PY" 06_crosspop.py
"$PY" 07_figures.py

echo "finished=$(date)"
ls -la "$RESULTS" "$FIGURES"
