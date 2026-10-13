#!/usr/bin/env bash
#SBATCH --job-name=vmh-scan
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/vmh-scan-%j.out
#SBATCH --error=logs/vmh-scan-%j.err

# Needs results/ and data/ from run_analysis.slurm (percell_scores.npz, colsum_vmh.txt,
# rowsum_vmh.txt). Scan A fits two logistic models per gene, so expect most of the wall
# time there; 8 cores gives about 16 GB, well above the dense arrays this builds.
set -euo pipefail
cd "$SLURM_SUBMIT_DIR"
source ./env.sh

for f in data/percell_scores.npz data/colsum_vmh.txt data/rowsum_vmh.txt; do
  [[ -s "$f" ]] || { echo "MISSING $f - run run_analysis.slurm first"; exit 1; }
done

echo "host=$(hostname)  started=$(date)"
"$PY" 08_gene_scan.py
echo "finished=$(date)"
