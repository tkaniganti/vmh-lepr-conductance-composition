#!/bin/bash
#SBATCH -J vmh-cac
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-cac-%j.out
#SBATCH -e logs/vmh-cac-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 41_cluster_axis_corrected.py
