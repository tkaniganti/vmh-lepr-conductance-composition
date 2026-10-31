#!/bin/bash
#SBATCH -J vmh-ulr
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-ulr-%j.out
#SBATCH -e logs/vmh-ulr-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 38_upper_level_recovery.py
