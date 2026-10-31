#!/bin/bash
#SBATCH -J vmh-lpo
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-lpo-%j.out
#SBATCH -e logs/vmh-lpo-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 36_lepr_positive_outside.py
