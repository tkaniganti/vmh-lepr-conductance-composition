#!/bin/bash
#SBATCH -J vmh-dbl
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-dbl-%j.out
#SBATCH -e logs/vmh-dbl-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 18_doublet_sensitivity.py
