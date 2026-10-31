#!/bin/bash
#SBATCH -J vmh-cov
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-cov-%j.out
#SBATCH -e logs/vmh-cov-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 53_effector_covariation.py
