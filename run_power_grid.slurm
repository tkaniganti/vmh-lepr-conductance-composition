#!/bin/bash
#SBATCH -J vmh-pgrid
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-pgrid-%j.out
#SBATCH -e logs/vmh-pgrid-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 51_power_grid.py
