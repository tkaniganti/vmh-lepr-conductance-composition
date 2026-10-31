#!/bin/bash
#SBATCH -J vmh-smn
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-smn-%j.out
#SBATCH -e logs/vmh-smn-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 55_scan_max_null.py
