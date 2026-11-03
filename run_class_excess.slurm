#!/bin/bash
#SBATCH -J vmh-ces
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-ces-%j.out
#SBATCH -e logs/vmh-ces-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 58_class_excess_specificity.py
