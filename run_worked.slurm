#!/bin/bash
#SBATCH -J vmh-wex
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-wex-%j.out
#SBATCH -e logs/vmh-wex-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 31_worked_example.py
