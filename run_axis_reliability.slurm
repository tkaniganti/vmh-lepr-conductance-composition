#!/bin/bash
#SBATCH -J vmh-rel
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 2
#SBATCH -t 00:30:00
#SBATCH -o logs/vmh-rel-%j.out
#SBATCH -e logs/vmh-rel-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 52_axis_reliability.py
