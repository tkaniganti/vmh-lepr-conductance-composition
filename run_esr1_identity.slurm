#!/bin/bash
#SBATCH --job-name=vmh-e11
#SBATCH --account=bio260082p
#SBATCH --partition=RM-shared
#SBATCH --qos=rm
#SBATCH --cpus-per-task=8
#SBATCH --time=04:00:00
#SBATCH --output=logs/vmh-e11-%j.out
#SBATCH --error=logs/vmh-e11-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
"$PY" 29_esr1_1_identity.py
