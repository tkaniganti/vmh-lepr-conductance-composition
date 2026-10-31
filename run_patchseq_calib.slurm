#!/bin/bash
#SBATCH -J vmh-psc
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -o logs/vmh-psc-%j.out
#SBATCH -e logs/vmh-psc-%j.err
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
D=data/patchseq/mini-atlas/data
"$PY" 42_patchseq_calibration.py \
  --counts $D/m1_patchseq_exon_counts.csv.gz \
  --ttypes $D/m1_patchseq_meta_data.csv \
  --ion-channels $D/ion-channel-genes-group-177.csv
