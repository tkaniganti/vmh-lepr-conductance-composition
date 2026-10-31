#!/bin/bash
#SBATCH -J vmh-pscw
#SBATCH -A bio260082p
#SBATCH -p RM-shared
#SBATCH -q rm
#SBATCH -c 8
#SBATCH -t 02:00:00
#SBATCH -o logs/vmh-pscw-%j.out
#SBATCH -e logs/vmh-pscw-%j.err

# Within-family Patch-seq calibration, Pvalb by default.
# Includes the positive control on the ephys ground-truth arm.

cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project
source ./env.sh
D=data/patchseq/mini-atlas/data
"$PY" 45_patchseq_within_family.py \
  --counts $D/m1_patchseq_exon_counts.csv.gz \
  --ttypes $D/m1_patchseq_meta_data.csv \
  --ephys  $D/m1_patchseq_ephys_features.csv \
  --ion-channels $D/ion-channel-genes-group-177.csv
