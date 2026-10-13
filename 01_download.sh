#!/usr/bin/env bash
# 01_download.sh - fetch all public inputs. RUN ON A LOGIN NODE (compute nodes have no network).
# Idempotent: existing non-empty files are skipped. Total ~0.5 GB.
set -euo pipefail
cd "$(dirname "$0")"
source ./env.sh
cd "$DATA"

GEO=https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn

fetch () {  # $1 = url, $2 = outfile
  if [[ -s "$2" ]]; then
    echo "  [skip] $2"
  else
    echo "  [get ] $2"
    curl -fsSL --retry 3 --max-time 1800 -o "$2.part" "$1"
    mv -f "$2.part" "$2"
  fi
}

echo "== GSE172204: Affinati VMH snRNA-seq (364 MB matrix) =="
B=$GEO/GSE172204/suppl
fetch "$B/GSE172204_matrix.mtx.gz"   GSE172204_matrix.mtx.gz
fetch "$B/GSE172204_features.tsv.gz" GSE172204_features.tsv.gz
fetch "$B/GSE172204_barcodes.tsv.gz" GSE172204_barcodes.tsv.gz
fetch "$B/GSE172204_metadata.csv.gz" GSE172204_metadata.csv.gz

echo "== GSE172206: SF1(Nr5a1) VMH TRAP-seq, Bead B1-B3 + Sup S1-S3 =="
fetch "$GEO/GSE172206/suppl/GSE172206_counts.csv.gz" Nr5a1_TRAP.csv.gz
fetch "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE172nnn/GSE172206/matrix/GSE172206_series_matrix.txt.gz" \
      Nr5a1_TRAP_series_matrix.txt.gz

echo "== GSE172200: Lepr VMH TRAP-seq, Bead/Sup x PBS/Leptin =="
fetch "$GEO/GSE172200/suppl/GSE172200_counts.csv.gz"   Lepr_TRAP.csv.gz
fetch "$GEO/GSE172200/suppl/GSE172200_metadata.csv.gz" Lepr_TRAP_meta.csv.gz

echo "== GSE172460: Rupp reporter-defined hypothalamic LepR snRNA-seq (50 MB) =="
B=$GEO/GSE172460/suppl
fetch "$B/GSE172460_matrix.mtx.gz"   GSE172460_matrix.mtx.gz
fetch "$B/GSE172460_features.tsv.gz" GSE172460_features.tsv.gz
fetch "$B/GSE172460_barcodes.tsv.gz" GSE172460_barcodes.tsv.gz
fetch "$B/GSE172460_metadata.csv.gz" GSE172460_metadata.csv.gz

echo
set +o pipefail
echo "Sanity check of the two label columns the analysis depends on:"
echo "-- GSE172204 VMH_cluster --"
zcat GSE172204_metadata.csv.gz | head -1
echo "-- GSE172460 Lepr_cluster counts --"
zcat GSE172460_metadata.csv.gz | awk -F, 'NR>1 && $NF!="" {print $NF}' | sort | uniq -c | sort -rn
echo "-- GSE172206 fraction assignment --"
zcat Nr5a1_TRAP_series_matrix.txt.gz | grep -E "^!Sample_title|^!Sample_description|fraction:" | cut -c1-300

echo "Done. Files in $DATA"
