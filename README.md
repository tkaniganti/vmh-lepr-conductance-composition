# VMH leptin-receptor conductance composition: reanalysis pipeline

Code accompanying *Conductance balance, not cell identity, as the basis of leptin-response
heterogeneity in VMH neurons* (Kaniganti, preprint).

Every reported statistic regenerates from the public accessions below. No bench data are
used; all inputs are deposited datasets from other laboratories.

## Data

| accession | contents | size |
|---|---|---|
| GSE172204 | Affinati et al. 2021, VMH single-nucleus RNA-seq, 10x trio plus `VMH_cluster` and `Conserved_cluster` metadata | 364 MB |
| GSE172206 | SF1 (*Nr5a1*) VMH TRAP-seq, three bead and three supernatant replicates | 722 KB |
| GSE172200 | *Lepr* VMH TRAP-seq, bead and supernatant by PBS and leptin, one sample each | 664 KB |
| GSE172460 | Rupp et al. 2023, reporter-defined hypothalamic LepR snRNA-seq, `Lepr_cluster` labels | 50 MB |

Three points that are easy to get wrong and are fixed in `config.py`:

1. The Rupp accession is **GSE172460**. `GSE172463` does not exist on GEO.
2. GSE172206 fraction assignment is **116053/054/055 bead, 116056/057/058 supernatant**,
   taken from `GSE172206_series_matrix.txt.gz` rather than inferred from column order.
3. GSE172200 is substantially arcuate in composition (*Pomc* enriched tenfold, *Agrp*
   sevenfold, *Nr5a1* 2.1-fold over its supernatant) and is used descriptively only.

## Running

```bash
bash 00_setup.sh            # LOGIN NODE, foreground, inside tmux: builds the conda env
bash 01_download.sh         # LOGIN NODE: compute nodes have no outbound network
sbatch run_analysis.slurm   # 02 through 07, RM-shared, about 30 min
sbatch run_scan.slurm       # 08, needs outputs from the above
# 09 to 19 run individually; 10, 11, 14 and 15 make a matrix pass and want a compute node
```

`env.sh` carries absolute paths, a conda prefix and a Slurm account specific to the machine
this was run on (Bridges-2 at the Pittsburgh Supercomputing Center). It is included for
provenance, not portability: edit `PROJECT`, `ENVPATH` and `ACCOUNT` before running elsewhere.

Environment constraints encoded in the scripts, which are specific to that system:
conda is activated by full prefix path rather than by name; Python is invoked as
`$ENVPATH/bin/python`; memory on RM-shared is requested by raising `--cpus-per-task`
because `--mem` is rejected; `-A` and `-q` are stated explicitly.

Count matrices are never loaded whole. Single `awk` passes over the gzipped matrix-market
files extract only the rows needed, and intermediates are cached in `data/`, so peak memory
stays near 1 to 2 GB.

## Scripts

Shared infrastructure:

| file | role |
|---|---|
| `config.py` | accessions, verified sample-to-fraction maps, gene panels, parameters |
| `vmhlib.py` | matrix streaming, control-matched module scores, bootstrap, permutation |
| `env.sh` | paths, conda activation, account |

Manuscript analyses:

| script | question | main output |
|---|---|---|
| `02_snrna_core.py` | population definition, detection, pseudobulk, per-cell module scores over 20 control-gene draws | `detection.json`, `cluster_scores.csv`, `percell_stats.json` |
| `03_trpc5_depth.py` | is the *Trpc5* gradient real or a capture-depth artifact | `trpc5_depth.json`, `trpc5_depth_null.csv` |
| `04_discreteness.py` | discrete subtypes or a continuum, with a power simulation | `discreteness.json` |
| `05_trap.py` | translatome cross-check with IP-over-supernatant enrichment | `trap_enrichment.csv`, `trap.json` |
| `06_crosspop.py` | placement against POMC and AgRP, contrasts computed within each control-gene draw | `crosspop_scores.csv`, `crosspop.json` |
| `08_gene_scan.py` | genome-wide scans for cluster structure and for correlates of effector balance | `gene_scan_cluster.csv`, `gene_scan_axis.csv` |
| `09_dispersion.py` | is the bidirectional population more dispersed than the arcuate references | `dispersion.json` |
| `10_conserved_classes.py` | variance partition between conserved classes against within the LepR class | `conserved_classes.json`, `conserved_mapping.csv` |
| `17_marker_classes.py` | do the transcriptomic classes reproduce the *in situ* marker composition | `marker_classes.json` |
| `18_doublet_sensitivity.py` | do the core results survive removal of author-flagged doublets | `doublet_sensitivity.json` |
| `19_manuscript_figures.py` | the four manuscript figures, every annotated statistic read from `results/` | `figures/manuscript_figures.pdf` |

`07_figures.py` produced an earlier figure set that the manuscript no longer uses. It is
retained because several analyses were developed against it.

## Scripts 11 to 16: a separate line of work

`11_ap_subunits.py` through `16_sigma1_trap.py` test adaptor protein complex subunit
composition across the same populations. That work is **not part of this manuscript** and
its results are reported elsewhere. The scripts are included because they share this
pipeline's infrastructure and cached intermediates, and because two of their findings bear
on the manuscript's methods: the median nucleus in which the AP-1 sigma position is detected
carries a single UMI across its three members, which is why per-cell compositional ratios are
not used anywhere in the manuscript, and `16` provides a worked example of the IP-over-
supernatant check applied to a composition rather than a single gene.

## Reproducibility notes

Control-gene selection for module scores is stochastic, so every module-based statistic is
computed over 20 independent draws and reported as mean with range rather than from a single
seed. Partial correlations residualise ranks on the rank of log depth; residualising on raw
depth leaves shared curvature that makes a depth-driven pair appear correlated.

The matrix-market header length is detected rather than assumed. An earlier version assumed
three header lines where the file has two, which silently dropped one matrix entry; the fix
is in `vmhlib.py` and `08_gene_scan.py`.

## Citation

See `CITATION.cff`. Licence: MIT (`LICENSE`).

## Variance partition and its controls

The structural result in the manuscript rests on scripts 20 to 25, added after the first
release. Scripts are listed in the order they were written; `20` and `21` are retained for
provenance but are superseded by `22` and `23`.

| script | question | main output |
|---|---|---|
| `20_variance_components.py` | mixed-model nested partition. **Superseded**: REML returned a between-class variance of exactly zero on every draw at six class levels, failing a positive control that a direct decomposition passes | `variance_components.json` |
| `21_variance_partition.py` | components of variance at each level separately. **Superseded** by `22`, since the three levels are nested and were not comparable as computed | `variance_partition.json` |
| `22_nested_partition.py` | proper two-level nested decomposition, with the detectability bound from injected cluster effects | `nested_partition.json` |
| `23_class_nulls.py` | per-class chance expectations, nuclei per cluster, leave-one-cluster-out, class-by-cluster membership | `class_nulls.json`, `class_membership.csv` |
| `24_mincluster_sensitivity.py` | do the per-class excesses survive excluding small clusters | `mincluster_sensitivity.json` |
| `25_supplementary_tables.py` | builds the six supplementary tables from the files above | `Supplementary_Tables.xlsx` |

Scripts 22 to 24 need a compute node and a few minutes each; slurm wrappers are provided.
`25` runs on a login node and needs `openpyxl`.

Two estimators appear in this project and should not be confused. Epsilon-squared from
Kruskal-Wallis, used in earlier scripts, rises with group count even in the absence of
structure and is reported only in Supplementary Table S2. The components-of-variance
estimator used in `22` to `24` is bias-corrected and truncated at zero, which is why
estimates can land exactly on zero.

The six supplementary tables are included under `supplementary/` so that reported values
can be checked without rerunning the pipeline.
