# VMH leptin-receptor conductance composition: reanalysis pipeline

Code accompanying *Cell identity does not resolve bidirectional leptin responses in VMH neurons: what single-nucleus data can and cannot measure* (Kaniganti, preprint).

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
provenance, not portability: edit `ENVPATH` and `ACCOUNT` before running elsewhere; `PROJECT` resolves itself.

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
| `19_manuscript_figures.py` | the five manuscript figures, every annotated statistic read from `results/` | `figures/manuscript_figures.pdf` |

`07_figures.py` produced an earlier figure set that the manuscript no longer uses. It is
retained because several analyses were developed against it.

## Scripts 11 to 16: a separate line of work

`11_ap_subunits.py` through `16_sigma1_trap.py` test adaptor protein complex subunit
composition across the same populations. That work is **not part of this manuscript** and
its results are reported elsewhere. The scripts are included because they share this
pipeline's infrastructure and cached intermediates, and because `16` provides a worked example of the IP-over-supernatant check applied to a composition rather than a single gene.

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
| `25_supplementary_tables.py` | builds the seven supplementary tables from the files above | `Supplementary_Tables.xlsx` |

Scripts 22 to 24 need a compute node and a few minutes each; slurm wrappers are provided.
`25` runs on a login node and needs `openpyxl`.

Two estimators appear in this project and should not be confused. Epsilon-squared from
Kruskal-Wallis, used in earlier scripts, rises with group count even in the absence of
structure and is reported only in Supplementary Table S2. The components-of-variance
estimator used in `22` to `24` is bias-corrected and truncated at zero, which is why
estimates can land exactly on zero.

The seven supplementary tables are included under `supplementary/` so that reported values
can be checked without rerunning the pipeline.

## Regenerating the supplementary tables

`25_supplementary_tables.py` reads from `results/` and does not compute anything itself, so
the files it needs must exist first. Tables S2 and S4 to S6 come from scripts 10 and 22 to
24; table S1 needs `18_doublet_sensitivity.py` and table S3 needs `09_dispersion.py`. Running
only 22 to 25 produces a workbook with those two sheets empty.

## Verifying a clean checkout

`run_verify.slurm` clones the released tag into a fresh directory, reruns scripts 22 to 25
against the cached downloads, and checks the results against the values reported in the
manuscript: five draw-averaged magnitudes, four claims stated as inequalities against each
class's own chance expectation, and five ordering and leave-one-out structures. It edits
paths in the slurm header and in `BASE` for the machine it runs on.


## A note on PROJECT

`env.sh` resolves `PROJECT` to the directory it lives in, falling back to
`SLURM_SUBMIT_DIR` and then the working directory, and it refuses to continue if the
result does not contain `vmhlib.py`. An inherited `PROJECT` is deliberately ignored, because that variable is in common use for other purposes and an inherited value silently redirected several runs during development.

`run_verify.slurm` writes its own log. It is not shipped in the release, because a log committed before a tag is cut necessarily describes a different tag than the one it ships in; the script regenerates it in about three minutes against whichever tag it names.

## Correctness and identity checks

Scripts 27 and 29 support claims made in the manuscript's Methods. The others are
exploratory and are not cited in the paper.

| script | question | main output |
|---|---|---|
| `26_insr_check.py` | is there a transcriptomic handle on insulin-responding VMH neurons | `insr_check.json` |
| `27_correctness_checks.py` | parser verified against the raw matrix; random-gene panels; full label shuffle | `correctness_checks.json` |
| `28_esr1_1_markers.py` | what distinguishes Esr1_1 from the rest of its class | `esr1_1_markers.csv` |
| `29_esr1_1_identity.py` | is Esr1_1 correctly assigned to the conserved Esr1 class | `esr1_1_identity.json` |
| `30_pomc_replication.py` | the unsupervised structure test applied to arcuate POMC, as calibration | `pomc_replication.json` |

Check 27 asks whether the pipeline is right rather than whether it reruns. The parser is
compared entry by entry against the deposited matrix-market file without the analysis
library; the effector panel is replaced by random genes to test whether the class-level
component is specific to conductances; and all labels are permuted to confirm the
decomposition reports nothing when the groupings carry no information. The pass criteria
are positional, comparing observed values against the distributions they should exceed,
with no threshold chosen by hand.

`vmhlib.pass_rows` keys its cache on the set of rows requested, not on the caller's name
for them. An earlier version keyed on the name alone, so a caller that reused a name for a
different set of genes silently received the previous extraction. That affected only
`27_correctness_checks.py`, whose random-panel loops reuse names across changing gene sets;
every other script uses one name per tag and was unaffected, which a cold-cache rerun of the
full pipeline confirmed.

Check 27 also runs the random-panel substitution within the leptin-receptor clusters, which
tests the circularity objection at the level the central claim is made, and reports observed
values as percentiles of their null distributions rather than against a chosen threshold.

`env.sh` no longer honours an inherited `PROJECT`. It always resolves to its own
directory, because that variable is in common use for other purposes and an inherited
value silently redirected several runs during development.

`25_supplementary_tables.py` now builds seven sheets; S7 reads `pomc_replication.json`.
