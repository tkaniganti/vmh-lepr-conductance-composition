# VMH LepR/SF1 leptin-effector reanalysis, Bridges-2 pipeline

Reanalysis of public data testing whether the bidirectional leptin response of SF1/LepR VMH
neurons reflects discrete effector-defined subtypes or a graded balance of co-expressed
conductances. All inputs are public; no bench data are used.

## Data

| accession | what | size |
|---|---|---|
| GSE172204 | Affinati VMH snRNA-seq, 10x trio + `VMH_cluster` metadata | 364 MB |
| GSE172206 | SF1(Nr5a1) VMH TRAP-seq, Bead B1-B3 + Sup S1-S3 | 722 KB |
| GSE172200 | Lepr VMH TRAP-seq, Bead/Sup x PBS/Leptin | 664 KB |
| GSE172460 | Rupp reporter-defined hypothalamic LepR snRNA-seq, `Lepr_cluster` labels | 50 MB |

Three corrections to the previous configuration, verified against GEO before writing:

1. The Rupp accession is **GSE172460**. `GSE172463` does not exist on GEO.
2. GSE172206 fraction assignment is **116053/054/055 = Bead, 116056/057/058 = Supernatant**
   (from `GSE172206_series_matrix.txt.gz`). The old positional slice `nr_hdr[2:5]` picked the
   three IP samples correctly, but by position rather than by design; it is now explicit in
   `config.py`.
3. Rupp raw label counts are L3.Nr5a1 424, L9.Pomc 231, L6.Agrp 189. The manuscript reports
   347 / 230 / 189, so `06_crosspop.py` prints the count at each filtering stage. Decide which
   filter the paper means and state it in Methods.

## Workflow

```bash
ssh kanigant@bridges2.psc.edu
tmux new -s vmh
cd /ocean/projects/bio260082p/kanigant/vmh-leptin-project

bash 00_setup.sh          # LOGIN NODE, foreground, small batches
bash 01_download.sh        # LOGIN NODE: compute nodes have no outbound network
sbatch run_analysis.slurm  # RM-shared, -c 8, about 16 GB, under 6 h
```

Interactive debugging:

```bash
srun -A bio260082p -p RM-shared -q rminteract -c 8 -t 02:00:00 --pty bash
source ./env.sh && "$PY" 03_trpc5_depth.py
```

Environment notes this pipeline respects:

- conda base is the shared system install; the hook is sourced and the env activated by full
  prefix path, never by name.
- Python is called as `$PY` (`$ENVPATH/bin/python`), never bare `python`.
- RM-shared rejects `--mem` (2000M per core cap), so memory is requested with
  `--cpus-per-task`. Eight cores is far more than the 1 to 2 GB these passes need.
- `-A` and `-q` are explicit, since two accounts are active and the default account does not
  carry every QOS.
- Scripts never overwrite with bare `>`, so `set -o noclobber` in `~/.bashrc` is harmless.
- The count matrices are never loaded whole. Single `awk` passes over the gzipped
  matrix-market files extract only the needed rows, and intermediates are cached in `data/`.

## Scripts

| script | what it does |
|---|---|
| `config.py` | accessions, verified sample-to-fraction maps, gene panels, parameters |
| `vmhlib.py` | matrix streaming, control-matched module scores, bootstrap, permutation |
| `02_snrna_core.py` | population definition, detection, pseudobulk, per-cell scores, KW over seeds |
| `03_trpc5_depth.py` | depth-controlled Trpc5 test: logistic regression, downsampling, empirical null |
| `04_discreteness.py` | GMM/BIC and dip test on the H-D axis, plus a power simulation |
| `05_trap.py` | TRAP CPM and IP-over-Sup enrichment for both series |
| `06_crosspop.py` | Rupp cross-population placement with bootstrap, permutation, leave-one-out |
| `07_figures.py` | all figures, with every annotated statistic read from results |

## What this changes in the manuscript

Fixes carried by `02`:

- the naive double-positive gate is reported separately for all VMH neurons and for LepR
  clusters, so the VMH-wide 344 can no longer be paired with the LepR denominator;
- labelled barcodes missing from `barcodes.tsv` and labelled nuclei with zero counts are
  counted and printed, which is where 6048/6049 and 2913/2914 came from;
- epsilon-squared is reported as mean and range over `N_SEEDS = 20` control-gene draws, so
  the ranges in the text have a stated origin.

New results that decide two of the manuscript's claims:

- `03` decides whether the Trpc5 gradient survives depth control. If the LR test, the
  downsampled chi-square, and the position of Trpc5 in the depth-matched null all hold, the
  claim can stay and the objection is answered in advance. If they do not, the gradient should
  come out of the paper, which costs little since the argument does not rest on it.
- `04` tests graded versus discrete directly rather than through cluster labels that were
  built without effector genes, and reports the size of split the pipeline could have
  detected. That bound belongs in the Limitations sentence in place of the current
  unquantified concession.
- `05` adds the enrichment the TRAP design provides, so "Kcnj11 is expressed" becomes a
  statement about enrichment over supernatant rather than bare IP abundance. It also prints
  the SF1-over-LepR Trpc5 ratio, which is the data-supported form of the Dhillon
  reconciliation currently argued by assertion in the Discussion.

## Still open

- Whether Allison et al. 2015 is actually reanalyzed. The scope note says yes, Methods says
  no. Nothing here uses it.
- The repository URL and a Zenodo DOI for the data availability statement.
