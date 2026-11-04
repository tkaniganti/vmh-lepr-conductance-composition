#!/usr/bin/env python3
"""
34_manifest.py - map every number the manuscript quotes to the object that produced it.

A review of v49 found seven numerical problems. Six were the same failure: correct
arithmetic applied to the wrong object. An observed value quoted where the text meant an
excess over chance. A twenty-draw mean in one paragraph and a single draw in another, both
called the same thing. Two genome-wide scans with different detection filters, both
reported as "the strongest correlate". Two independently drawn nulls, 300 shuffles and 200,
producing excesses that differed by 0.1 and read as a thresholding effect.

None of those is a computational error and none would fail a verification run, because the
code was right throughout. What was wrong was the pairing between a sentence and a value.

This script makes that pairing explicit. MANIFEST below declares, for every figure the
Results quote, the value as written, the file and key path it comes from, and the procedure
that produced it. The script then checks each declaration against the deposited results and
reports three things: declarations whose stored value no longer matches, declarations whose
key path has gone missing, and stored quantities that nothing in the manifest claims.

The maintenance cost is real: a new number in Results needs a new MANIFEST line. That is
the same cost as the figure-callout audit, which has caught enough to justify itself.

What a passing run does NOT establish. This verifies that a quoted number equals the value
at a declared key. It cannot verify that the key is the right one. A declaration pointing at
a plausibly named field whose contents are the wrong side of a filter passes whenever the
two happen to agree. detection.json contains exactly that case: n_vmh_labelled and
n_vmh_with_counts both hold 6,048, the first named as though it were the pre-filter count,
so a declaration naming either passes and only one is sourced correctly.

The procedure field carries the sourcing claim and nothing checks it. Write it as a
description of what the field contains, not a restatement of its name, and read a passing
run as evidence that the numbers match rather than that the claims are sourced.

Usage:
    python3 34_manifest.py                 check the manifest against results/
    python3 34_manifest.py --emit          print the manifest as a table for the supplement
    python3 34_manifest.py --orphans       list stored values nothing claims
    python3 34_manifest.py --scope         which producers the verification reaches

Output: printed, plus results/manifest_check.json when run without flags.
"""
import os, sys, json, csv
import config as cfg

R = cfg.RESULTS
TOL = 0.005          # a quoted figure may be rounded; this is the allowed gap


# ---------------------------------------------------------------------------
# MANIFEST
#
# (label, quoted, file, key path, procedure, locations)
#
# locations is optional and lists the sections that quote this figure, so one
# declaration serves several quoted locations rather than being repeated. A key path
# declared twice is reported as a duplicate: two declarations of one quantity can drift
# apart while both still pass.
#
# "quoted" is the number as it appears in Results, rounded as written.
# "procedure" is the thing that was wrong six times out of seven: which of several
# similar quantities this is. Be specific. "mean over 20 control-gene draws" and
# "single draw, seed = RANDOM_SEED" are different procedures even where both are
# correct.
# ---------------------------------------------------------------------------
MANIFEST = [
    # ---- 2.1 detection and gating
    ("VMH nuclei analysed", 6048, "detection.json", ["n_vmh_with_counts"],
     "VMH-labelled barcodes that also have a matrix column; n_vmh_labelled holds the pre-filter 6,049", ["2.1", "2.2"]),
    ("neurons analysed", 21582, "detection.json", ["n_neurons"],
     "Neuron_cluster non-empty AND present in barcodes.tsv; n_neurons_labelled holds the pre-filter 21,585"),
    ("Lepr-cluster nuclei analysed", 2913, "detection.json", ["n_lepr_with_counts"],
     "VMH_cluster beginning Lepr AND present in barcodes.tsv", ["2.1", "2.2", "2.5", "2.6"]),
    ("Lepr-cluster nuclei labelled", 2914, "detection.json", ["n_lepr_labelled"],
     "NOT QUOTED. VMH_cluster beginning Lepr in the metadata, before the barcode filter; the manuscript quotes 2,913, the post-filter count"),
    ("VMH nuclei labelled", 6049, "detection.json", ["n_vmh_labelled"],
     "VMH_cluster set in the metadata, before the barcode filter"),
    ("neurons labelled", 21585, "detection.json", ["n_neurons_labelled"],
     "Neuron_cluster set in the metadata, before the barcode filter"),
    ("VMH barcodes absent", 1, "detection.json", ["n_vmh_labelled_not_in_barcodes"],
     "labelled VMH rows with no column in barcodes.tsv"),
    ("neuron barcodes absent", 3, "detection.json", ["n_neurons_not_in_barcodes"],
     "labelled neuron rows with no column in barcodes.tsv"),
    ("Nr5a1 detected, all VMH", 24.2, "detection.json", ["vmh_detect", "Nr5a1"],
     "percent of VMH nuclei with at least one UMI"),
    ("Lepr detected, all VMH", 14.4, "detection.json", ["vmh_detect", "Lepr"],
     "percent of VMH nuclei with at least one UMI"),
    ("Slc17a6 detected, all VMH", 76.6, "detection.json", ["vmh_detect", "Slc17a6"],
     "percent of VMH nuclei with at least one UMI"),
    ("Nr5a1 detected, Lepr clusters", 33.5, "detection.json", ["lepr_detect", "Nr5a1"],
     "percent of Lepr-cluster nuclei with at least one UMI"),
    ("Lepr detected, Lepr clusters", 23.5, "detection.json", ["lepr_detect", "Lepr"],
     "percent of Lepr-cluster nuclei with at least one UMI"),
    ("gate retained, n", 316, "detection.json", ["naive_gate", "lepr_clusters"],
     "Nr5a1 and Lepr both detected, full set"),
    ("gate retained, pct", 10.8, "detection.json", ["naive_gate", "lepr_pct"],
     "percent of 2,913, full set"),

    # ---- 2.1 doublet profile
    ("flagged nuclei", 461, "doublet_sensitivity.json", ["n_lepr_doublet_flagged"],
     "Conserved_cluster beginning Doublets, within Lepr clusters"),
    ("gate, doublet-free n", 282, "doublet_sensitivity.json",
     ["without_doublets", "naive_gate_lepr"], "gate applied after excluding flagged"),
    ("median depth, flagged", 13095, "doublet_sensitivity.json",
     ["doublet_profile", "median_depth_flagged"], "median total counts"),
    ("median depth, other", 10837, "doublet_sensitivity.json",
     ["doublet_profile", "median_depth_other"], "median total counts"),
    ("Lepr detected, flagged", 18.2, "doublet_sensitivity.json",
     ["doublet_profile", "detection_pct_flagged", "Lepr"], "percent of flagged nuclei"),
    ("Lepr detected, other", 24.5, "doublet_sensitivity.json",
     ["doublet_profile", "detection_pct_other", "Lepr"], "percent of unflagged nuclei"),
    ("Nr5a1 detected, flagged", 33.2, "doublet_sensitivity.json",
     ["doublet_profile", "detection_pct_flagged", "Nr5a1"], "percent of flagged nuclei"),
    ("Lepr outside the Lepr clusters", 5.9, "doublet_sensitivity.json",
     ["doublet_profile", "lepr_pct_outside_lepr_clusters"],
     "pooled across the five non-Lepr VMH classes"),

    # ---- 2.2 per-cell effect sizes
    ("eps2, combined axis", 0.023, "percell_stats.json", ["HminusD", "eps2_mean"],
     "MEAN over 20 control-gene draws, not a single draw"),
    ("eps2, hyperpolarising", 0.0075, "percell_stats.json", ["HYPER", "eps2_mean"],
     "mean over 20 draws"),
    ("eps2, depolarising", 0.0228, "percell_stats.json", ["DEPOL", "eps2_mean"],
     "mean over 20 draws"),
    ("largest p, hyperpolarising", 0.0495, "percell_stats.json", ["HYPER", "p_max"],
     "MAXIMUM p across 20 draws: below 0.05 means every draw was significant"),

    # ---- 2.2 conserved classes
    ("conserved Lepr class size", 1857, "conserved_classes.json",
     ["class_sizes", "Lepr"], "all nuclei in the conserved Lepr class", ["2.2", "2.x within-class"]),
    ("pct of Lepr-cluster nuclei in the class", 75.5, "conserved_classes.json",
     ["pct_lepr_nuclei_in_lepr_class"],
     "DENOMINATOR is Lepr-cluster nuclei carrying a non-doublet conserved "
     "assignment (2,453), not 2,913 and not the 2,452 doublet-free set"),

    # ---- 2.2 nested partition
    ("cluster-in-class component", 4.73, "nested_partition.json",
     ["nested_all_vmh", "pct_cluster_in_class", "mean"], "mean over 20 draws"),
    ("nuclei in the partition", 5441, "nested_partition.json", ["n_nuclei"],
     "NOT QUOTED. VMH nuclei carrying a non-doublet conserved class; an intermediate the manuscript does not state"),

    # ---- 2.2 per-class excesses. Observed and excess are DIFFERENT quantities.
    ("Esr1 excess", 15.65, "class_nulls.json", ["by_class", "Esr1", "excess_over_null"],
     "EXCESS: observed minus its own null"),
    ("Esr1 observed", 15.898, "class_nulls.json",
     ["by_class", "Esr1", "observed_pct", "mean"], "OBSERVED, before the null is removed"),
    ("Nfib excess", 6.29, "class_nulls.json", ["by_class", "Nfib", "excess_over_null"],
     "EXCESS, null from 300 shuffles in script 23"),
    ("Nfib observed", 6.654, "class_nulls.json",
     ["by_class", "Nfib", "observed_pct", "mean"], "observed"),
    ("Lepr excess", -0.03, "class_nulls.json", ["by_class", "Lepr", "excess_over_null"],
     "EXCESS, negative because the estimator is bias-corrected"),
    ("Lepr observed", 0.095, "class_nulls.json",
     ["by_class", "Lepr", "observed_pct", "mean"],
     "OBSERVED: do not quote this where the text means the excess"),
    ("Lepr null", 0.122, "class_nulls.json", ["by_class", "Lepr", "null_pct", "mean"],
     "chance expectation for this class's k and n"),
    ("Foxp2 excess", 0.13, "class_nulls.json", ["by_class", "Foxp2", "excess_over_null"],
     "excess"),
    ("Dlk1 excess", 2.94, "class_nulls.json", ["by_class", "Dlk1", "excess_over_null"],
     "excess"),
    ("Fezf1 excess", 2.69, "class_nulls.json", ["by_class", "Fezf1", "excess_over_null"],
     "excess. Numerically close to the class component that earlier versions reported and unrelated to it; that component is no longer reported"),

    # ---- 2.2 threshold series. Its nulls are drawn separately from script 23's.
    ("Nfib at threshold 0", 6.19, "mincluster_sensitivity.json",
     ["by_class", "Nfib", "0", "excess"],
     "SERIES null, 200 shuffles in script 24; differs from the 6.29 above because the "
     "nulls are independent draws, not because of thresholding"),
    ("Nfib at threshold 20", 6.89, "mincluster_sensitivity.json",
     ["by_class", "Nfib", "20", "excess"], "series, 200 shuffles"),

    # ---- 2.4 power grid (51_power_grid.py).
    #
    # NON-SATURATED CELL PERCENTAGES ARE NOT DECLARED, and must not be quoted in the
    # manuscript either. Each cell is 100 simulations with its own seed, so anything
    # not at a true 0% or 100% moves between runs. Reseeding per cell on 29 Oct moved
    # balanced_2 at 2.0 SD from 100 to 99, minority_2 at 2.0 from 93 to 98,
    # sohn_3_skew at 2.5 from 100 to 97, and the dip at balanced_2 2.5 from 99 to 100.
    # The qualitative pattern did not move at all: every threshold below was identical
    # across both runs. So the THRESHOLDS are the reproducible quantity and the cell
    # percentages are not.
    #
    # The three configurations that recover only above 3 SD (sohn_3, sohn_3_asym,
    # uneven_3) are likewise not declared: thresholds_80pct stores them as strings
    # flagged "mechanism check only", and they are not detectability estimates.
    ("recovery threshold, balanced split", 2.0, "power_grid.json",
     ["thresholds_80pct", "balanced_2", "bic"],
     "lowest separation in SD units at which BIC selected k >= 2 in at least 80 of 100 "
     "simulations. Reproduced identically across two independent seedings"),
    ("recovery threshold, 18% minority", 2.0, "power_grid.json",
     ["thresholds_80pct", "minority_2", "bic"],
     "same, single 18% minority group; equal to the balanced threshold, which is how "
     "2.4 establishes that group size is not the constraint"),
    ("recovery threshold, skewed three-group", 2.5, "power_grid.json",
     ["thresholds_80pct", "sohn_3_skew", "bic"],
     "same, 18/18/64 displaced -0.2 and +0.8. The cell that locates the edge"),
    ("dip threshold, balanced split", 2.5, "power_grid.json",
     ["thresholds_80pct", "balanced_2", "dip"],
     "the dip test needs more separation than mixture selection, because a 50/50 "
     "equal-variance mixture is unimodal below 2 sigma"),
    ("widening at the resolving separation", 1.8, "power_grid.json",
     ["widening_if_resolvable", "by_separation", "5.0", "ratio_to_observed"],
     "ratio of the mixture SD to the observed SD for the sohn_3 configuration at the "
     "separation where it first recovers. Law of total variance over two 18% groups "
     "at +/- sep/2 and a 64% CENTRE GROUP AT ZERO, every component carrying the "
     "observed SD of 0.4515. THE CENTRE COMPONENT IS THE POINT: omitting it gives "
     "about 2.7 instead of 1.8, and two independent reconstructions from the "
     "manuscript's own figures landed 50% apart for exactly that reason, which is "
     "why this is computed in 51_power_grid.py rather than left to the reader. "
     "Single-draw quantity: observed_sd comes from one control-gene draw"),
    ("nuclei in the power grid", 2913, "power_grid.json", ["n"],
     "NOT QUOTED. Lepr-cluster nuclei the grid injects into"),
    ("observed SD of the axis", 0.4515, "power_grid.json", ["observed_sd"],
     "NOT QUOTED. Every separation is expressed in these units"),
    # ---- 2.8 TRAP. The BH correction is computed in 05_trap.py as of v1.1.6; before
    # that it existed nowhere in the project and the q values had no producer.
    ("BH family, panel", 12, "trap.json", ["bh_families", "n_panel_tested"],
     "panel genes with a finite Welch p, of 22 in HYPER and DEPOL"),
    ("BH family, panel total", 22, "trap.json", ["bh_families", "n_panel_total"],
     "HYPER and DEPOL only; MARKERS is not part of the family"),
    ("BH family, measured", 18, "trap.json", ["bh_families", "n_measured_tested"],
     "every gene tested with a finite Welch p, the wider family"),
    ("Nr5a1 enrichment", 3.84, "trap_enrichment.csv", ["Nr5a1", "SF1_log2_IP_over_Sup"],
     "log2 immunoprecipitate over supernatant, SF1 fraction"),
    ("Kcnj5 enrichment", 1.15, "trap_enrichment.csv", ["Kcnj5", "SF1_log2_IP_over_Sup"],
     "log2 IP over Sup"),
    ("Kcnj5 q", 0.0264, "trap_enrichment.csv", ["Kcnj5", "SF1_q_bh_panel"],
     "BH over the 12-gene panel family; manuscript rounds to 0.026"),
    ("Kcnj3 q", 0.0284, "trap_enrichment.csv", ["Kcnj3", "SF1_q_bh_panel"],
     "BH over the panel family; manuscript rounds to 0.028"),
    ("Hcn1 q", 0.047, "trap_enrichment.csv", ["Hcn1", "SF1_q_bh_panel"],
     "BH over the panel family"),
    ("Abcc8 q", 0.047, "trap_enrichment.csv", ["Abcc8", "SF1_q_bh_panel"],
     "BH over the panel family"),
    ("Pde3b q", 0.0264, "trap_enrichment.csv", ["Pde3b", "SF1_q_bh_panel"],
     "depleted, not enriched; BH over the panel family"),
    ("Gck enrichment", 1.44, "trap_enrichment.csv", ["Gck", "SF1_log2_IP_over_Sup"],
     "log2 IP over Sup"),

    # ---- 2.6 effector-panel cluster ranks. gene_scan_cluster.csv is truncated at
    # TOP_N = 300, so only Hcn1 and Trpc5 appear there; every other panel gene must be
    # read from gene_scan.json. Reading absence from the CSV as absence of structure is
    # what put Abcc8 on the wrong list in an earlier draft.
    ("Hcn1 cluster rank", 243, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Hcn1", "rank"],
     "rank by LR for cluster given log depth, of 12,361 tested"),
    ("Trpc5 cluster rank", 275, "gene_scan.json", ["scanA_trpc5", "rank"],
     "rank by LR for cluster given log depth"),
    ("Trpc5 percentile", 97.8, "gene_scan.json", ["scanA_trpc5", "percentile"],
     "percentile of the same ranking"),
    ("Abcc8 cluster rank", 1068, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Abcc8", "rank"],
     "rank by LR for cluster given log depth, of 12,361 tested; structured at "
     "q = 0.0, and quoted in 2.6 alongside the five panel-gene percentiles"),
    ("Abcc8 percentile", 91.4, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Abcc8", "percentile"], "percentile"),
    ("Kcnj3 percentile", 93.0, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Kcnj3", "percentile"], "percentile"),
    ("Kcnj5 percentile", 93.2, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Kcnj5", "percentile"], "percentile"),
    ("Kcnj6 percentile", 82.7, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Kcnj6", "percentile"],
     "percentile; sets the LOWER bound of the range stated in 2.6, which is "
     "why it must be declared alongside the other four"),
    ("Kcnj9 percentile", 85.9, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Kcnj9", "percentile"], "percentile"),
    ("Pde3b q", 0.5528, "gene_scan.json",
     ["scanA_effector_panel_ranks", "Pde3b", "q_bh"],
     "no cluster structure after depth correction"),
    ("genes tested, cluster scan", 12361, "gene_scan.json", ["n_genes_tested"],
     "adequately detected genes in scan A"),
    ("genes tested, axis scan", 12182, "gene_scan.json", ["scanB_n_tested"],
     "scan B, after excluding panel and control genes"),
    ("genes at 5% FDR, axis scan", 277, "gene_scan.json", ["scanB_n_sig_q05"],
     "scan B significant at q < 0.05"),

    # ---- 2.5 within-population scan
    ("strongest partial correlation", 0.1371, "gene_scan_axis.csv", ["Gm42418", "partial_rho"],
     "scan B of script 08, transcripts detected in 5% to 95% of Lepr nuclei; Rbfox1 at "
     "0.143 is excluded by the upper bound and comes from script 32, a sensitivity check"),

    # ---- 2.2 cluster-level axis, corrected measure (28 Oct 2026). None of these was
    # ever declared, which is why nothing caught that they came from cluster_scores.csv
    # (log1p pseudobulk, no z-scoring, no control correction) while the text called them
    # control-corrected per-gene SD scores. Recomputed by 41_cluster_axis_corrected.py
    # from per-nucleus module scores averaged within cluster over 20 draws.
    ("axis range, lowest cluster", -0.3215, "cluster_axis_corrected.json",
     ["axis_range_all_clusters", 0],
     "min over the 24 VMH cluster means of per-nucleus HminusD; per-gene SD units, "
     "NOT log-CP10k; the cluster is Esr1_2"),
    ("axis range, highest cluster", 0.1752, "cluster_axis_corrected.json",
     ["axis_range_all_clusters", 1], "max over the same 24 cluster means; Nfib_3"),
    ("Lepr span, lowest", -0.0925, "cluster_axis_corrected.json",
     ["axis_range_lepr_clusters", 0],
     "min over the nine Lepr-cluster means; negative, so these clusters are not "
     "confined to the hyperpolarizing side as earlier versions stated"),
    ("Lepr span, highest", 0.1582, "cluster_axis_corrected.json",
     ["axis_range_lepr_clusters", 1], "max over the same nine cluster means"),
    ("Esr1_1 axis value", 0.0925, "cluster_axis_corrected.json", ["esr1_1"],
     "the Esr1_1 cluster mean; the other three Esr1 clusters are at -0.3215, -0.2835 "
     "and -0.2000, stored under by_cluster"),
    ("depolarizing panel spread", 0.4477, "cluster_axis_corrected.json",
     ["DEPOL_spread", "span"],
     "max minus min of the DEPOL module score across the 24 cluster means. A panel "
     "score, not a position on the axis; it runs -0.1740 to 0.2737"),
    ("hyperpolarizing panel spread", 0.1684, "cluster_axis_corrected.json",
     ["HYPER_spread", "span"],
     "max minus min of the HYPER module score across the same means; -0.0717 to 0.0967"),

    # ---- 2.2 class contrasts, Figure 2C (declared 28 Oct 2026). Section 2.2 rests on
    # these now that the class component is withdrawn, and not one was declared before
    # tonight. The manuscript quotes them at three decimals; Nfib and Dlk1 differ only
    # at the fourth, so they are separate paths with separate values.
    ("Lepr vs Esr1 contrast", 0.1364, "conserved_classes.json",
     ["lepr_class_vs_others", "Lepr_vs_Esr1", "diff_mean_over_seeds"],
     "mean over 20 control-gene draws of the conserved Lepr class minus the "
     "ventrolateral Esr1 class on the balance axis; positive is more hyperpolarizing. "
     "Sign consistent across all 20 draws"),
    ("Lepr vs Fezf1 contrast", -0.0733, "conserved_classes.json",
     ["lepr_class_vs_others", "Lepr_vs_Fezf1", "diff_mean_over_seeds"],
     "same construction; negative, so Lepr is more depolarizing than Fezf1"),
    ("Lepr vs Nfib contrast", -0.0644, "conserved_classes.json",
     ["lepr_class_vs_others", "Lepr_vs_Nfib", "diff_mean_over_seeds"],
     "same construction; quoted as -0.064 in 2.2, distinct from the Dlk1 contrast "
     "only at the fourth decimal"),
    ("Lepr vs Dlk1 contrast", -0.0643, "conserved_classes.json",
     ["lepr_class_vs_others", "Lepr_vs_Dlk1", "diff_mean_over_seeds"],
     "same construction; also quoted as -0.064 in 2.2"),
    ("Lepr vs Foxp2 contrast", 0.0384, "conserved_classes.json",
     ["lepr_class_vs_others", "Lepr_vs_Foxp2", "diff_mean_over_seeds"],
     "same construction; the smallest of the five, and the only other positive one"),

    # ---- 2.2 reliability and 2.5's ceiling (52_axis_reliability.py, 29 Oct 2026).
    # Every DERIVED value here is computed in that script from the emitted inputs, not
    # written in: the ceiling is sqrt of the count-split reliability, the bound
    # conversions are the injected levels over that square root, and the cluster
    # reliabilities come from the variance split with cluster sizes read from
    # cluster_axis_corrected.json. A key holding 0.511 would be checkable only against
    # itself; computed from 0.261 it makes this manifest check the chain.
    ("count split, half correlation", 0.1501, "axis_reliability.json",
     ["count_split", "half_rho_mean"],
     "mean Spearman between two halves of each nucleus's molecules, 20 binomial "
     "thinnings; range 0.1272 to 0.1766, which is narrower than independent draws "
     "because the halves come from the same counts"),
    ("per-nucleus reliability", 0.261, "axis_reliability.json",
     ["count_split", "spearman_brown_full_depth"],
     "Spearman-Brown on the half correlation. The quantity 2.2 quotes and 2.3, 2.4 "
     "and 2.5 are read against; it appears in five places in the text and this is "
     "the key all five refer to"),
    ("panel split, lower end", -0.0131, "axis_reliability.json",
     ["panel_split", "half_rho_min"],
     "lowest of 20 stratified half-panel correlations. The range CONTAINS ZERO, which "
     "is why the text gives a range here and a point estimate for the count split"),
    ("panel split, upper end", 0.0623, "axis_reliability.json",
     ["panel_split", "half_rho_max"], "highest of the same 20"),
    ("attenuation ceiling", 0.5109, "axis_reliability.json",
     ["derived", "attenuation_ceiling"],
     "sqrt of the per-nucleus reliability: the largest correlation this axis can show "
     "with a perfectly measured transcript. 2.5's whole argument is the comparison of "
     "the observed 0.143 against this"),
    ("variance split, capture noise", 73.9, "axis_reliability.json",
     ["derived", "variance_split_pct", "capture_noise"],
     "percent of per-nucleus variance that is capture noise, as (1 - reliability)"),
    ("variance split, between cluster", 2.8, "axis_reliability.json",
     ["derived", "variance_split_pct", "between_cluster"],
     "size-weighted variance of the 24 cluster means as a percent of the total. NOT "
     "the 3.57% nested cluster-within-class component, which is a different object"),
    ("variance split, within cluster", 23.3, "axis_reliability.json",
     ["derived", "variance_split_pct", "within_cluster"],
     "the remainder; counts AGAINST a cluster mean's reliability, not toward it"),
    ("reliability, smallest cluster", 0.644, "axis_reliability.json",
     ["derived", "cluster_mean_reliability", "smallest_cluster", "reliability"],
     "cluster mean as an estimate of that cluster's value relative to others, at the "
     "smallest of the 24. The size is read from the data, so this follows the "
     "distribution rather than a literal in the text"),
    ("smallest cluster, n", 63, "axis_reliability.json",
     ["derived", "cluster_mean_reliability", "smallest_cluster", "n"],
     "nuclei in the smallest VMH cluster, from cluster_axis_corrected.json"),
    ("reliability at 200 nuclei", 0.8517, "axis_reliability.json",
     ["derived", "cluster_mean_reliability", "at_200", "reliability"], "same, at n=200"),
    ("reliability, conserved class", 0.9816, "axis_reliability.json",
     ["derived", "cluster_mean_reliability", "conserved_class", "reliability"],
     "same, at the 1,857 nuclei of the conserved Lepr class"),
    ("bound conversion, 0.10", 0.196, "axis_reliability.json",
     ["derived", "detectability_bound_conversion", "underlying_sd_units", 0],
     "2.3's injected 0.10 of the OBSERVED per-nucleus spread, expressed in standard "
     "deviations of the underlying quantity: 0.10 over sqrt(reliability)"),
    ("bound conversion, 0.20", 0.391, "axis_reliability.json",
     ["derived", "detectability_bound_conversion", "underlying_sd_units", 1],
     "the same for the 0.20 level"),
    ("reliability floor control", -0.0081, "axis_reliability.json",
     ["shuffled_floor"],
     "two INDEPENDENTLY permuted count matrices; must be near zero. An earlier version "
     "permuted once and scored that matrix twice, returning 0.92 against a ceiling of "
     "0.93, which is what the ceiling measures and not a floor"),
    ("reliability ceiling control", 0.9288, "axis_reliability.json",
     ["control_gene_draw_ceiling"],
     "same genes and counts under two control-gene draws: the most any split can reach"),

    # ---- 2.2 and 2.3 random-panel controls (27_correctness_checks.py).
    # These five values were quoted in the manuscript with NOTHING declared from this
    # file, so the scope check could not see that 27_correctness_checks.py was absent
    # from run_verify.slurm: it derives its source list from declarations, and a file
    # nothing declares from is invisible to it however much the paper takes from there.
    # Methods also attributed these controls to the nested decomposition, which reads
    # as script 22; script 22 has no random-panel control at all.
    ("random panels, class component mean", 0.6758, "correctness_checks.json",
     ["check_B_random_panels", "class_mean"],
     "mean class component over 100 detectability-matched random panels, scored by the "
     "same decomposition as the effector panel. Both arms carry the same upward bias, "
     "which is why the comparison holds where neither value stands alone"),
    ("random panels, class component median", 0.4535, "correctness_checks.json",
     ["check_B_random_panels", "class_median"], "median of the same 100"),
    ("observed class component, this control", 1.5946, "correctness_checks.json",
     ["check_B_random_panels", "observed_class"],
     "the effector panel's class component in the SINGLE control-gene draw used for "
     "this control; not the 20-draw value, and not comparable to it"),
    ("random panels, within-class mean", 0.3936, "correctness_checks.json",
     ["check_B2_within_class_random_panels", "random_mean"],
     "same substitution scored within the nine Lepr clusters, where the effector panel "
     "returns zero"),
    ("random panels, within-class median", 0.2677, "correctness_checks.json",
     ["check_B2_within_class_random_panels", "random_median"],
     "median within-class excess over 100 random panels matched to the effector panel "
     "for detectability, scored within the nine leptin-receptor clusters. Fixed by "
     "rng_b2 at RANDOM_SEED + 2000 with control-gene draws at RANDOM_SEED + 900 + i, "
     "so a rerun reproduces it; without that seed a legitimate rerun would move the "
     "third decimal for a reason no check can tell from staleness. The manuscript "
     "quoted 0.272 for many versions, which this pipeline does not produce: the key "
     "did not exist until 29 Oct and the origin of 0.272 is unknown. It passed "
     "44_manuscript_check because 0.2677 rounds to 0.27 and the tolerance absorbed a "
     "stale value rather than a prose rounding"),

    # ---- 2.3 and 2.7, the three load-bearing values the triage surfaced
    ("Nfib nuclei at the 20-nucleus threshold", 807, "mincluster_sensitivity.json",
     ["by_class", "Nfib", "20", "n_nuclei"],
     "nuclei retained in Nfib after cells below 20 are dropped; the denominator that "
     "makes the Nfib-versus-Lepr comparison readable on both sides"),
    ("Hcn1 TRAP p", 0.019, "trap_enrichment.csv", ["Hcn1", "SF1_welch_p_unadjusted"],
     "Welch p for Hcn1, IP over supernatant, SF1 fraction; unadjusted, with the BH q "
     "declared separately"),

    # ---- 2.9 Esr1_1 identity (29_esr1_1_identity.py).
    # Nine values the section quotes, none declared until 30 Oct. Like
    # 27_correctness_checks.py, this script was absent from run_verify.slurm and the
    # scope check could not report that, because it derives its source list from
    # declarations and nothing declared from this file. Found by the completeness
    # check only after its consumer search was fixed to stop counting a script's
    # mentions of its own output as consumption.
    ("Esr1_1 conserved-assignment agreement", 85.4, "esr1_1_identity.json",
     ["test1_annotation", "focus_matching_pct"],
     "246 of 288 nuclei of the cluster carry the matching conserved assignment. The "
     "denominator is the cluster AS ANNOTATED, which includes the one nucleus the "
     "deposited Conserved_cluster column marks Doublets_2; tests 2, 4 and 5 below work "
     "on 287 because that nucleus is dropped at the class assignment. The two counts "
     "are different populations by design, not a discrepancy, and the claim here is "
     "about annotation agreement across the cluster. Quoted in prose from the two "
     "counts until the ratio was emitted on 30 Oct"),
    ("Esr1_1 nuclei, as annotated", 288, "esr1_1_identity.json",
     ["test1_annotation", "focus_n_total"],
     "NOT QUOTED. The denominator of the agreement ratio above"),
    ("Esr1_1 pseudobulk r to its own class", 0.9868, "esr1_1_identity.json",
     ["test2_nearest_class", "correlations", "Esr1"],
     "transcriptome-wide correlation of the cluster's pseudobulk profile with the "
     "conserved Esr1 class; quoted as 0.987"),
    ("Esr1_1 pseudobulk r to next nearest", 0.981, "esr1_1_identity.json",
     ["test2_nearest_class", "correlations", "Fezf1"],
     "the next nearest class. The margin over the assigned class is 0.006, which the "
     "section states as its own limitation"),
    ("Satb2 detection in Esr1_1", 71.1, "esr1_1_identity.json",
     ["test5_markers", "Satb2", "FOCUS"],
     "percent of the cluster's nuclei detecting Satb2"),
    ("Satb2 detection in the Esr1 class", 23.4, "esr1_1_identity.json",
     ["test5_markers", "Satb2", "Esr1"], "the same for the rest of the class"),
    ("Esr1 detection in Esr1_1", 37.6, "esr1_1_identity.json",
     ["test5_markers", "Esr1", "FOCUS"],
     "the cluster is Esr1-low relative to its classmates, which the section reports "
     "rather than explains"),
    ("Esr1 detection in the Esr1 class", 77.2, "esr1_1_identity.json",
     ["test5_markers", "Esr1", "Esr1"], "the same for the rest of the class"),
    ("doublet rate in Esr1_1", 0.35, "esr1_1_identity.json",
     ["test3_doublets", "focus_pct"],
     "quoted as 0.3%, which is Python's round-half-to-even on 0.35; a reader rounding "
     "half away from zero gets 0.4%. The producers round to even and 39 and 44 were "
     "switched to half-up, so the two halves of the pipeline disagree at a 5 boundary"),
    ("doublet rate in the rest of the class", 0.21, "esr1_1_identity.json",
     ["test3_doublets", "others_pct"], "quoted as 0.2%"),

    # ---- 2.2 effector covariation (53_effector_covariation.py).
    ("covariation contrast, cluster level", 0.0638, "effector_covariation.json",
     ["conditions", "all genes, all clusters", "contrast"],
     "mean correlation between genes of the same panel minus the mean between panels, "
     "over depth-residualized pseudobulk profiles across the 24 clusters. THE PRIMARY "
     "CONDITION: all 22 genes, all clusters. The other three are pre-declared "
     "sensitivities and the strongest of them should not lead"),
    ("covariation percentile, cluster level", 96.5, "effector_covariation.json",
     ["conditions", "all genes, all clusters", "percentile_of_null"],
     "where that contrast falls among 200 detectability-matched gene sets. 94.0 "
     "without Pde3b and 98.0 with the three sub-10% genes removed, both declared below"),
    ("covariation contrast, per nucleus", 0.0, "effector_covariation.json",
     ["per_nucleus", "contrast"],
     "the same statistic on per-nucleus scores; quoted as below 0.0001"),
    ("covariation null p95, per nucleus", 0.0546, "effector_covariation.json",
     ["per_nucleus", "null_contrast_p95"],
     "95th percentile of the per-nucleus null; quoted as 0.055"),
    ("covariation at the matched loading", 0.5168, "effector_covariation.json",
     ["per_nucleus", "injected_at_matched_loading"],
     "per-nucleus contrast after injecting a common factor sized to reproduce the "
     "observed CLUSTER-level contrast, solved by bisection from a permuted baseline. "
     "Matched to the effect size, NOT asserted to be the same factor. Quoted as 0.517"),
    ("covariation percentile, without Pde3b", 94.0, "effector_covariation.json",
     ["conditions", "without Pde3b", "percentile_of_null"],
     "pre-declared sensitivity; falls inside the conventional band, which 2.2 states"),
    ("covariation percentile, detection >= 10%", 98.0, "effector_covariation.json",
     ["conditions", "detection >= 10%, all clusters", "percentile_of_null"],
     "pre-declared sensitivity, dropping Hcn4, Kcnj11 and Kcnj5"),

    # ---- Discussion, per-nucleus distribution (43_axis_distribution.py).
    # SINGLE CONTROL-GENE DRAW. The file records draws = 1, and every value below is
    # from that one draw, where 2.2 and 2.3 report ranges over 20. The manuscript does
    # not say so. Its sd, 0.4515, is the same quantity power_grid.json stores as
    # observed_sd, so a change to one must move the other.
    ("axis minimum", -1.78, "axis_distribution.json", ["min"],
     "lowest per-nucleus score across the 2,913 Lepr-cluster nuclei; the Discussion "
     "quotes the range as -1.78 to 1.77"),
    ("axis maximum", 1.7685, "axis_distribution.json", ["max"],
     "highest; quoted as 1.77"),
    ("fraction below zero", 47.0, "axis_distribution.json", ["pct_below_zero"],
     "quoted as 47.0%. Zero is where the two panels sit equally far above their own "
     "matched controls, not a membrane potential, so this describes the distribution "
     "and does not assign nuclei to response classes"),
    ("mean, in SD units from zero", 0.0199, "axis_distribution.json",
     ["mean_in_sd_units_from_zero"], "quoted as 0.02 standard deviations from zero"),
    ("axis mean", 0.009, "axis_distribution.json", ["mean"],
     "NOT QUOTED. The numerator of the figure above"),

    # ---- ranges that the old section-number pattern hid.
    # Until 30 Oct, 39_undeclared_numbers.py dropped every token shaped
    # ^\\d{1,2}\\.\\d{1,2}$ as a section number. That is a systematic bias rather than a
    # random one: a RANGE ENDPOINT is far likelier than a mean to be a short decimal,
    # so the check was blind to ranges specifically. In both cases below the central
    # value was declared and its range was not, and the range is what a reader checks.
    ("cluster-in-class component, low draw", 3.7085, "nested_partition.json",
     ["nested_all_vmh", "pct_cluster_in_class", "min"],
     "lowest of the 20 control-gene draws; quoted as 2.80 beside the declared mean "
     "of 3.57"),
    ("cluster-in-class component, high draw", 5.6571, "nested_partition.json",
     ["nested_all_vmh", "pct_cluster_in_class", "max"], "highest of the same 20"),
    ("Trpc5 detection, lowest cluster", 30.1, "trpc5_depth.json",
     ["downsampling", "pct_detected_raw_min"],
     "lowest per-cluster Trpc5 detection across the nine Lepr clusters, Lepr_4. The "
     "RANGE is emitted rather than its endpoints declared from the per-cluster dict: "
     "a later run moving Lepr_5 from 30.2 to 30.0 would leave 30.1 stored and the "
     "quoted floor wrong, and the declaration would still pass"),
    ("Trpc5 detection, highest cluster", 68.0, "trpc5_depth.json",
     ["downsampling", "pct_detected_raw_max"], "highest, Lepr_2"),

    # ---- the class level after the estimator correction (2 Nov 2026).
    # nested() used expected-mean-square coefficients for a balanced design at
    # both levels; see its docstring in 22_nested_partition.py. The cluster term
    # was understated and the class term inflated. Everything below is from the
    # corrected function, and the class component is now reported rather than
    # withheld because the reason for withholding it was the bias, not the design.
    ("class component", 1.6293, "nested_partition.json",
     ["nested_all_vmh", "pct_class", "mean"],
     "share of per-nucleus variance at the conserved-class level, mean over 20 "
     "control-gene draws. The superseded formula gave 2.6884; about 40% of that "
     "was cluster variance leaking through the wrong coefficient"),
    ("class component, low draw", 0.836, "nested_partition.json",
     ["nested_all_vmh", "pct_class", "min"], "lowest of the 20 draws"),
    ("class component, high draw", 2.3907, "nested_partition.json",
     ["nested_all_vmh", "pct_class", "max"], "highest of the 20"),
    ("class component, null 95th", 2.142, "nested_partition.json",
     ["class_component_null", "null_p95"],
     "95th percentile of 200 permutations that reassign whole CELLS to classes, "
     "preserving cells per class and every cell size. The observed 1.6293 does "
     "not exceed it"),
    ("class component, percentile of null", 92.0, "nested_partition.json",
     ["class_component_null", "percentile_of_null"],
     "where the observed value falls in that null. Quoted rather than a verdict, "
     "because 92 is below the conventional line and close to it, which is a "
     "different situation from a value sitting inside its null"),
    ("random panels reaching the component", 0.08, "correctness_checks.json",
     ["check_B_random_panels", "fraction_of_panels_at_or_above_observed"],
     "fraction of 100 detectability-matched arbitrary panels producing a class "
     "component at least as large as the effector panel's. Under the superseded "
     "formula this was 0.03 and the check passed; it now fails"),
    ("random panels, class 95th", 2.0521, "correctness_checks.json",
     ["check_B_random_panels", "class_p95"],
     "95th percentile of that random distribution, above the observed 1.5946"),

    # ---- contrast specificity (57_contrast_specificity.py).
    ("contrast specificity, joint", 24, "contrast_specificity.json",
     ["joint_largest_of_five"],
     "of 100 detectability-matched panels, how many reproduced the LARGEST of the "
     "five class contrasts. Criterion fixed before the run, so that one panel "
     "passing on one contrast cannot be read as passing on all five"),
    ("contrast specificity, Esr1", 5, "contrast_specificity.json",
     ["panels_reaching_observed", "Lepr_vs_Esr1"],
     "signed: a panel counts only if its difference is at least as large in the "
     "same direction. Two of five contrasts at 5 is close to what five tries give "
     "by chance, so the lowest two are the least reproduced rather than specific"),
    ("contrast specificity, Fezf1", 5, "contrast_specificity.json",
     ["panels_reaching_observed", "Lepr_vs_Fezf1"], "same"),
    ("contrast specificity, Nfib", 10, "contrast_specificity.json",
     ["panels_reaching_observed", "Lepr_vs_Nfib"],
     "22 without the depth correction, which was applied after the uncorrected "
     "result was seen; Nfib is the one class with appreciably lower median depth"),
    ("contrast specificity, Dlk1", 25, "contrast_specificity.json",
     ["panels_reaching_observed", "Lepr_vs_Dlk1"], "same"),
    ("contrast specificity, Foxp2", 31, "contrast_specificity.json",
     ["panels_reaching_observed", "Lepr_vs_Foxp2"], "same"),
    ("Lepr vs Esr1, single draw", 0.1423, "contrast_specificity.json",
     ["observed_single_draw", "Lepr_vs_Esr1"],
     "the contrast on the one draw the substitution test scores against, NOT the "
     "20-draw mean of 0.136 that 2.2 reports alongside it"),
    ("Lepr vs Esr1, labels permuted", 0.002, "contrast_specificity.json",
     ["control_zero_permuted_labels", "Lepr_vs_Esr1"],
     "zero control: permuting class labels collapses the largest contrast from "
     "0.1423 to this, so the differences do belong to the class partition even "
     "though they are not specific to effector genes"),


    # ---- the full label shuffle (27_correctness_checks.py, check_C).
    # The third null on the class component and the only one it clears. 2.2 quoted
    # 0.082 and 0.330 for the cluster term from the superseded estimator; both
    # moved. Declared at the new key paths, which did not exist in the manifest.
    ("full shuffle, cluster term mean", 0.1085, "correctness_checks.json",
     ["check_C_full_shuffle", "cluster_in_class", "shuffle_mean"],
     "200 shuffles of class and cluster labels across all nuclei, destroying every "
     "grouping. Was 0.082 under the superseded formula"),
    ("full shuffle, cluster term 95th", 0.438, "correctness_checks.json",
     ["check_C_full_shuffle", "cluster_in_class", "shuffle_p95"],
     "95th percentile of the same; was 0.330. No shuffle reached the observed "
     "4.6514 on that single draw"),
    ("full shuffle, class term mean", 0.023, "correctness_checks.json",
     ["check_C_full_shuffle", "class", "shuffle_mean"],
     "the same null at the class level, which the manuscript did not previously "
     "report. The class component clears this one and not the other two: it is "
     "distinguishable from no structure at all, but neither from an arbitrary "
     "regrouping of the same cells nor from an arbitrary gene panel"),
    ("full shuffle, class term 95th", 0.1446, "correctness_checks.json",
     ["check_C_full_shuffle", "class", "shuffle_p95"],
     "95th percentile; the observed 1.5946 is far above it, and no shuffle "
     "reached the observed value"),


    # ---- per-class specificity (58_class_excess_specificity.py).
    # Judged against the 95th percentile, the same rule applied to the class
    # component and the contrasts. Two of the four classes above the floor are
    # specific and two are not; the threshold was NOT moved after the counts were
    # seen. Foxp2 and Lepr sit at zero, where a one-sided count is vacuous, so they
    # are untestable this way rather than unanimous at 100 of 100.
    #
    # THE TESTED STATISTIC IS NOT THE REPORTED ONE. This is the raw one-way
    # component on a single depth-residualized draw; 2.3's figures are excesses over
    # per-class nulls averaged over draws. Esr1 reads 16.474 here and 15.65 there.
    # The counts transfer; the values do not.
    ("class specificity, Esr1", 1, "class_excess_specificity.json",
     ["panels_reaching_observed", "Esr1"],
     "of 100 detectability-matched panels. Esr1's component rests on the single "
     "Esr1_1 cell, so this is specificity of a single-cell effect"),
    ("class specificity, Nfib", 4, "class_excess_specificity.json",
     ["panels_reaching_observed", "Nfib"],
     "the one class carrying DISTRIBUTED structure that is also specific to the "
     "effector panel; this is what the title's within-class clause rests on"),
    ("class specificity, Fezf1", 6, "class_excess_specificity.json",
     ["panels_reaching_observed", "Fezf1"],
     "does not clear the 95th percentile. Fezf1's component rests on Lepr_7"),
    ("class specificity, Dlk1", 20, "class_excess_specificity.json",
     ["panels_reaching_observed", "Dlk1"],
     "the other distributed-structure class, and not specific. Named rather than "
     "averaged with the rest"),
    ("class specificity, Esr1 observed", 16.474, "class_excess_specificity.json",
     ["observed_per_class", "Esr1"],
     "NOT QUOTED as a result. The raw single-draw component the count above is "
     "scored against, recorded so the two statistics are not confused"),


    # ---- the between-cluster term, corrected 3 Nov 2026.
    # THREE ERRORS, all inflating it. The cluster means came from
    # cluster_axis_corrected.json (24 VMH clusters, 6,048 nuclei) while the total
    # came from this script's own axis over the 2,913 Lepr-cluster nuclei, so the
    # ratio divided one population by another; the two also came from different
    # control-gene draws; and the averaging sizes were the 24 clusters'. Means are
    # now computed from the script's own axis over its own nine clusters.
    #
    # The raw weighted variance of cluster means also carries their sampling noise.
    # With W = total - true, E[raw] = true + k*W/N, so
    # true = (raw - k*total/N)/(1 - k/N). Unlike the raw statistic this can return
    # a negative value under no structure, which is what makes a zero control mean
    # anything.
    ("between-cluster variance, raw", 0.006689, "axis_reliability.json",
     ["derived", "variance_split_pct", "between_raw"],
     "NOT QUOTED. Size-weighted variance of the nine Lepr cluster means before the "
     "noise correction"),
    ("between-cluster variance, corrected", 0.006039, "axis_reliability.json",
     ["derived", "variance_split_pct", "between_noise_corrected"],
     "NOT QUOTED as a variance; quoted as 2.8% of the total"),
    ("between control, zero injection", -5.2e-05, "axis_reliability.json",
     ["derived", "variance_split_pct", "between_control", 0, "recovered_mean"],
     "NOT QUOTED. The control that the raw statistic could not pass: under no "
     "between-cluster structure the corrected estimator returns zero, at 1.4 "
     "standard errors below it"),
    ("between control, realized at 0.005", 0.004182, "axis_reliability.json",
     ["derived", "variance_split_pct", "between_control", 2, "realized_mean"],
     "NOT QUOTED. Recovery is judged against what nine clusters actually realize, "
     "which is about 85% of the nominal injected parameter, not against the "
     "parameter itself. Recovered 0.004154 against this, within 0.1 SE"),
    ("reliability, median cluster", 0.8963, "axis_reliability.json",
     ["derived", "cluster_mean_reliability", "median_cluster", "reliability"],
     "at 301 nuclei, the median of the nine Lepr clusters"),
    ("reliability, largest cluster", 0.9487, "axis_reliability.json",
     ["derived", "cluster_mean_reliability", "largest_cluster", "reliability"],
     "at 644 nuclei"),


    # ---- confidence intervals and the two Trpc5 population scores in 2.9.
    # These have been in the manuscript since well before 3 Nov and were never
    # declared. They are quoted correctly and the values are stored; the gap is
    # simply that the declaration pass covered the point estimates and not the
    # intervals beside them. Found by reading the sentence, not by any check:
    # a value that is quoted and undeclared is invisible to the manifest, and
    # the undeclared listing shows it among a hundred others.
    ("POMC balance axis, CI low", -0.161, "crosspop_scores.csv",
     ["POMC", "HminusD_ci_lo"],
     "bootstrap lower bound on the arcuate POMC balance-axis score of -0.097"),
    ("POMC balance axis, CI high", -0.0349, "crosspop_scores.csv",
     ["POMC", "HminusD_ci_hi"], "upper bound; the text rounds this to -0.035"),
    ("POMC Trpc5 score", 0.9422, "crosspop_scores.csv", ["POMC", "Trpc5_mean"],
     "arcuate POMC Trpc5 module score in log-normalized units; quoted as 0.942. "
     "NOT a fold change, for the reason given for the 0.656 difference"),
    ("POMC Trpc5 CI low", 0.8343, "crosspop_scores.csv", ["POMC", "Trpc5_ci_lo"],
     "quoted as 0.834"),
    ("POMC Trpc5 CI high", 1.0563, "crosspop_scores.csv", ["POMC", "Trpc5_ci_hi"],
     "quoted as 1.056"),


    ("AgRP Trpc5 score", 0.2862, "crosspop_scores.csv", ["AgRP", "Trpc5_mean"],
     "the other half of the comparison the 0.656 difference is computed from; "
     "quoted as 0.286, and like the POMC value it is a log-normalized module "
     "score, not a fold change"),
    ("AgRP Trpc5 CI low", 0.2066, "crosspop_scores.csv", ["AgRP", "Trpc5_ci_lo"],
     "quoted as 0.207"),
    ("AgRP Trpc5 CI high", 0.3743, "crosspop_scores.csv", ["AgRP", "Trpc5_ci_hi"],
     "quoted as 0.374"),


    # ---- the VMH_SF1 row of crosspop_scores.csv, quoted in 2.9.
    # The last two of the five intervals the paper states. Found by grepping the
    # document for the interval pattern, which is the only way a quoted-and-
    # undeclared value can surface: the manifest reports on what it was told, so
    # a value nobody declared is outside its comparison entirely.
    ("VMH depolarizing module", -0.0431, "crosspop_scores.csv",
     ["VMH_SF1", "DEPOL_mean"],
     "the module that carried the VMH population's position; quoted as -0.043, "
     "and negative only in this population"),
    ("VMH depolarizing module, CI low", -0.0709, "crosspop_scores.csv",
     ["VMH_SF1", "DEPOL_ci_lo"], "quoted as -0.071"),
    ("VMH depolarizing module, CI high", -0.0149, "crosspop_scores.csv",
     ["VMH_SF1", "DEPOL_ci_hi"], "quoted as -0.015"),
    ("VMH Trpc5 score", 0.4614, "crosspop_scores.csv",
     ["VMH_SF1", "Trpc5_mean"],
     "intermediate between POMC at 0.9422 and AgRP at 0.2862; quoted as 0.461"),
    ("VMH Trpc5 CI low", 0.4022, "crosspop_scores.csv",
     ["VMH_SF1", "Trpc5_ci_lo"], "quoted as 0.402"),
    ("VMH Trpc5 CI high", 0.522, "crosspop_scores.csv",
     ["VMH_SF1", "Trpc5_ci_hi"], "quoted as 0.522"),
    ("VMH SF1 nuclei, cross-population", 424, "crosspop_scores.csv",
     ["VMH_SF1", "n"],
     "after the UMI filter, the population placed against POMC and AgRP"),


    # ---- the expected-mean-square coefficients Methods quotes.
    # They were computed, returned by nested(), and then dropped: the
    # nested_counts dict copied three keys and not these. They became
    # load-bearing when the 3 Nov reorganization moved the sentence that quotes
    # them into a paragraph-opening position, which the position triage counts.
    # The ratio between them is the whole reason the class term cannot be formed
    # by subtracting MS cluster-within-class, so they carry the argument.
    ("EMS coefficient c1", 61.502, "nested_partition.json",
     ["nested_counts", "ems_c1"],
     "the coefficient on the cluster component in E(MS cluster-within-class); "
     "quoted as 61.5"),
    ("EMS coefficient c2", 253.915, "nested_partition.json",
     ["nested_counts", "ems_c2"],
     "the coefficient on the cluster component in E(MS class); quoted as 253.9. "
     "The superseded formula subtracted MS cluster-within-class directly, which "
     "assumes this equals c1; here they differ by a factor of four"),
    ("EMS coefficient c3", 861.823, "nested_partition.json",
     ["nested_counts", "ems_c3"],
     "NOT QUOTED. The divisor on the class component, unchanged by the "
     "correction: the superseded n0a was already identically equal to it"),

]


def load(fn):
    p = os.path.join(R, fn)
    if not os.path.exists(p):
        return None
    if fn.endswith(".csv"):
        return list(csv.DictReader(open(p)))
    return json.load(open(p))


# Columns a CSV may use as its row key, in the order tried. Until 3 Nov this was
# hardcoded to "gene", which was true of the only two CSVs anything declared
# from. crosspop_scores.csv keys on "population", so every lookup against it
# returned None and the declaration came back "value not found" even though the
# value was plainly in the row. A resolver that silently matches nothing looks
# exactly like a missing value.
ROW_KEYS = ("gene", "population", "cluster", "name", "id")


def fetch(obj, path, fn):
    """Walk a key path. For a CSV the first path element is the row key."""
    if fn.endswith(".csv"):
        key, col = path[0], path[1]
        row = None
        for kc in ROW_KEYS:
            row = next((r for r in obj if r.get(kc) == key), None)
            if row is not None:
                break
        if row is None:
            return None
        if col not in row:
            raise KeyError(f"{fn}: row {key!r} has no column {col!r}; "
                           f"columns are {sorted(row)}")
        return float(row[col])
    cur = obj
    for k in path:
        if isinstance(cur, list):
            cur = cur[int(k)]
        else:
            cur = cur[k]
    return cur


def check():
    ok, bad, missing = [], [], []
    cache = {}
    seen_paths = {}
    for entry in MANIFEST:
        label, quoted, fn, path, proc = entry[:5]
        locs = entry[5] if len(entry) > 5 else []
        key = (fn, ".".join(map(str, path)))
        if key in seen_paths:
            print(f"  DUPLICATE PATH  {key[0]}:{key[1]}")
            print(f"      declared as '{seen_paths[key]}' and again as '{label}'")
            print("      merge them into one declaration with both locations\n")
        else:
            seen_paths[key] = label
        if fn not in cache:
            cache[fn] = load(fn)
        obj = cache[fn]
        if obj is None:
            missing.append((label, fn, ".".join(map(str, path)), "file not found"))
            continue
        try:
            got = fetch(obj, path, fn)
        except (KeyError, TypeError, IndexError, ValueError) as e:
            missing.append((label, fn, ".".join(map(str, path)), f"key error: {e}"))
            continue
        if got is None:
            missing.append((label, fn, ".".join(map(str, path)), "value not found"))
            continue
        try:
            agree = abs(float(got) - float(quoted)) <= max(TOL, abs(quoted) * 0.004)
        except (TypeError, ValueError):
            agree = str(got) == str(quoted)
        (ok if agree else bad).append((label, quoted, got, fn,
                                       ".".join(map(str, path)), proc))

    print(f"{len(MANIFEST)} declarations: {len(ok)} agree, {len(bad)} differ, "
          f"{len(missing)} unresolvable\n")
    if bad:
        print("DIFFER (the manuscript and the results no longer agree):")
        for label, q, g, fn, path, proc in bad:
            print(f"  {label:<40} manuscript {q}   stored {g}")
            print(f"      {fn}:{path}")
            print(f"      procedure: {proc}\n")
    if missing:
        print("UNRESOLVABLE (key path gone, or file not produced):")
        for label, fn, path, why in missing:
            print(f"  {label:<40} {fn}:{path}   {why}")
    if not bad and not missing:
        print("Every quoted figure matches the object the manifest declares for it.")

    json.dump({"n_declared": len(MANIFEST), "n_agree": len(ok),
               "n_differ": len(bad), "n_unresolvable": len(missing),
               "differ": [{"label": b[0], "manuscript": b[1], "stored": b[2],
                           "file": b[3], "key": b[4], "procedure": b[5]} for b in bad],
               "unresolvable": [{"label": m[0], "file": m[1], "key": m[2],
                                 "why": m[3]} for m in missing]},
              open(os.path.join(R, "manifest_check.json"), "w"), indent=2)
    return 1 if (bad or missing) else 0


def emit():
    """Print the manifest as a table for the supplement."""
    print("| quantity | value | source | key | quoted in | procedure |")
    print("|---|---|---|---|---|---|")
    for entry in MANIFEST:
        label, quoted, fn, path, proc = entry[:5]
        locs = ", ".join(entry[5]) if len(entry) > 5 and entry[5] else "-"
        print(f"| {label} | {quoted} | `{fn}` | `{'.'.join(map(str, path))}` "
              f"| {locs} | {proc} |")


def orphans():
    """Stored scalars that no manifest line claims. Candidates for the next Results draft."""
    claimed = {(e[2], ".".join(map(str, e[3]))) for e in MANIFEST}
    for fn in sorted({e[2] for e in MANIFEST if e[2].endswith('.json')}):
        obj = load(fn)
        if obj is None:
            continue
        found = []

        def walk(o, path):
            if isinstance(o, dict):
                for k, v in o.items():
                    walk(v, path + [k])
            elif isinstance(o, (int, float)) and not isinstance(o, bool):
                key = ".".join(path)
                if (fn, key) not in claimed:
                    found.append((key, o))

        walk(obj, [])
        if found:
            print(f"\n{fn}: {len(found)} stored values not claimed by the manifest")
            for k, v in found[:12]:
                print(f"    {k:<58} {v}")
            if len(found) > 12:
                print(f"    ... and {len(found) - 12} more")



def scope():
    """Which producers of quoted numbers does the verification reach?"""
    import glob, re as _re

    declared = sorted({e[2] for e in MANIFEST})
    # 34_manifest.py names every declared source file, and the loose second
    # pattern below would make it self-match as their writer. Exclude it.
    scripts = [sc for sc in sorted(glob.glob("[0-9][0-9]_*.py"))
               if os.path.basename(sc) != os.path.basename(__file__)]

    # which script writes each declared source file
    writer = {}
    for f in declared:
        stem = os.path.basename(f)
        for sc in scripts:
            body = open(sc).read()
            # Match a WRITE, not a mention: readers name the same file.
            # A write may name the file in a loop header rather than at the write
            # itself, so also accept the filename appearing in a tuple of outputs
            # alongside a DictWriter or open() call anywhere in the script.
            w = _re.search(
                rf'(json\.dump|to_csv|DictWriter|savez|open)\s*\([^)]*{_re.escape(stem)}'
                rf'|{_re.escape(stem)}[^)\n]*["\']w["\']', body)
            if not w and _re.search(rf'["\']{_re.escape(stem)}["\']', body) and (
                    "DictWriter" in body or "json.dump" in body or "savez" in body):
                w = True
            if w:
                writer.setdefault(f, sc)
                break

    try:
        verify = open("run_verify.slurm").read()
    except FileNotFoundError:
        print("run_verify.slurm not found; cannot report scope")
        return

    def expensive(sc):
        b = open(sc).read()
        return any(k in b for k in ("pass_rows", "pass_totals",
                                   "extract_lepr_block", "matrix.mtx"))

    covered, gap = [], []
    for f in declared:
        sc = writer.get(f)
        if sc is None:
            gap.append((f, "no producer found", False))
        elif sc in verify:
            covered.append((f, sc, expensive(sc)))
        else:
            gap.append((f, sc, expensive(sc)))

    n_decl = {f: sum(1 for e in MANIFEST if e[2] == f) for f in declared}

    print("Verification scope, derived from the manifest\n")
    print(f"  {len(declared)} source files hold the {len(MANIFEST)} quoted numbers\n")
    if covered:
        print("REACHED by run_verify.slurm:")
        for f, sc, exp in covered:
            print(f"  {f:<32} {sc:<28} {n_decl[f]:>2} declarations")
    if gap:
        print("\nNOT REACHED (a quoted number nothing verifies):")
        for f, sc, exp in gap:
            tag = "  [matrix pass]" if exp else ""
            print(f"  {f:<32} {sc:<28} {n_decl[f]:>2} declarations{tag}")
        cheap = [g for g in gap if not g[2] and g[1] != "no producer found"]
        dear = [g for g in gap if g[2]]
        print()
        if cheap:
            print(f"  cheap to add: {', '.join(sorted({g[1] for g in cheap}))}")
        if dear:
            print(f"  turns verification into a partial rerun: "
                  f"{', '.join(sorted({g[1] for g in dear}))}")
            print("  measured cost: 02 through 07 together 24 min cold, 08 6 min,\n"
                  "  18 2.5 min, so one full path is about 36 min cold and 12 warm")
        print("\n  labels affected:")
        for f, sc, exp in gap:
            labs = [e[0] for e in MANIFEST if e[2] == f]
            shown = " | ".join(labs[:4]) + (f" | and {len(labs) - 4} more"
                                            if len(labs) > 4 else "")
            print(f"    {f}: {shown}")
        print("\n  A declared number with no producer in the verification is a number"
              "\n  the archive reports but never checks. Extend run_verify.slurm or"
              "\n  remove the declaration.")
        print("\nNOTE: this compares declarations against the verification wrapper. A"
              "\nnumber quoted in the manuscript but never declared is invisible to it.")
        return 1
    else:
        print("\nEvery producer of a quoted number is reached.")
        print("\nNOTE: this compares declarations against the verification wrapper. A"
              "\nnumber quoted in the manuscript but never declared is invisible to it.")
    return 0


if __name__ == "__main__":
    if "--scope" in sys.argv:
        sys.exit(scope())
    elif "--emit" in sys.argv:
        emit()
    elif "--orphans" in sys.argv:
        orphans()
    else:
        sys.exit(check())
