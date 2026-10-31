#!/usr/bin/env python3
"""
55_scan_max_null.py - what does the largest of 12,182 partial correlations look like
when nothing is there?

WHY. Section 2.5 now argues from where the maximum sits. It reports a largest partial
correlation of 0.137 (0.143 with the upper detection bound removed) across 12,182 genes,
and compares that against the ceiling the axis's reliability permits, 0.511. The upper half
of that comparison has a reference; the lower half does not.

The 0.036 figure the section gives is a PER-TEST significance threshold. The quantity being
interpreted is a MAXIMUM over 12,182 tests, and the maximum of many null draws is far above
any per-test threshold by construction. So a reader asking whether 0.143 is simply the
largest of twelve thousand noise correlations has no answer in the text.

It is also the one claim in this manuscript reported against a threshold rather than as a
percentile of its own null. The per-class excesses, the covariation contrast, the
random-panel controls and the Trpc5 comparison are all percentiles; Methods says values are
reported that way rather than against a chosen threshold.

WHAT THIS DOES. Shuffles the per-nucleus axis across nuclei, recomputes every gene's partial
correlation exactly as 08_gene_scan.py does, and records the maximum absolute value. Repeats.
That gives the distribution of the maximum under no association, and 0.137 and 0.143 become
percentiles of it.

MATCHING THE OBSERVED COMPUTATION EXACTLY, from scan B of 08_gene_scan.py:
  - genes detected in 5% to 95% of Lepr-cluster nuclei (MIN_DETECT_FRAC, MAX_DETECT_FRAC)
  - resid_rank: rank the values, then remove log10 depth linearly from the ranks
  - partial rho: Spearman between the residualized ranks of the gene and of the axis
  Any departure from this measures a different quantity and the comparison is void.

WHAT MUST BE TRUE IF THIS IS RIGHT, stated before running:
  - the shuffled maximum must be WELL ABOVE the per-test threshold of 0.036, because it is
    a maximum over 12,182 draws. If it comes back near 0.036 the shuffle is not destroying
    what it should, or the recomputation is not matching the observed one
  - the shuffled MEDIAN gene correlation must be near zero
  - and the observed scan rerun through this code must reproduce 0.137, which is the
    positive control: if it does not, the null is not describing the same statistic

BOTH OUTCOMES ARE REPORTABLE. A shuffled maximum well below 0.137 sharpens the section: the
signal is real but small and the largest of it is ambient RNA. A shuffled maximum near 0.137
strengthens the claim differently: nothing in the transcriptome exceeds what shuffling
produces.

Output: results/scan_max_null.json
"""
import json
import os
import warnings
from collections import defaultdict

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "smn"
MIN_DETECT_FRAC = 0.05
MAX_DETECT_FRAC = 0.95
N_PERM = 50


def main():
    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    hd = np.asarray(z["HminusD"], float)
    cols = [int(c) for c in z["cols"]]
    cidx = {c: i for i, c in enumerate(cols)}
    n = len(cols)
    print(f">> {n} Lepr-cluster nuclei", flush=True)

    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_grp", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    depth = np.array([float(colsum[c]) for c in cols])
    logdepth = np.log10(np.where(depth <= 0, 1, depth))

    # every row with any count; the detection filter is applied after extraction, as in 08
    need = {int(r) for r in rowsum if rowsum[r] > 0}
    print(f">> extracting {len(need)} rows", flush=True)
    ent = V.pass_rows(PREFIX, need, TAG, "smn")
    expr = defaultdict(lambda: np.zeros(n))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            expr[int(r)][cidx[c]] = float(v)

    # Scan B excludes the panel genes AND their matched control genes, which define
    # the axis; scan A does not. That is the whole of the 12,361 against 12,182
    # difference the manuscript carries between 2.6 and 2.5. The first version of this
    # script omitted the exclusion and returned a maximum of 0.2735 against the
    # reported 0.137: the largest correlate was a control gene, which is a gene chosen
    # for resembling a panel gene in detection and then subtracted from the axis.
    panel_rows, ctrl_rows = V.choose_controls(
        name2row, V.read_two_col(V.dpath("rowsum_vmh.txt")), n, cfg.RANDOM_SEED)
    excluded = {r for rs in list(panel_rows.values()) + list(ctrl_rows.values())
                for r in rs}
    excluded |= {name2row[g] for p_ in cfg.PANELS for g in cfg.PANELS[p_]
                 if g in name2row}
    lo, hi = MIN_DETECT_FRAC * n, MAX_DETECT_FRAC * n
    # ORDER MATTERS, and the first version got it backwards. Scan B builds keep from
    # the detection bounds alone (line 119 of 08_gene_scan.py) and then skips excluded
    # rows inside the loop, so its count is keep minus the excluded rows that fall
    # WITHIN the bounds. Excluding first and filtering second removes excluded rows
    # that would never have passed anyway, and gave 12,149 against the reported 12,182:
    # 33 panel or control genes lying outside 5% to 95% detection.
    keep_all = [r for r in expr if lo <= (expr[r] > 0).sum() <= hi]
    keep = [r for r in keep_all if r not in excluded]
    print(f">> {len(keep_all)} genes pass detection; {len(keep_all) - len(keep)} of "
          f"them are panel or control genes, leaving {len(keep)}")
    print(f">> {len(excluded)} panel and control rows excluded, as scan B does")
    print(f">> section 2.5 reports 12,182", flush=True)

    # exactly scan B's construction
    def resid_rank(x):
        rx = stats.rankdata(x)
        b = np.polyfit(logdepth, rx, 1)
        return rx - (b[0] * logdepth + b[1])

    mat = np.vstack([resid_rank(expr[r]) for r in keep])

    def max_abs_rho(axis_vals):
        yr = resid_rank(axis_vals)
        yc = yr - yr.mean()
        den_y = np.sqrt((yc ** 2).sum())
        xc = mat - mat.mean(axis=1, keepdims=True)
        den_x = np.sqrt((xc ** 2).sum(axis=1))
        with np.errstate(invalid="ignore", divide="ignore"):
            rho = (xc @ yc) / (den_x * den_y)
        rho = rho[np.isfinite(rho)]
        return float(np.max(np.abs(rho))), float(np.median(rho))

    obs_max, obs_med = max_abs_rho(hd)
    print(f">> observed: max |partial rho| {obs_max:.4f}, median {obs_med:+.4f}",
          flush=True)
    print(f"   section 2.5 reports 0.137 within the detection bounds", flush=True)

    rng = np.random.default_rng(cfg.RANDOM_SEED)
    nulls, meds = [], []
    for i in range(N_PERM):
        m_, md = max_abs_rho(rng.permutation(hd))
        nulls.append(m_)
        meds.append(md)
        if (i + 1) % 10 == 0:
            print(f"   {i + 1}/{N_PERM} permutations, running max "
                  f"{np.mean(nulls):.4f}", flush=True)
    nulls = np.array(nulls)

    out = {"n_nuclei": n, "n_genes": len(keep), "n_permutations": N_PERM,
           "detection_bounds": [MIN_DETECT_FRAC, MAX_DETECT_FRAC],
           "construction": ("rank, remove log10 depth linearly from the ranks, then "
                            "Spearman between residualized ranks; identical to scan B "
                            "of 08_gene_scan.py"),
           "observed_max_abs_rho": round(obs_max, 4),
           "observed_median_rho": round(obs_med, 4),
           "null_max_mean": round(float(nulls.mean()), 4),
           "null_max_p5": round(float(np.percentile(nulls, 5)), 4),
           "null_max_p95": round(float(np.percentile(nulls, 95)), 4),
           "null_max_min": round(float(nulls.min()), 4),
           "null_max_max": round(float(nulls.max()), 4),
           "null_median_gene_rho_mean": round(float(np.mean(meds)), 5),
           "percentile_of_observed": round(float((nulls < obs_max).mean() * 100), 1),
           "per_test_threshold": round(float(1.96 / np.sqrt(n - 3)), 4)}

    print(f"\n{'=' * 68}")
    print(f"  observed maximum            {out['observed_max_abs_rho']}")
    print(f"  null maximum, mean          {out['null_max_mean']}  "
          f"(5th-95th {out['null_max_p5']} to {out['null_max_p95']})")
    print(f"  per-test threshold          {out['per_test_threshold']}  "
          f"(what 2.5 currently quotes, for a SINGLE test)")
    print(f"  observed sits at the        {out['percentile_of_observed']}th percentile "
          f"of the null maximum")

    # The first version demanded 3x the per-test threshold, a multiplier chosen
    # with no basis. A maximum over ~12,000 correlated draws runs about 2x.
    ctrl_a = out["null_max_mean"] > 1.5 * out["per_test_threshold"]
    ctrl_b = abs(out["null_median_gene_rho_mean"]) < 0.01
    out["controls_pass"] = bool(ctrl_a and ctrl_b)
    print(f"\n  null maximum well above the per-test threshold: {ctrl_a}")
    print(f"  shuffled median gene correlation near zero:     {ctrl_b} "
          f"({out['null_median_gene_rho_mean']})")
    if not out["controls_pass"]:
        print("\n  A control failed. The shuffle or the recomputation is not matching")
        print("  the observed statistic, and nothing above is readable.")
    elif out["percentile_of_observed"] > 95:
        print("\n  The observed maximum exceeds what shuffling produces. 2.5 can say")
        print("  the signal is real but small, and that the largest of it is ambient")
        print("  RNA, which is sharper than the current wording.")
    else:
        print("\n  The observed maximum is within the range shuffling produces. 2.5's")
        print("  claim becomes stronger and simpler: nothing in the transcriptome")
        print("  exceeds what a shuffled axis yields.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "scan_max_null.json"), "w"), indent=2)
    print("\nWritten to results/scan_max_null.json")


if __name__ == "__main__":
    main()
