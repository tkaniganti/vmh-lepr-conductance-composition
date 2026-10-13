#!/usr/bin/env python3
"""
03_trpc5_depth.py - is the Trpc5 detection gradient across LepR clusters real, or depth?

The manuscript argues that single-nucleus per-gene channel calls are unreliable because
detection tracks capture depth, then keeps one per-gene detection claim (Trpc5, 30% to 68%
across LepR clusters) as its strongest single-gene result. Detection rate is exactly the
statistic most sensitive to per-cluster differences in library size, so the claim has to
survive a depth control before it can be stated.

Three tests:
  1. Logistic regression  detected ~ C(cluster) + log10(total UMI), with a likelihood-ratio
     test of the cluster term against the depth-only model.
  2. Downsampling every nucleus to a common depth, then re-testing detection by chi-square.
  3. An empirical null: the same cluster LR statistic computed for ~200 expression-matched
     control genes, so the Trpc5 statistic is read against genes with no expected biology.

Output: results/trpc5_depth.json, results/trpc5_depth_null.csv
"""
import os, json, csv, math
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2 as chi2dist, chi2_contingency
import config as C
import vmhlib as V

PREFIX, TAG = "GSE172204", "vmh"
N_NULL_GENES = 200
DOWNSAMPLE_SEED = 7


def lr_cluster_effect(detected, cluster, logdepth):
    """LR statistic and df for the cluster term, over and above depth."""
    Xd = sm.add_constant(logdepth.reshape(-1, 1), has_constant="add")
    m0 = sm.Logit(detected, Xd).fit(disp=0, maxiter=200)
    D = pd.get_dummies(cluster, drop_first=True).to_numpy(dtype=float)
    X1 = np.column_stack([Xd, D])
    m1 = sm.Logit(detected, X1).fit(disp=0, maxiter=200)
    lr = 2 * (m1.llf - m0.llf)
    df = D.shape[1]
    return float(lr), int(df), float(chi2dist.sf(lr, df)), m1


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    col2group, _, _ = V.write_col2group(PREFIX, meta, bc2col, "VMH_cluster", TAG)
    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = g

    lepr_cols = sorted([c for c in colsum if col2cl[c].startswith("Lepr")])
    cluster = np.array([col2cl[c] for c in lepr_cols])
    depth = np.array([colsum[c] for c in lepr_cols], dtype=float)
    logdepth = np.log10(depth)
    print(f">> {len(lepr_cols)} LepR nuclei; median depth {np.median(depth):.0f}")
    bycl = pd.DataFrame({"cluster": cluster, "depth": depth}).groupby("cluster")["depth"]
    print(bycl.median().round(0).to_string())

    # ------------------------------------------------------------ gene of interest
    t5 = name2row["Trpc5"]
    ent = V.pass_rows(PREFIX, {t5}, TAG, "trpc5")
    counts = defaultdict(float)
    for l in open(ent):
        r, c, v = l.split()
        counts[int(c)] = float(v)
    det = np.array([1.0 if counts.get(c, 0) > 0 else 0.0 for c in lepr_cols])

    lr, df, p, m1 = lr_cluster_effect(det, cluster, logdepth)
    depth_coef = float(m1.params[1])
    print(f">> Trpc5: LR(cluster | depth) = {lr:.1f}, df = {df}, p = {p:.3g}; "
          f"depth coefficient = {depth_coef:.2f}")

    # ------------------------------------------------------------ 2. downsampling
    rng = np.random.default_rng(DOWNSAMPLE_SEED)
    target = int(np.percentile(depth, 10))     # common depth, 10th percentile
    ds_det = []
    for c, d in zip(lepr_cols, depth):
        k = counts.get(c, 0.0)
        if k <= 0:
            ds_det.append(0.0)
            continue
        # probability of retaining at least one of k molecules when thinning to target/d
        keep = min(1.0, target / d)
        ds_det.append(1.0 if rng.binomial(int(k), keep) > 0 else 0.0)
    ds_det = np.array(ds_det)
    tbl = []
    for k in sorted(set(cluster)):
        m = cluster == k
        tbl.append([int(ds_det[m].sum()), int(m.sum() - ds_det[m].sum())])
    ds_chi2, ds_p, ds_dof, _ = chi2_contingency(np.array(tbl))

    raw_pct = {k: round(100 * float(det[cluster == k].mean()), 1) for k in sorted(set(cluster))}
    ds_pct = {k: round(100 * float(ds_det[cluster == k].mean()), 1) for k in sorted(set(cluster))}
    print(f">> downsampled to {target} UMI: chi2 = {ds_chi2:.1f}, p = {ds_p:.3g}")

    # ------------------------------------------------------------ 3. empirical null
    mean_expr = {r: rowsum.get(r, 0.0) / len(colsum) for r in rowsum}
    t5_mean = mean_expr.get(t5, 0.0)
    cand = sorted((r for r in mean_expr if r != t5 and mean_expr[r] > 0),
                  key=lambda r: abs(math.log1p(mean_expr[r]) - math.log1p(t5_mean)))
    panel_rows = {name2row[g] for p in ("HYPER", "DEPOL", "MARKERS")
                  for g in C.PANELS[p] if g in name2row}
    null_rows = [r for r in cand if r not in panel_rows][:N_NULL_GENES]
    ent_n = V.pass_rows(PREFIX, set(null_rows), TAG, "null_depth")
    ncount = defaultdict(dict)
    for l in open(ent_n):
        r, c, v = l.split()
        ncount[int(r)][int(c)] = float(v)

    sym = {v: k for k, v in name2row.items()}
    null_rec = []
    for r in null_rows:
        d = np.array([1.0 if ncount[r].get(c, 0) > 0 else 0.0 for c in lepr_cols])
        if d.sum() < 20 or d.sum() > len(d) - 20:
            continue
        try:
            l_, df_, p_, _ = lr_cluster_effect(d, cluster, logdepth)
        except Exception:
            continue
        null_rec.append({"gene": sym.get(r, str(r)), "lr": round(l_, 2), "p": p_})
    null_lr = np.array([x["lr"] for x in null_rec]) if null_rec else np.array([np.nan])
    pct_rank = float((null_lr < lr).mean() * 100) if null_rec else float("nan")
    print(f">> empirical null: {len(null_rec)} depth-matched genes, "
          f"median LR {np.nanmedian(null_lr):.1f}; Trpc5 exceeds {pct_rank:.1f}% of them")

    out = {
        "n_lepr_nuclei": len(lepr_cols),
        "median_depth_by_cluster": {k: float(v) for k, v in bycl.median().items()},
        "logistic_cluster_given_depth": {"LR": round(lr, 2), "df": df, "p": p,
                                         "depth_coefficient": round(depth_coef, 3)},
        "downsampling": {"target_umi": target, "chi2": round(float(ds_chi2), 1),
                         "p": float(ds_p), "dof": int(ds_dof),
                         "pct_detected_raw": raw_pct, "pct_detected_downsampled": ds_pct},
        "empirical_null": {"n_genes": len(null_rec),
                           "median_LR": float(np.nanmedian(null_lr)),
                           "p90_LR": float(np.nanpercentile(null_lr, 90)),
                           "trpc5_percentile": pct_rank},
        "interpretation": ("Trpc5 cluster structure survives depth control if the LR p is "
                           "small, the downsampled chi-square remains significant, and the "
                           "Trpc5 LR sits far above the depth-matched null distribution. "
                           "If it does not, the gradient should be dropped from the paper."),
    }
    json.dump(out, open(os.path.join(C.RESULTS, "trpc5_depth.json"), "w"), indent=2)
    with open(os.path.join(C.RESULTS, "trpc5_depth_null.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "lr", "p"])
        w.writeheader(); w.writerows(null_rec)
    print("\nWritten to", C.RESULTS)


if __name__ == "__main__":
    main()
