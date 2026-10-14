#!/usr/bin/env python3
"""
23_class_nulls.py - per-class chance expectations, per-cluster sampling, leave-one-out.

The per-class cluster variances reported by 22_nested_partition.py are not comparable
without their own nulls, because the chance expectation of this estimator depends on the
number of clusters and the number of nuclei per cluster, and both differ by class. This
computes, for every conserved class:

  the observed cluster-level variance, averaged over the control-gene draws
  the chance expectation and 95th percentile from shuffling cluster labels within that class
  the excess over chance, which is the comparable quantity
  nuclei per cluster, which is the sampling the estimator actually had
  leave-one-cluster-out values, which say whether the result is distributed or carried by one

Run for all six classes rather than only the two of interest, so the leptin-receptor class
is described the same way as every other.

Output: results/class_nulls.json, results/class_membership.csv
"""
import os, csv, json, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "cn"
N_NULL = 300


def one_way(y, g):
    df = pd.DataFrame({"y": np.asarray(y, float), "g": np.asarray(g)})
    n, k = len(df), df.g.nunique()
    if k < 2 or n <= k:
        return None
    gm = df.y.mean()
    ssb = sum(len(x) * (x.y.mean() - gm) ** 2 for _, x in df.groupby("g"))
    ssw = sum(((x.y - x.y.mean()) ** 2).sum() for _, x in df.groupby("g"))
    msb, msw = ssb / (k - 1), ssw / (n - k)
    ni = df.groupby("g").size().to_numpy(float)
    n0 = (n - (ni ** 2).sum() / n) / (k - 1)
    return 100 * max((msb - msw) / n0, 0.0) / (max((msb - msw) / n0, 0.0) + msw)


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["Conserved_class"] = ("" if (not raw or raw.startswith("Doublets"))
                                else raw.split("_")[0])
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "Conserved_class", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    c2c = {}
    for l in open(cg):
        a, b = l.split()
        c2c[int(a)] = b
    cols = sorted(c for c in colsum if c in c2c)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cls = np.array([c2c[c] for c in cols])
    clu = np.array([(meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols])
    ld = np.log10(np.array([colsum[c] for c in cols], float))

    # class-by-cluster membership table, so cross-family membership is visible
    memb = defaultdict(lambda: defaultdict(int))
    for a, b in zip(cls, clu):
        memb[a][b] += 1
    with open(os.path.join(cfg.RESULTS, "class_membership.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["conserved_class", "vmh_cluster", "n_nuclei"])
        for a in sorted(memb):
            for b, n in sorted(memb[a].items(), key=lambda t: -t[1]):
                w.writerow([a, b, n])

    seeds = list(range(cfg.RANDOM_SEED, cfg.RANDOM_SEED + cfg.N_SEEDS))
    seed_sets, union = [], set()
    for s in seeds:
        pr, cr = V.choose_controls(name2row, rowsum, len(cols), s)
        seed_sets.append((pr, cr))
        union |= {r for rs_ in list(pr.values()) + list(cr.values()) for r in rs_}
    ent = V.pass_rows(PREFIX, union, TAG, "union")
    expr = defaultdict(dict)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in c2c:
            expr[c][int(r)] = float(v)

    scores = []
    for pr, cr in seed_sets:
        y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        scores.append(y - np.polyval(np.polyfit(ld, y, 1), ld))

    rng = np.random.default_rng(cfg.RANDOM_SEED)
    out = {"n_nuclei": len(cols), "n_seeds": len(seeds), "n_null_shuffles": N_NULL,
           "note": ("chance expectation is computed separately for each class, because it "
                    "depends on the number of clusters and on nuclei per cluster, both of "
                    "which differ by class; the comparable quantity is excess over chance")}
    per = {}
    for cl in sorted(set(cls.tolist())):
        m = cls == cl
        labs = clu[m]
        obs = [one_way(y[m], labs) for y in scores]
        obs = [v for v in obs if v is not None]
        y0 = scores[0][m]
        nulls, lab = [], labs.copy()
        for _ in range(N_NULL):
            rng.shuffle(lab)
            v = one_way(y0, lab)
            if v is not None:
                nulls.append(v)
        sizes = pd.Series(labs).value_counts()
        loo = {}
        for drop in sorted(set(labs)):
            keep = labs != drop
            v = one_way(y0[keep], labs[keep])
            loo[drop] = round(v, 3) if v is not None else None
        per[cl] = {
            "n_nuclei": int(m.sum()), "n_clusters": int(len(set(labs))),
            "nuclei_per_cluster_mean": round(float(sizes.mean()), 1),
            "nuclei_per_cluster_median": int(sizes.median()),
            "observed_pct": {"mean": round(float(np.mean(obs)), 4),
                             "min": round(float(np.min(obs)), 4),
                             "max": round(float(np.max(obs)), 4),
                             "n_draws_at_zero": int(sum(1 for v in obs if v == 0.0))},
            "null_pct": {"mean": round(float(np.mean(nulls)), 4),
                         "p95": round(float(np.percentile(nulls, 95)), 4)},
            "excess_over_null": round(float(np.mean(obs) - np.mean(nulls)), 4),
            "leave_one_cluster_out": loo,
            "loo_min": round(min(v for v in loo.values() if v is not None), 3),
            "loo_max": round(max(v for v in loo.values() if v is not None), 3)}
        print(f"   {cl:<7} obs {per[cl]['observed_pct']['mean']:7.3f}%  "
              f"null {per[cl]['null_pct']['mean']:6.3f}%  "
              f"excess {per[cl]['excess_over_null']:7.3f}%  "
              f"n/clu {per[cl]['nuclei_per_cluster_mean']:6.1f}  "
              f"LOO {per[cl]['loo_min']:.2f} to {per[cl]['loo_max']:.2f}", flush=True)
    out["by_class"] = per
    json.dump(out, open(os.path.join(cfg.RESULTS, "class_nulls.json"), "w"), indent=2)
    print("\nWritten to results/class_nulls.json and results/class_membership.csv")


if __name__ == "__main__":
    main()
