#!/usr/bin/env python3
"""
21_variance_partition.py - nested variance partition by components of variance.

Replaces the mixed-model approach in 20_variance_components.py, which failed its own
positive control: REML returned an exact zero for the between-class variance on all twenty
control-gene draws, while a direct decomposition of the same data gives 2.6%. Six class
levels are too few for REML to identify a variance component, and the optimiser parks at
the zero boundary. Depth is not the cause: removing it linearly leaves the class effect
slightly larger, and the mixed model returns zero with no depth covariate at all.

What this computes instead is a single estimator applied at both levels, which is what the
comparison in the manuscript requires. For a grouping with k groups and n observations:

    sigma_between = max((MS_between - MS_within) / n_bar, 0)
    share         = sigma_between / (sigma_between + MS_within)

reported alongside eta-squared from the same sums of squares. No optimiser is involved, so
there is no boundary behaviour, and the two levels are directly comparable because they are
computed the same way.

Levels:
  between the six cross-species conserved VMH classes (5,441 nuclei)
  between the VMH clusters inside the conserved Lepr class (11 clusters, 1,857 nuclei)

Depth is removed linearly from the score before partitioning, and the whole thing runs over
the N_SEEDS control-gene draws and is reported as mean with range.

Output: results/variance_partition.json
"""
import os, json, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "vp"
CLASS_COL, CLUSTER_COL = "Conserved_class", "VMH_cluster"
FOCUS_CLASS = "Lepr"


def components(y, grp):
    """Components of variance for one grouping. Returns share, eta-squared, k, n."""
    df = pd.DataFrame({"y": np.asarray(y, float), "g": np.asarray(grp)})
    n, k = len(df), df.g.nunique()
    if k < 2:
        return None
    grand = df.y.mean()
    ssb = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in df.groupby("g"))
    ssw = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in df.groupby("g"))
    msb, msw = ssb / (k - 1), ssw / (n - k)
    sb = max((msb - msw) / (n / k), 0.0)
    return (100 * sb / (sb + msw), ssb / (ssb + ssw), k, n)


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m[CLASS_COL] = "" if (not raw or raw.startswith("Doublets")) else raw.split("_")[0]

    col2group, _, _ = V.write_col2group(PREFIX, meta, bc2col, CLASS_COL, TAG)
    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    c2c = {}
    for l in open(col2group):
        a, b = l.split()
        c2c[int(a)] = b
    cols = sorted(c for c in colsum if c in c2c)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cls = np.array([c2c[c] for c in cols])
    clu = np.array([(meta[col2bc[c]].get(CLUSTER_COL) or "NA").strip() for c in cols])
    logdepth = np.log10(np.array([colsum[c] for c in cols], float))
    focus = cls == FOCUS_CLASS
    print(f">> {len(cols)} nuclei; {len(set(cls))} conserved classes; "
          f"{len(set(clu[focus]))} clusters inside the {FOCUS_CLASS} class "
          f"({int(focus.sum())} nuclei)", flush=True)

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

    rows = {"between_conserved_classes": [], "between_clusters_within_lepr_class": [],
            "between_all_clusters": []}
    for seed, (pr, cr) in zip(seeds, seed_sets):
        y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        y = y - np.polyval(np.polyfit(logdepth, y, 1), logdepth)   # depth removed linearly
        rows["between_conserved_classes"].append(components(y, cls))
        rows["between_clusters_within_lepr_class"].append(components(y[focus], clu[focus]))
        rows["between_all_clusters"].append(
            components(y, np.array([f"{a}|{b}" for a, b in zip(cls, clu)])))
        print(f"   seed {seed}: class {rows['between_conserved_classes'][-1][0]:.3f}%, "
              f"within-{FOCUS_CLASS} {rows['between_clusters_within_lepr_class'][-1][0]:.3f}%",
              flush=True)

    out = {"n_nuclei": len(cols), "n_seeds": len(seeds),
           "estimator": ("components of variance, sigma_between = max((MSB - MSW)/n_bar, 0); "
                         "share = sigma_between / (sigma_between + MSW); depth removed "
                         "linearly from the score before partitioning"),
           "note": ("the mixed-model estimator in 20_variance_components.py failed its "
                    "positive control at six class levels and is superseded by this")}
    for key, recs in rows.items():
        a = np.array([[r[0], r[1]] for r in recs if r], float)
        out[key] = {"k": recs[0][2], "n": recs[0][3],
                    "pct_mean": round(float(a[:, 0].mean()), 4),
                    "pct_min": round(float(a[:, 0].min()), 4),
                    "pct_max": round(float(a[:, 0].max()), 4),
                    "eta2_mean": round(float(a[:, 1].mean()), 4),
                    "eta2_min": round(float(a[:, 1].min()), 4),
                    "eta2_max": round(float(a[:, 1].max()), 4)}

    json.dump(out, open(os.path.join(cfg.RESULTS, "variance_partition.json"), "w"), indent=2)
    print("\n" + json.dumps({k: v for k, v in out.items()
                             if k.startswith("between")}, indent=2))
    print("\nWritten to", os.path.join(cfg.RESULTS, "variance_partition.json"))


if __name__ == "__main__":
    main()
