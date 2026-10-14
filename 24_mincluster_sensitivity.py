#!/usr/bin/env python3
"""
24_mincluster_sensitivity.py - does within-class structure survive dropping tiny clusters?

Median cluster size inside most conserved classes is 4 to 7 nuclei, against means of 46 to
66, so a between-cluster variance component is effectively driven by the one or two large
clusters in each class. This refits the per-class component after excluding clusters below a
minimum size, at several thresholds, and reports what survives.

The question is whether the positive control holds. If the Nfib excess persists when only
its substantial clusters are retained, the control rests on distributed structure. If it
collapses, the control rests on a handful of clusters and the manuscript must say so.

Output: results/mincluster_sensitivity.json
"""
import os, json, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "mc"
THRESHOLDS = [0, 10, 20, 50]
N_NULL = 200


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
    sb = max((msb - msw) / n0, 0.0)
    return 100 * sb / (sb + msw)


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
    out = {"thresholds": THRESHOLDS, "n_seeds": len(seeds), "by_class": {}}
    for cl in sorted(set(cls.tolist())):
        m = cls == cl
        labs = clu[m]
        sizes = pd.Series(labs).value_counts()
        rec = {}
        for thr in THRESHOLDS:
            keep_lab = set(sizes[sizes >= thr].index)
            sel = np.isin(labs, list(keep_lab))
            if len(keep_lab) < 2 or sel.sum() < 30:
                rec[str(thr)] = {"status": "too few clusters or nuclei",
                                 "n_clusters": len(keep_lab), "n_nuclei": int(sel.sum())}
                continue
            obs = [one_way(y[m][sel], labs[sel]) for y in scores]
            obs = [v for v in obs if v is not None]
            y0 = scores[0][m][sel]
            lab2, nulls = labs[sel].copy(), []
            for _ in range(N_NULL):
                rng.shuffle(lab2)
                v = one_way(y0, lab2)
                if v is not None:
                    nulls.append(v)
            rec[str(thr)] = {"n_clusters": len(keep_lab), "n_nuclei": int(sel.sum()),
                             "observed_mean": round(float(np.mean(obs)), 4),
                             "null_mean": round(float(np.mean(nulls)), 4),
                             "excess": round(float(np.mean(obs) - np.mean(nulls)), 4)}
            print(f"   {cl:<7} thr>={thr:<3} k={len(keep_lab):2d} n={int(sel.sum()):5d}  "
                  f"obs {rec[str(thr)]['observed_mean']:7.3f}%  "
                  f"excess {rec[str(thr)]['excess']:7.3f}%", flush=True)
        out["by_class"][cl] = rec
    json.dump(out, open(os.path.join(cfg.RESULTS, "mincluster_sensitivity.json"), "w"),
              indent=2)
    print("\nWritten to results/mincluster_sensitivity.json")


if __name__ == "__main__":
    main()
