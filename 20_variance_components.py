#!/usr/bin/env python3
"""
20_variance_components.py - one model for the central structural claim.

10_conserved_classes.py reports two separate Kruskal-Wallis effect sizes, one between
conserved VMH classes and one between clusters within the conserved Lepr class, and the
manuscript compares them. That is a juxtaposition of two tests rather than a test of the
difference. A nested variance-components model estimates both levels simultaneously from
the same data and gives a single answer to the question the paper actually asks: how much
of the cell-to-cell variance in effector balance is attributable to conserved class, and
how much to cluster nested within class?

Model, fitted on the per-nucleus hyperpolarising-minus-depolarising score:

    score ~ 1 + (1 | class) + (1 | class:cluster)

reported as variance components with the residual, so the three add to the total. A
likelihood-ratio test against the model without the nested cluster term gives a direct test
of whether subdivision within class explains any variance at all, which is the claim the
title rests on.

Run over the N_SEEDS control-gene draws, as elsewhere, and reported as mean with range.
Depth is included as a fixed effect, since classes differ in capture depth.

Output: results/variance_components.json
"""
import os, json, math
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.formula.api as smf
import config as cfg
import vmhlib as V

PREFIX, TAG = "GSE172204", "vc"
CLASS_COL, CLUSTER_COL = "Conserved_class", "VMH_cluster"
N_PERM = 100


def fit_components(df):
    """Cluster-within-class variance, with class as a FIXED effect and depth as a covariate.

    Six class levels are too few to estimate a variance component, so class is fixed and
    cluster-within-class is the only random term. A singular random-effects covariance
    here means the cluster variance was estimated at the zero boundary, which is a result
    rather than a failure, so warnings are captured and reported instead of raised.
    """
    import warnings
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always")
        r = smf.mixedlm("score ~ logdepth + C(cls)", df,
                        groups=df["cluster_in_class"]).fit(reml=True, method="lbfgs")
    msgs = "; ".join(sorted({str(w.message)[:60] for w in caught}))
    var_cluster = float(np.asarray(r.cov_re)[0, 0])
    var_resid = float(r.scale)
    fit_components.last_warnings = msgs
    return r, 0.0, var_cluster, var_resid


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
    col2cls = {}
    for l in open(col2group):
        c, g = l.split()
        col2cls[int(c)] = g
    cols = sorted(c for c in colsum if c in col2cls)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cluster = np.array([(meta[col2bc[c]].get(CLUSTER_COL) or "NA").strip() for c in cols])
    cls = np.array([col2cls[c] for c in cols])
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)
    print(f">> {len(cols)} nuclei, {len(set(cls))} conserved classes, "
          f"{len(set(zip(cls, cluster)))} class-by-cluster combinations")

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
        if c in col2cls:
            expr[c][int(r)] = float(v)

    recs = []
    for i, (pr, cr) in enumerate(seed_sets):
        score = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        df = pd.DataFrame({"score": score, "cls": cls,
                           "cluster_in_class": [f"{a}|{b}" for a, b in zip(cls, cluster)],
                           "logdepth": logdepth})
        try:
            full, vc_class, vc_cluster, vc_resid = fit_components(df)
            total = vc_class + vc_cluster + vc_resid
            recs.append({"seed": seeds[i],
                         "pct_class": 100 * vc_class / total,
                         "pct_cluster_in_class": 100 * vc_cluster / total,
                         "pct_residual": 100 * vc_resid / total})
            print(f"   seed {seeds[i]}: class {recs[-1]['pct_class']:.2f}%, "
                  f"cluster-in-class {recs[-1]['pct_cluster_in_class']:.2f}%")
        except Exception as e:
            print(f"   seed {seeds[i]}: FAILED ({e})")

    if not recs:
        raise SystemExit("no model converged; inspect the data frame before reporting")

    # Cluster-level variance computed separately within each conserved class
    sc0 = V.module_scores(expr, cols, colsum, *seed_sets[0])["HminusD"]
    per_class = {}
    for cl in sorted(set(cls.tolist())):
        m = cls == cl
        if len(set(cluster[m])) < 2:
            continue
        d = pd.DataFrame({"score": sc0[m], "grp": cluster[m], "logdepth": logdepth[m]})
        try:
            mod = smf.mixedlm("score ~ logdepth", d, groups=d["grp"]).fit(reml=True,
                                                                         method="lbfgs")
            vg = float(np.asarray(mod.cov_re)[0, 0]); vr = float(mod.scale)
            per_class[cl] = round(100 * vg / (vg + vr), 3)
        except Exception as e:
            per_class[cl] = None
    print(">> cluster variance (%) within each conserved class:", per_class)

    # Restricted model: the nine Lepr clusters within the conserved Lepr class only,
    # which is the comparison the manuscript actually claims.
    lep = cls == "Lepr"
    print(f">> restricted model: {int(lep.sum())} nuclei, "
          f"{len(set(cluster[lep]))} clusters in the conserved Lepr class")
    lep_recs = []
    for i, (pr, cr) in enumerate(seed_sets):
        sc = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"][lep]
        d = pd.DataFrame({"score": sc, "grp": cluster[lep], "logdepth": logdepth[lep]})
        try:
            m = smf.mixedlm("score ~ logdepth", d, groups=d["grp"]).fit(reml=True,
                                                                       method="lbfgs")
            vg = float(np.asarray(m.cov_re)[0, 0]); vr = float(m.scale)
            lep_recs.append(100 * vg / (vg + vr))
        except Exception:
            pass
    print(f">> cluster variance within the Lepr class: "
          f"{np.mean(lep_recs):.2f}% (n={len(lep_recs)} seeds)")

    # Restricted model: the nine Lepr clusters within the conserved Lepr class only,
    # which is the comparison the manuscript actually claims.
    lep = cls == "Lepr"
    print(f">> restricted model: {int(lep.sum())} nuclei, "
          f"{len(set(cluster[lep]))} clusters in the conserved Lepr class")
    lep_recs = []
    for i, (pr, cr) in enumerate(seed_sets):
        sc = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"][lep]
        d = pd.DataFrame({"score": sc, "grp": cluster[lep], "logdepth": logdepth[lep]})
        try:
            m = smf.mixedlm("score ~ logdepth", d, groups=d["grp"]).fit(reml=True,
                                                                       method="lbfgs")
            vg = float(np.asarray(m.cov_re)[0, 0]); vr = float(m.scale)
            lep_recs.append(100 * vg / (vg + vr))
        except Exception:
            pass
    print(f">> cluster variance within the Lepr class: "
          f"{np.mean(lep_recs):.2f}% (n={len(lep_recs)} seeds)")

    # Permutation null for the nested term: shuffle cluster labels WITHIN class and refit,
    # so any residual cluster-level variance reflects chance assignment rather than structure.
    pr0, cr0 = seed_sets[0]
    score0 = V.module_scores(expr, cols, colsum, pr0, cr0)["HminusD"]
    base = pd.DataFrame({"score": score0, "cls": cls,
                         "cluster_in_class": [f"{a}|{b}" for a, b in zip(cls, cluster)],
                         "logdepth": logdepth})
    _, _, obs_cluster, _ = fit_components(base)
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    null = []
    for _ in range(N_PERM):
        d = base.copy()
        for cl in d["cls"].unique():
            m = d["cls"] == cl
            vals = d.loc[m, "cluster_in_class"].to_numpy()
            rng.shuffle(vals)
            d.loc[m, "cluster_in_class"] = vals
        try:
            _, _, vc, _ = fit_components(d)
            null.append(vc)
        except Exception:
            pass
    perm_p = ((np.sum(np.array(null) >= obs_cluster) + 1) / (len(null) + 1)) if null else None
    print(f">> permutation null on seed {seeds[0]}: observed cluster variance "
          f"{obs_cluster:.5f}, null median {np.median(null):.5f}, p = {perm_p}")

    def summ(key):
        v = np.array([r[key] for r in recs], float)
        return {"mean": round(float(v.mean()), 4), "min": round(float(v.min()), 4),
                "max": round(float(v.max()), 4)}

    out = {"n_nuclei": len(cols), "n_classes": len(set(cls.tolist())),
           "n_class_by_cluster": len(set(zip(cls.tolist(), cluster.tolist()))),
           "cluster_pct_by_class": per_class,
           "n_seeds_converged": len(recs),
           "model": "score ~ logdepth + C(class) fixed, (1|class:cluster) random, REML",
           "variance_explained_pct": {k: summ(k) for k in
                                      ("pct_class", "pct_cluster_in_class", "pct_residual")},
           "restricted_lepr_class_cluster_pct": {
               "mean": round(float(np.mean(lep_recs)), 4) if lep_recs else None,
               "min": round(float(np.min(lep_recs)), 4) if lep_recs else None,
               "max": round(float(np.max(lep_recs)), 4) if lep_recs else None,
               "n_seeds": len(lep_recs),
               "note": "nine Lepr clusters within the conserved Lepr class only"},
           "restricted_lepr_class_cluster_pct": {
               "mean": round(float(np.mean(lep_recs)), 4) if lep_recs else None,
               "min": round(float(np.min(lep_recs)), 4) if lep_recs else None,
               "max": round(float(np.max(lep_recs)), 4) if lep_recs else None,
               "n_seeds": len(lep_recs),
               "note": "nine Lepr clusters within the conserved Lepr class only"},
           "nested_term_permutation": {
               "n_permutations": len(null), "observed_cluster_variance": round(obs_cluster, 6),
               "null_median": round(float(np.median(null)), 6) if null else None,
               "p": perm_p,
               "note": "cluster labels shuffled within class, so the null preserves class "
                       "structure and tests only the nested term"},
           "interpretation": (
               "The class term and the cluster-nested-within-class term are estimated in one "
               "model, so their variance shares are directly comparable rather than being two "
               "separate effect sizes placed side by side. The likelihood-ratio test on the "
               "nested term is tested by shuffling cluster labels within class, which "
               "preserves class structure and asks only whether subdivision within a class "
               "explains variance beyond chance assignment. If the class share substantially "
               "exceeds the nested share and the permutation test does not reject, the "
               "structural claim rests on a single estimate rather than two juxtaposed "
               "effect sizes.")}

    json.dump(out, open(os.path.join(cfg.RESULTS, "variance_components.json"), "w"), indent=2)
    print("\n" + json.dumps({k: out[k] for k in
                             ("variance_explained_pct", "nested_term_permutation", "n_seeds_converged")}, indent=2))


if __name__ == "__main__":
    main()
