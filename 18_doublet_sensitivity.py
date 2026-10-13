#!/usr/bin/env python3
"""
18_doublet_sensitivity.py - do the core results survive removal of author-flagged doublets?

The Affinati metadata assigns 461 of the 2,913 LepR-cluster nuclei (15.8%) to Doublets_1,
Doublets_2 or Doublets_3 in the Conserved_cluster column. Doublets inflate apparent
co-expression of opposing conductances, which is the central observation of this analysis,
so every headline number has to be shown with and without them.

Recomputed both ways, in one pass:
  detection rates and the naive double-positive gate
  per-cell module scores and Kruskal-Wallis epsilon-squared over N_SEEDS control-gene draws
  unsupervised structure on the H-D axis (Gaussian mixture by BIC, Hartigan dip test)
  the depth-controlled Trpc5 test (logistic regression of detection on cluster given log depth)

If the doublet-free numbers match, the full set stays primary and this is a supplementary
table. If they move, the doublet-free set becomes primary and the paper says so.

Output: results/doublet_sensitivity.json
"""
import os, json, math
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import kruskal, chi2 as chi2dist
from sklearn.mixture import GaussianMixture
import config as C
import vmhlib as V

try:
    import diptest
    HAVE_DIP = True
except Exception:
    HAVE_DIP = False

PREFIX = "GSE172204"
MARKERS = ["Nr5a1", "Lepr", "Slc17a6"]


def fit_mixture(x, kmax=4, seed=0):
    x = np.asarray(x, float).reshape(-1, 1)
    out = []
    for k in range(1, kmax + 1):
        gm = GaussianMixture(n_components=k, covariance_type="full",
                             random_state=seed, n_init=5).fit(x)
        out.append((k, float(gm.bic(x))))
    best = min(out, key=lambda r: r[1])
    return {"bic_by_k": {str(k): round(b, 1) for k, b in out}, "best_k": best[0]}


def lr_cluster_given_depth(y, cluster, logdepth):
    Xd = sm.add_constant(logdepth.reshape(-1, 1), has_constant="add")
    m0 = sm.Logit(y, Xd).fit(disp=0, maxiter=200)
    D = pd.get_dummies(cluster, drop_first=True).to_numpy(dtype=float)
    m1 = sm.Logit(y, np.column_stack([Xd, D])).fit(disp=0, maxiter=200)
    lr = 2 * (m1.llf - m0.llf)
    return float(lr), int(D.shape[1]), float(chi2dist.sf(max(lr, 0.0), D.shape[1]))


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)

    # doublet status from the deposited Conserved_cluster annotation
    is_doublet = {}
    for bc, m in meta.items():
        raw = (m.get("Conserved_cluster") or "").strip()
        is_doublet[bc] = raw.startswith("Doublets")

    col2group, _, _ = V.write_col2group(PREFIX, meta, bc2col, "VMH_cluster", "dbl")
    cs, rs = V.pass_totals(PREFIX, col2group, "dbl")
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = g
    cols = sorted(c for c in colsum if c in col2cl)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    dbl = np.array([is_doublet.get(col2bc.get(c, ""), False) for c in cols])
    clu = np.array([col2cl[c] for c in cols])
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)
    lepr_mask = np.array([g.startswith("Lepr") for g in clu])
    print(f">> {len(cols)} VMH nuclei; LepR {int(lepr_mask.sum())}; "
          f"doublet-flagged among LepR {int((lepr_mask & dbl).sum())} "
          f"({100 * (lepr_mask & dbl).sum() / lepr_mask.sum():.1f}%)")

    # ---------------------------------------------------- rows: markers + panels + controls
    seeds = list(range(C.RANDOM_SEED, C.RANDOM_SEED + C.N_SEEDS))
    seed_sets, union = [], set()
    for s in seeds:
        pr, cr = V.choose_controls(name2row, rowsum, len(cols), s)
        seed_sets.append((pr, cr))
        union |= {r for rs_ in list(pr.values()) + list(cr.values()) for r in rs_}
    union |= {name2row[g] for g in MARKERS + ["Trpc5"] if g in name2row}
    print(f">> extracting {len(union)} rows in one pass")
    ent = V.pass_rows(PREFIX, union, "dbl", "union")
    expr = defaultdict(dict)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in col2cl:
            expr[c][int(r)] = float(v)

    scores = [V.module_scores(expr, cols, colsum, pr, cr) for pr, cr in seed_sets]

    t5row = name2row["Trpc5"]
    t5 = np.array([1.0 if expr.get(c, {}).get(t5row, 0.0) > 0 else 0.0 for c in cols])
    mk = {g: np.array([1.0 if expr.get(c, {}).get(name2row[g], 0.0) > 0 else 0.0
                       for c in cols]) for g in MARKERS if g in name2row}

    # ---------------------------------------------------------------- both versions
    out = {"n_vmh": len(cols), "n_lepr": int(lepr_mask.sum()),
           "n_lepr_doublet_flagged": int((lepr_mask & dbl).sum()),
           "pct_lepr_doublet_flagged": round(100 * float((lepr_mask & dbl).mean() /
                                                         lepr_mask.mean()), 1)}

    for label, keep in (("with_doublets", np.ones(len(cols), bool)),
                        ("without_doublets", ~dbl)):
        sel = lepr_mask & keep
        vmh_sel = keep
        res = {"n_lepr": int(sel.sum()), "n_vmh": int(vmh_sel.sum())}

        res["detection_pct_vmh"] = {g: round(100 * float(mk[g][vmh_sel].mean()), 1)
                                    for g in mk}
        res["detection_pct_lepr"] = {g: round(100 * float(mk[g][sel].mean()), 1) for g in mk}
        if "Nr5a1" in mk and "Lepr" in mk:
            both = (mk["Nr5a1"] > 0) & (mk["Lepr"] > 0)
            res["naive_gate_lepr"] = int((both & sel).sum())
            res["naive_gate_lepr_pct"] = round(100 * float((both & sel).sum() / sel.sum()), 1)

        lab_sel = clu[sel]
        groups_present = sorted(set(lab_sel.tolist()))
        kw = defaultdict(list)
        for s in scores:
            for lab in ("HYPER", "DEPOL", "HminusD"):
                arr = s[lab][sel]
                g = [arr[lab_sel == k] for k in groups_present if (lab_sel == k).sum() > 1]
                H, _ = kruskal(*g)
                n = int(sum(len(x) for x in g))
                kw[lab].append(V.eps2_from_kruskal(float(H), len(g), n))
        res["eps2"] = {lab: {"mean": round(float(np.mean(v)), 4),
                             "min": round(float(np.min(v)), 4),
                             "max": round(float(np.max(v)), 4)}
                       for lab, v in kw.items()}

        hd = scores[0]["HminusD"][sel]
        res["mixture"] = fit_mixture(hd)
        if HAVE_DIP:
            d, p = diptest.diptest(hd)
            res["dip"] = {"D": round(float(d), 5), "p": round(float(p), 4)}

        try:
            lr, df, p = lr_cluster_given_depth(t5[sel], lab_sel, logdepth[sel])
            res["trpc5_depth_controlled"] = {"LR": round(lr, 2), "df": df, "p": p,
                                             "detect_pct": round(100 * float(t5[sel].mean()), 1)}
        except Exception as e:
            res["trpc5_depth_controlled"] = {"error": str(e)}

        out[label] = res
        print(f"   {label}: n={res['n_lepr']}, eps2 H-D={res['eps2']['HminusD']['mean']}, "
              f"best k={res['mixture']['best_k']}, "
              f"Trpc5 LR={res['trpc5_depth_controlled'].get('LR')}")

    a, b = out["with_doublets"], out["without_doublets"]
    out["comparison"] = {
        "eps2_HminusD_shift": round(b["eps2"]["HminusD"]["mean"] -
                                    a["eps2"]["HminusD"]["mean"], 4),
        "best_k_unchanged": a["mixture"]["best_k"] == b["mixture"]["best_k"],
        "trpc5_still_significant": (b["trpc5_depth_controlled"].get("p") is not None and
                                    b["trpc5_depth_controlled"]["p"] < 0.05),
        "gate_pct_shift": (round(b.get("naive_gate_lepr_pct", 0) -
                                 a.get("naive_gate_lepr_pct", 0), 1))}
    out["interpretation"] = (
        "If the doublet-free values match the full-set values, the full set stays primary "
        "and this becomes a supplementary table. If epsilon-squared, the selected number of "
        "mixture components, or the Trpc5 result changes, the doublet-free set becomes the "
        "primary analysis and the manuscript reports it as such.")

    json.dump(out, open(os.path.join(C.RESULTS, "doublet_sensitivity.json"), "w"), indent=2)
    print("\n" + json.dumps({k: out[k] for k in
                             ("n_lepr_doublet_flagged", "pct_lepr_doublet_flagged",
                              "with_doublets", "without_doublets", "comparison")}, indent=2))


if __name__ == "__main__":
    main()
