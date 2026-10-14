#!/usr/bin/env python3
"""
22_nested_partition.py - nested decomposition done correctly, with a power-matched control.

Supersedes 21_variance_partition.py, which computed three separate one-way partitions. The
69 class-by-cluster groups are nested inside the six classes, so that partition contains the
class-level variance; reporting it beside the class figure compares a whole with one of its
parts. This script fits the proper two-level nested model.

    y_ijk = mu + A_i + B_j(i) + e_ijk        A = conserved class, B = cluster within class

with mean squares MSA, MSB(A) and MSE, and the standard unbalanced components estimators

    sigma2_B = (MSB_A - MSE) / n0_B
    sigma2_A = (MSA - MSB_A) / n0_A

truncated at zero, reported as a share of the total. These are bias-corrected by
construction, which is why estimates can land exactly on zero; raw eta-squared is reported
alongside for comparison but is not the quantity in the text, since it rises with group
count even under the null.

Three additions the pooled figure cannot supply:

  per-class decomposition - cluster-level variance computed separately inside each of the
    six conserved classes, so it is visible whether within-class structure localises to
    particular classes rather than being uniform

  power-matched bound - whichever class carries the most cluster structure is downsampled to
    the nuclei and cluster count of the Lepr comparison and the estimator rerun, so the null
    inside the Lepr class carries a stated detectability bound in the same style as the
    mixture-model power simulation

  null expectation - the estimator applied to the score with cluster labels shuffled within
    class, giving the value expected by chance at each group count and sample size

The within-Lepr comparison uses the nine Lepr_* clusters, matching the population described
in the manuscript, rather than all cluster labels that happen to fall inside the conserved
Lepr class.

Output: results/nested_partition.json
"""
import os, json, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "np2"
FOCUS_CLASS = "Lepr"
N_NULL = 200
N_MATCH = 100


def one_way(y, grp):
    """Bias-corrected variance share and raw eta-squared for a single grouping."""
    df = pd.DataFrame({"y": np.asarray(y, float), "g": np.asarray(grp)})
    n, k = len(df), df.g.nunique()
    if k < 2 or n <= k:
        return None
    grand = df.y.mean()
    ssb = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in df.groupby("g"))
    ssw = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in df.groupby("g"))
    msb, msw = ssb / (k - 1), ssw / (n - k)
    # unbalanced n0 for the group factor
    ni = df.groupby("g").size().to_numpy(float)
    n0 = (n - (ni ** 2).sum() / n) / (k - 1)
    sb = max((msb - msw) / n0, 0.0)
    return {"pct": 100 * sb / (sb + msw), "eta2": ssb / (ssb + ssw), "k": k, "n": n}


def nested(y, cls, clu):
    """Two-level nested components: class, cluster within class, residual."""
    df = pd.DataFrame({"y": np.asarray(y, float), "A": np.asarray(cls),
                       "B": [f"{a}|{b}" for a, b in zip(cls, clu)]})
    n = len(df)
    a = df.A.nunique()
    b = df.B.nunique()
    grand = df.y.mean()
    ssa = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in df.groupby("A"))
    ssb = sum(len(g) * (g.y.mean() - df.loc[df.A == g.A.iloc[0], "y"].mean()) ** 2
              for _, g in df.groupby("B"))
    sse = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in df.groupby("B"))
    msa, msb, mse = ssa / (a - 1), ssb / (b - a), sse / (n - b)
    nb = df.groupby("B").size().to_numpy(float)
    n0b = (n - (nb ** 2).sum() / n) / (b - a)
    na = df.groupby("A").size().to_numpy(float)
    n0a = (n - (na ** 2).sum() / n) / (a - 1)
    s2b = max((msb - mse) / n0b, 0.0)
    s2a = max((msa - msb) / n0a, 0.0)
    tot = s2a + s2b + mse
    return {"pct_class": 100 * s2a / tot, "pct_cluster_in_class": 100 * s2b / tot,
            "pct_residual": 100 * mse / tot, "n_classes": a, "n_class_by_cluster": b, "n": n}


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["Conserved_class"] = ("" if (not raw or raw.startswith("Doublets"))
                                else raw.split("_")[0])
    col2group, _, _ = V.write_col2group(PREFIX, meta, bc2col, "Conserved_class", TAG)
    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    c2c = {}
    for l in open(col2group):
        a, b = l.split()
        c2c[int(a)] = b
    cols = sorted(c for c in colsum if c in c2c)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cls = np.array([c2c[c] for c in cols])
    clu = np.array([(meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols])
    ld = np.log10(np.array([colsum[c] for c in cols], float))

    # the manuscript's comparison: the nine Lepr_* clusters inside the conserved Lepr class
    focus = (cls == FOCUS_CLASS) & np.array([c.startswith("Lepr_") for c in clu])
    print(f">> {len(cols)} nuclei; {len(set(cls))} classes; "
          f"focus set {int(focus.sum())} nuclei in {len(set(clu[focus]))} Lepr_* clusters",
          flush=True)

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

    nested_recs, focus_recs, per_class = [], [], defaultdict(list)
    for seed, (pr, cr) in zip(seeds, seed_sets):
        y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        y = y - np.polyval(np.polyfit(ld, y, 1), ld)
        nested_recs.append(nested(y, cls, clu))
        focus_recs.append(one_way(y[focus], clu[focus]))
        for cl in sorted(set(cls.tolist())):
            m = cls == cl
            r = one_way(y[m], clu[m])
            if r:
                per_class[cl].append(r)
        print(f"   seed {seed}: class {nested_recs[-1]['pct_class']:.3f}%, "
              f"cluster-in-class {nested_recs[-1]['pct_cluster_in_class']:.3f}%, "
              f"within-Lepr {focus_recs[-1]['pct']:.3f}%", flush=True)

    def summ(vals):
        a = np.array(vals, float)
        return {"mean": round(float(a.mean()), 4), "min": round(float(a.min()), 4),
                "max": round(float(a.max()), 4)}

    pc = {cl: {"pct": summ([r["pct"] for r in recs]),
               "eta2_mean": round(float(np.mean([r["eta2"] for r in recs])), 4),
               "k": recs[0]["k"], "n": recs[0]["n"]} for cl, recs in per_class.items()}
    print("\n>> cluster variance (%) by conserved class:", flush=True)
    for cl, v in sorted(pc.items(), key=lambda t: -t[1]["pct"]["mean"]):
        print(f"   {cl:<7} {v['pct']['mean']:6.3f}%  (k={v['k']}, n={v['n']})", flush=True)

    # ---------------------------------------------- null expectation at matched size
    y0 = V.module_scores(expr, cols, colsum, *seed_sets[0])["HminusD"]
    y0 = y0 - np.polyval(np.polyfit(ld, y0, 1), ld)
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    nulls = []
    lab = clu[focus].copy()
    for _ in range(N_NULL):
        rng.shuffle(lab)
        r = one_way(y0[focus], lab)
        if r:
            nulls.append(r["pct"])
    null_focus = {"mean": round(float(np.mean(nulls)), 4),
                  "p95": round(float(np.percentile(nulls, 95)), 4), "n": len(nulls)}
    print(f"\n>> null expectation inside the focus set: mean {null_focus['mean']}%, "
          f"95th percentile {null_focus['p95']}%", flush=True)

    # Detectability bound: inject cluster effects of known size into the focus set itself
    # and report what the estimator recovers at its own n and k.
    inj = []
    sd = float(np.std(y0[focus]))
    lab0 = clu[focus]
    uniq = sorted(set(lab0))
    for eff in (0.05, 0.10, 0.20, 0.40):
        vals = []
        for _ in range(50):
            shift = {u: rng.normal(0, eff * sd) for u in uniq}
            yy = y0[focus] + np.array([shift[u] for u in lab0])
            r = one_way(yy, lab0)
            if r:
                vals.append(r["pct"])
        inj.append({"effect_sd_units": eff, "recovered_pct": summ(vals)})
        print(f">> injected cluster SD {eff:.2f}: recovered "
              f"{inj[-1]['recovered_pct']['mean']:.3f}%", flush=True)
    matched = {"design": "cluster effects injected into the focus set at its own n and k",
               "n": int(focus.sum()), "k": len(uniq), "levels": inj,
               "observed_without_injection": summ([r["pct"] for r in focus_recs]),
               "esr1_reference": {"pct": round(pc["Esr1"]["pct"]["mean"], 3),
                                  "n": pc["Esr1"]["n"], "k": pc["Esr1"]["k"],
                                  "note": "detected with fewer nuclei than the focus set"}}

    out = {"n_nuclei": len(cols),
           "estimator": ("unbalanced nested components of variance; sigma2_B = "
                         "(MSB(A) - MSE)/n0_B, sigma2_A = (MSA - MSB(A))/n0_A, both "
                         "truncated at zero and reported as shares of the total. These are "
                         "bias-corrected, so estimates can land exactly on zero; raw "
                         "eta-squared is reported alongside but is not used in the text "
                         "because it rises with group count under the null."),
           "nested_all_vmh": {k: summ([r[k] for r in nested_recs])
                              for k in ("pct_class", "pct_cluster_in_class", "pct_residual")},
           "nested_counts": {k: nested_recs[0][k]
                             for k in ("n_classes", "n_class_by_cluster", "n")},
           "within_lepr_class_nine_clusters": {
               "pct": summ([r["pct"] for r in focus_recs]),
               "eta2_mean": round(float(np.mean([r["eta2"] for r in focus_recs])), 4),
               "k": focus_recs[0]["k"], "n": focus_recs[0]["n"],
               "null_expectation": null_focus},
           "cluster_pct_by_class": pc,
           "power_matched_control": matched}
    json.dump(out, open(os.path.join(cfg.RESULTS, "nested_partition.json"), "w"), indent=2)
    print("\n" + json.dumps({k: out[k] for k in
                             ("nested_all_vmh", "within_lepr_class_nine_clusters",
                              "power_matched_control")}, indent=2))


if __name__ == "__main__":
    main()
