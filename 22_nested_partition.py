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


def _ems_coefficients(df):
    """c1, c2, c3 for an UNBALANCED two-level nested design.

    E(MSE)=s2e, E(MSB(A))=s2e+c1*s2B, E(MSA)=s2e+c2*s2B+c3*s2A
    (Sokal & Rohlf, Biometry, Box 10.6; Lynch & Walsh 1998, Eq. 18.32).
    """
    n = len(df)
    a = df.A.nunique()
    b = df.B.nunique()
    t_within = t_grand = t_cls = 0.0
    for _, gi in df.groupby("A"):
        n_i = len(gi)
        t_cls += n_i ** 2 / n
        for _, gij in gi.groupby("B"):
            n_ij = len(gij)
            t_within += n_ij ** 2 / n_i
            t_grand += n_ij ** 2 / n
    return ((n - t_within) / (b - a),
            (t_within - t_grand) / (a - 1),
            (n - t_cls) / (a - 1), a, b, n)


def nested(y, cls, clu):
    """Two-level nested components: class, cluster within class, residual.

    BOTH COEFFICIENTS WERE WRONG BEFORE 2 NOV 2026 and are fixed here.

    Lower level: n0b divided each squared cell size by the GRAND total rather
    than by that cell's own class total, so the divisor was too large and
    pct_cluster_in_class was understated.

    Upper level: s2a was (MSA - MSB) / n0a, which assumes c2 == c1. That holds
    only under balance; here they differ by close to an order of magnitude, so
    cluster variance leaked into the class term and inflated it. n0a itself was
    right and equals c3.

    Why the balanced control in 38_upper_level_recovery.py passed: under balance
    c2 DOES equal c1, so a clean balanced control is exactly what a formula valid
    only under balance produces. It could not have failed. The Limitations
    paragraph read that as the imbalance being unfixable; the estimator was wrong.
    nested_positive_control() below therefore runs on the REAL design.
    """
    df = pd.DataFrame({"y": np.asarray(y, float), "A": np.asarray(cls),
                       "B": [f"{a}|{b}" for a, b in zip(cls, clu)]})
    c1, c2, c3, a, b, n = _ems_coefficients(df)
    grand = df.y.mean()
    ssa = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in df.groupby("A"))
    ssb = sum(len(g) * (g.y.mean() - df.loc[df.A == g.A.iloc[0], "y"].mean()) ** 2
              for _, g in df.groupby("B"))
    sse = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in df.groupby("B"))
    msa, msb, mse = ssa / (a - 1), ssb / (b - a), sse / (n - b)
    s2b_raw = (msb - mse) / c1
    s2b = max(s2b_raw, 0.0)
    s2a_raw = (msa - mse - c2 * s2b) / c3
    s2a = max(s2a_raw, 0.0)
    tot = s2a + s2b + mse
    return {"pct_class": 100 * s2a / tot, "pct_cluster_in_class": 100 * s2b / tot,
            "pct_residual": 100 * mse / tot, "n_classes": a,
            "n_class_by_cluster": b, "n": n,
            "ems_c1": round(c1, 3), "ems_c2": round(c2, 3), "ems_c3": round(c3, 3),
            "s2a_before_truncation": round(s2a_raw, 6),
            "s2b_before_truncation": round(s2b_raw, 6)}


def class_component_null(y_draws, cls, clu, seed=0, n_perm=200):
    """Chance expectation for the class component.

    Every other component in this manuscript is reported as a percentile of its
    own null and the Methods say so; the class component was withheld rather
    than tested, so it had none. Reporting it against zero while everything else
    is reported against a null is the asymmetry that made the transcript scan
    untenable.

    Construction: reassign whole CELLS to classes at random, preserving the
    number of cells per class and every cell's size. That breaks the class
    grouping while leaving the nested design, the cell sizes and the
    within-cell structure exactly as observed, so the null asks whether cells
    group by class more than an arbitrary grouping of the same cells would.
    """
    cls = np.asarray(cls)
    clu = np.asarray(clu)
    cell = np.array([f"{a}|{b}" for a, b in zip(cls, clu)])
    cells = sorted(set(cell))
    cell_class = {c: c.split("|")[0] for c in cells}
    labels = [cell_class[c] for c in cells]

    rng = np.random.default_rng(seed)
    obs, null = [], []
    for y in y_draws:
        obs.append(nested(y, cls, clu)["pct_class"])
    for _ in range(n_perm):
        perm = list(rng.permutation(labels))
        remap = dict(zip(cells, perm))
        # cell identity is carried in clu so the permuted classes cannot collide
        null.append(nested(y_draws[0], np.array([remap[c] for c in cell]),
                           cell)["pct_class"])
    null = np.array(null)
    o = float(np.mean(obs))
    return {"observed_pct_class_mean": round(o, 4),
            "n_draws": len(y_draws), "n_permutations": int(n_perm),
            "null_mean": round(float(null.mean()), 4),
            "null_p95": round(float(np.percentile(null, 95)), 4),
            "null_max": round(float(null.max()), 4),
            "percentile_of_null": round(float((null < o).mean() * 100), 1),
            "exceeds_null_p95": bool(o > np.percentile(null, 95)),
            "note": ("whole cells reassigned to classes at random, preserving "
                     "cells per class and every cell size; the observed value is "
                     "the mean over control-gene draws and the null is built on "
                     "the first draw")}


def nested_positive_control(cls, clu, seed=0, reps=300):
    """Inject KNOWN components into the REAL design and check recovery.

    Untruncated, so max(., 0) cannot mask an upward bias at zero injection.
    """
    cls = np.asarray(cls)
    clu = np.asarray(clu)
    df = pd.DataFrame({"A": cls, "B": [f"{x}|{z}" for x, z in zip(cls, clu)]})
    c1, c2, c3, a, b, n = _ems_coefficients(df)
    out = {"design": {"n": int(n), "n_classes": int(a), "n_cells": int(b),
                      "c1": round(c1, 3), "c2": round(c2, 3), "c3": round(c3, 3)},
           "levels": []}
    cells = list(df.B)
    for true_a, true_b in ((0.0, 0.04), (0.0289, 0.04), (0.0, 0.12)):
        ga, gb = [], []
        for rep in range(reps):
            r = np.random.default_rng(seed + rep)
            ae = {i: r.normal(0, np.sqrt(true_a)) for i in set(cls)}
            be = {j: r.normal(0, np.sqrt(true_b)) for j in set(cells)}
            yy = (np.array([ae[i] for i in cls]) + np.array([be[j] for j in cells])
                  + r.normal(0, 1.0, n))
            d2 = df.assign(y=yy)
            grand = yy.mean()
            ssa = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in d2.groupby("A"))
            ssb = sum(len(g) * (g.y.mean()
                                - d2.loc[d2.A == g.A.iloc[0], "y"].mean()) ** 2
                      for _, g in d2.groupby("B"))
            sse = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in d2.groupby("B"))
            msa, msb, mse = ssa / (a - 1), ssb / (b - a), sse / (n - b)
            s2b = (msb - mse) / c1
            gb.append(s2b)
            ga.append((msa - mse - c2 * s2b) / c3)
        # Report the standard error beside the mean. At 30 reps all three
        # non-zero rows came back 6 to 8% high, which is within noise per row
        # but is three of three in one direction. With the class component now
        # reportable rather than withheld, a residual bias of that size matters,
        # so the mean needs its own uncertainty attached.
        out["levels"].append({"injected_s2A": true_a, "injected_s2B": true_b,
                              "recovered_s2A": round(float(np.mean(ga)), 5),
                              "recovered_s2A_se": round(
                                  float(np.std(ga, ddof=1) / np.sqrt(len(ga))), 5),
                              "recovered_s2B": round(float(np.mean(gb)), 5),
                              "recovered_s2B_se": round(
                                  float(np.std(gb, ddof=1) / np.sqrt(len(gb))), 5)})
    out["passes"] = bool(all(
        abs(l["recovered_s2A"] - l["injected_s2A"]) < 0.01
        and abs(l["recovered_s2B"] - l["injected_s2B"]) < 0.02
        for l in out["levels"]))
    return out


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
    y_draws = []   # kept for class_component_null, which needs the scores themselves
    for seed, (pr, cr) in zip(seeds, seed_sets):
        y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        y = y - np.polyval(np.polyfit(ld, y, 1), ld)
        y_draws.append(y)
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
           "estimator": ("unbalanced nested components of variance. "
                         "sigma2_B = (MSB(A) - MSE)/c1 and "
                         "sigma2_A = (MSA - MSE - c2*sigma2_B)/c3, with the "
                         "expected-mean-square coefficients computed for this design "
                         "(Sokal & Rohlf Box 10.6; Lynch & Walsh Eq. 18.32). Both are "
                         "truncated at zero and reported as shares of the total, so "
                         "estimates can land exactly on zero; raw eta-squared is "
                         "reported alongside but is not used in the text because it "
                         "rises with group count under the null. BEFORE 2 Nov 2026 "
                         "this used sigma2_A = (MSA - MSB(A))/n0_A, valid only under "
                         "balance, and a lower-level divisor formed with the grand "
                         "total rather than each class's own; the class term was "
                         "inflated and the cluster-within-class term understated."),
           "estimator_positive_control": nested_positive_control(cls, clu),
           "class_component_null": class_component_null(y_draws, cls, clu),
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
