#!/usr/bin/env python3
"""
10_conserved_classes.py - effector balance across the conserved VMH classes.

Affinati et al. reduced their 24 mouse VMH clusters to six cross-species conserved classes,
with Lepr expression marking a single dorsomedial class. The nine Lepr clusters are therefore
subdivisions of one conserved class, which means the Kruskal-Wallis in 02 partitions variance
*within* a class rather than between classes. This script uses the `Conserved_cluster`
annotation deposited in the same metadata to ask the better-posed question.

Three outputs:
  1. The mapping: which VMH clusters fall in which conserved class, confirming that the Lepr
     clusters are nested inside one class.
  2. Between-class comparison: effector balance across the six conserved classes, with the
     Lepr class tested against each other class, over N_SEEDS control-gene draws, with sign
     consistency as the robustness criterion.
  3. A variance partition: epsilon-squared between conserved classes against epsilon-squared
     between the nine Lepr clusters within the Lepr class. If the first greatly exceeds the
     second, effector balance is a property of conserved class identity and not of
     subdivision within it, which is a positive statement where the within-class result
     alone was only a small number.

Depth is controlled by repeating every comparison on scores residualized on log10 library
size, since conserved classes differ in capture depth.

One pass over the matrix for the union of rows needed across all seeds, so this is faster
than the per-seed passes in 02.

Output: results/conserved_classes.json, conserved_class_scores.csv, conserved_mapping.csv
"""
import os, csv, json, math, gzip
from collections import defaultdict
import numpy as np
from scipy.stats import kruskal
import config as C
import vmhlib as V

PREFIX, TAG = "GSE172204", "cons"
CLASS_COL = "Conserved_class"   # derived prefix of Conserved_cluster
CLUSTER_COL = "VMH_cluster"


def resid_on_depth(x, logdepth):
    b = np.polyfit(logdepth, x, 1)
    return x - (b[0] * logdepth + b[1])


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)

    # Conserved_cluster holds 17 labels; the six conserved classes of the paper are its
    # prefixes (Lepr, Nfib, Fezf1, Esr1, Dlk1, Foxp2). Doublets_* are excluded.
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["Conserved_class"] = "" if (not raw or raw.startswith("Doublets")) else raw.split("_")[0]

    cols0 = list(next(iter(meta.values())).keys())
    if CLASS_COL not in cols0:
        raise SystemExit(f"{CLASS_COL} not in metadata; columns = {cols0}")

    # ------------------------------------------------------------------ 1. mapping
    cross = defaultdict(lambda: defaultdict(int))
    for m in meta.values():
        k = (m.get(CLUSTER_COL) or "").strip()
        c = (m.get(CLASS_COL) or "").strip()
        if k and c:
            cross[c][k] += 1
    classes = sorted(cross)
    print(">> conserved classes and the VMH clusters nested in them:")
    maprows = []
    for c in classes:
        tot = sum(cross[c].values())
        members = sorted(cross[c], key=lambda k: -cross[c][k])
        n_lepr = sum(v for k, v in cross[c].items() if k.startswith("Lepr"))
        print(f"   {c}: n={tot}, Lepr nuclei={n_lepr}, clusters={members}")
        maprows.append({"conserved_class": c, "n": tot, "n_lepr_nuclei": n_lepr,
                        "vmh_clusters": ";".join(members)})
    with open(os.path.join(C.RESULTS, "conserved_mapping.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(maprows[0].keys()))
        w.writeheader(); w.writerows(maprows)

    lepr_classes = [c for c in classes
                    if sum(v for k, v in cross[c].items() if k.startswith("Lepr")) > 0]
    lepr_class = max(lepr_classes,
                     key=lambda c: sum(v for k, v in cross[c].items() if k.startswith("Lepr")))
    n_lepr_total = sum(sum(v for k, v in cross[c].items() if k.startswith("Lepr"))
                       for c in classes)
    n_in_class = sum(v for k, v in cross[lepr_class].items() if k.startswith("Lepr"))
    print(f">> Lepr nuclei are {100 * n_in_class / n_lepr_total:.1f}% inside class "
          f"'{lepr_class}' ({len(lepr_classes)} classes contain any)")

    # ------------------------------------------------------------- scores by class
    col2group, n_lab, n_unmatched = V.write_col2group(
        PREFIX, meta, bc2col, CLASS_COL, TAG)
    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = g
    cols = sorted(c for c in colsum if c in col2cl)
    cls = np.array([col2cl[c] for c in cols])
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)
    print(f">> {len(cols)} nuclei with a conserved class and counts "
          f"({n_unmatched} labelled barcodes absent from barcodes.tsv)")

    # cluster label per column, for the within-class partition
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    clu = np.array([(meta[col2bc[c]].get(CLUSTER_COL) or "NA").strip() for c in cols])

    # one pass for the union of rows across all seeds
    seeds = list(range(C.RANDOM_SEED, C.RANDOM_SEED + C.N_SEEDS))
    seed_sets = []
    union = set()
    for seed in seeds:
        pr, cr = V.choose_controls(name2row, rowsum, len(cols), seed)
        seed_sets.append((pr, cr))
        union |= {r for rs_ in list(pr.values()) + list(cr.values()) for r in rs_}
    print(f">> extracting {len(union)} gene rows in one pass")
    ent = V.pass_rows(PREFIX, union, TAG, "union")
    expr = defaultdict(dict)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in col2cl:
            expr[c][int(r)] = float(v)

    scores = []
    for (pr, cr), seed in zip(seed_sets, seeds):
        scores.append(V.module_scores(expr, cols, colsum, pr, cr))
        print(f"   scored seed {seed}")

    # ------------------------------------------------- 2. between-class comparison
    out = {"conserved_classes": classes, "lepr_class": lepr_class,
           "pct_lepr_nuclei_in_lepr_class": round(100 * n_in_class / n_lepr_total, 1),
           "n_seeds": C.N_SEEDS,
           "class_sizes": {c: int((cls == c).sum()) for c in classes},
           "median_depth_by_class": {c: float(np.median(depth[cls == c])) for c in classes}}

    rows = []
    for c in classes:
        m = cls == c
        rec = {"conserved_class": c, "n": int(m.sum()), "is_lepr_class": c == lepr_class}
        for lab in ("HYPER", "DEPOL", "HminusD"):
            across = np.array([s[lab][m].mean() for s in scores])
            mean, lo, hi = V.boot_ci(scores[0][lab][m], seed=1)
            rec[f"{lab}_mean"] = round(mean, 4)
            rec[f"{lab}_ci_lo"] = round(lo, 4)
            rec[f"{lab}_ci_hi"] = round(hi, 4)
            rec[f"{lab}_seed_range"] = round(float(across.max() - across.min()), 4)
        rows.append(rec)
    rows.sort(key=lambda r: -r["HminusD_mean"])
    with open(os.path.join(C.RESULTS, "conserved_class_scores.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    pair = {}
    for c in classes:
        if c == lepr_class:
            continue
        a, b = cls == lepr_class, cls == c
        diffs = np.array([float(s["HminusD"][a].mean() - s["HminusD"][b].mean())
                          for s in scores])
        obs, p = V.perm_test(scores[0]["HminusD"][a], scores[0]["HminusD"][b], seed=2)
        pair[f"{lepr_class}_vs_{c}"] = {
            "diff_mean_over_seeds": round(float(diffs.mean()), 4),
            "diff_range": round(float(diffs.max() - diffs.min()), 4),
            "perm_p_seed1": p,
            "sign_consistent_across_seeds": bool(np.all(np.sign(diffs) == np.sign(diffs[0]))),
            "abs_mean_over_range": (round(abs(float(diffs.mean())) / float(diffs.max() - diffs.min()), 2)
                                    if diffs.max() > diffs.min() else None)}
    out["lepr_class_vs_others"] = pair

    # --------------------------------------------------------- 3. variance partition
    part = {}
    for tag, use_resid in (("raw", False), ("depth_residualized", True)):
        between, within = [], []
        for s in scores:
            x = resid_on_depth(s["HminusD"], logdepth) if use_resid else s["HminusD"]
            g = [x[cls == c] for c in classes if (cls == c).sum() > 1]
            H, _ = kruskal(*g)
            between.append(V.eps2_from_kruskal(float(H), len(g), int(sum(len(y) for y in g))))
            m = cls == lepr_class
            sub = [x[m][clu[m] == k] for k in sorted(set(clu[m])) if (clu[m] == k).sum() > 1]
            if len(sub) > 1:
                H2, _ = kruskal(*sub)
                within.append(V.eps2_from_kruskal(float(H2), len(sub),
                                                  int(sum(len(y) for y in sub))))
        part[tag] = {
            "eps2_between_conserved_classes": {
                "mean": round(float(np.mean(between)), 4),
                "min": round(float(np.min(between)), 4),
                "max": round(float(np.max(between)), 4),
                "label": C.eps2_label(float(np.mean(between)))},
            "eps2_between_clusters_within_lepr_class": {
                "mean": round(float(np.mean(within)), 4),
                "min": round(float(np.min(within)), 4),
                "max": round(float(np.max(within)), 4),
                "label": C.eps2_label(float(np.mean(within)))} if within else None,
            "ratio_between_over_within": (round(float(np.mean(between)) / float(np.mean(within)), 2)
                                          if within and np.mean(within) > 0 else None)}
    out["variance_partition"] = part

    out["interpretation"] = (
        "If epsilon-squared between conserved classes substantially exceeds that between "
        "clusters within the Lepr class, effector balance is a property of conserved class "
        "identity rather than of subdivision within a class, and the small within-class "
        "value reported previously is the expected consequence of partitioning inside one "
        "class. If the two are comparable, effector balance does not track the conserved "
        "parcellation either and the graded reading stands at both levels.")

    json.dump(out, open(os.path.join(C.RESULTS, "conserved_classes.json"), "w"), indent=2)
    print("\n" + json.dumps({k: out[k] for k in
                             ("class_sizes", "lepr_class_vs_others", "variance_partition")},
                            indent=2))
    print("\nWritten to", C.RESULTS)


if __name__ == "__main__":
    main()
