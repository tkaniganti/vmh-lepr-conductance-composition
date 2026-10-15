#!/usr/bin/env python3
"""
26_insr_check.py - is there a transcriptomic handle on insulin-responding VMH neurons?

Sohn et al. (2016) describe insulin-inhibited SF1 neurons as a subpopulation distinct from
the leptin-responding ones, engaging K_ATP uniformly. Comparing within-population variance
in those cells against the leptin-receptor population would be informative, but only if the
insulin-responsive group can be defined the way the leptin-receptor group was: by cluster
identity rather than by transcript detection.

This asks whether that is possible. Three things decide it:

  1. Detection rate of Insr across VMH neurons. A receptor detected at a rate comparable to
     Lepr (14.4%) but spread evenly across clusters marks no population; one concentrated in
     particular clusters might.

  2. Concentration across clusters. Reported as the coefficient of variation of per-cluster
     detection rates and as the ratio of highest to lowest, with Lepr and Nr5a1 as reference
     points: Lepr is concentrated by construction, Slc17a6 is not.

  3. Whether any deposited cluster or conserved class is Insr-defined. If the atlas authors
     found an insulin-receptor population they would have named it.

Also reports Irs1, Irs2, Pik3cb and Abcc8/Kcnj11 for context, since the insulin pathway in
these cells runs PI3K p110-beta to K_ATP, and the K_ATP subunits are already known from the
main analysis to be measurable (Abcc8) and not (Kcnj11).

Nothing here is used in the VMH manuscript. It is a feasibility check for a separate
question.

Output: results/insr_check.json
"""
import os, json, warnings
from collections import defaultdict
import numpy as np
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "insr"
GENES = ["Insr", "Lepr", "Nr5a1", "Slc17a6", "Irs1", "Irs2", "Pik3cb", "Abcc8", "Kcnj11"]


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["Conserved_class"] = ("" if (not raw or raw.startswith("Doublets"))
                                else raw.split("_")[0])

    missing = [g for g in GENES if g not in name2row]
    present = [g for g in GENES if g in name2row]
    if missing:
        print(f">> not in the feature table: {missing}", flush=True)

    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "Conserved_class", TAG)
    cs, _ = V.pass_totals(PREFIX, cg, TAG)
    colsum = V.read_two_col(cs)
    c2c = {}
    for l in open(cg):
        a, b = l.split()
        c2c[int(a)] = b
    cols = sorted(c for c in colsum if c in c2c)
    colset = set(cols)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    clu = {c: (meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols}
    cls = {c: c2c[c] for c in cols}
    print(f">> {len(cols)} VMH nuclei, {len(set(clu.values()))} clusters, "
          f"{len(set(cls.values()))} conserved classes", flush=True)

    rows = {name2row[g]: g for g in present}
    ent = V.pass_rows(PREFIX, set(rows), TAG, "insr")
    detected = defaultdict(set)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in colset and float(v) > 0:
            detected[rows[int(r)]].add(c)

    out = {"n_nuclei": len(cols), "genes_absent_from_feature_table": missing, "genes": {}}
    for g in present:
        hit = detected[g]
        overall = 100 * len(hit) / len(cols)
        by_clu = {}
        for c in cols:
            k = clu[c]
            by_clu.setdefault(k, [0, 0])
            by_clu[k][1] += 1
            if c in hit:
                by_clu[k][0] += 1
        pct = {k: 100 * a / b for k, (a, b) in by_clu.items() if b >= 30}
        vals = np.array(list(pct.values()), float)
        top = sorted(pct.items(), key=lambda t: -t[1])[:3]
        bot = sorted(pct.items(), key=lambda t: t[1])[:3]
        cvv = float(vals.std() / vals.mean()) if vals.mean() > 0 else None
        ratio = (float(vals.max() / vals.min()) if vals.min() > 0 else None)
        out["genes"][g] = {
            "detection_pct_all_vmh": round(overall, 2),
            "n_clusters_scored": len(pct),
            "cv_across_clusters": round(cvv, 3) if cvv is not None else None,
            "max_over_min_ratio": round(ratio, 2) if ratio else "lowest cluster at zero",
            "highest_clusters": [(k, round(v, 1)) for k, v in top],
            "lowest_clusters": [(k, round(v, 1)) for k, v in bot]}
        print(f"   {g:<8} {overall:5.2f}% of nuclei | CV across clusters "
              f"{out['genes'][g]['cv_across_clusters']} | top {top[0][0]} {top[0][1]:.1f}%"
              f" | bottom {bot[0][0]} {bot[0][1]:.1f}%", flush=True)

    named = sorted({c.split("_")[0] for c in clu.values() if c != "NA"})
    out["cluster_name_prefixes"] = named
    out["any_insulin_named_population"] = any(
        p.lower().startswith(("insr", "ins")) for p in named)
    print(f"\n>> deposited cluster name prefixes: {', '.join(named)}")
    print(f">> any insulin-receptor-named population: "
          f"{out['any_insulin_named_population']}")

    ref = out["genes"].get("Lepr", {}).get("cv_across_clusters")
    ins = out["genes"].get("Insr", {}).get("cv_across_clusters")
    if ref and ins:
        out["verdict"] = (
            "Insr is less concentrated across clusters than Lepr "
            f"(CV {ins} against {ref}); no cluster-level handle"
            if ins < ref else
            f"Insr concentration (CV {ins}) is comparable to or above Lepr ({ref}); "
            "worth examining which clusters carry it")
        print(f"\n>> {out['verdict']}")

    json.dump(out, open(os.path.join(cfg.RESULTS, "insr_check.json"), "w"), indent=2)
    print("\nWritten to results/insr_check.json")


if __name__ == "__main__":
    main()
