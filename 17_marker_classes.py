#!/usr/bin/env python3
"""
17_marker_classes.py - do the transcriptomic classes match the RNAscope map?

The RNAscope atlas gives the composition of the recorded SF1/LepR population directly:
within SF1 cells that are also LepR positive, roughly 79% are Gck positive and 81% Bdnf
positive, but only about 31% Esr1 positive. If the conserved Lepr class is the transcriptomic
counterpart of that population, it should rank high for Gck and Bdnf and low for Esr1 against
the other VMH classes.

That is a prediction about ORDERING ACROSS CLASSES, not about absolute co-detection. Droplet
dropout puts per-nucleus co-detection far below the in-situ values (Lepr itself is detected in
23.5% of nuclei inside Lepr clusters), so the two assays cannot be compared on their face.
Detection percentages are reported for completeness and must not be read as the RNAscope
fractions measured again.

Three groupings: the six conserved VMH classes, the 24 VMH clusters, and the Rupp LepR
populations. Expression is per-group pseudobulk log-CP10k with a bootstrap CI over nuclei,
and each marker also gets a detection logistic on group given log10 depth, since classes
differ in capture depth.

Output: results/marker_classes_<tag>.csv, marker_classes.json
"""
import os, csv, json, math
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2 as chi2dist
import config as C
import vmhlib as V

MARKERS = ["Gck", "Bdnf", "Esr1", "Nr5a1", "Lepr", "Slc17a6"]
# fractions of the recorded SF1-and-LepR population carrying each marker, from the
# fluorescent triple-label series; the prediction is on rank, not on these values
RNASCOPE_WITHIN_SF1_LEPR = {"Gck": 0.791, "Bdnf": 0.807, "Esr1": 0.311}
N_BOOT = 2000


def boot_ci(x, seed=0):
    rng = np.random.default_rng(seed)
    x = np.asarray(x, float)
    m = rng.choice(x, (N_BOOT, x.size), replace=True).mean(1)
    return float(x.mean()), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def lr_group_given_depth(y, logdepth, group):
    Xd = sm.add_constant(logdepth.reshape(-1, 1), has_constant="add")
    m0 = sm.Logit(y, Xd).fit(disp=0, maxiter=200)
    D = pd.get_dummies(group, drop_first=True).to_numpy(dtype=float)
    m1 = sm.Logit(y, np.column_stack([Xd, D])).fit(disp=0, maxiter=200)
    lr = 2 * (m1.llf - m0.llf)
    return float(lr), int(D.shape[1]), float(chi2dist.sf(max(lr, 0.0), D.shape[1]))


def run(prefix, tag, label_col, relabel, keep, focus):
    print(f"\n=== {prefix} ({tag}) ===")
    feats, name2row = V.load_features(prefix)
    bcs, bc2col = V.load_barcodes(prefix)
    meta = V.load_metadata(prefix)
    if label_col == "Conserved_class":
        for m in meta.values():
            raw = (m.get("Conserved_cluster") or "").strip()
            m["Conserved_class"] = ("" if (not raw or raw.startswith("Doublets"))
                                    else raw.split("_")[0])

    present = [g for g in MARKERS if g in name2row]
    missing = [g for g in MARKERS if g not in name2row]
    if missing:
        print("   not in features:", missing)
    rows_needed = {name2row[g] for g in present}

    col2group, _, _ = V.write_col2group(prefix, meta, bc2col, label_col, f"{tag}_mk", keep=keep)
    cs, rs = V.pass_totals(prefix, col2group, f"{tag}_mk")
    colsum = V.read_two_col(cs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = relabel(g)
    cols = sorted(c for c in colsum if col2cl.get(c))
    grp = np.array([col2cl[c] for c in cols])
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)
    groups = sorted(set(grp.tolist()))
    cidx = {c: i for i, c in enumerate(cols)}
    print("   n:", {g: int((grp == g).sum()) for g in groups})

    ent = V.pass_rows(prefix, rows_needed, f"{tag}_mk", "markers17")
    row2gene = {name2row[g]: g for g in present}
    counts = {g: np.zeros(len(cols)) for g in present}
    for l in open(ent):
        r, c, v = l.split()
        g = row2gene.get(int(r))
        i = cidx.get(int(c))
        if g is not None and i is not None:
            counts[g][i] = float(v)

    tab, rank_summary = [], {}
    for g in present:
        ln = np.log1p(counts[g] / depth * C.NORM_TARGET)
        per = {}
        for q in groups:
            m = grp == q
            mean, lo, hi = boot_ci(ln[m], seed=3)
            per[q] = mean
            tab.append({"group": q, "gene": g, "n": int(m.sum()),
                        "log_cp10k_mean": round(mean, 4),
                        "ci_lo": round(lo, 4), "ci_hi": round(hi, 4),
                        "detect_pct": round(100 * float((counts[g][m] > 0).mean()), 1),
                        "median_depth": int(np.median(depth[m]))})
        order = sorted(groups, key=lambda q: -per[q])
        try:
            lr, df, p = lr_group_given_depth((counts[g] > 0).astype(float), logdepth, grp)
        except Exception:
            lr, df, p = None, None, None
        rank_summary[g] = {
            "ranking_high_to_low": order,
            "focus_rank": (order.index(focus) + 1) if focus in order else None,
            "n_groups": len(order),
            "focus_mean": round(per.get(focus, float("nan")), 4),
            "detection_LR_group_given_depth": None if lr is None else round(lr, 2),
            "detection_p": p}
        if focus in order:
            print(f"   {g:<9} {focus} ranks {order.index(focus) + 1}/{len(order)} "
                  f"(mean {per[focus]:.3f}); order: {', '.join(order)}")

    with open(os.path.join(C.RESULTS, f"marker_classes_{tag}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(tab[0].keys()))
        w.writeheader(); w.writerows(tab)

    verdict = None
    if all(k in rank_summary for k in ("Gck", "Bdnf", "Esr1")) and focus:
        n = rank_summary["Gck"]["n_groups"]
        hi_half = lambda r: r is not None and r <= n / 2
        verdict = {
            "prediction": ("Lepr-equivalent class high for Gck and Bdnf, low for Esr1, "
                           "matching the RNAscope composition of the recorded population"),
            "Gck_rank": rank_summary["Gck"]["focus_rank"],
            "Bdnf_rank": rank_summary["Bdnf"]["focus_rank"],
            "Esr1_rank": rank_summary["Esr1"]["focus_rank"],
            "of_groups": n,
            "supported": bool(hi_half(rank_summary["Gck"]["focus_rank"]) and
                              hi_half(rank_summary["Bdnf"]["focus_rank"]) and
                              not hi_half(rank_summary["Esr1"]["focus_rank"]))}

    return {"dataset": prefix, "tag": tag, "focus_group": focus,
            "n_per_group": {g: int((grp == g).sum()) for g in groups},
            "median_depth": {g: float(np.median(depth[grp == g])) for g in groups},
            "markers": rank_summary, "verdict": verdict,
            "rnascope_reference_within_SF1_LepR": RNASCOPE_WITHIN_SF1_LEPR,
            "note": ("Detection percentages are dropout-limited and are NOT comparable to the "
                     "RNAscope fractions; the prediction is on ranking across groups.")}


def main():
    out = {}
    out["conserved"] = run("GSE172204", "cons6", "Conserved_class",
                           relabel=lambda g: g, keep=lambda g: bool(g), focus="Lepr")
    out["clusters24"] = run("GSE172204", "vmh24mk", "VMH_cluster",
                            relabel=lambda g: g.split("_")[0], keep=lambda g: bool(g),
                            focus="Lepr")
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    out["rupp"] = run("GSE172460", "ruppmk", "Lepr_cluster",
                      relabel=lambda g: inv.get(g, g), keep=lambda g: bool(g),
                      focus="VMH_SF1")
    json.dump(out, open(os.path.join(C.RESULTS, "marker_classes.json"), "w"), indent=2)
    print("\n== verdicts ==")
    for k, v in out.items():
        print(k, json.dumps(v["verdict"]))


if __name__ == "__main__":
    main()
