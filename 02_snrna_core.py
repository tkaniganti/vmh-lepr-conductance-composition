#!/usr/bin/env python3
"""
02_snrna_core.py - population definition, detection, pseudobulk, per-cell module scores.

Fixes relative to the previous pipeline:
  * the naive double-positive gate is reported separately for all VMH neurons and for
    LepR clusters, so the 344 (VMH-wide) is never paired with the LepR denominator;
  * labelled barcodes absent from barcodes.tsv are counted and reported, which is the
    source of the 6048/6049 and 2913/2914 discrepancies;
  * epsilon-squared is computed over N_SEEDS control-gene draws and reported as
    mean with range, so the ranges in the manuscript have a stated origin.

Outputs: results/detection.json, cluster_scores.csv, percell_stats.json,
         trpc5_detection.csv, and data/percell_scores.npz for downstream scripts.
"""
import os, json, csv, math, gzip
from collections import defaultdict
import numpy as np
from scipy.stats import kruskal, chi2_contingency
import config as C
import vmhlib as V

PREFIX = "GSE172204"
TAG = "vmh"


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)

    col2group, n_lab, n_unmatched = V.write_col2group(
        PREFIX, meta, bc2col, "VMH_cluster", TAG)
    print(f">> VMH-labelled barcodes: {n_lab}; not found in barcodes.tsv: {n_unmatched}")

    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = g
    cols = sorted(colsum)                      # columns with at least one count
    zero_count_cols = sorted(set(col2cl) - set(colsum))
    print(f">> columns with counts: {len(cols)}; labelled but all-zero: {len(zero_count_cols)}")

    clu = np.array([col2cl[c] for c in cols])
    lepr_clusters = sorted({g for g in clu if g.startswith("Lepr")})

    # ---------------------------------------------------------------- detection
    marker_rows = {g: name2row[g] for g in ("Nr5a1", "Lepr", "Slc17a6") if g in name2row}
    ent = V.pass_rows(PREFIX, set(marker_rows.values()), TAG, "markers")
    row2gene = {v: k for k, v in marker_rows.items()}
    pos = defaultdict(set)
    for l in open(ent):
        r, c, v = l.split()
        if float(v) > 0:
            pos[row2gene[int(r)]].add(int(c))

    all_neuron_cols = set()
    for bc, m in meta.items():
        if (m.get("Neuron_cluster") or "").strip() and bc in bc2col:
            all_neuron_cols.add(bc2col[bc])
    vmh_cols = set(col2cl)
    lepr_cols = {c for c, g in col2cl.items() if g.startswith("Lepr")}

    def rate(cset, g):
        return round(100 * len(cset & pos[g]) / len(cset), 1)

    def dpos(cset):
        return len(cset & pos["Nr5a1"] & pos["Lepr"])

    det = {
        "n_total_cells": len(meta),
        "n_neurons": len(all_neuron_cols),
        "n_vmh_labelled": len(vmh_cols),
        "n_vmh_with_counts": len(cols),
        "n_vmh_labelled_not_in_barcodes": n_unmatched,
        "n_lepr_labelled": len(lepr_cols),
        "vmh_detect": {g: rate(vmh_cols, g) for g in marker_rows},
        "lepr_detect": {g: rate(lepr_cols, g) for g in marker_rows},
        "naive_gate": {
            "all_neurons": dpos(all_neuron_cols),
            "vmh_all": dpos(vmh_cols),
            "lepr_clusters": dpos(lepr_cols),
            "lepr_pct": round(100 * dpos(lepr_cols) / len(lepr_cols), 1),
        },
        "cluster_sizes": {g: int((clu == g).sum()) for g in sorted(set(clu))},
    }
    json.dump(det, open(os.path.join(C.RESULTS, "detection.json"), "w"), indent=2)
    print(">> naive gate: VMH-wide {vmh_all}, LepR clusters {lepr_clusters} "
          "({lepr_pct}%)".format(**det["naive_gate"]))

    # ---------------------------------------------------------------- pseudobulk
    panel_genes = [g for p in ("HYPER", "DEPOL", "MARKERS") for g in C.PANELS[p]]
    panel_rows_all = {name2row[g]: g for g in panel_genes if g in name2row}
    ent_p = V.pass_rows(PREFIX, set(panel_rows_all), TAG, "panels")
    gsum = defaultdict(lambda: defaultdict(float))
    totals = defaultdict(float)
    expr = defaultdict(dict)
    for l in open(ent_p):
        r, c, v = l.split()
        r, c, v = int(r), int(c), float(v)
        if c in col2cl:
            expr[c][r] = v
            gsum[col2cl[c]][panel_rows_all[r]] += v
    for c in cols:
        totals[col2cl[c]] += colsum[c]

    def ln(cl, g):
        return math.log1p(gsum[cl].get(g, 0.0) / totals[cl] * C.NORM_TARGET)

    rows = []
    for k in sorted(totals):
        H = float(np.mean([ln(k, g) for g in C.PANELS["HYPER"] if g in name2row]))
        D = float(np.mean([ln(k, g) for g in C.PANELS["DEPOL"] if g in name2row]))
        rows.append({"cluster": k, "is_lepr": k.startswith("Lepr"),
                     "HYPER": round(H, 4), "DEPOL": round(D, 4), "HminusD": round(H - D, 4),
                     **{g: round(ln(k, g), 4) for g in C.KEY_CHANNELS if g in name2row}})
    rows.sort(key=lambda r: -r["HminusD"])
    with open(os.path.join(C.RESULTS, "cluster_scores.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader(); w.writerows(rows)

    # ------------------------------------------------- per-cell scores over seeds
    lepr_mask = np.array([g.startswith("Lepr") for g in clu])
    lepr_cols_ord = [c for c, m in zip(cols, lepr_mask) if m]
    lepr_lab = clu[lepr_mask]

    per_seed = defaultdict(list)
    score_store = None
    for seed in range(C.RANDOM_SEED, C.RANDOM_SEED + C.N_SEEDS):
        panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(cols), seed)
        need = {r for rs in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs}
        missing = need - set(panel_rows_all)
        if missing:
            e = V.pass_rows(PREFIX, need, TAG, f"seed{seed}")
            for l in open(e):
                r, c, v = l.split()
                c = int(c)
                if c in col2cl:
                    expr[c][int(r)] = float(v)
        sc = V.module_scores(expr, cols, colsum, panel_rows, ctrl_rows)
        if seed == C.RANDOM_SEED:
            score_store = {k: v[lepr_mask] for k, v in sc.items()}
        for lab in ("HYPER", "DEPOL", "HminusD"):
            arr = sc[lab][lepr_mask]
            groups = [arr[lepr_lab == k] for k in lepr_clusters]
            Hst, p = kruskal(*groups)
            n = int(sum(g.size for g in groups))
            per_seed[lab].append({"H": float(Hst), "p": float(p),
                                  "eps2": V.eps2_from_kruskal(float(Hst), len(groups), n),
                                  "n": n})
        print(f"   seed {seed}: eps2 H-D = {per_seed['HminusD'][-1]['eps2']:.4f}")

    kw = {}
    for lab, recs in per_seed.items():
        e = np.array([r["eps2"] for r in recs])
        kw[lab] = {"n_seeds": C.N_SEEDS, "n": recs[0]["n"],
                   "eps2_mean": round(float(e.mean()), 4),
                   "eps2_min": round(float(e.min()), 4),
                   "eps2_max": round(float(e.max()), 4),
                   "label": C.eps2_label(float(e.mean())),
                   "p_max": max(r["p"] for r in recs)}

    # ---------------------------------------------------------------- Trpc5 (raw)
    t5 = name2row["Trpc5"]
    e5 = V.pass_rows(PREFIX, {t5}, TAG, "trpc5")
    t5pos = set()
    for l in open(e5):
        r, c, v = l.split()
        if float(v) > 0:
            t5pos.add(int(c))
    detrows, tbl = [], []
    for k in lepr_clusters:
        idx = [c for c in lepr_cols_ord if col2cl[c] == k]
        p = sum(1 for c in idx if c in t5pos)
        detrows.append({"cluster": k, "pos": p, "total": len(idx),
                        "pct": round(100 * p / len(idx), 1)})
        tbl.append([p, len(idx) - p])
    chi2, pv, dof, _ = chi2_contingency(np.array(tbl))
    kw["trpc5_chi2_uncorrected"] = {"chi2": round(float(chi2), 1), "p": float(pv), "dof": int(dof),
                                    "note": "not depth-corrected; see 03_trpc5_depth.py"}
    with open(os.path.join(C.RESULTS, "trpc5_detection.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["cluster", "pos", "total", "pct"])
        w.writeheader(); w.writerows(detrows)

    json.dump(kw, open(os.path.join(C.RESULTS, "percell_stats.json"), "w"), indent=2)
    np.savez(V.dpath("percell_scores.npz"),
             cols=np.array(lepr_cols_ord), labels=lepr_lab,
             **{k: v for k, v in score_store.items()})
    print("\n== per-cell effect sizes (over %d control-gene seeds) ==" % C.N_SEEDS)
    print(json.dumps({k: v for k, v in kw.items() if k != "trpc5_chi2_uncorrected"}, indent=2))


if __name__ == "__main__":
    main()
