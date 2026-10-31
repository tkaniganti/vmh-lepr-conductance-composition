#!/usr/bin/env python3
"""
31_worked_example.py - trace the module-score calculation from raw counts to a reported value.

Every number in the effector-composition sections rests on one quantity: the difference
between a hyperpolarising and a depolarising module score, each corrected against
expression-matched controls. The calculation involves a normalisation, a log transform, a
per-gene z-scoring across nuclei, two averaging steps and two subtractions, and the units
of the result follow from the z-scoring rather than from the log. A reader cannot check any
reported value without knowing that chain.

This reproduces it on a small, named subset: a few genes in a few nuclei, printing every
intermediate. The final axis value for those nuclei is recomputed by the pipeline's own
function and compared, so the worked example is demonstrably the same calculation and not
a parallel reimplementation that might drift from it.

Output: results/worked_example.json, and a table for the supplementary note.
"""
import os, json, math, warnings
from collections import defaultdict
import numpy as np
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "wex"
N_NUCLEI = 5                     # nuclei to show in full
SHOW_GENES = 3                   # panel genes to show per module


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""

    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_grp", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cols = sorted(c for c in colsum if c in col2bc
                  and (meta[col2bc[c]].get("VMH_cluster") or "").strip())
    clu = {c: meta[col2bc[c]]["VMH_cluster"].strip() for c in cols}
    print(f">> {len(cols)} VMH nuclei", flush=True)

    # the pipeline's own control draw, first seed, so this is the real calculation
    panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(cols), cfg.RANDOM_SEED)
    need = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs_}
    ent = V.pass_rows(PREFIX, need, TAG, "union")
    expr = defaultdict(dict)
    colset = set(cols)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in colset:
            expr[c][int(r)] = float(v)

    row2name = {v: k for k, v in name2row.items()}
    scores = V.module_scores(expr, cols, colsum, panel_rows, ctrl_rows)

    # pick nuclei spanning the axis so the example is not all near zero
    order = np.argsort(scores["HminusD"])
    pick = [cols[order[i]] for i in
            np.linspace(0, len(cols) - 1, N_NUCLEI).astype(int)]

    # ---------------------------------------------------------------- step by step
    # recompute the per-gene z-scores exactly as vmhlib does, to show the intermediates
    allrows = sorted(need)
    ln = {}
    for r in allrows:
        v = np.fromiter((math.log1p(expr.get(c, {}).get(r, 0.0) / colsum[c] * cfg.NORM_TARGET)
                         for c in cols), dtype=float, count=len(cols))
        sd = v.std()
        ln[r] = {"mean": float(v.mean()), "sd": float(sd),
                 "z": (v - v.mean()) / (sd if sd > 0 else 1.0)}

    idx = {c: i for i, c in enumerate(cols)}
    out = {"n_nuclei": len(cols), "seed": cfg.RANDOM_SEED,
           "norm_target": cfg.NORM_TARGET, "nuclei": []}

    print("\n" + "=" * 78)
    print("WORKED EXAMPLE: raw count to effector-balance axis")
    print("=" * 78)
    for c in pick:
        i = idx[c]
        rec = {"barcode": col2bc[c], "cluster": clu[c],
               "total_counts": colsum[c], "genes": {}}
        print(f"\nNucleus {col2bc[c]}  ({clu[c]}, {colsum[c]:,.0f} total counts)")
        for mod in ("HYPER", "DEPOL"):
            print(f"  {mod} panel:")
            for r in panel_rows[mod][:SHOW_GENES]:
                g = row2name.get(r, str(r))
                raw = expr.get(c, {}).get(r, 0.0)
                cp10k = raw / colsum[c] * cfg.NORM_TARGET
                lg = math.log1p(cp10k)
                z = ln[r]["z"][i]
                rec["genes"][g] = {"raw": raw, "cp10k": round(cp10k, 3),
                                   "log1p": round(lg, 4),
                                   "gene_mean": round(ln[r]["mean"], 4),
                                   "gene_sd": round(ln[r]["sd"], 4),
                                   "z": round(float(z), 4)}
                print(f"    {g:<9} raw {raw:>6.0f}  CP10k {cp10k:>8.2f}  "
                      f"log1p {lg:>7.4f}  (gene mean {ln[r]['mean']:.4f}, "
                      f"sd {ln[r]['sd']:.4f})  z {z:>+7.4f}")
            n_more = len(panel_rows[mod]) - SHOW_GENES
            if n_more > 0:
                print(f"    ... and {n_more} more panel genes")
            pan = float(np.mean([ln[r]["z"][i] for r in panel_rows[mod]]))
            ctl = float(np.mean([ln[r]["z"][i] for r in ctrl_rows[mod]]))
            print(f"    mean z over {len(panel_rows[mod])} panel genes      "
                  f"{pan:>+8.4f}")
            print(f"    mean z over {len(ctrl_rows[mod])} control genes    "
                  f"{ctl:>+8.4f}")
            print(f"    module score = panel - control          {pan - ctl:>+8.4f}")
            rec[f"{mod}_panel_mean_z"] = round(pan, 4)
            rec[f"{mod}_control_mean_z"] = round(ctl, 4)
            rec[f"{mod}_score"] = round(pan - ctl, 4)
        axis = rec["HYPER_score"] - rec["DEPOL_score"]
        rec["axis_hand"] = round(axis, 4)
        rec["axis_pipeline"] = round(float(scores["HminusD"][i]), 4)
        rec["agree"] = abs(rec["axis_hand"] - rec["axis_pipeline"]) < 1e-3
        print(f"  AXIS = HYPER - DEPOL = {rec['HYPER_score']:+.4f} - "
              f"{rec['DEPOL_score']:+.4f} = {axis:+.4f}")
        print(f"  pipeline value for this nucleus: {rec['axis_pipeline']:+.4f}  "
              f"{'(agrees)' if rec['agree'] else '(DISAGREES)'}")
        out["nuclei"].append(rec)

    print("\n" + "-" * 78)
    print("Units: each gene is centred and scaled across nuclei before averaging, so a")
    print("module score is in per-gene standard deviations, not in counts or log-fold")
    print("units. An axis value of +0.25 means the hyperpolarising contrast exceeds the")
    print("depolarising contrast by a quarter of a per-gene standard deviation. No fold")
    print("or percentage reading of these values is available.")
    print("-" * 78)

    out["all_agree"] = all(r["agree"] for r in out["nuclei"])
    out["units_note"] = ("per-gene z-scores averaged within panel and control sets; "
                         "module score is panel minus control; axis is HYPER minus DEPOL; "
                         "units are per-gene standard deviations across nuclei")
    json.dump(out, open(os.path.join(cfg.RESULTS, "worked_example.json"), "w"), indent=2)
    print(f"\nall nuclei agree with the pipeline: {out['all_agree']}")
    print("Written to results/worked_example.json")


if __name__ == "__main__":
    main()
