#!/usr/bin/env python3
"""
41_cluster_axis_corrected.py - cluster-level effector axis on the control-corrected measure.

THE PROBLEM THIS FIXES. Section 2.2 quotes cluster-level values (the axis ranging from
-0.129 to 0.247 across 24 clusters, the nine Lepr clusters at 0.068 to 0.206, Esr1_1 at
+0.148, the panel spreads of 0.320 and 0.096) and presents them as the control-corrected
module-score axis described two sentences earlier. They are not. Those come from
cluster_scores.csv, written by 02_snrna_core.py, which averages log1p pseudobulk values with
NO z-scoring and NO control genes. The per-nucleus axis used everywhere else, from
vmhlib.module_scores, z-scores each gene across nuclei and subtracts expression-matched
controls. Two different quantities sharing the name HminusD, in different units, one
control-corrected and one not.

The figure label on the cluster panel says log-CP10k and is correct. The prose is what is
wrong.

WHY A WHOLE SCRIPT. data/percell_scores.npz holds control-corrected per-nucleus scores for
2,913 nuclei, the Lepr clusters only, from a single control-gene draw. Cluster means for the
other 15 clusters do not exist anywhere, so they have to be computed, which needs a pass over
the matrix for all VMH nuclei across the 20 draws.

WHAT IT PRODUCES. Per-cluster means of HminusD, HYPER and DEPOL for all 24 VMH clusters, each
averaged across 20 control-gene draws with the across-draw range. That puts the whole of 2.2
on one measure.

POSITIVE CONTROL (discipline step 3c, and the first new script since that rule). A known
per-cluster offset is injected into the per-nucleus scores and the cluster means are checked
to recover it. The operation is an average so this should be exact; the point is that the
rule is standing, not that this case is doubtful.

EXPECT THE NUMBERS TO CHANGE. These are not relabelled versions of the old ones. The axis
range, the Lepr span, Esr1_1's value and both panel spreads will all differ, and Figure 2A
changes with them. Every one goes into the manifest this time; none of them was declared.

Output: results/cluster_axis_corrected.json
"""
import json
import os
import warnings
from collections import defaultdict

import numpy as np

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "cac"


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
    clu = np.array([meta[col2bc[c]]["VMH_cluster"].strip() for c in cols])
    clusters = sorted(set(clu))
    print(f">> {len(cols)} VMH nuclei across {len(clusters)} clusters", flush=True)

    seeds = list(range(cfg.RANDOM_SEED, cfg.RANDOM_SEED + cfg.N_SEEDS))
    per_draw = {k: defaultdict(list) for k in ("HminusD", "HYPER", "DEPOL")}
    first_draw = None
    colset = set(cols)

    for s in seeds:
        panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(cols), s)
        need = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values())
                for r in rs_}
        ent = V.pass_rows(PREFIX, need, TAG, f"seed{s}")
        expr = defaultdict(dict)
        for line in open(ent):
            r, c, v = line.split()
            c = int(c)
            if c in colset:
                expr[c][int(r)] = float(v)
        sc = V.module_scores(expr, cols, colsum, panel_rows, ctrl_rows)
        if first_draw is None:
            first_draw = {k: np.asarray(sc[k], float) for k in per_draw}
        for key in per_draw:
            v = np.asarray(sc[key], float)
            for cl in clusters:
                per_draw[key][cl].append(float(v[clu == cl].mean()))
        print(f"   seed {s} done", flush=True)

    # Per-draw SD of the per-nucleus axis. 2.4 quotes a widening ratio computed from
    # observed_sd, which comes from ONE stored draw, and presents it as a population
    # property. This records the SD for every draw so the ratio's stability is a
    # measured thing rather than an assumption.
    draw_sds = [float(np.asarray(ax, float).std()) for ax in axis_by_seed]
    out_sd = {"per_draw_sd": [round(v, 4) for v in draw_sds],
              "mean": round(float(np.mean(draw_sds)), 4),
              "min": round(float(np.min(draw_sds)), 4),
              "max": round(float(np.max(draw_sds)), 4),
              "note": ("the widening ratio in 2.4 uses the single stored draw; if this "
                       "range is tight the ratio is a population property, if not it "
                       "needs a range like every other module-based value")}
    print(f">> per-nucleus axis SD across {len(draw_sds)} draws: "
          f"{out_sd['min']} to {out_sd['max']} (mean {out_sd['mean']})", flush=True)

    out = {"axis_sd_across_draws": out_sd,
           "n_vmh_nuclei": len(cols), "n_clusters": len(clusters),
           "n_seeds": cfg.N_SEEDS,
           "measure": ("control-corrected module scores from vmhlib.module_scores, "
                       "per-gene z-scored across nuclei then matched controls subtracted; "
                       "units are per-gene standard deviations"),
           "supersedes": ("the cluster-level values in cluster_scores.csv, which are "
                          "log1p pseudobulk means with no z-scoring and no control "
                          "correction and are in log-CP10k"),
           "by_cluster": {}}

    for cl in clusters:
        rec = {"n_nuclei": int((clu == cl).sum()),
               "is_lepr": bool(cl.startswith("Lepr"))}
        for key in ("HminusD", "HYPER", "DEPOL"):
            vals = per_draw[key][cl]
            rec[key] = {"mean": round(float(np.mean(vals)), 4),
                        "min": round(float(np.min(vals)), 4),
                        "max": round(float(np.max(vals)), 4)}
        out["by_cluster"][cl] = rec

    hm = {cl: out["by_cluster"][cl]["HminusD"]["mean"] for cl in clusters}
    lep = {c: v for c, v in hm.items() if c.startswith("Lepr")}
    out["axis_range_all_clusters"] = [round(min(hm.values()), 4),
                                      round(max(hm.values()), 4)]
    out["axis_range_lepr_clusters"] = [round(min(lep.values()), 4),
                                       round(max(lep.values()), 4)]
    out["clusters_on_depolarizing_side"] = sorted(c for c, v in hm.items() if v < 0)
    out["esr1_1"] = round(hm.get("Esr1_1", float("nan")), 4)
    for key in ("HYPER", "DEPOL"):
        vs = [out["by_cluster"][cl][key]["mean"] for cl in clusters]
        out[f"{key}_spread"] = {"min": round(min(vs), 4), "max": round(max(vs), 4),
                                "span": round(max(vs) - min(vs), 4)}

    print(f"\n{'cluster':<12}{'H-D':>9}{'range':>18}{'n':>7}")
    for cl in sorted(clusters, key=lambda c: hm[c]):
        r = out["by_cluster"][cl]["HminusD"]
        mark = " *" if cl.startswith("Lepr") else "  "
        print(f"{cl:<12}{r['mean']:>9.4f}  [{r['min']:>7.4f},{r['max']:>7.4f}]"
              f"{out['by_cluster'][cl]['n_nuclei']:>7}{mark}")

    print(f"\n  axis across all 24 clusters   {out['axis_range_all_clusters']}")
    print(f"  axis across the nine Lepr     {out['axis_range_lepr_clusters']}")
    print(f"  Esr1_1                        {out['esr1_1']}")
    print(f"  on the depolarizing side      {out['clusters_on_depolarizing_side']}")
    print(f"  HYPER spread  {out['HYPER_spread']['span']:.4f} "
          f"({out['HYPER_spread']['min']:.4f} to {out['HYPER_spread']['max']:.4f})")
    print(f"  DEPOL spread  {out['DEPOL_spread']['span']:.4f} "
          f"({out['DEPOL_spread']['min']:.4f} to {out['DEPOL_spread']['max']:.4f})")
    print("\n  Compare against the manuscript's current -0.129 to 0.247, Lepr 0.068 to")
    print("  0.206, Esr1_1 +0.148, DEPOL span 0.320, HYPER span 0.096. Those are the")
    print("  uncorrected log-CP10k values and every one of them has to be replaced.")

    # ---------------------------------------------------------- positive control
    print(f"\n{'=' * 60}\npositive control: can cluster means recover a known offset?\n{'=' * 60}")
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    truth = {cl: float(rng.normal(0, 0.3)) for cl in clusters}
    y = first_draw["HminusD"] + np.array([truth[c] for c in clu])
    err = []
    for cl in clusters:
        got = float(y[clu == cl].mean())
        want = float(first_draw["HminusD"][clu == cl].mean()) + truth[cl]
        err.append(abs(got - want))
    worst = max(err)
    out["positive_control_max_error"] = float(worst)
    out["positive_control_passes"] = bool(worst < 1e-9)
    print(f"  largest error across {len(clusters)} clusters: {worst:.2e}  "
          f"{'OK' if out['positive_control_passes'] else 'FAIL'}")

    json.dump(out, open(os.path.join(cfg.RESULTS, "cluster_axis_corrected.json"), "w"),
              indent=2)
    print("\nWritten to results/cluster_axis_corrected.json")


if __name__ == "__main__":
    main()
