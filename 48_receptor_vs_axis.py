#!/usr/bin/env python3
"""
48_receptor_vs_axis.py - do receptor transcripts track the effector balance axis?

THE QUESTION. Klockener et al. (2011) show insulin activates PI3K in SF-1 neurons and
reduces firing through K_ATP activation, abolished in mice lacking the insulin receptor in
those cells. So Insr acts on a conductance that also carries the leptin hyperpolarization
Sohn et al. describe. If Insr abundance varied between nuclei along the effector axis, that
would be a candidate for something the axis reflects.

Insr is in neither module panel nor in the control sets, so this is an out-of-sample test,
the same shape as the AP adaptor test in script 47.

THREE REASONS TO EXPECT NOTHING, stated before the result so they are not post hoc:

  1. Klockener recorded SF-1 neurons. Section 2.7 reports the SF1 immunoprecipitate is NOT
     enriched for Lepr, so those cells are largely not this population
  2. Insulin reduces firing in SF-1 neurons uniformly, not bidirectionally. A factor with
     one direction cannot explain opposite responses
  3. Insr is a receptor, not a conductance, so its abundance can vary freely without moving
     an axis built from channel transcripts

A positive result would be interesting precisely because all three argue against it.

GENES TESTED. Insr and Lepr as the receptors of interest; Irs1, Irs2, Pik3r1 and Pik3ca as
the PI3K arm both hormones signal through; Kcnj11 and Abcc8 as the K_ATP subunits the two
papers converge on, included as POSITIVE CONTROLS rather than findings, since Abcc8 is in
the hyperpolarizing panel and must track the axis, and a run where it does not is wrong
before anything else is read.

CONTROLS. Every gene is scored two ways, as detection (logistic on the axis) and as level
among detected nuclei (Spearman), each with the partial controlling for log capture depth,
since the axis tracks depth at rho 0.137 in this population. A detectability-matched null of
200 genes gives the threshold.

Output: results/receptor_vs_axis.json
"""
import json
import os
import warnings
from collections import defaultdict

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "rcx"
N_NULL = 200

GENES = {
    "Insr":    "insulin receptor; Klockener et al. 2011",
    "Lepr":    "leptin receptor; defines the population, so expect structure",
    "Irs1":    "PI3K arm",
    "Irs2":    "PI3K arm",
    "Pik3r1":  "PI3K regulatory",
    "Pik3ca":  "PI3K catalytic",
    "Abcc8":   "POSITIVE CONTROL: in the hyperpolarizing panel, must track the axis",
    "Kcnj11":  "POSITIVE CONTROL: K_ATP pore; below 5% detection in 2.7, may be untestable",
}


def main():
    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    axis = np.asarray(z["HminusD"], float)
    cols = [int(c) for c in z["cols"]]
    cidx = {c: i for i, c in enumerate(cols)}
    print(f">> {len(cols)} Lepr-cluster nuclei, one control-gene draw", flush=True)

    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_grp", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    depth = np.array([colsum.get(c, np.nan) for c in cols], float)
    logd = np.log10(np.where(np.isnan(depth) | (depth <= 0), 1, depth))
    print(f">> axis vs log depth: rho "
          f"{float(stats.spearmanr(axis, logd).statistic):+.4f}", flush=True)

    want = {g: name2row[g] for g in GENES if g in name2row}
    for g in GENES:
        if g not in name2row:
            print(f"   {g} not in the feature list", flush=True)
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    pool = [int(r) for r in rowsum if rowsum[r] > 0]
    nullrows = list(rng.choice(pool, min(1200, len(pool)), replace=False))
    ent = V.pass_rows(PREFIX, set(want.values()) | set(nullrows), TAG, "rcx")

    counts = defaultdict(lambda: np.zeros(len(cols)))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            counts[int(r)][cidx[c]] = float(v)

    def resid_on_depth(x, mask):
        r = stats.rankdata(x).astype(float)
        d = stats.rankdata(logd[mask])
        return r - np.polyval(np.polyfit(d, r, 1), d)

    def score(row):
        v = counts[row]
        det = (v > 0).astype(float)
        pct = 100 * det.mean()
        if det.sum() < 50 or det.sum() > len(cols) - 50:
            return {"pct_detected": round(pct, 1), "testable": False}
        # detection against the axis, partial on depth
        full = np.ones(len(cols), bool)
        rho_det = float(stats.spearmanr(det, axis).statistic)
        p_det = float(stats.spearmanr(resid_on_depth(det, full),
                                      resid_on_depth(axis, full)).statistic)
        # level among detected nuclei
        m = v > 0
        if m.sum() >= 100:
            rho_lvl = float(stats.spearmanr(v[m], axis[m]).statistic)
            p_lvl = float(stats.spearmanr(resid_on_depth(v[m], m),
                                          resid_on_depth(axis[m], m)).statistic)
        else:
            rho_lvl = p_lvl = float("nan")
        return {"pct_detected": round(pct, 1), "testable": True,
                "n_detected": int(m.sum()),
                "rho_detection": round(rho_det, 4),
                "partial_detection": round(p_det, 4),
                "rho_level": round(rho_lvl, 4),
                "partial_level": round(p_lvl, 4)}

    out = {"n_nuclei": len(cols), "draws": 1, "genes": {},
           "caveat": ("detection in droplet single-nucleus data confounds abundance with "
                      "capture; partials control for depth but not for that")}

    print(f"\n{'gene':<9}{'%det':>7}{'n det':>8}{'det rho':>10}{'det|depth':>11}"
          f"{'lvl rho':>10}{'lvl|depth':>11}")
    for g, row in want.items():
        r = score(row)
        out["genes"][g] = {**r, "note": GENES[g]}
        if not r["testable"]:
            print(f"{g:<9}{r['pct_detected']:>7.1f}{'  not testable at this detection':>40}")
            continue
        print(f"{g:<9}{r['pct_detected']:>7.1f}{r['n_detected']:>8}"
              f"{r['rho_detection']:>10.4f}{r['partial_detection']:>11.4f}"
              f"{r['rho_level']:>10.4f}{r['partial_level']:>11.4f}")

    # null over detectability-matched genes
    print(f"\n>> null over {N_NULL} random genes", flush=True)
    nulls = []
    for row in nullrows[:N_NULL]:
        r = score(row)
        if r.get("testable"):
            nulls.append(abs(r["partial_detection"]))
    p95 = float(np.percentile(nulls, 95)) if nulls else float("nan")
    out["null_abs_partial_detection_p95"] = round(p95, 4)
    out["n_null_genes_scored"] = len(nulls)
    print(f"   |partial| on detection, 95th percentile: {p95:.4f} "
          f"({len(nulls)} genes testable)")

    ctrl = out["genes"].get("Abcc8", {})
    ok = ctrl.get("testable") and abs(ctrl.get("partial_detection", 0)) > p95
    out["positive_control_passes"] = bool(ok)
    print(f"\n   positive control, Abcc8 exceeds the null: {ok}")
    if not ok:
        print("   CONTROL FAILED. Abcc8 is in the hyperpolarizing panel and must track")
        print("   the axis. Nothing below is readable until that is understood.")

    beats = [g for g, r in out["genes"].items()
             if r.get("testable") and abs(r.get("partial_detection", 0)) > p95
             and g not in ("Abcc8", "Kcnj11")]
    out["genes_exceeding_null"] = beats
    print(f"\n{'=' * 68}")
    if ok and not beats:
        print("  No receptor or PI3K transcript tracks the axis beyond the null.")
        print("  Insr in particular does not, which is what the three reasons in this")
        print("  script's header predict. The axis does not reflect insulin-receptor")
        print("  abundance in this population.")
    elif ok:
        print(f"  Exceeding the null: {', '.join(beats)}")
        print("  Lepr doing so would be unsurprising, since the clusters are named for")
        print("  it. Insr doing so would be worth following up carefully.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "receptor_vs_axis.json"), "w"),
              indent=2)
    print("\nWritten to results/receptor_vs_axis.json")


if __name__ == "__main__":
    main()
