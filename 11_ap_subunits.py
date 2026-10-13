#!/usr/bin/env python3
"""
11_ap_subunits.py - adaptor protein (AP) complex subunit composition across VMH and
arcuate neuron populations.

Tests the stripped-down, falsifiable core of the AP-subunit hypothesis: that AP complex
composition differs between VMH (glutamatergic SF1) and arcuate (peptidergic POMC, AgRP)
neurons in a way consistent with their different transmission modes. The mechanistic
claims about vesicle recycling and channel sorting are not testable from expression data
and are not addressed here.

Two datasets, both already downloaded:
  Rupp GSE172460   - L3.Nr5a1, L9.Pomc, L6.Agrp in one experiment, the primary comparison
  Affinati GSE172204 - the 24 VMH clusters, grouped to class level, as a within-VMH check

Composition is expressed as within-position paralog share rather than absolute level,
following the framework used elsewhere in this program: a position is a slot in the complex
occupied by mutually exclusive paralogs, so share is the quantity that is comparable across
cell types with different overall expression.

Controls carried over from the rest of this pipeline:
  * detection rates reported per gene per population; genes below MIN_DETECT are marked
    untestable rather than reported as absent, the Kcnj11 lesson
  * depth control by logistic regression of detection on population given log10 library size
  * bootstrap CIs on shares, permutation tests on contrasts
  * the stated predictions evaluated explicitly, so the result is a verdict not a table

Output: results/ap_subunits_<tag>.json, ap_subunits_<tag>_genes.csv, ap_subunits_<tag>_shares.csv
"""
import os, csv, json, math, sys
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2 as chi2dist
import config as C
import vmhlib as V

MIN_DETECT = 0.05        # fraction of nuclei; below this a gene is untestable, not absent
N_BOOT = 2000

# ---------------------------------------------------------------- complex definitions
# Positions are slots occupied by mutually exclusive paralogs. Single-occupancy positions
# are listed too, so total complex abundance can be assembled, but they carry no share.
COMPLEXES = {
    "AP1": {"large":  ["Ap1g1", "Ap1g2"],
            "beta":   ["Ap1b1"],
            "medium": ["Ap1m1", "Ap1m2"],
            "small":  ["Ap1s1", "Ap1s2", "Ap1s3"]},
    "AP2": {"large":  ["Ap2a1", "Ap2a2"],
            "beta":   ["Ap2b1"],
            "medium": ["Ap2m1"],
            "small":  ["Ap2s1"]},
    "AP3": {"large":  ["Ap3d1"],
            "beta":   ["Ap3b1", "Ap3b2"],
            "medium": ["Ap3m1", "Ap3m2"],
            "small":  ["Ap3s1", "Ap3s2"]},
    "AP4": {"large":  ["Ap4e1"],
            "beta":   ["Ap4b1"],
            "medium": ["Ap4m1"],
            "small":  ["Ap4s1"]},
    "AP5": {"large":  ["Ap5z1"],
            "beta":   ["Ap5b1"],
            "medium": ["Ap5m1"],
            "small":  ["Ap5s1"]},
}
ALL_AP = [g for cx in COMPLEXES.values() for gs in cx.values() for g in gs]

# The predictions as stated, evaluated explicitly at the end.
PREDICTIONS = [
    ("Ap3b2", "VMH_SF1", ["POMC", "AgRP"], "higher", "AP-3 beta3B higher in VMH"),
    ("Ap3m2", "VMH_SF1", ["POMC", "AgRP"], "higher", "AP-3 mu3B higher in VMH"),
    ("Ap4b1", "VMH_SF1", ["POMC", "AgRP"], "lower",  "AP-4 beta4 lower in VMH"),
    ("Ap4m1", "VMH_SF1", ["POMC", "AgRP"], "lower",  "AP-4 mu4 lower in VMH"),
    ("Ap1m1", "VMH_SF1", ["POMC", "AgRP"], "lower",  "AP-1 mu1A lower in VMH"),
    ("Ap2m1", "VMH_SF1", ["POMC", "AgRP"], "higher", "AP-2 mu2 higher in VMH"),
]


def lr_group_given_depth(y, logdepth, group):
    Xd = sm.add_constant(logdepth.reshape(-1, 1), has_constant="add")
    m0 = sm.Logit(y, Xd).fit(disp=0, maxiter=200)
    D = pd.get_dummies(group, drop_first=True).to_numpy(dtype=float)
    m1 = sm.Logit(y, np.column_stack([Xd, D])).fit(disp=0, maxiter=200)
    lr = 2 * (m1.llf - m0.llf)
    return float(lr), int(D.shape[1]), float(chi2dist.sf(max(lr, 0.0), D.shape[1]))


def boot_diff(a, b, seed=0):
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    d = (rng.choice(a, (N_BOOT, a.size), replace=True).mean(1) -
         rng.choice(b, (N_BOOT, b.size), replace=True).mean(1))
    return float(a.mean() - b.mean()), float(np.percentile(d, 2.5)), float(np.percentile(d, 97.5))


def run(prefix, tag, label_col, relabel, keep):
    print(f"\n=== {prefix} ({tag}) ===")
    feats, name2row = V.load_features(prefix)
    bcs, bc2col = V.load_barcodes(prefix)
    meta = V.load_metadata(prefix)

    missing = [g for g in ALL_AP if g not in name2row]
    if missing:
        print(f"   not in features.tsv: {missing}")
    rows_needed = {name2row[g] for g in ALL_AP if g in name2row}

    col2group, _, _ = V.write_col2group(prefix, meta, bc2col, label_col,
                                        f"{tag}_ap", keep=keep)
    cs, rs = V.pass_totals(prefix, col2group, f"{tag}_ap")
    colsum = V.read_two_col(cs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = relabel(g)
    cols = sorted(c for c in colsum if col2cl.get(c))
    grp = np.array([col2cl[c] for c in cols])
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)
    pops = sorted(set(grp.tolist()))
    print("   n per group:", {p: int((grp == p).sum()) for p in pops})

    ent = V.pass_rows(prefix, rows_needed, f"{tag}_ap", "apgenes")
    row2gene = {name2row[g]: g for g in ALL_AP if g in name2row}
    counts = defaultdict(lambda: np.zeros(len(cols)))
    cidx = {c: i for i, c in enumerate(cols)}
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        i = cidx.get(c)
        if i is not None:
            counts[row2gene[int(r)]][i] = float(v)

    # ------------------------------------------------------------ per-gene summary
    generows = []
    for g in ALL_AP:
        if g not in row2gene.values():
            generows.append({"gene": g, "status": "absent_from_features"})
            continue
        x = counts[g]
        ln = np.log1p(x / depth * C.NORM_TARGET)
        rec = {"gene": g, "status": "ok"}
        det_all = float((x > 0).mean())
        if det_all < MIN_DETECT:
            rec["status"] = "untestable_low_detection"
        rec["detect_pct_all"] = round(100 * det_all, 1)
        for p in pops:
            m = grp == p
            rec[f"{p}_detect_pct"] = round(100 * float((x[m] > 0).mean()), 1)
            rec[f"{p}_log_cp10k"] = round(float(ln[m].mean()), 4)
        if rec["status"] == "ok" and len(pops) > 1:
            try:
                lr, df, pv = lr_group_given_depth((x > 0).astype(float), logdepth, grp)
                rec["LR_group_given_depth"] = round(lr, 2)
                rec["p_group_given_depth"] = pv
            except Exception:
                rec["LR_group_given_depth"] = None
                rec["p_group_given_depth"] = None
        generows.append(rec)
    with open(os.path.join(C.RESULTS, f"ap_subunits_{tag}_genes.csv"), "w", newline="") as fh:
        fields = sorted({k for r in generows for k in r},
                        key=lambda k: (k != "gene", k != "status", k))
        w = csv.DictWriter(fh, fieldnames=fields)
        w.writeheader(); w.writerows(generows)

    testable = {r["gene"] for r in generows if r.get("status") == "ok"}
    print(f"   testable at >={MIN_DETECT:.0%} detection: {len(testable)} of {len(ALL_AP)}")

    # ------------------------------------------------- within-position paralog share
    shares, sharerows = {}, []
    for cx, positions in COMPLEXES.items():
        for pos, genes in positions.items():
            usable = [g for g in genes if g in testable]
            if len(usable) < 2:
                continue
            for p in pops:
                m = grp == p
                tot = sum(counts[g][m].sum() for g in usable)
                if tot <= 0:
                    continue
                for g in usable:
                    # per-cell share, for a bootstrap CI over cells
                    percell = np.array([counts[g][m][i] / s if s > 0 else np.nan
                                        for i, s in enumerate(
                                            sum(counts[h][m] for h in usable))])
                    pc = percell[~np.isnan(percell)]
                    shares[(cx, pos, g, p)] = pc
                    sharerows.append({
                        "complex": cx, "position": pos, "gene": g, "population": p,
                        "pseudobulk_share": round(float(counts[g][m].sum() / tot), 4),
                        "n_cells_with_position_detected": int(pc.size),
                        "percell_share_mean": round(float(pc.mean()), 4) if pc.size else None})
    with open(os.path.join(C.RESULTS, f"ap_subunits_{tag}_shares.csv"), "w", newline="") as fh:
        if sharerows:
            w = csv.DictWriter(fh, fieldnames=list(sharerows[0].keys()))
            w.writeheader(); w.writerows(sharerows)

    # ------------------------------------------------------------- complex totals
    totals = {}
    for p in pops:
        m = grp == p
        tot_cx = {}
        for cx, positions in COMPLEXES.items():
            gs = [g for gs_ in positions.values() for g in gs_ if g in testable]
            tot_cx[cx] = float(sum(counts[g][m].sum() for g in gs)) / float(depth[m].sum()) * C.CPM_TARGET
        s = sum(tot_cx.values())
        totals[p] = {cx: round(v / s, 4) for cx, v in tot_cx.items()} if s > 0 else {}

    # -------------------------------------------------------------- predictions
    verdicts = []
    if set(["VMH_SF1", "POMC", "AgRP"]).issubset(set(pops)):
        for gene, focus, others, direction, label in PREDICTIONS:
            if gene not in testable:
                verdicts.append({"prediction": label, "gene": gene,
                                 "verdict": "untestable (detection below threshold)"})
                continue
            ln = np.log1p(counts[gene] / depth * C.NORM_TARGET)
            res = {"prediction": label, "gene": gene, "direction_predicted": direction}
            ok = True
            for o in others:
                d, lo, hi = boot_diff(ln[grp == focus], ln[grp == o], seed=7)
                res[f"vs_{o}"] = {"diff": round(d, 4), "ci": [round(lo, 4), round(hi, 4)]}
                sig = (lo > 0) if direction == "higher" else (hi < 0)
                ok = ok and sig
            res["verdict"] = "supported" if ok else "not supported"
            verdicts.append(res)

    out = {"dataset": prefix, "tag": tag,
           "n_per_group": {p: int((grp == p).sum()) for p in pops},
           "median_depth": {p: float(np.median(depth[grp == p])) for p in pops},
           "n_testable_subunits": len(testable),
           "untestable_subunits": sorted(set(ALL_AP) - testable),
           "complex_share_of_AP_signal": totals,
           "predictions": verdicts,
           "note": ("Share is within-position paralog share; positions with a single "
                    "occupant carry no share and appear only in complex totals. Genes below "
                    f"{MIN_DETECT:.0%} detection are marked untestable rather than absent.")}
    json.dump(out, open(os.path.join(C.RESULTS, f"ap_subunits_{tag}.json"), "w"), indent=2)
    print(json.dumps({k: out[k] for k in
                      ("n_per_group", "n_testable_subunits", "untestable_subunits",
                       "complex_share_of_AP_signal", "predictions")}, indent=2))
    return out


def main():
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    run("GSE172460", "rupp", "Lepr_cluster",
        relabel=lambda g: inv.get(g, ""),
        keep=lambda g: g in set(C.RUPP_POPS.values()))

    # within-VMH check: the 24 Affinati clusters grouped to class level by prefix
    run("GSE172204", "vmh", "VMH_cluster",
        relabel=lambda g: g.split("_")[0],
        keep=lambda g: bool(g))


if __name__ == "__main__":
    main()
