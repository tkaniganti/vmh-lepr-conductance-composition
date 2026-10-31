#!/usr/bin/env python3
"""
50_camp_pka_vs_axis.py - does the cAMP-PKA module track the effector balance axis?

THE PREDICTION, STATED BEFORE THE RUN. This is expected to be NULL. Section 2.5 scanned
12,182 transcripts individually against this axis and found nothing above about 2% of the
variance. A module score is largely an average of its members' correlations, so for a set to
track the axis while no member does, the members would have to covary with it weakly AND in
a coordinated direction. The scan makes that unlikely. This is written here so that a null
result reads as a prediction borne out rather than a rationalisation afterwards, and so that
a positive result has to survive the scepticism it would deserve.

WHY TEST IT ANYWAY. cAMP-PKA is the one pathway with existing evidence of non-uniformity in
this dataset: 2.7 reports Pde3b depleted in the SF1 immunoprecipitate at -0.82, q = 0.026.
It is also the post-translational mechanism the Discussion names as invisible to steady-state
mRNA, and the zebrafish case the Discussion cites works by phosphorylation converging on a
single population. If any module were going to carry information the individual scan missed,
this is the one with a reason behind it.

A NOTE ON WHAT A POSITIVE WOULD AND WOULD NOT MEAN. PKA acts by phosphorylation. Transcript
abundance of its subunits says nothing about kinase activity, and the Discussion's claim is
precisely that the determinant is a post-translational state. So a module correlation here
would be consistent with the mechanism without being evidence for it, and the manuscript
must not report it as support.

PDE3B IS EXCLUDED FROM THE MODULE. It sits in the hyperpolarizing panel, so including it
would correlate the axis partly with itself. Its own correlation is reported separately and
labelled circular.

CONTROLS, both gating:
  depth        the axis tracks capture depth at rho 0.137 here, and module scores built from
               detection-prone genes inherit that. Only the depth-partial is a result
  random sets  200 detectability-matched gene sets of the same size, scored identically.
               With 2,913 nuclei a partial of 0.04 means nothing without this

Output: results/camp_pka_vs_axis.json
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
PREFIX, TAG = "GSE172204", "cpk"
N_NULL = 200
MIN_DETECT = 0.05

# cAMP production, degradation, and the kinase. Pde3b excluded: it is in the
# hyperpolarizing panel and would correlate the axis with itself.
MODULE = {
    "adenylyl cyclase": ["Adcy1", "Adcy2", "Adcy3", "Adcy5", "Adcy6", "Adcy8", "Adcy9"],
    "phosphodiesterase": ["Pde4a", "Pde4b", "Pde4d", "Pde1a", "Pde1b", "Pde10a"],
    "PKA": ["Prkaca", "Prkacb", "Prkar1a", "Prkar1b", "Prkar2a", "Prkar2b"],
    "effectors": ["Creb1", "Rapgef3", "Rapgef4"],
}
CIRCULAR = ["Pde3b"]


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
    depth = np.array([float(colsum.get(c, 0)) for c in cols])
    logd = np.log10(np.where(depth <= 0, 1, depth))
    print(f">> axis vs log depth: rho "
          f"{float(stats.spearmanr(axis, logd).statistic):+.4f}", flush=True)

    flat = [g for gs in MODULE.values() for g in gs] + CIRCULAR
    present = {g: name2row[g] for g in flat if g in name2row}
    absent = [g for g in flat if g not in name2row]
    if absent:
        print(f">> not in the feature list: {', '.join(absent)}", flush=True)

    rng = np.random.default_rng(cfg.RANDOM_SEED)
    pool = [int(r) for r in rowsum if rowsum[r] > 0]
    nullpool = list(rng.choice(pool, min(2500, len(pool)), replace=False))
    ent = V.pass_rows(PREFIX, set(present.values()) | set(nullpool), TAG, "cpk")

    counts = defaultdict(lambda: np.zeros(len(cols)))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            counts[int(r)][cidx[c]] = float(v)

    def zscore(row):
        v = np.log1p(counts[row] / np.where(depth <= 0, 1, depth) * cfg.NORM_TARGET)
        sd = v.std()
        return (v - v.mean()) / (sd if sd > 0 else 1.0)

    def partial(x):
        rx = stats.rankdata(x).astype(float)
        ra = stats.rankdata(axis).astype(float)
        rd = stats.rankdata(logd).astype(float)
        fx = rx - np.polyval(np.polyfit(rd, rx, 1), rd)
        fa = ra - np.polyval(np.polyfit(rd, ra, 1), rd)
        return float(stats.spearmanr(x, axis).statistic), \
            float(stats.spearmanr(fx, fa).statistic)

    out = {"n_nuclei": len(cols), "draws": 1, "min_detect": MIN_DETECT,
           "prediction": ("null, because the 2.5 scan found no single transcript above "
                          "about 2% of the variance and a module score is largely an "
                          "average of its members"),
           "caveat": ("PKA acts by phosphorylation; subunit transcript abundance says "
                      "nothing about kinase activity, so a correlation here would be "
                      "consistent with the mechanism without being evidence for it"),
           "genes": {}, "submodules": {}}

    # per gene, so a module result can be traced to its members
    usable = {}
    print(f"\n{'gene':<10}{'%det':>7}{'rho raw':>10}{'rho|depth':>11}")
    for g, row in sorted(present.items()):
        pct = 100 * (counts[row] > 0).mean()
        if pct < MIN_DETECT * 100:
            out["genes"][g] = {"pct_detected": round(pct, 1), "testable": False}
            print(f"{g:<10}{pct:>7.1f}{'  below detection floor':>28}")
            continue
        raw, par = partial(zscore(row))
        out["genes"][g] = {"pct_detected": round(pct, 1), "testable": True,
                           "rho_raw": round(raw, 4), "rho_partial": round(par, 4),
                           "circular": g in CIRCULAR}
        if g not in CIRCULAR:
            usable[g] = row
        tag = "  (circular, in the panel)" if g in CIRCULAR else ""
        print(f"{g:<10}{pct:>7.1f}{raw:>10.4f}{par:>11.4f}{tag}")

    # submodules and the whole module
    print(f"\n{'module':<22}{'genes':>7}{'rho raw':>10}{'rho|depth':>11}")
    for name, gs in list(MODULE.items()) + [("ALL cAMP-PKA", [g for gs in MODULE.values()
                                                              for g in gs])]:
        rows = [usable[g] for g in gs if g in usable]
        if len(rows) < 2:
            print(f"{name:<22}{len(rows):>7}{'  too few usable genes':>28}")
            out["submodules"][name] = {"n_genes": len(rows), "testable": False}
            continue
        sc = np.mean([zscore(r) for r in rows], axis=0)
        raw, par = partial(sc)
        out["submodules"][name] = {"n_genes": len(rows), "testable": True,
                                   "rho_raw": round(raw, 4),
                                   "rho_partial": round(par, 4)}
        print(f"{name:<22}{len(rows):>7}{raw:>10.4f}{par:>11.4f}")

    # null: detectability-matched sets of the same size as the full module
    k = out["submodules"].get("ALL cAMP-PKA", {}).get("n_genes", 0)
    print(f"\n>> null over {N_NULL} random sets of {k} genes", flush=True)
    nulls = []
    cand = [r for r in nullpool if (counts[r] > 0).mean() >= MIN_DETECT]
    for _ in range(N_NULL):
        if len(cand) < k or k < 2:
            break
        pick = list(rng.choice(cand, k, replace=False))
        sc = np.mean([zscore(r) for r in pick], axis=0)
        nulls.append(abs(partial(sc)[1]))
    p95 = float(np.percentile(nulls, 95)) if nulls else float("nan")
    out["null_abs_partial_p95"] = round(p95, 4)
    out["n_null_sets"] = len(nulls)
    print(f"   |partial| 95th percentile: {p95:.4f} ({len(nulls)} sets, "
          f"{len(cand)} genes above the floor)")

    beats = [n for n, v in out["submodules"].items()
             if v.get("testable") and abs(v["rho_partial"]) > p95]
    out["submodules_exceeding_null"] = beats
    print(f"\n{'=' * 68}")
    if not beats:
        print("  The cAMP-PKA module does not track the axis. Predicted, for the reason")
        print("  in this script's header: the 2.5 scan already excluded single-transcript")
        print("  correlates, and a module average inherits that.")
        print("  Seventh negative. Nothing measured predicts position on this axis.")
    else:
        print(f"  Exceeding the null: {', '.join(beats)}")
        print("  Against the prediction stated in the header, so treat it sceptically:")
        print("  check which members carry it, whether they are detection-prone, and")
        print("  whether the result survives a second control-gene draw of the axis.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "camp_pka_vs_axis.json"), "w"),
              indent=2)
    print("\nWritten to results/camp_pka_vs_axis.json")


if __name__ == "__main__":
    main()
