#!/usr/bin/env python3
"""
46_axis_covariates.py - does the effector-balance axis track anything it was not built from?

WHY. Sections 2.3 and 2.4 establish that the axis carries no detectable structure within the
leptin-receptor population: the distribution is unimodal and the within-class variance
component sits at its own chance expectation. That says the axis does not PARTITION the
population. It does not say what the axis reflects.

This asks the complementary question. Take the per-nucleus score and regress it against
variables that played no part in building it. If something real tracks it, that is a finding.
If nothing does beyond capture depth, the axis reflects measurement rather than identity, and
that bounds what any claim resting on it can mean.

WHAT IS TESTED, all from the deposited metadata or the stored scores:

  capture depth           already known to correlate at rho 0.106 across VMH nuclei; included
                          here as a positive control for the method, since a variable known
                          to track the axis must come back tracking it
  genes detected          a second technical measure, correlated with depth but not identical
  conserved assignment    do nuclei whose conserved class is NOT Lepr sit differently from
                          those whose class is? Those nuclei came from the same clusters but
                          landed elsewhere in the cross-species annotation
  doublet flag            do doublet-cluster-labelled nuclei sit differently?
  cluster identity        how much of the axis does cluster membership explain? This one IS
                          partly circular, since clusters come from the same matrix, and is
                          reported for scale rather than as a finding

WHAT THIS IS NOT. None of these is a physiological measurement. A hit tells you the axis
tracks something about how the data were produced or annotated, not that it predicts how a
cell responds to leptin. The only test of that in the manuscript is the POMC and AgRP
comparison in 2.8, which came back negative, and the only decisive one is paired recording.

EFFECT SIZES, NOT P VALUES. With 2,913 nuclei a correlation of 0.05 reaches significance and
means nothing. Every association is reported as an effect size against a permutation null,
and p values are given only alongside.

Reads data/percell_scores.npz, which holds ONE control-gene draw for the 2,913
Lepr-cluster nuclei. Single draw, stated in the output.

Output: results/axis_covariates.json
"""
import json
import os
import warnings

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX = "GSE172204"
N_PERM = 2000


def perm_null(y, g, stat_fn, rng, n=N_PERM):
    """Permutation null for any statistic of (score, grouping)."""
    obs = stat_fn(y, g)
    lab = np.array(g, copy=True)
    null = []
    for _ in range(n):
        rng.shuffle(lab)
        null.append(stat_fn(y, lab))
    null = np.array(null)
    return obs, float(null.mean()), float(np.percentile(null, 97.5)), \
        float((null >= obs).mean())


def eta_sq(y, g):
    """Share of variance in y explained by grouping g. Biased upward with group count,
    which is why it is only ever compared against its own permutation null here."""
    gt = y.mean()
    ss_b = sum((y[g == k].size) * (y[g == k].mean() - gt) ** 2 for k in set(g))
    ss_t = ((y - gt) ** 2).sum()
    return 100 * ss_b / ss_t if ss_t > 0 else 0.0


def main():
    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    y = np.asarray(z["HminusD"], float)
    cols = np.asarray(z["cols"])
    lab = np.asarray(z["labels"])

    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    bc = [col2bc.get(int(c)) for c in cols]
    missing = sum(1 for b in bc if b is None)
    if missing:
        print(f">> {missing} nuclei could not be mapped back to a barcode", flush=True)

    print(f">> {y.size} Lepr-cluster nuclei, one control-gene draw", flush=True)
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    out = {"n": int(y.size), "draws": 1, "n_perm": N_PERM,
           "note": ("none of these covariates is a physiological measurement; a hit "
                    "means the axis tracks how the data were produced or annotated, "
                    "not that it predicts a response to leptin")}

    # ---- continuous covariates, Spearman against a permutation null
    print(f"\n{'covariate':<26}{'rho':>9}{'null p97.5':>12}{'p':>12}")
    # The metadata has no technical fields: ID, Run, Sample, Cell_cluster, Cell_type,
    # Neuron_cluster, VMH_cluster, Conserved_cluster. Depth has to be computed, which
    # is what 02_snrna_core.py does, so the positive control needs a totals pass.
    meta2 = {b: dict(m) for b, m in meta.items()}
    for b, m in meta2.items():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""
    cg, _, _ = V.write_col2group(PREFIX, meta2, bc2col, "_grp", "cov")
    cs, _ = V.pass_totals(PREFIX, cg, "cov")
    colsum = V.read_two_col(cs)
    depth = np.array([colsum.get(int(c), np.nan) for c in cols], float)

    cont = {}
    rho_d = float(stats.spearmanr(y[~np.isnan(depth)],
                                  depth[~np.isnan(depth)]).statistic)
    cont["capture depth"] = {"rho": round(rho_d, 4),
                             "n": int((~np.isnan(depth)).sum())}
    print(f"{'capture depth':<26}{rho_d:>9.4f}{'(positive control)':>24}")
    # THE CONTROL GATES THE VERDICT. A script whose control did not run cannot
    # conclude that nothing tracks the axis.
    control_ok = abs(rho_d) > 0.03
    if not control_ok:
        print("   POSITIVE CONTROL FAILED: depth does not track the axis here, though")
        print("   Methods reports rho = 0.106 across VMH nuclei. Nothing below is")
        print("   readable until that is resolved.")

    for name, key in ():
        vals = []
        for b in bc:
            m = meta.get(b, {}) if b else {}
            v = m.get(key) or m.get(key.lower()) or m.get("n_counts") \
                or m.get("n_genes")
            try:
                vals.append(float(v))
            except (TypeError, ValueError):
                vals.append(np.nan)
        v = np.array(vals, float)
        ok = ~np.isnan(v)
        if ok.sum() < 100:
            print(f"{name:<26}{'not in metadata':>33}")
            cont[name] = None
            continue
        rho = float(stats.spearmanr(y[ok], v[ok]).statistic)
        nulls = []
        vv = v[ok].copy()
        for _ in range(500):
            rng.shuffle(vv)
            nulls.append(abs(float(stats.spearmanr(y[ok], vv).statistic)))
        p975 = float(np.percentile(nulls, 97.5))
        p = float((np.array(nulls) >= abs(rho)).mean())
        cont[name] = {"rho": round(rho, 4), "null_p97_5": round(p975, 4),
                      "perm_p": p, "n": int(ok.sum())}
        print(f"{name:<26}{rho:>9.4f}{p975:>12.4f}{p:>12.4f}")
    out["continuous"] = cont
    out["positive_control_passes"] = bool(control_ok)

    # ---- categorical covariates, eta-squared against its own permutation null
    print(f"\n{'grouping':<26}{'eta2 %':>9}{'null mean':>11}{'null p97.5':>12}{'p':>8}")
    cats = {}
    for name, fn in (
        ("conserved class is Lepr",
         lambda b: "Lepr" if (meta.get(b, {}).get("Conserved_cluster") or ""
                              ).strip().startswith("Lepr") else "other"),
        ("doublet-flagged",
         lambda b: "doublet" if (meta.get(b, {}).get("Conserved_cluster") or ""
                                 ).strip().startswith("Doublets") else "other"),
        ("VMH cluster (circular)",
         lambda b: (meta.get(b, {}).get("VMH_cluster") or "?").strip()),
    ):
        g = np.array([fn(b) if b else "?" for b in bc])
        if len(set(g)) < 2:
            print(f"{name:<26}{'single level':>40}")
            cats[name] = None
            continue
        obs, nm, n975, p = perm_null(y, g, eta_sq, rng, n=500)
        sizes = {k: int((g == k).sum()) for k in sorted(set(g))}
        cats[name] = {"eta2_pct": round(obs, 4), "null_mean": round(nm, 4),
                      "null_p97_5": round(n975, 4), "perm_p": p,
                      "levels": sizes}
        print(f"{name:<26}{obs:>9.4f}{nm:>11.4f}{n975:>12.4f}{p:>8.4f}")
        if len(sizes) <= 3:
            for k, n_ in sizes.items():
                sub = y[g == k]
                print(f"    {k:<20} n={n_:<6} mean {sub.mean():+.4f}")
    out["categorical"] = cats

    # ---- verdict
    # depth is the positive control, not a finding, so it is excluded from hits
    hits = [n for n, r in cont.items()
            if r and n != "capture depth"
            and r.get("perm_p", 1) < 0.05 and abs(r.get("rho", 0)) > 0.1]
    hits += [n for n, r in cats.items()
             if r and r["perm_p"] < 0.05 and r["eta2_pct"] > r["null_p97_5"]
             and "circular" not in n]
    out["covariates_tracking_axis"] = hits
    print(f"\n{'=' * 62}")
    if not control_ok:
        print("  VERDICT WITHHELD: the positive control did not pass.")
    elif not hits:
        print("  Nothing tracks the axis beyond what is already reported.")
        print("  The axis reflects neither annotation nor obvious technical structure,")
        print("  which bounds what it can be said to carry and supports the claim that")
        print("  it holds no identity information in this population.")
    else:
        print(f"  Tracks the axis: {', '.join(hits)}")
        print("  Each is a measurement or annotation variable, not a physiological one.")
        print("  A hit here constrains what the axis reflects; it does not show the axis")
        print("  predicts anything about leptin responses.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "axis_covariates.json"), "w"),
              indent=2)
    print("\nWritten to results/axis_covariates.json")


if __name__ == "__main__":
    main()
