#!/usr/bin/env python3
"""
47_ap_share_vs_axis.py - does adaptor composition track the effector balance axis?

THE QUESTION. The Discussion names channel trafficking and surface availability as a
determinant that steady-state mRNA cannot report. AP adaptor complexes are the trafficking
machinery, so if within-position paralog share varies along the balance axis, the data are
at least consistent with trafficking differences between nuclei that differ in effector
composition.

WHY THIS IS NOT THE CLUSTERING 2.4 REFUSES. Nothing is partitioned and no groups are named.
The axis is built from channel transcripts; AP subunits appear in neither panel nor in the
control sets, so this is an out-of-sample test of whether the axis predicts anything.

WHAT IT CANNOT SHOW. Adaptor abundance is not adaptor activity, and surface expression is
set post-translationally. A correlation here would mean two transcript sets covary. It
would not establish differential trafficking, and the paper must not say that it does.

DESIGN, AND THE CHOICE BEHIND IT. Within-position paralog share is defined only for nuclei
where the position is detected, and coverage differs sharply by position: AP3 beta 79.7%,
AP1 small 61.4%, AP1 medium 32.8% of the 2,913 Lepr-cluster nuclei. Requiring all six
positions would give a clean AP1-versus-AP3 axis on a subset that is deeply sequenced BY
CONSTRUCTION, and capture depth already tracks the balance axis at rho 0.137, so any
correlation found that way would be partly depth. Each position is therefore scored
separately on whichever nuclei carry it, and six correlations are reported rather than one
composite. Less tidy, and it does not smuggle a depth confound into the answer.

CONTROLS, both of which gate the result:

  depth        each position's subset is compared against the full set for capture depth.
               A position detected only in deep nuclei will inherit the axis-depth
               correlation, so the partial correlation controlling for depth is reported
               alongside the raw one, and the raw one should not be read on its own.
  random sets  detectability-matched gene pairs scored as if they were a position, 200 of
               them, giving the rate at which any pair correlates with the axis at this
               sample size. With up to 2,323 nuclei a rho of 0.05 reaches significance and
               means nothing.

POSITIONS, taken from COMPLEXES in 11_ap_subunits.py, not redefined here. Single-occupant
positions carry no share and are excluded, which is that script's convention: AP1 beta
(Ap1b1) and AP3 large (Ap3d1) are therefore absent, so AP1 and AP3 contribute three
scorable positions each.

Output: results/ap_share_vs_axis.json
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
PREFIX, TAG = "GSE172204", "apx"
N_RANDOM = 200

# from COMPLEXES in 11_ap_subunits.py; single-occupant positions carry no share
POSITIONS = {
    ("AP1", "large"):  ["Ap1g1", "Ap1g2"],
    ("AP1", "medium"): ["Ap1m1", "Ap1m2"],
    ("AP1", "small"):  ["Ap1s1", "Ap1s2", "Ap1s3"],
    ("AP3", "beta"):   ["Ap3b1", "Ap3b2"],
    ("AP3", "medium"): ["Ap3m1", "Ap3m2"],
    ("AP3", "small"):  ["Ap3s1", "Ap3s2"],
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
    logd = np.log10(np.where(np.isnan(depth), 1, depth))

    rho_axis_depth = float(stats.spearmanr(axis, logd).statistic)
    print(f">> axis vs log depth over these nuclei: rho {rho_axis_depth:+.4f}", flush=True)

    # one pass for every gene we will need: the AP paralogs plus a random-pair pool
    ap_genes = sorted({g for gs in POSITIONS.values() for g in gs if g in name2row})
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    pool_rows = [r for r in rowsum if rowsum[r] > 0]
    extra = list(rng.choice(pool_rows, min(3000, len(pool_rows)), replace=False))
    need = {name2row[g] for g in ap_genes} | {int(r) for r in extra}
    ent = V.pass_rows(PREFIX, need, TAG, "apx")

    counts = defaultdict(lambda: np.zeros(len(cols)))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset:
            counts[int(r)][cidx[c]] = float(v)

    def share_and_corr(genes_rows, label):
        """Within-position share for nuclei where the position is detected."""
        mat = np.vstack([counts[r] for r in genes_rows])
        tot = mat.sum(0)
        have = tot > 0
        if have.sum() < 100:
            return None
        sh = mat[0][have] / tot[have]          # share of the FIRST paralog listed
        ax = axis[have]
        ld = logd[have]
        rho = float(stats.spearmanr(sh, ax).statistic)
        # partial, controlling for depth: residualise both on rank of log depth
        rs_ = stats.rankdata(ld)
        def resid(x):
            x = stats.rankdata(x).astype(float)
            b = np.polyfit(rs_, x, 1)
            return x - np.polyval(b, rs_)
        prho = float(stats.spearmanr(resid(sh), resid(ax)).statistic)
        return {"n": int(have.sum()), "pct_of_nuclei": round(100 * have.mean(), 1),
                "median_depth_subset": round(float(np.median(depth[have])), 1),
                "median_depth_all": round(float(np.median(depth)), 1),
                "rho_raw": round(rho, 4), "rho_partial_depth": round(prho, 4),
                "share_of": label}

    out = {"n_nuclei": len(cols), "draws": 1,
           "axis_vs_depth_rho": round(rho_axis_depth, 4),
           "positions": {},
           "caveat": ("adaptor abundance is not adaptor activity; a correlation here "
                      "means two transcript sets covary and does not establish "
                      "differential trafficking")}

    print(f"\n{'position:paralog':<22}{'n':>7}{'%':>6}{'rho raw':>10}"
          f"{'rho|depth':>11}{'depth s/a':>14}")
    for (cx, pos), gs in POSITIONS.items():
        rows = [name2row[g] for g in gs if g in name2row]
        if len(rows) < 2:
            continue
        # Every occupant gets its own share, not just the first. A three-way position
        # has three shares, and Ap1s2 or Ap1s3 could track the axis where Ap1s1 does
        # not. Scoring only gs[0] was a limitation of the first version.
        present = [g for g in gs if g in name2row]
        for lead in range(len(present)):
            ordered = [name2row[present[lead]]] + \
                      [name2row[g] for i, g in enumerate(present) if i != lead]
            r = share_and_corr(ordered, present[lead])
            key = f"{cx}_{pos}:{present[lead]}"
            out["positions"][key] = r
            if r is None:
                print(f"{key:<22}{'too few nuclei':>34}")
                continue
            print(f"{key:<22}{r['n']:>7}{r['pct_of_nuclei']:>6.1f}"
                  f"{r['rho_raw']:>10.4f}{r['rho_partial_depth']:>11.4f}"
                  f"{r['median_depth_subset']/r['median_depth_all']:>13.2f}x")

    # ---- control: how often does any detectability-matched pair correlate?
    print(f"\n>> random-pair control, {N_RANDOM} pairs", flush=True)
    nulls = []
    ex = [int(r) for r in extra if counts[int(r)].sum() > 0]
    for _ in range(N_RANDOM):
        pair = list(rng.choice(ex, 2, replace=False))
        r = share_and_corr(pair, "random")
        if r:
            nulls.append(abs(r["rho_partial_depth"]))
    nulls = np.array(nulls)
    p95 = float(np.percentile(nulls, 95)) if nulls.size else float("nan")
    out["random_pair_abs_partial_rho_p95"] = round(p95, 4)
    out["n_random_pairs_scored"] = int(nulls.size)
    print(f"   |rho| partial, 95th percentile of random pairs: {p95:.4f}")

    real = {k: v for k, v in out["positions"].items() if v}
    beats = [k for k, v in real.items() if abs(v["rho_partial_depth"]) > p95]
    out["positions_exceeding_random"] = beats
    print(f"\n{'=' * 66}")
    if not beats:
        print("  No position's share tracks the axis beyond what arbitrary gene pairs do.")
        print("  Adaptor composition carries no information about effector balance in")
        print("  this population, which is a further bound on what the axis reflects.")
    else:
        print(f"  Exceeding the random-pair null: {', '.join(beats)}")
        print("  Worth following up, and worth stating carefully: this is covariation")
        print("  between two transcript sets, not evidence of differential trafficking.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "ap_share_vs_axis.json"), "w"),
              indent=2)
    print("\nWritten to results/ap_share_vs_axis.json")


if __name__ == "__main__":
    main()
