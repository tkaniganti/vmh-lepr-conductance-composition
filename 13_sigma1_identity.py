#!/usr/bin/env python3
"""
13_sigma1_identity.py - does the AP-1 sigma composition gradient track neurotransmitter
identity rather than nucleus of origin or developmental lineage?

12_sigma1_test.py established a controlled, replicated gradient in sigma1A share:
AgRP 0.17, POMC 0.36, VMH SF1 0.55, surviving depth matching and equal-n downsampling, and
replicating across the six conserved VMH classes. The question is what orders it.

A lineage explanation fits badly, because AgRP and POMC neurons share progenitors yet sit at
opposite ends, and because the six VMH classes from one compartment span half the
cross-nuclear range. A neurotransmitter explanation fits the ordering exactly: GABAergic
(AgRP), mixed (POMC), glutamatergic (VMH SF1).

Three tests:
  1. Per-cell, among nuclei where the sigma1 position is detected: partial Spearman between
     a glutamatergic-minus-GABAergic index (Slc17a6 against Gad1/Gad2) and the Ap1s1 share,
     with log10 depth regressed out of both, and a permutation p.
  2. The same correlation computed for the leading paralog of every other testable position,
     so sigma1 is ranked rather than reported alone.
  3. The decisive one: the same test run WITHIN each population or conserved class, where
     between-group identity differences are removed. If sigma1 composition still tracks
     residual neurotransmitter character inside a single nominally glutamatergic class, the
     organizing variable is transmission mode and not nucleus or lineage.

A binned pseudobulk view is also written for plotting.

Output: results/sigma1_identity_<tag>.json, sigma1_identity_<tag>_bins.csv
"""
import os, csv, json, math
from collections import defaultdict
import numpy as np
from scipy import stats
import config as C
import vmhlib as V

N_PERM = 2000
N_BINS = 5
MIN_CELLS_WITHIN = 60

POSITIONS = {
    ("AP1", "large"):  ["Ap1g1", "Ap1g2"],
    ("AP1", "small"):  ["Ap1s1", "Ap1s2", "Ap1s3"],
    ("AP1", "medium"): ["Ap1m1", "Ap1m2"],
    ("AP2", "large"):  ["Ap2a1", "Ap2a2"],
    ("AP3", "beta"):   ["Ap3b1", "Ap3b2"],
    ("AP3", "medium"): ["Ap3m1", "Ap3m2"],
    ("AP3", "small"):  ["Ap3s1", "Ap3s2"],
}
FOCUS = ("AP1", "small")
GLUT = ["Slc17a6"]
GABA = ["Gad1", "Gad2"]
MIN_DETECT = 0.05


def resid(x, z):
    b = np.polyfit(z, x, 1)
    return x - (b[0] * z + b[1])


def partial_spearman(x, y, z, seed=0):
    """Spearman between x and y with z regressed out, plus a permutation p.

    z is rank-transformed before residualizing. Regressing ranks on raw z leaves shared
    curvature, which makes a pair driven purely by z look correlated; on ranks of z the
    relationship is linear and the confound is removed.
    """
    zr = stats.rankdata(z).astype(float)
    rx = resid(stats.rankdata(x).astype(float), zr)
    ry = resid(stats.rankdata(y).astype(float), zr)
    rx = (rx - rx.mean()) / (rx.std() + 1e-12)
    ry = (ry - ry.mean()) / (ry.std() + 1e-12)
    obs = float(np.dot(rx, ry) / rx.size)
    rng = np.random.default_rng(seed)
    hits = 0
    perm = ry.copy()
    for _ in range(N_PERM):
        rng.shuffle(perm)
        if abs(float(np.dot(rx, perm) / rx.size)) >= abs(obs):
            hits += 1
    return obs, (hits + 1) / (N_PERM + 1)


def run(prefix, tag, label_col, relabel, keep):
    print(f"\n=== {prefix} ({tag}) ===")
    feats, name2row = V.load_features(prefix)
    bcs, bc2col = V.load_barcodes(prefix)
    meta = V.load_metadata(prefix)
    apgenes = sorted({g for gs in POSITIONS.values() for g in gs if g in name2row})
    ntgenes = [g for g in GLUT + GABA if g in name2row]
    allg = apgenes + ntgenes
    rows_needed = {name2row[g] for g in allg}

    col2group, _, _ = V.write_col2group(prefix, meta, bc2col, label_col, f"{tag}_ap", keep=keep)
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
    cidx = {c: i for i, c in enumerate(cols)}
    print("   n:", {p: int((grp == p).sum()) for p in pops})

    ent = V.pass_rows(prefix, rows_needed, f"{tag}_ap", "apgenes_nt")
    row2gene = {name2row[g]: g for g in allg}
    counts = {g: np.zeros(len(cols)) for g in allg}
    for l in open(ent):
        r, c, v = l.split()
        g = row2gene.get(int(r))
        i = cidx.get(int(c))
        if g is not None and i is not None:
            counts[g][i] = float(v)

    ln = {g: np.log1p(counts[g] / depth * C.NORM_TARGET) for g in allg}
    glut = np.mean([ln[g] for g in GLUT if g in ln], axis=0)
    gaba = np.mean([ln[g] for g in GABA if g in ln], axis=0)
    nt_index = glut - gaba
    print(f"   neurotransmitter index: mean by group "
          f"{ {p: round(float(nt_index[grp == p].mean()), 3) for p in pops} }")

    testable = {g for g in apgenes if (counts[g] > 0).mean() >= MIN_DETECT}

    def share_vector(genes, mask=None):
        """Per-cell share of genes[0] among `genes`, and the mask of cells where defined."""
        tot = sum(counts[g] for g in genes)
        ok = tot > 0
        if mask is not None:
            ok = ok & mask
        return (counts[genes[0]][ok] / tot[ok]), ok

    # ------------------------------------------------- 1,2. per-cell, all positions
    posres = []
    for (cx, pos), genes in POSITIONS.items():
        usable = [g for g in genes if g in testable]
        if len(usable) < 2:
            continue
        sh, ok = share_vector(usable)
        if ok.sum() < 50:
            continue
        rho, p = partial_spearman(nt_index[ok], sh, logdepth[ok], seed=5)
        posres.append({"complex": cx, "position": pos, "lead_paralog": usable[0],
                       "n_cells": int(ok.sum()), "partial_rho": round(rho, 4), "perm_p": p})
    posres.sort(key=lambda r: -abs(r["partial_rho"]))
    rank = next((i + 1 for i, r in enumerate(posres)
                 if (r["complex"], r["position"]) == FOCUS), None)
    print(f"   AP1-small ranks {rank} of {len(posres)} positions by |partial rho|")

    out = {"dataset": prefix, "tag": tag,
           "n_per_group": {p: int((grp == p).sum()) for p in pops},
           "nt_index_by_group": {p: round(float(nt_index[grp == p].mean()), 4) for p in pops},
           "positions": posres,
           "AP1_small_rank_by_abs_rho": rank}

    # ------------------------------------------------------ 3. within-group version
    sg = [g for g in POSITIONS[FOCUS] if g in testable]
    within = {}
    if len(sg) >= 2:
        for p in pops:
            m = grp == p
            sh, ok = share_vector(sg, mask=m)
            if ok.sum() < MIN_CELLS_WITHIN:
                within[p] = {"n_cells": int(ok.sum()), "status": "too few cells"}
                continue
            rho, pv = partial_spearman(nt_index[ok], sh, logdepth[ok], seed=6)
            within[p] = {"n_cells": int(ok.sum()), "partial_rho": round(rho, 4), "perm_p": pv}
            print(f"   within {p}: n={int(ok.sum())}, rho={rho:+.3f}, p={pv:.4f}")
        vals = [v["partial_rho"] for v in within.values() if "partial_rho" in v]
        out["within_group"] = {
            "by_group": within,
            "n_groups_tested": len(vals),
            "n_positive": int(sum(1 for v in vals if v > 0)),
            "mean_rho": round(float(np.mean(vals)), 4) if vals else None,
            "all_same_sign": bool(vals and all(np.sign(v) == np.sign(vals[0]) for v in vals))}

        # --------------------------------------------------- binned pseudobulk view
        binrows = []
        for scope, mask in [("all", np.ones(len(cols), bool))] + \
                           [(p, grp == p) for p in pops]:
            idx = np.where(mask)[0]
            if idx.size < N_BINS * 20:
                continue
            q = np.quantile(nt_index[idx], np.linspace(0, 1, N_BINS + 1))
            q[-1] += 1e-9
            for b in range(N_BINS):
                sel = idx[(nt_index[idx] >= q[b]) & (nt_index[idx] < q[b + 1])]
                if sel.size == 0:
                    continue
                tot = sum(counts[g][sel].sum() for g in sg)
                if tot <= 0:
                    continue
                rec = {"scope": scope, "bin": b + 1, "n_cells": int(sel.size),
                       "nt_index_mean": round(float(nt_index[sel].mean()), 4),
                       "median_depth": int(np.median(depth[sel]))}
                for g in sg:
                    rec[f"{g}_share"] = round(float(counts[g][sel].sum() / tot), 4)
                binrows.append(rec)
        if binrows:
            with open(os.path.join(C.RESULTS, f"sigma1_identity_{tag}_bins.csv"),
                      "w", newline="") as fh:
                w = csv.DictWriter(fh, fieldnames=list(binrows[0].keys()))
                w.writeheader(); w.writerows(binrows)
            allbins = [r for r in binrows if r["scope"] == "all"]
            if len(allbins) >= 3:
                rho_b, p_b = stats.spearmanr([r["nt_index_mean"] for r in allbins],
                                             [r[f"{sg[0]}_share"] for r in allbins])
                out["binned_trend_all_cells"] = {
                    "n_bins": len(allbins), "spearman_rho": round(float(rho_b), 4),
                    "p": float(p_b),
                    "share_by_bin": [r[f"{sg[0]}_share"] for r in allbins],
                    "nt_index_by_bin": [r["nt_index_mean"] for r in allbins]}

    out["interpretation"] = (
        "A positive partial rho pooled across groups shows only that glutamatergic cells "
        "carry more sigma1A, which group identity alone could produce. The decisive field "
        "is within_group: if the correlation holds inside single populations, with the same "
        "sign in most of them, neurotransmitter character orders sigma1 composition "
        "independently of nucleus or class identity. If within-group correlations are near "
        "zero while the pooled one is positive, the gradient is a between-population "
        "property and the neurotransmitter reading is not established.")

    json.dump(out, open(os.path.join(C.RESULTS, f"sigma1_identity_{tag}.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items() if k != "interpretation"}, indent=2))
    return out


def main():
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    run("GSE172460", "rupp", "Lepr_cluster",
        relabel=lambda g: inv.get(g, ""), keep=lambda g: g in set(C.RUPP_POPS.values()))
    run("GSE172204", "vmh", "VMH_cluster",
        relabel=lambda g: g.split("_")[0], keep=lambda g: bool(g))


if __name__ == "__main__":
    main()
