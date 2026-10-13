#!/usr/bin/env python3
"""
14_sigma1_paralogs.py - what drives the sigma1 composition gradient, and what covaries
with it once neurotransmitter identity is excluded?

12 established the gradient and 13 ruled out transmission mode as its organizer. The
detection logistics in 12 pointed somewhere specific: in the Rupp populations, Ap1s1
detection did NOT differ across populations after depth control (LR 2.67, p = 0.26) while
Ap1s2 and Ap1s3 both did strongly (p around 2e-13). So the share shift may be the two minor
paralogs falling rather than sigma1A rising, which is a different claim with a different
mechanism behind it.

Four analyses:
  A. Decomposition. Absolute expression per paralog by group with bootstrap CIs, beside the
     shares, so "sigma1A higher in VMH" and "sigma1B/C lower in VMH" can be told apart.
  B. Paralog covariation. Partial Spearman between each paralog pair within each group,
     depth-residualized. If Ap1s2 and Ap1s3 move together and independently of Ap1s1, that
     points to a shared regulatory program over the two minor paralogs.
  C. Genome-wide search for correlates of the sigma1A share, computed WITHIN groups so
     between-population structure cannot drive it, then averaged across groups weighted by
     cell count. This is the step that could name the organizing variable 13 failed to find.
  D. Abundance control. Is the share itself a function of how much of the position was
     captured in a cell? If so, the per-cell share is partly a detection artifact.

Partial correlations residualize on the RANK of log depth, which is the corrected form; see
13 for why regressing on raw depth leaves a confound.

Output: results/sigma1_paralogs_<tag>.json, _decomposition.csv, _correlates.csv
"""
import os, csv, json, math
from collections import defaultdict
import numpy as np
from scipy import stats
import config as C
import vmhlib as V

SIGMA = ["Ap1s1", "Ap1s2", "Ap1s3"]
N_BOOT = 2000
N_PERM = 2000
MIN_DETECT_FRAC = 0.05
MAX_DETECT_FRAC = 0.95
MIN_CELLS = 100
TOP_N = 200


def resid(x, z):
    b = np.polyfit(z, x, 1)
    return x - (b[0] * z + b[1])


def partial_spearman(x, y, z):
    zr = stats.rankdata(z).astype(float)
    rx = resid(stats.rankdata(x).astype(float), zr)
    ry = resid(stats.rankdata(y).astype(float), zr)
    sx, sy = rx.std(), ry.std()
    if sx <= 0 or sy <= 0:
        return 0.0
    return float(np.dot((rx - rx.mean()) / sx, (ry - ry.mean()) / sy) / rx.size)


def perm_p(x, y, z, obs, seed=0):
    rng = np.random.default_rng(seed)
    zr = stats.rankdata(z).astype(float)
    rx = resid(stats.rankdata(x).astype(float), zr)
    ry = resid(stats.rankdata(y).astype(float), zr)
    rx = (rx - rx.mean()) / (rx.std() + 1e-12)
    ry = (ry - ry.mean()) / (ry.std() + 1e-12)
    hits, p = 0, ry.copy()
    for _ in range(N_PERM):
        rng.shuffle(p)
        if abs(float(np.dot(rx, p) / rx.size)) >= abs(obs):
            hits += 1
    return (hits + 1) / (N_PERM + 1)


def boot_mean(x, seed=0):
    rng = np.random.default_rng(seed)
    m = rng.choice(np.asarray(x, float), (N_BOOT, len(x)), replace=True).mean(1)
    return float(np.mean(x)), float(np.percentile(m, 2.5)), float(np.percentile(m, 97.5))


def extract_block(prefix, tag, cols):
    """One pass over the matrix keeping only the requested columns. Header detected."""
    out = V.dpath(f"block_{tag}.txt")
    colfile = V.dpath(f"blockcols_{tag}.txt")
    with open(colfile, "w") as fh:
        fh.write("\n".join(str(c) for c in sorted(cols)) + "\n")
    if not os.path.exists(out):
        V.sh("zcat " + V.dpath(f"{prefix}_matrix.mtx.gz") + " | awk '"
             "NR==FNR{keep[$1]=1; next} /^%/{next} !dims{dims=1; next} "
             "{ if($2 in keep) print $1, $2, $3 }' " + colfile + " - > " + out)
    return out


def run(prefix, tag, label_col, relabel, keep):
    print(f"\n=== {prefix} ({tag}) ===")
    feats, name2row = V.load_features(prefix)
    bcs, bc2col = V.load_barcodes(prefix)
    meta = V.load_metadata(prefix)
    sym = {v: k for k, v in name2row.items()}

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
    n = len(cols)
    print("   n:", {p: int((grp == p).sum()) for p in pops})

    block = extract_block(prefix, f"{tag}_blk", set(cols))
    nz = defaultdict(int)
    with open(block) as fh:
        for line in fh:
            r, c, v = line.split()
            if float(v) > 0:
                nz[int(r)] += 1
    lo, hi = MIN_DETECT_FRAC * n, MAX_DETECT_FRAC * n
    keeprows = sorted([r for r, k in nz.items() if lo <= k <= hi])
    srows = [name2row[g] for g in SIGMA if g in name2row]
    need = sorted(set(keeprows) | set(srows))
    kidx = {r: i for i, r in enumerate(need)}
    print(f"   {len(nz)} genes with signal; {len(keeprows)} in detection window")

    X = np.zeros((len(need), n), dtype=np.float32)
    with open(block) as fh:
        for line in fh:
            r, c, v = line.split()
            i = kidx.get(int(r))
            j = cidx.get(int(c))
            if i is not None and j is not None:
                X[i, j] = float(v)
    ln = np.log1p(X / depth * C.NORM_TARGET)

    sv = {g: X[kidx[name2row[g]]] for g in SIGMA if g in name2row}
    slv = {g: ln[kidx[name2row[g]]] for g in SIGMA if g in name2row}
    tot = sum(sv[g] for g in sv)
    detected = tot > 0

    # ------------------------------------------------------------ A. decomposition
    dec = []
    for p in pops:
        m = grp == p
        for g in sv:
            mean, lo_, hi_ = boot_mean(slv[g][m], seed=1)
            md = m & detected
            dec.append({"population": p, "gene": g, "n_cells": int(m.sum()),
                        "log_cp10k_mean": round(mean, 4),
                        "log_cp10k_ci_lo": round(lo_, 4), "log_cp10k_ci_hi": round(hi_, 4),
                        "detect_pct": round(100 * float((sv[g][m] > 0).mean()), 1),
                        "share_pseudobulk": round(float(sv[g][m].sum() /
                                                        max(tot[m].sum(), 1e-9)), 4)})
    with open(os.path.join(C.RESULTS, f"sigma1_paralogs_{tag}_decomposition.csv"),
              "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(dec[0].keys()))
        w.writeheader(); w.writerows(dec)

    # ------------------------------------------------------- B. paralog covariation
    cov = {}
    for p in pops:
        m = (grp == p) & detected
        if m.sum() < MIN_CELLS:
            continue
        pair = {}
        for a in range(len(SIGMA)):
            for b in range(a + 1, len(SIGMA)):
                ga, gb = SIGMA[a], SIGMA[b]
                if ga in slv and gb in slv:
                    r = partial_spearman(slv[ga][m], slv[gb][m], logdepth[m])
                    pair[f"{ga}_vs_{gb}"] = {"rho": round(r, 4),
                                             "perm_p": perm_p(slv[ga][m], slv[gb][m],
                                                              logdepth[m], r, seed=2)}
        cov[p] = {"n_cells": int(m.sum()), "pairs": pair}

    # --------------------------------------------- C. within-group genome-wide scan
    share = np.zeros(n)
    share[detected] = sv["Ap1s1"][detected] / tot[detected]
    per_gene = defaultdict(list)
    for p in pops:
        m = (grp == p) & detected
        if m.sum() < MIN_CELLS:
            continue
        y = share[m]
        z = logdepth[m]
        for r in keeprows:
            if r in srows:
                continue
            x = ln[kidx[r]][m]
            if x.std() <= 0:
                continue
            per_gene[r].append((int(m.sum()), partial_spearman(x, y, z)))
        print(f"   scanned {p} (n={int(m.sum())})")

    corr = []
    for r, vals in per_gene.items():
        if len(vals) < max(2, len(pops) - 1):
            continue
        w_ = np.array([v[0] for v in vals], float)
        rh = np.array([v[1] for v in vals], float)
        mean_rho = float(np.average(rh, weights=w_))
        t = mean_rho * math.sqrt(max(w_.sum() - 3, 1) / max(1 - mean_rho ** 2, 1e-12))
        corr.append({"gene": sym.get(r, str(r)), "row": r,
                     "mean_within_rho": round(mean_rho, 4),
                     "n_groups": len(vals),
                     "n_same_sign": int(max((rh > 0).sum(), (rh < 0).sum())),
                     "p_approx": float(2 * stats.t.sf(abs(t), max(int(w_.sum()) - 3, 1)))})
    if corr:
        pv = np.array([c["p_approx"] for c in corr])
        q = (stats.false_discovery_control(pv, method="bh")
             if hasattr(stats, "false_discovery_control") else np.full(pv.size, np.nan))
        for c, qq in zip(corr, q):
            c["q_bh"] = float(qq)
        corr.sort(key=lambda c: -abs(c["mean_within_rho"]))
        with open(os.path.join(C.RESULTS, f"sigma1_paralogs_{tag}_correlates.csv"),
                  "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=["gene", "row", "mean_within_rho", "n_groups",
                                               "n_same_sign", "p_approx", "q_bh"])
            w.writeheader(); w.writerows(corr[:TOP_N])

    # ------------------------------------------------------------ D. abundance check
    ab = {}
    for p in pops:
        m = (grp == p) & detected
        if m.sum() < MIN_CELLS:
            continue
        r1 = partial_spearman(tot[m], share[m], logdepth[m])
        r2 = float(stats.spearmanr(depth[m], share[m]).statistic)
        ab[p] = {"rho_share_vs_position_counts_given_depth": round(r1, 4),
                 "rho_share_vs_depth": round(r2, 4),
                 "median_position_counts": float(np.median(tot[m]))}

    out = {"dataset": prefix, "tag": tag,
           "n_per_group": {p: int((grp == p).sum()) for p in pops},
           "n_cells_with_position_detected": int(detected.sum()),
           "paralog_covariation": cov,
           "abundance_check": ab,
           "n_genes_scanned": len(corr),
           "top_within_group_correlates": [
               {k: c[k] for k in ("gene", "mean_within_rho", "n_groups",
                                  "n_same_sign", "q_bh")} for c in corr[:25]],
           "n_correlates_q05": int(sum(1 for c in corr if c.get("q_bh", 1) < 0.05))}
    json.dump(out, open(os.path.join(C.RESULTS, f"sigma1_paralogs_{tag}.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k != "top_within_group_correlates"}, indent=2))
    print("top correlates:", [c["gene"] for c in corr[:15]])
    return out


def main():
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    run("GSE172460", "rupp", "Lepr_cluster",
        relabel=lambda g: inv.get(g, ""), keep=lambda g: g in set(C.RUPP_POPS.values()))
    run("GSE172204", "vmh", "VMH_cluster",
        relabel=lambda g: g.split("_")[0], keep=lambda g: bool(g))


if __name__ == "__main__":
    main()
