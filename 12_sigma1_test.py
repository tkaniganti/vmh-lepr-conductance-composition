#!/usr/bin/env python3
"""
12_sigma1_test.py - does AP-1 sigma composition really differ across hypothalamic neuron
populations, or is the gradient an artifact of depth and sample size?

11_ap_subunits.py computed within-position paralog shares but tested none of them. The
AP1-small position (Ap1s1 sigma1A, Ap1s2 sigma1B, Ap1s3 sigma1C) showed the largest spread
of any position in either dataset, with the sigma1A share rising from AgRP to POMC to
VMH SF1 and replicating across the Affinati VMH classes. That is a point estimate with no
error bars, from populations differing 1.5-fold in median depth, for three paralogs of
different abundance, which is precisely the configuration that produced the first, wrong
reading of Trpc5.

Controls, all of which the claim must survive:
  1. bootstrap CIs on each share, resampling cells within population
  2. permutation test of composition difference, shuffling population labels across cells,
     using total variation distance between each population's share vector and the pooled one
  3. the same test for every other testable position, so sigma1 is ranked among them rather
     than reported alone
  4. depth quantile matching across populations, since share is sensitive to differential
     dropout between paralogs of different abundance
  5. equal-n downsampling to the smallest population
  6. per-subunit detection logistic regression on population given log10 depth

Output: results/sigma1_<tag>.json, sigma1_<tag>_positions.csv, sigma1_<tag>_shares.csv
"""
import os, csv, json, math
from collections import defaultdict
import numpy as np
import pandas as pd
import statsmodels.api as sm
from scipy.stats import chi2 as chi2dist
import config as C
import vmhlib as V

N_BOOT = 2000
N_PERM = 2000
N_DOWNSAMPLE = 200
MIN_DETECT = 0.05

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


def shares_from(counts, genes, idx):
    tot = sum(counts[g][idx].sum() for g in genes)
    if tot <= 0:
        return None
    return np.array([counts[g][idx].sum() / tot for g in genes])


def tvd_statistic(counts, genes, grp, pops, cell_idx):
    """Mean total variation distance between each population's share and the pooled share."""
    pooled = shares_from(counts, genes, cell_idx)
    if pooled is None:
        return None
    acc, n = 0.0, 0
    for p in pops:
        sub = cell_idx[grp[cell_idx] == p]
        s = shares_from(counts, genes, sub)
        if s is None:
            continue
        acc += 0.5 * np.abs(s - pooled).sum()
        n += 1
    return acc / n if n else None


def perm_test_composition(counts, genes, grp, pops, cell_idx, seed=0):
    obs = tvd_statistic(counts, genes, grp, pops, cell_idx)
    if obs is None:
        return None, None
    rng = np.random.default_rng(seed)
    g = grp[cell_idx].copy()
    hits = 0
    for _ in range(N_PERM):
        rng.shuffle(g)
        shuffled = grp.copy()
        shuffled[cell_idx] = g
        v = tvd_statistic(counts, genes, shuffled, pops, cell_idx)
        if v is not None and v >= obs:
            hits += 1
    return float(obs), (hits + 1) / (N_PERM + 1)


def boot_share(counts, genes, idx, seed=0):
    rng = np.random.default_rng(seed)
    out = []
    for _ in range(N_BOOT):
        pick = rng.choice(idx, size=idx.size, replace=True)
        s = shares_from(counts, genes, pick)
        if s is not None:
            out.append(s)
    if not out:
        return None
    A = np.array(out)
    return A.mean(0), np.percentile(A, 2.5, axis=0), np.percentile(A, 97.5, axis=0)


def quantile_match(grp, pops, depth, cell_idx, n_bins=10, seed=11):
    edges = np.quantile(depth[cell_idx], np.linspace(0, 1, n_bins + 1))
    edges[-1] += 1e-9
    rng = np.random.default_rng(seed)
    keep = []
    for b in range(n_bins):
        inbin = {p: cell_idx[(grp[cell_idx] == p) &
                             (depth[cell_idx] >= edges[b]) &
                             (depth[cell_idx] < edges[b + 1])] for p in pops}
        k = min(v.size for v in inbin.values())
        if k == 0:
            continue
        for p in pops:
            keep.extend(rng.choice(inbin[p], size=k, replace=False).tolist())
    return np.array(sorted(keep))


def lr_group_given_depth(y, logdepth, group):
    Xd = sm.add_constant(logdepth.reshape(-1, 1), has_constant="add")
    m0 = sm.Logit(y, Xd).fit(disp=0, maxiter=200)
    D = pd.get_dummies(group, drop_first=True).to_numpy(dtype=float)
    m1 = sm.Logit(y, np.column_stack([Xd, D])).fit(disp=0, maxiter=200)
    lr = 2 * (m1.llf - m0.llf)
    return float(lr), int(D.shape[1]), float(chi2dist.sf(max(lr, 0.0), D.shape[1]))


def run(prefix, tag, label_col, relabel, keep):
    print(f"\n=== {prefix} ({tag}) ===")
    feats, name2row = V.load_features(prefix)
    bcs, bc2col = V.load_barcodes(prefix)
    meta = V.load_metadata(prefix)
    allgenes = sorted({g for gs in POSITIONS.values() for g in gs if g in name2row})
    rows_needed = {name2row[g] for g in allgenes}

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
    allidx = np.arange(len(cols))
    print("   n:", {p: int((grp == p).sum()) for p in pops},
          " median depth:", {p: int(np.median(depth[grp == p])) for p in pops})

    ent = V.pass_rows(prefix, rows_needed, f"{tag}_ap", "apgenes")   # cached from 11
    row2gene = {name2row[g]: g for g in allgenes}
    cidx = {c: i for i, c in enumerate(cols)}
    counts = {g: np.zeros(len(cols)) for g in allgenes}
    for l in open(ent):
        r, c, v = l.split()
        # the cached entries file may hold more rows than this script asks for
        # (11_ap_subunits.py caches all AP genes), so skip anything unrequested
        g = row2gene.get(int(r))
        i = cidx.get(int(c))
        if g is not None and i is not None:
            counts[g][i] = float(v)

    testable = {g for g in allgenes if (counts[g] > 0).mean() >= MIN_DETECT}
    print(f"   testable subunits: {sorted(testable)}")

    # ------------------------------------------- 2,3. all positions, permutation ranked
    posrows = []
    for (cx, pos), genes in POSITIONS.items():
        usable = [g for g in genes if g in testable]
        if len(usable) < 2:
            posrows.append({"complex": cx, "position": pos, "paralogs": ";".join(genes),
                            "status": "untestable"})
            continue
        obs, p = perm_test_composition(counts, usable, grp, pops, allidx, seed=3)
        sh = {q: shares_from(counts, usable, allidx[grp == q]) for q in pops}
        spread = {usable[i]: round(float(max(sh[q][i] for q in pops) -
                                         min(sh[q][i] for q in pops)), 4)
                  for i in range(len(usable))}
        posrows.append({"complex": cx, "position": pos, "paralogs": ";".join(usable),
                        "status": "ok", "tvd": round(obs, 4), "perm_p": p,
                        "max_share_spread": max(spread.values()),
                        "spread_by_paralog": ";".join(f"{k}={v}" for k, v in spread.items())})
    posrows.sort(key=lambda r: -(r.get("tvd") or -1))
    with open(os.path.join(C.RESULTS, f"sigma1_{tag}_positions.csv"), "w", newline="") as fh:
        fields = ["complex", "position", "paralogs", "status", "tvd", "perm_p",
                  "max_share_spread", "spread_by_paralog"]
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader(); w.writerows(posrows)
    rank = next((i + 1 for i, r in enumerate(posrows)
                 if (r["complex"], r["position"]) == FOCUS), None)
    print(f"   AP1-small ranks {rank} of {len(posrows)} positions by TVD")

    # ---------------------------------------------------- 1. bootstrap CIs on sigma1
    sg = [g for g in POSITIONS[FOCUS] if g in testable]
    shrows = []
    if len(sg) >= 2:
        for q in pops:
            idx = allidx[grp == q]
            b = boot_share(counts, sg, idx, seed=5)
            s = shares_from(counts, sg, idx)
            for i, g in enumerate(sg):
                shrows.append({"population": q, "gene": g,
                               "share": round(float(s[i]), 4),
                               "boot_mean": round(float(b[0][i]), 4),
                               "ci_lo": round(float(b[1][i]), 4),
                               "ci_hi": round(float(b[2][i]), 4),
                               "n_cells": int(idx.size)})
        with open(os.path.join(C.RESULTS, f"sigma1_{tag}_shares.csv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(shrows[0].keys()))
            w.writeheader(); w.writerows(shrows)

    # ------------------------------------------------- 4. depth matched, 5. equal n
    out = {"dataset": prefix, "tag": tag,
           "n_per_population": {p: int((grp == p).sum()) for p in pops},
           "median_depth": {p: float(np.median(depth[grp == p])) for p in pops},
           "testable_subunits": sorted(testable),
           "position_rank_of_AP1_small_by_tvd": rank,
           "n_positions_tested": sum(1 for r in posrows if r["status"] == "ok")}

    if len(sg) >= 2:
        qidx = quantile_match(grp, pops, depth, allidx)
        obs_q, p_q = perm_test_composition(counts, sg, grp, pops, qidx, seed=4)
        out["depth_quantile_matched"] = {
            "n_per_population": int(qidx.size / len(pops)),
            "tvd": round(obs_q, 4) if obs_q else None, "perm_p": p_q,
            "shares": {q: ([round(float(x), 4) for x in _s]
                           if (_s := shares_from(counts, sg, qidx[grp[qidx] == q])) is not None
                           else None)
                       for q in pops},
            "paralogs": sg}

        rng = np.random.default_rng(C.RANDOM_SEED)
        nmin = min(int((grp == p).sum()) for p in pops)
        ds = defaultdict(list)
        for _ in range(N_DOWNSAMPLE):
            pick = np.concatenate([rng.choice(allidx[grp == p], nmin, replace=False)
                                   for p in pops])
            for q in pops:
                s = shares_from(counts, sg, pick[grp[pick] == q])
                if s is not None:
                    ds[q].append(s)
        out["equal_n"] = {"n_per_population": int(nmin), "n_draws": N_DOWNSAMPLE,
                          "paralogs": sg,
                          "share_median": {q: [round(float(x), 4)
                                               for x in np.median(np.array(v), axis=0)]
                                           for q, v in ds.items()},
                          "share_ci": {q: [[round(float(a), 4), round(float(b_), 4)]
                                           for a, b_ in zip(
                                               np.percentile(np.array(v), 2.5, axis=0),
                                               np.percentile(np.array(v), 97.5, axis=0))]
                                       for q, v in ds.items()}}

        # ------------------------------------------- 6. detection logistic per subunit
        det = {}
        for g in sg:
            y = (counts[g] > 0).astype(float)
            try:
                lr, df, pv = lr_group_given_depth(y, logdepth, grp)
                det[g] = {"detect_pct": round(100 * float(y.mean()), 1),
                          "LR_group_given_depth": round(lr, 2), "df": df, "p": pv}
            except Exception as e:
                det[g] = {"error": str(e)}
        out["detection_logistic"] = det

        tv = next(r for r in posrows if (r["complex"], r["position"]) == FOCUS)
        out["AP1_small_overall"] = {"tvd": tv["tvd"], "perm_p": tv["perm_p"],
                                    "shares": {q: [round(float(x), 4) for x in
                                                   shares_from(counts, sg, allidx[grp == q])]
                                               for q in pops}}

    out["interpretation"] = (
        "The sigma1 gradient stands only if the permutation p is small, AP1-small ranks at "
        "or near the top of the positions by TVD, the bootstrap CIs for the extreme "
        "populations do not overlap, and both the depth-matched and equal-n versions "
        "preserve the ordering. Failure of the depth-matched test specifically would mean "
        "the gradient reflects differential dropout among paralogs of different abundance.")

    json.dump(out, open(os.path.join(C.RESULTS, f"sigma1_{tag}.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("interpretation", "testable_subunits")}, indent=2))
    return out


def main():
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    run("GSE172460", "rupp", "Lepr_cluster",
        relabel=lambda g: inv.get(g, ""), keep=lambda g: g in set(C.RUPP_POPS.values()))
    run("GSE172204", "vmh", "VMH_cluster",
        relabel=lambda g: g.split("_")[0], keep=lambda g: bool(g))


if __name__ == "__main__":
    main()
