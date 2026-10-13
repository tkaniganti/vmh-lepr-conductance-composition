#!/usr/bin/env python3
"""
15_sigma1_group_scan.py - what covaries with sigma1A share when the unit of analysis is the
population rather than the cell?

14 showed the per-cell route is dead: the median cell with the AP1-small position detected
carries exactly ONE UMI across all three paralogs, so its per-cell share is 1 or 0 by chance.
Every per-cell sigma1 analysis is therefore uninformative, and the apparent anti-correlation
between paralogs is arithmetic, not co-regulation.

Pseudobulk is the measurement that survives, because counts are summed over hundreds of cells
before the ratio is taken. That is what the gradient in 12 rests on. This script makes the
group the unit: for every deposited grouping it computes the sigma1A share and the expression
of every gene, then correlates them across groups.

  Affinati GSE172204: 24 VMH clusters (and the 6 conserved classes, reported separately)
  Rupp GSE172460:     all Lepr_cluster populations, not just the three used earlier

Two filters, because n is 13 to 24 groups and per-dataset FDR alone would be either empty or
credulous:
  * a permutation null, shuffling the share across groups, which says how many hits of a
    given |rho| to expect by chance
  * replication, keeping only genes that agree in sign and rank highly in BOTH datasets,
    which is a far stronger filter than significance in either

Group median depth and group size are carried as covariates, since both could drive a share.

Output: results/sigma1_groupscan_<tag>.csv, sigma1_groupscan.json
"""
import os, csv, json, math
from collections import defaultdict
import numpy as np
from scipy import stats
import config as C
import vmhlib as V

SIGMA = ["Ap1s1", "Ap1s2", "Ap1s3"]
MIN_GROUP_CELLS = 40
MIN_GROUPS_DETECTING = 0.8     # gene must be expressed in this fraction of groups
N_PERM = 2000
TOP_N = 300


def extract_block(prefix, tag, cols):
    out = V.dpath(f"block_{tag}.txt")
    colfile = V.dpath(f"blockcols_{tag}.txt")
    with open(colfile, "w") as fh:
        fh.write("\n".join(str(c) for c in sorted(cols)) + "\n")
    if not os.path.exists(out):
        V.sh("zcat " + V.dpath(f"{prefix}_matrix.mtx.gz") + " | awk '"
             "NR==FNR{keep[$1]=1; next} /^%/{next} !dims{dims=1; next} "
             "{ if($2 in keep) print $1, $2, $3 }' " + colfile + " - > " + out)
    return out


def partial_spearman_groups(x, y, z):
    """Spearman between x and y with z (group median depth) regressed out of the ranks."""
    zr = stats.rankdata(z).astype(float)
    def r(v):
        rv = stats.rankdata(v).astype(float)
        b = np.polyfit(zr, rv, 1)
        return rv - (b[0] * zr + b[1])
    a, b_ = r(x), r(y)
    if a.std() <= 0 or b_.std() <= 0:
        return 0.0
    return float(np.dot((a - a.mean()) / a.std(), (b_ - b_.mean()) / b_.std()) / a.size)


def build(prefix, tag, label_col, relabel=lambda g: g):
    feats, name2row = V.load_features(prefix)
    bcs, bc2col = V.load_barcodes(prefix)
    meta = V.load_metadata(prefix)
    sym = {v: k for k, v in name2row.items()}
    col2group, _, _ = V.write_col2group(prefix, meta, bc2col, label_col,
                                        f"{tag}_grp", keep=lambda g: bool(g))
    cs, rs = V.pass_totals(prefix, col2group, f"{tag}_grp")
    colsum = V.read_two_col(cs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = relabel(g)
    cols = sorted(c for c in colsum if col2cl.get(c))
    block = extract_block(prefix, f"{tag}_grp_blk", set(cols))

    gsum = defaultdict(lambda: defaultdict(float))      # group -> row -> counts
    gtot = defaultdict(float)
    gcells = defaultdict(list)
    for c in cols:
        gtot[col2cl[c]] += colsum[c]
        gcells[col2cl[c]].append(colsum[c])
    ndet = defaultdict(set)
    with open(block) as fh:
        for line in fh:
            r, c, v = line.split()
            c = int(c)
            g = col2cl.get(c)
            if g is None:
                continue
            val = float(v)
            if val > 0:
                gsum[g][int(r)] += val
                ndet[int(r)].add(g)

    groups = sorted([g for g in gtot if len(gcells[g]) >= MIN_GROUP_CELLS])
    srows = [name2row[g] for g in SIGMA if g in name2row]
    share = []
    for g in groups:
        tot = sum(gsum[g].get(r, 0.0) for r in srows)
        share.append(gsum[g].get(srows[0], 0.0) / tot if tot > 0 else np.nan)
    share = np.array(share)
    depth = np.array([float(np.median(gcells[g])) for g in groups])
    size = np.array([len(gcells[g]) for g in groups], float)
    expr = {}
    need = int(MIN_GROUPS_DETECTING * len(groups))
    for r, gs in ndet.items():
        if r in srows or len(gs & set(groups)) < need:
            continue
        expr[r] = np.array([math.log1p(gsum[g].get(r, 0.0) / gtot[g] * C.NORM_TARGET)
                            for g in groups])
    return {"groups": groups, "share": share, "depth": depth, "size": size,
            "expr": expr, "sym": sym, "n_cells": {g: len(gcells[g]) for g in groups}}


def scan(d, tag):
    groups, share, depth, size = d["groups"], d["share"], d["depth"], d["size"]
    ok = ~np.isnan(share)
    share, depth, size = share[ok], depth[ok], size[ok]
    groups = [g for g, k in zip(groups, ok) if k]
    print(f"   {len(groups)} groups, share {share.min():.3f} to {share.max():.3f}")
    print(f"   share vs group median depth: rho={stats.spearmanr(depth, share).statistic:+.3f}"
          f"   vs group size: rho={stats.spearmanr(size, share).statistic:+.3f}")

    rows = []
    for r, v in d["expr"].items():
        v = v[ok]
        rows.append({"gene": d["sym"].get(r, str(r)), "row": r,
                     "rho": round(partial_spearman_groups(v, share, depth), 4)})
    obs = np.array([abs(x["rho"]) for x in rows])

    rng = np.random.default_rng(C.RANDOM_SEED)
    nulls = []
    idx = np.arange(share.size)
    for _ in range(min(N_PERM, 500)):
        rng.shuffle(idx)
        sp = share[idx]
        vals = [abs(partial_spearman_groups(d["expr"][r][ok], sp, depth))
                for r in list(d["expr"])[:400]]
        nulls.append(np.percentile(vals, 99))
    thresh = float(np.mean(nulls))
    rows.sort(key=lambda x: -abs(x["rho"]))
    for x in rows:
        x["exceeds_perm_99th"] = bool(abs(x["rho"]) > thresh)
    with open(os.path.join(C.RESULTS, f"sigma1_groupscan_{tag}.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "row", "rho", "exceeds_perm_99th"])
        w.writeheader(); w.writerows(rows[:TOP_N])
    return {"tag": tag, "n_groups": len(groups), "groups": groups,
            "share_by_group": {g: round(float(s), 4) for g, s in zip(groups, share)},
            "n_cells_by_group": {g: d["n_cells"][g] for g in groups},
            "share_vs_depth_rho": round(float(stats.spearmanr(depth, share).statistic), 4),
            "share_vs_size_rho": round(float(stats.spearmanr(size, share).statistic), 4),
            "perm_99th_threshold": round(thresh, 4),
            "n_genes": len(rows),
            "n_exceeding": int(sum(1 for x in rows if x["exceeds_perm_99th"])),
            "top30": [{"gene": x["gene"], "rho": x["rho"]} for x in rows[:30]],
            "_rows": {x["gene"]: x["rho"] for x in rows}}


def main():
    out = {}
    print("=== Affinati, 24 VMH clusters ===")
    d1 = build("GSE172204", "vmh24", "VMH_cluster")
    out["vmh24"] = scan(d1, "vmh24")

    print("\n=== Rupp, all Lepr_cluster populations ===")
    d2 = build("GSE172460", "rupp13", "Lepr_cluster")
    out["rupp13"] = scan(d2, "rupp13")

    # ------------------------------------------------------------ replication filter
    a, b = out["vmh24"].pop("_rows"), out["rupp13"].pop("_rows")
    shared = set(a) & set(b)
    rep = [{"gene": g, "rho_vmh": a[g], "rho_rupp": b[g],
            "same_sign": bool(np.sign(a[g]) == np.sign(b[g]) and a[g] != 0),
            "min_abs": round(min(abs(a[g]), abs(b[g])), 4)}
           for g in shared]
    rep = [r for r in rep if r["same_sign"]]
    rep.sort(key=lambda r: -r["min_abs"])
    with open(os.path.join(C.RESULTS, "sigma1_groupscan_replicated.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "rho_vmh", "rho_rupp", "same_sign", "min_abs"])
        w.writeheader(); w.writerows(rep[:TOP_N])
    out["replication"] = {
        "n_genes_in_both": len(shared),
        "n_same_sign": len(rep),
        "expected_same_sign_by_chance": round(len(shared) * 0.5, 1),
        "top30": [{"gene": r["gene"], "rho_vmh": r["rho_vmh"], "rho_rupp": r["rho_rupp"]}
                  for r in rep[:30]]}

    json.dump(out, open(os.path.join(C.RESULTS, "sigma1_groupscan.json"), "w"), indent=2)
    print("\n== replication ==")
    print(json.dumps(out["replication"], indent=2))
    print("\nshare by group (Affinati):", out["vmh24"]["share_by_group"])
    print("share by group (Rupp):", out["rupp13"]["share_by_group"])


if __name__ == "__main__":
    main()
