#!/usr/bin/env python3
"""
28_esr1_1_markers.py - what distinguishes Esr1_1 from the rest of its conserved class?

The per-class survey found that within-class variance in effector composition reaches 15.65%
in the conserved Esr1 class and collapses to 0.38% when Esr1_1 is removed, while removal of
any other cluster leaves it between 16.5% and 19.5%. The separation survives thresholding and
grows when small clusters are excluded, so it is not an artifact of sparse clusters. Nothing
in the VMH manuscript addresses what underlies it.

This is a two-group comparison inside one conserved class: Esr1_1 against the other Esr1_*
clusters in that class. Restricting to one class controls for everything that differs between
classes, so a gene that comes out here differs within a group the atlas treats as one type.

Three things are reported:

  effector panel first, since the separation was defined on conductance composition. If the
    15.65% is carried by one or two channels, those are named here.

  genome-wide, by Mann-Whitney on log-CP10k with Benjamini-Hochberg correction, restricted to
    genes detected in at least MIN_DETECT_PCT of nuclei in either group so that the test is
    not dominated by near-zero transcripts.

  depth as a confound. Clusters differ in capture depth, so median depth per group is printed
    and the genome-wide test is repeated on depth-matched subsamples. A gene that survives
    only the unmatched test is reported as depth-sensitive rather than as a marker.

Ambient-RNA artifacts are flagged, not removed: Gm42418 and similar rRNA-overlapping loci
appear at the top of many such tests and were the strongest hit in the within-Lepr scan.

This is exploratory and is not part of the VMH manuscript.

Output: results/esr1_1_markers.csv, results/esr1_1_markers.json
"""
import os, csv, json, warnings
from collections import defaultdict
import numpy as np
from scipy import stats
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "e1"
FOCUS_CLASS, FOCUS_CLUSTER = "Esr1", "Esr1_1"
MIN_DETECT_PCT = 10.0
N_MATCH_DRAWS = 20
AMBIENT_FLAGS = ("Gm42418", "Gm26917", "AY036118", "Lars2", "Malat1")


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m = len(p)
    prev = 1.0
    for rank, i in enumerate(o[::-1]):
        k = m - rank
        prev = min(prev, p[i] * m / k)
        q[i] = prev
    return q


def main():
    feats, name2row = V.load_features(PREFIX)
    row2name = {v: k for k, v in name2row.items()}
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        raw = (m.get("Conserved_cluster") or "").strip()
        m["Conserved_class"] = ("" if (not raw or raw.startswith("Doublets"))
                                else raw.split("_")[0])
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "Conserved_class", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    c2c = {}
    for l in open(cg):
        a, b = l.split()
        c2c[int(a)] = b
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cols = [c for c in sorted(colsum) if c2c.get(c) == FOCUS_CLASS]
    clu = {c: (meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols}

    # group A: the focus cluster. group B: the other Esr1_* clusters in this class.
    A = [c for c in cols if clu[c] == FOCUS_CLUSTER]
    B = [c for c in cols if clu[c].startswith(FOCUS_CLASS + "_") and clu[c] != FOCUS_CLUSTER]
    if len(A) < 20 or len(B) < 20:
        raise SystemExit(f"too few nuclei: {FOCUS_CLUSTER} n={len(A)}, others n={len(B)}")
    dA = np.array([colsum[c] for c in A], float)
    dB = np.array([colsum[c] for c in B], float)
    print(f">> {FOCUS_CLUSTER}: {len(A)} nuclei, median depth {np.median(dA):.0f}", flush=True)
    print(f">> other {FOCUS_CLASS} clusters ({len(set(clu[c] for c in B))}): "
          f"{len(B)} nuclei, median depth {np.median(dB):.0f}", flush=True)
    print(f"   depth ratio {np.median(dA)/np.median(dB):.2f}", flush=True)

    keep = {r for r, t in rowsum.items() if t > 0}
    ent = V.pass_rows(PREFIX, keep, TAG, "all")
    expr = defaultdict(dict)
    colset = set(A) | set(B)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in colset:
            expr[int(r)][c] = float(v)

    def cp10k(row, cells, depths):
        d = expr.get(row, {})
        return np.log1p(np.array([d.get(c, 0.0) for c in cells]) / depths * 1e4)

    # depth-matched subsamples, drawn once and reused for every gene
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    lo, hi = np.percentile(np.concatenate([dA, dB]), [10, 90])
    poolA = [i for i, d in enumerate(dA) if lo <= d <= hi]
    poolB = [i for i, d in enumerate(dB) if lo <= d <= hi]
    nm = min(len(poolA), len(poolB))
    matched = [(rng.choice(poolA, nm, replace=False), rng.choice(poolB, nm, replace=False))
               for _ in range(N_MATCH_DRAWS)]
    print(f">> depth-matched comparisons: {nm} nuclei per group x {N_MATCH_DRAWS} draws\n",
          flush=True)

    rows = []
    panel = {g for p in cfg.PANELS for g in cfg.PANELS[p]}
    for r in sorted(keep):
        g = row2name.get(r)
        if not g:
            continue
        d = expr.get(r, {})
        detA = 100 * sum(1 for c in A if d.get(c, 0) > 0) / len(A)
        detB = 100 * sum(1 for c in B if d.get(c, 0) > 0) / len(B)
        if max(detA, detB) < MIN_DETECT_PCT:
            continue
        xa, xb = cp10k(r, A, dA), cp10k(r, B, dB)
        try:
            u, p = stats.mannwhitneyu(xa, xb, alternative="two-sided")
        except ValueError:
            continue
        pm = []
        for ia, ib in matched:
            try:
                pm.append(stats.mannwhitneyu(xa[ia], xb[ib], alternative="two-sided")[1])
            except ValueError:
                pass
        rows.append({"gene": g, "in_effector_panel": g in panel,
                     "mean_log_cp10k_focus": round(float(xa.mean()), 4),
                     "mean_log_cp10k_others": round(float(xb.mean()), 4),
                     "diff": round(float(xa.mean() - xb.mean()), 4),
                     "detect_pct_focus": round(detA, 1), "detect_pct_others": round(detB, 1),
                     "p": p,
                     "p_depth_matched_median": (round(float(np.median(pm)), 6) if pm else None),
                     "ambient_flag": g in AMBIENT_FLAGS})
    q = bh([r["p"] for r in rows])
    for r, qq in zip(rows, q):
        r["q"] = round(float(qq), 6)
        r["survives_depth_matching"] = (r["p_depth_matched_median"] is not None
                                        and r["p_depth_matched_median"] < 0.05)
    rows.sort(key=lambda r: (-abs(r["diff"]), r["q"]))
    print(f">> {len(rows)} genes tested at >= {MIN_DETECT_PCT}% detection in either group")
    sig = [r for r in rows if r["q"] < 0.05]
    dm = [r for r in sig if r["survives_depth_matching"]]
    print(f">> {len(sig)} at q < 0.05, {len(dm)} of those surviving depth matching\n")

    print(">> effector panel genes:")
    for r in [x for x in rows if x["in_effector_panel"]][:12]:
        print(f"   {r['gene']:<9} diff {r['diff']:+.3f}  q {r['q']:.4g}  "
              f"depth-matched {'yes' if r['survives_depth_matching'] else 'no'}")

    print("\n>> largest differences surviving correction and depth matching:")
    for r in dm[:20]:
        flag = "  [ambient?]" if r["ambient_flag"] else ""
        print(f"   {r['gene']:<12} diff {r['diff']:+.3f}  q {r['q']:.3g}  "
              f"detect {r['detect_pct_focus']:.0f}% vs {r['detect_pct_others']:.0f}%{flag}")

    with open(os.path.join(cfg.RESULTS, "esr1_1_markers.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)
    json.dump({"focus": FOCUS_CLUSTER, "focus_n": len(A), "others_n": len(B),
               "other_clusters": sorted(set(clu[c] for c in B)),
               "median_depth_focus": float(np.median(dA)),
               "median_depth_others": float(np.median(dB)),
               "n_tested": len(rows), "n_q05": len(sig), "n_depth_matched": len(dm),
               "top_depth_matched": [r["gene"] for r in dm[:30]],
               "effector_panel_hits": [r["gene"] for r in dm if r["in_effector_panel"]]},
              open(os.path.join(cfg.RESULTS, "esr1_1_markers.json"), "w"), indent=2)
    print("\nWritten to results/esr1_1_markers.csv and .json")


if __name__ == "__main__":
    main()
