#!/usr/bin/env python3
"""
09_dispersion.py - is the VMH SF1/LepR population more heterogeneous in effector balance
than arcuate POMC and AgRP?

The placement analysis in 06 compares means. The question the recordings pose is about
spread: a population whose effector balance is broadly distributed should produce mixed
leptin response signs, while narrowly distributed populations should produce consistent
ones. POMC and AgRP have literature ground truth for consistent responses, so they are the
comparison this test needs.

Statistic: Brown-Forsythe (Levene centred on the median) on the per-cell H-D score, with
four controls, because every obvious confound inflates spread:

  1. depth       - residualize on log10 library size, and separately quantile-match depth
                   across populations, since capture noise scales with depth
  2. cell number - downsample all populations to the smallest (n = 189) over many draws,
                   since a larger sample explores more of its own distribution
  3. substructure - repeat within deposited neuron subclusters, since a population that is
                   two cell types is dispersed for an uninteresting reason
  4. gene choice  - an empirical null of random expression-matched pseudo-panels, since any
                   particular gene set can produce a spread difference by chance

Everything is repeated over the N_SEEDS control-gene draws, and sign consistency across
draws is the robustness criterion, as in 06.

Inputs:  cached Rupp intermediates from 06 (no new matrix pass)
Output:  results/dispersion.json, results/dispersion_percell.csv
"""
import os, csv, json, math
from collections import defaultdict
import numpy as np
from scipy import stats
import config as C
import vmhlib as V

PREFIX, TAG = "GSE172460", "rupp"
LABEL_COL = "Lepr_cluster"
SUBCLUSTER_COL = "Neuron_cluster"
N_DOWNSAMPLE = 500
N_NULL_PANELS = 200
MIN_SUBCLUSTER_N = 25


def brown_forsythe(groups):
    """Levene's test centred on the median."""
    dev = [np.abs(g - np.median(g)) for g in groups]
    F, p = stats.f_oneway(*dev)
    return float(F), float(p)


def resid_on_depth(x, logdepth):
    b = np.polyfit(logdepth, x, 1)
    return x - (b[0] * logdepth + b[1])


def quantile_match(pop_idx, depth, n_bins=10):
    """Equal-size sample from each population within each pooled depth decile."""
    edges = np.quantile(depth, np.linspace(0, 1, n_bins + 1))
    edges[-1] += 1e-9
    keep = {p: [] for p in pop_idx}
    for b in range(n_bins):
        inbin = {p: [i for i in idx if edges[b] <= depth[i] < edges[b + 1]]
                 for p, idx in pop_idx.items()}
        k = min(len(v) for v in inbin.values())
        if k == 0:
            continue
        rng = np.random.default_rng(1000 + b)
        for p, v in inbin.items():
            keep[p].extend(rng.choice(v, size=k, replace=False).tolist())
    return keep


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    wanted = set(C.RUPP_POPS.values())
    col2group, _, _ = V.write_col2group(
        PREFIX, meta, bc2col, LABEL_COL, TAG, keep=lambda g: g in wanted)
    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = g
    cols = sorted(colsum)
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    name = np.array([inv[col2cl[c]] for c in cols])
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)

    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    sub = np.array([(meta[col2bc[c]].get(SUBCLUSTER_COL) or "NA").strip() for c in cols])

    pops = ["POMC", "AgRP", "VMH_SF1"]
    pop_idx = {p: np.where(name == p)[0] for p in pops}
    print(">> n per population:", {p: len(v) for p, v in pop_idx.items()})
    print(">> median depth:", {p: int(np.median(depth[v])) for p, v in pop_idx.items()})

    # ---------------------------------------------------- scores, reusing cached entries
    expr = defaultdict(dict)
    per_seed, ctrl_pool = [], set()
    for seed in range(C.RANDOM_SEED, C.RANDOM_SEED + C.N_SEEDS):
        panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(cols), seed)
        need = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs_}
        ctrl_pool |= {r for rs_ in ctrl_rows.values() for r in rs_}
        e = V.pass_rows(PREFIX, need, TAG, f"seed{seed}")   # cached from 06
        for l in open(e):
            r, c, v = l.split()
            c = int(c)
            if c in col2cl:
                expr[c][int(r)] = float(v)
        per_seed.append((panel_rows, ctrl_rows,
                         V.module_scores(expr, cols, colsum, panel_rows, ctrl_rows)))
        print(f"   scored seed {seed}")

    def hd(i):
        return per_seed[i][2]["HminusD"]

    # ---------------------------------------------------------------- 1. depth residual
    out = {"n_per_population": {p: int(len(v)) for p, v in pop_idx.items()},
           "median_depth": {p: float(np.median(depth[v])) for p, v in pop_idx.items()},
           "n_seeds": C.N_SEEDS}

    sds, bf, ratios = defaultdict(list), [], defaultdict(list)
    for i in range(C.N_SEEDS):
        r = resid_on_depth(hd(i), logdepth)
        groups = [r[pop_idx[p]] for p in pops]
        for p, g in zip(pops, groups):
            sds[p].append(float(g.std(ddof=1)))
        F, pv = brown_forsythe(groups)
        bf.append({"F": round(F, 3), "p": pv})
        ratios["VMH_over_POMC"].append(sds["VMH_SF1"][-1] / sds["POMC"][-1])
        ratios["VMH_over_AgRP"].append(sds["VMH_SF1"][-1] / sds["AgRP"][-1])

    def summ(v):
        v = np.asarray(v, float)
        return {"mean": round(float(v.mean()), 4), "min": round(float(v.min()), 4),
                "max": round(float(v.max()), 4)}

    out["depth_residualized"] = {
        "sd_by_population": {p: summ(v) for p, v in sds.items()},
        "brown_forsythe_F": summ([b["F"] for b in bf]),
        "brown_forsythe_p_max": max(b["p"] for b in bf),
        "sd_ratio": {k: summ(v) for k, v in ratios.items()},
        "vmh_more_dispersed_than_both_in_all_seeds": bool(
            all(a > 1 for a in ratios["VMH_over_POMC"]) and
            all(a > 1 for a in ratios["VMH_over_AgRP"])),
    }

    # ---------------------------------------------------- 2. equal n, 3. depth matched
    rng = np.random.default_rng(C.RANDOM_SEED)
    nmin = min(len(v) for v in pop_idx.values())
    r0 = resid_on_depth(hd(0), logdepth)
    ds_ratio = {"VMH_over_POMC": [], "VMH_over_AgRP": []}
    ds_bf_p = []
    for _ in range(N_DOWNSAMPLE):
        pick = {p: rng.choice(pop_idx[p], size=nmin, replace=False) for p in pops}
        g = {p: r0[pick[p]] for p in pops}
        ds_ratio["VMH_over_POMC"].append(g["VMH_SF1"].std(ddof=1) / g["POMC"].std(ddof=1))
        ds_ratio["VMH_over_AgRP"].append(g["VMH_SF1"].std(ddof=1) / g["AgRP"].std(ddof=1))
        ds_bf_p.append(brown_forsythe([g[p] for p in pops])[1])
    out["equal_n"] = {
        "n_per_population": int(nmin), "n_draws": N_DOWNSAMPLE,
        "sd_ratio_median": {k: round(float(np.median(v)), 4) for k, v in ds_ratio.items()},
        "sd_ratio_ci": {k: [round(float(np.percentile(v, 2.5)), 4),
                            round(float(np.percentile(v, 97.5)), 4)]
                        for k, v in ds_ratio.items()},
        "pct_draws_bf_p_lt_05": round(100 * float(np.mean(np.array(ds_bf_p) < 0.05)), 1)}

    qm = quantile_match(pop_idx, depth)
    gq = {p: r0[np.array(qm[p])] for p in pops}
    Fq, pq = brown_forsythe([gq[p] for p in pops])
    out["depth_quantile_matched"] = {
        "n_per_population": int(len(qm["POMC"])),
        "sd_by_population": {p: round(float(gq[p].std(ddof=1)), 4) for p in pops},
        "brown_forsythe": {"F": round(Fq, 3), "p": pq}}

    # ------------------------------------------------------------- 4. within subcluster
    within = {}
    for p in pops:
        parts = []
        for s in sorted(set(sub[pop_idx[p]])):
            idx = pop_idx[p][sub[pop_idx[p]] == s]
            if idx.size >= MIN_SUBCLUSTER_N:
                parts.append(r0[idx] - np.median(r0[idx]))
        within[p] = np.concatenate(parts) if parts else np.array([])
        out.setdefault("subclusters_used", {})[p] = int(len(parts))
    if all(v.size > 0 for v in within.values()):
        Fw, pw = brown_forsythe([within[p] for p in pops])
        out["within_subcluster"] = {
            "sd_by_population": {p: round(float(within[p].std(ddof=1)), 4) for p in pops},
            "brown_forsythe": {"F": round(Fw, 3), "p": pw},
            "note": "each subcluster median-centred before pooling"}
    else:
        out["within_subcluster"] = "insufficient subcluster sizes"

    # ------------------------------------------------------------ 5. random-panel null
    pool = sorted(ctrl_pool)
    nH = len([r for r in per_seed[0][0]["HYPER"]])
    nD = len([r for r in per_seed[0][0]["DEPOL"]])
    null_ratio = {"VMH_over_POMC": [], "VMH_over_AgRP": []}
    rng2 = np.random.default_rng(99)
    for _ in range(N_NULL_PANELS):
        pick = rng2.permutation(pool)
        pr = {"HYPER": list(pick[:nH]), "DEPOL": list(pick[nH:nH + nD])}
        cr = {"HYPER": list(pick[nH + nD: nH + nD + 100]),
              "DEPOL": list(pick[nH + nD + 100: nH + nD + 200])}
        if len(cr["DEPOL"]) < 10:
            break
        s = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        s = resid_on_depth(s, logdepth)
        a = s[pop_idx["VMH_SF1"]].std(ddof=1)
        null_ratio["VMH_over_POMC"].append(a / s[pop_idx["POMC"]].std(ddof=1))
        null_ratio["VMH_over_AgRP"].append(a / s[pop_idx["AgRP"]].std(ddof=1))
    obs = {k: float(np.mean(v)) for k, v in ratios.items()}
    out["random_panel_null"] = {
        "n_panels": len(null_ratio["VMH_over_POMC"]),
        "pool_size": len(pool),
        "observed_ratio": {k: round(v, 4) for k, v in obs.items()},
        "null_median": {k: round(float(np.median(v)), 4) for k, v in null_ratio.items()},
        "observed_percentile": {
            k: round(100 * float(np.mean(np.asarray(null_ratio[k]) < obs[k])), 1)
            for k in null_ratio}}

    out["interpretation"] = (
        "The claim is supported only if the VMH SD ratio exceeds 1 against both arcuate "
        "populations in every control-gene draw, survives equal-n downsampling and depth "
        "quantile matching, persists within subclusters, and sits high against the "
        "random-panel null. A ratio near 1 anywhere means effector balance is no more "
        "dispersed in VMH SF1/LepR neurons and this direction should be abandoned.")

    json.dump(out, open(os.path.join(C.RESULTS, "dispersion.json"), "w"), indent=2)
    with open(os.path.join(C.RESULTS, "dispersion_percell.csv"), "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["population", "subcluster", "depth", "HminusD", "HminusD_depth_resid"])
        for i, c in enumerate(cols):
            w.writerow([name[i], sub[i], int(depth[i]), round(float(hd(0)[i]), 5),
                        round(float(r0[i]), 5)])
    print("\n" + json.dumps({k: out[k] for k in
                             ("depth_residualized", "equal_n", "depth_quantile_matched",
                              "within_subcluster", "random_panel_null")}, indent=2))


if __name__ == "__main__":
    main()
