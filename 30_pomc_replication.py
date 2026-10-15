#!/usr/bin/env python3
"""
30_pomc_replication.py - does the within-population null replicate in POMC/LepR neurons?

The VMH result rests on recordings from one laboratory's doctoral work. Biglari et al. (2021)
report the same three-way pattern in the POMC/LepR intersection with published, intersectionally
targeted electrophysiology: leptin excited 40%, inhibited 20% and left 40% of these neurons
unaffected. The Rupp reporter-defined atlas (GSE172460) contains that same intersection as a
molecularly defined population, so the within-population analysis can be repeated there.

This is a replication attempt, not a comparison. The POMC population in that atlas holds a few
hundred cells against 1,850 in the VMH comparison, so any null here is weaker by construction
and the detectability bound must be reported alongside it, not inferred from the VMH bound.

What runs:

  0. Inventory first. Print what cluster labels the deposited metadata actually provides for
     the POMC population. A nested partition needs two levels; if no substructure is deposited,
     the partition is skipped rather than invented, since clustering de novo and then testing
     for structure in the same cells is circular.

  1. Discreteness. Gaussian mixture BIC and Hartigan's dip on the per-cell effector-balance
     score, as for the VMH population. Needs no sub-labels.

  2. Detectability at this n. Two-component mixtures injected at known separations, to state
     what a split of a given size would have shown at this sample size. The VMH analysis
     reported that a two-standard-deviation split was recoverable and a one-standard-deviation
     split was not; the equivalent numbers here will be worse and are computed, not assumed.

  3. Within-population scan. Partial correlation of every adequately detected transcript with
     the balance score, residualising ranks on the rank of log depth, with BH correction. The
     VMH scan returned no transcript above about 2% of variance.

  4. Nested partition, only if the metadata provide sub-labels for this population.

Reading the result: a null here replicates the VMH finding in a population whose heterogeneity
is independently published, which strengthens it. Structure here would mean the VMH null is
population-specific, and the smaller n would then need addressing directly, since the weaker
dataset having found something the stronger did not is the harder result to interpret.

Output: results/pomc_replication.json
"""
import os, json, warnings
from collections import defaultdict, Counter
import numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172460", "pomc"
POP_FIELD, POP_MATCH = "Lepr_cluster", "POMC"
SPLITS = [0.5, 1.0, 1.5, 2.0, 3.0]
N_SIM = 100
MIN_DETECT_PCT = 10.0


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = np.empty_like(p)
    m, prev = len(p), 1.0
    for rank, i in enumerate(o[::-1]):
        prev = min(prev, p[i] * m / (m - rank))
        q[i] = prev
    return q


def main():
    feats, name2row = V.load_features(PREFIX)
    row2name = {v: k for k, v in name2row.items()}
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    out = {}

    # ---------------------------------------------------------------- 0
    print(">> inventory of deposited labels", flush=True)
    fields = set()
    for m in list(meta.values())[:50]:
        fields |= set(m.keys())
    print(f"   metadata fields: {sorted(fields)}", flush=True)
    pop_bcs = [b for b, m in meta.items()
               if POP_MATCH.lower() in (m.get(POP_FIELD) or "").lower()]
    if not pop_bcs:
        raise SystemExit(f"no cells matched {POP_FIELD} ~ {POP_MATCH}; "
                         f"check the field name against the printed list")
    labels = Counter((meta[b].get(POP_FIELD) or "").strip() for b in pop_bcs)
    print(f"   {len(pop_bcs)} cells matching '{POP_MATCH}' in {POP_FIELD}")
    for k, v in labels.most_common():
        print(f"     {k:<28} {v}")
    sub_field = None
    for f in sorted(fields):
        if f in (POP_FIELD,):
            continue
        vals = Counter((meta[b].get(f) or "").strip() for b in pop_bcs)
        vals.pop("", None)
        if 2 <= len(vals) <= 20:
            print(f"   candidate sub-label field '{f}': {len(vals)} values "
                  f"{dict(list(vals.most_common(6)))}")
            if sub_field is None and len(vals) >= 2:
                sub_field = f
    out["inventory"] = {"n_cells": len(pop_bcs),
                        "population_labels": dict(labels),
                        "candidate_sub_label_field": sub_field}

    cols = sorted(bc2col[b] for b in pop_bcs if b in bc2col)
    col2bc = {bc2col[b]: b for b in pop_bcs if b in bc2col}
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, POP_FIELD, TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    cols = [c for c in cols if c in colsum]
    ld = np.log10(np.array([colsum[c] for c in cols], float))
    print(f"\n>> {len(cols)} cells with counts; median depth "
          f"{np.median([colsum[c] for c in cols]):.0f}\n", flush=True)

    # module scores over the control-gene draws, as elsewhere
    seeds = list(range(cfg.RANDOM_SEED, cfg.RANDOM_SEED + cfg.N_SEEDS))
    seed_sets, union = [], set()
    for s in seeds:
        pr, cr = V.choose_controls(name2row, rowsum, len(cols), s)
        seed_sets.append((pr, cr))
        union |= {r for rs_ in list(pr.values()) + list(cr.values()) for r in rs_}
    ent = V.pass_rows(PREFIX, union, TAG, "union")
    expr = defaultdict(dict)
    colset = set(cols)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in colset:
            expr[c][int(r)] = float(v)
    scores = []
    for pr, cr in seed_sets:
        y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
        scores.append(y - np.polyval(np.polyfit(ld, y, 1), ld))
    y0 = scores[0]

    # ---------------------------------------------------------------- 1
    print(">> TEST 1: discreteness of the effector-balance score", flush=True)
    x = y0.reshape(-1, 1)
    bics = {}
    for k in (1, 2, 3, 4):
        try:
            bics[k] = float(GaussianMixture(k, random_state=cfg.RANDOM_SEED,
                                            n_init=5).fit(x).bic(x))
        except Exception:
            pass
    best = min(bics, key=bics.get)
    try:
        import diptest
        dstat, dp = diptest.diptest(y0)
    except Exception:
        dstat, dp = None, None
    print(f"   BIC {({k: round(v,1) for k,v in bics.items()})}, best k = {best}")
    print(f"   dip test D = {dstat}, p = {dp}")
    out["test1_discreteness"] = {"bic": {str(k): round(v, 2) for k, v in bics.items()},
                                 "best_k": best, "dip_D": dstat, "dip_p": dp}

    # ---------------------------------------------------------------- 2
    print("\n>> TEST 2: detectability at this sample size", flush=True)
    rng = np.random.default_rng(cfg.RANDOM_SEED)
    sd = float(np.std(y0))
    lev = []
    for d in SPLITS:
        hits = 0
        for _ in range(N_SIM):
            half = len(y0) // 2
            sim = y0.copy()
            idx = rng.permutation(len(sim))[:half]
            sim[idx] += d * sd
            try:
                b = {k: GaussianMixture(k, random_state=0, n_init=3)
                     .fit(sim.reshape(-1, 1)).bic(sim.reshape(-1, 1)) for k in (1, 2)}
                hits += b[2] < b[1]
            except Exception:
                pass
        lev.append({"split_sd": d, "pct_detected": round(100 * hits / N_SIM, 1)})
        print(f"   split {d:.1f} SD: detected in {lev[-1]['pct_detected']:.0f}% of simulations")
    out["test2_detectability"] = {"n_cells": len(cols), "levels": lev,
                                  "note": "equivalent VMH figures were 0% at 1 SD and 100% at 2 SD "
                                          "with 1,850 nuclei; this sample is smaller"}

    # ---------------------------------------------------------------- 3
    print("\n>> TEST 3: transcripts covarying with effector balance", flush=True)
    panel_rows = {r for p in cfg.PANELS for g in cfg.PANELS[p] if g in name2row
                  for r in [name2row[g]]}
    ctrl_rows = {r for _, cr in seed_sets for rs_ in cr.values() for r in rs_}
    excluded = panel_rows | ctrl_rows
    keep = {r for r, t in rowsum.items() if t > 0} - excluded
    ent2 = V.pass_rows(PREFIX, keep, TAG, "scan")
    vals = defaultdict(dict)
    for l in open(ent2):
        r, c, v = l.split()
        c = int(c)
        if c in colset:
            vals[int(r)][c] = float(v)
    ry = stats.rankdata(y0)
    rd = stats.rankdata(ld)
    res_y = ry - np.polyval(np.polyfit(rd, ry, 1), rd)
    recs = []
    for r, d in vals.items():
        g = row2name.get(r)
        if not g:
            continue
        v = np.array([d.get(c, 0.0) for c in cols])
        if 100 * (v > 0).mean() < MIN_DETECT_PCT:
            continue
        rv = stats.rankdata(np.log1p(v / np.array([colsum[c] for c in cols]) * 1e4))
        res_v = rv - np.polyval(np.polyfit(rd, rv, 1), rd)
        rho, p = stats.pearsonr(res_y, res_v)
        recs.append({"gene": g, "rho": float(rho), "p": float(p)})
    if recs:
        q = bh([r["p"] for r in recs])
        for r, qq in zip(recs, q):
            r["q"] = float(qq)
        recs.sort(key=lambda r: -abs(r["rho"]))
        sig = [r for r in recs if r["q"] < 0.05]
        print(f"   {len(recs)} transcripts tested, {len(sig)} at q < 0.05")
        print(f"   largest |partial rho| = {abs(recs[0]['rho']):.3f} ({recs[0]['gene']}), "
              f"variance explained {100*recs[0]['rho']**2:.1f}%")
        for r in recs[:8]:
            print(f"     {r['gene']:<12} rho {r['rho']:+.3f}  q {r['q']:.3g}")
        out["test3_scan"] = {"n_tested": len(recs), "n_q05": len(sig),
                             "max_abs_rho": round(abs(recs[0]["rho"]), 4),
                             "max_pct_variance": round(100 * recs[0]["rho"] ** 2, 2),
                             "top": [{k: (round(v, 4) if isinstance(v, float) else v)
                                      for k, v in r.items()} for r in recs[:20]]}

    # ---------------------------------------------------------------- 4
    if sub_field:
        print(f"\n>> TEST 4: partition using deposited sub-labels '{sub_field}'", flush=True)
        sub = np.array([(meta[col2bc[c]].get(sub_field) or "NA").strip() for c in cols])
        keep_m = sub != "NA"
        k = len(set(sub[keep_m]))
        if k >= 2:
            df_y = y0[keep_m]
            grp = sub[keep_m]
            n = len(df_y)
            gm = df_y.mean()
            ssb = sum(((grp == u).sum()) * (df_y[grp == u].mean() - gm) ** 2
                      for u in set(grp))
            ssw = sum(((df_y[grp == u] - df_y[grp == u].mean()) ** 2).sum() for u in set(grp))
            msb, msw = ssb / (k - 1), ssw / (n - k)
            ni = np.array([(grp == u).sum() for u in set(grp)], float)
            n0 = (n - (ni ** 2).sum() / n) / (k - 1)
            sb = max((msb - msw) / n0, 0.0)
            pct = 100 * sb / (sb + msw)
            print(f"   {k} sub-groups, {n} cells: between-group variance {pct:.3f}%")
            out["test4_partition"] = {"field": sub_field, "k": k, "n": n,
                                      "pct_between": round(pct, 4)}
    else:
        print("\n>> TEST 4 skipped: no deposited sub-labels for this population. "
              "Clustering de novo and then testing the same cells would be circular.")
        out["test4_partition"] = None

    json.dump(out, open(os.path.join(cfg.RESULTS, "pomc_replication.json"), "w"), indent=2)
    print("\nWritten to results/pomc_replication.json")


if __name__ == "__main__":
    main()
