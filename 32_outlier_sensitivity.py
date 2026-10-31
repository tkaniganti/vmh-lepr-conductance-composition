#!/usr/bin/env python3
"""
32_outlier_sensitivity.py - do rare-gene outliers drive the per-nucleus results?

The worked example (script 31) showed that a single UMI in a low-abundance panel gene can
produce an extreme per-nucleus score. Kcnj11 is detected in 3.1% of VMH nuclei, so its mean
log-CP10k across nuclei is 0.0174 with a standard deviation of 0.1033; one count gives a
z-score above +8, and a nucleus carrying it scored +1.99 on the effector axis against a
cluster-level range of -0.129 to 0.247.

Averaging over thirteen panel genes and then over nuclei in a cluster dilutes this, which is
why the pseudobulk values are stable. The per-nucleus analyses have no such dilution. This
asks whether the two of them that matter are affected:

  2.4  discreteness. Gaussian mixture BIC and Hartigan's dip on the per-nucleus axis. A
       heavy right tail from rare-gene counts could in principle create apparent structure,
       or mask it by widening the distribution.

  2.5  the within-population transcript scan. Partial correlations between each transcript
       and the axis. A handful of extreme axis values gives a few nuclei large leverage.

Four conditions, each rerunning both tests on the leptin-receptor population:

  as reported      the per-nucleus axis as used in the manuscript
  trimmed          the most extreme TRIM_PCT% of nuclei by |axis| removed
  winsorised       axis values clipped to the 1st and 99th percentiles
  rare genes out   panel genes detected in fewer than MIN_DETECT% of nuclei dropped from
                   both modules and the axis recomputed, which removes the cause rather
                   than the symptom

If the reported conclusions hold in all four, the outliers are a property of the data that
does not drive the result. If they do not, the manuscript needs to say which condition it
is reporting and why.

Output: results/outlier_sensitivity.json
"""
import os, json, warnings
from collections import defaultdict
import numpy as np
from scipy import stats
from sklearn.mixture import GaussianMixture
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "outl"
TRIM_PCT = 1.0          # per tail
MIN_DETECT = 10.0       # panel genes below this detection rate are "rare"
SCAN_MIN_DETECT = 10.0  # transcripts below this are not scanned
N_SIM_NULL = 200


def discreteness(y, seed):
    x = np.asarray(y, float).reshape(-1, 1)
    bic = {}
    for k in (1, 2, 3, 4):
        try:
            bic[k] = float(GaussianMixture(k, random_state=seed, n_init=5).fit(x).bic(x))
        except Exception:
            pass
    best = min(bic, key=bic.get) if bic else None
    try:
        import diptest
        D, p = diptest.diptest(np.asarray(y, float))
    except Exception:
        D, p = None, None
    return {"bic": {str(k): round(v, 2) for k, v in bic.items()},
            "best_k": best, "dip_D": D, "dip_p": p}


def main():
    feats, name2row = V.load_features(PREFIX)
    row2name = {v: k for k, v in name2row.items()}
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_grp", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    allcols = sorted(c for c in colsum if c in col2bc
                     and (meta[col2bc[c]].get("VMH_cluster") or "").strip())
    clu = {c: meta[col2bc[c]]["VMH_cluster"].strip() for c in allcols}
    cols = [c for c in allcols if clu[c].startswith("Lepr")]
    print(f">> {len(cols)} leptin-receptor nuclei of {len(allcols)} VMH", flush=True)

    panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(allcols), cfg.RANDOM_SEED)
    need = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs_}
    ent = V.pass_rows(PREFIX, need, TAG, "union")
    expr = defaultdict(dict)
    colset = set(allcols)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in colset:
            expr[c][int(r)] = float(v)

    # detection rate of each panel gene, which is what makes a gene "rare"
    det = {}
    for mod in ("HYPER", "DEPOL"):
        for r in panel_rows[mod]:
            hits = sum(1 for c in allcols if expr.get(c, {}).get(r, 0.0) > 0)
            det[r] = 100 * hits / len(allcols)
    rare = {r for r, d in det.items() if d < MIN_DETECT}
    print(f">> panel genes below {MIN_DETECT}% detection: "
          f"{sorted(row2name.get(r, r) for r in rare)}", flush=True)
    for r in sorted(det, key=lambda r: det[r]):
        print(f"   {row2name.get(r, r):<9} {det[r]:>5.1f}%"
              f"{'   RARE' if r in rare else ''}", flush=True)

    base = V.module_scores(expr, allcols, colsum, panel_rows, ctrl_rows)["HminusD"]
    idx = {c: i for i, c in enumerate(allcols)}
    y0 = np.array([base[idx[c]] for c in cols])

    # panel with rare genes removed
    pr2 = {m: [r for r in panel_rows[m] if r not in rare] for m in panel_rows}
    dropped = {m: len(panel_rows[m]) - len(pr2[m]) for m in panel_rows}
    y_norare = None
    if all(len(v) >= 3 for v in pr2.values()):
        s2 = V.module_scores(expr, allcols, colsum, pr2, ctrl_rows)["HminusD"]
        y_norare = np.array([s2[idx[c]] for c in cols])

    lo, hi = np.percentile(y0, [TRIM_PCT, 100 - TRIM_PCT])
    keep = (y0 >= lo) & (y0 <= hi)
    conditions = {
        "as_reported": y0,
        "trimmed": y0[keep],
        "winsorised": np.clip(y0, lo, hi),
    }
    if y_norare is not None:
        conditions["rare_genes_removed"] = y_norare

    out = {"n_nuclei": len(cols),
           "panel_detection_pct": {row2name.get(r, str(r)): round(d, 1)
                                   for r, d in det.items()},
           "rare_genes": sorted(row2name.get(r, str(r)) for r in rare),
           "n_dropped_per_module": dropped,
           "trim_pct_per_tail": TRIM_PCT,
           "axis_summary": {}, "discreteness": {}, "scan": {}}

    print("\n>> per-nucleus axis, by condition", flush=True)
    for lab, y in conditions.items():
        out["axis_summary"][lab] = {
            "n": int(len(y)), "mean": round(float(y.mean()), 4),
            "sd": round(float(y.std()), 4),
            "min": round(float(y.min()), 4), "max": round(float(y.max()), 4),
            "pct_above_1": round(100 * float((np.abs(y) > 1).mean()), 2)}
        s = out["axis_summary"][lab]
        print(f"   {lab:<20} n={s['n']:<5} sd={s['sd']:<7} "
              f"range {s['min']:+.3f} to {s['max']:+.3f}   "
              f"{s['pct_above_1']}% beyond |1|", flush=True)

    print("\n>> TEST 2.4: discreteness, by condition", flush=True)
    for lab, y in conditions.items():
        d = discreteness(y, cfg.RANDOM_SEED)
        out["discreteness"][lab] = d
        print(f"   {lab:<20} best k = {d['best_k']}   "
              f"dip D = {d['dip_D']}, p = {d['dip_p']}", flush=True)

    # ---------------------------------------------------------------- 2.5 scan
    print("\n>> TEST 2.5: strongest transcript correlate, by condition", flush=True)
    keep_rows = {r for r, t in rowsum.items() if t > 0}
    panel_and_ctrl = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values())
                      for r in rs_}
    scan_rows = keep_rows - panel_and_ctrl
    ent2 = V.pass_rows(PREFIX, scan_rows, TAG, "scan")
    vals = defaultdict(dict)
    cset = set(cols)
    for l in open(ent2):
        r, c, v = l.split()
        c = int(c)
        if c in cset:
            vals[int(r)][c] = float(v)
    depth = np.array([colsum[c] for c in cols], float)
    rd = stats.rankdata(np.log10(depth))

    for lab, y in conditions.items():
        if len(y) != len(cols):
            out["scan"][lab] = {"skipped": "trimming changes n; scan run on full set only"}
            print(f"   {lab:<20} skipped (n changed)", flush=True)
            continue
        ry = stats.rankdata(y)
        res_y = ry - np.polyval(np.polyfit(rd, ry, 1), rd)
        best = None
        for r, d in vals.items():
            g = row2name.get(r)
            if not g:
                continue
            v = np.array([d.get(c, 0.0) for c in cols])
            if 100 * (v > 0).mean() < SCAN_MIN_DETECT:
                continue
            rv = stats.rankdata(np.log1p(v / depth * cfg.NORM_TARGET))
            res_v = rv - np.polyval(np.polyfit(rd, rv, 1), rd)
            rho, p = stats.pearsonr(res_y, res_v)
            if best is None or abs(rho) > abs(best[1]):
                best = (g, float(rho), float(p))
        out["scan"][lab] = {"top_gene": best[0], "partial_rho": round(best[1], 4),
                            "p": best[2],
                            "pct_variance": round(100 * best[1] ** 2, 2)} if best else None
        if best:
            print(f"   {lab:<20} {best[0]:<12} rho {best[1]:+.3f}  "
                  f"({100*best[1]**2:.1f}% of variance)", flush=True)

    # ---------------------------------------------------------------- verdict
    ks = {lab: out["discreteness"][lab]["best_k"] for lab in conditions}
    same_k = len(set(ks.values())) == 1
    rhos = [v["partial_rho"] for v in out["scan"].values()
            if isinstance(v, dict) and "partial_rho" in v]
    out["verdict"] = {
        "discreteness_unchanged": bool(same_k),
        "best_k_by_condition": ks,
        "max_abs_partial_rho_across_conditions": round(max(abs(r) for r in rhos), 4)
        if rhos else None,
        "note": ("if best k is the same in every condition and no condition produces a "
                 "transcript correlate of appreciable size, the rare-gene outliers are a "
                 "property of the data that does not drive the reported conclusions")}
    print(f"\n>> best k by condition: {ks}")
    print(f">> discreteness conclusion unchanged across conditions: {same_k}")

    json.dump(out, open(os.path.join(cfg.RESULTS, "outlier_sensitivity.json"), "w"), indent=2)
    print("\nWritten to results/outlier_sensitivity.json")


if __name__ == "__main__":
    main()
