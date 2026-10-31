#!/usr/bin/env python3
"""
49_cellstate_vs_axis.py - do per-nucleus state measures track the effector balance axis?

WHY THESE AND NOT GENES. Section 2.5 already scanned 12,182 transcripts against this axis
within this population and found nothing above about 2% of the variance. Any single gene I
could nominate has been tested and failed. What the scan could NOT see is quantities that
are not individual transcripts, and that is what this measures:

  mitochondrial fraction   share of counts from mt- genes. A standard per-cell index of
                           metabolic state, and the Discussion names metabolic state acting
                           through the K_ATP conductance as a candidate determinant. Not a
                           gene, so never scanned
  ribosomal fraction       share from Rps/Rpl genes. Translational load, and a second state
                           index that behaves differently from the mitochondrial one
  genes detected           transcriptional complexity, distinct from total depth: two nuclei
                           at equal depth can differ in how many genes carry a count
  total counts             capture depth itself, the POSITIVE CONTROL. Known to track the
                           axis at rho 0.137 in this population, so it must come back
                           tracking it, and a run where it does not is wrong

THE CONFOUND THAT MATTERS. All four correlate with depth, two of them by construction, and
the axis tracks depth. So the raw correlation is uninterpretable for every measure here and
only the partial controlling for log depth is reported as a result. Mitochondrial fraction in
particular falls as depth rises in most droplet data, so a raw negative correlation with the
axis would be expected from depth alone.

NULL. Two hundred random gene SETS of matched size, each scored as a fraction of total
counts exactly as the real sets are, giving the rate at which any such fraction correlates
with the axis. Without it, 2,913 nuclei make a partial of 0.04 look meaningful.

A POSITIVE RESULT WOULD BE THE FIRST OF THE DAY and would bear directly on the Discussion's
claim that the determinant is a state rather than an identity. A null is a fifth negative and
strengthens the same claim from the other side.

Output: results/cellstate_vs_axis.json
"""
import json
import os
import re
import warnings
from collections import defaultdict

import numpy as np
from scipy import stats

import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "cst"
N_NULL = 200


def main():
    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    axis = np.asarray(z["HminusD"], float)
    cols = [int(c) for c in z["cols"]]
    cidx = {c: i for i, c in enumerate(cols)}
    print(f">> {len(cols)} Lepr-cluster nuclei, one control-gene draw", flush=True)

    feats, name2row = V.load_features(PREFIX)
    row2name = {r: g for g, r in name2row.items()}
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    for m in meta.values():
        m["_grp"] = "VMH" if (m.get("VMH_cluster") or "").strip() else ""
    cg, _, _ = V.write_col2group(PREFIX, meta, bc2col, "_grp", TAG)
    cs, rs = V.pass_totals(PREFIX, cg, TAG)
    colsum, rowsum = V.read_two_col(cs), V.read_two_col(rs)
    depth = np.array([float(colsum.get(c, 0)) for c in cols])
    logd = np.log10(np.where(depth <= 0, 1, depth))

    mt = sorted({r for g, r in name2row.items() if re.match(r"^mt-", g, re.I)})
    rb = sorted({r for g, r in name2row.items() if re.match(r"^Rp[sl]\d", g)})
    print(f">> {len(mt)} mitochondrial genes, {len(rb)} ribosomal genes", flush=True)
    if not mt:
        print("   WARNING: no mt- genes matched; the feature list may use another "
              "prefix. Check with: grep -i '^mt' features", flush=True)

    rng = np.random.default_rng(cfg.RANDOM_SEED)
    pool = [int(r) for r in rowsum if rowsum[r] > 0]
    nullsets = [list(rng.choice(pool, max(len(mt), 5), replace=False))
                for _ in range(N_NULL)]
    need = set(mt) | set(rb) | {r for s in nullsets for r in s}
    ent = V.pass_rows(PREFIX, need, TAG, "cst")

    counts = defaultdict(lambda: np.zeros(len(cols)))
    ngenes = np.zeros(len(cols))
    colset = set(cols)
    for line in open(ent):
        r, c, v = line.split()
        c = int(c)
        if c in colset and float(v) > 0:
            counts[int(r)][cidx[c]] = float(v)

    # genes detected needs every row, not just the ones extracted, so it comes from
    # the deposited per-nucleus figure if present; otherwise it is skipped rather
    # than computed from this partial extraction, which would be wrong.
    have_ngene = False

    def partial(x):
        """Spearman of x against the axis, controlling for log depth, both ranked."""
        rx = stats.rankdata(x).astype(float)
        ra = stats.rankdata(axis).astype(float)
        rd = stats.rankdata(logd).astype(float)
        fx = rx - np.polyval(np.polyfit(rd, rx, 1), rd)
        fa = ra - np.polyval(np.polyfit(rd, ra, 1), rd)
        return float(stats.spearmanr(x, axis).statistic), \
            float(stats.spearmanr(fx, fa).statistic)

    def frac(rows):
        tot = np.zeros(len(cols))
        for r in rows:
            tot += counts[r]
        with np.errstate(divide="ignore", invalid="ignore"):
            f = np.where(depth > 0, tot / depth, np.nan)
        return f

    out = {"n_nuclei": len(cols), "draws": 1,
           "n_mt_genes": len(mt), "n_ribo_genes": len(rb),
           "caveat": ("every measure here correlates with capture depth and the axis "
                      "tracks depth at rho 0.137, so only the depth-partial is a result"),
           "measures": {}}

    print(f"\n{'measure':<24}{'rho raw':>10}{'rho|depth':>12}{'median':>12}")
    for name, vals in (("capture depth (control)", depth),
                       ("mitochondrial fraction", frac(mt) if mt else None),
                       ("ribosomal fraction", frac(rb) if rb else None)):
        if vals is None:
            print(f"{name:<24}{'gene set empty':>34}")
            continue
        ok = ~np.isnan(vals)
        raw, par = partial(np.where(ok, vals, np.nanmedian(vals)))
        out["measures"][name] = {"rho_raw": round(raw, 4),
                                 "rho_partial_depth": round(par, 4),
                                 "median": round(float(np.nanmedian(vals)), 6)}
        print(f"{name:<24}{raw:>10.4f}{par:>12.4f}{np.nanmedian(vals):>12.5f}")

    # ---- null over random gene sets of matched size
    print(f"\n>> null over {N_NULL} random gene sets of {max(len(mt), 5)} genes",
          flush=True)
    nulls = []
    for s in nullsets:
        f = frac(s)
        if np.isfinite(f).sum() < 100:
            continue
        _, par = partial(np.where(np.isfinite(f), f, 0.0))
        nulls.append(abs(par))
    p95 = float(np.percentile(nulls, 95)) if nulls else float("nan")
    out["null_abs_partial_p95"] = round(p95, 4)
    out["n_null_sets"] = len(nulls)
    print(f"   |partial| 95th percentile: {p95:.4f} ({len(nulls)} sets)")

    ctrl = out["measures"].get("capture depth (control)", {})
    ok = abs(ctrl.get("rho_raw", 0)) > 0.1
    out["positive_control_passes"] = bool(ok)
    print(f"\n   positive control, depth tracks the axis: {ok} "
          f"(rho {ctrl.get('rho_raw')})")
    if not ok:
        print("   CONTROL FAILED. Nothing below is readable.")

    beats = [k for k, v in out["measures"].items()
             if "control" not in k and abs(v["rho_partial_depth"]) > p95]
    out["measures_exceeding_null"] = beats
    print(f"\n{'=' * 68}")
    if ok and not beats:
        print("  No state measure tracks the axis once depth is controlled.")
        print("  Metabolic and translational state, as indexed here, carry no")
        print("  information about effector balance. A fifth negative, and consistent")
        print("  with the Discussion's position that the determinant is not visible in")
        print("  steady-state transcript abundance at all.")
    elif ok:
        print(f"  Exceeding the null: {', '.join(beats)}")
        print("  The first positive. Treat carefully: a fraction of counts is not a")
        print("  measurement of metabolic state, and the Discussion's claim is about")
        print("  state at the time of recording, which no fixed tissue reports.")

    json.dump(out, open(os.path.join(cfg.RESULTS, "cellstate_vs_axis.json"), "w"),
              indent=2)
    print("\nWritten to results/cellstate_vs_axis.json")


if __name__ == "__main__":
    main()
