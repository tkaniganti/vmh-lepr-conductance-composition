#!/usr/bin/env python3
"""
08_gene_scan.py - two genome-wide scans over the 2,913 LepR nuclei.

Scan A, cluster structure given depth.
  For every adequately detected gene, the likelihood-ratio statistic for cluster identity
  over and above log10 library size, the same test 03_trpc5_depth.py applies to Trpc5.
  This places Trpc5 and the effector panels in the genome-wide ranking rather than against
  200 arbitrarily matched controls, and shows what actually separates the LepR clusters.

Scan B, covariation with the effector axis within the population.
  Partial Spearman correlation between each gene and the per-cell H-D score, with log depth
  regressed out of both. Panel and control genes are excluded, so the scan is not circular.
  This is the only analysis here that addresses within-population heterogeneity, which is
  what the recordings sample. It is hypothesis-generating: anything it surfaces needs the
  depth control and an orthogonal measurement before it enters the manuscript.

Inputs:  data/GSE172204_matrix.mtx.gz, data/col2group_vmh.txt, data/percell_scores.npz
Outputs: results/gene_scan_cluster.csv, gene_scan_axis.csv, gene_scan.json
"""
import os, json, csv, math, gzip
import numpy as np
from scipy import stats
import config as C
import vmhlib as V

PREFIX = "GSE172204"
MIN_DETECT_FRAC = 0.05      # genes detected in fewer than 5% of LepR nuclei are skipped
MAX_DETECT_FRAC = 0.95
TOP_N = 300                 # rows written per scan
L2_C = 1e6                  # effectively unpenalised logistic regression


# ---------------------------------------------------------------- logistic helpers
def _loglik(y, p):
    p = np.clip(p, 1e-12, 1 - 1e-12)
    return float(np.sum(y * np.log(p) + (1 - y) * np.log(1 - p)))


def _fit(X, y):
    from sklearn.linear_model import LogisticRegression
    try:
        m = LogisticRegression(C=np.inf, solver="lbfgs", max_iter=300)
        m.fit(X, y)
    except (TypeError, ValueError):
        m = LogisticRegression(C=L2_C, solver="lbfgs", max_iter=300)
        m.fit(X, y)
    return _loglik(y, m.predict_proba(X)[:, 1])


def lr_given_depth(y, Xdepth, Xfull, df):
    ll0 = _fit(Xdepth, y)
    ll1 = _fit(Xfull, y)
    lr = 2 * (ll1 - ll0)
    return float(lr), float(stats.chi2.sf(max(lr, 0.0), df))


# ---------------------------------------------------------------- matrix pass
def extract_lepr_block(lepr_cols):
    """One pass over the matrix, keeping only LepR columns. Header is detected, not assumed."""
    out = V.dpath("lepr_block.txt")
    colfile = V.dpath("lepr_cols.txt")
    with open(colfile, "w") as fh:
        fh.write("\n".join(str(c) for c in sorted(lepr_cols)) + "\n")
    if not os.path.exists(out):
        # skip every %-comment line, then the single dimensions line
        V.sh(
            "zcat " + V.dpath(f"{PREFIX}_matrix.mtx.gz") + " | awk '"
            "NR==FNR{keep[$1]=1; next} "
            "/^%/{next} "
            "!dims{dims=1; next} "
            "{ if($2 in keep) print $1, $2, $3 }' " + colfile + " - > " + out)
    return out


def count_header_lines():
    n_comment = 0
    with gzip.open(V.dpath(f"{PREFIX}_matrix.mtx.gz"), "rt") as fh:
        for line in fh:
            if line.startswith("%"):
                n_comment += 1
            else:
                return n_comment + 1, line.strip()
    return n_comment + 1, ""


def main():
    feats, name2row = V.load_features(PREFIX)
    sym = {v: k for k, v in name2row.items()}

    n_header, dims_line = count_header_lines()
    print(f">> matrix header lines: {n_header}  (dims: {dims_line})")
    # header length is detected here and in vmhlib.py; nothing assumes a fixed count

    z = np.load(V.dpath("percell_scores.npz"), allow_pickle=True)
    cols = np.asarray(z["cols"], dtype=np.int64)
    labels = np.asarray(z["labels"])
    hd = np.asarray(z["HminusD"], float)
    col_index = {c: i for i, c in enumerate(cols)}
    n = cols.size
    print(f">> {n} LepR nuclei, {len(set(labels.tolist()))} clusters")

    colsum = V.read_two_col(V.dpath("colsum_vmh.txt"))
    depth = np.array([colsum[c] for c in cols], float)
    logdepth = np.log10(depth)

    block = extract_lepr_block(set(cols.tolist()))

    # ------------------------------------------------- build detection and expression
    # Two passes over the extracted block: first count nonzeros per gene to decide which
    # genes to keep, then fill dense arrays only for those. Keeps peak memory modest.
    nz = {}
    with open(block) as fh:
        for line in fh:
            r, c, v = line.split()
            if float(v) > 0:
                nz[int(r)] = nz.get(int(r), 0) + 1
    lo, hi = MIN_DETECT_FRAC * n, MAX_DETECT_FRAC * n
    keep = sorted([r for r, k in nz.items() if lo <= k <= hi])
    kidx = {r: i for i, r in enumerate(keep)}
    print(f">> {len(nz)} genes with any signal; {len(keep)} kept "
          f"({MIN_DETECT_FRAC:.0%} to {MAX_DETECT_FRAC:.0%} detection)")

    det = np.zeros((len(keep), n), dtype=np.int8)
    expr = np.zeros((len(keep), n), dtype=np.float32)
    with open(block) as fh:
        for line in fh:
            r, c, v = line.split()
            r = int(r)
            i = kidx.get(r)
            if i is None:
                continue
            j = col_index.get(int(c))
            if j is None:
                continue
            val = float(v)
            if val > 0:
                det[i, j] = 1
                expr[i, j] = math.log1p(val / colsum[int(c)] * C.NORM_TARGET)
    print(f">> matrices built: det {det.shape}, expr {expr.shape}")

    # ------------------------------------------------------------------- scan A
    import pandas as pd
    D = pd.get_dummies(labels, drop_first=True).to_numpy(dtype=float)
    Xd = np.column_stack([np.ones(n), logdepth])
    Xf = np.column_stack([Xd, D])
    df_cluster = D.shape[1]

    rowsA = []
    for i, r in enumerate(keep):
        y = det[i].astype(float)
        try:
            lr, p = lr_given_depth(y, Xd, Xf, df_cluster)
        except Exception:
            continue
        rowsA.append({"gene": sym.get(r, str(r)), "row": r,
                      "detect_pct": round(100 * float(y.mean()), 1),
                      "LR_cluster_given_depth": round(lr, 2), "p": p})
        if (i + 1) % 1000 == 0:
            print(f"   scan A {i + 1}/{len(keep)}")
    rowsA.sort(key=lambda d: -d["LR_cluster_given_depth"])
    pA = np.array([d["p"] for d in rowsA])
    qA = stats.false_discovery_control(pA, method="bh") if hasattr(
        stats, "false_discovery_control") else np.full(pA.size, np.nan)
    for d, q in zip(rowsA, qA):
        d["q_bh"] = float(q)
    for rank, d in enumerate(rowsA, 1):
        d["rank"] = rank

    # ------------------------------------------------------------------- scan B
    panel_rows, ctrl_rows = V.choose_controls(
        name2row, V.read_two_col(V.dpath("rowsum_vmh.txt")), n, C.RANDOM_SEED)
    excluded = {r for rs in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs}
    excluded |= {name2row[g] for p in C.PANELS for g in C.PANELS[p] if g in name2row}

    def resid_rank(x):
        rx = stats.rankdata(x)
        b = np.polyfit(logdepth, rx, 1)
        return rx - (b[0] * logdepth + b[1])

    y_res = resid_rank(hd)
    y_res = (y_res - y_res.mean()) / (y_res.std() + 1e-12)

    rowsB = []
    for i, r in enumerate(keep):
        if r in excluded:
            continue
        x_res = resid_rank(expr[i])
        sd = x_res.std()
        if sd <= 0:
            continue
        x_res = (x_res - x_res.mean()) / sd
        rho = float(np.dot(x_res, y_res) / n)
        t = rho * math.sqrt(max(n - 3, 1) / max(1 - rho ** 2, 1e-12))
        p = float(2 * stats.t.sf(abs(t), n - 3))
        rowsB.append({"gene": sym.get(r, str(r)), "row": r,
                      "partial_rho": round(rho, 4), "p": p,
                      "detect_pct": round(100 * float(det[i].mean()), 1)})
        if (i + 1) % 2000 == 0:
            print(f"   scan B {i + 1}/{len(keep)}")
    pB = np.array([d["p"] for d in rowsB])
    qB = stats.false_discovery_control(pB, method="bh") if hasattr(
        stats, "false_discovery_control") else np.full(pB.size, np.nan)
    for d, q in zip(rowsB, qB):
        d["q_bh"] = float(q)
    rowsB.sort(key=lambda d: -abs(d["partial_rho"]))

    # ------------------------------------------------------------------- output
    def where(rows, gene, key):
        for i, d in enumerate(rows, 1):
            if d["gene"] == gene:
                return {"rank": i, "of": len(rows), "percentile": round(100 * (1 - i / len(rows)), 1),
                        key: d[key], "q_bh": round(d["q_bh"], 4)}
        return None

    panel_genes = [g for p in ("HYPER", "DEPOL") for g in C.PANELS[p]]
    summary = {
        "n_nuclei": int(n),
        "n_genes_tested": len(rowsA),
        "matrix_header_lines": n_header,
        "scanA_trpc5": where(rowsA, "Trpc5", "LR_cluster_given_depth"),
        "scanA_top20": [d["gene"] for d in rowsA[:20]],
        "scanA_effector_panel_ranks": {
            g: where(rowsA, g, "LR_cluster_given_depth") for g in panel_genes},
        "scanB_n_tested": len(rowsB),
        "scanB_top20": [{"gene": d["gene"], "rho": d["partial_rho"], "q": round(d["q_bh"], 4)}
                        for d in rowsB[:20]],
        "scanB_n_sig_q05": int(sum(1 for d in rowsB if d["q_bh"] < 0.05)),
        "note": ("Scan B is hypothesis-generating. A hit must survive the depth control in "
                 "03 and an orthogonal measurement before it enters the manuscript."),
    }

    with open(os.path.join(C.RESULTS, "gene_scan_cluster.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["rank", "gene", "row", "detect_pct",
                                           "LR_cluster_given_depth", "p", "q_bh"])
        w.writeheader(); w.writerows(rowsA[:TOP_N])
    with open(os.path.join(C.RESULTS, "gene_scan_axis.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "row", "partial_rho", "p", "q_bh",
                                           "detect_pct"])
        w.writeheader(); w.writerows(rowsB[:TOP_N])
    json.dump(summary, open(os.path.join(C.RESULTS, "gene_scan.json"), "w"), indent=2)
    print("\n" + json.dumps({k: summary[k] for k in
                             ("scanA_trpc5", "scanA_top20", "scanB_top20",
                              "scanB_n_sig_q05")}, indent=2))
    print("\nWritten to", C.RESULTS)


if __name__ == "__main__":
    main()
