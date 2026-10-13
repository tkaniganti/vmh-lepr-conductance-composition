"""vmhlib.py - shared loaders and scoring used by 02-06.

Design: the 10x matrix is never loaded whole. Single awk passes over the gzipped
matrix-market file extract only what is needed, and intermediates are cached in data/.
Peak memory stays near 1-2 GB, so RM-shared with -c 8 is ample.
"""
import os, gzip, csv, json, math, random, subprocess
from collections import defaultdict
import numpy as np
import config as C


def sh(cmd):
    subprocess.run(cmd, shell=True, check=True, executable="/bin/bash")


def dpath(*p):
    return os.path.join(C.DATA, *p)


# ---------------------------------------------------------------- annotations
def load_features(prefix):
    feats = [l.strip().strip('"') for l in gzip.open(dpath(f"{prefix}_features.tsv.gz"), "rt")]
    # features.tsv may be 1 or 3 columns; take the gene symbol column if tabbed
    feats = [f.split("\t")[1] if "\t" in f else f for f in feats]
    return feats, {g: i + 1 for i, g in enumerate(feats)}


def load_barcodes(prefix):
    bcs = [l.strip().strip('"') for l in gzip.open(dpath(f"{prefix}_barcodes.tsv.gz"), "rt")]
    return bcs, {bc: i + 1 for i, bc in enumerate(bcs)}


def load_metadata(prefix):
    rows = list(csv.DictReader(gzip.open(dpath(f"{prefix}_metadata.csv.gz"), "rt")))
    return {r["ID"]: r for r in rows}


def write_col2group(prefix, meta, bc2col, group_col, tag, keep=None):
    """Map matrix column -> group label for cells with a non-empty group_col.

    Returns (path, n_labelled, n_unmatched). n_unmatched is the count of labelled
    barcodes absent from barcodes.tsv; it is reported rather than silently dropped,
    which is where the 6048/6049 and 2913/2914 discrepancies came from.
    """
    f = dpath(f"col2group_{tag}.txt")
    n_lab = n_unmatched = 0
    lines = []
    for bc, m in meta.items():
        g = (m.get(group_col) or "").strip()
        if not g:
            continue
        if keep and not keep(g):
            continue
        n_lab += 1
        c = bc2col.get(bc)
        if c is None:
            n_unmatched += 1
            continue
        lines.append(f"{c}\t{g}")
    with open(f, "w") as fh:
        fh.write("\n".join(lines) + "\n")
    return f, n_lab, n_unmatched


# ---------------------------------------------------------------- matrix passes
def pass_totals(prefix, col2group, tag):
    """Per-column library size and per-row totals, restricted to the labelled columns."""
    cs, rs = dpath(f"colsum_{tag}.txt"), dpath(f"rowsum_{tag}.txt")
    if os.path.exists(cs) and os.path.exists(rs):
        return cs, rs
    sh(f"zcat {dpath(prefix + '_matrix.mtx.gz')} | awk 'NR==FNR{{cl[$1]=$2; next}} "
       f"/^%/{{next}} !dims{{dims=1; next}} {{ if($2 in cl){{colsum[$2]+=$3; rowsum[$1]+=$3}} }} "
       f"END{{ for(c in colsum) print c,colsum[c] > \"{cs}\"; "
       f"for(r in rowsum) print r,rowsum[r] > \"{rs}\" }}' {col2group} -")
    return cs, rs


def pass_rows(prefix, rows, tag, name):
    """Extract all matrix entries for a set of gene rows."""
    out = dpath(f"entries_{tag}_{name}.txt")
    if os.path.exists(out):
        return out
    rf = dpath(f"rows_{tag}_{name}.txt")
    with open(rf, "w") as fh:
        fh.write("\n".join(str(r) for r in sorted(rows)) + "\n")
    sh(f"zcat {dpath(prefix + '_matrix.mtx.gz')} | awk 'NR==FNR{{need[$1]=1; next}} "
       f"/^%/{{next}} !dims{{dims=1; next}} {{ if($1 in need) print $1,$2,$3 }}' {rf} - > {out}")
    return out


def read_two_col(path, cast=float):
    d = {}
    for l in open(path):
        a, b = l.split()
        d[int(a)] = cast(b)
    return d


# ---------------------------------------------------------------- module scores
def choose_controls(name2row, rowsum, n_cells, seed):
    """AddModuleScore-style control genes matched on binned mean expression."""
    rng = random.Random(seed)
    mean = {r: rowsum.get(r, 0.0) / n_cells for r in range(1, len(name2row) + 1)}
    expressed = sorted([r for r in mean if mean[r] > 0], key=lambda r: mean[r])
    binsz = max(1, len(expressed) // C.N_EXPR_BINS)
    bin_of = {r: min(i // binsz, C.N_EXPR_BINS - 1) for i, r in enumerate(expressed)}
    members = defaultdict(list)
    for r in expressed:
        members[bin_of[r]].append(r)
    panel_rows = {p: [name2row[g] for g in gs
                      if g in name2row and mean.get(name2row[g], 0) > 0]
                  for p, gs in C.PANELS.items() if p in ("HYPER", "DEPOL")}
    allpanel = {r for rs in panel_rows.values() for r in rs}
    ctrl = {}
    for p, rs in panel_rows.items():
        # equalise control-set size across panels rather than letting it follow panel size
        per = max(1, int(round(C.N_CTRL_GENES / len(rs))))
        s = set()
        for r in rs:
            cand = [x for x in members[bin_of[r]] if x not in allpanel and x not in s]
            rng.shuffle(cand)
            s.update(cand[:per])
        s = sorted(s)
        rng.shuffle(s)
        ctrl[p] = sorted(s[:C.N_CTRL_GENES])
    return panel_rows, ctrl


def module_scores(expr, cols, colsum, panel_rows, ctrl_rows):
    """Control-corrected per-cell module scores. expr: {col: {row: raw count}}."""
    allrows = sorted({r for rs in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs})
    ln = {}
    for r in allrows:
        v = np.fromiter((math.log1p(expr.get(c, {}).get(r, 0.0) / colsum[c] * C.NORM_TARGET)
                         for c in cols), dtype=float, count=len(cols))
        sd = v.std()
        ln[r] = (v - v.mean()) / (sd if sd > 0 else 1.0)
    out = {}
    for p in panel_rows:
        pan = np.mean([ln[r] for r in panel_rows[p]], axis=0)
        ctl = np.mean([ln[r] for r in ctrl_rows[p]], axis=0)
        out[p] = pan - ctl
    out["HminusD"] = out["HYPER"] - out["DEPOL"]
    return out


def eps2_from_kruskal(H, n_groups, n_total):
    return (H - n_groups + 1) / (n_total - n_groups)


def boot_ci(x, n_boot=None, seed=0):
    n_boot = n_boot or C.N_BOOT
    rng = np.random.default_rng(seed)
    x = np.asarray(x, dtype=float)
    means = rng.choice(x, size=(n_boot, x.size), replace=True).mean(axis=1)
    return float(x.mean()), float(np.percentile(means, 2.5)), float(np.percentile(means, 97.5))


def perm_test(a, b, n_perm=None, seed=0):
    """Two-sided label-permutation test on the difference in means."""
    n_perm = n_perm or C.N_PERM
    rng = np.random.default_rng(seed)
    a, b = np.asarray(a, float), np.asarray(b, float)
    obs = a.mean() - b.mean()
    pool = np.concatenate([a, b])
    na = a.size
    hits = 0
    for _ in range(n_perm):
        rng.shuffle(pool)
        if abs(pool[:na].mean() - pool[na:].mean()) >= abs(obs):
            hits += 1
    return float(obs), (hits + 1) / (n_perm + 1)
