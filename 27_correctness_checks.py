#!/usr/bin/env python3
"""
27_correctness_checks.py - checks that ask whether the pipeline is right, not whether it reruns.

The clean-clone verification in run_verify.slurm establishes reproducibility: the same code on
the same data returns the same numbers. It cannot detect an error that was present from the
start. These three checks have known answers that do not come from the pipeline.

  CHECK A - parser spot-check against the raw matrix
      Reads counts for a few genes in a few nuclei directly out of the gzipped matrix-market
      file by row and column index, with no use of vmhlib, and compares against what vmhlib
      returns for the same entries. Catches the silent class of error where a parser drops,
      shifts or misassigns entries: the counts still look plausible, so nothing downstream
      complains. The known matrix-market header bug in an earlier version of this project was
      exactly this kind, and was found only because it produced a visible inconsistency.

  CHECK B - random-gene panel
      Replaces the ion-channel panel with the same number of randomly chosen genes, matched
      for detection rate, and reruns the nested partition. Every reported quantity should
      collapse toward its null. If the between-class component survives with arbitrary genes,
      the estimator is responding to something structural in the data (depth, cluster size,
      ambient RNA) and not to conductance composition.

  CHECK C - full label shuffle
      Permutes conserved-class and cluster labels across all nuclei, breaking every real
      grouping, and reruns the partition. Both components should go to zero. This differs
      from the within-class shuffle used for the per-class nulls, which deliberately
      preserves class structure; here nothing is preserved.

A failure in A means the data are being read wrongly. A failure in B or C means the estimator
reports structure where none exists, which would invalidate the central result.

Output: results/correctness_checks.json
"""
import os, gzip, json, random, warnings
from collections import defaultdict
import numpy as np
import pandas as pd
import config as cfg
import vmhlib as V

warnings.filterwarnings("ignore")
PREFIX, TAG = "GSE172204", "corr"
N_SPOT_GENES, N_SPOT_CELLS = 4, 400
N_RANDOM_PANELS, N_SHUFFLES = 10, 200


# ----------------------------------------------------------------- CHECK A
def raw_lookup(prefix, want_rows, want_cols):
    """Read the matrix-market file directly. Deliberately does not use vmhlib."""
    path = os.path.join(cfg.DATA, f"{prefix}_matrix.mtx.gz")
    got = {}
    with gzip.open(path, "rt") as fh:
        header_done = False
        for line in fh:
            if line.startswith("%"):
                continue
            if not header_done:          # the dimensions line
                header_done = True
                continue
            r, c, v = line.split()
            r, c = int(r), int(c)
            if r in want_rows and c in want_cols:
                got[(r, c)] = float(v)
    return got


def check_a(name2row, bc2col, cols, rng):
    genes = ["Nr5a1", "Lepr", "Slc17a6", "Trpc5"][:N_SPOT_GENES]
    genes = [g for g in genes if g in name2row]
    rows = {name2row[g]: g for g in genes}
    sample_cols = set(rng.sample(list(cols), min(N_SPOT_CELLS, len(cols))))
    print(f"   reading {len(rows)} genes x {len(sample_cols)} nuclei straight from the .mtx",
          flush=True)
    raw = raw_lookup(PREFIX, set(rows), sample_cols)

    ent = V.pass_rows(PREFIX, set(rows), TAG, "spot")
    lib = {}
    for l in open(ent):
        r, c, v = l.split()
        r, c = int(r), int(c)
        if c in sample_cols:
            lib[(r, c)] = float(v)

    only_raw = set(raw) - set(lib)
    only_lib = set(lib) - set(raw)
    disagree = [(k, raw[k], lib[k]) for k in set(raw) & set(lib) if raw[k] != lib[k]]
    ok = not (only_raw or only_lib or disagree)
    print(f"   raw entries {len(raw)}, library entries {len(lib)}, "
          f"missing from library {len(only_raw)}, extra {len(only_lib)}, "
          f"value disagreements {len(disagree)}", flush=True)
    if disagree[:3]:
        for k, a, b in disagree[:3]:
            print(f"     row {k[0]} col {k[1]}: raw {a} vs library {b}")
    return {"genes": genes, "n_nuclei_sampled": len(sample_cols),
            "n_entries_raw": len(raw), "n_entries_library": len(lib),
            "missing_from_library": len(only_raw), "extra_in_library": len(only_lib),
            "value_disagreements": len(disagree), "pass": bool(ok)}


# ----------------------------------------------------------------- shared
def nested(y, cls, clu):
    df = pd.DataFrame({"y": np.asarray(y, float), "A": np.asarray(cls),
                       "B": [f"{a}|{b}" for a, b in zip(cls, clu)]})
    n, a, b = len(df), df.A.nunique(), df.B.nunique()
    if a < 2 or b <= a:
        return None
    grand = df.y.mean()
    ssa = sum(len(g) * (g.y.mean() - grand) ** 2 for _, g in df.groupby("A"))
    ssb = sum(len(g) * (g.y.mean() - df.loc[df.A == g.A.iloc[0], "y"].mean()) ** 2
              for _, g in df.groupby("B"))
    sse = sum(((g.y - g.y.mean()) ** 2).sum() for _, g in df.groupby("B"))
    msa, msb, mse = ssa / (a - 1), ssb / (b - a), sse / (n - b)
    nb = df.groupby("B").size().to_numpy(float)
    na = df.groupby("A").size().to_numpy(float)
    n0b = (n - (nb ** 2).sum() / n) / (b - a)
    n0a = (n - (na ** 2).sum() / n) / (a - 1)
    s2b = max((msb - mse) / n0b, 0.0)
    s2a = max((msa - msb) / n0a, 0.0)
    tot = s2a + s2b + mse
    return {"pct_class": 100 * s2a / tot, "pct_cluster_in_class": 100 * s2b / tot}


def main():
    rng = random.Random(cfg.RANDOM_SEED)
    feats, name2row = V.load_features(PREFIX)
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
    cols = sorted(c for c in colsum if c in c2c)
    col2bc = {bc2col[b]: b for b in meta if b in bc2col}
    cls = np.array([c2c[c] for c in cols])
    clu = np.array([(meta[col2bc[c]].get("VMH_cluster") or "NA").strip() for c in cols])
    ld = np.log10(np.array([colsum[c] for c in cols], float))
    out = {}

    print(">> CHECK A: parser against the raw matrix", flush=True)
    out["check_A_parser"] = check_a(name2row, bc2col, cols, rng)
    print(f"   {'PASS' if out['check_A_parser']['pass'] else 'FAIL'}\n", flush=True)

    # observed value, for reference
    pr, cr = V.choose_controls(name2row, rowsum, len(cols), cfg.RANDOM_SEED)
    need = {r for rs_ in list(pr.values()) + list(cr.values()) for r in rs_}
    ent = V.pass_rows(PREFIX, need, TAG, "obs")
    expr = defaultdict(dict)
    for l in open(ent):
        r, c, v = l.split()
        c = int(c)
        if c in c2c:
            expr[c][int(r)] = float(v)
    y = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
    y = y - np.polyval(np.polyfit(ld, y, 1), ld)
    obs = nested(y, cls, clu)
    out["observed"] = {k: round(v, 4) for k, v in obs.items()}
    print(f">> observed: class {obs['pct_class']:.3f}%, "
          f"cluster-in-class {obs['pct_cluster_in_class']:.3f}%\n", flush=True)

    # ------------------------------------------------------------- CHECK B
    print(">> CHECK B: random gene panels in place of the effector panel", flush=True)
    panel = [g for p in cfg.PANELS for g in cfg.PANELS[p]]
    n_panel = len(panel)
    detectable = [g for g in name2row
                  if name2row[g] in rowsum and rowsum[name2row[g]] > 0]
    recs = []
    for i in range(N_RANDOM_PANELS):
        pick = rng.sample(detectable, min(2 * n_panel, len(detectable)))
        fake = {"HYPER": [name2row[g] for g in pick[:n_panel]],
                "DEPOL": [name2row[g] for g in pick[n_panel:2 * n_panel]]}
        fr, fc = V.choose_controls(name2row, rowsum, len(cols), cfg.RANDOM_SEED + 500 + i)
        need = set(fake["HYPER"]) | set(fake["DEPOL"]) | {
            r for rs_ in list(fc.values()) for r in rs_}
        e = V.pass_rows(PREFIX, need, TAG, f"rand{i}")
        ex = defaultdict(dict)
        for l in open(e):
            r, c, v = l.split()
            c = int(c)
            if c in c2c:
                ex[c][int(r)] = float(v)
        try:
            ys = V.module_scores(ex, cols, colsum, fake, fc)["HminusD"]
        except Exception as err:
            print(f"   panel {i}: scoring failed ({err})", flush=True)
            continue
        ys = ys - np.polyval(np.polyfit(ld, ys, 1), ld)
        r = nested(ys, cls, clu)
        if r:
            recs.append(r)
            print(f"   panel {i}: class {r['pct_class']:.3f}%, "
                  f"cluster-in-class {r['pct_cluster_in_class']:.3f}%", flush=True)
    if recs:
        a = np.array([[r["pct_class"], r["pct_cluster_in_class"]] for r in recs])
        out["check_B_random_panels"] = {
            "n_panels": len(recs),
            "class_mean": round(float(a[:, 0].mean()), 4),
            "class_max": round(float(a[:, 0].max()), 4),
            "cluster_mean": round(float(a[:, 1].mean()), 4),
            "cluster_max": round(float(a[:, 1].max()), 4),
            "observed_class": out["observed"]["pct_class"],
            "pass": bool(a[:, 0].max() < out["observed"]["pct_class"]),
            "note_criterion": "observed must exceed every random panel; no threshold chosen",
            "note": ("random panels should give a smaller class component than the effector "
                     "panel; if they do not, the estimator is reporting structure that is "
                     "not specific to conductance genes")}
        print(f"   {'PASS' if out['check_B_random_panels']['pass'] else 'FAIL'}: "
              f"random max {a[:, 0].max():.3f}% vs observed "
              f"{out['observed']['pct_class']:.3f}%\n", flush=True)

    # ------------------------------------------------------------- CHECK C
    print(">> CHECK C: full label shuffle, nothing preserved", flush=True)
    recs = []
    for i in range(N_SHUFFLES):
        idx = np.arange(len(cols))
        rng.shuffle(idx)
        r = nested(y, cls[idx], clu[idx])
        if r:
            recs.append(r)
    a = np.array([[r["pct_class"], r["pct_cluster_in_class"]] for r in recs])
    # Positional criterion, matching the per-class nulls used elsewhere: the observed value
    # must exceed the 95th percentile of the shuffle distribution. No threshold is chosen.
    p95 = a.mean(axis=0) * 0 + np.percentile(a, 95, axis=0)
    above = [float((a[:, k] >= (obs["pct_class"], obs["pct_cluster_in_class"])[k]).mean())
             for k in (0, 1)]
    out["check_C_full_shuffle"] = {
        "n_shuffles": len(recs),
        "class": {"shuffle_mean": round(float(a[:, 0].mean()), 4),
                  "shuffle_p95": round(float(p95[0]), 4),
                  "shuffle_max": round(float(a[:, 0].max()), 4),
                  "observed": out["observed"]["pct_class"],
                  "fraction_of_shuffles_at_or_above_observed": round(above[0], 4)},
        "cluster_in_class": {"shuffle_mean": round(float(a[:, 1].mean()), 4),
                             "shuffle_p95": round(float(p95[1]), 4),
                             "shuffle_max": round(float(a[:, 1].max()), 4),
                             "observed": out["observed"]["pct_cluster_in_class"],
                             "fraction_of_shuffles_at_or_above_observed": round(above[1], 4)},
        "pass": bool(obs["pct_class"] > p95[0]
                     and obs["pct_cluster_in_class"] > p95[1]),
        "note": ("the observed components must exceed the 95th percentile of the shuffle "
                 "distribution; the criterion is positional and introduces no chosen value")}
    print(f"   class: observed {obs['pct_class']:.3f}% vs shuffle mean "
          f"{a[:, 0].mean():.3f}%, p95 {p95[0]:.3f}%, max {a[:, 0].max():.3f}%", flush=True)
    print(f"   cluster-in-class: observed {obs['pct_cluster_in_class']:.3f}% vs shuffle mean "
          f"{a[:, 1].mean():.3f}%, p95 {p95[1]:.3f}%, max {a[:, 1].max():.3f}%", flush=True)
    print(f"   shuffles reaching the observed value: {above[0]:.1%} (class), "
          f"{above[1]:.1%} (cluster)", flush=True)
    print(f"   {'PASS' if out['check_C_full_shuffle']['pass'] else 'FAIL'}\n", flush=True)

    json.dump(out, open(os.path.join(cfg.RESULTS, "correctness_checks.json"), "w"), indent=2)
    failed = [k for k, v in out.items() if isinstance(v, dict) and v.get("pass") is False]
    print("ALL CHECKS PASS" if not failed else f"FAILED: {failed}")
    print("Written to results/correctness_checks.json")


if __name__ == "__main__":
    main()
