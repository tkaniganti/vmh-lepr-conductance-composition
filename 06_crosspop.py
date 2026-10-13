#!/usr/bin/env python3
"""
06_crosspop.py - cross-population validation on the Rupp reporter-defined LepR atlas.

This is the analysis the previous bundle had no code for: the R figure scripts read
results/crosspop_scores.csv and crosspop_percell_stats.csv, and nothing generated them.

Accession: GSE172460 (verified present on GEO as a 10x trio plus metadata).
config.py's previous GSE172463 does not exist.

Population labels come from the metadata column `Lepr_cluster`. Raw label counts in that
file are L3.Nr5a1 424, L9.Pomc 231, L6.Agrp 189. The manuscript reports 347 / 230 / 189,
so this script prints the label count, the count after requiring a non-empty neuron
annotation, and the count after requiring at least one UMI, which is where the difference
is likely to sit. Record which filter the paper means and state it in Methods.

Validation gate: POMC must score more depolarizing than AgRP, and must carry higher Trpc5,
before the VMH SF1 population is placed. The script exits non-zero if that fails.

Output: results/crosspop_scores.csv, crosspop_percell_stats.csv, crosspop.json
"""
import os, csv, json, sys
from collections import defaultdict
import numpy as np
import config as C
import vmhlib as V

PREFIX, TAG = "GSE172460", "rupp"
LABEL_COL = "Lepr_cluster"


def main():
    feats, name2row = V.load_features(PREFIX)
    bcs, bc2col = V.load_barcodes(PREFIX)
    meta = V.load_metadata(PREFIX)
    if LABEL_COL not in next(iter(meta.values())):
        sys.exit(f"{LABEL_COL} not in metadata; columns = {list(next(iter(meta.values())))}")

    wanted = set(C.RUPP_POPS.values())
    col2group, n_lab, n_unmatched = V.write_col2group(
        PREFIX, meta, bc2col, LABEL_COL, TAG, keep=lambda g: g in wanted)
    print(f">> labelled barcodes in the three populations: {n_lab}; "
          f"not found in barcodes.tsv: {n_unmatched}")
    raw_counts = defaultdict(int)
    for m in meta.values():
        g = (m.get(LABEL_COL) or "").strip()
        if g in wanted:
            raw_counts[g] += 1
    print(">> raw label counts:", dict(raw_counts))

    cs, rs = V.pass_totals(PREFIX, col2group, TAG)
    colsum = V.read_two_col(cs)
    rowsum = V.read_two_col(rs)
    col2cl = {}
    for l in open(col2group):
        c, g = l.split()
        col2cl[int(c)] = g
    cols = sorted(colsum)
    print(">> counts after requiring >=1 UMI:",
          {g: sum(1 for c in cols if col2cl[c] == g) for g in sorted(wanted)})

    # pull the rows needed for panels and controls, over seeds, then score
    per_seed_scores = []
    trpc5_row = name2row.get("Trpc5")
    expr = defaultdict(dict)
    for seed in range(C.RANDOM_SEED, C.RANDOM_SEED + C.N_SEEDS):
        panel_rows, ctrl_rows = V.choose_controls(name2row, rowsum, len(cols), seed)
        need = {r for rs_ in list(panel_rows.values()) + list(ctrl_rows.values()) for r in rs_}
        need.add(trpc5_row)
        new = need - set().union(*[set(d) for d in expr.values()]) if expr else need
        if new:
            e = V.pass_rows(PREFIX, need, TAG, f"seed{seed}")
            for l in open(e):
                r, c, v = l.split()
                c = int(c)
                if c in col2cl:
                    expr[c][int(r)] = float(v)
        per_seed_scores.append(V.module_scores(expr, cols, colsum, panel_rows, ctrl_rows))
        print(f"   scored seed {seed}")

    grp = np.array([col2cl[c] for c in cols])
    inv = {v: k for k, v in C.RUPP_POPS.items()}
    name = np.array([inv[g] for g in grp])

    # Trpc5 as log-normalized expression, not detection
    t5 = np.array([np.log1p(expr.get(c, {}).get(trpc5_row, 0.0) / colsum[c] * C.NORM_TARGET)
                   for c in cols])

    def axis(lab, seed_i=0):
        return per_seed_scores[seed_i][lab]

    rows, percell = [], []
    for pop in ("POMC", "AgRP", "VMH_SF1"):
        m = name == pop
        rec = {"population": pop, "cluster": C.RUPP_POPS[pop], "n": int(m.sum())}
        for lab in ("HYPER", "DEPOL", "HminusD"):
            across = np.array([axis(lab, i)[m].mean() for i in range(C.N_SEEDS)])
            mean, lo, hi = V.boot_ci(axis(lab)[m], seed=1)
            rec[f"{lab}_mean"] = round(mean, 4)
            rec[f"{lab}_ci_lo"] = round(lo, 4)
            rec[f"{lab}_ci_hi"] = round(hi, 4)
            rec[f"{lab}_seed_range"] = round(float(across.max() - across.min()), 4)
        mean, lo, hi = V.boot_ci(t5[m], seed=1)
        rec.update({"Trpc5_mean": round(mean, 4), "Trpc5_ci_lo": round(lo, 4),
                    "Trpc5_ci_hi": round(hi, 4)})
        rows.append(rec)
        for v_hd, v_t5 in zip(axis("HminusD")[m], t5[m]):
            percell.append({"group": pop, "HminusD": round(float(v_hd), 5),
                            "Trpc5": round(float(v_t5), 5)})

    # ------------------------------------------------- contrasts, within each seed
    def vals(lab, seed_i):
        return t5 if lab == "Trpc5" else per_seed_scores[seed_i][lab]

    def contrast(a, b, lab):
        ma, mb = name == a, name == b
        n_eff = 1 if lab == "Trpc5" else C.N_SEEDS
        diffs = np.array([float(vals(lab, i)[ma].mean() - vals(lab, i)[mb].mean())
                          for i in range(n_eff)])
        obs, p = V.perm_test(vals(lab, 0)[ma], vals(lab, 0)[mb], seed=2)
        out_ = {"diff_seed1": round(float(obs), 4), "perm_p_seed1": p,
                "diff_mean_over_seeds": round(float(diffs.mean()), 4),
                "diff_min": round(float(diffs.min()), 4),
                "diff_max": round(float(diffs.max()), 4),
                "diff_range": round(float(diffs.max() - diffs.min()), 4),
                "sign_consistent_across_seeds": bool(
                    np.all(np.sign(diffs) == np.sign(diffs[0])) and diffs[0] != 0),
                "n_seeds": int(n_eff)}
        out_["abs_mean_over_range"] = (
            round(abs(out_["diff_mean_over_seeds"]) / out_["diff_range"], 2)
            if out_["diff_range"] > 0 else None)
        return out_

    stats = {
        "POMC_vs_AgRP_Trpc5": contrast("POMC", "AgRP", "Trpc5"),
        "POMC_vs_AgRP_HminusD": contrast("POMC", "AgRP", "HminusD"),
        "VMHSF1_vs_POMC_HminusD": contrast("VMH_SF1", "POMC", "HminusD"),
        "VMHSF1_vs_AgRP_HminusD": contrast("VMH_SF1", "AgRP", "HminusD"),
        "VMHSF1_vs_POMC_Trpc5": contrast("VMH_SF1", "POMC", "Trpc5"),
    }

    # leave-one-out over panel genes
    loo = {}
    for panel in ("HYPER", "DEPOL"):
        genes = [g for g in C.PANELS[panel] if g in name2row]
        keep = 0
        for drop in genes:
            pr, cr = V.choose_controls(name2row, rowsum, len(cols), C.RANDOM_SEED)
            pr = {p: [r for r in rs_ if r != name2row[drop]] for p, rs_ in pr.items()}
            sc = V.module_scores(expr, cols, colsum, pr, cr)["HminusD"]
            if sc[name == "VMH_SF1"].mean() > sc[name == "POMC"].mean():
                keep += 1
        loo[panel] = {"n_genes": len(genes), "n_preserving_VMH_most_hyperpolarizing": keep}

    gate_ok = (stats["POMC_vs_AgRP_Trpc5"]["diff_seed1"] > 0 and
               stats["POMC_vs_AgRP_HminusD"]["diff_mean_over_seeds"] < 0)
    out = {"accession": C.ACCESSIONS["Rupp_LepRb"], "label_column": LABEL_COL,
           "n_control_gene_seeds": C.N_SEEDS,
           "raw_label_counts": dict(raw_counts),
           "n_after_umi_filter": {r["population"]: r["n"] for r in rows},
           "contrasts": stats, "leave_one_out": loo,
           "validation_gate_passed": bool(gate_ok),
           "note": ("sign_consistent_across_seeds and diff_range carry the robustness claim; "
                    "the permutation p is from the first seed only and does not account "
                    "for control-gene selection.")}

    for nm, data, fields in (
            ("crosspop_scores.csv", rows, list(rows[0].keys())),
            ("crosspop_percell_stats.csv", percell, ["group", "HminusD", "Trpc5"])):
        with open(os.path.join(C.RESULTS, nm), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=fields)
            w.writeheader(); w.writerows(data)
    json.dump(out, open(os.path.join(C.RESULTS, "crosspop.json"), "w"), indent=2)
    print(json.dumps(out["contrasts"], indent=2))
    if not gate_ok:
        sys.exit("VALIDATION GATE FAILED: POMC did not score more depolarizing with higher "
                 "Trpc5 than AgRP. Do not place the VMH population until this is resolved.")


if __name__ == "__main__":
    main()
