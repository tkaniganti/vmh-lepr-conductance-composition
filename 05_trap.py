#!/usr/bin/env python3
"""
05_trap.py - translatome cross-check, with the IP-over-Sup enrichment the design provides.

The previous version used only the bead (IP) columns. A CPM value in an IP does not
establish cell-type expression, because an unenriched transcript can sit in the IP as
background. Both Affinati TRAP series deposit the matching supernatant fraction, which
is what the comparison needs:

  GSE172206 (SF1/Nr5a1): Bead B1-B3 = 116053/4/5, Sup S1-S3 = 116056/7/8   (n = 3 per fraction)
  GSE172200 (Lepr):      66853 Bead/PBS, 66854 Sup/PBS, 66855 Bead/Leptin, 66856 Sup/Leptin
                         (n = 1 per cell-by-treatment combination; point estimates only)

Fraction assignment is read from config, which was verified against the series matrix
rather than inferred from column position.

Output: results/trap_cpm.csv, results/trap_enrichment.csv, results/trap.json
"""
import os, csv, gzip, json
import numpy as np
from scipy import stats
import config as C
import vmhlib as V


def load_cpm(path, wanted):
    """Return {gene: {sample: CPM}} for the named sample columns, duplicate symbols summed."""
    with gzip.open(path, "rt") as f:
        r = csv.reader(f)
        hdr = [h.strip().strip('"') for h in next(r)]
        gi = hdr.index("gene_name")
        missing = [s for s in wanted if s not in hdr]
        if missing:
            raise SystemExit(f"columns {missing} not in {os.path.basename(path)}; header = {hdr}")
        ci = {s: hdr.index(s) for s in wanted}
        raw, tot = {}, {s: 0.0 for s in wanted}
        for row in r:
            g = row[gi].strip('"')
            d = raw.setdefault(g, {s: 0.0 for s in wanted})
            for s, i in ci.items():
                v = float(row[i])
                d[s] += v
                tot[s] += v
    return {g: {s: d[s] / tot[s] * C.CPM_TARGET for s in wanted} for g, d in raw.items()}


def main():
    nr_cols = C.TRAP_NR5A1["bead"] + C.TRAP_NR5A1["sup"]
    nr = load_cpm(V.dpath("Nr5a1_TRAP.csv.gz"), nr_cols)
    lp_cols = [C.TRAP_LEPR[k] for k in ("bead_pbs", "sup_pbs", "bead_leptin", "sup_leptin")]
    lp = load_cpm(V.dpath("Lepr_TRAP.csv.gz"), lp_cols)

    genes = C.KEY_CHANNELS + ["Lepr", "Nr5a1", "Slc17a6", "Gck", "Pomc", "Agrp"]
    cpm_rows, enr_rows = [], []
    for g in genes:
        nb = np.array([nr.get(g, {}).get(s, 0.0) for s in C.TRAP_NR5A1["bead"]])
        ns = np.array([nr.get(g, {}).get(s, 0.0) for s in C.TRAP_NR5A1["sup"]])
        lbp = lp.get(g, {}).get(C.TRAP_LEPR["bead_pbs"], 0.0)
        lsp = lp.get(g, {}).get(C.TRAP_LEPR["sup_pbs"], 0.0)
        lbl = lp.get(g, {}).get(C.TRAP_LEPR["bead_leptin"], 0.0)
        lsl = lp.get(g, {}).get(C.TRAP_LEPR["sup_leptin"], 0.0)
        cpm_rows.append({"gene": g,
                         "SF1_Bead_mean": round(float(nb.mean()), 1),
                         "SF1_Sup_mean": round(float(ns.mean()), 1),
                         "Lepr_Bead_PBS": round(lbp, 1), "Lepr_Sup_PBS": round(lsp, 1),
                         "Lepr_Bead_Leptin": round(lbl, 1), "Lepr_Sup_Leptin": round(lsl, 1)})
        # SF1 has replicates, so a test is possible; Lepr does not, so report the ratio only.
        try:
            t, p = stats.ttest_ind(np.log2(nb + 1), np.log2(ns + 1), equal_var=False)
            p = float(p)
        except Exception:
            t, p = float("nan"), float("nan")
        enr_rows.append({
            "gene": g,
            "SF1_log2_IP_over_Sup": round(float(np.log2((nb.mean() + 1) / (ns.mean() + 1))), 2),
            "SF1_welch_p_unadjusted": None if np.isnan(p) else round(p, 4),
            "Lepr_PBS_log2_IP_over_Sup": round(float(np.log2((lbp + 1) / (lsp + 1))), 2),
            "Lepr_Leptin_log2_IP_over_Sup": round(float(np.log2((lbl + 1) / (lsl + 1))), 2),
            "Lepr_leptin_vs_pbs_log2_IP": round(float(np.log2((lbl + 1) / (lbp + 1))), 2),
            "n_Lepr_replicates": 1,
        })

    for name, rows in (("trap_cpm.csv", cpm_rows), ("trap_enrichment.csv", enr_rows)):
        with open(os.path.join(C.RESULTS, name), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
            w.writeheader(); w.writerows(rows)

    # The SF1-versus-LepR Trpc5 contrast is the reconciliation of Dhillon (SF1 depolarized)
    # with the bidirectional observation in the SF1-and-LepR intersection, so report it directly.
    def gv(d, g, s):
        return d.get(g, {}).get(s, 0.0)
    t5_sf1 = float(np.mean([gv(nr, "Trpc5", s) for s in C.TRAP_NR5A1["bead"]]))
    t5_lepr = gv(lp, "Trpc5", C.TRAP_LEPR["bead_pbs"])
    summary = {
        "fraction_assignment_source": "GSE172206_series_matrix.txt.gz, fraction characteristic",
        "trpc5_SF1_IP_cpm": round(t5_sf1, 1),
        "trpc5_Lepr_IP_PBS_cpm": round(t5_lepr, 1),
        "trpc5_SF1_over_Lepr": round(t5_sf1 / max(t5_lepr, 1e-9), 2),
        "caveat": ("Lepr series has one sample per cell-by-treatment combination, so Lepr "
                   "enrichment and leptin effects are point estimates with no replication."),
    }
    json.dump(summary, open(os.path.join(C.RESULTS, "trap.json"), "w"), indent=2)
    print(json.dumps(summary, indent=2))
    print("\nWritten to", C.RESULTS)


if __name__ == "__main__":
    main()
