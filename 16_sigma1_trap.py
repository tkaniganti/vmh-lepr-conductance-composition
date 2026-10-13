#!/usr/bin/env python3
"""
16_sigma1_trap.py - does the AP-1 sigma composition hold in ribosome-associated mRNA?

The gradient established in 12 rests entirely on nuclear counts. Two claims in this project
already died when the TRAP supernatant fractions were brought in: Kcnj11 "clearly expressed"
(unenriched, 0.07 log2) and the Trpc5 SF1-over-LepR ratio (unenriched in both fractions). A
composition claim resting on one modality is exposed in exactly the same way.

The SF1 series is the clean one: three bead and three supernatant replicates, verified
assignment, and a working immunoprecipitation (Nr5a1 at 3.84 log2). The Lepr series has one
sample per cell-by-treatment combination and is arcuate-contaminated (Pomc enriched tenfold),
so it is reported as point estimates with that caveat and carries no inference.

Four questions:
  1. Are the sigma paralogs detectable in the translatome at all, and at what level?
  2. Is the per-replicate sigma1A share in the SF1 bead fraction consistent with the nuclear
     estimate for VMH SF1 (0.545 in Rupp L3.Nr5a1, 0.610 in the Affinati Lepr class)?
  3. Does the share differ between bead and supernatant? The supernatant is everything else
     in the same dissection, so a shift means the composition is a property of SF1 neurons
     rather than of the tissue.
  4. Is AP1-small the position with the largest bead-versus-supernatant shift, as it was the
     largest across populations in the nuclear data?

Per-replicate shares are the unit of analysis, so n = 3 per fraction is respected rather than
pooled into a single ratio.

Output: results/sigma1_trap.json, sigma1_trap_genes.csv, sigma1_trap_shares.csv
"""
import os, csv, gzip, json
import numpy as np
from scipy import stats
import config as C
import vmhlib as V

POSITIONS = {
    ("AP1", "large"):  ["Ap1g1", "Ap1g2"],
    ("AP1", "small"):  ["Ap1s1", "Ap1s2", "Ap1s3"],
    ("AP1", "medium"): ["Ap1m1", "Ap1m2"],
    ("AP2", "large"):  ["Ap2a1", "Ap2a2"],
    ("AP3", "beta"):   ["Ap3b1", "Ap3b2"],
    ("AP3", "medium"): ["Ap3m1", "Ap3m2"],
    ("AP3", "small"):  ["Ap3s1", "Ap3s2"],
}
FOCUS = ("AP1", "small")
# nuclear estimates for the VMH SF1 population, for the cross-modality comparison
SNRNA_REFERENCE = {"Rupp L3.Nr5a1": 0.5448, "Affinati Lepr class": 0.6103}
MIN_CPM = 1.0        # below this a gene is reported but not used for inference


def load_cpm(path, samples):
    with gzip.open(path, "rt") as f:
        r = csv.reader(f)
        hdr = [h.strip().strip('"') for h in next(r)]
        gi = hdr.index("gene_name")
        missing = [s for s in samples if s not in hdr]
        if missing:
            raise SystemExit(f"columns {missing} missing from {os.path.basename(path)}")
        ci = {s: hdr.index(s) for s in samples}
        raw, tot = {}, {s: 0.0 for s in samples}
        for row in r:
            g = row[gi].strip('"')
            d = raw.setdefault(g, {s: 0.0 for s in samples})
            for s, i in ci.items():
                v = float(row[i]); d[s] += v; tot[s] += v
    return {g: {s: d[s] / tot[s] * C.CPM_TARGET for s in samples} for g, d in raw.items()}


def per_replicate_share(cpm, genes, samples, lead=0):
    """Share of genes[lead] within the position, computed separately in each replicate."""
    out = []
    for s in samples:
        vals = [cpm.get(g, {}).get(s, 0.0) for g in genes]
        t = sum(vals)
        out.append(vals[lead] / t if t > 0 else np.nan)
    return np.array(out)


def main():
    bead, sup = C.TRAP_NR5A1["bead"], C.TRAP_NR5A1["sup"]
    nr = load_cpm(V.dpath("Nr5a1_TRAP.csv.gz"), bead + sup)
    lp_samples = [C.TRAP_LEPR[k] for k in ("bead_pbs", "sup_pbs", "bead_leptin", "sup_leptin")]
    lp = load_cpm(V.dpath("Lepr_TRAP.csv.gz"), lp_samples)

    allg = sorted({g for gs in POSITIONS.values() for g in gs}) + ["Nr5a1", "Lepr"]

    # ------------------------------------------------------------ 1. per-gene levels
    generows = []
    for g in allg:
        b = np.array([nr.get(g, {}).get(s, 0.0) for s in bead])
        s_ = np.array([nr.get(g, {}).get(s, 0.0) for s in sup])
        try:
            t, p = stats.ttest_ind(np.log2(b + 1), np.log2(s_ + 1), equal_var=False)
            p = float(p)
        except Exception:
            p = float("nan")
        generows.append({
            "gene": g,
            "SF1_bead_cpm": round(float(b.mean()), 2),
            "SF1_bead_sd": round(float(b.std(ddof=1)), 2) if b.size > 1 else None,
            "SF1_sup_cpm": round(float(s_.mean()), 2),
            "log2_IP_over_Sup": round(float(np.log2((b.mean() + 1) / (s_.mean() + 1))), 3),
            "welch_p": None if np.isnan(p) else round(p, 4),
            "above_min_cpm": bool(b.mean() >= MIN_CPM),
            "Lepr_bead_PBS_cpm": round(lp.get(g, {}).get(C.TRAP_LEPR["bead_pbs"], 0.0), 2),
            "Lepr_sup_PBS_cpm": round(lp.get(g, {}).get(C.TRAP_LEPR["sup_pbs"], 0.0), 2)})
    with open(os.path.join(C.RESULTS, "sigma1_trap_genes.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(generows[0].keys()))
        w.writeheader(); w.writerows(generows)

    sigma = POSITIONS[FOCUS]
    detectable = {r["gene"]: r["above_min_cpm"] for r in generows}
    print(">> sigma paralogs in the SF1 bead fraction (CPM):",
          {g: next(r["SF1_bead_cpm"] for r in generows if r["gene"] == g) for g in sigma})
    if not all(detectable.get(g) for g in sigma):
        print("   WARNING: at least one sigma paralog is below "
              f"{MIN_CPM} CPM; shares below are unreliable")

    # ----------------------------------------------- 2,3. shares, per replicate
    b_share = per_replicate_share(nr, sigma, bead)
    s_share = per_replicate_share(nr, sigma, sup)
    t, p = stats.ttest_ind(b_share, s_share, equal_var=False)
    sharerows = []
    for i, g in enumerate(sigma):
        bs = per_replicate_share(nr, sigma, bead, lead=i)
        ss = per_replicate_share(nr, sigma, sup, lead=i)
        sharerows.append({"gene": g, "fraction": "SF1_bead",
                          "share_mean": round(float(np.nanmean(bs)), 4),
                          "share_sd": round(float(np.nanstd(bs, ddof=1)), 4),
                          "replicates": [round(float(x), 4) for x in bs]})
        sharerows.append({"gene": g, "fraction": "SF1_sup",
                          "share_mean": round(float(np.nanmean(ss)), 4),
                          "share_sd": round(float(np.nanstd(ss, ddof=1)), 4),
                          "replicates": [round(float(x), 4) for x in ss]})
    for key, samp in (("Lepr_bead_PBS", C.TRAP_LEPR["bead_pbs"]),
                      ("Lepr_sup_PBS", C.TRAP_LEPR["sup_pbs"]),
                      ("Lepr_bead_Leptin", C.TRAP_LEPR["bead_leptin"])):
        for i, g in enumerate(sigma):
            v = per_replicate_share(lp, sigma, [samp], lead=i)
            sharerows.append({"gene": g, "fraction": key,
                              "share_mean": round(float(v[0]), 4) if not np.isnan(v[0]) else None,
                              "share_sd": None, "replicates": [round(float(v[0]), 4)]})
    with open(os.path.join(C.RESULTS, "sigma1_trap_shares.csv"), "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["gene", "fraction", "share_mean", "share_sd",
                                           "replicates"])
        w.writeheader(); w.writerows(sharerows)

    # ------------------------------------------------- 4. specificity across positions
    posrows = []
    for (cx, pos), genes in POSITIONS.items():
        bs = per_replicate_share(nr, genes, bead)
        ss = per_replicate_share(nr, genes, sup)
        lead_cpm = float(np.mean([nr.get(genes[0], {}).get(s, 0.0) for s in bead]))
        if np.isnan(bs).any() or np.isnan(ss).any():
            continue
        tt, pp = stats.ttest_ind(bs, ss, equal_var=False)
        posrows.append({"complex": cx, "position": pos, "lead_paralog": genes[0],
                        "bead_share": round(float(bs.mean()), 4),
                        "sup_share": round(float(ss.mean()), 4),
                        "shift": round(float(bs.mean() - ss.mean()), 4),
                        "welch_p": round(float(pp), 4),
                        "lead_cpm_bead": round(lead_cpm, 1)})
    posrows.sort(key=lambda r: -abs(r["shift"]))
    rank = next((i + 1 for i, r in enumerate(posrows)
                 if (r["complex"], r["position"]) == FOCUS), None)

    out = {
        "fraction_assignment": {"bead": bead, "sup": sup,
                                "source": "GSE172206_series_matrix.txt.gz"},
        "ip_works": {"Nr5a1_log2_IP_over_Sup":
                     next(r["log2_IP_over_Sup"] for r in generows if r["gene"] == "Nr5a1")},
        "sigma_cpm_bead": {g: next(r["SF1_bead_cpm"] for r in generows if r["gene"] == g)
                           for g in sigma},
        "all_sigma_above_min_cpm": bool(all(detectable.get(g) for g in sigma)),
        "sigma1A_share": {
            "SF1_bead_mean": round(float(np.nanmean(b_share)), 4),
            "SF1_bead_replicates": [round(float(x), 4) for x in b_share],
            "SF1_sup_mean": round(float(np.nanmean(s_share)), 4),
            "bead_vs_sup_welch_p": round(float(p), 4),
            "snRNA_reference": SNRNA_REFERENCE},
        "position_shift_ranking": posrows,
        "AP1_small_rank_by_shift": rank,
        "lepr_series_caveat": ("one sample per cell-by-treatment combination and the fraction "
                               "is arcuate-contaminated (Pomc enriched tenfold), so Lepr "
                               "values are descriptive only"),
        "interpretation": (
            "The nuclear gradient is corroborated if the sigma paralogs are detectable, the "
            "SF1 bead sigma1A share falls near the nuclear estimate for VMH SF1, and the "
            "bead share differs from the supernatant share, which would show the composition "
            "belongs to SF1 neurons rather than the surrounding tissue. A bead share equal to "
            "the supernatant share means the composition is a tissue property, not a "
            "cell-type one, and the claim must be narrowed accordingly.")}

    json.dump(out, open(os.path.join(C.RESULTS, "sigma1_trap.json"), "w"), indent=2)
    print(json.dumps({k: v for k, v in out.items()
                      if k not in ("interpretation", "lepr_series_caveat")}, indent=2))


if __name__ == "__main__":
    main()
