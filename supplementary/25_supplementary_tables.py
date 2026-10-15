#!/usr/bin/env python3
"""
25_supplementary_tables.py - build the six supplementary tables the manuscript cites.

Every value is read from results/ rather than typed, so the workbook regenerates with the
analysis. Nothing here is computed; this is a formatting step.

  S1  doublet sensitivity, full set against doublet-free        doublet_sensitivity.json
  S2  rank-based effect sizes superseded by the nested model    conserved_classes.json
  S3  dispersion of the balance axis across populations         dispersion.json
  S4  class-by-cluster membership                               class_membership.csv
  S5  per-class excess across minimum-cluster-size thresholds   mincluster_sensitivity.json
  S6  observed and chance values underlying the excesses        class_nulls.json

Output: results/Supplementary_Tables.xlsx and a plain-text copy of each sheet as CSV.
"""
import os, csv, json
import openpyxl
from openpyxl.styles import Font, Alignment, PatternFill
from openpyxl.utils import get_column_letter
import config as cfg

R = cfg.RESULTS
HEAD = Font(name="Arial", bold=True, size=10, color="FFFFFF")
BODY = Font(name="Arial", size=10)
TITLE = Font(name="Arial", bold=True, size=11)
NOTE = Font(name="Arial", size=9, italic=True, color="555555")
FILL = PatternFill("solid", fgColor="44546A")


def load(fn):
    p = os.path.join(R, fn)
    if not os.path.exists(p):
        print(f"   !! missing {fn}; sheet will be written empty")
        return None
    if fn.endswith(".json"):
        return json.load(open(p))
    with open(p) as fh:
        return list(csv.DictReader(fh))


def sheet(wb, name, title, note, header, rows):
    ws = wb.create_sheet(name)
    ws["A1"] = title
    ws["A1"].font = TITLE
    ws["A2"] = note
    ws["A2"].font = NOTE
    ws.append([])
    ws.append(header)
    for c in range(1, len(header) + 1):
        cell = ws.cell(row=4, column=c)
        cell.font, cell.fill = HEAD, FILL
        cell.alignment = Alignment(horizontal="center", wrap_text=True)
    for r in rows:
        ws.append(r)
    for r in ws.iter_rows(min_row=5):
        for c in r:
            c.font = BODY
    for c in range(1, len(header) + 1):
        width = max([len(str(header[c - 1]))] +
                    [len(str(r[c - 1])) for r in rows if c - 1 < len(r)] + [8])
        ws.column_dimensions[get_column_letter(c)].width = min(width + 3, 46)
    ws.freeze_panes = "A5"
    return ws


def main():
    wb = openpyxl.Workbook()
    wb.remove(wb.active)

    # ---------------------------------------------------------------- S1
    d = load("doublet_sensitivity.json")
    rows = []
    if d:
        a, b = d["with_doublets"], d["without_doublets"]
        rows = [
            ["LepR nuclei", a["n_lepr"], b["n_lepr"]],
            ["VMH nuclei", a["n_vmh"], b["n_vmh"]],
            ["Nr5a1 detected in VMH (%)", a["detection_pct_vmh"]["Nr5a1"],
             b["detection_pct_vmh"]["Nr5a1"]],
            ["Lepr detected in VMH (%)", a["detection_pct_vmh"]["Lepr"],
             b["detection_pct_vmh"]["Lepr"]],
            ["Slc17a6 detected in VMH (%)", a["detection_pct_vmh"]["Slc17a6"],
             b["detection_pct_vmh"]["Slc17a6"]],
            ["Double-positive gate retained (n)", a["naive_gate_lepr"], b["naive_gate_lepr"]],
            ["Double-positive gate retained (%)", a["naive_gate_lepr_pct"],
             b["naive_gate_lepr_pct"]],
            ["Effect size, combined axis", a["eps2"]["HminusD"]["mean"],
             b["eps2"]["HminusD"]["mean"]],
            ["Effect size, depolarising module", a["eps2"]["DEPOL"]["mean"],
             b["eps2"]["DEPOL"]["mean"]],
            ["Effect size, hyperpolarising module", a["eps2"]["HYPER"]["mean"],
             b["eps2"]["HYPER"]["mean"]],
            ["Mixture components selected by BIC", a["mixture"]["best_k"],
             b["mixture"]["best_k"]],
            ["Dip test p", a["dip"]["p"], b["dip"]["p"]],
            ["Trpc5 depth-controlled LR", a["trpc5_depth_controlled"]["LR"],
             b["trpc5_depth_controlled"]["LR"]],
            ["Trpc5 depth-controlled p", a["trpc5_depth_controlled"]["p"],
             b["trpc5_depth_controlled"]["p"]]]
    sheet(wb, "S1 doublets",
          "Supplementary Table S1. Sensitivity to doublet-annotated nuclei",
          "461 of 2,913 Lepr-cluster nuclei (15.8%) carry doublet assignments in the "
          "deposited Conserved_cluster annotation. Full-set values are reported in the "
          "manuscript.",
          ["Quantity", "Full set", "Doublet-free"], rows)

    # ---------------------------------------------------------------- S2
    d = load("conserved_classes.json")
    rows = []
    if d:
        for tag in ("raw", "depth_residualized"):
            p = d["variance_partition"][tag]
            for key, lab in (("eps2_between_conserved_classes", "between conserved classes"),
                             ("eps2_between_clusters_within_lepr_class",
                              "between clusters within the Lepr class")):
                v = p[key]
                rows.append([lab, tag.replace("_", " "), v["mean"], v["min"], v["max"]])
    sheet(wb, "S2 rank effect sizes",
          "Supplementary Table S2. Rank-based effect sizes, superseded by the nested model",
          "Epsilon-squared from Kruskal-Wallis, computed separately at each level. These "
          "are not comparable between levels and are reported here only because an earlier "
          "version of this analysis used them; the nested components of variance in the "
          "manuscript supersede them.",
          ["Level", "Depth handling", "Mean", "Min across draws", "Max across draws"], rows)

    # ---------------------------------------------------------------- S3
    d = load("dispersion.json")
    rows = []
    if d:
        sd = d["depth_residualized"]["sd_by_population"]
        for pop in sorted(sd):
            v = sd[pop]
            rows.append([pop.replace("_", " "), v["mean"], v["min"], v["max"]])
        r = d.get("equal_n", {}).get("sd_ratio_median", {})
        for k, v in r.items():
            rows.append([k.replace("_", " "), v, "", ""])
    sheet(wb, "S3 dispersion",
          "Supplementary Table S3. Dispersion of the effector-balance axis by population",
          "Standard deviation of the depth-residualised balance axis, with the range across "
          "20 control-gene draws. Ratios are medians over 500 equal-n downsampling draws.",
          ["Population or ratio", "Mean", "Min", "Max"], rows)

    # ---------------------------------------------------------------- S4
    d = load("class_membership.csv")
    rows = [[r["conserved_class"], r["vmh_cluster"], int(r["n_nuclei"])] for r in d] if d else []
    sheet(wb, "S4 class membership",
          "Supplementary Table S4. Class-by-cluster membership",
          "Membership of a conserved class is not restricted to similarly named clusters. "
          "Nuclei with a doublet assignment are excluded.",
          ["Conserved class", "VMH cluster", "Nuclei"], rows)

    # ---------------------------------------------------------------- S5
    d = load("mincluster_sensitivity.json")
    rows = []
    if d:
        for cl in sorted(d["by_class"]):
            for thr in d["thresholds"]:
                rec = d["by_class"][cl].get(str(thr), {})
                if "excess" in rec:
                    rows.append([cl, thr, rec["n_clusters"], rec["n_nuclei"],
                                 rec["observed_mean"], rec["null_mean"], rec["excess"]])
                else:
                    rows.append([cl, thr, rec.get("n_clusters", ""),
                                 rec.get("n_nuclei", ""), "", "", "not estimable"])
    sheet(wb, "S5 cluster size thresholds",
          "Supplementary Table S5. Within-class excess across minimum-cluster-size thresholds",
          "Clusters below the threshold are excluded and the chance expectation is "
          "recomputed. A variance component across two or three groups is barely estimable, "
          "so values at thresholds that reduce a class to that few clusters cannot be read "
          "as evidence either way.",
          ["Class", "Minimum cluster size", "Clusters retained", "Nuclei retained",
           "Observed (%)", "Chance (%)", "Excess (%)"], rows)

    # ---------------------------------------------------------------- S6
    d = load("class_nulls.json")
    rows = []
    if d:
        bc = d["by_class"]
        for cl in sorted(bc, key=lambda c: -bc[c]["excess_over_null"]):
            v = bc[cl]
            rows.append([cl, v["n_clusters"], v["n_nuclei"],
                         v["nuclei_per_cluster_mean"], v["nuclei_per_cluster_median"],
                         v["observed_pct"]["mean"], v["observed_pct"]["min"],
                         v["observed_pct"]["max"], v["observed_pct"]["n_draws_at_zero"],
                         v["null_pct"]["mean"], v["null_pct"]["p95"],
                         v["excess_over_null"], v["loo_min"], v["loo_max"]])
    sheet(wb, "S6 per-class nulls",
          "Supplementary Table S6. Observed and chance values underlying the per-class excesses",
          "The chance expectation is computed separately for each class by shuffling cluster "
          "labels among the same nuclei, because it depends on cluster count and nuclei per "
          "cluster. Leave-one-out columns give the range of the estimate when each cluster "
          "in turn is removed, on the first control-gene draw.",
          ["Class", "Clusters", "Nuclei", "Nuclei per cluster (mean)",
           "Nuclei per cluster (median)", "Observed mean (%)", "Observed min (%)",
           "Observed max (%)", "Draws at zero", "Chance mean (%)", "Chance 95th pct (%)",
           "Excess (%)", "Leave-one-out min (%)", "Leave-one-out max (%)"], rows)

    # ---------------------------------------------------------------- S7
    d = load("pomc_replication.json")
    rows = []
    if d:
        inv = d.get("inventory", {})
        t1 = d.get("test1_discreteness", {})
        t3 = d.get("test3_scan", {})
        t4 = d.get("test4_partition") or {}
        rows.append(["Cells analysed", inv.get("n_cells", ""), ""])
        for k, v in (t1.get("bic") or {}).items():
            rows.append([f"Mixture BIC, k = {k}", v, ""])
        rows.append(["Components selected by BIC", t1.get("best_k", ""), ""])
        rows.append(["Dip test D", t1.get("dip_D", ""), ""])
        rows.append(["Dip test p", t1.get("dip_p", ""), ""])
        for lv in (d.get("test2_detectability", {}).get("levels") or []):
            rows.append([f"Split of {lv['split_sd']} SD detected in",
                         f"{lv['pct_detected']}%", "of 100 simulations"])
        rows.append(["Transcripts tested", t3.get("n_tested", ""), ""])
        rows.append(["Transcripts at q < 0.05", t3.get("n_q05", ""), ""])
        rows.append(["Largest partial correlation", t3.get("max_abs_rho", ""),
                     (t3.get("top") or [{}])[0].get("gene", "")])
        rows.append(["Largest variance explained",
                     f"{t3.get('max_pct_variance', '')}%", "uncorrected maximum"])
        if t4:
            rows.append([f"Between-subgroup variance ({t4.get('field','')}, k = {t4.get('k','')})",
                         f"{t4.get('pct_between','')}%", ""])
        else:
            rows.append(["Between-subgroup variance", "not computed",
                         "no deposited sub-labels for this population"])
    sheet(wb, "S7 POMC calibration",
          "Supplementary Table S7. Unsupervised structure in arcuate POMC neurons",
          "The same test applied to a population known from intersectional targeting to "
          "contain functionally distinct subpopulations (Biglari et al., 2021). It likewise "
          "selects a single component, which calibrates the detectability bound against a "
          "real case rather than only against simulation. The bound here is weaker than in "
          "the VMH comparison because the sample is smaller.",
          ["Quantity", "Value", "Note"], rows)

    out = os.path.join(R, "Supplementary_Tables.xlsx")
    wb.save(out)
    print("Written", out)
    for ws in wb.worksheets:
        print(f"   {ws.title}: {ws.max_row - 4} rows")


if __name__ == "__main__":
    main()
